module Lattices
    use const
    use errors, only: stop_all
    use globals
    use Lattice_Data
    implicit none

    contains

    !Setup the lattice - both in real space and kspace if desired.
    subroutine Setup_Lattice()
        implicit none
        character(len=*), parameter :: t_r='Setup_Lattice'

        if(LatticeDim.eq.2) then
            if(tSquareLatt) then
                if(CellShape.eq.1) then
                    !Original tilted cell
!                    if(tTiltedLattice) then
                    call Setup2DLattice_Tilt()
                else
                    !Square or tilted cell
                    call Setup2DLattice_Square()
                endif
            else
                call stop_all(t_r,'Other lattice geometries not coded up yet')
            endif
        endif
            
        if((tDiag_KSpace.or.tProjectHFKPnts.or.tKSpaceOrbs).and.(.not.allocated(KPnts))) then
            call setup_kspace()
        endif

    end subroutine
    
    !Setup the kpoint mesh, and other things needed to work in kspace
    !See C. Gros, Z. Phys. B - Condensed Matter 86, 359-365 (1992) for details for 1 unit cell.
    !Alternatively, for multi-band models, we need an additional unitary matrix, which we take to be the unit matrix 
    !psi_(n,k) = 1/sqrt(N) \sum_(R,m) exp(ik.R) U_nm phi_m(r-R)
    !where psi_(n,k) are kspace orbitals, R is the unit cell lattice vectors, r is the basis function distance vector, k is kpoint.
    !U_nm is an arbitrary unitary matrix, which determines the mixing between the bands in a given kpoint
    subroutine setup_kspace()
        implicit none
        integer :: SS_Period,i,k,ind_1,ind_2,j,nKPnts_x,nKPnts_y,indx,indy,kpnt,n
        real(dp) :: PrimLattVec(LatticeDim),phase,ddot,r,r2,kpntx,kpnty,PrimLattVec_2(LatticeDim),phase1,phase5
        real(dp) :: tmp1,tmp2,tmp3,tmp4
        complex(dp) :: val,val2
        complex(dp) , allocatable :: temp(:,:),ham_temp(:,:),Random_CorrPot(:,:)
        character(len=*), parameter :: t_r='setup_kspace'

        !SS_Period is the size of the supercell repeating unit (e.g. the coupling correlation potential)
        !This will be equivalent to the number of bands per kpoint
        SS_Period = nImp
        if(mod(nSites,SS_Period).ne.0) call stop_all(t_r,'Lattice dimensions not consistent')
        nKPnts = nSites/SS_Period
        allocate(KPnts(LatticeDim,nKPnts))
        KPnts(:,:) = zero  
        allocate(RecipLattVecs(LatticeDim,LatticeDim))
        RecipLattVecs(:,:) = zero
            
        !Create k-space mesh
        if(LatticeDim.eq.1) then
            !Define the reciprocal lattice vector
            if(mod(nKPnts,2).eq.0) then
                if(tPeriodic) then
                    !We want to use a Gamma centered mesh
                    !I.e. include BZ boundary
                    tShift_Mesh = .false.
                else
                    !We want a Monkhort-Pack mesh
                    tShift_Mesh = .true.
                endif
            else
                !Odd number of k-points
                if(tPeriodic) then
                    !We want to use a Monkhorst pack mesh
                    tShift_Mesh = .true. 
                else
                    !Gamma-centered mesh (i.e. include BZ boundary)
                    tShift_Mesh = .false.
                endif
            endif
!            tShift_Mesh = .true. 
            if(tShift_Mesh) then
                write(6,"(A)") "Using a Monkhorst-pack kpoint mesh of 1st BZ - Symmetric k-points around Gamma point"
            else
                write(6,"(A)") "Using a gamma-centered kpoint mesh of 1st BZ - BZ boundary included"
            endif

            BZVol = 2.0_dp*pi/real(SS_Period,dp)
            write(6,"(A,G21.14)") "Brillouin zone volume: ",BZVol

            RecipLattVecs(1,1) = 2.0_dp*pi/real(SS_Period,dp)
            !Just use equally spaced mesh starting at -pi/SS_Period, and working our way across
            do k = 1,nKPnts
                KPnts(1,k) = -RecipLattVecs(1,1)/2.0_dp + (k-1)*RecipLattVecs(1,1)/nKPnts
                if(tShift_Mesh) then
                    !Shift kpoint mesh by half the kpoint spacing
                    KPnts(1,k) = KPnts(1,k) + RecipLattVecs(1,1)/(2.0_dp*real(nKPnts,dp))
                endif
            enddo
        elseif(LatticeDim.eq.2) then
            if(CellShape.eq.1) then
                call stop_all(t_r,'Cannot setup k-space - impurity tiling is not same as direct lattice. Use square lattice.')
            endif

            !Check that the number of kpoints is a square number
            if(abs(sqrt(real(nKPnts,dp))-real(nint(sqrt(real(nKPnts,dp))),dp)).gt.1.0e-8_dp) then
                call stop_all(t_r,'Number of kpoints is not a square number')
            endif
            nKPnts_x = nint(sqrt(real(nKPnts,dp)))
            nKPnts_y = nint(sqrt(real(nKPnts,dp)))
            write(6,"(A,I8)") "Number of kpoints in each dimension: ",nKPnts_x

            !We should have already defined whether we have periodic or antiperiodic boundary conditions
            if(mod(nKPnts_x,2).eq.0) then
                !This is correct for even kpoint meshes
                if(tPeriodic) tShift_Mesh = .false.
                if(.not.tPeriodic) tShift_Mesh = .true.
            else
                !This is correct for odd kpoint meshes
                if(tPeriodic) tShift_Mesh = .true. 
                if(.not.tPeriodic) tShift_Mesh = .false.
            endif

            if(tShift_Mesh) then
                write(6,"(A)") "Using a Monkhorst-pack kpoint mesh of 1st BZ - Symmetric k-points around Gamma point"
            else
                write(6,"(A)") "Using a gamma-centered kpoint mesh of 1st BZ - BZ boundary included"
            endif

            if(CellShape.eq.2) then
                !Reciprocal lattice vector in the x direction, defining a square grid
                RecipLattVecs(1,1) = 2.0_dp*pi/real(nImp_x,dp)
                RecipLattVecs(2,1) = 0.0_dp
                !Reciprocal lattice vector in the y-direction
                RecipLattVecs(1,2) = 0.0_dp
                RecipLattVecs(2,2) = 2.0_dp*pi/real(nImp_y,dp)
                BZVol = (2.0_dp*pi/real(nImp_x,dp))*(2.0_dp*pi/real(nImp_y,dp))
            elseif(CellShape.eq.3) then
                RecipLattVecs(1,1) = 2.0_dp*pi/real(nImp_x,dp)
                RecipLattVecs(2,1) = zero 
                RecipLattVecs(1,2) = -2.0_dp*pi/real(nImp_x,dp)
                RecipLattVecs(2,2) = 2.0_dp*pi/real(nImp_y,dp)
                BZVol = (2.0_dp*pi/real(nImp_x,dp))*(2.0_dp*pi/real(nImp_y,dp))
            endif

            !Reciprocal lattice set up correctly?
            tmp1 = zero
            tmp2 = zero
            tmp3 = zero
            tmp4 = zero
            do l=1,LatticeDim
                tmp1 = tmp1 + LatticeVec(l,1)*RecipLattVecs(l,1)
                tmp2 = tmp2 + LatticeVec(l,2)*RecipLattVecs(l,2)
                tmp3 = tmp3 + LatticeVec(l,1)*RecipLattVecs(l,2)
                tmp4 = tmp4 + LatticeVec(l,2)*RecipLattVecs(l,1)
            enddo
            if(abs(tmp1-2.0_dp*pi).gt.1.0e-8_dp) call stop_all(t_r,'Reciprocal lattice wrong')
            if(abs(tmp2-2.0_dp*pi).gt.1.0e-8_dp) call stop_all(t_r,'Reciprocal lattice wrong')
            if(abs(tmp3).gt.1.0e-8_dp) call stop_all(t_r,'Reciprocal lattice wrong')
            if(abs(tmp4).gt.1.0e-8_dp) call stop_all(t_r,'Reciprocal lattice wrong')

            !Now, to define the kpoints in each dimension
            do i = 1,nKPnts_x

                !Which x-coordinate do we have for this y cut of kpoints
                kpntx = -RecipLattVecs(1,1)/2.0_dp + (i-1)*RecipLattVecs(1,1)/nKPnts_x
                if(tShift_Mesh) kpntx = kpntx + RecipLattVecs(1,1)/(2.0_dp*real(nKPnts_x,dp))

                do j = 1,nKPnts_y
                    kpnt = (i-1)*nKPnts_x + j

                    kpnty = RecipLattVecs(2,2)/2.0_dp - (j-1)*RecipLattVecs(2,2)/nKPnts_y
                    if(tShift_Mesh) kpnty = kpnty - RecipLattVecs(2,2)/(2.0_dp*real(nKPnts_y,dp))

                    KPnts(1,kpnt) = kpntx
                    KPnts(2,kpnt) = kpnty

                enddo
            enddo
        else
            !Quoi?
            call stop_all(t_r,'Error here')
        endif
        write(6,"(A,L1)") "tShift_Mesh: ",tShift_Mesh

!        if(tWriteOut) then
        write(6,"(A)") "Writing out kpoint mesh: "
        do k=1,nKPnts
            if(LatticeDim.eq.1) then
                write(6,*) "KPnt ",k,KPnts(1,k),KPnts(1,k)/RecipLattVecs(1,1)
            else
                write(6,"(A,I5,2F13.5)") "KPnt ",k,KPnts(:,k)
            endif
        enddo
!        endif

        !Setup rotation matrix from site basis to k-space
        !First index r, second k
        allocate(RtoK_Rot(nSites,nSites))
        RtoK_Rot(:,:) = zzero
        !Run though all kpoints
        do k = 1,nKPnts
            !Construct rotation
            ind_1 = ((k-1)*nImp) + 1
            ind_2 = nImp*k
            do n = 0,nImp-1 !Run over bands in a given kpoint

                if(LatticeDim.eq.1) then
                    do j = 1,iImpRepeats
                        site = (j-1)*nImp + n + 1
                        PrimLattVec(1) = real((j-1)*nImp,dp)    !The real-space translation to the cell
                        phase = ddot(LatticeDim,KPnts(:,k),1,PrimLattVec,1)
                        RtoK_Rot(site,ind_1+n) = exp(cmplx(zero,phase,dp))/sqrt(real(nKPnts,dp))
                    enddo
                else
                    do kx = 0,iImpRepeats_x-1
                        PrimLattVec(1) = iImpRepeats_x*LatticeVector(1,1) + iImpRepeats_x*LatticeVector(1,2)
                        do ky = 0,iImpRepeats_y-1
                            PrimLattVec(2) = iImpRepeats_y*LatticeVector(2,1) + iImpRepeats_y*LatticeVector(2,2)
                            !Run through all cells, with PrimLattVec being the coordinate of the cell
                            !Choose the site of the cell = n+1
                            site = (j-1)*nImp + n + 1
                            phase = ddot(LatticeDim,KPnts(:,k),1,PrimLattVec,1)
                            RtoK_Rot(site,ind_1+n) = exp(cmplx(zero,phase,dp))/sqrt(real(nKPnts,dp))
                        enddo
                    enddo
                endif
!                do i = 1,nSites
!                    if(LatticeDim.eq.1) then
!                        if(mod(i-1,nImp).ne.n) cycle    !This is the unit rotation between n and m
!                        PrimLattVec(1) = real(i-1-mod(i-1,nImp))    !The real-space translation to the cell
!!                        PrimLattVec(1) = real(i-1,dp)               !The real-space translation to this site
!                    else
!                        call SiteIndexToLatCoord_2DSquare(i,indx,indy)
!                        !Since these indices are 1 indexed, we need to make them 0 indexed
!                        !The *site* displacement vectors are indx-1 and indy-1. 
!                        !The *cell* displacement vectors are indx-1-mod(indx-1,nImp_x) and the y version.
!                        if(((mod(indx-1,nImp_x)*nImp_x)+mod(indy-1,nImp_y)).ne.n) cycle
!                        PrimLattVec(1) = real(indx - 1 - mod(indx-1,nImp_x),dp)
!                        PrimLattVec(2) = real(indy - 1 - mod(indy-1,nImp_y),dp)
!!                        PrimLattVec(1) = real(indx - 1,dp)
!!                        PrimLattVec(2) = real(indy - 1,dp)
!!                        write(6,"(A,I5,A,2F10.4)") "Lattice site: ",i, " has coordinate: ",PrimLattVec(:)
!                    endif
!                    phase = ddot(LatticeDim,KPnts(:,k),1,PrimLattVec,1)
!!                    RtoK_Rot(i,ind_1+mod(i,nImp)) = exp(dcmplx(zero,phase))/sqrt(real(nKPnts,dp))
!                    RtoK_Rot(i,ind_1+n) = exp(dcmplx(zero,phase))/sqrt(real(nKPnts,dp))
!                enddo
            enddo
        enddo

        if(tWriteOut) call writematrixcomp(RtoK_Rot,'RtoK_Rot',.true.)

!        if(tCheck) then
        if(.true.) then
            !Now, is RtoK_Rot unitary?
            allocate(temp(nSites,nSites))
            !Check unitarity of matrix
            call ZGEMM('C','N',nSites,nSites,nSites,zone,RtoK_Rot,nSites,RtoK_Rot,nSites,zzero,temp,nSites) 
            do i = 1,nSites
                do j = 1,nSites
                    if((i.eq.j).and.(abs(temp(i,j)-zone).gt.1.0e-7_dp)) then
                        call writematrixcomp(temp,'Identity?',.false.)
                        write(6,*) "i,j: ",i,j,temp(i,j)
                        call stop_all(t_r,'Rotation matrix not unitary 1')
                    elseif((i.ne.j).and.(abs(temp(j,i)).gt.1.0e-7_dp)) then
                        call writematrixcomp(temp,'Identity?',.false.)
                        write(6,*) "i,j: ",i,j,temp(i,j)
                        call stop_all(t_r,'Rotation matrix not unitary 2')
                    endif
                enddo
            enddo
            !Try other way...
            call ZGEMM('N','C',nSites,nSites,nSites,zone,RtoK_Rot,nSites,RtoK_Rot,nSites,zzero,temp,nSites) 
            do i = 1,nSites
                do j = 1,nSites
                    if((i.eq.j).and.(abs(temp(i,j)-zone).gt.1.0e-7_dp)) then
                        call writematrixcomp(temp,'Identity?',.true.)
                        write(6,*) "i,j: ",i,j,temp(i,j)
                        call stop_all(t_r,'Rotation matrix not unitary 3')
                    elseif((i.ne.j).and.(abs(temp(j,i)).gt.1.0e-7_dp)) then
                        call writematrixcomp(temp,'Identity?',.true.)
                        write(6,*) "i,j: ",i,j,temp(j,i)
                        call stop_all(t_r,'Rotation matrix not unitary 4')
                    endif
                enddo
            enddo
            write(6,*) "Rotation matrix unitary... :)"
            deallocate(temp)
        endif

        !Right, just as a sanity check, see if this rotation block diagonalizes the hopping matrix
        allocate(ham_temp(nSites,nSites))
!        call writematrix(h0,'h0',.true.)
        do i = 1,nSites
            do j = 1,nSites
                ham_temp(j,i) = cmplx(h0(j,i),zero,dp)
            enddo
        enddo
        !Add a random correlation potential
        allocate(Random_Corrpot(nImp,nImp))
        Random_CorrPot(:,:) = zzero
        do i = 1,nImp
            do j = i,nImp
                if(i.eq.j) then
                    call random_number(r)
                    Random_CorrPot(j,i) = cmplx(r,zero,dp)
                else
                    call random_number(r)
                    call random_number(r2)
                    Random_CorrPot(j,i) = cmplx(r,r2,dp)
                endif
            enddo
        enddo
        call MakeBlockHermitian(Random_CorrPot,nImp)
        call add_localpot_comp_inplace(ham_temp,Random_CorrPot,tAdd=.true.)
!        call writematrixcomp(ham_temp,'h0 with corrpot',.false.)
        !Check hermitian
        do i = 1,nSites
            do j = 1,nSites
                if(i.eq.j) then
                    if(abs(aimag(ham_temp(i,i))).gt.1.0e-8_dp) then
                        call stop_all(t_r,'Matrix not diagonal real')
                    endif
                else
                    if(abs(ham_temp(j,i)-conjg(ham_temp(i,j))).gt.1.0e-8_dp) then
                        call stop_all(t_r,'Matrix not hermitian')
                    endif
                endif
            enddo
        enddo
        allocate(temp(nSites,nSites))
        call ZGEMM('C','N',nSites,nSites,nSites,zone,RtoK_Rot,nSites,ham_temp,nSites,zzero,temp,nSites)
        call ZGEMM('N','N',nSites,nSites,nSites,zone,temp,nSites,RtoK_Rot,nSites,zzero,ham_temp,nSites)
        !ham_temp is now in kspace
        !Zero the diagonal blocks
        do k = 1,nKPnts
            ind_1 = ((k-1)*nImp) + 1
            ind_2 = nImp*k
            ham_temp(ind_1:ind_2,ind_1:ind_2) = zzero
        enddo
        !Now check to see if there are any non-zero elements. If so, the unitary rotation is not correct
        do i = 1,nSites
            do j = 1,nSites
                if(abs(ham_temp(j,i)).gt.1.0e-7_dp) then
                    call writematrixcomp(ham_temp,'zero matrix',.true.)
                    write(6,*) "i,j: ",j,i
                    write(6,*) "ham in kspace: ",ham_temp(j,i)
                    call stop_all(t_r,'kspace rotations not correctly set up. ' &
                        //'One-electron hamiltonian is not block diagonal in kspace')
                endif
            enddo
        enddo
        deallocate(temp,Random_Corrpot,ham_temp)
        write(6,"(A)") "k-space rotations correctly block diagonalize general one-electron matrix"

    end subroutine setup_kspace
        
    
    !Setup arrays and indices needed to define a 2D square, non-tilted lattice
    !CellShape = 2 means square cell
    !CellShape = 3 means
!                      -*-
!                    -*-*-
!                    -*-
    subroutine Setup2DLattice_Square()
        implicit none
        integer :: nSitesOrig,nSites_x_low,nSites_x_high,nBelowSitesChange
        integer :: i,j,k,nAboveSitesChange,PhaseChange,StartInd
        integer :: indx,indy,indx_folded,indy_folded,index_2
        character(len=*), parameter :: t_r='Setup2DLattice_Square'

        write(6,"(A)") "Setting up a square lattice..."
        if(CellShape.eq.2) then
            write(6,*) "Simple cubic unit/supercell setup"
        elseif(CellShape.eq.3) then
            write(6,*) "Tilted lattice cell setup"
        else
            call stop_all(t_r,'Cell shape not recognized')
        endif
        if(LatticeDim.ne.2) call stop_all(t_r,'Should only be in here with 2D lattice')

        nImp_x = nint(sqrt(real(nImp,dp)))
        nImp_y = nImp_x
        if((sqrt(real(nImp,dp))-nint(sqrt(real(nImp,dp)))).gt.1.0e-8_dp) then
            !Constrain impurity cluster to have the same geometry unit cell as the 
            !entire lattice.
            call stop_all(t_r,'Number of impurity sites is not a square number.')
        endif

        LatticeVector(:,:) = zero
        if(CellShape.eq.2) then
            LatticeVector(1,1) = one*real(nImp_x,dp)
            LatticeVector(2,2) = one*real(nImp_y,dp)
            LatticeVector(:,:) = LatticeVector(:,:)
        elseif(CellShape.eq.3) then
            LatticeVector(1,1) = one*real(nImp_y,dp)
            LatticeVector(2,1) = one*real(nImp_x,dp)
            LatticeVector(2,2) = one*real(nImp_x,dp)
        endif
        !TODO: Check that the area of the unit cell is = nImp


        !Right, now how many lattice sites should be have. The constraints are that 
        !   o We need nSites to be a square number. We also need nSites/nImp to be a square number (number of kpoints).
        !   o We want nSites_x to be an even integer multiple of nImp_x
        !   o We want the total number to be as close as possible to the original value

        !First, calculate the two values satisfying 1 and 2, both above and below the original value
        nSitesOrig = nSites

        i = int(sqrt(real(nSites,dp))/real(nImp_x,dp))
        nSites_x_low = nImp_x * (i - mod(i,2))  !This is essentially giving the number of kpoints in each dimension.
        nSites_x_high = nImp_x * (i+1 + mod(i+1,2))
!        write(6,*) "nImp_x: ",nImp_x
!        write(6,*) "nSites_x_low: ",nSites_x_low
!        write(6,*) "nSites_x_high: ",nSites_x_high

        nBelowSitesChange = abs(nSites_x_low**2 - nSites)
        nAboveSitesChange = abs(nSites_x_high**2 - nSites)
        if(nBelowSitesChange.lt.nAboveSitesChange) then
            !Move number of sites down to nearest square number
            nSites_x = nSites_x_low
        else
            !Move number of sites down to nearest square number
            nSites_x = nSites_x_high   
        endif

        nSites_y = nSites_x   
        nSites = nSites_x**2

        if(nSites.ne.nSitesOrig) then
            write(6,"(A)") "Total number of sites changed to ensure square number and commensurate with impurity cluster"
        endif
        write(6,"(A,I9)") "Total number of sites: ",nSites

        !Now, set up ImpSites array, which gives the indices of the impurity cluster
        allocate(ImpSites(nImp))
        ImpSites = 0

        do i = 1,nImp_x
            do j = 1,nImp_x
                call LatCoordToSiteIndex_2DSquare(i,j,ImpSites(j+((i-1)*nImp_x)) 
            enddo
        enddo
        write(6,"(A)") "Impurity site indices defined as "
        do i = 1,nImp
            write(6,"(2I7)") i,ImpSites(i)
        enddo

        !TODO: Now, work out whether we want periodic or antiperiodic boundary conditions
!        tPeriodic = .true. 
!        tAntiPeriodic = .false. 

        if(tPeriodic) then
            write(6,"(A)") "PERIODIC boundary conditions chosen to ensure a closed shell fermi surface"
        else
            write(6,"(A)") "ANTIPERIODIC boundary conditions chosen to ensure a closed shell fermi surface"
        endif

        !Now, how many impurity "repeats" will there be once they are striped through the space?
        iImpRepeats = nint(real(nSites,dp)/real(nImp,dp))
        iImpRepeats_x = nint(real(nSites_x,dp)/real(nImp_x,dp))
        iImpRepeats_y = nint(real(nSites_y,dp)/real(nImp_y,dp))
        write(6,"(A,I6)") "Number of copied of correlation potential striped through space: ",iImpRepeats
        if((iImpRepeats_x*iImpRepeats_y).ne.iImpRepeats) call stop_all(t_r,'Error here')

        !TODO: Check that this has the same order as simply the site index
        allocate(StripedImpIndices(nImp,iImpRepeats))
        StripedImpIndices(:,:) = 0

!        write(6,*) "Impurity site indices: "
        do kx = 0,iImpRepeats_x-1
            dispx = iImpRepeats_x*LatticeVector(1,1) + iImpRepeats_x*LatticeVector(1,2)
            do ky = 0,iImpRepeats_y-1
                dispy = iImpRepeats_y*LatticeVector(2,1) + iImpRepeats_y*LatticeVector(2,2)
                do i = 1,nImp
                    !What is the displaced coordinate?
                    if(abs(indx-real(nint(indx),dp)).gt.1.0e-8_dp) call stop_all(t_r,'Error here x')
                    if(abs(indy-real(nint(indy),dp)).gt.1.0e-8_dp) call stop_all(t_r,'Error here y')
                    call FindDisplacedIndex_2DSquare(ImpSites(i),dispx,dispy,StripedImpIndices(i,(kx*iImpRepeats_x)+ky+1),PhaseChange,.false.)
                enddo
            enddo
        enddo

!        StartInd = 1    !This labels the first index of this impurity cluster
!        do k = 1,iImpRepeats
!            do i = 1,nImp_x
!                do j = 1,nImp_y
!                    call FindDisplacedIndex_2DSquare(StartInd,i-1,j-1,StripedImpIndices(((i-1)*nImp_x)+j,k),PhaseChange)
!
!!                    write(6,*) k,((i-1)*nImp_x)+j,StripedImpIndices(((i-1)*nImp_x)+j,k)
!
!                    !Check that we never have any boundary conditions since we should never have left the supercell
!                    if(PhaseChange.ne.1) call stop_all(t_r,'Should not be leaving supercell, so should not be changing phase!')
!                enddo
!            enddo
!
!            if(k.lt.iImpRepeats) then
!                !We want to move onto the next supercell index. 
!                !What are the lattice coordinates here?
!                call SiteIndexToLatCoord_2DSquare(StartInd,indx,indy)
!!                write(6,"(A,3I7)") "Current supercell index: ",StartInd,indx,indy
!                !Move indy nImp_y further down
!                indy = indy + nImp_y
!!                write(6,"(A,I7)") "Increasing y coordinate to: ",indy
!                !Map back in to the supercell
!                IndX_folded = py_mod(indx,nSites_X)
!                IndY_folded = py_mod(indy,nSites_Y)
!                if(IndX_folded.eq.0) IndX_folded = nSites_X
!                if(IndY_folded.eq.0) IndY_folded = nSites_Y
!                !Convert this back into a site index
!                call LatCoordToSiteIndex_2DSquare(IndX_folded,IndY_folded,Index_2)
!!                write(6,"(A,2I7)") "Mapped back into the supercell, such that index is now: ",indx_folded,indy_folded
!                if(Index_2.lt.StartInd) then
!                    !We have gone backwards. Move to the next column
!                    IndX_folded = IndX_folded + nImp_x
!!                    write(6,"(A,I7)") "We have actually mapped onto a previous supercell. Move to next column. New x coordinate: ",indx_folded
!                    if(IndX_folded.gt.nSites_x) call stop_all(t_r,'Error indexing')
!                    !Now, convert this back into a new lattice index
!                    call LatCoordToSiteIndex_2DSquare(IndX_folded,IndY_folded,Index_2)
!                    if(Index_2.lt.StartInd) call stop_all(t_r,'Error in indexing')
!                    if(Index_2.gt.nSites) call stop_all(t_r,'Error in indexing')
!                    StartInd = Index_2
!                else
!                    StartInd = Index_2
!                endif
!            endif
!        enddo


        !Check that the first repeat is the same as the indices we have already worked out for the main impurity
        do i = 1,nImp
            if(StripedImpIndices(i,1).ne.ImpSites(i)) call stop_all(t_r,'Something wrong here')
        enddo

    end subroutine Setup2DLattice_Square

    subroutine LatCoordToSiteIndex_2DSquare(Ind_X,Ind_Y,LatIndex)
        implicit none
        integer, intent(in) :: Ind_X,Ind_Y
        integer, intent(out) :: LatIndex
        character(len=*), parameter :: t_r='LatCoordToSiteIndex_2DSquare'

        if((Ind_X.lt.1).or.(Ind_X.gt.nSites_x)) call stop_all(t_r,'This routine must take coordinates inside the supercell')
        if((Ind_Y.lt.1).or.(Ind_Y.gt.nSites_y)) call stop_all(t_r,'This routine must take coordinates inside the supercell')

        LatIndex = (Ind_X-1)*nSites_y + Ind_Y

        if((LatIndex.lt.1).or.(LatIndex.gt.nSites)) call stop_all(t_r,'This routine must take an index inside the supercell')

    end subroutine LatCoordToSiteIndex_2DSquare

    !Find the lattice coordinate from the site index
    !The LatIndex has to be between 1 and nSites
    !The indexing is row major like a coordinate, i.e (1,1) (2,1)
    !                                                 (2,1) (2,2)
    subroutine SiteIndexToLatCoord_2DSquare(LatIndex,Ind_X,Ind_Y)
        implicit none
        integer, intent(in) :: LatIndex
        integer, intent(out) :: Ind_X,Ind_Y
        character(len=*), parameter :: t_r='SiteIndexToLatCoord_2DSquare'

        if((LatIndex.lt.1).or.(LatIndex.gt.nSites)) call stop_all(t_r,'This routine must take an index inside the supercell')

        Ind_Y = mod(LatIndex,nSites_y)
        if(Ind_Y.eq.0) Ind_Y = nSites_y
!        write(6,*) "IndY: ",Ind_Y
!        write(6,*) "Ind_X: ",((LatIndex - Ind_Y) / nSites_x) + 1,LatIndex - Ind_Y
        Ind_X = ((LatIndex - Ind_Y) / nSites_x) + 1

    end subroutine SiteIndexToLatCoord_2DSquare

    !Find the index of a lattice site displaced by a vector given by (DeltaX,DeltaY), and also calculate any phase change if necessary
    !In this scheme, the lattice is indexed by increasing in the y direction first (column major)
    subroutine FindDisplacedIndex_2DSquare(LatIndIn,DeltaX,DeltaY,LatIndOut,PhaseChange,tAllowWrap)
        implicit none
        integer, intent(in) :: LatIndIn,DeltaX,DeltaY
        integer, intent(out) :: LatIndOut
        integer, intent(out) :: PhaseChange
        logical, intent(in), optional :: tAllowWrap
        !local
        integer :: Ind_X,Ind_Y,IndX_folded,IndY_folded,flips_x,flips_y
        logical :: tAllowWrap_
        character(len=*), parameter :: t_r='FindDisplacedIndex_2DSquare'

        if(present(tAllowWrap)) then
            tAllowWrap_ = tAllowWrap
        else
            tAllowWrap_ = .true.
        endif

!        write(6,*) "Calling Find displaced index. Original ind: ",LatIndIn
!        write(6,*) "Displacement vector: ",DeltaX,DeltaY
        LatIndOut = LatIndIn
        PhaseChange = 1

        call SiteIndexToLatCoord_2DSquare(LatIndOut,Ind_X,Ind_Y)

!        write(6,*) "Original index coordinate: ",Ind_X,Ind_Y

        if(CellShape.eq.2) then
            Ind_X = Ind_X + DeltaX
            Ind_Y = Ind_Y + DeltaY

            !Now, we may need to map back into the supercell
            IndX_folded = py_mod(Ind_X,nSites_X)
            IndY_folded = py_mod(Ind_Y,nSites_Y)
            if(IndX_folded.eq.0) IndX_folded = nSites_X
            if(IndY_folded.eq.0) IndY_folded = nSites_Y

            if(tAntiPeriodic) then
                !Has this wrapped an even or odd number of times?
                if(Ind_X.ge.1) then
                    flips_x = abs(mod(int(real(Ind_X-1,dp)/real(nSites_x,dp)),2))
                else
                    flips_x = abs(mod(int(real(Ind_X-nSites_x,dp)/real(nSites_x,dp)),2))
                endif
                if(Ind_Y.ge.1) then
                    flips_y = abs(mod(int(real(Ind_Y-1,dp)/real(nSites_y,dp)),2))
                else
                    flips_y = abs(mod(int(real(Ind_Y-nSites_y,dp)/real(nSites_y,dp)),2))
                endif
                if(flips_x.eq.0) then
                    flips_x = 1
                else
                    flips_x = -1
                endif
                if(flips_y.eq.0) then
                    flips_y = 1
                else
                    flips_y = -1
                endif
            else
                !With PBCs, it doesn't actually matter how many flips we use, since we will never change phase
                flips_x = 1
                flips_y = 1
            endif

            Ind_X = IndX_folded
            Ind_Y = IndY_folded

            !What is the overall phase (assume same BCs on x and y directions)
            PhaseChange = flips_x * flips_y

        elseif(CellShape.eq.3)

            !First, translate one by one in the x direction
            do i = 1,abs(DeltaX)
                Ind_X = Ind_X + sign(1,DeltaX)
                Ind_Y = Ind_Y - sign(1,DeltaX)

                if(Ind_Y.eq.0) then
                    Ind_Y = Ind_Y + nSites_y
                    yFlip = yFlip + 1
                elseif(Ind_Y.eq.(nSites_y+1)) then
                    Ind_Y = Ind_Y - nSites_y
                    yFlip = yFlip + 1
                endif
                if(Ind_X.eq.0) then
                    Ind_X = Ind_X + nSites_x
                    xFlip = xFlip + 1
                elseif(Ind_X.eq.(nSites_x+1)) then
                    Ind_X = Ind_X - nSites_x
                    xFlip = xFlip + 1
                endif
            enddo
            do i = 1,abs(DeltaY)
                Ind_Y = Ind_Y + sign(1,DeltaY)
                if(Ind_Y.eq.(nSites_y+1)) then
                    Ind_Y = Ind_Y - nSites_y
                    yFlip = yFlip + 1
                elseif(Ind_Y.eq.0) then
                    Ind_Y = Ind_Y + nSites_y
                    yFlip = yFlip + 1
                endif
            enddo
            if(.not.tAllowFlips_) then
                if((xFlip+yFlip).ne.0) call stop_all(t_r,'Moved out of supercell')
            endif

            if(tPeriodic) then
                PhaseChange = 1
            else
                !APBCs (same in both directions
                PhaseChange = mod(yFlip+xFlip,2)
                if(PhaseChange.eq.0) then
                    PhaseChange = 1
                else
                    PhaseChange = -1
                endif
            endif
        endif

        !Now, translate back into a lattice index
        call LatCoordToSiteIndex_2DSquare(Ind_X,Ind_Y,LatIndOut)

!        write(6,*) "Final lattice index: ",LatIndOut
!        write(6,*) "phase: ",PhaseChange

    end subroutine FindDisplacedIndex_2DSquare
    
    !Setup the lattice for the tilted unit cell
    subroutine Setup2DLattice_Tilt()
        implicit none
        real(dp) :: dWidth
        integer :: TDLat_Width,x,y,dx,dy,ci,cj,site_imp
        integer :: i,j,k
        character(len=*), parameter :: t_r='Setup2DLattice'

        !Work out an appropriate width, and an actual nSites
        dWidth = sqrt(nSites*2.0_dp)
        TDLat_Width = 2 * nint(dWidth/2.0_dp)

        !There are two lattices, an x, y lattice, and an i, j lattice. These have had their axes rotated by 45 degrees.
        !2DLat_Ni is the width of each lattice (there are two interlocking in the (i,j) representation
        TDLat_Ni = TDLat_Width / 2
        !2DLat_Nj is the width of the lattice in the (x,y) representation
        TDLat_Nj = TDLat_Width

        !Actual nSites. 
        nSites = TDLat_Ni * TDLat_Nj

        write(6,*) "Updated number of sites in the 2D hubbard model will be: ",nSites
        
        if(mod(TDLat_Width/2,2).eq.0) then
            !Use anti-periodic boundary conditions
            !HF ground state only unique if Width=2*odd_number (direct PBC)
            !or if Width=2*even_number (anti-PBC)
            tAntiperiodic = .true.
            tPeriodic = .false.
        else
            tAntiperiodic = .false.
            tPeriodic = .true.
        endif
        if(tPeriodic) then
            write(6,*) "Periodic boundary conditions now in use"
        else
            write(6,*) "Anti-Periodic boundary conditions now in use"
        endif

        !Now to set up the impurity
        if(nImp.eq.1) then
            nImp_x = 1
            nImp_y = 1
        elseif(nImp.eq.2) then
            nImp_x = 1
            nImp_y = 2
        elseif(nImp.eq.4) then
            nImp_x = 2
            nImp_y = 2
        else
            call stop_all(t_r,'Cannot deal with impurities > 4')
        endif

        !Find the x,y coordinates for the middle of the array. This will be used to
        !define the corner site of the impurity
        call ij2xy(TDLat_Ni/2,TDLat_Nj/2,x,y)

        !Setup the impurity space, and how to tile the impurity through the space.
        !This creates the matrices TD_Imp_Lat and TD_Imp_Phase
        !If the correlation potential is a matrix of nImp x nImp, then TD_Imp_Lat
        !gives the index of that correlation potential which corresponds to the tiled
        !correlation potential through the space.
        !(TD_Imp_Lat(site,site)-1)/nImp + 1 gives the impurity index of that site. 

        call MakeVLocIndices()

        !Create the impurity cluster
        allocate(ImpSites(nImp))
        ImpSites = 0
        do dx = 0,nImp_x-1
            do dy = 0,nImp_y-1
                call xy2ij(x+dx,y+dy,ci,cj)
                !Remember - our sites are 1 indexed
                site_imp = ci + TDLat_Ni*cj + 1     !No need to take mods. We certainly shouldn't exceed the bounds of the ij lattice
                !write(6,*) "***",dx,dy,site_imp,(TD_Imp_Lat(site_imp,site_imp)),((TD_Imp_Lat(site_imp,site_imp)-1)/nImp) + 1
                ImpSites(((TD_Imp_Lat(site_imp,site_imp)-1)/nImp) + 1) = site_imp
            enddo
        enddo

        write(6,*) "Impurity sites defined as: ",ImpSites(:)

        !We now want to define a mapping, from the standard site indexing to an impurity ordering of the sites, such that
        ! Perm_indir(site) = Imp_site_index
        ! Perm_dir(Imp_site_index) = site       The first index maps you onto the original impurity sites
        ! Therefore, Perm_indir(site)%nImp = impurity site it maps to which repeat of the impurity you are on

        !In the impurity ordering (indirect), the impurity sites are first, followed by the repetitions of the striped
        !impurity space.
        !The direct space is the normal lattice ordering
        allocate(Perm_indir(nSites))
        allocate(Perm_dir(nSites))
        Perm_indir(:) = 0
        Perm_dir(:) = 0
        Perm_dir(1:nImp) = ImpSites(:)
        k = nImp+1
        loop: do i = 1,nSites
            do j = 1,nImp 
                if(i.eq.ImpSites(j)) then
                    cycle loop
                endif
            enddo
            Perm_dir(k) = i
            k = k + 1
        enddo loop
        if(k.ne.nSites+1) call stop_all(t_r,"Error here")

        do i=1,nSites
            Perm_indir(Perm_dir(i)) = i
        enddo
        !write(6,*) "Perm_dir: ",Perm_dir(:)
        !write(6,*) "Perm_indir: ",Perm_indir(:)

    end subroutine Setup2DLattice_Tilt

    !Find the indices for any diagonal operator, and the matrices for ...
    subroutine MakeVLocIndices()
        implicit none
        integer :: site,iv,i,j,xb,yb,iv_target,dx,dy,phase,ij,Conn_x,Conn_y,k
        integer, allocatable :: Imp_Connections(:,:)
        character(len=*), parameter :: t_r='MakeVLocIndices'

        if(LatticeDim.ne.2) call stop_all(t_r,'Should only be in here with 2D lattice')
        
        allocate(TD_Imp_Lat(nSites,nSites))
        allocate(TD_Imp_Phase(nSites,nSites))
        allocate(Imp_Connections(2,nImp))

        TD_Imp_Lat(:,:) = 0
        TD_Imp_Phase(:,:) = 0

!        write(6,*) "nImp_x: ",nImp_X
!        write(6,*) "nImp_y: ",nImp_y
!        write(6,*) "Li: ",TDLat_Ni
!        write(6,*) "Lj: ",TDLat_Nj
        do site = 0,nSites-1

            !Find ij indices
            call site2ij(site,i,j)
!            write(6,*) "site: ",site,"ij: ",i,j
            !Map i, j indices to x, y
            call ij2xy(i,j,xb,yb)
            !Which impurity copy are they in
            !All sites with the same (xb,yb) are in the same impurity copy
            xb = xb - py_mod(xb,nImp_x)
            yb = yb - py_mod(yb,nImp_y)

!            write(6,*) "site: ",site,"xb, yb: ",xb,yb

            Imp_Connections(:,:) = 0
            k = 1
            do dx = 0,nImp_x-1
                do dy = 0,nImp_y-1
                    if((mod(nImp_x,2).eq.1).and.(mod((xb/nImp_x),2).eq.1)) then
                        !Odd number of lattice points in the x direction of the impurity plaquette *and*
                        !this is an odd plaquette number
                        Imp_Connections(1,k) = nImp_X - dx - 1
                    else
                        Imp_Connections(1,k) = dx
                    endif
                    if((mod(nImp_y,2).eq.1).and.(mod((yb/nImp_y),2).eq.1)) then
                        Imp_Connections(2,k) = nImp_Y - dy - 1
                    else
                        Imp_Connections(2,k) = dy
                    endif
                    k = k+1
                enddo
            enddo

            iv = 0
            iv_target = -1
                    
!            write(6,*) "Delta: ",Imp_Connections(:,:)
            !Now run through the impurity plaquette of the site (xb,yb)
            !Find which impurity site this site maps on to (iv_target). 
            !All sites should map onto an impurity site
            k = 1
            do dx = 0,nImp_x-1
                do dy = 0,nImp_y-1

                    Conn_x = Imp_Connections(1,k)
                    Conn_y = Imp_Connections(2,k)

                    !write(6,*) "x, y: ",xb+Conn_x,yb+Conn_y
                    call xy2ij(xb+Conn_x,yb+Conn_y,i,j)
                    !write(6,*) "i, j: ",i,j,py_mod(i,TDLat_Ni),TDLat_Ni*py_mod(j,TDLat_Nj), &
                        !py_mod(i,TDLat_Ni) + TDLat_Ni*py_mod(j,TDLat_Nj)
                    call ij2site(i,j,ij)
!                    write(6,*) "checking: ",ij,i,j,site
                    if(ij.eq.site) then
                        !Does it map on to this impurity site?
                        if(iv_target.ne.-1) call stop_all(t_r,'iv target already found?')
                        iv_target = iv
!                        write(6,*) "Setting iv_target: ",iv_target
                    endif
                    iv = iv + 1
                    k = k + 1
                enddo
            enddo
            if(iv_target.eq.-1) call stop_all(t_r,'iv_target not set')

            iv = 0
            k = 1
            do dx = 0,nImp_x-1
                do dy = 0,nImp_y-1

                    Conn_x = Imp_Connections(1,k)
                    Conn_y = Imp_Connections(2,k)

                    call xy2ij(xb+Conn_x,yb+Conn_y,i,j)
                    call ij2site(i,j,ij)
                    
                    phase = 1
                    if(tAntiPeriodic.and.((i.lt.0).or.(i.ge.TDLat_Ni))) phase = -phase
                    if(tAntiPeriodic.and.((j.lt.0).or.(j.ge.TDLat_Nj))) phase = -phase
                    if(TD_Imp_Phase(ij+1,site+1).ne.0) then
                        call stop_all(t_r,'We should not have written here yet')
                    endif
                    TD_Imp_Lat(ij+1,site+1) = iv + nImp*iv_target + 1
                    TD_Imp_Phase(ij+1,site+1) = phase

                    iv = iv + 1
                    k = k + 1
                enddo
            enddo

        enddo
!        write(6,*) "impurity lattice matrix: "
!        do i = 1,nSites
!            do j = 1,nSites
!                write(6,"(I4)", advance='no') TD_Imp_Lat(i,j)
!            enddo
!            write(6,*) 
!        enddo
!        write(6,*) "***"
!        write(6,*) "Impurity phase mastrix: "
!        do i = 1,nSites
!            do j = 1,nSites
!                write(6,"(I4)", advance='no') TD_Imp_Phase(i,j)
!            enddo
!            write(6,*) 
!        enddo
        deallocate(Imp_Connections)
    end subroutine MakeVLocIndices
    

end module Lattices
