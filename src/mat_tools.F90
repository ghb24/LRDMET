!tools for matrix construction, manipulation and HF solver for DMET
module mat_tools
    use const
    use errors, only: stop_all, warning
    use globals
    implicit none

    interface DiagOneEOp
        module procedure DiagOneEOp_r
        module procedure DiagOneEOp_z
    end interface

    interface CheckRealSpaceTransInvar
        module procedure CheckRealSpaceTransInvar_r
        module procedure CheckRealSpaceTransInvar_z
    end interface

    contains

    !Setup arrays and indices needed to define a 2D square, non-tilted lattice
    subroutine Setup2DLattice_Square()
        implicit none
        integer :: nSitesOrig,nSites_x_low,nSites_x_high,nBelowSitesChange
        integer :: i,j,k,nAboveSitesChange,PhaseChange,StartInd
        integer :: indx,indy,indx_folded,indy_folded,index_2
        character(len=*), parameter :: t_r='Setup2DLattice_Square'

        write(6,"(A)") "Setting up a non-tilted square lattice..."

        if(LatticeDim.ne.2) call stop_all(t_r,'Should only be in here with 2D lattice')

        nImp_x = nint(sqrt(real(nImp,dp)))
        nImp_y = nImp_x
        if((sqrt(real(nImp,dp))-nint(sqrt(real(nImp,dp)))).gt.1.0e-8_dp) then
            !Constrain impurity cluster to have the same geometry unit cell as the 
            !entire lattice.
            call stop_all(t_r,'Number of impurity sites is not a square number.')
        endif

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

        nSites_y = nSites_x     !Non-tilted square lattice
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
                ImpSites(j+((i-1)*nImp_x)) = j + ((i-1)*nSites_x)
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
        write(6,"(A,I6)") "Number of copied of correlation potential striped through space: ",iImpRepeats

        allocate(StripedImpIndices(nImp,iImpRepeats))
        StripedImpIndices(:,:) = 0

!        write(6,*) "Impurity site indices: "
        StartInd = 1    !This labels the first index of this impurity cluster
        do k = 1,iImpRepeats
            do i = 1,nImp_x
                do j = 1,nImp_y
                    call FindDisplacedIndex_2DSquare(StartInd,i-1,j-1,StripedImpIndices(((i-1)*nImp_x)+j,k),PhaseChange)

!                    write(6,*) k,((i-1)*nImp_x)+j,StripedImpIndices(((i-1)*nImp_x)+j,k)

                    !Check that we never have any boundary conditions since we should never have left the supercell
                    if(PhaseChange.ne.1) call stop_all(t_r,'Should not be leaving supercell, so should not be changing phase!')
                enddo
            enddo

            if(k.lt.iImpRepeats) then
                !We want to move onto the next supercell index. 
                !What are the lattice coordinates here?
                call SiteIndexToLatCoord_2DSquare(StartInd,indx,indy)
!                write(6,"(A,3I7)") "Current supercell index: ",StartInd,indx,indy
                !Move indy nImp_y further down
                indy = indy + nImp_y
!                write(6,"(A,I7)") "Increasing y coordinate to: ",indy
                !Map back in to the supercell
                IndX_folded = py_mod(indx,nSites_X)
                IndY_folded = py_mod(indy,nSites_Y)
                if(IndX_folded.eq.0) IndX_folded = nSites_X
                if(IndY_folded.eq.0) IndY_folded = nSites_Y
                !Convert this back into a site index
                call LatCoordToSiteIndex_2DSquare(IndX_folded,IndY_folded,Index_2)
!                write(6,"(A,2I7)") "Mapped back into the supercell, such that index is now: ",indx_folded,indy_folded
                if(Index_2.lt.StartInd) then
                    !We have gone backwards. Move to the next column
                    IndX_folded = IndX_folded + nImp_x
!                    write(6,"(A,I7)") "We have actually mapped onto a previous supercell. Move to next column. New x coordinate: ",indx_folded
                    if(IndX_folded.gt.nSites_x) call stop_all(t_r,'Error indexing')
                    !Now, convert this back into a new lattice index
                    call LatCoordToSiteIndex_2DSquare(IndX_folded,IndY_folded,Index_2)
                    if(Index_2.lt.StartInd) call stop_all(t_r,'Error in indexing')
                    if(Index_2.gt.nSites) call stop_all(t_r,'Error in indexing')
                    StartInd = Index_2
                else
                    StartInd = Index_2
                endif
            endif
        enddo


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
    subroutine FindDisplacedIndex_2DSquare(LatIndIn,DeltaX,DeltaY,LatIndOut,PhaseChange)
        implicit none
        integer, intent(in) :: LatIndIn,DeltaX,DeltaY
        integer, intent(out) :: LatIndOut
        integer, intent(out) :: PhaseChange
        !local
        integer :: Ind_X,Ind_Y,IndX_folded,IndY_folded,flips_x,flips_y

!        write(6,*) "Calling Find displaced index. Original ind: ",LatIndIn
!        write(6,*) "Displacement vector: ",DeltaX,DeltaY
        LatIndOut = LatIndIn
        PhaseChange = 1

        call SiteIndexToLatCoord_2DSquare(LatIndOut,Ind_X,Ind_Y)

!        write(6,*) "Original index coordinate: ",Ind_X,Ind_Y

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
    
    !Make mean-field real-space hubbard matrix
    subroutine make_hop_mat()
        use DetTools, only: tospat
        use utils, only: get_free_unit
        implicit none
        integer :: i,j,k,x,y,li,lj,ij_link,site,dx,dy,iunit,iphase
        real(dp) :: phase,t
        real(dp), allocatable :: temp(:)
        logical :: exists
        character(len=*), parameter :: t_r='make_hop_mat'

        h0(:,:) = zero  
        t = -1.0_dp
        if(tReadSystem) then
            !Read in the hopping matrix
            inquire(file='CoreHam.dat',exist=exists)
            if(.not.exists) call stop_all(t_r,'Core Hamiltonian file cannot be found')
            iunit = get_free_unit()
            open(iunit,file='CoreHam.dat',status='old',action='read')
            if(tUHF) then
                allocate(temp(nSites*2))
                !Ordered as alpha, beta
                do i = 1,nSites*2
                    read(iunit,*) temp(:)

                    !Now split into respective components
                    do j = 1,nSites*2
                        if(mod(i,2).eq.1) then
                            !i is alpha spin-orbital

                            if((mod(j,2).eq.0).and.(abs(temp(j)).gt.1.0e-8_dp)) then
                                call stop_all(t_r,'Coupling between different spin types in core hamiltonian??')
                            elseif(mod(j,2).eq.1) then
                                !j is also alpha
                                h0(tospat(j),tospat(i)) = temp(j)
                            endif

                        else
                            !i is beta spin-orbital
                            if((mod(j,2).eq.1).and.(abs(temp(j)).gt.1.0e-8_dp)) then
                                call stop_all(t_r,'Coupling between different spin types in core hamiltonian??')
                            elseif(mod(j,2).eq.0) then
                                !j is also beta 
                                h0_b(tospat(j),tospat(i)) = temp(j)
                            endif
                        endif
                    enddo
                enddo
            else
                do i=1,nSites
                    read(iunit,*) h0(:,i)
                enddo
            endif
            close(iunit)
            !Also read in lattice matrix with correlation potential
            inquire(file='FinalFock.dat',exist=exists)
            if(.not.exists) call stop_all(t_r,'Lattice hamiltonian file cannot be found')
            iunit = get_free_unit()
            open(iunit,file='FinalFock.dat',status='old',action='read')
            if(tUHF) then
                do i=1,nSites*2
                    read(iunit,*) temp(:)
                    !Now split into respective components
                    do j = 1,nSites*2
                        if(mod(i,2).eq.1) then
                            !i is alpha spin-orbital

                            if((mod(j,2).eq.0).and.(abs(temp(j)).gt.1.0e-8_dp)) then
                                call stop_all(t_r,'Coupling between different spin types in one-e hamiltonian??')
                            elseif(mod(j,2).eq.1) then
                                !j is also alpha
                                h0v(tospat(j),tospat(i)) = temp(j)
                            endif

                        else
                            !i is beta spin-orbital
                            if((mod(j,2).eq.1).and.(abs(temp(j)).gt.1.0e-8_dp)) then
                                call stop_all(t_r,'Coupling between different spin types in one-e hamiltonian??')
                            elseif(mod(j,2).eq.0) then
                                !j is also beta 
                                h0v_b(tospat(j),tospat(i)) = temp(j)
                            endif
                        endif
                    enddo
                enddo

                deallocate(temp)
            else
                do i=1,nSites
                    read(iunit,*) h0v(:,i)
                enddo
            endif
            close(iunit)

        elseif(tSingFiss) then
            !Read in parameters for a two-band lattice
            call Read_TwoBandLatt()

        elseif(LatticeDim.eq.1) then
            !Tridiagonal matrix
            do i=1,nSites-1
                h0(i,i+1) = t
                h0(i+1,i) = t
            enddo

            if(tPeriodic) then
                h0(1,nSites) = t
                h0(nSites,1) = t
            elseif(tAntiPeriodic) then
                h0(1,nSites) = -t
                h0(nSites,1) = -t
            endif

            if(tChemPot.and.tAnderson) then
                h0(1,1) = h0(1,1) - (U/2.0_dp)
            endif

        elseif(LatticeDim.eq.2) then

            if(tTiltedLattice) then
                !Tilted lattices. No k-space coded up

                do site = 0,nSites-1
                    call site2ij(site,i,j)
                    !Now for coordinates on the original tilted lattice. This is the one where distances are not distorted
                    call ij2xy(i,j,x,y)
                    !write(6,*) "site: ",site," i,j: ",i,j,"x,y: ",x,y
                    !call flush(6)

                    do k = 1,4
                        !Move in all four possible directions
                        if(k.eq.1) then
                            dx = 1
                            dy = 0
                        elseif(k.eq.2) then
                            dx = 0
                            dy = 1
                        elseif(k.eq.3) then
                            dx = 0
                            dy = -1
                        elseif(k.eq.4) then
                            dx = -1
                            dy = 0
                        endif
                        !Move in the appropriate direction, and see where we wrap around to, by transforming into the i,j representation, where boundary conditions are easier
                        call xy2ij(x+dx,y+dy,li,lj)
                        phase = 1.0_dp
                        if(tAntiPeriodic) then
                            !We have to multiply the phase by a factor of -1 for every boundary condition we wrap around
                            if((li.lt.0).or.(li.ge.TDLat_Ni)) phase = -phase
                            if((lj.lt.0).or.(lj.ge.TDLat_Nj)) phase = -phase
                        endif

                        !Take our position modulo the axes
                        li = py_mod(li,TDLat_Ni) 
                        lj = py_mod(lj,TDLat_Nj) 
                        ij_link = li + TDLat_Ni*lj  !Now find the site number of the connecing lattice
                        h0(site+1,ij_link+1) = t * phase
                        h0(ij_link+1,site+1) = t * phase
                    enddo
                enddo

                !Change the ordering of the hopping hamiltonian, such that the impurity sites are defined firust, and then repeated in the order of the impurty tiling
!                write(6,*) "Hopping matrix in lattice ordering: "
!                do i = 1,nSites
!                    do j = 1,nSites
!                        write(6,"(F7.3)",advance='no') h0(j,i)
!                    enddo
!                    write(6,*) ""
!                enddo
                call Mat_to_imp_order(h0)
            else
                !Non-tilted lattice - square
                do i=1,nSites
                    !Run through all sites
                    do k = 1,4
                        !Move in all four possible directions
                        if(k.eq.1) then
                            dx = 1
                            dy = 0
                        elseif(k.eq.2) then
                            dx = 0
                            dy = 1
                        elseif(k.eq.3) then
                            dx = 0
                            dy = -1
                        elseif(k.eq.4) then
                            dx = -1
                            dy = 0
                        endif

                        !Find the index of the displaced lattice site
                        call FindDisplacedIndex_2DSquare(i,dx,dy,j,iphase)

                        if((h0(i,j).ne.zero).and.(abs(h0(i,j)-(t*iphase)).gt.1.0e-8_dp)) then
                            call stop_all(t_r,'Error in constructing h0 matrix')
                        endif
                        if((h0(j,i).ne.zero).and.(abs(h0(j,i)-(t*iphase)).gt.1.0e-8_dp)) then
                            call stop_all(t_r,'Error in constructing h0 matrix')
                        endif

                        h0(i,j) = t*iphase
                        h0(j,i) = t*iphase
                    enddo
                enddo

            endif
        else
            call stop_all(t_r,'Higher dimensional models not coded up')
        endif

    end subroutine make_hop_mat
            
    subroutine Read_TwoBandLatt()
        use utils, only: get_free_unit
        implicit none
        logical :: exists
        integer :: iunit,irow,icol,ios
        real(dp) :: val
        character(len=*), parameter :: t_r='Read_TwoBandLatt'
        
        write(6,"(A)") "Reading in global hopping integrals..."
    
        inquire(file='Tfort.dat',exist=exists)
        if(.not.exists) call stop_all(t_r,'T matrix file cannot be found')
        iunit = get_free_unit()
        open(iunit,file='Tfort.dat',status='old',action='read')
        icol = 1
        irow = 1
        do while(.true.)
            !i  = irow +  4*(icellrow-1) + 108*(icol-1) + 108*4*(icellcol-1)
            read(iunit,*,iostat=ios) val 
            if(ios.gt.0) call stop_all(t_r,'Error reading integrals')
            if(ios.lt.0) exit   !EOF

            h0(irow,icol) = val

            irow = irow + 1
            if(irow.gt.nSites) then
                irow = 1
                icol = icol + 1
            endif

        enddo
        close(iunit)

        !Now read in the local coulomb and exchange integrals
        if(.not.allocated(J_Ints)) allocate(J_Ints(nImp,nImp))
        if(.not.allocated(X_Ints)) allocate(X_Ints(nImp,nImp))
        J_Ints(:,:) = zero
        X_Ints(:,:) = zero

        write(6,"(A)") "Reading in local coulomb integrals..."

        inquire(file='Ufort.dat',exist=exists)
        if(.not.exists) call stop_all(t_r,'U matrix file cannot be found')
        open(iunit,file='Ufort.dat',status='old',action='read')
        icol = 1
        irow = 1
        do while(.true.)
            read(iunit,*,iostat=ios) val
            if(ios.gt.0) call stop_all(t_r,'Error reading integrals')
            if(ios.lt.0) exit   !EOF

            if((icol.le.nImp).and.(irow.le.nImp)) then
                !Store the integral
                J_Ints(irow,icol) = val
            endif
            irow = irow + 1
            if(irow.gt.nSites) then
                irow = 1
                icol = icol + 1
            endif
        enddo
        close(iunit)
        
        write(6,"(A)") "Reading in local exchange integrals..."

        inquire(file='Xfort.dat',exist=exists)
        if(.not.exists) call stop_all(t_r,'X matrix file cannot be found')
        open(iunit,file='Xfort.dat',status='old',action='read')
        icol = 1
        irow = 1
        do while(.true.)
            read(iunit,*,iostat=ios) val
            if(ios.gt.0) call stop_all(t_r,'Error reading integrals')
            if(ios.lt.0) exit   !EOF

            if((icol.le.nImp).and.(irow.le.nImp)) then
                !Store the integral
                X_Ints(irow,icol) = val
            endif
            irow = irow + 1
            if(irow.gt.nSites) then
                irow = 1
                icol = icol + 1
            endif
        enddo
        close(iunit)

    end subroutine Read_TwoBandLatt

    !Permute the ordering of a real matrix in the lattice ordering, such that it turns into a matrix
    !in the impurity ordering, such that the impurities are defined first.
    subroutine Mat_to_imp_order(h)
        implicit none
        real(dp) , intent(inout) :: h(nSites,nSites)
        real(dp) , allocatable :: temp(:,:)
        integer :: i
        
        allocate(temp(nSites,nSites))
        temp(:,:) = zero

        !Permute the columns
        do i = 1,nSites
            temp(:,i) = h(:,Perm_dir(i))
        enddo

        h(:,:) = zero
        !Permute the rows, and overwrite original matrix
        do i = 1,nSites
            h(i,:) = temp(Perm_dir(i),:)
        enddo
        deallocate(temp)

    end subroutine Mat_to_imp_order

    subroutine Mat_to_lattice_order(h)
        implicit none
        real(dp), intent(inout) :: h(nSites,nSites)
        real(dp) , allocatable :: temp(:,:)
        integer :: i

        allocate(temp(nSites,nSites))
        temp(:,:) = zero
        do i = 1,nSites
            temp(:,i) = h(:,Perm_indir(i))
        enddo
        h(:,:) = 0.0_dp
        do i = 1,nSites
            h(i,:) = temp(Perm_indir(i),:)
        enddo
        deallocate(temp)

    end subroutine Mat_to_lattice_order
    
    subroutine Mat_to_imp_order_comp(h)
        implicit none
        complex(dp) , intent(inout) :: h(nSites,nSites)
        complex(dp) , allocatable :: temp(:,:)
        integer :: i
        
        allocate(temp(nSites,nSites))
        temp(:,:) = zzero
        !Permute the columns
        do i = 1,nSites
            temp(:,i) = h(:,Perm_dir(i))
        enddo

        h(:,:) = zzero
        !Permute the rows, and overwrite original matrix
        do i = 1,nSites
            h(i,:) = temp(Perm_dir(i),:)
        enddo
        deallocate(temp)

    end subroutine Mat_to_imp_order_comp

    subroutine Mat_to_lattice_order_comp(h)
        implicit none
        complex(dp), intent(inout) :: h(nSites,nSites)
        complex(dp) , allocatable :: temp(:,:)
        integer :: i

        allocate(temp(nSites,nSites))
        temp(:,:) = zzero
        do i = 1,nSites
            temp(:,i) = h(:,Perm_indir(i))
        enddo
        h(:,:) = zzero
        do i = 1,nSites
            h(i,:) = temp(Perm_indir(i),:)
        enddo
        deallocate(temp)

    end subroutine Mat_to_lattice_order_comp

    !Care needed with these functions. Python 'floors' the result, whereas fortran will truncate towards zero. Therefore
    !the results will be the same when x .ge. y, but not otherwise.
    !Also, the modulo function in python always returns the same sign as the divisor, and also when dividing a negative number also truncates to -inf.
    !See http://stackoverflow.com/questions/1907565/c-python-different-behaviour-of-the-modulo-operation

    !Function to return the same value as the python modulo function
    !Will be the same if both arguments are positive
    elemental function py_mod(n,M) result(res)
        implicit none
        integer, intent(in) :: n,M
        integer :: res
        res = mod(mod(n,M) + M,M)
    end function py_mod

    pure subroutine xy2ij(x,y,i,j)
        implicit none
        integer, intent(in) :: x,y
        integer, intent(out) :: i,j

        if(x.ge.y) then
            i = (x-y)/2
        else
            if(mod(x-y,2).eq.0) then
                i = (x-y)/2
            else
                i = ( (x-y)/2 ) - 1
            endif
        endif
        j = x+y

    end subroutine xy2ij
                
    pure subroutine ij2xy(i,j,x,y)
        implicit none
        integer, intent(in) :: i,j
        integer, intent(out) :: x,y

        if(j.ge.0) then
            x = i + (j/2) + py_mod(j,2)
            y = (j/2) - i
        else
            if(mod(j,2).eq.0) then
                x = i + (j/2)
                y = (j/2) - i
            else
                x = i + (j/2) - 1 + py_mod(j,2)
                y = (j/2) - i
            endif
        endif

    end subroutine ij2xy

    pure subroutine ij2site(i,j,site)
        implicit none
        integer , intent(in) :: i,j
        integer , intent(out) :: site

        site = py_mod(i,TDLat_Ni) + TDLat_Ni*py_mod(j,TDLat_Nj)
    end subroutine ij2site
    pure subroutine site2ij(site,i,j)
        implicit none
        integer , intent(in) :: site
        integer , intent(out) :: i,j

        i = py_mod(site,TDLat_Ni)
        j = site/TDLat_Ni
        if(site.lt.0) then 
            if(mod(site,TDLat_Ni).ne.0) then
                j = site/TDLat_Ni - 1
            endif
        endif
    end subroutine site2ij
    
    !Diagonalizes the core hamiltonian, and stores the allowed CS occupations in terms of
    !number of fully occupied sites.
    subroutine find_occs()
        implicit none
        real(dp), allocatable :: Work(:),h0Eigenvecs(:,:),h0Eigenvals(:)
        integer :: lWork,info,i,j,k
        integer :: iImp,nHoppingsImp,nHoppingsEnv
        character(len=*), parameter :: t_r="find_occs"

!        call writematrix(h0,'Core hamil',.false.)

        if(allocated(allowed_occs)) deallocate(allowed_occs)
        if(tHalfFill) then
            N_occs = 1
            allocate(allowed_occs(N_occs))
            allowed_occs(1) = nSites/2 
        elseif(nElecFill.ne.0) then
            !We have a specific inputted filling fraction
            N_occs = 1
            allocate(allowed_occs(N_occs))
            allowed_occs(1) = nElecFill/2
        else
            if(tUHF) call stop_all(t_r,'Cannot run through all occupations with UHF yet - fix me')
            allocate(h0Eigenvecs(nSites,nSites))
            allocate(h0Eigenvals(nSites))
            h0Eigenvecs = 0.0_dp
            !First, diagonalize one-body hamiltonian
            h0Eigenvecs(:,:) = h0(:,:)
            h0Eigenvals(:) = 0.0_dp
            allocate(Work(1))
            lWork=-1
            info=0
            call dsyev('N','U',nSites,h0Eigenvecs,nSites,h0Eigenvals,Work,lWork,info)
            if(info.ne.0) call stop_all(t_r,'Workspace queiry failed')
            lwork=int(work(1))+1
            deallocate(work)
            allocate(work(lwork))
            call dsyev('N','U',nSites,h0Eigenvecs,nSites,h0Eigenvals,Work,lWork,info)
            if(info.ne.0) call stop_all(t_r,'Diag failed')
            deallocate(work)
!            call writevector(h0Eigenvals,'Core eigenvals')
            !Now count the number of allowed occupations corresponding to degenerate CS mean-field solutions
            !This counts number of totally occupied sites rather than electrons
            N_occs = 0 
            i = 1
            do while(i.le.nSites)   !i wants to always be pointing at the first element in a degenerate set
                j = i + 1
    !            write(6,*) i
                do while((j.le.nSites).and.(abs(h0Eigenvals(min(nSites,j))-h0Eigenvals(i)).lt.1.0e-10_dp))
                    !This loop will count the degenerate space
                    j = j + 1   
                enddo
                !j should not point at the end of the degenerate set which starts at i and ends at j-1
                !Now we have a CS number of sites
                if((j-1.lt.((nSites/2)+1)).and.(j-1.gt.int(0.2_dp*(nSites/2)))) then
                    !Constraint so that we only look at occupations from 0.2 x half filling to half filling
                    N_occs = N_occs + 1
                endif
                i = j    !Skip to end of degenerate set
            enddo

            write(6,"(A,I6)") "Number of allowed unique closed shell orbital configurations: ",N_occs

            allocate(allowed_occs(N_occs))
            allowed_occs(:) = 0
            !Now actually fill the array with these allowed occupations
            k = N_occs   !The slot for the occupation number
            i = 1
            do while(i.le.nSites)
                j = i + 1
                do while((j.le.nSites).and.(abs(h0Eigenvals(min(nSites,j))-h0Eigenvals(i)).lt.1.0e-10_dp))
                    !This loop will count the degenerate space
                    j = j + 1   
                enddo
                !j should now point at the end of the degenerate set which starts at i and ends at j-1
                !Now we have a CS number of sites
                if((j-1.lt.((nSites/2)+1)).and.(j-1.gt.int(0.2_dp*(nSites/2)))) then
                    !Constraint so that we only look at occupations from 0.2 x half filling to half filling
                    !Fill things up from low occupation to large, so low occupations go into final slot
                    if(tRampDownOcc) then
                        allowed_occs(k) = j - 1
                    else
                        allowed_occs(N_occs-k+1) = j - 1
                    endif
                    k = k - 1
                endif
                i = j    !Skip to end of degenerate set
            enddo

            write(6,*) "AllowedOccs: ",allowed_occs(:)

            if(k.ne.0) call stop_all(t_r,"Error in calculating occupations")
            deallocate(h0Eigenvecs,h0Eigenvals)
        endif
        if(mod(nSites,nImp).ne.0) call stop_all(t_r,'Number of sites should be factor of impurity size')

        !Count the hopping parameters available
        nHoppingsImp = 0
        nHoppingsEnv = 0
        do iImp=1,nImp
            do i=1,nImp
                !This counts hopping elements within the impurity cluster, which is defined as the first nImp sites
                if(abs(h0(i,iImp)).gt.1.0e-8_dp) nHoppingsImp = nHoppingsImp + 1
            enddo
            do i=nImp+1,nSites
                !This counts hopping elements from impurity to the environment, which is defined as the indices after the nImp sites
                if(abs(h0(i,iImp)).gt.1.0e-8_dp) nHoppingsEnv = nHoppingsEnv + 1
            enddo
        enddo
        write(6,"(A,I6,A,I6,A)") "Connections of impurity sites via t: ",nHoppingsImp," within imp, ",nHoppingsEnv," to env"

    end subroutine find_occs

    !Add the local potential striped across the core hamiltonian 
    !(only possible with translational invariance)
    !If tAdd, then the correlation potential is added to the core potential, otherwise it is subtracted
    subroutine add_localpot(core,core_v,CorrPot,tAdd,core_b,core_v_b,CorrPot_b)
        implicit none
        real(dp) , intent(in) :: core(nSites,nSites)
        real(dp) , intent(out) :: core_v(nSites,nSites)
        real(dp) , intent(in) :: CorrPot(nImp,nImp)
        logical , intent(in), optional :: tAdd
        real(dp) , intent(in), optional :: core_b(nSites,nSites)
        real(dp) , intent(out), optional :: core_v_b(nSites,nSites)
        real(dp) , intent(in), optional :: CorrPot_b(nImp,nImp)
        real(dp), allocatable :: temp(:,:)
        integer :: i,j,k,a,b,indi,indj
        logical :: tAdd_
        character(len=*) , parameter :: t_r='add_localpot'

        if(tReadSystem) then
            return
        endif
            
        core_v(:,:) = 0.0_dp
        if(tUHF) then
            if(.not.(present(core_b).and.present(core_v_b).and.present(CorrPot_b))) then
                call stop_all(t_r,'If using UHF, should be passing these through')
            endif
        endif

        if(present(tAdd)) then
            tAdd_ = tAdd
        else
            tAdd_ = .true.
        endif

        if(LatticeDim.eq.1) then
            !Construct new hamiltonian which is block diagonal in the local potential
            do k=0,(nSites/nImp)-1
                do i=1,nImp
                    do j=1,nImp
                        if(tAdd_) then
                            core_v((k*nImp)+i,(k*nImp)+j)=CorrPot(i,j)
                        else
                            core_v((k*nImp)+i,(k*nImp)+j)=-CorrPot(i,j)
                        endif
                    enddo
                enddo
            enddo

            !Add this to the original mean-field hamiltonian
            do i=1,nSites
                do j=1,nSites
                    core_v(i,j) = core_v(i,j) + core(i,j)
                enddo
            enddo

            if(tUHF) then
                do k=0,(nSites/nImp)-1
                    do i=1,nImp
                        do j=1,nImp
                            if(tAdd_) then
                                core_v_b((k*nImp)+i,(k*nImp)+j)=CorrPot_b(i,j)
                            else
                                core_v_b((k*nImp)+i,(k*nImp)+j)=-CorrPot_b(i,j)
                            endif
                        enddo
                    enddo
                enddo

                !Add this to the original mean-field hamiltonian
                do i=1,nSites
                    do j=1,nSites
                        core_v_b(i,j) = core_v_b(i,j) + core_b(i,j)
                    enddo
                enddo
            endif

        elseif(LatticeDim.eq.2) then
            !2D lattices.

            if(tTiltedLattice) then
               
                allocate(temp(nSites,nSites))
                temp(:,:) = core(:,:)
                call Mat_to_lattice_order(temp)
                
                !Add correlation potential to Core_v
                Core_v(:,:) = temp(:,:)

                !Tile through space
                do i = 1,nSites
                    do j = 1,nSites
                        !TD_Imp_Lat gives the element of the v_loc which should be added here
                        !(Row major)
                        if(TD_Imp_Lat(j,i).ne.0) then

                            !Convert these into the actual values of each dimension
                            b = mod(TD_Imp_Lat(j,i)-1,nImp) + 1
                            a = ((TD_Imp_Lat(j,i)-1)/nImp) + 1
                            !write(6,*) TD_Imp_Lat(j,i),
                            !write(6,*) "a: ",a
                            !write(6,*) "b: ",b

                            if(tAdd_) then
                                Core_v(j,i) = Core_v(j,i) + CorrPot(a,b)*TD_Imp_Phase(j,i)
                            else
                                Core_v(j,i) = Core_v(j,i) - CorrPot(a,b)*TD_Imp_Phase(j,i)
                            endif
                        endif
                    enddo
                enddo

                !Transform both core_v and core back to the impurity ordering
                call Mat_to_imp_order(Core_v)

                if(tUHF) then
                    temp(:,:) = core_b(:,:)
                    call Mat_to_lattice_order(temp)
                    
                    !Add correlation potential to Core_v
                    Core_v_b(:,:) = temp(:,:)

                    !Tile through space
                    do i = 1,nSites
                        do j = 1,nSites
                            !TD_Imp_Lat gives the element of the v_loc which should be added here
                            !(Row major)
                            if(TD_Imp_Lat(j,i).ne.0) then
                                !Convert these into the actual values of each dimension
                                b = mod(TD_Imp_Lat(j,i)-1,nImp) + 1
                                a = ((TD_Imp_Lat(j,i)-1)/nImp) + 1
                                !write(6,*) TD_Imp_Lat(j,i),
                                !write(6,*) "a: ",a
                                !write(6,*) "b: ",b

                                if(tAdd_) then
                                    Core_v_b(j,i) = Core_v_b(j,i) + CorrPot_b(a,b)*TD_Imp_Phase(j,i)
                                else
                                    Core_v_b(j,i) = Core_v_b(j,i) - CorrPot_b(a,b)*TD_Imp_Phase(j,i)
                                endif
                            endif
                        enddo
                    enddo

                    !Transform both core_v and core back to the impurity ordering
                    call Mat_to_imp_order(Core_v_b)

                endif

                deallocate(temp)
            else
                !Square lattice

                core_v(:,:) = core(:,:)
                if(tUHF) core_v_b(:,:) = core_b(:,:)

                do k = 1,iImpRepeats
                    do i = 1,nImp
                        Indi = StripedImpIndices(i,k)
                        do j = 1,nImp
                            Indj = StripedImpIndices(j,k)

                            if(tAdd_) then
                                Core_v(indj,indi) = Core_v(indj,indi) + CorrPot(j,i)
                            else
                                Core_v(indj,indi) = Core_v(indj,indi) - CorrPot(j,i)
                            endif

                            if(tUHF) then
                                if(tAdd_) then
                                    core_v_b(indj,indi) = core_v_b(indj,indi) + CorrPot_b(j,i)
                                else
                                    core_v_b(indj,indi) = core_v_b(indj,indi) - CorrPot_b(j,i)
                                endif
                            endif
                        enddo
                    enddo
                enddo

            endif
        endif

    end subroutine add_localpot
    
    
    !Add a *complex* local potential striped across a *real* hamiltonian 
    !(only possible with translational invariance)
    !If tAdd, then the correlation potential is added to the core potential, otherwise it is subtracted
    subroutine add_localpot_comp_inplace(core_v,CorrPot,tAdd)
        implicit none
        complex(dp) , intent(inout) :: core_v(nSites,nSites)
        complex(dp) , intent(in) :: CorrPot(nImp,nImp)
        logical , intent(in), optional :: tAdd
        integer :: i,j,k,a,b,indi,indj
        logical :: tAdd_

        if(tReadSystem) then
            return
        endif

        if(present(tAdd)) then
            tAdd_ = tAdd
        else
            tAdd_ = .true.
        endif

        if(LatticeDim.eq.1) then
            !Construct new hamiltonian which is block diagonal in the local potential
            do k=0,(nSites/nImp)-1
                do i=1,nImp
                    do j=1,nImp
                        if(tAdd_) then
                            core_v((k*nImp)+i,(k*nImp)+j)=core_v((k*nImp)+i,(k*nImp)+j) + CorrPot(i,j)
                        else
                            core_v((k*nImp)+i,(k*nImp)+j)=core_v((k*nImp)+i,(k*nImp)+j) - CorrPot(i,j)
                        endif
                    enddo
                enddo
            enddo
        elseif(LatticeDim.eq.2) then

            if(tTiltedLattice) then

                call Mat_to_lattice_order_comp(core_v)

                do i = 1,nSites
                    do j = 1,nSites
                        if(TD_Imp_Lat(j,i).ne.0) then
                            !Convert these into the actual values of each dimension
                            b = mod(TD_Imp_Lat(j,i)-1,nImp) + 1
                            a = ((TD_Imp_Lat(j,i)-1)/nImp) + 1
                            !write(6,*) TD_Imp_Lat(j,i),
                            !write(6,*) "a: ",a
                            !write(6,*) "b: ",b

                            if(tAdd_) then
                                Core_v(j,i) = Core_v(j,i) + CorrPot(a,b)*TD_Imp_Phase(j,i)
                            else
                                Core_v(j,i) = Core_v(j,i) - CorrPot(a,b)*TD_Imp_Phase(j,i)
                            endif
                        endif
                    enddo
                enddo

                !Transform both core_v and core back to the impurity ordering
                call Mat_to_imp_order_comp(Core_v)
            else
                !Square lattice
                
                do k = 1,iImpRepeats
                    do i = 1,nImp
                        Indi = StripedImpIndices(i,k)
                        do j = 1,nImp

                            Indj = StripedImpIndices(j,k)

                            if(tAdd_) then
                                Core_v(indj,indi) = Core_v(indj,indi) + CorrPot(j,i)
                            else
                                Core_v(indj,indi) = Core_v(indj,indi) - CorrPot(j,i)
                            endif

                        enddo
                    enddo
                enddo


            endif

        endif

    end subroutine add_localpot_comp_inplace

    !Add a *complex* local potential striped across a *real* hamiltonian 
    !(only possible with translational invariance)
    !If tAdd, then the correlation potential is added to the core potential, otherwise it is subtracted
    subroutine add_localpot_comp(core,core_v,CorrPot,tAdd)
        implicit none
        complex(dp) , intent(in) :: core(nSites,nSites)
        complex(dp) , intent(out) :: core_v(nSites,nSites)
        complex(dp) , intent(in) :: CorrPot(nImp,nImp)
        logical , intent(in), optional :: tAdd
        complex(dp), allocatable :: temp(:,:)
        integer :: i,j,k,a,b,indi,indj
        logical :: tAdd_

        if(tReadSystem) then
!            core_v(:,:) = core_v(:,:)
            return
        endif

        if(present(tAdd)) then
            tAdd_ = tAdd
        else
            tAdd_ = .true.
        endif

        if(LatticeDim.eq.1) then
            !Construct new hamiltonian which is block diagonal in the local potential
            core_v = zzero
            do k=0,(nSites/nImp)-1
                do i=1,nImp
                    do j=1,nImp
                        if(tAdd_) then
                            core_v((k*nImp)+i,(k*nImp)+j)=CorrPot(i,j)
                        else
                            core_v((k*nImp)+i,(k*nImp)+j)=-CorrPot(i,j)
                        endif
                    enddo
                enddo
            enddo

            !Add this to the original mean-field hamiltonian
            do i=1,nSites
                do j=1,nSites
                    core_v(i,j) = core_v(i,j) + core(i,j)
                enddo
            enddo
        elseif(LatticeDim.eq.2) then

            if(tTiltedLattice) then

                allocate(temp(nSites,nSites))
                temp(:,:) = core(:,:)
                call Mat_to_lattice_order_comp(temp)

                Core_v(:,:) = temp(:,:)
                !Tile through space
                do i = 1,nSites
                    do j = 1,nSites
                        !TD_Imp_Lat gives the element of the v_loc which should be added here
                        !(Row major)
                        if(TD_Imp_Lat(j,i).ne.0) then

                            !Convert these into the actual values of each dimension
                            b = mod(TD_Imp_Lat(j,i)-1,nImp) + 1
                            a = ((TD_Imp_Lat(j,i)-1)/nImp) + 1
                            !write(6,*) TD_Imp_Lat(j,i),
                            !write(6,*) "a: ",a
                            !write(6,*) "b: ",b

                            if(tAdd_) then
                                Core_v(j,i) = Core_v(j,i) + CorrPot(a,b)*TD_Imp_Phase(j,i)
                            else
                                Core_v(j,i) = Core_v(j,i) - CorrPot(a,b)*TD_Imp_Phase(j,i)
                            endif
                        endif
                    enddo
                enddo

                !Transform both core_v and core back to the impurity ordering
                call Mat_to_imp_order_comp(Core_v)
                deallocate(temp)
            else
                !2D Square lattice
                core_v(:,:) = core(:,:)

                do k = 1,iImpRepeats
                    do i = 1,nImp
                        Indi = StripedImpIndices(i,k)
                        do j = 1,nImp
                            Indj = StripedImpIndices(j,k)

                            if(tAdd_) then
                                Core_v(indj,indi) = Core_v(indj,indi) + CorrPot(j,i)
                            else
                                Core_v(indj,indi) = Core_v(indj,indi) - CorrPot(j,i)
                            endif

                        enddo
                    enddo
                enddo


            endif
        endif

    end subroutine add_localpot_comp

    !Run a full HF, including mean-field on-site repulsion term in the fock matrix
    !These are stored in FullHFOrbs and FullHFEnergies
    subroutine run_true_hf()
        use DetTools, only: GetHFAntisymInt_spinorb
        implicit none
        real(dp) :: HEl,PDiff,fockel
        real(dp), allocatable :: Work(:),OccOrbs_HF(:,:),PMatrix_old(:,:),PMatrix(:,:)
        real(dp), allocatable :: fock(:,:),temp(:,:),h0HF(:,:)
        integer :: i,lWork,info,ex(2,2),j,nIter
        logical :: tFailedSCF
        character(len=*), parameter :: t_r='run_true_hf'

        write(6,"(A)") "Constructing full HF solution. DMET will start from core hamiltonian solution."

        !Construct fock matrix
        !The fock matrix is just the core hamiltonian (without the fitted potential) + diag(1/2 U * rdm(i,i)) on the diagonals
        allocate(fock(nSites,nSites))
        fock(:,:) = h0(:,:) !Core hamiltonian
        if(.not.tAnderson) then
            do i=1,nSites
                !Include the on-site repulsion
                fock(i,i) = fock(i,i) + U * 0.5_dp * (real(NEl,dp)/real(nSites,dp))
            enddo
!        elseif(tChemPot) then
!            fock(1,1) = fock(1,1) - (U/2.0_dp)
        endif
        
        if(allocated(FullHFOrbs)) then
            deallocate(FullHFOrbs,FullHFEnergies)
        endif
        allocate(FullHFOrbs(nSites,nSites)) !The orbitals from the diagonalization of the fock matrix
        allocate(FullHFEnergies(nSites))
        !call writematrix(fock,'fock',.true.)

        !Now just diagonalise this fock matrix, rather than use diis
        !First, diagonalize one-body hamiltonian
        FullHFOrbs(:,:) = Fock(:,:)
        FullHFEnergies(:) = 0.0_dp
        allocate(Work(1))
        lWork=-1
        info=0
        call dsyev('V','U',nSites,FullHFOrbs,nSites,FullHFEnergies,Work,lWork,info)
        if(info.ne.0) call stop_all(t_r,'Workspace queiry failed')
        lwork=int(work(1))+1
        deallocate(work)
        allocate(work(lwork))
        call dsyev('V','U',nSites,FullHFOrbs,nSites,FullHFEnergies,Work,lWork,info)
        if(info.ne.0) call stop_all(t_r,'Diag failed')
        deallocate(work)

        PDiff = 1.0_dp
        allocate(OccOrbs_HF(nSites,nOcc))
        OccOrbs_HF(:,:) = FullHFOrbs(:,1:nOcc)

        allocate(PMatrix_old(nSites,nSites))
        allocate(PMatrix(nSites,nSites))

        !Calculate initial trial P matrix:
        call dgemm('N','T',nSites,nSites,nOcc,1.0_dp,OccOrbs_HF,nSites,OccOrbs_HF,nSites,0.0_dp,PMatrix_old,nSites)
        !call writevector(FullHFEnergies(1:10),'Initial HF eigenvalues')

        tFailedSCF = .false.
        nIter = 0
        do while(PDiff.gt.1.0e-8_dp)
            nIter = nIter + 1
            if(nIter.gt.1000) then
                call warning(t_r,'Failed to converge SCF. Exiting calculation of true HF orbitals...')
                tFailedSCF = .true.
                exit
            endif
            FullHFOrbs(:,:) = h0(:,:)
            if(.not.tAnderson) then
                do i=1,nSites
                    FullHFOrbs(i,i) = FullHFOrbs(i,i) + PMatrix_old(i,i)*U
                enddo
            else
                !In the anderson model, we only have an on-site repulsion on the first site
                FullHFOrbs(1,1) = FullHFOrbs(1,1) + PMatrix_old(1,1)*U
            endif
            FullHFEnergies(:) = 0.0_dp
            allocate(Work(1))
            lWork=-1
            info=0
            call dsyev('V','U',nSites,FullHFOrbs,nSites,FullHFEnergies,Work,lWork,info)
            if(info.ne.0) call stop_all(t_r,'Workspace queiry failed')
            lwork=int(work(1))+1
            deallocate(work)
            allocate(work(lwork))
            call dsyev('V','U',nSites,FullHFOrbs,nSites,FullHFEnergies,Work,lWork,info)
            if(info.ne.0) call stop_all(t_r,'Diag failed')
            deallocate(work)
        
            OccOrbs_HF(:,:) = FullHFOrbs(:,1:nOcc)

            !Create new PMatrix
            call dgemm('N','T',nSites,nSites,nOcc,1.0_dp,OccOrbs_HF,nSites,OccOrbs_HF,nSites,0.0_dp,PMatrix,nSites)

            PDiff = 0.0_dp
            do i=1,nSites
                do j=1,nSites
                    PDiff = PDiff + abs(PMatrix(i,j)-PMatrix_old(i,j))
                enddo
            enddo
            PMatrix_old(:,:) = PMatrix(:,:)
            !write(6,*) "PDiff: ",PDiff
            !call writevector(FullHFEnergies(1:10),'HF eigenvalues')
        enddo
        deallocate(PMatrix,PMatrix_old,OccOrbs_HF)

        !write(6,*) "Full HF Orbs: "
        !call writematrix(FullHFOrbs,'FullHFOrbs',.true.)
            
        if(.not.tFailedSCF) then
            write(6,*) "nOCC", nOcc
            write(6,*) "*True* Fock eigenvalues around fermi level: "
            do i=max(1,nOcc-3),nOcc
                write(6,*) FullHFEnergies(i),"*"
            enddo
            do i=nOcc+1,min(nSites,nOcc+3)
                write(6,*) FullHFEnergies(i)
            enddo

            !Convert core hamiltonian into HF basis
            allocate(temp(nSites,nSites))
            allocate(h0HF(nSites,nSites))
!            if(tChemPot.and.tAnderson) h0(1,1) = h0(1,1) - U/2.0_dp
            call dgemm('t','n',nSites,nSites,nSites,1.0_dp,FullHFOrbs,nSites,h0,nSites,0.0_dp,temp,nSites)
            call dgemm('n','n',nSites,nSites,nSites,1.0_dp,temp,nSites,FullHFOrbs,nSites,0.0_dp,h0HF,nSites)
!            if(tChemPot.and.tAnderson) h0(1,1) = h0(1,1) + U/2.0_dp
            deallocate(temp)

            if(.true.) then
                !Generate fock eigenvalues and see if they are the same
                do i=1,nSites
                    fockel = h0HF(i,i)
                    do j=1,nel
                        ex(1,1) = j
                        ex(1,2) = i*2
                        ex(2,1) = j
                        ex(2,2) = i*2
                        HEl = GetHFAntisymInt_spinorb(ex,FullHFOrbs)
                        fockel = fockel + HEl
                    enddo
                    !write(6,*) "Fock eigenvalue calculated: ",i,fockel
                    if(abs(fockel-FullHFEnergies(i)).gt.1.0e-7_dp) then
                        call stop_all(t_r,'HF solution not correct - fock eigenvalues do not agree')
                    endif
                enddo
            endif

            !Now calculate HF energy:
            HFEnergy = 0.0_dp
            do i=1,nOcc
                HFEnergy = HFEnergy + h0HF(i,i)*2.0_dp
            enddo
            do i=1,nel
                do j=1,nel
                    ex(1,1) = i
                    ex(1,2) = j
                    ex(2,1) = i
                    ex(2,2) = j
                    HEl = GetHFAntisymInt_spinorb(ex,FullHFOrbs)

                    HFEnergy = HFEnergy + 0.5_dp*HEl 
                enddo
            enddo
            write(6,*) "HF energy from core hamiltonian: ",HFEnergy

            HFEnergy = 0.0_dp
            do i=1,nOcc
                HFEnergy = HFEnergy + 2.0_dp*FullHFEnergies(i)
            enddo
            do i=1,nel
                do j=1,nel
                    ex(1,1) = i
                    ex(1,2) = j
                    ex(2,1) = i
                    ex(2,2) = j
                    HEl = GetHFAntisymInt_spinorb(ex,FullHFOrbs)
                    HFEnergy = HFEnergy - HEl*0.5_dp
                enddo
            enddo
            write(6,*) "HF energy from fock eigenvalues: ",HFEnergy

            deallocate(h0HF)
        else
            deallocate(FullHFOrbs,FullHFEnergies)
            call run_hf(0)
        endif

        deallocate(fock)

    end subroutine run_true_hf

    !Read orbitals from file
    subroutine read_orbitals()
        use utils, only: get_free_unit
        implicit none
        integer :: i,j,iunit
        real(dp), allocatable :: OccOrbs(:,:),temp(:,:)
        real(dp) :: Occupancy
        logical :: exists
        character(len=*), parameter :: t_r="read_orbitals"

        write(6,*) "Reading HF orbitals from disk..."

        if(.not.allocated(HFOrbs)) allocate(HFOrbs(nSites,nSites)) !The orbitals
        HFOrbs(:,:) = zero  
        HFEnergies(:) = zero

        if(tUHF) then
            if(.not.allocated(HFOrbs_b)) allocate(HFOrbs_b(nSites,nSites))
            HFOrbs_b(:,:) = zero
            HFEnergies_b(:) = zero
            inquire(file='FinalOrbitals_a.dat',exist=exists)
            if(.not.exists) call stop_all(t_r,'Alpha orbital file cannot be found')
            iunit = get_free_unit()
            open(iunit,file='FinalOrbitals_a.dat',status='old',action='read')
            do i=1,nSites
                read(iunit,*) HFOrbs(i,:)
            enddo
            close(iunit)
            inquire(file='FinalOrbitals_b.dat',exist=exists)
            if(.not.exists) call stop_all(t_r,'Beta orbital file cannot be found')
            open(iunit,file='FinalOrbitals_b.dat',status='old',action='read')
            do i=1,nSites
                read(iunit,*) HFOrbs_b(i,:)
            enddo
            close(iunit)

            !Eigenvalues
            inquire(file='FinalEigenvalue_a.dat',exist=exists)
            if(.not.exists) call stop_all(t_r,'Alpha eigenvalue file cannot be found')
            open(iunit,file='FinalEigenvalue_a.dat',status='old',action='read')
            read(iunit,*) HFEnergies(:)
            close(iunit)
            inquire(file='FinalEigenvalue_b.dat',exist=exists)
            if(.not.exists) call stop_all(t_r,'Beta eigenvalue file cannot be found')
            open(iunit,file='FinalEigenvalue_b.dat',status='old',action='read')
            read(iunit,*) HFEnergies_b(:)
            close(iunit)
        else
            inquire(file='FinalOrbitals.dat',exist=exists)
            if(.not.exists) call stop_all(t_r,'Orbital file cannot be found')
            iunit = get_free_unit()
            open(iunit,file='FinalOrbitals.dat',status='old',action='read')
            do i=1,nSites
                read(iunit,*) HFOrbs(i,:)
            enddo
            close(iunit)
            
            inquire(file='FinalEigenvalue.dat',exist=exists)
            if(.not.exists) call stop_all(t_r,'Eigenvalue file cannot be found')
            open(iunit,file='FinalEigenvalue.dat',status='old',action='read')
            read(iunit,*) HFEnergies(:)
            close(iunit)
        endif

        write(6,*) "nOCC", nOcc
        write(6,*) "Fock eigenvalues around fermi level: "
        do i=max(1,nOcc-7),nOcc
            if(tUHF) then
                write(6,*) HFEnergies(i),HFEnergies_b(i),"*"
            else
                write(6,*) HFEnergies(i),"*"
            endif
        enddo
        do i=nOcc+1,min(nSites,nOcc+7)
            if(tUHF) then
                write(6,*) HFEnergies(i),HFEnergies_b(i)
            else
                write(6,*) HFEnergies(i)
            endif
        enddo

        if(tCheck) then
            !Check that these orbitals are correct
            allocate(temp(nSites,nSites))
            call DGEMM('n','n',nSites,nSites,nSites,one,h0v,nSites,HFOrbs,nSites,zero,temp,nSites)

            do i = 1,nSites
                !Loop over orbitals
                do j = 1,nSites
                    !Loop over components of orbitals
                    if(abs((temp(j,i)/HFEnergies(i)) - HFOrbs(j,i)).gt.1.0e-8_dp) then
                        write(6,*) "Orbital: ",i
                        write(6,*) "Component: ",j
                        write(6,*) temp(j,i)/HFEnergies(i),HFOrbs(j,i),abs((temp(j,i)/HFEnergies(i)) - HFOrbs(j,i))
                        call stop_all(t_r,'Error in reading orbitals')
                    endif
                enddo
            enddo
            if(tUHF) then
                call DGEMM('n','n',nSites,nSites,nSites,one,h0v_b,nSites,HFOrbs_b,nSites,zero,temp,nSites)

                do i = 1,nSites
                    !Loop over orbitals
                    do j = 1,nSites
                        !Loop over components of orbitals
                        if(abs((temp(j,i)/HFEnergies_b(i)) - HFOrbs_b(j,i)).gt.1.0e-8_dp) then
                            write(6,*) "Orbital: ",i
                            write(6,*) "Component: ",j
                            write(6,*) temp(j,i)/HFEnergies_b(i),HFOrbs_b(j,i),abs((temp(j,i)/HFEnergies_b(i)) - HFOrbs_b(j,i))
                            call stop_all(t_r,'Error in reading beta orbitals')
                        endif
                    enddo
                enddo
            endif

            deallocate(temp)
            write(6,"(A)") "Orbitals read in correctly"
        endif
        
        !Now calculate the density matrix from the calculation based on double occupancy of the lowest lying nOcc orbitals
        !First, extract the occupied orbitals. Since eigenvalues are ordered in increasing order, these will be the first nOcc
        if(tUHF) then
            Occupancy = one
        else
            Occupancy = 2.0_dp
        endif
        allocate(OccOrbs(nSites,nOcc))
        OccOrbs(:,:) = HFOrbs(:,1:nOcc)
        !Now construct the density matrix in the original AO basis. The eigenvectors are given as AO x MO, so we want to contract out the
        !MO contributions in order to get the 1DM in the AO basis.
        call DGEMM('N','T',nSites,nSites,nOcc,Occupancy,OccOrbs,nSites,OccOrbs,nSites,zero,MeanFieldDM,nSites)
        if(tUHF) then
            OccOrbs(:,:) = HFOrbs_b(:,1:nOcc)
            call DGEMM('N','T',nSites,nSites,nOcc,Occupancy,OccOrbs,nSites,OccOrbs,nSites,zero,MeanFieldDM_b,nSites)
        endif

        ChemPot = (HFEnergies(nOcc) + HFEnergies(nOcc+1))/2.0_dp  !Chemical potential is half way between HOMO and LUMO
        HLGap = HFEnergies(nOcc+1)-HFEnergies(nOcc)   !one body HOMO-LUMO Gap
        if(tUHF) then
            ChemPot_b = (HFEnergies_b(nOcc) + HFEnergies_b(nOcc+1))/2.0_dp  !Chemical potential is half way between HOMO and LUMO
            HLGap_b = HFEnergies_b(nOcc+1)-HFEnergies_b(nOcc)   !one body HOMO-LUMO Gap
        endif

        write(6,"(A,F20.10)") "One-electron eigenvalue band-gap: ",HLGap
        if(tUHF) write(6,"(A,F20.10)") "One-electron beta eigenvalue band-gap: ",HLGap_b
        if(HLGap.lt.1.0e-6_dp) then
            write(6,"(A,G15.5,A)") "Warning. HL gap is: ",HLGap," Possible failure in assigning orbitals in degenerate set."
        endif

        deallocate(OccOrbs)

    end subroutine read_orbitals

    !Run a HF calculation on the entire system. In this case, it just consists of just diagonalizing the system rather than iterative DIIS (add later)
    subroutine run_hf(it)
        implicit none
        integer, intent(in) :: it
        real(dp), allocatable :: fock(:,:),OccOrbs(:,:),fock_b(:,:)
        integer :: i
        real(dp) :: Occupancy
        !complex(dp), allocatable :: k_ham(:,:)
        !character(len=*), parameter :: t_r="run_hf"

        !Construct fock matrix
        !The fock matrix is just the core hamiltonian (with the fitted potential) + diag(1/2 U * rdm(i,i)) on the diagonals
        allocate(fock(nSites,nSites))
        if(tUHF) allocate(fock_b(nSites,nSites))
        if(it.eq.0) then
            fock(:,:) = h0(:,:)
            if(tUHF) fock_b(:,:) = h0_b(:,:)
        else
            fock(:,:) = h0v(:,:)
            if(tUHF) fock_b(:,:) = h0v_b(:,:)
        endif
                    
!        call writematrix(h0v,'h0v',.true.)

        if(allocated(HFOrbs)) then
            !If HFOrbs is already allocated, then we want to maximize the overlap between orbitals of different iterations
        else
            allocate(HFOrbs(nSites,nSites)) !The orbitals from the diagonalization of the fock matrix
            if(tUHF) allocate(HFOrbs_b(nSites,nSites))
        endif
        !Don't include the 2 electron terms in the mean field. Therefore fock matrix is core hamiltonian
        !The diagonal on-site repulsion will just be mopped up by the correlation potential
!        do i=1,nSites
!            !Include the on-site repulsion
!            fock(i,i) = fock(i,i) + U * 0.5_dp * (NEl/real(nSites))
!        enddo

        !Now just diagonalise this fock matrix, rather than use diis
        !First, diagonalize one-body hamiltonian
        HFOrbs(:,:) = Fock(:,:)
        if(tChemPot.and.tAnderson) HFOrbs(1,1) = HFOrbs(1,1) + U/2.0_dp
        HFEnergies(:) = zero  
        if(tDiag_kspace.and.(LatticeDim.eq.1)) then
            call DiagOneEOp(HFOrbs,HFEnergies,nImp,nSites,tDiag_kspace)
        else
            call DiagOneEOp(HFOrbs,HFEnergies,nImp,nSites,.false.)
        endif
        if(tUHF) then
            HFOrbs_b(:,:) = Fock_b(:,:)
            HFEnergies_b(:) = zero
            if(tDiag_kspace.and.(LatticeDim.eq.1)) then
                call DiagOneEOp(HFOrbs_b,HFEnergies_b,nImp,nSites,tDiag_kspace)
            else
                call DiagOneEOp(HFOrbs_b,HFEnergies_b,nImp,nSites,.false.)
            endif
        endif
            
        write(6,*) "nOCC", nOcc
        write(6,*) "Fock eigenvalues around fermi level: "
        do i=max(1,nOcc-7),nOcc
            if(tUHF) then
                write(6,*) HFEnergies(i),HFEnergies_b(i),"*"
            else
                write(6,*) HFEnergies(i),"*"
            endif
        enddo
        do i=nOcc+1,min(nSites,nOcc+7)
            if(tUHF) then
                write(6,*) HFEnergies(i),HFEnergies_b(i)
            else
                write(6,*) HFEnergies(i)
            endif
        enddo

        if(tUHF) then
            Occupancy = one
        else
            Occupancy = 2.0_dp
        endif

        !Now calculate the density matrix from the calculation based on double occupancy of the lowest lying nOcc orbitals
        !First, extract the occupied orbitals. Since eigenvalues are ordered in increasing order, these will be the first nOcc
        allocate(OccOrbs(nSites,nOcc))
        OccOrbs(:,:) = HFOrbs(:,1:nOcc)
        !Now construct the density matrix in the original AO basis. The eigenvectors are given as AO x MO, so we want to contract out the
        !MO contributions in order to get the 1DM in the AO basis.
        call DGEMM('N','T',nSites,nSites,nOcc,Occupancy,OccOrbs,nSites,OccOrbs,nSites,zero,MeanFieldDM,nSites)

        if(tUHF) then
            OccOrbs(:,:) = HFOrbs_b(:,1:nOcc)
            !Now construct the density matrix in the original AO basis. The eigenvectors are given as AO x MO, so we want to contract out the
            !MO contributions in order to get the 1DM in the AO basis.
            call DGEMM('N','T',nSites,nSites,nOcc,Occupancy,OccOrbs,nSites,OccOrbs,nSites,0.0_dp,MeanFieldDM_b,nSites)
        endif


        ChemPot = (HFEnergies(nOcc) + HFEnergies(nOcc+1))/2.0_dp  !Chemical potential is half way between HOMO and LUMO
        HLGap = HFEnergies(nOcc+1)-HFEnergies(nOcc)   !one body HOMO-LUMO Gap
        if(tUHF) then
            ChemPot_b = (HFEnergies_b(nOcc) + HFEnergies_b(nOcc+1))/2.0_dp  !Chemical potential is half way between HOMO and LUMO
            HLGap_b = HFEnergies_b(nOcc+1)-HFEnergies_b(nOcc)   !one body HOMO-LUMO Gap
        endif

        if(HLGap.lt.1.0e-6_dp) then
            write(6,"(A,G15.5,A)") "Warning. HL gap is: ",HLGap," Possible failure in assigning orbitals in degenerate set."
        endif

        deallocate(Fock,OccOrbs)

    end subroutine run_hf

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
            if(tTiltedLattice) then
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

            !Reciprocal lattice vector in the x direction, defining a square grid
            RecipLattVecs(1,1) = 2.0_dp*pi/real(nImp_x,dp)
            RecipLattVecs(2,1) = 0.0_dp
            !Reciprocal lattice vector in the y-direction
            RecipLattVecs(1,2) = 0.0_dp
            RecipLattVecs(2,2) = 2.0_dp*pi/real(nImp_y,dp)
            BZVol = (2.0_dp*pi/real(nImp_x,dp))*(2.0_dp*pi/real(nImp_y,dp))

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
                write(6,*) "KPnt ",k,KPnts(1,k),KPnts(1,k)*(nImp/pi)
            else
                write(6,"(A,I5,4F13.5)") "KPnt ",k,KPnts(:,k),KPnts(:,k)*(nImp_x/pi)
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

                do i = 1,nSites
                    if(LatticeDim.eq.1) then
                        if(mod(i-1,nImp).ne.n) cycle    !This is the unit rotation between n and m
                        PrimLattVec(1) = real(i-1-mod(i-1,nImp))    !The real-space translation to the cell
!                        PrimLattVec(1) = real(i-1,dp)               !The real-space translation to this site
                    else
                        call SiteIndexToLatCoord_2DSquare(i,indx,indy)
                        !Since these indices are 1 indexed, we need to make them 0 indexed
                        !The *site* displacement vectors are indx-1 and indy-1. 
                        !The *cell* displacement vectors are indx-1-mod(indx-1,nImp_x) and the y version.
                        if(((mod(indx-1,nImp_x)*nImp_x)+mod(indy-1,nImp_y)).ne.n) cycle
                        PrimLattVec(1) = real(indx - 1 - mod(indx-1,nImp_x),dp)
                        PrimLattVec(2) = real(indy - 1 - mod(indy-1,nImp_y),dp)
!                        PrimLattVec(1) = real(indx - 1,dp)
!                        PrimLattVec(2) = real(indy - 1,dp)
!                        write(6,"(A,I5,A,2F10.4)") "Lattice site: ",i, " has coordinate: ",PrimLattVec(:)
                    endif
                    phase = ddot(LatticeDim,KPnts(:,k),1,PrimLattVec,1)
!                    RtoK_Rot(i,ind_1+mod(i,nImp)) = exp(dcmplx(zero,phase))/sqrt(real(nKPnts,dp))
                    RtoK_Rot(i,ind_1+n) = exp(dcmplx(zero,phase))/sqrt(real(nKPnts,dp))
                enddo
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
        
    !Compute k-space hamiltonian blocks from the content of h0v
    subroutine CreateKHamBlocks(k_Hams)
        implicit none
        complex(dp), intent(out) :: k_Hams(nImp,nImp,nKPnts)

        complex(dp), allocatable :: ztemp(:,:),k_Ham_tmp(:,:),CompHam(:,:)
        integer :: i,j,kPnt,ind_1,ind_2

        k_Hams(:,:,:) = zzero
        allocate(ztemp(nSites,nImp))
        allocate(k_Ham_tmp(nImp,nImp))
        allocate(CompHam(nSites,nSites))
        do i = 1,nSites
            do j = 1,nSites
                CompHam(j,i) = dcmplx(h0v(j,i),zero)
            enddo
        enddo
        do kPnt = 1,nKPnts
            ind_1 = ((kPnt-1)*nImp) + 1
            ind_2 = nImp*kPnt
            call ZGEMM('N','N',nSites,nImp,nSites,zone,CompHam,nSites,RtoK_Rot(:,ind_1:ind_2),nSites,zzero,ztemp,nSites)
            call ZGEMM('C','N',nImp,nImp,nSites,zone,RtoK_Rot(:,ind_1:ind_2),nSites,ztemp,nSites,zzero,k_Ham_tmp,nImp)
            k_Hams(:,:,kPnt) = k_Ham_tmp(:,:)
        enddo
        deallocate(k_Ham_tmp,ztemp,CompHam)
    end subroutine CreateKHamBlocks

    !Find the projection of the final one-electron orbitals onto each kpoint
    subroutine ProjectHFontoK()
        implicit none
        complex(dp), allocatable :: TempHF_Comp(:,:)
        integer :: i,j

        if(.not.allocated(HFtoKOrbs)) then
            allocate(HFtoKOrbs(nSites,nSites))
            if(tUHF) allocate(HFtoKOrbs_b(nSites,nSites))
        endif

        allocate(TempHF_Comp(nSites,nSites))
        TempHF_Comp(:,:) = zzero
        do i = 1,nSites
            do j = 1,nSites
                TempHF_Comp(j,i) = dcmplx(HFOrbs(j,i),zero)
            enddo
        enddo

        call ZGEMM('C','N',nSites,nSites,nSites,zone,RtoK_Rot,nSites,TempHF_Comp,nSites,zzero,HFtoKOrbs,nSites)
        if(tUHF) then
            TempHF_Comp(:,:) = zzero
            do i = 1,nSites
                do j = 1,nSites
                    TempHF_Comp(j,i) = dcmplx(HFOrbs_b(j,i),zero)
                enddo
            enddo
            call ZGEMM('C','N',nSites,nSites,nSites,zone,RtoK_Rot,nSites,TempHF_Comp,nSites,zzero,HFtoKOrbs_b,nSites)
        endif

        deallocate(TempHF_Comp)

        if(tWriteOut) call writematrixcomp(HFtoKOrbs,'HFtoKOrbs(k,i)',.false.)

    end subroutine ProjectHFontoK

    !This is for the end of a calculation, to get the kspace 1e orbitals
    subroutine GetKSpaceOrbs()
        use sort_mod, only: sort_d_i
        implicit none
        complex(dp), allocatable :: CompHam(:,:),k_Ham(:,:),ztemp(:,:),cWork(:)
        complex(dp), allocatable :: TempSchmidtBasis(:,:),ztemp2(:,:)
        real(dp), allocatable :: Work(:)
        integer :: i,j,k,SS_Period,lWork,ind_1,ind_2,info
        character(len=*), parameter :: t_r='GetKSpaceOrbs'

        write(6,"(A)") "Final diagonalization of the hamiltonian in kspace to get the complex orbitals"

        allocate(CompHam(nSites,nSites))
        do i = 1,nSites
            do j = 1,nSites
                CompHam(j,i) = dcmplx(h0v(j,i),zero)
            enddo
        enddo

        SS_Period = nImp    !The length of the supercell
        
        !Project this hamiltonian into kspace and diagonalize each block at a time
        allocate(k_Ham(SS_Period,SS_Period))
        if(.not.allocated(k_vecs)) then
            allocate(k_vecs(SS_Period,nSites))
        endif
        if(.not.allocated(k_HFEnergies)) allocate(k_HFEnergies(nSites))
        k_vecs(:,:) = zzero
        
        !Space for diagonalization
        allocate(ztemp(nSites,SS_Period))
        lwork = max(1,2*SS_Period-1)
        allocate(cWork(lWork))
        allocate(Work(max(1,3*SS_Period-2)))

        !Run though all kpoints
        do k = 1,nKPnts
            ind_1 = ((k-1)*SS_Period) + 1
            ind_2 = SS_Period*k
            !We now have the rotation matrix for the bands on this kpoint.
            !Rotate the hamiltonian into this basis
            call ZGEMM('N','N',nSites,SS_Period,nSites,zone,CompHam,nSites,RtoK_Rot(:,ind_1:ind_2),nSites,zzero,ztemp,nSites)
            call ZGEMM('C','N',SS_Period,SS_Period,nSites,zone,RtoK_Rot(:,ind_1:ind_2),nSites,ztemp,nSites,zzero,k_Ham,SS_Period)

            !Diagonalize this k-pure hamiltonian
            info = 0
            call ZHEEV('V','U',SS_Period,k_Ham,SS_Period,k_HFEnergies(ind_1:ind_2),cWork,lWork,Work,info)
            if(info.ne.0) call stop_all(t_r,'Diag failed')

            if(tWriteOut) then
                write(6,*) "For kpoint: ",k," Eigenvalues are:"
                do i = 0,SS_Period-1
                    write(6,*) k_HFEnergies(ind_1+i)
                enddo
                write(6,*) "Eigenvectors: "
                call writematrixcomp(k_Ham,'Eigenvec',.true.)
            endif

            k_vecs(:,ind_1:ind_2) = k_Ham(:,:)
        enddo
        deallocate(work,cWork,ztemp,k_Ham,CompHam)

        !We now have pure k-eigenvectors, as ( pure k-component : orbital number )

        !This map gives the index of the orbitals in terms of energy
        !KVec_EMapping(i) = index of ith lowest energy orbital
        if(.not.allocated(KVec_EMapping)) allocate(KVec_EMapping(nSites))
        do i = 1,nSites
            KVec_EMapping(i) = i
        enddo
        allocate(Work(nSites))
        Work(:) = k_HFEnergies(:)
        call sort_d_i(Work,KVec_EMapping,nSites)
        deallocate(Work)

        if(.not.allocated(KVec_InvEMap)) allocate(KVec_InvEMap(nSites))
        !KVec_InvEMapping takes as input the k orbital number, and returns its energetic order 
        do i = 1,nSites
            KVec_InvEMap(KVec_EMapping(i)) = i
        enddo

        if(tWriteOut) then
            call writevector(HFEnergies,'HFEnergies')
            call writevector(k_HFEnergies,'k_HFEnergies')
            call writevectorint(KVec_EMapping,'kVec_EMapping')
            call writevectorint(KVec_InvEMap,'kVec_InvEMap')
        endif

        do i = 1,nSites
            if(abs(HFEnergies(i)-k_HFEnergies(KVec_EMapping(i))).gt.1.0e-7) then
                call stop_all(t_r,'k-space HF energies do not match up with real space ones')
            endif
        enddo
        
        if(allocated(k_HFtoSchmidtTransform)) deallocate(k_HFtoSchmidtTransform)
        allocate(k_HFtoSchmidtTransform(nSites,nSites))

        allocate(ztemp(nSites,nSites))
        allocate(ztemp2(nSites,nSites))
        ztemp(:,:) = zzero
        !Construct full block diagonal representation of k-basis MOs
        do k = 1,nKPnts
            ind_1 = ((k-1)*SS_Period) + 1
            ind_2 = SS_Period*k
            ztemp(ind_1:ind_2,ind_1:ind_2) = k_vecs(:,ind_1:ind_2)
        enddo

        !First rotate k-space eigenvectors into AO basis
        call ZGEMM('N','N',nSites,nSites,nSites,zone,RtoK_Rot,nSites,ztemp,nSites,zzero,ztemp2,nSites)
        !Now rotate eigenvectors in AO basis into schmidt basis
        allocate(TempSchmidtBasis(nSites,nSites))
        TempSchmidtBasis(:,:) = zzero
        do j = 1,nSites
            do i = 1,nSites
                TempSchmidtBasis(i,j) = dcmplx(FullSchmidtBasis(i,j),zero)
            enddo
        enddo
        call ZGEMM('C','N',nSites,nSites,nSites,zone,TempSchmidtBasis,nSites,ztemp2,nSites,zzero,k_HFtoSchmidtTransform,nSites)
        deallocate(ztemp,ztemp2,TempSchmidtBasis)

        if(tWriteOut) call writematrixcomp(k_HFtoSchmidtTransform,'k_HFtoSchmidtTransform',.false.)

    end subroutine GetKSpaceOrbs
                
    !Given the lower triangle, make this matrix hermitian
    subroutine MakeBlockHermitian(Block,length)
        implicit none
        integer, intent(in) :: length
        complex(dp), intent(inout) :: Block(length,length)
        integer :: i,j

        do i = 1,nImp
            do j = i+1,nImp
                Block(i,j) = dconjg(Block(j,i))
            enddo
        enddo

    end subroutine MakeBlockHermitian
    
    !Add a non-local periodic set of real-space lattice parameters to a complex inplace hamiltonian
    subroutine add_Nonlocal_comp_inplace(ham,NonLocCoup,NonLocCoupLength,tAdd)
        implicit none
        integer, intent(in) :: NonLocCoupLength
        complex(dp), intent(inout) :: ham(nSites,nSites)
        complex(dp), intent(in) :: NonLocCoup(NonLocCoupLength)
        logical , intent(in), optional :: tAdd
        complex(dp), allocatable :: PartialBlock(:,:),FullBlocks(:,:,:),Block_T(:,:)
        real(dp) :: phase
        integer :: k,i,j,ind_1,ind_2,ind,ii,ind_v,iFullBlocks,ind_r,ind_c,ind_block
        character(len=*), parameter :: t_r='add_Nonlocal_comp_inplace'
            
        if(tPeriodic) then
            phase = 1.0_dp
        elseif(tAntiPeriodic) then
            phase = -1.0_dp
        endif

        if(iNonLocBlocks.eq.0) then
            !Use all blocks
            if(tOddFullNonlocCoups) then
                !iFullBlocks gives the number of full coupling unit cells, without the potentially special one which isn't repeated
                iFullBlocks = ((nSites/nImp)-2)/2
            else
                iFullBlocks = ((nSites/nImp)-1)/2
            endif
        else
            iFullBlocks = iNonLocBlocks
        endif
        allocate(FullBlocks(nImp,nImp,iFullBlocks))
        FullBlocks(:,:,:) = zzero
        if(tOddFullNonlocCoups) then
            !PartialBlock is the skew-symmetric block which is the furthest coupling if optimizing all, and there are an even number of kpoints
            allocate(PartialBlock(nImp,nImp))
            PartialBlock(:,:) = zzero
        endif

        !TODO: In the future, this can just use the var_to_couplingind routine
        ind = 0
        do k = 1,iFullBlocks
            do i = 1,nImp
                do j = 1,nImp
                    ind = ind + 1
                    if(ind.gt.NonLocCoupLength) then
                        write(6,*) "i,j,k: ",i,j,k
                        write(6,*) "ind: ",ind
                        stop 'Indexing error'
                    endif
                    call var_to_couplingind(ind,nImp,ind_r,ind_c,ind_block)
                    if(ind_block.ne.k) then
                        write(6,*) "k,i,j: ",k,i,j
                        write(6,*) "ind: ",ind
                        write(6,*) "ind_block, ind_r,ind_c: ",ind_block,ind_r,ind_c
                        stop 'error here: ind_block wrong'
                    endif
                    if(ind_r.ne.j) then
                        write(6,*) "k,i,j: ",k,i,j
                        write(6,*) "ind: ",ind
                        write(6,*) "ind_block, ind_r,ind_c: ",ind_block,ind_r,ind_c
                        stop 'error here 2: ind_r wrong'
                    endif
                    if(ind_c.ne.i) then
                        write(6,*) "k,i,j: ",k,i,j
                        write(6,*) "ind: ",ind
                        write(6,*) "ind_block, ind_r,ind_c: ",ind_block,ind_r,ind_c
                        stop 'error here 3: ind_c wrong'
                    endif
                    FullBlocks(j,i,k) = NonLocCoup(ind)
                enddo
            enddo
        enddo

        if(tOddFullNonlocCoups) then
            do i = 1,nImp
                do j = i+1,nImp
                    ind = ind + 1
                    if(ind.gt.NonLocCoupLength) stop 'Indexing error'
                    call var_to_couplingind(ind,nImp,ind_r,ind_c,ind_block)
                    if(ind_block.ne.iFullBlocks+1) then
                        write(6,*) "k,i,j: ",k,i,j
                        write(6,*) "ind: ",ind
                        write(6,*) "ind_block, ind_r,ind_c: ",ind_block,ind_r,ind_c
                        stop 'error here 4: ind_block special wrong'
                    endif
                    if(ind_r.ne.j) then
                        write(6,*) "k,i,j: ",k,i,j
                        write(6,*) "ind: ",ind
                        write(6,*) "ind_block, ind_r,ind_c: ",ind_block,ind_r,ind_c
                        stop 'error here 5: ind_r special wrong'
                    endif
                    if(ind_c.ne.i) then
                        write(6,*) "k,i,j: ",k,i,j
                        write(6,*) "ind: ",ind
                        write(6,*) "ind_block, ind_r,ind_c: ",ind_block,ind_r,ind_c
                        stop 'error here 6: ind_c special wrong'
                    endif
                    PartialBlock(j,i) = NonLocCoup(ind)
                    PartialBlock(i,j) = -NonLocCoup(ind)
                enddo
            enddo
        endif

        allocate(Block_T(nImp,nImp))
        do k = 1,nKPnts
            ind_1 = ((k-1)*nImp) + 1
            ind_2 = nImp*k

            do i = 1,iFullBlocks
                !Blocks to the right
                if(ind_1+(i*nImp).le.nSites) then
                    ham(ind_1:ind_2,ind_1+(i*nImp):ind_2+(i*nImp)) = FullBlocks(:,:,i)
                else
                    ham(ind_1:ind_2,ind_1+(i*nImp)-nSites:ind_2+(i*nImp)-nSites) = phase*FullBlocks(:,:,i)
                endif

                !Blocks to the left are flipped
                call TransposeBlock_z(FullBlocks(:,:,i),Block_T(:,:),nImp)
                if(ind_1-(i*nImp).ge.1) then
                    ham(ind_1:ind_2,ind_1-(i*nImp):ind_2-(i*nImp)) = Block_T(:,:)
                else
                    ham(ind_1:ind_2,ind_1-(i*nImp)+nSites:ind_2-(i*nImp)+nSites) = phase*Block_T(:,:)
                endif
            enddo

            if(tOddFullNonlocCoups) then
                !Add in the partial block. This is only used once. This should only really occur with APBCs.
                !...otherwise block would be symmetric?
                if(ind_1+((iFullBlocks+1)*nImp).le.nSites) then
                    ham(ind_1:ind_2,ind_1+((iFullBlocks+1)*nImp):ind_2+((iFullBlocks+1)*nImp)) = PartialBlock(:,:)
                else
                    call TransposeBlock_z(PartialBlock(:,:),Block_T(:,:),nImp)
                    ham(ind_1:ind_2,ind_1+((iFullBlocks+1)*nImp)-nSites:ind_2+((iFullBlocks+1)*nImp)-nSites) = Block_T(:,:)
                endif
            endif
        enddo
        deallocate(Block_T,FullBlocks)
        if(tOddFullNonlocCoups) deallocate(PartialBlock)

        if(tCheck) then
            !Check for hermiticity
            do i = 1,nSites
                do j = i+1,nSites
                    if(abs(ham(i,j)-dconjg(ham(j,i))).gt.1.0e-8_dp) then
                        write(6,*) "i,j: ",i,j
                        write(6,*) "ham(i,j): ",ham(i,j)
                        write(6,*) "ham(j,i): ",ham(j,i)
                        call stop_all(t_r,'hamiltonian not hermitian')
                    endif
                enddo
            enddo
        endif

    end subroutine add_Nonlocal_comp_inplace

    !swap the off-diagonal matrix elements (and cc if necessary)
    subroutine TransposeBlock_z(Block,BlockT,nSize)
        implicit none
        integer, intent(in) :: nSize
        complex(dp), intent(in) :: Block(nSize,nSize)
        complex(dp), intent(out) :: BlockT(nSize,nSize)
        integer :: i,j

        BlockT(:,:) = zzero
        do i = 1,nSize
            do j = 1,nSize
                if(i.eq.j) then
                    !Leave diagonals alone
                    BlockT(i,j) = Block(i,j)
                else
                    BlockT(i,j) = dconjg(Block(j,i))
                endif
            enddo
        enddo

    end subroutine TransposeBlock_z

    !Given a variable i, find the block, row and column within the block it corresponds to, in the first blocked row of the matrix
    subroutine var_to_couplingind(i,blocksize,ind_r,ind_c,ind_block)
        implicit none
        integer, intent(in) :: i,blocksize
        integer, intent(out) :: ind_r,ind_c,ind_block
        integer :: j,k,l,ind
        logical :: tFound

        if(tOddFullNonlocCoups.and.(i.gt.(((nSites/blocksize)-2)/2)*(blocksize**2))) then
            ind_block = ((nSites/blocksize)-2)/2 + 1
            j = i - (((nSites/blocksize)-2)/2)*(blocksize**2)
            ind = 0
            tFound = .false.
            loop: do k = 1,blocksize
                do l = k+1,nImp
                    ind = ind + 1
                    if(ind.eq.j) then
                        ind_r = l
                        ind_c = k
                        tFound = .true.
                        exit loop
                    endif
                enddo
            enddo loop
            if(.not.tFound) then
                write(6,*) "i: ",i
                write(6,*) "j: ",j
                write(6,*) "previous elements: ",(((nSites/blocksize)-2)/2)*(blocksize**2)
                stop 'error in indexing'
            endif
        else
            !This should be easier.
            ind_block = (i-1)/(blocksize**2) + 1
            j = mod(i,blocksize**2)
            k = mod(j,blocksize)
            if(k.eq.0) then
                ind_r = blocksize
            else
                ind_r = k
            endif
            if(j.eq.0) then
                !end element of block
                ind_c = blocksize
            else
!                write(6,*) "j: ",j
!                write(6,*) "mod(blocksize-j,blocksize): ",mod(blocksize-j,blocksize)
                ind_c = (j-1)/blocksize + 1
!                ind_c = (j+mod(blocksize-j,blocksize))/blocksize
            endif
        endif

    end subroutine var_to_couplingind
    
  !Find the index of the coupling vector from the block, column and row in the block
    subroutine couplingind_to_var(ind_r,ind_c,ind_block,blocksize,ind)
        implicit none
        integer, intent(in) :: ind_r,ind_c,ind_block,blocksize
        integer, intent(out) :: ind

        if(tOddFullNonlocCoups) then
            if(ind_block.gt.(((nSites/nImp)-2)/2)) then
                !ind_block > iFullBlocks
                !We now need to consider the special block.
                if(ind_c.ge.ind_r) stop 'Sending in wrong half indices for half-block'
                ind = (ind_block-1)*(nImp**2) + (nImp*(nImp-1))/2 -   &
                    ((nImp-ind_c)*(nImp-ind_c+1))/2 + ind_r - ind_c
            else
                ind = (ind_block-1)*(nImp**2) + (nImp*(ind_c-1)) + ind_r
            endif
        else
            !simples
            ind = (ind_block-1)*(nImp**2) + (nImp*(ind_c-1)) + ind_r
        endif

    end subroutine couplingind_to_var


    !Add coupling from the impurity sites to the other sites in the lattice
    !This is done to maintain the translational symmetry of the original impurity unit cell
    !If a cluster calculation, then the lattice coupling is seperate for each impurity of the unit cell
    !If optimizing the eigenvalues, then this just returns the hamiltonian from the eigenvalues, while
    !maintaining the same eigenvectors of the original lattice hamiltonian
    subroutine AddPeriodicImpCoupling_RealSpace(Ham,nLat,iCoupLength,SS_Period,Couplings)
        implicit none
        integer, intent(in) :: nLat,iCoupLength,SS_Period
        real(dp), intent(inout) :: Ham(nLat,nLat)
        real(dp), intent(in) :: Couplings(iCoupLength,SS_Period)
        real(dp), allocatable :: temp(:),temp2(:,:)
        complex(dp), allocatable :: ctemp(:,:),cham(:,:)
        integer :: CoupInd,i,j,k,offset,ImpCoup
        real(dp) :: mu
        logical, parameter :: tOneCouplingSet=.true.
        character(len=*), parameter :: t_r='AddPeriodicImpCoupling_RealSpace'

        if(tOptGF_EVals) then
            if(nImp.gt.1) call stop_all(t_r,'This is not set up for > 1 imp - bug GHB')
            if(tDiag_KSpace) then
                !Here, the eigenvalues are the positive/negative eigenvalues.
                !Use these to construct the lattice hamtiltonian
                if(tConstrainKSym.and.(.not.tConstrainphsym)) then
                    if(tShift_Mesh.and.(iCoupLength.ne.(nSites/2))) call stop_all(t_r,'Wrong # of params')
                    if((.not.tShift_Mesh).and.(iCoupLength.ne.((nSites/2)+1))) call stop_all(t_r,'Wrong # of params')
                elseif(tConstrainKSym.and.tConstrainphsym) then
                    if(tShift_Mesh.and.(iCoupLength.ne.(nSites/4))) call stop_all(t_r,'Wrong # of params')
                    if((.not.tShift_Mesh).and.(iCoupLength.ne.(((nSites/2)+1)/2))) call stop_all(t_r,'Wrong # of params')
                elseif(tConstrainphsym.and.(.not.tConstrainksym)) then
                    if(iCoupLength.ne.(nSites/2)) call stop_all(t_r,'Wrong # of params')
                else
                    if(iCoupLength.ne.nSites) call stop_all(t_r,'Wrong # of params')
                endif

                allocate(ctemp(nSites,nSites))
                allocate(cham(nSites,nSites))
                cham(:,:) = zzero
                if(tConstrainKSym.and.(.not.tConstrainphsym)) then
                    do i = 1,nSites/2
                        cham(i,i) = dcmplx(Couplings(i,1),zero)
                    enddo
                    if(tShift_Mesh) then
                        !No gamma point. All k-points symmetric
                        do i = 1,nSites/2
                            cham(i+(nSites/2),i+(nSites/2)) = dcmplx(Couplings((nSites/2)-i+1,1),zero)
                        enddo
                    else
                        !Put in the gamma point
                        cham((nSites/2)+1,(nSites/2)+1) = dcmplx(Couplings((nSites/2)+1,1),zero)
                        !Now, add in the other kpoints (not point on BZ boundary at start)
                        do i = 2,nSites/2
                            cham(i+(nSites/2),i+(nSites/2)) = dcmplx(Couplings((nSites/2)-i+2,1),zero)
                        enddo
                    endif
                elseif(tConstrainKSym.and.tConstrainphsym) then
                    mu = U/2.0_dp
                    !First, copy the eigenvalues to the particle states
                    do i = 1,iCoupLength
                        cham(i,i) = dcmplx(Couplings(i,1),zero)
                    enddo
                    !Now for the hole states
                    do i = 1,iCoupLength
                        cham(iCoupLength+i,iCoupLength+i) = dcmplx(2.0_dp*mu-Couplings(iCoupLength-i+1,1),zero)
                    enddo
                    !Now for ksym
                    if(tShift_Mesh) then
                        !No gamma point. All k-points symmetric
                        do i = 1,nSites/2
                            cham(i+(nSites/2),i+(nSites/2)) = cham((nSites/2)-i+1,(nSites/2)-i+1)
                        enddo
                    else
                        !Put in the gamma point (actually, this has already been put in)
                        cham((nSites/2)+1,(nSites/2)+1) = cham((nSites/2)+1,(nSites/2)+1)
                        !Now, add in the other kpoints (not point on BZ boundary at start)
                        do i = 2,nSites/2
                            cham(i+(nSites/2),i+(nSites/2)) = cham((nSites/2)-i+2,(nSites/2)-i+2)
                        enddo
                    endif
                else
                    do i = 1,nSites
                        cham(i,i) = dcmplx(Couplings(i,1),zero)
                    enddo
                endif
                call ZGEMM('N','N',nSites,nSites,nSites,zone,RtoK_Rot,nSites,cham,nSites,zzero,ctemp,nSites)
                call ZGEMM('N','C',nSites,nSites,nSites,zone,ctemp,nSites,RtoK_Rot,nSites,zzero,cham,nSites)
                deallocate(ctemp)

                !What about the phase factor? This will be needed when we go to > 1 impurity
                if(nImp.gt.1) call stop_all(t_r,'Cannot deal with > 1 impurty, as we will have to deal with the phase factor')
                !write(6,*) cham(1,:)

                do i = 1,nSites
                    do j = 1,nSites
                        if(abs(aimag(cham(j,i))).gt.1.0e-7_dp) then
                            write(6,*) cham(j,i)
                            call stop_all(t_r,'why is this not real? - Losing k-inversion symmetry')
                        endif
                        ham(j,i) = real(cham(j,i),dp)
                    enddo
                enddo
                deallocate(cham)
            else
                !Keep things real
                ham(:,:) = zero
                do i = 1,nSites
                    ham(i,i) = Couplings(i,1)
                enddo
                allocate(temp2(nSites,nSites))
                call DGEMM('N','N',nSites,nSites,nSites,one,HFOrbs,nSites,ham,nSites,zero,temp2,nSites)
                call DGEMM('N','T',nSites,nSites,nSites,one,temp2,nSites,HFOrbs,nSites,zero,ham,nSites)
                deallocate(temp2)
            endif
        else

            if(LatticeDim.eq.2) call stop_all(t_r,'Lattice coupling not set up for 2D lattices yet')

            if(SS_Period.gt.nImp) then
                call stop_all(t_r,'Coupling cannot connect to more sites than in the impurity unit cell')
            endif

    !        !Since the couplings need to maintain periodicity, we need them to apply to sites forwards and backwards.
    !        !Antiperiodic BCs should mean that we flip sign if we go around the lattice.

            if(tOneCouplingSet) then
                !Here, we only have one coupling (all columns of Couplings array should be the same)
                !This is inspired by assuming the supercell changes in a hubbard model with long-range couplings
                !This will actually be identical to below in the case of a single impurity
                ! *-*-B-C-0-C-B-A
                ! *-*-A-B-C-0-C-B
                ! B-A-*-*-B-C-0-C
                ! C-B-*-*-A-B-C-0
                ! 0-C-B-A-*-*-B-C
                ! C-0-C-B-*-*-A-B
                ! B-C-0-C-B-A-*-*
                ! A-B-C-0-C-B-*-*
                !As below, the top right diagonal will be negative with APBCs

                if(iCoupLength.gt.((nLat-SS_Period)/2)) then
                    call stop_all(t_r,'Impurity coupling should not extend beyond half the lattice in each dimension')
                endif

                do i = 1,nLat
                    offset = mod(i,SS_Period)
                    do j = i+offset+1,i+iCoupLength
                        if(j.gt.nLat) exit
                        Ham(i,j) = Couplings(j-i,1)
                        Ham(j,i) = Couplings(j-i,1)
                    enddo
                enddo

                !Now for the terms which wrap around the chain
                do i = 1,iCoupLength
                    do j = iCoupLength-i+1,1,-1
                        if(tPeriodic) then
                            Ham(i,nLat-j+1) = Couplings(j+i-1,1)
                            Ham(nLat-j+1,i) = Couplings(j+i-1,1)
                        elseif(tAntiPeriodic) then
                            Ham(i,nLat-j+1) = - Couplings(j+i-1,1)
                            Ham(nLat-j+1,i) = - Couplings(j+i-1,1)
                        else
                            !OBCs
                            continue
                        endif
                    enddo
                enddo
            else

                !Max number of couplings given by (nSites-nImp)/(2 x nImp)
                if(iCoupLength.gt.((nLat-SS_Period)/(2*SS_Period))) then
                    call stop_all(t_r,'Impurity coupling should not extend beyond' &
                        //' half the number of cell repeats in each dimension')
                endif

                !Example below for nLat=18, iCoupLength=3, SS_Period=2, Couplings= (/A,B,C  a,b,c/)
                ! *-*-A-0-B-0-C-0-0-0-0-0-C-0-B-0-A-0
                ! *-*-0-a-0-b-0-c-0-0-0-0-0-c-0-b-0-a
                ! A-0-*-*-A-0-B-0-C-0-0-0-0-0-C-0-B-0
                ! 0-a-*-*-0-a-0-b-0-c-0-0-0-0-0-c-0-b
                ! B-0-C-0-*-*-A-0-B-0-C-0-0-0-0-0-C-0
                ! 0-b-0-c-*-*-0-a-0-b-0-c-0-0-0-0-0-c
                ! C-0-B-0-A-0-*-*-A-0-B-0-C-0-0-0-0-0
                ! 0-c-0-b-0-a-*-*-0-a-0-b-0-c-0-0-0-0
                ! 0-0-C-0-B-0-A-0-*-*-A-0-B-0-C-0-0-0
                ! 0-0-0-c-0-b-0-a-*-*-0-a-0-b-0-c-0-0
                !                     *-*-A-0-B-0-C-0
                !                     *-*-0-a-0-b-0-c
                !                         *-*-A-0-B-0
                !                         *-*-0-a-0-b
                !                             *-*-A-0
                !                             *-*-0-a
                !                                 *-*
                !                                 *-*

                !A coupling of iCoupLength=1, SS_Period=1, Couplings=(/1/) will give normal lattice hamiltonian
                !If APBCs, then the corner couplings will be negative

                !First, simply put in the direct couplings to neighbours in one direction
                do i = 1,nLat
                    ImpCoup = mod(i,SS_Period)    !Which impurity coupling are we considering
                    if(ImpCoup.eq.0) ImpCoup = SS_Period
                    do j = 1,iCoupLength
                        if((i+(j*SS_Period)).gt.nLat) exit    !Do not go "around" the ring. We will deal with this afterwards, since this
                                                                        !may require different boundary conditions
                        Ham(i,i+(j*SS_Period)) = Couplings(j,ImpCoup)
                        Ham(i+(j*SS_Period),i) = Couplings(j,ImpCoup)
                    enddo
                enddo

                !Now, also put in the coupling in the other direction
                do i = 1,iCoupLength*SS_Period
                    ImpCoup = mod(i,SS_Period)
                    if(ImpCoup.eq.0) ImpCoup = SS_Period

                    CoupInd = iCoupLength + 1
                    do j = nLat - (iCoupLength*SS_Period) + i, nLat, SS_Period
                        CoupInd = CoupInd - 1   !Start with last coupling value
                        if(tAntiPeriodic) then
                            Ham(i,j) = - Couplings(CoupInd,ImpCoup)
                            Ham(j,i) = - Couplings(CoupInd,ImpCoup)
                        elseif(tPeriodic) then
                            Ham(i,j) = Couplings(CoupInd,ImpCoup)
                            Ham(j,i) = Couplings(CoupInd,ImpCoup)
                        else
                            !OBCs
                            continue
                        endif
                    enddo
                enddo
                    
    !            call writematrix(Ham,'Ham_realspace_BeforeInterchange',.true.)
                
    !            !OPTIONAL!
    !            !Attempt an alternative tiling, essentially coupling the 'other' bath sites
    !            !This means that we should change the order of the couplings in each supercell block of the lattice
    !            allocate(temp(SS_Period))
    !            do i = 1,nLat   !Run through every lattice site
    !                offset=mod(i,SS_Period) + 1
    !                do j = i+offset,nLat,SS_Period 
    !                    temp(:) = Ham(i,j:j+SS_Period-1)
    !                    do k = 0,SS_Period-1
    !                        !Now, reverse the order
    !                        Ham(i,j+SS_Period-1-k) = temp(k+1)
    !                        Ham(j+SS_Period-1-k,i) = temp(k+1)
    !                    enddo
    !                enddo
    !            enddo
    !            deallocate(temp)

            endif
        endif

    end subroutine AddPeriodicImpCoupling_RealSpace

!Take a real-space hamiltonian on the lattice, and check that it obeys translational symmetry of the unit cells,
!whose size is given by the impurity size. Use through interface for real and imaginary hamiltonians.
    subroutine CheckRealSpaceTransInvar_z(ham)
        implicit none
        complex(dp), intent(in) :: ham(nSites,nSites)
        integer :: i,j,ind_1,ind_2,phase,LatInd_UC,k,ind_kx,ind_ky,transvec_x,transvec_y
        integer :: Curr_UCx,Curr_UCy,Curr_UC,p,Index_2,k_trans,PhaseChange
        character(len=*), parameter :: t_r='CheckRealSpaceTransInvar_z'

        !Check for translational invariance of the system
        !That is: is each supercell row equivalent when translated to the supercell ajacent to it?
        if(LatticeDim.eq.1) then
            do j = 1,nSites
                do i = 1,nSites
                    ind_1 = i
                    ind_2 = j
                    !Translate a supercell down
                    ind_1 = i + nImp
                    ind_2 = j + nImp
                    phase = 1
                    if(ind_1.gt.nSites) then
                        ind_1 = ind_1 - nSites
                        if(tAntiPeriodic) phase = -1 * phase
                    endif
                    if(ind_2.gt.nSites) then
                        ind_2 = ind_2 - nSites
                        if(tAntiPeriodic) phase = -1 * phase
                    endif
                    if(abs(ham(i,j)-(cmplx(real(phase,dp),zero)*ham(ind_1,ind_2))).gt.1.0e-7_dp) then
                        write(6,"(2I6,A,2I6)") i,j," -> ",ind_1,ind_2
                        write(6,*) "Phase change: ",phase
                        write(6,*) ham(i,j),phase*ham(ind_1,ind_2)
                        call writematrixcomp(ham,'real space ham',.true.)
                        call stop_all(t_r,'Translational symmetry not maintained in real-space hamiltonian')
                    endif
                enddo
            enddo
        else
            !2D. Remember that lattice is organised as  1 3
            !                                           2 4
            do i = 1,nImp_x
                do j = 1,nImp_y
                    !Loop through first cell
                    call LatCoordToSiteIndex_2DSquare(i,j,LatInd_UC)
                    !This has coordinates (i,j) and a site index LatInd_UC 
                    do k = 1,nSites
                        !Look at all the coupling from the sites in the first cell, to all others in the system.
                        !This matrix element is (LatInd_UC,k)
                        if(abs(ham(LatInd_UC,k)-conjg(ham(k,LatInd_UC))).gt.1.0e-8_dp) then
                            call stop_all(t_r,'Hamiltonian matrix not hermitian')
                        endif
                        !What is the translation vector to site k?
                        call SiteIndexToLatCoord_2DSquare(k,ind_kx,ind_ky)
                        transvec_x = ind_kx - i
                        transvec_y = ind_ky - j
                        !(transvec_x,transvec_y) gives translation vector from site in unit cell to k
                        Curr_UCx = i !This is the chosen site of the first unit cell
                        Curr_UCy = j
                        Curr_UC = LatInd_UC
                        do p = 1,iImpRepeats-1
                            !Is this the same coupling as the one between the unit cell
                            !p away?
                            !Find the index for (i,j) translated by p unit cells. 
                            !See the Setup2DLattice_Square for how this is done 
                            Curr_UCy = Curr_UCy + nImp_y
                            Curr_UCy = py_mod(Curr_UCy,nSites_y)
                            if(Curr_UCy.eq.0) Curr_UCy = nSites_y
                            call LatCoordToSiteIndex_2DSquare(Curr_UCx,Curr_UCy,Index_2)
                            if(Index_2.lt.Curr_UC) then
                                Curr_UCx = Curr_UCx + nImp_x
                                if(Curr_UCx.gt.nSites_x) call stop_all(t_r,'error')
                                call LatCoordToSiteIndex_2DSquare(Curr_UCx,Curr_UCy,Index_2)
                                if(Index_2.lt.LatInd_UC) call stop_all(t_r,'error')
                                if(Index_2.gt.nSites) call stop_all(t_r,'error')
                            endif
                            Curr_UC = Index_2

                            !Curr_UC index and coordinates (Curr_UCx,Curr_UCy) is (i,j) displaced by p unit cells.
                            !Translate by (transvec_x,transvec_y)
                            call FindDisplacedIndex_2DSquare(Curr_UC,transvec_x,transvec_y,k_trans,PhaseChange)
                            !ham(LatInd_UC,k) should = PhaseChange * ham(Curr_UC,k_trans)
                            if(abs(ham(LatInd_UC,k)-(cmplx(real(PhaseChange,dp),zero)*ham(Curr_UC,k_trans))).gt.1.0e-7_dp) then
                                write(6,"(2I6,A,2I6)") LatInd_UC,k," -> ",Curr_UC,k_trans
                                write(6,*) "Phase change: ",PhaseChange
                                write(6,*) ham(LatInd_UC,k),cmplx(real(PhaseChange,dp),zero)*ham(Curr_UC,k_trans)
                                call writematrixcomp(ham,'real space ham',.true.)
                                call stop_all(t_r,'Translational symmetry not maintained in real-space hamiltonian')
                            endif
                        enddo
                    enddo
                enddo
            enddo
        endif
        write(6,"(A)") "Lattice system appropriately periodic in real space"

    end subroutine CheckRealSpaceTransInvar_z

!Take a real-space hamiltonian on the lattice, and check that it obeys translational symmetry of the unit cells,
!whose size is given by the impurity size. Use through interface for real and imaginary hamiltonians.
    subroutine CheckRealSpaceTransInvar_r(ham)
        implicit none
        real(dp), intent(in) :: ham(nSites,nSites)
        integer :: i,j,ind_1,ind_2,phase,LatInd_UC,k,ind_kx,ind_ky,transvec_x,transvec_y
        integer :: Curr_UCx,Curr_UCy,Curr_UC,p,Index_2,k_trans,PhaseChange
        character(len=*), parameter :: t_r='CheckRealSpaceTransInvar_r'

        !Check for translational invariance of the system
        !That is: is each supercell row equivalent when translated to the supercell ajacent to it?
        if(LatticeDim.eq.1) then
            do j = 1,nSites
                do i = 1,nSites
                    ind_1 = i
                    ind_2 = j
                    !Translate a supercell down
                    ind_1 = i + nImp
                    ind_2 = j + nImp
                    phase = 1
                    if(ind_1.gt.nSites) then
                        ind_1 = ind_1 - nSites
                        if(tAntiPeriodic) phase = -1 * phase
                    endif
                    if(ind_2.gt.nSites) then
                        ind_2 = ind_2 - nSites
                        if(tAntiPeriodic) phase = -1 * phase
                    endif
                    if(abs(ham(i,j)-(phase*ham(ind_1,ind_2))).gt.1.0e-7_dp) then
                        write(6,"(2I6,A,2I6)") i,j," -> ",ind_1,ind_2
                        write(6,*) "Phase change: ",phase
                        write(6,*) ham(i,j),phase*ham(ind_1,ind_2)
                        call writematrix(ham,'real space ham',.true.)
                        call stop_all(t_r,'Translational symmetry not maintained in real-space hamiltonian')
                    endif
                enddo
            enddo
        else
            !2D. Remember that lattice is organised as  1 3
            !                                           2 4
            do i = 1,nImp_x
                do j = 1,nImp_y
                    !Loop through first cell
                    call LatCoordToSiteIndex_2DSquare(i,j,LatInd_UC)
                    !This has coordinates (i,j) and a site index LatInd_UC 
                    do k = 1,nSites
                        !Look at all the coupling from the sites in the first cell, to all others in the system.
                        !This matrix element is (LatInd_UC,k)
                        if(abs(ham(LatInd_UC,k)-ham(k,LatInd_UC)).gt.1.0e-8_dp) then
                            call stop_all(t_r,'Hamiltonian matrix not hermitian')
                        endif
                        !What is the translation vector to site k?
                        call SiteIndexToLatCoord_2DSquare(k,ind_kx,ind_ky)
                        transvec_x = ind_kx - i
                        transvec_y = ind_ky - j
                        !(transvec_x,transvec_y) gives translation vector from site in unit cell to k
                        Curr_UCx = i !This is the chosen site of the first unit cell
                        Curr_UCy = j
                        Curr_UC = LatInd_UC
                        do p = 1,iImpRepeats-1
                            !Is this the same coupling as the one between the unit cell
                            !p away?
                            !Find the index for (i,j) translated by p unit cells. 
                            !See the Setup2DLattice_Square for how this is done 
                            Curr_UCy = Curr_UCy + nImp_y
                            Curr_UCy = py_mod(Curr_UCy,nSites_y)
                            if(Curr_UCy.eq.0) Curr_UCy = nSites_y
                            call LatCoordToSiteIndex_2DSquare(Curr_UCx,Curr_UCy,Index_2)
                            if(Index_2.lt.Curr_UC) then
                                Curr_UCx = Curr_UCx + nImp_x
                                if(Curr_UCx.gt.nSites_x) call stop_all(t_r,'error')
                                call LatCoordToSiteIndex_2DSquare(Curr_UCx,Curr_UCy,Index_2)
                                if(Index_2.lt.LatInd_UC) call stop_all(t_r,'error')
                                if(Index_2.gt.nSites) call stop_all(t_r,'error')
                            endif
                            Curr_UC = Index_2

                            !Curr_UC index and coordinates (Curr_UCx,Curr_UCy) is (i,j) displaced by p unit cells.
                            !Translate by (transvec_x,transvec_y)
                            call FindDisplacedIndex_2DSquare(Curr_UC,transvec_x,transvec_y,k_trans,PhaseChange)
                            !ham(LatInd_UC,k) should = PhaseChange * ham(Curr_UC,k_trans)
                            if(abs(ham(LatInd_UC,k)-(PhaseChange*ham(Curr_UC,k_trans))).gt.1.0e-7_dp) then
                                write(6,"(2I6,A,2I6)") LatInd_UC,k," -> ",Curr_UC,k_trans
                                write(6,*) "Phase change: ",PhaseChange
                                write(6,*) ham(LatInd_UC,k),phase*ham(Curr_UC,k_trans)
                                call writematrix(ham,'real space ham',.true.)
                                call stop_all(t_r,'Translational symmetry not maintained in real-space hamiltonian')
                            endif
                        enddo
                    enddo
                enddo
            enddo
        endif
        write(6,"(A)") "Lattice system appropriately periodic in real space"

    end subroutine CheckRealSpaceTransInvar_r

    !Routine to diagonalize a 1-electron, real operator, with the lattice periodicity.
    !the SS_Period is the size of the supercell repeating unit (e.g. the coupling correlation potential)
    !This will be equivalent to the number of bands per kpoint
    !Ham returns the eigenvectors (in real space).
    subroutine DiagOneEOp_r(Ham,Vals,SS_Period,nLat,tKSpace_Diag,tRealVectors)
        use sort_mod, only: sort_d_a_c
        implicit none
        integer, intent(in) :: nLat
        logical, intent(in) :: tKSpace_Diag
        integer, intent(in) :: SS_Period
        real(dp), intent(inout) :: Ham(nLat,nLat)
        real(dp), intent(out) :: Vals(nLat)
        logical, intent(in), optional :: tRealVectors
        real(dp), allocatable :: Work(:),r_vecs_real(:,:)
        real(dp) :: PrimLattVec(LatticeDim),SiteVec(LatticeDim),Expo,DDOT,phase
        complex(dp), allocatable :: RotMat(:,:),temp(:,:),CompHam(:,:),RotHam(:,:),cWork(:)
        complex(dp), allocatable :: CompHam_2(:,:),ztemp(:,:),k_Ham(:,:),k_vecs(:,:)
        complex(dp), allocatable :: r_vecs(:,:),Vec_temp(:),Vec_temp2(:)
        integer :: lWork,info,i,j,k,kSpace_ind,ki,kj,bandi,bandj,Ind_i,Ind_j,ind_1,ind_2
        integer :: ii,jj,xb,yb,impy,impx,impsite,MappedInd,MappedInd2
        real(dp), allocatable :: Couplings(:,:)
        logical :: tWrapped
        character(len=*), parameter :: t_r='DiagOneEOp_r'
        
        write(6,*) "Entered DiagOneEOp_r..."

        if(present(tRealVectors)) then
            if(.not.tRealVectors) then
                call stop_all(t_r,'Cannot have complex eigenvectors with real hamiltonian!')
            endif
        endif

        if(tKSpace_Diag.and.(LatticeDim.eq.1)) then
            if(tTiltedLattice.and.(LatticeDim.eq.2)) then
                call stop_all(t_r,'Cannot do k-space diagonalizations with tilted lattice '     &
                  //'- impurity site tiling is not same as direct lattice')
            endif

            !Testing of the long range coupling code
!            if(Ham(1,1).gt.0.001_dp) then
!                call writematrix(Ham,'Ham_realspace',.true.)
!                allocate(Couplings(3,SS_Period))
!                do i = 1,3
!                    Couplings(i,:) = -10.0_dp/real(i,dp)
!                enddo
!                call AddPeriodicImpCoupling_RealSpace(Ham,nLat,3,SS_Period,Couplings)
!                deallocate(Couplings)
!                call writematrix(Ham,'Ham_realspace',.true.)
!            endif
            if(tCheck) then
                if(nLat.eq.nSites) call CheckRealSpaceTransInvar(Ham)
            endif

            allocate(CompHam(nLat,nLat))
            do i=1,nLat
                do j=1,nLat
                    CompHam(j,i) = dcmplx(Ham(j,i),zero)
                enddo
            enddo

            allocate(RotMat(nLat,SS_Period))
            allocate(k_Ham(SS_Period,SS_Period))
            allocate(ztemp(nLat,SS_Period))
            allocate(k_vecs(SS_Period,nLat))
            k_vecs(:,:) = zzero
            
            !Space for diagonalization
            lwork = max(1,2*SS_Period-1)
            allocate(cWork(lWork))
            allocate(Work(max(1,3*SS_Period-2)))

            !Run though all kpoints
            do k = 1,nKPnts
                RotMat(:,:) = zzero     !Rot mat will be the rotation into the specfic kpoint of interest
                ind_1 = ((k-1)*SS_Period) + 1
                ind_2 = SS_Period*k
                RotMat(:,:) = RtoK_Rot(:,ind_1:ind_2)
                !We now have the rotation matrix for the bands on this kpoint.
                !Rotate the hamiltonian into this basis
                call ZGEMM('N','N',nLat,SS_Period,nLat,zone,CompHam,nLat,RotMat,nLat,zzero,ztemp,nLat)
                call ZGEMM('C','N',SS_Period,SS_Period,nLat,zone,RotMat,nLat,ztemp,nLat,zzero,k_Ham,SS_Period)

                !Diagonalize this k-pure hamiltonian
                info = 0
                call ZHEEV('V','U',SS_Period,k_Ham,SS_Period,Vals(ind_1:ind_2),cWork,lWork,Work,info)
                if(info.ne.0) call stop_all(t_r,'Diag failed')

                if(tWriteOut) then
                    write(6,*) "For kpoint: ",k," Eigenvalues are:"
                    do i = 0,SS_Period-1
                        write(6,*) Vals(ind_1+i)
                    enddo
                    write(6,*) "Eigenvectors: "
                    call writematrixcomp(k_Ham,'Eigenvec',.true.)
                endif

                k_vecs(:,ind_1:ind_2) = k_Ham(:,:)
            enddo
            
            deallocate(RotMat,k_Ham,ztemp,cWork,Work)

            !Now, rotate the k-space vectors back into r-space
            allocate(r_vecs(nLat,nLat))
            r_vecs(:,:) = zzero
            do k = 1,nKPnts
                ind_1 = ((k-1)*SS_Period) + 1
                ind_2 = SS_Period*k
                call ZGEMM('N','N',nLat,SS_Period,SS_Period,zone,RtoK_Rot(:,ind_1:ind_2),nLat,  &
                    k_vecs(:,ind_1:ind_2),SS_Period,zzero,r_vecs(:,ind_1:ind_2),nLat)
            enddo
            if(tWriteOut) call writematrixcomp(r_vecs,'r_vecs',.true.)

            if(tCheck) then
                !Do these satisfy the original eigenvalue problem?
                allocate(cWork(nLat))
                do i = 1,nLat
!                    write(6,*) "Checking eigenvector: ",i
                    call ZGEMV('N',nLat,nLat,zone,CompHam,nLat,r_vecs(:,i),1,zzero,cWork,1)
                    do j = 1,nLat
                        if(abs(cWork(j)-(Vals(i)*r_vecs(j,i))).gt.1.0e-8_dp) then
                            call stop_all(t_r,'Eigensystem not computed correctly in real basis')
                        endif
                    enddo
                enddo
                deallocate(cWork)
                !write(6,*) "Eigensystem correctly computed and transformed to real space"
            endif

!            write(6,*) "Before sorting..."
!            call writevector(Vals,'Vals')
!            call writematrixcomp(r_vecs,'r_vecs',.true.)
            
            !Order the vectors, such that the are in order of increasing eigenvalue
            call sort_d_a_c(Vals,r_vecs,nSites,nSites)
            
            !write(6,*) "After sorting..."
            !call writevector(Vals,'Vals')
            !call writematrixcomp(r_vecs,'r_vecs',.true.)
            
            if(tCheck) then
                !Do these satisfy the original eigenvalue problem?
                allocate(cWork(nLat))
                do i = 1,nLat
                    call ZGEMV('N',nLat,nLat,zone,CompHam,nLat,r_vecs(:,i),1,zzero,cWork,1)
                    do j = 1,nLat
                        if(abs(cWork(j)-(Vals(i)*r_vecs(j,i))).gt.1.0e-8_dp) then
                            call stop_all(t_r,'Eigensystem not computed correctly in real basis')
                        endif
                    enddo
                enddo
                deallocate(cWork)
                !write(6,*) "Eigensystem correctly computed and transformed to real space"
            endif
            
            !Degenerate sets?
            !It is *impossible* to have non-complex k-pure eigenvectors in the presence of degeneracies.
            !Within the degenerate set, we must rotate them together, losing the k-label. Boo.
            !HACK!!
            !Just rotate them together in equal quantities. We should only ever have mixtures of +-K, which should
            !be taken as the positive and negative LC. e^ik.r = cos k.r + i sin k.r     We want these trig functions
            i = 1
            do while(i.le.(nLat-1))
                j = i + 1
                do while(abs(Vals(i)-Vals(j)).lt.1.0e-8_dp)
                    j = j + 1
                    if(j.gt.nLat) exit
                enddo
                if((j-i).gt.1) then
                    !Degeneracy
                    !Originally, (for 1D systems?), if we had a degeneracy of two, we could assume that 
                    if((j-i).eq.2) then
                        !it was an equal linear combination of the two vectors
!                        if((j-i).gt.2) then
!                            !More than two fold degenerate. How do we rotate these??
!                            call stop_all(t_r,'Cannot handle more than 2x degeneracy in the k-space hamiltonian')
!                        endif
                        allocate(Vec_temp(nLat))
                        allocate(Vec_temp2(nLat))
                        !Degeneracy is between i and i + 1 (or i and j - 1)
                        if((i+1).ne.(j-1)) call stop_all(t_r,'Indexing error')
                        !Rotate these vectors together in positive and negative linear combinations
                        Vec_temp(:) = (r_vecs(:,i) + r_vecs(:,i+1)) / sqrt(2.0_dp)
                        Vec_temp2(:) = (r_vecs(:,i) - r_vecs(:,i+1)) / sqrt(2.0_dp)
                        r_vecs(:,i) = Vec_temp(:)
                        r_vecs(:,i+1) = Vec_temp2(:)
                        deallocate(Vec_temp,Vec_temp2)
                    else
                        !Brute force search for unitary rotation to make these vectors real
                        write(6,"(A,I6,A)") "Degenerate set of: ",j-i,' vectors.'
                        call FindUnitaryRotToRealVecs(r_vecs(:,i:j-1),nLat,j-i)
                    endif
                endif
                i = j 
            enddo

            allocate(r_vecs_real(nLat,nLat))
            r_vecs_real(:,:) = zero
            !Now, find the appropriate phase, such that the rotation will make the r_vecs real.
            !Apply the inverse of this rotation to the k_vecs, such that we end up with a complex set of
            !k-vectors (ordered by k-point), and real set of r_vecs (Ordered by energy).
            do i = 1,nLat   !Run through eigenvectors
                phase = zero
                if(tWriteOut) write(6,*) "Rotating eigenvector : ",i
                do j = 1,nLat
                    if((abs(aimag(r_vecs(j,i))).gt.1.0e-9_dp).and.(abs(r_vecs(j,i)).gt.1.0e-7_dp)) then
                        !Find the phase factor for this eigenvector
                        phase = atan(aimag(r_vecs(j,i))/real(r_vecs(j,i)))
                        if(tWriteOut) write(6,*) "Eigenvector: ",i,j,phase,r_vecs(j,i) * exp(dcmplx(0.0_dp,-phase))
                        exit
                    endif
                enddo
                !The phase should be the same for all components of the eigenvector
                r_vecs(:,i) = r_vecs(:,i) * exp(dcmplx(zero,-phase))
                do j = 1,nLat
                    if(abs(aimag(r_vecs(j,i))).gt.1.0e-6) then
                        write(6,*) "Error rotating component: ",j
                        write(6,*) phase,r_vecs(j,i)
                        call stop_all(t_r,'Eigenvectors not rotated correctly - degeneracies?')
                    endif
                    r_vecs_real(j,i) = real(r_vecs(j,i),dp)
                enddo
            enddo

            !Check again that these rotated r_vecs are correct eigenfunctions...
            if(tCheck) then
                !Do these satisfy the original eigenvalue problem?
                allocate(Work(nLat))
                do i = 1,nLat
                    call DGEMV('N',nLat,nLat,one,Ham,nLat,r_vecs_real(:,i),1,zero,Work,1)
                    do j = 1,nLat
                        if(abs(Work(j)-(Vals(i)*r_vecs_real(j,i))).gt.1.0e-7_dp) then
                            call stop_all(t_r,'Eigensystem not computed correctly in real real basis')
                        endif
                    enddo
                enddo
                deallocate(Work)
                !write(6,*) "Eigensystem correctly computed and transformed to real real space"
            endif

            Ham(:,:) = r_vecs_real(:,:)
            deallocate(CompHam,r_vecs_real)
        else
            !Normal real space diagonalization
            Vals(:) = 0.0_dp
            allocate(Work(1))
            lWork=-1
            info=0
            call dsyev('V','L',nLat,Ham,nLat,Vals,Work,lWork,info)
            if(info.ne.0) call stop_all(t_r,'Workspace query failed')
            lwork=int(work(1))+1
            deallocate(work)
            allocate(work(lwork))
            call dsyev('V','L',nLat,Ham,nLat,Vals,Work,lWork,info)
            if(info.ne.0) call stop_all(t_r,'Diag failed')
            deallocate(work)
!            call writevector(Vals,'Eigenvalues')
        endif

    end subroutine DiagOneEOp_r
                    
    !Find a unitary rotation of the vectors in Vecs, such that they are made to be real
    subroutine FindUnitaryRotToRealVecs(Vecs,nLat,nVecs)
        use MinAlgos
        implicit none
        integer, intent(in) :: nLat,nVecs
        complex(dp), intent(inout) :: Vecs(nLat,nVecs)
        !Local variables
        integer :: nParams
        real(dp), allocatable :: vars(:)    !The independent variables for the X matrix
        real(dp), allocatable :: step(:),var(:)
        complex(dp), allocatable :: UnitMat(:,:),X_Mat(:,:),NewVecs(:,:),VecPair(:,:)
        integer :: maxf,iprint,nloop,iquad,ierr,i,j,ind,iter,k,l
        real(dp) :: stopcr,simp,FinalErr,FuncValue,r,r2,rhoend,rhobeg,initialval,FinalErr_2
        real(dp) :: FinalErr_3,ax,bx,cx,phase
        logical :: tfirst,PosSuc
        logical, allocatable :: tExclude(:)
        logical, parameter :: TwoByTwoRots = .true. 
        logical, parameter :: tSimplex = .false.
        character(len=*), parameter :: t_r='FindUnitaryRotToRealVecs'

        if(TwoByTwoRots) then
            !Just use simplex to do a line-search - silly

            !First, just rotate the orbitals on their own. If they can be made real, then exclude them from the set.
            allocate(tExclude(nVecs))
            tExclude(:) = .false.
            allocate(vars(1))
            do i = 1,nVecs
                nParams = 1
                vars(1) = zero
                !Find bracket
                call MinComplexOrbs(vars,FinalErr,1,nLat,1,Vecs(:,i))
                if(FinalErr.lt.1.0e-7_dp) then
                    tExclude(i) = .true.
                    cycle
                endif
                ax = vars(1)
                vars(1) = 2*pi  !This should be the same value
                call MinComplexOrbs(vars,FinalErr_2,1,nLat,1,Vecs(:,i))
                if(abs(FinalErr-FinalErr_2).gt.1.0e-7_dp) call stop_all(t_r,'Error here')
                vars(1) = zero
                do j = 1,9
                    !Coarse search
                    vars(1) = vars(1) + (2.0_dp*pi)*real(j,dp)/10.0_dp
                    call MinComplexOrbs(vars,FinalErr_3,1,nLat,1,Vecs(:,i))
                    if(FinalErr_3.lt.FinalErr_2) then
                        FinalErr_2 = FinalErr_3
                        bx = vars(1)
                    endif
                enddo
                if(FinalErr_3.lt.FinalErr) then
                    !Bracketed

                    ax = zero
                    cx = 2.0_dp*pi
                    FinalErr = golden(ax,bx,cx,MinComplexOrbs,1.0e-7_dp,vars(1),nLat,Vecs(:,i)) 
!                    call uobyqa_2(nParams,vars,rhobeg,rhoend,iprint,maxf,FinalErr,MinComplexOrbs,nLat,1,Vecs(:,i))
                    if(abs(FinalErr).lt.1.0e-6_dp) then
                        !This vector can be rotated real on its own
                        !Don't rotate it here though - it will be done later on.
                        tExclude(i) = .true.
                    endif
                endif
            enddo

            !Now, attempt to rotate all pairs by a fixed linear combination
            allocate(VecPair(nLat,2))
            do i = 1,nVecs
                if(tExclude(i)) cycle
                do j = i+1,nVecs
                    if(tExclude(j)) cycle
                    VecPair(:,1) = (Vecs(:,i) + Vecs(:,j)) / sqrt(2.0_dp)
                    VecPair(:,2) = (Vecs(:,i) - Vecs(:,j)) / sqrt(2.0_dp)

                    !Can VecPair vectors be made to be real?
                    phase = zero
                    do l = 1,nLat
                        if((abs(aimag(VecPair(l,1))).gt.1.0e-5_dp).and.(abs(VecPair(l,1)).gt.1.0e-5_dp)) then
                            phase = atan(aimag(VecPair(l,1))/real(VecPair(l,1)))
                            exit
                        endif
                    enddo
                    VecPair(:,1) = VecPair(:,1)*exp(cmplx(zero,-phase,dp))
                    !Check the vector is real
                    FinalErr = zero
                    do l = 1,nLat
                        FinalErr = FinalErr + abs(aimag(VecPair(l,1)))
                    enddo
                    if(FinalErr.lt.1.0e-7_dp) then
                        write(6,"(A,2I6)") "Via trivial rotation, two eigenvectors made real: ",i,j
                        Vecs(:,i) = VecPair(:,1)

                        phase = zero
                        do l = 1,nLat
                            if((abs(aimag(VecPair(l,2))).gt.1.0e-5_dp).and.(abs(VecPair(l,2)).gt.1.0e-5_dp)) then 
                                phase = atan(aimag(VecPair(l,2))/real(VecPair(l,2)))
                                exit
                            endif
                        enddo
                        VecPair(:,2) = VecPair(:,2)*exp(cmplx(zero,-phase,dp))
                        FinalErr = zero
                        do l = 1,nLat
                            FinalErr = FinalErr + abs(aimag(VecPair(l,2)))
                        enddo
                        if(FinalErr.gt.1.0e-6_dp) then
                            call stop_all(t_r,'Cannot have one rotation allowed?')
                        endif
                        Vecs(:,j) = VecPair(:,2)
                        tExclude(i) = .true.
                        tExclude(j) = .true.
                    endif

                enddo
            enddo
            deallocate(vars)
                        
            nParams = 3
            allocate(vars(nParams))
            allocate(step(nParams))
            allocate(var(nParams))
            if(tWriteOut) then
                iprint = 50
            else
                iprint = -1
            endif
            iter = 0

            do while(.true.)
                iter = iter + 1
                do i = 1,nVecs
                    if(tExclude(i)) cycle
                    do j = i+1,nVecs
                        if(tExclude(j)) cycle

                        write(6,"(A,2I7)") "Rotating together vectors: ",i,j

                        VecPair(:,1) = Vecs(:,i)
                        VecPair(:,2) = Vecs(:,j)

                        vars(:) = zero  !Initial rotation angle
                        step(:) = 0.05_dp
                        !Fit a quadratic surface to be sure a minimum has been found
                        iquad = 1
                        !As function value is being evaluated in real(dp), it should
                        !be accurate to about 15 dp. If we set simp = 1.d-6, we should
                        !get about 9dp accuracy in fitting the surface
                        simp = 5.0e-7_dp
                        !Set value for stopping criterion. Stopping occurs when the 
                        !standard deviation of the values of the objective function at
                        !the points of the current simplex < stopcr
                        stopcr = 1.0e-7_dp
                        nloop = 2*nParams !This is how often the stopping criterion is checked
                        call MinComplexOrbPair(vars,InitialVal,nParams,nLat,VecPair)
                        write(6,"(A,G25.13)") "Initial function value for this pair: ",InitialVal
                        if(InitialVal.lt.1.0e-7_dp) then
                            write(6,"(A)") "Skipping optimization of this pair - pair are already real!"
                            cycle
                        endif

                        !Set max no of function evaluations. Default = 50*variables, print every 25
                        maxf = 5000*nParams

                        !Now call minim to do the work
                        tfirst = .true.
                        do while(.true.)
                            call minim_3(vars, step, nParams, FinalErr, maxf, iprint, stopcr, nloop, &
                                iquad, simp, var, MinComplexOrbPair, ierr, nLat, VecPair)

                            if(ierr.eq.0) exit
                            if(.not.tFirst) exit
                            tFirst = .false.    !We have found a minimum, but run again with a small number of iterations to check it is stable
                            maxf = 3*nParams
                        enddo

                        !On output, Err is the minimized objective function, var contains the diagonal of the inverse of the information matrix, whatever that is
                        if(ierr.eq.0) then
                            write(6,"(A)") "Simplex optimization successful."
                            write(6,"(A,F14.7)") "Minimized residual: ",FinalErr
                        elseif(ierr.eq.4) then
                            call stop_all(t_r,'nloop < 1')
                        elseif(ierr.eq.3) then
                            call stop_all(t_r,'nParams < 1')
                        elseif(ierr.eq.2) then
                            call warning(t_r,'Information matrix is not +ve semi-definite')
                            write(6,"(A,F14.7)") "Final residual: ",FinalErr
                        elseif(ierr.eq.1) then
                            call warning(t_r,'Max number of Simplex function evaluations reached.')
                            write(6,"(A,F14.7)") "Final residual: ",FinalErr
                        endif

                        !Update vectors
                        if(((InitialVal-FinalErr).gt.(InitialVal*0.1_dp)).or.((InitialVal-FinalErr).gt.0.2)) then
                            !Only update if change is more than 10% of initial error
                            write(6,"(A)") "Updating vector pair."
                            write(6,"(A,3F15.7)") "Rotation angle: ",vars(1),cos(vars(1)),sin(vars(1))
                            Vecs(:,i) = exp(cmplx(zero,vars(2),dp))*(cos(vars(1))*VecPair(:,1) + sin(vars(1))*VecPair(:,2))
                            Vecs(:,j) = exp(cmplx(zero,vars(3),dp))*(-sin(vars(1))*VecPair(:,1) + cos(vars(1))*VecPair(:,2))
                        else
                            write(6,"(A)") "Function has not decreased materially. Not updating."
                        endif

                        !Check if each vector is real. If it is, then add to the exclude list.
                        FinalErr = zero
                        FinalErr_2 = zero
                        do k = 1,nLat
                            FinalErr = FinalErr + abs(aimag(Vecs(k,i)))
                            FinalErr_2 = FinalErr_2 + abs(aimag(Vecs(k,j)))
                        enddo
                        if(FinalErr.lt.1.0e-6) tExclude(i) = .true.
                        if(FinalErr_2.lt.1.0e-6) tExclude(j) = .true.

                    enddo
                enddo
                !Check for convergence
                !Now, how complex are these vectors?
                FuncValue = zero
                do i = 1,nVecs
                    do j = 1,nLat
                        FuncValue = FuncValue + abs(aimag(Vecs(j,i)))
                    enddo
                enddo
                if(FuncValue.gt.(1.0e-6_dp*nVecs)) then
!                    write(6,*) "Sum of imaginary components of degenerate orbitals: ",FuncValue
                    write(6,*) "After macroiteration ",iter," objective function is: ",FuncValue
                else
                    exit    !Convergence successful
                endif

            enddo

            deallocate(step,var,VecPair,vars,tExclude)

        else
            !Try to optimize full unitary matrix all in one go

            !How many free parameters do we have?
            nParams = nVecs*(nVecs-1) + nVecs

            allocate(vars(nParams))
            !Starting values for anti-hermitian rotation matrix
            vars(:) = zero
            maxf = 5000*nParams

            if(tSimplex) then
                allocate(step(nParams))
                allocate(var(nParams))
                step(:) = 0.05_dp
                !Set max no of function evaluations. Default = 50*variables, print every 25
                if(tWriteOut) then
                    iprint = 50
                else
                    iprint = -1
                endif
                !Set value for stopping criterion. Stopping occurs when the 
                !standard deviation of the values of the objective function at
                !the points of the current simplex < stopcr
                stopcr = 1.0e-6_dp
                nloop = 2*nParams !This is how often the stopping criterion is checked
                !Fit a quadratic surface to be sure a minimum has been found
                iquad = 1
                !As function value is being evaluated in real(dp), it should
                !be accurate to about 15 dp. If we set simp = 1.d-6, we should
                !get about 9dp accuracy in fitting the surface
                simp = 1.0e-6_dp
            else
                !Brent
                rhobeg = 0.05_dp
                rhoend = 1.0e-6_dp
                iprint = 2
            endif

            iter = 0
            do while(.true.)
                iter = iter + 1
                    
                if(tSimplex) then
                    !Now call minim to do the work
                    tfirst = .true.
                    do while(.true.)
                        call minim_2(vars, step, nParams, FinalErr, maxf, iprint, stopcr, nloop, &
                            iquad, simp, var, MinComplexOrbs, ierr, nLat, nVecs, Vecs)

                        if(ierr.eq.0) exit
                        if(.not.tFirst) exit
                        tFirst = .false.    !We have found a minimum, but run again with a small number of iterations to check it is stable
                        maxf = 3*nParams
                    enddo
                else
                    call uobyqa_2(nParams,vars,rhobeg,rhoend,iprint,maxf,FinalErr,MinComplexOrbs,nLat,nVecs,Vecs)
                    ierr = 0
                endif

                if(FinalErr.lt.1.0e-6_dp) then
                    !Converged
                    !On output, Err is the minimized objective function, var contains the diagonal of the inverse of the information matrix, whatever that is
                    if(ierr.eq.0) then
                        write(6,"(A)") "Minimization successful."
                        write(6,"(A,F14.7)") "Minimized residual: ",FinalErr
                    elseif(ierr.eq.4) then
                        call stop_all(t_r,'nloop < 1')
                    elseif(ierr.eq.3) then
                        call stop_all(t_r,'nParams < 1')
                    elseif(ierr.eq.2) then
                        call warning(t_r,'Information matrix is not +ve semi-definite')
                        write(6,"(A,F14.7)") "Final residual: ",FinalErr
                    elseif(ierr.eq.1) then
                        call warning(t_r,'Max number of Simplex function evaluations reached.')
                        write(6,"(A,F14.7)") "Final residual: ",FinalErr
                    endif

                    !Now, do a final rotation of the orbitals, and check that they are real
                    allocate(X_Mat(nVecs,nVecs))
                    X_Mat(:,:) = zzero

                    do i = 1,nVecs
                        !The first values are the real diagonals
                        X_Mat(i,i) = cmplx(vars(i),zero,dp)
                    enddo

                    !Now, for the hermitian off-diagonals
                    ind = nVecs+1
                    do i = 1,nVecs
                        do j = i+1,nVecs
                            X_Mat(j,i) = cmplx(vars(ind),vars(ind+1),dp)
                            X_Mat(i,j) = conjg(X_Mat(j,i))
                            ind = ind + 2
                        enddo
                    enddo
                    if((ind-1).ne.nParams) call stop_all(t_r,'Error in indexing')
                    allocate(UnitMat(nVecs,nVecs))
                    call FindUnitaryMatrix(nVecs,X_Mat,UnitMat)
                    
                    !Now, rotate the vectors into a new basis!
                    allocate(NewVecs(nLat,nVecs))
                    call ZGEMM('N','N',nLat,nVecs,nVecs,zone,Vecs,nLat,UnitMat,nVecs,zzero,NewVecs,nLat)

                    !Now, how complex are these vectors?
                    FuncValue = zero
                    do i = 1,nVecs
                        do j = 1,nLat
                            FuncValue = FuncValue + abs(aimag(NewVecs(j,i)))
                        enddo
                    enddo
                    if(FuncValue.gt.1.0e-6_dp) then
                        write(6,*) "Sum of imaginary components of degenerate orbitals: ",FuncValue
                        call stop_all(t_r,'Orbitals still complex')
                    endif
                    !Update vectors
                    Vecs(:,:) = NewVecs(:,:)

                    deallocate(NewVecs,X_Mat,UnitMat)
                    exit
                else
                    !We know that strictly real eigenvectors are possible. Start from a new random matrix and search again!
                    write(6,*) "Restarting from random unitary matrix"
                    do i = 1,nParams
                        call random_number(r)
                        if(r.gt.0.5_dp) then
                            call random_number(r2)
                            vars(i) = r2*pi
                        else
                            call random_number(r2)
                            vars(i) = -r2*pi
                        endif
                    enddo
                endif
            enddo
            if(tSimplex) deallocate(step,var)
            deallocate(vars)
        endif

    end subroutine FindUnitaryRotToRealVecs
    
    !Given the parameters for a complex anti-hermitian matrix, calculate the imaginary component of the
    !resultant orbitals once rotated into the basis via the unitary transform.
    subroutine MinComplexOrbPair(vars,FuncValue,nParams,nLat,VecPair)
        implicit none
        integer, intent(in) :: nLat,nParams
        real(dp), intent(in) :: vars(nParams)
        real(dp), intent(out) :: FuncValue
        complex(dp), intent(in) :: VecPair(nLat,2)
        !Local variables
        complex(dp), allocatable :: NewVecs(:,:)
        integer :: i,j
        integer, parameter :: iPwr = 1
        character(len=*), parameter :: t_r='MinComplexOrbPair'

        if(nParams.ne.3) call stop_all(t_r,'Error')

        !Now, rotate the vectors into a new basis!
        allocate(NewVecs(nLat,2))

!        write(6,*) vars(1),vars(2),exp(cmplx(zero,vars(2),dp))
        NewVecs(:,1) = exp(cmplx(zero,vars(2),dp))*(cos(vars(1))*VecPair(:,1) + sin(vars(1))*VecPair(:,2))
        NewVecs(:,2) = exp(cmplx(zero,vars(3),dp))*(-sin(vars(1))*VecPair(:,1) + cos(vars(1))*VecPair(:,2))

        !Now, how complex are these vectors?
        FuncValue = zero
        do i = 1,2    
            do j = 1,nLat
                FuncValue = FuncValue + abs(aimag(NewVecs(j,i)))**iPwr
            enddo
        enddo

        deallocate(NewVecs)

    end subroutine MinComplexOrbPair

    !Given the parameters for a complex anti-hermitian matrix, calculate the imaginary component of the
    !resultant orbitals once rotated into the basis via the unitary transform.
    subroutine MinComplexOrbs(vars,FuncValue,nParams,nLat,nVecs,Vecs)
        implicit none
        integer, intent(in) :: nLat,nVecs,nParams
        real(dp), intent(in) :: vars(nParams)
        real(dp), intent(out) :: FuncValue
        complex(dp), intent(in) :: Vecs(nLat,nVecs)
        !Local variables
        complex(dp), allocatable :: X_Mat(:,:),UnitMat(:,:)
        complex(dp), allocatable :: NewVecs(:,:)
        integer :: i,ind,j,b
        integer, parameter :: iPwr = 1
        character(len=*), parameter :: t_r='MinComplexOrbs'

        if(nVecs.eq.1) then
            !Just a global phase factor
            if(nParams.ne.1) call stop_all(t_r,'error')
            FuncValue = zero
            do j = 1,nLat
                FuncValue = FuncValue + abs(aimag(exp(cmplx(zero,-vars(1),dp))*Vecs(j,1)))**iPwr
            enddo
        else
            allocate(X_Mat(nVecs,nVecs))
            X_Mat(:,:) = zzero

            do i = 1,nVecs
                !The first values are the real diagonals
                X_Mat(i,i) = cmplx(vars(i),zero,dp)
            enddo

            !Now, for the hermitian off-diagonals
            ind = nVecs+1
            do i = 1,nVecs
                do j = i+1,nVecs
                    X_Mat(j,i) = cmplx(vars(ind),vars(ind+1),dp)
                    X_Mat(i,j) = conjg(X_Mat(j,i))
                    ind = ind + 2
                enddo
            enddo
            if((ind-1).ne.nParams) call stop_all(t_r,'Error in indexing')

            allocate(UnitMat(nVecs,nVecs))
            call FindUnitaryMatrix(nVecs,X_Mat,UnitMat)

            !Now, rotate the vectors into a new basis!
            allocate(NewVecs(nLat,nVecs))
            call ZGEMM('N','N',nLat,nVecs,nVecs,zone,Vecs,nLat,UnitMat,nVecs,zzero,NewVecs,nLat)

            !Now, how complex are these vectors?
            FuncValue = zero
            do i = 1,nVecs
                do j = 1,nLat
                    FuncValue = FuncValue + abs(aimag(NewVecs(j,i)))**iPwr
                enddo
            enddo

            deallocate(NewVecs,X_Mat,UnitMat)
        endif

    end subroutine MinComplexOrbs
        
    !Find the unitary matrix resulting from a hermitian matrix X_Mat as exp(-iX)
    subroutine FindUnitaryMatrix(nVecs,X_Mat,UnitMat)
        implicit none
        integer, intent(in) :: nVecs
        complex(dp), intent(in) :: X_Mat(nVecs,nVecs)
        complex(dp), intent(out) :: UnitMat(nVecs,nVecs)
        complex(dp), allocatable :: EigenVecs(:,:)
        integer :: i,j,b,lWork,info
        complex(dp), allocatable :: cWork(:),ztemp(:,:)
        real(dp), allocatable :: EValues(:),Work(:)
        character(len=*), parameter :: t_r='FindUnitaryMatrix'

        if(tCheck) then
            !Check it is hermitian
            do i = 1,nVecs
                do j = 1,nVecs
                    if(abs(X_Mat(i,j)-conjg(X_Mat(j,i))).gt.1.0e-8_dp) then
                        call stop_all(t_r,'Matrix not hermitian')
                    endif
                enddo
                if(abs(aimag(X_Mat(i,i))).gt.1.0e-8_dp) then
                    call stop_all(t_r,'Matrix not hermitian')
                endif
            enddo
        endif

        !Diagonalize
        allocate(EigenVecs(nVecs,nVecs))
        EigenVecs(:,:) = X_Mat(:,:)
        allocate(EValues(nVecs))
        EValues(:) = zero
        allocate(cWork(1))
        allocate(Work(max(1,3*nVecs-2)))
        lWork = -1
        info = 0
!        write(6,*) "nVecs: ",nVecs
!        call writematrixcomp(EigenVecs,'X_Mat',.false.)
!        write(6,*) "EValues: ",EValues
        call zheev('V','U',nVecs,EigenVecs,nVecs,EValues,cWork,lWork,Work,info)
        if(info.ne.0) call stop_all(t_r,'Workspace query failed')
        lwork=int(abs(work(1)))+2 
        if(lwork.lt.max(1,2*nVecs-1)) lwork = max(1,2*nVecs-1)
        deallocate(cwork)
        allocate(cwork(lwork))
!        write(6,*) "lWork: ",lWork
!        write(6,*) "info: ",info
        call zheev('V','U',nVecs,EigenVecs,nVecs,EValues,cWork,lWork,Work,info)
        if(info.ne.0) then
            write(6,*) "info: ", info
            call stop_all(t_r,'Diag failed')
        endif
        deallocate(work,cwork)

        !Create diagonal complex matrix with -i*evalues on diagonal
        !Multiply this straight away by V*
        do b = 1,nVecs
            do i = 1,nVecs
                UnitMat(i,b) = exp(cmplx(zero,-EValues(i),dp))*conjg(EigenVecs(b,i))
            enddo
        enddo
        !Now the final rotation back into the original space
        deallocate(EValues)
        allocate(ztemp(nVecs,nVecs))
        call ZGEMM('N','N',nVecs,nVecs,nVecs,zone,EigenVecs,nVecs,UnitMat,nVecs,zzero,ztemp,nVecs)
        UnitMat(:,:) = ztemp(:,:)

        if(tCheck) then
            !Verify that this matrix is unitary
            call ZGEMM('C','N',nVecs,nVecs,nVecs,zone,UnitMat,nVecs,UnitMat,nVecs,zzero,ztemp,nVecs)
            !ztemp should be the unit matrix
            do i = 1,nVecs
                ztemp(i,i) = ztemp(i,i) - zone
            enddo
            !Check this is now the zero matrix
            do i = 1,nVecs
                do j = 1,nVecs
                    if(abs(ztemp(j,i)).gt.1.0e-7_dp) then
                        call stop_all(t_r,'Not created unitary matrix')
                    endif
                enddo
            enddo
            !For completeness, check the other way around
            call ZGEMM('N','C',nVecs,nVecs,nVecs,zone,UnitMat,nVecs,UnitMat,nVecs,zzero,ztemp,nVecs)
            !ztemp should be the unit matrix
            do i = 1,nVecs
                ztemp(i,i) = ztemp(i,i) - zone
            enddo
            !Check this is now the zero matrix
            do i = 1,nVecs
                do j = 1,nVecs
                    if(abs(ztemp(j,i)).gt.1.0e-7_dp) then
                        call stop_all(t_r,'Not created unitary matrix 2')
                    endif
                enddo
            enddo
        endif
        deallocate(ztemp,EigenVecs)
    end subroutine FindUnitaryMatrix
    
    !Routine to diagonalize a 1-electron, complex operator, with the lattice periodicity.
    !the SS_Period is the size of the supercell repeating unit (e.g. the coupling correlation potential)
    !This will be equivalent to the number of bands per kpoint
    !Ham returns the eigenvectors (in real space).
    subroutine DiagOneEOp_z(Ham,Vals,SS_Period,nLat,tKSpace_Diag,tRealVectors)
        use sort_mod, only: sort_d_a_c,sort_d_a_c_a_r
        implicit none
        integer, intent(in) :: nLat
        logical, intent(in) :: tKSpace_Diag
        integer, intent(in) :: SS_Period
        complex(dp), intent(inout) :: Ham(nLat,nLat)
        real(dp), intent(out) :: Vals(nLat)
        logical, intent(in), optional :: tRealVectors
        real(dp), allocatable :: Work(:),r_vecs_real(:,:),kPnts_tmp(:,:)
        real(dp) :: phase,kPnt_(LatticeDim),kPnt_1(LatticeDim),kPnt_2(LatticeDim),kPnt_3(LatticeDim)
        complex(dp), allocatable :: RotMat(:,:),CompHam(:,:),cWork(:),vec_temp(:)
        complex(dp), allocatable :: ztemp(:,:),k_Ham(:,:),k_vecs(:,:),vec_temp2(:)
        complex(dp), allocatable :: r_vecs(:,:)
        integer, allocatable :: PickedDegenerateOrbs(:)
        integer :: lWork,info,i,j,k,l,kSpace_ind,ind_1,ind_2,LowerDegenInd,UpperDegenInd
        integer :: nDegenOcc,nDegenOrbs,nDegenKPnts,ind,ind_found
        logical :: tRealVectors_,tSame
        logical :: tNetZeroMomDet
        character(len=*), parameter :: t_r='DiagOneEOp_z'

        write(6,*) "Entered DiagOneEOp_z..."

        if(present(tRealVectors)) then
            tRealVectors_ = tRealVectors
        else
            tRealVectors_ = .true. 
        endif

        if(tKSpace_Diag) then
            if(tRealVectors_.and.(LatticeDim.eq.2)) then
                call stop_all(t_r,'Cannot provide real eigenvectors for 2D systems. Sorry!')
            endif
            if(tTiltedLattice.and.(LatticeDim.eq.2)) then
                call stop_all(t_r,'Cannot do k-space diagonalizations - impurity site tiling is not same as direct lattice')
            endif
            if((LatticeDim.eq.2).and.(tShift_Mesh)) then
                tNetZeroMomDet = .true.
            else
                tNetZeroMomDet = .false.
            endif

            if(tCheck) then
                if(nLat.eq.nSites) call CheckRealSpaceTransInvar(Ham)
            endif

            allocate(CompHam(nLat,nLat))
            CompHam(:,:) = Ham(:,:)

            allocate(RotMat(nLat,SS_Period))
            allocate(k_Ham(SS_Period,SS_Period))
            allocate(ztemp(nLat,SS_Period))
            allocate(k_vecs(SS_Period,nLat))
            k_vecs(:,:) = zzero
            
            !Space for diagonalization
            lwork = max(1,2*SS_Period-1)
            allocate(cWork(lWork))
            allocate(Work(max(1,3*SS_Period-2)))

            !Run though all kpoints
            do k = 1,nKPnts
                RotMat(:,:) = zzero     !Rot mat will be the rotation into the specfic kpoint of interest
                ind_1 = ((k-1)*SS_Period) + 1
                ind_2 = SS_Period*k
                RotMat(:,:) = RtoK_Rot(:,ind_1:ind_2)
                !We now have the rotation matrix for the bands on this kpoint.
                !Rotate the hamiltonian into this basis
                call ZGEMM('N','N',nLat,SS_Period,nLat,zone,CompHam,nLat,RotMat,nLat,zzero,ztemp,nLat)
                call ZGEMM('C','N',SS_Period,SS_Period,nLat,zone,RotMat,nLat,ztemp,nLat,zzero,k_Ham,SS_Period)

                !Diagonalize this k-pure hamiltonian
                info = 0
                call ZHEEV('V','U',SS_Period,k_Ham,SS_Period,Vals(ind_1:ind_2),cWork,lWork,Work,info)
                if(info.ne.0) call stop_all(t_r,'Diag failed')

                if(tWriteOut) then
                    write(6,*) "For kpoint: ",k," Eigenvalues are:"
                    do i = 0,SS_Period-1
                        write(6,*) Vals(ind_1+i)
                    enddo
                    write(6,*) "Eigenvectors: "
                    call writematrixcomp(k_Ham,'Eigenvec',.true.)
                endif

                k_vecs(:,ind_1:ind_2) = k_Ham(:,:)
            enddo
            
            deallocate(RotMat,k_Ham,ztemp,cWork,Work)

            !Now, rotate the k-space vectors back into r-space
            allocate(r_vecs(nLat,nLat))
            r_vecs(:,:) = zzero
            do k = 1,nKPnts
                ind_1 = ((k-1)*SS_Period) + 1
                ind_2 = SS_Period*k
                call ZGEMM('N','N',nLat,SS_Period,SS_Period,zone,RtoK_Rot(:,ind_1:ind_2),nLat,  &
                    k_vecs(:,ind_1:ind_2),SS_Period,zzero,r_vecs(:,ind_1:ind_2),nLat)
            enddo
            if(tWriteOut) call writematrixcomp(r_vecs,'r_vecs',.true.)

            if(tCheck) then
                !Do these satisfy the original eigenvalue problem?
                allocate(cWork(nLat))
                do i = 1,nLat
!                    write(6,*) "Checking eigenvector: ",i
                    call ZGEMV('N',nLat,nLat,zone,CompHam,nLat,r_vecs(:,i),1,zzero,cWork,1)
                    do j = 1,nLat
                        if(abs(cWork(j)-(Vals(i)*r_vecs(j,i))).gt.1.0e-8_dp) then
                            call stop_all(t_r,'Eigensystem not computed correctly in real basis')
                        endif
                    enddo
                enddo
                deallocate(cWork)
                !write(6,*) "Eigensystem correctly computed and transformed to real space"
            endif

!            write(6,*) "Before sorting..."
!            call writevector(Vals,'Vals')
!            call writematrixcomp(r_vecs,'r_vecs',.true.)
            
            if(tNetZeroMomDet) then
                !Order according to energy, *but* keep the kpoint label!
                if((LatticeDim.ne.2).and.tShift_Mesh) call stop_all(t_r,'Error here')
                allocate(kPnts_tmp(LatticeDim,nKPnts))
                KPnts_tmp(:,:) = KPnts(:,:) 
                call sort_d_a_c_a_r(Vals,r_vecs,kPnts_tmp,nSites,nSites,LatticeDim)

                !Is there a degeneracy around nOcc? If so, we need to ensure that we order the kpoints such that
                !the net momentum of the occupied orbitals in the degenerate set = 0.
                j = nOcc + 1
                do while(abs(Vals(nOcc)-Vals(j)).lt.1.0e-8_dp)
                    j = j + 1
                    if(j.gt.nLat) exit
                enddo
                UpperDegenInd = j - 1
                j = nOcc - 1
                do while(abs(Vals(nOcc)-Vals(j)).lt.1.0e-8_dp)
                    j = j - 1
                    if(j.lt.1) exit
                enddo
                LowerDegenInd = j + 1
                if(UpperDegenInd-nOcc.gt.0) then
                    !Unfortunately, we have a degeneracy at the fermi surface.
                    !Choose a determinant which has a net zero momentum, and hope for the best.
                    !At least then it will have the right density.
                    nDegenOcc = nOcc - LowerDegenInd + 1            !Number of occupied pairs to have
                    nDegenOrbs = UpperDegenInd - LowerDegenInd + 1  !Number of degenerate orbitals
                    !Count the number of unique kpoints
                    nDegenKPnts = 0
                    do i = LowerDegenInd,UpperDegenInd
                        tSame = .false.
                        do j = LowerDegenInd,i-1
                            tSame = .true.
                            do l = 1,LatticeDim
                                tSame = tSame.and.(abs(KPnts_tmp(l,i)-KPnts_tmp(l,j)).lt.1.0e-8_dp)
                            enddo
                            if(tSame) exit
                        enddo
                        !If tSame then we have found a kpoint already sampled
                        if(.not.tSame) nDegenKPnts = nDegenKPnts + 1
                    enddo
                    write(6,*) "Number of degenerate kpoints at the fermi surface = ",nDegenKPnts

                    allocate(PickedDegenerateOrbs(nDegenOcc))
                    PickedDegenerateOrbs(:) = 0 !The index of the orbitals to pick

                    if(tShift_Mesh) then
                        !We should always have all four equivalent kpoints occupied.
                        !Therefore, just find them all, such that the net contribution is always zero
                        ind = LowerDegenInd - 1
                        ind_found = 0
                        do while(ind_found.lt.nDegenOcc)
                            ind = ind + 1
                            if(ind.gt.UpperDegenInd) call stop_all(t_r,'Indexing problem')
                            kPnt_(:) = KPnts_tmp(:,ind)
                            ind_found = ind_found + 1
                            if(ind_found.gt.nDegenOcc) call stop_all(t_r,'Indexing error')
                            PickedDegenerateOrbs(ind_found) = ind

                            !Now, find the three other equivalent degenerate orbitals to this one
                            kPnt_1(1) = kPnt_(2)
                            kPnt_1(2) = -kPnt_(1)
                            !Search for a degenerate orbital at this kpoint, ensuring that it is not an orbital we have picked already
                            do j = LowerDegenInd,UpperDegenInd
                                tSame = .true.
                                do l = 1,LatticeDim
                                    tSame = tSame.and.(abs(KPnts_tmp(l,j)-KPnt_1(l)).lt.1.0e-8_dp)
                                enddo
                                if(tSame) then
                                    !Now check it is not an orbital we have already picked
                                    do l = 1,ind_found
                                        if(PickedDegenerateOrbs(l).eq.j) then
                                            tSame = .false.
                                            exit
                                        endif
                                    enddo
                                endif
                                if(tSame) exit
                            enddo
                            if(.not.tSame) call stop_all(t_r,'Could not find equivalent kpoint')
                            ind_found = ind_found + 1
                            if(ind_found.gt.nDegenOcc) call stop_all(t_r,'Indexing error')
                            PickedDegenerateOrbs(ind_found) = j
                            
                            kPnt_2(1) = -kPnt_(2)
                            kPnt_2(2) = kPnt_(1)
                            !Search for a degenerate orbital at this kpoint, ensuring that it is not an orbital we have picked already
                            do j = LowerDegenInd,UpperDegenInd
                                tSame = .true.
                                do l = 1,LatticeDim
                                    tSame = tSame.and.(abs(KPnts_tmp(l,j)-KPnt_2(l)).lt.1.0e-8_dp)
                                enddo
                                if(tSame) then
                                    !Now check it is not an orbital we have already picked
                                    do l = 1,ind_found
                                        if(PickedDegenerateOrbs(l).eq.j) then
                                            tSame = .false.
                                            exit
                                        endif
                                    enddo
                                endif
                                if(tSame) exit
                            enddo
                            if(.not.tSame) call stop_all(t_r,'Could not find equivalent kpoint')
                            ind_found = ind_found + 1
                            if(ind_found.gt.nDegenOcc) call stop_all(t_r,'Indexing error')
                            PickedDegenerateOrbs(ind_found) = j

                            kPnt_3(1) = -kPnt_(1)
                            kPnt_3(2) = -kPnt_(2)
                            !Search for a degenerate orbital at this kpoint, ensuring that it is not an orbital we have picked already
                            do j = LowerDegenInd,UpperDegenInd
                                tSame = .true.
                                do l = 1,LatticeDim
                                    tSame = tSame.and.(abs(KPnts_tmp(l,j)-KPnt_3(l)).lt.1.0e-8_dp)
                                enddo
                                if(tSame) then
                                    !Now check it is not an orbital we have already picked
                                    do l = 1,ind_found
                                        if(PickedDegenerateOrbs(l).eq.j) then
                                            tSame = .false.
                                            exit
                                        endif
                                    enddo
                                endif
                                if(tSame) exit
                            enddo
                            if(.not.tSame) call stop_all(t_r,'Could not find equivalent kpoint')
                            ind_found = ind_found + 1
                            if(ind_found.gt.nDegenOcc) call stop_all(t_r,'Indexing error')
                            PickedDegenerateOrbs(ind_found) = j

                            kPnt_(:) = kPnt_(:) + kPnt_1(:) + kPnt_2(:) + kPnt_3(:)
                            do l = 1,LatticeDim
                                if(abs(kPnt_(l)).gt.1.0e-8_dp) call stop_all(t_r,'Not found translationally symmetric set')
                            enddo

                        enddo

                        do i = 1,nDegenOcc
                            do j = 1,nDegenOcc
                                if((i.ne.j).and.(PickedDegenerateOrbs(i).eq.PickedDegenerateOrbs(j))) then
                                    call stop_all(t_r,'Error here')
                                endif
                            enddo
                        enddo

!                        write(6,*) "List of orbitals and kpoints: "
!                        do i = 1,nSites
!                            write(6,"(I5,3F20.10)") i,Vals(i),kPnts_tmp(:,i)/pi
!                        enddo
!                        write(6,*) "Reordering vector: "
!                        do i = 1,nDegenOcc
!                            write(6,*) i,PickedDegenerateOrbs(i)
!                        enddo

                        !Now, we have a list of PickedDegenerateOrbs, of length nDegenOcc.
                        !These want to be the first ones in the degenerate set.
                        allocate(cWork(nSites))
                        do i = LowerDegenInd,LowerDegenInd+nDegenOcc-1
                            !Swap index PickedDegenOrbs(i) with i
                            !write(6,"(A,I8,A,I5)") "Swapping orbital ",PickedDegenerateOrbs(i)," and ",i
!                            write(6,*) "Swapping orbital ",PickedDegenerateOrbs(i-LowerDegenInd+1)," and ",i
                            cWork(:) = r_vecs(:,i)
                            phase = Vals(i)
                            kPnt_(:) = kPnts_tmp(:,i)
                            r_vecs(:,i) = r_vecs(:,PickedDegenerateOrbs(i-LowerDegenInd+1))
                            Vals(i) = Vals(PickedDegenerateOrbs(i-LowerDegenInd+1))
                            kPnts_tmp(:,i) = kPnts_tmp(:,PickedDegenerateOrbs(i-LowerDegenInd+1))
                            r_vecs(:,PickedDegenerateOrbs(i-LowerDegenInd+1)) = cWork(:)
                            Vals(PickedDegenerateOrbs(i-LowerDegenInd+1)) = phase
                            kPnts_tmp(:,PickedDegenerateOrbs(i-LowerDegenInd+1)) = kPnt_(:)
                            do j = i-LowerDegenInd+2,nDegenOcc
                                !Run through rest of swaps
                                if(PickedDegenerateOrbs(j).eq.i) then
                                    !i has swapped with PickedDegenerateOrbs(i)
!                                    write(6,*) "Slot ",j," changing to pick ",PickedDegenerateOrbs(i-LowerDegenInd+1)," rather than ",PickedDegenerateOrbs(j)
                                    PickedDegenerateOrbs(j) = PickedDegenerateOrbs(i-LowerDegenInd+1)
                                    exit
                                endif
                            enddo
!                            write(6,*) "After first swap, degenerate space is: "
!                            write(6,*) "List of orbitals and kpoints: "
!                            do j = LowerDegenInd,UpperDegenInd
!                                write(6,"(I5,3F20.10)") j,Vals(j),kPnts_tmp(:,j)/pi
!                            enddo
                        enddo
                        deallocate(cWork)
                    else
                        call stop_all(t_r,'Must be using Monkhorst-Pack mesh')
                    endif   !tShift_Mesh

                    deallocate(PickedDegenerateOrbs)
                endif   !Open shell at Fermi surface
                kPnt_(:) = zero
                do i = 1,nOcc
                    kPnt_(:) = kPnt_(:) + kPnts_tmp(:,i)
                enddo
                do l = 1,LatticeDim
                    if(abs(kPnt_(l)).gt.1.0e-8_dp) call stop_all(t_r,'Non zero momentum determinant chosen')
                enddo
                kPnt_(:) = zero
                do i = nOcc+1,nSites
                    kPnt_(:) = kPnt_(:) + kPnts_tmp(:,i)
                enddo
                do l = 1,LatticeDim
                    if(abs(kPnt_(l)).gt.1.0e-8_dp) call stop_all(t_r,'Non zero momentum virtual space')
                enddo
                deallocate(kPnts_tmp)
            else
                !Order the vectors, such that the are in order of increasing eigenvalue
                call sort_d_a_c(Vals,r_vecs,nSites,nSites)
            endif

            
            !write(6,*) "After sorting..."
            !call writevector(Vals,'Vals')
            !call writematrixcomp(r_vecs,'r_vecs',.true.)
            
            if(tCheck) then
                !Do these satisfy the original eigenvalue problem?
                allocate(cWork(nLat))
                do i = 1,nLat
                    call ZGEMV('N',nLat,nLat,zone,CompHam,nLat,r_vecs(:,i),1,zzero,cWork,1)
                    do j = 1,nLat
                        if(abs(cWork(j)-(Vals(i)*r_vecs(j,i))).gt.1.0e-8_dp) then
                            call stop_all(t_r,'Eigensystem not computed correctly in real basis')
                        endif
                    enddo
                enddo
                deallocate(cWork)
                !write(6,*) "Eigensystem correctly computed and transformed to real space"
            endif

            if(tRealVectors_) then
                !Degenerate sets?
                !It is *impossible* to have non-complex k-pure eigenvectors in the presence of degeneracies.
                !Within the degenerate set, we must rotate them together, losing the k-label. Boo.
                !HACK!!
                !Just rotate them together in equal quantities. We should only ever have mixtures of +-K, which should
                !be taken as the positive and negative LC. e^ik.r = cos k.r + i sin k.r     We want these trig functions
                allocate(Vec_temp(nLat))
                allocate(Vec_temp2(nLat))
                i = 1
                do while(i.le.(nLat-1))
                    j = i + 1
                    do while(abs(Vals(i)-Vals(j)).lt.1.0e-8)
                        j = j + 1
                        if(j.gt.nLat) exit
                    enddo
                    if((j-i).gt.1) then
                        !Degeneracy
                        if((j-i).gt.2) then
                            !More than two fold degenerate. How do we rotate these??
                            call stop_all(t_r,'Cannot handle more than 2x degeneracy in the k-space hamiltonian')
                        endif
                        !Degeneracy is between i and i + 1 (or i and j - 1)
                        if((i+1).ne.(j-1)) call stop_all(t_r,'Indexing error')
                        !Rotate these vectors together in positive and negative linear combinations
                        Vec_temp(:) = (r_vecs(:,i) + r_vecs(:,i+1)) / sqrt(2.0_dp)
                        Vec_temp2(:) = (r_vecs(:,i) - r_vecs(:,i+1)) / sqrt(2.0_dp)
                        r_vecs(:,i) = Vec_temp(:)
                        r_vecs(:,i+1) = Vec_temp2(:)
                    endif
                    i = j 
                enddo
                deallocate(Vec_temp,Vec_temp2)

                allocate(r_vecs_real(nLat,nLat))
                r_vecs_real(:,:) = zero
                !Now, find the appropriate phase, such that the rotation will make the r_vecs real.
                !Apply the inverse of this rotation to the k_vecs, such that we end up with a complex set of
                !k-vectors (ordered by k-point), and real set of r_vecs (Ordered by energy).
                do i = 1,nLat   !Run through eigenvectors
                    phase = zero
                    if(tWriteOut) write(6,*) "Rotating eigenvector : ",i
                    do j = 1,nLat
                        if((abs(aimag(r_vecs(j,i))).gt.1.0e-9_dp).and.(abs(r_vecs(j,i)).gt.1.0e-7_dp)) then
                            !Find the phase factor for this eigenvector
                            phase = atan(aimag(r_vecs(j,i))/real(r_vecs(j,i)))
                            if(tWriteOut) write(6,*) "Eigenvector: ",i,j,phase,r_vecs(j,i) * exp(dcmplx(0.0_dp,-phase))
                            exit
                        endif
                    enddo
                    !The phase should be the same for all components of the eigenvector
                    r_vecs(:,i) = r_vecs(:,i) * exp(dcmplx(zero,-phase))
                    do j = 1,nLat
                        if(abs(aimag(r_vecs(j,i))).gt.1.0e-6) then
                            write(6,*) "Error rotating component: ",j
                            write(6,*) phase,r_vecs(j,i)
                            call stop_all(t_r,'Eigenvectors not rotated correctly - degeneracies?')
                        endif
                        r_vecs_real(j,i) = real(r_vecs(j,i),dp)
                    enddo
                enddo

                !Check again that these rotated r_vecs are correct eigenfunctions...
!                if(tCheck) then
!                    !Do these satisfy the original eigenvalue problem?
!                    allocate(Work(nLat))
!                    do i = 1,nLat
!                        call DGEMV('N',nLat,nLat,one,Ham,nLat,r_vecs_real(:,i),1,zero,Work,1)
!                        do j = 1,nLat
!                            if(abs(Work(j)-(Vals(i)*r_vecs_real(j,i))).gt.1.0e-8_dp) then
!                                call stop_all(t_r,'Eigensystem not computed correctly in real real basis')
!                            endif
!                        enddo
!                    enddo
!                    deallocate(Work)
!                    !write(6,*) "Eigensystem correctly computed and transformed to real real space"
!                endif

                Ham(:,:) = cmplx(r_vecs_real(:,:),zero,dp)
                deallocate(CompHam,r_vecs_real,r_vecs)
            else
                !No need to maintain real eigenvectors
                
                !Check that these r_vecs are correct eigenfunctions...
                if(tCheck) then
                    !Do these satisfy the original eigenvalue problem?
                    allocate(vec_temp(nLat))
                    do i = 1,nLat
                        call ZGEMV('N',nLat,nLat,zone,Ham,nLat,r_vecs(:,i),1,zzero,vec_temp,1)
                        do j = 1,nLat
                            if(abs(vec_temp(j)-(Vals(i)*r_vecs(j,i))).gt.1.0e-8_dp) then
                                call stop_all(t_r,'Eigensystem not computed correctly in real basis')
                            endif
                        enddo
                    enddo
                    deallocate(vec_temp)
                    !write(6,*) "Eigensystem computed correctly"
                endif

                Ham(:,:) = r_vecs(:,:)
            endif
        else
            if(tRealVectors_) then
                call stop_all(t_r,'Cannot provide complex wavefunction without kpoint symmetry and hope for real vectors!')
            endif
            !Normal real space diagonalization
            Vals(:) = 0.0_dp
            allocate(cWork(1))
            allocate(Work(max(1,3*nLat-2)))
            lWork = -1
            info = 0
            call zheev('V','L',nLat,Ham,nLat,Vals,cWork,lWork,Work,info)
            if(info.ne.0) call stop_all(t_r,'Workspace queiry failed')
            lwork=int(abs(work(1)))+1
            deallocate(cwork)
            allocate(cwork(lwork))
            call zheev('V','L',nLat,Ham,nLat,Vals,cWork,lWork,Work,info)
            if(info.ne.0) call stop_all(t_r,'Diag failed')
            deallocate(work,cwork)
!            call writevector(Vals,'Eigenvalues')
        endif

    end subroutine DiagOneEOp_z

    !Diagonalize a lattice matrix, with optional k-independent additional terms,
    !in kspace or real space.
    subroutine DiagNonHermLatticeHam(LatVals,LatVecs_L,LatVecs_R,ham,cham,k_ind_mat)
        use sort_mod_c_a_c_a_c, only: Order_zgeev_vecs 
        use sort_mod, only: Orthonorm_zgeev_vecs
        implicit none
        complex(dp), intent(out) :: LatVals(nSites)
        complex(dp), intent(out) :: LatVecs_R(nSites,nSites),LatVecs_L(nSites,nSites)
        real(dp), intent(in), optional :: ham(nSites,nSites)
        complex(dp), intent(in), optional :: cham(nSites,nSites)
        complex(dp), intent(in), optional :: k_ind_mat(nImp,nImp)   !An optinoal matrix to stripe through the space

        complex(dp), allocatable :: k_Ham(:,:),RVec(:,:),LVec(:,:),W_Vals(:)
        real(dp), allocatable :: Work(:)
        complex(dp), allocatable :: ztemp(:,:),ham_temp(:,:),cWork(:)
        integer :: i,j,kPnt,ind_1,ind_2,lwork,info
        character(len=*), parameter :: t_r='DiagLatticeHam'

        if((.not.present(ham)).and.(.not.present(cham))) then
            call stop_all(t_r,'Neither real or complex hamiltonian passed in')
        endif

        if(tDiag_kspace) then
            !Diagonalize in the k-space supercell
            !Self-energy is then added to each kpoint
            allocate(k_Ham(nImp,nImp))
            allocate(RVec(nImp,nImp))
            allocate(LVec(nImp,nImp))
            allocate(W_Vals(nImp))
            allocate(Work(max(1,2*nImp)))
            allocate(ztemp(nSites,nImp))

            allocate(ham_temp(nSites,nSites))
            if(present(ham)) then
                do i = 1,nSites
                    do j = 1,nSites
                        ham_temp(j,i) = cmplx(ham(j,i),zero,dp)
                    enddo
                enddo
            else
                ham_temp(:,:) = cham(:,:)
            endif

            do kPnt = 1,nKPnts
                ind_1 = ((kPnt-1)*nImp) + 1
                ind_2 = nImp*kPnt

                !Transform one-electron hamiltonian into this k-basis
                call ZGEMM('N','N',nSites,nImp,nSites,zone,Ham_temp,nSites,RtoK_Rot(:,ind_1:ind_2),nSites,zzero,ztemp,nSites)
                call ZGEMM('C','N',nImp,nImp,nSites,zone,RtoK_Rot(:,ind_1:ind_2),nSites,ztemp,nSites,zzero,k_Ham,nImp)

                if(present(k_ind_mat)) then
                    !Include a k-independent, complex matrix to all kpoints
                    k_Ham(:,:) = k_Ham(:,:) + k_ind_mat(:,:)
                endif

                !Now diagonalize this 
                allocate(cWork(1))
                lwork = -1
                info = 0
                call zgeev('V','V',nImp,k_Ham,nImp,W_Vals,LVec,nImp,RVec,nImp,cWork,lWork,Work,info)
                if(info.ne.0) call stop_all(t_r,'Workspace query failed')
                lwork = int(abs(cWork(1)))+1
                deallocate(cWork)
                allocate(cWork(lWork))
                call zgeev('V','V',nImp,k_Ham,nImp,W_Vals,LVec,nImp,RVec,nImp,cWork,lWork,Work,info)
                if(info.ne.0) call stop_all(t_r,'Diagonalization of 1-electron GF failed')
                deallocate(cWork)

                !Rotate eigenvalues and eigenvectors back into r-space
                call ZGEMM('N','N',nSites,nImp,nImp,zone,RtoK_Rot(:,ind_1:ind_2),nSites,LVec,nImp,  &
                    zzero,LatVecs_L(:,ind_1:ind_2),nSites)
                call ZGEMM('N','N',nSites,nImp,nImp,zone,RtoK_Rot(:,ind_1:ind_2),nSites,RVec,nImp,  &
                    zzero,LatVecs_R(:,ind_1:ind_2),nSites)
                LatVals(ind_1:ind_2) = W_Vals(:)

            enddo
            deallocate(Ham_temp,ztemp,Work,W_Vals,LVec,RVec,k_Ham)
        else
            !Diagonalize in r-space
            allocate(Ham_temp(nSites,nSites))
            if(present(ham)) then
                do i = 1,nSites
                    do j = 1,nSites
                        Ham_temp(j,i) = cmplx(ham(j,i),zero,dp)
                    enddo
                enddo
            else
                Ham_temp(:,:) = cham(:,:)
            endif
            if(present(k_ind_mat)) then
                !Stripe the complex self-energy through the AO one-electron hamiltonian
                call add_localpot_comp_inplace(Ham_temp,k_ind_mat)
            endif

            !Now, diagonalize the resultant non-hermitian one-electron hamiltonian
            LatVals = zzero
            LatVecs_L = zzero
            LatVecs_R = zzero
            allocate(Work(max(1,2*nSites)))
            allocate(cWork(1))
            lwork = -1
            info = 0
            call zgeev('V','V',nSites,Ham_temp,nSites,LatVals,LatVecs_L,nSites,LatVecs_R,nSites,cWork,lWork,Work,info)
            if(info.ne.0) call stop_all(t_r,'Workspace query failed')
            lwork = int(abs(cWork(1)))+1
            deallocate(cWork)
            allocate(cWork(lWork))
            call zgeev('V','V',nSites,Ham_temp,nSites,LatVals,LatVecs_L,nSites,LatVecs_R,nSites,cWork,lWork,Work,info)
            if(info.ne.0) call stop_all(t_r,'Diag of H + SE failed')
            deallocate(work,cWork,Ham_temp)

        endif
    
        !zgeev does not order the eigenvalues in increasing magnitude for some reason. Ass.
        !This will order the eigenvectors according to increasing *REAL* part of the eigenvalues
        call Order_zgeev_vecs(LatVals,LatVecs_L,LatVecs_R)
        !call writevectorcomp(W_Vals,'Eigenvalues ordered')
        !Now, bi-orthogonalize sets of vectors in degenerate sets, and normalize all L and R eigenvectors against each other.
        call Orthonorm_zgeev_vecs(nSites,LatVals,LatVecs_L,LatVecs_R)

    end subroutine DiagNonHermLatticeHam


    !Convert a real-space symmetric operator into a k-space operator.
    !In: The operator in real space. Size = nSuperCell x nSupercell
    !   However, it assumes periodicity, such that only the first unit cell, and its connections to the other unit cells
    !   is referenced. I.e. R_Ham(1:nUnitCell,:). In the future, this should be changed so that only this is called.
    !Out: The operator in k-space. In this, the operator is block-diagonal, with each block being of size nUnitCell x nUnitCell.
    subroutine Convert1DtoKSpace(R_Ham,nSuperCell,nUnitCell,k_Ham)
        use sort_mod, only: sort_real
        implicit none
        integer, intent(in) :: nSuperCell,nUnitCell
        real(dp), intent(in) :: R_Ham(nSuperCell,nSuperCell)
        complex(dp), intent(out) :: k_Ham(nSuperCell,nSuperCell)
        complex(dp) :: KPntHam(nUnitCell,nUnitCell)
        integer :: k,i,j,a,b,nKpnts,cell_start,cell_end,k_start,k_end
        complex(dp) :: phase
        real(dp) :: KPnt_val
        real(dp), allocatable :: K_Vals(:),Orbs(:,:),W(:),Work(:),Vals(:),KPnts(:)
        complex(dp), allocatable :: cWork(:)
        integer :: lWork,info
        character(len=*), parameter :: t_r='Convert1DtoKSpace'

        if(tWriteOut) then
            write(6,*) "Converting real space operator to k-space: "
        endif
        if(tAntiPeriodic) call stop_all(t_r,'Cannot convert to k-space with anti-periodic boundary conditions')
        if(.not.tPeriodic) call stop_all(t_r,'Need periodic boundary conditions to convert to k-space')

        k_Ham(:,:) = dcmplx(0.0_dp,0.0_dp)

        !Number of kpoints = nSupercell/nUnitCell
        if(mod(nSupercell,nUnitCell).ne.0) call stop_all(t_r,'Not integer number of unit cells in supercell')
        nKpnts = nSupercell / nUnitCell
        if(mod(nKpnts,2).eq.1) then
            call stop_all(t_r,'For some reason, I am not getting hermitian operators with odd numbers of kpoints. " &
     &           //"Debug this routine.')
        endif

        allocate(KPnts(nKpnts))
        !Allocate values for k-vectors -> just 1D to start with here
!        write(6,*) "Number of k-points: ",nKpnts

!        if(mod(nKpnts,2).eq.0) then
            !Number of k-points even. Just use equally spaced mesh starting at -pi/a, and working our way across
            do k = 1,nKPnts
                KPnts(k) = -pi/real(nUnitCell,dp) + (k-1)*(2.0_dp*pi/nKpnts)/real(nUnitCell,dp)
            enddo
!        else
!            !Ensure that there is a kpoint at the Gamma point, and then equally spaced (don't sample BZ boundary)
!            do k = 1,nKPnts
!                KPnts(k) = -pi/real(nUnitCell,dp) + k*(2.0_dp*pi/(nKpnts+1))/real(nUnitCell,dp)
!            enddo
!        endif
!        call writevector(KPnts,'KPoint values')

        do k = 1,nKpnts
            !Create each block
            KPntHam(:,:) = R_Ham(1:nUnitCell,1:nUnitCell)   !The operator of the unit cell
            !write(6,*) "KPnt ",k,0,R_Ham(1:nUnitCell,1:nUnitCell)
            KPnt_val = KPnts(k)
            !write(6,*) "KPnt_val: ",KPnt_val

            do i = 1,nKpnts-1   !Real space translation lattice vectors to the (i+1)th cell
                cell_start = i*nUnitCell + 1
                cell_end = (i+1)*nUnitCell
                phase = exp(dcmplx(0.0_dp,KPnt_val*real(nUnitCell*i,dp))) !Phase between cell 1 and cell i+1
                !Add to the current kpoint, the opertor over the translated sites x by the phase change to them
                KPntHam(:,:) = KPntHam(:,:) + phase*R_Ham(1:nUnitCell,cell_start:cell_end)
                !write(6,*) "KPnt ",k,i,phase,phase*R_Ham(1:nUnitCell,cell_start:cell_end)
            enddo
                
            !TEST: Check hermiticity at all points
            do a = 1,nUnitCell
                do b = a,nUnitCell
                    if(abs(KPntHam(b,a)-conjg(KPntHam(a,b))).gt.1.0e-9_dp) then
                        write(6,*) a,b,KPntHam(a,b),KPntHam(b,a)
                        call writematrix(R_Ham(1:nUnitCell,:),'Coupling in real space',.true.)
                        call writematrix(R_Ham(:,:),'R_Ham',.true.)
                        call writematrixcomp(KPntHam,'Hamiltonian at kpoint',.true.)
                        call stop_all(t_r,'k-space operator not hermitian. Does input operator have periodicity')
                    endif
                enddo
            enddo

            !Now add this k-point to the full k-space operator
            k_start = (k-1)*nUnitCell + 1
            k_end = k*nUnitCell

            k_Ham(k_start:k_end,k_start:k_end) = KPntHam(:,:)
        enddo

        if(tWriteOut) then
            write(6,*) "Diagonal part of k-space operator: "
            do i = 1,nSuperCell
                write(6,*) i,k_Ham(i,i)
            enddo
        endif

        if(tCheck) then
            !TEST: Diagonalize real-space hamiltonian and check that eigenvalues are the same
            !We should now have the eigenvalues of the hamiltonian
            !Diagonalize the real-space hamiltonian and check that we have got this
            allocate(Orbs(nSuperCell,nSuperCell))
            Orbs(:,:) = R_Ham(:,:)
            allocate(W(nSuperCell))
            W(:) = 0.0_dp
            allocate(Work(1))
            lWork=-1
            info=0
            call dsyev('V','L',nSuperCell,Orbs,nSuperCell,W,Work,lWork,info)
            if(info.ne.0) call stop_all(t_r,'Workspace queiry failed')
            lwork=int(work(1))+1
            deallocate(work)
            allocate(work(lwork))
            call dsyev('V','L',nSuperCell,Orbs,nSuperCell,W,Work,lWork,info)
            if(info.ne.0) call stop_all(t_r,'Diag failed')
            deallocate(work)

            allocate(K_Vals(nSuperCell))
            if(nUnitCell.eq.1) then
                do i = 1,nSuperCell
                    K_Vals(i) = real(k_Ham(i,i),dp)
                enddo
            else
                lWork = max(1,2*nUnitCell-1)
                allocate(cWork(lWork))
                allocate(Work(max(1,3*nUnitCell-2)))
                allocate(Vals(nUnitCell))
                do k = 1,nKPnts
                    !Diagonalize block
                    k_start = (k-1)*nUnitCell + 1
                    k_end = k*nUnitCell
                    KPntHam(:,:) = k_Ham(k_start:k_end,k_start:k_end)
                    !Hermitian matrix diagonalization
                    call writematrixcomp(KPntHam,'KPntHam',.true.)
                    call ZHEEV('N','U',nUnitCell,KPntHam,nUnitCell,Vals,cWork,lWork,Work,info)
                    if(info.ne.0) call stop_all(t_r,'Diag Failed')
                    do j = 1,nUnitCell
                        K_Vals(k_start+j-1) = Vals(j)
                    enddo
                enddo
                deallocate(cWork,Work,Vals)
            endif

            !Sort values
            call sort_real(K_Vals,nSuperCell)

            do i = 1,nSuperCell
                !write(6,*) i,K_Vals(i),W(i)
                if(abs(K_Vals(i)-W(i)).gt.1.0e-9_dp) then
                    write(6,*) i,K_Vals(i),W(i)
                    call stop_all(t_r,'Conversion to k-space failed')
                endif
            enddo
            deallocate(K_Vals,W,Orbs)
        endif

        deallocate(KPnts)

    end subroutine Convert1DtoKSpace

    !The error metric used for the fitting of the self-energy in order to match the greens functions
    !IN: SE is the guess for the self-energy *correction* over the impurity sites (packed form)
    !    HL_GF is the set of greens functions for the DMET calculation
    !OUT: GF_Diff is the difference between the High-level greens functions, and the NI GF with the self-energy correction added (In packed form)
    subroutine GFErr(se,GF_Diff,HL_GF,Omega)
        implicit none
        complex(dp), intent(in) :: se(nVarSE)    !The guess for the self-energy *correction* (packed)
        complex(dp), intent(in) :: HL_GF(nImp,nImp) !The DMET calculated greens functions over all impurity sites
        real(dp), intent(in) :: Omega
        complex(dp), intent(out) :: GF_Diff(nVarSE) 
        complex(dp) :: ni_GFs(nImp,nImp),GF_Diff_unpacked(nImp,nImp)

        !Add se to h0v_SE, diagonalize and construct the non-interacting solutions for all impurity sites
        call mkgf(se,ni_GFs,Omega)
        GF_Diff_unpacked(:,:) = ni_GFs(:,:) - HL_GF(:,:)
        call ToCompPacked(nImp,GF_Diff_unpacked,GF_Diff)

    end subroutine GFErr

    !Add se to h0v_SE, diagonalize and construct the non-interacting greens functions between all impurity sites (unpacked)
    !This should really be diagonalized in k-space
    subroutine mkgf(se,ni_GFs,Omega)
        use sort_mod_c_a_c_a_c, only: Order_zgeev_vecs 
        use sort_mod, only: Orthonorm_zgeev_vecs
        implicit none
        complex(dp), intent(in) :: se(nVarSE)
        real(dp), intent(in) :: Omega
        complex(dp), intent(out) :: ni_GFs(nImp,nImp)
        complex(dp) :: se_unpacked(nImp,nImp)
        integer :: lWork,info,i,pertsite,pertBra
        complex(dp), allocatable :: AO_Ham(:,:),W_Vals(:),RVec(:,:),LVec(:,:),cWork(:)
        !complex(dp) :: NI_Cre(nImp,nImp),NI_Ann(nImp,nImp),NI_GF_Check(nImp,nImp)
        !complex(dp), allocatable :: HF_Ann_Ket(:,:),HF_Cre_Ket(:,:)
        real(dp), allocatable :: Work(:)
        character(len=*), parameter :: t_r='mkgf'

        call FromCompPacked(nImp,se,se_unpacked)

        !Now, stripe the (-)new self energy through the space
        allocate(AO_Ham(nSites,nSites))
        AO_Ham(:,:) = zzero
        call add_localpot_comp(h0v_SE,AO_Ham,se_unpacked,tAdd=.false.)

        !Now, diagonalize the resultant non-hermitian one-electron hamiltonian
        allocate(W_Vals(nSites))
        allocate(RVec(nSites,nSites))
        allocate(LVec(nSites,nSites))
        RVec = zzero
        LVec = zzero
        W_Vals = zzero
        allocate(Work(max(1,2*nSites)))
        allocate(cWork(1))
        lwork = -1
        info = 0
        call zgeev('V','V',nSites,AO_Ham,nSites,W_Vals,LVec,nSites,RVec,nSites,cWork,lWork,Work,info)
        if(info.ne.0) call stop_all(t_r,'Workspace query failed')
        lwork = int(abs(cWork(1)))+1
        deallocate(cWork)
        allocate(cWork(lWork))
        call zgeev('V','V',nSites,AO_Ham,nSites,W_Vals,LVec,nSites,RVec,nSites,cWork,lWork,Work,info)
        if(info.ne.0) call stop_all(t_r,'Diag of H - SE failed')
        deallocate(work,cWork,AO_Ham)

        !zgeev does not order the eigenvalues in increasing magnitude for some reason. Ass.
        !This will order the eigenvectors according to increasing *REAL* part of the eigenvalues
        call Order_zgeev_vecs(W_Vals,LVec,RVec)
        !call writevectorcomp(W_Vals,'Eigenvalues ordered')
        !Now, bi-orthogonalize sets of vectors in degenerate sets, and normalize all L and R eigenvectors against each other.
        call Orthonorm_zgeev_vecs(nSites,W_Vals,LVec,RVec)

        ni_gfs(:,:) = zzero
        do pertsite = 1,nImp
            do pertBra = 1,nImp
                do i = 1,nSites
                    ni_GFs(pertsite,pertBra) = ni_GFs(pertsite,pertBra) +   &
                        RVec(pertBra,i)*dconjg(LVec(pertsite,i))/(dcmplx(Omega,dDelta)-W_Vals(i))
                enddo
            enddo
        enddo

!        allocate(HF_Ann_Ket(nOcc,nImp))
!        allocate(HF_Cre_Ket(nOcc+1:nSites,nImp))
!        HF_Ann_Ket(:,:) = zzero
!        HF_Cre_Ket(:,:) = zzero
!        NI_Ann(:,:) = zzero
!        NI_Cre(:,:) = zzero
!        do pertsite = 1,nImp
!            do i = 1,nOcc
!                HF_Ann_Ket(i,pertsite) = dconjg(LVec(pertsite,i))/(dcmplx(Omega,dDelta)-W_Vals(i))
!            enddo
!            do a = nOcc+1,nSites
!                HF_Cre_Ket(a,pertsite) = dconjg(LVec(pertsite,a))/(dcmplx(Omega,dDelta)-W_Vals(a))
!            enddo
!            do pertBra = 1,nImp
!                do i = 1,nOcc
!                    NI_Ann(pertsite,pertBra) = NI_Ann(pertsite,pertBra) + RVec(pertBra,i)*HF_Ann_Ket(i,pertsite)
!                    !write(6,*) "mkgf: ",i,NI_Ann(pertsite,pertBra),RVec(pertBra,i),HF_Ann_Ket(i,pertsite)
!                enddo
!                do a = nOcc+1,nSites
!                    NI_Cre(pertsite,pertBra) = NI_Cre(pertsite,pertBra) + RVec(pertBra,a)*HF_Cre_Ket(a,pertsite)
!                enddo
!            enddo
!        enddo
!        NI_GF_Check(:,:) = NI_Cre(:,:) + NI_Ann(:,:)
!
!        do pertsite = 1,nImp
!            do pertBra = 1,nImp
!                if(abs(NI_GF_Check(pertBra,pertsite)-ni_GFs(pertBra,pertsite)).gt.1.0e-8_dp) then
!                    call stop_all(t_r,'NI GFs not consistent')
!                endif
!            enddo
!        enddo
!        deallocate(HF_Ann_Ket,HF_Cre_Ket)
!
!        call writematrixcomp(se_unpacked,'Delta SE',.true.)
!        write(6,*) "Full NI GF: ",NI_GF_Check(:,:)
!        write(6,*) "mkgf NI GF: ",ni_GFs(:,:)

        deallocate(W_Vals,RVec,LVec)
    end subroutine mkgf
    
    !The error metric used for the fitting of the vloc in order to match the RDMs
    !The error metric is not actually calculated, but can be considered as the squared sum of the elements in the
    !returned matrix. The matrix is then the gradients in each direction, and the jacobian is made numerically.
    !IN: vloc over impurity sites in triangular packed form. This is the *correction* to the correlation potential
    !OUT: Error matrix between the systems (just difference over all embedded sys) (triangular packed)
    subroutine RDMErr(v,ErrMat_packed)
        implicit none
        real(dp), intent(in) :: v(nImpCombs)
        real(dp), intent(out) :: ErrMat_packed(EmbCombs)
        real(dp) :: ErrMat_unpacked(EmbSize,EmbSize)
        real(dp) :: MFRdm(EmbSize,EmbSize)

        call mkrdm(v,MFRdm)
        ErrMat_unpacked(:,:) = MFRdm(:,:) - HL_1RDM(:,:)
!        call writematrix(MFRdm,'MFRdm',.true.)
        call ToTriangularPacked(EmbSize,ErrMat_unpacked,ErrMat_packed) 
    end subroutine RDMErr

    !Construct rdm over the embedded system from the fock + potential v (over impurity sites)
    !This will add the v to the impurity sites on the fock matrix, before diagonalizing, and rotating back to the RDM in the embedded system
    !with the original meanfield+vloc occupation numbers
    subroutine mkrdm(v,rdm)
        implicit none
        real(dp), intent(in) :: v(nImpCombs)
        real(dp), intent(out) :: rdm(EmbSize,EmbSize)
        real(dp) :: EValues(EmbSize),EVectors(EmbSize,EmbSize)
        real(dp) :: temp(EmbSize,EmbSize),temp2(EmbSize,EmbSize)
        integer :: i

        !Take the diagonalized system from mkorb, and construct the RDM in the basis of meanfield solution,
        call mkorb(v,EValues,EVectors)

!        call writevector(Evalues,'Evalues')
!        call writematrix(EVectors,'EVectors',.true.)

        !Now transform the occupation numbers from the embedded system into the mean-field+new_vloc natural orbital embedded basis
        ! U,diag(),U^T. This is what we want to match.
        !Create diagonal matrix
        temp(:,:) = 0.0_dp
        do i=1,EmbSize
            temp(i,i) = MFEmbOccs(i)
        enddo
!        call writematrix(temp,'MFEmbOccs',.true.)
        call DGEMM('N','N',EmbSize,EmbSize,EmbSize,1.0_dp,EVectors,EmbSize,temp,EmbSize,0.0_dp,temp2,EmbSize)
        call DGEMM('N','T',EmbSize,EmbSize,EmbSize,1.0_dp,temp2,EmbSize,EVectors,EmbSize,0.0_dp,rdm,EmbSize)

    end subroutine mkrdm

    !Take a potential over the impurity system (in triangular form), and add it to the fock matrix in the embedded system (only over the impurity sites) 
    !Diagonalize this embedded system and return the eigenvalues and vectors
    subroutine mkorb(v,EValues,EVectors)
        implicit none
        real(dp), intent(in) :: v(nImpCombs)
        real(dp), intent(out) :: EValues(EmbSize),EVectors(EmbSize,EmbSize)
        real(dp) :: vloc_unpacked(nImp,nImp)
        real(dp), allocatable :: work(:)
        integer :: lWork,info
        character(len=*), parameter :: t_r='mkorb'

        !Unpack potential
        call FromTriangularPacked(nImp,v,vloc_unpacked)
        !Add the potential over the impurity sites to the fock matrix over the entire embedded system
        EVectors(:,:) = Emb_Fock(:,:)
        EVectors(1:nImp,1:nImp) = EVectors(1:nImp,1:nImp) + vloc_unpacked(:,:)

        !EVectors now contains Fock + the vlocal over the impurity sites
        !Now diagonalize
        allocate(Work(1))
        lWork=-1
        info=0
        call dsyev('V','U',EmbSize,EVectors,EmbSize,EValues,Work,lWork,info)
        if(info.ne.0) call stop_all(t_r,'Workspace queiry failed')
        lwork=int(work(1))+1
        deallocate(work)
        allocate(work(lwork))
        call dsyev('V','U',EmbSize,EVectors,EmbSize,EValues,Work,lWork,info)
        if(info.ne.0) call stop_all(t_r,'Diag failed')
        deallocate(work)

    end subroutine mkorb

    !Pack a 2D complex matrix to a 1D array.
    !To start, assume nothing about the array
    subroutine ToCompPacked(n,unpacked,packed)
        implicit none
        integer, intent(in) :: n
        complex(dp), intent(in) :: unpacked(n,n)
        complex(dp), intent(out) :: packed(nVarSE)
        integer :: i,j,k
        character(len=*), parameter :: t_r='ToCompPacked'

        packed(:) = zzero
        if(iSE_Constraints.eq.1) then
            k = 1
            do i = 1,n
                do j = 1,n
                    packed(k) = unpacked(j,i)
                    k = k + 1
                enddo
            enddo
        else
            k = 1
            do i = 1,n
                do j = 1,i
                    packed(k) = unpacked(j,i)
                    if((i.ne.j).and.(abs(unpacked(j,i)-dconjg(unpacked(i,j))).gt.1.0e-8_dp)) then
                        write(6,*) "i,j ", i,j, unpacked(j,i),unpacked(i,j)
                        call stop_all(t_r,'Off-diagonal hermiticity in self-energy lost.')
                    endif
                    k = k + 1
                enddo
            enddo
        endif

    end subroutine ToCompPacked

    subroutine FromCompPacked(n,packed,unpacked)
        implicit none
        integer, intent(in) :: n
        complex(dp), intent(in) :: packed(nVarSE)
        complex(dp), intent(out) :: unpacked(n,n)
        integer :: i,j,k

        unpacked(:,:) = zzero
        if(iSE_Constraints.eq.1) then
            k = 1
            do i = 1,n
                do j = 1,n
                    unpacked(j,i) = packed(k)
                    k = k + 1
                enddo
            enddo
        else
            !Store upper triangle of packed hermitian matrix
            k = 1
            do i = 1,n
                do j = 1,i
                    unpacked(j,i) = packed(k)
                    if(i.ne.j) then
                        unpacked(i,j) = dconjg(packed(k))
                    endif
                    k = k + 1
                enddo
            enddo
        endif

    end subroutine FromCompPacked

    !Routine to triangular pack a matrix and return it in 'Packed'.
    pure subroutine ToTriangularPacked(Length,Unpacked,Packed)
        implicit none
        integer, intent(in) :: Length
        real(dp) , intent(in) :: Unpacked(Length,Length)
        real(dp) , intent(out) :: Packed((Length*(Length+1))/2)
        integer :: i,j,k

        Packed(:) = 0.0_dp
        k=1
        do i=1,Length
            do j=1,i
                Packed(k) = Unpacked(i,j)
                k=k+1
            enddo
        enddo

    end subroutine ToTriangularPacked

    pure subroutine FromTriangularPacked(Length,Packed,Unpacked)
        implicit none
        integer, intent(in) :: Length
        real(dp) , intent(out) :: Unpacked(Length,Length)
        real(dp) , intent(in) :: Packed((Length*(Length+1))/2)
        integer :: i,j,k
        k=1
        do i=1,Length
            do j=1,i
                Unpacked(i,j) = Packed(k)
                Unpacked(j,i) = Packed(k)
                k=k+1
            enddo
        enddo
    end subroutine FromTriangularPacked

    subroutine WriteMatrixcomp(mat,matname,tOneLine)
        implicit none
        complex(dp), intent(in) :: mat(:,:)
        character(len=*), intent(in) :: matname
        integer :: i,j
        logical :: tOneLine

        write(6,*) "Writing out matrix: ",trim(matname)
        write(6,"(A,I7,A,I7)") "Size: ",size(mat,1)," by ",size(mat,2)
        do i=1,size(mat,1)
            do j=1,size(mat,2)
                if(tOneLine) then
                    write(6,"(2G18.7)",advance='no') mat(i,j)
                else
                    write(6,"(2I6,2G18.7)") i,j,mat(i,j)
                endif
            enddo
            write(6,*)
        enddo
    end subroutine WriteMatrixcomp

    subroutine WriteMatrix(mat,matname,tOneLine)
        implicit none
        real(dp), intent(in) :: mat(:,:)
        character(len=*), intent(in) :: matname
        integer :: i,j
        logical :: tOneLine

        write(6,*) "Writing out matrix: ",trim(matname)
        write(6,"(A,I7,A,I7)") "Size: ",size(mat,1)," by ",size(mat,2)
        do i=1,size(mat,1)
            do j=1,size(mat,2)
                if(tOneLine) then
                    write(6,"(G25.10)",advance='no') mat(i,j)
                else
                    write(6,"(2I6,G25.10)") i,j,mat(i,j)
                endif
            enddo
            write(6,*)
        enddo
    end subroutine WriteMatrix

    subroutine WriteVector(vec,vecname)
        implicit none
        real(dp), intent(in) :: vec(:)
        character(len=*), intent(in) :: vecname
        integer :: i

        write(6,*) "Writing out vector: ",trim(vecname)
        write(6,"(A,I7,A,I7)") "Size: ",size(vec,1)
        do i=1,size(vec,1)
!            write(6,"(G25.10)",advance='no') vec(i)
            write(6,"(G25.10)") vec(i)
        enddo
        write(6,*)
    end subroutine WriteVector
    
    subroutine WriteVectorcomp(vec,vecname)
        implicit none
        complex(dp), intent(in) :: vec(:)
        character(len=*), intent(in) :: vecname
        integer :: i

        write(6,*) "Writing out vector: ",trim(vecname)
        write(6,"(A,I7,A,I7)") "Size: ",size(vec,1)
        do i=1,size(vec,1)
!            write(6,"(G25.10)",advance='no') vec(i)
            write(6,"(2G25.10)") vec(i)
        enddo
        write(6,*)
    end subroutine WriteVectorcomp
    
    subroutine WriteVectorInt(vec,vecname)
        implicit none
        integer, intent(in) :: vec(:)
        character(len=*), intent(in) :: vecname
        integer :: i

        write(6,*) "Writing out vector: ",trim(vecname)
        write(6,"(A,I7,A,I7)") "Size: ",size(vec,1)
        do i=1,size(vec,1)
!            write(6,"(G25.10)",advance='no') vec(i)
            write(6,"(I12)") vec(i)
        enddo
        write(6,*)
    end subroutine WriteVectorInt

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !*****************************************************************************
    function znrm2 ( n, x, incx )
    !*****************************************************************************
    !
    !! SCNRM2 returns the euclidean norm of a complex(kind=dp) vector.
    !
    !
    !  Discussion:
    !
    !    SCNRM2 := sqrt ( sum ( conjg ( x(1:n) ) * x(1:n) ) )
    !            = sqrt ( dot_product ( x(1:n), x(1:n) ) )
    !
    !  Parameters:
    !
    !    Input, integer N, the number of entries in the vector.
    !
    !    Input, complex(kind=dp) X(*), the vector.
    !
    !    Input, integer INCX, the increment between successive entries of X.
    !
    !    Output, real(kind=dp) SCNRM2, the norm of the vector.
    !
      implicit none
    !
      integer(ip), intent(in)      :: incx
      integer(ip)                  :: ix
      integer(ip), intent(in)      :: n
      real(kind=dp)                :: norm
     !real(kind=dp), parameter     :: one = 1.0_dp
      real(kind=dp)                :: scale
      real(kind=dp)                :: znrm2
      real(kind=dp)                :: ssq
      real(kind=dp)                :: temp
      complex(kind=dp), intent(in) :: x(*)
     
    !
      if ( n < 1 .or. incx < 1 ) then

        norm  = zero

      else

        scale = zero
        ssq = one

        do ix = 1, 1 + ( n - 1 ) * incx, incx
          if ( real(x(ix), dp) /= zero ) then  
            temp = abs ( real(x(ix), dp) )   
            if ( scale < temp ) then
              ssq = one + ssq * ( scale / temp )**2
              scale = temp
            else
              ssq = ssq + ( temp / scale )**2
            end if
          end if

          if ( aimag ( x(ix) ) /= zero ) then
            temp = abs ( aimag ( x(ix) ) )
            if ( scale < temp ) then
              ssq = one + ssq * ( scale / temp )**2
              scale = temp
            else
              ssq = ssq + ( temp / scale )**2
            end if

          end if

        end do

        norm  = scale * sqrt ( ssq )

      end if

      znrm2 = norm

      return
    end function znrm2

end module mat_tools
