!tools for matrix construction, manipulation and HF solver for DMET
module mat_tools
    use const
    use errors, only: stop_all, warning
    use globals
    implicit none

    contains

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
        integer :: i,j,k,x,y,li,lj,ij_link,site,dx,dy,iunit
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

        elseif(LatticeDim.eq.2) then

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
!            write(6,*) "Hopping matrix in lattice ordering: "
!            do i = 1,nSites
!                do j = 1,nSites
!                    write(6,"(F7.3)",advance='no') h0(j,i)
!                enddo
!                write(6,*) ""
!            enddo
            call Mat_to_imp_order(h0)
        else
            call stop_all(t_r,'Higher dimensional models not coded up')
        endif

!        if(tUHF) then
!            !TODO: THIS IS JUST FOR TESTING!!
!            h0_b(:,:) = h0(:,:)
!        endif

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
        integer :: i,j,k,a,b
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
        integer :: i,j,k,a,b
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
        integer :: i,j,k,a,b
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
        elseif(tChemPot) then
            fock(1,1) = fock(1,1) - (U/2.0_dp)
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
            if(tChemPot.and.tAnderson) h0(1,1) = h0(1,1) - U/2.0_dp
            call dgemm('t','n',nSites,nSites,nSites,1.0_dp,FullHFOrbs,nSites,h0,nSites,0.0_dp,temp,nSites)
            call dgemm('n','n',nSites,nSites,nSites,1.0_dp,temp,nSites,FullHFOrbs,nSites,0.0_dp,h0HF,nSites)
            if(tChemPot.and.tAnderson) h0(1,1) = h0(1,1) + U/2.0_dp
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
        HFEnergies(:) = zero  
        if(tDiag_kspace) then
            call DiagOneEOp(HFOrbs,HFEnergies,nImp,nSites,tDiag_kspace)
        else
            call DiagOneEOp(HFOrbs,HFEnergies,nImp,nSites,.false.)
        endif
        if(tUHF) then
            HFOrbs_b(:,:) = Fock_b(:,:)
            HFEnergies_b(:) = zero
            if(tDiag_kspace) then
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
        call DGEMM('N','T',nSites,nSites,nOcc,Occupancy,OccOrbs,nSites,OccOrbs,nSites,0.0_dp,MeanFieldDM,nSites)

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
    subroutine setup_kspace()
        implicit none
        integer :: SS_Period,i,k,ind_1,ind_2,j
        real(dp) :: PrimLattVec(LatticeDim),phase,ddot
        complex(dp) , allocatable :: temp(:,:)
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
                    tShift_Mesh = .false.
                else
                    !We want a Monkhort-Pack mesh
                    tShift_Mesh = .true.
                endif
            else
                if(tPeriodic) then
                    !We want to use a Gamma centered mesh
                    tShift_Mesh = .false.
                else
                    !Monkhorst-Pack mesh
                    tShift_Mesh = .true.
                endif
            endif
            if(tShift_Mesh) then
                write(6,"(A)") "Using a Monkhorst-pack kpoint mesh of 1st BZ - no gamma point"
            else
                write(6,"(A)") "Using a gamma-centered kpoint mesh of 1st BZ - gamma point included"
            endif

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
            call stop_all(t_r,'Cannot do k-space diagonalizations - impurity site tiling is not same as direct lattice')
            !Below is commented out code to ostensibly create 2D kpoint mesh
!            if(SS_Period.le.2) then
!                !1/2 impurity. Tilted lattice with real space lattice vector (1,1) and (1,-1)
!                !Always have 2 sites per direct unit cell
!                RecipLattVecs(:,1) = pi
!                RecipLattVecs(1,2) = pi
!                RecipLattVecs(2,2) = -pi
!
!                tShift = .false.
!
!                !Is it a regular grid?
!                kPerDim = nint(sqrt(real(nKPnts,dp)))
!                if(abs(real(kPerDim**2,dp)-real(nKPnts,dp)).gt.1.0e-8_dp) then
!                    write(6,*) "k-points per dimension: ",sqrt(real(nKPnts,dp))
!                    write(6,*) "Bands (needs to be square number for uniform mesh): ",nKPnts
!                    call stop_all(t_r,'Cannot do k-space diagonalization with non-uniform mesh (currently?)')
!                endif
!
!                do i = 1,kPerDim
!                    do j = 1,kPerDim
!                        kpnt = (i-1)*kPerDim + j
!
!                        KPnts(1,kpnt) = -pi    !Start at RHS of tilted FBZ
!                        KPnts(:,kpnt) = KPnts(:,kpnt) + (i-1)*RecipLattVecs(:,1)/kPerDim + (j-1)*RecipLattVecs(:,2)/kPerDim
!                        if(tShift) then
!                            KPnts(:,kpnt) = KPnts(:,kpnt) + RecipLattVecs(:,1)/(2.0_dp*real(kPerDim,dp)) + &
!                                RecipLattVecs(:,2)/(2.0_dp*real(kPerDim,dp))
!                        endif
!                        write(6,*) "KPnt ",kpnt,KPnts(:,kpnt)
!                    enddo
!                enddo
!!            else
!!                call stop_all(t_r,'Cannot deal with > 2 impurities with k-space diagonalization atm')
!            endif
        else
            !Quoi?
            call stop_all(t_r,'Error here')
        endif

        if(tWriteOut) then
            write(6,*) "Writing out kpoint mesh: "
            do k=1,nKPnts
                write(6,*) "KPnt ",k,KPnts(:,k)
            enddo
        endif
            
        !Setup rotation matrix from site basis to k-space
        !First index r, second k
        allocate(RtoK_Rot(nSites,nSites))
        RtoK_Rot(:,:) = zzero
        !Run though all kpoints
        do k = 1,nKPnts
            !Construct rotation
            ind_1 = ((k-1)*nImp) + 1
            ind_2 = nImp*k
            do i = 1,nSites
                if(LatticeDim.eq.1) then
                    PrimLattVec(1) = real(i-1,dp)   !The real-space translation to this site
                else
                    call stop_all(t_r,'Error')
                endif
                phase = ddot(LatticeDim,KPnts(:,k),1,PrimLattVec,1)
                RtoK_Rot(i,ind_1+mod(i,nImp)) = exp(dcmplx(zero,phase))/sqrt(real(nKPnts,dp))
            enddo
        enddo

        if(tWriteOut) call writematrixcomp(RtoK_Rot,'RtoK_Rot',.true.)

        if(tCheck) then
            !Now, is RtoK_Rot unitary?
            allocate(temp(nSites,nSites))
            !Check unitarity of matrix
            call ZGEMM('C','N',nSites,nSites,nSites,zone,RtoK_Rot,nSites,RtoK_Rot,nSites,zzero,temp,nSites) 
            do i = 1,nSites
                do j = 1,nSites
                    if((i.eq.j).and.(abs(temp(i,j)-zone).gt.1.0e-7_dp)) then
                        write(6,*) "i,j: ",i,j
                        call writematrixcomp(temp,'Identity?',.true.)
                        call stop_all(t_r,'Rotation matrix not unitary')
                    elseif((i.ne.j).and.(abs(temp(j,i)).gt.1.0e-7_dp)) then
                        write(6,*) "i,j: ",i,j
                        call writematrixcomp(temp,'Identity?',.true.)
                        call stop_all(t_r,'Rotation matrix not unitary 2')
                    endif
                enddo
            enddo
            !Try other way...
            call ZGEMM('N','C',nSites,nSites,nSites,zone,RtoK_Rot,nSites,RtoK_Rot,nSites,zzero,temp,nSites) 
            do i = 1,nSites
                do j = 1,nSites
                    if((i.eq.j).and.(abs(temp(i,j)-zone).gt.1.0e-7_dp)) then
                        write(6,*) "i,j: ",i,j
                        call writematrixcomp(temp,'Identity?',.true.)
                        call stop_all(t_r,'Rotation matrix not unitary')
                    elseif((i.ne.j).and.(abs(temp(j,i)).gt.1.0e-7_dp)) then
                        write(6,*) "i,j: ",i,j
                        call writematrixcomp(temp,'Identity?',.true.)
                        call stop_all(t_r,'Rotation matrix not unitary 2')
                    endif
                enddo
            enddo
            write(6,*) "Rotation matrix unitary... :)"
            deallocate(temp)
        endif

    end subroutine setup_kspace

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

    !Routine to diagonalize a 1-electron, real operator, with the lattice periodicity.
    !the SS_Period is the size of the supercell repeating unit (e.g. the coupling correlation potential)
    !This will be equivalent to the number of bands per kpoint
    !Ham returns the eigenvectors (in real space).
    subroutine DiagOneEOp(Ham,Vals,SS_Period,nLat,tKSpace_Diag)
        implicit none
        integer, intent(in) :: nLat
        logical, intent(in) :: tKSpace_Diag
        integer, intent(in) :: SS_Period
        real(dp), intent(inout) :: Ham(nLat,nLat)
        real(dp), intent(out) :: Vals(nLat)
        real(dp), allocatable :: Work(:),r_vecs_real(:,:)
        real(dp) :: PrimLattVec(LatticeDim),SiteVec(LatticeDim),Expo,DDOT,phase
        complex(dp), allocatable :: RotMat(:,:),temp(:,:),CompHam(:,:),RotHam(:,:),cWork(:)
        complex(dp), allocatable :: CompHam_2(:,:),ztemp(:,:),k_Ham(:,:),k_vecs(:,:)
        complex(dp), allocatable :: r_vecs(:,:)
        integer :: lWork,info,i,j,k,kSpace_ind,ki,kj,bandi,bandj,Ind_i,Ind_j,ind_1,ind_2
        integer :: ii,jj,xb,yb,impy,impx,impsite
        character(len=*), parameter :: t_r='DiagOneEOp'

        if(tKSpace_Diag) then
            if(LatticeDim.eq.2) then
                call stop_all(t_r,'Cannot do k-space diagonalizations - impurity site tiling is not same as direct lattice')
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
                    call ZGEMV('N',nLat,nLat,zone,CompHam,nLat,r_vecs(:,i),1,zzero,cWork,1)
                    do j = 1,nLat
                        if(abs(cWork(j)-(Vals(i)*r_vecs(j,i))).gt.1.0e-8_dp) then
                            call stop_all(t_r,'Eigensystem not computed correctly in real basis')
                        endif
                    enddo
                enddo
                deallocate(cWork)
                write(6,*) "Eigensystem correctly computed and transformed to real space"
            endif

            write(6,*) "Before sorting..."
            call writevector(Vals,'Vals')
            call writematrixcomp(r_vecs,'r_vecs',.true.)
            
            !Order the vectors, such that the are in order of increasing eigenvalue
            call sort_d_a_c(Vals,r_vecs,nSites,nSites)
            
            write(6,*) "After sorting..."
            call writevector(Vals,'Vals')
            call writematrixcomp(r_vecs,'r_vecs',.true.)
            
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
                write(6,*) "Eigensystem correctly computed and transformed to real space"
            endif

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
                        !temp(j,i) = temp(j,i) * exp(dcmplx(0.0_dp,-phase))
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

            !Degenerate sets?

            !Check again that these rotated r_vecs are correct eigenfunctions...
            if(tCheck) then
                !Do these satisfy the original eigenvalue problem?
                allocate(Work(nLat))
                do i = 1,nLat
                    call DGEMV('N',nLat,nLat,one,Ham,nLat,r_vecs_real(:,i),1,zero,Work,1)
                    do j = 1,nLat
                        if(abs(Work(j)-(Vals(i)*r_vecs_real(j,i))).gt.1.0e-8_dp) then
                            call stop_all(t_r,'Eigensystem not computed correctly in real real basis')
                        endif
                    enddo
                enddo
                deallocate(Work)
                write(6,*) "Eigensystem correctly computed and transformed to real real space"
            endif

            Ham(:,:) = r_vecs_real(:,:)
            deallocate(CompHam,r_vecs_real)
            return









                
!            allocate(RotMat(SS_Period,nLat))
!            do k = 1,nKPnts
!                !Construct the rotation matrix for each kpoint
!                !First index are the bands at that kpoint. Second index in the real-space basis
!                RotMat(:,:) = zzero
!
!                do i = 1,nLat
!                    if(LatticeDim.eq.1) then
!                        TransVec(1) = real(i-1,dp)  !Translation vector to this unit cell
!                    endif
!                    do j = 1,SS_Period  !Bands per kpoint
!
!                        Expo = ddot(LatticeDim,KPnts(:,k),1,TransVec(:),1)
!                        RotMat(j,i) = (1.0_dp/sqrt(real(nLat,dp))) * exp(dcmplx(0.0_dp,Expo))
!                    enddo
!                enddo
!
!
!            enddo

!            call writematrix(Ham,'Real space matrix',.true.)

            !Do incredibly naievly to start with by creating the entire rotation matrix and rotating
            !First index in k-space, second in r-space
            allocate(RotMat(nLat,nLat))
            RotMat(:,:) = zzero

            do i = 1,nLat
                if(LatticeDim.eq.1) then
                    !TransVec(1) = real(i-1,dp)  !/real(SS_Period,dp)
                    PrimLattVec(1) = real((i-1)/SS_Period,dp)*real(SS_Period,dp)  !All sites in the same lattice translation have the same primitive lattice vector
                    SiteVec(1) = real(mod(i-1,SS_Period),dp)    !Vector within the lattice to this basis function
                    impsite = mod(i,SS_Period)
                else
                    !which impurity site do they belong to?
                    call site2ij(i-1,ii,jj)
                    call ij2xy(ii,jj,xb,yb)
                    impx = py_mod(xb,nImp_x)
                    impy = py_mod(yb,nImp_y)
                    impsite = impy*nImp_y + impx
                    write(6,*) "Site located at impurity location: ",i,impx,impy,impsite
                endif
!                write(6,*) "Cell Translation: ",PrimLattVec(:)
!                write(6,*) "Site Vector: ",SiteVec(:)

                do k = 1,nKPnts

                    do j = 1,SS_Period  !Bands per kpoint

                        kSpace_ind = SS_Period*(k-1) + j

                        !Ensure that we only want to loop over bands corresponding to this lattice repeat, not all sites
                        if(mod(j,SS_Period).ne.impsite) cycle

                        Expo = ddot(LatticeDim,KPnts(:,k),1,PrimLattVec(:),1)
!                        phase = ddot(LatticeDim,KPnts(:,k),1,SiteVec(:),1)

                        RotMat(kSpace_ind,i) = (1.0_dp/sqrt(real(nLat/SS_Period,dp))) *     &
                            exp(dcmplx(0.0_dp,Expo)) !* exp(dcmplx(0.0_dp,phase))
                    enddo
                enddo
            enddo

            allocate(temp(nLat,nLat))
            !Check unitarity of matrix
            call ZGEMM('C','N',nLat,nLat,nLat,zone,RotMat,nLat,RotMat,nLat,zzero,temp,nLat) 
            do i = 1,nLat
                do j = 1,nLat
                    if((i.eq.j).and.(abs(temp(i,j)-zone).gt.1.0e-7_dp)) then
                        write(6,*) "i,j: ",i,j
                        call writematrixcomp(temp,'Identity?',.true.)
                        call stop_all(t_r,'Rotation matrix not unitary')
                    elseif((i.ne.j).and.(abs(temp(j,i)).gt.1.0e-7_dp)) then
                        write(6,*) "i,j: ",i,j
                        call writematrixcomp(temp,'Identity?',.true.)
                        call stop_all(t_r,'Rotation matrix not unitary 2')
                    endif
                enddo
            enddo
            !Try other way...
            call ZGEMM('N','C',nLat,nLat,nLat,zone,RotMat,nLat,RotMat,nLat,zzero,temp,nLat) 
            do i = 1,nLat
                do j = 1,nLat
                    if((i.eq.j).and.(abs(temp(i,j)-zone).gt.1.0e-7_dp)) then
                        write(6,*) "i,j: ",i,j
                        call writematrixcomp(temp,'Identity?',.true.)
                        call stop_all(t_r,'Rotation matrix not unitary')
                    elseif((i.ne.j).and.(abs(temp(j,i)).gt.1.0e-7_dp)) then
                        write(6,*) "i,j: ",i,j
                        call writematrixcomp(temp,'Identity?',.true.)
                        call stop_all(t_r,'Rotation matrix not unitary 2')
                    endif
                enddo
            enddo
            write(6,*) "Rotation matrix unitary... :)"

            allocate(CompHam(nLat,nLat))
            allocate(CompHam_2(nLat,nLat))
            allocate(RotHam(nLat,nLat))
            !Store the hamiltonian in complex form, so we can act on it
            do i = 1,nLat
                do j = 1,nLat
                    CompHam(j,i) = dcmplx(Ham(j,i),0.0_dp)
                    CompHam_2(j,i) = dcmplx(Ham(j,i),0.0_dp)
                enddo
            enddo

            !Now, simply transform the hamiltonian to k-space brute force
            call ZGEMM('N','N',nLat,nLat,nLat,zone,RotMat,nLat,CompHam,nLat,zzero,temp,nLat)
            call ZGEMM('N','C',nLat,nLat,nLat,zone,temp,nLat,RotMat,nLat,zzero,RotHam,nLat)

            if(tWriteOut) call writematrixcomp(RotHam,'k-space hamiltonian',.true.)
            
            !Is this now block diagonal?
            do ki = 1,nKPnts
                do bandi = 1,SS_Period
                    Ind_i = SS_Period*(ki-1) + bandi
                    do kj = 1,nKpnts
                        do bandj = 1,SS_Period
                            Ind_j = SS_Period*(kj-1) + bandj

                            if((ki.ne.kj).and.(abs(RotHam(Ind_j,Ind_i)).gt.1.0e-7_dp)) then
                                call stop_all(t_r,'Translational symmetry not conserved')
                            endif
                            if((Ind_j.eq.Ind_i).and.(aimag(RotHam(Ind_j,Ind_i)).gt.1.0e-7_dp)) then
                                call stop_all(t_r,'k-space hamiltonian not hermitian')
                            elseif((Ind_j.ne.Ind_i).and.(abs(RotHam(Ind_j,Ind_i)-dconjg(RotHam(Ind_i,Ind_j))).gt.1.0e-7_dp)) then
                                call stop_all(t_r,'k-space hamiltonian not hermitian')
                            endif
                        enddo
                    enddo
                enddo
            enddo

            !Now, bizarrely, just diagonalize the whole thing!
            CompHam(:,:) = RotHam(:,:)
            lWork = max(1,2*nLat-1)
            allocate(cWork(lWork))
            allocate(Work(max(1,3*nLat-2)))
            info = 0
            call ZHEEV('V','U',nLat,RotHam,nLat,Vals,cWork,lWork,Work,info)
            if(info.ne.0) call stop_all(t_r,'Diag failed')
            deallocate(cWork,Work)
                        
            if(tWriteOut) call writematrixcomp(RotHam,'k-space Eigenvectors',.true.)

            !Are the eigenvectors meaningful?
            allocate(cWork(nLat))
            do i = 1,nLat
                call ZGEMV('N',nLat,nLat,zone,CompHam,nLat,RotHam(:,i),1,zzero,cWork,1)
                do j = 1,nLat
                    if(abs(cWork(j)-(Vals(i)*RotHam(j,i))).gt.1.0e-8_dp) then
                        call stop_all(t_r,'Eigensystem not computed correctly')
                    endif
                enddo
            enddo
            deallocate(cWork)
    
            !Now, rotate eigenvectors back into real-space basis
            !call ZGEMM('C','N',nLat,nLat,nLat,zone,RotMat,nLat,RotHam,nLat,zzero,temp,nLat)
            !call ZGEMM('N','N',nLat,nLat,nLat,zone,temp,nLat,RotMat,nLat,zzero,RotHam,nLat)
            do i = 1,nLat
                !Rotate each eigenvector in turn
                call ZGEMV('C',nLat,nLat,zone,RotMat,nLat,RotHam(:,i),1,zzero,temp(:,i),1)
            enddo

            if(tWriteOut) then
                !Temp now contains the (complex) eigenvectors in r-space
                call writevector(Vals,'Eigenvalues')
                call writematrixcomp(temp,'r-space Eigenvectors',.true.)
            endif
            
            !Are the eigenvectors meaningful in real space compared to original hamiltonian?
            allocate(cWork(nLat))
            do i = 1,nLat
                call ZGEMV('N',nLat,nLat,zone,CompHam_2,nLat,temp(:,i),1,zzero,cWork,1)
                do j = 1,nLat
                    if(abs(cWork(j)-(Vals(i)*temp(j,i))).gt.1.0e-8_dp) then
                        call stop_all(t_r,'Eigensystem not computed correctly in real basis')
                    endif
                enddo
            enddo
            deallocate(cWork)

            do i = 1,nLat
                phase = 0.0_dp
                if(tWriteOut) write(6,*) "Rotating eigenvector : ",i
                do j = 1,nLat
                    if((abs(aimag(temp(j,i))).gt.1.0e-9_dp).and.(abs(temp(j,i)).gt.1.0e-7_dp)) then
                        !Find the phase factor for this eigenvector
                        phase = atan(aimag(temp(j,i))/real(temp(j,i)))
                        if(tWriteOut) write(6,*) "Eigenvector: ",i,j,phase,temp(j,i) * exp(dcmplx(0.0_dp,-phase))
                        !temp(j,i) = temp(j,i) * exp(dcmplx(0.0_dp,-phase))
                        !exit
                    endif
                enddo
!                if(j.le.nLat) then
!                    write(6,*) "phase computed from component: ",j,abs(temp(j,i))
!                    write(6,*) "phase computed to be: ",phase
!                else
!                    write(6,*) "No imaginary component in eigenvector found - no phase being applied"
!                    if(phase.ne.0.0_dp) call stop_all(t_r,'Error here')
!                endif
!
!                !Apply the same phase to all components of the eigenvector
!                do j = 1,nLat
!                    if(abs(phase-(atan(aimag(temp(j,i))/real(temp(j,i))))).gt.1.0e-8_dp) then
!                        write(6,*) "Eigenvector: ",i," component: ",j
!                        write(6,*) phase,atan(aimag(temp(j,i))/real(temp(j,i))),abs(phase-(atan(aimag(temp(j,i))/real(temp(j,i)))))
!                        write(6,*) temp(j,i)
!                        call stop_all(t_r,'A consistent phase cannot be found for this eigenvector')
!                    endif
!                    temp(j,i) = temp(j,i) * exp(dcmplx(0.0_dp,-phase))
!                enddo
            enddo

            !call stop_all(t_r,'stop')
            
            if(tWriteOut) call writematrixcomp(temp,'rotated r-space Eigenvectors',.true.)
            
            !Are the eigenvectors still meaningful?
            allocate(cWork(nLat))
            !Get back the r-space hamiltonian
            do i = 1,nLat
                do j = 1,nLat
                    CompHam(j,i) = dcmplx(Ham(j,i),0.0_dp)
                enddo
            enddo
            do i = 1,nLat
                call ZGEMV('N',nLat,nLat,zone,CompHam,nLat,temp(:,i),1,zzero,cWork,1)
                do j = 1,nLat
                    if(abs(cWork(j)-(Vals(i)*temp(j,i))).gt.1.0e-8_dp) then
                        call writevectorcomp(cWork,'H x vec')
                        call writevectorcomp(Vals(i)*temp(:,i),'RHS')
                        call stop_all(t_r,'Eigensystem not computed correctly 2')
                    endif
                enddo
            enddo
            deallocate(cWork)

            !Are the eigenvectors real now?
            Ham(:,:) = zero
            do i = 1,nLat
                do j = 1,nLat
                    if(aimag(temp(j,i)).gt.1.0e-8_dp) then
                        call writematrixcomp(temp,'Eigenvectors',.true.)
                        call stop_all(t_r,'Eigenvectors complex...')
                    endif
                    Ham(j,i) = real(temp(j,i),dp)
                enddo
            enddo

            deallocate(CompHam,temp,RotMat,RotHam)
        else
            !Normal real space diagonalization
            Vals(:) = 0.0_dp
            allocate(Work(1))
            lWork=-1
            info=0
            call dsyev('V','L',nLat,Ham,nLat,Vals,Work,lWork,info)
            if(info.ne.0) call stop_all(t_r,'Workspace queiry failed')
            lwork=int(work(1))+1
            deallocate(work)
            allocate(work(lwork))
            call dsyev('V','L',nLat,Ham,nLat,Vals,Work,lWork,info)
            if(info.ne.0) call stop_all(t_r,'Diag failed')
            deallocate(work)
!            call writevector(Vals,'Eigenvalues')
        endif

    end subroutine DiagOneEOp

    !Convert a real-space symmetric operator into a k-space operator.
    !In: The operator in real space. Size = nSuperCell x nSupercell
    !   However, it assumes periodicity, such that only the first unit cell, and its connections to the other unit cells
    !   is referenced. I.e. R_Ham(1:nUnitCell,:). In the future, this should be changed so that only this is called.
    !Out: The operator in k-space. In this, the operator is block-diagonal, with each block being of size nUnitCell x nUnitCell.
    subroutine Convert1DtoKSpace(R_Ham,nSuperCell,nUnitCell,k_Ham)
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
