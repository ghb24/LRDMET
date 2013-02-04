module solvers
    use const
    use errors, only: stop_all
    use globals
    use mat_tools, only: WriteVector,WriteMatrix
    implicit none

    contains

    !Transform the integrals over the embedded system, and solve for the energy and 1RDM on impurity.
    !tConverged means that the DMET is converged, and this is the last run, so calculate the 2RDM for other properties
    subroutine SolveSystem(tCreate2RDM)
        use utils, only: get_free_unit 
        use DetToolsData, only: FCIDetList,nFCIDet
        implicit none
        logical, intent(in) :: tCreate2RDM
        integer :: pSpaceDim,i,iunit,j,k,l!,Pivots(4),info
        character(len=256) :: cmd
        character(len=128) :: cmd3
        character(len=73) :: cmd2
        character(len=67) :: cmd1
        character(len=6) :: StrPSpace
        real(dp) :: Emb_nElec!,DD_Response,ZerothH(4,4)
!        real(dp), allocatable :: FullH1(:,:),LR_State(:)
!        real(dp) ::DDOT,Overlap
        real(dp) :: Check2eEnergy,trace
        real(dp), allocatable :: HL_2RDM_temp(:,:)
        character(len=*), parameter :: t_r="SolveSystem"

        !Calculate the number of electrons in the embedded system 
        Emb_nElec = 0.0_dp
        do i=1,EmbSize
            Emb_nElec = Emb_nElec + Emb_MF_DM(i,i)
        enddo
        write(6,"(A,I5)") "Number of electrons in full system: ",NEl
        Elec = nint(Emb_nElec)
        write(6,"(A,F10.7,I5)") "Number of electrons in embedded system: ",Emb_nElec,Elec

        if(nSys.gt.EmbSize) call stop_all(t_r,"Error in determining basis")

        if(.not.tCompleteDiag) then
            call WriteFCIDUMP()

            !Solve with Geralds FCI code
            pSpaceDim = 200
            if(U.gt.4.0) pSpaceDim = 400
            if(EmbSize.lt.8) pSpaceDim = 100
            if(EmbSize.ge.12) pSpaceDim = 1500

            write(StrPSpace,'(I6)') pSpaceDim
            
            if(U.gt.6.0) then
                cmd1 = "fci --subspace-dimension=12 --basis=Input --method='FCI' --pspace="
            else
                cmd1 = "fci --subspace-dimension=12 --basis=CoreH --method='FCI' --pspace="
            endif

            do i=1,6
                if(StrPSpace(i:i).ne.' ') exit
            enddo
            if(i.eq.7) call stop_all(t_r,'Error constructing system call')

            cmd2 = trim(cmd1)//trim(adjustl(StrPSpace(i:6)))
            if(tCreate2RDM) then
                cmd3 = " --thr-var=1e-12 --diis-block-size=208333 --save-rdm1='FCI1RDM' " &
                //"--fci-vec='FCIVec' --save-rdm2='FCI2RDM' 'FCIDUMP' > FCI.out"
            else
                cmd3 = " --thr-var=1e-12 --diis-block-size=208333 --save-rdm1='FCI1RDM' " &
                //"--fci-vec='FCIVec' 'FCIDUMP' > FCI.out"
            endif

            cmd = adjustl(cmd2)//trim(adjustl(cmd3))
!            write(6,*) "Command to call FCI code: ",cmd

            !Run FCI program
            call system(cmd)

            !Extract energy value to new file and read this is
            call system("grep 'FCI STATE 1 ENERGY' FCI.out | awk '{print$5}' > FCI.ene")
            HL_Energy = 0.0_dp
            iunit=get_free_unit()
            open(iunit,file='FCI.ene',status='old')
            read(iunit,*) HL_Energy
            close(iunit)
            if(HL_Energy.eq.0.0_dp) call stop_all(t_r,"FCI energy is 0")

            !read in the 1 (and 2) RDMs
            if(allocated(HL_1RDM)) deallocate(HL_1RDM)
            allocate(HL_1RDM(EmbSize,EmbSize))
            HL_1RDM(:,:) = 0.0_dp
            iunit=get_free_unit()
            open(iunit,file='FCI1RDM',status='old')
            read(iunit,*) cmd   !First line is a header
            do i=1,EmbSize
                !We can read in in this order since it should be symmetric
                read(iunit,*) HL_1RDM(:,i)
            enddo
            close(iunit)
!            call writematrix(HL_1RDM,'FCI 1RDM',.true.)

            if(tCreate2RDM) then
                write(6,*) "Creating 2RDM over impurity system..."
                if(allocated(HL_2RDM)) deallocate(HL_2RDM)
                allocate(HL_2RDM(EmbSize,EmbSize,EmbSize,EmbSize))  !< i^+ j^+ k l >
                allocate(HL_2RDM_temp(EmbSize*EmbSize,EmbSize*EmbSize))  !< i^+ j k^+ l >
                HL_2RDM(:,:,:,:) = 0.0_dp
                HL_2RDM_temp(:,:) = 0.0_dp
                iunit = get_free_unit()
                open(iunit,file='FCI2RDM',status='old')
                read(iunit,*) cmd   !First line is a header
                do i=1,EmbSize*EmbSize
                    read(iunit,*) HL_2RDM_temp(:,i)
                enddo
                close(iunit)
                !We now have a matrix, but we want to change it to the format we want.
                !HL_2RDM(i,j,k,l) = G_ij^kl = < i^+ k^+ l j >
                do i=1,EmbSize
                    do j=1,EmbSize
                        do k=1,EmbSize
                            do l=1,EmbSize

                                HL_2RDM(i,j,k,l) = HL_2RDM(i,j,k,l) + HL_2RDM_temp(((k-1)*EmbSize)+l,((i-1)*EmbSize)+j)

                                if(j.eq.k) then
                                    HL_2RDM(i,j,k,l) = HL_2RDM(i,j,k,l) - HL_1RDM(i,l)
                                endif

                            enddo
                        enddo
                    enddo
                enddo
                deallocate(HL_2RDM_temp)

            endif

        else
            !Do a complete diagonalization
            !Do not need to write FCIDUMP, since would only read it back in...
            call CompleteDiag(tCreate2RDM)
        endif
            
        write(6,"(A,F20.10)") "Embedded system energy is: ",HL_Energy

        !Now calculate the seperate 1-body and 2-body contribution to the energy of the embedded system
        !The one-body contribution can be calculated from Tr[h^T x FCI_RDM] = \sum_ij h_ij * FCI_RDM_ij
        !where h^T is the core hamiltonian in the embedded basis with the correlation potential over all but the impurity (the 1e integrals)
        One_ElecE = 0.0_dp
        do j=1,EmbSize
            do i=1,EmbSize
                if((i.eq.1).and.(j.eq.1).and.(tChemPot)) then
                    !Include the chemical potential in the one-electron hamiltonian
                    One_ElecE = One_ElecE + (Emb_h0v(i,j)-(U/2.0_dp))*HL_1RDM(i,j)
                else
                    One_ElecE = One_ElecE + Emb_h0v(i,j)*HL_1RDM(i,j)
                endif
            enddo
        enddo

        !The two electron contribution to the embedded system energy is just the FCI result - the 1e energy
        Two_ElecE = HL_Energy - One_ElecE

        if(tCreate2RDM) then
            !Do some tests to make sure we have the right 2RDM
            do i=1,EmbSize
                do j=1,EmbSize
                    do k=1,EmbSize
                        do l=1,EmbSize
                            if(abs(HL_2RDM(i,j,k,l)-HL_2RDM(k,l,i,j)).gt.1.0e-7_dp) then
                                write(6,*) "RDM(i,j,k,l): ",HL_2RDM(i,j,k,l)
                                write(6,*) "RDM(k,l,i,j): ",HL_2RDM(k,l,i,j)
                                call stop_all(t_r,'2RDM not symmetric')
                            endif
                            if(abs(HL_2RDM(i,j,k,l)-HL_2RDM(j,i,l,k)).gt.1.0e-7_dp) then
                                write(6,*) "RDM(i,j,k,l): ",HL_2RDM(i,j,k,l)
                                write(6,*) "RDM(j,i,l,k): ",HL_2RDM(j,i,l,k)
                                call stop_all(t_r,'2RDM not symmetric')
                            endif
                            if(abs(HL_2RDM(i,j,k,l)-HL_2RDM(l,k,j,i)).gt.1.0e-7_dp) then
                                write(6,*) "RDM(i,j,k,l): ",HL_2RDM(i,j,k,l)
                                write(6,*) "RDM(l,k,j,i): ",HL_2RDM(l,k,j,i)
                                call stop_all(t_r,'2RDM not symmetric')
                            endif
                        enddo
                    enddo
                enddo
            enddo

            !Check that 2RDM is correct by explicitly calculating the two-electron contribution to the FCI energy from the 2RDM
            !We only need to check the iiii components over the impurity sites ONLY
            Check2eEnergy = 0.0_dp
            if(tAnderson) then
                Check2eEnergy = Check2eEnergy + U*HL_2RDM(1,1,1,1)
            else
                do i=1,nImp
                    Check2eEnergy = Check2eEnergy + U*HL_2RDM(i,i,i,i)
                enddo
            endif
            Check2eEnergy = Check2eEnergy / 2.0_dp
            if(abs(Check2eEnergy-Two_ElecE).gt.1.0e-7_dp) then
                write(6,*) "Check2eEnergy: ",Check2eEnergy
                write(6,*) "Two_ElecE: ",Two_ElecE
                call stop_all(t_r,'2RDM calculated incorrectly')
            endif

            !Also check that trace condition is satisfied
            trace = 0.0_dp
            do i=1,EmbSize
                do j=1,EmbSize
                    trace = trace + HL_2RDM(i,i,j,j)
                enddo
            enddo
            if(abs(trace-real(elec*(elec-1),dp)).gt.1.0e-8_dp) then
                write(6,*) "trace: ",trace
                write(6,*) "elec*(elec-1): ",elec*(elec-1)
                call stop_all(t_r,'2RDM trace condition incorrect')
            endif

            !Now check that it reduced correctly down to the 1RDM
            do i=1,EmbSize
                do j=1,EmbSize
                    trace = 0.0_dp
                    do k=1,EmbSize
                        trace = trace + HL_2RDM(j,i,k,k)
                    enddo
                    if(abs(trace-((elec-1)*HL_1RDM(i,j))).gt.1.0e-8_dp) then
                        write(6,*) "Reduced 2RDM component: ",trace
                        write(6,*) "(elec-1)*HL_1RDM(i,j): ",(elec-1)*HL_1RDM(i,j)
                        call stop_all(t_r,'2RDM does not reduce to 1RDM appropriately')
                    endif
                enddo
            enddo
        endif

        !We only want to calculate the energy over the impurity site, along with the coupling to the bath.
        !Calculate one-electron energy contributions only over the impurity
        One_ElecE_Imp = 0.0_dp
        do j=1,nImp
            do i=1,nImp
                if(tChemPot.and.(i.eq.1).and.(j.eq.1)) then
                    One_ElecE_Imp = One_ElecE_Imp + (Emb_h0v(i,j)-U/2.0_dp)*HL_1RDM(i,j)
                else
                    One_ElecE_Imp = One_ElecE_Imp + Emb_h0v(i,j)*HL_1RDM(i,j)
                endif
            enddo
        enddo
        One_ElecE_Imp = One_ElecE_Imp/real(nImp)    !For energy per impurity

        !Two electron energy terms are not contained in core hamiltonian, or expressed over the bath, 
        !so just divided total 2e contribution by number of impurities to get the 2e contrib
        Two_ElecE_Imp = Two_ElecE/real(nImp)

        !There is also finally an interaction term between the bath and impurity which we can calculate
        !Model this as the h0v matrix with added 1/2 correlation potential so that there is no double counting between the interactions with the bath-impurity
        CoupE_Imp = 0.0_dp
        do j=nImp+1,EmbSize
            do i=1,nImp
                CoupE_Imp = CoupE_Imp + (Emb_h0v(i,j) + 0.5_dp*Emb_CorrPot(i,j))*HL_1RDM(i,j)
            enddo
        enddo
        CoupE_Imp = CoupE_Imp / real(nImp)

        TotalE_Imp = One_ElecE_Imp + Two_ElecE_Imp + CoupE_Imp

        write(6,"(A,F10.6)") "One electron energy per impurity:     ",One_ElecE_Imp
        write(6,"(A,F10.6)") "Two electron energy per impurity:     ",Two_ElecE_Imp
        write(6,"(A,F10.6)") "Coupling energy to bath per impurity: ",CoupE_Imp

        write(6,"(A,2F10.6)") "Total energy per impurity site:       ",TotalE_Imp

        !The target filling is the original filling of the mean field RDM per site
        !We know that the filling should be uniform due to the translational symmetry
        Targetfilling_Imp = 0.0_dp
        do i=1,nSites
            Targetfilling_Imp = Targetfilling_Imp + MeanFieldDM(i,i)
        enddo
        Targetfilling_Imp = Targetfilling_Imp / (2.0_dp*real(nSites))

        !Now what is the actually filling on each impurity from the embedded system
        Actualfilling_Imp = 0.0_dp
        do i=1,nImp
            Actualfilling_Imp = Actualfilling_Imp + HL_1RDM(i,i)
        enddo
        Actualfilling_Imp = Actualfilling_Imp / (2.0_dp*real(nImp))
        Fillingerror = Actualfilling_Imp - Targetfilling_Imp

        write(6,"(A,F15.7)") "Target filling per site: ",Targetfilling_Imp
        write(6,"(A,F15.7)") "Actual filling per site: ",Actualfilling_Imp
        write(6,"(A,F15.7)") "Filling error  per site: ",Fillingerror

!        if(tDebug) call writematrix(HL_1RDM,'hl_1rdm',.true.)
        if(tWriteOut) then
            call writematrix(HL_1RDM,'hl_1rdm',.true.)
        endif

    end subroutine SolveSystem

    !Calculate the FCI result for this hubbard model
    subroutine DiagFullSystem()
        use utils, only: get_free_unit
        implicit none
        integer :: pSpaceDim,i,j,iunit
        character(len=256) :: cmd
        character(len=128) :: cmd3
        character(len=73) :: cmd2
        character(len=67) :: cmd1
        character(len=6) :: StrPSpace
        character(len=*), parameter :: t_r="SolveFullSystem"

        write(6,"(A,I5)") "Number of electrons in full system: ",NEl
        if(.not.tCompleteDiag) then
        
            iunit = get_free_unit()
            open(iunit,file='FCIDUMP',status='unknown')
            write(iunit,*) "&FCI NORB=",nSites," ,NELEC=",NEl," ,MS2=0,"
            write(iunit,"(A)",advance='no') "ORBSYM= 1,"
            do i=2,nSites-1
                write(iunit,"(A)",advance='no') "1,"
            enddo
            write(iunit,"(A)") "1,"
            write(iunit,"(A)") "ISYM=1"
            write(iunit,"(A)") "&END"

            !Just define diagonal 2 electron contribution over impurity sites
            if(tAnderson) then
                write(iunit,"(F16.12,4I8)") U,1,1,1,1
            else
                do i=1,nSites
                    write(iunit,"(F16.12,4I8)") U,i,i,i,i
                enddo
            endif

            !Now for 1electron contributions
            do i=1,nSites
                do j=1,i
                    if(abs(h0(i,j)).gt.1.0e-10_dp) then
                        write(iunit,"(F16.12,4I8)") h0(i,j),i,j,0,0
                    endif
                enddo
            enddo
            !Core energy is zero
            write(iunit,"(F16.12,4I8)") 0.0_dp,0,0,0,0
            close(iunit)

            !Solve with Geralds FCI code
            pSpaceDim = 200
            if(U.gt.4.0) pSpaceDim = 400
            if(nSites.lt.8) pSpaceDim = 100
            if(nSites.ge.12) pSpaceDim = 1500

            write(StrPSpace,'(I6)') pSpaceDim
            
            if(U.gt.6.0) then
                cmd1 = "fci --subspace-dimension=12 --basis=Input --method='FCI' --pspace="
            else
                cmd1 = "fci --subspace-dimension=12 --basis=CoreH --method='FCI' --pspace="
            endif

            do i=1,6
                if(StrPSpace(i:i).ne.' ') exit
            enddo
            if(i.eq.7) call stop_all(t_r,'Error constructing system call')

            cmd2 = trim(cmd1)//trim(adjustl(StrPSpace(i:6)))
            cmd3 = " --thr-var=1e-12 --diis-block-size=208333 --save-rdm1='FCI1RDM' " &
            //"--fci-vec='FCIVec' 'FCIDUMP' > FCI.out"

            cmd = adjustl(cmd2)//trim(adjustl(cmd3))
!            write(6,*) "Command to call FCI code: ",cmd

            !Run FCI program
            call system(cmd)

            !Extract energy value to new file and read this is
            call system("grep 'FCI STATE 1 ENERGY' FCI.out | awk '{print$5}' > FCI.ene")
            HL_Energy = 0.0_dp
            iunit=get_free_unit()
            open(iunit,file='FCI.ene',status='old')
            read(iunit,*) HL_Energy
            close(iunit)
            if(HL_Energy.eq.0.0_dp) call stop_all(t_r,"FCI energy is 0")

            write(6,"(A,F20.10)") "Complete system energy is: ",HL_Energy

        else
            write(6,*) "Cannot currently diagonalize full system with complete diagonalizer. Use DAVIDSON."
            write(6,*) "Skipping diagonalization of full system (though won't be difficult to code up)."
        endif

    end subroutine DiagFullSystem

!Generate all determinants in the FCI space and do complete diagonalization
    subroutine CompleteDiag(tCreate2RDM)
        use DetToolsData
        implicit none
        logical, intent(in) :: tCreate2RDM
        integer :: OrbPairs,UMatSize,umatind
        real(dp), allocatable :: work(:)
        integer :: lwork,info,i,j
        character(len=*), parameter :: t_r='CompleteDiag'

        !First, allocate and fill the umat and tmat for the FCI space
        OrbPairs = (EmbSize*(EmbSize+1))/2
        UMatSize = (OrbPairs*(OrbPairs+1))/2
        write(6,*) "Allocating memory to store 2 electron integrals: ",UMatSize
        if(allocated(UMat)) deallocate(UMat)
        allocate(UMat(UMatSize))
        UMat(:) = 0.0_dp
        if(tAnderson) then
            umat(umatind(1,1,1,1)) = U
        else
            do i=1,nImp
                umat(umatind(i,i,i,i)) = U
            enddo
        endif
        if(allocated(tmat)) deallocate(tmat)
        allocate(tmat(EmbSize,EmbSize))
        tmat(:,:) = 0.0_dp
        do i=1,EmbSize
            do j=1,EmbSize
                if(abs(Emb_h0v(i,j)).gt.1.0e-10_dp) then
                    tmat(i,j) = Emb_h0v(i,j)
                endif
            enddo
        enddo
        if(tChemPot) then
            tmat(1,1) = tmat(1,1) - U/2.0_dp
        endif
        if(tWriteOut) then
            call writematrix(tmat,'tmat',.true.)
        endif

        !Now generate all determinants in the active space
        if(allocated(FCIDetList)) deallocate(FCIDetList)
        call GenDets(Elec,EmbSize,.false.,.false.,.false.)
        !FCIDetList now stores a list of all the determinants
        write(6,"(A,I14)") "Number of determinants in FCI space: ",nFCIDet
        write(6,"(A,F14.6,A)") "Allocating memory for the hamiltonian: ",real((nFCIDet**2)*8,dp)/1048576.0_dp," Mb"
        if(allocated(FullHamil)) deallocate(FullHamil)
        if(allocated(Spectrum)) deallocate(Spectrum)
        allocate(Spectrum(nFCIDet))
        allocate(FullHamil(nFCIDet,nFCIDet))
        FullHamil(:,:) = 0.0_dp
        !Construct the hamiltonian - slow
        do i=1,nFCIDet
            do j=1,nFCIDet
                call GetHElement(FCIDetList(:,i),FCIDetList(:,j),Elec,FullHamil(i,j))
            enddo
        enddo
!        call writematrix(FullHamil(1:nFCIDet,1:nFCIDet),'FCI hamil',.true.)

        !Diagonalize
        allocate(Work(1))
        lWork=-1
        info=0
        call dsyev('V','U',nFCIDet,FullHamil,nFCIDet,Spectrum,Work,lWork,info)
        if(info.ne.0) call stop_all(t_r,'Workspace queiry failed')
        lwork=int(work(1))+1
        deallocate(work)
        allocate(work(lwork))
        call dsyev('V','U',nFCIDet,FullHamil,nFCIDet,Spectrum,Work,lWork,info)
        if(info.ne.0) call stop_all(t_r,'Diag failed')
        deallocate(work)

        call writevector(Spectrum(1:(min(10,nFCIDet))),'FCI Spectrum')
            
        HL_Energy = Spectrum(1) !GS eigenvalue
            
        if(allocated(HL_1RDM)) deallocate(HL_1RDM)
        allocate(HL_1RDM(EmbSize,EmbSize))
        HL_1RDM(:,:) = 0.0_dp

        call FindFull1RDM(1,1,HL_1RDM)

        if(tWriteOut) call writematrix(HL_1RDM,'HL_1RDM',.true.)

        if(tCreate2RDM) then
            if(allocated(HL_2RDM)) deallocate(HL_2RDM)
            allocate(HL_2RDM(EmbSize,EmbSize,EmbSize,EmbSize))
            HL_2RDM(:,:,:,:) = 0.0_dp
            call FindFull2RDM(1,1,HL_2RDM)
        endif

    end subroutine CompleteDiag

    subroutine FindFull1RDM(StateBra,StateKet,RDM)
        use DetToolsData, only: FCIDetList,nFCIDet
        implicit none
        real(dp) , intent(out) :: RDM(EmbSize,EmbSize)
        integer , intent(in) :: StateBra,StateKet
        integer :: Ex(2),gtid,i,j,k,IC,iGetExcitLevel
        logical :: tSign
        character(len=*), parameter :: t_r='FindFull1RDM'

        RDM(:,:) = 0.0_dp

        do i=1,nFCIDet
            do j=1,nFCIDet
                IC = iGetExcitLevel(FCIDetList(:,i),FCIDetList(:,j),Elec)
                if(IC.eq.1) then
                    !Connected by a single
                    Ex(1) = 1
                    call GetExcitation(FCIDetList(:,i),FCIDetList(:,j),Elec,Ex,tSign)
                    if(tSign) then
                        RDM(gtid(Ex(1)),gtid(Ex(2))) = RDM(gtid(Ex(1)),gtid(Ex(2))) -   &
                            FullHamil(i,StateBra)*FullHamil(j,StateKet)
                    else
                        RDM(gtid(Ex(1)),gtid(Ex(2))) = RDM(gtid(Ex(1)),gtid(Ex(2))) +   &
                            FullHamil(i,StateBra)*FullHamil(j,StateKet)
                    endif
                elseif(IC.eq.0) then
                    !Same det
                    if(i.ne.j) call stop_all(t_r,'Error here')
                    do k=1,Elec
                        RDM(gtid(FCIDetList(k,i)),gtid(FCIDetList(k,i))) = RDM(gtid(FCIDetList(k,i)),gtid(FCIDetList(k,i))) &
                            + FullHamil(i,StateBra)*FullHamil(j,StateKet)
                    enddo
                endif

            enddo
        enddo

        do i=1,EmbSize
            do j=1,EmbSize
                if(abs(RDM(i,j)-RDM(j,i)).gt.1.0e-7_dp) call stop_all(t_r,'1RDM not symmetric')
            enddo
        enddo

    end subroutine FindFull1RDM

    !Find spin-integrated 2RDM very inefficiently
    !According to Helgakker <0|e_pqrs|0>
    !Done by running through all N^2 determinant pairs
    subroutine FindFull2RDM(StateBra,StateKet,RDM)
        use DetToolsData, only: FCIDetList,nFCIDet
        implicit none
        real(dp) , intent(out) :: RDM(EmbSize,EmbSize,EmbSize,EmbSize)
        integer , intent(in) :: StateBra,StateKet
        integer :: Ex(2,2),gtid,i,j,k,IC,iGetExcitLevel,kel,lel,l,temp
        logical :: tSign
        character(len=*), parameter :: t_r='FindFull2RDM'

        RDM(:,:,:,:) = 0.0_dp

        do i=1,nFCIDet
            do j=1,nFCIDet
                IC = iGetExcitLevel(FCIDetList(:,i),FCIDetList(:,j),Elec)
                if(IC.eq.2) then
                    !Connected by a double
                    Ex(1,1) = 2
                    call GetExcitation(FCIDetList(:,i),FCIDetList(:,j),Elec,Ex,tSign)
                    if(mod(Ex(1,1),2).ne.mod(Ex(1,2),2)) then
                        !We have a mixed spin excitation
                        !Ensure that the spin of i is the same as the spin of b
                        !If its not, then reverse a and b and flip the sign
                        if(mod(Ex(1,1),2).ne.mod(Ex(2,2),2)) then
                            temp = Ex(2,2)
                            Ex(2,2) = Ex(2,1)
                            Ex(2,1) = temp
                            tSign = .not.tSign
                        endif
                    endif
                    if(tSign) then
                        if(mod(Ex(1,1),2).eq.mod(Ex(1,2),2)) then
                            !same spin excitation
                            RDM(gtid(Ex(2,1)),gtid(Ex(1,2)),gtid(Ex(2,2)),gtid(Ex(1,1))) =  &
                                RDM(gtid(Ex(2,1)),gtid(Ex(1,2)),gtid(Ex(2,2)),gtid(Ex(1,1))) +  &
                                FullHamil(i,StateBra)*FullHamil(j,StateKet)
                            RDM(gtid(Ex(2,2)),gtid(Ex(1,1)),gtid(Ex(2,1)),gtid(Ex(1,2))) =  &
                                RDM(gtid(Ex(2,2)),gtid(Ex(1,1)),gtid(Ex(2,1)),gtid(Ex(1,2))) +  &
                                FullHamil(i,StateBra)*FullHamil(j,StateKet)
                            RDM(gtid(Ex(2,2)),gtid(Ex(1,2)),gtid(Ex(2,1)),gtid(Ex(1,1))) =  &
                                RDM(gtid(Ex(2,2)),gtid(Ex(1,2)),gtid(Ex(2,1)),gtid(Ex(1,1))) -  &
                                FullHamil(i,StateBra)*FullHamil(j,StateKet)
                            RDM(gtid(Ex(2,1)),gtid(Ex(1,1)),gtid(Ex(2,2)),gtid(Ex(1,2))) =  &
                                RDM(gtid(Ex(2,1)),gtid(Ex(1,1)),gtid(Ex(2,2)),gtid(Ex(1,2))) -  &
                                FullHamil(i,StateBra)*FullHamil(j,StateKet)

                        else
                            !Mixed spin excitation
                            !i has the same spin as b
                            RDM(gtid(Ex(2,1)),gtid(Ex(1,2)),gtid(Ex(2,2)),gtid(Ex(1,1))) = &
                                RDM(gtid(Ex(2,1)),gtid(Ex(1,2)),gtid(Ex(2,2)),gtid(Ex(1,1))) + &
                                FullHamil(i,StateBra)*FullHamil(j,StateKet)
                            RDM(gtid(Ex(2,2)),gtid(Ex(1,1)),gtid(Ex(2,1)),gtid(Ex(1,2))) = &
                                RDM(gtid(Ex(2,2)),gtid(Ex(1,1)),gtid(Ex(2,1)),gtid(Ex(1,2))) + &
                                FullHamil(i,StateBra)*FullHamil(j,StateKet)

                        endif

                    else
                        if(mod(Ex(1,1),2).eq.mod(Ex(1,2),2)) then
                            !same spin excitation
                            RDM(gtid(Ex(2,1)),gtid(Ex(1,2)),gtid(Ex(2,2)),gtid(Ex(1,1))) =  &
                                RDM(gtid(Ex(2,1)),gtid(Ex(1,2)),gtid(Ex(2,2)),gtid(Ex(1,1))) -  &
                                FullHamil(i,StateBra)*FullHamil(j,StateKet)
                            RDM(gtid(Ex(2,2)),gtid(Ex(1,1)),gtid(Ex(2,1)),gtid(Ex(1,2))) =  &
                                RDM(gtid(Ex(2,2)),gtid(Ex(1,1)),gtid(Ex(2,1)),gtid(Ex(1,2))) -  &
                                FullHamil(i,StateBra)*FullHamil(j,StateKet)
                            RDM(gtid(Ex(2,2)),gtid(Ex(1,2)),gtid(Ex(2,1)),gtid(Ex(1,1))) =  &
                                RDM(gtid(Ex(2,2)),gtid(Ex(1,2)),gtid(Ex(2,1)),gtid(Ex(1,1))) +  &
                                FullHamil(i,StateBra)*FullHamil(j,StateKet)
                            RDM(gtid(Ex(2,1)),gtid(Ex(1,1)),gtid(Ex(2,2)),gtid(Ex(1,2))) =  &
                                RDM(gtid(Ex(2,1)),gtid(Ex(1,1)),gtid(Ex(2,2)),gtid(Ex(1,2))) +  &
                                FullHamil(i,StateBra)*FullHamil(j,StateKet)
                        else
                            !Mixed spin excitation
                            !i has the same spin as b
                            RDM(gtid(Ex(2,1)),gtid(Ex(1,2)),gtid(Ex(2,2)),gtid(Ex(1,1))) = &
                                RDM(gtid(Ex(2,1)),gtid(Ex(1,2)),gtid(Ex(2,2)),gtid(Ex(1,1))) - &
                                FullHamil(i,StateBra)*FullHamil(j,StateKet)
                            RDM(gtid(Ex(2,2)),gtid(Ex(1,1)),gtid(Ex(2,1)),gtid(Ex(1,2))) = &
                                RDM(gtid(Ex(2,2)),gtid(Ex(1,1)),gtid(Ex(2,1)),gtid(Ex(1,2))) - &
                                FullHamil(i,StateBra)*FullHamil(j,StateKet)
                        endif
                    endif

                elseif(IC.eq.1) then
                    !Connected by a single
                    Ex(1,1) = 1
                    call GetExcitation(FCIDetList(:,i),FCIDetList(:,j),Elec,Ex,tSign)
                    do k=1,elec
                        kel = gtid(FCIDetList(k,i))

                        if(FCIDetList(k,i).ne.Ex(1,1)) then
                            if(tSign) then
                                if(mod(Ex(1,1),2).eq.mod(FCIDetList(k,i),2)) then
                                    !k also same spin
                                    RDM(kel,gtid(Ex(1,1)),gtid(Ex(2,1)),kel) = &
                                        RDM(kel,gtid(Ex(1,1)),gtid(Ex(2,1)),kel) + &
                                        FullHamil(i,StateBra)*FullHamil(j,StateKet)

                                    RDM(gtid(Ex(2,1)),gtid(Ex(1,1)),kel,kel) = &
                                        RDM(gtid(Ex(2,1)),gtid(Ex(1,1)),kel,kel) - &
                                        FullHamil(i,StateBra)*FullHamil(j,StateKet)

                                    RDM(kel,kel,gtid(Ex(2,1)),gtid(Ex(1,1))) = &
                                        RDM(kel,kel,gtid(Ex(2,1)),gtid(Ex(1,1))) - &
                                        FullHamil(i,StateBra)*FullHamil(j,StateKet)

                                    RDM(gtid(Ex(2,1)),kel,kel,gtid(Ex(1,1))) = &
                                        RDM(gtid(Ex(2,1)),kel,kel,gtid(Ex(1,1))) + &
                                        FullHamil(i,StateBra)*FullHamil(j,StateKet)

                                else
                                    !k opposite spin
                                    RDM(gtid(Ex(2,1)),gtid(Ex(1,1)),kel,kel) = &
                                        RDM(gtid(Ex(2,1)),gtid(Ex(1,1)),kel,kel) - &
                                        FullHamil(i,StateBra)*FullHamil(j,StateKet)
                                    RDM(kel,kel,gtid(Ex(2,1)),gtid(Ex(1,1))) = &
                                        RDM(kel,kel,gtid(Ex(2,1)),gtid(Ex(1,1))) - &
                                        FullHamil(i,StateBra)*FullHamil(j,StateKet)
                                endif


                            else
                                if(mod(Ex(1,1),2).eq.mod(FCIDetList(k,i),2)) then
                                    !k also same spin
                                    RDM(kel,gtid(Ex(1,1)),gtid(Ex(2,1)),kel) = &
                                        RDM(kel,gtid(Ex(1,1)),gtid(Ex(2,1)),kel) - &
                                        FullHamil(i,StateBra)*FullHamil(j,StateKet)

                                    RDM(gtid(Ex(2,1)),gtid(Ex(1,1)),kel,kel) = &
                                        RDM(gtid(Ex(2,1)),gtid(Ex(1,1)),kel,kel) + &
                                        FullHamil(i,StateBra)*FullHamil(j,StateKet)

                                    RDM(kel,kel,gtid(Ex(2,1)),gtid(Ex(1,1))) = &
                                        RDM(kel,kel,gtid(Ex(2,1)),gtid(Ex(1,1))) + &
                                        FullHamil(i,StateBra)*FullHamil(j,StateKet)

                                    RDM(gtid(Ex(2,1)),kel,kel,gtid(Ex(1,1))) = &
                                        RDM(gtid(Ex(2,1)),kel,kel,gtid(Ex(1,1))) - &
                                        FullHamil(i,StateBra)*FullHamil(j,StateKet)
                                else
                                    !k opposite spin
                                    RDM(gtid(Ex(2,1)),gtid(Ex(1,1)),kel,kel) = &
                                        RDM(gtid(Ex(2,1)),gtid(Ex(1,1)),kel,kel) + &
                                        FullHamil(i,StateBra)*FullHamil(j,StateKet)
                                    RDM(kel,kel,gtid(Ex(2,1)),gtid(Ex(1,1))) = &
                                        RDM(kel,kel,gtid(Ex(2,1)),gtid(Ex(1,1))) + &
                                        FullHamil(i,StateBra)*FullHamil(j,StateKet)
                                endif


                            endif

                        endif
                    enddo
                elseif(IC.eq.0) then
                    !Same det
                    Ex(1,1) = 0
                    if(i.ne.j) call stop_all(t_r,'Error here')
                    do k=1,Elec
                        kel = gtid(FCIDetList(k,i))
                        do l=k+1,Elec
                            lel = gtid(FCIDetList(l,i))
                            if(FCIDetList(k,i).eq.FCIDetList(l,i)) cycle

                            if(mod(FCIDetList(l,i),2).eq.mod(FCIDetList(k,i),2)) then
                                RDM(kel,kel,lel,lel) = RDM(kel,kel,lel,lel) + &
                                    FullHamil(i,StateBra)*FullHamil(j,StateKet)
                                RDM(lel,lel,kel,kel) = RDM(lel,lel,kel,kel) + &
                                    FullHamil(i,StateBra)*FullHamil(j,StateKet)
                                RDM(lel,kel,kel,lel) = RDM(lel,kel,kel,lel) - &
                                    FullHamil(i,StateBra)*FullHamil(j,StateKet)
                                RDM(kel,lel,lel,kel) = RDM(kel,lel,lel,kel) - &
                                    FullHamil(i,StateBra)*FullHamil(j,StateKet)
                            else
                                RDM(kel,kel,lel,lel) = RDM(kel,kel,lel,lel) + &
                                    FullHamil(i,StateBra)*FullHamil(j,StateKet)
                                RDM(lel,lel,kel,kel) = RDM(lel,lel,kel,kel) + &
                                    FullHamil(i,StateBra)*FullHamil(j,StateKet)
                            endif

                        enddo
                    enddo
                endif

            enddo
        enddo


    end subroutine FindFull2RDM

    !This routine is horrifically written! Rewrite!
    subroutine CompleteDiag_old()
        implicit none
        real(dp) , allocatable :: work(:)
        integer :: lWork,info
        character(len=*), parameter :: t_r='CompleteDiag_old'

        write(6,*) "Entering complete diagonalization"

        if(allocated(FullHamil)) deallocate(FullHamil)
        if(allocated(Spectrum)) deallocate(Spectrum)

        !Det 1: 1_alpha 1_beta
        !Det 2: 1_alpha 2_beta
        !Det 3: 2_alpha 1_beta
        !Det 4: 2_alpha 2_beta
        allocate(FullHamil(4,4))
        allocate(Spectrum(4))

        FullHamil(:,:) = 0.0_dp
        !Diagonals first
        FullHamil(1,1) = 2.0_dp*Emb_h0v(1,1) + U
        FullHamil(2,2) = Emb_h0v(1,1) + Emb_h0v(2,2)
        FullHamil(3,3) = Emb_h0v(1,1) + Emb_h0v(2,2)
        FullHamil(4,4) = 2.0_dp*Emb_h0v(2,2) 

        FullHamil(1,2) = Emb_h0v(1,2)
        FullHamil(1,3) = - Emb_h0v(1,2)

        FullHamil(2,4) = Emb_h0v(1,2)
        FullHamil(3,4) = - Emb_h0v(1,2)
        
        FullHamil(2,1) = Emb_h0v(1,2)
        FullHamil(4,1) = - Emb_h0v(1,2)

        FullHamil(4,2) = Emb_h0v(1,2)
        FullHamil(4,3) = - Emb_h0v(1,2)

        allocate(Work(1))
        lWork=-1
        info=0
        call dsyev('V','U',4,FullHamil,4,Spectrum,Work,lWork,info)
        if(info.ne.0) call stop_all(t_r,'Workspace queiry failed')
        lwork=int(work(1))+1
        deallocate(work)
        allocate(work(lwork))
        call dsyev('V','U',4,FullHamil,4,Spectrum,Work,lWork,info)
        if(info.ne.0) call stop_all(t_r,'Diag failed')
        deallocate(work)

        call writevector(Spectrum,'FCI Spectrum')
            
        HL_Energy = Spectrum(1) !GS eigenvalue
            
        if(allocated(HL_1RDM)) deallocate(HL_1RDM)
        allocate(HL_1RDM(EmbSize,EmbSize))
        HL_1RDM(:,:) = 0.0_dp

        call FindFull1RDM_old(1,1,HL_1RDM)

        call writematrix(HL_1RDM,'HL_1RDM',.true.)

    end subroutine CompleteDiag_old

    subroutine FindFull1RDM_old(StateBra,StateKet,RDM)
        implicit none
        real(dp) , intent(out) :: RDM(2,2)
        integer , intent(in) :: StateBra,StateKet

        RDM(:,:) = 0.0_dp

        !First, diagonal contributions
        RDM(1,1) = RDM(1,1) + FullHamil(1,StateBra)*FullHamil(1,StateKet)*2.0_dp !alpha alpha and beta beta component of det 1
        RDM(1,1) = RDM(1,1) + FullHamil(2,StateBra)*FullHamil(2,StateKet)    
        RDM(2,2) = RDM(2,2) + FullHamil(2,StateBra)*FullHamil(2,StateKet)
        RDM(1,1) = RDM(1,1) + FullHamil(3,StateBra)*FullHamil(3,StateKet)
        RDM(2,2) = RDM(2,2) + FullHamil(3,StateBra)*FullHamil(3,StateKet)
        RDM(2,2) = RDM(2,2) + FullHamil(4,StateBra)*FullHamil(4,StateKet)*2.0_dp !alpha alpha and beta beta component of det 4

        !Now the four off diagonal contributions    1 -> 2
        RDM(2,1) = RDM(2,1) + FullHamil(2,StateBra)*FullHamil(1,StateKet)
        RDM(1,2) = RDM(1,2) + FullHamil(1,StateBra)*FullHamil(2,StateKet)   ! 2 -> 1

        ! 1 -> 3 (Note minus sign)
        RDM(2,1) = RDM(2,1) - FullHamil(3,StateBra)*FullHamil(1,StateKet)   
        RDM(1,2) = RDM(1,2) - FullHamil(1,StateBra)*FullHamil(3,StateKet)

        ! 2 -> 4
        RDM(2,1) = RDM(2,1) + FullHamil(4,StateBra)*FullHamil(2,StateKet)
        RDM(1,2) = RDM(1,2) + FullHamil(2,StateBra)*FullHamil(4,StateKet)

        !Finally 3 -> 4 (also minus sign)
        RDM(2,1) = RDM(2,1) - FullHamil(4,StateBra)*FullHamil(3,StateKet)
        RDM(1,2) = RDM(1,2) - FullHamil(3,StateBra)*FullHamil(4,StateKet)

    end subroutine FindFull1RDM_old

    subroutine WriteFCIDUMP()
        use utils, only: get_free_unit
        implicit none
        integer :: iunit,i,j

        iunit = get_free_unit()
        open(iunit,file='FCIDUMP',status='unknown')
        write(iunit,*) "&FCI NORB=",EmbSize," ,NELEC=",Elec," ,MS2=0,"
        write(iunit,"(A)",advance='no') "ORBSYM= 1,"
        do i=2,EmbSize-1
            write(iunit,"(A)",advance='no') "1,"
        enddo
        write(iunit,"(A)") "1,"
        write(iunit,"(A)") "ISYM=1"
        write(iunit,"(A)") "&END"

        !Just define diagonal 2 electron contribution over impurity sites
        if(tAnderson) then
            write(iunit,"(F16.12,4I8)") U,1,1,1,1
        else
            do i=1,nImp
                write(iunit,"(F16.12,4I8)") U,i,i,i,i
            enddo
        endif

        !Now for 1electron contributions
        do i=1,EmbSize
            do j=1,i
                if(tChemPot.and.(i.eq.1).and.(j.eq.1)) then
                    write(iunit,"(F16.12,4I8)") Emb_h0v(i,j)-(U/2.0_dp),i,j,0,0
                elseif(abs(Emb_h0v(i,j)).gt.1.0e-10_dp) then
                    write(iunit,"(F16.12,4I8)") Emb_h0v(i,j),i,j,0,0
                endif
            enddo
        enddo
        !Core energy is zero
        write(iunit,"(F16.12,4I8)") 0.0_dp,0,0,0,0
        close(iunit)

    end subroutine WriteFCIDUMP



end module solvers
