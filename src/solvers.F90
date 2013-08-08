module solvers
    use const
    use errors, only: stop_all
    use globals
    use mat_tools, only: WriteVector,WriteMatrix,writematrixcomp
    implicit none

    contains

    !Transform the integrals over the embedded system, and solve for the energy and 1RDM on impurity.
    !tConverged means that the DMET is converged, and this is the last run, so calculate the 2RDM for other properties
    subroutine SolveSystem(tCreate2RDM)
        use utils, only: get_free_unit 
        use DetToolsData, only: FCIDetList,nFCIDet
        implicit none
        logical, intent(in) :: tCreate2RDM
        integer :: pSpaceDim,i,iunit,j,k,l,ios,i_spin,j_spin!,Pivots(4),info
        character(len=256) :: cmd
        character(len=128) :: cmd3
        character(len=73) :: cmd2
        character(len=67) :: cmd1
        character(len=6) :: StrPSpace
        real(dp) :: Emb_nElec,Hel,AFOrder!,DD_Response,ZerothH(4,4)
!        real(dp), allocatable :: FullH1(:,:),LR_State(:)
!        real(dp) ::DDOT,Overlap
        real(dp) :: Check2eEnergy,trace
        logical :: exists
        real(dp), allocatable :: HL_2RDM_temp(:,:),temp2rdm(:,:,:,:)
        real(dp), allocatable :: CoreH(:,:),W(:),work(:),temp(:,:)
        integer :: lWork,info,a,b,c,d
        character(len=*), parameter :: t_r="SolveSystem"

        !Calculate the number of electrons in the embedded system 
        write(6,"(A,I5)") "Number of electrons in full system: ",NEl
        Emb_nElec = 0.0_dp
        do i=1,EmbSize
            Emb_nElec = Emb_nElec + Emb_MF_DM(i,i)
        enddo
        if(tUHF) then
            write(6,"(A,F10.7)") "Number of electrons in alpha spin channel: ",Emb_nElec
            do i=1,EmbSize
                Emb_nElec = Emb_nElec + Emb_MF_DM_b(i,i)
            enddo
        endif
        Elec = nint(Emb_nElec)
        write(6,"(A,F10.7,I5)") "Number of total electrons in embedded system: ",Emb_nElec,Elec

        if(tCompleteDiag.or.tNonDirDavidson) then
            !Do a complete diagonalization, or solve with in-built non-direct davidson diagonalizer
            !Do not need to write FCIDUMP, since would only read it back in...
            call CompleteDiag(tCreate2RDM)
        elseif(tFCIQMC) then
            !Solve with FCIQMC
            call WriteFCIDUMP()

            !Ensure there is a correct input file
            inquire(file='input.neci',exist=exists)
            if(.not.exists) call stop_all(t_r,'Cannot find neci input file: input.neci')

            !Change number of electrons in input file
            cmd1 = "sed -e 's/XXX/'"
            write(StrPSpace,'(I6)') Elec
            do i=1,6
                if(StrPSpace(i:i).ne.' ') exit
            enddo
            if(i.eq.7) call stop_all(t_r,'Error constructing input file call to NECI')
            cmd2 = trim(adjustr(cmd1))//trim(adjustl(StrPSpace(i:6)))
            cmd3 = "'/g' input.neci > input.tmp"
            cmd = trim(adjustl(cmd2))//trim(adjustl(cmd3))
            !write(6,*) cmd
            call system(cmd)
            inquire(file='input.tmp',exist=exists)
            if(.not.exists) call stop_all(t_r,'Intermediate input.tmp file not found')
            !Overwrite initial input file
            call rename('input.tmp','input.neci')

            if(nNECICores.eq.0) then
                !Serial neci run
                cmd = "./neci.x input.neci > neci.out"
                write(6,"(A,A)") "Calling serial fciqmc code with system call: ",cmd
            else
                !Parallel neci run
                cmd1 = "mpirun -np "
                write(StrPSpace,'(I6)') nNECICores
                do i=1,6
                    if(StrPSpace(i:i).ne.' ') exit
                enddo
                if(i.eq.7) call stop_all(t_r,'Error constructing system call to NECI')
                cmd2 = cmd1//trim(adjustl(StrPSpace(i:6)))
                cmd3 = " ./neci.x input.neci > neci.out"
                cmd = adjustl(cmd2)//trim(adjustl(cmd3))
                write(6,"(A,A)") "Calling parallel fciqmc code with system call: ",cmd
            endif
            call system(cmd)

            !TODO: Check here whether FCIQMC calculation was successful or not
            call system("grep '*TOTAL ENERGY*' neci.out | awk '{print$9}' > FCI.ene")

            !TODO: How to calculate the high-level energy? Projected energy? Shift? Density matrix? input option?
            !Assume initially that we get it from the density matrix for consistency
            !Extract energy value to new file and read this is
            HL_Energy = 0.0_dp
            iunit=get_free_unit()
            open(iunit,file='FCI.ene',status='old')
            read(iunit,*) HL_Energy
            close(iunit)
            if(HL_Energy.eq.0.0_dp) call stop_all(t_r,"FCI energy is 0")

            !TODO: Read in 2RDMs, so we can check the energy, and calculate 2-particle stuff
            !read in the 1 (and 2?) RDMs
            if(allocated(HL_1RDM)) deallocate(HL_1RDM)
            allocate(HL_1RDM(EmbSize,EmbSize))
            HL_1RDM(:,:) = 0.0_dp
            iunit=get_free_unit()
            inquire(file='OneRDM',exist=exists)
            if(.not.exists) call stop_all(t_r,'"OneRDM" file not found after NECI calculation')
            open(iunit,file='OneRDM',status='old')
            ios = 0
            do while(ios.eq.0)
                read(iunit,*,iostat=ios) i_spin,j_spin,Hel
                if(ios.eq.0) then
                    if(mod(i_spin,2).eq.0) then
                        if(mod(j_spin,2).eq.0) then
                            HL_1RDM(i_spin/2,j_spin/2) = 2.0_dp*Hel
                            HL_1RDM(j_spin/2,i_spin/2) = 2.0_dp*Hel
                        else
                            HL_1RDM(i_spin/2,(j_spin+1)/2) = 2.0_dp*Hel
                            HL_1RDM((j_spin+1)/2,i_spin/2) = 2.0_dp*Hel
                        endif
                    else
                        if(mod(j_spin,2).eq.0) then
                            HL_1RDM((i_spin+1)/2,j_spin/2) = 2.0_dp*Hel
                            HL_1RDM(j_spin/2,(i_spin+1)/2) = 2.0_dp*Hel
                        else
                            HL_1RDM((i_spin+1)/2,(j_spin+1)/2) = 2.0_dp*Hel
                            HL_1RDM((j_spin+1)/2,(i_spin+1)/2) = 2.0_dp*Hel
                        endif
                    endif
                endif
            enddo
            if(ios.gt.0) call stop_all(t_r,'Error reading 1RDM')
            close(iunit)

            if(tCoreH_EmbBasis) then
                !We have read in the 1RDM in the CoreH basis - rotate back to
                !original basis
                !Diagonalize this original coreH basis
                allocate(CoreH(EmbSize,EmbSize))
                allocate(W(EmbSize))
                CoreH(:,:) = Emb_h0v(:,:)
                if(tChemPot) CoreH(1,1) = CoreH(1,1) - (U/2.0_dp)
                W(:) = 0.0_dp

                !Diagonalize
                allocate(work(1))
                lWork=-1
                info = 0
                call dsyev('V','U',EmbSize,CoreH,EmbSize,W,Work,lWork,info)
                if(info.ne.0) call stop_all(t_r,'Workspace query failed')
                lwork = int(work(1))+1
                deallocate(work)
                allocate(work(lwork))
                call dsyev('V','U',EmbSize,CoreH,EmbSize,W,Work,lWork,info)
                if(info.ne.0) call stop_all(t_r,'Diag failed')
                deallocate(work)

                !Now, rotate the RDM
                allocate(temp(EmbSize,EmbSize))
                call dgemm('n','n',EmbSize,EmbSize,EmbSize,1.0_dp,CoreH,EmbSize,HL_1RDM,EmbSize,0.0_dp,temp,EmbSize)
                call dgemm('n','t',EmbSize,EmbSize,EmbSize,1.0_dp,temp,EmbSize,CoreH,EmbSize,0.0_dp,HL_1RDM,EmbSize)
                deallocate(temp,CoreH,W)
            endif

        else
            !Solve with call to FCI code
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

            if(tCoreH_EmbBasis) then
                !We have read in the 1RDM in the CoreH basis - rotate back to
                !original basis
                !Diagonalize this original coreH basis
                allocate(CoreH(EmbSize,EmbSize))
                allocate(W(EmbSize))
                CoreH(:,:) = Emb_h0v(:,:)
                if(tChemPot) CoreH(1,1) = CoreH(1,1) - (U/2.0_dp)
                W(:) = 0.0_dp

                !Diagonalize
                allocate(work(1))
                lWork=-1
                info = 0
                call dsyev('V','U',EmbSize,CoreH,EmbSize,W,Work,lWork,info)
                if(info.ne.0) call stop_all(t_r,'Workspace query failed')
                lwork = int(work(1))+1
                deallocate(work)
                allocate(work(lwork))
                call dsyev('V','U',EmbSize,CoreH,EmbSize,W,Work,lWork,info)
                if(info.ne.0) call stop_all(t_r,'Diag failed')
                deallocate(work)

                !Now, rotate the RDM
                allocate(temp(EmbSize,EmbSize))
                call dgemm('n','n',EmbSize,EmbSize,EmbSize,1.0_dp,CoreH,EmbSize,HL_1RDM,EmbSize,0.0_dp,temp,EmbSize)
                call dgemm('n','t',EmbSize,EmbSize,EmbSize,1.0_dp,temp,EmbSize,CoreH,EmbSize,0.0_dp,HL_1RDM,EmbSize)
                deallocate(temp)

                if(tCreate2RDM) then
                    !Also transform the 2RDM back to the original basis
                    allocate(temp2rdm(EmbSize,EmbSize,EmbSize,EmbSize))
                    temp2rdm(:,:,:,:) = 0.0_dp

                    do a = 1,EmbSize
                        do b = 1,EmbSize
                            do c = 1,EmbSize
                                do d = 1,EmbSize
                                    hel = 0.0_dp
                                    do i = 1,EmbSize
                                        do j = 1,EmbSize
                                            do k = 1,EmbSize
                                                do l = 1,EmbSize
                                                    hel = hel + CoreH(a,i)*CoreH(b,j)*CoreH(c,k)*CoreH(d,l)*HL_2RDM(i,j,k,l)
                                                enddo
                                            enddo
                                        enddo
                                    enddo
                                    temp2rdm(a,b,c,d) = hel
                                enddo
                            enddo
                        enddo
                    enddo

                    HL_2RDM(:,:,:,:) = temp2rdm(:,:,:,:)
                    deallocate(temp2rdm)

                endif
                deallocate(W,CoreH)
            endif   !tCoreH_EmbBasis
            
            iunit=get_free_unit()
            inquire(file='FCI1RDM',exist=exists)
            if(exists) then
                open(unit=iunit,file='FCI1RDM',status='old')
                close(iunit,status='delete')
            endif
            inquire(file='FCI2RDM',exist=exists)
            if(exists) then
                open(unit=iunit,file='FCI2RDM',status='old')
                close(iunit,status='delete')
            endif
            inquire(file='FCIVec',exist=exists)
            if(exists) then
                open(unit=iunit,file='FCIVec',status='old')
                close(iunit,status='delete')
            endif
            inquire(file='FCI.out',exist=exists)
            if(exists) then
                open(unit=iunit,file='FCI.out',status='old')
                close(iunit,status='delete')
            endif
            inquire(file='FCI.ene',exist=exists)
            if(exists) then
                open(unit=iunit,file='FCI.ene',status='old')
                close(iunit,status='delete')
            endif

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
                    if(tUHF) One_ElecE = One_ElecE + (Emb_h0v_b(i,j)-(U/2.0_dp))*HL_1RDM_b(i,j)
                else
                    One_ElecE = One_ElecE + Emb_h0v(i,j)*HL_1RDM(i,j)
                    if(tUHF) One_ElecE = One_ElecE + Emb_h0v_b(i,j)*HL_1RDM_b(i,j)
                endif
            enddo
        enddo

        !The two electron contribution to the embedded system energy is just the FCI result - the 1e energy
        Two_ElecE = HL_Energy - One_ElecE

        if(tCreate2RDM) then

            if(tUHF) call stop_all(t_r,'2RDM not created with UHF')
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
                    if(tUHF) One_ElecE_Imp = One_ElecE_Imp + (Emb_h0v_b(i,j)-U/2.0_dp)*HL_1RDM_b(i,j)
                else
                    One_ElecE_Imp = One_ElecE_Imp + Emb_h0v(i,j)*HL_1RDM(i,j)
                    if(tUHF) One_ElecE_Imp = One_ElecE_Imp + Emb_h0v_b(i,j)*HL_1RDM_b(i,j)
                endif
            enddo
        enddo
        One_ElecE_Imp = One_ElecE_Imp/real(nImp)    !For energy per impurity

        !Two electron energy terms are not contained in core hamiltonian, or expressed over the bath, 
        !so just divided total 2e contribution by number of impurities to get the 2e contrib
        !Note, we can only do this since there is no correlated bath??
        Two_ElecE_Imp = Two_ElecE/real(nImp)

        !There is also finally an interaction term between the bath and impurity which we can calculate
        !Model this as the h0v matrix with added 1/2 correlation potential so that there is no double counting between the interactions with the bath-impurity
        CoupE_Imp = 0.0_dp
        do j=nImp+1,EmbSize
            do i=1,nImp
                CoupE_Imp = CoupE_Imp + (Emb_h0v(i,j) + 0.5_dp*Emb_CorrPot(i,j))*HL_1RDM(i,j)
                if(tUHF) CoupE_Imp = CoupE_Imp + (Emb_h0v_b(i,j) + 0.5_dp*Emb_CorrPot_b(i,j))*HL_1RDM_b(i,j)
!                CoupE_Imp = CoupE_Imp + (Emb_h0v(i,j) + 0.5_dp*v_loc(i,j))*HL_1RDM(i,j)
            enddo
        enddo
        CoupE_Imp = CoupE_Imp / real(nImp)

        TotalE_Imp = One_ElecE_Imp + Two_ElecE_Imp + CoupE_Imp

        write(6,"(A,F10.6)") "One electron energy per impurity:     ",One_ElecE_Imp
        write(6,"(A,F10.6)") "Two electron energy per impurity:     ",Two_ElecE_Imp
        write(6,"(A,F10.6)") "Coupling energy to bath per impurity: ",CoupE_Imp
        write(6,"(A,2F10.6)") "Total energy per impurity site:       ",TotalE_Imp
        call flush(6)

        !The target filling is the original filling of the mean field RDM per site
        !We know that the filling should be uniform due to the translational symmetry
        Targetfilling_Imp = 0.0_dp
        do i=1,nSites
            Targetfilling_Imp = Targetfilling_Imp + MeanFieldDM(i,i)
            if(tUHF) Targetfilling_Imp = Targetfilling_Imp + MeanFieldDM_b(i,i)
        enddo
        Targetfilling_Imp = Targetfilling_Imp / (2.0_dp*real(nSites))

        !Now what is the actually filling on each impurity from the embedded system
        Actualfilling_Imp = 0.0_dp
        do i=1,nImp
            Actualfilling_Imp = Actualfilling_Imp + HL_1RDM(i,i)
            if(tUHF) Actualfilling_Imp = Actualfilling_Imp + HL_1RDM_b(i,i)
        enddo
        Actualfilling_Imp = Actualfilling_Imp / (2.0_dp*real(nImp))
        Fillingerror = Actualfilling_Imp - Targetfilling_Imp

        write(6,"(A,F15.7)") "Target filling per site: ",Targetfilling_Imp
        write(6,"(A,F15.7)") "Actual filling per site: ",Actualfilling_Imp
        write(6,"(A,F15.7)") "Filling error  per site: ",Fillingerror
        write(6,"(A)") ""
        write(6,"(A)") "Impurity RDMs: "
        call writematrix(HL_1RDM(1:nImp,1:nImp),'hl_1RDM_Imp',.true.)
        if(tUHF) call writematrix(HL_1RDM_b(1:nImp,1:nImp),'hl_1RDM_Imp_b',.true.)

        if(tUHF) then
            call writematrix(HL_1RDM(1:nImp,1:nImp)-HL_1RDM_b(1:nImp,1:nImp),'Spin Density',.true.)
            AFOrder = zero
            do i = 1,nImp
                AFOrder = AFOrder + abs(HL_1RDM(i,i)-HL_1RDM_b(i,i))
            enddo
            AFOrder = AFOrder/nImp
            write(6,"(A,F15.7)") "AFM order: ",AFOrder
        endif
        call flush(6)

!        if(tDebug) call writematrix(HL_1RDM,'hl_1rdm',.true.)
        if(tWriteOut) then
            call writematrix(HL_1RDM,'hl_1rdm',.true.)
            if(tUHF) call writematrix(HL_1RDM_b,'hl_1rdm_b',.true.)
        endif

    end subroutine SolveSystem

    !Calculate the FCI result for this model in the full space
    subroutine DiagFullSystem()
        use utils, only: get_free_unit,append_ext,append_ext_real
        use DetBitOps, only: SQOperator
        use DetToolsData
        use DetTools, only: umatind,GenDets,GetHElement
        implicit none
        integer :: pSpaceDim,i,j,iunit,pertsitealpha,pertsitebeta,UMatSize,OrbPairs
        integer :: lWork,info,ilut,pertsite
        logical :: tParity
        complex(dp) :: DDRes,GFRes,GFRes_h,GFRes_p
        real(dp) :: ddot,Overlap,Omega,mu
        real(dp), allocatable :: V0(:),Work(:),V0_Ann(:),V0_Cre(:),Spectrum_Np1(:)
        real(dp), allocatable :: Spectrum_Nm1b(:),FullHamil_Np1(:,:),FullHamil_Nm1b(:,:)
        character(len=64) :: filename,filename2
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
                    if(tChemPot.and.(i.eq.1).and.(j.eq.1)) then
                        write(iunit,"(F16.12,4I8)") h0(1,1)-(U/2.0_dp),1,1,0,0
                    elseif(abs(h0(i,j)).gt.1.0e-10_dp) then
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

            write(6,*) "Performing complete diagonalization of full system...SLOW!"
            !First, allocate and fill the umat and tmat for the FCI space
            OrbPairs = (nSites*(nSites+1))/2
            UMatSize = (OrbPairs*(OrbPairs+1))/2
            write(6,*) "Allocating memory to store 2 electron integrals: ",UMatSize
            if(allocated(UMat)) deallocate(UMat)
            allocate(UMat(UMatSize))
            UMat(:) = 0.0_dp
            if(tAnderson) then
                umat(umatind(1,1,1,1)) = U
            else
                do i=1,nSites
                    umat(umatind(i,i,i,i)) = U
                enddo
            endif
            if(allocated(tmat)) deallocate(tmat)
            allocate(tmat(nSites,nSites))
            tmat(:,:) = 0.0_dp
            do i=1,nSites
                do j=1,nSites
                    if(abs(h0(i,j)).gt.1.0e-10_dp) then
                        tmat(i,j) = h0(i,j)
                    endif
                enddo
            enddo
            if(tChemPot) then
                tmat(1,1) = tmat(1,1) - U/2.0_dp
            endif

            !Now generate all determinants in the active space
            if(allocated(FCIDetList)) deallocate(FCIDetList)
            call GenDets(NEl,nSites,.true.,.true.,.true.)
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
                    call GetHElement(FCIDetList(:,i),FCIDetList(:,j),NEl,FullHamil(i,j))
                enddo
            enddo

            write(6,*) "Diagonalizing full system..."
            call flush(6)
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

            write(6,*) "Diagonalization complete..."
            write(6,*) "True ground state energy: ",Spectrum(1)
            call flush(6)

            if(tDDResponse) then
                write(6,*) "Calculating the exact *LOCAL* density-density response function..."

                iunit = get_free_unit()
                call append_ext_real('Exact_DDResponse',U,filename)
                if(.not.tHalfFill) then
                    call append_ext(filename,nOcc,filename2)
                else
                    filename2 = filename
                endif
                open(unit=iunit,file=filename2,status='unknown')
                write(iunit,"(A)") "# Frequency   DD_Response(Re)   DD_Response(Im)   "

                !Find V|0>
                allocate(V0(nFCIDet))
                V0(:) = 0.0_dp

                pertsite = 1    !The AO site that the perturbation takes place at
        
                pertsitealpha = 2*pertsite-1
                pertsitebeta = 2*pertsite

                do i = 1,nFCIDet
                    if(btest(FCIBitList(i),pertsitealpha-1)) then
                        V0(i) = V0(i) + FullHamil(i,1)
                    endif
                    if(btest(FCIBitList(i),pertsitebeta-1)) then
                        V0(i) = V0(i) + FullHamil(i,1)
                    endif
                enddo

                Omega = Start_Omega
                do while((Omega.lt.max(Start_Omega,End_Omega)+1.0e-5_dp).and.(Omega.gt.min(Start_Omega,End_Omega)-1.0e-5_dp))

                    DDRes = dcmplx(0.0_dp,0.0_dp)
                    do i = 2,nFCIDet

                        Overlap = ddot(nFCIDet,V0,1,FullHamil(:,i),1)
                        Overlap = Overlap**2
                        DDRes = DDRes + dcmplx(-Overlap,0.0_dp)/dcmplx(Spectrum(i)-Spectrum(1)-Omega,-dDelta)
                    enddo

                    write(iunit,"(3G22.10)") Omega,real(DDRes),-aimag(DDRes)

                    Omega = Omega + Omega_Step
                enddo

                close(iunit)
                deallocate(V0)

            endif
            if(tChargedResponse) then
                write(6,*) "Calculating the exact *LOCAL* greens function..."
                
                iunit = get_free_unit()
                call append_ext_real('Exact_GFResponse',U,filename)
                if(.not.tHalfFill) then
                    call append_ext(filename,nOcc,filename2)
                else
                    filename2 = filename
                endif
                open(unit=iunit,file=filename2,status='unknown')
                write(iunit,"(A)") "# Frequency   GF_Response(Re)   GF_Response(Im)   GF_part(Re)  GF_Part(Im)" &
                    //" GF_hole(Re)    GF_hole(Im)    "

                !Find V|0>
                allocate(V0_Cre(nNp1FCIDet))
                allocate(V0_Ann(nNm1bFCIDet))
                V0_Cre(:) = 0.0_dp
                V0_Ann(:) = 0.0_dp

                pertsite = 1
        
                pertsitealpha = 2*pertsite-1
                do i = 1,nFCIDet
                    if(.not.btest(FCIBitList(i),pertsitealpha-1)) then
                        ilut = FCIBitList(i)
                        call SQOperator(ilut,pertsitealpha,tParity,.false.)
                        do j = 1,nNp1FCIDet
                            if(Np1BitList(j).eq.ilut) then
                                if(tParity) then
                                    V0_Cre(j) = -FullHamil(i,1)
                                else
                                    V0_Cre(j) = FullHamil(i,1)
                                endif
                                exit
                            endif
                        enddo
                        if(j.ge.(nNp1FCIDet+1)) call stop_all(t_r,'Could not find corresponding determinant')
                    else
                        ilut = FCIBitList(i)
                        call SQOperator(ilut,pertsitealpha,tParity,.true.)
                        do j = 1,nNm1bFCIDet
                            if(Nm1bBitList(j).eq.ilut) then
                                if(tParity) then
                                    V0_Ann(j) = -FullHamil(i,1)
                                else
                                    V0_Ann(j) = FullHamil(i,1)
                                endif
                                exit
                            endif
                        enddo
                        if(j.gt.(nNm1bFCIDet+1)) call stop_all(t_r,'Could not find corresponding determinant')
                    endif
                enddo

                !Now find full eigenbasis for N+1_alpha and N-1_beta space
                allocate(Spectrum_Np1(nNp1FCIDet))
                allocate(FullHamil_Np1(nNp1FCIDet,nNp1FCIDet))
                FullHamil_Np1(:,:) = 0.0_dp
                !Construct the hamiltonian - slow
                do i=1,nNp1FCIDet
                    do j=1,nNp1FCIDet
                        call GetHElement(Np1FCIDetList(:,i),Np1FCIDetList(:,j),NEl+1,FullHamil_Np1(i,j))
                    enddo
                enddo

                write(6,*) "Diagonalizing N+1 system..."
                call flush(6)
                !Diagonalize
                allocate(Work(1))
                lWork=-1
                info=0
                call dsyev('V','U',nNp1FCIDet,FullHamil_Np1,nNp1FCIDet,Spectrum_Np1,Work,lWork,info)
                if(info.ne.0) call stop_all(t_r,'Workspace quiery failed')
                lwork=int(work(1))+1
                deallocate(work)
                allocate(work(lwork))
                call dsyev('V','U',nNp1FCIDet,FullHamil_Np1,nNp1FCIDet,Spectrum_Np1,Work,lWork,info)
                if(info.ne.0) call stop_all(t_r,'Diag failed')
                deallocate(work)

                allocate(Spectrum_Nm1b(nNm1bFCIDet))
                allocate(FullHamil_Nm1b(nNm1bFCIDet,nNm1bFCIDet))
                FullHamil_Nm1b(:,:) = 0.0_dp
                !Construct the hamiltonian - slow
                do i=1,nNm1bFCIDet
                    do j=1,nNm1bFCIDet
                        call GetHElement(Nm1bFCIDetList(:,i),Nm1bFCIDetList(:,j),NEl-1,FullHamil_Nm1b(i,j))
                    enddo
                enddo

                write(6,*) "Diagonalizing N-1 system..."
                call flush(6)
                !Diagonalize
                allocate(Work(1))
                lWork=-1
                info=0
                call dsyev('V','U',nNm1bFCIDet,FullHamil_Nm1b,nNm1bFCIDet,Spectrum_Nm1b,Work,lWork,info)
                if(info.ne.0) call stop_all(t_r,'Workspace quiery failed')
                lwork=int(work(1))+1
                deallocate(work)
                allocate(work(lwork))
                call dsyev('V','U',nNm1bFCIDet,FullHamil_Nm1b,nNm1bFCIDet,Spectrum_Nm1b,Work,lWork,info)
                if(info.ne.0) call stop_all(t_r,'Diag failed')
                deallocate(work)

                write(6,*) "Constructing exact greens function..."
                call flush(6)

                Omega = Start_Omega
                do while((Omega.lt.max(Start_Omega,End_Omega)+1.0e-5_dp).and.(Omega.gt.min(Start_Omega,End_Omega)-1.0e-5_dp))
            
                    if(.not.tAnderson) then
                        !In the hubbard model, apply a chemical potential of U/2
                        mu = U/2.0_dp
                    else
                        mu = 0.0_dp
                    endif

                    GFRes_p = dcmplx(0.0_dp,0.0_dp)
                    GFRes_h = dcmplx(0.0_dp,0.0_dp)
                    do i = 1,nNp1FCIDet

                        Overlap = ddot(nNp1FCIDet,V0_Cre,1,FullHamil_Np1(:,i),1)
                        Overlap = Overlap**2
                        GFRes_p = GFRes_p + dcmplx(Overlap,0.0_dp)/dcmplx(Omega+mu-Spectrum_Np1(i)+Spectrum(1),dDelta)
                    enddo
                    do i = 1,nNm1bFCIDet

                        Overlap = ddot(nNm1bFCIDet,V0_Ann,1,FullHamil_Nm1b(:,i),1)
                        Overlap = Overlap**2
                        GFRes_h = GFRes_h + dcmplx(Overlap,0.0_dp)/dcmplx(Omega+mu+Spectrum_Nm1b(i)-Spectrum(1),dDelta)
                    enddo

                    GFRes = GFRes_p + GFRes_h

                    write(iunit,"(7G22.10)") Omega,real(GFRes),-aimag(GFRes),real(GFRes_p),-aimag(GFRes_p), &
                        real(GFRes_h),-aimag(GFRes_h)

                    Omega = Omega + Omega_Step
                enddo

                close(iunit)
                deallocate(V0_Cre,V0_Ann,FullHamil_Nm1b,FullHamil_Np1,Spectrum_Nm1b,Spectrum_Np1)

            endif

            deallocate(FCIBitList,Nm1FCIDetList,Nm1BitList,FullHamil,Spectrum,FCIDetList,umat,tmat)
            deallocate(Np1FCIDetList,Np1BitList,Nm1bFCIDetList,Nm1bBitList,Np1bFCIDetList,Np1bBitList)

        endif

    end subroutine DiagFullSystem

!Generate all determinants in the FCI space and do complete diagonalization
    subroutine CompleteDiag(tCreate2RDM)
        use DetToolsData
        use Davidson, only: Real_NonDir_Davidson
        use DetTools, only: GetHElement,GenDets
        use utils, only: get_free_unit
        implicit none
        logical, intent(in) :: tCreate2RDM
        integer :: Nmax,ierr
        real(dp), allocatable :: work(:),CompressHam(:)
        integer, allocatable :: IndexHam(:)
        integer :: lwork,info,i,j,iSize,iunit_tmp
        logical :: texist
        character(len=*), parameter :: t_r='CompleteDiag'

        call CreateIntMats()

        !Now generate all determinants in the active space
        if(allocated(FCIDetList)) deallocate(FCIDetList)
        if(allocated(FCIBitList)) deallocate(FCIBitList)
        call GenDets(Elec,EmbSize,.false.,.true.,.false.)
        !FCIDetList now stores a list of all the determinants
        write(6,"(A,I14)") "Number of determinants in FCI space: ",nFCIDet
        if(tCompressedMats) then
            write(6,"(A)") "Compressing hamiltonian matrix..."
            if(tReadMats) then
                inquire(file='CompressHam_N',exist=texist)
                if(.not.texist) then
                    call CountSizeCompMat(FCIDetList(:,:),Elec,nFCIDet,Nmax,FCIBitList(:))
                else
                    iunit_tmp = get_free_unit()
                    open(iunit_tmp,file='CompressHam_N',status='old')
                    read(iunit_tmp,*) isize
                    Nmax = isize
                endif
            else
                if(iHamSize_N.eq.0) then
                    call CountSizeCompMat(FCIDetList(:,:),Elec,nFCIDet,Nmax,FCIBitList(:))
                    !iHamSize_N = Nmax  Cannot strictly do this.
                else
                    Nmax = iHamSize_N
                endif
            endif
            write(6,"(A,2I14)") "Size of Compressed/Full Hamiltonians: ",Nmax,nFCIDet**2
            write(6,"(A,F14.6,A)") "Allocating memory for compressed hamiltonian: ",real(Nmax*16,dp)/1048576.0_dp," Mb"
            write(6,"(A,F14.6,A)") "Saving in compression of: ",(real(nFCIDet**2 - Nmax*2,dp))*8/1048576.0_dp," Mb"
            allocate(CompressHam(Nmax),stat=ierr)
            allocate(IndexHam(Nmax),stat=ierr)
            if(ierr.ne.0) call stop_all(t_r,'Allocation failed')
            if(tReadMats.and.texist) then
                if(isize.ne.Nmax) call stop_all(t_r,'Error here')
                do i = 1,iSize
                    read(iunit_tmp,*) CompressHam(i),IndexHam(i)
                enddo
                close(iunit_tmp)
            else
                call StoreCompMat(FCIDetList(:,:),Elec,nFCIDet,Nmax,CompressHam,IndexHam,FCIBitList(:))
            endif
        else
            write(6,"(A,F14.6,A)") "Allocating memory for the hamiltonian: ",real((nFCIDet**2)*8,dp)/1048576.0_dp," Mb"
            if(allocated(FullHamil)) deallocate(FullHamil)
            allocate(FullHamil(nFCIDet,nFCIDet))
            FullHamil(:,:) = 0.0_dp
            !Construct the hamiltonian - slow
            do i=1,nFCIDet
                do j=1,nFCIDet
                    call GetHElement(FCIDetList(:,i),FCIDetList(:,j),Elec,FullHamil(i,j))
                    !write(6,*) i,j,FullHamil(i,j)
                enddo
            enddo
        endif
!        call writematrix(FullHamil(1:nFCIDet,1:nFCIDet),'FCI hamil',.true.)

        !Diagonalize
        if(allocated(HL_Vec)) deallocate(HL_Vec)
        allocate(HL_Vec(nFCIDet))
        if(tNonDirDavidson) then
            write(6,*) "Solving for ground state with non-direct davidson diagonalizer..."
            if(tCompressedMats) then
                call Real_NonDir_Davidson(nFCIDet,HL_Energy,HL_Vec,.false.,Nmax,CompressMat=CompressHam,IndexMat=IndexHam)
                if(tWriteMats) then
                    if(.not.(tReadMats.and.texist)) then
                        !Do not write out if we have already got them!
                        write(6,*) "Writing out matrices"
                        iunit_tmp = get_free_unit()
                        open(iunit_tmp,file='CompressHam_N',status='unknown')
                        write(iunit_tmp,*) Nmax
                        do i = 1,Nmax
                            write(iunit_tmp,*) CompressHam(i),IndexHam(i)
                        enddo
                        close(iunit_tmp)
                    endif
                endif
                deallocate(CompressHam,IndexHam)
            else
                call Real_NonDir_Davidson(nFCIDet,HL_Energy,HL_Vec,.false.,1,Mat=FullHamil)
                deallocate(FullHamil)
            endif
        else
            if(allocated(Spectrum)) deallocate(Spectrum)
            allocate(Spectrum(nFCIDet))
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

            !call writevector(Spectrum(1:(min(10,nFCIDet))),'FCI Spectrum')
            HL_Energy = Spectrum(1) !GS eigenvalue
            HL_Vec(:) = FullHamil(:,1)
            deallocate(Spectrum,FullHamil)
        endif

        write(6,*) "FCI energy: ",HL_Energy
        call flush(6)
            
        if(allocated(HL_1RDM)) deallocate(HL_1RDM)
        allocate(HL_1RDM(EmbSize,EmbSize))
        HL_1RDM(:,:) = zero
        if(tUHF) then
            if(allocated(HL_1RDM_b)) deallocate(HL_1RDM_b)
            allocate(HL_1RDM_b(EmbSize,EmbSize))
            HL_1RDM_b(:,:) = zero
        endif

        if(tReadMats) then
            inquire(file='RDM_N',exist=texist)
            if(texist) then
                iunit_tmp = get_free_unit()
                open(iunit_tmp,file='RDM_N',status='old')
                do i = 1,EmbSize
                    do j = 1,EmbSize
                        read(iunit_tmp,*) HL_1RDM(i,j)
                    enddo
                enddo
                close(iunit_tmp)
                if(tUHF) then
                    iunit_tmp = get_free_unit()
                    open(iunit_tmp,file='RDM_N_b',status='old')
                    do i = 1,EmbSize
                        do j = 1,EmbSize
                            read(iunit_tmp,*) HL_1RDM_b(i,j)
                        enddo
                    enddo
                    close(iunit_tmp)
                endif
            else
                if(tUHF) then
                    call FindFull1RDM(1,1,.true.,HL_1RDM,RDM_Beta=HL_1RDM_b)
                else
                    call FindFull1RDM(1,1,.true.,HL_1RDM)
                endif
            endif
        else
            if(tUHF) then
                call FindFull1RDM(1,1,.true.,HL_1RDM,RDM_Beta=HL_1RDM_b)
            else
                call FindFull1RDM(1,1,.true.,HL_1RDM)
            endif
        endif
        if(tWriteMats) then
            if(.not.(tReadMats.and.texist)) then
                !Do not write them out again if we have just read them in!
                write(6,*) "Writing out RDM"
                iunit_tmp = get_free_unit()
                open(iunit_tmp,file='RDM_N',status='unknown')
                do i = 1,EmbSize
                    do j = 1,EmbSize
                        write(iunit_tmp,*) HL_1RDM(i,j)
                    enddo
                enddo
                close(iunit_tmp)
                if(tUHF) then
                    open(iunit_tmp,file='RDM_N_b',status='unknown')
                    do i = 1,EmbSize
                        do j = 1,EmbSize
                            write(iunit_tmp,*) HL_1RDM_b(i,j)
                        enddo
                    enddo
                    close(iunit_tmp)
                endif
            endif
        endif

        if(tWriteOut) then
            call writematrix(HL_1RDM,'HL_1RDM',.true.)
            if(tUHF) call writematrix(HL_1RDM_b,'Beta HL_1RDM',.true.)
        endif

        if(tCreate2RDM) then
            if(allocated(HL_2RDM)) deallocate(HL_2RDM)
            allocate(HL_2RDM(EmbSize,EmbSize,EmbSize,EmbSize))
            HL_2RDM(:,:,:,:) = zero
            if(tUHF) then
                call stop_all(t_r,'Cannot calculate 2RDM with UHF currently')
            else
                call FindFull2RDM(1,1,.true.,HL_2RDM)
            endif
        endif

        if(allocated(FCIBitList)) deallocate(FCIBitList)

    end subroutine CompleteDiag


    !This subroutine will fill the integral arrays, tmat and umat,
    !used by the sltcnd routines, from the contents of Emb_h0v, and
    !put U over the impurity sites.
    !If tComp = .true., then the 1-electron array tmat_comp is filled, and the elements complex
    !Umat is real either way
    subroutine CreateIntMats(tComp)
        use DetToolsData
        use DetTools, only: umatind
        implicit none
        logical, intent(in), optional :: tComp
        logical :: tComp_
        integer :: OrbPairs,UMatSize,j,i

        if(present(tComp)) then
            tComp_ = tComp
        else
            !Assume that we want real 1-electron matrices if not specified otherwise.
            tComp_ = .false.
        endif

        !First, allocate and fill the umat and tmat for the FCI space
        OrbPairs = (EmbSizeSpin*(EmbSizeSpin+1))/2
        UMatSize = (OrbPairs*(OrbPairs+1))/2
        write(6,*) "Allocating memory to store 2 electron integrals: ",UMatSize
        if(allocated(UMat)) deallocate(UMat)
        allocate(UMat(UMatSize))
        UMat(:) = zero
        if(tAnderson) then
            umat(umatind(1,1,1,1)) = U
            if(tUHF) then
                !mixed spin
                umat(umatind(2,2,2,2)) = U
                umat(umatind(1,2,1,2)) = U
                umat(umatind(2,1,2,1)) = U
            endif
        elseif(tSingFiss) then
            !Store coulomb integrals
            do i = 1,nImp
                do j = 1,nImp
                    umat(umatind(i,j,i,j)) = J_Ints(i,j)
                enddo
            enddo
            !Store exchange integrals
            do i = 1,nImp
                do j = 1,nImp
                    if(i.eq.j) cycle    !The diagonal terms are included in coulomb part
                    umat(umatind(i,j,j,i)) = X_Ints(i,j)
                enddo
            enddo
        else
            do i=1,nImp
                umat(umatind(i,i,i,i)) = U
            enddo
            if(tUHF) then
                do i = nImp+1,nImp*2
                    umat(umatind(i,i,i,i)) = U
                enddo
                do i = 1,nImp*2,2
                    !abab and baba spins
                    umat(umatind(i,i+1,i,i+1)) = U
                    umat(umatind(i+1,i,i+1,i)) = U
                enddo
            endif
        endif

        !EmbSizeSpin = EmbSize for RHF, or EmbSize*2 for UHF
        if(tComp_) then
            if(allocated(tmat_comp)) deallocate(tmat_comp)
            allocate(tmat_comp(EmbSizeSpin,EmbSizeSpin))
            tmat_comp(:,:) = zzero
        else
            if(allocated(tmat)) deallocate(tmat)
            allocate(tmat(EmbSizeSpin,EmbSizeSpin))
            tmat(:,:) = zero
        endif
        do i=1,EmbSizeSpin
            do j=1,EmbSizeSpin
                if(tUHF) then
                    if(mod(i,2).ne.mod(j,2)) cycle
                    if(mod(i,2).eq.1) then
                        !alpha orbital
                        if(tComp_) then
                            tmat_comp(i,j) = dcmplx(Emb_h0v((i-1)/2 + 1,(j-1)/2 + 1),0.0_dp)
                        else
                            tmat(i,j) = Emb_h0v((i-1)/2 + 1,(j-1)/2 + 1)
                        endif
                    else
                        if(tComp_) then
                            tmat_comp(i,j) = dcmplx(Emb_h0v_b(i/2,j/2),0.0_dp)
                        else
                            tmat(i,j) = Emb_h0v_b(i/2,j/2)
                        endif
                    endif
                else
                    if(tComp_) then
                        tmat_comp(i,j) = dcmplx(Emb_h0v(i,j),0.0_dp)
                    else
                        tmat(i,j) = Emb_h0v(i,j)
                    endif
                endif
            enddo
        enddo
        if(tChemPot) then
            if(tComp_) then
                tmat_comp(1,1) = tmat_comp(1,1) - dcmplx(U/2.0_dp,0.0_dp)
                if(tUHF) tmat_comp(2,2) = tmat_comp(2,2) - dcmplx(U/2.0_dp)
            else
                tmat(1,1) = tmat(1,1) - U/2.0_dp
                if(tUHF) tmat(2,2) = tmat(2,2) - U/2.0_dp
            endif
        endif
        if(tWriteOut) then
            if(tComp_) then
                call writematrixcomp(tmat_comp,'tmat_comp',.true.)
            else
                call writematrix(tmat,'tmat',.true.)
            endif
        endif
    end subroutine CreateIntMats

    subroutine FindFull1RDM(StateBra,StateKet,tGroundState,RDM,RDM_Beta)
        use DetToolsData, only: FCIDetList,nFCIDet,FCIBitList
        use DetTools, only: iGetExcitLevel,GetExcitation,tospat
        use DetBitOps, only: FindBitExcitLevel
        implicit none
        real(dp) , intent(out) :: RDM(EmbSize,EmbSize)
        real(dp), intent(out), optional :: RDM_Beta(EmbSize,EmbSize)
        integer , intent(in) :: StateBra,StateKet
        logical , intent(in) :: tGroundState
        real(dp), pointer :: Bra(:),Ket(:)
        integer :: Ex(2),i,j,k,IC
        logical :: tSign
        character(len=*), parameter :: t_r='FindFull1RDM'

        if(tUHF.and.(.not.present(RDM_Beta))) then
            call stop_all(t_r,'Beta spin channel RDM not present')
        endif

        RDM(:,:) = zero 
        if(tUHF) RDM_Beta(:,:) = zero

        if(tGroundState) then
            Bra => HL_Vec(1:nFCIDet)
            Ket => HL_Vec(1:nFCIDet)
        else
            Bra => FullHamil(1:nFCIDet,StateBra)
            Ket => FullHamil(1:nFCIDet,StateKet)
        endif

        if(allocated(FCIBitList)) then
            do i=1,nFCIDet
                do j=i,nFCIDet
                    IC = FindBitExcitLevel(FCIBitList(i),FCIBitList(j))
                    if(IC.eq.1) then
                        !Connected by a single
                        Ex(1) = 1
                        call GetExcitation(FCIDetList(:,i),FCIDetList(:,j),Elec,Ex,tSign)
                        if(tSign) then
                            if(tUHF) then
                                if(mod(Ex(1),2).eq.1) then
                                    RDM(tospat(Ex(1)),tospat(Ex(2))) = RDM(tospat(Ex(1)),tospat(Ex(2))) -   &
                                        Bra(i)*Ket(j)
                                else
                                    RDM_Beta(tospat(Ex(1)),tospat(Ex(2))) = RDM_Beta(tospat(Ex(1)),tospat(Ex(2))) -   &
                                        Bra(i)*Ket(j)
                                endif
                            else
                                RDM(tospat(Ex(1)),tospat(Ex(2))) = RDM(tospat(Ex(1)),tospat(Ex(2))) -   &
                                    Bra(i)*Ket(j)
                            endif
                        else
                            if(tUHF) then
                                if(mod(Ex(1),2).eq.1) then
                                    RDM(tospat(Ex(1)),tospat(Ex(2))) = RDM(tospat(Ex(1)),tospat(Ex(2))) +   &
                                        Bra(i)*Ket(j)
                                else
                                    RDM_Beta(tospat(Ex(1)),tospat(Ex(2))) = RDM_Beta(tospat(Ex(1)),tospat(Ex(2))) +   &
                                        Bra(i)*Ket(j)
                                endif
                            else
                                RDM(tospat(Ex(1)),tospat(Ex(2))) = RDM(tospat(Ex(1)),tospat(Ex(2))) +   &
                                    Bra(i)*Ket(j)
                            endif
                        endif
                        !And also do other way round
                        Ex(1) = 1
                        call GetExcitation(FCIDetList(:,j),FCIDetList(:,i),Elec,Ex,tSign)
                        if(tSign) then
                            if(tUHF) then
                                if(mod(Ex(1),2).eq.1) then
                                    RDM(tospat(Ex(1)),tospat(Ex(2))) = RDM(tospat(Ex(1)),tospat(Ex(2))) -   &
                                        Bra(j)*Ket(i)
                                else
                                    RDM_Beta(tospat(Ex(1)),tospat(Ex(2))) = RDM_Beta(tospat(Ex(1)),tospat(Ex(2))) -   &
                                        Bra(j)*Ket(i)
                                endif
                            else
                                RDM(tospat(Ex(1)),tospat(Ex(2))) = RDM(tospat(Ex(1)),tospat(Ex(2))) -   &
                                    Bra(j)*Ket(i)
                            endif
                        else
                            if(tUHF) then
                                if(mod(Ex(1),2).eq.1) then
                                    RDM(tospat(Ex(1)),tospat(Ex(2))) = RDM(tospat(Ex(1)),tospat(Ex(2))) +   &
                                        Bra(j)*Ket(i)
                                else
                                    RDM_Beta(tospat(Ex(1)),tospat(Ex(2))) = RDM_Beta(tospat(Ex(1)),tospat(Ex(2))) +   &
                                        Bra(j)*Ket(i)
                                endif
                            else
                                RDM(tospat(Ex(1)),tospat(Ex(2))) = RDM(tospat(Ex(1)),tospat(Ex(2))) +   &
                                    Bra(j)*Ket(i)
                            endif
                        endif
                    elseif(IC.eq.0) then
                        !Same det
                        if(i.ne.j) call stop_all(t_r,'Error here')
                        do k=1,Elec
                            if(tUHF) then
                                if(mod(FCIDetList(k,i),2).eq.1) then
                                    RDM(tospat(FCIDetList(k,i)),tospat(FCIDetList(k,i))) =  &
                                        RDM(tospat(FCIDetList(k,i)),tospat(FCIDetList(k,i))) &
                                        + Bra(i)*Ket(j)
                                else
                                    RDM_Beta(tospat(FCIDetList(k,i)),tospat(FCIDetList(k,i))) =     &
                                        RDM_Beta(tospat(FCIDetList(k,i)),tospat(FCIDetList(k,i))) &
                                        + Bra(i)*Ket(j)
                                endif
                            else
                                RDM(tospat(FCIDetList(k,i)),tospat(FCIDetList(k,i))) =  &
                                    RDM(tospat(FCIDetList(k,i)),tospat(FCIDetList(k,i))) &
                                    + Bra(i)*Ket(j)
                            endif
                        enddo
                    endif

                enddo
                if(mod(i,25000).eq.0) then
                    write(6,"(A,2I12)") "Finding 1RDM (stupidly...) ",i,nFCIDet
                    call flush(6)
                endif
            enddo
        else
            do i=1,nFCIDet
                do j=1,nFCIDet
                    IC = iGetExcitLevel(FCIDetList(:,i),FCIDetList(:,j),Elec)
                    if(IC.eq.1) then
                        !Connected by a single
                        Ex(1) = 1
                        call GetExcitation(FCIDetList(:,i),FCIDetList(:,j),Elec,Ex,tSign)
                        if(tSign) then
                            if(tUHF) then
                                if(mod(Ex(1),2).eq.1) then
                                    RDM(tospat(Ex(1)),tospat(Ex(2))) = RDM(tospat(Ex(1)),tospat(Ex(2))) -   &
                                        Bra(i)*Ket(j)
                                else
                                    RDM_Beta(tospat(Ex(1)),tospat(Ex(2))) = RDM_Beta(tospat(Ex(1)),tospat(Ex(2))) -   &
                                        Bra(i)*Ket(j)
                                endif
                            else
                                RDM(tospat(Ex(1)),tospat(Ex(2))) = RDM(tospat(Ex(1)),tospat(Ex(2))) -   &
                                    Bra(i)*Ket(j)
                            endif
                        else
                            if(tUHF) then
                                if(mod(Ex(1),2).eq.1) then
                                    RDM(tospat(Ex(1)),tospat(Ex(2))) = RDM(tospat(Ex(1)),tospat(Ex(2))) +   &
                                        Bra(i)*Ket(j)
                                else
                                    RDM_Beta(tospat(Ex(1)),tospat(Ex(2))) = RDM_Beta(tospat(Ex(1)),tospat(Ex(2))) +   &
                                        Bra(i)*Ket(j)
                                endif
                            else
                                RDM(tospat(Ex(1)),tospat(Ex(2))) = RDM(tospat(Ex(1)),tospat(Ex(2))) +   &
                                    Bra(i)*Ket(j)
                            endif
                        endif
                    elseif(IC.eq.0) then
                        !Same det
                        if(i.ne.j) call stop_all(t_r,'Error here')
                        do k=1,Elec
                            if(tUHF) then
                                if(mod(FCIDetList(k,i),2).eq.1) then
                                    RDM(tospat(FCIDetList(k,i)),tospat(FCIDetList(k,i))) =  &
                                        RDM(tospat(FCIDetList(k,i)),tospat(FCIDetList(k,i))) &
                                        + Bra(i)*Ket(j)
                                else
                                    RDM_Beta(tospat(FCIDetList(k,i)),tospat(FCIDetList(k,i))) =     &
                                        RDM_Beta(tospat(FCIDetList(k,i)),tospat(FCIDetList(k,i))) &
                                        + Bra(i)*Ket(j)
                                endif
                            else
                                RDM(tospat(FCIDetList(k,i)),tospat(FCIDetList(k,i))) =  &
                                    RDM(tospat(FCIDetList(k,i)),tospat(FCIDetList(k,i))) &
                                    + Bra(i)*Ket(j)
                            endif
                        enddo
                    endif

                enddo
            enddo
        endif

        do i=1,EmbSize
            do j=1,EmbSize
                if(abs(RDM(i,j)-RDM(j,i)).gt.1.0e-7_dp) call stop_all(t_r,'1RDM not symmetric')
            enddo
        enddo
        if(tUHF) then
            do i=1,EmbSize
                do j=1,EmbSize
                    if(abs(RDM_Beta(i,j)-RDM_Beta(j,i)).gt.1.0e-7_dp) call stop_all(t_r,'1RDM beta space not symmetric')
                enddo
            enddo
        endif

    end subroutine FindFull1RDM

    !Find spin-integrated 2RDM very inefficiently
    !According to Helgakker <0|e_pqrs|0>
    !Done by running through all N^2 determinant pairs
    subroutine FindFull2RDM(StateBra,StateKet,tGroundState,RDM)
        use DetToolsData, only: FCIDetList,nFCIDet
        use DetTools, only: iGetExcitLevel,GetExcitation,gtid
        implicit none
        real(dp) , intent(out) :: RDM(EmbSize,EmbSize,EmbSize,EmbSize)
        integer , intent(in) :: StateBra,StateKet
        logical , intent(in) :: tGroundState
        real(dp), pointer :: Bra(:),Ket(:)
        integer :: Ex(2,2),i,j,k,IC,kel,lel,l,temp
        logical :: tSign
        character(len=*), parameter :: t_r='FindFull2RDM'
        
        RDM(:,:,:,:) = zero

        if(tGroundState) then
            Bra => HL_Vec(1:nFCIDet)
            Ket => HL_Vec(1:nFCIDet)
        else
            Bra => FullHamil(1:nFCIDet,StateBra)
            Ket => FullHamil(1:nFCIDet,StateKet)
        endif

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
                                Bra(i)*Ket(j)
                            RDM(gtid(Ex(2,2)),gtid(Ex(1,1)),gtid(Ex(2,1)),gtid(Ex(1,2))) =  &
                                RDM(gtid(Ex(2,2)),gtid(Ex(1,1)),gtid(Ex(2,1)),gtid(Ex(1,2))) +  &
                                Bra(i)*Ket(j)
                            RDM(gtid(Ex(2,2)),gtid(Ex(1,2)),gtid(Ex(2,1)),gtid(Ex(1,1))) =  &
                                RDM(gtid(Ex(2,2)),gtid(Ex(1,2)),gtid(Ex(2,1)),gtid(Ex(1,1))) -  &
                                Bra(i)*Ket(j)
                            RDM(gtid(Ex(2,1)),gtid(Ex(1,1)),gtid(Ex(2,2)),gtid(Ex(1,2))) =  &
                                RDM(gtid(Ex(2,1)),gtid(Ex(1,1)),gtid(Ex(2,2)),gtid(Ex(1,2))) -  &
                                Bra(i)*Ket(j)

                        else
                            !Mixed spin excitation
                            !i has the same spin as b
                            RDM(gtid(Ex(2,1)),gtid(Ex(1,2)),gtid(Ex(2,2)),gtid(Ex(1,1))) = &
                                RDM(gtid(Ex(2,1)),gtid(Ex(1,2)),gtid(Ex(2,2)),gtid(Ex(1,1))) + &
                                Bra(i)*Ket(j)
                            RDM(gtid(Ex(2,2)),gtid(Ex(1,1)),gtid(Ex(2,1)),gtid(Ex(1,2))) = &
                                RDM(gtid(Ex(2,2)),gtid(Ex(1,1)),gtid(Ex(2,1)),gtid(Ex(1,2))) + &
                                Bra(i)*Ket(j)

                        endif

                    else
                        if(mod(Ex(1,1),2).eq.mod(Ex(1,2),2)) then
                            !same spin excitation
                            RDM(gtid(Ex(2,1)),gtid(Ex(1,2)),gtid(Ex(2,2)),gtid(Ex(1,1))) =  &
                                RDM(gtid(Ex(2,1)),gtid(Ex(1,2)),gtid(Ex(2,2)),gtid(Ex(1,1))) -  &
                                Bra(i)*Ket(j)
                            RDM(gtid(Ex(2,2)),gtid(Ex(1,1)),gtid(Ex(2,1)),gtid(Ex(1,2))) =  &
                                RDM(gtid(Ex(2,2)),gtid(Ex(1,1)),gtid(Ex(2,1)),gtid(Ex(1,2))) -  &
                                Bra(i)*Ket(j)
                            RDM(gtid(Ex(2,2)),gtid(Ex(1,2)),gtid(Ex(2,1)),gtid(Ex(1,1))) =  &
                                RDM(gtid(Ex(2,2)),gtid(Ex(1,2)),gtid(Ex(2,1)),gtid(Ex(1,1))) +  &
                                Bra(i)*Ket(j)
                            RDM(gtid(Ex(2,1)),gtid(Ex(1,1)),gtid(Ex(2,2)),gtid(Ex(1,2))) =  &
                                RDM(gtid(Ex(2,1)),gtid(Ex(1,1)),gtid(Ex(2,2)),gtid(Ex(1,2))) +  &
                                Bra(i)*Ket(j)
                        else
                            !Mixed spin excitation
                            !i has the same spin as b
                            RDM(gtid(Ex(2,1)),gtid(Ex(1,2)),gtid(Ex(2,2)),gtid(Ex(1,1))) = &
                                RDM(gtid(Ex(2,1)),gtid(Ex(1,2)),gtid(Ex(2,2)),gtid(Ex(1,1))) - &
                                Bra(i)*Ket(j)
                            RDM(gtid(Ex(2,2)),gtid(Ex(1,1)),gtid(Ex(2,1)),gtid(Ex(1,2))) = &
                                RDM(gtid(Ex(2,2)),gtid(Ex(1,1)),gtid(Ex(2,1)),gtid(Ex(1,2))) - &
                                Bra(i)*Ket(j)
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
                                        Bra(i)*Ket(j)

                                    RDM(gtid(Ex(2,1)),gtid(Ex(1,1)),kel,kel) = &
                                        RDM(gtid(Ex(2,1)),gtid(Ex(1,1)),kel,kel) - &
                                        Bra(i)*Ket(j)

                                    RDM(kel,kel,gtid(Ex(2,1)),gtid(Ex(1,1))) = &
                                        RDM(kel,kel,gtid(Ex(2,1)),gtid(Ex(1,1))) - &
                                        Bra(i)*Ket(j)

                                    RDM(gtid(Ex(2,1)),kel,kel,gtid(Ex(1,1))) = &
                                        RDM(gtid(Ex(2,1)),kel,kel,gtid(Ex(1,1))) + &
                                        Bra(i)*Ket(j)

                                else
                                    !k opposite spin
                                    RDM(gtid(Ex(2,1)),gtid(Ex(1,1)),kel,kel) = &
                                        RDM(gtid(Ex(2,1)),gtid(Ex(1,1)),kel,kel) - &
                                        Bra(i)*Ket(j)
                                    RDM(kel,kel,gtid(Ex(2,1)),gtid(Ex(1,1))) = &
                                        RDM(kel,kel,gtid(Ex(2,1)),gtid(Ex(1,1))) - &
                                        Bra(i)*Ket(j)
                                endif


                            else
                                if(mod(Ex(1,1),2).eq.mod(FCIDetList(k,i),2)) then
                                    !k also same spin
                                    RDM(kel,gtid(Ex(1,1)),gtid(Ex(2,1)),kel) = &
                                        RDM(kel,gtid(Ex(1,1)),gtid(Ex(2,1)),kel) - &
                                        Bra(i)*Ket(j)

                                    RDM(gtid(Ex(2,1)),gtid(Ex(1,1)),kel,kel) = &
                                        RDM(gtid(Ex(2,1)),gtid(Ex(1,1)),kel,kel) + &
                                        Bra(i)*Ket(j)

                                    RDM(kel,kel,gtid(Ex(2,1)),gtid(Ex(1,1))) = &
                                        RDM(kel,kel,gtid(Ex(2,1)),gtid(Ex(1,1))) + &
                                        Bra(i)*Ket(j)

                                    RDM(gtid(Ex(2,1)),kel,kel,gtid(Ex(1,1))) = &
                                        RDM(gtid(Ex(2,1)),kel,kel,gtid(Ex(1,1))) - &
                                        Bra(i)*Ket(j)
                                else
                                    !k opposite spin
                                    RDM(gtid(Ex(2,1)),gtid(Ex(1,1)),kel,kel) = &
                                        RDM(gtid(Ex(2,1)),gtid(Ex(1,1)),kel,kel) + &
                                        Bra(i)*Ket(j)
                                    RDM(kel,kel,gtid(Ex(2,1)),gtid(Ex(1,1))) = &
                                        RDM(kel,kel,gtid(Ex(2,1)),gtid(Ex(1,1))) + &
                                        Bra(i)*Ket(j)
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
                                    Bra(i)*Ket(j)
                                RDM(lel,lel,kel,kel) = RDM(lel,lel,kel,kel) + &
                                    Bra(i)*Ket(j)
                                RDM(lel,kel,kel,lel) = RDM(lel,kel,kel,lel) - &
                                    Bra(i)*Ket(j)
                                RDM(kel,lel,lel,kel) = RDM(kel,lel,lel,kel) - &
                                    Bra(i)*Ket(j)
                            else
                                RDM(kel,kel,lel,lel) = RDM(kel,kel,lel,lel) + &
                                    Bra(i)*Ket(j)
                                RDM(lel,lel,kel,kel) = RDM(lel,lel,kel,kel) + &
                                    Bra(i)*Ket(j)
                            endif

                        enddo
                    enddo
                endif

            enddo
        enddo


    end subroutine FindFull2RDM

    subroutine WriteFCIDUMP()
        use utils, only: get_free_unit
        implicit none
        integer :: iunit,i,j,k,l,alpha,lWork,info,twoESize
        real(dp) :: hel
        real(dp), allocatable :: CoreH(:,:),W(:),work(:)
        character(len=*), parameter :: t_r='WriteFCIDUMP'

        !Open & write header
        iunit = get_free_unit()
        open(iunit,file='FCIDUMP',status='unknown')
        write(iunit,"(A,I9,A,I9,A)") "&FCI NORB=",EmbSize,", NELEC=",Elec,", MS2=0,"
        write(iunit,"(A)",advance='no') "ORBSYM= 1,"
        do i=2,EmbSize-1
            write(iunit,"(A)",advance='no') "1,"
        enddo
        write(iunit,"(A)") "1,"
        write(iunit,"(A)") "ISYM=1"
        write(iunit,"(A)") "&END"

        if(.not.tCoreH_EmbBasis) then
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
        else
            !Transform to the core hamiltonian basis and write out in this basis
            write(6,"(A)") "Transforming to core hamiltonian basis before writing out FCIDUMP..."
            !First, construct the core hamiltonian
            allocate(CoreH(EmbSize,EmbSize))
            allocate(W(EmbSize))
            CoreH(:,:) = Emb_h0v(:,:)
            if(tChemPot) CoreH(1,1) = CoreH(1,1) - (U/2.0_dp)
            W(:) = 0.0_dp

            !Diagonalize
            allocate(work(1))
            lWork=-1
            info = 0
            call dsyev('V','U',EmbSize,CoreH,EmbSize,W,Work,lWork,info)
            if(info.ne.0) call stop_all(t_r,'Workspace query failed')
            lwork = int(work(1))+1
            deallocate(work)
            allocate(work(lwork))
            call dsyev('V','U',EmbSize,CoreH,EmbSize,W,Work,lWork,info)
            if(info.ne.0) call stop_all(t_r,'Diag failed')
            deallocate(work)

            !Write out 2 electron integrals, transforming them into the coreH basis
            if(tAnderson) then
                twoESize = 1
            else
                twoESize = nImp
            endif
            do i = 1,EmbSize
                do j = 1,EmbSize
                    do k = 1,EmbSize
                        do l = 1,EmbSize
                            hel = 0.0_dp
                            do alpha = 1,twoESize
                                hel = hel + CoreH(alpha,i)*CoreH(alpha,j)*CoreH(alpha,k)*CoreH(alpha,l)
                            enddo
                            hel = hel * U
                            if(abs(hel).gt.1.0e-14_dp) then
                                write(iunit,"(F16.12,4I8)") hel,i,j,k,l
                            endif
                        enddo
                    enddo
                enddo
            enddo

            do i = 1,EmbSize
                write(iunit,"(F16.12,4I8)") W(i),i,i,0,0
            enddo

            deallocate(W,CoreH)
        endif

        !Core energy is zero
        write(iunit,"(F16.12,4I8)") 0.0_dp,0,0,0,0
        close(iunit)

    end subroutine WriteFCIDUMP

    !Find logical size of desired compressed hamiltonian
    subroutine CountSizeCompMat(DetList,Elec,nDet,Nmax,BitDetList)
        use DetTools, only: GetHElement
        implicit none
        integer, intent(in) :: Elec,nDet
        integer, intent(out) :: Nmax
        integer, intent(in) :: DetList(Elec,nDet)
        integer, intent(in), optional :: BitDetList(nDet)
        integer :: i,j
        real(dp) :: Elem

        Nmax = nDet + 1 !For storage of diagonal elements

        if(present(BitDetList)) then
            do i = 1,nDet
                do j = i+1,nDet
                    call GetHElement(DetList(:,i),DetList(:,j),Elec,Elem,ilutnI=BitDetList(i),ilutnJ=BitDetList(j))
                    if(abs(Elem).ge.CompressThresh) then
                        Nmax = Nmax + 2     !Due to both sides of matrix
                    endif
                enddo
                if(mod(i,25000).eq.0) then
                    write(6,"(A,2I12)") "Counting Mat: ",i,nDet
                    call flush(6)
                endif
            enddo
        else
            do i = 1,nDet
                do j = i+1,nDet
                    call GetHElement(DetList(:,i),DetList(:,j),Elec,Elem)
                    if(abs(Elem).ge.CompressThresh) then
                        Nmax = Nmax + 2     !Due to both sides of matrix
                    endif
                enddo
                if(mod(i,25000).eq.0) then
                    write(6,"(A,2I12)") "Counting Mat: ",i,nDet
                    call flush(6)
                endif
            enddo
        endif

    end subroutine CountSizeCompMat

    !Create compressed matrix form of hamiltonian from basis elements
    !DetList
    subroutine StoreCompMat(DetList,Elec,nDet,Nmax,sa,ija,BitDetList)
        use DetTools, only: GetHElement
        implicit none
        integer, intent(in) :: Elec,nDet,Nmax
        integer, intent(in) :: DetList(Elec,nDet)
        real(dp), intent(out) :: sa(Nmax)
        integer, intent(out) :: ija(Nmax)
        integer, intent(in), optional :: BitDetList(nDet)
        integer :: i,j,k
        real(dp) :: Elem
        character(len=*), parameter :: t_r='StoreCompMat'

        sa(:) = 0.0_dp  !Compressed matrix
        ija(:) = 0      !Index to compressed matrix

        !Store diagonals
        do j = 1,nDet
            if(present(BitDetList)) then
                call GetHElement(DetList(:,j),DetList(:,j),Elec,sa(j),ilutnI=BitDetList(j),ilutnJ=BitDetList(j))
            else
                call GetHElement(DetList(:,j),DetList(:,j),Elec,sa(j))
            endif
        enddo
        ija(1) = nDet + 2
        k = nDet + 1

        if(present(BitDetList)) then
            do i = 1,nDet
                do j = 1,nDet
                    if(i.eq.j) cycle
                    call GetHElement(DetList(:,i),DetList(:,j),Elec,Elem,ilutnI=BitDetList(i),ilutnJ=BitDetList(j))
                    if(abs(Elem).ge.CompressThresh) then
                        k = k + 1
                        if(k.gt.Nmax) call stop_all(t_r,'Compressed array sizes too small')
                        sa(k) = Elem
                        ija(k) = j
                    endif
                enddo
                ija(i+1) = k + 1
                if(mod(i,25000).eq.0) then
                    write(6,"(A,2I12)") "Building Mat: ",i,nDet
                    call flush(6)
                endif
            enddo
        else
            do i = 1,nDet
                do j = 1,nDet
                    if(i.eq.j) cycle
                    call GetHElement(DetList(:,i),DetList(:,j),Elec,Elem)
                    if(abs(Elem).ge.CompressThresh) then
                        k = k + 1
                        if(k.gt.Nmax) call stop_all(t_r,'Compressed array sizes too small')
                        sa(k) = Elem
                        ija(k) = j
                    endif
                enddo
                ija(i+1) = k + 1
                if(mod(i,25000).eq.0) then
                    write(6,"(A,2I12)") "Building Mat: ",i,nDet
                    call flush(6)
                endif
            enddo
        endif

    end subroutine StoreCompMat
    
    !Create compressed matrix form of hamiltonian from basis elements
    !DetList
    subroutine StoreCompMat_comp(DetList,Elec,nDet,Nmax,sa,ija,BitDetList)
        use DetTools, only: GetHElement_comp
        implicit none
        integer, intent(in) :: Elec,nDet,Nmax
        integer, intent(in) :: DetList(Elec,nDet)
        complex(dp), intent(out) :: sa(Nmax)
        integer, intent(out) :: ija(Nmax)
        integer, intent(in), optional :: BitDetList(nDet)
        integer :: i,j,k
        complex(dp) :: Elem
        character(len=*), parameter :: t_r='StoreCompMat_comp'

        sa(:) = zzero  !Compressed matrix
        ija(:) = 0      !Index to compressed matrix

        !Store diagonals
        do j = 1,nDet
            if(present(BitDetList)) then
                call GetHElement_comp(DetList(:,j),DetList(:,j),Elec,sa(j),ilutnI=BitDetList(j),ilutnJ=BitDetList(j))
            else
                call GetHElement_comp(DetList(:,j),DetList(:,j),Elec,sa(j))
            endif
        enddo
        ija(1) = nDet + 2
        k = nDet + 1

        if(present(BitDetList)) then
            do i = 1,nDet
                do j = 1,nDet
                    if(i.eq.j) cycle
                    call GetHElement_comp(DetList(:,i),DetList(:,j),Elec,Elem,ilutnI=BitDetList(i),ilutnJ=BitDetList(j))
                    if(abs(Elem).ge.CompressThresh) then
                        k = k + 1
                        if(k.gt.Nmax) call stop_all(t_r,'Compressed array sizes too small')
                        sa(k) = Elem
                        ija(k) = j
                    endif
                enddo
                ija(i+1) = k + 1
                if(mod(i,25000).eq.0) then
                    write(6,"(A,2I12)") "Building Mat: ",i,nDet
                    call flush(6)
                endif
            enddo
        else
            do i = 1,nDet
                do j = 1,nDet
                    if(i.eq.j) cycle
                    call GetHElement_comp(DetList(:,i),DetList(:,j),Elec,Elem)
                    if(abs(Elem).ge.CompressThresh) then
                        k = k + 1
                        if(k.gt.Nmax) call stop_all(t_r,'Compressed array sizes too small')
                        sa(k) = Elem
                        ija(k) = j
                    endif
                enddo
                ija(i+1) = k + 1
                if(mod(i,25000).eq.0) then
                    write(6,"(A,2I12)") "Building Mat: ",i,nDet
                    call flush(6)
                endif
            enddo
        endif

    end subroutine StoreCompMat_comp

end module solvers
