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

        !Set the correlation potential to zero over the impurity, since we do not need to include the fitted potential over the 
        !orbitals that we are going to do the high-level calculation on.
        !Emb_h0v is the core hamiltonian in the embedded basis, with correlation potential (not over the impurity sites)
        if(allocated(Emb_h0v)) deallocate(Emb_h0v)
        allocate(Emb_h0v(EmbSize,EmbSize))
        Emb_h0v(:,:) = Emb_CorrPot(:,:)
        Emb_h0v(1:nImp,1:nImp) = 0.0_dp     !Set correlation potential over the impurity sites to zero
        Emb_h0v(:,:) = Emb_h0v(:,:) + Emb_h0(:,:)    !Add the embedded core hamiltonian over all sites

        !Calculate the number of electrons in the embedded system 
        Emb_nElec = 0.0_dp
        do i=1,EmbSize
            Emb_nElec = Emb_nElec + Emb_MF_DM(i,i)
        enddo
        write(6,"(A,I5)") "Number of electrons in full system: ",NEl
        Elec = nint(Emb_nElec)
        write(6,"(A,F10.7,I5)") "Number of electrons in embedded system: ",Emb_nElec,Elec

        if(nSys.gt.EmbSize) call stop_all(t_r,"Error in determining basis")

        if(tGSFCI) then
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

            write(6,"(A,F20.10)") "Embedded system energy is: ",HL_Energy

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

                                HL_2RDM(i,j,k,l) = HL_2RDM(i,j,k,l) + HL_2RDM_temp(((i-1)*EmbSize)+j,((k-1)*EmbSize)+l)

                                if(j.eq.k) then
                                    HL_2RDM(i,j,k,l) = HL_2RDM(i,j,k,l) - HL_1RDM(i,l)
                                endif

                            enddo
                        enddo
                    enddo
                enddo
                deallocate(HL_2RDM_temp)

            endif

        elseif(tCompleteDiag) then
            !Do a complete diagonalization
            !Do not need to write FCIDUMP, since would only read it back in...
            call CompleteDiag(tCreate2RDM)
        endif

        call writematrix(HL_1RDM,'HL_1RDM',.true.)

        !Now calculate the seperate 1-body and 2-body contribution to the energy of the embedded system
        !The one-body contribution can be calculated from Tr[h^T x FCI_RDM] = \sum_ij h_ij * FCI_RDM_ij
        !where h^T is the core hamiltonian in the embedded basis with the correlation potential over all but the impurity (the 1e integrals)
        One_ElecE = 0.0_dp
        do j=1,EmbSize
            do i=1,EmbSize
                One_ElecE = One_ElecE + Emb_h0v(i,j)*HL_1RDM(i,j)
            enddo
        enddo

        !The two electron contribution to the embedded system energy is just the FCI result - the 1e energy
        Two_ElecE = HL_Energy - One_ElecE

        if(tCreate2RDM) then

            if(tCompleteDiag) then
                write(6,*) "FCI determinant GS: "
                do i=1,nFCIDet
                    write(6,*) FCIDetList(:,i),FullHamil(i,1)
                enddo
            endif

            !Write out the 2 electron RDM
            do i=1,EmbSize
                do j=1,EmbSize
                    do k=1,EmbSize
                        do l=1,EmbSize
                            write(6,*) "RDM: ",i,j,k,l,HL_2RDM(i,j,k,l)
                        enddo
                    enddo
                enddo
            enddo


            !Check that 2RDM is correct by explicitly calculating the two-electron contribution to the FCI energy from the 2RDM
            !We only need to check the iiii components over the impurity sites ONLY
            Check2eEnergy = 0.0_dp
            do i=1,nImp
                Check2eEnergy = Check2eEnergy + U*HL_2RDM(i,i,i,i)
            enddo
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
                    trace = trace + HL_2RDM(i,j,j,i)
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
                        trace = trace + HL_2RDM(i,k,k,j)
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
                One_ElecE_Imp = One_ElecE_Imp + Emb_h0v(i,j)*HL_1RDM(i,j)
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

!        if(tMFResponse) then
!            !Calculate the perturbed response
!
!            !Calculate the vector of length 4, H^(1) Psi^(0)
!            allocate(FullH1(2,2))
!            FullH1(:,:) = Emb_h1(:,:) + Emb_Pert(:,:) 
!            FullH1(:,:) = Lambda * FullH1(:,:)  !Add in the pertubation strength
!            if(abs(FullH1(1,2)-FullH1(2,1)).gt.1.0e-7) call stop_all(t_r,'Perturbation not hermitian')
!            allocate(LR_State(4))
!            LR_State(:) = 0.0_dp    !The FCI linear response to the perturbations
!
!            !Need to consider all singles (and diagonals) of all determinants
!            !Det 1 first
!            !Diagonals
!            LR_State(1) = LR_State(1) + 2.0_dp*FullH1(1,1)*FullHamil(1,1)
!            !1 -> 2
!            LR_State(2) = LR_State(2) + FullH1(1,2)*FullHamil(1,1)
!            !1 -> 3
!            LR_State(3) = LR_State(3) - FullH1(1,2)*FullHamil(1,1)
!
!            !Det 2
!            !Diagonals
!            LR_State(2) = LR_State(2) + FullH1(1,1)*FullHamil(2,1) + FullH1(2,2)*FullHamil(2,1)
!            !2 -> 1
!            LR_State(1) = LR_State(1) + FullH1(2,1)*FullHamil(2,1)
!            !2 -> 4
!            LR_State(4) = LR_State(4) + FullH1(1,2)*FullHamil(2,1)
!
!            !Det 3
!            !Diagonals
!            LR_State(3) = LR_State(3) + FullH1(1,1)*FullHamil(3,1) + FullH1(2,2)*FullHamil(3,1)
!            !3 -> 1
!            LR_State(1) = LR_State(1) - FullH1(2,1)*FullHamil(3,1)
!            !3 -> 4
!            LR_State(4) = LR_State(4) - FullH1(1,2)*FullHamil(3,1)
!
!            !Det 4
!            !Diagonals
!            LR_State(4) = LR_State(4) + 2.0_dp*FullH1(2,2)*FullHamil(4,1)
!            !4 -> 2
!            LR_State(2) = LR_State(2) + FullH1(2,1)*FullHamil(4,1)
!            !4 -> 3
!            LR_State(3) = LR_State(3) - FullH1(2,1)*FullHamil(4,1)
!
!            !Now, orthogonalize this state against the original Psi_0
!            Overlap = DDOT(4,LR_State(1:4),1,FullHamil(1:4,1),1)
!            LR_State(:) = LR_State(:) - Overlap*FullHamil(:,1)
!            LR_State(:) = -LR_State(:)
!
!            !Construct H^(0) - (E_0 + Omega)
!            ZerothH(:,:) = 0.0_dp
!
!            ZerothH(1,1) = 2.0_dp*Emb_h0v(1,1) + U - (HL_Energy + Omega)
!            ZerothH(2,2) = Emb_h0v(1,1) + Emb_h0v(2,2) - (HL_Energy + Omega)
!            ZerothH(3,3) = Emb_h0v(1,1) + Emb_h0v(2,2) - (HL_Energy + Omega)
!            ZerothH(4,4) = 2.0_dp*Emb_h0v(2,2) - (HL_Energy + Omega)
!
!            ZerothH(1,2) = Emb_h0v(1,2)
!            ZerothH(1,3) = - Emb_h0v(1,2)
!
!            ZerothH(2,4) = Emb_h0v(1,2)
!            ZerothH(3,4) = - Emb_h0v(1,2)
!            
!            ZerothH(2,1) = Emb_h0v(1,2)
!            ZerothH(4,1) = - Emb_h0v(1,2)
!
!            ZerothH(4,2) = Emb_h0v(1,2)
!            ZerothH(4,3) = - Emb_h0v(1,2)
!
!            !Now, solve the linear equation Ax = b, where b is LR_State, A is ZerothH and x will be the linear response wavefunction
!            call DGESV(4,1,ZerothH,4,Pivots,LR_State,4,info)
!            if(info.ne.0) call stop_all(t_r,'Error with solving linear system')
!
!            !LR_State is now the coefficients of Psi^(1)
!            !Now, apply the perturbation again, and project back onto the zeroth order wavefunction (still in FullHamil)
!            !Emb_Pert is the one-electron operator (not diagonal)
!            DD_Response = 0.0_dp
!            
!            !Det 1
!            !Diagonals
!            DD_Response = DD_Response + FullHamil(1,1)*2.0_dp*Emb_Pert(1,1)*LR_State(1)
!            !1 -> 2
!            DD_Response = DD_Response + FullHamil(2,1)*Emb_Pert(1,2)*LR_State(1)
!            !1 -> 3
!            DD_Response = DD_Response - FullHamil(3,1)*Emb_Pert(1,2)*LR_State(1)
!
!            !Det 2
!            !Diagonals
!            DD_Response = DD_Response + FullHamil(2,1)*Emb_Pert(1,1)*LR_State(2) + FullHamil(2,1)*Emb_Pert(2,2)*LR_State(2)
!            !2 -> 1
!            DD_Response = DD_Response + FullHamil(1,1)*Emb_Pert(2,1)*LR_State(2)
!            !2 -> 4
!            DD_Response = DD_Response + FullHamil(4,1)*Emb_Pert(1,2)*LR_State(2)
!
!            !Det 3
!            !Diagonals
!            DD_Response = DD_Response + FullHamil(3,1)*Emb_Pert(1,1)*LR_State(3) + FullHamil(3,1)*Emb_Pert(2,2)*FullHamil(3,1)
!            !3 -> 1
!            DD_Response = DD_Response - FullHamil(1,1)*Emb_Pert(2,1)*LR_State(3)
!            !3 -> 4
!            DD_Response = DD_Response - FullHamil(4,1)*Emb_Pert(1,2)*LR_State(3)
!
!            !Det 4
!            DD_Response = DD_Response + FullHamil(4,1)*2.0_dp*Emb_Pert(2,2)*LR_State(4)
!            !4 -> 2
!            DD_Response = DD_Response + FullHamil(2,1)*Emb_Pert(2,1)*LR_State(4)
!            !4 -> 3
!            DD_Response = DD_Response - FullHamil(3,1)*Emb_Pert(2,1)*LR_State(4)
!
!            write(6,*) "NEW DD_Response: ",DD_Response,Omega
!            deallocate(LR_State,FullH1)
!
!
!!            allocate(FullH1(2,2))
!!            FullH1(:,:) = Emb_h1(:,:) + Emb_Pert(:,:) 
!!            FullH1(:,:) = Lambda * FullH1(:,:)  !Add in the pertubation strength
!!            allocate(LR_State(4))
!!            LR_State(:) = 0.0_dp    !The FCI linear response to the perturbations
!!            allocate(temp(EmbSize,EmbSize))
!!            do j=2,4
!!                call FindFull1RDM(j,1,Temp_RDM)
!!!                call writematrix(Temp_RDM,'Non-hermit RDM',.true.)
!!                !Now contract this with the perturbations to get the denominator 
!!                call DGEMM('T','N',EmbSize,EmbSize,EmbSize,1.0_dp,FullH1,EmbSize,Temp_RDM,EmbSize,0.0_dp,temp,EmbSize)
!!                E1 = 0.0_dp
!!                do i=1,EmbSize
!!                    E1 = E1 + temp(i,i)
!!                enddo
!!                E1 = E1/(Spectrum(j)-Spectrum(1)-Omega)
!!
!!                !Hermitian conjugate
!!                call FindFull1RDM(1,j,Temp_RDM)
!!                !Now contract this with the perturbations to get the denominator 
!!                call DGEMM('T','N',EmbSize,EmbSize,EmbSize,1.0_dp,FullH1,EmbSize,Temp_RDM,EmbSize,0.0_dp,temp,EmbSize)
!!                E2 = 0.0_dp
!!                do i=1,EmbSize
!!                    E2 = E2 + temp(i,i)
!!                enddo
!!                E2 = E2/(Spectrum(j)-Spectrum(1)+Omega)
!!
!!                LR_State(j) = - E1 - E2
!!            enddo
!!
!!            !LR_State should be normalized - check this
!!            !Actually, it should be intermediately normalized, therefore will not sum squares to one. SHould be ok...
!!
!!            !First order wavefunction is now \sum_i=0,N [ci |Psi_i] where Psi_i are the eigenfunctions of the zeroth order hamiltonian
!!            !To calculate the density-density response function, we need to go through each component of Psi^(1), expand the eigenstate into determinants,
!!            DD_Response = 0.0_dp
!!            do i=2,4    !do H^(0) eigenstate components of Psi^(1)
!!                do j=1,4    !do determinant coefficient j for eigenstate i
!!                    do k=1,4    !do determinant coefficient k for eigenstate 1 (i.e. GS)
!!
!!                        if(((j.eq.1).and.(k.eq.2)).or.((j.eq.2).and.(k.eq.1))) then
!!                            !Transistion 1 -> 2
!!                            DD_Response = DD_Response + LR_State(i)*FullHamil(j,i)*FullHamil(k,1)*Emb_Pert(1,2)
!!                        elseif(((j.eq.1).and.(k.eq.3)).or.((j.eq.3).and.(k.eq.1))) then
!!                            !Transition 1 -> 3
!!                            DD_Response = DD_Response - LR_State(i)*FullHamil(j,i)*FullHamil(k,1)*Emb_Pert(1,2)
!!                        elseif(((j.eq.2).and.(k.eq.4)).or.((j.eq.4).and.(k.eq.2))) then
!!                            !Transition 2 -> 4
!!                            DD_Response = DD_Response + LR_State(i)*FullHamil(j,i)*FullHamil(k,1)*Emb_Pert(1,2)
!!                        elseif(((j.eq.3).and.(k.eq.4)).or.((k.eq.3).and.(j.eq.4))) then
!!                            !Transition 3 -> 4
!!                            DD_Response = DD_Response - LR_State(i)*FullHamil(j,i)*FullHamil(k,1)*Emb_Pert(1,2)
!!                        elseif((j.eq.k).and.(j.eq.1)) then
!!                            !Diagonal terms
!!                            DD_Response = DD_Response + LR_State(i)*FullHamil(j,i)*FullHamil(k,1)*Emb_Pert(1,1)
!!                        elseif((j.eq.k).and.(j.eq.2)) then
!!                            DD_Response = DD_Response + LR_State(i)*FullHamil(j,i)*FullHamil(k,1)*Emb_Pert(1,2)
!!                        elseif((j.eq.k).and.(j.eq.3)) then
!!                            DD_Response = DD_Response + LR_State(i)*FullHamil(j,i)*FullHamil(k,1)*Emb_Pert(1,2)
!!                        elseif((j.eq.k).and.(j.eq.4)) then
!!                            DD_Response = DD_Response + LR_State(i)*FullHamil(j,i)*FullHamil(k,1)*Emb_Pert(2,2)
!!                        endif
!!
!!                    enddo
!!                enddo
!!            enddo
!!
!!            write(6,"(A,F8.4,A,G20.10)") "Density-Density response for frequency ",Omega, " is ",DD_Response
!!
!!            DD_Response = 0.0_dp
!!            do i=2,4
!!                call FindFull1RDM(1,i,Temp_RDM)
!!                Temp_RDM(:,:) = Temp_RDM(:,:)*LR_State(i)
!!                do j=1,2
!!                    do k=1,2
!!                        DD_Response = DD_Response + Temp_RDM(j,k)*Emb_Pert(j,k)
!!                    enddo
!!                enddo
!!            enddo
!!            write(6,"(A,F8.4,A,G20.10)") "Density-Density RDM response for frequency ",Omega, " is ",DD_Response
!
!!            !Now, construct the density matrix from the linear response
!!            Temp_RDM(:,:) = 0.0_dp
!!            !First, diagonal contributions
!!            Temp_RDM(1,1) = Temp_RDM(1,1) + LR_State(1)*LR_State(1)*2.0_dp !alpha alpha and beta beta component of det 1
!!            Temp_RDM(1,1) = Temp_RDM(1,1) + LR_State(2)*LR_State(2)    
!!            Temp_RDM(2,2) = Temp_RDM(2,2) + LR_State(2)*LR_State(2)
!!            Temp_RDM(1,1) = Temp_RDM(1,1) + LR_State(3)*LR_State(3)
!!            Temp_RDM(2,2) = Temp_RDM(2,2) + LR_State(3)*LR_State(3)
!!            Temp_RDM(2,2) = Temp_RDM(2,2) + LR_State(4)*LR_State(4)*2.0_dp !alpha alpha and beta beta component of det 4
!!
!!            !Now the four off diagonal contributions    1 -> 2
!!            Temp_RDM(2,1) = Temp_RDM(2,1) + LR_State(2)*LR_State(1)
!!            Temp_RDM(1,2) = Temp_RDM(1,2) + LR_State(1)*LR_State(2)   ! 2 -> 1
!!
!!            ! 1 -> 3 (Note minus sign)
!!            Temp_RDM(2,1) = Temp_RDM(2,1) - LR_State(3)*LR_State(1)   
!!            Temp_RDM(1,2) = Temp_RDM(1,2) - LR_State(1)*LR_State(3)
!!
!!            ! 2 -> 4
!!            Temp_RDM(2,1) = Temp_RDM(2,1) + LR_State(4)*LR_State(2)
!!            Temp_RDM(1,2) = Temp_RDM(1,2) + LR_State(2)*LR_State(4)
!!
!!            !Finally 3 -> 4 (also minus sign)
!!            Temp_RDM(2,1) = Temp_RDM(2,1) - LR_State(4)*LR_State(3)
!!            Temp_RDM(1,2) = Temp_RDM(1,2) - LR_State(3)*LR_State(4)
!!
!!            !Finally, contract *this* RDM with the perturbation, to hopefully get the coupling
!!            call DGEMM('T','N',EmbSize,EmbSize,EmbSize,1.0_dp,FullH1,EmbSize,Temp_RDM,EmbSize,0.0_dp,temp,EmbSize)
!!            E2 = 0.0_dp
!!            do i=1,EmbSize
!!                E2 = E2 + temp(i,i)
!!            enddo
!!
!!            write(6,"(A,F10.6)") "Couping to linear response function: ",E2
!!            TotalE_Imp = TotalE_Imp + E2
!!
!!            deallocate(FullH1,LR_State,temp)
!
!        endif

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
        call writematrix(HL_1RDM,'hl_1rdm',.true.)

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
        if(tGSFCI) then
        
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
            do i=1,nSites
                write(iunit,"(F16.12,4I8)") U,i,i,i,i
            enddo

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
        do i=1,nImp
            umat(umatind(i,i,i,i)) = U
        enddo
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

        !Now generate all determinants in the active space
        if(allocated(FCIDetList)) deallocate(FCIDetList)
        call GenDets(Elec,EmbSize)
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
        !call writematrix(FullHamil(1:nFCIDet,1:nFCIDet),'FCI hamil',.true.)

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

        call writevector(Spectrum,'FCI Spectrum')
            
        HL_Energy = Spectrum(1) !GS eigenvalue
            
        if(allocated(HL_1RDM)) deallocate(HL_1RDM)
        allocate(HL_1RDM(EmbSize,EmbSize))
        HL_1RDM(:,:) = 0.0_dp

        call FindFull1RDM(1,1,HL_1RDM)

        call writematrix(HL_1RDM,'HL_1RDM',.true.)

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
                    Ex(1) = 0
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
                        !mixed spin excitation
                        !We need to ensure the same electron have the same spin
                        if(mod(Ex(1,1),2).ne.mod(Ex(2,1),2)) then
                            !We need to switch a and b (or i and j)
                            temp = Ex(2,2)
                            Ex(2,2) = Ex(2,1)
                            Ex(2,1) = temp
                            tSign = .not.tSign
                        endif
                    endif
                    if(gtid(Ex(1,1)).eq.1.and.gtid(Ex(2,1)).eq.1.and.gtid(Ex(1,2)).eq.2.and.gtid(Ex(2,2)).eq.2) then
                        write(6,*) "contrib (11,22): ",i,j,tSign,FullHamil(i,StateBra)*FullHamil(j,StateKet)
                    endif
                    if(tSign) then
                        RDM(gtid(Ex(1,1)),gtid(Ex(2,1)),gtid(Ex(1,2)),gtid(Ex(2,2))) =  &
                            RDM(gtid(Ex(1,1)),gtid(Ex(2,1)),gtid(Ex(1,2)),gtid(Ex(2,2))) -  &
                            FullHamil(i,StateBra)*FullHamil(j,StateKet)
                        RDM(gtid(Ex(1,2)),gtid(Ex(2,2)),gtid(Ex(1,1)),gtid(Ex(2,1))) =  &
                            RDM(gtid(Ex(1,2)),gtid(Ex(2,2)),gtid(Ex(1,1)),gtid(Ex(2,1))) -  &
                            FullHamil(i,StateBra)*FullHamil(j,StateKet)
!                        RDM(gtid(Ex(1,1)),gtid(Ex(1,2)),gtid(Ex(2,2)),gtid(Ex(2,1))) =  &
!                            RDM(gtid(Ex(1,1)),gtid(Ex(1,2)),gtid(Ex(2,2)),gtid(Ex(2,1))) -  &
!                            FullHamil(i,StateBra)*FullHamil(j,StateKet)
!                        RDM(gtid(Ex(1,2)),gtid(Ex(1,1)),gtid(Ex(2,2)),gtid(Ex(2,1))) =  &
!                            RDM(gtid(Ex(1,2)),gtid(Ex(1,1)),gtid(Ex(2,2)),gtid(Ex(2,1))) +  &
!                            FullHamil(i,StateBra)*FullHamil(j,StateKet)
!                        RDM(gtid(Ex(1,1)),gtid(Ex(1,2)),gtid(Ex(2,1)),gtid(Ex(2,2))) =  &
!                            RDM(gtid(Ex(1,1)),gtid(Ex(1,2)),gtid(Ex(2,1)),gtid(Ex(2,2))) +  &
!                            FullHamil(i,StateBra)*FullHamil(j,StateKet)
!                        RDM(gtid(Ex(1,2)),gtid(Ex(1,1)),gtid(Ex(2,1)),gtid(Ex(2,2))) =  &
!                            RDM(gtid(Ex(1,2)),gtid(Ex(1,1)),gtid(Ex(2,1)),gtid(Ex(2,2))) -  &
!                            FullHamil(i,StateBra)*FullHamil(j,StateKet)
                    else
                        RDM(gtid(Ex(1,1)),gtid(Ex(2,1)),gtid(Ex(1,2)),gtid(Ex(2,2))) =  &
                            RDM(gtid(Ex(1,1)),gtid(Ex(2,1)),gtid(Ex(1,2)),gtid(Ex(2,2))) +  &
                            FullHamil(i,StateBra)*FullHamil(j,StateKet)
                        RDM(gtid(Ex(1,2)),gtid(Ex(2,2)),gtid(Ex(1,1)),gtid(Ex(2,1))) =  &
                            RDM(gtid(Ex(1,2)),gtid(Ex(2,2)),gtid(Ex(1,1)),gtid(Ex(2,1))) +  &
                            FullHamil(i,StateBra)*FullHamil(j,StateKet)
!                        RDM(gtid(Ex(1,1)),gtid(Ex(1,2)),gtid(Ex(2,2)),gtid(Ex(2,1))) =  &
!                            RDM(gtid(Ex(1,1)),gtid(Ex(1,2)),gtid(Ex(2,2)),gtid(Ex(2,1))) +  &
!                            FullHamil(i,StateBra)*FullHamil(j,StateKet)
!                        RDM(gtid(Ex(1,2)),gtid(Ex(1,1)),gtid(Ex(2,2)),gtid(Ex(2,1))) =  &
!                            RDM(gtid(Ex(1,2)),gtid(Ex(1,1)),gtid(Ex(2,2)),gtid(Ex(2,1))) -  &
!                            FullHamil(i,StateBra)*FullHamil(j,StateKet)
!                        RDM(gtid(Ex(1,1)),gtid(Ex(1,2)),gtid(Ex(2,1)),gtid(Ex(2,2))) =  &
!                            RDM(gtid(Ex(1,1)),gtid(Ex(1,2)),gtid(Ex(2,1)),gtid(Ex(2,2))) -  &
!                            FullHamil(i,StateBra)*FullHamil(j,StateKet)
!                        RDM(gtid(Ex(1,2)),gtid(Ex(1,1)),gtid(Ex(2,1)),gtid(Ex(2,2))) =  &
!                            RDM(gtid(Ex(1,2)),gtid(Ex(1,1)),gtid(Ex(2,1)),gtid(Ex(2,2))) +  &
!                            FullHamil(i,StateBra)*FullHamil(j,StateKet)
                    endif
                elseif(IC.eq.1) then
                    !Connected by a single
                    Ex(1,1) = 1
                    call GetExcitation(FCIDetList(:,i),FCIDetList(:,j),Elec,Ex,tSign)
                    do k=1,elec
                        kel = gtid(FCIDetList(k,i))

                        if(FCIDetList(k,i).ne.Ex(1,1)) then
                            if(kel.eq.1.and.gtid(Ex(2,1)).eq.2.and.gtid(Ex(1,1)).gt.2) then
                                write(6,*) "contrib (11,22): ",i,j,tSign,FullHamil(i,StateBra)*FullHamil(j,StateKet)
                            endif
                            if(tSign) then
                                !RDM(gtid(Ex(2,1)),kel,kel,gtid(Ex(1,1))) = &
                                !    RDM(gtid(Ex(2,1)),kel,kel,gtid(Ex(1,1))) - &
                                !    FullHamil(i,StateBra)*FullHamil(j,StateKet)
                                RDM(kel,kel,gtid(Ex(2,1)),gtid(Ex(1,1))) = &
                                    RDM(kel,kel,gtid(Ex(2,1)),gtid(Ex(1,1))) - &
                                    FullHamil(i,StateBra)*FullHamil(j,StateKet)
                                RDM(gtid(Ex(2,1)),gtid(Ex(1,1)),kel,kel) = &
                                    RDM(gtid(Ex(2,1)),gtid(Ex(1,1)),kel,kel) - &
                                    FullHamil(i,StateBra)*FullHamil(j,StateKet)

!                                RDM(gtid(Ex(1,1)),kel,gtid(Ex(2,1)),kel) = &
!                                    RDM(gtid(Ex(1,1)),kel,gtid(Ex(2,1)),kel) - &
!                                    FullHamil(i,StateBra)*FullHamil(j,StateKet)
!                                RDM(kel,gtid(Ex(1,1)),gtid(Ex(2,1)),kel) = &
!                                    RDM(kel,gtid(Ex(1,1)),gtid(Ex(2,1)),kel) + &
!                                    FullHamil(i,StateBra)*FullHamil(j,StateKet)
!                                RDM(gtid(Ex(1,1)),kel,kel,gtid(Ex(2,1))) = &
!                                    RDM(gtid(Ex(1,1)),kel,kel,gtid(Ex(2,1))) + &
!                                    FullHamil(i,StateBra)*FullHamil(j,StateKet)
!                                RDM(kel,gtid(Ex(1,1)),kel,gtid(Ex(2,1))) = &
!                                    RDM(kel,gtid(Ex(1,1)),kel,gtid(Ex(2,1))) - &
!                                    FullHamil(i,StateBra)*FullHamil(j,StateKet)
                            else
                                !RDM(gtid(Ex(2,1)),kel,kel,gtid(Ex(1,1))) = &
                                !    RDM(gtid(Ex(2,1)),kel,kel,gtid(Ex(1,1))) + &
                                !    FullHamil(i,StateBra)*FullHamil(j,StateKet)
                                RDM(kel,kel,gtid(Ex(2,1)),gtid(Ex(1,1))) = &
                                    RDM(kel,kel,gtid(Ex(2,1)),gtid(Ex(1,1))) + &
                                    FullHamil(i,StateBra)*FullHamil(j,StateKet)
                                RDM(gtid(Ex(2,1)),gtid(Ex(1,1)),kel,kel) = &
                                    RDM(gtid(Ex(2,1)),gtid(Ex(1,1)),kel,kel) + &
                                    FullHamil(i,StateBra)*FullHamil(j,StateKet)
!
!                                RDM(gtid(Ex(1,1)),kel,gtid(Ex(2,1)),kel) = &
!                                    RDM(gtid(Ex(1,1)),kel,gtid(Ex(2,1)),kel) + &
!                                    FullHamil(i,StateBra)*FullHamil(j,StateKet)
!                                RDM(kel,gtid(Ex(1,1)),gtid(Ex(2,1)),kel) = &
!                                    RDM(kel,gtid(Ex(1,1)),gtid(Ex(2,1)),kel) - &
!                                    FullHamil(i,StateBra)*FullHamil(j,StateKet)
!                                RDM(gtid(Ex(1,1)),kel,kel,gtid(Ex(2,1))) = &
!                                    RDM(gtid(Ex(1,1)),kel,kel,gtid(Ex(2,1))) - &
!                                    FullHamil(i,StateBra)*FullHamil(j,StateKet)
!                                RDM(kel,gtid(Ex(1,1)),kel,gtid(Ex(2,1))) = &
!                                    RDM(kel,gtid(Ex(1,1)),kel,gtid(Ex(2,1))) + &
!                                    FullHamil(i,StateBra)*FullHamil(j,StateKet)
                            endif
                        endif
                    enddo
                elseif(IC.eq.0) then
                    !Same det
                    Ex(1,1) = 0
                    if(i.ne.j) call stop_all(t_r,'Error here')
                    do k=1,Elec
                        kel = gtid(FCIDetList(k,i))
                        do l=1,Elec
                            lel = gtid(FCIDetList(l,i))
                            if(FCIDetList(k,i).eq.FCIDetList(l,i)) cycle
                            if(kel.eq.1.and.lel.eq.2) then
                                write(6,*) "contrib (11,22): ",i,FullHamil(i,StateBra)*FullHamil(j,StateKet)
                            endif

                            RDM(kel,kel,lel,lel) = RDM(kel,kel,lel,lel) + &
                                FullHamil(i,StateBra)*FullHamil(j,StateKet)

                        enddo
                    enddo
                endif

            enddo
        enddo

!        do i=1,EmbSize
!            do j=1,EmbSize
!                do k=1,EmbSize
!                    do l=1,EmbSize
!                        if(abs(RDM(i,j,k,l)+RDM(j,i,k,l)).gt.1.0e-7_dp) then
!                            write(6,*) "RDM(i,j,k,l): ",RDM(i,j,k,l)
!                            write(6,*) "RDM(j,i,k,l): ",RDM(j,i,k,l)
!                            call stop_all(t_r,'2RDM not symmetric')
!                        endif
!                        if(abs(RDM(i,j,k,l)+RDM(i,j,l,k)).gt.1.0e-7_dp) then
!                            write(6,*) "RDM(i,j,k,l): ",RDM(i,j,k,l)
!                            write(6,*) "RDM(i,j,l,k): ",RDM(i,j,l,k)
!                            call stop_all(t_r,'2RDM not symmetric')
!                        endif
!                        if(abs(RDM(i,j,k,l)-RDM(j,i,l,k)).gt.1.0e-7_dp) then
!                            write(6,*) "RDM(i,j,k,l): ",RDM(i,j,k,l)
!                            write(6,*) "RDM(j,i,l,k): ",RDM(j,i,l,k)
!                            call stop_all(t_r,'2RDM not symmetric')
!                        endif
!                        if(abs(RDM(i,j,k,l)-RDM(k,l,i,j)).gt.1.0e-7_dp) then
!                            write(6,*) "RDM(i,j,k,l): ",RDM(i,j,k,l)
!                            write(6,*) "RDM(k,l,i,j): ",RDM(k,l,i,j)
!                            call stop_all(t_r,'2RDM not symmetric')
!                        endif
!                        if(abs(RDM(i,j,k,l)+RDM(l,k,i,j)).gt.1.0e-7_dp) then
!                            write(6,*) "RDM(i,j,k,l): ",RDM(i,j,k,l)
!                            write(6,*) "RDM(l,k,i,j): ",RDM(l,k,i,j)
!                            call stop_all(t_r,'2RDM not symmetric')
!                        endif
!                        if(abs(RDM(i,j,k,l)+RDM(k,l,j,i)).gt.1.0e-7_dp) then
!                            write(6,*) "RDM(i,j,k,l): ",RDM(i,j,k,l)
!                            write(6,*) "RDM(k,l,j,i): ",RDM(k,l,j,i)
!                            call stop_all(t_r,'2RDM not symmetric')
!                        endif
!                        if(abs(RDM(i,j,k,l)-RDM(l,k,j,i)).gt.1.0e-7_dp) then
!                            write(6,*) "RDM(i,j,k,l): ",RDM(i,j,k,l)
!                            write(6,*) "RDM(l,k,j,i): ",RDM(l,k,j,i)
!                            call stop_all(t_r,'2RDM not symmetric')
!                        endif
!                    enddo
!                enddo
!            enddo
!        enddo

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
        do i=1,nImp
            write(iunit,"(F16.12,4I8)") U,i,i,i,i
        enddo

        !Now for 1electron contributions
        do i=1,EmbSize
            do j=1,i
                if(abs(Emb_h0v(i,j)).gt.1.0e-10_dp) then
                    write(iunit,"(F16.12,4I8)") Emb_h0v(i,j),i,j,0,0
                endif
            enddo
        enddo
        !Core energy is zero
        write(iunit,"(F16.12,4I8)") 0.0_dp,0,0,0,0
        close(iunit)

    end subroutine WriteFCIDUMP



end module solvers
