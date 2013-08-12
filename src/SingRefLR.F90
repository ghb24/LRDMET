module SingRefLR
    use const
    use timing
    use errors, only: stop_all,warning
    use mat_tools, only: WriteVector,WriteMatrix,WriteVectorInt,WriteMatrixComp,WriteVectorComp,znrm2
    use matrixops, only: z_inv
    use globals

    contains

    subroutine CorrNI_MomGF()
        use utils, only: get_free_unit,append_ext_real
        use LRSolvers, only: WriteKVecHeader,GetNextkVal
        implicit none
        real(dp) :: k_val(LatticeDim),mu,Omega
        integer :: unit_a,unit_b,k,i,ind_1,ind_2,SS_Period
        complex(dp) :: ni_lr,ni_lr_ann,ni_lr_cre,ni_lr_b,ni_lr_ann_b,ni_lr_cre_b
        logical :: tFinishedk
        character(len=64) :: filename,filename_b
        character(len=*), parameter :: t_r='CorrNI_MomGF'

        write(6,"(A)") "Calculating non-interacting single-particle k-dependent greens function, including GS correlation potential"
        
        !Open file
        unit_a = get_free_unit()
        call append_ext_real('CorrNI_k_depGF',U,filename)
        open(unit_a,file=filename,status='unknown')
        write(unit_a,"(A)") "# 1.Omega  2.Re[GF]  3.Im[GF]  4.Im[GF^-]  5.Im[GF^+]"
        if(tUHF) then
            !If doing UHF, open file for beta results
            unit_b = get_free_unit()
            call append_ext_real('CorrNI_LocalGF_b',U,filename_b)
            open(unit_a,file=filename_b,status='unknown')
            write(unit_b,"(A)") "# 1.Omega  2.Re[GF]  3.Im[GF]  4.Im[GF^-]  5.Im[GF^+]"
        endif

        if(.not.tAnderson) then
            !In the hubbard model, apply a chemical potential of U/2
            mu = U/2.0_dp
        else
            mu = 0.0_dp
        endif

        SS_Period = nImp
            
        k = 0
        do while(.true.)
            call GetNextkVal(k,tFinishedk)
            if(tFinishedk) exit

            k_val = KPnts(:,k)
            call WriteKVecHeader(unit_a,k_val)
            if(tUHF) call WriteKVecHeader(unit_b,k_val)

            Omega = Start_Omega
            do while((Omega.lt.max(Start_Omega,End_Omega)+1.0e-5_dp).and.(Omega.gt.min(Start_Omega,End_Omega)-1.0e-5_dp))

                ni_lr_ann = zzero
                ni_lr_cre = zzero

                !Run through the corresponding k-space orbitals
                !ind_1 and ind_2 define the indices of the orbitals corresponding to this kpoint
                ind_1 = ((k-1)*SS_Period) + 1
                ind_2 = SS_Period*k

                do i = ind_1,ind_2  !Run through eigenvectors corresopnding to this kpoint
                    if(KVec_InvEMap(i).lt.nOcc) then
                        !This corresponds to an occupied orbtials
                        !Strictly speaking, we want to sum the weight of the vector over the kpoint,
                        !but since the vector is normalized, we can assume that this is one
                        ni_lr_ann = ni_lr_ann + zone/(dcmplx(Omega+mu,dDelta)-k_HFEnergies(i))
                    else
                        !Virtual orbital
                        ni_lr_cre = ni_lr_cre + zone/(dcmplx(Omega+mu,dDelta)-k_HFEnergies(i))
                    endif
                enddo

!                do i = 1,nOcc
!                    ni_lr_ann = ni_lr_ann + (HFtoKOrbs(k,i)*dconjg(HFtoKOrbs(k,i))) / &
!                        (dcmplx(Omega+mu,dDelta)-HFEnergies(i))
!                enddo
!                do a = nOcc+1,nSites
!                    ni_lr_cre = ni_lr_cre + (HFtoKOrbs(k,a)*dconjg(HFtoKOrbs(k,a))) / &
!                        (dcmplx(Omega+mu,dDelta)-HFEnergies(a))
!                enddo
                ni_lr = ni_lr_ann + ni_lr_cre 
                if(tUHF) then
                    call stop_all(t_r,'Not set up yet')
                    ni_lr_ann_b = zzero
                    ni_lr_cre_b = zzero
!                    do i = 1,nOcc
!                        ni_lr_ann_b = ni_lr_ann_b + (HFtoKOrbs_b(k,i)*dconjg(HFtoKOrbs_b(k,i))) / &
!                            (dcmplx(Omega+mu,dDelta)-HFEnergies_b(i))
!                    enddo
!                    do a = nOcc+1,nSites
!                        ni_lr_cre_b = ni_lr_cre_b + (HFtoKOrbs_b(k,a)*dconjg(HFtoKOrbs_b(k,a))) / &
!                            (dcmplx(Omega+mu,dDelta)-HFEnergies_b(a))
!                    enddo
                    ni_lr_b = ni_lr_ann_b + ni_lr_cre_b
                endif

                !Write out
                write(unit_a,"(5G22.10)") Omega,real(ni_lr),-aimag(ni_lr),-aimag(ni_lr_ann),-aimag(ni_lr_cre)
                if(tUHF) then
                    write(unit_b,"(5G22.10)") Omega,real(ni_lr_b),-aimag(ni_lr_b),-aimag(ni_lr_ann_b),-aimag(ni_lr_cre_b)
                endif

                Omega = Omega + Omega_Step
            enddo

            !End of k-point
            write(unit_a,"(A)") ""
            write(unit_a,"(A)") ""
            if(tUHF) then
                write(unit_b,"(A)") ""
                write(unit_b,"(A)") ""
            endif

        enddo   !End k

        close(unit_a)
        if(tUHF) close(unit_b)

    end subroutine CorrNI_MomGF

    subroutine CorrNI_LocalGF()
        use utils, only: get_free_unit,append_ext_real
        implicit none
        real(dp) :: Omega,mu
        integer :: unit_a,unit_b,i,a,nImp_GF,imp
        complex(dp), allocatable :: ni_lr_ann(:),ni_lr_cre(:),ni_lr(:),ni_lr_ann_b(:),ni_lr_cre_b(:),ni_lr_b(:)
        complex(dp) :: Tr_ni_lr,Tr_ni_lr_b
        character(64) :: filename,filename_b

        write(6,"(A)") "Calculating non-interacting single-particle local greens function including GS correlation potential"
        if(tUHF) write(6,"(A)") "Doing this for both alpha and beta greens functions..."

        !Open file
        unit_a = get_free_unit()
        call append_ext_real('CorrNI_LocalGF',U,filename)
        open(unit_a,file=filename,status='unknown')
        write(unit_a,"(A)") "# 1.Omega  2.Re[GF]  3.Im[GF]  4.Im[GF^-]  5.Im[GF^+]"
        if(tUHF) then
            !If doing UHF, open file for beta results
            unit_b = get_free_unit()
            call append_ext_real('CorrNI_LocalGF_b',U,filename_b)
            open(unit_a,file=filename_b,status='unknown')
            write(unit_b,"(A)") "# 1.Omega  2.Re[GF]  3.Im[GF]  4.Im[GF^-]  5.Im[GF^+]"
        endif
        
        if(.not.tAnderson) then
            !In the hubbard model, apply a chemical potential of U/2
            mu = U/2.0_dp
        else
            mu = 0.0_dp
        endif

        if(tAllImp_LR) then
            nImp_GF = nImp
        else
            nImp_GF = 1
        endif

        allocate(ni_lr_ann(nImp_GF))
        allocate(ni_lr_cre(nImp_GF))
        allocate(ni_lr_ann_b(nImp_GF))
        allocate(ni_lr_cre_b(nImp_GF))
        allocate(ni_lr(nImp_GF))
        allocate(ni_lr_b(nImp_GF))

        Omega = Start_Omega
        do while((Omega.lt.max(Start_Omega,End_Omega)+1.0e-5_dp).and.(Omega.gt.min(Start_Omega,End_Omega)-1.0e-5_dp))

            ni_lr_ann(:) = zzero
            ni_lr_cre(:) = zzero
            do i = 1,nOcc
                do imp = 1,nImp_GF
                    ni_lr_ann(imp) = ni_lr_ann(imp) + dcmplx(HFOrbs(imp,i)**2,zero) / &
                        (dcmplx(Omega+mu,dDelta)-HFEnergies(i))
                enddo
            enddo
            do a = nOcc+1,nSites
                do imp = 1,nImp_GF
                    ni_lr_cre(imp) = ni_lr_cre(imp) + dcmplx(HFOrbs(imp,a)**2,zero) / &
                        (dcmplx(Omega+mu,dDelta)-HFEnergies(a))
                enddo
            enddo
            ni_lr(:) = ni_lr_ann(:) + ni_lr_cre(:)
            if(tUHF) then
                ni_lr_ann_b(:) = zzero
                ni_lr_cre_b(:) = zzero
                do i = 1,nOcc
                    do imp = 1,nImp_GF
                        ni_lr_ann_b(imp) = ni_lr_ann_b(imp) + dcmplx(HFOrbs_b(imp,i)**2,zero) / &
                            (dcmplx(Omega+mu,dDelta)-HFEnergies_b(i))
                    enddo
                enddo
                do a = nOcc+1,nSites
                    do imp = 1,nImp_GF
                        ni_lr_cre_b(imp) = ni_lr_cre_b(imp) + dcmplx(HFOrbs_b(imp,a)**2,zero) / &
                            (dcmplx(Omega+mu,dDelta)-HFEnergies_b(a))
                    enddo
                enddo
                ni_lr_b(:) = ni_lr_ann_b(:) + ni_lr_cre_b(:)
            endif

            Tr_ni_lr = sum(ni_lr) / nImp_GF
            Tr_ni_lr_b = sum(ni_lr_b) / nImp_GF

            !Write out
            write(unit_a,"(3G22.10)",advance='no') Omega,real(Tr_ni_lr),-aimag(Tr_ni_lr)
            do i = 1,nImp_GF-1
                write(unit_a,"(2G22.10)",advance='no') -aimag(ni_lr_ann(i)),-aimag(ni_lr_cre(i))
            enddo
            write(unit_a,"(2G22.10)") -aimag(ni_lr_ann(nImp_GF)),-aimag(ni_lr_cre(nImp_GF))

            if(tUHF) then
                write(unit_b,"(5G22.10)") Omega,real(Tr_ni_lr_b),-aimag(Tr_ni_lr_b)
                do i = 1,nImp_GF-1
                    write(unit_b,"(2G22.10)",advance='no') -aimag(ni_lr_ann_b(i)),-aimag(ni_lr_cre_b(i))
                enddo
                write(unit_b,"(2G22.10)") -aimag(ni_lr_ann_b(nImp_GF)),-aimag(ni_lr_cre_b(nImp_GF))
            endif

            Omega = Omega + Omega_Step
        enddo
        close(unit_a)
        if(tUHF) close(unit_b)
        deallocate(ni_lr_ann,ni_lr_ann_b,ni_lr_cre,ni_lr_cre_b,ni_lr,ni_lr_b)

    end subroutine CorrNI_LocalGF

    subroutine CorrNI_LocalDD()
        use utils, only: get_free_unit,append_ext_real
        implicit none
        real(dp) :: Omega,EDiff
        integer :: unit_a,i,a
        complex(dp) :: ni_lr
        character(64) :: filename

        write(6,"(A)") "Calculating non-interacting two-particle local greens function including GS correlation potential"
        write(6,"(A)") "Doing this for combined alpha and beta spins..."

        !Open file
        unit_a = get_free_unit()
        call append_ext_real('CorrNI_LocalDD',U,filename)
        open(unit_a,file=filename,status='unknown')
        write(unit_a,"(A)") "# 1.Omega  2.Re[GF]  3.Im[GF]  4.Im[GF^-]  5.Im[GF^+]"
        
        Omega = Start_Omega
        do while((Omega.lt.max(Start_Omega,End_Omega)+1.0e-5_dp).and.(Omega.gt.min(Start_Omega,End_Omega)-1.0e-5_dp))


            ni_lr = zzero
            do i = 1,nOcc
                do a = nOcc+1,nSites
                    EDiff = HFEnergies(a)-HFEnergies(i)
                    ni_lr = ni_lr + dcmplx((HFOrbs(1,i)*HFOrbs(1,a))**2,0.0_dp) / &
                        (dcmplx(Omega-EDiff,dDelta))
                    if(tUHF) then
                        EDiff = HFEnergies_b(a)-HFEnergies_b(i)
                        ni_lr = ni_lr + dcmplx((HFOrbs_b(1,i)*HFOrbs_b(1,a))**2,0.0_dp) / &
                            (dcmplx(Omega-EDiff,dDelta))
                    endif
                enddo
            enddo
            if(.not.tUHF) ni_lr = ni_lr * 2.0_dp

            !Write out
            write(unit_a,"(3G22.10)") Omega,real(ni_lr),-aimag(ni_lr)

            Omega = Omega + Omega_Step
        enddo
        close(unit_a)

    end subroutine CorrNI_LocalDD


    subroutine NonInteractingLR()
        use utils, only: get_free_unit,append_ext_real,append_ext
        use DetTools, only: tospat
        implicit none
        integer :: ov_space,virt_start,i,a,a_spat,i_spat,ai_ind,iunit
        integer :: highbound,pertsite 
        real(dp) :: Omega,EDiff
        complex(dp) :: ResponseFn,ResponseFnPosW
        real(dp), allocatable :: transitions(:,:)   !(ov_space,2)   !1 = transition frequencies, 2 = moments
        character(len=64) :: filename,filename2
        !character(len=*), parameter :: t_r='NonInteractingLR'

        write(6,*) "Calculating the non-interacting linear response function"

        !Assume that the perturbation is local to site 1
        !They should all be the same in the NI limit anyway
        pertsite = 1

        !First, just enumerate transitions
        ov_space =2*nOcc*(nSites-nOcc)
        virt_start = (2*nOcc)+1
        allocate(transitions(ov_space,2))
        transitions(:,:) = 0.0_dp
        do i=1,nel
            do a=virt_start,2*nSites
                if(mod(i,2).ne.mod(a,2)) cycle      !Only want same spin excitations 
                ai_ind = ov_space_spinind(a,i)
                i_spat = tospat(i)
                a_spat = tospat(a)

                transitions(ai_ind,1) = FullHFEnergies(a_spat)-FullHFEnergies(i_spat)
                !Now calculate the moment
                transitions(ai_ind,2) = (FullHFOrbs(pertsite,i_spat)*FullHFOrbs(pertsite,a_spat))**2
                !write(6,*) "perturb_mo",i_spat,a_spat,FullHFOrbs(pertsite,i_spat)*FullHFOrbs(pertsite,a_spat)
                !write(6,*) "i,a,e_i,e_a: ",i_spat,a_spat,FullHFEnergies(i_spat),FullHFEnergies(a_spat)
!                write(6,*) "i, a: ",i,a
!                write(6,*) "moment: ",(FullHFOrbs(pertsite,i_spat)*FullHFOrbs(pertsite,a_spat))**2
                !write(6,*) Transitions(ai_ind,1),Transitions(ai_ind,2)
            enddo
        enddo
        call sort_real2(transitions,ov_space,2)

        call append_ext_real('NonInt_Transitions',U,filename)
        if(.not.tHalfFill) then
            !Also append occupation of lattice to the filename
            call append_ext(filename,nOcc,filename2)
        else
            filename2 = filename
        endif

        iunit = get_free_unit()
        open(unit=iunit,file=filename2,status='unknown')
        write(iunit,"(A)") "#Excitation     Transition_Frequency       Transition_Moment"
        do i=1,ov_space
            write(iunit,"(I8,2G22.12)") i,transitions(i,1),transitions(i,2)
        enddo
        close(iunit)
        write(6,*) "First 10 non-interacting transition frequencies: "
        highbound = min(ov_space,10)
        call writevector(transitions(1:highbound,1),'transition frequencies')
        deallocate(transitions)

        write(6,*) "Writing non-interacting linear response function to disk..."

        call append_ext_real('NonInt_DDResponse',U,filename)
        if(.not.tHalfFill) then
            !Also append occupation of lattice to the filename
            call append_ext(filename,nOcc,filename2)
        else
            filename2 = filename
        endif
        open(unit=iunit,file=filename2,status='unknown')
        write(iunit,"(A)") "# Frequency     DD_LinearResponse    DD_LinearResponse_PosW"
!        iunit2 = get_free_unit()
!        call append_ext_real('NonInt_DDResponse_posW',U,filename2)
!        if(.not.tHalfFill) then
!            !Also append occupation of lattice to the filename
!            call append_ext(filename2,nOcc,filename3)
!        else
!            filename3 = filename2
!        endif
!        open(unit=iunit2,file=filename3,status='unknown')
!        write(iunit2,"(A)") "# Frequency     DD_LinearResponse"

        Omega = Start_Omega
        do while((Omega.lt.max(Start_Omega,End_Omega)+1.0e-5_dp).and.(Omega.gt.min(Start_Omega,End_Omega)-1.0e-5_dp))

            ResponseFn = dcmplx(0.0_dp,0.0_dp)
            ResponseFnPosW = dcmplx(0.0_dp,0.0_dp) !Only positive frequency
            do i=1,nel
                do a=virt_start,2*nSites
                    if(mod(i,2).ne.mod(a,2)) cycle      !Only want same spin excitations 
                    i_spat = tospat(i)
                    a_spat = tospat(a)

                    EDiff = FullHFEnergies(a_spat)-FullHFEnergies(i_spat)
                    ResponseFn = ResponseFn + dcmplx((FullHFOrbs(pertsite,a_spat)* &
                        FullHFOrbs(pertsite,i_spat))**2,0.0_dp) /    &
                        (dcmplx(Omega,dDelta)-dcmplx(EDiff,0.0_dp))
                    ResponseFnPosW = ResponseFnPosW + dcmplx((FullHFOrbs(pertsite,a_spat)* &
                        FullHFOrbs(pertsite,i_spat))**2,0.0_dp) / &
                        (dcmplx(Omega,dDelta)-dcmplx(EDiff,0.0_dp))
                    ResponseFn = ResponseFn - dcmplx((FullHFOrbs(pertsite,a_spat)* &
                        FullHFOrbs(pertsite,i_spat))**2,0.0_dp) /    &
                        (dcmplx(Omega,dDelta)+dcmplx(EDiff,0.0_dp))
                enddo
            enddo
            ResponseFn = ResponseFn*Lambda
            write(iunit,"(5G25.10)") Omega,real(ResponseFn),-aimag(ResponseFn),real(ResponseFnPosW),-aimag(ResponseFnPosW)
!            write(iunit2,*) Omega,real(ResponseFnPosW),-aimag(ResponseFnPosW)

            Omega = Omega + Omega_Step

        enddo
        close(iunit)
!        close(iunit2)

    end subroutine NonInteractingLR

    subroutine TDA_LR()
        use utils, only: get_free_unit,append_ext_real,append_ext
        use DetToolsData, only: tmat,umat
        use DetTools, only: tospat,GetHFAntisymInt_spinorb,GetExcitation,GetHFInt_spinorb,umatind,GetHElement
        implicit none
        integer :: ov_space,virt_start,ierr,i,j,n,m,nj_ind,mi_ind,ex(2,2)
        integer :: m_spat,i_spat,lwork,info,k,l,orbpairs,umatsize,ai_ind,a
        integer :: state,iunit,a_spat,highbound,pertsite
        logical :: tSign
        integer, allocatable :: detHF(:),detR(:),detL(:)
        real(dp) :: HEl1,Omega
        complex(dp) :: ResponseFn
        real(dp), allocatable :: A_mat(:,:),W(:),Work(:),temp(:,:),Residues(:)
        real(dp), allocatable :: DM(:,:),DM_conj(:,:),DM_AO(:,:),DM_AO_conj(:,:)
        character(64) :: filename,filename2
        character(len=*), parameter :: t_r='TDA_LR'

        write(6,*) "Calculating the linear response function via the Tamm-Dancoff approximation"

        !Assume that the perturbation is local to site 1
        pertsite = 1

        ov_space =2*nOcc*(nSites-nOcc)
        virt_start = (2*nOcc)+1

        if(.false.) then
            !Temporarily create & store TMAT & UMAT
            if(allocated(tmat)) deallocate(tmat)
            if(allocated(umat)) deallocate(umat)
            allocate(tmat(nSites,nSites))
            allocate(temp(nSites,nSites))
            call dgemm('t','n',nSites,nSites,nSites,1.0_dp,FullHFOrbs,nSites,h0,nSites,0.0_dp,temp,nSites)
            call dgemm('n','n',nSites,nSites,nSites,1.0_dp,temp,nSites,FullHFOrbs,nSites,0.0_dp,tmat,nSites)
            deallocate(temp)
            if(tChemPot) then
                tmat(1,1) = tmat(1,1) - U/2.0_dp
            endif
            OrbPairs = (nSites*(nSites+1))/2
            umatsize = (OrbPairs*(OrbPairs+1))/2 
            allocate(umat(umatsize))
            do i=1,nSites
                do j=1,nSites
                    do k=1,nSites
                        do l=1,nSites
                            ex(1,1) = i*2
                            ex(1,2) = j*2
                            ex(2,1) = k*2
                            ex(2,2) = l*2
                            umat(umatind(i,j,k,l)) = GetHFInt_spinorb(ex,FullHFOrbs)
                        enddo
                    enddo
                enddo
            enddo

            !calculate matrix brute force to check it is right
            allocate(A_mat(ov_space,ov_space),stat=ierr)    !This is the singles hamiltonian
            if(ierr.ne.0) call stop_all(t_r,'alloc error')
            A_mat(:,:) = 0.0_dp

            allocate(detL(nel))
            allocate(detR(nel))
            allocate(detHF(nel))
            do k=1,nel
                detHF(k) = k
            enddo

            do j=1,nel
                do n=virt_start,2*nSites
                    if(mod(j,2).ne.mod(n,2)) cycle      !Only want same spin excitations for j and n
                    detL(:) = detHF(:)
                    do k=1,nel
                        if(detL(k).eq.j) then
                            detL(k) = n
                            exit
                        endif
                    enddo
                    call sort_int(detL,nel)
                    nj_ind = ov_space_spinind(n,j)

                    do i=1,nel
                        do m=virt_start,2*nSites
                            if(mod(i,2).ne.mod(m,2)) cycle      !Only want same spin excitations for i and m
                            detR(:) = detHF(:)
                            do k=1,nel
                                if(detR(k).eq.i) then
                                    detR(k) = m
                                    exit
                                endif
                            enddo
                            call sort_int(detR,nel)
                            mi_ind = ov_space_spinind(m,i)

                            call GetHElement(detL,detR,nel,HEl1)
                            A_mat(nj_ind,mi_ind) = HEl1
                        enddo
                    enddo
                enddo
            enddo

            !Remove the HF energy from the diagonals
            do i=1,ov_space
                A_mat(i,i) = A_mat(i,i) - HFEnergy
            enddo

            !Check that A is hermition 
            do i=1,ov_space
                do j=1,ov_space
                    if(abs(A_mat(i,j)-A_mat(j,i)).gt.1.0e-7) then
                        call stop_all(t_r,'A not hermitian')
                    endif
                enddo
            enddo

            !Diagonalize A
            allocate(Work(1))
            if(ierr.ne.0) call stop_all(t_r,"alloc err")
            allocate(W(ov_space))
            W(:)=0.0_dp
            lWork=-1
            info=0
            call dsyev('V','U',ov_space,A_mat,ov_space,W,Work,lWork,info)
            if(info.ne.0) call stop_all(t_r,'workspace query failed')
            lwork=int(work(1))+1
            deallocate(work)
            allocate(work(lwork))
            call dsyev('V','U',ov_space,A_mat,ov_space,W,Work,lwork,info)
            if (info.ne.0) call stop_all(t_r,"Diag failed")
            deallocate(work)

            write(6,*) "The first 10 transition frequencies are: "
            highbound = min(ov_space,10)
            call writevector(W(1:highbound),'Transition frequencies')

            deallocate(A_Mat,W,umat,tmat,detL,detR,detHF)

        endif
        
        allocate(A_mat(ov_space,ov_space),stat=ierr)    !This is the singles hamiltonian
        if(ierr.ne.0) call stop_all(t_r,'alloc error')
        A_mat(:,:) = 0.0_dp

        !First, construct A and B, which correspond to:
        !  <HF| [a*_i a_m [H, a*_n a_j]] |HF> = A
        do j=1,nel
            ex(1,2) = j     !second index in integral
            do n=virt_start,2*nSites
                if(mod(j,2).ne.mod(n,2)) cycle      !Only want same spin excitations for j and n
                nj_ind = ov_space_spinind(n,j)
                ex(2,2) = n !4th index in integral
                do i=1,nel
                    ex(2,1) = i !3rd index in integral
                    do m=virt_start,2*nSites
                        if(mod(i,2).ne.mod(m,2)) cycle      !Only want same spin excitations for i and m
                        mi_ind = ov_space_spinind(m,i)
                        ex(1,1) = m !First index in integral

                        !Calculate the antisymmetrized integral, < m j || i n > for A_mat and < m n || i j > for B_mat
                        HEl1 = GetHFAntisymInt_spinorb(ex,FullHFOrbs)
                        A_mat(mi_ind,nj_ind) = HEl1
                    enddo
                enddo
            enddo
        enddo

        !Now add the diagonal part to A
        do i=1,nel
            do m=virt_start,2*nSites
                if(mod(i,2).ne.mod(m,2)) cycle  !Only want same spin excitations
                mi_ind = ov_space_spinind(m,i)

                m_spat = tospat(m)
                i_spat = tospat(i)

                A_mat(mi_ind,mi_ind) = A_mat(mi_ind,mi_ind) + (FullHFEnergies(m_spat)-FullHFEnergies(i_spat))
            enddo
        enddo

        !Check that A is hermition 
        do i=1,ov_space
            do j=1,ov_space
                if(abs(A_mat(i,j)-A_mat(j,i)).gt.1.0e-7) then
                    call stop_all(t_r,'A not hermitian')
                endif
            enddo
        enddo

        !Diagonalize A
        allocate(Work(1))
        if(ierr.ne.0) call stop_all(t_r,"alloc err")
        allocate(W(ov_space))
        W(:)=0.0_dp
        lWork=-1
        info=0
        call dsyev('V','U',ov_space,A_mat,ov_space,W,Work,lWork,info)
        if(info.ne.0) call stop_all(t_r,'workspace query failed')
        lwork=int(work(1))+1
        deallocate(work)
        allocate(work(lwork))
        call dsyev('V','U',ov_space,A_mat,ov_space,W,Work,lwork,info)
        if (info.ne.0) call stop_all(t_r,"Diag failed")
        deallocate(work)

        write(6,*) "The first 10 transition frequencies are: "
        highbound = min(ov_space,10)
        call writevector(W(1:highbound),'Transition frequencies')
        write(6,*) "Calculating TDA Transition moments..."
        call flush(6)

        !call writematrix(A_mat,'evecs',.true.)

        !We now have our eigenvalues and eigenvectors. Calculate the response functions
        !First, calculate the residues
        allocate(Residues(ov_space))
        Residues(:) = 0.0_dp
        allocate(DM(nSites,nSites))
        allocate(DM_conj(nSites,nSites))
        allocate(DM_AO(nSites,nSites))
        allocate(DM_AO_conj(nSites,nSites))
        allocate(detL(nel))
        allocate(detHF(nel))
        do k=1,nel
            detHF(k) = k
        enddo
        allocate(temp(nSites,nSites))

        do i=1,nel
            do a=virt_start,2*nSites
                if(mod(i,2).ne.mod(a,2)) cycle  !Only want same spin excitations
                ai_ind = ov_space_spinind(a,i)
                detL(:) = detHF(:)
                do k=1,nel
                    if(detL(k).eq.i) then
                        detL(k) = a
                    endif
                enddo
                call sort_int(detL(:),nel)
                i_spat = tospat(i)
                a_spat = tospat(a)

                !Calculate permutation
                ex(1,1) = 1
                call GetExcitation(detL,detHF,nel,ex,tSign)

                DM(:,:) = 0.0_dp
                DM_conj(:,:) = 0.0_dp
                if(tSign) then
                    DM(a_spat,i_spat) = -1.0_dp
                    DM_conj(i_spat,a_spat) = -1.0_dp
                else
                    DM(a_spat,i_spat) = 1.0_dp
                    DM_conj(i_spat,a_spat) = 1.0_dp
                endif

                !Transfer to AO basis
                call dgemm('n','n',nSites,nSites,nSites,1.0_dp,FullHFOrbs,nSites,DM,nSites,0.0_dp,temp,nSites)
                call dgemm('n','t',nSites,nSites,nSites,1.0_dp,temp,nSites,FullHFOrbs,nSites,0.0_dp,DM_AO,nSites)
                !Do it also for the conjugate density matrix
                call dgemm('n','n',nSites,nSites,nSites,1.0_dp,FullHFOrbs,nSites,DM_conj,nSites,0.0_dp,temp,nSites)
                call dgemm('n','t',nSites,nSites,nSites,1.0_dp,temp,nSites,FullHFOrbs,nSites,0.0_dp,DM_AO_conj,nSites)

                !Extract pertsite,pertsite component
                !Now run over all states
                do state=1,ov_space
                    Residues(state) = Residues(state) + DM_AO(pertsite,pertsite)*DM_AO_conj(pertsite,pertsite)* &
                        A_mat(ai_ind,state)*A_mat(ai_ind,state)
!                    if(abs(A_mat(ai_ind,state)).gt.0.5_dp) then
!                        write(6,*) "i, a: ",i,a
!                        write(6,*) "moment: ",(FullHFOrbs(pertsite,i_spat)*FullHFOrbs(pertsite,a_spat))**2
!                        write(6,*) "Actually: ",DM_AO(pertsite,pertsite)*DM_AO_conj(pertsite,pertsite)
!                        write(6,*) "coefficient: ",A_mat(ai_ind,state)
!                        write(6,*) "ai_ind: ",ai_ind
!                        write(6,*) DM_AO(pertsite,pertsite),DM_AO_conj(pertsite,pertsite)
!                        write(6,*) FullHFOrbs(pertsite,i_spat)*FullHFOrbs(pertsite,a_spat)
!                    endif
                enddo
            enddo
        enddo
        deallocate(detHF,detL,temp,DM,DM_conj,DM_AO,DM_AO_conj)
        Residues(:) = Residues(:)*Lambda
        
        iunit = get_free_unit()
        call append_ext_real('TDA_Transitions',U,filename)
        if(.not.tHalfFill) then
            !Also append occupation of lattice to the filename
            call append_ext(filename,nOcc,filename2)
        else
            filename2 = filename
        endif
        open(unit=iunit,file=filename2,status='unknown')
        write(iunit,"(A)") "#Excitation     Transition_Frequency       Transition_Moment"
        do i=1,ov_space
            write(iunit,"(I8,2G22.12)") i,W(i),Residues(i)
        enddo
        close(iunit)

        write(6,*) "Writing TDA linear response function to disk..."

        call append_ext_real('TDA_DDResponse',U,filename)
        if(.not.tHalfFill) then
            !Also append occupation of lattice to the filename
            call append_ext(filename,nOcc,filename2)
        else
            filename2 = filename
        endif
        open(unit=iunit,file=filename2,status='unknown')
        write(iunit,"(A)") "# Frequency     DD_LinearResponse"

        Omega = Start_Omega
        do while((Omega.lt.max(Start_Omega,End_Omega)+1.0e-5_dp).and.(Omega.gt.min(Start_Omega,End_Omega)-1.0e-5_dp))

            ResponseFn = dcmplx(0.0_dp,0.0_dp)
            do i=1,ov_space
                ResponseFn = ResponseFn + ((dcmplx(Residues(i),0.0_dp)/(dcmplx(Omega,dDelta)-dcmplx(W(i),0.0_dp))) - &
                    (dcmplx(Residues(i),0.0_dp)/(dcmplx(Omega,dDelta)+dcmplx(W(i),0.0_dp))))
            enddo
            write(iunit,*) Omega,real(ResponseFn),-aimag(ResponseFn)

            Omega = Omega + Omega_Step

        enddo
        close(iunit)

        deallocate(A_Mat,W,Residues)

    end subroutine TDA_LR
    
    !Set up the RPA equations and solve them.
    !Finally, create the density-density linear response function from the resulting excitations/deexcitations
    !This is done in the spin-orbital space
    subroutine RPA_LR()
        use utils, only: get_free_unit,append_ext_real,append_ext
        use matrixops, only: d_inv
        use DetTools, only: tospat,GetHFAntisymInt_spinorb
        implicit none
        integer :: ov_space,virt_start,ierr,j,ex(2,2),ex2(2,2),n,i,m,nj_ind,mi_ind,info,lwork
        integer :: m_spat,i_spat,StabilitySize,mu,j_spat,ai_ind,iunit,a,excit,highbound
        integer :: pertsite
        real(dp) :: HEl1,HEl2,X_norm,Y_norm,norm,Energy_stab,DMEl1,DMEl2,Omega
        complex(dp) :: ResponseFn
        real(dp), allocatable :: A_mat(:,:),B_mat(:,:),Stability(:,:),StabilityCopy(:,:),W(:),Work(:)
        real(dp), allocatable :: S_half(:,:),temp(:,:),temp2(:,:),W2(:),X_stab(:,:),Y_stab(:,:)
        real(dp), allocatable :: trans_moment(:),AOMO_Spin(:,:),DM(:,:)
        character(64) :: filename,filename2
        character(len=*), parameter :: t_r='RPA_LR'

        if(tUHF) call stop_all(t_r,'Not sure this routine works with UHF')

        !Assume that the perturbation is local to site 1. This condition may want to be changed in the future
        pertsite = 1

        ov_space = 2*nOcc*(nSites-nOcc)
        virt_start = (2*nOcc)+1
        allocate(A_mat(ov_space,ov_space),stat=ierr)
        allocate(B_mat(ov_space,ov_space),stat=ierr)
        if(ierr.ne.0) call stop_all(t_r,'alloc error')
        A_mat(:,:) = 0.0_dp
        B_mat(:,:) = 0.0_dp

        !First, construct A and B, which correspond to:
        !  <HF| [a*_i a_m [H, a*_n a_j]] |HF> = A
        ! -<HF| [a*_i a_m [H, a*_j a_n]] |HF> = B
        do j=1,nel
            ex(1,2) = j     !second index in integral
            ex2(2,2) = j
            do n=virt_start,2*nSites
                if(mod(j,2).ne.mod(n,2)) cycle      !Only want same spin excitations for j and n
                nj_ind = ov_space_spinind(n,j)
                ex(2,2) = n !4th index in integral
                ex2(1,2) = n
                do i=1,nel
                    ex(2,1) = i !3rd index in integral
                    ex2(2,1) = i
                    do m=virt_start,2*nSites
                        if(mod(i,2).ne.mod(m,2)) cycle      !Only want same spin excitations for i and m
                        mi_ind = ov_space_spinind(m,i)
                        ex(1,1) = m !First index in integral
                        ex2(1,1) = m

                        !Calculate the antisymmetrized integral, < m j || i n > for A_mat and < m n || i j > for B_mat
                        HEl1 = GetHFAntisymInt_spinorb(ex,FullHFOrbs)
                        HEl2 = GetHFAntisymInt_spinorb(ex2,FullHFOrbs)
                        A_mat(mi_ind,nj_ind) = HEl1
                        B_mat(mi_ind,nj_ind) = HEl2
                    enddo
                enddo
            enddo
        enddo

        !Now add the diagonal part to A
        do i=1,nel
            do m=virt_start,2*nSites
                if(mod(i,2).ne.mod(m,2)) cycle  !Only want same spin excitations
                mi_ind = ov_space_spinind(m,i)

                m_spat = tospat(m)
                i_spat = tospat(i)

                A_mat(mi_ind,mi_ind) = A_mat(mi_ind,mi_ind) + (FullHFEnergies(m_spat)-FullHFEnergies(i_spat))
            enddo
        enddo

        !Check that A is hermition and B is symmetric
        do i=1,ov_space
            do j=1,ov_space
                if(abs(B_mat(i,j)-B_mat(j,i)).gt.1.0e-7_dp) then
                    call stop_all(t_r,'B not symmetric')
                endif
                if(abs(A_mat(i,j)-A_mat(j,i)).gt.1.0e-7) then
                    call stop_all(t_r,'A not hermitian')
                endif
            enddo
        enddo

        !Calculate here via direct diagonalization of the stability matrix
        write(6,*) "Calculating RPA from stability matrix"
        !call flush(6)

        !Stability = ( A  B  )
        !            ( B* A* )
        !Assume all integrals real to start with
        StabilitySize=2*ov_space
        allocate(Stability(StabilitySize,StabilitySize),stat=ierr)
        Stability(:,:)=0.0_dp
        Stability(1:ov_space,1:ov_space)=A_mat(1:ov_space,1:ov_space)
        Stability(ov_space+1:StabilitySize,1:ov_space)=B_mat(1:ov_space,1:ov_space)
        Stability(1:ov_space,ov_space+1:StabilitySize)=B_mat(1:ov_space,1:ov_space)
        Stability(ov_space+1:StabilitySize,ov_space+1:StabilitySize)=A_mat(1:ov_space,1:ov_space)

        !Now diagonalize
        !Find optimal space
        allocate(StabilityCopy(StabilitySize,StabilitySize))
        StabilityCopy(:,:)=Stability(:,:)
        allocate(W(StabilitySize),stat=ierr)    !Eigenvalues of stability matrix
        allocate(Work(1))
        if(ierr.ne.0) call stop_all(t_r,"alloc err")
        W(:)=0.0_dp
        lWork=-1
        info=0
        call dsyev('V','U',StabilitySize,Stability,StabilitySize,W,Work,lWork,info)
        if(info.ne.0) call stop_all(t_r,'workspace query failed')
        lwork=int(work(1))+1
        deallocate(work)
        allocate(work(lwork))
        call dsyev('V','U',StabilitySize,Stability,StabilitySize,W,Work,lwork,info)
        if (info.ne.0) call stop_all(t_r,"Diag failed")
        deallocate(work)

        do i=1,StabilitySize
            if(W(i).lt.0.0_dp) then
                write(6,*) i,W(i)
                call warning(t_r,"HF solution not stable. RPA failed. Recompute HF.")
                deallocate(Stability,W,StabilityCopy,A_mat,B_mat)
                return
            endif
        enddo
        write(6,"(A)") "Stability matrix positive definite. HF solution is minimum. RPA stable"

        !Now compute S^(1/2), and transform into original basis
        allocate(S_half(StabilitySize,StabilitySize))
        S_half(:,:)=0.0_dp
        do i=1,StabilitySize
            S_half(i,i)=sqrt(W(i))
        enddo
        allocate(temp(StabilitySize,StabilitySize))
        call dgemm('n','n',StabilitySize,StabilitySize,StabilitySize,1.0_dp,Stability,StabilitySize,    &
            S_half,StabilitySize,0.0_dp,temp,StabilitySize)
        call dgemm('n','t',StabilitySize,StabilitySize,StabilitySize,1.0_dp,temp,StabilitySize,Stability,   &
            StabilitySize,0.0_dp,S_half,StabilitySize)
        !S_half is now S^1/2 in the original basis

        !Check this by squaring it.
        call dgemm('n','t',StabilitySize,StabilitySize,StabilitySize,1.0_dp,S_half,StabilitySize,S_half,    &
            StabilitySize,0.0_dp,temp,StabilitySize)
        do i=1,StabilitySize
            do j=1,StabilitySize
                if(abs(StabilityCopy(i,j)-temp(i,j)).gt.1.0e-7) then
                    call stop_all(t_r,'S^1/2 not calculated correctly in original basis')
                endif
            enddo
        enddo

        temp(:,:)=0.0_dp
        do i=1,ov_space
            temp(i,i)=1.0_dp
        enddo
        do i=ov_space+1,StabilitySize
            temp(i,i)=-1.0_dp
        enddo
        allocate(temp2(StabilitySize,StabilitySize))
        call dgemm('n','n',StabilitySize,StabilitySize,StabilitySize,1.0_dp,S_half,StabilitySize,temp,  &
            StabilitySize,0.0_dp,temp2,StabilitySize)
        call dgemm('n','n',StabilitySize,StabilitySize,StabilitySize,1.0_dp,temp2,StabilitySize,S_half, &
            StabilitySize,0.0_dp,temp,StabilitySize)
        !Now diagonalize temp = S^(1/2) (1 0 \\ 0 -1 ) S^(1/2)

        lWork=-1
        allocate(W2(StabilitySize))
        allocate(Work(1))
        W2(:)=0.0_dp
        call dsyev('V','U',StabilitySize,temp,StabilitySize,W2,Work,lWork,info)
        if(info.ne.0) call stop_all(t_r,'workspace query failed')
        lwork=int(work(1))+1
        deallocate(work)
        allocate(work(lwork))
        call dsyev('V','U',StabilitySize,temp,StabilitySize,W2,Work,lwork,info)
        if (info.ne.0) call stop_all(t_r,"Diag failed")
        deallocate(work)
!            call writevector(W2,'Excitation energies')
        ! temp now holds the eigenvectors X~ Y~
        ! W2 runs over StabilitySize eigenvalues (ov_space*2). Therefor we expect redundant pairs of +-W2, corresponding
        ! to pairs of eigenvectors (X^v Y^v) and (X^v* Y^v*) (Same in real spaces).
        do i=1,ov_space
            !This they are listed in order of increasing eigenvalue, we should be able to easily check that they pair up
            if(abs(W2(i)+W2(StabilitySize-i+1)).gt.1.0e-7_dp) then
                write(6,*) i,StabilitySize-i+1, W2(i), W2(StabilitySize-i+1), abs(W2(i)-W2(StabilitySize-i+1))
                call stop_all(t_r,"Excitation energy eigenvalues do not pair")
            endif
        enddo

        !We actually have everything we need for the energy already now. However, calculate X and Y too.
        !Now construct (X Y) = S^(-1/2) (X~ Y~)
        !First get S^(-1/2) in the original basis
        S_half(:,:)=0.0_dp
        do i=1,StabilitySize
            S_half(i,i)=-sqrt(W(i))
        enddo
        call dgemm('n','n',StabilitySize,StabilitySize,StabilitySize,1.0_dp,Stability,StabilitySize,S_half, &
            StabilitySize,0.0_dp,temp2,StabilitySize)
        call dgemm('n','t',StabilitySize,StabilitySize,StabilitySize,1.0_dp,temp2,StabilitySize,Stability,  &
            StabilitySize,0.0_dp,S_half,StabilitySize)
        !S_half is now S^(-1/2) in the original basis

        !Now multiply S^(-1/2) (X~ y~)
        call dgemm('n','n',StabilitySize,StabilitySize,StabilitySize,1.0_dp,S_half,StabilitySize,temp,  &
            StabilitySize,0.0_dp,temp2,StabilitySize)

        !Check that eigenvectors are also paired.
        !Rotations among degenerate sets will screw this up though
!            do i=1,ov_space
!                write(6,*) "Eigenvectors: ",i,StabilitySize-i+1,W2(i),W2(StabilitySize-i+1)
!                do j=1,StabilitySize
!                    write(6,*) j,temp2(j,i),temp2(j,StabilitySize-i+1)
!                enddo
!            enddo
!            call writematrix(temp2,'X Y // Y X',.true.)
        !temp2 should now be a matrix of (Y X)
!                                            (X Y)
!           This is the other way round to normal, but due to the fact that our eigenvalues are ordered -ve -> +ve
!           TODO: Are the signs of this matrix correct?
        allocate(X_stab(ov_space,ov_space)) !First index is (m,i) compound index. Second is the eigenvector index.
        allocate(Y_stab(ov_space,ov_space))
        X_stab(:,:)=0.0_dp
        Y_stab(:,:)=0.0_dp
        !Put the eigenvectors corresponding to *positive* eigenvalues into the X_stab and Y_stab arrays.
        X_stab(1:ov_space,1:ov_space)=temp2(1:ov_space,ov_space+1:StabilitySize)
        Y_stab(1:ov_space,1:ov_space)=-temp2(ov_space+1:StabilitySize,ov_space+1:StabilitySize)
        deallocate(temp2)

        !Normalize the eigenvectors appropriately
        do mu=1,ov_space
            norm=0.0_dp
            Y_norm = 0.0_dp
            X_norm = 0.0_dp
            do i=1,ov_space
                norm = norm + X_stab(i,mu)*X_stab(i,mu) - Y_stab(i,mu)*Y_stab(i,mu)
                Y_norm = Y_norm + Y_stab(i,mu)*Y_stab(i,mu)
                X_norm = X_norm + X_stab(i,mu)*X_stab(i,mu)
            enddo
            if(norm.le.0.0_dp) then
                write(6,*) "Norm^2 for vector ",mu," is: ",norm
                call stop_all(t_r,'norm undefined')
            endif
            norm = sqrt(norm)
            do i=1,ov_space
                X_stab(i,mu) = X_stab(i,mu)/norm
                Y_stab(i,mu) = Y_stab(i,mu)/norm
            enddo
            if(Y_norm.gt.X_norm/2.0_dp) then
                write(6,*) "Warning: hole amplitudes large for excitation: ",mu,    &
                    " Quasi-boson approximation breaking down."
                write(6,*) "Norm of X component: ",X_norm
                write(6,*) "Norm of Y component: ",Y_norm
            endif
        enddo
!            call writematrix(X_stab,'X',.true.)

        !Now check orthogonality 
        !call Check_XY_orthogonality(X_stab,Y_stab)

!            call writevector(W2,'Stab_eigenvalues')

        !Now check that we satisfy the original RPA equations
        !For the *positive* eigenvalue space (since we have extracted eigenvectors corresponding to this), check that:
!           ( A B ) (X) = E_v(X )
!           ( B A ) (Y)      (-Y)
        deallocate(temp)
        allocate(temp(StabilitySize,ov_space))
        temp(1:ov_space,1:ov_space) = X_stab(:,:)
        temp(ov_space+1:StabilitySize,1:ov_space) = Y_stab(:,:)
        allocate(temp2(StabilitySize,ov_space))
        call dgemm('n','n',StabilitySize,ov_space,StabilitySize,1.0_dp,StabilityCopy,StabilitySize,temp,    &
            StabilitySize,0.0_dp,temp2,StabilitySize)
        do i=1,ov_space
            do j=1,ov_space
                if(abs(temp2(j,i)-(W2(i+ov_space)*X_stab(j,i))).gt.1.0e-6_dp) then
                    write(6,*) i,j,temp2(j,i),(W2(i+ov_space)*X_stab(j,i)),W2(i+ov_space)
                    call stop_all(t_r,"RPA equations not satisfied for positive frequencies in X matrix")
                endif
            enddo
        enddo
        do i=1,ov_space
            do j=1,ov_space
                if(abs(temp2(j+ov_space,i)-(-W2(i+ov_space)*Y_stab(j,i))).gt.1.0e-6_dp) then
                    write(6,*) i,j,temp2(j+ov_space,i),(-W2(i+ov_space)*Y_stab(j,i)),-W2(i+ov_space)
                    call stop_all(t_r,"RPA equations not satisfied for positive frequencies in Y matrix")
                endif
            enddo
        enddo
        deallocate(temp,temp2)

        !Is is also satisfied the other way around?
        !Check that we also satisfy (still for the *positive* eigenvalues):
!           ( A B ) (Y) = -E_v(Y )
!           ( B A ) (X)       (-X)
        allocate(temp(StabilitySize,ov_space))
        temp(1:ov_space,1:ov_space) = Y_stab(:,:)
        temp(ov_space+1:StabilitySize,1:ov_space) = X_stab(:,:)
        allocate(temp2(StabilitySize,ov_space))
        call dgemm('n','n',StabilitySize,ov_space,StabilitySize,1.0_dp,StabilityCopy,StabilitySize,temp,    &
            StabilitySize,0.0_dp,temp2,StabilitySize)
        do i=1,ov_space
            do j=1,ov_space
                if(abs(temp2(j,i)-(-W2(i+ov_space)*Y_stab(j,i))).gt.1.0e-6_dp) then
                    call stop_all(t_r,"RPA equations not satisfied for negative frequencies in X matrix")
                endif
            enddo
        enddo
        do i=1,ov_space
            do j=1,ov_space
                if(abs(temp2(j+ov_space,i)-(W2(i+ov_space)*X_stab(j,i))).gt.1.0e-6_dp) then
                    call stop_all(t_r,"RPA equations not satisfied for negative frequencies in Y matrix")
                endif
            enddo
        enddo
        deallocate(temp,temp2)
        !TODO: Finally, check that we satisfy eq. 1 in the Scuseria paper for X and Y defined for positive eigenvalues...
        do i=ov_space+1,StabilitySize
            do j=1,StabilitySize
                StabilityCopy(i,j)=-StabilityCopy(i,j)
            enddo
        enddo
        !Stability copy is now (A B // -B -A)
        allocate(temp(StabilitySize,ov_space))
        allocate(temp2(StabilitySize,ov_space))
        temp=0.0_dp
        temp(1:ov_space,1:ov_space) = X_stab(:,:)
        temp(ov_space+1:StabilitySize,1:ov_space) = Y_stab(:,:)
        call dgemm('n','n',StabilitySize,ov_space,StabilitySize,1.0_dp,StabilityCopy,StabilitySize,temp,    &
            StabilitySize,0.0_dp,temp2,StabilitySize)
        do i=1,ov_space
            do j=1,ov_space
                if(abs(temp2(j,i)-(W2(i+ov_space)*X_stab(j,i))).gt.1.0e-7) then
                    call stop_all(t_r,"RPA equations not satisfied for X")
                endif
            enddo
        enddo
        do i=1,ov_space
            do j=ov_space+1,StabilitySize
                if(abs(temp2(j,i)-(W2(i+ov_space)*Y_stab(j-ov_space,i))).gt.1.0e-7) then
                    call stop_all(t_r,"RPA equations not satisfied for Y")
                endif
            enddo
        enddo
        deallocate(temp,temp2)

        !Now calculate energy, in two different ways:
        !1. -1/2 Tr[A] + 1/2 sum_v E_v(positive)
        Energy_stab=0.0_dp
        do i=1,ov_space
            Energy_stab = Energy_stab + W2(ov_space+i) - A_mat(i,i)
        enddo
        Energy_stab = Energy_stab/2.0_dp

        write(6,"(A,G25.10)") "Full RPA energy from stability analysis (plasmonic RPA-TDA excitation energies): ",  &
            Energy_stab

        Energy_stab = 0.0_dp
        !E = 0.25 * Tr[BZ] where Z = Y X^-1

        allocate(temp2(ov_space,ov_space))
        temp2(:,:) = 0.0_dp
        !Find X^-1 
        call d_inv(X_stab,temp2)
!            call writematrix(temp2,'X^-1',.true.)
        allocate(temp(ov_space,ov_space))
        !Find Z (temp)
        call dgemm('n','n',ov_space,ov_space,ov_space,1.0_dp,Y_stab,ov_space,temp2,ov_space,0.0_dp,temp,ov_space)
        !Find BZ (temp2)
        call dgemm('n','n',ov_space,ov_space,ov_space,1.0_dp,B_mat,ov_space,temp,ov_space,0.0_dp,temp2,ov_space)
        !Take trace of BZ
        do i=1,ov_space
            Energy_stab = Energy_stab + temp2(i,i)
        enddo
        Energy_stab = Energy_stab/2.0_dp
        write(6,"(A,G25.10)") "Full RPA energy from stability analysis (Ring-CCD: 1/2 Tr[BZ]): ",Energy_stab

        Energy_stab = 0.0_dp
        do i=1,ov_space
            Y_norm = 0.0_dp
            do j=1,ov_space
                Y_norm = Y_norm + Y_stab(j,i)**2
            enddo
            Energy_stab = Energy_stab - W2(i+ov_space)*Y_norm
        enddo
        write(6,"(A,G25.10)") "Full RPA energy from stability analysis (Y-matrix): ",Energy_stab

        !Now, calculate the response functions for expectation value A and perturbation V
        !This is (for positive frequencies):
        !\sum_nu (<0|[A,Q_nu^+]|0><0|[Q_nu,V]|0> / (omega - W_nu)) - (<0|[V,Q_nu^+]|0><0|[Q_nu,A]|0> / (omega + W_nu))
        !Calculate the transition moments first, for a density-density response at site pertsite 
        allocate(trans_moment(ov_space))    
        trans_moment(:) = 0.0_dp

        !Construct an MO-AO orbital rotation matrix for spin-orbitals
        allocate(AOMO_Spin(nSites*2,nSites*2))
        AOMO_Spin(:,:) = 0.0_dp
        do i=1,nSites*2
            do j=1,nSites*2
                i_spat = tospat(i)
                j_spat = tospat(j)
                AOMO_Spin(j,i) = FullHFOrbs(tospat(j),tospat(i))
            enddo
        enddo
        allocate(DM(nSites*2,nSites*2))
        deallocate(temp)
        allocate(temp(nSites*2,nSites*2))

        write(6,*) "Calculating RPA Transition moments..."
        call flush(6)

        !Calculate <|[V,Q_nu^+]|0><0|[Q_nu,V]|0> and store for each nu
        do i=1,nel
            do a=virt_start,2*nSites
                if(mod(i,2).ne.mod(a,2)) cycle      !Only want same spin excitations 
                ai_ind = ov_space_spinind(a,i)      !This is the index in the array

                DM(:,:) = 0.0_dp
                if(mod(i,2).eq.1) then
                    !Alpha -> alpha transition.
                    !Parity is -1
                    DM(i,a) = -1.0_dp
                else
                    !Beta -> beta transition
                    !Parity is 1
                    DM(i,a) = 1.0_dp
                endif
                !Now, rotate this density matrix into the AO basis
                call dgemm('n','n',nSites*2,nSites*2,nSites*2,1.0_dp,AOMO_Spin,nSites*2,DM,nSites*2,0.0_dp,temp,nSites*2)
                call dgemm('n','t',nSites*2,nSites*2,nSites*2,1.0_dp,temp,nSites*2,AOMO_Spin,nSites*2,0.0_dp,DM,nSites*2)

                !Extract the pertsite,pertsite component, and sum the two spins
!                pertsite_alpha = pertsite*2 - 1
!                pertsite_beta = pertsite*2 
                DMEl1 = DM(pertsite*2-1,pertsite*2-1) + DM(pertsite*2,pertsite*2)

                !Now do <D_i^a|a_a^+ a_i|D_0> element, which is the one on the other side
                DM(:,:) = 0.0_dp
                if(mod(i,2).eq.1) then
                    DM(a,i) = -1.0_dp
                else
                    DM(i,a) = 1.0_dp
                endif
                !Now, rotate this density matrix into the AO basis
                call dgemm('n','n',nSites*2,nSites*2,nSites*2,1.0_dp,AOMO_Spin,nSites*2,DM,nSites*2,0.0_dp,temp,nSites*2)
                call dgemm('n','t',nSites*2,nSites*2,nSites*2,1.0_dp,temp,nSites*2,AOMO_Spin,nSites*2,0.0_dp,DM,nSites*2)

                !Extract the pertsite,pertsite component, and sum the two spins
                DMEl2 = DM(pertsite*2-1,pertsite*2-1) + DM(pertsite*2,pertsite*2)

                do excit=1,ov_space
                    trans_moment(excit) = trans_moment(excit) + ((X_stab(ai_ind,excit)*DMEl1 - Y_stab(ai_ind,excit)*DMEl2)*  &
                        (X_stab(ai_ind,excit)*DMEl2 - Y_stab(ai_ind,excit)*DMEl1))/4.0_dp   !Divide by 4 since each commutator is /2
                enddo
            enddo
        enddo
        trans_moment(:) = trans_moment(:)*Lambda

        iunit = get_free_unit()
        call append_ext_real('RPA_Transitions',U,filename)
        if(.not.tHalfFill) then
            !Also append occupation of lattice to the filename
            call append_ext(filename,nOcc,filename2)
        else
            filename2 = filename
        endif
        open(unit=iunit,file=filename2,status='unknown')
        write(iunit,"(A)") "#Excitation     Transition_Frequency       Transition_Moment"
        do i=1,ov_space
            write(iunit,"(I8,2G22.12)") i,W2(ov_space+i),trans_moment(i)
        enddo
        close(iunit)
        write(6,*) "First 10 RPA transition frequencies: "
        highbound = min(ov_space,10)
        call writevector(W2(ov_space+1:ov_space+highbound),'transition frequencies')

        write(6,*) "Writing RPA linear response function to disk..."

        call append_ext_real('RPA_DDResponse',U,filename)
        if(.not.tHalfFill) then
            !Also append occupation of lattice to the filename
            call append_ext(filename,nOcc,filename2)
        else
            filename2 = filename
        endif
        open(unit=iunit,file=filename2,status='unknown')
        write(iunit,"(A)") "# Frequency     DD_LinearResponse"

        Omega = Start_Omega
        do while((Omega.lt.max(Start_Omega,End_Omega)+1.0e-5_dp).and.(Omega.gt.min(Start_Omega,End_Omega)-1.0e-5_dp))

            ResponseFn = dcmplx(0.0_dp,0.0_dp)
            do i=1,ov_space
                ResponseFn = ResponseFn + ((dcmplx(trans_moment(i),0.0_dp)/    &
                    (dcmplx(Omega,dDelta)-dcmplx(W2(ov_space+i),0.0_dp))) &
                    - (dcmplx(trans_moment(i),0.0_dp)/(dcmplx(Omega,dDelta)+dcmplx(W2(ov_space+i),0.0_dp))))
            enddo
            write(iunit,*) Omega,real(ResponseFn),-aimag(ResponseFn)

            Omega = Omega + Omega_Step

        enddo
        close(iunit)

        deallocate(W2,W,temp,temp2,StabilityCopy,Stability,A_mat,B_Mat,trans_moment,S_half)
        deallocate(X_stab,Y_stab,AOMO_Spin,DM)

    end subroutine RPA_LR

    !Only want to consider single excitation space, consisting of i -> a
    !First list all alpha excitations, then beta excitations
    !Within each spin-type, it is virtual fast
    integer function ov_space_spinind(a,i)
        use DetTools, only: tospat
        implicit none
        integer, intent(in) :: i,a
        integer :: a_spat,i_spat,nVirt_spat

        if(mod(i,2).ne.mod(a,2)) ov_space_spinind = -1  !*Should* be easy to see where this goes wrong

        !Convert to spatial. Index the virtual excitations starting at 1
        a_spat = tospat(a-NEl)    !Runs from 1 -> number of spatial virtual orbitals
        i_spat = tospat(i)        !Runs from 1 -> nOcc
        nVirt_spat = nSites - nOcc

        if(mod(i,2).eq.1) then
            !It is an alpha -> alpha transition
            !These are indexed first
            ov_space_spinind = (i_spat-1)*nVirt_spat + a_spat
        else
            !It is a beta -> beta transition
            !Add on the entire set of alpha -> alpha transitions which come first
            ov_space_spinind = (i_spat-1)*nVirt_spat + a_spat + (nVirt_spat*nOcc)
        endif
    end function ov_space_spinind

    !dPsi/dLambda for Static DD response (i.e. omega -> 0)
    subroutine StaticMF_DD()
        real(dp) , allocatable :: Orbs(:,:),Energies(:),Work(:),PertHamil(:,:),PertOrbs(:,:)
        real(dp) , allocatable :: TempRDM(:,:),PertBath(:),GSBath(:),PertDM(:,:)
        real(dp) :: StaticResponse,Overlap,dStep,DDOT,PertNorm
        character(len=*), parameter :: t_r='StaticMF_DD'
        integer :: lWork, info,i,j,pertsite

        !Assume that the perturbation is local to site 1
        pertsite = 1
        
        allocate(Orbs(nSites,nSites))
        allocate(Energies(nSites))
        allocate(TempRDM(nSites,nSites))
        allocate(GSBath(nSites))
        Orbs(:,:) = h0(:,:)
        Energies(:) = 0.0_dp
        allocate(Work(1))
        lWork=-1
        info=0
        call dsyev('V','U',nSites,Orbs,nSites,Energies,Work,lWork,info)
        if(info.ne.0) call stop_all(t_r,'Workspace queiry failed')
        lwork=int(work(1))+1
        deallocate(work)
        allocate(work(lwork))
        call dsyev('V','U',nSites,Orbs,nSites,Energies,Work,lWork,info)
        if(info.ne.0) call stop_all(t_r,'Diag failed')
        deallocate(work)
        
        !Determine GS bath orbital
        call DGEMM('N','T',nSites,nSites,nOcc,2.0_dp,Orbs(:,1:nOcc),nSites,Orbs(:,1:nOcc),nSites,0.0_dp,TempRDM,nSites)
        GSBath(:) = 0.0_dp
        GSBath(nImp+1:nSites) = TempRDM(nImp+1:nSites,1)
        PertNorm = DDOT(nSites,GSBath(:),1,GSBath(:),1)
        GSBath(:) = GSBath(:) / sqrt(PertNorm)
        deallocate(TempRDM)
        
        dStep = 0.01

        do while(.true.)

            dStep = dStep/2.0_dp

            if(dStep.lt.1.0e-8_dp) exit

            allocate(PertBath(nSites))
            allocate(PertHamil(nSites,nSites))
            allocate(PertDM(nSites,nSites))
            PertHamil(:,:) = h0(:,:)
            PertHamil(pertsite,pertsite) = PertHamil(pertsite,pertsite) + dStep
            Energies(:) = 0.0_dp
            allocate(Work(1))
            lWork=-1
            info=0
            call dsyev('V','U',nSites,PertHamil,nSites,Energies,Work,lWork,info)
            if(info.ne.0) call stop_all(t_r,'Workspace queiry failed')
            lwork=int(work(1))+1
            deallocate(work)
            allocate(work(lwork))
            call dsyev('V','U',nSites,PertHamil,nSites,Energies,Work,lWork,info)
            if(info.ne.0) call stop_all(t_r,'Diag failed')
            deallocate(work)

            !Determine perturbed bath orbital
            !Calc RDM
            call DGEMM('N','T',nSites,nSites,nOcc,2.0_dp,PertHamil(:,1:nOcc),nSites,    &
                PertHamil(:,1:nOcc),nSites,0.0_dp,PertDM,nSites)
            PertBath(:) = 0.0_dp
            PertBath(nImp+1:nSites) = PertDM(nImp+1:nSites,1)
            PertNorm = DDOT(nSites,PertBath(:),1,PertBath(:),1)
            PertBath(:) = PertBath(:) / sqrt(PertNorm)

            !Calculate derivative
            PertBath(:) = PertBath(:) - GSBath(:)
            PertBath(:) = PertBath(:) / dStep

            call writevector(PertBath,'Static Pert Bath')

            allocate(PertOrbs(nSites,nOcc))
            do i=1,nOcc
                PertOrbs(:,i) = (PertHamil(:,i) - Orbs(:,i))/dStep
            enddo

            !Matrix elements are now PertHamil(pertsite,i)
            StaticResponse = 0.0_dp
            do i=1,nOcc
                Overlap = 0.0_dp
                do j=1,nSites
                    Overlap = Overlap + PertOrbs(j,i)*Orbs(j,i)
                enddo
                write(6,*) "Overlap: ",i,Overlap
                StaticResponse = StaticResponse + Overlap*Lambda*PertOrbs(pertsite,i)*Orbs(pertsite,i)
            enddo

            write(6,*) "Mean field static resonse = ",StaticResponse,dStep

            deallocate(PertOrbs,PertHamil,PertBath,PertDM)

        enddo

        deallocate(Orbs,Energies,GSBath)

                    
    end subroutine StaticMF_DD


    !Just calculate the response at the mean-field level to compare to as U -> 0
    subroutine non_interactingLR()
        implicit none
        real(dp) :: MFDD_Response,EDiff,Omega
        integer :: n,a,pertsite

!Assume perturbation acts at site 1, and this is a local greens function. They should all be the same in the non-interacting limit
        pertsite = 1    

        Omega = Start_Omega
        do while((Omega.lt.max(Start_Omega,End_Omega)+1.0e-5_dp).and.(Omega.gt.min(Start_Omega,End_Omega)-1.0e-5_dp))

            MFDD_Response = 0.0_dp

            do n=1,nOcc
                do a=nOcc+1,nSites
                    EDiff = HFEnergies(a)-HFEnergies(n)
                    MFDD_Response = MFDD_Response + ((HFOrbs(pertsite,n)*HFOrbs(pertsite,a))**2)*Lambda*  &
                        ((1.0_dp/(EDiff-Omega))+(1.0_dp/(EDiff+Omega)))
                enddo
            enddo

            write(6,*) "Mean-field DD response: ",Omega,MFDD_Response

            Omega = Omega + Omega_Step

        enddo

    end subroutine non_interactingLR
            
!    !Calculate density density response to perturbation of frequency omega at site pertsite 
!    subroutine calc_mf_dd_response()
!        implicit none
!        integer :: i,x,a,j
!        real(dp) :: CheckOrthog,DDOT    !,StepSize
!        real(dp) , allocatable :: temp(:,:),Pert(:,:),NormB0(:)
!        character(len=*) , parameter :: t_r='calc_mf_dd_response'
!
!        if(allocated(ResponseBasis)) deallocate(ResponseBasis)
!        allocate(ResponseBasis(nSites,2))
!
!        write(6,*) "Perturbation response for orbital: ",pertsite
!        write(6,*) "Frequency of perturbation: ",Omega
!        write(6,*) "Strength of perturbation: ",Lambda
!
!        if(nImp.gt.1) call stop_all(t_r,"Response not yet coded up for > 1 impurity site")
!
!        !The response operator is (\sum_{ia} V_ai |phi_a><psi_i| + V_ia|phi_i><phi_a|) / omega - (e_a-e_i)
!        !where V_ai = <phi_a|a_pertsite^+ a_pertsite|phi_i>
!
!        !The vector corresponding to this perturbation is calculated from the impurity to the environment sites
!        !Therefore, it is (c_imp,env_a)^(1) = <orb_imp| response operator | orb_env_a>
!
!
!        ResponseBasis(:,:) = 0.0_dp !Response over impurity sites = 0
!        do x=nImp+1,nSites
!            do i=1,nOcc
!                do a=nOcc+1,nSites
!                    !This is <phi_a| a_pertsite^+ a_pertsite |phi_i> * <imp|phi_a><phi_i|x>/omega-(e_a-e_i)
!                    ResponseBasis(x,2) = ResponseBasis(x,2) + (HFOrbs(pertsite,a)*HFOrbs(pertsite,i)*HFOrbs(1,a)*HFOrbs(x,i)/ &
!                        ((HFEnergies(a)-HFEnergies(i)) - omega))
!                    ResponseBasis(x,2) = ResponseBasis(x,2) + (HFOrbs(pertsite,i)*HFOrbs(pertsite,a)*HFOrbs(1,a)*HFOrbs(x,i)/ &
!                        ((HFEnergies(a)-HFEnergies(i)) + omega))
!                enddo
!            enddo
!        enddo
!
!        !Analytically calculate new bath orbital
!        !Renormalize the change in the first order bath orbital, so that it overall noramlized (to 1st order)
!        ResponseBasis(:,2) = ResponseBasis(:,2) / sqrt(ZerothBathNorm)
!
!!        !Add the newly normalized zeroth order orbital - do we need to do that if we just want the first-order change?
!!
!!!        ResponseBasis(:,2) = ResponseBasis(:,2) + EmbeddedBasis(:,2)*   &
!!!            (1.0_dp - DDOT(nSites,EmbeddedBasis(:,2),1,ResponseBasis(:,2),1)/ZerothBathNorm)
!!
!!
!!
!!
!!       !Numerically differentiate
!!        StepSize = 0.0001
!!
!!        ResponseBasis(:,2) = ResponseBasis(:,2) * StepSize 
!!
!!        ResponseBasis(1:nSites,2) = ResponseBasis(1:nSites,2) + EmbeddedBasis(:,2)  !Add original bath orbital
!!
!!        ResponseBasis(1:nSites,2) = ResponseBasis(1:nSites,2) - EmbeddedBasis(:,2)
!!        ResponseBasis(:,2) = ResponseBasis(:,2) / StepSize
!
!        call writevector(ResponseBasis(:,2),'ResponseBasis')
!
!        CheckOrthog = DDOT(nSites,ResponseBasis(:,2),1,ResponseBasis(:,2),1)
!        write(6,*) "norm: ",CheckOrthog
!
!        !ResponseBasis is now the bath orbital for first order change in the MF solution
!        !It should be orthogonal to the original bath orbital 
!        !However, since we have got a misture of the first and second order orbital in the solution, we have to project out the first
!        !order bath orbital from the original bath
!        !B^(0)/norm[B^(0)] * (1 - <B^(0)|B^(1)>/<B^(0)|B^(0)>
!        !We have to 'unnormalize' the states
!        CheckOrthog = DDOT(nSites,EmbeddedBasis(:,2)*sqrt(ZerothBathNorm),1,ResponseBasis(:,2)*sqrt(ZerothBathNorm),1)
!        CheckOrthog = 1.0_dp - CheckOrthog/ZerothBathNorm
!        allocate(NormB0(nSites))
!        NormB0(:) = EmbeddedBasis(:,2)*CheckOrthog
!        !only *now* can we correctly check for orthogonality
!        CheckOrthog = DDOT(nSites,NormB0(:),1,ResponseBasis,1)
!        write(6,*) "Projection against other bath: ",CheckOrthog
!        deallocate(NormB0)
!
!        !Add the impurity orbital to zero. We don't want to include impurity -> impurity or impurity -> bath^(0) coupling 
!        ResponseBasis(:,1) = 0.0_dp
!!        ResponseBasis(1,1) = 1.0_dp
!
!        !Now calculate the one-electron perturbations
!        !The standard 1-electron perturbation is 1/2 Lambda a_pertsite^+ a_pertsite.
!        !We calculate this first in the HF basis, and then transform into the zeroth-order embedding basis
!        if(allocated(Pert)) deallocate(Pert)
!        allocate(Pert(nSites,nSites))
!        allocate(temp(EmbSize,nSites))
!
!        if(allocated(Emb_Pert)) deallocate(Emb_Pert)
!        allocate(Emb_Pert(EmbSize,EmbSize))
!        if(allocated(Emb_h1)) deallocate(Emb_h1)
!        allocate(Emb_h1(EmbSize,EmbSize))
!
!        Pert(:,:) = 0.0_dp
!        do i=1,nSites
!            do j=1,nSites
!                Pert(i,j) = HFOrbs(pertsite,i)*HFOrbs(pertsite,j)
!            enddo
!        enddo
!        !Transform into embedding basis
!        call DGEMM('T','N',EmbSize,nSites,nSites,1.0_dp,EmbeddedBasis,nSites,Pert,nSites,0.0_dp,temp,EmbSize)
!        call DGEMM('N','N',EmbSize,EmbSize,nSites,1.0_dp,temp,EmbSize,EmbeddedBasis,nSites,0.0_dp,Emb_Pert,EmbSize)
!
!        !Now we need to calculate H^(1)
!        !Transform h0 into the embedding basis
!        !We want C^(1)T h0 C^(0) + C^(0)T h0 C^(1)
!        call DGEMM('T','N',EmbSize,nSites,nSites,1.0_dp,ResponseBasis,nSites,h0,nSites,0.0_dp,temp,EmbSize)
!        call DGEMM('N','N',EmbSize,EmbSize,nSites,1.0_dp,temp,EmbSize,EmbeddedBasis,nSites,0.0_dp,Emb_h1,EmbSize)
!
!        call DGEMM('T','N',EmbSize,nSites,nSites,1.0_dp,EmbeddedBasis,nSites,h0,nSites,0.0_dp,temp,EmbSize)
!        call DGEMM('N','N',EmbSize,EmbSize,nSites,1.0_dp,temp,EmbSize,ResponseBasis,nSites,1.0_dp,Emb_h1,EmbSize)
!
!        call writematrix(Emb_h1,'Emb_H1',.true.)
!
!        !We now have the perturbations delta H and V in the embedding basis.
!        deallocate(temp,Pert)
!
!    end subroutine calc_mf_dd_response

end module SingRefLR
