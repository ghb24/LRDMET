Program RealHub

use const
use timing
use errors, only: stop_all
use globals
use LinearResponse
use Solvers
use fitting
use mat_tools
implicit none

real :: start_time
logical, parameter :: tDebug = .false.

call init_calc()

call run_DMETcalc()

    contains
    
    !These want to be input parameters at some point
    subroutine set_defaults()
        implicit none

        tLR_DMET = .true. 
        tConstructFullSchmidtBasis = .true. 
        tMFResponse = .false. 
        tHalfFill = .true. 
        nSites = 12  
        LatticeDim = 1
        nImp = 2
        StartU = 1.0_dp
        EndU = 1.1_dp
        UStep = 0.2_dp
        tPeriodic = .false.
        tAntiPeriodic = .true. 
        tRampDownOcc = .true.
!        tCompleteDiag = tMFResponse
        tCompleteDiag = .true. 
        tGSFCI = .not.tCompleteDiag 
        Start_Omega = 0.0_dp
        End_Omega = 5.0_dp
        Omega_Step = 0.001_dp
        tDumpFCIDUMP = .true.
        tDiagFullSystem = .false.

    end subroutine set_defaults

    subroutine init_calc()
    use report, only: environment_report
    implicit none
!    logical :: exists

    write(6,"(A)") "***  Starting real-space hubbard calculation  ***"
    call cpu_time(start_time)
    call environment_report()

!    inquire(file='input.hub',exist=exists)
!    if(exists) then
!        call read_input()
!    else
        call set_defaults()
!    endif

    end subroutine init_calc

    subroutine run_DMETcalc()
        implicit none
        real(dp) :: ElecPerSite,FillingFrac,VarVloc,ErrRDM,mean_vloc
        integer :: i,it,Occ
        character(len=*), parameter :: t_r="run_DMETcalc"

        !Set up initial conditions, i.e. starting potential
        allocate(v_loc(nImp,nImp))
        v_loc = 0.0_dp
        allocate(h0(nSites,nSites))     !The core hamiltonian
        allocate(h0v(nSites,nSites))    !The core hamiltonian with the local potential
        allocate(HFEnergies(nSites))    !The fock energies
        allocate(MeanFieldDM(nSites,nSites))    !DM from mean-field
        HFEnergies = 0.0_dp
        h0 = 0.0_dp
        h0v = 0.0_dp

        U=StartU

        do while((U.lt.max(StartU,EndU)+1.0e-5_dp).and.(U.gt.min(StartU,EndU)-1.0e-5_dp))
            !We do a mean-field calculation on all U's on this range. (Only some later on will we do additional calculations (param sweeps)

            write(6,*) "Running DMET calculation with U = ",U

            !Calculate the core hamiltonian based on the hopping matrix of the hubbard model in real space
            call make_hop_mat()

            !Diagonalize the mean-field hamiltonian
            !Get occupations with unique GS
            call find_occs()

            do Occ=1,N_Occs

                !Loop from low occupation numbers to high
                !These occupations refer to number of closed shell orbitals, so total electrons is 2 x nOcc
                nOcc = allowed_occs(Occ)    !Number of occupied orbitals in this iteration
                NEl = 2*nOcc    !Total electrons in system
                ElecPerSite = (2.0_dp*nOcc)/real(nSites,dp) !Average number of electrons per site
                FillingFrac = ElecPerSite/2.0_dp    !A filling fraction of 1/2 is half filling (duh)

                !Write out some stats
                write(6,*) 
                write(6,"(A,F8.3,A,I5,A,I5,A)") "Electrons per site:   ",ElecPerSite," (in ", &
                    nOcc," doubly occupied orbitals on ",nSites," sites)"
                write(6,"(A,F10.5)")          "Filling Fraction:     ",FillingFrac
                write(6,"(A,F8.3)")           "Hubbard U:            ",U
                write(6,"(A,I5,A)")           "Embedded system size: ",nImp," sites"
                if(LatticeDim.eq.1) then
                write(6,"(A,I5,A)")           "1D hubbard lattice of ",nSites," sites"
                endif
                write(6,*) 
                
                if(tDiagFullSystem) then
                    !Find all eigenvalues of system
                    call DiagFullSystem()
                endif

                !Calculate full hf, including mean-field on-site repulsion (which is included in correlation potential in DMET
                call run_true_hf()

                if(tDumpFCIDUMP) then
                    call DumpFCIDUMP()
                endif

                !Calculate single reference linear response - non-interacting, TDA and RPA
                call SR_LinearResponse()

                !At this point, we have h0, U and a set of system sites (the first nImp indices), as well as a local potential
                do it=1,150

                    !Write out header file for the convergence
                    write(6,*) "Iteration: ",it

                    !Do 150 microiterations to converge the DMET for this occupation number
                    call add_localpot(h0,h0v,v_loc)

                    !Now run a HF calculation by constructing and diagonalizing the fock matrix
                    !This will also return the RDM in the AO basis
                    call run_hf(it)

                    !Construct the embedded basis
                    if(tConstructFullSchmidtBasis) then
                        call ConstructFullSchmidtBasis()
                    else
                        call CalcEmbedding()
                    endif

                    call writematrix(EmbeddedBasis,'EmbeddedBasis',.true.)
                    
                    !Now transform the 1 electron quantities into the embedded basis
                    !This should be exactly correct, i.e. we can now diagonalize the fock matrix in this basis
                    !to get the same result. We could also check that the number of electrons adds to the correct thing
                    call Transform1e()
                    
!                    if(tLR_DMET) then
!                        !Find the perturbation for the linear response in the schmidt basis of the zeroth order perturbation
!                        call FindSchmidtPert()
!                    endif

                    !Construct the two electron integrals in the system, and solve embedded system with high-level method
                    call SolveSystem(.true.)

                    if(tLR_DMET) then
                        call MR_LinearResponse()
                    endif
                    
                    call stop_all('end','end')

                    !Fit new potential
                    !vloc_change (global) is updated in here to reflect the optimal change
                    !VarVloc is a meansure of the change in the potential
                    !ErrRDM is a measure of the initial difference in the RDMs
                    call Fit_vloc(VarVloc,ErrRDM)

                    if(tDebug) call writematrix(vloc_change,'vloc_change',.true.)

                    !Mean vloc is actually for the old vloc for consistency with Geralds code
                    mean_vloc = 0.0_dp
                    do i=1,nImp
                        mean_vloc = mean_vloc + v_loc(i,i)
                    enddo
                    mean_vloc = mean_vloc/real(nImp)

                    !Write out stats:
                    !   Iter    E/Site  d[V]    ERR[RDM]    ERR[Filling]    mean[corr_pot]      Some RDM stuff...?
                    write(6,*) it,TotalE_Imp,VarVloc,ErrRDM,FillingError,mean_vloc

                    !Update vloc
                    v_loc(:,:) = v_loc(:,:) + vloc_change(:,:)

                    if(VarVloc.lt.1.0e-8) then
                        write(6,"(A)") "...correlation potential converged" 
                        exit
                    endif

                enddo

                if(it.eq.151) call stop_all(t_r,'DMET Convergence failed')

                deallocate(HFOrbs)

                !Potentially run FCI again now to get correlation functions from 2RDMs?

!                print "Eigenvalues of h0 (large system):"
!                nPrint = 5
!                iMin = max(0,nOcc-nPrint)
!                iMax = min(nOcc+nPrint,len(LatticeHfResult.EigenValues))
!                print "  iOrb:   %s" % (" ".join(["%10s"%("%i%s"%(o,[":2e","  "][o>=nOcc])) for o in range(iMin,iMax)]))
!                print "  ew(i): [%s]" % (" ".join(["%10.5f"%o for o in LatticeHfResult.EigenValues[iMin:iMax]]))
!                #print "EwSets += [(%.5f,%i,%r)]" % (U, nOcc,LatticeHfResult.EigenValues)
!                ELargePerSite = dot(LatticeHfResult.Rdm.flatten(), .5*(LatticeHfResult.Fock + h0).flatten())
!                print "Final energy per site:            %13.8f" % (HL["EnergyPerSite"])
!                print "Final 1e energy per site:         %10.5f  (%8.5f per hop)" % (HL["EImp1perSite"], nImp*HL["EImp1perSite"]/nHoppingsImp)
!                print "Final 2e energy per site:         %10.5f" % HL["EImp2perSite"]
!                print "Final imp-bath energy per site:   %10.5f  (%8.5f per hop)" % (HL["EBathImpPerSite"], nImp*HL["EBathImpPerSite"]/nHoppingsEnv)
!                print "Final filling on impurity:        %10.5f  (%+9.2f%% error)" % (HL["ActualFilling"], 100.*(HL["ActualFilling"] - Filling)/Filling)
!                print "Final U^2 * <na nb>:              %13.8f" %(U**2*HL["<na.nb>"])
!
!                print "!Diag1RDM: unique elements: ", " ".join(list(set(["%.4f" % o for o in Diag1RDM])))
!                AllResults.append( FResult(HL["EnergyPerSite"], nOcc, LatticeHfResult.Mu, LatticeHfResult.HlGap, HL["<na.nb>"], HL["EImp1perSite"], HL["EImp2perSite"], HL["EBathImpPerSite"], HL["nElec"], nImp, L,
!                    Notes = "(sym-brk: %.2e  dfil: %+.4e  dD: %.2e  dV: %.2e)" % (max(Diag1RDM)-min(Diag1RDM), HL["FillingError"], ErrRdm, VarVloc)) )
!                PrintResults(AllResults)

                write(6,"(A,F10.4,A,G20.10)") "FINAL energy per site for U=",U,' is: ',TotalE_Imp

                !Set potential for the next occupation number, or wipe it?
                v_loc(:,:) = 0.0_dp
                
            enddo   !Loop over occupations

            U=U+UStep   !Increment U

        enddo   !Loop over U values

    end subroutine run_DMETcalc
                
    !Construct full embedding basis, along with orthogonal core and virtual set, and check orthonormality of orbital space
    subroutine ConstructFullSchmidtBasis()
        implicit none
        real(dp), allocatable :: ProjOverlap(:,:),ProjOverlapEVals(:),Work(:)
        real(dp), allocatable :: RotOccOrbs(:,:),ImpurityOrbs(:,:),ProjOverlapVirt(:,:)
        real(dp), allocatable :: OverlapVirt(:,:),VirtSpace(:,:),temp(:,:) 
        real(dp) :: norm,DDOT,Overlap
        integer :: lwork,info,i,j
        character(len=*), parameter :: t_r='ConstructFullSchmidtBasis'
        integer :: nbath

        write(6,*) "Constructing full schmidt basis of HF det, including core + virt spaces"

        !HFOrbs defines HF basis.
        !Construct full projected overlap basis and diagonalize to see redundancies
        !Just do this for the occupied orbitals, and check orthogonality to the core HF and virtual orbitals
        allocate(ProjOverlap(nOcc,nOcc))
        call DGEMM('T','N',nOcc,nOcc,nImp,1.0_dp,HFOrbs(1:nImp,1:nOcc),nImp,HFOrbs(1:nImp,1:nOcc),nImp,0.0_dp,ProjOverlap,nOcc)

        !Diagonalize this
        allocate(ProjOverlapEVals(nOcc))
        ProjOverlapEVals(:)=0.0_dp
        allocate(Work(1))
        lWork=-1
        info=0
        call dsyev('V','U',nOcc,ProjOverlap,nOcc,ProjOverlapEVals,Work,lWork,info)
        if(info.ne.0) call stop_all(t_r,'Workspace queiry failed')
        lwork=int(work(1))+1
        deallocate(work)
        allocate(work(lwork))
        call dsyev('V','U',nOcc,ProjOverlap,nOcc,ProjOverlapEVals,Work,lWork,info)
        if(info.ne.0) call stop_all(t_r,'Diag failed')
        deallocate(work)

        call writevector(ProjOverlapEVals,'Projected overlap eigenvalues')

        !We should only have nImp non-zero eigenvalues
        nbath = 0
        do i=1,nOcc
            if(abs(ProjOverlapEVals(i)).gt.1.0e-7_dp) then
                nbath = nbath + 1
            endif
        enddo
        if(nbath.gt.nImp) call stop_all(t_r,'error here')

        !Now rotate original occupied orbitals into this new orthogonal basis
        allocate(RotOccOrbs(nSites,nOcc))
        call DGEMM('N','N',nSites,nOcc,nOcc,1.0_dp,HFOrbs(1:nSites,1:nOcc),nSites,ProjOverlap,nOcc,0.0_dp,RotOccOrbs,nSites)

        !RotOccOrbs now represents the rotated occupied orbitals into the bath basis. 
        !Only the last nImp orbitals will have any coupling to the impurity on them.
        !These RotOccOrbs constitute a legitamate HF wavefunction, are orthonormal to all other orbitals. Just simple rotation.
!        call writematrix(RotOccOrbs,'Occupied Orbs in schmidt basis',.true.)

        !Construct bath states by projecting out component on impurity and renormalizing
        !Assume last nImp states are the states with overlap with impurity only
        !Also normalize these orbitals
        do i=nOcc,nOcc-nImp+1,-1
            RotOccOrbs(1:nImp,i) = 0.0_dp
            norm = DDOT(nSites,RotOccOrbs(:,i),1,RotOccOrbs(:,i),1)
            norm = sqrt(norm)
            RotOccOrbs(:,i) = RotOccOrbs(:,i)/norm
        enddo

        !These states are now the bath states.
!        call writematrix(RotOccOrbs(:,nOcc-nImp+1:nOcc),'Bath orbitals',.true.)

        allocate(ImpurityOrbs(nSites,nImp))
        ImpurityOrbs(:,:) = 0.0_dp
        do i=1,nImp
            ImpurityOrbs(i,i) = 1.0_dp
        enddo

        !We now have all the orbitals. Which are orthogonal to which?
        write(6,*) "All impurity/bath orbitals orthogonal by construction"
        do i=nOcc,nOcc-nImp+1,-1
            do j=nOcc,nOcc-nImp+1,-1
                Overlap = DDOT(nSites,RotOccOrbs(:,i),1,RotOccOrbs(:,j),1)
                if(i.eq.j) then
                    if(abs(Overlap-1.0_dp).gt.1.0e-7_dp) then
                        call stop_all(t_r,'bath orbitals not normalized set')
                    endif
                else
                    if(abs(Overlap).gt.1.0e-7_dp) then
                        call stop_all(t_r,'bath orbitals not orthogonal set')
                    endif
                endif
            enddo
        enddo
        write(6,*) "All bath/bath orbitals orthonormal"

        !Now, consider the rotated core orbitals
        !Are they orthonormal to the bath orbitals (they are othogonal to impurity, since they have no component on them)
        !They must also be orthogonal to the bath orbitals, since they never had any component on the impurity, 
        !which is the only bit which has changed in constructing the bath (norm doesn't affect)
        do i=1,nOcc-nImp
            do j=nOcc,nOcc-nImp+1,-1
                !Overlap of core (i) with bath (j)
                Overlap = DDOT(nSites,RotOccOrbs(:,i),1,RotOccOrbs(:,j),1)
                if(abs(Overlap).gt.1.0e-7_dp) then
                    call stop_all(t_r,'bath orbitals with core not orthogonal set')
                endif
            enddo
        enddo

        !However, the virtual space is *not* orthogonal to the embedded system (though it is wrt the core).
        !Now create orthogonal set of orbitals from the virtual space. There will now be a redundancy.
        !Calculate the overlap of the virtual space with a projection onto the embedded system.
        !From the diagonalization of this, we expect exactly *nImp* non-zero eigenvalues, which are the redundant orbitals.
        !Remove these, and the rest are the now non-canonical virtual orbital space.
        allocate(ProjOverlapVirt(nOcc,nOcc))

        !This array is used to calculate the overlap of each virtual space function with each impurity
        !function
        allocate(OverlapVirt(2*nImp,nOcc))
        !The first nImp correspond to the impurity orbital, and the next two correspond to the bath orbitals
        OverlapVirt(:,:) = 0.0_dp
        do i=nOcc+1,nSites  !run through virtual space
            do j=1,nImp     !run through impurity orbitals
                OverlapVirt(j,i-nOcc) = HFOrbs(j,i)
            enddo
        enddo

        !Now calculate overlap with bath orbitals
        call DGEMM('T','N',nImp,nOcc,nSites,1.0_dp,RotOccOrbs(:,nOcc-nImp+1:nOcc),nSites,HFOrbs(:,nOcc+1:nOcc), &
            nSites,0.0_dp,OverlapVirt(nImp+1:2*nImp,1:nOcc),2*nImp-nImp)

        !Combine overlaps to get full projected overlap matrix
        call DGEMM('T','N',nOcc,nOcc,2*nImp,1.0_dp,OverlapVirt,2*nImp,OverlapVirt,2*nImp,0.0_dp,ProjOverlapVirt,nOcc)

        !Diagonalize this
        deallocate(ProjOverlapEVals)
        allocate(ProjOverlapEVals(nOcc))
        ProjOverlapEVals(:)=0.0_dp
        allocate(Work(1))
        lWork=-1
        info=0
        call dsyev('V','U',nOcc,ProjOverlapVirt,nOcc,ProjOverlapEVals,Work,lWork,info)
        if(info.ne.0) call stop_all(t_r,'Workspace queiry failed')
        lwork=int(work(1))+1
        deallocate(work)
        allocate(work(lwork))
        call dsyev('V','U',nOcc,ProjOverlapVirt,nOcc,ProjOverlapEVals,Work,lWork,info)
        if(info.ne.0) call stop_all(t_r,'Diag failed')
        deallocate(work)

!        call writevector(ProjOverlapEVals,'virtual overlap eigenvalues')

        nbath = 0   !Count the number of virtual functions which span the space of the embedded system
        do i=1,nOcc
            if(abs(ProjOverlapEVals(i)).gt.1.0e-7_dp) then
                nbath = nbath + 1
            endif
        enddo
        if(nbath.ne.nImp) then
            call stop_all(t_r,'Virtual space redundancy not as expected')
        endif

        !Now rotate orbitals such that they are orthogonal, while deleting the redundant ones
        !Assume the last nImp are redundant
        allocate(VirtSpace(nSites,nOcc-nImp))
        call DGEMM('N','N',nSites,nOcc-nImp,nOcc,1.0_dp,HFOrbs(:,nOcc+1:nSites),nSites,ProjOverlapVirt(1:nOcc,1:nOcc-nImp),nOcc, &
            0.0_dp,VirtSpace,nSites)

        !Check now that the virtual space is now orthogonal to all occupied HF orbitals, as well as the impurity and bath sites

        !First against the impurity sites
        do i=1,nImp
            do j=1,nOcc-nImp
                Overlap = DDOT(nSites,ImpurityOrbs(:,i),1,VirtSpace(:,j),1)
                if(abs(Overlap).gt.1.0e-7_dp) then
                    call stop_all(t_r,'virtual orbitals not orthogonal to impurity orbitals')
                endif
            enddo
        enddo

        do i=1,nOcc
            do j=1,nOcc-nImp
                Overlap = DDOT(nSites,RotOccOrbs(:,i),1,VirtSpace(:,j),1)
                if(abs(Overlap).gt.1.0e-7_dp) then
                    if(i.gt.(nOcc-nImp)) then
                        call stop_all(t_r,'virtual orbitals not orthogonal to bath orbitals')
                    else
                        call stop_all(t_r,'virtual orbitals not orthogonal to core orbitals')
                    endif
                endif
            enddo
        enddo

        if(allocated(FullSchmidtBasis)) deallocate(FullSchmidtBasis)
        allocate(FullSchmidtBasis(nSites,nSites))   ! (Atomicbasis,SchmidtBasis)
        FullSchmidtBasis(:,:) = 0.0_dp
        FullSchmidtBasis(:,1:nOcc-nImp) = RotOccOrbs(:,1:nOcc-nImp)    !The first orbitals consist of core  
        FullSchmidtBasis(:,nOcc-nImp+1:nOcc) = ImpurityOrbs(:,:)    !Impurity Orbs
        FullSchmidtBasis(:,nOcc+1:nOcc+nImp) = RotOccOrbs(:,nOcc-nImp+1:nOcc)   !Bath
        FullSchmidtBasis(:,nOcc+nImp+1:nSites) = VirtSpace(:,:)     !Virtual space

        !call writematrix(FullSchmidtBasis,'FullSchmidtBasis',.true.)

        !Construct unitary basis transformation matrix from HF to embedded basis
        if(allocated(HFtoSchmidtTransform)) deallocate(HFtoSchmidtTransform)
        allocate(HFtoSchmidtTransform(nSites,nSites))
        call DGEMM('T','N',nSites,nSites,nSites,1.0_dp,HFOrbs,nSites,FullSchmidtBasis,nSites,0.0_dp,HFtoSchmidtTransform,nSites)

        allocate(temp(nSites,nSites))
        !Check that this operator is unitary
        !call DGEMM('N','T',nSites,nSites,nSites,1.0_dp,HFtoSchmidtTransform,nSites,HFtoSchmidtTransform,nSites,0.0_dp,temp,nSites)
        !call writematrix(temp,'Test of unitarity of HF to Schmidt Transform',.true.)
        !do i=1,nSites
        !    do j=1,nSites
        !        if((i.ne.j).and.(abs(temp(i,j)).gt.1.0e-7_dp)) then
        !            call stop_all(t_r,'Transformation matrix not unitary')
        !        elseif((i.eq.j).and.(abs(temp(i,j)-1.0_dp).gt.1.0e-7)) then
        !            call stop_all(t_r,'Transformation matrix not unitary')
        !        endif
        !    enddo
        !enddo

        !Use this to rotate the fock operator into this new basis
        if(allocated(FockSchmidt)) deallocate(FockSchmidt)
        allocate(FockSchmidt(nSites,nSites))
        !Set up FockSchmidt to temperarily to be the HF basis fock operator (i.e. diagonal)
        FockSchmidt(:,:) = 0.0_dp
        do i=1,nSites
            FockSchmidt(i,i) = HFEnergies(i)
        enddo
        call DGEMM('T','N',nSites,nSites,nSites,1.0_dp,HFtoSchmidtTransform,nSites,FockSchmidt,nSites,0.0_dp,temp,nSites)
        call DGEMM('N','N',nSites,nSites,nSites,1.0_dp,temp,nSites,HFtoSchmidtTransform,nSites,0.0_dp,FockSchmidt,nSites)
        deallocate(temp)
        !FockSchmidt has been overwritten with the fock matrix in the schmidt basis
!        call writematrix(FockSchmidt,'Fock in schmidt basis',.true.)

        nSys = nImp !Fix this here
        
        if(allocated(EmbeddedBasis)) deallocate(EmbeddedBasis)
        allocate(EmbeddedBasis(nSites,nImp+nSys))
        EmbeddedBasis(:,:) = 0.0_dp
        EmbeddedBasis(:,1:nImp) = ImpurityOrbs(:,:)
        EmbeddedBasis(:,nImp+1:2*nImp) = FullSchmidtBasis(:,nOcc+1:nOcc+nImp)

        EmbSize = nImp+nSys !This is the total size of the embedded system with which to do the high-level calculation on 
        
        !Calculate some paramters which will be used later, which define the size of triangular packed arrays over the impurity sites, or
        !the entire embedding sites.
        nImpCombs = (nImp*(nImp+1))/2
        EmbCombs = (EmbSize*(EmbSize+1))/2

        deallocate(ProjOverlap,ProjOverlapEVals,RotOccOrbs,ImpurityOrbs,ProjOverlapVirt,OverlapVirt,VirtSpace)

    end subroutine ConstructFullSchmidtBasis


    !Transform the one-electron integrals into the embedded basis by unitary transformation
    subroutine Transform1e()
        implicit none
        real(dp), allocatable :: CorrPotential(:,:),temp(:,:)
!        real(dp), allocatable :: FockPotential(:,:)
        integer :: i,j
!        character(len=*), parameter :: t_r="Transform1e"
        
        if(allocated(Emb_Fock)) deallocate(Emb_Fock)
        if(allocated(Emb_h0)) deallocate(Emb_h0)                 !Core hamiltonian
        if(allocated(Emb_MF_DM)) deallocate(Emb_MF_DM)        !Mean field RDM
        if(allocated(Emb_FockPot)) deallocate(Emb_FockPot)      !The fock potential transforms h0v into the fock matrix in the AO basis
        if(allocated(Emb_CorrPot)) deallocate(Emb_CorrPot)      !The correlation potential is the block diagonal v_pot in the AO basis

        !Transform all of these quantities into the embedded basis
        allocate(Emb_h0(EmbSize,EmbSize))
        allocate(Emb_MF_DM(EmbSize,EmbSize))
        allocate(Emb_FockPot(EmbSize,EmbSize))
        allocate(Emb_CorrPot(EmbSize,EmbSize))
        allocate(Emb_Fock(EmbSize,EmbSize))     !For hub, this is just the core + corrPot

        !Rotate them
        allocate(temp(EmbSize,nSites))
        !Core hamiltonian

!        call writematrix(EmbeddedBasis,'EmbBasis',.true.)

        call DGEMM('T','N',EmbSize,nSites,nSites,1.0_dp,EmbeddedBasis,nSites,h0,nSites,0.0_dp,temp,EmbSize)
        call DGEMM('N','N',EmbSize,EmbSize,nSites,1.0_dp,temp,EmbSize,EmbeddedBasis,nSites,0.0_dp,Emb_h0,EmbSize)

        if(tDebug) call writematrix(Emb_h0,'Embedded Core hamil',.true.)

        !Mean field RDM
        call DGEMM('T','N',EmbSize,nSites,nSites,1.0_dp,EmbeddedBasis,nSites,MeanFieldDM,nSites,0.0_dp,temp,EmbSize)
        call DGEMM('N','N',EmbSize,EmbSize,nSites,1.0_dp,temp,EmbSize,EmbeddedBasis,nSites,0.0_dp,Emb_MF_DM,EmbSize)
        
        if(tDebug) call writematrix(Emb_MF_DM,'Embedded RDM',.true.)

        !Since we don't actually store the fock potential (the diagonal repulsion/2 electron terms), construct it here in the AO basis
        !This is the bit of the fock matrix not in the core hamiltonian
!        allocate(FockPotential(nSites,nSites))
!        FockPotential(:,:) = 0.0_dp
!        do i=1,nSites
!            !Include the on-site repulsion
!            FockPotential(i,i) = FockPotential(i,i) + U * 0.5_dp * (NEl/real(nSites))
!        enddo
!        call DGEMM('T','N',EmbSize,nSites,nSites,1.0_dp,EmbeddedBasis,nSites,FockPotential,nSites,0.0_dp,temp,EmbSize)
!        call DGEMM('N','N',EmbSize,EmbSize,nSites,1.0_dp,temp,EmbSize,EmbeddedBasis,nSites,0.0_dp,Emb_FockPot,EmbSize)
!        deallocate(FockPotential)
        Emb_FockPot(:,:)=0.0_dp     !We don't include the 2e terms in the fock matrix here, instead, leaving them to be captured by the correlation potential
!        call writematrix(Emb_FockPot,"Embedded Fock potential",.true.)

        !We also do not store the "Correlation potential", which is the potential which is added to the fock matrix to make the DMs match
        allocate(CorrPotential(nSites,nSites))
        CorrPotential(:,:) = 0.0_dp
        !It is just h0v - h0
        do i=1,nSites
            do j=1,nSites
                CorrPotential(i,j) = h0v(i,j) - h0(i,j)
            enddo
        enddo
        call DGEMM('T','N',EmbSize,nSites,nSites,1.0_dp,EmbeddedBasis,nSites,CorrPotential,nSites,0.0_dp,temp,EmbSize)
        call DGEMM('N','N',EmbSize,EmbSize,nSites,1.0_dp,temp,EmbSize,EmbeddedBasis,nSites,0.0_dp,Emb_CorrPot,EmbSize)
        deallocate(CorrPotential)
        deallocate(temp)

        if(tDebug) call writematrix(Emb_CorrPot,"Embedded correlation potential",.true.)

        !The 2e terms left out of fock matrix everywhere, although included basically in CorrPot
        Emb_Fock(:,:) = Emb_h0(:,:) + Emb_CorrPot(:,:) ! + Emb_FockPot  
        
    end subroutine Transform1e

    !Calculate the embedded basis, and transform the 1e operators into this new basis
    subroutine CalcEmbedding()
        implicit none
        real(dp) :: ImpurityOverlap(nImp,nImp),OverlapEVs(nImp),DDOT
        real(dp), allocatable :: RDMonImp(:,:),Work(:),SminHalf(:,:),temp(:,:)
        integer :: info,lWork,i,j,nDelete
        character(len=*), parameter :: t_r="CalcEmbedding"

        !Take density matrix over impurity: nSites-nImp x nImp
        allocate(RDMonImp(nSites-nImp,nImp))
        do i=nImp+1,nSites
            do j=1,nImp
                RDMonImp(i-nImp,j) = MeanFieldDM(i,j)
            enddo
        enddo

        !For response normalizationi
        ZerothBathNorm = 0.0_dp
        ZerothBathNorm = DDOT(nSites-nImp,RDMonImp(1:nSites-nImp,1),1,RDMonImp(1:nSites-nImp,1),1)

!        call writematrix(RDMonImp,'RDMonImp',.true.)

        !Now, we want to Lowdin orthogonalize these orbitals, i.e. |psi_alpha> = |psi_beta> S_beta,alpha^(-1/2)

        !Find overlap of this density over the impurity (contract out nSites-nImp)
        call DGEMM('T','N',nImp,nImp,nSites-nImp,1.0_dp,RDMonImp,nSites-nImp,RDMonImp,nSites-nImp,0.0_dp,ImpurityOverlap,nImp)
!        call writematrix(ImpurityOverlap,'s',.true.)
        !Now find the eigenbasis of this overlap in order to raise it to the power of -1/2
        !Diagonalize the system over the impurity sites
        allocate(Work(1))
        lWork=-1
        info=0
        call dsyev('V','U',nImp,ImpurityOverlap,nImp,OverlapEVs,Work,lWork,info)
        if(info.ne.0) call stop_all(t_r,'Workspace queiry failed')
        lwork=int(work(1))+1
        deallocate(work)
        allocate(work(lwork))
        call dsyev('V','U',nImp,ImpurityOverlap,nImp,OverlapEVs,Work,lWork,info)
        if(info.ne.0) call stop_all(t_r,'Diag failed')
        deallocate(work)
        if(nImp.eq.1) then
            if(abs(ImpurityOverlap(1,1)-1.0_dp).gt.1.0e-7_dp) call stop_all(t_r,'Diag error')
        endif
        !call writevector(OverlapEVs,'S_EVs')

        nSys = 0    !nSys is the number of orbitals in the bath having removed linear dependencies
        nDelete = 0
        do i=1,nImp
            if(OverlapEVs(i).lt.1.0e-10_dp) then
                write(6,*) "Warning: Bath basis linearly dependent"
                write(6,*) "Overlap eigenvalue: ",OverlapEVs(i)
                nDelete = nDelete + 1
            else
                nSys = nSys + 1
            endif
        enddo

        if(nSys.ne.nImp) write(6,*) "Bath orbitals removed due to linear dependency: ",nDelete
        write(6,*) "Total bath orbitals: ",nSys
        if(nDelete.ne.(nImp-nSys)) call stop_all(t_r,'nDelete.ne.(nImp-nSys)')

        !Now, we need to construct the full embedded basis including orbitals on the impurity sites
        !The orbitals to remove if necessary will be first.
        !Construct S^(-1/2), and then rotate back into the original basis
        allocate(SminHalf(nSys,nSys))
        SminHalf = 0.0_dp
        do i=1,nSys
            SminHalf(i,i) = OverlapEVs(i+nDelete)**(-0.5_dp)
        enddo
!        call writematrix(ImpurityOverlap(1:nImp,nDelete+1:nImp),'ImpOverlap',.true.)
        !Now rotate back into original basis as U D U^T (Having removed the first nDelete orbitals from the Impurity overlap eigenvectors
        allocate(temp(nImp,nSys))
        call DGEMM('N','N',nImp,nSys,nSys,1.0_dp,ImpurityOverlap(1:nImp,nDelete+1:nImp),nImp,SminHalf,nSys,0.0_dp,temp,nImp)
        call DGEMM('N','T',nImp,nImp,nSys,1.0_dp,temp,nImp,ImpurityOverlap(1:nImp,nDelete+1:nImp),nImp,0.0_dp,SminHalf,nImp)
        deallocate(temp)
        !SminHalf in now in the original basis
!        call writematrix(SminHalf,'SminHalf',.true.)

        !Now, we need to multiply the original bath orbitals extracted from the RDM on the impurity by this matrix to complete the Lowdin orthogonalization
        allocate(temp(nSites-nImp,nImp))
        call DGEMM('N','N',nSites-nImp,nImp,nImp,1.0_dp,RDMonImp,nSites-nImp,SminHalf,nImp,0.0_dp,temp,nSites-nImp)
        !temp now holds the orthogonalized bath orbitals
!        call writematrix(temp,'Orthog bath',.true.)

        !Now include the impurity sites along with the bath orbitals to make the final embedded basis
        !The system orbitals are represented by the first nImp orbitals, and the bath orbitals are the last nSys
        if(allocated(EmbeddedBasis)) deallocate(EmbeddedBasis)
        allocate(EmbeddedBasis(nSites,nImp+nSys))
        EmbeddedBasis(:,:) = 0.0_dp
        do i=1,nImp
            EmbeddedBasis(i,i) = 1.0_dp
        enddo

        if(nImp.ne.nSys) call stop_all(t_r,'You need to check in the code here to make sure we have the right bath orbitals')
        !Here, we only want to add orbitals on the bath which haven't been removed due to being linearly dependent
        !Which orbitals in temp are the removed linear dependent ones? TODO: Check this is being done correctly
        !Currently, we are just taking the first nSys of them
        EmbSize = nImp+nSys !This is the total size of the embedded system with which to do the high-level calculation on 
        EmbeddedBasis(nImp+1:nSites,nImp+1:EmbSize) = temp(:,1:nSys)
        
        !Calculate some paramters which will be used later, which define the size of triangular packed arrays over the impurity sites, or
        !the entire embedding sites.
        nImpCombs = (nImp*(nImp+1))/2
        EmbCombs = (EmbSize*(EmbSize+1))/2

        deallocate(temp,SminHalf)

    end subroutine CalcEmbedding

    subroutine DumpFCIDUMP()
        use utils, only: get_free_unit
        implicit none
        integer :: iunit,i,j,k,l,A,B,ex(2,2)
        real(dp) :: hel,GetHFInt_spinorb
        real(dp), allocatable :: temp(:,:),h0HF(:,:)
        
        !Calculate hopping matrix in MO basis
        allocate(temp(nSites,nSites))
        allocate(h0HF(nSites,nSites))
        call dgemm('t','n',nSites,nSites,nSites,1.0_dp,FullHFOrbs,nSites,h0,nSites,0.0_dp,temp,nSites)
        call dgemm('n','n',nSites,nSites,nSites,1.0_dp,temp,nSites,FullHFOrbs,nSites,0.0_dp,h0HF,nSites)
        deallocate(temp)

        iunit = get_free_unit()
        open(unit=iunit,file='FCIDUMP',status='unknown')
        write(iunit,'(2A6,I3,A7,I3,A5,I2,A)') '&FCI ','NORB=',nSites,'NELEC=',NEl,',MS2=',0,','
        WRITE(iunit,'(A9)',advance='no') 'ORBSYM='
        do i=1,nSites
            write(iunit,'(I1,A1)',advance='no') 1,','
        enddo
        write(iunit,*) ""
        WRITE(iunit,'(A7,I1)') 'ISYM=',1
        WRITE(iunit,'(A5)') '&END'
        do i=1,nSites
            do j=1,nSites
                A=(i*(i-1))/2+j
                DO k=1,nSites
                    DO l=1,nSites
                        B=(k*(k-1))/2+l

                        !IF(B.lt.A) CYCLE
                        !IF((i.lt.j).and.(k.lt.l)) CYCLE
                        !IF((i.gt.j).and.(k.lt.l)) CYCLE

                        ex(1,1) = 2*i
                        ex(1,2) = 2*k
                        ex(2,1) = 2*j
                        ex(2,2) = 2*l
                        hel = GetHFInt_spinorb(ex,FullHFOrbs)
                        if(abs(hel).gt.1.0e-8_dp) then
                            WRITE(iunit,'(1X,G20.14,4I3)') hel,i,j,k,l
                        endif
                    enddo
                enddo
            enddo
        enddo
        do i=1,nSites
            do j=1,i
                if(abs(h0HF(i,j)).gt.1.0e-8_dp) then
                    WRITE(iunit,'(1X,G20.14,4I3)') h0HF(i,j),i,j,0,0
                endif
            enddo
        enddo
        do i=1,nSites
            WRITE(iunit,'(1X,G20.14,4I3)') FullHFEnergies(i),i,0,0,0
        enddo
        WRITE(iunit,'(1X,G20.14,4I3)') 0.0_dp,0,0,0,0
        close(iunit)
        deallocate(h0HF)

        !Now dump in the AO basis
        open(unit=iunit,file='AO_FCIDUMP',status='unknown')
        write(iunit,'(2A6,I3,A7,I3,A5,I2,A)') '&FCI ','NORB=',nSites,'NELEC=',NEl,',MS2=',0,','
        WRITE(iunit,'(A9)',advance='no') 'ORBSYM='
        do i=1,nSites
            write(iunit,'(I1,A1)',advance='no') 1,','
        enddo
        write(iunit,*) ""
        WRITE(iunit,'(A7,I1)') 'ISYM=',1
        WRITE(iunit,'(A5)') '&END'
        do i=1,nSites
            write(iunit,'(1X,G20.14,4I3)') U,i,i,i,i
        enddo
        do i=1,nSites
            do j=1,i
                if(abs(h0(i,j)).gt.1.0e-8_dp) then
                    WRITE(iunit,'(1X,G20.14,4I3)') h0(i,j),i,j,0,0
                endif
            enddo
        enddo
        WRITE(iunit,'(1X,G20.14,4I3)') 0.0_dp,0,0,0,0
        close(iunit)

    end subroutine DumpFCIDUMP

End Program RealHub
