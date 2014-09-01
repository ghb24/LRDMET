module T1SCF
    use const
    implicit none

    contains

    !TODO: Check all feeds through correctly and that there are no more diagonalizations
    !of the lattice

    !Rotate the previous iterations Slater determinant into the new Slater determinant
    !over the entire lattice..
    !This is achieved by rotating the imp blocks of the T1 matrix into the site basis, 
    !before applying the rotation to all sites.
    subroutine RotDet_T1()
        implicit none


    end subroutine RotDet_T1

    !Find the T1 amplitudes of exp(T) which will rotate a HF wavefunction over the embedded
    !space to a new determinant which maximally overlaps with the HL wavefunction.
    !This will involve a HF calculation first (which will give back the entangled orbs?)
    !The T1 amplitudes of the rotation are then constructed from < HL | a^+ i | HF >.
    !Store these T1 amplitudes over the impurity + bath space in a global structure.
    subroutine GetT1Overlap()
        implicit none
        
        !Diagonalize embedded one-electron hamiltonian
        if(.not.allocated(Emb_h0)) call stop_all(t_r,'Emb_h0 not allocated')
        allocate(Emb_MFOrbs(EmbSize,EmbSize))
        allocate(EmbVals(EmbSize))
        Emb_MFOrbs(:,:) = Emb_h0(:,:)

        call DiagOneEOp(Emb_MFOrbs,EmbVals,nImp,EmbSize,.false.,.true.)

        !We now have a SD over the embedded space which corresponds to the equivalent
        !det over the lattice
        !We now need to find the nImp^2 T1 amplitudes
        !First, transform the HL wavefunction into the fock basis?
        !


    end subroutine GetT1Overlap

end module T1SCF
