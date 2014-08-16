module LatticesData
    use const
    implicit none
    save

    real(dp), allocatable :: LatticeVector(:,:)
    integer :: CellShape    !Flag for the type of lattice. 1=Orig tilted Lattice, 2 = Square cell, 3=Tilted
    logical :: tSquareLatt
    integer :: iImpRepeats  !Cells in the Supercell, or the number of repeats of the impurites striped through the system
    integer, allocatable :: ImpSites(:)             !The list of site indices for the impurity sites
    integer, allocatable :: StripedImpIndices(:,:)  !The index of the impurity repeats (nImp,iImpRepeats)

end module
