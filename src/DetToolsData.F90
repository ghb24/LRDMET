module DetToolsData
    use const
    implicit none

    integer , allocatable :: FCIDetList(:,:)    !FCI determinant list
    integer , allocatable :: FCIBitList(:)      !Bit representations of FCI list
    integer :: nFCIDet  !Number of determinants in active space
    real(dp), allocatable :: UMat(:)
    real(dp), allocatable :: TMat(:,:)
    complex(dp), allocatable :: TMat_Comp(:,:)  !Used when asking for complex matrix elements

    !Data for the coupling of the different electron number active spaces
    integer, allocatable :: Nm1FCIDetList(:,:)  !N-1 electron determinant list
    integer, allocatable :: Nm1BitList(:)
    integer :: nNm1FCIDet   !Number of determinants in N-1 list

    integer, allocatable :: Np1FCIDetList(:,:)  !N+1 electron determinant list
    integer, allocatable :: Np1BitList(:)
    integer :: nNp1FCIDet   !Number of determinants in N+1 list

    !If tSplitMS is used, then the above lists only refer to the Ms=1/2 lists, and the following arrays correspond to the beta lists
    integer, allocatable :: Nm1bFCIDetList(:,:)  !N-1 electron determinant list
    integer, allocatable :: Nm1bBitList(:)
    integer :: nNm1bFCIDet   !Number of determinants in N-1 list

    integer, allocatable :: Np1bFCIDetList(:,:)  !N+1 electron determinant list
    integer, allocatable :: Np1bBitList(:)
    integer :: nNp1bFCIDet   !Number of determinants in N+1 list

    integer :: ECoupledSpace !Total number of determinants in the N,N-1 and N+1 spaces

    integer :: DetListStorage   !Number of integers worth of data in the storage of all det info

end module DetToolsData
