module DetToolsData
    use const
    implicit none

    integer , allocatable :: FCIDetList(:,:)    !FCI determinant list
    integer , allocatable :: FCIBitList(:)      !Bit representations of FCI list
    integer :: nFCIDet  !Number of determinants in active space
    real(dp), allocatable :: UMat(:)
    real(dp), allocatable :: TMat(:,:)

    !Data for the coupling of the different electron number active spaces
    integer, allocatable :: Nm1FCIDetList(:,:)  !N-1 electron determinant list
    integer, allocatable :: Nm1BitList(:)
    integer :: nNm1FCIDet   !Number of determinants in N-1 list

    integer, allocatable :: Np1FCIDetList(:,:)  !N+1 electron determinant list
    integer, allocatable :: Np1BitList(:)
    integer :: nNp1FCIDet   !Number of determinants in N+1 list

    integer :: ECoupledSpace !Total number of determinants in the N,N-1 and N+1 spaces

end module DetToolsData
