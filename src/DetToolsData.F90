module DetToolsData
    use const
    implicit none

    integer , allocatable :: FCIDetList(:,:)    !FCI determinant list
    integer :: nFCIDet  !Number of determinants in active space
    real(dp), allocatable :: UMat(:)
    real(dp), allocatable :: TMat(:,:)

end module DetToolsData
