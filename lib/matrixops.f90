MODULE matrixops

  IMPLICIT NONE

!  INTERFACE det2
!    MODULE PROCEDURE r_det
!    MODULE PROCEDURE c_det
!  END INTERFACE

  INTERFACE sum2
    MODULE PROCEDURE d_sum2
  END INTERFACE

  INTERFACE mat_inv
    MODULE PROCEDURE d_inv
    MODULE PROCEDURE z_inv2
    !MODULE PROCEDURE z_inv
  END INTERFACE

  INTERFACE det
!    MODULE PROCEDURE s_det2
    MODULE PROCEDURE d_det2
!    MODULE PROCEDURE c_det2
    MODULE PROCEDURE z_det2
    ! Using LAPACK's LU
    !       sgetrf            Computes an LU factorization of a
    !       dgetrf            general matrix, using partial pivoting
    !       cgetrf            with row interchanges.
    !       zgetrf
  END INTERFACE
  INTERFACE write_mat
    MODULE PROCEDURE r_write_mat1
    MODULE PROCEDURE r_write_mat
    MODULE PROCEDURE c_write_mat
    MODULE PROCEDURE i_write_mat
    MODULE PROCEDURE l_write_mat
  END INTERFACE

CONTAINS

  SUBROUTINE z_inv2(mat,matinv,dimen)
    use const
    integer, intent(in) :: dimen
    complex(dp), intent(in) :: mat(dimen,dimen)
    complex(dp), dimension(dimen,dimen), intent(out) :: matinv
    integer, dimension(dimen) :: ipiv
    integer :: msize,nsize,lwork,info
    complex(dp), allocatable :: cWork(:)

    msize=dimen
    nsize=dimen
    if(msize.ne.nsize) stop 'error in z_inv2'
    matinv(:,:) = mat(:,:)

    info=0
    call ZGETRF(msize,nsize,matinv,nsize,ipiv,info)
    IF (INFO /= 0) then
        write(6,*) "info: ",info
        STOP 'Error with z_inv2 1'
    endif
    allocate(cWork(1))
    lwork = -1
    call ZGETRI(msize,matinv,msize,ipiv,cwork,lwork,info)
    if(info.ne.0) stop 'error with workspace query in z_inv2'
    lWork = int(real(cWork(1))) + 1
    deallocate(cWork)
    allocate(cWork(lWork))
    call ZGETRI(msize,matinv,msize,ipiv,cwork,lwork,info)
    if(info.ne.0) stop 'error with inversion in z_inv2'
    deallocate(cWork)

  end subroutine z_inv2




  SUBROUTINE z_inv (mat,matinv)
    use const
    
    CHARACTER,PARAMETER :: trans='N'
    COMPLEX(dp), INTENT(IN) :: mat(:,:)
!    COMPLEX(KIND=KIND(1.d0)), DIMENSION(SIZE(mat(:,1)),SIZE(mat(1,:))) ::&
    COMPLEX(dp), DIMENSION(SIZE(mat,1),SIZE(mat,2)) ::&
        & matinv, matdum
    INTEGER, DIMENSION(SIZE(mat,1)) :: IPIV
    INTEGER :: nsize,msize,i,INFO

    msize=SIZE(mat,1)
    nsize=SIZE(mat,2)
    matdum=mat
    matinv=zzero
    DO i=1,msize
       matinv(i,i)=zone
    END DO
    INFO=0 
    !write(6,*) matdum,msize,nsize,info
    CALL zGETRF(msize,nsize,matdum,nsize,IPIV,INFO)
    IF (INFO /= 0) then
        write(6,*) "info: ",info
        STOP 'Error with z_inv 1'
    endif
    INFO=0 
    CALL zGETRS(trans,msize,nsize,matdum,nsize,IPIV,matinv,msize,INFO)
    IF (INFO /= 0) STOP 'Error with z_inv 2'
    
  END SUBROUTINE z_inv

  SUBROUTINE d_inv (mat,matinv,dimen)
    
    use const
    implicit none
    integer, intent(in) :: dimen
    CHARACTER,PARAMETER :: trans='N'
    REAL(dp), INTENT(IN) :: mat(dimen,dimen)
!    REAL(KIND=KIND(1.d0)), DIMENSION(SIZE(mat(:,1)),SIZE(mat(1,:))) ::&
    REAL(dp), DIMENSION(dimen,dimen) ::&
        & matinv, matdum
    INTEGER, DIMENSION(dimen) :: IPIV
    INTEGER :: nsize,msize,i,INFO

    msize=dimen
    nsize=dimen
    matdum=mat
    matinv=zero
    DO i=1,msize
       matinv(i,i)=one
    END DO
    INFO=-1
    CALL dGETRF(msize,nsize,matdum,nsize,IPIV,INFO)
!    IF (INFO /= 0) STOP 'Error with d_inv 1'
    IF (INFO /= 0) THEN
      WRITE(*,*) 'Warning from d_inv 1', INFO
      STOP
    END IF
    INFO=-1
    CALL dGETRS(trans,msize,nsize,matdum,nsize,IPIV,matinv,msize,INFO)
    IF (INFO /= 0) STOP 'Error with d_inv 2'
    
  END SUBROUTINE d_inv

  FUNCTION d_sum2(A) RESULT(sum2)

    REAL(KIND=KIND(0.d0)) :: sum2, A(:)
    REAL(KIND=KIND(0.d0)), DIMENSION(SIZE(a(:))) :: Orig
    INTEGER :: N
!    INTEGER, DIMENSION(SIZE(A(:))) :: indx
    REAL(KIND=KIND(0.d0)) :: dnrm2

    N=SIZE(A(:))

    IF (N==0) THEN
      sum2=1.d0
    ELSE
      Orig=A
      sum2=dnrm2(N,Orig,1)
    END IF
  END FUNCTION d_sum2


  FUNCTION z_det2(A) RESULT(det)

    IMPLICIT NONE
    COMPLEX(KIND=KIND(0.d0)) :: det, A(:,:)
    COMPLEX(KIND=KIND(0.d0)), DIMENSION(SIZE(A,1),SIZE(A,2)) :: B
    INTEGER :: j, N, err, p
    INTEGER, DIMENSION(SIZE(A,1)) :: indx
    N=SIZE(A(:,1))

    IF (N==0) THEN
      det=(1.d0,0.d0)
    ELSE
      !CALL zgemt('N',N,N,(1.d0,0.d0),A,N,B,N)
      B=A
      CALL zgetrf(N,N,B,N,indx,err)
!      IF (err /= 0) STOP 'Error with det z_det2'
      IF (err /= 0) THEN 
         write(*,*) 'det->0'
         det=0.d0
      ENDIF
      p=0
      det=(1.d0,0.d0)
      DO j=1, N
    det=det*B(j,j)/SQRT(REAL(j,KIND=KIND(1.d0)))
    IF (indx(j)/=j) THEN
!	  p=p+1
      det=-det
    END IF
      END DO
!      det=det*(-1)**p
    END IF
  END FUNCTION z_det2
  FUNCTION d_det2(A) RESULT(det)

    IMPLICIT NONE
    REAL(KIND=KIND(0.d0)) :: det, A(:,:)
    REAL(KIND=KIND(0.d0)), DIMENSION(SIZE(A,1),SIZE(a,2)) :: B
    INTEGER :: j, N, err, p
    INTEGER, DIMENSION(SIZE(A,1)) :: indx
    N=SIZE(A(:,1))

    IF (N==0) THEN
      det=1.d0
    ELSE
      B=A
      CALL dgetrf(N,N,B,N,indx,err)
      IF (err /= 0) THEN 
    WRITE(*,*) 'Number ', err
    STOP 'Error with d_det2'
      END IF
!	  
!THEN
!STOP 'Not tested yet'
!     ELSE
!
!      END IF
      p=0
      det=1.d0
      DO j=1, N
    det=det*B(j,j)/SQRT(REAL(j,KIND=KIND(1.d0)))
    IF (indx(j)/=j) THEN
      p=p+1
!	  det=-det
    END IF
      END DO
      IF (MOD(p,2)==1) det=-det
!      det=det*(-1)**p
!      A=Orig
    END IF
  END FUNCTION d_det2
  SUBROUTINE r_write_mat(A,to)
    IMPLICIT NONE
    REAL(KIND=KIND(0.d0)) :: A(:,:)
    INTEGER :: i, j, N, to
    N=SIZE(A(1,:))
    WRITE(to,'(A)',advance='no') '#------------'
    DO i=2, N-1
      WRITE(to,'(A)',advance='no') '-------------'
    END DO
    WRITE(to,'(A)',advance='yes') '-------------'
    DO i=1, SIZE(A(:,1))
      DO j=1, N-1
    WRITE(to,'(ES14.5E3)',advance='no') A(i,j)
      END DO
      j=N
      WRITE(to,'(ES14.5E3)',advance='yes')  A(i,j)
    END DO
    WRITE(to,'(A)',advance='no') '#------------'
    DO i=2, N-1
      WRITE(to,'(A)',advance='no') '-------------'
    END DO
    WRITE(to,'(A)',advance='yes') '-------------'
  END SUBROUTINE r_write_mat
  SUBROUTINE c_write_mat(A,to)
    IMPLICIT NONE
    COMPLEX(KIND=KIND(0.d0)) :: A(:,:)
    INTEGER :: i, j, N, to
    N=SIZE(A(:,1))
    !  DO i=1, N-1
    !    WRITE(to,'(A)',advance='no') '-------------'
    !  END DO
    !  WRITE(to,'(A)',advance='yes') '-------------'
    DO i=1, SIZE(A(1,:))
      DO j=1, N-1
    WRITE(to,'(A,ES12.5,ES13.5,A)',advance='no') '(',A(i,j),') '
      END DO
      j=N
      WRITE(to,'(A,ES12.5,ES13.5,A)',advance='yes')  '(',A(i,j),') '
    END DO
    DO i=1, N-1
      WRITE(to,'(A)',advance='no') '-------------'
    END DO
    WRITE(to,'(A)',advance='yes') '-------------'
  END SUBROUTINE c_write_mat
  SUBROUTINE i_write_mat(A,to)
    IMPLICIT NONE
    INTEGER :: A(:,:)
    INTEGER :: i, j, N, to
    N=SIZE(A(:,1))
    DO i=1, N-1
      WRITE(to,'(A)',advance='no') '-------------'
    END DO
    WRITE(to,'(A)',advance='yes') '-------------'
    DO i=1, SIZE(A(1,:))
      DO j=1, N-1
    WRITE(to,'(i10)',advance='no') A(i,j)
      END DO
      j=N
      WRITE(to,'(i10)',advance='yes') A(i,j)
    END DO
    DO i=1, N-1
      WRITE(to,'(A)',advance='no') '-------------'
    END DO
    WRITE(to,'(A)',advance='yes') '-------------'
  END SUBROUTINE i_write_mat
  SUBROUTINE l_write_mat(A,to)
    IMPLICIT NONE
    LOGICAL :: A(:,:)
    INTEGER :: i, j, N, to
    N=SIZE(A(:,1))
    DO i=1, N-1
      WRITE(to,'(A)',advance='no') '-------------'
    END DO
    WRITE(to,'(A)',advance='yes') '-------------'
    DO i=1, SIZE(A(1,:))
      DO j=1, N-1
    WRITE(to,'(l6)',advance='no') A(i,j)
      END DO
      j=N
      WRITE(to,'(l6)',advance='yes') A(i,j)
    END DO
    DO i=1, N-1
      WRITE(to,'(A)',advance='no') '-------------'
    END DO
    WRITE(to,'(A)',advance='yes') '-------------'
  END SUBROUTINE l_write_mat
  SUBROUTINE r_write_mat1(A,to)
    IMPLICIT NONE
    REAL(KIND=KIND(0.d0)) :: A(:)
    INTEGER :: j, N, to
    N=SIZE(A(:))
    DO j=1, N-1
      WRITE(to,'(ES13.5)',advance='no') A(j)
    END DO
    j=N
    WRITE(to,'(ES13.5)',advance='yes')  A(j)
  END SUBROUTINE r_write_mat1

END MODULE matrixops
