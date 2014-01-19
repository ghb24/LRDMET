MODULE MinAlgos    
use const, only: dp
use Globals, only: nImp

IMPLICIT NONE
PRIVATE
PUBLIC :: minim, uobyqa
!minim: Simplex (no gradients required) (Bottom of file)
!uobyqa: Powell (no gradients required) (First routines)

CONTAINS

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% uobyqa.f %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
SUBROUTINE uobyqa(n, x, rhobeg, rhoend, iprint, maxfun, f, calfun, n_dum, g_dum, tMataxis_dum, tGrid, Freqs, Weights)

INTEGER, INTENT(IN)        :: n
REAL (dp), INTENT(IN OUT)  :: x(:)
REAL (dp), INTENT(IN)      :: rhobeg
REAL (dp), INTENT(IN)      :: rhoend
INTEGER, INTENT(IN)        :: iprint
INTEGER, INTENT(IN)        :: maxfun
REAL(dp), INTENT(OUT)      :: f
INTEGER, INTENT(IN)        :: n_dum
COMPLEX(dp), INTENT(IN)    :: g_dum(nImp,nImp,n_dum)
LOGICAL, INTENT(IN)        :: tMataxis_dum
LOGICAL, INTENT(IN)        :: tGrid
REAL(dp), INTENT(IN)       :: Freqs(n_dum)
REAL(dp), INTENT(IN)       :: Weights(n_dum)

INTERFACE
  SUBROUTINE calfun(x, f, n_dum, n, g_dum, tMataxis_dum, tGrid, Freqs, Weights, dJac)
    use const, only: dp
    use Globals, only: nImp
    IMPLICIT NONE
    INTEGER, INTENT(IN)    :: n, n_dum  !n_dum is number of frequency points
    REAL (dp), INTENT(IN)  :: x(n)
    REAL (dp), INTENT(OUT) :: f
    complex(dp), intent(in) :: g_dum(nImp,nImp,n_dum)
    logical, intent(in) :: tMataxis_dum
    logical, intent(in) :: tGrid
    real(dp), intent(in), optional :: Freqs(n_dum)
    real(dp), intent(in), optional :: Weights(n_dum)
    real(dp), intent(out), optional :: dJac(n)
  END SUBROUTINE calfun
END INTERFACE

!     This subroutine seeks the least value of a function of many variables,
!     by a trust region method that forms quadratic models by interpolation.
!     The algorithm is described in "UOBYQA: unconstrained optimization by
!     quadratic approximation" by M.J.D. Powell, Report DAMTP 2000/NA14,
!     University of Cambridge. The arguments of the subroutine are as follows.

!     N must be set to the number of variables and must be at least two.
!     Initial values of the variables must be set in X(1),X(2),...,X(N). They
!       will be changed to the values that give the least calculated F.
!     RHOBEG and RHOEND must be set to the initial and final values of a trust
!       region radius, so both must be positive with RHOEND<=RHOBEG. Typically
!       RHOBEG should be about one tenth of the greatest expected change to a
!       variable, and RHOEND should indicate the accuracy that is required in
!       the final values of the variables.
!     The value of IPRINT should be set to 0, 1, 2 or 3, which controls the
!       amount of printing.  Specifically, there is no output if IPRINT=0 and
!       there is output only at the return if IPRINT=1.  Otherwise, each new
!       value of RHO is printed, with the best vector of variables so far and
!       the corresponding value of the objective function. Further, each new
!       value of F with its variables are output if IPRINT=3.
!     MAXFUN must be set to an upper bound on the number of calls of CALFUN.
!     The array W will be used for working space. Its length must be at least
!       ( N**4 + 8*N**3 + 23*N**2 + 42*N + max [ 2*N**2 + 4, 18*N ] ) / 4.

!     SUBROUTINE CALFUN (N,X,F) must be provided by the user. It must set F to
!     the value of the objective function for the variables X(1),X(2),...,X(N).

!   Added, ghb24: calfun is the name of the external routine to calculate the function.

INTEGER  :: npt

!     Partition the working space array, so that different parts of it can be
!     treated separately by the subroutine that performs the main calculation.

npt = (n*n + 3*n + 2) / 2
CALL uobyqb(n, x, rhobeg, rhoend, iprint, maxfun, f, npt, calfun, n_dum, g_dum, tMataxis_dum, tGrid, Freqs, Weights)
RETURN
END SUBROUTINE uobyqa

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% uobyqb.f %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SUBROUTINE uobyqb(n, x, rhobeg, rhoend, iprint, maxfun, f, npt, calfun, n_dum, g_dum, tMataxis_dum, tGrid, Freqs, Weights)

INTEGER, INTENT(IN)        :: n
REAL (dp), INTENT(IN OUT)  :: x(:)
REAL (dp), INTENT(IN)      :: rhobeg
REAL (dp), INTENT(IN)      :: rhoend
INTEGER, INTENT(IN)        :: iprint
INTEGER, INTENT(IN)        :: maxfun
REAL(dp), INTENT(OUT)      :: f     !Final value of function
INTEGER, INTENT(IN)        :: npt
INTEGER, INTENT(IN)        :: n_dum
COMPLEX(dp), INTENT(IN)    :: g_dum(nImp,nImp,n_dum)
LOGICAL, INTENT(IN)        :: tMataxis_dum
LOGICAL, INTENT(IN)        :: tGrid
REAL(dp), INTENT(IN) :: Freqs(n_dum)
REAL(dp), INTENT(IN) :: Weights(n_dum)

INTERFACE
  SUBROUTINE calfun(x, f, n_dum, n, g_dum, tMataxis_dum, tGrid, Freqs, Weights, dJac)
    use const, only: dp
    use Globals, only: nImp
    IMPLICIT NONE
    INTEGER, INTENT(IN)    :: n, n_dum  !n_dum is number of frequency points
    REAL (dp), INTENT(IN)  :: x(n)
    REAL (dp), INTENT(OUT) :: f
    complex(dp), intent(in) :: g_dum(nImp,nImp,n_dum)
    logical, intent(in) :: tMataxis_dum
    logical, intent(in) :: tGrid
    real(dp), intent(in), optional :: Freqs(n_dum)
    real(dp), intent(in), optional :: Weights(n_dum)
    real(dp), intent(out), optional :: dJac(n)
  END SUBROUTINE calfun
END INTERFACE

! The following arrays were previously passed as arguments:

REAL (dp)  :: xbase(n), xopt(n), xnew(n), xpt(npt,n), pq(npt-1)
REAL (dp)  :: pl(npt,npt-1), h(n,n), g(n), d(n), vlag(npt), w(npt)

!     The arguments N, X, RHOBEG, RHOEND, IPRINT and MAXFUN are identical to
!       the corresponding arguments in SUBROUTINE UOBYQA.

!     NPT is set by UOBYQA to (N*N+3*N+2)/2 for the above dimension statement.
!     XBASE will contain a shift of origin that reduces the contributions from
!       rounding errors to values of the model and Lagrange functions.
!     XOPT will be set to the displacement from XBASE of the vector of
!       variables that provides the least calculated F so far.
!     XNEW will be set to the displacement from XBASE of the vector of
!       variables for the current calculation of F.
!     XPT will contain the interpolation point coordinates relative to XBASE.
!     PQ will contain the parameters of the quadratic model.
!     PL will contain the parameters of the Lagrange functions.
!     H will provide the second derivatives that TRSTEP and LAGMAX require.
!     G will provide the first derivatives that TRSTEP and LAGMAX require.
!     D is reserved for trial steps from XOPT, except that it will contain
!       diagonal second derivatives during the initialization procedure.
!     VLAG will contain the values of the Lagrange functions at a new point X.
!     The array W will be used for working space.

REAL (dp)  :: half = 0.5_dp, one = 1.0_dp, tol = 0.01_dp, two = 2.0_dp
REAL (dp)  :: zero = 0.0_dp
REAL (dp)  :: ddknew, delta, detrat, diff, distest, dnorm, errtol, estim
REAL (dp)  :: evalue, fbase, fopt, fsave, ratio, rho, rhosq, sixthm
REAL (dp)  :: sum, sumg, sumh, temp, tempa, tworsq, vmax, vquad, wmult
INTEGER    :: i, ih, ip, iq, iw, j, jswitch, k, knew, kopt, ksave, ktemp
INTEGER    :: nf, nftest, nnp, nptm

!     Set some constants.

nnp = n + n + 1
nptm = npt - 1
nftest = MAX(maxfun,1)

!     Initialization. NF is the number of function calculations so far.

rho = rhobeg
rhosq = rho * rho
nf = 0
DO  i = 1, n
  xbase(i) = x(i)
  xpt(1:npt,i) = zero
END DO
pl(1:npt,1:nptm) = zero

!     The branch to label 120 obtains a new value of the objective function
!     and then there is a branch back to label 50, because the new function
!     value is needed to form the initial quadratic model. The least function
!     value so far and its index are noted below.

50 x(1:n) = xbase(1:n) + xpt(nf+1,1:n)
GO TO 150

70 IF (nf == 1) THEN
  fopt = f
  kopt = nf
  fbase = f
  j = 0
  jswitch = -1
  ih = n
ELSE
  IF (f < fopt) THEN
    fopt = f
    kopt = nf
  END IF
END IF

!     Form the gradient and diagonal second derivatives of the initial
!     quadratic model and Lagrange functions.

IF (nf <= nnp) THEN
  jswitch = -jswitch
  IF (jswitch > 0) THEN
    IF (j >= 1) THEN
      ih = ih + j
      IF (w(j) < zero) THEN
        d(j) = (fsave+f-two*fbase) / rhosq
        pq(j) = (fsave-f) / (two*rho)
        pl(1,ih) = -two / rhosq
        pl(nf-1,j) = half / rho
        pl(nf-1,ih) = one / rhosq
      ELSE
        pq(j) = (4.0D0*fsave-3.0D0*fbase-f) / (two*rho)
        d(j) = (fbase+f-two*fsave) / rhosq
        pl(1,j) = -1.5D0 / rho
        pl(1,ih) = one / rhosq
        pl(nf-1,j) = two / rho
        pl(nf-1,ih) = -two / rhosq
      END IF
      pq(ih) = d(j)
      pl(nf,j) = -half / rho
      pl(nf,ih) = one / rhosq
    END IF
    
!     Pick the shift from XBASE to the next initial interpolation point
!     that provides diagonal second derivatives.
    
    IF (j < n) THEN
      j = j + 1
      xpt(nf+1,j) = rho
    END IF
  ELSE
    fsave = f
    IF (f < fbase) THEN
      w(j) = rho
      xpt(nf+1,j) = two * rho
    ELSE
      w(j) = -rho
      xpt(nf+1,j) = -rho
    END IF
  END IF
  IF (nf < nnp) GO TO 50
  
!     Form the off-diagonal second derivatives of the initial quadratic model.
  
  ih = n
  ip = 1
  iq = 2
END IF
ih = ih + 1
IF (nf > nnp) THEN
  temp = one / (w(ip)*w(iq))
  tempa = f - fbase - w(ip) * pq(ip) - w(iq) * pq(iq)
  pq(ih) = (tempa - half*rhosq*(d(ip)+d(iq))) * temp
  pl(1,ih) = temp
  iw = ip + ip
  IF (w(ip) < zero) iw = iw + 1
  pl(iw,ih) = -temp
  iw = iq + iq
  IF (w(iq) < zero) iw = iw + 1
  pl(iw,ih) = -temp
  pl(nf,ih) = temp
  
!     Pick the shift from XBASE to the next initial interpolation point
!     that provides off-diagonal second derivatives.
  
  ip = ip + 1
END IF
IF (ip == iq) THEN
  ih = ih + 1
  ip = 1
  iq = iq + 1
END IF
IF (nf < npt) THEN
  xpt(nf+1,ip) = w(ip)
  xpt(nf+1,iq) = w(iq)
  GO TO 50
END IF

!     Set parameters to begin the iterations for the current RHO.

sixthm = zero
delta = rho
80 tworsq = (two*rho) ** 2
rhosq = rho * rho

!     Form the gradient of the quadratic model at the trust region centre.

90 knew = 0
ih = n
DO  j = 1, n
  xopt(j) = xpt(kopt,j)
  g(j) = pq(j)
  DO  i = 1, j
    ih = ih + 1
    g(i) = g(i) + pq(ih) * xopt(j)
    IF (i < j) g(j) = g(j) + pq(ih) * xopt(i)
    h(i,j) = pq(ih)
  END DO
END DO

!     Generate the next trust region step and test its length.  Set KNEW
!     to -1 if the purpose of the next F will be to improve conditioning,
!     and also calculate a lower bound on the Hessian term of the model Q.

CALL trstep(n, g, h, delta, tol, d, evalue)
temp = zero
DO i = 1, n
  temp = temp + d(i)**2
END DO
dnorm = MIN(delta,SQRT(temp))
errtol = -one
IF (dnorm < half*rho) THEN
  knew = -1
  errtol = half * evalue * rho * rho
  IF (nf <= npt+9) errtol = zero
  GO TO 290
END IF

!     Calculate the next value of the objective function.

130 DO  i = 1, n
  xnew(i) = xopt(i) + d(i)
  x(i) = xbase(i) + xnew(i)
END DO
150 IF (nf >= nftest) THEN
  IF (iprint > 0) WRITE(*, 5000)
  GO TO 420
END IF
nf = nf + 1
CALL calfun(x, f, n_dum, n, g_dum, tMataxis_dum, tGrid, Freqs, Weights)
IF (iprint == 3) THEN
  WRITE(*, 5100) nf, f, x(1:n)
END IF
IF (nf <= npt) GO TO 70
IF (knew == -1) GO TO 420

!     Use the quadratic model to predict the change in F due to the step D,
!     and find the values of the Lagrange functions at the new point.

vquad = zero
ih = n
DO  j = 1, n
  w(j) = d(j)
  vquad = vquad + w(j) * pq(j)
  DO  i = 1, j
    ih = ih + 1
    w(ih) = d(i) * xnew(j) + d(j) * xopt(i)
    IF (i == j) w(ih) = half * w(ih)
    vquad = vquad + w(ih) * pq(ih)
  END DO
END DO
DO  k = 1, npt
  temp = zero
  DO  j = 1, nptm
    temp = temp + w(j) * pl(k,j)
  END DO
  vlag(k) = temp
END DO
vlag(kopt) = vlag(kopt) + one

!     Update SIXTHM, which is a lower bound on one sixth of the greatest
!     third derivative of F.

diff = f - fopt - vquad
sum = zero
DO  k = 1, npt
  temp = zero
  DO  i = 1, n
    temp = temp + (xpt(k,i)-xnew(i)) ** 2
  END DO
  temp = SQRT(temp)
  sum = sum + ABS(temp*temp*temp*vlag(k))
END DO
sixthm = MAX(sixthm, ABS(diff)/sum)

!     Update FOPT and XOPT if the new F is the least value of the objective
!     function so far.  Then branch if D is not a trust region step.

fsave = fopt
IF (f < fopt) THEN
  fopt = f
  xopt(1:n) = xnew(1:n)
END IF
ksave = knew
IF (knew <= 0) THEN
  
!     Pick the next value of DELTA after a trust region step.
  
  IF (vquad >= zero) THEN
    IF (iprint > 0) WRITE(*, 5200)
    GO TO 420
  END IF
  ratio = (f-fsave) / vquad
  IF (ratio <= 0.1D0) THEN
    delta = half * dnorm
  ELSE IF (ratio <= 0.7D0) THEN
    delta = MAX(half*delta,dnorm)
  ELSE
    delta = MAX(delta, 1.25D0*dnorm, dnorm+rho)
  END IF
  IF (delta <= 1.5D0*rho) delta = rho
  
!     Set KNEW to the index of the next interpolation point to be deleted.
  
  ktemp = 0
  detrat = zero
  IF (f >= fsave) THEN
    ktemp = kopt
    detrat = one
  END IF
  DO  k = 1, npt
    sum = zero
    DO  i = 1, n
      sum = sum + (xpt(k,i)-xopt(i)) ** 2
    END DO
    temp = ABS(vlag(k))
    IF (sum > rhosq) temp = temp * (sum/rhosq) ** 1.5D0
    IF (temp > detrat .AND. k /= ktemp) THEN
      detrat = temp
      ddknew = sum
      knew = k
    END IF
  END DO
  IF (knew == 0) GO TO 290
END IF

!     Replace the interpolation point that has index KNEW by the point XNEW,
!     and also update the Lagrange functions and the quadratic model.

DO  i = 1, n
  xpt(knew,i) = xnew(i)
END DO
temp = one / vlag(knew)
DO  j = 1, nptm
  pl(knew,j) = temp * pl(knew,j)
  pq(j) = pq(j) + diff * pl(knew,j)
END DO
DO  k = 1, npt
  IF (k /= knew) THEN
    temp = vlag(k)
    DO  j = 1, nptm
      pl(k,j) = pl(k,j) - temp * pl(knew,j)
    END DO
  END IF
END DO

!     Update KOPT if F is the least calculated value of the objective function.
!     Then branch for another trust region calculation.  The case KSAVE > 0
!     indicates that a model step has just been taken.

IF (f < fsave) THEN
  kopt = knew
  GO TO 90
END IF
IF (ksave > 0) GO TO 90
IF (dnorm > two*rho) GO TO 90
IF (ddknew > tworsq) GO TO 90

!     Alternatively, find out if the interpolation points are close
!     enough to the best point so far.

290 DO  k = 1, npt
  w(k) = zero
  DO  i = 1, n
    w(k) = w(k) + (xpt(k,i)-xopt(i)) ** 2
  END DO
END DO
320 knew = -1
distest = tworsq
DO  k = 1, npt
  IF (w(k) > distest) THEN
    knew = k
    distest = w(k)
  END IF
END DO

!     If a point is sufficiently far away, then set the gradient and Hessian
!     of its Lagrange function at the centre of the trust region, and find
!     half the sum of squares of components of the Hessian.

IF (knew > 0) THEN
  ih = n
  sumh = zero
  DO  j = 1, n
    g(j) = pl(knew,j)
    DO  i = 1, j
      ih = ih + 1
      temp = pl(knew,ih)
      g(j) = g(j) + temp * xopt(i)
      IF (i < j) THEN
        g(i) = g(i) + temp * xopt(j)
        sumh = sumh + temp * temp
      END IF
      h(i,j) = temp
    END DO
    sumh = sumh + half * temp * temp
  END DO
  
!     If ERRTOL is positive, test whether to replace the interpolation point
!     with index KNEW, using a bound on the maximum modulus of its Lagrange
!     function in the trust region.
  
  IF (errtol > zero) THEN
    w(knew) = zero
    sumg = zero
    DO  i = 1, n
      sumg = sumg + g(i) ** 2
    END DO
    estim = rho * (SQRT(sumg)+rho*SQRT(half*sumh))
    wmult = sixthm * distest ** 1.5D0
    IF (wmult*estim <= errtol) GO TO 320
  END IF
  
!     If the KNEW-th point may be replaced, then pick a D that gives a large
!     value of the modulus of its Lagrange function within the trust region.
!     Here the vector XNEW is used as temporary working space.
  
  CALL lagmax(n, g, h, rho, d, xnew, vmax)
  IF (errtol > zero) THEN
    IF (wmult*vmax <= errtol) GO TO 320
  END IF
  GO TO 130
END IF
IF (dnorm > rho) GO TO 90

!     Prepare to reduce RHO by shifting XBASE to the best point so far,
!     and make the corresponding changes to the gradients of the Lagrange
!     functions and the quadratic model.

IF (rho > rhoend) THEN
  ih = n
  DO  j = 1, n
    xbase(j) = xbase(j) + xopt(j)
    DO  k = 1, npt
      xpt(k,j) = xpt(k,j) - xopt(j)
    END DO
    DO  i = 1, j
      ih = ih + 1
      pq(i) = pq(i) + pq(ih) * xopt(j)
      IF (i < j) THEN
        pq(j) = pq(j) + pq(ih) * xopt(i)
        DO  k = 1, npt
          pl(k,j) = pl(k,j) + pl(k,ih) * xopt(i)
        END DO
      END IF
      DO  k = 1, npt
        pl(k,i) = pl(k,i) + pl(k,ih) * xopt(j)
      END DO
    END DO
  END DO
  
!     Pick the next values of RHO and DELTA.
  
  delta = half * rho
  ratio = rho / rhoend
  IF (ratio <= 16.0D0) THEN
    rho = rhoend
  ELSE IF (ratio <= 250.0D0) THEN
    rho = SQRT(ratio) * rhoend
  ELSE
    rho = 0.1D0 * rho
  END IF
  delta = MAX(delta,rho)
  IF (iprint >= 2) THEN
    IF (iprint >= 3) WRITE(*, 5300)
    WRITE(*, 5400) rho, nf
    WRITE(*, 5500) fopt, xbase(1:n)
  END IF
  GO TO 80
END IF

!     Return from the calculation, after another Newton-Raphson step, if
!     it is too short to have been tried before.

IF (errtol >= zero) GO TO 130
420 IF (fopt <= f) THEN
  DO  i = 1, n
    x(i) = xbase(i) + xopt(i)
  END DO
  f = fopt
END IF
IF (iprint >= 1) THEN
  WRITE(*, 5600) nf
  WRITE(*, 5500) f, x(1:n)
END IF
RETURN

5000 FORMAT (/T5, 'Return from UOBYQA because CALFUN has been',  &
    ' called MAXFUN times')
5100 FORMAT (/T5, 'Function number',i6,'    F =', g18.10,  &
    '    The corresponding X is:'/ (t3, 5g15.6))
5200 FORMAT (/T5, 'Return from UOBYQA because a trust',  &
    ' region step has failed to reduce Q')
5300 FORMAT (' ')
5400 FORMAT (/T5, 'New RHO =', g11.4, '     Number of function values =',i6)
5500 FORMAT (T5, 'Least value of F =', g23.15,  &
    '         The corresponding X is:'/ (t3, 5g15.6))
5600 FORMAT (/T5, 'At the return from UOBYQA',  &
    '     Number of function values =', i6)
END SUBROUTINE uobyqb

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% trstep.f %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SUBROUTINE trstep(n, g, h, delta, tol, d, evalue)

INTEGER, INTENT(IN)        :: n
REAL (dp), INTENT(IN)      :: g(:)
REAL (dp), INTENT(IN OUT)  :: h(:,:)
REAL (dp), INTENT(IN)      :: delta
REAL (dp), INTENT(IN)      :: tol
REAL (dp), INTENT(OUT)     :: d(:)
REAL (dp), INTENT(OUT)     :: evalue

!     N is the number of variables of a quadratic objective function, Q say.
!     G is the gradient of Q at the origin.
!     H is the Hessian matrix of Q.  Only the upper triangular and diagonal
!       parts need be set.  The lower triangular part is used to store the
!       elements of a Householder similarity transformation.
!     DELTA is the trust region radius, and has to be positive.
!     TOL is the value of a tolerance from the open interval (0,1).
!     D will be set to the calculated vector of variables.

!     EVALUE will be set to the least eigenvalue of H if and only if D is a
!     Newton-Raphson step.  Then EVALUE will be positive, but otherwise it
!     will be set to zero.

!     Let MAXRED be the maximum of Q(0)-Q(D) subject to ||D|| <= DELTA,
!     and let ACTRED be the value of Q(0)-Q(D) that is actually calculated.
!     We take the view that any D is acceptable if it has the properties

!             ||D|| <= DELTA  and  ACTRED <= (1-TOL)*MAXRED.

!     The calculation of D is done by the method of Section 2 of the paper
!     by MJDP in the 1997 Dundee Numerical Analysis Conference Proceedings,
!     after transforming H to tridiagonal form.

!     The arrays GG, TD, TN, W, PIV and Z will be used for working space.
REAL (dp)  :: gg(n), td(n), tn(n), w(n), piv(n), z(n)

REAL (dp)  :: delsq, dhd, dnorm, dsq, dtg, dtz, gam, gnorm, gsq, hnorm
REAL (dp)  :: par, parl, parlest, paru, paruest, phi, phil, phiu, pivksv
REAL (dp)  :: pivot, posdef, scale, shfmax, shfmin, shift, slope, sum
REAL (dp)  :: tdmin, temp, tempa, tempb, wsq, wwsq, wz, zsq
INTEGER    :: i, iterc, j, jp, k, kp, kpp, ksav, ksave, nm
REAL (dp)  :: one = 1.0_dp, two = 2.0_dp, zero = 0.0_dp

!     Initialization.

delsq = delta * delta
evalue = zero
nm = n - 1
DO  i = 1, n
  d(i) = zero
  td(i) = h(i,i)
  DO  j = 1, i
    h(i,j) = h(j,i)
  END DO
END DO

!     Apply Householder transformations to obtain a tridiagonal matrix that
!     is similar to H, and put the elements of the Householder vectors in
!     the lower triangular part of H.  Further, TD and TN will contain the
!     diagonal and other nonzero elements of the tridiagonal matrix.

DO  k = 1, nm
  kp = k + 1
  sum = zero
  IF (kp < n) THEN
    kpp = kp + 1
    DO  i = kpp, n
      sum = sum + h(i,k) ** 2
    END DO
  END IF
  IF (sum == zero) THEN
    tn(k) = h(kp,k)
    h(kp,k) = zero
  ELSE
    temp = h(kp,k)
    tn(k) = SIGN(SQRT(sum+temp*temp),temp)
    h(kp,k) = -sum / (temp+tn(k))
    temp = SQRT(two/(sum+h(kp,k)**2))
    DO  i = kp, n
      w(i) = temp * h(i,k)
      h(i,k) = w(i)
      z(i) = td(i) * w(i)
    END DO
    wz = zero
    DO  j = kp, nm
      jp = j + 1
      DO  i = jp, n
        z(i) = z(i) + h(i,j) * w(j)
        z(j) = z(j) + h(i,j) * w(i)
      END DO
      wz = wz + w(j) * z(j)
    END DO
    wz = wz + w(n) * z(n)
    DO  j = kp, n
      td(j) = td(j) + w(j) * (wz*w(j)-two*z(j))
      IF (j < n) THEN
        jp = j + 1
        DO  i = jp, n
          h(i,j) = h(i,j) - w(i) * z(j) - w(j) * (z(i)-wz*w(i))
        END DO
      END IF
    END DO
  END IF
END DO

!     Form GG by applying the similarity transformation to G.

gsq = zero
DO  i = 1, n
  gg(i) = g(i)
  gsq = gsq + g(i) ** 2
END DO
gnorm = SQRT(gsq)
DO  k = 1, nm
  kp = k + 1
  sum = zero
  DO  i = kp, n
    sum = sum + gg(i) * h(i,k)
  END DO
  DO  i = kp, n
    gg(i) = gg(i) - sum * h(i,k)
  END DO
END DO

!     Begin the trust region calculation with a tridiagonal matrix by
!     calculating the norm of H. Then treat the case when H is zero.

hnorm = ABS(td(1)) + ABS(tn(1))
tdmin = td(1)
tn(n) = zero
DO  i = 2, n
  temp = ABS(tn(i-1)) + ABS(td(i)) + ABS(tn(i))
  hnorm = MAX(hnorm,temp)
  tdmin = MIN(tdmin,td(i))
END DO
IF (hnorm == zero) THEN
  IF (gnorm == zero) GO TO 420
  scale = delta / gnorm
  DO  i = 1, n
    d(i) = -scale * gg(i)
  END DO
  GO TO 380
END IF

!     Set the initial values of PAR and its bounds.

parl = MAX(zero, -tdmin, gnorm/delta-hnorm)
parlest = parl
par = parl
paru = zero
paruest = zero
posdef = zero
iterc = 0

!     Calculate the pivots of the Cholesky factorization of (H+PAR*I).

160 iterc = iterc + 1
ksav = 0
piv(1) = td(1) + par
k = 1
170 IF (piv(k) > zero) THEN
  piv(k+1) = td(k+1) + par - tn(k) ** 2 / piv(k)
ELSE
  IF (piv(k) < zero .OR. tn(k) /= zero) GO TO 180
  ksav = k
  piv(k+1) = td(k+1) + par
END IF
k = k + 1
IF (k < n) GO TO 170
IF (piv(k) >= zero) THEN
  IF (piv(k) == zero) ksav = k
  
!     Branch if all the pivots are positive, allowing for the case when
!     G is zero.
  
  IF (ksav == 0 .AND. gsq > zero) GO TO 250
  IF (gsq == zero) THEN
    IF (par == zero) GO TO 380
    paru = par
    paruest = par
    IF (ksav == 0) GO TO 210
  END IF
  k = ksav
END IF

!     Set D to a direction of nonpositive curvature of the given tridiagonal
!     matrix, and thus revise PARLEST.

180 d(k) = one
IF (ABS(tn(k)) <= ABS(piv(k))) THEN
  dsq = one
  dhd = piv(k)
ELSE
  temp = td(k+1) + par
  IF (temp <= ABS(piv(k))) THEN
    d(k+1) = SIGN(one,-tn(k))
    dhd = piv(k) + temp - two * ABS(tn(k))
  ELSE
    d(k+1) = -tn(k) / temp
    dhd = piv(k) + tn(k) * d(k+1)
  END IF
  dsq = one + d(k+1) ** 2
END IF
190 IF (k > 1) THEN
  k = k - 1
  IF (tn(k) /= zero) THEN
    d(k) = -tn(k) * d(k+1) / piv(k)
    dsq = dsq + d(k) ** 2
    GO TO 190
  END IF
  d(1:k) = zero
END IF
parl = par
parlest = par - dhd / dsq

!     Terminate with D set to a multiple of the current D if the following
!     test suggests that it suitable to do so.

210 temp = paruest
IF (gsq == zero) temp = temp * (one-tol)
IF (paruest > zero .AND. parlest >= temp) THEN
  dtg = DOT_PRODUCT( d(1:n), gg(1:n) )
  scale = -SIGN(delta/SQRT(dsq),dtg)
  d(1:n) = scale * d(1:n)
  GO TO 380
END IF

!     Pick the value of PAR for the next iteration.

240 IF (paru == zero) THEN
  par = two * parlest + gnorm / delta
ELSE
  par = 0.5D0 * (parl+paru)
  par = MAX(par,parlest)
END IF
IF (paruest > zero) par = MIN(par,paruest)
GO TO 160

!     Calculate D for the current PAR in the positive definite case.

250 w(1) = -gg(1) / piv(1)
DO  i = 2, n
  w(i) = (-gg(i)-tn(i-1)*w(i-1)) / piv(i)
END DO
d(n) = w(n)
DO  i = nm, 1, -1
  d(i) = w(i) - tn(i) * d(i+1) / piv(i)
END DO

!     Branch if a Newton-Raphson step is acceptable.

dsq = zero
wsq = zero
DO  i = 1, n
  dsq = dsq + d(i) ** 2
  wsq = wsq + piv(i) * w(i) ** 2
END DO
IF (par /= zero .OR. dsq > delsq) THEN
  
!     Make the usual test for acceptability of a full trust region step.
  
  dnorm = SQRT(dsq)
  phi = one / dnorm - one / delta
  temp = tol * (one+par*dsq/wsq) - dsq * phi * phi
  IF (temp >= zero) THEN
    scale = delta / dnorm
    DO  i = 1, n
      d(i) = scale * d(i)
    END DO
    GO TO 380
  END IF
  IF (iterc >= 2 .AND. par <= parl) GO TO 380
  IF (paru > zero .AND. par >= paru) GO TO 380
  
!     Complete the iteration when PHI is negative.
  
  IF (phi < zero) THEN
    parlest = par
    IF (posdef == one) THEN
      IF (phi <= phil) GO TO 380
      slope = (phi-phil) / (par-parl)
      parlest = par - phi / slope
    END IF
    slope = one / gnorm
    IF (paru > zero) slope = (phiu-phi) / (paru-par)
    temp = par - phi / slope
    IF (paruest > zero) temp = MIN(temp,paruest)
    paruest = temp
    posdef = one
    parl = par
    phil = phi
    GO TO 240
  END IF
  
!     If required, calculate Z for the alternative test for convergence.
  
  IF (posdef == zero) THEN
    w(1) = one / piv(1)
    DO  i = 2, n
      temp = -tn(i-1) * w(i-1)
      w(i) = (SIGN(one,temp)+temp) / piv(i)
    END DO
    z(n) = w(n)
    DO  i = nm, 1, -1
      z(i) = w(i) - tn(i) * z(i+1) / piv(i)
    END DO
    wwsq = zero
    zsq = zero
    dtz = zero
    DO  i = 1, n
      wwsq = wwsq + piv(i) * w(i) ** 2
      zsq = zsq + z(i) ** 2
      dtz = dtz + d(i) * z(i)
    END DO
    
!     Apply the alternative test for convergence.
    
    tempa = ABS(delsq-dsq)
    tempb = SQRT(dtz*dtz+tempa*zsq)
    gam = tempa / (SIGN(tempb,dtz)+dtz)
    temp = tol * (wsq+par*delsq) - gam * gam * wwsq
    IF (temp >= zero) THEN
      DO  i = 1, n
        d(i) = d(i) + gam * z(i)
      END DO
      GO TO 380
    END IF
    parlest = MAX(parlest,par-wwsq/zsq)
  END IF
  
!     Complete the iteration when PHI is positive.
  
  slope = one / gnorm
  IF (paru > zero) THEN
    IF (phi >= phiu) GO TO 380
    slope = (phiu-phi) / (paru-par)
  END IF
  parlest = MAX(parlest,par-phi/slope)
  paruest = par
  IF (posdef == one) THEN
    slope = (phi-phil) / (par-parl)
    paruest = par - phi / slope
  END IF
  paru = par
  phiu = phi
  GO TO 240
END IF

!     Set EVALUE to the least eigenvalue of the second derivative matrix if
!     D is a Newton-Raphson step. SHFMAX will be an upper bound on EVALUE.

shfmin = zero
pivot = td(1)
shfmax = pivot
DO  k = 2, n
  pivot = td(k) - tn(k-1) ** 2 / pivot
  shfmax = MIN(shfmax,pivot)
END DO

!     Find EVALUE by a bisection method, but occasionally SHFMAX may be
!     adjusted by the rule of false position.

ksave = 0
350 shift = 0.5D0 * (shfmin+shfmax)
k = 1
temp = td(1) - shift

360 IF (temp > zero) THEN
  piv(k) = temp
  IF (k < n) THEN
    temp = td(k+1) - shift - tn(k) ** 2 / temp
    k = k + 1
    GO TO 360
  END IF
  shfmin = shift
ELSE
  IF (k < ksave) GO TO 370
  IF (k == ksave) THEN
    IF (pivksv == zero) GO TO 370
    IF (piv(k)-temp < temp-pivksv) THEN
      pivksv = temp
      shfmax = shift
    ELSE
      pivksv = zero
      shfmax = (shift*piv(k) - shfmin*temp) / (piv(k)-temp)
    END IF
  ELSE
    ksave = k
    pivksv = temp
    shfmax = shift
  END IF
END IF
IF (shfmin <= 0.99D0*shfmax) GO TO 350
370 evalue = shfmin

!     Apply the inverse Householder transformations to D.

380 nm = n - 1
DO  k = nm, 1, -1
  kp = k + 1
  sum = zero
  DO  i = kp, n
    sum = sum + d(i) * h(i,k)
  END DO
  DO  i = kp, n
    d(i) = d(i) - sum * h(i,k)
  END DO
END DO

!     Return from the subroutine.

420 RETURN
END SUBROUTINE trstep

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% lagmax.f %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SUBROUTINE lagmax(n, g, h, rho, d, v, vmax)

INTEGER, INTENT(IN)     :: n
REAL (dp), INTENT(IN)   :: g(:)
REAL (dp), INTENT(OUT)  :: h(:,:)
REAL (dp), INTENT(IN)   :: rho
REAL (dp), INTENT(OUT)  :: d(:)
REAL (dp), INTENT(OUT)  :: v(:)
REAL (dp), INTENT(OUT)  :: vmax

!     N is the number of variables of a quadratic objective function, Q say.
!     G is the gradient of Q at the origin.
!     H is the symmetric Hessian matrix of Q. Only the upper triangular and
!       diagonal parts need be set.
!     RHO is the trust region radius, and has to be positive.
!     D will be set to the calculated vector of variables.
!     The array V will be used for working space.
!     VMAX will be set to |Q(0)-Q(D)|.

!     Calculating the D that maximizes |Q(0)-Q(D)| subject to ||D|| <= RHO
!     requires of order N**3 operations, but sometimes it is adequate if
!     |Q(0)-Q(D)| is within about 0.9 of its greatest possible value. This
!     subroutine provides such a solution in only of order N**2 operations,
!     where the claim of accuracy has been tested by numerical experiments.

REAL (dp)  :: half = 0.5_dp, one = 1.0_dp, zero = 0.0_dp
REAL (dp)  :: dd, dhd, dlin, dsq, gd, gg, ghg, gnorm, halfrt, hmax, ratio
REAL (dp)  :: scale, sum, sumv, temp, tempa, tempb, tempc, tempd, tempv
REAL (dp)  :: vhg, vhv, vhw, vlin, vmu, vnorm, vsq, vv, wcos, whw, wsin, wsq
INTEGER    :: i, j, k

!     Preliminary calculations.

halfrt = SQRT(half)

!     Pick V such that ||HV|| / ||V|| is large.

hmax = zero
DO  i = 1, n
  sum = zero
  DO  j = 1, n
    h(j,i) = h(i,j)
    sum = sum + h(i,j) ** 2
  END DO
  IF (sum > hmax) THEN
    hmax = sum
    k = i
  END IF
END DO
DO  j = 1, n
  v(j) = h(k,j)
END DO

!     Set D to a vector in the subspace spanned by V and HV that maximizes
!     |(D,HD)|/(D,D), except that we set D=HV if V and HV are nearly parallel.
!     The vector that has the name D at label 60 used to be the vector W.

vsq = zero
vhv = zero
dsq = zero
DO  i = 1, n
  vsq = vsq + v(i) ** 2
  d(i) = DOT_PRODUCT( h(i,1:n), v(1:n) )
  vhv = vhv + v(i) * d(i)
  dsq = dsq + d(i) ** 2
END DO
IF (vhv*vhv <= 0.9999D0*dsq*vsq) THEN
  temp = vhv / vsq
  wsq = zero
  DO  i = 1, n
    d(i) = d(i) - temp * v(i)
    wsq = wsq + d(i) ** 2
  END DO
  whw = zero
  ratio = SQRT(wsq/vsq)
  DO  i = 1, n
    temp = DOT_PRODUCT( h(i,1:n), d(1:n) )
    whw = whw + temp * d(i)
    v(i) = ratio * v(i)
  END DO
  vhv = ratio * ratio * vhv
  vhw = ratio * wsq
  temp = half * (whw-vhv)
  temp = temp + SIGN(SQRT(temp**2+vhw**2),whw+vhv)
  DO  i = 1, n
    d(i) = vhw * v(i) + temp * d(i)
  END DO
END IF

!     We now turn our attention to the subspace spanned by G and D. A multiple
!     of the current D is returned if that choice seems to be adequate.

gg = zero
gd = zero
dd = zero
dhd = zero
DO  i = 1, n
  gg = gg + g(i) ** 2
  gd = gd + g(i) * d(i)
  dd = dd + d(i) ** 2
  sum = DOT_PRODUCT( h(i,1:n), d(1:n) )
  dhd = dhd + sum * d(i)
END DO
temp = gd / gg
vv = zero
scale = SIGN(rho/SQRT(dd),gd*dhd)
DO  i = 1, n
  v(i) = d(i) - temp * g(i)
  vv = vv + v(i) ** 2
  d(i) = scale * d(i)
END DO
gnorm = SQRT(gg)
IF (gnorm*dd <= 0.5D-2*rho*ABS(dhd) .OR. vv/dd <= 1.0D-4) THEN
  vmax = ABS(scale*(gd + half*scale*dhd))
  GO TO 170
END IF

!     G and V are now orthogonal in the subspace spanned by G and D.  Hence
!     we generate an orthonormal basis of this subspace such that (D,HV) is
!     negligible or zero, where D and V will be the basis vectors.

ghg = zero
vhg = zero
vhv = zero
DO  i = 1, n
  sum = DOT_PRODUCT( h(i,1:n), g(1:n) )
  sumv = DOT_PRODUCT( h(i,1:n), v(1:n) )
  ghg = ghg + sum * g(i)
  vhg = vhg + sumv * g(i)
  vhv = vhv + sumv * v(i)
END DO
vnorm = SQRT(vv)
ghg = ghg / gg
vhg = vhg / (vnorm*gnorm)
vhv = vhv / vv
IF (ABS(vhg) <= 0.01D0*MAX(ABS(ghg),ABS(vhv))) THEN
  vmu = ghg - vhv
  wcos = one
  wsin = zero
ELSE
  temp = half * (ghg-vhv)
  vmu = temp + SIGN(SQRT(temp**2+vhg**2),temp)
  temp = SQRT(vmu**2+vhg**2)
  wcos = vmu / temp
  wsin = vhg / temp
END IF
tempa = wcos / gnorm
tempb = wsin / vnorm
tempc = wcos / vnorm
tempd = wsin / gnorm
DO  i = 1, n
  d(i) = tempa * g(i) + tempb * v(i)
  v(i) = tempc * v(i) - tempd * g(i)
END DO

!     The final D is a multiple of the current D, V, D+V or D-V. We make the
!     choice from these possibilities that is optimal.

dlin = wcos * gnorm / rho
vlin = -wsin * gnorm / rho
tempa = ABS(dlin) + half * ABS(vmu+vhv)
tempb = ABS(vlin) + half * ABS(ghg-vmu)
tempc = halfrt * (ABS(dlin)+ABS(vlin)) + 0.25D0 * ABS(ghg+vhv)
IF (tempa >= tempb .AND. tempa >= tempc) THEN
  tempd = SIGN(rho,dlin*(vmu+vhv))
  tempv = zero
ELSE IF (tempb >= tempc) THEN
  tempd = zero
  tempv = SIGN(rho,vlin*(ghg-vmu))
ELSE
  tempd = SIGN(halfrt*rho,dlin*(ghg+vhv))
  tempv = SIGN(halfrt*rho,vlin*(ghg+vhv))
END IF
DO  i = 1, n
  d(i) = tempd * d(i) + tempv * v(i)
END DO
vmax = rho * rho * MAX(tempa,tempb,tempc)
170 RETURN
END SUBROUTINE lagmax


!   BELOW ARE THE ROUTINES FOR SIMPLEX


SUBROUTINE minim(p, step, nop, func, maxfn, iprint, stopcr, nloop, iquad,  &
                 simp, var, functn, ifault, n_dum, g_dum, tMataxis_dum, tGrid, Freqs, Weights)

!     A PROGRAM FOR FUNCTION MINIMIZATION USING THE SIMPLEX METHOD.

!     FOR DETAILS, SEE NELDER & MEAD, THE COMPUTER JOURNAL, JANUARY 1965

!     PROGRAMMED BY D.E.SHAW,
!     CSIRO, DIVISION OF MATHEMATICS & STATISTICS
!     P.O. BOX 218, LINDFIELD, N.S.W. 2070

!     WITH AMENDMENTS BY R.W.M.WEDDERBURN
!     ROTHAMSTED EXPERIMENTAL STATION
!     HARPENDEN, HERTFORDSHIRE, ENGLAND

!     Further amended by Alan Miller
!     CSIRO Division of Mathematical & Information Sciences
!     Private Bag 10, CLAYTON, VIC. 3169

!     Fortran 90 conversion by Alan Miller, June 1995
!     Alan.Miller @ vic.cmis.csiro.au
!     Latest revision - 5 December 1999

!     ARGUMENTS:-
!     P()     = INPUT, STARTING VALUES OF PARAMETERS
!               OUTPUT, FINAL VALUES OF PARAMETERS
!     STEP()  = INPUT, INITIAL STEP SIZES
!     NOP     = INPUT, NO. OF PARAMETERS, INCL. ANY TO BE HELD FIXED
!     FUNC    = OUTPUT, THE FUNCTION VALUE CORRESPONDING TO THE FINAL
!                 PARAMETER VALUES.
!     maxfn     = INPUT, THE MAXIMUM NO. OF FUNCTION EVALUATIONS ALLOWED.
!               Say, 20 times the number of parameters, NOP.
!     IPRINT  = INPUT, PRINT CONTROL PARAMETER
!                 < 0 NO PRINTING
!                 = 0 PRINTING OF PARAMETER VALUES AND THE FUNCTION
!                     VALUE AFTER INITIAL EVIDENCE OF CONVERGENCE.
!                 > 0 AS FOR IPRINT = 0 PLUS PROGRESS REPORTS AFTER EVERY
!                     IPRINT EVALUATIONS, PLUS PRINTING FOR THE INITIAL SIMPLEX.
!     STOPCR  = INPUT, STOPPING CRITERION.
!               The criterion is applied to the standard deviation of
!               the values of FUNC at the points of the simplex.
!     NLOOP   = INPUT, THE STOPPING RULE IS APPLIED AFTER EVERY NLOOP
!               FUNCTION EVALUATIONS.   Normally NLOOP should be slightly
!               greater than NOP, say NLOOP = 2*NOP.
!     IQUAD   = INPUT, = 1 IF FITTING OF A QUADRATIC SURFACE IS REQUIRED
!                      = 0 IF NOT
!               N.B. The fitting of a quadratic surface is strongly
!               recommended, provided that the fitted function is
!               continuous in the vicinity of the minimum.   It is often
!               a good indicator of whether a premature termination of
!               the search has occurred.
!     SIMP    = INPUT, CRITERION FOR EXPANDING THE SIMPLEX TO OVERCOME
!               ROUNDING ERRORS BEFORE FITTING THE QUADRATIC SURFACE.
!               The simplex is expanded so that the function values at
!               the points of the simplex exceed those at the supposed
!               minimum by at least an amount SIMP.
!     VAR()   = OUTPUT, CONTAINS THE DIAGONAL ELEMENTS OF THE INVERSE OF
!               THE INFORMATION MATRIX.
!     FUNCTN  = INPUT, NAME OF THE USER'S SUBROUTINE - ARGUMENTS (P,FUNC)
!               WHICH RETURNS THE FUNCTION VALUE FOR A GIVEN SET OF
!               PARAMETER VALUES IN ARRAY P.
!****     FUNCTN MUST BE DECLARED EXTERNAL IN THE CALLING PROGRAM.
!     IFAULT  = OUTPUT, = 0 FOR SUCCESSFUL TERMINATION
!                 = 1 IF MAXIMUM NO. OF FUNCTION EVALUATIONS EXCEEDED
!                 = 2 IF INFORMATION MATRIX IS NOT +VE SEMI-DEFINITE
!                 = 3 IF NOP < 1
!                 = 4 IF NLOOP < 1

!     N.B. P, STEP AND VAR (IF IQUAD = 1) MUST HAVE DIMENSION AT LEAST NOP
!          IN THE CALLING PROGRAM.

!*****************************************************************************

INTEGER, INTENT(IN)        :: nop, maxfn, iprint, nloop, iquad
INTEGER, INTENT(OUT)       :: ifault
REAL (dp), INTENT(IN)      :: stopcr, simp
REAL (dp), INTENT(IN OUT)  :: p(:), step(:)
REAL (dp), INTENT(OUT)     :: var(:), func

!Additional information passed in, just to get to the function evaluator
integer, intent(in) :: n_dum
complex(dp), intent(in) :: g_dum(nImp,nImp,n_dum)
logical, intent(in) :: tMataxis_dum
logical, intent(in) :: tGrid
real(dp), intent(in) :: Freqs(n_dum)
real(dp), intent(in) :: Weights(n_dum)

INTERFACE
  SUBROUTINE functn(p, func, n_dum, nop, g_dum, tMataxis_dum, tGrid, Freqs, Weights, dJac)
    use const, only: dp
    use Globals, only: nImp
    IMPLICIT NONE
!    INTEGER, PARAMETER     :: dp = SELECTED_REAL_KIND(12, 60)
    integer, intent(in) :: n_dum,nop
    REAL (dp), INTENT(IN)  :: p(nop)
    REAL (dp), INTENT(OUT) :: func
    complex(dp), intent(in) :: g_dum(nImp,nImp,n_dum)
    logical, intent(in) :: tMataxis_dum
    logical, intent(in) :: tGrid
    real(dp), intent(in), optional :: Freqs(n_dum)
    real(dp), intent(in), optional :: Weights(n_dum)
    real(dp), intent(out), optional :: dJac(nop)
  END SUBROUTINE functn
END INTERFACE

!     Local variables

REAL (dp)   :: g(nop+1,nop), h(nop+1), pbar(nop), pstar(nop), pstst(nop), &
               aval(nop), pmin(nop), temp(nop), bmat(nop*(nop+1)/2),  &
               vc(nop*(nop+1)/2), ymin, rmax, hstst, a0, hmin, test, hmean, &
               hstd, hstar, hmax, savemn

REAL (dp), PARAMETER :: zero = 0._dp, half = 0.5_dp, one = 1._dp, two = 2._dp
INTEGER     :: i, i1, i2, iflag, ii, ij, imax, imin, irank, irow, j, j1, jj, &
               k, l, loop, nap, neval, nmore, np1, nullty

!     A = REFLECTION COEFFICIENT, B = CONTRACTION COEFFICIENT, AND
!     C = EXPANSION COEFFICIENT.

REAL (dp), PARAMETER :: a = 1._dp, b = 0.5_dp, c = 2._dp

!     SET LOUT = LOGICAL UNIT NO. FOR OUTPUT

INTEGER, PARAMETER :: lout = 6

!     IF PROGRESS REPORTS HAVE BEEN REQUESTED, PRINT HEADING

IF (iprint > 0) then
    WRITE (lout,5000) iprint
    call flush(lout)
endif

!     CHECK INPUT ARGUMENTS

ifault = 0
IF (nop <= 0) ifault = 3
IF (nloop <= 0) ifault = 4
IF (ifault /= 0) RETURN

!     SET NAP = NO. OF PARAMETERS TO BE VARIED, I.E. WITH STEP /= 0

nap = COUNT(step /= zero)
neval = 0
loop = 0
iflag = 0

!     IF NAP = 0 EVALUATE FUNCTION AT THE STARTING POINT AND RETURN

IF (nap <= 0) THEN
  CALL functn(p,func, n_dum, nop, g_dum, tMataxis_dum, tGrid, Freqs, Weights)
  RETURN
END IF

!     SET UP THE INITIAL SIMPLEX

20 g(1,:) = p
irow = 2
DO i = 1, nop
  IF (step(i) /= zero) THEN
    g(irow,:) = p
    g(irow,i) = p(i) + step(i)
    irow = irow + 1
  END IF
END DO

np1 = nap + 1
DO i = 1, np1
  p = g(i,:)
  CALL functn(p,h(i), n_dum,nop,g_dum,tMataxis_dum, tGrid, Freqs, Weights)
  neval = neval + 1
  IF (iprint > 0) THEN
    WRITE (lout,5100) neval, h(i), p
    call flush(lout)
  END IF
END DO

!     START OF MAIN CYCLE.

!     FIND MAX. & MIN. VALUES FOR CURRENT SIMPLEX (HMAX & HMIN).

Main_loop: DO
  loop = loop + 1
  imax = 1
  imin = 1
  hmax = h(1)
  hmin = h(1)
  DO i = 2, np1
    IF (h(i) > hmax) THEN
      imax = i
      hmax = h(i)
    ELSE
      IF (h(i) < hmin) THEN
        imin = i
        hmin = h(i)
      END IF
    END IF
  END DO

!     FIND THE CENTROID OF THE VERTICES OTHER THAN P(IMAX)

  pbar = zero
  DO i = 1, np1
    IF (i /= imax) THEN
      pbar = pbar + g(i,:)
    END IF
  END DO
  pbar = pbar / nap

!     REFLECT MAXIMUM THROUGH PBAR TO PSTAR,
!     HSTAR = FUNCTION VALUE AT PSTAR.

  pstar = a * (pbar - g(imax,:)) + pbar
  CALL functn(pstar,hstar,n_dum,nop,g_dum,tmataxis_dum, tGrid, Freqs, Weights)
  neval = neval + 1
  IF (iprint > 0) THEN
    IF (MOD(neval,iprint) == 0) then
        WRITE (lout,5100) neval, hstar, pstar
        call flush(lout)
    endif
  END IF

!     IF HSTAR < HMIN, REFLECT PBAR THROUGH PSTAR,
!     HSTST = FUNCTION VALUE AT PSTST.

  IF (hstar < hmin) THEN
    pstst = c * (pstar - pbar) + pbar
    CALL functn(pstst,hstst, n_dum,nop,g_dum,tMataxis_dum, tGrid, Freqs, Weights)
    neval = neval + 1
    IF (iprint > 0) THEN
      IF (MOD(neval,iprint) == 0) then
          WRITE (lout,5100) neval, hstst, pstst
          call flush(lout)
      endif
    END IF

!     IF HSTST < HMIN REPLACE CURRENT MAXIMUM POINT BY PSTST AND
!     HMAX BY HSTST, THEN TEST FOR CONVERGENCE.

    IF (hstst >= hmin) THEN   ! REPLACE MAXIMUM POINT BY PSTAR & H(IMAX) BY HSTAR.
      g(imax,:) = pstar
      h(imax) = hstar
    ELSE
      g(imax,:) = pstst
      h(imax) = hstst
    END IF
    GO TO 250
  END IF

!     HSTAR IS NOT < HMIN.
!     TEST WHETHER IT IS < FUNCTION VALUE AT SOME POINT OTHER THAN
!     P(IMAX).   IF IT IS REPLACE P(IMAX) BY PSTAR & HMAX BY HSTAR.

  DO i = 1, np1
    IF (i /= imax) THEN
      IF (hstar < h(i)) THEN  ! REPLACE MAXIMUM POINT BY PSTAR & H(IMAX) BY HSTAR.
        g(imax,:) = pstar
        h(imax) = hstar
        GO TO 250
      END IF
    END IF
  END DO

!     HSTAR > ALL FUNCTION VALUES EXCEPT POSSIBLY HMAX.
!     IF HSTAR <= HMAX, REPLACE P(IMAX) BY PSTAR & HMAX BY HSTAR.

  IF (hstar <= hmax) THEN
    g(imax,:) = pstar
    hmax = hstar
    h(imax) = hstar
  END IF

!     CONTRACTED STEP TO THE POINT PSTST,
!     HSTST = FUNCTION VALUE AT PSTST.

  pstst = b * g(imax,:) + (one-b) * pbar
  CALL functn(pstst,hstst, n_dum,nop,g_dum,tMataxis_dum, tGrid, Freqs, Weights)
  neval = neval + 1
  IF (iprint > 0) THEN
    IF (MOD(neval,iprint) == 0) then
        WRITE (lout,5100) neval, hstst, pstst
        call flush(lout)
    endif
  END IF

!     IF HSTST < HMAX REPLACE P(IMAX) BY PSTST & HMAX BY HSTST.

  IF (hstst <= hmax) THEN
    g(imax,:) = pstst
    h(imax) = hstst
    GO TO 250
  END IF

!     HSTST > HMAX.
!     SHRINK THE SIMPLEX BY REPLACING EACH POINT, OTHER THAN THE CURRENT
!     MINIMUM, BY A POINT MID-WAY BETWEEN ITS CURRENT POSITION AND THE
!     MINIMUM.

  DO i = 1, np1
    IF (i /= imin) THEN
      DO j = 1, nop
        IF (step(j) /= zero) g(i,j) = (g(i,j) + g(imin,j)) * half
        p(j) = g(i,j)
      END DO
      CALL functn(p,h(i), n_dum,nop,g_dum,tMataxis_dum, tGrid, Freqs, Weights)
      neval = neval + 1
      IF (iprint > 0) THEN
        IF (MOD(neval,iprint) == 0) then
            WRITE (lout,5100) neval, h(i), p
            call flush(lout)
        endif
      END IF
    END IF
  END DO

!     IF LOOP = NLOOP TEST FOR CONVERGENCE, OTHERWISE REPEAT MAIN CYCLE.

  250 IF (loop < nloop) CYCLE Main_loop

!     CALCULATE MEAN & STANDARD DEVIATION OF FUNCTION VALUES FOR THE
!     CURRENT SIMPLEX.

  hmean = SUM( h(1:np1) ) / np1
  hstd = SUM( (h(1:np1) - hmean) ** 2 )
  hstd = SQRT(hstd / np1)

!     IF THE RMS > STOPCR, SET IFLAG & LOOP TO ZERO AND GO TO THE
!     START OF THE MAIN CYCLE AGAIN.

  IF (hstd > stopcr .AND. neval <= maxfn) THEN
    iflag = 0
    loop = 0
    CYCLE Main_loop
  END IF

!     FIND THE CENTROID OF THE CURRENT SIMPLEX AND THE FUNCTION VALUE THERE.

  DO i = 1, nop
    IF (step(i) /= zero) THEN
      p(i) = SUM( g(1:np1,i) ) / np1
    END IF
  END DO
  CALL functn(p,func, n_dum,nop,g_dum,tMataxis_dum, tGrid, Freqs, Weights)
  neval = neval + 1
  IF (iprint > 0) THEN
    IF (MOD(neval,iprint) == 0) then
        WRITE (lout,5100) neval, func, p
        call flush(lout)
    endif
  END IF

!     TEST WHETHER THE NO. OF FUNCTION VALUES ALLOWED, maxfn, HAS BEEN
!     OVERRUN; IF SO, EXIT WITH IFAULT = 1.

  IF (neval > maxfn) THEN
    ifault = 1
    IF (iprint < 0) RETURN
    WRITE (lout,5200) maxfn
    WRITE (lout,5300) hstd
    WRITE (lout,5400) p
    WRITE (lout,5500) func
    call flush(lout)
    RETURN
  END IF

!     CONVERGENCE CRITERION SATISFIED.
!     IF IFLAG = 0, SET IFLAG & SAVE HMEAN.
!     IF IFLAG = 1 & CHANGE IN HMEAN <= STOPCR THEN SEARCH IS COMPLETE.

  IF (iprint >= 0) THEN
    WRITE (lout,5600)
    WRITE (lout,5400) p
    WRITE (lout,5500) func
    call flush(lout)
  END IF

  IF (iflag == 0 .OR. ABS(savemn-hmean) >= stopcr) THEN
    iflag = 1
    savemn = hmean
    loop = 0
  ELSE
    EXIT Main_loop
  END IF

END DO Main_loop

IF (iprint >= 0) THEN
  WRITE (lout,5700) neval
  WRITE (lout,5800) p
  WRITE (lout,5900) func
  call flush(lout)
END IF
IF (iquad <= 0) RETURN

!------------------------------------------------------------------

!     QUADRATIC SURFACE FITTING

IF (iprint >= 0) WRITE (lout,6000)

!     EXPAND THE FINAL SIMPLEX, IF NECESSARY, TO OVERCOME ROUNDING
!     ERRORS.

hmin = func
nmore = 0
DO i = 1, np1
  DO
    test = ABS(h(i)-func)
    IF (test < simp) THEN
      DO j = 1, nop
        IF (step(j) /= zero) g(i,j) = (g(i,j)-p(j)) + g(i,j)
        pstst(j) = g(i,j)
      END DO
      CALL functn(pstst,h(i), n_dum,nop,g_dum,tMataxis_dum, tGrid, Freqs, Weights)
      nmore = nmore + 1
      neval = neval + 1
      IF (h(i) >= hmin) CYCLE
      hmin = h(i)
      IF (iprint >= 0) WRITE (lout,5100) neval, hmin, pstst
    ELSE
      EXIT
    END IF
  END DO
END DO

!     FUNCTION VALUES ARE CALCULATED AT AN ADDITIONAL NAP POINTS.

DO i = 1, nap
  i1 = i + 1
  pstar = (g(1,:) + g(i1,:)) * half
  CALL functn(pstar,aval(i), n_dum, nop,g_dum,tMataxis_dum, tGrid, Freqs, Weights)
  nmore = nmore + 1
  neval = neval + 1
END DO

!     THE MATRIX OF ESTIMATED SECOND DERIVATIVES IS CALCULATED AND ITS
!     LOWER TRIANGLE STORED IN BMAT.

a0 = h(1)
DO i = 1, nap
  i1 = i - 1
  i2 = i + 1
  DO j = 1, i1
    j1 = j + 1
    pstst = (g(i2,:) + g(j1,:)) * half
    CALL functn(pstst,hstst, n_dum,nop,g_dum,tMataxis_dum, tGrid, Freqs, Weights)
    nmore = nmore + 1
    neval = neval + 1
    l = i * (i-1) / 2 + j
    bmat(l) = two * (hstst + a0 - aval(i) - aval(j))
  END DO
END DO

l = 0
DO i = 1, nap
  i1 = i + 1
  l = l + i
  bmat(l) = two * (h(i1) + a0 - two*aval(i))
END DO

!     THE VECTOR OF ESTIMATED FIRST DERIVATIVES IS CALCULATED AND
!     STORED IN AVAL.

DO i = 1, nap
  i1 = i + 1
  aval(i) = two * aval(i) - (h(i1) + 3._dp*a0) * half
END DO

!     THE MATRIX Q OF NELDER & MEAD IS CALCULATED AND STORED IN G.

pmin = g(1,:)
DO i = 1, nap
  i1 = i + 1
  g(i1,:) = g(i1,:) - g(1,:)
END DO

DO i = 1, nap
  i1 = i + 1
  g(i,:) = g(i1,:)
END DO

!     INVERT BMAT

CALL syminv(bmat, nap, bmat, temp, nullty, ifault, rmax)
IF (ifault == 0) THEN
  irank = nap - nullty
ELSE                                 ! BMAT not +ve definite
                                     ! Resume search for the minimum
  IF (iprint >= 0) WRITE (lout,6100)
  ifault = 2
  IF (neval > maxfn) RETURN
  WRITE (lout,6200)
  step = half * step
  GO TO 20
END IF

!     BMAT*A/2 IS CALCULATED AND STORED IN H.

DO i = 1, nap
  h(i) = zero
  DO j = 1, nap
    IF (j <= i) THEN
      l = i * (i-1) / 2 + j
    ELSE
      l = j * (j-1) / 2 + i
    END IF
    h(i) = h(i) + bmat(l) * aval(j)
  END DO
END DO

!     FIND THE POSITION, PMIN, & VALUE, YMIN, OF THE MINIMUM OF THE
!     QUADRATIC.

ymin = DOT_PRODUCT( h(1:nap), aval(1:nap) )
ymin = a0 - ymin
DO i = 1, nop
  pstst(i) = DOT_PRODUCT( h(1:nap), g(1:nap,i) )
END DO
pmin = pmin - pstst
IF (iprint >= 0) THEN
  WRITE (lout,6300) ymin, pmin
  WRITE (lout,6400)
  call flush(lout)
END IF

!     Q*BMAT*Q'/2 IS CALCULATED & ITS LOWER TRIANGLE STORED IN VC

DO i = 1, nop
  DO j = 1, nap
    h(j) = zero
    DO k = 1, nap
      IF (k <= j) THEN
        l = j * (j-1) / 2 + k
      ELSE
        l = k * (k-1) / 2 + j
      END IF
      h(j) = h(j) + bmat(l) * g(k,i) * half
    END DO
  END DO

  DO j = i, nop
    l = j * (j-1) / 2 + i
    vc(l) = DOT_PRODUCT( h(1:nap), g(1:nap,j) )
  END DO
END DO

!     THE DIAGONAL ELEMENTS OF VC ARE COPIED INTO VAR.

j = 0
DO i = 1, nop
  j = j + i
  var(i) = vc(j)
END DO
IF (iprint < 0) RETURN
WRITE (lout,6500) irank
CALL print_tri_matrix(vc, nop, lout)

WRITE (lout,6600)
CALL syminv(vc, nap, bmat, temp, nullty, ifault, rmax)

!     BMAT NOW CONTAINS THE INFORMATION MATRIX

WRITE (lout,6700)
CALL print_tri_matrix(bmat, nop, lout)

ii = 0
ij = 0
DO i = 1, nop
  ii = ii + i
  IF (vc(ii) > zero) THEN
    vc(ii) = one / SQRT(vc(ii))
  ELSE
    vc(ii) = zero
  END IF
  jj = 0
  DO j = 1, i - 1
    jj = jj + j
    ij = ij + 1
    vc(ij) = vc(ij) * vc(ii) * vc(jj)
  END DO
  ij = ij + 1
END DO

WRITE (lout,6800)
ii = 0
DO i = 1, nop
  ii = ii + i
  IF (vc(ii) /= zero) vc(ii) = one
END DO
CALL print_tri_matrix(vc, nop, lout)

!     Exit, on successful termination.

WRITE (lout,6900) nmore
RETURN

5000 FORMAT (' Progress Report every',i4,' function evaluations'/  &
             ' EVAL.   FUNC.VALUE.          PARAMETER VALUES')
5100 FORMAT (/' ', i4, '  ', g12.5, '  ', 5G11.4, 3(/t22, 5G11.4))
5200 FORMAT (' No. of function evaluations > ',i5)
5300 FORMAT (' RMS of function values of last simplex =', g14.6)
5400 FORMAT (' Centroid of last simplex =',4(/' ', 6G13.5))
5500 FORMAT (' Function value at centroid =', g14.6)
5600 FORMAT (/' EVIDENCE OF CONVERGENCE')
5700 FORMAT (/' Minimum found after',i5,' function evaluations')
5800 FORMAT (' Minimum at',4(/' ', 6G13.6))
5900 FORMAT (' Function value at minimum =', g14.6)
6000 FORMAT (/' Fitting quadratic surface about supposed minimum'/)
6100 FORMAT (/' MATRIX OF ESTIMATED SECOND DERIVATIVES NOT +VE DEFN.'/ &
             ' MINIMUM PROBABLY NOT FOUND'/)
6200 FORMAT (/t11, 'Search restarting'/)
6300 FORMAT (' Minimum of quadratic surface =',g14.6,' at',4(/' ', 6G13.5))
6400 FORMAT (' IF THIS DIFFERS BY MUCH FROM THE MINIMUM ESTIMATED ',      &
             'FROM THE MINIMIZATION,'/' THE MINIMUM MAY BE FALSE &/OR THE ',  &
             'INFORMATION MATRIX MAY BE INACCURATE'/)
6500 FORMAT (' Rank of information matrix =',i3/ &
             ' Inverse of information matrix:-')
6600 FORMAT (/' If the function minimized was -LOG(LIKELIHOOD),'/         &
             ' this is the covariance matrix of the parameters.'/         &
             ' If the function was a sum of squares of residuals,'/       &
             ' this matrix must be multiplied by twice the estimated ',   &
             'residual variance'/' to obtain the covariance matrix.'/)
6700 FORMAT (' INFORMATION MATRIX:-'/)
6800 FORMAT (/' CORRELATION MATRIX:-')
6900 FORMAT (/' A further',i4,' function evaluations have been used'/)

END SUBROUTINE minim




SUBROUTINE syminv(a, n, c, w, nullty, ifault, rmax)

!     ALGORITHM AS7, APPLIED STATISTICS, VOL.17, 1968.

!     ARGUMENTS:-
!     A()    = INPUT, THE SYMMETRIC MATRIX TO BE INVERTED, STORED IN
!                LOWER TRIANGULAR FORM
!     N      = INPUT, ORDER OF THE MATRIX
!     C()    = OUTPUT, THE INVERSE OF A (A GENERALIZED INVERSE IF C IS
!                SINGULAR), ALSO STORED IN LOWER TRIANGULAR FORM.
!                C AND A MAY OCCUPY THE SAME LOCATIONS.
!     W()    = WORKSPACE, DIMENSION AT LEAST N.
!     NULLTY = OUTPUT, THE RANK DEFICIENCY OF A.
!     IFAULT = OUTPUT, ERROR INDICATOR
!                 = 1 IF N < 1
!                 = 2 IF A IS NOT +VE SEMI-DEFINITE
!                 = 0 OTHERWISE
!     RMAX   = OUTPUT, APPROXIMATE BOUND ON THE ACCURACY OF THE DIAGONAL
!                ELEMENTS OF C.  E.G. IF RMAX = 1.E-04 THEN THE DIAGONAL
!                ELEMENTS OF C WILL BE ACCURATE TO ABOUT 4 DEC. DIGITS.

!     LATEST REVISION - 1 April 1985

!***************************************************************************

REAL (dp), INTENT(IN OUT) :: a(:), c(:), w(:)
INTEGER, INTENT(IN)       :: n
INTEGER, INTENT(OUT)      :: nullty, ifault
REAL (dp), INTENT(OUT)    :: rmax

REAL (dp), PARAMETER :: zero = 0._dp, one = 1._dp
INTEGER              :: i, icol, irow, j, jcol, k, l, mdiag, ndiag, nn, nrow
REAL (dp)            :: x

nrow = n
ifault = 1
IF (nrow > 0) THEN
  ifault = 0

!     CHOLESKY FACTORIZATION OF A, RESULT IN C

  CALL chola(a, nrow, c, nullty, ifault, rmax, w)
  IF (ifault == 0) THEN

!     INVERT C & FORM THE PRODUCT (CINV)'*CINV, WHERE CINV IS THE INVERSE
!     OF C, ROW BY ROW STARTING WITH THE LAST ROW.
!     IROW = THE ROW NUMBER, NDIAG = LOCATION OF LAST ELEMENT IN THE ROW.

    nn = nrow * (nrow+1) / 2
    irow = nrow
    ndiag = nn
    10 IF (c(ndiag) /= zero) THEN
      l = ndiag
      DO i = irow, nrow
        w(i) = c(l)
        l = l + i
      END DO
      icol = nrow
      jcol = nn
      mdiag = nn

      30 l = jcol
      x = zero
      IF (icol == irow) x = one / w(irow)
      k = nrow
      40 IF (k /= irow) THEN
        x = x - w(k) * c(l)
        k = k - 1
        l = l - 1
        IF (l > mdiag) l = l - k + 1
        GO TO 40
      END IF

      c(l) = x / w(irow)
      IF (icol == irow) GO TO 60
      mdiag = mdiag - icol
      icol = icol - 1
      jcol = jcol - 1
      GO TO 30
    END IF ! (c(ndiag) /= zero)

    l = ndiag
    DO j = irow, nrow
      c(l) = zero
      l = l + j
    END DO

    60 ndiag = ndiag - irow
    irow = irow - 1
    IF (irow /= 0) GO TO 10
  END IF
END IF
RETURN

END SUBROUTINE syminv



SUBROUTINE chola(a, n, u, nullty, ifault, rmax, r)

!     ALGORITHM AS6, APPLIED STATISTICS, VOL.17, 1968, WITH
!     MODIFICATIONS BY A.J.MILLER

!     ARGUMENTS:-
!     A()    = INPUT, A +VE DEFINITE MATRIX STORED IN LOWER-TRIANGULAR
!                FORM.
!     N      = INPUT, THE ORDER OF A
!     U()    = OUTPUT, A LOWER TRIANGULAR MATRIX SUCH THAT U*U' = A.
!                A & U MAY OCCUPY THE SAME LOCATIONS.
!     NULLTY = OUTPUT, THE RANK DEFICIENCY OF A.
!     IFAULT = OUTPUT, ERROR INDICATOR
!                 = 1 IF N < 1
!                 = 2 IF A IS NOT +VE SEMI-DEFINITE
!                 = 0 OTHERWISE
!     RMAX   = OUTPUT, AN ESTIMATE OF THE RELATIVE ACCURACY OF THE
!                DIAGONAL ELEMENTS OF U.
!     R()    = OUTPUT, ARRAY CONTAINING BOUNDS ON THE RELATIVE ACCURACY
!                OF EACH DIAGONAL ELEMENT OF U.

!     LATEST REVISION - 1 April 1985

!***************************************************************************

REAL (dp), INTENT(IN)   :: a(:)
INTEGER, INTENT(IN)     :: n
REAL (dp), INTENT(OUT)  :: u(:)
INTEGER, INTENT(OUT)    :: nullty, ifault
REAL (dp), INTENT(OUT)  :: rmax, r(:)

!     ETA SHOULD BE SET EQUAL TO THE SMALLEST +VE VALUE SUCH THAT
!     1._dp + ETA IS CALCULATED AS BEING GREATER THAN 1._dp IN THE ACCURACY
!     BEING USED.

REAL (dp), PARAMETER :: eta = EPSILON(1.0_dp), zero = 0._dp
INTEGER              :: i, icol, irow, j, k, l, m
REAL (dp)            :: rsq, w

ifault = 1
IF (n > 0) THEN
  ifault = 2
  nullty = 0
  rmax = eta
  r(1) = eta
  j = 1
  k = 0

!     FACTORIZE COLUMN BY COLUMN, ICOL = COLUMN NO.

  DO  icol = 1, n
    l = 0

!     IROW = ROW NUMBER WITHIN COLUMN ICOL

    DO  irow = 1, icol
      k = k + 1
      w = a(k)
      IF (irow == icol) rsq = (w*eta) ** 2
      m = j
      DO  i = 1, irow
        l = l + 1
        IF (i == irow) EXIT
        w = w - u(l) * u(m)
        IF (irow == icol) rsq = rsq + (u(l)**2*r(i)) ** 2
        m = m + 1
      END DO

      IF (irow == icol) EXIT
      IF (u(l) /= zero) THEN
        u(k) = w / u(l)
      ELSE
        u(k) = zero
        IF (ABS(w) > ABS(rmax*a(k))) GO TO 60
      END IF
    END DO

!     END OF ROW, ESTIMATE RELATIVE ACCURACY OF DIAGONAL ELEMENT.

    rsq = SQRT(rsq)
    IF (ABS(w) > 5.*rsq) THEN
      IF (w < zero) RETURN
      u(k) = SQRT(w)
      r(i) = rsq / w
      IF (r(i) > rmax) rmax = r(i)
    ELSE
      u(k) = zero
      nullty = nullty + 1
    END IF

    j = j + icol
  END DO
  ifault = zero
END IF
60 RETURN

END SUBROUTINE chola



SUBROUTINE print_tri_matrix(a, n, lout)

INTEGER, INTENT(IN)    :: n, lout
REAL (dp), INTENT(IN)  :: a(:)

!     Local variables
INTEGER  :: i, ii, i1, i2, l

l = 1
DO l = 1, n, 6
  ii = l * (l-1) / 2
  DO i = l, n
    i1 = ii + l
    ii = ii + i
    i2 = MIN(ii,i1+5)
    WRITE (lout,'(1X,6G13.5)') a(i1:i2)
  END DO
END DO
RETURN

END SUBROUTINE print_tri_matrix

END MODULE MinAlgos    

!PROGRAM tminim
!
!! Use minim to maximize the objective function:
!
!! 1/{1 + (x-y)^2} + sin(pi.y.z/2) + exp[-{(x+z)/y - 2}^2]
!
!! with respect to x, y and z.
!! We will actually minimize its negative.
!! The maximum occurs at x = y = z = +/-sqrt(4n+1) for any integer n.
!
!USE Nelder_Mead
!IMPLICIT NONE
!INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(14, 60)
!REAL (dp)          :: object, p(3), simp, step(3), stopcr, var(3)
!INTEGER            :: ier, iprint, iquad, maxf, nloop, nop
!LOGICAL            :: first
!
!INTERFACE
!  SUBROUTINE objfun(p, func)
!    IMPLICIT NONE
!    INTEGER, PARAMETER     :: dp = SELECTED_REAL_KIND(14, 60)
!    REAL (dp), INTENT(IN)  :: p(:)
!    REAL (dp), INTENT(OUT) :: func
!  END SUBROUTINE objfun
!END INTERFACE
!
!! Set up starting values & step sizes.
!! Parameters p(1), p(2), p(3) are x, y & z.
!
!p(1) = 0._dp
!p(2) = 1._dp
!p(3) = 2._dp
!step = 0.4_dp
!nop = 3
!
!! Set max. no. of function evaluations = 250, print every 10.
!
!maxf = 250
!iprint = 10
!
!! Set value for stopping criterion.   Stopping occurs when the
!! standard deviation of the values of the objective function at
!! the points of the current simplex < stopcr.
!
!stopcr = 1.d-04
!nloop = 6
!
!! Fit a quadratic surface to be sure a minimum has been found.
!
!iquad = 1
!
!! As function value is being evaluated in REAL (dp), it
!! should be accurate to about 15 decimals.   If we set simp = 1.d-6,
!! we should get about 9 dec. digits accuracy in fitting the surface.
!
!simp = 1.d-6
!
!! Now call MINIM to do the work.
!
!first = .true.
!DO
!  CALL minim(p, step, nop, object, maxf, iprint, stopcr, nloop,   &
!             iquad, simp, var, objfun, ier)
!
!! If ier > 0, try a few more function evaluations.
!
!  IF (ier == 0) EXIT
!  IF (.NOT. first) STOP
!  first = .false.
!  maxf = 100
!END DO
!
!! Successful termination.
!
!WRITE(*, 900) object, p
!900 FORMAT(' Success !'/' Objective function = ', f12.6/ ' at: ', 3f12.6)
!WRITE(*, 910) var
!910 FORMAT(' Elements of var = ', 3f12.6)
!
!STOP
!
!END PROGRAM tminim
!
!
!
!SUBROUTINE objfun(p, func)
!
!! This is the subroutine which the user must write.
!! Remember that we are minimizing the negative of the function we
!! really want to maximize.
!
!IMPLICIT NONE
!INTEGER, PARAMETER     :: dp = SELECTED_REAL_KIND(14, 60)
!REAL (dp), INTENT(IN)  :: p(:)
!REAL (dp), INTENT(OUT) :: func
!
!!     Local variables
!
!REAL (dp), PARAMETER :: half = 0.5_dp, one = 1._dp, pi = 3.14159265357989_dp, &
!                        two = 2._dp
!REAL (dp)            :: x, y, z
!
!x = p(1)
!y = p(2)
!z = p(3)
!func = -one/(one + (x-y)**2) - SIN(half*pi*y*z) - EXP(-((x+z)/y - two)**2)
!RETURN
!END SUBROUTINE objfun
