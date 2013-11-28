module Continuation
    use const
    use errors, only: stop_all,warning
    use globals
    implicit none

    contains

        !A simple function (from numerical recipies!) for rational function interpolation
        !Fit the imaginary time data to a rational function, which theoretically can capture poles on the real axis.
        !See Journal of Low Temperature Physcis, Vol 29, No 3/4, page 179 (1977)

        !Given arrays xa and ya, each of length n, and given a value of x, this routine returns a value of y and an accuracy estimate dy.
        !The value returned is that of the diagonal rational function, evaluated at x, which passes through the n points (xa_i,ya_i), i=1,n
        subroutine ratint(xa,ya,n,x,y,dy)
            integer, intent(in) :: n
            complex(dp), intent(in) :: xa(n),ya(n),x
            complex(dp), intent(out) :: y,dy
            integer, parameter :: nmax=10000    !The largest expected value of n
            integer :: i,m,ns
            real(dp) :: hh,h_r
            complex(dp) :: dd,h,t,w
            complex(dp), allocatable :: c(:),d(:)
            character(len=*), parameter :: t_r='ratint'

            allocate(c(n))
            allocate(d(n))
            c(:) = zzero
            d(:) = zzero

            ns = 1
            hh = abs(x-xa(1))
            do i = 1,n
                h_r = real(abs(x-xa(i)),dp)
                if (h_r.eq.zero) then
                    y = ya(i)
                    dy = zzero
                    return
                elseif (h_r.lt.hh) then
                    ns = i
                    hh = h
                endif
                c(i) = ya(i)
                d(i) = ya(i) + tiny(h_r)
            enddo
            y = ya(ns)
            ns = ns - 1
            do m = 1,n-1
                do i = 1,n-m
                    w = c(i+1)-d(i)
                    h = xa(i+m)-x
                    t = (xa(i)-x)*d(i)/h
                    dd = t - c(i+1)
                    if(abs(dd).eq.zero) call stop_all(t_r,'Interpolating function has a plot at the requested frequency')
                    dd = w / dd
                    d(i) = c(i+1) * dd
                    c(i) = t * dd
                enddo
                if(2*ns.lt.n-m) then
                    dy = c(ns+1)
                else
                    dy = d(ns)
                    ns = ns - 1
                endif
                y = y + dy
            enddo

            deallocate(c,d)
        end subroutine ratint


!!       initially take tau to simply be the half-life  tau_step = log(2.0_dp)/End_Omega        !Assume that the spacing is given by the half-life of the highest frequency excitation.
!        !A very simple a naive Laplace transform to convert from frequency domain data to imaginary time data.
!        !This uses simple quadrature for the integration. I'm sure it can be done better!
!        subroutine LaplaceTrans(n_w,F_w,n_tau,tau_step,tau_end,F_tau)
!            integer, intent(in) :: n_w                !Number of frequency points.
!            complex(dp), intent(in) :: F_w(n_w)       !Function in frequency domain. (i.e. the greens function)
!            integer, intent(in) :: n_tau
!            real(dp), intent(in) :: tau_step, tau_end
!            real(dp), intent(out) :: F_tau(n_tau)    !Function returned in imaginary time domain.
!
!            integer :: i
!            real(dp) :: tau
!            character(len=*), parameter :: t_r='LaplaceTrans'
!
!            tau = 0.0_dp
!            i = 0
!            do while(tau.lt.tau_end) 
!                i = i + 1
!                if(i.gt.n_tau) call stop_all(t_r,'Too many imaginary time slices')
!
!                tau = tau + tau_step
!            enddo
!
!
!        end subroutine LaplaceTrans

end module Continuation
