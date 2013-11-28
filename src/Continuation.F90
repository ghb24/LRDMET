module Continuation
    use const
    use errors, only: stop_all,warning
    use globals
    use LinearResponse, only: GetNextOmega
    implicit none

    contains

        subroutine FromMatsubaraToRealFreq(n_Im,n_Re,Fn_Im,Fn_Re,Fn_Re_Err)
            implicit none
            integer, intent(in) :: n_Im,n_Re
            complex(dp), intent(in) :: Fn_Im(nImp,nImp,n_Im)
            complex(dp), intent(out) :: Fn_Re(nImp,nImp,n_Re)
            complex(dp), intent(out) :: Fn_Re_Err(nImp,nImp,n_Re)

            integer :: el_1,el_2,i,OmegaVal
            real(dp) :: Omega
            complex(dp) :: RealFreqPoint
            complex(dp), allocatable :: ImData(:),ImFreqPoints(:)

            Fn_Re(:,:,:) = zzero
            Fn_Re_Err(:,:,:) = zzero

            allocate(ImData(n_Im))
            allocate(ImFreqPoints(n_Im))
            !Fill array with matsubara frequency points
            ImFreqPoints(:) = zzero
            OmegaVal = 0
            do while(.true.)
                call GetNextOmega(Omega,OmegaVal,tMatbrAxis=.true.)
                if(OmegaVal.lt.0) exit
                ImFreqPoints(OmegaVal) = dcmplx(zero,Omega)
            enddo

            do el_1 = 1,nImp
                do el_2 = el_1,nImp

                    ImData(:) = Fn_Im(el_1,el_2,:)

                    OmegaVal = 0
                    do while(.true.)
                        !Cycle through the *real* frequency values
                        call GetNextOmega(Omega,OmegaVal,tMatbrAxis=.false.)
                        if(OmegaVal.lt.0) exit

                        !Which real frequency point do we want
                        RealFreqPoint = dcmplx(Omega,dDelta)

                        !Simple fit to Pade approximant
                        call ratint(ImFreqPoints,ImData,n_Im,RealFreqPoint,Fn_Re(el_1,el_2,OmegaVal),Fn_Re_Err(el_1,el_2,OmegaVal))

                    enddo

                    if(el_1.ne.el_2) then
                        !Copy the data to the other side to ensure off-diagonal hermiticity
                        do i = 1,n_Re
                            Fn_Re(el_2,el_1,i) = dconjg(Fn_Re(el_1,el_2,i))
                            Fn_Re_Err(el_2,el_1,i) = dconjg(Fn_Re_Err(el_1,el_2,i))
                        enddo
                    endif

                enddo
            enddo

            deallocate(ImFreqPoints,ImData)

        end subroutine FromMatsubaraToRealFreq

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
