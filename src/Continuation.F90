module Continuation
    use const
    use errors, only: stop_all,warning
    use globals
    implicit none

    contains

!       initially take tau to simply be the half-life  tau_step = log(2.0_dp)/End_Omega        !Assume that the spacing is given by the half-life of the highest frequency excitation.
        !A very simple a naive Laplace transform to convert from frequency domain data to imaginary time data.
        !This uses simple quadrature for the integration. I'm sure it can be done better!
        subroutine LaplaceTrans(n_w,F_w,n_tau,tau_step,tau_end,F_tau)
            integer, intent(in) :: n_w                !Number of frequency points.
            complex(dp), intent(in) :: F_w(n_w)       !Function in frequency domain. (i.e. the greens function)
            integer, intent(in) :: n_tau
            real(dp), intent(in) :: tau_step, tau_end
            real(dp), intent(out) :: F_tau(n_tau)    !Function returned in imaginary time domain.

            integer :: i
            real(dp) :: tau
            character(len=*), parameter :: t_r='LaplaceTrans'

            tau = 0.0_dp
            i = 0
            do while(tau.lt.tau_end) 
                i = i + 1
                if(i.gt.n_tau) call stop_all(t_r,'Too many imaginary time slices')

                tau = tau + tau_step
            enddo


        end subroutine LaplaceTrans

end module Continuation
