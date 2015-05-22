program matsusum
    implicit none
    integer    :: n=1000
    complex(8) :: zsum,matsu,ztmp,MatsuGF,zsum2
    integer    :: k,l,m,ntau,i
    real(8)    :: tau,scal,mu,beta,epsk,pi,tauval

    zsum=0.
    mu=1.0
    epsk=1.0
    tau=0.001d0
    beta=1000.d0

    ntau = int(beta/tau)

    pi=dacos(-1.d0)

    write(6,*) "Nyquist frequency: ",pi*n/beta

!    write(*,*) 'PI : ', pi
    open(77,file='MatsuGF',status='unknown')

    do k = -n,n
        matsu = CMPLX(0.d0,(pi*(2.0*k-1.0))/beta)
        MatsuGF = 1.d0 / (matsu + cmplx(mu - epsk,0.) ) 
        write(77,*) aimag(matsu),real(MatsuGF),aimag(MatsuGF)
    enddo
    close(77)

    open(77,file='Naive_TauGF',status='unknown')
    do i = 0,ntau
        tauval = i*tau

        zsum = cmplx(0.0,0.0)
        zsum2 = cmplx(0.0,0.0)
        do k = -n,n
            matsu = CMPLX(0.d0,(pi*(2.0*k-1.0))/beta)
            MatsuGF = 1.d0 / (matsu + cmplx(mu - epsk,0.) ) 
            zsum = zsum + MatsuGF*exp(-matsu*tauval)
            zsum2 = zsum2 + (MatsuGF - (1.0/matsu))*exp(-matsu*tauval)
        enddo
        zsum = zsum / beta
        zsum2 = (zsum2 / beta) - 0.5    !Remove the tail from the discrete sum,
                                        !and analytically FT it

        write(77,*) tauval,real(zsum),aimag(zsum),real(zsum2),aimag(zsum2)
    enddo
    close(77)



!     scal=abs(matsu)*tau

!     ztmp=MatsuGF / beta

!     ztmp = ztmp - 1.0/matsu

!     zsum = zsum + ztmp*exp(cmplx(0.0,matsu*tau))   !cmplx(cos(scal),sin(scal))

!     write(40,*) k,real(ztmp),aimag(ztmp)
!     write(41,*) k,real(zsum),aimag(zsum)

!    end do

!    write(*,*)'summation : ', zsum + 0.5
!    write(*,*)'FERMI     : ', 1.d0/(exp(beta*(epsk-mu))+1.d0)


end program

