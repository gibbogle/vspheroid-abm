! To determine checkpoint duration Tcp as a function of initial lesion number n0.
! The formulation comes from Curtis (1986).
! Tcp = Kcp.AUC(n0,Tcp)
! where AUC(n0,T1) = area under curve of n(t) from t=0 to t=T1,
! and:
! n(t) = (n0.e^(-et))/(1 + n0.(1-e^(-et))/e)
!
! The parameters are e = ePL and Kcp.
! Use time in hours.

module Tcp_mod

use real_kind_mod

implicit none

real(REAL_KIND) :: ePL, e, Kcp

contains

!--------------------------------------------------------------------------
! Time unit = hour, convert Tcp to seconds
! Note: Krepair, Kmisrepair are Curtis's epsilon_PL, epsilon_2PL
!--------------------------------------------------------------------------
subroutine makeTCP(Tcp, np, epsilon_PL, epsilon_2PL, Kcp_set)
real(REAL_KIND) :: Tcp(0:np)
integer :: np
real(REAL_KIND) :: epsilon_PL, epsilon_2PL(2), Kcp_set
real(REAL_KIND) :: n0, T1, auc, epsilon_ratio
integer :: k, i, nb
real(REAL_KIND) :: Tb, alpha

ePL = epsilon_PL    ! /h
e = epsilon_PL/sum(epsilon_2PL)
Kcp = Kcp_set

Tcp(0) = 0
do i = 1,np
    n0 = i
    !write(*,'(a,f6.0)') 'n0: ',n0
    do k = 1,10000
        T1 = k*0.01
        auc = getAUC(n0,T1)
!        write(*,'(a,3f8.3)') 'T1, auc: ',T1,auc,Kcp*auc
        if (T1 > Kcp*auc) exit
    enddo
    Tcp(i) = T1
!    write(*,'(a,f6.0,f6.2)') 'n0, Tcp: ',n0,T1
enddo
! Blend with linear variation
nb = 2/Kcp + 2
Tb = Tcp(nb)
do i = 1,nb
    alpha = real(i)/nb
    Tcp(i) = alpha*Tb
enddo
Tcp = 3600*Tcp  ! hours -> seconds

end subroutine

!--------------------------------------------------------------------------
! Lumping 3 kinds of lethal damage together
!--------------------------------------------------------------------------
function getAUC(n0, T1) result(area)
real(REAL_KIND) :: n0, T1, area
real(REAL_KIND) :: dt, t, n, np, et, darea
integer :: it, nt

dt = 0.01
nt = T1/dt + 0.5
area = 0
np = n0
do it = 1,nt
    t = it*dt
    et = exp(-ePL*t)
    n = n0*et/(1 + n0*(1-et)/e)
    darea = dt*(n + np)/2.
    area = area + darea
!    write(*,*) it,darea,area
    np = n
enddo
end function

end module

!program main
!use Tcp_mod
!implicit none
!
!real(REAL_KIND) :: n0, T1, auc
!integer :: k, i, np, nb
!real(REAL_KIND) :: Tb, alpha, Tcp(0:1000)
!
!np = 100
!ePL = 0.5     ! /h
!e = 9
!Kcp = 0.13
!dt = 0.001
!
!do i = 1,np
!    n0 = i
!    !write(*,'(a,f6.0)') 'n0: ',n0
!    do k = 1,10000
!        T1 = k*0.01
!        auc = getAUC(n0,T1)
!!        write(*,'(a,3f8.3)') 'T1, auc: ',T1,auc,Kcp*auc
!        if (T1 > Kcp*auc) exit
!    enddo
!    Tcp(i) = T1
!!    write(*,'(a,f6.0,f6.2)') 'n0, Tcp: ',n0,T1
!enddo
!! Blend with linear variation
!nb = 2/Kcp + 2
!Tb = Tcp(nb)
!do i = 1,nb
!    alpha = real(i)/nb
!    Tcp(i) = alpha*Tb
!enddo
!Tcp(0) = 0
!do i = 0,np
!    write(*,'(a,i3,f6.2)') 'blended: i, Tcp: ',i,Tcp(i)
!enddo
!end program
