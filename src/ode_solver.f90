! To use ODE solver on OGL metabolism
!----------------------------------------------------------------------------------
! Note: The value of spcrad was first determined by writing out the value computed in rkc.
! Later it was just determined by trial, then made into a run parameter.
!----------------------------------------------------------------------------------
double precision function spcrad(neqn,t,y)
!DEC$ ATTRIBUTES DLLEXPORT :: spcrad
use global
integer :: neqn
double precision :: t, y(neqn)
spcrad = spcrad_value
end function

!----------------------------------------------------------------------------------
!----------------------------------------------------------------------------------
module ode_solver

use chemokine
use metabolism
use rkc_90

implicit none

!real(REAL_KIND) :: work_rkc(8+5*3*2)
logical :: chemo_active(2*MAX_CHEMO)    ! flags necessity to solve for the constituent

!type(metabolism_type), pointer :: mp_rkc
!real(REAL_KIND) :: Cex_rkc(3)

! These global variables are a problem for parallelisation of f_rkc !!!!!!!!!!!!!!!!!!!

!integer, parameter :: N_EC_O = 20
!integer, parameter :: N_EC_G = 20
!integer, parameter :: N_EC_L = 20
!real(REAL_KIND) :: EC_O_array(N_EC_O)
!real(REAL_KIND) :: EC_G_array(N_EC_G)
!real(REAL_KIND) :: EC_L_array(N_EC_L)
!
!real(REAL_KIND),allocatable :: IC_O_lookup(:,:,:)
!real(REAL_KIND),allocatable :: IC_G_lookup(:,:,:)
!real(REAL_KIND),allocatable :: IC_L_lookup(:,:,:)
!real(REAL_KIND),allocatable :: dMdt_O_lookup(:,:,:)
!real(REAL_KIND),allocatable :: dMdt_G_lookup(:,:,:)
!real(REAL_KIND),allocatable :: dMdt_L_lookup(:,:,:)

contains

!----------------------------------------------------------------------------------
!----------------------------------------------------------------------------------
subroutine OGLSolver(kcell,tstart,dt,ok)
integer :: kcell
real(REAL_KIND) :: tstart, dt
logical :: ok
integer :: ichemo, k, ict, neqn, i
real(REAL_KIND) :: t, tend
real(REAL_KIND) :: Cin(3)
real(REAL_KIND) :: timer1, timer2
type(cell_type), pointer :: cp
! Variables for RKC
integer :: info(4), idid
real(REAL_KIND) :: rtol, atol(1)
type(rkc_comm) :: comm_rkc(1)
real(REAL_KIND) :: work_rkc(8+5*3*2)

cp => cell_list(kcell)
!if (kcell == 8000) write(*,'(a,i6,6f8.4)') 'OGLSolver: Cin,Cex: ',kcell,cp%Cin(1:3),cp%Cex(1:3)
kcell_now = kcell
! Try using icase for kcell
Cin(1:3) = cp%Cin(1:3)
neqn = 3

!mp_rkc => cp%metab
!Cex_rkc = cp%Cex(1:3)

info(1) = 1
info(2) = 1		! = 1 => use spcrad() to estimate spectral radius, != 1 => let rkc do it
info(3) = 1
info(4) = 0
rtol = 1d-2
atol = rtol

idid = 0
t = tstart
tend = t + dt
call rkc(comm_rkc(1),neqn,f_rkc_OGL,Cin,t,tend,rtol,atol,info,work_rkc,idid,kcell)	! ict)
if (idid /= 1) then
	write(logmsg,*) 'Solver: Failed at t = ',t,' with idid = ',idid
	call logger(logmsg)
	ok = .false.
	return
endif
cp%Cin(1:3) = Cin

!call get_metab_rates(ict,mp,Cin)
!if (kcell_now == 1) then
!	write(*,'(a,6f8.4)') 'after OGLSolver: Cin,Cex: ',cp%Cin(1:3),cp%Cex(1:3)
!	write(*,'(a,3e12.3)') 'rates: ',mp%O_rate,mp%G_rate,mp%L_rate
!endif
end subroutine

!----------------------------------------------------------------------------------
!----------------------------------------------------------------------------------
subroutine f_rkc_OGL(neqn,t,y,dydt,icase)
integer :: neqn, icase
real(REAL_KIND) :: t, y(neqn), dydt(neqn)
integer :: k, kk, i, ichemo, ict, Nlivecells
real(REAL_KIND) :: dCsum, dCdiff, dCreact, vol_cm3, Cex, Cin(3)
real(REAL_KIND) :: C, membrane_kin, membrane_kout, membrane_flux, area_factor, Cbnd, rate
type(cell_type), pointer :: cp
type(metabolism_type), pointer :: mp
!real(REAL_KIND) :: A, d, dX, dV, Kd, KdAVX
!type(metabolism_type), pointer :: mp
real(REAL_KIND) :: average_volume = 1.2
logical :: use_average_volume = .true.

!if (kcell_now == kcell_test) write(*,'(a,3f9.6)') 'f_rkc_OGL: y: ',y(1:3)

!ict = icase
cp => cell_list(icase)
mp => cp%metab
ict = cp%celltype
area_factor = (average_volume)**(2./3.)
vol_cm3 = Vcell_cm3		!!!!!!!!!!!!!! temporary !!!!!!!!!!!!
!call get_metab_rates(ict,mp_rkc,y)
call get_metab_rates(ict,mp,y)
!if (dbug) write(*,'(a,3e12.3)') 'rates: ',mp_rkc%O_rate,mp_rkc%G_rate,mp_rkc%L_rate
do ichemo = 1,3
	! First process IC reactions
	k = ichemo
	C = y(k)
    membrane_kin = chemo(ichemo)%membrane_diff_in
    membrane_kout = chemo(ichemo)%membrane_diff_out
!	membrane_flux = area_factor*(membrane_kin*Cex_rkc(ichemo) - membrane_kout*C)
	membrane_flux = area_factor*(membrane_kin*cp%Cex(ichemo) - membrane_kout*C)
    if (ichemo == OXYGEN) then
		rate = mp%O_rate
    elseif (ichemo == GLUCOSE) then	! 
		rate = mp%G_rate
    elseif (ichemo == LACTATE) then
		rate = -mp%L_rate
    endif
	dCreact = (membrane_flux - rate)/vol_cm3
	dydt(k) = dCreact
	if (dbug) write(*,'(a,i4,5e12.3)') 'dydt: ',ichemo,dCreact,membrane_flux,rate,C,cp%Cex(ichemo)
	if (isnan(dydt(k))) then
		write(nflog,*) 'f_rkc_OGL: dydt isnan: ',ichemo,dydt(k)
		write(*,*) 'f_rkc_OGL: dydt isnan: ',ichemo,dydt(k)
		stop
	endif	
enddo
end subroutine

subroutine testOGL
real(REAL_KIND) :: dt
integer :: ichemo, it, ityp, kcell
real(REAL_KIND) :: Kin, Kout, tstart, area_factor, CexG, CinG, CexO, CinO, H, fPDK, G, a_G, a_O, b_G, b_O
type(cell_type), pointer :: cp
logical :: ok

kcell = 1
cp => cell_list(kcell)
ityp = cp%celltype
tstart = 0
area_factor = (1.2)**(2./3.)
cp%Cex(1:3) = [0.1,2.5,5.0]
cp%Cin(1:3) = cp%Cex(1:3)
cp%metab%HIF1 = 0.3
cp%metab%PDK1 = 1.0
do it = 1,10
	dt = it*1
	call OGLSolver(kcell,tstart,dt,ok)	
	do ichemo = OXYGEN,LACTATE
		Kin = chemo(ichemo)%membrane_diff_in
		Kout = chemo(ichemo)%membrane_diff_out
		cp%dMdt(ichemo) = area_factor*(Kin*cp%Cex(ichemo) - Kout*cp%Cin(ichemo))
	enddo
	write(*,'(a,f6.1,3f8.5,3e11.3)') 'Cin: ',dt,cp%Cin(1:3),cp%dMdt(1:3)
enddo
! Now vary Cex to account for H
H = cp%metab%HIF1
G = cp%Cex(2)
a_G = (1 + K_Hb(ityp)*H)* G**Hill_N_G / (G**Hill_N_G + Hill_Km_G**Hill_N_G)
fPDK = cp%metab%PDK1
a_O = fPDK*cp%Cex(1)/(Hill_Km_O2 + cp%Cex(1))
CexG = (a_G/(1-a_G))*Hill_Km_G
CexO = (a_O/(1-a_O))*Hill_Km_O2
write(*,'(a,2f8.4)') 'a_G, CexG: ',a_G,CexG
write(*,'(a,2f8.4)') 'a_O, CexO: ',a_O,CexO
cp%Cex(2) = CexG
cp%Cex(1) = CexO
cp%metab%HIF1 = 0.0
cp%metab%PDK1 = 1.0
do it = 1,10
	dt = it*1
	call OGLSolver(kcell,tstart,dt,ok)	
	do ichemo = OXYGEN,LACTATE
		Kin = chemo(ichemo)%membrane_diff_in
		Kout = chemo(ichemo)%membrane_diff_out
		cp%dMdt(ichemo) = area_factor*(Kin*cp%Cex(ichemo) - Kout*cp%Cin(ichemo))
	enddo
	! Reverse H transform:
	G = cp%Cin(2)
	b_G = G/((1 + K_Hb(ityp)*H)*(G**Hill_N_G + Hill_Km_G**Hill_N_G))
	CinG = (b_G/(1-b_G))*Hill_Km_G
	b_O = cp%Cin(1)/(fPDK*(cp%Cin(1) + Hill_Km_O2))
	CinO = (b_O/(1-b_O))*Hill_Km_O2
	write(*,'(a,f6.1,2f8.5,3e11.3)') 'Cin,CinG,dMdt: ',dt,cp%Cin(2),CinG,cp%dMdt(2)
	write(*,'(a,f6.1,2f8.5,3e11.3)') 'Cin,CinO,dMdt: ',dt,cp%Cin(1),CinO,cp%dMdt(1)
enddo
end subroutine

!----------------------------------------------------------------------------------
!----------------------------------------------------------------------------------
!subroutine CreateTables
!type(cell_type) :: testcell
!type(cell_type), pointer :: cp
!
!if (allocated(IC_O_lookup)) deallocate(IC_O_lookup)
!if (allocated(IC_G_lookup)) deallocate(IC_G_lookup)
!if (allocated(IC_L_lookup)) deallocate(IC_L_lookup)
!if (allocated(dMdt_O_lookup)) deallocate(dMdt_O_lookup)
!if (allocated(dMdt_G_lookup)) deallocate(dMdt_G_lookup)
!if (allocated(dMdt_L_lookup)) deallocate(dMdt_L_lookup)
!allocate(IC_O_lookup(N_EC_O,N_EC_G,N_EC_L))
!allocate(IC_G_lookup(N_EC_O,N_EC_G,N_EC_L))
!allocate(IC_L_lookup(N_EC_O,N_EC_G,N_EC_L))
!allocate(dMdt_O_lookup(N_EC_O,N_EC_G,N_EC_L))
!
!cp => testcell
!O_max = 0.18
!G_max = 5.5
!L_max = 50
!d_O = O_max/(N_EC_O - 1)
!d_G = G_max/(N_EC_G - 1)
!d_L = L_max/(N_EC_L - 1)
!do i = 1,N_EC_O
!	EC_O_array(i) = (i - 1)*d_O
!enddo
!do i = 1,N_EC_G
!	EC_G_array(i) = (i - 1)*d_G
!enddo
!do i = 1,N_EC_L
!	EC_L_array(i) = (i - 1)*d_L
!enddo
!
!do i_O = 1,N__EC_O
!do i_G = 1,N__EC_G
!do i_L = 1,N__EC_L
!	C_O = EC_O_array(i_O)
!	C_G = EC_G_array(i_G)
!	C_L = EC_L_array(i_L)
!	
!	
!	call optimiser(ityp, r_G, C_L, r_Pm, mp, x, y, C_P)
!	f_G_lookup(i_r_G,i_C_L,i_r_Pm) = x
!	f_P_lookup(i_r_G,i_C_L,i_r_Pm) = y
!	C_P_lookup(i_r_G,i_C_L,i_r_Pm) = C_P
!!	write(*,'(a,4f8.4,e12.3)') 'x, y, C_L, C_P: ',x,y,C_L,C_P
!!	write(*,'(a,4e12.3)') 'A, I, O, L: ',mp%A_rate, mp%I_rate, mp%O_rate, mp%L_rate
!	r_P = 2*(1-x)*r_G - V*(K1*C_P - K2*C_L)
!	r_A = 2*(1-x)*r_G + f_PA*(1-y)*r_P
!	r_I = x*r_G + y*r_P
!	r_O = f_PO*r_P
!	r_L = V*(K1*C_P - K2*C_L)
!!	write(*,'(a,4e12.3)') 'r_A, r_I: ',r_A,r_I,r_O,r_L
!enddo
!enddo
!enddo
!
!end subroutine

end module
