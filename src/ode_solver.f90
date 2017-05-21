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

!real(REAL_KIND) :: dt_rkc
contains

!----------------------------------------------------------------------------------
! Try to estimate SS Cin from Cex, assuming that other variables change slowly,
! i.e. treating C_G, C_L as approx constant for determining SS C_O2, then
! use RKC to update these concs.
!----------------------------------------------------------------------------------
function get_CO2_SS(kcell) result(C)
integer :: kcell
real(REAL_KIND) :: C
!real(REAL_KIND) :: tstart, dt
!logical :: ok
!integer :: ichemo, k, neqn, i
!real(REAL_KIND) :: t, tend
!real(REAL_KIND) :: timer1, timer2
! Variables for RKC
!integer :: info(4), idid
!real(REAL_KIND) :: rtol, atol(1)
!type(rkc_comm) :: comm_rkc(1)
!real(REAL_KIND) :: work_rkc(8+5*3*2)
integer :: ict
real(REAL_KIND) :: Cin(3), area_factor = 1.2**(2./3.)
real(REAL_KIND) :: C0, C1, Cex, dC, r0, r1, a, b, Kmem	! these are OXYGEN
type(cell_type), pointer :: cp
type(metabolism_type), pointer :: mp

!write(*,'(a,2i6,6f8.4)') 'OGLSolver: Cin,Cex: ',istep,kcell,cp%Cin(1:3),cp%Cex(1:3)
cp => cell_list(kcell)
mp => cp%metab
ict = cp%celltype
kcell_now = kcell
Cin(1:3) = cp%Cin(1:3)
Cex = cp%Cex(1)
dC = min(Cex/10, 0.0001)
call get_metab_rates(ict,mp,Cin)
r0 = mp%O_rate
C1 = Cin(1) - dC
Cin(1) = C1
call get_metab_rates(ict,mp,Cin)
r1 = mp%O_rate
a = r0
b = (r1 - r0)/dC
! linearising gives O2 consumption rate r(C) = a + b(Cex - C)
Kmem = area_factor*chemo(OXYGEN)%membrane_diff_in
! equating membrane flux = Kmem(Cex - C) to linearised r(C) = a + b(Cex - C)
C = Cex - a/(Kmem - b)
end function

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
type(metabolism_type), pointer :: mp
! Variables for RKC
integer :: info(4), idid
real(REAL_KIND) :: rtol, atol(1)
type(rkc_comm) :: comm_rkc(1)
real(REAL_KIND) :: work_rkc(8+5*3*2)
logical :: use_SS_CO2 = .false.

!write(*,'(a,2i6,6f8.4)') 'OGLSolver: Cin,Cex: ',istep,kcell,cp%Cin(1:3),cp%Cex(1:3)
cp => cell_list(kcell)
kcell_now = kcell

mp => cp%metab
call Set_f_GP(cp%celltype,mp)	! dropped argument C from monolayer_m version

if (use_SS_CO2) then
	neqn = 2
	cp%Cin(1) = get_CO2_SS(kcell)
	Cin(1) = cp%Cin(2)		! glucose and lactate become 1 and 2 for RKC
	Cin(2) = cp%Cin(3)
!	write(*,'(a,2e12.3)') 'OGLsolver: cp%Cex(1), cp%Cin(1): ',cp%Cex(1),cp%Cin(1)
else
	neqn = 3
	Cin(1:3) = cp%Cin(1:3)
endif

info(1) = 1
info(2) = 1		! = 1 => use spcrad() to estimate spectral radius, != 1 => let rkc do it
info(3) = 1
info(4) = 0
rtol = 1d-3		! this is probably too large
atol = rtol

idid = 0
t = tstart
tend = t + dt
call rkc(comm_rkc(1),neqn,f_rkc_OGL,Cin,t,tend,rtol,atol,info,work_rkc,idid,kcell)	! ! Using icase to pass kcell
if (idid /= 1) then
	write(logmsg,*) 'Solver: Failed at t = ',t,' with idid = ',idid
	call logger(logmsg)
	ok = .false.
	return
endif
!write(*,'(a,2e12.3)') 'After rkc: Cin: ',Cin(1:2)
if (use_SS_CO2) then
	cp%Cin(2:3) = Cin(1:2)	! with use_SS_CO2 it should never happen that Cin(1) > Cex(1)
else
	if (Cin(1) > cp%Cex(1)) then
		Cin(1) = cp%Cex(1)
	!	stop
	endif
	cp%Cin(1:3) = Cin
endif

! This next step is to ensure that cp has the current rates.  May not be needed.
ict = cp%celltype
mp => cp%metab
call get_metab_rates(ict,mp,cp%Cin(1:3))

end subroutine

!----------------------------------------------------------------------------------
! C_O2, C_G are the intracellular concentrations
! When C_O2 < CO_H both f_G and f_P are reduced, to 0 when C_O2 <= CO_L
! When C_G < CG_H f_G is reduced, to 0 when C_G <= CG_L
! When both concentrations are below the H threshold the reduction factors are 
! multiplied for f_G
!----------------------------------------------------------------------------------
subroutine Set_f_GP(ityp,mp)
integer :: ityp
type(metabolism_type), pointer :: mp
!real(REAL_KIND) :: C(:)
real(REAL_KIND) :: C_O2, C_G, Ofactor, Gfactor, ATPfactor
real(REAL_KIND) :: CO_L, CG_L, f_ATP_H, f_ATP_L, f
real(REAL_KIND) :: alfa = 0.2
logical :: use_ATP = .true.

if (use_ATP) then
	f_ATP_L = f_ATPg(ityp)
	f_ATP_H = f_ATPramp(ityp)*f_ATPg(ityp)
	f = mp%A_rate/r_A_norm
	if (f > f_ATP_H) then
		ATPfactor = 1
	elseif (f < f_ATP_L) then
		ATPfactor = 0
	else
		ATPfactor = (f - f_ATP_L)/(f_ATP_H - f_ATP_L)
		ATPfactor = smoothstep(ATPfactor)
	endif
!	write(*,'(a,2e12.3)') 'A_rate, ATPfactor: ',mp%A_rate,ATPfactor
	mp%f_G = alfa*ATPfactor*f_G_norm + (1-alfa)*mp%f_G
	mp%f_P = alfa*ATPfactor*f_P_norm + (1-alfa)*mp%f_P
else
!	C_O2 = C(1)
!	C_G = C(N1D+2)
	CO_L = 0.8*CO_H(ityp)
	CG_L = 0.8*CG_H(ityp)
	Ofactor = 1
	Gfactor = 1

	if (C_O2 < CO_L) then
		Ofactor = 0
	elseif (C_O2 < CO_H(ityp)) then
		Ofactor = (C_O2 - CO_L)/(CO_H(ityp) - CO_L)
	endif
	if (C_G < CG_L) then
		Gfactor = 0
	elseif (C_G < CG_H(ityp)) then
		Gfactor = (C_G - CG_L)/(CG_H(ityp) - CG_L)
	endif
	mp%f_G = Gfactor*Ofactor*f_G_norm
	mp%f_P = Ofactor*f_P_norm
!	write(*,'(a,6f8.4)') 'G, Ofactor: ',Gfactor, Ofactor, C_G, C_O2, mp%f_G, mp%f_P
endif
end subroutine

!----------------------------------------------------------------------------------
!----------------------------------------------------------------------------------
function smoothstep(x) result(f)
real(REAL_KIND) :: x, f
f = x*x*(3 - 2*x)
end function

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
real(REAL_KIND) :: average_volume = 1.2
logical :: use_average_volume = .true.

!write(*,'(a,3f9.6)') 'f_rkc_OGL: y: ',y(1:3)

!ict = icase
cp => cell_list(icase)
mp => cp%metab
ict = cp%celltype
area_factor = (average_volume)**(2./3.)
vol_cm3 = Vcell_cm3		!!!!!!!!!!!!!! temporary !!!!!!!!!!!!
if (neqn == 2) then
	Cin(1) = cp%Cin(1)
	Cin(2) = y(1)
	Cin(3) = y(2)
else
	Cin = y
endif
call get_metab_rates(ict,mp,Cin)
!write(*,'(a,3e12.3)') 'rates: ',mp%O_rate,mp%G_rate,mp%L_rate
do i = 1,neqn
	! First process IC reactions
	if (neqn == 2) then
		ichemo = i+1
	else
		ichemo = i
	endif
	C = y(i)
    membrane_kin = chemo(ichemo)%membrane_diff_in
    membrane_kout = chemo(ichemo)%membrane_diff_out
!	membrane_flux = area_factor*(membrane_kin*Cex_rkc(ichemo) - membrane_kout*C)
	membrane_flux = area_factor*(membrane_kin*cp%Cex(ichemo) - membrane_kout*C)
    if (ichemo == OXYGEN) then
		rate = mp%O_rate
		if (rate < 0) then
			write(*,*) 'O_rate < 0: ',rate
			stop
		endif
    elseif (ichemo == GLUCOSE) then	! 
		rate = mp%G_rate
    elseif (ichemo == LACTATE) then
		rate = -mp%L_rate
    endif
	dCreact = (membrane_flux - rate)/vol_cm3
!	if (ichemo == OXYGEN) then
!		dCreact = min(dCreact, (cp%Cex(ichemo) - C)/dt_rkc)
!	endif
	dydt(i) = dCreact
!	write(*,'(a,i4,5e12.3)') 'dydt: ',ichemo,dCreact,membrane_flux,rate,C,cp%Cex(ichemo)
	if (isnan(dydt(i))) then
		write(nflog,'(a,i2,4e12.3)') 'f_rkc_OGL: dydt isnan: ',ichemo,dydt(i),C,cp%Cex(ichemo),rate
		write(*,'(a,i2,4e12.3)') 'f_rkc_OGL: dydt isnan: ',ichemo,dydt(i),C,cp%Cex(ichemo),rate
		stop
	endif	
enddo
end subroutine

!----------------------------------------------------------------------------------
!----------------------------------------------------------------------------------
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
