! To test simple cell metabolism model
! Concentration of ATP varies by < 10%  https://en.wikipedia.org/wiki/Glycolysis#Intermediates_for_other_pathways

! Units:
!     time				s = seconds
!     distance			cm
!     volume			cm^3
!     mass				micromole = 10^-6 mol = mumol
!     flux				mumol/s
!     concentration		mumol/cm^3 = mM
!
! Need to comment out 'use chemokine' when used in the test program metab.exe
!
! Question: How do the results of this model translate into cell rate of volume growth?
!----------------------------------------------------------------------------------------------------------------
module metabolism
use real_kind_mod
use global
use chemokine
implicit none

! From spheroid-abm, the max rates of consumption of oxygen and glucose are:
!   oxygen:  6.25e-17 moles/cell/s
!   glucose: 6.80e-17 moles/cell/s
! We work with mumol/cell/sec, and convert these values by scaling by 1.0e6, to give
!   oxygen:  6.25e-11 mumol/cell/s
!   glucose: 6.80e-11 mumol/cell/s

real(REAL_KIND) :: Hill_Km_O2
real(REAL_KIND) :: Hill_N_O2
real(REAL_KIND) :: Hill_Km_G					! Hill Km for dependence of glycolysis rate on glucose
real(REAL_KIND) :: Hill_N_G						! Hill N for dependence of glycolysis rate on glucose 
real(REAL_KIND) :: Hill_Km_P(MAX_CELLTYPES)     ! Hill Km for dependence of pyruvate oxidation rate on pyruvate
real(REAL_KIND) :: Hill_N_P						! Hill N for dependence of pyruvate oxidation rate on pyruvate
real(REAL_KIND) :: K_H1(MAX_CELLTYPES)     ! HIF-1 k1
real(REAL_KIND) :: K_H2(MAX_CELLTYPES)     ! HIF-1 k2
real(REAL_KIND) :: K_Hb(MAX_CELLTYPES)     ! HIF-1 kb
real(REAL_KIND) :: K_PDK(MAX_CELLTYPES)    ! K_PDK
real(REAL_KIND) :: PDKmin(MAX_CELLTYPES)   ! PDKmin
real(REAL_KIND) :: C_O2_norm(MAX_CELLTYPES)
real(REAL_KIND) :: C_G_norm(MAX_CELLTYPES)
real(REAL_KIND) :: C_L_norm(MAX_CELLTYPES)
real(REAL_KIND) :: K_PL(MAX_CELLTYPES)     ! P -> L
real(REAL_KIND) :: K_LP(MAX_CELLTYPES)     ! L -> P
real(REAL_KIND) :: Apoptosis_rate(MAX_CELLTYPES)     

real(REAL_KIND) :: f_G_norm, f_P_norm, r_P_norm, r_G_norm, r_A_norm, r_I_norm, C_P_norm
real(REAL_KIND) :: G_maxrate, O2_maxrate

real(REAL_KIND) :: A_rate_base(MAX_CELLTYPES)	! total rate of production of ATP under full nutrition

real(REAL_KIND), parameter :: base_O_rate = 0.0e-11 

real(REAL_KIND),allocatable :: f_G_lookup(:,:,:)
real(REAL_KIND),allocatable :: f_P_lookup(:,:,:)
real(REAL_KIND),allocatable :: C_P_lookup(:,:,:)

contains

!--------------------------------------------------------------------------
! This is the normalised rate
!--------------------------------------------------------------------------
function get_glycosis_rate(ityp, H, G) result(rate)
integer :: ityp
real(REAL_KIND) :: H, G, rate

rate = G_maxrate* (1 + K_Hb(ityp)*H)* G**Hill_N_G / (G**Hill_N_G + Hill_Km_G**Hill_N_G)
end function

!!--------------------------------------------------------------------------
!!--------------------------------------------------------------------------
!function get_Poxidation_rate(ityp, PDK, P) result(rate)
!integer :: ityp
!real(REAL_KIND) :: PDK, P, rate
!
!rate = PDK*P**Hill_N_P/(P**Hill_N_P + Hill_Km_P**Hill_N_P)
!end function


!--------------------------------------------------------------------------
! Currently this sets up parameters for type 1 cells only.
!--------------------------------------------------------------------------
subroutine SetupMetabolism
integer :: ityp, it, N_P, N_O2
real(REAL_KIND) :: f_Gn, f_Pn, f_PO, f_PA, MM_O2, MM_P, V, K1, K2, Km_P, Km_O2, C_P, C_L, r_Gn, r_Pn, r_Ln, r_An, r_In
type(metabolism_type), pointer :: mp

Hill_Km_O2 = chemo(OXYGEN)%MM_C0
Hill_N_O2 = chemo(OXYGEN)%Hill_N
Hill_Km_G = chemo(GLUCOSE)%MM_C0
Hill_N_G = chemo(GLUCOSE)%Hill_N
Hill_N_P = 1
N_PO = 3
ityp = 1
V = Vcell_cm3
!C_O2_norm = 0.02
!C_G_norm = 0.5
!C_L_norm = 3.0

metabolic%HIF1 = 0
metabolic%PDK1 = 1

Km_O2 = chemo(OXYGEN)%MM_C0
N_O2 = chemo(OXYGEN)%Hill_N
O2_maxrate = chemo(OXYGEN)%max_cell_rate
G_maxrate = chemo(GLUCOSE)%max_cell_rate

do ityp = 1,1
	MM_O2 = f_MM(C_O2_norm(ityp),Km_O2,N_O2)
    mp => metabolic(ityp)
	K1 = K_PL(ityp)
	K2 = K_LP(ityp)
	N_P = 1
	Km_P = Hill_Km_P(ityp)
	f_Gn = N_GI(ityp)
	f_Pn = N_PI(ityp)
	f_PO = N_PO(ityp)
	f_PA = N_PA(ityp)
	C_L = C_L_norm(ityp)
	C_P = C_L

	do it = 1,10
		MM_P = f_MM(C_P,Km_P,N_P)
		r_Gn = get_glycosis_rate(ityp,0.0d0,C_G_norm(ityp))
		r_Ln = 2*(1-f_Gn)*r_Gn - MM_P*MM_O2*(O2_maxrate-base_O_rate)/(f_PO*(1-f_Pn))
		r_Pn = 2*(1-f_Gn)*r_Gn - r_Ln
		r_An = 2*(1-f_Gn)*r_Gn + f_PA*(1-f_Pn)*r_Pn
		r_In = f_Gn*r_Gn + f_Pn*r_Pn
!		write(*,'(a,5e11.3)') 'G,P,L,A,I rates: ',r_Gn,r_Pn,r_Ln,r_An,r_In
		C_P = (r_Ln/V + K2*C_L)/K1
!		write(*,'(a,3f8.4)') 'Km_P,MM_P,C_P: ',Km_P,MM_P,C_P
		if (C_P <= 0) then
			write(*,*) 'Error: C_P < 0: decrease N_PI and/or N_GI'
			stop
		endif
	enddo

	f_G_norm = f_Gn
	f_P_norm = f_Pn
	r_G_norm = r_Gn
	r_P_norm = r_Pn
	r_A_norm = r_An
	r_I_norm = r_In
	C_P_norm = C_P
!	ATPg(ityp) = f_ATPg(ityp)*r_A_norm
	ATPs(ityp) = f_ATPs(ityp)*r_A_norm
	mp%I_rate_max = r_I_norm
enddo

end subroutine

!!--------------------------------------------------------------------------
!!--------------------------------------------------------------------------
!function get_I2Divide(cp) result(I2div)
!type(cell_type), pointer :: cp
!real(REAL_KIND) :: I2div
!integer :: ityp
!
!ityp = cp%celltype
!I2div = cp%divide_time*metabolic(ityp)%I_rate_max
!end function

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
function get_HIF1steadystate(ityp,C_O) result(H)
integer :: ityp
real(REAL_KIND) :: C_O, H

H = exp(-K_H1(ityp)*C_O)
end function

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
subroutine analyticSetHIF1(ityp, C_O, H, dt)
integer :: ityp
real(REAL_KIND) :: C_O, H, dt
real(REAL_KIND) :: a, b, c, e, H0

e = K_H1(ityp)*C_O
if (e > 100) then
    H = 0
    return
endif
a = K_H2(ityp)
b = exp(e)
H0 = H
c = 1 - b*H0
H = (1 - c*exp(-a*b*dt))/b
end subroutine

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
subroutine analyticSetPDK1(ityp, H, P, dt)
integer :: ityp
real(REAL_KIND) :: P, H, dt
real(REAL_KIND) :: a, b, c, d, P0

a = K_PDK(ityp)
d = 1 - PDKmin(ityp)
b = (1 - d*H)
P0 = P
c = P0 - b
P = b + c*exp(-a*dt)
end subroutine

!--------------------------------------------------------------------------
! The fraction of pyruvate that is converted to acetyl-CoA depends on the
! rate of glycolysis.  The normalised rate has a maximum value of 1 under
! normoxic conditions. This corresponds to the minimum pyruvate oxidation fraction
! (equivalently, the maximum lactate fraction).  The minimum pyruvate oxidation
! fraction is the fraction of pyruvate that is directed to the Krebs cycle
! when both glucose and oxygen are in ample supply.  
! In fact this is the nominal minimum, and the fraction can be even less 
! if hypoxia elevates glycosis rate.
! The basic idea is that over a range of glycolysis rate, the total rate of
! production of ATP is a constant.
! The total is the sum of ATP produced by glycolysis and by pyruvate oxidation.
! If the glycolysis rate is reduced, the intermediate production rate can be maintained 
! by increasing the pyruvate oxidation fraction (reducing lactate production), to the limit
! of fraction = 1.  Further reductions in glycolysis will then reduce ATP production.
! REVISED
! Glycolysis:
! A fraction N_GI goes to make intermediates, I_rate = N_GI*G_rate
! the remainder makes pyruvate, PP_rate = 2*(1 - N_GI)*G_rate
! and A_rate = PP_rate
! Pyruvate oxidation:
! A fraction N_PI goes to make intermediates via the TCA cycle, I_rate = N_PI*P_rate
! the remainder (1 - N_PI) is fully oxidised, producing N_PA ATP/mole, (N_PA = 18)
! and A_rate from pyruvate oxidation = N_PA*(1 - N_PI)*P_rate
! Bill:
! The intermediates used for anabolism are downstream of acetyl-CoA (including acetyl-CoA itself).
! Yes, PDK1 reduces pyruvate utilisation (r) but this is upstream of acetyl-CoA. So the formalism 
! needs to have the intermediates coming from acetyl-CoA rather than pyruvate itself.
!--------------------------------------------------------------------------
subroutine get_metab_rates(ityp, mp, Cin)
integer :: ityp
type(metabolism_type), pointer :: mp
real(REAL_KIND) :: Cin(:)

!if (Cin(GLUCOSE) == 0) then
!	write(*,*) 'Glucose concentration = 0'
!	stop
!endif

call f_metab(ityp,mp,Cin(OXYGEN),Cin(GLUCOSE),Cin(LACTATE))
end subroutine

!--------------------------------------------------------------------------
! This is used in the monolayer model, where all cells have the same IC
! concentrations, and therefore the same metabolism.
! We may not need G_rate, PP_rate, P_rate
! We update only by cell type here, because in the monolayer all cells of the
! same type have the same metabolic state.
! NOT USED YET
!--------------------------------------------------------------------------
!subroutine update_metabolism
!type(metabolism_type), pointer :: mp
!integer :: kcell, ityp
!real(REAL_KIND) :: Itotal, I2Divide
!type(cell_type), pointer :: cp
!
!do ityp = 1,Ncelltypes
!	mp => metabolic(ityp)
!!	call get_metab_rates(ityp, mp, Caverage(1:MAX_CHEMO))
!enddo
!
!do kcell = 1,nlist
!    if (colony_simulation) then
!        cp => ccell_list(kcell)
!    else
!        cp => cell_list(kcell)
!    endif
!    if (cp%state == DEAD) cycle
!    Itotal = cp%metab%Itotal	! The only fields we don't want to change are Itotal and I2Divide, which are cell-specific
!    I2Divide = cp%metab%I2Divide
!    cp%metab = metabolic(cp%celltype)
!    cp%metab%Itotal = Itotal
!    cp%metab%I2Divide = I2Divide
!enddo
!
!end subroutine


!--------------------------------------------------------------------------
! Use:
!	r_G for dG/dt, the rate of glycolysis = rate of glucose consumption
!	f_G for N_GI, the fraction of dG/dt that goes to intermediates
!	r_P for dP/dt, rate of utilisation of pyruvate
!	f_P for N_PI, the fraction of dP/dt that goes to intermediates
!	r_A for dA/dt, rate of ATP production
!	r_I for dI/dt, rate of intermediates production
!
!	r_G_norm for dG/dt under normal conditions, no nutrient constraints, H = 1
!	f_G_norm for the f_G under normal conditions (upper bound of f_G)
!	r_P_norm for dP/dt under normal conditions
!	f_P_norm for f_P under normal conditions (upper bound of f_P)
!	r_A_norm for r_A under normal conditions
!	r_I_norm for r_I under normal conditions
!	alpha = r_P as a fraction of dP/dt under normal conditions
!
!	C_G for IC glucose concentration
!	C_O2 for IC oxygen concentration
!	C_P for IC pyruvate concentration
!	C_L for IC lactate concentration
!	
! r_P = 2*(1 - f_G)*r_G + V*(K2*C_L - K1*C_P - dC_P/dt) = fPDK*r_P_max*MM(O2)*MM(C_P)
! with the constraint that C_P >= 0
! Steady-state approach may not be feasible, because if there is 
! plenty of glucose but O2 is very low, r_G will be high but r_P will tend
! towards 0.  This must lead to an increase in C_P.
!--------------------------------------------------------------------------
subroutine f_metab(ityp, mp, C_O2_, C_G_, C_L_)
integer :: ityp
type(metabolism_type), pointer :: mp
real(REAL_KIND) :: C_O2_, C_G_, C_L_
real(REAL_KIND) :: C_O2, C_G, C_L
real(REAL_KIND) :: r_G, fPDK
real(REAL_KIND) :: f_G, f_P, r_P, r_A, r_I, r_L, f_PO, f_PA
real(REAL_KIND) :: K1, K2, C_P
real(REAL_KIND) :: r_GP, r_GA, r_PA, r_Pm, V, Km_O2, Km_P, c1, c2, a, b, c, d, e, MM_P, MM_O2
real(REAL_KIND) :: r_Pm_base, r_A_target, f_G_new, C_P_max, dA
integer :: N_O2, N_P, it
real(REAL_KIND) :: r_G_threshold = 1.0e-12
real(REAL_KIND) :: dA_threshold = 1.0e-16

C_O2 = max(0.0,C_O2_)
C_G = max(0.0,C_G_)
C_L = max(0.0,C_L_)

N_O2 = Hill_N_O2
Km_O2 = Hill_Km_O2
N_P = 1
Km_P = Hill_Km_P(ityp)
V = Vcell_cm3		! should be actual cell volume cp%V
f_PO = N_PO(ityp)
f_PA = N_PA(ityp)
K1 = K_PL(ityp)
K2 = K_LP(ityp)
f_G = f_G_norm
f_P = f_P_norm

mp%G_rate = get_glycosis_rate(ityp,mp%HIF1,C_G)
r_G = mp%G_rate
!if (r_G == 0) then
!	write(*,'(a,3e12.3)') 'r_G = 0: G_maxrate,mp%HIF1,C_G: ',G_maxrate,mp%HIF1,C_G
!	stop
!endif
fPDK = mp%PDK1
MM_O2 = f_MM(C_O2,Km_O2,N_O2)
r_Pm_base = fPDK*MM_O2*(O2_maxrate-base_O_rate)/f_PO	! note that MM_P is not here, since it varies it is added as needed

if (dbug) write(*,'(a,5e12.3)') 'r_G,C_L,r_Pm: ',r_G,C_L,r_Pm_base,fPDK,MM_O2
if (r_G < r_G_threshold) then
	r_G = max(0.0d0,r_G)
	f_G = 0
	f_P = f_P_norm*r_G/r_G_threshold
	r_GA = 2*r_G
!	r_P = r_GA - V*(K1*C_P - K2*C_L) = MM(C_P)*r_Pm_base/(1-f_P) = r_Pm_base*C_P/((1-f_P)*(Km_P + C_P)) 
	e = r_GA + V*K2*C_L
	a = V*K1
	b = r_Pm_base/(1-f_P) - e + V*K1*Km_P
	c = -e*Km_P
	d = sqrt(b*b - 4*a*c)
	C_P = (-b + d)/(2*a)
	r_P = r_GA - V*(K1*C_P - K2*C_L)	
	r_PA = f_PA*(1-f_P)*r_P
else

	r_A_target = r_A_norm
	it = 0
	do
		it = it + 1
		if (it > 30) then
			write(*,*) 'f_metab: it,kcell: ',it,kcell_now
			write(*,'(a,3f9.6)') 'C O2,G,L: ',C_O2, C_G, C_L
			write(*,'(a,2f9.6)') 'HIF1,fPDK: ',mp%HIF1,fPDK
			stop
		endif
		if (dbug) write(*,'(a,i2,e12.3)') 'r_A_target: ',it,r_A_target
		c1 = r_Pm_base*(f_PA + 1/(1-f_P))
		c2 = r_A_target + V*K2*C_L
		a = V*K1
		b = V*K1*Km_P + c1 - c2
		c = -Km_P*c2
	!	C_P = solve_C_P(a,b,c)	! solves a*x^2 + b*x + c = 0
		d = sqrt(b*b - 4*a*c)
		C_P = (-b + d)/(2*a)
		r_GA = r_A_target - f_PA*r_Pm_base*C_P/(Km_P + C_P)
		r_P = r_GA - V*(K1*C_P - K2*C_L)
		r_P = max(0.0,r_P)
		C_P_max = (r_GA + V*K2*C_L)/(V*K1)
	!	if (C_P > C_P_max) then
	!		write(*,*) 'C_P_max is exceeded!!!!'
	!	endif
		if (dbug) write(*,'(a,4f8.4)') 'f_P, f_G, C_P, C_P_max: ',f_P, f_G, C_P, C_P_max
		if (dbug) write(*,'(a,2e12.3)') 'r_P, r_GA: ',r_P,r_GA
		dA = r_GA - 2*(1-f_G)*r_G
		if (dA < dA_threshold) then
			r_PA = f_PA*(1-f_P)*r_P
			r_A = r_GA + r_PA
			if (dbug) write(*,'(a,3e12.3)') 'r_GA,r_PA,r_A: ',r_GA,r_PA,r_A/r_A_norm
			if (dbug) write(*,*)
			exit
		else
			if (dbug) write(*,'(a,e12.3)') 'dA: ',dA
			f_G_new = max(0.0d0,f_G - dA/(2*r_G))
			dA = dA - 2*(f_G - f_G_new)*r_G
			if (f_G_new == 0) dA = 1.01*dA
			if (dbug) write(*,'(a,2e12.3)') 'f_G_new, dA: ',f_G_new, dA
			if (dbug) write(*,*)
			f_G = f_G_new
			r_A_target = r_A_target - dA
			cycle
		endif
	enddo
!if (kcell_now == 8000 .and. C_O2 == 0) then
!	write(*,'(a,3i6,2f8.4)') 'kcell,istep,it,C_O2,C_G: ',kcell_now,istep,it,C_O2,C_G
!endif
	! I think this has no function when ATPg = ATPs
!	if (r_A < ATPg(ityp)) then
!		f_G = 0
!		if (ATPg(ityp) > ATPs(ityp)) then
!			f_P = f_P_norm*(r_A - ATPs(ityp))/(ATPg(ityp) - ATPs(ityp))
!		else
!			r_P = 0
!		endif
!		r_GA = 2*r_G
!		e = r_GA + V*K2*C_L
!		a = V*K1
!		b = r_Pm_base/(1-f_P) - e + V*K1*Km_P
!		c = -e*Km_P
!		d = sqrt(b*b - 4*a*c)
!		C_P = (-b + d)/(2*a)
!		r_P = r_GA - V*(K1*C_P - K2*C_L)
!		r_P = max(0.0,r_P)
!		r_PA = f_PA*(1-f_P)*r_P
!	endif

endif

mp%A_rate = r_GA + r_PA				! production
mp%I_rate = f_G*r_G + f_P*r_P		! production
mp%P_rate = r_P						! utilisation
mp%O_rate = f_PO*r_P*(1-f_P)		! consumption
mp%L_rate = V*(K1*C_P - K2*C_L)		! production
mp%C_P = C_P
! Add base rate correction
mp%O_rate = mp%O_rate + MM_O2*base_O_rate	! TRY REMOVING THIS
mp%f_G = f_G
mp%f_P = f_P
if (dbug) then
	write(*,'(a,6e10.3)') 'r_G,A,P,I,L,O2: ',mp%G_rate,mp%A_rate,mp%P_rate,mp%I_rate,mp%L_rate,mp%O_rate
endif
end subroutine

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
function solve_C_P(a,b,c) result(x)
real(REAL_KIND) :: a, b, c, x
real(REAL_KIND) :: d

d = b*b - 4*a*c
if (d < 0) then
	write(*,*) 'Error: solve_C_P: a,b,c,d: ',a,b,c,d
!	write(*,'(a,3e12.3)') 'a,b,c: ',a,b,c
!	write(*,'(a,e12.3)') '-b/2a: ',-b/(2*a)
	x = 0
	return
else
	d = sqrt(d)
endif
x = (-b + d)/(2*a)
if (x < 0) then
	write(*,*) 'solve_C_P: x < 0: ',x
	stop
endif
end function

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
function f_MM(C,Km,N) result(v)
real(REAL_KIND) :: C, Km, v
integer :: N

v = C**N/(Km**N + C**N)
end function

end module

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!program main
!use metabolism
!integer :: nt = 10
!real(REAL_KIND) :: Hss, dGdt, dOdt
!
!open(nfout,file='metab.out',status='replace')
!call setup
!C_O = 0.01
!C_G = 3.0
!C_L = 1.0
!C_P = C_L*K_LP/K_PL
!M_I = 0
!
!Hss = get_HIF1steadystate(C_O)
!write(*,'(a,2e12.3)') 'steady-state HIF-1: ',C_O,Hss
!H = 1.0*Hss
!H = min(H, 1.0)
!H = max(H, 0.0)
!dGdt = -0.00004    ! test rate of change of ambient glucose
!dOdt = -0.0000004    ! test rate of change of ambient O2
!write(nfout,'(a)') '   hour C_G C_O C_L C_P H G_rate PP_rate P_rate L_rate A_rate I_rate O_rate'
!do istep = 1,60*24
!    C_G = C_G + dGdt*dt
!    C_G = max(C_G,0.0)
!    C_G = min(C_G,5.5)
!!    C_O = C_O + dOdt*dt
!!    C_O = max(C_O,0.0)
!!    C_O = min(C_O,0.10)
!    call timestep(nt)
!!    write(*,'(a,e12.3)') 'M_I: ', M_I
!enddo
!!write(*,*) 'Hss: ',Hss
!end program
