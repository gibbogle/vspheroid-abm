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
module FCN_mod
use real_kind_mod
implicit none

real(REAL_KIND) :: FCN_a1, FCN_b1, FCN_c1
real(REAL_KIND) :: FCN_r_Pm_base, FCN_fPn, FCN_Km_P, FCN_fGn, FCN_r_G, FCN_N_GA, FCN_V, FCN_K1, FCN_K2, FCN_C_L

end module

!----------------------------------------------------------------------------------------------------------------
module metabolism
use real_kind_mod
use global
#if .not. EXCEL
use chemokine
#endif
implicit none

integer, parameter :: FGP_SOLVER_MAXATP_TANDEM = 1
integer, parameter :: FGP_SOLVER_MAXATP_STAGED = 2
integer, parameter :: FGP_SOLVER_SURVIVAL_STAGED = 3

! From spheroid-abm, the max rates of consumption of oxygen and glucose are: 
!   oxygen:  6.25e-17 moles/cell/s
!   glucose: 6.80e-17 moles/cell/s
! We work with mumol/cell/sec, and convert these values by scaling by 1.0e6, to give
!   oxygen:  6.25e-11 mumol/cell/s
!   glucose: 6.80e-11 mumol/cell/s

real(REAL_KIND) :: Hill_Km_O2
real(REAL_KIND) :: Hill_N_O2
real(REAL_KIND) :: Hill_Km_G	! Hill Km for dependence of glycolysis rate on glucose
real(REAL_KIND) :: Hill_N_G		! Hill N for dependence of glycolysis rate on glucose 
real(REAL_KIND) :: Hill_Km_P    ! Hill Km for dependence of pyruvate oxidation rate on pyruvate
real(REAL_KIND) :: Hill_N_P		! Hill N for dependence of pyruvate oxidation rate on pyruvate
real(REAL_KIND) :: K_H1			! HIF-1 k1
real(REAL_KIND) :: K_H2			! HIF-1 k2
real(REAL_KIND) :: K_Hb			! HIF-1 kb
real(REAL_KIND) :: K_PDK		! K_PDK
real(REAL_KIND) :: PDKmin		! PDKmin
real(REAL_KIND) :: C_O2_norm
real(REAL_KIND) :: C_G_norm
real(REAL_KIND) :: C_L_norm
real(REAL_KIND) :: O2_baserate
real(REAL_KIND) :: G_baserate
real(REAL_KIND) :: K_PL			! P -> L
real(REAL_KIND) :: K_LP			! L -> P
real(REAL_KIND) :: r_P_norm, r_G_norm, r_A_norm, r_I_norm, C_P_norm, C_A_norm
real(REAL_KIND) :: G_maxrate, O2_maxrate
integer :: fgp_solver

!real(REAL_KIND) :: A_rate_base(MAX_CELLTYPES)	! total rate of production of ATP under full nutrition
!real(REAL_KIND), parameter :: base_O_rate = 0.0e-11 
!real(REAL_KIND),allocatable :: f_G_lookup(:,:,:)
!real(REAL_KIND),allocatable :: f_P_lookup(:,:,:)
!real(REAL_KIND),allocatable :: C_P_lookup(:,:,:)
!real(REAL_KIND) :: A_fract
!real(REAL_KIND) :: r_G_threshold = 1.0e-17		! possibly should be a GUI parameter, since it influences f_P
!real(REAL_KIND) :: dA_threshold = 1.0e-18

real(REAL_KIND), parameter :: r_H = 1.21e-4
real(REAL_KIND), parameter :: d_H = 2.3e-3
real(REAL_KIND), parameter :: Km_H = 20*0.18/160	! Kelly had 2 mmHg, Bill suggests 20
#if EXCEL
real(REAL_KIND), parameter :: PI = 4*atan(1.0)
#endif
logical :: use_Kelly = .false.

contains

!--------------------------------------------------------------------------
! This is NOT the normalised rate
! H = HIF1
! G = glucose concentration
! TRY:
!	a variable Km, with n = 1 
!	double (parallel) processes
!	boost rate at low C_G
!--------------------------------------------------------------------------
function get_glycosis_rate(H, G) result(rate)
integer :: ityp
real(REAL_KIND) :: H, G, rate
real(REAL_KIND) :: metab

metab = glucose_metab(G)
rate = G_maxrate* (1 + K_Hb*H)*metab
end function

!--------------------------------------------------------------------------
! Currently this sets up parameters for type 1 cells only.
!--------------------------------------------------------------------------
subroutine SetupMetabolism(mp,ok)
type(metabolism_type), pointer :: mp
logical :: ok
integer :: ityp, it, N_P, N_O2
real(REAL_KIND) :: f_Gn, f_Pn, f_PO, f_PA, MM_O2, MM_P, V, K1, K2, Km_P, C_P, C_L, r_Gn, r_Pn, r_Ln, r_An, r_In
real(REAL_KIND) :: Km_O2_factor = 1
real(REAL_KIND) :: average_volume = 1.2
real(REAL_KIND) :: H, fPDK, rP1, rP2, rP3

if (noSS) then 
    write(nflog,*) 'SetupMetabolism: not using SS solver'
else
    write(nflog,*) 'SetupMetabolism: using SS solver'
endif
ok = .true.
V = Vcell_cm3*average_volume
!mp => phase_metabolic(1)

Hill_Km_O2 = chemo(OXYGEN)%MM_C0
Hill_N_O2 = chemo(OXYGEN)%Hill_N
Hill_Km_G = chemo(GLUCOSE)%MM_C0
Hill_N_G = chemo(GLUCOSE)%Hill_N
Hill_N_P = 1
Hill_Km_O2 = Km_O2_factor*Hill_Km_O2
N_O2 = chemo(OXYGEN)%Hill_N
O2_maxrate = chemo(OXYGEN)%max_cell_rate
G_maxrate = chemo(GLUCOSE)%max_cell_rate

MM_O2 = f_MM(C_O2_norm,Hill_Km_O2,N_O2)
K1 = K_PL
K2 = K_LP
N_P = 1
Km_P = Hill_Km_P
f_Gn = f_G_norm  ! input parameters
f_Pn = f_P_norm
f_PO = N_PO
f_PA = N_PA
C_L = C_L_norm
C_P = (K2/K1)*C_L

H = get_HIF1steadystate(C_O2_norm)
call analyticSetPDK1(H, fPDK, 1.0d10)
mp%HIF1 = H
mp%PDK1 = fPDK

if (chemo(LACTATE)%used) then
	do it = 1,10	! new method, set f_G, f_P
		MM_P = f_MM(C_P,Km_P,N_P)
		r_Gn = get_glycosis_rate(H,C_G_norm)
		r_Ln = N_GA*(1-f_Gn)*r_Gn - fPDK*MM_P*MM_O2*O2_maxrate/(f_PO*(1-f_Pn))
		r_Pn = N_GA*(1-f_Gn)*r_Gn - r_Ln
		r_An = N_GA*(1-f_Gn)*r_Gn + f_PA*(1-f_Pn)*r_Pn
!		r_In = f_Gn*r_Gn + f_Pn*r_Pn
		r_In = f_Gn*r_Gn*N_GI + f_Pn*r_Pn*N_PI
		write(nflog,'(a,5e11.3)') 'G,P,L,A,I rates: ',r_Gn,r_Pn,r_Ln,r_An,r_In
		C_P = (r_Ln/V + K2*C_L)/K1
		if (C_P <= 0) then
			write(*,*) 'Error: C_P < 0: decrease N_PI and/or N_GI'
			stop
		endif
	enddo
else
	! We need r_O2 = O2_maxrate = f_PO*2*(1 - f_Gn)*(1 - f_Pn)*G_maxrate  
	write(nflog,'(a,2f8.4)') 'O2_maxrate/G_maxrate, f_PO*2*(1 - f_Gn)*(1 - f_Pn): ',O2_maxrate/G_maxrate, f_PO*2*(1 - f_Gn)*(1 - f_Pn)
	f_Pn = 1 - O2_maxrate/(G_maxrate*f_PO*2*(1 - f_Gn))
	write(nflog,'(a,f8.3,a,f8.3)') 'reset f_Pn to: ',f_Pn,'  f_PO: ',f_PO
	r_Gn = get_glycosis_rate(0.0d0,C_G_norm)
	r_Ln = 0
	r_Pn = N_GA*(1-f_Gn)*r_Gn
	r_An = N_GA*(1-f_Gn)*r_Gn + f_PA*(1-f_Pn)*r_Pn
!	r_In = f_Gn*r_Gn + f_Pn*r_Pn
	r_In = f_Gn*r_Gn*N_GI + f_Pn*r_Pn*N_PI
	C_P = 0
    f_G_norm = f_Gn
    f_P_norm = f_Pn
endif

r_G_norm = r_Gn
r_P_norm = r_Pn
r_A_norm = r_An
r_I_norm = r_In
C_P_norm = C_P
C_A_norm = C_P_norm		!#################  Just because we have no idea! #################
mp%f_G = f_Gn
mp%f_P = f_Pn
mp%C_P = C_P_norm
mp%C_A = C_A_norm
ATPg = f_ATPg*r_A_norm
ATPs = f_ATPs*r_A_norm
mp%I_rate_max = r_I_norm
mp%A_rate = r_A_norm
mp%I_rate = r_I_norm
mp%G_rate = r_G_norm
mp%O_rate = f_PO*(1 - f_Pn)*r_Pn
O2_baserate = 0.0*O2_maxrate
G_baserate = 0.0*G_maxrate
write(nflog,'(a,3e12.3)') 'f_G_norm,f_P_norm,ATPg: ',f_G_norm,f_P_norm,ATPg
end subroutine

!!--------------------------------------------------------------------------
!!--------------------------------------------------------------------------
!function get_I2Divide(cp) result(I2div)
!type(cell_type), pointer :: cp
!real(REAL_KIND) :: I2div
!integer :: ityp
!
!ityp = cp%celltype
!I2div = cp%divide_time*metabolic%I_rate_max
!end function

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
function get_HIF1steadystate1(C_O) result(H)
real(REAL_KIND) :: C_O, H
real(REAL_KIND) :: b

if (use_Kelly) then
	b = 1 + d_H*C_O/(r_H*(Km_H + C_O))
	H = 1/b
else
	H = exp(-K_H1*C_O)
endif
end function

!--------------------------------------------------------------------------
! Use K_H1 for the exponent, K_H2 for the rate coefficient
!--------------------------------------------------------------------------
function get_HIF1steadystate(C_O) result(H)
real(REAL_KIND) :: C_O, H
real(REAL_KIND) :: x
real(REAL_KIND) :: C_O_max = 0.18

x = min(C_O/C_O_max,1.0)
H = (1-x)**K_H1
!write(*,*) 'get_HIF1steadystate: K_H1, K_H2, K_HB: ', K_H1, K_H2, K_Hb, K_PDK, PDKmin
!write(*,*) 'C_O, x, H: ',C_O, x, H
end function

!--------------------------------------------------------------------------
!With the Kelly2008 model:
!  a = rH
!  b = 1 + dH*C_O/(rH*(Km + C_O))
!--------------------------------------------------------------------------
subroutine analyticSetHIF1z(C_O, H, dt)
real(REAL_KIND) :: C_O, H, dt
real(REAL_KIND) :: a, b, c, ee, H0

if (use_Kelly) then
	a = r_H
	b = 1 + d_H*C_O/(r_H*(Km_H + C_O))
else
	ee = K_H1*C_O
	if (ee > 100) then
		H = 0
		return
	endif
	a = K_H2
	b = exp(ee)
endif
H0 = H
c = 1 - b*H0
H = (1 - c*exp(-a*b*dt))/b
end subroutine

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
subroutine analyticSetHIF1(C_O, H, dt)
integer :: ityp
real(REAL_KIND) :: C_O, H, dt
real(REAL_KIND) :: x, H0, Heq
real(REAL_KIND) :: C_O_max = 0.18

x = min(C_O/C_O_max,1.0)
Heq = (1-x)**K_H1
H0 = H
H = Heq + (H0 - Heq)*exp(-K_H2*dt)
end subroutine

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
subroutine analyticSetPDK1(H, P, dt)
integer :: ityp
real(REAL_KIND) :: P, H, dt
real(REAL_KIND) :: a, b, c, d, P0

a = K_PDK
d = 1 - PDKmin
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
! the remainder makes pyruvate, PP_rate = N_GA*(1 - N_GI)*G_rate
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
subroutine get_metab_rates(mp, Cin)
integer :: ityp
type(metabolism_type), pointer :: mp
real(REAL_KIND) :: Cin(:)

if (noSS) then
    call f_metab_noSS(mp,Cin(OXYGEN),Cin(GLUCOSE),Cin(LACTATE),Cin(4))
else
    call f_metab(mp,Cin(OXYGEN),Cin(GLUCOSE),Cin(LACTATE))
endif
end subroutine

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
! r_P = N_GA*(1 - f_G)*r_G + V*(K2*C_L - K1*C_P - dC_P/dt) = fPDK*r_P_max*MM(O2)*MM(C_P)
! with the constraint that C_P >= 0
! Steady-state approach may not be feasible, because if there is 
! plenty of glucose but O2 is very low, r_G will be high but r_P will tend
! towards 0.  This must lead to an increase in C_P.
!--------------------------------------------------------------------------

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
subroutine f_metab(mp, C_O2_, C_G_, C_L_)
integer :: ityp
type(metabolism_type), pointer :: mp
real(REAL_KIND) :: C_O2_, C_G_, C_L_
real(REAL_KIND) :: C_O2, C_G, C_L
real(REAL_KIND) :: r_G, fPDK
real(REAL_KIND) :: f_G, f_P, r_P, r_A, r_I, r_L, f_PO, f_PA
real(REAL_KIND) :: K1, K2, C_P
real(REAL_KIND) :: r_GP, r_GA, r_PA, r_Pm, V, Km_O2, Km_P, q, a, b, c, d, e, MM_P, MM_O2, MM_G, Km_GO
real(REAL_KIND) :: r_GI, r_PI, r_O2
real(REAL_KIND) :: r_Pm_base
integer :: N_O2, N_P, it
real(REAL_KIND) :: average_volume = 1.2

C_O2 = max(0.0,C_O2_)
C_G = max(0.0,C_G_)
C_L = max(0.0,C_L_)

N_O2 = Hill_N_O2
Km_O2 = Hill_Km_O2
N_P = 1
Km_P = Hill_Km_P
V = Vcell_cm3*average_volume		! should be actual cell volume cp%V
f_PO = N_PO
f_PA = N_PA
K1 = K_PL
K2 = K_LP
f_G = mp%f_G
f_P = mp%f_P

mp%G_rate = get_glycosis_rate(mp%HIF1,C_G)
r_G = mp%G_rate
fPDK = mp%PDK1
MM_O2 = f_MM(C_O2,Km_O2,N_O2)
if (chemo(LACTATE)%used) then

	a = V*K1
	b = fPDK*O2_maxrate*MM_O2/(f_PO*(1 - f_P)) + Km_P*V*K1 - N_GA*(1 - f_G)*r_G - V*K2*C_L
	c = -Km_P*(N_GA*(1 - f_G)*r_G + V*K2*C_L)
	d = sqrt(b*b - 4*a*c) 
	C_P = (-b + d)/(2*a)
	
	r_O2 = fPDK*O2_maxrate*MM_O2*C_P/(Km_P + C_P)
	r_P = r_O2/(f_PO*(1 - f_P))
	r_L = N_GA*(1 - f_G)*r_G - r_P
	
	mp%A_rate = (1-f_G)*r_G*N_GA + (1-f_P)*r_P*N_PA	! production
!	mp%I_rate = f_G*r_G + f_P*r_P					! production
	mp%I_rate = f_G*r_G*N_GI + f_P*r_P*N_PI			! production
	mp%P_rate = r_P									! utilisation
	mp%O_rate = r_O2								! consumption
	mp%L_rate = r_L									! production
	mp%C_P = C_P
else    ! not corrected for revised treatment of N_GI, N_PI etc
	r_P = fPDK*MM_O2*N_GA*(1 - f_G)*r_G
	r_GI = r_G - r_P/2
	r_GA = r_P
	! Need to adjust f_P to maintain r_A at r_A_norm = N_GA*(1-f_Gn)*r_Gn + N_PA*(1-f_Pn)*r_Pn
	f_P = 1 - (r_A_norm - N_GA*(1-f_G)*r_G)/(N_PA*r_P)
	f_P = max(f_P,0.0)
	f_P = min(f_P,1.0)
	r_PI = f_P*r_P
	r_PA = f_PA*(1 - f_P)*r_P
	r_O2 = f_PO*(1 - f_P)*r_P
	mp%A_rate = r_GA + r_PA
	mp%I_rate = r_GI + r_PI
	mp%P_rate = r_P
	mp%O_rate = r_O2
	mp%L_rate = 0
	mp%C_P = 0
endif
	
! Add base rate correction
mp%O_rate = mp%O_rate + O2_baserate
mp%G_rate = mp%G_rate + G_baserate
!write(nflog,'(a,i8,2e12.3)') 'G_rate: ',istep,mp%HIF1,mp%G_rate
end subroutine

!--------------------------------------------------------------------------
! This method does not assume steady-state for C_P, instead the d.e. for 
! C_P is solved.
!--------------------------------------------------------------------------
subroutine f_metab_noSS(mp, C_O2_, C_G_, C_L_, C_P_)
integer :: ityp
type(metabolism_type), pointer :: mp
real(REAL_KIND) :: C_O2_, C_G_, C_L_, C_P_
real(REAL_KIND) :: C_O2, C_G, C_L, C_P
real(REAL_KIND) :: r_G, fPDK
real(REAL_KIND) :: f_G, f_P, r_P, r_A, r_I, r_L, f_PO, f_PA
real(REAL_KIND) :: K1, K2
real(REAL_KIND) :: r_GP, r_GA, r_PA, r_Pm, V, Km_O2, Km_P, q, a, b, c, d, e, MM_P, MM_O2, MM_G, Km_GO
real(REAL_KIND) :: r_GI, r_PI, r_O2
real(REAL_KIND) :: r_Pm_base
integer :: N_O2, N_P, it
real(REAL_KIND) :: average_volume = 1.2

C_O2 = max(0.0,C_O2_)
C_G = max(0.0,C_G_)
C_L = max(0.0,C_L_)
C_P = max(0.0,C_P_)

N_O2 = Hill_N_O2
Km_O2 = Hill_Km_O2
N_P = 1
Km_P = Hill_Km_P
V = Vcell_cm3*average_volume		! should be actual cell volume cp%V
f_PO = N_PO
f_PA = N_PA
K1 = K_PL
K2 = K_LP
f_G = mp%f_G
f_P = mp%f_P

mp%G_rate = get_glycosis_rate(mp%HIF1,C_G)
r_G = mp%G_rate
fPDK = mp%PDK1
MM_O2 = f_MM(C_O2,Km_O2,N_O2)

	r_O2 = fPDK*O2_maxrate*MM_O2*C_P/(Km_P + C_P)
	r_P = r_O2/(f_PO*(1 - f_P))
	r_L = N_GA*(1 - f_G)*r_G - r_P

	mp%A_rate = (1-f_G)*r_G*N_GA + (1-f_P)*r_P*N_PA	! production
!	mp%I_rate = f_G*r_G + f_P*r_P					! production
	mp%I_rate = f_G*r_G*N_GI + f_P*r_P*N_PI			! production
	mp%P_rate = r_P									! utilisation
	mp%O_rate = r_O2								! consumption
	mp%L_rate = r_L									! production
	mp%C_P = C_P
	
! Add base rate correction
mp%O_rate = mp%O_rate + O2_baserate
mp%G_rate = mp%G_rate + G_baserate
end subroutine

!--------------------------------------------------------------------------
! Only for no lactate
! NOT CORRECTED
!--------------------------------------------------------------------------
subroutine f_metab_noL(mp, C_O2_, C_G_, C_L_)
integer :: ityp
type(metabolism_type), pointer :: mp
real(REAL_KIND) :: C_O2_, C_G_, C_L_
real(REAL_KIND) :: C_O2, C_G, C_L
real(REAL_KIND) :: r_G, fPDK
real(REAL_KIND) :: f_G, f_P, r_P, r_A, r_I, r_L, f_PO, f_PA
real(REAL_KIND) :: K1, K2, C_P
real(REAL_KIND) :: r_GP, r_GA, r_PA, r_Pm, V, Km_O2, Km_P, q, a, b, c, d, e, MM_P, MM_O2
real(REAL_KIND) :: r_GI, r_PI, r_O2
real(REAL_KIND) :: r_Pm_base
integer :: N_O2, N_P, it

C_O2 = max(0.0,C_O2_)
C_G = max(0.0,C_G_)
C_L = max(0.0,C_L_)

N_O2 = Hill_N_O2
Km_O2 = Hill_Km_O2
N_P = 1
Km_P = Hill_Km_P
f_PO = N_PO
f_PA = N_PA
f_G = mp%f_G
f_P = mp%f_P

MM_O2 = f_MM(C_O2,Km_O2,N_O2)
mp%G_rate = MM_O2*get_glycosis_rate(mp%HIF1,C_G)
r_G = mp%G_rate
fPDK = mp%PDK1

r_G = r_G
r_P = fPDK*N_GA*(1 - f_G)*r_G
r_GI = f_G*r_G		! This implies glucose disappears when fPDK < 1
r_PI = f_P*r_P
r_GA = r_P
r_PA = f_PA*(1 - f_P)*r_P

! Now f_G and f_P are computed before the time step solve
if (.false.) then

if (r_GA + r_PA < ATPg) then	! adjust f_P to maintain ATPg
	! r_GA+r_PA = r_P*(f_P + f_PA*(1 - f_P)) = ATPg
	! => f_P = (ATPg/r_P - f_PA)/(1 - f_PA)
	f_P = (ATPg/r_P - f_PA)/(1 - f_PA)
	f_P = min(f_P,1.0)
	if (f_P < 0) then
		f_P = 0
		f_G = 1 - ATPg/(fPDK*N_GA*r_G*(1 + f_PA))
		f_G = max(f_G,0.0)
	endif
!	write(nflog,'(a,2f8.3)') 'f_P, f_G: ',f_P,f_G
!	stop
	r_P = fPDK*N_GA*(1 - f_G)*r_G
	r_GI = r_G - r_P/2
	r_GA = r_P
	r_PI = f_P*r_P
	r_PA = f_PA*(1 - f_P)*r_P
endif

endif

r_O2 = f_PO*(1 - f_P)*r_P
mp%A_rate = r_GA + r_PA
mp%I_rate = r_GI + r_PI
mp%P_rate = r_P
mp%O_rate = r_O2
mp%L_rate = 0
mp%G_rate = r_GI + 0.5*r_P	! this ensures mass conservation when fPDK < 1
mp%C_P = 0
mp%f_G = f_G
mp%f_P = f_P
! Add base rate correction
mp%O_rate = mp%O_rate + O2_baserate
mp%G_rate = mp%G_rate + G_baserate
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

!----------------------------------------------------------------------------------
subroutine Set_f_GP(mp,C)
type(metabolism_type), pointer :: mp
real(REAL_KIND) :: C(:)

if (fgp_solver == FGP_SOLVER_MAXATP_TANDEM) then
	call Set_f_GP_tandem(mp,C)
elseif (fgp_solver == FGP_SOLVER_MAXATP_STAGED) then
	call Set_f_GP_maxATP(mp,C)
elseif (fgp_solver == FGP_SOLVER_SURVIVAL_STAGED) then
	call Set_f_GP_survival(mp,C)
endif
end subroutine

!----------------------------------------------------------------------------------
! Assumes: 
!	always use_ATP
!	always use_lactate
!
! mp%A_rate = (1-f_G)*r_G*N_GA + (1-f_P)*r_P*N_PA	! production
! mp%I_rate = f_G*r_G*N_GI + f_P*r_P*N_PI		    ! production
! mp%P_rate = r_P								    ! utilisation
! mp%O_rate = f_PO*r_P*(1-f_P)					    ! consumption
! mp%L_rate = V*(K1*C_P - K2*C_L)				    ! production
!
! Note: moved from ODE_diffuse 5/2/18
!----------------------------------------------------------------------------------
subroutine Set_f_GP_survival(mp,C)
integer :: ityp
type(metabolism_type), pointer :: mp
real(REAL_KIND) :: C(:)
real(REAL_KIND) :: C_O2, C_G, C_L, Ofactor, Gfactor, ATPfactor, x, fmin
real(REAL_KIND) :: CO_L, CG_L, f_ATP_H, f_ATP_L, f, f_G_H, f_G_L
real(REAL_KIND) :: f_G, f_P, r_P, r_G, r_GI, r_GA, r_PI, r_PA, Km, C1, C2, r_L, r_A
integer :: N_O2, N_P, n
real(REAL_KIND) :: Km_O2, Km_P, V, f_PO, f_PA, K1, K2, MM_O2, fPDK, r_Pm_base, r_Ag, f_P_min, r_A_thresh
real(REAL_KIND) :: a1, b1, a2, b2, a3, b3, c3, d, fraction, C_P, r_O2
!type(metabolism_type) :: metab_save
!type(metabolism_type), pointer :: mp_save
real(REAL_KIND) :: alfa = 0.1
real(REAL_KIND) :: average_volume = 1.2
logical :: f_P_first = .true.	! First f_P is determined, then f_G (Note: f_P_first = .false. is not used)

if (N1D == 0) then
	! Cex for 3D case - vspheroid
	C_O2 = C(1)
	C_G = C(2)
	C_L = C(3)
else
	! Cex for 1D case - monolayer
	C_O2 = C(1)
	C_G = C(N1D+2)
	C_L = C(2*N1D+3)
endif

N_O2 = Hill_N_O2
Km_O2 = Hill_Km_O2
N_P = 1
Km_P = Hill_Km_P
V = Vcell_cm3*average_volume		! should be actual cell volume cp%V 
f_PO = N_PO
f_PA = N_PA
K1 = K_PL
K2 = K_LP
r_Ag = ATPg
r_G = get_glycosis_rate(mp%HIF1,C_G)
fPDK = mp%PDK1
MM_O2 = f_MM(C_O2,Km_O2,N_O2)
r_Pm_base = fPDK*MM_O2*O2_maxrate/f_PO	! note that MM_P is not here, since it varies it is added as needed

r_A_thresh = r_Ag + 0.5*(r_A_norm - r_Ag)    ! TESTING

if (f_P_first) then
	mp%f_G = f_G_norm
	mp%f_P = 0
	call f_metab(mp,C_O2,C_G,C_L)
!	write(nflog,'(a,4e12.3)') '(1) A_rate, r_Ag, f_ATPg, r_A_norm: ',mp%A_rate, r_Ag, f_ATPg, r_A_norm
	if (mp%A_rate < r_Ag) then
		! Need to set f_P = 0, reduce f_G below f_G_norm
		mp%f_G = 0
		call f_metab(mp,C_O2,C_G,C_L)
!		write(nflog,'(a,2e12.3)') '(2) A_rate, r_Ag: ',mp%A_rate, r_Ag
		if (mp%A_rate < r_Ag) then
!			write(nflog,*) 'Set f_P=0, f_G=0'
			! Need to set f_P = 0, f_G = 0
			f_P = 0
			f_G = 0
		else
			! Need to set f_P = 0, find level of f_G to maintain ATPg	*** (1)		! THIS IS OK
!			write(nflog,*) 'Set f_P=0, find f_G'
			f_P = 0
			a1 = 1./(N_GA*r_G*(1 + f_PA))
			a2 = a1*(r_Ag - V*f_PA*K2*C_L)
			b2 = a1*V*f_PA*K1
			a3 = N_GA*r_G*b2
			b3 = r_Pm_base*f_PA + N_GA*r_G*a2 + Km_P*N_GA*r_G*b2 - r_Ag
			c3 = Km_P*(N_GA*r_G*a2 - r_Ag)
			d = sqrt(b3*b3 - 4*a3*c3)
			C_P = (-b3 + d)/(2*a3)
			f_G = 1 - (a2 + b2*C_P)
		endif
		mp%f_G = f_G
		mp%f_P = f_P
		r_O2 = r_Pm_base*N_PO*C_P/(Km_P + C_P)
		r_P = r_O2/(N_PO*(1 - f_P))
!		write(*,'(a,4e12.3)') 'r_P, C_P, f_P, f_G: ',r_P, C_P, f_P, f_G
!		call f_metab(mp,C_O2,C_G,C_L)	! for checking f_G, f_P
!		if (mp%A_rate < r_Ag .and. mp%I_rate > 0) then !(mp%f_G > 0 .or. mp%f_P > 0)) then 
!			write(nflog,'(a,5e11.3)') 'A: r_A, r_Ag, f_G, f_P: ',mp%A_rate, r_Ag, mp%f_G, mp%f_P, mp%I_rate
!		endif
	else
		mp%f_P = f_P_norm
		call f_metab(mp,C_O2,C_G,C_L)
		! Need to set f_G = f_G_norm, find level of f_P to maintain ATPg *** (2)
		f_G = f_G_norm
		r_GA = N_GA*(1 - f_G)*r_G
		a1 = (r_GA - r_Ag)/(f_PA*V*K1)
		b1 = (r_GA + V*K2*C_L)/(V*K1)
		a2 = a1*((r_Ag - r_GA)/f_PA - r_Pm_base)
		b2 = r_Pm_base*b1 - (Km_P + b1)*(r_Ag - r_GA)/f_PA
		! z = 1/(1-f_P) = b2/a2, -> f_P = 1 - a2/b2
		f_P_min = 1 - a2/b2     ! this is f_P needed to reach r_Ag 
		if (mp%A_rate < r_Ag) then
			mp%f_G = f_G
			mp%f_P = f_P_min
		elseif (mp%A_rate < r_A_thresh) then
		    fraction = (mp%A_rate - r_Ag)/(r_A_thresh - r_Ag)
		    mp%f_P = min(f_P_min + fraction*(f_P_norm - f_P_min),f_P_norm)
			mp%f_G = f_G_norm
		else
!			write(nflog,*) 'Set f_G=f_Gn, f_P=f_Pn'
			! Need to set f_P = f_P_norm, f_G = f_G_norm
			mp%f_P = f_P_norm
			mp%f_G = f_G_norm
		endif
!		call f_metab(mp,C_O2,C_G,C_L)	! for checking f_G, f_P
!		if (mp%A_rate < r_Ag .and. mp%I_rate > 0) then !(mp%f_G > 0 .or. mp%f_P > 0)) then
!			write(nflog,'(a,5e11.3)') 'B: r_A, r_Ag, f_G, f_P: ',mp%A_rate, r_Ag, mp%f_G, mp%f_P, mp%I_rate
!		endif
	endif
else    ! f_G first
	mp%f_G = 0
	mp%f_P = f_P_norm
	call f_metab(mp,C_O2,C_G,C_L)
!	write(nflog,'(a,4e12.3)') '(1) A_rate, r_Ag, f_ATPg, r_A_norm: ',mp%A_rate, r_Ag, f_ATPg, r_A_norm
	if (mp%A_rate < r_Ag) then
		! Need to set f_G = 0, reduce f_P below f_P_norm
		mp%f_P = 0
		call f_metab(mp,C_O2,C_G,C_L)
!		write(nflog,'(a,2e12.3)') '(2) A_rate, r_Ag: ',mp%A_rate, r_Ag
		if (mp%A_rate < r_Ag) then
!			write(nflog,*) 'Set f_P=0, f_G=0'
			! Need to set f_P = 0, f_G = 0
			mp%f_P = 0
			mp%f_G = 0
		else
			! Need to set f_G = 0, find level of f_P to maintain ATPg	*** (1)
!			write(nflog,*) 'Set f_G=0, find f_P'
			f_G = 0
			r_GA = N_GA*(1 - f_G)*r_G
			a1 = (r_GA - r_Ag)/(f_PA*V*K1)
			b1 = (r_GA + V*K2*C_L)/(V*K1)
			a2 = a1*((r_Ag - r_GA)/f_PA - r_Pm_base)
			b2 = r_Pm_base*b1 - (Km_P + b1)*(r_Ag - r_GA)/f_PA
			! z = 1/(1-f_P) = b2/a2, -> f_P = 1 - a2/b2
			mp%f_G = f_G
			mp%f_P = 1 - a2/b2		
		endif
	else
		mp%f_G = f_G_norm
		call f_metab(mp,C_O2,C_G,C_L)
!		write(nflog,'(a,2e12.3)') '(3) A_rate, r_Ag: ',mp%A_rate, r_Ag
		if (mp%A_rate < r_Ag) then
!			write(nflog,*) 'Set f_P=f_Pn, find f_G'
			! Need to set f_P = f_P_norm, find level of f_G to maintain ATPg *** (2)
			f_P = f_P_norm
			a1 = N_GA*r_G/(f_PA*(1-f_P))
			b1 = (r_Ag - N_GA*r_G)/(f_PA*(1-f_P))
			a2 = (-N_GA*r_G - a1)/(V*K1)
			b2 = (N_GA*r_G + V*K2*C_L - b1)/(V*K1)
			a3 = a1*a2
			b3 = (a1*(Km_P + b2) + a2*b1 - r_Pm_base*a2/(1 - f_P))
			c3 = b1*(Km_P + b2) - r_Pm_base*b2/(1 - f_P)
			! a3.y^2 + b3.y + c3 = 0 where y = f_G
			d = sqrt(b3*b3 - 4*a3*c3)
			mp%f_P = f_P
			mp%f_G = (-b3 + d)/(2*a3)			
		else
!			write(nflog,*) 'Set f_G=f_Gn, f_P=f_Pn'
			! Need to set f_P = f_P_norm, f_G = f_G_norm
			mp%f_P = f_P_norm
			mp%f_G = f_G_norm
		endif
	endif
    f_PA = N_PA
    r_G = mp%G_rate
    r_P = mp%P_rate
endif
return

f_ATP_L = f_ATPg
f_ATP_H = min(f_ATPramp*f_ATPg,1.0)
f = mp%A_rate/r_A_norm
if (f > f_ATP_H) then
	ATPfactor = 1
elseif (f < f_ATP_L) then
	ATPfactor = 0
else
	ATPfactor = (f - f_ATP_L)/(f_ATP_H - f_ATP_L)
endif
x = mp%G_rate/r_G_norm
C1 = 1
C2 = 0
n = 1
if (x > C1) then
	Gfactor = 1
elseif (x < C2) then
	Gfactor = 0
else
	x = (x-C2)/(C1-C2)
	Gfactor = x**n
endif
C_O2 = C(1)
C1 = 0.10	! C_O2_norm
C2 = 0
n = 5
fmin = 0.3
if (C_O2 > C1) then
	Ofactor = 1 - fmin
elseif (C_O2 < C2) then
	Ofactor = 0
else
	x = (C_O2 - C2)/(C1 - C2)
	Ofactor = x**n*(1-fmin)
endif
Ofactor = Ofactor + fmin
f_P = Ofactor*f_P_norm
f_G = Ofactor*f_G_norm
ATPfactor = Gfactor*ATPfactor		!!!!!! TRY THIS !!!!!

mp%f_P = alfa*ATPfactor*f_P + (1-alfa)*mp%f_P
mp%f_G = alfa*ATPfactor*f_G + (1-alfa)*mp%f_G
write(nflog,'(a,4e12.3)') 'mp%A_rate/r_A_norm,ATPfactor,mp%f_G,mp%f_P: ',mp%A_rate/r_A_norm,ATPfactor,mp%f_G,mp%f_P
end subroutine

!----------------------------------------------------------------------------------
subroutine Set_f_GP_maxATP(mp,C)
integer :: ityp
type(metabolism_type), pointer :: mp
real(REAL_KIND) :: C(:)
real(REAL_KIND) :: C_O2, C_G, C_L
real(REAL_KIND) :: f_G, f_P, r_P, r_G, r_O2
integer :: N_O2, N_P
real(REAL_KIND) :: Km_O2, Km_P, V, K1, K2, MM_O2, fPDK, r_Pm_base, r_Ag, r_An
real(REAL_KIND) :: a1, b1, a2, b2, a3, b3, c3, d, C_P, x, y, z, zz
real(REAL_KIND) :: average_volume = 1.2
logical :: f_P_first = .true.	! First f_P is determined, then f_G (Note: f_P_first = .false. is not used)

if (N1D == 0) then
	! Cex for 3D case - vspheroid
	C_O2 = C(1)
	C_G = C(2)
	C_L = C(3)
else
	! Cex for 1D case - monolayer
	C_O2 = C(1)
	C_G = C(N1D+2)
	C_L = C(2*N1D+3)
endif

N_O2 = Hill_N_O2
Km_O2 = Hill_Km_O2
Km_P = Hill_Km_P
V = Vcell_cm3*average_volume		! should be actual cell volume cp%V 
K1 = K_PL
K2 = K_LP
r_Ag = ATPg
r_An = r_A_norm

r_G = get_glycosis_rate(mp%HIF1,C_G)
fPDK = mp%PDK1
MM_O2 = f_MM(C_O2,Km_O2,N_O2)
r_Pm_base = fPDK*MM_O2*O2_maxrate/N_PO	! note that MM_P is not here, since it varies it is added as needed

f_G = f_G_norm
mp%f_G = f_G
mp%f_P = 0
call f_metab(mp,C_O2,C_G,C_L)
if (mp%A_rate < r_An) then
	f_P = 0
	! Now evaluate f_G < f_G_norm
	a1 = 1/(N_GA*r_G*(1 + N_PA))
	a2 = (r_An - V*K2*C_L*N_PA)*a1
	b2 = V*K1*N_PA*a1
	a3 = N_GA*r_G*b2
	b3 = r_Pm_base*N_PA + N_GA*r_G*a2 + Km_P*N_GA*r_G*b2 - r_An
	c3 = Km_P*(N_GA*r_G*a2 - r_An)
	d = sqrt(b3*b3 - 4*a3*c3)
	C_P = (-b3 + d)/(2*a3)
	y = a2 + b2*C_P
	f_G = max(0.0, 1-y)
else
	! leave f_G = f_G_norm, evaluate f_P < f_P_norm
	y = 1 - f_G
	zz = (r_An - y*r_G*N_GA)/(r_Pm_base*N_PA)
	z = Km_P*zz/(1 - zz)
	x = r_Pm_base*z/((Km_P + z)*(y*r_G*N_GA - V*(k1*z - K2*C_L)))
	f_P = 1.0 - x
	C_P = z
!	write(*,'(a,5e12.3)') 'zz,z,x,f_P,C_P: ',zz,z,x,f_P,C_P
	r_O2 = r_Pm_base*N_PO*C_P/(Km_P + C_P)
	r_P = r_O2/(N_PO*(1 - f_P))
endif

#if 0
if (f_P > f_P_norm) then
	f_P = f_P_norm
endif
if (f_P < 0) then
	f_P = 0
	! Now evaluate f_G
	a1 = 1/(N_GA*r_G*(1 + N_PA))
	a2 = (r_An - V*K2*C_L*N_PA)*a1
	b2 = V*K1*N_PA*a1
	a3 = N_GA*r_G*b2
	b3 = r_Pm_base*N_PA + N_GA*r_G*a2 + Km_P*N_GA*r_G*b2 - r_An
	c3 = Km_P*(N_GA*r_G*a2 - r_An)
	d = sqrt(b3*b3 - 4*a3*c3)
	C_P = (-b3 + d)/(2*a3)
	y = a2 + b2*C_P
	f_G = max(1-y, 0.0)
endif
#endif
!write(*,'(a,5e12.3)') 'r_P, r_O2, C_P, f_P, f_G: ',r_P, r_O2, C_P, f_P, f_G
mp%f_G = f_G
mp%f_P = f_P
end subroutine

!----------------------------------------------------------------------------------
subroutine Set_f_GP_tandem(mp,C)
integer :: ityp
type(metabolism_type), pointer :: mp
real(REAL_KIND) :: C(:)
real(REAL_KIND) :: C_O2, C_G, C_L
real(REAL_KIND) :: f_G, f_P, fGn, fPn, r_P, r_G, r_O2
integer :: N_O2, N_P, n, i
real(REAL_KIND) :: Km_O2, Km_P, V, K1, K2, MM_O2, fPDK, r_Pm_base, r_An
real(REAL_KIND) :: a1, b1, c1, a2, b2, c2, a3, b3, c3, e1, e2, r(3), w, z, C_P, val
real(REAL_KIND) :: average_volume = 1.2
logical :: f_P_first = .true.	! First f_P is determined, then f_G (Note: f_P_first = .false. is not used)

integer :: ierr, Ndim, k, ifail
real(REAL_KIND) :: x(1), Rtol, f(1), x0, x1, dx, f0, f1, dfdx

if (N1D == 0) then
	! Cex for 3D case - vspheroid
	C_O2 = C(1)
	C_G = C(2)
	C_L = C(3)
else
	! Cex for 1D case - monolayer
	C_O2 = C(1)
	C_G = C(N1D+2)
	C_L = C(2*N1D+3)
endif

N_O2 = Hill_N_O2
Km_O2 = Hill_Km_O2
Km_P = Hill_Km_P
V = Vcell_cm3*average_volume		! should be actual cell volume cp%V 
K1 = K_PL
K2 = K_LP
r_An = r_A_norm

r_G = get_glycosis_rate(mp%HIF1,C_G)
fPDK = mp%PDK1
MM_O2 = f_MM(C_O2,Km_O2,N_O2)
r_Pm_base = fPDK*MM_O2*O2_maxrate/N_PO	! note that MM_P is not here, since it varies it is added as needed

fGn = f_G_norm
fPn = f_P_norm
f_G = fGn
f_P = fPn
mp%f_G = f_G
mp%f_P = f_P
call f_metab(mp,C_O2,C_G,C_L)
if (mp%A_rate < r_An) then
	e1 = r_An/(r_Pm_base*N_PA)
	e2 = r_G*N_GA/(r_Pm_base*N_PA)
	a1 = (e2 - e1)*Km_P/(fGn*e2)
	b1 = (e2 - e1 + 1)/(fGn*e2)
	c1 = Km_P
	! a1, b1, c1 look OK, give w = 1 for steady-state C_P
!	z = 0.0568
!	w = (a1 + b1*z)/(c1 + z)
!	write(*,*) 'z,w: ',z,w
!	w = 1
!	val = r_Pm_base*z/((1 - w*fPn)*(Km_P+z)) - (1-w*fGn)*r_G*N_GA + V*(K1*z - K2*C_L)
!	write(*,*) 'val: ',val
	! This equation is OK
		
	rtol = 1.0d-6
!	x0 = 0.06
	x0 = mp%C_P
	dx = x0/1000
	do k = 1,100
!		x(1) = x0
!		call FCN(1,X,f,IFAIL)
!		f0 = f(1)
		w = (a1 + b1*x0)/(c1 + x0)
		w = min(w,1.0)
		f0 = r_Pm_base*x0/((1 - w*fPn)*(Km_P+x0)) - (1-w*fGn)*r_G*N_GA + V*(K1*x0 - K2*C_L)
!		x(1) = x0 + dx
!		call FCN(1,X,f,IFAIL)
!		f1 = f(1)
		x1 = x0 + dx
		w = (a1 + b1*x1)/(c1 + x1)
		w = min(w,1.0)
		f1 = r_Pm_base*x1/((1 - w*fPn)*(Km_P+x1)) - (1-w*fGn)*r_G*N_GA + V*(K1*x1 - K2*C_L)
		dfdx = (f1 - f0)/dx
		x0 = x0 - f0/dfdx
		if (abs(f0/dfdx) < rtol) exit
	enddo
	z = x1

#if 0
	FCN_a1 = a1
	FCN_b1 = b1
	FCN_c1 = c1
	FCN_r_Pm_base = r_Pm_base
	FCN_fPn = fPn
	FCN_Km_P = Km_P
	FCN_fGn = fGn
	FCN_r_G = r_G
	FCN_N_GA = N_GA
	FCN_V = V
	FCN_K1 = K1
	FCN_K2 = K2
	FCN_C_L = C_L
	Ndim = 1
	rtol = 1.0d-5
	x(1) = 0.06	! initial guess
	CALL NLEQ1E(FCN,dummyjac,Ndim,x,rtol,ierr)
	z = x(1)
#endif

	w = (a1 + b1*z)/(c1 + z)
!	write(*,*) 'k, f1, z, w: ',k,f1,z,w
	w = max(w,0.0)
	f_G = w*fGn
	f_P = w*fPn
endif
mp%f_G = f_G
mp%f_P = f_P

!write(*,'(a,5e12.3)') 'r_P, r_O2, C_P, f_P, f_G: ',r_P, r_O2, C_P, f_P, f_G
end subroutine

#if 0
!---------------------------------------------------------------------------
subroutine FCN(N,X,F,IFAIL)
use FCN_mod
integer :: N, IFAIL
real(8) :: X(N),F(N)
real(8) :: z, w

z = x(1)
w = (FCN_a1 + FCN_b1*z)/(FCN_c1 + z)
F(1) = FCN_r_Pm_base*z/((1 - w*FCN_fPn)*(FCN_Km_P+z)) - (1-w*FCN_fGn)*FCN_r_G*FCN_N_GA + FCN_V*(FCN_K1*z - FCN_K2*FCN_C_L)
end subroutine

!---------------------------------------------------------------------------
subroutine dummyjac(N,LDJAC,X,DFDX,IFAIL)
integer :: N, LDJAC, IFAIL
real(8) :: X(N), DFDX(LDJAC,N)
end subroutine

!----------------------------------------------------------------------------------
subroutine Set_f_GP_noSS(mp,C,C_P)
integer :: ityp
type(metabolism_type), pointer :: mp
real(REAL_KIND) :: C(:), C_P
real(REAL_KIND) :: C_O2, C_G, C_L, Ofactor, Gfactor, ATPfactor, x, fmin
real(REAL_KIND) :: CO_L, CG_L, f_ATP_H, f_ATP_L, f, f_G_H, f_G_L
real(REAL_KIND) :: f_G, f_P, r_P, r_G, r_GI, r_GA, r_PI, r_PA, Km, C1, C2, r_L, r_A
integer :: N_O2, N_P, n
real(REAL_KIND) :: Km_O2, Km_P, V, f_PO, f_PA, K1, K2, MM_O2, fPDK, r_Pm_base, r_Ag
real(REAL_KIND) :: a1, b1, a2, b2, a3, b3, c3, d
type(metabolism_type) :: metab_save
type(metabolism_type), pointer :: mp_save
real(REAL_KIND) :: alfa = 0.1
logical :: f_P_first = .true.	! First f_P is determined, then f_G

if (N1D == 0) then
	! Cex for 3D case - vspheroid
	C_O2 = C(1)
	C_G = C(2)
	C_L = C(3)
else
	! Cex for 1D case - monolayer
	C_O2 = C(1)
	C_G = C(N1D+2)
	C_L = C(2*N1D+3)
endif

N_O2 = Hill_N_O2
Km_O2 = Hill_Km_O2
N_P = 1
Km_P = Hill_Km_P
V = Vcell_cm3		! should be actual cell volume cp%V 
f_PO = N_PO
f_PA = N_PA
K1 = K_PL
K2 = K_LP
r_Ag = ATPg
r_G = get_glycosis_rate(mp%HIF1,C_G)
fPDK = mp%PDK1
MM_O2 = f_MM(C_O2,Km_O2,N_O2)
r_Pm_base = fPDK*MM_O2*O2_maxrate/f_PO	! note that MM_P is not here, since it varies it is added as needed

if (f_P_first) then
	mp%f_G = f_G_norm
	mp%f_P = 0
	call f_metab_noSS(mp,C_O2,C_G,C_L,C_P)
!	write(nflog,'(a,4e12.3)') '(1) A_rate, r_Ag, f_ATPg, r_A_norm: ',mp%A_rate, r_Ag, f_ATPg, r_A_norm
	if (mp%A_rate < r_Ag) then
		! Need to set f_P = 0, reduce f_G below f_G_norm
		mp%f_G = 0
		call f_metab_noSS(mp,C_O2,C_G,C_L,C_P)
!		write(nflog,'(a,2e12.3)') '(2) A_rate, r_Ag: ',mp%A_rate, r_Ag
		if (mp%A_rate < r_Ag) then
!			write(nflog,*) 'Set f_P=0, f_G=0'
			! Need to set f_P = 0, f_G = 0
			mp%f_P = 0
			mp%f_G = 0
		else
			! Need to set f_P = 0, find level of f_G to maintain ATPg	*** (1)		! THIS IS OK
!			write(nflog,*) 'Set f_P=0, find f_G'
			f_P = 0
			a1 = N_GA*r_G/(f_PA*(1-f_P))
			b1 = (r_Ag - N_GA*r_G)/(f_PA*(1-f_P))
			a2 = (-N_GA*r_G - a1)/(V*K1)
			b2 = (N_GA*r_G + V*K2*C_L - b1)/(V*K1)
			a3 = a1*a2
			b3 = (a1*(Km_P + b2) + a2*b1 - r_Pm_base*a2/(1 - f_P))
			c3 = b1*(Km_P + b2) - r_Pm_base*b2/(1 - f_P)
			! a3.y^2 + b3.y + c3 = 0 where y = f_G
			d = sqrt(b3*b3 - 4*a3*c3)
			mp%f_P = f_P
			mp%f_G = (-b3 + d)/(2*a3)
		endif
	else
		mp%f_P = f_P_norm
		call f_metab_noSS(mp,C_O2,C_G,C_L,C_P)
!		write(nflog,'(a,2e12.3)') '(3) A_rate, r_Ag: ',mp%A_rate, r_Ag
		if (mp%A_rate < r_Ag) then
!			write(nflog,*) 'Set f_G=f_Gn, find f_P'
			! Need to set f_G = f_G_norm, find level of f_P to maintain ATPg *** (2)
			f_G = f_G_norm
			r_GA = N_GA*(1 - f_G)*r_G
			a1 = (r_GA - r_Ag)/(f_PA*V*K1)
			b1 = (r_GA + V*K2*C_L)/(V*K1)
			a2 = a1*((r_Ag - r_GA)/f_PA - r_Pm_base)
			b2 = r_Pm_base*b1 - (Km_P + b1)*(r_Ag - r_GA)/f_PA
			! z = 1/(1-f_P) = b2/a2, -> f_P = 1 - a2/b2
			mp%f_G = f_G
			mp%f_P = 1 - a2/b2
		else
!			write(nflog,*) 'Set f_G=f_Gn, f_P=f_Pn'
			! Need to set f_P = f_P_norm, f_G = f_G_norm
			mp%f_P = f_P_norm
			mp%f_G = f_G_norm
		endif
	endif
else    ! f_G first
	mp%f_G = 0
	mp%f_P = f_P_norm
	call f_metab_noSS(mp,C_O2,C_G,C_L,C_P)
!	write(nflog,'(a,4e12.3)') '(1) A_rate, r_Ag, f_ATPg, r_A_norm: ',mp%A_rate, r_Ag, f_ATPg, r_A_norm
	if (mp%A_rate < r_Ag) then
		! Need to set f_G = 0, reduce f_P below f_P_norm
		mp%f_P = 0
		call f_metab_noSS(mp,C_O2,C_G,C_L,C_P)
!		write(nflog,'(a,2e12.3)') '(2) A_rate, r_Ag: ',mp%A_rate, r_Ag
		if (mp%A_rate < r_Ag) then
!			write(nflog,*) 'Set f_P=0, f_G=0'
			! Need to set f_P = 0, f_G = 0
			mp%f_P = 0
			mp%f_G = 0
		else
			! Need to set f_G = 0, find level of f_P to maintain ATPg	*** (1)
!			write(nflog,*) 'Set f_G=0, find f_P'
			f_G = 0
			r_GA = N_GA*(1 - f_G)*r_G
			a1 = (r_GA - r_Ag)/(f_PA*V*K1)
			b1 = (r_GA + V*K2*C_L)/(V*K1)
			a2 = a1*((r_Ag - r_GA)/f_PA - r_Pm_base)
			b2 = r_Pm_base*b1 - (Km_P + b1)*(r_Ag - r_GA)/f_PA
			! z = 1/(1-f_P) = b2/a2, -> f_P = 1 - a2/b2
			mp%f_G = f_G
			mp%f_P = 1 - a2/b2			
		endif
	else
		mp%f_G = f_G_norm
		call f_metab_noSS(mp,C_O2,C_G,C_L,C_P)
!		write(nflog,'(a,2e12.3)') '(3) A_rate, r_Ag: ',mp%A_rate, r_Ag
		if (mp%A_rate < r_Ag) then
!			write(nflog,*) 'Set f_P=f_Pn, find f_G'
			! Need to set f_P = f_P_norm, find level of f_G to maintain ATPg *** (2)
			f_P = f_P_norm
			a1 = N_GA*r_G/(f_PA*(1-f_P))
			b1 = (r_Ag - N_GA*r_G)/(f_PA*(1-f_P))
			a2 = (-N_GA*r_G - a1)/(V*K1)
			b2 = (N_GA*r_G + V*K2*C_L - b1)/(V*K1)
			a3 = a1*a2
			b3 = (a1*(Km_P + b2) + a2*b1 - r_Pm_base*a2/(1 - f_P))
			c3 = b1*(Km_P + b2) - r_Pm_base*b2/(1 - f_P)
			! a3.y^2 + b3.y + c3 = 0 where y = f_G
			d = sqrt(b3*b3 - 4*a3*c3)
			mp%f_P = f_P
			mp%f_G = (-b3 + d)/(2*a3)			
		else
!			write(nflog,*) 'Set f_G=f_Gn, f_P=f_Pn'
			! Need to set f_P = f_P_norm, f_G = f_G_norm
			mp%f_P = f_P_norm
			mp%f_G = f_G_norm
		endif
	endif
    f_PA = N_PA
    r_G = mp%G_rate
    r_P = mp%P_rate
!    write(*,'(a,4e12.3)') 'f_G,f_P: ',mp%f_G,mp%f_P
!    write(*,'(a,4e12.3)') 'A_rate,I_rate: ',2*(1-mp%f_G)*r_G, f_PA*(1-mp%f_P)*r_P,mp%f_G*r_G,mp%f_P*r_P
endif
!write(nflog,'(a,2e12.3)') 'f_P,f_G: ',mp%f_P,mp%f_G
!write(nflog,*)
end subroutine

!--------------------------------------------------------------------------
! Test metab rates with very low O2, holding glucose and lactate constant.
!--------------------------------------------------------------------------
subroutine testmetab1
type(metabolism_type), target :: metab
type(metabolism_type), pointer :: mp
integer :: ityp, i
real(REAL_KIND) :: C_O2, C_G, C_L

metab = cell_list(1)%metab
mp => metab
ityp = 1
C_G = 100
C_L = 18
C_O2 = 0.1
mp%HIF1 = 1
mp%PDK1 = 0.6
write(nflog,*) 'i  C_O2  C_L   A_rate   I_rate   P_rate   O_rate   HIF1 PDK1   C_P'
do i = 1,100
	call f_metab(mp, C_O2, C_G, C_L)
	write(nflog,'(i6,9e12.3)') i,C_O2,C_L,mp%A_rate,mp%I_rate,mp%P_rate,mp%O_rate,mp%HIF1,mp%PDK1,mp%C_P
	C_O2 = 0.8*C_O2
!	C_O2 = C_O2 - 0.00005
enddo
stop
end subroutine

!--------------------------------------------------------------------------
! Test metab rates with a range of O2, holding glucose and lactate constant.
! Note: lactate conc C_L has no effect.  Is this correct?
!--------------------------------------------------------------------------
subroutine testmetab2
type(metabolism_type), target :: metab
type(metabolism_type), pointer :: mp
integer :: ityp, i, j
real(REAL_KIND) :: C_O2, C_G, C_L, C_P
real(REAL_KIND) :: C(3*N1D+3)
logical :: single_case = .false.

metab = cell_list(1)%metab
mp => metab
ityp = 1

write(*,*) 'ATPg: ',ATPg
write(nflog,*) 'ATPg: ',ATPg
write(nflog,'(a)') '     i        C_O2         C_L      G_rate      A_rate      I_rate      P_rate      L_rate      O_rate        HIF1        PDK1'
if (single_case) then
! To check low-glucose metabolism
C_G = 0.005
C_L = 3.6
C_P = 0.06
C_O2 = 0.15
C(1) = C_O2
C(N1D+2) = C_G
C(2*N1D+3) = C_L
mp%HIF1 = get_HIF1steadystate(C_O2)
mp%PDK1 = 1 - (1 - PDKmin)*mp%HIF1  ! steady-state PDK1?
call set_f_GP_noSS(mp,C,C_P)
call f_metab_noSS(mp, C_O2, C_G, C_L, C_P)
write(nflog,'(i6,10e12.3)') i,C_O2,C_L,mp%G_rate,mp%A_rate,mp%I_rate,mp%P_rate,mp%L_rate,mp%O_rate,mp%HIF1,mp%PDK1
stop
endif

C_G = 5.5
C_P = 0.01
j = 1
!do j = 1,11
	C_L = (j-1)*0.3
	C_O2 = 0.15
!	write(*,*) 'ATPg: ',ATPg
!	write(nflog,*) 'ATPg: ',ATPg
!	write(nflog,'(a)') '     i        C_G         C_L      G_rate      A_rate      I_rate      P_rate      L_rate      O_rate        HIF1        PDK1'
	do i = 1,50
		C(1) = C_O2
		C(N1D+2) = C_G
		C(2*N1D+3) = C_L
		mp%HIF1 = get_HIF1steadystate(C_O2)
		mp%PDK1 = 1 - (1 - PDKmin)*mp%HIF1  ! steady-state PDK1?
		call set_f_GP_noSS(mp,C,C_P)
		call f_metab_noSS(mp, C_O2, C_G, C_L, C_P)
		write(nflog,'(i6,10e12.3)') i,C_G,C_L,mp%G_rate,mp%A_rate,mp%I_rate,mp%P_rate,mp%L_rate,mp%O_rate,mp%HIF1,mp%PDK1
		C_G = 0.8*C_G
	enddo
!enddo
stop
end subroutine

!--------------------------------------------------------------------------
! Test metab rates with a range of glucose, holding oxygen and lactate constant.
!--------------------------------------------------------------------------
subroutine testmetab3
type(metabolism_type), target :: metab
type(metabolism_type), pointer :: mp
integer :: ityp, i, j
real(REAL_KIND) :: C_O2, C_G, C_L
real(REAL_KIND) :: C(3*N1D+3)

metab = cell_list(1)%metab
mp => metab
ityp = 1
C_O2 = 0.15
C_L = 0
!do j = 1,11
!	C_L = (j-1)*0.3
	C_G = 5.5
	write(nflog,*)
	write(nflog,'(a)') '     i        C_G          C_L      G_rate      A_rate      I_rate      P_rate      L_rate      O_rate        HIF1         f_G         f_P'
	do i = 1,100
		C(1) = C_O2
		C(N1D+2) = C_G
		C(2*N1D+3) = C_L
		mp%HIF1 = get_HIF1steadystate(C_O2)
		mp%PDK1 = 1 - (1 - PDKmin)*mp%HIF1
		call set_f_GP(mp,C)
		call f_metab(mp, C_O2, C_G, C_L)
		write(nflog,'(i6,11e12.3)') i,C_G,C_L,mp%G_rate,mp%A_rate,mp%I_rate,mp%P_rate,mp%L_rate,mp%O_rate,mp%HIF1, mp%f_G, mp%f_P
		C_G = 0.8*C_G
	enddo
!enddo
stop
end subroutine
#endif

!----------------------------------------------------------------------------------
! Computes metabolism rate as a fraction of the maximum cell rate
! Use the "soft landing" option for Hill_N = 1 if MM_threshold = 0
!----------------------------------------------------------------------------------
function O2_metab(C) result(metab)
integer :: ichemo
real(REAL_KIND) :: C
real(REAL_KIND) :: metab

ichemo = OXYGEN
if (ichemo == OXYGEN) then
	if (chemo(ichemo)%Hill_N == 2) then
		if (C > 0) then
			metab = C*C/(chemo(ichemo)%MM_C0*chemo(ichemo)%MM_C0 + C*C)
		else
			metab = 0
		endif
	else
!		if (MM_THRESHOLD > 0) then
!			if (C > ODEdiff%C1_soft) then
!				metab = (C-ODEdiff%deltaC_soft)/(chemo(ichemo)%MM_C0 + C - ODEdiff%deltaC_soft)
!			elseif (C > 0) then
!				metab = ODEdiff%k_soft*C*C
!			else
!				metab = 0
!			endif
!		else
			if (C > 0) then
				metab = C/(chemo(ichemo)%MM_C0 + C)
			else
				metab = 0
			endif
!		endif
	endif
endif
end function

!----------------------------------------------------------------------------------
! Computes metabolism rate as a fraction of the maximum cell rate
!----------------------------------------------------------------------------------
function glucose_metab(C) result(metab)
real(REAL_KIND) :: C, metab
real(REAL_KIND) :: Kmin, Kmax, Km1
real(REAL_KIND) :: Vmax1, Vmax2, Km2, n1, n2
real(REAL_KIND) :: fV = 0.6
real(REAL_KIND) :: fK = 0.08
real(REAL_KIND) :: fboost = 2
real(REAL_KIND) :: Cboost = 0.1
logical :: variable_Km = .false.
logical :: double_Km = .false.
logical :: use_boost = .false.

if (C == 0) then
	metab = 0
	return
endif
if (use_boost) then

elseif (double_Km) then
	Km1 = Hill_Km_G
	Km2 = fK*Km1
	n1 = Hill_N_G
	n2 = 1
	metab = fV*C**n1/(Km1**n1 + C**n1) + (1 - fV)*C**n2/(Km2**n2 + C**n2)
elseif (variable_Km) then
	Kmax = Hill_Km_G	! These are completely arbitrary values
	Kmin = Kmax/15
	Km1 = 1*Kmin
	metab = C*(Km1 + C)/(Kmin*Km1 + Kmax*C + C*(Km1 + C))
else
	metab = C**Hill_N_G /(C**Hill_N_G + Hill_Km_G**Hill_N_G)
endif
end function

#if EXCEL
!--------------------------------------------------------------------------------------
! Determine real roots r(:) of the cubic equation:
! x^3 + a.x^2 + b.x + c = 0
! If there is one real root, n=1 and the root is r(1)
! If there are three distinct real roots, n=3 and the roots are r(1), r(2), r(3)
! If there is a repeated root, n=2 and the single root is r(1), the repeated root is r(2)
!--------------------------------------------------------------------------------------
subroutine cubic_roots(a, b, c, r, n)
real(REAL_KIND) :: a, b, c, r(3)
integer :: n
real(REAL_KIND) :: QQ, RR, theta, R2, Q3, AA, BB

QQ = (a*a - 3*b)/9
RR = (2*a*a*a - 9*a*b + 27*c)/54
Q3 = QQ*QQ*QQ
R2 = RR*RR
if (R2 < Q3) then
	n = 3
	theta = acos(RR/sqrt(Q3))
	r(1) = -2*sqrt(QQ)*cos(theta/3) - a/3
	r(2) = -2*sqrt(QQ)*cos((theta+2*PI)/3) - a/3
	r(3) = -2*sqrt(QQ)*cos((theta-2*PI)/3) - a/3
else
	n = 1
	AA = -sign(1.d0,RR)*(abs(RR) + sqrt(R2 - Q3))**(1.d0/3.d0)
	if (AA == 0) then
		BB = 0
	else
		BB = QQ/AA
	endif
	r(1) = AA + BB - a/3
endif
end subroutine
#endif

end module


