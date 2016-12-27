!---------------------------------------------------------------------------------------------------------------
! Force-based simulation
!
! UNITS: (OLD)
!    distance um
!    time     seconds
! Units:
!     time				s = seconds
!     distance			cm
!     volume			cm^3
!     mass				micromole = 10^-6 mol = mumol
!     flux				mumol/s
!     concentration		mumol/cm^3 = mM 
!---------------------------------------------------------------------------------------------------------------
module vspheroid

use global
use fmotion
use nbr
use chemokine
use react_diff
use transfer
use cellstate
use packer
use colony
use Tcp_mod

#include "../src/version.h"

implicit none

contains

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine ArrayInitialisation(ok)
logical :: ok

call RngInitialisation
! These are deallocated here instead of in subroutine wrapup so that when a simulation run ends 
! it will still be possible to view the cell distributions and chemokine concentration fields.
if (allocated(cell_list)) deallocate(cell_list)
if (allocated(grid)) deallocate(grid)
if (allocated(perm_index)) deallocate(perm_index)
if (allocated(Cextra_all)) deallocate(Cextra_all)
if (allocated(Caverage)) deallocate(Caverage)
if (allocated(Cflux)) deallocate(Cflux)
if (allocated(Cflux_prev)) deallocate(Cflux_prev)
if (allocated(stencil)) deallocate(stencil)
if (allocated(gaplist)) deallocate(gaplist)
call logger('did deallocation')
allocate(cell_list(MAX_NLIST))
allocate(gaplist(max_ngaps))
allocate(grid(NX,NY,NZ))
allocate(perm_index(MAX_NLIST))
allocate(Cextra_all(NX,NY,NZ,NCONST))
allocate(Caverage(NX,NY,NZ,NCONST))
allocate(Cflux(NX,NY,NZ,NCONST))
allocate(Cflux_prev(NX,NY,NZ,NCONST))
allocate(stencil(NX,NY,NZ))

ncells = 0
nlist = 0
ngaps = 0
ok = .true.

end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine ReadCellParams(ok)
logical :: ok
integer :: i, idrug, nmetab, im, ichemo
integer :: itestcase, ictype, ishow_progeny
integer :: iuse_oxygen, iuse_glucose, iuse_tracer, iuse_drug, iuse_metab, iV_depend, iV_random, iuse_FD
integer :: iuse_extra, iuse_relax, iuse_par_relax, iuse_gd_all
real(REAL_KIND) :: days, percent, fluid_fraction, d_layer, sigma(MAX_CELLTYPES), Vsite_cm3, bdry_conc, spcrad_value, d_n_limit
real(REAL_KIND) :: anoxia_tag_hours, anoxia_death_hours, aglucosia_tag_hours, aglucosia_death_hours
integer :: iuse_drop, isaveprofiledata, isaveslicedata, iusecellcycle, iusemetabolism
integer :: iconstant, ioxygengrowth, iglucosegrowth, ioxygendeath, iglucosedeath
logical :: use_metabolites
character*(12) :: drug_name
character*(1) :: numstr
type(cycle_parameters_type),pointer :: ccp

ok = .true.

open(nfcell,file=inputfile,status='old')

read(nfcell,'(a)') header
if (header(1:3) == 'GUI') then
	gui_run_version = header
	header = 'DD/MM/YYYY header_string'
else
	read(nfcell,*) gui_run_version				! program run version number
endif
read(nfcell,*) dll_run_version				! DLL run version number
read(nfcell,*) NX							! size of fine grid
NX = 33
read(nfcell,*) initial_count				! initial number of tumour cells
read(nfcell,*) divide_time_median(1)
read(nfcell,*) divide_time_shape(1)
read(nfcell,*) divide_time_median(2)
read(nfcell,*) divide_time_shape(2)
read(nfcell,*) iV_depend
read(nfcell,*) iV_random
read(nfcell,*) days							! number of days to simulate
read(nfcell,*) d_n_limit					! possible limit on diameter or number of cells
diam_count_limit = d_n_limit
read(nfcell,*) DELTA_T						! time step size (sec)
read(nfcell,*) n_substeps					! # of time step subdivisions for drug solving
read(nfcell,*) NXB							! size of coarse grid = NYB
read(nfcell,*) NZB							! size of coarse grid
read(nfcell,*) DELTA_X						! grid size (um)
read(nfcell,*) a_separation
read(nfcell,*) a_force
read(nfcell,*) c_force
read(nfcell,*) x0_force
read(nfcell,*) x1_force
read(nfcell,*) kdrag
read(nfcell,*) frandom
read(nfcell,*) NT_CONC						! number of subdivisions of DELTA_T for diffusion computation
read(nfcell,*) Nmm3							! number of cells/mm^3
read(nfcell,*) fluid_fraction				! fraction of the (non-necrotic) tumour that is fluid
read(nfcell,*) medium_volume0				! initial total volume (medium + spheroid) (cm^3)
read(nfcell,*) d_layer						! thickness of the unstirred layer around the spheroid (cm)
read(nfcell,*) Vdivide0						! nominal cell volume multiple for division
read(nfcell,*) dVdivide						! variation about nominal divide volume
read(nfcell,*) MM_THRESHOLD					! O2 concentration threshold Michaelis-Menten "soft-landing" (uM)
read(nfcell,*) anoxia_threshold			    ! O2 threshold for anoxia (uM)
read(nfcell,*) anoxia_tag_hours				! anoxia time leading to tagging to die by anoxia (h)
read(nfcell,*) anoxia_death_hours			! time after tagging to death by anoxia (h)
read(nfcell,*) aglucosia_threshold			! glucose threshold for aglucosia (uM)
read(nfcell,*) aglucosia_tag_hours			! aglucosia time leading to tagging to die by aglucosia (h)
read(nfcell,*) aglucosia_death_hours		! time after tagging to death by aglucosia (h)
read(nfcell,*) itestcase                    ! test case to simulate
read(nfcell,*) seed(1)						! seed vector(1) for the RNGs
read(nfcell,*) seed(2)						! seed vector(2) for the RNGs
read(nfcell,*) ncpu_input					! for GUI just a placeholder for ncpu, used only when execute parameter ncpu = 0
read(nfcell,*) Ncelltypes					! maximum number of cell types in the spheroid
do ictype = 1,Ncelltypes
	read(nfcell,*) percent
	celltype_fraction(ictype) = percent/100
enddo
read(nfcell,*) NT_GUI_OUT					! interval between GUI outputs (timesteps) 
read(nfcell,*) show_progeny                 ! if != 0, the number of the cell to show descendents of

read(nfcell,*) iuse_oxygen		! chemo(OXYGEN)%used
read(nfcell,*) ioxygengrowth
chemo(OXYGEN)%controls_growth = (ioxygengrowth == 1)
read(nfcell,*) ioxygendeath
chemo(OXYGEN)%controls_death = (ioxygendeath == 1)
read(nfcell,*) chemo(OXYGEN)%diff_coef
read(nfcell,*) chemo(OXYGEN)%medium_diff_coef
read(nfcell,*) chemo(OXYGEN)%membrane_diff_in
read(nfcell,*) chemo(OXYGEN)%membrane_diff_out
Vsite_cm3 = 2.0e-9		! this is from spheroid-abm, and the scaling here is for consistency with spheroid-abm
chemo(OXYGEN)%membrane_diff_in = chemo(OXYGEN)%membrane_diff_in*Vsite_cm3/60		! /min -> /sec
chemo(OXYGEN)%membrane_diff_out = chemo(OXYGEN)%membrane_diff_out*Vsite_cm3/60		! /min -> /sec
read(nfcell,*) chemo(OXYGEN)%bdry_conc
read(nfcell,*) iconstant
chemo(OXYGEN)%constant = (iconstant == 1)
read(nfcell,*) chemo(OXYGEN)%max_cell_rate
chemo(OXYGEN)%max_cell_rate = chemo(OXYGEN)%max_cell_rate*1.0e6					! mol/cell/s -> mumol/cell/s
read(nfcell,*) chemo(OXYGEN)%MM_C0
read(nfcell,*) chemo(OXYGEN)%Hill_N
read(nfcell,*) iuse_glucose		!chemo(GLUCOSE)%used
read(nfcell,*) iglucosegrowth
chemo(GLUCOSE)%controls_growth = (iglucosegrowth == 1)
read(nfcell,*) iglucosedeath
chemo(GLUCOSE)%controls_death = (iglucosedeath == 1)
read(nfcell,*) chemo(GLUCOSE)%diff_coef
read(nfcell,*) chemo(GLUCOSE)%medium_diff_coef
read(nfcell,*) chemo(GLUCOSE)%membrane_diff_in
read(nfcell,*) chemo(GLUCOSE)%membrane_diff_out
chemo(GLUCOSE)%membrane_diff_in = chemo(GLUCOSE)%membrane_diff_in*Vsite_cm3/60		! /min -> /sec
chemo(GLUCOSE)%membrane_diff_out = chemo(GLUCOSE)%membrane_diff_out*Vsite_cm3/60	! /min -> /sec
read(nfcell,*) chemo(GLUCOSE)%bdry_conc
read(nfcell,*) iconstant
chemo(GLUCOSE)%constant = (iconstant == 1)
read(nfcell,*) chemo(GLUCOSE)%max_cell_rate
chemo(GLUCOSE)%max_cell_rate = chemo(GLUCOSE)%max_cell_rate*1.0e6					! mol/cell/s -> mumol/cell/s
read(nfcell,*) chemo(GLUCOSE)%MM_C0
read(nfcell,*) chemo(GLUCOSE)%Hill_N

read(nfcell,*) chemo(LACTATE)%diff_coef
read(nfcell,*) chemo(LACTATE)%medium_diff_coef
read(nfcell,*) chemo(LACTATE)%membrane_diff_in
chemo(LACTATE)%membrane_diff_in = chemo(LACTATE)%membrane_diff_in*Vsite_cm3/60	! /min -> /sec
read(nfcell,*) chemo(LACTATE)%membrane_diff_out
chemo(LACTATE)%membrane_diff_out = chemo(LACTATE)%membrane_diff_out*Vsite_cm3/60	! /min -> /sec
read(nfcell,*) chemo(LACTATE)%bdry_conc
chemo(LACTATE)%bdry_conc = max(0.001,chemo(LACTATE)%bdry_conc)
read(nfcell,*) chemo(LACTATE)%max_cell_rate
chemo(LACTATE)%max_cell_rate = chemo(LACTATE)%max_cell_rate*1.0e6					! mol/cell/s -> mumol/cell/s
read(nfcell,*) chemo(LACTATE)%MM_C0
read(nfcell,*) chemo(LACTATE)%Hill_N

read(nfcell,*) iuse_tracer		!chemo(TRACER)%used
read(nfcell,*) chemo(TRACER)%diff_coef
read(nfcell,*) chemo(TRACER)%medium_diff_coef
read(nfcell,*) chemo(TRACER)%membrane_diff_in
chemo(TRACER)%membrane_diff_in = chemo(TRACER)%membrane_diff_in*Vsite_cm3/60		! /min -> /sec
chemo(TRACER)%membrane_diff_out = chemo(TRACER)%membrane_diff_in
read(nfcell,*) chemo(TRACER)%bdry_conc
read(nfcell,*) iconstant
chemo(TRACER)%constant = (iconstant == 1)
read(nfcell,*) chemo(TRACER)%max_cell_rate
read(nfcell,*) chemo(TRACER)%MM_C0
read(nfcell,*) chemo(TRACER)%Hill_N

read(nfcell,*) LQ(1)%alpha_H
read(nfcell,*) LQ(1)%beta_H
read(nfcell,*) LQ(1)%OER_am
read(nfcell,*) LQ(1)%OER_bm
read(nfcell,*) LQ(1)%K_ms
read(nfcell,*) LQ(1)%death_prob
read(nfcell,*) LQ(1)%growth_delay_factor
read(nfcell,*) LQ(1)%growth_delay_N
read(nfcell,*) LQ(2)%alpha_H
read(nfcell,*) LQ(2)%beta_H
read(nfcell,*) LQ(2)%OER_am
read(nfcell,*) LQ(2)%OER_bm
read(nfcell,*) LQ(2)%K_ms
read(nfcell,*) LQ(2)%death_prob
read(nfcell,*) LQ(2)%growth_delay_factor
read(nfcell,*) LQ(2)%growth_delay_N
read(nfcell,*) iuse_gd_all
use_radiation_growth_delay_all = (iuse_gd_all == 1)
read(nfcell,*) iusecellcycle
use_cell_cycle = (iusecellcycle == 1)
ccp => cc_parameters
call ReadCellCycleParameters(nfcell,ccp)
read(nfcell,*) iusemetabolism
use_metabolism = (iusemetabolism == 1)
call ReadMetabolismParameters(nfcell)
read(nfcell,*) O2cutoff(1)
read(nfcell,*) O2cutoff(2)
read(nfcell,*) O2cutoff(3)
read(nfcell,*) hypoxia_threshold
read(nfcell,*) growthcutoff(1)
read(nfcell,*) growthcutoff(2)
read(nfcell,*) growthcutoff(3)
read(nfcell,*) spcrad_value
read(nfcell,*) iuse_extra
read(nfcell,*) iuse_relax
read(nfcell,*) iuse_par_relax
read(nfcell,*) iuse_FD
read(nfcell,*) iuse_drop
read(nfcell,*) Ndrop
read(nfcell,*) alpha_shape
read(nfcell,*) beta_shape
read(nfcell,*) isaveprofiledata
read(nfcell,*) saveprofile%filebase
read(nfcell,*) saveprofile%dt
read(nfcell,*) saveprofile%nt
read(nfcell,*) isaveslicedata
read(nfcell,*) saveslice%filebase
read(nfcell,*) saveslice%dt
read(nfcell,*) saveslice%nt

read(nfcell,*) Ndrugs_used
call ReadDrugData(nfcell)

if (use_events) then
	call ReadProtocol(nfcell)
	use_treatment = .false.
endif

close(nfcell)

if (chemo(OXYGEN)%Hill_N /= 1 .and. chemo(OXYGEN)%Hill_N /= 2) then
	call logger('Error: OXYGEN_HILL_N must be 1 or 2')
	ok = .false.
	return
endif
if (chemo(GLUCOSE)%Hill_N /= 1 .and. chemo(GLUCOSE)%Hill_N /= 2) then
	call logger('Error: GLUCOSE_HILL_N must be 1 or 2')
	ok = .false.
	return
endif
DELTA_X = 1.0e-4*DELTA_X	! um -> cm
dxf = DELTA_X
dxb = NRF*dxf

!NX = 33		! fixed for now

NY = NX
NZ = NX
NYB = NXB
Kdrag = 1.0e5*Kdrag
MM_THRESHOLD = MM_THRESHOLD/1000					! uM -> mM
anoxia_threshold = anoxia_threshold/1000			! uM -> mM
aglucosia_threshold = aglucosia_threshold/1000		! uM -> mM
O2cutoff = O2cutoff/1000							! uM -> mM
hypoxia_threshold = hypoxia_threshold/1000			! uM -> mM
chemo(OXYGEN)%used = (iuse_oxygen == 1)
chemo(GLUCOSE)%used = (iuse_glucose == 1)
chemo(TRACER)%used = (iuse_tracer == 1)
chemo(OXYGEN)%MM_C0 = chemo(OXYGEN)%MM_C0/1000		! uM -> mM
chemo(GLUCOSE)%MM_C0 = chemo(GLUCOSE)%MM_C0/1000	! uM -> mM
chemo(LACTATE)%MM_C0 = chemo(LACTATE)%MM_C0/1000	! uM -> mM
if (.not.chemo(OXYGEN)%used) then
    chemo(OXYGEN)%controls_growth = .false.
    chemo(OXYGEN)%controls_death = .false.
endif
if (.not.chemo(GLUCOSE)%used) then
    chemo(GLUCOSE)%controls_growth = .false.
    chemo(GLUCOSE)%controls_death = .false.
endif

mitosis_duration = ccp%T_M(1)  ! seconds
LQ(:)%growth_delay_factor = 60*60*LQ(:)%growth_delay_factor	! hours -> seconds
t_anoxia_limit = 60*60*anoxia_tag_hours				! hours -> seconds
anoxia_death_delay = 60*60*anoxia_death_hours		! hours -> seconds
t_aglucosia_limit = 60*60*aglucosia_tag_hours		! hours -> seconds
aglucosia_death_delay = 60*60*aglucosia_death_hours	! hours -> seconds
nsteps = days*24*3600./DELTA_T
write(logmsg,*) 'nsteps: ',nsteps
call logger(logmsg)

divide_dist(1:2)%class = LOGNORMAL_DIST
divide_time_median(1:2) = 60*60*divide_time_median(1:2)		! hours -> seconds
sigma(1:2) = log(divide_time_shape(1:2))
!divide_dist%p1 = log(divide_time_mean/exp(sigma*sigma/2))	
divide_dist(1:2)%p1 = log(divide_time_median(1:2))	
divide_dist(1:2)%p2 = sigma(1:2)
divide_time_mean(1:2) = exp(divide_dist(1:2)%p1 + 0.5*divide_dist(1:2)%p2**2)	! mean = median.exp(sigma^2/2)
write(logmsg,'(a,24e12.4)') 'shape, sigma: ',divide_time_shape(1:2),sigma(1:2)
call logger(logmsg)
write(logmsg,'(a,4e12.4)') 'Median, mean divide time: ',divide_time_median(1:2)/3600,divide_time_mean(1:2)/3600
call logger(logmsg)

use_V_dependence = (iV_depend == 1)
saveprofile%active = (isaveprofiledata == 1)
saveprofile%it = 1
saveprofile%dt = 60*saveprofile%dt		! mins -> seconds
saveslice%active = (isaveslicedata == 1)
saveslice%it = 1
saveslice%dt = 60*saveslice%dt			! mins -> seconds

use_dropper = (iuse_drop == 1)

randomise_initial_volume = (iV_random == 1)
use_extracellular_O2 = (iuse_extra == 1)

open(nfout,file=outputfile,status='replace')
write(nfout,'(a,a)') 'GUI version: ',gui_run_version
write(nfout,'(a,a)') 'DLL version: ',dll_run_version
write(nfout,*)

write(nflog,*)
write(nflog,'(a,a)') 'GUI version: ',gui_run_version
write(nflog,'(a,a)') 'DLL version: ',dll_run_version
write(nflog,*)

open(nfres,file='vspheroid_ts.out',status='replace')
write(nfres,'(a)') 'data info GUI_version DLL_version &
istep hour vol_mm3 diam_um Ncells(1) Ncells(2) &
Nanoxia_dead(1) Nanoxia_dead(2) Naglucosia_dead(1) Naglucosia_dead(2) NdrugA_dead(1) NdrugA_dead(2) &
NdrugB_dead(1) NdrugB_dead(2) Nradiation_dead(1) Nradiation_dead(2) &
Ntagged_anoxia(1) Ntagged_anoxia(2) Ntagged_aglucosia(1) Ntagged_aglucosia(2) Ntagged_drugA(1) Ntagged_drugA(2) &
Ntagged_drugB(1) Ntagged_drugB(2) Ntagged_radiation(1) Ntagged_radiation(2) &
f_hypox(1) f_hypox(2) f_hypox(3) &
f_clonohypox(1) f_clonohypox(2) f_clonohypox(3) &
f_growth(1) f_growth(2) f_growth(3) &
f_necrot plating_efficiency(1) plating_efficiency(2) &
medium_oxygen medium_glucose medium_lactate &
medium_drugA medium_drugA_metab1 medium_drugA_metab2 medium_drugB medium_drugB_metab1 medium_drugB_metab2 &
bdry_oxygen bdry_glucose bdry_lactate &
bdry_drugA bdry_drugA_metab1 bdry_drugA_metab2 bdry_drugB bdry_drugB_metab1 bdry_drugB_metab2 &
doubling_time glycolysis_rate pyruvate_oxidation_rate ATP_rate intermediates_rate Ndivided pyruvate_oxidised_fraction'

write(logmsg,*) 'Opened nfout: ',outputfile
call logger(logmsg)
call DetermineKd
ok = .true.

end subroutine

!-----------------------------------------------------------------------------------------
! The cell cycle parameters include the parameters for radiation damage and repair, 
! and for the associated checkpoint duration limits Tcp(:).
! Time unit = hour
!-----------------------------------------------------------------------------------------
subroutine ReadCellCycleParameters(nf, ccp)
integer :: nf
type(cycle_parameters_type),pointer :: ccp

read(nf,*) ccp%T_G1(1)
read(nf,*) ccp%T_G1(2)
read(nf,*) ccp%T_S(1)
read(nf,*) ccp%T_S(2)
read(nf,*) ccp%T_G2(1)
read(nf,*) ccp%T_G2(2)
read(nf,*) ccp%T_M(1)
read(nf,*) ccp%T_M(2)
read(nf,*) ccp%G1_mean_delay(1)
read(nf,*) ccp%G1_mean_delay(2)
read(nf,*) ccp%G2_mean_delay(1)
read(nf,*) ccp%G2_mean_delay(2)
read(nf,*) ccp%eta_PL
read(nf,*) ccp%eta_L(1)
read(nf,*) ccp%eta_L(2)
read(nf,*) ccp%Krepair_base
read(nf,*) ccp%Krepair_max
read(nf,*) ccp%Kmisrepair(1)
read(nf,*) ccp%Kmisrepair(2)
read(nf,*) ccp%Kcp

ccp%T_G1 = 3600*ccp%T_G1                    ! hours -> seconds
ccp%T_S = 3600*ccp%T_S
ccp%T_G2 = 3600*ccp%T_G2
ccp%T_M = 3600*ccp%T_M
ccp%G1_mean_delay = 3600*ccp%G1_mean_delay
ccp%G2_mean_delay = 3600*ccp%G2_mean_delay

ccp%Pk_G1 = 1./ccp%G1_mean_delay    ! /sec
ccp%Pk_G2 = 1./ccp%G2_mean_delay    ! /sec
call logger('makeTCP')
call makeTCP(ccp%tcp,NTCP,ccp%Krepair_base,ccp%Kmisrepair,ccp%Kcp) ! set checkpoint repair time limits 
call logger('did makeTCP')
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine ReadMetabolismParameters(nf)
integer :: nf
integer :: ityp

do ityp = 1,Ncelltypes
	read(nf,*) N_GA(ityp)
	read(nf,*) N_PA(ityp)
	read(nf,*) N_GI(ityp)
	read(nf,*) N_PI(ityp)
	read(nf,*) N_PO(ityp)
	read(nf,*) K_H1(ityp)
	read(nf,*) K_H2(ityp)
	read(nf,*) K_HB(ityp)
	read(nf,*) K_PDK(ityp)
	read(nf,*) PDKmin(ityp)
	read(nf,*) C_O2_norm(ityp)
	read(nf,*) C_G_norm(ityp)
	read(nf,*) C_L_norm(ityp)
!	read(nf,*) f_ATPg(ityp)
	read(nf,*) f_ATPs(ityp)
	read(nf,*) K_PL(ityp)
	read(nf,*) K_LP(ityp)
	read(nf,*) Hill_Km_P(ityp)
	read(nf,*) Apoptosis_rate(ityp)
enddo
PDKmin(:) = 0.3
Hill_N_P = 1
Hill_Km_P = Hill_Km_P/1000	! uM -> mM
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine ReadDrugData(nf)
integer :: nf
integer :: idrug, im, ictyp, ival
character*(16) :: drugname

if (allocated(drug)) then
	deallocate(drug)
endif
allocate(drug(Ndrugs_used))
do idrug = 1,Ndrugs_used
	read(nf,'(a)') drug(idrug)%classname
	if (drug(idrug)%classname == 'TPZ') then
		drug(idrug)%drugclass = TPZ_CLASS
	elseif (drug(idrug)%classname == 'DNB') then
		drug(idrug)%drugclass = DNB_CLASS
	endif
	drug(idrug)%use_metabolites = .true.	! currently by default all drugs use 2 metabolites
	drug(idrug)%nmetabolites = 2		
    do im = 0,2			! 0 = parent, 1 = metab_1, 2 = metab_2 
		read(nf,'(a)') drugname
		if (im == 0) then
			drug(idrug)%name = drugname
		endif
		read(nf,*) drug(idrug)%diff_coef(im)
		read(nf,*) drug(idrug)%medium_diff_coef(im)
		read(nf,*) drug(idrug)%membrane_diff_in(im)
		read(nf,*) drug(idrug)%membrane_diff_out(im)
		read(nf,*) drug(idrug)%halflife(im)
		drug(idrug)%membrane_diff_in(im) = drug(idrug)%membrane_diff_in(im)*Vsite_cm3/60	! /min -> /sec
		drug(idrug)%membrane_diff_out(im) = drug(idrug)%membrane_diff_out(im)*Vsite_cm3/60	! /min -> /sec
		do ictyp = 1,ncelltypes
            read(nf,*) drug(idrug)%Kmet0(ictyp,im)
            read(nf,*) drug(idrug)%C2(ictyp,im)
            read(nf,*) drug(idrug)%KO2(ictyp,im)
            read(nf,*) drug(idrug)%Vmax(ictyp,im)
            read(nf,*) drug(idrug)%Km(ictyp,im)
            read(nf,*) drug(idrug)%Klesion(ictyp,im)
            read(nf,*) drug(idrug)%kill_O2(ictyp,im)
            read(nf,*) drug(idrug)%kill_drug(ictyp,im)
            read(nf,*) drug(idrug)%kill_duration(ictyp,im)
            read(nf,*) drug(idrug)%kill_fraction(ictyp,im)
            read(nf,*) drug(idrug)%SER_max(ictyp,im)
            read(nf,*) drug(idrug)%SER_Km(ictyp,im)
            read(nf,*) drug(idrug)%SER_KO2(ictyp,im)
            read(nf,*) drug(idrug)%n_O2(ictyp,im)
            read(nf,*) drug(idrug)%death_prob(ictyp,im)
            read(nf,*) ival
            drug(idrug)%kills(ictyp,im) = (ival == 1)
            read(nf,*) ival
            drug(idrug)%kill_model(ictyp,im) = ival
            read(nf,*) ival
            drug(idrug)%sensitises(ictyp,im) = (ival == 1)
            drug(idrug)%Kmet0(ictyp,im) = drug(idrug)%Kmet0(ictyp,im)/60					! /min -> /sec
            drug(idrug)%KO2(ictyp,im) = 1.0e-3*drug(idrug)%KO2(ictyp,im)					! um -> mM
            drug(idrug)%kill_duration(ictyp,im) = 60*drug(idrug)%kill_duration(ictyp,im)	! min -> sec
		enddo
    enddo
    write(nflog,*) 'drug: ',idrug,drug(idrug)%classname,'  ',drug(idrug)%name
enddo
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
!subroutine getIndices(indx, idrug, nmetab)
!integer :: indx, idrug, nmetab
!
!if (indx == 1) then
!	idrug = TPZ_DRUG
!	nmetab = 2
!elseif (indx == 2) then
!	idrug = DNB_DRUG
!	nmetab = 2
!else
!	idrug = 0
!	nmetab = 0
!endif
!end subroutine

!-----------------------------------------------------------------------------------------
! Skip lines until the 'PROTOCOL' line
!-----------------------------------------------------------------------------------------
subroutine ReadProtocol(nf)
integer :: nf
integer :: itime, ntimes, kevent, ichemo, idrug, im
character*(64) :: line
character*(16) :: drugname
character*(1)  :: numstr
real(REAL_KIND) :: t, dt, vol, conc, O2conc, O2flush, dose, O2medium
type(event_type) :: E

write(logmsg,*) 'ReadProtocol'
call logger(logmsg)
chemo(TRACER+1:)%used = .false.
do
	read(nf,'(a)') line
	if (trim(line) == 'PROTOCOL') exit
enddo
read(nf,*) ntimes
if (ntimes == 0) then
	Nevents = 0
	return
endif
if (allocated(event)) deallocate(event)
allocate(event(2*ntimes))
kevent = 0
do itime = 1,ntimes
	read(nf,'(a)') line
	if (trim(line) == 'DRUG') then
		kevent = kevent + 1
		event(kevent)%etype = DRUG_EVENT
		read(nf,'(a)') line
		drugname = trim(line)
		do idrug = 1,ndrugs_used
			if (drugname == drug(idrug)%name) then
				ichemo = DRUG_A + 3*(idrug-1)
				exit
			endif
		enddo
		! Need to copy drug(idrug) parameters to chemo(ichemo) 
		call CopyDrugParameters(idrug,ichemo)
		read(nf,*) t
		read(nf,*) dt
		read(nf,*) vol
		read(nf,*) O2conc
		read(nf,*) O2flush
		read(nf,*) conc
		event(kevent)%time = t
		event(kevent)%ichemo = ichemo
		event(kevent)%idrug = idrug
		event(kevent)%volume = vol
		event(kevent)%conc = conc
		event(kevent)%O2conc = O2conc
		event(kevent)%dose = 0
		chemo(ichemo)%used = .true.
		write(nflog,'(a,i3,2f8.3)') 'define DRUG_EVENT: volume, O2conc: ',kevent,event(kevent)%volume,event(kevent)%O2conc
		if (drug(idrug)%use_metabolites) then
			do im = 1,drug(idrug)%nmetabolites
				chemo(ichemo+im)%used = .true.
			enddo
		endif

		kevent = kevent + 1
		event(kevent)%etype = MEDIUM_EVENT
		event(kevent)%time = t + dt
		event(kevent)%ichemo = 0
		event(kevent)%volume = medium_volume0
		event(kevent)%conc = 0
		event(kevent)%O2medium = O2flush
		event(kevent)%dose = 0
		write(nflog,'(a,i3,2f8.3)') 'define MEDIUM_EVENT: volume: ',kevent,event(kevent)%volume,event(kevent)%O2medium
	elseif (trim(line) == 'MEDIUM') then
		kevent = kevent + 1
		event(kevent)%etype = MEDIUM_EVENT
		read(nf,*) t
		read(nf,*) vol
		read(nf,*) O2medium
		event(kevent)%time = t
		event(kevent)%volume = vol	
		event(kevent)%ichemo = 0
		event(kevent)%O2medium = O2medium
		event(kevent)%dose = 0
		write(nflog,'(a,i3,2f8.3)') 'define MEDIUM_EVENT: volume: ',kevent,event(kevent)%volume,event(kevent)%O2medium
	elseif (trim(line) == 'RADIATION') then
		kevent = kevent + 1
		event(kevent)%etype = RADIATION_EVENT
		read(nf,*) t
		read(nf,*) dose
		event(kevent)%time = t
		event(kevent)%dose = dose	
		event(kevent)%ichemo = 0
		event(kevent)%volume = 0
		event(kevent)%conc = 0
	endif
enddo
Nevents = kevent
! Set events not done
! convert time from hours to seconds
! convert volume from uL to cm^3  NO LONGER - now using cm^3 everywhere
do kevent = 1,Nevents
	event(kevent)%done = .false.
	event(kevent)%time = event(kevent)%time*60*60
!	event(kevent)%V = event(kevent)%V*1.0e-3
	E = event(kevent)
!	write(*,'(a,i3,f8.0,2i3,3f8.4)') 'event: ',kevent,E%time,E%etype,E%ichemo,E%V,E%conc,E%dose
enddo
! Check that events are sequential
do kevent = 1,Nevents-1
	if (event(kevent)%time >= event(kevent+1)%time) then
		write(logmsg,*) 'Error: non-sequential event: ',kevent,event(kevent)%time
		call logger(logmsg)
		stop
	endif
enddo
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine CopyDrugParameters(idrug,ichemo)
integer :: idrug,ichemo
integer :: im, im1, im2
character*(1) :: numstr

im1 = 0
chemo(ichemo)%name = drug(idrug)%name
if (drug(idrug)%use_metabolites) then
	do im = 1,drug(idrug)%nmetabolites
		chemo(ichemo+im)%used = .true.
		chemo(ichemo+im)%name = trim(chemo(ichemo)%name) // '_metab'
		write(numstr,'(i1)') im
		chemo(ichemo+im)%name = trim(chemo(ichemo+im)%name) // numstr
	enddo
	im2 = 2
else
	im2 = 0
endif
do im = im1, im2
	chemo(ichemo+im)%diff_coef = drug(idrug)%diff_coef(im)
	chemo(ichemo+im)%medium_diff_coef = drug(idrug)%medium_diff_coef(im)
	chemo(ichemo+im)%membrane_diff_in = drug(idrug)%membrane_diff_in(im)
	chemo(ichemo+im)%membrane_diff_out = drug(idrug)%membrane_diff_out(im)
	chemo(ichemo+im)%halflife = drug(idrug)%halflife(im)
!	chemo(ichemo+im)%medium_dlayer = d_layer
	chemo(ichemo+im)%decay = (chemo(ichemo+im)%halflife > 0)
	if (chemo(ichemo+im)%decay) then
		chemo(ichemo+im)%decay_rate = DecayRate(chemo(ichemo+im)%halflife)
	else
		chemo(ichemo+im)%decay_rate = 0
	endif
enddo
end subroutine

!-----------------------------------------------------------------------------------------
! d~ = d - (R1+R2)
! d_detach is the value of d~ at which V -> 0, i.e. it depends on R1+R2
!-----------------------------------------------------------------------------------------
subroutine setup(ncpu, infile, outfile, ok)
integer :: ncpu
character*(128) :: infile, outfile
logical :: ok
integer :: kcell
real(REAL_KIND) :: k_v, tgrowth(MAX_CELLTYPES)
type(cycle_parameters_type),pointer :: ccp
ok = .true.
initialized = .false.
par_zig_init = .false.

inputfile = infile
outputfile = outfile
call logger("ReadCellParams")
call ReadCellParams(ok)
if (.not.ok) return
call logger("did ReadCellParams")

if (ncpu == 0) then
	ncpu = ncpu_input
endif
Mnodes = ncpu
write(logmsg,*) 'ncpu: ',ncpu 
call logger(logmsg)

#if defined(OPENMP) || defined(_OPENMP)
    call logger("OPENMP defined")
    call omp_initialisation(ok)
    if (.not.ok) return
#else
    call logger("OPENMP NOT defined")
    if (Mnodes > 1) then
        write(logmsg,'(a)') 'No OpenMP, using one thread only'
        call logger(logmsg)
        Mnodes = 1
    endif
#endif

NY = NX
NZ = NX
Ndim(1) = NX
Ndim(2) = NY
Ndim(3) = NZ
call ArrayInitialisation(ok)
if (.not.ok) return
call logger('did ArrayInitialisation')

!call create_stencil

chemo(LACTATE)%used = use_metabolism
call SetupChemo
call logger('did SetupChemo')
if (use_metabolism) then
	call SetupMetabolism
endif

! Assume that Raverage is for cells midway in growth, i.e. between 0.8 and 1.6, at 1.2
! as a fraction of the average cell volume, (or, equivalently, between 1.0 and 2.0, at 1.5)
! and that actual divide volume Vdivide is: Vdivide0-dVdivide < Vdivide < Vdivide0+dVdivide
! Note that 2.0/1.5 = 1.6/1.2
! Growth occurs at a constant rate for (divide_time - mitosis_duration)
! This is option A.
! Try a different option: B
Vdivide0 = (2.0/1.5)*(4.*PI/3.)*Raverage**3
dVdivide = dVdivide*Vdivide0	
Rdivide0 = Raverage*(2.0/1.5)**(1./3.)
d_nbr_limit = 1.2*2*Rdivide0	! 1.5 is an arbitrary choice - was 1.2

!test_growthrate = Vdivide0/(2*(divide_time_median(1) - mitosis_duration))	! um3/sec
write(logmsg,'(a,2e12.3)') 'Raverage,Vdivide0: ',Raverage,Vdivide0
call logger(logmsg)

! New cell cycle formulation - need a value for max (unconstrained) growth rate
use_volume_method = .not.use_cell_cycle
if (use_cell_cycle .and. .not.use_volume_based_transition) then
    use_constant_growthrate = .true.
endif
! Growth occurs during G1, S and G2, not in checkpoints
ccp => cc_parameters
tgrowth = ccp%T_G1 + ccp%T_S + ccp%T_G2
max_growthrate = Vdivide0/(2*tgrowth)

call setup_force_parameters

call make_jumpvec
call PlaceCells(ok)
call logger('did PlaceCells')
!do kcell = 1,nlist
!	write(nflog,'(i6,3f8.1)') kcell,1.0e4*cp%centre(:,1)
!enddo
alpha_v = 0.25
epsilon = 7.5
es_e = 1
shift = -6
Dfactor = 1
sqr_es_e = sqrt(es_e)

ntries = 1
t_fmover = 0
ndtotal = 0
ndtotal_last = 0
total_dMdt = 0

Nradiation_tag = 0
Ndrug_tag = 0
Nanoxia_tag = 0
Naglucosia_tag = 0
Nradiation_dead = 0
Ndrug_dead = 0
Nanoxia_dead = 0
Naglucosia_dead = 0

ndoublings = 0
doubling_time_sum = 0
ncells_mphase = 0

is_dropped = .false.

k_v = 2/alpha_v - sqr_es_e + sqrt(es_e - shift/epsilon)
k_detach = k_v*alpha_v/2
!write(*,'(a,f8.4)') 'k_detach: ',k_detach

call setup_nbrlists(ok)
if (.not.ok) return
call logger('did setup_nbrlists')

call setup_react_diff(ok)
! For simple testing...
!cell_list(:)%Cin(OXYGEN) = chemo(OXYGEN)%bdry_conc
!cell_list(:)%Cin(GLUCOSE) = chemo(GLUCOSE)%bdry_conc
call logger('did setup_react_diff')

call SetInitialGrowthRate
limit_stop = .false.
medium_change_step = .false.		! for startup 
call logger('completed Setup')

nspeedtest = 10000
!call testOGL
!stop

end subroutine

!-----------------------------------------------------------------------------------------
! SN30000 CT model
! ----------------
! The rate of cell killing, which becomes a cell death probability rate, is inferred from
! the cell kill experiment.
! The basic assumption is that the rate of killing depends on the drug metabolism rate.
! There are five models:
! kill_model = 1:
!   killing rate = c = Kd.dM/dt
! kill_model = 2:
!   killing rate = c = Kd.Ci.dM/dt
! kill_model = 3:
!   killing rate = c = Kd.(dM/dt)^2
! kill_model = 4:
!   killing rate = c = Kd.Ci
! kill_model = 5:
!   killing rate = c = Kd.Ci^2
! where dM/dt = F(O2).kmet0.Ci
! In the kill experiment both O2 and Ci are held constant:
! O2 = CkillO2, Ci = Ckill
! In this case c is constant and the cell population N(t) is given by:
! N(t) = N(0).exp(-ct), i.e.
! c = -log(N(T)/N(0))/T where T is the duration of the experiment
! N(T)/N(0) = 1 - f, where f = kill fraction
! kill_model = 1:
!   c = Kd.F(CkillO2).kmet0.Ckill => Kd = -log(1-f)/(T.F(CkillO2).kmet0.Ckill)
! kill_model = 2:
!   c = Kd.F(CkillO2).kmet0.Ckill^2 => Kd = -log(1-f)/(T.F(CkillO2).kmet0.Ckill^2)
! kill_model = 3:
!   c = Kd.(F(CkillO2).kmet0.Ckill)^2 => Kd = -log(1-f)/(T.(F(CkillO2).kmet0.Ckill)^2)
! kill_model = 4:
!   c = Kd.Ckill => Kd = -log(1-f)/(T.Ckill)
! kill_model = 5:
!   c = Kd.Ckill^2 => Kd = -log(1-f)/(T.Ckill^2)
!-----------------------------------------------------------------------------------------
subroutine DetermineKd
real(REAL_KIND) :: C2, KO2, n_O2, Kmet0, kmet 
real(REAL_KIND) :: f, T, Ckill, Ckill_O2, Kd
integer :: idrug, ictyp, im, kill_model

do idrug = 1,ndrugs_used
!	if (idrug == 1 .and. .not.chemo(TPZ_DRUG)%used) cycle
!	if (idrug == 2 .and. .not.chemo(DNB_DRUG)%used) cycle
	do ictyp = 1,Ncelltypes
		do im = 0,2
			if (drug(idrug)%kills(ictyp,im)) then
				C2 = drug(idrug)%C2(ictyp,im)
				KO2 = drug(idrug)%KO2(ictyp,im)
				n_O2 = drug(idrug)%n_O2(ictyp,im)
				Kmet0 = drug(idrug)%Kmet0(ictyp,im)
				kill_model = drug(idrug)%kill_model(ictyp,im)
				Ckill_O2 = drug(idrug)%kill_O2(ictyp,im)
				f = drug(idrug)%kill_fraction(ictyp,im)
				T = drug(idrug)%kill_duration(ictyp,im)
				Ckill = drug(idrug)%kill_drug(ictyp,im)
				kmet = (1 - C2 + C2*(KO2**n_O2)/(KO2**n_O2 + Ckill_O2**n_O2))*Kmet0
				if (kill_model == 1) then
					Kd = -log(1-f)/(T*kmet*Ckill)
				elseif (kill_model == 2) then
					Kd = -log(1-f)/(T*kmet*Ckill**2)
				elseif (kill_model == 3) then
					Kd = -log(1-f)/(T*(kmet*Ckill)**2)
				elseif (kill_model == 4) then
					Kd = -log(1-f)/(T*Ckill)
				elseif (kill_model == 5) then
					Kd = -log(1-f)/(T*Ckill**2)
				endif
				drug(idrug)%Kd(ictyp,im) = Kd
			endif
!			if (idrug == 1) then
!				TPZ%Kd(i) = Kd
!			elseif (idrug == 2) then
!				DNB%Kd(i,im) = Kd
!			endif
		enddo
	enddo
enddo

end subroutine


!-----------------------------------------------------------------------------------------
! Cells are initially placed in a regular rectangular grid pattern, spacing between cell 
! centres equal to 2*Raverage
!-----------------------------------------------------------------------------------------
subroutine PlaceCells(ok)
logical :: ok
integer :: ix, iy, iz, kcell, i, site(3), irad, lastid, nblob
real(REAL_KIND) :: Radius, cellradius, d, r2lim, r2, rad(3), rsite(3), centre(3), ave(3)
logical, allocatable :: occup(:,:,:)

!blobcentre = DELTA_X*[(NX+1)/2,(NY+1)/2,(NZ+1)/2]
blobcentre = DELTA_X*[(NX-1)/2,(NY-1)/2,(NZ-1)/2]
write(nflog,'(a,3f8.4)') 'PlaceCells: blobcentre: ',blobcentre
if (use_packer) then
	nblob = initial_count
	cellradius = Raverage
	call SelectCellLocations(nblob, cellradius)
	ave = 0
	do i = 1,nblob
		call get_cellcentre(i,centre)
		ave = ave + centre
	enddo
	ave = ave/nblob
	write(nflog,'(a,3f8.4)') 'Average cell centre: ',ave
	write(nflog,'(a,3f8.4)') 'blobcentre: ',blobcentre
!	write(*,*) 'Cell centres:'
	kcell = 0
	do i = 1,nblob
		call get_cellcentre(i,centre)
		centre = centre - ave + blobcentre
		kcell = kcell + 1
		rsite = centre
!		write(*,'(i6,3f8.4)') kcell,rsite
		call AddCell(kcell,rsite)
	enddo
	call FreeCellLocations
else
	d = 2.2*Raverage
	Radius = (3.0*initial_count/(4.0*PI))**(1./3.)	! approx initial blob radius, scaled by /d
	write(nflog,*) 'blobcentre, d: ',blobcentre,d,Radius
	irad = Radius + 2
	allocate(occup(-irad:irad,-irad:irad,-irad:irad))
	occup = .false.
	r2lim = 0.97*Radius*Radius
	lastID = 0
!	write(*,'(a,3f8.4)') 'blobcentre: ',blobcentre
!	write(*,*) 'Cell centres:'
	kcell = 0
	do ix = -irad, irad
		do iy = -irad, irad
			do iz = -irad, irad
				site = (/ix,iy,iz/)
				rad = site
				r2 = dot_product(rad,rad)
				if (r2 > r2lim) cycle
				kcell = kcell + 1
				rsite = blobcentre + d*site
				call AddCell(kcell,rsite)
				occup(ix,iy,iz) = .true.
			enddo
		enddo
	enddo
	if (kcell > initial_count) then
		write(logmsg,*) 'Cell count already exceeds specified number: ',kcell,initial_count
		call logger(logmsg)
		ok = .false.
		return
	endif
	! Now add cells to make the count up to the specified initial_count
	if (kcell < initial_count) then
		call AddBdryCells(kcell, occup, irad, d, blobcentre)
		kcell = initial_count
	endif
	deallocate(occup)
endif
nlist = kcell
ncells = kcell
end subroutine

!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------
subroutine AddCell(kcell, rsite)
integer :: kcell
real(REAL_KIND) :: rsite(3)
integer :: ityp, k, kpar = 0
real(REAL_KIND) :: v(3), c(3), R1, R2, V0, Tdiv, Vdiv, p(3), R
type(cell_type), pointer :: cp
type(cycle_parameters_type),pointer :: ccp
	
cp => cell_list(kcell)
ccp => cc_parameters
cp%ID = kcell
cp%state = ALIVE
cp%generation = 1
cp%celltype = random_choice(celltype_fraction,Ncelltypes,kpar)
ityp = cp%celltype
if (ityp == 2) then
	write(logmsg,*) 'So far only cell type 1 is implemented'
	call logger(logmsg)
	stop
endif
Ncells_type(ityp) = Ncells_type(ityp) + 1
cp%Iphase = .true.
cp%nspheres = 1

V0 = Vdivide0/2
cp%divide_volume = get_divide_volume(ityp, V0, Tdiv)
cp%divide_time = Tdiv
cp%dVdt = max_growthrate(ityp)
if (use_volume_method) then
    !cp%divide_volume = Vdivide0
    if (initial_count == 1) then
	    cp%V = 0.9*cp%divide_volume
    else
    !	cp%V = (0.5 + 0.49*par_uni(kpar))*cp%divide_volume
        if (randomise_initial_volume) then
	        cp%V = cp%divide_volume*0.5*(1 + par_uni(kpar))
        else
	        cp%V = cp%divide_volume/1.6
        endif
    endif
    cp%t_divide_last = 0    ! not correct
else
    cp%NL1 = 0
    cp%NL2 = 0
    ! Need to assign phase, volume to complete phase, current volume
    ! Simplest way to assign phase fractions is in proportion to phase durations
    ! (exclude M-phase)
    p(1) = ccp%T_G1(ityp)
    p(2) = ccp%T_S(ityp)
    p(3) = ccp%T_G2(ityp)
    p = p/sum(p)
    k = random_choice(p,3,kpar)
    if (k == 1) then
        cp%phase = G1_phase
        cp%G1_flag = .false.
        R = par_uni(kpar)
        cp%G1_time = R*ccp%T_G1(ityp)
        cp%V = V0 + max_growthrate(ityp)*R*ccp%T_G1(ityp)
        cp%G1_V = V0 + max_growthrate(ityp)*ccp%T_G1(ityp)
        cp%t_divide_last = cp%G1_time - ccp%T_G1(ityp)
    elseif (k == 2) then
        cp%phase = S_phase
        R = par_uni(kpar)
        cp%S_time = R*ccp%T_S    (ityp)
        cp%V = V0 + max_growthrate(ityp)*(ccp%T_G1(ityp) + R*ccp%T_S(ityp))
        cp%S_V = V0 + max_growthrate(ityp)*(ccp%T_G1(ityp) + ccp%T_S(ityp))
        cp%t_divide_last = cp%S_time - ccp%T_S(ityp) - ccp%T_G1(ityp) - ccp%G1_mean_delay(ityp)
    elseif (k == 3) then
        cp%phase = G2_phase
        R = par_uni(kpar)
        cp%G2_time = R*ccp%T_G2(ityp)
        cp%V = V0 + max_growthrate(ityp)*(ccp%T_G1(ityp) + ccp%T_S(ityp) + R*ccp%T_G2(ityp))
        cp%S_V = V0 + max_growthrate(ityp)*(ccp%T_G1(ityp) + ccp%T_S(ityp) + ccp%T_G2(ityp))
        cp%t_divide_last = cp%G2_time - ccp%T_G2(ityp) - ccp%T_S(ityp) - ccp%T_G1(ityp) - ccp%G1_mean_delay(ityp)
    endif
endif
cp%radius(1) = (3*cp%V/(4*PI))**(1./3.)
cp%centre(:,1) = rsite
cp%site = rsite/DELTA_X + 1
cp%d = 0
cp%birthtime = 0
!cp%growthrate = test_growthrate
!cp2%divide_volume = get_divide_volume()
cp%d_divide = (3*cp%divide_volume/PI)**(1./3.)
cp%mitosis = 0

cp%drug_tag = .false.
cp%radiation_tag = .false.
cp%anoxia_tag = .false.
cp%aglucosia_tag = .false.
cp%growth_delay = .false.
cp%G2_M = .false.
cp%p_rad_death = 0

cp%t_anoxia = 0
cp%t_aglucosia = 0

call get_random_vector3(v)	! set initial axis direction
cp%d = 0.1*small_d
c = cp%centre(:,1)
cp%centre(:,1) = c + (cp%d/2)*v
cp%centre(:,2) = c - (cp%d/2)*v
cp%nbrs = 0
cp%Cex = Caverage(1,1,1,:)	! initially the concentrations are the same everywhere
cp%Cin = cp%Cex
cp%CFSE = generate_CFSE(1.d0)

cp%ndt = ndt

if (use_metabolism) then
	cp%metab = metabolic(1)
	cp%growth_rate_factor = get_growth_rate_factor()
	cp%ATP_rate_factor = get_ATP_rate_factor()
	if (cp%metab%A_rate == 0) then
		write(*,*) 'A_rate = 0'
		stop
	endif
endif
end subroutine

!--------------------------------------------------------------------------------
! Add cells at the boundary to bring the total count from k up to initial_count
! (1) Make a list of all boundary sites (sites in contact with an OUTSIDE site)
! (2) Iteratively traverse the list to select the adjacent OUTSIDE site closest 
! to the centre.
!--------------------------------------------------------------------------------
subroutine AddBdryCells(klast, occup, irad, d, centre)
integer :: klast, irad
logical :: occup(-irad:irad,-irad:irad,-irad:irad)
real(REAL_KIND) :: d, centre(3)
integer :: ix, iy, iz
integer :: kcell, i, kb, site(3), nbsite(3), nbt, kbmin, imin
integer, allocatable :: sitelist(:,:)
real(REAL_KIND) :: r2, r2min, rad(3), rsite(3)

nbt = 0
do ix = -irad, irad
	do iy = -irad, irad
		do iz = -irad, irad
			if (.not.occup(ix,iy,iz)) cycle
			site = [ix,iy,iz]
			do i = 1,27
				if (i == 14) cycle
				nbsite = site + jumpvec(:,i)
				if (.not.occup(nbsite(1),nbsite(2),nbsite(3))) then
					nbt = nbt+1
					exit
				endif
			enddo
		enddo
	enddo
enddo

allocate(sitelist(3,nbt))

nbt = 0
do ix = -irad, irad
	do iy = -irad, irad
		do iz = -irad, irad
			if (.not.occup(ix,iy,iz)) cycle
			site = [ix,iy,iz]
			do i = 1,27
				if (i == 14) cycle
				nbsite = site + jumpvec(:,i)
				if (.not.occup(nbsite(1),nbsite(2),nbsite(3))) then
					nbt = nbt+1
					sitelist(:,nbt) = site
					exit
				endif
			enddo
		enddo
	enddo
enddo	

do kcell = klast+1,initial_count
	r2min = 1.0e10
	do kb = 1,nbt
		site = sitelist(:,kb)
		do i = 1,27
			if (i == 14) cycle
			nbsite = site + jumpvec(:,i)
			if (.not.occup(nbsite(1),nbsite(2),nbsite(3))) then
				rad = d*site
				r2 = dot_product(rad,rad)
				if (r2 < r2min) then
					kbmin = kb
					imin = i
					r2min = r2
				endif
			endif
		enddo
	enddo
	site = sitelist(:,kbmin) + jumpvec(:,imin)
	rsite = centre + d*site
	call AddCell(kcell,rsite)
	occup(site(1),site(2),site(3)) = .true.
enddo

deallocate(sitelist)
		
end subroutine

!-----------------------------------------------------------------------------------------
! Simulate through a full time step: DELTA_T 
!-----------------------------------------------------------------------------------------
subroutine simulate_step(res) BIND(C)
!DEC$ ATTRIBUTES DLLEXPORT :: simulate_step  
use, intrinsic :: iso_c_binding
integer(c_int) :: res
integer :: kcell, hour, kpar=0
real(REAL_KIND) :: radiation_dose, dt, volume, maxarea, t_ave_double
integer :: i, k, nit, irepeat, nrepeat, nt_diff, it_diff, ncells0, nhypoxic(3)
integer :: nshow = 100
integer :: Nhop, nt_hour, nt_nbr
integer :: nvars, ns
real(REAL_KIND) :: dxc, ex_conc(35*O2_BY_VOL+1)		! just for testing
real(REAL_KIND) :: t0, t1, tmover, tgrower, ave_cin(0:2), ave_cex(0:2)
real(REAL_KIND), parameter :: FULL_NBRLIST_UPDATE_HOURS = 4
type(cell_type), pointer :: cp
logical :: ok, done, changed

if (kcell_test > 0) then
	cp => cell_list(kcell_test)
else
	cp => cell_list(1)
endif
!Nhop = 10*(30/DELTA_T)
Nhop = 1
nt_hour = 3600/DELTA_T
nt_nbr = nt_hour*FULL_NBRLIST_UPDATE_HOURS
tmover = 0
tgrower = 0
if (Ncells == 0) then
	call logger('Ncells = 0')
    res = 2
    return
endif
if (limit_stop) then
	call logger('Spheroid size limit reached')
	res = 6
	return
endif

if (istep == -100) then
	tnow = istep*DELTA_T
!	call make_colony_distribution(tnow)
!	stop
endif

if (use_dropper .and. .not.is_dropped .and. ncells > ndrop) then
	call drop_blob
	is_dropped = .true.
endif

!call getVolume(volume,maxarea)
blobradius = getRadius()
blobcentre = getCentre()

istep = istep + 1
t_simulation = (istep-1)*DELTA_T	! seconds
radiation_dose = 0
if (use_events) then
	call ProcessEvent(radiation_dose)
endif
if (radiation_dose > 0) then
	write(logmsg,'(a,f6.1)') 'Radiation dose: ',radiation_dose
	call logger(logmsg)
endif
call SetupChemomap
!dt = DELTA_T/NT_CONC
! the idea is to accumulate time steps until DELTA_T is reached 
call make_perm_index(ok)
if (.not.ok) then
	call logger('make_perm_index error')
	res = 4
	return
endif
!if (ncells > nshow) write(*,*) 'start moving'
if (istep == 1) then
	nrepeat = 1
else
	nrepeat = 1
endif
!nrepeat = 60
!dt = DELTA_T/nrepeat
ndt = ndt + 1
!write(*,*) 'mover - grower: ndt: ',ndt
do irepeat = 1,nrepeat
t_fmover = 0
nit = 0
done = .false.
do while (.not.done)
	nit = nit + 1
	t0 = mytimer()
	call fmover(dt,done,ok)
	t1 = mytimer()
	tmover = tmover + t1 - t0
	if (.not.ok) then
		call logger('fmover error')
		res = 1
		return
	endif
	tnow = istep*DELTA_T + t_fmover
	ncells0 = ncells

	call GrowCells(radiation_dose,dt,changed,ok)
	if (.not.ok) then
		call logger('grower error')
		res = 3
		return
	endif
	tgrower = tgrower + mytimer() - t1
	radiation_dose = 0
	if (changed) then
		call make_perm_index(ok)
	endif
	t_fmover = t_fmover + dt
enddo
enddo
!write(*,'(a,2f8.2)') 't mover, grower: ',tmover,tgrower
t0 = mytimer()
call update_all_nbrlists
!write(*,'(a,f8.2)') 't update_nbr_lists: ',mytimer() - t0
!if (mod(istep,Nhop) == 0) then
!	! determine cell death and tagging for death 
!	call setup_grid_cells
!	call update_all_nbrlists
!	if (Ncells > 0) then
!		write(logmsg,'(a,5i8)') 'istep,Ncells,Nsteps,Nit,ndt: ',istep,Ncells,Nsteps,Nit,ndt
!		call logger(logmsg)
!	endif
!endif

call make_grid_flux_weights

if (ngaps > 2000 .or. mod(istep,nt_nbr) == 0) then
	if (ngaps > 2000) then
		call squeezer
	endif
	call setup_nbrlists(ok)
	if (.not.ok) then
		res = 5
		return
	endif
endif

if (use_metabolism) then
	call update_HIF1(DELTA_T)
endif

!write(*,'(a,3e12.3)') 'simulate_step: Cin: ',cp%Cin(1:3)
! Reaction-diffusion system
if (medium_change_step .or. chemo(DRUG_A)%present) then
	nt_diff = n_substeps
else
	nt_diff = 1
endif
dt = DELTA_T/nt_diff
call setup_grid_cells
call update_all_nbrlists
do it_diff = 1,nt_diff
	call diff_solver(dt,ok)
	if (.not.ok) then
		res = 7
		return
	endif
enddo
medium_change_step = .false.
call set_bdry_conc()

if (saveprofile%active) then
	if (istep*DELTA_T >= saveprofile%it*saveprofile%dt) then
		call WriteProfileData
		saveprofile%it = saveprofile%it + 1
		if (saveprofile%it > saveprofile%nt) then
			saveprofile%active = .false.
		endif
	endif
endif
if (saveslice%active) then
	if (istep*DELTA_T >= saveslice%it*saveslice%dt) then
		call WriteSliceData
		saveslice%it = saveslice%it + 1
		if (saveslice%it > saveslice%nt) then
			saveslice%active = .false.
		endif
	endif
endif

if (mod(istep,nt_hour) == 0) then
	if (ndoublings > 0) then
	    t_ave_double = doubling_time_sum/(ndoublings*3600)
	else
	    t_ave_double = 0
	endif
	write(logmsg,'(a,4i8,a,f5.2)') 'istep, hour, Ncells, ntagged: ',istep,istep/nt_hour,Ncells,Ndrug_tag(1,1),' t_double: ',t_ave_double
	call logger(logmsg)
	if (dbug) write(*,*) 'DRUG_A present: ',chemo(DRUG_A:DRUG_A+2)%present
!	write(*,'(7e11.3)') Caverage(NX/2,NY/2,11:17,DRUG_A)
!	write(*,'(7e11.3)') Caverage(NX/2,NY/2,11:17,DRUG_A+1)
!	write(*,'(7e11.3)') Caverage(NX/2,NY/2,11:17,DRUG_A+2)
	if (chemo(DRUG_A)%present) then
		call get_average_concs(DRUG_A,ave_cin, ave_cex)
		if (ave_cin(1) < 1.0e-6) then
			call ClearDrug(DRUG_A)
			chemo(DRUG_A:DRUG_A+2)%present = .false.
		endif
	endif
	if (chemo(DRUG_B)%present) then
		call get_average_concs(DRUG_B,ave_cin,ave_cex)
		if (ave_cin(1) < 1.0e-6) then
			call ClearDrug(DRUG_B)
			chemo(DRUG_B:DRUG_B+2)%present = .false.
		endif
	endif
	if (.not.use_TCP) then
		call get_concdata(nvars, ns, dxc, ex_conc)
	!	write(*,'(a,3f8.4)') 'cell #1: ',cell_list(1)%Cex(1),cell_list(1)%Cin(1),cell_list(1)%Cex(1)-cell_list(1)%Cin(1) 
	endif
!	call write_bdryconcs
    call showcells
endif
res = 0
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine update_HIF1(dt)
real(REAL_KIND) :: dt
real(REAL_KIND) :: HIF1, PDK1
integer :: kcell, ityp
type(metabolism_type), pointer :: mp
type(cell_type), pointer :: cp

do kcell = 1,nlist
	cp => cell_list(kcell)
	if (cp%state == DEAD) cycle
	mp => cp%metab
	HIF1 = mp%HIF1
	ityp = cp%celltype
	call analyticSetHIF1(ityp,cp%Cin(OXYGEN),HIF1,dt)
	mp%HIF1 = HIF1
	PDK1 = mp%PDK1
	call analyticSetPDK1(ityp,HIF1,PDK1,dt)
	mp%PDK1 = PDK1
enddo
end subroutine

!-----------------------------------------------------------------------------------------
! Find average drug concs in cells
!-----------------------------------------------------------------------------------------
subroutine get_average_concs(ichemo,cin,cex)
integer :: ichemo
real(REAL_KIND) :: cin(0:2), cex(0:2)
type(cell_type), pointer :: cp
integer :: kcell, n

n = 0
cin = 0
cex = 0
do kcell = 1,nlist
	cp => cell_list(kcell)
	if (cp%state == DEAD) cycle
	n = n+1
	cin = cin + cp%Cin(ichemo:ichemo+2)
	cex = cex + cp%Cex(ichemo:ichemo+2)
enddo
cin = cin/n
cex = cex/n
if (dbug) then
	write(*,'(a,3e12.3)') 'Average IC drug concs: ',cin
	write(*,'(a,3e12.3)') 'Average EC drug concs: ',cex
endif
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine showcell(kcell)
integer :: kcell
real(REAL_KIND) :: Vn   ! to normalise volumes
type(cell_type), pointer :: cp

Vn = Vdivide0/1.6
cp => cell_list(kcell)
write(nflog,'(a,i6,4e12.3)') 'kcell, volume, divide_volume, dVdt, divide_time: ', &
                kcell, cp%V/Vn, cp%divide_volume/Vn, cp%dVdt/Vn, cp%divide_time
!write(*,'(a,i6,4e12.3)') 'kcell, volume, divide_volume, dVdt, divide_time: ', &
!                kcell, cp%V/Vn, cp%divide_volume/Vn, cp%dVdt/Vn, cp%divide_time
end subroutine

!-----------------------------------------------------------------------------------------
! Average volumes etc
!-----------------------------------------------------------------------------------------
subroutine showcells
integer :: kcell, n
real(REAL_KIND) :: Vsum,divVsum,dVdtsum,divtsum
real(REAL_KIND) :: Vn   ! to normalise volumes
type(cell_type), pointer :: cp

Vn = Vdivide0/1.6
Vsum=0
divVsum=0
dVdtsum=0
divtsum=0
n=0
do kcell = 1,nlist
    cp => cell_list(kcell)
    if (cp%state == DEAD) cycle
    n = n+1
    Vsum = Vsum + cp%V/Vn
    divVsum = divVsum + cp%divide_volume/Vn
    dVdtsum = dVdtsum + cp%dVdt/Vn
    divtsum = divtsum + cp%divide_time
enddo
!write(nflog,'(a,4e12.3)') 'ave volume, divide_volume, dVdt, divide_time: ', Vsum/n,divVsum/n,dVdtsum/n,divtsum/n
!write(*,'(a,4e12.3)') 'ave volume, divide_volume, dVdt, divide_time: ', Vsum/n,divVsum/n,dVdtsum/n,divtsum/n
end subroutine

!----------------------------------------------------------------------------------
!----------------------------------------------------------------------------------
subroutine squeezer
integer :: igap, kcell, nlist0

if (ngaps == 0) return
write(*,*) 'squeezer: ngaps: ',ngaps
!do igap = 1,ngaps
!	kcell = gaplist(igap)
!	if (
!	cell_list(kcell) = cell_list(nlist)
!	nlist = nlist - 1
!enddo
nlist0 = nlist
nlist = 0
do kcell = 1,nlist0
	if (cell_list(kcell)%state == DEAD) cycle
	nlist = nlist + 1
	cell_list(nlist) = cell_list(kcell)
enddo
write(nflog,*) '==squeezer== nlist0, nlist: ',nlist0, nlist
ngaps = 0
end subroutine

!----------------------------------------------------------------------------------
!----------------------------------------------------------------------------------
subroutine ProcessEvent(radiation_dose)
real(REAL_KIND) :: radiation_dose
integer :: kevent, ichemo, idrug, im, nmetab
real(REAL_KIND) :: V, C(MAX_CHEMO)
type(event_type) :: E

!write(logmsg,*) 'ProcessEvent'
!call logger(logmsg)
do kevent = 1,Nevents
	E = event(kevent)
	if (t_simulation >= E%time .and. .not.E%done) then
		write(nflog,'(a,i3,2f8.0,i3,2f10.4)') 'Event: ',E%etype,t_simulation,E%time,E%ichemo,E%volume,E%conc
		if (E%etype == RADIATION_EVENT) then
			radiation_dose = E%dose
			write(logmsg,'(a,f8.0,f8.3)') 'RADIATION_EVENT: time, dose: ',t_simulation,E%dose
			call logger(logmsg)
		elseif (E%etype == MEDIUM_EVENT) then
			write(logmsg,'(a,f8.0,f8.3,2f8.4)') 'MEDIUM_EVENT: time, volume, O2medium: ',t_simulation,E%volume,E%O2medium
			call logger(logmsg)
			C = 0
			C(OXYGEN) = E%O2medium
			chemo(OXYGEN)%bdry_conc = E%O2medium
			C(GLUCOSE) = chemo(GLUCOSE)%bdry_conc
			V = E%volume
			call MediumChange(V,C)
		elseif (E%etype == DRUG_EVENT) then
			C = 0
			C(OXYGEN) = E%O2conc
			chemo(OXYGEN)%bdry_conc = E%O2conc
			C(GLUCOSE) = chemo(GLUCOSE)%bdry_conc
			ichemo = E%ichemo
			idrug = E%idrug
			C(ichemo) = E%conc
			V = E%volume
			write(logmsg,'(a,2f8.3)') 'DRUG_EVENT: volume, conc: ',E%volume,E%conc
			call logger(logmsg)
			! set %present
			chemo(ichemo)%present = .true.
			chemo(ichemo)%bdry_conc = 0
			nmetab = drug(idrug)%nmetabolites
			do im = 1,nmetab
				if (chemo(ichemo + im)%used) then
					chemo(ichemo + im)%present = .true.
					chemo(ichemo + im)%bdry_conc = 0
				endif
			enddo
			call MediumChange(V,C)
		endif
		event(kevent)%done = .true.
	endif
enddo
end subroutine

!-----------------------------------------------------------------------------------------
! If the volume removed is Vr, the fraction of constituent mass that is retained
! in the medium is (Vm - Vr)/Vm.  The calculation does not apply to oxygen.
! Usually Vr = Ve.
! With grid, need to exclude gridcells that contain cells and those interior gridcells.
! The total mass in the external gridcells is determined, and a fraction of it may be retained 
!-----------------------------------------------------------------------------------------
subroutine MediumChange(Ve,Ce)
real(REAL_KIND) :: Ve, Ce(:)
real(REAL_KIND) :: R, Vm_old, Vm_new, Vr, Vkeep, Vblob, fkeep
integer :: kcell, site(3), siteb(3), ixb, iyb, izb, izb_1, izb_2, ichemo, ic
!integer, allocatable :: zinrng(:,:,:)
real(REAL_KIND), allocatable :: exmass(:), exconc(:)
real(REAL_KIND), pointer :: Cave(:,:,:)

write(nflog,*) 'MediumChange: ',Ve,Ce
allocate(ngcells(NXB,NYB,NZB))
!allocate(zinrng(2,NXB,NYB))
allocate(exmass(MAX_CHEMO))
allocate(exconc(MAX_CHEMO))
exmass = 0
ngcells = 0
Vblob = 0
! need to identify gridcells containing cells, to exclude from the calculation
do kcell = 1,nlist
	if (cell_list(kcell)%state == DEAD) cycle
	site = cell_list(kcell)%site	! this gives the location in the fine grid
	call getSiteb(site,siteb)
	ngcells(siteb(1),siteb(2),siteb(3)) = ngcells(siteb(1),siteb(2),siteb(3)) + 1
enddo
do ixb = 1,NXB
	do iyb = 1,NYB
		! Need to include all internal empty gridcells (necrotic)
		! First find first and last non-empty gridcells
		izb_1 = 0
		izb_2 = 0
		do izb = 1,NZB
			if (ngcells(ixb,iyb,izb) > 0) then
				if (izb_1 == 0) izb_1 = izb
				izb_2 = izb
			endif
		enddo
		if (izb_1 /= 0) then
			Vblob = Vblob + izb_2 - izb_1 + 1
		endif
		do izb = 1,NZB
			if (izb_1 /= 0 .and. (izb >= izb_1 .and. izb <= izb_2)) cycle	! blob gridcells
			do ichemo = 1,MAX_CHEMO
				if (.not.chemo(ichemo)%used) cycle
				exmass(ichemo) = exmass(ichemo) + dxb3*chemo(ichemo)%Cave_b(ixb,iyb,izb)
			enddo
		enddo
!		zinrng(:,ixb,iyb) = [izb_1,izb_2]
	enddo
enddo
Vblob = dxb3*Vblob

R = getRadius()
Vblob = (4./3.)*PI*R**3			! cm3
Vm_old = total_volume - Vblob	! this is the initial external (medium) volume
Vr = min(Vm_old,Ve)				! this is the amount of medium volume that is exchanged
Vkeep = Vm_old - Vr
fkeep = Vkeep/Vm_old			! fraction of initial medium volume that is kept. i.e. fraction of exmass(:)
Vm_new = Ve + Vkeep				! new medium volume
total_volume = Vm_new + Vblob	! new total volume
!write(*,'(6f8.4)') total_volume,Vblob,Vm_old,Vkeep,fkeep,Vm_new
! Concentrations in the gridcells external to the blob (those with no cells) are set
! to the values of the mixture of old values and added medium. 
do ichemo = 1,MAX_CHEMO
	if (.not.chemo(ichemo)%used) cycle
	if (Ce(ichemo) == 0) then
		exconc(ichemo) = 0
	else
		exconc(ichemo) = (Ve*Ce(ichemo) + fkeep*exmass(ichemo))/Vm_new
	endif
	write(nflog,'(a,i2,6f8.4)') 'MediumChange: exconc: ',ichemo,Ce(ichemo),Ve,fkeep,exmass(ichemo),Vm_new,exconc(ichemo)
enddo
do ichemo = 1,MAX_CHEMO
	if (.not.chemo(ichemo)%used) cycle
	write(nflog,*) 'Setting Cave_b: ',ichemo,exconc(ichemo)
	chemo(ichemo)%Cave_b(:,:,:) = exconc(ichemo)
	chemo(ichemo)%Cprev_b(:,:,:) = exconc(ichemo)
	chemo(ichemo)%medium_Cbnd = exconc(ichemo)
	Caverage(:,:,:,ichemo) = exconc(ichemo)		! This is approximately OK, because the medium volume is so much greater than the blob volume
	! Try this
	Cflux(:,:,:,ichemo) = 0
	chemo(ichemo)%Fcurr_b = 0
	call update_Cex(ichemo)
enddo
medium_change_step = .true.

! Use interpolation to determine concentrations on the boundary of the fine grid
do ic = 1,nchemo
	ichemo = chemomap(ic)
	Cave => Caverage(:,:,:,ichemo)
	call interpolate_Cave(ichemo, Cave, chemo(ichemo)%Cave_b)
enddo
deallocate(ngcells)
!deallocate(zinrng)
deallocate(exmass)
deallocate(exconc)
end subroutine

!-----------------------------------------------------------------------------------------
! Determine which gridcell in the big grid (siteb(:)) the fine gridcell site(:) falls in.
! Note that %site(:) and %centre(:,:) always refer to the fine grid, because it is
! assumed that cells never move outside the fine grid.
!-----------------------------------------------------------------------------------------
subroutine getSiteb(site,siteb)
integer :: site(3), siteb(3)
integer :: ixb1, iyb1, ixb, iyb, izb

ixb1 = (NXB+1)/2 - (NX-1)/(2*NRF)
ixb = ixb1 + (site(1)-1)/NRF
iyb1 = (NYB+1)/2 - (NY-1)/(2*NRF)
iyb = iyb1 + (site(2)-1)/NRF
izb = 1 + (site(3)-1)/NRF
siteb = [ixb, iyb, izb]
end subroutine

!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------
subroutine get_dimensions(NX_dim, NY_dim, NZ_dim, nsteps_dim, deltat, maxchemo, nextra, cused, dfraction, deltax) BIND(C)
!DEC$ ATTRIBUTES DLLEXPORT :: get_dimensions
use, intrinsic :: iso_c_binding
integer(c_int) :: NX_dim,NY_dim,NZ_dim,nsteps_dim, maxchemo, nextra
real(c_double) :: deltat, dfraction, deltax
logical(c_bool) :: cused(*)
integer :: ichemo

NX_dim = NX
NY_dim = NY
NZ_dim = NZ
nsteps_dim = nsteps
deltat = DELTA_T
deltax = DELTA_X
maxchemo = MAX_CHEMO
nextra = N_EXTRA
do ichemo = 1,MAX_CHEMO
	cused(ichemo+1) = .true.	!chemo(ichemo)%used
enddo
cused(1) = .true.			! CFSE
cused(MAX_CHEMO+2) = .true.	! Growth rate
dfraction = 1.0	!2*cell_radius/DELTA_X
end subroutine

!--------------------------------------------------------------------------------
! TC = tumour cell
!--------------------------------------------------------------------------------
subroutine get_scene(nTC_list,TC_list,dropped,droppedcentre) BIND(C)
!DEC$ ATTRIBUTES DLLEXPORT :: get_scene 
use, intrinsic :: iso_c_binding
integer(c_int) :: nTC_list
type(celldata_type) :: TC_list(*)
logical(c_bool) :: dropped
real(c_double) :: droppedcentre(3)
integer :: kcell, nspheres, is
real(REAL_KIND) :: centre(3)
type(cell_type), pointer :: cp

dropped = is_dropped
droppedcentre = 0
nTC_list = 0
do kcell = 1,nlist
	cp => cell_list(kcell)
	if (cp%state == DEAD) cycle
	if (cp%Iphase) then
		nspheres = 1
	else
		nspheres = 2
	endif
	do is = 1,nspheres
		nTC_list = nTC_list + 1
		TC_list(nTC_list)%tag = nTC_list
		TC_list(nTC_list)%radius = 10000*cp%radius(is)		! cm -> um
		centre = 10000*cp%centre(:,is)	! cm -> um
		droppedcentre = droppedcentre + centre
		! rotate to make z the vertical axis on screen
		TC_list(nTC_list)%centre(1) = centre(2)
		TC_list(nTC_list)%centre(2) = centre(3)	! invert in z-dirn 
		TC_list(nTC_list)%centre(3) = centre(1)
		TC_list(nTC_list)%celltype = cp%celltype
		TC_list(nTC_list)%status = 0
		if (cp%mitosis > 0) then
			TC_list(nTC_list)%status = 1
		endif
!		write(nflog,'(2i6,4e12.3)') nTC_list,TC_list(nTC_list)%tag,TC_list(nTC_list)%radius,TC_list(nTC_list)%centre
	enddo
enddo
droppedcentre = droppedcentre/nTC_list
!write(nflog,*) 'nTC_list: ',nTC_list
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine get_DLL_build_version(version_array,array_len) BIND(C) 
!DEC$ ATTRIBUTES DLLEXPORT :: get_dll_build_version
use, intrinsic :: iso_c_binding
character(c_char) :: version_array(12)
integer(c_int) :: array_len
integer :: k

dll_version = DLL_BUILD_VERSION
gui_version = GUI_BUILD_VERSION
do k = 1,12
	version_array(k) = dll_version(k:k)
	if (version_array(k) == ' ') then
		version_array(k) = char(0)
		array_len = k
		exit
	endif
enddo
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine Execute(ncpu,infile_array,inbuflen,outfile_array,outbuflen,centre) BIND(C)
!subroutine Execute() BIND(C)
!DEC$ ATTRIBUTES DLLEXPORT :: execute
use, intrinsic :: iso_c_binding
!use :: seq_mod
character(c_char) :: infile_array(128), outfile_array(128)
integer(c_int) :: ncpu, inbuflen, outbuflen
real(c_double) :: centre(*)
character*(128) :: infile, outfile
integer :: i, res
logical :: ok, isopen

!call test_seq
!stop

infile = ''
do i = 1,inbuflen
	infile(i:i) = infile_array(i)
enddo
outfile = ''
do i = 1,outbuflen
	outfile(i:i) = outfile_array(i)
enddo

awp_0%is_open = .false.
awp_1%is_open = .false.
par_zig_init = .false.
logfile = 'vspheroid.log'
inquire(unit=nflog,OPENED=isopen)
if (.not.isopen) then
    open(nflog,file=logfile,status='replace')
endif
#if defined(OPENMP) || defined(_OPENMP)
    write(logmsg,'(a)') 'Executing with OpenMP'
	call logger(logmsg)
#endif

write(logmsg,*) 'inputfile:  ', infile
call logger(logmsg)
write(logmsg,*) 'outputfile: ', outfile 
call logger(logmsg)
if (use_TCP) then
	write(nflog,*) 'call connecter'
	call connecter(ok)
	if (.not.ok) then
		call logger('Failed to make TCP connections')
		return
	endif
endif
call Setup(ncpu,infile,outfile,ok)
if (ok) then
	clear_to_send = .true.
	simulation_start = .true.
	istep = 0
	res = 0
else
	call logger('=== Setup failed ===')
	res = 1
	stop
endif
execute_t1 = wtime()
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine DisableTCP
!DEC$ ATTRIBUTES DLLEXPORT :: disableTCP
!DEC$ ATTRIBUTES STDCALL, REFERENCE, MIXED_STR_LEN_ARG, ALIAS:"DISABLETCP" :: disableTCP

use_TCP = .false.   ! because this is called from spheroid_main()	
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine Connection(awp,port,error)
TYPE(winsockport) :: awp
integer :: port, error
integer :: address = 0
!!!character*(64) :: ip_address = "127.0.0.1"C      ! need a portable way to make a null-terminated C string
character*(64) :: host_name = "localhost"

if (.not.winsock_init(1)) then
    call logger("winsock_init failed")
    stop
endif

awp%handle = 0
awp%host_name = host_name
awp%ip_port = port
awp%protocol = IPPROTO_TCP
call Set_Winsock_Port (awp,error)

if (.not.awp%is_open) then
    write(nflog,*) 'Error: connection: awp not open: ',port
else
    write(nflog,*) 'connection: awp open: ',port, error
endif
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine Connecter(ok)
logical :: ok
integer :: error

! Main connection
ok = .true.
error = 0
call Connection(awp_0,TCP_PORT_0,error)
if (awp_0%handle < 0 .or. error /= 0) then
    write(logmsg,'(a)') 'TCP connection to TCP_PORT_0 failed'
    call logger(logmsg)
    ok = .false.
    return
endif
if (.not.awp_0%is_open) then
	write(logmsg,'(a)') 'No connection to TCP_PORT_0'
    call logger(logmsg)
    ok = .false.
    return
endif
write(logmsg,'(a)') 'Connected to TCP_PORT_0  '
call logger(logmsg)

if (use_CPORT1) then
	call connection(awp_1,TCP_PORT_1,error)
	if (awp_1%handle < 0 .or. error /= 0) then
		write(logmsg,'(a)') 'TCP connection to TCP_PORT_1 failed'
		call logger(logmsg)
		ok = .false.
		return
	endif
	if (.not.awp_1%is_open) then
		write(logmsg,'(a)') 'No connection to TCP_PORT_1'
		call logger(logmsg)
		ok = .false.
		return
	endif
	write(logmsg,'(a)') 'Connected to TCP_PORT_1  '
	call logger(logmsg)
endif
! Allow time for completion of the connection
call sleeper(2)
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine terminate_run(res) BIND(C)
!DEC$ ATTRIBUTES DLLEXPORT :: terminate_run 
use, intrinsic :: iso_c_binding
integer(c_int) :: res
character*(8), parameter :: quit = '__EXIT__'
integer :: error, i

call logger('terminate_run')
call Wrapup

if (res == 0) then
	call logger(' Execution successful!')
elseif (res == -1) then
	call logger(' Execution stopped')
elseif (res == 2) then
	call logger(' No more live cells')
elseif (res == 6) then
	call logger(' Spheroid size limit reached')
elseif (res == 1) then
	call logger(' === Execution failed === ERROR in fmover')
elseif (res == 3) then
	call logger(' === Execution failed === ERROR in GrowCells')
elseif (res == 4) then
	call logger(' === Execution failed === ERROR in make_perm_index')
elseif (res == 5) then
	call logger(' === Execution failed === ERROR in setup_nbrlists')
elseif (res == 7) then
	call logger(' === Execution failed === ERROR in diff_solver')
endif
write(logmsg,'(a,f10.2)') 'Execution time (min): ',(wtime() - execute_t1)/60
call logger(logmsg)

!close(nflog)

if (use_TCP) then
	if (stopped) then
	    call winsock_close(awp_0)
	    if (use_CPORT1) call winsock_close(awp_1)
	else
	    call winsock_send(awp_0,quit,8,error)
	    call winsock_close(awp_0)
		if (use_CPORT1) then
			call winsock_send(awp_1,quit,8,error)
			call winsock_close(awp_1)
		endif
	endif
endif
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine Wrapup
integer :: ierr
logical :: isopen

call logger('doing wrapup ...')
ierr = 0
!call logger('deallocated all arrays')

! Close all open files
inquire(unit=nfout,OPENED=isopen)
if (isopen) then
	close(nfout)
	call logger('closed nfout')
endif
inquire(nfres,OPENED=isopen)
if (isopen) close(nfres)
call logger('closed files')

if (par_zig_init) then
	call par_zigfree
endif
call logger('freed par_zig')
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine testresults
integer :: k1, k2
real(REAL_KIND) :: R1, R2, c1(3), c2(3), v(3), d2, d

write(*,*) 'testresults'
do k1 = 1,nlist
	if (cell_list(k1)%state == DEAD) cycle
	R1 = cell_list(k1)%radius(1)
	c1 = cell_list(k1)%centre(:,1)
	write(*,'(a,i4,3f6.2)') 'k1, c: ',k1,c1
	do k2 = k1+1, nlist
		if (cell_list(k2)%state == DEAD) cycle
		R2 = cell_list(k2)%radius(1)
		c2 = cell_list(k2)%centre(:,1)
		v = c2 - c1
		d2 = dot_product(v,v)
		d = sqrt(d2)
		write(*,'(2i4,f8.2)') k1,k2,d
	enddo
enddo

end subroutine


end module