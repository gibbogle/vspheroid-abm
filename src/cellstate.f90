! Cancer cell state development

module cellstate
use global
use chemokine
use metabolism
use nbr
use cycle_mod

implicit none

integer :: kcell_dividing = 0
real(REAL_KIND) :: total_dose_time

contains

!-----------------------------------------------------------------------------------------
! Need to initialize site and cell concentrations when a cell divides and when there is
! cell death.
!-----------------------------------------------------------------------------------------
subroutine GrowCells(dt,changed, ok)
real(REAL_KIND) :: dt
logical :: changed, ok
type(cell_type),pointer :: cp

!call logger('GrowCells')
!tnow = istep*DELTA_T
ok = .true.
!if (dose > 0) then
!	call Irradiation(dose, ok)
!	if (.not.ok) return
!endif
changed = .false.
!call grower(dt,changed,OK)
call new_grower(dt,changed,OK)
if (.not.ok) return
if (use_death) then
	call CellDeath(dt,changed, ok)
	if (.not.ok) return
endif
!if (use_migration) then
!	call CellMigration(ok)
!	if (.not.ok) return
!endif
!cp => cell_list(kcell_test)
!write(nflog,*) 'ndivided: ',ndivided
!write(nflog,'(a,2i6,2f8.3,e12.3)') 'kcell_test: ',kcell_test,cp%phase,tnow/3600,cp%mitosis,cp%dVdt
end subroutine

!-----------------------------------------------------------------------------------------
! The O2 concentration to use with cell kcell is either the intracellular concentration,
! or is use_extracellular_O2, the corresponding extracellular concentration
!-----------------------------------------------------------------------------------------
subroutine getO2conc(cp, C_O2)
type(cell_type), pointer :: cp
real(REAL_KIND) :: C_O2

if (use_extracellular_O2 .and. istep > 1) then		! fix 30/04/2015
	C_O2 = cp%Cex(OXYGEN)
else
	C_O2 = cp%Cin(OXYGEN)
endif
end subroutine

!-----------------------------------------------------------------------------------------
! The glucose concentration to use with cell kcell is either the intracellular concentration,
! or if use_extracellular_O2 (!), the corresponding extracellular concentration
!-----------------------------------------------------------------------------------------
subroutine getGlucoseconc(cp, C_glucose)
type(cell_type), pointer :: cp
real(REAL_KIND) :: C_glucose

if (use_extracellular_O2 .and. istep > 1) then
	C_glucose = cp%Cex(GLUCOSE)
else
	C_glucose = cp%Cin(GLUCOSE)
endif
end subroutine

!-----------------------------------------------------------------------------------------
! Irradiate cells with dose (and duration tmin to be added)
!-----------------------------------------------------------------------------------------
subroutine Irradiation(dose,ok)
real(REAL_KIND) :: dose
logical :: ok
integer :: kcell, site(3), iv, ityp, n, kpar=0
real(REAL_KIND) :: C_O2, SER, SER_OER(2), p_death, p_recovery, R, kill_prob, tmin
type(cell_type), pointer :: cp
type(cycle_parameters_type), pointer :: ccp

ok = .true.
if (dose == 0) then
	return
endif
call logger('Irradiation')
!tnow = istep*DELTA_T	! seconds
if (use_volume_method) then
    do kcell = 1,nlist
        if (colony_simulation) then
            cp => ccell_list(kcell)
        else
            cp => cell_list(kcell)
        endif
	    if (cp%state == DEAD) cycle
	    if (cp%radiation_tag) cycle	! we do not tag twice (yet)
	    ityp = cp%celltype
	    call getO2conc(cp,C_O2)
	    ! Compute sensitisation SER
	    if (Ndrugs_used > 0) then
	        SER = getSER(cp,C_O2)
	    else
    	    SER = 1.0
    	endif
!	    SER = 1.0
!	    do idrug = 1,Ndrugs_used
!		    ichemo = 4 + 3*(idrug-1)
!		    if (.not.chemo(ichemo)%present) cycle
!		    do im = 0,2
!			    ichemo = 4 + 3*(idrug-1) + im
!			    if (drug(idrug)%sensitises(ityp,im)) then
!				    Cs = cp%Cin(ichemo)	! concentration of drug/metabolite in the cell
!				    SER_max0 = drug(idrug)%SER_max(ityp,im)
!				    SER_Km = drug(idrug)%SER_Km(ityp,im)
!				    SER_KO2 = drug(idrug)%SER_KO2(ityp,im)
!				    SERmax = (Cs*SER_max0 + SER_Km)/(Cs + SER_Km)
!				    SER = SER*(C_O2 + SER_KO2*SERmax)/(C_O2 + SER_KO2)
!			    endif
!		    enddo
!	    enddo
	    call get_kill_probs(ityp,dose,C_O2,SER,p_recovery,p_death)
	    kill_prob = 1 - p_recovery
	    R = par_uni(kpar)
	    if (R < kill_prob) then
		    cp%radiation_tag = .true.	! volume_method
		    Nradiation_tag(ityp) = Nradiation_tag(ityp) + 1
		    cp%p_rad_death = p_death
		    if (LQ(ityp)%growth_delay_N > 0 .and. cp%Iphase) then
			    cp%growth_delay = .true.
			    cp%dt_delay = LQ(ityp)%growth_delay_factor*dose
			    cp%N_delayed_cycles_left = LQ(ityp)%growth_delay_N
		    else
			    cp%growth_delay = .false.
		    endif
	    elseif (use_radiation_growth_delay_all .and. LQ(ityp)%growth_delay_N > 0) then
		    cp%growth_delay = .true.
		    cp%dt_delay = LQ(ityp)%growth_delay_factor*dose
		    cp%N_delayed_cycles_left = LQ(ityp)%growth_delay_N
	    else
		    cp%growth_delay = .false.
	    endif
    enddo
else
	n = 0
    tmin = 1.0      ! for now...
    do kcell = 1,nlist
        if (colony_simulation) then
            cp => ccell_list(kcell)
        else
            cp => cell_list(kcell)
        endif
	    if (cp%state == DEAD) cycle
	    ityp = cp%celltype
	    ccp => cc_parameters(ityp)
	    call getO2conc(cp,C_O2)
	    ! Compute sensitisation SER 
	    if (Ndrugs_used > 0) then
	        SER = getSER(cp,C_O2)
	    else
    	    SER = 1.0
    	endif
        SER_OER(1) = SER*(LQ(ityp)%OER_am*C_O2 + LQ(ityp)%K_ms)/(C_O2 + LQ(ityp)%K_ms)      ! OER_alpha
        SER_OER(2) = SER*(LQ(ityp)%OER_bm*C_O2 + LQ(ityp)%K_ms)/(C_O2 + LQ(ityp)%K_ms)      ! OER_beta
        SER_OER = 1
        call radiation_damage(cp, ccp, dose, SER_OER(1), tmin)
	    if (cp%irrepairable) then	! irrepairable damage 
			if (.not.cp%radiation_tag) then
				n = n+1
				cp%radiation_tag = .true.	! tagged, but not DYING yet
			    Nradiation_tag(ityp) = Nradiation_tag(ityp) + 1
!				write(nflog,*) 'Tagged with IRL: ',istep,cp%N_IRL,Nradiation_tag(ityp),kcell,cp%phase
			endif	
			if (cp%phase == S_phase .or. cp%phase == M_phase) then	! kill immediately, otherwise dies at mitosis
				call CellDies(kcell,.false.)
!				Nradiation_dead(ityp) = Nradiation_dead(ityp) + 1
!				write(nflog,*) 'S-phase - killed'
			endif
		endif
    enddo
    call check_radiation
    write(logmsg,*) '# irrepairable: ',n
    call logger(logmsg)
endif
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine check_radiation
type(cell_type), pointer :: cp
type(cycle_parameters_type), pointer :: ccp
integer :: kcell, n
real(REAL_KIND) :: aveCh1, aveCh2, aveIRL, avePL

n = 0
aveCh1 = 0 
aveCh2 = 0 
aveIRL = 0 
avePL = 0
do kcell = 1,nlist
    cp => cell_list(kcell)
    if (cp%state == DEAD .or. cp%state == DYING) cycle
	n = n+1
	aveCh1 = aveCh1 + cp%N_Ch1
	aveCh2 = aveCh2 + cp%N_Ch2
	aveIRL = aveIRL + cp%N_IRL
	avePL = avePL + cp%N_PL
enddo
write(nflog,*) 'check_radiation: n: ',n
write(nflog,'(a,f9.3)') 'aveCh1: ',aveCh1/n
write(nflog,'(a,f9.3)') 'aveCh2: ',aveCh2/n
write(nflog,'(a,f9.3)') 'aveIRL: ',aveIRL/n
write(nflog,'(a,f9.3)') 'avePL: ',avePL/n
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
function getSER(cp,C_O2) result(SER)
type(cell_type), pointer :: cp
real(REAL_KIND) ::C_O2, SER
real(REAL_KIND) :: Cs							! concentration of radiosensitising drug
real(REAL_KIND) :: SER_max0, SER_Km, SER_KO2	! SER parameters of the drug
real(REAL_KIND) :: SERmax						! max sensitisation at the drug concentration
integer :: ityp, idrug, ichemo, im

ityp = cp%celltype
SER = 1.0
do idrug = 1,Ndrugs_used
    ichemo = 4 + 3*(idrug-1)
    if (.not.chemo(ichemo)%present) cycle
    do im = 0,2
	    ichemo = 4 + 3*(idrug-1) + im
	    if (drug(idrug)%sensitises(ityp,im)) then
		    Cs = cp%Cin(ichemo)	! concentration of drug/metabolite in the cell
		    SER_max0 = drug(idrug)%SER_max(ityp,im)
		    SER_Km = drug(idrug)%SER_Km(ityp,im)
		    SER_KO2 = drug(idrug)%SER_KO2(ityp,im)
		    SERmax = (Cs*SER_max0 + SER_Km)/(Cs + SER_Km)
		    SER = SER*(C_O2 + SER_KO2*SERmax)/(C_O2 + SER_KO2)
	    endif
    enddo
enddo
end function

!-----------------------------------------------------------------------------------------
! A cell that receives a dose of radiation either recovers completely before reaching 
! mitosis or retains damage that has a probability of causing cell death during mitosis.
! A damaged cell that does not die at this point passes the damage on to the progeny cells.
! The probability of complete recovery = p_recovery = p_r
! The probability of death for a damaged cell at mitosis = p_death = p_d
! To ensure that the short-term death probability is consistent with the previous
! LQ formulation, we require p_d(1-p_r) = kill_prob as previously calculated.
! If p_d is determined (currently it is fixed), then 1-p_r = kill_prob/p_d,
! therefore p_r = 1 - kill_prob/p_d
!-----------------------------------------------------------------------------------------
subroutine get_kill_probs(ityp,dose,C_O2,SER,p_recovery,p_death)
integer :: ityp
real(REAL_KIND) :: dose, C_O2, p_recovery, p_death
real(REAL_KIND) :: SER						! sensitisation
real(REAL_KIND) :: OER_alpha_d, OER_beta_d, expon, kill_prob_orig

OER_alpha_d = dose*(LQ(ityp)%OER_am*C_O2 + LQ(ityp)%K_ms)/(C_O2 + LQ(ityp)%K_ms)
OER_beta_d = dose*(LQ(ityp)%OER_bm*C_O2 + LQ(ityp)%K_ms)/(C_O2 + LQ(ityp)%K_ms)
!expon = LQ(ityp)%alpha_H*OER_alpha_d + LQ(ityp)%beta_H*OER_alpha_d**2		! 07/08/2015

OER_alpha_d = OER_alpha_d*SER
OER_beta_d = OER_beta_d*SER

expon = LQ(ityp)%alpha_H*OER_alpha_d + LQ(ityp)%beta_H*OER_beta_d**2
p_recovery = exp(-expon)	! = SF
p_death = LQ(ityp)%death_prob

!kill_prob_orig = 1 - exp(-expon)
!call get_pdeath(ityp,dose,C_O2,p_death)
!p_recovery = 1 - kill_prob_orig/p_death
end subroutine

!-----------------------------------------------------------------------------------------
! This is the probability of death at time of division of cell that received a dose of 
! radiation and did not recover.
! In general one would expect this to depend of damage, i.e. on dose and C_O2, but
! for now it is a constant for a cell type
!-----------------------------------------------------------------------------------------
subroutine get_pdeath(ityp,dose,C_O2,p_death)
integer :: ityp
real(REAL_KIND) :: dose, C_O2, p_death

p_death = LQ(ityp)%death_prob
end subroutine

!-----------------------------------------------------------------------------------------
! Cells move to preferable nearby sites.
! For now this is turned off - need to formulate a sensible interpretation of "preferable"
!-----------------------------------------------------------------------------------------
subroutine CellMigration(ok)
logical :: ok
integer :: kcell, j, indx, site0(3), site(3), jmax
real(REAL_KIND) :: C0(MAX_CHEMO), C(MAX_CHEMO), v0, v, vmax, d0, d

call logger('CellMigration is not yet implemented')
ok = .false.
return

end subroutine


!-----------------------------------------------------------------------------------------
! Cells can be tagged to die, or finally die of anoxia or glucosia, or they can be tagged 
! for death at division time if the drug is effective.
!-----------------------------------------------------------------------------------------
subroutine CellDeath(dt,changed,ok)
real(REAL_KIND) :: dt
logical :: changed, ok
integer :: kcell, nlist0, site(3), i, ichemo, idrug, im, ityp, n, killmodel, kpar=0 
real(REAL_KIND) :: C_O2, C_glucose, kmet, Kd, dMdt, Cdrug, n_O2, kill_prob, dkill_prob, death_prob, survival_prob, u
!logical :: anoxia_death, aglucosia_death
real(REAL_KIND) :: delayed_death_prob
logical :: flag
type(cell_type), pointer :: cp
type(drug_type), pointer :: dp
real(REAL_KIND) :: p_tag = 0.3
type(cycle_parameters_type), pointer :: ccp
!logical :: allow_anoxia = .true.

!write(nflog,*) 'CellDeath'
ok = .true.
if (colony_simulation) then
    return
endif

!tnow = istep*DELTA_T	! seconds
!anoxia_death = chemo(OXYGEN)%controls_death
!aglucosia_death = chemo(GLUCOSE)%controls_death
nlist0 = nlist
flag = .false.
n = 0
do kcell = 1,nlist
    cp => cell_list(kcell)
	if (cp%state == DEAD) cycle
	ityp = cp%celltype
	ccp => cc_parameters(ityp)
	if (cp%state == DYING) then
		delayed_death_prob = ccp%apoptosis_rate*dt/3600
	    if (par_uni(kpar) < delayed_death_prob) then	
			call CellDies(kcell,.true.)
			changed = .true.
!			NATP_dead(ityp) = NATP_dead(ityp) + 1		! note that death from anoxia here means metabolic death 
		endif
		cycle
	endif
	if (cp%irrepairable) n = n+1
	call getO2conc(cp,C_O2)
	if (use_metabolism) then
		if (cp%metab%A_rate < ATPs) then
			cp%ATP_tag = .true.
			NATP_tag(ityp) = NATP_tag(ityp) + 1
			cp%dVdt = 0
			call CellDies(kcell,.false.)
!			cp%state = DYING
!			Ndying(ityp) = Ndying(ityp) + 1
			cycle
		endif
	else
	
#if 0
		if (cp%anoxia_tag) then
!			if (tnow >= cp%t_anoxia_die) then
			if (tnow >= cp%t_anoxia_die .and. par_uni(kpar) < p_tag*dt/DELTA_T) then
				call CellDies(kcell,cp)
				changed = .true.
				Nanoxia_dead(ityp) = Nanoxia_dead(ityp) + 1
				cycle
			endif
		else
			if (anoxia_death .and. C_O2 < anoxia_threshold) then
				cp%t_anoxia = cp%t_anoxia + dt
				if (cp%t_anoxia > t_anoxia_limit) then
					cp%anoxia_tag = .true.						! tagged to die later
					cp%t_anoxia_die = tnow + anoxia_death_delay	! time that the cell will die
					Nanoxia_tag(ityp) = Nanoxia_tag(ityp) + 1
				endif
			else
				cp%t_anoxia = 0
			endif
		endif
		call getGlucoseconc(cp,C_glucose)
		if (cp%aglucosia_tag) then
			if (tnow >= cp%t_aglucosia_die) then
				call CellDies(kcell,cp)
				changed = .true.
				Naglucosia_dead(ityp) = Naglucosia_dead(ityp) + 1
				cycle
			endif
		else
			if (aglucosia_death .and. C_O2 < aglucosia_threshold) then
				cp%t_aglucosia = cp%t_aglucosia + dt
				if (cp%t_aglucosia > t_aglucosia_limit) then
					cp%aglucosia_tag = .true.						    ! tagged to die later
					cp%t_aglucosia_die = tnow + aglucosia_death_delay	! time that the cell will die
					Naglucosia_tag(ityp) = Naglucosia_tag(ityp) + 1
				endif
			else
				cp%t_aglucosia = 0
			endif
		endif
#endif
	endif
	do idrug = 1,ndrugs_used	
		ichemo = DRUG_A + 3*(idrug-1)	
		if (.not.chemo(ichemo)%present) cycle
		if (cp%drug_tag(idrug)) cycle	! don't tag more than once
		dp => drug(idrug)
		kill_prob = 0
		death_prob = 0
		survival_prob = 1
		do im = 0,2
			if (.not.dp%kills(ityp,im)) cycle
			Cdrug = cp%Cin(ichemo + im)
			if (Cdrug == 0) cycle
!			if (.not.flag) write(nflog,*) 'idrug,im,kcell: ',idrug,im,kcell
			killmodel = dp%kill_model(ityp,im)		! could use %drugclass to separate kill modes 
			Kd = dp%Kd(ityp,im)
			n_O2 = dp%n_O2(ityp,im)
!			Cdrug = 0.54
			kmet = (1 - dp%C2(ityp,im) + dp%C2(ityp,im)*dp%KO2(ityp,im)**n_O2/(dp%KO2(ityp,im)**n_O2 + C_O2**n_O2))*dp%Kmet0(ityp,im)
			if (.not.flag .and. C_O2 < 0.01) write(nflog,'(a,6e12.3)') 'C_O2,Kmet0,kmet,Kd,Cdrug,dt: ',C_O2,dp%Kmet0(ityp,im),kmet,Kd,Cdrug,dt
			dMdt = kmet*Cdrug
			call getDrugKillProb(killmodel,Kd,dMdt,Cdrug,dt,dkill_prob)
			survival_prob = survival_prob*(1 - dkill_prob)
			if (.not.flag .and. C_O2 < 0.01) write(nflog,'(a,2e12.3)') 'dkill_prob: ',dkill_prob,survival_prob
			death_prob = max(death_prob,dp%death_prob(ityp,im))
		enddo
		kill_prob = 1 - survival_prob
		if (kill_prob == 0) cycle
		if (.not.flag .and. idrug==1) then
			total_dose_time = total_dose_time + dt
			write(nflog,'(a,2f7.4,2f8.1)') 'kill_prob, SF: ',kill_prob,1-kill_prob,dt,total_dose_time
		endif
		u = par_uni(kpar)
	    if (u < kill_prob) then
			cp%p_drug_death(idrug) = death_prob
			cp%drug_tag(idrug) = .true.
            Ndrug_tag(idrug,ityp) = Ndrug_tag(idrug,ityp) + 1
		endif
		flag = .true.
	enddo
enddo
!if (n > 0) then
!	write(logmsg,*) '# of irrepairable left: ',n
!	call logger(logmsg)
!endif
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine getDrugKillProb(kill_model,Kd,dMdt,Cdrug,dt,dkill_prob)
integer :: kill_model
real(REAL_KIND) :: Kd, dMdt, Cdrug, dt, dkill_prob
real(REAL_KIND) :: SF, dtstep, kill_prob, c
integer :: Nsteps, istep

if (kill_model == 1) then
	c = Kd*dMdt
elseif (kill_model == 2) then
	c = Kd*dMdt*Cdrug
elseif (kill_model == 3) then
	c = Kd*dMdt**2
elseif (kill_model == 4) then
	c = Kd*Cdrug
elseif (kill_model == 5) then
	c = Kd*Cdrug**2
endif
SF = exp(-c*dt)
dkill_prob = 1 - SF
end subroutine

!-----------------------------------------------------------------------------------------
! If the dying cell site is less than a specified fraction f_migrate of the blob radius,
! the site migrates towards the blob centre.
! %indx -> 0
! If the site is on the boundary, it is removed from the boundary list, and %indx -> OUTSIDE_TAG
! The cell contents should be released into the site.
!-----------------------------------------------------------------------------------------
subroutine CellDies(kcell,now)
integer :: kcell
logical :: now
type(cell_type), pointer :: cp
integer :: ityp, idrug

!write(*,*) 'CellDies: ',kcell
cp => cell_list(kcell)
ityp = cp%celltype
if (.not.now) then
	cp%state = DYING
	Ndying(ityp) = Ndying(ityp) + 1
!	write(nflog,*) 'CellDies: DYING: ',kcell,Ndying(ityp)
	return
endif
if (cp%state == DYING) then
	Ndying(ityp) = Ndying(ityp) - 1
endif
Ncells = Ncells - 1
Ncells_type(ityp) = Ncells_type(ityp) - 1
Ndead(ityp) = Ndead(ityp) + 1
!write(nflog,*) 'CellDies: DEAD: ',kcell,Ndead(ityp)
if (cp%ATP_tag) then
	NATP_tag(ityp) = NATP_tag(ityp) - 1
	NATP_dead(ityp) = NATP_dead(ityp) + 1
endif
do idrug = 1,ndrugs_used
	if (cp%drug_tag(idrug)) then
		Ndrug_tag(idrug,ityp) = Ndrug_tag(idrug,ityp) - 1
		Ndrug_dead(idrug,ityp) = Ndrug_dead(idrug,ityp) + 1
	endif
enddo
if (cp%radiation_tag) then
	Nradiation_tag(ityp) = Nradiation_tag(ityp) - 1
	Nradiation_dead(ityp) = Nradiation_dead(ityp) + 1
endif
ngaps = ngaps + 1
if (ngaps > max_ngaps) then
    write(logmsg,'(a,i6,i6)') 'CellDies: ngaps > max_ngaps: ',ngaps,max_ngaps
    call logger(logmsg)
    stop
endif
gaplist(ngaps) = kcell
cp%state = DEAD

end subroutine

!-----------------------------------------------------------------------------------------
! Cell growth, death and division are handled here.  Division occurs when cell volume 
! exceeds the divide volume. 
! As the cell grows we need to adjust both Cin and Cex to maintain mass conservation.
! GROWTH DELAY
! When a cell has received a dose of radiation (or possibly drug - not yet considered)
! the cycle time is increased by an amount that depends on the dose.  The delay may be
! transmitted to progeny cells.
! 
! NOTE: now the medium concentrations are not affected by cell growth
!-----------------------------------------------------------------------------------------
subroutine new_grower(dt, changed, ok)
real(REAL_KIND) :: dt
logical :: changed, ok
integer :: nlivecells
integer, allocatable :: livelist(:)
integer :: k, kcell, idrug, ityp, n
logical :: divide, drug_killed(2), radiation_killed
integer, parameter :: MAX_DIVIDE_LIST = 10000
integer :: ndivide, divide_list(MAX_DIVIDE_LIST)
integer, parameter :: MAX_DRUG_LIST = 10000
integer :: ndrug_kill(2), drug_kill_list(2,MAX_DRUG_LIST)
integer, parameter :: MAX_RADIATION_LIST = 10000
integer :: nradiation_kill, radiation_kill_list(MAX_RADIATION_LIST)
type(cell_type), pointer :: cp

ok = .true.
changed = .false.
ndivide = 0
!tnow = istep*DELTA_T !+ t_fmover
!if (colony_simulation) write(*,*) 'grower: ',nlist0,use_volume_method,tnow

allocate(livelist(nlist))
nlivecells = 0
do kcell = 1,nlist
	if (colony_simulation) then
	    cp => ccell_list(kcell)
	else
    	cp => cell_list(kcell)
    endif
	if (cp%state == DEAD) cycle
	nlivecells = nlivecells + 1
	livelist(nlivecells) = kcell
enddo

ndivide = 0
ndrug_kill = 0
nradiation_kill = 0
call omp_set_num_threads(Mnodes)
!$omp parallel do private(kcell, divide, drug_killed, radiation_killed, idrug)
do k = 1,nlivecells
	kcell = livelist(k)
	call omp_grower(dt, kcell, divide, drug_killed, radiation_killed)
	if (divide) then
		ndivide = ndivide + 1
		if (ndivide > MAX_DIVIDE_LIST) then
		    write(logmsg,*) 'Error: growcells: MAX_DIVIDE_LIST exceeded: ',MAX_DIVIDE_LIST
		    call logger(logmsg)
		    ok = .false.
!		    return
		endif
		divide_list(ndivide) = kcell
	endif
	do idrug = 1,2
		if (drug_killed(idrug)) then
			ndrug_kill(idrug) = ndrug_kill(idrug) + 1
			drug_kill_list(idrug,ndrug_kill(idrug)) = kcell
		endif
	enddo
	if (radiation_killed) then
		nradiation_kill = nradiation_kill + 1
		radiation_kill_list(nradiation_kill) = kcell
	endif
enddo
!$omp end parallel do

do idrug = 1,2
	do k = 1,ndrug_kill(idrug)
		changed = .true.
		kcell = drug_kill_list(idrug,k)
		if (colony_simulation) then
			cp => ccell_list(kcell)
		else
    		cp => cell_list(kcell)
		endif
		ityp = cp%celltype
		call CellDies(kcell,.false.)
!		Ndrug_dead(idrug,ityp) = Ndrug_dead(idrug,ityp) + 1
	enddo
enddo

do k = 1,nradiation_kill
	changed = .true.
	kcell = radiation_kill_list(k)
	if (colony_simulation) then
		cp => ccell_list(kcell)
	else
		cp => cell_list(kcell)
	endif
	ityp = cp%celltype
	call CellDies(kcell,.false.)
!	Nradiation_dead(ityp) = Nradiation_dead(ityp) + 1
enddo

do k = 1,ndivide
	changed = .true.
	kcell = divide_list(k)
	if (colony_simulation) then
	    cp => ccell_list(kcell)
	else
    	cp => cell_list(kcell)
    endif
	call divider(kcell, ok)
	if (.not.ok) return
enddo
deallocate(livelist)
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine omp_grower(dt, kcell, divide, drug_killed, radiation_killed)
real(REAL_KIND) :: dt
integer :: kcell
logical :: divide, drug_killed(2), radiation_killed

integer :: k, ityp, idrug, prev_phase, kpar=0
type(cell_type), pointer :: cp
type(cycle_parameters_type), pointer :: ccp
real(REAL_KIND) :: rr(3), c(3), rad, d_desired, R, pdeath, mitosis_duration, deviation, rad0
logical :: mitosis_entry, in_mitosis, tagged
logical :: dbug

dbug = (kcell == -977)
divide = .false.
drug_killed = .false.
radiation_killed = .false.

	if (colony_simulation) then
	    cp => ccell_list(kcell)
	else
    	cp => cell_list(kcell)
    endif
	ityp = cp%celltype
	ccp => cc_parameters(ityp)
	mitosis_entry = .false.
	mitosis_duration = ccp%T_M
	in_mitosis = .false.
!	tagged = cp%anoxia_tag .or. cp%aglucosia_tag .or. (cp%state == DYING)
!	if (tagged) then
	if (cp%state == DYING) then
		cp%dVdt = 0
		return
	endif
	if (use_volume_method) then
!        if (colony_simulation) then
!            write(*,'(a,i6,L2,2e12.3)') 'kcell: ',kcell,cp%Iphase,cp%V,cp%divide_volume
!        endif
	    if (cp%Iphase) then
		    call growcell(cp,dt)
		    if (cp%V > cp%divide_volume) then	! time to enter mitosis
!    	        mitosis_entry = .true.
!				cp%Iphase = .false.
	            in_mitosis = .true.
!				cp%mitosis = 0
!				cp%t_start_mitosis = tnow
				
				cp%Iphase = .false.
				cp%mitosis = 0
				cp%t_start_mitosis = tnow
				ncells_mphase = ncells_mphase + 1
				cp%nspheres = 2
				call get_random_vector3(rr)	! set initial axis direction
				cp%d = 0.1*small_d
				c = cp%centre(:,1)
				cp%centre(:,1) = c + (cp%d/2)*rr
				cp%centre(:,2) = c - (cp%d/2)*rr
				cp%d_divide = 2.0**(2./3)*cp%radius(1)
	        endif
	    else
	        in_mitosis = .true.
	    endif
	else
		if (dbug) write(*,*) 'phase: ',cp%phase
	    prev_phase = cp%phase
!	    if (cp%dVdt == 0) then
!			write(nflog,*) 'dVdt=0: kcell, phase: ',kcell,cp%phase 
!		endif
	    if (cp%dVdt > 0) then
			call timestep(cp, ccp, dt)
		endif
!        if (.not.cp%radiation_tag .and.(cp%NL2(1) > 0 .or. cp%NL2(2) > 0)) then	! irrepairable damage
!			! For now, tag for death at mitosis
!			cp%radiation_tag = .true.	! new_grower
!		    Nradiation_tag(ityp) = Nradiation_tag(ityp) + 1
!		endif
        if (cp%phase >= M_phase) then
            if (prev_phase == Checkpoint2) then		! this is mitosis entry
                if (.not.cp%radiation_tag .and. cp%N_PL > 0) then		! lesions still exist, no time for repair
					cp%radiation_tag = .true.
				    Nradiation_tag(ityp) = Nradiation_tag(ityp) + 1
!				    write(nflog,*) 'Tagged with N_PL: ',istep,cp%N_PL,Nradiation_tag(ityp)
				endif
				
				cp%Iphase = .false.
				cp%nspheres = 2
				cp%mitosis = 0
				cp%t_start_mitosis = tnow
				ncells_mphase = ncells_mphase + 1
				call get_random_vector3(rr)	! set initial axis direction
				rrsum = rrsum + rr
				cp%d = 0.1*small_d
				c = cp%centre(:,1)
				cp%centre(:,1) = c + (cp%d/2)*rr
				cp%centre(:,2) = c - (cp%d/2)*rr
				cp%d_divide = 2.0**(2./3)*cp%radius(1)
				
				! For cells with Ch1 or Ch2, check for death
				if (cp%N_Ch1 > 0 .and. .not.cp%radiation_tag) then
					pdeath = 1 - ccp%psurvive_Ch1**cp%N_Ch1
					R = par_uni(kpar)
					if (R < pdeath) then				
						cp%radiation_tag = .true.
						Nradiation_tag(ityp) = Nradiation_tag(ityp) + 1
!						write(nflog,*) 'Tagged with N_Ch1: ',cp%N_Ch1,Nradiation_tag(ityp)
					endif
				endif
				if (cp%N_Ch2 > 0 .and. .not.cp%radiation_tag) then
					pdeath = 1 - ccp%psurvive_Ch2**cp%N_Ch2
					R = par_uni(kpar)
					if (R < pdeath) then				
						cp%radiation_tag = .true.
						Nradiation_tag(ityp) = Nradiation_tag(ityp) + 1
!						write(nflog,*) 'Tagged with N_Ch2: ',cp%N_Ch2,Nradiation_tag(ityp)
					endif
				endif
            endif
            in_mitosis = .true.
        endif
		if (cp%phase < Checkpoint2 .and. cp%phase /= Checkpoint1) then
		    call growcell(cp,dt)
		endif	
	endif
	
	if (in_mitosis) then
		do idrug = 1,ndrugs_used
			if (cp%drug_tag(idrug)) then
				drug_killed(idrug) = .true.
				return
			endif
		enddo

		cp%mitosis = (tnow - cp%t_start_mitosis)/mitosis_duration
		d_desired = max(cp%mitosis*cp%d_divide,small_d)
		rr = cp%centre(:,2) - cp%centre(:,1)
		cp%d = sqrt(dot_product(rr,rr))
		deviation = (cp%d - d_desired)/d_desired
!		rr = rr/cp%d	! axis direction
		c = (cp%centre(:,1) + cp%centre(:,2))/2
		cp%site = c/DELTA_X + 1
!		call cubic_solver(d_desired,cp%V,rad)
        if (cp%d > cp%d_divide) then
            divide = .true.
        else
		    call cubic_solver(d_desired,cp%V,rad)
		    cp%radius = rad
!		    if (cp%d > 0.75*cp%d_divide/2 .and. abs(deviation) > 0.2) then
!		        write(logmsg,'(a,i8,f8.4)') 'deviation: ',kcell,deviation
!			    call logger(logmsg)
!		    endif
!		    if (cp%d > 2*rad) then	! completion of mitosis - note that this overrides cp%phase 
!			    write(logmsg,'(a,i6,2e12.3,f7.3)') 'divides: d > 2*rad: ',kcell,cp%d,rad,cp%mitosis
!			    call logger(logmsg)
!			    write(nflog,'(a,6e12.3)') 'centre(1),centre(2): ',cp%centre(:,1),cp%centre(:,2)
!			    divide = .true.
!		    endif
		endif
			
		if (use_volume_method) then
			if (cp%growth_delay) then
				if (cp%G2_M) then
					if (tnow > cp%t_growth_delay_end) then
						cp%G2_M = .false.
					else
!						cycle
						return
					endif
				else
					cp%t_growth_delay_end = tnow + cp%dt_delay
					cp%G2_M = .true.
!					cycle
					return
				endif
			endif
			! try moving death prob test to here
			if (cp%radiation_tag) then
				R = par_uni(kpar)
				if (R < cp%p_rad_death) then
					radiation_killed = .true.
					return
				endif
			endif		
		else
		    ! Check for cell death by radiation lesions
			if (cp%radiation_tag) then
				radiation_killed = .true.
				return
			endif		        
		endif
		
!		write(*,*) 'istep, kcell, mitosis: ',istep,kcell,cp%mitosis
        if (cp%mitosis >= 1) then
			divide = .true.
		endif
	endif
	
end subroutine


!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine old_growcell(cp, dt, c_rate, r_mean)
type(cell_type), pointer :: cp
real(REAL_KIND) :: dt, c_rate, r_mean
real(REAL_KIND) :: Cin_0(NCONST), Cex_0(NCONST)		! at some point NCONST -> MAX_CHEMO
real(REAL_KIND) :: dVdt,  Vin_0, dV, metab_O2, metab_glucose, metab, dVdt_new
logical :: glucose_growth
integer :: C_option = 1	! we must use this

glucose_growth = chemo(GLUCOSE)%controls_growth
Cin_0 = cp%Cin
metab_O2 = O2_metab(Cin_0(OXYGEN))	! Note that currently growth depends only on O2
metab_glucose = glucose_metab(Cin_0(GLUCOSE))
if (glucose_growth) then
	metab = metab_O2*metab_glucose
else
	metab = metab_O2
endif
dVdt = get_dVdt(cp,metab)
if (suppress_growth) then	! for checking solvers
	dVdt = 0
endif
Cex_0 = cp%Cex
cp%dVdt = dVdt
Vin_0 = cp%V
dV = dVdt*dt
cp%V = Vin_0 + dV
if (C_option == 1) then
	! Calculation based on transfer of an extracellular volume dV with constituents, i.e. holding extracellular concentrations constant
	cp%Cin = (Vin_0*Cin_0 + dV*Cex_0)/(Vin_0 + dV)
elseif (C_option == 2) then
	! Calculation based on change in volumes without mass transfer of constituents
	cp%Cin = Vin_0*Cin_0/(Vin_0 + dV)
endif
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine growcell(cp, dt)
type(cell_type), pointer :: cp
real(REAL_KIND) :: dt
real(REAL_KIND) :: Cin_0(NCONST), Cex_0(NCONST)		! at some point NCONST -> MAX_CHEMO
real(REAL_KIND) :: dVdt,  Vin_0, dV, metab_O2, metab_glucose, metab
logical :: oxygen_growth, glucose_growth
integer :: ityp
integer :: C_option = 1	! we must use this

!if (use_cell_cycle .and. .not.(cp%phase == G1_phase .or. cp%phase == S_phase .or. cp%phase == G2_phase)) then
!	write(*,*) 'no growth - phase'
!	return
!endif
ityp = cp%celltype
if (colony_simulation) then		! FIX THIS LATER
    metab = 1
    dVdt = get_dVdt(cp,metab)
else
	if (use_metabolism) then
		cp%metab%Itotal = cp%metab%Itotal + dt*cp%metab%I_rate
		! need to set cp%dVdt from cp%metab%I_rate
!		dVdt = max_growthrate(ityp)*cp%metab%I_rate/cp%metab%I_rate_max
!		dVdt = cp%growth_rate_factor*dVdt
!		cp%V = cp%V + cp%dVdt*dt
		metab = cp%metab%I_rate/cp%metab%I_rate_max
		dVdt = get_dVdt(cp,metab)
		Cin_0 = cp%Cin
!		if (kcell_now == 1) then
!			write(*,*) 'cp%growth_rate_factor: ',cp%growth_rate_factor
!			write(*,'(a,4e12.3)') 'dVdt: ',max_growthrate(ityp),cp%metab%I_rate/cp%metab%I_rate_max,cp%dVdt,get_dVdt(cp,metab)
!			write(*,'(a,i3,4e12.3)') 'phase,V,Vdivide,I..: ', cp%phase,cp%V,cp%divide_volume,cp%metab%Itotal,cp%metab%I2Divide
!		endif
	else
		oxygen_growth = chemo(OXYGEN)%controls_growth
		glucose_growth = chemo(GLUCOSE)%controls_growth
		Cin_0 = cp%Cin
		metab_O2 = O2_metab(Cin_0(OXYGEN))	! Note that currently growth depends only on O2
		metab_glucose = glucose_metab(Cin_0(GLUCOSE))
		if (oxygen_growth .and. glucose_growth) then
			metab = metab_O2*metab_glucose
		elseif (oxygen_growth) then
			metab = metab_O2
		elseif (glucose_growth) then
			metab = metab_glucose
		endif
		dVdt = get_dVdt(cp,metab)
		if (suppress_growth) then	! for checking solvers
			dVdt = 0
		endif
	endif
    Cex_0 = cp%Cex
endif
cp%dVdt = dVdt
Vin_0 = cp%V
dV = dVdt*dt
if (use_cell_cycle .and. .not.(cp%phase == G1_phase .or. cp%phase == S_phase .or. cp%phase == G2_phase)) then
    dV = 0
endif
cp%V = Vin_0 + dV
if (.not.colony_simulation) then
    if (C_option == 1) then
        ! Calculation based on transfer of an extracellular volume dV with constituents, i.e. holding extracellular concentrations constant
        cp%Cin = (Vin_0*Cin_0 + dV*Cex_0)/(Vin_0 + dV)
    elseif (C_option == 2) then
        ! Calculation based on change in volumes without mass transfer of constituents
        cp%Cin = Vin_0*Cin_0/(Vin_0 + dV)
    endif
endif
end subroutine

!-----------------------------------------------------------------------------------------
! Need to account for time spent in mitosis.  Because there is no growth during mitosis,
! the effective divide_time must be reduced by mitosis_duration.
! Note that TERMINAL_MITOSIS is the only option.
!-----------------------------------------------------------------------------------------
function get_dVdt(cp, metab) result(dVdt)
type(cell_type), pointer :: cp
real(REAL_KIND) :: metab, dVdt
integer :: ityp
real(REAL_KIND) :: r_mean, c_rate, mitosis_duration
type(cycle_parameters_type), pointer :: ccp

ityp = cp%celltype
ccp => cc_parameters(ityp)
mitosis_duration = ccp%T_M
if (use_cell_cycle) then
!    if (cp%phase == G1_phase .or. cp%phase == S_phase .or. cp%phase == G2_phase) then
        if (use_constant_growthrate) then
            dVdt = max_growthrate(ityp)
        else
            dVdt = metab*max_growthrate(ityp)
        endif
!    else
!        dVdt = 0
!    endif
else
    if (use_V_dependence) then
	    if (use_constant_divide_volume) then
		    dVdt = metab*log(2.0)*cp%V/(cp%divide_time - mitosis_duration)
	    else
		    c_rate = log(2.0)/(divide_time_mean(ityp) - mitosis_duration)
		    dVdt = c_rate*cp%V*metab
	    endif
    else
	    if (use_constant_divide_volume) then
		    dVdt = 0.5*metab*Vdivide0/(cp%divide_time  - mitosis_duration)
	    else
		    ityp = cp%celltype
		    r_mean = 0.5*Vdivide0/(divide_time_mean(ityp) - mitosis_duration)
		    dVdt = r_mean*metab
	    endif
    endif
endif
dVdt = dVdt/cp%fg   ! apply growth factor here instead of in cycle
!if (cp%dVdt == 0) then
!    write(nflog,*) 'get_dVdt: = 0'
!    write(*,*) 'get_dVdt: dVdt = 0'
!    write(*,'(a,i6,4e12.3)') 'kcell, metab, Vdivide0: ',kcell_now, metab, Vdivide0,cp%Cin(OXYGEN),cp%Cin(GLUCOSE)
!    stop
!endif
end function

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine SetInitialGrowthRate
integer :: kcell, ityp
real(REAL_KIND) :: C_O2, C_glucose, metab, metab_O2, metab_glucose, dVdt
logical :: oxygen_growth, glucose_growth
type(cell_type), pointer :: cp

write(nflog,*) 'SetInitialGrowthRate'
oxygen_growth = chemo(OXYGEN)%controls_growth
glucose_growth = chemo(GLUCOSE)%controls_growth
do kcell = 1,nlist
    if (colony_simulation) then
        cp => ccell_list(kcell)
    else
        cp => cell_list(kcell)
    endif
	if (cp%state == DEAD) cycle
	if (use_metabolism) then
		metab = 1
	else
		C_O2 = chemo(OXYGEN)%bdry_conc
		C_glucose = chemo(GLUCOSE)%bdry_conc
		if (oxygen_growth .and. glucose_growth) then
			metab_O2 = O2_metab(C_O2)
			metab_glucose = glucose_metab(C_glucose)
			metab = metab_O2*metab_glucose
		elseif (oxygen_growth) then
			metab_O2 = O2_metab(C_O2)
			metab = metab_O2
		elseif (glucose_growth) then
			metab_glucose = glucose_metab(C_glucose)
			metab = metab_glucose
		endif
	endif
	dVdt = get_dVdt(cp,metab)
	if (suppress_growth) then	! for checking solvers
		dVdt = 0
	endif
	cp%dVdt = dVdt
enddo
end subroutine

!-----------------------------------------------------------------------------------------
! A single cell is replaced by two.
!-----------------------------------------------------------------------------------------
subroutine old_divider(kcell1, ok)
integer :: kcell1
logical :: ok
integer :: kcell2, ityp, nbrs0
real(REAL_KIND) :: r(3), c(3), cfse0, cfse1, V0, Tdiv, gfactor
type(cell_type), pointer :: cp1, cp2

!write(*,*) 'divider:'
!write(logmsg,*) 'divider: ',kcell1 
!call logger(logmsg)
ok = .true.
tnow = istep*DELTA_T
cp1 => cell_list(kcell1)
if (ngaps > 0) then
    kcell2 = gaplist(ngaps)
    ngaps = ngaps - 1
else
	nlist = nlist + 1
	if (nlist > MAX_NLIST) then
		ok = .false.
		call logger('Error: Maximum number of cells MAX_NLIST has been exceeded.  Increase and rebuild.')
		return
	endif
	kcell2 = nlist
endif
ncells = ncells + 1
ityp = cp1%celltype
ncells_type(ityp) = ncells_type(ityp) + 1
ncells_mphase = ncells_mphase - 1
cp2 => cell_list(kcell2)
cp1%state = ALIVE
cp1%generation = cp1%generation + 1
V0 = cp1%V/2
cp1%V = V0
cp1%site = cp1%centre(:,1)/DELTA_X + 1
cp1%d = 0
cp1%birthtime = tnow
!cp1%divide_volume = get_divide_volume1()
cp1%divide_volume = get_divide_volume(ityp,V0,Tdiv,gfactor)
cp1%divide_time = Tdiv
cp1%d_divide = (3*cp1%divide_volume/PI)**(1./3.)
cp1%mitosis = 0
cp1%t_divide_last = tnow
cfse0 = cp1%CFSE
cp1%CFSE = generate_CFSE(cfse0/2)
cfse1 = cfse0 - cp1%CFSE

cp1%drug_tag = .false.
!cp1%anoxia_tag = .false.
!cp1%t_anoxia = 0
!cp1%aglucosia_tag = .false.
!cp1%t_aglucosia = 0

if (cp1%growth_delay) then
	cp1%N_delayed_cycles_left = cp1%N_delayed_cycles_left - 1
	cp1%growth_delay = (cp1%N_delayed_cycles_left > 0)
!	write(*,*) 'growth_delay cell divides: ',kcell1,kcell2,cp1%N_delayed_cycles_left
endif
cp1%G2_M = .false.

nbrs0 = cp1%nbrs
cp1%nbrs = nbrs0 + 1
cp1%nbrlist(cp1%nbrs)%indx = kcell2
cp1%nbrlist(cp1%nbrs)%contact = .false.
cp1%nbrlist(cp1%nbrs)%contact(1,1) = .true.

cp2%ID = cp1%ID
cp2%celltype = cp1%celltype
cp2%state = ALIVE
cp2%generation = cp1%generation
cp2%V = V0
cp2%radius(1) = cp1%radius(2)
cp2%centre(:,1) = cp1%centre(:,2)
cp2%site = cp2%centre(:,1)/DELTA_X + 1
cp2%d = 0
cp2%birthtime = tnow
!cp2%divide_volume = get_divide_volume1()
cp2%divide_volume = get_divide_volume(ityp,V0,Tdiv,gfactor)
cp2%divide_time = Tdiv
cp2%d_divide = (3*cp2%divide_volume/PI)**(1./3.)
cp2%mitosis = 0
cp2%t_divide_last = tnow
cp2%CFSE = cfse1

cp2%ID = cp1%ID
cp2%p_rad_death = cp1%p_rad_death
cp2%radiation_tag = cp1%radiation_tag
if (cp2%radiation_tag) then
	Nradiation_tag(ityp) = Nradiation_tag(ityp) + 1
endif
cp2%drug_tag = .false.
!cp2%anoxia_tag = .false.
!cp2%t_anoxia = 0
!cp2%aglucosia_tag = .false.
!cp2%t_aglucosia = 0

cp2%growth_delay = cp1%growth_delay
if (cp2%growth_delay) then
	cp2%dt_delay = cp1%dt_delay
	cp2%N_delayed_cycles_left = cp1%N_delayed_cycles_left
endif
cp2%G2_M = .false.

cp2%Cin = cp1%Cin
cp2%Cex = cp1%Cex

!allocate(cp2%nbrlist(MAX_NBRS))
! Need to set up neighbour lists for both cp1 and cp2 <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< DO THIS
! Use neighbour lists of neighbours
! First simple test:
cp2%nbrlist(1:nbrs0) = cp1%nbrlist(1:nbrs0)
cp2%nbrs = nbrs0 + 1
cp2%nbrlist(cp2%nbrs)%indx = kcell1
cp2%nbrlist(cp2%nbrs)%contact = .false.
cp2%nbrlist(cp2%nbrs)%contact(1,1) = .true.
cp1%Iphase = .true.
cp1%nspheres = 1
cp2%Iphase = .true.
cp2%nspheres = 1
cp2%ndt = cp1%ndt
call update_nbrlist(kcell1)
call update_nbrlist(kcell2)
! Note: any cell that has kcell1 in its nbrlist now needs to add kcell2 to the nbrlist.
end subroutine

!-----------------------------------------------------------------------------------------
! New version
!-----------------------------------------------------------------------------------------
subroutine divider(kcell1, ok)
integer :: kcell1
logical :: ok
integer :: kcell2, ityp, nbrs0
real(REAL_KIND) :: r(3), c(3), cfse0, cfse2, V0, Tdiv, gfactor
type(cell_type), pointer :: cp1, cp2
type(cycle_parameters_type), pointer :: ccp

!write(logmsg,*) 'divider: ',kcell1 
!call logger(logmsg)
ok = .true.
if (colony_simulation) then
    cp1 => ccell_list(kcell1)
else
	cp1 => cell_list(kcell1)
endif
!write(nflog,'(a,i6,4e12.3)') 'divider: ',kcell1,cp1%mitosis,cp1%V,cp1%radius(1:2)
!write(nflog,'(a,i6,f6.2)') 'divided in time: ',kcell1,(tnow - cp1%t_divide_last)/3600
if (ngaps > 0) then
    kcell2 = gaplist(ngaps)
    ngaps = ngaps - 1
else
	nlist = nlist + 1
	if (nlist > MAX_NLIST) then
		ok = .false.
		call logger('Error: Maximum number of cells MAX_NLIST has been exceeded.  Increase and rebuild.')
		return
	endif
	kcell2 = nlist
endif
ncells = ncells + 1
ityp = cp1%celltype
ccp => cc_parameters(ityp)
ncells_type(ityp) = ncells_type(ityp) + 1
ncells_mphase = ncells_mphase - 1
if (colony_simulation) then
    cp2 => ccell_list(kcell2)
else
	cp2 => cell_list(kcell2)
endif

cp1%state = ALIVE
cp1%generation = cp1%generation + 1
V0 = cp1%V/2
cp1%V = V0
cp1%site = cp1%centre(:,1)/DELTA_X + 1
cp1%d = 0
cp1%birthtime = tnow
!cp1%divide_volume = get_divide_volume1()
cp1%divide_volume = get_divide_volume(ityp,V0,Tdiv,gfactor)
cp1%divide_time = Tdiv
cp1%fg = gfactor
cp1%d_divide = (3*cp1%divide_volume/PI)**(1./3.)
if (use_metabolism) then	! Fraction of I needed to divide = fraction of volume needed to divide
	cp1%metab%I2Divide = get_I2Divide(cp1)
	cp1%metab%Itotal = 0
endif
cp1%mitosis = 0
cfse0 = cp1%CFSE
cp1%CFSE = generate_CFSE(cfse0/2)
cfse2 = cfse0 - cp1%CFSE

cp1%drug_tag = .false.
!cp1%anoxia_tag = .false.
!cp1%t_anoxia = 0
!cp1%aglucosia_tag = .false.
!cp1%t_aglucosia = 0

if (cp1%growth_delay) then
	cp1%N_delayed_cycles_left = cp1%N_delayed_cycles_left - 1
	cp1%growth_delay = (cp1%N_delayed_cycles_left > 0)
endif
!if (use_metabolism .and. cp1%dVdt == 0) then
!	cp1%dVdt = (cp1%metab%I_rate/cp1%metab%I_rate_max)*max_growthrate(ityp)
!endif
!cp1%G1_time = tnow + (max_growthrate(ityp)/cp1%dVdt)*cp1%fg*ccp%T_G1    ! time spend in G1 varies inversely with dV/dt
cp1%G1_time = tnow + (max_growthrate(ityp)/cp1%dVdt)*ccp%T_G1    ! time spend in G1 varies inversely with dV/dt

if (.not.colony_simulation) then
    nbrs0 = cp1%nbrs
    cp1%nbrs = nbrs0 + 1
    cp1%nbrlist(cp1%nbrs)%indx = kcell2
    cp1%nbrlist(cp1%nbrs)%contact = .false.
    cp1%nbrlist(cp1%nbrs)%contact(1,1) = .true.
endif
cp1%G2_M = .false.
cp1%Iphase = .true.
cp1%nspheres = 1
cp1%phase = G1_phase

ndoublings = ndoublings + 1
doubling_time_sum = doubling_time_sum + tnow - cp1%t_divide_last
cp1%t_divide_last = tnow

! Second cell
cp2 = cp1

! These are the variations from cp1
cp2%divide_volume = get_divide_volume(ityp,V0,Tdiv,gfactor)
cp2%divide_time = Tdiv
cp2%fg = gfactor
if (use_metabolism) then	! Fraction of I needed to divide = fraction of volume needed to divide
	cp2%metab%I2Divide = get_I2Divide(cp2)
	cp2%metab%Itotal = 0
!	cp2%growth_rate_factor = get_growth_rate_factor()
!	cp2%ATP_rate_factor = get_ATP_rate_factor()
endif
cp2%CFSE = cfse2
if (cp2%radiation_tag) then
	Nradiation_tag(ityp) = Nradiation_tag(ityp) + 1
endif
if (.not.colony_simulation) then
    cp2%d_divide = (3*cp2%divide_volume/PI)**(1./3.)
    cp2%radius(1) = cp1%radius(2)
    cp2%centre(:,1) = cp1%centre(:,2)
    cp2%site = cp2%centre(:,1)/DELTA_X + 1
    cp2%nbrlist(cp2%nbrs)%indx = kcell1
    call update_nbrlist(kcell1)
    call update_nbrlist(kcell2)
endif
! Note: any cell that has kcell1 in its nbrlist now needs to add kcell2 to the nbrlist.
ndivided = ndivided + 1
end subroutine

!----------------------------------------------------------------------------------
! This is used to adjust a cell's growth rate to introduce some random variability
! into divide times.
! dr = 0 suppresses variability
!----------------------------------------------------------------------------------
!function get_growth_rate_factor() result(r)
!real(REAL_KIND) :: r
!real(REAL_KIND) :: dr = 0.0	! was 0.2
!integer :: kpar = 0
!
!r = 1 + (par_uni(kpar) - 0.5)*dr
!end function

!----------------------------------------------------------------------------------
! This is used to adjust a cell's ATP thresholds to introduce some random variability
! into transitions to no-growth and death.
! Too much variability?
! dr = 0 suppresses variability - not needed in spheroid
!----------------------------------------------------------------------------------
function get_ATP_rate_factor() result(r)
real(REAL_KIND) :: r
real(REAL_KIND) :: dr = 0.0	! was 0.2
integer :: kpar = 0

r = 1 + (par_uni(kpar) - 0.5)*dr
end function

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
function get_divide_volume1() result(vol)
real(REAL_KIND) :: vol
integer :: kpar = 0
real(REAL_KIND) :: U

U = par_uni(kpar)
vol = Vdivide0 + dVdivide*(2*U-1)
end function

!-----------------------------------------------------------------------------------------
! Generate a random value for CFSE from a distribution with mean = average
! In the simplest case we can allow a uniform distribution about the average.
! Multiplying factor in the range (1-a, 1+a)
! Better to make it a Gaussian distribution: 
!  = average*(1+s*R)
! where R = N(0,1), s = std deviation
!-----------------------------------------------------------------------------------------
real(REAL_KIND) function generate_CFSE(average)
real(REAL_KIND) :: average, std
integer :: kpar = 0
real(REAL_KIND) :: R

! Uniform distribution
!R = par_uni(kpar)
!generate_CFSE = (1 - a + 2*a*R)*average
! Gaussian distribution
R = par_rnor(kpar)	! N(0,1)
generate_CFSE = (1 + CFSE_std*R)*average
end function

!-----------------------------------------------------------------------------------------
! Cap volume: http://mathworld.wolfram.com/SphericalCap.html
! http://en.wikipedia.org/wiki/Cubic_function
! We find x1, the first real root
! dc = centre separation distance (um)
! V = constant volume (um3)
!-----------------------------------------------------------------------------------------
subroutine cubic_solver(dc,V,R)
real(REAL_KIND) :: dc, V, R
real(REAL_KIND) :: a, b, c, d
real(REAL_KIND) :: DD0, DD1, CC, u1, x1

a = 2.0
b = 1.5*dc
c = 0
d = -(dc**3/8 + 3*V/(2*PI))

DD0 = b**2 - 3*a*c
DD1 = 2*b**3 - 9*a*b*c + 27*a**2*d
CC = (DD1 + sqrt(DD1**2 - 4*DD0**3))/2
if (CC < 0) then
	CC = -(-CC)**(1.d0/3)
else
	CC = CC**(1.d0/3)
endif

u1 = 1
R = -(1/(3*a))*(b + u1*CC + DD0/(u1*CC))

!write(*,*) 'R: ',R
!write(*,*) 'cubic: ',a*R**3 + b*R**2 + c*R + d
end subroutine

end module
