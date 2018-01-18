! Cell cycle from Basse2002

module cycle_mod
use real_kind_mod
use global

implicit none

integer, parameter :: G1_phase    = 1
integer, parameter :: Checkpoint1 = 2
integer, parameter :: S_phase     = 3
integer, parameter :: G2_phase    = 4
integer, parameter :: Checkpoint2 = 5
integer, parameter :: M_phase     = 6
integer, parameter :: dividing    = 7

logical :: use_volume_based_transition = .false.
real(REAL_KIND) :: starvation_arrest_threshold = 5
real(REAL_KIND) :: max_arrest_time = 6*3600
!integer, parameter :: NTCP = 200

contains

!--------------------------------------------------------------------------
! Damage from radiation dose following Curtis1986
! dose is the Gy received, tmin is the duration (min) over which it is delivered.
! Krepair, Kmisrepair are Curtis's epsilon_PL, epsilon_2PL
! As modified according to Bill Wilson:
! There are three classes of lesions:
!	%NL1 = number of potentially lethal (i.e. repairable) lesions = Curtis's n_PL
!	%NL2(1) = number of lethal lesions of type b
!	%NL2(2) = number of lethal lesions of type c
! During the radiation exposure, damage of different kinds occurs at a rate
! that is given by the product dose intensity (Gy/h), eta, and SER_OER.
! eta_PL and eta_L(:) are the rate constants for PL and L lesions
! SER_OER is the product of the SER and OER, where
! SER = Sensitivity Enhancement Ratio which is drug dependent
! OER = Oxygen Enhancement Ratio which is a function of intracellular O2 conc.
! (Note: OER for PL uses OER_alpha, for L uses OER_beta)
! Over the duration of the exposure repair is also occurring (this will be
! insignificant if the duration is very short).
!--------------------------------------------------------------------------
subroutine radiation_damage(cp, ccp, dose, SER_OER, tmin)
type(cell_type), pointer :: cp
type(cycle_parameters_type), pointer :: ccp
real(REAL_KIND) :: dose, SER_OER(2), tmin
real(REAL_KIND) :: dDdt, rnL(2), rnPL, drnPLdt, drnLdt(2), dtmin, dthour, fraction, Krepair, Kmissum, dr, R
integer :: nt, it, n, dn, k, ityp, kpar=0

dtmin = 0.01
dthour = dtmin/60
nt = max(tmin/dtmin, 1.0)
dDdt = dose/(nt*dthour)     ! Gy/h
rnPL = cp%NL1
rnL = cp%NL2
Kmissum = sum(ccp%Kmisrepair)
ityp = cp%celltype
if (cp%phase < S_phase) then
    Krepair = ccp%Krepair_base
elseif (cp%phase > S_phase) then
    Krepair = ccp%Krepair_max
else
    if (use_volume_based_transition) then   ! fraction = fraction of passage through S phase
        fraction = (cp%V - cp%G1_V)/(cp%S_V - cp%G1_V)
    else
        fraction = 1 - (cp%S_time - tnow)/ccp%T_S(ityp)
    endif
    Krepair = ccp%Krepair_base + fraction*(ccp%Krepair_max - ccp%Krepair_base)
endif
do it = 1,nt
    drnPLdt = ccp%eta_PL*SER_OER(1)*dDdt - Krepair*rnPL - Kmissum*rnPL**2
    rnPL = rnPL + drnPLdt*dthour
    drnLdt = ccp%eta_L*SER_OER(2)*dDdt + ccp%Kmisrepair*rnPL**2
    rnL = rnL + drnLdt*dthour
enddo
n = rnPL
dr = rnPL - n
R = par_uni(kpar)
dn = 0
if (R < dr) dn = 1
cp%NL1 = cp%NL1 + n + dn
do k = 1,2
    n = rnL(k)
    dr = rnL(k) - n
    R = par_uni(kpar)
    dn = 0
    if (R < dr) dn = 1
    cp%NL2(k) = cp%NL2(k) + n + dn
enddo
end subroutine

!--------------------------------------------------------------------------
! Phase transitions are now based on cell volume, cp%V, to allow for delay
! when growth is slowed by starvation of oxygen and/or glucose.
! Note that the volumes required for the transitions (cp%G1_V,..)  never change.
!--------------------------------------------------------------------------
subroutine timestep(cp, ccp, dt)
type(cell_type), pointer :: cp
type(cycle_parameters_type), pointer :: ccp
real(REAL_KIND) :: dt
integer :: phase, ityp, kpar=0
real(REAL_KIND) :: cf_O2, cf_glucose, pcp_O2, pcp_glucose, pcp_starvation, R
logical :: switch

phase = cp%phase
if (cp%dVdt == 0) then
	if (phase == G1_phase .or. phase == S_phase .or. phase == G2_phase) return
endif
ityp = cp%celltype
if (.not.use_metabolism) then
	if (.not.colony_simulation .and. (phase == Checkpoint1)) then    ! check for starvation arrest
		cf_O2 = cp%Cin(OXYGEN)/anoxia_threshold
		pcp_O2 = getPcp_release(cf_O2,dt)
		cf_glucose = cp%Cin(GLUCOSE)/aglucosia_threshold
		pcp_glucose = getPcp_release(cf_glucose,dt)
		pcp_starvation = pcp_O2*pcp_glucose
		if (pcp_starvation == 0) then
			return
		elseif (pcp_starvation < 1) then
			R = par_uni(kpar)
			if (R < pcp_starvation) return
		endif
	endif	
endif
if (phase == G1_phase) then
    if (use_volume_based_transition) then
        switch = (cp%V > cp%G1_V)
    else
        switch = (tnow > cp%G1_time)
    endif
    if (switch) then
        cp%phase = Checkpoint1
        cp%G1_flag = .false.
        cp%G1S_flag = .false.
        cp%G1S_time = tnow + ccp%Tcp(cp%NL1)
    endif
elseif (phase == Checkpoint1) then  ! this checkpoint combines the release from G1 delay and the G1S repair check
    if (.not.cp%G1_flag) then
        R = par_uni(kpar)
        cp%G1_flag = (R < ccp%Pk_G1(ityp)*dt)
    endif
    cp%G1S_flag = (cp%NL1 == 0 .or. tnow > cp%G1S_time)
    if (use_metabolism) then
		cp%G1S_flag = cp%G1S_flag .and. (cp%metab%A_rate > ATPg(ityp))
	endif
    if (cp%G1_flag .and. cp%G1S_flag) then
        cp%phase = S_phase
! Note: now %I_rate has been converted into equivalent %dVdt, to simplify code
!        if (use_metabolism) then
!	        cp%S_time = tnow + (cp%metab%I_rate_max/cp%metab%I_rate)*ccp%T_S(ityp)
!	    else
	        cp%S_time = tnow + (max_growthrate(ityp)/cp%dVdt)*ccp%T_S(ityp)
!	    endif
    endif
elseif (phase == S_phase) then
    if (use_volume_based_transition) then
        switch = (cp%V > cp%S_V)
    else
        switch = (tnow > cp%S_time)
    endif
    if (switch) then
        cp%phase = G2_phase
! Note: now %I_rate has been converted into equivalent %dVdt, to simplify code
!        if (use_metabolism) then
!	        cp%G2_time = tnow + (cp%metab%I_rate_max/cp%metab%I_rate)*ccp%T_G2(ityp)
!	    else
			cp%G2_time = tnow + (max_growthrate(ityp)/cp%dVdt)*ccp%T_G2(ityp)
!		endif
    endif
elseif (phase == G2_phase) then
    if (use_volume_based_transition) then
        switch = (cp%V > cp%G2_V)
    else
! Note: now %I_rate has been converted into equivalent %dVdt, to simplify code
!		if (use_metabolism) then
!			switch = (tnow > cp%G2_time .and. cp%metab%Itotal > cp%metab%I2divide) ! try this to prevent volumes decreasing 
!		else
			switch = (tnow > cp%G2_time .and. cp%V > cp%divide_volume) ! try this to prevent volumes decreasing 
!		endif
    endif
    if (switch) then
        cp%phase = Checkpoint2
        cp%G2_flag = .false.
        cp%G2M_flag = .false.
        cp%G2M_time = tnow + ccp%Tcp(cp%NL1)
    endif
elseif (phase == Checkpoint2) then ! this checkpoint combines the release from G2 delay and the G2M repair check
    if (.not.cp%G2_flag) then
        R = par_uni(kpar)
        cp%G2_flag = (R < ccp%Pk_G2(ityp)*dt)
    endif
    cp%G2M_flag = (cp%NL1 == 0 .or. tnow > cp%G2M_time)
    if (use_metabolism) then
		cp%G2M_flag = cp%G2M_flag .and. (cp%metab%A_rate > ATPg(ityp))
	endif
    if (cp%G2_flag .and. cp%G2M_flag) then
        cp%phase = M_phase
        cp%M_time = tnow + ccp%T_M (ityp)   
    endif
elseif (phase == M_phase) then
    if (tnow > cp%M_time) then
        cp%phase = dividing
!        cp%doubling_time = tnow
    endif
endif    
if (cp%NL1 > 0) then
    call radiation_repair(cp, ccp, dt)
endif
end subroutine

!--------------------------------------------------------------------------
! We want:
!    Prel = 1 for Fc >= Fmax
!    Prel = 0 for Fc <= 1
!    Prel to vary between 0 and 1 as Fc varies between 1 and Fmax
! For intermediate Fc values, make the mean duration of arrest an
! exponentially-distributed random variable with mean Tarr, then in a time
! step of length dt the probability of release is Prel = dt/Tarr.
! For Prel to vary between 0 and 1 we would need Tarr to vary between 
! infinity and dt.  Since infinity is not an option, instead we can specify
! a maximum arrest time, e.g. Tarr = Tmax (= 10h, say).  This means that as 
! Fc approaches 1, the minimum value of Prel = dt/Tmax, which for dt = 10 min 
! and Tmax = 10h gives Prel = 1/60.
!--------------------------------------------------------------------------
function getPcp_release(fc,dt) result(pcp)
real(REAL_KIND) :: fc, dt, pcp
real(REAL_KIND) :: mu

if (fc <= 1) then
    pcp = 0
elseif (fc >= starvation_arrest_threshold) then
    pcp = 1
else
    mu = dt + (max_arrest_time - dt)*(starvation_arrest_threshold - fc)/(starvation_arrest_threshold - 1)
    pcp = dt/mu
endif
end function

!--------------------------------------------------------------------------
! Time unit = hour
! This may need to be changed, because it implicitly assumes that no more
! than one repair and one misrepair of each type can occur within a time step.
! The fix would be to subdivide the time step, as in the damage subroutine.
!--------------------------------------------------------------------------
subroutine radiation_repair(cp, ccp, dt)
type(cell_type), pointer :: cp
type(cycle_parameters_type), pointer :: ccp
real(REAL_KIND) :: dt
integer :: i, ityp, kpar=0
real(REAL_KIND) :: dthour, fraction, Krepair, R

dthour = dt/3600    ! sec -> hour
ityp = cp%celltype
if (cp%phase < S_phase) then
    Krepair = ccp%Krepair_base
elseif (cp%phase > S_phase) then
    Krepair = ccp%Krepair_max
else
    if (use_volume_based_transition) then   ! fraction = fraction of passage through S phase
        fraction = (cp%V - cp%G1_V)/(cp%S_V - cp%G1_V)
    else
        fraction = 1 - (cp%S_time - tnow)/ccp%T_S(ityp)
    endif
    Krepair = ccp%Krepair_base + fraction*(ccp%Krepair_max - ccp%Krepair_base)
endif
R = par_uni(kpar)
if (R < cp%NL1*Krepair*dthour) then
    cp%NL1 = cp%NL1 - 1
    if (cp%NL1 == 0) return
endif
do i = 1,2
    R = par_uni(kpar)
    if (R < cp%NL1**2*ccp%Kmisrepair(i)*dthour) then
        cp%NL1 = cp%NL1 - 1
        cp%NL2(i) = cp%NL2(i) + 1
        if (cp%NL1 == 0) return
    endif
enddo
end subroutine

end module

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!program main
!use cycle_mod
!use Tcp_mod
!implicit none
!integer :: istep, icell, nd
!real(REAL_KIND) :: tnow, Td
!real(REAL_KIND) :: DELTA_T, hours, dose, tmin
!integer :: Ncells, Nsteps, i
!type(cell_type), allocatable, target :: cell_list(:)
!type(cell_type), pointer :: cp
!type(cycle_parameters_type), target :: cc_parameters
!type(cycle_parameters_type), pointer :: ccp
!logical :: all_divided
!
!ccp => cc_parameters
!ccp%eta_L = 0.14                    ! /Gy
!ccp%eta_PL = 0.6                    ! /Gy
!ccp%T_G1 = 12                       ! h
!ccp%T_S = 9                         ! h
!ccp%T_G2 = 1                        ! h
!ccp%T_M = 0.5                       ! h
!ccp%G1_mean_delay = 3               ! h
!ccp%G2_mean_delay = 2               ! h
!ccp%Pk_G1 = 1./ccp%G1_mean_delay    ! /h
!ccp%Pk_G2 = 1./ccp%G2_mean_delay    ! /h
!ccp%Krepair = 0.1
!ccp%Kmisrepair(1:3) = 0.01852 = (1./3)*(0.5/9)
!ccp%Kcp = 0.13
!
!DELTA_T = 0.05   ! hours
!Ncells = 100
!allocate(cell_list(Ncells))
!
!! test radiation dose
!cp => cell_list(1)
!dose = 5
!tmin = 5
!do i = 1,10
!    cp%NL1 = 0
!    cp%NL2 = 0
!    call radiation_damage(cp, ccp, dose, tmin)
!    write(*,*) 'radiation damage: ',cp%NL1,cp%NL2
!enddo
!stop
!
!write(*,*) 'makeTCP'
!call makeTCP(ccp%tcp,NTCP,ccp%epsilon_PL,ccp%epsilon_2PL,ccp%Kcp) ! set checkpoint repair time limits 
!write(*,*) 'did makeTCP'
!
!hours = 10*24
!Nsteps = hours/DELTA_T
!tnow = 0
!do icell = 1,Ncells
!    cp => cell_list(icell)
!    cp%phase = G1_phase
!    cp%G1_flag = .false.
!    cp%G1_time = tnow + ccp%T_G1
!    cp%NL1 = 100
!    cp%NL2 = 0
!enddo
!do istep = 1,Nsteps
!    all_divided = .true.
!    do icell = 1,Ncells
!        cp => cell_list(icell)
!        if (cp%phase == divided) cycle
!        all_divided = .false.
!        tnow = istep*DELTA_T
!        call timestep(cp, ccp, tnow, DELTA_T)
!    enddo
!    if (all_divided) then
!        write(*,*) 'All cells have divided'
!        exit
!    endif
!enddo
!do icell = 1,Ncells
!    cp => cell_list(icell)
!    write(*,'(a,i6,f6.2,5i6)') 'Cell: ',icell,cp%doubling_time, cp%phase, cp%NL1, cp%NL2
!enddo
!stop
!
!Td = 0
!nd = 0
!do icell = 1,Ncells
!    if (cell_list(icell)%phase == divided) then
!        nd = nd+1
!        Td = Td + cell_list(icell)%doubling_time
!    endif
!enddo
!write(*,'(a,i6,f6.2)') 'Average doubling time: ',nd,Td/nd
!end program

