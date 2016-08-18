! To determine distribution of colony size

module colony

use global
use cellstate
implicit none

!integer, parameter :: n_colony_days=10
integer :: nmax

contains

!---------------------------------------------------------------------------------------------------
! Simulate fate of cells grown with no nutrient constraints.
! The only determinants of the colony size for a cell are (considering radiation only):
! volume
! divide_time_mean(ityp)
! radiation_tag
! p_rad_death
! growth_delay
! G2_M
! dt_delay
! t_growth_delay_end
! N_delayed_cycles_left
! The new method simply continues the simulation from where it ended, for 10 days.
!---------------------------------------------------------------------------------------------------
subroutine make_colony_distribution(n_colony_days,dist,ddist,ndist) bind(C)
!DEC$ ATTRIBUTES DLLEXPORT :: make_colony_distribution
use, intrinsic :: iso_c_binding
real(c_double) :: n_colony_days, dist(*), ddist
integer(c_int) :: ndist
real(REAL_KIND) :: V0, dVdt, dt, t, tend
real(REAL_KIND) :: tnow_save
integer :: kcell, ityp, n, idist, ncycmax, ntot, nlist_save
type (cell_type), pointer :: cp

write(logmsg,*) 'make_colony_distribution: nlist, ndist: ',nlist,ndist
call logger(logmsg)
colony_simulation = .true.
nlist_save = nlist
tnow_save = tnow
ncycmax = 24*3600*n_colony_days/divide_time_mean(1) + 3
nmax = 2**ncycmax
allocate(ccell_list(nmax))
dist(1:ndist) = 0
ntot = 0
tend = tnow + n_colony_days*24*3600    ! plate for 10 days
do kcell = 1, nlist_save
	cp => cell_list(kcell)
	if (cp%state == DEAD) cycle
	ityp = cp%celltype
	
	! Now simulate colony growth from a single cell
	tnow = tnow_save
	call make_colony(kcell,tend,n)
	if (mod(kcell,100) == 0) then
	    write(logmsg,*) 'cell: n: ',kcell,n
	    call logger(logmsg)
	endif
	ntot = ntot + n
	idist = n/ddist + 1
	dist(idist) = dist(idist) + 1
enddo 
dist(1:ndist) = dist(1:ndist)/sum(dist(1:ndist))
!write(logmsg,'(a,2i8,f8.1)') 'Colony size distribution: ', nlist_save,ntot,real(ntot)/nlist_save
!call logger(logmsg)
write(nflog,'(a,2i8,f8.1)') 'Colony size distribution: ', nlist_save,ntot,real(ntot)/nlist_save
do idist = 1,ndist
!	write(logmsg,'(i4,a,i4,f6.3)') int((idist-1)*ddist),'-',int(idist*ddist),dist(idist)
!    call logger(logmsg)
	write(nflog,'(i4,a,i4,f6.3)') int((idist-1)*ddist),'-',int(idist*ddist),dist(idist)
enddo
deallocate(ccell_list)
colony_simulation = .false.
nlist = nlist_save
tnow = tnow_save
end subroutine

!---------------------------------------------------------------------------------------------------
! The cell is at the point of division - possibly G2_M (arrested at G2/M checkpoint)
! For now only radiation tagging is handled
! Growth rate dVdt (mean) is used only to estimate the time of next division
!---------------------------------------------------------------------------------------------------
subroutine make_colony(kcell,tend,n)
integer :: kcell, n
real(REAL_KIND) :: tend, dt 
integer :: icell, ityp, nlist0, kpar=0
real(REAL_KIND) :: V0, Tdiv0, r_mean, c_rate, dVdt, Tmean, R
logical :: changed, ok
type (cell_type), pointer :: cp

!write(*,'(a,i6,2f8.0)') 'make_colony: ',kcell,tnow,tend
ccell_list(1) = cell_list(kcell)
ccell_list(1)%anoxia_tag = .false.
ccell_list(1)%aglucosia_tag = .false.
ccell_list(1)%drug_tag = .false.
ityp = ccell_list(1)%celltype
!ccell_list(1)%divide_volume = get_divide_volume(ityp,V0,Tdiv0)
dt = DELTA_T
ngaps = 0
!write(*,*) 'dVdt: ',ccell_list(1)%dVdt
nlist = 1
do while (tnow < tend)
	tnow = tnow + dt
    call colony_grower(dt,changed,ok)
!    write(*,*) 'tnow, nlist: ',tnow,nlist
enddo
n = 0
do icell = 1,nlist
	if (ccell_list(icell)%state /= DEAD) n = n+1
enddo
end subroutine

!-----------------------------------------------------------------------------------------
subroutine colony_grower(dt, changed, ok)
real(REAL_KIND) :: dt
logical :: changed, ok
integer :: k, kcell, nlist0, ityp, idrug, prev_phase, kpar=0
type(cell_type), pointer :: cp
type(cycle_parameters_type), pointer :: ccp
real(REAL_KIND) :: R
integer :: ndivide, divide_list(1000)
logical :: drugkilled
logical :: mitosis_entry, in_mitosis, divide

ok = .true.
changed = .false.
ccp => cc_parameters
nlist0 = nlist
ndivide = 0
!tnow = istep*DELTA_T !+ t_fmover
!if (colony_simulation) write(*,*) 'grower: ',nlist0,use_volume_method,tnow
do kcell = 1,nlist0
	kcell_now = kcell
	if (colony_simulation) then
	    cp => ccell_list(kcell)
	else
    	cp => cell_list(kcell)
    endif
	if (cp%state == DEAD) cycle
	ityp = cp%celltype
	divide = .false.
	mitosis_entry = .false.
	in_mitosis = .false.
	if (use_volume_method) then
!        if (colony_simulation) then
!            write(*,'(a,i6,L2,2e12.3)') 'kcell: ',kcell,cp%Iphase,cp%V,cp%divide_volume
!        endif
	    if (cp%Iphase) then
		    call growcell(cp,dt)
		    if (cp%V > cp%divide_volume) then	! time to enter mitosis
    	        mitosis_entry = .true.
	        endif
	    else
	        in_mitosis = .true.
	    endif
	else
	    prev_phase = cp%phase
!		write(*,'(a,i6,i2,4e12.3)') 'kcell, phase (a): ',kcell,cp%phase,tnow,cp%G2_time,cp%V,cp%dVdt
        call timestep(cp, ccp, dt)
        if (cp%phase >= M_phase) then
            if (prev_phase == Checkpoint2) then
                mitosis_entry = .true.
            else
                in_mitosis = .true.
            endif
        endif
		if (cp%phase < Checkpoint2 .and. cp%phase /= Checkpoint1) then
		    call growcell(cp,dt)
		endif	
	endif
	if (mitosis_entry) then
!	    write(*,*) 'mitosis_entry: ',tnow
		drugkilled = .false.
		do idrug = 1,ndrugs_used
			if (cp%drug_tag(idrug)) then
				call CellDies(kcell,cp)
				changed = .true.
				Ndrug_dead(idrug,ityp) = Ndrug_dead(idrug,ityp) + 1
				drugkilled = .true.
				exit
			endif
		enddo
		if (drugkilled) cycle
			
		if (use_volume_method) then
			if (cp%growth_delay) then
				if (cp%G2_M) then
					if (tnow > cp%t_growth_delay_end) then
						cp%G2_M = .false.
					else
						cycle
					endif
				else
					cp%t_growth_delay_end = tnow + cp%dt_delay
					cp%G2_M = .true.
					cycle
				endif
			endif
			! try moving death prob test to here
			if (cp%radiation_tag) then
				R = par_uni(kpar)
				if (R < cp%p_rad_death) then
					call CellDies(kcell,cp)
					changed = .true.
					Nradiation_dead(ityp) = Nradiation_dead(ityp) + 1
					cycle
				endif
			endif		
		else
		    ! Check for cell death by radiation lesions
		    ! For simplicity: 
		    !   cell death occurs only at mitosis entry
		    !   remaining L1 lesions and L2c misrepair (non-reciprocal translocation) are treated the same way
		    !   L2a and L2b are treated as non-fatal
		    if (cp%NL1 > 0 .or. cp%NL2(2) > 0) then
				call CellDies(kcell,cp)
				changed = .true.
				Nradiation_dead(ityp) = Nradiation_dead(ityp) + 1
				cycle
			endif		        
		endif
		
		cp%Iphase = .false.
		cp%mitosis = 0
		cp%t_start_mitosis = tnow
!		ncells_mphase = ncells_mphase + 1
	elseif (in_mitosis) then
!	    write(*,*) 'in_mitosis: ',tnow
		cp%mitosis = (tnow - cp%t_start_mitosis)/mitosis_duration
        if (cp%mitosis >= 1) then
			divide = .true.
		endif
	endif
	if (divide) then
		ndivide = ndivide + 1
		divide_list(ndivide) = kcell
	endif
enddo
do k = 1,ndivide
	changed = .true.
	kcell = divide_list(k)
	if (colony_simulation) then
	    cp => ccell_list(kcell)
	else
    	cp => cell_list(kcell)
    endif
	call colony_divider(kcell, ok)
	if (.not.ok) return
enddo
end subroutine

!-----------------------------------------------------------------------------------------
subroutine colony_divider(kcell1, ok)
integer :: kcell1
logical :: ok
integer :: kcell2, ityp, nbrs0
real(REAL_KIND) :: r(3), c(3), cfse0, cfse2, V0, Tdiv
type(cell_type), pointer :: cp1, cp2

!write(*,*) 'divider:'
!write(logmsg,*) 'divider: ',kcell1 
!call logger(logmsg)
ok = .true.
!tnow = istep*DELTA_T
if (colony_simulation) then
    cp1 => ccell_list(kcell1)
else
	cp1 => cell_list(kcell1)
endif
if (ngaps > 0) then
    kcell2 = gaplist(ngaps)
    ngaps = ngaps - 1
else
	nlist = nlist + 1
	if (nlist > nmax) then
		ok = .false.
		call logger('Error: Maximum number of cells nmax has been exceeded.  Increase and rebuild.')
		return
	endif
	kcell2 = nlist
endif
ncells = ncells + 1
if (colony_simulation) then
    cp2 => ccell_list(kcell2)
else
	cp2 => cell_list(kcell2)
endif
ityp = cp1%celltype
ncells_type(ityp) = ncells_type(ityp) + 1
ncells_mphase = ncells_mphase - 1
!cp2 => cell_list(kcell2)

cp1%state = ALIVE
cp1%generation = cp1%generation + 1
V0 = cp1%V/2
cp1%V = V0
cp1%birthtime = tnow
cp1%divide_volume = get_divide_volume(ityp,V0,Tdiv)
cp1%divide_time = Tdiv
cp1%mitosis = 0
cfse0 = cp1%CFSE
cp1%CFSE = generate_CFSE(cfse0/2)
cfse2 = cfse0 - cp1%CFSE

cp1%drug_tag = .false.
cp1%anoxia_tag = .false.
cp1%t_anoxia = 0
cp1%aglucosia_tag = .false.
cp1%t_aglucosia = 0

if (cp1%growth_delay) then
	cp1%N_delayed_cycles_left = cp1%N_delayed_cycles_left - 1
	cp1%growth_delay = (cp1%N_delayed_cycles_left > 0)
endif
cp1%G2_M = .false.
cp1%Iphase = .true.
cp1%phase = G1_phase

ndoublings = ndoublings + 1
doubling_time_sum = doubling_time_sum + tnow - cp1%t_divide_last
cp1%t_divide_last = tnow

! Second cell
!cell_list(kcell2) = cell_list(kcell1)
cp2 = cp1

! These are the variations from cp1
cp2%divide_volume = get_divide_volume(ityp,V0,Tdiv)
cp2%divide_time = Tdiv
cp2%CFSE = cfse2
if (cp2%radiation_tag) then
	Nradiation_tag(ityp) = Nradiation_tag(ityp) + 1
endif
!if (colony_simulation) write(*,'(a,i6,2e12.3)') 'new cell: ',kcell2,cp2%V,cp2%divide_volume
end subroutine

end module
