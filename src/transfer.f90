! Transferring results to the GUI

module transfer

use global
use chemokine
use metabolism
use cycle_mod
use envelope

implicit none

real(REAL_KIND), allocatable :: SFlookup(:)
real(REAL_KIND) :: dp_survive

contains

!-----------------------------------------------------------------------------------------
! Removing references to anoxia, aglucosia
! Adding non-viable, nogrow, ndead
! Note: different variable naming
!    Vmonolayer                                        Vspheroid
!    ----------                                        ---------
!    EC (Caverage at medium bottom, at cell layer)
!    caverage (IC concentration - all the same)        cmedium (average in the medium outside the blob)
!    cmedium (Cmediumave, average medium)              cbdry (average at blob boundary)
!-----------------------------------------------------------------------------------------


!--------------------------------------------------------------------------------
! Need to compare growth rate with a fraction of average growth rate
! Should ngrowth(i) include cells at checkpoints?
! nogrow is the count of viable cells that are not growing
!--------------------------------------------------------------------------------
subroutine getGrowthCount(ngrowth, nogrow, nphase,nmutations,nclono)
integer :: ngrowth(3), nogrow(:), nphase(:), nmutations, nclono
integer :: kcell, i, ityp, iphase, nch1, nch2, ip
real(REAL_KIND) :: r_mean(2), ps1, ps2, ps, rclono, alfa, SF
real(REAL_KIND) :: N50 = 5.64	! number of divisions to get 50 cells
type(cell_type), pointer :: cp
type(cycle_parameters_type), pointer :: ccp

if (use_cell_cycle) then
    r_mean = max_growthrate
else
    r_mean = Vdivide0/(2*divide_time_mean)
endif
ngrowth = 0
nogrow = 0
nphase = 0
nmutations = 0
rclono = 0
do kcell = 1,nlist
	cp => cell_list(kcell)
	ityp = cp%celltype
	if (cp%state == DEAD) cycle
	if (.not.(cp%state == DYING) .and. cp%metab%A_rate < ATPg) then
		nogrow(ityp) = nogrow(ityp) + 1
	endif
	do i = 1,3
		if (cp%dVdt < growthcutoff(i)*r_mean(ityp) .or. &
		    cp%phase == G1_checkpoint .or. cp%phase == S_checkpoint .or. cp%phase == G2_checkpoint) then
		    ngrowth(i) = ngrowth(i) + 1
!		    write(*,'(a,3i6,3e12.3)') 'getGrowthCount: ',kcell,ityp,i,cell_list(kcell)%dVdt,growthcutoff(i),r_mean(ityp)
		endif
	enddo
	iphase = min(cp%phase,M_phase)
	nphase(iphase) = nphase(iphase) + 1
	if (iphase == S_phase) then
	    if (.not.cp%arrested) nphase(7) = nphase(7) + 1
	endif
	nmutations = nmutations + cp%N_Ch1 + cp%N_Ch2
	ccp => cc_parameters(ityp)
	if (.not.cp%state == DYING) then
		nch1 = cp%N_Ch1
		if (nch1 > 0) then
			ps1 = ccp%psurvive_Ch1**nch1
		else
			ps1 = 1
		endif
		nch2 = cp%N_Ch2
		if (nch2 > 0) then
			ps2 = ccp%psurvive_Ch2**nch2
		else
			ps2 = 1
		endif
		ps = ps1*ps2
		if (ps == 1) then
			rclono = rclono + 1
		else
!			rclono = rclono + ps**N50
            ip = ps/dp_survive
            alfa = ps - ip*dp_survive
            SF = (1-alfa)*SFlookup(ip) + alfa*SFlookup(ip+1)
		endif
	endif
enddo
nclono = rclono + 0.5
end subroutine

!--------------------------------------------------------------------------------
! From the probability of survival of division, estimate the probability of
! growing a colony with more than Nlim cells.
!--------------------------------------------------------------------------------
subroutine clonogenic(ityp, p_survive, Nlim, SF)
real(REAL_KIND) :: p_survive, SF
integer :: ityp, Nlim
real(REAL_KIND) :: R, p_d
integer :: nc, irun, nruns, idiv, ndiv, n, ntemp, i, kpar = 0
logical :: ended, dbug
integer :: n_colony_days = 10

p_d = 1 - p_survive
dbug = .false.
nruns = 10000
ndiv = 24*3600*n_colony_days/divide_time_mean(ityp) + 1
nc = 0
do irun = 1,nruns
	n = 1
	ended = .false.
	do idiv = 1,ndiv
		ntemp = n
		do i = 1,ntemp
!			call random_number(R)
            R = par_uni(kpar)
			if (R < p_d) then
				n = n - 1
				if (n == 0) then
					ended = .true.
					exit
				endif
			else
				n = n + 1
			endif
	        if (n >= Nlim) then
	            nc = nc + 1
	            ended = .true.
	            exit
	        endif
		enddo
	    if (dbug) then
	        write(*,*) 'irun, n: ',irun,n
	    endif
		if (ended) exit
	enddo
!	if (n == 0) then
!		kbox = 0
!	else
!		if (boxsize == 1) then
!			kbox = n
!		else
!			kbox = n/boxsize + 1
!		endif
!	endif
!	if (kbox > nboxes) then
!		write(*,*) 'Error: irun, n, kbox, nboxes: ',irun, n, kbox, nboxes
!		stop
!	endif
!	ncol(kbox) = ncol(kbox) + 1
enddo
SF = nc/real(nruns)
end subroutine

!--------------------------------------------------------------------------------
! Generate and store a lookup table for SF as a function of p (division survival prob).
!--------------------------------------------------------------------------------
subroutine GenerateSFlookup(ityp)
integer :: ityp
integer :: np_survive
integer :: Nlim = 50
integer :: i
real(REAL_KIND) :: p_survive

write(nflog,*)
write(nflog,*) 'GenerateSFlookup: ityp: ',ityp
dp_survive = 0.01
np_survive = 100
if (allocated(SFlookup)) deallocate(SFlookup)
allocate(SFlookup(0:np_survive))
do i = 0,np_survive
    p_survive = i*dp_survive
    if (i == 0) then
        SFlookup(i) = 0
    elseif (i == np_survive) then
        SFlookup(i) = 1
    else
        call clonogenic(ityp, p_survive, Nlim, SFlookup(i))
    endif
    write(nflog,'(i4,2f8.5)') i,p_survive,SFlookup(i)
enddo
end subroutine

!-----------------------------------------------------------------------------------------
! Note: the only requirement is that the array has the SAVE attribute
!-----------------------------------------------------------------------------------------
!subroutine test_array(array) bind(C)
subroutine test_array(narr1, narr2, cptr) bind(C)
!DEC$ ATTRIBUTES DLLEXPORT :: test_array
use, intrinsic :: iso_c_binding
!TYPE, BIND(C) :: PASS
!   INTEGER (C_INT) :: lenf
!   TYPE (C_PTR)    :: f
!	real(c_double) :: f(*)
!END TYPE PASS
!TYPE (PASS), INTENT(INOUT) :: array
integer(c_int) :: narr1, narr2
TYPE (C_PTR)    :: cptr
!real(c_double), allocatable, save :: a(:,:)
real(c_double), save :: a(4,3)
integer :: i1, i2
integer :: n1=4, n2=3

!allocate(a(n1,n2))
do i1 = 1,n1
	do i2 = 1,n2
		a(i1,i2) = i1 + 10*i2
	enddo
enddo

narr1 = n1
narr2 = n2
cptr = c_loc(a)

end subroutine

!-----------------------------------------------------------------------------------------
! Distances are still all cm 
! Added ixyz to pass the appropriate Caverage array index
!-----------------------------------------------------------------------------------------
subroutine get_fielddata(axis, fraction, fdata, ixyz, res) bind(C)
!DEC$ ATTRIBUTES DLLEXPORT :: get_fielddata
use, intrinsic :: iso_c_binding
integer(c_int) :: axis, ixyz, res
real(c_double) :: fraction
type (fielddata_type)    :: fdata
type(celldata_type), save :: cdata(4000)
type(cell_type), pointer :: cp
real(REAL_KIND) :: x, y, z, c(3), r, rad, csum(3), bcentre(3)
integer :: i, kcell, nc

fdata%NX = NX
fdata%NY = NY
fdata%NZ = NZ
fdata%NCONST = NCONST
fdata%DX = DELTA_X
fdata%conc_ptr = c_loc(Caverage)
fdata%cell_ptr = c_loc(cdata)

!write(nflog,*) 'get_fielddata: axis: ',axis
! Find blob centre
csum = 0
nc = 0
do kcell = 1,nlist
	cp => cell_list(kcell)
	if (cp%state == DEAD) cycle
	nc = nc+1
	csum = csum + cp%centre(:,1)
!	write(nflog,'(a,i6,3f8.4)') 'cell centre: ',kcell,cp%centre(:,1) 
enddo
bcentre = csum/nc

nc = 0
! start by fixing on central plane, ignoring fraction
if (axis == X_AXIS) then
!	x = ((NX-1)/2)*DELTA_X
	x = bcentre(1)	
	ixyz = x/DELTA_X		! approx?
elseif (axis == Y_AXIS) then
!	y = ((NY-1)/2)*DELTA_X 
	y = bcentre(2)	
	ixyz = y/DELTA_X
elseif (axis == Z_AXIS) then
!	z = ((NZ-1)/2)*DELTA_X 
	z = bcentre(3)	
	ixyz = z/DELTA_X
endif
do kcell = 1,nlist
	cp => cell_list(kcell)
	if (cp%state == DEAD) cycle
	do i = 1,cp%nspheres
		c = cp%centre(:,i)
		r = cp%radius(i)
		if (axis == X_AXIS) then
			if (c(1) + r < x .or. c(1) - r > x) cycle
			rad = sqrt(r**2-(c(1)-x)**2)
			nc = nc + 1
			cdata(nc)%radius = rad
			cdata(nc)%centre(1:2) = [c(2),(NX-1.1)*DELTA_X - c(3)]		! always use centre(1:2) to store the 2D coordinates
		elseif (axis == Y_AXIS) then
			if (c(2) + r < y .or. c(2) - r > y) cycle
			rad = sqrt(r**2-(c(2)-y)**2)
			nc = nc + 1
			cdata(nc)%radius = rad
			cdata(nc)%centre(1:2) = [c(1),(NX-1.1)*DELTA_X - c(3)]		! invert z = c(3).  1.1 corrects for cell-wall interpenetration
		elseif (axis == Z_AXIS) then
			if (c(3) + r < z .or. c(3) - r > z) cycle
			rad = sqrt(r**2-(c(3)-z)**2)
			nc = nc + 1
			cdata(nc)%radius = rad
			cdata(nc)%centre(1:2) = [c(1),c(2)]		! always use centre(1:2) to store the 2D coordinates
		endif
		! status
!		if (cp%anoxia_tag .or. cp%state == DYING) then
!			cdata(nc)%status = 2	! tagged to die of anoxia
!		elseif (cp%aglucosia_tag) then
!			cdata(nc)%status = 4	! tagged to die of aglucosia
!		elseif (cp%radiation_tag) then
!			cdata(nc)%status = 10
!			write(nflog,*) 'Tagged to die from radiation: ',kcell
!		elseif (cp%drug_tag(1)) then
!			cdata(nc)%status = 11
!		elseif (cp%drug_tag(1)) then
!			cdata(nc)%status = 12
!		elseif (cp%Cin(OXYGEN) < hypoxia_threshold) then
!			cdata(nc)%status = 1	! radiobiological hypoxia
!		elseif (cp%mitosis > 0) then
!			cdata(nc)%status = 3	! in mitosis
!		else
!			cdata(nc)%status = 0
!		endif
		cdata(nc)%status = getstatus(cp)
		cdata(nc)%celltype = cp%celltype
	enddo
enddo
fdata%ncells = nc
res = 0
end subroutine

!--------------------------------------------------------------------------------
! A growing cell (r_A > ATPg) has status 0 for celltype 1, 10 for celltype 2
! Otherwise:
!
! icolourscheme = 0
!   1 = not growing (ATPs < r_A < ATPg
!   2 = DYING
!   3 = in mitosis
!--------------------------------------------------------------------------------
function getstatus(cp) result(status)
type(cell_type), pointer :: cp
real(REAL_KIND) :: A_rate
integer :: status

A_rate = cp%metab%A_rate
if (icolourscheme == 0) then
    if (cp%phase >= M_phase) then
        status = 3
    elseif (A_rate >= ATPg) then    ! happy growing cell
	    if (cp%celltype == 1) then
		    status = 0
	    else
		    status = 10
	    endif
    elseif (cp%state == DYING) then
        status = 2
    elseif (A_rate < ATPg) then
        status = 1
    endif
elseif (icolourscheme == 1 .or. icolourscheme == 2) then
    if (cp%state == DYING) then
        if (cp%ATP_tag) then                            ! ATP death
            status = 4
        else                                            ! treatment death
            status = 5
        endif
    elseif (cp%Cin(OXYGEN) < hypoxia_threshold) then    ! hypoxic
        status = 6
    elseif (A_rate < ATPg) then                         ! not growing
        status = 3
    else                                                ! growing
        if (cp%phase == S_phase) then
            status = 1                                  ! S-phase
        elseif (cp%phase == M_phase) then
            status = 2                                  ! M-phase
        elseif (cp%celltype == 1) then
            status = 0                                  ! type 1
        else
            status = 10                                 ! type 2
        endif
    endif   
endif
end function

!--------------------------------------------------------------------------------
! A growing cell has status 0 for celltype 1, 10 for celltype 2
!--------------------------------------------------------------------------------
function old_getstatus(cp) result(status)
type(cell_type), pointer :: cp
integer :: status

!if (cp%anoxia_tag) then
!	status = 2	! tagged to die of anoxia
!elseif (cp%aglucosia_tag) then
!	status = 4	! tagged to die of aglucosia
if (cp%state == DYING) then
    status = 2
elseif (cp%radiation_tag) then
	status = 11
elseif (cp%drug_tag(1)) then
	status = 12
elseif (cp%drug_tag(2)) then
	status = 13
elseif (cp%Cin(OXYGEN) < hypoxia_threshold) then
	status = 1	! radiobiological hypoxia
!elseif (cp%mitosis > 0) then
elseif (cp%V > 0.9*cp%divide_volume) then  ! just a surrogate for mitosis
	status = 3	! in mitosis
else
	if (cp%celltype == 1) then
		status = 0
	else
		status = 10
	endif
endif
end function

!--------------------------------------------------------------------------------
! This version computes concentrations in a spherical shell of thickness dr.
!--------------------------------------------------------------------------------
subroutine WriteProfileData
integer :: ns
real(REAL_KIND), allocatable :: ex_conc(:,:)
real(REAL_KIND) :: dr = 20.e-4	! 20um -> cm
real(REAL_KIND) :: dx, xyz0(3), dxyz(3), r2, r
integer :: ntot, ir, nr, ichemo, kcell, site(3)
integer :: i, ic, nc
integer, parameter :: max_shells = 100
integer :: cnt(max_shells)
integer :: icmap(0:MAX_CHEMO+N_EXTRA)		! maps ichemo -> ic
character*(16) :: title(1:MAX_CHEMO+N_EXTRA)
character*(128) :: filename
character*(6) :: mintag
type(cell_type), pointer :: cp

!write(nflog,*) 'WriteProfileData: MAX_CHEMO,N_EXTRA: ',MAX_CHEMO,N_EXTRA
! First find the centre of the blob
dx = DELTA_X
xyz0 = 0
ntot = 0
do kcell = 1,nlist
	cp => cell_list(kcell)
	if (cp%state == DEAD) cycle
	ntot = ntot+1
	xyz0 = xyz0 + cp%centre(:,1)
enddo
if (ntot == 0) return
xyz0 = xyz0/ntot	! this is the blob centre

! Set up icmap
ic = 0
do ichemo = 0,MAX_CHEMO+N_EXTRA
	if (ichemo == CFSE) then
		if (.not.chemo(ichemo)%used) cycle
		ic = ic + 1
		title(ic) = 'CFSE'
	    icmap(ichemo) = ic
	elseif (ichemo <= MAX_CHEMO) then
		if (.not.chemo(ichemo)%used) cycle
		ic = ic + 1
		title(ic) = chemo(ichemo)%name
	    icmap(ichemo) = ic
	elseif (ichemo == GROWTH_RATE) then
		ic = ic + 1
		title(ic) = 'Growth_rate'
	    icmap(ichemo) = ic
	elseif (ichemo == CELL_VOLUME) then
		ic = ic + 1
		title(ic) = 'Cell_volume'
	    icmap(ichemo) = ic
!	elseif (ichemo == O2_BY_VOL) then
!		ic = ic + 1
!		title(ic) = 'Cell_O2xVol'
!	    icmap(ichemo) = ic
	endif
enddo
nc = ic
allocate(ex_conc(nc,max_shells))
ex_conc = 0

! Now look at shells at dr spacing. 
nr = 0
cnt = 0
do kcell = 1,nlist
	cp => cell_list(kcell)
	if (cp%state == DEAD) cycle
	dxyz = cp%centre(:,1) - xyz0
	r2 = 0
	do i = 1,3
		r2 = r2 + dxyz(i)**2
	enddo
	r = sqrt(r2)
	ir = r/dr + 1
	nr = max(ir,nr)
	cnt(ir) = cnt(ir) + 1
	do ichemo = 0,MAX_CHEMO+N_EXTRA
		ic = icmap(ichemo)
		if (ichemo == CFSE .and. chemo(ichemo)%used) then
			ex_conc(ic,ir) = ex_conc(ic,ir) + cp%cfse
		elseif (ichemo <= MAX_CHEMO .and. chemo(ichemo)%used) then
!			if (cp%conc(ichemo) > 10) then
!			write(*,'(a,3i6,e12.3)') 'bad conc: ',kcell,ichemo,ic,cp%conc(ichemo)
!			stop
!			endif
			ex_conc(ic,ir) = ex_conc(ic,ir) + cp%Cin(ichemo)
		elseif (ichemo == GROWTH_RATE) then
			ex_conc(ic,ir) = ex_conc(ic,ir) + cp%dVdt
		elseif (ichemo == CELL_VOLUME) then
			ex_conc(ic,ir) = ex_conc(ic,ir) + cp%V
		endif
	enddo
enddo

! Average
do ir = 1,nr
	do ic = 1,nc
		ex_conc(ic,ir) = ex_conc(ic,ir)/cnt(ir)
	enddo
enddo			
	
! Remove time tag from the filename for download from NeSI
write(mintag,'(i6)') int(istep*DELTA_T/60)
filename = saveprofile%filebase
filename = trim(filename)//'_'
filename = trim(filename)//trim(adjustl(mintag))
filename = trim(filename)//'min.dat'
open(nfprofile,file=filename,status='replace')
write(nfprofile,'(a,a)') 'GUI version: ',gui_run_version
write(nfprofile,'(a,a)') 'DLL version: ',dll_run_version
write(nfprofile,'(i6,a)') int(istep*DELTA_T/60),' minutes'
write(nfprofile,'(i6,a)') nr,' shells'
write(nfprofile,'(f6.2,a)') 10000*dr,' dr (um)'
write(nfprofile,'(32a16)') title(1:nc)
do ir = 1,nr
	write(nfprofile,'(32(e12.3,4x))') ex_conc(1:nc,ir)
enddo
close(nfprofile)
deallocate(ex_conc)
!write(nflog,*) 'did WriteProfileData: ',filename
end subroutine

!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------
subroutine WriteSliceData
character*(128) :: filename
character*(6) :: mintag
integer, parameter :: NSMAX = 200
integer :: iz, nslices
real(REAL_KIND) :: dzslice
real(REAL_KIND) :: rad(NSMAX,2)

! Remove time tag from the filename for download from NeSI
write(mintag,'(i6)') int(istep*DELTA_T/60)
filename = saveslice%filebase
filename = trim(filename)//'_'
filename = trim(filename)//trim(adjustl(mintag))
filename = trim(filename)//'min.dat'
open(nfslice,file=filename,status='replace')
write(nfslice,'(a,a)') 'GUI version: ',gui_run_version
write(nfslice,'(a,a)') 'DLL version: ',dll_run_version
write(nfslice,'(i6,a)') int(istep*DELTA_T/60),' minutes'
call getSlices(nslices,dzslice,NSMAX,rad)
write(nfslice,'(e12.3,a)') dzslice,'  slice thickness (cm)'
write(nfslice,'(i4,a)') nslices,'  number of slices'
do iz = 1,nslices
    write(nfslice,'(i4,3e12.3)') iz,rad(iz,:),rad(iz,2)-rad(iz,1)
enddo
close(nfslice)
end subroutine

!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------
subroutine WriteFACSData
character*(128) :: filename
character*(6) :: mintag
!integer :: ic, nc, ichem(2)
type(cell_type), pointer :: cp
real(REAL_KIND) :: hour, val, facs_data(32)
integer :: k, kcell, iextra, ichemo, ivar, nvars, var_index(32)
character*(24) :: var_name(32)

hour = istep*DELTA_T/3600.
nvars = 1	! CFSE
var_index(nvars) = 0
var_name(nvars) = 'CFSE'
do ichemo = 1,MAX_CHEMO
	if (.not.chemo(ichemo)%used) cycle
	nvars = nvars + 1
	var_index(nvars) = ichemo
	var_name(nvars) = chemo(ichemo)%name
enddo
do iextra = 1,N_EXTRA-2
	nvars = nvars + 1
	var_index(nvars) = MAX_CHEMO + iextra
	if (iextra == 1) then
		var_name(nvars) = 'GROWTH_RATE'
	elseif (iextra == 2) then
		var_name(nvars) = 'CELL_VOLUME'
!	elseif (iextra == 3) then
!		var_name(nvars) = 'O2_BY_VOL'
	endif
enddo

!nc = 0
!if (chemo(DRUG_A)%present) then
!	nc = nc + 1
!	ichem(nc) = DRUG_A
!endif
!if (chemo(DRUG_B)%present) then
!	nc = nc + 1
!	ichem(nc) = DRUG_B
!endif
! Remove time tag from the filename for download from NeSI
write(mintag,'(i6)') int(istep*DELTA_T/60)
filename = saveFACS%filebase
filename = trim(filename)//'_'
filename = trim(filename)//trim(adjustl(mintag))
filename = trim(filename)//'min.dat'
open(nfFACS,file=filename,status='replace')
write(nfFACS,'(a,a,2a12,a,i8,a,f8.2)') trim(header),' ',gui_run_version, dll_run_version, ' step: ',istep, ' hour: ',hour
do ivar = 1,nvars
	write(nfFACS,'(a14,$)') trim(var_name(ivar))
enddo
write(nfFACS,*)
do kcell = 1,nlist
	cp => cell_list(kcell)
	if (cp%state == DEAD) cycle
!	write(nfFACS,'(2e14.5)') (cp%conc(ichem(ic)),ic=1,nc)
	k = 0
	do ivar = 1,nvars
		ichemo = var_index(ivar)
		if (ichemo == 0) then
			val = cell_list(kcell)%CFSE
		elseif (ichemo <= MAX_CHEMO) then
			val = cell_list(kcell)%Cin(ichemo)
		elseif (ichemo == GROWTH_RATE) then
			val = cell_list(kcell)%dVdt
		elseif (ichemo == CELL_VOLUME) then
			val = cell_list(kcell)%V
!		elseif (ichemo == O2_BY_VOL) then
!			val = cell_list(kcell)%volume*cell_list(kcell)%conc(OXYGEN)
		endif
		k = k+1
		facs_data(k) = val
	enddo
	write(nfFACS,'(32e14.5)') facs_data(1:nvars)
enddo
close(nfFACS)
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine get_constituents(nvars,cvar_index,nvarlen,name_array,narraylen) BIND(C)
!DEC$ ATTRIBUTES DLLEXPORT :: get_constituents
use, intrinsic :: iso_c_binding
character(c_char) :: name_array(0:*)
integer(c_int) :: nvars, cvar_index(0:*), nvarlen, narraylen
integer :: ivar, k, ichemo
character*(24) :: name
character(c_char) :: c

write(nflog,*) 'get_cell_constituents'
nvarlen = 24
ivar = 0
k = ivar*nvarlen
cvar_index(ivar) = 0	! CFSE
name = 'CFSE'
call copyname(name,name_array(k),nvarlen)
do ichemo = 1,MAX_CHEMO
	if (.not.chemo(ichemo)%used) cycle
	ivar = ivar + 1
	k = ivar*nvarlen
	cvar_index(ivar) = ichemo
	name = chemo(ichemo)%name
	write(nflog,*) 'get_cell_constituents: ',ichemo,name
	call copyname(name,name_array(k),nvarlen)
enddo
ivar = ivar + 1
k = ivar*nvarlen
cvar_index(ivar) = GROWTH_RATE
name = 'Growth rate'
call copyname(name,name_array(k),nvarlen)
ivar = ivar + 1
k = ivar*nvarlen
cvar_index(ivar) = CELL_VOLUME
name = 'Cell volume'
call copyname(name,name_array(k),nvarlen)
ivar = ivar + 1
k = ivar*nvarlen
cvar_index(ivar) = O2_BY_VOL
name = 'Cell O2xVol'
call copyname(name,name_array(k),nvarlen)
ivar = ivar + 1
k = ivar*nvarlen
cvar_index(ivar) = CYCLE_PHASE
name = 'Cycle phase'
call copyname(name,name_array(k),nvarlen)
nvars = ivar + 1
write(nflog,*) 'did get_cell_constituents'
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine copyname(name,name_array,n)
character*(*) :: name
character :: name_array(*)
integer :: n
integer :: k

do k = 1,n
	name_array(k) = name(k:k)
enddo
end subroutine

!-----------------------------------------------------------------------------------------
! Get number of live cells
!-----------------------------------------------------------------------------------------
subroutine get_nFACS(n) BIND(C)
!DEC$ ATTRIBUTES DLLEXPORT :: get_nfacs
use, intrinsic :: iso_c_binding
integer(c_int) :: n
integer :: k, kcell

n = 0
do kcell = 1,nlist
	if (cell_list(kcell)%state == DEAD) cycle
	n = n+1
enddo
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine get_FACS(facs_data, vmin, vmax, scale_drug_by_vol) BIND(C)
!DEC$ ATTRIBUTES DLLEXPORT :: get_facs
use, intrinsic :: iso_c_binding
real(c_double) :: val, facs_data(*), vmin(*), vmax(*)
logical(c_bool),value :: scale_drug_by_vol
integer :: k, kcell, iextra, ichemo, ivar, nvars, var_index(32)
real(REAL_KIND) :: cfse_min, Vcell_pL

write(nflog,*) 'get_FACS'
nvars = 1	! CFSE
var_index(nvars) = 0
do ichemo = 1,MAX_CHEMO
	if (.not.chemo(ichemo)%used) cycle
	nvars = nvars + 1
	var_index(nvars) = ichemo
enddo
do iextra = 1,N_EXTRA-1
	nvars = nvars + 1
	var_index(nvars) = MAX_CHEMO + iextra
enddo
vmin(1:MAX_CHEMO+N_EXTRA) = 1.0e10
vmax(1:MAX_CHEMO+N_EXTRA) = -1.0e10
cfse_min = 1.0e20
k = 0
do kcell = 1,nlist
	if (cell_list(kcell)%state == DEAD) cycle
	Vcell_pL = cm3_pL*cell_list(kcell)%V
	do ivar = 1,nvars
		ichemo = var_index(ivar)
		if (ichemo == 0) then
			val = cell_list(kcell)%CFSE
			cfse_min = min(val,cfse_min)
		elseif (ichemo <= MAX_CHEMO) then
			val = cell_list(kcell)%Cin(ichemo)
			if (ichemo >= DRUG_A) then
			    val = val*Vcell_pL
			endif
		elseif (ichemo == GROWTH_RATE) then
			val = cm3_pL*cell_list(kcell)%dVdt
		elseif (ichemo == CELL_VOLUME) then
			val = Vcell_pL
		elseif (ichemo == O2_BY_VOL) then
			val = Vcell_pL*cell_list(kcell)%Cin(OXYGEN)
		elseif (ichemo == CYCLE_PHASE) then
			val = cell_list(kcell)%phase
		endif
		k = k+1
		facs_data(k) = val
		vmin(ichemo+1) = min(val,vmin(ichemo+1))
		vmax(ichemo+1) = max(val,vmax(ichemo+1))
	enddo
enddo
end subroutine

!-----------------------------------------------------------------------------------------
! nhisto is the number of histogram boxes
! vmax(ivar) is the maximum value for variable ivar
! Probably need to adjust vmax() to a roundish value
!
! Compute 3 distributions: 1 = both cell types
!                          2 = type 1
!                          3 = type 2
! Stack three cases in vmax() and histo_data()
!-----------------------------------------------------------------------------------------
subroutine get_histo(nhisto, histo_data, vmin, vmax, histo_data_log, vmin_log, vmax_log, scale_drug_by_vol) BIND(C)
!DEC$ ATTRIBUTES DLLEXPORT :: get_histo
use, intrinsic :: iso_c_binding
integer(c_int),value :: nhisto
real(c_double) :: vmin(*), vmax(*), histo_data(*)
real(c_double) :: vmin_log(*), vmax_log(*), histo_data_log(*)
logical(c_bool),value :: scale_drug_by_vol
real(REAL_KIND) :: val, val_log, Vcell_pL
integer :: n(3), i, ih, k, kcell, ict, ichemo, ivar, nvars, var_index(32)
integer,allocatable :: cnt(:,:,:)
real(REAL_KIND),allocatable :: dv(:,:), valmin(:,:), valmax(:,:)
integer,allocatable :: cnt_log(:,:,:)
integer :: phase_cnt(3,7), ivar_phase
real(REAL_KIND),allocatable :: dv_log(:,:), valmin_log(:,:), valmax_log(:,:)
type(cell_type), pointer :: cp
!real(REAL_KIND) :: vmin_log(100), vmax_log(100)
!real(REAL_KIND),allocatable :: histo_data_log(:)

write(nflog,*) 'get_histo: histogram using volume scaling for drugs: ',scale_drug_by_vol
nvars = 1	! CFSE
var_index(nvars) = 0
do ichemo = 1,MAX_CHEMO
	if (.not.chemo(ichemo)%used) cycle
	nvars = nvars + 1
	var_index(nvars) = ichemo
enddo
nvars = nvars + 1
var_index(nvars) = GROWTH_RATE
nvars = nvars + 1
var_index(nvars) = CELL_VOLUME
nvars = nvars + 1
var_index(nvars) = O2_BY_VOL
nvars = nvars + 1
var_index(nvars) = CYCLE_PHASE
ivar_phase = nvars

allocate(cnt(3,nvars,nhisto))
allocate(dv(3,nvars))
allocate(valmin(3,nvars))
allocate(valmax(3,nvars))
allocate(cnt_log(3,nvars,nhisto))
allocate(dv_log(3,nvars))
allocate(valmin_log(3,nvars))
allocate(valmax_log(3,nvars))
!allocate(histo_data_log(10000))
cnt = 0
valmin = 1.0e10
valmax = -1.0e10
cnt_log = 0
valmin_log = 1.0e10
valmax_log = -1.0e10
n = 0
do kcell = 1,nlist
	cp => cell_list(kcell)
	if (cp%state == DEAD) cycle
	ict = cp%celltype
	Vcell_pL = cm3_pL*cp%V
	do ivar = 1,nvars
		ichemo = var_index(ivar)
		if (ichemo == 0) then
			val = cp%CFSE
		elseif (ichemo <= MAX_CHEMO) then
			val = cp%Cin(ichemo)
		elseif (ichemo == GROWTH_RATE) then
			val = cm3_pL*cp%dVdt
		elseif (ichemo == CELL_VOLUME) then
			val = Vcell_pL
		elseif (ichemo == O2_BY_VOL) then
			val = cp%Cin(OXYGEN)*Vcell_pL
		elseif (ichemo == CYCLE_PHASE) then
			val = cp%phase
		endif
		if (scale_drug_by_vol .and. (ichemo >= DRUG_A .and. ichemo <= DRUG_B + 2)) then
			val = val*Vcell_pL
		endif
		valmax(ict+1,ivar) = max(valmax(ict+1,ivar),val)	! cell type 1 or 2
		valmax(1,ivar) = max(valmax(1,ivar),val)			! both
		if (val <= 1.0e-8) then
			val_log = -8
		else
			val_log = log10(val)
		endif
		valmin_log(ict+1,ivar) = min(valmin_log(ict+1,ivar),val_log)	! cell type 1 or 2
		valmin_log(1,ivar) = min(valmin_log(1,ivar),val_log)			! both
		valmax_log(ict+1,ivar) = max(valmax_log(ict+1,ivar),val_log)	! cell type 1 or 2
		valmax_log(1,ivar) = max(valmax_log(1,ivar),val_log)			! both
	enddo
	n(ict+1) = n(ict+1) + 1
	n(1) = n(1) + 1
enddo
do ivar = 1,nvars
	ichemo = var_index(ivar)
	if (ichemo == CELL_VOLUME) then
		valmin(:,ivar) = Vcell_pL*0.75
		valmin_log(:,ivar) = log10(Vcell_pL*0.75)
	else
		valmin(:,ivar) = 0
	endif
enddo

dv = (valmax - valmin)/nhisto
!write(nflog,*) 'dv'
!write(nflog,'(e12.3)') dv
dv_log = (valmax_log - valmin_log)/nhisto
!write(nflog,*) 'dv_log'
!write(nflog,'(e12.3)') dv_log
phase_cnt = 0
do kcell = 1,nlist
	cp => cell_list(kcell)
	if (cp%state == DEAD) cycle
	ict = cp%celltype
	Vcell_pL = cm3_pL*cp%V
	do ivar = 1,nvars
		ichemo = var_index(ivar)
		if (ichemo == 0) then
			val = cp%CFSE
		elseif (ichemo <= MAX_CHEMO) then
			val = cp%Cin(ichemo)
		elseif (ichemo == GROWTH_RATE) then
			val = cm3_pL*cp%dVdt
		elseif (ichemo == CELL_VOLUME) then
			val = Vcell_pL
		elseif (ichemo == O2_BY_VOL) then
			val = cp%Cin(OXYGEN)*Vcell_pL
		elseif (ichemo == CYCLE_PHASE) then
			val = cp%phase
			phase_cnt(1,cp%phase) = phase_cnt(1,cp%phase) + 1
			phase_cnt(ict+1,cp%phase) = phase_cnt(ict+1,cp%phase) + 1
		endif
		if (scale_drug_by_vol .and. (ichemo >= DRUG_A .and. ichemo <= DRUG_B + 2)) then
			val = val*Vcell_pL
		endif
		k = (val-valmin(1,ivar))/dv(1,ivar) + 1
		k = min(k,nhisto)
		k = max(k,1)
		cnt(1,ivar,k) = cnt(1,ivar,k) + 1
		k = (val-valmin(ict+1,ivar))/dv(ict+1,ivar) + 1
		k = min(k,nhisto)
		k = max(k,1)
		cnt(ict+1,ivar,k) = cnt(ict+1,ivar,k) + 1
		if (val <= 1.0e-8) then
			val_log = -8
		else
			val_log = log10(val)
		endif
		k = (val_log-valmin_log(1,ivar))/dv_log(1,ivar) + 1
		k = min(k,nhisto)
		k = max(k,1)
		cnt_log(1,ivar,k) = cnt_log(1,ivar,k) + 1
		k = (val_log-valmin_log(ict+1,ivar))/dv_log(ict+1,ivar) + 1
		k = min(k,nhisto)
		k = max(k,1)
		cnt_log(ict+1,ivar,k) = cnt_log(ict+1,ivar,k) + 1
	enddo
enddo

do i = 1,3
	if (n(i) == 0) then
		vmin((i-1)*nvars+1:i*nvars) = 0
		vmax((i-1)*nvars+1:i*nvars) = 0
		histo_data((i-1)*nvars*nhisto+1:i*nhisto*nvars) = 0
		vmin_log((i-1)*nvars+1:i*nvars) = 0
		vmax_log((i-1)*nvars+1:i*nvars) = 0
		histo_data_log((i-1)*nvars*nhisto+1:i*nhisto*nvars) = 0
	else
		do ivar = 1,nvars
			vmin((i-1)*nvars+ivar) = valmin(i,ivar)
			vmax((i-1)*nvars+ivar) = valmax(i,ivar)
			do ih = 1,nhisto
				k = (i-1)*nvars*nhisto + (ivar-1)*nhisto + ih
				histo_data(k) = (100.*cnt(i,ivar,ih))/n(i)
			enddo
			vmin_log((i-1)*nvars+ivar) = valmin_log(i,ivar)
			vmax_log((i-1)*nvars+ivar) = valmax_log(i,ivar)
			do ih = 1,nhisto
				k = (i-1)*nvars*nhisto + (ivar-1)*nhisto + ih
				histo_data_log(k) = (100.*cnt_log(i,ivar,ih))/n(i)
			enddo
		enddo
	endif
	! Now overwrite CYCLE_PHASE data
	ivar = ivar_phase
	do ih = 1,nhisto
		k = (i-1)*nvars*nhisto + (ivar-1)*nhisto + ih
		if (ih <= 7) then
			histo_data(k) = (100.*phase_cnt(i,ih))/n(i)
		else
			histo_data(k) = 0
		endif
	enddo
enddo
deallocate(cnt)
deallocate(dv)
deallocate(valmin)
deallocate(valmax)
deallocate(cnt_log)
deallocate(dv_log)
deallocate(valmin_log)
deallocate(valmax_log)
!deallocate(histo_data_log)
end subroutine

!--------------------------------------------------------------------------------
! Returns the distribution of cell volume.
! nv is passed from the GUI
! Min divide volume = Vdivide0 - dVdivide
! therefore Vmin = (Vdivide0 - dVdivide)/2
! Max divide volume = Vmax = Vdivide0 + dVdivide
! dv = (Vmax - Vmin)/nv
! v0 = Vmin + dv/2
!--------------------------------------------------------------------------------
subroutine get_volprob(nv, v0, dv, prob) BIND(C)
!DEC$ ATTRIBUTES DLLEXPORT :: get_volprob
use, intrinsic :: iso_c_binding
integer(c_int) :: nv
real(c_double) :: v0, dv, prob(*)
integer :: n, kcell, k
real(REAL_KIND) :: v, Vmin, Vmax

!call logger('get_volprob')
Vmin = (Vdivide0 - dVdivide)/2
Vmax = Vdivide0 + dVdivide
dv = (Vmax - Vmin)/nv
v0 = Vmin + dv/2
prob(1:nv) = 0
n = 0
do kcell = 1,nlist
	if (cell_list(kcell)%state == DEAD) cycle
	v = cell_list(kcell)%V
	k = (v - Vmin)/dv + 1
	k = min(k,nv)
	k = max(k,1)
	prob(k) = prob(k) + 1
	n = n+1
enddo
prob(1:nv) = prob(1:nv)/n
end subroutine

!--------------------------------------------------------------------------------
! Returns the distribution of intracellular O2 level
! nv is passed from the GUI
!--------------------------------------------------------------------------------
subroutine get_oxyprob(nv, v0, dv, prob) BIND(C)
!DEC$ ATTRIBUTES DLLEXPORT :: get_oxyprob
use, intrinsic :: iso_c_binding
integer(c_int) :: nv
real(c_double) :: v0, dv, prob(*)
integer :: n, kcell, k
real(REAL_KIND) :: v, Vmin, Vmax, O2max

!call logger('get_oxyprob')
Vmin = 0
Vmax = chemo(OXYGEN)%bdry_conc
v0 = Vmin
dv = (Vmax - Vmin)/nv
prob(1:nv) = 0
n = 0
O2max = 0
do kcell = 1,nlist
	if (cell_list(kcell)%state == DEAD) cycle
	v = cell_list(kcell)%Cin(OXYGEN)
	O2max = max(O2max,v)
	k = (v - Vmin)/dv + 1
	k = min(k,nv)
	k = max(k,1)
	prob(k) = prob(k) + 1
	n = n+1
enddo
prob(1:nv) = prob(1:nv)/n
end subroutine

!--------------------------------------------------------------------------------
! Estimate average concentrations in the medium.
! In fact the average is over the coarse grid points that lie outside the blob.
!--------------------------------------------------------------------------------
subroutine getMediumConc(cmedium, cbdry)
real(REAL_KIND) :: cmedium(:), cbdry(:)
real(REAL_KIND) :: cntr(3), rng(3), radius, r2, d2
integer :: x, y, z
integer :: ixb, iyb, izb, nsum, ic, ichemo
logical :: bdry = .true.


cbdry(:) = chemo(:)%medium_Cbnd
call getBlobCentreRange(cntr,rng,radius)
r2 = radius**2
cmedium = 0
nsum = 0
do ixb = 1,NXB
	x = (ixb-1)*dxb
	do iyb = 1,NYB
		y = (iyb-1)*dxb
		do izb = 1,NZB
			z = (izb-1)*dxb
			d2 = (x - cntr(1))**2 + (y - cntr(2))**2 + (z - cntr(3))**2
			if (d2 > r2) then
				nsum = nsum + 1
				do ic = 1,nchemo
					ichemo = chemomap(ic)
					cmedium(ichemo) = cmedium(ichemo) + chemo(ichemo)%cave_b(ixb,iyb,izb)
				enddo
			endif
		enddo
	enddo
enddo
cmedium = cmedium/nsum
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
!subroutine getNecroticFraction(necrotic_fraction, totvol_cm3)
!real(REAL_KIND) :: necrotic_fraction, totvol_cm3
!real(REAL_KIND) :: dz, cellvol_cm3, dvol
!!real(REAL_KIND) :: necrotic_vol_cm3
!
!!call getNecroticVolume(necrotic_vol_cm3)
!!necrotic_fraction = necrotic_vol_cm3/vol_cm3
!dz = 2*Raverage
!cellvol_cm3 = Ncells/(1000.*Nmm3)	! Nmm3 = calibration target # of cells/mm3
!dvol = totvol_cm3-cellvol_cm3
!necrotic_fraction = dvol/totvol_cm3
!!write(*,'(a,i6,3e12.3,f6.3)') 'getNecroticFraction: ',Ncells,cellvol_cm3,totvol_cm3,dvol,necrotic_fraction
!if (necrotic_fraction < 0.005) necrotic_fraction = 0
!end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine getDiamVol(diam_cm,vol_cm3)
real(REAL_KIND) :: diam_cm, vol_cm3
real(REAL_KIND) :: area_cm2, radius, cntr(3), rng(3)
!call getBlobCentreRange(cntr,rng,radius)
!diam_cm = 2*radius
!vol_cm3 = (4*PI/3)*radius**3
call getVolume(5,vol_cm3, area_cm2)
diam_cm = 2*sqrt(area_cm2/PI)
end subroutine


!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine AverageCellMetabolism(r_G,r_P,r_A,r_I,P_util)
real(REAL_KIND) :: r_G,r_P,r_A,r_I,P_util
integer :: kcell, n, ityp
type(metabolism_type), pointer :: mp

r_G = 0
r_P = 0
r_A = 0
r_I = 0
P_util = 0
n = 0
do kcell = 1,nlist
	if (cell_list(kcell)%state == DEAD) cycle
	n = n+1
	mp => cell_list(kcell)%metab
	r_G = r_G + mp%G_rate/r_G_norm
	r_P = r_P + mp%P_rate/r_P_norm
	r_A = r_A + mp%A_rate/r_A_norm
	r_I = r_I + mp%I_rate/r_I_norm
	if (mp%G_rate > 0) then
		P_util = P_util + mp%P_rate/(2*(1-mp%f_G)*mp%G_rate)
	endif
enddo
r_G = r_G/n
r_P = r_P/n
r_A = r_A/n
r_I = r_I/n
P_util = p_util/n
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
function rint(i) result(r)
integer :: i
real(REAL_KIND) :: r
r = i
end function



!--------------------------------------------------------------------------------
! Compute total uptake rate for a constituent
!--------------------------------------------------------------------------------
subroutine sum_dMdt(ichemo)
integer :: ichemo
integer :: kcell, Nc
real(REAL_KIND) :: asum
type(cell_type), pointer :: cp

Nc = 0
asum = 0
do kcell = 1,nlist
	cp => cell_list(kcell)
	if (cp%state == DEAD) cycle
	Nc = Nc + 1
	asum = asum + cp%dMdt(ichemo)
!	if (ichemo==GLUCOSE .and. cp%dMdt(ichemo) < 1.0e-15) then
!		write(*,*) 'sum_dMdt: glucose small dMdt: ',kcell,cp%dMdt(ichemo)
!		stop
!	endif
enddo
total_dMdt = total_dMdt + asum
!write(*,'(a,2i6,2e12.3)') 'sum_dMdt: ',ichemo,Nc,asum,total_dMdt*3600
end subroutine

!!--------------------------------------------------------------------------------
!!--------------------------------------------------------------------------------
!subroutine getNviable(Nviable, Nlive)
!integer :: Nviable(:), Nlive(:)
!integer :: kcell, ityp, idrug
!logical :: tag
!type(cell_type), pointer :: cp
!
!Nviable = 0
!Nlive = 0
!do kcell = 1,nlist
!	cp => cell_list(kcell)
!	if (cp%state == DEAD) cycle
!    ityp = cp%celltype
!    Nlive(ityp) = Nlive(ityp) + 1
!!	if (cp%anoxia_tag .or. cp%aglucosia_tag .or. cp%radiation_tag .or. cp%state == DYING) cycle
!	if (cp%radiation_tag .or. cp%state == DYING) cycle
!    tag = .false.
!    do idrug = 1,ndrugs_used
!		if (cp%drug_tag(idrug)) tag = .true.
!	enddo
!	if (tag) cycle
!	Nviable(ityp) = Nviable(ityp) + 1
!enddo	
!end subroutine

!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------
subroutine getHypoxicCount(nhypoxic)
integer :: nhypoxic(3)
integer :: kcell, i

nhypoxic = 0
do kcell = 1,nlist
	if (cell_list(kcell)%state == DEAD) cycle
	do i = 1,3
		if (cell_list(kcell)%Cin(OXYGEN) < O2cutoff(i)) then
			nhypoxic(i) = nhypoxic(i) + 1
!			write(*,'(a,i6,2e12.3)') 'hypoxic: ',kcell,cell_list(kcell)%Cin(OXYGEN),O2cutoff(i)
		endif
	enddo
enddo
end subroutine

!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------
subroutine getClonoHypoxicCount(nclonohypoxic)
integer :: nclonohypoxic(3)
integer :: kcell, i, idrug
logical :: tagged

nclonohypoxic = 0
do kcell = 1,nlist
	if (cell_list(kcell)%state == DEAD) cycle
	if (cell_list(kcell)%state == DYING) cycle
!	if (cell_list(kcell)%anoxia_tag) cycle
!	if (cell_list(kcell)%aglucosia_tag) cycle
	if (cell_list(kcell)%radiation_tag) cycle
	tagged = .false.
	do idrug = 1,MAX_DRUGTYPES
	    if (cell_list(kcell)%drug_tag(idrug)) tagged = .true.
	enddo
	if (tagged) cycle
	do i = 1,3
		if (cell_list(kcell)%Cin(OXYGEN) < O2cutoff(i)) then
			nclonohypoxic(i) = nclonohypoxic(i) + 1
		endif
	enddo
enddo
end subroutine


!--------------------------------------------------------------------------------
! Returns all the extracellular concentrations along a line through the blob centre.
! Together with CFSE, growth rate (dVdt), cell volume,...
! Store the constituent profiles one after the other.
! The challenge is to find the appropriate projection/interpolation from the cells
! to the equi-spaced points on the line.
! Method:
! The (estimated) blob centre and range provides a line through the centre
! parallel to the X axis between the extremes of the blob in this direction.
! If the centre is c(3) and the range is rng(3), then the line is:
! x1 <= x <= x2  where x1 = c(1) - rng(1)/2, x2 = c(1) + rng(1)/2
! y = c(2)
! z = c(3)
! First the variable values need to be determined at all grid pts that enclose this line.
! In the case of extracellular concentrations, the grid pt values are used.
! For the non-concentration variables, i.e. those connected with cells, the
! grid pt values must be estimated by averaging over nearby cells.
! The relevant grid cells are those with "lower-left" corner pts given by:
! ixcnr = x1/dx + 1 ... x2/dx + 1
! iycnr = c(2)/dx + 1
! izcnr = c(3)/dx + 1
!--------------------------------------------------------------------------------
subroutine get_concdata1(nvars, ns, dxc, ex_conc) BIND(C)
!DEC$ ATTRIBUTES DLLEXPORT :: get_concdata
use, intrinsic :: iso_c_binding
integer(c_int) :: nvars, ns
real(c_double) :: dxc, ex_conc(0:*)
real(REAL_KIND) ::  dx, x1, x2, y0, z0, x, y, z
real(REAL_KIND) :: cntr(3), rng(3), radius
real(REAL_KIND) :: alfax, alfay, alfaz
integer :: nxgpts, ix1, ix2, iy0, iz0, kx, ky, kz
integer :: ix, iy, iz, k, ks, ichemo, kcell, offset
integer :: gcell(200)
integer, allocatable ::  ngc(:)
real(REAL_KIND), allocatable :: ctemp(:,:,:,:)
type(cell_type), pointer :: cp

!call logger('get_concdata')
dx = DELTA_X

call getBlobCentreRange(cntr,rng,radius)

x1 = cntr(1) - rng(1)/2
x2 = cntr(1) + rng(1)/2
y0 = cntr(2)
z0 = cntr(3)
ix1 = x1/dx + 1
ix2 = x2/dx + 1
iy0 = y0/dx + 1
iz0 = z0/dx + 1
nxgpts = ix2 - ix1 + 2
nvars = 1 + MAX_CHEMO + N_EXTRA
allocate(ngc(NX))
allocate(ctemp(nxgpts,2,2,0:nvars-1))

ichemo = OXYGEN
write(nflog,'(a,40f8.5)') 'EC O2: ',(Caverage(ix,iy0,iz0,ichemo),ix=ix1,ix2+1)

! First calculate averages at the enclosing grid pts
ctemp = 0
ngc = 0
do kx = 1,nxgpts
	do ky = 1,2
		do kz = 1,2
			ix = ix1 + kx - 1
			iy = iy0 + ky - 1
			iz = iz0 + kz - 1
			! Find the cells associated with this gridpt to compute cell-related variables
			call get_gridptcells(ix,iy,iz,ngc(ix),gcell)
			if (ngc(ix) > 0) then
				do k = 1,ngc(ix)
					kcell = gcell(k)
					cp => cell_list(kcell)
					ctemp(kx,ky,kz,0) = ctemp(kx,ky,kz,0) + cp%CFSE
					ctemp(kx,ky,kz,GROWTH_RATE) = ctemp(kx,ky,kz,GROWTH_RATE) + cp%dVdt
					ctemp(kx,ky,kz,CELL_VOLUME) = ctemp(kx,ky,kz,CELL_VOLUME) + cp%V
					ctemp(kx,ky,kz,O2_BY_VOL) = ctemp(kx,ky,kz,O2_BY_VOL) + cp%Cin(OXYGEN)*cp%V
				enddo
				ctemp(kx,ky,kz,0:nvars-1) = ctemp(kx,ky,kz,0:nvars-1)/ngc(ix)
			endif
			! For constituents simply use the concentration fields
			do ichemo = 1,MAX_CHEMO
				ctemp(kx,ky,kz,ichemo) = Caverage(ix,iy,iz,ichemo)
			enddo
		enddo
	enddo
enddo

ns = 20
dxc = (x2-x1)/(ns-1)

! Now need to interpolate values at ns points along the line from (x1,y0,z0) to (x2,y0,z0)
! alfax, alfay, alfaz each lies in (0,1)
y = y0
z = z0
iy = y/dx + 1
iz = z/dx + 1
alfay = (y - (iy-1)*dx)/dx
alfaz = (z - (iz-1)*dx)/dx
do ks = 1,ns
	x = x1 + (ks-1)*dxc
	ix = x/dx + 1
	alfax = (x - (ix-1)*dx)/dx
	kx = ix - ix1 + 1
!	if (kx >= nxgpts) then
!		write(*,*) 'bad kx: ',kx,nxgpts
!		write(*,'(a,i4,4f8.4,2i4)') 'ks,x,x1,x2,dxc,ix,kx: ',ks,x,x1,x2,dxc,ix,kx
!		stop
!	endif
	do ichemo = 0, nvars-1
		offset = ichemo*ns
		k = offset - 1 + ks
		ex_conc(k) = (1-alfax)*(1-alfay)*(1-alfaz)*ctemp(kx,1,1,ichemo) + &
					 (1-alfax)*alfay*(1-alfaz)*ctemp(kx,2,1,ichemo) + &
					 (1-alfax)*(1-alfay)*alfaz*ctemp(kx,1,2,ichemo) + &
					 (1-alfax)*alfay*alfaz*ctemp(kx,2,2,ichemo) + &
					 alfax*(1-alfay)*(1-alfaz)*ctemp(kx+1,1,1,ichemo) + &
					 alfax*alfay*(1-alfaz)*ctemp(kx+1,2,1,ichemo) + &
					 alfax*(1-alfay)*alfaz*ctemp(kx+1,1,2,ichemo) + &
					 alfax*alfay*alfaz*ctemp(kx+1,2,2,ichemo)
!		if (ichemo == OXYGEN) then
!		    write(nflog,'(a,5i6,f8.3)') 'get_concdata: ',ichemo,ix,ks,ns,k,ex_conc(k)
!		endif
!			if (ex_conc(k) > 0.18) then
!				write(*,'(2i4,2f8.4)') kp,k,x,ex_conc(k)
!				 write(*,'(3i4,4f8.4)') kx,1,1,(1-alfax),(1-alfay),(1-alfaz),ctemp(kx,1,1,ichemo)
!				 write(*,'(3i4,4f8.4)') kx,2,1,(1-alfax),alfay,(1-alfaz),ctemp(kx,2,1,ichemo)
!				 write(*,'(3i4,4f8.4)') kx,1,2,(1-alfax),(1-alfay),alfaz,ctemp(kx,1,2,ichemo)
!				 write(*,'(3i4,4f8.4)') kx,2,2,(1-alfax),alfay,alfaz,ctemp(kx,2,2,ichemo)
!				 write(*,'(3i4,4f8.4)') kx+1,1,1,alfax,(1-alfay),(1-alfaz),ctemp(kx+1,1,1,ichemo)
!				 write(*,'(3i4,4f8.4)') kx+1,2,1,alfax,alfay,(1-alfaz),ctemp(kx+1,2,1,ichemo)
!				 write(*,'(3i4,4f8.4)') kx+1,1,2,alfax,(1-alfay),alfaz,ctemp(kx+1,1,2,ichemo)
!				 write(*,'(3i4,4f8.4)') kx+1,2,2,alfax,alfay,alfaz,ctemp(kx+1,2,2,ichemo)
!				stop
!			endif
!		endif
	enddo
enddo
!write(nflog,'(10f8.3)') ex_conc(0:nvars*ns-1)
deallocate(ctemp)
deallocate(ngc)
end subroutine

!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------
subroutine get_concdata(axis, nvars, ns, dxc, ex_conc) BIND(C)
!DEC$ ATTRIBUTES DLLEXPORT :: get_concdata
use, intrinsic :: iso_c_binding
integer(c_int) :: axis, nvars, ns
real(c_double) :: dxc, ex_conc(0:*)
real(REAL_KIND) ::  dx, x1, x2, y1, y2, z1, z2, x0, y0, z0, x, y, z
real(REAL_KIND) :: cntr(3), rng(3), radius
real(REAL_KIND) :: alfax, alfay, alfaz
integer :: nxgpts, ix1, ix2, iy1, iy2, iz1, iz2, ix0, iy0, iz0, kx, ky, kz
integer :: ix, iy, iz, k, ks, ichemo, kcell, offset
integer :: gcell(200)
integer, allocatable ::  ngc(:)
real(REAL_KIND), allocatable :: ctemp(:,:,:,:)
type(cell_type), pointer :: cp

!call logger('get_concdata')
dx = DELTA_X
ns = 20

call getBlobCentreRange(cntr,rng,radius)

if (axis == 1) then
	x1 = cntr(1) - rng(1)/2
	x2 = cntr(1) + rng(1)/2
	y0 = cntr(2)
	z0 = cntr(3)
	ix1 = x1/dx + 1
	ix2 = x2/dx + 1
	iy0 = y0/dx + 1
	iz0 = z0/dx + 1
	nxgpts = ix2 - ix1 + 2
	nvars = 1 + MAX_CHEMO + N_EXTRA
	allocate(ngc(NX))
	allocate(ctemp(nxgpts,2,2,0:nvars-1))

	ichemo = OXYGEN

	! First calculate averages at the enclosing grid pts
	ctemp = 0
	ngc = 0
	do kx = 1,nxgpts
		do ky = 1,2
			do kz = 1,2
				ix = ix1 + kx - 1
				iy = iy0 + ky - 1
				iz = iz0 + kz - 1
				! Find the cells associated with this gridpt to compute cell-related variables
				call get_gridptcells(ix,iy,iz,ngc(ix),gcell)
				if (ngc(ix) > 0) then
					do k = 1,ngc(ix)
						kcell = gcell(k)
						cp => cell_list(kcell)
						ctemp(kx,ky,kz,0) = ctemp(kx,ky,kz,0) + cp%CFSE
						ctemp(kx,ky,kz,GROWTH_RATE) = ctemp(kx,ky,kz,GROWTH_RATE) + cp%dVdt
						ctemp(kx,ky,kz,CELL_VOLUME) = ctemp(kx,ky,kz,CELL_VOLUME) + cp%V
						ctemp(kx,ky,kz,O2_BY_VOL) = ctemp(kx,ky,kz,O2_BY_VOL) + cp%Cin(OXYGEN)*cp%V
					enddo
					ctemp(kx,ky,kz,0:nvars-1) = ctemp(kx,ky,kz,0:nvars-1)/ngc(ix)
				endif
				! For constituents simply use the concentration fields
				do ichemo = 1,MAX_CHEMO
					ctemp(kx,ky,kz,ichemo) = Caverage(ix,iy,iz,ichemo)
				enddo
			enddo
		enddo
	enddo

	dxc = (x2-x1)/(ns-1)

	! Now need to interpolate values at ns points along the line from (x1,y0,z0) to (x2,y0,z0)
	! alfax, alfay, alfaz each lies in (0,1)
	y = y0
	z = z0
	iy = y/dx + 1
	iz = z/dx + 1
	alfay = (y - (iy-1)*dx)/dx
	alfaz = (z - (iz-1)*dx)/dx
	do ks = 1,ns
		x = x1 + (ks-1)*dxc
		ix = x/dx + 1
		alfax = (x - (ix-1)*dx)/dx
		kx = ix - ix1 + 1
	!	if (kx >= nxgpts) then
	!		write(*,*) 'bad kx: ',kx,nxgpts
	!		write(*,'(a,i4,4f8.4,2i4)') 'ks,x,x1,x2,dxc,ix,kx: ',ks,x,x1,x2,dxc,ix,kx
	!		stop
	!	endif
		do ichemo = 0, nvars-1
			offset = ichemo*ns
			k = offset - 1 + ks
			ex_conc(k) = (1-alfax)*(1-alfay)*(1-alfaz)*ctemp(kx,1,1,ichemo) + &
						 (1-alfax)*alfay*(1-alfaz)*ctemp(kx,2,1,ichemo) + &
						 (1-alfax)*(1-alfay)*alfaz*ctemp(kx,1,2,ichemo) + &
						 (1-alfax)*alfay*alfaz*ctemp(kx,2,2,ichemo) + &
						 alfax*(1-alfay)*(1-alfaz)*ctemp(kx+1,1,1,ichemo) + &
						 alfax*alfay*(1-alfaz)*ctemp(kx+1,2,1,ichemo) + &
						 alfax*(1-alfay)*alfaz*ctemp(kx+1,1,2,ichemo) + &
						 alfax*alfay*alfaz*ctemp(kx+1,2,2,ichemo)
		enddo
	enddo
	!write(nflog,'(10f8.3)') ex_conc(0:nvars*ns-1)
	deallocate(ctemp)
	deallocate(ngc)
elseif (axis == 2) then
	y1 = cntr(2) - rng(2)/2
	y2 = cntr(2) + rng(2)/2
	x0 = cntr(1)
	z0 = cntr(3)
	iy1 = y1/dx + 1
	iy2 = y2/dx + 1
	ix0 = x0/dx + 1
	iz0 = z0/dx + 1
	nxgpts = iy2 - iy1 + 2
	nvars = 1 + MAX_CHEMO + N_EXTRA
	allocate(ngc(NX))
	allocate(ctemp(2,nxgpts,2,0:nvars-1))

	ichemo = OXYGEN
	write(nflog,'(a,40f8.5)') 'EC y O2: ',(Caverage(ix0,iy,iz0,ichemo),iy=iy1,iy2+1)

	! First calculate averages at the enclosing grid pts
	ctemp = 0
	ngc = 0
	do ky = 1,nxgpts
		do kx = 1,2
			do kz = 1,2
				iy = iy1 + ky - 1
				ix = ix0 + kx - 1
				iz = iz0 + kz - 1
				! Find the cells associated with this gridpt to compute cell-related variables
				call get_gridptcells(ix,iy,iz,ngc(iy),gcell)
				if (ngc(ix) > 0) then
					do k = 1,ngc(iy)
						kcell = gcell(k)
						cp => cell_list(kcell)
						ctemp(kx,ky,kz,0) = ctemp(kx,ky,kz,0) + cp%CFSE
						ctemp(kx,ky,kz,GROWTH_RATE) = ctemp(kx,ky,kz,GROWTH_RATE) + cp%dVdt
						ctemp(kx,ky,kz,CELL_VOLUME) = ctemp(kx,ky,kz,CELL_VOLUME) + cp%V
						ctemp(kx,ky,kz,O2_BY_VOL) = ctemp(kx,ky,kz,O2_BY_VOL) + cp%Cin(OXYGEN)*cp%V
					enddo
					ctemp(kx,ky,kz,0:nvars-1) = ctemp(kx,ky,kz,0:nvars-1)/ngc(ix)
				endif
				! For constituents simply use the concentration fields
				do ichemo = 1,MAX_CHEMO
					ctemp(kx,ky,kz,ichemo) = Caverage(ix,iy,iz,ichemo)
				enddo
			enddo
		enddo
	enddo

	dxc = (y2-y1)/(ns-1)

	! Now need to interpolate values at ns points along the line from (x0,y1,z0) to (x0,y2,z0)
	! alfax, alfay, alfaz each lies in (0,1)
	x = x0
	z = z0
	ix = x/dx + 1
	iz = z/dx + 1
	alfax = (x - (ix-1)*dx)/dx
	alfaz = (z - (iz-1)*dx)/dx
	do ks = 1,ns
		y = y1 + (ks-1)*dxc
		iy = y/dx + 1
		alfay = (y - (iy-1)*dx)/dx
		ky = iy - iy1 + 1
		do ichemo = 0, nvars-1
			offset = ichemo*ns
			k = offset - 1 + ks
			ex_conc(k) = (1-alfay)*(1-alfax)*(1-alfaz)*ctemp(1,ky,1,ichemo) + &
						 (1-alfay)*alfax*(1-alfaz)*ctemp(2,ky,1,ichemo) + &
						 (1-alfay)*(1-alfax)*alfaz*ctemp(1,ky,2,ichemo) + &
						 (1-alfay)*alfax*alfaz*ctemp(2,ky,2,ichemo) + &
						 alfay*(1-alfax)*(1-alfaz)*ctemp(1,ky+1,1,ichemo) + &
						 alfay*alfax*(1-alfaz)*ctemp(2,ky+1,1,ichemo) + &
						 alfay*(1-alfax)*alfaz*ctemp(1,ky+1,2,ichemo) + &
						 alfay*alfax*alfaz*ctemp(2,ky+1,2,ichemo)
		enddo
	enddo
	deallocate(ctemp)
	deallocate(ngc)
elseif (axis == 3) then
	z1 = cntr(3) - rng(3)/2
	z2 = cntr(3) + rng(3)/2
	y0 = cntr(2)
	x0 = cntr(1)
	iz1 = z1/dx + 1
	iz2 = z2/dx + 1
	iy0 = y0/dx + 1
	ix0 = x0/dx + 1
	nxgpts = iz2 - iz1 + 2
	nvars = 1 + MAX_CHEMO + N_EXTRA
	allocate(ngc(NX))
	allocate(ctemp(2,2,nxgpts,0:nvars-1))

	ichemo = OXYGEN
	write(nflog,'(a,40f8.5)') 'EC z O2: ',(Caverage(ix0,iy0,iz,ichemo),iz=iz1,iz2+1)

	! First calculate averages at the enclosing grid pts
	ctemp = 0
	ngc = 0
	do kz = 1,nxgpts
		do ky = 1,2
			do kx = 1,2
				iz = iz1 + kz - 1
				iy = iy0 + ky - 1
				ix = ix0 + kx - 1
				! Find the cells associated with this gridpt to compute cell-related variables
				call get_gridptcells(ix,iy,iz,ngc(iz),gcell)
				if (ngc(iz) > 0) then
					do k = 1,ngc(iz)
						kcell = gcell(k)
						cp => cell_list(kcell)
						ctemp(kx,ky,kz,0) = ctemp(kx,ky,kz,0) + cp%CFSE
						ctemp(kx,ky,kz,GROWTH_RATE) = ctemp(kx,ky,kz,GROWTH_RATE) + cp%dVdt
						ctemp(kx,ky,kz,CELL_VOLUME) = ctemp(kx,ky,kz,CELL_VOLUME) + cp%V
						ctemp(kx,ky,kz,O2_BY_VOL) = ctemp(kx,ky,kz,O2_BY_VOL) + cp%Cin(OXYGEN)*cp%V
					enddo
					ctemp(kx,ky,kz,0:nvars-1) = ctemp(kx,ky,kz,0:nvars-1)/ngc(ix)
				endif
				! For constituents simply use the concentration fields
				do ichemo = 1,MAX_CHEMO
					ctemp(kx,ky,kz,ichemo) = Caverage(ix,iy,iz,ichemo)
				enddo
			enddo
		enddo
	enddo

	dxc = (z2-z1)/(ns-1)

	! Now need to interpolate values at ns points along the line from (x0,y0,z1) to (x0,y0,z2)
	! alfax, alfay, alfaz each lies in (0,1)
	y = y0
	x = x0
	iy = y/dx + 1
	ix = x/dx + 1
	alfay = (y - (iy-1)*dx)/dx
	alfax = (x - (ix-1)*dx)/dx
	do ks = 1,ns
		z = z1 + (ks-1)*dxc
		iz = z/dx + 1
		alfaz = (z - (iz-1)*dx)/dx
		kz = iz - iz1 + 1
		do ichemo = 0, nvars-1
			offset = ichemo*ns
			k = offset - 1 + ks
			ex_conc(k) = (1-alfaz)*(1-alfay)*(1-alfax)*ctemp(1,1,kz,ichemo) + &
						 (1-alfaz)*alfay*(1-alfax)*ctemp(2,1,kz,ichemo) + &
						 (1-alfaz)*(1-alfay)*alfax*ctemp(1,2,kz,ichemo) + &
						 (1-alfaz)*alfay*alfax*ctemp(2,2,kz,ichemo) + &
						 alfaz*(1-alfay)*(1-alfax)*ctemp(1,1,kz+1,ichemo) + &
						 alfaz*alfay*(1-alfax)*ctemp(2,1,kz+1,ichemo) + &
						 alfaz*(1-alfay)*alfax*ctemp(1,2,kz+1,ichemo) + &
						 alfaz*alfay*alfax*ctemp(2,2,kz+1,ichemo)
		enddo
	enddo
	deallocate(ctemp)
	deallocate(ngc)
endif
end subroutine

!--------------------------------------------------------------------------------
! Returns all the intracellular concentrations along a line through the blob centre.
! Need to find all cells within a yz-tube about the centreline, then average over
! blocks of width dx.  A yz-tube has specified radius tube_radius, centred on centreline.
! Together with CFSE, growth rate (dVdt), cell volume,...
! Store the constituent profiles one after the other.
!--------------------------------------------------------------------------------
subroutine get_IC_concdata1(nvars, ns, dxc, ic_conc) BIND(C)
!DEC$ ATTRIBUTES DLLEXPORT :: get_ic_concdata
use, intrinsic :: iso_c_binding
integer(c_int) :: nvars, ns
real(c_double) :: dxc, ic_conc(0:*)
integer :: k, ks, ichemo, kcell, offset
real(REAL_KIND) :: cntr(3), rng(3), radius, tube_radius, tube_radius2
real(REAL_KIND) :: x1, x2, y0, z0, xdiff, xmin, xmax, r2
real(REAL_KIND), allocatable :: csum(:,:)
integer, allocatable :: ncsum(:)
type(cell_type), pointer :: cp

!call logger('get_IC_concdata')
nvars = 1 + MAX_CHEMO + N_EXTRA
tube_radius = DELTA_X
tube_radius2 = tube_radius*tube_radius

call getBlobCentreRange(cntr,rng,radius)

x1 = cntr(1) - rng(1)/2
x2 = cntr(1) + rng(1)/2
y0 = cntr(2)
z0 = cntr(3)

ns = 20
dxc = (x2-x1)/ns
!dxc = DELTA_X
!ns = (x2-x1)/dxc
xdiff = (x2-x1) - ns*dxc
xmin = x1 + xdiff/2
xmax = x2 - xdiff/2
! x range for block ks is xmin + ks.dx -> xmin + (ks+1)dx for ks: 0,ns-1
allocate(csum(ns,0:nvars-1))
allocate(ncsum(ns))
csum = 0
ncsum = 0
do kcell = 1,nlist
	if (cell_list(kcell)%state == DEAD) cycle
	cp => cell_list(kcell)
	if (cp%centre(1,1) < xmin .or. cp%centre(1,1) > xmax) cycle
    r2 = (cp%centre(2,1) - y0)**2 + (cp%centre(3,1) - z0)**2
    if (r2 < tube_radius2) then
        ks = (cp%centre(1,1) - xmin)/dxc + 1
        ncsum(ks) = ncsum(ks) + 1
        csum(ks,CFSE) = csum(ks,CFSE) + cp%CFSE
        csum(ks,1:MAX_CHEMO) = csum(ks,1:MAX_CHEMO) + cp%Cin
        csum(ks,GROWTH_RATE) = csum(ks,GROWTH_RATE) + cp%dVdt
        csum(ks,CELL_VOLUME) = csum(ks,CELL_VOLUME) + cm3_pL*cp%V
        csum(ks,O2_BY_VOL) = csum(ks,O2_BY_VOL) + cp%Cin(OXYGEN)*cm3_pL*cp%V
    endif
enddo
do ks = 1,ns
	do ichemo = 0, nvars-1
		offset = ichemo*ns
		k = offset - 1 + ks
		if (ncsum(ks) > 0) then
            ic_conc(k) = csum(ks,ichemo)/ncsum(ks)
        else
            ic_conc(k) = 0
        endif
    enddo
enddo

end subroutine

!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------
subroutine get_IC_concdata(axis, nvars, ns, dxc, ic_conc) BIND(C)
!DEC$ ATTRIBUTES DLLEXPORT :: get_ic_concdata
use, intrinsic :: iso_c_binding
integer(c_int) :: axis, nvars, ns
real(c_double) :: dxc, ic_conc(0:*)
integer :: k, ks, ichemo, kcell, offset
real(REAL_KIND) :: cntr(3), rng(3), radius, tube_radius, tube_radius2
real(REAL_KIND) :: x1, x2, y1, y2, z1, z2, x0, y0, z0, xdiff, ydiff, zdiff, xmin, xmax, ymin, ymax, zmin, zmax, r2
real(REAL_KIND), allocatable :: csum(:,:)
integer, allocatable :: ncsum(:)
type(cell_type), pointer :: cp

!call logger('get_IC_concdata')
nvars = 1 + MAX_CHEMO + N_EXTRA
tube_radius = DELTA_X
tube_radius2 = tube_radius*tube_radius
ns = 20

call getBlobCentreRange(cntr,rng,radius)

if (axis == 1) then
	x1 = cntr(1) - rng(1)/2
	x2 = cntr(1) + rng(1)/2
	y0 = cntr(2)
	z0 = cntr(3)

	dxc = (x2-x1)/ns
	!dxc = DELTA_X
	!ns = (x2-x1)/dxc
	xdiff = (x2-x1) - ns*dxc
	xmin = x1 + xdiff/2
	xmax = x2 - xdiff/2
	! x range for block ks is xmin + ks.dx -> xmin + (ks+1)dx for ks: 0,ns-1
	allocate(csum(ns,0:nvars-1))
	allocate(ncsum(ns))
	csum = 0
	ncsum = 0
	do kcell = 1,nlist
		if (cell_list(kcell)%state == DEAD) cycle
		cp => cell_list(kcell)
		if (cp%centre(1,1) < xmin .or. cp%centre(1,1) > xmax) cycle
		r2 = (cp%centre(2,1) - y0)**2 + (cp%centre(3,1) - z0)**2
		if (r2 < tube_radius2) then
			ks = (cp%centre(1,1) - xmin)/dxc + 1
			ncsum(ks) = ncsum(ks) + 1
			csum(ks,CFSE) = csum(ks,CFSE) + cp%CFSE
			csum(ks,1:MAX_CHEMO) = csum(ks,1:MAX_CHEMO) + cp%Cin
			csum(ks,GROWTH_RATE) = csum(ks,GROWTH_RATE) + cm3_pL*cp%dVdt
			csum(ks,CELL_VOLUME) = csum(ks,CELL_VOLUME) + cm3_pL*cp%V
			csum(ks,O2_BY_VOL) = csum(ks,O2_BY_VOL) + cp%Cin(OXYGEN)*cm3_pL*cp%V
		endif
	enddo
	do ks = 1,ns
		do ichemo = 0, nvars-1
			offset = ichemo*ns
			k = offset - 1 + ks
			if (ncsum(ks) > 0) then
				ic_conc(k) = csum(ks,ichemo)/ncsum(ks)
			else
				ic_conc(k) = 0
			endif
		enddo
	enddo
elseif (axis == 2) then
	y1 = cntr(2) - rng(2)/2
	y2 = cntr(2) + rng(2)/2
	x0 = cntr(1)
	z0 = cntr(3)

	dxc = (y2-y1)/ns
	ydiff = (y2-y1) - ns*dxc
	ymin = y1 + ydiff/2
	ymax = y2 - ydiff/2
	! y range for block ks is ymin + ks.dx -> ymin + (ks+1)dx for ks: 0,ns-1
	allocate(csum(ns,0:nvars-1))
	allocate(ncsum(ns))
	csum = 0
	ncsum = 0
	do kcell = 1,nlist
		if (cell_list(kcell)%state == DEAD) cycle
		cp => cell_list(kcell)
		if (cp%centre(2,1) < ymin .or. cp%centre(2,1) > ymax) cycle
		r2 = (cp%centre(1,1) - x0)**2 + (cp%centre(3,1) - z0)**2
		if (r2 < tube_radius2) then
			ks = (cp%centre(2,1) - ymin)/dxc + 1
			ncsum(ks) = ncsum(ks) + 1
			csum(ks,CFSE) = csum(ks,CFSE) + cp%CFSE
			csum(ks,1:MAX_CHEMO) = csum(ks,1:MAX_CHEMO) + cp%Cin
			csum(ks,GROWTH_RATE) = csum(ks,GROWTH_RATE) + cp%dVdt
			csum(ks,CELL_VOLUME) = csum(ks,CELL_VOLUME) + cm3_pL*cp%V
			csum(ks,O2_BY_VOL) = csum(ks,O2_BY_VOL) + cp%Cin(OXYGEN)*cm3_pL*cp%V
		endif
	enddo
	do ks = 1,ns
		do ichemo = 0, nvars-1
			offset = ichemo*ns
			k = offset - 1 + ks
			if (ncsum(ks) > 0) then
				ic_conc(k) = csum(ks,ichemo)/ncsum(ks)
			else
				ic_conc(k) = 0
			endif
		enddo
	enddo
elseif (axis == 3) then
	z1 = cntr(3) - rng(3)/2
	z2 = cntr(3) + rng(3)/2
	y0 = cntr(2)
	x0 = cntr(1)

	dxc = (z2-z1)/ns
	zdiff = (z2-z1) - ns*dxc
	zmin = z1 + zdiff/2
	zmax = z2 - zdiff/2
	! z range for block ks is zmin + ks.dx -> zmin + (ks+1)dx for ks: 0,ns-1
	allocate(csum(ns,0:nvars-1))
	allocate(ncsum(ns))
	csum = 0
	ncsum = 0
	do kcell = 1,nlist
		if (cell_list(kcell)%state == DEAD) cycle
		cp => cell_list(kcell)
		if (cp%centre(3,1) < zmin .or. cp%centre(3,1) > zmax) cycle
		r2 = (cp%centre(2,1) - y0)**2 + (cp%centre(1,1) - x0)**2
		if (r2 < tube_radius2) then
			ks = (cp%centre(3,1) - zmin)/dxc + 1
			ncsum(ks) = ncsum(ks) + 1
			csum(ks,CFSE) = csum(ks,CFSE) + cp%CFSE
			csum(ks,1:MAX_CHEMO) = csum(ks,1:MAX_CHEMO) + cp%Cin
			csum(ks,GROWTH_RATE) = csum(ks,GROWTH_RATE) + cp%dVdt
			csum(ks,CELL_VOLUME) = csum(ks,CELL_VOLUME) + cm3_pL*cp%V
			csum(ks,O2_BY_VOL) = csum(ks,O2_BY_VOL) + cp%Cin(OXYGEN)*cm3_pL*cp%V
		endif
	enddo
	do ks = 1,ns
		do ichemo = 0, nvars-1
			offset = ichemo*ns
			k = offset - 1 + ks
			if (ncsum(ks) > 0) then
				ic_conc(k) = csum(ks,ichemo)/ncsum(ks)
			else
				ic_conc(k) = 0
			endif
		enddo
	enddo
endif

end subroutine

!--------------------------------------------------------------------------------
! Returns all the extracellular concentrations along a line through the blob centre.
! Together with CFSE, growth rate (dVdt), cell volume,...
! Store the constituent profiles one after the other.
! The challenge is to find the appropriate projection/interpolation from the cells
! to the equi-spaced points on the line.
!--------------------------------------------------------------------------------
subroutine get_concdata2(nvars, ns, dx, ex_conc) BIND(C)
!DEC$ ATTRIBUTES DLLEXPORT :: get_concdata
use, intrinsic :: iso_c_binding
integer(c_int) :: nvars, ns
real(c_double) :: dx, ex_conc(0:*)
real(REAL_KIND) :: cbnd, cmin = 1.0e-6
integer :: ix, iy, iz, i, ic, k, ichemo, kcell, x, y, z, x1, x2, offset
integer :: gcell(100)
integer, allocatable ::  ngc(:)
real(REAL_KIND), allocatable :: ctemp(:,:)
type(cell_type), pointer :: cp

!call logger('get_concdata')
nvars = 1 + MAX_CHEMO + N_EXTRA
allocate(ngc(NX))
allocate(ctemp(NX,0:nvars-1))

dx = DELTA_X

iy = NY/2+1
iz = NZ/2+1
!write(*,*) 'get_concdata: iy,iz: ',iy,iz

! First calculate averages at the grid pts

x1 = 0
x2 = 0
ctemp = 0
ngc = 0
do ix = 1,NX
	call get_gridptcells(ix,iy,iz,ngc(ix),gcell)
	if (ngc(ix) == 0) cycle
	if (x1 == 0) x1 = ix
	do k = 1,ngc(ix)
		kcell = gcell(k)
		cp => cell_list(kcell)
		ctemp(ix,0) = ctemp(ix,0) + cp%CFSE
		ctemp(ix,1:MAX_CHEMO) = ctemp(ix,1:MAX_CHEMO) + cp%Cex(1:MAX_CHEMO)
		ctemp(ix,GROWTH_RATE) = ctemp(ix,GROWTH_RATE) + cp%dVdt
		ctemp(ix,CELL_VOLUME) = ctemp(ix,CELL_VOLUME) + cp%V
		ctemp(ix,O2_BY_VOL) = ctemp(ix,O2_BY_VOL) + cp%Cin(OXYGEN)*cp%V
	enddo
	ctemp(ix,:) = ctemp(ix,:)/ngc(ix)
	x2 = ix
enddo
do ix = x1,x2
	if (ngc(ix) == 0) then
		ctemp(ix,:) = ctemp(ix-1,:)
	endif
enddo
!write(*,'(10i6)') ngc(:)
!write(*,'(10f7.3)') ctemp(:,GLUCOSE)
ns = x2 - x1 + 1
!write(*,*) 'x1, x2, ns: ',x1,x2,ns

do ichemo = 0,nvars-1
	offset = ichemo*ns
	k = offset - 1
	do x = x1, x2
		k = k + 1
		ex_conc(k) = ctemp(x,ichemo)
!		if (ichemo == 1) write(*,*) x,k,ex_conc(k)
	enddo
enddo

deallocate(ctemp)
deallocate(ngc)

#if 0
rng(:,1) = Centre(:) - (adrop*Radius + 2)
rng(:,2) = Centre(:) + (adrop*Radius + 2)
!rng(axis,:) = Centre(axis) + fraction*Radius
y = Centre(2) + 0.5
z = Centre(3) + 0.5

! First need to establish the range of x that is inside the blob: (x1,x2)
x1 = 0
do x = rng(1,1),rng(1,2)
	kcell = occupancy(x,y,z)%indx(1)
	if (kcell <= OUTSIDE_TAG) then
		if (x1 == 0) then
			cycle
		else
			exit
		endif
	elseif (x1 == 0) then
		x1 = x
	endif
	x2 = x
enddo

ns = x2 - x1 + 1 
do ichemo = 0,nvars-1
	offset = ichemo*ns
	k = offset - 1
	do x = x1, x2
		k = k + 1
		kcell = occupancy(x,y,z)%indx(1)
		if (kcell <= OUTSIDE_TAG) then
			ex_conc(k) = 0
			cycle
		endif
		i = ODEdiff%ivar(x,y,z)
        if (ichemo == 0) then	! CFSE
			if (kcell > 0) then
				ex_conc(k) = cell_list(kcell)%CFSE
			else
				ex_conc(k) = 0
			endif
       elseif (ichemo <= MAX_CHEMO) then
			if (chemo(ichemo)%used) then
				if (i > 0) then
					ex_conc(k) = allstate(i,ichemo)	
				else
					ex_conc(k) = 0
				endif
			else
				ex_conc(k) = 0
			endif
        elseif (ichemo == GROWTH_RATE) then	
			if (kcell > 0) then
				ex_conc(k) = cell_list(kcell)%dVdt
			else
				ex_conc(k) = 0
			endif
        elseif (ichemo == CELL_VOLUME) then	
			if (kcell > 0) then
				ex_conc(k) = Vcell_pL*cell_list(kcell)%V
			else
				ex_conc(k) = 0
			endif
        elseif (ichemo == O2_BY_VOL) then	
			if (kcell > 0) then
				ex_conc(k) = allstate(i,OXYGEN)*Vcell_pL*cell_list(kcell)%V
			else
				ex_conc(k) = 0
			endif
		endif
    enddo
enddo
#endif
! Add concentrations at the two boundaries 
! At ns=1, at at ns=ns+1
!ns = ns+1
!do ic = 1,nvars
!	ichemo = ic - 1
!	if (ichemo == 0) then	! CFSE
!		k = ic	
!		ex_conc(k) = 0
!        k = (ns-1)*nvars + ic
!        ex_conc(k) = 0	
!	elseif (ichemo <= MAX_CHEMO) then
!		k = ic
!		if (chemo(ichemo)%used) then
!			ex_conc(k) = BdryConc(ichemo,t_simulation)
!		else
!			ex_conc(k) = 0
!		endif      
!		k = (ns-1)*nvars + ic
!		if (chemo(ichemo)%used) then
!			ex_conc(k) = BdryConc(ichemo,t_simulation)
!		else
!			ex_conc(k) = 0
!		endif      
!    elseif (ichemo == MAX_CHEMO+1) then	! growth rate
!		k = ic
!		ex_conc(k) = 0
!! 			write(nflog,'(a,4i6,f8.3)') 'Growth rate: ',x,1,i,k,ex_conc(k)
!        k = (ns-1)*nvars + ic
!        ex_conc(k) = 0
!! 			write(nflog,'(a,4i6,f8.3)') 'Growth rate: ',x,ns,i,k,ex_conc(k)
!    elseif (ichemo == MAX_CHEMO+2) then	! cell volume
!		k = ic
!		ex_conc(k) = 0
!! 			write(nflog,'(a,4i6,f8.3)') 'Growth rate: ',x,1,i,k,ex_conc(k)
!        k = (ns-1)*nvars + ic
!        ex_conc(k) = 0
!! 			write(nflog,'(a,4i6,f8.3)') 'Growth rate: ',x,ns,i,k,ex_conc(k)
!    endif
!enddo
end subroutine

!--------------------------------------------------------------------------------
! Generates a list of cells with centres within the "cube of influence" of the
! grid pt (ix,iy,iz).
! The boundaries of the grid pt (ix,iy,iz) "cube of influence" (dotted lines) are:
!	x in (ix-1)*dxf - dxf/2, (ix-1)*dxf + dxf/2
!	y in (iy-1)*dxf - dyf/2, (iy-1)*dyf + dyf/2
!	z in (iz-1)*dxf - dzf/2, (iz-1)*dzf + dzf/2
! (respecting grid boundaries)
!--------------------------------------------------------------------------------
subroutine get_gridptcells(ix,iy,iz,ngc,gcell)
integer :: ix,iy,iz,ngc,gcell(:)
integer :: idx, idy, idz, ixx, iyy, izz, nc, k, kcell
real(REAL_KIND) :: xmax, ymax, zmax, x1, x2, y1, y2, z1, z2, c(3)
type(cell_type), pointer :: cp

xmax = (NX-1)*dxf
ymax = (NY-1)*dxf
zmax = (NZ-1)*dxf
x1 = max(0.0,(ix-1.5)*dxf)
x2 = min(xmax,x1 + dxf)
y1 = max(0.0,(iy-1.5)*dxf)
y2 = min(ymax,y1 + dxf)
z1 = max(0.0,(iz-1.5)*dxf)
z2 = min(zmax,z1 + dxf)
ngc = 0
do idx = -1,0
	ixx = ix + idx
	if (ixx < 1 .or. ixx == NX) cycle
	do idy = -1,0
		iyy = iy + idy
		if (iyy < 1 .or. iyy == NY) cycle
		do idz = -1,0
			izz = iz + idz
			if (izz < 1 .or. izz == NZ) cycle
			nc = grid(ixx,iyy,izz)%nc
			if (nc == 0) cycle
			do k = 1,nc
				kcell = grid(ixx,iyy,izz)%cell(k)
				cp => cell_list(kcell)
				if (cp%state == DEAD) cycle
				if (cp%nspheres == 1) then
					c = cp%centre(:,1)
				else
					c = 0.5*(cp%centre(:,1) + cp%centre(:,2))
				endif
				if ((c(1) < x1 .or. c(1) > x2) .or. (c(2) < y1 .or. c(2) > y2) .or. (c(3) < z1 .or. c(3) > z2)) cycle
				ngc = ngc + 1
				gcell(ngc) = kcell
			enddo
		enddo
	enddo
enddo
end subroutine

!-----------------------------------------------------------------------------------------
! Estimate necrotic volume by determining the range in the three axis directions
! Try to improve the method so as to handle "no-hit" cases.
! The idea is to use a sheaf of chords in each axis direction instead of a single one,
! e.g. a 5x5 = 25 sheaf.
!-----------------------------------------------------------------------------------------
subroutine getNecroticVolume1(vol)
real(REAL_KIND) :: vol
real(REAL_KIND) :: rmin(3), rmax(3), cblob(3), cb(3), c(3), r, smin, smax, R3, Rn(3), Rdiv, Rsum
integer :: kcell, i, j, k, js, ix0, iy0, iz0
integer :: ixx, iyy, izz, axis, axis1, axis2, ixyz, rng(3,2), dx, dy, dz, dxm, dym, dzm, n
logical :: overlap, hit1, hit2
type(cell_type), pointer :: cp
real(REAL_KIND) :: fac = 1.2
integer :: rsheaf = 2

! First estimate blob centre location
rmin = 1.0e10
rmax = 0
do kcell = 1,nlist
	if (cell_list(kcell)%state == DEAD) cycle
	cp => cell_list(kcell)
	do k = 1,cp%nspheres
		rmin = min(rmin,cp%centre(:,k)-cp%radius(k))
		rmax = max(rmax,cp%centre(:,k)+cp%radius(k))
	enddo
enddo
cblob(:) = (rmin(:) + rmax(:))/2	! approximate blob centre
!write(*,'(a,3f8.4,a,3f8.4)') 'cblob: ',cblob,'  size: ',rmax-rmin
ix0 = cblob(1)/DELTA_X + 1
iy0 = cblob(2)/DELTA_X + 1
iz0 = cblob(3)/DELTA_X + 1
!write(*,*) 'getNecroticVolume:'
!write(*,'(a,3f8.4)') 'centre: ',cblob
!write(*,'(a,3i4)') 'centre gridsite: ',ix,iy,iz
! Check if any cells in this gridcell overlap with C(:)
! A cell is treated here as a cube with side = cell diameter
overlap = .false.
do k = 1,grid(ix0,iy0,iz0)%nc
	kcell = grid(ix0,iy0,iz0)%cell(k)
	cp => cell_list(kcell)
	do js = 1,cp%nspheres
		r = fac*cp%radius(js)
		c = cp%centre(:,js)
!		write(*,'(a,i6,3(2f8.4,2x))') 'kcell,bnds: ',kcell,c(1)-r,c(1)+r,c(2)-r,c(2)+r,c(3)-r,c(3)+r
		if ((cblob(1) >= c(1)-r .and. cblob(1) <= c(1)+r) .and. &
			(cblob(2) >= c(2)-r .and. cblob(2) <= c(2)+r) .and. &
			(cblob(3) >= c(3)-r .and. cblob(3) <= c(3)+r)) then
			overlap = .true.
			exit
		endif
	enddo
	if (overlap) exit
enddo
if (overlap) then
	vol = 0
	return
endif
! There is a gap in the centre, find the gap extent in each axis direction

do axis = 1,3
	rng(1,:) = [ix0, ix0]
	rng(2,:) = [iy0, iy0]
	rng(3,:) = [iz0, iz0]
	if (axis == 1) then
		rng(1,1) = 1
		rng(1,2) = NX-1
		axis1 = 2
		axis2 = 3
		dxm = 0
		dym = rsheaf
		dzm = rsheaf
	elseif (axis == 2) then
		rng(2,1) = 1
		rng(2,2) = NY-1
		axis1 = 1
		axis2 = 3
		dxm = rsheaf
		dym = 0
		dzm = rsheaf
	elseif (axis == 3) then
		rng(3,1) = 1
		rng(3,2) = NZ-1
		axis1 = 1
		axis2 = 2
		dxm = rsheaf
		dym = rsheaf
		dzm = 0
	endif
!	write(*,'(a,i2,3(2i3,2x))') 'axis: ',axis,((rng(i,j),j=1,2),i=1,3)
	! X direction.  
	!write(*,*) 'seeking rmin, rmax in X direction'
	
	Rsum = 0
	n = 0
	do dx = -dxm,dxm
	do dy = -dym,dym
	do dz = -dzm,dzm
	
	cb(1) = cblob(1) + 0.5*dx*Raverage
	cb(2) = cblob(2) + 0.5*dy*Raverage
	cb(3) = cblob(3) + 0.5*dz*Raverage

	smin = 1.0e10
	smax = -1.0e10
	hit1 = .false.
	hit2 = .false.
	do ixx = rng(1,1),rng(1,2)
	do iyy = rng(2,1),rng(2,2)
	do izz = rng(3,1),rng(3,2)
!		write(*,*) 'ixx,iyy,izz,nc: ',ixx,iyy,izz,grid(ixx,iyy,izz)%nc
		do k = 1,grid(ixx,iyy,izz)%nc
			kcell = grid(ixx,iyy,izz)%cell(k)
			cp => cell_list(kcell)
!			write(*,'(a,2i6,3f8.4)') 'k,kcell,cp%centre: ',k,kcell,cp%centre(:,1)
			do js = 1,cp%nspheres
!				r = fac*cp%radius(js)
				r = fac*Raverage
				c = cp%centre(:,js)
				if ((cb(axis1) > c(axis1)-r .and. cb(axis1) < c(axis1)+r) .and. &	! chord passes through the cube enclosing the cell
					(cb(axis2) > c(axis2)-r .and. cb(axis2) < c(axis2)+r)) then
	!				write(*,'(3(2f8.4,2x))') c(1)-r,c(1)+r,c(2)-r,c(2)+r,c(3)-r,c(3)+r
					if (smax < c(axis)+r .and. cb(axis) > c(axis)+r) then
						smax = c(axis)+r
	!					write(*,'(a,f8.4)') 'smax: ',smax
						hit1 = .true.
					endif
					if (smin > c(axis)-r .and. cb(axis) < c(axis)-r) then
						smin = c(axis)-r
	!					write(*,'(a,f8.4)') 'smin: ',smin
						hit2 = .true.
					endif
				endif
			enddo
		enddo
	enddo
	enddo
	enddo
	if (hit1 .and. hit2 .and. smax < smin) then
		n = n + 1
		Rsum = Rsum + smin - smax
!		write(*,'(a,2i4,3f8.4)') 'axis,n,smin,smax,Rsum: ',axis,n,smin,smax,Rsum
	endif
	
	enddo
	enddo
	enddo
	
	if (Rsum == 0) then
		vol = 0
		return
	endif
	Rn(axis) = Rsum/(2*n)
!	write(*,'(a,2i4,2f8.4)') 'axis done: ',axis,n,Rsum,Rn(axis)
enddo
	
R3 = (1./3.)*(Rn(1)**3 + Rn(2)**3 + Rn(3)**3)
Rdiv = (3*Vdivide0/(4*PI))**(1./3.)
!write(*,'(a,e12.3)') 'Rdiv: ',Rdiv
!write(*,'(a,3e12.3)') 'Rn: ',Rn
if (Rn(1) < 0 .or. Rn(2) < 0 .or. Rn(3) < 0) then
	write(logmsg,*) 'Rn < 0'
	call logger(logmsg)
	stop
endif
if (R3 > 0) then
	R = R3**(1./3)
	if (R < 4*Rdiv) then
		R3 = 0
		vol = 0
	endif
elseif (R3 == 0) then
	R = 0
else
	write(logmsg,*) 'Error: getNecroticVolume: R3 < 0: ',R3
	call logger(logmsg)
	stop
endif
vol = (4./3.)*PI*R3 
!write(*,'(a,e12.3)') 'vol: ',vol
end subroutine


!-----------------------------------------------------------------------------------------
! EC values should be average in the medium.
! DRUG_A = TPZ_DRUG e.g. SN30000
! DRUG_B = DNB_DRUG e.g. PR104A 
!-----------------------------------------------------------------------------------------
subroutine get_values(nvars,varID,ysim)
!DEC$ ATTRIBUTES DLLEXPORT :: get_values
integer :: nvars
character*(24) :: varID(nvars)
real(REAL_KIND) :: ysim(nvars)
integer :: ivar, ityp
integer :: Nviable(MAX_CELLTYPES)	!, Nlive(MAX_CELLTYPES)
real(REAL_KIND) :: plate_eff(MAX_CELLTYPES)

call MakeMediumCave
do ivar = 1,nvars
	if (varID(ivar) == 'OXYGEN_EC') then
		ysim(ivar) = chemo(OXYGEN)%medium_Cave
	elseif (varID(ivar) == 'GLUCOSE_EC') then
		ysim(ivar) = chemo(GLUCOSE)%medium_Cave
	elseif (varID(ivar) == 'LACTATE_EC') then
		ysim(ivar) = chemo(LACTATE)%medium_Cave
	elseif (varID(ivar) == 'SN30000_EC') then
		ysim(ivar) = chemo(DRUG_A)%medium_Cave
	elseif (varID(ivar) == 'SN30000_METAB1_EC') then
		ysim(ivar) = chemo(DRUG_A+1)%medium_Cave
	elseif (varID(ivar) == 'SN30000_METAB2_EC') then
		ysim(ivar) = chemo(DRUG_A+2)%medium_Cave
	elseif (varID(ivar) == 'PR104A_EC') then
		ysim(ivar) = chemo(DRUG_B)%medium_Cave
	elseif (varID(ivar) == 'PR104A_METAB1_EC') then
		ysim(ivar) = chemo(DRUG_B+1)%medium_Cave
	elseif (varID(ivar) == 'PR104A_METAB2_EC') then
		ysim(ivar) = chemo(DRUG_B+2)%medium_Cave
	elseif (varID(ivar) == 'OXYGEN_IC') then
!		ysim(ivar) = Caverage(MAX_CHEMO+OXYGEN)
		ysim(ivar) = 0
	elseif (varID(ivar) == 'GLUCOSE_IC') then
!		ysim(ivar) = Caverage(MAX_CHEMO+GLUCOSE)
		ysim(ivar) = 0
	elseif (varID(ivar) == 'LACTATE_IC') then
!		ysim(ivar) = Caverage(MAX_CHEMO+LACTATE)
		ysim(ivar) = 0
	elseif (varID(ivar) == 'SN30000_IC') then
!		ysim(ivar) = Caverage(MAX_CHEMO+DRUG_A)
		ysim(ivar) = 0
	elseif (varID(ivar) == 'SN30000_METAB1_IC') then
!		ysim(ivar) = Caverage(MAX_CHEMO+DRUG_A+1)
		ysim(ivar) = 0
	elseif (varID(ivar) == 'SN30000_METAB2_IC') then
!		ysim(ivar) = Caverage(MAX_CHEMO+DRUG_A+2)
		ysim(ivar) = 0
	elseif (varID(ivar) == 'PR104A_IC') then
!		ysim(ivar) = Caverage(MAX_CHEMO+DRUG_A)
		ysim(ivar) = 0
	elseif (varID(ivar) == 'PR104A_METAB1_IC') then
!		ysim(ivar) = Caverage(MAX_CHEMO+DRUG_A+1)
		ysim(ivar) = 0
	elseif (varID(ivar) == 'PR104A_METAB2_IC') then
!		ysim(ivar) = Caverage(MAX_CHEMO+DRUG_A+2)
		ysim(ivar) = 0
	elseif (varID(ivar) == 'NCELLS') then
		ysim(ivar) = Ncells			! for now, total live cells
	elseif (varID(ivar) == 'PE') then
		call getNviable	!(Nviable, Nlive)
		do ityp = 1,Ncelltypes
!			if (Nlive(ityp) > 0) then
!				plate_eff(ityp) = real(Nviable(ityp))/Nlive(ityp)
			if (Ncells_type(ityp) > 0) then
				plate_eff(ityp) = real(Nviable(ityp))/Ncells_type(ityp)
			else
				plate_eff(ityp) = 0
			endif
		enddo
		ysim(ivar) = plate_eff(1)	! for now, just type 1 cells
	elseif (varID(ivar) == 'RADIATION') then
		ysim(ivar) = -1
	else
		write(*,*) 'varID is not in the list of possible IDs: ',varID(ivar)
		stop
	endif
enddo
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine MakeMediumCave
integer :: ichemo

do ichemo = 1,MAX_CHEMO
	if (chemo(ichemo)%used) then
		chemo(ichemo)%medium_Cave = sum(chemo(ichemo)%Cave_b)/(NXB*NYB*NZB)
	endif
enddo
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine get_summary(summaryData,i_hypoxia_cutoff,i_growth_cutoff) BIND(C)
!DEC$ ATTRIBUTES DLLEXPORT :: get_summary
use, intrinsic :: iso_c_binding
real(c_double) :: summaryData(*)
integer(c_int) :: i_hypoxia_cutoff,i_growth_cutoff
!integer :: Nviable(MAX_CELLTYPES)
integer :: nhypoxic(3), nclonohypoxic(3), ngrowth(3), nogrow(MAX_CELLTYPES), nphase(6+1), nmutations, nclono, &
    medium_oxygen, medium_glucose, medium_lactate, medium_drug(2,0:2), &
    bdry_oxygen, bdry_glucose, bdry_lactate, bdry_drug(2,0:2)
integer :: TNradiation_dead, TNdrug_dead(2),  TNdead, TNviable, TNnonviable, TNATP_dead, TNnogrow, &
           Ntagged_ATP(MAX_CELLTYPES), Ntagged_radiation(MAX_CELLTYPES), Ntagged_drug(2,MAX_CELLTYPES), &
           TNtagged_ATP, TNtagged_radiation, TNtagged_drug(2)
integer :: ityp, i, im, idrug
real(REAL_KIND) :: diam_um,  npmm3, Tplate_eff
real(REAL_KIND) :: hypoxic_fraction, clonohypoxic_fraction(3), growth_fraction, nogrow_fraction, viable_fraction, clono_fraction, phase_fraction(6+1), rmutations
real(REAL_KIND) :: diam_cm, vol_cm3, vol_mm3, hour, necrotic_fraction, doubling_time, plate_eff(MAX_CELLTYPES)
real(REAL_KIND) :: volume_cm3(5), maxarea(5), diameter_um(5)
real(REAL_KIND) :: cmedium(MAX_CHEMO), cbdry(MAX_CHEMO)
real(REAL_KIND) :: r_G, r_P, r_A, r_I, P_utilisation

hour = istep*DELTA_T/3600.
call getDiamVol(diam_cm,vol_cm3)
vol_mm3 = vol_cm3*1000				! volume in mm^3
diam_um = diam_cm*10000
npmm3 = Ncells/vol_mm3

Ntagged_ATP(:) = NATP_tag(:)				! number currently dying of ATP starvation
Ntagged_radiation(:) = Nradiation_tag(:)	! number currently tagged by radiation
Ntagged_drug(1,:) = Ndrug_tag(1,:)			! number currently tagged by drugA
Ntagged_drug(2,:) = Ndrug_tag(2,:)			! number currently tagged by drugA

TNtagged_ATP = sum(Ntagged_ATP(1:Ncelltypes))
TNtagged_radiation = sum(Ntagged_radiation(1:Ncelltypes))
TNtagged_drug(1) = sum(Ntagged_drug(1,1:Ncelltypes))
TNtagged_drug(2) = sum(Ntagged_drug(2,1:Ncelltypes))

TNATP_dead = sum(NATP_dead(1:Ncelltypes))
TNradiation_dead = sum(Nradiation_dead(1:Ncelltypes))
TNdrug_dead(1) = sum(Ndrug_dead(1,1:Ncelltypes))
TNdrug_dead(2) = sum(Ndrug_dead(2,1:Ncelltypes))

call getNviable	!(Nviable, Nlive)
TNviable = sum(Nviable(1:Ncelltypes))
TNnonviable = sum(Ndying(1:Ncelltypes))
TNdead = sum(Ndead(1:Ncelltypes))
viable_fraction = TNviable/real(Ncells)

call getHypoxicCount(nhypoxic)
hypoxic_fraction = nhypoxic(i_hypoxia_cutoff)/real(Ncells)
call getClonoHypoxicCount(nclonohypoxic)
if (TNviable > 0) then
	clonohypoxic_fraction = nclonohypoxic(:)/real(TNviable)
else
	clonohypoxic_fraction = 0
endif		
!call getGrowthCount(ngrowth, nogrow, nphase)
call getGrowthCount(ngrowth,nogrow,nphase,nmutations,nclono)
TNnogrow = sum(nogrow(:))
nogrow_fraction = TNnogrow/real(TNviable)
growth_fraction = ngrowth(i_growth_cutoff)/real(Ncells)
phase_fraction = nphase/real(Ncells)
clono_fraction = nclono/real(Ncells)
rmutations = (nmutations*1000000.)/Ncells
!call getNecroticFraction(necrotic_fraction,vol_cm3)
!necrotic_percent = 100.*necrotic_fraction
do ityp = 1,Ncelltypes
!	if (Nlive(ityp) > 0) then
!		plate_eff(ityp) = real(Nviable(ityp))/Nlive(ityp)
	if (Ncells_type(ityp) > 0) then
		plate_eff(ityp) = real(Nviable(ityp))/Ncells_type(ityp)
	else
		plate_eff(ityp) = 0
	endif
enddo
plate_eff = 100.*plate_eff
Tplate_eff = 0
do ityp = 1,Ncelltypes
	Tplate_eff = Tplate_eff + plate_eff(ityp)*celltype_fraction(ityp)
enddo

if (use_metabolism) then
	! Get normalised metabolism state variables, averaged over cells
	call AverageCellMetabolism(r_G,r_P,r_A,r_I,P_utilisation)
else
	r_G = 0
	r_P = 0
	r_A = 0
	r_I = 0
	P_utilisation = 0
endif

call getMediumConc(cmedium, cbdry)

call getVolumes(volume_cm3,maxarea)
do i = 1,5
	diameter_um(i) = 10000*(6*volume_cm3(i)/PI)**(1./3)
enddo
necrotic_fraction = volume_cm3(1)/volume_cm3(5)

if (ndoublings > 0) then
    doubling_time = doubling_time_sum/(3600*ndoublings)
else
    doubling_time = 0
endif

summaryData(1:57) = [ rint(istep), rint(Ncells), rint(TNviable), rint(TNnonviable), &
	rint(TNATP_dead), rint(TNdrug_dead(1)), rint(TNdrug_dead(2)), rint(TNradiation_dead), rint(TNdead), &
    rint(TNtagged_ATP), rint(TNtagged_drug(1)), rint(TNtagged_drug(2)), rint(TNtagged_radiation), &
	diam_um, vol_mm3, &
	100*viable_fraction, 100*hypoxic_fraction, 100*clonohypoxic_fraction(i_hypoxia_cutoff), &
	100*growth_fraction, 100*nogrow_fraction, 100*necrotic_fraction, 100*clono_fraction, &
	Tplate_eff, npmm3, &
	cmedium(OXYGEN), cmedium(GLUCOSE), cmedium(LACTATE), cmedium(DRUG_A:DRUG_A+2), cmedium(DRUG_B:DRUG_B+2), &
	cbdry(OXYGEN), cbdry(GLUCOSE), cbdry(LACTATE), cbdry(DRUG_A:DRUG_A+2), cbdry(DRUG_B:DRUG_B+2), &
	doubling_time, r_G, r_P, r_A, r_I, rint(ndoublings), P_utilisation, &
	100*phase_fraction(1:7), rmutations]
write(nfres,'(a,a,2a12,i8,7e12.4,21i7,42e12.4)') trim(header),' ',gui_run_version, dll_run_version, &
	istep, hour, vol_mm3, diameter_um, &
	Ncells_type(1:2), TNviable, TNnonviable, &
	NATP_dead(1:2), Ndrug_dead(1,1:2), Ndrug_dead(2,1:2), Nradiation_dead(1:2), TNdead, &
    Ntagged_ATP(1:2), Ntagged_drug(1,1:2), Ntagged_drug(2,1:2), Ntagged_radiation(1:2), &
	viable_fraction, hypoxic_fraction, clonohypoxic_fraction(i_hypoxia_cutoff), &
	growth_fraction, nogrow_fraction, necrotic_fraction, &
	clono_fraction, plate_eff(1:2), &
	cmedium(OXYGEN), cmedium(GLUCOSE), cmedium(LACTATE), cmedium(DRUG_A:DRUG_A+2), cmedium(DRUG_B:DRUG_B+2), &
	cbdry(OXYGEN), cbdry(GLUCOSE), cbdry(LACTATE), cbdry(DRUG_A:DRUG_A+2), cbdry(DRUG_B:DRUG_B+2), &
	doubling_time, r_G, r_P, r_A, r_I, real(ndoublings), P_utilisation, &
	100*phase_fraction(1:7), rmutations
		
call sum_dMdt(GLUCOSE)

if (diam_count_limit > LIMIT_THRESHOLD) then
	if (Ncells > diam_count_limit) limit_stop = .true.
elseif (diam_count_limit > 0) then
	if (diam_um > diam_count_limit) limit_stop = .true.
endif

ndoublings = 0
doubling_time_sum = 0

end subroutine

!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------
subroutine getNviable	!(Nviable, Nlive)
!integer :: Nviable(:) 
integer :: Nlive(MAX_CELLTYPES)
integer :: kcell, ityp, idrug
logical :: tag
type(cell_type), pointer :: cp

Nviable = 0
Nlive = 0
do kcell = 1,nlist
	cp => cell_list(kcell)
	if (cp%state == DEAD) cycle
    ityp = cp%celltype
    Nlive(ityp) = Nlive(ityp) + 1
!	if (cp%anoxia_tag .or. cp%aglucosia_tag .or. cp%radiation_tag .or. cp%state == DYING) cycle
!    tag = .false.
!    do idrug = 1,ndrugs_used
!		if (cell_list(kcell)%drug_tag(idrug)) tag = .true.
!	enddo
!	if (tag) cycle
	if (cp%state == DYING) cycle
	Nviable(ityp) = Nviable(ityp) + 1
enddo
if (Nlive(1) /= Ncells_type(1)) then
	write(*,*) 'Error: getNviable: Nlive /= Ncells_type: ',Nlive(1),Ncells_type(1)
	stop
endif
if (Nviable(1) /= Ncells_type(1) - Ndying(1)) then
	write(*,'(a,3i8)') 'Error: getNviable: Nviable /= Ncells_type - Ndying: ',Nviable(1),Ncells_type(1),Ndying(1)
	stop
endif
end subroutine

!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------
subroutine set_colourscheme(icset) bind(C)
!DEC$ ATTRIBUTES DLLEXPORT :: set_colourscheme
use, intrinsic :: iso_c_binding
integer(c_int),value :: icset

icolourscheme = icset
end subroutine

end module

#if 0
!--------------------------------------------------------------------------------
! Need to compare growth rate with a fraction of average growth rate
!--------------------------------------------------------------------------------
subroutine getGrowthCount(ngrowth)
integer :: ngrowth(3)
integer :: kcell, i, ityp
real(REAL_KIND) :: r_mean(2)
type(cell_type), pointer :: cp

r_mean(1:2) = Vdivide0/(2*divide_time_mean(1:2))
ngrowth = 0
do kcell = 1,nlist
    cp => cell_list(kcell)
	if (cp%state == DEAD) cycle
	if (cp%dVdt == 0) cycle
	ityp = cp%celltype
	do i = 1,3
		if (cp%dVdt < growthcutoff(i)*r_mean(ityp)) then
		    ngrowth(i) = ngrowth(i) + 1
		endif
	enddo
enddo
end subroutine

#endif