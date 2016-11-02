! Motion based on forces
module fmotion

use global
implicit none

real(REAL_KIND) :: fwall

contains

!-----------------------------------------------------------------------------------------
! Start with a very simple dropping, just place the blob slightly above the bottom of the well.
!-----------------------------------------------------------------------------------------
subroutine drop_blob
integer :: kcell, isph, kmin
real(REAL_KIND) :: blob_zmin
type(cell_type), pointer :: cp

blob_zmin = 1.0e10
do kcell = 1,nlist
	cp => cell_list(kcell)
	if (cp%state == DEAD) cycle
	do isph = 1,cp%nspheres
		if (cp%centre(3,isph)-cp%radius(isph) < blob_zmin) then
			blob_zmin = cp%centre(3,isph)-cp%radius(isph)
			kmin = kcell
		endif
	enddo
enddo
write(nflog,*) 'blob_zmin: ',blob_zmin
write(nflog,*) 'min x,y: ',cell_list(kmin)%centre(1:2,1)
do kcell = 1,nlist
	cp => cell_list(kcell)
	if (cp%state == DEAD) cycle
	do isph = 1,cp%nspheres
		cp%centre(3,isph) = cp%centre(3,isph) - blob_zmin
	enddo
enddo
end subroutine

!-----------------------------------------------------------------------------------------
! In each time step we want to find cell positions that make all forces = 0, or
! that minimises sum of squares of forces.  We want to minimise the number of force
! computations needed.
! The time step size dt_move is set dynamically here.
!-----------------------------------------------------------------------------------------
subroutine fmover(dt_move, done, ok)
real(REAL_KIND) :: dt_move
logical :: done, ok
integer :: k1, kcell, kpar, nd, nr, nc(0:8), kfrom(0:8), kto(0:8), tMnodes
real(REAL_KIND), allocatable :: force(:,:,:)
real(REAL_KIND) :: fmax, dx1(3), dx2(3)
type(cell_type), pointer :: cp1
logical :: penetrated

done = .false.
ok = .true.
dt_move = ndt*delta_tmove
if (ncells <= 1 .and. cell_list(1)%Iphase) then
	dt_move = DELTA_T	! one time step of DELTA_T
	done = .true.
	return
endif
allocate(force(3,ncells,2))
call forces(force,fmax,penetrated,ok)
if (.not.ok) return

if (fmax > 0) then
do
	if (dt_move*fmax/kdrag > delta_max .or. penetrated) then
		ndt = ndt/2
!		write(*,'(a,5e12.3,i4)') 'fmover: -ndt: ',dt_move,fmax,kdrag,dt_move*fmax/kdrag,delta_min,ndt
		if (ndt <= 1) then
			ndt = 1
			dt_move = ndt*delta_tmove
			exit
		endif
!		cycle
	elseif (dt_move*fmax/kdrag < delta_min) then
		ndt = ndt + 1
!		write(*,'(a,5e12.3,i4)') 'fmover: +ndt: ',dt_move,fmax,kdrag,dt_move*fmax/kdrag,delta_min,ndt
		if (ndt >= ndt_max) then
			ndt = ndt_max
			dt_move = ndt*delta_tmove
			exit
		endif
!		cycle
	endif
	exit
enddo
endif
if (dt_move + t_fmover > DELTA_T) then
	dt_move = DELTA_T - t_fmover
	done = .true.
endif

if (ncells < 500) then
	tMnodes = 1
else
	tMnodes = Mnodes
endif
if (ncells <= tMnodes) then
    do kpar = 0,tMnodes-1
        if (kpar < ncells) then
            kfrom(kpar) = kpar+1
            kto(kpar) = kpar+1
        else
            kfrom(kpar) = 1
            kto(kpar) = 0
        endif
    enddo
else
    nd = ncells/tMnodes
    nr = ncells - nd*tMnodes
    do kpar = 0,tMnodes-1
        kfrom(kpar) = kpar*nd + 1
        kto(kpar) = (kpar+1)*nd
        if (kpar == tMnodes-1) kto(kpar) = ncells
    enddo
endif

do kpar = 0,tMnodes-1
    do k1 = kfrom(kpar),kto(kpar)
!do k1 = 1,ncells
	kcell = perm_index(k1)
	cp1 => cell_list(kcell)
	if (cp1%Iphase) then
		dx1 = dt_move*force(:,k1,1)/kdrag
		cp1%centre(:,1) = cp1%centre(:,1) + dx1
	else
		dx1 = dt_move*force(:,k1,1)/kdrag
		dx2 = dt_move*force(:,k1,2)/kdrag
		cp1%centre(:,1) = cp1%centre(:,1) + dx1
		cp1%centre(:,2) = cp1%centre(:,2) + dx2
	endif
enddo
enddo
!omp end parallel do
deallocate(force)
ok = .true.
end subroutine

!-----------------------------------------------------------------------------------------
! Trying a modification:
! all cells move through dt - can this work?
!-----------------------------------------------------------------------------------------
subroutine tried_fmover(dt, done, ok)
real(REAL_KIND) :: dt
logical :: done, ok
integer :: k1, kcell, kpar, nd, nr, nc(0:8), kfrom(0:8), kto(0:8), tMnodes, cell_ndt
real(REAL_KIND), allocatable :: force(:,:,:)
real(REAL_KIND) :: fmax, dx1(3), dx2(3), dt_move, t_moved
type(cell_type), pointer :: cp1
real(REAL_KIND),allocatable :: cell_fmax(:)

!write(*,*) 'fmover: '
done = .false.
ok = .true.
!dt_move = ndt*delta_tmove
if (ncells <= 1 .and. cell_list(1)%Iphase) then
!	dt_move = DELTA_T	! one time step of DELTA_T
	done = .true.
	return
endif
allocate(force(3,ncells,2))
allocate(cell_fmax(ncells))
cell_fmax = 0
!call forces(force,fmax,cell_fmax,ok)
if (.not.ok) return

if (fmax == 0) return

!if (ncells < 500) then
!	tMnodes = 1
!else
!	tMnodes = Mnodes
!endif
!if (ncells <= tMnodes) then
!    do kpar = 0,tMnodes-1
!        if (kpar < ncells) then
!            kfrom(kpar) = kpar+1
!            kto(kpar) = kpar+1
!        else
!            kfrom(kpar) = 1
!            kto(kpar) = 0
!        endif
!    enddo
!else
!    nd = ncells/tMnodes
!    nr = ncells - nd*tMnodes
!    do kpar = 0,tMnodes-1
!        kfrom(kpar) = kpar*nd + 1
!        kto(kpar) = (kpar+1)*nd
!        if (kpar == tMnodes-1) kto(kpar) = ncells
!    enddo
!endif

!do kpar = 0,tMnodes-1
!    do k1 = kfrom(kpar),kto(kpar)
do k1 = 1,ncells
	kcell = perm_index(k1)
	cp1 => cell_list(kcell)
	cell_ndt = cp1%ndt + 1      ! now ndt is by cell
    dt_move = cell_ndt*delta_tmove
    t_moved = 0
    done = .false.
    do while (.not.done)
	    if (dt_move*cell_fmax(kcell)/kdrag > delta_max) then
		    cell_ndt = cell_ndt/2
		    if (cell_ndt <= 1) then
			    cell_ndt = 1
			    dt_move = cell_ndt*delta_tmove
    !			exit
		    endif
	    elseif (dt_move*cell_fmax(kcell)/kdrag < delta_min) then
		    cell_ndt = cell_ndt + 1
		    if (cell_ndt >= ndt_max) then
			    cell_ndt = ndt_max
			    dt_move = cell_ndt*delta_tmove
    !			exit
		    endif
	    endif
        if (dt_move + t_moved > dt) then
	        dt_move = dt - t_moved
	        done = .true.
        endif
        cp1%ndt = cell_ndt

	    if (cp1%Iphase) then
		    dx1 = dt_move*force(:,k1,1)/kdrag
		    cp1%centre(:,1) = cp1%centre(:,1) + dx1
	    else
		    dx1 = dt_move*force(:,k1,1)/kdrag
		    dx2 = dt_move*force(:,k1,2)/kdrag
		    cp1%centre(:,1) = cp1%centre(:,1) + dx1
		    cp1%centre(:,2) = cp1%centre(:,2) + dx2
	    endif
	    t_moved = t_moved + dt_move
    enddo
enddo

!!omp end parallel do
deallocate(force)
deallocate(cell_fmax)
ok = .true.
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine forces(force,fmax,penetrated,ok)
real(REAL_KIND) :: fmax, force(:,:,:)
logical :: penetrated,ok
integer :: k1, kcell, kpar, nd, nr, nc(0:8), kfrom(0:8), kto(0:8), tMnodes
real(REAL_KIND) :: F(3,2), r(3), amp, fsum, fwall_prev
real(REAL_KIND),allocatable :: cell_fmax(:)
logical :: badforce, pen

if (ncells < 100) then
	tMnodes = 1
else
	tMnodes = Mnodes
endif
if (ncells <= tMnodes) then
    do kpar = 0,tMnodes-1
        if (kpar < ncells) then
            kfrom(kpar) = kpar+1
            kto(kpar) = kpar+1
        else
            kfrom(kpar) = 1
            kto(kpar) = 0
        endif
    enddo
else
    nd = ncells/tMnodes
    nr = ncells - nd*tMnodes
    do kpar = 0,tMnodes-1
 !       nc(kpar) = nd
 !       if (kpar == tMnodes-1) nc(kpar) = nd + nr
        kfrom(kpar) = kpar*nd + 1
        kto(kpar) = (kpar+1)*nd
        if (kpar == tMnodes-1) kto(kpar) = ncells
    enddo
endif
!write(*,'(15i4)') ncells,Mnodes,kfrom(0:Mnodes-1),kto(0:Mnodes-1)
allocate(cell_fmax(ncells))
cell_fmax = 0
fsum = 0
fwall_prev = fwall
fwall = 0
badforce = .false.
call omp_set_num_threads(tMnodes)
penetrated = .false.
!$omp parallel do private(k1,kcell,F,amp,ok)
!do k1 = 1,ncells
!	write(*,*) 'Threads, max: ',omp_get_num_threads(),omp_get_max_threads()
do kpar = 0,tMnodes-1
    do k1 = kfrom(kpar),kto(kpar)
	    kcell = perm_index(k1)
	    if (kcell == 0) then
	        write(*,*) 'forces: kpar,k1,kcell: ',kpar,k1,kcell
	        stop
	    endif
	    call get_cell_force(kcell,F,pen,ok)
!	    if (kcell == 75) write(*,'(a,3e12.3)') 'F: (a): ',F(:,1)
	    if (.not.ok) then
			badforce = .true.
			exit
		endif
		if (pen) then
		    penetrated = .true.
		endif
	    call get_random_dr(r)
	    F(:,1) = F(:,1) + frandom*r + fwall_prev*[0,0,1]
!	    if (kcell == 1) write(*,*) 'F: (b): ',F(:,1)
	    force(:,k1,1) = F(:,1)
	    amp = sqrt(dot_product(F(:,1),F(:,1)))
	    cell_fmax(k1) = max(cell_fmax(k1),amp)
!	    if (amp > 0) write(*,*) 'amp: ',k1,amp
	    if (.not.cell_list(kcell)%Iphase) then
			call get_random_dr(r)
			F(:,2) = F(:,2) + frandom*r + fwall_prev*[0,0,1]
		    force(:,k1,2) = F(:,2)
		    cell_fmax(k1) = max(cell_fmax(k1),sqrt(dot_product(F(:,2),F(:,2))))
        endif
    enddo
enddo
!omp end parallel do
fmax = maxval(cell_fmax(:))
!if (ncells > 50e10 ) then
!	write(*,*) 'fsum: ',fsum
!	write(*,'(5e12.3)') cell_fmax
!	if (fmax == 0) then
!		write(*,*) 'fmax = 0'
!		stop
!	endif
!endif
deallocate(cell_fmax)
ok = .not.badforce
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine get_cell_force(kcell,F,penetrated,ok)
integer :: kcell
real(REAL_KIND) :: F(3,2)
logical :: penetrated, ok
integer :: k2, knbr	
integer :: isphere1, isphere2, nspheres1, nspheres2
real(REAL_KIND) :: R1, c1(3), R2, c2(3), v(3), d, d_hat, dF
logical :: incontact, Mphase
type(cell_type), pointer :: cp1, cp2
!integer :: nvar	! nv = ncells + # of Mphase cells

ok = .true.
F = 0
cp1 => cell_list(kcell)
!write(*,*) 'get_cell_force: ',kcell,' nbrs: ',cp1%nbrs
Mphase = .not.cp1%Iphase
if (Mphase) then
	nspheres1 = 2
else
	nspheres1 = 1
endif
F = 0
do k2 = 1,cp1%nbrs
	knbr = cp1%nbrlist(k2)%indx
	cp2 => cell_list(knbr)
	if (cp2%Iphase) then
		nspheres2 = 1
	else
		nspheres2 = 2
	endif
	do isphere1 = 1,nspheres1
		R1 = cp1%radius(isphere1)
		c1 = cp1%centre(:,isphere1)
!		if (isphere1 == 1) then
!			iv = cp1%varindex
!		else
!			iv = cp1%varindex + 1
!		endif
		do isphere2 = 1,nspheres2
			R2 = cp2%radius(isphere2)
			c2 = cp2%centre(:,isphere2)
			v = c1 - c2
			d = sqrt(dot_product(v,v))
			v = v/d
			incontact = cp1%nbrlist(k2)%contact(isphere1,isphere2)
			dF = get_force(kcell,R1,R2,d,incontact,penetrated,ok)	! returns magnitude and sign of force 
!			if (kcell == 75) then
!				write(*,'(a,3i6,e12.3)') 'dF: ',knbr,isphere1,isphere2,dF
!				write(*,'(6e12.3)') cp1%centre,cp2%centre
!			endif
			if (isnan(dF)) then
				write(*,'(a,2i4,4e12.3)') 'dF: ',kcell,knbr,R1,R2,d,dF
				stop
			endif
!			if (dF /= 0) then
!				write(*,'(a,2i4,4e12.3)') 'dF: ',kcell,knbr,R1,R2,d,dF
!				stop
!			endif
            if (.not.ok) then
                write(nflog,*) kcell,knbr,isphere1,isphere2
                write(nflog,'(a,i6,3e12.3)') 'nbrlist: ',kcell,c1
                write(nflog,'(10i6)') cp1%nbrlist(1:cp1%nbrs)%indx
                write(nflog,'(a,i6,3e12.3)') 'nbrlist: ',knbr,c2
                write(nflog,'(10i6)') cp2%nbrlist(1:cp2%nbrs)%indx
!                write(*,'(a,3f8.3)') 'cell 10: ',cell_list(10)%centre(:,1)
!                write(*,'(a,3f8.3)') 'cell 30: ',cell_list(30)%centre(:,1)
!                write(*,'(a,3f8.3)') 'cell 44: ',cell_list(44)%centre(:,1)
                return
            endif
			if (isphere1 == 1) F(:,isphere1) = F(:,isphere1) + dF*v
			if (Mphase .and. isphere1 == 2) F(:,isphere1) = F(:,isphere1) + dF*v
		enddo
	enddo
enddo
if (Mphase) then	! compute force between separating spheres, based on deviation from desired d = cp1%d
	c1 = cp1%centre(:,1)
	c2 = cp1%centre(:,2)
	v = c1 - c2
	d = sqrt(dot_product(v,v))
	v = v/d
	d_hat = max(cp1%mitosis*cp1%d_divide,small_d)
	dF = get_separating_force(d,d_hat)
	F(:,1) = F(:,1) + dF*v	! dF acts in the direction of v on sphere #1 when dF > 0
	F(:,2) = F(:,2) - dF*v	! dF acts in the direction of -v on sphere #2 when dF > 0 
endif

! Check for wall force (bottom of the well) NOT USED
! For now the bottom is a plane
! Treat wall as another cell
do isphere1 = 1,nspheres1
	R1 = cp1%radius(isphere1)
	c1 = cp1%centre(:,isphere1)
	d = c1(3)
	incontact = (d < R1)
	v = [0, 0, 1]
	if (R1 < 1.5*d) then
		dF = get_force(kcell,R1,R1,2*d,incontact,penetrated,ok)	! returns magnitude and sign of force 
		if (dF > 0) then
			fwall = max(fwall,dF)
		endif
!		if (dF < 0) then
!			dF = wall_attraction_factor*dF
!		endif
!		if (isphere1 == 1) F(:,isphere1) = F(:,isphere1) + dF*v
!		if (Mphase .and. isphere1 == 2) F(:,isphere1) = F(:,isphere1) + dF*v
	endif
enddo
end subroutine
!-----------------------------------------------------------------------------------------
! Each cell has a varindex, and for Mphase the second sphere is varindex+1
!-----------------------------------------------------------------------------------------
subroutine forces1(force,fmax,ok)
real(REAL_KIND) :: fmax, force(:,:,:)
logical :: ok
integer :: k1, k2, kcell, knbr, iv(2)
integer :: isphere1, isphere2, nspheres1, nspheres2
real(REAL_KIND) :: R1, c1(3), R2, c2(3), v(3), d, d_hat, dF, F(3,2)
logical :: incontact, Mphase, use_force, penetrated
type(cell_type), pointer :: cp1, cp2
integer :: nvar	! nv = ncells + # of Mphase cells

ok = .true.
if (ncells == 1 .and. cell_list(1)%Iphase) then
	force(:,1,:) = 0
	return
endif
!write(logmsg,*) 'forces: ',istep,ncells,ncells_mphase
!call logger(logmsg)
use_force = .true.
fmax = 0
do k1 = 1,ncells
	kcell = perm_index(k1)
	cp1 => cell_list(kcell)
	Mphase = .not.cp1%Iphase
	if (cp1%Iphase) then
		nspheres1 = 1
!		iv(1) = cp1%varindex
	else
		nspheres1 = 2
!		iv(1) = cp1%varindex
!		iv(2) = cp1%varindex + 1
	endif
	F = 0
	if (use_force) then
	do k2 = 1,cp1%nbrs
		knbr = cp1%nbrlist(k2)%indx
		cp2 => cell_list(knbr)
		if (cp2%Iphase) then
			nspheres2 = 1
		else
			nspheres2 = 2
		endif
		do isphere1 = 1,nspheres1
			R1 = cp1%radius(isphere1)
			c1 = cp1%centre(:,isphere1)
			!if (isphere1 == 1) then
			!	iv = cp1%varindex
			!else
			!	iv = cp1%varindex + 1
			!endif
			do isphere2 = 1,nspheres2
				R2 = cp2%radius(isphere2)
				c2 = cp2%centre(:,isphere2)
				v = c1 - c2
				d = sqrt(dot_product(v,v))
				v = v/d
				incontact = cp1%nbrlist(k2)%contact(isphere1,isphere2)
				dF = get_force(kcell,R1,R2,d,incontact,penetrated,ok)	! returns magnitude and sign of force
                if (.not.ok) then
                    stop
                endif
				if (isphere1 == 1) F(:,isphere1) = F(:,isphere1) + dF*v
				if (Mphase .and. isphere1 == 2) F(:,isphere1) = F(:,isphere1) + dF*v
			enddo
		enddo
	enddo
	endif
	if (Mphase) then	! compute force between separating spheres, based on deviation from desired d = cp1%d
		c1 = cp1%centre(:,1)
		c2 = cp1%centre(:,2)
		v = c1 - c2
		d = sqrt(dot_product(v,v))
		v = v/d
		d_hat = max(cp1%mitosis*cp1%d_divide,small_d)
		dF = get_separating_force(d,d_hat)
!		write(*,'(a,3e12.3)') 'dF: ',d,d_hat,dF
		F(:,1) = F(:,1) + dF*v	! dF acts in the direction of v on sphere #1 when dF > 0
		F(:,2) = F(:,2) - dF*v	! dF acts in the direction of -v on sphere #2 when dF > 0
	endif
	force(:,k1,1) = F(:,1)
	fmax = max(fmax,sqrt(dot_product(F(:,1),F(:,1))))
	if (Mphase) then
		force(:,k1,2) = F(:,2)
!		write(*,'(i4,3f8.3,4x,3f8.3)') kcell,F(1,:),F(2,:)
		fmax = max(fmax,sqrt(dot_product(F(:,2),F(:,2))))
	endif
	
enddo

end subroutine

!-----------------------------------------------------------------------------------------
! This force combines cell-cell repulsion and cell-cell attraction, which depends upon
! the contact history, is use_hysteresis = .true.
! F > 0 = repulsion
! F < 0 = attraction
!-----------------------------------------------------------------------------------------
function get_force(kcell,R1,R2,d,incontact,penetrated,ok) result(F)
integer :: kcell
real(REAL_KIND) :: R1, R2, d
logical :: incontact, penetrated, ok
real(REAL_KIND) :: F
real(REAL_KIND) :: x, dx0, e
logical :: allow_penetration = .true.

ok = .true.
penetrated = .false.
if (allow_penetration) then
    dx0 = 0.01*x0_force 
    e = (x0_force + dx0)*(a_force/(dx0*(x1_force - x0_force + dx0)) + b_force)
else
    dx0 = 0
endif
x = d/(R1+R2)
!write(*,'(a,4e12.3)') 'get_force: ',R1,R2,d,x
if (use_hysteresis) then
	if ((incontact .and. x > xcross2_force) .or. &
		(.not.incontact .and. x > xcross1_force)) then
		F = 0
		return
	endif
	if (x < x0_force) then
		write(logmsg,'(a,i6,3e12.3,2f8.3,a,L2)') 'Error: get_force: x < x0: ',kcell,R1,R2,d,x,x0_force,' incontact: ',incontact
		call logger(logmsg)
		F = 0
		ok = .false.
		return
	!	stop
	else
		F = a_force/((x-x0_force)*(x1_force-x)) + b_force
	endif
else
	if (x > xcross2_force) then
		F = 0
		return
	endif
	if (x < x0_force + dx0) then
	    if (allow_penetration) then
	        F = e/x
!	        penetrated = .true.     ! this may not be necessary
!	        write(*,'(a,4e12.3)') 'penetrated: ',x,x0_force + dx0,e,F
	    else
    		write(logmsg,'(a,i6,3e12.3,2f8.3,a,L2)') 'Error: get_force: x < x0: ',kcell,R1,R2,d,x,x0_force,' incontact: ',incontact
	    	call logger(logmsg)
		    F = 0
	    	ok = .false.
		    return
		endif
	!	stop
	else
		F = a_force/((x-x0_force)*(x1_force-x)) + b_force
	endif
endif
end function

!-----------------------------------------------------------------------------------------
! This force tries to maintain the desired separation distance d^ of the two parts of a
! dividing cell.  
! The force F > 0, tending to push the two "spheres" apart, when d < d^, and
! the force F < 0, tending to push them together, when d > d^.
! The original setup used distances in um, so to keep forces the same convert separation
! distance to um.
!-----------------------------------------------------------------------------------------
function get_separating_force(d,d_hat) result(F)
real(REAL_KIND) :: d, d_hat
real(REAL_KIND) :: F
real(REAL_KIND) :: dum

dum = (d_hat - d)*1.0e4		! cm -> um
F = a_separation*dum**3
end function

!-----------------------------------------------------------------------------------------
! The force is expressed as a function of x = d/(R1+R2)
! where d is the centre-centre distance, R1 and R2 are the two cell radii.
! The cell-cell force parameters are:
! a_force		scales the basic force function 1/((x-x0)(x1-x))
! c_force		magnitude of maximum attractive force, i.e. minimum of function at x = (x0 + x1)/2
! x0_force		location of lower asymptote
! x1_force		location of upper asymptote
! xhat_force	location of zero crossing
!
! F(x) = a/((x-x0)(x1-x)) - c - 4a/(x1-x0)^2
! the minimum value -c occurs at x = (x0+x1)/2
! the relationship between the parameters is:
! xhat = (x0+x1)/2 + (x1-x0)/2.(c(x1-x0)^2/(c(x1-x0)^2 + 4a))^(1/2)
!-----------------------------------------------------------------------------------------
subroutine setup_force_parameters
integer :: i
real(REAL_KIND) :: dx, delta, x, F

!a_separation = 1
!a_force = 1.0
!c_force = 0.25	! 1.0
!x0_force = 0.80
!x1_force = 1.3
write(logmsg,'(a,f6.2)') 'force parameters: a_separation: ',a_separation
call logger(logmsg)

dx = x1_force - x0_force
b_force = -c_force - 4*a_force/dx**2
delta = dx**2 + 4*a_force/b_force
xcross1_force = (x0_force+x1_force)/2 - 0.5*sqrt(delta)
xcross2_force = (x0_force+x1_force)/2 + 0.5*sqrt(delta)
write(logmsg,'(a,4f6.2)') 'force parameters: a_force,b_force,x0_force,x1_force: ',a_force,b_force,x0_force,x1_force
call logger(logmsg)
write(logmsg,*) 'force parameters: xcross:',xcross2_force
call logger(logmsg)

!dt_min = 1.0
!delta_min = 0.02e-4		! um -> cm
!delta_max = 0.30e-4		! um -> cm
dt_min = 2.0
delta_min = 0.05e-4		! um -> cm
delta_max = 0.50e-4		! um -> cm
delta_tmove = dt_min
ndt = 10

! Test forces
!dx = (1.1*xcross2_force - 1.05*x0_force)/100
!write(nflog,'(3f8.3)') 1.1*xcross2_force,1.05*x1_force,dx
!do i = 1,100
!	x = 1.05*x0_force + i*dx
!	if (x > xcross2_force) then
!		F = 0
!	else
!		F = a_force/((x-x0_force)*(x1_force-x)) + b_force
!	endif
!	write(nflog,*) i,x,F
!enddo
!stop	

end subroutine


end module


