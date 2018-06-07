module nbr

use global

implicit none

contains

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine getquads(cp,ix1,ix2,iy1,iy2,iz1,iz2,bcentre)
type(cell_type), pointer :: cp
integer :: ix1,ix2,iy1,iy2,iz1,iz2
real(REAL_KIND) :: bcentre(3)
integer :: nspheres
real(REAL_KIND) :: c(3)
real(REAL_KIND) :: dtol = 4*Raverage

!write(nflog,'(a,3f8.4)') 'bcentre: ',bcentre
!write(nflog,'(a,3f8.4)') 'dtol: ',dtol
!write(nflog,'(a,3f8.4)') 'bcentre+dtol: ',bcentre+dtol
!write(nflog,'(a,3f8.4)') 'bcentre-dtol: ',bcentre-dtol
!write(nflog,'(a,4f8.4)') 'c1: ',cp%centre(:,1), cp%radius(1)
nspheres = cp%nspheres
c = cp%centre(:,1)
!if (nspheres == 2) write(nflog,*) 'c2: ',cp%centre(:,2),cp%radius(2)
! for x
if (c(1) < bcentre(1) - dtol) then
	ix1 = 1
	ix2 = 1
elseif (c(1) > bcentre(1) + dtol) then
	ix1 = 2
	ix2 = 2
else
	ix1 = 1
	ix2 = 2
endif
if (nspheres == 2) then
	c = cp%centre(:,2)
	if (c(1) < bcentre(1) - dtol) then
		ix1 = min(ix1,1)
	elseif (c(1) > bcentre(1) + dtol) then
		ix2 = max(ix2,2)
	else
		ix1 = 1
		ix2 = 2
	endif
endif
!write(nflog,*) 'ix1,ix2: ',ix1,ix2
! for y
if (c(2) < bcentre(2) - dtol) then
	iy1 = 1
	iy2 = 1
elseif (c(2) > bcentre(2) + dtol) then
	iy1 = 2
	iy2 = 2
else
	iy1 = 1
	iy2 = 2
endif
if (nspheres == 2) then
	c = cp%centre(:,2)
	if (c(2) < bcentre(2) - dtol) then
		iy1 = min(iy1,1)
	elseif (c(2) > bcentre(2) + dtol) then
		iy2 = max(iy2,2)
	else
		iy1 = 1
		iy2 = 2
	endif
endif
!write(nflog,*) 'iy1,iy2: ',iy1,iy2
! for z
if (c(3) < bcentre(3) - dtol) then
	iz1 = 1
	iz2 = 1
elseif (c(3) > bcentre(3) + dtol) then
	iz1 = 2
	iz2 = 2
else
	iz1 = 1
	iz2 = 2
endif
if (nspheres == 2) then
	c = cp%centre(:,2)
	if (c(3) < bcentre(3) - dtol) then
		iz1 = min(iz1,1)
	elseif (c(3) > bcentre(3) + dtol) then
		iz2 = max(iz2,2)
	else
		iz1 = 1
		iz2 = 2
	endif
endif
!write(nflog,*) 'iz1,iz2: ',iz1,iz2
!stop
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine setup_nbrlists(ok)
logical :: ok
type(cell_type), pointer :: cp1, cp2
real(REAL_KIND) :: R1, c1(3), R2, c2(3), v(3), d2, d, Rsum, dfactor, c(3), bcentre(3)
integer :: kcell, k2, k, isphere1, nspheres1, isphere2, nspheres2, nbrs, i, j
integer :: nbrlist(1000), t(1000)
real(REAL_KIND) :: nbr_d(1000)
logical :: near, incontact, contact(2,2)
integer :: ix1, ix2, iy1, iy2, iz1, iz2, ix, iy, iz, kount(2,2,2), listmax
integer, allocatable :: clist(:,:,:,:)
logical :: dbug = .false.

listmax = max(1000,nlist/2)
allocate(clist(2,2,2,listmax))
bcentre = blobcentre0
!write(*,*) 'bcentre: ',bcentre
kount = 0
do kcell = 1,nlist
	cp1 => cell_list(kcell)
	if (cp1%state == DEAD) cycle
    call getquads(cp1,ix1,ix2,iy1,iy2,iz1,iz2,bcentre)
!	write(nflog,'(7i6)') kcell,ix1,ix2,iy1,iy2,iz1,iz2
	do ix = ix1,ix2
		do iy = iy1, iy2
			do iz = iz1, iz2
				kount(ix,iy,iz) = kount(ix,iy,iz) + 1
				if (kount(ix,iy,iz) > listmax) then
					write(*,*) 'Error: setup_nbrlists: kount > listmax: ',ix,iy,iz,kount(ix,iy,iz)
					write(nflog,*) 'Error: setup_nbrlists: kount > listmax: ',ix,iy,iz,kount(ix,iy,iz)
					stop
				endif
				clist(ix,iy,iz,kount(ix,iy,iz)) = kcell
			enddo
		enddo
	enddo
enddo
!write(*,*) 'kount: ',kount

do kcell = 1,nlist
	cp1 => cell_list(kcell)
	if (cp1%state == DEAD) cycle
    call getquads(cp1,ix1,ix2,iy1,iy2,iz1,iz2,bcentre)
    nspheres1 = cp1%nspheres
	nbrs = 0
	do ix = ix1,ix2
		do iy = iy1,iy2
			do iz = iz1,iz2
				do j = 1,kount(ix,iy,iz)
					k2 = clist(ix,iy,iz,j)
					if (k2 == kcell) cycle
					cp2 => cell_list(k2)
					nspheres2 = cp2%nspheres
					near = .false.
					do isphere1 = 1,nspheres1
!						R1 = cp1%radius(isphere1)
						c1 = cp1%centre(:,isphere1)
						do isphere2 = 1,nspheres2
!							R2 = cp2%radius(isphere2)
							c2 = cp2%centre(:,isphere2)
!							Rsum = R1 + R2
							v = c2 - c1
							d2 = dot_product(v,v)
							d = sqrt(d2)
							if (d < d_nbr_limit) then
								near = .true.
								exit
							endif
						enddo
						if (near) exit
					enddo
					if (near) then
						nbrs = nbrs + 1
						if (nbrs > 1000) then
							write(logmsg,'(a,5i6)') 'Size of nbrlist exceeded: istep, kcell, ncells: ',istep,kcell,ncells,nbrs,1000
							call logger(logmsg)
							write(nflog,'(10i6)') nbrlist(:)
							ok = .false.
							return
						endif
						nbrlist(nbrs) = k2
						nbr_d(nbrs) = d
					endif
				enddo
			enddo
		enddo
	enddo
!   Note: replacing method of creating nbrlist with one that orders the neighbours by distance, and keeps the closest.
	do i = 1,nbrs
		t(i) = i
	enddo
	if (nbrs == 0) then
!		write(*,*) 'setup_nbrlists: nbrs = 0: kcell: ',kcell
	else
		call qqsort(nbr_d,nbrs,t)     ! sort in increasing order
	endif
	! Now the ith nearest is nbrlist(t(i))
	! set cp1%nbrlist(:) as the closest #
	cp1%nbrs = min(nbrs,NBR_LIST_MAX)
	do i = 1,cp1%nbrs
		cp1%nbrlist(i)%indx = nbrlist(t(i))
	enddo
!	write(nflog,*) 'cell, nbrs: ',kcell,cp1%nbrs
enddo
deallocate(clist)
!stop
end subroutine

!-----------------------------------------------------------------------------------------
! Dumb version to start
! This takes a lot of time when Ncells gets large: about 40% of time when Ncells = 25k
!-----------------------------------------------------------------------------------------
subroutine setup_nbrlists1(ok)
logical :: ok
type(cell_type), pointer :: cp1, cp2
real(REAL_KIND) :: R1, c1(3), R2, c2(3), v(3), d2, d, Rsum, dfactor
integer :: kcell, k2, k, isphere1, nspheres1, isphere2, nspheres2, nbrs, i
!type(neighbour_type) :: nbrlist(1000)
integer :: nbrlist(1000), t(1000)
real(REAL_KIND) :: nbr_d(1000)
logical :: near, incontact, contact(2,2)
logical :: dbug = .false.

!call logger('start setup_nbrlists')
do kcell = 1,nlist
	cp1 => cell_list(kcell)
	if (cp1%state == DEAD) cycle
	!if (cp1%Iphase) then
	!	nspheres1 = 1
	!else
	!	nspheres1 = 2
	!endif
    nspheres1 = cp1%nspheres
	nbrs = 0
	do k2 = 1,nlist
		if (k2 == kcell) cycle
		cp2 => cell_list(k2)
		if (cp2%state == DEAD) cycle
		!if (cp2%Iphase) then
		!	nspheres2 = 1
		!else
		!	nspheres2 = 2
		!endif
        nspheres2 = cp2%nspheres
		near = .false.
		incontact = .false.
		contact = .false.
		do isphere1 = 1,nspheres1
			R1 = cp1%radius(isphere1)
			c1 = cp1%centre(:,isphere1)
			do isphere2 = 1,nspheres2
				R2 = cp2%radius(isphere2)
				c2 = cp2%centre(:,isphere2)
				Rsum = R1 + R2
				v = c2 - c1
				d2 = dot_product(v,v)
				d = sqrt(d2)
				if (use_hysteresis) then
					if (d < d_nbr_limit) then
						!if (dbug) write(*,'(5f8.2)') R1,R2,Rsum,d,d_nbr_limit
						dfactor = 1
						do k = 1,cp1%nbrs
							if (cp1%nbrlist(k)%indx == k2) then
								if (cp1%nbrlist(k)%contact(isphere1,isphere2)) dfactor = k_detach
								exit
							endif
						enddo
						near = .true.
						if (d < dfactor*Rsum) then
							incontact = .true.
							contact(isphere1,isphere2) = .true.
						endif
					endif
				else
					if (d < d_nbr_limit) then
						near = .true.
						exit
					endif
				endif
			enddo
			if (near) exit
		enddo
		if (near) then
			nbrs = nbrs + 1
!			if (nbrs > MAX_NBRS) then
			if (nbrs > 1000) then
				write(logmsg,'(a,5i6)') 'Size of nbrlist exceeded: istep, kcell, ncells: ',istep,kcell,ncells,nbrs,1000
				call logger(logmsg)
!				write(nflog,'(10i6)') nbrlist(:)%indx
				write(nflog,'(10i6)') nbrlist(:)
				ok = .false.
				return
			endif
!			nbrlist(nbrs)%indx = k2
!			nbrlist(nbrs)%contact = contact
			nbrlist(nbrs) = k2
			nbr_d(nbrs) = d
		endif
	enddo
!	cp1%nbrs = nbrs
!	cp1%nbrlist(1:nbrs) = nbrlist(1:nbrs)
!   Note: replacing method of creating nbrlist with one that orders the neighbours by distance, and keeps the closest.
	do i = 1,nbrs
		t(i) = i
	enddo
	if (nbrs == 0) then
!		write(*,*) 'setup_nbrlists: nbrs = 0: kcell: ',kcell
	else
		call qqsort(nbr_d,nbrs,t)     ! sort in increasing order
	endif
	! Now the ith nearest is nbrlist(t(i))
	! set cp1%nbrlist(:) as the closest #
	cp1%nbrs = min(nbrs,NBR_LIST_MAX)
	do i = 1,cp1%nbrs
		cp1%nbrlist(i)%indx = nbrlist(t(i))
	enddo

enddo
!call logger('end setup_nbrlists')
ok = .true.
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine update_all_nbrlists
integer :: kcell, max_nbrs

!write(nflog,*) 'update_all_nbrlists'
max_nbrs = 0
do kcell = 1,nlist
	if (cell_list(kcell)%state == DEAD) cycle
	call update_nbrlist(kcell)
	max_nbrs = max(max_nbrs,cell_list(kcell)%nbrs)
enddo
!write(nflog,*) 'max_nbrs: ',max_nbrs
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine update_nbrlist(kcell)
integer :: kcell
type (cell_type), pointer :: cp1, cp2, cp3
integer :: nbrlist(1000), t(1000)
real(REAL_KIND) :: nbr_d(1000)
integer :: nbrs, ns1, k, k2, k3, knbr, i
real(REAL_KIND) :: c1(3,2), d
logical :: inlist

cp1 => cell_list(kcell)
if (cp1%Iphase) then
    ns1 = 1
    c1(:,1) = cp1%centre(:,1)
else
    ns1 = 2
    c1(:,1) = cp1%centre(:,1)
    c1(:,2) = cp1%centre(:,2)
endif
nbrs = 0
do k2 = 1,cp1%nbrs
    knbr = cp1%nbrlist(k2)%indx
    if (knbr == kcell) cycle
    cp2 => cell_list(knbr)
	if (cp2%state == DEAD) cycle
    call get_dist(c1,ns1,cp2,d)
    if (d < d_nbr_limit) then
        call add_to_list(knbr,d,nbrlist,nbr_d,nbrs)
    endif
    do k3 = 1,cp2%nbrs
        knbr = cp2%nbrlist(k3)%indx
        if (knbr == kcell) cycle
        cp3 => cell_list(knbr)
		if (cp3%state == DEAD) cycle
        call get_dist(c1,ns1,cp3,d)
        if (d < d_nbr_limit) then
            call add_to_list(knbr,d,nbrlist,nbr_d,nbrs)
        endif
    enddo
enddo
if (nbrs == 0) return
! sort nbrlist(:) by nbr_d
do i = 1,nbrs
    t(i) = i
enddo
call qqsort(nbr_d,nbrs,t)     ! sort in increasing order
! Now the ith nearest is nbrlist(t(i))
! set cp1%nbrlist(:) as the closest #
cp1%nbrs = min(nbrs,NBR_LIST_MAX)
do i = 1,cp1%nbrs
	cp1%nbrlist(i)%indx = nbrlist(t(i))
enddo
! Need to ensure that if knbr is in nbrlist then kcell is in the list for knbr
do i = 1,cp1%nbrs
	knbr = cp1%nbrlist(i)%indx
    cp2 => cell_list(knbr)
	if (cp2%state == DEAD) cycle
	inlist = .false.
    do k2 = 1,cp2%nbrs
		if (kcell == cp2%nbrlist(k2)%indx) then
			inlist = .true.
			exit
		endif
	enddo
	if (.not.inlist) then
		cp2%nbrs = cp2%nbrs + 1
		cp2%nbrlist(cp2%nbrs)%indx = kcell
	endif
enddo
	
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine get_dist(c1,ns1,cp2,d)
integer :: ns1
type (cell_type), pointer :: cp2
real(REAL_KIND) :: c1(3,2), d
real(REAL_KIND) :: r2min, r(3)
integer :: is1, is2, ns2

if (cp2%Iphase) then
    ns2 = 1
else
    ns2 = 2
endif
r2min = 1.0e10
do is1 = 1,ns1
    do is2 = 1,ns2
        r = c1(:,is1) - cp2%centre(:,is2)
        r2min = min(r2min,dot_product(r,r))
    enddo
enddo
d = sqrt(r2min)
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine add_to_list(knbr,d,nbrlist,nbr_d,nbrs)
integer :: knbr, nbrs, nbrlist(:)
real(REAL_KIND) :: d,nbr_d(:)
integer :: k

do k = 1,nbrs
    if (nbrlist(k) == knbr) return
enddo
nbrs = nbrs + 1
nbrlist(nbrs) = knbr
nbr_d(nbrs) = d
end subroutine

!--------------------------------------------------------------------------------
!     NON-RECURSIVE STACK VERSION OF QUICKSORT FROM N.WIRTH'S PASCAL
!     BOOK, 'ALGORITHMS + DATA STRUCTURES = PROGRAMS'.
!     SINGLE PRECISION, ALSO CHANGES THE ORDER OF THE ASSOCIATED ARRAY T.
!--------------------------------------------------------------------------------
SUBROUTINE qqsort(a, n, t)
IMPLICIT NONE

INTEGER, INTENT(IN)    :: n
REAL(REAL_KIND), INTENT(INOUT)    :: a(n)
INTEGER, INTENT(INOUT) :: t(n)

!     Local Variables

INTEGER         :: i, j, k, l, r, s, stackl(15), stackr(15), ww
REAL(REAL_KIND) :: w, x

s = 1
stackl(1) = 1
stackr(1) = n

!     KEEP TAKING THE TOP REQUEST FROM THE STACK UNTIL S = 0.

10 CONTINUE
l = stackl(s)
r = stackr(s)
s = s - 1

!     KEEP SPLITTING A(L), ... , A(R) UNTIL L >= R.

20 CONTINUE
i = l
j = r
k = (l+r) / 2
x = a(k)

!     REPEAT UNTIL I > J.

DO
  DO
    IF (a(i).LT.x) THEN                ! Search from lower end
      i = i + 1
      CYCLE
    ELSE
      EXIT
    END IF
  END DO

  DO
    IF (x.LT.a(j)) THEN                ! Search from upper end
      j = j - 1
      CYCLE
    ELSE
      EXIT
    END IF
  END DO

  IF (i.LE.j) THEN                     ! Swap positions i & j
    w = a(i)
    ww = t(i)
    a(i) = a(j)
    t(i) = t(j)
    a(j) = w
    t(j) = ww
    i = i + 1
    j = j - 1
    IF (i.GT.j) EXIT
  ELSE
    EXIT
  END IF
END DO

IF (j-l.GE.r-i) THEN
  IF (l.LT.j) THEN
    s = s + 1
    stackl(s) = l
    stackr(s) = j
  END IF
  l = i
ELSE
  IF (i.LT.r) THEN
    s = s + 1
    stackl(s) = i
    stackr(s) = r
  END IF
  r = j
END IF

IF (l.LT.r) GO TO 20
IF (s.NE.0) GO TO 10

RETURN
END SUBROUTINE qqsort

end module
