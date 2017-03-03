! To test solution of reaction-diffusion on a rectangular grid, using fgmres with ITSOL ILUK preconditioning
! This version uses two rectangular grids, coarse and fine.
! The fine grid is large enough to contain the spheroid blob completely, and has grid resolution
! such that a grid-cell can hold a few tumour cells, e.g. dxf = 20 -> 30 um
! The fine grid is embedded in a coarse grid, which has a resolution dxb that is a multiple of dxf,
! e.g. dxb = 4*dxf. Coarse grid points coincide with points in the fine grid, and in particular
! the points on the boundary of the fine grid coincides with a line of the course grid.
! The size of the coarse grid is set to match the medium volume, if not the shape.
! 
! The idea behind the solution method is that solving on the coarse grid provides the boundary
! values for solution on the fine grid.
!
! Quasi-steady-state method:
! -------------------------
! Time-dependence is handled in the coarse solution, using the IMEX 2-SBDF scheme.  The most recent 
! rates of uptake-secretion of constituents by cells are used to compute grid flux F().
! After the Cave values on the coarse grid have been determined, the boundary concentrations on
! the fine grid are estimated by interpolation, and the steady-state solution on the fine grid
! is computed.  It is this solution that provides the flux values for the next coarse grid time step.
! 
! Units:
!     time				s = seconds
!     distance			cm
!     volume			cm^3
!     mass				micromole = 10^-6 mol = mumol
!     flux				mumol/s
!     concentration		mumol/cm^3 = mM
!
! Parallel issues:
! https://software.intel.com/en-us/articles/threading-fortran-applications-for-parallel-performance-on-multi-core-systems/
!
module react_diff

use real_kind_mod
use global
use omp_lib
use sparse_map
use par_zig_mod
use chemokine
use continuum
use metabolism
use ode_solver

use, intrinsic :: iso_c_binding

implicit none

contains

!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
subroutine setup_react_diff(ok)
logical :: ok
integer :: ix, iy, iz, ic, maxnz, ichemo
real(REAL_KIND) :: C0
character*(10) :: emapfile, bmapfile
real(REAL_KIND), pointer :: Cprev(:,:,:), Fprev(:,:,:), Fcurr(:,:,:)
real(REAL_KIND), pointer :: Cave_b(:,:,:), Cprev_b(:,:,:), Fprev_b(:,:,:), Fcurr_b(:,:,:)
logical :: zero

write(nflog,*) 'setup_react_diff: NX,NY,NZ: ',NX,NY,NZ
dxf = DELTA_X
dx3 = dxf*dxf*dxf
nrow = (NX-2)*(NY-2)*(NZ-1)		! embedded map, Dirichlet conditions on ix=1,NX, iy=1,NY, iz=NZ
maxnz = MAX_CHEMO*nrow
if (allocated(amap)) deallocate(amap)
if (allocated(ja)) deallocate(ja)
if (allocated(ia)) deallocate(ia)
allocate(amap(maxnz,0:3))
allocate(ja(maxnz))
allocate(ia(nrow+1))

dxb3 = dxb*dxb*dxb
nrow_b = NXB*NYB*NZB
maxnz = MAX_CHEMO*nrow_b
if (allocated(amap_b)) deallocate(amap_b)
if (allocated(ja_b)) deallocate(ja_b)
if (allocated(ia_b)) deallocate(ia_b)
allocate(amap_b(maxnz,0:3))
allocate(ja_b(maxnz))
allocate(ia_b(nrow_b+1))

emapfile = ''
write(emapfile,'(a,i2.0,a)') 'emap',NX,'.dat'
write(nflog,*) 'setup_react_diff: ',emapfile
call make_sparse_emap(emapfile,.true.,ok)
if (.not.ok) return
write(nflog,*) 'made emapfile: ',emapfile

bmapfile = ''
write(bmapfile,'(a,i2.0,a)') 'bmap',NXB,'.dat'
write(nflog,*) 'setup_react_diff: ',bmapfile
call make_sparse_map(bmapfile,.false.,ok)
if (.not.ok) return
write(nflog,*) 'made bmapfile: ',bmapfile

call make_grid_flux_weights(ok)
do ichemo = 1,MAX_CHEMO
	if (.not.chemo(ichemo)%used) cycle
	if (ichemo > TRACER) then
		chemo(ichemo)%bdry_conc = 0
	endif
	Cextra_all(:,:,:,ichemo) = chemo(ichemo)%bdry_conc
	Caverage(:,:,:,ichemo) = chemo(ichemo)%bdry_conc
	chemo(ichemo)%medium_Cbnd = chemo(ichemo)%bdry_conc
	write(nflog,*) 'ichemo, Caverage: ',ichemo,chemo(ichemo)%bdry_conc
!	write(*,*) 'call update_Cex_Cin: ',ichemo
	call update_Cex_Cin_const(ichemo)	! initialise with SS values
	Fcurr => Cflux(:,:,:,ichemo)
	call getF_const(ichemo,Fcurr)	! estimates all fine grid pt fluxes Cflux(:,:,:,:) from Cextra(:,:,:,:)
! Sum fine grid fluxes to initialise Fcurr_b, Fprev_b
	C0 = chemo(ichemo)%bdry_conc
	Cprev => chemo(ichemo)%Cprev
	Fprev => chemo(ichemo)%Fprev
	Cprev_b => chemo(ichemo)%Cprev_b
	Fprev_b => chemo(ichemo)%Fprev_b
	Fcurr_b => chemo(ichemo)%Fcurr_b
	Cave_b => chemo(ichemo)%Cave_b
	
	Cave_b = C0
	Cprev_b = C0
	Cprev = C0
	Fprev = Cflux(:,:,:,ichemo)
	call makeF_b(ichemo,Fprev_b,Fprev,DELTA_T,zero)
	Fcurr_b = Fprev_b
enddo
ok = .true.
end subroutine

!-------------------------------------------------------------------------------------------
! The flux values on the coarse grid are derived from the values on the fine grid by
! appropriate summation.  It's important that the total flux is consistent.
!-------------------------------------------------------------------------------------------
subroutine makeF_b(ichemo,F_b,F,dt,zero)
integer, parameter :: dN = NRF/2
real(REAL_KIND) :: F_b(:,:,:), F(:,:,:), dt
logical :: zero
real(REAL_KIND) :: Fsum_b, Fsum, Fmin, dF
real(REAL_KIND) :: wt(-dN:dN,-dN:dN,-dN:dN)
real(REAL_KIND) :: xfac, yfac, zfac
integer :: idx, idy, idz, xb0, yb0, idxb, idyb, zb1
integer :: ixb, iyb, izb, ix0, iy0, iz0, ix, iy, iz
integer :: ichemo

xb0 = (NXB+1)/2			! these must all be integer - check NX,NY,NZ,NXB,NYB,NZB,NRF earlier
idxb = (NX-1)/(2*NRF)
yb0 = (NYB+1)/2
idyb = (NY-1)/(2*NRF)
zb1 = (NZ-1)/NRF + 1

do idx = -dN,dN
do idy = -dN,dN
do idz = -dN,dN
	if (idx > -dN .and. idx < dN) then
		xfac = 1
	else
		xfac = 0.5
	endif
	if (idy > -dN .and. idy < dN) then
		yfac = 1
	else
		yfac = 0.5
	endif
	if (idz > -dN .and. idz < dN) then
		zfac = 1
	else
		zfac = 0.5
	endif
	wt(idx,idy,idz) = xfac*yfac*zfac
enddo
enddo
enddo

Fsum_b = 0
F_b = 0
do ixb = xb0-idxb,xb0+idxb
	ix0 = (ixb - xb0)*NRF + (NX+1)/2
	do iyb = yb0-idyb,yb0+idyb
		iy0 = (iyb - yb0)*NRF + (NY+1)/2
		do izb = 1,zb1
			iz0 = (izb - 1)*NRF + 1
			do idx = -dN,dN
				ix = ix0 + idx
				if (ix < 1 .or. ix > NX) cycle
				do idy = -dN,dN
					iy = iy0 + idy
					if (iy < 1 .or. iy > NY) cycle
					do idz = -dN,dN
						iz = iz0 + idz
						if (iz < 1 .or. iz > NZ) cycle
						dF = wt(idx,idy,idz)*F(ix,iy,iz)
						F_b(ixb,iyb,izb) = F_b(ixb,iyb,izb) + dF
						Fsum_b = Fsum_b + dF
					enddo
				enddo
			enddo
		enddo
	enddo
enddo
zero = (Fsum_b == 0)
if (dbug .and. ichemo >= DRUG_A) write(*,'(a,i2,e12.3)') 'flux makeF_b: DRUG_A sum: ',ichemo,Fsum_b
end subroutine

!-------------------------------------------------------------------------------------------
! This version is for the embedded grid.
! Given boundary concentrations and grid point fluxes, we need to solve for the 
! steady-state concentration field.
!-------------------------------------------------------------------------------------------
subroutine make_csr_SS(a, ichemo, Cave, Fcurr, rhs)
integer :: ichemo
real(REAL_KIND) :: a(:), Cave(:,:,:), Fcurr(:,:,:), rhs(:)
integer :: k, ix, iy, iz, krow, kcol, nc, idx, idy, idz, ixx, iyy, izz, n, ncsum
integer :: nc_max = 10	! just a wild guess but not a bad one
real(REAL_KIND) :: Ktissue, Kmedium, Kdiff, alfa, Kr, Vex, Cbdry, Kdiff_sum
logical, save :: first = .true.

Ktissue = chemo(ichemo)%diff_coef
Kmedium = chemo(ichemo)%medium_diff_coef

krow = 0
do k = 1,nnz
	if (k == ia(krow+1)) krow = krow+1
	kcol = ja(k)
	a(k) = amap(k,0)
enddo

Kdiff_sum = 0
k = 0
do ix = 2,NX-1
	do iy = 2,NY-1
		do iz = 1,NZ-1
			krow = (ix-2)*(NY-2)*(NZ-1) + (iy-2)*(NZ-1) + iz
			
			! This is very crude!!!!!!!!!!!!!!! 
			! nc = # of cells in the gridcell with LL corner = (ix,iy,iz), not centred on the grid pt
			! nc = grid(ix,iy,iz)%nc
			! Try making it symmetrical about the grid pt
			ncsum = 0
			n = 0
			do idx = -1,0
				ixx = ix + idx
				do idy = -1,0
					iyy = iy + idy
					do idz = -1,0
						izz = iz + idz
						if (izz == 0) cycle
						ncsum = ncsum + grid(ixx,iyy,izz)%nc
						n = n+1
					enddo
				enddo
			enddo
			nc = ncsum/n
			
			alfa = min(nc,nc_max)/nc_max
!			Kdiff = Kdiff*(1 - chemo(ichemo)%diff_reduction_factor*alfa)
			! Kdiff should range between Kmedium and Ktissue as nc goes from 0 to nc_max
			Kdiff = (1-alfa)*Kmedium + alfa*Ktissue
!			Kdiff = Kmedium
!			if (nc > 0) then
!				write(*,'(a,2i4,f8.4)') 'Kdiff reduction: ',nc,nc_max,chemo(ichemo)%diff_reduction_factor*min(nc,nc_max)/nc_max
!				k = k+1
!				Kdiff_sum = Kdiff_sum + Kdiff
!			endif
			Kr = 1/(dxf*Kdiff)
			rhs(krow) = -Kr*Fcurr(ix,iy,iz)
		enddo
	enddo
enddo
if (k > 0) write(nflog,'(a,i4,e12.3)') 'Kdiff ave: ',ichemo,Kdiff_sum/k
ix = 2
do iy = 2,NY-1
	do iz = 1,NZ-1
		krow = (ix-2)*(NY-2)*(NZ-1) + (iy-2)*(NZ-1) + iz
		Cbdry = Cave(1,iy,iz)
		rhs(krow) = rhs(krow) + Cbdry
	enddo
enddo
ix = NX-1
do iy = 2,NY-1
	do iz = 1,NZ-1
		krow = (ix-2)*(NY-2)*(NZ-1) + (iy-2)*(NZ-1) + iz
		Cbdry = Cave(NX,iy,iz)
		rhs(krow) = rhs(krow) + Cbdry		
	enddo
enddo
iy = 2
do ix = 2,NX-1
	do iz = 1,NZ-1
		krow = (ix-2)*(NY-2)*(NZ-1) + (iy-2)*(NZ-1) + iz
		Cbdry = Cave(ix,1,iz)
		rhs(krow) = rhs(krow) + Cbdry
	enddo
enddo
iy = NY-1
do ix = 2,NX-1
	do iz = 1,NZ-1
		krow = (ix-2)*(NY-2)*(NZ-1) + (iy-2)*(NZ-1) + iz
		Cbdry = Cave(ix,NY,iz)
		rhs(krow) = rhs(krow) + Cbdry
	enddo
enddo
iz = NZ-1
do ix = 2,NX-1
	do iy = 2,NY-1
		krow = (ix-2)*(NY-2)*(NZ-1) + (iy-2)*(NZ-1) + iz
		Cbdry = Cave(ix,iy,NZ)
		rhs(krow) = rhs(krow) + Cbdry
	enddo
enddo
first = .false.
end subroutine

!-------------------------------------------------------------------------------------------
! This is the version for the coarse grid.
! Use variable numbering (ix,iy) -> k = (ix-1)*NY + iy
! Since the equation is now dM/dt = ... need to derive different expression for Kr
! Need to add treatment of top boundary for O2
! Need to add decay!  For now add it explicitly in solver
!-------------------------------------------------------------------------------------------
subroutine make_csr_b(a_b, ichemo, dt, Cave_b, Cprev_b, Fcurr_b, Fprev_b, rhs, zero)
integer :: ichemo
real(REAL_KIND) :: dt, a_b(:), Cave_b(:,:,:), Cprev_b(:,:,:), Fcurr_b(:,:,:), Fprev_b(:,:,:), rhs(:)
logical :: zero
integer :: ixb, iyb, izb, k, i, krow, kcol, nc
real(REAL_KIND) :: Kdiff, Ktissue, Kmedium, Kr, Cbdry, Fsum
integer, parameter :: m = 3

zero = .true.
!Kdiff = chemo(ichemo)%medium_diff_coef
Ktissue = chemo(ichemo)%diff_coef
Kmedium = chemo(ichemo)%medium_diff_coef
Fsum = 0
krow = 0
Kdiff = Kmedium
do k = 1,nnz_b
	if (k == ia_b(krow+1)) krow = krow+1
	kcol = ja_b(k)
	if (amap_b(k,0) == 2*m) then
    	Kr = dxb*dxb/Kdiff
		a_b(k) = 3*Kr/(2*dt) + 2*m
	else
		a_b(k) = amap_b(k,0)
	endif
	if (ichemo == OXYGEN) then
		if (amap_b(k,3) == NZB .and. kcol == krow-1) then
			if (a_b(k) /= -2) then
				write(nflog,*) 'Error in OXYGEN bdry adjustment'
				stop
			endif
			a_b(k) = -1
		endif
	endif
enddo

! To make OMP useable, need to ensure that values of krow concurrently evaluated are well separated.
! Using ixb as the loop variable should suffice.
!!$omp parallel do private(iyb, izb, krow, Kr)
do ixb = 1,NXB
	do iyb = 1,NYB
		do izb = 1,NZB
			krow = (ixb-1)*NYB*NZB + (iyb-1)*NZB + izb
!			if (Fcurr_b(ixb,iyb,izb) > 0) then
!			    Kdiff = Ktissue
!			else
!			    Kdiff = Kmedium
!			endif
!		    Fsum = Fsum + Fcurr_b(ixb,iyb,izb)
			Kr = dxb*dxb/Kdiff
			rhs(krow) = Kr*((-2*Fcurr_b(ixb,iyb,izb) + Fprev_b(ixb,iyb,izb))/dxb3 + (1./(2*dt))*(4*Cave_b(ixb,iyb,izb) - Cprev_b(ixb,iyb,izb)))
			if (rhs(krow) /= 0) zero = .false.
		enddo
	enddo
enddo
!!$omp end parallel do

if (ichemo == OXYGEN) then
	Cbdry = chemo(ichemo)%bdry_conc
	izb = NZB
	do ixb = 1,NXB
		do iyb = 1,NYB
			krow = (ixb-1)*NYB*NZB + (iyb-1)*NZB + izb
			rhs(krow) = rhs(krow) + Cbdry		
		enddo
	enddo
endif
!write(nflog,*) 'make_csr_b: Fsum: ',ichemo,Fsum
end subroutine

!-------------------------------------------------------------------------------------------
! Update Cextra(:,:,:) and thereby cp%Cex(ichemo) and cp%Cin(ichemo) from Caverage(:,:,:,ichemo), 
!-------------------------------------------------------------------------------------------
subroutine update_Cex_Cin_const(ichemo)
integer :: ichemo
integer :: kcell
real(REAL_KIND) :: Kin, Kout, dC, dCdt, dMdt
type(cell_type), pointer :: cp
real(REAL_KIND), pointer :: Cextra(:,:,:)

!write(*,*) 'update_Cex_Cin_const: ichemo: ',ichemo

Cextra => Caverage(:,:,:,ichemo)		! currently using the average concentration!

Kin = chemo(ichemo)%membrane_diff_in
Kout = chemo(ichemo)%membrane_diff_out
do kcell = 1,nlist
	cp => cell_list(kcell)
	if (cp%state == DEAD .or. cp%state == DYING) cycle
	call extra_concs_const(kcell, Cextra, cp%Cex(ichemo))
	cp%Cin(ichemo) = getCin_SS(kcell,ichemo,cp%V,cp%Cex(ichemo))
enddo
end subroutine

!-------------------------------------------------------------------------------------------
! Update Cextra(:,:,:) and thereby cp%Cex(ichemo) and cp%Cin(ichemo) from Caverage(:,:,:,ichemo), 
!-------------------------------------------------------------------------------------------
subroutine update_Cex_Cin_dCdt_const(ichemo, dt)
integer :: ichemo
real(REAL_KIND) :: dt
integer :: kcell, ix, iy, iz
real(REAL_KIND) :: alfa(3), Clast, Kin, Kout, dC, dCexdt, dMdt
type(cell_type), pointer :: cp
real(REAL_KIND), pointer :: Cextra(:,:,:)
real(REAL_KIND), pointer :: Cprev(:,:,:)

!write(*,*) 'update_Cex_Cin_dCdt_const: ',ichemo

Cextra => Caverage(:,:,:,ichemo)		! currently using the average concentration!
Cprev => chemo(ichemo)%Cprev
Kin = chemo(ichemo)%membrane_diff_in
Kout = chemo(ichemo)%membrane_diff_out
!$omp parallel do private(cp, ix, iy, iz, alfa, Clast, dCexdt)
do kcell = 1,nlist
	cp => cell_list(kcell)
	if (cp%state == DEAD .or. cp%state == DYING) cycle
	call grid_interp(kcell, alfa)
	ix = cp%site(1)
	iy = cp%site(2)
	iz = cp%site(3)
	cp%Cex(ichemo) = (1-alfa(1))*(1-alfa(2))*(1-alfa(3))*Cextra(ix,iy,iz)  &
        + (1-alfa(1))*alfa(2)*(1-alfa(3))*Cextra(ix,iy+1,iz)  &
        + (1-alfa(1))*alfa(2)*alfa(3)*Cextra(ix,iy+1,iz+1)  &
        + (1-alfa(1))*(1-alfa(2))*alfa(3)*Cextra(ix,iy,iz+1)  &
        + alfa(1)*(1-alfa(2))*(1-alfa(3))*Cextra(ix+1,iy,iz)  &
        + alfa(1)*alfa(2)*(1-alfa(3))*Cextra(ix+1,iy+1,iz)  &
        + alfa(1)*alfa(2)*alfa(3)*Cextra(ix+1,iy+1,iz+1)  &
        + alfa(1)*(1-alfa(2))*alfa(3)*Cextra(ix+1,iy,iz+1)
	Clast = (1-alfa(1))*(1-alfa(2))*(1-alfa(3))*Cprev(ix,iy,iz)  &
        + (1-alfa(1))*alfa(2)*(1-alfa(3))*Cprev(ix,iy+1,iz)  &
        + (1-alfa(1))*alfa(2)*alfa(3)*Cprev(ix,iy+1,iz+1)  &
        + (1-alfa(1))*(1-alfa(2))*alfa(3)*Cprev(ix,iy,iz+1)  &
        + alfa(1)*(1-alfa(2))*(1-alfa(3))*Cprev(ix+1,iy,iz)  &
        + alfa(1)*alfa(2)*(1-alfa(3))*Cprev(ix+1,iy+1,iz)  &
        + alfa(1)*alfa(2)*alfa(3)*Cprev(ix+1,iy+1,iz+1)  &
        + alfa(1)*(1-alfa(2))*alfa(3)*Cprev(ix+1,iy,iz+1)
    dCexdt = (cp%Cex(ichemo) - Clast)/dt
	cp%Cin(ichemo) = getCin(kcell,ichemo,cp%V,cp%Cex(ichemo),dCexdt)
enddo
!$omp end parallel do
end subroutine

!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
subroutine update_Cin_const_SS(ichemo)
integer :: ichemo
integer :: kcell, ix, iy, iz
real(REAL_KIND) :: alfa(3), cmax
type(cell_type), pointer :: cp
real(REAL_KIND), pointer :: Cextra(:,:,:)

!write(*,*) 'update_Cin_const_SS: ',ichemo

!write(*,*) 'update_Cex_Cin_dCdt_const: ',ichemo

Cextra => Caverage(:,:,:,ichemo)		! currently using the average concentration!
cmax = 0
!$omp parallel do private(cp, ix, iy, iz, alfa)
do kcell = 1,nlist
	cp => cell_list(kcell)
	if (cp%state == DEAD .or. cp%state == DYING) cycle
	call grid_interp(kcell, alfa)
	ix = cp%site(1)
	iy = cp%site(2)
	iz = cp%site(3)
	cp%Cex(ichemo) = (1-alfa(1))*(1-alfa(2))*(1-alfa(3))*Cextra(ix,iy,iz)  &
        + (1-alfa(1))*alfa(2)*(1-alfa(3))*Cextra(ix,iy+1,iz)  &
        + (1-alfa(1))*alfa(2)*alfa(3)*Cextra(ix,iy+1,iz+1)  &
        + (1-alfa(1))*(1-alfa(2))*alfa(3)*Cextra(ix,iy,iz+1)  &
        + alfa(1)*(1-alfa(2))*(1-alfa(3))*Cextra(ix+1,iy,iz)  &
        + alfa(1)*alfa(2)*(1-alfa(3))*Cextra(ix+1,iy+1,iz)  &
        + alfa(1)*alfa(2)*alfa(3)*Cextra(ix+1,iy+1,iz+1)  &
        + alfa(1)*(1-alfa(2))*alfa(3)*Cextra(ix+1,iy,iz+1)
	cp%Cin(ichemo) = getCin_SS(kcell,ichemo,cp%V,cp%Cex(ichemo))
!	if (ichemo == OXYGEN) then
!	    cmax = max(cmax,cp%Cex(ichemo))
!	endif
!    if (ichemo == OXYGEN .and. kcell == 4593) then
!        write(nflog,'(a,3i4)') 'update_Cin_const_SS: O2: kcell=4593: ix,iy,iz: ',ix,iy,iz
!        write(nflog,'(a,e12.3,2x,3e12.3)') 'Cex, alfa: ',cp%Cex(ichemo),alfa
!    endif
enddo
!$omp end parallel do
!if (ichemo == OXYGEN) then
!    write(nflog,*) 'Cextra: ix,iy: ',17,17
!    write(nflog,'(10f8.4)') Cextra(17,17,:)
!    write(nflog,'(a,e12.3)') 'cp%Cex max (O2): ',cmax
!endif
end subroutine

!-------------------------------------------------------------------------------------------
! Updates Cex for all constituents, Cin for all except drugs.
!-------------------------------------------------------------------------------------------
subroutine update_Cex(ichemo)
integer :: ichemo
real(REAL_KIND) :: dt
integer :: kcell, ix, iy, iz
real(REAL_KIND) :: alfa(3)
type(cell_type), pointer :: cp
type(metabolism_type), pointer :: mp
real(REAL_KIND), pointer :: Cextra(:,:,:)

if (.not.chemo(ichemo)%present) return
Cextra => Caverage(:,:,:,ichemo)		! currently using the average concentration!
!if (ichemo == DRUG_A .and. Cextra(NX/2,NY/2,NZ/2) == 0) then
!	write(*,*) 'update_Cex: DRUG_A: Cextra = 0'
!	stop
!endif

!if (ichemo == DRUG_A+1) then
!	write(*,*) 'update_Cex: ichemo: ', ichemo
!	do ix = NX/2-1,NX/2+1
!		write(*,'(10f7.4)') Cextra(ix,NY/2,:)
!	enddo
!endif
!!$omp parallel do private(cp, ix, iy, iz, alfa)
do kcell = 1,nlist
	cp => cell_list(kcell)
	if (cp%state == DEAD .or. cp%state == DYING) cycle
	call grid_interp(kcell, alfa)
	ix = cp%site(1)
	iy = cp%site(2)
	iz = cp%site(3)
	cp%Cex(ichemo) = (1-alfa(1))*(1-alfa(2))*(1-alfa(3))*Cextra(ix,iy,iz)  &
        + (1-alfa(1))*alfa(2)*(1-alfa(3))*Cextra(ix,iy+1,iz)  &
        + (1-alfa(1))*alfa(2)*alfa(3)*Cextra(ix,iy+1,iz+1)  &
        + (1-alfa(1))*(1-alfa(2))*alfa(3)*Cextra(ix,iy,iz+1)  &
        + alfa(1)*(1-alfa(2))*(1-alfa(3))*Cextra(ix+1,iy,iz)  &
        + alfa(1)*alfa(2)*(1-alfa(3))*Cextra(ix+1,iy+1,iz)  &
        + alfa(1)*alfa(2)*alfa(3)*Cextra(ix+1,iy+1,iz+1)  &
        + alfa(1)*(1-alfa(2))*alfa(3)*Cextra(ix+1,iy,iz+1)
    if (ichemo < DRUG_A) then
		cp%Cin(ichemo) = getCin_SS(kcell,ichemo,cp%V,cp%Cex(ichemo))
	endif
!	if (ichemo == DRUG_A+1) write(*,'(i6,e12.3)') kcell,cp%Cex(ichemo)
enddo
!!$omp end parallel do
end subroutine

!-------------------------------------------------------------------------------------------
! Updates Cin and dMdt for drugs and metabolites.
!-------------------------------------------------------------------------------------------
subroutine integrate_Cin(dt)
real(REAL_KIND) :: dt
integer :: kcell
type(cell_type), pointer :: cp

!$omp parallel do
do kcell = 1,nlist
	cp => cell_list(kcell)
	if (cp%state == DEAD .or. cp%state == DYING) cycle
	if (chemo(DRUG_A)%present) then
		call integrate_cell_Cin(DRUG_A,kcell,dt)
	endif
	if (chemo(DRUG_B)%present) then
		call integrate_cell_Cin(DRUG_B,kcell,dt)
	endif
enddo
!$omp end parallel do
end subroutine

!-------------------------------------------------------------------------------------------
! The system quickly equilibrates.
! Would it be best to find the steady-state solution analytically?
! 3 simultaneous quadratics... but this approach is wrong, I think.
! Updates Cin and dMdt by integrating ODEs for drug and metabolites.
!-------------------------------------------------------------------------------------------
subroutine integrate_cell_Cin(ichemo_parent, kcell, dt)
integer :: ichemo_parent, kcell
real(REAL_KIND) :: dt
type(cell_type), pointer :: cp
integer :: ichemo, ictyp, idrug, im, it
real(REAL_KIND) :: Vin, Cex(0:2), Cin(0:2), dCdt(0:2), dMdt(0:2), dtt, alfa
real(REAL_KIND) :: CO2
real(REAL_KIND) :: Kin(0:2), Kout(0:2), decay_rate(0:2)
real(REAL_KIND) :: CC2, KO2(0:2), Kmet0(0:2), Vmax(0:2), Km(0:2)
real(REAL_KIND) :: K1(0:2), K2(0:2)
real(REAL_KIND) :: a0,a1,a2,b0,b1,b2,c0,c1,c2,d0,d1,d2,a,b,c,d,y0,y1,y2
integer :: n_O2(0:2)
type(drug_type), pointer :: dp
integer :: nt = 20
real(REAL_KIND) :: tol = 1.0e-12
logical :: use_analytic = .false.

dtt = dt/nt
cp => cell_list(kcell)
ictyp = cp%celltype
Vin = cp%V
CO2 = cp%Cin(OXYGEN)
!if (kcell == 1) then
!	write(*,'(a,i4,3e12.3)') 'integrate_cell_Cin: kcell,CO2,Vin,dt: ',kcell,CO2,Vin,dt
!	write(*,'(a,3e12.3)') 'Cin: ',cp%Cin(DRUG_A:DRUG_A+2)
!	write(*,'(a,3e12.3)') 'Cex: ',cp%Cex(DRUG_A:DRUG_A+2)
!endif
idrug = (ichemo_parent - DRUG_A)/3 + 1
dp => drug(idrug)
do im = 0,2
	ichemo = ichemo_parent + im
	Cin(im) = cp%Cin(ichemo)
	Cex(im) = cp%Cex(ichemo)
	Kin(im) = chemo(ichemo)%membrane_diff_in
	Kout(im) = chemo(ichemo)%membrane_diff_out
	decay_rate(im) = chemo(ichemo)%decay_rate
	Kmet0(im) = dp%Kmet0(ictyp,im)
	CC2 = dp%C2(ictyp,im)
	KO2(im) = dp%KO2(ictyp,im)
	n_O2(im) = dp%n_O2(ictyp,im)
	Kmet0(im) = dp%Kmet0(ictyp,im)
	Km(im) = dp%Km(ictyp,im)
	Vmax(im) = dp%Vmax(ictyp,im)
	K1(im) = (1 - CC2 + CC2*KO2(im)**n_O2(im)/(KO2(im)**n_O2(im) + CO2**n_O2(im)))
!	if (kcell == 1) then
!		write(*,'(2i2,3e12.3)') im,ichemo,Kin(im),Kout(im),Kd(im)
!		write(*,'(i2,3e12.3)') im,Kmet0(im),C2(im),KO2(im)
!		write(*,'(2i2,3e12.3)') im,N_O2(im),Km(im),Vmax(im),K1(im)
!	endif
enddo
if (use_analytic) then
	a0 = Kin(0)*Cex(0)/Vin
	a1 = Kin(1)*Cex(1)/Vin
	a2 = Kin(2)*Cex(2)/Vin
	b0 = Kout(0)/Vin + decay_rate(0)
	b1 = Kout(1)/Vin + decay_rate(1)
	b2 = Kout(2)/Vin + decay_rate(2)
	c0 = K1(0)*Kmet0(0)
!	write(*,*) 'K1,Kmet0: ',K1(0),Kmet0(0)
	c1 = K1(1)*Kmet0(1)
	c2 = K1(2)*Kmet0(2)
	d0 = K1(0)*Vmax(0)
	d1 = K1(1)*Vmax(1)
	d2 = K1(2)*Vmax(2)
	a = b0 + c0
	b = c0*Km(0) + d0 - a0 + b0*Km(0)
	c = -a0*Km(0)
	d = sqrt(b*b - 4*a*c)
!	write(*,*) 'a0,b0,c0: ',a0,b0,c0
!	write(*,*) 'a,b,c,d: ',a,b,c,d
	y0 = (-b + d)/(2*a)
!	write(*,*) 'y0: ',y0
	
	a1 = a1 + (a0 - b0*y0)
	a = b1 + c1
	b = c1*Km(1) + d1 - a1 + b1*Km(1)
	c = -a1*Km(1)
	d = sqrt(b*b - 4*a*c)
!	write(*,*) 'a1,b1,c1: ',a1,b1,c1
!	write(*,*) 'a,b,c,d: ',a,b,c,d
	y1 = (-b + d)/(2*a)
!	write(*,*) 'y1: ',y1
	
	a2 = a2 + (a1 - b1*y1)
	a = b2 + c2 
	b = c2*Km(2) + d2 - a2 + b2*Km(2)
	c = -a2*Km(2)
	d = sqrt(b*b - 4*a*c)
	y2 = (-b + d)/(2*a)
!	write(*,'(a,3e12.4)') 'y: ',y0,y1,y2
	Cin = [y0,y1,y2]
else
	cp%dMdt(ichemo_parent:ichemo_parent+2) = 0
	do it = 1,nt
		K2(0) = Kmet0(0) + Vmax(0)/(Km(0) + Cin(0))
		K2(1) = Kmet0(1) + Vmax(1)/(Km(1) + Cin(1))
		K2(2) = Kmet0(2) + Vmax(2)/(Km(2) + Cin(2))
		dCdt(0) = Kin(0)*Cex(0)/Vin - Kout(0)*Cin(0)/Vin - decay_rate(0)*Cin(0) - K1(0)*K2(0)*Cin(0)
		dCdt(1) = Kin(1)*Cex(1)/Vin - Kout(1)*Cin(1)/Vin - decay_rate(1)*Cin(1) - K1(1)*K2(1)*Cin(1) + K1(0)*K2(0)*Cin(0)
		dCdt(2) = Kin(2)*Cex(2)/Vin - Kout(2)*Cin(2)/Vin - decay_rate(2)*Cin(2) - K1(2)*K2(2)*Cin(2) + K1(1)*K2(1)*Cin(1)
	!	if (kcell == 1) then
	!		write(*,'(a,i4,6e11.4)') 'dCdt: ',it,Cin,dCdt
	!		write(*,'(5e14.6)') Kin(0),Kout(0),decay_rate(0),K1(0),K2(0)
	!		write(*,'(5e14.6)') Cin(0),Cex(0),dCdt(0),Vin,dtt
	!	endif
		do im = 0,2
			Cin(im) = max(0.0, Cin(im) + dCdt(im)*dtt)
		enddo
	enddo
endif
do im = 0,2
	cp%dMdt(ichemo_parent+im) = Kin(im)*Cex(im) - Kout(im)*Cin(im)
	cp%Cin(ichemo_parent+im) = Cin(im)
enddo
!if (kcell == 1) write(*,'(a,3e12.3)') 'Cin: ',Cin(0:2)
!if (kcell == 1) write(*,'(a,2f7.4,3e12.3)') 'integrate_Cin: ',Cex(0),Cin(0),Kin(0),Kout(0),cp%dMdt(ichemo_parent)

end subroutine

!-------------------------------------------------------------------------------------------
! Estimate total flux values associated with each fine grid pt, consistent with the current
! average concentrations Cave in the neighbourhood of the grid pt.  The grid pt flux is
! determined from the uptake fluxes of nearby cells, and a positive computed flux F 
! corresponds to a reduction in grid pt concentration.
! In an iterative procedure, the average grid pt concentrations C are used to estimate
! Cex and Cin for each nearby cell, which then provides an estimate of the effective
! extracellular concentration associated with each grid pt.
! The membrane flux is the equilibrium value, i.e. 
! rate of consumption is balanced by the rate of mass transport across the membrane.  
! Note: It might make sense to compute all the cell concentrations at the same time,
! since the interpolations are identical for each constituent, at a given cell.
! The current best estimate of the extracellular concentrations at the fine grid pts
! is stored in Cextra(:,:,:,ichemo)
! This is initialised to Caverage(:,:,:,:)
! Not necessarily SS
!-------------------------------------------------------------------------------------------
subroutine getF_const(ichemo, Cflux_const)
integer :: ichemo
real(REAL_KIND) :: Cflux_const(:,:,:)
real(REAL_KIND) :: Kin, Kout, total_flux, zmax
integer :: kcell, kcellmax
type(cell_type), pointer :: cp

!write(*,*) 'getF_const: ',ichemo,nlist
Kin = chemo(ichemo)%membrane_diff_in
Kout = chemo(ichemo)%membrane_diff_out

! Compute cell fluxes cp%dMdt
total_flux = 0
zmax = 0
!!$omp parallel do private(cp)
do kcell = 1,nlist
	cp => cell_list(kcell)
	if (cp%state == DEAD) cycle
	if (cp%state == DYING) then
		cp%dMdt(ichemo) = 0
	else
		cp%dMdt(ichemo) = Kin*cp%Cex(ichemo) - Kout*cp%Cin(ichemo)
	endif
!	total_flux = total_flux + cp%dMdt(ichemo)
enddo
!!$omp end parallel do
!write(nflog,'(a,i6,e12.3)') 'kcellmax, zmax: ',kcellmax,zmax

! Estimate grid pt flux values F
call make_grid_flux(ichemo,Cflux_const)
end subroutine

!-------------------------------------------------------------------------------------------
! 1-alfa is the amount of the previous flux 
!-------------------------------------------------------------------------------------------
subroutine make_grid_flux(ichemo,Cflux_const)
integer :: ichemo
real(REAL_KIND) :: Cflux_const(:,:,:)
integer :: kcell, k, cnr(3,8), ix, iy, iz, j, n(3)
real(REAL_KIND) :: wt(8), fsum(3), avec(3,2)
type(cell_type), pointer :: cp
real(REAL_KIND) :: alpha_flux = 0.2		! was 0.3

Cflux_const = 0
fsum = 0
n = 0
avec = 0
do kcell = 1,nlist
	cp => cell_list(kcell)
	if (cp%state == DEAD .or. cp%state == DYING) cycle
	cnr = cp%cnr
	wt = cp%wt
	ix = cp%centre(1,1)/DELTA_X + 1
	if (ix < (NX+1)/2) then
		j = 1
!		do k = 1,8
!			if (cnr(1,k) > (NX+1)/2) then
!				write(*,*) 'Bad cnr: ', ix,k,cnr(1,k),cp%centre(1,1)
!				stop
!			endif
!		enddo
	else
		j = 3
!		do k = 1,8
!			if (cnr(1,k) < (NX+1)/2) then
!				write(*,*) 'Bad cnr: ', ix,k,cnr(1,k),cp%centre(1,1)
!				write(*,*) 'Bad cnr'
!				stop
!			endif
!		enddo
	endif
	n(j) = n(j) + 1
	fsum(j) = fsum(j) + cp%dMdt(ichemo)
	avec(j,1:2) = avec(j,1:2) + cp%Cex(1:2)
!	wt = 1./8	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! try this
	do k = 1,8
		Cflux_const(cnr(1,k),cnr(2,k),cnr(3,k)) = Cflux_const(cnr(1,k),cnr(2,k),cnr(3,k)) + cp%dMdt(ichemo)*wt(k)
	enddo
enddo
!write(*,'(a,2f8.5)') 'Average L concs: ',avec(1,1)/n(1),avec(1,2)/n(1)
!write(*,'(a,2f8.5)') 'Average R concs: ',avec(3,1)/n(3),avec(3,2)/n(3)
!write(*,'(a,i2,3i6,3e11.3)') 'make_grid_flux: fsum(a): ',ichemo,n,fsum
!Cflux_const(:,:,:) = alpha_flux*Cflux_const(:,:,:) + (1-alpha_flux)*Cflux_prev(:,:,:,ichemo) 
fsum = 0
n = 0
do ix = 1,NX
	if (ix < (NX+1)/2) then
		j = 1
	elseif (ix == (NX+1)/2) then
		j = 2
	else
		j = 3
	endif
	do iy = 1,NY
		do iz = 1,NZ
			if (Cflux_const(ix,iy,iz) /= 0) then
				fsum(j) = fsum(j) + Cflux_const(ix,iy,iz)
				Cflux_const(ix,iy,iz) = alpha_flux*Cflux_const(ix,iy,iz) + (1-alpha_flux)*Cflux_prev(ix,iy,iz,ichemo)
			endif
		enddo
	enddo
	n(j) = n(j) + 1
enddo
if (dbug .and. ichemo >= DRUG_A) then
	write(*,*) 'make_grid_flux: total flux: ',ichemo,sum(Cflux_const)
endif
!write(*,'(a,i2,3i6,3e11.3)') 'make_grid_flux: fsum(b): ',ichemo,n,fsum
end subroutine

!-------------------------------------------------------------------------------------------
! Called from simulate_step() after cells have been moved.
! Note that cp%site is the grid pt at the lower left (in 3D) of the cell centre.  
! This implies that the surrounding grid pts all have incremented indices.
!-------------------------------------------------------------------------------------------
subroutine make_grid_flux_weights(ok)
logical :: ok
integer :: ix, iy, iz, kcell, k
real(REAL_KIND) :: c(3)
type(cell_type), pointer :: cp

!!$omp parallel do private(cp,ix,iy,iz) 
do kcell = 1,nlist
	cp => cell_list(kcell)
	if (cp%state == DEAD) cycle
	c = cp%centre(:,1)
	ix = c(1)/DELTA_X + 1
	iy = c(2)/DELTA_X + 1
	iz = c(3)/DELTA_X + 1
	if (ix < 1 .or. ix > NX-1 .or. iy < 1 .or. iy > NY-1 .or. iz < 1 .or. iz > NZ-1) then
		write(logmsg,*) 'make_grid_flux_weights: blob too big, cell outside grid: ',kcell,ix,iy,iz
		call logger(logmsg)
		ok = .false.
		return
	endif
	cp%site = [ix, iy, iz]
	cp%cnr(:,1) = [ix, iy, iz]
	cp%cnr(:,2) = [ix, iy+1, iz]
	cp%cnr(:,3) = [ix, iy, iz+1]
	cp%cnr(:,4) = [ix, iy+1, iz+1]
	cp%cnr(:,5) = [ix+1, iy, iz]
	cp%cnr(:,6) = [ix+1, iy+1, iz]
	cp%cnr(:,7) = [ix+1, iy, iz+1]
	cp%cnr(:,8) = [ix+1, iy+1, iz+1]
!	write(*,*) 'make_grid_flux_weights: ',kcell,cp%state
	call grid_flux_weights(kcell, cp%cnr, cp%wt)
	if (kcell == -1) then
		write(*,'(a,3f8.5)') 'centre, cnr, wt: ',cp%centre(:,1)
		do k = 1,8
			write(*,'(3i4,f8.5)') cp%cnr(:,k),cp%wt(k)
		enddo
		write(*,'(a,f8.5)') 'sum wt: ',sum(cp%wt)
	endif
enddo
!!$omp end parallel do
end subroutine

!-------------------------------------------------------------------------------------------
! This valid for steady-state with any constituent, with decay, when Cex is known
! In fact it is currently valid for oxygen and glucose.
! This is effectively the same as getCin() with dCexdt = 0
!-------------------------------------------------------------------------------------------
function getCin_SS(kcell,ichemo, Vin, Cex) result(Cin)
integer :: kcell, ichemo
real(REAL_KIND) :: Vin, Cex, Cin, Cex_t, Cin_t
real(REAL_KIND) :: Kin, Kout, Kd, Kmax, VKdecay, C0, a, b, c, D, r(3)
integer :: i, n, ictyp, idrug, im
real(REAL_KIND) :: CO2, C_parent, C_metab1
real(REAL_KIND) :: C2(0:2), KO2(0:2), n_O2(0:2), Kmet0(0:2), K1(0:2)
real(REAL_KIND) :: K2, Km, Vmax	!,K1
type(drug_type), pointer :: dp

!write(*,*) 'getCin_SS: istep,ichemo,Vin,Cex: ',istep,ichemo,Vin,Cex
if (Cex < 0) then
	Cex = 0
	Cin = 0
	return
endif
ictyp = cell_list(kcell)%celltype
CO2 = cell_list(kcell)%Cin(OXYGEN)
Kin = chemo(ichemo)%membrane_diff_in
Kout = chemo(ichemo)%membrane_diff_out
Kd = chemo(ichemo)%decay_rate
if (ichemo <= LACTATE) then
	Kmax = chemo(ichemo)%max_cell_rate
	VKdecay = Vin*Kd
	C0 = chemo(ichemo)%MM_C0
	if (chemo(ichemo)%Hill_N == 2) then
		b = C0*C0
		a = (Kmax + b*VKdecay - Kin*Cex)/(Kout+VKdecay)
		c = -b*Cex*Kin/(Kout + VKdecay)
		call cubic_roots(a,b,c,r,n)
		if (n == 1) then
			Cin = r(1)
		else
			n = 0
			do i = 1,3
				if (r(i) > 0) then
					n = n+1
					Cin = r(i)
				endif
			enddo
			if (n > 1) then
				write(nflog,*) 'getCin_SS: two roots > 0: ',r
				stop
			endif
		endif
	elseif (chemo(ichemo)%Hill_N == 1) then
		b = (Kmax + Kout*C0 + 2*VKdecay*C0 - Kin*Cex)/(Kout + VKdecay)
		c = -C0*Cex*Kin/(Kout + VKdecay)
		D = sqrt(b*b - 4*c)
		Cin = (D - b)/2
	endif
elseif (ichemo == TRACER) then
	
else	! parent drug or drug metabolite
	idrug = (ichemo - DRUG_A)/3 + 1
	dp => drug(idrug)
	im = ichemo - DRUG_A - 3*(idrug-1)
	Kmet0(im) = dp%Kmet0(ictyp,im)
	C2(im) = dp%C2(ictyp,im)
	KO2(im) = dp%KO2(ictyp,im)
	n_O2(im) = dp%n_O2(ictyp,im)
	K1(im) = (1 - C2(im) + C2(im)*KO2(im)**n_O2(im)/(KO2(im)**n_O2(im) + CO2**n_O2(im)))*Kmet0(im)
	if (im == 0) then		! parent
		Km = dp%Km(ictyp,0)
		Vmax = dp%Vmax(ictyp,0)
		K2 = K1(0)*Vmax/Kmet0(0)
		if (K2 /= 0) then	!quadratic: a.x^2 + b.x + c = 0
			a = K1(0) + Kd + Kout/Vin
			b = a*Km + K2 - Kin*Cex/Vin
			c = -Kin*Cex*Km/Vin
			b = b/a
			c = c/a
			D = sqrt(b*b - 4*c)
			Cin = (D - b)/2
		else				! linear: a.x + b = 0
			a = K1(0) + Kd + Kout/Vin
			b = -Kin*Cex/Vin
			Cin = -b/a
		endif
	elseif (im == 1) then	! metab1
		CO2 = cell_list(kcell)%Cin(OXYGEN)
		C_parent = cell_list(kcell)%Cin(ichemo-1)
		Cin = (K1(0)*C_parent + Kin*Cex/Vin)/(K1(1) + Kd + Kout/Vin)
	elseif (im == 2) then	! metab2
		CO2 = cell_list(kcell)%Cin(OXYGEN)
		C_metab1 = cell_list(kcell)%Cin(ichemo-1)
		Cin = (K1(1)*C_metab1 + Kin*Cex/Vin)/(K1(2) + Kd + Kout/Vin)
	endif
endif
end function

!-------------------------------------------------------------------------------------------
! Use an approximation to dCin/dt derived from dCex/dt to arrive at a better estimate of Cin and flux
!-------------------------------------------------------------------------------------------
function getCin(kcell, ichemo, Vin, Cex, dCexdt) result(Cin)
integer :: kcell, ichemo
real(REAL_KIND) :: Vin, Cex, dCexdt, Cin
real(REAL_KIND) :: Kin, Kout, Kd, Kmax, VKdecay, dCdt, delta, C0, a, b, c, D, r(3)
integer :: i, n, ictyp, idrug, im
real(REAL_KIND) :: CO2, C_parent, C_metab1
real(REAL_KIND) :: C2(0:2), KO2(0:2), n_O2(0:2), Kmet0(0:2), K1(0:2)
real(REAL_KIND) :: Km, Vmax, K2
type(drug_type), pointer :: dp
type(cell_type), pointer :: cp

if (Cex < 0) then
	Cex = 0
	Cin = 0
	return
endif
cp => cell_list(kcell)
ictyp = cp%celltype
CO2 = cp%Cin(OXYGEN)
Kin = chemo(ichemo)%membrane_diff_in
Kout = chemo(ichemo)%membrane_diff_out
Kd = chemo(ichemo)%decay_rate
if (ichemo <= GLUCOSE) then
	Kmax = chemo(ichemo)%max_cell_rate
	VKdecay = Vin*Kd
	C0 = chemo(ichemo)%MM_C0
	if (chemo(ichemo)%Hill_N == 2) then
		dCdt = dCexdt*(Kin/Vin)/(Kout/Vin + Kd + 2*Kmax*cp%Cin(ichemo)*C0**2/(C0**2+cp%Cin(ichemo)**2))
		delta = Vin*dCdt
		b = C0*C0
		a = (Kmax + b*VKdecay - (Kin*Cex-delta))/(Kout+VKdecay)
		c = -b*(Kin*Cex-delta)/(Kout + VKdecay)
		call cubic_roots(a,b,c,r,n)
		if (n == 1) then
			Cin = r(1)
		else
			n = 0
			do i = 1,3
				if (r(i) > 0) then
					n = n+1
					Cin = r(i)
				endif
			enddo
			if (n > 1) then
				write(nflog,*) 'getCin_SS: two roots > 0: ',r
				stop
			endif
		endif
	elseif (chemo(ichemo)%Hill_N == 1) then
		dCdt = dCexdt*(Kin/Vin)/(Kout/Vin + Kd + Kmax*C0/(C0 + cp%Cin(ichemo)))
		delta = Vin*dCdt
		b = (Kmax + Kout*C0 + 2*VKdecay*C0 - (Kin*Cex-delta))/(Kout + VKdecay)
		c = -C0*(Kin*Cex-delta)/(Kout + VKdecay)
		D = sqrt(b*b - 4*c)
		Cin = (D - b)/2
	endif
elseif (ichemo == TRACER) then
	stop
else	! parent drug or drug metabolite
	idrug = (ichemo - DRUG_A)/3 + 1
	dp => drug(idrug)
	im = ichemo - DRUG_A - 3*(idrug-1)
	Kmet0(im) = dp%Kmet0(ictyp,im)
	C2(im) = dp%C2(ictyp,im)
	KO2(im) = dp%KO2(ictyp,im)
	n_O2(im) = dp%n_O2(ictyp,im)
	K1(im) = (1 - C2(im) + C2(im)*KO2(im)**n_O2(im)/(KO2(im)**n_O2(im) + CO2**n_O2(im)))*Kmet0(im)
	if (im == 0) then		! parent
		Km = dp%Km(ictyp,0)
		Vmax = dp%Vmax(ictyp,0)
		K2 = K1(0)*Vmax/Kmet0(0)
		dCdt = (Kin*dCexdt/Vin)/(Kout/Vin + Kd + K1(0) + K2*Km/(Km + CO2)**2)
		if (K2 /= 0) then	!quadratic: a.x^2 + b.x + c = 0
			a = K1(0) + Kd + Kout/Vin
			b = a*Km + K2 - (Kin*Cex/Vin - dCdt)
			c = -Km*(Kin*Cex/Vin - dCdt)
			b = b/a
			c = c/a
			D = sqrt(b*b - 4*c)
			Cin = (D - b)/2
		else				! linear: a.x + b = 0
			a = K1(0) + Kd + Kout/Vin
			b = -(Kin*Cex/Vin - dCdt)
			Cin = -b/a
		endif
		Cin = max(Cin,0.0)
	elseif (im == 1) then	! metab1
		C_parent = cp%Cin(ichemo-1)
		dCdt = (Kin*dCexdt/Vin + K1(0)*cp%dCdt(ichemo-1))/(Kout/Vin + Kd + K1(1))
		Cin = (K1(0)*C_parent + (Kin*Cex/Vin - dCdt))/(K1(1) + Kd + Kout/Vin)
		Cin = max(Cin,0.0)
	elseif (im == 2) then	! metab2
		C_metab1 = cp%Cin(ichemo-1)
		dCdt = (Kin*dCexdt/Vin + K1(1)*cp%dCdt(ichemo-1))/(Kout/Vin + Kd + K1(2))
		Cin = (K1(1)*C_metab1 + (Kin*Cex/Vin - dCdt))/(K1(2) + Kd + Kout/Vin)
		Cin = max(Cin,0.0)
	endif
endif
end function

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

!--------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------
subroutine unpack_csr(a)
real(REAL_KIND) :: a(:)
integer :: k, irow, icol
real(REAL_KIND) :: asum

write(nflog,*) 'nrow: ',nrow
allocate(afull(nrow,nrow))
irow = 0
do k = 1,nnz
	if (k == ia(irow+1)) irow = irow+1
	icol = ja(k)
	afull(irow,icol) = a(k)
enddo
!do irow = 1,nrow
!	asum = 0
!	do icol = 1,nrow
!		asum = asum + afull(irow,icol)*9
!	enddo
!!	write(*,'(i6,2f8.3)') irow,asum,rhs(irow)
!enddo
!write(*,*) 'rhs:'
!write(*,'(25f7.2)') rhs(:)
end subroutine

!--------------------------------------------------------------------------------------
! Cex values on the boundary of the fine grid are interpolated from Cave_b on the coarse grid
! The coarse grid indices that coincide with the boundary of the fine grid are:
!	ixb = xb0-idxb, xb0+idxb = xb1, xb2
!	iyb = yb0-idyb, yb0+idyb = yb1, yb2
!	izb = zb1
! (see below)
! There are 5 faces.
!
! Average over the faces to obtain an approximation to the average blob boundary conc.
!
! Interpolate over the whole fine grid?
!--------------------------------------------------------------------------------------
subroutine interpolate_Cave(ichemo, Cave, Cave_b)
integer :: ichemo
real(REAL_KIND) :: Cave(:,:,:), Cave_b(:,:,:)
integer :: idx, idy, idz, xb0, yb0, idxb, idyb, xb1, xb2, yb1, yb2, zb1, zb2
integer :: ixb, iyb, izb, ix0, iy0, iz0, ix, iy, iz, nsum, ncsum
real(REAL_KIND) :: ax, ay, az, alpha_conc, asum(2), Cnew, csum, finegrid_Cbnd, cbnd(20000)
real(REAL_KIND) :: cbndx(2)
real(REAL_KIND) :: alpha = 0.5
logical :: use_bdryaverage = .false.

alpha_conc = alpha
if (use_bdryaverage) write(nflog,*) 'interpolate_Cave: use_bdryaverage: ',ichemo
xb0 = (NXB+1)/2			
idxb = (NX-1)/(2*NRF)
xb1 = xb0 - idxb
xb2 = xb0 + idxb
yb0 = (NYB+1)/2
idyb = (NY-1)/(2*NRF)
yb1 = yb0 - idyb
yb2 = yb0 + idyb
zb1 = 1
zb2 = (NZ-1)/NRF + 1

!write(nflog,'(a,6i4)') 'interpolate_Cave: xb1,xb2,...: ',xb1,xb2,yb1,yb2,zb1,zb2
!write(nflog,*) 'Cave_b(14,14,9): ',Cave_b(14,14,9)
csum = 0
ncsum = 0

! Left and right faces: ixb = xb1, xb2  
!                       iyb = yb1..yb2
!						izb = zb1..zb2
asum = 0
nsum = 0
do iyb = yb1,yb2-1
	iy0 = (iyb - yb0)*NRF + (NY+1)/2
	do izb = zb1,zb2-1
		iz0 = (izb - 1)*NRF + 1
		do idy = 0,4
			iy = iy0 + idy
			ay = (4. - idy)/4.
			do idz = 0,4
				iz = iz0 + idz
				az = (4. - idz)/4.
				ixb = xb1
				ix = (ixb - xb0)*NRF + (NX+1)/2
				nsum = nsum + 1
				Cnew =     ay*az*Cave_b(ixb, iyb,   izb) + &
								(1-ay)*az*Cave_b(ixb, iyb+1, izb) + &
								ay*(1-az)*Cave_b(ixb, iyb,   izb+1) + &
							(1-ay)*(1-az)*Cave_b(ixb, iyb+1, izb+1)
				Cave(ix,iy,iz) = alpha_conc*Cnew + (1-alpha_conc)*Cave(ix,iy,iz)
				asum(1) = asum(1) + Cave(ix,iy,iz)
				csum = csum + Cave(ix,iy,iz)
				ncsum = ncsum + 1
				cbnd(ncsum) = Cave(ix,iy,iz)
				if (ichemo == OXYGEN .and. iy == NY/2 .and. iz == NZ/2) then
					cbndx(1) = Cave(ix,iy,iz)
!					write(*,*) 'ix: cbndx(1): ',ix,cbndx(1)
				endif
				
				ixb = xb2
				ix = (ixb - xb0)*NRF + (NX+1)/2
				Cnew =     ay*az*Cave_b(ixb,iyb,  izb) + &
								(1-ay)*az*Cave_b(ixb,iyb+1,izb) + &
								ay*(1-az)*Cave_b(ixb,iyb,  izb+1) + &
							(1-ay)*(1-az)*Cave_b(ixb,iyb+1,izb+1)
				Cave(ix,iy,iz) = alpha_conc*Cnew + (1-alpha_conc)*Cave(ix,iy,iz)
				asum(2) = asum(2) + Cave(ix,iy,iz)
				csum = csum + Cave(ix,iy,iz)
				ncsum = ncsum + 1
				cbnd(ncsum) = Cave(ix,iy,iz)
				if (ichemo == OXYGEN .and. iy == NY/2 .and. iz == NZ/2) then
					cbndx(2) = Cave(ix,iy,iz)
!					write(*,*) 'ix: cbndx(2): ',ix,cbndx(2)
				endif
			enddo
		enddo
	enddo
enddo
!write(*,*) 'Average left and right faces: ',asum/nsum

! Back and front faces: iyb = yb1, yb2
!                       ixb = xb1..xb2
!						izb = zb1..zb2
asum = 0
nsum = 0
do ixb = xb1,xb2-1
	ix0 = (ixb - xb0)*NRF + (NX+1)/2
	do izb = zb1,zb2-1
		iz0 = (izb - 1)*NRF + 1
		do idx = 0,4
			ix = ix0 + idx
			ax = (4. - idx)/4.
			do idz = 0,4
				iz = iz0 + idz
				az = (4. - idz)/4.
				iyb = yb1
				iy = (iyb - yb0)*NRF + (NY+1)/2
				nsum = nsum + 1
				Cnew =     ax*az*Cave_b(ixb,   iyb, izb) + &
					   (1-ax)*az*Cave_b(ixb+1, iyb, izb) + &
					   ax*(1-az)*Cave_b(ixb,   iyb, izb+1) + &
				   (1-ax)*(1-az)*Cave_b(ixb+1, iyb, izb+1)
				Cave(ix,iy,iz) = alpha_conc*Cnew + (1-alpha_conc)*Cave(ix,iy,iz)
				asum(1) = asum(1) + Cave(ix,iy,iz)
				csum = csum + Cave(ix,iy,iz)
				ncsum = ncsum + 1
				cbnd(ncsum) = Cave(ix,iy,iz)

				iyb = yb2
				iy = (iyb - yb0)*NRF + (NY+1)/2
				Cnew =     ax*az*Cave_b(ixb,   iyb, izb) + &
					   (1-ax)*az*Cave_b(ixb+1, iyb, izb) + &
					   ax*(1-az)*Cave_b(ixb,   iyb, izb+1) + &
				   (1-ax)*(1-az)*Cave_b(ixb+1, iyb, izb+1)
				Cave(ix,iy,iz) = alpha_conc*Cnew + (1-alpha_conc)*Cave(ix,iy,iz)
				asum(2) = asum(2) + Cave(ix,iy,iz)
				csum = csum + Cave(ix,iy,iz)
				ncsum = ncsum + 1
				cbnd(ncsum) = Cave(ix,iy,iz)
			enddo
		enddo
	enddo
enddo
!write(*,*) 'Average back and front faces: ',asum/nsum

! Top face: izb = zb2
!           ixb = xb1..xb2
!			iyb = yb1..yb2
nsum = 0
asum = 0
do ixb = xb1,xb2-1
	ix0 = (ixb - xb0)*NRF + (NX+1)/2
!	write(*,*) 'ix0: ',ix0
	do iyb = yb1,yb2-1
		iy0 = (iyb - yb0)*NRF + (NY+1)/2
!		write(*,*) 'iy0: ',iy0
		do idx = 0,4
			ix = ix0 + idx
			ax = (4. - idx)/4.
			do idy = 0,4
				iy = iy0 + idy
				ay = (4. - idy)/4.
				izb = zb2
				iz = (izb - 1)*NRF + 1
				nsum = nsum + 1
				Cnew =		ax*ay*Cave_b(ixb,   iyb,   izb) + &
						(1-ax)*ay*Cave_b(ixb+1, iyb,   izb) + &
						ax*(1-ay)*Cave_b(ixb,   iyb+1, izb) + &
					(1-ax)*(1-ay)*Cave_b(ixb+1, iyb+1, izb)
				Cave(ix,iy,iz) = alpha_conc*Cnew + (1-alpha_conc)*Cave(ix,iy,iz)
				csum = csum + Cave(ix,iy,iz)
				ncsum = ncsum + 1
				cbnd(ncsum) = Cave(ix,iy,iz)
			enddo
		enddo
	enddo
enddo
chemo(ichemo)%fine_grid_cbnd = csum/ncsum		! average concentration on the boundary of the fine grid.
if (dbug .and. ichemo >= DRUG_A) then
    write(*,'(a,i2,i6,e12.3)') 'finegrid_Cbnd (bdry of fine grid): ',ichemo,ncsum,chemo(ichemo)%fine_grid_cbnd
endif
! Try putting average conc everywhere on the boundary to see if it still
! results in asymetry of the profile at low O2
if (use_bdryaverage) then
	write(*,'(a,i2,f8.5)') 'finegrid bdry average: ',ichemo,chemo(ichemo)%fine_grid_cbnd
	do iyb = yb1,yb2-1
		iy0 = (iyb - yb0)*NRF + (NY+1)/2
		do izb = zb1,zb2-1
			iz0 = (izb - 1)*NRF + 1
			do idy = 0,4
				iy = iy0 + idy
				do idz = 0,4
					iz = iz0 + idz
					ixb = xb1
					ix = (ixb - xb0)*NRF + (NX+1)/2
					Cave(ix,iy,iz) = chemo(ichemo)%fine_grid_cbnd
					ixb = xb2
					ix = (ixb - xb0)*NRF + (NX+1)/2
					Cave(ix,iy,iz) = chemo(ichemo)%fine_grid_cbnd
				enddo
			enddo
		enddo
	enddo
	do ixb = xb1,xb2-1
		ix0 = (ixb - xb0)*NRF + (NX+1)/2
		do izb = zb1,zb2-1
			iz0 = (izb - 1)*NRF + 1
			do idx = 0,4
				ix = ix0 + idx
				do idz = 0,4
					iz = iz0 + idz
					iyb = yb1
					iy = (iyb - yb0)*NRF + (NY+1)/2
					Cave(ix,iy,iz) = chemo(ichemo)%fine_grid_cbnd
					iyb = yb2
					iy = (iyb - yb0)*NRF + (NY+1)/2
					Cave(ix,iy,iz) = chemo(ichemo)%fine_grid_cbnd
				enddo
			enddo
		enddo
	enddo
	do ixb = xb1,xb2-1
		ix0 = (ixb - xb0)*NRF + (NX+1)/2
		do iyb = yb1,yb2-1
			iy0 = (iyb - yb0)*NRF + (NY+1)/2
			do idx = 0,4
				ix = ix0 + idx
				do idy = 0,4
					iy = iy0 + idy
					izb = zb2
					iz = (izb - 1)*NRF + 1
					Cave(ix,iy,iz) = chemo(ichemo)%fine_grid_cbnd
				enddo
			enddo
		enddo
	enddo				
endif

!write(nflog,'(a,i4)') 'interpolate_Cave: ichemo, Cave: ',ichemo
!write(nflog,'(10f8.2)') Cave(NX/2,NY/2,:)
end subroutine

!--------------------------------------------------------------------------------------
! Estimate the total mass of a constituent.
! Add the extracellular contributions from Cave_b and the intracellular contributions
! from cell_list.
! Note that there is double counting here because of the micky-mouse way that Cex is
! treated inside the blob.  Not important right now for our purposes.
! Since V is cm3 and C is muM, mass is mumols
!--------------------------------------------------------------------------------------
subroutine getMass(ichemo,mass)
integer :: ichemo
real(REAL_KIND) :: mass, cmass
integer :: ixb, iyb, izb, kcell
type(cell_type),pointer :: cp

mass = 0
do ixb = 1,NXB
	do iyb = 1,NYB
		do izb = 1,NZB
			mass = mass + dxb3*chemo(ichemo)%Cave_b(ixb,iyb,izb)
		enddo
	enddo
enddo
cmass = 0
do kcell = 1,nlist
	cp => cell_list(kcell)
	if (cp%state == DEAD) cycle
	cmass = cmass + cp%V*cp%Cin(ichemo)
enddo
if (dbug .and. ichemo >= DRUG_A) write(*,'(a,i2,2e12.3)') 'medium,cell mass: ',ichemo,mass,cmass
mass = mass + cmass
end subroutine

!--------------------------------------------------------------------------------------
! Use the solution for Cave to estimate the average blob boundary concentration: %medium_Cbnd
! by averaging over a sphere with radius = blob radius (cm) centred at the blob centre
! centre_b -> centre_f on the fine grid. (unless this changes significantly)
! To generate a random point on the sphere, it is necessary only to generate two random numbers, z between -R and R, phi between 0 and 2 pi, each with a uniform distribution
! To find the latitude (theta) of this point, note that z=Rsin(theta), so theta=sin-1(z/R); its longitude is (surprise!) phi.
! In rectilinear coordinates, x=Rcos(theta)cos(phi), y=Rcos(theta)sin(phi), z=Rsin(theta)= (surprise!) z.
! Needless to say, (x,y,z) are not independent, but are rather constrained by x2+y2+z2=R2.
! medium_cbnd is not used anywhere, just for info.
!--------------------------------------------------------------------------------------
subroutine set_bdry_conc
integer :: kpar = 0
real(REAL_KIND) :: rad, xb0, yb0, zb0, x, y, z, p(3), phi, theta, c(MAX_CHEMO), csum(MAX_CHEMO)
real(REAL_KIND) :: xf0, yf0, zf0    ! centre on fine grid
real(REAL_KIND) :: pmax(3), cmax
integer :: i, ic, ichemo, ix, iy, iz, k, n = 100

!call SetRadius(Nsites)
rad = blobradius
xf0 = blobcentre(1)
yf0 = blobcentre(2)
zf0 = blobcentre(3)
csum = 0
cmax = 0
do
	z = -rad + 2*rad*par_uni(kpar)
	phi = 2*PI*par_uni(kpar)
	theta = asin(z/rad)
	x = xf0 + rad*cos(theta)*cos(phi)
	y = yf0 + rad*cos(theta)*sin(phi)
	z = zf0 + z
	ix = x/DXF + 1
	iy = y/DXF + 1
	iz = z/DXF + 1
	if (ix < 1 .or. ix+1 > NX) cycle
	if (iy < 1 .or. iy+1 > NY) cycle
	if (iz < 1 .or. iz+1 > NZ) cycle
	k = k+1
	p = [x, y, z]
	call getConc_f(p,c)
	csum = csum + c
	if (c(1) > cmax) then
	    cmax = c(1)
	    pmax = p
	endif
	if (k == n) exit
enddo
do ic = 1,nchemo
	ichemo = chemomap(ic)
	chemo(ichemo)%medium_Cbnd = csum(ichemo)/n
enddo
end subroutine

!--------------------------------------------------------------------------------------
! Interpolate to obtain concentrations at p(:) = (x,y,z) from Caverage(:,:,:,:)
! Copied from spheroid_abm, which uses the coarse grid solution chemo(ichemo)%Cave_b(:,:,:)  
! Here we want to use the fine grid, which has more resolution.
!--------------------------------------------------------------------------------------
subroutine getConc_f(cb,c)
real(REAL_KIND) :: cb(3), c(:)
integer :: ix, iy, iz, grid(3), i
real(REAL_KIND) :: alfa(3)
real(REAL_KIND), pointer :: Cextra(:,:,:,:)

ix = cb(1)/DXF + 1
iy = cb(2)/DXF + 1
iz = cb(3)/DXF + 1
grid = [ix, iy, iz]
do i = 1,3
	alfa(i) = (cb(i) - (grid(i)-1)*DXF)/DXF
enddo

c = 0
Cextra => Caverage(:,:,:,:)
c(:) = (1-alfa(1))*(1-alfa(2))*(1-alfa(3))*Cextra(ix,iy,iz,:)  &
		+ (1-alfa(1))*alfa(2)*(1-alfa(3))*Cextra(ix,iy+1,iz,:)  &
		+ (1-alfa(1))*alfa(2)*alfa(3)*Cextra(ix,iy+1,iz+1,:)  &
		+ (1-alfa(1))*(1-alfa(2))*alfa(3)*Cextra(ix,iy,iz+1,:)  &
		+ alfa(1)*(1-alfa(2))*(1-alfa(3))*Cextra(ix+1,iy,iz,:)  &
		+ alfa(1)*alfa(2)*(1-alfa(3))*Cextra(ix+1,iy+1,iz,:)  &
		+ alfa(1)*alfa(2)*alfa(3)*Cextra(ix+1,iy+1,iz+1,:)  &
		+ alfa(1)*(1-alfa(2))*alfa(3)*Cextra(ix+1,iy,iz+1,:)
end subroutine

!--------------------------------------------------------------------------------------
! See D:\ITSOL\tests\input
!
! Time for 20 iterations with NX=1000
! ILUtype=1 33 sec
! ILUtype=2 is hopelessly slow with the parameters specified
! ILUtype=3 56 sec
! ILUtype=4 197 sec
!
! Note: Solving for each constituent in parallel with OpenMP like this is possible
!       only if the constituent reactions are independent
! For drugs, the three related constituents are solved for sequentially in the fine grid.
! First drug, then metab1 (using drug results), then metab2 (using metab1 results).
!-------------------------------------------------------------------------------------- 
subroutine diff_solver(dt,ok)
real(REAL_KIND) :: dt
logical :: ok
integer :: i, k, k1, ix, iy, iz, irow, icol, kc, ic, icc, it, ix0, ix1, ix2
integer :: ixb, iyb, izb, izb0
integer :: ichemo, ierr, nfill, iters, maxits, im_krylov, nltz
real(REAL_KIND) :: R, tol, tol_b, asum, t, Vex_curr, Vex_next, Vin_curr, Vin_next, fdecay
real(REAL_KIND) :: Csum, dCsum, msum, mass(MAX_CHEMO), Cmin, masslimit
real(REAL_KIND), allocatable :: x(:), rhs(:)
real(REAL_KIND), pointer :: Cave(:,:,:), Cprev(:,:,:), Fprev(:,:,:), Fcurr(:,:,:)
real(REAL_KIND), pointer :: Cave_b(:,:,:), Cprev_b(:,:,:), Fprev_b(:,:,:), Fcurr_b(:,:,:)
real(REAL_KIND), allocatable :: a(:), a_b(:)
real(REAL_KIND) :: alpha_cave = 0.2
real(REAL_KIND) :: dCtol = 1.0e-4
real(REAL_KIND) :: massmin = 1.0e-8
integer :: ILUtype = 1
integer :: nfinemap, finemap(MAX_CHEMO)
integer :: im, im1, im2, ichemof, iy0, iz0
logical :: zeroF(MAX_CHEMO), zeroC(MAX_CHEMO)
logical :: done
logical :: do_fine = .true.
logical :: use_const = .true.
logical :: test_symmetry = .false.
real(REAL_KIND) :: tdiff, tmetab, t0, t1

if (dbug) write(*,'(a,L2,f8.2)') 'diff_solver: medium_change_step,dt: ',medium_change_step,dt
ok = .true.
do ichemo = 1,MAX_CHEMO
	if (.not.chemo(ichemo)%used) cycle
	call getMass(ichemo,mass(ichemo))
	if (chemo(ichemo)%present) then
		if (chemo(DRUG_A)%present .and. (ichemo == DRUG_A+1 .or. ichemo == DRUG_A+2)) cycle
		if (chemo(DRUG_B)%present .and. (ichemo == DRUG_B+1 .or. ichemo == DRUG_B+2)) cycle
		masslimit = 1.0e-12
		if (mass(ichemo) < masslimit) then
			chemo(ichemo)%present = .false.
			write(nflog,*) 'Set present = false: ',ichemo 
			if (dbug) write(*,*) 'Set present = false: ',ichemo 
		endif
	else
		if (mass(ichemo) > 1.0e-12) then
			chemo(ichemo)%present = .true.
		endif
	endif
enddo
if (dbug) write(*,'(a,3e12.3)') 'DRUG_A mass: ',mass(DRUG_A:DRUG_A+2)
call SetupChemomap

nfinemap = 0
do ic = 1,nchemo
	ichemo = chemomap(ic)
	if (ichemo <= TRACER) then
		nfinemap = nfinemap + 1
		finemap(nfinemap) = ichemo
	elseif (mod(ichemo-DRUG_A,3) == 0) then
		nfinemap = nfinemap + 1
		finemap(nfinemap) = ichemo		! idrug = (ichemo-DRUG_A)/3 + 1
	endif
enddo		

nfill = 1	! Level of fill for ILUK preconditioner
tol = 1.0d-6
tol_b = 1.0d-6
im_krylov = 60	! dimension of Krylov subspace in (outer) FGMRES
maxits = 50

! Compute all steady-state grid point fluxes in advance from Cextra(:,:,:,:): Cflux(:,:,:,:)

t0 = mytimer()
!$omp parallel do private(Cave, Fcurr, Cave_b, Cprev_b, Fprev_b, Fcurr_b, a_b, x, rhs, ix, iy, iz, ixb, iyb, izb, it, done, ichemo, icc, k, dCsum, msum, iters, ierr)
do ic = 1,nchemo
	ichemo = chemomap(ic)
	if (chemo(ichemo)%constant) cycle
	!write(*,'(a,i2)') 'coarse grid: ichemo: ',ichemo
	ichemo_curr = ichemo
	icc = ichemo - 1
	allocate(rhs(nrow_b))
	allocate(x(nrow_b))
	allocate(a_b(MAX_CHEMO*nrow_b))
	Fcurr => Cflux(:,:,:,ichemo)
	Cave_b => chemo(ichemo)%Cave_b
	Cprev_b => chemo(ichemo)%Cprev_b
	Fprev_b => chemo(ichemo)%Fprev_b
	Fcurr_b => chemo(ichemo)%Fcurr_b
	Cave => Caverage(:,:,:,ichemo)
	Fprev_b = Fcurr_b
	call makeF_b(ichemo,Fcurr_b, Fcurr, dt,zeroF(ichemo))
    izb0 = 5    ! as in spheroid-abm 
!	if (ichemo == OXYGEN) then
!	    write(nflog,*) 'total flux_f: O2: ',sum(Fcurr(:,:,:))
!	    write(nflog,*) 'total flux_b: O2: ',sum(Fcurr_b(:,:,:))
!	    write(nflog,*) 'Cave_b: O2:'
!	    write(nflog,'(10e12.3)') Cave_b(NXB/2,:,izb0)
!	    write(nflog,*) 'Cave_f: O2:'
!	    write(nflog,'(10e12.3)') Cave(NX/2,:,NX/2)
!	endif
!	if (ichemo == GLUCOSE) then
!	    write(nflog,*) 'Cave_b: glucose: ixb,..,izb: ',NXB/2,izb0
!	    write(nflog,'(10e12.3)') Cave_b(NXB/2,:,izb0)
!	endif

	call make_csr_b(a_b, ichemo, dt, Cave_b, Cprev_b, Fcurr_b, Fprev_b, rhs, zeroC(ichemo))		! coarse grid

	! Solve Cave_b(t+dt) on coarse grid
	!----------------------------------
	call itsol_create_matrix(icc,nrow_b,nnz_b,a_b,ja_b,ia_b,ierr)
	!write(nflog,*) 'itsol_create_matrix: ierr: ',ierr
	if (ierr /= 0) then
		ok = .false.
	endif
		
	if (ILUtype == 1) then
		call itsol_create_precond_ILUK(icc,nfill,ierr)
	!	write(nflog,*) 'itsol_create_precond_ILUK: ierr: ',ierr 
	elseif (ILUtype == 2) then
		call itsol_create_precond_VBILUK(icc,nfill,ierr)
	!	write(*,*) 'itsol_create_precond_VBILUK: ierr: ',ierr 
	elseif (ILUtype == 3) then
		call itsol_create_precond_ILUT(icc,nfill,tol_b,ierr)
	!	write(*,*) 'itsol_create_precond_ILUT: ierr: ',ierr 
	elseif (ILUtype == 4) then
		call itsol_create_precond_ARMS(icc,nfill,tol_b,ierr)
	!	write(*,*) 'itsol_create_precond_ARMS: ierr: ',ierr 
	endif
	if (ierr /= 0) then
		ok = .false.
	endif

	do izb = 1,NZB
		do iyb = 1,NYB
			do ixb = 1,NXB
				k = (ixb-1)*NYB*NZB + (iyb-1)*NZB + izb
				x(k) = Cave_b(ixb,iyb,izb)		! initial guess
			enddo
		enddo
	enddo
	if (.not.zeroC(ichemo)) then
	!	write(nflog,*) 'call itsol_solve_fgmr_ILU'
		call itsol_solve_fgmr_ILU(icc, rhs, x, im_krylov, maxits, tol_b, iters, ierr)
	!	write(nflog,*) 'itsol_solve_fgmr_ILU: Cave_b: ierr, iters: ',ierr,iters
		if (ierr /= 0) then
			ok = .false.
		endif
	else
		write(nflog,*) 'no solve, zeroC: ',ichemo
	endif
	call itsol_free_precond_ILU(icc, ierr)
!	write(nflog,*) 'did itsol_free_precond_ILU'
	call itsol_free_matrix(icc, ierr)
!	write(nflog,*) 'did itsol_free_matrix'

	Cprev_b = Cave_b
	fdecay = exp(-chemo(ichemo)%decay_rate*dt)		!1 - chemo(ichemo)%decay_rate*dt
	msum = 0
	nltz = 0
	do izb = 1,NZB
		do iyb = 1,NYB
			do ixb = 1,NXB
				k = (ixb-1)*NYB*NZB + (iyb-1)*NZB + izb
				msum = msum + x(k)*dxb3		! this sums the mass of constituent in mumols
				if (x(k) < 0) then
					nltz = nltz+1
!					write(*,'(a,4i3,e12.3)') 'Cave_b < 0: ',ichemo,ixb,iyb,izb,x(k)
					x(k) = 0
				endif 
				Cave_b(ixb,iyb,izb) = fdecay*x(k)
			enddo
		enddo
	enddo
	if (dbug .and. nltz > 0) write(*,*) 'coarse grid n < zero: ',ichemo,nltz
	! interpolate Cave_b on fine grid boundary
	call interpolate_Cave(ichemo, Cave, Cave_b)
	deallocate(a_b, x, rhs)
enddo
!$omp end parallel do

if (test_symmetry) then
	write(*,*) 'test_symmetry = .true.'
	! Make Caverage and Cflux for Oxygen symmetric about ix = (NX+1)/2
	ichemo = OXYGEN
	Cave => Caverage(:,:,:,ichemo)
	Fcurr => Cflux(:,:,:,ichemo)
	ix0 = (NX+1)/2
	do iy = 1,NY
		do iz = 1,NZ
			do ix1 = 1,ix0-1
				ix2 = NX - ix1 + 1
!				Cave(ix1,iy,iz) = (Cave(ix1,iy,iz) + Cave(ix2,iy,iz))/2
!				Cave(ix2,iy,iz) = Cave(ix1,iy,iz)
				Fcurr(ix1,iy,iz) = (Fcurr(ix1,iy,iz) + Fcurr(ix2,iy,iz))/2
				Fcurr(ix2,iy,iz) = Fcurr(ix1,iy,iz)
			enddo
		enddo
	enddo
endif
	
!$omp parallel do private(Cave, Cprev, Fprev, Fcurr, a, x, rhs, ix, iy, iz, ixb, iyb, izb, it, done, ichemo, icc, k, Csum, dCsum, msum, iters, ichemof, im, im1, im2, ierr)
do ic = 1,nfinemap
	ichemof = finemap(ic)
	if (chemo(ichemof)%constant) cycle
	allocate(rhs(nrow))
	allocate(x(nrow))
	allocate(a(MAX_CHEMO*nrow))
	im1 = 0
	if (ichemof <= TRACER) then
		im2 = 0
	else
		im2 = 2
	endif
	do im = im1, im2	! this is to ensure that the sequence is parent -> metab1 -> metab2
		ichemo = ichemof + im
		if (.not.chemo(ichemo)%present) cycle
		ichemo_curr = ichemo
		icc = ichemo - 1
		Cave => Caverage(:,:,:,ichemo)
		Cprev => chemo(ichemo)%Cprev
		Fcurr => Cflux(:,:,:,ichemo)
		
		if (dbug) write(*,*) 'fine grid: ichemo: ',ichemo
!		write(nflog,*) 'fine grid: ichemo: ',ichemo
!		write(nflog,'(a,i4)') 'Cave: '
!		write(nflog,'(10f8.2)') Cave(NX/2,NY/2,:)

!		if (mass(ichemo) > massmin) then	! remove this check for low levels

		call make_csr_SS(a, ichemo, Cave, Fcurr, rhs)	! fine grid - note: using the same flux values as the Cave_b solution!
		
		! Solve Cave(t+dt) steady-state on fine grid
		!-------------------------------------------
		call itsol_create_matrix(icc,nrow,nnz,a,ja,ia,ierr)
		!write(*,*) 'itsol_create_matrix: ierr: ',ierr
		if (ierr /= 0) then
			ok = .false.
		endif
		if (ILUtype == 1) then
			call itsol_create_precond_ILUK(icc,nfill,ierr)
		!	write(*,*) 'itsol_create_precond_ILUK: ierr: ',ierr 
		elseif (ILUtype == 2) then
			call itsol_create_precond_VBILUK(icc,nfill,ierr)
		!	write(*,*) 'itsol_create_precond_VBILUK: ierr: ',ierr 
		elseif (ILUtype == 3) then
			call itsol_create_precond_ILUT(icc,nfill,tol,ierr)
		!	write(*,*) 'itsol_create_precond_ILUT: ierr: ',ierr 
		elseif (ILUtype == 4) then
			call itsol_create_precond_ARMS(icc,nfill,tol,ierr)
		!	write(*,*) 'itsol_create_precond_ARMS: ierr: ',ierr 
		endif
		if (ierr /= 0) then
			ok = .false.
		endif

		do iz = 1,NZ-1
			do iy = 2,NY-1
				do ix = 2,NX-1
					k = (ix-2)*(NY-2)*(NZ-1) + (iy-2)*(NZ-1) + iz
					x(k) = Cave(ix,iy,iz)		! initial guess
				enddo
			enddo
		enddo
		
		if (.not.zeroC(ichemo)) then
!			write(*,*) 'call itsol_solve_fgmr_ILU'
!			write(nflog,*) 'before fgmr: istep,ichemo: ',istep,ichemo
!			write(nflog,*) 'x:'
!			write(nflog,'(10e12.3)') x(1:100)
!			write(nflog,*) 'rhs:'
!			write(nflog,'(10e12.3)') rhs(1:100)
			call itsol_solve_fgmr_ILU(icc,rhs, x, im_krylov, maxits, tol, iters, ierr)
!			write(*,*) 'itsol_solve_fgmr_ILU: Cave: ierr, iters: ',ierr,iters
!			write(nflog,*) 'x:'
!			write(nflog,'(10e12.3)') x(1:100)
			if (ierr /= 0) then
				write(logmsg,*) 'itsol_solve_fgmr_ILU failed with err: ',ierr
				call logger(logmsg)
				ok = .false.
			endif
		else
			write(nflog,*) 'no solve, zeroC: ',ichemo
		endif
		call itsol_free_precond_ILU(icc,ierr)
		call itsol_free_matrix(icc,ierr)
		
!		dCsum = 0
		Cmin = 100
		Csum = 0
		nltz = 0
		do iz = 1,NZ-1
			do iy = 2,NY-1
				do ix = 2,NX-1
					k = (ix-2)*(NY-2)*(NZ-1) + (iy-2)*(NZ-1) + iz
					if (x(k) < 0) then
						nltz = nltz+1
!						write(*,'(a,4i3,e12.3)') 'fine grid x(k)<0 : ',ichemo,ix,iy,iz,x(k)
						x(k) = 0
					endif
!					Cave(ix,iy,iz) = max(0.0,alpha_cave*x(k) + (1-alpha_cave)*Cave(ix,iy,iz))
					Cave(ix,iy,iz) = alpha_cave*x(k) + (1-alpha_cave)*Cave(ix,iy,iz)
					Csum = Csum + Cave(ix,iy,iz)
					Cmin = min(Cave(ix,iy,iz),Cmin)
!						dCsum = dCsum + abs((Cave(ix,iy,iz) - Cprev(ix,iy,iz))/Cave(ix,iy,iz))
				enddo
			enddo
		enddo
		if (dbug .and. nltz > 0 .and. ichemo >= DRUG_A) write(*,*) 'zzzzzzzzzzzzz fine grid n < zero: zzzzzzzzzzzzzzzz: ',ichemo,nltz
		if (dbug .and. ichemo == DRUG_A+1) then
			write(nflog,'(a,i2,3e12.3)') 'ichemo: min, ave, mid: ',ichemo,Cmin,Csum/(NX*NY*NZ),Cave(NX/2,NY/2,NZ/2)
			write(nflog,'(40e11.3)') Cave(NX/2,NY/2,:)
		endif
		
!		endif
		
		Cflux_prev(:,:,:,ichemo) = Fcurr
		
!		if (use_SS) then
!			call update_Cin_const_SS(ichemo)				! true steady-state
!			call getF_const(ichemo, Fcurr)
!		elseif (use_integration) then       ! uses SS for O2, glucose. Now assume use_integration always
			if (.not.use_metabolism) then
				! This updates Cex and Cin
				call update_Cex(ichemo)
			endif
!		else
!			call update_Cex_Cin_dCdt_const(ichemo,dt)		! quasi steady-state
!			call getF_const(ichemo, Fcurr)
!		endif
		Cprev = Cave
							
	enddo
	deallocate(a, x, rhs)
	
enddo
!$omp end parallel do

t1 = mytimer()
tdiff = t1 - t0

! This updates Cex, Cin and dMdt, grid fluxes
if (use_metabolism) then
	call update_IC
	if (Ncells > nspeedtest) then
		call speedtest
		nspeedtest = nspeedtest + 10000
	endif
endif

tmetab = mytimer() - t1
!write(logmsg,'(a,3f8.4)') 'tdiff,tmetab,tmetab/(tdiff+tmetab): ',tdiff,tmetab,tmetab/(tdiff+tmetab)
!call logger(logmsg)

!if (use_integration) then	! always
if (chemo(DRUG_A)%present .or. chemo(DRUG_B)%present) then
	! solve for Cin and dMdt for drug + metabolites by integrating them together
	call integrate_Cin(dt)
endif
do ic = 1,nchemo
	ichemo = chemomap(ic)
	Fcurr => Cflux(:,:,:,ichemo)
	if (ichemo < DRUG_A) then
		if (.not.use_metabolism) then
			call getF_const(ichemo, Fcurr)	! this computes dMdt (SS cell fluxes) then calls make_grid_flux
		endif
	else
		! use the cell fluxes dMdt previously computed in integrate_Cin
		call make_grid_flux(ichemo,Fcurr)
	endif
!	if (ichemo == OXYGEN) then
!		iy0 = (NY+1)/2
!		iz0 = (NZ+1)/2
!		write(*,'(7e11.3)') Fcurr(:,iy0,iz0)
!	endif
enddo
!endif

end subroutine

!----------------------------------------------------------------------------------
!----------------------------------------------------------------------------------
subroutine CheckDrugConcs
integer :: ndrugs_present, drug_present(3*MAX_DRUGTYPES), drug_number(3*MAX_DRUGTYPES)
integer :: idrug, iparent, im, kcell, ichemo, i
type(cell_type), pointer :: cp

ndrugs_present = 0
drug_present = 0
do idrug = 1,ndrugs_used
	iparent = DRUG_A + 3*(idrug-1)
	if (chemo(iparent)%present) then		! simulation with this drug has started
	    do im = 0,2
	        ichemo = iparent + im
	        ndrugs_present = ndrugs_present + 1
	        drug_present(ndrugs_present) = ichemo
	        drug_number(ndrugs_present) = idrug
	    enddo
	endif
enddo

do kcell = 1,nlist
	if (cell_list(kcell)%state == DEAD) cycle
    cp => cell_list(kcell)
	do i = 1,ndrugs_present
	    ichemo = drug_present(i)
	    idrug = drug_number(i)
	    if (cp%Cin(ichemo) > Cthreshold) drug_gt_cthreshold(idrug) = .true.
	    if (cp%Cex(ichemo) > Cthreshold) drug_gt_cthreshold(idrug) = .true.
	enddo
enddo
    
end subroutine

!-------------------------------------------------------------------------------------- 
! Computes Cex, metabolic rates, updates SS Cin, dMdt for each cell, then grid fluxes
! Although the solution time is only 2 sec, this has been shown to be sufficient to 
! reach steady state.
!-------------------------------------------------------------------------------------- 
subroutine update_IC
integer :: kcell, ic, ichemo, ix, iy, iz, it, nt=10
real(REAL_KIND) :: alfa(3), dCdt(3), Kin, Kout, dtt, area_factor, vol_cm3, tstart, y(3), rate(3)
type(cell_type), pointer :: cp
type(metabolism_type), pointer :: mp
real(REAL_KIND), pointer :: Cextra(:,:,:), Fcurr(:,:,:)
real(REAL_KIND) :: average_volume = 1.2
logical :: ok

!write(*,*) 'update_IC: '
area_factor = (average_volume)**(2./3.)

!$omp parallel do private(cp, alfa, ix, iy, iz, ic, ichemo, Cextra, tstart, dtt, Kin, Kout, ok)
do kcell = 1,nlist
	cp => cell_list(kcell)
	if (cp%state == DEAD .or. cp%state == DYING) cycle
	call grid_interp(kcell, alfa)
	ix = cp%site(1)
	iy = cp%site(2)
	iz = cp%site(3)
	do ic = 1,nchemo
		ichemo = chemomap(ic)
		Cextra => Caverage(:,:,:,ichemo)		! currently using the average concentration!	
		cp%Cex(ichemo) = &
			  (1-alfa(1))*(1-alfa(2))*(1-alfa(3))*Cextra(ix,iy,iz)  &
			+ (1-alfa(1))*alfa(2)*(1-alfa(3))*Cextra(ix,iy+1,iz)  &
			+ (1-alfa(1))*(1-alfa(2))*alfa(3)*Cextra(ix,iy,iz+1)  &
			+ (1-alfa(1))*alfa(2)*alfa(3)*Cextra(ix,iy+1,iz+1)  &
			+ alfa(1)*(1-alfa(2))*(1-alfa(3))*Cextra(ix+1,iy,iz)  &
			+ alfa(1)*alfa(2)*(1-alfa(3))*Cextra(ix+1,iy+1,iz)  &
			+ alfa(1)*(1-alfa(2))*alfa(3)*Cextra(ix+1,iy,iz+1) &
			+ alfa(1)*alfa(2)*alfa(3)*Cextra(ix+1,iy+1,iz+1)
	enddo
	tstart = 0
	dtt = 2
	call OGLSolver(kcell,tstart,dtt,ok)
!	do while (cp%Cin(OXYGEN) > cp%Cex(OXYGEN))
!		tstart = tstart + dtt
!		dtt = 0.5
!		call OGLSolver(kcell,tstart,dtt,ok)
!	enddo
	ichemo = GLUCOSE
	if (cp%Cin(ichemo) > cp%Cex(ichemo)) then
		write(nflog,*) 'Cin > Cex: ',kcell,ichemo,cp%Cin(ichemo),cp%Cex(ichemo)
	endif
	do ichemo = OXYGEN,LACTATE
		Kin = chemo(ichemo)%membrane_diff_in
		Kout = chemo(ichemo)%membrane_diff_out
		cp%dMdt(ichemo) = area_factor*(Kin*cp%Cex(ichemo) - Kout*cp%Cin(ichemo))
	enddo
enddo
!$omp end parallel do

do ichemo = OXYGEN,LACTATE
	Fcurr => Cflux(:,:,:,ichemo)
	call make_grid_flux(ichemo,Fcurr)
enddo
end subroutine

!-------------------------------------------------------------------------------------- 
! The aim is to see how time to solve OGL for all cells varies with num_threads
!-------------------------------------------------------------------------------------- 
subroutine speedtest
type(cell_type), allocatable :: cell_save(:)
integer :: kcell, nthreads, ntmax, npr
real(REAL_KIND) :: dtt, tstart, t0
type(cell_type), pointer :: cp
logical :: ok

allocate(cell_save(nlist))
do kcell = 1,nlist
	cell_save(kcell) = cell_list(kcell)
enddo

npr = omp_get_num_procs()
ntmax = omp_get_max_threads()
do nthreads = 1,npr
	do kcell = 1,nlist
		cell_list(kcell) = cell_save(kcell)
	enddo
	call omp_set_num_threads(nthreads)
	t0 = mytimer()
!$omp parallel do private(tstart, dtt, ok)
	do kcell = 1,nlist
		cp => cell_list(kcell)
		if (cp%state == DEAD .or. cp%state == DYING) cycle
		tstart = 0
		dtt = 2
		call OGLSolver(kcell,tstart,dtt,ok)	
	enddo
!$omp end parallel do
	write(logmsg,'(a,i2,e12.4)') 'OGLsolver speed test: nthreads,time: ',nthreads,mytimer() - t0
	call logger(logmsg)
enddo
deallocate(cell_save)
call omp_set_num_threads(ntmax)
end subroutine

!-------------------------------------------------------------------------------------- 
!-------------------------------------------------------------------------------------- 
subroutine ClearDrug(ichemo)
integer :: ichemo
integer :: kcell, im
type(cell_type),pointer :: cp

write(logmsg,*) 'ClearDrug: ',ichemo
call logger(logmsg)
do im = 0,2
	Caverage(:,:,:,ichemo+im) = 0
	Cflux(:,:,:,ichemo+im) = 0
	chemo(ichemo+im)%Cave_b = 0
	chemo(ichemo+im)%Cprev_b = 0
	chemo(ichemo+im)%Fprev_b = 0
	chemo(ichemo+im)%Fcurr_b = 0
enddo
do kcell = 1,nlist
    cp => cell_list(kcell)
	if (cp%state == DEAD .or. cp%state == DYING) cycle
	cp%Cin(ichemo:ichemo+2) = 0
	cp%Cex(ichemo:ichemo+2) = 0
enddo
end subroutine

!-------------------------------------------------------------------------------------- 
! To test fine grid solver with constant boundary condition for oxygen.
! Initially set Cave to constant values everywhere for oxygen and glucose
! (and lactate if use_metabolism).
! The glucose level can remain constant.  The solver changes only interior values on
! the grid, therefore bdry values do not need to be reset.
! Need to set up Cflux(), leave it fixed.
! Simplest to set a constant Cflux value for grid points within a range of distances 
! from the centre.
! The solver finds the steady state concentration field.
! The conclusion is that the solver generates a symmetrical solution from symmetrical
! boundary and source conditions.
!-------------------------------------------------------------------------------------- 
subroutine test_finesolver
integer :: i, k, k1, ix, iy, iz, irow, icol, kc, ic, icc, it
integer :: ixb, iyb, izb, izb0
integer :: ichemo, ierr, nfill, iters, maxits, im_krylov, nltz
real(REAL_KIND) :: R, tol, tol_b, asum, t, Vex_curr, Vex_next, Vin_curr, Vin_next, fdecay
real(REAL_KIND) :: Csum, dCsum, msum, mass(MAX_CHEMO), Cmin
real(REAL_KIND), allocatable :: x(:), rhs(:)
real(REAL_KIND), pointer :: Cave(:,:,:), Cprev(:,:,:), Fprev(:,:,:), Fcurr(:,:,:)
real(REAL_KIND), allocatable :: a(:)
real(REAL_KIND) :: alpha_cave = 0.2
integer :: ILUtype = 1
logical :: ok
real(REAL_KIND) :: xx, yy, zz, R1, R2, R12, R22, c(3), d2, constant_flux, df, factor

write(*,*) 'test_finesolver:'
nfill = 1	! Level of fill for ILUK preconditioner
tol = 1.0d-6
im_krylov = 60	! dimension of Krylov subspace in (outer) FGMRES 
maxits = 50
!!$omp parallel do private(Cave, Cprev, Fprev, Fcurr, a, x, rhs, ix, iy, iz, ixb, iyb, izb, it, done, ichemo, icc, k, Csum, dCsum, msum, iters, ichemof, im, im1, im2, ierr)
if (.not.use_metabolism) then
	write(*,*) 'test_finesolver needs use_metabolism = true'
	stop
endif
ichemo = OXYGEN
ichemo_curr = ichemo
icc = ichemo - 1
Cave => Caverage(:,:,:,ichemo)
Cprev => chemo(ichemo)%Cprev
Fcurr => Cflux(:,:,:,ichemo)

Cave = 0.05
constant_flux = 4.4e-10

! Set grid pt fluxes
c(1) = (NX-1)*DELTA_X/2
c(2) = (NY-1)*DELTA_X/2
c(3) = (NZ-1)*DELTA_X/2
R1 = 5*DELTA_X
R2 = 10*DELTA_X
R12 = R1*R1
R22 = R2*R2

df = 0.002
do ix = 1,NX
	if (ix < (NX+1)/2) then
		factor = 1 - df
	elseif (ix > (NX+1)/2) then
		factor = 1 + df
	else
		factor = 1
	endif
	do iy = 1,NY
		do iz = 1,NZ
			xx = (ix-1)*DELTA_X
			yy = (iy-1)*DELTA_X
			zz = (iz-1)*DELTA_X
			d2 = (xx-c(1))**2 + (yy-c(2))**2 + (zz-c(3))**2
			if (d2 >= R12 .and. d2 <= R22) then
				Fcurr(ix,iy,iz) = factor*constant_flux
			else
				Fcurr(ix,iy,iz) = 0
			endif
		enddo
	enddo
enddo

do it = 1,40

Cprev = Cave

!do ic = 1,nfinemap
!	ichemof = finemap(ic)
!	if (chemo(ichemof)%constant) cycle
	allocate(rhs(nrow))
	allocate(x(nrow))
	allocate(a(MAX_CHEMO*nrow))
!	im1 = 0
!	if (ichemof <= TRACER) then
!		im2 = 0
!	else
!		im2 = 2
!	endif
!	do im = im1, im2	! this is to ensure that the sequence is parent -> metab1 -> metab2
!		ichemo = ichemof + im
!		if (.not.chemo(ichemo)%present) cycle
		
		if (dbug) write(*,*) 'fine grid: ichemo: ',ichemo
!		write(nflog,*) 'fine grid: ichemo: ',ichemo
!		write(nflog,'(a,i4)') 'Cave: '
!		write(nflog,'(10f8.2)') Cave(NX/2,NY/2,:)

!		if (mass(ichemo) > massmin) then	! remove this check for low levels

		call make_csr_SS(a, ichemo, Cave, Fcurr, rhs)	! fine grid - note: using the same flux values as the Cave_b solution!
		
		! Solve Cave(t+dt) steady-state on fine grid
		!-------------------------------------------
		call itsol_create_matrix(icc,nrow,nnz,a,ja,ia,ierr)
		!write(*,*) 'itsol_create_matrix: ierr: ',ierr
		if (ierr /= 0) then
			ok = .false.
		endif
		if (ILUtype == 1) then
			call itsol_create_precond_ILUK(icc,nfill,ierr)
		!	write(*,*) 'itsol_create_precond_ILUK: ierr: ',ierr 
		elseif (ILUtype == 2) then
			call itsol_create_precond_VBILUK(icc,nfill,ierr)
		!	write(*,*) 'itsol_create_precond_VBILUK: ierr: ',ierr 
		elseif (ILUtype == 3) then
			call itsol_create_precond_ILUT(icc,nfill,tol,ierr)
		!	write(*,*) 'itsol_create_precond_ILUT: ierr: ',ierr 
		elseif (ILUtype == 4) then
			call itsol_create_precond_ARMS(icc,nfill,tol,ierr)
		!	write(*,*) 'itsol_create_precond_ARMS: ierr: ',ierr 
		endif
		if (ierr /= 0) then
			ok = .false.
		endif

		do iz = 1,NZ-1
			do iy = 2,NY-1
				do ix = 2,NX-1
					k = (ix-2)*(NY-2)*(NZ-1) + (iy-2)*(NZ-1) + iz
					x(k) = Cave(ix,iy,iz)		! initial guess
				enddo
			enddo
		enddo
		
!		if (.not.zeroC(ichemo)) then
!			write(*,*) 'call itsol_solve_fgmr_ILU'
!			write(nflog,*) 'before fgmr: istep,ichemo: ',istep,ichemo
!			write(nflog,*) 'x:'
!			write(nflog,'(10e12.3)') x(1:100)
!			write(nflog,*) 'rhs:'
!			write(nflog,'(10e12.3)') rhs(1:100)
			call itsol_solve_fgmr_ILU(icc,rhs, x, im_krylov, maxits, tol, iters, ierr)
!			write(*,*) 'itsol_solve_fgmr_ILU: Cave: ierr, iters: ',ierr,iters
!			write(nflog,*) 'x:'
!			write(nflog,'(10e12.3)') x(1:100)
			if (ierr /= 0) then
				write(logmsg,*) 'itsol_solve_fgmr_ILU failed with err: ',ierr
				call logger(logmsg)
				ok = .false.
			endif
!		else
!			write(nflog,*) 'no solve, zeroC: ',ichemo
!		endif
		call itsol_free_precond_ILU(icc,ierr)
		call itsol_free_matrix(icc,ierr)
		
!		dCsum = 0
		Cmin = 100
		Csum = 0
		nltz = 0
		do iz = 1,NZ-1
			do iy = 2,NY-1
				do ix = 2,NX-1
					k = (ix-2)*(NY-2)*(NZ-1) + (iy-2)*(NZ-1) + iz
					if (x(k) < 0) then
						nltz = nltz+1
!						write(*,'(a,4i3,e12.3)') 'fine grid x(k)<0 : ',ichemo,ix,iy,iz,x(k)
						x(k) = 0
					endif
!					Cave(ix,iy,iz) = max(0.0,alpha_cave*x(k) + (1-alpha_cave)*Cave(ix,iy,iz))
					Cave(ix,iy,iz) = alpha_cave*x(k) + (1-alpha_cave)*Cave(ix,iy,iz)
					Csum = Csum + Cave(ix,iy,iz)
					Cmin = min(Cave(ix,iy,iz),Cmin)
!						dCsum = dCsum + abs((Cave(ix,iy,iz) - Cprev(ix,iy,iz))/Cave(ix,iy,iz))
				enddo
			enddo
		enddo
		
!		endif
		
		Cflux_prev(:,:,:,ichemo) = Fcurr
		
		if (.not.use_metabolism) then
			! This updates Cex and Cin
			call update_Cex(ichemo)
		endif
		Cprev = Cave
							
!	enddo
	deallocate(a, x, rhs)
	
!enddo
!!$omp end parallel do

write(*,'(10f7.4)') Cave(:,17,17)

enddo

stop

end subroutine

!-------------------------------------------------------------------------------------- 
!-------------------------------------------------------------------------------------- 
subroutine test_coarsesolver
integer :: i, k, k1, ix, iy, iz, irow, icol, kc, ic, icc, it
integer :: ixb, iyb, izb, ixb0, iyb0, izb0
integer :: ichemo, ierr, nfill, iters, maxits, im_krylov, nltz
real(REAL_KIND) :: R, tol, tol_b, asum, t, Vex_curr, Vex_next, Vin_curr, Vin_next, fdecay
real(REAL_KIND) :: Csum, dCsum, msum, mass(MAX_CHEMO), Cmin, dt
real(REAL_KIND), allocatable :: x(:), rhs(:)
real(REAL_KIND), pointer :: Cave(:,:,:), Cprev(:,:,:), Fprev(:,:,:), Fcurr(:,:,:)
real(REAL_KIND), pointer :: Cave_b(:,:,:), Cprev_b(:,:,:), Fprev_b(:,:,:), Fcurr_b(:,:,:)
real(REAL_KIND), allocatable :: a(:),a_b(:)
real(REAL_KIND) :: alpha_cave = 0.3
integer :: ILUtype = 1
logical :: zeroF(MAX_CHEMO), zeroC(MAX_CHEMO)
logical :: ok
real(REAL_KIND) :: xx, yy, zz, R1, R2, R12, R22, c(3), d2, constant_flux

constant_flux = 1.0e-8
nfill = 1	! Level of fill for ILUK preconditioner
tol_b = 1.0d-6
im_krylov = 60	! dimension of Krylov subspace in (outer) FGMRES
maxits = 50
dt = DELTA_T

ichemo = OXYGEN
ichemo_curr = ichemo
icc = ichemo - 1
zeroC(ichemo) = .false.
allocate(rhs(nrow_b))
allocate(x(nrow_b))
allocate(a_b(MAX_CHEMO*nrow_b))
Fcurr => Cflux(:,:,:,ichemo)
Cave_b => chemo(ichemo)%Cave_b
Cprev_b => chemo(ichemo)%Cprev_b
Fprev_b => chemo(ichemo)%Fprev_b
Fcurr_b => chemo(ichemo)%Fcurr_b
Cave => Caverage(:,:,:,ichemo)

! Set initial Cave_b
Cave_b = 0.18
Cprev_b = 0.18

ixb0 = (NXB+1)/2
iyb0 = (NYB+1)/2
izb0 = 5    ! as in spheroid-abm 

! Create constant flux on the coarse grid
Fcurr_b = 0
do ixb = ixb0-2,ixb0+2
do iyb = iyb0-2,iyb0+2
do izb = izb0-2,izb0+2
	Fcurr_b(ixb,iyb,izb) = constant_flux
enddo
enddo
enddo

write(*,*) 'test_coarsesolver:'
do it = 1,20
write(*,*) 'it: ',it
Fprev_b = Fcurr_b
!call makeF_b(ichemo,Fcurr_b, Fcurr, dt,zeroF(ichemo))
!	if (ichemo == OXYGEN) then
!	    write(nflog,*) 'total flux_f: O2: ',sum(Fcurr(:,:,:))
!	    write(nflog,*) 'total flux_b: O2: ',sum(Fcurr_b(:,:,:))
!	    write(nflog,*) 'Cave_b: O2:'
!	    write(nflog,'(10e12.3)') Cave_b(NXB/2,:,izb0)
!	    write(nflog,*) 'Cave_f: O2:'
!	    write(nflog,'(10e12.3)') Cave(NX/2,:,NX/2)
!	endif
!	if (ichemo == GLUCOSE) then
!	    write(nflog,*) 'Cave_b: glucose: ixb,..,izb: ',NXB/2,izb0
!	    write(nflog,'(10e12.3)') Cave_b(NXB/2,:,izb0)
!	endif

call make_csr_b(a_b, ichemo, dt, Cave_b, Cprev_b, Fcurr_b, Fprev_b, rhs, zeroC(ichemo))		! coarse grid

! Solve Cave_b(t+dt) on coarse grid
!----------------------------------
call itsol_create_matrix(icc,nrow_b,nnz_b,a_b,ja_b,ia_b,ierr)
!write(nflog,*) 'itsol_create_matrix: ierr: ',ierr
if (ierr /= 0) then
	ok = .false.
endif
	
if (ILUtype == 1) then
	call itsol_create_precond_ILUK(icc,nfill,ierr)
!	write(nflog,*) 'itsol_create_precond_ILUK: ierr: ',ierr 
elseif (ILUtype == 2) then
	call itsol_create_precond_VBILUK(icc,nfill,ierr)
!	write(*,*) 'itsol_create_precond_VBILUK: ierr: ',ierr 
elseif (ILUtype == 3) then
	call itsol_create_precond_ILUT(icc,nfill,tol_b,ierr)
!	write(*,*) 'itsol_create_precond_ILUT: ierr: ',ierr 
elseif (ILUtype == 4) then
	call itsol_create_precond_ARMS(icc,nfill,tol_b,ierr)
!	write(*,*) 'itsol_create_precond_ARMS: ierr: ',ierr 
endif
if (ierr /= 0) then
	ok = .false.
endif

do izb = 1,NZB
	do iyb = 1,NYB
		do ixb = 1,NXB
			k = (ixb-1)*NYB*NZB + (iyb-1)*NZB + izb
			x(k) = Cave_b(ixb,iyb,izb)		! initial guess
		enddo
	enddo
enddo
!if (.not.zeroC(ichemo)) then
!	write(nflog,*) 'call itsol_solve_fgmr_ILU'
	call itsol_solve_fgmr_ILU(icc, rhs, x, im_krylov, maxits, tol_b, iters, ierr)
!	write(nflog,*) 'itsol_solve_fgmr_ILU: Cave_b: ierr, iters: ',ierr,iters
	if (ierr /= 0) then
		ok = .false.
	endif
!else
!	write(nflog,*) 'no solve, zeroC: ',ichemo
!endif
call itsol_free_precond_ILU(icc, ierr)
!	write(nflog,*) 'did itsol_free_precond_ILU'
call itsol_free_matrix(icc, ierr)
!	write(nflog,*) 'did itsol_free_matrix'

Cprev_b = Cave_b
fdecay = exp(-chemo(ichemo)%decay_rate*dt)		!1 - chemo(ichemo)%decay_rate*dt
msum = 0
nltz = 0
do izb = 1,NZB
	do iyb = 1,NYB
		do ixb = 1,NXB
			k = (ixb-1)*NYB*NZB + (iyb-1)*NZB + izb
			msum = msum + x(k)*dxb3		! this sums the mass of constituent in mumols
			if (x(k) < 0) then
				nltz = nltz+1
!					write(*,'(a,4i3,e12.3)') 'Cave_b < 0: ',ichemo,ixb,iyb,izb,x(k)
				x(k) = 0
			endif 
			Cave_b(ixb,iyb,izb) = fdecay*x(k)
		enddo
	enddo
enddo
if (dbug .and. nltz > 0) write(*,*) 'coarse grid n < zero: ',ichemo,nltz
! interpolate Cave_b on fine grid boundary
!call interpolate_Cave(ichemo, Cave, Cave_b)

write(*,'(10f7.4)') Cave_b(:,iyb0,izb0)

enddo
deallocate(a_b, x, rhs)

stop
end subroutine

end module

