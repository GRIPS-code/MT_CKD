module co2_continuum
use, intrinsic :: iso_fortran_env, only: real64
use netcdf
use utils
implicit none


public :: carbon_dioxide_continuum


contains


subroutine frnco2(v1c, v2c, dvc, nptc, c, tave, v1ss, v2ss, v1abs, v2abs, path)

  real(kind=real64), intent(out) :: v1c, v2c, dvc
  integer, intent(out) :: nptc
  real(kind=real64), dimension(:), intent(inout) :: c ![cm3].
  real(kind=real64), intent(in) :: tave
  real(kind=real64), intent(out) :: v1ss, v2ss
  real(kind=real64), intent(in) :: v1abs, v2abs
  character(len=*), intent(in) :: path

  integer :: dimid, i, i1, i2, j, ncid, npts, varid
  real(kind=real64) :: dvs, tcor, trat, v1s, v2s
  real(kind=real64), dimension(:), allocatable :: s
  real(kind=real64), dimension(1196:1220) :: tdep_bandhead

  !Read data from netcdf file.
  call nc_check(nf90_open(path, nf90_nowrite, ncid))
  call nc_check(nf90_inq_varid(ncid, "tdep_bandhead", varid))
  call nc_check(nf90_get_var(ncid, varid, tdep_bandhead))
  call nc_check(nf90_inq_varid(ncid, "bfco2", varid))
  call nc_check(nf90_get_att(ncid, varid, "wavenumber_lower_bound", v1s))
  call nc_check(nf90_get_att(ncid, varid, "wavenumber_upper_bound", v2s))
  call nc_check(nf90_get_att(ncid, varid, "wavenumber_resolution", dvs))
  call nc_check(nf90_inq_dimid(ncid, "n3", dimid))
  call nc_check(nf90_inquire_dimension(ncid, dimid, len=npts))
  allocate(s(npts))
  call nc_check(nf90_get_var(ncid, varid, s))
  call nc_check(nf90_close(ncid))

  trat = tave/246._real64

  dvc = dvs
  v1ss = v1s
  v2ss = v2s
  v1c = v1abs - dvc
  v2c = v2abs + dvc

  if (v1c .lt. v1s) then
    i1 = -1
  else
    i1 = int((v1c - v1s)/dvs + 0.01_real64)
  endif

  v1c = v1s + dvs*real(i1 - 1)
  i2 = int((v2c - v1s)/dvs + 0.01_real64)
  nptc = i2 - i1 + 3
  if (nptc .gt. npts) nptc = npts + 4
  v2c = v1c + dvs*real(nptc - 1)

  do j = 1, nptc
    i = i1 + (j - 1)
    c(j) = 0._real64
    if ((i .ge. 1) .and. (i .le. npts)) then
      tcor = 1._real64
      if (i .ge. 1196 .and. i .le. 1220) tcor = (trat)**tdep_bandhead(i)
      c(j) = tcor*s(i) ![cm3].
    endif
  enddo
  deallocate(s)
end subroutine frnco2


!> @brief Carbon dixoide continuum.
subroutine carbon_dioxide_continuum(n_absrb, jrad, wk, grid, tave, pave, path, absrb)

  integer, intent(in) :: n_absrb, jrad
  real(kind=real64), intent(in) :: wk !< Carbon dioxide number density [cm-3].
  real(kind=real64), intent(in) :: tave !< Temperature [K].
  real(kind=real64), intent(in) :: pave !< Pressure [mb].
  type(SpectralGrid), intent(in) :: grid !< Spectral grid.
  character(len=*), intent(in) :: path !< Path to dataset.
  real(kind=real64), dimension(:), intent(inout) :: absrb !< Extinction [cm-1].

  integer :: ncid, varid, nptc, j, jfac, ist, last
  real(kind=real64) :: wco2, v1c, v2c, dvc, v1ss, v2ss, vj, cfac
  real(kind=real64), dimension(n_absrb) :: fco2
  real(kind=real64), dimension(500) :: xfacco2
  real(kind=real64), dimension(6000) :: c

  if (grid%xn .gt. -20._real64 .and. grid%x0 .lt. 10000._real64) then
    call nc_check(nf90_open(path, nf90_nowrite, ncid))
    call nc_check(nf90_inq_varid(ncid, "xfacco2", varid))
    call nc_check(nf90_get_var(ncid, varid, xfacco2))
    call nc_check(nf90_close(ncid))
    fco2(:) = 0._real64 ![cm3].
    wco2 = wk*(pave/p0)*(t0/tave)*1.0e-20_real64 ![cm-3].
    call frnco2(v1c, v2c, dvc, nptc, fco2, tave, v1ss, v2ss, grid%x0, grid%xn, path)
    do j = 1, nptc
      vj = v1c + dvc*real(j - 1)
      cfac = 1._real64 !unitless
      if (vj .ge. 2000._real64 .and. vj .le. 2998._real64) then
        jfac = int((vj - 1998._real64)/2._real64 + 0.00001_real64)
        cfac = xfacco2(jfac)
      endif
      fco2(j) = cfac*fco2(j) ![cm3].
      c(j) = fco2(j)*wco2 !unitless.
      if (jrad .eq. 1) c(j) = c(j)*radfn(vj, tave/radcn2) ![cm-1].
    enddo
    call pre_xint(v1ss, v2ss, grid%x0, grid%dx, grid%n, ist, last)
    call xint(v1c, v2c, dvc, c, 1._real64, grid%x0, grid%dx, absrb, ist, last)
  endif
end subroutine carbon_dioxide_continuum


end module co2_continuum
