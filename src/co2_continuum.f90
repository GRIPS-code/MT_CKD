module co2_continuum
use netcdf
implicit none


public :: frnco2


contains


subroutine frnco2(v1c, v2c, dvc, nptc, c, tave, v1ss, v2ss, v1abs, v2abs, path)

  real, intent(out) :: v1c, v2c, dvc
  integer, intent(out) :: nptc
  real, dimension(:), allocatable, intent(inout) :: c
  real, intent(in) :: tave
  real, intent(out) :: v1ss, v2ss
  real, intent(in) :: v1abs, v2abs
  character(len=*), intent(in) :: path

  integer :: dimid, err, i, i1, i2, j, ncid, npts, varid
  real :: dvs, tcor, trat, v1s, v2s
  real, dimension(:), allocatable :: s
  real, dimension(1196:1220) :: tdep_bandhead

!     temparature dependence coefficients for wavenumbers between 2386
!     and 2434. computed based on (line-coupled) continuum coefficients
!     at 250k and 296k, set to unity at t_eff (determined by invariance
!     of calculations in this region for iasi low pwv cases).
  data (tdep_bandhead(i),i=1196,1220)/ &
  & 1.44e-01, 3.61e-01, 5.71e-01, 7.63e-01, 8.95e-01, &
  & 9.33e-01, 8.75e-01, 7.30e-01, 5.47e-01, 3.79e-01, &
  & 2.55e-01, 1.78e-01, 1.34e-01, 1.07e-01, 9.06e-02, &
  & 7.83e-02, 6.83e-02, 6.00e-02, 5.30e-02, 4.72e-02, &
  & 4.24e-02, 3.83e-02, 3.50e-02, 3.23e-02, 3.01e-02/

  !Read data from netcdf file.
  err = nf90_open(path, nf90_nowrite, ncid)
  err = nf90_inq_varid(ncid, "bfco2", varid)
  err = nf90_get_att(ncid, varid, "wavenumber_lower_bound", v1s)
  err = nf90_get_att(ncid, varid, "wavenumber_upper_bound", v2s)
  err = nf90_get_att(ncid, varid, "wavenumber_resolution", dvs)
  err = nf90_inq_dimid(ncid, "n3", dimid)
  err = nf90_inquire_dimension(ncid, dimid, len=npts)
  allocate(s(npts))
  err = nf90_get_var(ncid, varid, s)
  err = nf90_close(ncid)

  trat = tave/246.

  dvc = dvs
  v1ss = v1s
  v2ss = v2s
  v1c = v1abs - dvc
  v2c = v2abs + dvc

  if (v1c .lt. v1s) then
    i1 = -1
  else
    i1 = (v1c - v1s)/dvs + 0.01
  endif

  v1c = v1s + dvs*real(i1 - 1)
  i2 = (v2c - v1s)/dvs + 0.01
  nptc = i2 - i1 + 3
  if (nptc .gt. npts) nptc = npts + 4
  v2c = v1c + dvs*real(nptc - 1)

  allocate(c(nptc))
  do j = 1, nptc
    i = i1 + (j - 1)
    c(j) = 0.
    if ((i .ge. 1) .and. (i .le. npts)) then
      tcor = 1.
      if (i .ge. 1196 .and. i .le. 1220) tcor = (trat)**tdep_bandhead(i)
      c(j) = tcor*s(i)
    endif
  enddo
  deallocate(s)
end subroutine frnco2


end module co2_continuum
