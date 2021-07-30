module h2o_continuum
use netcdf
implicit none


public :: sl296, sl260, frn296


contains


subroutine sl296(v1c, v2c, dvc, nptc, c, v1ss, v2ss, v1abs, v2abs, path)

  real, intent(out) :: v1c, v2c, dvc
  integer, intent(out) :: nptc
  real, dimension(:), allocatable, intent(inout) :: c
  real, intent(out) :: v1ss, v2ss
  real, intent(in) :: v1abs, v2abs
  character(len=*), intent(in) :: path

  integer :: dimid, err, i, i1, i2, j, ncid, npts, varid
  real :: dvs, v1s, v2s
  real, dimension(:), allocatable :: s

  !Read data from netcdf file.
  err = nf90_open(path, nf90_nowrite, ncid)
  err = nf90_inq_varid(ncid, "bs296", varid)
  err = nf90_get_att(ncid, varid, "wavenumber_lower_bound", v1s)
  err = nf90_get_att(ncid, varid, "wavenumber_upper_bound", v2s)
  err = nf90_get_att(ncid, varid, "wavenumber_resolution", dvs)
  err = nf90_inq_dimid(ncid, "n2", dimid)
  err = nf90_inquire_dimension(ncid, dimid, len=npts)
  allocate(s(npts))
  err = nf90_get_var(ncid, varid, s)
  err = nf90_close(ncid)

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
    if ((i .lt. 1) .or. (i .gt. npts)) continue
    c(j) = s(i)
  enddo
  deallocate(s)
end subroutine sl296


subroutine sl260(v1c, v2c, dvc, nptc, c, v1ss, v2ss, v1abs, v2abs, path)

  real, intent(out) :: v1c, v2c, dvc
  integer, intent(out) :: nptc
  real, dimension(:), allocatable, intent(inout) :: c
  real, intent(out) :: v1ss, v2ss
  real, intent(in) :: v1abs, v2abs
  character(len=*), intent(in) :: path

  integer :: dimid, err, i, i1, i2, j, ncid, npts, varid
  real :: dvs, v1s, v2s
  real, dimension(:), allocatable :: s

  !Read data from netcdf file.
  err = nf90_open(path, nf90_nowrite, ncid)
  err = nf90_inq_varid(ncid, "bs260", varid)
  err = nf90_get_att(ncid, varid, "wavenumber_lower_bound", v1s)
  err = nf90_get_att(ncid, varid, "wavenumber_upper_bound", v2s)
  err = nf90_get_att(ncid, varid, "wavenumber_resolution", dvs)
  err = nf90_inq_dimid(ncid, "n2", dimid)
  err = nf90_inquire_dimension(ncid, dimid, len=npts)
  allocate(s(npts))
  err = nf90_get_var(ncid, varid, s)
  err = nf90_close(ncid)

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
    if ((i .lt. 1) .or. (i .gt. npts)) continue
    c(j) = s(i)
  enddo
  deallocate(s)
end subroutine sl260


subroutine frn296(v1c, v2c, dvc, nptc, c, v1ss, v2ss, v1abs, v2abs, path)

  real, intent(out) :: v1c, v2c, dvc
  integer, intent(out) :: nptc
  real, dimension(:), allocatable, intent(inout) :: c
  real, intent(out) :: v1ss, v2ss
  real, intent(in) :: v1abs, v2abs
  character(len=*), intent(in) :: path

  integer :: dimid, err, i, i1, i2, j, ncid, npts, varid
  real :: dvs, v1s, v2s
  real, dimension(:), allocatable :: s

  !Read data from netcdf file.
  err = nf90_open(path, nf90_nowrite, ncid)
  err = nf90_inq_varid(ncid, "bfh2o", varid)
  err = nf90_get_att(ncid, varid, "wavenumber_lower_bound", v1s)
  err = nf90_get_att(ncid, varid, "wavenumber_upper_bound", v2s)
  err = nf90_get_att(ncid, varid, "wavenumber_resolution", dvs)
  err = nf90_inq_dimid(ncid, "n2", dimid)
  err = nf90_inquire_dimension(ncid, dimid, len=npts)
  allocate(s(npts))
  err = nf90_get_var(ncid, varid, s)
  err = nf90_close(ncid)

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
      c(j) = s(i)
    endif
  enddo
  deallocate(s)
end subroutine frn296


end module h2o_continuum
