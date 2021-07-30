module o3_continuum
use netcdf
implicit none


public :: xo3chp, o3hht0, o3hht1, o3hht2, o3hhuv


contains


subroutine xo3chp(v1c, v2c, dvc, nptc, c0, c1, c2, v1ss, v2ss, v1abs, v2abs, path)

  real, intent(out) :: v1c, v2c, dvc
  integer, intent(out) :: nptc
  real, dimension(:), allocatable, intent(inout) :: c0, c1, c2
  real, intent(out) :: v1ss, v2ss
  real, intent(in) :: v1abs, v2abs
  character(len=*), intent(in) :: path

  integer :: dimid, err, i, i1, i2, j, ncid, npts, varid
  real :: dvs, vj, v1s, v2s
  real, dimension(:), allocatable :: x, y, z

  !Read data from netcdf file.
  err = nf90_open(path, nf90_nowrite, ncid)
  err = nf90_inq_varid(ncid, "x_o2", varid)
  err = nf90_get_att(ncid, varid, "wavenumber_lower_bound", v1s)
  err = nf90_get_att(ncid, varid, "wavenumber_upper_bound", v2s)
  err = nf90_get_att(ncid, varid, "wavenumber_resolution", dvs)
  err = nf90_inq_dimid(ncid, "n7", dimid)
  err = nf90_inquire_dimension(ncid, dimid, len=npts)
  allocate(x(npts), y(npts), z(npts))
  err = nf90_get_var(ncid, varid, x)
  err = nf90_inq_varid(ncid, "y_o2", varid)
  err = nf90_get_var(ncid, varid, y)
  err = nf90_inq_varid(ncid, "z_o2", varid)
  err = nf90_get_var(ncid, varid, z)
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

  allocate(c0(nptc), c1(nptc), c2(nptc))
  do j = 1, nptc
    i = i1 + (j - 1)
    if ((i .lt. 1) .or. (i .gt. npts)) then
      c0(j) = 0.
      c1(j) = 0.
      c2(j) = 0.
    else
!            remove radiation field from diffuse ozone
      vj = v1c + dvc*real(j - 1)
      c0(j) = x(i)/vj
      c1(j) = y(i)/vj
      c2(j) = z(i)/vj
    endif
  enddo
  deallocate(x, y, z)
end subroutine xo3chp


subroutine o3hht0(v1c, v2c, dvc, nptc, c, v1ss, v2ss, v1abs, v2abs, path)

  real, intent(out) :: v1c, v2c, dvc
  integer, intent(out) :: nptc
  real, dimension(:), allocatable, intent(inout) :: c
  real, intent(out) :: v1ss, v2ss
  real, intent(in) :: v1abs, v2abs
  character(len=*), intent(in) :: path

  integer :: dimid, err, i, i1,i2, j, ncid, npts, varid
  real :: dvs, vj, v1s, v2s
  real, dimension(:), allocatable :: s

  !Read data from netcdf file.
  err = nf90_open(path, nf90_nowrite, ncid)
  err = nf90_inq_varid(ncid, "o3_hh0", varid)
  err = nf90_get_att(ncid, varid, "wavenumber_lower_bound", v1s)
  err = nf90_get_att(ncid, varid, "wavenumber_upper_bound", v2s)
  err = nf90_get_att(ncid, varid, "wavenumber_resolution", dvs)
  err = nf90_inq_dimid(ncid, "n8", dimid)
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
    vj = v1c + dvc*real(j - 1)
    c(j) = s(i)/vj
  enddo
  deallocate(s)
end subroutine o3hht0


subroutine o3hht1(v1c, v2c, dvc, nptc, c, v1abs, v2abs, path)

  real, intent(out) :: v1c, v2c, dvc
  integer, intent(out) :: nptc
  real, dimension(:), allocatable, intent(inout) :: c
  real, intent(in) :: v1abs, v2abs
  character(len=*), intent(in) :: path

  integer :: dimid, err, i, i1, i2, j, ncid, npts, varid
  real :: dvs, v1s, v2s
  real, dimension(:), allocatable :: s

  !Read data from netcdf file.
  err = nf90_open(path, nf90_nowrite, ncid)
  err = nf90_inq_varid(ncid, "o3_hh1", varid)
  err = nf90_get_att(ncid, varid, "wavenumber_lower_bound", v1s)
  err = nf90_get_att(ncid, varid, "wavenumber_upper_bound", v2s)
  err = nf90_get_att(ncid, varid, "wavenumber_resolution", dvs)
  err = nf90_inq_dimid(ncid, "n8", dimid)
  err = nf90_inquire_dimension(ncid, dimid, len=npts)
  allocate(s(npts))
  err = nf90_get_var(ncid, varid, s)
  err = nf90_close(ncid)

  dvc = dvs
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
end subroutine o3hht1


subroutine o3hht2(v1c, v2c, dvc, nptc, c, v1abs, v2abs, path)

  real, intent(out) :: v1c, v2c, dvc
  integer, intent(out) :: nptc
  real, dimension(:), allocatable, intent(inout) :: c
  real, intent(in) :: v1abs, v2abs
  character(len=*), intent(in) :: path

  integer :: dimid, err, i, i1, i2, j, ncid, npts, varid
  real :: dvs, v1s, v2s
  real, dimension(:), allocatable :: s

  !Read data from netcdf file.
  err = nf90_open(path, nf90_nowrite, ncid)
  err = nf90_inq_varid(ncid, "o3_hh2", varid)
  err = nf90_get_att(ncid, varid, "wavenumber_lower_bound", v1s)
  err = nf90_get_att(ncid, varid, "wavenumber_upper_bound", v2s)
  err = nf90_get_att(ncid, varid, "wavenumber_resolution", dvs)
  err = nf90_inq_dimid(ncid, "n8", dimid)
  err = nf90_inquire_dimension(ncid, dimid, len=npts)
  allocate(s(npts))
  err = nf90_get_var(ncid, varid, s)
  err = nf90_close(ncid)

  dvc = dvs
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
end subroutine o3hht2


subroutine o3hhuv(v1c, v2c, dvc, nptc, c, v1ss, v2ss, v1abs, v2abs, path)

  real, intent(out) :: v1c, v2c, dvc
  integer, intent(out) :: nptc
  real, dimension(:), allocatable, intent(inout) :: c
  real, intent(out) :: v1ss, v2ss
  real, intent(in) :: v1abs, v2abs
  character(len=*), intent(in) :: path

  integer :: dimid, err, i, i1, i2, j, ncid, npts, varid
  real :: dvs, vj, v1s, v2s
  real, dimension(:), allocatable :: s

  !Read data from netcdf file.
  err = nf90_open(path, nf90_nowrite, ncid)
  err = nf90_inq_varid(ncid, "o3_huv", varid)
  err = nf90_get_att(ncid, varid, "wavenumber_lower_bound", v1s)
  err = nf90_get_att(ncid, varid, "wavenumber_upper_bound", v2s)
  err = nf90_get_att(ncid, varid, "wavenumber_resolution", dvs)
  err = nf90_inq_dimid(ncid, "n9", dimid)
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
    vj = v1c + dvc*real(j - 1)
    c(j) = s(i)/vj
  enddo
  deallocate(s)
end subroutine o3hhuv


end module o3_continuum
