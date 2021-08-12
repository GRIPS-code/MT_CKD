module o3_continuum
use, intrinsic :: iso_fortran_env, only: real64
use netcdf
use utils
implicit none


public :: ozone_cw_continuum, ozone_hh_continuum, ozone_uv_continuum


contains


subroutine xo3chp(v1c, v2c, dvc, nptc, c0, c1, c2, v1ss, v2ss, v1abs, v2abs, path)

  real(kind=real64), intent(out) :: v1c, v2c, dvc
  integer, intent(out) :: nptc
  real(kind=real64), dimension(:), intent(inout) :: c0 ![cm3].
  real(kind=real64), dimension(:), intent(inout) :: c1 ![cm3 K-1].
  real(kind=real64), dimension(:), intent(inout) :: c2 ![cm3 K-2].
  real(kind=real64), intent(out) :: v1ss, v2ss
  real(kind=real64), intent(in) :: v1abs, v2abs
  character(len=*), intent(in) :: path

  integer :: dimid, i, i1, i2, j, ncid, npts, varid
  real(kind=real64) :: dvs, vj, v1s, v2s
  real(kind=real64), dimension(:), allocatable :: x, y, z

  !Read data from netcdf file.
  call nc_check(nf90_open(path, nf90_nowrite, ncid))
  call nc_check(nf90_inq_varid(ncid, "x_o3", varid))
  call nc_check(nf90_get_att(ncid, varid, "wavenumber_lower_bound", v1s))
  call nc_check(nf90_get_att(ncid, varid, "wavenumber_upper_bound", v2s))
  call nc_check(nf90_get_att(ncid, varid, "wavenumber_resolution", dvs))
  call nc_check(nf90_inq_dimid(ncid, "n7", dimid))
  call nc_check(nf90_inquire_dimension(ncid, dimid, len=npts))
  allocate(x(npts), y(npts), z(npts))
  call nc_check(nf90_get_var(ncid, varid, x))
  call nc_check(nf90_inq_varid(ncid, "y_o3", varid))
  call nc_check(nf90_get_var(ncid, varid, y))
  call nc_check(nf90_inq_varid(ncid, "z_o3", varid))
  call nc_check(nf90_get_var(ncid, varid, z))
  call nc_check(nf90_close(ncid))

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

  do j = 1, nptc
    i = i1 + (j - 1)
    if ((i .lt. 1) .or. (i .gt. npts)) then
      c0(j) = 0.
      c1(j) = 0.
      c2(j) = 0.
    else
!            remove radiation field from diffuse ozone
      vj = v1c + dvc*real(j - 1)
      c0(j) = x(i)/vj ![cm3].
      c1(j) = y(i)/vj ![cm3 K-1].
      c2(j) = z(i)/vj ![cm3 K-2].
    endif
  enddo
  deallocate(x, y, z)
end subroutine xo3chp


subroutine o3hht0(v1c, v2c, dvc, nptc, c, v1ss, v2ss, v1abs, v2abs, path)

  real(kind=real64), intent(out) :: v1c, v2c, dvc
  integer, intent(out) :: nptc
  real(kind=real64), dimension(:), intent(inout) :: c ![cm3].
  real(kind=real64), intent(out) :: v1ss, v2ss
  real(kind=real64), intent(in) :: v1abs, v2abs
  character(len=*), intent(in) :: path

  integer :: dimid, i, i1,i2, j, ncid, npts, varid
  real(kind=real64) :: dvs, vj, v1s, v2s
  real(kind=real64), dimension(:), allocatable :: s

  !Read data from netcdf file.
  call nc_check(nf90_open(path, nf90_nowrite, ncid))
  call nc_check(nf90_inq_varid(ncid, "o3_hh0", varid))
  call nc_check(nf90_get_att(ncid, varid, "wavenumber_lower_bound", v1s))
  call nc_check(nf90_get_att(ncid, varid, "wavenumber_upper_bound", v2s))
  call nc_check(nf90_get_att(ncid, varid, "wavenumber_resolution", dvs))
  call nc_check(nf90_inq_dimid(ncid, "n8", dimid))
  call nc_check(nf90_inquire_dimension(ncid, dimid, len=npts))
  allocate(s(npts))
  call nc_check(nf90_get_var(ncid, varid, s))
  call nc_check(nf90_close(ncid))

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

  do j = 1, nptc
    i = i1 + (j - 1)
    c(j) = 0.
    if ((i .lt. 1) .or. (i .gt. npts)) continue
    vj = v1c + dvc*real(j - 1)
    c(j) = s(i)/vj ![cm3].
  enddo
  deallocate(s)
end subroutine o3hht0


subroutine o3hht1(v1c, v2c, dvc, nptc, c, v1abs, v2abs, path)

  real(kind=real64), intent(out) :: v1c, v2c, dvc
  integer, intent(out) :: nptc
  real(kind=real64), dimension(:), intent(inout) :: c ![K-1].
  real(kind=real64), intent(in) :: v1abs, v2abs
  character(len=*), intent(in) :: path

  integer :: dimid, i, i1, i2, j, ncid, npts, varid
  real(kind=real64) :: dvs, v1s, v2s
  real(kind=real64), dimension(:), allocatable :: s

  !Read data from netcdf file.
  call nc_check(nf90_open(path, nf90_nowrite, ncid))
  call nc_check(nf90_inq_varid(ncid, "o3_hh1", varid))
  call nc_check(nf90_get_att(ncid, varid, "wavenumber_lower_bound", v1s))
  call nc_check(nf90_get_att(ncid, varid, "wavenumber_upper_bound", v2s))
  call nc_check(nf90_get_att(ncid, varid, "wavenumber_resolution", dvs))
  call nc_check(nf90_inq_dimid(ncid, "n8", dimid))
  call nc_check(nf90_inquire_dimension(ncid, dimid, len=npts))
  allocate(s(npts))
  call nc_check(nf90_get_var(ncid, varid, s))
  call nc_check(nf90_close(ncid))

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

  do j = 1, nptc
    i = i1 + (j - 1)
    c(j) = 0.
    if ((i .lt. 1) .or. (i .gt. npts)) continue
    c(j) = s(i) ![K-1].
  enddo
  deallocate(s)
end subroutine o3hht1


subroutine o3hht2(v1c, v2c, dvc, nptc, c, v1abs, v2abs, path)

  real(kind=real64), intent(out) :: v1c, v2c, dvc
  integer, intent(out) :: nptc
  real(kind=real64), dimension(:), intent(inout) :: c ![K-2].
  real(kind=real64), intent(in) :: v1abs, v2abs
  character(len=*), intent(in) :: path

  integer :: dimid, i, i1, i2, j, ncid, npts, varid
  real(kind=real64) :: dvs, v1s, v2s
  real(kind=real64), dimension(:), allocatable :: s

  !Read data from netcdf file.
  call nc_check(nf90_open(path, nf90_nowrite, ncid))
  call nc_check(nf90_inq_varid(ncid, "o3_hh2", varid))
  call nc_check(nf90_get_att(ncid, varid, "wavenumber_lower_bound", v1s))
  call nc_check(nf90_get_att(ncid, varid, "wavenumber_upper_bound", v2s))
  call nc_check(nf90_get_att(ncid, varid, "wavenumber_resolution", dvs))
  call nc_check(nf90_inq_dimid(ncid, "n8", dimid))
  call nc_check(nf90_inquire_dimension(ncid, dimid, len=npts))
  allocate(s(npts))
  call nc_check(nf90_get_var(ncid, varid, s))
  call nc_check(nf90_close(ncid))

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

  do j = 1, nptc
    i = i1 + (j - 1)
    c(j) = 0.
    if ((i .lt. 1) .or. (i .gt. npts)) continue
    c(j) = s(i) ![K-2].
  enddo
  deallocate(s)
end subroutine o3hht2


subroutine o3hhuv(v1c, v2c, dvc, nptc, c, v1ss, v2ss, v1abs, v2abs, path)

  real(kind=real64), intent(out) :: v1c, v2c, dvc
  integer, intent(out) :: nptc
  real(kind=real64), dimension(:), intent(inout) :: c ![cm3].
  real(kind=real64), intent(out) :: v1ss, v2ss
  real(kind=real64), intent(in) :: v1abs, v2abs
  character(len=*), intent(in) :: path

  integer :: dimid, i, i1, i2, j, ncid, npts, varid
  real(kind=real64) :: dvs, vj, v1s, v2s
  real(kind=real64), dimension(:), allocatable :: s

  !Read data from netcdf file.
  call nc_check(nf90_open(path, nf90_nowrite, ncid))
  call nc_check(nf90_inq_varid(ncid, "o3_huv", varid))
  call nc_check(nf90_get_att(ncid, varid, "wavenumber_lower_bound", v1s))
  call nc_check(nf90_get_att(ncid, varid, "wavenumber_upper_bound", v2s))
  call nc_check(nf90_get_att(ncid, varid, "wavenumber_resolution", dvs))
  call nc_check(nf90_inq_dimid(ncid, "n9", dimid))
  call nc_check(nf90_inquire_dimension(ncid, dimid, len=npts))
  allocate(s(npts))
  call nc_check(nf90_get_var(ncid, varid, s))
  call nc_check(nf90_close(ncid))

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

  do j = 1, nptc
    i = i1 + (j - 1)
    c(j) = 0.
    if ((i .lt. 1) .or. (i .gt. npts)) continue
    vj = v1c + dvc*real(j - 1)
    c(j) = s(i)/vj ![cm3].
  enddo
  deallocate(s)
end subroutine o3hhuv


subroutine ozone_cw_continuum(jrad, wk, grid, tave, path, absrb)

  integer, intent(in) :: jrad
  real(kind=real64), intent(in) :: wk !< Ozone number density [cm-3].
  real(kind=real64), intent(in) :: tave !< Temperature [K].
  type(SpectralGrid), intent(in) :: grid !< Spectral grid.
  character(len=*), intent(in) :: path !< Path to dataset.
  real(kind=real64), dimension(:), intent(inout) :: absrb !< Extinction [cm-1].

  integer :: npto3, j, ist, last
  real(kind=real64) :: wo3, v1c, v2c, dvc, v1ss, v2ss, dt, vj
  real(kind=real64), dimension(5150):: cch0, cch1, cch2

  if (grid%xn .gt. 8920._real64 .and. grid%x0 .le. 24665._real64) then
    cch0(:) = 0._real64 ![cm3].
    cch1(:) = 0._real64 ![cm3 K-1].
    cch2(:) = 0._real64 ![cm3 K-2].
    wo3 = wk*1.0e-20_real64 ![cm-3].
    call xo3chp(v1c, v2c, dvc, npto3, cch0, cch1, cch2, v1ss, v2ss, grid%x0, grid%xn, path)
    dt = tave - 273.15_real64 ![K].
    do j = 1, npto3
      cch0(j) = (cch0(j) + (cch1(j) + cch2(j)*dt)*dt)*wo3 !unitless.
      vj = v1c + dvc*real(j - 1)
      if (jrad .eq. 1) cch0(j) = cch0(j)*radfn(vj, tave/radcn2) ![cm-1].
    enddo
    call pre_xint(v1ss, v2ss, grid%x0, grid%dx, grid%n, ist, last)
    call xint(v1c, v2c, dvc, cch0, 1._real64, grid%x0, grid%dx, absrb, ist, last)
  endif
end subroutine ozone_cw_continuum


subroutine ozone_hh_continuum(n_absrb, jrad, wk, grid, tave, path, absrb)

  integer, intent(in) :: n_absrb, jrad
  real(kind=real64), intent(in) :: wk !< Ozone number density [cm-3].
  real(kind=real64), intent(in) :: tave !< Temperature [K].
  type(SpectralGrid), intent(in) :: grid !< Spectral grid.
  character(len=*), intent(in) :: path !< Path to dataset.
  real(kind=real64), dimension(:), intent(inout) :: absrb !< Extinction [cm-1].

  integer :: npto3, npt1, npt2, j, i, i_fix, ist, last
  real(kind=real64) :: wo3, tc, v1c, v2c, dvc, v1ss, v2ss, v1t1, v2t1, dvt1, v1t2, &
                       v2t2, dvt2, vj
  real(kind=real64), dimension(6000) :: c
  real(kind=real64), dimension(n_absrb) :: c0, ct1, ct2, absbsv

  if (grid%xn .gt. 27370._real64 .and. grid%x0 .lt. 40800._real64) then
    c0(:) = 0._real64 ![cm3].
    ct1(:) = 0._real64 ![K-1].
    ct2(:) = 0._real64 ![K-2].
    wo3 = wk*1.e-20_real64 ![cm-3].
    tc = tave - 273.15_real64 ![K].
    call o3hht0(v1c, v2c, dvc, npto3, c0, v1ss, v2ss, grid%x0, grid%xn, path)
    call o3hht1(v1t1, v2t1, dvt1, npt1, ct1, grid%x0, grid%xn, path)
    call o3hht2(v1t2, v2t2, dvt2, npt2, ct2, grid%x0, grid%xn, path)
    do j = 1, npto3
      c(j) = c0(j)*wo3 !unitless.
      vj = v1c + dvc*real(j - 1) ![cm-1].
      if (jrad .eq. 1) c(j) = c(j)*radfn(vj, tave/radcn2) ![cm-1].
      c(j) = c(j)*(1. + ct1(j)*tc + ct2(j)*tc*tc) ![cm-1].
    enddo

    !Save non-Hartley Huggins optical depth contribution to
    !prevent double counting for wavenumber region beyond
    !40800 cm-1.
    if (vj .gt. 40815._real64 .and. grid%xn .gt. 40800._real64) then
      i_fix = int((40800._real64 - grid%x0)/grid%dx + 1.001_real64)
      do i = i_fix, grid%n
        absbsv(i) = absrb(i) ![cm-1].
      enddo
    endif

    !Combine Hartley Huggins with previous optical depths.
    call pre_xint(v1ss, v2ss, grid%x0, grid%dx, grid%n, ist, last)
    call xint(v1c, v2c, dvc, c, 1._real64, grid%x0, grid%dx, absrb, ist, last)

    !If V2 > 40800 cm-1, replace points with previously
    !saved values (non-Hartley Huggins contribution)
    if (vj .gt. 40815._real64 .and. grid%xn .gt. 40800._real64) then
      do i = i_fix, grid%n
        absrb(i) = absbsv(i) ![cm-1].
      enddo
    endif
  endif
end subroutine ozone_hh_continuum


subroutine ozone_uv_continuum(n_absrb, jrad, wk, grid, tave, path, absrb)

  integer :: n_absrb, jrad
  real(kind=real64), intent(in) :: wk !< Ozone number density [cm-3].
  real(kind=real64), intent(in) :: tave !< Temperature [K].
  type(SpectralGrid), intent(in) :: grid !< Spectral grid.
  character(len=*), intent(in) :: path !< Path to dataset.
  real(kind=real64), dimension(:), intent(inout) :: absrb !< Extinction [cm-1].

  integer :: npto3, j, i, i_fix, ist, last
  real(kind=real64) :: wo3, v1c, v2c, dvc, v1ss, v2ss, vj
  real(kind=real64), dimension(6000) :: c
  real(kind=real64), dimension(n_absrb) :: c0, absbsv

  if (grid%xn .gt. 40800._real64 .and. grid%x0 .lt. 54000._real64) then
    c0(:) = 0._real64 ![cm3].
    wo3 = wk ![cm-3].
    call o3hhuv(v1c, v2c, dvc, npto3, c0, v1ss, v2ss, grid%x0, grid%xn, path)
    do j = 1, npto3
      c(j) = c0(j)*wo3 !unitless.
      vj = v1c + dvc*real(j - 1) ![cm-1].
      if (jrad .eq. 1) c(j) = c(j)*radfn(vj, tave/radcn2) ![cm-1].
    enddo

    !Save non-Hartley Huggins UV optical depth contribution to
    !prevent double counting for wavenumber region before
    !40800 cm-1.
    if (grid%x0 .lt. 40800._real64) then
      i_fix = int((40800._real64 - grid%x0)/grid%dx + 1.001_real64)
      do i = 1, i_fix - 1
        absbsv(i) = absrb(i) ![cm-1].
      enddo
    endif

    !Combine UV Hartley Huggins with previous optical depths
    call pre_xint(v1ss, v2ss, grid%x0, grid%dx, grid%n, ist, last)
    call xint(v1c, v2c, dvc, c, 1._real64, grid%x0, grid%dx, absrb, ist, last)

    !If V1 < 40800 cm-1, replace points with previously
    !saved values (non-Hartley Huggins UV contribution)
    if (grid%x0 .lt. 40800._real64) then
      do i = 1, i_fix - 1
        absrb(i) = absbsv(i) ![cm-1].
      enddo
    endif
  endif
end subroutine ozone_uv_continuum


end module o3_continuum
