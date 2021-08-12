module h2o_continuum
use, intrinsic :: iso_fortran_env, only: real64
use netcdf
use utils
implicit none


public :: water_vapor_self_continuum, water_vapor_foreign_continuum


contains


subroutine sl296(v1c, v2c, dvc, nptc, c, v1ss, v2ss, v1abs, v2abs, path)

  real(kind=real64), intent(out) :: v1c, v2c, dvc
  integer, intent(out) :: nptc
  real(kind=real64), dimension(:), intent(inout) :: c ![cm3].
  real(kind=real64), intent(out) :: v1ss, v2ss
  real(kind=real64), intent(in) :: v1abs, v2abs
  character(len=*), intent(in) :: path

  integer :: dimid, i, i1, i2, j, ncid, npts, varid
  real(kind=real64) :: dvs, v1s, v2s
  real(kind=real64), dimension(:), allocatable :: s

  !Read data from netcdf file.
  call nc_check(nf90_open(path, nf90_nowrite, ncid))
  call nc_check(nf90_inq_varid(ncid, "bs296", varid))
  call nc_check(nf90_get_att(ncid, varid, "wavenumber_lower_bound", v1s))
  call nc_check(nf90_get_att(ncid, varid, "wavenumber_upper_bound", v2s))
  call nc_check(nf90_get_att(ncid, varid, "wavenumber_resolution", dvs))
  call nc_check(nf90_inq_dimid(ncid, "n2", dimid))
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
    c(j) = s(i) ![cm3].
  enddo
  deallocate(s)
end subroutine sl296


subroutine sl260(v1c, v2c, dvc, nptc, c, v1ss, v2ss, v1abs, v2abs, path)

  real(kind=real64), intent(out) :: v1c, v2c, dvc
  integer, intent(out) :: nptc
  real(kind=real64), dimension(:), intent(inout) :: c ![cm3].
  real(kind=real64), intent(out) :: v1ss, v2ss
  real(kind=real64), intent(in) :: v1abs, v2abs
  character(len=*), intent(in) :: path

  integer :: dimid, i, i1, i2, j, ncid, npts, varid
  real(kind=real64) :: dvs, v1s, v2s
  real(kind=real64), dimension(:), allocatable :: s

  !Read data from netcdf file.
  call nc_check(nf90_open(path, nf90_nowrite, ncid))
  call nc_check(nf90_inq_varid(ncid, "bs260", varid))
  call nc_check(nf90_get_att(ncid, varid, "wavenumber_lower_bound", v1s))
  call nc_check(nf90_get_att(ncid, varid, "wavenumber_upper_bound", v2s))
  call nc_check(nf90_get_att(ncid, varid, "wavenumber_resolution", dvs))
  call nc_check(nf90_inq_dimid(ncid, "n2", dimid))
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
    c(j) = s(i) ![cm3].
  enddo
  deallocate(s)
end subroutine sl260


subroutine frn296(v1c, v2c, dvc, nptc, c, v1ss, v2ss, v1abs, v2abs, path)

  real(kind=real64), intent(out) :: v1c, v2c, dvc
  integer, intent(out) :: nptc
  real(kind=real64), dimension(:), intent(inout) :: c ![cm3].
  real(kind=real64), intent(out) :: v1ss, v2ss
  real(kind=real64), intent(in) :: v1abs, v2abs
  character(len=*), intent(in) :: path

  integer :: dimid, i, i1, i2, j, ncid, npts, varid
  real(kind=real64) :: dvs, v1s, v2s
  real(kind=real64), dimension(:), allocatable :: s

  !Read data from netcdf file.
  call nc_check(nf90_open(path, nf90_nowrite, ncid))
  call nc_check(nf90_inq_varid(ncid, "bfh2o", varid))
  call nc_check(nf90_get_att(ncid, varid, "wavenumber_lower_bound", v1s))
  call nc_check(nf90_get_att(ncid, varid, "wavenumber_upper_bound", v2s))
  call nc_check(nf90_get_att(ncid, varid, "wavenumber_resolution", dvs))
  call nc_check(nf90_inq_dimid(ncid, "n2", dimid))
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
    if ((i .ge. 1) .and. (i .le. npts)) then
      c(j) = s(i) ![cm3].
    endif
  enddo
  deallocate(s)
end subroutine frn296


!> @brief Water vapor self continuum.
subroutine water_vapor_self_continuum(n_absrb, jrad, wk, wtot, grid, &
                                      tave, pave, path, h2o_grid, csh2o, &
                                      absrb, cself)

  integer, intent(in) :: n_absrb, jrad
  real(kind=real64), intent(in) :: wk !< Water vapor number density [cm-3].
  real(kind=real64), intent(in) :: wtot !< Total air density [cm-3].
  type(SpectralGrid), intent(in) :: grid !< Spectral grid.
  real(kind=real64), intent(in) :: tave !< Temperature [K].
  real(kind=real64), intent(in) :: pave !< Pressure [mb].
  character(len=*), intent(in) :: path !< Path to dataset.
  type(SpectralGrid), intent(out) :: h2o_grid
  real(kind=real64), dimension(:), intent(inout) :: csh2o !< [cm3].
  real(kind=real64), dimension(:), intent(inout) :: absrb !< Extinction [cm-1].
  real(kind=real64), dimension(:), intent(inout) :: cself !< Extinction [cm-1].

  integer :: j, nptc, ist, last
  real(kind=real64) :: h2o_fac, rself, v1c, v2c, dvc, v1ss, v2ss, vj, sh2o, tfac
  real(kind=real64), dimension(n_absrb) :: sh2ot0, sh2ot1

  h2o_fac = wk/wtot !unitless.
  rself = h2o_fac*(pave/p0)*(t0/tave)*1.e-20_real64 !unitless.

  if (grid%xn .gt. -20.0_real64 .and. grid%x0 .lt. 20000._real64) then
    sh2ot0(:) = 0._real64 ![cm3].
    sh2ot1(:) = 0._real64 ![cm3].
    call sl296(v1c, v2c, dvc, nptc, sh2ot0, v1ss, v2ss, grid%x0, grid%xn, path)
    call sl260(v1c, v2c, dvc, nptc, sh2ot1, v1ss, v2ss, grid%x0, grid%xn, path)
    tfac = (tave - t0)/(260._real64 - t0) !unitless
    do j = 1, nptc
      vj = v1c + dvc*real(j - 1) ![cm-1].
      sh2o = 0.
      if (sh2ot0(j) .gt. 0._real64) then
        sh2o = sh2ot0(j)*(sh2ot1(j)/sh2ot0(j))**tfac ![cm3].
      endif
      cself(j) = wk*sh2o*rself !unitless
      h2o_grid%x0 = v1c
      h2o_grid%xn = v2c
      h2o_grid%dx = dvc
      h2o_grid%n = nptc
      csh2o(j) = 1.e-20_real64*sh2o ![cm3].
      if (jrad .eq. 1) cself(j) = cself(j)*radfn(vj, tave/radcn2) ![cm-1].
    enddo

    !Interpolate to total optical depth grid.
    call pre_xint(v1ss, v2ss, grid%x0, grid%dx, grid%n, ist, last)
    call xint(v1c, v2c, dvc, cself, 1._real64, grid%x0, grid%dx, absrb, ist, last)
  endif
end subroutine water_vapor_self_continuum


!> @brief Water vapor foreign continuum.
subroutine water_vapor_foreign_continuum(n_absrb, icflg, jrad, wk, wtot, grid, &
                                         tave, pave, cself, path, cfh2o, absrb)

  integer, intent(in) :: n_absrb, icflg, jrad
  real(kind=real64), intent(in) :: wk !< Water vapor number density [cm-3].
  real(kind=real64), intent(in) :: wtot !< Total air density [cm-3].
  type(SpectralGrid), intent(in) :: grid !< Spectral grid.
  real(kind=real64), intent(in) :: tave !< Temperature [K].
  real(kind=real64), intent(in) :: pave !< Pressure [mb].
  real(kind=real64), dimension(:), intent(in) :: cself !< Extinction [cm-1].
  character(len=*), intent(in) :: path !< Path to dataset.
  real(kind=real64), dimension(:), intent(inout) :: cfh2o !< [cm3].
  real(kind=real64), dimension(:), intent(inout) :: absrb !< Extinction [cm-1].

  integer :: ncid, varid, n_1, n_2, nptc, j, jfac, ist, last
  integer, parameter :: ipts = 5050
  real(kind=real64) :: f0, v0f1, hwsq1, beta1, c_1, c_2, beta2, v1c, v2c, dvc, v1ss, &
                       v2ss, vj, fscal, vdelsq1, vdelmsq1, vf1, vmf1, vf2, h2o_fac, &
                       rfrgn, rfrgn_aj, c_f
  real(kind=real64), dimension(n_absrb) :: fh2o
  real(kind=real64), dimension(-1:61) :: xfac_rhu
  real(kind=real64), dimension(6000) :: c
  real(kind=real64), dimension(ipts) :: cfrgn_aj

  if (grid%xn .gt. -20.0_real64 .and. grid%x0 .lt. 20000._real64) then
    call nc_check(nf90_open(path, nf90_nowrite, ncid))
    call nc_check(nf90_inq_varid(ncid, "xfac_rhu", varid))
    call nc_check(nf90_get_var(ncid, varid, xfac_rhu))
    call nc_check(nf90_close(ncid))
    h2o_fac = wk/wtot !unitless.
    rfrgn = (1._real64 - h2o_fac)*(pave/p0)*(t0/tave)*1.e-20_real64 !unitless
    rfrgn_aj = h2o_fac*(pave/p0)*(t0/tave)*1.e-20_real64 !unitless
    fh2o(:) = 0.
    f0 = 0.06
    v0f1 = 255.67
    hwsq1 = 240.**2
    beta1 = 57.83
    c_1 = -0.42
    n_1 = 8
    c_2 = 0.3
    beta2 = 630.
    n_2 = 8
    call frn296(v1c, v2c, dvc, nptc, fh2o, v1ss, v2ss, grid%x0, grid%xn, path)
    do j = 1, nptc
      vj = v1c + dvc*real(j - 1) ![cm-1].
      if (vj .le. 600._real64) then
        jfac = int((vj + 10._real64)/10._real64 + 0.00001_real64)
        fscal = xfac_rhu(jfac)
      else
        vdelsq1 = (vj - v0f1)**2
        vdelmsq1 = (vj + v0f1)**2
        vf1 = ((vj - v0f1)/beta1)**n_1
        vmf1 = ((vj + v0f1)/beta1)**n_1
        vf2 = (vj/beta2)**n_2
        fscal = 1. + &
                (f0 + c_1*((hwsq1/(vdelsq1 + hwsq1 + vf1)) + &
                (hwsq1/(vdelmsq1 + hwsq1 + vmf1))))/ &
                (1. + c_2*vf2)
      endif
      fh2o(j) = fh2o(j)*fscal ![cm3].
      c_f = wk*fh2o(j) !unitless
      cfh2o(j) = 1.e-20_real64*fh2o(j) ![cm3].
      if (jrad .eq. 1) c_f = c_f*radfn(vj, tave/radcn2) ![cm-1].
      c(j) = c_f*rfrgn ![cm-1].
      cfrgn_aj(j) = c_f*rfrgn_aj ![cm-1].
    enddo
    call pre_xint(v1ss, v2ss, grid%x0, grid%dx, grid%n, ist, last)
    call xint(v1c, v2c, dvc, c, 1._real64, grid%x0, grid%dx, absrb, ist, last)
    if (icflg .eq. 1) then
      do j = 1, nptc
        c(j) = cself(j) - cfrgn_aj(j) ![cm-1].
      enddo
      call pre_xint(v1ss, v2ss, grid%x0, grid%dx, grid%n, ist, last)
      call xint(v1c, v2c, dvc, c, 1._real64, grid%x0, grid%dx, absrb, ist, last)
    endif
  endif
end subroutine water_vapor_foreign_continuum


end module h2o_continuum
