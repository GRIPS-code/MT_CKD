module n2_continuum
use, intrinsic :: iso_fortran_env, only: real64
use netcdf
use utils
implicit none


public :: nitrogen_continuum, nitrogen_fundamental_continuum, nitrogen_overtone_continuum


contains


subroutine xn2_r(v1c, v2c, dvc, nptc, c, fo2, tave, v1ss, v2ss, v1abs, v2abs, path)

  real(kind=real64), intent(out) :: v1c, v2c, dvc
  integer, intent(out) :: nptc
  real(kind=real64), dimension(:), intent(inout) :: c, fo2
  real(kind=real64), intent(in) :: tave
  real(kind=real64), intent(out) :: v1ss, v2ss
  real(kind=real64), intent(in) :: v1abs, v2abs
  character(len=*), intent(in) :: path

  integer :: dimid, i, i1, i2, j, ncid, nptb, npts, varid
  real(kind=real64) :: dvb, dvs, sf_t, tfac, v1b, v2b, v1s, v2s, xn2, xo2
  real(kind=real64), dimension(:), allocatable :: c_220, c_296, sf_220, sf_296

  !Read data from netcdf file.
  call nc_check(nf90_open(path, nf90_nowrite, ncid))
  call nc_check(nf90_inq_varid(ncid, "ct_296", varid))
  call nc_check(nf90_get_att(ncid, varid, "wavenumber_lower_bound", v1s))
  call nc_check(nf90_get_att(ncid, varid, "wavenumber_upper_bound", v2s))
  call nc_check(nf90_get_att(ncid, varid, "wavenumber_resolution", dvs))
  call nc_check(nf90_inq_dimid(ncid, "n4", dimid))
  call nc_check(nf90_inquire_dimension(ncid, dimid, len=npts))
  allocate(c_296(npts))
  call nc_check(nf90_get_var(ncid, varid, c_296))
  call nc_check(nf90_inq_varid(ncid, "sf_296", varid))
  allocate(sf_296(npts))
  call nc_check(nf90_get_var(ncid, varid, sf_296))
  call nc_check(nf90_inq_varid(ncid, "ct_220", varid))
  call nc_check(nf90_get_att(ncid, varid, "wavenumber_lower_bound", v1b))
  call nc_check(nf90_get_att(ncid, varid, "wavenumber_upper_bound", v2b))
  call nc_check(nf90_get_att(ncid, varid, "wavenumber_resolution", dvb))
  call nc_check(nf90_inq_dimid(ncid, "n4", dimid))
  call nc_check(nf90_inquire_dimension(ncid, dimid, len=nptb))
  allocate(c_220(npts))
  call nc_check(nf90_get_var(ncid, varid, c_220))
  call nc_check(nf90_inq_varid(ncid, "sf_220", varid))
  allocate(sf_220(npts))
  call nc_check(nf90_get_var(ncid, varid, sf_220))
  call nc_check(nf90_close(ncid))

  xo2 = 0.21_real64
  xn2 = 0.79_real64

  tfac = (tave - 296._real64)/(220._real64 - 296._real64) !unitless.

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
    if ((i .lt. 1) .or. (i .gt. npts)) continue
    c(j) = c_296(i)*((c_220(i)/c_296(i))**tfac) !unitless.
    sf_t = sf_296(i)*((sf_220(i)/sf_296(i))**tfac) !unitless.
    fo2(j) = (sf_t - 1._real64)*xn2/xo2 !unitless.
  enddo
  deallocate(c_220, c_296, sf_220, sf_296)
end subroutine xn2_r


subroutine n2_ver_1(v1c, v2c, dvc, nptc, c, c1, c2, t, v1ss, v2ss, v1abs, v2abs, path)

  real(kind=real64), intent(out) :: v1c, v2c, dvc
  integer, intent(out) :: nptc
  real(kind=real64), dimension(:), intent(inout) :: c, c1, c2
  real(kind=real64), intent(in) :: t
  real(kind=real64), intent(out) :: v1ss, v2ss
  real(kind=real64), intent(in) :: v1abs, v2abs
  character(len=*), intent(in) :: path

  integer :: dimid, i, i1, i2, j, ncid, npts, varid
  real(kind=real64) :: a_o2, dvs, vj, v1s, v2s, xt_lin, xtfac
  real(kind=real64), dimension(:), allocatable :: xn2_272, xn2_228, a_h2o

  !Read data from netcdf file.
  call nc_check(nf90_open(path, nf90_nowrite, ncid))
  call nc_check(nf90_inq_varid(ncid, "xn2_272", varid))
  call nc_check(nf90_get_att(ncid, varid, "wavenumber_lower_bound", v1s))
  call nc_check(nf90_get_att(ncid, varid, "wavenumber_upper_bound", v2s))
  call nc_check(nf90_get_att(ncid, varid, "wavenumber_resolution", dvs))
  call nc_check(nf90_inq_dimid(ncid, "n5", dimid))
  call nc_check(nf90_inquire_dimension(ncid, dimid, len=npts))
  allocate(xn2_272(npts), xn2_228(npts), a_h2o(npts))
  call nc_check(nf90_get_var(ncid, varid, xn2_272))
  call nc_check(nf90_inq_varid(ncid, "xn2_228", varid))
  call nc_check(nf90_get_var(ncid, varid, xn2_228))
  call nc_check(nf90_inq_varid(ncid, "a_h2o", varid))
  call nc_check(nf90_get_var(ncid, varid, a_h2o))
  call nc_check(nf90_close(ncid))

  xtfac = ((1._real64/t) - (1._real64/272._real64))/((1._real64/228._real64) - &
          (1._real64/272._real64)) !unitless.]
  xt_lin = (t - 272._real64)/(228._real64 - 272._real64) !unitless.
  a_o2  = 1.294_real64 - 0.4545_real64*t/296._real64 !unitless.

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
    if ((i .lt. 1) .or. (i .gt. npts)) continue
    vj = v1c + dvc*real(j - 1)
    if ((xn2_272(i) .gt. 0._real64) .and. (xn2_228(i) .gt. 0._real64)) then
       c(j) = xn2_272(i)*(xn2_228(i)/xn2_272(i))**xtfac
    else
      c(j) = xn2_272(i) + (xn2_228(i) - xn2_272(i))*xt_lin
    endif
    c(j) = c(j)/vj !unitless
    c1(j) = a_o2*c(j) !unitless
    c2(j) = (9._real64/7._real64)*a_h2o(i)*c(j) !unitless.
  enddo
  deallocate(xn2_272, xn2_228, a_h2o)
end subroutine n2_ver_1


subroutine n2_overtone1(v1c, v2c, dvc, nptc, c, v1ss, v2ss, v1abs, v2abs, path)

  real(kind=real64), intent(out) :: v1c, v2c, dvc ![cm-1].
  integer, intent(out) :: nptc
  real(kind=real64), dimension(:), intent(inout) :: c ![1.e20 cm3 amagat-1].
  real(kind=real64), intent(out) :: v1ss, v2ss
  real(kind=real64), intent(in) :: v1abs, v2abs
  character(len=*), intent(in) :: path

  integer :: dimid, i, i1, i2, j, ncid, npts, varid
  real(kind=real64) :: dvs, vj, v1s, v2s
  real(kind=real64), dimension(:), allocatable :: xn2

  !Read data from netcdf file.
  call nc_check(nf90_open(path, nf90_nowrite, ncid))
  call nc_check(nf90_inq_varid(ncid, "xn2", varid))
  call nc_check(nf90_get_att(ncid, varid, "wavenumber_lower_bound", v1s))
  call nc_check(nf90_get_att(ncid, varid, "wavenumber_upper_bound", v2s))
  call nc_check(nf90_get_att(ncid, varid, "wavenumber_resolution", dvs))
  call nc_check(nf90_inq_dimid(ncid, "n6", dimid))
  call nc_check(nf90_inquire_dimension(ncid, dimid, len=npts))
  allocate(xn2(npts))
  call nc_check(nf90_get_var(ncid, varid, xn2))
  call nc_check(nf90_close(ncid))

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
    c(j) = 0.
    if ((i .lt. 1) .or. (i .gt. npts)) continue
    vj = v1c + dvc*real(j - 1)
    c(j) = xn2(i) ![cm-1].
    c(j) = c(j)/vj !unitless
  enddo
  deallocate(xn2)
end subroutine n2_overtone1


subroutine nitrogen_continuum(n_absrb, jrad, wn2, grid, &
                              tave, pave, x_vmr_n2, x_vmr_o2, x_vmr_h2o, path, absrb)

  integer, intent(in) :: n_absrb, jrad
  real(kind=real64), intent(in) :: wn2 !< Nitrogen number density [cm-3].
  real(kind=real64), intent(in) :: tave !< Temperature [K].
  real(kind=real64), intent(in) :: pave !< Pressure [mb].
  real(kind=real64), intent(in) :: x_vmr_n2, x_vmr_o2, x_vmr_h2o
  type(SpectralGrid), intent(in) :: grid !< Spectral grid.
  character(len=*), intent(in) :: path !< Path to dataset.
  real(kind=real64), dimension(:), intent(inout) :: absrb !< Extinction [cm-1].

  integer :: nptc, j, ist, last
  real(kind=real64) :: a_h2o, tau_fac, v1c, v2c, dvc, v1ss, v2ss, vj
  real(kind=real64), dimension(n_absrb) :: c0, c1
  real(kind=real64), dimension(6000) :: c

  if (grid%xn .gt. -10._real64 .and. grid%x0 .lt. 350._real64) then
    c0(:) = 0._real64 !unitless.
    c1(:) = 0._real64 !unitless.
    a_h2o = 1._real64 !unitless.
    tau_fac = (wn2/loschmidt)*(pave/p0)*(t_273/tave) !unitless.
    call xn2_r(v1c, v2c, dvc, nptc, c0, c1, tave, v1ss, v2ss, grid%x0, grid%xn, path)
    do j = 1, nptc
      vj = v1c + dvc*real(j - 1) ![cm-1].
      c(j) = tau_fac*c0(j)*(x_vmr_n2 + c1(j)*x_vmr_o2 + a_h2o*x_vmr_h2o) !unitless.
      if (jrad .eq. 1) c(j) = c(j)*radfn(vj, tave/radcn2) ![cm-1].
    enddo
    call pre_xint(v1ss, v2ss, grid%x0, grid%dx, grid%n, ist, last)
    call xint(v1c, v2c, dvc, c, 1._real64, grid%x0, grid%dx, absrb, ist, last)
  endif
end subroutine nitrogen_continuum


subroutine nitrogen_fundamental_continuum(jrad, wn2, grid, &
                                          tave, pave, x_vmr_n2, x_vmr_o2, x_vmr_h2o, &
                                          path, absrb)

  integer, intent(in) :: jrad
  real(kind=real64), intent(in) :: wn2 !< Nitrogen number density [cm-3].
  real(kind=real64), intent(in) :: tave !< Temperature [K].
  real(kind=real64), intent(in) :: pave !< Pressure [mb].
  real(kind=real64), intent(in) :: x_vmr_n2, x_vmr_o2, x_vmr_h2o
  type(SpectralGrid), intent(in) :: grid !< Spectral grid.
  character(len=*), intent(in) :: path !< Path to dataset.
  real(kind=real64), dimension(:), intent(inout) :: absrb !< Extinction [cm-1].

  integer :: nptc, j, ist, last
  real(kind=real64) :: tau_fac, v1c, v2c, dvc, v1ss, v2ss, vj
  real(kind=real64), dimension(5150) :: cn0, cn1, cn2 ![cm3 amg-1].
  real(kind=real64), dimension(6000) :: c

  if (grid%xn .gt. 2001.77_real64 .and. grid%x0 .lt. 2897.59_real64) then
    cn0(:) = 0._real64 !unitless.
    cn1(:) = 0._real64 !unitless.
    cn2(:) = 0._real64 !unitless.
    tau_fac = (wn2/loschmidt)*(pave/p0)*(t_273/tave) !unitless.
    call n2_ver_1(v1c, v2c, dvc, nptc, cn0, cn1, cn2, tave, v1ss, v2ss, grid%x0, grid%xn, path)
    do j = 1, nptc
      vj = v1c + dvc*real(j - 1) ![cm-1].
      c(j) = tau_fac*(x_vmr_n2*cn0(j) + x_vmr_o2*cn1(j) + x_vmr_h2o*cn2(j)) !unitless.
      if (jrad .eq. 1) c(j) = c(j)*radfn(vj, tave/radcn2) ![cm-1].
    enddo
    call pre_xint(v1ss, v2ss, grid%x0, grid%dx, grid%n, ist, last)
    call xint(v1c, v2c, dvc, c, 1._real64, grid%x0, grid%dx, absrb, ist, last)
  endif
end subroutine nitrogen_fundamental_continuum


subroutine nitrogen_overtone_continuum(n_absrb, jrad, wn2, grid, &
                                       tave, pave, x_vmr_n2, x_vmr_o2, x_vmr_h2o, &
                                       path, absrb)

  integer, intent(in):: n_absrb, jrad
  real(kind=real64), intent(in) :: wn2 !< Nitrogen number density [cm-3].
  type(SpectralGrid), intent(in) :: grid !< Spectral grid.
  real(kind=real64), intent(in) :: tave !< Temperature [K].
  real(kind=real64), intent(in) :: pave !< Pressure [mb].
  real(kind=real64), intent(in) :: x_vmr_n2 !< N2 volume mixing ratio [mol mol-1].
  real(kind=real64), intent(in) :: x_vmr_o2 !< O2 volume mixing ratio [mol mol-1].
  real(kind=real64), intent(in) :: x_vmr_h2o !< H2O volume mixing ratio [mol mol-1].
  character(len=*), intent(in) :: path !< Path to the input dataset.
  real(kind=real64), dimension(:), intent(inout) :: absrb !< Extinction [cm-1].

  integer :: nptc, j, ist, last
  real(kind=real64) :: tau_fac, v1c, v2c, dvc, v1ss, v2ss, vj
  real(kind=real64), dimension(n_absrb) :: c0
  real(kind=real64), dimension(6000) :: c

  if (grid%xn .gt. 4340._real64 .and. grid%x0 .lt. 4910._real64) then
    c0(:) = 0._real64 !unitless.
    tau_fac = (wn2/loschmidt)*(pave/p0)*(t_273/tave)*(x_vmr_n2 + x_vmr_o2 + x_vmr_h2o) !unitless.
    call n2_overtone1(v1c, v2c, dvc, nptc, c0, v1ss, v2ss, grid%x0, grid%xn, path)
    do j = 1, nptc
      vj = v1c + dvc*real(j - 1) !unitless.
      c(j) = tau_fac*c0(j) !unitless.
      if (jrad .eq. 1) c(j) = c(j)*radfn(vj, tave/radcn2) ![cm-1].
    enddo
    call pre_xint(v1ss, v2ss, grid%x0, grid%dx, grid%n, ist, last)
    call xint(v1c, v2c, dvc, c, 1._real64, grid%x0, grid%dx, absrb, ist, last)
  endif
end subroutine nitrogen_overtone_continuum


end module n2_continuum
