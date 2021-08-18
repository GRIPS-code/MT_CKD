module o2_continuum
use, intrinsic :: iso_fortran_env, only: real64
use netcdf
use utils
implicit none


public :: oxygen_fundamental_continuum, oxygen_nir1_continuum, oxygen_nir2_continuum, &
          oxygen_nir3_continuum, oxygen_visible_continuum, oxygen_herzberg_continuum, &
          oxygen_uv_continuum


contains


subroutine o2_ver_1(v1c, v2c, dvc, nptc, c, t, v1ss, v2ss, v1abs, v2abs, path)

  real(kind=real64), intent(out) :: v1c, v2c, dvc
  integer, intent(out) :: nptc
  real(kind=real64), dimension(:), allocatable, intent(inout) :: c ![cm3].
  real(kind=real64), intent(in) :: t
  real(kind=real64), intent(out) :: v1ss, v2ss
  real(kind=real64), intent(in) :: v1abs, v2abs
  character(len=*), intent(in) :: path

  integer :: dimid, i, i1, i2, j, ncid, npts, varid
  real(kind=real64) :: dvs, factor, vj, v1s, v2s, xktfac
  real(kind=real64), dimension(:), allocatable :: xo2, xo2t

  !Read data from netcdf file.
  call nc_check(nf90_open(path, nf90_nowrite, ncid))
  call nc_check(nf90_inq_varid(ncid, "o2_f", varid))
  call nc_check(nf90_get_att(ncid, varid, "wavenumber_lower_bound", v1s))
  call nc_check(nf90_get_att(ncid, varid, "wavenumber_upper_bound", v2s))
  call nc_check(nf90_get_att(ncid, varid, "wavenumber_resolution", dvs))
  call nc_check(nf90_inq_dimid(ncid, "n10", dimid))
  call nc_check(nf90_inquire_dimension(ncid, dimid, len=npts))
  allocate(xo2(npts), xo2t(npts))
  call nc_check(nf90_get_var(ncid, varid, xo2))
  call nc_check(nf90_inq_varid(ncid, "o2_t", varid))
  call nc_check(nf90_get_var(ncid, varid, xo2t))
  call nc_check(nf90_close(ncid))

  xktfac = (1._real64/296._real64) - (1._real64/t) ![K-1].
  factor = (1.e20_real64/loschmidt) ![cm3].

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

  allocate(c(nptc))
  do j = 1, nptc
    i = i1 + (j - 1)
    c(j) = 0.
    if ((i .lt. 1) .or. (i .gt. npts)) continue
    vj = v1c + dvc*real(j - 1)
    c(j) = factor*xo2(i)*exp(xo2t(i)*xktfac)/vj ![cm3].
  enddo
  deallocate(xo2, xo2t)
end subroutine o2_ver_1


subroutine o2inf1(v1c, v2c, dvc, nptc, c, v1ss, v2ss, v1abs, v2abs, path)

  real(kind=real64), intent(out) :: v1c, v2c, dvc
  integer, intent(out) :: nptc
  real(kind=real64), dimension(:), allocatable, intent(inout) :: c !unitless.
  real(kind=real64), intent(out) :: v1ss, v2ss
  real(kind=real64), intent(in) :: v1abs, v2abs
  character(len=*), intent(in) :: path

  integer :: dimid, i, i1, i2, j, ncid, npts, varid
  real(kind=real64) :: dvs, vj, v1s, v2s
  real(kind=real64), dimension(:), allocatable :: xo2inf1

  !Read data from netcdf file.
  call nc_check(nf90_open(path, nf90_nowrite, ncid))
  call nc_check(nf90_inq_varid(ncid, "o2_inf1", varid))
  call nc_check(nf90_get_att(ncid, varid, "wavenumber_lower_bound", v1s))
  call nc_check(nf90_get_att(ncid, varid, "wavenumber_upper_bound", v2s))
  call nc_check(nf90_get_att(ncid, varid, "wavenumber_resolution", dvs))
  call nc_check(nf90_inq_dimid(ncid, "n11", dimid))
  call nc_check(nf90_inquire_dimension(ncid, dimid, len=npts))
  allocate(xo2inf1(npts))
  call nc_check(nf90_get_var(ncid, varid, xo2inf1))
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

  allocate(c(nptc))
  do j = 1, nptc
    i = i1 + (j - 1)
    c(j) = 0._real64
    if (i .lt. 1 .or. i .gt. npts) continue
    vj = v1c + dvc*real(j - 1)
    c(j) = xo2inf1(i)/vj !unitless.
  enddo
  deallocate(xo2inf1)
end subroutine o2inf1


subroutine o2inf2(v1c, v2c, dvc, nptc, c, v1ss, v2ss, v1abs, v2abs)

  real(kind=real64), intent(out) :: v1c, v2c, dvc
  integer, intent(out) :: nptc
  real(kind=real64), dimension(:), allocatable, intent(inout) :: c ![cm3].
  real(kind=real64), intent(out) :: v1ss, v2ss
  real(kind=real64), intent(in) :: v1abs, v2abs

  integer :: j
  real(kind=real64) :: damp1, damp2, dvs, dv1, dv2, hw1, hw2, o2inf, s1, s2, v1_osc, &
                       v2_osc, vj, v1s, v2s

  v1_osc = 9375._real64
  hw1 = 58.96_real64
  v2_osc = 9439._real64
  hw2 = 45.04_real64
  s1 = 1.166e-04_real64
  s2 = 3.086e-05_real64

  v1s = 9100._real64
  v2s = 11000._real64
  dvs = 2._real64
  dvc = dvs
  v1ss = v1s
  v2ss = v2s

  v1c = v1abs - dvc
  v2c = v2abs + dvc

  ! the following lines prevent a possible problem that can only occur if
  ! v2abs-v1abs >> 2000 cm-1 (i.e. in the standalone continuum code).
  if (v1c .lt. v1s) v1c = v1s - 2._real64*dvs
  if (v2c .gt. v2s) v2c = v2s + 2._real64*dvs
  nptc = int((v2c - v1c)/dvc + 3.01_real64)
  v2c = v1c + dvc*real(nptc - 1)

  allocate(c(nptc))
  do j = 1, nptc
    c(j) = 0._real64
    vj = v1c + dvc*real(j - 1)
    if ((vj .gt. v1s) .and. (vj .lt. v2s)) then
      dv1 = vj - v1_osc
      dv2 = vj - v2_osc
      if (dv1 .lt. 0._real64) then
        damp1 = exp(dv1/176.1_real64)
      else
        damp1 = 1._real64
      endif
      if (dv2 .lt. 0._real64) then
        damp2 = exp(dv2/176.1_real64)
      else
        damp2 = 1._real64
      endif
      o2inf = 0.31831_real64*(((s1*damp1/hw1)/(1._real64 + (dv1/hw1)**2)) &
              + ((s2*damp2/hw2)/(1._real64 + (dv2/hw2)**2)))*1.054_real64
      c(j) = o2inf/vj
    endif
  enddo
end subroutine o2inf2


subroutine o2inf3(v1c, v2c, dvc, nptc, c, v1ss, v2ss, v1abs, v2abs, path)

  real(kind=real64), intent(out) :: v1c, v2c, dvc
  integer, intent(out) :: nptc
  real(kind=real64), dimension(:), allocatable, intent(inout) :: c !unitless.
  real(kind=real64), intent(out) :: v1ss, v2ss
  real(kind=real64), intent(in) :: v1abs, v2abs
  character(len=*), intent(in) :: path

  integer :: dimid, i, i1, i2, j, ncid, npts, varid
  real(kind=real64) :: dvs, vj, v1s, v2s
  real(kind=real64), dimension(:), allocatable :: xo2inf3

  !Read data from netcdf file.
  call nc_check(nf90_open(path, nf90_nowrite, ncid))
  call nc_check(nf90_inq_varid(ncid, "o2_inf3", varid))
  call nc_check(nf90_get_att(ncid, varid, "wavenumber_lower_bound", v1s))
  call nc_check(nf90_get_att(ncid, varid, "wavenumber_upper_bound", v2s))
  call nc_check(nf90_get_att(ncid, varid, "wavenumber_resolution", dvs))
  call nc_check(nf90_inq_dimid(ncid, "n12", dimid))
  call nc_check(nf90_inquire_dimension(ncid, dimid, len=npts))
  allocate(xo2inf3(npts))
  call nc_check(nf90_get_var(ncid, varid, xo2inf3))
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

  allocate(c(nptc))
  do j = 1, nptc
    i = i1 + (j - 1)
    c(j) = 0._real64
    if ((i .lt. 1) .or. (i .gt. npts)) continue
    vj = v1c + dvc*real(j - 1)
    c(j) = xo2inf3(i)/vj !unitless.
  enddo
  deallocate(xo2inf3)
end subroutine o2inf3


subroutine o2_vis(v1c, v2c, dvc, nptc, c, v1ss, v2ss, v1abs, v2abs, path)

  real(kind=real64), intent(out) :: v1c, v2c, dvc
  integer, intent(out) :: nptc
  real(kind=real64), dimension(:), allocatable, intent(inout) :: c ![cm3].
  real(kind=real64), intent(out) :: v1ss, v2ss
  real(kind=real64), intent(in) :: v1abs, v2abs
  character(len=*), intent(in) :: path

  integer :: dimid, i, i1, i2, j, ncid, npts, varid
  real(kind=real64) :: factor, vj
  real(kind=real64) :: v1s, v2s, dvs
  real(kind=real64), dimension(:), allocatable :: s

  !Read data from netcdf file.
  call nc_check(nf90_open(path, nf90_nowrite, ncid))
  call nc_check(nf90_inq_varid(ncid, "o2_invis", varid))
  call nc_check(nf90_get_att(ncid, varid, "wavenumber_lower_bound", v1s))
  call nc_check(nf90_get_att(ncid, varid, "wavenumber_upper_bound", v2s))
  call nc_check(nf90_get_att(ncid, varid, "wavenumber_resolution", dvs))
  call nc_check(nf90_inq_dimid(ncid, "n13", dimid))
  call nc_check(nf90_inquire_dimension(ncid, dimid, len=npts))
  allocate(s(npts))
  call nc_check(nf90_get_var(ncid, varid, s))
  call nc_check(nf90_close(ncid))

  factor = 1./((loschmidt*1.e-20*(55.*273./296.)**2)*89.5) ![cm3].

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

  allocate(c(nptc))
  do j = 1, nptc
    i = i1 + (j - 1)
    c(j) = 0.
    if ((i .lt. 1) .or. (i .gt. npts)) continue
    vj = v1c + dvc*real(j - 1)
    c(j) = factor*s(i)/vj ![cm3].
  enddo
  deallocate(s)
end subroutine o2_vis


subroutine o2herz(v1c, v2c, dvc, nptc, c, t, p, v1ss, v2ss, v1abs, v2abs)

  real(kind=real64), intent(out) :: v1c, v2c, dvc
  integer, intent(out) :: nptc
  real(kind=real64), dimension(:), allocatable, intent(inout) :: c ![cm3].
  real(kind=real64), intent(in) :: t, p
  real(kind=real64), intent(out) :: v1ss, v2ss
  real(kind=real64), intent(in) :: v1abs, v2abs

  integer :: i, i1, i2, j
  real(kind=real64) :: dvs, herz, vj, v1s

  v1s = 36000._real64
  v1ss = v1s
  v2ss = 99999._real64

  dvs = 10._real64
  dvc = dvs

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
  v2c = v1c + dvs*real(nptc - 1)

  allocate(c(nptc))
  do j = 1, nptc
    i = i1 + (j - 1)
    c(j) = 0._real64
    if (i .lt. 1) continue
    vj = v1c + dvc*real(j - 1)
    call hertda(herz, vj)
    call herprs(herz, t, p)
    c(j) = herz/vj ![cm3].
  enddo
end subroutine o2herz


subroutine hertda(herz, v)

  real(kind=real64), intent(out) :: herz
  real(kind=real64), intent(in) :: v

  real(kind=real64) :: corr, yratio

  herz = 0._real64
  if (v .le. 36000._real64) return
  corr = 0._real64
  if (v .le. 40000._real64) corr = ((40000._real64 - v)/4000._real64)*7.917e-7_real64
  yratio = v/48811._real64
  herz = 6.884e-4_real64*(yratio)*exp(-69.738_real64*(log(yratio))**2) - corr
end subroutine hertda


subroutine herprs(herz, t, p)

  real(kind=real64), intent(inout) :: herz
  real(kind=real64), intent(in) :: t, p

  herz = herz*(1._real64 + 0.83_real64*(p/p0)*(273.16_real64/t))
end subroutine herprs


!> @brief O2 continuum in the far-uv is from lu et al. (2010).
subroutine o2fuv(v1c, v2c, dvc, nptc, c, v1ss, v2ss, v1abs, v2abs, path)

  real(kind=real64), intent(out) :: v1c, v2c, dvc
  integer, intent(out) :: nptc
  real(kind=real64), dimension(:), allocatable, intent(inout) :: c ![cm3].
  real(kind=real64), intent(out) :: v1ss, v2ss
  real(kind=real64), intent(in) :: v1abs, v2abs
  character(len=*), intent(in) :: path

  integer :: dimid, i, i1, i2, j, ncid, npts, varid
  real(kind=real64) :: dvs, vj, v1s, v2s
  real(kind=real64), dimension(:), allocatable  :: s

  !Read data from netcdf file.
  call nc_check(nf90_open(path, nf90_nowrite, ncid))
  call nc_check(nf90_inq_varid(ncid, "o2_infuv", varid))
  call nc_check(nf90_get_att(ncid, varid, "wavenumber_lower_bound", v1s))
  call nc_check(nf90_get_att(ncid, varid, "wavenumber_upper_bound", v2s))
  call nc_check(nf90_get_att(ncid, varid, "wavenumber_resolution", dvs))
  call nc_check(nf90_inq_dimid(ncid, "n14", dimid))
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
    i1 = int((v1c - v1s)/dvs + 1.e-5_real64)
  endif
  v1c = v1s + dvs*real(i1 - 1)
  i2 = int((v2c - v1s)/dvs + 1.e-5_real64)

  nptc = i2 - i1 + 3
  if (nptc .gt. npts) nptc = npts + 4
  v2c = v1c + dvs*real(nptc - 1)

  allocate(c(nptc))
  do j = 1, nptc
    i = i1 + (j - 1)
    c(j) = 0.
    if ((i .lt. 1) .or. (i .gt. npts)) continue
    vj = v1c + dvc*real(j - 1)
    c(j) = s(i)/vj ![cm3].
  enddo
  deallocate(s)
end subroutine o2fuv


subroutine oxygen_fundamental_continuum(jrad, wk, grid, tave, pave, path, absrb)

  integer, intent(in) :: jrad
  real(kind=real64), intent(in) :: wk !< Oxygen number density [cm-3].
  real(kind=real64), intent(in) :: tave !< Temperature [K].
  real(kind=real64), intent(in) :: pave !< Pressure [mb].
  type(SpectralGrid), intent(in) :: grid !< Spectral grid.
  character(len=*), intent(in) :: path !< Path to dataset.
  real(kind=real64), dimension(:), intent(inout) :: absrb !< Extinction [cm-1].

  integer :: nptc, j, ist, last
  real(kind=real64) :: tau_fac, v1c, v2c, dvc, v1ss, v2ss, vj
  real(kind=real64), dimension(:), allocatable :: c, c0

  if (grid%xn .gt. 1340._real64 .and. grid%x0 .lt. 1850._real64) then
    tau_fac = wk*1.e-20_real64*(pave/p0)*(t_273/tave) ![cm-3].
    call o2_ver_1(v1c, v2c, dvc, nptc, c0, tave, v1ss, v2ss, grid%x0, grid%xn, path)
    allocate(c(nptc))
    do j = 1, nptc
      vj = v1c + dvc*real(j - 1) ![cm-1].
      c(j) = tau_fac*c0(j) !unitless.
      if (jrad .eq. 1) c(j) = c(j)*radfn(vj, tave/radcn2) ![cm-1].
    enddo
    call pre_xint(v1ss, v2ss, grid%x0, grid%dx, grid%n, ist, last)
    call xint(v1c, v2c, dvc, c, 1._real64, grid%x0, grid%dx, absrb, ist, last)
    deallocate(c, c0)
  endif
end subroutine oxygen_fundamental_continuum


subroutine oxygen_nir1_continuum(jrad, wk, grid, &
                                 tave, pave, x_vmr_o2, x_vmr_n2, x_vmr_h2o, &
                                 path, absrb)

  integer, intent(in) :: jrad
  real(kind=real64), intent(in) :: wk !< Oxygen number density [cm-3].
  real(kind=real64), intent(in) :: tave !< Temperature [K].
  real(kind=real64), intent(in) :: pave !< Pressure [mb].
  real(kind=real64), intent(in) :: x_vmr_o2, x_vmr_n2, x_vmr_h2o
  type(SpectralGrid), intent(in) :: grid !< Spectral grid.
  character(len=*), intent(in) :: path !< Path to dataset.
  real(kind=real64), dimension(:), intent(inout) :: absrb !< Extinction [cm-1].

  integer :: nptc, j, ist, last
  real(kind=real64) :: a_o2, a_n2, a_h2o, tau_fac, v1c, v2c, dvc, v1ss, v2ss, vj
  real(kind=real64), dimension(:), allocatable :: c, c0

  if (grid%xn .gt. 7536._real64 .and. grid%x0 .lt. 8500._real64) then
    a_o2  = 1._real64/0.446_real64
    a_n2  = 0.3_real64/0.446_real64
    a_h2o = 1._real64
    tau_fac = (wk/loschmidt)*(pave/p0)*(t_273/tave)* &
              (a_o2*x_vmr_o2 + a_n2*x_vmr_n2 + a_h2o*x_vmr_h2o) !unitless
    call o2inf1(v1c, v2c, dvc, nptc, c0, v1ss, v2ss, grid%x0, grid%xn, path)
    allocate(c(nptc))
    do j = 1, nptc
      c(j) = tau_fac*c0(j) !unitless.
      vj = v1c + dvc*real(j - 1) ![cm-1].
      if (jrad .eq. 1) c(j) = c(j)*radfn(vj, tave/radcn2) ![cm-1].
    enddo
    call pre_xint(v1ss, v2ss, grid%x0, grid%dx, grid%n, ist, last)
    call xint(v1c, v2c, dvc, c, 1._real64, grid%x0, grid%dx, absrb, ist, last)
    deallocate(c, c0)
  endif
end subroutine oxygen_nir1_continuum


subroutine oxygen_nir2_continuum(jrad, wk, wtot, grid, tave, pave, absrb)

  integer, intent(in) :: jrad
  real(kind=real64), intent(in) :: wk !< Oxygen number density [cm-3].
  real(kind=real64), intent(in) :: wtot !< Total air number density [cm-3].
  real(kind=real64), intent(in) :: tave !< Temperture [K].
  real(kind=real64), intent(in) :: pave !< Pressure [mb].
  type(SpectralGrid), intent(in) :: grid !< Spectral grid.
  real(kind=real64), dimension(:), intent(inout) :: absrb !< Extinction [cm-1].

  integer :: nptc, j, ist, last
  real(kind=real64) :: v1c, v2c, dvc, v1ss, v2ss, wo2, adjwo2, vj
  real(kind=real64), dimension(:), allocatable :: c, c0

  if (grid%xn .gt. 9100._real64 .and. grid%x0 .lt. 11000._real64) then
    call o2inf2(v1c, v2c, dvc, nptc, c0, v1ss, v2ss, grid%x0, grid%xn)
    wo2 = wk*1.e-20_real64*(pave/p0)*(t0/tave) ![cm-3].
    adjwo2 = (wk/wtot)*(1._real64/0.209_real64)*wo2 ![cm-3].
    allocate(c(nptc))
    do j = 1, nptc
      c(j) = c0(j)*adjwo2 !unitless.
      vj = v1c + dvc*real(j - 1) ![cm-1].
      if (jrad .eq. 1) c(j) = c(j)*radfn(vj, tave/radcn2) ![cm-1].
    enddo
    call pre_xint(v1ss, v2ss, grid%x0, grid%dx, grid%n, ist, last)
    call xint(v1c, v2c, dvc, c, 1._real64, grid%x0, grid%dx, absrb, ist, last)
    deallocate(c, c0)
  endif
end subroutine oxygen_nir2_continuum


subroutine oxygen_nir3_continuum(jrad, wk, grid, tave, pave, path, absrb)

  integer, intent(in) :: jrad
  real(kind=real64), intent(in) :: wk !< Oxygen number density [cm-3].
  real(kind=real64), intent(in) :: tave !< Temperature [K].
  real(kind=real64), intent(in) :: pave !< Pressure [mb].
  type(SpectralGrid), intent(in) :: grid !< Spectral grid.
  character(len=*), intent(in) :: path !< Path to dataset.
  real(kind=real64), dimension(:), intent(inout) :: absrb !< Extinction [cm-1].

  integer :: nptc, j, ist, last
  real(kind=real64) :: tau_fac, v1c, v2c, dvc, v1ss, v2ss, vj
  real(kind=real64), dimension(:), allocatable :: c, c0

  if (grid%xn .gt. 12961.5_real64 .and. grid%x0 .lt. 13221.5_real64) then
    tau_fac = (wk/loschmidt)*(pave/p0)*(t_273/tave) !unitless
    call o2inf3(v1c, v2c, dvc, nptc, c0, v1ss, v2ss, grid%x0, grid%xn, path)
    allocate(c(nptc))
    do j = 1, nptc
      c(j) = tau_fac*c0(j) !unitless.
      vj = v1c + dvc*real(j - 1) ![cm-1].
      if (jrad .eq. 1) c(j) = c(j)*radfn(vj, tave/radcn2) ![cm-1].
    enddo
    call pre_xint(v1ss, v2ss, grid%x0, grid%dx, grid%n, ist, last)
    call xint(v1c, v2c, dvc, c, 1._real64, grid%x0, grid%dx, absrb, ist, last)
    deallocate(c, c0)
  endif
end subroutine oxygen_nir3_continuum


subroutine oxygen_visible_continuum(jrad, wk, wtot, grid, tave, pave, path, absrb)

  integer, intent(in) :: jrad
  real(kind=real64), intent(in) :: wk !< Oxygen number density [cm-3].
  real(kind=real64), intent(in) :: wtot !< Total air number density [cm-3].
  real(kind=real64), intent(in) :: tave !< Temperature [K].
  real(kind=real64), intent(in) :: pave !< Pressure [mb].
  type(SpectralGrid), intent(in) :: grid !< Spectral grid.
  character(len=*), intent(in) :: path !< Path to dataset.
  real(kind=real64), dimension(:), intent(inout) :: absrb !< Extinction [cm-1].

  integer :: nptc, j, ist, last
  real(kind=real64) :: wo2, chio2, adjwo2, v1c, v2c, dvc, v1ss, v2ss, vj
  real(kind=real64), dimension(:), allocatable :: c, c0

  if (grid%xn .gt. 15000._real64 .and. grid%x0 .lt. 29870._real64) then
    wo2 = wk*1.e-20_real64*(pave/p0)*(t_273/tave) ![cm-3].
    chio2 = wk/wtot !unitless.
    adjwo2 = chio2*wo2 ![cm-3].
    call o2_vis(v1c, v2c, dvc, nptc, c0, v1ss, v2ss, grid%x0, grid%xn, path)
    allocate(c(nptc))
    do j = 1, nptc
      c(j) = c0(j)*adjwo2 !unitless.
      vj = v1c + dvc*real(j - 1) ![cm-1].
      if (jrad .eq. 1) c(j) = c(j)*radfn(vj, tave/radcn2) ![cm-1].
    enddo
    call pre_xint(v1ss, v2ss, grid%x0, grid%dx, grid%n, ist, last)
    call xint(v1c, v2c, dvc, c, 1._real64, grid%x0, grid%dx, absrb, ist, last)
    deallocate(c, c0)
  endif
end subroutine oxygen_visible_continuum


subroutine oxygen_herzberg_continuum(jrad, wk, grid, tave, pave, absrb)

  integer, intent(in) :: jrad
  real(kind=real64), intent(in) :: wk !< Oxygen number density [cm-3].
  real(kind=real64), intent(in) :: tave !< Temperature [K].
  real(kind=real64), intent(in) :: pave !< Pressure [mb].
  type(SpectralGrid), intent(in) :: grid !< Spectral grid.
  real(kind=real64), dimension(:), intent(inout) :: absrb !< Extinction [cm-1].

  integer :: nptc, j, ist, last
  real(kind=real64) :: wo2, v1c, v2c, dvc, v1ss, v2ss, vj
  real(kind=real64), dimension(:), allocatable :: c, c0

  if (grid%xn .gt. 36000._real64) then
    wo2 = wk*1.e-20_real64 ![cm-3].
    call o2herz(v1c, v2c, dvc, nptc, c0, tave, pave, v1ss, v2ss, grid%x0, grid%xn)
    allocate(c(nptc))
    do j = 1, nptc
      c(j) = c0(j)*wo2 !unitless.
      vj = v1c + dvc*real(j - 1)
      if (jrad .eq. 1) c(j) = c(j)*radfn(vj, tave/radcn2) ![cm-1].
    enddo
    call pre_xint(v1ss, v2ss, grid%x0, grid%dx, grid%n, ist, last)
    call xint(v1c, v2c, dvc, c, 1._real64, grid%x0, grid%dx, absrb, ist, last)
    deallocate(c, c0)
  endif
end subroutine oxygen_herzberg_continuum


subroutine oxygen_uv_continuum(jrad, wk, grid, tave, path, absrb)

  integer, intent(in) :: jrad
  real(kind=real64), intent(in) :: wk !< Oxygen number density [cm-3].
  real(kind=real64), intent(in) :: tave !< Temperature [K].
  type(SpectralGrid), intent(in) :: grid !< Spectral grid.
  character(len=*), intent(in) :: path !< Path to dataset.
  real(kind=real64), dimension(:), intent(inout) :: absrb !< Extinction [cm-1].

  integer :: nptc, j, ist, last
  real(kind=real64) :: wo2, v1c, v2c, dvc, v1ss, v2ss, vj
  real(kind=real64), dimension(:), allocatable :: c, c0

  if (grid%xn .gt. 56740._real64) then
    wo2 = wk*1.e-20_real64 ![cm-3].
    call o2fuv(v1c, v2c, dvc, nptc, c0, v1ss, v2ss, grid%x0, grid%xn, path)
    allocate(c(nptc))
    do j = 1, nptc
      vj = v1c + dvc*real(j - 1)
      c(j) = c0(j)*wo2 !unitless.
      if (jrad .eq. 1) c(j) = c(j)*radfn(vj, tave/radcn2) ![cm-1].
    enddo
    call pre_xint(v1ss, v2ss, grid%x0, grid%dx, grid%n, ist, last)
    call xint(v1c, v2c, dvc, c, 1._real64, grid%x0, grid%dx, absrb, ist, last)
    deallocate(c, c0)
  endif
end subroutine oxygen_uv_continuum


end module o2_continuum
