module n2_continuum
use netcdf
implicit none


public :: xn2_r, n2_ver_1, n2_overtone1


contains


subroutine xn2_r(v1c, v2c, dvc, nptc, c, fo2, tave, v1ss, v2ss, v1abs, v2abs, path)

  real, intent(out) :: v1c, v2c, dvc
  integer, intent(out) :: nptc
  real, dimension(:), allocatable, intent(inout) :: c, fo2
  real, intent(in) :: tave
  real, intent(out) :: v1ss, v2ss
  real, intent(in) :: v1abs, v2abs
  character(len=*), intent(in) :: path

  integer :: dimid, err, i, i1, i2, j, ncid, nptb, npts, varid
  real :: dvb, dvs, sf_t, tfac, v1b, v2b, v1s, v2s, xn2, xo2
  real, dimension(:), allocatable :: c_220, c_296, sf_220, sf_296

!
!     Model used:
!      Borysow, A, and L. Frommhold, "Collision-induced
!         rototranslational absorption spectra of N2-N2
!         pairs for temperatures from 50 to 300 K", The
!         Astrophysical Journal, 311, 1043-1057, 1986.
!
!     Updated 2004/09/22 based on:
!
!      Boissoles, J., C. Boulet, R.H. Tipping, A. Brown and Q. Ma,
!         Theoretical CAlculations of the Translation-Rotation
!         Collision-Induced Absorption in N2-N2, O2-O2 and N2-O2 Pairs,
!         J.Quant. Spec. Rad. Transfer, 82,505 (2003).
!
!     The scale factors are reported to account for the efect of o2-o2
!     and n2-o2 collision induced effects.
!     The values for scale factor values (sf296) for 296K are based on
!     linear interpolation of Boissoles at al. values at 250K and 300K

  !Read data from netcdf file.
  err = nf90_open(path, nf90_nowrite, ncid)
  err = nf90_inq_varid(ncid, "ct_296", varid)
  err = nf90_get_att(ncid, varid, "wavenumber_lower_bound", v1s)
  err = nf90_get_att(ncid, varid, "wavenumber_upper_bound", v2s)
  err = nf90_get_att(ncid, varid, "wavenumber_resolution", dvs)
  err = nf90_inq_dimid(ncid, "n4", dimid)
  err = nf90_inquire_dimension(ncid, dimid, len=npts)
  allocate(c_296(npts))
  err = nf90_get_var(ncid, varid, c_296)
  err = nf90_inq_varid(ncid, "sf_296", varid)
  allocate(sf_296(npts))
  err = nf90_get_var(ncid, varid, sf_296)
  err = nf90_inq_varid(ncid, "ct_220", varid)
  err = nf90_get_att(ncid, varid, "wavenumber_lower_bound", v1b)
  err = nf90_get_att(ncid, varid, "wavenumber_upper_bound", v2b)
  err = nf90_get_att(ncid, varid, "wavenumber_resolution", dvb)
  err = nf90_inq_dimid(ncid, "n4", dimid)
  err = nf90_inquire_dimension(ncid, dimid, len=nptb)
  allocate(c_220(npts))
  err = nf90_get_var(ncid, varid, c_220)
  err = nf90_inq_varid(ncid, "sf_220", varid)
  allocate(sf_220(npts))
  err = nf90_get_var(ncid, varid, sf_220)
  err = nf90_close(ncid)

  xo2 = 0.21
  xn2 = 0.79

  tfac = (tave - 296.)/(220. - 296.)

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

!*******  absorption coefficient in units of cm-1 amagat-2
  allocate(c(nptc), fo2(nptc))
  do j = 1, nptc
    i = i1 + (j - 1)
    c(j) = 0.
    if ((i .lt. 1) .or. (i .gt. npts)) continue
    c(j) = c_296(i)*((c_220(i)/c_296(i))**tfac)
    sf_t = sf_296(i)*((sf_220(i)/sf_296(i))**tfac)
!        correct for incorporation of air mixing ratios in sf
!        fo2 is now ~ the ratio of alpha(n2-o2)/alpha(n2-n2)
!        eq's 7 and 8 in the boissoles paper.
    fo2(j) = (sf_t - 1.)*(xn2)/(xo2)
  enddo
  deallocate(c_220, c_296, sf_220, sf_296)
end subroutine xn2_r


subroutine n2_ver_1(v1c, v2c, dvc, nptc, c, c1, c2, t, v1ss, v2ss, v1abs, v2abs, path)

  real, intent(out) :: v1c, v2c, dvc
  integer, intent(out) :: nptc
  real, dimension(:), allocatable, intent(inout) :: c, c1, c2
  real, intent(in) :: t
  real, intent(out) :: v1ss, v2ss
  real, intent(in) :: v1abs, v2abs
  character(len=*), intent(in) :: path

  integer :: dimid, err, i, i1, i2, j, ncid, npts, varid
  real :: a_o2, dvs, vj, v1s, v2s, xt_lin, xtfac
  real, dimension(:), allocatable :: xn2_272, xn2_228, a_h2o

  !Read data from netcdf file.
  err = nf90_open(path, nf90_nowrite, ncid)
  err = nf90_inq_varid(ncid, "xn2_272", varid)
  err = nf90_get_att(ncid, varid, "wavenumber_lower_bound", v1s)
  err = nf90_get_att(ncid, varid, "wavenumber_upper_bound", v2s)
  err = nf90_get_att(ncid, varid, "wavenumber_resolution", dvs)
  err = nf90_inq_dimid(ncid, "n5", dimid)
  err = nf90_inquire_dimension(ncid, dimid, len=npts)
  allocate(xn2_272(npts), xn2_228(npts), a_h2o(npts))
  err = nf90_get_var(ncid, varid, xn2_272)
  err = nf90_inq_varid(ncid, "xn2_228", varid)
  err = nf90_get_var(ncid, varid, xn2_228)
  err = nf90_inq_varid(ncid, "a_h2o", varid)
  err = nf90_get_var(ncid, varid, a_h2o)
  err = nf90_close(ncid)

!     Nitrogen Collision Induced Fundamental

!     Lafferty, W.J., A.M. Solodov,A. Weber, W.B. Olson and J._M.
!        Hartmann, Infrared collision-induced absorption by N2 near 4.3
!        microns for atmospheric applications: measurements and
!         emprirical modeling, Appl. Optics, 35, 5911-5917, (1996).
!
!     mt_ckd_2.8: Coefficients for N2-H2O relative efficiency determined
!     by Mlawer and Alvarado based primarily on measurements from
!     Baranov and Lafferty (2012). These coefficients were derived
!     simultaneously with water vapor foreign and self continuum
!     coefficients from 1800-2600 cm-1 using IASI measurements as in
!     Alvarado et al. (2012).
!
  xtfac = ((1./t) - (1./272.))/((1./228.) - (1./272.))
  xt_lin = (t - 272.)/(228. - 272.)
!
!     a_o2  represents the relative broadening efficiency of o2
  a_o2  = 1.294 - 0.4545*t/296.

!     a_h2o represents the relative broadening efficiency of h2o.  It
!     has spectral dependence and is stored on same grid as xn2.

!     The absorption coefficients from the Lafferty et al. reference
!     are for pure nitrogen (absorber and broadener)
!
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

  allocate(c(nptc), c1(nptc), c2(nptc))
  do j = 1, nptc
    i = i1 + (j - 1)
    c(j) = 0.
    if ((i .lt. 1) .or. (i .gt. npts)) continue
    vj = v1c + dvc*real(j - 1)
    if ((xn2_272(i) .gt. 0.) .and. (xn2_228(i) .gt. 0.)) then
!           logarithmic interpolation in reciprical of temperature
       c(j) = xn2_272(i)*(xn2_228(i)/xn2_272(i))**xtfac
    else
!           linear interpolation  (xn2_272 or xn2_228 = 0 to get here)
      c(j) = xn2_272(i) + (xn2_228(i) - xn2_272(i))*xt_lin
    endif
!     the radiation field is removed with 1/vj
    c(j) = c(j)/vj
    c1(j) = a_o2*c(j)
!     the factor 9/7 is a modification to a preliminary formulation
!     of n2-h2o absorption.
    c2(j) = (9./7.)*a_h2o(i)*c(j)
  enddo
  deallocate(xn2_272, xn2_228, a_h2o)
end subroutine n2_ver_1


subroutine n2_overtone1(v1c, v2c, dvc, nptc, c, v1ss, v2ss, v1abs, v2abs, path)

  real, intent(out) :: v1c, v2c, dvc
  integer, intent(out) :: nptc
  real, dimension(:), allocatable, intent(inout) :: c
  real, intent(out) :: v1ss, v2ss
  real, intent(in) :: v1abs, v2abs
  character(len=*), intent(in) :: path

  integer :: dimid, err, i, i1, i2, j, ncid, npts, varid
  real :: dvs, vj, v1s, v2s
  real, dimension(:), allocatable :: xn2

  !Read data from netcdf file.
  err = nf90_open(path, nf90_nowrite, ncid)
  err = nf90_inq_varid(ncid, "xn2", varid)
  err = nf90_get_att(ncid, varid, "wavenumber_lower_bound", v1s)
  err = nf90_get_att(ncid, varid, "wavenumber_upper_bound", v2s)
  err = nf90_get_att(ncid, varid, "wavenumber_resolution", dvs)
  err = nf90_inq_dimid(ncid, "n6", dimid)
  err = nf90_inquire_dimension(ncid, dimid, len=npts)
  allocate(xn2(npts))
  err = nf90_get_var(ncid, varid, xn2)
  err = nf90_close(ncid)

!     nitrogen collision induced first overtone

!     shapiro and gush (1966) modified by mlawer and gombos (2015).

!     the absorption coefficients are for pure nitrogen (absorber and
!     broadener.
!
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
!
  allocate(c(nptc))
  do j = 1, nptc
    i = i1 + (j - 1)
    c(j) = 0.
    if ((i .lt. 1) .or. (i .gt. npts)) continue
    vj = v1c + dvc*real(j - 1)
    c(j) = xn2(i)
!     the radiation field is removed with 1/vj
    c(j) = c(j)/vj
  enddo
  deallocate(xn2)
end subroutine n2_overtone1


end module n2_continuum
