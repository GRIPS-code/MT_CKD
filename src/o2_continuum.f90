module o2_continuum
use, intrinsic :: iso_fortran_env, only: real64
use netcdf
implicit none


public :: o2_ver_1, o2_vis, o2fuv, o2herz, o2inf1, o2inf2, o2inf3


contains


subroutine o2_ver_1(v1c, v2c, dvc, nptc, c, t, v1ss, v2ss, v1abs, v2abs, path)

  real(kind=real64), intent(out) :: v1c, v2c, dvc
  integer, intent(out) :: nptc
  real(kind=real64), dimension(:), intent(inout) :: c
  real(kind=real64), intent(in) :: t
  real(kind=real64), intent(out) :: v1ss, v2ss
  real(kind=real64), intent(in) :: v1abs, v2abs
  character(len=*), intent(in) :: path

  integer :: dimid, err, i, i1, i2, j, ncid, npts, varid
  real(kind=real64) :: dvs, factor, vj, v1s, v2s, xktfac, xlosmt
  real(kind=real64), dimension(:), allocatable :: xo2, xo2t

!     oxygen collision induced fundamental
!     f. thibault, v. menoux, r. le doucen, l. rosenman, j.-m. hartmann,
!     and ch. boulet
!     infrared collision-induced absorption by o2 near 6.4 microns for
!     atmospheric applications: measurements and emprirical modeling,
!     appl. optics, 35, 5911-5917, (1996).

  !Read data from netcdf file.
  err = nf90_open(path, nf90_nowrite, ncid)
  err = nf90_inq_varid(ncid, "o2_f", varid)
  err = nf90_get_att(ncid, varid, "wavenumber_lower_bound", v1s)
  err = nf90_get_att(ncid, varid, "wavenumber_upper_bound", v2s)
  err = nf90_get_att(ncid, varid, "wavenumber_resolution", dvs)
  err = nf90_inq_dimid(ncid, "n10", dimid)
  err = nf90_inquire_dimension(ncid, dimid, len=npts)
  allocate(xo2(npts), xo2t(npts))
  err = nf90_get_var(ncid, varid, xo2)
  err = nf90_inq_varid(ncid, "o2_t", varid)
  err = nf90_get_var(ncid, varid, xo2t)
  err = nf90_close(ncid)

  xlosmt = 2.68675e+19
  xktfac = (1./296.) - (1./t)
!
!     correct formulation for consistency with lblrtm:
  factor = (1.e+20/xlosmt)
!
!     a factor of 0.21, the mixing ratio of oxygen, in the thibault et
!     al. formulation is not included here.  this factor is in the
!     column amount.

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
  nptc = i2-i1+3
  if (nptc .gt. npts) nptc = npts + 4
  v2c = v1c + dvs*real(nptc - 1)

  do j = 1, nptc
    i = i1 + (j - 1)
    c(j) = 0.
    if ((i .lt. 1) .or. (i .gt. npts)) continue
    vj = v1c + dvc*real(j - 1)
!     the radiation field is removed with 1/vj
    c(j) = factor*xo2(i)*exp(xo2t(i)*xktfac)/vj
  enddo
  deallocate(xo2, xo2t)
end subroutine o2_ver_1


subroutine o2inf1(v1c, v2c, dvc, nptc, c, v1ss, v2ss, v1abs, v2abs, path)

  real(kind=real64), intent(out) :: v1c, v2c, dvc
  integer, intent(out) :: nptc
  real(kind=real64), dimension(:), intent(inout) :: c
  real(kind=real64), intent(out) :: v1ss, v2ss
  real(kind=real64), intent(in) :: v1abs, v2abs
  character(len=*), intent(in) :: path

  integer :: dimid, err, i, i1, i2, j, ncid, npts, varid
  real(kind=real64) :: dvs, vj, v1s, v2s
  real(kind=real64), dimension(:), allocatable :: xo2inf1

!        o2 continuum formulated by mate et al. over the spectral region
!        7550-8486 cm-1:  "absolute intensities for the o2 1.27 micron
!        continuum absorption", b. mate, c. lugez, g.t. fraser, and
!        w.j. lafferty, j. geophys. res., 104, 30,585-30,590, 1999.
!
!        the units of these continua coefficients are
!         1 / (amagat_o2*amagat_air).
!        also, refer to the paper "observed  atmospheric
!        collision induced absorption in near infrared oxygen bands",
!        mlawer, clough, brown, stephen, landry, goldman, & murcray,
!        journal of geophysical research (1997).

  !Read data from netcdf file.
  err = nf90_open(path, nf90_nowrite, ncid)
  err = nf90_inq_varid(ncid, "o2_inf1", varid)
  err = nf90_get_att(ncid, varid, "wavenumber_lower_bound", v1s)
  err = nf90_get_att(ncid, varid, "wavenumber_upper_bound", v2s)
  err = nf90_get_att(ncid, varid, "wavenumber_resolution", dvs)
  err = nf90_inq_dimid(ncid, "n11", dimid)
  err = nf90_inquire_dimension(ncid, dimid, len=npts)
  allocate(xo2inf1(npts))
  err = nf90_get_var(ncid, varid, xo2inf1)
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

  do j = 1, nptc
    i = i1 + (j - 1)
    c(j) = 0.
    if ((i .lt. 1) .or. (i .gt. npts)) continue
    vj = v1c + dvc*real(j - 1)
    c(j) = xo2inf1(i)/vj
  enddo
  deallocate(xo2inf1)
end subroutine o2inf1


subroutine o2inf2(v1c, v2c, dvc, nptc, c, v1ss, v2ss, v1abs, v2abs)

  real(kind=real64), intent(out) :: v1c, v2c, dvc
  integer, intent(out) :: nptc
  real(kind=real64), dimension(:), intent(inout) :: c
  real(kind=real64), intent(out) :: v1ss, v2ss
  real(kind=real64), intent(in) :: v1abs, v2abs

  integer :: j
  real(kind=real64) :: damp1, damp2, dvs, dv1, dv2, hw1, hw2, o2inf, s1, s2, v1_osc, &
                       v2_osc, vj, v1s, v2s

  v1_osc = 9375.
  hw1 = 58.96
  v2_osc = 9439.
  hw2 = 45.04
  s1 = 1.166e-04
  s2 = 3.086e-05

  v1s = 9100.
  v2s = 11000.
  dvs = 2.
  dvc = dvs
  v1ss = v1s
  v2ss = v2s

  v1c = v1abs - dvc
  v2c = v2abs + dvc

  ! the following lines prevent a possible problem that can only occur if
  ! v2abs-v1abs >> 2000 cm-1 (i.e. in the standalone continuum code).
  if (v1c .lt. v1s) v1c = v1s - 2.*dvs
  if (v2c .gt. v2s) v2c = v2s + 2.*dvs
  nptc = (v2c - v1c)/dvc + 3.01
  v2c = v1c + dvc*real(nptc - 1)

  do j = 1, nptc
    c(j) = 0.
    vj = v1c + dvc*real(j - 1)
    if ((vj .gt. v1s) .and. (vj .lt. v2s)) then
      dv1 = vj - v1_osc
      dv2 = vj - v2_osc
      if (dv1 .lt. 0.0) then
        damp1 = exp(dv1/176.1)
      else
        damp1 = 1.0
      endif
      if (dv2 .lt. 0.0) then
        damp2 = exp(dv2/176.1)
      else
        damp2 = 1.0
      endif
      o2inf = 0.31831*(((s1*damp1/hw1)/(1. + (dv1/hw1)**2)) &
              + ((s2*damp2/hw2)/(1. + (dv2/hw2)**2)))*1.054
      c(j) = o2inf/vj
    endif
  enddo
end subroutine o2inf2


subroutine o2inf3(v1c, v2c, dvc, nptc, c, v1ss, v2ss, v1abs, v2abs, path)

  real(kind=real64), intent(out) :: v1c, v2c, dvc
  integer, intent(out) :: nptc
  real(kind=real64), dimension(:), intent(inout) :: c
  real(kind=real64), intent(out) :: v1ss, v2ss
  real(kind=real64), intent(in) :: v1abs, v2abs
  character(len=*), intent(in) :: path

  integer :: dimid, err, i, i1, i2, j, ncid, npts, varid
  real(kind=real64) :: dvs, vj, v1s, v2s
  real(kind=real64), dimension(:), allocatable :: xo2inf3

!
!        O2 A-band continuum formulated by Mlawer based on solar FTS measurements.
!        See Payne et al. (2020) "Absorption Coefficient (ABSCO) Tables for the
!        Orbiting Carbon Observatories: Version 5.1", JQSRT
!
!        Units of these coefficients are 1 / (amagat_O2*amagat_air)
!        Spectral range of coefficients is 12961.5 - 13221.5 cm-1.
!
  !Read data from netcdf file.
  err = nf90_open(path, nf90_nowrite, ncid)
  err = nf90_inq_varid(ncid, "o2_inf3", varid)
  err = nf90_get_att(ncid, varid, "wavenumber_lower_bound", v1s)
  err = nf90_get_att(ncid, varid, "wavenumber_upper_bound", v2s)
  err = nf90_get_att(ncid, varid, "wavenumber_resolution", dvs)
  err = nf90_inq_dimid(ncid, "n12", dimid)
  err = nf90_inquire_dimension(ncid, dimid, len=npts)
  allocate(xo2inf3(npts))
  err = nf90_get_var(ncid, varid, xo2inf3)
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

  do j = 1, nptc
    i = i1 + (j - 1)
    c(j) = 0.
    if ((i .lt. 1) .or. (i .gt. npts)) continue
    vj = v1c + dvc*real(j - 1)
    c(j) = xo2inf3(i)/vj
  enddo
  deallocate(xo2inf3)
end subroutine o2inf3


subroutine o2_vis(v1c, v2c, dvc, nptc, c, v1ss, v2ss, v1abs, v2abs, path)

  real(kind=real64), intent(out) :: v1c, v2c, dvc
  integer, intent(out) :: nptc
  real(kind=real64), dimension(:), intent(inout) :: c
  real(kind=real64), intent(out) :: v1ss, v2ss
  real(kind=real64), intent(in) :: v1abs, v2abs
  character(len=*), intent(in) :: path

  integer :: dimid, err, i, i1, i2, j, ncid, npts, varid
  real(kind=real64) :: factor, vj, xlosmt
  real(kind=real64) :: v1s, v2s, dvs
  real(kind=real64), dimension(:), allocatable :: s

!
!        o2 continuum formulated by greenblatt et al. over the spectral
!        region 8797-29870 cm-1:  "absorption coefficients of oxygen
!        between 330 and 1140 nm, g.d. greenblatt, j.j. orlando, j.b.
!        burkholder and a.r. ravishabkara,  j. geophys. res., 95,
!        18577-18582, 1990.
!
!        the units conversion  is to (cm^2/molec)/atm(o2)
!
!      these are the conditions reported in the paper by greenblatt et
!      al. for the spectrum of fig. 1.
!
!     conditions:  55 atm.; 296 k; 89.5 cm path
!
  !Read data from netcdf file.
  err = nf90_open(path, nf90_nowrite, ncid)
  err = nf90_inq_varid(ncid, "o2_invis", varid)
  err = nf90_get_att(ncid, varid, "wavenumber_lower_bound", v1s)
  err = nf90_get_att(ncid, varid, "wavenumber_upper_bound", v2s)
  err = nf90_get_att(ncid, varid, "wavenumber_resolution", dvs)
  err = nf90_inq_dimid(ncid, "n13", dimid)
  err = nf90_inquire_dimension(ncid, dimid, len=npts)
  allocate(s(npts))
  err = nf90_get_var(ncid, varid, s)
  err = nf90_close(ncid)

  xlosmt = 2.68675e+19
  factor = 1./((xlosmt*1.e-20*(55.*273./296.)**2)*89.5)

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
  nptc = i2 - i1+3
  if (nptc .gt. npts) nptc = npts + 4
  v2c = v1c + dvs*real(nptc - 1)

  do j = 1, nptc
    i = i1 + (j - 1)
    c(j) = 0.
    if ((i .lt. 1) .or. (i .gt. npts)) continue
    vj = v1c + dvc*real(j - 1)
    c(j) = factor*s(i)/vj
  enddo
  deallocate(s)
end subroutine o2_vis


subroutine o2herz(v1c, v2c, dvc, nptc, c, t, p, v1ss, v2ss, v1abs, v2abs)

  real(kind=real64), intent(out) :: v1c, v2c, dvc
  integer, intent(out) :: nptc
  real(kind=real64), dimension(:), intent(inout) :: c
  real(kind=real64), intent(in) :: t, p
  real(kind=real64), intent(out) :: v1ss, v2ss
  real(kind=real64), intent(in) :: v1abs, v2abs

  integer :: i, i1, i2, j
  real(kind=real64) :: dvs, herz, vj, v1s

  v1s = 36000.
  v1ss = v1s
  v2ss = 99999.

  dvs = 10.
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
!         IF (NPTC.GT.NPTS) NPTC=NPTS+4
!        mja, 10-27-2011 - this seems to be redundant as the
!        Herzberg O2 continuum is a function, not block data
  v2c = v1c + dvs*real(nptc - 1)

  do j = 1, nptc
    i = i1 + (j - 1)
    c(j) = 0.
    if (i .lt. 1) continue
    vj = v1c + dvc*real(j - 1)
    call hertda(herz, vj)
    call herprs(herz, t, p)
    c(j) = herz/vj
  enddo
end subroutine o2herz


subroutine hertda(herz, v)

  real(kind=real64), intent(out) :: herz
  real(kind=real64), intent(in) :: v

  real(kind=real64) :: corr, yratio
!
!     HERZBERG O2 ABSORPTION
!     HALL,1987 PRIVATE COMMUNICATION, BASED ON:
!
!     REF. JOHNSTON, ET AL., JGR,89,11661-11665,1984
!          NICOLET, 1987 (RECENT STUDIES IN ATOMIC
!                         & MOLECULAR PROCESSES,
!                         PLENUM PUBLISHING CORP, NY 1987)
!
!     AND YOSHINO, ET AL., 1988 (PREPRINT OF "IMPROVED ABSORPTION
!          CROSS SECTIONS OF OXYGEN IN THE WAVELENGTH REGION 205-240NM
!          OF THE HERZBERG CONTINUUM")
!
!     **** NOTE:  CROSS SECTION AT 0 PRESSURE  ***
!     THE PRESSURE DEPENDENT TERM IS IN SUBROUTINE HERPRS
!
!C    COMMON /CNSTNS/ PI,CA,DEG,GCAIR,BIGNUM,BIGEXP
!
  herz = 0.0
  if (v .le. 36000.00) return
!
!     EXTRAPOLATE SMOOTHLY THROUGH THE HERZBERG BAND REGION
!     NOTE: HERZBERG BANDS ARE NOT CORRECTLY INCLUDED
!
  corr = 0.
!      IF (V.LE.40000.) CORR = ((40000.-V)/4000.)*7.917E-27
!     factor of 1.e-20 removed; put in front factor
  if (v .le. 40000.) corr = ((40000. - v)/4000.)*7.917e-07
!
!     UNITS ARE (CM2)
!
!     HALL'S NEW HERZBERG  (LEAST SQRS FIT, LN(P))
!
!     YRATIO=2048.7/WL(I)  ****IN ANGSTOMS****
!           =.20487/WN(I)     IN MICRONS
!           =WCM(I)/48811.0   IN CM-1
!
  yratio = v/48811.0
!     HERZ = 6.884E-24*(YRATIO)*EXP(-69.738*( LOG(YRATIO))**2)-CORR
!     factor of 1.e-20 removed; put in front factor
  herz = 6.884e-04*(yratio)*exp(-69.738*(log(yratio))**2) - corr
end subroutine hertda


subroutine herprs(herz, t, p)

  real(kind=real64), intent(inout) :: herz
  real(kind=real64), intent(in) :: t, p
!
!     CORRECT THE HERZBERG CONTINUUM CROSS SECTION FOR PRESSURE
!     DEPENDENCE; BASED ON SHARDANAND, JQRST, 18, 525-530, 1977.
!                 FOR UN2| BROADENING
!                 AND YOSHINO ET AL 1988 FOR UO2| BROADENING
!
!     PO2= PARTIAL PRESSURE OF O2
!     PN2= PARTIAL PRESSURE OF N2; BN2=.45*BO2
!
!     DATA BO2 / 1.72E-3 /
!
!     Changed in Herzberg continuum pressure,
!     Reference:
!     "Atmospheric Propagation in the UV, Visible, IR and MM-wave
!     Region and Related Systems Aspects".
!     G.P. Anderson,F.X. Kneizys, E.P. Shettle, L.W. Abreu,
!     J.H. Chetwynd, R.E. Huffman, and L.A. Hall; Conference
!     Proceedings No. 454 of the Advisory Group for Aerospace
!     Research & Development; 1990.
!     NOTE:  THE HERZBERG CONTINUUM OBEYS BEER'S LAW
!            OPTICAL DEPTH(TOTAL)=SUM OVER LAYER O.D.(I)
!
!     BO2= RATIO OF SIGMA(O2-O2)/(SIGMA(O2)) * 760(TORR)*.2095
!     BN2=.45*BO2= RATIO OF SIGMA(O2-N2)/(SIGMA(O2)) * 760(TORR)*.78
!
!     BO2*760*(.2095+.45*.78) = .73 , AS BELOW
!
!     Changed Herzberg continuum pressure (see above reference)
!
!     BO2*760*(.2095+.45*.78) = .83 , AS BELOW
  herz = herz*(1.+.83*(p/1013.)*(273.16/t))
end subroutine herprs


!> @brief O2 continuum in the far-uv is from lu et al. (2010).
subroutine o2fuv(v1c, v2c, dvc, nptc, c, v1ss, v2ss, v1abs, v2abs, path)

  real(kind=real64), intent(out) :: v1c, v2c, dvc
  integer, intent(out) :: nptc
  real(kind=real64), dimension(:), intent(inout) :: c
  real(kind=real64), intent(out) :: v1ss, v2ss
  real(kind=real64), intent(in) :: v1abs, v2abs
  character(len=*), intent(in) :: path

  integer :: dimid, err, i, i1, i2, j, ncid, npts, varid
  real(kind=real64) :: dvs, vj, v1s, v2s
  real(kind=real64), dimension(:), allocatable  :: s

  !Read data from netcdf file.
  err = nf90_open(path, nf90_nowrite, ncid)
  err = nf90_inq_varid(ncid, "o2_infuv", varid)
  err = nf90_get_att(ncid, varid, "wavenumber_lower_bound", v1s)
  err = nf90_get_att(ncid, varid, "wavenumber_upper_bound", v2s)
  err = nf90_get_att(ncid, varid, "wavenumber_resolution", dvs)
  err = nf90_inq_dimid(ncid, "n14", dimid)
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
    i1 = (v1c - v1s)/dvs + 1.e-5
  endif
  v1c = v1s + dvs*real(i1 - 1)

  i2 = (v2c - v1s)/dvs + 1.e-5

  nptc = i2 - i1 + 3
  if (nptc .gt. npts) nptc = npts + 4
  v2c = v1c + dvs*real(nptc - 1)

  do j = 1, nptc
    i = i1 + (j - 1)
    c(j) = 0.
    if ((i .lt. 1) .or. (i .gt. npts)) continue
    vj = v1c + dvc*real(j - 1)
    c(j) = s(i)/vj
  enddo
  deallocate(s)
end subroutine o2fuv


end module o2_continuum
