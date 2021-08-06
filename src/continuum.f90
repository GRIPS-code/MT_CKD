module continuum
use, intrinsic :: iso_fortran_env, only: real64
use utils, only: pre_xint, radfn, xint
use co2_continuum, only: frnco2
use h2o_continuum, only: sl296, sl260, frn296
use n2_continuum, only: xn2_r, n2_ver_1, n2_overtone1
use o2_continuum, only: o2_ver_1, o2_vis, o2fuv, o2herz, o2inf1, o2inf2, o2inf3
use o3_continuum, only: xo3chp, o3hht0, o3hht1, o3hht2, o3hhuv
implicit none


public :: contnm


contains


!contains the continuum data which is interpolated into the array absrb.
subroutine contnm(icflg, nmol, nptabs, jrad, xcnt, pave, tave, wbroad, v1abs, &
                  v2abs, dvabs, v1, v2, wk, npth, v1h, dvh, absrb, path, csh2o, cfh2o)

  integer, intent(in) :: icflg, nmol, nptabs, jrad
  real(kind=real64), dimension(7), intent(in) :: xcnt
  real(kind=real64), intent(in) :: pave, tave, wbroad, v1abs, v2abs, dvabs, v1, v2
  real(kind=real64), dimension(:), intent(in) :: wk
  integer, intent(out) :: npth
  real(kind=real64), intent(out) :: v1h, dvh
  real(kind=real64), dimension(:), intent(inout) :: absrb
  character(len=*), intent(in) :: path
  real(kind=real64), dimension(:), intent(inout) :: csh2o, cfh2o

  character(len=18) :: hvrcnt
  integer :: i, i_fix, j, jfac, m, nptabsc, nptc, ist, last, n_1, n_2, npto3, npt1, npt2
  integer, parameter :: ipts = 5050
  integer, parameter :: n_absrb = 5050
  integer, dimension(17) :: fscdid
  real(kind=real64) :: rhoave, xkt, amagat, wn2, wtot, x_vmr_h2o, x_vmr_o2, x_vmr_n2, &
                       v1absc, v2absc, dvabsc, h2o_fac, xself, rself, xfrgn, rfrgn, &
                       rfrgn_aj, tfac, v1c, v2c, dvc, v1ss, v2ss, vj, sh2o, f0, &
                       v0f1, hwsq1, beta1, c_1, c_2, beta2, fscal, vdelsq1, vdelmsq1, &
                       vf1, vmf1, vf2, c_f, xco2c, wco2, cfac, xo3cn, wo3, dt, tc, &
                       v1t1, v1t2, dvt1, v2t1, v2t2, dvt2, xo2cn, tau_fac, a_o2, a_n2, &
                       a_h2o, adjwo2, chio2, xn2cn, xrayl, conv_cm2mol, wo2, &
                       vrayleigh, xvrayleigh, ray_ext
  real(kind=real64), parameter :: p0 = 1013. ![mb].
  real(kind=real64), parameter :: t0 = 296. ![K].
  real(kind=real64), parameter :: xlosmt = 2.68675e19
  real(kind=real64), parameter :: radcn2 = 1.4387752 ![cm K].
  real(kind=real64), dimension(ipts) :: cself, cfrgn_aj
  real(kind=real64), dimension(n_absrb) :: sh2ot0, sh2ot1, fh2o, fco2, c0, c1, ct1, ct2, &
                                           absbsv
  real(kind=real64), dimension(6000) :: c
  real(kind=real64), dimension(500) :: xfacco2
  real(kind=real64), dimension(5150):: cch0, cch1, cch2, cn0, cn1, cn2
  real(kind=real64), dimension(-1:61) :: xfac_rhu

  data xfacco2/ &
        1.0000,0.9998,0.9997,0.9996,0.9995,0.9994,0.9992,0.9991, &
        0.9991,0.9990,0.9990,0.9989,0.9988,0.9988,0.9987,0.9986, &
        0.9985,0.9984,0.9983,0.9982,0.9981,0.9980,0.9979,0.9978, &
        0.9976,                                                  &
        0.9975,0.9973,0.9972,0.9970,0.9969,0.9967,0.9965,0.9963, &
        0.9961,0.9958,0.9956,0.9954,0.9951,0.9948,0.9946,0.9943, &
        0.9940,0.9936,0.9933,0.9929,0.9926,0.9922,0.9918,0.9913, &
        0.9909,                                                  &
        0.9904,0.9899,0.9894,0.9889,0.9884,0.9878,0.9872,0.9866, &
        0.9859,0.9853,0.9846,0.9838,0.9831,0.9823,0.9815,0.9806, &
        0.9798,0.9789,0.9779,0.9770,0.9759,0.9749,0.9738,0.9727, &
        0.9716,                                                  &
        0.9704,0.9691,0.9679,0.9666,0.9652,0.9638,0.9624,0.9609, &
        0.9594,0.9578,0.9562,0.9546,0.9529,0.9511,0.9493,0.9475, &
        0.9456,0.9436,0.9417,0.9396,0.9375,0.9354,0.9332,0.9310, &
        0.9287,                                                  &
        0.9264,0.9240,0.9216,0.9191,0.9166,0.9140,0.9114,0.9087, &
        0.9060,0.9032,0.9004,0.8976,0.8947,0.8917,0.8887,0.8857, &
        0.8827,0.8796,0.8764,0.8732,0.8700,0.8668,0.8635,0.8602, &
        0.8568,                                                  &
        0.8534,0.8500,0.8466,0.8432,0.8397,0.8362,0.8327,0.8292, &
        0.8257,0.8221,0.8186,0.8151,0.8115,0.8080,0.8044,0.8009, &
        0.7973,0.7938,0.7903,0.7868,0.7833,0.7799,0.7764,0.7730, &
        0.7697,                                                  &
        0.7663,0.7630,0.7597,0.7565,0.7533,0.7502,0.7471,0.7441, &
        0.7411,0.7382,0.7354,0.7326,0.7298,0.7272,0.7246,0.7221, &
        0.7197,0.7173,0.7150,0.7129,0.7108,0.7088,0.7068,0.7050, &
        0.7033,                                                  &
        0.7016,0.7001,0.6986,0.6973,0.6961,0.6949,0.6939,0.6930, &
        0.6921,0.6914,0.6908,0.6903,0.6899,0.6897,0.6895,0.6895, &
        0.6895,0.6895,0.6895,0.6895,0.6908,0.7014,0.7121,0.7227, &
        0.7552,                                                  &
        0.8071,0.8400,0.9012,0.9542,1.0044,1.0330,1.0554,1.0766, &
        1.0967,1.1160,1.1346,1.1525,1.1700,1.1869,1.2035,1.2196, &
        1.2354,1.2509,1.2662,1.2811,1.2958,1.3103,1.3245,1.3386, &
        1.3525,                                                  &
        1.3661,1.3796,1.3930,1.4062,1.4193,1.4322,1.4449,1.4576, &
        1.4701,1.4825,1.4949,1.5070,1.5191,1.5311,1.5430,1.5548, &
        1.5550,1.5550,1.5550,1.5550,1.5550,1.5550,1.5550,1.5550, &
        1.5550,                                                  &
        1.5550,1.5550,1.5550,1.5550,1.5550,1.5550,1.5550,1.5550, &
        1.5550,1.5550,1.5550,1.5549,1.5547,1.5543,1.5539,1.5532, &
        1.5525,1.5516,1.5506,1.5494,1.5481,1.5467,1.5452,1.5435, &
        1.5417,                                                  &
        1.5397,1.5377,1.5355,1.5332,1.5308,1.5282,1.5255,1.5228, &
        1.5199,1.5169,1.5137,1.5105,1.5072,1.5037,1.5002,1.4966, &
        1.4929,1.4890,1.4851,1.4811,1.4771,1.4729,1.4686,1.4643, &
        1.4599,                                                  &
        1.4555,1.4509,1.4463,1.4417,1.4370,1.4322,1.4274,1.4225, &
        1.4176,1.4126,1.4076,1.4025,1.3974,1.3923,1.3872,1.3820, &
        1.3768,1.3716,1.3663,1.3611,1.3558,1.3505,1.3452,1.3400, &
        1.3347,                                                  &
        1.3294,1.3241,1.3188,1.3135,1.3083,1.3030,1.2978,1.2926, &
        1.2874,1.2822,1.2771,1.2720,1.2669,1.2618,1.2568,1.2518, &
        1.2468,1.2419,1.2370,1.2322,1.2274,1.2227,1.2180,1.2133, &
        1.2087,                                                  &
        1.2041,1.1996,1.1952,1.1907,1.1864,1.1821,1.1778,1.1737, &
        1.1695,1.1654,1.1614,1.1575,1.1536,1.1497,1.1460,1.1422, &
        1.1386,1.1350,1.1314,1.1280,1.1246,1.1212,1.1179,1.1147, &
        1.1115,                                                  &
        1.1084,1.1053,1.1024,1.0994,1.0966,1.0938,1.0910,1.0883, &
        1.0857,1.0831,1.0806,1.0781,1.0757,1.0734,1.0711,1.0688, &
        1.0667,1.0645,1.0624,1.0604,1.0584,1.0565,1.0546,1.0528, &
        1.0510,                                                  &
        1.0493,1.0476,1.0460,1.0444,1.0429,1.0414,1.0399,1.0385, &
        1.0371,1.0358,1.0345,1.0332,1.0320,1.0308,1.0296,1.0285, &
        1.0275,1.0264,1.0254,1.0244,1.0235,1.0226,1.0217,1.0208, &
        1.0200,                                                  &
        1.0192,1.0184,1.0177,1.0170,1.0163,1.0156,1.0150,1.0143, &
        1.0137,1.0132,1.0126,1.0121,1.0116,1.0111,1.0106,1.0101, &
        1.0097,1.0092,1.0088,1.0084,1.0081,1.0077,1.0074,1.0070, &
        1.0067,                                                  &
        1.0064,1.0061,1.0058,1.0055,1.0053,1.0050,1.0048,1.0046, &
        1.0043,1.0041,1.0039,1.0037,1.0036,1.0034,1.0032,1.0030, &
        1.0029,1.0027,1.0026,1.0025,1.0023,1.0022,1.0021,1.0020, &
        1.0019,                                                  &
        1.0018,1.0017,1.0016,1.0015,1.0014,1.0014,1.0013,1.0012, &
        1.0011,1.0011,1.0010,1.0010,1.0009,1.0009,1.0008,1.0007, &
        1.0006,1.0005,1.0004,1.0003,1.0002,1.0001,1.0000,1.0000, &
        1.0000/

  data (xfac_rhu(i), i = -1, 61) / &
        0.7620,0.7840, &
        0.7820,0.7840,0.7620,0.7410,0.7970, &
        0.9140,0.9980,0.9830,0.9330,0.8850, &
        0.8420,0.8070,0.8000,0.8010,0.8100, &
        0.8090,0.8320,0.8180,0.7970,0.8240, &
        0.8640,0.8830,0.8830,0.8470,0.8380, &
        0.8660,0.9410,1.0400,1.0680,1.1410, &
        1.0800,1.0340,1.1550,1.0990,1.0270, &
        0.9500,0.8950,0.8150,0.7830,0.7700, &
        0.7000,0.7650,0.7750,0.8500,0.9000, &
        0.9050,0.9540,1.0200,1.0200,1.0250, &
        1.0200,1.1000,1.1250,1.1200,1.1110, &
        1.1370,1.1600,1.1490,1.1070,1.0640, &
        1.0450/

  hvrcnt = '$Revision: 32919 $'
  rhoave = (pave/p0)*(t0/tave)
  xkt = tave/radcn2 ![cm-1].
  amagat = (pave/p0)*(273./tave) ![amg = 44.615 mol m-3]
  wtot = wbroad
  do m = 1, nmol
    wtot = wtot + wk(m) ![cm-2].
  enddo
  x_vmr_h2o = wk(1)/wtot !unitless.
  x_vmr_o2  = wk(7)/wtot !unitless.
  x_vmr_n2  = 1. - x_vmr_h2o - x_vmr_o2 !unitless.
  wn2 = x_vmr_n2*wtot ![cm-2].

  fscdid(:) = 0
  xself = xcnt(1)
  xfrgn = xcnt(2)
  xco2c = xcnt(3)
  xo3cn = xcnt(4)
  xo2cn = xcnt(5)
  xn2cn = xcnt(6)
  xrayl = xcnt(7)

  !H2O continuum derivatives are computed w.r.t. ln(q)
  !dqh2o must be returned with the radiation field included
  if (icflg .gt. 0) then
    wn2 = 0.
    do j = 1, ipts
      cself(j) = 0.
      cfrgn_aj(j) = 0.
    enddo
    v1absc = v1abs
    v2absc = v2abs
    dvabsc = dvabs
    nptabsc = nptabs
  endif

  !CLOUD EFFECTIVE OPTICAL DEPTH  FROM "in_lblrtm_cld" file
  if (fscdid(4) .eq. 5) then
!   call cld_od(v1c, v2c, dvc, nptc, c_cld, layer, xkt)
!   !Radiation field
!   if (jrad .eq. 1) then
!     do j = 1, nptc
!       vj = v1c + real(j - 1)*dvc
!       c_cld(j) = c_cld(j)*RADFN(VJ,XKT)
!     enddo
!   endif
!   !Interpolate to total optical depth grid
!   call xint(v1c, v2c, dvc, c_cld, 1._real64, v1abs, dvabs, absrb, 1, nptabs)
  endif

  !WATER VAPOR
  h2o_fac = wk(1)/wtot !unitless.
  rself = h2o_fac*rhoave*1.e-20*xself
  rfrgn = (1. - h2o_fac)*rhoave*1.e-20*xfrgn
  rfrgn_aj = h2o_fac*rhoave*1.e-20*xfrgn

  !CORRECTION TO THE WATER VAPOR CONTINUUM    mt_ckd_2.4   Nov 2008
  !     The following modifications to the water vapor continuum arise
  !     from new analyses of ARM measurements in the microwave and far-IR
  !     regions. Analyses of measurements in the microwave are based
  !     primarily on the two-channel MWR (23.8 and 31.4 GHz) at SGP,
  !     with supporting evidence from 150 GHz MWRHF measurements during
  !     the COPS campaign and from 170 GHz GVRP measurements at SGP (V. H.
  !     Payne, E. J. Mlawer and S. A. Clough). Measurements in the far-IR
  !     were from the AERO_ext at the NSA site, in the time surrounding
  !     and including the RHUBC-I campaign (J. Delamere and S. A. Clough).
  !                             SELF
  if ((v2 .gt. -20.0) .and. (v1 .lt. 20000.) .and. xself .gt. 0.) then
    sh2ot0(:) = 0.
    sh2ot1(:) = 0.
    call sl296(v1c, v2c, dvc, nptc, sh2ot0, v1ss, v2ss, v1abs, v2abs, path)
    call sl260(v1c, v2c, dvc, nptc, sh2ot1, v1ss, v2ss, v1abs, v2abs, path)
    !Loop calculating self continuum optical depth
    tfac = (tave - t0)/(260. - t0)

    !MT_CKD_3.5  All previous IR corrections now included in stored coefficients
    !rather than correction functions.
    do j = 1, nptc
      vj = v1c + dvc*real(j - 1)
      sh2o = 0.
      if (sh2ot0(j) .gt. 0.) then
        sh2o = sh2ot0(j)*(sh2ot1(j)/sh2ot0(j))**tfac
      endif
      cself(j) = wk(1)*(sh2o*rself)
      v1h = v1c
      dvh = dvc
      npth = nptc
      csh2o(j) = 1.e-20*sh2o*xself ![cm3].
      if (jrad .eq. 1) cself(j) = cself(j)*radfn(vj, xkt)
    enddo

    !Interpolate to total optical depth grid
    call pre_xint(v1ss, v2ss, v1abs, dvabs, nptabs, ist, last)
    call xint(v1c, v2c, dvc, cself, 1._real64, v1abs, dvabs, absrb, ist, last)
  endif

!                             FOREIGN
  if ((v2 .gt. -20.0) .and. (v1 .lt. 20000.) .and. xfrgn .gt. 0.) then
    fh2o(:) = 0.
    !correction to foreign continuum   mt_ckd_2.4  nov 2008   sac
    f0 = 0.06
    v0f1 = 255.67
    hwsq1 = 240.**2
    beta1 = 57.83
    c_1 = -0.42
    n_1 = 8

    c_2 = 0.3
    beta2 = 630.
    n_2 = 8
    !mt_ckd_2.8    March 2016     Mlawer and Alvarado
    !Extensive changes to the foreign continuum were made for
    !mt_ckd_2.8.  Based on measurements by Baranov and Lafferty
    !(2012) from 850-1160 cm-1 and Mondelain et al.(2014) at
    !4255 cm-1, a revised MT_CKD foreign continuum formulation
    !was derived and has been implemented in window regions >
    !4000 cm-1 (blended with previous coefficients in transition
    !regions between windows and bands). For 1800-3000 cm-1,
    !this formulation guided the spectral shape of the foreign
    !coefficients, but the values were reduced to obtain
    !agreement with measurements in this window by Baranov and
    !Lafferty (2012) and IASI measurements from 1900-2150 cm-1.
    !Coefficients in this region were derived simultaneously
    !with N2-H2O CIA coefficients and water vapor self continuum
    !coefficents.
    call frn296(v1c, v2c, dvc, nptc, fh2o, v1ss, v2ss, v1abs, v2abs, path)
    do j = 1, nptc
      vj = v1c + dvc*real(j - 1)
      if (vj .le. 600.) then
        jfac = (vj + 10.)/10. + 0.00001
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
      fh2o(j) = fh2o(j)*fscal
      c_f = wk(1)*fh2o(j)
      cfh2o(j) = 1.e-20*fh2o(j)*xfrgn
      if (jrad .eq. 1) c_f = c_f*radfn(vj, xkt)
      c(j) = c_f*rfrgn
      cfrgn_aj(j) = c_f*rfrgn_aj
    enddo
    call pre_xint(v1ss, v2ss, v1abs, dvabs, nptabs, ist, last)
    call xint(v1c, v2c, dvc, c, 1._real64, v1abs, dvabs, absrb, ist, last)
    if (icflg .eq. 1) then
      do j = 1, nptc
        c(j) = cself(j) - cfrgn_aj(j)
      enddo
      call pre_xint(v1ss, v2ss, v1abs, dvabs, nptabs, ist, last)
      call xint(v1c, v2c, dvc, c, 1._real64, v1abs, dvabs, absrb, ist, last)
    endif
  endif

  !CARBON DIOXIDE
  if ((v2 .gt. -20.0) .and. (v1 .lt. 10000.) .and. xco2c .gt. 0) then
    fco2(:) = 0.
    wco2 = wk(2)*rhoave*1.0e-20*xco2c
    call frnco2(v1c, v2c, dvc, nptc, fco2, tave, v1ss, v2ss, v1abs, v2abs, path)
    do j = 1, nptc
      vj = v1c + dvc*real(j - 1)
      !**mt_ckd_2.0      11 July 2007    sac
      !             This continuum differs from mt_ck_1.3 in that an entirely
      !             new co2 continuum has been developed based on the line
      !             coupling parameters from Hartmann's group as distributed
      !             with hitran.  This continuum must be used with lblrtm_v11
      !             and spectral line parameters including Hartmann's line
      !             parameters for co2.
      !             Based on recent validation studies, a scaling of the
      !             continuum for v3 is required to achieve an acceptable
      !             result at 2385 cm-1, the 'bandhead' of v3.
      !             Clough et al., presentation at EGU 2007
      !*** mt_ckd_2.5  Adjustment to the original scaling made.
      !                (temperature dependence of continuum also added)
      cfac = 1.
      if (vj .ge. 2000. .and. vj .le. 2998.) then
        jfac = (vj - 1998.)/2. + 0.00001
        cfac = xfacco2(jfac)
      endif
      fco2(j) = cfac*fco2(j)
      c(j) = fco2(j)*wco2
      if (jrad.eq.1) c(j) = c(j)*radfn(vj, xkt)
    enddo
    call pre_xint(v1ss, v2ss, v1abs, dvabs, nptabs, ist, last)
    call xint(v1c, v2c, dvc, c, 1._real64, v1abs, dvabs, absrb, ist, last)
  endif

  !DIFFUSE OZONE
  if (v2 .gt. 8920. .and. v1 .le. 24665. .and. xo3cn .gt. 0.) then
    cch0(:) = 0.
    cch1(:) = 0.
    cch2(:) = 0.

    wo3 = wk(3)*1.0e-20*xo3cn
    call xo3chp(v1c, v2c, dvc, npto3, cch0, cch1, cch2, v1ss, v2ss, v1abs, v2abs, path)
    dt = tave - 273.15

    do j = 1, npto3
      cch0(j) = (cch0(j) + (cch1(j) + cch2(j)*dt)*dt)*wo3
      vj = v1c + dvc*real(j - 1)
      if (jrad .eq. 1) cch0(j) = cch0(j)*radfn(vj, xkt)
    enddo
    call pre_xint(v1ss, v2ss, v1abs, dvabs, nptabs, ist, last)
    call xint(v1c, v2c, dvc, cch0, 1._real64, v1abs, dvabs, absrb, ist, last)
  endif

  if (v2 .gt. 27370. .and. v1 .lt. 40800. .and.xo3cn .gt. 0.) then
    c0(:) = 0.
    ct1(:) = 0.
    ct2(:) = 0.

    wo3 = wk(3)*1.e-20*xo3cn
    tc = tave - 273.15
    call o3hht0(v1c, v2c, dvc, npto3, c0, v1ss, v2ss, v1abs, v2abs, path)
    call o3hht1(v1t1, v2t1, dvt1, npt1, ct1, v1abs, v2abs, path)
    call o3hht2(v1t2, v2t2, dvt2, npt2, ct2, v1abs, v2abs, path)

    do j = 1, npto3
      c(j) = c0(j)*wo3
      vj = v1c + dvc*real(j - 1)
      if (jrad .eq. 1) c(j) = c(j)*radfn(vj, xkt)
      c(j) = c(j)*(1. + ct1(j)*tc + ct2(j)*tc*tc)
    enddo

    !Save non-Hartley Huggins optical depth contribution to
    !prevent double counting for wavenumber region beyond
    !40800 cm-1.
    if ((vj .gt. 40815.) .and. (v2 .gt. 40800) .and. xo3cn .gt. 0.) then
      i_fix = (40800. - v1abs)/dvabs + 1.001
      do i = i_fix, nptabs
        absbsv(i) = absrb(i)
      enddo
    endif

    !Combine Hartley Huggins with previous optical depths
    call pre_xint(v1ss, v2ss, v1abs, dvabs, nptabs, ist, last)
    call xint(v1c, v2c, dvc, c, 1._real64, v1abs, dvabs, absrb, ist, last)

    !If V2 > 40800 cm-1, replace points with previously
    !saved values (non-Hartley Huggins contribution)
    if ((vj .gt. 40815.) .and. (v2 .gt. 40800) .and. xo3cn .gt. 0.) then
      do i = i_fix, nptabs
        absrb(i) = absbsv(i)
      enddo
    endif
  endif

  if (v2 .gt. 40800. .and. v1 .lt. 54000. .and. xo3cn .gt. 0.) then
    c0(:) = 0.

    wo3 = wk(3)*xo3cn
    call o3hhuv(v1c, v2c, dvc, npto3, c0, v1ss, v2ss, v1abs, v2abs, path)

    do j = 1, npto3
      c(j) = c0(j)*wo3
      vj = v1c + dvc*real(j - 1)
      if (jrad .eq. 1) c(j) = c(j)*radfn(vj, xkt)
    enddo

    !Save non-Hartley Huggins UV optical depth contribution to
    !prevent double counting for wavenumber region before
    !40800 cm-1.
    if (v1 .lt. 40800) then
      i_fix = (40800. - v1abs)/dvabs + 1.001
      do i = 1, i_fix - 1
        absbsv(i) = absrb(i)
      enddo
    endif

    !Combine UV Hartley Huggins with previous optical depths
    call pre_xint(v1ss, v2ss, v1abs, dvabs, nptabs, ist, last)
    call xint(v1c, v2c, dvc, c, 1._real64, v1abs, dvabs, absrb, ist, last)

    !If V1 < 40800 cm-1, replace points with previously
    !saved values (non-Hartley Huggins UV contribution)
    if (v1 .lt. 40800) then
      do i = 1, i_fix - 1
        absrb(i) = absbsv(i)
      enddo
    endif
  endif

  !O2 OXYGEN COLLISION INDUCED FUNDAMENTAL
  !version_1 of the Oxygen Collision Induced Fundamental
  !     F. Thibault, V. Menoux, R. Le Doucen, L. Rosenman, J.-M. Hartmann,
  !        and Ch. Boulet,
  !        Infrared collision-induced absorption by O2 near 6.4 microns
  !        for atmospheric applications: measurements and emprirical
  !        modeling, Appl. Optics, 35, 5911-5917, (1996).
  if ((v2 .gt. 1340.0) .and. (v1 .lt. 1850.) .and. xo2cn .gt. 0.) then
    c0(:) = 0.
    tau_fac = xo2cn*wk(7)*1.e-20*amagat

    !Wk(7) is the oxygen column amount in units of molec/cm2
    !amagat is in units of amagats (air)

    !The temperature correction is done in the subroutine o2_ver_
    call o2_ver_1(v1c, v2c, dvc, nptc, c0, tave, v1ss, v2ss, v1abs, v2abs, path)

    !c0 are the oxygen absorption coefficients at temperature
    !tave
    !   - these absorption coefficients are in units of
    !     [(cm^2/molec) 10^20)]/(cm-1  amagat)
    !   - cm-1 in the denominator arises through the removal
    !     of the radiation field
    !   - for this case, an amagat is interpreted as one
    !     loschmidt of air (273K)
    do j = 1, nptc
      vj = v1c + dvc*real(j - 1)
      c(j) = tau_fac*c0(j)
      if (jrad .eq. 1) c(j) = c(j)*radfn(vj, xkt)
    enddo
    call pre_xint(v1ss, v2ss, v1abs, dvabs, nptabs, ist, last)
    call xint(v1c, v2c, dvc, c, 1._real64, v1abs, dvabs, absrb, ist, last)
  endif

  !O2 Collision Induced
  !O2 continuum formulated by Mate et al. over the spectral region
  !7550-8486 cm-1:  "Absolute Intensities for the O2 1.27 micron
  !continuum absorption", B. Mate, C. Lugez, G.T. Fraser, and
  !W.J. Lafferty, J. Geophys. Res., 104, 30,585-30,590, 1999.

  !Units of these coefficients are 1 / (amagat_O2*amagat_air)

  !Also, refer to the paper "Observed  Atmospheric
  !Collision Induced Absorption in Near Infrared Oxygen Bands",
  !Mlawer, Clough, Brown, Stephen, Landry, Goldman, & Murcray,
  !Journal of Geophysical Research (1998).
  if ((v2 .gt. 7536.0) .and. (v1 .lt. 8500.) .and. xo2cn .gt. 0.) then
    c0(:) = 0.
    a_o2  = 1./0.446
    a_n2  = 0.3/0.446
    a_h2o = 1.
    tau_fac = xo2cn*(wk(7)/xlosmt)*amagat* &
              (a_o2*x_vmr_o2 + a_n2*x_vmr_n2 + a_h2o*x_vmr_h2o)
    call o2inf1(v1c, v2c, dvc, nptc, c0, v1ss, v2ss, v1abs, v2abs, path)
    do j = 1, nptc
      c(j) = tau_fac*c0(j)
      vj = v1c + dvc*real(j - 1)
      if (jrad .eq. 1) c(j) = c(j)*radfn(vj, xkt)
    enddo
    call pre_xint(v1ss, v2ss, v1abs, dvabs, nptabs, ist, last)
    call xint(v1c, v2c, dvc, c, 1._real64, v1abs, dvabs, absrb, ist, last)
  endif

  !O2 continuum formulated by Mlawer et al. over the spectral
  !region 9100-11000 cm-1. Refer to the paper "Observed
  !Atmospheric Collision Induced Absorption in Near Infrared
  !Oxygen Bands", Mlawer, Clough, Brown, Stephen, Landry, Goldman,
  !& Murcray, Journal of Geophysical Research (1998).
  if ((v2 .gt. 9100.0) .and. (v1 .lt. 11000.) .and. xo2cn .gt. 0.) then
    c0(:) = 0.
    call o2inf2(v1c, v2c, dvc, nptc, c0, v1ss, v2ss, v1abs, v2abs)
    wo2 = xo2cn*wk(7)*1.e-20*rhoave
    adjwo2 = (wk(7)/wtot)*(1./0.209)*wo2
    do j = 1, nptc
      c(j) = c0(j)*adjwo2
      vj = v1c + dvc*real(j - 1)
      if (jrad .eq. 1) c(j) = c(j)*radfn(vj, xkt)
    enddo
    call pre_xint(v1ss, v2ss, v1abs, dvabs, nptabs, ist, last)
    call xint(v1c, v2c, dvc, c, 1._real64, v1abs, dvabs, absrb, ist, last)
  endif

  !O2 A-band continuum formulated by Mlawer based on solar FTS measurements.
  if ((v2 .gt. 12961.5) .and. (v1 .lt. 13221.5) .and. xo2cn .gt. 0.) then
    c0(:) = 0.
    tau_fac = xo2cn*(wk(7)/xlosmt)*amagat
    call o2inf3(v1c, v2c, dvc, nptc, c0, v1ss, v2ss, v1abs, v2abs, path)
    do j = 1, nptc
      c(j) = tau_fac*c0(j)
      vj = v1c + dvc*real(j - 1)
      if (jrad .eq. 1) c(j) = c(j)*radfn(vj, xkt)
    enddo
    call pre_xint(v1ss, v2ss, v1abs, dvabs, nptabs, ist, last)
    call xint(v1c, v2c, dvc, c, 1._real64, v1abs, dvabs, absrb, ist, last)
  endif

  !O2 continuum formulated by Greenblatt et al. over the spectral
  !region 8797-29870 cm-1:  "Absorption Coefficients of Oxygen
  !Between 330 and 1140 nm, G.D. Greenblatt, J.J. Orlando, J.B.
  !Burkholder, and A.R. Ravishabkara,  J. Geophys. Res., 95,
  !18577-18582, 1990.

  !The units conversion to (cm^2/molec)/atm(o2)  has been done in
  !subroutine o2_vis.
  if ((v2 .gt. 15000.0) .and. (v1 .lt. 29870.) .and. xo2cn .gt. 0.) then
    c0(:) = 0.
    wo2 = wk(7)*1.e-20*((pave/1013.)*(273./tave))*xo2cn
    chio2 =  wk(7)/wtot
    adjwo2 = chio2*wo2
    call o2_vis(v1c, v2c, dvc, nptc, c0, v1ss, v2ss, v1abs, v2abs, path)
    do j = 1, nptc
      c(j) = c0(j)*adjwo2
      vj = v1c + dvc*real(j - 1)
      if (jrad .eq. 1) c(j) = c(j)*radfn(vj, xkt)
    enddo
    call pre_xint(v1ss, v2ss, v1abs, dvabs, nptabs, ist, last)
    call xint(v1c, v2c, dvc, c, 1._real64, v1abs, dvabs, absrb, ist, last)
  endif

  if (v2 .gt. 36000.0 .and. xo2cn .gt. 0.) then
    c0(:) = 0.
    wo2 = wk(7)*1.e-20*xo2cn
    call o2herz(v1c, v2c, dvc, nptc, c0, tave, pave, v1ss, v2ss, v1abs, v2abs)
    do j = 1, nptc
      c(j) = c0(j)*wo2
      vj = v1c + dvc*real(j - 1)
      if (jrad .eq. 1) c(j) = c(j)*radfn(vj, xkt)
    enddo
    call pre_xint(v1ss, v2ss, v1abs, dvabs, nptabs, ist, last)
    call xint(v1c, v2c, dvc, c, 1._real64, v1abs, dvabs, absrb, ist, last)
  endif

  !O2 in far-UV (large portion is referred to as Schumann-Runge continuum)
  !Only calculate if V2 > 56740. cm-1.
  if (v2 .gt. 56740.0 .and. xo2cn .gt. 0.) then
    c0(:) = 0.
    wo2 = wk(7)*1.e-20*xo2cn
    call o2fuv(v1c, v2c, dvc, nptc, c0, v1ss, v2ss, v1abs, v2abs, path)
    do j = 1, nptc
      vj = v1c + dvc*real(j - 1)
      c(j) = c0(j)*wo2
      if (jrad .eq. 1) c(j) = c(j)*radfn(vj, xkt)
    enddo
    call pre_xint(v1ss, v2ss, v1abs, dvabs, nptabs, ist, last)
    call xint(v1c, v2c, dvc, c, 1._real64, v1abs, dvabs, absrb, ist, last)
  endif

  !NITROGEN CONTINUA
  !NITROGEN COLLISION INDUCED PURE ROTATION BAND
  !        Model used:
  !         Borysow, A, and L. Frommhold, "Collision-induced
  !            rototranslational absorption spectra of N2-N2
  !            pairs for temperatures from 50 to 300 K", The
  !            Astrophysical Journal, 311, 1043-1057, 1986.

  !     Uodated 2004/09/22 based on:

  !      Boissoles, J., C. Boulet, R.H. Tipping, A. Brown and Q. Ma,
  !         Theoretical Calculations of the Translation-Rotation
  !         Collision-Induced Absorption in N2-N2, O2-O2 and N2-O2 Pairs,
  !         J.Quant. Spec. Rad. Transfer, 82,505 (2003).

  !         The temperature dependence between the two reference
  !         temperatures has been assumed the same as that for the
  !         original continuum.

  !        THIS NITROGEN CONTINUUM IS IN UNITS OF 1./(CM AMAGAT^2)

  !        Only calculate if V2 > -10. cm-1 and V1 <  350. cm-1
  if ((v2 .gt. -10.0) .and. (v1 .lt. 350.) .and. xn2cn .gt. 0.) then
    c0(:) = 0.
    c1(:) = 0.
    !The following puts WXN2 units in 1./(CM AMAGAT)
    !     c1(j) represents the relative broadening efficiency of o2
    !     a_h2o represents the relative broadening efficiency of h2o
    a_h2o = 1.

    !correct formulation for consistency with LBLRTM (per molec/cm^2)
    tau_fac = xn2cn*(wn2/xlosmt)*amagat
    call xn2_r(v1c, v2c, dvc, nptc, c0, c1, tave, v1ss, v2ss, v1abs, v2abs, path)
    !c1 is  ~ the ratio of alpha(n2-o2)/alpha(n2-n2)
    !Eq's 7 and 8 in the Boissoles paper.

    do j = 1, nptc
      vj = v1c + dvc*real(j - 1)
      c(j) = tau_fac*c0(j)*(x_vmr_n2 + c1(j)*x_vmr_o2 + a_h2o*x_vmr_h2o)
      if (jrad .eq. 1) c(j) = c(j)*radfn(vj, xkt)
    enddo
    call pre_xint(v1ss, v2ss, v1abs, dvabs, nptabs, ist, last)
    call xint(v1c, v2c, dvc, c, 1._real64, v1abs, dvabs, absrb, ist, last)
  endif

  !NITROGEN COLLISION INDUCED FUNDAMENTAL

  !        version_1 of the Nitrogen Collision Induced Fundamental

  !        Lafferty, W.J., A.M. Solodov,A. Weber, W.B. Olson and
  !        J._M. Hartmann, Infrared collision-induced absorption by
  !        N2 near 4.3 microns for atmospheric applications:
  !        Measurements and emprirical modeling, Appl. Optics, 35,
  !        5911-5917, (1996).

  !        Only calculate if V2 > 2001.77 cm-1 and V1 < 2897.59 cm-1.
  !        Bounds for original implmentation were from Lafferty (2085-
  !        2670 cm-1).  N2-N2 coefficients were analytically extended
  !        to 2000-2900 cm-1 to allow implementation of N2-H2O, which
  !        has been assumed to have a wider spectral shape to allow
  !        agreement with  Baranov and Lafferty (2012).
  if ((v2 .gt. 2001.77) .and. (v1 .lt. 2897.59) .and. xn2cn .gt. 0.) then
    cn0 = 0.
    cn1 = 0.
    cn2 = 0.

    !           The absorption coefficients from the Lafferty et al.
    !           reference are for pure nitrogen (absorber and broadener).

    !     correct formulation for consistency with LBLRTM (per molec/cm^2)
    tau_fac = xn2cn*(wn2/xlosmt)*amagat

    !           Wn2 is in units of molec/cm2
    !           amagat is in units of amagats (air)

    !           The temperature correction of the absorption coefficient and the
    !           adjustments for relative broadening efficiency are done in
    !           subroutine n2_ver_1:
    call n2_ver_1 (v1c, v2c, dvc, nptc, cn0, cn1, cn2, tave, v1ss, v2ss, v1abs, v2abs, path)

    !           cn0 are the nitrogen absorption coefficients at
    !           temperature tave
    !              - these absorption coefficients are in units of
    !                   [(cm^2/molec) 10^20)]/(cm-1  amagat)
    !              - cm-1 in the denominator arises through the removal
    !                   of the radiation field
    !              - for this case, an amagat is interpreted as one
    !                   loschmidt of air (273K)
    !           cn1 are N2 absorption coefficents with collision partner O2
    !           cn2 are N2 absorption coefficents with collision partner H2O
    do j = 1, nptc
      vj = v1c + dvc*real(j - 1)
      c(j) = tau_fac*(x_vmr_n2*cn0(j) + x_vmr_o2*cn1(j) + x_vmr_h2o*cn2(j))
      if (jrad .eq. 1) c(j) = c(j)*radfn(vj, xkt)
    enddo
    call pre_xint(v1ss, v2ss, v1abs, dvabs, nptabs, ist, last)
    call xint(v1c, v2c, dvc, c, 1._real64, v1abs, dvabs, absrb, ist, last)
  endif

  !NITROGEN COLLISION INDUCED FIRST OVERTONE

  !        version_1 of the Nitrogen Collision Induced First Overtone

  !        Shapiro and Gush (1966) modified by Mlawer and Gombos (2015)
  !        based on comparisons with measurements from SGP solar FTS
  !        (TCCON network).
  if ((v2 .gt. 4340.0) .and. (v1 .lt. 4910.) .and. xn2cn .gt. 0.) then
    c0(:) = 0.

    !        All species are assumed to have the same broadening efficiency.
    !           a_o2  represents the relative broadening efficiency of o2
    !           a_h2o represents the relative broadening efficiency of h2o
    a_o2  = 1.
    a_h2o = 1.

    !     correct formulation for consistency with LBLRTM (per molec/cm^2)
    tau_fac = xn2cn*(wn2/xlosmt)*amagat*(x_vmr_n2 + a_o2*x_vmr_o2 + a_h2o*x_vmr_h2o)

    !           Wn2 is in units of molec/cm2
    !           amagat is in units of amagats (air)

    !           The absorption coefficients are assumed to have no
    !           temperature dependence.
    call n2_overtone1(v1c, v2c, dvc, nptc, c0, v1ss, v2ss, v1abs, v2abs, path)

    !           c0 are the nitrogen absorption coefficients
    !              - these absorption coefficients are in units of
    !                   [(cm^2/molec) 10^20)]/(cm-1  amagat)
    !              - cm-1 in the denominator arises through the removal
    !                   of the radiation field
    !              - for this case, an amagat is interpreted as one
    !                   loschmidt of air (273K)
    do j = 1, nptc
      vj = v1c + dvc*real(j - 1)
      c(j) = tau_fac*c0(j)
      if (jrad .eq. 1) c(j) = c(j)*radfn(vj, xkt)
    enddo
    call pre_xint(v1ss, v2ss, v1abs, dvabs, nptabs, ist, last)
    call xint(v1c, v2c, dvc, c, 1._real64, v1abs, dvabs, absrb, ist, last)
  endif

  !Rayleigh Scattering calculation
  !     (sac)
  !     The effects of Rayleigh scattering are also included in module
  !     lbllow.f with the aerosol/cloud properties.  In the case that
  !     that lbllow.f is selected (fscdid(4) .ne. 0), the decision has been
  !     made to include this effect with the scattering processes,
  !     even though it is molecular in nature.  Otherwise the effects
  !     of Rayleigh scattering are included here.

  !     The formulation, adopted from MODTRAN_3.5 (using approximation
  !     of Shettle et al., (Appl Opt, 2873-4, 1980) with depolarization
  !     = 0.0279, output in km-1 for T=273K & P=1 ATM) is as follows:

  !     The rayleigh extinction coefficient (scattering out of the direct
  !     beam), ray_ext, can be defined as

  !         ray_ext = (vrayleigh**4/(9.38076E18-1.08426E09*vrayleigh**2))
  !     *        *wmol_tot*conv_cm2mol

  !     where vrayleigh is the wavenumber value, wmol_tot is the total
  !     molecular amount in the layer, and conv_cm2mol is the conversion
  !     factor derived by multiplying air density (2.68675E19 mol/cm3)
  !     at 273 K with the number of km per cm (1.e-5 km/cm).

  !     The total layer amount of all molecules is calculated above as
  !     WTOT. For numerical purposes a factor of 1.e-20  has been
  !     included in conv_cm2mol and the same factor has been included
  !     in the air density in the denominator as well.
  !     In addition, a factor of 10,000 (cm-1) has been
  !     divided out of vrayleigh. Finally, the radiation field is
  !     excluded, so xvrayleigh**4 is replaced by xvrayleigh**3. When
  !     JRAD=1, the radiation field is put in by multiplying the
  !     absorption by xvrayleigh.

  !     Rayleigh scattering in the direct beam is only calculated for
  !     model runs > 820 cm-1.
  if ((fscdid(4) .eq. 0 .or. fscdid(4) .eq. 5) .and. v2 .ge. 820. &
      .and. xrayl .gt. 0.) then
    conv_cm2mol = xrayl*1.e-20/(2.68675e-1*1.e5)
    do i = 1, nptabs
      vrayleigh = v1abs + (i - 1)*dvabs
      xvrayleigh = vrayleigh/1.e4
      ray_ext = (xvrayleigh**3/(9.38076e2 - 10.8426*xvrayleigh**2))*(wtot*conv_cm2mol)
      if (jrad .eq. 1) then
        ray_ext = ray_ext*xvrayleigh
      elseif (jrad .eq. 0) then
        ray_ext = ray_ext*xvrayleigh/radfn(vrayleigh,xkt)
      endif
      absrb(i) = absrb(i) + ray_ext
    enddo
  endif
end subroutine contnm


!subroutine cld_od(V1C,V2C,DVC,NPTC,C,layer,xkt)
!   COMMON /ABSORB/ V1ABS,V2ABS,DVABS,NPTABS,ABSRB(n_absrb)

!   COMMON /IFIL/ IRD,IPR,IPU,NOPR,NFHDRF,NPHDRF,NFHDRL,NPHDRL,       &
!   &              NLNGTH,KFILE,KPANEL,LINFIL,NFILE,IAFIL,IEXFIL,      &
!   &              NLTEFL,LNFIL4,LNGTH4

!   parameter (n_lyr=200,n_cld=500)

!   COMMON /cld_rd/ n_freq, n_align, v_cloud_freq(n_cld),             &
!   &                                      cloudodlayer(n_lyr,n_cld)
!   DIMENSION C(*)

!   logical EX
!   character*55 in_cld_file
!   dimension i_layer(n_lyr), pres_layer_dum(n_lyr), v_cntnm(n_lyr)

!   data in_cld_file /'in_lblrtm_cld'/
!   data dvs /5./

!     Read in TES cloud effective optical depth file
!  if (layer .eq. 1) then

!     open (35,FILE=in_cld_file,STATUS='OLD')
!     read (35,*) n_freq
!     read (35,*) (v_cloud_freq(j),j=1,n_freq)

!     write (ipr,*)
!     write (ipr,*)                                                  &
!     &           '** iaersl=5; Cloud Information from "in_cld_file" **'
!     write (ipr,'(" n_freq = ",i5)') n_freq
!     write (ipr,'(5x,10f10.4)') (v_cloud_freq(j),j=1,n_freq)

!     read (35,*) n_layer

!     write (ipr,'(" n_layer = ",i5)') n_layer

!     do l =1,n_layer
!        read (35,*) i_layer(l), pres_layer_dum(l)
!        read (35,*) (cloudodlayer(l,j),j=1,n_freq)
!        write (ipr,'(i5,f12.5)') i_layer(l), pres_layer_dum(l)
!        write (ipr,'(5x,10f10.4)') (cloudodlayer(l,j),j=1,n_freq)
!     enddo
!     close (35)
!  endif

!        Generated output continuum grid
!  DVC = DVS
!  V1C = V1ABS-10.
!  V2C = V2ABS+10.
!  NPTC = ((V2C - V1C)/DVC) + 1
!  do J = 1, NPTC
!     v_cntnm(j) = V1C+DVC* REAL(J-1)
!  enddo

!        Linearly interpolate TES cloud effective OD onto continuum grid
!  ilo = 1
!  do j=1, NPTC
!     IF (v_cntnm(j).LE.v_cloud_freq(1))      THEN
!        C(j) = cloudodlayer(layer,1)
!        GO TO 10
!     ELSE IF (v_cntnm(j).GT.v_cloud_freq(n_freq)) THEN
!        C(j) = cloudodlayer(layer,n_freq)
!        GO TO 10
!     END IF

!     do i=ilo, n_freq
!        IF (v_cntnm(j).LE.v_cloud_freq(i))  THEN
!           v_m=(cloudodlayer(layer,i)-cloudodlayer(layer,i-1))/ (         &
!              v_cloud_freq(i)-v_cloud_freq(i-1))
!           C(j) = cloudodlayer(layer,i-1)+ (v_cntnm(j)-v_cloud_freq(i-1))*&
!              v_m
!           ilo = i-1
!           GO TO 10
!        END IF
!     enddo

!10    continue

!     if (v_cntnm(j).eq.0.) then
!        C(j) = 0.
!     else
!        C(j) = C(j)/RADFN(v_cntnm(j),XKT)
!     endif

!     if (ilo.lt.1) ilo = 1
!   enddo
!end subroutine cld_od


end module continuum
