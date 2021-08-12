module continuum
use, intrinsic :: iso_fortran_env, only: real64
use utils, only: pre_xint, radfn, SpectralGrid, xint, t0, t_273, p0, loschmidt, radcn2
use co2_continuum, only: carbon_dioxide_continuum
use h2o_continuum, only: water_vapor_self_continuum, water_vapor_foreign_continuum
use n2_continuum, only: nitrogen_continuum, nitrogen_fundamental_continuum, &
                        nitrogen_overtone_continuum
use o2_continuum, only: oxygen_fundamental_continuum, oxygen_nir1_continuum, &
                        oxygen_nir2_continuum, oxygen_nir3_continuum, &
                        oxygen_visible_continuum, oxygen_herzberg_continuum, &
                        oxygen_uv_continuum
use o3_continuum, only: ozone_cw_continuum, ozone_hh_continuum, ozone_uv_continuum
implicit none


public :: contnm


integer, parameter, public :: h2o = 1
integer, parameter, public :: co2 = 2
integer, parameter, public :: o3 = 3
integer, parameter, public :: o2 = 4
integer, parameter, public :: n2 = 5
integer, parameter, public :: ar = 6


contains


!contains the continuum data which is interpolated into the array absrb.
subroutine contnm(pave, tave, vmr, path_length, grid, icflg, jrad, xcnt, &
                  h2o_grid, absrb, path, csh2o, cfh2o)

  real(kind=real64), intent(in) :: pave !< Pressure [mb].
  real(kind=real64), intent(in) :: tave !< Temperature [K].
  real(kind=real64), dimension(:), intent(in) :: vmr !< Volume mixing ratio in dry air [mol mol-1].
  real(kind=real64), intent(in) :: path_length !< Path length [cm].
  type(SpectralGrid), intent(in) :: grid
  integer, intent(in) :: icflg, jrad
  real(kind=real64), dimension(7), intent(in) :: xcnt
  type(SpectralGrid), intent(out) :: h2o_grid
  real(kind=real64), dimension(:), intent(inout) :: absrb
  character(len=*), intent(in) :: path
  real(kind=real64), dimension(:), intent(inout) :: csh2o, cfh2o

  integer :: j, nptabsc
  integer, parameter :: ipts = 5050
  integer, parameter :: n_absrb = 5050
  integer, dimension(17) :: fscdid
  real(kind=real64) :: wn2, x_vmr_h2o, x_vmr_o2, x_vmr_n2, &
                       v1absc, v2absc, dvabsc, xself, xfrgn, xco2c, xo3cn, xo2cn, &
                       xn2cn, xrayl
  real(kind=real64), dimension(ipts) :: cself, cfrgn_aj

  real(kind=real64) :: air_density, dry_air_density, total_density
  real(kind=real64), dimension(size(vmr)) :: density

  air_density = loschmidt*(pave/p0)*(t_273/tave) ![cm-3].
  dry_air_density = air_density*(1. - vmr(h2o)) ![cm-3].
  density(:) = vmr(:)*dry_air_density ![cm-3].
  total_density = sum(density) ![cm-3].
  x_vmr_h2o = density(h2o)/total_density !unitless (different from vmr(h2o)?).
  x_vmr_o2  = density(o2)/total_density !unitless (different from vmr(o2)?).
  x_vmr_n2  = 1. - x_vmr_h2o - x_vmr_o2 !unitless (different from vmr(n2)?).
  wn2 = x_vmr_n2*total_density ![cm-3].

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
    v1absc = grid%x0
    v2absc = grid%xn
    dvabsc = grid%dx
    nptabsc = grid%n
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

  !Water vapor
  if (xself .gt. 0.) then
    call water_vapor_self_continuum(n_absrb, jrad, density(h2o), total_density, &
                                    grid, tave, pave, path, h2o_grid, &
                                    csh2o, absrb, cself)
  endif
  if (xfrgn .gt. 0) then
    call water_vapor_foreign_continuum(n_absrb, icflg, jrad, density(h2o), &
                                       total_density, grid, tave, pave, &
                                       cself, path, cfh2o, absrb)
  endif

  !Carbon dioxide
  if (xco2c .gt. 0.) then
    call carbon_dioxide_continuum(n_absrb, jrad, density(co2), grid, tave, pave, &
                                  path, absrb)
  endif

  !DIFFUSE OZONE
  if (xo3cn .gt. 0.) then
    call ozone_cw_continuum(jrad, density(o3), grid, tave, path, absrb)
    call ozone_hh_continuum(n_absrb, jrad, density(o3), grid, &
                            tave, path, absrb)
    call ozone_uv_continuum(n_absrb, jrad, density(o3), grid, &
                            tave, path, absrb)
  endif

  !OXYGEN
  if (xo2cn .gt. 0.) then
    call oxygen_fundamental_continuum(n_absrb, jrad, density(o2), grid, &
                                      tave, pave, path, absrb)
    call oxygen_nir1_continuum(n_absrb, jrad, density(o2), grid, &
                               tave, pave, x_vmr_o2, x_vmr_n2, x_vmr_h2o, &
                               path, absrb)
    call oxygen_nir2_continuum(n_absrb, jrad, density(o2), total_density, &
                               grid, tave, pave, absrb)
    call oxygen_nir3_continuum(n_absrb, jrad, density(o2), grid, &
                               tave, pave, path, absrb)
    call oxygen_visible_continuum(n_absrb, jrad, density(o2), total_density, &
                                  grid, tave, pave, path, absrb)
    call oxygen_herzberg_continuum(n_absrb, jrad, density(o2), grid, &
                                   tave, pave, absrb)
    call oxygen_uv_continuum(n_absrb, jrad, density(o2), grid, &
                             tave, path, absrb)
  endif

  !NITROGEN CONTINUA
  if (xn2cn .gt. 0.) then
    call nitrogen_continuum(n_absrb, jrad, wn2, grid, &
                            tave, pave, x_vmr_n2, x_vmr_o2, x_vmr_h2o, path, absrb)
    call nitrogen_fundamental_continuum(jrad, wn2, grid, &
                                        tave, pave, x_vmr_n2, x_vmr_o2, x_vmr_h2o, &
                                        path, absrb)
    call nitrogen_overtone_continuum(n_absrb, jrad, wn2, grid, &
                                     tave, pave, x_vmr_n2, x_vmr_o2, x_vmr_h2o, &
                                     path, absrb)
  endif

  !Rayleigh Scattering calculation
  if (xrayl .gt. 0.) then
    call rayleigh_scattering(fscdid(4), jrad, total_density, grid, &
                             tave, absrb)
  endif
  absrb(:) = absrb(:)*path_length
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


subroutine rayleigh_scattering(fscdid, jrad, wtot, grid, tave, absrb)

  integer, intent(in) :: fscdid, jrad
  real(kind=real64), intent(in) :: wtot, tave
  type(SpectralGrid), intent(in) :: grid
  real(kind=real64), dimension(:), intent(inout) :: absrb

  integer :: i
  real(kind=real64) :: conv_cm2mol, vrayleigh, xvrayleigh, ray_ext

  if ((fscdid .eq. 0 .or. fscdid .eq. 5) .and. grid%xn .ge. 820._real64) then
    conv_cm2mol = 1.e-20_real64/(2.68675e-1_real64*1.e5_real64)
    do i = 1, grid%n
      vrayleigh = grid%x0 + real(i - 1)*grid%dx
      xvrayleigh = vrayleigh/1.e4_real64
      ray_ext = (xvrayleigh**3/(9.38076e2_real64 - 10.8426_real64*xvrayleigh**2))* &
                (wtot*conv_cm2mol)
      if (jrad .eq. 1) then
        ray_ext = ray_ext*xvrayleigh
      elseif (jrad .eq. 0) then
        ray_ext = ray_ext*xvrayleigh/radfn(vrayleigh, tave/radcn2)
      endif
      absrb(i) = absrb(i) + ray_ext
    enddo
  endif
end subroutine rayleigh_scattering


end module continuum
