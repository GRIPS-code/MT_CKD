!This program calculates the continuum optical depth for an homogeneous layer.
program mt_ckd
use, intrinsic :: iso_fortran_env, only: input_unit, output_unit, real64
use netcdf
use utils, only: nc_check, radfn
use continuum, only: contnm
implicit none


integer :: i, icflg, nmol, nptabs, npth, jrad, ncid, varid, err
integer, dimension(1) :: dimid
real(kind=real64), dimension(7) :: xcnt
real(kind=real64), dimension(5050) :: wavenumber, absrb, csh2o, cfh2o, csh2or, cfh2or
real(kind=real64), dimension(60) :: wk
real(kind=real64) :: pave, tave, vmrh2o, xlength, wtot, w_dry, wa, wn2, wbroad, &
                     v1abs, v2abs, dvabs, xkt, vi, v1h, dvh, radfld, v1, v2
!character(len=8), dimension(60) :: hmolid
!character(len=8), parameter :: holn2 = " OTHER  "
real(kind=real64), parameter :: radcn2 = 1.4387752
real(kind=real64), parameter :: alosmt = 2.6867775e19
character(len=120) :: path

!data hmolid/"  H2O   ", "  CO2   ", "   O3   ", "  N2O   ", &
!            "   CO   ", "  CH4   ", "   O2   ", 53*"        "/

path = "/home/rlm/MT_CKD/mt_ckd/mt-ckd.nc"
icflg = -999
do i = 1, 7
  xcnt(i) = 1._real64 !Controls which molecules are calculated:
enddo

do i = 1, 5050
  absrb(i) = 0._real64
enddo

do i = 1, 60
  wk(i) = 0._real64
enddo

write(output_unit, *)
write(output_unit, *) "This version is limited to 5000 values"

!
!   THE FOLLOWING QUANTITIES MUST BE SPECIFIED: 
!
!          PRESSURE                   PAVE (MB)
!
!          TEMPERATURE                TAVE ( K)
!
!          COLUMN AMOUNT
!            NITROGEN                 WN2    (MOLEC/CM**2)
!            OXYGEN                   WK(7)  (MOLEC/CM**2)
!            CARBON DIOXIDE           WK(2)  (MOLEC/CM**2)
!            WATER VAPOR              WK(1)  (MOLEC/CM**2)
!
!          NUMBER OF MOLECULES        NMOL
!
!          BEGINNING WAVENUMBER       V1ABS (CM-1)
!
!          ENDING WAVENUMBER          V2ABS (CM-1)
!
!          SAMPLING INTERVAL          DVABS (CM-1)
!
!          NUMBER OF VALUES           NPTABS
!
!
!   THE RESULTS ARE IN ARRAY ABSORB
!
!   NOTE THAT FOR AN ATMOSPHERI! LAYER: 
!
!            WTOT   = ALOSMT * (PAVE/1013) * (273/TAVE) * (PATH LENGTH)
!
!            WBROAD = the column amount for all species not explicitly provided
!
!            WK(M)  = (VOLUME MIXING RATIO) * (COLUMN OF DRY AIR)
!
!
!
!   THE FOLLOWING IS AN EXAMPLE FOR A ONE CM PATH (SEE CNTNM.OPTDPT FOR RESULTS)
!

pave = 1013. ![mb].
tave = 260. ![K].
vmrh2o = 0.01 ![mol mol-1]?
xlength = 1. ![cm].

write(output_unit, *)
write(output_unit, *) "*** USE phys_constsor this program, vmr_h2o is taken "  &
                      //" with respect to the total column ***"
write(output_unit, *)
!write(output_unit, *) " read: pressure (mb)  if negative use default values"
!read(input_unit, *) press_rd
!if (press_rd .gt. 0.) then
!  pave = press_rd
!  write(output_unit, *) " read:   temperature (K)"
!  read(input_unit, *) tave
!  write(output_unit, *) " read:   path length (cm)"
!  read(input_unit, *) xlength
!  write(output_unit, *) " read:   vmr h2o"
!  read(input_unit, *) vmrh2o
!endif
write(output_unit, *) "Pressure [mb], Temperature [K], Path Length [cm], VMR H2O"
write(output_unit, "(1x,f13.6,f17.4,f18.4,f12.8)") pave, tave, xlength, vmrh2o

!It may be preferable to specifiy the water column directly!
wtot = alosmt*(pave/1013.)*(273./tave)*xlength ![cm-2].
w_dry = wtot*(1. - vmrh2o) !dry air [cm-2].
wa = 0.009*w_dry !argon [cm-2].
wn2 = 0.78*w_dry !nitrogen [cm-2].
wk(7) = 0.21*w_dry !oxygen [cm-2].
wk(2) = 345.e-06*w_dry !carbon dioxide [cm-2].
!wk(3) = 3.75e-6*w_dry !ozone [cm-2].
wk(3) = 0. !ozone [cm-2].
if (abs(vmrh2o - 1.) .lt. 1.e-5) then
  wk(1) = wtot !water vapor [cm-2].
else
  wk(1) = vmrh2o*w_dry !water vapor [cm-2].
endif
wbroad = wn2 + wa !nitrogen plus argon [cm-2].
nmol = 7

!Spectral grid.
v1abs = 0. ![cm-1].
v2abs = 10000. ![cm-1].
dvabs = 2. ![cm-1].
nptabs = int(1. + (v2abs - v1abs)/dvabs)
do i = 1, nptabs
  absrb(i) = 0.
enddo

!write out inputs
!write(ipr, 970) pave, tave
!write(ipr, 975) (hmolid(i), i=1, 7), holn2
!write(ipr, 980) (wk(i), i=1, 7), wbroad
!write(ipu, 970) pave, tave
!write(ipu, 975) (hmolid(i), i=1, 7), holn2
!write(ipu, 980) (wk(i), i=1, 7), wbroad
!70 format(/, 29x, "P(MB)", 7X, "T(K)", //, 23x, 0P, F12.3, F9.2)
!75 format(/, 9X, "MOLECULAR AMOUNTS (MOL/CM**2) BY LAYER ", //, 8(1X, A6, 3X))
!80 format(/, 1P, 8E10.3, //)

jrad = 1
v1 = v1abs
v2 = v2abs
call contnm(icflg, nmol, nptabs, jrad, xcnt, pave, tave, wbroad, v1abs, &
            v2abs, dvabs, v1, v2, wk, npth, v1h, dvh, absrb, path, csh2o, cfh2o)

call nc_check(nf90_create("out.nc", ior(nf90_clobber, nf90_netcdf4), ncid))
call nc_check(nf90_def_dim(ncid, "wavenumber", nptabs, dimid(1)))
call nc_check(nf90_def_var(ncid, "wavenumber", nf90_double, dimid(1:1), varid))
call nc_check(nf90_put_att(ncid, varid, "units", "cm-1"))
wavenumber(:) = 0._real64
do i = 1, nptabs
  wavenumber(i) = v1abs + real(i - 1)*dvabs
enddo
call nc_check(nf90_put_var(ncid, varid, wavenumber(1:nptabs)))
call nc_check(nf90_def_var(ncid, "absrb", nf90_double, dimid(1:1), varid))
call nc_check(nf90_put_att(ncid, varid, "units", ""))
call nc_check(nf90_put_var(ncid, varid, absrb(1:nptabs)))


wavenumber(:) = 0._real64
xkt = tave/radcn2
do i = 1, npth
  vi = v1h + real(i - 1)*dvh
  if (vi .ge. v1abs .and. vi .le. v2abs) then
    radfld = radfn(vi, xkt)
    csh2or(i) = csh2o(i)*radfld
    cfh2or(i) = cfh2o(i)*radfld
  endif
  wavenumber(i) = vi
enddo


call nc_check(nf90_def_dim(ncid, "wavenumber2", npth, dimid(1)))
call nc_check(nf90_def_var(ncid, "wavenumber2", nf90_double, dimid(1:1), varid))
call nc_check(nf90_put_att(ncid, varid, "units", "cm-1"))
call nc_check(nf90_put_var(ncid, varid, wavenumber(1:npth)))

call nc_check(nf90_def_var(ncid, "csh2o", nf90_double, dimid(1:1), varid))
call nc_check(nf90_put_att(ncid, varid, "units", "cm3"))
call nc_check(nf90_put_att(ncid, varid, "longname", &
              "water vapor self continuum coefficient without radiation term"))
call nc_check(nf90_put_var(ncid, varid, csh2o(1:npth)))

call nc_check(nf90_def_var(ncid, "cfh2o", nf90_double, dimid(1:1), varid))
call nc_check(nf90_put_att(ncid, varid, "units", "cm3"))
call nc_check(nf90_put_att(ncid, varid, "longname", &
              "water vapor foreign continuum coefficient without radiation term"))
call nc_check(nf90_put_var(ncid, varid, cfh2o(1:npth)))

call nc_check(nf90_def_var(ncid, "csh2or", nf90_double, dimid(1:1), varid))
call nc_check(nf90_put_att(ncid, varid, "units", "cm2"))
call nc_check(nf90_put_att(ncid, varid, "longname", &
              "water vapor self continuum coefficient with radiation term"))
call nc_check(nf90_put_var(ncid, varid, csh2or(1:npth)))

call nc_check(nf90_def_var(ncid, "cfh2or", nf90_double, dimid(1:1), varid))
call nc_check(nf90_put_att(ncid, varid, "units", "cm2"))
call nc_check(nf90_put_att(ncid, varid, "longname", &
              "water vapor foreign continuum coefficient with radiation term"))
call nc_check(nf90_put_var(ncid, varid, cfh2or(1:npth)))

call nc_check(nf90_close(ncid))
end program mt_ckd
