!This program calculates the continuum optical depth for an homogeneous layer.
program mt_ckd
use, intrinsic :: iso_fortran_env, only: input_unit, output_unit, real64
use netcdf
use utils, only: nc_check, radfn, SpectralGrid
use continuum, only: contnm, h2o, co2, o3, o2, n2, ar
implicit none


integer :: i, icflg, jrad, ncid, varid
integer, dimension(1) :: dimid
real(kind=real64), dimension(7) :: xcnt
real(kind=real64), dimension(:), allocatable :: wavenumber, absrb, csh2o, cfh2o, csh2or, cfh2or
real(kind=real64) :: pave ! Pressure [mb].
real(kind=real64) :: tave ! Temperature [K].
real(kind=real64) :: path_length !Path length [cm-1].
real(kind=real64) :: xkt, vi, radfld
real(kind=real64), parameter :: radcn2 = 1.4387752_real64
character(len=120) :: path
real(kind=real64), dimension(6) :: vmr ! Volume mixing ratio in dry air [mol mol-1].
type(SpectralGrid) :: grid, h2o_grid


!User input.
path = "/home/rlm/MT_CKD/mt_ckd/new-file.nc"
pave = 1013. ![mb].
tave = 260. ![K].
vmr(h2o) = 0.01_real64 ![mol mol-1].
vmr(co2) = 3.45e-4_real64 ![mol mol-1].
vmr(o3) = 3.75e-6_real64 ![mol mol-1].
!vmr(o3) = 0._real64 ![mol mol-1].
vmr(o2) = 0.21_real64 ![mol mol-1].
vmr(n2) = 0.78_real64 ![mol mol-1].
vmr(ar) = 0.009_real64 ![mol mol-1].
path_length = 1. ![cm].
call grid%construct(0._real64, 1.e5_real64, 2._real64) !Spectral grid.

!Initialization.
icflg = -999
xcnt(:) = 1._real64 !Controls which molecules are calculated:
absrb(:) = 0._real64
write(output_unit, *)
write(output_unit, *) "This version is limited to 5000 values"
write(output_unit, *) "Pressure [mb], Temperature [K], Path Length [cm], VMR H2O"
write(output_unit, "(1x,f13.6,f17.4,f18.4,f12.8)") pave, tave, path_length, vmr(h2o)
jrad = 1

allocate(absrb(grid%n))
call contnm(pave, tave, vmr, path_length, grid, icflg, jrad, xcnt, &
            h2o_grid, absrb, path, csh2o, cfh2o)

!Write out the data.
call nc_check(nf90_create("out.nc", ior(nf90_clobber, nf90_netcdf4), ncid))
call nc_check(nf90_def_dim(ncid, "wavenumber", grid%n, dimid(1)))
call nc_check(nf90_def_var(ncid, "wavenumber", nf90_double, dimid(1:1), varid))
call nc_check(nf90_put_att(ncid, varid, "units", "cm-1"))
allocate(wavenumber(grid%n))
wavenumber(:) = 0._real64
do i = 1, grid%n
  wavenumber(i) = grid%x0 + real(i - 1)*grid%dx
enddo
call nc_check(nf90_put_var(ncid, varid, wavenumber))
call nc_check(nf90_def_var(ncid, "absrb", nf90_double, dimid(1:1), varid))
call nc_check(nf90_put_att(ncid, varid, "units", ""))
call nc_check(nf90_put_var(ncid, varid, absrb))
deallocate(wavenumber, absrb)

if (xcnt(1) .gt. 0._real64) then
  allocate(wavenumber(h2o_grid%n), csh2or(h2o_grid%n))
  wavenumber(:) = 0._real64
  xkt = tave/radcn2
  do i = 1, h2o_grid%n
    vi = h2o_grid%x0 + real(i - 1)*h2o_grid%dx
    if (vi .ge. grid%x0 .and. vi .le. grid%xn) then
      radfld = radfn(vi, xkt)
      csh2or(i) = csh2o(i)*radfld
    endif
    wavenumber(i) = vi
  enddo
  call nc_check(nf90_def_dim(ncid, "wavenumber2", h2o_grid%n, dimid(1)))
  call nc_check(nf90_def_var(ncid, "wavenumber2", nf90_double, dimid(1:1), varid))
  call nc_check(nf90_put_att(ncid, varid, "units", "cm-1"))
  call nc_check(nf90_put_var(ncid, varid, wavenumber))
  call nc_check(nf90_def_var(ncid, "csh2o", nf90_double, dimid(1:1), varid))
  call nc_check(nf90_put_att(ncid, varid, "units", "cm3"))
  call nc_check(nf90_put_att(ncid, varid, "longname", &
                "water vapor self continuum coefficient without radiation term"))
  call nc_check(nf90_put_var(ncid, varid, csh2o))
  call nc_check(nf90_def_var(ncid, "csh2or", nf90_double, dimid(1:1), varid))
  call nc_check(nf90_put_att(ncid, varid, "units", "cm2"))
  call nc_check(nf90_put_att(ncid, varid, "longname", &
                "water vapor self continuum coefficient with radiation term"))
  call nc_check(nf90_put_var(ncid, varid, csh2or))
  deallocate(csh2o, wavenumber, csh2or)
endif
if (xcnt(2) .gt. 0._real64) then
  allocate(wavenumber(h2o_grid%n), cfh2or(h2o_grid%n))
  wavenumber(:) = 0._real64
  xkt = tave/radcn2
  do i = 1, h2o_grid%n
    vi = h2o_grid%x0 + real(i - 1)*h2o_grid%dx
    if (vi .ge. grid%x0 .and. vi .le. grid%xn) then
      radfld = radfn(vi, xkt)
      cfh2or(i) = cfh2o(i)*radfld
    endif
    wavenumber(i) = vi
  enddo
  call nc_check(nf90_def_dim(ncid, "wavenumber3", h2o_grid%n, dimid(1)))
  call nc_check(nf90_def_var(ncid, "wavenumber3", nf90_double, dimid(1:1), varid))
  call nc_check(nf90_put_att(ncid, varid, "units", "cm-1"))
  call nc_check(nf90_put_var(ncid, varid, wavenumber))
  call nc_check(nf90_def_var(ncid, "cfh2o", nf90_double, dimid(1:1), varid))
  call nc_check(nf90_put_att(ncid, varid, "units", "cm3"))
  call nc_check(nf90_put_att(ncid, varid, "longname", &
                "water vapor foreign continuum coefficient without radiation term"))
  call nc_check(nf90_put_var(ncid, varid, cfh2o))
  call nc_check(nf90_def_var(ncid, "cfh2or", nf90_double, dimid(1:1), varid))
  call nc_check(nf90_put_att(ncid, varid, "units", "cm2"))
  call nc_check(nf90_put_att(ncid, varid, "longname", &
                "water vapor foreign continuum coefficient with radiation term"))
  call nc_check(nf90_put_var(ncid, varid, cfh2or))
  deallocate(cfh2o, wavenumber, cfh2or)
endif
call nc_check(nf90_close(ncid))

end program mt_ckd
