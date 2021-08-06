module utils
use, intrinsic :: iso_fortran_env, only: error_unit, real64
use netcdf
implicit none


public :: nc_check, pre_xint, radfn, xint


contains


!> @brief Exits the program if a netCDF error occurs.
subroutine nc_check(err)

  integer, intent(in) :: err !< NetCDF return code.

  if (err .ne. nf90_noerr) then
    write(error_unit, *) "Error: "//nf90_strerror(err)
    stop 1
  endif
end subroutine nc_check


subroutine pre_xint(v1ss, v2ss, v1abs, dvabs, nptabs, ist, last)

  real(kind=real64), intent(in) :: v1ss, v2ss, v1abs, dvabs
  integer, intent(in) :: nptabs
  integer, intent(out) :: ist, last

  integer :: nbnd_v1c, nbnd_v2c
!   set up needed variables for call to xint
!   output variables
!     nptabs_loc - number of values to be processed in xint
!     ist - index of first value to be processed in xint

  nbnd_v1c = 2 + (v1ss - v1abs)/dvabs + 1.e-5
  ist = max(1, nbnd_v1c)

  nbnd_v2c = 1 + (v2ss - v1abs)/dvabs + 1.e-5
  last = min(nptabs, nbnd_v2c)
end subroutine pre_xint


subroutine xint(v1a, v2a, dva, a, afact, vft, dvr3, r3, n1r3, n2r3)

  real(kind=real64), intent(in) :: v1a, v2a, dva
  real(kind=real64), dimension(:), intent(in) :: a
  real(kind=real64), intent(in) :: afact, vft, dvr3
  real(kind=real64), dimension(:), intent(inout) :: r3
  integer, intent(in) :: n1r3, n2r3
!     this subroutine interpolates the a array stored
!     from v1a to v2a in increments of dva using a multiplicative
!     factor afact, into the r3 array from location n1r3 to n2r3 in
!     increments of dvr3

  integer :: i, ilo, ihi, j
  real(kind=real64) :: b, b1, b2, c, conti, p, recdva, vi, vj

  recdva = 1./dva
  ilo = (v1a + dva - vft)/dvr3 + 1. + 0.999
  ilo = max(ilo, n1r3)
  ihi = (v2a - dva - vft)/dvr3 + 0.999
  ihi = min(ihi,n2r3)

  do i = ilo, ihi
    vi = vft + dvr3*real(i - 1)
    j = (vi - v1a)*recdva + 1.001
    vj = v1a + dva*float(j - 1)
    p = recdva*(vi - vj)
    c = (3. - 2.*p)*p*p
    b = 0.5*p*(1. - p)
    b1 = b*(1. - p)
    b2 = b*p
    conti = -a(j-1)*b1 + a(j)*(1. - c + b2) + a(j+1)*(c + b1) - a(j+2)*b2
    r3(i) = r3(i) + conti*afact
  enddo
end subroutine xint


!> @brief Calculates the radiation term for the line shape.
function radfn(vi, xkt) result(r)

  real(kind=real64), intent(in) :: vi !< Wavenumber [cm-1].
  real(kind=real64), intent(in) :: xkt !< Temperature divided by the second
                                       !! radiation constant [cm-1].
  real(kind=real64) :: r !< Radiation term [cm-1].

  real(kind=real64) :: expvkt, xvi, xviokt

  xvi = vi
  if (xkt .gt. 0.0) then
    xviokt = xvi/xkt
    if (xviokt .le. 0.01) then
      r = 0.5*xviokt*xvi
    elseif (xviokt .le. 10.0) then
      expvkt = exp(-xviokt)
      r = xvi*(1. - expvkt)/(1. + expvkt)
    else
      r = xvi
    endif
  else
    r = xvi
  endif
end function radfn


end module utils
