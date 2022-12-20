!-------------------------------------------
! MODULE: Set precision
! real32  -> simple precision
! real64  -> double precision
! real128 -> quadruple precision
!-------------------------------------------
module precision
   use ISO_FORTRAN_ENV
   implicit none
   integer, parameter :: wp = real64

end module precision