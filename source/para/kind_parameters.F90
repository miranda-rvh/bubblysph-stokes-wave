module kind_parameters
!! This module defines the precision level - single or double - according
!! to flags set in the Makefile
  use iso_c_binding

  implicit none

  private
  public :: rkind,ikind,dkind

#if single
  integer ,parameter :: rkind =c_float 
#else
  integer ,parameter :: rkind =c_double 
#endif  
  integer ,parameter :: dkind =c_double

  integer ,parameter :: ikind =c_int    

end module
