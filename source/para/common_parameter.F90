module common_parameter
!! This module defines the maximum size of all arrays, number of dimensions,
!! and expected stencil sizes for different kernels
  use kind_parameters
  implicit none 

  !! Maximum number of particles (my workstation memory is limited to ~1.8e7
  integer(ikind) ,parameter :: npar=17500000 
  integer(ikind) ,parameter :: nparb=2000000 !! Max number of bubbles 
  
  ! dims is # of dimensions, must=2 or 3
  ! nplink is max allowable number of neighbours (will throw seg fault if breached)

#if dim3
  integer(ikind) ,parameter :: dims = 3
#if kernel==2
  integer(ikind) ,parameter :: nplink=150   !! SPH Wendland 3D 
#else
  integer(ikind) ,parameter :: nplink=300  ! Quintic kernel requires more neighbours
#endif
#else  
  integer(ikind) ,parameter :: dims = 2
#if kernel==2
  integer(ikind) ,parameter :: nplink=50    !! SPH Wendland 2D
#else
  integer(ikind) ,parameter :: nplink=100  ! Quintic kernel requires more neighbours
#endif
#endif

  integer(ikind) ,parameter :: nclink=600 !! Max number of particles in a cell...

  !! Define some common numbers
  real(rkind), parameter :: pi=3.14159265358979323846d0
  real(rkind), parameter :: one=1.0_rkind
  real(rkind), parameter :: zero=0.0_rkind
  real(rkind), parameter :: half=0.5_rkind
  real(rkind), parameter :: quarter=0.25_rkind
  real(rkind), parameter :: two = 2.0_rkind
  real(rkind), parameter :: three = 3.0_rkind  
  real(rkind), parameter :: four = 4.0_rkind  
  real(rkind), parameter :: five = 5.0_rkind    
  real(rkind), parameter :: onethird = one/three      
  real(rkind), parameter :: twothirds = two/three      
  real(rkind), parameter :: twofifths = two/five
  real(rkind), parameter :: threefifths = three/five
  real(rkind), parameter :: fivethirds = five/three
  real(rkind), parameter :: twoninths = two/9.0_rkind    
  real(rkind), parameter :: vsmall = epsilon(one)



end module common_parameter
