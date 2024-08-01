      module global_variables 
      use kind_parameters      
      implicit none

      !! Dimensionless governing parameters
      real(rkind) :: Re,Ma,beta,Fr,We

      integer(ikind) :: itime, itmax
      integer(ikind) :: n_par_fw,n_par_w
      integer(ikind) :: i_open_domain

      real(rkind) :: dx, dy, h, time, dt,  dt_out
      real(rkind) :: grx, gry
      real(rkind) :: xb_min, xb_max, yb_min, yb_max, xl
      real(rkind), dimension(:), allocatable :: xp, yp ,up, vp,wp,zp, p, vol

      !! jack's boundary condition framework
      real(rkind),dimension(:,:),allocatable, target :: b_node,b_edge
      integer(ikind),dimension(:),allocatable,target :: b_type,b_periodic_parent
      integer(ikind) :: nb_patches
     
      logical :: slip_walls
     
      end module global_variables 
