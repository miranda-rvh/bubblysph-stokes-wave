module common_2d

  use iso_c_binding
  use common_parameter
  use kind_parameters      

  implicit none
  
  !! Governing parameters
  real(rkind) :: Re,Ma,beta,Fr,We,Maa

 
  integer(ikind) ::  i_open_domain  
  
  !! Particle properties
  real(rkind), dimension(:), allocatable, target :: p,conc,a_out,vort,vol,h
  real(rkind),dimension(:,:),allocatable, target :: rp, up
  real(rkind), dimension(:,:), allocatable :: rpo, upo
  real(rkind), dimension(:,:,:), allocatable :: grad_vel,kgcm
  real(rkind), allocatable, dimension(:,:) :: grad_p  
  integer(ikind), dimension(:), allocatable, target :: n_surf,c_func
  logical( c_bool ), dimension(:), allocatable, target :: inbin
  real(rkind),allocatable,dimension(:) :: ppe_smooth
  real(rkind), bind(c) :: dx, dy, dv, av_conc
  integer(ikind), dimension(:), allocatable, target :: irelation,vrelation
  real(rkind),dimension(:),allocatable, target :: dP_mp
  real(rkind),dimension(:,:),allocatable, target :: surf_norm,grad_conc
  real(rkind), dimension(:,:), allocatable,target :: ushift
  
  !! VOF scheme for bubbly flows
  real(rkind),dimension(:),allocatable :: a0,a00   !! Vol frac of liquid (phase 0)
  real(rkind),dimension(:,:),allocatable :: M0 !! momentum transfer to phase 0 (from all bubble groups)
  real(rkind),dimension(:),allocatable :: nu_BI !! bubble induced turbulent viscosity
  real(rkind),dimension(:),allocatable,target :: eps_srs,tke
  
  !! LES stuff
  real(rkind),dimension(:),allocatable,target :: I_LM,I_MM  
  real(rkind) :: delta_F
  
  !! Misc
  real(rkind) :: t_impact
   
  
  !! Control flags useful for debugging
  logical :: output_everytime ! Output every time step
  logical :: fick_shifting  ! Use shifting
  logical :: advection_terms  ! Include advection terms
  logical :: output_mirrors  ! Include mirror particles in output files
  logical :: lagrangian    ! Lagrangian or Eulerian  


  integer(ikind) n_inbin
  
  !! important physical parameters
  real(rkind), dimension(:),allocatable, target :: grav
  
  !! Profiling
  real(rkind) :: tTot,tffa,tppeb,tppes,tt1,tt2  


  !! jack's boundary condition framework
  real(rkind),dimension(:,:),allocatable, target :: b_node,b_edge
  integer(ikind),dimension(:),allocatable,target :: b_type,b_periodic_parent
  real(rkind),dimension(:),allocatable,target :: b_corner_angle
  integer(ikind),dimension(:),allocatable,target :: p_inflow
  integer(ikind),dimension(:,:),allocatable,target :: rp_inflow  
  real(rkind) :: U_inflow,U_inflow0
  integer(ikind) :: nb_patches,np_inflow
  integer(ikind),dimension(:),allocatable,target :: p_free
  integer(ikind) :: n_free
  
  !! Wavemaker X-position and X-velocity
  real(rkind) :: X_wmaker,U_wmaker
  real(rkind) :: x_target

  !! linear solver and PPE stuff
  real(rkind), dimension(:), allocatable, target :: lhs_mat,rhs_vec
  integer(ikind), dimension(:), allocatable, target :: lhs_index
  integer(ikind) :: LS_iters,nnz_old
  real(rkind) :: LS_residual

  !! time steps
  real(rkind), bind(c) :: dt
  real(rkind) :: dt_coef_advection,dt_coef_viscous,dt_coef_acoustic! time-step coefficients (ie. cfl number)
  real(rkind) :: umax
  real(rkind) :: time, dt2, dt_out,next_dump,next_dump2,time_max
  integer(ikind), bind(c) :: itime,itmax,n_dump,n_dump2
  logical :: output_this_step

  real(rkind), bind(c) :: D,r0
  integer(ikind) :: nmirror,n_count,nbsr,n_diff
  integer(ikind) :: n_par_fwm,n_par_fw,n_par_w
  integer(ikind), bind(c) :: nmirror_esti,n_par_fw_old,maxneighbours,maxneighbours_old

  real(rkind), bind(c) :: xmin,xmax,ymin,ymax

  !! Cell numbers and link lists
  integer(ikind), bind(c) :: ncx,ncy,nct,ncz
  integer(ikind),dimension(:),allocatable :: ic_count
  integer(ikind),dimension(:,:),allocatable :: ic_link

  !! Particle link lists and derivative operator arrays
  integer(ikind),dimension(:),allocatable,target :: ij_count
  integer(ikind),dimension(:,:),allocatable,target :: ij_link
  real(rkind),dimension(:,:),allocatable :: ij_w_L
  real(rkind),dimension(:,:,:),allocatable :: ij_w_G
  
!! Bubbly flows stuff -------------------------------------------
  !! Link lists for bubbles
  integer(ikind),dimension(:),allocatable,target :: ijb_count,ijb_linearlink,ijb_firstindex
  integer(ikind),dimension(:,:),allocatable,target :: ijb_link
  !! Liquid "DENSITY" at bubble locations
  real(rkind),dimension(:),allocatable,target :: bnfs
  integer(ikind),dimension(500) :: i_bubble_trace
  real(rkind),dimension(:),allocatable :: deformation_distance,mean_dissrate
  
  !! Position and velocity
  real(rkind),dimension(:,:),allocatable,target :: rb,ub,ubturb,u_l2b   
  !! Bubble size
  real(rkind),dimension(:),allocatable :: radb,b_out
  real(rkind),dimension(:),allocatable :: bubble_age,bubble_LE  !! age and life expectancy...
  integer(ikind),dimension(:),allocatable :: bubble_EoL  !! Flag if in end of life
  !! Numbers etc
  real(rkind) :: Rhinze
  integer(ikind) :: n_bub,n_bub_m !! number of bubbles, number of bubbles plus mirrors
  !! Relations for boundary conditions  
  integer(ikind), dimension(:), allocatable, target :: irelation_b,vrelation_b
  !! Destroyed bubbles and bubbles-in-bin
  integer(ikind),dimension(:),allocatable :: b_freeindices
  integer(ikind) :: b_nfree
  logical,dimension(:),allocatable :: b_inbin
  
 
!! ---------------------------------------------------------------

  real(rkind), bind(c) :: h0,uno_cell_size,sup_size,sup_size_2,cell_size

  real(rkind), bind(c) :: ss,ss2

  real(rkind), bind(c) ::  xb_min,xb_max,yb_min,yb_max,zb_min,zb_max
  real(rkind), bind(c) :: xmint,xmaxt,ymint,ymaxt,zmint,zmaxt                    

  real(rkind), bind(c) :: ad_7,ad_7h

end module common_2d
