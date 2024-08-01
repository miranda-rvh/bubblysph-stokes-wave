module input
  use common_parameter
  use common_2d
  use cell_griddin
  use bubble_acoustics
  implicit none

  private 
  public :: getdata

contains
!! ------------------------------------------------------------------------------------------------
  subroutine getdata
    real(rkind) :: re_time,p_dist,error_mag,dummy_real,tmp_real
    real(rkind),dimension(dims) :: tmp
    real(rkind) :: x,y,z,fpios,pios,beltrami_mag,dz
    integer(ikind) :: i,j,ib
    real(rkind),dimension(:,:),allocatable :: rptmp,uptmp
    integer(ikind) :: n_par_fwtmp,nz,n_par_wtmp,ii,dummy_int
    logical :: keepgoing

    !! Check compiler flag option compatibility...
#ifdef in3d
#ifndef dim3
    write(6,*) "Compiled for 2D simulation with 3D input files. Incompatible. Stopping."
    stop
#endif
#endif

    !! ====================================================
    !! read in data on run parameters
    allocate(grav(dims));grav=zero


    open(11,file='control.in')
    
    !! Read header
    read(11,*)
    read(11,*)    
    
    !! Main non-dimensional parameters
    read(11,*)    
    read(11,*) Re
    read(11,*)
    
    read(11,*)        
    read(11,*) Ma
    read(11,*)
    
    read(11,*) 
    read(11,*) Maa
    read(11,*) 
    
    read(11,*)
    read(11,*) beta  
    read(11,*)
    
    read(11,*)
    read(11,*) Fr
    read(11,*)
    
    read(11,*)        
    read(11,*) grav(1),grav(2)
    read(11,*)
    
    read(11,*)
    read(11,*) We
    read(11,*)
    
    read(11,*)
    read(11,*) time,time_max
    read(11,*)
    itime = 0
    
    read(11,*)
    read(11,*) dt_out    
    read(11,*)
    
    read(11,*)
    read(11,*) dt_coef_advection,dt_coef_viscous,dt_coef_acoustic    
    read(11,*)    
    
    close(11)
    
    !! Load ipart metadata
    open(unit=13,file='IPART')
    
    read(13,*) keepgoing
    !! Check whether compilation and setup data match in terms of slip/no-slip
#ifdef slip
    if(.not.keepgoing) then
       write(6,*) "Warning, setup is for no slip, but sph compiled for slip. Stopping."
       stop
    end if
#else
    if(keepgoing) then
       write(6,*) "Warning, setup is for slip, but sph compiled for no slip. Stopping."
       stop
    end if
#endif     


         
    !! Some SPH and simulation controls     
    read(13,*) i_open_domain
    read(13,*) n_par_fw
    read(13,*) n_par_w
    read(13,*) h0
    read(13,*) dx
#ifdef in3d
    read(13,*) dy,dz
#else
    read(13,*) dy
#endif    
    read(13,*) nz
    read(13,*) zb_max
#ifdef in3d
    zb_max = half*zb_max;zb_min = -zb_max
    write(6,*) "Loading 3D input files, 3D dimension domain sizes::"
    write(6,*) "zb_min,zb_max",zb_min,zb_max
#endif    
    
    !! Re-scale gravity vector with Fr^2 
    grav(:) = grav(:)/Fr/Fr

    !! Set the size of the implicit (SPH) filter    
    delta_F = h0!dx !! Width of SPH (grid) filter        
     
    !! estimate the number of mirror particles - this is useful to reduce costs associated
    !! with memory allocation later
    nmirror_esti=int(0.5*n_par_fw)
    n_par_fw_old = 1
    maxneighbours_old = 1
    nnz_old = 1
    
    !! Space for ij weights and kernel gradient correction
    allocate(ij_w_G(1,1,1),ij_w_L(1,1),kgcm(1,1,1))

    !! read in data on boundaries
    open(10,file='IBOUND')
    read(10,*) nb_patches
    allocate(b_node(nb_patches,2),b_edge(nb_patches,2))
    allocate(b_type(nb_patches),b_periodic_parent(nb_patches))
    do i=1,nb_patches
       read(10,*) b_node(i,:)
       read(10,*) b_edge(i,:)
       read(10,*) b_type(i)
       if(b_type(i)==2)then !periodic boundary patch
          read(10,*) b_periodic_parent(i)
       end if
       if(b_type(i)==3)then ! inflow patch
          read(10,*) U_inflow
          U_inflow0=U_inflow
       end if      
    end do
    close(10)
    write(6,*) "Boundaries read in. There are ",nb_patches," boundary patches"

    !! Initialise wavemaker position (regardless of whether we have a wavemaker)    
    X_wmaker=b_node(1,1)
    U_wmaker=zero  

    !! Various checks to ensure the boundaries are correctly specified
    call check_boundaries

    !! Create an X-Y aligned box which contains the domain
    xb_min = minval(b_node(:,1))
    xb_max = maxval(b_node(:,1))
    yb_min = minval(b_node(:,2))
    yb_max = maxval(b_node(:,2))
    write(6,*) "xb_min, xb_max=",xb_min,xb_max
    write(6,*) "yb_min, yb_max=",yb_min,yb_max

    !! Initialise the time-step to something small, just in case
    dt=1.0d-8

    !! Allocate space and read in particle data, either from...
    allocate(rp(npar,dims));rp=zero      !! Position
    allocate(up(npar,dims));up=zero      !! Velocity
    allocate(p(npar),a_out(npar))        !! Pressure and spare scalar for output
    allocate(vol(npar),h(npar))           !! Volume and smoothing length
    allocate(n_surf(npar))                !! Free surface flag
    allocate(inbin(npar))                 !! Bin flag
    allocate(a0(npar));a0=one              !! Liquid volume fraction    

    !! a DUMP file
#ifdef restart
    open(14,file='DUMP')       
    read(14,*)dummy_int,dummy_real,dummy_int
    read(14,*)n_par_fwm,n_par_fw,n_par_w,np_inflow
    read(14,*)xb_min,xb_max,yb_min,yb_max, dx, dy, av_conc
    print*,'CASE RESTART FROM DUMP FILE'
    do i=1,n_par_fw
#ifdef dim3
       read(14,*)rp(i,1), rp(i,2),rp(i,3)
#else
       read(14,*)rp(i,1), rp(i,2)
#endif               
    enddo
    if(np_inflow/=0)then
       allocate(p_inflow(np_inflow+2))
       do i=1,np_inflow
          read(14,*) p_inflow(i)
       end do
    end if       
    close(14)
    !! Or the input file
#else
!    open(13,file='IPART')
    do i=1,n_par_fw
#ifdef in3d
       read(13,*)rp(i,1),rp(i,2),rp(i,3),up(i,1),up(i,2),up(i,3),dv
       dv = dv*dz
#else    
       read(13,*)rp(i,1), rp(i,2), up(i,1), up(i,2), dv
#endif       
    enddo    
    close(13)       
#endif

    !! Initialise all fluid particles to free surface if open domain
    n_surf(:) = 0
    if(i_open_domain==0)then
       n_surf(1:n_par_fw) = 0
    else
       n_surf(1:n_par_fw) = 1
    end if
    

    !! Build a 3D domain (overwriting a chunk of the above) if required
    !! but only if not loading 3d input files
#ifndef in3d
#if dim3   
    !! Adjust the particle volume
    dv = dx*dx*dx
    !! Determine the domain length in Z
    zb_min = (-half*(nz-1) - half)*dx;zb_max = (half*(nz-1) + half)*dx
    write(6,*) "zb_min,zb_max",zb_min,zb_max
    !! Either load in from a DUMP file...
#ifdef restart
    open(36,file='DUMP')
    read(36,*) dummy_int,dummy_real,dummy_int
    read(36,*) dummy_int,n_par_fw,dummy_int,dummy_int
    read(36,*) dummy_real,dummy_real,dummy_real,dummy_real,dummy_real,dummy_real,dummy_real
    write(6,*) "Three dimensional simulation, with ",n_par_fw,"particles!!!!!"    
    do i=1,n_par_fw
       read(36,*) rp(i,1),rp(i,2),rp(i,3),up(i,1),up(i,2),up(i,3),p(i)

       !! Sometimes we want to overwrite the velocity field 
!        up(i,:) = init_velocity_field(rp(i,1),rp(i,2),rp(i,3))
    end do      

    close(36)
#else    
    !! or make nz copies of the 2D domain
    allocate(rptmp(npar,dims),uptmp(npar,dims))     
    !! Make boundary particles first... 
    j=0
    do ii=1,nz
       do i=1,n_par_w
          j=j+1
          rptmp(j,:) = rp(i,:)
          rptmp(j,3) = zb_min + half*dx + dble(ii-1)*dx
          uptmp(j,:) = up(i,:)      
       end do
    end do
    n_par_wtmp = j
    !! Now make fluid particles...
    do ii=1,nz
       do i=n_par_w+1,n_par_fw
          j=j+1
          rptmp(j,:) = rp(i,:)
          rptmp(j,3) = zb_min + half*dx + dble(ii-1)*dx
          
          !! Copy velocity from 2D slice
          uptmp(j,:) = up(i,:)
          
          !! Prescribe new velocity field
!          uptmp(j,:) = init_velocity_field(rptmp(j,1),rptmp(j,2),rptmp(j,3))

       end do
    end do
    n_par_fwtmp = j
    n_par_w=n_par_wtmp
    !! estimate the number of mirror particles
    nmirror_esti = j*2 + 5*n_par_fw
write(6,*) "Estimated mirrors is ",nmirror_esti  
    n_par_fw = n_par_fwtmp

    if(n_par_fw<=1000) nmirror = n_par_fw*8
    do i=1,n_par_fw
       rp(i,:)=rptmp(i,:)
       up(i,:)=uptmp(i,:)     
     end do
    deallocate(rptmp,uptmp)     
      
    write(6,*) "Three dimensional simulation, with ",n_par_fw,"particles!"    
#endif
#endif    
#endif

    !! In some cases, we may want to prescribe an initial pressure field 
    !! only really for Taylor-Green 
    if(.false.)then  
       do i=1,n_par_fw
          p(i) = quarter*(cos(two*rp(i,1)) + cos(two*rp(i,2)))
       end do
    end if
   
    !! This section initialises bubbles
#ifdef bubbles    
    !! Initialise acoustic arrays
    Rhinze = half*We**(-threefifths)                  !! Hinze-radius estimate
    call create_acoustic_array
 
    allocate(rb(nparb,dims),ub(nparb,dims));ub=zero   !! Bubble velocity and position
    allocate(ubturb(nparb,dims));ubturb=zero     !! Turbulent fluctuations for bubble interactions    
    allocate(u_l2b(nparb,dims));u_l2b = zero !! Used to store previous fluid velocity at bubble location
    allocate(radb(nparb));radb = Rhinze         !! initialise bubble radii to Hinze radii for now
    allocate(b_out(nparb));b_out = zero          !! Spare field for outputting
    !! When impact surface, how long to coallesce, has it impacted surface?  
    allocate(bubble_age(nparb),bubble_LE(nparb),bubble_EoL(nparb))
    bubble_age=zero;bubble_LE=1.0d10;bubble_EoL=0
    t_impact=-one !! initialise t_impact as -ve...
    
    !! For bubble breakup model
    allocate(deformation_distance(nparb));deformation_distance = zero
    allocate(mean_dissrate(nparb));mean_dissrate = zero
   
    !! Initialise bubbles at start (if any)
    n_bub = 0! 64*64*64!10970!21940
    if(.false.)then
       !! Initial random bubble placement - normal distribution - rand + Box-Muller transform
       !! this gives a sort of Gaussian cloud of bubbles - lengthscales hard-coded
       do i=1,n_bub
          keepgoing = .true.
          do while(keepgoing)
             call random_number(dummy_real);call random_number(tmp_real)   
             rb(i,1) = 0.05*cos(two*pi*tmp_real)*sqrt(-two*log(dummy_real))
             if(rb(i,1)>=xb_min.and.rb(i,1)<=xb_max) keepgoing = .false.
          end do 
          keepgoing = .true.         
          do while(keepgoing)
             call random_number(dummy_real);call random_number(tmp_real)   
             rb(i,2) = 0.05*cos(two*pi*tmp_real)*sqrt(-two*log(dummy_real)) 
             if(rb(i,2)>=yb_min.and.rb(i,2)<=yb_max) keepgoing = .false.
          end do
#if dim3       
          keepgoing = .true.         
          do while(keepgoing)
             call random_number(dummy_real);call random_number(tmp_real)   
             rb(i,3) = 0.02*cos(two*pi*tmp_real)*sqrt(-two*log(dummy_real))
             if(rb(i,3)>=zb_min.and.rb(i,3)<=zb_max) keepgoing = .false.
          end do       
#endif       
       end do
    else if(.false.)then  !! Random within bounds and uniformly distributed
       do i=1,n_bub
          dummy_real = rand();rb(i,1) = xb_min + (xb_max-xb_min)*dummy_real
          dummy_real = rand();rb(i,2) = yb_min + (yb_max-yb_min)*dummy_real
#if dim3
          dummy_real = rand();rb(i,3) = zb_min + (zb_max-zb_min)*dummy_real                   
#endif          
          !! Zero bubble velocity
          ub(i,:)=zero;u_l2b(i,:) = ub(i,:) 

          !! Prescribe initial velocity to match liquid?          
          ub(i,:) = init_velocity_field(rb(i,1),rb(i,2),rb(i,3))
          u_l2b(i,:) = init_velocity_field(rb(i,1),rb(i,2),rb(i,3))                              

          radb(i)=0.3d0*dx
          ub(i,2)=0.1d0;rb(i,1)=zero;rb(i,2)=0.80d0;rb(i,3)=zero;radb(i)=0.003d0 !! Single bubble
       end do 
    end if
                
    !! Book-keeping: Initialise free index list, bubbles in bin
    allocate(b_inbin(nparb));b_inbin = .false.
    allocate(b_freeindices(nparb));b_freeindices = 0
    b_nfree = 0
    !! We may wish to plot the trajectories of individual bubbles. Do this through bubble traces.
    do i=1,500   !! Initialise the bubbble traces
       i_bubble_trace(i)=i*10
    end do
#endif    
    
    !! Back to liquid. Here initialise arrays used in sub-resolution-scale turbulence modelling
    allocate(I_LM(n_par_fw),I_MM(n_par_fw));I_LM=1.0d-16;I_MM=1.d-8 !! tracked properties for Lagrangian averaging
    allocate(eps_srs(npar));eps_srs=zero       !! sub-resolution turb dissipation rate
  
    !! Setup which particles are inflow particles...
#ifndef restart
    call set_inflow_particles
    n_free = 0
#endif       
    allocate(p_free(npar))
       
 
    !! Initial volumes and smoothing lengths
    vol(:)=dv
    h(:)=h0

    !! Set support radius and square of support radius
    !! also set ad_7 and ad_7h, which are kernel normalisation factors
#if kernel==1
    ss = three;ss2=ss*ss
    sup_size     = ss*h0                      !! QUINTIC SPLINE
    sup_size_2 = ss2*h0*h0
#if dim3
    ad_7 = one/120.0d0/pi/h0/h0/h0  !! 3D kernel
#else
    ad_7=7.0d0/478.0d0/pi/h0/h0    !! 2D kernel
#endif       
#elif kernel==2
    ss = two;ss2=ss*ss
    sup_size     = ss*h0                      !! WENDLAND C2
    sup_size_2 = ss2*h0*h0
#if dim3
    ad_7 = 21.0/16.0d0/pi/h0/h0/h0  !! 3D kernel
#else
    ad_7=7.0d0/four/pi/h0/h0    !! 2D kernel
#endif       
#endif    
    ad_7h=ad_7/h0

    !! Cell size  -  increased by factor of 1.26 allows volume fractions to approx 50%
    cell_size = ss*h0*1.26d0
    uno_cell_size = one/cell_size
 

    !! initially, all particles are in the domain. If they leave for
    !! some reason, or something goes wrong, a particle will be placed 'in the bin'
    inbin(:) = .false.
    
    
    !! open some files for outputting
    open(unit=92,file='./data_directory/TIME')          !! Contains time,dt,max(u) etc every step
    open(unit=93,file='./data_directory/LS')            !! Contains iters and residual for each PPE solve
    open(unit=94,file='./data_directory/TIME_OUT')      !! Contains particle, bubble etc counts for outputs
    
    !! Force first output before any calculations
    output_everytime = .true.    
    
    return
  end subroutine getdata
!! ------------------------------------------------------------------------------------------------
  subroutine check_boundaries
    !! This subroutine checks the boundary conditions read in from
    !! IBOUND for Jack's boundary condition framework. 
    real(rkind),dimension(2) :: tmp
    real(rkind) :: error_mag,tmp1,tmp2,p_dist
    integer(ikind) :: i,j,n_fail,ibm1,ibp1
    logical :: b_check_pass

    ! check 1: node(i) + edge(i) = node(i+1) ??         [Connectivity]
    b_check_pass = .true.
    n_fail = 0
    do i=1,nb_patches
       if(i==nb_patches) then
          j=1
       else
          j=i+1
       end if
       tmp(:) = b_node(i,:) + b_edge(i,:) - b_node(j,:)
       error_mag = sqrt(dot_product(tmp,tmp))
       write(6,*) i,j,tmp
       if(error_mag>three*vsmall)then       
          n_fail = n_fail + 1
          write(6,*) "Check 1 failed: b_nodes",i,"and ",j,"not connected"
          b_check_pass = .false.
       end if
    end do
    if(b_check_pass)then
       write(6,*) "Check 1 passed: all nodes and edges connected within tolerance"
    end if

    !check 2: sum(edges) = 0 ??                [a closed loop]
    tmp(:) = zero
    do i=1,nb_patches
       tmp(:) = tmp(:) + b_edge(i,:)
       write(6,*) tmp(:)
    end do
    error_mag = sqrt(dot_product(tmp,tmp))
!    if(error_mag>1.1d0*vsmall)then
    if(tmp(1)>three*vsmall.or.tmp(2)>three*vsmall)then
       n_fail = n_fail + 1
       write(6,*) "Check 2 failed: edges sum to",tmp,vsmall
    else
       write(6,*) "Check 2 passed: edges sum to zero"
    end if

    ! checks 3: periodic patches matching ??        [length of periodic patches]
    b_check_pass = .true.
    do i=1,nb_patches
       if(b_type(i)==2) then !periodic
          tmp = b_edge(i,:) + b_edge(b_periodic_parent(i),:)
          error_mag = sqrt(dot_product(tmp,tmp))
          if(abs(error_mag)>1.0d-10)then
             n_fail = n_fail + 1
             write(6,*) "Check 3 failed: periodic patches ",i,b_periodic_parent(i),"different length/direction"
             write(6,*) error_mag
             b_check_pass = .false.
          end if
       end if
    end do
    if(b_check_pass)then
       write(6,*) "Check 3 passed: all (if any) periodic patches well matched"
    end if
    if(n_fail/=0)then
       stop
    end if

    ! check 4: calculate angles at corners...           [sum of internal angles]
    allocate(b_corner_angle(nb_patches))
    error_mag = zero
    write(6,*) "For information, the corner angles are:"
    do i=1,nb_patches  ! for each corner
       ibm1 = mod(i+nb_patches-2,nb_patches)+1
       ibp1 = mod(i,nb_patches)+1
       ! find the perpendicular distance (signed) between edge ibm1 and node ibp1
       p_dist = b_edge(ibm1,2)*b_node(ibp1,1) - b_edge(ibm1,1)*b_node(ibp1,2) + &
         b_node(i,1)*b_node(ibm1,2) - b_node(i,2)*b_node(ibm1,1)
       p_dist = -one*p_dist/sqrt(dot_product(b_edge(ibm1,:),b_edge(ibm1,:)))
       ! find the magnitude of edge i
       tmp1 = sqrt(dot_product(b_edge(i,:),b_edge(i,:)))
       ! find the angle
       if(p_dist<=zero.and.abs(p_dist)>abs(tmp1))then
          p_dist = -tmp1
       else if(p_dist>zero.and.abs(p_dist)>abs(tmp1))then
          p_dist = tmp1
       end if
!       p_dist = min(p_dist,tmp1) ! if p_dist=1.0000000002(eg) tmp2=NaN, hence this line
       tmp2 = asin(p_dist/tmp1)
       b_corner_angle(i) = tmp2
       write(6,*) ibm1,i,tmp2
       error_mag = error_mag + tmp2
    end do
    write(6,*) "Check 4: Angles sum to : ",error_mag
    return
  end subroutine check_boundaries
!! ------------------------------------------------------------------------------------------------
  subroutine set_inflow_particles

    integer(ikind) :: i,j,ib,ibp1,np_inflow_esti
    real(rkind) :: p_dist,tmp,alph_inflow
    real(rkind),dimension(2) :: nrm 
    integer(ikind),dimension(:),allocatable :: p_inflow_tmp
    real(rkind),dimension(:,:),allocatable :: rp_inflow_tmp
    logical,dimension(:),allocatable :: p_inflow_mask

    ! I assume there is only 1 inflow patch - if >1, need to make changes
    ! I also assume that the order I find particles will be in the order they appear
    ! along the boundary patch... (have to be careful here)
    n_diff = 0  !! Counter to check the number of particles in domain is ~constant
    np_inflow = 0
    alph_inflow = one*dx  ! this is for cartesian particle spacing
    do ib = 1,nb_patches
       if(b_type(ib)/=3) cycle ! Only for the inflow patches
       ibp1 = mod(ib,nb_patches)+1
       ! calculate the boundary normal
       tmp = sqrt(dot_product(b_edge(ib,:),b_edge(ib,:)))
       nrm(1) = -one*b_edge(ib,2)/tmp
       nrm(2) = b_edge(ib,1)/tmp

       ! allocate p_inflow
       np_inflow_esti = 2*100*int(tmp/dx) + 2! *2 for safety!!
       allocate(p_inflow_tmp(np_inflow_esti));p_inflow_tmp(:)=0
       allocate(rp_inflow_tmp(np_inflow_esti,dims));rp_inflow_tmp(:,:)=zero
       
       do i=1,n_par_fw ! look through all particles
          ! find the perpendicular distance (signed) between particle i
          ! and boundary patch ib.
          p_dist = b_edge(ib,2)*rp(i,1) - b_edge(ib,1)*rp(i,2) + &
               b_node(ibp1,1)*b_node(ib,2) - b_node(ibp1,2)*b_node(ib,1)
          p_dist = -one*p_dist/sqrt(dot_product(b_edge(ib,:),b_edge(ib,:)))
          if(p_dist<0.0) cycle ! discard if i is the wrong side of ib
          if(p_dist>alph_inflow) cycle ! discard if it is too far from ib
          ! find position along the boundary patch which is closest to particle i
          if(abs(b_edge(ib,1))>abs(b_edge(ib,2)))then
             tmp = (rp(i,1) - b_node(ib,1) - p_dist*nrm(1))/b_edge(ib,1)
          else
             tmp = (rp(i,2) - b_node(ib,2) - p_dist*nrm(2))/b_edge(ib,2)
          end if
          if(tmp<zero.or.tmp>one) cycle  ! discard if it is not between ends of ib
          if(i<=n_par_w.and.p_dist>1.0d-10) cycle  ! if boundary, and not on the corner/ b_node, cycle

          ! if still here, it should be labelled an inflow particle!!
          np_inflow = np_inflow + 1
          p_inflow_tmp(np_inflow) = i
          rp_inflow_tmp(np_inflow,:) = rp(i,:)
       end do

       !! Look through inflow particles, and order along patch...
       allocate(p_inflow(np_inflow))
       allocate(rp_inflow(np_inflow,dims));rp_inflow(1:np_inflow,:)=rp_inflow_tmp(1:np_inflow,:)
       p_inflow(1:np_inflow) = p_inflow_tmp(1:np_inflow)
       
       deallocate(p_inflow_tmp,rp_inflow_tmp)

    end do
    return
  end subroutine set_inflow_particles
!! ------------------------------------------------------------------------------------------------
#ifdef dim3
  function init_velocity_field(x,y,z) result(vel)
     real(rkind),intent(in) :: x,y,z
#else
  function init_velocity_field(x,y) result(vel)
     real(rkind),intent(in) :: x,y
     real(rkind) :: z
#endif  
     real(rkind),dimension(dims) :: vel
     real(rkind) :: fpios,pios,beltrami_mag

#ifndef dim3
     z=one
#endif     
     
     !! This function takes in coordinates and returns and velocity field.
     !! It can be used to set initial flow states for problems like
     !! Taylor-Green, or to zero a velocity field. Just change the .true.
     !! and .false. flags.
     
     
     !! Beltrami flow - valid for a 2pi triply periodic box
     if(.false.)then
        beltrami_mag = (four*sqrt(two)/three/sqrt(three))
        fpios=five*pi/6.0d0;pios=pi/6.0d0        
        vel(1) = beltrami_mag*(sin(two*x-fpios)*cos(two*y-pios)*sin(two*z) - cos(two*z-fpios)*sin(two*x-pios)*sin(two*y))
        vel(2) = beltrami_mag*(sin(two*y-fpios)*cos(two*z-pios)*sin(two*x) - cos(two*x-fpios)*sin(two*y-pios)*sin(two*z))
#ifdef dim3        
        vel(3) = beltrami_mag*(sin(two*z-fpios)*cos(two*x-pios)*sin(two*y) - cos(two*y-fpios)*sin(two*z-pios)*sin(two*x))
#endif        
     end if     
     
     !! Taylor-Green initial velocity field
     if(.true.)then
       vel(1) = sin(two*pi/3.0)*(two/sqrt(three))*sin(x)*cos(y)*cos(z)
       vel(2) = sin(-two*pi/3.0)*(two/sqrt(three))*cos(x)*sin(y)*cos(z)       
#ifdef dim3
       vel(3) = zero     
#endif       
     end if
     
     !! Zero field
     if(.true.)then
        vel(:) = zero
     end if         
  
  end function init_velocity_field  
!! ------------------------------------------------------------------------------------------------  
end module input
