module bubble_boundaries
  !! This module contains routines to impose boundary conditions on the bubbles, and also 
  !! to create and destroy bubbles.
  use kind_parameters
  use common_parameter
  use common_2d

  implicit none
  real(rkind),dimension(2) :: norm_be
  real(rkind) :: perp_dist
  integer(ikind),allocatable,dimension(:) :: n_near_patch,n_near_corner
  integer(ikind),allocatable,dimension(:,:) :: b_neighbour,bc_neighbour
  real(rkind),allocatable,dimension(:,:,:) :: b_neigh_pos
  logical :: LDC

#ifdef bubbles
  public :: create_mirror_bubbles,mirror_velocities_b,replace_escapees_b
#endif  
contains
#ifdef bubbles
!! ------------------------------------------------------------------------------------------------
  subroutine create_mirror_bubbles
    !! Create mirrors of bubbles to apply wall or periodic BCs, directly following the approach for
    !! SPH particles
    real(rkind),dimension(2) :: trans_bp,bp,tang_be
    real(rkind) tmp,tmp_angle,m_be
    integer(ikind) :: i,j,ib,imp,imp_temp,k
    integer(ikind) :: ibp,ibm1,ibb,ibp1


    !! Insert or remove bubbles for inflow/outflow
    call insert_or_remove_bubbles
    
    !! Create boxes of bubbles close to each boundary patch (except non-boundaries)
    allocate(n_near_patch(nb_patches),n_near_corner(nb_patches))
    allocate(b_neighbour(nb_patches,n_bub),bc_neighbour(nb_patches,n_bub))
    allocate(b_neigh_pos(nb_patches,n_bub,2))
    call create_boundary_boxes_b
    
    !! Allocate arrays for mirror-parent relations
    allocate(irelation_b(n_bub+1:5*n_bub))
    allocate(vrelation_b(n_bub+1:5*n_bub))
    imp = 0 ! Initially there are no mirror bubbles

    !! Loop over each patch, and create mirrors as required
    do ib = 1,nb_patches
       ibp1 = mod(ib,nb_patches)+1

       ! OPTION 0: A non-boundary, do nothing
       if(b_type(ib)==0) cycle 

       ! OPTION 1: A wall boundary
       if(b_type(ib)==1) then
          norm_be(:) = bound_norm(ib)
          do j=1,n_near_patch(ib) ! loop over bubbles near patch ib
             i=b_neighbour(ib,j)
             bp = b_neigh_pos(ib,j,:) ! the position of bubble relative to boundary
             if(bp(2)<=vsmall) cycle  ! if wrong side of bound

             ! if the corner angle is less than 0 and a mirror of i inteferes with corner ib
             ! allow 0.1dx so mirror bubbles don't come too close (only for one bit of corner)
             if(b_corner_angle(ib)<zero)then
                tmp_angle = 0.5*(pi + b_corner_angle(ib))
                m_be = sqrt(dot_product(b_edge(ib,:),b_edge(ib,:)))
                if(bp(2)+0.1d0*dx>=bp(1)*tan(tmp_angle)*m_be) cycle
             end if
             ! if the corner angle is less than 0 and a mirror of i inteferes with corner ibp1
             if(b_corner_angle(ibp1)<zero)then
                tmp_angle = 0.5*(pi + b_corner_angle(ibp1))
                m_be = sqrt(dot_product(b_edge(ib,:),b_edge(ib,:)))
                if(bp(2)>=(one-bp(1))*tan(tmp_angle)*m_be) cycle
             end if

             ! if still in loop, i needs a mirror k
             imp = imp + 1
             k = n_bub + imp   ! index of the mirror bubble
             irelation_b(k) = i ! parent of the mirror
             vrelation_b(k) = ib !mirror-parent velocity relationship

             ! position of the mirror
             rb(k,1:2) = rb(i,1:2) - two*bp(2)*norm_be(:)

             ! velocity of the mirror is *-1.0
             tang_be(1)=norm_be(2);tang_be(2)=-norm_be(1)
             ub(k,1:2) = -ub(i,1:2)  
#ifdef wavemaker
             if(ib==nb_patches) then
                ub(k,1) = ub(k,1) + two*U_wmaker
             end if
#endif                      
#if dim3
             rb(k,3) = rb(i,3);ub(k,3)=-ub(i,3)
#endif              
          end do
       end if

       ! OPTION 2: A periodic boundary
       if(b_type(ib)==2)then
          ! we are looking for bubbles near boundary patch b_periodic_parent(ib)
          ibp = b_periodic_parent(ib)  !identify parent patch for mirror bubbles
          m_be = sqrt(dot_product(b_edge(ib,:),b_edge(ib,:)))    
          norm_be(:) = bound_norm(ib)
          ! loop over every bubble near patch ibp
          do j=1,n_near_patch(ibp)
             i=b_neighbour(ibp,j)
             bp(:) = b_neigh_pos(ibp,j,:)
             if(bp(2)<=zero) cycle ! if perp_dist negative (or 0)

             ! if still in loop, i needs a mirror k
             imp = imp + 1
             k = n_bub + imp   ! index of the mirror bubble
             irelation_b(k) = i ! parent of the mirror
             vrelation_b(k) = ib !mirror-parent velocity relationship
             
             ! position the mirror (relative to patches now, to account for non-parallel patches)
             trans_bp(:) = -one*bp(1)*b_edge(ib,:) - bp(2)*norm_be(:)
             rb(k,1:2) = b_node(ib,:) + b_edge(ib,:) + trans_bp(:)
#if dim3
             rb(k,3) = rb(i,3)
#endif 
             ! velocity of the mirror is unchanged (patches must be parallel)
             ub(k,:) = ub(i,:)
          end do
       end if

       ! OPTION 3: An inflow patch (prescribed U_in) modified!
       if(b_type(ib)==3)then
          !! DO ABSOLUTELY NOTHING FOR NOW
 
       end if

       ! OPTION 4: An outflow patch
       if(b_type(ib)==4)then
          !! DO ABSOLUTELY NOTHING AGAIN
       end if
    end do

    ! STEP 6: Create mirrors in the corners...
    do ib = 1,nb_patches !loop over the corners       
       ibm1 = mod(ib+nb_patches-2,nb_patches)+1
       
       ! OPTION 1: One of the patches is a non-boundary, do nothing
       if(b_type(ib)==0.or.b_type(ibm1)==0) cycle

       ! OPTION 2: The two patches are co-linear (cross product =0), do nothing
       tmp = b_edge(ibm1,1)*b_edge(ib,2) - b_edge(ibm1,2)*b_edge(ib,1)
       if(abs(tmp)<1d-10) cycle

       ! OPTION 3: a double wall corner
       if(b_type(ib)==1.and.b_type(ibm1)==1)then   ! a wall corner
          ! Check to see if the corner needs doing...
          if(b_corner_angle(ib)<=zero) cycle

          do j=1,n_near_corner(ib) ! loop over paticles near corner
             i = bc_neighbour(ib,j)
             
             trans_bp(:) = rb(i,1:2)-b_node(ib,:)
             perp_dist = sqrt(dot_product(trans_bp,trans_bp))

             ! create a mirror bubble
             imp = imp + 1
             k = n_bub + imp
             irelation_b(k) = i ! parent of the mirror
             vrelation_b(k) = -ib !mirror-parent velocity relationship

             ! position the mirror
             rb(k,1:2) = rb(i,1:2) - two*trans_bp(:)
#if dim3
             rb(k,3) = rb(i,3)
#endif 

             ub(k,:) = ub(i,:)              ! don't reverse the velocity
#ifdef wavemaker
             if(ib==1) then !First corner
                ub(k,1) = two*U_wmaker - ub(i,1)
             end if
#endif               
           end do
       end if

       ! OPTION 4: A double periodic corner 
       if(b_type(ib)==2.and.b_type(ibm1)==2)then
          ibp = b_periodic_parent(ib)  !identify parent patch for mirror bubbles
          do j=1,n_near_corner(ibp) ! look for bubbles near b_node(ibp)
             i=bc_neighbour(ibp,j)

             ! create a mirror bubble
             imp = imp + 1
             k = n_bub + imp
             irelation_b(k) = i ! parent of the mirror
             vrelation_b(k) = -ib !mirror-parent velocity relationship

             ! position the mirror
             trans_bp(:) = b_node(ib,:) - b_node(ibp,:)
             rb(k,1:2) = rb(i,1:2) + trans_bp(:)
#if dim3
             rb(k,3) = rb(i,3)
#endif 
             ! don't change the velocity
             ub(k,:) = ub(i,:)
          end do
       end if

       ! OPTION 5: The two patches are of different type
       if(b_type(ib)/=b_type(ibm1))then
          ! OPTION 5a. One of the two is a wall, we will look near this first...
          if(b_type(ib)==1.or.b_type(ibm1)==1)then 
             if(b_type(ibm1)==1)then
                ibb = ibm1   !! Ibb is the wall patch
             else
                ibb = ib
             end if
             imp_temp = imp
             do j=n_bub+1,n_bub + imp_temp  ! look through existing mirror bubbles near wall
                perp_dist = bound_dist(j,ibb)
                if(perp_dist>sup_size+vsmall) cycle  ! cycle if it's not near the patch
                if(perp_dist<vsmall) cycle ! cycle if it's actually aligned with the patch
                norm_be(:) = bound_norm(ibb)  ! find the normal...
                tmp = position_along_patch(j,ibb)            
                ! is it off the end (the correct end) of the patch?
                if(tmp<=zero.and.ibb==ib.or.tmp>=one.and.ibb==ibm1) then
                   imp = imp + 1
                   k = n_bub + imp
                   ! parent of the mirror is the parent of the mirror in which it is mirrored!
                   i = irelation_b(j)
                   irelation_b(k) = i
                   vrelation_b(k) = -j !mirror-parent velocity relationship

                   ! position the mirror
                   rb(k,1:2) = rb(j,1:2) - two*perp_dist*norm_be(:)

                   tang_be(1)=norm_be(2);tang_be(2)=-norm_be(1)
                   ub(k,1:2) =-ub(j,1:2)  
#if dim3
                   rb(k,3) = rb(i,3);ub(k,3)=-ub(i,3)
#endif 
                end if
             end do
          else
             ! OPTION 5b: Neither of the two is a wall...(TBCompleted)  
          end if
       end if
    end do
    
    !! Make mirrors in third dimension - periodic only. 
    !! N.B. The third dimension is done last, and includes mirrors of mirror bubbles created for
    !! the first two dimensions.
#if dim3    
    n_bub_m = n_bub + imp
    do i=1,n_bub_m
       if(b_inbin(i)) cycle
       if(rb(i,3)<=zb_min+sup_size)then !! It's near the minZ bound
          imp = imp + 1
          k = n_bub + imp
          if(i>n_bub)then
             irelation_b(k) = irelation_b(i)
             if(vrelation_b(i)<-nb_patches) then  !! wall-periodic corner
                vrelation_b(k)=9999
             else if(vrelation_b(i)<0) then !! wall-wall or periodic-periodic corner
                vrelation_b(k)=9997
             else if(vrelation_b(i)<=nb_patches) then !! periodic, wall or inflow patch
                if(b_type(vrelation_b(i))==1)then ! wall
                   vrelation_b(k)=9996
                else
                   vrelation_b(k)=9995
                end if
             end if
          else
             irelation_b(k) = i
             vrelation_b(k)=9998
          end if
#ifdef zwall
          ub(k,:) = -ub(i,:)    
          rb(k,:)=rb(i,:);rb(k,3)= two*zb_min - rb(i,3)
#else          
          ub(k,:) = ub(i,:)           
          rb(k,:)=rb(i,:);rb(k,3)=rb(i,3) + zb_max-zb_min
#endif          
       end if
       if(rb(i,3)>=zb_max-sup_size)then !! It's near the maxZ bound
          imp = imp + 1
          k = n_bub + imp
          if(i>n_bub)then
             irelation_b(k) = irelation_b(i)
             if(vrelation_b(i)<-nb_patches) then  !! wall-periodic corner
                vrelation_b(k)=9999
             else if(vrelation_b(i)<0) then !! wall-wall or periodic-periodic corner
                vrelation_b(k)=9997
             else if(vrelation_b(i)<=nb_patches) then !! periodic, wall or inflow patch
                if(b_type(vrelation_b(i))==1)then ! wall
                   vrelation_b(k)=9996
                else
                   vrelation_b(k)=9995
                end if
             end if
          else
             irelation_b(k) = i
             vrelation_b(k)=9998
          end if
#ifdef zwall
          ub(k,:) = -ub(i,:)    
          rb(k,:)=rb(i,:);rb(k,3)= two*zb_max - rb(i,3)
#else          
          ub(k,:) = ub(i,:) 
          rb(k,:)=rb(i,:);rb(k,3)=rb(i,3) - zb_max+zb_min
#endif          
       end if
    end do    
#endif    
    
    
    !! Update np and nmirror
    n_bub_m = n_bub + imp
    
    !! Copy bubble sizes and other properties from parents
    !$omp parallel do private(i)
    do j=n_bub+1,n_bub_m
       i=irelation_b(j)
       radb(j) = radb(i)    
    end do
    !$omp end parallel do

    deallocate(n_near_patch,n_near_corner)
    deallocate(b_neighbour,bc_neighbour)
    deallocate(b_neigh_pos)

    return
  end subroutine create_mirror_bubbles
!! ------------------------------------------------------------------------------------------------
  subroutine replace_escapees_b
    !! Check whether any bubbles have escaped. For bubbles which have passed through a wall,
    !! reflect their position and velocity. For bubbles which have passed through a periodic 
    !! boundary, insert them in the corresponding periodic parent patch.
    real(rkind),dimension(2) :: trans_bp
    real(rkind) tmp,m_be,u_norm
    integer(ikind) :: i,ib,ibp,n_escapees

    n_escapees = 0
    do ib = 1,nb_patches
       if(b_type(ib)==0) cycle  ! if it's a non-boundary, do nothing
       if(b_type(ib)==1.or.b_type(ib)==2) then   ! it's a wall or periodic
          ! loop over every bubble
          norm_be(:) = bound_norm(ib)  !! N.B. wall normals point into fluid
          ibp = b_periodic_parent(ib)
          m_be = sqrt(dot_product(b_edge(ib,:),b_edge(ib,:)))
          
          
          !$omp parallel do private(perp_dist,tmp,trans_bp,u_norm) reduction(+:n_escapees)
          do i=1,n_bub
             if(.not.b_inbin(i))then
                ! find the perpendicular distance between i and the boundary edge
                perp_dist = bound_dist(i,ib) 
   
                ! for bubbles close but on the wrong side of the boundary
                if(perp_dist<zero.and.abs(perp_dist)<=sup_size+vsmall) then  
                   tmp = position_along_patch(i,ib)
                   ! if i is close to the boundary patch but beyond the ends, cycle
                   ! it needs replacing, are we at a wall or a periodic?
                   if(b_type(ib)==1) then ! WALL
                      if(tmp>=zero.and.tmp<=one) then
                         ! put bubble back in (reflect in boundary)
                         rb(i,1:2) = rb(i,1:2) - two*perp_dist*norm_be(:)
                         u_norm=dot_product(ub(i,1:2),norm_be(:))   !! normal component of velocity
                         if(u_norm<=zero)then
                            ub(i,1:2) = ub(i,1:2) + two*u_norm*norm_be(:) !! reflect normal component..
                         end if                         
                         write(6,*) "replaced bubble which had passed through wall",i,ib,rb(i,:),perp_dist
                         n_escapees = n_escapees + 1
                      end if
                   else ! PERIODIC
                      if(tmp<zero.or.tmp>one) then  !! if unexpectedly about to leave via a corner
                         trans_bp(:) = -one*tmp*b_edge(ibp,:) - perp_dist*bound_norm(ibp)
                         rb(i,1:2) = b_node(ibp,:) + b_edge(ibp,:) + trans_bp(:) 
                         n_escapees = n_escapees + 1               
                      else
                         ! reposition the bubble relative to the boundary patches, which
                         ! allows for non-parallel patch-parents (eg. bent Poiseuille flow)
                         trans_bp(:) = -one*tmp*b_edge(ibp,:) - perp_dist*bound_norm(ibp)
                         rb(i,1:2) = b_node(ibp,:) + b_edge(ibp,:) + trans_bp(:)
                         n_escapees = n_escapees + 1
                      end if
                   end if
                end if
             end if
          end do
          !$omp end parallel do
       end if
    end do
    
#if dim3    
    !! 3rd dimension 
    !$omp parallel do
    do i=1,n_bub
       if(.not.b_inbin(i)) then    
          if(rb(i,3)<zb_min)then
#ifdef zwall
             rb(i,3) = two*zb_min - rb(i,3);ub(i,3) = -ub(i,3)          
#else          
             rb(i,3) = rb(i,3) + zb_max - zb_min
#endif
          end if
          if(rb(i,3)>=zb_max)then
#ifdef zwall
             rb(i,3) = two*zb_max - rb(i,3);ub(i,3) = -ub(i,3)
#else
             rb(i,3) = rb(i,3) - zb_max + zb_min
#endif
          end if
       end if
    end do
    !$omp end parallel do
#endif    
    
  end subroutine replace_escapees_b
!! ------------------------------------------------------------------------------------------------
  subroutine insert_or_remove_bubbles
    !! A dummy of the inflow/outflow insertion/removal routine, followed by calculations of 
    !! entrainment rates or plume sources to insert bubbles into the simulation.
    use sphtools
    real(rkind) :: tmp,tmp2
    integer(ikind) :: i,j,ib,n_new,steps_per_nb,n_tries
    real(rkind) :: nb_per_step,nb_remainder,Rh
    real(rkind) :: Rb0,Ab0,Vb0,E_bc_avail,V_bc_avail,Rmin,Rmax
    logical :: keep_going,make_bubble
    
    
    do ib = 1,nb_patches
       if(b_type(ib)<=2) cycle ! Do nothing for non, wall or periodic patches
       if(b_type(ib)==4) then !if outflow, we will deal with it first        
          !! DO NOTHING FOR NOW
       else if(b_type(ib)==3) then ! if it's an inflow boundary, deal with it second
          !! DO NOTHING FOR NOW
       end if
    end do
    
if(.false.)then
    !! A bubble source within the flow?
    !! How many bubbles?
    nb_per_step = 1989.4*dt!6366.2*dt
    if(nb_per_step<one) then                !! If <1 bubble per step, find steps per new bubble,
       steps_per_nb = floor(one/nb_per_step)   !! and only do new bubble if mod(itime,spnb) is zero
       if(mod(itime,steps_per_nb)==0) then
          n_new = 1
       else
          n_new = 0
       end if
    else              !! Set n_new = floor or ceiling of nb_per_step, with choice set randomly (weighted by nb_remainder)
       n_new = floor(nb_per_step)
       nb_remainder = nb_per_step - floor(nb_per_step)
       call random_number(tmp2)
       if(tmp2<nb_remainder) n_new = n_new + 1
    end if

    if(n_new/=0)then
       do i=1,n_new
          !! Create a new bubble either by increasing n_bub or using a free index
          if(b_nfree==0) then
             n_bub = n_bub + 1
             j=n_bub !! j is index of new bubble
          else
             j = b_freeindices(b_nfree)
             b_nfree = b_nfree - 1
          end if
          b_inbin(j) = .false.
!          call random_number(tmp);call random_number(tmp2)
          tmp = rand();tmp2=rand()
          tmp = tmp*0.02;tmp2 = tmp2*2.0*pi ! tmp and tmp2 are radius and angle
          rb(j,1) = tmp*cos(tmp2)  !! Bubbles are created in circular disk radius 0.02, center -0.1,0,0
#ifdef dim3
          rb(j,3) = tmp*sin(tmp2);ub(j,3)=zero
#endif         
          rb(j,2) = 0.09d0 ! 0.09

          keep_going = .true.
          do while(keep_going)
             tmp = rand();tmp2=rand()
             radb(j) = min(0.001d0 + 0.001*cos(two*pi*tmp2)*sqrt(-two*log(tmp)),dx)
             if(radb(j)>=0.0001) keep_going = .false.
          end do

          
          ub(j,1) = zero;ub(j,2) = 0.2d0!0.00063333d0!0.2d0     
          bubble_EoL(j) = 0
          bubble_age(j) = zero
          bubble_LE(j) = 1.0d10
           
       end do
    end if    
       
!    write(6,*) "added ",n_new,"bubbles this step"
else if(.true.)then
!! Bubbles generated by turbulence

#ifdef srs
    !! Loop over all liquid particles
    do i=1,n_par_fw
       if(inbin(i)) cycle
       !! Check whether the particle is at the free surface and the FS isn't saturated
       if(n_surf(i)==1.and.a0(i)>0.7d0) then
          !! Check whether the SGS dissipation rate is sufficient
          if(eps_srs(i)>0.2d0) then 
             
             !! Set the amount of energy available for bubble creation in this time-step
             E_bc_avail = 0.01*beta*a0(i)*eps_srs(i)*dv*dt
                
             !! Set the volume available to be filled with bubbles
             V_bc_avail = (one/ad_7)*(one/0.7 - one/a0(i))
             if(V_bc_avail>zero)then
                
             !! Loop as long as there's energy to make more bubbles
             keep_going = .true.
             n_new = 0
             n_tries = 0
             do while(keep_going)
                make_bubble = .true.
                n_tries = n_tries + 1
                
                !! Set max and min bubbles sizes (made dimensionless by Hinze scale)
                !! N.B. Hinze scale appears here (based on We & Re, set in input module)
                !! but only for making outputs dimensionless. It does not enter the physics
                !! of the model.
                Rh = Rhinze
                Rmax = dx/Rh!4.0!0.4*dx/Rh  !! Set max bubble size = particle spacing
                Rmin = 5.0/1.0d2!0.1       !! set min bubble size relative to max bubble size
                
                !! Randomly asign a bubble size between Rmin and Rmax
                !! Choice of pdf for random number
!                Rb0 = deane_and_stokes_spectrum_rand(Rmin,Rmax)
                Rb0 = flat_spectrum_rand(Rmin,Rmax)                     
!                Rb0 = linear_spectrum_rand(Rmin,Rmax)
                Rb0 = Rb0*Rh       !! Rescale by Rh to give absolute size
               
                !! Bubble surface area and volume
                Ab0 = 4.0*pi*Rb0*Rb0         
                Vb0 = 4.0*pi*Rb0*Rb0*Rb0/3.0  
                   
                !! Only make the bubble if there is enough energy remaining
                if(E_bc_avail<Ab0/We)then
                   nb_remainder = E_bc_avail/(Ab0/We) !! A bit of randomness so we don't penalize large bubbles
                   if(rand()>nb_remainder) then
                      make_bubble = .false.
                   end if
                end if
                   
                !! Only make the bubble if it won't result in too small liquid volume fraction
                if(V_bc_avail<=Vb0) then
                   nb_remainder = V_bc_avail/Vb0 !! A bit of randomness to sometimes allow bigger bubbles
                   if(rand()>nb_remainder)then
                      make_bubble = .false.
                   end if
                end if
                                                
                !! If checks were passed, make the bubble                               
                if(make_bubble) then  
                   n_new = n_new + 1              
                   E_bc_avail = E_bc_avail - Ab0/We  !! Reduce available energy                      
                   V_bc_avail = V_bc_avail - Vb0     !! Reduce available volume
                   call make_new_bubble(i,Rb0)
                end if

                !! Check whether we have capacity for more bubbles 
                if(E_bc_avail<=4.0*pi*Rh*Rh*Rmin*Rmin/We) then  !! Energy 
                   keep_going = .false.
                end if
                if(V_bc_avail<=zero) then  !! Volume
                   keep_going = .false.
                end if
                if(n_new>=10.or.n_tries>20) then  !! Limit number of new bubbles and number of tries
                   keep_going = .false.
                end if
                if(n_tries-n_new>5)then !! if X failed attempts
                   keep_going = .false.
                end if
                                      
             end do
             end if
          end if
       end if
    end do
    
    if(itime.eq.10) then
    do i=1,n_par_fw
       !! Check whether the particle is at the free surface and the FS isn't saturated
       if(n_surf(i)==1) then       
          j=i
       end if
    end do

!    call make_new_bubble(j,Rhinze)    
    
    end if
    
#endif    
end if 
      
     
      
    return
  end subroutine insert_or_remove_bubbles
!! ------------------------------------------------------------------------------------------------
  subroutine make_new_bubble_at_position(x,y,z,radius)
    !! Create a new bubble of a specified size and position
    real(rkind),intent(in) :: x,y,z,radius    !! Size of bubble
    integer(ikind) :: j
    real(rkind) :: tmp
    
    !! Use a new index or a free index?
    if(b_nfree==0) then 
       n_bub = n_bub + 1
       j=n_bub
    else
       j = b_freeindices(b_nfree) 
       b_nfree = b_nfree - 1
    end if
    
    !! It is not in the bin
    b_inbin(j) = .false.
    
    !! Copy position and velocity
    rb(j,1) = x;rb(j,2) = y;rb(j,3) = z;ub(j,:)=zero
    
    !! allocate size
    radb(j) = radius
    
    !! Set some age flags
    bubble_EoL(j) = 0        !! New bubble is not in End-of-Life
    bubble_age(j) = zero    !! New bubble has age zero
    bubble_LE(j) = 1.0d10    !! Life expectancy
    
    !! The bubble is initially undeformed::
    deformation_distance(j) = zero
    mean_dissrate(j) = zero
    
    !! Set the liquid velocity at bubble location (for use in virtual mass calculation)
    u_l2b(j,:) = zero
  
    return
  end subroutine make_new_bubble_at_position
!! ------------------------------------------------------------------------------------------------  
  subroutine make_new_bubble(i,radius)
    use bubble_acoustics
    !! Make a new bubble at SPH particle i, with specified size.
    integer(ikind),intent(in) :: i      !! Index of SPH particle "creating" this bubble
    real(rkind),intent(in) :: radius    !! Size of bubble
    integer(ikind) :: j
    real(rkind) :: tmp,tbub
    
    !! Use a new index or a free index?
    if(b_nfree==0) then 
       n_bub = n_bub + 1
       j=n_bub
    else
       j = b_freeindices(b_nfree) 
       b_nfree = b_nfree - 1
    end if
    
    !! It is not in the bin
    b_inbin(j) = .false.
    
    !! Copy position and velocity
    rb(j,:) = rp(i,:);ub(j,:) = up(i,:)
    
    !! allocate size
    radb(j) = radius
    
    !! Set some age flags
    bubble_EoL(j) = 0        !! New bubble is not in End-of-Life
    bubble_age(j) = dt*rand() !! Bubble was created sometime in last SPH step
    bubble_LE(j) = 1.0d10    !! Life expectancy
    
    !! The bubble is initially undeformed::
    deformation_distance(j) = zero
    mean_dissrate(j) = zero
    
    !! Add some random noise to the position ::
    tmp = rand();rb(j,1) = rb(j,1) + dx*(tmp-half)
    tmp = rand();rb(j,2) = rb(j,2) + dx*(tmp-half)
#ifdef dim3
    tmp = rand();rb(j,3) = rb(j,3) + dx*(tmp-half)                                            
#endif    

    !! Add some random noise to the velocity
!    tmp = rand();ub(j,1) = ub(j,1)*(one + 0.1*(tmp-half))
!    tmp = rand();ub(j,2) = ub(j,2)*(one + 0.1*(tmp-half))
#ifdef dim3
!    tmp = rand();ub(j,3) = ub(j,3)*(one + 0.1*(tmp-half))
#endif    

    !! Set the liquid velocity at bubble location (for use in virtual mass calculation)
    u_l2b(j,:) = up(i,:)
    
#ifdef acoustics
    tbub = time - bubble_age(j)
    call add_pressure_from_bubble(radius,rb(j,:),tbub)
#endif    
  
    return
  end subroutine make_new_bubble    
!! ------------------------------------------------------------------------------------------------
  subroutine create_boundary_boxes_b
    !! Loops through all boundary patches and corners, and builds lists of bubbles close to each
    !! patch and corner, and their position relative to that patch.
    integer(ikind) :: i,ib
    real(rkind) :: tmp

    n_near_patch(:) = 0
    n_near_corner(:) = 0
    do ib = 1,nb_patches
       if(b_type(ib)==0) cycle ! if it's a non-boundary, i don't care
       norm_be(:) = bound_norm(ib)
       do i=1,n_bub
          if(b_inbin(i)) cycle       
          ! find the perpendicular distance between i and the boundary edge
          perp_dist = bound_dist(i,ib)
          ! if perp_dist too big, cycle
          if(abs(perp_dist)>sup_size+vsmall) cycle
          ! find the normal to the boundary (pointing inwards (fluid-wards))
          tmp = position_along_patch(i,ib)
          ! if i is close to the line but beyond the ends, cycle
          if(tmp<zero.or.tmp>one) cycle
          n_near_patch(ib) = n_near_patch(ib) + 1 ! how many in box?
          b_neighbour(ib,n_near_patch(ib)) = i   ! index of this within box
          b_neigh_pos(ib,n_near_patch(ib),1) = tmp  ! position relative to bnode
          b_neigh_pos(ib,n_near_patch(ib),2) = perp_dist ! position relative to bnode
       end do
       do i=1,n_bub
          perp_dist = sqrt((rb(i,1) - b_node(ib,1))**2 + (rb(i,2)-b_node(ib,2))**2)
          if(perp_dist>sup_size+vsmall) cycle
          ! it is near a corner!!
          n_near_corner(ib) = n_near_corner(ib) + 1
          bc_neighbour(ib,n_near_corner(ib)) = i
       end do
    end do

    return
  end subroutine create_boundary_boxes_b
!! ------------------------------------------------------------------------------------------------
  function bound_dist(i,ib) result(p_dist)
    integer(ikind),intent(in) :: i,ib
    real(rkind) :: p_dist
    integer(ikind) :: ibp1

    ! the ib+1'th patch is 1 if ib = nb_patches
    ibp1 = mod(ib,nb_patches)+1

    ! find the perpendicular distance (signed) between bubble i
    ! and boundary patch ib.
    p_dist = b_edge(ib,2)*rb(i,1) - b_edge(ib,1)*rb(i,2) + &
         b_node(ibp1,1)*b_node(ib,2) - b_node(ibp1,2)*b_node(ib,1)
    p_dist = -one*p_dist/sqrt(dot_product(b_edge(ib,:),b_edge(ib,:)))

  end function bound_dist
!! ------------------------------------------------------------------------------------------------
  function bound_norm(ib) result(nrm)
    ! returns the unit normal vector to the boundary edge
    integer(ikind),intent(in) :: ib
    real(rkind) :: tmp1
    real(rkind),dimension(2) :: nrm

    tmp1 = sqrt(dot_product(b_edge(ib,:),b_edge(ib,:)))  !! edge vector magnitude

    nrm(1) = -one*b_edge(ib,2)/tmp1           !! unit normal
    nrm(2) = b_edge(ib,1)/tmp1
  end function bound_norm
!! ------------------------------------------------------------------------------------------------
  function position_along_patch(i,ib) result(tmp)
    !! Returns the distance along the boundary patch ib of bubble i:
    !! rb(i,:) = b_node(ib,:) + tmp*b_edge(ib,:) + perp_dist*norm_be(:)
    !!
    !! N.B. This function assumes that perp_dist and norm_be have already been obtained for
    !! patch ib.
    integer(ikind),intent(in) :: i,ib
    real(rkind) :: tmp
    
    ! find position along the boundary patch which is closest to bubble i
    if(abs(b_edge(ib,1))>abs(b_edge(ib,2)))then
       tmp = (rb(i,1) - b_node(ib,1) - perp_dist*norm_be(1))/b_edge(ib,1)
    else
       tmp = (rb(i,2) - b_node(ib,2) - perp_dist*norm_be(2))/b_edge(ib,2)
    end if
             
  end function position_along_patch
!! ------------------------------------------------------------------------------------------------
   subroutine mirror_velocities_b
   ! set mirror velocities (ie. after viscous, before ppe setup)
     integer(ikind) i,j,ib,ibm1
     real(rkind),dimension(dims) :: rij
     real(rkind),dimension(2) :: tang_be
     real(rkind) :: tmp
     
     !! j is the mirror bubble under consideration
     !! i is the parent of j
     
     do j=n_bub + 1,n_bub_m
        i=irelation_b(j)
                
        !! NON-CORNER
        if(vrelation_b(j)>0.and.vrelation_b(j)<999) then 
        
           !! WALL
           if(abs(b_type(vrelation_b(j)))==1)then 
              rij = half*(rb(i,:)-rb(j,:));tmp=sqrt(dot_product(rij,rij))
              tang_be(1)=rij(2)/tmp;tang_be(2)=-rij(1)/tmp
              ub(j,1:2) =  - ub(i,1:2)
#if dim3
              ub(j,3)=-ub(i,3)
#endif 
#ifdef wavemaker
              if(vrelation_b(j)==nb_patches) then
                 ub(j,1) = ub(i,1) + two*U_wmaker
              end if
#endif
           else if(b_type(vrelation_b(j))==2) then ! periodic
              ub(j,:) = ub(i,:)  ! no change in velocity (assuming parallel periodic patches)
           else if(b_type(vrelation_b(j))==3) then ! inflow, u* is U_inflow
              ub(j,:) = ub(i,:)        
           else if(b_type(vrelation_b(j))==4)then ! outflow
              ub(j,:) = ub(i,:) !! No change in velocity (du/dn = 0)
           end if
        end if
        
        !! CORNER (wall-wall or periodic-periodic)
        if(vrelation_b(j)<0.and.vrelation_b(j)>-n_bub)then
        
           ib = -vrelation_b(j)
           ibm1 = mod(ib+nb_patches-2,nb_patches)+1      

           !! WALL WALL CORNER
           if(b_type(ib)==1.and.b_type(ibm1)==1) then 

              ub(j,:) = ub(i,:)  ! don't reverse velocity
#ifdef wavemaker
              if(ib==1) then !! First corner
                 ub(j,1) = two*U_wmaker - ub(i,1)
              endif
#endif

           !! PERIODIC PERIODIC CORNER
           else if(b_type(ib)==2.and.b_type(ibm1)==2) then 
              ub(j,:) = ub(i,:)   ! assuming parallel patches
           end if
           
        !! WALL-OTHER CORNER   
        else if(vrelation_b(j)<-n_bub) then ! wall-periodic or wall-inflow corner
           ub(j,:) = -ub(-vrelation_b(j),:) 
        end if
        
#ifdef dim3        
        !! Third dimension mirrors
        if(vrelation_b(j)==9999)then ! parent is wall-periodic corner  
#ifdef zwall
           ub(j,:) = ub(i,:)
#else        
           ub(j,:) = -ub(i,:)
#endif           
        end if
        if(vrelation_b(j)==9998)then ! parent is fluid
#ifdef zwall
           ub(j,:) = -ub(i,:)
#else
           ub(j,:) = ub(i,:)
#endif           
        end if
        if(vrelation_b(j)==9997)then ! parent is wall-wall or periodic-periodic corner
#ifdef zwall
           ub(j,:) = -ub(i,:)
#else        
           ub(j,:) = ub(i,:)
#endif           
        end if
        if(vrelation_b(j)==9996)then ! parent is wall    
#ifdef zwall
           ub(j,:) = ub(i,:)
#else
           ub(j,:) = -ub(i,:)
#endif           
        end if
        if(vrelation_b(j)==9995)then ! parent is periodic,inflow  
#ifdef zwall
           ub(j,:) = -ub(i,:)
#else
           ub(j,:) = ub(i,:)
#endif           
        end if        
#endif
     end do

   return
   end subroutine mirror_velocities_b
#endif
!! ------------------------------------------------------------------------------------------------   
end module bubble_boundaries
