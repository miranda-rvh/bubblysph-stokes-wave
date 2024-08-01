module mirror_boundaries
  use kind_parameters
  use common_parameter
  use common_2d

  implicit none
  real(rkind),dimension(2) :: norm_be
  real(rkind) :: perp_dist
  integer(ikind),allocatable,dimension(:) :: n_near_patch,n_near_corner
  integer(ikind),allocatable,dimension(:,:) :: b_neighbour,bc_neighbour
  real(rkind),allocatable,dimension(:,:,:) :: b_neigh_pos

  public :: create_mirror_particles,angle_between_patches,mirror_velocities 
contains

! intype=1 gives uniform inflow/outflow. intype=2 gives SS Poiseuille inflow/outflow profiles.
#define intype 1 
  subroutine create_mirror_particles
!! -- NB:
!! For a mirror particle j with true-parent i, irelation(j)=i

!! The flag vrelation(j) indicates how the mirroring took place:
!! 0<vrelation(j)<=nb_patches 
!! If j is a wall mirror, vrelation(j) = ib (index of wall patch j is near)
!! If j is a periodic mirror, vrelation(j) = ib (the index of the periodic patch j is near)
!! If j is an inflow mirror, vrelation(j) = ib (index of inflow patch j is near)

!! -nb_patches<=vrelation(j)<0 
!! If j is a double wall corner, vrelation(j) = -ib (the index of the boundary node (corner) it is near)
!! If j is a double periodic corner, vrelation(j) = -ib (the index of the boundary node (corner) it is near)

!! If j is a wall-periodic corner, vrelation(j) = -i  (here i=index of mirror parent)

!! The variable imp counts up the mirror particles as we progress through different types of boundary patch
    real(rkind),dimension(2) :: trans_bp,bp,tang_be
    real(rkind) tmp,tmp_angle,m_be
    integer(ikind) :: i,j,ib,imp,imp_temp,k,l
    integer(ikind) :: ibp,ibm1,ibb,ibp1

    ! STEP 1: Replace any particles which have escaped
    call replace_escapees

    ! STEP 2: Insert or remove particles for inflow/outflow
    call insert_or_remove_particles
    !! NOTE: I will modify this subroutine for inflow/outflow in due course...
    
    ! STEP 3: Create boxes of particles close to each boundary patch (except non-boundaries)
    allocate(n_near_patch(nb_patches),n_near_corner(nb_patches))
    allocate(b_neighbour(nb_patches,n_par_fw),bc_neighbour(nb_patches,n_par_fw))
    allocate(b_neigh_pos(nb_patches,n_par_fw,2))
    call create_boundary_boxes
    
    ! STEP 4: Allocate arrays for mirror-parent relations
    nmirror_esti = 2*n_par_fw  ! adjust estimate for max number of mirrors (as n_par_fw increases)
    if(n_par_fw<=1000) nmirror_esti = n_par_fw*8
    allocate(irelation(n_par_fw+1:n_par_fw+nmirror_esti))
    allocate(vrelation(n_par_fw+1:n_par_fw+nmirror_esti))
    allocate(dP_mp(n_par_fw+1:n_par_fw+nmirror_esti))
    imp = 0 ! Initially there are no mirror particles

    ! STEP 5: Loop over each patch, and create mirrors as required
    do ib = 1,nb_patches
       ibp1 = mod(ib,nb_patches)+1

       ! OPTION 0: A non-boundary, do nothing
       if(b_type(ib)==0) cycle 

       ! OPTION 1: A wall boundary
       if(b_type(ib)==1) then
          norm_be(:) = bound_norm(ib)
          do j=1,n_near_patch(ib) ! loop over particles near patch ib
             i=b_neighbour(ib,j)
             bp = b_neigh_pos(ib,j,:) ! the position of particle relative to boundary
             if(bp(2)<=vsmall) cycle  ! if wrong side of bound

             ! if the corner angle is less than 0 and a mirror of i inteferes with corner ib
             ! allow 0.1dx so mirror particles don't come too close (only for one bit of corner)
             if(b_corner_angle(ib)<zero)then
                tmp_angle = half*(pi + b_corner_angle(ib))
                m_be = sqrt(dot_product(b_edge(ib,:),b_edge(ib,:)))
                if(bp(2)+0.1d0*dx>=bp(1)*tan(tmp_angle)*m_be) cycle
             end if
             ! if the corner angle is less than 0 and a mirror of i inteferes with corner ibp1
             if(b_corner_angle(ibp1)<zero)then
                tmp_angle = half*(pi + b_corner_angle(ibp1))
                m_be = sqrt(dot_product(b_edge(ib,:),b_edge(ib,:)))
                if(bp(2)>=(one-bp(1))*tan(tmp_angle)*m_be) cycle
             end if

             ! if still in loop, i needs a mirror k
             imp = imp + 1
             k = n_par_fw + imp   ! index of the mirror particle
             irelation(k) = i ! parent of the mirror
             vrelation(k) = ib !mirror-parent velocity relationship

             ! position of the mirror
             rp(k,1:2) = rp(i,1:2) - two*bp(2)*norm_be(:)

#ifdef slip
             up(k,1:2) = reverse_normal_velocity(norm_be(1),norm_be(2),up(i,1),up(i,2))
#else
             ! velocity of the mirror is *-1.0
             tang_be(1)=norm_be(2);tang_be(2)=-norm_be(1)
             up(k,1:2) = -up(i,1:2)           
#endif

#ifdef wavemaker
             if(ib==nb_patches) then !! If wavemaker, and the final patch
                up(k,1) = up(k,1) + two*U_wmaker  !! Adjust the velocity of the mirror
             end if
#endif
             
#if dim3
             rp(k,3) = rp(i,3);up(k,3)=-up(i,3)
#endif              
             ! set the pressure difference between mirror and parent
             dP_mp(k) = dot_product(grav,rp(k,:)-rp(i,:))
          end do
       end if

       ! OPTION 2: A periodic boundary
       if(b_type(ib)==2)then
          ! we are looking for particles near boundary patch b_periodic_parent(ib)
          ibp = b_periodic_parent(ib)  !identify parent patch for mirror particles
          m_be = sqrt(dot_product(b_edge(ib,:),b_edge(ib,:)))    
          norm_be(:) = bound_norm(ib)
          ! loop over every particle near patch ibp
          do j=1,n_near_patch(ibp)
             i=b_neighbour(ibp,j)
             bp(:) = b_neigh_pos(ibp,j,:)
             if(bp(2)<=zero) cycle ! if perp_dist negative (or 0)

             ! if still in loop, i needs a mirror k
             imp = imp + 1
             k = n_par_fw + imp   ! index of the mirror particle
             irelation(k) = i ! parent of the mirror
             vrelation(k) = ib !mirror-parent velocity relationship
             
             ! position the mirror (relative to patches now, to account for non-parallel patches)
             trans_bp(:) = -one*bp(1)*b_edge(ib,:) - bp(2)*norm_be(:)
             rp(k,1:2) = b_node(ib,:) + b_edge(ib,:) + trans_bp(:)
#if dim3
             rp(k,3) = rp(i,3)
#endif 
             ! velocity of the mirror is unchanged (patches must be parallel)
             up(k,:) = up(i,:)

             ! set the pressure difference between mirror and parent
             dP_mp(k) = zero
          end do
       end if

       ! OPTION 3: An inflow patch (prescribed U_in) modified!
       if(b_type(ib)==3)then
          ! find the normal to the boundary (pointing inwards (fluid-wards))
          norm_be(:) = bound_norm(ib)
          do j=1,np_inflow  ! For every inflow particle (in the domain)
             i=p_inflow(j) 

             do l=1,3 ! create 3 mirrors, with spacings dx               
                imp = imp + 1
                k = n_par_fw + imp
                irelation(k) = i ! parent of the mirror
                vrelation(k) = ib !mirror-parent velocity relationship

                ! position of the mirror (replicate band sup_size wide near patch)
                rp(k,1:2) = rp(i,1:2) - real(l)*dx*norm_be(:)

                ! Velocity of mirror particle
                up(k,1:2) = norm_be(:)*U_inflow
#if dim3
                rp(k,3) = rp(i,3);up(k,3)=zero ! No velocity at inflow in 3rd dimension
#endif 
                ! set the pressure difference between mirror and parent
                dP_mp(k) = dot_product(grav,rp(k,:)-rp(i,:))
             end do
          end do
          
          do j=1,n_near_patch(ib)
             i=b_neighbour(ib,j)
             bp = b_neigh_pos(ib,j,:) ! the position of particle relative to boundary
             if(bp(2)<=vsmall) cycle  ! if wrong side of bound
             !! if still in loop, flag as noupdate
           
          end do
       end if

       ! OPTION 4: An outflow patch
       if(b_type(ib)==4)then
          norm_be(:) = bound_norm(ib)
          do j=1,n_near_patch(ib) ! loop over particles near patch ib
             i=b_neighbour(ib,j)
             bp = b_neigh_pos(ib,j,:) ! the position of particle relative to boundary
             if(bp(2)<=vsmall) cycle  ! if wrong side of bound
             ! if still in loop, i needs a mirror k
             imp = imp + 1
             k = n_par_fw + imp   ! index of the mirror particle
             irelation(k) = i ! parent of the mirror
             vrelation(k) = ib !mirror-parent velocity relationship

             rp(k,1:2) = rp(i,1:2) - sup_size*norm_be(:)!two*bp(2)*norm_be(:) !
!             rp(k,1:2) = rp(i,1:2) - two*bp(2)*norm_be(:) !             
             up(k,:) = up(i,:)   ! velocity of the mirror is unchanged
#if dim3
             rp(k,3) = rp(i,3)
#endif 
             ! set the pressure difference between mirror and parent
             dP_mp(k) = dot_product(grav,rp(k,:)-rp(i,:))

          end do
       end if
    end do

    !! Make sure mirrors of solid particles have the same velocity as their parents...
    do j=n_par_fw+1,n_par_fw+imp
       i=irelation(j)
       tmp = dot_product(up(i,:),up(i,:))
       if(i<=n_par_w.and.tmp/=zero) then  !! if mirroring a moving solid particle
          up(j,:)=up(i,:)
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
             
             trans_bp(:) = rp(i,1:2)-b_node(ib,:)
             perp_dist = sqrt(dot_product(trans_bp,trans_bp))

             ! create a mirror particle
             imp = imp + 1
             k = n_par_fw + imp
             irelation(k) = i ! parent of the mirror
             vrelation(k) = -ib !mirror-parent velocity relationship

             ! position the mirror
             rp(k,1:2) = rp(i,1:2) - two*trans_bp(:)
#if dim3
             rp(k,3) = rp(i,3)
#endif 

#ifdef slip
             up(k,1:2) = -up(i,1:2)
#ifdef dim3
             up(k,3) = up(i,3)
#endif             
#else             
             up(k,:) = up(i,:)              ! don't reverse the velocity
#endif             

#ifdef wavemaker
             if(ib.eq.1) then !! If it's the first corner (wall-wavemaker)
                up(k,1) = two*U_wmaker - up(i,1)   !! Adjust velocities accordingly
             endif 
#endif

             ! set the pressure difference between mirror and parent
             dP_mp(k) = dot_product(grav,rp(k,:)-rp(i,:))
          end do
       end if

       ! OPTION 4: A double periodic corner 
       if(b_type(ib)==2.and.b_type(ibm1)==2)then
          ibp = b_periodic_parent(ib)  !identify parent patch for mirror particles
          do j=1,n_near_corner(ibp) ! look for particles near b_node(ibp)
             i=bc_neighbour(ibp,j)

             ! create a mirror particle
             imp = imp + 1
             k = n_par_fw + imp
             irelation(k) = i ! parent of the mirror
             vrelation(k) = -ib !mirror-parent velocity relationship

             ! position the mirror
             trans_bp(:) = b_node(ib,:) - b_node(ibp,:)
             rp(k,1:2) = rp(i,1:2) + trans_bp(:)
#if dim3
             rp(k,3) = rp(i,3)
#endif 
             ! don't change the velocity
             up(k,:) = up(i,:)

             ! set the pressure difference between mirror and parent
             dP_mp(k) = zero

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
             do j=n_par_fw+1,n_par_fw + imp_temp  ! look through existing mirror particles near wall
                perp_dist = bound_dist(j,ibb)
                if(perp_dist>sup_size+vsmall) cycle  ! cycle if it's not near the patch
                if(perp_dist<vsmall) cycle ! cycle if it's actually aligned with the patch
                norm_be(:) = bound_norm(ibb)  ! find the normal...
                tmp = position_along_patch(j,ibb)            
                ! is it off the end (the correct end) of the patch?
                if(tmp<=zero.and.ibb==ib.or.tmp>=one.and.ibb==ibm1) then
                   imp = imp + 1
                   k = n_par_fw + imp
                   ! parent of the mirror is the parent of the mirror in which it is mirrored!
                   i = irelation(j)
                   irelation(k) = i
                   vrelation(k) = -j !mirror-parent velocity relationship

                   ! position the mirror
                   rp(k,1:2) = rp(j,1:2) - two*perp_dist*norm_be(:)

#ifdef slip
                   up(k,1:2) = reverse_normal_velocity(norm_be(1),norm_be(2),up(j,1),up(j,2))
#else

                   tang_be(1)=norm_be(2);tang_be(2)=-norm_be(1)
                   up(k,1:2) = -up(j,1:2)  
#endif
#if dim3
                   rp(k,3) = rp(i,3);up(k,3)=-up(i,3)
#endif 
                   
                   ! set the pressure difference between mirror and parent
                   dP_mp(k) = dP_mp(j) ! diff between parent-mirror and parent-mirror-parent
                   dP_mp(k) = dP_mp(k) + dot_product(grav,rp(k,:)-rp(j,:)) !gravity
                end if
             end do
          else
             ! OPTION 5b: Neither of the two is a wall...(TBCompleted)  
          end if
       end if
    end do
    
    !! Make mirrors in third dimension - periodic only. 
    !! N.B. The third dimension is done last, and includes mirrors of mirror particles created for
    !! the first two dimensions.
#if dim3    
    nmirror = imp
    n_par_fwm = n_par_fw + nmirror
    do i=1,n_par_fwm

       if(inbin(i)) cycle ! Don't mirror particles in the bin
       if(rp(i,3)<=zb_min+sup_size)then !! It's near the minZ bound
          imp = imp + 1
          k = n_par_fw + imp
          if(i>n_par_fw)then
             irelation(k) = irelation(i)
             dp_mp(k) = dp_mp(i)
             if(vrelation(i)<-nb_patches) then  !! wall-periodic corner
                vrelation(k)=9999
             else if(vrelation(i)<0) then !! wall-wall or periodic-periodic corner
                vrelation(k)=9997
             else if(vrelation(i)<=nb_patches) then !! periodic, wall or inflow patch
                if(b_type(vrelation(i))==1)then ! wall
                   vrelation(k)=9996
                else
                   vrelation(k)=9995
                end if
             end if
          else
             irelation(k) = i
             dp_mp(k)=zero
             vrelation(k)=9998
          end if
#ifdef zwall
          up(k,:) = -up(i,:)
          rp(k,:) = rp(i,:);rp(k,3) = two*zb_min - rp(i,3)
#else          
          up(k,:) = up(i,:) 
          rp(k,:)=rp(i,:);rp(k,3)=rp(i,3) + zb_max-zb_min
#endif          
       end if
       if(rp(i,3)>=zb_max-sup_size)then !! It's near the maxZ bound
          imp = imp + 1
          k = n_par_fw + imp
          if(i>n_par_fw)then
             irelation(k) = irelation(i)
             dp_mp(k) = dp_mp(i)
             if(vrelation(i)<-nb_patches) then  !! wall-periodic corner
                vrelation(k)=9999
             else if(vrelation(i)<0) then !! wall-wall or periodic-periodic corner
                vrelation(k)=9997
             else if(vrelation(i)<=nb_patches) then !! periodic, wall or inflow patch
                if(b_type(vrelation(i))==1)then ! wall
                   vrelation(k)=9996
                else
                   vrelation(k)=9995
                end if
             end if
          else
             irelation(k) = i
             dp_mp(k)=zero
             vrelation(k)=9998
          end if
#ifdef zwall
          up(k,:) = -up(i,:)
          rp(k,:) = rp(i,:);rp(k,3) = two*zb_max - rp(i,3)
#else          
          up(k,:) = up(i,:) 
          rp(k,:)=rp(i,:);rp(k,3)=rp(i,3) - zb_max+zb_min
#endif          
       end if
    end do    
#endif    

   
    ! STEP 7: Update n_par_fwm and nmirror
    nmirror = imp
    n_par_fwm = n_par_fw + nmirror

    !! Copy the particle volumes
    !$omp parallel do private(i)
    do j=n_par_fw+1,n_par_fwm
       i=irelation(j)
       vol(j)=vol(i)
       h(j) = h(i)
    end do
    !$omp end parallel do

    
    deallocate(n_near_patch,n_near_corner)
    deallocate(b_neighbour,bc_neighbour)
    deallocate(b_neigh_pos)
    
    return
  end subroutine create_mirror_particles
!! ------------------------------------------------------------------------------------------------
  subroutine replace_escapees
    !! Check whether any particles have escaped. For particles which have passed through a wall,
    !! reflect their position and velocity. For particles which have passed through a periodic 
    !! boundary, insert them in the corresponding periodic parent patch.
    real(rkind),dimension(2) :: trans_bp
    real(rkind) tmp,m_be,u_norm
    integer(ikind) :: i,ib,ibp,n_escapees
    n_escapees = 0

    do ib = 1,nb_patches
       if(b_type(ib)==0) cycle  ! if it's a non-boundary, do nothing
       if(b_type(ib)==1.or.b_type(ib)==2) then   ! it's a wall or periodic
          ! loop over every particle
          norm_be(:) = bound_norm(ib)  !! N.B. wall normals point into fluid
          ibp = b_periodic_parent(ib)
          m_be = sqrt(dot_product(b_edge(ib,:),b_edge(ib,:)))
          
          
          !$omp parallel do private(perp_dist,tmp,trans_bp,u_norm) reduction(+:n_escapees)
          do i=n_par_w+1,n_par_fw
             if(.not.inbin(i)) then ! don't replace particles in bin
                ! find the perpendicular distance between i and the boundary edge
                perp_dist = bound_dist(i,ib) 

                ! for particles close but on the wrong side of the boundary
                if(perp_dist<-vsmall.and.abs(perp_dist)<=sup_size+vsmall) then  
                   tmp = position_along_patch(i,ib)
                   ! if i is close to the boundary patch but beyond the ends, cycle
                   ! it needs replacing, are we at a wall or a periodic?
                   if(b_type(ib)==1) then ! WALL
                      if(tmp>=zero.and.tmp<=one) then 
                         ! put particle back in (reflect in boundary)
                         rp(i,1:2) = rp(i,1:2) - two*perp_dist*norm_be(:)
                         u_norm=dot_product(up(i,1:2),norm_be(:))   !! normal component of velocity
                         if(u_norm<=zero)then
                            up(i,1:2) = up(i,1:2) + two*u_norm*norm_be(:) !! reflect normal component..
                         end if                         
                         write(6,*) "replaced particle which had passed through wall",i,ib,rp(i,:),perp_dist
                         n_escapees = n_escapees + 1
                      end if
                   else ! PERIODIC
                      if(tmp<zero.or.tmp>one) then  !! if unexpectedly about to leave via a corner
                         trans_bp(:) = -1.0*tmp*b_edge(ibp,:) - perp_dist*bound_norm(ibp)
                         rp(i,1:2) = b_node(ibp,:) + b_edge(ibp,:) + trans_bp(:) 
                         n_escapees = n_escapees + 1               
                      else
                         ! reposition the particle relative to the boundary patches, which
                         ! allows for non-parallel patch-parents (eg. bent Poiseuille flow)
                         trans_bp(:) = -1.0*tmp*b_edge(ibp,:) - perp_dist*bound_norm(ibp)
                         rp(i,1:2) = b_node(ibp,:) + b_edge(ibp,:) + trans_bp(:)
                         n_escapees = n_escapees + 1
                      end if
                   end if
                   ! set the o-arrays, as we're doing inflow/outflow after o-array allocation
                   rpo(i,:) = rp(i,:) - up(i,:)*dt
                   upo(i,:) = up(i,:)
                end if
             end if
          end do
          !$omp end parallel do
       end if
    end do
    
#if dim3    
    !! 3rd dimensions
    !$omp parallel do
    do i=1,n_par_fw
       if(.not.inbin(i)) then
          if(rp(i,3)<zb_min)then
#ifdef zwall
             rp(i,3) = two*zb_min - rp(i,3);up(i,3) = -up(i,3)   !! If walls, reflect and then reverse wallnormal u
#else          
             rp(i,3) = rp(i,3) + zb_max - zb_min  !! If periodic, just translate 
#endif             
             rpo(i,:) = rp(i,:) - up(i,:)*dt;upo(i,:) = up(i,:)
          end if
          if(rp(i,3)>=zb_max)then
#ifdef zwall
             rp(i,3) = two*zb_max - rp(i,3);up(i,3) = -up(i,3)
#else
             rp(i,3) = rp(i,3) - zb_max + zb_min
#endif             
             rpo(i,:) = rp(i,:) - up(i,:)*dt;upo(i,:) = up(i,:)
          end if
       end if
    end do
    !$omp end parallel do
#endif    
    
  end subroutine replace_escapees
!! ------------------------------------------------------------------------------------------------
  subroutine insert_or_remove_particles
    real(rkind) alph_inflow,m_be
    integer(ikind) :: i,j,ib,k,n_left,n_joined

    ! If we want a non-constant U_inflow, set it here, for example:
    !U_inflow = 0.75 + 0.375*sin(real(itime/200))    or
    !U_inflow = time
    n_left=0;n_joined=0
    do ib = 1,nb_patches
       if(b_type(ib)<=2) cycle ! Do nothing for non, wall or periodic patches
       if(b_type(ib)==4) then !if outflow, we will deal with it first        
          norm_be(:) = bound_norm(ib)             
          do i=n_par_w+1,n_par_fw
             perp_dist = bound_dist(i,ib) ! how far is i from ib?
             if(perp_dist>0.or.perp_dist<=-sup_size) cycle  !! if it's in the domain or too far out, neglect
             inbin(i) = .true.
             n_inbin = n_inbin + 1
             rp(i,:) = -1.0d10*dx
             n_free = n_free + 1   ! add 1 to the number of free indices for particles
             p_free(n_free) = i    ! store the free index in the list of free indices
             n_left=n_left+1
          end do
       else if(b_type(ib)==3) then ! if it's an inflow boundary, deal with it second
          alph_inflow = 1.1d0*dx
          norm_be(:) = bound_norm(ib)
          m_be = sqrt(dot_product(b_edge(ib,:),b_edge(ib,:)))
          ! loop through the inflow particles
          do i=1,np_inflow
             j = p_inflow(i) ! j is the i'th inflow particle
             perp_dist = bound_dist(j,ib) ! how far is j from ib?
             if(perp_dist>alph_inflow) then
                n_joined=n_joined+1
                ! j has moved away from ib, so create a new particle
                if(n_free==0) then   ! if no free indices, increment n_par_fw
                   n_par_fw = n_par_fw + 1
                   k=n_par_fw
                else 
                   k = p_free(n_free)  ! use last free index
                   n_free = n_free - 1  ! subtract 1 from the number of free indices for particles
                   inbin(k) = .false.    ! take it out of the bin...
                   n_inbin = n_inbin - 1
                end if
                p_inflow(i) = k   ! mark the new particle as inflow
                ! position the new particle the correct distance along ib
!                tmp = (real(i)-half)*dx/m_be  !! This for jet

!                rp(k,1:2) = b_node(ib,:) + tmp*b_edge(ib,:) + (perp_dist-dx)*norm_be(:)
!#if dim3
!                rp(k,3) = rp(j,3)
!#endif                
                rp(k,:) = rp_inflow(i,:)
                
                ! give the new particle a velocity
                up(k,1:2) = norm_be(:)*U_inflow

                ! set the o-arrays, as we're doing inflow/outflow after o-array allocation
                rpo(k,:) = rp(k,:) - up(k,:)*dt
                upo(k,:) = up(k,:)
             end if
          end do
       end if
    end do
    n_diff = n_diff + n_joined - n_left
!    write(1,*) itime,n_diff
    return
  end subroutine insert_or_remove_particles
!! ------------------------------------------------------------------------------------------------
  subroutine create_boundary_boxes
    !! Loops through all boundary patches and corners, and builds lists of particles close to each
    !! patch and corner, and their position relative to that patch.
    integer(ikind) :: i,ib
    real(rkind) :: tmp

    n_near_patch(:) = 0
    n_near_corner(:) = 0
    do ib = 1,nb_patches
       if(b_type(ib)==0) cycle ! if it's a non-boundary, i don't care
       norm_be(:) = bound_norm(ib)
       do i=1,n_par_fw
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
       do i=n_par_w+1,n_par_fw
          perp_dist = sqrt((rp(i,1) - b_node(ib,1))**2 + (rp(i,2)-b_node(ib,2))**2)
          if(perp_dist>sup_size+vsmall) cycle
          ! it is near a corner!!
          n_near_corner(ib) = n_near_corner(ib) + 1
          bc_neighbour(ib,n_near_corner(ib)) = i
       end do
    end do

    return
  end subroutine create_boundary_boxes
!! ------------------------------------------------------------------------------------------------
  function bound_dist(i,ib) result(p_dist)
    integer(ikind),intent(in) :: i,ib
    real(rkind) :: p_dist
    integer(ikind) :: ibp1

    ! the ib+1'th patch is 1 if ib = nb_patches
    ibp1 = mod(ib,nb_patches)+1

    ! find the perpendicular distance (signed) between particle i
    ! and boundary patch ib.
    p_dist = b_edge(ib,2)*rp(i,1) - b_edge(ib,1)*rp(i,2) + &
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

    nrm(1) = -1.0*b_edge(ib,2)/tmp1           !! unit normal
    nrm(2) = b_edge(ib,1)/tmp1
  end function bound_norm
!! ------------------------------------------------------------------------------------------------
  function position_along_patch(i,ib) result(tmp)
    !! Returns the distance along the boundary patch ib of particle i:
    !! rp(i,:) = b_node(ib,:) + tmp*b_edge(ib,:) + perp_dist*norm_be(:)
    !!
    !! N.B. This function assumes that perp_dist and norm_be have already been obtained for
    !! patch ib.
    integer(ikind),intent(in) :: i,ib
    real(rkind) :: tmp
    
    ! find position along the boundary patch which is closest to particle i
    if(abs(b_edge(ib,1))>abs(b_edge(ib,2)))then
       tmp = (rp(i,1) - b_node(ib,1) - perp_dist*norm_be(1))/b_edge(ib,1)
    else
       tmp = (rp(i,2) - b_node(ib,2) - perp_dist*norm_be(2))/b_edge(ib,2)
    end if
             
  end function position_along_patch
!! ------------------------------------------------------------------------------------------------
  function angle_between_patches(a,b) result(Rmat)
     !! Returns a rotation matrix to rotate a vector aligned with patch a into alignment with patch b (I think...)
     integer(ikind),intent(in) :: a,b
     real(rkind), dimension(2,2) :: Rmat
     real(rkind) tmp_angle

     tmp_angle = dot_product(b_edge(a,:),b_edge(b,:))/ &
        sqrt(dot_product(b_edge(a,:),b_edge(a,:))*  &
        dot_product(b_edge(b,:),b_edge(b,:)))
     tmp_angle = acos(tmp_angle) - pi
     Rmat(1,1) =  cos(tmp_angle)
     Rmat(1,2) = -sin(tmp_angle)
     Rmat(2,1) = sin(tmp_angle)
     Rmat(2,2) = cos(tmp_angle)
   end function angle_between_patches
!! ------------------------------------------------------------------------------------------------
   subroutine mirror_velocities
   ! set mirror velocities (ie. after viscous, before ppe setup)
     integer(ikind) i,j,ib,ibm1
     real(rkind),dimension(dims) :: rij
     real(rkind),dimension(2) :: tang_be
     real(rkind) :: tmp
     
     !! j is the mirror particle under consideration
     !! i is the parent of j
     
     do j=n_par_fw + 1,n_par_fwm
        i=irelation(j)
        
        
        !! NON-CORNER
        if(vrelation(j)>0.and.vrelation(j)<999) then 
        
           !! WALL
           if(abs(b_type(vrelation(j)))==1)then 
              rij = half*(rp(i,:)-rp(j,:));tmp=sqrt(dot_product(rij,rij))
              tang_be(1)=rij(2)/tmp;tang_be(2)=-rij(1)/tmp

#ifdef slip 
              up(j,1:2) = reverse_normal_velocity(rij(1)/tmp,rij(2)/tmp,up(i,1),up(i,2))
#else
              up(j,1:2) = - up(i,1:2)
#endif             
#if dim3
              up(j,3)=-up(i,3)
#endif 

#ifdef wavemaker
              if(vrelation(j)==nb_patches) then !! If it's the final patch
                 up(j,1) = up(j,1) + two*U_wmaker
              end if
#endif

           else if(b_type(vrelation(j))==2) then ! periodic
              up(j,:) = up(i,:)  ! no change in velocity (assuming parallel periodic patches)
           else if(b_type(vrelation(j))==3) then ! inflow, u* is U_inflow
              up(j,:) = up(i,:)        
           else if(b_type(vrelation(j))==4)then ! outflow
              up(j,:) = up(i,:) !! No change in velocity (du/dn = 0)
           end if
        end if
        
        !! CORNER (wall-wall or periodic-periodic)
        if(vrelation(j)<0.and.vrelation(j)>-n_par_fw)then
        
           ib = -vrelation(j)
           ibm1 = mod(ib+nb_patches-2,nb_patches)+1      

           !! WALL WALL CORNER
           if(b_type(ib)==1.and.b_type(ibm1)==1) then 

#ifdef slip
                 up(j,1:2) = -up(j,1:2)
#ifdef dim3
                 up(j,3) = up(j,3)
#endif                 
#else
                 up(j,:) = up(i,:)  ! don't reverse velocity
#endif                 

#ifdef wavemaker
                 if(ib==1) then !! If the first corner, wavemaker-wall
                    up(j,1) = two*U_wmaker - up(i,1)
                 end if
#endif

           !! PERIODIC PERIODIC CORNER
           else if(b_type(ib)==2.and.b_type(ibm1)==2) then 
              up(j,:) = up(i,:)   ! assuming parallel patches
           end if
           
        !! WALL-OTHER CORNER   
        else if(vrelation(j)<-n_par_fw) then ! wall-periodic or wall-inflow corner
#ifdef slip
           rij = half*(rp(-vrelation(j),:)-rp(j,:));tmp=sqrt(dot_product(rij,rij))
           up(j,1:2) = reverse_normal_velocity(rij(1)/tmp,rij(2)/tmp,up(-vrelation(j),1),up(-vrelation(j),2))
#else
           up(j,:) = -up(-vrelation(j),:) 
#endif           
        end if
        
#ifdef dim3        
        !! Third dimension mirrors
        if(vrelation(j)==9999)then ! parent is wall-periodic corner  
#ifdef zwall
           up(j,:) = up(i,:)
#else        
           up(j,:) = -up(i,:)
#endif           
        end if
        if(vrelation(j)==9998)then ! parent is fluid
#ifdef zwall
           up(j,:) = -up(i,:)
#else
           up(j,:) = up(i,:)
#endif           
        end if
        if(vrelation(j)==9997)then ! parent is wall-wall or periodic-periodic corner
#ifdef zwall
           up(j,:) = -up(i,:)
#else        
           up(j,:) = up(i,:)
#endif           
        end if
        if(vrelation(j)==9996)then ! parent is wall    
#ifdef zwall
           up(j,:) = up(i,:)
#else
           up(j,:) = -up(i,:)
#endif           
        end if
        if(vrelation(j)==9995)then ! parent is periodic,inflow  
#ifdef zwall
           up(j,:) = -up(i,:)
#else
           up(j,:) = up(i,:)
#endif           
        end if        
#endif
     end do

     !! Make sure mirrors of solid particles have the same velocity as their parents...
     do j=n_par_fw+1,n_par_fwm
        i=irelation(j)
        tmp = dot_product(up(i,:),up(i,:))
        if(i<=n_par_w.and.tmp/=zero) then  !! if mirroring a moving solid particle
           up(j,:)=up(i,:)
        end if
     end do
   return
   end subroutine mirror_velocities
!! ------------------------------------------------------------------------------------------------
  function reverse_normal_velocity(n_n1,n_n2,u_or1,u_or2) result(u_out)
    ! returns the unit normal vector to the boundary edge
    real(rkind),intent(in) :: n_n1,n_n2,u_or1,u_or2
    real(rkind),dimension(2) :: u_out
    real(rkind),dimension(2) :: xv,yv,n_t,u_tmp,n_n,u_or
    real(rkind),dimension(2,2) :: Qij,QijT
    
    !! Pass scalars to vectors
    n_n(1) = n_n1;n_n(2) = n_n2
    u_or(1) = u_or1;u_or(2) = u_or2
    
    !! Initialise unit vectors
    xv = (/1,0/);yv = (/0,1/)
    
    !! Build tangent
    n_t(1) = -n_n(2);n_t(2) = n_n(1)

    !! Build rotation matrix
    Qij(1,1) = dot_product(xv,n_n);Qij(1,2) = dot_product(xv,n_t)
    Qij(2,1) = dot_product(yv,n_n);Qij(2,2) = dot_product(yv,n_t)

    !! Build reverse matrix
    QijT = transpose(Qij)
    
    !! Rotate into FoR with wall
    u_tmp = matmul(Qij,u_or)
    
    !! Reverse normal component
    u_tmp(1) = -u_tmp(1)
    
    !! Rotate back into x-y FoR
    u_out = matmul(QijT,u_tmp)
 
  end function reverse_normal_velocity
!! ------------------------------------------------------------------------------------------------   
end module mirror_boundaries
