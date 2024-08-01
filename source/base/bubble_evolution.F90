module bubble_evolution
  !! This module contains routines which evaluate the bubble density, update the bubble velocity
  !! an model sub-resolution fluctuations, free-surface interaction, and breakup.
  !! N.B. bubble entrainment is handled in bubble_boundaries.
  use kind_parameters
  use common_parameter
  use common_2d
  use calculation_gradient
  use bubble_boundaries
  use omp_lib
  implicit none
  
  private
#ifdef bubbles  
  public ::  update_bubble_density,put_bubbles_in_bin,calculate_bubble_velocities
#endif  

  real(rkind),dimension(:,:),allocatable :: vel_l2b,snorm_l2b
  real(rkind),dimension(:,:,:),allocatable :: gradvel_l2b
  real(rkind),dimension(:),allocatable :: vol_l2b,tke_l2b,eps_srs_l2b
  real(rkind),dimension(:,:),allocatable :: M0_bubbles  

contains
!! ------------------------------------------------------------------------------------------------
#ifdef bubbles
  subroutine update_bubble_density
     !! This routine evaluates the liquid volume fraction by interpolating the bubble volumes
     !! onto the SPH particle positions
     use mirror_boundaries
     use sphtools
     integer(ikind) :: i,j,k
     real(rkind),dimension(dims) :: rij,gradw
     real(rkind) :: qq,tmp_vf,bubble_volume,Vgas,kcoeff,wqq,rad,tmpa,tmpv
     real(rkind),dimension(:),allocatable :: tmp_smoothed
     
     !! Loop over all SPH particles         
     !$omp parallel do private(tmp_vf,k,j,rij,qq,wqq,bubble_volume,Vgas,kcoeff,gradw,rad)
     do i=1,n_par_fw

        !! Loop over bubble neighbours to find bubble volumes
        tmp_vf = zero
        kcoeff = dv/vol(i)   !! previous Dilation coefficient
        do k=1,ijb_count(i)
           j = ijb_linearlink(ijb_firstindex(i)+k-1)           
           rij = rp(i,:) - rb(j,:)
          
           rad = sqrt(dot_product(rij,rij));qq=rad/h(i)
#if dim3
           bubble_volume = (four/three)*pi*radb(j)**three
#else           
           bubble_volume = radb(j)**two
#endif
           wqq = wab(qq)
           
           !! Sum of bubble volumes interpolated to SPH particles
           tmp_vf = tmp_vf + wqq*bubble_volume
        end do
        
        !! Vgas is the volume of bubbles associated with this SPH particle
        Vgas = kcoeff*dv*tmp_vf
        
        !! ARTIFICIAL LIMIT: Set Vgas <= half*dv, limiting a0<=1/3. (clearly violates volume conservation)
        if(Vgas>half*dv) Vgas = half*dv
        
        !! The total volume of this SPH particle is now liquid + gas volume
        vol(i) = dv + Vgas
        
        !! The liquid volume fraction is liquid volume divided by total volume
        a0(i) = dv/vol(i)
     end do
     !$omp end parallel do
     
     !! Adjust the smoothing length to account for the volume expansion
     !$omp parallel do
     do i=1,n_par_fw
#if dim3
        h(i) = h0*(a0(i))**(-one/three)
#else
        h(i) = h0*(a0(i))**(-one/two)
#endif        
      end do
     !$omp end parallel do
            
     
     !! Mirror the volume fractions, volumes, and smoothing lengths.
     !$OMP PARALLEL DO PRIVATE(j)
     do i=n_par_fw+1,n_par_fwm
        j=irelation(i)
        a0(i) = a0(j)
        vol(i) = vol(j)
        h(i) = h(j)
     end do
     !$OMP END PARALLEL DO     
     
     !! Number of bubbles and volume of bubbles entrained
     tmpa=zero
     j=0
     !$omp parallel do reduction(+:j,tmpa)
     do i=1,n_bub
        if(.not.b_inbin(i)) then
           j=j+1
           tmpa=tmpa + radb(i)**three                
        end if
     end do
     !$omp end parallel do
     tmpa = tmpa*four*pi/three     
     write(186,*) time,tmpa,dble(j)
     flush(186)
     
                      
     
     return
  end subroutine update_bubble_density
!! ------------------------------------------------------------------------------------------------
  subroutine calculate_bubble_velocities
     !! Interpolate properties from liquid to bubble, calculate fluctuations, evolve bubble 
     !! velocities one step (many substeps depending on bubble size), calculate break up
     !! and surface interactions etc...
     use sphtools
     integer(ikind) :: i,j,k
     real(rkind),dimension(dims) :: rij,urel,fdrag,fvm,tmp_M0,ub0,fb,dudt_b,fl,tmp_vec,fsurf,tmpm
     real(rkind) :: qq,wqq,bubble_volume,tmp_wsum,coeff_vm,urelmag,rad,tmpv,tmpa
     real(rkind) :: Cd,Reb,wx,wy,wz,Cl,phimax,dt_bubble,tmp_nu_BI,temp_L,temp,tmp,tmp2
     real(rkind) :: kcoeff,b_nearsurf,snorm_mag,We_b,ubsn,fsmag,Tc
     integer(ikind),dimension(:),allocatable :: n_neighbours_l2b,n_fsneighbours_l2b
     real(rkind),dimension(:),allocatable :: tmp_smoothed
           
     !! Interpolate various properties from liquid to bubbles
     !! N.B. all properties interpolated from liquid to bubbles have suffix l2b     
     allocate(vel_l2b(n_bub_m,dims));vel_l2b=zero
     allocate(gradvel_l2b(n_bub_m,dims,dims));gradvel_l2b = zero
     allocate(vol_l2b(n_bub_m),tke_l2b(n_bub_m));vol_l2b=zero;tke_l2b=zero
     allocate(snorm_l2b(n_bub_m,dims));snorm_l2b=zero
     allocate(n_neighbours_l2b(n_bub_m));n_neighbours_l2b=0
     allocate(n_fsneighbours_l2b(n_bub_m));n_fsneighbours_l2b=0
     allocate(eps_srs_l2b(n_bub_m));eps_srs_l2b =zero
     !$omp parallel do private(k,j,rij,qq,wqq) reduction(+:vel_l2b,gradvel_l2b,vol_l2b,tke_l2b &
     !$omp ,n_neighbours_l2b,n_fsneighbours_l2b,snorm_l2b,eps_srs_l2b)
     do i=1,n_par_fwm   
        do k=1,ijb_count(i)
           j = ijb_linearlink(ijb_firstindex(i)+k-1)           

           rij = rp(i,:)-rb(j,:);qq = sqrt(dot_product(rij,rij))/h(i)
           wqq = wab(qq)*vol(i)
           
           !! SPH particle concentration (for normalisation)
           vol_l2b(j) = vol_l2b(j) + wqq

           !! Liquid velocity and velocity gradient           
           vel_l2b(j,:) = vel_l2b(j,:) + wqq*up(i,:)
           gradvel_l2b(j,:,:) = gradvel_l2b(j,:,:) + wqq*grad_vel(:,:,i)   

           !! sub-grid Turbulent kinetic energy
           tke_l2b(j) = tke_l2b(j) + wqq*tke(i) 
           
           !! sub-grid turbulent dissipation rate
           eps_srs_l2b(j) = eps_srs_l2b(j) + wqq*eps_srs(i)          
           
           !! Total and free-surface neighbours     
           n_neighbours_l2b(j) = n_neighbours_l2b(j) + 1
           n_fsneighbours_l2b(j) = n_fsneighbours_l2b(j) + n_surf(i)    
           
           !! Surface normals
           snorm_l2b(j,:) = snorm_l2b(j,:) + wqq*surf_norm(i,:)                   
        end do
     end do
     !$omp end parallel do

     !! De-allocate memory no longer required
     if(allocated(tke)) deallocate(tke)
     allocate(bnfs(n_bub_m));bnfs = zero     !! bnfs is a measure of bubble-near-free-surface        
            
     !! Put bubbles in bin if [short of SPH neighbours && all neighbours are free surface]
     !! Also evaluate bnfs
     !$omp parallel do
     do i=1,n_bub
        !! Bubble-near-free-surface (0.353d0 is magnitude of surface normal at plane surface)
        bnfs(i) = sqrt(dot_product(snorm_l2b(i,:),snorm_l2b(i,:)))/0.353d0/vol_l2b(i)
        if(bnfs(i)<=1.0d-4.or.i_open_domain==0) bnfs(i)=vsmall
        if(n_neighbours_l2b(i)<=20.and.n_fsneighbours_l2b(i)==n_neighbours_l2b(i)) then
           if(.not.b_inbin(i))then
              b_inbin(i) = .true.           !! Put in the bin
              rb(i,:) = -1.0d4   !! Put it out of the way
              ub(i,:) = zero 
           end if          
        end if

     end do
     !$omp end parallel do  

     !! Normalise with volume all properties interpolated from liquid to bubbles
     !$omp parallel do private(snorm_mag)
     do i=1,n_bub_m
        vel_l2b(i,:) = vel_l2b(i,:)/vol_l2b(i)
        gradvel_l2b(i,:,:) = gradvel_l2b(i,:,:)/vol_l2b(i)      
        tke_l2b(i) = tke_l2b(i)/vol_l2b(i)      
        eps_srs_l2b(i) = eps_srs_l2b(i)/vol_l2b(i)
        
        !! Except surface normals, which become unit vectors
        snorm_mag = max(sqrt(dot_product(snorm_l2b(i,:),snorm_l2b(i,:))),vsmall)
        snorm_l2b(i,:) = snorm_l2b(i,:)/snorm_mag           
     end do
     !$omp end parallel do
          
    
     !! Lift coefficient
     Cl = half
         
     !! Evaluate the sub-grid fluctuations "felt" or "seen" by the bubbles
     call calc_subgrid_fluctuations
     
     !! Output some tracer positions if required. Used for some bubble plume stuff
     if(.false.) call bubble_traces          
    
     !! Loop over all bubbles, evaluate forces and update velocities
     allocate(M0_bubbles(n_bub_m,dims));M0_bubbles=zero
     !$omp parallel do private(urel,fdrag,bubble_volume,fvm,ub0,coeff_vm,fb,urelmag &
     !$omp ,Cd,Reb,wx,wy,wz,fl,tmp_vec,b_nearsurf,We_b,fsurf,ubsn,fsmag &
     !$omp ,Tc,dt_bubble,j,k)
     do i=1,n_bub     
        if(.not.b_inbin(i)) then                            
          
           !! Increase the age of all bubbles
           bubble_age(i) = bubble_age(i) + dt         
        
           !! Keep track of the mean dissipation rate (of turbulent energy in the liquid) experienced
           !! by the bubble over its lifetime...
           mean_dissrate(i) = mean_dissrate(i)*(bubble_age(i)-dt)/bubble_age(i) + eps_srs_l2b(i)*dt/bubble_age(i)        
        
           !! Initial bubble velocity
           ub0 = ub(i,:)                      
           
           !! Time-derivative of liquid velocity in FoR moving with bubble
           !! N.B. u_l2b holds the liquid velocity at the last (SPH) time-step. vel_l2b holds same at this step.
           dudt_b(:) = (vel_l2b(i,:) - u_l2b(i,:))/dt           
                                
           !! Modify "felt" velocity with the subgrid fluctuations           
           vel_l2b(i,:) = vel_l2b(i,:) + ubturb(i,:)

b_out(i) = sqrt(dot_product(ubturb(i,:),ubturb(i,:)))     


#if dim3
           !! Bubble volume                     
           bubble_volume = (four/three)*pi*radb(i)**three
           coeff_vm = half          

           !! Components of vorticity
           wx = gradvel_l2b(i,3,2)-gradvel_l2b(i,2,3)
           wy = gradvel_l2b(i,1,3)-gradvel_l2b(i,3,1)
           wz = gradvel_l2b(i,2,1)-gradvel_l2b(i,1,2)        
#else
           bubble_volume = radb(i)**two
           coeff_vm = four
        
           !! Single non-zero component of vorticity
           wz = gradvel_l2b(i,2,1)-gradvel_l2b(i,1,2)        
#endif     

           !! Set bubble integration time as fraction of SPH time
           !! Sort of empirical at present. In due course need to calculate stability criteria based on
           !! stiffness of system and set it that way.
           j = ceiling(dx/radb(i))
           dt_bubble = dt/dble(j)


           !! Loop for j small steps
           do k=1,j                
                      
              !! Relative (felt) velocity
              urel = vel_l2b(i,:) - ub(i,:)            
              urelmag = sqrt(dot_product(urel,urel))
              
              
              !! Evaluate bubble Reynolds number and drag coefficient
              Reb = two*radb(i)*urelmag*Re
              if(Reb>1000.0) then
                 Cd = 0.44
              else
                 Cd = (24.0d0/max(Reb,0.0001))*(one+0.15d0*Reb**0.687d0)
              end if
              
              !! Evaluate deformation velocity
              ubsn = 8.2d0*(two*eps_srs_l2b(i)*radb(i))**twothirds - 6.0/We/radb(i)
              if(ubsn>zero)then
                 ubsn = sqrt(ubsn)
              else
                 ubsn = -sqrt(-ubsn)
              end if

              !! Integrate deformation velocity to obtain deformation distance 
              if(ubsn>zero) then
                 deformation_distance(i) = deformation_distance(i) + ubsn*dt_bubble  
              end if                                   
              !! The deformation can't be less than zero (spherical)
              deformation_distance(i) = max(deformation_distance(i),zero)

             

              !! Drag force, volumes and VM coefficient, vorticity and lift force
#if dim3
              fdrag = (half*pi*radb(i)**two)*beta*urel*urelmag*Cd       
 
              !! Slip velocity cross vorticity
              tmp_vec(1) = urel(2)*wz - urel(3)*wy
              tmp_vec(2) = urel(3)*wx - urel(1)*wz
              tmp_vec(3) = urel(1)*wy - urel(2)*wx        
#else
              fdrag = half*radb(i)*beta*urel*urelmag*Cd 
     
              !! Slip velocity cross vorticity
              tmp_vec(1) = urel(2)*wz
              tmp_vec(2) = urel(1)*wz
#endif     

              !! Lift force
              fl = beta*Cl*bubble_volume*tmp_vec(:)

              !! Buoyancy force
              fb = (one-beta)*bubble_volume*(grav)
          
              !! Partial virtual mass force
              fvm = beta*(coeff_vm)*bubble_volume*dudt_b(:)

           
              !! THIS SECTION ON FREE SURFACE INTERACTION 
              if(bubble_EoL(i)==0.and.bnfs(i)>0.01d0) then  !! For bubbles not yet interacting with FS

                 !! Relative normal velocity and bubble Weber number
                 ubsn = dot_product(urel,snorm_l2b(i,:))
                 We_b = We*radb(i)*urelmag**two

                 !! Non-zero timescale to prevent bubbles being marked for bursting at instant of creation
                 !! Also gives a measure of the timescale for decelleration
                 Tc = dx/urelmag

                 !! Mark bubbles as in fscontact if required
                 if(bnfs(i)>0.1d0.and.bubble_age(i)>Tc) then 
                    bubble_EoL(i)=1               
                    bubble_LE(i) = bubble_age(i) + Tc !! Expect it to remain for another Tc time units                   
                 end if
              end if
              
              !! For bubbles in contact
              if(bubble_EoL(i)==1) then

                 !! Evaluate the relative velocity component normal to the surface              
                 ubsn = dot_product(urel(:),snorm_l2b(i,:))
              
                 !! The free-surface force 
                 fsmag = ubsn*(one+coeff_vm*beta)*bubble_volume/max(bubble_LE(i)-bubble_age(i),dt_bubble) - & 
                                  dot_product(fb+fdrag+fvm+fl,snorm_l2b(i,:))
                 fsurf = fsmag*snorm_l2b(i,:)
                 ubsn = half*(one+erf(two*log(five*bnfs(i))))   !! Smoothing function
                 fsurf = fsurf*ubsn
                 
                 !! If no longer in contact
                 if(bnfs(i)<=0.1) then
                    bubble_EoL(i) = 0      !! A new lease of life, if it has moved away from FS...
                 end if
              else
                 fsurf = zero
              end if                                 
              !! END OF FREE SURFACE INTERACTION !! -----------------------------------------------
           
              !! Update bubble velocity
              ub(i,:) = ub(i,:) + dt_bubble*(fb+fdrag+fvm+fl+fsurf)/((one+coeff_vm*beta)*bubble_volume)      
        
           end do

           !! Calculate final virtual mass force
           fvm = fvm - beta*coeff_vm*bubble_volume*(ub(i,:)-ub0(:))/dt

           !! Total momentum exchange between bubble and liquid
           M0_bubbles(i,:) = -(fdrag + fvm + fl)/beta !! Divide by beta because liquid momentum eqn is u, not rho*u
              
           !! Store fluid velocity at current bubble position for next time step
           u_l2b(i,:) = vel_l2b(i,:)

        end if   !! End of inbin if
     end do     
     !$omp end parallel do
     
     !! Check whether any bubble has NaN velocity. If so, put it in the bin so it doesn't
     !! ruin the SPH simulation
     call put_bubbles_in_bin_velocity_only
     
     !! Mirror momentum exchanges
     !$omp parallel do private(i)
     do j=n_bub+1,n_bub_m
        i=irelation_b(j);M0_bubbles(j,:)=M0_bubbles(i,:)
     end do
     !$omp end parallel do

     !! Interpolate the momentum exchange back onto the liquid particles
     !$omp parallel do private(tmp_M0,k,j,rij,qq,kcoeff)
     do i=1,n_par_fw
        tmp_M0 = zero
        kcoeff = dv/vol(i)
        do k=1,ijb_count(i)
           j = ijb_linearlink(ijb_firstindex(i)+k-1)           

           rij = rp(i,:) - rb(j,:)
           qq = sqrt(dot_product(rij,rij))/h(i)
           tmp_M0(:) = tmp_M0(:) + wab(qq)*M0_bubbles(j,:)
        end do
        M0(i,:) = kcoeff*tmp_M0(:)*vol(i)
     end do
     !$omp end parallel do

     !! Re-mirror bubble velocities.
     call mirror_velocities_b         

     !! Interpolate the bubble induced turbulence back onto the liquid particles
     !$omp parallel do private(tmp_nu_BI,k,j,rij,qq,kcoeff,bubble_volume)
     do i=1,n_par_fw
        tmp_nu_BI = zero
        kcoeff = dv/vol(i)
        do k=1,ijb_count(i)
           j = ijb_linearlink(ijb_firstindex(i)+k-1)           
           rij = rp(i,:) - rb(j,:);qq = sqrt(dot_product(rij,rij))/h(i)

#if dim3
           bubble_volume = (four/three)*pi*radb(j)**three
#else           
           bubble_volume = radb(j)**two
#endif
 
           !! Roughly equiv. to Sato & Sekoguchi 1975.
           tmp_nu_BI = tmp_nu_BI + kcoeff*wab(qq) &
                                   *radb(j)*two &
                                   *sqrt(dot_product(vel_l2b(j,:)-ub(j,:),vel_l2b(j,:)-ub(j,:))) &
                                   *vol(i) ! VOLUME TERM                                 
        end do
        nu_BI(i) =  0.6*(one-a0(i))*tmp_nu_BI*Re
     end do
     !$omp end parallel do   
            
     !! Determine whether any bubbles break up            
     call bubble_breakup
          
     !! Deallocate memory     
     deallocate(vel_l2b,M0_bubbles,gradvel_l2b,vol_l2b,tke_l2b,snorm_l2b)
     deallocate(n_neighbours_l2b,n_fsneighbours_l2b) 
     deallocate(eps_srs_l2b)      
    

     return
  end subroutine calculate_bubble_velocities
!! ------------------------------------------------------------------------------------------------  
  subroutine put_bubbles_in_bin
     !! Determine if a bubble is out of the domain, and if so, put it in the bin.
     !! This routine will only be called prior to neighbour finding (as it resets
     !! b_nfree and does a full recount of free indices.
     integer(ikind) :: i,n_death,j
     logical :: bad_bubble
     real(rkind) :: tmp,tmp2,bubble_extra_time
     
     !! Initialise counters, then loop over all bubbles
     b_nfree = 0
     n_death = 0
     do i=1,n_bub
        !! Initial assumption that all bubbles are OK.
        bad_bubble = .false.
        
        !! Test if outside x-y domain 
        if(rb(i,1)<xb_min-sup_size.or.rb(i,1)>xb_max+sup_size.or. &
           rb(i,2)<yb_min-sup_size.or.rb(i,2)>yb_max+sup_size)then      
           bad_bubble = .true.
        end if
        
        !! Test if outside z domain
#ifdef dim3
!#ifdef zwall
        if(rb(i,3)<zb_min-sup_size.or.rb(i,3)>zb_max+sup_size) then !Outside z if walls in z 
           bad_bubble = .true.
        end if
!#endif     
#endif   
        !! Test if bubble has exceeded life expectancy + extra time.
        !! N.B. "bubble_extra_time" is the persistence time predicted by the model from 
        !! S. Poulain et al. (2018) JFM 851:636–671. doi: 10.1017/jfm.2018.471.        
        bubble_extra_time = sqrt(radb(i))*(We**(3.0/4.0))*sqrt(Fr)*(one/Re)*500.0
        if(bubble_EoL(i)==1.and.bubble_age(i)>bubble_LE(i)+bubble_extra_time)then 
           n_death = n_death + 1
           bad_bubble = .true.               
        end if
        
        !! Test if position or velocity is NaN
        do j=1,dims
           if(isnan(rb(i,j))) bad_bubble = .true.
           if(isnan(ub(i,j))) bad_bubble = .true.
        end do

        !! Put in bin if any tests failed
        if(bad_bubble) then
           call bubble_in_bin(i)
           b_nfree = b_nfree + 1           
        end if
               
     end do
        
     !! Output number of burst bubbles/time step   
     write(189,*) time,dble(n_death)/max(dt,vsmall)
     flush(189)        
     
     return
  end subroutine put_bubbles_in_bin
!! ------------------------------------------------------------------------------------------------
  subroutine put_bubbles_in_bin_velocity_only
     !! Determine if a bubble is out of the domain, and if so, put it in the bin.
     !! Only testing the velocity
     integer(ikind) :: i,n_death,j
     logical :: bad_bubble
     
     b_nfree = 0
     n_death = 0
     do i=1,n_bub
        !! Initial assumption that all bubbles are OK.
        bad_bubble = .false.
        
        !! Test if position or velocity is NaN
        do j=1,dims
           if(isnan(ub(i,j))) bad_bubble = .true.
           if(isnan(vel_l2b(i,j))) bad_bubble = .true.
           if(isnan(ubturb(i,j))) bad_bubble = .true.
        end do

        !! Put in bin if any tests failed
        if(bad_bubble) then
           call bubble_in_bin(i)
           b_nfree = b_nfree + 1    
           M0_bubbles(i,:) = zero  
           vel_l2b(i,:) = zero    
        end if
        
        
     end do
                 
     return
  end subroutine put_bubbles_in_bin_velocity_only  
!! ------------------------------------------------------------------------------------------------
  subroutine bubble_in_bin(i)
     !! This switches an element of b_inbin and resets bubble parameters
     integer(ikind),intent(in) :: i
     b_inbin(i) = .true.           !! Put in the bin
     b_freeindices(b_nfree+1) = i    !! Store the free index
     rb(i,:) = -1.0d2   !! Put it out of the way
     ub(i,:) = zero      !! Zero some velocities
     ubturb(i,:) = zero            
     u_l2b(i,:) = zero
     return
  end subroutine bubble_in_bin 
!! ------------------------------------------------------------------------------------------------
  subroutine calc_subgrid_fluctuations
     !! Follows method of Beuer & Hoppe 2017, which generalises Pozorski & Apte 2009
     !! N.B. rnd1 and rnd2 take several uses here.
     integer(ikind) :: i,j
     real(rkind) :: rnd1,rnd2,sigma_srs,T_srs,C0,urel_mag,b_coeff2
     real(rkind),dimension(dims) :: zeta,urel,utvec
     real(rkind),dimension(dims,dims) :: Wmat,Emat
     real(rkind) :: T_para,T_perp,urel_mag2,yplus

     !! Constants
     C0=one
     b_coeff2 = one !! beta in Beuer & Hoppe

     !$omp parallel do private(rnd1,rnd2,sigma_srs,T_srs,urel,urel_mag,T_para,T_perp &
     !$omp ,Wmat,Emat,urel_mag2,utvec,yplus)
     do i=1,n_bub 
        if(.not.b_inbin(i)) then          
           
           !! Standard deviation of subgrid fluctuations
           sigma_srs = sqrt(two*tke_l2b(i)/three)
           
           !! sub-grid fluctuation time-scales
           T_srs = C0*delta_F/max(sigma_srs,vsmall)
           
           !! Relative velocity and its magnitude
           urel(:) = vel_l2b(i,:) - ub(i,:)
           urel_mag2 = dot_product(urel,urel)  !squared
           urel_mag = sqrt(urel_mag2)
           
           !! Parallel and perpendicular timescales...
           rnd1 = b_coeff2*urel_mag2/max((two*tke_l2b(i)/three),vsmall)
           T_para = T_srs/sqrt(one + rnd1)
           T_perp = T_srs/sqrt(one + four*rnd1)
           
           !! Construct W, the covariance matrix
           rnd1 = sigma_srs*sqrt(one - exp(-two*dt/T_para))
           rnd2 = sigma_srs*sqrt(one - exp(-two*dt/T_perp))
           do j=1,dims
              Wmat(j,:) = (rnd1-rnd2)*urel(j)*urel(:)/max(urel_mag2,vsmall)
              Wmat(j,j) = Wmat(j,j) + rnd2
           end do           
           
           !! Construct E, the exponential of the drift matrix
           rnd1 = exp(-dt/T_para)
           rnd2 = exp(-dt/T_perp)
           do j=1,dims
              Emat(j,:) = (rnd1-rnd2)*urel(j)*urel(:)/max(urel_mag2,vsmall)
              Emat(j,j) = Emat(j,j) + rnd2
           end do
           
           !! A dims-dimensional normally distributed random vector (mean=0, sigma=1)                      
           do j=1,dims
              rnd1 = max(vsmall,rand());rnd2=rand()
              zeta(j) = cos(two*pi*rnd2)*sqrt(-two*log(rnd1))
           end do
                   
           !! Update turbulent velocity fluctuation
           utvec = ubturb(i,:)
           ubturb(i,:) = matmul(Emat(:,:),utvec(:)) + matmul(Wmat(:,:),zeta(:))      
              
        end if
        
     end do
     !$omp end parallel do
      
     return
  end subroutine calc_subgrid_fluctuations  
!! ------------------------------------------------------------------------------------------------
  subroutine bubble_traces
     !! If required, this routine stores the size,position, and relative velocity
     !! of a set of 499 bubbles at each time-step. Only used for some of the visualisations
     !! of the Fraga plume.
     integer(ikind) :: i,j,fnum
     
     !! Loop over each trace
     do i=1,499
     
        !! What is the file-name/number (starts at 200)
        fnum = 200+i
        
        !! Output the position of this bubble (if it exists, and is within n_bub)
        j=i_bubble_trace(i)
        if(.not.b_inbin(j))then
        if(n_bub>j) then
           write(fnum,*) time,radb(j)/Rhinze,rb(j,:),vel_l2b(j,:)+ubturb(j,:)-ub(j,:)
      
           !! Flush
           flush(fnum)

        end if
        end if
             

        !! Change the index if inbin (it should be the next new bubble)
!        if(b_inbin(j)) then
!           i_bubble_trace(i) = b_freeindices(b_nfree)
!        end if

     end do
          
      
     return
  end subroutine bubble_traces
!! ------------------------------------------------------------------------------------------------  
  subroutine bubble_breakup
     !! Determine whether a bubble breaks, and if so, the child sizes. Then break the bubble.
     integer(ikind) :: i,n_bub_tmp,n_breakup
     real(rkind) :: new_vol_frac
     real(rkind) :: l,Vb,Vmin
     
     !! Loop over all bubbles
     n_breakup = 0
     n_bub_tmp = n_bub
     do i=1,n_bub_tmp
        !! Ignore those in the bin
        if(.not.b_inbin(i)) then          

           !! Consider bubbles with deformation distance>bubble radius
           if(deformation_distance(i)>radb(i)) then !! Bubble breakup occurs..
              n_breakup = n_breakup + 1    
                           
              l = quarter*((12.0/8.2/We)**threefifths)*(mean_dissrate(i)**(-twofifths))/(two*radb(i))
              
!              if(l>2.0**(-6.0/45.0)) then
              if(l<0.90) then       !! Avoid issues with Lambda close to the critical value...
                 call breakup_spectrum(l,0.001_rkind,new_vol_frac) 
                 call split_bubble(i,new_vol_frac)   
              end if

           end if
          
        end if
     end do
     
     !! Output the breakup rate
     write(188,*) time,dble(n_breakup)/dt
     flush(188)
  
     return
  end subroutine bubble_breakup  
!! ------------------------------------------------------------------------------------------------
  subroutine breakup_spectrum(l,Vmin,nvf)
     !! Calculates the breakup spectrum based on a minimum volume ratio (no tiny tiny fragments)
     !! and Lambda. Follows the model in Martinez-Bazan 2010.
     real(rkind),intent(in) :: l    !! Lambda in Martinez-Bazan model
     real(rkind),intent(in) :: Vmin
     real(rkind),intent(out) :: nvf  !! New size fraction (in terms of volume)
     real(rkind) :: ll,Vmax,delta_V,Vi,sumV_tot,tmp
     integer(ikind) :: i,nn
     real(rkind),dimension(:),allocatable :: sumV
     logical :: keep_going
     
     !! Params
     Vmax = half
     ll = l**fivethirds
     
     !! Evaluate the normalisation of the p.d.f. and store numerical integral in sumV
     nn = 101
     allocate(sumV(nn));sumV=zero
     delta_V = (Vmax - Vmin)/dble(nn-1)
     do i=2,nn
        Vi = Vmin + (i-1)*(Vmax-Vmin)/dble(nn-1)   !! Vi is Volume of i-th interval
        tmp = (Vi**(-twothirds))*((one-Vi)**(-twothirds))*(Vi**twoninths-ll)*((one-Vi)**twoninths-ll) !! non-normalised pdf
        tmp = max(zero,tmp) !! pdf cannot be negative
        sumV(i) = sumV(i-1) + delta_V*tmp
     end do
     sumV_tot = sumV(nn)
     sumV = sumV/sumV_tot !! Normalise so sumV now contains the c.d.f.
        
     !! Pick a random number        
     tmp = rand()

     keep_going = .true.
     i=0
     !! Numerical inversion of the c.d.f.
     do while(keep_going)
        i=i+1
        if(sumV(i)>=tmp) then   !! When sumV>tmp, we have passed the point on inverse c.d.f.
           keep_going = .false.
           Vi = Vmin + (i-1)*(Vmax-Vmin)/dble(nn-1)   !! Vi is Volume of i-th interval
        end if
     end do             
     
     !! Set the new size fraction to Vi
     nvf = Vi
                        
  
     return
  end subroutine breakup_spectrum   
!! ------------------------------------------------------------------------------------------------  
  subroutine split_bubble(i,nvf)
     !! Takes the index of a parent bubble i, and the child volume fraction.
     !! creates 2 child bubbles, one with index i, and one with a different index.
     integer(ikind),intent(in) :: i
     real(rkind),intent(in) :: nvf
     integer(ikind) :: j
     real(rkind) :: tmp,r_c1,r_c2,r_p

     !! Radius of parent bubble
     r_p = radb(i)     
  
     !! radii of child bubbles
     r_c1 = r_p*nvf**onethird
     r_c2 = r_p*(one-nvf)**onethird          
    
     !! make child 1 ---------------------------------------    
     !! Use a new index or a free index?
     if(b_nfree==0) then 
        n_bub = n_bub + 1
        j=n_bub
     else
        j = b_freeindices(b_nfree) 
        b_nfree = b_nfree - 1
     end if
     
     !! Position, bin-status, size
     b_inbin(j) = .false.
     rb(j,:) = rb(i,:);ub(j,:) = ub(i,:)   
     radb(j) = r_c1
    
     !! Set some age flags and reset integral properties
     bubble_EoL(j) = 0;bubble_age(j) = zero;bubble_LE(j) = 1.0d10
     deformation_distance(j) = zero;mean_dissrate(i)=zero
    
     !! Add some random noise to the position ::
     tmp = rand();rb(j,1) = rb(j,1) + half*(tmp-half)*radb(i)
     tmp = rand();rb(j,2) = rb(j,2) + half*(tmp-half)*radb(i)
#ifdef dim3
     tmp = rand();rb(j,3) = rb(j,3) + half*(tmp-half)*radb(i)
#endif    

     !! Set the liquid velocity at bubble location (for use in virtual mass calculation)
     u_l2b(j,:) = u_l2b(i,:)
     ubturb(j,:)=ubturb(i,:)
     
     
     !! make child 2 --------------------------------------
     
     !! (using index of parent bubble)
     radb(i) = r_c2
     bubble_EoL(i) = 0;bubble_age(i) = zero;bubble_LE(i) = 1.0d10;        
     deformation_distance(i) = zero
     mean_dissrate(i)=zero
            
     !! Write relative sizes to file for plotting breakup spectrum
     write(181,*) nvf
     write(181,*) one-nvf
     flush(181)  

  
     return
  end subroutine split_bubble
!! ------------------------------------------------------------------------------------------------    
#endif  
!! ------------------------------------------------------------------------------------------------
end module bubble_evolution
