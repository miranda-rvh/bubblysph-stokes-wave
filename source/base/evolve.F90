module evolve
  !! This module contains routines to carry out a single step, calculate the time-step, and 
  !! evaluate any statistics or secondary properties we might want to output.
  use kind_parameters
  use common_parameter
  use common_2d
  use calculation_gradient
  use neighbour_finding
  use mirror_boundaries
  use freesurf
  use shifting
  use predictor_mod
  use sphtools
  use omp_lib
  use bubble_evolution
  use bubble_boundaries  
  use output
  implicit none

  private
  public :: step

contains
!! ------------------------------------------------------------------------------------------------
  subroutine step
    !! Carry out one step of the projection method.
    
    integer(ikind) :: i

    !! Some hard-coded control flags, largely useful for debugging.
    lagrangian = .true.                    !! Lagrangian or Eulerian.
    fick_shifting = .true.                 !! Use shifting (always yes!, unless debugging)
    output_everytime = .false.             !! output every time-step (useful for debugging)
    output_mirrors = .false.               !! Include mirror particles in output (for debugging)
    advection_terms = .false.              !! Include u_{ps}.grad(u).


    !! Pass position & velocity into o variables for storage
    allocate(rpo(3*n_par_fw,dims),upo(3*n_par_fw,dims))
    upo(1:n_par_fw,:) = up(1:n_par_fw,:)
    rpo(1:n_par_fw,:) = rp(1:n_par_fw,:)
    !! ------------------------------------------------------------------------
    
    !! Calculate the value of the time-step
    call set_time_step
    !! ------------------------------------------------------------------------            

    !! Initial advection to r*
    if(lagrangian)then 
       !$OMP PARALLEL DO
       do i=n_par_w+1, n_par_fw   !! N.B. counter - we don't move boundary particles...
          rp(i,:) = rpo(i,:) + up(i,:)*dt
       end do
       !$OMP END PARALLEL DO
    end if   
#ifdef bubbles
    !$omp parallel do
    do i=1,n_bub
       rb(i,:) = rb(i,:) + ub(i,:)*dt    !! Update bubble positions
    end do
    !$omp end parallel do
#endif        

   !! Modify boundary bits if using the wavemaker.
#ifdef wavemaker
    call wavemaker_velocity
    U_wmaker = 0.1*sin(pi*time)
    X_wmaker = X_wmaker + U_wmaker*dt
    b_node(1,1) = X_wmaker
    b_node(nb_patches,1) = X_wmaker
#endif 
    !! ------------------------------------------------------------------------       

    !! Put in bin, build boundaries and find neighbours
#ifdef bubbles
    call replace_escapees_b          !! Replace any bubbles which have left domain (if periodic)
    call put_bubbles_in_bin          !! Any bubbles in trouble - put them in the bubble-bin
#endif
    if(lagrangian.or.itime==1)then
       call put_in_bin               !! Any SPH particles in trouble - put them in the bin.
    end if
#ifdef bubbles    
    call create_mirror_bubbles       !! Create mirrors of the bubbles at boundaries
#endif     
    if(lagrangian.or.itime==1)then
       call create_mirror_particles  !! Create mirrors of the SPH particles at boundaries
       call particle_allocation      !! Build neighbour lists (both SPH and bubble)
    end if
    !! ------------------------------------------------------------------------        

    !! Evaluate the liquid volume fraction, and adjust SPH volumes and smoothing lengths
    allocate(a00(n_par_fwm));a00=zero    
    a00(1:n_par_fwm) = a0(1:n_par_fwm)   !! a00 is the volume fraction at the start of the time step 
#ifdef bubbles
    !! Update the bubble volume fractions...
    call update_bubble_density
#endif       
    !! ------------------------------------------------------------------------    

    !! Evaluate weights for gradient and Laplacian, identify free surface particles
    call kernel_calculation  !! Find weights for gradients        
    call locate_free_surf_part !! Free surface identification and various bits
    !! ------------------------------------------------------------------------    

        
    !! Calculate the shifting velocity
    allocate(ushift(n_par_fwm,dims));ushift=zero
    if(lagrangian)then
       if(fick_shifting)then
          call new_fick_shift
       else
          ushift(1:n_par_fw,:) = zero
          deallocate(grad_conc)          
       end if
    else  !! Eulerian SPH
       ushift(1:n_par_fw,:) = -up(1:n_par_fw,:)  
       deallocate(grad_conc)
    end if
    !! ------------------------------------------------------------------------        
    
    !! Prediction step: Calculate the viscous forces and update the velocities   
    call prediction_step               
    call mirror_velocities     !! re-do velocities in mirrors in prep for PPE
    !! ------------------------------------------------------------------------        
    
    !! Setup and solve the PPE
    call ppe_solve_in_house
    !! ------------------------------------------------------------------------        

    !! Calculate pressure gradient
    allocate(grad_p(n_par_fw,dims))   
    call calc_pressure_gradient(grad_p) !! see comments in subroutine for details, contains grad(a0*p) - g
    
    
    !! Calculate the final velocity
    !$OMP PARALLEL DO
    do i=n_par_w+1,n_par_fw
      up(i,:) = up(i,:) - grad_p(i,:)*dt/a0(i)  !! a0*u^{n+1} = a0*u^{*} - dt*grad_p 

      !! For spray particles - those with few neighbours - only gravity applies
      if(ij_count(i)<=2) up(i,:) = upo(i,:) + grav(:)*dt 
    enddo
    !$OMP END PARALLEL DO
    deallocate(grad_p)      
    
    !! Modify velocities as required
    call adjust_velocity_walls ! Remove component of velocity towards the wall. Required for low viscosity violent flows
#ifdef wavemaker
    call damping_zone !! Damping zone used by wavemaker calculations
#endif    
    !! ------------------------------------------------------------------------            

    !! Calculate vorticity if required for output
    if(output_this_step) call vorticity_calc
    !! ------------------------------------------------------------------------        
    
    !! Particle advection - move particles to final positions
    if(lagrangian)then
       !$OMP PARALLEL DO
       do i=n_par_w+1,n_par_fw    !! Don't move boundary particles
          rp(i,:) = rpo(i,:) + half*(upo(i,:) + up(i,:))*dt + ushift(i,:)*dt 
       end do
       !$OMP END PARALLEL DO
    end if
    !! ------------------------------------------------------------------------        
    
    !! Book-keeping: deallocate linking lists and other things
    if(lagrangian)then
       deallocate(irelation,dP_mp,vrelation)
       deallocate(ij_count,ij_link)
#ifdef bubbles
       deallocate(ijb_count,ijb_linearlink,ijb_firstindex)
       deallocate(irelation_b,vrelation_b)
       deallocate(bnfs)
#endif             
    end if
    deallocate(ushift)
    if(allocated(surf_norm)) deallocate(surf_norm)
    deallocate(rpo,upo,a00)
    !! ------------------------------------------------------------------------        

    return
  end subroutine step
!! ------------------------------------------------------------------------------------------------
  subroutine set_time_step
    !! This routine sets the time-step based on advective, viscous and acoustic criteria
    real(rkind) :: u_abs
    real(rkind) :: dt_viscous, dt_advection,dt_acoustic
    integer(ikind) :: i

    ! Constraint 1: CFL condition - coeff*h/max(|u|)
    umax = 1.0d-10
    do i=n_par_w+1,n_par_fw
       if(inbin(i)) cycle
       u_abs = sqrt(dot_product(up(i,:),up(i,:)))
       umax = max(umax,u_abs)
    end do
    dt_advection = dt_coef_advection*dx/umax

    ! Constraint 2: Viscous condition
    dt_viscous = dt_coef_viscous*Re*dx*dx

    ! Constraint 3: sound speed condition if required...
#ifdef isotcomp
    dt_acoustic = dt_coef_acoustic*Ma*dx
#else
    dt_acoustic = 1.0d10
#endif

    !! Choose the smallest       
    dt = min(dt_acoustic,min(dt_advection,dt_viscous))    
        
    return
  end subroutine set_time_step
!! ------------------------------------------------------------------------------------------------
  subroutine vorticity_calc
     !! Calculate the vorticity in preparation for outputting
     integer(ikind) :: i,j,k
     real(rkind) :: tmp,vort_tmp,vx,vy,vz

     allocate(vort(n_par_fwm));vort=zero

#if dim3
        !! 3D - magnitude of vorticity vector
        !$OMP PARALLEL DO PRIVATE(k,j,tmp,vz,vy,vx)
        do i=1,n_par_fw
           vz=zero;vy=zero;vx=zero
           do k=1,ij_count(i)
              j=ij_link(k,i)
              tmp = (up(i,2)-up(j,2))*ij_w_G(1,k,i) - (up(i,1)-up(j,1))*ij_w_G(2,k,i)
              vz = vz - tmp
              tmp = (up(i,3)-up(j,3))*ij_w_G(2,k,i) - (up(i,2)-up(j,2))*ij_w_G(3,k,i)
              vx = vx - tmp
              tmp = (up(i,1)-up(j,1))*ij_w_G(3,k,i) - (up(i,3)-up(j,3))*ij_w_G(1,k,i)                            
              vy = vy - tmp
           end do
           vort(i) = sqrt(vx*vx + vy*vy + vz*vz)
        end do
        !$OMP END PARALLEL DO
#else        
        !! 2D - Z component (other components are zero)
        !$OMP PARALLEL DO PRIVATE(k,j,tmp,vort_tmp)
        do i=1,n_par_fw
           vort_tmp=zero
           do k=1,ij_count(i)
              j=ij_link(k,i)
              tmp = (up(i,2)-up(j,2))*ij_w_G(1,k,i) - (up(i,1)-up(j,1))*ij_w_G(2,k,i)
              vort_tmp = vort_tmp - tmp
           end do
           vort(i) = vort_tmp
        end do
        !$OMP END PARALLEL DO
#endif    
  
     return
  end subroutine vorticity_calc  
!! ------------------------------------------------------------------------------------------------
  subroutine wavemaker_velocity
#ifdef wavemaker  
     integer(ikind) :: i, nspec, n,nn
     real(rkind) :: depth
     real(rkind) :: wavenumber, kd
     real(rkind) :: fp,  AN, om, ce, cg, f1, f2, df, ff
     real(rkind) :: Sp_sum, vel_sum, f(100), Sp(100), ak(100), apm(100)
     real(rkind) :: stroke2(100), tfocus, phi(100),ga,rr,gam,sig
     real(rkind) :: arg1,tmp1,tmp2
     real(rkind) :: target_wnumber

     !! Depth
     depth= 0.5d0!2.8d0 ! deep water depth

     !! Frequency of peak in spectrum
     target_wnumber = 2.0d0*pi
       
     !! Desired focal position
     x_target = 3.0d0

     !! Gravity
     ga = abs(grav(2))

     !! Desired amplitude at focal point    
     AN = 0.1d0!0.115d0!0.268d0 !Amplitude at focal point
     
     !! Find angular frequency of peak and peak frequency
     wavenumber = target_wnumber
     om = sqrt(wavenumber*ga*tanh(wavenumber*depth))
     fp = om/2.0d0/pi
     
!write(6,*) wavenumber,fp
     ce = om/wavenumber  !! Celerity
     cg = ce*(1d0+2d0*wavenumber*depth/sinh(2d0*wavenumber*depth))/2d0      !! Group celerity?

     !! upper and lower frequencies, number of intervals, and frequency increments
     f1=fp*0.5d0
     f2=fp*3.0d0
     nspec=100
     df= (f2-f1)/real(nspec)
     Sp_sum=0d0

     !! Loop over all frequency increments
     do n=1,nspec
        f(n) = f1 + df*real(n-1) + df/2d0
        ff = fp/f(n)
          
        !!Calculate power spectrum (JONSWAP)
        gam = 2.0d0!3.3d0

        if(f(n).gt.fp)then
           sig = 0.09d0
        else
           sig = 0.07d0
        endif
        rr = (f(n)-fp)/(fp*sig)
        rr = rr**2
        rr = exp(-rr/2d0)
          
        !! Sp is the power for the n-th component   
        Sp(n)= (ff**5d0)*exp(-1.25*ff**4d0)*(gam**rr)

        !! Sum of the power spectrum           
        Sp_sum = Sp_sum+ Sp(n)*df

        !! Determine the wavenumber of this frequency component
        om=2d0*pi*f(n)  !! angular frequency
        wavenumber=om*om/ga
        do nn=1,100
           wavenumber=om*om/(ga*tanh(wavenumber*depth))
        end do          
        ak(n)=wavenumber
     end do
     
     
     !! Loop over all components and calculate stroke of each component
     !! Note that for certain conditions, sums of hyperbolic trig get problematic --> inf/NaN
     do n=1,nspec

        !! amplitude of n-th component   (eqn 14 in LSR)
        apm(n) = AN*Sp(n)*df/Sp_sum 

        !! Calculate stroke amplitude for n-th component (eqn 17 in LSR)
        stroke2(n) = (sinh(2d0*ak(n)*depth)+2d0*ak(n)*depth)/2d0
        stroke2(n) = stroke2(n)/(cosh(2d0*ak(n)*depth)-1d0)
        stroke2(n) = apm(n)*stroke2(n)
!write(6,*) apm(n),stroke2(n)
     end do

     !! Calculate the focal time
     tfocus= 2d0*x_target/cg  
      
     !! Sum the velocities of individual components to obtain the paddle velocity
     vel_sum = 0d0
     do n=1,nspec
        om= 2d0*pi*f(n)
        phi(n)=om*tfocus-ak(n)*x_target     !! phi shifts each component so they all focus at x_target,tfocus
        vel_sum = vel_sum + stroke2(n)*om*cos(-om*time+phi(n))
     end do

     U_wmaker = vel_sum
!write(6,*) "time",time,"U_wmaker",U_wmaker     

     !! Write the position and velocity to output
     write(174,*) time,X_wmaker + U_wmaker*dt,U_wmaker
     flush(174)
     
#endif
     return
  end subroutine wavemaker_velocity
!! ------------------------------------------------------------------------------------------------
  subroutine damping_zone
#ifdef wavemaker
     real(rkind) :: damp_cutoff,pos_along_domain,pos_within_damping,damping_coeff  
     integer(ikind) :: i

     !! Damp the final 1-X proportion of the domain
     damp_cutoff = 0.75
     
     !! Loop over particles
     !$omp parallel do private(pos_along_domain)
     do i=1,n_par_fw
     
        !! Find position along the domain (normalised)
        pos_along_domain = (rp(i,1)-xb_min)/(xb_max-xb_min)  !! Proportion of distance across domain       
        
        !! Damp velocities if required
        if(pos_along_domain>=damp_cutoff) then !! Final X% of domain, damp
        
           !! pos_within_damping runs linearly from zero to one between start and end of damping zone
           pos_within_damping = (pos_along_domain - damp_cutoff)/(one-damp_cutoff)

           !! Some exponential function for damping
!           damping_coeff = exp(-3.0d0*pos_within_damping)
           damping_coeff = half - half*erf(5.0d0*(pos_within_damping-half))
!           damping_coeff = exp(-10.0d0*pos_within_damping**4.0d0)
           
           !! Scale the velocity
           up(i,:) = up(i,:)*damping_coeff
                     
        end if
     end do
     !$omp end parallel do
#endif
     return
  end subroutine damping_zone
!! ------------------------------------------------------------------------------------------------  
end module evolve
