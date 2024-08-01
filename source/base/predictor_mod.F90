module predictor_mod
  !! This module calculates the intermediate velocity u* = u^n + dt*(viscous + interphase)
  !! It also performs the following functions:
  !!   > Evaluate effective viscosity based on sub-resolution turbulence closure model
  !!   > Evaluate velocity gradients, which are used in other routines
  !!   > Evaluate dissipation rate used in bubble calculations, and total dissipation 
  !!     rate which is output to a file
  !!   > Apply any external body forces (with a hard-coded switch)
  
  !! Set this flag to determine the closure model.
#define srs_closure 3
  !! srs_closure = 0   - standard smagorinsky
  !! srs_closure = 1   - Germano, with Shepard averaging  
  !! srs_closure = 2   - Germano, with Lagrangian averaging
  !! srs_closure = 3   - Mixed-Scale-Model (by Sagaut)
  use kind_parameters
  use common_parameter
  use common_2d
  use sphtools
  use calculation_gradient
  use mirror_boundaries
  use omp_lib

  implicit none

  real(rkind),dimension(:,:),allocatable :: acc
  real(rkind),dimension(dims) :: uij,acc_tmp,norm_temp
  real(rkind),dimension(:),allocatable :: norm_srt,mu_srs
  real(rkind),dimension(:,:,:),allocatable :: srt  

  private
  public :: prediction_step

contains
!! ------------------------------------------------------------------------------------------------
  subroutine prediction_step
    use bubble_evolution  
    integer(ikind) :: i,j,k,kk
    real(rkind) :: mu_ij,a0_ij,x,y,z,fpios,pios,rad,qq,wqq,temp_L,temp
    real(rkind) :: aatau,tmp_acc,phase_coeff,tmp_diss,tmp_vol
    real(rkind),dimension(:),allocatable :: mu_eff,diss_rate_tmp
    real(rkind),dimension(dims) :: rij
    
    real(rkind) :: diss,diss2

    !! Calculate the velocity gradient 
    !! (technically only needed when bubbles or srs or advection_terms or not(lagrangian))
    allocate(grad_vel(dims,dims,n_par_fwm))
    call grad_operator(grad_vel(1,:,:),up(:,1))
    call grad_operator(grad_vel(2,:,:),up(:,2))
#if dim3
    call grad_operator(grad_vel(3,:,:),up(:,3))
#endif   


#ifdef srs
    !! Calculate strain-rate tensor norm
    allocate(norm_srt(n_par_fwm),srt(dims,dims,n_par_fwm))
    !$OMP PARALLEL DO 
    do i=1,n_par_fw
       srt(:,:,i) = half*(grad_vel(:,:,i) + transpose(grad_vel(:,:,i)))
#if dim3
       norm_srt(i) = sqrt(two*(srt(1,1,i)**two + two*srt(1,2,i)**two + srt(2,2,i)**two &
                  + two*srt(1,3,i)**two + two*srt(2,3,i)**two + srt(3,3,i)**two))
#else
       norm_srt(i) = sqrt(two*(srt(1,1,i)**two + srt(1,2,i)**two + srt(2,1,i)**two + srt(2,2,i)**two))
#endif       

    end do
    !$OMP END PARALLEL DO
    !! Copy strain rate tensor and norm to mirrors (no inversion...)
    !$OMP PARALLEL DO PRIVATE(i)
    do j=n_par_fw+1,n_par_fwm
       i=irelation(j)
       srt(:,:,j) = srt(:,:,i);norm_srt(j)=norm_srt(i)
    end do 
    !$OMP END PARALLEL DO
    
    !! Obtain the sub-resolution-scale viscosity from the appropriate model
    !! These routines are defined lower down the module
#if srs_closure==0    
    call static_smag      
#elif srs_closure==1
    call dynamic_smag_local      
#elif srs_closure==2
    call dynamic_smag_lagrangian
#elif srs_closure==3
    call mixed_scale_model
#endif    

   
    !! Calculate the sub-resolution turbulent kinetic energy, only used for bubble fluctuations
#ifdef bubbles
    call calc_subgrid_tke
#endif       

    !! Calculate the dissipation rate and effective viscosity, EXCLUDING bubble induced contributions.
    allocate(diss_rate_tmp(n_par_fwm))
    allocate(mu_eff(n_par_fwm));mu_eff=zero               
    diss = zero
    !$OMP PARALLEL DO reduction(+:diss)
    do i=1,n_par_fw
       mu_eff(i) = one + mu_srs(i)
       
       !! SGS dissipation rate
       diss_rate_tmp(i)=(mu_srs(i)/Re)*norm_srt(i)**two

       !! Summing to get global measure of dissipation
       diss = diss + a00(i)*vol(i)*(mu_eff(i)/Re)*norm_srt(i)**two                    
    end do
    !$OMP END PARALLEL DO
    write(198,*) time,diss !<--- writing volume integrated dissipation rate to file
    flush(198)
    deallocate(norm_srt,srt,mu_srs)    
    
    
    !! Copy effective viscosity and sugrid dissipation rate to neighbours
    !$OMP PARALLEL DO PRIVATE(i)
    do j=n_par_fw+1,n_par_fwm 
       i=irelation(j);mu_eff(j)=mu_eff(i);diss_rate_tmp(j) = diss_rate_tmp(i)
    end do
    !$OMP END PARALLEL DO
    
    !! Smooth the dissipation rate for use in predicting bubble entrainment rates
    eps_srs = zero
    !$OMP PARALLEL DO PRIVATE(k,j,rij,rad,qq,wqq,tmp_vol,tmp_diss)
    do i=1,n_par_fw
       tmp_diss = zero;
       tmp_vol = zero
       do k=1,ij_count(i)
          j=ij_link(k,i)
          
          !! Smooth eps_srs for use in bubble model
          rij = rp(i,:)-rp(j,:);rad = sqrt(dot_product(rij,rij));qq=rad/h(i)
          wqq = wab(qq)
          tmp_diss = tmp_diss + wqq*diss_rate_tmp(j)
          tmp_vol = tmp_vol + wqq
         
       end do
       eps_srs(i) = tmp_diss/max(tmp_vol,epsilon(tmp_vol))
a_out(i) = eps_srs(i)       
    end do 
    !$omp end parallel do    
    deallocate(diss_rate_tmp)           

    !! Copy dissipation rate to mirrors
    !$OMP PARALLEL DO PRIVATE(i)
    do j=n_par_fw+1,n_par_fwm 
          i=irelation(j);eps_srs(j)=eps_srs(i)
    end do
    !$OMP END PARALLEL DO    
#else
    !! If no srs model, then zero the dissipation rate, and set the effective viscosity to unity
    eps_srs = zero
    allocate(mu_eff(n_par_fwm));mu_eff=zero               
    mu_eff(1:n_par_fwm)=one
#endif
   
   
#ifdef bubbles
    !! Update bubble velocities, and calculate inter-phase momentum exchange and bubble induced viscosity
    allocate(nu_BI(n_par_fw));nu_BI=zero
    allocate(M0(n_par_fw,dims));M0=zero
    call calculate_bubble_velocities   
    
    !! Adjust effective viscosity
    !$omp parallel do
    do i=1,n_par_fw
       mu_eff(i) = mu_eff(i) + nu_BI(i)
    end do
    !$omp end parallel do
    !! Copy to mirrors (again)
    !$omp parallel do private(i)
    do j=n_par_fw+1,n_par_fwm 
          i=irelation(j);mu_eff(j)=mu_eff(i)
    end do
    !$omp end parallel do
    deallocate(nu_BI)
#endif   
   
    !! Calculate acceleration due to viscous term
    allocate(acc(n_par_fw,dims))
    acc(:,:)=zero
    diss = zero
    !$OMP PARALLEL DO PRIVATE(k,j,uij,mu_ij,acc_tmp,a0_ij,rij,rad,qq,wqq) &
    !$omp reduction(+:diss)
    do i=1,n_par_fw
       acc_tmp(:) = zero
       tmp_diss = zero
       do k=1,ij_count(i)
          j=ij_link(k,i)
          uij(:) = up(i,:) - up(j,:)         
          mu_ij = mu_eff(i)*mu_eff(j)/(mu_eff(i) + mu_eff(j)) ! Harmonic mean (dropped a factor of 2)
          a0_ij = a00(i)*a00(j)/(a00(i) + a00(j)) ! Harmonic mean of vol frac (dropped a factor of 2)
      
          acc_tmp(:) = acc_tmp(:) + a0_ij*mu_ij*uij(:)*ij_w_L(k,i)     
       end do
   
       !! Finalise viscous acceleration. Nb. the 4 comes from the 2 Harmonic means, the 1/Re from the non-dimensionalisation.
       acc(i,:) = 4.0d0*acc_tmp(:)/Re 
            
       !! A secondary calculation of dissipation rate, just for outputting
       diss = diss + (4.0d0/Re)*vol(i)*sqrt(dot_product(acc_tmp(:),acc_tmp(:)))  
       

    end do
    !$OMP END PARALLEL DO
    write(199,*) time,diss   !<-- another file to store this second calculation of volume integrated dissipation rate
    flush(199)
    deallocate(mu_eff)       
   
    !! Add external forcing if required
    !! These routines are defined lower down
    if(.false.) call kolmogorov_forcing
    if(.false.) call linear_forcing
    if(.false.) call antuono_forcing        

    !! Add the advection terms if required
    if(advection_terms.or..not.lagrangian) then 
       !$OMP PARALLEL DO
       do i=1,n_par_fw             !! N.B. should be additional grad alpha term in here, but we will neglect it.
          acc(i,1) = acc(i,1) + dot_product(ushift(i,:),grad_vel(1,:,i))         
          acc(i,2) = acc(i,2) + dot_product(ushift(i,:),grad_vel(2,:,i))          
#if dim3
          acc(i,3) = acc(i,3) + dot_product(ushift(i,:),grad_vel(3,:,i))
#endif          
       end do
       !$OMP END PARALLEL DO
    end if     
#ifdef bubbles
    !! Add the interphase momentum exchange
    !$omp parallel do
    do i=1,n_par_fw
       acc(i,:) = acc(i,:) + M0(i,:)
    end do
    !$omp end parallel do
    deallocate(M0)
#endif    

    !! Update the velocity to find U*  
    !$omp parallel do
    do i=1+n_par_w,n_par_fw
       !!a0^n+1 u* = a0^n u^n + dt*acc    
       up(i,:) = (a00(i)*up(i,:) + acc(i,:)*dt)/a0(i)                  
    end do
    !$omp end parallel do
     

    if(allocated(grad_vel)) deallocate(grad_vel)
    deallocate(acc)
  end subroutine prediction_step
!! ------------------------------------------------------------------------------------------------
  subroutine static_smag
     !! Static smagorinsky model
     use sphtools
     integer(ikind) :: i
     real(rkind) :: C_smag,yplus


     allocate(mu_srs(n_par_fw))
     C_smag=0.18!12 !! Limit max C_smag=half

     !$OMP PARALLEL DO PRIVATE(yplus)
     do i=1,n_par_fw
        mu_srs(i) = Re*((C_smag*delta_F)**two)*norm_srt(i)        
               
     end do
     !$OMP END PARALLEL DO                 

     return
  end subroutine static_smag 
!! ------------------------------------------------------------------------------------------------
  subroutine dynamic_smag_local
     !! Dynamic smagorinsky model with local averaging
     !! Follows Kirkby 2014, Vremen, Geurts & Kuerten 1997
     use sphtools
     integer(ikind) :: i,j,k
     real(rkind),dimension(:,:,:),allocatable :: uu,srt_norm_srt
     real(rkind),dimension(dims) :: u_filt,rij
     real(rkind),dimension(dims,dims) :: uu_filt,srt_filt,sns_filt,L_tens,M_tens
     real(rkind),dimension(dims,dims) :: ufilt_ufilt
     real(rkind) :: norm_srt_filt
     real(rkind) :: qq,rad,temp,c_tmp,a_coef,delta_G,C_smag
     real(rkind),dimension(:),allocatable :: c_smag_sq

     allocate(mu_srs(n_par_fw))
     delta_G = sup_size*two  !! Top hat test filter width of stencil
     a_coef = delta_G/delta_F  !! Ratio of test filter to SPH (grid) filter widths...

     !! Step 1: allocate variables to be filtered
     allocate(uu(n_par_fwm,dims,dims),srt_norm_srt(n_par_fwm,dims,dims))
     uu(1:n_par_fwm,1,1)=up(1:n_par_fwm,1)**two;uu(1:n_par_fwm,1,2)=up(1:n_par_fwm,1)*up(1:n_par_fwm,2)
     uu(1:n_par_fwm,2,1)=uu(1:n_par_fwm,1,2);uu(1:n_par_fwm,2,2)=up(1:n_par_fwm,2)**two
#if dim3
     uu(1:n_par_fwm,1,3)=up(1:n_par_fwm,1)*up(1:n_par_fwm,3);uu(1:n_par_fwm,2,3)=up(1:n_par_fwm,2)*up(1:n_par_fwm,3)
     uu(1:n_par_fwm,3,3)=up(1:n_par_fwm,3)*up(1:n_par_fwm,3);uu(1:n_par_fwm,3,1)=uu(1:n_par_fwm,1,3)
     uu(1:n_par_fwm,3,2)=uu(1:n_par_fwm,2,3)
#endif     
     !$OMP PARALLEL DO
     do i=1,n_par_fwm
         srt_norm_srt(i,:,:)=norm_srt(i)*srt(:,:,i)
     end do
     !$OMP END PARALLEL DO

     !! Step 2: for each particle, find the filtered vars...
     allocate(c_smag_sq(n_par_fwm));c_smag_sq=zero
     !$OMP PARALLEL DO PRIVATE(k,j,rij,rad,qq,temp,u_filt,c_tmp,uu_filt,srt_filt, &
     !$OMP norm_srt_filt,sns_filt,M_tens,L_tens,ufilt_ufilt)
     do i=1,n_par_fw
        u_filt=up(i,:);c_tmp=one!wab(zero)
        uu_filt=uu(i,:,:)
        srt_filt=srt(:,:,i);norm_srt_filt=norm_srt(i)
        sns_filt=srt_norm_srt(i,:,:)
        do k=1,ij_count(i)
           j=ij_link(k,i)
           rij = rp(i,:)-rp(j,:);rad=sqrt(dot_product(rij,rij));qq=rad/h(i)
           temp = wab(qq)
           u_filt(:) = u_filt(:) + temp*up(j,:)
           uu_filt(:,:) = uu_filt(:,:) + temp*uu(j,:,:)
           c_tmp = c_tmp + temp
           srt_filt(:,:)=srt_filt(:,:) + temp*srt(:,:,j)
           norm_srt_filt = norm_srt_filt + temp*norm_srt(j)
           sns_filt(:,:)=sns_filt(:,:) + temp*srt_norm_srt(j,:,:)                 
        end do
        u_filt = u_filt/c_tmp
        uu_filt = uu_filt/c_tmp
        srt_filt = srt_filt/c_tmp
        norm_srt_filt = norm_srt_filt/c_tmp
        sns_filt = sns_filt/c_tmp

        !! Step 3: construct the L and M tensors
        ufilt_ufilt(1,1)=u_filt(1)**two;ufilt_ufilt(1,2)=u_filt(1)*u_filt(2)
        ufilt_ufilt(2,1)=u_filt(2)*u_filt(1);ufilt_ufilt(2,2)=u_filt(2)**two
#if dim3
        ufilt_ufilt(1,3)=u_filt(1)*u_filt(3);ufilt_ufilt(2,3)=u_filt(2)*u_filt(3)
        ufilt_ufilt(3,3)=u_filt(3)*u_filt(3);ufilt_ufilt(3,1)=u_filt(3)*u_filt(1)
        ufilt_ufilt(3,2)=u_filt(3)*u_filt(2)
#endif        
        L_tens(:,:) = ufilt_ufilt(:,:) - uu_filt(:,:)
        M_tens(:,:) = a_coef*norm_srt_filt*srt_filt(:,:) - sns_filt(:,:)

        !! Step 4: Set C_smag
        temp=zero;rad=zero
        do j=1,dims
           do k=1,dims
              temp = temp + M_tens(j,k)*L_tens(j,k)      !! temp is L*M
              rad = rad + M_tens(j,k)**two               !! rad is M*M
           end do
        end do
        C_smag_sq(i) = max(zero,temp/(two*rad*delta_G**two))
     end do
     !$OMP END PARALLEL DO
     deallocate(uu,srt_norm_srt)

     !! Local averaging of C_smag
     !$OMP PARALLEL DO PRIVATE(j)
     do i=n_par_fw+1,n_par_fwm
        j=irelation(i);c_smag_sq(i)=c_smag_sq(j)
     end do
     !$OMP END PARALLEL DO
     !$OMP PARALLEL DO PRIVATE(k,j,rij,rad,qq,c_tmp,C_smag)
     do i=1,n_par_fw
        c_tmp=wab(zero)
        C_smag=zero
        do k=1,ij_count(i)
           j=ij_link(k,i)
           rij = rp(i,:)-rp(j,:);rad=sqrt(dot_product(rij,rij))
           qq = wab(rad/h(i))
           c_tmp = c_tmp + qq   !! Kernel sum
           C_smag = C_smag + c_smag_sq(j)*qq  !! SPH smoothing of c_smag_sq
        end do
        C_smag=min(sqrt(C_smag/c_tmp),half) !! Limit max C_smag=half
        mu_srs(i) = Re*((C_smag*delta_F)**two)*norm_srt(i)        
     end do
     !$OMP END PARALLEL DO
     deallocate(c_smag_sq)

     return
  end subroutine dynamic_smag_local
!! ------------------------------------------------------------------------------------------------
  subroutine dynamic_smag_lagrangian
     !! Dynamic smagorinsky model with Lagrangian averaging
     !! C. Meneveau, T.S. Lund and W.H. Cabot J. Fluid Mech., 319 (1996), pp. 353-385
     use sphtools
     integer(ikind) :: i,j,k
     real(rkind),dimension(:,:,:),allocatable :: uu,srt_norm_srt
     real(rkind),dimension(dims) :: u_filt,rij
     real(rkind),dimension(dims,dims) :: uu_filt,srt_filt,sns_filt,L_tens,M_tens
     real(rkind),dimension(dims,dims) :: ufilt_ufilt
     real(rkind) :: norm_srt_filt
     real(rkind) :: qq,rad,temp,c_tmp,a_coef,delta_G,unovT,C_smag
     real(rkind),dimension(:),allocatable :: I_LM0,I_MM0
     real(rkind),dimension(dims) :: gradILM,gradIMM

     allocate(mu_srs(n_par_fw))
     delta_G = sup_size*two  !! Top hat test filter width of stencil
     a_coef = delta_G/delta_F  !! Ratio of test filter to SPH (grid) filter widths...

     !! If advective terms, prepare to find I_LM,I_MM gradients...
     if(advection_terms.or..not.lagrangian)then
!        allocate(I_LM0(n_par_fwm),I_MM0(n_par_fwm))
!        I_LM0(1:n_par_fw)=I_LM(1:n_par_fw);I_MM0(1:n_par_fw)=I_MM(1:n_par_fw)
!        !$OMP PARALLEL DO PRIVATE(j)
!        do i=n_par_fw+1,n_par_fwm
!           j=irelation(i)
!           I_LM0(i)=I_LM(j);I_MM0(i)=I_MM(j)
!        end do
!        !$OMP END PARALLEL DO
     end if

     !! allocate variables to be filtered
     allocate(uu(n_par_fwm,dims,dims),srt_norm_srt(n_par_fwm,dims,dims))
     uu(1:n_par_fwm,1,1)=up(1:n_par_fwm,1)**two;uu(1:n_par_fwm,1,2)=up(1:n_par_fwm,1)*up(1:n_par_fwm,2)
     uu(1:n_par_fwm,2,1)=uu(1:n_par_fwm,1,2);uu(1:n_par_fwm,2,2)=up(1:n_par_fwm,2)**two
#if dim3
     uu(1:n_par_fwm,1,3)=up(1:n_par_fwm,1)*up(1:n_par_fwm,3);uu(1:n_par_fwm,2,3)=up(1:n_par_fwm,2)*up(1:n_par_fwm,3)
     uu(1:n_par_fwm,3,3)=up(1:n_par_fwm,3)*up(1:n_par_fwm,3);uu(1:n_par_fwm,3,1)=uu(1:n_par_fwm,1,3)
     uu(1:n_par_fwm,3,2)=uu(1:n_par_fwm,2,3)
#endif       
     !$OMP PARALLEL DO
     do i=1,n_par_fwm
         srt_norm_srt(i,:,:)=norm_srt(i)*srt(:,:,i)
     end do
     !$OMP END PARALLEL DO

     !! Step 2: for each particle, find the filtered vars...
     !$OMP PARALLEL DO PRIVATE(k,j,rij,rad,qq,temp,u_filt,c_tmp,uu_filt,srt_filt, &
     !$OMP norm_srt_filt,sns_filt,M_tens,L_tens,ufilt_ufilt,gradILM,gradIMM,unovT,C_smag) 
     do i=1,n_par_fw
        u_filt=up(i,:);c_tmp=one!wab(zero)
        uu_filt=uu(i,:,:)
        srt_filt=srt(:,:,i);norm_srt_filt=norm_srt(i)
        sns_filt=srt_norm_srt(i,:,:)
        do k=1,ij_count(i)
           j=ij_link(k,i)
           rij = rp(i,:)-rp(j,:);rad=sqrt(dot_product(rij,rij));qq=rad/h(i)
           temp = wab(qq)
           u_filt(:) = u_filt(:) + temp*up(j,:)
           uu_filt(:,:) = uu_filt(:,:) + temp*uu(j,:,:)
           c_tmp = c_tmp + temp
           srt_filt(:,:)=srt_filt(:,:) + temp*srt(:,:,j)
           norm_srt_filt = norm_srt_filt + temp*norm_srt(j)
           sns_filt(:,:)=sns_filt(:,:) + temp*srt_norm_srt(j,:,:)  
        end do
        u_filt = u_filt/c_tmp
        uu_filt = uu_filt/c_tmp
        srt_filt = srt_filt/c_tmp
        norm_srt_filt = norm_srt_filt/c_tmp
        sns_filt = sns_filt/c_tmp

        !! Step 3: construct the L and M tensors
        ufilt_ufilt(1,1)=u_filt(1)**two;ufilt_ufilt(1,2)=u_filt(1)*u_filt(2)
        ufilt_ufilt(2,1)=u_filt(2)*u_filt(1);ufilt_ufilt(2,2)=u_filt(2)**two
#if dim3
        ufilt_ufilt(1,3)=u_filt(1)*u_filt(3);ufilt_ufilt(2,3)=u_filt(2)*u_filt(3)
        ufilt_ufilt(3,3)=u_filt(3)*u_filt(3);ufilt_ufilt(3,1)=u_filt(3)*u_filt(1)
        ufilt_ufilt(3,2)=u_filt(3)*u_filt(2)
#endif          
        L_tens(:,:) = ufilt_ufilt(:,:) - uu_filt(:,:)
        M_tens(:,:) = a_coef*norm_srt_filt*srt_filt(:,:) - sns_filt(:,:)

        !! L_ij*M_ij = temp; M_ij*M_ij = rad
        temp=zero;rad=zero
        do j=1,dims
           do k=1,dims
              temp = temp + M_tens(j,k)*L_tens(j,k)    !! temp is L*M
              rad = rad + M_tens(j,k)**two              !! rad is M*M
           end do
        end do
        rad = rad*two*delta_G**two

        
!        unovT = (two*abs(I_LM(i))**0.25)/delta_G
        unovT = (abs(I_LM(i)*I_MM(i))**0.125d0)/delta_G/1.5d0
        I_LM(i) = (I_LM(i) + dt*unovT*temp)/(one+dt*unovT)  !! Implicit backwards Euler
        I_MM(i) = (I_MM(i) + dt*unovT*rad)/(one+dt*unovT)

        !! gradILM,gradIMM terms
        if(advection_terms)then        !! Explicit advection terms...
!           gradILM=zero;gradIMM=zero
!           do k=1,ij_count(i)
!              j=ij_link(k,i)
!              gradILM=gradILM + (I_LM0(j)-I_LM0(i))*ij_w_G(:,k,i)               
!              gradIMM=gradIMM + (I_MM0(j)-I_MM0(i))*ij_w_G(:,k,i)               
!           end do        
!           I_LM(i) = I_LM(i) + dt*dot_product(ushift(i,:),matmul(kgcm(:,:,i),gradILM))
!           I_MM(i) = I_MM(i) + dt*dot_product(ushift(i,:),matmul(kgcm(:,:,i),gradIMM))
        end if        

        !! Step 4: Set C_smag and mu_srs
        C_smag = sqrt(max(I_LM(i),zero)/max(I_MM(i),vsmall))
        mu_srs(i) = Re*((C_smag*delta_F)**two)*norm_srt(i)        
     end do
     !$OMP END PARALLEL DO
     deallocate(uu,srt_norm_srt)

!     if(advection_terms) deallocate(I_LM0,I_MM0)

     return
  end subroutine dynamic_smag_lagrangian
!! ------------------------------------------------------------------------------------------------
  subroutine mixed_scale_model
     !! Mixed scale model, follows Sagaut, Lubin PhD thesis and Lubin 2006
     use sphtools
     integer(ikind) :: i,j,k
     real(rkind),dimension(dims) :: u_filt,rij
     real(rkind) :: qq,rad,temp,c_tmp,a_coef,qsquared,C_msm

     allocate(mu_srs(n_par_fw))
     a_coef = half !! Weighting for geometric average
     C_msm = 0.06   !! mixed scale model coefficient
     
     !$OMP PARALLEL DO PRIVATE(k,j,rij,rad,qq,temp,u_filt,c_tmp,qsquared)
     do i=1,n_par_fw
        !! Filter the velocity with a Shephard filter
        u_filt=zero;c_tmp=zero
        do k=1,ij_count(i)
           j=ij_link(k,i)
           rij = rp(i,:)-rp(j,:);rad=sqrt(dot_product(rij,rij));qq=rad/h(i)
           temp = wab(qq)
           u_filt(:) = u_filt(:) + temp*up(j,:)
           c_tmp = c_tmp + temp
        end do
        u_filt = u_filt/c_tmp
 
        !! find just small scales
        u_filt = up(i,:) - u_filt(:)
        
        !! TKE at cut-off
        qsquared = half*dot_product(u_filt,u_filt)

        !! MSM model
        mu_srs(i) = Re*C_msm*(delta_F**(one+a_coef))*(norm_srt(i)**half*a_coef)*qsquared**(half-half*a_coef)

     end do
     !$OMP END PARALLEL DO
     return
  end subroutine mixed_scale_model  
!! ------------------------------------------------------------------------------------------------
  subroutine calc_subgrid_tke
     !! Estimate the subgrid TKE from the double filtered velocity.
     use sphtools
     integer(ikind) :: i,j,k
     real(rkind),dimension(dims) :: u_filt,rij
     real(rkind) :: qq,rad,temp,c_tmp
    
     allocate(tke(n_par_fwm));tke=zero
    
     !$OMP PARALLEL DO PRIVATE(k,j,rij,rad,qq,temp,u_filt,c_tmp)
     do i=1,n_par_fw
        !! Filter the velocity with a Shephard filter
        u_filt=zero;c_tmp=zero
        do k=1,ij_count(i)
           j=ij_link(k,i)
           rij = rp(i,:)-rp(j,:);rad=sqrt(dot_product(rij,rij));qq=rad/h(i)
           temp = wab(qq)
           u_filt(:) = u_filt(:) + temp*up(j,:)
           c_tmp = c_tmp + temp
        end do
        u_filt = u_filt/c_tmp
 
        !! subtract from the (implicitly filtered) velocity
        u_filt = up(i,:) - u_filt(:)
        
        !! TKE at cut-off
        tke(i) = half*dot_product(u_filt,u_filt)

     end do
     !$OMP END PARALLEL DO
     
     !$omp parallel do private(j)
     do i=n_par_fw+1,n_par_fwm
        j=irelation(i)
        tke(i)=tke(j)
     end do
     !$omp end parallel do
     
     return
  end subroutine calc_subgrid_tke  
!! ------------------------------------------------------------------------------------------------
  subroutine kolmogorov_forcing
     integer(ikind) :: i
     
     !! This will need modifying depending on the domain size. May need magnitude modifying too.
     !$omp parallel do
     do i=1,n_par_fw
        acc(i,1) = acc(i,1) + sin(four*rp(i,2))   !! Kolmogorov forcing!!!
     end do
     !$omp end parallel do
     
     return
  end subroutine kolmogorov_forcing
!! ------------------------------------------------------------------------------------------------
  subroutine linear_forcing
     !! Routine to add linear forcing term proportional to the velocity, as in XX (need to find ref)
     integer(ikind) :: i
     real(rkind),dimension(dims) :: meanvel
     real(rkind) :: totenergy,totvol     
       
     !! Evaluate the total energy in the domain, and the mean velocity
     !$omp parallel do reduction(+:totenergy,meanvel,totvol)   
     do i=1,n_par_fw
        totenergy = totenergy + vol(i)*dot_product(up(i,:),up(i,:))
        meanvel = meanvel + vol(i)*up(i,:)
        totvol = totvol + vol(i)
     end do
     !$OMP END PARALLEL DO
     meanvel = meanvel/totvol
   
     !! Add forcing to acceleration term
     !! N.B. the mean velocity is removed from the forcing, to prevent drift.
     !$omp parallel do
     do i=1,n_par_fw
        acc(i,:) = acc(i,:) + one*totvol*(up(i,:)-meanvel)/totenergy
     end do
     !$omp end parallel do       
         
     return
  end subroutine linear_forcing
!! ------------------------------------------------------------------------------------------------
  subroutine antuono_forcing
     !! Add the forcing used in Antuono et al. Physics of Fluids 33 (2021) 015102. 
     !! doi:doi:10.1063/5.0034568.
     !! Results in a generalised Beltrami flow at t=1
     integer(ikind) :: i
     real(rkind) :: x,y,z
     real(rkind) :: ramp,fpios,pios
     real(rkind),dimension(dims) :: acc_tmp

     if(time<0.1d0) then   
        ramp = 10.0*time
     else if(time<0.9d0) then
        ramp = one
     else if(time<one) then
        ramp = 10.0d0*(one-time)
     else
        ramp = zero
     end if
     ramp = -(four*sqrt(two)/three/sqrt(three))*ramp*three*16.0d0*pi*pi*1.25d-3!/Re
     fpios=five*pi/6.0d0;pios=pi/6.0d0
     !$omp parallel do private(acc_tmp,x,y,z)
     do i=1,n_par_fw
         x=two*pi*rp(i,1);y=two*pi*rp(i,2);z=two*pi*rp(i,3)
         acc_tmp(1) = sin(two*x-fpios)*cos(two*y-pios)*sin(two*z) - cos(two*z-fpios)*sin(two*x-pios)*sin(two*y)
         acc_tmp(2) = sin(two*y-fpios)*cos(two*z-pios)*sin(two*x) - cos(two*x-fpios)*sin(two*y-pios)*sin(two*z)
         acc_tmp(3) = sin(two*z-fpios)*cos(two*x-pios)*sin(two*y) - cos(two*y-fpios)*sin(two*z-pios)*sin(two*x) 
         acc(i,:) = acc(i,:) + ramp*acc_tmp(:)
     end do
     !$omp end parallel do       
       
     return
  end subroutine antuono_forcing
!! ------------------------------------------------------------------------------------------------
end module predictor_mod
