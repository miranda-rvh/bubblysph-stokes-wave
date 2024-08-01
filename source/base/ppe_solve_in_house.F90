subroutine ppe_solve_in_house
  !! This routine calculates the LHS and RHS of the PPE (or Helmholtz equation)
  !! and then calls an iterative solver to obtain the pressure.
  use kind_parameters
  use common_parameter
  use common_2d
  use omp_lib
  use linear_solver
  use calculation_gradient
  implicit none
  real(rkind) :: temp,temp_L,temp_D
  integer(ikind) :: itermax
  real(rkind) :: tol  
  integer ind,i,j,k,num,ij,n_t,nnz
  real(rkind) :: h_coh
  real(rkind),dimension(dims) :: uij,grada0,dvw,grad1
  real(rkind) :: rhs_vec_tmp,lhs_diag
  integer(ikind), dimension(:),allocatable :: ij_num,lhs_row_firstindex
  
  !! For calculating initial pressure field
  real(rkind),dimension(:,:,:),allocatable :: grad_vel0  
  real(rkind),dimension(:,:),allocatable :: udotgradu
  real(rkind),dimension(:),allocatable :: zero_rhs
  real(rkind),dimension(dims) :: rij
  real(rkind) :: Gam_sch_xx,Gam_sch_yy,Gam_sch_zz,Gam_sch
  
  
  !! Construct the RHS: (a0/dt)*div u 
  allocate(rhs_vec(n_par_fw));rhs_vec=zero
  temp=zero
  !$OMP PARALLEL DO PRIVATE(k,j,uij,rhs_vec_tmp,dvw) reduction(+:temp)
  do i=1,n_par_fw
     rhs_vec_tmp = zero
     do k=1,ij_count(i)
        j=ij_link(k,i)

        dvw = matmul(kgcm(:,:,i),ij_w_G(:,k,i))  !! Divergence weights...
   
        uij(:) = up(j,:) - up(i,:)               !! ij velocity difference         
        rhs_vec_tmp=rhs_vec_tmp + dot_product(uij(:),dvw(:))  !! ij contribution to divergence
          
      
     enddo
     rhs_vec(i) = rhs_vec_tmp*a0(i)/dt   
     
     !! Summing all elements so we can subtract mean
     temp = temp + rhs_vec(i) 
  enddo
  !$OMP END PARALLEL DO
   
  !! Remove DC from rhs_vec if closed domain flow
  if(i_open_domain==0)then
     temp = temp/dble(n_par_fw)
     !$OMP PARALLEL DO
     do i=1,n_par_fw
        rhs_vec(i)=rhs_vec(i)-temp     
     end do
     !$OMP END PARALLEL DO
  end if
  
  !! For free surface flows, we need to solve an initial Poisson equation to obtain a pressure field
  !! for the first step, to feed in to the Helmholtz equation. The RHS of this Poisson equation is zero,
  !! but that does not mean the discrete vector rhs_vec will be zero, as some of the boundary conditions
  !! on the LHS manifest in non-zero elements on the rhs_vec. Hence, we take a copy of the rhs_vec before
  !! adding boundary conditions, which is then subtracted later to obtain the zero RHS.
  if(itime==1)then 
     allocate(zero_rhs(n_par_fw));zero_rhs(:) = rhs_vec(:)
  end if

  !! Determine the first index of every row in the LHS, and also the total size of the LHS
  !! This is the only part which is done in serial, but it is FAST, and by doing it first,
  !! we can do ALL the rest in parallel.  
  allocate(lhs_row_firstindex(n_par_fw+1))
  lhs_row_firstindex(1)=n_par_fw+2
  num=n_par_fw+1
  do i=1,n_par_fw
     lhs_row_firstindex(i+1) = lhs_row_firstindex(i) + ij_count(i)
  end do
  nnz = lhs_row_firstindex(n_par_fw+1)
  
  !! Size of ppe and allocation...(only when required...)
  !! Determine whether we need to re-allocate lhs_mat and lhs_index, as allocation of
  !! these large arrays is expensive. Only do it when they grow.
  !! nnz is initialised to zero in the input routine, so lhs is automatically allocated
  !! at the first time step.
  if(nnz>nnz_old) then !! size of ppe has increased, reallocate lhs_mat,lhs_index
     if(allocated(lhs_mat))deallocate(lhs_mat,lhs_index)
     allocate(lhs_mat(nnz),lhs_index(nnz))
     nnz_old = nnz
  end if
  lhs_index=0;lhs_mat=zero
  
  !! ij_num is a temporary array which provides a reverse link to help building boundaries efficiently.
  allocate(ij_num(n_par_fw));ij_num=0

  !! This is the main loop used to build the LHS
  !$omp parallel do private(num,k,j,ij_num,temp_D,temp_L,lhs_diag,ij,n_t, &
  !$omp grada0,Gam_sch_xx,Gam_sch_yy,Gam_sch_zz,rij,Gam_sch,grad1)
  do i=1,n_par_fw   ! for all particles 
     
     !! lhs_index(1:n_par_fw) contains the index of the first element in each row.
     lhs_index(i) = lhs_row_firstindex(i)
     !! First embedded loop builds lhs_index so we can track where things go...
     !! We also evaluate the gradient of the volume fraction, for later.
     num = lhs_index(i) - 1
     grada0=zero
#ifdef schwaiger
     Gam_sch_xx=zero;Gam_sch_yy=zero;Gam_sch_zz=zero
     grad1=zero
#endif          
     do k=1,ij_count(i)   ! for all neighbours
        j=ij_link(k,i)     ! now considering particles i and j
        if(j<=n_par_fw)then
           num=num+1          
           ij_num(j)=num
           lhs_index(num)=j
        endif
        if(j>n_par_fw)then
           num=num+1
           ij = irelation(j)
           ij_num(ij)=num
           lhs_index(num) = ij
        end if
#ifdef bubbles        
        !! Build the gradient of log(a0)
        grada0 = grada0 + (log(a0(j))-log(a0(i)))*ij_w_G(:,k,i)
#endif    

#ifdef schwaiger
        !! Construct Gamma for the Schwaiger operator  
        rij = (rp(i,:)-rp(j,:))
        Gam_sch_xx = Gam_sch_xx - ij_w_L(k,i)*rij(1)*rij(1)
        Gam_sch_yy = Gam_sch_yy - ij_w_L(k,i)*rij(2)*rij(2)
#ifdef dim3
        Gam_sch_zz = Gam_sch_zz - ij_w_L(k,i)*rij(3)*rij(3)                
#endif        
#endif    

#ifdef schwaiger
       !! Construct the kernel gradient sum and store in grad1
       grad1 = grad1 + ij_w_G(:,k,i)
#endif       
    
     end do
     grada0 = matmul(kgcm(:,:,i),grada0)  !! Kernel-corrected (log of liquid fraction) gradient..

#ifdef schwaiger
     !! Finalise Gamma, and store tr(Gamma)^{-1}/n.
     !! N.B. ij_w_L contains a coefficient of two, which is cancelled here by having two in the numerator.
     Gam_sch_xx = two/(Gam_sch_xx+vsmall)
     Gam_sch_yy = two/(Gam_sch_yy+vsmall)
#ifdef dim3
     Gam_sch_zz = two/(Gam_sch_zz+vsmall)          
     Gam_sch = (Gam_sch_xx + Gam_sch_yy + Gam_sch_zz)/three
#else   
     Gam_sch = (Gam_sch_xx + Gam_sch_yy)/two
#endif     
#endif   

     !! The second embedded loop is used to populate lhs_mat
     lhs_diag=zero
     do k=1,ij_count(i)
        j=ij_link(k,i)    
        temp_L = -ij_w_L(k,i)  !! i-j interaction for Laplacian (Morris)
        temp_D = -dot_product(grada0,ij_w_G(:,k,i))  !! i-j interaction for (grad(ln(a0))).grad(a0*p)
#ifdef schwaiger
        temp_L = -Gam_sch*ij_w_L(k,i) !! Scale Morris operator by Gam_sch
        !! Add Gam_sch-scaled gradient part of Schwaiger operator
        temp_D = temp_D - Gam_sch*two*dot_product(grad1,matmul(kgcm(:,:,i),ij_w_G(:,k,i)))
#endif        
        if(i/=j)then
           lhs_diag=lhs_diag - temp_L - temp_D   ! Contribution to diagonal of lhs_mat
           if(j>n_par_fw)then    ! If j is mirror, then
              ij=irelation(j)    ! look for parent of j, ij
              if(i/=ij)then
                 n_t = ij_num(ij)
                 lhs_mat(n_t)=lhs_mat(n_t)+temp_L+temp_D   ! Contribution to lhs_mat(i,ij)
              else
                 lhs_diag=lhs_diag+temp_L+temp_D       ! Contribution to lhs_mat(i,ij) if ij=i 
              end if
              ! Building Neumann pressure boundary condition: multiply the pressure increment by a0...
              rhs_vec(i)=rhs_vec(i)-a0(i)*dP_mp(j)*(temp_L+temp_D)
           elseif(j<=n_par_fw)then   ! if j is not a mirror
              n_t = ij_num(j)   
              lhs_mat(n_t)=lhs_mat(n_t)+temp_L+temp_D    ! Contribution to lhs_mat(i,j)
           end if
        end if
     end do
     !! Build homogeneous (a0*P=0) Dirichlet conditions. Setting the Diagonal=HUGE and the RHS element=non-HUGE
     !! can only be satisfied by the solution=0.
     if(n_surf(i)/=0)then
        lhs_diag=-1.0d30 
     end if
     
     !! Add the diagonal contribution to the main matrix array
     lhs_mat(i)=lhs_mat(i)+lhs_diag
    
  end do
  !$omp end parallel do
  !! Final element of lhs_index  
  lhs_index(n_par_fw+1)=lhs_row_firstindex(n_par_fw+1)
  deallocate(ij_num,lhs_row_firstindex)  

  !! ==========================================================================
  !! This section only happens at the first time step
  !! ----------------------------------------------------------------  
  !! Poisson equation for the initial pressure field...  
  !! grad^{2}P = -div.(u.grad(u))
  if(itime==1)then
  
     !! Remove the divu* /dt component from zero_rhs (leaving only the BCs)
     !$omp parallel do
     do i=1,n_par_fw
        zero_rhs(i) = rhs_vec(i)-zero_rhs(i)
     end do
  
     !! Evaluate initial velocity gradient
     allocate(grad_vel0(dims,dims,n_par_fw));grad_vel0=zero
     call grad_operator(grad_vel0(1,:,:),upo(:,1))
     call grad_operator(grad_vel0(2,:,:),upo(:,2))
#if dim3
     call grad_operator(grad_vel0(3,:,:),upo(:,3))
#endif      
    
     !! Evaluate u.grad(u)
     allocate(udotgradu(dims,n_par_fwm))
     !$omp parallel do
     do i=1,n_par_fw
        udotgradu(1,i) = dot_product(upo(i,:),grad_vel0(1,:,i))
        udotgradu(2,i) = dot_product(upo(i,:),grad_vel0(2,:,i))
#if dim3        
        udotgradu(3,i) = dot_product(upo(i,:),grad_vel0(3,:,i))                
#endif        
     end do
     !$omp end parallel do
         
     !! Copy to mirrors
     !$omp parallel do private(i)
     do j=n_par_fw+1,n_par_fwm
        i=irelation(j);udotgradu(:,j)=udotgradu(:,i)
     end do
     !$omp end parallel do
     
     !! Add Div. (u.grad(u)) to zero_rhs
     !$OMP PARALLEL DO PRIVATE(k,j,uij,dvw)
     do i=1,n_par_fw
        do k=1,ij_count(i)
           j=ij_link(k,i)

           dvw = matmul(kgcm(:,:,i),ij_w_G(:,k,i))  !! Divergence weights...
      
           uij(:) = udotgradu(:,i) - udotgradu(:,j) !! Difference in udotgradu (-ve because want -div.(u.grad(u)))
           zero_rhs(i) = zero_rhs(i) + dot_product(uij(:),dvw(:))  !! ij contribution to divergence     
        enddo
     enddo
     !$OMP END PARALLEL DO     
       
     tol = 1e-6
     itermax = 1000  
     do i=1,n_par_fw
        p(i)=zero
     end do
     
     if(.true..and.i_open_domain==1)then   
        !$OMP PARALLEL DO PRIVATE(j)
        do i=1,n_par_fw
           !! Smooth the PPE
           if(n_surf(i)==0.and.ppe_smooth(i)/=one)then !! internal but near fs  
              zero_rhs(i) = zero_rhs(i)*ppe_smooth(i)                
              do j=lhs_index(i),lhs_index(i+1)-1
                 lhs_mat(j) = lhs_mat(j)*ppe_smooth(i)
              end do
           end if
        end do
        !$OMP END PARALLEL DO     
     end if
     !! Solve the PPE. (including scaling to reduce iterations)
     temp = 1000.0d0
     zero_rhs = zero_rhs*temp     
     call SPARSE_PBICGSTAB(lhs_mat,lhs_index,P,zero_rhs,n_par_fw,nnz,itermax,tol,LS_iters,LS_residual)
     P = P/temp
     write(6,*) "Initial pressure",minval(p(1:n_par_fw)),maxval(p(1:n_par_fw)),LS_iters,LS_residual    
     deallocate(zero_rhs,udotgradu,grad_vel0)
  end if
  !! --------------------------------------------------------------------------
  !! ==========================================================================
  
  !! Adding isothermal compressibility: the term (1/roL*c*c*dt)d(al p)/dt 
#ifdef isotcomp
  !$OMP PARALLEL DO PRIVATE(h_coh)
  do i=1,n_par_fw
     h_coh = (Ma/dt)**two
     lhs_mat(i) = lhs_mat(i) - h_coh*a0(i)
     rhs_vec(i) = rhs_vec(i) - h_coh*a00(i)*p(i)
  end do
  !$OMP END PARALLEL DO
#endif
  
  !! PPE SMOOTHING for low viscosity free surface flows
  if(.true..and.i_open_domain==1)then
     !$omp parallel do private(j)
     do i=1,n_par_fw
        if(n_surf(i)==0.and.ppe_smooth(i)/=one)then !! internal but near fs
           rhs_vec(i) = rhs_vec(i)*ppe_smooth(i) 
           do j=lhs_index(i),lhs_index(i+1)-1
              lhs_mat(j) = lhs_mat(j)*ppe_smooth(i)
           end do
        end if
     end do
     !$omp end parallel do
     deallocate(ppe_smooth)
  end if
  
  !! Set the tolerance and maximum iterations we want  
#if single  
  tol = 1e-6_rkind
#else
  tol = 1e-8_rkind
#endif    
  itermax = 500
  
  !! Zero the solution vector
  p(1:n_par_fw) = zero
  
  !! Solve the linear system. For closed domains there is no Dirichlet BC, and we need a different solver
  if(i_open_domain==0)then
     call SPARSE_PBICGSTAB_nonull(lhs_mat,lhs_index,p,rhs_vec,n_par_fw,nnz,itermax,tol,LS_iters,LS_residual)
  else
  
     !! Scale the RHS of the PPE depending on the max velocity, then reverse scale the resulting pressure.
     !! This reduces the required iterations for very low velocity flows.
     temp = one + one/(umax + 0.001d0)**two
     rhs_vec = rhs_vec*temp
     call SPARSE_PBICGSTAB(lhs_mat,lhs_index,P,rhs_vec,n_par_fw,nnz,itermax,tol,LS_iters,LS_residual)
     P = P/temp
  end if

  !! Check whether the solver failed.     
  if(isnan(LS_residual))then
     write(6,*) itime,"was NaN. Stopping."
     stop
  end if

  !! Free up space     
  deallocate(rhs_vec)
   
  !! Solver returns p(i) = pressure*a0, so need to divide by a0 to get pressure
  !$OMP PARALLEL DO
  do i=1,n_par_fw
     p(i)=p(i)/a0(i)
  end do
  !$OMP END PARALLEL DO  
    
end subroutine ppe_solve_in_house
