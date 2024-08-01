module freesurf
  use kind_parameters
  use common_2d
  use common_parameter
  use sphtools
  use omp_lib
  implicit none

  real(rkind) :: beta_div_r
  real(rkind),dimension(dims) :: rij
  real(rkind),dimension(:),allocatable :: div_r
  private
  public :: locate_free_surf_part

  !! Module contains routines to locate the free surface

contains
!! ------------------------------------------------------------------------------------------------
  subroutine locate_free_surf_part 
    !! Calculate some things about the particle distribution, then mark free surface particles
    integer(ikind) :: i
   
    allocate(div_r(n_par_fw));div_r=zero

    !! find div.r, conc and surf_norm, grad conc 
    call calc_part_dist_vars
 
    !! Free surface threshold (dimension dependent)
#if dim3
    beta_div_r=2.3_rkind   !! 2.5 for wendland
#else 
    beta_div_r=1.5_rkind
#endif    

    !! Indentify free surface particles for open and closed domains
    if(i_open_domain==1)then
       !$OMP PARALLEL DO 
       do i=1,n_par_fw
          ! Free surface identity criterion
          if(div_r(i)>beta_div_r) then 
             n_surf(i)=0
          else
             n_surf(i)=1
          end if
       end do
       !$OMP END PARALLEL DO
    else  !! Closed domain,
       n_surf(:)=0             
    end if

    !! For mirrors, if parent is FS, mirror is FS
    !$OMP PARALLEL DO
    do i=n_par_fw+1,n_par_fwm  
       n_surf(i) = n_surf(irelation(i))
       surf_norm(i,:) = surf_norm(irelation(i),:)
    end do
    !$OMP END PARALLEL DO

    !! smooth the surface normals
    if(i_open_domain==1)then
       call smooth_normals
    end if

    !! A variable to smooth the PPE near free surfaces
    if(i_open_domain==1)then
       if(.not.allocated(ppe_smooth)) allocate(ppe_smooth(n_par_fw))
       ppe_smooth(:) = one
       do i =1,n_par_fw
          if(div_r(i)>beta_div_r+0.32_rkind) cycle
          if(n_surf(i)==1) then
             ppe_smooth(i) = zero
          else
             !! Option 1 works with Wendland, option 2 with Quintic??
             ppe_smooth(i) = half*(one - cos((div_r(i)-beta_div_r)*pi/(0.32_rkind)))
!             ppe_smooth(i) = half*(one - cos((div_r(i)-beta_div_r)*pi/(0.2d0))) !*beta_div_r
           end if
       end do  
    end if  

    deallocate(div_r)
  end subroutine locate_free_surf_part
!! ------------------------------------------------------------------------------------------------
  subroutine calc_part_dist_vars
    !! Calculate variables related to particle distribution: div.r, surf_norm, conc, grad_conc
    !! used for free surface identification and shifting
    integer(ikind) :: i,j,k
    real(rkind) :: qq,rad,temp,wtmp
    real(rkind) :: conc_tmp,div_r_tmp,kcoeff
    real(rkind),dimension(dims) :: surf_norm_tmp,sum_ft_tmp
    real(rkind) :: nt,Rt,Wab1,fab,conc_sum
   
    !! initial allocation
!    if(allocated(surf_norm))deallocate(surf_norm)    
    if(.not.allocated(conc)) allocate(conc(n_par_fwm))
    if(.not.allocated(surf_norm)) allocate(surf_norm(n_par_fwm,dims))
    allocate(grad_conc(n_par_fw,dims))

    !! Parameters for the tensile instability correction
    nt=four
    Rt=quarter
    qq=dx/h0  ! JRCK put this in (so h/dx can vary) (previously qq=1.3 was hard-coded)
    Wab1=Rt*Wab(qq)**(-nt)

    !! calculate div.r, normals and concentration, concentration gradient 
    !$OMP PARALLEL DO PRIVATE(conc_tmp,div_r_tmp,surf_norm_tmp,k,j,rij,rad,qq,temp,wtmp,sum_ft_tmp,kcoeff)
    do i=1,n_par_fw
       conc_tmp =zero;div_r_tmp=zero;surf_norm_tmp=zero
       sum_ft_tmp =zero
       kcoeff = dv/vol(i)
       do k=1,ij_count(i)
          j=ij_link(k,i) 
          rij(:) = rp(i,:) - rp(j,:)
          rad = sqrt(dot_product(rij,rij));qq = rad/h(i)
          wtmp = Wab(qq)*kcoeff
          
          !! Concentration
          conc_tmp = conc_tmp + wtmp*vol(j)


          !! Grad-C interaction
          fab=one + wtmp*wtmp*wtmp*wtmp*wab1     

          div_r_tmp=div_r_tmp - dot_product(rij,ij_w_G(:,k,i))
          surf_norm_tmp(:)=surf_norm_tmp(:) - ij_w_G(:,k,i)

          sum_ft_tmp(:)=sum_ft_tmp(:)+ij_w_G(:,k,i)*fab

       end do
       conc(i) = conc_tmp
       div_r(i)=div_r_tmp
       surf_norm(i,:) = surf_norm_tmp(:)*h(i) ! normalise so |n|=0.504989 for planar surface cartesian distribution for any h
       grad_conc(i,:) = sum_ft_tmp(:)
    end do
    !$OMP END PARALLEL DO
    conc(1:n_par_fw) = conc(1:n_par_fw)
    
    
  end subroutine calc_part_dist_vars
!! ------------------------------------------------------------------------------------------------
  subroutine smooth_normals
    !! Shepard filter the surface normal vectors
    integer(ikind) :: i,j,k
    real(rkind) :: rad,qq,kcoeff
    real(rkind),dimension(dims) :: sn_temp_tmp
    real(rkind),dimension(:,:),allocatable :: sn_temp
    allocate(sn_temp(n_par_fwm,dims))

    !$OMP PARALLEL DO PRIVATE(sn_temp_tmp,k,j,rij,rad,qq,kcoeff)
    do i=1,n_par_fw
       sn_temp_tmp=zero
       kcoeff = dv/vol(i)
       do k=1,ij_count(i)
          j=ij_link(k,i)
          if(j<=n_par_fw) then ! there are no surf_norm for j in mirror
             rij(:) = rp(i,:) - rp(j,:)
                
             rad = sqrt(dot_product(rij,rij))
             qq = rad/h(i)
             sn_temp_tmp(:) = sn_temp_tmp(:) + kcoeff*Wab(qq)*surf_norm(j,:)*vol(j)
          end if
       end do
       sn_temp(i,:) = sn_temp_tmp(:)
    end do
    !$OMP END PARALLEL DO

    !! Normalise with concentration
    !$OMP PARALLEL DO
    do i=1,n_par_fw
       surf_norm(i,:) = sn_temp(i,:)/conc(i)
    end do
    !$OMP END PARALLEL DO
    !! Copy to mirrors
    !$omp parallel do
    do i=n_par_fw+1,n_par_fwm
       surf_norm(i,:) = surf_norm(irelation(i),:)
    end do
    !$omp end parallel do
    deallocate(sn_temp)
  end subroutine smooth_normals 
!! ------------------------------------------------------------------------------------------------
end module freesurf
