subroutine kernel_calculation
  use kind_parameters
  use common_parameter
  use common_2d
  use sphtools
  implicit none

  integer i,j,k

  real(rkind) rad,qq,temp,temp2
  real(rkind) rr2_tmp,kgcoeff
  real(rkind),dimension(dims) :: rij,gradw
  real(rkind),dimension(dims,dims) :: kgcf_tmp,kgcm_tmp
  
  !! Hard-coded compiler flag to specify whether to correct kernel gradients
  !! KGC==0 is plain SPH
  !! KGC==1 is Bonet & Lok correction tensor
#define KGC 1

  !! Calculate the kernel correction - Bonet and Lok tensor..
  !! Calculate the weights for gradients and laplacians..
  
  !! Allocation of memory is expensive, so only re-allocate if n_par_fw or maxneighbours has increased
  maxneighbours = maxval(ij_count(1:n_par_fw)) 
  if(n_par_fw>n_par_fw_old)then
     deallocate(kgcm)
     allocate(kgcm(dims,dims,n_par_fw));kgcm = zero !! kernel gradient correction matrix
  end if  
  if(maxneighbours>maxneighbours_old.or.n_par_fw>n_par_fw_old) then
     deallocate(ij_w_L,ij_w_G)
     allocate(ij_w_L(maxneighbours,n_par_fw));ij_w_L=zero
     allocate(ij_w_G(dims,maxneighbours,n_par_fw));ij_w_G=zero
     maxneighbours_old = maxneighbours
     n_par_fw_old = n_par_fw
  end if
  

  !$OMP PARALLEL DO PRIVATE(k,j,rij,rr2_tmp,rad,qq,gradw,temp,kgcf_tmp,temp2,kgcoeff,kgcm_tmp) &
  !$OMP SHARED(ij_w_G,ij_w_L)
  do i=1,n_par_fw !loop over fluid particles
     kgcf_tmp(:,:)=zero
     !! Dilation factor as h and vol increase...
     kgcoeff = dv*h0/(h(i)*vol(i))
     do k=1,ij_count(i)
        j=ij_link(k,i)  ! Particles j within interaction distance of i
        rij(:) = rp(i,:) - rp(j,:)     
       
        rr2_tmp = max(dot_product(rij,rij),vsmall)  ! /0 mollification here, then max() is called fewer times
       
        rad=sqrt(rr2_tmp);
        qq = rad/h(i)
        gradw(:) = kgcoeff*rij(:)*fac(qq)/rad   !uncorrected kernel gradient
        ij_w_G(:,k,i) = gradw(:)*vol(j)
        ij_w_L(k,i) = two*dot_product(rij,ij_w_G(:,k,i))/rr2_tmp                
#if KGC==1
        kgcf_tmp(1,:) = kgcf_tmp(1,:) - vol(j)*gradw(:)*rij(1)
        kgcf_tmp(2,:) = kgcf_tmp(2,:) - vol(j)*gradw(:)*rij(2)  
#if dim3
        kgcf_tmp(3,:) = kgcf_tmp(3,:) - vol(j)*gradw(:)*rij(3)  
#endif        
#endif        
     end do
#if KGC==1
#if dim3
!     kgcm(:,:,i) = matinv3(kgcf_tmp(:,:))
     call M33INV(kgcf_tmp,kgcm_tmp)
!write(6,*) "me",i,kgcm(1,1,i),kgcm(1,2,i),kgcm(1,3,i)
!write(6,*) "nasa",i,kgcm_tmp(1,1),kgcm_tmp(1,2),kgcm_tmp(1,3)
     kgcm(1,1,i) = kgcm_tmp(1,1)
     kgcm(1,2,i) = kgcm_tmp(1,2)
     kgcm(1,3,i) = kgcm_tmp(1,3)
     kgcm(2,1,i) = kgcm_tmp(2,1)
     kgcm(2,2,i) = kgcm_tmp(2,2)
     kgcm(2,3,i) = kgcm_tmp(2,3)
     kgcm(3,1,i) = kgcm_tmp(3,1)
     kgcm(3,2,i) = kgcm_tmp(3,2)
     kgcm(3,3,i) = kgcm_tmp(3,3)        
     

     
#else     
     !! And the inverse for the 2D case
     temp = one/((kgcf_tmp(1,1)*kgcf_tmp(2,2)-kgcf_tmp(1,2)*kgcf_tmp(2,1))+1.0d-15)
     temp2 = kgcf_tmp(1,1)
     kgcm(1,1,i) = temp*kgcf_tmp(2,2)
     kgcm(1,2,i) = -temp*kgcf_tmp(1,2)
     kgcm(2,1,i) = -temp*kgcf_tmp(2,1)
     kgcm(2,2,i) = temp*temp2

#endif     
#else
     !! Not correcting, set correction matrices to identity matrix
#if dim3
     kgcm(:,:,i)=zero
     kgcm(1,1,i)=one;kgcm(2,2,i)=one;kgcm(3,3,i)=one
#else
     kgcm(1,1,i)=one;kgcm(2,2,i)=one;kgcm(1,2,i)=zero;kgcm(2,1,i)=zero
#endif     
#endif
  end do
  !$OMP END PARALLEL DO
  

end subroutine kernel_calculation

