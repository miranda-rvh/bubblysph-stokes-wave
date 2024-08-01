module calculation_gradient
  !! Two routines for gradient calculations
  use kind_parameters
  use common_parameter
  use common_2d
  use mirror_boundaries
  use sphtools
  use omp_lib
  implicit none

  private
  public :: calc_pressure_gradient,grad_operator
contains
!! ------------------------------------------------------------------------------------------------
  subroutine calc_pressure_gradient (grad_p)
  !! Calculate grad(p) - grav
    integer(ikind) :: i,j,k
    real(rkind) :: temp
    real(rkind),dimension(dims) :: grad_p_tmp
    real(rkind), dimension(:,:), intent(inout) :: grad_p

    ! set the pressure in mirror particles using dP_mp
    !$OMP PARALLEL DO PRIVATE(j)
    do i=n_par_fw+1,n_par_fwm
       j = irelation(i)
       P(i) = P(j) + dP_mp(i)    !! May need modifying with volume fraction
    end do
    !$OMP END PARALLEL DO
    
    grad_p=zero

    ! The form of the gradient approximation: Eqn 5, Lind et al.,JCP,2011
    !$OMP PARALLEL DO PRIVATE(k,j,temp,grad_p_tmp)
    do i=1,n_par_fw
       grad_p_tmp(:) = zero
       do k=1,ij_count(i)
          j=ij_link(k,i)
          temp = (a0(j)*P(j) - a0(i)*P(i))
          grad_p_tmp(:) = grad_p_tmp(:) + temp*ij_w_G(:,k,i)
       end do
       grad_p(i,:) = matmul(kgcm(:,:,i),grad_p_tmp(:)) - a0(i)*grav(:)  !! Including gravity here...   
    end do
    !$OMP END PARALLEL DO

  end subroutine calc_pressure_gradient
!! ------------------------------------------------------------------------------------------------
  subroutine grad_operator(grad,pha)
    !! Calculate the gradient of a field pha
    integer(ikind) :: i, j, k
    real(rkind) :: temp
    real(rkind), dimension(:), intent(in) :: pha
    real(rkind), dimension(:,:), intent(inout) ::grad
    real(rkind), dimension(dims) :: grad_tmp

    grad(:,:)=zero
    !$OMP PARALLEL DO PRIVATE(k,j,temp,grad_tmp)
    do i=1,n_par_fw
       grad_tmp(:)=zero
       do k=1,ij_count(i)
          j=ij_link(k,i)
          temp = pha(j) - pha(i)
          grad_tmp(:) = grad_tmp(:) + temp*ij_w_G(:,k,i)
       end do
       grad(:,i) = matmul(kgcm(:,:,i),grad_tmp(:))     
    end do
    !$OMP END PARALLEL DO

  end subroutine grad_operator
!! ------------------------------------------------------------------------------------------------
end module calculation_gradient
