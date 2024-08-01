module shifting
  !! Shift particles, and apply some little velocity corrections near walls 

  use kind_parameters
  use common_parameter
  use common_2d
  use calculation_gradient
  use sphtools
  implicit none

!! Set the type of shifting coefficient. 0 gives the 0.25h^2 type, 1 gives the \propto U type
!! see King & Lind (2021) JNNFM for details.
#define shift_type 0


  private 
  public :: new_fick_shift,adjust_velocity_walls
contains
!! ------------------------------------------------------------------------------------------------
  subroutine new_fick_shift
    integer(ikind) :: i
    real(rkind) :: temp, dCds,maxU,absU,dCds2
    real(rkind), dimension(dims) :: surf_t,surf_n,surf_t2,tmp_vec
    real(rkind) :: Diff_coeff,meanshift

    ! Find the maximum particle speed if required (only really for v. low Re flows).
#if shift_type==1    
    maxU = vsmall
    !$omp parallel do private(temp) reduction(max:maxU)
    do i=1,n_par_fw
       temp = sqrt(dot_product(up(i,:),up(i,:)))
       maxU = temp
    end do
    !$omp end parallel do
#endif    

    !! Loop over all particles
    !$OMP PARALLEL DO PRIVATE(diff_coeff,absU,temp,surf_n,surf_t,surf_t2,dCds,dCds2,tmp_vec)
    do i=1+n_par_w,n_par_fw
    
       !! Set the shifting coefficient
#if shift_type==0
#ifdef wavemaker
       diff_coeff = quarter*h(i)*h(i)
!       if(rp(i,1).le.50.0d0*dx) diff_coeff = diff_coeff*2.0d0 !! Make shifting coefficient bigger near wavemaker
#else
       diff_coeff = quarter*h(i)*h(i)
#endif       
#else
       !! Alternative shift coeff, useful for low Reynolds number flows
       absU = sqrt(dot_product(up(i,:),up(i,:)))
       temp = 0.2*exp(-absU/(0.2*maxU))  !! 0.2 for jet!!!
       diff_coeff = four*h(i)*dt*(absU + temp*maxU)  !4.0
#endif     
    
       !! Evaluate the shifting vector
       tmp_vec = zero
       ! OPTION 1: A free surface particle - remove shifting normal to surface  
       if(n_surf(i) == 1.and.i_open_domain==1) then   !! If using n_surf in closed domain, ignore for shifting...

          ! Unit surface normal and tangent
          temp=sqrt(dot_product(surf_norm(i,:),surf_norm(i,:)))
          surf_n(:)=surf_norm(i,:)/max(temp,vsmall) ! stabilise for div by 0
#if dim3          
          surf_t(1) = -surf_n(2)   !! non-parallel with n
          surf_t(2) = surf_n(1)
          surf_t(3) = surf_n(3)

          !! surf_t2 = surf_n x surf_t
          surf_t2(1) = surf_n(2)*surf_t(3) - surf_t(2)*surf_n(3)  !! orthogonal to n
          surf_t2(2) = surf_n(3)*surf_t(1) - surf_t(3)*surf_n(1)
          surf_t2(3) = surf_n(1)*surf_t(2) - surf_t(1)*surf_n(2)

          surf_t(1) = surf_n(2)*surf_t2(3) - surf_t2(2)*surf_n(3)  !! orthogonal to n and t2
          surf_t(2) = surf_n(3)*surf_t2(1) - surf_t2(3)*surf_n(1)
          surf_t(3) = surf_n(1)*surf_t2(2) - surf_t2(1)*surf_n(2)

          dCds = dot_product(surf_t(:),grad_conc(i,:))
          dCds2 = dot_product(surf_t2(:),grad_conc(i,:))
          tmp_vec(:) = -Diff_coeff*(dCds*surf_t(:) + dCds2*surf_t2(:))
#else
          surf_t(1) = -surf_n(2);surf_t(2) = surf_n(1)
          dCds = dot_product(surf_t(:),grad_conc(i,:))
          tmp_vec(:) = -Diff_coeff*dCds*surf_t(:)
#endif          
          
       ! OPTION 2: An ordinary internal particle, shift according to Fick's law
       else
          tmp_vec(:) = -Diff_coeff*grad_conc(i,:)
       endif

       ! Set limit on distance which can be shifted (per time step) < 0.1dx (Xu recommends, esp. for violent flows)
       temp=sqrt(dot_product(tmp_vec,tmp_vec))
       if(temp>0.1d0*h(i)) then
             tmp_vec(:) = 0.1d0*h(i)*tmp_vec(:)/temp
       end if
       
       !! Convert shifting vector into shifting velocity (divide by dt)
       ushift(i,:) = tmp_vec(:)/max(dt,vsmall)
    enddo
    !$OMP END PARALLEL DO
    
    
    meanshift = zero
    !$omp parallel do reduction(+:meanshift)
    do i=1,n_par_fw
       meanshift = meanshift + dot_product(ushift(i,:),ushift(i,:))
    end do
    !$omp end parallel do
    meanshift = sqrt(meanshift/dble(n_par_fw))
    
!    write(711,*) time,meanshift
!    flush(711)
    
    !! This was allocated and evaluated in free-surf module, used here, no longer needed
    deallocate(grad_conc)

    return
  end subroutine new_fick_shift
!! ------------------------------------------------------------------------------------------------    
  subroutine adjust_velocity_walls
    !! Finds all fluid particles within 0.5dx of a wall, and reverses the component of
    !! velocity normal to the wall
    integer(ikind) :: i, j, k,jj
    real(rkind) :: temp1,u_norm
    integer(ikind), dimension(:), allocatable :: imirror_close

    real(rkind) :: r_min2
    real(rkind),dimension(dims) :: uij,wall_norm

    allocate(imirror_close(n_par_fw))
    imirror_close=0_ikind

    ! Find the particles close to wall mirrors
    !$OMP PARALLEL DO PRIVATE(i,k,j,r_min2)
    do jj=n_par_fw+1,n_par_fwm !! Look through all mirrors
       i=irelation(jj)  !! and then look through their parents...
       if(i>n_par_w) then !! but only their fluid parents...
          r_min2=dx*dx ! cut-off distance is dx
          do k=1,ij_count(i)  !! Look through their neighbours        
             j=ij_link(k,i)

             !! For particles close to wall mirrors, mark them as such
             if(j>n_par_fw .and. dot_product(rp(i,:)-rp(j,:),rp(i,:)-rp(j,:))< r_min2) then
                if(vrelation(j)>0.and.vrelation(j)<999)then ! edges
                   if(b_type(vrelation(j))==1) then             ! walls
                   imirror_close(i)=j
                   end if
                end if
                if(vrelation(j)==-1.or.vrelation(j)<-2)then  ! corners
                   imirror_close(i)=j
                end if
             end if
          end do
       end if
    end do
    !$OMP END PARALLEL DO

    ! Correct the velocity
    !$OMP PARALLEL DO PRIVATE(i,j,wall_norm,temp1,uij,u_norm)
    do jj=n_par_fw+1,n_par_fwm
       i = irelation(jj)
       if(imirror_close(i)/=0)then   ! if it is close to a wall mirror
          j=imirror_close(i)

	  ! calculate the unit vector between i and the wall mirror
          wall_norm(:) = rp(i,:)-rp(j,:)
          temp1=sqrt(dot_product(wall_norm,wall_norm))
          wall_norm(:) = wall_norm(:)/temp1

          !! relative velocity
          uij(:) = up(i,:) - up(j,:)

	  ! calculate velocity component towards the mirror
          u_norm=dot_product(uij,wall_norm)  

          !! remove particle velocity component toward mirror
          if(u_norm<=zero)then
             up(i,:)=up(i,:) - u_norm*wall_norm(:)
          end if
       end if
    end do
    !$OMP END PARALLEL DO

    deallocate(imirror_close)


  end subroutine adjust_velocity_walls
!! ------------------------------------------------------------------------------------------------  
end module shifting
