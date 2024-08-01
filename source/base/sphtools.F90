module sphtools
  !! This module contains:
  !!   > SPH kernels 
  !!   > 3x3 matrix inversion
  !!   > specialised random number generators
  use kind_parameters
  use common_parameter
  use common_2d
  implicit none

contains
#if kernel==1  
!! QUINTIC SPLINE
!! ------------------------------------------------------------------------------------------------  
  function fac(qq) result(factemp)
    real(rkind), intent(in) :: qq
    real(rkind) ::  factemp
    integer(ikind) :: iselect
    real(rkind) ::  q3, q2, q1, q34, q24, q14

!    factemp = 0.0d0
    !!  fac1=factemp/rad       !kernel gradient along i<->j direction
    iselect = floor(qq)
    if(iselect==0)then
       q3=three-qq
       q2=two-qq
       q1=one-qq
       q34=q3**four
       q24=q2**four
       q14=q1**four
       factemp=ad_7h*(-five*q34+30.0*q24-75.0*q14) 
    elseif(iselect==1)then
       q3=three-qq
       q2=two-qq
       q34=q3**four
       q24=q2**four
       factemp=ad_7h*(-five*q34+30.0*q24)
    elseif(iselect==2)then
       q3=three-qq
       q34=q3**four
       factemp=ad_7h*(-five*q34)    
    else
       factemp=zero
    endif
  end function fac
!! ------------------------------------------------------------------------------------------------  
  function Wab(qq) result(fval)
    real(rkind) ::  qq
    real(rkind) :: fval
    real(rkind) ::  q3, q2, q1, q35, q25, q15

    q3=three-qq
    q2=two-qq
    q1=one-qq
    q35=q3**five
    q25=q2**five
    q15=q1**five
    if(qq>=zero .and. qq<one)then
       fval=ad_7*(q35-6.0_rkind*q25+15_rkind*q15)
    elseif(qq>=one .and. qq<two)then
       fval=ad_7*(q35-6.0_rkind*q25)
    elseif(qq>=two .and. qq<three)then
       fval=ad_7*q35
    else!if(qq>=three)then
       fval=zero
    endif
  end function Wab
!! ------------------------------------------------------------------------------------------------  
#elif kernel==2
!! WENDLAND C2 KERNEL
!! ------------------------------------------------------------------------------------------------   
  function wab(qq) result(factemp) 
    real(rkind), intent(in) :: qq
    real(rkind) ::  factemp
    real(rkind) :: qq4
    qq4 = one - half*qq;qq4=qq4**four
    if(qq<=two) then
       factemp = ad_7*(two*qq + one)*qq4
    else
       factemp = zero    
    end if
  end function wab 
!! ------------------------------------------------------------------------------------------------                
  function fac(qq) result(factemp)      
    real(rkind), intent(in) :: qq
    real(rkind) ::  factemp
    real(rkind) :: qq4,qq3,qq2
    qq2 = qq*qq;qq3=qq*qq2;qq4=qq*qq3


    if(qq<=two) then
       factemp = ad_7h*(0.625_rkind*qq4 - 3.75_rkind*qq3 + 7.5_rkind*qq2 - five*qq) 
    else
       factemp = zero    
    end if
  end function fac
!! ------------------------------------------------------------------------------------------------
#endif
!! ------------------------------------------------------------------------------------------------
  function matinv3(A) result(B)
    !! Performs a direct calculation of the inverse of a 3×3 matrix.
    real(rkind), intent(in) :: A(3,3)   !! Matrix
    real(rkind)             :: B(3,3)   !! Inverse matrix
    real(rkind)             :: detinv

    ! Calculate the inverse determinant of the matrix
    detinv = A(1,1)*A(2,2)*A(3,3)  &
           - A(1,1)*A(2,3)*A(3,2)  &
           - A(1,2)*A(2,1)*A(3,3)  &
           + A(1,2)*A(2,3)*A(3,1)  &
           + A(1,3)*A(2,1)*A(3,2)  &
           - A(1,3)*A(2,2)*A(3,1)
    if(abs(detinv).le.vsmall) then
       detinv = one/(detinv+vsmall)
    else
       detinv = one/detinv
    end if
                   

    ! Calculate the inverse of the matrix
    B(1,1) = +detinv * (A(2,2)*A(3,3) - A(2,3)*A(3,2))
    B(2,1) = -detinv * (A(2,1)*A(3,3) - A(2,3)*A(3,1))
    B(3,1) = +detinv * (A(2,1)*A(3,2) - A(2,2)*A(3,1))
    B(1,2) = -detinv * (A(1,2)*A(3,3) - A(1,3)*A(3,2))
    B(2,2) = +detinv * (A(1,1)*A(3,3) - A(1,3)*A(3,1))
    B(3,2) = -detinv * (A(1,1)*A(3,2) - A(1,2)*A(3,1))
    B(1,3) = +detinv * (A(1,2)*A(2,3) - A(1,3)*A(2,2))
    B(2,3) = -detinv * (A(1,1)*A(2,3) - A(1,3)*A(2,1))
    B(3,3) = +detinv * (A(1,1)*A(2,2) - A(1,2)*A(2,1))
        
    
  end function
!! ------------------------------------------------------------------------------------------------ 
  SUBROUTINE M33INV (A, AINV)
  IMPLICIT NONE

      real(rkind), DIMENSION(3,3), INTENT(IN)  :: A
      real(rkind), DIMENSION(3,3), INTENT(OUT) :: AINV
      LOGICAL :: OK_FLAG

      real(rkind), PARAMETER :: EPS = 1.0D-10
      real(rkind) :: DET
      real(rkind), DIMENSION(3,3) :: COFACTOR


      DET =   A(1,1)*A(2,2)*A(3,3)  &
            - A(1,1)*A(2,3)*A(3,2)  &
            - A(1,2)*A(2,1)*A(3,3)  &
            + A(1,2)*A(2,3)*A(3,1)  &
            + A(1,3)*A(2,1)*A(3,2)  &
            - A(1,3)*A(2,2)*A(3,1)

      IF (ABS(DET) .LE. EPS) THEN
         AINV = 0.0D0
         OK_FLAG = .FALSE.
         RETURN
      END IF

      COFACTOR(1,1) = +(A(2,2)*A(3,3)-A(2,3)*A(3,2))
      COFACTOR(1,2) = -(A(2,1)*A(3,3)-A(2,3)*A(3,1))
      COFACTOR(1,3) = +(A(2,1)*A(3,2)-A(2,2)*A(3,1))
      COFACTOR(2,1) = -(A(1,2)*A(3,3)-A(1,3)*A(3,2))
      COFACTOR(2,2) = +(A(1,1)*A(3,3)-A(1,3)*A(3,1))
      COFACTOR(2,3) = -(A(1,1)*A(3,2)-A(1,2)*A(3,1))
      COFACTOR(3,1) = +(A(1,2)*A(2,3)-A(1,3)*A(2,2))
      COFACTOR(3,2) = -(A(1,1)*A(2,3)-A(1,3)*A(2,1))
      COFACTOR(3,3) = +(A(1,1)*A(2,2)-A(1,2)*A(2,1))

      AINV = TRANSPOSE(COFACTOR) / DET

      OK_FLAG = .TRUE.

      RETURN

      END SUBROUTINE M33INV


  
!! ------------------------------------------------------------------------------------------------
  function deane_and_stokes_spectrum_rand(ymin,ymax) result(y)
     !! Generate a random number with Deane & Stokes distribution in range [ymin,ymax]
     !! (slope -3/2 for y<1, -10/3 for y>1)
     real(rkind),intent(in) :: ymin,ymax  !! Args
     real(rkind) :: y    !! Output
     real(rkind) :: a,b,f0,X0,x !! Internals
     
     ! Exponents for y<=1 and y>1
     a=-three/two;b=-1.0d1/three
     
     !! Normalisation for p.d.f.
     f0 = (one/(a+one))*(one-ymin**(a+one)) + (one/(b+one))*(ymax**(b+one) - one)

     !! c.d.f at y=1
     X0 = (one/(a+one)/f0)*(one-ymin**(a+one))
     
     !! Generate uniform distributed random number
     x = rand()
     
     !! Evaluate y from inverse c.d.f of x
     if(x<=X0) then
        y = (ymin**(a+one) + f0*(a+one)*x)**(one/(a+one))
     else
        y = (f0*(b+one)*(x-X0) + one)**(one/(b+one))
     end if
   
  end function deane_and_stokes_spectrum_rand
!! ------------------------------------------------------------------------------------------------  
  function linear_spectrum_rand(ymin,ymax) result(y)
     !! Generate a random number with linear distribution within the range [ymin,ymax]  
     real(rkind),intent(in) :: ymin,ymax  !! Args
     real(rkind) :: y    !! Output
     real(rkind) :: x,min4,diff4
     real(rkind) :: a,b,c
     
     !! Slope
     a=0.2d0
     
     !! Integral from xmin to xmax (for normalisation)
     c= half*a*ymax*ymax - half*a*ymin*ymin
                     
     !! Generate uniform distributed random number
     x = rand()
          
     !! Evaluate y from inverse c.d.f. of x
     y= sqrt((x + half*a*ymin*ymin)*two/C/a)
   
  end function linear_spectrum_rand  
!! ------------------------------------------------------------------------------------------------  
  function cubed_spectrum_rand(ymin,ymax) result(y)
     !! Generate a random number with cubic distribion within the range [ymin,ymax]
     real(rkind),intent(in) :: ymin,ymax  !! Args
     real(rkind) :: y    !! Output
     real(rkind) :: x,min4,diff4,xtest,cdf
     integer(ikind) :: i
    
     !! powers of 4...
     min4 = ymin - quarter*ymin**four
     diff4 = ymax - quarter*ymax**four - min4

     !! Uniformly distributed random number
     xtest = rand()

     !! Numerically invert the cdf
     cdf = zero     
     i=0
     do while (xtest>=cdf.and.i<=101)
        i = i+1
        y = ymin + (ymax-ymin)*dble(i-1)/100.0        
        cdf = (y - quarter*y**four - min4)/diff4
     end do
     
   
  end function cubed_spectrum_rand    
!! ------------------------------------------------------------------------------------------------  
  function flat_spectrum_rand(ymin,ymax) result(y)
     !! Generate a random number uniformly distributed within the range [ymin,ymax]
     real(rkind),intent(in) :: ymin,ymax  !! Args
     real(rkind) :: y    !! Output
    
     !! Generate uniform distributed random number
     y = rand()

     !! Scale to range
     y = ymin + (ymax - ymin)*y
   
  end function flat_spectrum_rand    
!! ------------------------------------------------------------------------------------------------  
end module sphtools
