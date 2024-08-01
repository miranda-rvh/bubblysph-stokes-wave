module bubble_acoustics
  !! This module contains routines to determine the acoustic signal at a given location
  use kind_parameters
  use common_parameter
  use common_2d

  implicit none
  real(rkind),dimension(dims) :: ra !! Position vector of microphone
  real(rkind) :: ab_min 
  real(rkind),dimension(:),allocatable :: pa,ta !! Time and pressure for acoustic signals
  integer(ikind) :: na_samp
  real(rkind) :: dta_samp
  real(rkind),parameter :: phiom = sqrt(3.0d0*1.3d0*1.0d5/1.0d3) !! from deane & Stokes
  integer(ikind) :: ia_lastout
  integer(ikind),parameter :: Nimages=100
  

contains
#ifdef bubbles
!! ------------------------------------------------------------------------------------------------
  subroutine create_acoustic_array
     !! Set up arrays to hold acoustic signal
     integer(ikind) :: i
     real(rkind) :: ommax,fmax
     
     
     !! Minimum bubble size
     ab_min = 5.0d0*Rhinze/1.0d2
         
     !! Maximum frequency:
     ommax = phiom/ab_min
     fmax = ommax/(two*pi)
     
     !! Acoustic sampling interval and number of samples
     dta_samp = one/(two*fmax) !! Satisfy Nyquist
     na_samp = floor(time_max/dta_samp) + 1
             
     !! Allocate arrays
     allocate(pa(na_samp),ta(na_samp))
     pa=zero
     ta=zero
     
     !! Initialise time data
     !$omp parallel do
     do i=1,na_samp
        ta(i) = dble(i-1)*dta_samp
     end do
     !$omp end parallel do
     
     
     !! Set the location of the microphone - hard-coded to the top-right corner for now
     ra(1) = xb_max
     ra(2) = yb_max
#ifdef dim3
     ra(3) = zero
#endif                   
     
     write(6,*) "Acoustic array sizes: fmax,dta_samp,na_samp"
     write(6,*) fmax,dta_samp,na_samp
  
     !! Initialise the output:
     ia_lastout = 0
     open(unit=301,file='data_directory/statistics/acoustic_signal.out')
  
     return
  end subroutine create_acoustic_array
!! ------------------------------------------------------------------------------------------------
  subroutine add_pressure_from_bubble(abub,rbub,tbub)
     !! This routine takes in a bubble size, position and time-of-creation, and adds the contribution
     !! of this bubble to the acoustic pressure array
     real(rkind),intent(in) :: abub,tbub
     real(rkind),dimension(dims),intent(in) :: rbub
     real(rkind),dimension(dims) :: rarbub
     real(rkind) :: r,om,T_decay,A,tend,ti,T_delay
     integer(ikind) :: i,istart,iend,j,k
     real(rkind) :: sign_flag,tmp
     
     !! Determine whether to have a positive or negative initial oscillation?
     tmp = rand()
     if(tmp.gt.half) then
        sign_flag = one
     else
        sign_flag = -one
     end if
     
     !! Evaluate angular frequency of oscillations
     om = phiom/abub

     !! Decay time-scale
     T_decay = one/(0.0023*om**(4.0/3.0))    
     
     !! Loop over all images of the bubble due to z-periodicity
#ifdef dim3
     k=-Nimages
     do j=1,2*Nimages+1
        k=k+1
#endif        
        !! Evaluate distance between bubble and microphone
        rarbub= rbub - ra
#ifdef dim3
        rarbub(3) = rarbub(3) + dble(k)*(zb_max-zb_min)
#endif        
        r = sqrt(dot_product(rarbub,rarbub))
        r=max(r,vsmall)
        
        !! Acoustic travel time between bubble and microphone
        T_delay = r*Maa
 
        !! Amplitude (add the scaling with distance here: more efficient)
!        A = abub/r
        A = exp(-0.1*(abub/Rhinze)**two)
        A = A/r
     
        !! Index immediately after tbub
        istart = floor((tbub+T_delay)/dta_samp) + 1
     
        !! End time, and last index
        tend = tbub + T_delay + 7.0d0*T_decay  !! exp(-tend/T_decay) ~= 0.001
        iend = floor(tend/dta_samp) + 1
         
        !! Loop over the interval during which the bubble produces sound
        do i=istart,iend
           ti = ta(i)-(tbub+T_delay) !! Time since tbub at index i

           !! Augment pressure
           pa(i) = pa(i) + A*exp(-ti/T_decay)*cos(om*ti)
        end do
#ifdef dim3
     end do
#endif
  
     return
  end subroutine add_pressure_from_bubble
!! ------------------------------------------------------------------------------------------------
  subroutine write_acoustic_file
     !! This routine writes the acoustic file for a given set of time-steps
     integer(ikind) :: iend,i
   
     ia_lastout = ia_lastout + 1
     
     iend = floor((time-dt)/dta_samp)
     iend = max(iend,ia_lastout)
     
     do i=ia_lastout,iend
        write(301,*) ta(i),pa(i)
     end do
     flush(301)
     
     ia_lastout = iend
     
     return
  end subroutine write_acoustic_file  
!! ------------------------------------------------------------------------------------------------
#endif
end module bubble_acoustics
