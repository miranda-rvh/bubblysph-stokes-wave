subroutine put_in_bin
  !! This routine loops over all particles and identifies any which are "in trouble"
  !! which it then puts in the "bin", so they are skipped by the remainder of the code.
  use kind_parameters
  use common_parameter
  use common_2d
  integer(ikind) :: i


  ! parallel loop over all fluid particles
  !$omp parallel do
  do i=n_par_w+1,n_par_fw
     ! is the particle out of the domain?
     if(rp(i,1)<xmint.or.rp(i,1)>xmaxt.or.rp(i,2)<ymint.or.rp(i,2)>ymaxt)then
        inbin(i) = .true.
     end if
#ifdef dim3
#ifdef zwall
    if(rp(i,3)<zmint.or.rp(i,3)>zmaxt) then  !! Is it out of the domain in Z?
       inbin(i) = .true.
    end if
#endif
#endif
     ! is the particle position NaN
     if(isnan(rp(i,1)).or.isnan(rp(i,2)))then
        inbin(i) = .true.
     end if
  end do
  !$omp end parallel do

  ! now (serial) loop to correct position of everything in the bin  
  n_inbin = 0
  n_free = 0
  do i=n_par_w+1,n_par_fw
     if(inbin(i))then
        ! re-locate the particle to somewhere out of the way...
        rp(i,1) = two*xmint - xmaxt
        rp(i,2) = two*ymint - ymaxt
        ! give it a zero velocity
        up(i,:) = zero
        n_inbin = n_inbin + 1      !! add 1 to the number in the bin

        !! Put the particle in the free particles array
        n_free = n_free + 1   ! add 1 to the number of free indices 
        p_free(n_free) = i    ! store the free index in the list of free indices       
     end if
  end do
  
  !! Check whether everything is in the bin... stop if so.
  if(n_inbin==n_par_fw-n_par_w.and.n_inbin/=0)then
     write(6,*) "All ",n_par_fw - n_par_w," liquid particles are in the bin. Stopping."
     stop
  end if
  
  return
end subroutine put_in_bin
!! ------------------------------------------------------------------------------------------------
