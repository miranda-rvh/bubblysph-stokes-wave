!!================================================================================================= 
!!                  BUBBLYSPH
!! Dr Jack King                                                                          
!!                                                                                               
!! Most relevant references                                                                      
!! ------------------------
!! (1) King et al. (2022) ArXiv:2206.01641 
!!     Bubbly-LES-SPH. Details of the turbulence and bubble modelling used
!!      
!! (2) King & Lind (2021) JNNFM 293:104556, available at ArXiv:2009.12245
!!     Viscoelastic SPH and a bit more SPH detail
!!
!! (3) Lind et al. (2012) JCP 231(4):1499-1523, doi: 10.1016/j.jcp.2011.10.027
!!     Incompressible SPH with shifting
!!
!!=================================================================================================
program bubblysph
  use kind_parameters
  use common_parameter
  use common_2d
  use input
  use omp_lib
  use evolve
  use output
  use bubble_acoustics

  implicit none
  real(rkind) :: ts_start,ts_end,t_run,t_per_dt,t_last_10
  integer(ikind) :: n_threads,narg
  character(len=200) :: argch

  !! Ensure we don't run non-slip walls with wavemaker
#ifdef wavemaker
#ifndef slip
  write(6,*) "Code compiled with wavemaker and no slip walls. Incompatible."
  write(6,*) "Please recompile either with slip, or without wavemaker."
  write(6,*) "Stopping."
  stop
#endif
#endif  
  
  
  !! Set the number of threads:
  narg = command_argument_count()
  if(narg/=0)then                         !! set threads from the argument (e.g. ./sph 4 runs on 4 threads)
     call get_command_argument(1,argch)
     argch = adjustl(argch)
     read(argch,938) n_threads
     938 format(i4)     
     call omp_set_num_threads(n_threads)
     write(6,*) "User specified ",n_threads,"threads"
  else                                    !! Or just use as many threads as set by environment
     !$omp parallel
     n_threads = omp_get_num_threads()
     !$omp end parallel
     write(6,*) "Automatically allocated",n_threads,"threads"
  end if

  ! Read in data and initialize the variables
  call getdata
  call setup_stats_files
  
  !! Set the total domain size to allow for mirror particles
  xmint=xb_min - sup_size*three*half
  xmaxt=xb_max + sup_size*three*half
  ymint=yb_min - sup_size*three*half
  ymaxt=yb_max + sup_size*three*half
#ifdef dim3
  zmint = zb_min - sup_size*three*half
  zmaxt = zb_max + sup_size*three*half
#endif  
  
 
  !! Initialise some counters and time-stamps
  nct = 0
  t_run = zero
  t_last_10 = zero
  next_dump = time;next_dump2 = time
  n_dump = -1;n_dump2= 0
  output_this_step = .true.
  
  !! START OF MAIN LOOP--------------------------------------------------------
  !! ==========================================================================
!  do while (itime<=itmax)
  do while (time<=time_max)

     !! Increment counter
     itime = itime+1

     ! profiling for screen output
     ts_start=omp_get_wtime()

     !! Build the cells if we need to (only at step 0, when nct=0)
     call build_cells

     !! Check whether to output data this step
     output_this_step = .false.
     if(time>next_dump.or.output_everytime)then
        next_dump = next_dump + dt_out
        n_dump = n_dump + 1
        output_this_step = .true.
     end if

     !! Carry out a single complete time step
     !! "step" is a subroutine in the evolve module, and is the heart of the code. 
     !! it's a good place to start.
     call step

     !! Output acoustics if required
#ifdef acoustics
     call write_acoustic_file
#endif   
     
     !! Increment time
     time = time + dt

     !! profiling for screen output
     ts_end=omp_get_wtime()
     t_run = t_run + ts_end - ts_start
     t_per_dt = t_run/dble(itime)
     t_last_10 = t_last_10 + ts_end - ts_start
     
     !! Output data to screen - this will look nice on a 24-line deep terminal.
     if(mod(itime,10)==0)then
        write(6,*)"itime,time,dt=", itime,time,dt
        write(6,*) "n_par_fwm,n_par_fw,nmirror",n_par_fwm,n_par_fw,nmirror
        write(6,*) "Linear solver:",LS_iters,"iterations,error:",LS_residual
        write(6,*) "Maximum velocity in domain: ",umax
        write(6,*) "There are ",n_inbin,"particles in the bin"
        !$OMP PARALLEL
        n_threads = omp_get_num_threads()
        call random_seed()
        !$OMP END PARALLEL       
        write(6,*) "Number of threads=",n_threads,"Run time=",t_run
        write(6,*) "run-time/dt=",t_per_dt,"Moving avg=",0.1*t_last_10
        t_last_10 = zero
#ifdef bubbles
        write(6,*) " "
        write(6,*) " "        
        write(6,*) "BUBBLES"
        write(6,*) "Number of bubbles in domain = ",n_bub-b_nfree
        write(6,*) " "
#else        
        write(6,'(/,/,/,/,A)') "  "
#endif        
        write(6,'(/,/,/,A)') " "
        write(6,'(/,/,/,/,/,/,A)') " "                       
     end if
     
     !! Output data to files if required
     if(output_this_step) then
        call output_fields
     end if     
     call output_statistics !! This happens EVERY step       

  


     !! Do a little de-allocation
     deallocate(conc) 
     
     !! Output some data EVERY step - only ever used if something goes wrong for debugging really...   
     write(92,*) time,dt,n_par_fw-n_inbin,maxval(abs(up(1:n_par_fw,1))) &
                 ,maxval(abs(up(1:n_par_fw,2))),maxval(abs(P(1:n_par_fw)))
     write(93,*) LS_iters,LS_residual
     flush(92);flush(93)

     

  end do
  !! END OF MAIN LOOP ---------------------------------------------------------
  !! ==========================================================================
  
  !! Output at final step
  call output_fields

  !! Close some files and deallocate some space. In practice, this isn't important as
  !! we rarely run the code until the end - usually stop it, run out of time, or
  !! it crashes!!
  close(92)  
  close(93)  
  close(94)  
  deallocate(rp,up,p,inbin)

  stop
end program bubblysph
!! ------------------------------------------------------------------------------------------------
subroutine build_cells
  use kind_parameters
  use common_parameter
  use common_2d
  use cell_griddin  
  
  !! Set the support size
  sup_size = ss*maxval(h(1:n_par_fw))
  sup_size_2 = sup_size*sup_size

  !! Build cells and cell neighbour lists, but only if this is the start
  if(nct==0) then

     !! Deallocate as required
     if(allocated(ic_count))then
        deallocate(ic_count,ic_link)
     end if
     
     !! Generate cells including space for mirror particles
     call cell_generator
#if dim3
     call cells_third_dimension
#endif      
     
     !! Cell neighbour lists
#if dim3
     call neighbour_cells3d
#else
     call neighbour_cells  
#endif

   end if

   !! Check if the smoothing length is too big...
   !! N.B. this becomes an issue if the liquid volume fraction gets too small. It could be avoided 
   !! by re-building cells, but in practice we limit the liquid volume fraction and it's not an issue.
   if(sup_size>cell_size) then
      write(6,*) "ERROR: sup_size>cell_size. Stopping."
      stop
   end if

   return
end subroutine build_cells
!! ------------------------------------------------------------------------------------------------
