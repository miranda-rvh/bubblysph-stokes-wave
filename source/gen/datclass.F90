program datgen
  !! This program generates the input data for sph3d: IPART, IBOUND and INDAT
  !! It reads in a case number, then generates the data.
  !! Details of each case are HARD-CODED.
  !! More cases can be added, just copy an existing one as a template.

  use kind_parameters
  use common_parameter
  use global_variables 
  implicit none

  real(rkind) :: x,y
  integer(ikind) :: nz
  real(rkind) :: dt_coef_advection,dt_coef_viscous,dt_coef_acoustic
  real(rkind) :: wave_slope,xmax_local

  integer ipart,itest,shift_coeff
  integer i,j,ires
  double precision h0,eta_fs,yl
  real(rkind) :: k_,h_, g_,ak_,a_,alpa,sgma,b_

  write(*,*) 'Cases: '
  write(*,*) '  case  1:  Still water/bubble column'
  write(*,*) '  case  2:  Lobovsky dam break'
  write(*,*) '  case  3:  Taylor-Green'
  write(*,*) '  case  4:  Poiseuille Flow/Channel flow'
  write(*,*) '  case  5:  3rd order Stokes Wave'
  write(*,*) '  case  6:  Focussed waves with a wavemaker'
  write(*,*) '  '
  write(*,*) 'Input test case number: '
  
  itest=5
  write(6,*) "HARDCODED TO TEST CASE 5"
!  read(*,*) itest

  select case (itest) 
!! ------------------------------------------------------------------------------------------------
!!  Below are several cases, each goes through the same setup process. Details are HARD CODED, 
!!  so you need to recompile each time you change a case.
!! ------------------------------------------------------------------------------------------------     
  case(1)
    !! Still water (also used for bubble column)


     !! Set domain size, resolution and flow type
     i_open_domain=1               !! 1 indicates a free surface, 0 indicates a closed domain. 
     slip_walls=.false.            !! Slip or no-slip walls
     h0=1.0d0;xl=1.0d0;yl=1.3d0    !! A lengthscale, and X and Y lengthscales for domain. Must be O(1)
     dx=xl/100.0;dy=dx             !! Particle spacing relative to lengthscale
     nz = 10                       !! Number of copies in 3rd dimension (if 3D)
     h=1.3d0*dx                    !! Initial smoothing length relative to particle spacing
     
     !! Build the domain boundaries
     nb_patches = 4                                             !! The number of boundary segments
     allocate(b_node(nb_patches,2),b_edge(nb_patches,2))        !! Allocate space
     allocate(b_type(nb_patches),b_periodic_parent(nb_patches))
     b_type(:) = (/ 1, 2, 0, 2/)             !! The type of patch (1=wall,2=periodic,0=none,3=inflow,4=outflow)
     b_periodic_parent(:) = (/ 3, 4, 1, 2/)  !! Parents for periodic patches
     b_node(1,:) = (/ -0.5d0*xl, 0.0d0 /)   !! X-Y coordinates of the boundary nodes.
     b_node(2,:) = (/ 0.5d0*xl, 0.0d0 /)
     b_node(3,:) = (/ 0.5d0*xl, yl /)
     b_node(4,:) = (/ -0.5d0*xl, yl /)
     call make_boundary_edge_vectors      !! <---- this routine makes a set of edge vectors which connect the nodes
     xb_min = minval(b_node(:,1));xb_max = maxval(b_node(:,1));yb_min = minval(b_node(:,2));yb_max = maxval(b_node(:,2))    
     call make_boundary_particles    !! <----- allocates xp,yp,up,vp and makes boundary particles for no slip walls
     ipart = n_par_w     

     !! Run through the domain in a square lattice and make liquid particles as desired.
     x = xb_min + 0.5d0*dx
     do while(x<xb_max)!-0.5*dx)
        y = yb_min + dx
        do while(y<h0)
           ipart=ipart+1;xp(ipart)= x;yp(ipart)= y;up(ipart)=0.0d0;vp(ipart)=0.0d0 
           y = y + dx
        enddo       
        x = x + dx
     enddo
     n_par_fw = ipart
!! ------------------------------------------------------------------------------------------------     
  case(2)
     !  Lobovsky dam break...
     !  ooooooooooooooooooooooooooooooo
     !  o                             o
     !  o                             o
     !  o                             o
     !  o                             o
     !  o                             o
     !  o                             o
     !  o                             o
     !  o                             o
     !  o                             o
     !  oxxxxxxxxxxxx                 o
     !  oxxxxxxxxxxxx                 o
     !  oxxxLIQUIDxxx                 o
     !  oxxxxxxxxxxxx                 o
     !  oxxxxxxxxxxxx                 o
     !  ooooooooooooooooooooooooooooooo 

     !! Set domain size, resolution and flow type
     i_open_domain=1               !! 1 indicates a free surface, 0 indicates a closed domain. 
     slip_walls=.false.            !! Slip or no-slip walls     
     h0=0.3d0;xl=1.61d0;yl=0.805d0 !! A lengthscale, and X and Y lengthscales for domain. Must be O(1)
     dx=xl/300.0;dy=dx             !! Particle spacing relative to lengthscale
     nz = 10                       !! Number of copies in 3rd dimension (if 3D)
     h=1.3d0*dx                    !! Initial smoothing length relative to particle spacing

     !! Build domain and boundary particles
     nb_patches = 4
     allocate(b_node(nb_patches,2),b_edge(nb_patches,2))
     allocate(b_type(nb_patches),b_periodic_parent(nb_patches))
     b_type(:) = (/ 1, 1, 0, 1 /) 
     b_periodic_parent(:) = (/ 3, 4, 1, 2 /)
     b_node(1,:) = (/ -0.6d0, 0.0d0 /) 
     b_node(2,:) = (/ xl-0.6d0, 0.0d0 /)
     b_node(3,:) = (/ xl-0.6d0, yl /)
     b_node(4,:) = (/ -0.6d0, yl /)
     call make_boundary_edge_vectors
     xb_min = minval(b_node(:,1));xb_max = maxval(b_node(:,1));yb_min = minval(b_node(:,2));yb_max = maxval(b_node(:,2))
     call make_boundary_particles
     ipart = n_par_w

     !! Liquid particles
     y=yb_min+dx
     do while(y<h0 + 8.0*dx)!yb_max-0.5*dx)
        x = xb_min + dx
        do while(x<xb_max-0.5*dx)
           if(x<0.0d0.and.y<h0) then
              ipart=ipart+1;xp(ipart)= x;yp(ipart)= y;up(ipart)=0.0d0;vp(ipart)=0.0d0        
           end if
           x = x + dx
        enddo
        y = y + dx
     enddo
     n_par_fw = ipart
!! ------------------------------------------------------------------------------------------------     
  case(3)
     !      Taylor-Green Vortices
     !
     !    ----------------------------
     !    |   |    ^        |    ^   |
     !    |   v    |        v    |   |
     !    |<--      <------       <--| 
     !    |                          |
     !    |-->       ------>      -->|
     !    |   |    ^        |    ^   |
     !    |   |    |        |    |   |
     !    |   |    |        |    |   |
     !    |   v    |        v    |   |
     !    |<--      <------       <--|
     !    |                          |
     !    |-->       ------>      -->|
     !    |   |    ^        |    ^   |
     !    |   v    |        v    |   |
     !     --------------------------


     !! Set domain size, resolution and flow type
     i_open_domain=0               !! 1 indicates a free surface, 0 indicates a closed domain. 
     slip_walls=.false.            !! Slip or no-slip walls     
     h0=2.0d0*pi;xl=2.0d0*pi;yl=2.0d0*pi    !! A lengthscale, and X and Y lengthscales for domain. Must be O(1)
     nz = 20                       !! Number of copies in 3rd dimension (if 3D)
     
     
     write(6,*) "Enter resolution (L/dx):"
     read(*,*) ires

     dx=xl/dble(ires);dy=dx              !! Particle spacing relative to lengthscale
     h=1.3d0*dx                    !! Initial smoothing length relative to particle spacing

     !! Build domain and boundary particles
     nb_patches = 4
     allocate(b_node(nb_patches,2),b_edge(nb_patches,2))
     allocate(b_type(nb_patches),b_periodic_parent(nb_patches))    
     b_type(:) = (/ 2, 2, 2, 2 /)  ! all periodic boundaries
     b_periodic_parent(:) = (/ 3, 4, 1, 2 /)
     b_node(1,:) = (/ -0.5d0*h0, -0.5d0*h0 /)
     b_node(2,:) = (/  0.5d0*h0, -0.5d0*h0 /)
     b_node(3,:) = (/  0.5d0*h0,  0.5d0*h0 /)
     b_node(4,:) = (/ -0.5d0*h0,  0.5d0*h0 /)
     call make_boundary_edge_vectors
     xb_min = minval(b_node(:,1));xb_max = maxval(b_node(:,1));yb_min = minval(b_node(:,2));yb_max = maxval(b_node(:,2))
     call make_boundary_particles
     ipart = n_par_w

     !! Make liquid particles
     y=yb_min+dy/2.0d0
     do while (y < yb_max)
        x=xb_min+dx/2.0d0
        do while(x < xb_max)
           ipart=ipart+1
           xp(ipart)=x
           yp(ipart)=y

           !! Taylor-Green velocity field
           up(ipart) = -cos(2.0*pi*x/h0)*sin(2.0*pi*y/h0)
           vp(ipart) = sin(2.0*pi*x/h0)*cos(2.0*pi*y/h0)
                      
           x=x+dx
        enddo
        y=y+dy
     enddo
     n_par_fw=ipart
!! ------------------------------------------------------------------------------------------------     
  case(4)

     !! POISEUILLE FLOW

     ! ooooooooooooooooooooooooooooooooooooooooo
     ! ooooooooooooooooooooooooooooooooooooooooo
     ! |                                        |
     ! |                                        |
     ! |          ----> pressure gradient       |
     ! |                                        |
     ! |                                        |
     ! |                periodic bounds in x    |
     ! |                                        |  
     ! ooooooooooooooooooooooooooooooooooooooooo
     ! ooooooooooooooooooooooooooooooooooooooooo


     !! Set domain size, resolution and flow type
     i_open_domain=0                !! 1 indicates a free surface, 0 indicates a closed domain. 
     slip_walls=.false.            !! Slip or no-slip walls     
     h0=1.0d0;xl=4.0d0*pi;yl=2.0d0  !! A lengthscale, and X and Y lengthscales for domain. Must be O(1)
     dx=xl/300.0;dy=yl/floor(yl/dx) !! Particle spacing relative to lengthscale
     nz = 75                        !! Number of copies in 3rd dimension (if 3D)
     h=1.3d0*max(dx,dy)             !! Initial smoothing length relative to particle spacing
     

     !! Build boundaries and boundary particles
     nb_patches = 4
     allocate(b_node(nb_patches,2),b_edge(nb_patches,2))
     allocate(b_type(nb_patches),b_periodic_parent(nb_patches))
     b_type(:) = (/ 1, 2, 1, 2/)
     b_periodic_parent(:) = (/ 3, 4, 1, 2/)
     b_node(1,:) = (/ -0.5d0*xl, 0.0d0 /)
     b_node(2,:) = (/ 0.5d0*xl, 0.0d0 /)
     b_node(3,:) = (/ 0.5d0*xl, yl /)
     b_node(4,:) = (/ -0.5d0*xl, yl /)
     call make_boundary_edge_vectors
     xb_min = minval(b_node(:,1));xb_max = maxval(b_node(:,1));yb_min = minval(b_node(:,2));yb_max = maxval(b_node(:,2))    
     call make_boundary_particles
     ipart = n_par_w     

      ! -- fluid particles-- 
     x = xb_min + 0.5*dx
     do while(x<xb_max)
        y = yb_min + dy
        do while(y<yb_max - 0.5*dy)
           ipart=ipart+1;xp(ipart)= x;yp(ipart)= y
           up(ipart) = y*(yl-y)*4.0d0/yl/yl          !! Poiseuille velocity field
           vp(ipart)=0.0d0
           y = y + dy
        enddo
        x = x + dx
     enddo
     n_par_fw = ipart
!! ------------------------------------------------------------------------------------------------     
  case(5)
    !! Moin breaking wave case
    ! 
    ! |oooo                                oooo|
    ! |    ooo                          ooo    |
    ! |       oo                      oo       |
    ! |         ooo                ooo         |
    ! |            oooo        oooo            | 
    ! |                oooooooo                | 
    ! |                                        |
    ! |     3rd order Stokes wave              |
    ! |                                        |
    ! |        Periodic in x                   |
    ! |                                        |
    ! |                                        |
    ! ------------------------------------------
    


     !! Additional parameter here for wave profile
     wave_slope = 0.35
     
     ak_ = 0.35
     k_ = 2.0d0*pi
     h_ = 0.5d0
     g_ = 1.0d0
     a_ = ak_/k_
     alpa = 1.0d0/tanh(k_*h_)
     sgma = sqrt(g_*k_*tanh(k_*h_)*(1.0 + (k_**2.0d0)*(a_**2.0d0)*(9.0/8.0*((alpa**2.0d0) - 1.0)*((alpa**2.0d0) - 1.0) &
     + (alpa**2.0d0))))
     b_ = a_*g_/sgma

     !! Set domain size, resolution and flow type
     i_open_domain=1               !! 1 indicates a free surface, 0 indicates a closed domain. 
     slip_walls=.false.            !! Slip or no-slip walls     
     h0=0.5d0;xl=1.0d0;yl=1.0d0    !! A lengthscale, and X and Y lengthscales for domain. Must be O(1)
 
     write(6,*) "Enter resolution (L/dx):"
     read(*,*) ires
     
     dx=xl/dble(ires);dy=dx             !! Particle spacing relative to lengthscale
     nz = 5                       !! Number of copies in 3rd dimension (if 3D)
     h=1.3d0*dx                    !! Initial smoothing length relative to particle spacing
     
     !! Build boundaries and boundary particles
     nb_patches = 4
     allocate(b_node(nb_patches,2),b_edge(nb_patches,2))
     allocate(b_type(nb_patches),b_periodic_parent(nb_patches))
     b_type(:) = (/ 1, 2, 0, 2 /) 
     b_periodic_parent(:) = (/ 3, 4, 1, 2 /)
     b_node(1,:) = (/ -0.5d0*xl, -0.5d0*yl /) 
     b_node(2,:) = (/ 0.5d0*xl, -0.5d0*yl /)
     b_node(3,:) = (/ 0.5d0*xl, 0.5d0*yl /)
     b_node(4,:) = (/ -0.5d0*xl, 0.5d0*yl /)

     call make_boundary_edge_vectors
     xb_min = minval(b_node(:,1));xb_max = maxval(b_node(:,1));yb_min = minval(b_node(:,2));yb_max = maxval(b_node(:,2))
     call make_boundary_particles
     ipart = n_par_w

     !! Liquid particles
     y=yb_min+0.5d0*dx
     do while(y<yb_max)
        x = xb_min + 0.5d0*dx
        do while(x<xb_max-0.1*dx)
           !! eta_fs gives the position of the free surface as a function of X
           !! eta_fs = wave_slope*cos(2.0*pi*x/xl) &
                 !! + 0.5*(wave_slope**2.0d0)*cos(4.0*pi*x/xl) &
                 !! + (3.0/8.0)*(wave_slope**3.0d0)*cos(6.0*pi*x/xl)
           !! eta_fs = eta_fs/(2.0d0*pi)
            eta_fs = (a_*cos(k_*x/xl)) + ak_*(1.0/4.0*alpa*(3.0*(alpa**2.0d0) - 1.0)*(a_**2.0d0)*k_*cos(2.0*k_*x/xl)) &
               + (ak_**2.0d0)*(-3.0/8.0*((alpa**4.0d0) - 3.0*(alpa**2.0d0) + 3.0)*(a_**3.0d0)*(k_**2.0d0)*cos(k_*x/xl) + &
               3.0/64.0*(8.0*(alpa**6.0d0) + ((alpa**2.0d0) - 1.0)*((alpa**2.0d0) - 1.0)) &
               *(a_**3.0d0)*(k_**2.0d0)*cos(3.0*k_*x/xl))
           if(y<=eta_fs) then
              ipart=ipart+1;xp(ipart)= x;yp(ipart)= y       
              !! up(ipart) = (wave_slope/sqrt(2.0*pi))*sqrt(1.0+wave_slope*wave_slope)*cos(2.0*pi*x/xl)*exp(2.0*pi*y/xl)
              
              up(ipart) = b_*cosh(k_*(y/xl + h_))/cosh(k_*h_)*k_*cos(k_*x/xl) + &
                  ak_*3.0*ak_*b_/(8.0*alpa)*((alpa**2.0d0) - 1.0)*((alpa**2.0d0) - 1.0)* &
                  cosh(2.0*k_*(y/xl + h_))*2.0*k_*cos(2.0*k_*x/xl)/cosh(2.0*k_*h_) &
                  + ak_*ak_*1.0/64.0*((alpa**2.0d0) - 1.0)*((alpa**2.0d0) + 3.0)* &
                  (9.*(alpa**2.0d0) - 13.0)*cosh(3.0*k_*(y/xl + h_))/cosh(3.0*k_*h_)*a_*a_*k_*k_*b_*3.0*k_*cos(3.0*k_*x/xl)

              !! vp(ipart) = (wave_slope/sqrt(2.0*pi))*sqrt(1.0+wave_slope*wave_slope)*sin(2.0*pi*x/xl)*exp(2.0*pi*y/xl)
              
              vp(ipart) = b_*k_*sinh(k_*(y/xl + h_))/cosh(k_*h_)*sin(k_*x/xl) & 
                  + ak_*3.0*ak_*b_/(8.0*alpa)*((alpa**2.0d0) - 1.0)*((alpa**2.0d0) - 1.0)* &
                  2.0*k_*sinh(2.0*k_*(y/xl + h_))*sin(2.0*k_*x/xl)/cosh(2.0*k_*h_) & 
                  + ak_*ak_*1.0/64.0*((alpa**2.0d0) - 1.0)*((alpa**2.0d0) + 3.0)* &
                  (9.0*(alpa**2.0d0) - 13.0)*3.0*k_*sinh(3.0*k_*(y/xl + h_))/cosh(3.0*k_*h_)*a_*a_*k_*k_*b_*sin(3.0*k_*x/xl)
           else
              !! Do nothing.
           end if
           x = x + dx
        enddo
        y = y + dx
     enddo
     n_par_fw = ipart
!! ------------------------------------------------------------------------------------------------  
  case(6)
    !! Focussed waves with a wavemaker
    ! 
    ! |oooooooooooooooooooooooooooooooooooooooo|
    ! |                                        |    Free surface initially flat.
    ! |                                        |    Requires compiler flags for 
    ! |                                        |    wavemaker
    ! |                                        |    
    ! |                                        |
    ! |                                        |
    ! |                                        |
    ! ------------------------------------------
    
    !! Might put in a beach...
    


     !! Additional parameter here for wave profile
     wave_slope = 0.0

     !! Set domain size, resolution and flow type
     i_open_domain=1                   !! 1 indicates a free surface, 0 indicates a closed domain. 
     slip_walls=.true.                 !! Slip or no-slip walls     
     h0=0.5d0;xl=12.0d0*h0;yl=2.0d0*h0 !! A lengthscale, and X and Y lengthscales for domain. Must be O(1)
     dx=h0/120.0;dy=dx                  !! Particle spacing relative to lengthscale
     nz = 5                            !! Number of copies in 3rd dimension (if 3D)
     h=1.3d0*dx                        !! Initial smoothing length relative to particle spacing
     
     !! Build boundaries and boundary particles
     nb_patches = 4
     allocate(b_node(nb_patches,2),b_edge(nb_patches,2))
     allocate(b_type(nb_patches),b_periodic_parent(nb_patches))
     b_type(:) = (/ 1,1, 0, 1 /)      
     b_node(1,:) = (/ 0.0d0, -0.5d0*yl /) 
     b_node(2,:) = (/ xl, -0.5d0*yl /)
     b_node(3,:) = (/ xl, 0.5d0*yl /)
     b_node(4,:) = (/ 0.0d0, 0.5d0*yl /)

     call make_boundary_edge_vectors
     xb_min = minval(b_node(:,1));xb_max = maxval(b_node(:,1));yb_min = minval(b_node(:,2));yb_max = maxval(b_node(:,2))
     call make_boundary_particles
     ipart = n_par_w

     !! Liquid particles
     y=yb_min+0.5d0*dx
     do while(y<yb_max)
        x = xb_min + 0.5d0*dx
        xmax_local = xb_max!min(xb_max,xb_max + (y-b_node(3,2))*(xb_max-b_node(2,1))/(b_node(3,2)-b_node(2,2)))
        do while(x<xmax_local)!xb_max)
           !! eta_fs gives the position of the free surface as a function of X
           eta_fs = 0.0d0
           if(y<=eta_fs) then           
              ipart=ipart+1;xp(ipart)= x;yp(ipart)= y       
              up(ipart) = 0.0d0
              vp(ipart) = 0.0d0
           else

           end if
           x = x + dx
        enddo
        y = y + dx
     enddo
     n_par_fw = ipart
!! ------------------------------------------------------------------------------------------------  

  end select

  !! Remainder of the program writes the case data to files
  
  !! Write particle metadata to file
  open(13,file='./IPART')  
  write(13,*) slip_walls
  write(13,*) i_open_domain,"        :open domain (0-closed,1-open)"
  write(13,*) n_par_fw,"        :n_par_fw"
  write(13,*) n_par_w,"        :n_par_w"
  write(13,*) h," :h"
  write(13,*) dx," :dx"
  write(13,*) dy," :dy"  
  write(13,*) nz
  write(13,*) nz*dx
  
  !! Main particle data
  do i=1,n_par_fw
     write(13,*) xp(i), yp(i), up(i), vp(i), dx*dy
  end do
  close(13)
  deallocate(xp, yp)
  deallocate(up, vp)


  !! Write a little to screen
  write(6,*) 'n_par_fw = ',n_par_fw,"n_par_w = ",n_par_w
  write(6,*) 'dx,dy',dx,dy,'h',h


  !! Write boundary data to file.
  open(21,file='./IBOUND')
  write(21,*) nb_patches,"                          :nb_patches"              
  do i=1,nb_patches
     write(21,*) b_node(i,:),"   :b_node",i
     write(21,*) b_edge(i,:),"   :b_edge",i
     write(21,*) b_type(i),"                    :b_type",i
     if(b_type(i)==2)then
        write(21,*) b_periodic_parent(i),"                :b_periodic_parent"
     end if
     if(b_type(i)==3)then
        write(21,*) 1.0d0,"             :U_inflow"  !! May want to change hard coded unit characteristc velocity
     end if
  end do
  close(21)
  deallocate(b_node,b_edge)
  deallocate(b_type,b_periodic_parent)
  
  call system('cp I* ../../.')

  write(6,*) 'END of DATCLASS'
  stop
end program datgen
!! ------------------------------------------------------------------------------------------------
   subroutine make_boundary_edge_vectors
     use kind_parameters
     use common_parameter
     use global_variables
     implicit none
     integer(ikind) ib,ibp1

     do ib = 1,nb_patches ! loop over all boundary patches
        ibp1 = mod(ib,nb_patches) + 1   
        b_edge(ib,:) = b_node(ibp1,:) - b_node(ib,:)  ! calculate b_edge
     end do

     return 
   end subroutine make_boundary_edge_vectors
!! ------------------------------------------------------------------------------------------------
   subroutine make_boundary_particles
     use kind_parameters
     use common_parameter
     use global_variables 

     implicit none
     integer(ikind) :: i,ipart,ib,ibm1
     real(rkind) :: x,y,m_be,tmp,tmp2
     integer(ikind) :: nround,iround

     
     ipart=0
     !! we only allocate memory when we need
     allocate(xp(npar), yp(npar))
     allocate(up(npar), vp(npar))

     if(slip_walls.eqv..false.) then

     !! Wall particles
     do ib=1,nb_patches  ! loop over all boundary patches
        ibm1 = mod(ib+nb_patches-2,nb_patches)+1
        if(abs(b_type(ib))==1)then  ! if it is a wall boundary patch
           m_be = dsqrt(dot_product(b_edge(ib,:),b_edge(ib,:)))
           x = b_node(ib,1);y = b_node(ib,2)
           tmp = 0.0
           do while(tmp<1.0-1.0d-10)   ! move along the patch in increments of dx
              ipart = ipart + 1
              xp(ipart) = x;yp(ipart) = y
              up(ipart) = 0.0d0
              vp(ipart) = 0.0d0
              tmp = tmp + dx/m_be  ! note, in future we should allow for dx/=dy here
              x = b_node(ib,1) + tmp*b_edge(ib,1)
              y = b_node(ib,2) + tmp*b_edge(ib,2)
           end do
        else if(abs(b_type(ibm1))==1)then  ! if the previous boundary patch was a wall (but ib is not a wall...)
           ipart = ipart + 1
           xp(ipart) = b_node(ib,1)  ! place a single particle at the node (to complete the previous wall)
           yp(ipart) = b_node(ib,2)
           m_be = dsqrt(dot_product(b_edge(ibm1,:),b_edge(ibm1,:)))
           up(ipart) = 0.0d0
           vp(ipart) = 0.0d0
        end if
     end do
     end if
     n_par_w=ipart    

    

     write(6,*) 'no. of solid boundary particles:',n_par_w
     return
   end subroutine make_boundary_particles
!! ------------------------------------------------------------------------------------------------
