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
  integer i,j
  double precision h0,eta_fs,yl
  double precision :: f0,sigma,depth,g,k,Li
  real(rkind) :: omega,theta,m1,n1,m,n,d,Lx,Ly,z,dz,zl,kd
  integer(ikind) :: nn
  real(rkind) :: omega_0,omega_1,omega_m,K2,K1,b1,b2,b3,beta1,beta2,beta3
  real(rkind) :: eta_fs_0,eta_fs_1,x_hat,y_hat,t_hat,z_hat,xx,yy,zz
  real(rkind) :: gamma1,gamma3,omega_2,Beta13_const,Beta13,Beta31_const,Beta31,Beta33
  real(rkind) :: b11,b13_const,b13,b31_const,b31,b33,Beta20,eps,eta_fs_2
  real(rkind) :: u_m_hat,v_n_hat,w_hat

  write(*,*) 'Cases: '
  write(*,*) '  case  1:  Monochromatic cross waves'
  write(*,*) '  case  2:  Unidirectional third order Stokes wave'  
  write(*,*) '  '
  write(*,*) 'Input test case number: '
  read(*,*) itest

  select case (itest) 
!! ------------------------------------------------------------------------------------------------
!!  Below are several cases, each goes through the same setup process. Details are HARD CODED, 
!!  so you need to recompile each time you change a case.
!! ------------------------------------------------------------------------------------------------     
  case(1)
     !! Monochromatic crossing waves.
     !! Initial fields come from theory provided by Sam Draycott.
  
     !! CROSSING ANGLE     
     theta  = 45*pi/180; ! angle, 0 is fully standing, 90 is fully travelling. 45 deg is equivalent to two waves at 90 deg.
     if(theta.eq.0.0d0) then
        write(6,*) "Theta cannot be zero. Go back to gen2D... stopping"
        stop
     end if

     !! Dimensional parameters::
     depth = 0.135
     g = 9.81

     !! Dimensionless parameters
     eps = 0.5 !! Wave slope ?  (=kA)     
     kd = pi!0.6281

     k = kd/depth
          
     !! With k now known, find sigma, f0 and L
     sigma = sqrt(k*g*tanh(kd))
     f0 = sigma/(2.0*pi)    
     Li = 2*pi/k  !! Uni-directional wavelength (dimensional ) _ THIS IS CHARACTERISTIC LENGTH SCALE
     
     write(6,*) "sigma",sigma
     
     !! Characteristic scales (only written for output, not actually used)
     write(6,*) "Characteristic scales for wave theory: L, U, T"
     write(6,*) Li,sqrt(g/k),Li/sqrt(g/k)
     write(6,*) "Characteristic scales for bubblysph code: L,U,T"
     write(6,*) Li,sqrt(g/k),Li/sqrt(g/k)
     write(6,*) "With these characteristic scales, we should set Fr=sqrt(1/2pi)"
     
     !! Dimensionless parameters now. Determining angular frequency and domain size
     omega= sigma/((g*k)**0.5) ! this will also change, re-defined later
     m1 = k*sin(theta)
     n1 = k*cos(theta)
     m = sin(theta)
     n = cos(theta)
     d = k*depth
     Lx = 2*pi/m1  !! Note wave people like z to be depth, but the bubbly-sph code takes y as depth. 
     Ly = 2*pi/n1  !! This input file still uses z as depth for wave calculations, then swaps for outputs...  
     
     !! Lx and Ly are DIMENSIONAL.
     write(6,*) "DOMAIN SIZE: ",Lx,Ly           
     !! Length-scales are then made dimensionless with Li

     !! Set domain size, resolution and flow type
     i_open_domain=1                   !! 1 indicates a free surface, 0 indicates a closed domain. 
     slip_walls=.true.                 !! Slip or no-slip walls     
     xl = Lx/Li                         !! Domain size in x,y, and (vertical) z
     yl = Ly/Li
     zl = 2.0d0*depth/Li
     dx = (Li/Li)/100.0d0         !! Spacing in terms of Li, made dimensionless
     dx = xl/nint(xl/dx)                  !! Particle spacing
     dy = yl/nint(yl/dx)               !! slightly adjusted to fit Y domain
     nz = 1     !! This parameter needs outputting, but serves no purpose for this code
     h=1.3d0*dx                        !! Initial smoothing length relative to particle spacing
     
     write(6,*) "Dimensionless domain size and spacings: xl,yl,zl; dx,dy"
     write(6,*) xl,yl,zl
     write(6,*) dx,dy

     !! Build boundaries and boundary particles
     nb_patches = 4
     allocate(b_node(nb_patches,2),b_edge(nb_patches,2))
     allocate(b_type(nb_patches),b_periodic_parent(nb_patches))
     b_type(:) = (/ 1, 2, 0, 2 /) 
     b_periodic_parent(:) = (/ 3, 4, 1, 2 /)
     b_node(1,:) = (/ -0.0d0, -0.5d0*zl /) 
     b_node(2,:) = (/ xl, -0.5d0*zl /)
     b_node(3,:) = (/ xl, 0.5d0*zl /)
     b_node(4,:) = (/ 0.0d0, 0.5d0*zl /)

     call make_boundary_edge_vectors
     xb_min = minval(b_node(:,1));xb_max = maxval(b_node(:,1));yb_min = minval(b_node(:,2));yb_max = maxval(b_node(:,2))
     call make_boundary_particles
     ipart = n_par_w

     !! Constants for wave solution (from Sam's code)
     omega_0 = sqrt(tanh(d))
     omega_1 = 0.0d0
     omega_m = sqrt(tanh(m*d))

     K2 = ((1.0+omega_m**4.0)*((2.0*m**2.0-2.0*n**2.0+1.0)*omega_0**(-3.0) &
        - 3.0*omega_0))/((1.0+omega_m**4.0) - m*(omega_m/omega_0)**2.0)
     K1 = K2/(1.0+omega_m**4.0)
     b1 = (1.0/8.0)*(3.0*omega_0**(-6.0)-omega_0**(-2.0))
     b2 = (1.0/8.0)*(3.0*omega_0**2.0 - omega_0**(-2.0)*(m**2.0-n**2.0)+omega_0*K2)
     b3 = (1.0/8.0)*(omega_0**2.0 - omega_0**(-2.0)*(m**2.0-n**2.0))
     Beta1 = (1.0/8.0)*(-omega_0**(-3.0)+omega_0)
     Beta2 = 3.0*(omega_0**(-7.0)-omega_0)/(16.0*cosh(2.0*d))
     Beta3 = K2/(16.0*cosh(2.0*m*d));   

     gamma1 = (m**2.0+9.0*n**2.0)**(0.5)
     gamma3 = (9.0*m**2.0+n**2.0)**(0.5)

     omega_2 = (1.0/32.0)*(6.0*omega_0**(-7.0)- 8.0*omega_0**(-3.0) - 6.0*omega_0 - 8.0*omega_0**5.0) &
         -(1.0/8.0)*m*omega_0**2.0*omega_m**2.0*K1 &
         +(1.0/16.0)*(omega_0**4.0-4.0*m**2.0+1.0)*K2 &
         +m**2.0*((1.0/32.0)*(3.0*omega_0**(-7.0) &
         - 2.0*omega_0**(-3.0) + 43.0*omega_0) - (1.0/8.0)*(m**2.0-n**2.0)*omega_0**(-3.0)) &
         + n**2.0*((1.0/32.0)*(3.0*omega_0**(-7.0) - 2.0*omega_0**(-3.0) + 5.0*omega_0) &
         - (1.0/8.0)*(m**2.0-n**2.0)*omega_0**(-3.0))


     !! splitting out the constant first to see if helps de-bug
     Beta13_const = (16.0*cosh(gamma1*d)*(gamma1*tanh(gamma1*d)-omega_0**2.0))**(-1.0)
     Beta13 = Beta13_const * &
         (   (-3.0*omega_0**(-7.0) +8.0*omega_0**(-3.0) - 3.0*omega_0+  2.0*omega_0**5.0) &
             + m**2.0 *(-6.0*omega_0**(-7.0) + 4.0*omega_0**(-3.0) - 10.0*omega_0) &
             + n**2.0*(6.0*omega**(-7.0) - 4.0*omega_0**(-3.0) - 2.0*omega_0) &
             + 4.0*n**2.0*(m**2.0-n**2.0)*omega_0**(-3.0)  )
    
     !! splitting out the constant first to see if helps de-bug (FIX on 3rd line
     !! of Beta31, missing omega_m**2)
     Beta31_const = (16.0*cosh(gamma3*d)*(gamma3*tanh(gamma3*d)-9.0*omega_0**2.0))**(-1.0)
     Beta31 =  Beta31_const * &
         (   (-9.0*omega_0**(-7.0) +64.0*omega_0**(-3.0) - 33.0*omega_0+  18.0*omega_0**5.0) &
             + 36.0*m*omega_0**2.0*omega_m**2.0*K1 &
             + 2.0*(3.0*omega_0**4.0 - 8.0*m**2.0-1.0)*K2 &
             + m**2.0*(-18.0*omega_0**(-7.0) +4.0*omega_0**(-3.0) - 30.0*omega_0) &
             +  4.0*m**2.0*(m**2.0-n**2.0)*omega_0**(-3.0) &
             + n**2.0*(18.0*omega_0**(-7.0) - 4.0*omega_0**(-3.0) + 2.0*omega_0)    )
                         

     Beta33 = ((128.0*cosh(3.0*d))**(-1.0))*(1.0+3.0*omega_0**4.0)*(9.0*omega_0**(-13.0) &
            - 22.0*omega_0**(-9.0) + 13.0*omega_0**(-5.0))


     b11 = (1.0/16.0)*(5.0*omega**(-4.0) - 4.0 + 4.0*omega_0**4.0) &
         + (1.0/8.0)*m*omega_0*omega_m**2.0*K1 &
         + (1.0/16.0)*(omega_0**3.0 + 2.0*m**2.0*omega_0**(-1.0) - omega_0**(-1.0))*K2 &
         + (1.0/32.0)*m**2.0 * ( ( 3.0*omega_0**(-8.0) - 2.0*omega_0**(-4.0) - 1.0) &
         - 4.0*(m**2.0-n**2.0)*omega_0**(-4.0)) &
         + (1.0/32.0)*n**2.0 * ( ( 3.0*omega_0**(-8.0) - 2.0*omega_0**(-4.0) - 1.0) &
         + 4.0*(m**2.0-n**2.0)*omega_0**(-4.0))


     !! splitting out the constant first to see if helps de-bug 
     b13_const = (16.0*(gamma1*tanh(gamma1*d)-omega_0**2.0))**(-1.0)

     b13 = (1.0/16.0)* (9.0*omega_0**(-4.0) - 6.0 + 2.0*omega_0**4.0) &
         - (1.0/16.0)*m**2.0*(3.0*omega_0**(-8.0) + 5.0) &
         + (1.0/16.0)*n**2.0*(3.0*omega_0**(-8.0) + 1.0) &
         + b13_const* &
         (   (-3.0*omega_0**(-6.0) + 8.0*omega_0**(-2.0) - 3.0*omega_0**2.0 + 2.0*omega_0**6.0) &
             + m**2.0*(-6.0*omega_0**(-6.0) + 4.0*omega_0**(-2.0) - 10.0*omega_0**2.0) &
             + n**2.0*(6.0*omega_0**(-6.0) - 4.0*omega_0**(-2.0) - 2.0*omega_0**2.0) &
             + 4.0*n**2.0*(m**2.0-n**2.0)*omega_0**(-2.0)  )


     !! splitting out the constant first to see if helps de-bug 
     !%b31_const = (3/16)*(gamma3*tanh(gamma3*d) - 9*omega_0**2); % is this (3/(16*gamma3 etc...)?
     b31_const = 3.0/(16.0*(gamma3*tanh(gamma3*d) - 9.0*omega_0**2.0))
     
     b31 = (1.0/16.0)*(21.0*omega_0**(-4.0) - 10.0 + 6.0*omega_0**4.0) &
         + (3.0/4.0)*m*omega_0*omega_m**2.0*K1 &
         -(1.0/16.0)*m**2.0*(3.0*omega_0**(-8.0) + 5.0) &
         + (1.0/16.0)*n**2.0*(3.0*omega_0**(-8.0)+1.0) &
         + (1.0/8.0)*(omega_0**3.0 - m**2.0*omega_0**(-1.0))*K2 &
         + b31_const* &
         (   (-9.0*omega_0**(-6.0) + 64.0*omega_0**(-2.0) - 33.0*omega_0**2.0 + 18.0*omega_0**6.0) &
             + 36.0*m*omega_0**3.0*omega_m**2.0*K1  &
             + 2.0*K2*(3.0*omega_0**5.0 -8.0*m**2.0*omega_0 - omega_0) &
             + n**2.0*(18.0*omega_0**(-6.0) - 4.0*omega_0**(-2.0) + 2.0*omega_0**2.0) &
             + m**2.0*(-18.0*omega_0**(-6.0) + 4.0*omega_0**(-2.0)  - 30.0*omega_0**2.0) &
             + 4.0*m**2.0*(m**2.0-n**2.0)*omega_0**(-2.0)  )


     b33 = (1.0/16.0)*(-3.0*omega_0**(-8.0) + 21.0*omega_0**(-4.0) - 15.0) &
         + ( (16.0*(tanh(3.0*d) - 3.0*omega_0**2.0))**(-1.0) * ((-27.0*omega_0**(-6.0) &
         + 66.0*omega_0**(-2.0) - 39.0*omega_0**2.0)))

     Beta20 = 0
     
     !! All required constants now built...

     !! Loop over a grid of x and y positions...
     y=0.5d0*dy
     do while(y<yl)
        x=0.5d0*dx
        do while(x<xl)
        
           !! Make x non-dimensional
           x_hat = k*x*Li
           y_hat = k*y*Li
           t_hat = 0.0d0

           !! First order solution
           eta_fs_0 = cos(n*y_hat)*cos(m*x_hat-t_hat)
           
           !! Second order solution
           eta_fs_1 = b1*cos(2.0*n*y_hat)*cos(2*(m*x_hat-t_hat)) + b2*cos(2.0*(m*x_hat-t_hat)) +b3*cos(2.0*n*y_hat)
        
           !! Third order solution
           eta_fs_2 = (b11*cos(n*y_hat) + b13*cos(3.0*n*y_hat))*cos(m*x_hat-t_hat) &
                      + (b31*cos(n*y_hat) + b33*cos(3.0*n*y_hat))*cos(3.0*(m*x_hat - t_hat))                 
eta_fs_2 = 0.0d0                      
                      
 
           !! Evaluate FS elevation 
           eta_fs = eta_fs_0 + eps*eta_fs_1 + 0.5*eps*eps*eta_fs_2           
           eta_fs = eta_fs*eps
           
           !! Re-scale FS elevation to make dimensionless
           eta_fs= eta_fs/(2.0*pi)


           !! Output counter to show some progress...           
!           write(6,*) x/xl,y/yl,eta_fs
           
           !! Set spacing in vertical direction so there's an integer number of particles...
           dz = (eta_fs+0.5*zl)/(dble(nint((eta_fs+0.5*zl)/(0.5*(dx+dy)) - 0.5))+0.5)
           
          
           !! Start at surface
           z = eta_fs
           do while(z>-0.5*zl + 0.2*dz)
              ipart = ipart + 1
              xp(ipart) = x;yp(ipart) = z;zp(ipart) = y !! N.B. here we swap back to y being vertical...

              !! Dimensionless depth (for the velocity)
              z_hat = z*2.0*pi 
              
              !! Scaled coords
              zz = z_hat + d
              yy = n*y_hat
              xx = m*x_hat-t_hat

              u_m_hat = eps*omega_0*(cosh(zz)/sinh(d))*cos(yy)*cos(xx) &
                      + (eps**2.0)*(2.0*Beta2*cosh(2.0*zz)*cos(2.0*yy)*cos(2.0*xx)+ &
                                 2.0*Beta3*cosh(2.0*m*zz)*cos(2.0*xx)) !&
!                      +0.5*(eps**3.0)*(Beta13*cosh(gamma1*zz)*cos(3.0*yy)*cos(xx) &
 !                                   +3.0*Beta31*cosh(gamma3*zz)*cos(yy)*cos(3.0*xx) &
 !                                   +3.0*Beta33*cosh(3.0*zz)*cos(3.0*yy)*cos(3.0*xx)   )

!write(6,*) u_m_hat

              v_n_hat = -eps*omega_0*(cosh(zz)/sinh(d))*sin(yy)*sin(xx) &
                        -eps**2.0*(2.0*Beta2*cosh(2.0*zz)*sin(2.0*yy)*sin(2.0*xx)) !&
!                        -0.5*eps**3.0*(Beta13*cosh(gamma1*zz)*sin(3.0*yy)*sin(xx) &
!                                      +3.0*Beta31*cosh(gamma3*zz)*sin(yy)*sin(3.0*xx) &
!                                      +3.0*Beta33*cosh(3.0*zz)*sin(3.0*yy)*sin(3.0*xx)   )

              w_hat = eps*omega_0*(sinh(zz)/sinh(d))*cos(yy)*sin(xx) &
                    + eps**2.0*(2.0*Beta2*sinh(2.0*zz)*cos(2.0*yy)*sin(2.0*xx) &
                                + 2.0*m*Beta3*sinh(2.0*m*zz)*sin(2.0*xx)) !&
!                    +0.5*eps**3.0*(gamma1*Beta13*sinh(gamma1*zz)*cos(3.0*yy)*sin(xx) &
!                                  +gamma3*Beta31*sinh(gamma3*zz)*cos(yy)*sin(3*xx) &
!                                  +3.0*Beta33*sinh(3.0*zz)*cos(3.0*yy)*sin(3.0*xx)   )


              !! Pass velocities back to arrays (note change of coord system again...)
              !! These velocities are made dimensionless by sqrt(g/k)
              up(ipart) = u_m_hat*m!*eps
              vp(ipart) = w_hat!*eps
              wp(ipart) = v_n_hat*n!*eps
              
              !! Make velocities dimensional
!              up(ipart) = up(ipart)*sqrt(g/k)
!              vp(ipart) = vp(ipart)*sqrt(g/k)              
!              wp(ipart) = wp(ipart)*sqrt(g/k)

              !! Make dimensionless through Froude speed
!              up(ipart) = up(ipart)/sqrt(g*depth)
!              vp(ipart) = vp(ipart)/sqrt(g*depth)              
!              wp(ipart) = wp(ipart)/sqrt(g*depth)
                        
              
              z = z - dz  
           end do
           x = x + dx
        end do
        y = y + dy     
     end do
     n_par_fw = ipart

     write(6,*) "New omega",omega_0 + eps*omega_1 + 0.5*(eps**2.0)*omega_2
!! ------------------------------------------------------------------------------------------------
  case(2)
     !! Monochromatic unidirectional wave
  
     !! Dimensional parameters::
     depth = 0.135
     g = 9.81

     !! Dimensionless parameters
     eps = 0.55 !! Wave slope ?  (=kA)     
     kd = pi!0.6281

     k = kd/depth
          
     !! With k now known, find sigma, f0 and L
     sigma = sqrt(k*g*tanh(kd))
     f0 = sigma/(2.0*pi)    
     Li = 2*pi/k  !! Uni-directional wavelength (dimensional ) _ THIS IS CHARACTERISTIC LENGTH SCALE
     
     write(6,*) "sigma",sigma
     
     !! Characteristic scales (only written for output, not actually used)
     write(6,*) "Characteristic scales for wave theory: L, U, T"
     write(6,*) Li,sqrt(g/k),Li/sqrt(g/k)
     write(6,*) "Characteristic scales for bubblysph code: L,U,T"
     write(6,*) Li,sqrt(g/k),Li/sqrt(g/k)
     write(6,*) "With these characteristic scales, set Fr=sqrt(1/2pi)"
     
     !! Dimensionless parameters now. Determining angular frequency and domain size
     omega= sigma/((g*k)**0.5) ! this will also change, re-defined later
     d = k*depth

     !! Set domain size, resolution and flow type
     i_open_domain=1                   !! 1 indicates a free surface, 0 indicates a closed domain. 
     slip_walls=.true.                 !! Slip or no-slip walls     
     xl = 1.0d0
!     yl = 0.06d0
     zl = 2.0d0*depth/Li
     dx = (Li/Li)/150.0d0         !! Spacing in terms of Li, made dimensionless
     dx = xl/nint(xl/dx)                  !! Particle spacing
     yl = dx*5.0d0
     dy = yl/nint(yl/dx)               !! slightly adjusted to fit Y domain
     nz = 1     !! This parameter needs outputting, but serves no purpose for this code
     h=1.3d0*dx                        !! Initial smoothing length relative to particle spacing
     
     write(6,*) "Dimensionless domain size and spacings: xl,yl,zl; dx,dy"
     write(6,*) xl,yl,zl
     write(6,*) dx,dy
     
     !! Build boundaries and boundary particles
     nb_patches = 4
     allocate(b_node(nb_patches,2),b_edge(nb_patches,2))
     allocate(b_type(nb_patches),b_periodic_parent(nb_patches))
     b_type(:) = (/ 1, 2, 0, 2 /) 
     b_periodic_parent(:) = (/ 3, 4, 1, 2 /)
     b_node(1,:) = (/ -0.5d0*xl, -0.5d0*zl /) 
     b_node(2,:) = (/ 0.5d0*xl, -0.5d0*zl /)
     b_node(3,:) = (/ 0.5d0*xl, 0.5d0*zl /)
     b_node(4,:) = (/ -0.5d0*xl, 0.5d0*zl /)

     call make_boundary_edge_vectors
     xb_min = minval(b_node(:,1));xb_max = maxval(b_node(:,1));yb_min = minval(b_node(:,2));yb_max = maxval(b_node(:,2))
     call make_boundary_particles
     ipart = n_par_w

     !! Loop over a grid of x and y positions...
     y=0.5d0*dy
     do while(y<yl)
        x=xb_min + 0.5d0*dx
        do while(x<xb_max)
        


           eta_fs = eps*cos(2.0*pi*x/xl) &
                  + 0.5*(eps**2.0d0)*cos(4.0*pi*x/xl) &
                  + (3.0/8.0)*(eps**3.0d0)*cos(6.0*pi*x/xl)
           eta_fs = eta_fs/(2.0d0*pi)  !! eta normalised by wavelength
                          
           !! Set spacing in vertical direction so there's an integer number of particles...
           dz = (eta_fs+0.5*zl)/(dble(nint((eta_fs+0.5*zl)/(0.5*(dx+dy)) - 0.5))+0.5)
                     
           !! Start at surface
           z = eta_fs
           do while(z>-0.5*zl + 0.2*dz)
              ipart = ipart + 1
              xp(ipart) = x;yp(ipart) = z;zp(ipart) = y !! N.B. here we swap back to y being vertical...


              up(ipart) = eps*sqrt(1.0+eps*eps)*cos(2.0*pi*x)*exp(2.0*pi*z)
              vp(ipart) = eps*sqrt(1.0+eps*eps)*sin(2.0*pi*x)*exp(2.0*pi*z)              
!              up(ipart) = up(ipart)/sqrt(2.0*pi)
!              vp(ipart) = vp(ipart)/sqrt(2.0*pi)    
                                   
              z = z - dz  
           end do
           x = x + dx
        end do
        y = y + dy     
     end do
     n_par_fw = ipart
        
!! ------------------------------------------------------------------------------------------------  

  end select

  !! Remainder of the program writes the case data to files
  
  !! Write particle data to file
  open(13,file='./IPART')
    
  write(13,*) slip_walls
  write(13,*) i_open_domain
  write(13,*) n_par_fw
  write(13,*) n_par_w
  write(13,*) h
  write(13,*) dx
  write(13,*) dy,dy
  write(13,*) nz
  write(13,*) yl  
  do i=1,n_par_fw
     write(13,*) xp(i), yp(i),zp(i), up(i), vp(i),wp(i), dx*dy
  end do
  close(13)
  deallocate(xp, yp,zp)
  deallocate(up, vp,wp)


  !! Write a little to screen
  write(6,*) 'n_par_fw = ',n_par_fw,"n_par_w = ",n_par_w
  write(6,*) 'dx,dy',dx,dy,'h',h

  !! Write control data to file


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
     allocate(xp(npar), yp(npar),zp(npar))
     allocate(up(npar), vp(npar),wp(npar))

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
