module output
  !! This module contains routines to output field data and statistics, both every time-step, and
  !! every dt_out.
  use kind_parameters
  use common_parameter
  use common_2d
  use sphtools
  implicit none
  integer(ikind) :: n_bub_out,nslice,ngrid,nsurf,n_out,nelev

  private
  public :: output_fields,output_statistics,setup_stats_files
contains
!! ------------------------------------------------------------------------------------------------
  subroutine output_statistics
     !! This routine calls routines for specific statistics. It is called every time-step.
  
     if(.true.) call global_energy_stats     
          
     if(.true.) call bubble_entrainment_stats     
  
     if(.false.) call plume_slice_data      
     
     if(.false.) call taylor_green_error  
     
     if(.true.) call fs_elev
  
  
     return
  end subroutine output_statistics
!! ------------------------------------------------------------------------------------------------
  subroutine output_fields
    !! This routine calls routines for specific field outputs - particles, slices, grid, surf,
    !! bubbles, and also does regular DUMP files for restarts.
    !! This routine is only called every dt_out
    integer(ikind) :: idp,nf,i,i_out,int_out

    call particle_output
    
    call bubble_output
    
    if(.false.) call channel_flow_statistics        
      
    if(.false.) call fraga_analysis    

#if dim3
    !! Interpolate (v. low order) to a 2d x-y grid for 3D simulations
    if(.true.) call slice_output

    !! Full 3D grid
    if(.false.) call grid_output !! SPH interpolation into a grid!

    !! Full 3D grid, velocity only (used to find K.E. spectrum)
    if(.false.) call grid_output_u !! MLS interpolation into a grid!
    
    !! Elevation on 2D grid
    if(.true.) call elevation_output    
    
    !! 3D grid of points near surface - used to get nice 3d free-surface plots
    if(.false.) call surf_output
#endif     
     
    !! write and push some time data
    !! time,#particles,itime,dt,#elements in slice, #elements in grid,# bubbles out,# elements in surface,#elements in welev
    write(94,*) time,n_out,itime,dt,nslice,ngrid,n_bub_out,nsurf,nelev
    flush(94)



    !!  Dump output for restarts
    open(100,file='DUMP')
    write(100,*)itime,time,n_dump
    write(100,*)n_par_fwm-n_inbin,n_par_fw-n_inbin,n_par_w,np_inflow
    write(100,*)xb_min,xb_max,yb_min,yb_max, dx, dy, av_conc
    do i=1,n_par_fw
       if(inbin(i)) cycle
       write(100,*) rp(i,:),up(i,:),p(i)
    enddo
    !! Indicate the inflow particles
    if(np_inflow/=0)then
       do i=1,np_inflow
          write(100,*) p_inflow(i)
       end do
    end if
#ifdef bubbles
    do i=1,n_bub
       write(100,*) rb(i,:),ub(i,:)
    end do
#endif     

    close(100)
    
    !! Deallocate properties no longer required.
    deallocate(vort)
          
    return
  end subroutine output_fields
!! ------------------------------------------------------------------------------------------------
  subroutine particle_output

    integer(ikind) :: idp,nf,i,i_out,int_out
    character(70) fnm1

    ! File naming stuff
    nf = 22
    idp = n_dump
    ! Create the name of the output file
    if( idp < 10 ) then 
       write(fnm1,'(A21,I1)') './data_directory/PART', idp
    else if( idp < 100 ) then 
       write(fnm1,'(A21,I2)') './data_directory/PART', idp
    else if( idp < 1000 ) then
       write(fnm1,'(A21,I3)') './data_directory/PART', idp
    else
       write(fnm1,'(A21,I4)') './data_directory/PART', idp
    end if

    ! Open the file
    open(nf+1,file=fnm1,status='unknown')

    ! Will we output mirror particles?
    if(output_mirrors) then
       n_out=n_par_fw + nmirror
    else
       n_out=n_par_fw 
    endif

    ! Output the data
    i_out = 0
    do i=1,n_out
       if(inbin(i)) cycle  ! don't waste space outputting parts in bin
       !! Output surface or colour function
       int_out = n_surf(i)
        
#if dim3
       if(abs(rp(i,3))>dx.and.n_surf(i)==0) cycle  !! This line means only plot particles at centreline (z=0), and in FS
       i_out = i_out + 1        
       write(nf+1,200) rp(i,1),rp(i,2),rp(i,3),up(i,1),up(i,2),up(i,3),P(i),conc(i),a_out(i),one-a0(i),int_out,vort(i) 
#else        
       i_out = i_out + 1
       write(nf+1,200) rp(i,1),rp(i,2),zero,up(i,1),up(i,2),zero,P(i),conc(i),a_out(i),one-a0(i),int_out,vort(i) 
#endif      
    end do
         
     
  200 FORMAT(' ',10(E11.5,1X),I2,1X,E11.5)       
    close(nf+1)
    n_out = i_out

   

     return
  end subroutine particle_output
!! ------------------------------------------------------------------------------------------------
  subroutine bubble_output

    integer(ikind) :: idp,nf,i,i_out,int_out
    character(70) fnm2

    nf = 22
    idp = n_dump

    ! Create the name of the output file
#ifdef bubbles
    if( idp < 10 ) then 
       write(fnm2,'(A23,I1)') './data_directory/BUBBLE', idp
    else if( idp < 100 ) then 
       write(fnm2,'(A23,I2)') './data_directory/BUBBLE', idp
    else if( idp < 1000 ) then
       write(fnm2,'(A23,I3)') './data_directory/BUBBLE', idp
    else
       write(fnm2,'(A23,I4)') './data_directory/BUBBLE', idp
    end if
#endif     

    !! Output the bubble positions, velocities and sizes...
#ifdef bubbles
    open(nf+1,file=fnm2,status='unknown')
    n_bub_out=n_bub;i_out = 0
    do i=1,n_bub_out
       if(b_inbin(i)) cycle  ! don't waste space outputting bubbles in bin     
       i_out = i_out + 1
#if dim3
       !! outputs pubble position (x,y,z), bubble velocity (u,v,w), bubble size, and a scalar
       write(nf+1,300) rb(i,1),rb(i,2),rb(i,3),ub(i,1),ub(i,2),ub(i,3),radb(i)/Rhinze,b_out(i)
#else        
       write(nf+1,300) rb(i,1),rb(i,2),zero,ub(i,1),ub(i,2),zero,radb(i)/Rhinze,b_out(i)
#endif      
  300 FORMAT(' ',8(E11.5,1X))  
    enddo
    close(nf+1)
    n_bub_out = i_out

#endif

     return
  end subroutine bubble_output
!! ------------------------------------------------------------------------------------------------
  subroutine slice_output
#if dim3
    integer idp,nf,i,n_out,i_out
    character(70) fnm1
    integer :: ii,jj,nxs,nys,kk,j,i0,j0,nss
    real(rkind),dimension(:),allocatable :: xs,ys,us,vs,ws,ps,vorts,alphs,surfs,vfs
    real(rkind),dimension(dims) :: rij
    real(rkind) :: rad,qq,wtmp,tta,ttb,kcoeff
  

    ! File naming stuff
    nf = 22
    idp = n_dump
  
    ! Create the name of the output file
    if( idp < 10 ) then 
       write(fnm1,'(A22,I1)') './data_directory/SLICE', idp
    else if( idp < 100 ) then 
       write(fnm1,'(A22,I2)') './data_directory/SLICE', idp
    else if( idp < 1000 ) then
       write(fnm1,'(A22,I3)') './data_directory/SLICE', idp
    else
       write(fnm1,'(A22,I4)') './data_directory/SLICE', idp
    end if
  
    ! Open the file
    open(nf+1,file=fnm1,status='unknown')
    
    !! Make slice grid
    if(i_open_domain==1) then
       ymax = maxval(rp(n_par_w+1:n_par_fw,2)+3.0*dx)
    else
       ymax = yb_max
    end if
  
    nxs = 2*ceiling((xb_max-xb_min)/dx)
    nys = 2*ceiling((ymax-yb_min)/dy)
    nslice=nxs*nys
    
    nss = ceiling(two*sup_size/dx)
     
    allocate(xs(nslice),ys(nslice))
    kk=0
    do ii=1,nxs
       do jj=1,nys
          kk=kk+1
          xs(kk) = xb_min + dble(ii-1)*(xb_max-xb_min)/dble(nxs-1)
          ys(kk) = yb_min + dble(jj-1)*(ymax-yb_min)/dble(nys-1)        
       end do
    end do
  

    !! Loop over all particles, then find a smallish range of slice-nodes to loop over
    allocate(us(nslice),vs(nslice),ws(nslice),ps(nslice),vorts(nslice),alphs(nslice),surfs(nslice),vfs(nslice))
    us=zero;vs=zero;ws=zero;ps=zero;vorts=zero;alphs=zero;surfs=zero;vfs=zero
!   !$omp parallel do private(i0,j0,ii,jj,kk,rij,rad,wtmp) reduction(+:us,vs,ws,ps,vorts,alphs)
    do i=1,n_par_fwm
       kcoeff = dv/vol(i)
       if(rp(i,3)>sup_size) cycle 
       i0 = floor(nxs*(rp(i,1)-xb_min)/(xb_max-xb_min));i0=min(i0,nxs);i0=max(i0,1)
       j0 = floor(nys*(rp(i,2)-yb_min)/(ymax-yb_min));j0=min(j0,nys);j0=max(j0,1)
       do ii=max(i0-nss,1),min(i0+nss,nxs)
          do jj=max(j0-nss,1),min(j0+nss,nys)
             kk = nys*(ii-1) + jj !! Linear index..
             rij = (/xs(kk),ys(kk),zero/) - rp(i,:)
             rad = sqrt(dot_product(rij,rij));qq=rad/h(i)
             wtmp = wab(qq)*kcoeff        

             us(kk) = us(kk) + dv*wtmp*up(i,1)
             vs(kk) = vs(kk) + dv*wtmp*up(i,2)
             ws(kk) = ws(kk) + dv*wtmp*up(i,3)
             ps(kk) = ps(kk) + dv*wtmp*p(i)                                 
             vorts(kk)=vorts(kk) + dv*wtmp*vort(i)
             alphs(kk)=alphs(kk) + dv*wtmp*a_out(i)
             surfs(kk) = surfs(kk) + dv*wtmp
             vfs(kk) = vfs(kk) + dv*wtmp*(one-a0(i))
          end do
       end do     
    end do
!    !$omp end parallel do

    !! Write to output
    nss = nslice
    nslice = 0
    do kk=1,nss   !! Loop over all points in slice
       if(surfs(kk)>1.0d-1) then
          nslice = nslice + 1
          write(nf+1,400) xs(kk),ys(kk),zero,us(kk),vs(kk),ws(kk),ps(kk),vorts(kk),alphs(kk),surfs(kk),vfs(kk)
       end if
    400 FORMAT(' ',11(E11.5,1X))       
    end do
  
    deallocate(xs,ys,us,vs,ws,ps,alphs,vorts,surfs,vfs)

    close(nf+1)
#endif
    return
  end subroutine slice_output
!! ------------------------------------------------------------------------------------------------
  subroutine grid_output
#if dim3
    integer idp,nf,i,n_out,i_out
    character(70) fnm1
    integer :: ii,jj,nxs,nys,kk,j,i0,j0,nss,ll,nzs,l0
    real(rkind),dimension(:),allocatable :: xs,ys,us,vs,ws,ps,vorts,alphs,zs
    real(rkind),dimension(dims) :: rij
    real(rkind) :: rad,qq,wtmp,tta,ttb,kcoeff
  

    ! File naming stuff
    nf = 22
    idp = n_dump

    ! Create the name of the output file
    if( idp < 10 ) then 
       write(fnm1,'(A22,I1)') './data_directory/WGRID', idp
    else if( idp < 100 ) then 
       write(fnm1,'(A22,I2)') './data_directory/WGRID', idp
    else if( idp < 1000 ) then
       write(fnm1,'(A22,I3)') './data_directory/WGRID', idp
    else
       write(fnm1,'(A22,I4)') './data_directory/WGRID', idp
    end if

    ! Open the file
    open(nf+1,file=fnm1,status='unknown')
  
    !! Make slice grid
    if(i_open_domain==1) then
       ymax = maxval(rp(n_par_w+1:n_par_fw,2))
    else
       ymax = yb_max
    end if
  
    nxs = ceiling((xb_max-xb_min)/dx)
    nys = ceiling((ymax-yb_min)/dx)
    nzs = ceiling((zb_max-zb_min)/dx)
    ngrid=nxs*nys*nzs
   
    nss = ceiling(two*sup_size/dx)
   
    allocate(xs(ngrid),ys(ngrid),zs(ngrid))
    kk=0
    do ii=1,nxs
       do jj=1,nys
          do ll=1,nzs
             kk=kk+1
             xs(kk) = xb_min + dble(ii-1)*(xb_max-xb_min)/dble(nxs-1)
             ys(kk) = yb_min + dble(jj-1)*(ymax-yb_min)/dble(nys-1)        
             zs(kk) = zb_min + dble(ll-1)*(zb_max-zb_min)/dble(nzs-1)
          end do
       end do
    end do
  

    !! Loop over all particles, then find a smallish range of slice-nodes to loop over
    allocate(us(ngrid),vs(ngrid),ws(ngrid),ps(ngrid),vorts(ngrid),alphs(ngrid))
    us=zero;vs=zero;ws=zero;ps=zero;vorts=zero;alphs=zero
!    !$omp parallel do private(i0,j0,ii,jj,kk,rij,rad,wtmp) reduction(+:us,vs,ws,ps,vorts,alphs)
    do i=1,n_par_fwm
       kcoeff = dv/vol(i)
       i0 = floor(nxs*(rp(i,1)-xb_min)/(xb_max-xb_min));i0=min(i0,nxs);i0=max(i0,1)
       j0 = floor(nys*(rp(i,2)-yb_min)/(ymax-yb_min));j0=min(j0,nys);j0=max(j0,1)
       l0 = floor(nzs*(rp(i,3)-zb_min)/(zb_max-zb_min));l0=min(l0,nzs);l0=max(l0,1)
       do ii=max(i0-nss,1),min(i0+nss,nxs)
          do jj=max(j0-nss,1),min(j0+nss,nys)
             do ll=max(l0-nss,1),min(l0+nss,nzs)
                kk = nzs*nys*(ii-1) + nzs*(jj-1) + ll !! Linear index..
                rij = (/xs(kk),ys(kk),zs(kk)/) - rp(i,:)
                rad = sqrt(dot_product(rij,rij));qq=rad/h(i)
                wtmp = wab(qq)*kcoeff        
  
                us(kk) = us(kk) + vol(i)*wtmp*up(i,1)
                vs(kk) = vs(kk) + vol(i)*wtmp*up(i,2)
                ws(kk) = ws(kk) + vol(i)*wtmp*up(i,3)
                ps(kk) = ps(kk) + vol(i)*wtmp*p(i)                                 
                vorts(kk)=vorts(kk) + vol(i)*wtmp*vort(i)
                alphs(kk)=alphs(kk) + vol(i)*wtmp*a_out(i)
             end do
          end do
       end do     
    end do
!    !$omp end parallel do

    !! Write to output
    do kk=1,ngrid   !! Loop over all points in slice
       write(nf+1,500) xs(kk),ys(kk),zs(kk),us(kk),vs(kk),ws(kk),ps(kk),vorts(kk),alphs(kk)
    500 FORMAT(' ',9(E11.5,1X))       
    end do
  
    deallocate(xs,ys,us,vs,ws,ps,alphs,vorts,zs)

    close(nf+1)
#endif
    return
  end subroutine grid_output
!! ------------------------------------------------------------------------------------------------
  subroutine elevation_output
#if dim3
    integer idp,nf,i,n_out,i_out
    character(70) fnm1
    integer :: ii,jj,nxs,nys,kk,j,i0,j0,nss,ll,nzs,l0
    real(rkind),dimension(:),allocatable :: xs,zs,us,vs,ws,es,vols
    real(rkind),dimension(dims) :: rij
    real(rkind) :: rad,qq,wtmp,tta,ttb,kcoeff
  

    ! File naming stuff
    nf = 22
    idp = n_dump

    ! Create the name of the output file
    if( idp < 10 ) then 
       write(fnm1,'(A22,I1)') './data_directory/WELEV', idp
    else if( idp < 100 ) then 
       write(fnm1,'(A22,I2)') './data_directory/WELEV', idp
    else if( idp < 1000 ) then
       write(fnm1,'(A22,I3)') './data_directory/WELEV', idp
    else
       write(fnm1,'(A22,I4)') './data_directory/WELEV', idp
    end if

    ! Open the file
    open(nf+1,file=fnm1,status='unknown')
  
    !! Grid size 
    nxs = ceiling((xb_max-xb_min)/dx)
    nzs = ceiling((zb_max-zb_min)/dx)
    nelev=nxs*nzs
   
    nss = ceiling(two*sup_size/dx)
   
    allocate(xs(nelev),zs(nelev))
    kk=0
    do ii=1,nxs
       do ll=1,nzs
          kk=kk+1
          xs(kk) = xb_min + dble(ii-1)*(xb_max-xb_min)/dble(nxs-1)
          zs(kk) = zb_min + dble(ll-1)*(zb_max-zb_min)/dble(nzs-1)
       end do
    end do
  

    !! Loop over all particles, then find a smallish range of slice-nodes to loop over
    allocate(us(nelev),vs(nelev),ws(nelev),es(nelev),vols(nelev))
    us=zero;vs=zero;ws=zero;es=zero;vols=zero
    do i=1,n_par_fwm
       if(n_surf(i).eq.1) then  !! Only consider free surface particles...
          kcoeff = dv/vol(i)
          i0 = floor(nxs*(rp(i,1)-xb_min)/(xb_max-xb_min));i0=min(i0,nxs);i0=max(i0,1)
          l0 = floor(nzs*(rp(i,3)-zb_min)/(zb_max-zb_min));l0=min(l0,nzs);l0=max(l0,1)
          do ii=max(i0-nss,1),min(i0+nss,nxs)
             do ll=max(l0-nss,1),min(l0+nss,nzs)
                kk = nzs*(ii-1) +  ll !! Linear index..
                rij = (/xs(kk),0.0_rkind,zs(kk)/) - rp(i,:)     !! Taking 2D SPH interpolation over surface
                rij(2) = 0.0d0
                rad = sqrt(dot_product(rij,rij));qq=rad/h(i)
                wtmp = wab(qq)*kcoeff            
               
                vols(kk) = vols(kk) + vol(i)*wtmp
                us(kk) = us(kk) + vol(i)*wtmp*up(i,1)
                vs(kk) = vs(kk) + vol(i)*wtmp*up(i,2)
                ws(kk) = ws(kk) + vol(i)*wtmp*up(i,3)
                es(kk) = es(kk) + vol(i)*wtmp*rp(i,2)                                 
             end do
          end do     
       end if
    end do
!    !$omp end parallel do

    !! Write to output
    do kk=1,nelev   !! Loop over all points in slice
       us(kk) = us(kk)/vols(kk)
       vs(kk) = vs(kk)/vols(kk)
       ws(kk) = ws(kk)/vols(kk)
       es(kk) = es(kk)/vols(kk)
       write(nf+1,500) xs(kk),zs(kk),us(kk),vs(kk),ws(kk),es(kk)
    500 FORMAT(' ',9(E11.5,1X))       
    end do
  
    deallocate(xs,zs,us,vs,ws,es,vols)

    close(nf+1)
#endif
    return
  end subroutine elevation_output  
!! ------------------------------------------------------------------------------------------------  
  subroutine grid_output_u
    use svd_lib  
#if dim3
    integer idp,nf,i,n_out,i_out
    character(70) fnm1
    integer :: ii,jj,nxs,nys,kk,j,i0,j0,nss,ll,nzs,l0
    real(rkind),dimension(:),allocatable :: xs,ys,us,vs,ws,zs,vols
    real(rkind),dimension(dims) :: rij
    real(rkind) :: rad,qq,wtmp,tta,ttb,kcoeff,ff1,oosqrt2,wij,x,y,z
    integer(ikind),dimension(:),allocatable :: grid_neighbourcount
    integer(ikind),dimension(:,:),allocatable :: grid_neighbourlist
    integer(ikind) :: nsize,nsize_local,i1,k
    real(dkind),dimension(:,:),allocatable :: amat
    real(dkind),dimension(:),allocatable :: xvec,wvec,psivec
    
    oosqrt2 = one/sqrt(2.0d0)
    
    nsize=20
    allocate(amat(nsize,nsize))
    allocate(xvec(nsize),wvec(nsize),psivec(nsize))
  

    ! File naming stuff
    nf = 22
    idp = n_dump

    ! Create the name of the output file
    if( idp < 10 ) then 
       write(fnm1,'(A22,I1)') './data_directory/UGRID', idp
    else if( idp < 100 ) then 
       write(fnm1,'(A22,I2)') './data_directory/UGRID', idp
    else if( idp < 1000 ) then
       write(fnm1,'(A22,I3)') './data_directory/UGRID', idp
    else
       write(fnm1,'(A22,I4)') './data_directory/UGRID', idp
    end if

    ! Open the file
    open(nf+1,file=fnm1,status='unknown')
  
    !! Make slice grid
    if(i_open_domain==1) then
       ymax = maxval(rp(n_par_w+1:n_par_fw,2))
    else
       ymax = yb_max
    end if
  
    nxs = ceiling((xb_max-xb_min)/dx)
    nys = ceiling((ymax-yb_min)/dx)
    nzs = ceiling((zb_max-zb_min)/dx)
    ngrid=nxs*nys*nzs
   
    nss = ceiling(two*sup_size/dx)
   
    allocate(xs(ngrid),ys(ngrid),zs(ngrid))
    kk=0
    do ii=1,nxs
       do jj=1,nys
          do ll=1,nzs
             kk=kk+1
             xs(kk) = xb_min + dble(ii-1)*(xb_max-xb_min)/dble(nxs-1)
             ys(kk) = yb_min + dble(jj-1)*(ymax-yb_min)/dble(nys-1)        
             zs(kk) = zb_min + dble(ll-1)*(zb_max-zb_min)/dble(nzs-1)
          end do
       end do
    end do
  

    !! Loop over all particles, then find a smallish range of slice-nodes to loop over
    allocate(us(ngrid),vs(ngrid),ws(ngrid),vols(ngrid))
    us=zero;vs=zero;ws=zero;vols=zero
!    !$omp parallel do private(i0,j0,ii,jj,kk,rij,rad,wtmp) reduction(+:us,vs,ws,vols)
    
    allocate(grid_neighbourcount(ngrid),grid_neighbourlist(nplink,ngrid))
    grid_neighbourcount=0;grid_neighbourlist=0
    do i=1,n_par_fwm
       kcoeff = dv/vol(i)
       i0 = floor(nxs*(rp(i,1)-xb_min)/(xb_max-xb_min));i0=min(i0,nxs);i0=max(i0,1)
       j0 = floor(nys*(rp(i,2)-yb_min)/(ymax-yb_min));j0=min(j0,nys);j0=max(j0,1)
       l0 = floor(nzs*(rp(i,3)-zb_min)/(zb_max-zb_min));l0=min(l0,nzs);l0=max(l0,1)
       do ii=max(i0-nss,1),min(i0+nss,nxs)
          do jj=max(j0-nss,1),min(j0+nss,nys)
             do ll=max(l0-nss,1),min(l0+nss,nzs)
                kk = nzs*nys*(ii-1) + nzs*(jj-1) + ll !! Linear index..
                
                
                rij = (/xs(kk),ys(kk),zs(kk)/) - rp(i,:)
                rad = sqrt(dot_product(rij,rij));qq=rad/h(i)

                !! Build the grid neighbours
                if(qq.le.2.0d0) then
                   grid_neighbourcount(kk) = grid_neighbourcount(kk) + 1
                   grid_neighbourlist(grid_neighbourcount(kk),kk) = i
                endif

             end do
          end do
       end do     
    end do
!    !$omp end parallel do

    !! Loop over all grid points, and build MLS interpolation
    !$omp parallel do private(k,j,xvec,wvec,amat,psivec,x,y,z,i1,nsize_local)
    do i=1,ngrid
    
       !! Loop over SPH particle neighbours
       amat = zero
       nsize_local = nsize
       do k=1,grid_neighbourcount(i)
          j=grid_neighbourlist(k,i)     !! Index of SPH particle neighbour
          
          !! Find x_{i} (grid) - x_{j} (particle)
          x = rp(j,1) - xs(i)
          y = rp(j,2) - ys(i)
          z = rp(j,3) - zs(i)                    
          
          !! Scale x,y,z for hermite
          rad = sqrt(x*x + y*y)/h0;qq=rad
          ff1 = Wab(qq)
          x=x/h0;y=y/h0;z=z/h0             
                    
          !! Build Vector of Taylor monomials
          xvec(1) = one
          xvec(2) = x
          xvec(3) = y
          xvec(4) = z
          xvec(5) = (1.0/2.0)*x*x
          xvec(6) = x*y
          xvec(7) = (1.0/2.0)*y*y
          xvec(8) = y*z
          xvec(9) = (1.0/2.0)*z*z
          xvec(10) = z*x
          xvec(11) = (1.0/6.0)*x*x*x
          xvec(12) = (1.0/2.0)*x*x*y
          xvec(13) = (1.0/2.0)*x*y*y
          xvec(14) = (1.0/6.0)*y*y*y
          xvec(15) = (1.0/2.0)*y*y*z
          xvec(16) = (1.0/2.0)*y*z*z
          xvec(17) = (1.0/6.0)*z*z*z
          xvec(18) = (1.0/2.0)*z*z*x
          xvec(19) = (1.0/2.0)*z*x*x
          xvec(20) = x*y*z 
                    
          !! Build vector of Basis functions
          !Re-scale x,y,z
          x=x*oosqrt2;y=y*oosqrt2;z=z*oosqrt2          
          wvec(1) = ff1*one
          wvec(2) = ff1*Hermite1(x)*oosqrt2
          wvec(3) = ff1*Hermite1(y)*oosqrt2
          wvec(4) = ff1*Hermite1(z)*oosqrt2                    
          wvec(5) = ff1*Hermite2(x)*half
          wvec(6) = ff1*Hermite1(x)*Hermite1(y)*half
          wvec(7) = ff1*Hermite2(y)*half
          wvec(8) = ff1*Hermite1(y)*Hermite1(z)*half
          wvec(9) = ff1*Hermite2(z)*half
          wvec(10)= ff1*Hermite1(z)*Hermite1(x)*half
          wvec(11)= ff1*Hermite3(x)*half*oosqrt2
          wvec(12)= ff1*Hermite2(x)*Hermite1(y)*half*oosqrt2
          wvec(13)= ff1*Hermite2(y)*Hermite1(x)*half*oosqrt2
          wvec(14)= ff1*Hermite3(y)*half*oosqrt2
          wvec(15)= ff1*Hermite2(y)*Hermite1(z)*half*oosqrt2
          wvec(16)= ff1*Hermite2(z)*Hermite1(y)*half*oosqrt2
          wvec(17)= ff1*Hermite3(z)*half*oosqrt2
          wvec(18)= ff1*Hermite2(z)*Hermite1(x)*half*oosqrt2
          wvec(19)= ff1*Hermite2(x)*Hermite1(z)*half*oosqrt2
          wvec(20)= ff1*Hermite1(x)*Hermite1(y)*Hermite1(z)*half*oosqrt2
          
          do i1=1,nsize
             amat(i1,:) = amat(i1,:) + xvec(i1)*wvec(:)   !! Contribution to LHS for this interaction
          end do   
                             
       end do
    
       
       !! Solve linear system to get Psi
       psivec = zero;psivec(1) = one
       call svd_solve(amat,nsize,psivec)       
       
       !! Second loop over SPH particle neighbours
       us(i)=zero;vs(i)=zero;ws(i)=zero
       do k=1,grid_neighbourcount(i)
          j=grid_neighbourlist(k,i)     !! Index of SPH particle neighbour       
       
          !! Find x_{i} (grid) - x_{j} (particle)
          x = rp(j,1) - xs(i)
          y = rp(j,2) - ys(i)
          z = rp(j,3) - zs(i)                    
          
          !! Scale x,y,z for hermite
          rad = sqrt(x*x + y*y)/h0;qq=rad
          ff1 = Wab(qq)
          x=x/h0;y=y/h0;z=z/h0             
                            
          !! Build vector of Basis functions
          !Re-scale x,y,z
          x=x*oosqrt2;y=y*oosqrt2;z=z*oosqrt2          
          wvec(1) = ff1*one
          wvec(2) = ff1*Hermite1(x)*oosqrt2
          wvec(3) = ff1*Hermite1(y)*oosqrt2
          wvec(4) = ff1*Hermite1(z)*oosqrt2                    
          wvec(5) = ff1*Hermite2(x)*half
          wvec(6) = ff1*Hermite1(x)*Hermite1(y)*half
          wvec(7) = ff1*Hermite2(y)*half
          wvec(8) = ff1*Hermite1(y)*Hermite1(z)*half
          wvec(9) = ff1*Hermite2(z)*half
          wvec(10)= ff1*Hermite1(z)*Hermite1(x)*half
          wvec(11)= ff1*Hermite3(x)*half*oosqrt2
          wvec(12)= ff1*Hermite2(x)*Hermite1(y)*half*oosqrt2
          wvec(13)= ff1*Hermite2(y)*Hermite1(x)*half*oosqrt2
          wvec(14)= ff1*Hermite3(y)*half*oosqrt2
          wvec(15)= ff1*Hermite2(y)*Hermite1(z)*half*oosqrt2
          wvec(16)= ff1*Hermite2(z)*Hermite1(y)*half*oosqrt2
          wvec(17)= ff1*Hermite3(z)*half*oosqrt2
          wvec(18)= ff1*Hermite2(z)*Hermite1(x)*half*oosqrt2
          wvec(19)= ff1*Hermite2(x)*Hermite1(z)*half*oosqrt2
          wvec(20)= ff1*Hermite1(x)*Hermite1(y)*Hermite1(z)*half*oosqrt2
          
          !! Dot-product with Psi to get weights.
          wij = dot_product(psivec,wvec)
          
          !! Evaluate properties at grid point i
          us(i) = us(i) + up(j,1)*wij          
          vs(i) = vs(i) + up(j,2)*wij
          ws(i) = ws(i) + up(j,3)*wij              
       
       end do
    
    
    end do
    !$omp end parallel do



    !! Write to output
    do kk=1,ngrid   !! Loop over all points in slice
!       us(kk)=us(kk)/vols(kk);vs(kk)=vs(kk)/vols(kk);ws(kk)=ws(kk)/vols(kk) !! Normalise by volume sum
       write(nf+1,500) xs(kk),ys(kk),zs(kk),us(kk),vs(kk),ws(kk)
    500 FORMAT(' ',6(E11.5,1X))       
    end do
  
    deallocate(xs,ys,us,vs,ws,zs,vols)
    deallocate(grid_neighbourcount,grid_neighbourlist)
    deallocate(amat,xvec,wvec,psivec)

    close(nf+1)
#endif
    return
  end subroutine grid_output_u
!! ------------------------------------------------------------------------------------------------
  subroutine surf_output
#if dim3
    integer idp,nf,i,n_out,i_out,nsgrid
    character(70) fnm1
    integer :: ii,jj,nxs,nys,kk,j,i0,j0,nss,ll,nzs,l0
    real(rkind),dimension(:),allocatable :: xs,ys,zs,srf,alp
    real(rkind),dimension(dims) :: rij
    real(rkind) :: rad,qq,wtmp,tta,ttb,kcoeff,dxgrid
  

    ! File naming stuff
    nf = 22
    idp = n_dump

    ! Create the name of the output file
    if( idp < 10 ) then 
       write(fnm1,'(A22,I1)') './data_directory/WSURF', idp
    else if( idp < 100 ) then 
       write(fnm1,'(A22,I2)') './data_directory/WSURF', idp
    else if( idp < 1000 ) then
       write(fnm1,'(A22,I3)') './data_directory/WSURF', idp
    else
       write(fnm1,'(A22,I4)') './data_directory/WSURF', idp
    end if

    ! Open the file
    open(nf+1,file=fnm1,status='unknown')
  
    !! Make slice grid
    ymax = yb_max

    dxgrid = 0.8*dx

    nxs = ceiling((xb_max-xb_min)/dxgrid)
    nys = ceiling((ymax-yb_min)/dxgrid)
    nzs = ceiling((zb_max-zb_min)/dxgrid)
    nsgrid=nxs*nys*nzs

    dxgrid = (xb_max-xb_min)/dble(nxs)
  
    !! Neighbour search size
    nss = ceiling(sup_size/dxgrid)
   
    allocate(xs(nsgrid),ys(nsgrid),zs(nsgrid))
    kk=0
    do ii=1,nxs
       do jj=1,nys
          do ll=1,nzs
             kk=kk+1
             xs(kk) = xb_min + dble(ii-1)*(xb_max-xb_min)/dble(nxs-1)
             ys(kk) = yb_min + dble(jj-1)*(ymax-yb_min)/dble(nys-1)        
             zs(kk) = zb_min + dble(ll-1)*(zb_max-zb_min)/dble(nzs-1)
          end do
       end do
    end do
  
    !! Loop over all particles, then find a smallish range of slice-nodes to loop over
    allocate(srf(nsgrid),alp(nsgrid))
    srf=zero;alp=zero
    !$omp parallel do private(i0,j0,l0,kcoeff,ll,qq,ii,jj,kk,rij,rad,wtmp) reduction(+:srf,alp)
    do i=1,n_par_fwm
       kcoeff = dv/vol(i)
       i0 = floor(nxs*(rp(i,1)-xb_min)/(xb_max-xb_min));i0=min(i0,nxs);i0=max(i0,1)
       j0 = floor(nys*(rp(i,2)-yb_min)/(ymax-yb_min));j0=min(j0,nys);j0=max(j0,1)
       l0 = floor(nzs*(rp(i,3)-zb_min)/(zb_max-zb_min));l0=min(l0,nzs);l0=max(l0,1)
       do ii=max(i0-nss,1),min(i0+nss,nxs)
          do jj=max(j0-nss,1),min(j0+nss,nys)
             do ll=max(l0-nss,1),min(l0+nss,nzs)
                kk = nzs*nys*(ii-1) + nzs*(jj-1) + ll !! Linear index..
                rij = (/xs(kk),ys(kk),zs(kk)/) - rp(i,:)
                rad = sqrt(dot_product(rij,rij));qq=rad/h(i)
                wtmp = wab(qq)*kcoeff        

                srf(kk) = srf(kk) + vol(i)*wtmp
                alp(kk) = alp(kk) + vol(i)*wtmp*a_out(i)             
             end do
          end do
       end do   
    end do
    !$omp end parallel do

    !! Write to output
    nsurf = 0
    do kk=1,nsgrid   !! Loop over all points in slice
      if(srf(kk)>5.0d-2.and.srf(kk)<=9.5d-1) then
           nsurf = nsurf+1
          write(nf+1,600) xs(kk),ys(kk),zs(kk),srf(kk),alp(kk)
       end if
    600 FORMAT(' ',5(E11.5,1X))       
    end do
  
    deallocate(xs,ys,srf,zs,alp)

    close(nf+1)
#endif
    return
  end subroutine surf_output
!! ------------------------------------------------------------------------------------------------  
  subroutine global_energy_stats
     !! Calculate and output total energy and related things
     integer(ikind) :: i
     real(rkind) :: tmp,wtmp
     real(rkind),dimension(dims) :: rij
     
     !! Total kinetic energy
     tmp=zero;wtmp = zero
     !$omp parallel do reduction(+:tmp)
     do i=1,n_par_fw
        if(conc(i).gt.0.75) then
           tmp = tmp + a0(i)*dot_product(up(i,:),up(i,:))*vol(i)/(4*pi*pi)
        end if
     end do    
     !$omp end parallel do
     tmp = half*tmp!/wtmp
     write(196,*) time,tmp
     flush(196)
     
     !! Mean velocity
     rij=zero;wtmp = zero
     !$omp parallel do reduction(+:rij,wtmp)
     do i=1,n_par_fw
        rij(:) = rij(:) + up(i,:)*vol(i)
        wtmp = wtmp + vol(i)
     end do    
     !$omp end parallel do
     rij = rij/wtmp
     write(197,*) time,rij
     flush(197)
     
     !! Graviational potential energy (for Breaking wave case)
     wtmp = zero
     !$omp parallel do reduction(+:wtmp) 
     do i=1,n_par_fw
        if(conc(i).gt.0.75) then     
           wtmp = wtmp + a0(i)*(rp(i,2)+(yb_max-yb_min)*0.5d0)*vol(i)
        end if
     end do    
     !$omp end parallel do
     write(195,*) time,wtmp
     flush(195)

     !! Surface energy of all bubbles       
#ifdef bubbles
     tmp = zero
     !$omp parallel do reduction(+:tmp)
     do i=1,n_bub
        tmp = tmp + four*pi*radb(i)**two
     end do
     !$omp end parallel do
     tmp = tmp/We
#endif        
     write(194,*) time,tmp
     flush(194)    
     
     return 
  end subroutine global_energy_stats
!! ------------------------------------------------------------------------------------------------
  subroutine plume_slice_data
     !! Output lots of data about bubble motion at a specific depth, for each bubble 
     !! as it passes that depth
     real(rkind) :: yslice,tmp,tmp_ste
     integer(ikind) :: i     

#ifdef bubbles
     !! Write lots of data about bubble motion at specific depth, for each bubble
     yslice = 0.5545d0       !D=0.081,z/D=5.5, z=0.4455, 1-z
     do i=1,n_bub
        tmp = rb(i,2)-yslice
        tmp_ste = rb(i,2)-yslice + ub(i,2)*dt
        tmp = tmp*tmp_ste !! Flag should be -ve only when a bubble is about to cross yslice
        if(tmp<=zero) then
           tmp_ste = sqrt(rb(i,1)**2 + rb(i,3)**2)
           !! time, ab/ah,x,z,ub*3,ul*3,ubturb*3
           write(173,*) time, radb(i)/Rhinze,rb(i,1),rb(i,3),ub(i,:),u_l2b(i,:),ubturb(i,:)
           flush(173)
        end if
     end do
#endif
     
     return
  end subroutine plume_slice_data     
!! ------------------------------------------------------------------------------------------------
  subroutine bubble_entrainment_stats
     !! Output number and volume of bubbles (since impact), and whether sub- or super- Hinze scale
     real(rkind) :: tmp
     integer(ikind) :: i,counta,countb
  
#ifdef bubbles
     !! Number of bubbles and total volume of air...  
     if(t_impact<zero.and.n_bub/=0) then
        t_impact = time           
     end if

     tmp = zero;counta=0
     !$omp parallel do reduction(+:tmp,counta)
     do i=1,n_bub
        if(.not.b_inbin(i))then
           tmp = tmp + (four/three)*pi*radb(i)**three
         
           !! Count super-Hinze bubbles
           if(radb(i)>Rhinze)then
              counta = counta + 1
           end if
        endif
     end do
     !$omp end parallel do
     tmp = tmp/(dble(n_par_fw)*dv)
        
     countb = n_bub-b_nfree - counta
        
     if(t_impact>=zero) then
        !! time, Vbubbles/Vliquid, total #bub, superHinze #, sub-Hinze # 
        write(172,*) time-t_impact,tmp,dble(n_bub-b_nfree),dble(counta),dble(countb)
        flush(172)
     endif
#endif         
     return
  end subroutine bubble_entrainment_stats
!! ------------------------------------------------------------------------------------------------
  subroutine channel_flow_statistics
     !! Calculate mean and fluctuating velocity for channel flow
     integer(ikind) :: i,j,k,nf
     real(rkind) :: qq,rad,wtmp
     real(rkind),dimension(dims) :: rij
     real(rkind) :: kcoeff
     real(rkind),dimension(100,dims) :: u_bar,u_bar_now
     real(rkind),dimension(100) :: yj,c_bar
     character(70) fnm1
             
     !! Mean flow over channel at time-instant...
     u_bar_now = zero
     c_bar = zero
     do j=1,100
        yj(j) = dble(j-1)*two/99.0d0
     end do
     nf = ceiling(two*h0/(yj(2)-yj(1)))+2
     !$omp parallel do private(k,j,qq,wtmp,rij,rad,kcoeff) reduction(+:c_bar,u_bar_now)
     do i=1,n_par_fwm   !! Loop over all particles
        kcoeff = dv/vol(i)
#if dim3           
        if(abs(rp(i,1))<=two*h(i).and.abs(rp(i,3))<=two*h(i)) then
#else
        if(abs(rp(i,1))<=two*h(i)) then              
#endif              
           k = 1+floor(99.0d0*rp(i,2)/two) !! find the index of point just below particle i
           k=max(1,k);k=min(100,k)
           do j = max(1,k-nf),min(k+nf,100)   !! Loop over 10 nearest nodes                
              rij = rp(i,:)-(/zero,yj(j),zero/)
              rad = sqrt(dot_product(rij,rij));qq = rad/h(i)
              wtmp = kcoeff*wab(qq)*vol(j)
              u_bar_now(j,:) = u_bar_now(j,:) + wtmp*up(i,:)
              c_bar(j) = c_bar(j) + wtmp
           end do           
        end if
     end do
     !$omp end parallel do
     do j=1,100
        u_bar_now(j,:)=u_bar_now(j,:)/c_bar(j)
     end do         
           
     ! Create the name of the output file
     n_dump = n_dump + 1;nf=99 
     if( n_dump < 10 ) then 
        write(fnm1,'(A4,I1)') 'mean', n_dump
     else if( n_dump < 100 ) then 
        write(fnm1,'(A4,I2)') 'mean', n_dump
     else if( n_dump < 1000 ) then
        write(fnm1,'(A4,I3)') 'mean', n_dump
     else
        write(fnm1,'(A4,I4)') 'mean', n_dump
     end if

     ! Open the file
     open(nf+1,file=fnm1,status='unknown')
     do i=1,100
        write(nf+1,*) yj(i),u_bar_now(i,:)
     end do
     close(nf+1)
       
  
     return
  end subroutine channel_flow_statistics
!! ------------------------------------------------------------------------------------------------
  subroutine fraga_analysis
     !! Evaluate velocity fields at slices at a couple of different depths for the Fraga plume case
     integer(ikind) :: i,j,k,nf
     real(rkind) :: x_loc1,qq,rad,wtmp
     real(rkind),dimension(dims) :: rij
     real(rkind) :: kcoeff
     character(70) fnm1
     real(rkind) :: yslice
     real(rkind),dimension(:),allocatable :: xj,zj
     integer(ikind) :: nn,mm,kx,kz,jx,jz
     real(rkind),dimension(:,:,:),allocatable :: uslice
     real(rkind),dimension(:,:),allocatable :: cslice       
  
  
#ifdef dim3        

     !! ---------------------------------------------------------------------------------
     !! Depth (or hieght above base) of slice
     yslice =0.45325d0       !D=0.081,z/D=6.75, z=0.54675, 1-z      
           
     !! Allocate and build xj,zj arrays (size nn*nn)
     mm = 50;nn = mm*2 + 1
     allocate(xj(nn),zj(nn))                 
     do j=1,nn
        xj(j) = dble(j-mm-1)*half/dble(mm)
        zj(j) = dble(j-mm-1)*half/dble(mm)              
     end do
           
     !! allocate velocity in slice
     allocate(uslice(nn,nn,dims),cslice(nn,nn))
     uslice = zero;cslice = zero

     !! How many grid nodes in a stencil...
     nf = ceiling(two*h0/(xj(2)-xj(1)))+2
           
     !! Loop over all particles
     !$omp parallel do private(kx,kz,jx,jz,rij,rad,wtmp,kcoeff) reduction(+:uslice,cslice)
     do i=1,n_par_fw !! Omit mirror particles, as we're interested in plume not boundary stuff   

        kcoeff = dv/vol(i)
              
        !! Only do for particles within 2h of yslice
        if(abs(rp(i,2))-yslice<=two*h(i))then   
              
           !! Indentify indices of a nearby grid node
           kx = floor(one + mm + two*rp(i,1)*dble(mm))
           kz = floor(one + mm + two*rp(i,3)*dble(mm))
                 
           !! Limit indices...
           kx=max(1,kx);kx=min(nn,kx);kz=max(1,kz);kz=min(nn,kz)
                 
           !! Loop over nearby nodes...
           do jx = max(1,kx-nf),min(kx+nf,nn)   
              do jz = max(1,kz-nf),min(kz+nf,nn)
                               
                 !! slice-grid-node particle distance and kernel                               
                 rij = rp(i,:)-(/xj(jx),yslice,zj(jz)/)
                 rad = sqrt(dot_product(rij,rij));qq = rad/h(i)
                 wtmp = kcoeff*wab(qq)*vol(j)
                    
                 !! Contribution to slice velocity and kernel sum
                 uslice(jx,jz,:) = uslice(jx,jz,:) + wtmp*up(i,:)
                 cslice(jx,jz) = cslice(jx,jz) + wtmp
              end do
           end do           
        end if
     end do
     !$omp end parallel do

     !! Normalise by kernel sum
     do jx=1,nn
        do jz=1,nn
           uslice(jx,jz,:) = uslice(jx,jz,:)/max(cslice(jx,jz),vsmall)  !! avoid NaN if zero conc...
        end do
     end do         
           
     ! Create the name of the output file
     n_dump = n_dump + 1;nf=99 
     if( n_dump < 10 ) then 
        write(fnm1,'(A4,I1)') 'u_r0', n_dump
     else if( n_dump < 100 ) then 
        write(fnm1,'(A4,I2)') 'u_r0', n_dump
     else if( n_dump < 1000 ) then
        write(fnm1,'(A4,I3)') 'u_r0', n_dump
     else
        write(fnm1,'(A4,I4)') 'u_r0', n_dump
     end if

     ! Open the file
     open(nf+1,file=fnm1,status='unknown')
     do jx=1,nn
        do jz=1,nn
           wtmp = sqrt(xj(jx)**two + zj(jz)**two)
           write(nf+1,*) wtmp,uslice(jx,jz,:)
        end do
     end do
     close(nf+1)
     deallocate(xj,zj,uslice,cslice)                 
        
     !! ---------------------------------------------------------------------------------
     !! Depth (or hieght above base) of slice
     yslice = 0.5545d0       !D=0.081,z/D=5.5, z=0.4455, 1-z
         
     !! Allocate and build xj,zj arrays (size nn*nn)
     mm = 50;nn = mm*2 + 1
     allocate(xj(nn),zj(nn))                 
     do j=1,nn
        xj(j) = dble(j-mm-1)*half/dble(mm)
        zj(j) = dble(j-mm-1)*half/dble(mm)              
     end do
           
     !! allocate velocity in slice
     allocate(uslice(nn,nn,dims),cslice(nn,nn))
     uslice = zero;cslice = zero

     !! How many grid nodes in a stencil...
     nf = ceiling(two*h0/(xj(2)-xj(1)))+2
           
     !! Loop over all particles
     !$omp parallel do private(kx,kz,jx,jz,rij,rad,wtmp,kcoeff) reduction(+:uslice,cslice)
     do i=1,n_par_fw !! Omit mirror particles, as we're interested in plume not boundary stuff   

        kcoeff = dv/vol(i)
              
        !! Only do for particles within 2h of yslice
        if(abs(rp(i,2))-yslice<=two*h(i))then   
             
           !! Indentify indices of a nearby grid node
           kx = floor(one + mm + two*rp(i,1)*dble(mm))
           kz = floor(one + mm + two*rp(i,3)*dble(mm))
                 
           !! Limit indices...
           kx=max(1,kx);kx=min(nn,kx);kz=max(1,kz);kz=min(nn,kz)
                
           !! Loop over nearby nodes...
           do jx = max(1,kx-nf),min(kx+nf,nn)   
              do jz = max(1,kz-nf),min(kz+nf,nn)
                            
                 !! slice-grid-node particle distance and kernel                               
                 rij = rp(i,:)-(/xj(jx),yslice,zj(jz)/)
                 rad = sqrt(dot_product(rij,rij));qq = rad/h(i)
                 wtmp = kcoeff*wab(qq)*vol(j)
                    
                 !! Contribution to slice velocity and kernel sum
                 uslice(jx,jz,:) = uslice(jx,jz,:) + wtmp*up(i,:)
                 cslice(jx,jz) = cslice(jx,jz) + wtmp
              end do
           end do           
        end if
     end do
     !$omp end parallel do

     !! Normalise by kernel sum
     do jx=1,nn
        do jz=1,nn
           uslice(jx,jz,:) = uslice(jx,jz,:)/max(cslice(jx,jz),vsmall)  !! avoid NaN if zero conc...
        end do
     end do         
           
     ! Create the name of the output file
     nf=99 
     if( n_dump < 10 ) then 
        write(fnm1,'(A4,I1)') 'u_r1', n_dump
     else if( n_dump < 100 ) then 
        write(fnm1,'(A4,I2)') 'u_r1', n_dump
     else if( n_dump < 1000 ) then
        write(fnm1,'(A4,I3)') 'u_r1', n_dump
     else
        write(fnm1,'(A4,I4)') 'u_r1', n_dump
     end if

     ! Open the file
     open(nf+1,file=fnm1,status='unknown')
     do jx=1,nn
        do jz=1,nn
           wtmp = sqrt(xj(jx)**two + zj(jz)**two)
           write(nf+1,*) wtmp,uslice(jx,jz,:)
        end do
     end do
     close(nf+1)
     deallocate(xj,zj,uslice,cslice)                 
#endif  
  
     return
  end subroutine fraga_analysis
!! ------------------------------------------------------------------------------------------------
  subroutine taylor_green_error
     !! L2 of velocity for Taylor Green flow
     real(rkind) :: tmpa,tmpb,tmpc,uexpon,U_c,ue,ve,umage
     integer(ikind) :: i
     
     tmpa=zero;tmpb=zero;tmpc=zero
     uexpon=exp(-two*(time)/Re)     
     !$omp parallel do private(ue,ve,umage,U_c) reduction(+:tmpa,tmpb,tmpc)
     do i=1,n_par_fw
        ue = uexpon*sin(rp(i,1))*cos(rp(i,2))
        ve = -uexpon*cos(rp(i,1))*sin(rp(i,2))
        umage = sqrt(ue*ue + ve*ve)               !! Exact
        U_c = sqrt(dot_product(up(i,:),up(i,:)))  !! Numerical
        tmpa = tmpa + (U_c-umage)**two  !! L2 of error
        tmpb = tmpb + umage**two  !! L2 of exact
        tmpc = tmpc + U_c**two    !! L2 of numerical
     end do
     !$omp end parallel do
     tmpa = sqrt(tmpa/n_par_fw)/sqrt(tmpb/n_par_fw)  !! L2 of error
     tmpc = sqrt(tmpc/n_par_fw)   !! L2 of velocity
     tmpb = sqrt(tmpb/n_par_fw)
     write(192,*) time,tmpa   
     flush(192)
    
     return
  end subroutine taylor_green_error
!! ------------------------------------------------------------------------------------------------  
  subroutine fs_elev
     !! Find max surface elevation at x_target
     integer(ikind) :: i,j
     real(rkind),dimension(10) :: fs_max,x_loc
     real(rkind) :: fs_mean
     integer(ikind) :: fs_sum
     
     !! set the locations to measure at: 0, 0.25, 0.5, 0.75, 1, 1.25 etc *x_target
     do i=1,10
        x_loc(i) = dble(i-1)*0.25*x_target
     end do
     
     fs_max = -1.0d5 !! Initialise big negative number
     fs_mean = zero
     !$omp parallel do private(j) reduction(max:fs_max) reduction(+:fs_mean,fs_sum)
     do i=1,n_par_fw
        if(n_surf(i).eq.1) then  !! Only search through free-surface particles
           fs_mean = fs_mean + rp(i,2) !! Augment mean
           fs_sum = fs_sum + 1
           do j=1,10
              if(abs(rp(i,1)-x_loc(j)).le.two*dx)then   !! Only particles in vertical slice near x_target
                 fs_max(j) = max(fs_max(j),rp(i,2))
              end if        
           end do
        end if
     end do
     !$omp end parallel do
     !! Evaluate mean free-surface height
     fs_mean = fs_mean/dble(fs_sum)
     
     !! Remove mean from max
     fs_max = fs_max - fs_mean
     
     !! Output
     write(175,*) time,fs_max
     flush(175)
    
     return
  end subroutine fs_elev 
!! ------------------------------------------------------------------------------------------------
  subroutine setup_stats_files

     !! Total kinetic energy  
     open(unit=196,file='./data_directory/statistics/total_ke')
     
     !! Mean dissipation rate - calculated from enstrophy
     open(unit=198,file='./data_directory/statistics/mean_diss_rate')
     
     !! Mean dissipation rate - calculated from viscous term
     open(unit=199,file='./data_directory/statistics/mean_diss_rate2')

     !! Mean velocity
     open(unit=197,file='./data_directory/statistics/mean_velocity')    

     !! Total surface energy of all bubbles
     open(unit=194,file='./data_directory/statistics/total_surface_energy')     
     
     !! Gravitational potential energy
     open(unit=195,file='./data_directory/statistics/total_gpe')
     
     !! L2 Error relative to analytic Taylor-Green solution
     open(unit=192,file='./data_directory/statistics/taylor_green_error_l2')     

     !! parent-child volume ratios for breakup events
     open(unit=181,file='./data_directory/statistics/bubble_parent_child_vol_ratios')         

     !! Bubble breakup rate
     open(unit=188,file='./data_directory/statistics/bubble_breakup_rate')         

     !! Bubble death rate
     open(unit=189,file='./data_directory/statistics/bubble_death_rate')      
     
     !! Volume of bubbles and number of bubbles
     open(unit=186,file='./data_directory/statistics/bubble_count')         
     
     !! Time since impact, vol bubbles, number of bubbles, number sub and super Hinze
     open(unit=172,file='./data_directory/statistics/bubble_count2')         
     
     !! time, size,position,velocity,liquid velocity,srs fluctuations of bubbles as they pass
     !! a specific y-position. Used in plume
     open(unit=173,file='./data_directory/statistics/bubble_velocity_at_depth')     
     
     !! time, wavemaker position, wavemaker velocity
     open(unit=174,file='./data_directory/statistics/wavemaker_X_U')         
     
     !! time, free surface elevation
     open(unit=175,file='./data_directory/statistics/free_surface_elev')
     
   
  
     return
  end subroutine setup_stats_files
!! ------------------------------------------------------------------------------------------------
!! ================================================================================================
  !! Univariate Hermite Polynomials (Physicists kind)
  function Hermite1(q) result(Hres)
     real(rkind),intent(in) :: q
     real(rkind) :: Hres
     Hres = 2.0d0*q
  end function Hermite1
  function Hermite2(q) result(Hres)
     real(rkind),intent(in) :: q
     real(rkind) :: Hres
     Hres = 4.0d0*q*q - 2.0d0
  end function Hermite2
  function Hermite3(q) result(Hres)
     real(rkind),intent(in) :: q
     real(rkind) :: Hres
     Hres = 8.0d0*q*q*q - 12.0d0*q
  end function Hermite3
  function Hermite4(q) result(Hres)
     real(rkind),intent(in) :: q
     real(rkind) :: Hres
     Hres = 16.0d0*q*q*q*q - 48.0d0*q*q + 12.0d0
  end function Hermite4
  
!! ------------------------------------------------------------------------------------------------
end module output
