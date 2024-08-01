module neighbour_finding
  !! Build neighbour lists for SPH and bubbles, by placing particles and bubbles
  !! in cells and then searching over neighbouring cells.
  use common_parameter
  use common_2d
  use sphtools
  use omp_lib
  implicit none

  private
  public :: particle_allocation

  integer(ikind), dimension(:), allocatable ::ist, ip, nc
  integer(ikind), dimension(:), allocatable :: cellpart
#ifdef bubbles
  integer(ikind), dimension(:), allocatable ::ist_b, ip_b, nc_b
  integer(ikind), dimension(:), allocatable :: cellpart_b
#endif
   
!! Note: This neighbour finding module is a bit simplified by the fact that we create mirror
!! particles for periodic boundaries. The additional cost & memory of creating mirrors for periodic
!! boundaries seems minimal, and it a) simplifies this module, and b) avoids having to check & 
!! modify rij for every kernel evaluation.  

contains
!! ------------------------------------------------------------------------------------------------
  subroutine particle_allocation
    !! This subroutine does some allocation, then calls routines to put particles in
    !! cells, reorder them, and look through cells to build neighbour lists..
    integer nct_temp
    

    !! allocate linking lists
    allocate(ij_count(n_par_fw),ij_link(nplink,n_par_fw))   
#ifdef bubbles
    allocate(ijb_count(n_par_fwm))          !! Note that mirror particles also hold a count of bubble neighbours
    allocate(ijb_firstindex(n_par_fwm+1))
#endif
    
    !! allocate cell lists
    nct_temp=nct+ncx+ncy+3
    allocate(ist(nct_temp+1),nc(nct_temp),ip(n_par_fw+nmirror_esti));nc=0
#ifdef bubbles
    allocate(ist_b(nct_temp+1),nc_b(nct_temp),ip_b(n_bub_m));nc_b=0
#endif    

    !! put particles in different cells, including mirror particles	
#if dim3
    call divide3d
#else    
    call divide
#endif    

    !! Do the same for bubbles
#ifdef bubbles
#if dim3
    call divide3d_b     
#else    
    call divide_b
#endif    
#endif

    !! Search through the neighbouring cells and build neighbour lists
    call neighboring_list_parallel
#ifdef bubbles
    call neighboring_list_parallel_b
#endif

    !! Free up memory again
    deallocate(ist,nc,ip,cellpart)
#ifdef bubbles
    deallocate(ist_b,nc_b,ip_b,cellpart_b)
#endif    

  end subroutine particle_allocation
!! ------------------------------------------------------------------------------------------------
  subroutine divide
    !! This subroutine divides the particles amongst the cells, and creates
    !! a list of particles ordered by their cell
    real(rkind) :: ddx,ddy
    integer(ikind) :: nc_ii_max, k, icell, jcell, ii,nini,nend
    nini=1;nend=n_par_fwm
    nc=0_ikind  ! nc(i,mkind)=0

    nc_ii_max = 0
    allocate(cellpart(nend-nini+1))

    !$omp parallel do private(ddx,ddy,icell,jcell,ii) reduction(+:nc) reduction(max:nc_ii_max)
    do k=nini, nend    ! loop over specified range of particles
       if(.not.inbin(k)) then  ! don't try and find a cell for particles in the bin.
          ddx = rp(k,1) - xmint    ! find x,y of particle if origin at BL corner of domain
          ddy = rp(k,2) - ymint

          icell = int( ddx *uno_cell_size) + 1       !! Cell label in x
          jcell = int( ddy *uno_cell_size) + 1       !! Cell label in y
          ii = icell + (jcell - 1)*ncx        !! ii is linear cell position in the matrix of cells
          if(ii > nct .or. ii < 1) then
             print *, "particles are out of domain"
             print *, "ii=", ii, "nct_t=", nct
             print *, "nct=", nct
             print *, "particle in trouble is ",k,ddx,ddy
             stop
          endif
          nc(ii) = nc(ii)+1     ! add 1 to the number of particles in this cell
          !! If the number of particles in a cell is greater than
          !! some upper limit then stop
          nc_ii_max = nc(ii)  ! ensure nc_ii_max is max(nc(ii) for all ii
          cellpart(k)=ii             !! the index of the cell in which particle i is
       end if
    end do
    !$omp end parallel do

    !! loop over all cells
    ist(1)=1_ikind
    do ii=1, nct
       ist(ii+1)=ist(ii)+nc(ii)   ! ist(ii) is the starting index for each cell
       nc(ii)=0_ikind             ! erase the nc in each cell here 
    end do

    !! index look-up: ip(j)=k where j is the cell-ordered index and k is the original index
    do k=nini,nend
       if(inbin(k))cycle  ! don't do this for particles in the bin
       ip(ist(cellpart(k))+nc(cellpart(k)))=k   
        nc(cellpart(k))=nc(cellpart(k))+1
    end do 

  end subroutine divide
!! ------------------------------------------------------------------------------------------------
  subroutine  neighboring_list_parallel
    !! This subroutine loops through every particle i, finds the cell ic containing
    !! particle i, then goes through every cell jc adjacent to ic to check for neighbours  
    !! The list of neighbouring cells is built at the start of the simulation in the cell_gridding
    !! module.
    integer(ikind) :: i,ic
    integer(ikind) :: jc,kc
    real(rkind) :: ss_local2

    ij_count=0
    ij_link=0 !neighbour list arrays
        
    !! Loop over all particles (in parallel)
    !$omp parallel do private(ic,jc,kc,ss_local2)
    do i=1,n_par_fw 
       if(.not.inbin(i))then
          !! Which cell are we in?
          ic = cellpart(i)  
          
          !! What is the local support size
          ss_local2 = ss2*h(i)*h(i)       

          !! Loop over all cells neighbouring cell ic
          do kc=1,ic_count(ic)
             jc = ic_link(ic,kc)

             !! Add any neighbours in cell jc to neighbour list of particle i
             call neighbour_particle_cell(i,rp(i,:),jc,ss_local2)
   
          end do
       end if   
    end do
    !$omp end parallel do

  end subroutine neighboring_list_parallel
!! ------------------------------------------------------------------------------------------------
  subroutine neighbour_particle_cell(ii,ri,jc,ss_local2)
    !! looks through particles in cell jc and checks whether they are neighbours of 
    !! particle ii which has position (vector) ri.
    integer(ikind),intent(in) :: ii,jc !! jc is cell
    real(rkind),dimension(:),intent(in) :: ri
    real(rkind),intent(in) :: ss_local2
    real(rkind) :: rr2tmp
    integer(ikind) :: j,jj,is,ie
    real(rkind),dimension(dims) :: rij
    
    if(nc(jc)/=0)then !! if the cell isn't empty
       is =ist(jc);ie=ist(jc+1)-1  !! Start and end indices of particles in cell jc
       do jj=is,ie       !! Loop over all particles in cell jc
          j=ip(jj)       !! j is regular index, jj is cell-ordered index

          rij(:) = ri(:)-rp(j,:)
          rr2tmp=dot_product(rij,rij)  !! Distance squared
!          if(j/=ii .and. rr2tmp<=sup_size_2)then ! Don't include self in neighbour list
          if(rr2tmp <= ss_local2)then          ! Include self (faster)    
             ij_count(ii)=ij_count(ii)+1;              !! Increment count
             ij_link(ij_count(ii),ii)=j                !! add to list  
          endif          
       end do
    end if
        
    return
  end subroutine neighbour_particle_cell
!! ------------------------------------------------------------------------------------------------
  subroutine divide3d
    !! This subroutine divides the particles amongst the cells, and creates
    !! a list of particles ordered by their cell  
    real(rkind) :: ddx,ddy,ddz
    integer(ikind) :: nc_ii_max, k, icell, jcell, ii,kcell,nsheet,nini,nend
    nsheet = ncx*ncy
    nini=1;nend=n_par_fwm
    nc=0_ikind

    nc_ii_max = 0
    allocate(cellpart(nend-nini+1))

    !! This loop is parallelised, but I don't expect the parallelism to be particularly efficient,
    !! because we don't know which cell each particle is in a priori, and the particles are not
    !! ordered by cell yet...
    !$omp parallel do private(ddx,ddy,ddz,icell,jcell,kcell,ii) reduction(+:nc) reduction(max:nc_ii_max)
    do k=nini,nend    ! loop over specified range of particles
       if(.not.inbin(k)) then  ! don't try and find a cell for particles in the bin.
          ddx = rp(k,1) - xmint    ! find position of particle w.r.t. origin at BL corner of domain
          ddy = rp(k,2) - ymint
          ddz = rp(k,3) - zmint

          icell = int( ddx *uno_cell_size) + 1       !! Cell label in x
          jcell = int( ddy *uno_cell_size) + 1       !! Cell label in y
          kcell = int( ddz *uno_cell_size) + 1       !! Cell label in z
          ii = icell + (jcell - 1)*ncx + (kcell - 1)*nsheet     !! ii is linear cell position in the matrix of cells
          if(ii > nct .or. ii < 1) then
             print *, "particles are out of domain",icell,jcell,kcell
             print *, "ii=", ii, "nct_t=", nct
             print *, "particle in trouble is ",k,ddx,ddy,ddz
             write(6,*) "in xyz FoR: ",rp(k,:)
             if(k>n_par_fw) write(6,*) "parent is",irelation(k),rp(irelation(k),:),"inbin",inbin(irelation(k))             
             stop
          endif
          nc(ii) = nc(ii)+1     ! add 1 to the number of particles in this cell
          !! If the number of particles in a cell is greater than  some upper limit then stop
          nc_ii_max = nc(ii)  ! ensure nc_ii_max is max(nc(ii) for all ii
          cellpart(k)=ii             !! the index of the cell in which particle i is
       end if
    end do
    !$omp end parallel do

    !! loop over all cells to find starting index of particles for each cell
    ist(1)=1_ikind
    do ii=1, nct
       ist(ii+1)=ist(ii)+nc(ii)   ! ist(ii) is the starting index for each cell
       nc(ii)=0_ikind             ! erase the nc in each cell here 
    end do

    !! index look-up: ip(j)=k where j is the cell-ordered index and k is the original index
    do k=nini,nend
       if(inbin(k))cycle  ! don't do this for particles in the bin
       ip(ist(cellpart(k))+nc(cellpart(k)))=k   
       nc(cellpart(k))=nc(cellpart(k))+1  
    end do 

  end subroutine divide3d  
!! ------------------------------------------------------------------------------------------------
#ifdef bubbles
  subroutine divide_b 
    !! This subroutine divides the particles amongst the cells, and creates
    !! a list of particles ordered by their cell
    real(rkind) :: ddx,ddy
    integer(ikind) :: nc_ii_max, k, icell, jcell, ii,nini,nend
    nini=1;nend=n_bub_m
    nc_b=0_ikind  

    nc_ii_max = 0
    allocate(cellpart_b(nend-nini+1))

    !$omp parallel do private(ddx,ddy,icell,jcell,ii) reduction(+:nc_b) reduction(max:nc_ii_max)
    do k=nini, nend    ! loop over specified range of particles
       if(.not.b_inbin(k)) then
          ddx = rb(k,1) - xmint    ! find x,y of particle if origin at BL corner of domain
          ddy = rb(k,2) - ymint

          icell = int( ddx *uno_cell_size) + 1       !! Cell label in x
          jcell = int( ddy *uno_cell_size) + 1       !! Cell label in y
          ii = icell + (jcell - 1)*ncx        !! ii is linear cell position in the matrix of cells
          if(ii > nct .or. ii < 1) then
             print *, "bubbles are out of domain"
             print *, "ii=", ii, "nct_t=", nct
             print *, "nct=", nct
             print *, "bubble in trouble is ",k,ddx,ddy
             stop
          endif
          nc_b(ii) = nc_b(ii)+1     ! add 1 to the number of particles in this cell
          !! If the number of particles in a cell is greater than
          !! some upper limit then stop
          nc_ii_max = nc_b(ii)  ! ensure nc_ii_max is max(nc(ii) for all ii
          cellpart_b(k)=ii             !! the index of the cell in which particle i is
       end if
    end do
    !$omp end parallel do

    !! loop over all cells
    ist_b(1)=1_ikind
    do ii=1, nct
       ist_b(ii+1)=ist_b(ii)+nc_b(ii)   ! ist(ii) is the starting index for each cell
       nc_b(ii)=0_ikind             ! erase the nc in each cell here 
    end do

    !! index look-up: ip(j)=k where j is the cell-ordered index and k is the original index
    do k=nini,nend
       if(b_inbin(k))cycle  ! don't do this for bubbles in the bin    
       ip_b(ist_b(cellpart_b(k))+nc_b(cellpart_b(k)))=k   
       nc_b(cellpart_b(k))=nc_b(cellpart_b(k))+1
    end do 
    
  end subroutine divide_b  
!! ------------------------------------------------------------------------------------------------  
  subroutine divide3d_b
    !! This subroutine divides the particles amongst the cells, and creates
    !! a list of particles ordered by their cell  
    real(rkind) :: ddx,ddy,ddz
    integer(ikind) :: nc_ii_max, k, icell, jcell, ii,kcell,nsheet,nini,nend
    nsheet = ncx*ncy
    nini=1;nend=n_bub_m
    nc_b=0_ikind;nc_ii_max = 0
    allocate(cellpart_b(nend-nini+1))

    !! This loop is parallelised, but I don't expect the parallelism to be particularly efficient,
    !! because we don't know which cell each particle is in a priori, and the particles are not
    !! ordered by cell yet...
    !$omp parallel do private(ddx,ddy,ddz,icell,jcell,kcell,ii) reduction(+:nc_b) reduction(max:nc_ii_max)
    do k=nini,nend    ! loop over specified range of particles
       if(.not.b_inbin(k))then
          ddx = rb(k,1) - xmint    ! find position of particle w.r.t. origin at BL corner of domain
          ddy = rb(k,2) - ymint
          ddz = rb(k,3) - zmint

          icell = int( ddx *uno_cell_size) + 1       !! Cell label in x
          jcell = int( ddy *uno_cell_size) + 1       !! Cell label in y
          kcell = int( ddz *uno_cell_size) + 1       !! Cell label in z
          ii = icell + (jcell - 1)*ncx + (kcell - 1)*nsheet     !! ii is linear cell position in the matrix of cells
          if(ii > nct .or. ii < 1) then
             print *, "bubbles are out of domain"
             print *, "ii=", ii, "nct_t=", nct
             print *, "bubble in trouble is ",k,ddx,ddy,ddz
             print *, "bubble numbers: n_bub,n_bub_m",n_bub,n_bub_m
          endif
          nc_b(ii) = nc_b(ii)+1     ! add 1 to the number of particles in this cell
          !! If the number of particles in a cell is greater than  some upper limit then stop
          nc_ii_max = nc_b(ii)  ! ensure nc_ii_max is max(nc(ii) for all ii
          cellpart_b(k)=ii             !! the index of the cell in which particle i is
       end if       
    end do
    !$omp end parallel do

    !! loop over all cells to find starting index of particles for each cell
    ist_b(1)=1_ikind
    do ii=1, nct
       ist_b(ii+1)=ist_b(ii)+nc_b(ii)   ! ist(ii) is the starting index for each cell
       nc_b(ii)=0_ikind             ! erase the nc in each cell here 
    end do

    !! index look-up: ip(j)=k where j is the cell-ordered index and k is the original index
    do k=nini,nend
       if(b_inbin(k))cycle  ! don't do this for bubbles in the bin    
       ip_b(ist_b(cellpart_b(k))+nc_b(cellpart_b(k)))=k   
       nc_b(cellpart_b(k))=nc_b(cellpart_b(k))+1  
    end do 

  end subroutine divide3d_b    
!! ------------------------------------------------------------------------------------------------
  subroutine  neighboring_list_parallel_b
    !! This subroutine loops through every particle i, finds the cell ic containing
    !! particle i, then goes through every cell jc adjacent to ic to check for neighbours  
    !! The list of neighbouring cells is built at the start of the simulation in the cell_gridding
    !! module.
    integer(ikind) :: i,ic
    integer(ikind) :: jc,kc
    real(rkind) :: ss_local2
    integer(ikind) :: linearlink_length

    !! Loop over all particles, and estimate a max number of bubble neighbours
    ijb_count = 0
    !$omp parallel do private(ic,kc,jc)
    do i=1,n_par_fwm
       if(.not.inbin(i))then    
          !! Which cell are we in?
          ic = cellpart(i)
       
          !! Loop over all cells neighbouring cell ic
          do kc=1,ic_count(ic)
             jc = ic_link(ic,kc)
          
             !! Add the number of bubbles in each cell
             ijb_count(i) = ijb_count(i) + nc_b(jc)             
          end do
       end if
    end do
       
    !! Locate the starting index for this neighbour-list
    !! This is the only part not in parallel
    ijb_firstindex(1) = 1
    do i=1,n_par_fwm
       ijb_firstindex(i+1) = ijb_firstindex(i) + ijb_count(i)
    end do
    linearlink_length = ijb_firstindex(n_par_fwm+1)
    
    !! Allocate the linear link list, flag and reset the count
    allocate(ijb_linearlink(linearlink_length));ijb_linearlink=0
    ijb_count=0
        
    !! Loop over all particles (in parallel)
    !$omp parallel do private(ic,jc,kc,ss_local2)
    do i=1,n_par_fwm 
       if(.not.inbin(i))then
          !! Which cell are we in?
          ic = cellpart(i)     
          
          !! What is the local support size
          ss_local2 = h(i)*h(i)*ss2    

          !! Loop over all cells neighbouring cell ic
          do kc=1,ic_count(ic)
             jc = ic_link(ic,kc)
          
             !! Add any neighbours in cell jc to neighbour list of particle i
             call neighbour_bubble_cell(i,rp(i,:),jc,ss_local2)
   
          end do
       end if     
    end do
    !$omp end parallel do
    
  end subroutine neighboring_list_parallel_b
!! ------------------------------------------------------------------------------------------------
  subroutine neighbour_bubble_cell(ii,ri,jc,ss_local2)
    !! looks through particles in cell jc and checks whether they are neighbours of 
    !! particle ii which has position (vector) ri.
    integer(ikind),intent(in) :: ii,jc !! jc is cell
    real(rkind),dimension(:),intent(in) :: ri
    real(rkind),intent(in) :: ss_local2    
    real(rkind) :: rr2tmp
    integer(ikind) :: j,jj,is,ie
    real(rkind),dimension(dims) :: rij
    
    if(nc_b(jc)/=0)then !! if the cell isn't empty
       is =ist_b(jc);ie=ist_b(jc+1)-1  !! Start and end indices of particles in cell jc
       do jj=is,ie       !! Loop over all particles in cell jc
          j=ip_b(jj)       !! j is regular index, jj is cell-ordered index

          rij(:) = ri(:)-rb(j,:)
          rr2tmp=dot_product(rij,rij)  !! Distance squared
          if(rr2tmp <= ss_local2)then          ! Include self (faster)    
             ijb_count(ii)=ijb_count(ii)+1;              !! Increment count
             ijb_linearlink(ijb_firstindex(ii)-1+ijb_count(ii))=j !! Add to linear link list
          endif          
       end do
    end if
    return
  end subroutine neighbour_bubble_cell  
#endif  
!! ------------------------------------------------------------------------------------------------
end module neighbour_finding
