module cell_griddin
  !! The routines in this module are used to generate cells and linked lists of neighbour-cells
  !! The routines are called only at the start of a simulation, and only from the main program
  
  !! For bubbly flows, they *may* be called during a simulation, as cell sizes require adjusting
  !! to account for particle volume changes.
  
  !! Periodicity is NOT controlled in this module. 
  use kind_parameters
  use common_parameter
  use common_2d
  implicit none

  private
  public :: cell_generator, cells_third_dimension, neighbour_cells,neighbour_cells3d

contains
!! ------------------------------------------------------------------------------------------------
  subroutine cell_generator

    !! determine number of cells in x, y directions and total
    ncx = int( (xmaxt-xmint)*uno_cell_size ) + 1 
    ncy = int( (ymaxt-ymint)*uno_cell_size ) + 1

    nct=ncx*ncy ! Number of cells in a XY sheet

  end subroutine cell_generator
!! ------------------------------------------------------------------------------------------------
  subroutine cells_third_dimension
    !! Multiply the 2D sheet of cells by the appropriate number
    !! of cells in the third dimension...
    
    !! How many cells in the third dimension: ncz
    ncz = int( (zmaxt-zmint)*uno_cell_size) + 1
    
    !! How many cells in total: nct = nct*ncz
    nct = nct*ncz

    return
  end subroutine cells_third_dimension
!! ------------------------------------------------------------------------------------------------  
  subroutine neighbour_cells
    !! For each cell, build a list of the 9 neighbouring cells (including self)
    use common_2d
    integer(ikind) :: ic,icx,icy,jc,nsheet
   
    nsheet = ncx*ncy
   
    !! Allocate the cell-link-lists
    allocate(ic_count(nct),ic_link(nct,9))
 
    !! Loop over all cells
    icx=0;icy=1
    do ic = 1,nct
       icx = icx + 1
       if(icx>ncx) then
          icx = 1
          icy = icy + 1
       end if
     
       ! This cell
       ic_count(ic)=1;ic_link(ic,1)=ic
       
       if(icx+1<=ncx)then
          !! East cell
          jc = ic+1
          ic_count(ic)=ic_count(ic)+1;ic_link(ic,ic_count(ic))=jc
          if(icy+1<=ncy)then
             !! NE cell
             jc = ic + ncx + 1    
             ic_count(ic)=ic_count(ic)+1;ic_link(ic,ic_count(ic))=jc
          end if          
          if(icy-1>=1)then
             !! SE cell
             jc = ic - ncx + 1
             ic_count(ic)=ic_count(ic)+1;ic_link(ic,ic_count(ic))=jc
          end if
       end if       
       if(icx-1>=1)then
          !! West cell
          jc = ic - 1
          ic_count(ic)=ic_count(ic)+1;ic_link(ic,ic_count(ic))=jc
          if(icy+1<=ncy)then
             !! NW cell
             jc = ic + ncx - 1
          ic_count(ic)=ic_count(ic)+1;ic_link(ic,ic_count(ic))=jc
          end if          
          if(icy-1>=1)then
             !! SW cell
             jc = ic - ncx - 1
          ic_count(ic)=ic_count(ic)+1;ic_link(ic,ic_count(ic))=jc
          end if         
       end if       
       if(icy+1<=ncy) then
          !! North cell
          jc = ic + ncx
          ic_count(ic)=ic_count(ic)+1;ic_link(ic,ic_count(ic))=jc
       end if      
       if(icy-1>=1)then
          !! South cell
          jc = ic - ncx
          ic_count(ic)=ic_count(ic)+1;ic_link(ic,ic_count(ic))=jc
       end if          
    end do
   
    return
  end subroutine neighbour_cells 
!! ------------------------------------------------------------------------------------------------     
  subroutine neighbour_cells3d
    !! For each cell, build a list of the 27 neighbouring cells (including itself)
    use common_2d
    integer(ikind) :: ic,icx,icy,icz,ics,jc,nsheet
    
    !! Cells per sheet?
    nsheet = ncx*ncy
  
    !! Allocate the cell-link-lists
    allocate(ic_count(nct),ic_link(nct,27))
 
    !! Loop over all cells
    icx=0;icy=1;icz=1
    do ic = 1,nct
       icx = icx + 1
       if(icx>ncx) then
          icx = 1
          icy = icy + 1
       end if
       if(icy>ncy) then
          icy = 1
          icz = icz + 1
       end if  


       !!! This sheet
    
       ! This cell
       ic_count(ic)=1;ic_link(ic,1)=ic       
       if(icx+1<=ncx)then          !! East cell
          jc = ic+1;ic_count(ic)=ic_count(ic)+1;ic_link(ic,ic_count(ic))=jc
          if(icy+1<=ncy)then             !! NE cell
             jc = ic + ncx + 1;ic_count(ic)=ic_count(ic)+1;ic_link(ic,ic_count(ic))=jc
          end if          
          if(icy-1>=1)then             !! SE cell
             jc = ic - ncx + 1;ic_count(ic)=ic_count(ic)+1;ic_link(ic,ic_count(ic))=jc
          end if
       end if       
       if(icx-1>=1)then          !! West cell
          jc = ic - 1;ic_count(ic)=ic_count(ic)+1;ic_link(ic,ic_count(ic))=jc
          if(icy+1<=ncy)then             !! NW cell
             jc = ic + ncx - 1;ic_count(ic)=ic_count(ic)+1;ic_link(ic,ic_count(ic))=jc
          end if          
          if(icy-1>=1)then             !! SW cell
             jc = ic - ncx - 1;ic_count(ic)=ic_count(ic)+1;ic_link(ic,ic_count(ic))=jc
          end if         
       end if       
       if(icy+1<=ncy) then          !! North cell
          jc = ic + ncx;ic_count(ic)=ic_count(ic)+1;ic_link(ic,ic_count(ic))=jc
       end if      
       if(icy-1>=1)then          !! South cell
          jc = ic - ncx;ic_count(ic)=ic_count(ic)+1;ic_link(ic,ic_count(ic))=jc
       end if  

       !! Previous sheet
       if(icz-1>=1)then
          ics = ic - nsheet      
          jc = ics;ic_count(ic)=ic_count(ic)+1;ic_link(ic,ic_count(ic))=jc
          if(icx+1<=ncx)then             !! East cell
             jc = ics+1;ic_count(ic)=ic_count(ic)+1;ic_link(ic,ic_count(ic))=jc
             if(icy+1<=ncy)then                !! NE cell
                jc = ics + ncx + 1;ic_count(ic)=ic_count(ic)+1;ic_link(ic,ic_count(ic))=jc
             end if          
             if(icy-1>=1)then                !! SE cell
                jc = ics - ncx + 1;ic_count(ic)=ic_count(ic)+1;ic_link(ic,ic_count(ic))=jc
             end if
          end if       
          if(icx-1>=1)then             !! West cell
             jc = ics - 1;ic_count(ic)=ic_count(ic)+1;ic_link(ic,ic_count(ic))=jc
             if(icy+1<=ncy)then                !! NW cell
                jc = ics + ncx - 1;ic_count(ic)=ic_count(ic)+1;ic_link(ic,ic_count(ic))=jc
             end if         
             if(icy-1>=1)then                !! SW cell
                jc = ics - ncx - 1;ic_count(ic)=ic_count(ic)+1;ic_link(ic,ic_count(ic))=jc
             end if         
          end if
          if(icy+1<=ncy) then             !! North cell
             jc = ics + ncx;ic_count(ic)=ic_count(ic)+1;ic_link(ic,ic_count(ic))=jc
          end if
          if(icy-1>=1)then             !! South cell
             jc = ics - ncx;ic_count(ic)=ic_count(ic)+1;ic_link(ic,ic_count(ic))=jc
          end if 
       end if
       
       !! Next sheet
       if(icz+1<=ncz)then
          ics = ic + nsheet      

          jc = ics;ic_count(ic)=ic_count(ic)+1;ic_link(ic,ic_count(ic))=jc
          if(icx+1<=ncx)then             !! East cell
             jc = ics+1;ic_count(ic)=ic_count(ic)+1;ic_link(ic,ic_count(ic))=jc
             if(icy+1<=ncy)then                !! NE cell
                jc = ics + ncx + 1;ic_count(ic)=ic_count(ic)+1;ic_link(ic,ic_count(ic))=jc
             end if          
             if(icy-1>=1)then                !! SE cell
                jc = ics - ncx + 1;ic_count(ic)=ic_count(ic)+1;ic_link(ic,ic_count(ic))=jc
             end if
          end if       
          if(icx-1>=1)then             !! West cell
             jc = ics - 1;ic_count(ic)=ic_count(ic)+1;ic_link(ic,ic_count(ic))=jc
             if(icy+1<=ncy)then                !! NW cell
                jc = ics + ncx - 1;ic_count(ic)=ic_count(ic)+1;ic_link(ic,ic_count(ic))=jc
             end if         
             if(icy-1>=1)then                !! SW cell
                jc = ics - ncx - 1;ic_count(ic)=ic_count(ic)+1;ic_link(ic,ic_count(ic))=jc
             end if         
          end if
          if(icy+1<=ncy) then             !! North cell
             jc = ics + ncx;ic_count(ic)=ic_count(ic)+1;ic_link(ic,ic_count(ic))=jc
          end if
          if(icy-1>=1)then             !! South cell
             jc = ics - ncx;ic_count(ic)=ic_count(ic)+1;ic_link(ic,ic_count(ic))=jc
          end if             
       end if                
               
    end do
      
    return
  end subroutine neighbour_cells3d   
!! ------------------------------------------------------------------------------------------------   
end module cell_griddin
