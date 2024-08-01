c    "Copyright 2009 Prof. Robert Dalrymple, Prof. M. Gomez Gesteira, Dr Benedict Rogers, 
c     Dr. Alejandro Crespo, Dr. Muthukumar Narayanaswamy, Dr Shan Zou, Dr Andrea Panizzo "
c
c    This file is part of SPHYSICS.
c
c    SPHYSICS is free software; you can redistribute it and/or modify
c    it under the terms of the GNU General Public License as published by
c    the Free Software Foundation; either version 3 of the License, or
c    (at your option) any later version.
c
c    SPHYSICS is distributed in the hope that it will be useful,
c    but WITHOUT ANY WARRANTY; without even the implied warranty of
c    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
c    GNU General Public License for more details.
c
c    You should have received a copy of the GNU General Public License
c    along with this program.  If not, see <http://www.gnu.org/licenses/>.


c       initial data for SPH models                                                                     72

      program PART2VTU_2D

      parameter(np_max = 2700000,i_PART_counter_max=20000)
      parameter(num_phase_max=9)
      character chartemp*40, name_orig*40
      character name_vtu*40, name_vtu2*12, name_vtu3*9
      character np_string3*3, np_string4*4, np_string5*5
      character np_string6*6, np_string7*7, np_string8*8
      character frame_string1*1, frame_string2*2, frame_string3*3
      character frame_string4*4, frame_string5*5, frame_string6*6
      character supp*4,supp3*3,supp2*2,supp1*1, zero_string
      character string1*100,string2*100,string3*100,string4*100
      character chartemp2*100
      CHARACTER(LEN=10) :: FMT,FMT1
      CHARACTER(LEN=1)  :: TAB,DQ
      
      real xp(np_max),yp(np_max),zp(np_max),up(np_max),wp(np_max)
      real vp(np_max),aBubbles(np_max)
      real p(np_max),rho(np_max),alpha(np_max),vort(np_max)
      real time(i_PART_counter_max), DT(i_PART_counter_max)
      integer nsurf(np_max),fltype(np_max)
      integer np_all(i_PART_counter_max), IT(i_PART_counter_max)
      integer iinp,ioutp
      real  DT1(i_PART_counter_max),DT2(i_PART_counter_max)  
      integer flag2d
      
      real np_phase_start(num_phase_max),np_phase(num_phase_max)
      real rho0_phase(num_phase_max),P0_phase(num_phase_max) 
      real gamma_phase(num_phase_max),viscos_phase(num_phase_max)
      real ST_coeff(num_phase_max),backgroundPressure
          
      TAB=CHAR(9)     
      FMT="(A)"
      FMT1="(2A)"
      DQ=CHAR(34)
           
c     %LOAD AND READ TIME,DT FROM FILE DT. THE FIRST HEADERLINE IN 
        open(unit=70,file='../data_directory/TIME_OUT',status='old')
        i_loop_finish = 0
        i_PART_counter = 0
        do while(i_loop_finish.eq.0)
          i_PART_counter = i_PART_counter + 1
          if(i_PART_counter.gt.i_PART_counter_max)then
            print*,'Number of entries in file DT exceeds max value'
            print*,'i_PART_counter.gt.i_PART_counter_max'
            print*,'Adjust i_PART_counter_max, i_PART_counter_max = ',
     &              i_PART_counter_max
            stop
          endif
          read(70,*,END = 76)time(i_PART_counter),
     &np_all(i_PART_counter),IT(i_PART_counter),DT(i_PART_counter)
           
          !Determine whether to exit loop
          if(i_loop_finish.eq.0)then
            i_loop_finish = i_loop_finish - 1
          endif
76          i_loop_finish = i_loop_finish + 1
          !print*
        enddo
        N_start = 1
        Nframes = i_PART_counter-1  !Why -2?
      write(6,*) "There are ",Nframes,"frames."


      !! JRCK addition...
!      write(6,*) "Enter starting frame"
!      read(*,*) N_start
      N_start=1
      
!      stop        
      ngrab = N_start-2
      do iframe=N_start,Nframes
              
c       % READ IN THE PART FILE FOR EACH FRAME
        ngrab=ngrab+1
        if(ngrab.lt.10) then
           write(supp1,'(i0)') ngrab
           name_orig='../data_directory/PART'//supp1
        end if
        if(ngrab.ge.10.and.ngrab.lt.100) then
           write(supp2,'(i2)') ngrab
           name_orig='../data_directory/PART'//supp2
        end if
        if(ngrab.ge.100.and.ngrab.lt.1000) then
           write(supp3,'(i3)') ngrab
           name_orig='../data_directory/PART'//supp3
        end if
        if(ngrab.ge.1000.and.ngrab.lt.10000) then
           write(supp,'(i4)') ngrab
           name_orig='../data_directory/PART'//supp
        end if
        write(supp,'(i4.4)') ngrab

        name_vtu ='../paraview_files/PART'//supp//'.vtu'
        print*, 'iframe, name_orig, name_vtu ',
     &              iframe, ' ',name_orig, name_vtu 
        
        iinp = 23
        ioutp = 24
        open(iinp,file=name_orig,status='old')
        open(ioutp,file=name_vtu,status='unknown')
c       % READ POSITION, VELOCITY, DENSITY, PRESSURE, MASS AND VORTICITY DATA FOR ALL PARTICLES                      
         npp=0 ! keeps track of number of particles
         np = np_all(iframe)!king(iframe+1)
         do i=1,np
             read(iinp,*,end=300) xp(i),zp(i),yp(i),up(i),wp(i),vp(i)
     &,p(i),rho(i),alpha(i),aBubbles(i),nsurf(i),vort(i)
            npp=npp+1
         enddo
300    np=npp

        close (iinp)
        print*,'np ',np
                                                                
              
       
201     format(a40)
202     format(a100)
203     format(a25,i7,a17,i7,a2)
211     format(a21)


c     % OUTPUT TO FILE IN VTU FORMAT 
        if(np.lt.1000)then       
          write(np_string3,'(i3.3)') np
        string4 = '  <Piece NumberOfPoints='//DQ//np_string3//DQ//' Numb
     &erOfCells='//DQ//np_string3//DQ//'>'
        elseif(np.lt.10000)then       
          write(np_string4,'(i4.4)') np
        string4 = '  <Piece NumberOfPoints='//DQ//np_string4//DQ//' Numb
     &erOfCells='//DQ//np_string4//DQ//'>'
        elseif(np.lt.100000)then       
          write(np_string5,'(i5.5)') np
        string4 = '  <Piece NumberOfPoints='//DQ//np_string5//DQ//' Numb
     &erOfCells='//DQ//np_string5//DQ//'>'
        elseif(np.lt.1000000)then       
          write(np_string6,'(i6.6)') np
        string4 = '  <Piece NumberOfPoints='//DQ//np_string6//DQ//' Numb
     &erOfCells='//DQ//np_string6//DQ//'>'
        elseif(np.lt.10000000)then       
          write(np_string7,'(i7.7)') np
        string4 = '  <Piece NumberOfPoints='//DQ//np_string7//DQ//' Numb
     &erOfCells='//DQ//np_string7//DQ//'>'
        elseif(np.lt.100000000)then       
          write(np_string8,'(i8.8)') np
        string4 = '  <Piece NumberOfPoints='//DQ//np_string8//DQ//' Numb
     &erOfCells='//DQ//np_string8//DQ//'>'
        else
          print*,'Too many particles for np_string'
          stop  
        endif
        !print*,'np_string, np ',np_string, np 
        string1 = '<?xml version='//DQ//'1.0'//DQ//'?>'
        string2 = '<VTKFile type= '//DQ//'UnstructuredGrid'//DQ//'  vers
     &ion= '//DQ//'0.1'//DQ//'  byte_order= '//DQ//'BigEndian'//DQ//'>'
        string3 = ' <UnstructuredGrid>'
        write(ioutp,211)string1
        write(ioutp,202)string2
        write(ioutp,202)string3
        write(ioutp,202)string4
              
c       % WRITE IN PRESSURE DATA
        string1 = '   <PointData Scalars='//DQ//'Pressure'//DQ//' Vector
     &s='//DQ//'Velocity'//DQ//'>'
        string2 = '    <DataArray type='//DQ//'Float32'//DQ//' Name='//D
     &Q//'Pressures'//DQ//' format='//DQ//'ascii'//DQ//'>'
        write(ioutp,202)string1
        write(ioutp,202)string2
        do ii=1,np
          write(ioutp,*)p(ii)
        enddo
        string3 = '    </DataArray>'
        write(ioutp,202)string3

c       % WRITE density DATA        
        string1 = '    <DataArray type='//DQ//'Float32'//DQ//' Name='//D
     &Q//'Density'//DQ//' format='//DQ//'ascii'//DQ//'>'
        write(ioutp,202)string1
        do ii=1,np
          write(ioutp,*)rho(ii)
        enddo
        string3 = '    </DataArray>'
        write(ioutp,202) string3

c       % WRITE alpha DATA        
        string1 = '    <DataArray type='//DQ//'Float32'//DQ//' Name='//D
     &Q//'Alpha'//DQ//' format='//DQ//'ascii'//DQ//'>'
        write(ioutp,202)string1
        do ii=1,np
          write(ioutp,*)alpha(ii)
        enddo
        string3 = '    </DataArray>'
        write(ioutp,202) string3
        
c       % WRITE aB DATA        
        string1 = '    <DataArray type='//DQ//'Float32'//DQ//' Name='//D
     &Q//'aB'//DQ//' format='//DQ//'ascii'//DQ//'>'
        write(ioutp,202)string1
        do ii=1,np
          write(ioutp,*)aBubbles(ii)
        enddo
        string3 = '    </DataArray>'
        write(ioutp,202) string3        

c       % WRITE vort DATA        
        string1 = '    <DataArray type='//DQ//'Float32'//DQ//' Name='//D
     &Q//'vort'//DQ//' format='//DQ//'ascii'//DQ//'>'
        write(ioutp,202)string1
        do ii=1,np
          write(ioutp,*)vort(ii)
        enddo
        string3 = '    </DataArray>'
        write(ioutp,202) string3

c       % WRITE U-VELOCITY DATA        
        string1 = '    <DataArray type='//DQ//'Float32'//DQ//' Name='//D
     &Q//'u'//DQ//' format='//DQ//'ascii'//DQ//'>'
        write(ioutp,202)string1
        do ii=1,np
          write(ioutp,*)up(ii)
        enddo
        string3 = '    </DataArray>'
        write(ioutp,202) string3

c       % WRITE W_VELOCITY DATA        
        string1 = '    <DataArray type='//DQ//'Float32'//DQ//' Name='//D
     &Q//'v'//DQ//' format='//DQ//'ascii'//DQ//'>'
        write(ioutp,202)string1
        do ii=1,np
          write(ioutp,*)wp(ii)
        enddo
        string3 = '    </DataArray>'
        write(ioutp,202) string3
        
c       % WRITE V_VELOCITY DATA        
        string1 = '    <DataArray type='//DQ//'Float32'//DQ//' Name='//D
     &Q//'w'//DQ//' format='//DQ//'ascii'//DQ//'>'
        write(ioutp,202)string1

        do ii=1,np
          write(ioutp,*)vp(ii)
        enddo
        string3 = '    </DataArray>'
        write(ioutp,202) string3        

c       % WRITE surface DATA
        string1 = '    <DataArray type='//DQ//'Int32'//DQ//' Name='//D
     &Q//'surface'//DQ//' format='//DQ//'ascii'//DQ//'>'
        write(ioutp,202)string1
        do ii=1,np
          write(ioutp,*)nsurf(ii)
        enddo
        string3 = '    </DataArray>'
        write(ioutp,202) string3
              
c       % WRITE VELOCITY DATA
        string1 = '    <DataArray type='//DQ//'Float32'//DQ//' Name='//D
     &Q//'Velocity'//DQ//' NumberOfComponents='//DQ//'3'//DQ//' format='
     &//DQ//'ascii'//DQ//'>'
        write(ioutp,202) string1
        do ii=1,np
          write(ioutp,*)up(ii),wp(ii),vp(ii)
        enddo
        string3 = '    </DataArray>'
        write(ioutp,202) string3
        string4 = '   </PointData>'
        write(ioutp,202) string4
              
c       % WRITE PARTICLE POSITION DATA
        string2 = '   <Points>'
        string1 = '    <DataArray type='//DQ//'Float32'//DQ//' NumberOfCo
     &omponents='//DQ//'3'//DQ//' format='//DQ//'ascii'//DQ//'>'
        write(ioutp,202) string2
        write(ioutp,202) string1
        do ii=1,np
          write(ioutp,*)xp(ii),zp(ii),yp(ii)
        enddo        
        string3 = '    </DataArray>'
        string2 = '   </Points>'
        write(ioutp,202) string3
        write(ioutp,202) string2
             
c       % WRITE CELL DATA. CELL IS OF TYPE VERTEX.        
        string2 = '   <Cells>'
        string1 = '    <DataArray type='//DQ//'Int32'//DQ//' Name='//DQ/
     &/'connectivity'//DQ//' format='//DQ//'ascii'//DQ//'>'
        write(ioutp,202) string2
        write(ioutp,202) string1
        do ii=1,np
          write(ioutp,*)ii-1
        enddo
        string3 = '    </DataArray>'
        write(ioutp,202) string3
        
        
        string1 = '    <DataArray type='//DQ//'Int32'//DQ//' Name='//DQ/
     &/'offsets'//DQ//' format='//DQ//'ascii'//DQ//'>'
        write(ioutp,202) string1
        do ii=1,np
          write(ioutp,*)ii
        enddo
        string3 = '    </DataArray>'
        write(ioutp,202) string3
        
        
        string1 = '    <DataArray type='//DQ//'Int32'//DQ//' Name='//DQ/
     &/'types'//DQ//' format='//DQ//'ascii'//DQ//'>'
        write(ioutp,202) string1
        do ii=1,np
          write(ioutp,*)1
        enddo
        string3 = '    </DataArray>'
        write(ioutp,202) string3
        
        
        string1 = '   </Cells>' 
        string2 = '  </Piece>'
        string3 = ' </UnstructuredGrid>'
        string4 = '</VTKFile>'
        write(ioutp,202) string1
        write(ioutp,202) string2
        write(ioutp,202) string3
        write(ioutp,202) string4
      enddo
      close(ioutp);
      stop
      end
