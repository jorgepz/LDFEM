
module vtkmodule

  implicit none

  contains

  subroutine vtkgen( vtkname, pointdataname, celldataname  &
                     , NNodPC , fileunit  &
                     , pointdata , celldata &
                     , vtkconec , vtknodes &
                     )

    use declarmodule
    
    implicit none

    ! ------------------------------------------------------------------
    character(len=8), intent(in)  :: vtkname ,pointdataname, celldataname
  
    integer                :: i, j, k                      ! Indexes
    integer  , intent(in)  :: NNodPC, fileunit

    real*8    ,intent(in)    :: vtknodes (nnod,3)     ! Nodes coordinates matrix
    integer ,intent(in)    :: vtkconec (ntet ,4)            ! Nodes coordinates matrix 
    integer                :: pointdatacomp, celldatacomp
    real*8     , intent(in)  :: pointdata(:,:), celldata(:,:)    ! Nodes coordinates matrix
    ! ------------------------------------------------------------------

!~ write(*,*) vtkname


    !
    open (unit = fileunit , file = vtkname//".vtk" )
    !
    write (fileunit,'(A)') "# vtk DataFile Version 2.0"
    write (fileunit,'(A)') vtkname
    write (fileunit,'(A)') "ASCII"
    write (fileunit,'(A)') "DATASET UNSTRUCTURED_GRID"
    write (fileunit,'(A)') ""
    write (fileunit,'(A,i10,A)') "POINTS", nnod ," float"
    !


    do i=1,nnod
      !
      write (fileunit,'(3e25.10)') ( vtknodes(i,j),j=1,3 )
      !
    end do
    !


    write (fileunit,'(A)') ""

    write (fileunit,'(A,i10,i10)') "CELLS", ntet ,ntet*(NNodPC+1)
    !
    do i= 1 , ntet
      !
      write (fileunit,'(i1,A,5i12)') NNodPC," ",( vtkconec(i,j)-1,j=1,NNodPC)
      !
    end do
  

    write (fileunit,'(A)') ""


    if ( NNodPC == 4) then
      !
      write (fileunit,'(A,i10)') "CELL_TYPES", ntet
      do i=1,ntet
        write (fileunit,'(A)') "10"
      end do

    else if ( NNodPC == 8) then
      
      write (fileunit,'(A,i10)') "CELL_TYPES", ntet
      do i=1,ntet
        write (fileunit,'(A)') "12"
      end do

    end if
  
    pointdatacomp = size(pointdata,2)

    !
    write (fileunit,'(A)') ""
    write (fileunit,'(A,i10)') "POINT_DATA", nnod
    write (fileunit,'(A,i1)') "SCALARS "//pointdataname//" float ", pointdatacomp
    write (fileunit,'(A)') "LOOKUP_TABLE default"
    !
    do i=1,nnod
      !
      write (fileunit,'(3e25.10)') (pointdata(i,j),j=1,pointdatacomp)
      !
    end do
    !

    !------------------------------
    celldatacomp = size( celldata , 2 )

    write (fileunit,'(A)') ""
    write (fileunit,'(A,i10)') "CELL_DATA", ntet
    write (fileunit,'(A,i1)') "SCALARS "//celldataname//" float ", celldatacomp
    write (fileunit,'(A)') "LOOKUP_TABLE default"
    !
    do i=1,ntet
      !
      do j=1,celldatacomp
        if (j<celldatacomp) then
          write (fileunit,'(e12.5)',advance='no') (celldata(i,j))
        else
          write (fileunit,'(e12.5)') (celldata(i,j))
        end if
      end do
      !
    end do
    !------------------------------
    close(fileunit)

    return
    
  end subroutine vtkgen
  
!~ ####################################################################    
!~ ####################################################################     
!~     Comienzo edición de Pablo - Creación de nuevas subrutinas de vtk  05-11-2014
!~ ####################################################################
!~ ####################################################################
  
  subroutine vtkgen_point( vtkname, pointdataname  &
                          , NNodPC , fileunit  &
                          , pointdata  &
                          , vtkconec , vtknodes &
                          )

    use declarmodule
    
    implicit none

    ! ------------------------------------------------------------------
    character(len=8), intent(in)  :: vtkname ,pointdataname
  
    integer                :: i, j, k                      ! Indexes
    integer  , intent(in)  :: NNodPC, fileunit

    real*8    ,intent(in)    :: vtknodes (nnod,3)     ! Nodes coordinates matrix
    integer ,intent(in)    :: vtkconec (ntet ,4)            ! Conectivity matrix 
    integer                :: pointdatacomp             
    real*8     , intent(in)  :: pointdata(:,:)    ! Point data
    ! ------------------------------------------------------------------

    !
    open (unit = fileunit , file = vtkname//".vtk" )
    !
    write (fileunit,'(A)') "# vtk DataFile Version 2.0"
    write (fileunit,'(A)') vtkname
    write (fileunit,'(A)') "ASCII"
    write (fileunit,'(A)') "DATASET UNSTRUCTURED_GRID"
    write (fileunit,'(A)') ""
    write (fileunit,'(A,i10,A)') "POINTS", nnod ," float"
    !


    do i=1,nnod
      !
      write (fileunit,'(3e25.10)') ( vtknodes(i,j),j=1,3 )
      !
    end do
    !


    write (fileunit,'(A)') ""

    write (fileunit,'(A,i10,i10)') "CELLS", ntet ,ntet*(NNodPC+1)
    !
    do i= 1 , ntet
      !
      write (fileunit,'(i1,A,5i12)') NNodPC," ",( vtkconec(i,j)-1,j=1,NNodPC)
      !
    end do
  

    write (fileunit,'(A)') ""


    if ( NNodPC == 4) then
      !
      write (fileunit,'(A,i10)') "CELL_TYPES", ntet
      do i=1,ntet
        write (fileunit,'(A)') "10"
      end do

    else if ( NNodPC == 8) then
      
      write (fileunit,'(A,i10)') "CELL_TYPES", ntet
      do i=1,ntet
        write (fileunit,'(A)') "12"
      end do

    end if
  
    pointdatacomp = size(pointdata,2)

    !
    write (fileunit,'(A)') ""
    write (fileunit,'(A,i10)') "POINT_DATA", nnod
    write (fileunit,'(A,i1)') "SCALARS "//pointdataname//" float ", pointdatacomp
    write (fileunit,'(A)') "LOOKUP_TABLE default"
    !
    do i=1,nnod
      !
      write (fileunit,'(3e25.10)') (pointdata(i,j),j=1,pointdatacomp)
      !
    end do
    !

    close(fileunit)

    
  end subroutine vtkgen_point
  
  
  subroutine vtkgen_cell( vtkname, celldataname  &
                     , NNodPC , fileunit  &
                     , celldata &
                     , vtkconec , vtknodes &
                     )

    use declarmodule
    
    implicit none

    ! ------------------------------------------------------------------
    character(len=8), intent(in)  :: vtkname , celldataname
  
    integer                :: i, j, k                      ! Indexes
    integer  , intent(in)  :: NNodPC, fileunit

    real*8    ,intent(in)    :: vtknodes (nnod,3)     ! Nodes coordinates matrix
    integer ,intent(in)    :: vtkconec (ntet ,4)            ! Conectivity matrix 
    integer                ::  celldatacomp
    real*8     , intent(in)  ::  celldata(:,:)    ! Data cell
    ! ------------------------------------------------------------------

    !
    open (unit = fileunit , file = vtkname//".vtk" )
    !
    write (fileunit,'(A)') "# vtk DataFile Version 2.0"
    write (fileunit,'(A)') vtkname
    write (fileunit,'(A)') "ASCII"
    write (fileunit,'(A)') "DATASET UNSTRUCTURED_GRID"
    write (fileunit,'(A)') ""
    write (fileunit,'(A,i10,A)') "POINTS", nnod ," float"
    !


    do i=1,nnod
      !
      write (fileunit,'(3f12.5)') ( vtknodes(i,j),j=1,3 )
      !
    end do
    !


    write (fileunit,'(A)') ""

    write (fileunit,'(A,i10,i10)') "CELLS", ntet ,ntet*(NNodPC+1)
    !
    do i= 1 , ntet
      !
      write (fileunit,'(i1,A,5i12)') NNodPC," ",( vtkconec(i,j)-1,j=1,NNodPC)
      !
    end do
  

    write (fileunit,'(A)') ""


    if ( NNodPC == 4) then
      !
      write (fileunit,'(A,i10)') "CELL_TYPES", ntet
      do i=1,ntet
        write (fileunit,'(A)') "10"
      end do

    else if ( NNodPC == 8) then
      
      write (fileunit,'(A,i10)') "CELL_TYPES", ntet
      do i=1,ntet
        write (fileunit,'(A)') "12"
      end do

    end if
  

    !------------------------------
    celldatacomp = size( celldata , 2 )

    write (fileunit,'(A)') ""
    write (fileunit,'(A,i10)') "CELL_DATA", ntet
    write (fileunit,'(A,i1)') "SCALARS "//celldataname//" float ", celldatacomp
    write (fileunit,'(A)') "LOOKUP_TABLE default"
    !
    do i=1,ntet
      !
      do j=1,celldatacomp
        if (j<celldatacomp) then
          write (fileunit,'(e12.5)',advance='no') (celldata(i,j))
        else
          write (fileunit,'(e12.5)') (celldata(i,j))
        end if
      end do
      !
    end do
    !------------------------------
    close(fileunit)

    
  end subroutine vtkgen_cell

!~ ####################################################################    
!~ ####################################################################     
!~     Fin edición de Pablo - Creación de nuevas subrutinas de vtk  05-11-2014
!~ ####################################################################
!~ ####################################################################

  subroutine grid2vtk( nodes , cells, grid_scalars )
    !
    !     % http://www.visualization.hpc.mil/wiki/VTK_Data_Formats_%26_Reading_in_Data
    !
    integer :: cells (:,:)
    !
    real*8 :: grid_scalars(:), nodes(:,:)
    integer :: nnodes, i,j, nelem, caras
    !
    nnodes = size(nodes,1)
    !

    !
    open (unit=4,file='dist_vtk.vtk')
    !
    write (4,'(A)') '# vtk DataFile Version 2.0'
    write (4,'(A)') 'dist_vtk'
    !
    write (4,'(A,/,A,/)') 'ASCII','DATASET POLYDATA'
    !
    write (4,'(A,i6,A)')  'POINTS ', nnodes, ' float'
    !
    do i = 1, nnodes
      !
      write (4,'(4f8.2)')  nodes(i,:)
!~       fprintf(fid,'%3.2f %3.2f %3.2f \n', nodes(i,:) ) ;
      !
    end do
    !
    nelem = size(cells,1)
    !
    caras = 4
    !
    write (4,'(/,A,i8,i8)') 'POLYGONS ', nelem , nelem*(caras+1)

    do i = 1, nelem
      !
      write (4,'(i3,4i9)')  caras, cells(i,:)-1
      !
    end do
    !
    !
    write (4,'(/,A,i6)') 'POINT_DATA ', nnodes
    !
    ! types
    write (4,'(A)') 'SCALARS regions float'
!~     write (4,'(A)') 'SCALARS regions int 1'
    write (4,'(A)') 'LOOKUP_TABLE default'
    !
    do i = 1, nnodes
      write (4,'(f10.5)')  grid_scalars(i)
    end do
    !
    close(4)
    !
    !
    !
  end subroutine grid2vtk

  
end module vtkmodule
