

module usefull
  !
  implicit none
  !
  contains
  !
  subroutine linspace(x, x_start, x_end, x_len )
    !

    real*8, dimension(:), intent (out) :: x
    real*8                             :: x_start, x_end, dx
    integer                          :: x_len, i

    dx         = ( x_end - x_start ) / ( x_len - 1 )
    x(1:x_len) = [ ( x_start + ((i-1)*dx) , i=1 , x_len ) ]

  end subroutine


  subroutine tranvoigsym(T,v)
  
    implicit none

    real*8, intent(in) :: T(3,3)  
    real*8,intent(out) :: v(6)

    v= 0.0
    
    v(1) = T(1,1)
    v(2) = T(2,2)
    v(3) = T(3,3)
      
    v(4) = ( T(2,3)+T(3,2) ) * 0.5 
    v(5) = ( T(1,3)+T(3,1) ) * 0.5 
    v(6) = ( T(1,2)+T(2,1) ) * 0.5 
      
    return
    
  end subroutine tranvoigsym


  subroutine sparmatvecmul( inds , valA , vec , prod )
  
    integer :: inds(:,:), nvals, sizevec, i,j, auxi, auxj
    real*8  :: valA(:), vec(:)
    real*8, allocatable,intent(out)  :: prod(:)

    sizevec = size( vec  , 1 )
    nvals   = size( valA , 1 )
    
    allocate ( prod(sizevec) )
    
    prod = 0.0d+0
    
    do i=1, nvals
      !
      auxi = inds( i , 1 )
      auxj = inds( i , 2 )
      prod( auxi ) = prod( auxi ) + valA(i) * vec( auxj ) 
      if ( auxi .ne. auxj ) then
        prod( auxj ) = prod( auxj ) + valA(i) * vec( auxi ) 
      end if       
      !
    end do
  
    return
    
  end subroutine sparmatvecmul


  !
  !
  ! ---------------------------------------------------
  subroutine plot(x,y,title,scalesflag)
    !
    ! Example:   call plot( real*8( [2,3,7,4]) , real*8([2,8,5,10]) )
    !
    real*8, intent(in), dimension(:) :: x, y
    integer                          :: size_x,size_y,i
    character, intent(in)            :: title(6)
    integer , intent(in)             :: scalesflag

    size_x = size(x)
    size_y = size(y)

    if ( size_x /= size_y ) then
      print *, "Array size mismatch"
    else
      !
      ! Writes file with data
      open(unit=8,file='data.dat')
      !
      do i=1,size_x
        write(8,'(e16.8,e16.8)') x(i), y(i)
      end do
      !
      close(unit=8)
      ! -------------------------------------

      ! --- Writes file with plot format ----
      open(unit=8,file='style.gnu')
      !
      write(8,'(A)') "set terminal png"
      write(8,'(A,6A,A)') "set output '",title(1:6),".png'"
      write(8,'(A)') "set grid"
      write(8,'(A)') "set xlabel 'Iterations'"
      write(8,'(A,6A,A)') "set ylabel '",title,"'"

      if (scalesflag==1) then
        write(8,'(A)') "set logscale x"
      else if (scalesflag==2) then
        write(8,'(A)') "set logscale y"
      end if
      !
      write(8,'(A,A,A)') "plot 'data.dat' with linespoints ls 7 lc 9"
!~       write(8,'(A,A,A)') "plot ",title," with linespoints ls 7 lc 9"
      !
      close(unit=8)
      ! -------------------------------------

    end if

    call system('gnuplot style.gnu')
    call system('rm data.dat')


  end subroutine plot


  subroutine trace(A,tr)
    !
    implicit none
    !
    real*8, intent(in)    :: A(3,3)
    integer               :: j
    real*8                :: tr
    !
    tr = 0.0d+0
    !
    do j = 1, 3
      !
      tr = tr + A(j,j)
      !
    end do
    !
  end subroutine trace


  subroutine traceComplex(A,tr)
    !
    implicit none
    !
    double complex, intent(in)    :: A(3,3)
    integer                       :: j
    double complex                :: tr
    !
    tr = 0.0d+0
    !
    do j = 1, 3
      !
      tr = tr + A(j,j)
      !
    end do
    !
  end subroutine traceComplex

end module
