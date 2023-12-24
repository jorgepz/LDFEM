
module elemmodule

  implicit none

  contains


  ! ======================================================================
  subroutine ShapeFun( elemtype , derivord , x, y , z , fun )

    ! --- modules ---
    use declarmodule

    implicit none

    ! --- declarations ---
    integer,intent(in)   :: elemtype, derivord
    real*8 , intent(in)  :: x,y,z
    real*8               :: A(4,4), Ainv(4,4), detA
    !
    integer              :: i, j, k , l    , ind               ! Indexes
    integer              :: auxint(12)
    !
    real*8,allocatable,intent(out)          :: fun(:,:)
    ! ---------------------
    
    if (elemtype == 1) then

      if (derivord == 0) then
      
        allocate(fun(4,1))
       
        fun(1,1) = x
        fun(2,1) = 1.0d+0 - x - y - z
        fun(3,1) = z
        fun(4,1) = y

      else if (derivord == 1 ) then

        allocate(fun(3,4))

        fun = 0.0d+0  
        fun(1,1) = 1.0d+0
        fun(1:3,2) = [ -1.0d+0 , -1.0d+0 , -1.0d+0 ]
        fun(3,3) = 1.0d+0
        fun(2,4) = 1.0d+0
    
      end if
      
    end if
      
    return

  end subroutine ShapeFun
  ! ======================================================================




  ! ======================================================================
  subroutine EpsMats ( EpsMat , deriv )

  ! --- modules ---
    use declarmodule
    
    implicit none

    ! --- declarations ---
    real*8 , intent(in)  :: deriv(:,:)
    real*8               :: EpsMat(6,12)
    integer :: i
    ! ---------------------

    EpsMat = 0.0d+0
    
    do i=1,4

      EpsMat (1 , (i-1)*3+1   ) = deriv(1,i)
      EpsMat (2 , (i-1)*3+2   ) = deriv(2,i)
      EpsMat (3 , (i-1)*3+3   ) = deriv(3,i)
      
      EpsMat (4 , (i-1)*3+2   ) = deriv(3,i)
      EpsMat (4 , (i-1)*3+3   ) = deriv(2,i)
      
      EpsMat (5 , (i-1)*3+1   ) = deriv(3,i)
      EpsMat (5 , (i-1)*3+3   ) = deriv(1,i)
      
      EpsMat (6 , (i-1)*3+1   ) = deriv(2,i)
      EpsMat (6 , (i-1)*3+2   ) = deriv(1,i)

    end do

    return

  end subroutine EpsMats
  ! ======================================================================





  ! ======================================================================
  subroutine BgrandeMats ( BgrandeMat , deriv , F )

  ! --- modules ---
    use declarmodule
    
    implicit none

    ! --- declarations ---
    real*8 , intent(in)  :: deriv(:,:), F(:,:)
    real*8               :: BgrandeMat(6,12)
    integer              :: i, j, k
    ! ---------------------

    BgrandeMat = 0.0d+0
    
    do k=1,4

      do i=1,3
        do j=1,3
          BgrandeMat ( i , (k-1)*3 + j   ) = deriv(i,k) * F(j,i)
        end do
      end do
            

      do j=1,3
        BgrandeMat ( 4 , (k-1)*3 + j   ) = deriv(2,k) * F(j,3) + deriv(3,k) * F(j,2)
        BgrandeMat ( 5 , (k-1)*3 + j   ) = deriv(1,k) * F(j,3) + deriv(3,k) * F(j,1)
        BgrandeMat ( 6 , (k-1)*3 + j   ) = deriv(1,k) * F(j,2) + deriv(2,k) * F(j,1)
      end do

    end do

    return

  end subroutine BgrandeMats
  ! ======================================================================




  ! ======================================================================
  subroutine BsCalc( elemtype )

    ! --- modules ---
    use declarmodule
    use usefull
    
    implicit none

    ! --- declarations ---
    integer,intent(in)   :: elemtype
    real*8               :: Xelem(4),Yelem(4),Zelem(4)
    real*8               :: A(4,4), Ainv(4,4), detA
    !
    integer              :: i, j, k , l    , ind               ! Indexes
    integer              :: auxint(12)
    !
    logical              :: okflag
    ! ---------------------
    
    if (elemtype == 1) then
      
      do i = 1,ntet

        ! find nodal coordinates
        Xelem = nodes(Cellsconec(i,1:4),1)
        Yelem = nodes(Cellsconec(i,1:4),2)
        Zelem = nodes(Cellsconec(i,1:4),3)

        A (1:4,1) = dble ( [1,1,1,1] )
        A (1:4,2) = dble ( Xelem )
        A (1:4,3) = dble ( Yelem )
        A (1:4,4) = dble ( Zelem )

        call M44INV ( A , Ainv , detA , okflag )

        ElemVols(i) = detA/6.0d+0 
        !
        if (detA<=0.0) then
          write(*,*) "Error volumen ",i, "negativoo!! ", ElemVols(i)
          write(*,*) "hay como: ", Ntet, "elementos en total"
!~           stop
        end if

        auxint = [ ((i-1)*12 + j , j=1,12) ]

        Ainv = transpose(Ainv)

        do j = 1,4

          ! ------------
          BsMat (1 , (i-1)*12 + (j-1)*3+1   ) = Ainv(j,2)
          BsMat (2 , (i-1)*12 + (j-1)*3+1 +1) = Ainv(j,3)
          BsMat (3 , (i-1)*12 + (j-1)*3+1 +2) = Ainv(j,4)

          BsMat (4 , (i-1)*12 + (j-1)*3+1 +1) = Ainv(j,4)
          BsMat (4 , (i-1)*12 + (j-1)*3+1 +2) = Ainv(j,3)

          BsMat (5 , (i-1)*12 + (j-1)*3+1   ) = Ainv(j,4)
          BsMat (5 , (i-1)*12 + (j-1)*3+1 +2) = Ainv(j,2)

          BsMat (6 , (i-1)*12 + (j-1)*3+1   ) = Ainv(j,3)
          BsMat (6 , (i-1)*12 + (j-1)*3+1 +1) = Ainv(j,2)
          ! ------------

          BMat (j,((i-1)*3+1):((i-1)*3+3)) = Ainv (j,2:4)

        end do  ! end do j
      end do  ! end do elem
    end if ! end if elemtye

    return

  end subroutine BsCalc
  ! ======================================================================




  ! ======================================================================
  subroutine GradHessCalc ( elemtype )

    ! --- modules ---
    use declarmodule
    use usefull

    implicit none

    ! --- declarations ---
    integer,intent(in)   :: elemtype
    !
    integer              :: i, j, k , l    , ind      ,elem         ! Indexes
    integer              :: auxint(12)
    !
    real*8               :: auxreal1, auxreal2        , ojo(3,3)  , auxtrace, caca(3,3)

    real*8               :: Xelemdef(4), Yelemdef(4), Zelemdef(4)
    real*8               :: Mk(3,4), BGe(4,3), Fk(3,3) , Ek(3,3) , Sk(3,3)    , invFK(3,3)
    real*8               :: Mij(3,3), Mkl(3,3)
    real*8               :: PsiHess (6,6)

    ! these were used for debugging
    ! real*8               :: GradMatE ( Nnodes*3 ), HessMatE(Nnodes*3,Nnodes*3), CellsconecE(Nelem,4)
    logical :: okflag

    real*8 :: voigMij(6), voigMkl(6)   , detFk
    ! ------------------------------

    ojo      = 0.0d+0
    ojo(1,1) = 1.0d+0
    ojo(2,2) = 1.0d+0
    ojo(3,3) = 1.0d+0

    ! psi hessian for saint venant kirchoff function
    PsiHess = 0.0d+0
    PsiHess (1,1:3) = ( shear / (1.0d+0-2.0d+0*nu) ) * 2.0d+0 * [ 1.0d+0-nu , nu     , nu ] 
    PsiHess (2,1:3) = ( shear / (1.0d+0-2.0d+0*nu) ) * 2.0d+0 * [ nu     , 1.0d+0-nu , nu ] 
    PsiHess (3,1:3) = ( shear / (1.0d+0-2.0d+0*nu) ) * 2.0d+0 * [ nu     , nu     , 1.0d+0-nu ] 
    PsiHess (4,4) = shear 
    PsiHess (5,5) = shear 
    PsiHess (6,6) = shear 
    ! -------------------


    if (elemtype == 1) then

      HessMat      = 0.0d+0
      GradMat      = 0.0d+0
      RightHand    = 0.0D+0
      !
      ! === compute grad and hess ===
      do elem = 1 , Ntet
        !
        Xelemdef = ldnodesdef ( Cellsconec(elem,1:4) , 1 )
        Yelemdef = ldnodesdef ( Cellsconec(elem,1:4) , 2 )
        Zelemdef = ldnodesdef ( Cellsconec(elem,1:4) , 3 )
        !

        Mk (1 , 1:4) = Xelemdef
        Mk (2 , 1:4) = Yelemdef
        Mk (3 , 1:4) = Zelemdef
        !

        BGe = BMat (1:4, ((elem-1)*3+1):(elem*3) )
        !
        Fk = matmul ( Mk , BGe )
        
        call M33INV (Fk , invFK , okflag , detFk )        

!~          if (elem==1) then
!~     write (*,*) "------"
!~     write (*,*) "Fk:", Fk
!~   end if
  
        if (detFk < 0.0) then
          !
          write(*,*) "Error - detFk negativo: ",detFk , " en elemento: ", elem
          !
        end if
      
        ! --- lagrange tensor ---
        Ek = 0.5d+0 * ( matmul( transpose(Fk) , Fk ) - ojo )


        call trace ( Ek , auxtrace )

        ! --- cosserat tensor ---
        Sk = lambda * auxtrace * ojo  +  2.0d+0 * shear * Ek

!~          if (elem==1) then
!~     write (*,*) "------"
!~     write (*,*) "Sk:", Sk
!~   end if

!~   if (elem==1) then
!~     write (*,*) "------"
!~     write (*,*) "Sk:", Sk
!~   end if

        !
!~                    write(*,'(A)') "Sk"
!~                    write(*,'(3f8.3)') Sk
        ! --- loop in coordinates ---
        do i = 1, 3
          !
          ! --- loop in nodes ---
          do j = 1, 4
            !
            Mij = transpose ( matmul ( &
                        transpose ( reshape( Fk(i,:) , [1,3])  ) , &
                        reshape( BGe (j,:) , [1,3])  ) )
            
            ! --- inner product --- 
            call tensinneprod ( Sk , Mij , auxreal1 )



            ! ----------------------
            ! --- Gradient ---------
            GradMat     ( (Cellsconec(elem,j)-1)*3+i ) = GradMat     ( (Cellsconec(elem,j)-1)*3+i ) &
              + auxreal1 * ElemVols(elem)
            ! ----------------------


            ! --- loop in coordinates ---
            do k=1,3
              !
              ! --- loop in nodes ---
              do l=1,4
                !
                Mkl = transpose ( matmul ( &
                        transpose ( reshape( Fk(k,:) , [1,3])  ) , &
                        reshape( BGe (l,:) , [1,3])  ) )

                !
                if ( i==k ) then
                  !
                  caca = transpose ( matmul (  transpose ( reshape( BGe(l,:) , [1,3])  ) , reshape( BGe (j,:) , [1,3])  ) )

                  call tensinneprod ( Sk , caca , auxreal1 )
                  !
                else
                  !
                  auxreal1 = 0.0d+0
                end if
                
                
                call trace ( Mkl , auxtrace )
                
                voigMkl = 0.0d+0
                voigMij = 0.0d+0
                
                call tranvoigsym(Mkl,voigMkl)
                call tranvoigsym(Mij,voigMij)
                                
                voigMij(4:6)  = voigMij(4:6)  * 2.0d+0
                voigMkl(4:6)  = voigMkl(4:6)  * 2.0d+0
                                
                auxreal2 = dot_product ( voigMkl , matmul (PsiHess , voigMij ) )

                ! ----------------------
                ! ---- Hessian  ---------
                HessMat   ( (Cellsconec(elem,j)-1)*3+i  , (Cellsconec(elem,l)-1)*3+k ) = &
                  HessMat ( (Cellsconec(elem,j)-1)*3+i  , (Cellsconec(elem,l)-1)*3+k ) + &
                  ( auxreal2 + auxreal1 )  * ElemVols(elem)
                  
                  
                ! ----------------------

              end do  ! end loop node l
              !
            end do  ! end loop coord k
            !
          end do ! end loop node j
          !
        end do  ! end loop coord i



        ! --- faces loads calculation ---
!~         if ( any ( FacesNodLoad( 1:NFacesNodLoad , 1 ) == elem ) ) then
!~         
!~           !
!~           i = minloc( abs( FacesNodLoad( 1:NFacesNodLoad , 1 ) - elem ) ,1 )
!~           
!~           ! area triang ---
!~           auxreal1 = 0.5D+0 * abs( M33DET( [ nodes( FacesNodLoad(i , 2:4 ) , 2 ) , &
!~                                              nodes( FacesNodLoad(i , 2:4 ) , 3 ) , &
!~                                              [1.0D+0 , 1.0D+0 , 1.0D+0  ]       ] ) )
!~ 
!~           write(*,*) "face",i, "area: ", auxreal1, "elem: ",FacesNodLoad( i,1), "nodes: ", FacesNodLoad( i,2:4)
!~ 
!~           !
!~           do j=1,3
!~             RightHand ( FacesNodLoad(i,j+1)*3-2 ) = RightHand ( FacesNodLoad(i,j+1)*3-2 ) + loadx*auxreal1 / 3.0d+0 
!~             RightHand ( FacesNodLoad(i,j+1)*3-1 ) = RightHand ( FacesNodLoad(i,j+1)*3-1 ) + loady*auxreal1 / 3.0d+0 
!~           end do ! loop on 3 nodes of face
!~           !
!~         end if ! end if of element with loaded face


        !
      end do ! end loop elements

!~       write ( *,'(A)') "--- Check RightHand ---"
!~       write ( *,'(3e12.5)') RightHand
!~       write ( *,'(A)') "-----------------------"

      ! --- gradient final ---
      GradMat = GradMat - RightHand

    end if ! end if elemtype

    return

  end subroutine GradHessCalc
  ! ======================================================




  ! ======================================================================
  ! ======================================================
  subroutine VonMissesCalc ( elemtype )

    ! --- modules ---

    use declarmodule
    use usefull

    implicit none

    ! --- declarations ---
    integer,intent(in)   :: elemtype
    !
    integer              :: i, j, k , l    , ind      ,elem         ! Indexes
    integer              :: auxint(12)
    !
    real*8               :: auxreal1, auxreal2        , ojo(3,3)  , auxtrace, caca(3,3)

    real*8               :: Xelemdef(4), Yelemdef(4), Zelemdef(4)
    real*8               :: Mk(3,4), BGe(4,3), Fk(3,3) , Ek(3,3) , Sk(3,3)    , invFK(3,3)
    real*8               :: Mij(3,3), Mkl(3,3)
    real*8               :: PsiHess (6,6)

    ! these were used for debugging
    ! real*8               :: GradMatE ( Nnodes*3 ), HessMatE(Nnodes*3,Nnodes*3), CellsconecE(Nelem,4)
    logical :: okflag

    real*8 :: voigMij(6), voigMkl(6)   , detFk
    ! ------------------------------

    ojo      = 0.0d+0
    ojo(1,1) = 1.0d+0
    ojo(2,2) = 1.0d+0
    ojo(3,3) = 1.0d+0

    ! psi hessian for saint venant kirchoff function
    PsiHess = 0.0d+0
    PsiHess (1,1:3) = ( shear / (1.0d+0-2.0d+0*nu) ) * 2.0d+0 * [ 1.0d+0-nu , nu     , nu ] 
    PsiHess (2,1:3) = ( shear / (1.0d+0-2.0d+0*nu) ) * 2.0d+0 * [ nu     , 1.0d+0-nu , nu ] 
    PsiHess (3,1:3) = ( shear / (1.0d+0-2.0d+0*nu) ) * 2.0d+0 * [ nu     , nu     , 1.0d+0-nu ] 
    PsiHess (4,4)   = shear 
    PsiHess (5,5)   = shear 
    PsiHess (6,6)   = shear 
    ! -------------------

    if (elemtype == 1) then

      HessMat      = 0.0d+0
!~       HessMatAntes = 0.0d+0
      GradMat      = 0.0d+0
      RightHand    = 0.0D+0
      !
      ! === compute grad and hess ===
      do elem = 1 , Ntet
        !
        Xelemdef = ldnodesdef ( Cellsconec(elem,1:4) , 1 )
        Yelemdef = ldnodesdef ( Cellsconec(elem,1:4) , 2 )
        Zelemdef = ldnodesdef ( Cellsconec(elem,1:4) , 3 )
        !
!~         write(*,*)  "nodes def:  "
!~         write(*,'(3f9.3)') transpose ( nodesdef( Cellsconec(elem,1:4) , 1:3 ) )

        Mk (1 , 1:4) = Xelemdef
        Mk (2 , 1:4) = Yelemdef
        Mk (3 , 1:4) = Zelemdef
        !
!~         write(*,'(A,/,4f9.3,/,4f9.3,/,4f9.3)') "Mk: ", Mk(1,:), Mk(2,:), Mk(3,:)

        BGe = BMat (1:4, ((elem-1)*3+1):(elem*3) )
        !
        Fk = matmul ( Mk , BGe )
        
        call M33INV (Fk , invFK , okflag , detFk )        
        
        if (detFk <0.0) then
          !
          write(*,*) "Error - detFk negativo: ",detFk , " en elemento: ", elem
!~           stop
          !
        end if
      
        ! --- lagrange tensor ---
        Ek = 0.5d+0 * ( matmul( transpose(Fk) , Fk ) - ojo )
        
        call trace ( Ek , auxtrace )

        ! --- cosserat tensor ---
        Sk = lambda * auxtrace * ojo  +  2.0d+0 * shear * Ek
        !
      end do ! end loop elements
      !
    end if ! end if elemtype

    return

  end subroutine VonMissesCalc
  ! ======================================================


      
end module elemmodule
