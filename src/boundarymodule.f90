
module boundarymodule

  implicit none

  contains

  ! =====================================================================
  subroutine boundarycond ( time )

    ! --- modules ---
    use declarmodule
    use usefull
    use depsmodule

    implicit none

    integer                             :: elem , i, j, k , l    , ind   , it            ! Indexes
  
    integer                             :: auxint
    integer                             :: auxsurf(3)
    real*8                              :: auxarea 
    real*8                              :: aux1, aux2
    real*8 :: Xnode(3), CenterRot(3), AxisRot(3), MatRot(3,3)

    real*8 , intent(in) :: time
    real*8 :: AngleRot, toleranciapresion

    Udiri = 0.0d+0
    
    toleranciapresion = 1.0D-30
    
    ! --- Impose nodal conditions ---
    !
    do i = 1, nnod
      aux1 = NodBCMat(i) ! # of BC on i-th node
      !
      if ( BCDefMat( int(aux1) , 1 ) == 0 ) then ! es una traslacion
        !
        do j = 1, 3
          !
          aux2   = BCDefMat( int(aux1) , j+1 )  ! type BC for j-th coordinate
          auxint = 3*(i-1) + j                  ! degree of freedom
          !
          if (aux2 == 0 ) then ! it is neumann dof
            !
            Fext(auxint) = Fext(auxint) + BCDefMat(int(aux1),j+3+1) * time
            !
          elseif ( aux2 == 1 ) then ! it is diri
            !
            Udiri(auxint) =  BCDefMat( int(aux1) , j+3+1 ) * time
            AuxNeumDOFs(auxint) = 0
            
            !
          end if
        end do
        !
      elseif ( BCDefMat( int(aux1) , 1) == 1 ) then
        ! en este caso (hecho para
        ! desplazamientos) se imponen desplazamientos del tipo de una rotacion.
        !
        ! AxisRot
        Xnode = ldnodesdef(i,:)
        CenterRot = [1.0d+0 , 0.5d+0 , 0.5d+0 ]
        !
        AngleRot = BCDefMat( int(aux1) , 2) * time
        
        MatRot      = 0.0d+0
        MatRot(1,1) = 1.0d+0
        MatRot(2,2) =  cos( AngleRot )
        MatRot(2,3) =  sin( AngleRot )
        MatRot(3,2) = -sin( AngleRot )
        MatRot(3,3) =  cos( AngleRot )
        !
        Udiri( [ 3*(i-1)+1,3*(i-1)+2,3*(i-1)+3] ) =  &
          matmul( MatRot ,  Xnode - CenterRot )  - ( Xnode - CenterRot)
        AuxNeumDOFs([ 3*(i-1)+1,3*(i-1)+2,3*(i-1)+3]) = [ 0, 0, 0 ]

      !imposed disp normal direction
      elseif ( BCDefMat( int(aux1) , 1) == 2 ) then  
        ! aca se impone un desplazamiento radial con sim de rev
        ! con centro en centerrot y eje  Oz
        ! anglerot en este caso seria el desplazamiento radial
        !
        !  AxisRot
        Xnode = ldnodesdef(i,:)
        CenterRot = [0.0d+0 , 0.0d+0 , 0.5d+0 ]
        !
        AngleRot = BCDefMat( int(aux1) , 2) * time
        !
        Udiri( [ 3*(i-1)+1,3*(i-1)+2,3*(i-1)+3] ) =  [ ( Xnode(1:2) - CenterRot(1:2) ) / &
          sqrt( dot_product( ( Xnode(1:2) - CenterRot(1:2) ), ( Xnode(1:2) - CenterRot(1:2) ) ) ) &
          * AngleRot , [0.0D+0] ]
        AuxNeumDOFs([ 3*(i-1)+1,3*(i-1)+2,3*(i-1)+3]) = [ 0, 0, 0 ]

        !
      end if
    end do
    ! ---------------------------------------
    
!~     stop  

    
    ! --- Impone condiciones superficiales ---
    
    ! recorro todos los triangulos
    do i = 1, ntri
      !
      ! aux1 es el numero de cc de cada triangulo
      aux1 = SurfBCMat(i)
      !
      if ( BCDefMat( int(aux1) , 1 ) == 0 ) then ! es una traslacion
        !
        do j = 1, 3
          
          aux2    = BCDefMat( int(aux1) , j+1 )  ! type

          ! nodes of the surface
          auxsurf = [3*(Surfconec(i,1)-1) + j , 3*(Surfconec(i,2)-1) + j , 3*(Surfconec(i,3)-1) + j ]
          
          ! condicion de presion
          if (aux2 == 0 ) then
            !       
!~             if ( abs( BCDefMat(int(aux1),j+3+1) ) > toleranciapresion ) then                  
              auxarea = 0.5D+0 * &
                        sqrt( abs( M33DET( [ nodes( Surfconec(i,:) , 1 ) , &
                                             nodes( Surfconec(i,:) , 2 ) , &
                                            [1.0D+0 , 1.0D+0 , 1.0D+0  ] ] ) )**2 &
                             +abs( M33DET( [ nodes( Surfconec(i,:) , 1 ) , &
                                             nodes( Surfconec(i,:) , 3 ) , &
                                            [1.0D+0 , 1.0D+0 , 1.0D+0  ] ] ) )**2 &
                             +abs( M33DET( [ nodes( Surfconec(i,:) , 2 ) , &
                                             nodes( Surfconec(i,:) , 3 ) , &
                                            [1.0D+0 , 1.0D+0 , 1.0D+0  ] ] ) )**2 )
              !
              ! add loads to load vector
              
              Fext(auxsurf) = Fext(auxsurf) + auxarea * time * BCDefMat(int(aux1),j+3+1) * &
                [ 1.0D+0/3.0D+0, 1.0D+0/3.0D+0, 1.0D+0/3.0D+0 ]
          !
          ! condicion de desplazamiento
          elseif ( aux2 == 1 ) then
            !
            Udiri(auxsurf) =  &
              [ BCDefMat(int(aux1),j+3+1) , BCDefMat(int(aux1),j+3+1) , BCDefMat(int(aux1),j+3+1) ] * time
            AuxNeumDOFs(auxsurf) = [0,0,0]
            !
          end if
        end do
        !
      elseif ( BCDefMat( int(aux1) , 1) == 1 ) then
        !
        do j=1,3 ! for en los 3 nodos
          ind = Surfconec(i,j)
          Xnode = ldnodesdef( ind , : )

          CenterRot = [1.0d+0 , 0.5d+0 , 0.5d+0 ]
          AngleRot = BCDefMat( int(aux1) , 2) * time
          !
          MatRot      = 0.0d+0
          MatRot(1,1) = 1.0d+0
          MatRot(2,2) =  cos( AngleRot )
          MatRot(2,3) =  sin( AngleRot )
          MatRot(3,2) = -sin( AngleRot )
          MatRot(3,3) =  cos( AngleRot )
          !
          Udiri( [ 3*(ind-1)+1,3*(ind-1)+2,3*(ind-1)+3] ) =  &
            matmul( MatRot ,  Xnode - CenterRot )   - ( Xnode - CenterRot)
          AuxNeumDOFs([ 3*(ind-1)+1,3*(ind-1)+2,3*(ind-1)+3]) = [ 0, 0, 0 ]
          !
        end do
        !

      !imposed disp normal direction
      elseif ( BCDefMat( int(aux1) , 1) == 2 ) then  
        !
        do j=1,3 ! for en los 3 nodos
          ind = Surfconec(i,j)

          Xnode = ldnodesdef(ind,:)
          CenterRot = [0.0d+0 , 0.0d+0 , 0.5d+0 ]
          !
          AngleRot = BCDefMat( int(aux1) , 2) * time
          !
          Udiri( [ 3*(ind-1)+1,3*(ind-1)+2,3*(ind-1)+3] ) =  [ ( Xnode(1:2) - CenterRot(1:2) ) / &
            sqrt( dot_product( ( Xnode(1:2) - CenterRot(1:2) ), ( Xnode(1:2) - CenterRot(1:2) ) ) ) &
            * AngleRot , [0.0D+0] ]
        if (Xnode(1)==0.05) then
            write(*,*) "Xnode ", Xnode
            write(*,*) "Udiri: ", Udiri( [ 3*(ind-1)+1,3*(ind-1)+2,3*(ind-1)+3] )
            stop
          end if
          AuxNeumDOFs([ 3*(ind-1)+1,3*(ind-1)+2,3*(ind-1)+3]) = [ 0, 0, 0 ]
  
        end do
        !
      end if

    end do
    !
  end subroutine boundarycond

end module boundarymodule
