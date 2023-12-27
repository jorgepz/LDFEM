!
module materialmodule

  implicit none

  contains

  ! ======================================================================
  ! Cosserat Tensor function
  ! ======================================================================
  subroutine cosserat ( E, FisPropElem , S )

    ! --- modules ---
    use declarmodule
    use usefull
    use depsmodule

    implicit none

    integer                             :: elem , i, j, k , l    , ind   , it, FisPropElem            ! Indexes
  
    real*8                              :: auxtrace, ojo(3,3)

    real*8                              :: time
    real*8 , intent(in)                 :: E(3,3)
    real*8 , intent(out)                :: S(3,3)
    logical                             :: ok_flag
    real*8                              :: C(3,3), CINV(3,3)
    real*8                              :: DETC

    ! identity
    ojo      = 0.0d+0
    ojo(1,1) = 1.0d+0
    ojo(2,2) = 1.0d+0
    ojo(3,3) = 1.0d+0
    
    
    if (constitutive_model == 1 ) then ! svk model
      !
      ! --- Materials ---
      young = MatsMat(FisPropElem,1)
      nu    = MatsMat(FisPropElem,2)
      !
      lambda  = young * nu / ( (1.0d+0 + nu) * (1.0d+0 - 2.0d+0 * nu) )
      shear   = young      / ( 2.0d+0 * (1.0d+0 + nu) )
      ! ---------------------------
!~       write(*,*) lambda,shear
!~       stop
      !
      call trace ( E , auxtrace )
      S = lambda * auxtrace * ojo  +  2.0d+0 * shear * E
      !
    elseif (constitutive_model == 2) then ! curnier model
      !
      lambda = MatsMat(FisPropElem,1)
      shear  = MatsMat(FisPropElem,2)
      !
      C = 2.0D+0 * E + ojo
      call M33INV ( C , CINV, OK_FLAG )
      DETC = M33DET(C)
      S = lambda * ( sqrt(DETC) -1.0D+0 )* CINV +  2.0d+0 * shear * E
    end if
    
  end subroutine cosserat
  ! ======================================================================



  ! ======================================================================
  ! Cosserat Tensor function
  ! ======================================================================
  subroutine cosseratincomp ( E , FisPropElem, p , S )

    ! --- modules ---
    use declarmodule
    use usefull
    use depsmodule

    implicit none

    integer                             :: elem , i, j, k , l    , ind   , it   , FisPropElem          ! Indexes
  
    real*8                              :: auxtrace, ojo(3,3)

    real*8                              :: time
    real*8 , intent(in)                 :: E(3,3), p
    real*8 , intent(out)                :: S(3,3)
    logical                             :: ok_flag
    real*8                              :: C(3,3), CINV(3,3)
    real*8                              :: DETC, DETF, trC
    real*8 :: c1
    
    ! identity
    ojo      = 0.0d+0
    ojo(1,1) = 1.0d+0
    ojo(2,2) = 1.0d+0
    ojo(3,3) = 1.0d+0
    
    
    if (constitutive_model == 6 ) then ! neohookean
      !
      c1 = MatsMat(FisPropElem,1)
      !
      C = 2.0D+0 * E + ojo
      call M33INV ( C , CINV, OK_FLAG )
      DETC = M33DET(C)
      DETF = sqrt(DETC)
      call trace( C, trC )
      
      S =  DETF * p * CINV  +  2.0 * DETF**(-2.0/3.0) * c1 * ( ojo - 1.0/3.0 * trC * CINV )
    end if
    
  end subroutine cosseratincomp
  ! ======================================================================


  ! ======================================================================
  ! Complex version of the Cosserat Tensor function
  ! ======================================================================
  subroutine cosseratcomp ( Ecomp, FisPropElem , Scomp )

    ! --- modules ---
    use declarmodule
    use usefull
    use depsmodule

    implicit none

    integer                             :: elem , i, j, k , l    , ind   , it, FisPropElem            ! Indexes
  
    real*8                              :: ojo(3,3)

    real*8                              :: time
    double complex :: DETC
    double complex :: auxtracecomp, auxdetcomp, Ccomp(3,3), CcompINV(3,3)
    double complex , intent(in)                 :: Ecomp(3,3)
    double complex , intent(out)                :: Scomp(3,3)
    logical :: ok_flag
    
    double complex :: residmat(3,3)

!~     double complex  :: 

    ! identity
    ojo      = 0.0d+0
    ojo(1,1) = 1.0d+0
    ojo(2,2) = 1.0d+0
    ojo(3,3) = 1.0d+0
    !
    
    if (constitutive_model == 1 ) then ! svk model
      !
      young = MatsMat(FisPropElem,1)
      nu    = MatsMat(FisPropElem,2)
      !
      lambda  = young * nu / ( (1.0d+0 + nu) * (1.0d+0 - 2.0d+0 * nu) )
      shear   = young      / ( 2.0d+0 * (1.0d+0 + nu) )
      !
      call traceComplex ( Ecomp , auxtracecomp )
      Scomp = lambda * auxtracecomp * ojo  +  2.0d+0 * shear * Ecomp
      !
    elseif (constitutive_model == 2) then ! curnier model
      !
      lambda = MatsMat(FisPropElem,1)
      shear  = MatsMat(FisPropElem,2)
      !
      Ccomp = 2.0D+0 * Ecomp + ojo
      call M33INVComp ( Ccomp , CcompINV, OK_FLAG)
      DETC = M33DETComp(Ccomp)
!~       residmat = matmul(Ccomp, CcompINV)-ojo
!~       write(*,*) "C", Ccomp
!~       write(*,*) "resid mat:", residmat
!~       write(*,*) "detc:", DETC
      Scomp = lambda * ( sqrt(DETC) -1.0D+0 )* CcompINV +  2.0d+0 * shear * Ecomp
      !
    end if

  end subroutine cosseratcomp
  ! ======================================================================
  


  ! ======================================================================
  ! Complex version of the Cosserat Tensor function
  ! ======================================================================
  subroutine cosseratcompIncomp ( Ecomp , FisPropElem, p, Scomp )

    ! --- modules ---
    use declarmodule
    use usefull
    use depsmodule

    implicit none

    integer                             :: elem , i, j, k , l    , ind   , it   , FisPropElem          ! Indexes
  
    real*8                              :: ojo(3,3), c1
real*8 :: p
    real*8                              :: time
    double complex :: DETC, DETF
    double complex :: auxtracecomp, auxdetcomp, Ccomp(3,3), CcompINV(3,3)
    double complex , intent(in)                 :: Ecomp(3,3)
    double complex , intent(out)                :: Scomp(3,3)
    logical :: ok_flag
    
    double complex :: residmat(3,3)

    ! identity
    ojo      = 0.0d+0
    ojo(1,1) = 1.0d+0
    ojo(2,2) = 1.0d+0
    ojo(3,3) = 1.0d+0
    !
    

    if (constitutive_model == 6 ) then ! neohookean
      !
      c1 = MatsMat(FisPropElem,1)
      !
      Ccomp = 2.0D+0 * Ecomp + ojo
      call M33INVComp ( Ccomp , CcompINV, OK_FLAG )
      DETC = M33DETComp(Ccomp)
      DETF = sqrt(DETC)
      
      call traceComplex ( Ccomp , auxtracecomp )
      
      Scomp =  DETF * p * CcompINV  +  DETF**(-2.0/3.0) * 2* c1 * ( ojo - 1.0/3.0 * auxtracecomp * CcompINV )
    end if
!~     write(*,'(6e12.3)') Ccomp


  end subroutine cosseratcompIncomp
  ! ======================================================================


  
  ! ======================================================================
  ! Hessian of strain energy function
  ! ======================================================================
  subroutine ElasticityTensors ( E, FisPropElem )

    ! --- modules ---
    use declarmodule
    use usefull
    use depsmodule

    implicit none

    integer                             :: elem , i, j, k , l    , ind   , it, FisPropElem            ! Indexes
  
    real*8                              :: auxtrace, ojo(3,3)

    real*8                              :: time, direction(3,3) , deltaesc, dSdDir(3,3)
    real*8 , intent(in)                 :: E(3,3)

    deltaesc = 1.0d-10
    
    do i = 1, 6
      
      direction      = 0.0d+0
      if (i <4)  direction(i,i) = 1.0d+0
      if (i==4)  direction(2,3) = 1.0d+0
      if (i==5)  direction(1,3) = 1.0d+0
      if (i==6)  direction(1,2) = 1.0d+0
      
      call CompDerivTensor ( E, FisPropElem , direction , deltaesc , dSdDir )

      do j = 1, 3
        ConsMat (i , j ) = dSdDir(j,j)
      end do
      ConsMat (i , 4 ) = dSdDir(2,3) * 0.5d+0
      ConsMat (i , 5 ) = dSdDir(1,3) * 0.5d+0
      ConsMat (i , 6 ) = dSdDir(1,2) * 0.5d+0
      
    end do
    !
  end subroutine ElasticityTensors
  ! ======================================================================




  ! ======================================================================
  ! Hessian of strain energy function
  ! ======================================================================
  subroutine ElasticityTensorsIncomp ( E , FisPropElem, p)

    ! --- modules ---
    use declarmodule
    use usefull
    use depsmodule

    implicit none
        integer                                 :: FisPropElem

    integer                             :: elem , i, j, k , l    , ind   , it            ! Indexes
  
    real*8                              :: auxtrace, ojo(3,3)

    real*8                              :: time, direction(3,3) , deltaesc, dSdDir(3,3)
    real*8 , intent(in)                 :: E(3,3), p

    deltaesc = 1.0d-10
    
    do i = 1, 6
      
      direction      = 0.0d+0
      if (i <4)  direction(i,i) = 1.0d+0
      if (i==4)  direction(2,3) = 1.0d+0
      if (i==5)  direction(1,3) = 1.0d+0
      if (i==6)  direction(1,2) = 1.0d+0
      
      call CompDerivTensorIncomp ( E , FisPropElem, p, direction , deltaesc , dSdDir )

      do j = 1, 3
        ConsMat (i , j ) = dSdDir(j,j)
      end do
      ConsMat (i , 4 ) = dSdDir(2,3) * 0.5d+0
      ConsMat (i , 5 ) = dSdDir(1,3) * 0.5d+0
      ConsMat (i , 6 ) = dSdDir(1,2) * 0.5d+0
      
    end do
    !
  end subroutine ElasticityTensorsIncomp
  ! ======================================================================


  ! ======================================================================
  ! Complex Step Directional derivative of tensor
  ! ======================================================================
  subroutine CompDerivTensor ( Ereal, FisPropElem, Dir, deltaesc, dfcsm )
    !
    implicit none
    !
    integer                                 :: FisPropElem
    !
    real*8     ,              intent(in)    :: Ereal(:,:) , Dir(:,:) , deltaesc
    real*8     ,              intent(out)   :: dfcsm(3,3)
    double complex                          :: Ecomp(3,3), Scomp(3,3)
    !

    Ecomp = Ereal * ( 1.0d+0 , 0.0d+0 )  + deltaesc * Dir * ( 0.0d+0 , 1.0d+0 )

    call cosseratcomp ( Ecomp, FisPropElem , Scomp )
    
    dfcsm = imagpart( Scomp ) / deltaesc

  end subroutine CompDerivTensor
  ! ======================================================================


  ! ======================================================================
  ! Complex Step Directional derivative of tensor
  ! ======================================================================
  subroutine CompDerivTensorIncomp ( Ereal, FisPropElem, p, Dir, deltaesc, dfcsm )
    !
    implicit none
    !
        integer                                 :: FisPropElem

    real*8     ,              intent(in)    :: Ereal(:,:) , Dir(:,:) , deltaesc, p
    real*8     ,              intent(out)   :: dfcsm(3,3)
    double complex                          :: Ecomp(3,3), Scomp(3,3)
    !

    Ecomp = Ereal * ( 1.0d+0 , 0.0d+0 )  + deltaesc * Dir * ( 0.0d+0 , 1.0d+0 )

    call cosseratcompIncomp ( Ecomp , FisPropElem, p , Scomp )
    
    dfcsm = imagpart( Scomp ) / deltaesc

  end subroutine CompDerivTensorIncomp
  ! ======================================================================



end module materialmodule
