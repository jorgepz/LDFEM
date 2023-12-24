! ==============================================================================
! -------    LDFEM: a Large Deformations Finite Element Method solver   --------
!
! version history:
! 0.1.7 : 2016.09.01 jpz cambios varios
! 0.1.6 : 2015.11.12 jpz
! 0.1.5 : 2015.01.23 jpz
! 0.1.4 : 2014.12.01 pcg
! 0.1.3 : 2014.08.24 jpz
! 0.1.2 : 2014.08.12 jpz
! 0.1.1 : 2014.08.11 jpz
!
! ==============================================================================

program ldfem

  ! ======================================================================
  ! --- solver modules  ---
  use declarmodule
  use vtkmodule
  use elemmodule
  use boundarymodule
  use materialmodule
  use usefull

  ! ======================================================================


  ! ======================================================================
  implicit none

  integer                             :: elem , i, j, k , l    , ind   , it,cn, aux3            ! Indexes
  !
  character(len=20)                   :: auxchar
  character(len=20)                   :: meshfilename, probname
  integer                             :: auxint
  integer                             :: auxsurf(3)
  real*8                              :: auxarea 
  !
  real*8                              :: dx1, dx2, dx3, L1, L2, L3, despnodo(1:3), H(1:3,1:3), a, b      , xi, wi
  !
  real*8                              :: F(3,3), invF(3,3), detF, E(3,3), S(3,3), P(3,3), Sig(3,3)
  
  integer                             :: counter, elemtype
  !
  integer                             :: latexflag, geomflag, condefflag , facei1isfix
  !
  integer                             :: elemind(NNodPE*3)
  real*8                              :: aux1, aux2, sumauxdet, auxdet1, auxdet2, auxdet3, auxdet4, tolauxdet
  real*8        ,allocatable          :: sigmas(:,:), resul(:) , resul2(:)
  !
  real*8                              :: Be(6,12), Ue(12),  matauxgeom(4,4)       , otramataux(12,12)
  real*8                              :: pe

  real*8                              :: voigCINV(6)
  
  integer                             :: linsolverflag

  integer                             :: tolits
  real*8                              :: tolgrad, normgrad, toldeltagrad, normrighthand

  real*8                              :: auxreal1, auxreal2        , ojo(3,3)

  real*8                              :: normgradhist(100), normresidload

  real*8        , allocatable         :: fun(:,:), deriv(:,:), funder(:,:)
  
  real*8        , allocatable         :: Kml(:,:), Kmg(:,:), elecoordmat(:,:), elecoordspa(:,:), Kgl(:,:) , Kgg(:,:)

  ! tensors calculation for compressible case
  real*8        , allocatable         :: AllCosMat(:,:),AllPioMat(:,:), AllSigMat(:,:), sigVM_Cau(:,:)

  real*8                              :: EpsMat(6,12) , BgrandeMat(6,12)   , voigS(6) , voigP(6) , voigSig(6)
  integer                             :: nodeselem(4), dofselem(12)
  real*8                              :: jacobianmat (3,3) , invjacobianmat (3,3), jacobiandet , DeeMat(6,6)
  real*8                              :: jacobianmatspa (3,3) , invjacobianmatspa (3,3), jacobiandetspa
  logical                             :: okflag, ex
  
  real*8 :: DispTotElem(12) ! disp of nodes from 1 elem   
  real*8 :: ConstMatPerturbed(6,6)
  
  real*8 :: detC, Cinv(3,3), C(3,3)
  
  real*8 , allocatable :: ones(:)
  real*8 :: kappa
  
  integer                             :: ok
  integer, allocatable                :: pivot(:)


  real*8 :: checkmat(4,4)

  real*8 :: ConsMatAnali(6,6)

  ! control node
  integer, allocatable                 :: ControlNode(:)
  real*8, allocatable                  :: ControlNodeCoord(:,:)
  real*8, allocatable                  :: DistancesToControlNode(:)
  integer, allocatable                 :: ControlElem(:)
  real*8, allocatable                  :: ControlElemCoord(:,:), ControlElemPesos(:,:)

  integer :: compresibleflag   ! 1 si es compresible 2 si es incompresible
  
  integer             :: indtime, Ntimes
  real*8,allocatable  :: times(:)
  real*8              :: deltatime, factorini, factorfin
  real :: timesolver1, timesolver2
  
  integer :: contadormaxHess
  ! ---------------------------------------------------------------------------------


  ! ======================================================================
  
  
  ! ======================================================================
  ! =========  0. Folders check  =========================================
  ! ======================================================================
  
  write(*,'(A)') "===================================="
  write(*,'(A)') "-------   Welcome to LDFEM   -------"

  call getarg(1,probname)

  INQUIRE(FILE = "output/", EXIST = ex)

  if ( .not. ex ) call system( "mkdir output" )
  
  INQUIRE(FILE = "output/"//trim(probname), EXIST = ex)

  if ( .not. ex ) then ! si no existe la crea
    call system( "mkdir output/" //trim(probname))
  else                 ! si existe la borra y la crea
    call system( "rm -R output/" //trim(probname))
    call system( "mkdir output/" //trim(probname))
  end if
  
  INQUIRE(FILE = "examples", EXIST = ex)

  if ( .not. ex ) write(*,*) "FALTA CARPETA input"
  
  INQUIRE(FILE = "examples/"//trim(probname), EXIST = ex)

  if ( .not. ex )  write(*,*) "FALTA CARPETA "//trim(probname)//" DENTRO DE LA CARPETA input"

  ! ======================================================================

  

  ! ======================================================================
  ! =========  1. Parameters Def/Reading  ================================
  ! ======================================================================

  ! ----   iteration parameters   ----------------------------
  tolits       = 20
  tolgrad      = 1.0d-8
  toldeltagrad = 1.0d-8
 
  ! --- quasicompressible formulation penalization factor
  kappa = 1000000.0
 
  ! identity ( equivalent to eye(3) )
  ojo = 0.0d+0
  ojo(1,1) = 1.0d+0
  ojo(2,2) = 1.0d+0
  ojo(3,3) = 1.0d+0
 
  ! ---------------------------------------------------------
  ! --------   Input parameters reading  --------------------
  
  call chdir("examples/" //trim(probname))

  open (unit = 11 , file = trim(probname)//".inp" )

  read (11,*) auxchar
  if ( auxchar /= '0.1.7' ) then
    write(*,'(A)') "Error version in input"
    stop
  end if

  read (11,*) meshfilename
  
  
  ! materials definitions
  read (11,*) constitutive_model, linsolverflag
  
  if ( constitutive_model > 5) then
    compresibleflag = 2
  else
    compresibleflag = 1
  end if

  ! time factors
  read (11,*) factorini, factorfin, Ntimes
  allocate( times(Ntimes) )
  if (Ntimes == 1) then
    write(*,*) "debe haber al menos dos tiempos"
    stop
  end if
  deltatime = (factorfin - factorini) / dble(Ntimes-1)
  call linspace( times , deltatime, deltatime, Ntimes )
  times(1) = factorini
  write(*,*) "Los factores de carga a resolver son: ",times
  ! ----------------

  ! control node
  read (11,*) ncontrolnodes
  
  allocate( ControlNodeCoord(ncontrolnodes,3) )
  
  do i =1, ncontrolnodes
    read (11,*) ( ControlNodeCoord(i,j), j=1 , 3 )
  end do
  ! ----------------

  ! materials definitions
  read (11,*) nmats

  allocate( MatsMat(nmats,4) )

  do i = 1, nmats
    read (11,*) auxreal1, ( MatsMat(i,j), j=1 , 4 )
  end do
  ! ----------------

  ! boundary conditions
  read (11,*) nBC

  allocate( BCDefMat(nBC,7) )

  BCDefMat = 0.0
  
  do i = 1, nBC
    read (11,*) auxreal1, ( BCDefMat(i,j), j=1 , 7 )
  end do
  ! ----------------
  
  
  
  ! scale factors
  read (11,*) nscalefactors

  allocate( ScaleFactors(nscalefactors) )

  read (11,*) ( ScaleFactors(j), j=1 , nscalefactors )
  ! ----------------
    

  close(11)
  ! -------  finished parameters reading -------------------


  ! -------  mesh file reading -----------------------------
  open (unit = 11 , file = meshfilename )
  !
  write(*,*) meshfilename
  read (11,*) nnod, ntet, ntri
  write(*,'(A,/,A,i6,/,A,i6,/,A,i6)') "Totales:", "  nodos:", nnod, "  tetra:", ntet, "  trian:", ntri
  !
  allocate( nodes(nnod,3) , Cellsconec(ntet,NNodPE) , Surfconec(ntri,3))
  nodes        = 0.0
  Cellsconec   = 0
  Surfconec    = 0
  !
  allocate( NodBCMat( nnod) , CellsBCMat( ntet ) , SurfBCMat(ntri), VolsMats(ntet) )
  !
  NodBCMat     = 0
  CellsBCMat   = 0
  SurfBCMat    = 0

  do i = 1, nnod
    read (11,*) ( nodes(i,j), j=1 , 3 ), auxreal1
    NodBCMat(i) = int( auxreal1 )
  end do
  !
  do i = 1, ntet
    read (11,*) Cellsconec(i,1), Cellsconec(i,2), Cellsconec(i,3), Cellsconec(i,4), auxreal1
    VolsMats(i) = int( auxreal1 )
  end do
  !
  do i = 1, ntri
    read (11,*) ( Surfconec(i,j),j=1,3 ), auxreal1
    SurfBCMat(i) = int( auxreal1 )
  end do
  !
  close(11)

  allocate( ones( 3*nnod ) )
  ones(:) = [ (1.0,i=1,3*nnod ) ]

  ! --------------------------------------------------
  ! Calcula cuales son los nodos mas cercanos a cada nodo de control insertado
  ! las coordenadas originales quedan guardadas en ControlElemCoord

  allocate( DistancesToControlNode(nnod) , ControlNode(ncontrolnodes) )

  write(*,'(A)')    "-------------------------------"
  do j=1, ncontrolnodes
    DistancesToControlNode =  0.0

    do i = 1, 3
!~       DistancesToControlNode(:) = DistancesToControlNode(:)  + &
!~         abs( real ( nodes(:,i)  ) - ControlNodeCoord(j,i) )
      DistancesToControlNode(:) = DistancesToControlNode(:)  + &
        abs( nodes(:,i) - ControlNodeCoord(j,i) )
    end do
    
    ControlNode(j) = minloc( DistancesToControlNode , 1)
    write(*,'(A,i10,A,i10)')    "Control node number ",j,"/", ncontrolnodes
    write(*,'(A,i10,A,3e12.3)') "Closest node to coordinates entered is: ", ControlNode(j), &
      " distance: ", DistancesToControlNode( ControlNode(j))
    write(*,'(A,3e12.3)') "Coordinates: ", nodes(ControlNode(j),:)
  end do
!~   write(*,'(A)')    "-------------------------------"
  ! --------------------------------------------------
  
  

  ! --------------------------------------------------
  ! calcula cuales son los elementos de control mas cercanos a cada nodo de contro
  ! utliza las coordenadas originales de cada nodo de control
  
  allocate( ControlElem(ncontrolnodes) ,  ControlElemPesos(ncontrolnodes,4) )
  ControlElem      = 0      
  ControlElemPesos = 0      
  
  tolauxdet   = 1.0D-6 ! tolerancia para volumenes relativos de tetraedros
  
  do cn = 1, ncontrolnodes
    i = 1
    aux3 = 1
    do while (i <= ntet .and. aux3 == 1)
      j = 1
      do while (j <= 4 .and. aux3 == 1)
        if ( Cellsconec(i,j) == ControlNode(cn) ) then          
          auxdet1 = M44DET( [ [1.0D+0,      ControlNodeCoord(cn , :)], &
                              [1.0D+0, nodes( Cellsconec(i,2) , : )], &
                              [1.0D+0, nodes( Cellsconec(i,3) , : )], &
                              [1.0D+0, nodes( Cellsconec(i,4) , : )]] )
                              
          auxdet2 = M44DET( [ [1.0D+0, nodes( Cellsconec(i,1) , : )], &
                              [1.0D+0,      ControlNodeCoord(cn , :)], &
                              [1.0D+0, nodes( Cellsconec(i,3) , : )], &
                              [1.0D+0, nodes( Cellsconec(i,4) , : )]] )
                              
          auxdet3 = M44DET( [ [1.0D+0, nodes( Cellsconec(i,1) , : )], &
                              [1.0D+0, nodes( Cellsconec(i,2) , : )], &
                              [1.0D+0,      ControlNodeCoord(cn , :)], &
                              [1.0D+0, nodes( Cellsconec(i,4) , : )]] )
                              
          auxdet4 = M44DET( [ [1.0D+0, nodes( Cellsconec(i,1) , : )], &
                              [1.0D+0, nodes( Cellsconec(i,2) , : )], &
                              [1.0D+0, nodes( Cellsconec(i,3) , : )], &
                              [1.0D+0,      ControlNodeCoord(cn , :)]] )
                              
          sumauxdet = auxdet1+auxdet2+auxdet3+auxdet4
!~           write(*,*) auxdet1,auxdet2,auxdet3,auxdet4
!~           write(*,*) sumauxdet
          if  ( ( auxdet1/sumauxdet >= 0 .and. auxdet2/sumauxdet >= 0 .and. &
                  auxdet3/sumauxdet >= 0 .and. auxdet4/sumauxdet >= 0 ) &
                .or. &
                ( auxdet1/sumauxdet >= -tolauxdet .or. auxdet2/sumauxdet >= -tolauxdet .or. &
                  auxdet3/sumauxdet >= -tolauxdet .or. auxdet4/sumauxdet >= -tolauxdet        ) ) then
            ControlElem(cn) = i
            aux3 = 0
            ControlElemPesos(cn,:) = [auxdet1/sumauxdet, auxdet2/sumauxdet, auxdet3/sumauxdet, auxdet4/sumauxdet]
          end if
        end if
        j = j+1
      end do
      i = i+1
    end do
    
    if (cn == 315 .or. cn==314 ) then
      write(*,*) aux3, cn, ControlElem(cn)
      if ( cn ==315) then
!~         ControlNodeCoord(k,:)
!~         stop
      end if
    end if

  end do
  ! --------------------------------------------------
   
  deallocate(DistancesToControlNode)

  ! ------------------------------------------  
  ! ======================================================================
  
  

  ! ======================================================================
  ! =========     2. Pre-process       ===================================
  ! ======================================================================

  allocate ( Fext(nnod*3) , Udiri(nnod*3) , AuxNeumDOFs(3*nnod) , Udirikm1(nnod*3 ) )

  Fext        = 0.0D+0
  Udiri       = 0.0D+0
  Udirikm1    = 0.0D+0
  AuxNeumDOFs = [(i,i=1,3*nnod)] 

  allocate ( NeumDOFs     (nnod*3)  , DiriDOFs     (nnod*3)  )
  allocate ( NeumSortDOFs (nnod*3)  , DiriSortDOFs (nnod*3)  )
  !
  NNeumDOFs    = 0
  NeumDOFs     = 0
  NeumSortDOFs = 0
  !
  NDiriDOFs    = 0
  DiriDOFs     = 0
  DiriSortDOFs = 0

  allocate( DispBC(3*nnod)   )
  DispBC           = 0.0D+0

  allocate( ldnodesdef(nnod,3) )
 
  ! coord nodes large deformations
  ldnodesdef = nodes

  ! ----- Boundary Conditions SubRoutine ------
  call boundarycond ( times( 1) )
  ! -------------------------------------------

  NnonzeroDiriDOFs = 0
  
  ! recorro los auxneum dofs para ver donde tengo diris impuestos
  do i = 1, 3*nnod
    if (AuxNeumDOFs(i) == 0 ) then !it is diri
      !
      NDiriDOFs           = NDiriDOFs + 1
      DiriDOFs(NDiriDOFs) = i
      DiriSortDOFs(i)     = NDiriDOFs
      !
      ! cuento los diris no cero
      if ( Udiri( i ) /= 0.0 ) NnonzeroDiriDOFs = NnonzeroDiriDOFs + 1
      !
        
    else
      ! 
      NNeumDOFs = NNeumDOFs + 1
      NeumDOFs(NNeumDOFs) = i  
      NeumSortDOFs(i)     = NNeumDOFs
      !
    end if
  end do
  

  write(*,*) "boundary condtitions:", " nneumdofs:", NNeumDOFs, " Ndiridofs:", NDiriDOFs


  if ( NnonzeroDiriDOFs > 0 ) then
    write(*,*) "Nnonzero diridofs:", NnonzeroDiriDOFs
    ! ------      ---------------
    allocate( NonzeroDiriDOFs ( NnonzeroDiriDOFs ) , NonzeroDiriVALs ( NnonzeroDiriDOFs ) ) 
    NnonzeroDiriDOFs = 0
    NonzeroDiriDOFs = 0
    NonzeroDiriVALs = 0.0D+0

    do i = 1, NDiriDOFs
      if ( Udiri( DiriDOFs(i) ) /= 0.0 ) then
        NnonzeroDiriDOFs = NnonzeroDiriDOFs + 1
        NonzeroDiriDOFs( NnonzeroDiriDOFs )  = DiriDOFs(i)
        NonzeroDiriVALs( NnonzeroDiriDOFs )  = Udiri( DiriDOFs(i) )
      end if
    end do
    ! ---------------------
  end if !  --- ( NnonzeroDiriDOFs > 0 ) ---


  
  DispBC           = 0.0D+0
  DispBC = Udiri    

  allocate( BCMat(nnod,6) )
  !  '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
  do i = 1, nnod
    auxsurf = [3*i-2 , 3*i-1 , 3*i ]
    BCMat(i,[1,2,3]) = DispBC(auxsurf)
    BCMat(i,[4,5,6]) = Fext  (auxsurf)
  end do
  
  call chdir("../../output/" //trim(probname))
  
  ! --- vtk generation ---
  call vtkgen &
    ( "inputvtk" ,"nodes_bc", "cellsind" &
    , NNodPE , 3   &
    , reshape( [ BCMat(:,1) , BCMat(:,2) , BCMat(:,3) , BCMat(:,4) , BCMat(:,5) , BCMat(:,6)  ] , [nnod,6] ) &
    , reshape( [ (dble(j),j=1,ntet) ] ,[ntet,1]) &
    , Cellsconec(1:ntet,1:4) , nodes(1:nnod,1:3) )
  ! ----------------------

  INQUIRE(FILE = "LD_"//trim(probname) , EXIST = ex)
  
  if ( .not. ex ) then
    !
    call system( "mkdir LD_"//trim(probname) )
    !
  end if

  call chdir("LD_"//trim(probname))
  ! ======================================================================


  ! ======================================================================
  ! =========  3. Process              ===================================
  ! ======================================================================

  allocate ( GradMat ( nnod*3 ) , ElemVols( ntet ) )

  !
  latexflag  = 0
  condefflag = 0
  geomflag   = 1 ! 1 cube geometry ; 2 pipe
  elemtype   = 1 ! 1 tethraedron   ; 2 cube
  ! --------------------------------------------------------------------

  if ( compresibleflag == 2 ) then
    allocate (  pTot(ntet)  )
    pTot = 0.0d+0
  end if


  ! ==================== LOOP ===================
  !
  allocate( DisplacNeum (NNeumDOFs) , DisplacTot (nnod*3) )
  allocate( pivot(NNeumDOFs) )
  allocate( resul(NNeumDOFs) )
  allocate( resul2(3*nnod) )
  allocate( MatDisplacTot(nnod,3)) !, SumaMatDisplacTot(nnod,3) )

  
!~ ####################################################################    
!~ ####################################################################     
!~     Comienzo edición de Pablo - Cálculo de tensores para el caso compresible 03-11-2014
!~ ####################################################################
!~ ####################################################################
  
  INQUIRE(FILE = "Stress/", EXIST = ex)

  if ( .not. ex ) call system( "mkdir Stress" )
  
  INQUIRE(FILE = "Stress/Cosserat", EXIST = ex)

  if ( .not. ex ) call system( "mkdir Stress/Cosserat" )
  
  INQUIRE(FILE = "Stress/Piola", EXIST = ex)

  if ( .not. ex ) call system( "mkdir Stress/Piola" )
  
  INQUIRE(FILE = "Stress/Cauchy", EXIST = ex)

  if ( .not. ex ) call system( "mkdir Stress/Cauchy" )
  
  INQUIRE(FILE = "Stress/Von_Misses", EXIST = ex)

  if ( .not. ex ) call system( "mkdir Stress/Von_Misses" )
  
  INQUIRE(FILE = "Deformed/", EXIST = ex)

  if ( .not. ex ) call system( "mkdir Deformed" )
  
  INQUIRE(FILE = "Deformed/Step_by_step/", EXIST = ex)

  if ( .not. ex ) call system( "mkdir Deformed/Step_by_step" )
  
  allocate( AllCosMat(6,ntet), AllPioMat(9,ntet), AllSigMat(6,ntet), sigVM_Cau(ntet,1) )
  

  call chdir("Deformed")
    open (unit = 14 , file = trim(probname)//"_ControlNode_Def_Position.dat" )
  call chdir("../")
  
!~   if (NNonZeroDiriDOFs>0) then
    open (unit = 15 , file = trim(probname)//"_Grad_Material.dat" )
    write(15,'(A)') "||                    Load|           Material Grad||"
    write(15,'(A)') "-----------------------------------------------------"
!~   end if
  
  call chdir("Stress")
    open (unit = 16 , file = trim(probname)//"_ControlElem_Cauchy.dat" )
    open (unit = 17 , file = trim(probname)//"_ControlElem_Cosserat.dat" )
    open (unit = 18 , file = trim(probname)//"_ControlElem_Piola.dat" )
    open (unit = 19 , file = trim(probname)//"_ControlElem_Von_Misses.dat" )
  call chdir("../")
  
  do j=1, ncontrolnodes
    aux3=ControlNode(j)
    if (j<ncontrolnodes) then
!~       write(14,'(A,i10,A,i10,A)',advance='no') "||                             Control Node", j &
!~       ,"(",aux3,")                           ||   "
    else
!~       write(14,'(A,i10,A,i10,A)') "||                             Control Node", j,"(",aux3,")                           ||"
    end if
  end do
  do j=1, ncontrolnodes
    if (j<ncontrolnodes) then
!~       write(14,'(A,A,A,A)',advance='no') "-------------------" , "-------------------------", "-------------------------"&
!~       , "-------------------------   "
    else
!~       write(14,'(A,A,A,A)') "-------------------" , "-------------------------", "-------------------------"&
!~       , "-------------------------"
    end if
  end do
  do j=1, ncontrolnodes
    if (j<ncontrolnodes) then
!~       write(14,'(A,A,A,A)',advance='no') "||           Load|" , "                       X|", "                       Y|"&
!~       , "                       Z||   "
    else
!~       write(14,'(A,A,A,A)') "||           Load|" , "                       X|", "                       Y|"&
!~       , "                       Z||"
    end if
  end do
  do j=1, ncontrolnodes
    if (j<ncontrolnodes) then
!~       write(14,'(A,A,A,A)',advance='no') "-------------------" , "-------------------------", "-------------------------"&
!~       , "-------------------------   "
    else
!~       write(14,'(A,A,A,A)') "-------------------" , "-------------------------", "-------------------------"&
!~       , "-------------------------"
    end if
  end do


    do cn=1, ncontrolnodes
      aux3 = ControlElem(cn)

      write(14,'(e10.3,i8,e28.15,e28.15,e28.15)') 0.0 , cn &
      , dot_product(ldnodesdef(Cellsconec(aux3,:),1), ControlElemPesos(cn,:) ) &
      , dot_product(ldnodesdef(Cellsconec(aux3,:),2), ControlElemPesos(cn,:) ) &
      , dot_product(ldnodesdef(Cellsconec(aux3,:),3), ControlElemPesos(cn,:) )

    end do

  
  do j=1, ncontrolnodes
    aux3 = ControlElem(j)
    if (j<ncontrolnodes) then
      write(16,'(A,i10,A,i10,A)',advance='no') "||                                                              Control Element"&
      , j,"(",aux3 ,")                                                                              ||   "
    else
      write(16,'(A,i10,A,i10,A)') "||                                                              Control Element"&
      , j,"(",aux3 ,")                                                                              ||"
    end if
  end do
  do j=1, ncontrolnodes
    if (j<ncontrolnodes) then
      write(16,'(A,A)',advance='no') "--------------------------------------------------------------------"&
,"-----------------------------------------------------------------------------------------------------------------   "
    else
      write(16,'(A,A)') "--------------------------------------------------------------------"&
,"-----------------------------------------------------------------------------------------------------------------"
    end if
  end do
  do j=1, ncontrolnodes
    if (j<ncontrolnodes) then
      write(16,'(A,A,A,A,A,A,A)',advance='no') "||           Load|" , "                        XX|", "                        YY|"&
  ,"                        ZZ|", "                        XY|", "                        YZ|", "                        XZ||   "
    else
      write(16,'(A,A,A,A,A,A,A)') "||           Load|" , "                        XX|", "                        YY|"&
  ,"                        ZZ|", "                        XY|", "                        YZ|", "                        XZ||"
    end if
  end do
  do j=1, ncontrolnodes
    if (j<ncontrolnodes) then
      write(16,'(A,A)',advance='no') "--------------------------------------------------------------------"&
,"-----------------------------------------------------------------------------------------------------------------   "
    else
      write(16,'(A,A)') "--------------------------------------------------------------------"&
,"-----------------------------------------------------------------------------------------------------------------"
    end if
  end do
  
  do j=1, ncontrolnodes
    aux3 = ControlElem(j)
    if (j<ncontrolnodes) then
      write(17,'(A,i10,A,i10,A)',advance='no') "||                                                              Control Element"&
      , j,"(",aux3 ,")                                                                              ||   "
    else
      write(17,'(A,i10,A,i10,A)') "||                                                              Control Element"&
      , j,"(",aux3 ,")                                                                              ||"
    end if
  end do
  do j=1, ncontrolnodes
    if (j<ncontrolnodes) then
      write(17,'(A,A)',advance='no') "--------------------------------------------------------------------"&
,"-----------------------------------------------------------------------------------------------------------------   "
    else
      write(17,'(A,A)') "--------------------------------------------------------------------"&
,"-----------------------------------------------------------------------------------------------------------------"
    end if
  end do
  do j=1, ncontrolnodes
    if (j<ncontrolnodes) then
      write(17,'(A,A,A,A,A,A,A)',advance='no') "||           Load|" , "                        XX|", "                        YY|"&
  ,"                        ZZ|", "                        XY|", "                        YZ|", "                        XZ||   "
    else
      write(17,'(A,A,A,A,A,A,A)') "||           Load|" , "                        XX|", "                        YY|"&
  ,"                        ZZ|", "                        XY|", "                        YZ|", "                        XZ||"
    end if
  end do
  do j=1, ncontrolnodes
    if (j<ncontrolnodes) then
      write(17,'(A,A)',advance='no') "--------------------------------------------------------------------"&
,"-----------------------------------------------------------------------------------------------------------------   "
    else
      write(17,'(A,A)') "--------------------------------------------------------------------"&
,"-----------------------------------------------------------------------------------------------------------------"
    end if
  end do

  
  do j=1, ncontrolnodes
    aux3 = ControlElem(j)
    if (j<ncontrolnodes) then
      write(18,'(A,A,i10,A,i10,A)',advance='no') "||                                                                         "&
      ,"                                         Control Element"&
, j,"(",aux3 ,")                                                                                                           ||   "
    else
      write(18,'(A,A,i10,A,i10,A)') "||                                                                                      "&
      ,"                            Control Element"&
  , j,"(",aux3 ,")                                                                                                           ||"
    end if
  end do
  do j=1, ncontrolnodes
    if (j<ncontrolnodes) then
      write(18,'(A,A,A)',advance='no') "--------------------------------------------------------------------------------------"&
       ,"----------------------------------------------------------------------------------------------------------------------"&
       ,"----------------------------------------------------------   "
    else
      write(18,'(A,A,A)') "--------------------------------------------------------------------------------------"&
       ,"----------------------------------------------------------------------------------------------------------------------"&
       ,"----------------------------------------------------------"
    end if
  end do
  do j=1, ncontrolnodes
    if (j<ncontrolnodes) then
      write(18,'(A,A,A,A,A,A,A,A,A,A)',advance='no') "||           Load|" , "                     XX-P0|"&
      , "                     YY-P1|", "                     ZZ-P2|", "                     XY-P3|"&
      , "                     YZ-P4|", "                     XZ-P5|", "                     YX-P6|"&
      , "                     ZY-P7|", "                     ZX-P8||   "
    else
      write(18,'(A,A,A,A,A,A,A,A,A,A)') "||           Load|" , "                     XX-P0|"&
      , "                     YY-P1|", "                     ZZ-P2|", "                     XY-P3|"&
      , "                     YZ-P4|", "                     XZ-P5|", "                     YX-P6|"&
      , "                     ZY-P7|", "                     ZX-P8||"
    end if
  end do
  do j=1, ncontrolnodes
    if (j<ncontrolnodes) then
      write(18,'(A,A,A)',advance='no') "--------------------------------------------------------------------------------------"&
       ,"----------------------------------------------------------------------------------------------------------------------"&
       ,"----------------------------------------------------------   "
    else
      write(18,'(A,A,A)') "--------------------------------------------------------------------------------------"&
       ,"----------------------------------------------------------------------------------------------------------------------"&
       ,"----------------------------------------------------------"
    end if
  end do
  
  do j=1, ncontrolnodes
    aux3 = ControlElem(j)
    if (j<ncontrolnodes) then
      write(19,'(A,i10,A,i10,A)',advance='no') "||  Control Element", j,"(", aux3 ,")   ||   " 
    else
      write(19,'(A,i10,A,i10,A)') "||  Control Element", j,"(", aux3 ,")   ||   "
    end if
  end do
  do j=1, ncontrolnodes
    if (j<ncontrolnodes) then
      write(19,'(A,A)',advance='no') "------------------" , "----------------------------   "
    else
      write(19,'(A,A)') "------------------" , "----------------------------"
    end if
  end do
  do j=1, ncontrolnodes
    if (j<ncontrolnodes) then
      write(19,'(A,A)',advance='no') "||           Load|", "                        VM||   "
    else
      write(19,'(A,A)') "||           Load|", "                        VM||   "
    end if
  end do
  do j=1, ncontrolnodes
    if (j<ncontrolnodes) then
      write(19,'(A,A)',advance='no') "------------------" , "----------------------------   "
    else
      write(19,'(A,A)') "------------------" , "----------------------------"
    end if
  end do
  
  contadormaxHess = 13 * 6 * ntet
  allocate( HessMatSparse( contadormaxHess ) , HessMatIndexes( contadormaxHess ,2) )

  if (linsolverflag == 1 ) allocate( Kmg( 3*nnod , 3*nnod ) )
  
  if ( NnonzeroDiriDOFs > 0 ) then
    allocate( HessMatSparseNeumNonzeroDiri ( contadormaxHess ) &
      , HessMatIndexesNeumNonzeroDiri ( contadormaxHess ,2) )
  end if

!~ ####################################################################    
!~ ####################################################################     
!~     Fin edición de Pablo - Cálculo de tensores para el caso compresible 03-11-2014
!~ ####################################################################
!~ ####################################################################

  allocate( Kml(12,12) , Kgl(12,12) , elecoordspa(3,4) , elecoordmat(3,4)  , funder(3,4) )

  write(*,'(7f8.2)') transpose(BCDefMat)
  Fext  = 0.0d+0
    








  ! ===================================================
  ! ------------ Time loop ---------------
  times_do: do indtime = 1 , Ntimes
  
    write ( * , '(/,A)' ) "-------  New time step: -------"
    write ( * , '(A,i4,A,e10.2,A,e10.2,A)' ) " indtime: ", indtime , &
      " factor: ", sum( times(1:indtime) ), " (delta: ", deltatime, ")"
    it            = 0
    normgrad      = 1.0d+10
    normrighthand = 1.0d+0

    call boundarycond ( times( indtime) )

    if ( NnonzeroDiriDOFs > 0 ) then
      NonzeroDiriVALs = 0.0D+0
      NnonzeroDiriDOFs = 0
      !
      do i = 1, NDiriDOFs
        if ( Udiri( DiriDOFs(i) ) /= 0.0 ) then
          NnonzeroDiriDOFs = NnonzeroDiriDOFs + 1
          NonzeroDiriVALs( NnonzeroDiriDOFs )  = Udiri( DiriDOFs(i) )
        end if
      end do
    end if !  --- ( NnonzeroDiriDOFs > 0 ) ---
  
    DispBC = Udiri 
    Udirikm1 = Udiri 
  
  
    call cpu_time(time1)
  
    ! ===================================================
    ! ============    Newton Raphson    =================
    !
    ! -------------  starts central loop -----------------------------
    write ( * , '(A)' ) "--- Starts Newton-Raphson -------- "
    write ( * , '(A)' ) "#it |    gradient norm     |  residual load norm  |"
    write ( * , '(A)' ) "---------------------------------------------------"
    !
    newtonraphson: do while ( normgrad/normrighthand > tolgrad  .and. it < tolits )
      !
      it = it +1
      contador = 0
      contadorHNeNzDi = 0
      GradMat        = 0.0D+0
      HessMatSparse  = 0.0D+0
      HessMatIndexes = 0

      if (NnonzeroDiriDOFs > 0 ) then
        HessMatSparseNeumNonzeroDiri  = 0.0D+0
        HessMatIndexesNeumNonzeroDiri = 0
      end if
      
      ! ===========================================================
      ! ===========   Grad hess assemb  ===========================
      !
      elements_do: do elem = 1, ntet
      
        nodeselem =  Cellsconec( elem ,1:4)

        call nods2dofs( nodeselem , 4 , 3 , dofselem )
      
!~         write(*,*) "elem:", elem , "dofs:" , dofselem
      
        Kml = 0.0d+0
        Kgl = 0.0d+0
        
        elecoordmat(1,1:4) = nodes( nodeselem , 1 )
        elecoordmat(2,1:4) = nodes( nodeselem , 2 )
        elecoordmat(3,1:4) = nodes( nodeselem , 3 )
  
        elecoordspa(1,1:4) = ldnodesdef( nodeselem , 1 )
        elecoordspa(2,1:4) = ldnodesdef( nodeselem , 2 )
        elecoordspa(3,1:4) = ldnodesdef( nodeselem , 3 )
      
        xi = 0.25d+0
        wi = 1.0d+0/6.0d+0
      
        call ShapeFun( 1 , 0 , xi, xi , xi , fun )
        call ShapeFun( 1 , 1 , xi, xi , xi , deriv )

        ! jacobiana
        jacobianmat    = matmul ( elecoordmat , transpose( deriv) )
        jacobianmatspa = matmul ( elecoordspa , transpose( deriv) )
        
        ! material coord
        call M33INV ( jacobianmat    , invjacobianmat    , okflag , jacobiandet    )
        ! spatial coord
        call M33INV ( jacobianmatspa , invjacobianmatspa , okflag , jacobiandetspa )
        
        ! --------------------------------------
        if ( jacobiandet <= 0.0 ) then
          write(*,*) "Error volumen ",i, "negativoo!! ", jacobiandet
          write(*,*) "hay como: ", Ntet, "elementos en total"
          stop
        end if
  
        funder = matmul ( transpose ( invjacobianmat ) , deriv )
  
        
        ! --- disp grad ---
        H = matmul( (elecoordspa - elecoordmat )  , transpose ( funder ) )
  
        ! --- deformation grad ---
        F = ojo + H

        ! --- Lagrange tensor ---
        E = 0.5D+0 * ( H + transpose(H) + matmul( transpose(H) , H ) )

        C = 2.0D+0 * E + ojo

        if ( compresibleflag == 1 ) then

          ! --- cosserat tensor ---
          call cosserat ( E, VolsMats(elem) , S )
  
          call ElasticityTensors(E, VolsMats(elem))
  
  
          
          call  BgrandeMats ( BgrandeMat , funder , F )
  
        elseif (compresibleflag == 2 ) then

          pe = pTot(elem)
          
          call M33INV ( C , CINV, OKFLAG , DETC)
          call tranvoigsym( CINV , voigCINV )
  
          ! --- cosserat tensor ---
          call cosseratIncomp ( E ,VolsMats(elem), pe , S )
          call ElasticityTensorsIncomp(E,VolsMats(elem),pe)
          call  BgrandeMats ( BgrandeMat , funder , F )
          ConstMatPerturbed = kappa * DETC * &
            matmul( reshape(voigCINV,[6,1]) , transpose( reshape(voigCINV,[6,1]) ) ) 
            
          ConsMat = ConsMat + ConstMatPerturbed

        end if
        ! --------------------------------------------------------------


        ! local tangent stifness matrices
        Kml        = matmul ( matmul ( transpose( BgrandeMat ) , ConsMat ) , BgrandeMat ) * jacobiandet * wi
        matauxgeom = matmul ( matmul ( transpose( funder     ) , S       ) , funder     ) * jacobiandet * wi

        ! --- construyo matriz geometrica local ---
        do i=1,4
          do j=1,4
            Kgl( (i-1)*3+1 , (j-1)*3+1 ) = matauxgeom(i,j)
            Kgl( (i-1)*3+2 , (j-1)*3+2 ) = matauxgeom(i,j)
            Kgl( (i-1)*3+3 , (j-1)*3+3 ) = matauxgeom(i,j)
          end do
        end do

        ! Recorro todos los grados de libertad de cada elemento y veo si son de neumann 
        localiK: do i = 1,12
          localjK: do j = i, 12
            !
            aux1 = min( dofselem(i)  , dofselem(j) )
            aux2 = max( dofselem(i)  , dofselem(j) )


            ! veo si los dos grados de libertad ESTAN en los de neumann
            ! si en los vectores de los neuman ordenados los dos grados de lib NO dan cero
            ! entonces esa entrada de la matriz hay que ponerla en la matriz global hessiana neuman neuman 
            if ( ( NeumSortDOFs(int(aux1)) .ne. 0 ) .and. (NeumSortDOFs(int(aux2)) .ne. 0 ) ) then
              contador                        = contador + 1
              HessMatIndexes ( contador , : ) = [ NeumSortDOFs(int(aux1)) , NeumSortDOFs(int(aux2)) ]
              HessMatSparse  ( contador     ) = Kgl(i,j) + Kml(i,j) 
              
!~             if ( ( aux1 == 4) .and. ( aux2==4) ) then
!~               write(*,*) "WTFFFFFFFFFFF"
!~               write( *,*) HessMatIndexes(contador,1)
!~               stop
!~             end if



            ! veo si el grado de libertad mas grande esta en neumann
            elseif ( ( NnonzeroDiriDOFs > 0 ) .and. (it == 1) ) then
              if ( any( NonzeroDiriDOFs .eq. int(aux2) ) ) then
                contadorHNeNzDi                        = contadorHNeNzDi + 1
                HessMatIndexesNeumNonzeroDiri ( contadorHNeNzDi , : ) = &
                  [ int(aux1) , int(aux2) ]
                HessMatSparseNeumNonzeroDiri  ( contadorHNeNzDi     ) = &
                  Kgl(i,j) + Kml(i,j) 

              end if
            end if
            !
          end do localjK
        end do localiK
        ! --------------------------------

        call tranvoigsym(S,voigS)
        if (compresibleflag == 1 ) then

          GradMat ( dofselem ) = GradMat ( dofselem ) + matmul( transpose ( BgrandeMat ) , voigS ) * jacobiandet * wi

        elseif (compresibleflag == 2 ) then

          GradMat ( dofselem ) = GradMat ( dofselem ) + matmul( transpose ( BgrandeMat ) , voigS + &
          ( kappa * ( M33DET(F) -1 )  - pe ) * M33DET(F) * voigCINV ) * jacobiandet * wi
        end if


      end do elements_do
      

!~ write(*,*) "min", minval( HessMatIndexes(1:contador,1))
!~               stop
      
      !
      if (contador >= contadormaxHess ) then
        write(*,*) "el contador de la hessmat supero el numero alocado"
        stop
      end if
      ! =============================================================

    
      DisplacNeum = 0.0d+0
      DisplacTot  = 0.0d+0
  
      GradMat = GradMat - Fext

      if ( it == 1 ) FextIt1 = Fext  

      if ( (it == 1) .and. ( NnonzeroDiriDOFs > 0 ) ) then

        resul2 = 0.0d+0
        if (indTime == 1 ) then
          call sparmatvecmul ( HessMatIndexesNeumNonzeroDiri( 1:contadorHNeNzDi , : ) &
            , HessMatSparseNeumNonzeroDiri(1:contadorHNeNzDi) , DispBC , resul2 )
        else
          call sparmatvecmul ( HessMatIndexesNeumNonzeroDiri( 1:contadorHNeNzDi , : ) &
            , HessMatSparseNeumNonzeroDiri(1:contadorHNeNzDi) , DispBC , resul2 )
        end if
        !
        GradMat = GradMat + resul2
        FextIt1 = FextIt1 + resul2

      end if

      DisplacNeum = - GradMat ( NeumDOFs(1:NNeumDOFs) )

      ! ===============================================================
      ! =======   Linear system resolution  ===========================
  
  
      if ( linsolverflag == 1 ) then ! dgesv lapack solver

        DisplacTot  = 0.0d+0

        pivot = 0.0
        Kmg = 0.0

        do i=1,contador
          if ( HessMatIndexes(i,1) == HessMatIndexes(i,2) ) then
            Kmg ( NeumDOFs(HessMatIndexes(i,1)) , NeumDOFs(HessMatIndexes(i,2)) ) = &
            Kmg ( NeumDOFs(HessMatIndexes(i,1)) , NeumDOFs(HessMatIndexes(i,2)) ) + HessMatSparse(i)
         else
            Kmg ( NeumDOFs(HessMatIndexes(i,1)) , NeumDOFs(HessMatIndexes(i,2)) ) = &
            Kmg ( NeumDOFs(HessMatIndexes(i,1)) , NeumDOFs(HessMatIndexes(i,2)) ) + HessMatSparse(i)
            Kmg ( NeumDOFs(HessMatIndexes(i,2)) , NeumDOFs(HessMatIndexes(i,1)) ) = &
            Kmg ( NeumDOFs(HessMatIndexes(i,1)) , NeumDOFs(HessMatIndexes(i,2)) )
         end if
        end do


        call DGESV( NNeumDOFs , 1 , Kmg ( NeumDOFs(1:NNeumDOFs) , NeumDOFs(1:NNeumDOFs) ) , NNeumDOFs ,&
                       pivot , DisplacNeum , NNeumDOFs, ok)



      elseif (linsolverflag == 2 ) then  ! hsl solver

        call cpu_time(timesolver1)

!~         n  = NNeumDOFs
!~         ne = contador

!~         matrix%n = n
!~         matrix%ne = ne
        
!~         ! Allocate arrays of appropriate sizes
!~         allocate( matrix%val(ne), matrix%row(ne), matrix%col(ne))
!~         allocate( sysb(n) , sysx(n) )
        
!~         ! Read matrix and right-hand side
!~         matrix%row = HessMatIndexes(1:ne,1)
!~         matrix%col = HessMatIndexes(1:ne,2)
!~         matrix%val = HessMatSparse(1:ne)


!~         if (indtime==1) then   
!~           open (unit = 40 , file = "t1.dat" )
!~           do i=1,ne
!~             write(40,'(i10,A,i10,A,e20.10)') HessMatIndexes(i,1), " ", HessMatIndexes(i,2), " ", HessMatSparse(i)
!~           end do
!~           close(40)
!~         elseif (indtime==2) then
!~           open (unit = 40 , file = "t2.dat" )
!~           do i=1,ne
!~             write(40,'(i10,A,i10,A,e20.10)') HessMatIndexes(i,1), " ", HessMatIndexes(i,2), " ", HessMatSparse(i)
!~           end do
!~           close(40)        
!~         end if

     
!~         sysb = DisplacNeum
      
!~         ! Initialize the structures
!~         call ma57_initialize(factors,control)
    
!~         !  thiiis sets warnings off . for example the duplicates warning
!~         control%wp = -1
        
!~         ! Analyse
!~         call ma57_analyse(matrix,factors,control,ainfo)
!~         !    write(*,*) "duplicates", ainfo%dup, "ainfoflag", ainfo%flag
!~         if (ainfo%flag<0) then
!~           write(6,'(a,i2)') &
!~             ' Failure of ma57_analyse with ainfo%flag=', ainfo%flag
!~           stop
!~         end if
      
!~         ! Factorize
!~         call ma57_factorize(matrix,factors,control,finfo)
!~         if (finfo%flag<0) then
!~           write(6,'(a,i2)') &
!~           ' Failure of ma57_factorize with finfo%flag=', finfo%flag
!~           stop
!~         end if
      
!~         ! Solve without refinement
!~         sysx = sysb
!~         call ma57_solve(matrix,factors,sysx,control,sinfo)
      
!~         ! Perform one refinement
!~         call ma57_solve(matrix,factors,sysx,control,sinfo,sysb)
      
!~         DisplacNeum = sysx
      
!~         ! Clean up
!~         deallocate( matrix%val, matrix%row, matrix%col, sysb, sysx )
      
!~         call ma57_finalize(factors,control,sysinfo)

        call cpu_time(timesolver2)
        
!~         write(*,'(1e12.3)')  timesolver2-timesolver1

      end if
      ! ===============================================================
  

      ! --------------------------------------------------------------------------

      DisplacTot( NeumDOFs(1:NNeumDOFs)) = DisplacNeum(1:NNeumDOFs)

      if ( (it == 1) .and. ( NnonzeroDiriDOFs > 0 ) ) then
        DisplacTot( NonzeroDiriDOFs(1:NnonzeroDiriDOFs) ) = DispBC( NonzeroDiriDOFs(1:NnonzeroDiriDOFs) )
      end if


! aca tendria q ir la presion

      do i =1,nnod
        MatDisplacTot(i,:) = DisplacTot( ((i-1)*3+1):(i*3) )
      end do


      if (compresibleflag == 1 ) then

! ========================== POOOOOOOOOOSSSSSSSS   PROCESSSOOOOOOO  ==================================

        elements_p_do: do elem = 1, ntet
        
          nodeselem =  Cellsconec( elem ,1:4)
  
          call nods2dofs( nodeselem , 4 , 3 , dofselem )
  
          DispTotElem = DisplacTot ( dofselem )
  !~         write(*,'(12e12.4)') DispTotElem
        
          Kml = 0.0d+0
          Kgl = 0.0d+0
          
          elecoordmat(1,1:4) = nodes( nodeselem , 1 )
          elecoordmat(2,1:4) = nodes( nodeselem , 2 )
          elecoordmat(3,1:4) = nodes( nodeselem , 3 )
    
          elecoordspa(1,1:4) = ldnodesdef( nodeselem , 1 )
          elecoordspa(2,1:4) = ldnodesdef( nodeselem , 2 )
          elecoordspa(3,1:4) = ldnodesdef( nodeselem , 3 )
        
          xi = 0.25d+0
          wi = 1.0d+0/6.0d+0
        
          call ShapeFun( 1 , 0 , xi, xi , xi , fun )
          call ShapeFun( 1 , 1 , xi, xi , xi , deriv )
  
          ! jacobiana
          jacobianmat    = matmul ( elecoordmat , transpose( deriv) )
          jacobianmatspa = matmul ( elecoordspa , transpose( deriv) )
          
          ! material coord
          call M33INV ( jacobianmat    , invjacobianmat    , okflag , jacobiandet    )
          ! spatial coord
          call M33INV ( jacobianmatspa , invjacobianmatspa , okflag , jacobiandetspa )
          
          ! --------------------------------------
          if ( jacobiandet <= 0.0 ) then
            write(*,*) "Error volumen ",i, "negativoo!! ", jacobiandet
            write(*,*) "hay como: ", Ntet, "elementos en total"
            stop
          end if
    
          funder = matmul ( transpose ( invjacobianmat ) , deriv )
          
          ! --- disp grad ---
          H = matmul( (elecoordspa - elecoordmat )  , transpose ( funder ) )
    
          ! --- deformation grad ---
          F = ojo + H
  
          ! --- Lagrange tensor ---
          E = 0.5D+0 * ( H + transpose(H) + matmul( transpose(H) , H ) )
  
          C = 2.0D+0 * E + ojo
          call  BgrandeMats ( BgrandeMat , funder , F )
  
          call M33INV ( C , CINV, OKFLAG , DETC)
          call tranvoigsym( CINV , voigCINV )
  
        end do elements_p_do

      
      elseif (compresibleflag == 2 ) then

        elements_pincomp_do: do elem = 1, ntet
        
          nodeselem =  Cellsconec( elem ,1:4)
  
          call nods2dofs( nodeselem , 4 , 3 , dofselem )
  
          DispTotElem = DisplacTot ( dofselem )
  !~         write(*,'(12e12.4)') DispTotElem
        
          Kml = 0.0d+0
          Kgl = 0.0d+0
          
          elecoordmat(1,1:4) = nodes( nodeselem , 1 )
          elecoordmat(2,1:4) = nodes( nodeselem , 2 )
          elecoordmat(3,1:4) = nodes( nodeselem , 3 )
    
          elecoordspa(1,1:4) = ldnodesdef( nodeselem , 1 )
          elecoordspa(2,1:4) = ldnodesdef( nodeselem , 2 )
          elecoordspa(3,1:4) = ldnodesdef( nodeselem , 3 )
        
          xi = 0.25d+0
          wi = 1.0d+0/6.0d+0
        
          call ShapeFun( 1 , 0 , xi, xi , xi , fun )
          call ShapeFun( 1 , 1 , xi, xi , xi , deriv )
  
          ! jacobiana
          jacobianmat    = matmul ( elecoordmat , transpose( deriv) )
          jacobianmatspa = matmul ( elecoordspa , transpose( deriv) )
          
          ! material coord
          call M33INV ( jacobianmat    , invjacobianmat    , okflag , jacobiandet    )
          ! spatial coord
          call M33INV ( jacobianmatspa , invjacobianmatspa , okflag , jacobiandetspa )
          
          ! --------------------------------------
          if ( jacobiandet <= 0.0 ) then
            write(*,*) "Error volumen ",i, "negativoo!! ", jacobiandet
            write(*,*) "hay como: ", Ntet, "elementos en total"
            stop
          end if
    
          funder = matmul ( transpose ( invjacobianmat ) , deriv )
          
          ! --- disp grad ---
          H = matmul( (elecoordspa - elecoordmat )  , transpose ( funder ) )
    
          ! --- deformation grad ---
          F = ojo + H
  
          ! --- Lagrange tensor ---
          E = 0.5D+0 * ( H + transpose(H) + matmul( transpose(H) , H ) )
  
          C = 2.0D+0 * E + ojo
          call  BgrandeMats ( BgrandeMat , funder , F )
  
          call M33INV ( C , CINV, OKFLAG , DETC)
          call tranvoigsym( CINV , voigCINV )
  
  
          if ( compresibleflag == 2 ) then
      !~    write(*,*) "detf", M33DET(F)
            pTot(elem) = kappa * ( M33DET(F) * &
            dot_product( &
              reshape( matmul( transpose(reshape(voigCINV,[6,1])) , BgrandeMat ) , [12] ) &
              , DispTotElem ) ) &
              + &
              kappa * ( M33DET(F) - 1 )
      !~         DisplacTot( ((i-1)*3+1):(i*3) )
             end if
        end do elements_pincomp_do
        
        write(*,'(A,e12.3)') "ptot(20)", pTot(10)

      end if

        
!~       SumaMatDisplacTot = SumaMatDisplacTot + MatDisplacTot
  
      ldnodesdef      = ldnodesdef + MatDisplacTot
          
      normgrad         = sqrt( dot_product ( GradMat(NeumDOFs(1:NNeumDOFs)) , GradMat(NeumDOFs(1:NNeumDOFs)) ))
      normrighthand    = sqrt( dot_product ( FextIt1(NeumDOFs(1:NNeumDOFs)) , FextIt1(NeumDOFs(1:NNeumDOFs)) ))
      normgradhist(it) = normgrad
  
      resul = 0.0d0
      call sparmatvecmul ( HessMatIndexes( 1:contador , : ) &
         , HessMatSparse(1:contador) , DisplacTot( NeumDOFs(1:NNeumDOFs) ) , resul )
  
      normresidload = &
        sqrt( dot_product ( resul + GradMat(NeumDOFs(1:NNeumDOFs)) , resul + GradMat(NeumDOFs(1:NNeumDOFs)) ) )

      write ( * , '(i3,A,e20.12,A,e20.12,A)' ) it ," | ", normgradhist(it) , " | " , normresidload ," | "
  
    end do newtonraphson
    ! --------------------------------

    call cpu_time(time2)
  
    write(*,*) "linear system solver: ", linsolverflag , "   time:", time2-time1, "seconds"
    write ( * , '(A)' ) "---   Loop end   --------"

    ! ----------------------------------------------------------------------------------------------
    ! Obtiene desplazamientos de nodos de control interpolando entre desplazamientos nodales del 
    ! elemento de control.
    do cn=1, ncontrolnodes
      aux3 = ControlElem(cn)

      write(14,'(e10.3,i8,e28.15,e28.15,e28.15)') sum( times(1:indTime)) , cn &
      , dot_product(ldnodesdef(Cellsconec(aux3,:),1), ControlElemPesos(cn,:) ) &
      , dot_product(ldnodesdef(Cellsconec(aux3,:),2), ControlElemPesos(cn,:) ) &
      , dot_product(ldnodesdef(Cellsconec(aux3,:),3), ControlElemPesos(cn,:) )

    end do
    ! ---------------------------------------------------------------------------------------------
    
    
    if (NNonZeroDiriDOFs>0) then
      write(15,'(A,e24.15,A,e24.15,A,i4)') "  ",times(indTime) ," " &
      ,  sqrt( dot_product( GradMat(NonzeroDiriDOFs),GradMat(NonzeroDiriDOFs)  ) ),"  ",it
    else
      write(15,'(A,e24.15,A,e24.15,A,i4)') "  ",times(indTime) ," " &
      ,  sqrt( dot_product( GradMat(NeumDOFs),GradMat(NeumDOFs)  ) ),"  ",it
    end if

    call chdir("Deformed/Step_by_step")
    
    if ( Ntimes <= 9999 ) then
      if     ( indtime < 10 ) then
        write( auxchar,'(A,A,i1)') "Def_", "000", indtime
      elseif   (  ( indtime >= 10 ) .and. ( indtime < 100 ) ) then
        write( auxchar,'(A,A,i2)') "Def_", "00", indtime
      elseif   ( ( indtime >= 100 ) .and. (indtime < 1000 ) ) then
        write( auxchar,'(A,A,i3)') "Def_", "0", indtime
      else
        write( auxchar,'(A,i4)') "Def_", indtime
      end if
    end if

    ! ---------------------------------------------
    ! Generacion archivo vtk
    if (compresibleflag == 1 ) then
      call vtkgen_point &
        (  auxchar ,"displacs" &
        , NNodPE , 3   &
        , ldnodesdef - nodes &
        , Cellsconec(1:ntet,1:4) , ldnodesdef &
        )


    else if (compresibleflag == 2 ) then
      call vtkgen &
        (  auxchar ,"displacs", "presiones" &
        , NNodPE , 3   &
        , ldnodesdef - nodes &
        , reshape( pTot ,[ntet,1] ) &
        , Cellsconec(1:ntet,1:4) , ldnodesdef &
        )
    end if
    
    call chdir("../../")
    ! ----------------------------------------------------
    
!~ ####################################################################    
!~ ####################################################################     
!~     Comienzo edición de Pablo - Cálculo de tensores para el caso compresible 03-11-2014
!~ ####################################################################
!~ ####################################################################
      
      elements_do2: do elem = 1, ntet
      
        nodeselem =  Cellsconec( elem ,1:4)

        call nods2dofs( nodeselem , 4 , 3 , dofselem )
      
        Kml = 0.0d+0
        Kgl = 0.0d+0
        
        elecoordmat(1,1:4) = nodes( nodeselem , 1 )
        elecoordmat(2,1:4) = nodes( nodeselem , 2 )
        elecoordmat(3,1:4) = nodes( nodeselem , 3 )
  
        elecoordspa(1,1:4) = ldnodesdef( nodeselem , 1 )
        elecoordspa(2,1:4) = ldnodesdef( nodeselem , 2 )
        elecoordspa(3,1:4) = ldnodesdef( nodeselem , 3 )
      
        xi = 0.25d+0
        wi = 1.0d+0/6.0d+0
      
        call ShapeFun( 1 , 0 , xi, xi , xi , fun )
        call ShapeFun( 1 , 1 , xi, xi , xi , deriv )

        ! jacobiana
        jacobianmat    = matmul ( elecoordmat , transpose( deriv) )
        jacobianmatspa = matmul ( elecoordspa , transpose( deriv) )
        
        ! material coord
        call M33INV ( jacobianmat    , invjacobianmat    , okflag , jacobiandet    )
        ! spatial coord
        call M33INV ( jacobianmatspa , invjacobianmatspa , okflag , jacobiandetspa )
  
        funder = matmul ( transpose ( invjacobianmat ) , deriv )
        
        ! --- disp grad ---
        H = matmul( (elecoordspa - elecoordmat )  , transpose ( funder ) )
  
        ! --- deformation grad ---
        F = ojo + H

        ! --- Lagrange tensor ---
        E = 0.5D+0 * ( H + transpose(H) + matmul( transpose(H) , H ) )

        C = 2.0D+0 * E + ojo

      if ( compresibleflag == 1 ) then

        ! --- cosserat tensor ---
        call cosserat ( E, VolsMats(elem) , S )
        AllCosMat(1,elem) = S(1,1)
        AllCosMat(2,elem) = S(2,2)
        AllCosMat(3,elem) = S(3,3)
        AllCosMat(4,elem) = S(1,2)
        AllCosMat(5,elem) = S(2,3)
        AllCosMat(6,elem) = S(1,3)
        
        ! --- piola tensor ---
        P = matmul( F , S )
        AllPioMat(1,elem) = P(1,1)
        AllPioMat(2,elem) = P(2,2)
        AllPioMat(3,elem) = P(3,3)
        AllPioMat(4,elem) = P(1,2)
        AllPioMat(5,elem) = P(2,3)
        AllPioMat(6,elem) = P(1,3)
        AllPioMat(7,elem) = P(2,1)
        AllPioMat(8,elem) = P(3,2)
        AllPioMat(9,elem) = P(3,1)
        
        ! --- cauchy tensor ---
        call M33INV ( F , invF , okflag , detF )
        Sig = 1.0D+0 / detF *  matmul( F , P )
        AllSigMat(1,elem) = Sig(1,1)
        AllSigMat(2,elem) = Sig(2,2)
        AllSigMat(3,elem) = Sig(3,3)
        AllSigMat(4,elem) = Sig(1,2)
        AllSigMat(5,elem) = Sig(2,3)
        AllSigMat(6,elem) = Sig(1,3)
        
        ! --- von misses - cauchy tensor---
        
        sigVM_Cau (elem,1) = sqrt( (Sig(1,1) +Sig(2,2) + Sig(3,3) )**2 &
        - 3*(Sig(1,1)*Sig(2,2) + Sig(1,1)*Sig(3,3)+ Sig(2,2)*Sig(3,3)   &
        - Sig(1,2)**2 - Sig(1,3)**2 - Sig(2,3)**2 ) )
        
        

      elseif (compresibleflag == 2 ) then

      end if

    end do elements_do2
    
    ! --- vtk generation ---
    
    call chdir("Stress/Cosserat")
    
    if ( Ntimes <= 9999 ) then
      if     ( indtime < 10 ) then
        write( auxchar,'(A,A,i1)') "Cos_", "000", indtime
      elseif   (  ( indtime >= 10 ) .and. ( indtime < 100 ) ) then
        write( auxchar,'(A,A,i2)') "Cos_", "00", indtime
      elseif   ( ( indtime >= 100 ) .and. ( indtime < 1000 ) ) then
        write( auxchar,'(A,A,i3)') "Cos_", "0", indtime
      else
        write( auxchar,'(A,i4)') "Cos_", indtime
      end if
    end if
    
    call vtkgen_cell &
    (  auxchar , "Cosserat" &
    , NNodPE , 3   &
    , transpose(AllCosMat) &
    , Cellsconec(1:ntet,1:4) , ldnodesdef &
    )
    
    call chdir("../")
    call chdir("Piola")
    
    if ( Ntimes <= 9999 ) then
      if     ( indtime < 10 ) then
        write( auxchar,'(A,A,i1)') "Pio_", "000", indtime
      elseif   (  ( indtime >= 10 ) .and. ( indtime < 100 ) ) then
        write( auxchar,'(A,A,i2)') "Pio_", "00", indtime
      elseif   ( ( indtime >= 100 ) .and. ( indtime < 1000 ) ) then
        write( auxchar,'(A,A,i3)') "Pio_", "0", indtime
      else
        write( auxchar,'(A,i4)') "Pio_", indtime
      end if
    end if
    
    
    call vtkgen_cell &
    (  auxchar , "Piolaaaa" &
    , NNodPE , 3   &
    , transpose(AllPioMat) &
    , Cellsconec(1:ntet,1:4) , ldnodesdef &
    )
    
    call chdir("../")
    call chdir("Cauchy")
    
    if ( Ntimes <= 9999 ) then
      if     ( indtime < 10 ) then
        write( auxchar,'(A,A,i1)') "Cau_", "000", indtime
      elseif   (  ( indtime >= 10 ) .and. ( indtime < 100 ) ) then
        write( auxchar,'(A,A,i2)') "Cau_", "00", indtime
      elseif   ( ( indtime >= 100 ) .and. ( indtime < 1000 ) ) then
        write( auxchar,'(A,A,i3)') "Cau_", "0", indtime
      else
        write( auxchar,'(A,i4)') "Cau_", indtime
      end if
    end if
    
    call vtkgen_cell &
    (  auxchar , "Cauchyyy" &
    , NNodPE , 3   &
    , transpose(AllSigMat) &
    , Cellsconec(1:ntet,1:4) , ldnodesdef &
    )
    
    call chdir("../")
    call chdir("Von_Misses")
    
    if ( Ntimes <= 9999 ) then
      if     ( indtime < 10 ) then
        write( auxchar,'(A,A,i1)') "VM_", "000", indtime
      elseif   (  ( indtime >= 10 ) .and. ( indtime < 100 ) ) then
        write( auxchar,'(A,A,i2)') "VM_", "00", indtime
      elseif   ( ( indtime >= 100 ) .and. ( indtime < 1000 ) ) then
        write( auxchar,'(A,A,i3)') "VM_", "0", indtime
      else
        write( auxchar,'(A,i4)') "VM_", indtime
      end if
    end if
    
    call vtkgen_cell &
    (  auxchar , "Von_Miss" &
    , NNodPE , 3   &
    , sigVM_Cau &
    , Cellsconec(1:ntet,1:4) , ldnodesdef &
    )
    
    
    call chdir("../../")
  
  do cn=1, ncontrolnodes
    aux3 = ControlElem(cn)
    if (cn<ncontrolnodes) then
      write(16,'(A,e15.5,A,e26.16,A,e26.16,A,e26.16,A,e26.16,A,e26.16,A,e26.16,A)',advance='no') "||"&
      ,times(indTime) ,"|"&
      , AllSigMat(1,aux3),"|", AllSigMat(2,aux3),"|", AllSigMat(3,aux3),"|"&
      , AllSigMat(4,aux3),"|", AllSigMat(5,aux3),"|", AllSigMat(6,aux3),"||   "  
    else
      write(16,'(A,e15.5,A,e26.16,A,e26.16,A,e26.16,A,e26.16,A,e26.16,A,e26.16,A)') "||"&
      ,times(indTime) ,"|"&
      , AllSigMat(1,aux3),"|", AllSigMat(2,aux3),"|", AllSigMat(3,aux3),"|"&
      , AllSigMat(4,aux3),"|", AllSigMat(5,aux3),"|", AllSigMat(6,aux3),"||" 
    end if
  end do
  
  
  do cn=1, ncontrolnodes
    aux3 = ControlElem(cn)
    if (cn<ncontrolnodes) then
      write(17,'(A,e15.5,A,e26.16,A,e26.16,A,e26.16,A,e26.16,A,e26.16,A,e26.16,A)',advance='no') "||"&
      ,times(indTime) ,"|"&
      , AllCosMat(1,aux3),"|", AllCosMat(2,aux3),"|", AllCosMat(3,aux3),"|"&
      , AllCosMat(4,aux3),"|", AllCosMat(5,aux3),"|", AllCosMat(6,aux3),"||   "   
    else
      write(17,'(A,e15.5,A,e26.16,A,e26.16,A,e26.16,A,e26.16,A,e26.16,A,e26.16,A)') "||",times(indTime) ,"|"&
      , AllCosMat(1,aux3),"|", AllCosMat(2,aux3),"|", AllCosMat(3,aux3),"|"&
      , AllCosMat(4,aux3),"|", AllCosMat(5,aux3),"|", AllCosMat(6,aux3),"||"    
    end if
  end do
  
  
  ! ----- escribo entradas de tensor de Piola -----
  do cn=1, ncontrolnodes
    aux3 = ControlElem(cn)
    if (cn<ncontrolnodes) then
      write(18,'(A,e15.5,A,e26.16,A,e26.16,A,e26.16,A,e26.16,A,e26.16,A,e26.16,A,e26.16,A,e26.16,A,e26.16,A)',advance='no')" "&
      ,times(indTime) ," "&
      , AllPioMat(1,aux3)," ", AllPioMat(2,aux3)," ", AllPioMat(3,aux3)," "&
      , AllPioMat(4,aux3)," ", AllPioMat(5,aux3)," ", AllPioMat(6,aux3)," "&
      , AllPioMat(7,aux3)," ", AllPioMat(8,aux3)," ", AllPioMat(9,aux3),"   "   
    else
      write(18,'(A,e15.5,A,e26.16,A,e26.16,A,e26.16,A,e26.16,A,e26.16,A,e26.16,A,e26.16,A,e26.16,A,e26.16,A)') " "&
      ,times(indTime) ," "&
      , AllPioMat(1,aux3)," ", AllPioMat(2,aux3)," ", AllPioMat(3,aux3)," "&
      , AllPioMat(4,aux3)," ", AllPioMat(5,aux3)," ", AllPioMat(6,aux3)," "&
      , AllPioMat(7,aux3)," ", AllPioMat(8,aux3)," ", AllPioMat(9,aux3)," "
    end if
  end do
  ! ---------------------------------------
  
  
  do cn=1, ncontrolnodes
    aux3 = ControlElem(cn)
    if (cn<ncontrolnodes) then
      write(19,'(A,e15.5,A,e26.16,A)',advance='no') "  ",times(indTime) ," "&
      , sigVM_Cau(aux3,1),"     "   
    else
      write(19,'(A,e15.5,A,e26.16,A)') "  ",times(indTime) ," "&
      , sigVM_Cau(aux3,1),"  "
    end if
  end do
  

  
!~ ####################################################################    
!~ ####################################################################     
!~     Fin edición de Pablo - Cálculo de tensores para el caso compresible 03-11-2014
!~ ####################################################################
!~ #################################################################### 

  end do times_do
  
  allocate( ControlElemCoord(3,ncontrolnodes) )
  
  do j=1,ncontrolnodes
    aux3 = ControlElem(j)
    ControlElemCoord(1,j) = sum(nodes(Cellsconec(aux3,:),1))/4.0 
    ControlElemCoord(2,j) = sum(nodes(Cellsconec(aux3,:),2))/4.0
    ControlElemCoord(3,j) = sum(nodes(Cellsconec(aux3,:),3))/4.0
  end do
  
  do j=1, ncontrolnodes
    if (j<ncontrolnodes) then
!~       write(14,'(A,A,A,A)',advance='no') "-------------------" , "-------------------------", "-------------------------"&
!~       , "-------------------------   "
    else
!~       write(14,'(A,A,A,A)') "-------------------" , "-------------------------", "-------------------------"&
!~       , "-------------------------"
    end if
  end do
  
!~   if (NNonZeroDiriDOFs>0) then
    write(15,'(A)') "-----------------------------------------------------"
!~   end if
  
  do j=1, ncontrolnodes
    if (j<ncontrolnodes) then
      write(16,'(A,A)',advance='no') "--------------------------------------------------------------------"&
,"-----------------------------------------------------------------------------------------------------------------   "
    else
      write(16,'(A,A)') "--------------------------------------------------------------------"&
,"-----------------------------------------------------------------------------------------------------------------"
    end if
  end do
  do j=1, ncontrolnodes
    if (j<ncontrolnodes) then
      write(16,'(A,e26.16,A,e26.16,A,e26.16,A,A)',advance='no') "|| Barycenter Coord:  X=", ControlElemCoord(1,j),",  Y="&
      , ControlElemCoord(2,j),",  Z=", ControlElemCoord(3,j),"                                                         "&
      ,"          ||   "  
    else
      write(16,'(A,e26.16,A,e26.16,A,e26.16,A,A)') "|| Barycenter Coord:  X=", ControlElemCoord(1,j),",  Y="&
      , ControlElemCoord(2,j),",  Z=", ControlElemCoord(3,j),"                                                           "&
      ,"        ||" 
    end if
  end do
  do j=1, ncontrolnodes
    if (j<ncontrolnodes) then
      write(16,'(A,A)',advance='no') "--------------------------------------------------------------------"&
,"-----------------------------------------------------------------------------------------------------------------   "
    else
      write(16,'(A,A)') "--------------------------------------------------------------------"&
,"-----------------------------------------------------------------------------------------------------------------"
    end if
  end do
  
  
  
  do j=1, ncontrolnodes
    if (j<ncontrolnodes) then
      write(17,'(A,A)',advance='no') "--------------------------------------------------------------------"&
,"-----------------------------------------------------------------------------------------------------------------   "
    else
      write(17,'(A,A)') "--------------------------------------------------------------------"&
,"-----------------------------------------------------------------------------------------------------------------"
    end if
  end do
  do j=1, ncontrolnodes
    if (j<ncontrolnodes) then
      write(17,'(A,e26.16,A,e26.16,A,e26.16,A,A)',advance='no') "|| Barycenter Coord:  X=", ControlElemCoord(1,j),",  Y="&
      , ControlElemCoord(2,j),",  Z=", ControlElemCoord(3,j),"                                                         "&
      ,"          ||   "  
    else
      write(17,'(A,e26.16,A,e26.16,A,e26.16,A,A)') "|| Barycenter Coord:  X=", ControlElemCoord(1,j),",  Y="&
      , ControlElemCoord(2,j),",  Z=", ControlElemCoord(3,j),"                                                           "&
      ,"        ||" 
    end if
  end do
  do j=1, ncontrolnodes
    if (j<ncontrolnodes) then
      write(17,'(A,A)',advance='no') "--------------------------------------------------------------------"&
,"-----------------------------------------------------------------------------------------------------------------   "
    else
      write(17,'(A,A)') "--------------------------------------------------------------------"&
,"-----------------------------------------------------------------------------------------------------------------"
    end if
  end do
  
  
  do j=1, ncontrolnodes
    if (j<ncontrolnodes) then
      write(18,'(A,A,A)',advance='no') "--------------------------------------------------------------------------------------"&
       ,"----------------------------------------------------------------------------------------------------------------------"&
       ,"----------------------------------------------------------   "
    else
      write(18,'(A,A,A)') "--------------------------------------------------------------------------------------"&
       ,"----------------------------------------------------------------------------------------------------------------------"&
       ,"----------------------------------------------------------"
    end if
  end do
  do j=1, ncontrolnodes
    if (j<ncontrolnodes) then
      write(18,'(A,e26.16,A,e26.16,A,e26.16,A,A)',advance='no') "|| Barycenter Coord:  X=", ControlElemCoord(1,j),",  Y="&
      , ControlElemCoord(2,j),",  Z=", ControlElemCoord(3,j),"                                                          "&
      ,"                                                                                          ||   "  
    else
      write(18,'(A,e26.16,A,e26.16,A,e26.16,A,A)') "|| Barycenter Coord:  X=", ControlElemCoord(1,j),",  Y="&
      , ControlElemCoord(2,j),",  Z=", ControlElemCoord(3,j),"                                            "&
      ,"                                                                                                        ||" 
    end if
  end do
  do j=1, ncontrolnodes
    if (j<ncontrolnodes) then
      write(18,'(A,A,A)',advance='no') "--------------------------------------------------------------------------------------"&
       ,"----------------------------------------------------------------------------------------------------------------------"&
       ,"----------------------------------------------------------   "
    else
      write(18,'(A,A,A)') "--------------------------------------------------------------------------------------"&
       ,"----------------------------------------------------------------------------------------------------------------------"&
       ,"----------------------------------------------------------"
    end if
  end do
  
  
  do j=1, ncontrolnodes
    if (j<ncontrolnodes) then
      write(19,'(A,A)',advance='no') "------------------" , "----------------------------   "
    else
      write(19,'(A,A)') "------------------" , "----------------------------"
    end if
  end do
  do j=1, ncontrolnodes
    if (j<ncontrolnodes) then
      write(19,'(A)',advance='no') "|| Barycenter Coord:                        ||   "
    else
      write(19,'(A)') "|| Barycenter Coord:                        ||"
    end if
  end do
  do j=1, ncontrolnodes
    if (j<ncontrolnodes) then
      write(19,'(A,e26.16,A)',advance='no') "||         X=",ControlElemCoord(1,j),"     ||   "
    else
      write(19,'(A,e26.16,A)') "||         X=",ControlElemCoord(1,j),"     ||"
    end if
  end do
  do j=1, ncontrolnodes
    if (j<ncontrolnodes) then
      write(19,'(A,e26.16,A)',advance='no') "||         Y=",ControlElemCoord(2,j),"     ||   "
    else
      write(19,'(A,e26.16,A)') "||         Y=",ControlElemCoord(2,j),"     ||"
    end if
  end do
  do j=1, ncontrolnodes
    if (j<ncontrolnodes) then
      write(19,'(A,e26.16,A)',advance='no') "||         Z=",ControlElemCoord(3,j),"     ||   "
    else
      write(19,'(A,e26.16,A)') "||         Z=",ControlElemCoord(3,j),"     ||"
    end if
  end do
  do j=1, ncontrolnodes
    if (j<ncontrolnodes) then
      write(19,'(A,A)',advance='no') "------------------" , "----------------------------   "
    else
      write(19,'(A,A)') "------------------" , "----------------------------"
    end if
  end do
  
  close(14)
!~   if (NNonZeroDiriDOFs>0) then
    close(15)
!~   end if
  close(16)
  close(17)
  close(18)
  close(19)
  ! ======================================================================
  ! ======================================================================



  ! ======================================================================
  ! =========   4. Pos Process    ========================================
  ! ======================================================================


  ! --- tensions! ---

!~   do i=1,ntet
!~     !
!~     elemind = 0
!~     !
!~     call nods2dofs ( Cellsconec(i,1:4) , 4 , 3 , elemind )
!~ 
!~     Be = BsMat(1:6 , ((i-1)*12+1):(i*12) )
!~ 
!~     Ue = DisplacTot(elemind)
!~ 
!~     sigmas(1:6,i) = matmul ( matmul ( DeeMat  , Be ) , Ue )
!~ 
!~     sigVM ( i) = sqrt( (sigmas(1,i) +sigmas(2,i) + sigmas(3,i) )**2 &
!~       - 3*(sigmas(1,i)*sigmas(2,i) + sigmas(3,i)*sigmas(2,i)+ sigmas(1,i)*sigmas(3,i)   &
!~             - sigmas(4,i)**2 - sigmas(5,i)**2 - sigmas(6,i)**2 ) )
!~ 
!~ !    sigVM ( i) = sigmas(1,i)
!~ 
!~   end do
  !
  
  
  ! =============================================================================
  ! =============================================================================
!~   allocate( UndeformedMat(nnod) )
!~     
!~   do i=1,nnod
!~     !
!~     if ( (BCMat(i,1)==2) .or. (BCMat(i,2)==2) .or. (BCMat(i,1)==3) ) then
!~       !
!~       auxint =  BCMat(i,1)*3**2 + BCMat(i,2)*3 + BCMat(i,3)
!~       !
!~       if (auxint <= 11) then
!~         !
!~         if (BCMat(i,2)==2) then
!~           !
!~           UndeformedMat(i)= auxint + 4
!~           !
!~         elseif (BCMat(i,1)==1) then
!~           !
!~           UndeformedMat(i)= auxint + 2
!~           !
!~         else
!~           !
!~           UndeformedMat(i)= auxint + 7
!~           !
!~         end if
!~         !
!~       else
!~         !
!~         UndeformedMat(i)= auxint
!~         !
!~       end if
!~       ! 
!~     else
!~       !
!~       UndeformedMat(i)= BCMat(i,1)*2**2 + BCMat(i,2)*2 + BCMat(i,3)
!~       !
!~     end if
!~     !
!~   end do

  ! =============================================================================
  ! =============================================================================

!~   allocate( DisplacTotX(nnod), DisplacTotY(nnod), DisplacTotZ(nnod) )
!~   DisplacTotX = ldnodesdef(:,1)-nodes(:,1)
!~   DisplacTotY = ldnodesdef(:,2)-nodes(:,2)
!~   DisplacTotZ = ldnodesdef(:,3)-nodes(:,3)

!~   do l=0,1
!~ 
!~     if     (l>99) then
!~        write( auxchar,'(i3)') l
!~     elseif (l>9 ) then
!~       write( auxchar,'(A,i2)') "0",l
!~     elseif (l>-1) then
!~       write( auxchar,'(A,i1)') "00",l
!~     end if
!~ 
!~     if     (l>9) then
!~       write( auxchar,'(i2)') l
!~     elseif (l>-1) then
!~       write( auxchar,'(A,i1)') "0",l
!~     end if
!~ 
!~     do i=1,Nnodes
!~ !~
!~       nodesdef(i,1:3) = nodes(i,1:3) + DisplacTot(((i-1)*3+1):(i*3))* l *0.5
!~ !~
!~     end do
!~     
!~     if (l==0) then
!~     ! --- vtk generation ---
!~ 
!~       ( "defvtk"//auxchar &



!~   call vtkgen &
!~     ( "undeform","nodsbcva", "cellsind" &
!~     , NNodPE , 3   &
!~     , reshape( [UndeformedMat(:), BCMat(:,4) , BCMat(:,5) , BCMat(:,6) ] , [nnod,4] ) &
!~     , reshape( [ (dble(j),j=1,ntet) ] ,[ntet,1] ) &
!~     , Cellsconec(1:ntet,1:4) , nodes &
!~     )


  
  call chdir("Deformed")
  ! --------------  Output nodes def ----------------------
  open (unit = 12 , file = trim(probname)//"_Deformed_nodes_position_.txt" )
  !
  write (12 , '(A)')       "---------------------------------------------------------------------------------------"
  write (12 , '(A)')       "Deformed nodes - last iteration for last step"
  write (12 , '(A)')       "---------------------------------------------------------------------------------------"
  write (12 , '(A,A,A,A)') "|      node|", "                       X|", "                       Y|", "                       Z|"
  write (12 , '(A)')       "---------------------------------------------------------------------------------------"
  do i =1,nnod
    !
    write (12 , '(A,i10,A,e24.12,A,e24.12,A,e24.12,A)') "|",i,"|", ldnodesdef(i,1),"|", ldnodesdef(i,2),"|", ldnodesdef(i,3),"|"
  end do
  !
  close(12)

  open (unit = 12 , file = trim(probname)//"_Nodes_displacement_.txt" )
  !
  write (12 , '(A)')       "---------------------------------------------------------------------------------------"
  write (12 , '(A)')       "Nodes displacement - last iteration for last step"
  write (12 , '(A)')       "---------------------------------------------------------------------------------------"
  write (12 , '(A,A,A,A)') "|      node|", "                   Xdesp|", "                   Ydesp|", "                   Zdesp|"
  write (12 , '(A)')       "---------------------------------------------------------------------------------------"
  do i =1,nnod
    !
    write (12 , '(A,i10,A,e24.12,A,e24.12,A,e24.12,A)') "|",i,"|", ldnodesdef(i,1) - nodes(i,1),"|", ldnodesdef(i,2) - nodes(i,2) &
                                                      ,"|", ldnodesdef(i,3) - nodes(i,3),"|"
  end do
  !
  close(12)
  ! -------------------------------------------------------
  call chdir("../")

  
  write(*,*) "All files written, ldfem finished. Go and open the vtks with Paraview!"

end program ldfem

! ==================================================================================================
! ==================================================================================================

subroutine nods2dofs( nods , n , degree , dofs )
  implicit none
  integer, intent(in) :: n, degree
  integer :: i, j
  integer, intent(in) :: nods(n)
  integer             :: dofs(n*degree)

  do i=1,n
    dofs( ((i-1)*degree+1):(i*degree) ) = [ ( (nods(i)-1)*degree+j,j=1,degree ) ]
  end do

end subroutine nods2dofs

subroutine tensinneprod( A , B , prod )
  implicit none
  real*8, intent(in) :: A(3,3), B(3,3)
  integer :: i, j
  real*8     :: prod

  prod = 0.0d+0

  do i=1,3
    do j=1,3
      prod = prod + A(i,j)*B(i,j)
    end do
  end do

end subroutine tensinneprod


subroutine unique(v,u)
  !
  implicit none
  !
  real*8, intent(in)    :: v(:)
  real*8    :: u(:)
  !
end subroutine unique
