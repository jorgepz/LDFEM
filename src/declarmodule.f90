
module declarmodule
  
  integer, parameter                   :: NNodPE  =   4  !& Number of nodes per element & 1x1

  ! Input params
  integer, save                        :: nnod           !& Number of nodes & 1x1
  integer, save                        :: ntet           !& Number of volume elements (assumed tetrahedrons) & 1x1
  integer, save                        :: ntri           !& Number of surface elements (assumed triangles) & 1x1
  integer, save                        :: nscalefactors  !& Scale factors quantity & 1x1
  integer, save                        :: ncontrolnodes  !& Scale factors quantity & 1x1
  
  ! Geometry 
  integer  ,save        , allocatable  :: Cellsconec(:,:)    !& Conectivity matrix of volume elements & ntet x 4 
  integer  ,save        , allocatable  :: Surfconec(:,:)     !& Conectivity matrix of surface elements & ntri x 3 
  real*8   ,save        , allocatable  :: nodes(:,:)         !& Nodes coord. matrix & nnod x 3
  real*8   ,save        , allocatable  :: ldnodesdef(:,:)    !& Deformed nodes coord. matrix & nnod x 3
  real*8   ,save        , allocatable  :: MatDisplacTot(:,:) !& Matrix with displacements in each NR step & nnod x 3

  ! Boundary Conditions 
  integer,save                         :: nBC            !& number of boundary conditions defined
  real*8,save           , allocatable  :: BCDefMat(:,:)  !& matrix with boundary conditions definitions
  integer,save          , allocatable  :: NodBCMat(:)    !& vector with nodes boundary conditions labels & nnod x 1
  integer,save          , allocatable  :: CellsBCMat(:)  !& vector with cells boundary conditions labels & nnod x 1
  integer,save          , allocatable  :: SurfBCMat(:)   !& vector with surfaces boundary conditions labels & nnod x 1

  real*8, save          , allocatable  :: Fext(:)        !& nodal loads imposed at current time and it & 3*nnod x 1
  real*8, save          , allocatable  :: Udiri(:)       !& nodal disps imposed at current time  & 3*nnod x 1
  real*8, save          , allocatable  :: Udirikm1(:)       !& nodal disps imposed at current time  & 3*nnod x 1
  real*8, save          , allocatable  :: FextIt1(:)     !& loads imposed at iteration 1               & 3*nnod x 1
  real*8, save          , allocatable  :: DisplacNeum(:) !&  
  real*8, save          , allocatable  :: DisplacTot(:)  !& 
  real*8, save          , allocatable  :: pTot(:)  !& 
  
  real*8, save          , allocatable  :: BCMat(:,:)    !&    & nnod x 6
  real*8, save          , allocatable  :: NonzeroDiriVALs(:)
  real*8, save          , allocatable  :: DispBC(:)

  integer               , allocatable  :: NeumDOFs(:) 
  integer               , allocatable  :: DiriDOFs (:)
  integer               , allocatable  :: AuxNeumDOFs(:)
  integer               , allocatable  :: NonzeroDiriDOFs(:)
  integer               , allocatable  :: NeumSortDOFs(:)  !& Array with DOFs for assemble of reduced Neuman-Neuman matrix & 3*nnod x 1 
  integer               , allocatable  :: DiriSortDOFs (:) !& 
  double precision,save , allocatable  :: UndeformedMat(:), DisplacTotX(:), DisplacTotY(:), DisplacTotZ(:)
  integer                              :: NNeumDOFs , NDiriDOFs, NnonzeroDiriDOFs


  double precision,save , allocatable  :: GradMat(:)
  double precision, save, allocatable  :: HessMat(:,:) , HessMatSparse(:), HessMatSparseNeumNonzeroDiri(:)
  integer         , save, allocatable  :: HessMatIndexes(:,:), HessMatIndexesNeumNonzeroDiri(:,:)
  integer         , save               :: contador !& contador Hessiana & estimado
  integer         , save               :: contadorHNeNzDi !& contador Hessiana NeumNonzeroDiri & estimado

  real*8,save           , allocatable  :: RightHand(:) , RightHandPunt(:), ScaleFactors(:)

  double precision,save , allocatable  :: ElemVols(:)
  

  real*8,save           , allocatable  :: BsMat(:,:) , BMat(:,:)


  ! Material properties
  integer,save                         :: nmats
  real*8,save                          :: young , nu , shear , lambda , ConsMat(6,6)
  real*8                               :: load
  integer,save          , allocatable  :: VolsMats(:) !& En la entrada i-esima tiene la propiedad física del volumen i-esimo

  real*8, save          , allocatable  :: MatsMat(:,:)
  real*8, save                         :: MatPars (10)

  integer , save :: solverflag, assembflag
  real, save     :: time1, time2

  integer, save :: constitutive_model

  ! gauss integration points
!~   real, save     :: xi, wi

end module declarmodule




!~   real*8   ,save        , allocatable  :: SumaMatDisplacTot(:,:) !:: 
!~   real*8   ,save        , allocatable  :: elecoordmat(:,:), elecoordspa(:,:)
!~   real*8   ,save        , allocatable  :: nodesdefvtkld(:,:) !& REVISAR Deformed nodes coord. matrix for paraview graphs & nnod x 3
