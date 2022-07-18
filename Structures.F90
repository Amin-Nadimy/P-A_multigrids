module Structures
    implicit none
    !IMPORTANT: THE ARRAY ARE ALWAYS DEFINED BY (COLUMN, ROW)
    !->This is structure is as follows:
    !The triangles are divided by their orientation, therefore we have two
    !kind of triangles -> Up
    !                  -> Down
    !Inside every triangle the nodes are stored as one dimension more of the array:
    !(Column, Row, Node)
    !Since we are in a triangle the codification of the nodes are not trivial, hence we
    !have defined it depending on the line that contains them (this is global,
    !it does not depend on the orientation of the triangle):
    !Numeration of the nodes
    !   For nodes in the sides:
    !   1º -> /
    !   2º -> _
    !   3º -> \
    !   4º -> Centered nodes (to allow working with CC and Mix, we in the code this will be denoted by
    !						 an "a" which is defined by a = size(Tri%UpUF,3) )
    !The X direction is the _ and the Y direction is the \ as a consequence our triangles numeration is the following:
    !
    !                           /\
    !                          /  \     <- Up (1,1,:)
    !                         /____\
    !                        /\    /\   <- Down (1,1,:)
    !                       /  \  /  \
    !         Up(1,2,:)->  /____\/____\ <- Up (2,2,:)
    !
    !**************Triangle coordinates*****************
    !                   Y Direction
    !                    \
    !                     \   angle = alpha
    !                      \_______ X Direction
    !
    !Numeration and nomenclature, for internal meshed triangles:
    !                Xp3          Xp1 __l2__ Xp2
    !                / \              \    /
    !           l1->/   \<-l3      l1->\  /<-l3
    !              /_____\              \/
    !           Xp1  l2   Xp2           Xp3
    !
    !The way to acces the next triangles is the following:
    !Up triangles are surrounded by down triangles and viceversa
    !Taking this into account we realise that we have to act different deppending of our triangle:
    !Example with Up triangles: Up(2,3) will have the following     -> Down (2,2) -> One column less
    !                                                               -> Down (2,3) -> Same coordinates
    !                                                               -> Down (1,2) -> One row and one column less
    !
    !Example with Down triangles: Down(2,2) will have the following -> Up (2,3) -> One column higher
    !                                                               -> Up (2,2) -> Same coordinates
    !                                                               -> Up (3,3) -> One row and one column more
    !
    !->EXPLANATION OF THE ARRAYS FOR TWO-DIMENSIONAL MULTIGRID WITH TRIANGLES
    !Ideas:
    !We have subdivided the triangles between up ones and down ones.
    !We need to store 4 diferent values for each node:
    !                                               - Values of our data in that node (U)
    !                                               - Right hand side of the equation (F)
    !                                               - Residual of our function in each node (Res)
    !                                               - An auxiliar that can store for example in problems time
    !                                                 dependency the previous value of (U)
    !
    !As we said previously we will access from and up triangle the down triangles. Therefore
    !in order to make our algorithms coherents, we don't want to distinguish between inner triangles
    !and triangles in contact with the boundaries. Hence, we will store our boundaries/overlap
    !with the down triangles.
    !Following the previous example we can see very easily that only the up triangles have contact
    !with the boundaries/overlap and as our stencil only have four points(next three triangles)
    !The example before presented would be stored as follows:
    !
    !   UpUF:                  DownUF: (size = number of down triangles in the largest column + 3)
    !                                    (If we want both to have the same numberation
    !                                    this matrix must start with 0,0)
    !       |0|f|f|                     |-|-|-|
    !       |u|0|f|                     |b|b|f|
    !       |u|u|0|                     |-|1|-|
    !
    !   f -> Second part of the equation
    !   b -> Boundaries/overlap
    !   u -> Values in the nodes
    !   - -> Empty value
    !NOTE FOR NODES IN THE MIDDLE OF THE SIDES:
    !In this situation we will only work with up triangles and the stencils requiered are reduced to three.
    !The data will be stored as follows:
    !For each side there will be one array therefore, there will be three arrays to store all the nodes.
    !The overlap will be stored in UpUF in the position required, that are the diagonal, the colum(0) and
    !in the last row, each variable will use only two of this places.
    !NOTE:As here we are using vectors, the overlap must be stored with negative value.
    !
    !->Xp1,Xp2,Xp3: 1D array(2). The first component is the x position and the second is the y position.
    !               Xp2 is the coordinate with the highest x position
    !               Xp3 is the coordinate with the highest y position unless Xp2 has both the highest x and y,
    !                   then Xp3 will contain the coordinate with the second highest y
    !               Xp1 is the oder
    !->SizeUp/Down is the length of the largest column/row of the triangle. This is done because
    !in order to iterate we go from column 1 to the one previous to the diagonal
    !and from row 1 to Size + 1
    !
    !StencilCC is an array of doubles(3,3,7). Contains for each node the molecule of the multigrid
    !StencilSD is an array of doubles(3,3,3). Contains for each node the molecule of the multigrid
    !StencilMix is an array of doubles(3,3,:). Contains for each node the molecule of the multigrid
    !			The first three positions are for SD type nodes.
    !
    !->Visited is a logical to store if that level has been visited already. Only useful for W and F cycles.

    type, public :: Triangle
        double precision,allocatable, dimension(:,:,:) :: UpUF
        double precision,allocatable, dimension(:,:,:) :: UpResAux
        double precision,allocatable, dimension(:,:,:) :: DownUF
        double precision,allocatable, dimension(:,:,:) :: DownResAux
        double precision,allocatable, dimension(:,:,:) :: StencilCC
        double precision,allocatable, dimension(:,:,:) :: StencilSD
        double precision,allocatable, dimension(:,:,:) :: StencilMix
        double precision, dimension(2,3) :: X
        Integer :: SizeUp,SizeDown
        Logical :: Visited

    end type Triangle

    type, public :: multigrid_scaling
      real, allocatable :: detwei(:), nx(:,:,:), sdetwei(:,:)
    end type multigrid_scaling
    !This structure is developed to be used in an array, where the position inside the array is the name
    !of the triangle
    !Xp1 :: double(2). X,Y coordinate of Xp1 vertex
    !Xp2 :: double(2). X,Y coordinate of Xp2 vertex
    !Xp3 :: double(2). X,Y coordinate of Xp3 vertex
    !Method :: Is a number that indicates which smoother will be applied to that triangle
    !v1 :: Integer. Number of pre-smoothers
    !v2 :: Integer. Number of post-smoothers
    !weight :: Is the weight that shoud be applied to the sleected smoother
    !Neig(3) :: integer(3).Name of the neighbour triangles
    !				Goes from 1 wich is the l1 side to 3 the l3 side
    !Normals(3,2) :: double. Normal vectors of Mpos triangle.(Side,(Xcomponent, Ycomponent))
    !k_coef :: double. K coefficient, for example heat diffusion coefficient, by default = 1. This is the same coefficient
    !				as the one stored in Kcoord. This one is only for plotting purposes
    !Dir(3) :: logical. True if the direction of the sides are the same or not
    !					This is important when you have to share information between the boundaries.
    !					If Dir is true in that side, means that the first element of your side, coincides
    !					with the first element of the other triangle. If not then the first match with the
    !					last of the other one.
    !Tri :: type(Triangle)(:)
    type, public :: Mesh
        double precision, dimension(2,3) :: X
        real, dimension(2) :: center
        real, dimension(3) :: dc_unele
        real, dimension(3) :: dc_str_ele
        real :: str_area
        real, allocatable :: detwei(:) !ngi
        real, allocatable :: sdetwei(:,:) !sngi,nface
        ! it saves structured snorm not un_ele
        real, allocatable :: snorm(:,:,:) !(sngi,ndim,nface)
        real, allocatable :: nx(:,:,:) !( ngi, ndim, nloc )

        integer :: method, v1, v2
        integer, dimension(3) :: Neig, fNeig
        integer, allocatable :: S_nodes(:,:)
        logical, dimension(3) :: Dir
        real :: k_coef
        real, allocatable:: t_overlap(:,:) ! (2**n * nloc, nface)
        real, allocatable:: t_overlap_old(:,:) ! (2**n * nloc, nface)
        real, allocatable:: u_overlap(:,:,:) ! (ndim, 2**n * nloc, nface)
        real, allocatable:: u_ele(:,:,:) ! (ndim, nloc, totele_str)
        integer :: region_id
        integer, allocatable :: s_ele(:,:) ! global surface element of my neig relative to my face (2**n_split, nface)

        !Both arrays have as dimension the number of levels
        type(Triangle), allocatable, dimension(:) :: Tri
        type(multigrid_scaling), allocatable, dimension (:) :: scaling_var
    end type Mesh



    type, public :: element_info
      integer, allocatable :: surf_ele(:,:), str_neig(:,:)
      real, allocatable :: x_loc(:,:)
    end type element_info



    type, public:: pointer_Mesh
      type(Mesh), pointer :: ptr
    end type

    ! @brief:: this is a class used semi-structured and multigrids
    type, public :: fields
      real, dimension(:,:,:), allocatable :: tnew, told, error
      real, allocatable :: residuale(:,:,:), RHS(:,:,:), source(:,:,:)
    end type fields




    ! @brief:: this is a class for global CSR mass and stiffness matrices
    ! g_iloc :: global row number
    ! g_iloc :: global column number
    type, public :: sparse
        integer, allocatable :: g_iloc(:)
        integer, allocatable :: g_jloc(:)
        doubleprecision, allocatable :: val(:)

    end type sparse





    ! @brief:: this is a class for global CSR mass and stiffness matrices
    type, public :: m_CSR
        integer :: ele_id
        double precision, allocatable :: values(:,:)
    end type m_CSR


    ! @brief:: it gives unstructured neighbours full data
    type, public :: neig_data
      integer, dimension(3)::  Nside
      integer, dimension(3) :: Npos
      integer, allocatable :: Nnodes(:)
    end type neig_data

    ! type, public :: Semi_Mesh
    !     double precision, dimension(2,3) :: X
    !     real, allocatable:: up(:,:,:), down(:,:,:), x_str(:,:,:)
    !     !Both arrays have as dimension the number of levels
    !     type(Triangle), allocatable, dimension(:) :: Tri
    ! end type Semi_Mesh


	!Stores 4 coordinates, and the K coeficient to apply to the triangles inside that square
    !Xc1 :: double(2). X,Y coordinate of X1 vertex
    !Xc2 :: double(2). X,Y coordinate of X2 vertex
    !Xc3 :: double(2). X,Y coordinate of X3 vertex
    !Xc3 :: double(2). X,Y coordinate of X3 vertex
    !K :: double. Coefficient to apply
    type, public :: Kcoord
        double precision, dimension(2) :: Xc1
        double precision, dimension(2) :: Xc2
        double precision, dimension(2) :: Xc3
        double precision, dimension(2) :: Xc4
        double precision :: k
    end type Kcoord

end module Structures
