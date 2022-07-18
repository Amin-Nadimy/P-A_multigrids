module LinearSolvers
	use structures
	use TriangleOps, only: printmatrix, ConstructorSD, exportMatrix
!	use Smoothers
	use MeshOps, only: getNeigDataMesh
	use generic, only: vector_by_Matrix, Matrix_by_vector, Matrix_by_Matrix, AreEqual
    implicit none

    private
!    public :: solver_TriSD!<-Introduce circular dependency
	public :: solver_gauss, solver_GaussP, solver_Thomas, fact_PLU, solver_PLU, solver_MeshCC
	public :: GenerateA_CC, GenerateB_CC, solver_TriCC, solver_Thomas_Thin, Solve2x2
	public :: GSsolver_MeshCC, GSsolver_MeshSD,getInverse, solver_BlockThomas, GSsolver_MeshMix
	!Subroutine 1:
	!This function solves a 2x2 system, explicitly, and therefore, faster.
    !mat(in) :: double(2,2)
    !b(in) :: double(2)
    !Returns :: double(2). Solution of the equation system.
    !Subroutine 2:
	!This function solves a 2x2 system, explicitly.
	!Faster than subroutine 1 for loops with the same mat.
    !mat(in) :: double(2,2)
    !b(in) :: double(2)
    !det(in) :: double. The determinant of the mat matrix.
    !Returns :: double(2). Solution of the equation system.
	interface Solve2x2
	    module procedure Solve2x2A
	    module procedure Solve2x2B
	end interface Solve2x2


contains

    !Solves a linear equation system by gauss without pivoting
    !The input matrix must be defined as (Column, Row)
    !mat(in) :: double(:,:)
    !b(in) :: double(:)
    !Returns :: double(:). Solution of the equation system.
    !Note all the input vectors and matrices must have the same dimension
    function solver_gauss(mat, b)
        implicit none
        !Global variables
        double precision, intent(in), dimension(:) :: b
        double precision, intent(in), dimension(:,:) :: mat
        !Local variables
        double precision, dimension(size(b,1)) :: solver_gauss, b2
        double precision, dimension(size(b,1),size(b,1)) :: mat2
        double precision :: piv, idet
        integer :: j, i, k

        solver_gauss = 0d0
        b2 = b
        mat2 = transpose(mat)

        !Si el numero de variables es 1 o 2 sabemos la solución y no hace falta calcularla:
        if(size(b2,1)==1)then
            solver_gauss(1)=b2(1)/mat2(1,1)
            return
        else if(size(b2,1)==2)then
            idet = 1.d0/(mat2(1,1)*mat2(2,2)-mat2(1,2)*mat2(2,1))
            solver_gauss(1)=idet*(mat2(2,2)*b2(1)-mat2(1,2)*b2(2))
            solver_gauss(2)=idet*(-mat2(2,1)*b2(1)+mat2(1,1)*b2(2))
            return
        end if

        ! Hacemos ceros bajo la diagonal:
        do i=1,size(b2,1)
            do j=i+1,size(b2,1)
                piv = mat2(j,i)/mat2(i,i)
                do k=i+1,size(b2,1)
                    mat2(j,k) = mat2(j,k) - piv*mat2(i,k)
                end do
                b2(j) = b2(j) - piv*b2(i)
            end do
        end do

        !Despejamos las incógnitas:
        do i=size(b2,1),1,-1
            solver_gauss(i) = b2(i)
            do j=i+1,size(b2,1)
                solver_gauss(i) = solver_gauss(i) - mat2(i,j)*solver_gauss(j)
            end do
            solver_gauss(i) = solver_gauss(i)/mat2(i,i)
        end do
    end function solver_gauss

	!This function solves a 2x2 system, explicitly, and therefore, faster.
    !mat(in) :: double(2,2)
    !b(in) :: double(2)
    !Returns :: double(2). Solution of the equation system.
	function Solve2x2A(mat,b)
        implicit none
        !Global variables
        double precision, intent(in), dimension(2) :: b
        double precision, intent(in), dimension(2,2) :: mat
        !Local variables
        double precision, dimension(2):: Solve2x2A
        double precision :: aux

        aux = mat(1,1)*mat(2,2) - mat(2,1) * mat(1,2)

        Solve2x2A(1) = ( b(1) * mat(2,2) - b(2) * mat(2,1) ) / aux
        Solve2x2A(2) = ( b(2) * mat(1,1) - b(1) * mat(1,2) ) / aux

	end function Solve2x2A


	!This function solves a 2x2 system, explicitly, and therefore, faster.
    !mat(in) :: double(2,2)
    !b(in) :: double(2)
    !det(in) :: double. The determinant of the mat matrix.
    !Returns :: double(2). Solution of the equation system.
	function Solve2x2B(mat,b, det)
        implicit none
        !Global variables
        double precision, intent(in), dimension(2) :: b
        double precision, intent(in), dimension(2,2) :: mat
        double precision :: det
        !Local variables
        double precision, dimension(2):: Solve2x2B

        Solve2x2B(1) = ( b(1) * mat(2,2) - b(2) * mat(2,1) ) / det
        Solve2x2B(2) = ( b(2) * mat(1,1) - b(1) * mat(1,2) ) / det

	end function Solve2x2B
    !Solves a linear equation system by gauss with pivoting
    !The input matrix must be defined as (Column, Row)
    !mat(in) :: double(:,:)
    !b(in) :: double(:)
    !Returns :: double(:). Solution of the equation system.
    !Note all the input vectors and matrices must have the same dimension
	function solver_GaussP(mat, b)
        implicit none
        !Global variables
        double precision, intent(in), dimension(:) :: b
        double precision, intent(in), dimension(:,:) :: mat
        !Local variables
		double precision :: aux, maxi
		integer :: i, k, n
		integer, dimension(1) :: rowmax
		double precision, dimension(size(b,1),size(b,1)) :: A
		double precision, dimension(size(b,1)) :: auxRow, solver_GaussP
	    n = size(b,1)

	    !If size==2 use a fastest method
		if (n==2) then
		    solver_GaussP = Solve2x2A(mat,b)
		    return
		end if
		!If mix type problem
!		if (n==4) then
!		    if (AreEqual(0d0,mat(4,4)).and.AreEqual(0d0,mat(2,1)).and.AreEqual(0d0,mat(3,1)).and.AreEqual(0d0,mat(3,2))&
!		    	.and.AreEqual(0d0,mat(1,2)).and.AreEqual(0d0,mat(1,3)).and.AreEqual(0d0,mat(2,3))) then
!		        solver_GaussP(4) = ( b(1)+b(2)+b(3)-b(4) )/( mat(4,1) + mat(4,2) + mat(4,3) )
!		        solver_GaussP(1) = b(1) + solver_GaussP(4) * mat(4,1)
!		        solver_GaussP(2) = b(2) + solver_GaussP(4) * mat(4,2)
!		        solver_GaussP(3) = b(3) + solver_GaussP(4) * mat(4,3)
!		    end if
!		end if
		!Prepare variables
		A = mat
		solver_GaussP = b

		do i=1,n

			!Pivot only if it is worthy
	        if (A(i,i)/maxval(abs(A(i,i:n)))<1d-3) then
			    maxi = 0d0
			    rowmax(1) = i
			    !Before calling the Gauss method we will seek for
			    !a better row to pivote with
				rowmax =  (i-1) + maxloc(abs(A(i,i:n)))
		        !Now we modify the A matrix and the right-hand side vector
		        auxRow =  A(:,i)
		        A(:,i) = A(:,rowmax(1))
		        A(:,rowmax(1)) = auxRow
		        maxi = solver_GaussP(i)
		        solver_GaussP(i) = solver_GaussP(rowmax(1))
		        solver_GaussP(rowmax(1)) = maxi
	        end if

		    !Gauss Method
		    do k = i+1, n
		        aux = (A(i,k)/A(i,i))
		        A(:,k) = A(:,k) - A(:,i) * aux
		        solver_GaussP(k) = solver_GaussP(k) - solver_GaussP(i) * aux
		    end do

		end do

		solver_GaussP = Backward(A,solver_GaussP)
	end function solver_GaussP

	!Returns the condition number of the input matrix
    !mat(in) :: double(:,:)
    !Returns :: double
	doubleprecision function getConditionNumber(mat)
        implicit none
        !Global variables
        double precision, intent(in), dimension(:,:) :: mat

        getConditionNumber = maxval(abs(mat))*maxval(abs(getInverse(mat)))

	end function getConditionNumber

    !Returns the inverse for the input triangle by
    !Gauss-Jordan with partial pivoting
    !mat(in) :: double(:,:)
    !Inverse(out) :: double(:). Inverse
    !Piv(out) :: double(:,:). The pivoting matrix
    !Note all the input vectors and matrices must have the same dimension
	subroutine getInverse_Piv(mat, Inverse, p)
        implicit none
        !Global variables
        double precision, intent(in), dimension(:,:) :: mat
        double precision, intent(out), dimension(:,:) :: p, Inverse
        !Local variables
		double precision :: aux
		integer :: i, k, n
		integer, dimension(1) :: rowmax
		double precision, dimension(size(mat,1)*2,size(mat,1)) :: A
		double precision, dimension(size(mat,1)*2) :: auxRow
		double precision, dimension(size(mat,1)) :: auxRow2

		!Prepare data
	    n = size(mat,1)
		A(1:n,:) = mat
		A(n+1:2*n,:) = 0d0
		!Introduce the identity matrix
		do i = 1, n
		    A(n+i,i) = 1d0
		end do

		P = 0d0
		do i = 1, n
		    P(i,i) = 1d0
		end do

		do i=1,n
		    rowmax(1) = i

		    !Before calling the Gauss method we will seek for
		    !a better row to pivote with
			rowmax =  (i-1) + maxloc(abs(A(i,i:n)))

	        !Now we modify the A matrix
	        auxRow =  A(:,i)
	        A(:,i) = A(:,rowmax(1))
	        A(:,rowmax(1)) = auxRow

			!And the Permutation matrix
	        auxRow2 =  P(:,i)
	        P(:,i) = P(:,rowmax(1))
	        P(:,rowmax(1)) = auxRow2

		    !Gauss Method
		    do k = i+1, n
		        A(:,k) = A(:,k) - A(:,i) * A(i,k)/A(i,i)
		    end do
		    !Normalize
			A(:,i) = A(:,i)/A(i,i)
		end do

		do i = n, 1,-1
		    !Gauss Method
		    do k = i-1, 1,-1
		        A(:,k) = A(:,k) - A(:,i) * A(i,k)
		    end do
		end do

		Inverse = A(size(mat,1)+1:size(A,1),:)

	end subroutine getInverse_Piv

    !Returns the inverse for the input triangle by
    !Gauss-Jordan with partial pivoting
    !mat(in) :: double(:,:)
    !Returns :: double(:). Inverse
    !Note all the input vectors and matrices must have the same dimension
	function getInverse(mat)
        Implicit none
        !Global variables
        double precision, intent(in), dimension(:,:) :: mat
       !Local variables
		double precision :: aux
		integer :: i, k, n
		integer, dimension(1) :: rowmax
		double precision, dimension(size(mat,1)*2,size(mat,2)) :: A
		double precision, dimension(size(mat,1),size(mat,2)) :: getInverse!, P
		double precision, dimension(size(mat,1)*2) :: auxRow
!		double precision, dimension(size(mat,1)) :: auxRow2

	    n = size(mat,1)
	    if (n==2) then!Faster formula
	    	aux = mat(1,1)*mat(2,2)-mat(1,2)*mat(2,1)
	        getInverse(1,1) = mat(2,2)/aux
	        getInverse(2,2) = mat(1,1)/aux
	        getInverse(1,2) = -mat(1,2)/aux
	        getInverse(2,1) = -mat(2,1)/aux
	        return
	    end if

		!Prepare data
		A(1:n,:) = mat
		A(n+1:2*n,:) = 0d0
		!Introduce the identity matrix
		do i = 1, n
		    A(n+i,i) = 1d0
		end do

		do i=1,n

			!Pivote only if it is worthy
			if (A(i,i)/maxval(abs(A(i,i:n)))<1d-3) then
			    rowmax(1) = i
			    !Before calling the Gauss method we will seek for
			    !a better row to pivote with
				rowmax =  (i-1) + maxloc(abs(A(i,i:n)))
		        !Now we modify the A matrix
		        auxRow =  A(:,i)
		        A(:,i) = A(:,rowmax(1))
		        A(:,rowmax(1)) = auxRow
			end if

		    !Gauss Method
		    do k = i+1, n
		        A(:,k) = A(:,k) - A(:,i) * A(i,k)/A(i,i)
		    end do
		    !Normalize
			A(:,i) = A(:,i)/A(i,i)
		end do

		do i = n, 1,-1
		    !Gauss Method
		    do k = i-1, 1,-1
		        A(:,k) = A(:,k) - A(:,i) * A(i,k)
		    end do
		end do

		getInverse = A(size(mat,1)+1:size(A,1),:)

	end function getInverse

	!Thomas solver for block-tridiagonal systems
	!lower(inout) :: double(:,:,:)
	!upper(inout) :: double(:,:,:)
	!central(inout) :: double(:,:,:)
	!b(inout) :: double(:,:)
	!returns :: double(:,:) Solution of the system
	!Note: First component is the position, and the other two are the block, for example
	!for 2x2 blocks: lower(:,2,2) and b(:,2)
	function solver_BlockThomas(lower,upper,central,b)
        implicit none
        !Global variables
        double precision, intent(inout), dimension(:,:) :: b
        double precision, intent(inout), dimension(:,:,:) :: lower,upper,central
        !Local variables
		integer :: n, j
		double precision, dimension(size(b,1), size(b,2)) :: solver_BlockThomas
		double precision, dimension(size(b,2),size(b,2)) :: aux
		!Prepare data
	    n = size(b,1)

		do j = 2, n
		    !Component of L
		    lower(j,:,:) = Matrix_by_Matrix(lower(j,:,:),getInverse(central(j-1,:,:)))
		    !Diagonal component of U
		    central(j,:,:) = central(j,:,:) - Matrix_by_Matrix(lower(j,:,:),upper(j-1,:,:))
		    !Right-hand side modification
		    b(j,:) = b(j,:) - matrix_by_vector(lower(j,:,:),b(j-1,:))
		end do

		!Solver
		solver_BlockThomas(n,:) = solver_GaussP(central(n,:,:), b(n,:))

		do j = n-1, 1, -1
			solver_BlockThomas(j,:) = solver_GaussP(central(j,:,:), b(j,:)- matrix_by_vector(upper(j,:,:),solver_BlockThomas(j+1,:)))
		end do

	end function solver_BlockThomas

	!Thomas solver for tridiagonal systems
	!mat(in) :: double(:,:)
	!b(in) :: double(:)
	!returns :: double(:) Solution of the system
	!NOTE: (Columns, Rows)
	function solver_Thomas(mat, b)
        implicit none
        !Global variables
        double precision, intent(inout), dimension(:) :: b
        double precision, intent(in), dimension(:,:) :: mat
        !Local variables
		integer :: i, n
		double precision, dimension(size(b,1),size(b,1)) :: A, L
		double precision, dimension(size(b,1)) :: solver_Thomas,aux
	    n = size(b,1)
		A = mat
		L = 0d0
		do i = 1, n
		    L(i,i) = 1d0
		end do

		do i = 2, n
		    !Component of L
		    A(i-1,i) = A(i-1,i) / A(i-1,i - 1)
		    !Diagonal component of U
		    A(i,i) = A(i,i) - A(i,i-1) * A(i-1,i)
		end do

		do i = 2, n
	        L(i-1,i) = A(i-1,i)
    		A(i-1,i) = 0
		end do

		aux = b
		!Forward solver
		do i = 2, n
		    aux(i) = (aux(i) - aux(i-1)*L(i-1,i))
		end do

		!Backward solver
		solver_Thomas(n) = aux(n)/A(n,n)
		do i = n-1, 1, -1
		    solver_Thomas(i) = (aux(i) - solver_Thomas(i+1)*A(i+1,i)) / A(i,i)
		end do

	end function solver_Thomas

	!Thomas solver for tridiagonal systems
	!However this method use a [:,3] matrix
	!A(inout) :: double(:,3)
	!b(inout) :: double(:)
	!returns :: double(:) Solution of the system
	!EXAMPLE:
	!Usual matrix:
	!	level 1->	1 2 0 0
	!	level 2->	3 4 5 0
	!	level 3->	0 6 7 8
	!	level 4->	0 0 9 10
	!This matrix must be transformed to this:
	!	level 1->	0 3 6 9
	!	level 2->	1 4 7 10
	!	level 3->	2 5 8 0
	function solver_Thomas_Thin(A, b)
        implicit none
        !Global variables
        double precision, intent(inout), dimension(:) :: b
        double precision, intent(inout), dimension(:,:) :: A
        !Local variables
		integer :: k, n, j
		double precision, dimension(size(b,1)) :: solver_Thomas_Thin
		!Prepare data
	    n = size(b,1)
		do j = 2, n
		    !Component of L
		    A(j,1) = A(j,1) / A(j - 1,2)
		    !Diagonal component of U
		    A(j,2) = A(j,2) - A(j-1,3) * A(j,1)
		end do

		!Forward solver
		do k = 2, n
		    b(k) = (b(k) - b(k-1)*A(k,1))
		end do

		!Backward solver
		solver_Thomas_Thin(n) = b(n)/A(n,2)
		do k = n-1, 1, -1
		    solver_Thomas_Thin(k) = (b(k) - solver_Thomas_Thin(k+1)*A(k,3)) / A(k,2)
		end do

	end function solver_Thomas_Thin

    !Solves a linear equation system by gauss with pivoting
    !The input matrix must be defined as (Column, Row)
    !P(in) :: double(:,:)
    !U(in) :: double(:,:). Matrix cotaining U and L
    !b(in) :: double(:)
    !Returns :: double(:). Solution of the equation system.
    !Note all the input vectors and matrices must have the same dimension
    !FUNCTION TO BE USED AFTER THE fact_PLU
	function solver_PLU(P,U,b)
        implicit none
        !Global variables
        double precision, intent(inout), dimension(:) :: b
        double precision, intent(in), dimension(:,:) :: P,U
        !Local variables
		double precision, dimension(size(b,1)) :: solver_PLU

	    solver_PLU = Backward(U, ForwardPLU(U,matmul(P,b)))

	end function solver_PLU

    !Does the PLU factorization
    !The input matrix must be defined as (Column, Row)
    !mat(in) :: double(:,:)
    !P(out) :: double(:) (Transposed, to save some calculations in solver_PLU)
    !U(out) :: double(:).(Contains the L and the U)
    !Note: all the input vectors and matrices must have the same dimension
    !	   In order to solve an equation system, use solver_PLU after this method
    !FUNCTION TO BE USED BEFORE THE solver_PLU ONLY!!
	subroutine fact_PLU(mat, P, U)
        implicit none
        !Global variables
        double precision, intent(in), dimension(:,:) :: mat
        double precision, intent(inout), dimension(:,:) :: P, U!, L
        !Local variables
		double precision :: aux, maxi
		integer :: i, j, k, n
		integer, dimension(1) :: rowmax
		double precision, dimension(size(mat,1)) :: auxRow
	    n = size(mat,1)
		U = mat

		P = 0d0
		do j = 1, n
		    P(j,j) = 1d0
		end do

		do i=1,n

			!Pivot only if it is worthy
	        if (U(i,i)/maxval(abs(U(i,i:n)))<1d-3) then
			    maxi = 0d0
			    rowmax(1) = i
			    !We will seek for a better row to pivote with
				rowmax =  (i-1) + maxloc(abs(U(i,i:n)))

		        !Now we modify the matrix
		        auxRow =  U(:,i)
		        U(:,i) = U(:,rowmax(1))
		        U(:,rowmax(1)) = auxRow

		        auxRow =  P(:,i)
		        P(:,i) = P(:,rowmax(1))
		        P(:,rowmax(1)) = auxRow
			end if

		    !LU method
		    do k = i+1, n
		        aux = (U(i,k)/U(i,i))
		        do j = i, n
		            U(j,k) = U(j,k) - U(j,i) * aux
		        end do
		        U(i,k) = aux
		    end do

		end do

		!Finally transpose because in Fortran we work
		!in (columns, row), not the other way around
		P = transpose(P)

	end subroutine fact_PLU
	!Solves the lower triangular system Ax = b
	!mat(in) :: double(:,:)
	!b(in) :: double(:)
	!Return :: double(size(b)). Solution of the system
	function Forward(mat,b)
		Implicit none
		!Global variables
		double precision, intent(in), dimension(:,:) :: mat
		double precision, intent(in), dimension(:) :: b
		!Local variables
		double precision, dimension(size(b,1)) :: Forward
		double precision :: aux
		integer :: i, j

		do i = 1, size(b,1)
			aux = 0d0
		    do j = 1, i-1
				aux = aux + forward(j)*mat(j,i)
		    end do
		    forward(i) = (b(i) - aux) / mat(i,i)
		end do
	end function Forward

	!Solves the upper triangular system Ax = b
	!mat(in) :: double(:,:)
	!b(in) :: double(:)
	!Return :: double(size(b)). Solution of the system
	function Backward(mat,b)
		Implicit none
		!Global variables
		double precision, intent(in), dimension(:,:) :: mat
		double precision, intent(in), dimension(:) :: b
		!Local variables
		double precision, dimension(size(b,1)) :: Backward
		double precision :: aux
		integer :: i, j, n

		n = size(b,1)
		do i = n, 1, -1
			aux = 0d0
		    do j = n, i+1, -1
				aux = aux + Backward(j)*mat(j,i)
		    end do
		    Backward(i) = (b(i) - aux) / mat(i,i)
		end do
	end function Backward

	!Solves the lower triangular system Ax = b
	!mat(in) :: double(:,:)
	!b(in) :: double(:)
	!Return :: double(size(b)). Solution of the system
	!THIS METHOD IS ONLY TO USE AFTER THE PLU,
	!BECAUSE THE DIAGONALS WILL BE TAKEN AS 1
	function ForwardPLU(mat,b)
		Implicit none
		!Global variables
		double precision, intent(in), dimension(:,:) :: mat
		double precision, intent(in), dimension(:) :: b
		!Local variables
		double precision, dimension(size(b,1)) :: ForwardPLU
		double precision :: aux
		integer :: i, j

		do i = 1, size(b,1)
			aux = 0d0
		    do j = 1, i-1
				aux = aux + ForwardPLU(j)*mat(j,i)
		    end do
		    ForwardPLU(i) = (b(i) - aux)
		end do
	end function ForwardPLU

	!!Just to solve the unstructured mesh!!
	!meshL(inout) :: type (mesh)(:)
	!level(in) :: integer. The level to solve. It MUST be the deepest
	!it(in) :: integer. Number of iterations. By default the half of the total number of triangles
	!It uses the iterative Gauss-seidel method to solve the COARSEST grid
	subroutine GSsolver_MeshCC(meshL, level, it)
		Implicit none
		!Global variables
		type(mesh), intent(inout), dimension(:) :: meshL
		integer, intent(in) :: level
		integer, optional, intent(in) :: it
		!Local variables
		integer :: k, iteration, m,a
		double precision :: aux
		!Prepare data
		iteration = size(meshL)
		if (present(it)) iteration = it
		a = size(meshL(1)%Tri(level)%UpUF,3)
		do m = 1, iteration

			do k = 1, size(meshL)
				aux = meshL(k)%Tri(level)%StencilCC(2,2,7)*meshL(k)%Tri(level)%UpUF(1,2,a)
			    !L1 side
		        if (meshL(k)%Neig(1)/=0) then
		            aux = aux + meshL(k)%Tri(level)%StencilCC(1,2,7)&
		            	*meshL(meshL(k)%Neig(1))%Tri(level)%UpUF(1,2,a)
		        end if
			    !L2 side
		        if (meshL(k)%Neig(2)/=0) then
		            aux = aux + meshL(k)%Tri(level)%StencilCC(2,3,7)&
		            	*meshL(meshL(k)%Neig(2))%Tri(level)%UpUF(1,2,a)
		        end if
		        !L3 side
		        if (meshL(k)%Neig(3)/=0) then
		            aux = aux + meshL(k)%Tri(level)%StencilCC(3,2,7)&
		            	*meshL(meshL(k)%Neig(3))%Tri(level)%UpUF(1,2,a)
		        end if

		        meshL(k)%Tri(level)%UpUF(1,2,a) = meshL(k)%Tri(level)%UpUF(1,2,a) &
		        		+ (meshL(k)%Tri(level)%UpUF(2,1,a)-aux)/meshL(k)%Tri(level)%StencilCC(2,2,7)
	        end do

		end do

	end subroutine GSsolver_MeshCC


	!!Just to solve the unstructured mesh!!
	!meshL(inout) :: type (mesh)(:)
	!level(in) :: integer. The level to solve. It MUST be the deepest
	subroutine solver_MeshCC(meshL, level)
		Implicit none
		!Global variables
		type(mesh), intent(inout), dimension(:) :: meshL
		integer, intent(in) :: level
		!Local variables
		double precision, dimension(size(meshL),size(meshL)) :: A
		double precision, dimension(size(meshL)) :: b
		integer ::  k, f
		!Prepare data
		A = 0d0
		f = size(meshL(1)%Tri(level)%UpUF,3)
		!Create A and b
		do k = 1, size(meshL)
		    A(k,k) = meshL(k)%Tri(level)%StencilCC(2,2,7)
		    b(k) = meshL(k)%Tri(level)%UpUF(2,1,f)
		    !L1 side
	        if (meshL(k)%Neig(1)/=0) then
	            A(meshL(k)%Neig(1),k) = meshL(k)%Tri(level)%StencilCC(1,2,7)
	        end if
		    !L2 side
	        if (meshL(k)%Neig(2)/=0) then
	            A(meshL(k)%Neig(2),k) = meshL(k)%Tri(level)%StencilCC(2,3,7)
	        end if
	        !L3 side
	        if (meshL(k)%Neig(3)/=0) then
	            A(meshL(k)%Neig(3),k) = meshL(k)%Tri(level)%StencilCC(3,2,7)
	        end if
		end do

		!Solve and save data
		b = solver_GaussP(A,b)
		do k = 1, size(meshL)
		    meshL(k)%Tri(level)%UpUF(1,2,f) = b(k)
		end do

	end subroutine solver_MeshCC
	!!Just to solve the unstructured mesh!!
	!meshL(inout) :: type (mesh)(:)
	!level(in) :: integer. The level to solve. It MUST be the deepest
	!it(in) :: integer. Number of iterations. By default the half of the total number of triangles
	!It uses the iterative Gauss-seidel method to solve the COARSEST grid
	subroutine GSsolver_MeshSD(meshL, level, it)
		Implicit none
		!Global variables
		type(mesh), intent(inout), dimension(:) :: meshL
		integer, intent(in) :: level
		integer, optional, intent(in) :: it
		!Local variables
		integer :: k, iteration, m
		integer :: inicio, fin, inc, Npos, Nside
	    integer, dimension(3) :: Ndata
	    integer, dimension(2) :: Offset
		double precision :: aux
		!Prepare data
		iteration = size(meshL)/2
		if (present(it)) iteration = it

		do m = 1, iteration
			do k = 1, size(meshL)
			    !L1 side
				call getNeigDataMesh(meshL, k, 1, level, inicio, fin, inc, Npos, Nside, Ndata,.false., Offset)
		        if (Npos > 0) then
		        	!Calculate residual
			        aux = meshL(k)%Tri(Level)%UpUF(2,1,1) -&
			        ( meshL(k)%Tri(Level)%StencilSD(2,2,4) * meshL(k)%Tri(Level)%UpUF(1,2,1)&
			        + meshL(k)%Tri(Level)%StencilSD(2,3,4) * meshL(k)%Tri(Level)%UpUF(1,2,2)&
			        + meshL(k)%Tri(Level)%StencilSD(3,2,4) * meshL(k)%Tri(Level)%UpUF(1,2,3)&
			        - meshL(k)%Tri(Level)%StencilSD(1,2,4) * meshL(Npos)%Tri(Level)%UpUF(1,2,Ndata(2))&
			        - meshL(k)%Tri(Level)%StencilSD(2,1,4) * meshL(Npos)%Tri(Level)%UpUF(1,2,Ndata(3)))

                    meshL(k)%Tri(level)%UpUF(1,2,1) = meshL(k)%Tri(level)%UpUF(1,2,1)&
                    	 + (aux)/meshL(k)%Tri(level)%StencilSD(2,2,4)
		        end if
			    !L2 side
		    	call getNeigDataMesh(meshL, k, 2, level, inicio, fin, inc, Npos, Nside, Ndata,.false., Offset)
		        if (Npos > 0) then

		        	aux = meshL(k)%Tri(level)%UpUF(2,1,2) -&
		        	( meshL(k)%Tri(level)%StencilSD(2,2,5) * meshL(k)%Tri(level)%UpUF(1,2,2)&
		        	+ meshL(k)%Tri(level)%StencilSD(2,1,5) * meshL(k)%Tri(level)%UpUF(1,2,1)&
		        	+ meshL(k)%Tri(level)%StencilSD(3,2,5) * meshL(k)%Tri(level)%UpUF(1,2,3)&
		        	- meshL(k)%Tri(level)%StencilSD(1,2,5) * meshL(Npos)%Tri(level)%UpUF(1,2,Ndata(2))&
		        	- meshL(k)%Tri(level)%StencilSD(2,3,5) * meshL(Npos)%Tri(level)%UpUF(1,2,Ndata(3)))

                    meshL(k)%Tri(level)%UpUF(1,2,2) = meshL(k)%Tri(level)%UpUF(1,2,2)&
                    	 + (aux)/meshL(k)%Tri(level)%StencilSD(2,2,5)
		        end if
		        !L3 side
		    	call getNeigDataMesh(meshL, k, 3, level, inicio, fin, inc, Npos, Nside, Ndata,.false., Offset)
		        if (Npos > 0) then
		        	aux = meshL(k)%Tri(level)%UpUF(2,1,3) -&
		        	( meshL(k)%Tri(level)%StencilSD(2,2,6) * meshL(k)%Tri(level)%UpUF(1,2,3)&
		        	+ meshL(k)%Tri(level)%StencilSD(1,2,6) * meshL(k)%Tri(level)%UpUF(1,2,1)&
		        	+ meshL(k)%Tri(level)%StencilSD(2,3,6) * meshL(k)%Tri(level)%UpUF(1,2,2)&
		        	- meshL(k)%Tri(level)%StencilSD(3,2,6) * meshL(Npos)%Tri(level)%UpUF(1,3,Ndata(2))&
		        	- meshL(k)%Tri(level)%StencilSD(2,1,6) * meshL(Npos)%Tri(level)%UpUF(2,3,Ndata(3)))

                    meshL(k)%Tri(level)%UpUF(1,2,3) = meshL(k)%Tri(level)%UpUF(1,2,3)&
                    	 + (aux)/meshL(k)%Tri(level)%StencilSD(2,2,6)
		        end if
	        end do
		end do

	end subroutine GSsolver_MeshSD

	!!Just to solve the unstructured mesh!!
	!meshL(inout) :: type (mesh)(:)
	!level(in) :: integer. The level to solve. It MUST be the deepest
	!it(in) :: integer. Number of iterations. By default the half of the total number of triangles
	!It uses the iterative Gauss-seidel method to solve the COARSEST grid
	subroutine GSsolver_MeshMix(meshL, level, it)
		Implicit none
		!Global variables
		type(mesh), intent(inout), dimension(:) :: meshL
		integer, intent(in) :: level
		integer, optional, intent(in) :: it
		!Local variables
		integer :: Mpos, iteration, m,Npos
		double precision, dimension(4,4) :: Mat
		double precision, dimension(4) :: b
		!Prepare data
		iteration = size(meshL)*5
		if (present(it)) iteration = it
		mat = 0d0

		do m = 1, iteration
			do Mpos = 1, size(meshL)
		    	!Mix type problem
				mat(1,1) = 1d0
				mat(2,2) = 1d0
				mat(3,3) = 1d0
				mat(1,4) = meshL(Mpos)%Tri(level)%StencilSD(2,2,1)
				mat(2,4) = meshL(Mpos)%Tri(level)%StencilSD(2,3,1)
				mat(3,4) = meshL(Mpos)%Tri(level)%StencilSD(3,2,1)
				mat(4,1) = meshL(Mpos)%Tri(level)%StencilCC(1,2,1)
				mat(4,2) = meshL(Mpos)%Tri(level)%StencilCC(2,3,2)
				mat(4,3) = meshL(Mpos)%Tri(level)%StencilCC(3,2,3)
				b = 0d0
			    !L1 side
				Npos = meshL(Mpos)%Neig(1)
		        if (Npos > 0) then
		        	!Calculate residual
		    		b(1) = meshL(Mpos)%Tri(level)%UpUF(2,1,1) &
		    			- (meshL(Mpos)%Tri(level)%UpUF(1,2,1) - (meshL(Npos)%Tri(level)%UpUF(1,2,4)&
		    			-meshL(Mpos)%Tri(level)%UpUF(1,2,4))*meshL(Mpos)%Tri(level)%StencilCC(1,2,1) )
		        end if
			    !L2 side
		    	Npos = meshL(Mpos)%Neig(2)
		        if (Npos > 0) then
					b(2) = meshL(Mpos)%Tri(level)%UpUF(2,1,2) &
		    			- (meshL(Mpos)%Tri(level)%UpUF(1,2,2) - (meshL(Npos)%Tri(level)%UpUF(1,2,4) &
		    			- meshL(Mpos)%Tri(level)%UpUF(1,2,4))*meshL(Mpos)%Tri(level)%StencilCC(2,3,2) )
		        end if
		        !L3 side
		        Npos = meshL(Mpos)%Neig(3)
		        if (Npos > 0) then
					b(3) = meshL(Mpos)%Tri(level)%UpUF(2,1,3) &
		    			- (meshL(Mpos)%Tri(level)%UpUF(1,2,3) - (meshL(Npos)%Tri(level)%UpUF(1,2,4) &
		    			- meshL(Mpos)%Tri(level)%UpUF(1,2,4))*meshL(Mpos)%Tri(level)%StencilCC(3,2,3) )
		        end if
				b(4) = meshL(Mpos)%Tri(level)%UpUF(2,1,4) - (meshL(Mpos)%Tri(level)%StencilSD(2,2,1)&
						*meshL(Mpos)%Tri(level)%UpUF(1,2,1) + meshL(Mpos)%Tri(level)%StencilSD(2,3,1)&
						*meshL(Mpos)%Tri(level)%UpUF(1,2,2) + meshL(Mpos)%Tri(level)%StencilSD(3,2,1)&
						*meshL(Mpos)%Tri(level)%UpUF(1,2,3))

				b = solver_GaussP(mat,b)
				meshL(Mpos)%Tri(level)%UpUF(1,2,:) = meshL(Mpos)%Tri(level)%UpUF(1,2,:) + b(:)
	        end do
		end do

	end subroutine GSsolver_MeshMix
	!Solves analytically triangle tri (Cell centered)
	!IMPORTANT: It works if there are no boundaries, or if they are zero. This means that this can be used
	!with all the grids except the finest one
	!Uses the Gauss solver with partial pivoting
	!Tri(inout) :: type(Triangle)
	!The numeration follows the L3 side. Starting with Xp1 as 1. As an example:
	!            /\
	!		    /  \
	!		   /_4__\
	!		  /\    /\  L3
	!		 /  \ 3/  \
	!	Xp1	/_1__\/_2__\
	subroutine solver_TriCC(Tri)
		Implicit none
		!Global variables
		type(Triangle), intent(inout) :: Tri
		double precision, dimension(Tri%SizeUp**2) :: b
		!Local variables
		integer :: j,i,k, a

		a = size(Tri%UpUF,3)

		if (Tri%SizeUp==1) then
			Tri%UpUF(1,2,a) = Tri%UpUF(2,1,a)/Tri%StencilCC(2,2,7)
	    else
		    !Firstly I solve the system
			b = GenerateB_CC(Tri)
			b = solver_GaussP(GenerateA_CC(Tri),b)

			!Now I introduce the data into Tri
			k = 1
			!Up triangles
			do i = 2, Tri%SizeUp+1
			    do j = 1, i-1
			        Tri%UpUF(j,i,a) = b(k)
			        k = k + 2
			    end do
			    k = k - 1
			end do
			!Down triangles
			k = 3
			do i = 2, Tri%SizeDown+1
			    do j = 1, i-1
			        Tri%DownUF(j,i,a) = b(k)
			        k = k + 2
			    end do
			    k = k + 1
			end do
		end if

	end subroutine solver_TriCC
	!Solves analytically triangle tri (Side centered)
	!IMPORTANT: It works if there are no boundaries, or if they are zero. This means that this can be used
	!with all the grids except the finest one
	!Uses the Gauss solver with partial pivoting
	!Tri(inout) :: type(Triangle)
	!The numeration follows the L2 side. Starting with Xp1 as 1. The nodes of the L1 and L3 side,
	!are one row, and the L2 nodes are the next row
	!            /\
	!		    /  \
	!		   /_3__\
	!	L1	  /\    /\   L3
	!		 /  1  2  \
	!	Xp1	/____\/____\
	!			L2
!	subroutine solver_TriSD(Tri)
!		Implicit none
!		!Global variables
!		type(Triangle), intent(inout) :: Tri
!		double precision, dimension(((Tri%SizeUp/2)*(Tri%SizeUp+1)-Tri%SizeUp)*3) :: b
!		!Local variables
!		integer :: j,i,k, TamRow
!
!		!Solve linear system
!		b = GenerateB_SD(Tri)
!		b = solver_GaussP(GenerateA_SD(Tri%Xp1,Tri%Xp2,Tri%Xp3,Tri%SizeUp*Tri%SizeUp),b)
!
!		!Introduce the solution into Tri
!		k = 1
!		!First L3 nodes
!		TamRow = Tri%SizeUp-1!Size of the elements of the L2 kind
!		do j = 1, Tri%SizeUp-1
!		    do i = j+1, Tri%SizeUp
!				Tri%UpUF(j,i,3) = b(k)
!				k = k + 2
!		    end do
!		    k = k + TamRow
!		    TamRow = TamRow - 1
!		end do
!		k = 2
!		!Secondly L1 nodes
!		TamRow = Tri%SizeUp-1!Size of the elements of the L2 kind
!		do j = 1, Tri%SizeUp-1
!		    do i = j+2, Tri%SizeUp+1
!				Tri%UpUF(j,i,1) = b(k)
!				k = k + 2
!		    end do
!		    k = k + TamRow
!		    TamRow = TamRow - 1
!		end do
!		k = Tri%SizeUp * 2 - 1
!		!Finally L2 nodes
!		TamRow = Tri%SizeUp * 2 - 4
!		do j = 2, Tri%SizeUp
!		    do i = j+1, Tri%SizeUp+1
!				Tri%UpUF(j,i,2) = b(k)
!				k = k + 1
!		    end do
!		    k = k + TamRow
!		    TamRow = TamRow - 2
!		end do
!
!	end subroutine solver_TriSD



	!Creates the matrix A for the equation system asociated to our stencil
	!Tri(in) :: Type(Triangle)
	!Returns :: double (:,:)
	!The numeration follows the L3 side. Starting with Xp1 as 1. As an example:
	!            /\
	!		    /  \
	!		   /_4__\
	!		  /\    /\  L3
	!		 /  \ 3/  \
	!	Xp1	/_1__\/_2__\
	function GenerateA_CC(Tri)
		Implicit none
		!Global variables
		type(Triangle), intent(in) :: Tri
		!Local variables
		double precision, dimension(Tri%SizeUp**2,Tri%SizeUp**2) :: GenerateA_CC
		integer ::TamCol, aux, k, m, aux2
		!Prepare data
		GenerateA_CC = 0d0

		!First Xp1
		GenerateA_CC(1,1) = Tri%StencilCC(2,2,1)
		GenerateA_CC(3,1) = Tri%StencilCC(3,2,1)

		TamCol = 3
		k = 2
		do aux = 2, Tri%SizeUp-1
		    !First L2
		    GenerateA_CC(k,k) = Tri%StencilCC(2,2,5)
		    GenerateA_CC(k+1,k) = Tri%StencilCC(1,2,5)
		    GenerateA_CC(k+1+TamCol,k) = Tri%StencilCC(3,2,5)
		    k = k + 1
		    !General situation, Up triangles
		    aux2 = k
		    k = k + 1
		    do m = 2, TamCol-3,2
		    	GenerateA_CC(k-1,k) = Tri%StencilCC(2,3,7)
		        GenerateA_CC(k,k) = Tri%StencilCC(2,2,7)
		        GenerateA_CC(k+1,k) = Tri%StencilCC(1,2,7)
		        GenerateA_CC(k+1+TamCol,k) = tri%StencilCC(3,2,7)
		        k = k + 2
		    end do
		    !General situation, down triangles
		    k = aux2
		    do m = 1, TamCol-2,2
		        GenerateA_CC(k-1,k) = tri%StencilCC(3,2,7)
		        GenerateA_CC(k,k) = Tri%StencilCC(2,2,7)
		        GenerateA_CC(k+1,k) = Tri%StencilCC(2,3,7)
		        GenerateA_CC(k+1-TamCol,k) = Tri%StencilCC(1,2,7)
		        k = k + 2
		    end do
		    k = k - 1
			!Finally L1
		    GenerateA_CC(k,k) = Tri%StencilCC(2,2,4)
		    GenerateA_CC(k-1,k) = Tri%StencilCC(2,3,4)
		    GenerateA_CC(k+1+TamCol,k) = Tri%StencilCC(3,2,4)
		    k = k + 1
		    TamCol = TamCol + 2
		end do
		!Finally L3 side
	    !First Xp2
	    GenerateA_CC(k,k) = Tri%StencilCC(2,2,2)
	    GenerateA_CC(k+1,k) = Tri%StencilCC(1,2,2)
	    k = k + 1
	    !General situation, Up triangles
	    aux2 = k
	    k = k + 1
	    do m = 2, TamCol-3,2
	    	GenerateA_CC(k-1,k) = Tri%StencilCC(2,3,6)
	        GenerateA_CC(k,k) = Tri%StencilCC(2,2,6)
	        GenerateA_CC(k+1,k) = Tri%StencilCC(1,2,6)
	        k = k + 2
	    end do
	    !General situation, down triangles
	    k = aux2
	    do m = 1, TamCol-2,2
	        GenerateA_CC(k-1,k) = tri%StencilCC(3,2,7)
	        GenerateA_CC(k,k) = Tri%StencilCC(2,2,7)
	        GenerateA_CC(k+1,k) = Tri%StencilCC(2,3,7)
	        GenerateA_CC(k+1-TamCol,k) = Tri%StencilCC(1,2,7)
	        k = k + 2
	    end do
	    k = k - 1
		!Finally Xp3
	    GenerateA_CC(k,k) = Tri%StencilCC(2,2,3)
	    GenerateA_CC(k-1,k) = Tri%StencilCC(2,3,3)
	end function GenerateA_CC
	!Creates the vector b for the equation system asociated to our stencil
	!Tri(in) :: Type(Triangle)
	!Returns :: double (:)
	!The numeration follows the L3 side. Starting with Xp1 as 1. As an example:
	!            /\
	!		    /  \
	!		   /_4__\
	!		  /\    /\  L3
	!		 /  \ 3/  \
	!	Xp1	/_1__\/_2__\
	function GenerateB_CC(Tri)
	    Implicit none
		!Global variables
		type(Triangle), intent(in) :: Tri
		!Local variables
		double precision, dimension(Tri%SizeUp**2) :: GenerateB_CC
		integer :: j,i,k, a
		k = 1
		!prepare data
		a = size(Tri%UpUF,3)
		!Up triangles
		do i = 2, Tri%SizeUp+1
		    do j = 1, i-1
		        GenerateB_CC(k) = Tri%UpUF(i,j,a)
		        k = k + 2
		    end do
		    k = k - 1
		end do
		!Down triangles
		k = 3
		do i = 2, Tri%SizeDown+1
		    do j = 1, i-1
		        GenerateB_CC(k) = Tri%DownUF(i,j,a)
		        k = k + 2
		    end do
		    k = k + 1
		end do
	end function GenerateB_CC
	!Returns the equivalent matrix to the stencil we use.
	!Xp1 (in) :: double(:). Vertex of the triangle
	!Xp2 (in) :: double(:). Vertex of the triangle
	!Xp3 (in) :: double(:). Vertex of the triangle
	!Trg (in) :: integer. Number of triangles
	!Return :: double(:,:)
	!The numeration follows the L2 side. Starting with Xp1 as 1. The nodes of the L1 and L3 side,
	!are one row, and the L2 nodes are the next row
	!            /\
	!		    /  \
	!		   /_3__\
	!	L1	  /\    /\   L3
	!		 /  1  2  \
	!	Xp1	/____\/____\
	!			L2
!	function GenerateA_SD(Xp1, Xp2, Xp3, Trg)
!	    Implicit none
!	    !Global variables
!	    double precision, intent(in), dimension(2) :: Xp1, Xp2, Xp3
!	    Integer, intent(in) :: Trg
!	    !Local variables
!	    double precision, allocatable, dimension(:,:) :: GenerateA_SD
!	    type (Triangle) :: Tri
!	    integer :: i,j,m
!
!		!Prepare Triangle structure
!		call ConstructorSD(Tri, Xp1, Xp2, Xp3, Trg)
!		Tri%UpUF = 0d0
!		allocate(GenerateA_SD(((Tri%SizeUp/2)*(Tri%SizeUp+1)-Tri%SizeUp)*3, ((Tri%SizeUp/2)*(Tri%SizeUp+1)-Tri%SizeUp)*3))
!		!***L2 manually***
!		!Xp1
!		m = 1
!		Tri%UpUF(1,2,3) = 1d0
!		call getResidualSD(Tri)
!		GenerateA_SD(m,:) = Triangle2vector(Tri)
!		m = m + 1
!		Tri%UpUF=0d0
!		!Interior points
!		do i = 3, Tri%SizeUp
!			Tri%UpUF(1,i,1) = 1d0
!			call getResidualSD(Tri)
!			GenerateA_SD(m,:) = Triangle2vector(Tri)
!			m = m + 1
!			Tri%UpUF=0d0
!
!			Tri%UpUF(1,i,3) = 1d0
!			call getResidualSD(Tri)
!			GenerateA_SD(m,:) = Triangle2vector(Tri)
!			m = m + 1
!			Tri%UpUF=0d0
!		end do
!		!Xp2
!		Tri%UpUF(1,Tri%SizeUp+1,1) = 1d0
!		call getResidualSD(Tri)
!		GenerateA_SD(m,:) = Triangle2vector(Tri)
!		m = m + 1
!		Tri%UpUF=0d0
!
!		!***General situation***
!		do j = 2, Tri%SizeUp-1
!		    !**First L2 type nodes**
!		    do i = j+1, Tri%SizeUp+1
!				Tri%UpUF(j,i,2) = 1d0
!				call getResidualSD(Tri)
!				GenerateA_SD(m,:) = Triangle2vector(Tri)
!				m = m + 1
!				Tri%UpUF=0d0
!		    end do
!
!		    !**Other nodes**
!		    !L1 side
!			Tri%UpUF(j,j+1,3) = 1d0
!			call getResidualSD(Tri)
!			GenerateA_SD(m,:) = Triangle2vector(Tri)
!			m = m + 1
!			Tri%UpUF=0d0
!
!		    do i = j+2, Tri%SizeUp
!				Tri%UpUF(j,i,1) = 1d0
!				call getResidualSD(Tri)
!				GenerateA_SD(m,:) = Triangle2vector(Tri)
!				m = m + 1
!				Tri%UpUF=0d0
!
!				Tri%UpUF(j,i,3) = 1d0
!				call getResidualSD(Tri)
!				GenerateA_SD(m,:) = Triangle2vector(Tri)
!				m = m + 1
!				Tri%UpUF=0d0
!		    end do
!		    !L3 side
!			Tri%UpUF(j,Tri%SizeUp+1,1) = 1d0
!			call getResidualSD(Tri)
!			GenerateA_SD(m,:) = Triangle2vector(Tri)
!			m = m + 1
!			Tri%UpUF=0d0
!		end do
!		!***Xp3 manually***
!		Tri%UpUF(Tri%SizeUp,Tri%SizeUp+1,2) = 1d0
!		call getResidualSD(Tri)
!		GenerateA_SD(m,:) = Triangle2vector(Tri)
!		!Finally change the sign
!		GenerateA_SD = - transpose(GenerateA_SD)
!	end function GenerateA_SD
	!Creates the vector b for the equation system asociated to our stencil
	!Tri(in) :: Type(Triangle)
	!Returns :: double (:)
	!The numeration follows the L2 side. Starting with Xp1 as 1. The nodes of the L1 and L3 side,
	!are one row, and the L2 nodes are the next row
	!            /\
	!		    /  \
	!		   /_3__\
	!	L1	  /\    /\   L3
	!		 /  1  2  \
	!	Xp1	/____\/____\
	!			L2
	function GenerateB_SD(Tri)
	    Implicit none
		!Global variables
		type(Triangle), intent(in) :: Tri
		!Local variables
		double precision, dimension(((Tri%SizeUp/2)*(Tri%SizeUp+1)-Tri%SizeUp)*3) :: GenerateB_SD
		integer :: j,i,k, TamRow
		k = 1
		!First L3 nodes
		TamRow = Tri%SizeUp-1!Size of the elements of the L2 kind
		do j = 1, Tri%SizeUp-1
		    do i = j+1, Tri%SizeUp
				GenerateB_SD(k) = Tri%UpUF(i,j,3)
				k = k + 2
		    end do
		    k = k + TamRow
		    TamRow = TamRow - 1
		end do
		k = 2
		!Secondly L1 nodes
		TamRow = Tri%SizeUp-1!Size of the elements of the L2 kind
		do j = 1, Tri%SizeUp-1
		    do i = j+2, Tri%SizeUp+1
				GenerateB_SD(k) = Tri%UpUF(i,j,1)
				k = k + 2
		    end do
		    k = k + TamRow
		    TamRow = TamRow - 1
		end do
		k = Tri%SizeUp * 2 - 1
		!Finally L2 nodes
		TamRow = Tri%SizeUp * 2 - 4
		do j = 2, Tri%SizeUp
		    do i = j+1, Tri%SizeUp+1
				GenerateB_SD(k) = Tri%UpUF(i,j,2)
				k = k + 1
		    end do
		    k = k + TamRow
		    TamRow = TamRow - 2
		end do
	end function GenerateB_SD
	!Function to use only in GenerateA_SD
	!Returns the Residual of the input Type(Triangle) as a vector
	function Triangle2vector(Tri)
		Implicit none
		!Global variables
		type (Triangle), intent(in) :: Tri
		!Local variables
		double precision, allocatable, dimension(:) :: Triangle2vector
		integer :: i, j, m

		allocate(Triangle2vector(((Tri%SizeUp/2)*(Tri%SizeUp+1)-Tri%SizeUp)*3))

		!***L2 manually***
		!Xp1
		m=1
		Triangle2vector(m) = Tri%UpResAux(1,2,3)
		m = m + 1
		!Interior points
		do i = 3, Tri%SizeUp
		    Triangle2vector(m) = Tri%UpResAux(1,i,1)
		    m = m + 1
		    Triangle2vector(m) = Tri%UpResAux(1,i,3)
		    m = m + 1
		end do
		!Xp2
		Triangle2vector(m) = Tri%UpResAux(1,Tri%SizeUp+1,1)
		m = m + 1
		!***General situation***
		do j = 2, Tri%SizeUp-1
		    !**First L2 type nodes**
		    do i = j+1, Tri%SizeUp+1
				Triangle2vector(m) = Tri%UpResAux(j,i,2)
				m = m + 1
		    end do
		    !**Other nodes**
		    !L1 side
		    Triangle2vector(m) = Tri%UpResAux(j,j+1,3)
		    m = m + 1
		    do i = j+2, Tri%SizeUp
			    Triangle2vector(m) = Tri%UpResAux(j,i,1)
		    	m = m + 1
			    Triangle2vector(m) = Tri%UpResAux(j,i,3)
		    	m = m + 1
		    end do
		    !L3 side
	 	    Triangle2vector(m) = Tri%UpResAux(j,Tri%SizeUp+1,1)
	    	m = m + 1
		end do
		!***Xp3 manually***
		Triangle2vector(m) = Tri%UpResAux(Tri%SizeUp,Tri%SizeUp+1,2)

	end function
end module LinearSolvers
