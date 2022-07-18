module Generic
    implicit none


	private
	public :: AreEqual, CheckVector, NumLoc, QuickSort, getTrapezoidArea, getTriangleArea
	public :: checkPointInsideSquare, vector_by_Matrix, Matrix_by_vector, Matrix_by_Matrix
	public :: getProcessor
	!Check whether two double are equal or not
	!
	!Subroutine 1:
	!X(in) :: double/real
	!Y(in) :: double/real
	!returns true if they are equal
	!
	!Subroutine 2:
	!X(in) :: double(:)
	!Y(in) :: double(:)
	!returns true if they are equal
    interface AreEqual
        module procedure AreEqual1
        module procedure AreEqual2
        module procedure AreEqual3
    end interface AreEqual

!This module includes all the generic subroutines or functiones needed by other algorithms
contains


	!Check whether two double are equal or not
	!X(in) :: double
	!Y(in) :: double
	!returns true if they are equal
	logical function AreEqual1(X, Y)
	    Implicit none
	    !Global variables
	    double precision, intent(in) :: X, Y
	    !Local variables
	    double precision :: eps
	    eps = epsilon(eps)

	    AreEqual1 = abs(X-Y) < eps !* max(abs(X),abs(Y))

	end function AreEqual1

	!Check whether two double are equal or not
	!X(in) :: double(:)
	!Y(in) :: double(:)
	!returns true if they are equal
	logical function AreEqual2(X, Y)
	    Implicit none
	    !Global variables
	    double precision, intent(in), dimension(:) :: X, Y
	    !Local variables
	    double precision :: eps
	    eps = epsilon(eps)
	    AreEqual2 = sqrt(dot_product(X-Y, X-Y)) < eps !* max(maxval(abs(X)),maxval(abs(Y)))
	end function AreEqual2

	!Check whether two double are equal or not
	!X(in) :: real
	!Y(in) :: real
	!returns true if they are equal
	logical function AreEqual3(X, Y)
	    Implicit none
	    !Global variables
	    real, intent(in) :: X, Y
	    !Local variables
	    real :: eps
	    eps = epsilon(eps)

	    AreEqual3 = abs(X-Y) < eps !* max(abs(X),abs(Y))

	end function AreEqual3
	!Checks if the vector contains the number
	!Vector(in) :: integer(:)
	!Num(in) :: integer
	!Return: true if founded
	logical function CheckVector(Vector, num)
	    Implicit none
	    !Global variables
	    integer, intent(in), dimension(:) :: Vector
	    integer, intent(in) :: Num
	    !Local variables
	    integer :: i
	    CheckVector = .false.
!	    if (any(vector==Num)) CheckVector = .true.
	    do i = 1, size(vector)
	        if (vector(i) == Num) then
	        	CheckVector = .true.
	        	return
	        end if
	    end do
	end function CheckVector

	!Multiplication of a matrix by a vector
	!vector(in) :: double(:)
	!mat(in) :: double(:,:)
	!Return: double(size(mat,1))
	!Note, the row size of mat must be the same that the size of vector
	function vector_by_Matrix(vector, mat)
	    Implicit none
	    !Global variables
	    double precision, intent(in), dimension(:) :: Vector
	    double precision, intent(in), dimension(:,:) :: mat
	    !Local variables
	    double precision, dimension(size(mat,1)) :: vector_by_Matrix
	    integer :: i

	    do i = 1, size(mat,1)
	        vector_by_Matrix(i) = dot_product(Vector,mat(i,:))
	    end do

	end function vector_by_Matrix

	!Multiplication of two vectors
	!vect1(in) :: double(:)
	!vect2(in) :: double(:)
	!Return: double
	doubleprecision function Vector_by_vector(vect1,vect2)
	    Implicit none
	    !Global variables
	    double precision, intent(in), dimension(:) :: Vect1, vect2
		!Local variables
		integer :: i
		Vector_by_vector = 0d0
		do i = 1, size(vect1,1)
		    Vector_by_vector = Vector_by_vector + vect1(i)*vect2(i)
		end do

	end function Vector_by_vector

	!Multiplication of a matrix by a vector
	!vector(in) :: double(:)
	!mat(in) :: double(:,:)
	!Return: double(size(mat,1))
	!Note, the row size of mat must be the same that the size of vector
	function Matrix_by_vector(mat,vector)
	    Implicit none
	    !Global variables
	    double precision, intent(in), dimension(:) :: Vector
	    double precision, intent(in), dimension(:,:) :: mat
	    !Local variables
	    double precision, dimension(size(mat,1)) :: Matrix_by_vector
	    integer :: i

	    do i = 1, size(mat,1)
	        Matrix_by_vector(i) = dot_product(mat(:,i),vector)
	    end do

	end function Matrix_by_vector

	!Multiplication of a matrix by a vector
	!mat1(in) :: double(:)
	!mat2(in) :: double(:,:)
	!Returns: Matrix_by_Matrix(j,i) = dot_product(mat1(:,i),mat2(j,:))
	function Matrix_by_Matrix(mat1,mat2)
	    Implicit none
	    !Global variables
	    double precision, intent(in), dimension(:,:) :: mat1, mat2
	    !Local variables
	    double precision, dimension(size(mat2,2),size(mat1,1)) :: Matrix_by_Matrix
	    integer :: i, j

	    do i = 1, size(mat1,1)
		    do j = 1, size(mat2,2)
		        Matrix_by_Matrix(j,i) = dot_product(mat1(:,i),mat2(j,:))

		    end do
	    end do

	end function Matrix_by_Matrix

	!Returns the position of a number inside the vector.
	!The seek goes from 1 to size(vector). If not found returns 0
	!Vector(in) :: Integer(:)
	!Number(in) :: Integer
	!Return: Position as an integer
	integer function NumLoc(Vector, Number)
	    Implicit none
	    !Global variables
	    integer, intent(in), dimension(:) :: vector
	    Integer, intent(in) :: number
	    !Local variables
	    integer :: i

	    do NumLoc = 1, size(vector)
	        if (Number == Vector(NumLoc)) return
	    end do
	    	NumLoc = 0
	end function NumLoc
	!Returns true if the poitn Xp falls inside the area
	!delimited by X1, X2, X3 and X4
	!X1(in) :: double(2)
	!X2(in) :: double(2)
	!X3(in) :: double(2)
	!X4(in) :: double(2)
	!Xp(in) :: double(2)
	!Return :: logical
	logical function checkPointInsideSquare(X1,X2,X3,X4,Xp)
	    Implicit none
	    !Global variables
	    double precision, intent(in), dimension(2) :: X1,X2,X3,X4, Xp
	    !Local variables
	    double precision :: A1,A2,A3,A4,Area, aux3, aux4
	    double precision, dimension(4) :: aux
	    double precision, dimension(2) :: aux2, Xp1,Xp2,Xp3,Xp4
	    logical :: ABC, ACD, ABD
    	!Prepare data
    	checkPointInsideSquare = .false.
    	Xp1 = X1
    	Xp2 = X2
    	Xp3 = X3
    	Xp4 = X4
!    	!Check clockwise/anticlokwise ordering
!		aux(1) = dot_product(Xp1,Xp1)
!		aux(2) = dot_product(Xp2,Xp2)
!		aux(3) = dot_product(Xp3,Xp3)
!		aux(4) = dot_product(Xp4,Xp4)
!		!Sort the list
!		call QuickSort(aux)
!		!Store the highest in Xp4
!		if (AreEqual(aux(4),dot_product(Xp3,Xp3))) then
!			aux2 = Xp4
!			Xp4 = Xp3
!			Xp3 = aux2
!		end if
!		if (AreEqual(aux(4),dot_product(Xp2,Xp2))) then
!			aux2 = Xp4
!			Xp4 = Xp2
!			Xp2 = aux2
!		end if
!		if (AreEqual(aux(4),dot_product(Xp1,Xp1))) then
!			aux2 = Xp4
!			Xp4 = Xp1
!			Xp1 = aux2
!		end if
!		!And the lowest in Xp1
!		if (AreEqual(aux(1),dot_product(Xp3,Xp3))) then
!			aux2 = Xp1
!			Xp1 = Xp3
!			Xp3 = aux2
!		end if
!		if (AreEqual(aux(1),dot_product(Xp2,Xp2))) then
!			aux2 = Xp1
!			Xp1 = Xp2
!			Xp2 = aux2
!		end if
!		if (AreEqual(aux(1),dot_product(Xp4,Xp4))) then
!			aux2 = Xp1
!			Xp1 = Xp4
!			Xp4 = aux2
!		end if
!		!Sort coordinates
!		!First cirtual triangle:
!		if (getTriangleDeterminant(Xp1,Xp2,Xp4)>0) then
!		    if (getTriangleDeterminant(Xp1,Xp2,Xp4)>0) then
!				aux2 = Xp3
!				Xp3 = Xp4
!				Xp4 = aux2
!		    end if
!		end if

!		!To order Xp3 and Xp2 we will check a point we are sure is inside, with an initial
!		!configuration. If it works it is well ordered. If not, we swap Xp3 with Xp2
!		 aux3 = abs(getTriangleDeterminant(Xp1,Xp2,Xp1)+getTriangleDeterminant(Xp2,Xp4,Xp1)+&
!		 	getTriangleDeterminant(Xp4,Xp3,Xp1) + getTriangleDeterminant(Xp3,Xp1,Xp1))
!
!		 aux4 = abs(getTriangleDeterminant(Xp1,Xp2,Xp1))+abs(getTriangleDeterminant(Xp2,Xp4,Xp1))+&
!		 	abs(getTriangleDeterminant(Xp4,Xp3,Xp1)) + abs(getTriangleDeterminant(Xp3,Xp1,Xp1))
!
!		!We swap values
!		if (.not.AreEqual(aux3,aux4)) then
!			aux2 = Xp2
!			Xp2 = Xp3
!			Xp3 = aux2
!		end if

		!Order clockwise/anti-clockwise
		ABC = getTriangleDeterminant(Xp1,Xp2,Xp4) > 0d0
		ACD = getTriangleDeterminant(Xp1,Xp4,Xp3) > 0d0
		ABD = getTriangleDeterminant(Xp1,Xp2,Xp3) > 0d0

		if (.not.ABC.and.ACD.and..not.ABD) then
			aux2 = Xp3
			Xp3 = Xp4
			Xp4 = aux2
		end if
		if (ABC.and..not.ACD.and..not.ABD) then
			aux2 = Xp2
			Xp2 = Xp4
			Xp4 = aux2
		end if
		if (ABC.and..not.ACD.and.ABD) then
			aux2 = Xp3
			Xp3 = Xp4
			Xp4 = aux2
		end if
		if (.not.ABC.and.ACD.and.ABD) then
			aux2 = Xp2
			Xp2 = Xp4
			Xp4 = aux2
		end if
		!Now that it is indeed clockwise ordered we check the input point
		 aux3 = abs(getTriangleDeterminant(Xp1,Xp2,Xp)+getTriangleDeterminant(Xp2,Xp4,Xp)+&
		 	getTriangleDeterminant(Xp4,Xp3,Xp) + getTriangleDeterminant(Xp3,Xp1,Xp))

		 aux4 = abs(getTriangleDeterminant(Xp1,Xp2,Xp))+abs(getTriangleDeterminant(Xp2,Xp4,Xp))+&
		 	abs(getTriangleDeterminant(Xp4,Xp3,Xp)) + abs(getTriangleDeterminant(Xp3,Xp1,Xp))
		if (AreEqual(aux3,aux4)) checkPointInsideSquare = .true.

	end function checkPointInsideSquare

	!Returns the determinant of the triangle, only to be used with checkPointInsideSquare
	real function getTriangleDeterminant(X1,X2,X3)
	    Implicit none
	    !Global variables
	    double precision, intent(in), dimension(2) :: X1,X2,X3

	   	getTriangleDeterminant = (X2(1)*X3(2) - X3(1)*X2(2)) - &
	   		(X1(1)*X3(2) - X3(1)*X1(2)) + (X1(1)*X2(2)-X2(1)*X1(2))

	end function getTriangleDeterminant
	!Returns the area of the triangle
	!X1(in) :: double(2)
	!X2(in) :: double(2)
	!X3(in) :: double(2)
	function getTriangleArea(X1,X2,X3)
	    Implicit none
	    !Global variables
	    double precision, intent(in), dimension(2) :: X1,X2,X3
    	!Local variables
    	double precision :: getTriangleArea
    	double precision, dimension(3) :: l

        l(1) = sqrt(dot_product(X1-X3, X1-X3))
        l(2) = sqrt(dot_product(X1-X2, X1-X2))
        l(3) = sqrt(dot_product(X2-X3, X2-X3))

    	getTriangleArea = 0.25d0*sqrt( (l(1)+l(2)-l(3))*(l(1)-l(2)+l(3))*(-l(1)+l(2)+l(3))*(l(1)+l(2)+l(3)) )

	end function getTriangleArea
	!Returns the area difined by 4 coordinates
	!X1(in) :: double(2)
	!X2(in) :: double(2)
	!X3(in) :: double(2)
	!X4(in) :: double(2)
	function getTrapezoidArea(X1,X2,X3,X4)
	    Implicit none
	    !Global variables
	    double precision, intent(in), dimension(2) :: X1,X2,X3,X4
    	!Local variables
    	double precision :: getTrapezoidArea
    	double precision, dimension(4) :: l

    	l(1) = sqrt(dot_product(X1-X2,X1-X2))
    	l(2) = sqrt(dot_product(X1-X3,X1-X3))
    	l(3) = sqrt(dot_product(X4-X2,X4-X2))
    	l(4) = sqrt(dot_product(X4-X3,X4-X3))
		!l(1) will be Lmin and l(4) Lmax
		call QuickSort (l)

		if(AreEqual(l(1),0d0)) then
			l(1) = l(2)
			l(2) = 0d0
		end if

		!Calculate area
		if (l(1)==l(2) .and. l(2)==l(3).and.l(3)==l(4)) then
			!Square
			getTrapezoidArea = l(1)**2
		else if (AreEqual(l(2),0d0)) then
				!Triangle
			   getTrapezoidArea = 0.25d0*sqrt( (l(1)+l(4)-l(3))*(l(1)-l(4)+l(3))*(-l(1)+l(4)+l(3))*(l(1)+l(4)+l(3)) )
		else
			!Trapezoid
			getTrapezoidArea = sqrt( (-l(4)+l(3)+l(1)+l(2))*(l(4)-l(3)+l(1)+l(2))*(l(4)-l(3)+l(1)-l(2))*(l(4)-l(3)-l(1)+l(2)) )
			getTrapezoidArea = getTrapezoidArea * ( l(1) + l(4) )/( 4*(l(4)-l(1)) )
		end if

	end function getTrapezoidArea

	!Returns the processor given a triangle
	!Mpos (in) :: Integer
	!NTriangles (in) :: Integer. Usually size(meshL)
	!Total (in) :: Integer. Number of processors
	integer function getProcessor(Mpos, NTriangles, total)
	    Implicit none
	    !Global variables
	    integer, intent(in) :: Mpos, NTriangles, total
	    !Local variables
	    integer :: aux

		aux = Mpos
		getProcessor = -1
		do while (aux>0)
		    aux = aux - NTriangles/total
		    getProcessor = getProcessor + 1
		end do
		if (getProcessor>=total-1) getProcessor = total-1
	end function getProcessor

!QuickSort algorithm, taken from http://rosettacode.org
!Sorts an array of doubles
!QuickSort(a), where a is double(:)
!Order from less to max
RECURSIVE SUBROUTINE QuickSort(a)

  DOUBLE PRECISION, INTENT(IN OUT) :: a(:)
  DOUBLE PRECISION :: split

  IF(size(a) > 1) THEN
     CALL Partition(a, split)
     CALL QuickSort(a(:split-1))
     CALL QuickSort(a(split:))
  END IF

END SUBROUTINE QuickSort

SUBROUTINE Partition(a, marker)

  DOUBLE PRECISION, INTENT(IN OUT) :: a(:)
  DOUBLE PRECISION, INTENT(OUT) :: marker
  DOUBLE PRECISION :: left, right, pivot, temp

  pivot = (a(1) + a(size(a))) / 2  ! Average of first and last elements to prevent quadratic
  left = 0                         ! behavior with sorted or reverse sorted data
  right = size(a) + 1

  DO WHILE (left < right)
     right = right - 1
     DO WHILE (a(right) > pivot)
        right = right-1
     END DO
     left = left + 1
     DO WHILE (a(left) < pivot)
        left = left + 1
     END DO
     IF (left < right) THEN
        temp = a(left)
        a(left) = a(right)
        a(right) = temp
     END IF
  END DO

  IF (left == right) THEN
     marker = left + 1
  ELSE
     marker = left
  END IF

END SUBROUTINE Partition


end module Generic
