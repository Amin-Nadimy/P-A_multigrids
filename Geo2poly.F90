module Geo2poly
	use poly_io
	use structures
	use strings
	use evaluate
	use generic
    implicit none

    type, public :: geometry
        double precision, allocatable, dimension(:,:) :: Point
		integer, allocatable, dimension(:,:) :: Line
		integer, allocatable, dimension(:,:) :: LineLoop
		integer, allocatable, dimension(:,:) :: Surface
    end type geometry

contains
	!Deallocates all the geo components
	!geo(inout) :: type(geometry)
	subroutine Geo_destructor(geo)
	    implicit none
	    !Global variables
	    type (geometry), intent(inout) :: geo

	    deallocate(geo%point,geo%line,geo%LineLoop,geo%Surface)

	end subroutine Geo_destructor

	!Reads the 2D geometrical data from a .geo file (but doesn't work if there are circles)
	!geo(inout) :: type(PSGL)
	!Io(out) :: integer. If /= 0 then a problem happend while reading the file
	!file(in) :: string with the relative path to the .geo file.
	!NOTE: ALL THE ARRAYS USED HERE HAVE A FIXED SIZE, THIS MAY BE A PROBLEM IN THE FUTURE
	subroutine Read_geo(geo, file, Io)
	    implicit none
	    !Global variables
	    type (geometry), intent(inout) :: geo
	    character (len=*), intent(in) :: file
	    integer, intent(inout) :: Io
	    !Local variables
	    character (len=1000) :: cadena, filex, path
	    character(len=100),dimension(100) :: parsed
	    character (len=100) :: delim
	    integer :: aux, point, line, lineloop, plane, i, k, p
	    integer, dimension(500) :: TableConversion
	    double precision :: auxR
	    !Prepare data
	    delim = "()=;{},"
		filex = trim(file)
		TableConversion = 0
		!Allocate (stimate values)
		allocate(geo%Point(2,1000),geo%Line(2,500), geo%LineLoop(100,100),geo%Surface(50,50))
		!Just to identify the numbers I have introduce
		geo%Point = huge(auxR)
		geo%Line = huge(aux)
		geo%LineLoop = huge(aux)
		geo%Surface = huge(aux)
		!Remove the .geo if exists
		call delall(filex, ".geo")
		!Add the .geo
		call insertstr(filex,".geo",len_trim(filex)+1)
		!Open file
		call getCWD (path)
		!To avoid using debug as working directory
		call delsubstr(path,"/Debug")
		path = trim(path)//"/"//trim(filex)
		open(1,file = trim(path), status="old", action = "read")
		point = 1
		line = 1
		lineloop = 1
		plane = 1
		!Read a line
		call readline(1,cadena,io)
		do while (Io == 0)
			!To lowercase
			cadena = lowercase(cadena)
			!Remove spaces
			call removesp(cadena)
			!Parse a line
			call parse(cadena, delim, parsed, aux)
			!Analize
			!Four posibilities
			!Obtain the points
			if (trim(parsed(1))=="point") then
				i = 3
				do while (len_trim(parsed(i))==0)
				    i = i + 1
				end do
				!Get (x)
				call value(trim(parsed(i)),auxR,Io)
				geo%Point(1,point) = auxR
				!Get (y)
				call value(trim(parsed(i+1)),auxR,Io)
				geo%Point(2,point) = auxR
				point = point + 1
			end if
			!Connect those points, to create lines
			if (trim(parsed(1))=="line") then
				i = 3
				do while (len_trim(parsed(i))==0)
				    i = i + 1
				end do
				!Get first point
				call value(trim(parsed(i)),aux,Io)
				geo%Line(1,line) = aux
				!Get last point
				call value(trim(parsed(i+1)),aux,Io)
				geo%Line(2,line) = aux
				line = line + 1
			end if
			!Connect lines, to create "surfaces"
			if (trim(parsed(1))=="lineloop") then
				!Get the position
				call value(trim(parsed(2)),aux,Io)
				TableConversion(lineloop) = aux
				i = 3
				do while (len_trim(parsed(i))==0)
				    i = i + 1
				end do
				!get a "line"
				call value(trim(parsed(i)),aux,Io)
				k = 1
				do
					!Save the "line"
				    geo%LineLoop(lineloop,k) = abs(aux)
				    i = i + 1
				    k = k + 1
				    !Until loop
					if (len_trim(parsed(i))==0) exit
					call value(parsed(i),aux,Io)
				end do
				lineloop = lineloop + 1
			end if
			!Create surfaces, taking into account holes
			if (trim(parsed(1))=="planesurface") then
				i = 3
				do while (len_trim(parsed(i))==0)
				    i = i + 1
				end do
				!get a "lineloop"
				call value(trim(parsed(i)),aux,Io)
				k = 1
				do
					p = 1
					!Use tableConversion
					do while (TableConversion(p)/=aux)
					    p = p + 1
					end do

					!Save the "Surface"
				    geo%Surface(plane,k) = p
				    i = i + 1
				    k = k + 1
				    if (len_trim(parsed(i))==0) exit
					call value(trim(parsed(i)),aux,Io)
				end do
				plane = plane + 1
				!Restart Io
				Io = 0
			end if

			!Read a line
			call readline(1,cadena,io)
		end do

	end subroutine Read_geo



	!Creates a poly archive to be used by acute, this file will be name as input.poly
	!geo(inout) :: type(geometry)
	subroutine CreatePoly(geo)
	    Implicit none
	    !Global variables
	    type (geometry), intent(inout) :: geo
	    !Local
	    character(len=500) :: path
		integer :: points, lines, holes, i, aux, k, aux2, j
		double precision :: auxR
		double precision, allocatable, dimension(:,:) :: hole_point
		!Number of points
		auxR = huge(auxR)
		aux = huge(aux)
		!Calculate number of points
		points = 1
		do while (.not.AreEqual(geo%point(1,points),auxR))
		    points = points + 1
		end do
		points = points - 1
		!Calculate number of Lines
		lines = 1
		do while (geo%Line(1,lines)/= aux)
		    lines = lines + 1
		end do
		lines = lines - 1
		!Calculate number of holes
		i = 1
		holes = 0
		do while (geo%Surface(i,1)/= aux)
			k = 2
		    do while (geo%Surface(i,k)/= aux)
		        k = k + 1
		        holes = holes + 1
		    end do
		    i = i + 1
		end do

		!Hole points
		allocate(hole_point(2,holes))

		!Look for holes
		j = 1
		i = 1
		do while (geo%Surface(i,1)/= aux)
			k = 2
		    do while (geo%Surface(i,k)/= aux)
		        !Hole detected, look for a point
		        !Firstly find the surface
		        aux2 = geo%LineLoop(geo%Surface(i,k),1)
		        !Get a point
		        aux2 = geo%Line(1,aux2)
		        !Store that point
		        hole_point(:,j) = geo%Point(:,aux2)
		        j = j + 1
		        k = k + 1
		    end do
		    i = i + 1
		end do
		!Open file
		call getCWD (path)
		!To avoid using debug as working directory
		call delsubstr(path,"/Debug")
		path = trim(path)//"/input.poly"
		call poly_write ( trim(path), points, geo%point(:,1:points), lines,&
			geo%line(:,1:points), holes, hole_point )

	end subroutine CreatePoly

	!Reads the information generated by the aCute software given a .poly file
	!meshList(out) :: type(Mesh)(:)
	!ios(out) :: integer. If /= 0 then a problem happend while reading the file
	!remove(in) :: logical. If true, then the .ele and the .node will be erased, as well as the .1.poly.
	!(Optional)filename(in) :: string with the relative path to the .poly file.
	!							By DEFAULT: input
	!(Optional)subiteration(in) :: integer. ACute creates different archives deppending on the iteration.
	!							By DEFAULT: 1 or none
	subroutine ReadPoly_NotIntel(meshL,ios, remove,filename, subiteration)
        Implicit none
	    !Global variables
		type(Mesh), intent(out), allocatable, dimension(:) :: meshL
		character (len = *), optional, intent(in) :: filename
		integer, intent(out) :: ios
		integer,optional, intent(in) :: subiteration
		logical, intent(in) :: remove
		!Local variables
		character(len=500) :: ele, cadena, node, filex, path
		character(len=2) :: delim
		character(len=100),dimension(30) :: parsed
		integer :: i, aux
		logical :: bol
		double precision :: daux
		double precision, allocatable, dimension(:,:) :: vertex !(2,:)(X,Y)
		!Prepare data
		delim = " "
		filex = "input"!.poly"

		if (present(filename)) then
		    filex = filename
		    !Remove the .poly
		    call delall(filex, ".poly")
		end if
		!Adapt to parallel
		call getCWD (path)
		!To avoid using debug as working directory
		call delsubstr(path,"/Debug")
		filex = trim(path)//"/"//trim(filex)

		if(present(subiteration))then
			call writenum(subiteration,cadena,"i5")
			!Create .ele link
			ele = trim(filex)//trim(cadena)//".ele"
			!Create .node link
			node = trim(filex)//trim(cadena)//".node"
		else
			!Create .ele link
			ele = trim(filex)//".1.ele"
			inquire( file=ele, exist=bol)
			if (.not.bol) then
			    ele = trim(filex)//".ele"
			end if

			!Create .node link
			node = trim(filex)//".1.node"
			inquire( file=node, exist=bol)
			if (.not.bol) then
			    node = trim(filex)//".node"
			end if

		end if
		!****At first i will store the coordinates****
		open(1,file = node, status="old", action = "read")

		call readline(1,cadena,ios)
		!Get the number of vertex
		call parse(cadena, delim, parsed, aux)
		!I am only interested in the first number, wich is the total number of vertex
		call value(parsed(1),aux,ios)
		allocate(vertex(2,aux))

		!Store all the coordinates
		do i = 1, aux
			!Now read line by line including the data
			call readline(1,cadena,ios)
			!Split cadena
			call parse(cadena, delim, parsed, aux)
			!We are only interested in the second number, wich is the X
			call value(parsed(2),daux,ios)
			vertex(1,i) = daux
			!And in the third, wich is the Y
			call value(parsed(3),daux,ios)
			vertex(2,i) = daux
		end do
		close(1)

		!****By reading the .ele file we will create the meshList****
		open(1,file = ele, status="old", action = "read")

		call readline(1,cadena,ios)
		!Get the number of triangles
		call parse(cadena, delim, parsed, aux)
		!I am only interested in the first number, wich is the total number of triangles
		call value(parsed(1),aux,ios)
		allocate(meshL(aux))
		!Store the triangles
		do i = 1, aux
			!Now read line by line including the data
			call readline(1,cadena,ios)
			!Get split cadena
			call parse(cadena, delim, parsed, aux)
			!We are interested in the second number
			call value(parsed(2),aux,ios)
			meshL(i)%Xp1=vertex(:,aux)
			!The third
			call value(parsed(3),aux,ios)
			meshL(i)%Xp2=vertex(:,aux)
			!And the forth
			call value(parsed(4),aux,ios)
			meshL(i)%Xp3=vertex(:,aux)
		end do
		close(1)

		!****Look for neighbours****
		!Despite we can use aCute to calculate the neighbours, and it is probably faster,
		!it will not tell us the direction between two connected triangles
		do i = 1, size(meshL,1)
		    do aux = i + 1, size(meshL,1)
		        call CheckNeig_NotIntel(meshL,i,aux)
		    end do
		end do

	end subroutine ReadPoly_NotIntel
!
	!Returns true if the input triangles of the mesh list are neighbours
	!ml(inout) :: type(Mesh)(:)
	!i(in) :: integer
	!j(in) :: integer
	!Return: logical
	subroutine CheckNeig_NotIntel(ml, i, j)
		Implicit none
		!Global variables
		type(Mesh), intent(inout), dimension(:) :: ml
		integer, intent(in) :: j,i
		!Local variables
		integer, dimension(4) :: Vertex
		logical :: oneMatch, twoMatch, threeMatch, oneMatch2, twoMatch2, threeMatch2
		!Initialize variables
		oneMatch = .false.
		twomatch = .false.
		threematch = .false.
		oneMatch2 = .false.
		twomatch2 = .false.
		threematch2 = .false.

		Vertex = 0
		!Vertex goes from Xp1 to Xp3 of both triangles
		!The forth is just to check that 2 vertex coincides
		!Which means that they do not have a pair

		!Check Xp1 of i triangle
		if (AreEqual(ml(i)%Xp1,ml(j)%Xp1)) then
			vertex(1) = 1
			oneMatch = .true.
			oneMatch2 = .true.
		else if (AreEqual(ml(i)%Xp1,ml(j)%Xp2)) then
			vertex(1) = 2
			oneMatch = .true.
			twoMatch2 = .true.
		else if (AreEqual(ml(i)%Xp1,ml(j)%Xp3)) then
			vertex(1) = 3
			oneMatch = .true.
			threeMatch2 = .true.
		end if

		!Check Xp2 of i triangle
		if (AreEqual(ml(i)%Xp2,ml(j)%Xp1)) then
			vertex(2) = 2
			twomatch = .true.
			oneMatch2 = .true.
		else if (AreEqual(ml(i)%Xp2,ml(j)%Xp2)) then
			vertex(2) = 5
			twomatch = .true.
			twoMatch2 = .true.
		else if (AreEqual(ml(i)%Xp2,ml(j)%Xp3)) then
			vertex(2) = 6
			twomatch = .true.
			threeMatch2 = .true.
		end if

		!Check Xp3 of i triangle
		if (AreEqual(ml(i)%Xp3,ml(j)%Xp1)) then
			vertex(3) = 3
			threematch = .true.
			oneMatch2 = .true.
		else if (AreEqual(ml(i)%Xp3,ml(j)%Xp2)) then
			vertex(3) = 6
			threematch = .true.
			twoMatch2 = .true.
		else if (AreEqual(ml(i)%Xp3,ml(j)%Xp3)) then
			vertex(3) = 9
			threematch = .true.
			threeMatch2 = .true.
		end if

		!Store neighbour
		if (oneMatch .and. threeMatch) then
			ml(i)%Neig(1)=j
			if (checkvector(vertex,1) .or.( CheckVector(vertex,2) .and. CheckVector(vertex,9)) ) then
			    ml(i)%Dir(1) = .true.
			end if
			threeMatch = .false.
	 	else if (oneMatch .and. twoMatch) then
			ml(i)%Neig(2)=j
			if (checkvector(vertex,1) .or.( CheckVector(vertex,6) .and. CheckVector(vertex,2)) ) then
			    ml(i)%Dir(2) = .true.
			end if
			onematch = .false.
		else if (threeMatch .and. twoMatch) then
			ml(i)%Neig(3)=j

			if (checkvector(vertex,9) .or.( CheckVector(vertex,6) .and. CheckVector(vertex,2)) ) then
			    ml(i)%Dir(3) = .true.
			end if
			twomatch = .false.
		end if

		!Store neighbour in neighbour
		if (oneMatch2 .and. threeMatch2) then
			ml(j)%Neig(1)=i
		 	if (OneMatch) then
		 		ml(j)%Dir(1) = ml(i)%Dir(1)
	 	    else if (TwoMatch) then
	 	    	ml(j)%Dir(1) = ml(i)%Dir(2)
 	        else if (ThreeMatch) then
 	        	ml(j)%Dir(1) = ml(i)%Dir(3)
		 	end if

	 	else if (oneMatch2 .and. twoMatch2) then
			ml(j)%Neig(2)=i
		 	if (OneMatch) then
		 		ml(j)%Dir(2) = ml(i)%Dir(1)
	 	    else if (TwoMatch) then
	 	    	ml(j)%Dir(2) = ml(i)%Dir(2)
 	        else if (ThreeMatch) then
 	        	ml(j)%Dir(2) = ml(i)%Dir(3)
		 	end if
		else if (threeMatch2 .and. twoMatch2) then
			ml(j)%Neig(3)=i
		 	if (OneMatch) then
		 		ml(j)%Dir(3) = ml(i)%Dir(1)
	 	    else if (TwoMatch) then
	 	    	ml(j)%Dir(3) = ml(i)%Dir(2)
 	        else if (ThreeMatch) then
 	        	ml(j)%Dir(3) = ml(i)%Dir(3)
		 	end if
		end if

	end subroutine

end module Geo2poly
