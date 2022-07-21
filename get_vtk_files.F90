module get_vtk_files
  use structured_meshgen
  implicit none

contains

  !@Brief> This subroutine generates VTU files based on XML format for unstructured grids.
  ! an example is used from https://people.math.sc.edu/Burkardt/data/vtu/vtu.html
  ! overall format comes from http://www.princeton.edu/~efeibush/viscourse/vtk.pdf
  subroutine get_vtu(x_all, tnew, error,analytical, totele_str, totele_unst, totele, nloc, itime, ndim, cell_type, solve_for)
    implicit none
    ! Global variables
    integer, intent(in):: itime, nloc, totele, ndim, cell_type, totele_str, totele_unst
    real, intent(in):: tnew(nloc,totele_str, totele_unst), x_all(ndim,nloc,totele),error(nloc,totele_str, totele_unst)
    real,intent(in) :: analytical(nloc,totele_str, totele_unst)
    character (len=20), intent(in):: solve_for

    !local variables
    character (len=25):: file_name
    character (len=10), allocatable:: str_coor(:,:,:)
    integer:: vtk_io, ele,s_ele, u_ele, i, nface
    allocate(str_coor(ndim+1,nloc,totele))

    do ele=1,totele
      do i=1,nloc
        write(str_coor(1,i,ele), '(F10.3)') x_all(1,i,ele)
        str_coor(1,i,ele)=sweep_blanks(str_coor(1,i,ele))
        write(str_coor(2,i,ele), '(F10.3)') x_all(2,i,ele)
        str_coor(2,i,ele)=sweep_blanks(str_coor(2,i,ele))
        write(str_coor(3,i,ele), '(F10.3)') 0.0
        str_coor(3,i,ele)=sweep_blanks(str_coor(3,i,ele))
      end do
    end do


    write(file_name, '(2a,i0,a)') trim(solve_for),'_', itime,'.vtu'
    open(unit=1, file=file_name)
      write(1,'(4(a,X))', advance="yes" ) '<VTKFile', 'type="UnstructuredGrid"', 'version="0.1"', 'byte_order="LittleEndian">'
      write(1,'(1(2X,a))', advance="yes" ) '<UnstructuredGrid>'
      write(1,'(4X,a,X,a,i0,a,X,a,i0,a,X)', advance="yes" ) '<Piece', 'NumberOfPoints="', nloc*totele,  '"',&
                                                                          'NumberOfCells="', totele, '">'
      write(1,'((6X,a,X,a))', advance="yes" ) '<PointData', 'Scalars="scalars">'

      !######################### values for your unknown #######################
      write(1,'((8X,a,X,a,X,a,a,a,X,a))', advance="yes" ) '<DataArray', 'type="Float32"', &
      'Name="', trim(solve_for),'"', 'Format="ascii">'

      do u_ele=1,totele_unst
        do s_ele =1,totele_str
        write(1,'(10x,F12.10,2x)',advance="yes") tnew(1,s_ele, u_ele)
        write(1,'(10x,F12.10,2x)',advance="yes") tnew(2,s_ele, u_ele)
        write(1,'(10x,F12.10,2x)',advance="yes") tnew(3,s_ele, u_ele)
      end do
    end do

      write(1,'((8X,a))', advance="yes" ) '</DataArray>'

      !###################### values of unkwon 2 if you have ####################
      write(1,'((8X,a,X,a,X,a,a,a,X,a))', advance="yes" ) '<DataArray', 'type="Float32"', &
      'Name="', 'error','"', 'Format="ascii">'

      do u_ele=1,totele_unst
        do s_ele =1,totele_str
        write(1,'(10x,F10.7,2x)',advance="yes") error(1,s_ele, u_ele)
        write(1,'(10x,F10.7,2x)',advance="yes") error(2,s_ele, u_ele)
        write(1,'(10x,F10.7,2x)',advance="yes") error(3,s_ele, u_ele)
      end do
    end do

      write(1,'((8X,a))', advance="yes" ) '</DataArray>'

      !###################### values of unkwon 2 if you have ####################
      write(1,'((8X,a,X,a,X,a,a,a,X,a))', advance="yes" ) '<DataArray', 'type="Float32"', &
      'Name="', 'analytical','"', 'Format="ascii">'

      do u_ele=1,totele_unst
        do s_ele =1,totele_str
        write(1,'(10x,F10.7,2x)',advance="yes") analytical(1,s_ele, u_ele)
        write(1,'(10x,F10.7,2x)',advance="yes") analytical(2,s_ele, u_ele)
        write(1,'(10x,F10.7,2x)',advance="yes") analytical(3,s_ele, u_ele)
      end do
    end do

      write(1,'((8X,a))', advance="yes" ) '</DataArray>'

      !####################### end if specifying unknown values #################
      write(1,'((6X,a))', advance="yes" ) '</PointData>'

      !############## point coordinates #########################################
      write(1,'((6X,a))', advance="yes" ) '<Points>'
      write(1,'((8X,a))', advance="yes" ) '<DataArray type="Float32" NumberOfComponents="3" Format="ascii">'
      do ele=1,totele
        do i=1,nloc
          write(1,'(10x, 3(a,x),2x)', advance="yes") adjustl(trim((str_coor(1,i,ele)))), adjustl(trim(str_coor(2,i,ele)))&
          ,adjustl(trim(str_coor(3,i,ele)))
        end do
      end do
      write(1,'((8X,a))', advance="yes" ) '</DataArray>'
      write(1,'((6X,a))', advance="yes" ) '</Points>'

      !############### cell connectivity points ################################
      write(1,'((6X,a))', advance="yes" ) '<Cells>'
      write(1,'((8X,a))', advance="yes" ) '<DataArray type="Int32" Name="connectivity" Format="ascii">'
      do ele=1,totele
        write(1,'(10x,i0,1x,i0,1x,i0)',advance="yes" ) ele*3-3,ele*3-2,ele*3-1
      end do
      write(1,'((8X,a))', advance="yes" ) '</DataArray>'

      !################ offsets ###############################################
      write(1,'((8X,a))', advance="yes" ) '<DataArray type="Int32" Name="offsets" Format="ascii">'
      write(1,'(10x,i0)',advance="no" ) nloc
      do ele=2,totele
        write(1,'(2x,i0)',advance="no" ) nloc*ele
      end do
      write(1,'(/,(8X,a))', advance="yes" ) '</DataArray>'

      !############################## cell types ###############################
      write(1,'((8X,a))', advance="yes" ) '<DataArray type="Int32" Name="types" Format="ascii">'
      write(1,'(10x,I0)', advance="no") cell_type
      do ele=2,totele
        write(1,'(x,I0)', advance="no") cell_type
      end do
      write(1,'(/,(8X,a))', advance="yes" ) '</DataArray>'

      write(1,'((6X,a))', advance="yes" ) '</Cells>'
      write(1,'((4X,a))', advance="yes" ) '</Piece>'
      write(1,'((2X,a))', advance="yes" ) '</UnstructuredGrid>'
      write(1,'((a))', advance="yes" ) '</VTKFile>'








    close(1)


  end subroutine get_vtu





















  !>@Brief
  ! cell_type is defined by vtk file formates, found in: (triangle =5)
  ! https://kitware.github.io/vtk-examples/site/VTKFileFormats/
  !cell_type is found in the link below (triangle=5)
  ! https://kitware.github.io/vtk-examples/site/VTKFileFormats/
  ! https://vtk.org/wp-content/uploads/2015/04/file-formats.pdf
  subroutine get_vtk(x_all, tnew, totele, nloc, itime, ndim, cell_type)

    ! Global variables
    integer, intent(in):: itime, nloc, totele, ndim, cell_type
    real, intent(in):: tnew(nloc,totele), x_all(ndim,nloc,totele)

    !local variables
    character (len=25):: file_name
    character (len=10), allocatable:: str_coor(:,:,:)
    integer:: vtk_io, ele, i, nface
    allocate(str_coor(ndim+1,nloc,totele))


    do ele=1,totele
      do i=1,nloc
        write(str_coor(1,i,ele), '(F10.3)') x_all(1,i,ele)
        str_coor(1,i,ele)=sweep_blanks(str_coor(1,i,ele))
        write(str_coor(2,i,ele), '(F10.3)') x_all(2,i,ele)
        str_coor(2,i,ele)=sweep_blanks(str_coor(2,i,ele))
        write(str_coor(3,i,ele), '(F10.3)') 0.0
        str_coor(3,i,ele)=sweep_blanks(str_coor(3,i,ele))
      end do
    end do


    write(file_name, '(a,i0,a)') 'VTK_', itime,'.vtk'
    open(unit=1, file=file_name)
      write(1,'(5(a,X))', advance="yes") '#', 'vtk', 'DataFile', 'Version', '2.0'
      write(1,'(3(a,X))', advance="yes") 'Unstructured', 'Grid', 'Example'
      write(1,'(a)', advance="yes") 'ASCII'
      write(1,'(2(a,X),//)', advance="yes") 'DATASET', 'UNSTRUCTURED_GRID'


      write(1,'(a,X,i0,x,a)', advance="yes") 'POINTS', nloc*totele, 'float'
      do ele=1,totele
        do i=1,nloc
          write(1,'(3(a,x),2x)', advance="no") adjustl(trim((str_coor(1,i,ele)))), adjustl(trim(str_coor(2,i,ele)))&
          ,adjustl(trim(str_coor(3,i,ele)))
        end do
      end do


      write(1,'(//,a,X,i0,x,i0)', advance="yes") 'CELLS', totele, totele*(nloc+1)
      do ele=1,totele
        write(1,'(i0,1x,i0,1x,i0,1x,i0)') nloc, ele*3-3,ele*3-2,ele*3-1
      end do


      write(1,'(//,a,X,i0)', advance="yes") 'CELL_TYPES', totele
      do ele=1,totele
        write(1,'(I0)') cell_type
      end do


      write(1,'(/,a,X,i0)', advance="yes") 'POINT_DATA', nloc*totele
      write(1,'(4(a,X))', advance="yes") 'SCALARS','scalars','float', '1'
      write(1,'(2(a,X))', advance="yes") 'LOOKUP_TABLE','default'
      do ele=1,totele
        write(1,'(F7.5,2x)',advance="no") tnew(1,ele)
        write(1,'(F7.5,2x)',advance="no") tnew(2,ele)
        write(1,'(F7.5,2x)',advance="no") tnew(3,ele)
      end do
    close(1)

    deallocate(str_coor)

  end subroutine get_vtk


  !cell_type is found in the link below (triangle=5)
  ! https://kitware.github.io/vtk-examples/site/VTKFileFormats/
  ! https://vtk.org/wp-content/uploads/2015/04/file-formats.pdf
  subroutine get_vtk_str_tri(tnew, totele, no_ele_row, no_ele_col, nface, &
                            dx, dy, nloc, itime, ndim, cell_type)

    ! External variables
    integer, intent(in):: itime, nloc, totele, ndim, cell_type, nface, no_ele_row, no_ele_col
    real, intent(in):: tnew(nloc,totele), dx, dy

    !local variables
    character (len=25):: file_name
    character (len=10), allocatable:: str_coor(:,:)
    integer:: vtk_io, ele, i
    real :: x_loc(ndim,nloc)
    allocate(str_coor(3,nloc))


    do ele=1,totele
      do i=1,nloc

      end do
    end do


    write(file_name, '(a,i0,a)') 'VTK_', itime,'.vtk'
    open(unit=1, file=file_name)
      write(1,'(5(a,X))', advance="yes") '#', 'vtk', 'DataFile', 'Version', '2.0'
      write(1,'(3(a,X))', advance="yes") 'Unstructured', 'Grid', 'Example'
      write(1,'(a)', advance="yes") 'ASCII'
      write(1,'(2(a,X),//)', advance="yes") 'DATASET', 'UNSTRUCTURED_GRID'


      write(1,'(a,X,i0,x,a)', advance="yes") 'POINTS', nloc*totele, 'float'
      do ele=1,totele
        do i=1,nloc

          call str_tri_X_nodes(ele, x_loc, ndim, nloc, dx, dy, no_ele_row, no_ele_col)
          write(str_coor(1,i), '(F10.3)') x_loc(1,i)
          str_coor(1,i)=sweep_blanks(str_coor(1,i))
          write(str_coor(2,i), '(F10.3)') x_loc(2,i)
          str_coor(2,i)=sweep_blanks(str_coor(2,i))
          write(str_coor(3,i), '(F10.3)') 0.0
          str_coor(3,i)=sweep_blanks(str_coor(3,i))

          write(1,'(3(a,x),2x)', advance="no") adjustl(trim((str_coor(1,i)))), adjustl(trim(str_coor(2,i)))&
          ,adjustl(trim(str_coor(3,i)))
        end do
      end do


      write(1,'(//,a,X,i0,x,i0)', advance="yes") 'CELLS', totele, totele*(nloc+1)
      do ele=1,totele
        write(1,'(i0,1x,i0,1x,i0,1x,i0)') nloc, ele*3-3,ele*3-2,ele*3-1
      end do


      write(1,'(//,a,X,i0)', advance="yes") 'CELL_TYPES', totele
      do ele=1,totele
        write(1,'(I0)') cell_type
      end do


      write(1,'(/,a,X,i0)', advance="yes") 'POINT_DATA', nloc*totele
      write(1,'(4(a,X))', advance="yes") 'SCALARS','scalars','float', '1'
      write(1,'(2(a,X))', advance="yes") 'LOOKUP_TABLE','default'
      do ele=1,totele
        write(1,'(F7.5,2x)',advance="no") tnew(1,ele)
        write(1,'(F7.5,2x)',advance="no") tnew(2,ele)
        write(1,'(F7.5,2x)',advance="no") tnew(3,ele)
      end do
    close(1)

    deallocate(str_coor)

  end subroutine get_vtk_str_tri




  character(30) function sweep_blanks(in_str)
    character(*), intent(in) :: in_str
    character(30) :: out_str
    character :: ch
    integer :: j

    out_str = " "
    do j=1, len_trim(in_str)
      ! get j-th char
      ch = in_str(j:j)
      if (ch .ne. " ") then
        out_str = trim(out_str) // ch
      endif
      sweep_blanks = out_str
    end do
  end function sweep_blanks

end module get_vtk_files
