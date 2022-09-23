module structured_meshgen
  use Msh2Tri
contains

  ! This subroutine generates structured rectangular mesh
  subroutine ele_info(totele, nface, face_ele, no_ele_row, row, row2, &
                           x_all, dx, dy, ndim, nloc, no_ele_col, col)
    ! ordering the face numbers: bottom face=1, right=1, left=3 and top=4
    ! row and row2 are row number associated with ele and ele2
    ! no_ele_row is total number of element in each row
    ! no_ele_col is total number of element in each column
    ! ndim is no of dimensions
    ! nloc is no of corner nodes
    ! face_ele(iface, ele) = given the face no iface and element no return the element next to
    ! the surface or if negative return the negative of the surface element number between element ele and face iface.

    implicit none
    integer, intent(in) :: no_ele_row, totele, no_ele_col, nloc, ndim, nface
    integer, intent(out) :: row, row2 , col
    integer:: face_ele(nface,totele), face_list_no(nface,totele), ele, iface

    real, intent(in) :: dx, dy
    real, intent(out) :: x_all(ndim,nloc,totele)

    do ele=1,totele
      row = ceiling(real(ele)/no_ele_row)
      col = ele-(no_ele_row*(row-1))

      ! corner node coordinates
      x_all(1,1,ele) = dx*(col-1); x_all(2,1,ele) = dy*(row-1)
      x_all(1,2,ele) = dx*col    ; x_all(2,2,ele) = dy*(row-1)
      x_all(1,3,ele) = dx*(col-1); x_all(2,3,ele) = dy*row
      x_all(1,4,ele) = dx*col    ; x_all(2,4,ele) = dy*row

      ! findin neighbiuring ele and face numbers
      do iface=1,nface
        if (iface==1) then
          face_ele(iface,ele) = ele - no_ele_row
          ! face_list_no(iface,ele) = 4
          ! row2 = ceiling(real(ele)/no_ele_row)
          ! if (row2.EQ.1) face_list_no(iface,ele) = -1*face_list_no(iface,ele)

        elseif ( iface==2 ) then
          face_ele(iface,ele) = ele - 1
          ! face_list_no(iface,ele) = 3
          row2 = ceiling(real(face_ele(iface,ele))/no_ele_row)
          if (row2.NE.row) then
            face_ele(iface,ele) = -1*face_ele(iface,ele)  !It is a boundary element located at the right side of the domain
            ! face_list_no(iface,ele) = -1*face_list_no(iface,ele)
          end if

        elseif ( iface==3 ) then
          face_ele(iface,ele) = ele +1   !It is a boundary element
          ! face_list_no(iface,ele) = 2
          row2 = ceiling(real(face_ele(iface,ele))/no_ele_row)
          if (row2.NE.row) then
            face_ele(iface,ele) = -1*face_ele(iface,ele)  !It is a boundary element located at the lest side of the domain
            ! face_list_no(iface,ele) = -1*face_list_no(iface,ele)
          end if

        elseif ( iface==4 ) then
          face_ele(iface,ele) = ele + no_ele_row
          ! face_list_no(iface,ele) = 1
          if (face_ele(iface,ele).GT.totele) then
            face_ele(iface,ele) = -1*face_ele(iface,ele)  !It is a boundary element located at the top of the domain
            ! face_list_no(iface,ele) = -1*face_list_no(iface,ele)
          end if
        end if
      end do
    end do
  end subroutine ele_info


  ! This subroutine generates structured triangular mesh
  subroutine tri_ele_info(totele, nface, face_ele, no_ele_row, &
                           x_all, dx, dy, ndim, nloc, no_ele_col)
    ! ordering the face numbers:
    !    |\        __1__
    !   2| \ 3     \   |
    !    |  \      3\  | 2
    !    |___\       \ |
    !      1          \|

    ! node numbering
    !    2
    !    |\        1___3
    !    | \       \   |
    !    |  \       \  |
    !    |___\       \ |
    !    3   1        \|
    !                  2
    ! row and row2 are row number associated with ele and ele2
    ! no_ele_row is total number of element in each row
    ! no_ele_col is total number of element in each column
    ! ndim is no of dimensions
    ! nloc is no of corner nodes
    ! face_ele(iface, ele) = given the face no iface and element no return the element next to
    ! the surface or if negative return the negative of the surface element number between element ele and face iface.

    implicit none
    ! External vbls
    integer, intent(in) :: no_ele_row, totele, no_ele_col, nloc, ndim, nface
    real, intent(in) :: dx, dy
    real, intent(out) :: x_all(ndim,nloc,totele)

    ! local vbls
    integer :: row, row2 , col
    integer:: face_ele(nface,totele), face_list_no(nface,totele), ele, iface


    do ele=1,totele
      row = ceiling(real(ele)/no_ele_row)
      col = ele-(no_ele_row*(row-1))

      ! corner node coordinates
      if ( MOD(ele,2)>0 ) then ! pointing up triangles
        x_all(1,1,ele) = dx*(col/2+1) ; x_all(2,1,ele) = dy*(row-1)
        x_all(1,2,ele) = dx*(col/2)   ; x_all(2,2,ele) = dy*row
        x_all(1,3,ele) = dx*(col/2)   ; x_all(2,3,ele) = dy*(row-1)
        ! x_all(1,4,ele) = dx*col    ; x_all(2,4,ele) = dy*row

        ! findin neighbiuring ele and face numbers
        do iface=1,nface
          if (iface==1) then
            face_ele(iface,ele) = ele - no_ele_row+1
            if (face_ele(iface,ele).LT.1) THEN
              face_ele(iface,ele) = 0!-1*face_ele(iface,ele)
            end if

          elseif ( iface==2 ) then
            face_ele(iface,ele) = ele - 1   !It is a boundary element
            ! face_list_no(iface,ele) = 2
            row2 = ceiling(real(face_ele(iface,ele))/no_ele_row)
            if (row2.NE.row) then
              face_ele(iface,ele) = 0!-1*face_ele(iface,ele)  !It is a boundary element located at the lest side of the domain
              ! face_list_no(iface,ele) = -1*face_list_no(iface,ele)
            end if

          elseif ( iface==3 ) then
            face_ele(iface,ele) = ele + 1
            ! face_list_no(iface,ele) = 3
            row2 = ceiling(real(face_ele(iface,ele))/no_ele_row)
            if (row2.NE.row) then
              face_ele(iface,ele) = 0!-1*face_ele(iface,ele)  !It is a boundary element located at the right side of the domain
              ! face_list_no(iface,ele) = -1*face_list_no(iface,ele)
            end if
          end if
        end do ! iface

      else ! pointing down triangles
        x_all(1,1,ele) = dx*(col/2-1) ; x_all(2,1,ele) = dy*row
        x_all(1,2,ele) = dx*(col/2)   ; x_all(2,2,ele) = dy*(row-1)
        x_all(1,3,ele) = dx*(col/2)   ; x_all(2,3,ele) = dy*row

        ! findin neighbiuring ele and face numbers
        do iface=1,nface
          if (iface==1) then
            face_ele(iface,ele) = ele + no_ele_row-1
            ! face_list_no(iface,ele) = 1
            if (face_ele(iface,ele).GT.totele) then
              face_ele(iface,ele) = 0!-1*face_ele(iface,ele)  !It is a boundary element located at the top of the domain
              ! face_list_no(iface,ele) = -1*face_list_no(iface,ele)
            end if

          elseif ( iface==2 ) then
            face_ele(iface,ele) = ele + 1
            ! face_list_no(iface,ele) = 3
            row2 = ceiling(real(face_ele(iface,ele))/no_ele_row)
            if (row2.NE.row) then
              face_ele(iface,ele) = 0!-1*face_ele(iface,ele)  !It is a boundary element located at the right side of the domain
              ! face_list_no(iface,ele) = -1*face_list_no(iface,ele)
            end if

          elseif ( iface==3 ) then
            face_ele(iface,ele) = ele - 1   !It is a boundary element
            ! face_list_no(iface,ele) = 2
            row2 = ceiling(real(face_ele(iface,ele))/no_ele_row)
            if (row2.NE.row) then
              face_ele(iface,ele) = 0!-1*face_ele(iface,ele)  !It is a boundary element located at the lest side of the domain
              ! face_list_no(iface,ele) = -1*face_list_no(iface,ele)
            end if
          end if
        end do ! iface
      end if
    end do ! ele
  end subroutine tri_ele_info


  ! This subroutine generates structured triangular mesh
  subroutine tri_ele_info2(totele, ele, iface, ele22, no_ele_row)
    ! ordering the face numbers:
    !    |\        __1__
    !   2| \ 3     \   |
    !    |  \      3\  | 2
    !    |___\       \ |
    !      1          \|

    ! node numbering
    !    2
    !    |\        1___3
    !    | \       \   |
    !    |  \       \  |
    !    |___\       \ |
    !    3   1        \|
    !                  2
    ! row and row2 are row number associated with ele and ele2
    ! no_ele_row is total number of element in each row
    ! no_ele_col is total number of element in each column
    ! ndim is no of dimensions
    ! nloc is no of corner nodes
    ! face_ele(iface, ele) = given the face no iface and element no return the element next to
    ! the surface or if negative return the negative of the surface element number between element ele and face iface.
    implicit none
    ! External vbls
    integer, intent(in) :: no_ele_row, totele
    ! local vbls
    integer :: row, row2 , col
    integer:: ele22, ele, iface


    row = ceiling(real(ele)/no_ele_row)
    col = ele-(no_ele_row*(row-1))
    ! corner node coordinates
    ! call str_tri_X_nodes(ele, x_loc, ndim, nloc, dx, dy, no_ele_row, no_ele_col)

    if ( MOD(ele,2)>0 ) then ! pointing up triangles
      ! findin neighbiuring ele and face numbers
        if (iface==1) then
          ele22 = ele - no_ele_row+1
          if (ele22.LT.1) THEN
            ele22 = 0!-1*face_ele(iface,ele)
          end if

        elseif ( iface==2 ) then
          ele22 = ele - 1   !It is a boundary element
          row2 = ceiling(real(ele22)/no_ele_row)
          if (row2.NE.row) then
            ele22 = 0!-1*face_ele(iface,ele)  !It is a boundary element located at the lest side of the domain
          end if

        elseif ( iface==3 ) then
          ele22 = ele + 1
          row2 = ceiling(real(ele22)/no_ele_row)
          if (row2.NE.row) then
            ele22 = 0!-1*face_ele(iface,ele)  !It is a boundary element located at the right side of the domain
          end if
        end if

    else ! pointing down triangles
      ! findin neighbiuring ele and face numbers
        if (iface==1) then
          ele22 = ele + no_ele_row-1
          if (ele22.GT.totele) then
            ele22 = 0!-1*face_ele(iface,ele)  !It is a boundary element located at the top of the domain
          end if

        elseif ( iface==2 ) then
          ele22 = ele + 1
          row2 = ceiling(real(ele22)/no_ele_row)
          if (row2.NE.row) then
            ele22 = 0!-1*face_ele(iface,ele)  !It is a boundary element located at the right side of the domain
          end if

        elseif ( iface==3 ) then
          ele22 = ele - 1   !It is a boundary element
          row2 = ceiling(real(ele22)/no_ele_row)
          if (row2.NE.row) then
            ele22 = 0!-1*face_ele(iface,ele)  !It is a boundary element located at the lest side of the domain
          end if
        end if
    end if
  end subroutine tri_ele_info2


  ! this subroutine returns x_loc only
  subroutine str_tri_X_nodes(ele, x_loc, ndim, nloc, dx, dy, no_ele_row, no_ele_col)
    implicit none
    ! external vbls
    integer, intent(in) :: ele, ndim, nloc, no_ele_row, no_ele_col
    real, intent(in) :: dx, dy
    !local vbls
    real, intent(inout) :: x_loc(ndim, nloc)
    integer:: row, col

    row = ceiling(real(ele)/no_ele_row)
    col = ele-(no_ele_row*(row-1))

    if ( MOD(ele,2)>0 ) then ! pointing up triangles
      x_loc(1,1) = dx*(col/2+1) ; x_loc(2,1) = dy*(row-1)
      x_loc(1,2) = dx*(col/2)   ; x_loc(2,2) = dy*row
      x_loc(1,3) = dx*(col/2)   ; x_loc(2,3) = dy*(row-1)

    else
      x_loc(1,1) = dx*(col/2-1) ; x_loc(2,1) = dy*row
      x_loc(1,2) = dx*(col/2)   ; x_loc(2,2) = dy*(row-1)
      x_loc(1,3) = dx*(col/2)   ; x_loc(2,3) = dy*row
    end if

  end subroutine str_tri_X_nodes

  ! This subroutine returns x_all and face_ele forunstructured triangular mesh
  ! subroutine unst_tri_ele_info(totele, nface, face_ele, x_all, ndim, nloc)
  !   ! ndim is no of dimensions
  !   ! nloc is no of corner nodes
  !   ! face_ele(iface, ele) = given the face no iface and element no return the element next to
  !   ! the surface or if negative return the negative of the surface element number between element ele and face iface.
  !
  !   implicit none
  !   integer :: totele, nloc, ndim, nface, ierr, ele, iface
  !   type(Mesh), allocatable, dimension(:) :: meshList
  !   type(Mesh):: meshList
  !   integer, allocatable:: face_ele(:,:), face_list_no(:,:)
  !   real, allocatable:: x_all(:,:,:)
  !
  !   call ReadMSH(meshList,'./test_msh',ierr)
  !   totele = size(meshlist)
  !   allocate(face_ele(nface,totele), face_list_no(nface,totele),x_all(ndim,nloc,totele))
  !   ! print*, meshList(1)%Xp1(2)
  !   do ele=1,totele
  !     ! corner node coordinates
  !       ! x_all(1,1,ele) = meshList(ele)%Xp1(1)   ; x_all(2,1,ele) = meshList(ele)%Xp1(2)
  !       ! x_all(1,2,ele) = meshList(ele)%Xp2(1)   ; x_all(2,2,ele) = meshList(ele)%Xp2(2)
  !       ! x_all(1,3,ele) = meshList(ele)%Xp3(1)   ; x_all(2,3,ele) = meshList(ele)%Xp3(2)
  !
  !       ! findin neighbiuring ele and face numbers
  !       do iface=1,nface
  !           ! face_ele(iface,ele) = meshList(ele)%Neig(iface)
  !       end do ! iface
  !   end do ! ele
  ! end subroutine unst_tri_ele_info



end module structured_meshgen
