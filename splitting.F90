module splitting
  use Msh2Tri

implicit none
contains



  ! brief:: this subroutine restricts the fine mehs into a coarser one
  subroutine restrictor(tracer, totele_unst, i_split, ilevel)
    implicit none
    ! global vbl
    integer, intent(in) ::  i_split, totele_unst,ilevel
    type(fields), dimension(:), INTENT(INOUT) :: tracer
    ! local vbl
    integer :: fin_ele(4),totele_str_dummy,coarse_ele, un_ele

    if ( ilevel>1 ) then

      totele_str_dummy = 4**(i_split)
      do un_ele=1,totele_unst
        do coarse_ele = 1,totele_str_dummy
          call element_conversion(fin_ele,coarse_ele,i_split)
          !################## now doing restriction #########################################

          tracer(ilevel)%RHS(1,coarse_ele,un_ele) = average(tracer(ilevel-1)%residuale(:,fin_ele(3),un_ele))
          tracer(ilevel)%RHS(2,coarse_ele,un_ele) = average(tracer(ilevel-1)%residuale(:,fin_ele(4),un_ele))
          tracer(ilevel)%RHS(3,coarse_ele,un_ele) = average(tracer(ilevel-1)%residuale(:,fin_ele(1),un_ele))
        end do
      end do
    end if
  end subroutine restrictor



  subroutine prolongator(tracer, totele_unst, i_split, ilevel)
    implicit none
    ! global vbl
    integer, intent(in) ::  i_split,ilevel,totele_unst
    type(fields), dimension(:), INTENT(INOUT) :: tracer
    ! local vbl
    integer :: fin_ele(4), un_ele, totele_str,coarse_ele

    totele_str = 4**(i_split)
    do un_ele = 1,totele_unst
      do coarse_ele = 1, totele_str
        call element_conversion(fin_ele,coarse_ele,i_split-1)

        !################## now doing restriction #########################################
        ! should be == tnew(iloc) + 0.5*error() + 0.5*error()
        ! basically, here error is coming from solved coarser grid
        !####################### for fin_ele(1) ##################################
        tracer(ilevel)%tnew(1,fin_ele(1),un_ele) = tracer(ilevel)%tnew(1,fin_ele(1),un_ele) &
                                                       + 0.5* tracer(ilevel)%error(3,coarse_ele,un_ele)&
                                                       + 0.5* tracer(ilevel)%error(1,coarse_ele,un_ele)
        tracer(ilevel)%tnew(2,fin_ele(1),un_ele) = tracer(ilevel)%tnew(2,fin_ele(1),un_ele) &
                                                       + 0.5* tracer(ilevel)%error(2,coarse_ele,un_ele)&
                                                       + 0.5* tracer(ilevel)%error(3,coarse_ele,un_ele)
        tracer(ilevel)%tnew(3,fin_ele(1),un_ele) = tracer(ilevel)%tnew(3,fin_ele(1),un_ele) &
                                                       + tracer(ilevel)%error(3,coarse_ele,un_ele)
        !####################### for fin_ele(2) ##################################
        tracer(ilevel)%tnew(1,fin_ele(2),un_ele) = tracer(ilevel)%tnew(1,fin_ele(2),un_ele) &
                                                       + tracer(ilevel)%error(2,fin_ele(1),un_ele)
        tracer(ilevel)%tnew(2,fin_ele(2),un_ele) = tracer(ilevel)%tnew(2,fin_ele(2),un_ele) &
                                                       + tracer(ilevel)%error(1,fin_ele(1),un_ele)
        tracer(ilevel)%tnew(3,fin_ele(2),un_ele) = tracer(ilevel)%tnew(3,fin_ele(2),un_ele) &
                                                       + 0.5* tracer(ilevel)%error(1,coarse_ele,un_ele)&
                                                       + 0.5* tracer(ilevel)%error(2,coarse_ele,un_ele)
        !####################### for fin_ele(3) ##################################
        tracer(ilevel)%tnew(1,fin_ele(3),un_ele) = tracer(ilevel)%tnew(1,fin_ele(3),un_ele) &
                                                       + tracer(ilevel)%error(1,coarse_ele,un_ele)
        tracer(ilevel)%tnew(2,fin_ele(3),un_ele) = tracer(ilevel)%tnew(2,fin_ele(3),un_ele) &
                                                       + tracer(ilevel)%error(3,fin_ele(2),un_ele)
        tracer(ilevel)%tnew(3,fin_ele(3),un_ele) = tracer(ilevel)%tnew(3,fin_ele(3),un_ele) &
                                                       + tracer(ilevel)%error(2,fin_ele(2),un_ele)
        !####################### for fin_ele(4) ##################################
        tracer(ilevel)%tnew(1,fin_ele(4),un_ele) = tracer(ilevel)%tnew(1,fin_ele(4),un_ele) &
                                                       + tracer(ilevel)%error(3,fin_ele(2),un_ele)
        tracer(ilevel)%tnew(2,fin_ele(4),un_ele) = tracer(ilevel)%tnew(2,fin_ele(4),un_ele) &
                                                       + tracer(ilevel)%error(2,coarse_ele,un_ele)
        tracer(ilevel)%tnew(3,fin_ele(4),un_ele) = tracer(ilevel)%tnew(3,fin_ele(4),un_ele) &
                                                       + tracer(ilevel)%error(1,fin_ele(2),un_ele)
      end do
    end do
  end subroutine prolongator



  !brief:: this subroutine returns local number of fine eles at fines multigrid level corresponding to the coarser coarse_ele
  ! n_split:: refers to no spliting of the coarser grid NOT the finer one
  subroutine element_conversion(fin_ele,coarse_ele,i_split)
    implicit none
    ! global vbl
    integer, intent(in) ::  coarse_ele,i_split
    ! local vbl
    integer :: counter,orientation,ipos, irow, rowx,tot_fine
    integer, intent(out) :: fin_ele(4)

    call get_str_info(i_split, coarse_ele, irow, ipos, orientation)
    tot_fine = 0
    rowx = 2**(i_split+1)*2-1

    select case(orientation)
    case(1) ! uptriangle    4
            !            1  2  3
        counter = 2
        do while (counter < irow*2)
          tot_fine = tot_fine + rowx
          rowx = rowx -2
          counter = counter +1
        end do

        fin_ele(1) = ipos * 2 -1 + tot_fine
        fin_ele(2) = fin_ele(1) +1
        fin_ele(3) = fin_ele(1) +2

        tot_fine = tot_fine + rowx
        fin_ele(4) = ipos * 2 -1 + tot_fine

      case(0) ! down triangle   3  2  1
              !                    4
        counter = 1
        do while (counter < irow*2)
          tot_fine = tot_fine + rowx
          rowx = rowx -2
          counter = counter +1
        end do

        fin_ele(3) = (ipos/2-1)*3 + ipos/2 + tot_fine +1
        fin_ele(2) = fin_ele(3) + 1
        fin_ele(1) = fin_ele(3) + 2
        fin_ele(4) = fin_ele(1) - rowx-2
    end select
  end subroutine element_conversion





  real function average(numbers) result(ave)
    implicit none
    real, intent(in) :: numbers(3)

    ave = sum(numbers(:))/3.
  end function average








  !@brief:: this sub, generates structured ele numbers located on the iface of each un_ele in semi-str mesh
  ! subroutine get_surface_ele(meshlist, n_split, totele_unst, nface)
  !   implicit none
  !   type(Mesh), dimension(:), INTENT(INOUT) :: meshList
  !   integer, intent(in) :: n_split, nface, totele_unst
  !   integer :: A, B, i, j, k, un_ele
  !
  !   do un_ele = 1,totele_unst
  !     ! for iface==1
  !     A = (un_ele-1)*4**n_split + 1
  !     do k = 2**n_split,1,-1
  !       meshList(un_ele)%s_ele(k,1) = A
  !       A=A+2
  !     end do
  !     ! for iface==2
  !     B = (un_ele -1)*4**n_split
  !     i=2*2**n_split -1 ; j=i
  !     do k = 2**n_split,1,-1
  !       meshList(un_ele)%s_ele(k,2) = B+i
  !       i=i+j-2
  !       j=j-2
  !     end do
  !     ! for iface==3
  !     i=1 ; j=1
  !     do k = 1,2**n_split
  !       meshList(un_ele)%s_ele(k,3) = B+i
  !       i=i+2*2**n_split-j
  !       j=j+2
  !     end do
  !   end do
  ! end subroutine get_surface_ele


! brief:: this subroutine return global surface elements of my neig
  subroutine get_surface_ele(meshlist, n_split, totele_unst, nface)
    implicit none
    type(Mesh), dimension(:), INTENT(INOUT) :: meshList
    integer, intent(in) :: n_split, nface, totele_unst
    integer :: A, B, i, j, k, un_ele, iface, Npos, Nside, lnodes2(2)


    do un_ele=1,totele_unst
      do iface=1,nface
        call getNeigDataMesh(meshlist, un_ele, iface, Npos, Nside, lnodes2)
        if ( Npos/=0 ) then
            if ( iface==1 ) then
              if ( meshList(un_ele)%Dir(1) ) then
                    if ( Nside==1 ) then
                      A = (Npos-1)*4**n_split + 1
                      do k = 2**n_split,1,-1
                        meshList(un_ele)%s_ele(k,1) = A
                        A=A+2
                      end do
                    elseif ( Nside==2 ) then
                      B = (Npos -1)*4**n_split
                      i=2*2**n_split -1 ; j=i
                      do k = 2**n_split,1,-1
                        meshList(un_ele)%s_ele(k,1) = B+i
                        i=i+j-2
                        j=j-2
                      end do
                    elseif ( Nside==3 ) then
                      B = (Npos -1)*4**n_split
                      i=1 ; j=1
                      do k = 1,2**n_split
                        meshList(un_ele)%s_ele(k,1) = B+i
                        i=i+2*2**n_split-j
                        j=j+2
                      end do
                    end if
              else
                    if ( Nside==1 ) then
                      A = (Npos-1)*4**n_split + 1
                      do k = 1,2**n_split
                        meshList(un_ele)%s_ele(k,1) = A
                        A=A+2
                      end do
                    elseif ( Nside==2 ) then
                      B = (Npos -1)*4**n_split
                      i=2*2**n_split -1 ; j=i
                      do k = 1, 2**n_split
                        meshList(un_ele)%s_ele(k,1) = B+i
                        i=i+j-2
                        j=j-2
                      end do
                    elseif ( Nside==3 ) then
                      B = (Npos -1)*4**n_split
                      i=1 ; j=1
                      do k = 2**n_split,1,-1
                        meshList(un_ele)%s_ele(k,1) = B+i
                        i=i+2*2**n_split-j
                        j=j+2
                      end do
                    end if


              end if ! if dir(iface)

            elseif ( iface==2 ) then
                if ( meshList(un_ele)%Dir(2) ) then
                      if ( Nside==1 ) then
                        A = (Npos-1)*4**n_split + 1
                        do k = 2**n_split,1,-1
                          meshList(un_ele)%s_ele(k,2) = A
                          A=A+2
                        end do
                      elseif ( Nside==2 ) then
                        B = (Npos -1)*4**n_split
                        i=2*2**n_split -1 ; j=i
                        do k = 1,2**n_split
                          meshList(un_ele)%s_ele(k,2) = B+i
                          i=i+j-2
                          j=j-2
                        end do
                      elseif ( Nside==3 ) then
                        B = (Npos -1)*4**n_split
                        i=1 ; j=1
                        do k = 2**n_split,1,-1
                          meshList(un_ele)%s_ele(k,2) = B+i
                          i=i+2*2**n_split-j
                          j=j+2
                        end do
                      end if
                else
                      if ( Nside==1 ) then
                        A = (Npos-1)*4**n_split + 1
                        do k = 1,2**n_split
                          meshList(un_ele)%s_ele(k,2) = A
                          A=A+2
                        end do
                      elseif ( Nside==2 ) then
                        B = (Npos -1)*4**n_split
                        i=2*2**n_split -1 ; j=i
                        do k = 2**n_split,1,-1
                          meshList(un_ele)%s_ele(k,2) = B+i
                          i=i+j-2
                          j=j-2
                        end do
                      elseif ( Nside==3 ) then
                        B = (Npos -1)*4**n_split
                        i=1 ; j=1
                        do k = 1,2**n_split
                          meshList(un_ele)%s_ele(k,2) = B+i
                          i=i+2*2**n_split-j
                          j=j+2
                        end do
                      end if


                end if ! if dir(iface)

            elseif ( iface==3 ) then
              if ( meshList(un_ele)%Dir(3) ) then
                    if ( Nside==1 ) then
                          A = (Npos-1)*4**n_split + 1
                          do k = 1,2**n_split
                            meshList(un_ele)%s_ele(k,3) = A
                            A=A+2
                          end do
                    elseif ( Nside==2 ) then
                          B = (Npos -1)*4**n_split
                          i=2*2**n_split -1 ; j=i
                          do k = 2**n_split,1,-1
                            meshList(un_ele)%s_ele(k,3) = B+i
                            i=i+j-2
                            j=j-2
                          end do
                    elseif ( Nside==3 ) then
                          B = (Npos -1)*4**n_split
                          i=1 ; j=1
                          do k = 1,2**n_split
                            meshList(un_ele)%s_ele(k,3) = B+i
                            i=i+2*2**n_split-j
                            j=j+2
                          end do
                    end if
              else
                    if ( Nside==1 ) then
                          A = (Npos-1)*4**n_split + 1
                          do k = 2**n_split,1,-1
                            meshList(un_ele)%s_ele(k,3) = A
                            A=A+2
                          end do
                    elseif ( Nside==2 ) then
                          B = (Npos -1)*4**n_split
                          i=2*2**n_split -1 ; j=i
                          do k = 1,2**n_split
                            meshList(un_ele)%s_ele(k,3) = B+i
                            i=i+j-2
                            j=j-2
                          end do
                    elseif ( Nside==3 ) then
                          B = (Npos -1)*4**n_split
                          i=1 ; j=1
                          do k = 2**n_split,1,-1
                            meshList(un_ele)%s_ele(k,3) = B+i
                            i=i+2*2**n_split-j
                            j=j+2
                          end do
                    end if


              end if ! if dir(iface)
            end if ! if iface
        end if ! if Npos/=0
      end do ! iface
    end do ! un_ele
    ! do un_ele = 1,totele_unst
    !   ! for iface==1
    !   A = (un_ele-1)*4**n_split + 1
    !   do k = 1,2**n_split
    !     meshList(un_ele)%s_ele(k,1) = A
    !     A=A+2
    !   end do
    !   ! for iface==2
    !   B = (un_ele -1)*4**n_split
    !   i=2*2**n_split -1 ; j=i
    !   do k = 1, 2**n_split
    !     meshList(un_ele)%s_ele(k,2) = B+i
    !     i=i+j-2
    !     j=j-2
    !   end do
    !   ! for iface==3
    !   i=1 ; j=1
    !   do k = 1,2**n_split
    !     meshList(un_ele)%s_ele(k,3) = B+i
    !     i=i+2*2**n_split-j
    !     j=j+2
    !   end do
    ! end do
  end subroutine get_surface_ele


  !brief>:: this subroutine retuns an 2D array of local structured boundary element numbers starting from face 1, 2 and 3
subroutine loc_surf_ele2(surf_ele, n_split)
  implicit none
  !global vbl
  integer, intent(in) :: n_split
  integer, allocatable, intent(inout) :: surf_ele(:,:)
  integer :: i, ele, counter
  if (allocated(surf_ele)) deallocate(surf_ele)
  allocate(surf_ele(2**n_split,3))
  !localvbl

  ! face 1
  i=2
  surf_ele(1,1) = 1
  do ele=2,2**n_split
    surf_ele(ele,1) = surf_ele(ele-1,1) + 2
  end do

  ! face 2 & 3
  surf_ele(1,3) = 1
  counter=surf_ele(ele-1,1)
  surf_ele(1,2) = counter
  do ele=2,2**n_split
    surf_ele(ele,2) = surf_ele(ele-1,2) + counter-2
    surf_ele(ele,3) = surf_ele(ele-1,2) +1
    counter = counter -2
  end do

end subroutine loc_surf_ele2



!brief>:: this subroutine retuns an 2D array of local structured boundary element numbers starting from face 1, 2 and 3
! it for multigrid
subroutine loc_surf_ele_multigrid(ele_info, ilevel, n_split)
implicit none
!global vbl
type(element_info), intent(inout) :: ele_info(:)
integer, intent(in) :: n_split, ilevel
integer :: i, ele, counter

! face 1
i=2
ele_info(ilevel)%surf_ele(1,1) = 1
do ele=2,2**n_split
  ele_info(ilevel)%surf_ele(ele,1) = ele_info(ilevel)%surf_ele(ele-1,1) + 2
end do

! face 2 & 3
ele_info(ilevel)%surf_ele(1,3) = 1
counter=ele_info(ilevel)%surf_ele(ele-1,1)
ele_info(ilevel)%surf_ele(1,2) = counter
do ele=2,2**n_split
  ele_info(ilevel)%surf_ele(ele,2) = ele_info(ilevel)%surf_ele(ele-1,2) + counter-2
  ele_info(ilevel)%surf_ele(ele,3) = ele_info(ilevel)%surf_ele(ele-1,2) +1
  counter = counter -2
end do

end subroutine loc_surf_ele_multigrid


!brief>:: this subroutine retuns an 1D array of local structured boundary element numbers starting from face 1, 2 and 3
!e.g. for n_split=2, it gives: 1,3,5,7,7,12,15,16,1,8,13,16
subroutine loc_surf_ele(surf_ele, n_split)
  implicit none
  !global vbl
  integer, intent(in) :: n_split
  integer, dimension(:), intent(inout) :: surf_ele
  !localvbl
  integer :: i, ele, counter,j, k

  ! face 1
  i=2
  surf_ele(1) = 1
  do ele=2,2**n_split
    surf_ele(ele) = surf_ele(ele-1) + 2
  end do

  counter=surf_ele(ele-1)
  surf_ele(ele) = surf_ele(ele-1)

  surf_ele(ele+2**n_split+1) = surf_ele(ele-1)+1

  do j=ele+1,ele+2**n_split-1
    surf_ele(j) = surf_ele(j-1) + counter-2
    counter = counter -2
  end do

  surf_ele(j) = surf_ele(1)

  do k=j+1,j+2**n_split-1
    surf_ele(k) = surf_ele(k-2**n_split-1)+ 1
  end do


  ! surf_ele(ele,3) = surf_ele(ele-1,2) +1
end subroutine loc_surf_ele




  !@brief> this subroutine recieves an unstructured mesh (meshlist) and splits each ele 'n' times
  ! and returns a semi_structured mesh
!   subroutine get_splitting2(meshlist, n, semi_meshlist)
!     ! External vbls
!     type(Mesh), dimension(:) :: meshList
!     integer, intent(in):: n
!
!     ! local vbls
!     ! type(Semi_Mesh), allocatable, DIMENSION(:), INTENT(OUT) :: semi_meshlist
!     integer :: Un_ele, ele, i,j,k, idim, ndim=2, nloc=3, iloc, no_ups_row,totele_up,totele_down
!     integer:: no_up_row1,i_rowele, irow, row_num, no_down_row, no_down_row1, totele_unst
!     integer:: total, next, current, totele_str
!     real, allocatable :: x_down(:,:,:), x_up(:,:,:), disp(:,:), x_str(:,:,:)
!
!     totele_unst = size(meshlist)
!     no_ups_row = 2**n ; no_up_row1 = 2**n
!     no_down_row = 2**n-1 ; no_down_row1 = 2**n-1
!     totele_up = no_ups_row *(no_ups_row+1)*0.5
!     totele_down = 4**n - totele_up!(no_ups_row-1) * (no_ups_row)*0.5
!     totele_str = totele_down+totele_up
!     allocate(x_down(ndim,nloc,totele_down), x_up(ndim,nloc,totele_up), disp(ndim, ndim))
!     allocate( semi_meshlist(size(meshlist)) )
!
!     do i=1,totele_unst
!       allocate(semi_meshlist(i)%up(ndim,nloc,totele_unst))
!       allocate(semi_meshlist(i)%down(ndim,nloc,totele_unst))
!       allocate(semi_meshlist(i)%x_str(ndim,nloc,totele_str))
!     end do
!
!
!     do Un_ele=1,totele_unst ! loop over unstr ele
!       no_ups_row = 2**n; no_up_row1 = 2**n
!       no_down_row = 2**n-1; no_down_row1 = 2**n-1
!       ! up triangels
!       ! it just defines the first ele in each unstr ele
!       do idim=1,ndim
!         x_up(idim,1,1)= (meshlist(un_ele)%X(idim,1) - meshlist(un_ele)%X(idim,3))/(2**n)&
!                                                     +meshlist(un_ele)%X(idim,3)
!         x_up(idim,2,1)= (meshlist(un_ele)%X(idim,2) - meshlist(un_ele)%X(idim,3))/(2**n)&
!                                                     +meshlist(un_ele)%X(idim,3)
!         x_up(idim,3,1)= meshlist(un_ele)%X(idim,3)
!       end do
!
!       ! vector-displacement of the side with points 1 and 3
!       disp(1,1) = ((x_up(1,1,1)-x_up(1,3,1)))
!       disp(1,2) = ((x_up(2,1,1)-x_up(2,3,1)))
!       ! vector-displacement of the side with points 1 and 2
!       disp(2,1) = ((x_up(1,2,1)-x_up(1,3,1)))
!       disp(2,2) = ((x_up(2,2,1)-x_up(2,3,1)))
!
!
!       row_num = 2**n
!       i_rowele = 2**n
!       ele=2
!       do while( row_num>1 )
!         do while ( ele<= i_rowele )
!           do idim=1,ndim
!             x_up(idim,:,ele) = x_up(idim,:,ele-1) + disp(1,idim)
!           end do
!           ele = ele +1
!         end do
!
!         ! calculate 1st ele from next row
!         do iloc=1,nloc
!           x_up(:,iloc,ele) = x_up(:,iloc,ele-no_ups_row) + disp(2,:)
!         end do
!         ele = ele+1
!         no_ups_row = no_ups_row -1
!         no_up_row1=no_up_row1-1
!         i_rowele = i_rowele + no_up_row1
!         row_num = row_num -1
!       end do ! end do while row_num
!
!
!       ! down triangels
!       ! defining the 1st down tri
!       do idim = 1,ndim
!         x_down(idim,1,1) = x_up(idim,2,1)
!         x_down(idim,2,1) = x_up(idim,1,1)
!         x_down(idim,3,1) = x_up(idim,2,1) + disp(1,idim)
!       end do
!
!       row_num = 2**n
!       i_rowele = 2**n -1
!       ele=2
!       do while (row_num>2)
!         do while ( ele<=i_rowele)
!           do idim=1,ndim
!             x_down(idim,:,ele) = x_down(idim,:,ele-1) + disp(1,idim)
!           end do
!           ele = ele +1
!         end do
!
!         ! calculate 1st ele from next row
!         do idim=1,ndim
!           x_down(idim,:,ele) = x_down(idim,:,ele-no_down_row) + disp(2,idim)
!         end do
!
!         ele = ele+1
!         no_down_row = no_down_row -1
!         no_down_row1=no_down_row1-1
!         i_rowele = i_rowele + no_down_row1
!         row_num = row_num -1
!
!       end do
!
!       semi_meshlist(un_ele)%down = x_down
!       semi_meshlist(un_ele)%up = x_up
!
!       irow = 2**n
!       i_rowele = 2*2**n-1
!       i=1 ; j=1 ; k=1
!       total= 2**(n+1) -1 ; next = total ; current = next
!
!       do while ( irow>=1)
!         do iloc=1,nloc
!           do idim=1,ndim
!             semi_meshlist(un_ele)%x_str(idim, iloc, i)= x_up(idim,iloc,j)
!           end do
!         end do
!         i=i+1 ; j=j+1
!         if ( i < total ) then
!           do iloc=1,nloc
!             do idim=1,ndim
!               semi_meshlist(un_ele)%x_str(idim,iloc,i) = x_down(idim,iloc,k)
!             end do
!           end do
!           i=i+1 ; k=k+1
!         else
!           irow = irow - 1
!           i_rowele = i_rowele - 2
!           next = current -2
!           total = total +next
!           current = next
!         end if
!       end do ! end while
!
! ! if ( un_ele==1 ) then
! !   print*, 'un_ele = ',un_ele
! !   print*, 'ups'
! !   do ele=1,totele_up
! !     print*, ele
! !     print*, x_up(1,:,ele)
! !     print*, x_up(2,:,ele)
! !   end do
! !
! !   print*, 'down'
! !   do ele=1,totele_down
! !     print*, ele
! !     print*, x_down(1,:,ele)
! !     print*, x_down(2,:,ele)
! !   end do
! !
! ! end if
!
!     end do !end loop over unstr ele
!
!     deallocate(x_up, x_down, disp)
!   end subroutine get_splitting2







!@brief> this subroutine gets number of splitting of an unstr ele and returns neighbours of
! each internal str ele after splitting
! We use it once to create an stencil for all unstr_ele.
! local face numbers:
!    |\        __1__
!   2| \ 3     \   |
!    |  \      3\  | 2
!    |___\       \ |
!      1          \|
subroutine get_str_neig(n, str_neig)
  !external vlbs
  integer, intent(in):: n

  !local vbls
  integer:: total, current, irow, ele, totele
  integer, allocatable, dimension(:), intent(inout):: str_neig(:,:)

  totele = 4**n
  if (allocated(str_neig)) deallocate(str_neig)
  allocate(str_neig(3,totele))

  total = 2**(n+1)-1; current=total
  irow = 2**n -1

  !ROW 1
  str_neig(1,1)=0 ; str_neig(2,1)=0 ; str_neig(3,1)=2 ! only ele 1
  ele=2

  do while ( ele <= total)
    str_neig(2,ele)= ele+1 ; str_neig(3,ele)=ele-1 ; str_neig(1,ele)=ele+total -1
    ele = ele+1
    str_neig(2,ele)=ele-1 ; str_neig(1,ele)=0 ; str_neig(3,ele)=ele+1
    ele = ele+1
  end do

  str_neig(3,ele-1)=0 ! last ele in the row

  do while (irow >= 1)
    total = total +current -2
    current = current -2

    str_neig(2,ele)=0 ; str_neig(3,ele)=ele+1 ; str_neig(1,ele)=ele - current -1 ! 1st ele in new row

    ele = ele+1
    do while ( ele <= total )
      str_neig(2,ele)= ele+1 ; str_neig(3,ele)=ele-1 ; str_neig(1,ele)=ele+current -1
      ele = ele+1

      str_neig(2,ele)=ele-1 ; str_neig(3,ele)=ele+1 ; str_neig(1,ele)=ele - current -1
      ele = ele+1
    end do
    str_neig(3,ele-1)=0 ! last ele in the row
    irow = irow -1
  end do

end subroutine get_str_neig






!@brief> this subroutine gets number of splitting of an unstr ele and returns neighbours of
! each internal str ele after splitting
! We use it once to create an stencil for all unstr_ele.
! it is for multigrid version
! local face numbers:
!    |\        __1__
!   2| \ 3     \   |
!    |  \      3\  | 2
!    |___\       \ |
!      1          \|
subroutine get_str_neig_multigrid(ele_info, ilevel, n)
  !external vlbs
  integer, intent(in):: n, ilevel
  type(element_info), allocatable :: ele_info(:)

  !local vbls
  integer:: total, current, irow, ele


  total = 2**(n+1)-1; current=total
  irow = 2**n -1

  !ROW 1
  ele_info(ilevel)%str_neig(1,1)=0 ; ele_info(ilevel)%str_neig(2,1)=0 ; ele_info(ilevel)%str_neig(3,1)=2 ! only ele 1
  ele=2

  do while ( ele <= total)
    ele_info(ilevel)%str_neig(2,ele)= ele+1; ele_info(ilevel)%str_neig(3,ele)=ele-1; ele_info(ilevel)%str_neig(1,ele)=ele+total-1
    ele = ele+1
    ele_info(ilevel)%str_neig(2,ele)=ele-1 ; ele_info(ilevel)%str_neig(1,ele)=0 ; ele_info(ilevel)%str_neig(3,ele)=ele+1
    ele = ele+1
  end do

  ele_info(ilevel)%str_neig(3,ele-1)=0 ! last ele in the row

  do while (irow >= 1)
    total = total +current -2
    current = current -2

    ! 1st ele in new row
    ele_info(ilevel)%str_neig(2,ele)=0; ele_info(ilevel)%str_neig(3,ele)=ele+1; ele_info(ilevel)%str_neig(1,ele)=ele-current-1

    ele = ele+1
    do while ( ele <= total )
      ele_info(ilevel)%str_neig(2,ele)= ele+1;ele_info(ilevel)%str_neig(3,ele)=ele-1;ele_info(ilevel)%str_neig(1,ele)=ele+current-1
      ele = ele+1

      ele_info(ilevel)%str_neig(2,ele)=ele-1;ele_info(ilevel)%str_neig(3,ele)=ele+1;ele_info(ilevel)%str_neig(1,ele)=ele-current -1
      ele = ele+1
    end do
    ele_info(ilevel)%str_neig(3,ele-1)=0 ! last ele in the row
    irow = irow -1
  end do

end subroutine get_str_neig_multigrid










!@brief> this subroutine gets number of splitting (n) and ele and returns
! ele row number and its position in the row
subroutine get_ups_info(n, ele, irow, ipos, sum_up, irow_ups)
  ! sum_up:: sum of up tris between row 1 and (irow-1)

  ! external vbls
  integer, intent(in):: n, ele
  ! local vbls
  integer, INTENT(INOUT):: irow, ipos, sum_up, irow_ups
  integer ::  ele2, row, i

  ele2 = ele
  i=1

    row = 2**n
    do while( ele2>row )
      ele2 = ele2 - row
      row = row - 1
      i = i+1
    end do
    irow = i
    ipos = ele2

    irow_ups = 2**n - (irow-1)
    ! this is arithmetic sum rule
    sum_up = ((irow-1) * ((2**(n+1)) + ((irow-1) -1)*(-1)))/2

  ! print*, 'up ele =',ele
  ! print*, irow, ipos,sum_up, irow_ups

end subroutine get_ups_info


subroutine get_downs_info(n, ele, irow, ipos, sum_down, irow_downs)
  ! sum_down:: sum of down tris between row 1 and (irow-1)

  ! external vbls
  integer, intent(in):: n, ele
  ! local vbls
  integer, INTENT(INOUT):: irow, ipos, sum_down, irow_downs
  integer ::  ele2, row, i

  ele2 = ele
  i=1

    row = 2**n
    do while( ele2>row )
      ele2 = ele2 - row
      row = row - 1
      i = i+1
    end do
    irow = i
    ipos = ele2

    irow_downs = 2**n - (irow-1)
    ! this is arithmetic sum rule
    sum_down = ((irow-1) * ((2**(n+1)) + ((irow-1) -1)*(-1)))/2

  ! print*, 'down ele=', ele
  ! print*, irow, ipos,sum_down, irow_downs

end subroutine get_downs_info


subroutine update_overlaps3(meshlist, tnew, told, t_bc, n_split, nface, totele_unst, nloc, str_neig)
  ! this subroutine updates overlaps of each unstructured ele
  ! local face numbers:
  !    |\
  !   3| \ 2
  !    |  \
  !    |___\
  !      1

  ! local node numbers
  !    2
  !    |\
  !    | \
  !    |  \
  !    |___\
  !    3   1
  !
  ! eternal vbls
  type(Mesh), dimension(:), INTENT(INOUT) :: meshList
  real, INTENT(IN) :: tnew(:,:,:), t_bc(:,:), told(:,:,:)
  integer, intent(in) :: n_split, nface, totele_unst, nloc, str_neig(:,:)

  ! local vbls
  integer :: str_ele, j, Nnodes(2), Npos, Nside, ele22, i_got_boundary, r_got_boundary
  integer :: iface, ipos, irow, u_iface, un_ele, orientation

  do un_ele=1,totele_unst
    do str_ele=1,4**n_split
      call get_str_info(n_split, str_ele, irow, ipos, orientation)
      do iface=1,nface
        ele22 = str_neig(iface,str_ele)
        i_got_boundary = (sign(1, -Npos) + 1 )/2
        r_got_boundary = real(i_got_boundary)
        if ( ele22==0 ) then

          if ( iface==1 ) then
            call getNeigDataMesh(meshlist, un_ele, 1, Npos, Nside, Nnodes)
            if ( Npos==0 ) then
              meshlist(un_ele)%t_overlap(:,1) = t_bc(:,1)
              meshlist(un_ele)%t_overlap_old(:,1) = t_bc(:,1)
            else
              if ( Nside==2 ) then
                if ( meshList(un_ele)%Dir(1) ) then
                  meshList(Npos)%t_overlap((2**n_split-(ipos/2+1) +1)*nloc-2:(2**n_split-(ipos/2+1) +1)*nloc,Nside)&
                  = tnew(:,str_ele, un_ele)
                  meshList(Npos)%t_overlap_old((2**n_split-(ipos/2+1) +1)*nloc-2:(2**n_split-(ipos/2+1) +1)*nloc,Nside)&
                  = told(:,str_ele, un_ele)
                else
                  meshList(Npos)%t_overlap((ipos/2+1)*nloc-2:(ipos/2+1)*nloc,Nside)= tnew(:,str_ele, un_ele)
                  meshList(Npos)%t_overlap_old((ipos/2+1)*nloc-2:(ipos/2+1)*nloc,Nside)= told(:,str_ele, un_ele)
                end if
              else
                if ( meshList(un_ele)%Dir(1) ) then
                  meshList(Npos)%t_overlap((ipos/2+1)*nloc-2:(ipos/2+1)*nloc,Nside) = tnew(:,str_ele, un_ele)
                  meshList(Npos)%t_overlap_old((ipos/2+1)*nloc-2:(ipos/2+1)*nloc,Nside) = told(:,str_ele, un_ele)
                else
                  meshList(Npos)%t_overlap((2**n_split-(ipos/2+1) +1)*nloc-2:(2**n_split-(ipos/2+1) +1)*nloc,Nside)&
                  = tnew(:,str_ele, un_ele)
                  meshList(Npos)%t_overlap_old((2**n_split-(ipos/2+1) +1)*nloc-2:(2**n_split-(ipos/2+1) +1)*nloc,Nside)&
                  = told(:,str_ele, un_ele)
                end if
              end if
            end if


          elseif ( iface==2 ) then
            call getNeigDataMesh(meshlist, un_ele, 3, Npos, Nside, Nnodes)
            if ( Npos==0 ) then
              meshlist(un_ele)%t_overlap(:,3) = t_bc(:,1)
              meshlist(un_ele)%t_overlap_old(:,3) = t_bc(:,1)
            else

              if ( Nside==2 ) then
                if ( meshList(un_ele)%Dir(3) ) then
                  meshList(Npos)%t_overlap((2**n_split-irow +1)*nloc-2,Nside) = tnew(1,str_ele, un_ele)
                  meshList(Npos)%t_overlap((2**n_split-irow +1)*nloc-1,Nside) = tnew(2,str_ele, un_ele)
                  meshList(Npos)%t_overlap((2**n_split-irow +1)*nloc  ,Nside) = tnew(3,str_ele, un_ele)

                  meshList(Npos)%t_overlap_old((2**n_split-irow +1)*nloc-2,Nside) = told(1,str_ele, un_ele)
                  meshList(Npos)%t_overlap_old((2**n_split-irow +1)*nloc-1,Nside) = told(2,str_ele, un_ele)
                  meshList(Npos)%t_overlap_old((2**n_split-irow +1)*nloc  ,Nside) = told(3,str_ele, un_ele)
                else
                  meshList(Npos)%t_overlap(irow*nloc-2,Nside) = tnew(1,str_ele, un_ele)
                  meshList(Npos)%t_overlap(irow*nloc-1,Nside) = tnew(2,str_ele, un_ele)
                  meshList(Npos)%t_overlap(irow*nloc  ,Nside) = tnew(3,str_ele, un_ele)

                  meshList(Npos)%t_overlap_old(irow*nloc-2,Nside) = told(1,str_ele, un_ele)
                  meshList(Npos)%t_overlap_old(irow*nloc-1,Nside) = told(2,str_ele, un_ele)
                  meshList(Npos)%t_overlap_old(irow*nloc  ,Nside) = told(3,str_ele, un_ele)
                end if
              else
                if ( meshList(un_ele)%Dir(3) ) then
                  meshList(Npos)%t_overlap(irow*nloc-2,Nside)= tnew(1,str_ele, un_ele)
                  meshList(Npos)%t_overlap(irow*nloc-1,Nside)= tnew(2,str_ele, un_ele)
                  meshList(Npos)%t_overlap(irow*nloc  ,Nside)= tnew(3,str_ele, un_ele)

                  meshList(Npos)%t_overlap_old(irow*nloc-2,Nside)= told(1,str_ele, un_ele)
                  meshList(Npos)%t_overlap_old(irow*nloc-1,Nside)= told(2,str_ele, un_ele)
                  meshList(Npos)%t_overlap_old(irow*nloc  ,Nside)= told(3,str_ele, un_ele)
                else
                  meshList(Npos)%t_overlap((2**n_split-irow +1)*nloc-2,Nside)= tnew(1,str_ele, un_ele)
                  meshList(Npos)%t_overlap((2**n_split-irow +1)*nloc-1,Nside)= tnew(2,str_ele, un_ele)
                  meshList(Npos)%t_overlap((2**n_split-irow +1)*nloc  ,Nside)= tnew(3,str_ele, un_ele)

                  meshList(Npos)%t_overlap_old((2**n_split-irow +1)*nloc-2,Nside)= told(1,str_ele, un_ele)
                  meshList(Npos)%t_overlap_old((2**n_split-irow +1)*nloc-1,Nside)= told(2,str_ele, un_ele)
                  meshList(Npos)%t_overlap_old((2**n_split-irow +1)*nloc  ,Nside)= told(3,str_ele, un_ele)
                end if
              end if
            end if


          else
            call getNeigDataMesh(meshlist, un_ele, 2, Npos, Nside, Nnodes)
            if ( Npos==0 ) then
              meshlist(un_ele)%t_overlap(:,2) = t_bc(:,1)
              meshlist(un_ele)%t_overlap_old(:,2) = t_bc(:,1)
            else

              if ( Nside==2 ) then
                if ( meshList(un_ele)%Dir(2) ) then
                  meshList(Npos)%t_overlap((irow-1)*nloc+1,Nside)= tnew(1,str_ele, un_ele)
                  meshList(Npos)%t_overlap((irow-1)*nloc+2,Nside)= tnew(2,str_ele, un_ele)
                  meshList(Npos)%t_overlap((irow-1)*nloc+3,Nside)= tnew(3,str_ele, un_ele)

                  meshList(Npos)%t_overlap_old((irow-1)*nloc+1,Nside)= told(1,str_ele, un_ele)
                  meshList(Npos)%t_overlap_old((irow-1)*nloc+2,Nside)= told(2,str_ele, un_ele)
                  meshList(Npos)%t_overlap_old((irow-1)*nloc+3,Nside)= told(3,str_ele, un_ele)
                else
                  meshList(Npos)%t_overlap((2**n_split-irow +1)*nloc-2,Nside)= tnew(1,str_ele, un_ele)
                  meshList(Npos)%t_overlap((2**n_split-irow +1)*nloc-1,Nside)= tnew(2,str_ele, un_ele)
                  meshList(Npos)%t_overlap((2**n_split-irow +1)*nloc  ,Nside)= tnew(3,str_ele, un_ele)

                  meshList(Npos)%t_overlap_old((2**n_split-irow +1)*nloc-2,Nside)= told(1,str_ele, un_ele)
                  meshList(Npos)%t_overlap_old((2**n_split-irow +1)*nloc-1,Nside)= told(2,str_ele, un_ele)
                  meshList(Npos)%t_overlap_old((2**n_split-irow +1)*nloc  ,Nside)= told(3,str_ele, un_ele)
                end if

              else
                if ( meshList(un_ele)%Dir(2) ) then
                  meshList(Npos)%t_overlap((2**n_split-irow +1)*nloc-2,Nside)= tnew(1,str_ele, un_ele)
                  meshList(Npos)%t_overlap((2**n_split-irow +1)*nloc-1,Nside)= tnew(2,str_ele, un_ele)
                  meshList(Npos)%t_overlap((2**n_split-irow +1)*nloc  ,Nside)= tnew(3,str_ele, un_ele)

                  meshList(Npos)%t_overlap_old((2**n_split-irow +1)*nloc-2,Nside)= told(1,str_ele, un_ele)
                  meshList(Npos)%t_overlap_old((2**n_split-irow +1)*nloc-1,Nside)= told(2,str_ele, un_ele)
                  meshList(Npos)%t_overlap_old((2**n_split-irow +1)*nloc  ,Nside)= told(3,str_ele, un_ele)
                else
                  meshList(Npos)%t_overlap(irow*nloc-2,Nside)= tnew(1,str_ele, un_ele)
                  meshList(Npos)%t_overlap(irow*nloc-1,Nside)= tnew(2,str_ele, un_ele)
                  meshList(Npos)%t_overlap(irow*nloc  ,Nside)= tnew(3,str_ele, un_ele)

                  meshList(Npos)%t_overlap_old(irow*nloc-2,Nside)= told(1,str_ele, un_ele)
                  meshList(Npos)%t_overlap_old(irow*nloc-1,Nside)= told(2,str_ele, un_ele)
                  meshList(Npos)%t_overlap_old(irow*nloc  ,Nside)= told(3,str_ele, un_ele)
                end if
              end if
            end if


          end if
        end if ! ele22
      end do ! iface
    end do ! str_ele
  end do ! un_ele


end subroutine update_overlaps3


subroutine update_overlaps2(meshlist,surf_ele, tnew, told, t_bc, n_split, nface, totele_unst, nloc, str_neig,x_loc)
  ! this subroutine updates overlaps of each unstructured ele
  ! local face numbers:
  !    |\
  !   3| \ 2
  !    |  \
  !    |___\
  !      1

  ! local node numbers
  !    2
  !    |\
  !    | \
  !    |  \
  !    |___\
  !    3   1
  !
  ! eternal vbls
  type(Mesh), dimension(:), INTENT(INOUT) :: meshList
  real, INTENT(IN) :: tnew(:,:,:), told(:,:,:)
  integer, intent(in) :: n_split, nface, totele_unst, nloc, str_neig(:,:), surf_ele(:,:)
  real, intent(inout) ::  t_bc(:), x_loc(:,:)

  ! local vbls
  integer :: str_ele, j, Nnodes(2), Npos, Nside, ele22
  integer :: iface, ipos, irow, u_iface, un_ele, orientation,i

  do un_ele=1,totele_unst
    do i=1,2**n_split
      str_ele = surf_ele(i,1) ! for iface==1
      call get_str_info(n_split, str_ele, irow, ipos, orientation)
      Npos = meshList(un_ele)%Neig(1)
      call get_splitting(meshList(un_ele)%X, n_split, str_ele, x_loc)
      if ( Npos==0 ) then
      t_bc(1) = boundary(x_loc(1,1) , x_loc(2,1))
      t_bc(2) = boundary(x_loc(1,3) , x_loc(2,3))
        meshlist(un_ele)%t_overlap((ipos/2)*nloc+1,1) = t_bc(1)
        meshlist(un_ele)%t_overlap((ipos/2)*nloc+3,1) = t_bc(2)

        meshlist(un_ele)%t_overlap_old((ipos/2)*nloc+1,1) = t_bc(1)
        meshlist(un_ele)%t_overlap_old((ipos/2)*nloc+3,1) = t_bc(2)

      else
        Nside = meshList(un_ele)%fNeig(1)
        if ( Nside==2 ) then
          if ( meshList(un_ele)%Dir(1) ) then
            meshList(Npos)%t_overlap((2**n_split-(ipos/2+1) +1)*nloc-2:(2**n_split-(ipos/2+1) +1)*nloc,Nside)&
            = tnew(:,str_ele, un_ele)
            meshList(Npos)%t_overlap_old((2**n_split-(ipos/2+1) +1)*nloc-2:(2**n_split-(ipos/2+1) +1)*nloc,Nside)&
            = told(:,str_ele, un_ele)
          else
            meshList(Npos)%t_overlap((ipos/2+1)*nloc-2:(ipos/2+1)*nloc,Nside)= tnew(:,str_ele, un_ele)
            meshList(Npos)%t_overlap_old((ipos/2+1)*nloc-2:(ipos/2+1)*nloc,Nside)= told(:,str_ele, un_ele)
          end if
        else
          if ( meshList(un_ele)%Dir(1) ) then
            meshList(Npos)%t_overlap((ipos/2+1)*nloc-2:(ipos/2+1)*nloc,Nside) = tnew(:,str_ele, un_ele)
            meshList(Npos)%t_overlap_old((ipos/2+1)*nloc-2:(ipos/2+1)*nloc,Nside) = told(:,str_ele, un_ele)
          else
            meshList(Npos)%t_overlap((2**n_split-(ipos/2+1) +1)*nloc-2:(2**n_split-(ipos/2+1) +1)*nloc,Nside)&
            = tnew(:,str_ele, un_ele)
            meshList(Npos)%t_overlap_old((2**n_split-(ipos/2+1) +1)*nloc-2:(2**n_split-(ipos/2+1) +1)*nloc,Nside)&
            = told(:,str_ele, un_ele)
          end if
        end if
      end if
    end do ! str_ele fir face 1


    do i=1,2**n_split ! iface==3
      str_ele = surf_ele(i,3)
      call get_str_info(n_split, str_ele, irow, ipos, orientation)
      Npos = meshList(un_ele)%Neig(3)
      call get_splitting(meshList(un_ele)%X, n_split, str_ele, x_loc)
      if ( Npos==0 ) then
            t_bc(1) = boundary(x_loc(1,2) , x_loc(2,2))
            t_bc(2) = boundary(x_loc(1,3) , x_loc(2,3))
            meshlist(un_ele)%t_overlap((irow-1)*nloc+2,3) = t_bc(1)
            meshlist(un_ele)%t_overlap((irow-1)*nloc+3,3) = t_bc(2)

            meshlist(un_ele)%t_overlap_old((irow-1)*nloc+2,3) = t_bc(1)
            meshlist(un_ele)%t_overlap_old((irow-1)*nloc+3,3) = t_bc(2)

      else
        Nside = meshList(un_ele)%fNeig(3)
        if ( Nside==2 ) then
          if ( meshList(un_ele)%Dir(3) ) then
            meshList(Npos)%t_overlap((2**n_split-irow +1)*nloc-2,Nside) = tnew(1,str_ele, un_ele)
            meshList(Npos)%t_overlap((2**n_split-irow +1)*nloc-1,Nside) = tnew(2,str_ele, un_ele)
            meshList(Npos)%t_overlap((2**n_split-irow +1)*nloc  ,Nside) = tnew(3,str_ele, un_ele)

            meshList(Npos)%t_overlap_old((2**n_split-irow +1)*nloc-2,Nside) = told(1,str_ele, un_ele)
            meshList(Npos)%t_overlap_old((2**n_split-irow +1)*nloc-1,Nside) = told(2,str_ele, un_ele)
            meshList(Npos)%t_overlap_old((2**n_split-irow +1)*nloc  ,Nside) = told(3,str_ele, un_ele)
          else
            meshList(Npos)%t_overlap(irow*nloc-2,Nside) = tnew(1,str_ele, un_ele)
            meshList(Npos)%t_overlap(irow*nloc-1,Nside) = tnew(2,str_ele, un_ele)
            meshList(Npos)%t_overlap(irow*nloc  ,Nside) = tnew(3,str_ele, un_ele)

            meshList(Npos)%t_overlap_old(irow*nloc-2,Nside) = told(1,str_ele, un_ele)
            meshList(Npos)%t_overlap_old(irow*nloc-1,Nside) = told(2,str_ele, un_ele)
            meshList(Npos)%t_overlap_old(irow*nloc  ,Nside) = told(3,str_ele, un_ele)
          end if
        else
          if ( meshList(un_ele)%Dir(3) ) then
            meshList(Npos)%t_overlap(irow*nloc-2,Nside)= tnew(1,str_ele, un_ele)
            meshList(Npos)%t_overlap(irow*nloc-1,Nside)= tnew(2,str_ele, un_ele)
            meshList(Npos)%t_overlap(irow*nloc  ,Nside)= tnew(3,str_ele, un_ele)

            meshList(Npos)%t_overlap_old(irow*nloc-2,Nside)= told(1,str_ele, un_ele)
            meshList(Npos)%t_overlap_old(irow*nloc-1,Nside)= told(2,str_ele, un_ele)
            meshList(Npos)%t_overlap_old(irow*nloc  ,Nside)= told(3,str_ele, un_ele)
          else
            meshList(Npos)%t_overlap((2**n_split-irow +1)*nloc-2,Nside)= tnew(1,str_ele, un_ele)
            meshList(Npos)%t_overlap((2**n_split-irow +1)*nloc-1,Nside)= tnew(2,str_ele, un_ele)
            meshList(Npos)%t_overlap((2**n_split-irow +1)*nloc  ,Nside)= tnew(3,str_ele, un_ele)

            meshList(Npos)%t_overlap_old((2**n_split-irow +1)*nloc-2,Nside)= told(1,str_ele, un_ele)
            meshList(Npos)%t_overlap_old((2**n_split-irow +1)*nloc-1,Nside)= told(2,str_ele, un_ele)
            meshList(Npos)%t_overlap_old((2**n_split-irow +1)*nloc  ,Nside)= told(3,str_ele, un_ele)
          end if
        end if
      end if
    end do


    do i=1,2**n_split
      str_ele = surf_ele(i,2)
      call get_str_info(n_split, str_ele, irow, ipos, orientation)
      npos = meshList(un_ele)%Neig(2)
      call get_splitting(meshList(un_ele)%X, n_split, str_ele, x_loc)
      if ( Npos==0 ) then
        t_bc(1) = boundary(x_loc(1,1) , x_loc(2,1))
        t_bc(2) = boundary(x_loc(1,2) , x_loc(2,2))
        meshlist(un_ele)%t_overlap((irow-1)*nloc+1,2) = t_bc(1)
        meshlist(un_ele)%t_overlap((irow-1)*nloc+2,2) = t_bc(2)

        meshlist(un_ele)%t_overlap_old((irow-1)*nloc+1,2) = t_bc(1)
        meshlist(un_ele)%t_overlap_old((irow-1)*nloc+2,2) = t_bc(2)

      else
        Nside = meshList(un_ele)%fNeig(2)
        if ( Nside==2 ) then
          if (  meshList(un_ele)%Dir(2) ) then
            meshList(Npos)%t_overlap((irow-1)*nloc+1,Nside)= tnew(1,str_ele, un_ele)
            meshList(Npos)%t_overlap((irow-1)*nloc+2,Nside)= tnew(2,str_ele, un_ele)
            meshList(Npos)%t_overlap((irow-1)*nloc+3,Nside)= tnew(3,str_ele, un_ele)

            meshList(Npos)%t_overlap_old((irow-1)*nloc+1,Nside)= told(1,str_ele, un_ele)
            meshList(Npos)%t_overlap_old((irow-1)*nloc+2,Nside)= told(2,str_ele, un_ele)
            meshList(Npos)%t_overlap_old((irow-1)*nloc+3,Nside)= told(3,str_ele, un_ele)
          else
            meshList(Npos)%t_overlap((2**n_split-irow +1)*nloc-2,Nside)= tnew(1,str_ele, un_ele)
            meshList(Npos)%t_overlap((2**n_split-irow +1)*nloc-1,Nside)= tnew(2,str_ele, un_ele)
            meshList(Npos)%t_overlap((2**n_split-irow +1)*nloc  ,Nside)= tnew(3,str_ele, un_ele)

            meshList(Npos)%t_overlap_old((2**n_split-irow +1)*nloc-2,Nside)= told(1,str_ele, un_ele)
            meshList(Npos)%t_overlap_old((2**n_split-irow +1)*nloc-1,Nside)= told(2,str_ele, un_ele)
            meshList(Npos)%t_overlap_old((2**n_split-irow +1)*nloc  ,Nside)= told(3,str_ele, un_ele)
          end if

        else
          if ( meshList(un_ele)%Dir(2) ) then
            meshList(Npos)%t_overlap((2**n_split-irow +1)*nloc-2,Nside)= tnew(1,str_ele, un_ele)
            meshList(Npos)%t_overlap((2**n_split-irow +1)*nloc-1,Nside)= tnew(2,str_ele, un_ele)
            meshList(Npos)%t_overlap((2**n_split-irow +1)*nloc  ,Nside)= tnew(3,str_ele, un_ele)

            meshList(Npos)%t_overlap_old((2**n_split-irow +1)*nloc-2,Nside)= told(1,str_ele, un_ele)
            meshList(Npos)%t_overlap_old((2**n_split-irow +1)*nloc-1,Nside)= told(2,str_ele, un_ele)
            meshList(Npos)%t_overlap_old((2**n_split-irow +1)*nloc  ,Nside)= told(3,str_ele, un_ele)
          else
            meshList(Npos)%t_overlap(irow*nloc-2,Nside)= tnew(1,str_ele, un_ele)
            meshList(Npos)%t_overlap(irow*nloc-1,Nside)= tnew(2,str_ele, un_ele)
            meshList(Npos)%t_overlap(irow*nloc  ,Nside)= tnew(3,str_ele, un_ele)

            meshList(Npos)%t_overlap_old(irow*nloc-2,Nside)= told(1,str_ele, un_ele)
            meshList(Npos)%t_overlap_old(irow*nloc-1,Nside)= told(2,str_ele, un_ele)
            meshList(Npos)%t_overlap_old(irow*nloc  ,Nside)= told(3,str_ele, un_ele)
          end if
        end if
      end if
    end do ! str_ele
  end do ! un_ele


end subroutine update_overlaps2



subroutine update_overlaps(meshlist,surf_ele, tnew, told, t_bc, i_split, nface, totele_unst, nloc, str_neig)
  ! this subroutine updates overlaps of each unstructured ele
  ! local face numbers:
  !    |\
  !   3| \ 2
  !    |  \
  !    |___\
  !      1

  ! local node numbers
  !    2
  !    |\
  !    | \
  !    |  \
  !    |___\
  !    3   1
  !
  ! eternal vbls
  type(Mesh), dimension(:), INTENT(INOUT) :: meshList
  real, INTENT(IN) :: tnew(:,:,:), told(:,:,:)
  integer, intent(in) :: i_split, nface, totele_unst, nloc, str_neig(:,:), surf_ele(:,:)
  real, intent(inout) ::  t_bc(:)

  ! local vbls
  real :: x_loc(2,3)
  integer :: str_ele, j, Nnodes(2), Npos, Nside, ele22
  integer :: iface, ipos, irow, u_iface, un_ele, orientation,i

  do un_ele=1,totele_unst
    do i=1,2**i_split
      str_ele = surf_ele(i,1) ! for iface==1
! print*, un_ele,i_split,  i, surf_ele(i,1)
      call get_str_info(i_split, str_ele, irow, ipos, orientation)
      Npos = meshList(un_ele)%Neig(1)
      call get_splitting(meshList(un_ele)%X, i_split, str_ele, x_loc)
      if ( Npos==0 ) then
      t_bc(1) = boundary(x_loc(1,1) , x_loc(2,1))
      t_bc(2) = boundary(x_loc(1,3) , x_loc(2,3))
        meshlist(un_ele)%t_overlap((ipos/2)*nloc+1,1) = t_bc(1)
        meshlist(un_ele)%t_overlap((ipos/2)*nloc+3,1) = t_bc(2)

        meshlist(un_ele)%t_overlap_old((ipos/2)*nloc+1,1) = t_bc(1)
        meshlist(un_ele)%t_overlap_old((ipos/2)*nloc+3,1) = t_bc(2)

      else
        Nside = meshList(un_ele)%fNeig(1)
        if ( Nside==2 ) then
          if ( meshList(un_ele)%Dir(1) ) then
            meshList(Npos)%t_overlap((2**i_split-(ipos/2+1) +1)*nloc-2:(2**i_split-(ipos/2+1) +1)*nloc,Nside)&
            = tnew(:,str_ele, un_ele)
            meshList(Npos)%t_overlap_old((2**i_split-(ipos/2+1) +1)*nloc-2:(2**i_split-(ipos/2+1) +1)*nloc,Nside)&
            = told(:,str_ele, un_ele)
          else
            meshList(Npos)%t_overlap((ipos/2+1)*nloc-2:(ipos/2+1)*nloc,Nside)= tnew(:,str_ele, un_ele)
            meshList(Npos)%t_overlap_old((ipos/2+1)*nloc-2:(ipos/2+1)*nloc,Nside)= told(:,str_ele, un_ele)
          end if
        else
          if ( meshList(un_ele)%Dir(1) ) then
            meshList(Npos)%t_overlap((ipos/2+1)*nloc-2:(ipos/2+1)*nloc,Nside) = tnew(:,str_ele, un_ele)
            meshList(Npos)%t_overlap_old((ipos/2+1)*nloc-2:(ipos/2+1)*nloc,Nside) = told(:,str_ele, un_ele)
          else
            meshList(Npos)%t_overlap((2**i_split-(ipos/2+1) +1)*nloc-2:(2**i_split-(ipos/2+1) +1)*nloc,Nside)&
            = tnew(:,str_ele, un_ele)
            meshList(Npos)%t_overlap_old((2**i_split-(ipos/2+1) +1)*nloc-2:(2**i_split-(ipos/2+1) +1)*nloc,Nside)&
            = told(:,str_ele, un_ele)
          end if
        end if
      end if
    end do ! str_ele fir face 1


    do i=1,2**i_split ! iface==3
      str_ele = surf_ele(i,3)
      call get_str_info(i_split, str_ele, irow, ipos, orientation)
      Npos = meshList(un_ele)%Neig(3)
      call get_splitting(meshList(un_ele)%X, i_split, str_ele, x_loc)
      if ( Npos==0 ) then
            t_bc(1) = boundary(x_loc(1,2) , x_loc(2,2))
            t_bc(2) = boundary(x_loc(1,3) , x_loc(2,3))
            meshlist(un_ele)%t_overlap((irow-1)*nloc+2,3) = t_bc(1)
            meshlist(un_ele)%t_overlap((irow-1)*nloc+3,3) = t_bc(2)

            meshlist(un_ele)%t_overlap_old((irow-1)*nloc+2,3) = t_bc(1)
            meshlist(un_ele)%t_overlap_old((irow-1)*nloc+3,3) = t_bc(2)

      else
        Nside = meshList(un_ele)%fNeig(3)
        if ( Nside==2 ) then
          if ( meshList(un_ele)%Dir(3) ) then
            meshList(Npos)%t_overlap((2**i_split-irow +1)*nloc-2,Nside) = tnew(1,str_ele, un_ele)
            meshList(Npos)%t_overlap((2**i_split-irow +1)*nloc-1,Nside) = tnew(2,str_ele, un_ele)
            meshList(Npos)%t_overlap((2**i_split-irow +1)*nloc  ,Nside) = tnew(3,str_ele, un_ele)

            meshList(Npos)%t_overlap_old((2**i_split-irow +1)*nloc-2,Nside) = told(1,str_ele, un_ele)
            meshList(Npos)%t_overlap_old((2**i_split-irow +1)*nloc-1,Nside) = told(2,str_ele, un_ele)
            meshList(Npos)%t_overlap_old((2**i_split-irow +1)*nloc  ,Nside) = told(3,str_ele, un_ele)
          else
            meshList(Npos)%t_overlap(irow*nloc-2,Nside) = tnew(1,str_ele, un_ele)
            meshList(Npos)%t_overlap(irow*nloc-1,Nside) = tnew(2,str_ele, un_ele)
            meshList(Npos)%t_overlap(irow*nloc  ,Nside) = tnew(3,str_ele, un_ele)

            meshList(Npos)%t_overlap_old(irow*nloc-2,Nside) = told(1,str_ele, un_ele)
            meshList(Npos)%t_overlap_old(irow*nloc-1,Nside) = told(2,str_ele, un_ele)
            meshList(Npos)%t_overlap_old(irow*nloc  ,Nside) = told(3,str_ele, un_ele)
          end if
        else
          if ( meshList(un_ele)%Dir(3) ) then
            meshList(Npos)%t_overlap(irow*nloc-2,Nside)= tnew(1,str_ele, un_ele)
            meshList(Npos)%t_overlap(irow*nloc-1,Nside)= tnew(2,str_ele, un_ele)
            meshList(Npos)%t_overlap(irow*nloc  ,Nside)= tnew(3,str_ele, un_ele)

            meshList(Npos)%t_overlap_old(irow*nloc-2,Nside)= told(1,str_ele, un_ele)
            meshList(Npos)%t_overlap_old(irow*nloc-1,Nside)= told(2,str_ele, un_ele)
            meshList(Npos)%t_overlap_old(irow*nloc  ,Nside)= told(3,str_ele, un_ele)
          else
            meshList(Npos)%t_overlap((2**i_split-irow +1)*nloc-2,Nside)= tnew(1,str_ele, un_ele)
            meshList(Npos)%t_overlap((2**i_split-irow +1)*nloc-1,Nside)= tnew(2,str_ele, un_ele)
            meshList(Npos)%t_overlap((2**i_split-irow +1)*nloc  ,Nside)= tnew(3,str_ele, un_ele)

            meshList(Npos)%t_overlap_old((2**i_split-irow +1)*nloc-2,Nside)= told(1,str_ele, un_ele)
            meshList(Npos)%t_overlap_old((2**i_split-irow +1)*nloc-1,Nside)= told(2,str_ele, un_ele)
            meshList(Npos)%t_overlap_old((2**i_split-irow +1)*nloc  ,Nside)= told(3,str_ele, un_ele)
          end if
        end if
      end if
    end do


    do i=1,2**i_split
      str_ele = surf_ele(i,2)
      call get_str_info(i_split, str_ele, irow, ipos, orientation)
      npos = meshList(un_ele)%Neig(2)
      call get_splitting(meshList(un_ele)%X, i_split, str_ele, x_loc)
      if ( Npos==0 ) then
        t_bc(1) = boundary(x_loc(1,1) , x_loc(2,1))
        t_bc(2) = boundary(x_loc(1,2) , x_loc(2,2))
        meshlist(un_ele)%t_overlap((irow-1)*nloc+1,2) = t_bc(1)
        meshlist(un_ele)%t_overlap((irow-1)*nloc+2,2) = t_bc(2)

        meshlist(un_ele)%t_overlap_old((irow-1)*nloc+1,2) = t_bc(1)
        meshlist(un_ele)%t_overlap_old((irow-1)*nloc+2,2) = t_bc(2)

      else
        Nside = meshList(un_ele)%fNeig(2)
        if ( Nside==2 ) then
          if (  meshList(un_ele)%Dir(2) ) then
            meshList(Npos)%t_overlap((irow-1)*nloc+1,Nside)= tnew(1,str_ele, un_ele)
            meshList(Npos)%t_overlap((irow-1)*nloc+2,Nside)= tnew(2,str_ele, un_ele)
            meshList(Npos)%t_overlap((irow-1)*nloc+3,Nside)= tnew(3,str_ele, un_ele)

            meshList(Npos)%t_overlap_old((irow-1)*nloc+1,Nside)= told(1,str_ele, un_ele)
            meshList(Npos)%t_overlap_old((irow-1)*nloc+2,Nside)= told(2,str_ele, un_ele)
            meshList(Npos)%t_overlap_old((irow-1)*nloc+3,Nside)= told(3,str_ele, un_ele)
          else
            meshList(Npos)%t_overlap((2**i_split-irow +1)*nloc-2,Nside)= tnew(1,str_ele, un_ele)
            meshList(Npos)%t_overlap((2**i_split-irow +1)*nloc-1,Nside)= tnew(2,str_ele, un_ele)
            meshList(Npos)%t_overlap((2**i_split-irow +1)*nloc  ,Nside)= tnew(3,str_ele, un_ele)

            meshList(Npos)%t_overlap_old((2**i_split-irow +1)*nloc-2,Nside)= told(1,str_ele, un_ele)
            meshList(Npos)%t_overlap_old((2**i_split-irow +1)*nloc-1,Nside)= told(2,str_ele, un_ele)
            meshList(Npos)%t_overlap_old((2**i_split-irow +1)*nloc  ,Nside)= told(3,str_ele, un_ele)
          end if

        else
          if ( meshList(un_ele)%Dir(2) ) then
            meshList(Npos)%t_overlap((2**i_split-irow +1)*nloc-2,Nside)= tnew(1,str_ele, un_ele)
            meshList(Npos)%t_overlap((2**i_split-irow +1)*nloc-1,Nside)= tnew(2,str_ele, un_ele)
            meshList(Npos)%t_overlap((2**i_split-irow +1)*nloc  ,Nside)= tnew(3,str_ele, un_ele)

            meshList(Npos)%t_overlap_old((2**i_split-irow +1)*nloc-2,Nside)= told(1,str_ele, un_ele)
            meshList(Npos)%t_overlap_old((2**i_split-irow +1)*nloc-1,Nside)= told(2,str_ele, un_ele)
            meshList(Npos)%t_overlap_old((2**i_split-irow +1)*nloc  ,Nside)= told(3,str_ele, un_ele)
          else
            meshList(Npos)%t_overlap(irow*nloc-2,Nside)= tnew(1,str_ele, un_ele)
            meshList(Npos)%t_overlap(irow*nloc-1,Nside)= tnew(2,str_ele, un_ele)
            meshList(Npos)%t_overlap(irow*nloc  ,Nside)= tnew(3,str_ele, un_ele)

            meshList(Npos)%t_overlap_old(irow*nloc-2,Nside)= told(1,str_ele, un_ele)
            meshList(Npos)%t_overlap_old(irow*nloc-1,Nside)= told(2,str_ele, un_ele)
            meshList(Npos)%t_overlap_old(irow*nloc  ,Nside)= told(3,str_ele, un_ele)
          end if
        end if
      end if
    end do ! str_ele
  end do ! un_ele


end subroutine update_overlaps


! this function defines the equation for boundary conditions
real function boundary(a,b) result(bc)
  implicit none
  real, intent(in) :: a,b
  bc = sin(a+b)
end function boundary


end module splitting
