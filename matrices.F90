module matrices
  use structured_meshgen
  use Structures
  use Msh2Tri
  use splitting
implicit none

contains


  ! !brief:: this subroutine generate mass and stiff stencil matrices
  ! subroutine gen_stcl(mass_stcl, stiff_stcl, diff_stcl, n,nx,k,detwei,nloc,ngi,ndim,dt)
  !   implicit none
  !   !global vbl
  !   integer , intent(in) :: nloc, ngi, ndim
  !   REAL, intent( in ) :: N(ngi,nloc), nx(ngi,ndim,nloc),detwei(ngi),dt,k
  !   real, dimension(ngi, ndim, nloc, nloc), intent(inout) :: stiff_stcl, diff_stcl
  !   real, dimension(ngi,nloc,nloc), intent(inout) :: mass_stcl
  !
  !   !local vbl
  !   integer :: idim, iloc, jloc, gi
  !
  !   do jloc=1,nloc
  !     do iloc=1,nloc
  !       do gi=1,ngi
  !         mass_stcl(gi, iloc, jloc) = n(gi,iloc) * n(gi,jloc)/dt
  !         do idim=1,ndim
  !           stiff_stcl(gi,idim, iloc, jloc) =  nx(gi,idim,iloc)*n(gi,jloc)*detwei(gi)
  !           diff_stcl(gi,idim,iloc, jloc) = k * nx(gi,idim,iloc)* nx(gi,idim,jloc) * detwei(gi)
  !         end do
  !       end do
  !     end do
  !   end do
  !
  ! end subroutine gen_stcl


  !@brief> this subroutine calculates the volume intergral of the diffusion term after integration by parts
  subroutine add_diffusion_vol(k, nx, ngi, ndim, nloc, detwei, tnew_xgi, diff_vol)
    implicit none
    !external vbls
    integer, intent(in) :: ngi, nloc, ndim
    real, intent(in) :: k
    real, dimension( ngi ), intent(in) :: detwei
    real, dimension( ngi, ndim ), intent(in) :: tnew_xgi
    real, dimension( ngi, ndim, nloc ), intent(in) :: nx
    !local vbls
    integer :: idim, iloc, gi
    real, intent(inout) :: diff_vol(:)

    diff_vol(:)=0.0
    do iloc=1,nloc
      do idim=1,ndim
        do gi=1,ngi
          diff_vol(iloc) = diff_vol(iloc) + k * nx(gi,idim,iloc) * tnew_xgi(gi,idim) * detwei(gi)
        end do
      end do ! quadrature
    end do ! iloc
! print*, diff_vol
  end subroutine add_diffusion_vol


  !@brief> this subroutine calculates the surface intergral of the diffusion term after integration by parts
  ! node1 and node2 are corner nodes of un_ele located on mface
  ! mf_center is center of mface of un_ele
  subroutine add_diffusion_surf(meshlist, un_ele, mface, iface, nface, k, n_split, n, ngi, sngi, ele22,sn2,&
                                nloc, sdetwei, tnew_sgi,tnew_sgi2, diff_surf, delta_x,sn,my_diff_surf,neig_diff_surf)
    implicit none
    !external vbls
    type(Mesh), dimension(:), INTENT(IN) :: meshList
    integer, intent(in) :: ngi, sngi,nloc, ele22, mface,n_split,un_ele,nface,iface
    real, intent(in) :: k
    real, dimension( sngi ), intent(in) :: sdetwei
    real, dimension( sngi ), intent(in) :: tnew_sgi,tnew_sgi2
    real, dimension( ngi, nloc ), intent(in) :: n
    real, intent(in) :: sn(sngi,nloc)
    real,intent(inout) :: my_diff_surf(nloc,nloc,nface),neig_diff_surf(nloc,nloc,nface),sn2(:,:)
    !local vbls
    real, intent(out) :: delta_x
    integer :: idim, iloc, gi,i, Npos, node1, node2, sgi,jloc
    real :: mf_center(2) !(x,y)
    real, intent(inout) :: diff_surf(:)

    if ( iface==1 ) then
      i=1
      node1=1
      node2=3
    elseif ( iface==2 ) then
      i=3
      node1=1
      node2=2
    else
      i=2
      node1=2
      node2=3
    end if
    Npos = meshList(un_ele)%Neig(i)

    ! ################### calculate delta x #########################################
    if ( ele22 /= 0 ) then
      delta_x = meshList(un_ele)%dc_str_ele(iface)
    elseif ( ele22==0 .and. Npos /=0) then
      delta_x = meshList(un_ele)%dc_unele(mface)/(2**n_split)
    elseif ( ele22==0 .and. Npos ==0) then
      mf_center(1) = (meshList(un_ele)%X(1,node1) + meshList(un_ele)%X(1,node2))/2
      mf_center(2) = (meshList(un_ele)%X(2,node1) + meshList(un_ele)%X(2,node2))/2

      delta_x = (meshList(un_ele)%center(1) - mf_center(1))**2 + (meshList(un_ele)%center(2) - mf_center(2))**2
      delta_x = sqrt(delta_x) / (2**n_split)
    end if

    ! ################ ed dx calculation ############################################
    ! do iloc=1,nloc
    !     diff_surf(iloc) = diff_surf(iloc) + (k/delta_x)*sum(sn(:,iloc) *  (tnew_sgi(:) - tnew_sgi2(:)) * sdetwei(:))
    ! end do ! iloc

  end subroutine add_diffusion_surf



  ! brief:: this subroutine calculates area of a semi_structured triangle within un_ele
  ! based on 3 corner nodes! and saves it in meshlist(un_ele)%area
  subroutine get_area(meshlist,un_ele, n_split)
    implicit none
    !global vbl
    type(Mesh), dimension(:), INTENT(INout) :: meshList
    integer, intent(in) :: un_ele, n_split
    ! local vbl
    real :: area
    area = meshList(un_ele)%X(1,1)*(meshList(un_ele)%X(2,2)-meshList(un_ele)%X(2,3)) - &
           meshList(un_ele)%X(2,1)*(meshList(un_ele)%X(1,2)-meshList(un_ele)%X(1,3)) + &
           meshList(un_ele)%X(1,2)*meshList(un_ele)%x(2,3) - meshList(un_ele)%x(2,2)*meshList(un_ele)%x(1,3)
    meshList(un_ele)%str_area = 0.5 * abs(area)/(4**n_split)
  end subroutine get_area




!brief>:: it converts a CSR matrix (sparse_matrix) into a normal matrix (global_matrix)
subroutine gen_global_matrix(totele, nloc, sparse_matrix, glob_matrix)
  implicit none
  ! external vbl
  type(sparse), intent(in) :: sparse_matrix
  integer :: totele, nloc
  !local vbl
  real, intent(inout) :: glob_matrix(:,:)
  integer :: L,i ,j, counter,k

  ! if (allocated(glob_matrix) ) deallocate(glob_matrix)
  ! allocate(glob_matrix(totele*nloc, totele*nloc))
  glob_matrix(:,:) = 0.0
    counter=0
    do i=1,size(sparse_matrix%g_iloc)
      if ( i<size(sparse_matrix%g_iloc) ) then
        L=sparse_matrix%g_iloc(i+1) - sparse_matrix%g_iloc(i)
        do j=1,L
          counter = counter +1
          glob_matrix(i, sparse_matrix%g_jloc(counter)) = sparse_matrix%val(counter)
        end do
      elseif (i==size(sparse_matrix%g_iloc)) then
        counter=counter+1
        do k=counter, size(sparse_matrix%val)
          glob_matrix(i, sparse_matrix%g_jloc(k)) = sparse_matrix%val(k)
        end do
      end if
    end do
end subroutine gen_global_matrix


!@brief> this subroutine takes CSR matrix and times it to told to make rhs (mass_matrix*told)
! (n*n) times (n*1) = (n*1)
subroutine csr_mul_array(sparse_matrix, array, result)
  implicit none
  ! external vbls
  type(sparse), intent(in) :: sparse_matrix
  real, allocatable, dimension(:), intent(in) :: array
  !local vbl
  real, allocatable, intent(inout) :: result(:)
  integer :: c1,c2,n, nloc, totele, ele

  nloc = size(sparse_matrix%g_jloc)/size(sparse_matrix%g_iloc)
  totele = size(sparse_matrix%g_jloc)/nloc
  c1=1
  c2=1
  result = 0.0
  do ele=1,size(sparse_matrix%g_iloc)
    do n=1,3
            result(c1) = result(c1) + sparse_matrix%val(c2) * array(sparse_matrix%g_jloc(c2))
            c2 = c2+1
    end do
    c1 = c1+1
  end do
end subroutine csr_mul_array



! !@brief> this subroutine takes CSR lhs and adds it to CSR lhs and saves it back to CSR lhs
! ! here because size of sparse_flux is largerthe lhs, everything will be saved in sparse_flux instead of sparse_lhs
! subroutine flux_plus_lhs(sparse_lhs, sparse_flux)
! implicit none
! ! external vbls
! type(sparse), intent(in) :: sparse_lhs
! type(sparse), intent(inout) :: sparse_flux
! !local vbl
! integer :: i,n
!
! do i=1,size(sparse_lhs%g_iloc)
!   i_f= sparse_flux%g_iloc(i)
!   i_lhs = sparse_lhs%i
!   result(i) = sum(sparse_matrix%val(n:n+2) * array(i))
! end do
! end subroutine flux_plus_lhs



!@brief> :: it gets told(nloc,totele) as its matrix and convertes it into an array
! of global_array(nloc*totele)
subroutine told_to_array(totele, nloc, matrix, global_array)
  implicit none
  !external vbl
  integer, intent(in) :: totele, nloc
  real, intent(in) :: matrix(:,:)
  !local vbl
  integer :: ele, iloc, counter
  real, allocatable, intent(out) :: global_array(:)
  allocate(global_array(totele*nloc))

  counter=1
  do ele=1,totele
    do iloc=1,nloc
      global_array(counter) = matrix(iloc,ele)
      counter = counter+1
    end do
  end do
end subroutine told_to_array


!@brief> :: it gets told(nloc,totele_str,totele_unst) as its matrix and convertes it into an array
! of global_array(nloc*totele)
subroutine told_to_array_semi(totele_str, totele_unst, nloc, matrix, global_array)
  implicit none
  !external vbl
  integer, intent(in) :: totele_unst, totele_str, nloc
  real, intent(in) :: matrix(:,:,:)
  !local vbl
  integer :: str, unstr, iloc, counter
  real, allocatable, intent(out) :: global_array(:)
  allocate(global_array(totele_unst*totele_str*nloc))

  counter=1
  do unstr=1,totele_unst
    do str=1,totele_str
      do iloc=1,nloc
        global_array(counter) = matrix(iloc,str, unstr)
        counter = counter+1
      end do
    end do
  end do
end subroutine told_to_array_semi


!"brief>:: it's opposite to told_to_array. It receives an array and turns it into (nloc,totele)
subroutine array_to_told(totele, nloc, array, told_style)
  implicit none
  !external vbl
  integer, intent(in) :: totele, nloc
  real, intent(in) :: array(:)
  !local vbl
  integer :: ele, iloc, counter
  real, allocatable, intent(out) :: told_style(:,:)
  allocate(told_style(nloc, totele))

  counter=1
  do ele=1,totele
    do iloc=1,nloc
      told_style(iloc,ele)=array(counter)
      counter= counter+1
    end do
  end do
end subroutine array_to_told



!"brief>:: it's opposite to told_to_array_semi. It receives an array and turns
! it into told(nloc,totele_str,totele_unst)
subroutine array_to_told_semi(totele_str,totele_unst, nloc, array, told_style)
  implicit none
  !external vbl
  integer, intent(in) :: totele_str,totele_unst, nloc
  real, intent(in) :: array(:)
  !local vbl
  integer :: str,unstr, iloc, counter
  real, allocatable, intent(out) :: told_style(:,:,:)
  allocate(told_style(nloc, totele_str,totele_unst))

  counter=1
  do unstr=1,totele_unst
    do str=1,totele_str
      do iloc=1,nloc
        told_style(iloc,str,unstr)=array(counter)
        counter= counter+1
      end do
    end do
  end do
end subroutine array_to_told_semi



subroutine matrix_to_array(matrix, array)
  implicit none
  !external vbl
  real, intent(in) :: matrix(:,:)
  ! local vbl
  real, allocatable, intent(out) :: array(:)
  integer :: no_rows, no_cols, i, j, counter

no_rows = size(matrix(:,1))
no_cols = size(matrix(1,:))
counter=1
do i=1,no_rows
  do j=1,no_cols
    array(counter) = matrix(i,j)
    counter = counter+1
  end do
end do
end subroutine matrix_to_array


!@brief>:: this subroutines initialises CSR matrix and put numbers for row & column
! sparxe_matrix%g_iloc :: refers to global row number in CSR format
! sparxe_matrix%g_jloc :: refers to global colum number in CSR format
! Amin in case of semi-structured instead on Npos/=0 I need to use overlaps or boundary conditions
subroutine make_sparse_matrix(sparse_matrix, totele, nloc)
  implicit none
  !external vbl
  ! real, intent(in) :: x_all(:,:,:)
  integer :: totele, nloc
  !local vbl
  type(sparse), intent(out) :: sparse_matrix
  integer :: i, ele, iloc, jloc, counter, counter2
  integer:: Npos, Nside, lnodes2(2)
  ! logical :: T, F
  ! T=.true.
  ! F=.false.

  if (allocated(sparse_matrix%g_iloc) ) deallocate(sparse_matrix%g_iloc)
  if (allocated(sparse_matrix%g_jloc) ) deallocate(sparse_matrix%g_jloc)
  if (allocated(sparse_matrix%val) ) deallocate(sparse_matrix%val)

  allocate(sparse_matrix%g_iloc(totele*nloc), sparse_matrix%g_jloc(totele*nloc*(nloc)))
  allocate(sparse_matrix%val(totele*nloc*(nloc)))
  counter=1
  counter2=1
  do ele=1,totele

    do iloc=1,nloc

      sparse_matrix%g_iloc(counter)= counter2
      counter = counter+1

      do jloc=1,nloc
        sparse_matrix%g_jloc(counter2)= glob_no(ele, nloc, jloc)
        counter2 = counter2+1
      end do ! ijloc

    end do ! iloc
  end do ! end ele
end subroutine make_sparse_matrix



! !@brief>:: this subroutines initialises CSR matrix and put numbers for row & column
! ! sparxe_matrix%g_iloc :: refers to global row number in CSR format
! ! sparxe_matrix%g_jloc :: refers to global colum number in CSR format
! ! Amin in case of semi-structured instead on Npos/=0 I need to use overlaps or boundary conditions
! subroutine make_sparse_matrix(sparse_matrix, meshList, totele, nloc)
!   implicit none
!   !external vbl
!   type(Mesh), intent(in), dimension(:) :: meshList
!   ! real, intent(in) :: x_all(:,:,:)
!   integer :: totele, nloc
!   !local vbl
!   type(sparse), intent(out) :: sparse_matrix
!   integer :: i, ele, iloc, jloc, counter, counter2
!   integer:: Npos, Nside, lnodes2(2)
!   ! logical :: T, F
!   ! T=.true.
!   ! F=.false.
!
!
!   allocate(sparse_matrix%g_iloc(totele*nloc), sparse_matrix%g_jloc(totele*nloc*(nloc)))
!   allocate(sparse_matrix%val(totele*nloc*(nloc)))
! print*, 'size g_iloc:', size(sparse_matrix%g_iloc),size(sparse_matrix%g_jloc),size(sparse_matrix%val)
!   counter=1
!   counter2=1
!   do ele=1,totele
!     do iloc=1,nloc
!       sparse_matrix%g_iloc(counter)= counter2
! PRINT*, COUNTER
!       counter = counter+1
!       do jloc=1,3
!         sparse_matrix%g_jloc(counter2)= glob_no(ele,nloc, jloc)
! print*, 'glob_no(ele,nloc, jloc):', counter2, glob_no(ele,nloc, jloc)
!         counter2 = counter2+1
!       end do ! ijloc
!
!       if ( iloc==1 ) then
!         call getNeigDataMesh(meshlist, ele, 1, Npos, Nside, lnodes2)
!         if ( Npos/=0 ) then
!           ! do i=1,nloc
!             ! if ( x_all(1,1,ele)==x_all(1,i,Npos) .and. x_all(2,1,ele)==x_all(2,i,Npos) ) then
!             !   sparse_matrix%g_jloc(counter2) = glob_no(Npos,nloc,i)
!             !   counter2 = counter2 +1
!             ! end if
!             if ( meshList(ele)%Dir(1) ) then
!               sparse_matrix%g_jloc(counter2) = glob_no(Npos,nloc,lnodes2(2))
! print*, 'glob_no(ele,nloc, jloc):', counter2, glob_no(ele,nloc, jloc)
!               counter2 = counter2 +1
!               sparse_matrix%g_jloc(counter2) = glob_no(Npos,nloc,lnodes2(1))
! print*, 'glob_no(ele,nloc, jloc):', counter2, glob_no(ele,nloc, jloc)
!               counter2 = counter2 +1
!             else
!               sparse_matrix%g_jloc(counter2) = glob_no(Npos,nloc,lnodes2(1))
! print*, 'glob_no(ele,nloc, jloc):', counter2, glob_no(ele,nloc, jloc)
!               counter2 = counter2 +1
!               sparse_matrix%g_jloc(counter2) = glob_no(Npos,nloc,lnodes2(2))
! print*, 'glob_no(ele,nloc, jloc):', counter2, glob_no(ele,nloc, jloc)
!               counter2 = counter2 +1
!             end if
!           ! end do
!         end if
!         call getNeigDataMesh(meshlist, ele, 2, Npos, Nside, lnodes2)
!
!         if ( Npos/=0 ) then
!           ! do i=1,nloc
!             ! if ( x_all(1,1,ele)==x_all(1,i,Npos) .and. x_all(2,1,ele)==x_all(2,i,Npos) ) then
!             !   sparse_matrix%g_jloc(counter2) = glob_no(Npos,nloc,i)
!             !   counter2 = counter2 +1
!             if ( meshList(ele)%Dir(2) ) then
!               sparse_matrix%g_jloc(counter2) = glob_no(Npos,nloc,lnodes2(2))
!               counter2 = counter2 +1
! print*, 'glob_no(ele,nloc, jloc):', counter2, glob_no(ele,nloc, jloc)
!               sparse_matrix%g_jloc(counter2) = glob_no(Npos,nloc,lnodes2(1))
!               counter2 = counter2 +1
! print*, 'glob_no(ele,nloc, jloc):', counter2, glob_no(ele,nloc, jloc)
!             else
!               sparse_matrix%g_jloc(counter2) = glob_no(Npos,nloc,lnodes2(1))
!               counter2 = counter2 +1
! print*, 'glob_no(ele,nloc, jloc):', counter2, glob_no(ele,nloc, jloc)
!               sparse_matrix%g_jloc(counter2) = glob_no(Npos,nloc,lnodes2(2))
!               counter2 = counter2 +1
!             print*, 'glob_no(ele,nloc, jloc):', counter2, glob_no(ele,nloc, jloc)
!             end if
!           ! end do
!         end if
!
!
!       elseif ( iloc==2 ) then
!         call getNeigDataMesh(meshlist, ele, 2, Npos, Nside, lnodes2)
!
!         if ( Npos/=0 ) then
!           ! do i=1,nloc
!             ! if ( x_all(1,2,ele)==x_all(1,i,Npos) .and. x_all(2,2,ele)==x_all(2,i,Npos) ) then
!             !   sparse_matrix%g_jloc(counter2) = glob_no(Npos,iloc,i)
!             !   counter2 = counter2 +1
!             ! end if
!             if ( meshList(ele)%Dir(2) ) then
!               sparse_matrix%g_jloc(counter2) = glob_no(Npos,nloc,lnodes2(2))
! print*, 'glob_no(ele,nloc, jloc):', counter2, glob_no(ele,nloc, jloc)
!               counter2 = counter2 +1
!               sparse_matrix%g_jloc(counter2) = glob_no(Npos,nloc,lnodes2(1))
! print*, 'glob_no(ele,nloc, jloc):', counter2, glob_no(ele,nloc, jloc)
!               counter2 = counter2 +1
!             else
!               sparse_matrix%g_jloc(counter2) = glob_no(Npos,nloc,lnodes2(1))
! print*, 'glob_no(ele,nloc, jloc):', counter2, glob_no(ele,nloc, jloc)
!               counter2 = counter2 +1
!               sparse_matrix%g_jloc(counter2) = glob_no(Npos,nloc,lnodes2(2))
! print*, 'glob_no(ele,nloc, jloc):', counter2, glob_no(ele,nloc, jloc)
!               counter2 = counter2 +1
!             end if
!           ! end do
!         end if
!         call getNeigDataMesh(meshlist, ele, 3, Npos, Nside, lnodes2)
!
!         if (  Npos/=0 ) then
!           ! do i=1,nloc
!             ! if ( x_all(1,2,ele)==x_all(1,i,Npos) .and. x_all(2,2,ele)==x_all(2,i,Npos) ) then
!             !   sparse_matrix%g_jloc(counter2) = glob_no(Npos,nloc,i)
!             !   counter2 = counter2 +1
!             ! end if
!             if ( meshList(ele)%Dir(3) ) then
!               sparse_matrix%g_jloc(counter2) = glob_no(Npos,nloc,lnodes2(2))
! print*, 'glob_no(ele,nloc, jloc):', counter2, glob_no(ele,nloc, jloc)
!               counter2 = counter2 +1
!               sparse_matrix%g_jloc(counter2) = glob_no(Npos,nloc,lnodes2(1))
! print*, 'glob_no(ele,nloc, jloc):', counter2, glob_no(ele,nloc, jloc)
!               counter2 = counter2 +1
!             else
!               sparse_matrix%g_jloc(counter2) = glob_no(Npos,nloc,lnodes2(1))
! print*, 'glob_no(ele,nloc, jloc):', counter2, glob_no(ele,nloc, jloc)
!               counter2 = counter2 +1
!               sparse_matrix%g_jloc(counter2) = glob_no(Npos,nloc,lnodes2(2))
! print*, 'glob_no(ele,nloc, jloc):', counter2, glob_no(ele,nloc, jloc)
!               counter2 = counter2 +1
!             end if
!           ! end do
!         end if
!
!       elseif ( iloc==3 ) then
!         call getNeigDataMesh(meshlist, ele, 1, Npos, Nside, lnodes2)
!
!         if ( Npos/=0 ) then
!           ! do i=1,nloc
!           !   if ( x_all(1,3,ele)==x_all(1,i,Npos) .and. x_all(2,3,ele)==x_all(2,i,Npos) ) then
!           !     sparse_matrix%g_jloc(counter2) = glob_no(Npos,iloc,i)
!           !     counter2 = counter2 +1
!           !   end if
!           ! end do
!           if ( meshList(ele)%Dir(1) ) then
!             sparse_matrix%g_jloc(counter2) = glob_no(Npos,nloc,lnodes2(2))
! print*, 'glob_no(ele,nloc, jloc):', counter2, glob_no(ele,nloc, jloc)
!             counter2 = counter2 +1
!             sparse_matrix%g_jloc(counter2) = glob_no(Npos,nloc,lnodes2(1))
! print*, 'glob_no(ele,nloc, jloc):', counter2, glob_no(ele,nloc, jloc)
!             counter2 = counter2 +1
!           else
!             sparse_matrix%g_jloc(counter2) = glob_no(Npos,nloc,lnodes2(1))
! print*, 'glob_no(ele,nloc, jloc):', counter2, glob_no(ele,nloc, jloc)
!             counter2 = counter2 +1
!             sparse_matrix%g_jloc(counter2) = glob_no(Npos,nloc,lnodes2(2))
! print*, 'glob_no(ele,nloc, jloc):', counter2, glob_no(ele,nloc, jloc)
!             counter2 = counter2 +1
!           end if
!         end if
!         call getNeigDataMesh(meshlist, ele, 3, Npos, Nside, lnodes2)
!
!         if ( Npos/=0 ) then
!           ! do i=1,nloc
!           !   if ( x_all(1,3,ele)==x_all(1,i,Npos) .and. x_all(2,3,ele)==x_all(2,i,Npos) ) then
!           !     sparse_matrix%g_jloc(counter2) = glob_no(Npos,nloc,i)
!           !     counter2 = counter2 +1
!           !   end if
!           ! end do
!           if ( meshList(ele)%Dir(3) ) then
!             sparse_matrix%g_jloc(counter2) = glob_no(Npos,nloc,lnodes2(2))
! print*, 'glob_no(ele,nloc, jloc):', counter2, glob_no(ele,nloc, jloc)
!             counter2 = counter2 +1
!             sparse_matrix%g_jloc(counter2) = glob_no(Npos,nloc,lnodes2(1))
! print*, 'glob_no(ele,nloc, jloc):', counter2, glob_no(ele,nloc, jloc)
!             counter2 = counter2 +1
!           else
!             sparse_matrix%g_jloc(counter2) = glob_no(Npos,nloc,lnodes2(1))
! print*, 'glob_no(ele,nloc, jloc):', counter2, glob_no(ele,nloc, jloc)
!             counter2 = counter2 +1
!             sparse_matrix%g_jloc(counter2) = glob_no(Npos,nloc,lnodes2(2))
! print*, 'glob_no(ele,nloc, jloc):', counter2, glob_no(ele,nloc, jloc)
!             counter2 = counter2 +1
!           end if
!         end if
!       end if
!     end do ! iloc
!
!   end do ! end ele
! end subroutine make_sparse_matrix



!@brief>:: this subroutines initialises CSR matrix for structured tri mesh and put numbers for row & column
! sparxe_matrix%g_iloc :: refers to global row number in CSR format
! sparxe_matrix%g_jloc :: refers to global colum number in CSR format
! Amin in case of semi-structured instead on Npos/=0 I need to use overlaps or boundary conditions
subroutine make_sparse_matrix_flux_str(sparse_matrix, no_ele_row, totele, nloc)
  implicit none
  !external vbl
  integer, intent(in) :: totele, nloc, no_ele_row
  !local vbl
  type(sparse), intent(out) :: sparse_matrix
  integer :: i, ele, iloc, jloc, counter, counter2
  integer:: Npos, Nside
  ! logical :: T, F
  ! T=.true.
  ! F=.false.


  allocate(sparse_matrix%g_iloc(totele*nloc), sparse_matrix%g_jloc(totele*nloc*(nloc+2*(nloc))))
  allocate(sparse_matrix%val(totele*nloc*(nloc+2*(nloc))))
  ! +2*nloc here is for 2 possible neighbours of the node

  ! sparse_matrix%g_iloc=0
  ! sparse_matrix%g_jloc=0
  sparse_matrix%val=0.0
  counter=1
  counter2=1
  do ele=1,totele
    do iloc=1,nloc
      sparse_matrix%g_iloc(counter)= counter2
      counter = counter+1
      ! first loop over 3 inner nloc nodes
      do jloc=1,3
        sparse_matrix%g_jloc(counter2)= glob_no(ele,nloc, jloc)
        counter2 = counter2+1
      end do ! ijloc
      ! now loop over neig ele nloc nodes
      if ( iloc==1 ) then
        call tri_ele_info2(totele, ele, 1, Npos, no_ele_row)
        if ( Npos/=0 ) then
              sparse_matrix%g_jloc(counter2) = glob_no(Npos,nloc,3)
              counter2 = counter2 +1
              sparse_matrix%g_jloc(counter2) = glob_no(Npos,nloc,1)
              counter2 = counter2 +1
              sparse_matrix%g_jloc(counter2) = glob_no(Npos, nloc, 2)
              counter2 = counter2 +1

        elseif ( Npos==0 ) then
          sparse_matrix%g_jloc(counter2) = glob_no(ele,nloc,1)
          counter2 = counter2 +1
          sparse_matrix%g_jloc(counter2) = glob_no(ele,nloc,3)
          counter2 = counter2 +1
          sparse_matrix%g_jloc(counter2) = glob_no(ele,nloc,2)
          counter2 = counter2 +1
        end if

        call tri_ele_info2(totele, ele, 3, Npos, no_ele_row)
        if ( Npos/=0 ) then
              sparse_matrix%g_jloc(counter2) = glob_no(Npos,nloc,2)
              counter2 = counter2 +1
              sparse_matrix%g_jloc(counter2) = glob_no(Npos,nloc,1)
              counter2 = counter2 +1
              sparse_matrix%g_jloc(counter2) = glob_no(Npos, nloc, 3)
              counter2 = counter2 +1
        elseif ( Npos==0 ) then
          sparse_matrix%g_jloc(counter2) = glob_no(ele,nloc,1)
          counter2 = counter2 +1
          sparse_matrix%g_jloc(counter2) = glob_no(ele,nloc,2)
          counter2 = counter2 +1
          sparse_matrix%g_jloc(counter2) = glob_no(ele,nloc,3)
          counter2 = counter2 +1
        end if


      elseif ( iloc==2 ) then
        call tri_ele_info2(totele, ele, 3, Npos, no_ele_row)
        if ( Npos/=0 ) then
              sparse_matrix%g_jloc(counter2) = glob_no(Npos,nloc,1)
              counter2 = counter2 +1
              sparse_matrix%g_jloc(counter2) = glob_no(Npos,nloc,2)
              counter2 = counter2 +1
              sparse_matrix%g_jloc(counter2) = glob_no(Npos, nloc, 3)
              counter2 = counter2 +1

        elseif ( Npos==0 ) then
          sparse_matrix%g_jloc(counter2) = glob_no(ele,nloc,2)
          counter2 = counter2 +1
          sparse_matrix%g_jloc(counter2) = glob_no(ele,nloc,1)
          counter2 = counter2 +1
          sparse_matrix%g_jloc(counter2) = glob_no(ele,nloc,3)
          counter2 = counter2 +1
        end if

        call tri_ele_info2(totele, ele, 2, Npos, no_ele_row)
        if ( Npos/=0 ) then
              sparse_matrix%g_jloc(counter2) = glob_no(Npos,nloc,3)
              counter2 = counter2 +1
              sparse_matrix%g_jloc(counter2) = glob_no(Npos,nloc,2)
              counter2 = counter2 +1
              sparse_matrix%g_jloc(counter2) = glob_no(Npos, nloc, 1)
              counter2 = counter2 +1
        elseif ( Npos==0 ) then
          sparse_matrix%g_jloc(counter2) = glob_no(ele,nloc,2)
          counter2 = counter2 +1
          sparse_matrix%g_jloc(counter2) = glob_no(ele,nloc,3)
          counter2 = counter2 +1
          sparse_matrix%g_jloc(counter2) = glob_no(ele,nloc,1)
          counter2 = counter2 +1
        end if


      elseif ( iloc==3 ) then
        call tri_ele_info2(totele, ele, 1, Npos, no_ele_row)
        if ( Npos/=0 ) then
              sparse_matrix%g_jloc(counter2) = glob_no(Npos,nloc,1)
              counter2 = counter2 +1
              sparse_matrix%g_jloc(counter2) = glob_no(Npos,nloc,3)
              counter2 = counter2 +1
              sparse_matrix%g_jloc(counter2) = glob_no(Npos, nloc, 2)
              counter2 = counter2 +1

        elseif ( Npos==0 ) then
          sparse_matrix%g_jloc(counter2) = glob_no(ele,nloc,3)
          counter2 = counter2 +1
          sparse_matrix%g_jloc(counter2) = glob_no(ele,nloc,1)
          counter2 = counter2 +1
          sparse_matrix%g_jloc(counter2) = glob_no(ele,nloc,2)
          counter2 = counter2 +1
        end if

        call tri_ele_info2(totele, ele, 2, Npos, no_ele_row)
        if ( Npos/=0 ) then
              sparse_matrix%g_jloc(counter2) = glob_no(Npos,nloc,2)
              counter2 = counter2 +1
              sparse_matrix%g_jloc(counter2) = glob_no(Npos,nloc,3)
              counter2 = counter2 +1
              sparse_matrix%g_jloc(counter2) = glob_no(Npos, nloc, 1)
              counter2 = counter2 +1
        elseif ( Npos==0 ) then
          sparse_matrix%g_jloc(counter2) = glob_no(ele,nloc,3)
          counter2 = counter2 +1
          sparse_matrix%g_jloc(counter2) = glob_no(ele,nloc,2)
          counter2 = counter2 +1
          sparse_matrix%g_jloc(counter2) = glob_no(ele,nloc,1)
          counter2 = counter2 +1
        end if
      end if
    end do ! iloc

  end do ! end ele
end subroutine make_sparse_matrix_flux_str





!@brief>:: this subroutines initialises CSR matrix for unstr mesh and put numbers for row & column
! sparxe_matrix%g_iloc :: refers to global row number in CSR format
! sparxe_matrix%g_jloc :: refers to global colum number in CSR format
! Amin in case of semi-structured instead on Npos/=0 I need to use overlaps or boundary conditions
subroutine make_sparse_matrix_flux(sparse_matrix, meshList, totele, nloc)
  implicit none
  !external vbl
  type(Mesh), intent(in), dimension(:) :: meshList
  integer, intent(in) :: totele, nloc
  !local vbl
  type(sparse), intent(out) :: sparse_matrix
  integer :: i, ele, iloc, jloc, counter, counter2
  integer:: Npos, Nside, lnodes2(2)
  ! logical :: T, F
  ! T=.true.
  ! F=.false.


  allocate(sparse_matrix%g_iloc(totele*nloc), sparse_matrix%g_jloc(totele*nloc*(nloc+2*(nloc))))
  allocate(sparse_matrix%val(totele*nloc*(nloc+2*(nloc))))
  ! +2*nloc here is for 2 possible neighbours of the node

  ! sparse_matrix%g_iloc=0
  ! sparse_matrix%g_jloc=0
  sparse_matrix%val=0.0
  counter=1
  counter2=1
  do ele=1,totele
    do iloc=1,nloc
      sparse_matrix%g_iloc(counter)= counter2
      counter = counter+1
      ! first loop over 3 inner nloc nodes
      do jloc=1,3
        sparse_matrix%g_jloc(counter2)= glob_no(ele,nloc, jloc)
        counter2 = counter2+1
      end do ! ijloc
      ! now loop over neig ele nloc nodes
      if ( iloc==1 ) then
        call getNeigDataMesh(meshlist, ele, 1, Npos, Nside, lnodes2)
        if ( Npos/=0 ) then
            if ( meshList(ele)%Dir(1) ) then
              sparse_matrix%g_jloc(counter2) = glob_no(Npos,nloc,lnodes2(2))
              counter2 = counter2 +1
              sparse_matrix%g_jloc(counter2) = glob_no(Npos,nloc,lnodes2(1))
              counter2 = counter2 +1
              !Amin :: here the 3rd Npos iloc should be added because we're looping over nloc
              do i=1,nloc
                if ( i/=lnodes2(1) .and. i/=lnodes2(2) ) then
                  sparse_matrix%g_jloc(counter2) = glob_no(Npos, nloc, i)
                  counter2 = counter2 +1
                  exit
                end if
              end do
            else
              sparse_matrix%g_jloc(counter2) = glob_no(Npos,nloc,lnodes2(1))
              counter2 = counter2 +1
              sparse_matrix%g_jloc(counter2) = glob_no(Npos,nloc,lnodes2(2))
              counter2 = counter2 +1
              do i=1,nloc
                if ( i/=lnodes2(1) .and. i/=lnodes2(2) ) then
                  sparse_matrix%g_jloc(counter2) = glob_no(Npos, nloc, i)
                  counter2 = counter2 +1
                  exit
                end if
              end do
            end if
        elseif ( Npos==0 ) then
          ! Amin, I just defined 3 extra nodes after total length of %g_jloc to represent
          ! domain boundary indecies.
          sparse_matrix%g_jloc(counter2) = glob_no(ele,nloc,1)
          counter2 = counter2 +1
          sparse_matrix%g_jloc(counter2) = glob_no(ele,nloc,3)
          counter2 = counter2 +1
          sparse_matrix%g_jloc(counter2) = glob_no(ele,nloc,2)
          counter2 = counter2 +1
        end if

        call getNeigDataMesh(meshlist, ele, 2, Npos, Nside, lnodes2)
        if ( Npos/=0 ) then
            if ( meshList(ele)%Dir(2) ) then
              sparse_matrix%g_jloc(counter2) = glob_no(Npos,nloc,lnodes2(2))
              counter2 = counter2 +1
              sparse_matrix%g_jloc(counter2) = glob_no(Npos,nloc,lnodes2(1))
              counter2 = counter2 +1
              do i=1,nloc
                if ( i/=lnodes2(1) .and. i/=lnodes2(2) ) then
                  sparse_matrix%g_jloc(counter2) = glob_no(Npos, nloc, i)
                  counter2 = counter2 +1
                  exit
                end if
              end do
            else
              sparse_matrix%g_jloc(counter2) = glob_no(Npos,nloc,lnodes2(1))
              counter2 = counter2 +1
              sparse_matrix%g_jloc(counter2) = glob_no(Npos,nloc,lnodes2(2))
              counter2 = counter2 +1
              do i=1,nloc
                if ( i/=lnodes2(1) .and. i/=lnodes2(2) ) then
                  sparse_matrix%g_jloc(counter2) = glob_no(Npos, nloc, i)
                  counter2 = counter2 +1
                  exit
                end if
              end do
            end if
        elseif ( Npos==0 ) then
          sparse_matrix%g_jloc(counter2) = glob_no(ele,nloc,1)
          counter2 = counter2 +1
          sparse_matrix%g_jloc(counter2) = glob_no(ele,nloc,2)
          counter2 = counter2 +1
          sparse_matrix%g_jloc(counter2) = glob_no(ele,nloc,3)
          counter2 = counter2 +1
        end if


      elseif ( iloc==2 ) then
        call getNeigDataMesh(meshlist, ele, 2, Npos, Nside, lnodes2)

        if ( Npos/=0 ) then
            if ( meshList(ele)%Dir(2) ) then
              sparse_matrix%g_jloc(counter2) = glob_no(Npos,nloc,lnodes2(2))
              counter2 = counter2 +1
              sparse_matrix%g_jloc(counter2) = glob_no(Npos,nloc,lnodes2(1))
              counter2 = counter2 +1
              do i=1,nloc
                if ( i/=lnodes2(1) .and. i/=lnodes2(2) ) then
                  sparse_matrix%g_jloc(counter2) = glob_no(Npos, nloc, i)
                  counter2 = counter2 +1
                  exit
                end if
              end do
            else
              sparse_matrix%g_jloc(counter2) = glob_no(Npos,nloc,lnodes2(1))
              counter2 = counter2 +1
              sparse_matrix%g_jloc(counter2) = glob_no(Npos,nloc,lnodes2(2))
              counter2 = counter2 +1
              do i=1,nloc
                if ( i/=lnodes2(1) .and. i/=lnodes2(2) ) then
                  sparse_matrix%g_jloc(counter2) = glob_no(Npos, nloc, i)
                  counter2 = counter2 +1
                  exit
                end if
              end do
            end if
          ! end do
        elseif ( Npos==0 ) then
          sparse_matrix%g_jloc(counter2) = glob_no(ele,nloc,2)
          counter2 = counter2 +1
          sparse_matrix%g_jloc(counter2) = glob_no(ele,nloc,1)
          counter2 = counter2 +1
          sparse_matrix%g_jloc(counter2) = glob_no(ele,nloc,3)
          counter2 = counter2 +1
        end if
        call getNeigDataMesh(meshlist, ele, 3, Npos, Nside, lnodes2)

        if (  Npos/=0 ) then
            if ( meshList(ele)%Dir(3) ) then
              sparse_matrix%g_jloc(counter2) = glob_no(Npos,nloc,lnodes2(2))
              counter2 = counter2 +1
              sparse_matrix%g_jloc(counter2) = glob_no(Npos,nloc,lnodes2(1))
              counter2 = counter2 +1
              do i=1,nloc
                if ( i/=lnodes2(1) .and. i/=lnodes2(2) ) then
                  sparse_matrix%g_jloc(counter2) = glob_no(Npos, nloc, i)
                  counter2 = counter2 +1
                  exit
                end if
              end do
            else
              sparse_matrix%g_jloc(counter2) = glob_no(Npos,nloc,lnodes2(1))
              counter2 = counter2 +1
              sparse_matrix%g_jloc(counter2) = glob_no(Npos,nloc,lnodes2(2))
              counter2 = counter2 +1
              do i=1,nloc
                if ( i/=lnodes2(1) .and. i/=lnodes2(2) ) then
                  sparse_matrix%g_jloc(counter2) = glob_no(Npos, nloc, i)
                  counter2 = counter2 +1
                  exit
                end if
              end do
            end if
          elseif ( Npos==0 ) then
            sparse_matrix%g_jloc(counter2) = glob_no(ele,nloc,2)
            counter2 = counter2 +1
            sparse_matrix%g_jloc(counter2) = glob_no(ele,nloc,3)
            counter2 = counter2 +1
            sparse_matrix%g_jloc(counter2) = glob_no(ele,nloc,1)
            counter2 = counter2 +1
        end if
      elseif ( iloc==3 ) then
        call getNeigDataMesh(meshlist, ele, 1, Npos, Nside, lnodes2)

        if ( Npos/=0 ) then
          if ( meshList(ele)%Dir(1) ) then
            sparse_matrix%g_jloc(counter2) = glob_no(Npos,nloc,lnodes2(2))
            counter2 = counter2 +1
            sparse_matrix%g_jloc(counter2) = glob_no(Npos,nloc,lnodes2(1))
            counter2 = counter2 +1
            do i=1,nloc
              if ( i/=lnodes2(1) .and. i/=lnodes2(2) ) then
                sparse_matrix%g_jloc(counter2) = glob_no(Npos, nloc, i)
                counter2 = counter2 +1
                exit
              end if
            end do
          else
            sparse_matrix%g_jloc(counter2) = glob_no(Npos,nloc,lnodes2(1))
            counter2 = counter2 +1
            sparse_matrix%g_jloc(counter2) = glob_no(Npos,nloc,lnodes2(2))
            counter2 = counter2 +1
            do i=1,nloc
              if ( i/=lnodes2(1) .and. i/=lnodes2(2) ) then
                sparse_matrix%g_jloc(counter2) = glob_no(Npos, nloc, i)
                counter2 = counter2 +1
                exit
              end if
            end do
          end if
        elseif ( Npos==0 ) then
          sparse_matrix%g_jloc(counter2) = glob_no(ele,nloc,3)
          counter2 = counter2 +1
          sparse_matrix%g_jloc(counter2) = glob_no(ele,nloc,1)
          counter2 = counter2 +1
          sparse_matrix%g_jloc(counter2) = glob_no(ele,nloc,2)
          counter2 = counter2 +1
        end if
        call getNeigDataMesh(meshlist, ele, 3, Npos, Nside, lnodes2)

        if ( Npos/=0 ) then
          if ( meshList(ele)%Dir(3) ) then
            sparse_matrix%g_jloc(counter2) = glob_no(Npos,nloc,lnodes2(2))
            counter2 = counter2 +1
            sparse_matrix%g_jloc(counter2) = glob_no(Npos,nloc,lnodes2(1))
            counter2 = counter2 +1
            do i=1,nloc
              if ( i/=lnodes2(1) .and. i/=lnodes2(2) ) then
                sparse_matrix%g_jloc(counter2) = glob_no(Npos, nloc, i)
                counter2 = counter2 +1
                exit
              end if
            end do
          else
            sparse_matrix%g_jloc(counter2) = glob_no(Npos,nloc,lnodes2(1))
            counter2 = counter2 +1
            sparse_matrix%g_jloc(counter2) = glob_no(Npos,nloc,lnodes2(2))
            counter2 = counter2 +1
            do i=1,nloc
              if ( i/=lnodes2(1) .and. i/=lnodes2(2) ) then
                sparse_matrix%g_jloc(counter2) = glob_no(Npos, nloc, i)
                counter2 = counter2 +1
                exit
              end if
            end do
          end if
        elseif ( Npos==0 ) then
          sparse_matrix%g_jloc(counter2) = glob_no(ele,nloc,3)
          counter2 = counter2 +1
          sparse_matrix%g_jloc(counter2) = glob_no(ele,nloc,2)
          counter2 = counter2 +1
          sparse_matrix%g_jloc(counter2) = glob_no(ele,nloc,1)
          counter2 = counter2 +1
        end if
      end if
    end do ! iloc

  end do ! end ele
end subroutine make_sparse_matrix_flux



!@brief>:: this subroutines initialises CSR matrix for unstr mesh and put numbers for row & column
! sparxe_matrix%g_iloc :: refers to global row number in CSR format
! sparxe_matrix%g_jloc :: refers to global colum number in CSR format
! Amin in case of semi-structured instead on Npos/=0 I need to use overlaps or boundary conditions
subroutine make_sparse_matrix_flux_semi(sparse_matrix, meshList, n_split, totele_unst, totele_str,nloc, nface,&
                                        str_neig)
  implicit none
  !external vbl
  type(Mesh), intent(in), dimension(:) :: meshList
  integer, intent(in) :: totele_unst, totele_str, nloc,n_split, nface
  integer, dimension(:,:), intent(in) :: str_neig
  !local vbl
  type(sparse), intent(out) :: sparse_matrix
  integer :: un_ele, str_ele, iloc, jloc, counter, counter2
  integer:: Npos, Nside, iface, irow, ipos, ele22, target_ele, lnodes2(2), orientation

  if (allocated(sparse_matrix%g_iloc) ) deallocate(sparse_matrix%g_iloc)
  if (allocated(sparse_matrix%g_jloc) ) deallocate(sparse_matrix%g_jloc)
  if (allocated(sparse_matrix%val) ) deallocate(sparse_matrix%val)

  allocate(sparse_matrix%g_iloc(totele_unst*totele_str*nloc))
  allocate(sparse_matrix%g_jloc(totele_unst*totele_str*nloc*(nloc+2*(nloc))))
  allocate(sparse_matrix%val(totele_unst*totele_str*nloc*(nloc+2*(nloc))))
  ! +2*nloc here is for 2 possible neighbours of the node

  sparse_matrix%val=0.0
  counter=1
  counter2=1
  do un_ele=1,totele_unst
    do str_ele=1,totele_str
      call get_str_info(n_split, str_ele, irow, ipos, orientation)
      ipos=int(ipos/2)+1
      do iloc=1,nloc
        sparse_matrix%g_iloc(counter)= counter2
        counter = counter+1
        do jloc=1,3   !loop over 3 inner nloc nodes
          sparse_matrix%g_jloc(counter2)= glob_no_semi(un_ele, str_ele, n_split, nloc, jloc)
          counter2 = counter2+1
        end do ! jloc

        if ( iloc==1 ) then
          do iface=1,nface
              if ( iface==1 ) then
                    call getNeigDataMesh(meshlist, un_ele, 1, Npos, Nside, lnodes2)
                    ele22 = str_neig(iface,str_ele) ! neighbour within str elements
                    if ( ele22/=0 ) then !str_ele is within un_ele
                          do jloc=1,nloc
                            sparse_matrix%g_jloc(counter2) = glob_no_semi(un_ele, ele22, n_split, nloc, jloc)
                            counter2 = counter2 +1
                          end do

                    elseif ( ele22==0 .and. Npos/=0 ) then !str_ele on bc of un_ele which un_ele is within domain
                          ! Npos = meshList(un_ele)%Neig(iface)
                          ! Nside = meshList(un_ele)%fNeig(iface)
                          target_ele = meshList(un_ele)%s_ele(ipos, 1)
                          do jloc=1,nloc
                            sparse_matrix%g_jloc(counter2) = (target_ele -1)*nloc + jloc
                            counter2 = counter2 +1
                          end do

                    elseif ( ele22==0 .and. Npos==0 ) then !str_ele on bc of un_ele which un_ele is on bc of domain
                          do jloc=1,nloc
                            sparse_matrix%g_jloc(counter2) = glob_no_semi(un_ele, str_ele, n_split, nloc, jloc) ! not sure what to put here. I just put str_ele
                            counter2 = counter2 +1
                          end do
                    end if !ele22

              elseif ( iface==3 ) then
                    call getNeigDataMesh(meshlist, un_ele, 2, Npos, Nside, lnodes2)
                    ele22 = str_neig(iface,str_ele) ! neighbour within str elements
                    if ( ele22/=0 ) then !str_ele is within un_ele
                          do jloc=1,nloc
                            sparse_matrix%g_jloc(counter2) = glob_no_semi(un_ele, ele22, n_split, nloc, jloc)
                            counter2 = counter2 +1
                          end do

                    elseif ( ele22==0 .and. Npos/=0 ) then !str_ele on bc of un_ele which un_ele is within domain
                          ! Npos = meshList(un_ele)%Neig(iface)
                          ! Nside = meshList(un_ele)%fNeig(iface)
                          target_ele = meshList(un_ele)%s_ele(irow, 2)
                          do jloc=1,nloc
                            sparse_matrix%g_jloc(counter2) = (target_ele -1)*nloc + jloc
                            counter2 = counter2 +1
                          end do

                    elseif ( ele22==0 .and. Npos==0 ) then !str_ele on bc of un_ele which un_ele is on bc of domain
                          do jloc=1,nloc
                            sparse_matrix%g_jloc(counter2) = glob_no_semi(un_ele, str_ele, n_split, nloc, jloc) ! not sure what to put here. I just put str_ele
                            counter2 = counter2 +1
                          end do
                    end if !ele22
              end if !iface
            end do !iface

        elseif ( iloc==2 ) then
          do iface=1,nface
              if ( iface==2 ) then
                    call getNeigDataMesh(meshlist, un_ele, 3, Npos, Nside, lnodes2)
                    ele22 = str_neig(iface,str_ele) ! neighbour within str elements
                    if ( ele22/=0 ) then !str_ele is within un_ele
                          do jloc=1,nloc
                            sparse_matrix%g_jloc(counter2) = glob_no_semi(un_ele, ele22, n_split, nloc, jloc)
                            counter2 = counter2 +1
                          end do

                    elseif ( ele22==0 .and. Npos/=0 ) then !str_ele on bc of un_ele which un_ele is within domain
                          ! Npos = meshList(un_ele)%Neig(iface)
                          ! Nside = meshList(un_ele)%fNeig(iface)
                          target_ele = meshList(un_ele)%s_ele(irow, 3)
                          do jloc=1,nloc
                            sparse_matrix%g_jloc(counter2) = (target_ele -1)*nloc + jloc
                            counter2 = counter2 +1
                          end do

                    elseif ( ele22==0 .and. Npos==0 ) then !str_ele on bc of un_ele which un_ele is on bc of domain
                          do jloc=1,nloc
                            sparse_matrix%g_jloc(counter2) = glob_no_semi(un_ele, str_ele, n_split, nloc, jloc) ! not sure what to put here. I just put str_ele
                            counter2 = counter2 +1
                          end do
                    end if !ele22

              elseif( iface==3 ) then
                    call getNeigDataMesh(meshlist, un_ele, 2, Npos, Nside, lnodes2)
                    ele22 = str_neig(iface,str_ele) ! neighbour within str elements
                    if ( ele22/=0 ) then !str_ele is within un_ele
                          do jloc=1,nloc
                            sparse_matrix%g_jloc(counter2) = glob_no_semi(un_ele, ele22, n_split, nloc, jloc)
                            counter2 = counter2 +1
                          end do

                    elseif ( ele22==0 .and. Npos/=0 ) then !str_ele on bc of un_ele which un_ele is within domain
                          ! Npos = meshList(un_ele)%Neig(iface)
                          ! Nside = meshList(un_ele)%fNeig(iface)
                          target_ele = meshList(un_ele)%s_ele(irow, 2)
                          do jloc=1,nloc
                            sparse_matrix%g_jloc(counter2) = (target_ele -1)*nloc + jloc
                            counter2 = counter2 +1
                          end do

                    elseif ( ele22==0 .and. Npos==0 ) then !str_ele on bc of un_ele which un_ele is on bc of domain
                          do jloc=1,nloc
                            sparse_matrix%g_jloc(counter2) = glob_no_semi(un_ele, str_ele, n_split, nloc, jloc) ! not sure what to put here. I just put str_ele
                            counter2 = counter2 +1
                          end do
                    end if !ele22
              end if
          end do !iface

        elseif ( iloc==3 ) then
          do iface=1,nface
              if ( iface==1 ) then
                    call getNeigDataMesh(meshlist, un_ele, 1, Npos, Nside, lnodes2)
                    ele22 = str_neig(iface,str_ele) ! neighbour within str elements
                    if ( ele22/=0 ) then !str_ele is within un_ele
                          do jloc=1,nloc
                            sparse_matrix%g_jloc(counter2) = glob_no_semi(un_ele, ele22, n_split, nloc, jloc)
                            counter2 = counter2 +1
                          end do

                    elseif ( ele22==0 .and. Npos/=0 ) then !str_ele on bc of un_ele which un_ele is within domain
                          ! Npos = meshList(un_ele)%Neig(iface)
                          ! Nside = meshList(un_ele)%fNeig(iface)
                          target_ele = meshList(un_ele)%s_ele(ipos, 1)
                          do jloc=1,nloc
                            sparse_matrix%g_jloc(counter2) = (target_ele -1)*nloc + jloc
                            counter2 = counter2 +1
                          end do

                    elseif ( ele22==0 .and. Npos==0 ) then !str_ele on bc of un_ele which un_ele is on bc of domain
                          do jloc=1,nloc
                            sparse_matrix%g_jloc(counter2) = glob_no_semi(un_ele, str_ele, n_split, nloc, jloc) ! not sure what to put here. I just put str_ele
                            counter2 = counter2 +1
                          end do
                    end if !ele22

              elseif (iface==2) then
                    call getNeigDataMesh(meshlist, un_ele, 3, Npos, Nside, lnodes2)
                    ele22 = str_neig(iface,str_ele) ! neighbour within str elements
                    if ( ele22/=0 ) then !str_ele is within un_ele
                          do jloc=1,nloc
                            sparse_matrix%g_jloc(counter2) = glob_no_semi(un_ele, ele22, n_split, nloc, jloc)
                            counter2 = counter2 +1
                          end do

                    elseif ( ele22==0 .and. Npos/=0 ) then !str_ele on bc of un_ele which un_ele is within domain
                          ! Npos = meshList(un_ele)%Neig(iface)
                          ! Nside = meshList(un_ele)%fNeig(iface)
                          target_ele = meshList(un_ele)%s_ele(irow, 3)
                          do jloc=1,nloc
                            sparse_matrix%g_jloc(counter2) = (target_ele -1)*nloc + jloc
                            counter2 = counter2 +1
                          end do
                    elseif ( ele22==0 .and. Npos==0 ) then !str_ele on bc of un_ele which un_ele is on bc of domain
                          do jloc=1,nloc
                            sparse_matrix%g_jloc(counter2) = glob_no_semi(un_ele, str_ele, n_split, nloc, jloc) ! not sure what to put here. I just put str_ele
                            counter2 = counter2 +1
                          end do
                    end if !ele22
              end if !iface
          end do !iface
        end if !iloc

      end do !iloc
    end do !totele_str
  end do !totele_unst
end subroutine make_sparse_matrix_flux_semi



!@brief>:: this subroutines initialises CSR matrix for unstr mesh and put numbers for row & column
! sparxe_matrix%g_iloc :: refers to global row number in CSR format
! sparxe_matrix%g_jloc :: refers to global colum number in CSR format
! Amin in case of semi-structured instead on Npos/=0 I need to use overlaps or boundary conditions
! subroutine make_sparse_compact(sparse_matrix, meshList, n_split, totele_unst, totele_str,nloc, nface,&
!                                         str_neig, str_ele, un_ele)
!   implicit none
!   !external vbl
!   type(Mesh), intent(in), dimension(:) :: meshList
!   integer, intent(in) :: totele_unst, totele_str, nloc,n_split, nface, un_ele, str_ele
!   integer, dimension(:,:), intent(in) :: str_neig
!   !local vbl
!   type(sparse2), intent(out) :: sparse_matrix
!   integer :: iloc, jloc, counter, counter2
!   integer:: Npos, Nside, iface, irow, ipos, ele22, target_ele, lnodes2(2)
!
!   if (allocated(sparse_matrix%g_iloc) ) deallocate(sparse_matrix%g_iloc)
!   if (allocated(sparse_matrix%g_jloc) ) deallocate(sparse_matrix%g_jloc)
!   if (allocated(sparse_matrix%val) ) deallocate(sparse_matrix%val)
!
!   allocate(sparse_matrix%g_iloc(nloc))
!   allocate(sparse_matrix%g_jloc((nloc+2*(nloc))))
!   allocate(sparse_matrix%val(nloc,(nloc+2*(nloc))))
!   ! +2*nloc here is for 2 possible neighbours of the node
!
!   sparse_matrix%val=0.0
!   counter=1
!   counter2=1
!   ! do un_ele=1,totele_unst
!   !   do str_ele=1,totele_str
!       call get_str_info(n_split, str_ele, irow, ipos, orientation)
!       ipos=int(ipos/2)+1
!       do iloc=1,nloc
!         sparse_matrix%g_iloc(counter)= counter2
!         counter = counter+1
!         do jloc=1,3   !loop over 3 inner nloc nodes
!           sparse_matrix%g_jloc(counter2)= glob_no_semi(un_ele, str_ele, n_split, nloc, jloc)
!           counter2 = counter2+1
!         end do ! jloc
!
!         if ( iloc==1 ) then
!           do iface=1,nface
!               if ( iface==1 ) then
!                     call getNeigDataMesh(meshlist, un_ele, 1, Npos, Nside, lnodes2)
!                     ele22 = str_neig(iface,str_ele) ! neighbour within str elements
!                     if ( ele22/=0 ) then !str_ele is within un_ele
!                           do jloc=1,nloc
!                             sparse_matrix%g_jloc(counter2) = glob_no_semi(un_ele, ele22, n_split, nloc, jloc)
!                             counter2 = counter2 +1
!                           end do
!
!                     elseif ( ele22==0 .and. Npos/=0 ) then !str_ele on bc of un_ele which un_ele is within domain
!                           ! Npos = meshList(un_ele)%Neig(iface)
!                           ! Nside = meshList(un_ele)%fNeig(iface)
!                           target_ele = meshList(un_ele)%s_ele(ipos, 1)
!                           do jloc=1,nloc
!                             sparse_matrix%g_jloc(counter2) = (target_ele -1)*nloc + jloc
!                             counter2 = counter2 +1
!                           end do
!
!                     elseif ( ele22==0 .and. Npos==0 ) then !str_ele on bc of un_ele which un_ele is on bc of domain
!                           do jloc=1,nloc
!                             sparse_matrix%g_jloc(counter2) = glob_no_semi(un_ele, str_ele, n_split, nloc, jloc) ! not sure what to put here. I just put str_ele
!                             counter2 = counter2 +1
!                           end do
!                     end if !ele22
!
!               elseif ( iface==3 ) then
!                     call getNeigDataMesh(meshlist, un_ele, 2, Npos, Nside, lnodes2)
!                     ele22 = str_neig(iface,str_ele) ! neighbour within str elements
!                     if ( ele22/=0 ) then !str_ele is within un_ele
!                           do jloc=1,nloc
!                             sparse_matrix%g_jloc(counter2) = glob_no_semi(un_ele, ele22, n_split, nloc, jloc)
!                             counter2 = counter2 +1
!                           end do
!
!                     elseif ( ele22==0 .and. Npos/=0 ) then !str_ele on bc of un_ele which un_ele is within domain
!                           ! Npos = meshList(un_ele)%Neig(iface)
!                           ! Nside = meshList(un_ele)%fNeig(iface)
!                           target_ele = meshList(un_ele)%s_ele(irow, 2)
!                           do jloc=1,nloc
!                             sparse_matrix%g_jloc(counter2) = (target_ele -1)*nloc + jloc
!                             counter2 = counter2 +1
!                           end do
!
!                     elseif ( ele22==0 .and. Npos==0 ) then !str_ele on bc of un_ele which un_ele is on bc of domain
!                           do jloc=1,nloc
!                             sparse_matrix%g_jloc(counter2) = glob_no_semi(un_ele, str_ele, n_split, nloc, jloc) ! not sure what to put here. I just put str_ele
!                             counter2 = counter2 +1
!                           end do
!                     end if !ele22
!               end if !iface
!             end do !iface
!
!         elseif ( iloc==2 ) then
!           do iface=1,nface
!               if ( iface==2 ) then
!                     call getNeigDataMesh(meshlist, un_ele, 3, Npos, Nside, lnodes2)
!                     ele22 = str_neig(iface,str_ele) ! neighbour within str elements
!                     if ( ele22/=0 ) then !str_ele is within un_ele
!                           do jloc=1,nloc
!                             sparse_matrix%g_jloc(counter2) = glob_no_semi(un_ele, ele22, n_split, nloc, jloc)
!                             counter2 = counter2 +1
!                           end do
!
!                     elseif ( ele22==0 .and. Npos/=0 ) then !str_ele on bc of un_ele which un_ele is within domain
!                           ! Npos = meshList(un_ele)%Neig(iface)
!                           ! Nside = meshList(un_ele)%fNeig(iface)
!                           target_ele = meshList(un_ele)%s_ele(irow, 3)
!                           do jloc=1,nloc
!                             sparse_matrix%g_jloc(counter2) = (target_ele -1)*nloc + jloc
!                             counter2 = counter2 +1
!                           end do
!
!                     elseif ( ele22==0 .and. Npos==0 ) then !str_ele on bc of un_ele which un_ele is on bc of domain
!                           do jloc=1,nloc
!                             sparse_matrix%g_jloc(counter2) = glob_no_semi(un_ele, str_ele, n_split, nloc, jloc) ! not sure what to put here. I just put str_ele
!                             counter2 = counter2 +1
!                           end do
!                     end if !ele22
!
!               elseif( iface==3 ) then
!                     call getNeigDataMesh(meshlist, un_ele, 2, Npos, Nside, lnodes2)
!                     ele22 = str_neig(iface,str_ele) ! neighbour within str elements
!                     if ( ele22/=0 ) then !str_ele is within un_ele
!                           do jloc=1,nloc
!                             sparse_matrix%g_jloc(counter2) = glob_no_semi(un_ele, ele22, n_split, nloc, jloc)
!                             counter2 = counter2 +1
!                           end do
!
!                     elseif ( ele22==0 .and. Npos/=0 ) then !str_ele on bc of un_ele which un_ele is within domain
!                           ! Npos = meshList(un_ele)%Neig(iface)
!                           ! Nside = meshList(un_ele)%fNeig(iface)
!                           target_ele = meshList(un_ele)%s_ele(irow, 2)
!                           do jloc=1,nloc
!                             sparse_matrix%g_jloc(counter2) = (target_ele -1)*nloc + jloc
!                             counter2 = counter2 +1
!                           end do
!
!                     elseif ( ele22==0 .and. Npos==0 ) then !str_ele on bc of un_ele which un_ele is on bc of domain
!                           do jloc=1,nloc
!                             sparse_matrix%g_jloc(counter2) = glob_no_semi(un_ele, str_ele, n_split, nloc, jloc) ! not sure what to put here. I just put str_ele
!                             counter2 = counter2 +1
!                           end do
!                     end if !ele22
!               end if
!           end do !iface
!
!         elseif ( iloc==3 ) then
!           do iface=1,nface
!               if ( iface==1 ) then
!                     call getNeigDataMesh(meshlist, un_ele, 1, Npos, Nside, lnodes2)
!                     ele22 = str_neig(iface,str_ele) ! neighbour within str elements
!                     if ( ele22/=0 ) then !str_ele is within un_ele
!                           do jloc=1,nloc
!                             sparse_matrix%g_jloc(counter2) = glob_no_semi(un_ele, ele22, n_split, nloc, jloc)
!                             counter2 = counter2 +1
!                           end do
!
!                     elseif ( ele22==0 .and. Npos/=0 ) then !str_ele on bc of un_ele which un_ele is within domain
!                           ! Npos = meshList(un_ele)%Neig(iface)
!                           ! Nside = meshList(un_ele)%fNeig(iface)
!                           target_ele = meshList(un_ele)%s_ele(ipos, 1)
!                           do jloc=1,nloc
!                             sparse_matrix%g_jloc(counter2) = (target_ele -1)*nloc + jloc
!                             counter2 = counter2 +1
!                           end do
!
!                     elseif ( ele22==0 .and. Npos==0 ) then !str_ele on bc of un_ele which un_ele is on bc of domain
!                           do jloc=1,nloc
!                             sparse_matrix%g_jloc(counter2) = glob_no_semi(un_ele, str_ele, n_split, nloc, jloc) ! not sure what to put here. I just put str_ele
!                             counter2 = counter2 +1
!                           end do
!                     end if !ele22
!
!               elseif (iface==2) then
!                     call getNeigDataMesh(meshlist, un_ele, 3, Npos, Nside, lnodes2)
!                     ele22 = str_neig(iface,str_ele) ! neighbour within str elements
!                     if ( ele22/=0 ) then !str_ele is within un_ele
!                           do jloc=1,nloc
!                             sparse_matrix%g_jloc(counter2) = glob_no_semi(un_ele, ele22, n_split, nloc, jloc)
!                             counter2 = counter2 +1
!                           end do
!
!                     elseif ( ele22==0 .and. Npos/=0 ) then !str_ele on bc of un_ele which un_ele is within domain
!                           ! Npos = meshList(un_ele)%Neig(iface)
!                           ! Nside = meshList(un_ele)%fNeig(iface)
!                           target_ele = meshList(un_ele)%s_ele(irow, 3)
!                           do jloc=1,nloc
!                             sparse_matrix%g_jloc(counter2) = (target_ele -1)*nloc + jloc
!                             counter2 = counter2 +1
!                           end do
!                     elseif ( ele22==0 .and. Npos==0 ) then !str_ele on bc of un_ele which un_ele is on bc of domain
!                           do jloc=1,nloc
!                             sparse_matrix%g_jloc(counter2) = glob_no_semi(un_ele, str_ele, n_split, nloc, jloc) ! not sure what to put here. I just put str_ele
!                             counter2 = counter2 +1
!                           end do
!                     end if !ele22
!               end if !iface
!           end do !iface
!         end if !iloc
!
!       end do !iloc
!   !   end do !totele_str
!   ! end do !totele_unst
! end subroutine make_sparse_compact



!@brief>:: adds val to sparse_matrix at position [g_iloc,g_jloc]
subroutine add_to_CSR(sparse_matrix, g_iloc, g_jloc, val)
  implicit none
  !external vbl
  integer, intent(in) :: g_iloc, g_jloc
  real, intent(in) :: val
  type(sparse), intent(inout) :: sparse_matrix

  !local vbl
  integer :: A, i

  A = sparse_matrix%g_iloc(g_iloc)
  do i = A, A+6
    if ( sparse_matrix%g_jloc(i) == g_jloc ) then
      sparse_matrix%val(i) = val
    end if
  end do
end subroutine add_to_CSR





!@brief>:: adds val to sparse_matrix at position [g_iloc,g_jloc]
subroutine add_to_CSR_flux(sparse_matrix, g_iloc, g_jloc, val)
  implicit none
  !external vbl
  integer, intent(in) :: g_iloc, g_jloc
  real, intent(in) :: val
  type(sparse), intent(inout) :: sparse_matrix

  !local vbl
  integer :: A, i

  A = sparse_matrix%g_iloc(g_iloc)
  do i = A, A+8
    if ( sparse_matrix%g_jloc(i) == g_jloc ) then
      sparse_matrix%val(i) =  sparse_matrix%val(i)+ val
    end if
  end do
end subroutine add_to_CSR_flux



! !@brief>:: adds val to sparse_matrix at position [g_iloc,g_jloc]
! subroutine add_to_CSR_flux(sparse_matrix, g_iloc, g_jloc, val)
!   implicit none
!   !external vbl
!   integer, intent(in) :: g_iloc, g_jloc
!   real, intent(in) :: val
!   type(sparse), intent(inout) :: sparse_matrix
!
!   !local vbl
!   integer :: A, i
!
!   A = sparse_matrix%g_iloc(g_iloc)
!   do i = A, A+5
!     if ( sparse_matrix%g_jloc(i) == g_jloc ) then
!       sparse_matrix%val(i) = val
!     end if
!   end do
! end subroutine add_to_CSR_flux



! ! @brief>> This function works for semi-structured rid only
! ! in :: unstructured and structured ele number
! ! out :: global node numbers
! subroutine loc_to_glob_semi(un_ele, str_ele, nloc, glob_i)
!   implicit none
!   ! external vbls
!   integer, intent(in) :: un_ele, str_ele, nloc
!   ! internal vbls
!   integer, intent(inout) :: glob_i(nloc)
!   integer :: j
!
!   j=0
!   do while (j<nloc)
!     glob_i(nloc-j) = un_ele * str_ele * nloc - j
!     j=j+1
!   end do
! end subroutine loc_to_glob_semi


! @brief> :: Coverts local node id to global node id
integer function glob_no_semi(un_ele, str_ele, n_split, nloc, iloc) result(glob_id)
  implicit none
  integer, intent(in) :: un_ele, str_ele, n_split, nloc, iloc
    glob_id = nloc*(4**n_split)*(un_ele -1) + (str_ele -1)*nloc + iloc
end function glob_no_semi





! @brief> :: Coverts local node id to global node id
integer function glob_no(ele, nloc, iloc) result(glob_id)
  implicit none
  integer, intent(in) :: ele, nloc, iloc
    glob_id = (ele-1) * nloc + iloc
end function glob_no


!@brief>:: it receives matrix2 and converts iy to an array, then returns
! dot product of matrix1 and the new array and the result
subroutine matrix_dot_array(matrix1, matrix2, result)
  implicit none
  ! external vlb
  real, intent(in) :: matrix1(:,:), matrix2(:,:)
  ! local vbl
  real, intent(inout) :: result(:)
  real, allocatable :: array(:)
  integer :: counter, no_cols, no_rows,i ,j

  no_rows = size(matrix1(1,:))
  no_cols = size(matrix1(1,:))
  allocate(array(size(result)))

  counter =1
  do i=1,no_rows
    do j=1,no_cols
      array(counter)=matrix2(i,j)
      counter = counter +1
    end do
  end do

  counter=1
  do i=1,no_rows
    do j=1,no_cols
      result(counter) = result(counter) + matrix1(i,j)*array(j)
    end do
    counter = counter +1
  end do

  deallocate(array)
end subroutine matrix_dot_array



subroutine print_my_matrix(matrix, row)
  implicit none
  !external vbl
  real, intent(in) :: matrix(:,:)
  integer, intent(in) :: row
  integer :: i

  do i=1,row
    print*, i, matrix(i,:)
  end do
end subroutine print_my_matrix


!brief> this subroutine removes duplicated and zero elements of an array
subroutine remove_dups(example, res)
    ! https://rosettacode.org/wiki/Remove_duplicate_elements#Fortran
    implicit none
    integer :: n ! size of the sample
    integer:: k  ! number of unique elements
    integer:: i, j
    integer, dimension(:), intent(inout) :: example
    integer, allocatable, dimension(:), intent(inout) :: res
    integer, allocatable, dimension(:) :: res2 ! dummy array

    n= size(example)
    if (allocated(res)) deallocate(res)
    allocate(res(n))
    res=0
    k = 1
    res(1) = example(1)
    outer: do i=2,size(example)
        do j=1,k
            if (res(j) == example(i)) then
                cycle outer
            end if
        end do
        k = k + 1
        res(k) = example(i)
    end do outer
    ! allocate(res(count(res2/=0)))
    ! res=res2(:count(res2/=0))
    ! deallocate(res2)

    !===============================================
    ! n= size(example)
    ! if (allocated(res)) deallocate(res)
    ! allocate(res(n))
    !
    ! k = 1
    ! res(1) = example(1)
    ! outer: do i=2,size(example)
    !     do j=1,k
    !         if (res(j) == example(i)) then
    !             cycle outer
    !         end if
    !     end do
    !     k = k + 1
    !     res(k) = example(i)
    ! end do outer
    ! ! allocate(res(count(res2/=0)))
    ! ! res=res2(:count(res2/=0))
    ! ! deallocate(res2)
end subroutine remove_dups




! Finds inverse of a matrix
subroutine FINDInv(matrix, inverse, n, errorflag)
  ! Returns the inverse of a matrix calculated by finding the LU
  ! decomposition.  Depends on LAPACK.
  !Subroutine to find the inverse of a square matrix
  !Author : Louisda16th a.k.a Ashwith J. Rego
  !Reference : Algorithm has been well explained in:
  !http://math.uww.edu/~mcfarlat/inverse.htm
  !http://www.tutor.ms.unimelb.edu.au/matrix/matrix_inverse.html
  !https://www.dreamincode.net/forums/topic/366231-FORTRAN-90%3A-Matrix-Inversion/

  IMPLICIT NONE
  !Declarations
  INTEGER, INTENT(IN) :: n
  INTEGER, INTENT(OUT) :: errorflag  !Return error status. -1 for error, 0 for normal
  REAL, INTENT(IN), DIMENSION(n,n) :: matrix  !Input matrix
  REAL, INTENT(OUT), DIMENSION(n,n) :: inverse !Inverted matrix

  LOGICAL :: FLAG = .TRUE.
  INTEGER :: i, j, k, l
  REAL :: m
  REAL, DIMENSION(n,2*n) :: augmatrix !augmented matrix

  !Augment input matrix with an identity matrix
  DO i = 1, n
      DO j = 1, 2*n
          IF (j <= n ) THEN
              augmatrix(i,j) = matrix(i,j)
          ELSE IF ((i+n) == j) THEN
              augmatrix(i,j) = 1
          Else
              augmatrix(i,j) = 0
          ENDIF
      END DO
  END DO

  !Reduce augmented matrix to upper traingular form
  DO k =1, n-1
      IF (augmatrix(k,k) == 0) THEN
          FLAG = .FALSE.
          DO i = k+1, n
              IF (augmatrix(i,k) /= 0) THEN
                  DO j = 1,2*n
                      augmatrix(k,j) = augmatrix(k,j)+augmatrix(i,j)
                  END DO
                  FLAG = .TRUE.
                  EXIT
              ENDIF
              IF (FLAG .EQV. .FALSE.) THEN
                  PRINT*, "Matrix is non - invertible"
                  inverse = 0
                  errorflag = -1
                  return
              ENDIF
          END DO
      ENDIF
      DO j = k+1, n
          m = augmatrix(j,k)/augmatrix(k,k)
          DO i = k, 2*n
              augmatrix(j,i) = augmatrix(j,i) - m*augmatrix(k,i)
          END DO
      END DO
  END DO

  !Test for invertibility
  DO i = 1, n
      IF (augmatrix(i,i) == 0) THEN
          PRINT*, "Matrix is non - invertible"
          inverse = 0
          errorflag = -1
          return
      ENDIF
  END DO

  !Make diagonal elements as 1
  DO i = 1 , n
      m = augmatrix(i,i)
      DO j = i , (2 * n)
             augmatrix(i,j) = (augmatrix(i,j) / m)
      END DO
  END DO

  !Reduced right side half of augmented matrix to identity matrix
  DO k = n-1, 1, -1
      DO i =1, k
      m = augmatrix(i,k+1)
          DO j = k, (2*n)
              augmatrix(i,j) = augmatrix(i,j) -augmatrix(k+1,j) * m
          END DO
      END DO
  END DO

  !store answer
  DO i =1, n
      DO j = 1, n
          inverse(i,j) = augmatrix(i,j+n)
      END DO
  END DO
  errorflag = 0
end subroutine FINDinv

end module matrices
