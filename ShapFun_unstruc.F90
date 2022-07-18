module ShapFun_unstruc
  use Structures
  use ShapFun

  contains

  ! Amin added snloc to all subroutines
  subroutine get_shape_funs_spec_ustr(n, nlx, nlx_lxx, nlxx, sn_orig, weight, nlx_nod, &
              nloc, snloc, sngi, ngi, ndim, nface,n_s_list_no, face_sn, face_sn2, face_snlx, &
              sweight, npoly, ele_type, totele, face_list_no )
      implicit none
      ! nloc=no of local nodes per element
      ! ngi=no of quadrature points
      ! ndim=no of dimensions.
      integer, intent(in) :: nloc, sngi, ngi, ndim, nface, n_s_list_no, snloc, totele
      integer, INTENT(INOUT) :: face_list_no( nface, totele)
      ! shape functions....
      ! if .not.got_shape_funs then get the shape functions else assume we have them already
      ! n, nlx are the volume shape functions and their derivatives in local coordinates.
      ! weight are the weights of the quadrature points.
      ! nlx_nod are the derivatives of the local coordinates at the nods.
      ! nlx_lxx = the 3rd order local derivatives at the nodes.
      ! face info:
      ! face_ele(iface, ele) = given the face no iface and element no return the element next to
      ! the surface or if negative return the negative of the surface element number between element ele and face iface.
      ! face_list_no(iface, ele) returns the possible origantation number which defines the numbering
      ! of the non-zeros of the nabouting element.
      real, intent(inout) :: n(ngi,nloc), nlx(ngi,ndim,nloc), nlxx(ngi,nloc), nlx_lxx(ngi,ndim,nloc), weight(ngi)
      real, intent(inout) :: nlx_nod(nloc,ndim,nloc), sn_orig(sngi,snloc)
      real, intent(inout) :: face_sn(sngi,nloc,nface), face_sn2(sngi,nloc,n_s_list_no)
      real, intent(inout) :: face_snlx(sngi,ndim,nloc,nface)
      ! npoly=order of polynomial in Cartesian space; ele_type=type of element including order of poly.
      integer, intent(in) :: npoly,ele_type
      real, allocatable :: sweight(:)

      ! form the shape functions...
      call get_shape_funs_with_faces_unstr(n, nlx, sn_orig, weight,  &
                 nloc, snloc, sngi, ngi, ndim, nface,n_s_list_no, face_sn, face_sn2, face_snlx, &
                 sweight, npoly,ele_type, totele, face_list_no)

      ! Calculate high order derivatives of shape functions nlx_lxx, nlxx
      ! and also calculate node-wise deriavtives of shape functions nlx_nod.
      ! Amin for higher order uncomment get_high_order_shape_funs below
      ! call get_high_order_shape_funs(n, nlx, nlx_lxx, nlxx, weight, nlx_nod, nloc, ngi, ndim, &
      !              npoly, ele_type)
  end subroutine get_shape_funs_spec_ustr


  subroutine get_shape_funs_with_faces_unstr(n, nlx, sn_orig, weight,  &
                   nloc, snloc, sngi, ngi, ndim, nface,n_s_list_no, face_sn, face_sn2, face_snlx,  &
                   sweight, npoly,ele_type, totele, face_list_no)
      implicit none
      ! nloc=no of local nodes per element
      ! ngi=no of quadrature points
      ! ndim=no of dimensions.
      ! ele_type= element type
      integer, intent(in) :: nloc, sngi, ngi, ndim, nface, n_s_list_no, snloc
      integer, intent(in) :: totele
      integer, INTENT(INOUT) :: face_list_no( nface, totele)
      ! shape functions....
      ! if .not.got_shape_funs then get the shape functions else assume we have them already
      ! n, nlx are the volume shape functions and their derivatives in local coordinates.
      ! weight are the weights of the quadrature points.
      ! nlx_nod are the derivatives of the local coordinates at the nods.
      ! nlx_lxx = the 3rd order local derivatives at the nodes.
      ! face info:
      ! face_ele(iface, ele) = given the face no iface and element no return the element next to
      ! the surface or if negative return the negative of the surface element number between element ele and face iface.
      ! face_list_no(iface, ele) returns the possible origantation number which defines the numbering
      ! of the non-zeros of the nabouting element.
      real, intent(inout) :: n(ngi,nloc), nlx(ngi,ndim,nloc), weight(ngi)
      real, intent(inout) :: face_sn(sngi,nloc,nface), face_sn2(sngi,nloc,n_s_list_no)
      real, intent(inout) :: face_snlx(sngi,ndim,nloc,nface), sn_orig(sngi,snloc)
      ! npoly=order of polynomial in Cartesian space; ele_type=type of element including order of poly.
      integer, intent(in) :: npoly, ele_type
      ! local variables...
      integer is_triangle_or_tet
      parameter(is_triangle_or_tet=100)
      integer sndim, suf_ngi, suf_ndim, suf_nloc, ipoly, IQADRA
      real, allocatable :: rdum1(:), rdum2(:), rdum3(:)
      real, allocatable :: sn2(:,:),snlx_orig(:,:,:),sweight(:)
      real, allocatable :: suf_n(:,:),suf_nlx(:,:,:),suf_weight(:)
      ! integer, allocatable :: face_list_no(:,:)

      ! allocate memory...
      allocate(rdum1(10000), rdum2(10000), rdum3(10000) )
      allocate(sn2(sngi,nloc),snlx_orig(sngi,ndim,snloc),sweight(sngi) )
      allocate(suf_n(sngi,nloc),suf_nlx(sngi,ndim,nloc),suf_weight(sngi) )

      ipoly=npoly
      IQADRA=IPOLY+1

      sndim=ndim-1
      if(ele_type < is_triangle_or_tet) then ! not triangle...
        call get_shape_funs(ngi,nloc,ndim,  &
                 weight,n,nlx, ipoly,iqadra, &
                 snloc, sngi, sndim, sweight,sn_orig,snlx_orig, .true.   )

      else
        ! allocate(rdum1(10000),rdum2(10000),rdum3(10000))
        ! gives a surface triangle with time slab in surface integral.
        ! Amin- tured .true. to .false.
        call get_shape_funs(ngi,nloc,ndim,  &
                  weight,n,nlx, ipoly,iqadra, &
                  snloc, sngi, sndim, sweight,sn_orig,snlx_orig, .false.   )
       call unstr_tri_surface_pointers_sn(nface, sngi, snloc, nloc, ndim, sndim, totele, n_s_list_no, &
                      sn_orig, snlx_orig, face_sn, face_snlx, face_sn2, face_list_no)
        ! return a surface tet...
        ! Amin uncomment the 6 lines below
        ! suf_ngi=NGI/IQADRA
        ! suf_ndim=ndim-1
        ! suf_nloc=NLOC/(IPOLY+1)
        ! call get_shape_funs(suf_ngi,suf_nloc,suf_ndim,  &
        !           suf_weight,suf_n,suf_nlx, ipoly,iqadra, &
        !           sngi, sndim, rdum1,rdum2,rdum3, .false.   )

      endif
  end subroutine get_shape_funs_with_faces_unstr



  ! ! !this subroutine is for untructured triangular meshes created by Gmsh and read in by Msh2Tri.F90
  subroutine unstr_tri_surface_pointers_sn(nface, sngi, snloc, nloc, ndim, sndim, totele, n_s_list_no, &
                   sn_orig, snlx_orig, face_sn, face_snlx, face_sn2, face_list_no, neigh_orientation)
      ! ********************************************************************************************************
      ! calculate face_list_no( iface, ele), face_sn(:,:,iface); face_snlx(:,:,:,iface), face_sn2(:,:,s_list_no)
      ! ********************************************************************************************************
      ! ****integers:
      ! nface=no of faces of each elemenet.
      ! sngi=no of surface quadrature points of the faces - this is set to the max no of all faces.
      ! snloc=no of nodes on a surface element.
      ! nloc=no of nodes within a volume element.
      ! ndim=no of dimensions - including time possibly.
      ! sndim=ndim-1 dimensions of the surface elements.
      ! totele=no of elements.
      ! n_s_list_no= no of different oriantations for the surface element numbering.
      !
      ! ****original surface element shape functions:
      ! sn_orig(sgi,siloc) = shape function at quadrature point sgi and surface node siloc
      ! snlx_orig(sgi,1,siloc) = shape function derivative (only one in 2D) at quadrature point sgi and surface node siloc
      !
      ! ****new surface element shape functions designed to be highly efficient:
      ! face_sn(sgi,iloc,iface) = shape function at quadrature point sgi and volume node iloc for face iface of element.
      ! face_snlx(sgi,1,iloc,iface) = shape function derivative at quadrature point sgi and volume node iloc for face iface of element.
      ! face_sn2(sgi,iloc,n_s_list_no) = shape function at quadrature point sgi but on the otherside of face and volume node iloc for face iface of element.
      ! This works from: sn2 =face_sn2(:,:,s_list_no) in which  s_list_no = face_list_no( iface, ele) .
      implicit none
      integer, intent(in) :: nface, sngi, snloc, nloc, ndim, sndim, totele, n_s_list_no
      real, intent(in) :: sn_orig(sngi,snloc), snlx_orig(sngi,sndim,snloc)
      integer, optional, intent(in) :: neigh_orientation
      real, intent(out) :: face_sn(sngi,nloc,nface), face_snlx(sngi,sndim,nloc,nface), face_sn2(sngi,nloc,n_s_list_no)
      integer, intent(out) :: face_list_no(nface,totele)
      ! local variables...
      integer iface,lnod1,lnod2,i_s_list_no

      face_sn(:,:,:)=0.0
      face_snlx(:,:,:,:)=0.0
      face_sn2(:,:,:)=0.0

      iface=1
      lnod1=1
      lnod2=3
      face_sn(:,lnod1,iface)=sn_orig(:,1)
      face_sn(:,lnod2,iface)=sn_orig(:,2)
      face_snlx(:,1,lnod1,iface)=snlx_orig(:,1,1)
      face_snlx(:,1,lnod2,iface)=snlx_orig(:,1,2)
      i_s_list_no = iface
      ! face_list_no(iface,:) = i_s_list_no


      iface=2
      lnod1=2
      lnod2=1
      face_sn(:,lnod1,iface)=sn_orig(:,1)
      face_sn(:,lnod2,iface)=sn_orig(:,2)
      face_snlx(:,1,lnod1,iface)=snlx_orig(:,1,1)
      face_snlx(:,1,lnod2,iface)=snlx_orig(:,1,2)
      i_s_list_no = iface
      ! face_list_no(iface,:) = i_s_list_no

      iface=3
      lnod1=3
      lnod2=2
      face_sn(:,lnod1,iface)=sn_orig(:,1)
      face_sn(:,lnod2,iface)=sn_orig(:,2)
      face_snlx(:,1,lnod1,iface)=snlx_orig(:,1,1)
      face_snlx(:,1,lnod2,iface)=snlx_orig(:,1,2)
      i_s_list_no = iface

  end subroutine unstr_tri_surface_pointers_sn


  subroutine get_unstr_sn2(face_sn2,nface,iface,Nside,sngi,nloc,snloc,sn_orig,snlx_orig)
    ! @AMIN :: I commented case section as all parts are identical
    !@AMIN :: I can delete if statement and use face_sn2(:,locn,iface,Nside)
    implicit none
    integer, intent(in):: iface, Nside, sngi, snloc, nloc, nface
    real, intent(in):: sn_orig(sngi,snloc), snlx_orig(sngi,1,snloc)
    real, intent(inout):: face_sn2(sngi,nloc,nface)
    integer:: locn1, locn2

    face_sn2=0.0
    ! select case(iface)
      ! case(1)
        if ( Nside==1 ) then
          locn1=3
          locn2=1
          face_sn2(:,locn1,iface) = sn_orig(:,1)
          face_sn2(:,locn2,iface) = sn_orig(:,2)

        else if ( Nside==2 ) then
          locn1=1
          locn2=2
          face_sn2(:,locn1,iface)=sn_orig(:,1)
          face_sn2(:,locn2,iface)=sn_orig(:,2)

        else if ( Nside==3 ) then
          locn1=2
          locn2=3
          face_sn2(:,locn1,iface)=sn_orig(:,1)
          face_sn2(:,locn2,iface)=sn_orig(:,2)
        end if

  end subroutine get_unstr_sn2


  ! subroutine get_semi_sn2(meshlist, un_ele, mface, sn_orig, sn2)
  !   ! external vbl
  !   type(Mesh), dimension(:), INTENT(IN) :: meshList
  !   integer, INTENT(IN) :: un_ele, mface
  !   real, INTENT(IN) :: sn_orig(:,:)
  !
  !   ! local vbls
  !   real, INTENT(INOUT) :: sn2(:,:)
  !
  !   if (meshList(un_ele)%Dir(mface)) then
  !     sn2(:,meshList(un_ele)%S_nodes(2,mface))=sn_orig(:,1)
  !     sn2(:,meshList(un_ele)%S_nodes(1,mface))=sn_orig(:,2)
  !   else
  !     sn2(:,meshList(un_ele)%S_nodes(1,mface))=sn_orig(:,1)
  !     sn2(:,meshList(un_ele)%S_nodes(2,mface))=sn_orig(:,2)
  !   end if
  !
  !
  ! end subroutine get_semi_sn2


  subroutine get_sn2_implicit(sn, snlx, sn2,&
                    ele22, iface, ipos, irow, meshlist, un_ele, sn_orig, face_sn2_str,face_snlx_str,face_sn_str)
    implicit none
    !external vbls
    integer, intent(in) :: ele22, iface, ipos, irow, un_ele
    real, intent(in) :: sn_orig(:,:), face_sn2_str(:,:,:), face_snlx_str(:,:,:,:), face_sn_str(:,:,:)
    type(Mesh), dimension(:), INTENT(IN) :: meshList
    real, intent(inout) :: sn(:,:), snlx(:,:,:), sn2(:,:)
    ! local vbl
    integer :: mface, sp

    sn   = face_sn_str(:,:,iface)
    snlx = face_snlx_str(:,:,:,iface)

    sn2=0
    if ( ele22==0 ) then
          if ( iface==1 ) then
            sp = ipos/2+1 ! position of ele along the un_iface
            mface = 1 ! un_iface number which str_ele is located on
          elseif ( iface==2 ) then
            sp = irow
            mface = 3
          elseif (iface==3 ) then
            sp = irow
            mface = 2
          end if
          call get_semi_sn2_implicit(meshlist, un_ele, mface, sn_orig, sn2)
    else if ( ele22 /= 0 ) then
      sn2  = face_sn2_str(:,:,iface)
    end if ! end if ele22=0

  end subroutine get_sn2_implicit


  subroutine get_semi_sn2_implicit(meshlist, un_ele, mface, sn_orig, sn2)
    ! external vbl
    type(Mesh), dimension(:), INTENT(IN) :: meshList
    integer, INTENT(IN) :: un_ele, mface
    real, INTENT(IN) :: sn_orig(:,:)
! Amin it seems that mface and Nside are vise versa. I mean mface=3 means Nside=3 I am not sure why
    ! local vbls
    real, INTENT(INOUT) :: sn2(:,:)

          if ( mface==1 ) then
            sn2(:,meshList(un_ele)%S_nodes(:,mface))=sn_orig(:,:)
            ! sn2(:,meshList(un_ele)%S_nodes(1,mface))=sn_orig(:,1)
            ! sn2(:,meshList(un_ele)%S_nodes(2,mface))=sn_orig(:,2)
          else
            sn2(:,meshList(un_ele)%S_nodes(1,mface))=sn_orig(:,2)
            sn2(:,meshList(un_ele)%S_nodes(2,mface))=sn_orig(:,1)
          end if

  end subroutine get_semi_sn2_implicit


  !brief:: this subroutine gives stencil of un_str ele
  subroutine get_un_ele_mass_stiff_diffvol(mass_stcl,stiff_stcl,diff_vol_stcl,&
                                    n,nx,detwei,k,nloc,ngi,ndim,dt,ml_ele)
    implicit none
    !global vbl
    integer , intent(in) :: nloc, ngi, ndim
    REAL, intent( in ) :: N(ngi,nloc), nx(ngi,ndim,nloc), dt
    REAL, pointer, intent( in ) :: detwei(:)
    real, intent(in) :: k
    real, dimension(:,:), intent(inout) :: mass_stcl
    real, dimension(:), intent(inout) :: ml_ele
    ! real, dimension(:,:,:), intent(inout) :: stiff_ele !,diff_vol
    real, dimension(:,:,:,:), intent(inout) :: stiff_stcl !(nloc,ngi,idim, nloc) ! the last nloc actually referse to the final tnew_loc(iloc)
    real, intent(inout) :: diff_vol_stcl(:,:,:,:) !(ngi,ndim,nloc)

    !local vbl
    integer :: idim, iloc, jloc, gi
    ! real, dimension(ngi,nloc) :: mass_ele

    ! ! in stiff_stcl locations of gi and idim might not make sense but they are driven from full expasion
    ! ! of the main code and the final results are correct
    do jloc=1,nloc
      ml_ele(jloc)=sum(n(:,jloc)*detwei(:))
      do iloc=1,nloc
        mass_stcl(iloc,jloc) = sum(n(:,iloc)*detwei(:)*n(:,jloc))
        do idim=1,ndim
          do gi=1,ngi
            stiff_stcl(iloc,gi,idim,jloc) = nx(gi,idim,jloc)*detwei(gi)*n(gi,iloc)
            diff_vol_stcl(gi,idim,iloc,jloc) = k * nx(gi,idim,iloc) * detwei(gi)*nx(gi,idim,jloc)
          end do
        end do
      end do
    end do

    !----------- here is old version of the code which also works fine ----------
    ! stiff_ele = 0.0
    ! mass_ele = 0.0
    ! diff_vol=0.0
    ! do iloc=1,nloc
      ! do gi=1,ngi
        ! ml_ele(iloc)=sum(n(:,iloc)*detwei(:))/dt ! lumped mass matrix in element (use later for solver).
        ! do idim=1,ndim
          ! lumped stiff matrix
          ! stiff_ele(gi,idim, iloc) =  nx(gi,idim,iloc)*detwei(gi)
          ! diff_vol(gi,idim,iloc) = k * nx(gi,idim,iloc) * detwei(gi)
        ! end do
        ! adding mass matrix
        ! mass_ele(gi,iloc) = n(gi,iloc)*detwei(gi)/dt

      ! end do ! ngi
    ! end do ! iloc

    !---------------------- get the diff_vol stencils -------------------------
    ! do iloc=1,nloc
    !   do gi=1,ngi
    !     do idim=1,ndim
    !       diff_vol_stcl(gi,idim,iloc) = k * nx(gi,idim,iloc) * detwei(gi)*n(gi,idim,iloc)
    !     end do
    !   end do ! ngi
    ! end do ! iloc



    ! ----------- Here is the expansion of do loops above ----------------------
    !---------------------- get the Mass_ele stencils --------------------------
    ! sum of this times by t_loc(:) gives mass_ele(1)
    ! mass_stcl(1,1) = sum(mass_ele(:,1)*n(:,1))
    ! mass_stcl(1,2) = sum(mass_ele(:,1)*n(:,2))
    ! mass_stcl(1,3) = sum(mass_ele(:,1)*n(:,3))
    !
    ! ! sum of this times by t_loc(:) gives mass_ele(2)
    ! mass_stcl(2,1) = sum(mass_ele(:,2)*n(:,1))
    ! mass_stcl(2,2) = sum(mass_ele(:,2)*n(:,2))
    ! mass_stcl(2,3) = sum(mass_ele(:,2)*n(:,3))
    !
    ! ! sum of this times by t_loc(:) gives mass_ele(3)
    ! mass_stcl(3,1) = sum(mass_ele(:,3)*n(:,1))
    ! mass_stcl(3,2) = sum(mass_ele(:,3)*n(:,2))
    ! mass_stcl(3,3) = sum(mass_ele(:,3)*n(:,3))
    !---------------------- get the stiff_ele stencils -------------------------
    ! stiff_stcl(nloc,ngi,idim, nloc) ! the last nloc actually referse to the final tnew_loc(iloc)
    ! below is the expansion of the loop above for stiff_ele stencil
    ! stiff_stcl(1,1,1,1) = nx(1,1,1)*detwei(1)*n(1,1)
    ! stiff_stcl(1,2,1,1) = nx(2,1,1)*detwei(2)*n(2,1)
    ! stiff_stcl(1,3,1,1) = nx(3,1,1)*detwei(3)*n(3,1)
    ! stiff_stcl(1,1,2,1) = nx(1,2,1)*detwei(1)*n(1,1)
    ! stiff_stcl(1,2,2,1) = nx(2,2,1)*detwei(2)*n(2,1)
    ! stiff_stcl(1,3,2,1) = nx(3,2,1)*detwei(3)*n(3,1)
    ! !
    ! stiff_stcl(2,1,1,1) = nx(1,1,1)*detwei(1)*n(1,2)
    ! stiff_stcl(2,2,1,1) = nx(2,1,1)*detwei(2)*n(2,2)
    ! stiff_stcl(2,3,1,1) = nx(3,1,1)*detwei(3)*n(3,2)
    ! stiff_stcl(2,1,2,1) = nx(1,2,1)*detwei(1)*n(1,2)
    ! stiff_stcl(2,2,2,1) = nx(2,2,1)*detwei(2)*n(2,2)
    ! stiff_stcl(2,3,2,1) = nx(3,2,1)*detwei(3)*n(3,2)
    ! !
    ! stiff_stcl(3,1,1,1) = nx(1,1,1)*detwei(1)*n(1,3)
    ! stiff_stcl(3,2,1,1) = nx(2,1,1)*detwei(2)*n(2,3)
    ! stiff_stcl(3,3,1,1) = nx(3,1,1)*detwei(3)*n(3,3)
    ! stiff_stcl(3,1,2,1) = nx(1,2,1)*detwei(1)*n(1,3)
    ! stiff_stcl(3,2,2,1) = nx(2,2,1)*detwei(2)*n(2,3)
    ! stiff_stcl(3,3,2,1) = nx(3,2,1)*detwei(3)*n(3,3)


  end subroutine get_un_ele_mass_stiff_diffvol



end module ShapFun_unstruc
