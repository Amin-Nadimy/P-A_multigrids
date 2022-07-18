module subs
contains


  ! This sub calculates the local corrds L1, L2, L3, L4 and
  ! weights at the quadrature points.
  SUBROUTINE TRIQUAold(L1, L2, L3, L4, WEIGHT, D3,NGI)
    ! This sub calculates the local corrds L1, L2, L3, L4 and
    ! weights at the quadrature points.
    ! If D3 it does this for 3Dtetrahedra elements else
    ! triangular elements.
    IMPLICIT NONE
    INTEGER , intent(in):: NGI
    LOGICAL , intent(in) :: D3
    REAL , dimension(ngi) , intent(inout) ::L1, L2, L3, L4, WEIGHT
    ! Local variables...
    REAL :: ALPHA,BETA
    REAL :: ALPHA1,BETA1
    REAL :: ALPHA2,BETA2
    real :: rsum
    INTEGER I

    IF(D3) THEN
      ! this is for a tetrahedra element...
      ! This is for one point.
      IF(NGI.EQ.1) THEN
        ! Degree of precision is 1
        DO I=1,NGI
          L1(I)=0.25
          L2(I)=0.25
          L3(I)=0.25
          L4(I)=0.25
          WEIGHT(I)=1.0
        END DO
      ENDIF

      IF(NGI.EQ.4) THEN
        ! Degree of precision is 2
        ALPHA=0.58541020
        BETA=0.13819660
        DO I=1,NGI
          L1(I)=BETA
          L2(I)=BETA
          L3(I)=BETA
          L4(I)=BETA
          WEIGHT(I)=0.25
        END DO
        L1(1)=ALPHA
        L2(2)=ALPHA
        L3(3)=ALPHA
        L4(4)=ALPHA
      ENDIF

      IF(NGI.EQ.5) THEN
        ! Degree of precision is 3
        L1(1)=0.25
        L2(1)=0.25
        L3(1)=0.25
        L4(1)=0.25
        WEIGHT(1)=-4./5.

        DO I=2,NGI
          L1(I)=1./6.
          L2(I)=1./6.
          L3(I)=1./6.
          L4(I)=1./6.
          WEIGHT(I)=9./20.
        END DO
        L1(2)=0.5
        L2(3)=0.5
        L3(4)=0.5
        L4(5)=0.5
      ENDIF

      IF(NGI.EQ.11) THEN
        ! Degree of precision is 4
        ALPHA=(1.+SQRT(5./14.))/4.0
        BETA =(1.-SQRT(5./14.))/4.0
        I=1
        L1(I)=0.25
        L2(I)=0.25
        L3(I)=0.25
        WEIGHT(I)=-6.*74.0/5625.0
        DO I=2,5
          L1(I)=1./14.
          L2(I)=1./14.
          L3(I)=1./14.
          WEIGHT(I)=6.*343./45000.
        END DO
        L1(2)=11./14.
        L2(3)=11./14.
        L3(4)=11./14.
        DO I=6,11
          L1(I)=ALPHA
          L2(I)=ALPHA
          L3(I)=ALPHA
          WEIGHT(I)=6.*56.0/2250.0
        END DO
        L3(6)=BETA
        L2(7)=BETA
        L2(8)=BETA
        L3(8)=BETA
        L1(9)=BETA
        L1(10)=BETA
        L3(10)=BETA
        L1(11)=BETA
        L2(11)=BETA
        ! ENDOF IF(NGI.EQ.11) THEN...
      ENDIF

      ! if (NGI == 15) then!Fith order quadrature
      !   ! Degree of precision is 5
      !   L1=(/0.2500000000000000, 0.0000000000000000, 0.3333333333333333, 0.3333333333333333, 0.3333333333333333, &
      !          0.7272727272727273, 0.0909090909090909, 0.0909090909090909, 0.0909090909090909, 0.4334498464263357, &
      !          0.0665501535736643, 0.0665501535736643, 0.0665501535736643, 0.4334498464263357, 0.4334498464263357/)
      !   L2=(/0.2500000000000000, 0.3333333333333333, 0.3333333333333333, 0.3333333333333333, 0.0000000000000000, &
      !          0.0909090909090909, 0.0909090909090909, 0.0909090909090909, 0.7272727272727273, 0.0665501535736643, &
      !          0.4334498464263357, 0.0665501535736643, 0.4334498464263357, 0.0665501535736643, 0.4334498464263357/)
      !   L3=(/0.2500000000000000, 0.3333333333333333, 0.3333333333333333, 0.0000000000000000, 0.3333333333333333, &
      !          0.0909090909090909, 0.0909090909090909, 0.7272727272727273, 0.0909090909090909, 0.0665501535736643, &
      !          0.0665501535736643, 0.4334498464263357, 0.4334498464263357, 0.4334498464263357, 0.0665501535736643/)
      !   !We divide the weights later by 6
      !   WEIGHT=(/0.1817020685825351, 0.0361607142857143, 0.0361607142857143, 0.0361607142857143, 0.0361607142857143, &
      !          0.0698714945161738, 0.0698714945161738, 0.0698714945161738, 0.0698714945161738, 0.0656948493683187, &
      !          0.0656948493683187, 0.0656948493683187, 0.0656948493683187, 0.0656948493683187, 0.0656948493683187/)
      ! end if

      ! if (NGI == 45) then!Eighth order quadrature, for bubble shape functions or P3
      !   !Obtained from: https://people.sc.fsu.edu/~jburkardt/datasets/quadrature_rules_tet/quadrature_rules_tet.html
      !   !Here to get triangle quadrature sets:
      !   !https://people.sc.fsu.edu/~jburkardt/datasets/quadrature_rules_tri/quadrature_rules_tri.html
      !   ! Degree of precision is 8
      !   L1=(/0.2500000000000000,0.6175871903000830,0.1274709365666390,0.1274709365666390,0.1274709365666390,0.9037635088221031,&
      !     0.0320788303926323,0.0320788303926323,0.0320788303926323,0.4502229043567190,0.0497770956432810,0.0497770956432810,&
      !     0.0497770956432810,0.4502229043567190,0.4502229043567190,0.3162695526014501,0.1837304473985499,0.1837304473985499,&
      !     0.1837304473985499,0.3162695526014501,0.3162695526014501,0.0229177878448171,0.2319010893971509,0.2319010893971509,&
      !     0.5132800333608811,0.2319010893971509,0.2319010893971509,0.2319010893971509,0.0229177878448171,0.5132800333608811,&
      !     0.2319010893971509,0.0229177878448171,0.5132800333608811,0.7303134278075384,0.0379700484718286,0.0379700484718286,&
      !     0.1937464752488044,0.0379700484718286,0.0379700484718286,0.0379700484718286,0.7303134278075384,0.1937464752488044,&
      !     0.0379700484718286,0.7303134278075384,0.1937464752488044/)
      !   L2=(/0.2500000000000000,0.1274709365666390,0.1274709365666390,0.1274709365666390,0.6175871903000830,0.0320788303926323,&
      !     0.0320788303926323,0.0320788303926323,0.9037635088221031,0.0497770956432810,0.4502229043567190,0.0497770956432810,&
      !     0.4502229043567190,0.0497770956432810,0.4502229043567190,0.1837304473985499,0.3162695526014501,0.1837304473985499,&
      !     0.3162695526014501,0.1837304473985499,0.3162695526014501,0.2319010893971509,0.0229177878448171,0.2319010893971509,&
      !     0.2319010893971509,0.5132800333608811,0.2319010893971509,0.0229177878448171,0.5132800333608811,0.2319010893971509,&
      !     0.5132800333608811,0.2319010893971509,0.0229177878448171,0.0379700484718286,0.7303134278075384,0.0379700484718286,&
      !     0.0379700484718286,0.1937464752488044,0.0379700484718286,0.7303134278075384,0.1937464752488044,0.0379700484718286,&
      !     0.1937464752488044,0.0379700484718286,0.7303134278075384/)
      !   L3=(/0.2500000000000000,0.1274709365666390,0.1274709365666390,0.6175871903000830,0.1274709365666390,0.0320788303926323,&
      !     0.0320788303926323,0.9037635088221031,0.0320788303926323,0.0497770956432810,0.0497770956432810,0.4502229043567190,&
      !     0.4502229043567190,0.4502229043567190,0.0497770956432810,0.1837304473985499,0.1837304473985499,0.3162695526014501,&
      !     0.3162695526014501,0.3162695526014501,0.1837304473985499,0.2319010893971509,0.2319010893971509,0.0229177878448171,&
      !     0.2319010893971509,0.2319010893971509,0.5132800333608811,0.5132800333608811,0.2319010893971509,0.0229177878448171,&
      !     0.0229177878448171,0.5132800333608811,0.2319010893971509,0.0379700484718286,0.0379700484718286,0.7303134278075384,&
      !     0.0379700484718286,0.0379700484718286,0.1937464752488044,0.1937464752488044,0.0379700484718286,0.7303134278075384,&
      !     0.7303134278075384,0.1937464752488044,0.0379700484718286/)
      !   !We divide the weights later by 6
      !   WEIGHT=(/-0.2359620398477557,0.0244878963560562,0.0244878963560562,0.0244878963560562,0.0244878963560562,0.0039485206398261,&
      !     0.0039485206398261,0.0039485206398261,0.0039485206398261,0.0263055529507371,0.0263055529507371,0.0263055529507371,&
      !     0.0263055529507371,0.0263055529507371,0.0263055529507371,0.0829803830550589,0.0829803830550589,0.0829803830550589,&
      !     0.0829803830550589,0.0829803830550589,0.0829803830550589,0.0254426245481023,0.0254426245481023,0.0254426245481023,&
      !     0.0254426245481023,0.0254426245481023,0.0254426245481023,0.0254426245481023,0.0254426245481023,0.0254426245481023,&
      !     0.0254426245481023,0.0254426245481023,0.0254426245481023,0.0134324384376852,0.0134324384376852,0.0134324384376852,&
      !     0.0134324384376852,0.0134324384376852,0.0134324384376852,0.0134324384376852,0.0134324384376852,0.0134324384376852,&
      !     0.0134324384376852,0.0134324384376852,0.0134324384376852/)
      ! end if

      DO I=1,NGI
        L4(I)=1.0-L1(I)-L2(I)-L3(I)
      END DO

      ! Now multiply by 1/6. to get weigts correct...
      DO I=1,NGI
        WEIGHT(I)=WEIGHT(I)/6.
      END DO
      ! ENDOF IF(D3) THEN...
    ENDIF

    IF(.NOT.D3) THEN
      ! 2-D TRAINGULAR ELEMENTS...
      IF(NGI.EQ.1) THEN
        ! LINEAR
        I=1
        L1(I)=1./3.
        L2(I)=1./3.
        WEIGHT(I)=1.0
      ENDIF

      IF(NGI.EQ.3) THEN
        ! QUADRASTIC
        DO I=1,NGI
          L1(I)=0.5
          L2(I)=0.5
          WEIGHT(I)=1.0/3.0
        END DO
        L1(2)=0.0
        L2(3)=0.0
      ENDIF

      IF(NGI.EQ.4) THEN
        ! CUBIC
        I=1
        L1(I)=1./3.
        L2(I)=1./3.
        WEIGHT(I)=-27./48.
        DO I=2,NGI
          L1(I)=0.2
          L2(I)=0.2
          WEIGHT(I)=25./48.
        END DO
        L1(1)=0.6
        L2(2)=0.6
      ENDIF

      IF(NGI.EQ.7) THEN
        ! QUNTIC
        ALPHA1=0.0597158717
        BETA1 =0.4701420641
        ALPHA2=0.7974269853
        BETA2 =0.1012865073
        I=1
        L1(I)=1./3.
        L2(I)=1./3.
        WEIGHT(I)=0.225
        DO I=2,4
          L1(I)=BETA1
          L2(I)=BETA1
          WEIGHT(I)=0.1323941527
        END DO
        L1(2)=ALPHA1
        L2(4)=ALPHA1
        DO I=5,7
          L1(I)=BETA2
          L2(I)=BETA2
          WEIGHT(I)=0.1259391805
        END DO
        L1(5)=ALPHA2
        L2(6)=ALPHA2
        ! ENDOF IF(NGI.EQ.7) THEN...
      ENDIF

      IF(NGI.EQ.14) THEN
        ! 5th order quadrature set...
        L1(1) = 6.943184420297371E-002
        L1(2) = 6.943184420297371E-002
        L1(3) = 6.943184420297371E-002
        L1(4) = 6.943184420297371E-002
        L1(5) = 6.943184420297371E-002
        L1(6) = 0.330009478207572
        L1(7) = 0.330009478207572
        L1(8) = 0.330009478207572
        L1(9) = 0.330009478207572
        L1(10) = 0.669990521792428
        L1(11) = 0.669990521792428
        L1(12) = 0.669990521792428
        L1(13) = 0.930568155797026
        L1(14) = 0.930568155797026
        ! local coord 1:
        L2(1) = 4.365302387072518E-002
        L2(2) = 0.214742881469342
        L2(3) = 0.465284077898513
        L2(4) = 0.715825274327684
        L2(5) = 0.886915131926301
        L2(6) = 4.651867752656094E-002
        L2(7) = 0.221103222500738
        L2(8) = 0.448887299291690
        L2(9) = 0.623471844265867
        L2(10) = 3.719261778493340E-002
        L2(11) = 0.165004739103786
        L2(12) = 0.292816860422638
        L2(13) = 1.467267513102734E-002
        L2(14) = 5.475916907194637E-002
        ! local coord 2:
        WEIGHT(1) = 1.917346464706755E-002
        WEIGHT(2) = 3.873334126144628E-002
        WEIGHT(3) = 4.603770904527855E-002
        WEIGHT(4) = 3.873334126144628E-002
        WEIGHT(5) = 1.917346464706755E-002
        WEIGHT(6) = 3.799714764789616E-002
        WEIGHT(7) = 7.123562049953998E-002
        WEIGHT(8) = 7.123562049953998E-002
        WEIGHT(9) = 3.799714764789616E-002
        WEIGHT(10) = 2.989084475992800E-002
        WEIGHT(11) = 4.782535161588505E-002
        WEIGHT(12) = 2.989084475992800E-002
        WEIGHT(13) = 6.038050853208200E-003
        WEIGHT(14) = 6.038050853208200E-003
        rsum=SUM(WEIGHT(1:NGI))
        WEIGHT(1:NGI)=WEIGHT(1:NGI)/RSUM
        ! ENDOF IF(NGI.EQ.14) THEN...
      ENDIF


      DO I=1,NGI
        L3(I)=1.0-L1(I)-L2(I)
      END DO
      ! ENDOF IF(.NOT.D3) THEN...
    ENDIF
    RETURN
  END subroutine TRIQUAold


  subroutine get_shape_funs_spec(n, nlx, nlx_lxx, nlxx, weight, nlx_nod, &
            nloc, sngi, ngi, ndim, nface,max_face_list_no, face_sn, face_sn2, face_snlx, face_sweigh, &
            npoly, ele_type )
    implicit none
    ! nloc=no of local nodes per element
    ! ngi=no of quadrature points
    ! ndim=no of dimensions.
    integer, intent(in) :: nloc, sngi, ngi, ndim, nface, max_face_list_no
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
    real, intent(inout) :: nlx_nod(nloc,ndim,nloc)
    real, intent(inout) :: face_sn(sngi,nloc,nface), face_sn2(sngi,nloc,max_face_list_no)
    real, intent(inout) :: face_snlx(sngi,ndim,nloc,nface), face_sweigh(sngi,nface)
    ! npoly=order of polynomial in Cartesian space; ele_type=type of element including order of poly.
    integer, intent(in) :: npoly,ele_type

    ! form the shape functions...
    call get_shape_funs_with_faces(n, nlx, weight,  &
               nloc, sngi, ngi, ndim, nface,max_face_list_no, face_sn, face_sn2, face_snlx, face_sweigh, &
               npoly,ele_type)

    ! Calculate high order derivatives of shape functions nlx_lxx, nlxx
    ! and also calculate node-wise deriavtives of shape functions nlx_nod.
    ! Amin for higher order uncomment get_high_order_shape_funs below
    ! call get_high_order_shape_funs(n, nlx, nlx_lxx, nlxx, weight, nlx_nod, nloc, ngi, ndim, &
    !              npoly, ele_type)
  end subroutine get_shape_funs_spec


  subroutine get_shape_funs_with_faces(n, nlx, weight,  &
                 nloc, sngi, ngi, ndim, nface,max_face_list_no, face_sn, face_sn2, face_snlx, face_sweigh, &
                 npoly,ele_type)
    implicit none
    ! nloc=no of local nodes per element
    ! ngi=no of quadrature points
    ! ndim=no of dimensions.
    ! ele_type= element type
    integer, intent(in) :: nloc, sngi, ngi, ndim, nface, max_face_list_no
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
    real, intent(inout) :: face_sn(sngi,nloc,nface), face_sn2(sngi,nloc,max_face_list_no)
    real, intent(inout) ::  face_snlx(sngi,ndim,nloc,nface), face_sweigh(sngi,nface)
    ! npoly=order of polynomial in Cartesian space; ele_type=type of element including order of poly.
    integer, intent(in) :: npoly, ele_type
    ! local variables...
    integer is_triangle_or_tet
    parameter(is_triangle_or_tet=100)
    integer sndim, suf_ngi, suf_ndim, suf_nloc, ipoly, IQADRA
    real, allocatable :: rdum1(:), rdum2(:), rdum3(:)
    real, allocatable :: sn(:,:),sn2(:,:),snlx(:,:,:),sweight(:)
    real, allocatable :: suf_n(:,:),suf_nlx(:,:,:),suf_weight(:)

    ! allocate memory...
    allocate(rdum1(10000), rdum2(10000), rdum3(10000) )
    allocate(sn(sngi,nloc),sn2(sngi,nloc),snlx(sngi,ndim,nloc),sweight(sngi) )
    allocate(suf_n(sngi,nloc),suf_nlx(sngi,ndim,nloc),suf_weight(sngi) )

    ipoly=npoly
    IQADRA=IPOLY+1

    sndim=ndim-1
    if(ele_type < is_triangle_or_tet) then ! not triangle...
      call get_shape_funs(ngi,nloc,ndim,  &
               weight,n,nlx, ipoly,iqadra, &
               sngi, sndim, sweight,sn,snlx, .true.   )
      print*, 'if'

    else
      print*, 'else'
      ! allocate(rdum1(10000),rdum2(10000),rdum3(10000))
      ! gives a surface triangle with time slab in surface integral.
      ! Amin- tured .true. to .false.
      call get_shape_funs(ngi,nloc,ndim,  &
                weight,n,nlx, ipoly,iqadra, &
                sngi, sndim, sweight,sn,snlx, .false.   )
      ! return a surface tet...
      suf_ngi=NGI/IQADRA
      suf_ndim=ndim-1
      suf_nloc=NLOC/(IPOLY+1)
      call get_shape_funs(suf_ngi,suf_nloc,suf_ndim,  &
                suf_weight,suf_n,suf_nlx, ipoly,iqadra, &
                sngi, sndim, rdum1,rdum2,rdum3, .false.   )

    endif
  end subroutine get_shape_funs_with_faces


  ! form volume and surface shape functions
  subroutine get_shape_funs(ngi,nloc,ndim,  &
               weight,n,nlx, ipoly,iqadra, &
               sngi, sndim, sweight,sn,snlx, with_time_slab   )
    ! ***************************************
    ! form volume and surface shape functions
    ! ***************************************
    IMPLICIT NONE
    INTEGER, intent(in) :: sngi, NGI, NLOC, ndim, sndim
    INTEGER, intent(in) :: IPOLY,IQADRA
    logical, intent(in) :: with_time_slab
    REAL, intent(inout) ::  n(ngi,nloc)
    REAL, intent(inout) ::  nlx(ngi,ndim,nloc)
    REAL, intent(inout) ::  sn(sngi, nloc)
    REAL, intent(inout) ::  snlx(sngi, sndim, nloc)
    real, intent(inout) :: WEIGHT(ngi), sweight(sngi)
    ! local variables...
    integer mloc, NDNOD, snloc_temp
    logical square
    real, allocatable :: m(:)

    allocate(m(10000))
    mloc=1
    SQUARE = .false.

    IF(SQUARE) THEN ! Square in up to 4D
      ! volumes
      ! AMIN turn on  interface_SPECTR
      ! call interface_SPECTR(NGI,NLOC,WEIGHT,N,NLX, ndim, IPOLY,IQADRA  )
      ! surfaces...
      NDNOD =INT((NLOC**(1./real(ndim) ))+0.1)
      snloc_temp=NDNOD**(ndim-1)
      sweight=0.0; sN=0.0; sNLX=0.0

      !AMIN turn on interface_SPECTR
      ! call interface_SPECTR(sNGI,sNLOC_temp,sweight,sN,sNLX, ndim-1, IPOLY,IQADRA  )

    ELSE
      ! volume tet plus time slab...
      call triangles_tets_with_time( nloc,ngi, ndim, n,nlx, weight, ipoly,iqadra, with_time_slab)
      ! surfaces...
      call triangles_tets_with_time( nloc,sngi, ndim-1, sn,snlx, sweight, ipoly,iqadra, with_time_slab)
    ENDIF
  end subroutine get_shape_funs


  subroutine triangles_tets_with_time(nloc,ngi,ndim,n,nlx,weight,ipoly,iqadra,with_time_slab)
    ! ****************************************************************************************
    implicit none
    integer, intent(in) :: nloc,ngi,ndim, ipoly,iqadra
    real, intent(inout) :: n(ngi,nloc), nlx(ngi,ndim,nloc)
    real, intent(inout) :: weight(ngi)
    logical, intent(in) :: with_time_slab
    ! local variables...
    integer nloc_space, ngi_space, ndim_space,   nloc_t, ngi_t, ndim_t
    real, allocatable :: l1(:), l2(:), l3(:), l4(:)
    real, allocatable :: weight_space(:), n_space(:,:), nlx_space(:,:,:)
    real, allocatable :: weight_t(:), n_t(:,:), nlx_t(:,:,:)

    if(with_time_slab) then ! form a space-time slab.

      nloc_t=nloc/(ipoly+1)
      NGI_T=IQADRA
      !ngi_t=ngi/(ipoly+2)
      ndim_t = 1

      nloc_space=nloc/nloc_t
      ngi_space=ngi/ngi_t
      ndim_space = ndim-1
      allocate(l1(ngi_space), l2(ngi_space), l3(ngi_space), l4(ngi_space) )
      allocate(weight_space(ngi_space), n_space(ngi_space,nloc_space), nlx_space(ngi_space,ndim_space,nloc_space) )

      allocate(weight_t(ngi_t), n_t(ngi_t,nloc_t), nlx_t(ngi_t,ndim_t,nloc_t) )

      ! triangles or tetrahedra...
      call SHATRInew(L1, L2, L3, L4, WEIGHT_space, &
            NLOC_space,NGI_space,ndim_space,  n_space,nlx_space)
      ! ! extend into time domain...
      !AMIN turn on interface_SPECTR
      ! call interface_SPECTR(NGI_t,NLOC_t,WEIGHT_t,N_t,NLX_t, ndim_t, IPOLY,IQADRA  )

      ! ! combine-space time...
      !AMIN tuen on make_space_time_shape_funs
      ! call make_space_time_shape_funs(nloc_space,ngi_space, ndim_space, n_space,nlx_space, weight_space, &
      !                                  nloc_t,ngi_t, ndim_t, n_t, nlx_t, weight_t, &
      !                                  nloc,ngi, ndim, n, nlx, weight )
    else ! just a triangle or tet without time slab...
      ! triangles or tetrahedra...
      allocate(l1(ngi), l2(ngi), l3(ngi), l4(ngi) )
      ! Amin for d3=3 turn .false. to .true.
      call TRIQUAold(L1, L2, L3, L4, WEIGHT, .false.,NGI)
      call SHATRInew(L1, L2, L3, L4, weight, nloc,ngi,ndim,  n,nlx)
    endif
  end subroutine triangles_tets_with_time


    ! Interface to SHATRIold using the new style variables
  SUBROUTINE SHATRInew(L1, L2, L3, L4, WEIGHT, NLOC,NGI,ndim,  N,NLX_ALL)
    ! Interface to SHATRIold using the new style variables
    IMPLICIT NONE
    INTEGER , intent(in) :: NLOC,NGI,ndim
    REAL , dimension(ngi), intent(in) :: L1, L2, L3, L4
    REAL , dimension(ngi), intent(inout) :: WEIGHT
    REAL , dimension(ngi, nloc ), intent(inout) ::N
    real, dimension (ngi,ndim,nloc), intent(inout) :: NLX_ALL

    call SHATRIold(L1, L2, L3, L4, WEIGHT, .false., &
           NLOC,NGI,  &
           N,NLX_ALL(:,1,:),NLX_ALL(:,2,:),NLX_ALL(:,3,:))
  end subroutine SHATRInew


  ! Work out the shape functions and there derivatives...
  SUBROUTINE SHATRIold(L1, L2, L3, L4, WEIGHT, D3,NLOC,NGI,N,NLX,NLY,NLZ)
    ! Work out the shape functions and there derivatives...
    IMPLICIT NONE
    INTEGER , intent(in) :: NLOC,NGI
    LOGICAL , intent(in) :: D3
    REAL , dimension(ngi), intent(in) :: L1, L2, L3, L4
    REAL , dimension(ngi), intent(inout) :: WEIGHT
    REAL , dimension(nloc, ngi ), intent(inout) ::N,NLX,NLY,NLZ
    ! Local variables...
    INTEGER ::  GI

    IF(.NOT.D3) THEN
      ! Assume a triangle...

      IF(NLOC.EQ.1) THEN
        Loop_Gi_Nloc1: DO GI=1,NGI
          N(1,GI)=1.0
          NLX(1,GI)=0.0
          NLY(1,GI)=0.0
        end DO Loop_Gi_Nloc1
      ELSE IF((NLOC.EQ.3).OR.(NLOC.EQ.4)) THEN
        Loop_Gi_Nloc3_4: DO GI=1,NGI
          N(1,GI)=L1(GI)
          N(2,GI)=L2(GI)
          N(3,GI)=L3(GI)

          NLX(1,GI)=1.0
          NLX(2,GI)=0.0
          NLX(3,GI)=-1.0

          NLY(1,GI)=0.0
          NLY(2,GI)=1.0
          NLY(3,GI)=-1.0
          IF(NLOC.EQ.4) THEN
            ! Bubble function...
            !alpha == 1 behaves better than the correct value of 27. See Osman et al. 2019
            N(4,GI)  =1. * L1(GI)*L2(GI)*L3(GI)
            NLX(4,GI)=1. * L2(GI)*(1.-L2(GI))-2.*L1(GI)*L2(GI)
            NLY(4,GI)=1. * L1(GI)*(1.-L1(GI))-2.*L1(GI)*L2(GI)
          ENDIF
        end DO Loop_Gi_Nloc3_4
      ELSE IF((NLOC.EQ.6).OR.(NLOC.EQ.7)) THEN
        Loop_Gi_Nloc_6_7: DO GI=1,NGI
          N(1,GI)=(2.*L1(GI)-1.)*L1(GI)
          N(2,GI)=(2.*L2(GI)-1.)*L2(GI)
          N(3,GI)=(2.*L3(GI)-1.)*L3(GI)

          N(4,GI)=4.*L1(GI)*L2(GI)
          N(5,GI)=4.*L2(GI)*L3(GI)
          N(6,GI)=4.*L1(GI)*L3(GI)

          ! nb L1+L2+L3+L4=1
          ! x-derivative...
          NLX(1,GI)=4.*L1(GI)-1.
          NLX(2,GI)=0.
          NLX(3,GI)=-4.*(1.-L2(GI))+4.*L1(GI) + 1.

          NLX(4,GI)=4.*L2(GI)
          NLX(5,GI)=-4.*L2(GI)
          NLX(6,GI)=4.*(1.-L2(GI))-8.*L1(GI)

          ! y-derivative...
          NLY(1,GI)=0.
          NLY(2,GI)=4.*L2(GI)-1.0
          NLY(3,GI)=-4.*(1.-L1(GI))+4.*L2(GI) + 1.

          NLY(4,GI)=4.*L1(GI)
          NLY(5,GI)=4.*(1.-L1(GI))-8.*L2(GI)
          NLY(6,GI)=-4.*L1(GI)
          IF(NLOC.EQ.7) THEN
            ! Bubble function...
            N(7,GI)  =L1(GI)*L2(GI)*L3(GI)
            NLX(7,GI)=L2(GI)*(1.-L2(GI))-2.*L1(GI)*L2(GI)
            NLY(7,GI)=L1(GI)*(1.-L1(GI))-2.*L1(GI)*L2(GI)
          ENDIF
        END DO Loop_Gi_Nloc_6_7
        ! ENDOF IF(NLOC.EQ.6) THEN...
      ! ELSE IF(NLOC==10) THEN ! Cubic triangle...
        ! get the shape functions for a cubic triangle...
        ! call shape_triangle_cubic( l1, l2, l3, l4, weight, d3, &
        !          nloc, ngi, &
        !          n, nlx, nly, nlz )

      ELSE ! has not found the element shape functions
        stop 811
      ENDIF
      ! ENDOF IF(.NOT.D3) THEN
    ENDIF


    IF(D3) THEN
      ! Assume a tet...
      ! This is for 5 point quadrature.
      IF((NLOC.EQ.10).OR.(NLOC.EQ.11)) THEN
        Loop_Gi_Nloc_10_11: DO GI=1,NGI
          !ewrite(3,*)'gi,L1(GI),L2(GI),L3(GI),L4(GI):',gi,L1(GI),L2(GI),L3(GI),L4(GI)
          N(1,GI)=(2.*L1(GI)-1.)*L1(GI)
          N(3,GI)=(2.*L2(GI)-1.)*L2(GI)
          N(5,GI)=(2.*L3(GI)-1.)*L3(GI)
          N(10,GI)=(2.*L4(GI)-1.)*L4(GI)

          !if(L1(GI).gt.-1.93) ewrite(3,*)'gi,L1(GI), L2(GI), L3(GI), L4(GI),N(1,GI):', &
          !                            gi,L1(GI), L2(GI), L3(GI), L4(GI),N(1,GI)
          !
          !
          N(2,GI)=4.*L1(GI)*L2(GI)
          N(6,GI)=4.*L1(GI)*L3(GI)
          N(7,GI)=4.*L1(GI)*L4(GI)

          N(4,GI) =4.*L2(GI)*L3(GI)
          N(9,GI) =4.*L3(GI)*L4(GI)
          N(8,GI)=4.*L2(GI)*L4(GI)
          ! nb L1+L2+L3+L4=1
          ! x-derivative...
          NLX(1,GI)=4.*L1(GI)-1.
          NLX(3,GI)=0.
          NLX(5,GI)=0.
          NLX(10,GI)=-4.*(1.-L2(GI)-L3(GI))+4.*L1(GI) + 1.
          !if(L1(GI).gt.-1.93) ewrite(3,*)'Nlx(1,GI):', &
          !     Nlx(1,GI)

          NLX(2,GI)=4.*L2(GI)
          NLX(6,GI)=4.*L3(GI)
          NLX(7,GI)=4.*(L4(GI)-L1(GI))

          NLX(4,GI) =0.
          NLX(9,GI) =-4.*L3(GI)
          NLX(8,GI)=-4.*L2(GI)

          ! y-derivative...
          NLY(1,GI)=0.
          NLY(3,GI)=4.*L2(GI)-1.0
          NLY(5,GI)=0.
          NLY(10,GI)=-4.*(1.-L1(GI)-L3(GI))+4.*L2(GI) + 1.

          NLY(2,GI)=4.*L1(GI)
          NLY(6,GI)=0.
          NLY(7,GI)=-4.*L1(GI)

          NLY(4,GI) =4.*L3(GI)
          NLY(9,GI) =-4.*L3(GI)
          NLY(8,GI)=4.*(1-L1(GI)-L3(GI))-8.*L2(GI)

          ! z-derivative...
          NLZ(1,GI)=0.
          NLZ(3,GI)=0.
          NLZ(5,GI)=4.*L3(GI)-1.
          NLZ(10,GI)=-4.*(1.-L1(GI)-L2(GI))+4.*L3(GI) + 1.

          NLZ(2,GI)=0.
          NLZ(6,GI)=4.*L1(GI)
          NLZ(7,GI)=-4.*L1(GI)

          NLZ(4,GI) =4.*L2(GI)
          NLZ(9,GI) =4.*(1.-L1(GI)-L2(GI))-8.*L3(GI)
          NLZ(8,GI)=-4.*L2(GI)
          IF(NLOC.EQ.11) THEN
            ! Bubble function...
            N(11,GI)  =L1(GI)*L2(GI)*L3(GI)*L4(GI)
            NLX(11,GI)=L2(GI)*L3(GI)*(1.-L2(GI)-L3(GI))-2.*L1(GI)*L2(GI)*L3(GI)
            NLY(11,GI)=L1(GI)*L3(GI)*(1.-L1(GI)-L3(GI))-2.*L1(GI)*L2(GI)*L3(GI)
            NLZ(11,GI)=L1(GI)*L2(GI)*(1.-L1(GI)-L2(GI))-2.*L1(GI)*L2(GI)*L3(GI)
          ENDIF
        end DO Loop_Gi_Nloc_10_11
        ! ENDOF IF(NLOC.EQ.10) THEN...
      ENDIF

      IF((NLOC.EQ.4).OR.(NLOC.EQ.5)) THEN
        Loop_Gi_Nloc_4_5: DO GI=1,NGI
          N(1,GI)=L1(GI)
          N(2,GI)=L2(GI)
          N(3,GI)=L3(GI)
          N(4,GI)=L4(GI)

          NLX(1,GI)=1.0
          NLX(2,GI)=0
          NLX(3,GI)=0
          NLX(4,GI)=-1.0

          NLY(1,GI)=0.0
          NLY(2,GI)=1.0
          NLY(3,GI)=0.0
          NLY(4,GI)=-1.0

          NLZ(1,GI)=0.0
          NLZ(2,GI)=0.0
          NLZ(3,GI)=1.0
          NLZ(4,GI)=-1.0
          IF(NLOC.EQ.5) THEN
            ! Bubble function ...
            !alpha == 50 behaves better than the correct value of 256. See Osman et al. 2019
            N(5,GI)  = 50. * L1(GI)*L2(GI)*L3(GI)*L4(GI)
            NLX(5,GI)= 50. * L2(GI)*L3(GI)*(1.-L2(GI)-L3(GI))-2.*L1(GI)*L2(GI)*L3(GI)
            NLY(5,GI)= 50. * L1(GI)*L3(GI)*(1.-L1(GI)-L3(GI))-2.*L1(GI)*L2(GI)*L3(GI)
            NLZ(5,GI)= 50. * L1(GI)*L2(GI)*(1.-L1(GI)-L2(GI))-2.*L1(GI)*L2(GI)*L3(GI)
          ENDIF
        end DO Loop_Gi_Nloc_4_5
      ENDIF
      IF(NLOC.EQ.1) THEN
        Loop_Gi_Nloc_1: DO GI=1,NGI
          N(1,GI)=1.0
          NLX(1,GI)=0.0
          NLY(1,GI)=0.0
          NLZ(1,GI)=0.0
        end DO Loop_Gi_Nloc_1
      ENDIF
      ! ENDOF IF(D3) THEN...
    ENDIF
    RETURN
  END SUBROUTINE SHATRIold

end module subs





program shatritest
  use subs
  implicit none
  integer :: nloc, ngi, ndim, ipoly,iqadra, sngi, sndim, nface, ele_type,max_face_list_no
  integer:: npoly
  logical:: with_time_slab
  real, dimension( : ), allocatable:: l1, l2, l3, l4, weight
  logical :: d3
  real, dimension( :, : ), allocatable :: n,face_sweigh
  real, dimension(:,:,:), allocatable :: NLX_ALL,nlx_nod, nlx_lxx
  real, allocatable :: face_snlx(:,:,:,:), nlxx(:,:)
  real, allocatable :: sn(:,:),snlx(:,:,:),sweight(:), face_sn(:,:,:),face_sn2(:,:,:)

  ngi = 3; nloc = 3; ndim=2; npoly=1; iqadra=2; with_time_slab=.false.
  sngi=2; nface=3; ele_type=101; max_face_list_no=3
  ipoly=npoly

  allocate(l1(ngi), L2(ngi), L3(ngi), L4(ngi), weight(ngi))
  allocate(n(nloc,ngi), nlx_all(ngi,ndim,nloc))
  allocate(sn(sngi,nloc),snlx(sngi,sndim,nloc),sweight(sngi), face_sn(sngi,nloc,nface))
  allocate(face_sn2(sngi,nloc,max_face_list_no),face_snlx(sngi,ndim,nloc,nface),face_sweigh(sngi,nface))
  allocate(nlx_lxx(ngi,ndim,nloc), nlx_nod(nloc,ndim,nloc),nlxx(ngi,nloc))


  call TRIQUAold(L1, L2, L3, L4, WEIGHT, D3,NGI)
  ! call shatri( l1, l2, l3, l4, weight, d3, nloc, ngi, n, nlx, nly, nlz )
  call get_shape_funs_spec(n, nlx_all, nlx_lxx, nlxx, weight, nlx_nod, &
            nloc, sngi, ngi, ndim, nface,max_face_list_no, face_sn, face_sn2, face_snlx, face_sweigh, &
            npoly, ele_type )
print*, nlx_all
print*, n

  deallocate(n, nlx_all, L1, L2, L3 ,L4)
end program shatritest
