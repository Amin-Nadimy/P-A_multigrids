module ShapFun
  use Structures
  use Msh2Tri
  contains



    ! this subroutine calculates your unknown & velocity @ local surface and @ Gaussian Quadrature points
    ! it works only for semi-structured
    ! get_loc_sgi return its firs line of vbls
    subroutine get_loc_sgi(tnew_sgi, tnew_sgi2, told_sgi, told_sgi2,usgi,usgi2,iface,str_ele,un_ele,&
                          u_loc, u_loc2, sn, sn2, tnew_loc, tnew_loc2, told_loc, told_loc2, nloc, ndim)
      implicit none
      ! external vbl
      integer, intent(in) :: nloc, ndim,iface,str_ele,un_ele
      real, dimension(:), intent(inout) :: tnew_sgi, tnew_sgi2, told_loc, told_loc2, told_sgi,told_sgi2
      real, dimension(:,:), intent(inout) :: u_loc2,usgi, usgi2, sn, sn2
      real,intent(in) :: tnew_loc2(:,:,:,:)
      real,pointer, intent(in) :: tnew_loc(:),u_loc(:,:)
      !local vbl
      integer :: iloc, idim

      ! if you want to have this in the code just repleace call get_loc_sgi in the code by the lines below
      usgi=0.0; usgi2=0.0; tnew_sgi=0.0; tnew_sgi2=0.0; told_sgi=0.0; told_sgi2=0.0!; xsgi=0.0

      do iloc=1,nloc ! use all of the nodes not just the surface nodes.
        do idim=1,ndim
          usgi(:,idim)  = usgi(:,idim)  + sn(:,iloc)*u_loc(idim,iloc)
          usgi2(:,idim) = usgi2(:,idim) + sn2(:,iloc)*u_loc2(idim,iloc)
          ! xsgi(:,idim)  = xsgi(:,idim)  + sn(:,iloc)*x_loc(idim,iloc)
        end do
        tnew_sgi  = tnew_sgi(:)  + sn(:,iloc)*tnew_loc(iloc)
        tnew_sgi2 = tnew_sgi2(:) + sn2(:,iloc)*tnew_loc2(iloc,str_ele,un_ele,iface)

        told_sgi  = told_sgi(:)  + sn(:,iloc)*told_loc(iloc)
        told_sgi2 = told_sgi2(:) + sn2(:,iloc)*told_loc2(iloc)
      end do

    end subroutine get_loc_sgi





    ! rectangular ele !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! this subroutine calculates n, nlx for rectangular element
  SUBROUTINE RE2DN4(LOWQUA,NGI,NLOC,MLOC,M,WEIGHT,N,NLX,NLY,SNGI,SNLOC,SWEIGH,SN,SNLX)
    ! use FLDebug
    IMPLICIT NONE
    ! NB might have to define surface elements for p and (u,v,w)
    ! in here as well.
    ! This subroutine defines the shape functions M and N and their
    ! derivatives at the Gauss points
    ! REAL M(1,NGI),WEIGHT(NGI),N(4,NGI),NLX(4,NGI),NLY(4,NGI)
    INTEGER, intent(in):: NGI,NLOC,MLOC
    REAL:: M(NGI,MLOC),WEIGHT(NGI)
    REAL:: N(NGI,NLOC),NLX(NGI,NLOC),NLY(NGI,NLOC)
    REAL:: POSI,TLY
    REAL:: LX(16),LY(16),LXP(4),LYP(4)
    REAL:: WEIT(16)
    INTEGER:: SNGI,SNLOC
    REAL ::SWEIGH(SNGI)
    REAL:: SN(SNGI,SNLOC),SNLX(SNGI,SNLOC)
    INTEGER:: P,Q,CORN,GPOI,ILOC,JLOC,NDGI
    LOGICAL:: LOWQUA,GETNDP
    INTEGER:: I
    ! NB LXP(I) AND LYP(I) ARE THE LOCAL X AND Y COORDS OF NODAL POINT I

    ! ewrite(3,*)'inside re2dn4, nloc,mloc,ngi',&
                    ! nloc,mloc,ngi

    LXP(1)=-1
    LYP(1)=-1

    LXP(2)=1
    LYP(2)=-1

    ! LXP(3)=1
    ! LYP(3)=1

    ! LXP(4)=-1
    ! LYP(4)=1

    LXP(3)=-1
    LYP(3)=1

    LXP(4)=1
    LYP(4)=1

    IF(NGI.EQ.4) THEN
      POSI=1.0/SQRT(3.0)
      LX(1)=-POSI
      LY(1)=-POSI
      LX(2)= POSI
      LY(2)= POSI

      do  Q=1,2! Was loop 23
        do  P=1,2! Was loop 24
          do  CORN=1,4! Was loop 25
            GPOI=(Q-1)*2 + P

            IF(MLOC.EQ.1)  M(GPOI,1)=1.
              WEIGHT(GPOI)=1.

              N(GPOI,CORN)=0.25*(1.+LXP(CORN)*LX(P))&
                          *(1.+LYP(CORN)*LY(Q))
              NLX(GPOI,CORN)=0.25*LXP(CORN)*(1.+LYP(CORN)*LY(Q))
              NLY(GPOI,CORN)=0.25*LYP(CORN)*(1.+LXP(CORN)*LX(P))
          end do ! Was loop 25
        end do ! Was loop 24
      end do ! Was loop 23
      ! ewrite(3,*) 'here 1'
      ! ewrite(3,*) 'N:',N
      ! ewrite(3,*) 'NLX:',NLX
      ! ewrite(3,*) 'NLY:',NLY
      ! Surface shape functions
      IF((SNGI.GT.1).AND.(SNLOC.GT.1)) THEN
         ! ewrite(3,*) '***************** SNGI=',SNGI
        do  P=1,2! Was loop 27
          do  CORN=1,2! Was loop 27
                       GPOI=P
                       SN(GPOI,CORN)=0.5*(1.+LXP(CORN)*LX(P))
                       SNLX(GPOI,CORN)=0.5*LXP(CORN)
                       SWEIGH(GPOI)=1.
            end do ! Was loop 27
        end do ! Was loop 27
      ENDIF
    ! IF(NGI.EQ.4) THEN ...
    ELSE
      NDGI =INT(SQRT(NGI+0.1) +0.1)
      ! ewrite(3,*) 'ndgi,ngi,sngi:',ndgi,ngi,sngi

      GETNDP=.FALSE.
      CALL LAGROT(WEIT,LX,NDGI,GETNDP)
      LY(1:NDGI) = LX(1:NDGI)
      ! ewrite(3,*) 'weit:',weit
      ! ewrite(3,*) 'lx:',lx

      do  Q=1,NDGI! Was loop 323
        do  P=1,NDGI! Was loop 324
          do  CORN=1,4! Was loop 325
            ! ewrite(3,*) 'q,p,corn:',q,p,corn
            GPOI=(Q-1)*NDGI + P
            IF(MLOC.EQ.1)  M(GPOI,1)=1.
            WEIGHT(GPOI)=WEIT(P)*WEIT(Q)
            ! ewrite(3,*) 'here1'
            N(GPOI,CORN)=0.25*(1.+LXP(CORN)*LX(P))&
                             *(1.+LYP(CORN)*LY(Q))
            ! ewrite(3,*) 'here2'
            NLX(GPOI,CORN)=0.25*LXP(CORN)*(1.+LYP(CORN)*LY(Q))
            NLY(GPOI,CORN)=0.25*LYP(CORN)*(1.+LXP(CORN)*LX(P))
            ! ewrite(3,*) 'here3'
          end do ! Was loop 325
        end do ! Was loop 324
      end do ! Was loop 323
      ! ewrite(3,*) 'here 1'
      ! ewrite(3,*) 'N:',N
      ! ewrite(3,*) 'NLX:',NLX
      ! ewrite(3,*) 'NLY:',NLY
      ! Surface shape functions
      ! ewrite(3,*) '***************** SNGI=',SNGI
      IF(SNGI.GT.0) THEN
        GETNDP=.FALSE.
        CALL LAGROT(WEIT,LX,SNGI,GETNDP)
        do  P=1,SNGI! Was loop 327
          do  CORN=1,2! Was loop 327
            GPOI=P
            SN(GPOI,CORN)=0.5*(1.+LXP(CORN)*LX(P))
            SNLX(GPOI,CORN)=0.5*LXP(CORN)
            SWEIGH(GPOI)=WEIT(P)
          end do ! Was loop 327
        end do ! Was loop 327
      ! ENDOF IF(SNGI.GT.0) THEN...
      ENDIF
    ! END OF IF(NGI.EQ.4) THEN ELSE ...
    ENDIF

    IF(MLOC.EQ.NLOC) THEN
      do  I=1,4! Was loop 2545
        do  CORN=1,4! Was loop 2545
          M(I,CORN)=N(I,CORN)
        end do ! Was loop 2545
      end do ! Was loop 2545
    ENDIF
    ! ewrite(3,*) 'in re2dn4.f here 2 ngi,sngi',ngi,sngi
    ! ewrite(3,*) 'N:',N
    ! ewrite(3,*) 'NLX:',NLX
    ! ewrite(3,*) 'NLY:',NLY
    ! END
  END SUBROUTINE RE2DN4


  ! This computes the weight and points for standard Gaussian quadrature.
  SUBROUTINE LAGROT(WEIT,QUAPOS,NDGI,GETNDP)
    ! use RGPTWE_module
    IMPLICIT NONE
    ! This computes the weight and points for standard Gaussian quadrature.
    ! IF(GETNDP) then get the POSITION OF THE NODES
    ! AND DONT BOTHER WITH THE WEITS.
    INTEGER:: NDGI
    REAL:: WEIT(NDGI),QUAPOS(NDGI)
    LOGICAL:: GETNDP
    LOGICAL:: WEIGHT
    INTEGER ::IG
    !real function...
    ! real :: RGPTWE
    !real, allocatable:: RGPTWE(:,:,:)
    !allocate(RGPTWE(ndgi,1,1))

    IF(.NOT.GETNDP) THEN
      WEIGHT=.TRUE.
      do IG=1,NDGI
        WEIT(IG)=RGPTWE(IG,NDGI,WEIGHT)
      END DO

      WEIGHT=.FALSE.
      do IG=1,NDGI
        QUAPOS(IG)=RGPTWE(IG,NDGI,WEIGHT)
      END DO
    ELSE
      IF(NDGI.EQ.1) THEN
        QUAPOS(1)=0.
      ELSE
        do IG=1,NDGI
          QUAPOS(IG)= -1+2.*REAL(IG-1)/REAL(NDGI-1)
        END DO
      ENDIF
    ENDIF
  END SUBROUTINE LAGROT


  ! this subroutine is for rectangular structured meshes only
  subroutine surface_pointers_sn(nface, sngi, snloc, nloc, ndim, sndim, totele, n_s_list_no, nele_x, nele_y, &
               sn_orig, snlx_orig, face_sn, face_snlx, face_sn2, face_list_no)
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
    ! nele_x = no of elements across in the x-direction.
    ! nele_y = no of elements across in the y-direction.
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
    integer, intent(in) :: nface, sngi, snloc, nloc, ndim, sndim, totele, n_s_list_no, nele_x, nele_y
    real, intent(in) :: sn_orig(sngi,snloc), snlx_orig(sngi,sndim,snloc)
    real, intent(out) :: face_sn(sngi,nloc,nface), face_snlx(sngi,sndim,nloc,nface), face_sn2(sngi,nloc,n_s_list_no)
    integer, intent(out) :: face_list_no(nface,totele)
    ! local variables...
    integer iface,lnod1,lnod2,i_s_list_no

    face_sn=0.0
    face_snlx=0.0
    face_sn2=0.0

    ! local node numbers:
    !  3   4
    !  1   2
    !
    ! face numbers:
    !    4
    !  2   3
    !    1

    iface=1
    lnod1=2
    lnod2=1
    face_sn(:,lnod1,iface)=sn_orig(:,1)
    face_sn(:,lnod2,iface)=sn_orig(:,2)
    face_snlx(:,1,lnod1,iface)=snlx_orig(:,1,1)
    face_snlx(:,1,lnod2,iface)=snlx_orig(:,1,2)
    i_s_list_no = iface
    lnod1=4
    lnod2=3
    face_sn2(:,lnod1,i_s_list_no)=sn_orig(:,1)
    face_sn2(:,lnod2,i_s_list_no)=sn_orig(:,2)
    face_list_no(iface,:) = i_s_list_no

    iface=2
    lnod1=1
    lnod2=3
    face_sn(:,lnod1,iface)=sn_orig(:,1)
    face_sn(:,lnod2,iface)=sn_orig(:,2)
    face_snlx(:,1,lnod1,iface)=snlx_orig(:,1,1)
    face_snlx(:,1,lnod2,iface)=snlx_orig(:,1,2)
    i_s_list_no = iface
    lnod1=2
    lnod2=4
    face_sn2(:,lnod1,i_s_list_no)=sn_orig(:,1)
    face_sn2(:,lnod2,i_s_list_no)=sn_orig(:,2)
    face_list_no(iface,:) = i_s_list_no

    iface=3
    lnod1=4
    lnod2=2
    face_sn(:,lnod1,iface)=sn_orig(:,1)
    face_sn(:,lnod2,iface)=sn_orig(:,2)
    face_snlx(:,1,lnod1,iface)=snlx_orig(:,1,1)
    face_snlx(:,1,lnod2,iface)=snlx_orig(:,1,2)
    i_s_list_no = iface
    lnod1=3
    lnod2=1
    face_sn2(:,lnod1,i_s_list_no)=sn_orig(:,1)
    face_sn2(:,lnod2,i_s_list_no)=sn_orig(:,2)
    face_list_no(iface,:) = i_s_list_no

    iface=4
    lnod1=3
    lnod2=4
    face_sn(:,lnod1,iface)=sn_orig(:,1)
    face_sn(:,lnod2,iface)=sn_orig(:,2)
    face_snlx(:,1,lnod1,iface)=snlx_orig(:,1,1)
    face_snlx(:,1,lnod2,iface)=snlx_orig(:,1,2)
    i_s_list_no = iface
    lnod1=1
    lnod2=2
    face_sn2(:,lnod1,i_s_list_no)=sn_orig(:,1)
    face_sn2(:,lnod2,i_s_list_no)=sn_orig(:,2)
    face_list_no(iface,:) = i_s_list_no
    !
    ! calculate face_list_no(nface,totele) =
    !  do jele=1,nele_y
    !    do iele=1,nele_x
    !       ele=(jele-1)*nele_x + iele
    !    end do
    !  end do
  end subroutine surface_pointers_sn

  ! triangles !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! This sub calculates the local corrds L1, L2, L3, L4 and
  !  weights at the quadrature points.
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
            L1(I) = 0.5
            L2(I) = 0.5
            WEIGHT(I)=1.0/3.0
          END DO
          L1(2)=0.
          L2(3)=0.
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

  ! Amin added snloc to all subroutines
  subroutine get_shape_funs_spec(n, nlx, nlx_lxx, nlxx, weight, nlx_nod, &
              nloc, snloc, sngi, ngi, ndim, nface,n_s_list_no, face_sn, face_sn2, face_snlx, &
              sweight, npoly, ele_type, totele)!, face_list_no )
      implicit none
      ! nloc=no of local nodes per element
      ! ngi=no of quadrature points
      ! ndim=no of dimensions.
      integer, intent(in) :: nloc, sngi, ngi, ndim, nface, n_s_list_no, snloc, totele
      ! integer, INTENT(INOUT) :: face_list_no( nface, totele)
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
      real, intent(inout) :: face_sn(sngi,nloc,nface), face_sn2(sngi,nloc,n_s_list_no)
      real, intent(inout) :: face_snlx(sngi,ndim,nloc,nface)
      ! npoly=order of polynomial in Cartesian space; ele_type=type of element including order of poly.
      integer, intent(in) :: npoly,ele_type
      real, allocatable :: sweight(:)

      ! form the shape functions...
      call get_shape_funs_with_faces(n, nlx, weight,  &
                 nloc, snloc, sngi, ngi, ndim, nface,n_s_list_no, face_sn, face_sn2, face_snlx, &
                 sweight, npoly,ele_type, totele)!, face_list_no)

      ! Calculate high order derivatives of shape functions nlx_lxx, nlxx
      ! and also calculate node-wise deriavtives of shape functions nlx_nod.
      ! Amin for higher order uncomment get_high_order_shape_funs below
      ! call get_high_order_shape_funs(n, nlx, nlx_lxx, nlxx, weight, nlx_nod, nloc, ngi, ndim, &
      !              npoly, ele_type)
  end subroutine get_shape_funs_spec


  subroutine get_shape_funs_with_faces(n, nlx, weight,  &
                   nloc, snloc, sngi, ngi, ndim, nface,n_s_list_no, face_sn, face_sn2, face_snlx,  &
                   sweight, npoly,ele_type, totele)!, face_list_no)
      implicit none
      ! nloc=no of local nodes per element
      ! ngi=no of quadrature points
      ! ndim=no of dimensions.
      ! ele_type= element type
      integer, intent(in) :: nloc, sngi, ngi, ndim, nface, n_s_list_no, snloc
      integer, intent(in) :: totele
      ! integer, INTENT(INOUT) :: face_list_no( nface, totele)
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
      real, intent(inout) :: face_snlx(sngi,ndim,nloc,nface)
      ! npoly=order of polynomial in Cartesian space; ele_type=type of element including order of poly.
      integer, intent(in) :: npoly, ele_type
      ! local variables...
      integer is_triangle_or_tet
      parameter(is_triangle_or_tet=100)
      integer sndim, suf_ngi, suf_ndim, suf_nloc, ipoly, IQADRA
      real, allocatable :: rdum1(:), rdum2(:), rdum3(:)
      real, allocatable :: sn_orig(:,:),sn2(:,:),snlx_orig(:,:,:),sweight(:)
      real, allocatable :: suf_n(:,:),suf_nlx(:,:,:),suf_weight(:)
      ! integer, allocatable :: face_list_no(:,:)

      ! allocate memory...
      allocate(rdum1(10000), rdum2(10000), rdum3(10000) )
      allocate(sn_orig(sngi,snloc),sn2(sngi,nloc),snlx_orig(sngi,ndim,snloc),sweight(sngi) )
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
        call tri_surface_pointers_sn(nface, sngi, snloc, nloc, ndim, sndim, totele, n_s_list_no, &
                       sn_orig, snlx_orig, face_sn, face_snlx, face_sn2)!, face_list_no)
        ! return a surface tet...
        ! Amin uncomment the 6 lines below
        ! suf_ngi=NGI/IQADRA
        ! suf_ndim=ndim-1
        ! suf_nloc=NLOC/(IPOLY+1)
        ! call get_shape_funs(suf_ngi,suf_nloc,suf_ndim,  &
        !           suf_weight,suf_n,suf_nlx, ipoly,iqadra, &
        !           sngi, sndim, rdum1,rdum2,rdum3, .false.   )

      endif
  end subroutine get_shape_funs_with_faces

  ! this subroutine is for triangular structured meshes only generated by tri_ele_info in the structured_meshgen.F90
  subroutine tri_surface_pointers_sn(nface, sngi, snloc, nloc, ndim, sndim, totele, n_s_list_no, &
                   sn_orig, snlx_orig, face_sn, face_snlx, face_sn2)!, face_list_no)
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
      ! no_ele_row = no of elements across in the x-direction.
      ! no_ele_col = no of elements across in the y-direction.
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
      real, intent(out) :: face_sn(sngi,nloc,nface), face_snlx(sngi,sndim,nloc,nface), face_sn2(sngi,nloc,n_s_list_no)
      ! integer, intent(out) :: face_list_no(nface,totele)
      ! local variables...
      integer iface,lnod1,lnod2,i_s_list_no

      face_sn(:,:,:)=0.0
      face_snlx=0.0!(:,:,:,:)
      face_sn2(:,:,:)=0.0

      ! local face numbers:
      !    |\        __1__
      !   2| \ 3     \   |
      !    |  \      3\  | 2
      !    |___\       \ |
      !      1          \|

      ! local node numbers
      !    2
      !    |\        1___3
      !    | \       \   |
      !    |  \       \  |
      !    |___\       \ |
      !    3   1        \|
      !                 2

      iface=1
      lnod1=1
      lnod2=3
      face_sn(:,lnod1,iface)=sn_orig(:,1)
      face_sn(:,lnod2,iface)=sn_orig(:,2)
      face_snlx(:,1,lnod1,iface)=snlx_orig(:,1,1)
      face_snlx(:,1,lnod2,iface)=snlx_orig(:,1,2)
      ! i_s_list_no = iface
      lnod1=3
      lnod2=1
      ! face_sn2(:,lnod1,i_s_list_no)=sn_orig(:,1)
      ! face_sn2(:,lnod2,i_s_list_no)=sn_orig(:,2)
      ! face_list_no(iface,:) = i_s_list_no
      face_sn2(:,lnod1,iface)=sn_orig(:,1)
      face_sn2(:,lnod2,iface)=sn_orig(:,2)


      iface=2
      lnod1=3
      lnod2=2
      face_sn(:,lnod1,iface)=sn_orig(:,1)
      face_sn(:,lnod2,iface)=sn_orig(:,2)
      face_snlx(:,1,lnod1,iface)=snlx_orig(:,1,1)
      face_snlx(:,1,lnod2,iface)=snlx_orig(:,1,2)
      ! i_s_list_no = iface
      lnod1=2
      lnod2=3
      ! face_sn2(:,lnod1,i_s_list_no)=sn_orig(:,1)
      ! face_sn2(:,lnod2,i_s_list_no)=sn_orig(:,2)
      ! face_list_no(iface,:) = i_s_list_no
      face_sn2(:,lnod1,iface)=sn_orig(:,1)
      face_sn2(:,lnod2,iface)=sn_orig(:,2)

      iface=3
      lnod1=2
      lnod2=1
      face_sn(:,lnod1,iface)=sn_orig(:,1)
      face_sn(:,lnod2,iface)=sn_orig(:,2)
      face_snlx(:,1,lnod1,iface)=snlx_orig(:,1,1)
      face_snlx(:,1,lnod2,iface)=snlx_orig(:,1,2)
      ! i_s_list_no = iface
      lnod1=1
      lnod2=2
      ! face_sn2(:,lnod1,i_s_list_no)=sn_orig(:,1)
      ! face_sn2(:,lnod2,i_s_list_no)=sn_orig(:,2)
      ! face_list_no(iface,:) = i_s_list_no
      face_sn2(:,lnod1,iface)=sn_orig(:,1)
      face_sn2(:,lnod2,iface)=sn_orig(:,2)

      ! calculate face_list_no(nface,totele) =
      !  do jele=1,no_ele_col
      !    do iele=1,no_ele_row
      !       ele=(jele-1)*no_ele_row + iele
      !    end do
      !  end do
  end subroutine tri_surface_pointers_sn


  ! form volume and surface shape functions
  subroutine get_shape_funs(ngi,nloc,ndim,  &
                 weight,n,nlx, ipoly,iqadra, &
                 snloc, sngi, sndim, sweight,sn_orig,snlx_orig, with_time_slab   )
      ! ***************************************
      ! form volume and surface shape functions
      ! ***************************************
      IMPLICIT NONE
      INTEGER, intent(in) :: sngi, NGI, NLOC, ndim, sndim, snloc
      INTEGER, intent(in) :: IPOLY,IQADRA
      logical, intent(in) :: with_time_slab
      REAL, intent(inout) ::  n(ngi,nloc)
      REAL, intent(inout) ::  nlx(ngi,ndim,nloc)
      REAL, intent(inout) ::  sn_orig(sngi, snloc)
      REAL, intent(inout) ::  snlx_orig(sngi, sndim, snloc)
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
        sweight=0.0; sN_orig=0.0; sNLX_orig=0.0

        !AMIN turn on interface_SPECTR
        ! call interface_SPECTR(sNGI,sNLOC_temp,sweight,sN,sNLX, ndim-1, IPOLY,IQADRA  )

      ELSE
        ! volume tet plus time slab...
        call triangles_tets_with_time( nloc,ngi, ndim, n,nlx, weight, ipoly,iqadra, with_time_slab)
        ! surfaces...
        ! AMIN I've changed nloc to snloc
        call triangles_tets_with_time( snloc,sngi, ndim-1, sn_orig,snlx_orig, sweight, ipoly,iqadra, with_time_slab)
      ENDIF
  end subroutine get_shape_funs


  subroutine triangles_tets_with_time( nloc,ngi, ndim, n,nlx, weight, ipoly,iqadra, with_time_slab)
      ! ****************************************************************************************
      implicit none
      integer, intent(in) :: nloc,ngi,ndim, ipoly,iqadra
      real, intent(inout) :: n(:,:), nlx(:,:,:)
      real, intent(inout) :: weight(:)
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
        ! extend into time domain...
        !AMIN turn on interface_SPECTR
        ! call interface_SPECTR(NGI_t,NLOC_t,WEIGHT_t,N_t,NLX_t, ndim_t, IPOLY,IQADRA  )

        ! combine-space time...
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
  SUBROUTINE SHATRInew(L1, L2, L3, L4, WEIGHT,NLOC,NGI,ndim,N,NLX)
      ! Interface to SHATRIold using the new style variables
      IMPLICIT NONE
      INTEGER , intent(in) :: NLOC,NGI,ndim
      REAL , dimension(ngi), intent(in) :: L1, L2, L3, L4
      REAL , dimension(ngi), intent(inout) :: WEIGHT
      REAL , dimension(:, : ), intent(inout) ::N
      real, dimension (:,:,:), intent(inout) :: NLX

      call SHATRIold(L1, L2, L3, L4, WEIGHT, .false., &
             NLOC,NGI,ndim,  &
             N,NLX)
  end subroutine SHATRInew



  SUBROUTINE SHATRIold(L1, L2, L3, L4, WEIGHT, D3,NLOC,NGI,ndim,N,NLX)
      ! Work out the shape functions and there derivatives...
      IMPLICIT NONE
      INTEGER , intent(in) :: NLOC,NGI,ndim
      LOGICAL , intent(in) :: D3
      REAL , dimension(ngi), intent(in) :: L1, L2, L3, L4
      REAL , dimension(ngi), intent(inout) :: WEIGHT
      REAL , dimension(:, : ), intent(inout) ::N
      real, dimension(:,:,:), INTENT(INOUT) :: NLX ! (ngi,ndim,nloc)
      ! Local variables...
      INTEGER ::  GI

      INTEGER :: P,CORN,GPOI
      REAL :: LXP(2), LX(2)

      IF(.NOT.D3) THEN
        ! Assume a triangle...

        IF(NLOC.EQ.1) THEN
          Loop_Gi_Nloc1: DO GI=1,NGI
            N(GI,1)=1.0
            NLX(GI,1,1)=0.0
            NLX(GI,2,1)=0.0
          end DO Loop_Gi_Nloc1
        ELSE IF((NLOC.EQ.3).OR.(NLOC.EQ.4)) THEN
          Loop_Gi_Nloc3_4: DO GI=1,NGI
            N(GI,1)=L1(GI)
            N(GI,2)=L2(GI)
            N(GI,3)=L3(GI)

            NLX(GI,1,1)=1.0
            NLX(GI,1,2)=0.0
            NLX(GI,1,3)=-1.0

            NLX(GI,2,1)=0.0
            NLX(GI,2,2)=1.0
            NLX(GI,2,3)=-1.0
            ! IF(NLOC.EQ.4) THEN
            !   ! Bubble function...
            !   !alpha == 1 behaves better than the correct value of 27. See Osman et al. 2019
            !   N(4,GI)  =1. * L1(GI)*L2(GI)*L3(GI)
            !   NLX(4,GI)=1. * L2(GI)*(1.-L2(GI))-2.*L1(GI)*L2(GI)
            !   NLY(4,GI)=1. * L1(GI)*(1.-L1(GI))-2.*L1(GI)*L2(GI)
            ! ENDIF
          end DO Loop_Gi_Nloc3_4
        ELSE IF((NLOC.EQ.6).OR.(NLOC.EQ.7)) THEN
          ! Loop_Gi_Nloc_6_7: DO GI=1,NGI
          !   N(1,GI)=(2.*L1(GI)-1.)*L1(GI)
          !   N(2,GI)=(2.*L2(GI)-1.)*L2(GI)
          !   N(3,GI)=(2.*L3(GI)-1.)*L3(GI)
          !
          !   N(4,GI)=4.*L1(GI)*L2(GI)
          !   N(5,GI)=4.*L2(GI)*L3(GI)
          !   N(6,GI)=4.*L1(GI)*L3(GI)
          !
          !   ! nb L1+L2+L3+L4=1
          !   ! x-derivative...
          !   NLX(1,GI)=4.*L1(GI)-1.
          !   NLX(2,GI)=0.
          !   NLX(3,GI)=-4.*(1.-L2(GI))+4.*L1(GI) + 1.
          !
          !   NLX(4,GI)=4.*L2(GI)
          !   NLX(5,GI)=-4.*L2(GI)
          !   NLX(6,GI)=4.*(1.-L2(GI))-8.*L1(GI)
          !
          !   ! y-derivative...
          !   NLY(1,GI)=0.
          !   NLY(2,GI)=4.*L2(GI)-1.0
          !   NLY(3,GI)=-4.*(1.-L1(GI))+4.*L2(GI) + 1.
          !
          !   NLY(4,GI)=4.*L1(GI)
          !   NLY(5,GI)=4.*(1.-L1(GI))-8.*L2(GI)
          !   NLY(6,GI)=-4.*L1(GI)
          !   IF(NLOC.EQ.7) THEN
          !     ! Bubble function...
          !     N(7,GI)  =L1(GI)*L2(GI)*L3(GI)
          !     NLX(7,GI)=L2(GI)*(1.-L2(GI))-2.*L1(GI)*L2(GI)
          !     NLY(7,GI)=L1(GI)*(1.-L1(GI))-2.*L1(GI)*L2(GI)
          !   ENDIF
          ! END DO Loop_Gi_Nloc_6_7
          ! ENDOF IF(NLOC.EQ.6) THEN...
        ELSE IF(NLOC==10) THEN ! Cubic triangle...
          ! get the shape functions for a cubic triangle...
          ! uncomment line below to activate shape_triangle_cubic
          ! call shape_triangle_cubic( l1, l2, l3, l4, weight, d3, &
          !          nloc, ngi, &
          !          n, nlx, nly, nlz )
          ! AMIN I've added this else if statement.
        else if ( ndim.EQ.1 ) then !surface integral in 1D
            LXP(1)=-1
            LXP(2)= 1
            LX(1)=-1.0/SQRT(3.0)
            Lx(2)=1.0/SQRT(3.0)
            do  P=1,NGI
               do  CORN=1,2! Was loop 327
                  GPOI=P
                  N(GPOI,CORN)=0.5*(1.+LXP(CORN)*LX(P))
                  NLX(GPOI,ndim,CORN)=0.5*LXP(CORN)
                  WEIGHT(GPOI)=1
               end do
            end do

        ELSE ! has not found the element shape functions
          stop 811
        ENDIF
        ! ENDOF IF(.NOT.D3) THEN
      ENDIF


      IF(D3) THEN
        ! ! Assume a tet...
        ! ! This is for 5 point quadrature.
        ! IF((NLOC.EQ.10).OR.(NLOC.EQ.11)) THEN
        !   Loop_Gi_Nloc_10_11: DO GI=1,NGI
        !     !ewrite(3,*)'gi,L1(GI),L2(GI),L3(GI),L4(GI):',gi,L1(GI),L2(GI),L3(GI),L4(GI)
        !     N(1,GI)=(2.*L1(GI)-1.)*L1(GI)
        !     N(3,GI)=(2.*L2(GI)-1.)*L2(GI)
        !     N(5,GI)=(2.*L3(GI)-1.)*L3(GI)
        !     N(10,GI)=(2.*L4(GI)-1.)*L4(GI)
        !
        !     !if(L1(GI).gt.-1.93) ewrite(3,*)'gi,L1(GI), L2(GI), L3(GI), L4(GI),N(1,GI):', &
        !     !                            gi,L1(GI), L2(GI), L3(GI), L4(GI),N(1,GI)
        !     !
        !     !
        !     N(2,GI)=4.*L1(GI)*L2(GI)
        !     N(6,GI)=4.*L1(GI)*L3(GI)
        !     N(7,GI)=4.*L1(GI)*L4(GI)
        !
        !     N(4,GI) =4.*L2(GI)*L3(GI)
        !     N(9,GI) =4.*L3(GI)*L4(GI)
        !     N(8,GI)=4.*L2(GI)*L4(GI)
        !     ! nb L1+L2+L3+L4=1
        !     ! x-derivative...
        !     NLX(1,GI)=4.*L1(GI)-1.
        !     NLX(3,GI)=0.
        !     NLX(5,GI)=0.
        !     NLX(10,GI)=-4.*(1.-L2(GI)-L3(GI))+4.*L1(GI) + 1.
        !     !if(L1(GI).gt.-1.93) ewrite(3,*)'Nlx(1,GI):', &
        !     !     Nlx(1,GI)
        !
        !     NLX(2,GI)=4.*L2(GI)
        !     NLX(6,GI)=4.*L3(GI)
        !     NLX(7,GI)=4.*(L4(GI)-L1(GI))
        !
        !     NLX(4,GI) =0.
        !     NLX(9,GI) =-4.*L3(GI)
        !     NLX(8,GI)=-4.*L2(GI)
        !
        !     ! y-derivative...
        !     NLY(1,GI)=0.
        !     NLY(3,GI)=4.*L2(GI)-1.0
        !     NLY(5,GI)=0.
        !     NLY(10,GI)=-4.*(1.-L1(GI)-L3(GI))+4.*L2(GI) + 1.
        !
        !     NLY(2,GI)=4.*L1(GI)
        !     NLY(6,GI)=0.
        !     NLY(7,GI)=-4.*L1(GI)
        !
        !     NLY(4,GI) =4.*L3(GI)
        !     NLY(9,GI) =-4.*L3(GI)
        !     NLY(8,GI)=4.*(1-L1(GI)-L3(GI))-8.*L2(GI)
        !
        !     ! z-derivative...
        !     NLZ(1,GI)=0.
        !     NLZ(3,GI)=0.
        !     NLZ(5,GI)=4.*L3(GI)-1.
        !     NLZ(10,GI)=-4.*(1.-L1(GI)-L2(GI))+4.*L3(GI) + 1.
        !
        !     NLZ(2,GI)=0.
        !     NLZ(6,GI)=4.*L1(GI)
        !     NLZ(7,GI)=-4.*L1(GI)
        !
        !     NLZ(4,GI) =4.*L2(GI)
        !     NLZ(9,GI) =4.*(1.-L1(GI)-L2(GI))-8.*L3(GI)
        !     NLZ(8,GI)=-4.*L2(GI)
        !     IF(NLOC.EQ.11) THEN
        !       ! Bubble function...
        !       N(11,GI)  =L1(GI)*L2(GI)*L3(GI)*L4(GI)
        !       NLX(11,GI)=L2(GI)*L3(GI)*(1.-L2(GI)-L3(GI))-2.*L1(GI)*L2(GI)*L3(GI)
        !       NLY(11,GI)=L1(GI)*L3(GI)*(1.-L1(GI)-L3(GI))-2.*L1(GI)*L2(GI)*L3(GI)
        !       NLZ(11,GI)=L1(GI)*L2(GI)*(1.-L1(GI)-L2(GI))-2.*L1(GI)*L2(GI)*L3(GI)
        !     ENDIF
        !   end DO Loop_Gi_Nloc_10_11
        !   ! ENDOF IF(NLOC.EQ.10) THEN...
        ! ENDIF
        !
        ! IF((NLOC.EQ.4).OR.(NLOC.EQ.5)) THEN
        !   Loop_Gi_Nloc_4_5: DO GI=1,NGI
        !     N(1,GI)=L1(GI)
        !     N(2,GI)=L2(GI)
        !     N(3,GI)=L3(GI)
        !     N(4,GI)=L4(GI)
        !
        !     NLX(1,GI)=1.0
        !     NLX(2,GI)=0
        !     NLX(3,GI)=0
        !     NLX(4,GI)=-1.0
        !
        !     NLY(1,GI)=0.0
        !     NLY(2,GI)=1.0
        !     NLY(3,GI)=0.0
        !     NLY(4,GI)=-1.0
        !
        !     NLZ(1,GI)=0.0
        !     NLZ(2,GI)=0.0
        !     NLZ(3,GI)=1.0
        !     NLZ(4,GI)=-1.0
        !     IF(NLOC.EQ.5) THEN
        !       ! Bubble function ...
        !       !alpha == 50 behaves better than the correct value of 256. See Osman et al. 2019
        !       N(5,GI)  = 50. * L1(GI)*L2(GI)*L3(GI)*L4(GI)
        !       NLX(5,GI)= 50. * L2(GI)*L3(GI)*(1.-L2(GI)-L3(GI))-2.*L1(GI)*L2(GI)*L3(GI)
        !       NLY(5,GI)= 50. * L1(GI)*L3(GI)*(1.-L1(GI)-L3(GI))-2.*L1(GI)*L2(GI)*L3(GI)
        !       NLZ(5,GI)= 50. * L1(GI)*L2(GI)*(1.-L1(GI)-L2(GI))-2.*L1(GI)*L2(GI)*L3(GI)
        !     ENDIF
        !   end DO Loop_Gi_Nloc_4_5
        ! ENDIF
        ! IF(NLOC.EQ.1) THEN
        !   Loop_Gi_Nloc_1: DO GI=1,NGI
        !     N(1,GI)=1.0
        !     NLX(1,GI)=0.0
        !     NLY(1,GI)=0.0
        !     NLZ(1,GI)=0.0
        !   end DO Loop_Gi_Nloc_1
        ! ENDIF
        ! ! ENDOF IF(D3) THEN...
        print*, 'deal with it later'
      ENDIF
      RETURN
  END SUBROUTINE SHATRIold


  !subroutine calculates det of J time weights @ quadrature points for rectangles
  subroutine det_nlx( x_loc, n, nlx, nx, detwei, weight, ndim, nloc, ngi, INV_JAC )
    ! ****************************************************
    ! This sub form the derivatives of the shape functions
    ! ****************************************************
    ! x_loc: spatial nodes.
    ! n, nlx, nlx_lxx: shape function and local derivatives of the shape functions (nlx_lxx is local grad of the local laplacian)- defined by shape functional library.
    ! nx, nx_lxx: derivatives of the shape functions.
    ! detwei, inv_jac: determinant at the quadrature pots and inverse of Jacobian at quadrature pts.
    ! ndim,nloc,ngi: no of dimensions, no of local nodes within an element, no of quadrature points.
    ! nlx_nod, nx_nod: same as nlx and nx but formed at the nodes not quadrature points.
    implicit none
    integer, intent( in ) :: ndim,nloc,ngi

    REAL, DIMENSION( ndim,nloc ), intent( in ) :: x_loc
    REAL, DIMENSION( ngi, nloc ), intent( in ) :: N
    REAL, DIMENSION( ngi, ndim, nloc ), intent( in ) :: nlx
    REAL, DIMENSION( ngi ), intent( in ) :: WEIGHT
    REAL, DIMENSION( ngi ), intent( inout ) :: DETWEI
    REAL, DIMENSION( ngi, ndim, nloc ), intent( inout ) :: nx
    REAL, DIMENSION( ngi, ndim, ndim ), intent( inout ):: INV_JAC
    ! Local variables
    REAL :: AGI, BGI, CGI, DGI, EGI, FGI, GGI, HGI, KGI, A11, A12, A13, A21, &
            A22, A23, A31, A32, A33, DETJ
    INTEGER :: GI, L, IGLX, ii

    if (ndim==2) then
     ! conventional:
      do  GI=1,NGI! Was loop 331

        AGI=0.
        BGI=0.

        CGI=0.
        DGI=0.

        do  L=1,NLOC! Was loop 79
          IGLX=L

          AGI=AGI+NLX(GI,1,L)*x_loc(1,L)
          BGI=BGI+NLX(GI,1,L)*x_loc(2,L)

          CGI=CGI+NLX(GI,2,L)*x_loc(1,L)
          DGI=DGI+NLX(GI,2,L)*x_loc(2,L)

        end do ! Was loop 79

        DETJ= AGI*DGI-BGI*CGI
        DETWEI(GI)=ABS(DETJ)*WEIGHT(GI)
        ! For coefficient in the inverse mat of the jacobian.
        A11= DGI /DETJ
        A21=-BGI /DETJ

        A12=-CGI /DETJ
        A22= AGI /DETJ

        INV_JAC( GI, 1,1 )= A11
        INV_JAC( GI, 1,2 )= A21

        INV_JAC( GI, 2,1 )= A12
        INV_JAC( GI, 2,2 )= A22

        do  L=1,NLOC! Was loop 373
          NX(GI,1,L)= A11*NLX(GI,1,L)+A12*NLX(GI,2,L)
          NX(GI,2,L)= A21*NLX(GI,1,L)+A22*NLX(GI,2,L)
        end do ! Was loop 373
      end do ! GI Was loop 331

      !jac(1) = AGI; jac(2) = DGI ; jac(3) = BGI ; jac(4) = EGI

    elseif ( ndim.eq.3 ) then
      do  GI=1,NGI! Was loop 331

        AGI=0.
        BGI=0.
        CGI=0.

        DGI=0.
        EGI=0.
        FGI=0.

        GGI=0.
        HGI=0.
        KGI=0.

        do  L=1,NLOC! Was loop 79
          IGLX=L
          !ewrite(3,*)'xndgln, x, nl:', &
          !     iglx, l, x(iglx), y(iglx), z(iglx), NLX(L,GI), NLY(L,GI), NLZ(L,GI)
          ! NB R0 does not appear here although the z-coord might be Z+R0.
          AGI=AGI+NLX(GI,1,L)*x_loc(1,IGLX)
          BGI=BGI+NLX(GI,1,L)*x_loc(2,IGLX)
          CGI=CGI+NLX(GI,1,L)*x_loc(3,IGLX)

          DGI=DGI+NLX(GI,2,L)*x_loc(1,IGLX)
          EGI=EGI+NLX(GI,2,L)*x_loc(2,IGLX)
          FGI=FGI+NLX(GI,2,L)*x_loc(3,IGLX)

          GGI=GGI+NLX(GI,3,L)*x_loc(1,IGLX)
          HGI=HGI+NLX(GI,3,L)*x_loc(2,IGLX)
          KGI=KGI+NLX(GI,3,L)*x_loc(3,IGLX)
        end do ! Was loop 79

        DETJ=AGI*(EGI*KGI-FGI*HGI)&
            -BGI*(DGI*KGI-FGI*GGI)&
            +CGI*(DGI*HGI-EGI*GGI)
        DETWEI(GI)=ABS(DETJ)*WEIGHT(GI)
        ! ewrite(3,*)'gi, detj, weight(gi)', gi, detj, weight(gi)
        ! rsum = rsum + detj
        ! rsumabs = rsumabs + abs( detj )
        ! For coefficient in the inverse mat of the jacobian.
        A11= (EGI*KGI-FGI*HGI) /DETJ
        A21=-(DGI*KGI-FGI*GGI) /DETJ
        A31= (DGI*HGI-EGI*GGI) /DETJ

        A12=-(BGI*KGI-CGI*HGI) /DETJ
        A22= (AGI*KGI-CGI*GGI) /DETJ
        A32=-(AGI*HGI-BGI*GGI) /DETJ

        A13= (BGI*FGI-CGI*EGI) /DETJ
        A23=-(AGI*FGI-CGI*DGI) /DETJ
        A33= (AGI*EGI-BGI*DGI) /DETJ

        INV_JAC( GI, 1,1 )= A11
        INV_JAC( GI, 2,1 )= A21
        INV_JAC( GI, 3,1 )= A31
            !
        INV_JAC( GI, 1,2 )= A12
        INV_JAC( GI, 2,2 )= A22
        INV_JAC( GI, 3,2 )= A32
            !
        INV_JAC( GI, 1,3 )= A13
        INV_JAC( GI, 2,3 )= A23
        INV_JAC( GI, 3,3 )= A33

        do  L=1,NLOC! Was loop 373
          NX(GI,1,L)= A11*NLX(GI,1,L)+A12*NLX(GI,2,L)+A13*NLX(GI,3,L)
          NX(GI,2,L)= A21*NLX(GI,1,L)+A22*NLX(GI,2,L)+A23*NLX(GI,3,L)
          NX(GI,3,L)= A31*NLX(GI,1,L)+A32*NLX(GI,2,L)+A33*NLX(GI,3,L)
        end do ! Was loop 373
      end do ! GI Was loop 331
    end if
  end subroutine det_nlx

  !subroutine calculates det of J time weights @ quadrature points for triangles
  subroutine tri_det_nlx( x_loc, n, nlx, nx, detwei, weight, ndim, nloc, ngi, INV_JAC )
    ! ****************************************************
    ! This sub form the derivatives of the shape functions
    ! ****************************************************
    ! x_loc: spatial nodes.
    ! n, nlx, nlx_lxx: shape function and local derivatives of the shape functions (nlx_lxx is local grad of the local laplacian)- defined by shape functional library.
    ! nx, nx_lxx: derivatives of the shape functions.
    ! detwei, inv_jac: determinant at the quadrature pots and inverse of Jacobian at quadrature pts.
    ! ndim,nloc,ngi: no of dimensions, no of local nodes within an element, no of quadrature points.
    ! nlx_nod, nx_nod: same as nlx and nx but formed at the nodes not quadrature points.
    implicit none
    integer, intent( in ) :: ndim,nloc,ngi

    REAL, DIMENSION( ndim,nloc ), intent( in ) :: x_loc
    REAL, DIMENSION( ngi, nloc ), intent( in ) :: N
    REAL, DIMENSION( ngi, ndim, nloc ), intent( in ) :: nlx
    REAL, DIMENSION( ngi ), intent( in ) :: WEIGHT
    REAL, DIMENSION( ngi ), intent( inout ) :: DETWEI
    REAL, DIMENSION( ngi, ndim, nloc ), intent( inout ) :: nx
    REAL, DIMENSION( ngi, ndim, ndim ), intent( inout ):: INV_JAC
    ! Local variables
    REAL :: AGI, BGI, CGI, DGI, EGI, FGI, GGI, HGI, KGI, A11, A12, A13, A21, &
            A22, A23, A31, A32, A33, DETJ
    INTEGER :: GI, L, IGLX, ii

    if (ndim==2) then
       ! conventional:
        do  GI=1,NGI! Was loop 331

          AGI=0.
          BGI=0.

          CGI=0.
          DGI=0.

          do  L=1,NLOC! Was loop 79
            IGLX=L

            AGI=AGI+NLX(GI,1,L)*x_loc(1,L)
            BGI=BGI+NLX(GI,1,L)*x_loc(2,L)

            CGI=CGI+NLX(GI,2,L)*x_loc(1,L)
            DGI=DGI+NLX(GI,2,L)*x_loc(2,L)

          end do ! Was loop 79

          DETJ= AGI*DGI-BGI*CGI
          DETWEI(GI)=0.5*ABS(DETJ)*WEIGHT(GI)
          ! For coefficient in the inverse mat of the jacobian.
          A11= DGI /DETJ
          A21=-CGI /DETJ

          A12=-BGI /DETJ
          A22= AGI /DETJ

          INV_JAC( GI, 1,1 )= A11
          INV_JAC( GI, 1,2 )= A21

          INV_JAC( GI, 2,1 )= A12
          INV_JAC( GI, 2,2 )= A22

          do  L=1,NLOC! Was loop 373
            NX(GI,1,L)= A11*NLX(GI,1,L)+A12*NLX(GI,2,L)
            NX(GI,2,L)= A21*NLX(GI,1,L)+A22*NLX(GI,2,L)
          end do ! Was loop 373
        end do ! GI Was loop 331

        !jac(1) = AGI; jac(2) = DGI ; jac(3) = BGI ; jac(4) = EGI

    elseif ( ndim.eq.3 ) then
      do  GI=1,NGI! Was loop 331

        AGI=0.
        BGI=0.
        CGI=0.

        DGI=0.
        EGI=0.
        FGI=0.

        GGI=0.
        HGI=0.
        KGI=0.

        do  L=1,NLOC! Was loop 79
          IGLX=L
          !ewrite(3,*)'xndgln, x, nl:', &
          !     iglx, l, x(iglx), y(iglx), z(iglx), NLX(L,GI), NLY(L,GI), NLZ(L,GI)
          ! NB R0 does not appear here although the z-coord might be Z+R0.
          AGI=AGI+NLX(GI,1,L)*x_loc(1,IGLX)
          BGI=BGI+NLX(GI,1,L)*x_loc(2,IGLX)
          CGI=CGI+NLX(GI,1,L)*x_loc(3,IGLX)

          DGI=DGI+NLX(GI,2,L)*x_loc(1,IGLX)
          EGI=EGI+NLX(GI,2,L)*x_loc(2,IGLX)
          FGI=FGI+NLX(GI,2,L)*x_loc(3,IGLX)

          GGI=GGI+NLX(GI,3,L)*x_loc(1,IGLX)
          HGI=HGI+NLX(GI,3,L)*x_loc(2,IGLX)
          KGI=KGI+NLX(GI,3,L)*x_loc(3,IGLX)
        end do ! Was loop 79

        DETJ=AGI*(EGI*KGI-FGI*HGI)&
            -BGI*(DGI*KGI-FGI*GGI)&
            +CGI*(DGI*HGI-EGI*GGI)
        DETWEI(GI)=ABS(DETJ)*WEIGHT(GI)
        ! ewrite(3,*)'gi, detj, weight(gi)', gi, detj, weight(gi)
        ! rsum = rsum + detj
        ! rsumabs = rsumabs + abs( detj )
        ! For coefficient in the inverse mat of the jacobian.
        A11= (EGI*KGI-FGI*HGI) /DETJ
        A21=-(DGI*KGI-FGI*GGI) /DETJ
        A31= (DGI*HGI-EGI*GGI) /DETJ

        A12=-(BGI*KGI-CGI*HGI) /DETJ
        A22= (AGI*KGI-CGI*GGI) /DETJ
        A32=-(AGI*HGI-BGI*GGI) /DETJ

        A13= (BGI*FGI-CGI*EGI) /DETJ
        A23=-(AGI*FGI-CGI*DGI) /DETJ
        A33= (AGI*EGI-BGI*DGI) /DETJ

        INV_JAC( GI, 1,1 )= A11
        INV_JAC( GI, 2,1 )= A21
        INV_JAC( GI, 3,1 )= A31
            !
        INV_JAC( GI, 1,2 )= A12
        INV_JAC( GI, 2,2 )= A22
        INV_JAC( GI, 3,2 )= A32
            !
        INV_JAC( GI, 1,3 )= A13
        INV_JAC( GI, 2,3 )= A23
        INV_JAC( GI, 3,3 )= A33

        do  L=1,NLOC! Was loop 373
          NX(GI,1,L)= A11*NLX(GI,1,L)+A12*NLX(GI,2,L)+A13*NLX(GI,3,L)
          NX(GI,2,L)= A21*NLX(GI,1,L)+A22*NLX(GI,2,L)+A23*NLX(GI,3,L)
          NX(GI,3,L)= A31*NLX(GI,1,L)+A32*NLX(GI,2,L)+A33*NLX(GI,3,L)
        end do ! Was loop 373
      end do ! GI Was loop 331
  end if
  end subroutine tri_det_nlx



  SUBROUTINE det_snlx_all( SNLOC, SNGI, SNDIM, ndim, XSL_ALL, SN, SNLX, SWEIGHT, SDETWE, SAREA, NORMXN_ALL, NORMX_ALL )
    !inv_jac )
    IMPLICIT NONE
    INTEGER, intent( in ) :: SNLOC, SNGI, SNDIM, ndim
    REAL, DIMENSION( NDIM, SNLOC ), intent( in ) :: XSL_ALL
    REAL, DIMENSION( SNGI, SNLOC ), intent( in ) :: SN
    REAL, DIMENSION( SNGI, SNDIM, SNLOC ), intent( in ) :: SNLX
    REAL, DIMENSION( SNGI ), intent( in ) :: SWEIGHT
    REAL, DIMENSION( SNGI ), intent( inout ) :: SDETWE
    REAL, intent( inout ) ::  SAREA
    REAL, DIMENSION( sngi, ndim ), intent( inout ) :: NORMXN_ALL
    REAL, DIMENSION( ndim ), intent( in ) :: NORMX_ALL
    !REAL, DIMENSION( NDIM,ndim ), intent( in ) :: inv_jac
    ! Local variables
    INTEGER :: GI, SL, IGLX
    REAL :: DXDLX, DXDLY, DYDLX, DYDLY, DZDLX, DZDLY
    REAL :: A, B, C, DETJ, RUB3=0.0!, RUB4

    SAREA=0.

    if (ndim==2) then
      DO GI=1,SNGI

        DXDLX=0.
        DXDLY=0.
        DYDLX=0.
        DYDLY=0.
        DZDLX=0.
        DZDLY=0.

        DO SL=1,SNLOC
          DXDLX=DXDLX + SNLX(GI,1,SL)*XSL_ALL(1,SL)
          ! DXDLY=DXDLY + SNLX(GI,2,SL)*XSL_ALL(1,SL)
          DYDLX=DYDLX + SNLX(GI,1,SL)*XSL_ALL(2,SL)
          ! DYDLY=DYDLY + SNLX(GI,2,SL)*XSL_ALL(2,SL)
          ! DZDLX=DZDLX + SNLX(GI,1,SL)*XSL_ALL(3,SL)
          ! DZDLY=DZDLY + SNLX(GI,2,SL)*XSL_ALL(3,SL)
        END DO

        A = DYDLX!*DZDLY - DYDLY*DZDLX
        B = DXDLX!*DZDLY - DXDLY*DZDLX
        ! C = DXDLX*DYDLY - DXDLY*DYDLX

        DETJ=SQRT( A**2 + B**2)! + C**2)
        ! inv_jac(1,1)=DXDLX; inv_jac(1,2)=DXDLY; inv_jac(1,3)=DXDLZ
        ! inv_jac(2,1)=DyDLX; inv_jac(2,2)=DyDLY; inv_jac(2,3)=DyDLZ
        ! inv_jac(3,1)=DzDLX; inv_jac(3,2)=DzDLY; inv_jac(3,3)=DzDLZ
        ! inv_jac=inv_jac/detj
        SDETWE(GI)=DETJ*SWEIGHT(GI)
        SAREA=SAREA+SDETWE(GI)

        ! Calculate the normal at the Gauss pts...
        ! Perform x-product. N=T1 x T2
        CALL NORMGI(NORMXN_ALL(GI,1),NORMXN_ALL(GI,2),rub3, &
                    DXDLX,DYDLX,DZDLX, DXDLY,DYDLY,1.0, &
                    NORMX_ALL(1),NORMX_ALL(2),1.0)
      END DO

    elseif ( ndim==3 ) then
      DO GI=1,SNGI

        DXDLX=0.
        DXDLY=0.
        DYDLX=0.
        DYDLY=0.
        DZDLX=0.
        DZDLY=0.

        DO SL=1,SNLOC
          DXDLX=DXDLX + SNLX(GI,1,SL)*XSL_ALL(1,SL)
          DXDLY=DXDLY + SNLX(GI,2,SL)*XSL_ALL(1,SL)
          DYDLX=DYDLX + SNLX(GI,1,SL)*XSL_ALL(2,SL)
          DYDLY=DYDLY + SNLX(GI,2,SL)*XSL_ALL(2,SL)
          DZDLX=DZDLX + SNLX(GI,1,SL)*XSL_ALL(3,SL)
          DZDLY=DZDLY + SNLX(GI,2,SL)*XSL_ALL(3,SL)
        END DO

        A = DYDLX*DZDLY - DYDLY*DZDLX
        B = DXDLX*DZDLY - DXDLY*DZDLX
        C = DXDLX*DYDLY - DXDLY*DYDLX

        DETJ=SQRT( A**2 + B**2 + C**2)
        ! inv_jac(1,1)=DXDLX; inv_jac(1,2)=DXDLY; inv_jac(1,3)=DXDLZ
        ! inv_jac(2,1)=DyDLX; inv_jac(2,2)=DyDLY; inv_jac(2,3)=DyDLZ
        ! inv_jac(3,1)=DzDLX; inv_jac(3,2)=DzDLY; inv_jac(3,3)=DzDLZ
        ! inv_jac=inv_jac/detj
        SDETWE(GI)=DETJ*SWEIGHT(GI)
        SAREA=SAREA+SDETWE(GI)

        ! Calculate the normal at the Gauss pts...
        ! Perform x-product. N=T1 x T2
        CALL NORMGI(NORMXN_ALL(GI,1),NORMXN_ALL(GI,2),NORMXN_ALL(GI,3), &
                    DXDLX,DYDLX,DZDLX, DXDLY,DYDLY,DZDLY, &
                    NORMX_ALL(1),NORMX_ALL(2),NORMX_ALL(3))
      END DO
    end if

    RETURN
  END SUBROUTINE det_snlx_all


  ! @brief> :: this subroutine calculates detwei and nx of structured elements based on detwei of parent un_ele
  subroutine semi_tri_det_nlx(meshlist, un_ele, x_loc_un_ele, n, nlx, weight, ndim,&
                              nloc, ngi, INV_JAC, n_split)
    implicit none
    integer, intent( in ) :: ndim,nloc,ngi, un_ele

    type(Mesh), dimension(:), intent(inout) :: meshList
    real, DIMENSION( ndim,nloc ), intent( in ) :: x_loc_un_ele
    REAL, DIMENSION( ngi, nloc ), intent( in ) :: N
    REAL, DIMENSION( ngi, ndim, nloc ), intent( in ) :: nlx
    REAL, DIMENSION( ngi ), intent( in ) :: WEIGHT
    integer, intent(in) :: n_split
    ! REAL, DIMENSION( ngi ), intent( inout ) :: DETWEI
    ! REAL, DIMENSION( ngi, ndim, nloc ), intent( inout ) :: nx
    REAL, DIMENSION( ngi, ndim, ndim ), intent( inout ):: INV_JAC

    call tri_det_nlx( x_loc_un_ele, n, nlx, meshList(un_ele)%nx, meshList(un_ele)%detwei,&
                                                        weight, ndim, nloc, ngi, INV_JAC )

    meshList(un_ele)%detwei = meshList(un_ele)%detwei / 4**n_split
    meshList(un_ele)%nx = meshList(un_ele)%nx*2**n_split
  end subroutine semi_tri_det_nlx



  ! @brief> :: this subroutine calculates detwei and nx of structured elements based on detwei of parent un_ele
  subroutine semi_tri_det_nlx_multigrid(meshlist, multi_levels, un_ele, x_loc_un_ele, n, nlx, weight, ndim,&
                              nloc, ngi, INV_JAC, n_split)
    implicit none
    integer, intent( in ) :: ndim,nloc,ngi, un_ele, multi_levels

    type(Mesh), dimension(:), intent(inout) :: meshList
    real, DIMENSION( ndim,nloc ), intent( in ) :: x_loc_un_ele
    REAL, DIMENSION( ngi, nloc ), intent( in ) :: N
    REAL, DIMENSION( ngi, ndim, nloc ), intent( in ) :: nlx
    REAL, DIMENSION( ngi ), intent( in ) :: WEIGHT
    integer, intent(in) :: n_split
    ! REAL, DIMENSION( ngi ), intent( inout ) :: DETWEI
    ! REAL, DIMENSION( ngi, ndim, nloc ), intent( inout ) :: nx
    REAL, DIMENSION( ngi, ndim, ndim ), intent( inout ):: INV_JAC
    ! local vbl
    integer :: ilevel

    do ilevel=1,multi_levels
      call tri_det_nlx( x_loc_un_ele, n, nlx, meshList(un_ele)%scaling_var(ilevel)%nx, meshList(un_ele)%scaling_var(ilevel)%detwei,&
      weight, ndim, nloc, ngi, INV_JAC )
      meshList(un_ele)%scaling_var(ilevel)%detwei = meshList(un_ele)%scaling_var(ilevel)%detwei / 4**(n_split-ilevel+1)
      meshList(un_ele)%scaling_var(ilevel)%nx = meshList(un_ele)%scaling_var(ilevel)%nx*2**(n_split-ilevel+1)
    end do
  end subroutine semi_tri_det_nlx_multigrid



! brief>:: It adds sdetwei of nested structured ele and snorm of parental un_ele to the meshlist class
subroutine semi_det_snlx(meshlist, n_split, nface, ndim, sndim, nloc, sngi,sweight, sn,face_sn,&
                                  face_snlx,un_ele)
  implicit none
!global vbl
type(Mesh), dimension(:), intent(inout) :: meshList
integer, intent(in) :: nface, ndim, nloc, sngi, sndim, n_split, un_ele
real, target, intent( in ) :: face_snlx(sngi,sndim,nloc,nface),face_sn(sngi,nloc,nface)
real, intent(in) :: sweight(sngi)
real, pointer :: sn(:,:)
!local vbl
integer :: iface, iloc, idim
real :: xsgi(sngi,ndim), norm(ndim),sarea, dummy(sngi,ndim), dummy2(sngi)


  do iface=1,nface
    xsgi=0.0
    sn => face_sn(:,:,iface)
    do iloc=1,nloc ! use all of the nodes not just the surface nodes.
      do idim=1,ndim
        xsgi(:,idim)  = xsgi(:,idim)  + sn(:,iloc)*sngl(meshList(un_ele)%X(idim,iloc))
      end do
    end do

    do idim=1,ndim
       norm(idim) = sum(xsgi(:,idim))/real(sngi) - sum(meshList(un_ele)%X(idim,:))/real(nloc)
    end do

    call det_snlx_all( nloc, sngi, sndim, ndim, sngl(meshList(un_ele)%X), face_sn(:,:,iface), face_snlx(:,:,:,iface)&
                      , sweight, meshList(un_ele)%sdetwei(:,iface), sarea, meshList(un_ele)%snorm(:,:,iface), norm )


    meshList(un_ele)%sdetwei(:,iface) = meshList(un_ele)%sdetwei(:,iface)/2**n_split
  end do
  ! I know it is not a good implementatino but I had to do it due to different face numbering between un_ele and str_ele
  dummy = meshList(un_ele)%snorm(:,:,2)
  meshList(un_ele)%snorm(:,:,2) = meshList(un_ele)%snorm(:,:,3)
  meshList(un_ele)%snorm(:,:,3) = dummy

  dummy2 = meshList(un_ele)%sdetwei(:,2)
  meshList(un_ele)%sdetwei(:,2) = meshList(un_ele)%sdetwei(:,3)
  meshList(un_ele)%sdetwei(:,3) = dummy2


end subroutine semi_det_snlx



! brief>:: It adds sdetwei of nested structured ele and snorm of parental un_ele to the meshlist class
subroutine semi_det_snlx_multigrid(meshlist, multi_levels, n_split, nface, ndim, sndim, nloc, sngi,sweight, sn,face_sn,&
                                  face_snlx,un_ele)
  implicit none
!global vbl
type(Mesh), dimension(:), intent(inout) :: meshList
integer, intent(in) :: nface, ndim, nloc, sngi, sndim, n_split, un_ele, multi_levels
real, target, intent( in ) :: face_snlx(sngi,sndim,nloc,nface),face_sn(sngi,nloc,nface)
real, intent(in) :: sweight(sngi)
real, pointer :: sn(:,:)
!local vbl
integer :: iface, iloc, idim, ilevel
real :: xsgi(sngi,ndim), norm(ndim),sarea, dummy(sngi,ndim), dummy2(sngi)


  do iface=1,nface
    xsgi=0.0
    sn => face_sn(:,:,iface)
    do iloc=1,nloc ! use all of the nodes not just the surface nodes.
      do idim=1,ndim
        xsgi(:,idim)  = xsgi(:,idim)  + sn(:,iloc)*sngl(meshList(un_ele)%X(idim,iloc))
      end do
    end do

    do idim=1,ndim
       norm(idim) = sum(xsgi(:,idim))/real(sngi) - sum(meshList(un_ele)%X(idim,:))/real(nloc)
    end do

    do ilevel=1,multi_levels
      call det_snlx_all( nloc, sngi, sndim, ndim, sngl(meshList(un_ele)%X), face_sn(:,:,iface), face_snlx(:,:,:,iface)&
                , sweight, meshList(un_ele)%scaling_var(ilevel)%sdetwei(:,iface), sarea, meshList(un_ele)%snorm(:,:,iface), norm )


      meshList(un_ele)%scaling_var(ilevel)%sdetwei(:,iface) = meshList(un_ele)%scaling_var(ilevel)%sdetwei(:,iface)&
                                                                                      /2**(n_split-ilevel+1)
    end do
  end do
  ! I know it is not a good implementatino but I had to do it due to different face numbering between un_ele and str_ele
  dummy = meshList(un_ele)%snorm(:,:,2)
  meshList(un_ele)%snorm(:,:,2) = meshList(un_ele)%snorm(:,:,3)
  meshList(un_ele)%snorm(:,:,3) = dummy

  dummy2 = meshList(un_ele)%sdetwei(:,2)
  meshList(un_ele)%sdetwei(:,2) = meshList(un_ele)%sdetwei(:,3)
  meshList(un_ele)%sdetwei(:,3) = dummy2


end subroutine semi_det_snlx_multigrid




!@brief>:: this subroutine return nx and position of structured ele in an unstructured ele
subroutine semi_get_nx_pos(irow, ipos, orientation, n_split, ele, nx, updown)
  implicit none
  ! global vbls
  integer, intent(in):: n_split, ele!, ngi, ndim, nloc
  REAL, DIMENSION( :,:,: ), intent( inout ) :: nx

  ! local vbls
  integer, INTENT(OUT):: irow, ipos
  integer :: orientation, updown


  nx = nx *((1-orientation)+orientation*(-1))
  call get_str_info(n_split, ele, irow, ipos, orientation)
  nx = nx *((1-orientation)+orientation*(-1))
  updown = orientation
  if ( orientation==0 ) updown=-1

end subroutine semi_get_nx_pos






  REAL FUNCTION RGPTWE(IG,ND,WEIGHT)
      IMPLICIT NONE
      !     NB If WEIGHT is TRUE in function RGPTWE then return the Gauss-pt weight
      !     else return the Gauss-pt.
      !     NB there are ND Gauss points we are looking for either the
      !     weight or the x-coord of the IG'th Gauss point.
      INTEGER IG,ND
      LOGICAL WEIGHT

        IF(WEIGHT) THEN
           GO TO (10,20,30,40,50,60,70,80,90,100) ND
           !     +++++++++++++++++++++++++++++++
           !     For N=1 +++++++++++++++++++++++
    10     CONTINUE
           RGPTWE=2.0
           GO TO 1000
           !     For N=2 +++++++++++++++++++++++
    20     CONTINUE
           RGPTWE=1.0
           GO TO 1000
           ! For N=3 +++++++++++++++++++++++
    30     CONTINUE
           GO TO (11,12,11) IG
    11     RGPTWE= 0.555555555555556
           GO TO 1000
    12     RGPTWE= 0.888888888888889
           GO TO 1000
           ! For N=4 +++++++++++++++++++++++
    40     CONTINUE
           GO TO (21,22,22,21) IG
    21     RGPTWE= 0.347854845137454
           GO TO 1000
    22     RGPTWE= 0.652145154862546
           GO TO 1000
           ! For N=5 +++++++++++++++++++++++
    50     CONTINUE
           GO TO (31,32,33,32,31) IG
    31     RGPTWE= 0.236926885056189
           GO TO 1000
    32     RGPTWE= 0.478628670499366
           GO TO 1000
    33     RGPTWE= 0.568888888888889
           GO TO 1000
           ! For N=6 +++++++++++++++++++++++
    60     CONTINUE
           GO TO (41,42,43,43,42,41) IG
    41     RGPTWE= 0.171324492379170
           GO TO 1000
    42     RGPTWE= 0.360761573048139
           GO TO 1000
    43     RGPTWE= 0.467913934572691
           GO TO 1000
           ! For N=7 +++++++++++++++++++++++
    70     CONTINUE
           GO TO (51,52,53,54,53,52,51) IG
    51     RGPTWE= 0.129484966168870
           GO TO 1000
    52     RGPTWE= 0.279705391489277
           GO TO 1000
    53     RGPTWE= 0.381830050505119
           GO TO 1000
    54     RGPTWE= 0.417959183673469
           GO TO 1000
           ! For N=8 +++++++++++++++++++++++
    80     CONTINUE
           GO TO (61,62,63,64,64,63,62,61) IG
    61     RGPTWE= 0.101228536290376
           GO TO 1000
    62     RGPTWE= 0.222381034453374
           GO TO 1000
    63     RGPTWE= 0.313706645877877
           GO TO 1000
    64     RGPTWE= 0.362683783378362
           GO TO 1000
           ! For N=9 +++++++++++++++++++++++
    90     CONTINUE
           GO TO (71,72,73,74,75,74,73,72,71) IG
    71     RGPTWE= 0.081274388361574
           GO TO 1000
    72     RGPTWE= 0.180648160694857
           GO TO 1000
    73     RGPTWE= 0.260610696402935
           GO TO 1000
    74     RGPTWE= 0.312347077040003
           GO TO 1000
    75     RGPTWE= 0.330239355001260
           GO TO 1000
           ! For N=10 +++++++++++++++++++++++
    100    CONTINUE
           GO TO (81,82,83,84,85,85,84,83,82,81) IG
    81     RGPTWE= 0.066671344308688
           GO TO 1000
    82     RGPTWE= 0.149451349150581
           GO TO 1000
    83     RGPTWE= 0.219086362515982
           GO TO 1000
    84     RGPTWE= 0.269266719309996
           GO TO 1000
    85     RGPTWE= 0.295524224714753
           !
    1000   CONTINUE
        ELSE
           GO TO (210,220,230,240,250,260,270,280,290,200) ND
           ! +++++++++++++++++++++++++++++++
           ! For N=1 +++++++++++++++++++++++ THE GAUSS POINTS...
    210    CONTINUE
           RGPTWE=0.0
           GO TO 2000
           ! For N=2 +++++++++++++++++++++++
    220    CONTINUE
           RGPTWE= 0.577350269189626
           GO TO 2000
           ! For N=3 +++++++++++++++++++++++
    230    CONTINUE
           GO TO (211,212,211) IG
    211    RGPTWE= 0.774596669241483
           GO TO 2000
    212    RGPTWE= 0.0
           GO TO 2000
           ! For N=4 +++++++++++++++++++++++
    240    CONTINUE
           GO TO (221,222,222,221) IG
    221    RGPTWE= 0.861136311594953
           GO TO 2000
    222    RGPTWE= 0.339981043584856
           GO TO 2000
           ! For N=5 +++++++++++++++++++++++
    250    CONTINUE
           GO TO (231,232,233,232,231) IG
    231    RGPTWE= 0.906179845938664
           GO TO 2000
    232    RGPTWE= 0.538469310105683
           GO TO 2000
    233    RGPTWE= 0.0
           GO TO 2000
           ! For N=6 +++++++++++++++++++++++
    260    CONTINUE
           GO TO (241,242,243,243,242,241) IG
    241    RGPTWE= 0.932469514203152
           GO TO 2000
    242    RGPTWE= 0.661209386466265
           GO TO 2000
    243    RGPTWE= 0.238619186083197
           GO TO 2000
           ! For N=7 +++++++++++++++++++++++
    270    CONTINUE
           GO TO (251,252,253,254,253,252,251) IG
    251    RGPTWE= 0.949107912342759
           GO TO 2000
    252    RGPTWE= 0.741531185599394
           GO TO 2000
    253    RGPTWE= 0.405845151377397
           GO TO 2000
    254    RGPTWE= 0.0
           GO TO 2000
           ! For N=8 +++++++++++++++++++++++
    280    CONTINUE
           GO TO (261,262,263,264,264,263,262,261) IG
    261    RGPTWE= 0.960289856497536
           GO TO 2000
    262    RGPTWE= 0.796666477413627
           GO TO 2000
    263    RGPTWE= 0.525532409916329
           GO TO 2000
    264    RGPTWE= 0.183434642495650
           GO TO 2000
           ! For N=9 +++++++++++++++++++++++
    290    CONTINUE
           GO TO (271,272,273,274,275,274,273,272,271) IG
    271    RGPTWE= 0.968160239507626
           GO TO 2000
    272    RGPTWE= 0.836031107326636
           GO TO 2000
    273    RGPTWE= 0.613371432700590
           GO TO 2000
    274    RGPTWE= 0.324253423403809
           GO TO 2000
    275    RGPTWE= 0.0
           GO TO 2000
           ! For N=10 +++++++++++++++++++++++
    200    CONTINUE
           GO TO (281,282,283,284,285,285,284,283,282,281) IG
    281    RGPTWE= 0.973906528517172
           GO TO 2000
    282    RGPTWE= 0.865063366688985
           GO TO 2000
    283    RGPTWE= 0.679409568299024
           GO TO 2000
    284    RGPTWE= 0.433395394129247
           GO TO 2000
    285    RGPTWE= 0.148874338981631
           !
    2000   CONTINUE
           IF(IG.LE.INT((ND/2)+0.1)) RGPTWE=-RGPTWE
        ENDIF
  END FUNCTION RGPTWE


  ! normal at GI
  SUBROUTINE NORMGI( NORMXN, NORMYN, NORMZN, &
                     DXDLX, DYDLX, DZDLX, DXDLY, DYDLY, DZDLY, &
                     NORMX, NORMY, NORMZ)
    ! Calculate the normal at the Gauss pts
    ! Perform x-product. N=T1 x T2
    implicit none
    REAL, intent( inout ) :: NORMXN, NORMYN, NORMZN
    REAL, intent( in )    :: DXDLX, DYDLX, DZDLX, DXDLY, DYDLY, DZDLY
    REAL, intent( in )    :: NORMX, NORMY, NORMZ
    ! Local variables
    REAL :: RN, SIRN

    CALL XPROD1( NORMXN, NORMYN, NORMZN, &
                 DXDLX, DYDLX, DZDLX, &
                 DXDLY, DYDLY, DZDLY )

    RN = SQRT( NORMXN**2 + NORMYN**2 + NORMZN**2 )

    SIRN = SIGN( 1.0 / RN, NORMXN * NORMX + NORMYN * NORMY + NORMZN * NORMZ )

    NORMXN = SIRN * NORMXN
    NORMYN = SIRN * NORMYN
    NORMZN = SIRN * NORMZN

    RETURN
  END SUBROUTINE NORMGI


  ! cross product
  SUBROUTINE XPROD1( AX, AY, AZ, &
                     BX, BY, BZ, &
                     CX, CY, CZ )
    implicit none
    REAL, intent( inout ) :: AX, AY, AZ
    REAL, intent( in )    :: BX, BY, BZ, CX, CY, CZ

    ! Perform cross product(x-product). a=b x c
    AX =    BY * CZ - BZ * CY
    AY = -( BX * CZ - BZ * CX )
    AZ =    BX * CY - BY * CX

    RETURN
  END subroutine XPROD1

end module ShapFun
