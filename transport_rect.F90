module transport_rect
  use ShapFun
  use structured_meshgen
  use matrices
  contains

  Subroutine trans_rec(CFL, ndim, direct_solver, njac_its, time, nits, no_ele_row,&
    no_ele_col, nface, nloc, snloc, ngi, sngi, n_s_list_no, u_x, u_y, with_stab )
    implicit none

    ! external vbls
    integer, intent(in) :: ndim, njac_its, nits, no_ele_col, no_ele_row, nface
    integer :: nloc, snloc, ngi, sngi, n_s_list_no
    real, intent(in) :: time, u_x, u_y, CFL
    logical, intent(in) :: direct_solver, with_stab

    !local vbls
    integer :: gi, iface, totele, ele, ele2, ele22
    integer :: s_list_no, s_gi, iloc, jloc, i
    integer:: itime, ntime, idim, sndim, nonodes, mloc, col, row, row2
    integer :: errorflag, i_got_boundary, jac_its, its
    integer:: npoly, ele_type, diff_coe_factor
    integer:: offset_x, offset_y

    real :: sarea, volume, dt, L, dx, dy, r_got_boundary, toler, a_coef
    real:: x_Length=100., y_length=100.

    logical:: with_time_slab, D3
    logical :: LOWQUA

    integer, allocatable :: face_ele(:,:), face_list_no(:,:)
    real, allocatable :: n(:,:), nlx(:,:,:), nx(:,:,:), M(:,:)
    real, allocatable :: nlx_lxx(:,:,:), nlxx(:,:), nlx_nod(:,:,:)
    real, allocatable :: weight(:), detwei(:), sdetwei(:)
    real, allocatable :: sn(:,:),sn2(:,:),snlx(:,:,:),sweight(:), s_cont(:,:)
    real, allocatable :: tnew_loc(:), tnew_loc2(:), told(:,:), tnew(:,:), tnew_nonlin(:,:)
    real, allocatable :: tnew_analytical(:,:)
    real, allocatable :: tnew_nonlin_loc(:)
    real, allocatable :: t_bc(:,:), u_bc(:,:,:)
    real, allocatable :: tnew_gi(:), tnew_xgi(:,:), tnew_sgi(:), tnew_sgi2(:), told_gi(:)
    real, allocatable :: face_sn(:,:,:), face_sn2(:,:,:), face_snlx(:,:,:,:)
    real, allocatable :: u_loc(:,:), u_ele(:,:,:), u_loc2(:,:), x_loc(:,:), ugi(:,:)
    real, allocatable :: x_all(:,:,:)
    real, allocatable :: xsgi(:,:), usgi(:,:), usgi2(:,:), income(:), snorm(:,:), norm(:)
    real, allocatable :: mass_ele(:,:), mat_loc(:,:), mass_ele_inv(:,:), mat_loc_inv(:,:)
    real, allocatable :: rhs_loc(:), ml_ele(:), rhs_jac(:), mass_t_new(:)
    real, allocatable :: SN_orig(:,:),SNLX_orig(:,:,:), inv_jac(:, :, : ), mat_diag_approx(:)
    real, allocatable :: a_star(:,:), p_star(:), rgi(:), diff_coe(:), told_loc(:), stab(:,:), mat_tnew(:), mass_told(:)
    real, allocatable:: L1(:), L2(:), L3(:), L4(:)
    real, allocatable:: NLX_ALL(:,:,:)

    sndim = ndim-1 ! dimensions of the surface element coordinates.
    toler = 0.00000000001
    totele = no_ele_row * no_ele_col
    dx = x_length/no_ele_row
    dy = y_length/no_ele_col
    dt = CFL*dx

    allocate(n( ngi, nloc ), nx( ngi, ndim, nloc ), nlx( ngi, ndim, nloc ), M(ngi,nloc))
    allocate(weight(ngi), detwei(ngi), sdetwei(sngi))
    allocate(face_list_no( nface, totele), face_ele(nface,totele))
    allocate(sn(sngi,nloc),sn2(sngi,nloc),snlx(sngi,sndim,nloc),sweight(sngi))
    allocate(SN_orig(sngi,snloc),SNLX_orig(sngi,sndim,snloc), s_cont(sngi,ndim) )
    allocate(u_loc(ndim,nloc), u_loc2(ndim,nloc), ugi(ngi,ndim), u_ele(ndim,nloc,totele))
    allocate(xsgi(sngi,ndim), usgi(sngi,ndim), usgi2(sngi,ndim), income(sngi), snorm(sngi,ndim), norm(ndim))
    allocate(tnew_loc(nloc), tnew_loc2(nloc), told(nloc,totele), tnew(nloc,totele), tnew_nonlin(nloc,totele))
    allocate(tnew_analytical(nloc,totele))
    allocate(tnew_nonlin_loc(nloc))
    allocate(t_bc(nloc,totele), u_bc(ndim,nloc,totele))
    allocate(tnew_gi(ngi), tnew_xgi(ngi,ndim), tnew_sgi(sngi), tnew_sgi2(sngi), told_gi(ngi))
    allocate(face_sn(sngi,nloc,nface), face_sn2(sngi,nloc,n_s_list_no), face_snlx(sngi,sndim,nloc,nface) )
    allocate(x_all(ndim,nloc,totele), x_loc(ndim,nloc))
    allocate(mass_ele(nloc,nloc), mat_loc_inv(nloc,nloc), rhs_loc(nloc), ml_ele(nloc))
    allocate(rhs_jac(nloc), mass_t_new(nloc), inv_jac(ngi, ndim, ndim ) )
    allocate(a_star(ngi,nloc), p_star(ngi), rgi(ngi), diff_coe(ngi), told_loc(nloc), stab(nloc, nloc))
    allocate(mat_tnew(nloc), mat_diag_approx(nloc), mass_told(nloc))

    t_bc(:,:)=0.0 ! this contains the boundary conditions just outside the domain
    u_bc(:,:,:)=0.0 ! this contains the boundary conditions on velocity just outside the domain
    ! 1D initial conditions
    tnew(:,:) = 0.0
    tnew(:,no_ele_row/5:no_ele_row/2) = 1.0 ! a bit off - is this correct -only correct for 1D?

    ! ! 2D initial conditions
    ! do i=1,no_ele_col/2
    !   tnew(:,i*no_ele_row+no_ele_row/5:i*no_ele_row+no_ele_row/2) = 1.0
    ! end do

    u_ele(:,:,:) = 0
    u_ele(1,:,1:totele) = u_x ! suggest this should be unity
    u_ele(2,:,1:totele) = u_y ! suggest this should be unity
    u_bc(:,:,:)=u_ele(:,:,:) ! this contains the boundary conditions on velocity just outside the domain
    diff_coe_factor = 1



    ! just for error analysis
    ! it moves the initial condition with offset to the place that it should be
    tnew_analytical(:,:) = 0.0
    offset_x = int(u_x * dt*ntime*no_ele_row/x_Length+1)
    offset_y = int(u_y * dt*ntime*no_ele_row/y_Length+1)

    ! 1D condition
    tnew_analytical(:,offset_x+no_ele_row/5: offset_x+no_ele_row/2)=1.

    ! ! 2D initial conditions
    ! do i=1,no_ele_col/2
    !   tnew_analytical(:,(i+offset_y)*no_ele_row+((no_ele_row/5)+offset_x):&
    !   (i+offset_y)*no_ele_row+((no_ele_row/2)+offset_x)) = 1.0
    ! end do



    LOWQUA=.true.
    call RE2DN4(LOWQUA,NGI,NLOC,MLOC,M,WEIGHT,N,NLX(:,1,:),NLX(:,2,:),  SNGI,SNLOC,sweight,SN_orig,SNLX_orig)
    ! call ele_info(totele, nface, face_ele, no_ele_row, row, row2, &
    !               x_all, dx, dy, ndim, nloc, no_ele_col, col)
    call surface_pointers_sn(nface, sngi, snloc, nloc, ndim, sndim, totele, n_s_list_no, no_ele_row, no_ele_col, &
                 sn_orig, snlx_orig, face_sn, face_snlx, face_sn2, face_list_no)

    do itime=1,ntime
      told = tnew ! prepare for next time step
      tnew_nonlin = tnew
      do its=1,nits
        tnew = tnew_nonlin ! for non-linearity
        do ele=1,totele
          ! volume integration
          x_loc(:,:)=x_all(:,:,ele) ! x_all contains the coordinates of the corner nodes
          call det_nlx( x_loc, n, nlx, nx, detwei, weight, ndim, nloc, ngi, inv_jac )
          !volume=sum(detwei)
      !      do gi=1,ngi
      !         print *,'gi=',gi
      !         print *,'nx(gi,1,:):',nx(gi,1,:)
      !         print *,'nx(gi,2,:):',nx(gi,2,:)
      !      end do
      !      stop 221
          u_loc(:,:) = u_ele(:,:,ele) ! this is the
          !AMIN changed to t_old
          tnew_loc(:) =  tnew(:,ele) ! this is the FEM element - 2nd index is the local node number
          told_loc(:) =  told(:,ele) ! this is the FEM element - 2nd index is the local node number

          ! calculates vel at gi sum = phi_1(gi)*c1 +phi_2(gi)*c2 +phi_3(gi)*c3 +phi_4(gi)*c4
          ! do idim=1,ndim
          !   do gi=1,ngi
          !      !ugi_x(gi,idim)=sum(nx(gi,idim,:)*u_loc(idim,:))
          !      ugi(gi,idim)=sum(n(gi,:)*u_loc(idim,:))
          !   end do
          ! end do
          ! calculates t at gi sum = phi_1(gi)*t1 +phi_2(gi)*t2 +phi_3(gi)*t3 +phi_4(gi)*t4
          do gi=1,ngi
            do idim=1,ndim
               ugi(gi,idim)=sum(n(gi,:)*u_loc(idim,:))
               tnew_xgi(gi,idim)=sum(nx(gi,idim,:)*tnew_loc(:))
            end do
            tnew_gi(gi)=sum(n(gi,:)*tnew_loc(:))
            told_gi(gi)=sum(n(gi,:)*told_loc(:))
            rgi(gi)=(tnew_gi(gi)-told_gi(gi))/dt + sum(ugi(gi,:)*tnew_xgi(gi,:))
      ! a_star
      ! eq 4
            ! a_coef=sum(ugi(gi,:)*txgi(gi,:))/max(toler, sum( txgi(gi,:)**2 ) )
            a_coef=rgi(gi)/max(toler, sum( tnew_xgi(gi,:)**2 ) )
            a_star(gi,:) = a_coef * tnew_xgi(gi,:)
            ! eq 14
            ! p_star(gi) =0.0
            ! do iloc=1,nloc
            !    p_star(gi) = max(p_star(gi), abs(sum( a_star(gi,:)*nx(gi,:,iloc) ))  )
            ! end do
            ! p_star(gi) = 0.25/max(toler, p_star(gi))
            ! eq 23 for P*
            p_star(gi) =0.0
            do idim=1,ndim
              do iloc=1,nloc
                 p_star(gi) = max(p_star(gi), abs(sum( a_star(gi,:)*(inv_jac(gi,:,idim)))  ))
              end do
            end do
            p_star(gi) = min(1/toler ,0.25/p_star(gi) )
            ! eq 18
            diff_coe(gi) = diff_coe_factor*(0.25*rgi(gi)**2 *p_star(gi) /max(toler, sum( tnew_xgi(gi,:)**2 ) ))
          end do
          rhs_loc=0.0 ! the local to element rhs vector.
          mass_ele=0.0 ! initialize mass matric to be zero.
          stab=0.0
          do iloc=1,nloc
            !inod = glob_no(iloc,ele)
            do jloc=1,nloc
              !jnod = glob_no(ele,jloc)

              do gi=1,ngi
                stab(iloc,jloc) = stab(iloc,jloc) + sum(diff_coe(gi)* nx(gi,:,jloc)* nx(gi,:,iloc))* detwei(gi)
                !M(inod,jnod) = M(inod,jnod) + ngi_3p_wei(g)*sh_func(iloc)*sh_func(jloc)*det_jac
                mass_ele(iloc,jloc) = mass_ele(iloc,jloc) + n(gi,iloc)*n(gi,jloc)*detwei(gi)
              end do ! quadrature
            end do ! ijloc

            ! Amin- How can it be lumped mass matrix as it does not have any mass_ele involved?
            ml_ele(iloc)=sum(n(:,iloc)*detwei(:)) ! lumped mass matrix in element (use later for solver).
            do gi=1,ngi
              do idim=1,ndim
                rhs_loc(iloc) = rhs_loc(iloc) + nx(gi,idim,iloc)*ugi(gi,idim)*tnew_gi(gi)*detwei(gi)
              end do
            end do ! quadrature
          end do ! iloc

          !Include the surface integral here:
          do iface = 1,nface
            ele22 = face_ele( iface, ele) ! if ele22 is negative assume on boundary
            i_got_boundary = (sign(1, -ele22) + 1 )/2
            r_got_boundary = real(i_got_boundary)
            ele2 = ele*i_got_boundary + ele22*(1-i_got_boundary)
            ! r_got_boundary=1.0 if we want to use the boundary conditions and have incomming velocity.
            ! r_got_boundary=0.0 not on boundary

            tnew_loc2(:)=tnew(:,ele2) * (1.0-r_got_boundary)     + t_bc(:,ele)  * r_got_boundary
            u_loc2(:,:)=u_ele(:,:,ele2)* (1.0-r_got_boundary)  + u_bc(:,:,ele)* r_got_boundary

            !Surface integral along an element
            ! need to work on this also
            s_list_no = face_list_no( iface, ele) ! correct
            sn(:,:)     = face_sn(:,:,iface)  ! correct
            snlx(:,:,:) = face_snlx(:,:,:,iface)
            sn2(:,:)    = face_sn2(:,:,s_list_no) ! correct

            usgi=0.0; usgi2=0.0; xsgi=0.0; tnew_sgi=0.0; tnew_sgi2=0.0
            do iloc=1,nloc ! use all of the nodes not just the surface nodes.
              do idim=1,ndim

                usgi(:,idim)  = usgi(:,idim)  + sn(:,iloc)*u_loc(idim,iloc)
                usgi2(:,idim) = usgi2(:,idim) + sn2(:,iloc)*u_loc2(idim,iloc)
                xsgi(:,idim)  = xsgi(:,idim)  + sn(:,iloc)*x_loc(idim,iloc)

              end do
              tnew_sgi(:)  = tnew_sgi(:)  + sn(:,iloc)*tnew_loc(iloc)
              tnew_sgi2(:) = tnew_sgi2(:) + sn2(:,iloc)*tnew_loc2(iloc)
            end do
            !        usgi=0.0
            !        usgi2=0.0

            ! this is the approximate normal direction...
            do idim=1,ndim
               norm(idim) = sum(xsgi(:,idim))/real(sngi) - sum(x_loc(idim,:))/real(nloc)
            end do
            ! start to do the integration
            call det_snlx_all( nloc, sngi, sndim, ndim, x_loc, sn, snlx, sweight, sdetwei, sarea, snorm, norm )
            !      print *,'iface,sum(sdetwei(:)):',iface,sum(sdetwei(:))

            ! income=1 if info is comming from neighbouring element.
            do s_gi=1,sngi
              income(s_gi)=0.5 + 0.5*sign(1.0, -sum(snorm(s_gi,:)*0.5*(usgi(s_gi,:)+usgi2(s_gi,:)))  )
            end do
            !        print *,'iface,income:',iface,income
            !        print *,'snorm(1,:):',snorm(1,:)
            !        print *,'snorm(2,:):',snorm(2,:)
            !        print *,'ele22, i_got_boundary, r_got_boundary, ele, ele2:',ele22, i_got_boundary, r_got_boundary, ele, ele2

            ! sloc_vec=0.0; sloc_vec2=0.0
            do idim=1,ndim
              s_cont(:,idim) = snorm(:,idim)*sdetwei(:) &
                          *( (1.-income(:))* usgi(:,idim)*tnew_sgi(:) + income(:)*usgi2(:,idim)*tnew_sgi2(:) )
            end do

            do iloc=1,nloc
              do idim=1,ndim
                rhs_loc(iloc)  = rhs_loc(iloc)  - sum( sn(:,iloc)*s_cont(:,idim) )
              end do
            end do
          end do ! iface

          if(with_stab) then
            mat_loc= mass_ele + dt*stab
          else
            mat_loc= mass_ele
          end if

          if(direct_solver) then
            ! inverse of the mass matric (nloc,nloc)
            ! AMIN I would say the passed matrix should be mass_ele not mat_loc based on the formula
            call FINDInv(mass_ele, mat_loc_inv, nloc, errorflag)
            !can be only this part instead
            ! do iloc=1,nloc
            !    tnew_nonlin(iloc,ele)=told_loc(iloc) + dt* sum( mat_loc_inv(iloc,:)* rhs_loc(:) ) &
            !                                  + told_loc(iloc)* dt* sum( mat_loc_inv(iloc,:) * stab(iloc,:) )
            ! end do
            do iloc=1,nloc
              mass_told(iloc)=sum( mat_loc(iloc,:) * told_loc(:) )
            end do
            do iloc=1,nloc
              tnew_nonlin(iloc,ele)=sum( mat_loc_inv(iloc,:)* (mass_told(:)+dt*rhs_loc(:)))
            end do

          else ! iterative solver
             do iloc=1,nloc
                rhs_jac(iloc)= sum( mass_ele(iloc,:) * told_loc(:) ) + dt*rhs_loc(iloc)
                mat_diag_approx(iloc) = ml_ele(iloc) + dt * stab(iloc,iloc)
             end do

             tnew_nonlin_loc(:) = tnew_nonlin(:,ele)
             do jac_its=1,njac_its ! jacobi iterations...
                do iloc=1,nloc
                   mat_tnew(iloc)= sum( mat_loc(iloc,:) * tnew_nonlin_loc(:) )
                end do
                do iloc=1,nloc
                   tnew_nonlin_loc(iloc)=  (mat_diag_approx(iloc)*tnew_nonlin_loc(iloc) - mat_tnew(iloc) + rhs_jac(iloc) )&
                                                                  / mat_diag_approx(iloc)
                end do
             end do
             tnew_nonlin(:,ele) = tnew_nonlin_loc(:)
          endif ! endof if(direct_solver) then else

        end do ! do ele=1,totele
        tnew=tnew_nonlin
      end do ! do its=1,nits
    end do ! do itime=1,ntime

! print*, sum(abs(tnew-tnew_analytical))*100/(totele*nloc)

    OPEN(unit=10, file='DG-rectangular_structured')

    ! for 1D results only
      do ele=1,totele
        write(10,*) x_all(1,1,ele), tnew(1,ele)
        write(10,*) x_all(1,2,ele), tnew(2,ele)
      end do

      ! for 2D results only
      ! do ele=1,totele
      !   write(10,*) x_all(1,1,ele), x_all(2,1,ele), tnew(1,ele)
      !   write(10,*) x_all(1,2,ele), x_all(2,2,ele), tnew(2,ele)
      !   write(10,*) x_all(1,3,ele), x_all(2,3,ele), tnew(3,ele)
      !   write(10,*) x_all(1,4,ele), x_all(2,4,ele), tnew(4,ele)
      ! end do
    close(10)

    OPEN(unit=11, file='DG-rectangular_structured_analytical')
    ! for 1D results only
    do ele=1,totele
      write(11,*) x_all(1,1,ele), tnew_analytical(1,ele)
      write(11,*) x_all(1,2,ele), tnew_analytical(2,ele)
      write(11,*) x_all(1,3,ele), tnew_analytical(3,ele)
      write(11,*) x_all(1,4,ele), tnew_analytical(4,ele)
    end do

      ! for 2D results only
      ! do ele=1,totele
      !   write(11,*) x_all(1,1,ele), x_all(2,1,ele), tnew_analytical(1,ele)
      !   write(11,*) x_all(1,2,ele), x_all(2,2,ele), tnew_analytical(2,ele)
      !   write(11,*) x_all(1,3,ele), x_all(2,3,ele), tnew_analytical(3,ele)
      !   write(11,*) x_all(1,4,ele), x_all(2,4,ele), tnew_analytical(4,ele)
      ! end do
    close(11)
      !       print *,'tnew:',tnew
      !       print *,'told:',told
     ! print *,' '
     ! do ele=1,totele
     !    print *,dx*real(ele-1),tnew(1,ele)
     !    print *,dx*real(ele),  tnew(2,ele)
     ! end do
     ! print *,' '
     ! do ele=1,totele
     !    print *,dx*real(ele-1),tnew(1,ele)
     !    print *,dx*real(ele),  tnew(2,ele)
     !    print *,dx*real(ele-1),told(1,ele)
     !    print *,dx*real(ele),  told(2,ele)
     ! end do

    deallocate(face_ele, face_list_no, n, nlx, nx, M, weight, detwei, sdetwei,&
               sn,sn2,snlx,sweight, s_cont, tnew_loc, tnew_loc2, told, tnew, tnew_nonlin,&
               t_bc, u_bc, tnew_gi, tnew_sgi, tnew_sgi2, told_gi, face_sn, face_sn2, face_snlx,&
               u_loc, u_ele, u_loc2, x_loc, ugi, x_all,&
               xsgi, usgi, usgi2, income, snorm, norm, &
               mass_ele, mat_loc_inv, rhs_loc, ml_ele, rhs_jac, mass_t_new,&
               SN_orig,SNLX_orig, inv_jac, mat_diag_approx,  &
               a_star, p_star, rgi, diff_coe, told_loc, stab, mat_tnew)


  end subroutine trans_rec
end module transport_rect
