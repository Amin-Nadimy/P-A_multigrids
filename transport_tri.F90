module transport_tri
  use Structures
  use ShapFun
  use structured_meshgen
  use matrices
  use get_vtk_files
  contains


  ! brief> this is for transport structured tri ele
  ! implicit method
  Subroutine str_implicit(CFL, ndim, direct_solver, njac_its, time, nits, no_ele_row,&
    no_ele_col, nface, dx, dy, nloc, snloc, ngi, sngi, n_s_list_no, u_x, u_y, with_stab, vtk_interval, cell_type )
    implicit none

    ! external vbls
    real, intent(in) :: CFL, time, dx, dy, u_x, u_y
    integer, intent(in) :: ndim, njac_its, nits, no_ele_col, no_ele_row, nface, nloc, snloc, ngi, sngi
    integer, intent(in) :: n_s_list_no, vtk_interval, cell_type
    logical, intent(in) :: direct_solver, with_stab

    ! local vbls
    integer :: gi, iface, totele, ele, ele2, ele22, errorflag, i_got_boundary, jac_its, its
    integer :: s_list_no, s_gi, iloc, jloc, i, npoly, ele_type, vtk_io
    integer:: itime, ntime, idim, sndim, nonodes, mloc, g_iloc, g_jloc, target_ele
    integer, allocatable :: face_ele2(:)!,face_ele(:,:), face_list_no(:,:)

    real :: sarea, volume, dt, L, r_got_boundary, toler, a_coef
    real :: t1, t2, t1_get_shape_funs_spec, t2_get_shape_funs_spec
    real :: t1_tri_ele_info2, t2_tri_ele_info2, time_tri_ele_info2
    real :: t1_tri_det_nlx, t2_tri_det_nlx, time_tri_det_nlx
    real :: t1_det_snlx_all, t2_det_snlx_all, time_det_snlx_all
    real :: lhs, mass_ele, stiff_ele, flux_ele

    character(10) :: start, finish
    logical:: with_time_slab, D3
    type(sparse) :: sparse_mass, sparse_lhs, sparse_flux

    real, allocatable :: n(:,:), nlx(:,:,:), nx(:,:,:), M(:,:)
    real, allocatable :: nlx_lxx(:,:,:), nlxx(:,:), nlx_nod(:,:,:)
    real, allocatable :: weight(:), detwei(:), sdetwei(:)
    real, allocatable :: sn(:,:),sn2(:,:),snlx(:,:,:),sweight(:), s_cont(:,:)
    real, allocatable :: tnew_loc(:), tnew_loc2(:), told(:,:), tnew(:,:), tnew_nonlin(:,:)
    real, allocatable :: tnew_nonlin_loc(:)
    real, allocatable :: t_bc(:,:), u_bc(:,:,:)
    real, allocatable :: tnew_gi(:), tnew_xgi(:,:), tnew_sgi(:), tnew_sgi2(:), told_gi(:)
    real, allocatable :: face_sn(:,:,:), face_sn2(:,:,:), face_snlx(:,:,:,:)
    real, allocatable :: u_loc(:,:), u_ele(:,:,:), u_loc2(:,:), x_loc(:,:), ugi(:,:)
    real, allocatable :: xsgi(:,:), usgi(:,:), usgi2(:,:), income(:), snorm(:,:), norm(:)
    real, allocatable :: mat_loc(:,:), mat_loc_inv(:,:), rhs_loc(:), ml_ele(:), rhs_jac(:), mass_t_new(:)
    real, allocatable :: SN_orig(:,:),SNLX_orig(:,:,:), inv_jac(:, :, : ), mat_diag_approx(:)
    real, allocatable :: a_star(:,:), p_star(:), rgi(:), diff_coe(:), told_loc(:), stab(:,:), mat_tnew(:), mass_told(:)
    real, allocatable :: glob_mass(:,:), glob_lhs(:,:), glob_told_array(:),rhs_array(:),inv_lhs(:,:)

    real, allocatable:: L1(:), L2(:), L3(:), L4(:), tnew_nonlin2(:)
    real, allocatable:: NLX_ALL(:,:,:), glob_flux(:,:), glob_flux2(:,:)



    totele = no_ele_row * no_ele_col *2
    sndim = ndim-1 ! dimensions of the surface element coordinates.
    mloc = 1
    toler = 0.00000000001
    ele_type = 3
    vtk_io = vtk_interval

    allocate(n( ngi, nloc ), nx( ngi, ndim, nloc ), nlx( ngi, ndim, nloc ), M(ngi,nloc))
    allocate(nlx_lxx(ngi,ndim,nloc), nlxx(ngi,nloc), nlx_nod(nloc,ndim,nloc))
    allocate(weight(ngi), detwei(ngi), sdetwei(sngi))
    allocate( face_ele2(nface))!, face_ele(nface,totele),face_list_no( nface, totele))
    allocate(sn(sngi,nloc),sn2(sngi,nloc),snlx(sngi,sndim,nloc))!,sweight(sngi))
    allocate(SN_orig(sngi,snloc),SNLX_orig(sngi,sndim,snloc), s_cont(sngi,ndim) )
    allocate(u_loc(ndim,nloc), u_loc2(ndim,nloc), ugi(ngi,ndim), u_ele(ndim,nloc,totele))
    allocate(xsgi(sngi,ndim), usgi(sngi,ndim), usgi2(sngi,ndim), income(sngi), snorm(sngi,ndim), norm(ndim))
    allocate(tnew_loc(nloc), tnew_loc2(nloc), told(nloc,totele), tnew(nloc,totele), tnew_nonlin(nloc,totele))
    allocate(tnew_nonlin_loc(nloc))
    allocate(t_bc(nloc,totele), u_bc(ndim,nloc,totele))
    allocate(tnew_gi(ngi), tnew_xgi(ngi,ndim), tnew_sgi(sngi), tnew_sgi2(sngi), told_gi(ngi))
    allocate(face_sn(sngi,nloc,nface), face_sn2(sngi,nloc,n_s_list_no), face_snlx(sngi,sndim,nloc,nface) )
    allocate(x_loc(ndim,nloc),tnew_nonlin2(nloc*totele))
    allocate(mat_loc_inv(nloc,nloc), rhs_loc(nloc), ml_ele(nloc))
    allocate(rhs_jac(nloc), mass_t_new(nloc), inv_jac(ngi, ndim, ndim ) )
    allocate(a_star(ngi,nloc), p_star(ngi), rgi(ngi), diff_coe(ngi), told_loc(nloc), stab(nloc, nloc))
    allocate(mat_tnew(nloc), mat_diag_approx(nloc), mass_told(nloc))
    allocate(glob_mass(totele*nloc, totele*nloc), glob_lhs(totele*nloc, totele*nloc))
    allocate(rhs_array(nloc*totele),inv_lhs(totele*nloc,totele*nloc))


    allocate(L1(ngi), L2(ngi), L3(ngi), L4(ngi))
    allocate(NLX_ALL(ngi,ndim,nloc),glob_flux(totele*nloc,totele*nloc), glob_flux2(totele*nloc,totele*nloc))

    ! CPU_TIME(t1)
    ! SYSTEM_CLOCK
    ! character(10) :: time
    ! call DATE_AND_TIME(TIME = time)
    ! print *, 'Hour:  ', time(1)
    ! print *, 'Minute:', time(2)
    ! print *, 'Second:', time(3)

    t_bc(:,:)=0.0 ! this contains the boundary condif ( itime==2 .AND. ele==2 .and. its==1 ) then

    u_bc(:,:,:)=0.0 ! this contains the boundary conditions on velocity just outside the domain
    tnew(:,:) = 0.0

    ! 1D initial conditions
    ! tnew(:,1:2) = 1.0 ! only correct for 1D?
    ! do ele=1,no_ele_col-50
    !   tnew(:,ele*no_ele_row+10:ele*no_ele_row+50) = 1.0 ! just for now, you can use the line above for 1D or lines below for 2D instead.
    ! end do

    ! 2D initial conditions
      tnew(:,1:no_ele_row/5+1) = 1.0
    do i=2,no_ele_col/5
      tnew(:,(i-1)*no_ele_row+1 : no_ele_row*i-(no_ele_row*4/5)) = 1.0
    end do

    u_ele(:,:,:) = 0
    u_ele(1,:,1:totele) = u_x ! suggest this should be unity
    u_ele(2,:,1:totele) = u_y ! suggest this should be unity
    u_bc(:,:,:)=u_ele(:,:,:) ! this contains the boundary conditions on velocity just outside the domain

    dt = CFL*dx
    ! dt = CFL/(u_x/dx+u_y/dy)
    ntime = time/dt


    with_time_slab =.false.
    D3=.false.
    ! IQADRA = ipoly+1
    npoly=1; ele_type= 101

    call TRIQUAold(L1, L2, L3, L4, WEIGHT, D3,NGI)

    call make_sparse_matrix(sparse_lhs, totele, nloc)
    call make_sparse_matrix(sparse_mass, totele, nloc)
    call make_sparse_matrix_flux_str(sparse_flux, no_ele_row, totele, nloc)
    call get_shape_funs_spec(n, nlx, nlx_lxx, nlxx, weight, nlx_nod, &
              nloc, snloc, sngi, ngi, ndim, nface,n_s_list_no, face_sn, face_sn2, face_snlx, &
              sweight, npoly, ele_type, totele)!, face_list_no)
    ! call CPU_TIME(t2_get_shape_funs_spec)

    ! time_tri_ele_info2 = 0.0
    ! time_tri_det_nlx = 0.0
    ! time_det_snlx_all = 0.0
    print*, 'totele = ',totele
    print*, 'ntime = ', ntime
    call CPU_TIME(t1)
    ! call DATE_AND_TIME(TIME = start)

    do itime=1,ntime
      !generating VTK files
      if ( vtk_io <= itime  ) then
        call get_vtk_str_tri(tnew, totele, no_ele_row, no_ele_col, nface, &
                                  dx, dy, nloc, itime, ndim, cell_type)
        vtk_io = vtk_io + vtk_interval
      end if

      told = tnew ! prepare for next time step
      tnew_nonlin = tnew
      do its=1,nits
        glob_flux2 =0.0
        sparse_mass%val=0.
        sparse_lhs%val=0.
        sparse_flux%val=0.0
        tnew = tnew_nonlin ! for non-linearity
        do ele=1,totele
          ! call tri_ele_info2(totele, ele, nface, face_ele2, no_ele_row, &
          !                   x_loc, dx, dy, ndim, nloc, no_ele_col)
          call str_tri_X_nodes(ele, x_loc, ndim, nloc, dx, dy, no_ele_row, no_ele_col)
          ! time_tri_ele_info2 = time_tri_ele_info2 + (t2_tri_ele_info2 - t1_tri_ele_info2)

          ! volume integration
          call tri_det_nlx( x_loc, n, nlx, nx, detwei, weight, ndim, nloc, ngi, inv_jac )
          ! time_tri_det_nlx = time_tri_det_nlx + (t2_tri_det_nlx - t1_tri_det_nlx)

          u_loc = u_ele(:,:,ele) ! this is the
          tnew_loc =  tnew(:,ele) ! this is the FEM element - 2nd index is the local node number
          told_loc =  told(:,ele) ! this is the FEM element - 2nd index is the local node number

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
            p_star(gi) =0.0
            do idim=1,ndim
              do iloc=1,nloc
                 p_star(gi) = max(p_star(gi), abs(sum( a_star(gi,:)*(inv_jac(gi,:,idim)))  ))
              end do
            end do
            p_star(gi) = min(1/toler ,0.25/p_star(gi) )
            ! eq 18
            diff_coe(gi) = 0.25*rgi(gi)**2 *p_star(gi) /max(toler, sum( tnew_xgi(gi,:)**2 ) )
          end do
          rhs_loc=0.0 ! the local to element rhs vector.
          stab=0.0
          do iloc=1,nloc
            g_iloc = glob_no(ele, nloc, iloc)
            do jloc=1,nloc
              g_jloc = glob_no(ele, nloc, jloc)
              mass_ele=0.0
              stiff_ele = 0.0
              lhs=0.0
              do gi=1,ngi
                stab(iloc,jloc) = stab(iloc,jloc) + sum(diff_coe(gi)* nx(gi,:,jloc)* nx(gi,:,iloc))* detwei(gi)
                mass_ele = mass_ele + n(gi,iloc)*n(gi,jloc)*detwei(gi)/dt
                do idim=1,ndim
                  stiff_ele = stiff_ele + nx(gi,idim,iloc)*ugi(gi,idim)*detwei(gi)*n(gi,jloc)
                end do
              end do ! quadrature
              call add_to_CSR(sparse_mass,g_iloc, g_jloc, mass_ele)
              lhs = mass_ele-stiff_ele
              call add_to_CSR(sparse_lhs, g_iloc, g_jloc, lhs)
            end do ! iloc
          end do ! iloc


          !Include the surface integral here:
          do iface = 1,nface
            call tri_ele_info2(totele, ele, iface, ele22, no_ele_row)
            i_got_boundary = (sign(1, -ele22) + 1 )/2 ! 1 means it is on the boundary
            r_got_boundary = real(i_got_boundary)
            ele2 = ele*i_got_boundary + ele22*(1-i_got_boundary)
            ! r_got_boundary=1.0 if we want to use the boundary conditions and have incomming velocity.
            ! r_got_boundary=0.0 not on boundary

            tnew_loc2(:)=tnew(:,ele2) * (1.0-r_got_boundary)   + t_bc(:,ele)  * r_got_boundary
            u_loc2(:,:)=u_ele(:,:,ele2)* (1.0-r_got_boundary)  + u_bc(:,:,ele)* r_got_boundary

            !Surface integral along an element
            sn   = face_sn(:,:,iface)
            snlx = face_snlx(:,:,:,iface)
            sn2  = face_sn2(:,:,iface)


            usgi=0.0; usgi2=0.0; xsgi=0.0; tnew_sgi=0.0; tnew_sgi2=0.0
            do iloc=1,nloc ! use all of the nodes not just the surface nodes.
              do idim=1,ndim

                usgi(:,idim)  = usgi(:,idim)  + sn(:,iloc)*u_loc(idim,iloc)
                usgi2(:,idim) = usgi2(:,idim) + sn2(:,iloc)*u_loc2(idim,iloc)
                xsgi(:,idim)  = xsgi(:,idim)  + sn(:,iloc)*x_loc(idim,iloc)

              end do
              tnew_sgi  = tnew_sgi(:)  + sn(:,iloc)*tnew_loc(iloc)
              tnew_sgi2 = tnew_sgi2(:) + sn2(:,iloc)*tnew_loc2(iloc)
            end do
            !        usgi=0.0
            !        usgi2=0.0

            ! this is the approximate normal direction...
            do idim=1,ndim
               norm(idim) = sum(xsgi(:,idim))/real(sngi) - sum(x_loc(idim,:))/real(nloc)
            end do
            ! start to do the integration
            call det_snlx_all( nloc, sngi, sndim, ndim, x_loc, sn, snlx, sweight, sdetwei, sarea, snorm, norm )

            ! income=1 if info is comming from neighbouring element.
            do s_gi=1,sngi
              income(s_gi)=0.5 + 0.5*sign(1.0, -sum(snorm(s_gi,:)*0.5*(usgi(s_gi,:)+usgi2(s_gi,:)))  )
            end do

            ! sloc_vec=0.0; sloc_vec2=0.0
            do idim=1,ndim
              s_cont(:,idim) = snorm(:,idim)*sdetwei(:) &
                          *( (1.-income(:))* usgi(:,idim) + income(:)*usgi2(:,idim) )
            end do

            target_ele = int( (1.-income(1))*ele + income(1)*ele22 )
            ! if ( target_ele==0 ) then ! it happens when ele is on the boundary
            !   target_ele=ele
            ! end if
            do iloc=1,nloc
              g_iloc=glob_no(ele,nloc,iloc)
              do jloc=1,nloc
                g_jloc=glob_no(target_ele, nloc, jloc)
                flux_ele =0.0
                do idim=1,ndim
                  glob_flux2(g_iloc,g_jloc) = glob_flux2(g_iloc,g_jloc) +sum( snorm(:,idim)*sdetwei(:)*sn(:,iloc)&
                                                                *((1.-income(:))*sn(:,jloc)*usgi(:,idim) &
                                                                +     income(:)*sn2(:,jloc)*usgi2(:,idim)))
                  flux_ele = flux_ele +sum( snorm(:,idim)*sdetwei(:)*sn(:,iloc) &
                                                                *((1.-income(:))*sn(:,jloc)*usgi(:,idim) &
                                                                +     income(:)*sn2(:,jloc)*usgi2(:,idim)))
                end do
                call add_to_CSR_flux(sparse_flux, g_iloc, g_jloc,  flux_ele)
              end do ! jloc
            end do ! iloc
          end do ! iface

          ! if (with_stab) then
            ! mat_loc= mass_ele + dt*stab
          ! else
          !   mat_loc= mass_ele
          ! end if

        end do ! do ele=1,totele
        call gen_global_matrix(totele, nloc, sparse_mass, glob_mass)
        call gen_global_matrix(totele, nloc, sparse_lhs, glob_lhs)
        call gen_global_matrix(totele, nloc, sparse_flux, glob_flux)
        call told_to_array(totele, nloc, told, glob_told_array)
        rhs_array = matmul(glob_mass, glob_told_array)
        glob_lhs = glob_lhs + glob_flux2 !glob_flux

        call FINDInv(glob_lhs, inv_lhs, nloc*totele, errorflag)
        tnew_nonlin2 = matmul(inv_lhs,rhs_array )
        call array_to_told(totele, nloc, tnew_nonlin2, tnew_nonlin)
        tnew=tnew_nonlin
      end do ! do its=1,nits
      print*, 'srtuctured ',itime
    end do ! do itime=1,ntime

    !runing time
    call CPU_TIME(t2)
    print*, 'cpu_time for time_loop = ', t2-t1
    print*, 'cpu_time tri_ele_info2 = ', time_tri_ele_info2
    print*, 'cpu_time tri_det_nlx = ', time_tri_det_nlx
    print*, ''
    print*, 'cpu_time det_snlx_all = ', time_det_snlx_all
    print*, ''


    deallocate(face_ele2, n, nlx, nx, M,&
               nlx_lxx, nlxx, nlx_nod, weight, detwei, sdetwei, sn,sn2,snlx,sweight, s_cont,&
               tnew_loc, tnew_loc2, told, tnew, tnew_nonlin, tnew_nonlin_loc, t_bc, u_bc,&
               tnew_gi, tnew_xgi, tnew_sgi, tnew_sgi2, told_gi, face_sn, face_sn2, face_snlx,&
               u_loc, u_ele, u_loc2, x_loc, ugi, xsgi, usgi, usgi2, income, snorm, norm,&
               mat_loc_inv, rhs_loc, ml_ele, rhs_jac, mass_t_new,&
               SN_orig,SNLX_orig, inv_jac, mat_diag_approx,tnew_nonlin2,&
               a_star, p_star, rgi, diff_coe, told_loc, stab, mat_tnew, mass_told,glob_flux,glob_flux2,&
               L1, L2, L3, L4, NLX_ALL, glob_mass, glob_lhs,glob_told_array,rhs_array,inv_lhs)
  end subroutine str_implicit



  ! brief> this is for transport structured trianggular elements
  ! explicit method
  Subroutine str_explicit(CFL, ndim, direct_solver, njac_its, time, nits, no_ele_row,&
    no_ele_col, nface, dx, dy, nloc, snloc, ngi, sngi, n_s_list_no, u_x, u_y, with_stab, vtk_interval, cell_type )
    implicit none

    ! external vbls
    real, intent(in) :: CFL, time, dx, dy, u_x, u_y
    integer, intent(in) :: ndim, njac_its, nits, no_ele_col, no_ele_row, nface, nloc, snloc, ngi, sngi
    integer, intent(in) :: n_s_list_no, vtk_interval, cell_type
    logical, intent(in) :: direct_solver, with_stab

    ! local vbls
    integer :: gi, iface, totele, ele, ele2, ele22, errorflag, i_got_boundary, jac_its, its
    integer :: s_list_no, s_gi, iloc, jloc, i, npoly, ele_type, vtk_io
    integer:: itime, ntime, idim, sndim, nonodes, mloc
    integer, allocatable :: face_ele2(:)!,face_ele(:,:), face_list_no(:,:)

    real :: sarea, volume, dt, L, r_got_boundary, toler, a_coef
    real :: t1, t2, t1_get_shape_funs_spec, t2_get_shape_funs_spec
    real :: t1_tri_ele_info2, t2_tri_ele_info2, time_tri_ele_info2
    real :: t1_tri_det_nlx, t2_tri_det_nlx, time_tri_det_nlx
    real :: t1_det_snlx_all, t2_det_snlx_all, time_det_snlx_all,face_nodes(2,3)

    character(10) :: start, finish
    logical:: with_time_slab, D3

    real, allocatable :: n(:,:), nlx(:,:,:), nx(:,:,:), M(:,:)
    real, allocatable :: nlx_lxx(:,:,:), nlxx(:,:), nlx_nod(:,:,:)
    real, allocatable :: weight(:), detwei(:), sdetwei(:)
    real, allocatable :: sn(:,:),sn2(:,:),snlx(:,:,:),sweight(:), s_cont(:,:)
    real, allocatable :: tnew_loc(:), tnew_loc2(:), told(:,:), tnew(:,:), tnew_nonlin(:,:)
    real, allocatable :: tnew_nonlin_loc(:)
    real, allocatable :: t_bc(:,:), u_bc(:,:,:)
    real, allocatable :: tnew_gi(:), tnew_xgi(:,:), tnew_sgi(:), tnew_sgi2(:), told_gi(:)
    real, allocatable :: face_sn(:,:,:), face_sn2(:,:,:), face_snlx(:,:,:,:)
    real, allocatable :: u_loc(:,:), u_ele(:,:,:), u_loc2(:,:), x_loc(:,:), ugi(:,:)
    real, allocatable :: xsgi(:,:), usgi(:,:), usgi2(:,:), income(:), snorm(:,:), norm(:)
    real, allocatable :: mass_ele(:,:), mat_loc(:,:), mat_loc_inv(:,:), rhs_loc(:), ml_ele(:), rhs_jac(:), mass_t_new(:)
    real, allocatable :: SN_orig(:,:),SNLX_orig(:,:,:), inv_jac(:, :, : ), mat_diag_approx(:)
    real, allocatable :: a_star(:,:), p_star(:), rgi(:), diff_coe(:), told_loc(:), stab(:,:), mat_tnew(:), mass_told(:)
    real, allocatable :: stnew_loc2(:,:)

    real, allocatable:: L1(:), L2(:), L3(:), L4(:)
    real, allocatable:: NLX_ALL(:,:,:)



    totele = no_ele_row * no_ele_col *2
    sndim = ndim-1 ! dimensions of the surface element coordinates.
    mloc = 1
    toler = 0.00000000001
    ele_type = 3
    vtk_io = vtk_interval

    allocate(n( ngi, nloc ), nx( ngi, ndim, nloc ), nlx( ngi, ndim, nloc ), M(ngi,nloc))
    allocate(nlx_lxx(ngi,ndim,nloc), nlxx(ngi,nloc), nlx_nod(nloc,ndim,nloc))
    allocate(weight(ngi), detwei(ngi), sdetwei(sngi))
    allocate( face_ele2(nface))!, face_ele(nface,totele),face_list_no( nface, totele))
    allocate(sn(sngi,snloc),sn2(sngi,snloc),snlx(sngi,sndim,nloc))!,sweight(sngi))
    allocate(SN_orig(sngi,snloc),SNLX_orig(sngi,sndim,snloc), s_cont(sngi,ndim) )
    allocate(u_loc(ndim,nloc), u_loc2(ndim,nloc), ugi(ngi,ndim), u_ele(ndim,nloc,totele))
    allocate(xsgi(sngi,ndim), usgi(sngi,ndim), usgi2(sngi,ndim), income(sngi), snorm(sngi,ndim), norm(ndim))
    allocate(tnew_loc(nloc), tnew_loc2(nloc), told(nloc,totele), tnew(nloc,totele), tnew_nonlin(nloc,totele))
    allocate(tnew_nonlin_loc(nloc))
    allocate(t_bc(nloc,totele), u_bc(ndim,nloc,totele))
    allocate(tnew_gi(ngi), tnew_xgi(ngi,ndim), tnew_sgi(sngi), tnew_sgi2(sngi), told_gi(ngi))
    allocate(face_sn(sngi,snloc,nface), face_sn2(sngi,snloc,n_s_list_no), face_snlx(sngi,sndim,nloc,nface) )
    allocate(x_loc(ndim,nloc))
    allocate(mass_ele(nloc,nloc), mat_loc_inv(nloc,nloc), rhs_loc(nloc), ml_ele(nloc))
    allocate(rhs_jac(nloc), mass_t_new(nloc), inv_jac(ngi, ndim, ndim ) )
    allocate(a_star(ngi,nloc), p_star(ngi), rgi(ngi), diff_coe(ngi), told_loc(nloc), stab(nloc, nloc))
    allocate(mat_tnew(nloc), mat_diag_approx(nloc), mass_told(nloc), stnew_loc2(snloc,nface))


    allocate(L1(ngi), L2(ngi), L3(ngi), L4(ngi))
    allocate(NLX_ALL(ngi,ndim,nloc))

    ! CPU_TIME(t1)
    ! SYSTEM_CLOCK
    ! character(10) :: time
    ! call DATE_AND_TIME(TIME = time)
    ! print *, 'Hour:  ', time(1)
    ! print *, 'Minute:', time(2)
    ! print *, 'Second:', time(3)

    t_bc(:,:)=0.0 ! this contains the boundary condif ( itime==2 .AND. ele==2 .and. its==1 ) then

    u_bc(:,:,:)=0.0 ! this contains the boundary conditions on velocity just outside the domain
    tnew(:,:) = 0.0
    face_nodes(1,1) = 1
    face_nodes(2,1) = 3
    face_nodes(1,2) = 3
    face_nodes(2,2) = 2
    face_nodes(1,3) = 2
    face_nodes(2,3) = 1


    ! 1D initial conditions
    ! tnew(:,1:2) = 1.0 ! only correct for 1D?
    ! do ele=1,no_ele_col-50
    !   tnew(:,ele*no_ele_row+10:ele*no_ele_row+50) = 1.0 ! just for now, you can use the line above for 1D or lines below for 2D instead.
    ! end do

    ! 2D initial conditions
      tnew(:,1:no_ele_row/5+1) = 1.0
    do i=2,no_ele_col
      tnew(:,(i-1)*no_ele_row+1 : no_ele_row*i-(no_ele_row*4/5)) = 1.0
    end do

    u_ele(:,:,:) = 0
    u_ele(1,:,1:totele) = u_x ! suggest this should be unity
    u_ele(2,:,1:totele) = u_y ! suggest this should be unity
    u_bc(:,:,:)=u_ele(:,:,:) ! this contains the boundary conditions on velocity just outside the domain

    dt = CFL*dx
    ! dt = CFL/(u_x/dx+u_y/dy)
    ntime = 100!time/dt


    with_time_slab =.false.
    D3=.false.
    ! IQADRA = ipoly+1
    npoly=1; ele_type= 101

    call TRIQUAold(L1, L2, L3, L4, WEIGHT, D3,NGI)

    ! call CPU_TIME(t1_get_shape_funs_spec)
    call get_shape_funs_spec(n, nlx, nlx_lxx, nlxx, weight, nlx_nod, &
              nloc, snloc, sngi, ngi, ndim, nface,n_s_list_no, face_sn, face_sn2, face_snlx, &
              sweight, npoly, ele_type, totele)!, face_list_no)


    ! call CPU_TIME(t2_get_shape_funs_spec)

    ! time_tri_ele_info2 = 0.0
    ! time_tri_det_nlx = 0.0
    ! time_det_snlx_all = 0.0
    print*, 'totele = ',totele
    print*, 'ntime = ', ntime
    call CPU_TIME(t1)
    ! call DATE_AND_TIME(TIME = start)

    do itime=1,ntime
      !generating VTK files
      if ( vtk_io <= itime  ) then
        call get_vtk_str_tri(tnew, totele, no_ele_row, no_ele_col, nface, &
                                  dx, dy, nloc, itime, ndim, cell_type)
        vtk_io = vtk_io + vtk_interval
      end if

      told = tnew ! prepare for next time step
      tnew_nonlin = tnew
      do its=1,nits
        tnew = tnew_nonlin ! for non-linearity
        do ele=1,totele
          ! call CPU_TIME(t1_tri_ele_info2)
          ! call tri_ele_info2(totele, ele, nface, face_ele2, no_ele_row, &
          !                   x_loc, dx, dy, ndim, nloc, no_ele_col)
          call str_tri_X_nodes(ele, x_loc, ndim, nloc, dx, dy, no_ele_row, no_ele_col)
          ! call CPU_TIME(t2_tri_ele_info2)
          ! time_tri_ele_info2 = time_tri_ele_info2 + (t2_tri_ele_info2 - t1_tri_ele_info2)

          ! volume integration
          ! call CPU_TIME(t1_tri_det_nlx)
          call tri_det_nlx( x_loc, n, nlx, nx, detwei, weight, ndim, nloc, ngi, inv_jac )
          ! call CPU_TIME(t2_tri_det_nlx)
          ! time_tri_det_nlx = time_tri_det_nlx + (t2_tri_det_nlx - t1_tri_det_nlx)

          !volume=sum(detwei)
      !      do gi=1,ngi
      !         print *,'gi=',gi
      !         print *,'nx(gi,1,:print):',nx(gi,1,:)
      !         print *,'nx(gi,2,:):',nx(gi,2,:)
      !      end do
      !      stop 221
          u_loc = u_ele(:,:,ele) ! this is the
          !AMIN changed to t_old
          tnew_loc =  tnew(:,ele) ! this is the FEM element - 2nd index is the local node number
          told_loc =  told(:,ele) ! this is the FEM element - 2nd index is the local node number

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
            ! rgi(gi)=(tnew_gi(gi)-told_gi(gi))/dt + sum(ugi(gi,:)*tnew_xgi(gi,:))
            ! ! a_star
            ! ! eq 4
            ! ! a_coef=sum(ugi(gi,:)*txgi(gi,:))/max(toler, sum( txgi(gi,:)**2 ) )
            ! a_coef=rgi(gi)/max(toler, sum( tnew_xgi(gi,:)**2 ) )
            ! a_star(gi,:) = a_coef * tnew_xgi(gi,:)
            ! ! eq 14
            ! ! p_star(gi) =0.0
            ! ! do iloc=1,nloc
            ! !    p_star(gi) = max(p_star(gi), abs(sum( a_star(gi,:)*nx(gi,:,iloc) ))  )
            ! ! end do
            ! ! p_star(gi) = 0.25/max(toler, p_star(gi))
            ! ! eq 23 for P*
            ! p_star(gi) =0.0
            ! do idim=1,ndim
            !   do iloc=1,nloc
            !      p_star(gi) = max(p_star(gi), abs(sum( a_star(gi,:)*(inv_jac(gi,:,idim)))  ))
            !   end do
            ! end do
            ! p_star(gi) = min(1/toler ,0.25/p_star(gi) )
            ! ! eq 18
            ! diff_coe(gi) = 0.25*rgi(gi)**2 *p_star(gi) /max(toler, sum( tnew_xgi(gi,:)**2 ) )
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

            ml_ele(iloc)=sum(n(:,iloc)*detwei(:)) ! lumped mass matrix in element (use later for solver).
            do gi=1,ngi
              do idim=1,ndim
                rhs_loc(iloc) = rhs_loc(iloc) + nx(gi,idim,iloc)*ugi(gi,idim)*tnew_gi(gi)*detwei(gi)
              end do
            end do ! quadrature
          end do ! iloc

          !Include the surface integral here:
          do iface = 1,nface
            call tri_ele_info2(totele, ele, iface, ele22, no_ele_row)

            i_got_boundary = (sign(1, -ele22) + 1 )/2 ! 1 means it is on the boundary
            r_got_boundary = real(i_got_boundary)
            ele2 = ele*i_got_boundary + ele22*(1-i_got_boundary)
            ! r_got_boundary=1.0 if we want to use the boundary conditions and have incomming velocity.
            ! r_got_boundary=0.0 not on boundary

            tnew_loc2(:)=tnew(:,ele2) * (1.0-r_got_boundary)   + t_bc(:,ele)  * r_got_boundary
            ! stnew_loc2(2,1) = tnew_loc2(1)
            ! stnew_loc2(1,1) = tnew_loc2(3)
            ! stnew_loc2(2,2) = tnew_loc2(3)
            ! stnew_loc2(1,2) = tnew_loc2(2)
            ! stnew_loc2(2,3) = tnew_loc2(2)
            ! stnew_loc2(1,3) = tnew_loc2(1)

            stnew_loc2(2,iface) = tnew_loc2(face_nodes(1,iface))
            stnew_loc2(1,iface) = tnew_loc2(face_nodes(2,iface))

            u_loc2(:,:)=u_ele(:,:,ele2)* (1.0-r_got_boundary)  + u_bc(:,:,ele)* r_got_boundary

            !Surface integral along an element
            ! s_list_no = face_list_no( iface, ele)
            sn   = face_sn(:,:,iface)
            snlx = face_snlx(:,:,:,iface)
            sn2  = face_sn2(:,:,iface)

            usgi=0.0; usgi2=0.0; xsgi=0.0; tnew_sgi=0.0; tnew_sgi2=0.0
            do iloc=1,snloc ! use all of the nodes not just the surface nodes.
              do idim=1,ndim
                usgi(:,idim)  = usgi(:,idim)  + sn(:,iloc)*u_loc(idim,face_nodes(iloc,iface))
                xsgi(:,idim)  = xsgi(:,idim)  + sn(:,iloc)*x_loc(idim,face_nodes(iloc,iface))
              end do
              tnew_sgi  = tnew_sgi(:)  + sn(:,iloc)*tnew_loc(face_nodes(iloc,iface))
            end do

            do iloc=1,snloc ! use all of the nodes not just the surface nodes.
              do idim=1,ndim
                usgi2(:,idim) = usgi2(:,idim) + sn2(:,iloc)*u_loc2(idim,iloc)
              end do
              tnew_sgi2 = tnew_sgi2(:) + sn2(:,iloc)*stnew_loc2(iloc,iface)
            end do


            ! this is the approximate normal direction...
            do idim=1,ndim
               norm(idim) = sum(xsgi(:,idim))/real(sngi) - sum(x_loc(idim,:))/real(nloc)
            end do
            ! start to do the integration
            ! call CPU_TIME(t1_det_snlx_all)
            call det_snlx_all( nloc, sngi, sndim, ndim, x_loc, sn, snlx, sweight, sdetwei, sarea, snorm, norm )
            ! call CPU_TIME(t2_det_snlx_all)
            ! time_det_snlx_all = time_det_snlx_all + (t2_det_snlx_all - t1_det_snlx_all)

                 ! print *,'iface,sum(sdetwei(:)):',iface,sum(sdetwei(:))
            ! income=1 if info is comming from neighbouring element.
            do s_gi=1,sngi
              income(s_gi)=0.5 + 0.5*sign(1.0, -sum(snorm(s_gi,:)*0.5*(usgi(s_gi,:)+usgi2(s_gi,:)))  )
            end do
                   ! print *,'iface,income:',iface,income
                   ! print *,'snorm(1,:):',snorm(1,:)
                   ! print *,'snorm(2,:):',snorm(2,:)
                   ! print *,'ele22, i_got_boundary, r_got_boundary, ele, ele2:',ele22, i_got_boundary, r_got_boundary, ele, ele2

            ! sloc_vec=0.0; sloc_vec2=0.0
            do idim=1,ndim
              s_cont(:,idim) = snorm(:,idim)*sdetwei(:) &
                          *( (1.-income(:))* usgi(:,idim)*tnew_sgi(:) + income(:)*usgi2(:,idim)*tnew_sgi2(:) )
            end do

            do iloc=1,snloc
              do idim=1,ndim
                rhs_loc(face_nodes(iloc,iface))=rhs_loc(face_nodes(iloc,iface))-sum(sn(:,iloc)*s_cont(:,idim))
              end do
            end do
          end do ! iface


          ! if (with_stab) then
            mat_loc= mass_ele + dt*stab
          ! else
          !   mat_loc= mass_ele
          ! end if

          ! if(direct_solver) then
          !   ! inverse of the mass matric (nloc,nloc)
          !   call FINDInv(mass_ele, mat_loc_inv, nloc, errorflag)
          !   !can be only this part instead
          !   ! do iloc=1,nloc
          !   !    tnew_nonlin(iloc,ele)=told_loc(iloc) + dt* sum( mat_loc_inv(iloc,:)* rhs_loc(:) ) &
          !   !                                  + told_loc(iloc)* dt* sum( mat_loc_inv(iloc,:) * stab(iloc,:) )
          !   ! end do
          !   do iloc=1,nloc
          !     mass_told(iloc)=sum( mat_loc(iloc,:) * told_loc(:) )
          !   end do
          !   do iloc=1,nloc
          !     tnew_nonlin(iloc,ele)=sum( mat_loc_inv(iloc,:)* (mass_told(:)+dt*rhs_loc(:)))
          !   end do
          !
          ! else ! iterative solver
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
          ! endif ! endof if(direct_solver) then else

        end do ! do ele=1,totele
        tnew=tnew_nonlin
      end do ! do its=1,nits
      print*, 'srtuctured ',itime
    end do ! do itime=1,ntime

    !runing time
    call CPU_TIME(t2)
    print*, 'cpu_time for time_loop = ', t2-t1
    ! ! call DATE_AND_TIME(TIME = finish)
    ! ! print*, 'finishing time = ', finish, '  starting time =  ',start
    ! !
    ! ! print*, 'cpu_time get_shape_funs_spec = ', t2_get_shape_funs_spec - t1_get_shape_funs_spec
    ! print*, 'cpu_time tri_ele_info2 = ', time_tri_ele_info2
    ! ! print*, ''
    ! !
    ! print*, 'cpu_time tri_det_nlx = ', time_tri_det_nlx
    ! print*, ''
    ! !
    ! print*, 'cpu_time det_snlx_all = ', time_det_snlx_all
    ! print*, ''


    ! ! OPEN(unit=10, file='1DG-FEM, test_triangle')
    ! !   do ele=1,totele
    ! !     write(10,*) x_all(1,1,ele), x_all(2,1,ele), tnew(1,ele)
    ! !     write(10,*) x_all(1,2,ele), x_all(2,2,ele), tnew(2,ele)
    ! !   end do
    ! ! close(10)
    ! ! local face numbers:
    ! !    |\        __1__
    ! !   2| \ 3     \   |
    ! !    |  \      3\  | 2
    ! !    |___\       \ |
    ! !      1          \|
    !
    ! ! local node numbers
    ! !    2
    ! !    |\        1___3
    ! !    | \       \   |
    ! !    |  \       \  |
    ! !    |___\       \ |
    ! !    3   1        \|
    ! !                 2
    ! OPEN(unit=10, file='DG_triangle')
    !   do ele=1,totele
    !     if ( MOD(ele,2).ne.0 ) then
    !       write(10,*) x_all(1,3,ele), x_all(2,3,ele), tnew(3,ele)
    !       write(10,*) x_all(1,1,ele), x_all(2,1,ele), tnew(1,ele)
    !     else
    !       write(10,*) x_all(1,2,ele), x_all(2,2,ele), tnew(2,ele)
    !       ! write(10,*) x_all(1,2,ele), x_all(2,2,ele), tnew(2,ele)
    !     end if
    !   end do
    !   ! do ele=1,totele
    !   !   if ( MOD(ele,2).ne.0 ) then
    !   !     write(10,*) x_all(1,2,ele), x_all(2,2,ele), tnew(2,ele)
    !   !   else
    !   !     write(10,*) x_all(1,1,ele), x_all(2,1,ele), tnew(1,ele)
    !   !     write(10,*) x_all(1,3,ele), x_all(2,3,ele), tnew(3,ele)
    !   !   end if
    !   ! end do
    ! close(10)
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
    deallocate(face_ele2, n, nlx, nx, M,&
               nlx_lxx, nlxx, nlx_nod, weight, detwei, sdetwei, sn,sn2,snlx,sweight, s_cont,&
               tnew_loc, tnew_loc2, told, tnew, tnew_nonlin, tnew_nonlin_loc, t_bc, u_bc,&
               tnew_gi, tnew_xgi, tnew_sgi, tnew_sgi2, told_gi, face_sn, face_sn2, face_snlx,&
               u_loc, u_ele, u_loc2, x_loc, ugi, xsgi, usgi, usgi2, income, snorm, norm,&
               mass_ele, mat_loc, mat_loc_inv, rhs_loc, ml_ele, rhs_jac, mass_t_new,&
               SN_orig,SNLX_orig, inv_jac, mat_diag_approx,&
               a_star, p_star, rgi, diff_coe, told_loc, stab, mat_tnew, mass_told,&
               L1, L2, L3, L4, NLX_ALL,stnew_loc2)
  end subroutine str_explicit
end module transport_tri
