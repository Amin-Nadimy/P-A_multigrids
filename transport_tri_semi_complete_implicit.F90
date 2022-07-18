module transport_tri_semi
  use ShapFun
  use structured_meshgen
  use ShapFun_unstruc
  use matrices
  use Msh2Tri
  use splitting
  use get_vtk_files
  contains


    ! brief this is for semi-str implicit method
    Subroutine Semi_implicit_iterative(CFL, ndim, direct_solver, njac_its, time, nits,&
      nface, dx, dy, nloc, snloc, ngi, sngi, n_s_list_no, u_x, u_y, with_stab, vtk_interval, cell_type)
      implicit none

      ! external variables
      real, intent(in) :: CFL, time, u_x, u_y, dx, dy
      integer, intent(in) :: ndim, nits, nface, nloc, ngi, sngi, snloc, njac_its
      integer, intent(in) :: vtk_interval, cell_type, n_s_list_no
      logical, intent(in) :: direct_solver, with_stab

      ! local vbls
      type(Mesh), allocatable, dimension(:) :: meshList
      type(sparse) :: sparse_mass, sparse_lhs, sparse_flux

      character (len=20):: file_name
      logical :: with_time_slab, D3

      integer :: gi, iface, ele, ele2, ele22, n_split, irow, ipos, sum_up, irow_ups
      integer :: totele_unst, un_ele, str_ele, totele_str
      integer :: s_list_no, s_gi, iloc, jloc, i,j,k, ntime, vtk_io
      integer :: itime, idim, sndim, nonodes, mloc, col, row, row2
      integer :: errorflag, i_got_boundary, jac_its, its
      integer :: npoly, ele_type, ierr, i_overlaps, target_loc
      integer :: Mpos, Npos, Mside, Nside, lnodes(2), lnodes2(2), totnodes, mface, sp
      integer :: g_iloc, g_jloc, target_ele, orientation

      real :: sarea, volume, dt, L, r_got_boundary, toler, a_coef
      real :: length_x, length_y
      real :: t1, t2, t1_get_shape_funs_spec, t2_get_shape_funs_spec
      real :: t1_tri_det_nlx, t2_tri_det_nlx, time_tri_det_nlx, i_tri_det_nlx
      real :: t1_det_snlx_all, t2_det_snlx_all, time_det_snlx_all, i_det_snlx_all
      real :: t1_ReadMSH, t2_ReadMSH, time_ReadMSH
      real :: t1_getNeigDataMesh, t2_getNeigDataMesh, time_getNeigDataMesh
      real :: t1_get_unstr_sn2, t2_get_unstr_sn2, time_get_unstr_sn2, str_x(2,3)
      real :: mass_ele, stiff_ele, flux_ele2, lhs, flux_ele

      integer, allocatable ::  face_ele(:,:), face_list_no(:,:), str_neig(:,:), res(:), result(:)

      real, allocatable :: n(:,:), nlx(:,:,:), nx(:,:,:), M(:,:)
      real, allocatable :: nlx_lxx(:,:,:), nlxx(:,:), nlx_nod(:,:,:)
      real, allocatable :: weight(:), detwei(:), sdetwei(:)
      real, allocatable :: sn(:,:),sn2(:,:),snlx(:,:,:), sweight(:), s_cont(:,:), sweight_str(:)
      real, allocatable :: tnew_loc(:), tnew_loc2(:), told(:,:,:), tnew(:,:,:), tnew_nonlin(:,:,:)
      real, allocatable :: tnew_nonlin_loc(:)
      real, allocatable :: t_bc(:,:), u_bc(:,:,:)
      real, allocatable :: tnew_gi(:), tnew_xgi(:,:), tnew_sgi(:), tnew_sgi2(:), told_gi(:)
      real, allocatable :: face_sn(:,:,:), face_sn2(:,:,:), face_snlx(:,:,:,:)
      real, allocatable :: face_sn_str(:,:,:), face_sn2_str(:,:,:), face_snlx_str(:,:,:,:)
      real, allocatable :: u_loc(:,:), u_ele(:,:,:), u_loc2(:,:), x_loc(:,:), ugi(:,:)
      real, allocatable :: xsgi(:,:), usgi(:,:), usgi2(:,:), income(:), snorm(:,:), norm(:)
      real, allocatable :: mat_loc(:,:), mat_loc_inv(:,:), rhs_loc(:), ml_ele(:), rhs_jac(:), mass_t_new(:)
      real, allocatable :: SN_orig(:,:),SNLX_orig(:,:,:), inv_jac(:, :, : ), mat_diag_approx(:)
      real, allocatable :: a_star(:,:), p_star(:), rgi(:), diff_coe(:), told_loc(:), stab(:,:), mat_tnew(:), mass_told(:)
      real, allocatable :: L1(:), L2(:), L3(:), L4(:), x_all_str(:,:,:), glob_told_array(:)
      real, allocatable :: glob_lhs(:,:), inv_lhs(:,:), rhs_array(:)
      real, allocatable :: tnew_nonlin2(:)

      ! xp(2): coordinates of nodes
      ! call ReadMSH(meshList,'./semi_mesh.msh',ierr, totnodes)
      ! call ReadMSH(meshList,'./irregular shape.msh',ierr, totnodes)
      ! call ReadMSH(meshList,'./new_semi_mesh.msh',ierr, totnodes)
      ! call ReadMSH(meshList,'./1_unele_test.msh',ierr, totnodes)
      ! call ReadMSH(meshList,'./gmsh_100.msh',ierr, totnodes)
      print*, 'starting reading the .msh file'
      call CPU_TIME(t1_ReadMSH)
      ! call ReadMSH(meshList,'./10_elements.msh',ierr, totnodes)
      call ReadMSH(meshList,'./1_unele_test.msh',ierr, totnodes)
      ! call ReadMSH(meshList,'./2_unele_test.msh',ierr, totnodes)
      ! call ReadMSH(meshList,'./4_elements.msh',ierr, totnodes)
      ! call ReadMSH(meshList,'./irregular.msh',ierr, totnodes)
      call CPU_TIME(t2_ReadMSH)
      time_ReadMSH = t2_ReadMSH - t1_ReadMSH
      print*, 'finished reading .msh file'

      totele_unst = size(meshList)

      n_split =1
      totele_str = 4**n_split

      call get_str_neig(n_split, str_neig)

      sndim = ndim-1 ! dimensions of the surface element coordinates.
      mloc=1
      toler = 0.00000000001
      ele_type = 3
      vtk_io=vtk_interval

      allocate(n( ngi, nloc ), nx( ngi, ndim, nloc ), nlx( ngi, ndim, nloc ), M(ngi,nloc))
      allocate(nlx_lxx(ngi,ndim,nloc), nlxx(ngi,nloc), nlx_nod(nloc,ndim,nloc))
      allocate(weight(ngi), detwei(ngi), sdetwei(sngi))
      allocate(face_list_no( nface, totele_str), face_ele(nface,totele_str))
      allocate(sn(sngi,nloc),sn2(sngi,nloc),snlx(sngi,sndim,nloc))!,sweight(sngi))
      allocate(SN_orig(sngi,snloc),SNLX_orig(sngi,sndim,snloc), s_cont(sngi,ndim) )
      allocate(u_loc(ndim,nloc), u_loc2(ndim,nloc), ugi(ngi,ndim))
      allocate(xsgi(sngi,ndim), usgi(sngi,ndim), usgi2(sngi,ndim), income(sngi), snorm(sngi,ndim), norm(ndim))
      allocate(tnew_loc(nloc), tnew_loc2(nloc), told(nloc,totele_str,totele_unst))
      allocate(tnew(nloc,totele_str,totele_unst), tnew_nonlin(nloc,totele_str,totele_unst))
      allocate(tnew_nonlin_loc(nloc))
      allocate(t_bc(2**n_split*nloc,4), u_bc(ndim,2**n_split*nloc,4)) ! it works just for semi_mesh.msh
      allocate(tnew_gi(ngi), tnew_xgi(ngi,ndim), tnew_sgi(sngi), tnew_sgi2(sngi), told_gi(ngi))
      allocate(face_sn(sngi,nloc,nface), face_sn2(sngi,nloc,n_s_list_no), face_snlx(sngi,sndim,nloc,nface) )
      allocate(face_sn_str(sngi,nloc,nface), face_sn2_str(sngi,nloc,n_s_list_no), face_snlx_str(sngi,sndim,nloc,nface))
      allocate(x_loc(ndim,nloc))
      allocate(mat_loc_inv(nloc,nloc), rhs_loc(nloc), ml_ele(nloc))
      allocate(rhs_jac(nloc), mass_t_new(nloc), inv_jac(ngi, ndim, ndim ) )
      allocate(a_star(ngi,nloc), p_star(ngi), rgi(ngi), diff_coe(ngi), told_loc(nloc), stab(nloc, nloc))
      allocate(mat_tnew(nloc), mat_diag_approx(nloc), mass_told(nloc))
      allocate(L1(ngi), L2(ngi), L3(ngi), L4(ngi), x_all_str(ndim,nloc,totele_unst*totele_str))
      allocate(inv_lhs(totele_str*totele_unst*nloc, totele_str*totele_unst*nloc))
      allocate(rhs_array(totele_str*totele_unst*nloc), tnew_nonlin2(totele_str*totele_unst*nloc))
      allocate(glob_lhs(totele_str*totele_unst*nloc,totele_str*totele_unst*nloc))

      do un_ele=1,size(meshlist)
        allocate(meshList(un_ele)%t_overlap(2**n_split*nloc, nface))
        allocate(meshList(un_ele)%t_overlap_old(2**n_split*nloc, nface))
        allocate(meshList(un_ele)%u_overlap(ndim, 2**n_split*nloc , nface))
        allocate(meshList(un_ele)%u_ele(ndim, nloc, 4**n_split))
        allocate(meshList(un_ele)%S_nodes(sngi,nface))
        allocate(meshList(un_ele)%s_ele(2**n_split,nface))
        meshList(un_ele)%t_overlap=0.0
        meshList(un_ele)%u_overlap(1,:,:)=u_x
        meshList(un_ele)%u_overlap(2,:,:)=u_y
        meshList(un_ele)%u_ele(:,:,:) = 0
        meshList(un_ele)%u_ele(1,:,:) = u_x ! suggest this should be unity
        meshList(un_ele)%u_ele(2,:,:) = u_y ! suggest this should be unity
        do iface=1,nface

          call getNeigDataMesh(meshlist, un_ele, iface, Npos, Nside, lnodes2)
          meshList(un_ele)%S_nodes(:,iface)=lnodes2
          meshList(un_ele)%fNeig(iface)=Nside
        end do
      end do

      ! meshList(1)%t_overlap=0.0
      ! meshList(2)%t_overlap=0.0

      t_bc(:,:)=0.0 ! this contains the boundary conditions just outside the domain
      u_bc(1,:,:)=u_x ! this contains the boundary conditions on velocity just outside the domain
      u_bc(2,:,:)=u_y
      tnew(:,:,:) = 0.0

      !!!!!!!!!!!!!!!!!!!!!!!!!!! new semi mesh ICs !!!!!!!!!!!!!!!!!!!!!!!!!!!!
      tnew(:,3,1)=1.0
      do un_ele=1,totele_unst
          if ( meshList(un_ele)%region_id == 13 ) then
            tnew(:,:,un_ele) = 1.
          end if
        end do

      dt = CFL*dx
      ! dt = CFL*dx/u_x + CFL*dy/u_y
      ntime = time/dt
      with_time_slab =.false.
      D3=.false.
      npoly=1; ele_type= 101
      call TRIQUAold(L1, L2, L3, L4, WEIGHT, D3,NGI)
      call get_surface_ele(meshlist, n_split, totele_unst, nface)
      call make_sparse_matrix(sparse_mass, totele_str*totele_unst, nloc)
      call make_sparse_matrix_flux_semi(sparse_lhs,meshList,n_split,totele_unst,totele_str,nloc,nface,str_neig)
      call make_sparse_matrix_flux_semi(sparse_flux,meshList,n_split,totele_unst,totele_str,nloc,nface,str_neig)
      call get_shape_funs_spec_ustr(n, nlx, nlx_lxx, nlxx, sn_orig, weight, nlx_nod, &
                  nloc, snloc, sngi, ngi, ndim, nface,n_s_list_no, face_sn, face_sn2, face_snlx, &
                  sweight, npoly, ele_type, totele_unst, face_list_no )
! i=1;j=9
! do gi=1,size(sparse_flux%g_jloc)/9
!   call remove_dups(sparse_flux%g_jloc(i:j), result)
!   sparse_flux%g_jloc(i:j)= result
!   print*, sparse_flux%g_jloc(i:j)
!   i=j+1
!   j=j+9
! end do

      call get_shape_funs_spec(n, nlx, nlx_lxx, nlxx, weight, nlx_nod, &
                nloc, snloc, sngi, ngi, ndim, nface,n_s_list_no, face_sn_str, face_sn2_str, face_snlx_str, &
                sweight_str, npoly, ele_type, totele_str)!, face_list_no)

      i=1 ! saving x_all structured triangles to be used for VTK generation
      do un_ele=1,totele_unst
        do str_ele=1,totele_str
          call get_splitting(meshList(un_ele)%X, n_split, str_ele, str_x)
          x_all_str(:,:,i) = str_x !(ndim,nloc,totele)
          i=i+1
        end do
      end do
      print*, 'totele = ', totele_str * totele_unst
      print*, 'ntime = ', ntime
      call CPU_TIME(t1)


      do itime=1,ntime
        ! generating VTK files
        if ( vtk_io <= itime ) then
          call get_vtk(x_all_str, tnew, totele_str*totele_unst, nloc, itime,ndim, cell_type)
          vtk_io = vtk_io + vtk_interval
        end if

        told = tnew ! prepare for next time step
        tnew_nonlin = tnew
        do its=1,nits
          tnew = tnew_nonlin ! for non-linearity
          call update_overlaps(meshlist, tnew, told, t_bc, n_split, nface, totele_unst, nloc, str_neig)
          sparse_mass%val=0.
          sparse_lhs%val=0.
          sparse_flux%val=0.

          do un_ele=1,totele_unst
            do str_ele=1,totele_str
              call get_str_info(n_split, str_ele, irow, ipos, orientation)
              call get_splitting(meshList(un_ele)%X, n_split, str_ele, x_loc)
              call tri_det_nlx( x_loc, n, nlx, nx, detwei, weight, ndim, nloc, ngi, inv_jac )

              !volume=sum(detwei)
              u_loc = meshList(un_ele)%u_ele(:,:,str_ele) ! this is the
              tnew_loc =  tnew(:,str_ele, un_ele) ! this is the FEM element - 2nd index is the local node number
              told_loc =  told(:,str_ele, un_ele) ! this is the FEM element - 2nd index is the local node number

              do gi=1,ngi
                do idim=1,ndim
                  ugi(gi,idim)=sum(n(gi,:)*u_loc(idim,:))
                  tnew_xgi(gi,idim)=sum(nx(gi,idim,:)*tnew_loc(:))
                end do
                tnew_gi(gi)=sum(n(gi,:)*tnew_loc(:))
                told_gi(gi)=sum(n(gi,:)*told_loc(:))
                rgi(gi)=(tnew_gi(gi)-told_gi(gi))/dt + sum(ugi(gi,:)*tnew_xgi(gi,:))
                ! a_starif(ele_type < is_triangle_or_tet) then
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
                g_iloc = glob_no_semi(un_ele, str_ele, n_split, nloc, iloc)
                do jloc=1,nloc
                  g_jloc = glob_no_semi(un_ele, str_ele, n_split, nloc, jloc)
                  mass_ele = 0.0
                  stiff_ele = 0.0
                  do gi=1,ngi
                    stab(iloc,jloc) = stab(iloc,jloc) + sum(diff_coe(gi)* nx(gi,:,jloc)* nx(gi,:,iloc))* detwei(gi)
                    mass_ele = mass_ele + n(gi,iloc)*n(gi,jloc)*detwei(gi)/dt
                    do idim=1,ndim
                      stiff_ele = stiff_ele + nx(gi,idim,iloc)*ugi(gi,idim)*detwei(gi)*n(gi,jloc)
                    end do
                  end do ! quadrature
                  call add_to_CSR(sparse_mass, g_iloc, g_jloc, mass_ele)
                  lhs = mass_ele-stiff_ele
                  call add_to_CSR_flux(sparse_lhs, g_iloc, g_jloc, lhs)
                end do ! jloc
              end do ! iloc



              do iface = 1,nface
                ele22 = str_neig(iface,str_ele)!face_ele( iface, ele) ! if ele22 is 0 assume on boundary
                i_got_boundary = (sign(1, -ele22) + 1 )/2
                r_got_boundary = real(i_got_boundary)
                ele2 = str_ele*i_got_boundary + ele22*(1-i_got_boundary)! if ele22=0, then ele2=ele, if ele22/=0, ele2=ele22

                sn   = face_sn_str(:,:,iface)
                snlx = face_snlx_str(:,:,:,iface)

                sn2=0
                if ( ele22==0 .and. iface==1 ) then
                  sp = ipos/2+1 ! position of ele along the un_iface
                  mface = 1 ! un_iface number which str_ele is located on
                  call get_semi_sn2_implicit(meshlist, un_ele, mface, sn_orig, sn2)

                else if ( ele22==0 .and. iface==2 ) then
                  sp = irow
                  mface = 3
                  call get_semi_sn2_implicit(meshlist, un_ele, mface, sn_orig, sn2)

                else if ( ele22==0 .and. iface==3 ) then
                  sp = irow
                  mface = 2
                  call get_semi_sn2_implicit(meshlist, un_ele, mface, sn_orig, sn2)

                else if ( ele22 /= 0 ) then
                  sn2  = face_sn2_str(:,:,iface)
                end if ! end if ele22=0

                tnew_loc2(:)= tnew(:,ele2, un_ele) * (1.0-r_got_boundary) &
                            + meshlist(un_ele)%t_overlap( (sp)*nloc-2:(sp)*nloc, mface ) * r_got_boundary

                u_loc2(:,:)= meshList(un_ele)%u_ele(:,:,ele2)* (1.0-r_got_boundary) &
                            + meshlist(un_ele)%u_overlap(:, sp*nloc-2:sp*nloc, mface ) * r_got_boundary
                            ! 2nd section of u_overlap should refer ro 3 ilocs from the un_ele u_bc which can be identified based on irow

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

                ! this is the approximate normal direction...
                do idim=1,ndim
                   norm(idim) = sum(xsgi(:,idim))/real(sngi) - sum(x_loc(idim,:))/real(nloc)
                end do

                call det_snlx_all( nloc, sngi, sndim, ndim, x_loc, sn, snlx, sweight, sdetwei, sarea, snorm, norm )
                ! income=1 if info is comming from neighbouring element.

                do s_gi=1,sngi
                  income(s_gi)=0.5 + 0.5*sign(1.0, -sum(snorm(s_gi,:)*0.5*(usgi(s_gi,:)+usgi2(s_gi,:)))  )
                end do

                do idim=1,ndim
                  s_cont(:,idim) = snorm(:,idim)*sdetwei(:) &
                              *( (1.-income(:))* usgi(:,idim)*tnew_sgi(:) + income(:)*usgi2(:,idim) )
                end do

                target_ele = int( (1.-income(1))*str_ele + income(1)*ele22 )
                Npos = meshList(un_ele)%Neig(mface)
                Nside = meshList(un_ele)%fNeig(mface)
                if ( target_ele==0 .and. Npos==0 ) then ! it happens when ele is on the boundary
                elseif ( target_ele==0 .and. Npos/=0) then
                  if ( iface==1 ) then
                    ipos = int(ipos/2)+1
                    target_ele = meshList(un_ele)%s_ele(ipos, 1)
                  elseif ( iface==2 ) then
                    target_ele = meshList(un_ele)%s_ele(irow, 3)
                  elseif ( iface==3 ) then
                    target_ele = meshList(un_ele)%s_ele(irow, 2)
                  end if
                elseif ( target_ele /= 0) then
                  target_ele = (un_ele-1)*4**n_split + target_ele
                end if

                do iloc=1,nloc
                  g_iloc = glob_no_semi(un_ele, str_ele, n_split, nloc, iloc)
                  do jloc=1,nloc
                    g_jloc = (target_ele-1)*nloc + jloc !glob_no_semi(un_ele, target_ele, n_split, nloc, jloc)
                    flux_ele = 0.0
                    do idim=1,ndim
                      ! glob_flux(g_iloc,g_jloc)=glob_flux(g_iloc,g_jloc) +sum( snorm(:,idim)*sdetwei(:)*sn(:,iloc) &
                      !                                                   *((1.-income(:))*sn(:,jloc)*usgi(:,idim) &
                      !                                                   +     income(:)*sn2(:,jloc)*usgi2(:,idim)))
                      flux_ele = flux_ele +sum( snorm(:,idim)*sdetwei(:)*sn(:,iloc) &
                                          *((1.-income(:))*sn(:,jloc)*usgi(:,idim) &
                                          +     income(:)*sn2(:,jloc)*usgi2(:,idim)))
                    end do ! idim
                    call add_to_CSR_flux(sparse_flux, g_iloc, g_jloc,  flux_ele)
                  end do ! jloc
                end do ! iloc
              end do ! iface

              ! if (with_stab) then
                ! mat_loc= mass_ele + dt*stab
              ! else
              !   stab=0
              !   mat_loc= mass_ele
              ! end if

            end do ! end do totele_str
          end do ! end totele_unst
          sparse_lhs%val = sparse_lhs%val + sparse_flux%val
          call gen_global_matrix(totele_str*totele_unst, nloc, sparse_lhs, glob_lhs)
          call told_to_array_semi(totele_str, totele_unst, nloc, told, glob_told_array)
          call csr_mul_array(sparse_mass, glob_told_array, rhs_array)
          call FINDInv(glob_lhs, inv_lhs, nloc*totele_str*totele_unst, errorflag)
           tnew_nonlin2 = matmul(inv_lhs,rhs_array )
          call array_to_told_semi(totele_str, totele_unst, nloc, tnew_nonlin2, tnew_nonlin)
          tnew=tnew_nonlin
        end do ! do its=1,nits
      print*, 'semi', itime
      end do ! do itime=1,ntime
      call CPU_TIME(t2)
      print*, 'cpu_time for time_loop = ', t2-t1

      deallocate(face_ele, face_list_no, str_neig, n, nlx, nx, M, nlx_lxx, nlxx, nlx_nod, weight, detwei, sdetwei,&
                 sn,sn2,snlx, sweight, s_cont, sweight_str, tnew_loc, tnew_loc2, told, tnew, tnew_nonlin,&
                 tnew_nonlin_loc, t_bc, u_bc, tnew_gi, tnew_xgi, tnew_sgi, tnew_sgi2, told_gi,&
                 face_sn, face_sn2, face_snlx, face_sn_str, face_sn2_str, face_snlx_str, u_loc, u_loc2, x_loc, ugi,&
                 xsgi, usgi, usgi2, income, snorm, norm, mat_loc_inv, rhs_loc, ml_ele, rhs_jac, mass_t_new,&
                 SN_orig,SNLX_orig, inv_jac, mat_diag_approx, a_star, p_star, rgi, diff_coe, told_loc, stab, mat_tnew, mass_told,&
                 L1, L2, L3, L4, x_all_str, inv_lhs, rhs_array,glob_told_array,glob_lhs,&
                 tnew_nonlin2, sparse_flux%g_iloc, sparse_flux%g_jloc, sparse_flux%val,&
                 sparse_lhs%g_iloc, sparse_lhs%g_jloc, sparse_lhs%val, sparse_mass%g_iloc, sparse_mass%g_jloc,sparse_mass%val)
    end subroutine Semi_implicit_iterative



  ! brief this is for semi-str complete implicit method
  Subroutine Semi_implicit_direct(CFL, ndim, direct_solver, njac_its, time, nits,&
    nface, dx, dy, nloc, snloc, ngi, sngi, n_s_list_no, u_x, u_y, with_stab, vtk_interval, cell_type)
    implicit none

    ! external variables
    real, intent(in) :: CFL, time, u_x, u_y, dx, dy
    integer, intent(in) :: ndim, nits, nface, nloc, ngi, sngi, snloc, njac_its
    integer, intent(in) :: vtk_interval, cell_type, n_s_list_no
    logical, intent(in) :: direct_solver, with_stab

    ! local vbls
    type(Mesh), allocatable, dimension(:) :: meshList
    type(sparse) :: sparse_mass, sparse_lhs, sparse_flux

    character (len=20):: file_name
    logical :: with_time_slab, D3

    integer :: gi, iface, ele, ele2, ele22, n_split, irow, ipos, sum_up, irow_ups
    integer :: totele_unst, un_ele, str_ele, totele_str
    integer :: s_list_no, s_gi, iloc, jloc, i,j,k, ntime, vtk_io
    integer :: itime, idim, sndim, nonodes, mloc, col!, row, row2
    integer :: errorflag, i_got_boundary, jac_its, its
    integer :: npoly, ele_type, ierr, i_overlaps, target_loc
    integer :: Mpos, Npos, Mside, Nside, lnodes(2), lnodes2(2), totnodes, mface, sp
    integer :: g_iloc, g_jloc, target_ele, orientation

    real :: sarea, volume, dt, L, r_got_boundary, toler, a_coef
    real :: length_x, length_y
    real :: t1, t2, t1_get_shape_funs_spec, t2_get_shape_funs_spec
    real :: t1_tri_det_nlx, t2_tri_det_nlx, time_tri_det_nlx, i_tri_det_nlx
    real :: t1_det_snlx_all, t2_det_snlx_all, time_det_snlx_all, i_det_snlx_all
    real :: t1_ReadMSH, t2_ReadMSH, time_ReadMSH
    real :: t1_getNeigDataMesh, t2_getNeigDataMesh, time_getNeigDataMesh
    real :: t1_get_unstr_sn2, t2_get_unstr_sn2, time_get_unstr_sn2, str_x(2,3)
    real :: mass_ele, stiff_ele, flux_ele2, lhs, flux_ele

    integer, allocatable ::  face_ele(:,:), face_list_no(:,:), str_neig(:,:)

    real, allocatable :: n(:,:), nlx(:,:,:), nx(:,:,:), M(:,:)
    real, allocatable :: nlx_lxx(:,:,:), nlxx(:,:), nlx_nod(:,:,:)
    real, allocatable :: weight(:), detwei(:), sdetwei(:)
    real, allocatable :: sn(:,:),sn2(:,:),snlx(:,:,:), sweight(:), s_cont(:,:), sweight_str(:)
    real, allocatable :: tnew_loc(:), tnew_loc2(:), told(:,:,:), tnew(:,:,:), tnew_nonlin(:,:,:)
    real, allocatable :: tnew_nonlin_loc(:)
    real, allocatable :: t_bc(:,:), u_bc(:,:,:)
    real, allocatable :: tnew_gi(:), tnew_xgi(:,:), tnew_sgi(:), tnew_sgi2(:), told_gi(:)
    real, allocatable :: face_sn(:,:,:), face_sn2(:,:,:), face_snlx(:,:,:,:)
    real, allocatable :: face_sn_str(:,:,:), face_sn2_str(:,:,:), face_snlx_str(:,:,:,:)
    real, allocatable :: u_loc(:,:), u_ele(:,:,:), u_loc2(:,:), x_loc(:,:), ugi(:,:)
    real, allocatable :: xsgi(:,:), usgi(:,:), usgi2(:,:), income(:), snorm(:,:), norm(:)
    real, allocatable :: mat_loc(:,:), mat_loc_inv(:,:), rhs_loc(:), ml_ele(:), rhs_jac(:), mass_t_new(:)
    real, allocatable :: SN_orig(:,:),SNLX_orig(:,:,:), inv_jac(:, :, : ), mat_diag_approx(:)
    real, allocatable :: a_star(:,:), p_star(:), rgi(:), diff_coe(:), told_loc(:), stab(:,:), mat_tnew(:), mass_told(:)
    real, allocatable :: L1(:), L2(:), L3(:), L4(:), x_all_str(:,:,:), glob_flux(:,:), glob_told_array(:)
    real, allocatable :: glob_lhs(:,:), inv_lhs(:,:), rhs_array(:), glob_mass(:,:)
    real, allocatable :: tnew_nonlin2(:)

    ! xp(2): coordinates of nodes
    ! call ReadMSH(meshList,'./semi_mesh.msh',ierr, totnodes)
    ! call ReadMSH(meshList,'./irregular shape.msh',ierr, totnodes)
    ! call ReadMSH(meshList,'./new_semi_mesh.msh',ierr, totnodes)
    ! call ReadMSH(meshList,'./1_unele_test.msh',ierr, totnodes)
    ! call ReadMSH(meshList,'./gmsh_100.msh',ierr, totnodes)
    print*, 'starting reading the .msh file'
    call CPU_TIME(t1_ReadMSH)
    ! call ReadMSH(meshList,'./10_elements.msh',ierr, totnodes)
    ! call ReadMSH(meshList,'./1_unele_test.msh',ierr, totnodes)
    ! call ReadMSH(meshList,'./2_unele_test.msh',ierr, totnodes)
    ! call ReadMSH(meshList,'./4_elements.msh',ierr, totnodes)
    call ReadMSH(meshList,'./irregular.msh',ierr, totnodes)
    call CPU_TIME(t2_ReadMSH)
    time_ReadMSH = t2_ReadMSH - t1_ReadMSH
    print*, 'finished reading .msh file'

    totele_unst = size(meshList)

    n_split =2
    totele_str = 4**n_split

    call get_str_neig(n_split, str_neig)

    sndim = ndim-1 ! dimensions of the surface element coordinates.
    mloc=1
    toler = 0.00000000001
    ele_type = 3
    vtk_io=vtk_interval

    allocate(n( ngi, nloc ), nx( ngi, ndim, nloc ), nlx( ngi, ndim, nloc ), M(ngi,nloc))
    allocate(nlx_lxx(ngi,ndim,nloc), nlxx(ngi,nloc), nlx_nod(nloc,ndim,nloc))
    allocate(weight(ngi), detwei(ngi), sdetwei(sngi))
    allocate(face_list_no( nface, totele_str), face_ele(nface,totele_str))
    allocate(sn(sngi,nloc),sn2(sngi,nloc),snlx(sngi,sndim,nloc))!,sweight(sngi))
    allocate(SN_orig(sngi,snloc),SNLX_orig(sngi,sndim,snloc), s_cont(sngi,ndim) )
    allocate(u_loc(ndim,nloc), u_loc2(ndim,nloc), ugi(ngi,ndim))
    allocate(xsgi(sngi,ndim), usgi(sngi,ndim), usgi2(sngi,ndim), income(sngi), snorm(sngi,ndim), norm(ndim))
    allocate(tnew_loc(nloc), tnew_loc2(nloc), told(nloc,totele_str,totele_unst))
    allocate(tnew(nloc,totele_str,totele_unst), tnew_nonlin(nloc,totele_str,totele_unst))
    allocate(tnew_nonlin_loc(nloc))
    allocate(t_bc(2**n_split*nloc,4), u_bc(ndim,2**n_split*nloc,4)) ! it works just for semi_mesh.msh
    allocate(tnew_gi(ngi), tnew_xgi(ngi,ndim), tnew_sgi(sngi), tnew_sgi2(sngi), told_gi(ngi))
    allocate(face_sn(sngi,nloc,nface), face_sn2(sngi,nloc,n_s_list_no), face_snlx(sngi,sndim,nloc,nface) )
    allocate(face_sn_str(sngi,nloc,nface), face_sn2_str(sngi,nloc,n_s_list_no), face_snlx_str(sngi,sndim,nloc,nface))
    allocate(x_loc(ndim,nloc))
    allocate(mat_loc_inv(nloc,nloc), rhs_loc(nloc), ml_ele(nloc))
    allocate(rhs_jac(nloc), mass_t_new(nloc), inv_jac(ngi, ndim, ndim ) )
    allocate(a_star(ngi,nloc), p_star(ngi), rgi(ngi), diff_coe(ngi), told_loc(nloc), stab(nloc, nloc))
    allocate(mat_tnew(nloc), mat_diag_approx(nloc), mass_told(nloc))
    allocate(L1(ngi), L2(ngi), L3(ngi), L4(ngi), x_all_str(ndim,nloc,totele_unst*totele_str))
    allocate(glob_flux(totele_str*totele_unst*nloc, totele_str*totele_unst*nloc))
    allocate(inv_lhs(totele_str*totele_unst*nloc, totele_str*totele_unst*nloc))
    allocate(rhs_array(totele_str*totele_unst*nloc), tnew_nonlin2(totele_str*totele_unst*nloc))
    allocate( glob_mass(totele_str*totele_unst*nloc,totele_str*totele_unst*nloc))
    allocate(glob_lhs(totele_str*totele_unst*nloc,totele_str*totele_unst*nloc))

    do un_ele=1,size(meshlist)
      allocate(meshList(un_ele)%t_overlap(2**n_split*nloc, nface))
      allocate(meshList(un_ele)%t_overlap_old(2**n_split*nloc, nface))
      allocate(meshList(un_ele)%u_overlap(ndim, 2**n_split*nloc , nface))
      allocate(meshList(un_ele)%u_ele(ndim, nloc, 4**n_split))
      allocate(meshList(un_ele)%S_nodes(sngi,nface))
      allocate(meshList(un_ele)%s_ele(2**n_split,nface))
      meshList(un_ele)%t_overlap=0.0
      meshList(un_ele)%u_overlap(1,:,:)=u_x
      meshList(un_ele)%u_overlap(2,:,:)=u_y
      meshList(un_ele)%u_ele(:,:,:) = 0
      meshList(un_ele)%u_ele(1,:,:) = u_x ! suggest this should be unity
      meshList(un_ele)%u_ele(2,:,:) = u_y ! suggest this should be unity
      do iface=1,nface

        call getNeigDataMesh(meshlist, un_ele, iface, Npos, Nside, lnodes2)
        meshList(un_ele)%S_nodes(:,iface)=lnodes2
        meshList(un_ele)%fNeig(iface)=Nside
      end do
    end do
    ! meshList(1)%t_overlap=0.0
    ! meshList(2)%t_overlap=0.0
    t_bc(:,:)=0.0 ! this contains the boundary conditions just outside the domain
    u_bc(1,:,:)=u_x ! this contains the boundary conditions on velocity just outside the domain
    u_bc(2,:,:)=u_y
    tnew(:,:,:) = 0.0

    !!!!!!!!!!!!!!!!!!!!!!!!!!! new semi mesh ICs !!!!!!!!!!!!!!!!!!!!!!!!!!!!
    do un_ele=1,totele_unst
      if ( meshList(un_ele)%region_id == 13 ) then
        tnew(:,:,un_ele) = 1.
      end if
    end do

    dt = CFL*dx
    ! dt = CFL*dx/u_x + CFL*dy/u_y
    ntime = time/dt
    with_time_slab =.false.
    D3=.false.
    npoly=1; ele_type= 101

    call TRIQUAold(L1, L2, L3, L4, WEIGHT, D3,NGI)
    call get_surface_ele(meshlist, n_split, totele_unst, nface)
    call make_sparse_matrix(sparse_mass, totele_str*totele_unst, nloc)
    call make_sparse_matrix_flux_semi(sparse_lhs,meshList,n_split,totele_unst,totele_str,nloc,nface,str_neig)
    call make_sparse_matrix_flux_semi(sparse_flux,meshList,n_split,totele_unst,totele_str,nloc,nface,str_neig)
    call get_shape_funs_spec_ustr(n, nlx, nlx_lxx, nlxx, sn_orig, weight, nlx_nod, &
                nloc, snloc, sngi, ngi, ndim, nface,n_s_list_no, face_sn, face_sn2, face_snlx, &
                sweight, npoly, ele_type, totele_unst, face_list_no )

    call get_shape_funs_spec(n, nlx, nlx_lxx, nlxx, weight, nlx_nod, &
              nloc, snloc, sngi, ngi, ndim, nface,n_s_list_no, face_sn_str, face_sn2_str, face_snlx_str, &
              sweight_str, npoly, ele_type, totele_str)!, face_list_no)
    ! ! saving x_all structured triangles to be used for VTK generation
    i=1
    do un_ele=1,totele_unst
      do str_ele=1,totele_str
        call get_splitting(meshList(un_ele)%X, n_split, str_ele, str_x)
        x_all_str(:,:,i) = str_x !(ndim,nloc,totele)
        i=i+1
      end do
    end do
    print*, 'totele = ', totele_str * totele_unst
    print*, 'ntime = ', ntime
    call CPU_TIME(t1)

    do itime=1,ntime

      ! generating VTK files
      if ( vtk_io <= itime ) then
        call get_vtk(x_all_str, tnew, totele_str*totele_unst, nloc, itime,ndim, cell_type)
        vtk_io = vtk_io + vtk_interval
      end if

      told = tnew ! prepare for next time step
      tnew_nonlin = tnew
      do its=1,nits
        tnew = tnew_nonlin ! for non-linearity
        call update_overlaps(meshlist, tnew, told, t_bc, n_split, nface, totele_unst, nloc, str_neig)
        sparse_mass%val=0.
        sparse_lhs%val=0.
        sparse_flux%val=0.

        do un_ele=1,totele_unst
          do str_ele=1,totele_str
            call get_str_info(n_split, str_ele, irow, ipos, orientation)
            call get_splitting(meshList(un_ele)%X, n_split, str_ele, x_loc)
            call tri_det_nlx( x_loc, n, nlx, nx, detwei, weight, ndim, nloc, ngi, inv_jac )
            !volume=sum(detwei)
            u_loc = meshList(un_ele)%u_ele(:,:,str_ele) ! this is the
            tnew_loc =  tnew(:,str_ele, un_ele) ! this is the FEM element - 2nd index is the local node number
            told_loc =  told(:,str_ele, un_ele) ! this is the FEM element - 2nd index is the local node number

            do gi=1,ngi
              do idim=1,ndim
                ugi(gi,idim)=sum(n(gi,:)*u_loc(idim,:))
                tnew_xgi(gi,idim)=sum(nx(gi,idim,:)*tnew_loc(:))
              end do
              tnew_gi(gi)=sum(n(gi,:)*tnew_loc(:))
              told_gi(gi)=sum(n(gi,:)*told_loc(:))
              rgi(gi)=(tnew_gi(gi)-told_gi(gi))/dt + sum(ugi(gi,:)*tnew_xgi(gi,:))
              ! a_starif(ele_type < is_triangle_or_tet) then
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
            ! mass_ele=0.0 ! initialize mass matric to be zero.
            stab=0.0
            do iloc=1,nloc
              g_iloc = glob_no_semi(un_ele, str_ele, n_split, nloc, iloc)
              do jloc=1,nloc
                g_jloc = glob_no_semi(un_ele, str_ele, n_split, nloc, jloc)
                mass_ele = 0.0
                stiff_ele = 0.0
                do gi=1,ngi
                  stab(iloc,jloc) = stab(iloc,jloc) + sum(diff_coe(gi)* nx(gi,:,jloc)* nx(gi,:,iloc))* detwei(gi)
                  ! mass_ele(iloc,jloc) = mass_ele(iloc,jloc) + n(gi,iloc)*n(gi,jloc)*detwei(gi)
                  mass_ele = mass_ele + n(gi,iloc)*n(gi,jloc)*detwei(gi)/dt
                  do idim=1,ndim
                    stiff_ele = stiff_ele + nx(gi,idim,iloc)*ugi(gi,idim)*detwei(gi)*n(gi,jloc)
                  end do
                end do ! quadrature
                call add_to_CSR(sparse_mass, g_iloc, g_jloc, mass_ele)
                lhs = mass_ele-stiff_ele
                call add_to_CSR_flux(sparse_lhs, g_iloc, g_jloc, lhs)
              end do ! jloc
            end do ! iloc


            do iface = 1,nface
              ele22 = str_neig(iface,str_ele)!face_ele( iface, ele) ! if ele22 is 0 assume on boundary
              i_got_boundary = (sign(1, -ele22) + 1 )/2
              r_got_boundary = real(i_got_boundary)
              ele2 = str_ele*i_got_boundary + ele22*(1-i_got_boundary)! if ele22=0, then ele2=ele, if ele22/=0, ele2=ele22

              sn   = face_sn_str(:,:,iface)
              snlx = face_snlx_str(:,:,:,iface)

              sn2=0
              if ( ele22==0 .and. iface==1 ) then
                sp = ipos/2+1 ! position of ele along the un_iface
                mface = 1 ! un_iface number which str_ele is located on
                call get_semi_sn2_implicit(meshlist, un_ele, mface, sn_orig, sn2)

              else if ( ele22==0 .and. iface==2 ) then
                sp = irow
                mface = 3
                call get_semi_sn2_implicit(meshlist, un_ele, mface, sn_orig, sn2)

              else if ( ele22==0 .and. iface==3 ) then
                sp = irow
                mface = 2
                call get_semi_sn2_implicit(meshlist, un_ele, mface, sn_orig, sn2)

              else if ( ele22 /= 0 ) then
                sn2  = face_sn2_str(:,:,iface)
              end if ! end if ele22=0

              tnew_loc2(:)= tnew(:,ele2, un_ele) * (1.0-r_got_boundary) &
                          + meshlist(un_ele)%t_overlap( (sp)*nloc-2:(sp)*nloc, mface ) * r_got_boundary

              u_loc2(:,:)= meshList(un_ele)%u_ele(:,:,ele2)* (1.0-r_got_boundary) &
                          + meshlist(un_ele)%u_overlap(:, sp*nloc-2:sp*nloc, mface ) * r_got_boundary
                          ! 2nd section of u_overlap should refer ro 3 ilocs from the un_ele u_bc which can be identified based on irow

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

              do idim=1,ndim
                s_cont(:,idim) = snorm(:,idim)*sdetwei(:) &
                            *( (1.-income(:))* usgi(:,idim)*tnew_sgi(:) + income(:)*usgi2(:,idim) )
              end do

              target_ele = int( (1.-income(1))*str_ele + income(1)*ele22 )
              Npos = meshList(un_ele)%Neig(mface)
              Nside = meshList(un_ele)%fNeig(mface)
              if ( target_ele==0 .and. Npos==0 ) then ! it happens when ele is on the boundary
              elseif ( target_ele==0 .and. Npos/=0) then
                        if ( iface==1 ) then
                              ipos = int(ipos/2)+1
                              target_ele = meshList(un_ele)%s_ele(ipos, 1)

                        elseif ( iface==2 ) then
                            target_ele = meshList(un_ele)%s_ele(irow, 3)

                          elseif ( iface==3 ) then
                              target_ele = meshList(un_ele)%s_ele(irow, 2)
                          end if

              elseif ( target_ele /= 0) then
                target_ele = (un_ele-1)*4**n_split + target_ele
              end if
              do iloc=1,nloc
                g_iloc = glob_no_semi(un_ele, str_ele, n_split, nloc, iloc)
                do jloc=1,nloc
                  g_jloc = (target_ele-1)*nloc + jloc !glob_no_semi(un_ele, target_ele, n_split, nloc, jloc)
                  flux_ele = 0.0
                  do idim=1,ndim
                    ! glob_flux(g_iloc,g_jloc)=glob_flux(g_iloc,g_jloc) +sum( snorm(:,idim)*sdetwei(:)*sn(:,iloc) &
                    !                                                   *((1.-income(:))*sn(:,jloc)*usgi(:,idim) &
                    !                                                   +     income(:)*sn2(:,jloc)*usgi2(:,idim)))
                    flux_ele = flux_ele +sum( snorm(:,idim)*sdetwei(:)*sn(:,iloc) &
                                        *((1.-income(:))*sn(:,jloc)*usgi(:,idim) &
                                        +     income(:)*sn2(:,jloc)*usgi2(:,idim)))
                  end do ! idim

                  call add_to_CSR_flux(sparse_flux, g_iloc, g_jloc,  flux_ele)
                end do ! jloc
              end do ! iloc
            end do ! iface

            ! if (with_stab) then
              ! mat_loc= mass_ele + dt*stab
            ! else
            !   stab=0
            !   mat_loc= mass_ele
            ! end if
          end do ! end do totele_str
        end do ! end totele_unst
        call gen_global_matrix(totele_str*totele_unst, nloc, sparse_mass, glob_mass)
        call gen_global_matrix(totele_str*totele_unst, nloc, sparse_lhs, glob_lhs)
        call gen_global_matrix(totele_str*totele_unst, nloc, sparse_flux, glob_flux)
        call told_to_array_semi(totele_str, totele_unst, nloc, told, glob_told_array)
        rhs_array = matmul(glob_mass, glob_told_array)
        glob_lhs = glob_lhs + glob_flux

        call FINDInv(glob_lhs, inv_lhs, nloc*totele_str*totele_unst, errorflag)
        tnew_nonlin2 = matmul(inv_lhs,rhs_array )
        call array_to_told_semi(totele_str, totele_unst, nloc, tnew_nonlin2, tnew_nonlin)
        tnew=tnew_nonlin
      end do ! do its=1,nits
    print*, 'semi', itime
    end do ! do itime=1,ntime

    call CPU_TIME(t2)
    print*, 'cpu_time for time_loop = ', t2-t1

    deallocate(face_ele, face_list_no, str_neig, n, nlx, nx, M, nlx_lxx, nlxx, nlx_nod, weight, detwei, sdetwei,&
               sn,sn2,snlx, sweight, s_cont, sweight_str, tnew_loc, tnew_loc2, told, tnew, tnew_nonlin,&
               tnew_nonlin_loc, t_bc, u_bc, tnew_gi, tnew_xgi, tnew_sgi, tnew_sgi2, told_gi,&
               face_sn, face_sn2, face_snlx, face_sn_str, face_sn2_str, face_snlx_str, u_loc, u_loc2, x_loc, ugi,&
               xsgi, usgi, usgi2, income, snorm, norm, mat_loc_inv, rhs_loc, ml_ele, rhs_jac, mass_t_new,&
               SN_orig,SNLX_orig, inv_jac, mat_diag_approx, a_star, p_star, rgi, diff_coe, told_loc, stab, mat_tnew, mass_told,&
               L1, L2, L3, L4, x_all_str, inv_lhs, rhs_array, glob_flux,glob_mass,glob_told_array,glob_lhs,&
               tnew_nonlin2, sparse_flux%g_iloc, sparse_flux%g_jloc, sparse_flux%val,&
               sparse_lhs%g_iloc, sparse_lhs%g_jloc, sparse_lhs%val, sparse_mass%g_iloc, sparse_mass%g_jloc,sparse_mass%val)
  end subroutine Semi_implicit_direct


  ! brief> this is for semi-structured explicit method
  Subroutine semi_explicit(CFL, ndim, direct_solver, njac_its, time, nits,&
    nface, dx, dy, nloc, snloc, ngi, sngi, n_s_list_no, u_x, u_y, with_stab, vtk_interval, cell_type)
    implicit none

    ! external variables
    real, intent(in) :: CFL, time, u_x, u_y, dx, dy
    integer, intent(in) :: ndim, nits, nface, nloc, ngi, sngi, snloc, njac_its
    integer, intent(in) :: vtk_interval, cell_type, n_s_list_no
    logical, intent(in) :: direct_solver, with_stab

    ! local vbls
    type(Mesh), allocatable, dimension(:) :: meshList
    ! type(Semi_Mesh), allocatable, dimension(:) :: semi_meshlist

    character (len=20):: file_name
    logical :: with_time_slab, D3

    integer :: gi, iface, totele, ele, ele2, ele22, n_split, irow, ipos, sum_up, irow_ups
    integer :: totele_unst, un_ele, str_ele, totele_str
    integer :: s_list_no, s_gi, iloc, jloc, i,j,k, ntime, vtk_io
    integer :: itime, idim, sndim, nonodes, mloc, col, row, row2
    integer :: errorflag, i_got_boundary, jac_its, its
    integer :: npoly, ele_type, ierr, i_overlaps, orientation
    integer :: Mpos, Npos, Mside, Nside, lnodes(2), lnodes2(2), totnodes, mface, sp

    real :: sarea, volume, dt, L, r_got_boundary, toler, a_coef
    real :: length_x, length_y
    real :: t1, t2, t1_get_shape_funs_spec, t2_get_shape_funs_spec
    real :: t1_tri_det_nlx, t2_tri_det_nlx, time_tri_det_nlx, i_tri_det_nlx
    real :: t1_det_snlx_all, t2_det_snlx_all, time_det_snlx_all, i_det_snlx_all
    real :: t1_ReadMSH, t2_ReadMSH, time_ReadMSH
    real :: t1_getNeigDataMesh, t2_getNeigDataMesh, time_getNeigDataMesh
    real :: t1_get_unstr_sn2, t2_get_unstr_sn2, time_get_unstr_sn2, str_x(2,3)

    integer, allocatable ::  face_ele(:,:), face_list_no(:,:), str_neig(:,:)

    real, allocatable :: n(:,:), nlx(:,:,:), nx(:,:,:), M(:,:)
    real, allocatable :: nlx_lxx(:,:,:), nlxx(:,:), nlx_nod(:,:,:)
    real, allocatable :: weight(:), detwei(:), sdetwei(:)
    real, allocatable :: sn(:,:),sn2(:,:),snlx(:,:,:), sweight(:), s_cont(:,:), sweight_str(:)
    real, allocatable :: tnew_loc(:), tnew_loc2(:), told(:,:,:), tnew(:,:,:), tnew_nonlin(:,:,:)
    real, allocatable :: tnew_nonlin_loc(:)
    real, allocatable :: t_bc(:,:), u_bc(:,:,:)
    real, allocatable :: tnew_gi(:), tnew_xgi(:,:), tnew_sgi(:), tnew_sgi2(:), told_gi(:)
    real, allocatable :: face_sn(:,:,:), face_sn2(:,:,:), face_snlx(:,:,:,:)
    real, allocatable :: face_sn_str(:,:,:), face_sn2_str(:,:,:), face_snlx_str(:,:,:,:)
    real, allocatable :: u_loc(:,:), u_ele(:,:,:), u_loc2(:,:), x_loc(:,:), ugi(:,:)
    real, allocatable :: xsgi(:,:), usgi(:,:), usgi2(:,:), income(:), snorm(:,:), norm(:)
    real, allocatable :: mass_ele(:,:), mat_loc(:,:), mat_loc_inv(:,:), rhs_loc(:), ml_ele(:), rhs_jac(:), mass_t_new(:)
    real, allocatable :: SN_orig(:,:),SNLX_orig(:,:,:), inv_jac(:, :, : ), mat_diag_approx(:)
    real, allocatable :: a_star(:,:), p_star(:), rgi(:), diff_coe(:), told_loc(:), stab(:,:), mat_tnew(:), mass_told(:)
    real, allocatable :: L1(:), L2(:), L3(:), L4(:), x_all_str(:,:,:)


    ! xp(2): coordinates of nodes
    ! call ReadMSH(meshList,'./semi_mesh.msh',ierr, totnodes)
    ! call ReadMSH(meshList,'./irregular shape.msh',ierr, totnodes)
    ! call ReadMSH(meshList,'./new_semi_mesh.msh',ierr, totnodes)
    ! call ReadMSH(meshList,'./gmsh_100.msh',ierr, totnodes)
    print*, 'starting reading the .msh file'
    call CPU_TIME(t1_ReadMSH)
    ! call ReadMSH(meshList,'./10_elements.msh',ierr, totnodes)
    ! call ReadMSH(meshList,'./1_unele_test.msh',ierr, totnodes)
    call ReadMSH(meshList,'./2_unele_test.msh',ierr, totnodes)
    call CPU_TIME(t2_ReadMSH)
    time_ReadMSH = t2_ReadMSH - t1_ReadMSH
    print*, 'finished reading .msh file'

    totele_unst = size(meshList)

    n_split = 2
    totele_str = 4**n_split

    call get_str_neig(n_split, str_neig)

    sndim = ndim-1 ! dimensions of the surface element coordinates.
    mloc=1
    toler = 0.00000000001
    ele_type = 3
    vtk_io=vtk_interval

    allocate(n( ngi, nloc ), nx( ngi, ndim, nloc ), nlx( ngi, ndim, nloc ), M(ngi,nloc))
    allocate(nlx_lxx(ngi,ndim,nloc), nlxx(ngi,nloc), nlx_nod(nloc,ndim,nloc))
    allocate(weight(ngi), detwei(ngi), sdetwei(sngi))
    allocate(face_list_no( nface, totele_str), face_ele(nface,totele_str))
    allocate(sn(sngi,nloc),sn2(sngi,nloc),snlx(sngi,sndim,nloc))!,sweight(sngi))
    allocate(SN_orig(sngi,snloc),SNLX_orig(sngi,sndim,snloc), s_cont(sngi,ndim) )
    allocate(u_loc(ndim,nloc), u_loc2(ndim,nloc), ugi(ngi,ndim))
    allocate(xsgi(sngi,ndim), usgi(sngi,ndim), usgi2(sngi,ndim), income(sngi), snorm(sngi,ndim), norm(ndim))
    allocate(tnew_loc(nloc), tnew_loc2(nloc), told(nloc,totele_str,totele_unst))
    allocate(tnew(nloc,totele_str,totele_unst), tnew_nonlin(nloc,totele_str,totele_unst))
    allocate(tnew_nonlin_loc(nloc))
    allocate(t_bc(2**n_split*nloc,4), u_bc(ndim,2**n_split*nloc,4)) ! it works just for semi_mesh.msh
    allocate(tnew_gi(ngi), tnew_xgi(ngi,ndim), tnew_sgi(sngi), tnew_sgi2(sngi), told_gi(ngi))
    allocate(face_sn(sngi,nloc,nface), face_sn2(sngi,nloc,n_s_list_no), face_snlx(sngi,sndim,nloc,nface) )
    allocate(face_sn_str(sngi,nloc,nface), face_sn2_str(sngi,nloc,n_s_list_no), face_snlx_str(sngi,sndim,nloc,nface))
    allocate(x_loc(ndim,nloc))
    allocate(mass_ele(nloc,nloc), mat_loc_inv(nloc,nloc), rhs_loc(nloc),ml_ele(nloc))
    allocate(rhs_jac(nloc), mass_t_new(nloc), inv_jac(ngi, ndim, ndim ) )
    allocate(a_star(ngi,nloc), p_star(ngi), rgi(ngi), diff_coe(ngi), told_loc(nloc), stab(nloc, nloc))
    allocate(mat_tnew(nloc), mat_diag_approx(nloc), mass_told(nloc))
    allocate(L1(ngi), L2(ngi), L3(ngi), L4(ngi), x_all_str(ndim,nloc,totele_unst*totele_str))

    do un_ele=1,size(meshlist)
      allocate(meshList(un_ele)%t_overlap(2**n_split*nloc, nface))
      allocate(meshList(un_ele)%t_overlap_old(2**n_split*nloc, nface))
      allocate(meshList(un_ele)%u_overlap(ndim, 2**n_split*nloc , nface))
      allocate(meshList(un_ele)%u_ele(ndim, nloc, 4**n_split))
      allocate(meshList(un_ele)%S_nodes(sngi,nface))
      meshList(un_ele)%t_overlap=0.0
      meshList(un_ele)%u_overlap(1,:,:)=u_x
      meshList(un_ele)%u_overlap(2,:,:)=u_y
      meshList(un_ele)%u_ele(:,:,:) = 0
      meshList(un_ele)%u_ele(1,:,:) = u_x ! suggest this should be unity
      meshList(un_ele)%u_ele(2,:,:) = u_y ! suggest this should be unity
      do iface=1,nface

        call getNeigDataMesh(meshlist, un_ele, iface, Npos, Nside, lnodes2)
        meshList(un_ele)%S_nodes(:,iface)=lnodes2
      end do
    end do

    meshList(1)%t_overlap=0.0
    ! meshList(2)%t_overlap=0.0

    t_bc(:,:)=0.0 ! this contains the boundary conditions just outside the domain
    u_bc(1,:,:)=u_x ! this contains the boundary conditions on velocity just outside the domain
    u_bc(2,:,:)=u_y
    tnew(:,:,:) = 0.0

    tnew(:,9,1)=1.0
    tnew(:,10,1)=1.0


    dt = CFL*dx
    ! dt = CFL*dx/u_x + CFL*dy/u_y
    ntime = time/dt
    with_time_slab =.false.
    D3=.false.
    npoly=1; ele_type= 101

    call TRIQUAold(L1, L2, L3, L4, WEIGHT, D3,NGI)

    call get_shape_funs_spec_ustr(n, nlx, nlx_lxx, nlxx, sn_orig, weight, nlx_nod, &
                nloc, snloc, sngi, ngi, ndim, nface,n_s_list_no, face_sn, face_sn2, face_snlx, &
                sweight, npoly, ele_type, totele_unst, face_list_no )

    call get_shape_funs_spec(n, nlx, nlx_lxx, nlxx, weight, nlx_nod, &
              nloc, snloc, sngi, ngi, ndim, nface,n_s_list_no, face_sn_str, face_sn2_str, face_snlx_str, &
              sweight_str, npoly, ele_type, totele_str)!, face_list_no)

    ! ! saving x_all structured triangles to be used for VTK generation
    i=1
    do un_ele=1,totele_unst
      do str_ele=1,totele_str
        call get_splitting(meshList(un_ele)%X, n_split, str_ele, str_x)
        x_all_str(:,:,i) = str_x !(ndim,nloc,totele)
        i=i+1
      end do
    end do
    print*, 'totele = ', totele_str * totele_unst
    print*, 'ntime = ', ntime
    call CPU_TIME(t1)


    do itime=1,ntime

      ! generating VTK files
      if ( vtk_io <= itime ) then
        call get_vtk(x_all_str, tnew, totele_str*totele_unst, nloc, itime,ndim, cell_type)
        vtk_io = vtk_io + vtk_interval
      end if

      told = tnew ! prepare for next time step
      tnew_nonlin = tnew
      do its=1,nits
        tnew = tnew_nonlin ! for non-linearity
        call update_overlaps(meshlist, tnew, told, t_bc, n_split, nface, totele_unst, nloc, str_neig)
        do un_ele=1,totele_unst

          do str_ele=1,totele_str
            call get_str_info(n_split, str_ele, irow, ipos, orientation)
            call get_splitting(meshList(un_ele)%X, n_split, str_ele, x_loc)
            ! x_loc = str_x ! x_all contains the coordinates of the corner nodes

            call tri_det_nlx( x_loc, n, nlx, nx, detwei, weight, ndim, nloc, ngi, inv_jac )

            !volume=sum(detwei)
            u_loc = meshList(un_ele)%u_ele(:,:,str_ele) ! this is the
            tnew_loc =  tnew(:,str_ele, un_ele) ! this is the FEM element - 2nd index is the local node number
            told_loc =  told(:,str_ele, un_ele) ! this is the FEM element - 2nd index is the local node number

            do gi=1,ngi
              do idim=1,ndim
                ugi(gi,idim)=sum(n(gi,:)*u_loc(idim,:))
                tnew_xgi(gi,idim)=sum(nx(gi,idim,:)*tnew_loc(:))
              end do
              tnew_gi(gi)=sum(n(gi,:)*tnew_loc(:))
              told_gi(gi)=sum(n(gi,:)*told_loc(:))
              rgi(gi)=(tnew_gi(gi)-told_gi(gi))/dt + sum(ugi(gi,:)*tnew_xgi(gi,:))
              ! a_starif(ele_type < is_triangle_or_tet) then
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
            mass_ele=0.0 ! initialize mass matric to be zero.
            stab=0.0
            do iloc=1,nloc
              !inod = glob_no(iloc,ele)
              do jloc=1,nloc
                !jnod = glob_no(ele,jloc)
                do gi=1,ngi
                  stab(iloc,jloc) = stab(iloc,jloc) + sum(diff_coe(gi)* nx(gi,:,jloc)* nx(gi,:,iloc))* detwei(gi)
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


            do iface = 1,nface
              ele22 = str_neig(iface,str_ele)!face_ele( iface, ele) ! if ele22 is 0 assume on boundary
              i_got_boundary = (sign(1, -ele22) + 1 )/2
              r_got_boundary = real(i_got_boundary)
              ele2 = str_ele*i_got_boundary + ele22*(1-i_got_boundary)! if ele22=0, then ele2=ele, if ele22/=0, ele2=ele22

              ! s_list_no = face_list_no( iface, ele)
              sn   = face_sn_str(:,:,iface)
              snlx = face_snlx_str(:,:,:,iface)

              sn2=0
              if ( ele22==0 .and. iface==1 ) then
                sp = ipos/2+1 ! position of ele along the un_iface
                mface = 1 ! un_iface number which str_ele is located on
                call get_semi_sn2(meshlist, un_ele, mface, sn_orig, sn2)

              else if ( ele22==0 .and. iface==2 ) then
                sp = irow
                mface = 3
                call get_semi_sn2(meshlist, un_ele, mface, sn_orig, sn2)

              else if ( ele22==0 .and. iface==3 ) then
                sp = irow
                mface = 2
                call get_semi_sn2(meshlist, un_ele, mface, sn_orig, sn2)

              else if ( ele22 /= 0 ) then
                sn2  = face_sn2_str(:,:,iface)
              end if ! end if ele22=0

              tnew_loc2(:)= tnew(:,ele2, un_ele) * (1.0-r_got_boundary) &
                          + meshlist(un_ele)%t_overlap( (sp)*nloc-2:(sp)*nloc, mface ) * r_got_boundary

              u_loc2(:,:)= meshList(un_ele)%u_ele(:,:,ele2)* (1.0-r_got_boundary) &
                          + meshlist(un_ele)%u_overlap(:, sp*nloc-2:sp*nloc, mface ) * r_got_boundary
                          ! 2nd section of u_overlap should refer ro 3 ilocs from the un_ele u_bc which can be identified based on irow

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

            if (with_stab) then
              mat_loc= mass_ele + dt*stab
            else
              stab=0
              mat_loc= mass_ele
            end if

            if(direct_solver) then
              ! inverse of the mass matric (nloc,nloc)
              call FINDInv(mass_ele, mat_loc_inv, nloc, errorflag)
              !can be only this part instead
              do iloc=1,nloc
                 tnew_nonlin(iloc,str_ele, un_ele)=told_loc(iloc) + dt* sum( mat_loc_inv(iloc,:)* rhs_loc(:) )! &
                                               !+ told_loc(iloc)* dt* sum( mat_loc_inv(iloc,:) * stab(iloc,:) )
              end do
              ! do iloc=1,nloc
              !   mass_told(iloc)=sum( mat_loc(iloc,:) * told_loc(:) )
              ! end do
              ! do iloc=1,nloc
              !   tnew_nonlin(iloc,str_ele, un_ele)=sum( mat_loc_inv(iloc,:)* (mass_told(:)+dt*rhs_loc(:)))
              ! end do

            else ! iterative solver
               do iloc=1,nloc
                  rhs_jac(iloc)= sum( mass_ele(iloc,:) * told_loc(:) ) !+ dt*rhs_loc(iloc)
                  mat_diag_approx(iloc) = ml_ele(iloc) - dt*rhs_loc(iloc) + dt * stab(iloc,iloc)
               end do

               tnew_nonlin_loc(:) = tnew_nonlin(:,str_ele, un_ele)
               do jac_its=1,njac_its ! jacobi iterations...
                  do iloc=1,nloc
                     mat_tnew(iloc)= sum( mat_loc(iloc,:) * tnew_nonlin_loc(:) )- dt*rhs_loc(iloc)
                  end do
                  do iloc=1,nloc
                     tnew_nonlin_loc(iloc)=  (mat_diag_approx(iloc)*tnew_nonlin_loc(iloc) - mat_tnew(iloc) + rhs_jac(iloc) )&
                                                                    / mat_diag_approx(iloc)
                  end do
               end do
               tnew_nonlin(:,str_ele, un_ele) = tnew_nonlin_loc(:)
            endif ! endof if(direct_solver) then else

          end do ! end do totele_str
        end do ! end totele_unst
        tnew=tnew_nonlin
      end do ! do its=1,nits
      print*, 'semi-structured ', itime
    end do ! do itime=1,ntime

    call CPU_TIME(t2)
    print*, 'cpu_time for time_loop = ', t2-t1

    deallocate(face_ele, face_list_no, str_neig, n, nlx, nx, M, nlx_lxx, nlxx, nlx_nod, weight, detwei, sdetwei,&
               sn,sn2,snlx, sweight, s_cont, sweight_str, tnew_loc, tnew_loc2, told, tnew, tnew_nonlin,&
               tnew_nonlin_loc, t_bc, u_bc, tnew_gi, tnew_xgi, tnew_sgi, tnew_sgi2, told_gi,&
               face_sn, face_sn2, face_snlx, face_sn_str, face_sn2_str, face_snlx_str, u_loc, u_loc2, x_loc, ugi,&
               xsgi, usgi, usgi2, income, snorm, norm, mass_ele, mat_loc, mat_loc_inv, rhs_loc, ml_ele, rhs_jac, mass_t_new,&
               SN_orig,SNLX_orig, inv_jac, mat_diag_approx, a_star, p_star, rgi, diff_coe, told_loc, stab, mat_tnew, mass_told,&
               L1, L2, L3, L4, x_all_str)
  end subroutine semi_explicit

end module transport_tri_semi
