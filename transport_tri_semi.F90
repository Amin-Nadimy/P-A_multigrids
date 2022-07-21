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
    ! solver:: defines type of solvers, 1 :: Jacobi
    !                                   2 :: Richardon
    Subroutine Semi_implicit_iterative(CFL, ndim, direct_solver, n_multigrid, time, nits, nface, dx, dy,&
              nloc, snloc, ngi, sngi, n_s_list_no, u_x, u_y, with_stab, vtk_interval, cell_type, solver, multi_levels,n_smooth)
      implicit none

      ! global variables
      real, intent(in) :: CFL, time, u_x, u_y, dx, dy
      integer, intent(in) :: ndim, nits, nface, nloc, ngi, sngi, snloc, n_multigrid, solver
      integer, intent(in) :: vtk_interval, cell_type, n_s_list_no, multi_levels,n_smooth
      logical, intent(in) :: direct_solver, with_stab

      ! local vbls
      type(Mesh), allocatable, dimension(:), target :: meshList
      type(pointer_Mesh), target :: mesh_pnt
      type(fields), allocatable, dimension(:), target :: tracer
      type(sparse) :: sparse_mass, sparse_lhs, sparse_flux
      type(element_info), allocatable, target :: ele_info(:)

      character (len=20):: file_name, solve_for
      logical :: with_time_slab, D3

      integer :: gi, iface, ele, ele2, n_split, irow, ipos, sum_up, irow_ups, i_split
      integer :: totele_unst, un_ele, str_ele, totele_str,updown
      integer :: s_list_no, s_gi, iloc, jloc, i,j, ntime, vtk_io,smooth
      integer :: itime, idim, sndim, nonodes, mloc, col, row, row2
      integer :: errorflag, i_got_boundary, multigrid, its, ilevel
      integer :: npoly, ele_type, ierr, i_overlaps, target_loc
      integer :: Mpos, Npos, Mside, Nside, lnodes(2), lnodes2(2), totnodes, mface, sp
      integer :: g_iloc, g_jloc, target_ele, orientation
      integer, pointer :: ele22

      real :: k, delta_x,omga

      real :: sarea, volume, dt, L
      real :: length_x, length_y, theta
      real :: t1, t2, t1_get_shape_funs_spec, t2_get_shape_funs_spec
      real :: t1_tri_det_nlx, t2_tri_det_nlx, time_tri_det_nlx, i_tri_det_nlx
      real :: t1_det_snlx_all, t2_det_snlx_all, time_det_snlx_all, i_det_snlx_all
      real :: t1_ReadMSH, t2_ReadMSH, time_ReadMSH
      real :: t1_getNeigDataMesh, t2_getNeigDataMesh, time_getNeigDataMesh
      real :: t1_get_unstr_sn2, t2_get_unstr_sn2, time_get_unstr_sn2, str_x(2,3),mat_diag(nloc,nloc)

      logical, parameter :: activate_advection = .false.
      logical, parameter :: activate_diffusion = .true.
      logical, parameter :: activate_time = .false.



      integer, allocatable ::  face_ele(:,:), face_list_no(:,:), res(:), result(:),surf_ele(:,:)
      integer, allocatable :: surf_ele2(:),fin_ele(:)

      real, allocatable :: n(:,:), nlx(:,:,:), nx(:,:,:), M(:,:)
      real, allocatable :: nlx_lxx(:,:,:), nlxx(:,:), nlx_nod(:,:,:)
      real, allocatable :: weight(:)
      real, allocatable :: sn2(:,:),snlx(:,:,:), sweight(:), sweight_str(:)
      real, allocatable :: t_bc(:), u_bc(:,:,:)
      real, allocatable :: tnew_gi(:), tnew_xgi(:,:), tnew_sgi(:), tnew_sgi2(:), told_gi(:),told_sgi(:),told_sgi2(:)
      real, allocatable :: face_sn(:,:,:), face_sn2(:,:,:), face_snlx(:,:,:,:)
      real, allocatable :: face_sn2_str(:,:,:), face_snlx_str(:,:,:,:)
      real, allocatable :: u_ele(:,:,:), x_loc(:,:), ugi(:,:)
      real, allocatable :: xsgi(:,:), usgi(:,:), usgi2(:,:), income(:), snorm(:,:), norm(:)
      real, allocatable :: mat_loc(:,:), ml_ele(:), rhs_jac(:), mass_t_new(:)
      real, allocatable :: SN_orig(:,:),SNLX_orig(:,:,:), inv_jac(:, :, : ), mat_diag_approx(:)
      real, allocatable :: mat_tnew(:), mass_told(:)
      real, allocatable :: L1(:), L2(:), L3(:), L4(:), x_all_str(:,:,:)!, glob_told_array(:)
      real, allocatable :: A_x(:), mass_ele_old(:),stiff_ele_old(:)
      real, allocatable :: stiff_ele_new(:), s_cont(:,:),flux_ele_new(:), mass_ele_new(:)
      real, allocatable :: flux_ele_old(:),s_cont_old(:,:), diff_vol(:), diff_surf(:),analytical(:,:,:)
      real, allocatable :: mass_ele(:,:),diff_vol2(:,:,:),stiff_ele(:,:,:), tnew_loc2(:,:,:,:)
      real, allocatable :: mass_stcl(:,:),stiff_stcl(:,:,:,:),diff_vol_stcl(:,:,:,:)
      real, allocatable :: my_diff_surf(:,:,:),neig_diff_surf(:,:,:),stiff_ele_new1(:,:), stiff_ele_old1(:,:)
      real, allocatable :: tnew_nonlin_loc2(:,:,:,:), jac_res1(:)
      real, allocatable, target :: face_sn_str(:,:,:), tnew_nonlin(:,:,:),diff_vol1(:,:)
      real, pointer :: sn(:,:), sdetwei(:), u_loc2(:,:), detwei(:),u_loc(:,:)
      real, pointer :: told_loc2(:),tnew_nonlin_loc(:),told_loc(:),tnew_loc(:)

      ! xp(2): coordinates of nodes
      ! call ReadMSH(meshList,'./semi_mesh.msh',ierr, totnodes)
      ! call ReadMSH(meshList,'./irregular shape.msh',ierr, totnodes)
      ! call ReadMSH(meshList,'./new_semi_mesh.msh',ierr, totnodes)
      ! call ReadMSH(meshList,'./1_unele_test.msh',ierr, totnodes)
      ! call ReadMSH(meshList,'./gmsh_100.msh',ierr, totnodes)
      print*, '---------------------------------------------------------'
      print*, '|       Reading the .msh file     |'
      call CPU_TIME(t1_ReadMSH)
      ! call ReadMSH(meshList,'./1_unele_test.msh',ierr, totnodes)
      ! call ReadMSH(meshList,'./2_unele_test.msh',ierr, totnodes)
      call ReadMSH(meshList,'./new_2_ele.msh',ierr, totnodes)
      ! call ReadMSH(meshList,'./irregular.msh',ierr, totnodes)
      ! call ReadMSH(meshList,'./semi_structured_mesh.msh',ierr, totnodes)
      ! call ReadMSH(meshList,'./P_structured_mesh',ierr, totnodes)
      ! call ReadMSH(meshList,'./untitled8192.msh',ierr, totnodes) ! n==2
      ! call ReadMSH(meshList,'./untitled2048.msh',ierr, totnodes) ! n==3
      ! call ReadMSH(meshList,'./untitled8.msh',ierr, totnodes)  ! n==7

      call CPU_TIME(t2_ReadMSH)
      time_ReadMSH = t2_ReadMSH - t1_ReadMSH
      print*, '|   Time for reading .msh file    |', time_ReadMSH

      totele_unst = size(meshList)
      theta = 1.
      n_split = 3
      ! orientation = 1 ! 1 is up-triangle and -1 is down-triangle
      totele_str = 4**n_split
      solve_for = 'Tracer'



      sndim = ndim-1 ! dimensions of the surface element coordinates.
      mloc=1
      ele_type = 3
      vtk_io=vtk_interval
      dt = CFL*dx
      ! dt = CFL*dx/u_x + CFL*dy/u_y
      ntime = 10!time/dt
      k = 1. !diffusion coeficient for the diffusion term, m^2/s for water diffuses into air
      with_time_slab =.false.
      D3=.false.
      npoly=1; ele_type= 101
      omga = 0.8

      allocate(n( ngi, nloc ), nx( ngi, ndim, nloc ), nlx( ngi, ndim, nloc ), M(ngi,nloc))
      allocate(nlx_lxx(ngi,ndim,nloc), nlxx(ngi,nloc), nlx_nod(nloc,ndim,nloc))
      allocate(weight(ngi))
      allocate(face_list_no( nface, totele_str), face_ele(nface,totele_str))
      allocate(sn2(sngi,nloc),snlx(sngi,sndim,nloc))
      allocate(SN_orig(sngi,snloc),SNLX_orig(sngi,sndim,snloc), s_cont(sngi,ndim),s_cont_old(sngi,ndim) )
      allocate(ugi(ngi,ndim))
      allocate(xsgi(sngi,ndim), usgi(sngi,ndim), usgi2(sngi,ndim), income(sngi), snorm(sngi,ndim), norm(ndim))
      ! allocate(tnew_nonlin(nloc,totele_str,totele_unst))
      allocate(t_bc(2), u_bc(ndim,2**n_split*nloc,4)) ! it works just for semi_mesh.msh
      allocate(tnew_gi(ngi),tnew_xgi(ngi,ndim),tnew_sgi(sngi),tnew_sgi2(sngi),told_gi(ngi),told_sgi(sngi),told_sgi2(sngi))
      allocate(face_sn(sngi,nloc,nface), face_sn2(sngi,nloc,n_s_list_no), face_snlx(sngi,sndim,nloc,nface) )
      allocate(face_sn_str(sngi,nloc,nface), face_sn2_str(sngi,nloc,n_s_list_no), face_snlx_str(sngi,sndim,nloc,nface))
      allocate(x_loc(ndim,nloc))
      allocate(ml_ele(nloc))
      allocate(rhs_jac(nloc), mass_t_new(nloc), inv_jac(ngi, ndim, ndim ), tnew_loc2(nloc,totele_str,totele_unst,nface) )
      ! allocate(tnew_nonlin_loc2(nloc,totele_str,totele_unst,nface))
      allocate(mat_tnew(nloc), mat_diag_approx(nloc), mass_told(nloc))
      allocate(L1(ngi), L2(ngi), L3(ngi), L4(ngi), x_all_str(ndim,nloc,totele_unst*totele_str))
      allocate(A_x(nloc),diff_vol(nloc),fin_ele(4))
      allocate(stiff_ele_old(nloc), diff_surf(nloc), jac_res1(nloc))
      allocate(stiff_ele_new(nloc),flux_ele_new(nloc), flux_ele_old(nloc), mass_ele_new(nloc),mass_ele_old(nloc))
      allocate(mass_ele(ngi,nloc),diff_vol2(ngi,ndim,nloc),stiff_ele(ngi,ndim,nloc),surf_ele(2**n_split,nface))
      allocate(surf_ele2(3*2**n_split),mass_stcl(nloc,nloc),stiff_stcl(nloc,ngi,ndim,nloc),diff_vol_stcl(ngi,ndim,nloc,nloc))
      allocate(diff_vol1(nloc,nloc),my_diff_surf(nloc,nloc,nface))
      allocate(neig_diff_surf(nloc,nloc,nface))
      allocate(stiff_ele_new1(nloc,nloc),stiff_ele_old1(nloc,nloc),analytical(nloc,totele_str,totele_unst))

      allocate(tracer(multi_levels), ele_info(multi_levels))
      do i=1,multi_levels
        totele_str = 4**(n_split-i+1)
        allocate(tracer(i)%tnew(nloc,totele_str,totele_unst))
        allocate(tracer(i)%told(nloc,totele_str,totele_unst))
        allocate(tracer(i)%RHS(nloc,totele_str,totele_unst))
        allocate(tracer(i)%residuale(nloc,totele_str,totele_unst))
        allocate(tracer(i)%error(nloc,totele_str,totele_unst))
        allocate(tracer(i)%source(nloc,totele_str,totele_unst))

        allocate(ele_info(i)%surf_ele(2**(n_split-i+1),nface))
        call loc_surf_ele_multigrid(ele_info, i, n_split-i+1)

        allocate(ele_info(i)%str_neig(3,4**(n_split-i+1)))
        call get_str_neig_multigrid(ele_info, i, n_split-i+1)

      end do

      do un_ele=1,size(meshlist)
        allocate(meshList(un_ele)%t_overlap(2**n_split*nloc, nface))
        allocate(meshList(un_ele)%t_overlap_old(2**n_split*nloc, nface))
        allocate(meshList(un_ele)%u_overlap(ndim, 2**n_split*nloc , nface))
        allocate(meshList(un_ele)%u_ele(ndim, nloc, 4**n_split))
        allocate(meshList(un_ele)%S_nodes(sngi,nface))
        allocate(meshList(un_ele)%s_ele(2**n_split,nface))
        allocate(meshList(un_ele)%detwei(ngi))
        allocate(meshList(un_ele)%sdetwei(sngi,nface))
        allocate(meshList(un_ele)%snorm(sngi,ndim,nface))
        allocate(meshList(un_ele)%nx(ngi,ndim,nloc))
        meshList(un_ele)%t_overlap=0.0
        meshList(un_ele)%u_overlap(1,:,:)=u_x
        meshList(un_ele)%u_overlap(2,:,:)=u_y
        meshList(un_ele)%u_ele(:,:,:) = 0
        meshList(un_ele)%u_ele(1,:,:) = u_x ! suggest this should be unity
        meshList(un_ele)%u_ele(2,:,:) = u_y ! suggest this should be unity

        allocate(meshList(un_ele)%scaling_var(multi_levels))
        do ilevel=1,multi_levels
          allocate(meshList(un_ele)%scaling_var(ilevel)%detwei(ngi))
          allocate(meshList(un_ele)%scaling_var(ilevel)%sdetwei(sngi,nface))
          allocate(meshList(un_ele)%scaling_var(ilevel)%nx(ngi,ndim,nloc))
        end do


        do iface=1,nface

          call getNeigDataMesh(meshlist, un_ele, iface, Npos, Nside, lnodes2)
          meshList(un_ele)%S_nodes(:,iface)=lnodes2
          meshList(un_ele)%fNeig(iface)=Nside

        end do
      end do

      ! meshList(1)%t_overlap=0.0
      ! meshList(2)%t_overlap=0.0
! just testing git2
      ! t_bc(:)=0.0 ! this contains the boundary conditions just outside the domain
      u_bc(1,:,:)=u_x ! this contains the boundary conditions on velocity just outside the domain
      u_bc(2,:,:)=u_y
      do ilevel=1,multi_levels
        tracer(ilevel)%tnew = 0.0
        tracer(ilevel)%error = 0.0
      end do
      !!!!!!!!!!!!!!!!!!!!!!!!!!! new semi mesh ICs !!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! tnew(:,3,1)=1.0

      do un_ele=1,totele_unst
        do iface=1,nface
          call getNeigDataMesh(meshlist, un_ele, iface, Npos, Nside, lnodes2)
          call get_d_center(MeshList,un_ele, Npos, iface, nface, nloc, n_split,ele_info(1)%str_neig)
        end do
        if ( meshList(un_ele)%region_id == 12 ) then
          tracer(1)%tnew(:,:,un_ele) = 0.
        end if
      end do


      ! call TRIQUAold(L1, L2, L3, L4, WEIGHT, D3,NGI)
        ! call get_surface_ele(meshlist, n_split, totele_unst, nface)
      ! call loc_surf_ele2(surf_ele, n_split)

      call get_shape_funs_spec_ustr(n, nlx, nlx_lxx, nlxx, sn_orig, weight, nlx_nod, &
                  nloc, snloc, sngi, ngi, ndim, nface,n_s_list_no, face_sn, face_sn2, face_snlx, &
                  sweight, npoly, ele_type, totele_unst, face_list_no )

      call get_shape_funs_spec(n, nlx, nlx_lxx, nlxx, weight, nlx_nod, &
                nloc, snloc, sngi, ngi, ndim, nface,n_s_list_no, face_sn_str, face_sn2_str, face_snlx_str, &
                sweight_str, npoly, ele_type, totele_str)!, face_list_no)


      i=1 ! saving x_all structured triangles to be used for VTK generation
      do un_ele=1,totele_unst
        do str_ele=1,totele_str
          call get_splitting(meshList(un_ele)%X, n_split, str_ele, x_loc)
          x_all_str(:,:,i) = x_loc !(ndim,nloc,totele)
          i=i+1

          do iloc=1,nloc
            analytical(iloc,str_ele,un_ele) =  boundary(x_loc(1,iloc),x_loc(2,iloc))
          end do

        end do

        call semi_tri_det_nlx_multigrid(meshlist, multi_levels, un_ele, sngl(meshList(un_ele)%X), n, nlx, weight, ndim,&
                                    nloc, ngi, INV_JAC, n_split)
        call semi_det_snlx_multigrid(meshlist, multi_levels,n_split, nface, ndim, sndim, nloc, sngi,sweight, sn,face_sn,&
                                          face_snlx,un_ele)
        ! call get_area(meshlist,un_ele,n_split)
      end do



      print*, '|   n_split =', n_split
      print*, '|   totele un_ele & totele',  totele_unst, totele_str * totele_unst
      print*, '|   ntime = ', ntime
      print*, '---------------------------------------------------------'
      call CPU_TIME(t1)

      do itime=1,ntime
        ! generating VTK files
        if ( vtk_io <= itime ) then
          totele_str = 4**(n_split)
          ! call get_vtu(x_all_str, tnew, totele_str*totele_unst, nloc, itime, ndim, cell_type, solve_for)
          call get_vtu(x_all_str, tracer(1)%tnew,tracer(1)%error, &
                                                analytical,str_ele, un_ele, totele_str*totele_unst,&
                                                nloc, itime, ndim, cell_type, solve_for)


          ! call get_vtk(x_all_str, tnew, totele_str*totele_unst, nloc, itime,ndim, cell_type)
          vtk_io = vtk_io + vtk_interval
        end if

        tracer(1)%told = tracer(1)%tnew
        tnew_nonlin = tracer(1)%tnew
        ! do its=1,nits
          do multigrid=1,n_multigrid
            do ilevel=1,multi_levels
              i_split = n_split-ilevel+1
              if (allocated(tnew_nonlin)) deallocate(tnew_nonlin)
              allocate(tnew_nonlin(nloc,4**(i_split),totele_unst))
              tracer(ilevel)%tnew = 0.0
              tnew_nonlin = 0.0

              call smoother

              call update_overlaps(meshlist,ele_info(ilevel)%surf_ele, tracer(ilevel)%tnew, tracer(ilevel)%told,&
                          t_bc, i_split, nface,totele_unst, nloc, ele_info(ilevel)%str_neig)

              call restrictor(tracer,totele_unst, i_split, ilevel) 

              call get_residual

            end do ! ilevel multigrids


            do ilevel = multi_levels-1,1,-1

            end do !ilevel
          end do ! solve_its
        print*, 'semi', itime
      end do ! do itime=1,ntime

      call CPU_TIME(t2)
      print*, '----------------------------------------------------------'
      print*, '|        cpu_time for time_loop = ', t2-t1, '     |'
      print*, '----------------------------------------------------------'



      nullify(mesh_pnt%ptr , sdetwei, tnew_nonlin_loc,told_loc2,u_loc2,told_loc,tnew_loc,sn,detwei,u_loc)
      do i=1,multi_levels
        deallocate(tracer(i)%tnew, tracer(i)%told, tracer(i)%RHS, tracer(i)%residuale, tracer(i)%error,tracer(i)%source)
      end do
      deallocate(face_ele, face_list_no, n, nlx, nx, M, nlx_lxx, nlxx, nlx_nod, weight,&
                 sn2,snlx, sweight, s_cont, sweight_str, tnew_nonlin,&
                 t_bc, u_bc, tnew_gi, tnew_xgi, tnew_sgi, tnew_sgi2, told_gi,&
                 face_sn, face_sn2, face_snlx, face_sn_str, face_sn2_str, face_snlx_str, x_loc, ugi,&
                 xsgi, usgi, usgi2, income, snorm, norm, ml_ele, rhs_jac, mass_t_new,&
                 SN_orig,SNLX_orig, inv_jac, mat_diag_approx, mat_tnew,mass_ele_old,stiff_ele_old,stiff_ele_new,&
                 L1, L2, L3, L4, x_all_str, A_x,mass_told,flux_ele_new,mass_ele_new,flux_ele_old,&
                 s_cont_old,diff_vol,diff_surf,fin_ele, jac_res1,&
                 mass_ele,diff_vol2,stiff_ele,mass_stcl,stiff_stcl,diff_vol_stcl,diff_vol1,&
                 my_diff_surf,neig_diff_surf,tnew_loc2,stiff_ele_new1,stiff_ele_old1,tnew_nonlin_loc2)

    contains



      !####################### cal A and b #####################################
      subroutine get_A_x(upto_date)
        ! upto_date defines if we use Gausidel or Jacobi we need Tnew updated or not
        ! if true, Gauss. if false, Jacobi
        logical, INTENT(IN) :: upto_date

        do iloc=1,nloc
          mass_ele_old(iloc) = 1/dt* sum(mass_stcl(iloc,:)*told_loc(:)) !sum(mass_ele(:,iloc)*told_gi(:))
          stiff_ele_old(iloc) = sum(stiff_ele_old1(iloc,:)*told_loc(:))
        end do

        if (upto_date) then
          do iloc=1,nloc
            mass_ele_new(iloc) = 1/dt* sum(mass_stcl(iloc,:)*tnew_nonlin_loc(:)) !sum(mass_ele(:,iloc)*tnew_gi(:))
            diff_vol(iloc) = sum(diff_vol1(iloc,:)*tnew_nonlin_loc(:))
            stiff_ele_new(iloc) = sum(stiff_ele_new1(iloc,:)*tnew_nonlin_loc(:))
            forall ( jloc=1:nloc, iface=1:nface)
              diff_surf(iloc) = diff_surf(iloc) + sum(my_diff_surf(iloc,:,iface) * tnew_nonlin_loc(:))   &
                                                - sum(Neig_diff_surf(iloc,:,iface) * tnew_nonlin_loc2(:,str_ele,un_ele,iface))
            end forall
          end do
        else
         do iloc=1,nloc
           mass_ele_new(iloc) = 1/dt* sum(mass_stcl(iloc,:)*tnew_loc(:)) !sum(mass_ele(:,iloc)*tnew_gi(:))
           diff_vol(iloc) = sum(diff_vol1(iloc,:)*tnew_loc(:))
           stiff_ele_new(iloc) = sum(stiff_ele_new1(iloc,:)*tnew_loc(:))
           forall ( jloc=1:nloc, iface=1:nface)
             diff_surf(iloc) = diff_surf(iloc) + sum(my_diff_surf(iloc,:,iface) * tnew_loc(:))   &
                                               - sum(Neig_diff_surf(iloc,:,iface) * tnew_nonlin_loc2(:,str_ele,un_ele,iface))
           end forall
         end do
       end if
       do iloc=1,nloc
         A_x(iloc) = theta*(mass_ele_new(iloc) - stiff_ele_new(iloc) + flux_ele_new(iloc) &
         + diff_vol(iloc) + diff_surf(iloc))&
         +(1.-theta)*(mass_ele_new(iloc))
       end do
      end subroutine get_A_x



      subroutine get_RHS ! b
        if ( ilevel==1 ) then

       do iloc=1,nloc
           tracer(ilevel)%source(iloc,str_ele,un_ele) = sum(mass_stcl(iloc,:)*tracer(ilevel)%source(:,str_ele,un_ele))

           tracer(ilevel)%RHS(iloc,str_ele,un_ele)=theta*(mass_ele_old(iloc)+tracer(ilevel)%source(iloc,str_ele,un_ele))&
           +(1.-theta)*(mass_ele_old(iloc) + stiff_ele_old(iloc)-flux_ele_old(iloc) &
           - diff_vol(iloc) - diff_surf(iloc) + tracer(ilevel)%source(iloc,str_ele,un_ele) )
         end do
       end if
       ! b(:,str_ele,un_ele) => tracer(ilevel)%RHS(iloc,str_ele,un_ele)
      end subroutine get_RHS



      subroutine get_diff_surf_stencl
        do iloc=1,nloc
          do jloc=1,nloc
            my_diff_surf(iloc,jloc,iface) =  (k/delta_x)*sum(sn(:,iloc) *sn(:,jloc) * sdetwei(:))
            Neig_diff_surf(iloc,jloc,iface) =(k/delta_x)*sum(sn(:,iloc) *sn2(:,jloc)* sdetwei(:))
          end do
        end do ! iloc
      end subroutine get_diff_surf_stencl



      subroutine get_diagonal()
        do iloc=1,nloc
          mat_diag_approx(iloc) = 1/dt* ml_ele(iloc) + diff_vol1(iloc,iloc) &!
          + sum(my_diff_surf(iloc,iloc,:))
        end do
      end subroutine get_diagonal



      !########################## Solvers ######################################
      subroutine solve_Jacobi
        ! D x^{i+1} = D x^i - A x^i +s
        do iloc = 1, nloc
          tnew_nonlin_loc(iloc) = tnew_loc(iloc) +   omga/mat_diag_approx(iloc) *(tracer(ilevel)%RHS(iloc,str_ele,un_ele)&
                                                                                  - A_x(iloc))
        end do
      end subroutine solve_Jacobi



      subroutine solve_Gauss_Seidel
        ! D x^{i+1} = D x^{i+1} - A x^{i+1} +s
        do iloc = 1, nloc
          tnew_nonlin_loc(iloc) = tnew_nonlin_loc(iloc) +   omga/mat_diag_approx(iloc) &
                                                        *(tracer(ilevel)%RHS(iloc,str_ele,un_ele)- A_x(iloc))
        end do
      end subroutine solve_Gauss_Seidel



      subroutine solve_Richardson
        ! implicit none
        ! ! x^{k+1} = x^k + omga(b-Ax^k)
        do iloc=1,nloc
          tnew_nonlin_loc(iloc) = tnew_nonlin_loc(iloc) + &
          omga * (tracer(ilevel)%RHS(iloc,str_ele,un_ele)-(mass_ele_new(iloc) - stiff_ele_new(iloc) + flux_ele_new(iloc)))
        end do
      end subroutine solve_Richardson



      subroutine get_local_residual
        do iloc=1,nloc
          tracer(ilevel)%residuale(iloc,str_ele,un_ele) = A_x(iloc)*tnew_nonlin_loc(iloc)&
                                                                  -tracer(ilevel)%RHS(iloc,str_ele,un_ele)
        end do
      end subroutine get_local_residual



      subroutine get_error
        do iloc=1,nloc
          ! analytical(iloc,str_ele,un_ele) =  boundary(x_loc(1,iloc),x_loc(2,iloc))
          tracer(ilevel)%error(iloc,str_ele,un_ele) = abs(tnew_nonlin_loc(iloc) - analytical(iloc,str_ele,un_ele))
        end do
      end subroutine get_error


      subroutine smoother
        if (allocated(tnew_nonlin_loc2)) deallocate(tnew_nonlin_loc2)
        allocate(tnew_nonlin_loc2(nloc,4**(i_split),totele_unst,nface))


        do smooth =1, n_smooth
          ! tnew = tnew_nonlin ! for non-linearity
          tracer(ilevel)%tnew = tnew_nonlin ! for non-linearity

          ! call update_overlaps3(meshlist, tnew, told, t_bc, n_split, nface, totele_unst, nloc, str_neig)
          ! call update_overlaps2(meshlist, surf_ele2, tnew, told, t_bc, n_split, nface, totele_unst, nloc, str_neig)
          ! call update_overlaps3(meshlist,surf_ele, tnew, told, t_bc, n_split, nface, totele_unst, nloc, str_neig)
          call update_overlaps(meshlist,ele_info(ilevel)%surf_ele, tracer(ilevel)%tnew, tracer(ilevel)%told,&
                            t_bc, i_split, nface,totele_unst, nloc, ele_info(ilevel)%str_neig)

          do un_ele=1,totele_unst
            mesh_pnt%ptr => meshList(un_ele)
            detwei => mesh_pnt%ptr%scaling_var(ilevel)%detwei
            nx = mesh_pnt%ptr%scaling_var(ilevel)%nx
            call get_un_ele_mass_stiff_diffvol(mass_stcl,stiff_stcl,diff_vol_stcl,n,nx,detwei,k,nloc,ngi,ndim,dt,ml_ele)
            totele_str = 4**(i_split)
            do str_ele=1,totele_str
              call semi_get_nx_pos(irow, ipos, orientation, i_split, str_ele, nx, updown)

              ! get_splitting splits an unstructured ele n_split times
              call get_splitting(mesh_pnt%ptr%X, i_split, str_ele, x_loc)

              ! #################################################################
              u_loc => mesh_pnt%ptr%u_ele(:,:,str_ele)
              tnew_loc =>  tracer(ilevel)%tnew(:,str_ele, un_ele)
              told_loc =>  tracer(ilevel)%told(:,str_ele, un_ele)

              do gi=1,ngi
                do idim=1,ndim
                  ugi(gi,idim)=sum(n(gi,:)*u_loc(idim,:))
                  ! tnew_xgi(gi,idim)=sum(nx(gi,idim,:)*tnew_loc(:))
                end do
                ! tnew_gi(gi)=sum(n(gi,:)*tnew_loc(:))
                ! told_gi(gi)=sum(n(gi,:)*told_loc(:))
              end do

              !#################### cal M and K ##################################
              stiff_ele_new = 0.0
              stiff_ele_new1 = 0.0
              stiff_ele_old = 0.0
              stiff_ele_old1 = 0.0
              mass_ele_new = 0.0
              mass_ele_old = 0.0
              diff_vol1=0.0
              do iloc=1,nloc
                tracer(ilevel)%source(iloc,str_ele,un_ele) = -2*k*boundary(x_loc(1,iloc) , x_loc(2,iloc))

                do jloc=1,nloc
                  stiff_ele_new1(iloc,jloc) = stiff_ele_new1(iloc,jloc) + &
                                    (sum(stiff_stcl(jloc,:,1,iloc)*ugi(:,1)) + sum(stiff_stcl(jloc,:,2,iloc)*ugi(:,2)))
                  stiff_ele_old1(iloc,jloc) = stiff_ele_old1(iloc,jloc) + &
                                    (sum(stiff_stcl(jloc,:,1,iloc)*ugi(:,1)) + sum(stiff_stcl(jloc,:,2,iloc)*ugi(:,2)))
                end do

                do idim=1,ndim
                  do jloc=1,nloc
                    diff_vol1(iloc,jloc) = diff_vol1(iloc,jloc) + sum(diff_vol_stcl(:,idim,iloc,jloc))
                  end do
                end do
              end do ! iloc
              stiff_ele_new1 = stiff_ele_new1*updown
              stiff_ele_old1 = stiff_ele_old1*updown

              !############################### surface integral ##################
              flux_ele_new=0.0
              flux_ele_old=0.0
              diff_surf=0.0
              my_diff_surf =0.0
              Neig_diff_surf =0.0
              do iface = 1,nface
                ele22 => ele_info(ilevel)%str_neig(iface,str_ele)!face_ele( iface, ele) ! if ele22 is 0 assume on boundary
                i_got_boundary = (sign(1, -ele22) + 1 )/2
                ele2 = str_ele*i_got_boundary + ele22*(1-i_got_boundary)! if ele22=0, then ele2=ele, if ele22/=0, ele2=ele22

                !######################## vol and surf shape functions ############
                sn => face_sn_str(:,:,iface)
                ! snlx = face_snlx_str(:,:,:,iface)
                sn2=0
                if ( ele22==0 ) then
                      if ( iface==1 ) then
                        sp = ipos/2+1 ! position of ele along the un_iface
                        mface = 1 ! un_iface number which str_ele is located on
                      elseif ( iface==2 ) then
                        sp = irow
                        mface = 3
                      else
                        sp = irow
                        mface = 2
                      end if
                      call get_semi_sn2_implicit(meshlist, un_ele, mface, sn_orig, sn2)
                else
                  sn2  = face_sn2_str(:,:,iface)
                end if ! end if ele22=0

                ! ###################################################################################
                 if ( i_got_boundary < 1e-1 ) then
                   tnew_loc2(:,str_ele,un_ele,iface) = tracer(ilevel)%tnew(:,ele2, un_ele)
                   tnew_nonlin_loc2(:,str_ele,un_ele,iface) = tnew_nonlin(:,ele2, un_ele)
                   told_loc2 => tracer(ilevel)%told(:,ele2, un_ele)
                   u_loc2 => mesh_pnt%ptr%u_ele(:,:,ele2)
                 else
                   tnew_loc2(:,str_ele,un_ele,iface) = mesh_pnt%ptr%t_overlap( (sp)*nloc-2:(sp)*nloc, mface )
                   tnew_nonlin_loc2(:,str_ele,un_ele,iface) = mesh_pnt%ptr%t_overlap( (sp)*nloc-2:(sp)*nloc, mface )
                   told_loc2 => mesh_pnt%ptr%t_overlap_old( (sp)*nloc-2:(sp)*nloc, mface )
                   u_loc2 => mesh_pnt%ptr%u_overlap(:, sp*nloc-2:sp*nloc, mface )
                 end if

                call get_loc_sgi(tnew_sgi, tnew_sgi2, told_sgi, told_sgi2,usgi,usgi2,iface,str_ele,un_ele,&
                                      u_loc, u_loc2, sn, sn2, tnew_loc, tnew_loc2, told_loc, told_loc2, nloc, ndim)

                ! ############################# get flux ############################################
                sdetwei => mesh_pnt%ptr%scaling_var(ilevel)%sdetwei(:,iface)
                snorm = updown * mesh_pnt%ptr%snorm(:,:,iface)
                ! income=1 if info is comming from neighbouring element.

                do s_gi=1,sngi
                  income(s_gi)=0.5 + 0.5*sign(1.0, -sum(snorm(s_gi,:)*0.5*(usgi(s_gi,:)+usgi2(s_gi,:)))  )
                end do

                do idim=1,ndim
                  s_cont(:,idim) = snorm(:,idim)*sdetwei(:) &
                        *( (1.-income(:))*usgi(:,idim)*tnew_sgi(:) + income(:)*tnew_sgi2(:)*usgi2(:,idim) )
                  s_cont_old(:,idim) = snorm(:,idim)*sdetwei(:) &
                        *( (1.-income(:))*usgi(:,idim)*told_sgi(:) + income(:)*told_sgi2(:)*usgi2(:,idim) )
                end do

                do iloc=1,nloc
                  do idim=1,ndim
                    flux_ele_new(iloc)  = flux_ele_new(iloc)  + sum( sn(:,iloc)*s_cont(:,idim) )
                    flux_ele_old(iloc)  = flux_ele_old(iloc)  + sum( sn(:,iloc)*s_cont(:,idim) )
                  end do
                end do

                !################################# end get flux ##########################################
                call add_diffusion_surf(meshlist, un_ele, mface, iface, nface, k, i_split, n, ngi, sngi, ele22,sn2,&
                                      nloc, sdetwei, tnew_sgi,tnew_sgi2, diff_surf, delta_x,sn,my_diff_surf,neig_diff_surf)
                call get_diff_surf_stencl
              end do ! iface

              !################################# solve the system ########################################
              tnew_nonlin_loc => tnew_nonlin(:,str_ele, un_ele)

              select case(solver)
                case(1)
                  call get_A_x(.false.)
                  call get_RHS
                  call get_diagonal
                  call solve_Jacobi
                case(2)
                  call solve_Richardson
                case(3)
                  call get_A_x(.true.)
                  call get_RHS
                  call get_diagonal
                  call solve_Gauss_Seidel
              end select

              call get_error
            end do ! end do totele_str
          end do ! end totele_unst
        end do ! smooth
      end subroutine smoother


      subroutine get_residual
        do un_ele=1,totele_unst
          mesh_pnt%ptr => meshList(un_ele)
          detwei => mesh_pnt%ptr%scaling_var(ilevel)%detwei
          nx = mesh_pnt%ptr%scaling_var(ilevel)%nx
          call get_un_ele_mass_stiff_diffvol(mass_stcl,stiff_stcl,diff_vol_stcl,n,nx,detwei,k,nloc,ngi,ndim,dt,ml_ele)
          totele_str = 4**(n_split-ilevel+1)
          do str_ele=1,totele_str
            call semi_get_nx_pos(irow, ipos, orientation, i_split, str_ele, nx, updown)

            ! get_splitting splits an unstructured ele n_split times
            call get_splitting(mesh_pnt%ptr%X, i_split, str_ele, x_loc)

            ! #################################################################
            u_loc => mesh_pnt%ptr%u_ele(:,:,str_ele)
            tnew_loc =>  tracer(ilevel)%tnew(:,str_ele, un_ele)
            told_loc =>  tracer(ilevel)%told(:,str_ele, un_ele)

            do gi=1,ngi
              do idim=1,ndim
                ugi(gi,idim)=sum(n(gi,:)*u_loc(idim,:))
                ! tnew_xgi(gi,idim)=sum(nx(gi,idim,:)*tnew_loc(:))
              end do
              ! tnew_gi(gi)=sum(n(gi,:)*tnew_loc(:))
              ! told_gi(gi)=sum(n(gi,:)*told_loc(:))
            end do

            !#################### cal M and K ##################################
            stiff_ele_new = 0.0
            stiff_ele_new1 = 0.0
            stiff_ele_old = 0.0
            stiff_ele_old1 = 0.0
            mass_ele_new = 0.0
            mass_ele_old = 0.0
            diff_vol1=0.0
            do iloc=1,nloc
              tracer(ilevel)%source(iloc,str_ele,un_ele) = -2*k*boundary(x_loc(1,iloc) , x_loc(2,iloc))

              do jloc=1,nloc
                stiff_ele_new1(iloc,jloc) = stiff_ele_new1(iloc,jloc) + &
                                  (sum(stiff_stcl(jloc,:,1,iloc)*ugi(:,1)) + sum(stiff_stcl(jloc,:,2,iloc)*ugi(:,2)))
                stiff_ele_old1(iloc,jloc) = stiff_ele_old1(iloc,jloc) + &
                                  (sum(stiff_stcl(jloc,:,1,iloc)*ugi(:,1)) + sum(stiff_stcl(jloc,:,2,iloc)*ugi(:,2)))
              end do

              do idim=1,ndim
                do jloc=1,nloc
                  diff_vol1(iloc,jloc) = diff_vol1(iloc,jloc) + sum(diff_vol_stcl(:,idim,iloc,jloc))
                end do
              end do
            end do ! iloc
            stiff_ele_new1 = stiff_ele_new1*updown
            stiff_ele_old1 = stiff_ele_old1*updown

            !############################### surface integral ##################
            flux_ele_new=0.0
            flux_ele_old=0.0
            diff_surf=0.0
            my_diff_surf =0.0
            Neig_diff_surf =0.0
            do iface = 1,nface
              ele22 => ele_info(ilevel)%str_neig(iface,str_ele)!face_ele( iface, ele) ! if ele22 is 0 assume on boundary
              i_got_boundary = (sign(1, -ele22) + 1 )/2
              ele2 = str_ele*i_got_boundary + ele22*(1-i_got_boundary)! if ele22=0, then ele2=ele, if ele22/=0, ele2=ele22

              !######################## vol and surf shape functions ############
              sn => face_sn_str(:,:,iface)
              ! snlx = face_snlx_str(:,:,:,iface)
              sn2=0
              if ( ele22==0 ) then
                    if ( iface==1 ) then
                      sp = ipos/2+1 ! position of ele along the un_iface
                      mface = 1 ! un_iface number which str_ele is located on
                    elseif ( iface==2 ) then
                      sp = irow
                      mface = 3
                    else
                      sp = irow
                      mface = 2
                    end if
                    call get_semi_sn2_implicit(meshlist, un_ele, mface, sn_orig, sn2)
              else
                sn2  = face_sn2_str(:,:,iface)
              end if ! end if ele22=0

              ! ###################################################################################
               if ( i_got_boundary < 1e-1 ) then
                 tnew_loc2(:,str_ele,un_ele,iface) = tracer(ilevel)%tnew(:,ele2, un_ele)
                 tnew_nonlin_loc2(:,str_ele,un_ele,iface) = tnew_nonlin(:,ele2, un_ele)
                 told_loc2 => tracer(ilevel)%told(:,ele2, un_ele)
                 u_loc2 => mesh_pnt%ptr%u_ele(:,:,ele2)
               else
                 tnew_loc2(:,str_ele,un_ele,iface) = mesh_pnt%ptr%t_overlap( (sp)*nloc-2:(sp)*nloc, mface )
                 tnew_nonlin_loc2(:,str_ele,un_ele,iface) = mesh_pnt%ptr%t_overlap( (sp)*nloc-2:(sp)*nloc, mface )
                 told_loc2 => mesh_pnt%ptr%t_overlap_old( (sp)*nloc-2:(sp)*nloc, mface )
                 u_loc2 => mesh_pnt%ptr%u_overlap(:, sp*nloc-2:sp*nloc, mface )
               end if

              call get_loc_sgi(tnew_sgi, tnew_sgi2, told_sgi, told_sgi2,usgi,usgi2,iface,str_ele,un_ele,&
                                    u_loc, u_loc2, sn, sn2, tnew_loc, tnew_loc2, told_loc, told_loc2, nloc, ndim)

              ! ############################# get flux ############################################
              sdetwei => mesh_pnt%ptr%scaling_var(ilevel)%sdetwei(:,iface)
              snorm = updown * mesh_pnt%ptr%snorm(:,:,iface)
              ! income=1 if info is comming from neighbouring element.

              do s_gi=1,sngi
                income(s_gi)=0.5 + 0.5*sign(1.0, -sum(snorm(s_gi,:)*0.5*(usgi(s_gi,:)+usgi2(s_gi,:)))  )
              end do

              do idim=1,ndim
                s_cont(:,idim) = snorm(:,idim)*sdetwei(:) &
                      *( (1.-income(:))*usgi(:,idim)*tnew_sgi(:) + income(:)*tnew_sgi2(:)*usgi2(:,idim) )
                s_cont_old(:,idim) = snorm(:,idim)*sdetwei(:) &
                      *( (1.-income(:))*usgi(:,idim)*told_sgi(:) + income(:)*told_sgi2(:)*usgi2(:,idim) )
              end do

              do iloc=1,nloc
                do idim=1,ndim
                  flux_ele_new(iloc)  = flux_ele_new(iloc)  + sum( sn(:,iloc)*s_cont(:,idim) )
                  flux_ele_old(iloc)  = flux_ele_old(iloc)  + sum( sn(:,iloc)*s_cont(:,idim) )
                end do
              end do

              !################################# end get flux ##########################################
              call add_diffusion_surf(meshlist, un_ele, mface, iface, nface, k, i_split, n, ngi, sngi, ele22,sn2,&
                                    nloc, sdetwei, tnew_sgi,tnew_sgi2, diff_surf, delta_x,sn,my_diff_surf,neig_diff_surf)
              call get_diff_surf_stencl
            end do ! iface

            !################################# solve the system ########################################
            tnew_nonlin_loc => tnew_nonlin(:,str_ele, un_ele)
            call get_A_x(.false.)
            call get_RHS
            tracer(ilevel)%residuale(:,str_ele, un_ele) = A_x(:) -  tracer(ilevel)%RHS(:,str_ele,un_ele)
          end do
        end do

      end subroutine get_residual

    end subroutine Semi_implicit_iterative
















    Subroutine Semi_implicit_iterative_P(CFL, ndim, direct_solver, nsolver_its, time, nits,&
          nface, dx, dy, nloc, snloc, ngi, sngi, n_s_list_no, u_x, u_y, with_stab, vtk_interval, cell_type)
          implicit none

          ! external variables
          real, intent(in) :: CFL, time, u_x, u_y, dx, dy
          integer, intent(in) :: ndim, nits, nface, nloc, ngi, sngi, snloc, nsolver_its
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
          integer :: errorflag, i_got_boundary, solver_its, its
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
          real :: mass_ele, stiff_ele, lhs, flux_ele,omega

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
          real, allocatable :: mat_loc(:,:), mat_loc_inv(:,:), ml_ele(:),  rhs_jac(:), mass_t_new(:)
          real, allocatable :: SN_orig(:,:),SNLX_orig(:,:,:), inv_jac(:, :, : ), mat_diag_approx(:)
          real, allocatable :: a_star(:,:), p_star(:), rgi(:), diff_coe(:), told_loc(:), stab(:,:), mat_tnew(:), mass_told(:)
          real, allocatable :: L1(:), L2(:), L3(:), L4(:), x_all_str(:,:,:), glob_told_array(:)
          real, allocatable :: glob_lhs(:,:), inv_lhs(:,:), rhs_array(:), glob_flux_ele(:,:)
          real, allocatable :: tnew_nonlin2(:), stiff_ele2(:), mass_ele2(:,:),b(:),A(:,:)
          real, allocatable :: tnew_nonlin_sgi(:), tnew_nonlin_sgi2(:), flux_ele2(:),tnew_nonlin_gi(:)

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
          omega=.5
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
          allocate(mat_loc_inv(nloc,nloc), ml_ele(nloc))
          allocate(rhs_jac(nloc), mass_t_new(nloc), inv_jac(ngi, ndim, ndim ) )
          allocate(a_star(ngi,nloc), p_star(ngi), rgi(ngi), diff_coe(ngi), told_loc(nloc), stab(nloc, nloc))
          allocate(mat_tnew(nloc), mat_diag_approx(nloc), mass_told(nloc))
          allocate(L1(ngi), L2(ngi), L3(ngi), L4(ngi), x_all_str(ndim,nloc,totele_unst*totele_str))
          allocate(inv_lhs(totele_str*totele_unst*nloc, totele_str*totele_unst*nloc))
          allocate(rhs_array(totele_str*totele_unst*nloc), tnew_nonlin2(totele_str*totele_unst*nloc))
          allocate(glob_lhs(totele_str*totele_unst*nloc,totele_str*totele_unst*nloc))
          allocate(glob_flux_ele(nloc,totele_str*totele_unst*nloc), stiff_ele2(nloc))
          allocate(mass_ele2(nloc,nloc), flux_ele2(nloc))
          allocate(b(totele_str*totele_unst*nloc), A(nloc,totele_str*totele_unst*nloc))
          allocate(tnew_nonlin_sgi(sngi),tnew_nonlin_sgi2(sngi),tnew_nonlin_gi(sngi))

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
          tnew(:,6:7,1)=1.0
          ! do un_ele=1,totele_unst
          !     if ( meshList(un_ele)%region_id == 13 ) then
          !       tnew(:,:,un_ele) = 1.
          !     end if
          !   end do

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
              do solver_its=1,nsolver_its ! jacobi iterations...
              tnew = tnew_nonlin ! for non-linearity
              call update_overlaps3(meshlist, tnew, told, t_bc, n_split, nface, totele_unst, nloc, str_neig)
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
                    !################### Petrov coef calculations ######################
                    do gi=1,ngi
                      do idim=1,ndim
                        ugi(gi,idim)=sum(n(gi,:)*u_loc(idim,:))
                        tnew_xgi(gi,idim)=sum(nx(gi,idim,:)*tnew_loc(:))
                      end do
                      tnew_gi(gi)=sum(n(gi,:)*tnew_loc(:))
                      told_gi(gi)=sum(n(gi,:)*told_loc(:))
                      ! tnew_nonlin_gi(gi) = sum(n(gi,:)*tnew_nonlin_loc(:))
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

                    !#################### cal M and K ##################################
                    stab=0.0
                    stiff_ele2=0.0
                    mass_ele2=0.0
                    do iloc=1,nloc
                      g_iloc = glob_no_semi(un_ele, str_ele, n_split, nloc, iloc)
                      do jloc=1,nloc
                        g_jloc = glob_no_semi(un_ele, str_ele, n_split, nloc, jloc)
                        mass_ele = 0.0
                        stiff_ele = 0.0
                        do gi=1,ngi
                          ! stab(iloc,jloc) = stab(iloc,jloc) + sum(diff_coe(gi)* nx(gi,:,jloc)* nx(gi,:,iloc))* detwei(gi)
                          mass_ele = mass_ele + n(gi,iloc)*n(gi,jloc)*detwei(gi)/dt
                          do idim=1,ndim
                            stiff_ele = stiff_ele + nx(gi,idim,iloc)*ugi(gi,idim)*detwei(gi)*n(gi,jloc)
                          end do
                        end do ! quadrature
                        mass_ele2(iloc,jloc)=mass_ele
                        call add_to_CSR(sparse_mass, g_iloc, g_jloc, mass_ele)

                        lhs = mass_ele-stiff_ele
                        call add_to_CSR_flux(sparse_lhs, g_iloc, g_jloc, lhs)
                      end do ! jloc
                      ml_ele(iloc)=sum(n(:,iloc)*detwei(:))/dt ! lumped mass matrix in element (use later for solver).

                      do gi=1,ngi
                        do idim=1,ndim
                          stiff_ele2(iloc) = stiff_ele2(iloc) + nx(gi,idim,iloc)*tnew_gi(gi)*ugi(gi,idim)*detwei(gi)
                        end do
                      end do ! quadrature
                    end do ! iloc


                    !############################### surface integral ##################
                    flux_ele=0.0
                    flux_ele2=0.0
                    do iface = 1,nface
                      ele22 = str_neig(iface,str_ele)!face_ele( iface, ele) ! if ele22 is 0 assume on boundary
                      i_got_boundary = (sign(1, -ele22) + 1 )/2
                      r_got_boundary = real(i_got_boundary)
                      ele2 = str_ele*i_got_boundary + ele22*(1-i_got_boundary)! if ele22=0, then ele2=ele, if ele22/=0, ele2=ele22

                      !######################## vol and surf shape functions ############
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


                      tnew_loc2(:)= tnew(:,ele2, un_ele) * (1.0-r_got_boundary) &
                                  + meshlist(un_ele)%t_overlap( (sp)*nloc-2:(sp)*nloc, mface ) * r_got_boundary

                      u_loc2(:,:)= meshList(un_ele)%u_ele(:,:,ele2)* (1.0-r_got_boundary) &
                                  + meshlist(un_ele)%u_overlap(:, sp*nloc-2:sp*nloc, mface ) * r_got_boundary
                                  ! 2nd section of u_overlap should refer to 3 ilocs from the un_ele u_bc which can be identified based on irow

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


                      !############### Finding neighbouring element ##################
                      target_ele = int( (1.-income(1))*str_ele + income(1)*ele22 )
                      Npos = meshList(un_ele)%Neig(mface)
                      Nside = meshList(un_ele)%fNeig(mface)

                      if ( target_ele==0 .and. Npos==0 ) then ! it happens when ele is on the boundary
                        tnew_nonlin_sgi=0.0
                        tnew_nonlin_sgi2=0.0
                        do iloc=1,nloc
                          tnew_nonlin_sgi = tnew_nonlin_sgi(:) + sn(:,iloc)*tnew_nonlin(iloc,str_ele,un_ele)
                          tnew_nonlin_sgi2 = tnew_nonlin_sgi
                        end do
                      elseif ( target_ele==0 .and. Npos/=0) then
                        if ( iface==1 ) then
                          ipos = int(ipos/2)+1
                          target_ele = meshList(un_ele)%s_ele(ipos, 1)
                        elseif ( iface==2 ) then
                          target_ele = meshList(un_ele)%s_ele(irow, 3)
                        elseif ( iface==3 ) then
                          target_ele = meshList(un_ele)%s_ele(irow, 2)
                        end if
                        tnew_nonlin_sgi=0.0
                        tnew_nonlin_sgi2=0.0
                        do iloc=1,nloc
                          tnew_nonlin_sgi = tnew_nonlin_sgi(:) + sn(:,iloc)*tnew_nonlin(iloc,str_ele,un_ele)
                          tnew_nonlin_sgi2 = tnew_nonlin_sgi2(:) + sn2(:,iloc)*tnew_nonlin(iloc,target_ele,Npos)
                        end do
                      elseif ( target_ele /= 0) then
                        target_ele = (un_ele-1)*4**n_split + target_ele
                        tnew_nonlin_sgi=0.0
                        tnew_nonlin_sgi2=0.0
                        do iloc=1,nloc
                          tnew_nonlin_sgi = tnew_nonlin_sgi(:) + sn(:,iloc)*tnew_nonlin(iloc,str_ele,un_ele)
                          tnew_nonlin_sgi2 = tnew_nonlin_sgi2(:) + sn2(:,iloc)*tnew_nonlin(iloc,target_ele,un_ele)
                        end do
                      end if

                      !################## calculate flux ##############################
                      do idim=1,ndim
                        s_cont(:,idim) = snorm(:,idim)*sdetwei(:) &
                        *( (1.-income(:))* usgi(:,idim)*tnew_sgi(:) + income(:)*usgi2(:,idim)*tnew_sgi2(:))
                      end do

                      do iloc=1,nloc
                        do idim=1,ndim
                          flux_ele2(iloc)  =  flux_ele2(iloc) + sum( sn(:,iloc)*s_cont(:,idim) )
                        end do
                      end do

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

                  !######################## Solver ################################
                  tnew_nonlin_loc(:) = tnew_nonlin(:,str_ele, un_ele)

                  do iloc=1,nloc
                    tnew_nonlin_loc(iloc) = tnew_nonlin_loc(iloc) + 1. * (sum(mass_ele2(iloc,:)*told_loc(:)) &
                                                      -(sum(mass_ele2(iloc,:)*tnew_nonlin_loc(:))-stiff_ele2(iloc)&
                                                      + flux_ele2(iloc)))
                  end do

                  tnew_nonlin(:,str_ele, un_ele) = tnew_nonlin_loc(:)
                    !#############################################################

                end do ! end do totele_str
              end do ! end totele_unst
              ! sparse_lhs%val = sparse_lhs%val + sparse_flux%val
              ! call gen_global_matrix(totele_str*totele_unst, nloc, sparse_lhs, glob_lhs)
              ! call told_to_array_semi(totele_str, totele_unst, nloc, told, glob_told_array)
              ! call csr_mul_array(sparse_mass, glob_told_array, rhs_array)
              ! tnew_nonlin2=0.0
              ! do solver_its=1,nsolver_its
              !     tnew_nonlin2 = tnew_nonlin2 + omega * (rhs_array - matmul(glob_lhs, tnew_nonlin2))
              ! end do

              ! call FINDInv(glob_lhs, inv_lhs, nloc*totele_str*totele_unst, errorflag)
              !  tnew_nonlin2 = matmul(inv_lhs,rhs_array )
              ! call array_to_told_semi(totele_str, totele_unst, nloc, tnew_nonlin2, tnew_nonlin)
              tnew=tnew_nonlin
            end do ! Gauss
            end do ! do its=1,nits
            print*, 'semi_P', itime
          end do ! do itime=1,ntime

          deallocate(face_ele, face_list_no, str_neig, n, nlx, nx, M, nlx_lxx, nlxx, nlx_nod, weight, detwei, sdetwei,&
                     sn,sn2,snlx, sweight, s_cont, sweight_str, tnew_loc, tnew_loc2, told, tnew, tnew_nonlin,&
                     tnew_nonlin_loc, t_bc, u_bc, tnew_gi, tnew_xgi, tnew_sgi, tnew_sgi2, told_gi,&
                     face_sn, face_sn2, face_snlx, face_sn_str, face_sn2_str, face_snlx_str, u_loc, u_loc2, x_loc, ugi,&
                     xsgi, usgi, usgi2, income, snorm, norm, mat_loc_inv, ml_ele, rhs_jac, mass_t_new,&
                     SN_orig,SNLX_orig, inv_jac, mat_diag_approx, a_star, p_star, rgi, diff_coe, told_loc, stab, mat_tnew,&
                     L1, L2, L3, L4, x_all_str, inv_lhs, rhs_array,glob_lhs,b,A,mass_told,tnew_nonlin_gi,&
                     tnew_nonlin2, sparse_flux%g_iloc, sparse_flux%g_jloc, sparse_flux%val,stiff_ele2,flux_ele2,mass_ele2,&
                     sparse_lhs%g_iloc, sparse_lhs%g_jloc, sparse_lhs%val, sparse_mass%g_iloc, sparse_mass%g_jloc,sparse_mass%val)
        end subroutine Semi_implicit_iterative_P










  ! brief this is for semi-str complete implicit method
  Subroutine Semi_implicit_direct(CFL, ndim, direct_solver, nsolver_its, time, nits,&
    nface, dx, dy, nloc, snloc, ngi, sngi, n_s_list_no, u_x, u_y, with_stab, vtk_interval, cell_type)
    implicit none

    ! external variables
    real, intent(in) :: CFL, time, u_x, u_y, dx, dy
    integer, intent(in) :: ndim, nits, nface, nloc, ngi, sngi, snloc, nsolver_its
    integer, intent(in) :: vtk_interval, cell_type, n_s_list_no
    logical, intent(in) :: direct_solver, with_stab

    ! local vbls
    type(Mesh), allocatable, dimension(:) :: meshList
    type(sparse) :: sparse_mass, sparse_lhs, sparse_flux

    character (len=20):: file_name
    logical :: with_time_slab, D3, diffusion

    integer :: gi, iface, ele, ele2, ele22, n_split, irow, ipos, sum_up, irow_ups
    integer :: totele_unst, un_ele, str_ele, totele_str
    integer :: s_list_no, s_gi, iloc, jloc, i,j, ntime, vtk_io
    integer :: itime, idim, sndim, nonodes, mloc, col!, row, row2
    integer :: errorflag, i_got_boundary, solver_its, its
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
    real :: mass_ele, stiff_ele, flux_ele, lhs
    real :: k

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
    real, allocatable :: tnew_nonlin2(:),diff_vol(:)

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

    n_split =1
    totele_str = 4**n_split
    diffusion = .true.
    k = 2*19**(-5) !diffusion coeficient for the diffusion term, m^2/s for water diffuses into air

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
    allocate(glob_lhs(totele_str*totele_unst*nloc,totele_str*totele_unst*nloc), diff_vol(nloc))

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
    ! tnew(:,3,1)=1.0

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
    call make_sparse_matrix_flux_semi(sparse_mass,meshList,n_split,totele_unst,totele_str,nloc,nface,str_neig)
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


    !#################### Time loop ##########################################
    do itime=1,ntime
      ! generating VTK files
      if ( vtk_io <= itime ) then
        call get_vtk(x_all_str, tnew, totele_str*totele_unst, nloc, itime,ndim, cell_type)
        vtk_io = vtk_io + vtk_interval
      end if

      told = tnew
      tnew_nonlin = tnew
      do its=1,nits
        tnew = tnew_nonlin ! for non-linearity
        call update_overlaps3(meshlist, tnew, told, t_bc, n_split, nface, totele_unst, nloc, str_neig)
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
                call add_to_CSR_flux(sparse_mass, g_iloc, g_jloc, mass_ele)
                lhs = mass_ele-stiff_ele
                call add_to_CSR_flux(sparse_lhs, g_iloc, g_jloc, lhs)
              end do ! jloc
            end do ! iloc

            if ( diffusion ) call add_diffusion_vol(k, nx, ngi, ndim, nloc, detwei, tnew_xgi, diff_vol)



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
                            *( (1.-income(:))* usgi(:,idim) + income(:)*usgi2(:,idim) )
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
! if ( itime==1 .and. its==1 ) then
  ! print*, 'sparse_flux%val'
  ! print*, sparse_flux%val
!   print*, 'mass_glob'
!   call print_my_matrix(glob_mass, totele_str*totele_unst*nloc)
!   print*, ''
!   print*, 'lhs_glob'
!   call print_my_matrix(glob_lhs, totele_str*totele_unst*nloc)
!   print*, ''
!   print*, 'flux_glob'
!   call print_my_matrix(glob_flux, totele_str*totele_unst*nloc)
!   print*, ''
! end if
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
               SN_orig,SNLX_orig, inv_jac, mat_diag_approx, a_star, p_star, rgi, diff_coe, told_loc, stab, mat_tnew,&
               L1, L2, L3, L4, x_all_str, inv_lhs, rhs_array, glob_flux,glob_mass,glob_told_array,glob_lhs,&
               tnew_nonlin2, sparse_flux%g_iloc, sparse_flux%g_jloc, sparse_flux%val,mass_told,&
               sparse_lhs%g_iloc, sparse_lhs%g_jloc, sparse_lhs%val, sparse_mass%g_iloc, sparse_mass%g_jloc,sparse_mass%val)
  end subroutine Semi_implicit_direct


  ! brief> this is for semi-structured explicit method
  Subroutine semi_explicit(CFL, ndim, direct_solver, nsolver_its, time, nits,&
    nface, dx, dy, nloc, snloc, ngi, sngi, n_s_list_no, u_x, u_y, with_stab, vtk_interval, cell_type)
    implicit none

    ! external variables
    real, intent(in) :: CFL, time, u_x, u_y, dx, dy
    integer, intent(in) :: ndim, nits, nface, nloc, ngi, sngi, snloc, nsolver_its
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
    integer :: errorflag, i_got_boundary, solver_its, its
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


    ! call ReadMSH(meshList,'./semi_mesh.msh',ierr, totnodes)
    ! call ReadMSH(meshList,'./irregular shape.msh',ierr, totnodes)
    ! call ReadMSH(meshList,'./new_semi_mesh.msh',ierr, totnodes)
    ! call ReadMSH(meshList,'./gmsh_100.msh',ierr, totnodes)
    print*, 'starting reading the .msh file'
    call CPU_TIME(t1_ReadMSH)
    ! call ReadMSH(meshList,'./10_elements.msh',ierr, totnodes)
    ! call ReadMSH(meshList,'./1_unele_test.msh',ierr, totnodes)
    ! call ReadMSH(meshList,'./2_unele_test.msh',ierr, totnodes)
    ! call ReadMSH(meshList,'./teseting_mesh.msh',ierr, totnodes)
    call ReadMSH(meshList,'./irregular.msh',ierr, totnodes)


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
    allocate(mass_ele(nloc,nloc), mat_loc_inv(nloc,nloc), rhs_loc(nloc), ml_ele(nloc))
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

    !!!!!!!!!!!!!!!!!!!!!!!!!!! new semi mesh ICs !!!!!!!!!!!!!!!!!!!!!!!!!!!!
    do un_ele=1,totele_unst
        if ( meshList(un_ele)%region_id == 13 ) then
          tnew(1,:,un_ele) = 1.
        end if
      end do
  ! tnew(:,7,1)=1.


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
        call update_overlaps3(meshlist, tnew, told, t_bc, n_split, nface, totele_unst, nloc, str_neig)
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

            !################### Petrov coef calculations ######################
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
! if ( itime==1 .and. its==1 ) then
!   print*, tnew_loc(:)
!   print*, tnew_xgi(1,:)
!   print*, ''
! end if
            !#################### cal M and K ##################################
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
! if ( itime==1 .and.its==1 ) then
! print 1, un_ele, str_ele
! 1 format(i2,4X,i2)
! call print_my_matrix(mass_ele, nloc)
! print*, ''
! end if
            !############################### surface integral ##################
            do iface = 1,nface
              ele22 = str_neig(iface,str_ele)!face_ele( iface, ele) ! if ele22 is 0 assume on boundary
              i_got_boundary = (sign(1, -ele22) + 1 )/2
              r_got_boundary = real(i_got_boundary)
              ele2 = str_ele*i_got_boundary + ele22*(1-i_got_boundary)! if ele22=0, then ele2=ele, if ele22/=0, ele2=ele22

              !######################## vol and surf shape functions ############
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

              ! ########### T and u @ quadrature points ########################
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

              !################## calculate flux ##############################
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
                            *( (1.-income(:))* usgi(:,idim)*tnew_sgi(:) + income(:)*usgi2(:,idim)*tnew_sgi2(:) )
              end do

              !######################## cal (K-F)Tnew ##########################
              do iloc=1,nloc
                do idim=1,ndim
                  rhs_loc(iloc)  = rhs_loc(iloc)  - sum( sn(:,iloc)*s_cont(:,idim) )
                end do
              end do
            end do ! iface

            !################## add Petrov if we have it on ####################
            ! if (with_stab) then
              mat_loc= mass_ele + dt*stab
            ! else
            !   stab=0
            !   mat_loc= mass_ele
            ! end if

            !##################### solve the system #############################
            ! if(direct_solver) then
              ! inverse of the mass matric (nloc,nloc)
              ! call FINDInv(mass_ele, mat_loc_inv, nloc, errorflag)
              !can be only this part instead
              ! do iloc=1,nloc
              !    tnew_nonlin(iloc,ele)=told_loc(iloc) + dt* sum( mat_loc_inv(iloc,:)* rhs_loc(:) ) &
              !                                  + told_loc(iloc)* dt* sum( mat_loc_inv(iloc,:) * stab(iloc,:) )
              ! end do
              ! do iloc=1,nloc
              !   mass_told(iloc)=sum( mat_loc(iloc,:) * told_loc(:) )
              ! end do
              ! do iloc=1,nloc
              !   tnew_nonlin(iloc,str_ele, un_ele)=sum( mat_loc_inv(iloc,:)* (mass_told(:)+dt*rhs_loc(:)))
              ! end do

            ! else ! iterative solver
               do iloc=1,nloc
                  rhs_jac(iloc)= sum( mass_ele(iloc,:) * told_loc(:) ) + dt*rhs_loc(iloc)
                  mat_diag_approx(iloc) = ml_ele(iloc) + dt * stab(iloc,iloc)
               end do

               tnew_nonlin_loc(:) = tnew_nonlin(:,str_ele, un_ele)
               do solver_its=1,nsolver_its ! jacobi iterations...
                  do iloc=1,nloc
                     mat_tnew(iloc)= sum( mat_loc(iloc,:) * tnew_nonlin_loc(:) )
                  end do
                  do iloc=1,nloc
                     tnew_nonlin_loc(iloc)=  (mat_diag_approx(iloc)*tnew_nonlin_loc(iloc) - mat_tnew(iloc) + rhs_jac(iloc) )&
                                                                    / mat_diag_approx(iloc)
                  end do
               end do
               tnew_nonlin(:,str_ele, un_ele) = tnew_nonlin_loc(:)
            ! endif ! endof if(direct_solver) then else
            !###################################################################


          end do ! end do totele_str
        end do ! end totele_unst
        tnew=tnew_nonlin
      end do ! do its=1,nits
      print*, 'semi-explicit ', itime
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
