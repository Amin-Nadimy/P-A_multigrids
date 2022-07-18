!>@Brief
! cell_type is defined by vtk file formates, found in: (triangle =5)
! https://kitware.github.io/vtk-examples/site/VTKFileFormats/

module amin
  use ShapFun
  use ShapFun_unstruc
  ! use structured_meshgen
  use matrices
  use Msh2Tri
  use splitting
  use get_vtk_files
  contains

  Subroutine diffusion(CFL, ndim, direct_solver, njac_its, time, nits,&
    nface, dx, dy, nloc, snloc, ngi, sngi, n_s_list_no, u_x, u_y, with_stab, vtk_interval, cell_type, diff_coef)
    implicit none

    ! income variables
    real, intent(in):: CFL, dx, dy, time, u_x, u_y, diff_coef
    integer, intent(in) :: ndim, nits, njac_its, nface, nloc, ngi, sngi, snloc, n_s_list_no
    integer, intent(in):: vtk_interval, cell_type
    logical, INTENT(IN) :: direct_solver, with_stab


    ! local vbl
    character (len=20):: file_name
    integer:: vtk_io, g_iloc, g_jloc, target_ele, target_iface

    integer :: gi, iface, totele, ele, ele2, ele22, n_split
    integer :: s_list_no, s_gi, iloc, jloc, i, j, ntime
    integer:: itime, idim, sndim, nonodes, mloc, col, row, row2
    integer :: errorflag, i_got_boundary, jac_its, its
    integer:: npoly, ele_type, ierr
    integer:: Mpos, Npos, Mside, Nside, lnodes(2), lnodes2(2), totnodes

    real :: sarea, volume, dt, L, r_got_boundary, toler, a_coef
    real:: length_x, length_y, mass_ele, stiff_ele, lhs, diff_ele
    real :: t1, t2, flux_ele, bc_diff
    real :: t1_ReadMSH, t2_ReadMSH, time_ReadMSH
    character(10) :: start, finish
    logical:: with_time_slab, D3

    type(sparse) :: sparse_mass, sparse_lhs, sparse_flux
    type(Mesh), allocatable, dimension(:) :: meshList
    type(neig_data), allocatable, dimension(:) :: neig_info

    character (len=10), allocatable:: str_coor(:,:,:)
    integer, allocatable :: face_ele(:,:), face_list_no(:,:)
    real, allocatable :: n(:,:), nlx(:,:,:), nx(:,:,:), M(:,:)
    real, allocatable :: nlx_lxx(:,:,:), nlxx(:,:), nlx_nod(:,:,:)
    real, allocatable :: weight(:), detwei(:), sdetwei(:)
    real, allocatable :: sn(:,:),sn2(:,:),snlx(:,:,:), sweight(:), s_cont(:,:)
    real, allocatable :: tnew_loc(:), tnew_loc2(:), told(:,:), tnew(:,:), tnew_nonlin(:,:)
    real, allocatable :: tnew_nonlin_loc(:), tnew_nonlin2(:)
    real, allocatable :: t_bc(:,:), u_bc(:,:,:)
    real, allocatable :: tnew_gi(:), tnew_xgi(:,:), tnew_sgi(:), tnew_sgi2(:), told_gi(:)
    real, allocatable :: face_sn(:,:,:), face_sn2(:,:,:), face_snlx(:,:,:,:)
    real, allocatable :: u_loc(:,:), u_ele(:,:,:), u_loc2(:,:), x_loc(:,:), ugi(:,:)
    real, allocatable :: x_all(:,:,:)
    real, allocatable :: xsgi(:,:), usgi(:,:), usgi2(:,:), income(:), snorm(:,:), norm(:)
    real, allocatable :: mat_loc(:,:), mat_loc_inv(:,:), rhs_loc(:), ml_ele(:), rhs_jac(:), mass_t_new(:)
    real, allocatable :: SN_orig(:,:),SNLX_orig(:,:,:), inv_jac(:, :, : ), mat_diag_approx(:)
    real, allocatable :: a_star(:,:), p_star(:), rgi(:), diff_coe(:), told_loc(:), stab(:,:), mat_tnew(:), mass_told(:)
    real, allocatable :: rhs_array(:), glob_lhs_matrix(:,:), glob_mass_matrix(:,:), glob_flux_matrix(:,:)
    real, allocatable :: glob_told_array(:), inv_lhs(:,:), glob_flux(:,:)
    real, allocatable:: L1(:), L2(:), L3(:), L4(:)
    real, allocatable:: NLX_ALL(:,:,:)

    ! xp(2): coordinates of nodes
    ! call ReadMSH(meshList,'./unstr_950272.msh',ierr, totnodes)
    ! call ReadMSH(meshList,'./semi_mesh.msh',ierr, totnodes)
    ! call ReadMSH(meshList,'./900_ele.msh',ierr, totnodes)
print*, 'starting reading the .msh file'
call CPU_TIME(t1_ReadMSH)
! call ReadMSH(meshList,'./10_elements.msh',ierr, totnodes)
! call ReadMSH(meshList,'./gmsh_100.msh',ierr, totnodes)
call ReadMSH(meshList,'./4_elements.msh',ierr, totnodes)
call CPU_TIME(t2_ReadMSH)
time_ReadMSH = t2_ReadMSH - t1_ReadMSH
print*, 'finished reading .msh file'

    totele = size(meshList)

    sndim = ndim-1 ! dimensions of the surface element coordinates.
    mloc=1
    toler = 0.00000000001
    ele_type = 3
    dt = CFL*dx
    length_x = totele/2
    vtk_io=vtk_interval

    allocate(n( ngi, nloc ), nx( ngi, ndim, nloc ), nlx( ngi, ndim, nloc ), M(ngi,nloc))
    allocate(nlx_lxx(ngi,ndim,nloc), nlxx(ngi,nloc), nlx_nod(nloc,ndim,nloc))
    allocate(weight(ngi), detwei(ngi), sdetwei(sngi))
    allocate(face_list_no( nface, totele), face_ele(nface,totele))
    allocate(sn(sngi,nloc),sn2(sngi,nloc),snlx(sngi,sndim,nloc))
    allocate(SN_orig(sngi,snloc),SNLX_orig(sngi,sndim,snloc), s_cont(sngi,ndim) )
    allocate(u_loc(ndim,nloc), u_loc2(ndim,nloc), ugi(ngi,ndim), u_ele(ndim,nloc,totele))
    allocate(xsgi(sngi,ndim), usgi(sngi,ndim), usgi2(sngi,ndim), income(sngi), snorm(sngi,ndim), norm(ndim))
    allocate(tnew_loc(nloc), tnew_loc2(nloc), told(nloc,totele), tnew(nloc,totele), tnew_nonlin(nloc,totele))
    allocate(tnew_nonlin_loc(nloc), tnew_nonlin2(nloc*totele))
    allocate(t_bc(nloc,totele), u_bc(ndim,nloc,totele))
    allocate(tnew_gi(ngi), tnew_xgi(ngi,ndim), tnew_sgi(sngi), tnew_sgi2(sngi), told_gi(ngi))
    allocate(face_sn(sngi,nloc,nface), face_sn2(sngi,nloc,n_s_list_no), face_snlx(sngi,sndim,nloc,nface) )
    allocate(x_all(ndim,nloc,totele), x_loc(ndim,nloc))
    allocate( mat_loc_inv(nloc,nloc), rhs_loc(nloc), ml_ele(nloc))
    allocate(rhs_jac(nloc), mass_t_new(nloc), inv_jac(ngi, ndim, ndim ) )
    allocate(a_star(ngi,nloc), p_star(ngi), rgi(ngi), diff_coe(ngi), told_loc(nloc), stab(nloc, nloc))
    allocate(mat_tnew(nloc), mat_diag_approx(nloc), mass_told(nloc))
    allocate(L1(ngi), L2(ngi), L3(ngi), L4(ngi))
    allocate(NLX_ALL(ngi,ndim,nloc))
    allocate(str_coor(ndim+1,nloc,totele))
    allocate(neig_info(totele), glob_flux(totele*nloc,totele*nloc))
    allocate(rhs_array(nloc*totele), inv_lhs(totele*nloc,totele*nloc))
    allocate( glob_mass_matrix(totele*nloc,totele*nloc), glob_lhs_matrix(totele*nloc,totele*nloc))
    allocate(glob_flux_matrix(totele*nloc,totele*nloc))



    t_bc(:,:)=0.0 ! this contains the boundary conditions just outside the domain
    u_bc(:,:,:)=0.0 ! this contains the boundary conditions on velocity just outside the domain
    tnew(:,:) = 0.0

    u_ele(:,:,:) = 0
    u_ele(1,:,1:totele) = u_x ! suggest this should be unity
    u_ele(2,:,1:totele) = u_y ! suggest this should be unity
    u_bc(:,:,:)=u_ele(:,:,:) ! this contains the boundary conditions on velocity just outside the domain
    ! dt = CFL/((u_ele(1,1,1)/dx)+(u_ele(2,1,1)/dy))
    dt = CFL*dx
    ntime= time/dt
print*, 'totele = ',totele
print*, 'ntime = ', ntime

    with_time_slab =.false.
    D3=.false.
    ! IQADRA = ipoly+1
    npoly=1; ele_type= 101

    call TRIQUAold(L1, L2, L3, L4, WEIGHT, D3,NGI)


    do ele=1,totele
      ! corner node coordinates
      x_all(1,1,ele) = meshList(ele)%X(1,1)   ; x_all(2,1,ele) = meshList(ele)%X(2,1)
      x_all(1,2,ele) = meshList(ele)%X(1,2)   ; x_all(2,2,ele) = meshList(ele)%X(2,2)
      x_all(1,3,ele) = meshList(ele)%X(1,3)   ; x_all(2,3,ele) = meshList(ele)%X(2,3)

      ! initial condition for 900 elements
      ! if ( meshList(ele)%region_id == 9 ) then
      !   tnew(:,ele) = 1.
      ! end if

      ! tnew(:,2) = 1.

      ! 100 elements
      ! tnew(:,2:3) = 1
      ! tnew(:,5) = 1
      ! tnew(:,7) = 1
      ! tnew(:,9) = 1
      ! tnew(:,11) = 1

      !10 elements
      ! tnew(:,1:2) = 1
      ! tnew(:,5) = 1
      ! tnew(:,9) = 1

      !4 eleents
      tnew(:,1)=1.
      tnew(:,3)=1.


      ! findin neighbiuring ele and face numbers
      do iface=1,nface
        face_ele(iface,ele) = meshList(ele)%Neig(iface)
        call getNeigDataMesh(meshlist, ele, iface, Npos, Nside, lnodes2)
        neig_info(ele)%Npos(iface) = Npos
        neig_info(ele)%Nside(iface) = Nside
        neig_info(ele)%Nnodes = lnodes2
      end do ! iface
    end do ! ele


    call make_sparse_matrix(sparse_lhs, totele, nloc)
    call make_sparse_matrix(sparse_mass, totele, nloc)
    call make_sparse_matrix_flux(sparse_flux, meshList, totele, nloc)

    call get_shape_funs_spec_ustr(n, nlx, nlx_lxx, nlxx, sn_orig, weight, nlx_nod, &
                nloc, snloc, sngi, ngi, ndim, nface,n_s_list_no, face_sn, face_sn2, face_snlx, &
                sweight, npoly, ele_type, totele, face_list_no )

    snlx_orig (:,1,1) = face_snlx(:,sndim,1,1) ! just copying from already calculated face_snlx to avoid doing it again
    snlx_orig (:,1,2) = face_snlx(:,sndim,3,1) ! it' copied from unstr_tri_surface_pointers_sn

    call CPU_TIME(t1)
    ! call DATE_AND_TIME(TIME = start)
    do itime=1,ntime

      ! generating VTK files
      if ( vtk_io <= itime  ) then
        call get_vtk(x_all, tnew, totele, nloc, itime,ndim, cell_type)
        vtk_io = vtk_io + vtk_interval
      end if

      told = tnew ! prepare for next time step
      tnew_nonlin = tnew
      do its=1,nits
        sparse_mass%val=0.
        sparse_lhs%val=0.
        ! glob_flux=0.
        sparse_flux%val=0.0
        tnew = tnew_nonlin ! for non-linearity
        do ele=1,totele
          ! volume integration
          x_loc=x_all(:,:,ele) ! x_all contains the coordinates of the corner nodes

          call tri_det_nlx( x_loc, n, nlx, nx, detwei, weight, ndim, nloc, ngi, inv_jac )

          !volume=sum(detwei)
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
          ! calculates t at gi sum = phi_1(gi)*t1 +phi_2(gi)*t2 +phi_3(gi)*t3 +phi_4(gi)*t4if(ele_type < is_triangle_or_tet) then
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
            diff_coe(gi) = 0.25*rgi(gi)**2 *p_star(gi) /max(toler, sum( tnew_xgi(gi,:)**2 ) )
          end do

          stab=0.0


          do iloc=1,nloc
            g_iloc = glob_no(ele, nloc, iloc)
            do jloc=1,nloc
              mass_ele=0.0
              stiff_ele = 0.0
              diff_ele =0.0
              lhs=0.0
              g_jloc = glob_no(ele, nloc, jloc)
              do gi=1,ngi
                stab(iloc,jloc) = stab(iloc,jloc) + sum(diff_coe(gi)* nx(gi,:,jloc)* nx(gi,:,iloc))* detwei(gi)
                mass_ele = mass_ele + n(gi,iloc)*n(gi,jloc)*detwei(gi)/dt
                do idim=1,ndim
                  stiff_ele = stiff_ele + nx(gi,idim,iloc)*ugi(gi,idim)*detwei(gi)*n(gi,jloc)
                  diff_ele  = diff_ele  + nx(gi,idim,iloc)*diff_coef*detwei(gi)*nx(gi,idim,jloc)
                end do
              end do ! quadrature
              call add_to_CSR(sparse_mass,g_iloc, g_jloc, mass_ele)
              lhs = mass_ele - stiff_ele - diff_ele
              call add_to_CSR(sparse_lhs,g_iloc, g_jloc, lhs)
            end do ! jloc
          end do ! iloc

          !Include the surface integral here:
          do iface = 1,nface
            ele22 = neig_info(ele)%Npos(iface) ! if ele22 is zero or negative assume on boundary
            i_got_boundary = (sign(1, -ele22) + 1 )/2 ! 1 means it is on the boundary
            r_got_boundary = real(i_got_boundary)
            ele2 = ele*i_got_boundary + ele22*(1-i_got_boundary)
            ! r_got_boundary=1.0 if we want to use the boundary conditions and have incomming velocity.
            ! r_got_boundary=0.0 not on boundary

            tnew_loc2(:)=tnew(:,ele2) * (1.0-r_got_boundary)   + t_bc(:,ele)  * r_got_boundary
            u_loc2(:,:)=u_ele(:,:,ele2)* (1.0-r_got_boundary)  + u_bc(:,:,ele)* r_got_boundary

            !Surface integral along an element
            sn   = face_sn(:,:,iface)             ! sn(sngi,nloc)
            snlx = face_snlx(:,:,:,iface)         ! slnlx(sngi, sndim, nloc)
            call get_unstr_sn2(face_sn2,nface,iface,neig_info(ele)%Nside(iface),sngi,nloc,snloc,SN_orig,snlx_orig)

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
            if ( target_ele==0 ) then ! it happens when ele is on the boundary
              target_ele=ele
            end if

            do iloc=1,nloc
              g_iloc=glob_no(ele,nloc,iloc)
              do jloc=1,nloc
                g_jloc=glob_no(target_ele, nloc, jloc)
                flux_ele = 0.0
                bc_diff = 0.0
                do idim=1,ndim
                  ! glob_flux(g_iloc,g_jloc)=glob_flux(g_iloc,g_jloc) +sum( snorm(:,idim)*sdetwei(:)*sn(:,iloc) &
                  !                                               *((1.-income(:))*sn(:,jloc)*usgi(:,idim) &
                  !                                               +     income(:)*sn2(:,jloc)*usgi2(:,idim)))

                  flux_ele = flux_ele +sum( snorm(:,idim)*sdetwei(:)*sn(:,iloc) &
                                                                *((1.-income(:))*sn(:,jloc)*usgi(:,idim) &
                                                                +     income(:)*sn2(:,jloc)*usgi2(:,idim)))
                  bc_diff = bc_diff +sum( snorm(:,idim)*sdetwei(:)*sn(:, iloc)*diff_coef &
                                                                *0.5*(snlx(:,sndim,jloc)))! deleted snlx2
                end do
                flux_ele = flux_ele + bc_diff
                call add_to_CSR_flux(sparse_flux, g_iloc, g_jloc,  flux_ele)
              end do ! jloc
            end do ! iloc
          end do ! iface


          ! ! if (with_stab) then
          !   mat_loc= mass_ele + dt*stab
          ! ! else
          ! !   mat_loc= mass_ele
          ! ! end if
          !
        end do ! do ele=1,totele


      call gen_global_matrix(totele, nloc, sparse_mass, glob_mass_matrix)
      call gen_global_matrix(totele, nloc, sparse_lhs, glob_lhs_matrix)
      call gen_global_matrix(totele, nloc, sparse_flux, glob_flux_matrix)
! if ( itime==1 .and. its==1 ) then
!   call print_my_matrix(glob_flux_matrix, totele*nloc)
! end if
      call told_to_array(totele, nloc, told, glob_told_array)
      rhs_array = matmul(glob_mass_matrix, glob_told_array)
      glob_lhs_matrix = glob_lhs_matrix + glob_flux_matrix

      call FINDInv(glob_lhs_matrix, inv_lhs, nloc*totele, errorflag)
      tnew_nonlin2 = matmul(inv_lhs,rhs_array )
      call array_to_told(totele, nloc, tnew_nonlin2, tnew_nonlin)
      tnew=tnew_nonlin
      end do ! do its=1,nits
print*, itime
    end do ! do itime=1,ntime

    !runing time
    call CPU_TIME(t2)
print*, 'cpu time_ReadMSH = ', time_ReadMSH
print*, ''
print*, 'cpu_time for time_loop = ', t2-t1



    deallocate (face_ele,face_list_no,n,nlx,nx,M,nlx_lxx,nlxx,nlx_nod,weight,detwei,sdetwei,&
                sn,sn2,snlx,sweight,s_cont,tnew_loc,tnew_loc2,told,tnew,tnew_nonlin,&
                tnew_nonlin_loc,tnew_nonlin2,t_bc,u_bc,tnew_gi,tnew_xgi,tnew_sgi,tnew_sgi2,told_gi,&
                face_sn,face_sn2,face_snlx,u_loc,u_ele,u_loc2,x_loc,ugi,x_all,xsgi,usgi,usgi2,&
                income,snorm,norm,mat_loc_inv,rhs_loc,ml_ele,rhs_jac,mass_t_new,SN_orig,&
                SNLX_orig,inv_jac,mat_diag_approx,a_star,p_star,rgi,diff_coe,told_loc,stab,mat_tnew,&
                mass_told,rhs_array,glob_lhs_matrix,glob_mass_matrix,glob_flux_matrix,glob_told_array,&
                inv_lhs,glob_flux,L1,L2,L3,L4,NLX_ALL)
  end subroutine diffusion

end module amin
