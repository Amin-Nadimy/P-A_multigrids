program Transport_equation
  use precision
  use Structures
  use Generic
  use strings
  use evaluate
  use Msh2Tri
  use transport_rect
  use transport_tri
  use transport_tri_unstr
  use amin
  use transport_tri_semi

  implicit none

  integer :: mode = 9
  select case(mode)
      case(1)
        call trans_rec(0.7, 2, .false., 10, 250., 2, 2*100, 1, 4, 4, 2, 4, 2, 4, 2*0.01428571, 0.0, .true.)

      case(2)
        call str_explicit(0.7, 2, .false., 10, 0.07*2, 2, 707, 354, 3, 0.1, 0.1, 3, 2, 3, 2, 3, 0.1, 0.0, .false.,1,5)

      case(3)
        call str_implicit(0.7, 2, .false., 10, 0.07*30, 2, 2, 1, 3, 0.1, 0.1, 3, 2, 3, 2, 3, 0.1, 0.0, .true.,1,5)
!=====================================================================================================================================
      case(4)
        call unstr_explicit(0.7, 2, .false., 10, 0.7*0.1*1, 2,3, 0.001 , 0.001, 3, 2, 3, 2, 3, 0.9, 0., .true., 1, 5)

      case(5)
        call unstr_implicit(.7, 2, .false., 1, .7*.1*2, &
                                              2, 3, 0.1, 0.1, 3, 2, 3, 2, 3, -.1, .1, .true., 1, 5)

      case(6) ! amin.F90
        call diffusion(.7, 2, .false., 1, .7*.1*30, 2, 3, 0.1, 0.1, 3, 2, 3, 2, 3, 1., 0.0, .true., 1, 5, .5)
!=====================================================================================================================================
      case(7)
        call semi_explicit(0.7, 2, .false., 10,  0.7*0.1*40, 2, 3, 0.1, 0.1, 3, 2, 3, 2, 3, 0.1, 0.1, .true., 1, 5)

      case(8)
        call Semi_implicit_direct(.7, 2, .false., 10,  .7*0.5*10, 2, 3, 0.5, 0.5, 3, 2, 3, 2, 3, 0.1, 0.1, .true., 1, 5)

      case(9)
        ! (CFL, ndim, direct_solver, n_multigrid, time, nits, nface, dx, dy,&
        !           nloc, snloc, ngi, sngi, n_s_list_no, u_x, u_y, with_stab, vtk_interval, cell_type, solver, multi_levels,n_smooth)
        call Semi_implicit_iterative(1., 2, .false., 10, .025*0.005*150, 2, 3, 0.000125, 0.000125, 3, 2, 3, 2, 3, 0.,&
         0., .false., 1, 5, 3, 2,4)

      case(10)
        call Semi_implicit_iterative_P(.7, 2, .false., 10, .7*0.1*40, 2, 3, 0.1, 0.1, 3, 2, 3, 2, 3, 0.1, 0.1, .true., 1, 5)
  end select
end program Transport_equation
