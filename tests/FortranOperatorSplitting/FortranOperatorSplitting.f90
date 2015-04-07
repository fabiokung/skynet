program operator_splitting

  use iso_c_binding
  use skynet
  implicit none

  interface

    integer(c_int) function check_results(num_zones, num_isotopes, num_times,   &
        p_total_heating, p_total_entropy_changes, p_final_mass_fractions)       &
        bind(C, name = 'CheckResults')

      use iso_c_binding
      implicit none

      integer(c_int), value :: num_zones, num_isotopes, num_times
      type(c_ptr), value :: p_total_heating, p_total_entropy_changes,           &
          p_final_mass_fractions

    end function

  end interface

  logical*1, parameter :: do_screening = .true.

  integer, parameter :: num_zones = 2
  integer, parameter :: num_isotopes = 13
  integer, parameter :: num_times = 11

  integer, target :: a(num_isotopes)
  integer, target :: z(num_isotopes)
  real*8, target :: composition(num_zones, num_isotopes)

  real*8 :: time(num_times)
  real*8 :: T9s(num_zones, num_times)
  real*8 :: rhos(num_zones, num_times)
  real*8, target :: T9(num_zones)
  real*8, target :: rho(num_zones)

  real*8, target :: heating(num_zones)
  real*8, target :: entropy_changes(num_zones)
  real*8, target :: total_heating(num_zones)
  real*8, target :: total_entropy_changes(num_zones)

  integer :: i, j, exit_status

  ! read a, z and initial composition
  open(unit = 1, file = "initial_compositions")

  ! skip first line
  read(1, *)

  ! read num_isotopes records
  do i = 1, num_isotopes
    read(1, '(5X,i4,i4,f5.1,f5.1)') a(i), z(i), composition(1, i),              &
        composition(2, i)
  end do

  close(1)


  ! read times, temperatures, and densities
  open(unit = 2, file = "temperature_density")

  ! skip first line
  read(2, *)

  ! read num_times records
  do i = 1, num_times
    read(2, '(e10.1,f7.2,f7.2,e10.1,e10.1)') time(i), T9s(1, i), T9s(2, i),     &
        rhos(1, i), rhos(2, i)
  end do

  close(2)


  ! init SkyNet
  call skynet_set_option_convergence_criterion(1);
  call skynet_set_option_mass_deviation_threshold(1.0d-6)
  call skynet_set_option_max_dt(1.0d-8)
  call skynet_set_option_do_screening(do_screening)

!  do i = 1, num_isotopes
!    write(6, '(f6.2,"  ",f6.2)') composition(1, i), composition(2, i)
!  end do

  T9(:) = T9s(:, 1)
  rho(:) = rhos(:, 1)
  call skynet_init(num_zones, num_isotopes, time(1), c_loc(T9), c_loc(rho),     &
      c_loc(a), c_loc(z), c_loc(composition))

  ! take steps
  do i = 1, num_zones
    total_heating(i) = 0.0d0
  end do

  do i = 1, num_times - 1
    T9(:) = T9s(:, i)
    rho(:) = rhos(:, i)
    call skynet_take_step(c_loc(T9), c_loc(rho), time(i + 1),  c_loc(heating),  &
        c_loc(entropy_changes), c_loc(composition))

!    write(6, *) "time = ", time(i + 1), ", heating = ", heating(1), ", ",       &
!        heating(2)

    do j = 1, num_zones
      total_heating(j) = total_heating(j) + heating(j)
      total_entropy_changes(j) = total_entropy_changes(j) + entropy_changes(j)
    end do
  end do

!  write(6, '("total heating: ",e10.3,e10.3)') total_heating(1), total_heating(2)

  ! clean up
  call skynet_finish()

  ! check the results with the same networks run in C++
  exit_status = check_results(num_zones, num_isotopes, num_times,               &
      c_loc(total_heating), c_loc(total_entropy_changes), c_loc(composition))

  call exit(exit_status)

end program operator_splitting
