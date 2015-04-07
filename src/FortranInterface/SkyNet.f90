! This module contains bindings to C functions in SkyNet so that the SkyNet
! operator split functions can be called from FORTRAN.

module skynet

  interface

    ! The skynet_set_option_* functions must be called BEFORE skynet_init, if
    ! they are called after skynet_init, they don't have any effect (until
    ! skynet_init is called again).
    !
    ! See the file src/Network/NetworkOptions.hpp for more information about the
    ! available options. Currently, not all options are setable via the Fortran
    ! interface.

    ! set the convergence criterion (default Mass)
    subroutine skynet_set_option_convergence_criterion(convergence_criterion)   &
        bind(C, name = 'SkyNetSetConvergenceCriterion')
      use iso_c_binding
      implicit none

      ! [in] set the convergence criterion for SkyNet:
      !   1: DeltaYByY
      !   2: Mass
      !   3: BothDeltaYByYAndMass
      integer(c_int), value :: convergence_criterion
    end subroutine skynet_set_option_convergence_criterion


    ! set the DeltaYByY convergence threshold (default 1.0E-10)
    subroutine skynet_set_option_delta_y_by_y_threshold(delta_y_by_y_threshold) &
        bind(C, name = 'SkyNetSetOptionDeltaYByYThreshold')
      use iso_c_binding
      implicit none

      ! [in]  with the DeltaYByY convergence criterion, the Newton-Raphson
      ! iteration is considered to have converged when the error is smaller than
      ! this threshold
      real(c_double), value :: delta_y_by_y_threshold
    end subroutine skynet_set_option_delta_y_by_y_threshold


    ! set the mass convergence criterion (default 1.0E-6)
    subroutine skynet_set_option_mass_deviation_threshold(                      &
        mass_deviation_threshold)                                               &
        bind(C, name = 'SkyNetSetOptionMassDeviationThreshold')
      use iso_c_binding
      implicit none

      ! [in] if the total mass differs from 1 by less than this threshold, then
      ! the Newton-Raphson step is considered converged
      real(c_double), value :: mass_deviation_threshold
    end subroutine skynet_set_option_mass_deviation_threshold


    ! set the minimum time step size (default 1.0E-16)
    subroutine skynet_set_option_min_dt(min_dt)                                 &
        bind(C, name = 'SkyNetSetOptionMinDt')
      use iso_c_binding
      implicit none

      ! [in] smallest allowed time step size
      real(c_double), value :: min_dt
    end subroutine skynet_set_option_min_dt


    ! set the maximum time step size (default 1.0E16)
    subroutine skynet_set_option_max_dt(max_dt)                                 &
        bind(C, name = 'SkyNetSetOptionMaxDt')
      use iso_c_binding
      implicit none

      ! [in] smallest allowed time step size
      real(c_double), value :: max_dt
    end subroutine skynet_set_option_max_dt


    ! set whether to do screening or not
    subroutine skynet_set_option_do_screening(do_screening)                     &
        bind(C, name = 'SkyNetSetOptionDoScreening')
      use iso_c_binding
      implicit none

      ! [in] whether to do screening
      logical(c_bool), value :: do_screening
    end subroutine skynet_set_option_do_screening

    ! initialize the nuclide library (this is done automatically by skynet_init,
    ! but with this function, one can only inititalize the nuclide library
    ! without creating any networks, which is needed to use NSE and the EOS)
    subroutine skynet_init_nuclide_library(num_isotopes, p_as, p_zs)            &
        bind(C, name = 'SkyNetInitNuclideLibrary')
      use iso_c_binding
      implicit none

      ! [in] number of isotopes to use in the network
      integer(c_int), value :: num_isotopes

      ! [in] pointer (ints) to the A values of the isotopes
      type(c_ptr), value :: p_as

      ! [in] pointer (ints) to the Z values of the isotopes
      type(c_ptr), value :: p_zs
    end subroutine skynet_init_nuclide_library


    ! initialize SkyNet and set initial mass fractions
    subroutine skynet_init(num_zones, num_isotopes, start_time, p_initial_T9s,  &
        p_initial_rhos, p_as, p_zs, p_initial_mass_fractions)                   &
        bind(C, name = 'SkyNetInit')
      use iso_c_binding
      implicit none

      ! [in] number of zones, i.e. copies of SkyNet to run
      integer(c_int), value :: num_zones

      ! [in] number of isotopes to use in the network
      integer(c_int), value :: num_isotopes

      ! [in] start time of the SkyNet evolution (in seconds)
      real(c_double), value :: start_time

      ! [in] pointer (doubles) to the initial temperature values (in Giga
      ! Kelvin) of the zones
      type(c_ptr), value :: p_initial_T9s

      ! [in] pointer (doubles) to the initial density values (in g / cm^3) of
      ! the zones
      type(c_ptr), value :: p_initial_rhos

      ! [in] pointer (ints) to the A values of the isotopes
      type(c_ptr), value :: p_as

      ! [in] pointer (ints) to the Z values of the isotopes
      type(c_ptr), value :: p_zs

      ! [in] pointer (doubles) to the initial mass fractions of all the zones,
      ! the FORTRAN array is initial_mass_fractions(num_zones, num_isotopes)
      type(c_ptr), value :: p_initial_mass_fractions
    end subroutine skynet_init


    ! evolve SkyNet to the given end time at constant temperature and density,
    ! returns the amount of heating that has occurred in this step in the zones
    subroutine skynet_take_step(p_T9s, p_rhos, end_time, p_specific_heating,    &
        p_entropy_changes, p_mass_fractions) bind(C, name = 'SkyNetTakeStep')
      use iso_c_binding
      implicit none

      ! [in] pointer (doubles) to the constant temperature values (in Giga
      ! Kelvin) of the zones
      type(c_ptr), value :: p_T9s

      ! [in] pointer (doubles) to the constant density values (in g / cm^3) of
      ! the zones
      type(c_ptr), value :: p_rhos

      ! [in] the end time for the SkyNet evolution step (in seconds)
      real(c_double), value :: end_time

      ! [out] pointer (doubles) to the amount of specific heating (in erg / g)
      ! that has occurred during this SkyNet step
      type(c_ptr), value :: p_specific_heating

      ! [out] pointer (doubles) to the change in entropy (in k_b / baryon)
      ! that has occurred during this SkyNet step
      type(c_ptr), value :: p_entropy_changes

      ! [out] pointer (doubles) to the mass fractions of all the zones at the
      ! end of the step, the FORTRAN array is
      ! mass_fractions(num_zones, num_isotopes)
      type(c_ptr), value :: p_mass_fractions
    end subroutine skynet_take_step


    ! finishes the SkyNet evolution in all zones and destroys the networks
    subroutine skynet_finish() bind(C, name = 'SkyNetFinish')
      use iso_c_binding
      implicit none
    end subroutine skynet_finish


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!                                                                       !!
    !!  NSE INTERFACE                                                        !!
    !!                                                                       !!
    !!  MUST CALL skynet_init_nuclide_library OR skynet_init FIRST           !!
    !!                                                                       !!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! calcualte NSE from temperature and density
    subroutine skynet_nse_from_t_rho(T9, rho, Ye, p_mass_fractions)             &
        bind(C, name = 'SkyNetNSECalcFromTemperatureAndDensity')
      use iso_c_binding
      implicit none

      ! [in] temperature (in Giga Kelvin)
      real(c_double), value :: T9

      ! [in] density (in gm / cm^3)
      real(c_double), value :: rho

      ! [in] electron faction
      real(c_double), value :: Ye

      ! [out] pointer (doubles) to the NSE mass fractions, the FORTRAN array is
      ! mass_fractions(num_isotopes)
      type(c_ptr), value :: p_mass_fractions
    end subroutine skynet_nse_from_t_rho


    ! calcualte NSE from temperature and entropy
    subroutine skynet_nse_from_t_s(T9, entropy, Ye, p_mass_fractions, p_rho)    &
        bind(C, name = 'SkyNetNSECalcFromTemperatureAndEntropy')
      use iso_c_binding
      implicit none

      ! [in] temperature (in Giga Kelvin)
      real(c_double), value :: T9

      ! [in] entropy (in k_B / baryon)
      real(c_double), value :: entropy

      ! [in] electron faction
      real(c_double), value :: Ye

      ! [out] pointer (doubles) to the NSE mass fractions, the FORTRAN array is
      ! mass_fractions(num_isotopes)
      type(c_ptr), value :: p_mass_fractions

      ! [out] the NSE density (in gm / cm^3)
      real(c_double), value :: p_rho
    end subroutine skynet_nse_from_t_s


    ! calcualte NSE from entropy and density
    subroutine skynet_nse_from_s_rho(entropy, rho, Ye, p_mass_fractions, p_T9)  &
        bind(C, name = 'SkyNetNSECalcFromEntropyAndDensity')
      use iso_c_binding
      implicit none

      ! [in] entropy (in k_B / baryon)
      real(c_double), value :: entropy

      ! [in] density (in gm / cm^3)
      real(c_double), value :: rho

      ! [in] electron faction
      real(c_double), value :: Ye

      ! [out] pointer (doubles) to the NSE mass fractions, the FORTRAN array is
      ! mass_fractions(num_isotopes)
      type(c_ptr), value :: p_mass_fractions

      ! [out] the NSE temperature (in Giga Kelvin)
      type(c_ptr), value :: p_T9
    end subroutine skynet_nse_from_s_rho


    ! calcualte NSE from entropy and density with a guess for the temperature
    subroutine skynet_nse_from_s_rho_with_t_guess(entropy, rho, Ye, T9_guess,   &
        p_mass_fractions, p_T9)                                                 &
        bind(C, name = 'SkyNetNSECalcFromEntropyAndDensityWithTemperatureGuess')
      use iso_c_binding
      implicit none

      ! [in] entropy (in k_B / baryon)
      real(c_double), value :: entropy

      ! [in] density (in gm / cm^3)
      real(c_double), value :: rho

      ! [in] electron faction
      real(c_double), value :: Ye

      ! [in] guess for the temperature (in Giga Kelvin)
      real(c_double), value :: T9_guess

      ! [out] pointer (doubles) to the NSE mass fractions, the FORTRAN array is
      ! mass_fractions(num_isotopes)
      type(c_ptr), value :: p_mass_fractions

      ! [out] the NSE temperature (in Giga Kelvin)
      type(c_ptr), value :: p_T9
    end subroutine skynet_nse_from_s_rho_with_t_guess


    ! calcualte NSE from internal energy and density
    subroutine skynet_nse_from_u_rho(u, rho, Ye, p_mass_fractions, p_T9)        &
        bind(C, name = 'SkyNetNSECalcFromInternalEnergyAndDensity')
      use iso_c_binding
      implicit none

      ! [in] specific internal energy (in erg / g)
      real(c_double), value :: u

      ! [in] density (in gm / cm^3)
      real(c_double), value :: rho

      ! [in] electron faction
      real(c_double), value :: Ye

      ! [out] pointer (doubles) to the NSE mass fractions, the FORTRAN array is
      ! mass_fractions(num_isotopes)
      type(c_ptr), value :: p_mass_fractions

      ! [out] the NSE temperature (in Giga Kelvin)
      type(c_ptr), value :: p_T9
    end subroutine skynet_nse_from_u_rho


    ! calcualte NSE from internal energy and density with a guess for the
    ! temperature
    subroutine skynet_nse_from_u_rho_with_t_guess(u, rho, Ye, T9_guess,         &
        p_mass_fractions, p_T9)                                                 &
        bind(C, name =                                                          &
            'SkyNetNSECalcFromInternalEnergyAndDensityWithTemperatureGuess')
      use iso_c_binding
      implicit none

      ! [in] specific internal energy (in erg / g)
      real(c_double), value :: u

      ! [in] density (in gm / cm^3)
      real(c_double), value :: rho

      ! [in] electron faction
      real(c_double), value :: Ye

      ! [in] guess for the temperature (in Giga Kelvin)
      real(c_double), value :: T9_guess

      ! [out] pointer (doubles) to the NSE mass fractions, the FORTRAN array is
      ! mass_fractions(num_isotopes)
      type(c_ptr), value :: p_mass_fractions

      ! [out] the NSE temperature (in Giga Kelvin)
      type(c_ptr), value :: p_T9
    end subroutine skynet_nse_from_u_rho_with_t_guess


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!                                                                       !!
    !!  EOS INTERFACE                                                        !!
    !!                                                                       !!
    !!  MUST CALL skynet_init_nuclide_library OR skynet_init FIRST           !!
    !!                                                                       !!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! call EOS given temperature and density
    subroutine skynet_eos_from_t_rho(T9, rho, p_mass_fractions, p_s, p_u,       &
        p_dSdT9) bind(C, name = 'SkyNetEOSCalcFromTAndRho')
      use iso_c_binding
      implicit none

      ! [in] temperature (in Giga Kelvin)
      real(c_double), value :: T9

      ! [in] density (in gm / cm^3)
      real(c_double), value :: rho

      ! [in] pointer (doubles) to the mass fractions, the FORTRAN array is
      ! mass_fractions(num_isotopes)
      type(c_ptr), value :: p_mass_fractions

      ! [out] the EOS entropy (in k_b / baryon)
      type(c_ptr), value :: p_s

      ! [out] the EOS specific internal energy (in erg / gm)
      type(c_ptr), value :: p_u

      ! [out] the EOS entropy derivative w.r.t. temperature
      ! (in k_b / baryon / GK)
      type(c_ptr), value :: p_dSdT9
    end subroutine skynet_eos_from_t_rho


    ! call EOS given entropy and density
    subroutine skynet_eos_from_s_rho(entropy, rho, p_mass_fractions, p_T9, p_u, &
        p_dSdT9) bind(C, name = 'SkyNetEOSCalcFromSAndRho')
      use iso_c_binding
      implicit none

      ! [in] entropy (in k_b / baryon)
      real(c_double), value :: entropy

      ! [in] density (in gm / cm^3)
      real(c_double), value :: rho

      ! [in] pointer (doubles) to the mass fractions, the FORTRAN array is
      ! mass_fractions(num_isotopes)
      type(c_ptr), value :: p_mass_fractions

      ! [out] the EOS temperature (in Giga Kelvin)
      type(c_ptr), value :: p_T9

      ! [out] the EOS specific internal energy (in erg / gm)
      type(c_ptr), value :: p_u

      ! [out] the EOS entropy derivative w.r.t. temperature
      ! (in k_b / baryon / GK)
      type(c_ptr), value :: p_dSdT9
    end subroutine skynet_eos_from_s_rho


    ! call EOS given entropy and density with a guess for the temperature
    subroutine skynet_eos_from_s_rho_with_t_guess(entropy, rho, T9_guess,       &
        p_mass_fractions, p_T9, p_u, p_dSdT9)                                        &
        bind(C, name = 'SkyNetEOSCalcFromSAndRhoWithTGuess')
      use iso_c_binding
      implicit none

      ! [in] entropy (in k_b / baryon)
      real(c_double), value :: entropy

      ! [in] density (in gm / cm^3)
      real(c_double), value :: rho

      ! [in] guess for the temperature (in Giga Kelvin)
      real(c_double), value :: T9_guess

      ! [in] pointer (doubles) to the mass fractions, the FORTRAN array is
      ! mass_fractions(num_isotopes)
      type(c_ptr), value :: p_mass_fractions

      ! [out] the EOS temperature (in Giga Kelvin)
      type(c_ptr), value :: p_T9

      ! [out] the EOS specific internal energy (in erg / gm)
      type(c_ptr), value :: p_u

      ! [out] the EOS entropy derivative w.r.t. temperature
      ! (in k_b / baryon / GK)
      type(c_ptr), value :: p_dSdT9
    end subroutine skynet_eos_from_s_rho_with_t_guess


    ! call EOS given temperature and entropy
    subroutine skynet_eos_from_t_s(T9, entropy, p_mass_fractions, p_rho, p_u,   &
        p_dSdT9) bind(C, name = 'SkyNetEOSCalcFromTAndS')
      use iso_c_binding
      implicit none

      ! [in] temperature (in Giga Kelvin)
      real(c_double), value :: T9

      ! [in] entropy (in k_b / baryon)
      real(c_double), value :: entropy

      ! [in] pointer (doubles) to the mass fractions, the FORTRAN array is
      ! mass_fractions(num_isotopes)
      type(c_ptr), value :: p_mass_fractions

      ! [out] the EOS density (in gm / cm^3)
      type(c_ptr), value :: p_rho

      ! [out] the EOS specific internal energy (in erg / gm)
      type(c_ptr), value :: p_u

      ! [out] the EOS entropy derivative w.r.t. temperature
      ! (in k_b / baryon / GK)
      type(c_ptr), value :: p_dSdT9
    end subroutine skynet_eos_from_t_s


    ! call EOS given temperature and entropy with a guess for the density
    subroutine skynet_eos_from_t_s_with_rho_guess(T9, entropy, rho_guess,       &
        p_mass_fractions, p_rho, p_u, p_dSdT9)                                  &
        bind(C, name = 'SkyNetEOSCalcFromTAndSWithRhoGuess')
      use iso_c_binding
      implicit none

      ! [in] temperature (in Giga Kelvin)
      real(c_double), value :: T9

      ! [in] entropy (in k_b / baryon)
      real(c_double), value :: entropy

      ! [in] guess for the density (in gm / cm^3)
      real(c_double), value :: rho_guess

      ! [in] pointer (doubles) to the mass fractions, the FORTRAN array is
      ! mass_fractions(num_isotopes)
      type(c_ptr), value :: p_mass_fractions

      ! [out] the EOS density (in gm / cm^3)
      type(c_ptr), value :: p_rho

      ! [out] the EOS specific internal energy (in erg / gm)
      type(c_ptr), value :: p_u

      ! [out] the EOS entropy derivative w.r.t. temperature
      ! (in k_b / baryon / GK)
      type(c_ptr), value :: p_dSdT9
    end subroutine skynet_eos_from_t_s_with_rho_guess

  end interface

end module skynet
