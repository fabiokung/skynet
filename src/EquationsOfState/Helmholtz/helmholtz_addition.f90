! wrapper to call helmeos on a vector of length n (but with constant density, abar and zbar)
subroutine evaluate_helm_eos(rho, abar, zbar, n, success_flag)

  ! read_helm_table must have been called already
  ! temp_row must also have been set

  implicit none
  include 'vector_eos.dek'
  real(8) :: rho, abar, zbar
  integer :: n, j
  logical :: success_flag

  real(8) :: ye

  if (n < 1) then
    stop 'n must be at least 1 in evaluate_helm_eos.'
  endif

  ye = max(1.0d-16,(1.0d0/abar) * zbar)

  if (ye * rho < 1.0d-12) then
    write(6,*) 'ye * rho < 1e-12'
    success_flag = .False.
    return
  endif

  jlo_eos = 1
  jhi_eos = n

  do j = jlo_eos, jhi_eos
    if (temp_row(j) < 1.0d3) then
      write(6,*) 'temperature < 1e3'
      success_flag = .False.
      return
    endif

    den_row(j)  = rho
    abar_row(j) = abar
    zbar_row(j) = zbar
  enddo

  call helmeos

  if (eosfail) then
    write(6,*) 'eosfail is true'
    success_flag = .False.
    return
  endif

  ! throw away ion contributions, they will be calculated in SkyNet
  ! since SkyNet does not include screening (yet), we also throw out the
  ! coulomb corrections, for consistency
  do j = jlo_eos, jhi_eos
!     ! with coulomb corrections
!     stot_row(j) = srad_row(j) + sele_row(j) + scou_row(j)

    ! without coulomb corrections
    stot_row(j) = srad_row(j) + sele_row(j)
    etot_row(j) = erad_row(j) + eele_row(j)
    dst_row(j) = dsradt_row(j) + dsept_row(j)
  enddo

  success_flag = .True.

end subroutine evaluate_helm_eos

! wrapper for read_helm_table that remembers whether the table has already been read
subroutine read_table(helm_table_path, str_length)
  logical :: table_read
  COMMON table_read

  integer, intent(in) :: str_length
  character(str_length), intent(in) :: helm_table_path

  call read_helm_table(helm_table_path)
  table_read = .true.

end subroutine read_table

! Wrapper for the eos routine
subroutine Helm_eos_wrap(T, rho, Abar, Zbar, electronPhotonEntropy,             &
    electronPhotonInternalEnergy, electronEta,                                  &
    electronPhotonEntropyDerivativeWithTemperature, successFlag)

  implicit none
  include 'vector_eos.dek'

  logical :: table_read
  COMMON table_read

  ! input parameters
  real*8 :: T      ! temperature in Kelvin
  real*8 :: rho    ! density in grams per cm^3
  real*8 :: Abar   ! average mass nubmer A
  real*8 :: Zbar   ! average atomic number Z

  ! output parameters

  ! specific entropy of electrons and photons in ergs per gram per Kelvin
  real*8 :: electronPhotonEntropy

  ! internal energy of electrons and photons in ergs per gram
  real*8 :: electronPhotonInternalEnergy

  ! Electron degeneracy parameter
  real*8 :: electronEta

  ! derivative of the specific entropy of electrons and photons
  ! in ergs per gram per Kelvin^2
  real*8 :: electronPhotonEntropyDerivativeWithTemperature

  logical :: successFlag ! true

  if (.not. table_read) then
    write(6,*) 'must read helm_table.dat before calling Helmholtz EOS'
    successFlag = .false.
    return
  endif

  temp_row(1) = T
  call evaluate_helm_eos(rho, Abar, Zbar, 1, successFlag)

  if (.not. successFlag) then
    electronPhotonEntropy = -1.0
    electronPhotonInternalEnergy = -1.0
    electronEta = -1.0

    electronPhotonEntropyDerivativeWithTemperature = 0.0
    return
  endif

  electronPhotonEntropy = stot_row(1)
  electronPhotonInternalEnergy = etot_row(1)
  electronEta = etaele_row(1)

  electronPhotonEntropyDerivativeWithTemperature = dst_row(1)

end subroutine Helm_eos_wrap
