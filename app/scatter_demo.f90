program scatter_demo
  use neutron_collision, only: collision_data, dp
  use iso_fortran_env, only: int64
  implicit none
  type(collision_data) :: data
  real(dp) :: energy, energy_out, direction(3), mu, sigma_total, sigma_elastic
  character(512) :: argument, directory
  integer :: n, i, ios, seed_size, seed_value
  integer, allocatable :: seed(:)

  energy = 1.e6_dp
  n = 1
  seed_value = 12345
  directory = 'data/u238'
  if (command_argument_count() > 4) error stop 'Usage: scatter_demo [E_eV] [N] [seed] [data_directory]'
  if (command_argument_count() >= 1) then
    call get_command_argument(1,argument)
    read(argument,*,iostat=ios) energy
    if (ios /= 0) error stop 'Invalid incident energy'
  end if
  if (command_argument_count() >= 2) then
    call get_command_argument(2,argument)
    read(argument,*,iostat=ios) n
    if (ios /= 0) error stop 'Invalid sample count'
  end if
  if (n < 1 .or. n > 10000000) error stop 'N must be in [1,10000000]'
  if (command_argument_count() >= 3) then
    call get_command_argument(3,argument)
    read(argument,*,iostat=ios) seed_value
    if (ios /= 0) error stop 'Invalid seed'
  end if
  if (seed_value < 0) error stop 'Seed must be nonnegative'
  if (command_argument_count() >= 4) call get_command_argument(4,directory)
  call random_seed(size=seed_size)
  allocate(seed(seed_size))
  do i = 1, seed_size
    seed(i) = int(modulo(int(seed_value,int64)+104729_int64*i, 2147483647_int64))
  end do
  call random_seed(put=seed)
  call data%load(trim(directory))
  call data%cross_sections(energy,sigma_total,sigma_elastic)
  print '(a)', '# U-238 ENDF/B-VII.1; elastic collision conditional on MT=2; stationary target'
  print '(a,es24.16)', '# E_in_eV = ', energy
  print '(a,es24.16)', '# AWR = ', data%mass_ratio()
  print '(a,es24.16)', '# sigma_total_barn = ', sigma_total
  print '(a,es24.16)', '# sigma_elastic_barn = ', sigma_elastic
  print '(a,es24.16)', '# mean_mu_cm = ', data%mean_cosine_cm(energy)
  print '(a)', 'E_out_eV,omega_x,omega_y,omega_z,mu_cm'
  do i = 1, n
    call data%scatter(energy,[0.0_dp,0.0_dp,1.0_dp],energy_out,direction,mu)
    print '(es24.16,4(",",es24.16))', energy_out,direction,mu
  end do
end program
