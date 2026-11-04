program test_collision
  use neutron_collision, only: collision_data, dp
  implicit none
  type(collision_data) :: data
  real(dp) :: energies(6), energy, energy_out, mu, a, alpha, omega(3), out(3), u(2)
  real(dp) :: st, se, norm, recoil, expected, integral, moment, dx, x, weight, cumulative
  real(dp) :: last_cdf, sum_mu, sum_ratio, dsdo
  integer :: i,j,k, n
  character(64) :: mode

  call data%load('data/u238')
  call get_command_argument(1,mode)
  if (trim(mode) == 'bad_direction') then
    call data%scatter(1.e6_dp,[0.0_dp,0.0_dp,0.0_dp],energy_out,out)
    error stop 'Invalid direction was accepted'
  else if (trim(mode) == 'bad_uniform') then
    call data%scatter(1.e6_dp,[0.0_dp,0.0_dp,1.0_dp],energy_out,out,uniforms=[1.0_dp,0.0_dp])
    error stop 'Invalid uniform was accepted'
  end if
  a = data%mass_ratio()
  call check(abs(a-236.0058_dp) < 1.e-10_dp,'AWR must come from evaluation, not integer 238')
  alpha = ((a-1)/(a+1))**2
  energies = [1.5e5_dp,1.e6_dp,1.1e6_dp,1.e7_dp,1.0000005e7_dp,1.9e7_dp]
  do k = 1,size(energies)
    energy = energies(k)
    call data%cross_sections(energy,st,se)
    call check(st >= se .and. se > 0,'cross section ordering')
    call check(abs(data%angular_cdf(energy,-1.0_dp)) < 1.e-13_dp,'CDF(-1)')
    call check(abs(data%angular_cdf(energy,1.0_dp)-1) < 1.e-13_dp,'CDF(1)')
    ! Composite Simpson quadrature checks normalization, first moment and dsigma/dOmega.
    n = 8192
    dx = 2.0_dp/n
    integral = 0
    moment = 0
    dsdo = 0
    last_cdf = -1
    do j = 0,n
      x = -1 + dx*j
      weight = 2
      if (mod(j,2) == 1) weight = 4
      if (j == 0 .or. j == n) weight = 1
      integral = integral + weight*data%angular_pdf(energy,x)
      moment = moment + weight*x*data%angular_pdf(energy,x)
      dsdo = dsdo + weight*data%differential_cross_section(energy,x)
      cumulative = data%angular_cdf(energy,x)
      call check(cumulative >= last_cdf-1.e-13_dp,'CDF monotonicity')
      last_cdf = cumulative
    end do
    call check(abs(integral*dx/3-1) < 1.e-8_dp,'PDF normalization')
    call check(abs(moment*dx/3-data%mean_cosine_cm(energy)) < 1.e-8_dp,'Legendre first moment')
    call check(abs(dsdo*dx/3*2*acos(-1.0_dp)-se) < 1.e-7_dp,'differential integrates to elastic')

    sum_mu = 0
    sum_ratio = 0
    n = 4000
    do i = 1,n
      ! Stratified inverse-CDF inputs avoid stochastic failures of moment checks.
      u = [(i-0.5_dp)/n,modulo(i*0.6180339887498949_dp,1.0_dp)]
      select case(mod(i,4))
      case(0)
        omega = [0.0_dp,0.0_dp,1.0_dp]
      case(1)
        omega = [0.0_dp,0.0_dp,-1.0_dp]
      case(2)
        omega = [1.0_dp,0.0_dp,0.0_dp]
      case default
        omega = [1.0_dp,2.0_dp,3.0_dp]/sqrt(14.0_dp)
      end select
      call data%scatter(energy,omega,energy_out,out,mu,u)
      norm = sum(out**2)
      call check(abs(norm-1) < 1.e-13_dp,'unit outgoing direction')
      call check(energy_out >= alpha*energy .and. energy_out <= energy,'kinematic energy bounds')
      expected = energy*(a*a+1+2*a*mu)/(a+1)**2
      call check(abs(energy_out-expected)/energy < 1.e-13_dp,'angle-energy correlation')
      expected = (1+a*mu)/sqrt(a*a+1+2*a*mu)
      call check(abs(dot_product(omega,out)-expected) < 1.e-13_dp,'laboratory angle')
      recoil = sum((sqrt(energy)*omega-sqrt(energy_out)*out)**2)/a
      call check(abs(energy_out+recoil-energy)/energy < 2.e-13_dp,'energy and momentum conservation')
      call check(abs(data%angular_cdf(energy,mu)-u(1)) < 1.e-11_dp,'inverse CDF residual')
      sum_mu = sum_mu + mu
      sum_ratio = sum_ratio + energy_out/energy
    end do
    call check(abs(sum_mu/n-data%mean_cosine_cm(energy)) < 2.e-5_dp,'sample angular mean')
    expected = (a*a+1+2*a*data%mean_cosine_cm(energy))/(a+1)**2
    call check(abs(sum_ratio/n-expected) < 2.e-7_dp,'sample energy mean')
  end do
  ! Reloading must replace, rather than leak or mix, data.
  call data%load('data/u238')
  print '(a)', 'PASS: ENDF parsing, six energies, PDF/CDF, moments, 24000 collisions, conservation, reload'
contains
  subroutine check(condition,message)
    logical,intent(in) :: condition
    character(*),intent(in) :: message
    if (.not. condition) then
      print '(a)', 'FAIL: '//message
      error stop 1
    end if
  end subroutine
end program
