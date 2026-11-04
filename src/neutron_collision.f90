! ENDF/B-VII.1 U-238: conditional elastic collision, stationary target.
module neutron_collision
  use iso_fortran_env, only: real64
  use ieee_arithmetic, only: ieee_is_finite
  implicit none
  private
  integer, parameter, public :: dp = real64
  real(dp), parameter :: pi = acos(-1.0_dp)
  ! Above the unresolved resonance region; below nonpositive high-E PDFs.
  real(dp), parameter :: supported_min = 1.5e5_dp, supported_max = 1.9e7_dp
  integer, parameter :: max_order = 64

  type :: cross_section_table
    real(dp), allocatable :: energy(:), sigma(:)
  end type

  type, public :: collision_data
    private
    real(dp) :: awr = 0.0_dp
    type(cross_section_table) :: total, elastic
    real(dp), allocatable :: angular_energy(:), coefficient(:,:)
  contains
    procedure :: load => load_data
    procedure :: cross_sections
    procedure :: angular_pdf
    procedure :: angular_cdf
    procedure :: differential_cross_section
    procedure :: mean_cosine_cm
    procedure :: scatter
    procedure :: mass_ratio
  end type

contains

  subroutine require(condition, message)
    logical, intent(in) :: condition
    character(*), intent(in) :: message
    if (.not. condition) then
      print '(a)', 'neutron_collision: '//message
      error stop 1
    end if
  end subroutine

  ! ENDF numbers have an implicit E: " 1.234567+6".
  function endf_real(field) result(value)
    character(*), intent(in) :: field
    real(dp) :: value
    character(:), allocatable :: text
    integer :: i, ios
    text = trim(adjustl(field))
    if (len(text) == 0) then
      value = 0.0_dp
      return
    end if
    if (scan(text, 'eEdD') == 0) then
      do i = len(text), 2, -1
        if (text(i:i) == '+' .or. text(i:i) == '-') then
          text = text(:i-1)//'E'//text(i:)
          exit
        end if
      end do
    end if
    read(text, *, iostat=ios) value
    call require(ios == 0, 'invalid ENDF numeric field: '//field)
    call require(ieee_is_finite(value), 'nonfinite ENDF value')
  end function

  subroutine record(unit, mf, mt, values)
    integer, intent(in) :: unit, mf, mt
    real(dp), intent(out) :: values(6)
    character(256) :: line
    integer :: ios, actual_mf, actual_mt, mat, i
    read(unit, '(a)', iostat=ios) line
    call require(ios == 0, 'truncated or unreadable ENDF section')
    read(line(67:75), '(i4,i2,i3)', iostat=ios) mat, actual_mf, actual_mt
    call require(ios == 0, 'invalid ENDF section identifier')
    call require(mat == 9237 .and. actual_mf == mf .and. actual_mt == mt, &
                 'expected U-238 MAT=9237 and requested MF/MT')
    do i = 1, 6
      values(i) = endf_real(line(11*i-10:11*i))
    end do
  end subroutine

  subroutine values_array(unit, mf, mt, values)
    integer, intent(in) :: unit, mf, mt
    real(dp), intent(out) :: values(:)
    real(dp) :: v(6)
    integer :: i, n
    do i = 1, size(values), 6
      call record(unit, mf, mt, v)
      n = min(6, size(values)-i+1)
      values(i:i+n-1) = v(:n)
    end do
  end subroutine

  subroutine open_input(path, unit)
    character(*), intent(in) :: path
    integer, intent(out) :: unit
    integer :: ios
    open(newunit=unit, file=path, status='old', action='read', iostat=ios)
    call require(ios == 0, 'cannot open '//path)
  end subroutine

  subroutine read_xs(path, mt, table, awr)
    character(*), intent(in) :: path
    integer, intent(in) :: mt
    type(cross_section_table), intent(out) :: table
    real(dp), intent(out) :: awr
    integer :: unit, n
    real(dp) :: v(6), interpolation(2)
    real(dp), allocatable :: pairs(:)
    call open_input(path, unit)
    call record(unit, 3, mt, v)
    call require(nint(v(1)) == 92238 .and. v(2) > 1, 'wrong isotope or mass ratio')
    awr = v(2)
    call record(unit, 3, mt, v)
    call require(nint(v(5)) == 1, 'only one MF3 interpolation region is supported')
    n = nint(v(6))
    call require(n >= 2 .and. n <= 2000000, 'invalid MF3 point count')
    call values_array(unit, 3, mt, interpolation)
    call require(nint(interpolation(1)) == n .and. nint(interpolation(2)) == 2, &
                 'MF3 requires ENDF linear-linear interpolation (INT=2)')
    allocate(pairs(2*n), table%energy(n), table%sigma(n))
    call values_array(unit, 3, mt, pairs)
    table%energy = pairs(1::2)
    table%sigma = pairs(2::2)
    call require(all(table%energy(2:) > table%energy(:n-1)), 'MF3 energies must increase')
    call require(all(table%sigma >= 0), 'negative cross section')
    call require(table%energy(1) <= supported_min .and. table%energy(n) >= supported_max, &
                 'cross section does not cover supported energies')
    close(unit)
  end subroutine

  subroutine load_data(self, directory)
    class(collision_data), intent(inout) :: self
    character(*), intent(in) :: directory
    integer :: unit, n, i, order, j
    real(dp) :: v(6), interpolation(2), other_awr, a(max_order), energy, mu
    call read_xs(directory//'/total.endf', 1, self%total, self%awr)
    call read_xs(directory//'/elastic.endf', 2, self%elastic, other_awr)
    call require(abs(self%awr-other_awr) < 1.e-8_dp, 'inconsistent isotope masses')
    call open_input(directory//'/angular.endf', unit)
    call record(unit, 4, 2, v)
    call require(nint(v(1)) == 92238 .and. abs(v(2)-self%awr) < 1.e-8_dp, 'wrong angular isotope')
    call require(nint(v(3)) == 0 .and. nint(v(4)) == 1, 'requires MF4 LVT=0, LTT=1 (Legendre)')
    call record(unit, 4, 2, v)
    call require(nint(v(3)) == 0 .and. nint(v(4)) == 2, 'requires anisotropic center-of-mass data')
    call record(unit, 4, 2, v)
    call require(nint(v(5)) == 1, 'only one angular energy interpolation region is supported')
    n = nint(v(6))
    call require(n >= 2 .and. n <= 100000, 'invalid angular energy count')
    call values_array(unit, 4, 2, interpolation)
    call require(nint(interpolation(1)) == n .and. nint(interpolation(2)) == 2, &
                 'angular coefficients require linear-linear interpolation')
    if (allocated(self%angular_energy)) deallocate(self%angular_energy, self%coefficient)
    allocate(self%angular_energy(n), self%coefficient(max_order,n))
    self%coefficient = 0
    do i = 1, n
      call record(unit, 4, 2, v)
      self%angular_energy(i) = v(2)
      order = nint(v(5))
      call require(order >= 0 .and. order <= max_order, 'unsupported Legendre order')
      call values_array(unit, 4, 2, self%coefficient(:order,i))
    end do
    close(unit)
    call require(all(self%angular_energy(2:) > self%angular_energy(:n-1)), 'angular energies must increase')
    call require(self%angular_energy(1) <= supported_min .and. self%angular_energy(n) >= supported_max, &
                 'angular data do not cover supported energies')
    ! Check every contributing node, including brackets of the supported range.
    ! Interpolation is a convex mixture of these densities. Never clip negatives.
    do i = 1, n
      if (i < n) then
        if (self%angular_energy(i+1) < supported_min) cycle
      end if
      if (i > 1) then
        if (self%angular_energy(i-1) >= supported_max) cycle
      end if
      energy = min(supported_max, max(supported_min, self%angular_energy(i)))
      call coefficients(self, energy, a)
      do j = 0, 8192
        mu = -1.0_dp + 2.0_dp*real(j,dp)/8192.0_dp
        call require(pdf_from_coefficients(a,mu) >= 0, 'nonpositive angular data in supported range')
      end do
    end do
  end subroutine

  subroutine check_energy(self, energy)
    class(collision_data), intent(in) :: self
    real(dp), intent(in) :: energy
    call require(allocated(self%angular_energy), 'load data before use')
    call require(ieee_is_finite(energy), 'nonfinite energy')
    call require(energy >= supported_min .and. energy <= supported_max, &
                 'energy outside supported interval [150000,19000000] eV')
  end subroutine

  integer function lower_index(grid, x) result(lo)
    real(dp), intent(in) :: grid(:), x
    integer :: hi, mid
    lo = 1
    hi = size(grid)
    do while (hi-lo > 1)
      mid = (lo+hi)/2
      if (x >= grid(mid)) then
        lo = mid
      else
        hi = mid
      end if
    end do
  end function

  function interpolate_xs(table, energy) result(sigma)
    type(cross_section_table), intent(in) :: table
    real(dp), intent(in) :: energy
    real(dp) :: sigma, t
    integer :: i
    i = lower_index(table%energy, energy)
    t = (energy-table%energy(i))/(table%energy(i+1)-table%energy(i))
    sigma = (1-t)*table%sigma(i) + t*table%sigma(i+1)
  end function

  subroutine cross_sections(self, energy, sigma_total, sigma_elastic)
    class(collision_data), intent(in) :: self
    real(dp), intent(in) :: energy
    real(dp), intent(out) :: sigma_total, sigma_elastic
    call check_energy(self, energy)
    sigma_total = interpolate_xs(self%total, energy)
    sigma_elastic = interpolate_xs(self%elastic, energy)
    call require(sigma_total > 0 .and. sigma_elastic <= sigma_total, 'inconsistent cross sections')
  end subroutine

  subroutine coefficients(self, energy, a)
    class(collision_data), intent(in) :: self
    real(dp), intent(in) :: energy
    real(dp), intent(out) :: a(max_order)
    real(dp) :: t
    integer :: i
    call check_energy(self, energy)
    i = lower_index(self%angular_energy, energy)
    t = (energy-self%angular_energy(i))/(self%angular_energy(i+1)-self%angular_energy(i))
    a = (1-t)*self%coefficient(:,i) + t*self%coefficient(:,i+1)
  end subroutine

  function legendre(mu) result(p)
    real(dp), intent(in) :: mu
    real(dp) :: p(0:max_order+1)
    integer :: l
    p(0) = 1
    p(1) = mu
    do l = 2, max_order+1
      p(l) = ((2*l-1)*mu*p(l-1)-(l-1)*p(l-2))/l
    end do
  end function

  function pdf_from_coefficients(a, mu) result(pdf)
    real(dp), intent(in) :: a(max_order), mu
    real(dp) :: pdf, p(0:max_order+1)
    integer :: l
    p = legendre(mu)
    pdf = 0.5_dp
    do l = 1, max_order
      pdf = pdf + 0.5_dp*(2*l+1)*a(l)*p(l)
    end do
  end function

  function cdf_from_coefficients(a, mu) result(cdf)
    real(dp), intent(in) :: a(max_order), mu
    real(dp) :: cdf, p(0:max_order+1)
    integer :: l
    p = legendre(mu)
    cdf = 0.5_dp*(mu+1)
    do l = 1, max_order
      cdf = cdf + 0.5_dp*a(l)*(p(l+1)-p(l-1))
    end do
  end function

  function angular_pdf(self, energy, mu) result(pdf)
    class(collision_data), intent(in) :: self
    real(dp), intent(in) :: energy, mu
    real(dp) :: pdf, a(max_order)
    call require(ieee_is_finite(mu), 'nonfinite cosine')
    call require(abs(mu) <= 1, 'cosine outside [-1,1]')
    call coefficients(self, energy, a)
    pdf = pdf_from_coefficients(a, mu)
    call require(pdf >= 0, 'negative angular density')
  end function

  function angular_cdf(self, energy, mu) result(cdf)
    class(collision_data), intent(in) :: self
    real(dp), intent(in) :: energy, mu
    real(dp) :: cdf, a(max_order)
    call require(ieee_is_finite(mu), 'nonfinite cosine')
    call require(abs(mu) <= 1, 'cosine outside [-1,1]')
    call coefficients(self, energy, a)
    cdf = cdf_from_coefficients(a, mu)
  end function

  ! Center-of-mass d sigma_el / d Omega, in barn/sr.
  function differential_cross_section(self, energy, mu) result(dsdo)
    class(collision_data), intent(in) :: self
    real(dp), intent(in) :: energy, mu
    real(dp) :: dsdo, sigma_total, sigma_elastic
    call self%cross_sections(energy, sigma_total, sigma_elastic)
    dsdo = sigma_elastic*self%angular_pdf(energy,mu)/(2*pi)
  end function

  function mean_cosine_cm(self, energy) result(g)
    class(collision_data), intent(in) :: self
    real(dp), intent(in) :: energy
    real(dp) :: g, a(max_order)
    call coefficients(self, energy, a)
    g = a(1)
  end function

  function mass_ratio(self) result(awr)
    class(collision_data), intent(in) :: self
    real(dp) :: awr
    call require(self%awr > 0, 'load data before use')
    awr = self%awr
  end function

  ! Uniforms may be supplied by a caller's RNG; otherwise use RANDOM_NUMBER.
  ! Output energy is eV, output direction is in the laboratory frame.
  subroutine scatter(self, energy_in, direction_in, energy_out, direction_out, mu_cm, uniforms)
    class(collision_data), intent(in) :: self
    real(dp), intent(in) :: energy_in, direction_in(3)
    real(dp), intent(out) :: energy_out, direction_out(3)
    real(dp), intent(out), optional :: mu_cm
    real(dp), intent(in), optional :: uniforms(2)
    real(dp) :: a(max_order), u(2), lo, hi, mu, mu_lab, phi, sine, denominator
    real(dp) :: w(3), tangent(3), bitangent(3), norm, axis(3)
    integer :: j, k
    call coefficients(self, energy_in, a)
    call require(all(ieee_is_finite(direction_in)), 'nonfinite direction')
    call require(maxval(abs(direction_in)) <= 1.000001_dp, 'input direction must be unit length')
    norm = sqrt(sum(direction_in**2))
    call require(abs(norm-1) < 1.e-10_dp, 'input direction must be unit length')
    w = direction_in/norm
    if (present(uniforms)) then
      u = uniforms
    else
      call random_number(u)
    end if
    call require(all(ieee_is_finite(u)), 'nonfinite random variate')
    call require(all(u >= 0) .and. all(u < 1), 'uniforms must be in [0,1)')
    ! Exact polynomial CDF, no angular mesh or histogram approximation.
    lo = -1
    hi = 1
    do j = 1, 50
      mu = 0.5_dp*(lo+hi)
      if (cdf_from_coefficients(a,mu) < u(1)) then
        lo = mu
      else
        hi = mu
      end if
    end do
    mu = 0.5_dp*(lo+hi)
    call require(pdf_from_coefficients(a,mu) >= 0, 'invalid angular density at sampled angle')
    denominator = sqrt(self%awr**2 + 1 + 2*self%awr*mu)
    energy_out = energy_in*(denominator/(self%awr+1))**2
    mu_lab = (1+self%awr*mu)/denominator
    mu_lab = max(-1.0_dp,min(1.0_dp,mu_lab)) ! Roundoff guard on a cosine only.
    sine = sqrt(max(0.0_dp,1-mu_lab**2))
    phi = 2*pi*u(2)
    ! Cross with the least-aligned Cartesian axis: stable also at both poles.
    k = minloc(abs(w), dim=1)
    axis = 0
    axis(k) = 1
    tangent = cross(axis,w)
    tangent = tangent/sqrt(sum(tangent**2))
    bitangent = cross(w,tangent)
    direction_out = mu_lab*w + sine*(cos(phi)*tangent + sin(phi)*bitangent)
    direction_out = direction_out/sqrt(sum(direction_out**2))
    if (present(mu_cm)) mu_cm = mu
  end subroutine

  function cross(a,b) result(c)
    real(dp), intent(in) :: a(3), b(3)
    real(dp) :: c(3)
    c = [a(2)*b(3)-a(3)*b(2), a(3)*b(1)-a(1)*b(3), a(1)*b(2)-a(2)*b(1)]
  end function
end module
