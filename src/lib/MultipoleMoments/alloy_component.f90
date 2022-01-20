module alloy_component_mod

   implicit none (type, external)
   private
   ! Module variables

   ! Public members
   integer, parameter, public :: N_MAX_POLES = 4
   integer, parameter, public :: MAX_DOS_MOMENT = 3

   public :: create_alloy_component
   public :: init_partial_waves
   public :: init_component_gf
   public :: init_multipole_moments
   public :: update_partial_waves
   public :: update_potential_parameters
   public :: find_potf_zeros
   public :: update_charge_density
   public :: update_multipole_moments
   public :: recalculate_core
   public :: LsfEquation
   public :: is_lsf_mode


contains

   !-----------------------------------------------------------------------------------
   !> Initialize multipole moments.
   !-----------------------------------------------------------------------------------
   subroutine init_multipole_moments(lmax_qlm, comp)
      !> `l`-cutoff for multipole moments
      integer, intent(in) :: lmax_qlm

      !> Alloy component object
      type(alloy_component), intent(inout) :: comp

      integer :: nlm

      comp%lmax_qlm = lmax_qlm

      nlm = (lmax_qlm + 1)**2

      allocate(comp%Qlm(nlm), source = 0.0_dp)
      allocate(comp%Qlm_new(nlm), source = 0.0_dp)
   end subroutine init_multipole_moments

   !-----------------------------------------------------------------------------------
   !> Calculate a set of partial waves for a given iz, isp.
   !! Logarithmic derivatives are also updated here.
   !-----------------------------------------------------------------------------------
   subroutine update_partial_waves(z_mesh, iz, isp, comp, rad_solver)
      use energy_mesh_mod, only: energy_mesh
      use radial_wavefunction_mod, only: c_radial_wavefunction, &
         create_radial_wavefunction, &
         free_radial_wavefunction
      use radial_function_mod, only: r_radial_function, &
         create_radial_function, &
         copy_radial_function, &
         free_radial_function, &
         integrate_square
      use rad_solver_base_mod, only: IRadSolver
      use debug_dump_info, only: dump_object
      type(energy_mesh), intent(in) :: z_mesh
      integer, intent(in) :: iz
      integer, intent(in) :: isp
      type(alloy_component), intent(inout), target :: comp
      class(IRadSolver), intent(in) :: rad_solver

      integer :: ir_max, lmax
      complex(kind = dp) :: z
      real(kind = dp) :: beta
      type(c_lm_radial_function), pointer :: p_wave
      type(c_radial_wavefunction) :: wf
      type(r_radial_function) :: rpot
      real(kind = dp) :: sgn
      complex(kind = dp) :: ctmp
      integer :: kappa, l, ierr, ir_ws
      integer :: ilm, ilm1, ilm2, ilm0, ilmp

      ASSERT(0 < isp .and. isp <= size(comp%pot%val, 2), &
      "Spin index is outside of the allowed range")
      ASSERT(0 < iz .and. iz <= z_mesh%nz_tot, &
      "Energy index is outside of the allowed range")

      ASSERT(allocated(comp%p_waves), &
      "Partial waves are not initialized")

      !--- Allocate a new partial wave if necessary.
      !    This is done to ensure that the functions are allocated
      !    only for (i_at, iz) actually used in calculations.
      if (.not. allocated(comp % p_waves(iz, isp) % val)) then
         call create_radial_function(comp % r_mesh, comp % lmax, &
            comp % p_waves(iz, isp))
      endif

      p_wave => comp%p_waves(iz, isp)

      call copy_radial_function(comp%pot, rpot)
      call create_radial_wavefunction(comp%r_mesh, wf)

      z = z_mesh%z(iz)

      ir_max = comp%ir_chd
      ir_ws = comp%ir_ps
      lmax = p_wave%lmax

      ! TODO: Define `radial_function` with a pointer array and make
      !       a view construction
      !--- Copy potential to a temporary function
      rpot%val = comp%pot%val(:, isp)
      rpot%ir_max = comp%pot%ir_max

      do l = 0, lmax
         ilm1 = l**2 + 1
         ilm2 = (l + 1)**2

         !--- Solve Dirac equation
         kappa = -l - 1

         call rad_solver%c_eval(z, kappa, ir_max, rpot, beta, wf, ierr)


         !--- Evaluate the log-derivative
         comp%log_ders(ilm1, iz, isp) = rad_solver%log_derivative(wf, ir_ws, l, rpot%val(ir_ws))

         do ilm = ilm1 + 1, ilm2
            comp%log_ders(ilm, iz, isp) = comp%log_ders(ilm1, iz, isp)
         enddo

         call rad_solver%normalize(wf, ir_ws, beta, &
            sign_correction = .True.)

         !
         !--- Distribute `l`-components over `lm`-function
         !
         ilm0 = p_wave%lm_map(ilm1)
         p_wave%val(:, ilm0) = wf%p%val

         do ilm = ilm1 + 1, ilm2
            ilmp = p_wave%lm_map(ilm)
            if(ilmp /= ilm0) then
               p_wave%val(:, ilmp) = p_wave%val(:, ilm0)
            endif
         enddo
      end do

      call free_radial_wavefunction(wf)
      call free_radial_function(rpot)
   end subroutine update_partial_waves

   !-----------------------------------------------------------------------------------
   !> Calculate potential functions and parameters for a given iz, isp.
   !-----------------------------------------------------------------------------------
   subroutine update_potential_parameters(iz, isp, ws_rad, scr_alpha, comp)
      use rad_solver_base_mod, only: IRadSolver

      !> Energy point index
      integer, intent(in) :: iz

      !> Spin channel
      integer, intent(in) :: isp

      !> Average WS radius (used for scaling)
      real(kind = dp), intent(in) :: ws_rad

      !> Screening parameters
      real(kind = dp), intent(in) :: scr_alpha(:)

      !> Alloy component
      type(alloy_component), intent(inout) :: comp

      integer :: l, nlm, ilm, ilm1, ilm2, ilm_w
      real(kind = dp) :: ws_ratio, ws_ratio_12, fac
      complex(kind = dp) :: logder, p_ws, potf, pdot12

      nlm = (comp%p_waves(iz, isp)%lmax + 1)**2
      ASSERT(size(scr_alpha, 1) >= nlm, &
      "Screening-parameter array is too small for the site `l`-cutoff")

      do l = 0, comp%p_waves(iz, isp)%lmax
         ilm1 = l**2 + 1
         ilm2 = (l + 1)**2

         ws_ratio = (ws_rad / comp%ps_rad)**(2 * l + 1)
         ws_ratio_12 = sqrt(ws_ratio)

         ! Here we assume that `P(z)` depends only on `l`!
         ! Logarithmic derivative at the potential sphere
         logder = comp%log_ders(ilm1, iz, isp)

         ! Partial wave at the potential sphere
         ilm_w = comp%p_waves(iz, isp)%lm_map(ilm1)
         p_ws = comp%p_waves(iz, isp)%val(comp%ir_ps, ilm_w)

         ! Potential function and its energy derivative
         potf = 2 * (2 * l + 1) * (logder + l + 1) / (logder - l) * ws_ratio

         fac = (2 * l + 1.0_dp) * sqrt(2.0_dp / comp%ps_rad) * comp%ps_rad * ws_ratio_12
         pdot12 = fac / p_ws / (logder - l)

         ! Screening transformation
         do ilm = ilm1, ilm2
            comp%potf(ilm, iz, isp) = potf / (1.0_dp - scr_alpha(ilm) * potf)
            comp%pdot12(ilm, iz, isp) = pdot12 / (1.0_dp - scr_alpha(ilm) * potf)
         enddo

         ! TODO: add a tail to the partial wave for `r > r_ws`.
         !       The tail is a properly normalized linear combination of functions K(r), J(r).
      enddo

   end subroutine update_potential_parameters

   !-----------------------------------------------------------------------------------
   !> Calculate potential functions and parameters for a given iz, isp.
   !-----------------------------------------------------------------------------------
   subroutine find_potf_zeros(isp, comp, ierr, e_cut, rad_solver)
      use constants, only: pi
      use radial_function_mod, only: r_radial_function, &
         create_radial_function, &
         free_radial_function, &
         integrate_square
      use radial_wavefunction_mod, only: r_radial_wavefunction, &
         create_radial_wavefunction, &
         free_radial_wavefunction
      use gf_poles_mod, only: gf_poles, EtaLogder

      use zero_finding_mod, only: solve_univariate
      use rad_solver_base_mod, only: IRadSolver
      use utils_mod, only: option

      !> Spin channel
      integer, intent(in) :: isp

      !> Alloy component
      type(alloy_component), intent(inout), target :: comp

      !> Error code
      integer, intent(out) :: ierr

      !> Cut-off energy after which the poles are discarded
      real(kind = dp), intent(in), optional :: e_cut
      class(IRadSolver), intent(in) :: rad_solver !< Radial solver

      real(kind = dp), parameter :: e_step = 0.5_dp
      real(kind = dp), parameter :: eta_tol = 1e-14_dp
      integer, parameter :: n_iter_max = 30

      type(gf_poles), pointer :: poles
      integer :: l_qn, l_prev, il, ilm1, ilm, nlm
      integer :: ip, ierr2
      integer :: k_qn
      type(r_radial_function) :: rpot
      type(r_radial_wavefunction) :: rad_wf
      real(kind = dp) :: e_nu, d_nu, eta_nu, e_cut_
      real(kind = dp) :: wf_norm, beta
      logical :: skip

      type(EtaLogder) :: eta_logder

      e_cut_ = option(100.0_dp, e_cut)

      ierr = 0
      poles => comp%poles(isp)

      call create_radial_function(comp%r_mesh, rpot)
      call create_radial_wavefunction(comp%r_mesh, rad_wf)

      rpot%val = comp%pot%val(:, isp)
      rpot%ir_max = comp%pot%ir_max

      nlm = size(poles%n_poles, 1)

      skip = .False.

      call eta_logder%init(0, comp%ir_ps, comp%ir_chd, rpot, rad_solver)

      ! TODO: Additional poles
      do ip = 1, 2

         ilm1 = 1
         l_prev = -1
         do ilm = 1, nlm
            l_qn = poles%l_poles(ilm)
            il = l_qn + 1

            if(l_prev /= l_qn) then
               l_prev = l_qn
               ilm1 = ilm

               call eta_logder%update(l_qn)

               !
               !--- Find such `e` that D(e) = -l - 1
               !
               d_nu = real(-l_qn - 1, kind = dp)
               eta_nu = 0.5_dp - atan(d_nu) / pi + &
                  real(comp%species%wf_nodes(il) + ip - 1, kind = dp)

               call solve_univariate(eta_logder, eta_nu, &
                  poles%e_poles(ilm, ip), e_step, eta_tol, &
                  n_iter_max, e_nu, ierr2)

               !
               !--- If not converged: another attempt with a different initial value
               !
               if(ierr2 == 1) then
                  call solve_univariate(eta_logder, eta_nu, &
                     0.0_dp, e_step, eta_tol, n_iter_max, e_nu, ierr2)

                  if(ierr2 == 1) then
                     ierr = il ! Fatal: Return `l` for which the pole is not found
                     return
                  endif
               endif

               if (ip == 1 .or. (ip > 1 .and. e_nu < e_cut_)) then
                  skip = .False.

                  poles%n_poles(ilm) = ip
                  poles%e_poles(ilm, ip) = e_nu
                  poles%r_poles(ilm, ip) = 1.0_dp

                  !
                  !--- Find the wavefunction
                  !
                  k_qn = -l_qn - 1

                  call rad_solver%eval(e_nu, k_qn, comp%ir_chd, rpot, beta, rad_wf, ierr2)

                  call rad_solver%normalize(rad_wf, comp%ir_ps, beta)

                  ! Store normalized WF
                  poles%wf_poles(ilm, ip)%val = rad_wf%p%val
               else
                  skip = .True.
               endif
            else
               !
               !--- Copy from `ilm1`
               !
               if (.not. skip) then
                  poles%n_poles(ilm) = poles%n_poles(ilm1)
                  poles%r_poles(ilm, ip) = poles%r_poles(ilm1, ip)
                  poles%e_poles(ilm, ip) = poles%e_poles(ilm1, ip)
                  poles%wf_poles(ilm, ip)%val = poles%wf_poles(ilm1, ip)%val
               endif
            endif
         enddo
      enddo

      call free_radial_wavefunction(rad_wf)
      call free_radial_function(rpot)
   end subroutine find_potf_zeros

   !-----------------------------------------------------------------------------------
   !> Update the charge density and the charge transfer.
   !-----------------------------------------------------------------------------------
   subroutine update_charge_density(comp, chd_schema)
      use constants, only: pi
      use chd_schema_mod, only: ChdSchema
      use energy_mesh_mod, only: energy_mesh
      use utils_mod, only: spin_factor

      !> Alloy component
      type(alloy_component), intent(inout), target :: comp

      !> Charge-density schema
      type(ChdSchema), intent(inout) :: chd_schema

      integer :: isp, nsp
      type(energy_mesh), pointer :: z_mesh
      real(kind = dp) :: sp_fac
      complex(kind = dp) :: z, z_wght
      integer :: l, lmax
      integer :: ilm, ilm_w, lm_i, lm_f
      integer :: ip, ir, iz
      complex(kind = dp) :: gf_rr, pw
      complex(kind = dp) :: pole_wght(N_MAX_POLES)
      real(kind = dp) :: phi, phidot, phi2dot, rmom(0:2)

      comp % chd_new % val = 0.0_dp

      nsp = comp%chd_new%nsp

      !
      !--- Set spin factor
      !
      sp_fac = spin_factor(nsp)

      !--- `l`-cutoff
      lmax = comp % lmax

      pole_wght = 0.0_dp

      !
      !--- Calculate the charge density
      !
      do isp = 1, nsp
         if (.not. chd_schema % is_spin_domain(isp)) cycle

         !
         !--- Real-axis contribution.
         !    We place it before the main contour to be able to broadcast
         !    the part for `proc_map(1, i_at)` immediately.
         !
         !    Note that since LMTO functions depend only on `l` we can
         !    sum moments over `m`.
         !
         !    Additional contributions (`aux_mom` and `core_chd`) are
         !    added only on processes associated with `iz = 1`.
         !
         if (chd_schema % atom_domain % proc_map(1, comp % i_at) == &
            chd_schema % atom_domain % rank) then

            do l = 0, lmax
               lm_i = l**2 + 1
               lm_f = (l + 1)**2

               rmom(0:2) = sum(comp%aux_mom(0:2, lm_i:lm_f, isp), 2) * sp_fac

               do ir = 1, comp%ir_chd
                  phi = comp%lmto_pars(isp)%phis(0, l)%val(ir)
                  phidot = comp%lmto_pars(isp)%phis(1, l)%val(ir)
                  phi2dot = comp%lmto_pars(isp)%phis(2, l)%val(ir)

                  comp%chd_new%val(ir, isp) = comp%chd_new%val(ir, isp) + &
                     rmom(0) * phi**2 + &
                     rmom(1) * 2.0_dp * phi * phidot + &
                     rmom(2) * (phidot**2 + phi * phi2dot)
               enddo
            enddo ! l

            !
            !--- Add core density
            !
            comp%chd_new%val(:, isp) = comp%chd_new%val(:, isp) + comp%core_chd%val(:, isp)
         endif

         z_mesh => comp%gf(isp)%z_mesh

         !
         !--- Main-contour contribution
         !
         do iz = 1, z_mesh%nzm
            if (.not. chd_schema % atom_domain % is_in_domain(iz, comp % i_at)) cycle

            z_wght = z_mesh%weights(iz) * sp_fac / pi
            z = z_mesh%z(iz)

            do ilm = 1, (lmax + 1)**2
               ilm_w = comp%p_waves(iz, isp)%lm_map(ilm)

               !
               !--- Pre-calculate pole weights
               !
               do ip = 1, comp%poles(isp)%n_poles(ilm)
                  pole_wght(ip) = 1.0_dp / &
                     (z - comp%poles(isp)%e_poles(ilm, ip))
               enddo

               do ir = 1, comp%ir_chd
                  !--- Diagonal GF
                  pw = comp%p_waves(iz, isp)%val(ir, ilm_w)
                  gf_rr = comp%gf(isp)%val(ilm, ilm, iz) * pw**2

                  !--- Pole contributions
                  do ip = 1, comp%poles(isp)%n_poles(ilm)
                     phi = comp%poles(isp)%wf_poles(ilm, ip)%val(ir)
                     gf_rr = gf_rr - phi**2 * pole_wght(ip)
                  enddo

                  comp%chd_new%val(ir, isp) = comp%chd_new%val(ir, isp) - &
                     aimag(z_wght * gf_rr)
               enddo ! ir
            enddo ! ilm
         enddo ! iz

         call chd_schema % reduce_z(comp % i_at, isp, comp % chd_new % val)
      enddo ! isp

      !
      !--- Update the charge transfer
      !
      comp%qtr_new = sum(comp%dos_mom(0, :, :)) * sp_fac - comp%nval
   end subroutine update_charge_density

   !-----------------------------------------------------------------------------------
   !> Update the multipole moments.
   !!
   !! Output values `Qlm_new` are meaningful only for `l > 0`.
   !! The monopole term is taken from `Qtr` calculated in `update_charge_density`.
   !-----------------------------------------------------------------------------------
   subroutine update_multipole_moments(ws_rad, comp, qlm_schema)
      use constants, only: pi
      use qlm_schema_mod, only: QlmSchema
      use radial_function_mod, only: &
         r_radial_function, &
         c_radial_function, &
         create_radial_function, &
         free_radial_function, &
         integrate
      use energy_mesh_mod, only: energy_mesh
      use spherical_harmonics_mod, only: gaunt_real_table
      use utils_mod, only: spin_factor

      !> WS radius, needed for proper scaling
      real(kind = dp), intent(in) :: ws_rad

      !> Alloy component object
      type(alloy_component), intent(inout) :: comp

      !> MPI schema for multipole moments
      type(QlmSchema), intent(inout) :: qlm_schema

      type(energy_mesh), pointer :: z_mesh
      complex(kind = dp), allocatable :: q_ll(:, :, :) ! Main contribution
      real(kind = dp), allocatable :: qp_ll(:, :) ! Pole contributions
      type(c_radial_function) :: temp
      type(r_radial_function) :: rtmp, rl_w

      integer :: nlm, nz, nsp, n_poles, lmax
      integer :: iz, isp, ip
      integer :: l, m, lm, l1, m1, lm1, l2, m2, lm2
      integer :: ilm_w1, ilm_w2
      real(kind = dp) :: fac, sp_fac
      real(kind = dp) :: gaunt
      complex(kind = dp) :: z, z_wght, gf_ll
      complex(kind = dp) :: cg_ff


      ASSERT(allocated(comp%Qlm), "Multipole moments are not initialized")

      !
      !--- Quick return if no multipole moments are considered
      !
      if(comp%lmax_qlm == 0) then
         comp%Qlm(1) = 0.0_dp
         return
      endif

      lmax = comp % lmax
      nlm = (lmax + 1)**2
      nz = size(comp%p_waves, 1)
      nsp = size(comp%p_waves, 2)

      sp_fac = spin_factor(nsp)

      n_poles = 0
      do isp = 1, nsp
         n_poles = max(n_poles, maxval(comp%poles(isp)%n_poles))
      enddo

      allocate(qp_ll(nlm, n_poles), source = 0.0_dp)
      allocate(q_ll(nlm, nlm, nz), source = (0.0_dp, 0.0_dp))

      call create_radial_function(comp%r_mesh, temp)
      call create_radial_function(comp%r_mesh, rl_w)
      call create_radial_function(comp%r_mesh, rtmp)

      comp%Qlm = 0.0_dp

      sp: do isp = 1, nsp
         if (.not. qlm_schema % atom_domain % is_spin_domain(isp)) cycle

         z_mesh => comp%gf(isp)%z_mesh

         lqn: do l = 0, comp%lmax_qlm
            !
            !--- (r / w)^{l}
            !
            !    r == r_{R}
            !
            rl_w%val = comp%r_mesh%r**l / ws_rad**l

            fac = sqrt(4 * pi) / (2 * l + 1)

            !
            !--- Partial-wave contributions from the poles
            !
            !    q^{p}_L1 = \delta_{L1, L2} \int_{0}^{s_R} dr (r / w)^{l} * phi_L1(r, e_p)**2
            !
            do lm1 = 1, nlm
               do ip = 1, comp%poles(isp)%n_poles(lm1)
                  rtmp%val = comp%poles(isp)%wf_poles(lm1, ip)%val * &
                     comp%poles(isp)%wf_poles(lm1, ip)%val * rl_w%val

                  qp_ll(lm1, ip) = integrate(rtmp, comp%ir_ps)
               enddo
            enddo

            !
            !--- Partial-wave contributions from the contour
            !
            !    q_{L1, L2}(z) = \int_{0}^{s_R} dr (r / w)^{l} *
            !                        phi_L1(r, z) * phi_L2(r, z)
            !
            do iz = 1, z_mesh%nzm
               if (.not. qlm_schema % atom_domain % is_in_domain(iz, comp % i_at)) cycle

               do lm1 = 1, nlm
                  do lm2 = 1, nlm
                     ilm_w1 = comp%p_waves(iz, isp)%lm_map(lm1)
                     ilm_w2 = comp%p_waves(iz, isp)%lm_map(lm2)

                     temp%val = comp%p_waves(iz, isp)%val(:, ilm_w1) * &
                        comp%p_waves(iz, isp)%val(:, ilm_w2) * &
                        rl_w%val

                     q_ll(lm1, lm2, iz) = integrate(temp, comp%ir_ps)
                  enddo
               enddo
            enddo

            mqn: do m = -l, l
               lm = l**2 + l + m + 1

               !
               !--- C(L, L1, L2) * \int_C dz Im{ Tr [phi_L1(z) * G(z) * phi_L2(z)] -
               !        \delta_{L1, L2} * [phi^{p}_{L1}]^{2} / (z - e_p) }
               !
               do iz = 1, z_mesh%nzm
                  if (.not. qlm_schema % atom_domain % is_in_domain(iz, comp % i_at)) cycle

                  z_wght = z_mesh%weights(iz) * sp_fac / pi
                  z = z_mesh%z(iz)

                  cg_ff = 0.0_dp
                  L_1: do l1 = 0, lmax
                     do m1 = -l1, l1
                        lm1 = l1**2 + l1 + m1 + 1

                        L_2: do l2 = 0, lmax
                           do m2 = -l2, l2
                              lm2 = l2**2 + l2 + m2 + 1

                              gaunt = gaunt_real_table(l1, l2, l, m1, m2, m)

                              ! Skip zero values
                              if(abs(gaunt) < 1e-13_dp) then
                                 cycle
                              endif

                              !
                              !--- G_{L1, L2}(z) * q_{L1, L2}(z)
                              !
                              gf_ll = comp%gf(isp)%val(lm1, lm2, iz) * q_ll(lm1, lm2, iz)

                              !
                              !--- \delta_{L1, L2} * q^{p}_{L1} / (z - e_p)
                              !
                              if(lm1 == lm2) then
                                 do ip = 1, comp%poles(isp)%n_poles(lm1)
                                    gf_ll = gf_ll - qp_ll(lm1, ip) / &
                                       (z - comp%poles(isp)%e_poles(lm1, ip))
                                 enddo
                              endif

                              cg_ff = cg_ff - gaunt * gf_ll
                           enddo  ! m2
                        enddo  L_2
                     enddo ! m1
                  enddo L_1

                  comp%Qlm(lm) = comp%Qlm(lm) + aimag(z_wght * cg_ff) * fac
               enddo ! iz

            enddo mqn
         enddo lqn

         if (qlm_schema % atom_domain % spin_domain % is_decomposed()) then
            call qlm_schema % reduce_z(comp % i_at, isp, comp % Qlm)
         endif
      enddo sp

      !--- No spin decomposition => Qlm is already summed over spins,
      !    sum-reduce over z-points only once
      if (.not. qlm_schema % atom_domain % spin_domain % is_decomposed()) then
         call qlm_schema % reduce_z(comp % i_at, 1, comp % Qlm)
      endif

      !
      !--- NB!!!!
      !    Real-axis contributions (from the auxiliary contour)
      !    are omitted because they do not matter once the system is
      !    close to convergence.
      !

      call free_radial_function(rtmp)
      call free_radial_function(rl_w)
      call free_radial_function(temp)
      deallocate(qp_ll)
   end subroutine update_multipole_moments


end module alloy_component_mod

