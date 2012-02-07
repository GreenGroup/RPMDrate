!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   RPMD - Ring polymer molecular dynamics simulations
!
!   Copyright (c) 2012 by Joshua W. Allen (jwallen@mit.edu)
!                         Yury V. Suleimanov (ysuleyma@mit.edu)
!                         William H. Green (whgreen@mit.edu)
!
!   Permission is hereby granted, free of charge, to any person obtaining a
!   copy of this software and associated documentation files (the "Software"),
!   to deal in the Software without restriction, including without limitation
!   the rights to use, copy, modify, merge, publish, distribute, sublicense,
!   and/or sell copies of the Software, and to permit persons to whom the
!   Software is furnished to do so, subject to the following conditions:
!
!   The above copyright notice and this permission notice shall be included in
!   all copies or substantial portions of the Software.
!
!   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
!   IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
!   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
!   THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
!   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
!   FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
!   DEALINGS IN THE SOFTWARE.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module system

    implicit none

    integer, parameter :: MAX_ATOMS = 100
    double precision :: dt
    double precision :: beta
    double precision :: mass(MAX_ATOMS)
    double precision :: pi = dacos(-1.0d0)

contains

    ! Advance the simluation by one time step using the velocity Verlet
    ! algorithm.
    ! Parameters:
    !   t - The current simulation time
    !   p - The momentum of each bead in each atom
    !   q - The position of each bead in each atom
    !   V - The potential of each bead
    !   dVdq - The force exerted on each bead in each atom
    !   Natoms - The number of atoms in the molecular system
    !   Nbeads - The number of beads to use per atom
    !   potential - A function that evaluates the potential and force for a given position
    ! Returns:
    !   t - The updated simulation time
    !   p - The updated momentum of each bead in each atom
    !   q - The updated position of each bead in each atom
    !   V - The updated potential of each bead
    !   dVdq - The updated force exerted on each bead in each atom
    !   result - A flag that indicates if the time step completed successfully (if zero) or that an error occurred (if nonzero)
    subroutine verlet_step(t, p, q, V, dVdq, Natoms, Nbeads, potential, result)

        implicit none
        external potential
        integer, intent(in) :: Natoms, Nbeads
        double precision, intent(inout) :: t, p(3,Natoms,Nbeads), q(3,Natoms,Nbeads)
        double precision, intent(inout) :: V(Nbeads), dVdq(3,Natoms,Nbeads)
        integer, intent(out) :: result

        integer :: i, j

        result = 0

        ! Update momentum (half time step)
        p = p - 0.5d0 * dt * dVdq

        ! Update position (full time step)
        if (Nbeads .eq. 1) then
            ! For a single bead, there are no free ring polymer terms to add,
            ! so we simply update the positions using the momentum, as in
            ! classical trajectories
            do i = 1, 3
                do j = 1, Natoms
                    q(i,j,1) = q(i,j,1) + p(i,j,1) * dt / mass(j)
                end do
            end do
        else
            ! For multiple beads, we update the positions and momenta for the
            ! harmonic free ring term in the Hamiltonian by transforming to
            ! and from normal mode space
            call free_ring_polymer_step(p, q, Natoms, Nbeads)
        end if

        ! Update potential and forces using new position
        call potential(q, V, dVdq, Natoms, Nbeads)

        ! Update momentum (half time step)
        p = p - 0.5d0 * dt * dVdq

        ! Update time
        t = t + dt

    end subroutine verlet_step

    ! Update the positions and momenta of each atom in each free ring polymer
    ! bead for the term in the Hamiltonian describing the harmonic free ring
    ! polymer interactions. This is most efficiently done in normal mode space;
    ! this function therefore uses fast Fourier transforms (from the FFTW3
    ! library) to transform to and from normal mode space.
    ! Parameters:
    !   p - The momentum of each bead in each atom
    !   q - The position of each bead in each atom
    !   Natoms - The number of atoms in the molecular system
    !   Nbeads - The number of beads to use per atom
    ! Returns:
    !   p - The updated momentum of each bead in each atom
    !   q - The updated position of each bead in each atom
    subroutine free_ring_polymer_step(p, q, Natoms, Nbeads)

        implicit none
        integer, intent(in) :: Natoms, Nbeads
        double precision, intent(inout) :: p(3,Natoms,Nbeads), q(3,Natoms,Nbeads)

        double precision :: poly(4,Nbeads)
        double precision :: beta_n, twown, pi_n, wk, wt, wm, cos_wt, sin_wt, p_new
        integer :: i, j, k

        ! Transform to normal mode space
        do i = 1, 3
            do j = 1, Natoms
                call rfft(p(i,j,:), Nbeads)
                call rfft(q(i,j,:), Nbeads)
            end do
        end do

        do j = 1, Natoms

            poly(1,1) = 1.0d0
            poly(2,1) = 0.0d0
            poly(3,1) = dt / mass(j)
            poly(4,1) = 1.0d0

            if (Nbeads .gt. 1) then
                beta_n = beta / Nbeads
                twown = 2.0d0 / beta_n
                pi_n = pi / Nbeads
                do k = 1, Nbeads / 2
                    wk = twown * dsin(k * pi_n)
                    wt = wk * dt
                    wm = wk * mass(j)
                    cos_wt = dcos(wt)
                    sin_wt = dsin(wt)
                    poly(1,k+1) = cos_wt
                    poly(2,k+1) = -wm*sin_wt
                    poly(3,k+1) = sin_wt/wm
                    poly(4,k+1) = cos_wt
                end do
                do k = 1, (Nbeads - 1) / 2
                    poly(1,Nbeads-k+1) = poly(1,k+1)
                    poly(2,Nbeads-k+1) = poly(2,k+1)
                    poly(3,Nbeads-k+1) = poly(3,k+1)
                    poly(4,Nbeads-k+1) = poly(4,k+1)
                end do
            end if

            do k = 1, Nbeads
                do i = 1, 3
                    p_new = p(i,j,k) * poly(1,k) + q(i,j,k) * poly(2,k)
                    q(i,j,k) = p(i,j,k) * poly(3,k) + q(i,j,k) * poly(4,k)
                    p(i,j,k) = p_new
                end do
            end do

        end do

        ! Transform back to Cartesian space
        do i = 1, 3
            do j = 1, Natoms
                call irfft(p(i,j,:), Nbeads)
                call irfft(q(i,j,:), Nbeads)
            end do
        end do

    end subroutine free_ring_polymer_step

    ! Compute the total energy of all ring polymers in the RPMD system.
    ! Parameters:
    !   q - The position of each bead in each atom
    !   Natoms - The number of atoms in the molecular system
    !   Nbeads - The number of beads to use per atom
    ! Returns:
    !   Ering - The total energy of all ring polymers
    subroutine get_ring_polymer_energy(q, Natoms, Nbeads, Ering)

        implicit none
        integer, intent(in) :: Natoms, Nbeads
        double precision, intent(in) :: q(3,Natoms,Nbeads)
        double precision, intent(out) :: Ering

        double precision :: wn, dx, dy, dz
        integer :: j, k

        Ering = 0.0d0
        wn = Nbeads / beta
        do j = 1, Natoms
            dx = q(1,j,1) - q(1,j,Nbeads)
            dy = q(2,j,1) - q(2,j,Nbeads)
            dz = q(3,j,1) - q(3,j,Nbeads)
            Ering = Ering + 0.5d0 * mass(j) * wn * wn * (dx * dx + dy * dy + dz * dz)
            do k = 2, Nbeads
                dx = q(1,j,k-1) - q(1,j,k)
                dy = q(2,j,k-1) - q(2,j,k)
                dz = q(3,j,k-1) - q(3,j,k)
                Ering = Ering + 0.5d0 * mass(j) * wn * wn * (dx * dx + dy * dy + dz * dz)
            end do
        end do

    end subroutine get_ring_polymer_energy

    ! Compute the total kinetic energy of the RPMD system.
    ! Parameters:
    !   p - The momentum of each bead in each atom
    !   Natoms - The number of atoms in the molecular system
    !   Nbeads - The number of beads to use per atom
    ! Returns:
    !   Ek - The kinetic energy of the system
    subroutine get_kinetic_energy(p, Natoms, Nbeads, Ek)

        implicit none
        integer, intent(in) :: Natoms, Nbeads
        double precision, intent(in) :: p(3,Natoms,Nbeads)
        double precision, intent(out) :: Ek

        integer :: i, j, k

        Ek = 0.0d0
        do i = 1, 3
            do j = 1, Natoms
                do k = 1, Nbeads
                    Ek = Ek + 0.5d0 * p(i,j,k) * p(i,j,k) / mass(j)
                end do
            end do
        end do

    end subroutine get_kinetic_energy

    ! Compute the center of mass position of the RPMD system.
    ! Parameters:
    !   q - The position of each bead in each atom
    !   Natoms - The number of atoms in the molecular system
    !   Nbeads - The number of beads to use per atom
    ! Returns:
    !   cm - The center of mass of the system
    subroutine get_center_of_mass(q, Natoms, Nbeads, cm)

        implicit none
        integer, intent(in) :: Natoms, Nbeads
        double precision, intent(in) :: q(3,Natoms,Nbeads)
        double precision, intent(out) :: cm(3)

        double precision :: total_mass
        integer :: i, j, k

        cm(:) = 0.0d0
        total_mass = sum(mass)
        do i = 1, 3
            do j = 1, Natoms
                do k = 1, Nbeads
                    cm(i) = cm(i) + q(i,j,k) * mass(j)
                end do
            end do
            cm(i) = cm(i) / total_mass
        end do

    end subroutine get_center_of_mass

    ! Compute the centroid position of each atom in the RPMD system.
    ! Parameters:
    !   q - The position of each bead in each atom
    !   Natoms - The number of atoms in the molecular system
    !   Nbeads - The number of beads to use per atom
    ! Returns:
    !   centroid - The centroid of each atom
    subroutine get_centroid(q, Natoms, Nbeads, centroid)

        implicit none
        integer, intent(in) :: Natoms, Nbeads
        double precision, intent(in) :: q(3,Natoms,Nbeads)
        double precision, intent(out) :: centroid(3,Natoms)

        integer :: i, j, k

        centroid(:,:) = 0.0d0
        do i = 1, 3
            do j = 1, Natoms
                do k = 1, Nbeads
                    centroid(i,j) = centroid(i,j) + q(i,j,k)
                end do
                centroid(i,j) = centroid(i,j) / Nbeads
            end do
        end do

    end subroutine get_centroid

    ! Compute the radius of gyration of each atom in the RPMD system. This is
    ! a useful quantity to check while debugging.
    ! Parameters:
    !   q - The position of each bead in each atom
    !   Natoms - The number of atoms in the molecular system
    !   Nbeads - The number of beads to use per atom
    ! Returns:
    !   R - The radius of gyration of each atom
    subroutine get_radius_of_gyration(q, Natoms, Nbeads, R)

        implicit none
        integer, intent(in) :: Natoms, Nbeads
        double precision, intent(in) :: q(3,Natoms,Nbeads)
        double precision, intent(out) :: R(Natoms)

        double precision :: centroid(3,Natoms), dx
        integer :: i, j, k

        call get_centroid(q, Natoms, Nbeads, centroid)

        R(:) = 0.0d0
        do j = 1, Natoms
            do i = 1, 3
                do k = 1, Nbeads
                    dx = q(i,j,k) - centroid(i,j)
                    R(j) = R(j) + dx * dx
                end do
            end do
            R(j) = sqrt(R(j) / Nbeads)
        end do

    end subroutine get_radius_of_gyration

end module system
