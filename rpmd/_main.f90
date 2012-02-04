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
