!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   RPMDrate - Bimolecular reaction rates via ring polymer molecular dynamics
!
!   Copyright (c) 2012 by Joshua W. Allen (jwallen@mit.edu)
!                         William H. Green (whgreen@mit.edu)
!                         Yury V. Suleimanov (ysuleyma@mit.edu, ysuleyma@princeton.edu)
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

! The ``transition_state`` module defines a dividing surface for a transition
! state, defined by a set of forming bonds and a set of breaking bonds. The
! attributes are:
!
! ================================= ============================================
! Attribute                         Description
! ================================= ============================================
! ``number_of_transition_states``   The number of equivalent transition states that define the dividing surface
! ``number_of_forming_bonds``       The number of bonds being formed by the reaction
! ``number_of_breaking_bonds``      The number of bonds being broken by the reaction
! ``forming_bonds``                 An array listing the pairs of indices of each forming bond in each transition state
! ``forming_bond_lengths``          An array listing the lengths of each forming bond in each transition state
! ``breaking_bonds``                An array listing the pairs of indices of each breaking bond in each transition state
! ``breaking_bond_lengths``         An array listing the lengths of each breaking bond in each transition state
! ================================= ============================================
!
module transition_state

    implicit none

    integer, parameter :: MAX_ATOMS = 100, MAX_BONDS = 16, MAX_TS = 16

    integer :: number_of_transition_states
    integer :: number_of_bonds

    integer :: forming_bonds(MAX_TS,MAX_BONDS,2)
    integer :: breaking_bonds(MAX_TS,MAX_BONDS,2)
    double precision :: forming_bond_lengths(MAX_TS,MAX_BONDS)
    double precision :: breaking_bond_lengths(MAX_TS,MAX_BONDS)

contains

    ! Return the value of the dividing surface function for each of the
    ! equivalent transition states that define the dividing surface.
    ! Parameters:
    !   position - A 3 x Natoms array of atomic positions
    !   Natoms - The number of atoms
    ! Returns:
    !   values - The value of the dividing surface function for each transition state
    subroutine evaluate_all(position, Natoms, values)

        implicit none
        integer, intent(in) :: Natoms
        double precision, intent(in) :: position(3,Natoms)
        double precision, dimension(:,:), intent(out) :: values

        integer :: m, n, atom1, atom2
        double precision :: Rx, Ry, Rz, R

        values(:,:) = 0.0d0

        do n = 1, number_of_transition_states

            do m = 1, number_of_bonds
                ! Forming bond
                atom1 = forming_bonds(n,m,1)
                atom2 = forming_bonds(n,m,2)
                Rx = position(1,atom1) - position(1,atom2)
                Ry = position(2,atom1) - position(2,atom2)
                Rz = position(3,atom1) - position(3,atom2)
                R = sqrt(Rx * Rx + Ry * Ry + Rz * Rz)
                values(n,m) = values(n,m) + (forming_bond_lengths(n,m) - R)
                ! Breaking bond
                atom1 = breaking_bonds(n,m,1)
                atom2 = breaking_bonds(n,m,2)
                Rx = position(1,atom1) - position(1,atom2)
                Ry = position(2,atom1) - position(2,atom2)
                Rz = position(3,atom1) - position(3,atom2)
                R = sqrt(Rx * Rx + Ry * Ry + Rz * Rz)
                values(n,m) = values(n,m) - (breaking_bond_lengths(n,m) - R)
            end do

        end do

    end subroutine evaluate_all

    ! Return the value of the transition state dividing surface function. This
    ! is the maximum of the individual values for each equivalent transition
    ! state, as determined by evaluate_all().
    ! Parameters:
    !   position - A 3 x Natoms array of atomic positions
    !   Natoms - The number of atoms
    ! Returns:
    !   s1 - The value of the dividing surface function
    subroutine value(position, Natoms, s1)

        implicit none
        integer, intent(in) :: Natoms
        double precision, intent(in) :: position(3,Natoms)
        double precision, intent(out) :: s1

        double precision :: values(number_of_transition_states,number_of_bonds)

        call evaluate_all(position, Natoms, values)

        s1 = minval(maxval(values, 1))

    end subroutine value

    ! Return the value of the gradient of the transition state dividing surface
    ! function.
    ! Parameters:
    !   position - A 3 x Natoms array of atomic positions
    !   Natoms - The number of atoms
    ! Returns:
    !   ds1 - The gradient of the dividing surface function
    subroutine gradient(position, Natoms, ds1)

        implicit none
        integer, intent(in) :: Natoms
        double precision, intent(in) :: position(3,Natoms)
        double precision, intent(out) :: ds1(3,Natoms)

        double precision :: values(number_of_transition_states,number_of_bonds)
        integer :: m, n, atom1, atom2
        double precision :: Rx, Ry, Rz, Rinv

        ds1(:,:) = 0.0d0

        call evaluate_all(position, Natoms, values)

        n = minloc(maxval(values, 1), 1)
        m = maxloc(values(n,:), 1)

        ! Forming bond
        atom1 = forming_bonds(n,m,1)
        atom2 = forming_bonds(n,m,2)
        Rx = position(1,atom1) - position(1,atom2)
        Ry = position(2,atom1) - position(2,atom2)
        Rz = position(3,atom1) - position(3,atom2)
        Rinv = 1.0/sqrt(Rx * Rx + Ry * Ry + Rz * Rz)
        ds1(1,atom1) = ds1(1,atom1) - Rx * Rinv
        ds1(2,atom1) = ds1(2,atom1) - Ry * Rinv
        ds1(3,atom1) = ds1(3,atom1) - Rz * Rinv
        ds1(1,atom2) = ds1(1,atom2) + Rx * Rinv
        ds1(2,atom2) = ds1(2,atom2) + Ry * Rinv
        ds1(3,atom2) = ds1(3,atom2) + Rz * Rinv

        ! Breaking bond
        atom1 = breaking_bonds(n,m,1)
        atom2 = breaking_bonds(n,m,2)
        Rx = position(1,atom1) - position(1,atom2)
        Ry = position(2,atom1) - position(2,atom2)
        Rz = position(3,atom1) - position(3,atom2)
        Rinv = 1.0/sqrt(Rx * Rx + Ry * Ry + Rz * Rz)
        ds1(1,atom1) = ds1(1,atom1) + Rx * Rinv
        ds1(2,atom1) = ds1(2,atom1) + Ry * Rinv
        ds1(3,atom1) = ds1(3,atom1) + Rz * Rinv
        ds1(1,atom2) = ds1(1,atom2) - Rx * Rinv
        ds1(2,atom2) = ds1(2,atom2) - Ry * Rinv
        ds1(3,atom2) = ds1(3,atom2) - Rz * Rinv

    end subroutine gradient

    ! Return the value of the Hessian of the transition state dividing surface
    ! function.
    ! Parameters:
    !   position - A 3 x Natoms array of atomic positions
    !   Natoms - The number of atoms
    ! Returns:
    !   d2s1 - The Hessian of the dividing surface function
    subroutine hessian(position, Natoms, d2s1)

        implicit none
        integer, intent(in) :: Natoms
        double precision, intent(in) :: position(3,Natoms)
        double precision, intent(out) :: d2s1(3,Natoms,3,Natoms)

        double precision :: values(number_of_transition_states,number_of_bonds)
        integer :: m, n, atom1, atom2
        double precision :: Rx, Ry, Rz, Rinv
        double precision :: dxx, dyy, dzz, dxy, dxz, dyz

        d2s1(:,:,:,:) = 0.0d0

        call evaluate_all(position, Natoms, values)

        n = minloc(maxval(values, 1), 1)
        m = maxloc(values(n,:), 1)

        ! Forming bond
        atom1 = forming_bonds(n,m,1)
        atom2 = forming_bonds(n,m,2)
        Rx = position(1,atom1) - position(1,atom2)
        Ry = position(2,atom1) - position(2,atom2)
        Rz = position(3,atom1) - position(3,atom2)
        Rinv = 1.0/sqrt(Rx * Rx + Ry * Ry + Rz * Rz)

        dxx = -(Ry * Ry + Rz * Rz) * (Rinv * Rinv * Rinv)
        dyy = -(Rz * Rz + Rx * Rx) * (Rinv * Rinv * Rinv)
        dzz = -(Rx * Rx + Ry * Ry) * (Rinv * Rinv * Rinv)
        dxy = Rx * Ry * (Rinv * Rinv * Rinv)
        dxz = Rx * Rz * (Rinv * Rinv * Rinv)
        dyz = Ry * Rz * (Rinv * Rinv * Rinv)

        d2s1(1,atom1,1,atom1) = d2s1(1,atom1,1,atom1) + dxx
        d2s1(1,atom1,2,atom1) = d2s1(1,atom1,2,atom1) + dxy
        d2s1(1,atom1,3,atom1) = d2s1(1,atom1,3,atom1) + dxz
        d2s1(2,atom1,1,atom1) = d2s1(2,atom1,1,atom1) + dxy
        d2s1(2,atom1,2,atom1) = d2s1(2,atom1,2,atom1) + dyy
        d2s1(2,atom1,3,atom1) = d2s1(2,atom1,3,atom1) + dyz
        d2s1(3,atom1,1,atom1) = d2s1(3,atom1,1,atom1) + dxz
        d2s1(3,atom1,2,atom1) = d2s1(3,atom1,2,atom1) + dyz
        d2s1(3,atom1,3,atom1) = d2s1(3,atom1,3,atom1) + dzz

        d2s1(1,atom1,1,atom2) = d2s1(1,atom1,1,atom2) - dxx
        d2s1(1,atom1,2,atom2) = d2s1(1,atom1,2,atom2) - dxy
        d2s1(1,atom1,3,atom2) = d2s1(1,atom1,3,atom2) - dxz
        d2s1(2,atom1,1,atom2) = d2s1(2,atom1,1,atom2) - dxy
        d2s1(2,atom1,2,atom2) = d2s1(2,atom1,2,atom2) - dyy
        d2s1(2,atom1,3,atom2) = d2s1(2,atom1,3,atom2) - dyz
        d2s1(3,atom1,1,atom2) = d2s1(3,atom1,1,atom2) - dxz
        d2s1(3,atom1,2,atom2) = d2s1(3,atom1,2,atom2) - dyz
        d2s1(3,atom1,3,atom2) = d2s1(3,atom1,3,atom2) - dzz

        d2s1(1,atom2,1,atom1) = d2s1(1,atom2,1,atom1) - dxx
        d2s1(1,atom2,2,atom1) = d2s1(1,atom2,2,atom1) - dxy
        d2s1(1,atom2,3,atom1) = d2s1(1,atom2,3,atom1) - dxz
        d2s1(2,atom2,1,atom1) = d2s1(2,atom2,1,atom1) - dxy
        d2s1(2,atom2,2,atom1) = d2s1(2,atom2,2,atom1) - dyy
        d2s1(2,atom2,3,atom1) = d2s1(2,atom2,3,atom1) - dyz
        d2s1(3,atom2,1,atom1) = d2s1(3,atom2,1,atom1) - dxz
        d2s1(3,atom2,2,atom1) = d2s1(3,atom2,2,atom1) - dyz
        d2s1(3,atom2,3,atom1) = d2s1(3,atom2,3,atom1) - dzz

        d2s1(1,atom2,1,atom2) = d2s1(1,atom2,1,atom2) + dxx
        d2s1(1,atom2,2,atom2) = d2s1(1,atom2,2,atom2) + dxy
        d2s1(1,atom2,3,atom2) = d2s1(1,atom2,3,atom2) + dxz
        d2s1(2,atom2,1,atom2) = d2s1(2,atom2,1,atom2) + dxy
        d2s1(2,atom2,2,atom2) = d2s1(2,atom2,2,atom2) + dyy
        d2s1(2,atom2,3,atom2) = d2s1(2,atom2,3,atom2) + dyz
        d2s1(3,atom2,1,atom2) = d2s1(3,atom2,1,atom2) + dxz
        d2s1(3,atom2,2,atom2) = d2s1(3,atom2,2,atom2) + dyz
        d2s1(3,atom2,3,atom2) = d2s1(3,atom2,3,atom2) + dzz

        ! Breaking bond
        atom1 = breaking_bonds(n,m,1)
        atom2 = breaking_bonds(n,m,2)
        Rx = position(1,atom1) - position(1,atom2)
        Ry = position(2,atom1) - position(2,atom2)
        Rz = position(3,atom1) - position(3,atom2)
        Rinv = 1.0/sqrt(Rx * Rx + Ry * Ry + Rz * Rz)

        dxx = (Ry * Ry + Rz * Rz) * (Rinv * Rinv * Rinv)
        dyy = (Rz * Rz + Rx * Rx) * (Rinv * Rinv * Rinv)
        dzz = (Rx * Rx + Ry * Ry) * (Rinv * Rinv * Rinv)
        dxy = -Rx * Ry * (Rinv * Rinv * Rinv)
        dxz = -Rx * Rz * (Rinv * Rinv * Rinv)
        dyz = -Ry * Rz * (Rinv * Rinv * Rinv)

        d2s1(1,atom1,1,atom1) = d2s1(1,atom1,1,atom1) + dxx
        d2s1(1,atom1,2,atom1) = d2s1(1,atom1,2,atom1) + dxy
        d2s1(1,atom1,3,atom1) = d2s1(1,atom1,3,atom1) + dxz
        d2s1(2,atom1,1,atom1) = d2s1(2,atom1,1,atom1) + dxy
        d2s1(2,atom1,2,atom1) = d2s1(2,atom1,2,atom1) + dyy
        d2s1(2,atom1,3,atom1) = d2s1(2,atom1,3,atom1) + dyz
        d2s1(3,atom1,1,atom1) = d2s1(3,atom1,1,atom1) + dxz
        d2s1(3,atom1,2,atom1) = d2s1(3,atom1,2,atom1) + dyz
        d2s1(3,atom1,3,atom1) = d2s1(3,atom1,3,atom1) + dzz

        d2s1(1,atom1,1,atom2) = d2s1(1,atom1,1,atom2) - dxx
        d2s1(1,atom1,2,atom2) = d2s1(1,atom1,2,atom2) - dxy
        d2s1(1,atom1,3,atom2) = d2s1(1,atom1,3,atom2) - dxz
        d2s1(2,atom1,1,atom2) = d2s1(2,atom1,1,atom2) - dxy
        d2s1(2,atom1,2,atom2) = d2s1(2,atom1,2,atom2) - dyy
        d2s1(2,atom1,3,atom2) = d2s1(2,atom1,3,atom2) - dyz
        d2s1(3,atom1,1,atom2) = d2s1(3,atom1,1,atom2) - dxz
        d2s1(3,atom1,2,atom2) = d2s1(3,atom1,2,atom2) - dyz
        d2s1(3,atom1,3,atom2) = d2s1(3,atom1,3,atom2) - dzz

        d2s1(1,atom2,1,atom1) = d2s1(1,atom2,1,atom1) - dxx
        d2s1(1,atom2,2,atom1) = d2s1(1,atom2,2,atom1) - dxy
        d2s1(1,atom2,3,atom1) = d2s1(1,atom2,3,atom1) - dxz
        d2s1(2,atom2,1,atom1) = d2s1(2,atom2,1,atom1) - dxy
        d2s1(2,atom2,2,atom1) = d2s1(2,atom2,2,atom1) - dyy
        d2s1(2,atom2,3,atom1) = d2s1(2,atom2,3,atom1) - dyz
        d2s1(3,atom2,1,atom1) = d2s1(3,atom2,1,atom1) - dxz
        d2s1(3,atom2,2,atom1) = d2s1(3,atom2,2,atom1) - dyz
        d2s1(3,atom2,3,atom1) = d2s1(3,atom2,3,atom1) - dzz

        d2s1(1,atom2,1,atom2) = d2s1(1,atom2,1,atom2) + dxx
        d2s1(1,atom2,2,atom2) = d2s1(1,atom2,2,atom2) + dxy
        d2s1(1,atom2,3,atom2) = d2s1(1,atom2,3,atom2) + dxz
        d2s1(2,atom2,1,atom2) = d2s1(2,atom2,1,atom2) + dxy
        d2s1(2,atom2,2,atom2) = d2s1(2,atom2,2,atom2) + dyy
        d2s1(2,atom2,3,atom2) = d2s1(2,atom2,3,atom2) + dyz
        d2s1(3,atom2,1,atom2) = d2s1(3,atom2,1,atom2) + dxz
        d2s1(3,atom2,2,atom2) = d2s1(3,atom2,2,atom2) + dyz
        d2s1(3,atom2,3,atom2) = d2s1(3,atom2,3,atom2) + dzz

    end subroutine hessian

    ! Return the value, gradient, and Hessian of the transition state dividing
    ! surface function.
    ! Parameters:
    !   position - A 3 x Natoms array of atomic positions
    !   Natoms - The number of atoms
    ! Returns:
    !   s1 - The value of the dividing surface function
    !   ds1 - The gradient of the dividing surface function
    !   d2s1 - The Hessian of the dividing surface function
    subroutine evaluate(position, Natoms, s1, ds1, d2s1)

        implicit none
        integer, intent(in) :: Natoms
        double precision, intent(in) :: position(3,Natoms)
        double precision, intent(out) :: s1, ds1(3,Natoms), d2s1(3,Natoms,3,Natoms)

        call value(position, Natoms, s1)
        call gradient(position, Natoms, ds1)
        call hessian(position, Natoms, d2s1)

    end subroutine evaluate

    ! Check that a given position is valid by ensuring that all breaking and
    ! forming bonds have lengths less than Rmax.
    ! Parameters:
    !   position - A 3 x Natoms array of atomic positions
    !   Natoms - The number of atoms
    !   Rmax - The maximum valid bond length in atomic units for breaking/forming bonds
    ! Returns:
    !   result - 0 if the position is valid, 1 if invalid
    subroutine check_for_valid_position(position, Natoms, Rmax, result)

        implicit none
        integer, intent(in) :: Natoms
        double precision, intent(in) :: position(3,Natoms)
        double precision, intent(in) :: Rmax
        integer, intent(inout) :: result

        double precision :: values(number_of_transition_states,number_of_bonds)
        integer :: m, n, atom1, atom2
        double precision :: Rx, Ry, Rz, R

        call evaluate_all(position, Natoms, values)

        n = minloc(maxval(values, 1), 1)

        do m = 1, number_of_bonds
            ! Forming bond
            atom1 = forming_bonds(n,m,1)
            atom2 = forming_bonds(n,m,2)
            Rx = position(1,atom1) - position(1,atom2)
            Ry = position(2,atom1) - position(2,atom2)
            Rz = position(3,atom1) - position(3,atom2)
            R = sqrt(Rx * Rx + Ry * Ry + Rz * Rz)
            if (R .gt. Rmax) then
                result = 1
                return
            end if
            ! Breaking bond
            atom1 = breaking_bonds(n,m,1)
            atom2 = breaking_bonds(n,m,2)
            Rx = position(1,atom1) - position(1,atom2)
            Ry = position(2,atom1) - position(2,atom2)
            Rz = position(3,atom1) - position(3,atom2)
            R = sqrt(Rx * Rx + Ry * Ry + Rz * Rz)
            if (R .gt. Rmax) then
                result = 1
                return
            end if
        end do

    end subroutine

    ! Check that the transition state dividing surface does not give the
    ! same value for equivalent transition states at a given position.
    ! Parameters:
    !   position - A 3 x Natoms array of atomic positions
    !   Natoms - The number of atoms
    ! Returns:
    !   result - 0 if the position is valid, 1 if invalid
    subroutine check_values(position, Natoms, result)

        implicit none
        integer, intent(in) :: Natoms
        double precision, intent(in) :: position(3,Natoms)
        integer, intent(inout) :: result

        double precision :: values(number_of_transition_states,number_of_bonds), max_value
        integer :: m, n, Neq

        call evaluate_all(position, Natoms, values)

        n = minloc(maxval(values, 1), 1)
        m = maxloc(values(n,:), 1)
        max_value = values(n,m)

        Neq = 0
        do n = 1, number_of_transition_states
            if (abs(values(n,m) - max_value) < 1.0e-10) Neq = Neq + 1
        end do

        if (Neq .ne. 1) result = 1

    end subroutine

end module transition_state

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! The ``reactants`` module defines a dividing surface for a set of bimolecular
! reactants, characterized by the distance at which the interaction of the
! reactant molecules becomes negligible. The attributes are:
!
! ========================= ====================================================
! Attribute                 Description
! ========================= ====================================================
! ``Nreactant1_atoms``      The number of atoms in the first reactant
! ``Nreactant2_atoms``      The number of atoms in the second reactant
! ``reactant1_atoms``       An array of indices for each atom in the first reactant
! ``reactant2_atoms``       An array of indices for each atom in the second reactant
! ``massfrac``              The mass fraction of each atom in each reactant
! ``Rinf``                  The distance at which the interaction of the reactants becomes negligible
! ========================= ====================================================
!
module reactants

    implicit none

    integer, parameter :: MAX_ATOMS = 100

    integer :: Nreactant1_atoms
    integer :: Nreactant2_atoms
    integer :: reactant1_atoms(MAX_ATOMS)
    integer :: reactant2_atoms(MAX_ATOMS)
    double precision :: massfrac(MAX_ATOMS)
    double precision :: Rinf

contains

    ! Return the center of mass of the first reactant.
    ! Parameters:
    !   position - A 3 x Natoms array of atomic positions
    !   Natoms - The number of atoms
    ! Returns:
    !   cm - The center of mass of the first reactant
    subroutine reactant1_center_of_mass(position, Natoms, cm)

        implicit none
        integer, intent(in) :: Natoms
        double precision, intent(in) :: position(3,Natoms)
        double precision, intent(out) :: cm(3)

        integer :: n, atom, i

        cm(:) = 0.0
        do n = 1, Nreactant1_atoms
            atom = reactant1_atoms(n)
            do i = 1, 3
                cm(i) = cm(i) + massfrac(atom) * position(i,atom)
            end do
        end do

    end subroutine reactant1_center_of_mass

    ! Return the center of mass of the second reactant.
    ! Parameters:
    !   position - A 3 x Natoms array of atomic positions
    !   Natoms - The number of atoms
    ! Returns:
    !   cm - The center of mass of the second reactant
    subroutine reactant2_center_of_mass(position, Natoms, cm)

        implicit none
        integer, intent(in) :: Natoms
        double precision, intent(in) :: position(3,Natoms)
        double precision, intent(out) :: cm(3)

        integer :: n, atom, i

        cm(:) = 0.0
        do n = 1, Nreactant2_atoms
            atom = reactant2_atoms(n)
            do i = 1, 3
                cm(i) = cm(i) + massfrac(atom) * position(i,atom)
            end do
        end do

    end subroutine reactant2_center_of_mass

    ! Return the value of the bimolecular reactants dividing surface function.
    ! Parameters:
    !   position - A 3 x Natoms array of atomic positions
    !   Natoms - The number of atoms
    ! Returns:
    !   s0 - The value of the reactants dividing surface function
    subroutine value(position, Natoms, s0)

        implicit none
        integer, intent(in) :: Natoms
        double precision, intent(in) :: position(3,Natoms)
        double precision, intent(out) :: s0

        double precision :: cm1(3), cm2(3)
        double precision :: Rx, Ry, Rz, R

        call reactant1_center_of_mass(position, Natoms, cm1)
        call reactant2_center_of_mass(position, Natoms, cm2)

        Rx = cm2(1) - cm1(1)
        Ry = cm2(2) - cm1(2)
        Rz = cm2(3) - cm1(3)
        R = sqrt(Rx * Rx + Ry * Ry + Rz * Rz)

        s0 = Rinf - R

    end subroutine value

    ! Return the gradient of the bimolecular reactants dividing surface
    ! function.
    ! Parameters:
    !   position - A 3 x Natoms array of atomic positions
    !   Natoms - The number of atoms
    ! Returns:
    !   ds0 - The gradient of the reactants dividing surface function
    subroutine gradient(position, Natoms, ds0)

        implicit none
        integer, intent(in) :: Natoms
        double precision, intent(in) :: position(3,Natoms)
        double precision, intent(out) :: ds0(3,Natoms)

        double precision :: cm1(3), cm2(3)
        double precision :: Rx, Ry, Rz, Rinv
        integer :: n, atom

        call reactant1_center_of_mass(position, Natoms, cm1)
        call reactant2_center_of_mass(position, Natoms, cm2)

        Rx = cm2(1) - cm1(1)
        Ry = cm2(2) - cm1(2)
        Rz = cm2(3) - cm1(3)
        Rinv = 1.0/sqrt(Rx * Rx + Ry * Ry + Rz * Rz)

        do n = 1, Nreactant1_atoms
            atom = reactant1_atoms(n)
            ds0(1,atom) = Rx * Rinv * massfrac(atom)
            ds0(2,atom) = Ry * Rinv * massfrac(atom)
            ds0(3,atom) = Rz * Rinv * massfrac(atom)
        end do
        do n = 1, Nreactant2_atoms
            atom = reactant2_atoms(n)
            ds0(1,atom) = -Rx * Rinv * massfrac(atom)
            ds0(2,atom) = -Ry * Rinv * massfrac(atom)
            ds0(3,atom) = -Rz * Rinv * massfrac(atom)
        end do

    end subroutine gradient

    ! Return the Hessian of the bimolecular reactants dividing surface function.
    ! Parameters:
    !   position - A 3 x Natoms array of atomic positions
    !   Natoms - The number of atoms
    ! Returns:
    !   d2s0 - The Hessian of the reactants dividing surface function
    subroutine hessian(position, Natoms, d2s0)

        implicit none
        integer, intent(in) :: Natoms
        double precision, intent(in) :: position(3,Natoms)
        double precision, intent(out) :: d2s0(3,Natoms,3,Natoms)

        double precision :: cm1(3), cm2(3)
        double precision :: Rx, Ry, Rz, Rinv
        double precision :: dxx, dyy, dzz, dxy, dxz, dyz, massfactor
        integer :: n1, atom1, n2, atom2

        call reactant1_center_of_mass(position, Natoms, cm1)
        call reactant2_center_of_mass(position, Natoms, cm2)

        Rx = cm2(1) - cm1(1)
        Ry = cm2(2) - cm1(2)
        Rz = cm2(3) - cm1(3)
        Rinv = 1.0/sqrt(Rx * Rx + Ry * Ry + Rz * Rz)

        dxx = -(Ry * Ry + Rz * Rz) * (Rinv * Rinv * Rinv)
        dyy = -(Rz * Rz + Rx * Rx) * (Rinv * Rinv * Rinv)
        dzz = -(Rx * Rx + Ry * Ry) * (Rinv * Rinv * Rinv)
        dxy = Rx * Ry * (Rinv * Rinv * Rinv)
        dxz = Rx * Rz * (Rinv * Rinv * Rinv)
        dyz = Ry * Rz * (Rinv * Rinv * Rinv)

        do n1 = 1, Nreactant1_atoms
            atom1 = reactant1_atoms(n1)
            do n2 = 1, Nreactant1_atoms
                atom2 = reactant1_atoms(n2)
                massfactor = massfrac(atom1) * massfrac(atom2)
                d2s0(1,atom1,1,atom2) = dxx * massfactor
                d2s0(1,atom1,2,atom2) = dxy * massfactor
                d2s0(1,atom1,3,atom2) = dxz * massfactor
                d2s0(2,atom1,1,atom2) = dxy * massfactor
                d2s0(2,atom1,2,atom2) = dyy * massfactor
                d2s0(2,atom1,3,atom2) = dyz * massfactor
                d2s0(3,atom1,1,atom2) = dxz * massfactor
                d2s0(3,atom1,2,atom2) = dyz * massfactor
                d2s0(3,atom1,3,atom2) = dzz * massfactor
            end do
            do n2 = 1, Nreactant2_atoms
                atom2 = reactant2_atoms(n2)
                massfactor = massfrac(atom1) * massfrac(atom2)
                d2s0(1,atom1,1,atom2) = -dxx * massfactor
                d2s0(1,atom1,2,atom2) = -dxy * massfactor
                d2s0(1,atom1,3,atom2) = -dxz * massfactor
                d2s0(2,atom1,1,atom2) = -dxy * massfactor
                d2s0(2,atom1,2,atom2) = -dyy * massfactor
                d2s0(2,atom1,3,atom2) = -dyz * massfactor
                d2s0(3,atom1,1,atom2) = -dxz * massfactor
                d2s0(3,atom1,2,atom2) = -dyz * massfactor
                d2s0(3,atom1,3,atom2) = -dzz * massfactor
            end do
        end do
        do n1 = 1, Nreactant2_atoms
            atom1 = reactant2_atoms(n1)
            do n2 = 1, Nreactant1_atoms
                atom2 = reactant1_atoms(n2)
                massfactor = massfrac(atom1) * massfrac(atom2)
                d2s0(1,atom1,1,atom2) = -dxx * massfactor
                d2s0(1,atom1,2,atom2) = -dxy * massfactor
                d2s0(1,atom1,3,atom2) = -dxz * massfactor
                d2s0(2,atom1,1,atom2) = -dxy * massfactor
                d2s0(2,atom1,2,atom2) = -dyy * massfactor
                d2s0(2,atom1,3,atom2) = -dyz * massfactor
                d2s0(3,atom1,1,atom2) = -dxz * massfactor
                d2s0(3,atom1,2,atom2) = -dyz * massfactor
                d2s0(3,atom1,3,atom2) = -dzz * massfactor
            end do
            do n2 = 1, Nreactant2_atoms
                atom2 = reactant2_atoms(n2)
                massfactor = massfrac(atom1) * massfrac(atom2)
                d2s0(1,atom1,1,atom2) = dxx * massfactor
                d2s0(1,atom1,2,atom2) = dxy * massfactor
                d2s0(1,atom1,3,atom2) = dxz * massfactor
                d2s0(2,atom1,1,atom2) = dxy * massfactor
                d2s0(2,atom1,2,atom2) = dyy * massfactor
                d2s0(2,atom1,3,atom2) = dyz * massfactor
                d2s0(3,atom1,1,atom2) = dxz * massfactor
                d2s0(3,atom1,2,atom2) = dyz * massfactor
                d2s0(3,atom1,3,atom2) = dzz * massfactor
            end do
        end do

    end subroutine hessian

    ! Return the value, gradient, and Hessian of the bimolecular reactants
    ! dividing surface function.
    ! Parameters:
    !   position - A 3 x Natoms array of atomic positions
    !   Natoms - The number of atoms
    ! Returns:
    !   s0 - The value of the dividing surface function
    !   ds0 - The gradient of the dividing surface function
    !   d2s0 - The Hessian of the dividing surface function
    subroutine evaluate(position, Natoms, s1, ds1, d2s1)

        implicit none
        integer, intent(in) :: Natoms
        double precision, intent(in) :: position(3,Natoms)
        double precision, intent(out) :: s1, ds1(3,Natoms), d2s1(3,Natoms,3,Natoms)

        call value(position, Natoms, s1)
        call gradient(position, Natoms, ds1)
        call hessian(position, Natoms, d2s1)

    end subroutine evaluate

end module reactants
