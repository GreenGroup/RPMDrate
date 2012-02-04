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

! Seed the pseudo-random number generator using the system clock. The
! implementation is based on an example in the gfortran documentation for
! the intrinsic random_seed() subroutine. Note that multiple calls to this
! subroutine do not reseed the pseudo-random number generator.
subroutine random_init()
    implicit none
    integer :: i, n, clock
    integer, dimension(:), allocatable :: seed

    call random_seed(size=n)
    allocate(seed(n))

    call system_clock(count=clock)

    seed = clock + 37 * (/ (i - 1, i = 1, n) /)
    call random_seed(put = seed)

    deallocate(seed)

end subroutine random_init

! Compute a pseudo-random number uniformly distributed in [0,1].
! Returns:
!   rn - The pseudo-random number
subroutine random(rn)
    implicit none
    double precision, intent(out) :: rn
    call random_number(rn)
end subroutine random

! Compute a pseudo-random number weighted by a standard normal distribution.
! The Marsaglia polar method is used to convert from a uniform distribution to
! a normal distribution.
! Returns:
!   rn - The pseudo-random number
subroutine randomn(rn)

    implicit none
    double precision, intent(out) :: rn
    integer :: iset
    double precision :: gset, u, v, S, fac

    save iset,gset

    data iset/0/

    ! The Marsaglia polar method
    if (iset .eq. 0) then
        S = 1.d0
        do while (S .ge. 1.d0 .or. S .eq. 0.d0)
            call random_number(u)
            call random_number(v)
            u = 2.d0 * u - 1.d0
            v = 2.d0 * v - 1.d0
            S = u * u + v * v
        end do
        fac = sqrt(-2 * log(S) / S)
        gset = u * fac
        rn = v * fac
        iset = 1
    else
        rn = gset
        iset = 0
    end if

end subroutine randomn
