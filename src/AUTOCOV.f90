 !  ==========================================================================  !
 !  G Martos Venturini                                                          !
 !  Copyright (C) 2020                                                          !
 !  --------------------------------------------------------------------------  !
 !  This program is free software; you can redistribute it and/or modify        !
 !  it under the terms of the GNU General Public License as published by        !
 !  the Free Software Foundation; either version 2 of the License, or           !
 !  (at your option) any later version.                                         !
 !                                                                              !
 !  This program is distributed in the hope that it will be useful,             !
 !  but WITHOUT ANY WARRANTY; without even the implied warranty of              !
 !  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the               !
 !  GNU General Public License for more details.                                !
 !  ==========================================================================  !

subroutine autocov_(e, n, answer)
  implicit none
  integer :: i, j, n
  double precision, dimension(n) :: answer
  double precision, dimension(n, 2) :: e
  do i = 1, n
    answer(i) = 0
  end do
    do i = 1, n
      do j = 1, (n-i+1)
           answer(i) = answer(i) + 2*e(j+i-1,1)*e(j,1) + e(j+i-1,1)*e(j,2) + e(j+i-1,1)*e(j,2) + 2*e(j+i-1,2)*e(j,2) 
      end do
      answer(i) = answer(i)/(6*(n))
    end do
end subroutine
