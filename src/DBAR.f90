 !  ==========================================================================  !
 !  Miguel de Carvalho                                                          !
 !  Copyright (C) 2017                                                          !
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

subroutine dbar_(y, l, k, answer)
  implicit none
  
  integer :: i, j, k, l, ll, n
  double precision, dimension(k + l - 1) :: answer, count, sum
  double precision, dimension(l, k) :: Y
  
  n = l + (k - 1)

  do i = 1, n
    count(i) = 0
    sum(i) = 0
  end do
  
  do ll = 2, (k + l)
    do i = 1, l
      do j = 1, k
        if(i + j == ll) then
          count(ll - 1) = count(ll - 1) + 1
          sum(ll - 1) = sum(ll - 1) + Y(i, j)
          answer(ll - 1) = sum(ll - 1) / count(ll - 1)
        end if
      end do
    end do
 end do
end subroutine
