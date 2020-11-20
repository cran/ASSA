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

subroutine yyi_(YA, YB, U, l, k, ansA, ansB )
  implicit none
  integer :: i, j, k, l, ll
  double precision, dimension(l, l) :: U
  double precision, dimension(l, k) :: YA, YB, ansA, ansB
 
do j = 1, k
 do i = 1, l
  do ll = 1, l
      ansA(i,j) = ansA(i,j) + U(i,ll)*YA(ll,j)
      ansB(i,j) = ansB(i,j) + U(i,ll)*YB(ll,j)
     end do
   end do
 end do
end subroutine
