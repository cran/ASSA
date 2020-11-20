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

subroutine yi_(YA, YB, u, l, k, ansA, ansB )
  implicit none
  
  integer :: i, j, k, l
  double precision, dimension(k) :: vA, vB
  double precision, dimension(l) :: u
  double precision, dimension(k, l) :: YA, YB
  double precision, dimension(l, k) :: ansA, ansB
  do j = 1, k
    vA(j) = 0
    vB(j) = 0
          do i = 1, l
          ansA(i,j) = 0
          ansB(i,j) = 0
      end do    
  end do
  
    do j = 1, k
      do i = 1, l
        if(u(i) > 0 ) then
          vA(j) = vA(j) + YA(j,i)*u(i)
          vB(j) = vB(j) + YB(j,i)*u(i)
        end if
        if(u(i) < 0 ) then
          vA(j) = vA(j) + YB(j,i)*u(i)
          vB(j) = vB(j) + YA(j,i)*u(i)
        end if
      end do
 end do
     do j = 1, k
      do i = 1, l
        if(u(i) > 0 ) then
          ansA(i,j) = ansA(i,j) + vA(j)*u(i)
          ansB(i,j) = ansB(i,j) + vB(j)*u(i)
        end if
        if(u(i) < 0 ) then
          ansA(i,j) = ansA(i,j) + vB(j)*u(i)
          ansB(i,j) = ansB(i,j) + vA(j)*u(i)
        end if
      end do
 end do
end subroutine
