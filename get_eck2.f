      subroutine get_eck(geom,eck)
      implicit real*8 (a-h,o-z)
      include 'param.inc'
      parameter (natoms = natom)
      common/masses/um(natoms),a0(3,natoms),utot
      dimension f(2,2),finv(2,2),transf(3,3),eig(2),tmp(2),f0(2,2)
      dimension eck(3,3),vec(3,3),geom(3,natoms),com(3)
      dimension x(natoms,3)
        do j = 1,3
            com(j) = 0.d0
        enddo
      do i = 1,natoms
         do j = 1,3
            com(j) = com(j) + um(i)*geom(j,i)
         enddo
      enddo
      vec = 0.
      do j = 1,3
         com(j) = com(j)/utot
         do i = 1,natoms
            g0 = geom(j,i)
            geom(j,i) = geom(j,i)-com(j)
            do k = 1,2
               vec(j,k) = vec(j,k) + um(i)*geom(j,i)*a0(k,i)
            enddo
         enddo
      enddo
      f = 0.
      do i = 1,2
         do j = 1,2
            do k = 1,3
               f(j,i) = f(j,i) + vec(k,j)*vec(k,i)
            enddo
         enddo
      enddo
      f0 = f
      call house(f,2,2,eig,tmp)
      do j = 1,2
         eig(j) = 1./sqrt(eig(j))
      enddo
      finv = 0.
      do i = 1,2
         do j = 1,2
            do k = 1,2
               finv(j,i) = finv(j,i)+eig(k)*f(j,k)*f(i,k)
            enddo
         enddo
      enddo
      eck = 0.
      do i = 1,2
         do j = 1,3
            do k = 1,2
               eck(j,i) = eck(j,i) + vec(j,k)*finv(k,i)
            enddo
         enddo
      enddo
      eck(1,3) = eck(2,1)*eck(3,2)-eck(3,1)*eck(2,2)
      eck(2,3) = eck(3,1)*eck(1,2)-eck(1,1)*eck(3,2)
      eck(3,3) = eck(1,1)*eck(2,2)-eck(2,1)*eck(1,2)
      x = 0.
      d = 0.
      do i = 1,natoms
         do j = 1,3
            do k = 1,3
               x(i,j) = x(i,j) + geom(k,i)*eck(k,j)
            enddo
         enddo
      enddo
      do i = 1,natoms
         do j = 1,3
	    geom(j,i) = x(i,j)
         enddo
      enddo
      return
      end
