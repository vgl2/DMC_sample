        subroutine get_x(r,psips)
        implicit real*8(a-h,o-z)
        parameter (nvib=3)
        parameter (ndim=9)
        dimension psips(ndim),r(nvib)
        psips(1) = 0.d0
        psips(2) = 0.d0
        psips(3) = 0.d0
        psips(4) = r(1)
        psips(5) = 0.d0
        psips(6) = 0.d0
        psips(7) = r(2)*dcos(r(3))
        psips(8) = r(2)*dsin(r(3))
        psips(9) = 0.d0
        return
        end subroutine
