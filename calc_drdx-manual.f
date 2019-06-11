        subroutine calc_drdx(coord,icart,iint,first)
        implicit real*8(a-h,o-z)
        parameter (natoms=3)
        parameter (nvib=3)
        parameter (nparam=500)
        common/wa/roh1(nparam),wave(nparam),feed1(nparam)
        common/pcoef/alpha,c2,re,rea,c1
        dimension coord(natoms*3),r(nvib),test(5)
        pi = dacos(-1.d0)
        dr = 0.0001
        do i = 1,5
            test(i) = 0.d0
        enddo
        do j = 1,5
            if (j.eq.1) then
                coord(icart) = coord(icart) - 2*dr
                call calc_r(1,coord(:),r(:))
                if (iint.eq.1) then
                    test(1) = r(1)
                elseif (iint.eq.2) then
                    test(1) = r(2)
                elseif (iint.eq.3) then
                    test(1) = r(3)
                endif
                coord(icart) = coord(icart) + 2*dr
            else if (j.eq.2) then
                coord(icart) = coord(icart) - dr
                call calc_r(1,coord(:),r(:))
                if (iint.eq.1) then
                    test(2) = r(1)
                elseif (iint.eq.2) then
                    test(2) = r(2)
                elseif (iint.eq.3) then
                    test(2) = r(3)
                endif
                coord(icart) = coord(icart) + dr
            else if (j.eq.3) then
                coord(icart) = coord(icart) 
                call calc_r(1,coord(:),r(:))
                if (iint.eq.1) then
                    test(3) = r(1)
                elseif (iint.eq.1) then
                    test(3) = r(2)
                elseif (iint.eq.3) then
                    test(3) = r(3)
                endif
                coord(icart) = coord(icart) 
            else if (j.eq.4) then  
                coord(icart) = coord(icart) + dr
                call calc_r(1,coord(:),r(:))
                if (iint.eq.1) then
                    test(4) = r(1)
                elseif (iint.eq.2) then
                    test(4) = r(2)
                elseif (iint.eq.3) then
                    test(4) = r(3)
                endif
                coord(icart) = coord(icart) - dr
            else if (j.eq.5) then  
                coord(icart) = coord(icart) + 2*dr
                call calc_r(1,coord(:),r(:))
                if (iint.eq.1) then
                    test(5) = r(1)
                elseif (iint.eq.2) then
                    test(5) = r(2)
                elseif (iint.eq.3) then
                    test(5) = r(3)
                endif
                coord(icart) = coord(icart) - 2*dr
            endif
        enddo
        first=(test(1)-(8*test(2))+(8*test(4))-(test(5)))/(12*dr)
        return
        end subroutine



        subroutine calc_drdx_sq(coord,icart,iint,second)
        implicit real*8(a-h,o-z)
        parameter (natoms=3)
        parameter (nvib=3)
        parameter (nparam=500)
        common/wa/roh1(nparam),wave(nparam),feed1(nparam)
        common/pcoef/alpha,c2,re,rea,c1
        dimension x(natoms,3),coord(natoms*3),r(nvib),test(5)
        pi = dacos(-1.d0)
        dr = 0.0001
        do i = 1,5
            test(i) = 0.d0
        enddo
        do j = 1,5
            if (j.eq.1) then
                coord(icart) = coord(icart) - 2*dr
                call calc_r(1,coord(:),r(:))
                if (iint.eq.1) then
                    test(1) = r(1)
                elseif (iint.eq.2) then
                    test(1) = r(2)
                elseif (iint.eq.3) then
                    test(1) = r(3)
                endif
                coord(icart) = coord(icart) + 2*dr
            else if (j.eq.2) then
                coord(icart) = coord(icart) - dr
                call calc_r(1,coord(:),r(:))
                if (iint.eq.1) then
                    test(2) = r(1)
                elseif (iint.eq.2) then
                    test(2) = r(2)
                elseif (iint.eq.3) then
                    test(2) = r(3)
                endif
                coord(icart) = coord(icart) + dr
            else if (j.eq.3) then
                coord(icart) = coord(icart) 
                call calc_r(1,coord(:),r(:))
                if (iint.eq.1) then
                    test(3) = r(1)
                elseif (iint.eq.2) then
                    test(3) = r(2)
                elseif (iint.eq.3) then
                    test(3) = r(3)
                endif
                coord(icart) = coord(icart) 
            else if (j.eq.4) then  
                coord(icart) = coord(icart) + dr
                call calc_r(1,coord(:),r(:))
                if (iint.eq.1) then
                    test(4) = r(1)
                elseif (iint.eq.2) then
                    test(4) = r(2)
                elseif (iint.eq.3) then
                    test(4) = r(3)
                endif
                coord(icart) = coord(icart) - dr
            else if (j.eq.5) then  
                coord(icart) = coord(icart) + 2*dr
                call calc_r(1,coord(:),r(:))
                if (iint.eq.1) then
                    test(5) = r(1)
                elseif (iint.eq.2) then
                    test(5) = r(2)
                elseif (iint.eq.3) then
                    test(5) = r(3)
                endif
                coord(icart) = coord(icart) - 2*dr
            endif
        enddo
        second=(-(test(1))+(16*test(2))-(30*test(3))+(16*test(4))-
     1  (test(5)))/(12*(dr)**2)
        return
        end subroutine
