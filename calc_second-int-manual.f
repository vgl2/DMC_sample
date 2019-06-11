        subroutine calc_second_psi(coord,sig_bend,nmono,second_mat)
        implicit real*8(a-h,o-z)
        parameter (natoms=6)
        parameter (nvib=6)
        parameter (nvibw=3)
        parameter (nparam=500)
        common/wa/roh1(nparam),wave(nparam),feed1(nparam)
        common/pcoef/alpha,c2,re,rea,c1
        dimension coord(natoms*3),test(5,5),r(nvibw,nmono),
     1  second(nvibw,nvibw,nmono),second_mat(nvib,nvib)
        pi = dacos(-1.d0)
        call calc_r(1,coord(:),nmono,r(:,:))
        dr = 0.0001
        do i = 1,5
            do j = 1,5
                test(i,j) = 0.d0
            enddo
        enddo
        do m = 1,nmono
            do k = 1,nvibw
                do l = k,nvibw
                    do i = 1,5
                        if (i.eq.1) then
                            r(k,m) = r(k,m) - 2*dr
                        else if (i.eq.2) then
                            r(k,m) = r(k,m) - dr
                        else if (i.eq.3) then
                            r(k,m) = r(k,m) 
                        else if (i.eq.4) then  
                            r(k,m) = r(k,m) + dr
                        else if (i.eq.5) then  
                            r(k,m) = r(k,m) + 2*dr
                        endif
                        do j = 1,5
                            if (j.eq.1) then
                                r(l,m) = r(l,m) + 2*dr
                                call calc_psi(1,r(:,:),sig_bend,nmono,
     1                          psi)
                                test(i,1) = psi
                                r(l,m) = r(l,m) - 2*dr
                            else if (j.eq.2) then
                                r(l,m) = r(l,m) + dr
                                call calc_psi(1,r(:,:),sig_bend,nmono,
     1                          psi)
                                test(i,2) = psi
                                r(l,m) = r(l,m) - dr
                            else if (j.eq.3) then
                                r(l,m) = r(l,m) 
                                call calc_psi(1,r(:,:),sig_bend,nmono,
     1                          psi)
                                test(i,3) = psi
                                r(l,m) = r(l,m) 
                            else if (j.eq.4) then  
                                r(l,m) = r(l,m) - dr
                                call calc_psi(1,r(:,:),sig_bend,nmono,
     1                          psi)
                                test(i,4) = psi
                                r(l,m) = r(l,m) + dr
                            else if (j.eq.5) then  
                                r(l,m) = r(l,m) - 2*dr
                                call calc_psi(1,r(:,:),sig_bend,nmono,
     1                          psi)
                                test(i,5) = psi
                                r(l,m) = r(l,m) + 2*dr
                            endif
                        enddo
                        if (i.eq.1) then
                            r(k,m) = r(k,m) + 2*dr
                        else if (i.eq.2) then
                            r(k,m) = r(k,m) + dr
                        else if (i.eq.3) then
                            r(k,m) = r(k,m) 
                        else if (i.eq.4) then  
                            r(k,m) = r(k,m) - dr
                        else if (i.eq.5) then  
                            r(k,m) = r(k,m) - 2*dr
                        endif
                    enddo
                    if (k.eq.l) then
                        second(k,l,m)=(-test(3,5)+(16*test(3,4))
     1                  -(30*(test(3,3)))+(16*(test(3,2)))-test(3,1))/
     1                  (12*dr**2)
                    else 
                        second(k,l,m)=((8*(test(4,5)+test(5,4)+test(1,2)
     1                  +test(2,1)))-(8*(test(2,5)+test(1,4)+test(4,1)
     1                  +test(5,2)))-(test(5,5)+test(1,1)-test(5,1)-
     1                  test(1,5))+(64*(test(2,4)+test(4,2)-test(4,4)
     1                  -test(2,2))))/(144.d0*dr**2)
                    endif
                    do i = 1,5
                        do j = 1,5
                            test(i,j) = 0.d0
                        enddo
                    enddo
                enddo
            enddo
            do k = 1,nvibw
                do l = 1,nvibw
                    second(l,k,m) = second(k,l,m)
                enddo
            enddo
        enddo
        do j =1 ,nmono
            if (j.eq.1) then
                do k = 1,nvibw
                    do l = 1,nvibw
                        second_mat(k,l) = second(k,l,j)
                    enddo
                enddo
            endif
            if (j.eq.2) then
                do k = 1,nvibw
                    do l = 1,nvibw
                        second_mat(k+nvibw,l+nvibw) = second(k,l,j)
                    enddo
                enddo
            endif
            if (j.eq.3) then
                do k = 1,nvibw
                    do l = 1,nvibw
                        second_mat(k+2*nvibw,l+2*nvibw) = second(k,l,j)
                    enddo
                enddo
            endif
            if (j.eq.4) then
                do k = 1,nvibw
                    do l = 1,nvibw
                        second_mat(k+3*nvibw,l+3*nvibw) = second(k,l,j)
                    enddo
                enddo
            endif
            if (j.eq.5) then
                do k = 1,nvibw
                    do l = 1,nvibw
                        second_mat(k+4*nvibw,l+4*nvibw) = second(k,l,j)
                    enddo
                enddo
            endif
            if (j.eq.6) then
                do k = 1,nvibw
                    do l = 1,nvibw
                        second_mat(k+5*nvibw,l+5*nvibw) = second(k,l,j)
                    enddo
                enddo
            endif
        enddo
        return
        end subroutine


        subroutine calc_first_psi(coord,iint,sig_bend,first)
        implicit real*8(a-h,o-z)
        parameter (natoms=3)
        parameter (nvib=3)
        parameter (nparam=500)
        common/wa/roh1(nparam),wave(nparam),feed1(nparam)
        common/pcoef/alpha,c2,re,rea,c1
        dimension coord(natoms*3),test(5),r(nvib),wf(nvib)
        pi = dacos(-1.d0)
        r(1) = sqrt((coord(4)-coord(1))**2 + (coord(5)-coord(2))**2 +
     1  (coord(6)-coord(3))**2)
        r(2) = sqrt((coord(7)-coord(1))**2 + (coord(8)-coord(2))**2 +
     1  (coord(9)-coord(3))**2)
        r(3) = dacos((((coord(4)-coord(1))*(coord(7)-coord(1)))+
     1  ((coord(5)-coord(2))*(coord(8)-coord(2)))+((coord(6)-coord(3))*
     1  (coord(9)-coord(3))))/(r(1)*r(2)))
        dr = 0.0001
        do j = 1,5
            test(j) = 0.d0
        enddo
        do j = 1,5
            if (j.eq.1) then
                r(iint) = r(iint) + 2*dr
                call calc_psi(1,r(:),sig_bend,wf(:),psi)
                test(1) = psi
                r(iint) = r(iint) - 2*dr
            else if (j.eq.2) then
                r(iint) = r(iint) + dr
                call calc_psi(1,r(:),sig_bend,wf(:),psi)
                test(2) = psi
                r(iint) = r(iint) - dr
            else if (j.eq.3) then
                r(iint) = r(iint)
                call calc_psi(1,r(:),sig_bend,wf(:),psi)
                test(3) = psi
                r(iint) = r(iint)
            else if (j.eq.4) then
                r(iint) = r(iint) - dr
                call calc_psi(1,r(:),sig_bend,wf(:),psi)
                test(4) = psi
                r(iint) = r(iint) + dr
            else if (j.eq.5) then
                r(iint) = r(iint) - 2*dr
                call calc_psi(1,r(:),sig_bend,wf(:),psi)
                test(5) = psi
                r(iint) = r(iint) + 2*dr
            endif
        enddo
        first = (-(test(1))+(8*(test(2)))-(8*(test(4)))+(test(5)))/
     1  (12*dr)
        return
        end subroutine

