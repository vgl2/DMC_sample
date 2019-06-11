        subroutine calc_r(n,coord,nmono,r)
        implicit real*8(a-h,o-z)
        parameter (nmax=25000)
        parameter (ndim=27)
        parameter (ndimw=9)
        parameter (nvibw=3)
        parameter (nvib=9)
        dimension coord(ndim,nmax),coord_mono(ndimw,nmono,nmax),
     1  r(nvibw,nmono,nmax)
        do i = 1,n
            ip = 0
            do j = 1,nmono
                do k = 1,ndimw
                    ip = ip + 1
                    coord_mono(k,j,i) = coord(ip,i)
                enddo
            enddo
        enddo
        do i = 1,n
            do j = 1,nmono
                r(1,j,i)=sqrt(((coord_mono(4,j,i)-coord_mono(1,j,i))**2)
     1          +((coord_mono(5,j,i)-coord_mono(2,j,i))**2)+
     1          ((coord_mono(6,j,i)-coord_mono(3,j,i))**2))
                r(2,j,i)=sqrt(((coord_mono(7,j,i)-coord_mono(1,j,i))**2)
     1          +((coord_mono(8,j,i)-coord_mono(2,j,i))**2)+
     1          ((coord_mono(9,j,i)-coord_mono(3,j,i))**2))
                r(3,j,i)=dacos((((coord_mono(4,j,i)-coord_mono(1,j,i))*
     1          (coord_mono(7,j,i)-coord_mono(1,j,i)))+
     1          ((coord_mono(5,j,i)-coord_mono(2,j,i))*
     1          (coord_mono(8,j,i)-coord_mono(2,j,i)))+
     1          ((coord_mono(6,j,i)-coord_mono(3,j,i))*
     1          (coord_mono(9,j,i)-coord_mono(3,j,i))))/(r(1,j,i)*
     1          r(2,j,i)))
            enddo
        enddo
        return
        end subroutine

        subroutine calc_psi(n,r,sig_bend,nmono,psi)
        implicit real*8(a-h,o-z)
        parameter (nmax=25000)
        parameter (nvibw=3)
        parameter (nparam=500)
        common/pcoef/alpha,c2,re,rea,c1
        common/wa1/roh1(nparam),wave_oh(nparam),feed_oh(nparam)
        common/wa2/rod1(nparam),wave_od(nparam),feed_od(nparam)
        dimension r(nvibw,nmono,nmax),wf(nvibw,nmono,nmax),psi(nmax)
        pi = dacos(-1.d0)
        do i = 1,n
            psi(i) = 1.d0
        enddo
        do i = 1,n
            do j = 1,nmono
                do k = 1,nvibw
                    if (k.ne.3) then
                        call splint(roh1(:),wave_oh(:),feed_oh(:),nparam
     1                  ,r(k,j,i),wf(k,j,i))
                        psi(i) = psi(i)*wf(k,j,i)
c                       psi(i) = psi(i)
                    else
                        wf(k,j,i) = (alpha/((sig_bend**2)*pi))
     1                  **(0.25)*exp(-0.5*alpha*(((r(k,j,i)-
     1                  rea)**2)/(sig_bend**2)))
                        psi(i) = psi(i)*wf(k,j,i)
                    endif
                enddo
            enddo
        enddo
        return
        end subroutine

        subroutine calc_grad(n,r,sig_bend,psi,nmono,dr,first_mat,
     1  second_mat)
        implicit real*8(a-h,o-z)
        parameter (nmax=25000)
        parameter (nvibw=3)
        parameter (nvib=9)
        parameter (nparam=500)
        common/pcoef/alpha,c2,re,rea,c1
        common/wa1/roh1(nparam),wave_oh(nparam),feed_oh(nparam)
        common/wa2/rod1(nparam),wave_od(nparam),feed_od(nparam)
        dimension r(nvibw,nmono,nmax),dr(nvibw,nmono,nmax),psi(nmax),
     1  first(nvibw,nmono,nmax),second(nvibw,nmono,nmax),
     1  first_mat(nvib,nmax),second_mat(nvib,nmax)
        do i = 1,n
            do j = 1,nmono
                do k = 1,nvibw
c       calculate first derivative to get the drift term
                    if (k.ne.3) then
                        call calc_first(r(:,:,i),first(k,j,i),
     1                  second(k,j,i),k,nmono,j,sig_bend)
                        dr(k,j,i) = (2*first(k,j,i))/psi(i)
c                    dr(j,i) = 0.d0
                    else
                        dr(k,j,i)=(-2*alpha*(r(k,j,i)-rea)/
     1                  (sig_bend**2))
                        first(k,j,i) = dr(k,j,i)*psi(i)/2.d0
                        second(k,j,i) = ((alpha**2*((r(k,j,i)-
     1                  rea)**2))-alpha)*psi(i)
                    endif
                enddo
            enddo
        enddo
        do i = 1,n
            ip = 0
            do j = 1,nmono
                do k = 1,nvibw
                    ip = ip + 1
                    first_mat(ip,i) = first(k,j,i)
                    second_mat(ip,i) = second(k,j,i)
                enddo
            enddo
        enddo
        return 
        end subroutine

        subroutine calc_t(n,coord,psi,sig_bend,nmono,first_mat,
     1  second_mat,t)
        implicit real*8(a-h,o-z)
        parameter (nmax=25000)
        parameter (nvib=9)
        parameter (nparam=500)
        parameter (ndim=27)
        common/pcoef/alpha,c2,re,rea,c1
        common/kin/weight(ndim),d(ndim)
        common /wa1/roh1(nparam),wave_oh(nparam),feed_oh(nparam)
        common /wa2/rod1(nparam),wave_od(nparam),feed_od(nparam)
        common/der/drdx(nvib,ndim,nmax),dr2dx2(nvib,ndim,nmax)
        dimension coord(ndim,nmax),psi(nmax),t(nmax),second(nvib,nvib),
     1  tot_deriv(ndim,nmax),tot_deriv1(ndim,nmax),tot_deriv2(ndim,nmax)
     1  ,first_mat(nvib,nmax),second_mat(nvib,nmax)
        do i = 1,n
            t(i) = 0.d0
            do k = 1,ndim
                tot_deriv1(k,i) = 0.d0
                tot_deriv(k,i) = 0.d0
                tot_deriv2(k,i) = 0.d0
            enddo
        enddo
        do i = 1,n
            do k =1,ndim
                do l = 1,nvib
                    do j = 1,nvib
                        if (l.eq.j) then
                            tot_deriv1(k,i) = tot_deriv1(k,i)+
     1                      (second_mat(l,i)*drdx(l,k,i)*drdx(j,k,i))
                        else
                            tot_deriv1(k,i) = tot_deriv1(k,i)+
     1                      (((first_mat(l,i)*drdx(l,k,i)*first_mat(j,i)
     1                      *drdx(j,k,i)))/(psi(i)))
                        endif
                    enddo
                enddo
                do l = 1,nvib
                    tot_deriv2(k,i) = tot_deriv2(k,i)+(first_mat(l,i)*
     1              dr2dx2(l,k,i))
                enddo
                tot_deriv(k,i) = tot_deriv1(k,i) + tot_deriv2(k,i)
            enddo
        enddo
        do i = 1,n
            do j =1,ndim
                t(i) = t(i) -(0.5*(tot_deriv(j,i)/weight(j)))
            enddo
        enddo
        do i = 1,n
            t(i) = t(i)/psi(i)
        enddo
        return
        end
            
C
c       calculate second derivatives        
        
        
        
        subroutine initial_roh_wave(sig)
        implicit real*8(a-h,o-z)
        parameter (nparam=500)
        common/wa1/roh1(nparam),wave_oh(nparam),feed_oh(nparam)
c        common/pcoef/alpha,c2,re,rea,c1
        dimension harm(nparam),wave2(nparam)
        open(unit=200,file='roh_morse.dat',status='old')
        do i = 1,nparam
            read(200,*) roh1(i),harm(i),wave2(i),wave_oh(i)
        enddo
        maxwave = maxloc(wave_oh,1)
        re = roh1(maxwave)
        do i =1 ,nparam
           roh1(i) = roh1(i) - re
        enddo 
        do i = 1,nparam
            roh1(i) = roh1(i)*sig
        enddo
        do i = 1,nparam
            roh1(i) = roh1(i) + re
        enddo
        dx = roh1(2)-roh1(1)
        tot_wave = 0.d0
        do i = 1,nparam
            tot_wave = tot_wave + wave_oh(i)**2*dx
        enddo
        do i = 1,nparam
            wave_oh(i) = wave_oh(i)/sqrt(tot_wave)
        enddo
        pt1 = (wave_oh(2)-wave_oh(1))/(roh1(2)-roh1(1))
        ptn =(wave_oh(nparam)-wave_oh(nparam-1))/(roh1(nparam)-
     1  roh1(nparam-1))
        call spline(roh1,wave_oh,nparam,pt1,ptn,feed_oh)
        return
        end         

        subroutine initial_rod_wave(sig)
        implicit real*8(a-h,o-z)
        parameter (nparam=500)
        common/wa2/rod1(nparam),wave_od(nparam),feed_od(nparam)
c        common/pcoef/alpha,c2,re,rea,c1
        dimension harm(nparam),wave2(nparam)
        open(unit=200,file='rod_morse.dat',status='old')
        do i = 1,nparam
            read(200,*) rod1(i),harm(i),wave2(i),wave_od(i)
        enddo
        maxwave = maxloc(wave_od,1)
        re = rod1(maxwave)
        do i =1 ,nparam
           rod1(i) = rod1(i) - re
        enddo 
        do i = 1,nparam
            rod1(i) = rod1(i)*sig
        enddo
        do i = 1,nparam
            rod1(i) = rod1(i) + re
        enddo
        dx = rod1(2)-rod1(1)
        tot_wave = 0.d0
        do i = 1,nparam
            tot_wave = tot_wave + wave_od(i)**2*dx
        enddo
        do i = 1,nparam
            wave_od(i) = wave_od(i)/sqrt(tot_wave)
        enddo
        pt1 = (wave_od(2)-wave_od(1))/(rod1(2)-rod1(1))
        ptn =(wave_od(nparam)-wave_od(nparam-1))/(rod1(nparam)-
     1  rod1(nparam-1))
        call spline(rod1,wave_od,nparam,pt1,ptn,feed_od)
        return
        end         

        subroutine calc_first(r,first,second,k,nmono,j,sig_bend)
        implicit real*8(a-h,o-z)
        parameter (nparam=500)
        parameter (nvibw=3)
        common /wa1/roh1(nparam),wave_oh(nparam),feed_oh(nparam)
        common /wa2/rod1(nparam),wave_od(nparam),feed_od(nparam)
        dimension test(5),r(nvibw,nmono)
        dr = 0.0001
        do i = 1,5
            test(i) = 0.d0
        enddo
        do i = 1,5
            if (i.eq.1) then
                r(k,j) = r(k,j) + 2*dr
                call calc_psi(1,r(:,:),sig_bend,nmono,psi)
                test(1) = psi
                r(k,j)= r(k,j) - 2*dr
            endif
            if (i.eq.2) then
                r(k,j) = r(k,j) + dr
                call calc_psi(1,r(:,:),sig_bend,nmono,psi)
                test(2) = psi
                r(k,j)= r(k,j) - dr
            endif
            if (i.eq.3) then
                r(k,j) = r(k,j) 
                call calc_psi(1,r(:,:),sig_bend,nmono,psi)
                test(3) = psi
                r(k,j)= r(k,j) 
            endif
            if (i.eq.4) then
                r(k,j) = r(k,j) - dr
                call calc_psi(1,r(:,:),sig_bend,nmono,psi)
                test(4) = psi
                r(k,j)= r(k,j) + dr
            endif
            if (i.eq.5) then
                r(k,j) = r(k,j) - 2*dr
                call calc_psi(1,r(:,:),sig_bend,nmono,psi)
                test(5) = psi
                r(k,j)= r(k,j) + 2*dr
            endif
        enddo
        first = ((-test(1))+(8*test(2))+(-8*test(4))+(test(5)))/
     1  (12*(dr))
        second = ((-test(1))+(16*test(2))+(-30*test(3))+(16*test(4))+
     1  (-test(5)))/(12*(dr**2))
            
        return
        end subroutine

