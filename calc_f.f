        program calc_f
        implicit real*8(a-h,o-z)
        parameter (nbin=100)
        parameter (nparam=500)
        common/wa/roh1(nparam),wave(nparam),feed1(nparam)
        dimension rval(nbin),box(nbin),psi(nbin),f_val(nbin),
     1  f_norm(nbin)
        open(unit=9,file='noimp-hist.dat',status='old')
        open(unit=10,file='f_val-morse.dat',status='unknown')
        pi = dacos(-1.d0)
        call initial_wave()
        do i = 1,nbin
            read(9,*) rval(i),box(i)
        enddo
        h_weight = 1.00782503223d0
        o_weight = 15.99491461957d0
c        weight =  (h_weight*o_weight)/(h_weight+o_weight)
c        weight = weight*1822.88852962d0
c        freq = 3887.5378889312051/219474.6
c        re = 1.810055141
c        do i = 1,nbin
c            psi(i) = ((weight*freq)/pi)**(0.25)*exp(-0.5*(weight*freq)*
c     1      ((rval(i)-re)**2))
c        enddo
        do i = 1,nbin
            call splint(roh1(:),wave(:),feed1(:),nparam,rval(i),psi(i))
        enddo
        do i = 1,nbin
            f_val(i) = psi(i)*box(i)
        enddo
        dx = rval(2)-rval(1)
        tot_int = 0.d0
        do i = 1,nbin
            tot_int = tot_int + dx*(f_val(i)**2)
        enddo
        do i = 1,nbin
            f_norm(i) = f_val(i)/sqrt(tot_int)
        enddo
        do i = 1,nbin
            write(10,*) rval(i),f_val(i),f_norm(i)
        enddo
        avg_r = 0.d0
        do i = 1,nbin
            avg_r = avg_r + rval(i)*f_norm(i)*dx
        enddo
        end program


        subroutine initial_wave()
        implicit real*8(a-h,o-z)
        parameter (nparam=500)
        common/wa/roh1(nparam),wave(nparam),feed1(nparam)
        dimension harm(nparam)
        open(unit=200,file='roh_morse.dat',status='old')
        do i = 1,nparam
            read(200,*) roh1(i),harm(i),wave(i)
        enddo
        pt1 = (wave(2)-wave(1))/(roh1(2)-roh1(1))
        ptn =(wave(nparam)-wave(nparam-1))/(roh1(nparam)-roh1(nparam-1))
        call spline(roh1,wave,nparam,pt1,ptn,feed1)
        return
        end         
