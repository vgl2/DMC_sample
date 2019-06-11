        program hist2
        implicit real*8 (a-h,o-z)
C       this program will take my coordinates and generate the
C       wavefunction. Victor Lee 7/12/2016
        parameter (nbin=125)
        parameter (nmax=25000)
        parameter (ndim=1)
        dimension x1(nmax),y1(nmax),z1(nmax),x2(nmax),y2(nmax),z2(nmax),
     1  box(nbin),r(nmax),nfam(nmax),rval(nbin),box_squared(nbin),
     1  standard_box(nbin),old_box(nbin),analytical_soln(nbin),
     1  analytical_squared(nbin),analytical_box(nbin),dr(nmax)   
        open(unit=9,file='descendant.dat',status='old')
        open(unit=10,file='plot.dat',status='unknown')
        omega = 4401.d0/219474.6
        weight =0.5*(1.0078250/(6.022140857e23*9.10938291e-28))
        force = weight*omega**2
        re = 1.40104295 
        alpha = sqrt(force*weight)
        do i = 1,nbin
            old_box(i) = 0.
            box(i) = 0.

        enddo
        do k = 1,11
            read(9,*) n 
            do i = 1,n
                read(9,*) x1(i),y1(i),z1(i),x2(i),y2(i),z2(i),nfam(i)
                r(i)=sqrt((x2(i)-x1(i))**2+(y2(i)-y1(i))**2+(z2(i)-
     1  z1(i))**2)
                dr(i) = r(i)-re
C                print *, r(i)
            enddo
            nzero = 0
            zero = 0
            ndescend = 0
            descend = 0
            do i = 1,n
                if (nfam(i).lt.1) then
                    nzero= nzero + 1
                endif
                ndescend = ndescend + nfam(i)
            enddo
            zero = nzero 
            descend = ndescend
            print *, zero/descend
            rmin = -3.
            rmax = 3.
            do j = 1,nbin
                old_box(j) = box(j)
            enddo
            if (k.gt.1) then
                call hist(dr,nbin,n,box,rval,rmin,rmax)
                do i = 1,nbin
                    box(i) = old_box(i) + box(i) 
                enddo
            endif
            enddo
C                call linspace(rmin,rmax,nbin,rval)
                box_sum = 0
                box_total = 0
                analytical_total = 0
                do i = 1,nbin
                    box_squared(i) = box(i)**2
                    box_total = box_total + box_squared(i)
                    box_sum = box_sum + box(i)
                enddo
                do i = 1,nbin
                    analytical_soln(i) = (alpha)**(0.25)*exp(-0.5*alpha*
     1  (rval(i)**2))
                enddo
                do i = 1,nbin
                    analytical_squared(i) = analytical_soln(i)**2
                    analytical_total = analytical_total + 
     1  analytical_squared(i)
                    analytical_sum = analytical_sum + box(i)
                enddo
            print *, sqrt(box_total),'sqrt box tot'
            print *, box_sum,'sum box tot'
            do i = 1,nbin
                standard_box(i) = box(i)/sqrt(box_total)
                analytical_box(i) =analytical_soln(i)/
     1  sqrt(analytical_total)
                write(10,*), rval(i),standard_box(i),analytical_box(i)
            enddo
            tot_box = 0
            do i = 1,nbin
                tot_box = tot_box + standard_box(i)**2
            enddo
            print *, tot_box, '1???'
        end 
            
