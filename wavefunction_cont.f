        subroutine wavefunction_cont(x,fam,n,rval,standard_box,
     1  dx,nbin)
c        program hist2
        implicit real*8 (a-h,o-z)
C       this program will take my coordinates and generate the
C       wavefunction. Victor Lee 7/12/2016
c        parameter (nbin=75)
        parameter (nmax=25000)
        parameter (ndim=1)
        dimension x(nmax),y1(nmax),nfam(nmax),rval(nbin),box(nbin),
     1  box_squared(nbin), standard_box(nbin),old_box(nbin),
     1  x_new(nmax),fam(nmax),psi_sq(nbin)
c        open(unit=9,file='descendant_cont_weight.dat',status='old')
c        open(unit=10,file='plot_descendant_cont.dat',status='unknown')
        do i = 1,nbin
            old_box(i) = 0.
            box(i) = 0.
        enddo
        descend = 0.d0
        do i =1,n
            descend = descend + fam(i)
        enddo
        rmin = 0.d0
        rmax = 4.d0
        do j = 1,nbin
            old_box(j) = box(j)
        enddo
        call hist_weight(x,nbin,n,box,rval,rmin,rmax,fam,
     1  descend,dx)
        do i = 1,nbin
            box(i) = old_box(i) + box(i) 
        enddo
        box_total = 0
        do i = 1,nbin
            box_squared(i) = box(i)**2
            box_total = box_total + box_squared(i)*dx
        enddo
        do i = 1,nbin
            standard_box(i) = box(i)/sqrt(box_total)
        enddo
        tot_box = 0
        do i = 1,nbin
            tot_box = tot_box + dx*standard_box(i)**2
        enddo
c        print *, tot_box, '1???'
        end 
            
