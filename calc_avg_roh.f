        program calc_roh
        implicit real*8(a-h,o-z)
        parameter (nwavetot=20)
        parameter (nmax=25000)
        parameter (ndim=9)
        parameter (nbin=100)
        dimension n(nwavetot),n0(nwavetot),time(nwavetot),roh(nmax),
     1  psips(ndim,nmax),weight(nmax),tot_box(nbin),box(nbin),rval(nbin)
        open(unit=8,file='wf.dat',status='old') 
        open(unit=9,file='wf-hist.dat',status='unknown')
        do i = 1,nbin
            tot_box(i) = 0.d0
            box(i)  = 0.d0
        enddo
        read(8,*) nwave
        do k = 1,nwave
            print *, k
            read(8,*) n(k),n0(k),time(k)
            do i = 1,n(k)
                read(8,*) (psips(j,i),j=1,ndim)
                weight(i) = 1.d0
            enddo
            do i = 1,n(k)
                roh(i) = sqrt(((psips(4,i)-psips(1,i))**2)+((psips(5,i)-
     1          psips(2,i))**2)+((psips(6,i)-psips(3,i))**2))
            enddo
            call wavefunction_cont(roh(:),weight(:),n(k),rval(:),box(:),
     1      dx,nbin)
            do i = 1,nbin
                tot_box(i) = tot_box(i) + box(i)
                box(i) = 0.d0
            enddo
        enddo
        do i =1,nbin
            tot_box(i) = tot_box(i)/dfloat(nwave)
        enddo
        do i = 1,nbin
            write(9,*) rval(i),tot_box(i)
        enddo
        end program
