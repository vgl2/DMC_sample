C        program hist
        subroutine hist(x,nbin,nwalk,box,p,xmin,xmax)
C       This program will be able to generate a histogram data-like file
C       on gnuplot.
        implicit real*8 (a-h,o-z)
        parameter (ndim=6)
        parameter (nmax=25000)
        dimension box(nbin),x(nmax),vector(nbin),p(nbin)
C        common/pts/psips(ndim,nmax)
C       Initialization       
        bin = nbin
        boxwidth = (xmax-xmin)/bin
        box = 0
C       look over walkers instead of bins here       
        do i = 1,nwalk
C            do j = 1,ndim
                ibin = int((x(i)-xmin)/boxwidth) + 1
                box(ibin) = box(ibin) + 1
            if (ibin.gt.nbin) then
                 print *, 'LOOK HERE', ibin
                 print *, x(i)
            endif
        enddo
C       Center the binwidths
C        p = boxwidth/2
C       placing the values on the center of the histogram
        do i = 1,nbin
            if (i.eq.1) then
                p(i) = xmin+(boxwidth/2.0)
            else
                p(i) = p(i-1) + boxwidth
            endif
        enddo
        end subroutine
