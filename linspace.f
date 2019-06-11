        subroutine linspace(xmin,xmax,nbin,vector)
C       This subroutine will take the minimum and maximum of the two
C       numbers that were given and make a equally spaced vector between
C       xmin and xmax based on the amount of nbin specified by the
C       user.
C       Victor Lee May 13,2016
        implicit real*8 (a-h,o-z)
        dimension vector(nbin)
C       Figure out the range of x
        delta_x = (xmax-xmin)/nbin
C        print *, xmax, 'max'
c        print *, xmin, 'min'
C        print *, delta_x, 'delta'
C        read(*,*)
C       loop over number of nbin
        do i = 1,nbin
                if (i.eq.1) then
                        vector(i) = xmin + delta_x
                else
                        vector(i) = vector(i-1) + delta_x
                endif
        enddo
        end subroutine 
