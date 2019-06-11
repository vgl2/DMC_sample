        program historam
C       This program will be able to take the results from the DMC
C       simulation and generate a proper histogram file. 
C       Victor Lee 5/16/16
        implicit real*8(a-h,o-z)
        double precision box_sum,norm_constant,norm_constant_check
        parameter(nbin=20)
        parameter(xmin=-5)
        parameter(xmax=5)
        dimension box(nbin),box_sq(nbin),vector(nbin),box_norm(nbin),
     1  step(nbin)
        open(unit=8,file='hist.dat', status='old')
C        open(unit=9, file='normwf.dat',status='new')
C       read in dmc results
        do i = 1,nbin
                read(8,*) step(i),box(i)
                box_sq(i) = box(i)**2
        enddo
        box_sum = sum(box_sq)
        norm_constant = sqrt(box_sum)
        do i = 1,nbin 
                box(i) = box(i)/norm_constant
                print *, box(i) ,'normalized'
        enddo
C        norm = 0
        end program 
