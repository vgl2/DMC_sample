        program average
C       this program will take the mean of the latter half of the
C       elements in an array and tell that number
        implicit real*8(a-h,o-z)
        real add,mean
        parameter (n=20000) !number of energy files
        dimension timestep(n),energy(n),pop(n)
        open (unit=9,file='energy.dat',status='old')
        do i = 1,n
                read(9,*) timestep(i),energy(i),pop(i) 
        enddo
        add = 0.d0
        do i = (n/2)+1,n
C                print *, energy(i),i
                add = add + energy(i)
C                print *, add
C                read(*,*)
        enddo
C        print *, add
        mean = add/(n/2)
        print *, mean
        end program
