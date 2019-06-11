        program descend
C       This program will be able to take the descendent weight of the
C       walkers to calculate properties of interest (<r>, <theta>, etc)
C       Victor Lee 7-11-16
        implicit real*8(a-h,o-z)
        parameter (natom=2)
        parameter (ndim=6)
        parameter (nmax=25000)
        parameter (nresult=20)
        dimension x1(nmax),y1(nmax),z1(nmax),x2(nmax),y2(nmax),z2(nmax),
     1  nfam(nmax),r(nmax),fam(nmax),avgr(nresult),e(nmax),avge(nresult)
        open(unit=9,file='descendant.dat',status='old')
        totavgr = 0
        totavge = 0
        totavgr_no = 0
        do k = 1,nresult
            read(9,*) n,time
C           print *, n
            ndescend = 0
            nzero = 0
            do i = 1,n
                read(9,*) x1(i),y1(i),z1(i),x2(i),y2(i),z2(i),nfam(i)
                r(i)=sqrt((x2(i)-x1(i))**2+(y2(i)-y1(i))**2+(z2(i)-
     1  z1(i))**2)
                e(i)=poten(r(i))
                if (nfam(i).eq.0) then
                    nzero = nzero + 1
                endif
                ndescend = ndescend + nfam(i)
            enddo
            avge(k) = 0
            avgr(k) = 0
            avgr_no= 0
            do i = 1,n
                fam(i) = nfam(i)
                avgr(k) = avgr(k) + ((fam(i)/ndescend)*r(i))  
                avgr_no = avgr_no + r(i)
                avge(k) = avge(k) + ((fam(i)/ndescend)*e(i)) 
            enddo
                print *, avgr_no/n ,'<r>'
                read(*,*) 
            totavgr = totavgr + avgr(k)
            totavge = totavge + avge(k)
        enddo
            print *, '<r> is', totavgr/nresult
            print *, '<e> is', (totavge/nresult)*219474.6
        end
        
        function poten(x)
        implicit real*8 (a-h,o-z)
C
        c = 4401.d0
        c = c/219474.6
        d = 1.0078250
        weight = d/(6.022140857e23*9.10938291e-28)
        red_weight = 0.5*weight
        re = 1.40104295d0
        poten =  0.5*red_weight*c**2*(x-re)**2
        return
        end
               
