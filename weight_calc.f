        program weight_calc
C       This program will make multiple continuous weighting dmc 
C       calculations on the same wavefunction with different seed 
C       values to get the average weight of the simulation.
        implicit real*8(a-h,o-z)
        parameter (nmax=25000)
        parameter (nwavetot = 50)
        parameter (nmaxtot = nmax*nwavetot)
        parameter (ndim = 27)
        parameter (natoms = 3)
        parameter (nvib=9)
        parameter (ndimw=9)
        parameter (nvibw=3)
        parameter (nmono=3)
        parameter (nparam=500)
        common /wa1/roh1(nparam),wave_oh(nparam),feed_oh(nparam)
        common /wa2/rod1(nparam),wave_od(nparam),feed_od(nparam)
        common/pcoef/alpha,c2,re,rea,c1
        dimension con_weight(nmax),n0(nwavetot),simtime(nwavetot),
     1  etrial(nwavetot),psips(ndim,nmax),n(nwavetot),dseed(nwavetot),
     1  freq_keep(nvibw),eq(3,natoms),y(ndim),req(nvibw)
        character(len=1024) psisq, holes, wavefunction, weight
        read(*,*) dt
        read(*,*) endtime
        read(*,*) psisq
        read(*,*) weight
        read(*,*) dseedi
        read(*,*) wavefunction
        read(*,*) holes
        open(unit=8,file=wavefunction,status='old')
        open(unit=300,file=psisq,status='unknown')
        open(unit=100,file=weight,status='unknown')
        pi = dacos(-1.d0)
        ntot = 0
        h_weight = 1.00782503223d0
        o_weight = 15.99491461957d0
        c1 = 1655.93078643302/219474.6
        c2 = 3887.5378889312051/219474.6
        open(unit=39,file='eq_x.dat',status='unknown')
        read(8,*) nwave
        wave = nwave
        call ggubs(dseedi,nwave,dseed)
        call cpu_time(start)
        tot_etrial = 0.d0
        read(39,*) ((eq(j,i),j=1,3),i=1,natoms)
        l = 0
        do i = 1,natoms
            do j = 1,3
                l = l + 1
                y(l) = eq(j,i)/0.52917721067
            enddo
        enddo
        call calc_r(1,y,nmono,req)
        re = req(1)
        rea = req(3)
        sig_stretch = 1.d0
        call initial_roh_wave(sig_stretch)
        call initial_rod_wave(sig_stretch)
        g = (1/(h_weight*(re**2)))+(1/(h_weight*(re**2)))+(1/(o_weight*
     1  (re**2)))+(1/(o_weight*(re**2)))-((2*dcos(rea))/(o_weight*re*
     1  re))
        g = (1/g)*1822.88852962d0
        alpha = g*c1
        do k = 1,nwave
            read(8,*) n(k),n0(k),simtime(k)
            do i = 1,n(k)
                read(8,*) (psips(j,i),j=1,ndim)
            enddo
c            if (k.eq.2) then
            call cont_weight(psips(:,:),dseed(k),n(k),n0(k),
     1      con_weight(:),etrial(k),dt,endtime,holes)
            print *, etrial(k)*219474.6,k
            write(100,*) n(k)
            do i = 1,n(k)
                write(100,*) con_weight(i)
            enddo
c            print *, etrial(k)*219474.6
c     1          nmax,nbin)
            tot_etrial = tot_etrial + etrial(k)
            call cpu_time(finish)
c            endif
        enddo
        avg_etrial = tot_etrial/wave
        st_avg_eref = 0.d0
        do k = 1,nwave
            st_avg_eref = st_avg_eref + (etrial(k)-avg_etrial)**2 
        enddo
        st_dev_eref = sqrt((1/dfloat(nwave))*st_avg_eref)
        print *, avg_etrial*219474.6,st_dev_eref*219474.6
c       Calculate <x> and <x^2>
        print *, finish-start,'cpu time'
        end program
