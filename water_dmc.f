        program qmc
        implicit real*8 (a-h,o-z)
C
C

        parameter (ndim =27)
        parameter (nmax = 25000)
        parameter (nparam=500)
        parameter (nvib=9)
        parameter (nvibw=3)
        parameter (nmono=3)
        common/pts/psips(ndim,nmax),sig_bend,t(nmax),v(nmax)
        common/kin/weight(ndim),d(ndim)
        common/pcoef/alpha,c2,re,rea,c1
        common/wa1/roh1(nparam),wave_oh(nparam),feed_oh(nparam)
        common/wa2/rod1(nparam),wave_od(nparam),feed_od(nparam)
        common/wav/psi_old(nmax),r(nvibw,nmono,nmax)
        common/der/drdx(nvib,ndim,nmax),dr2dx2(nvib,ndim,nmax)
        dimension step(nmax*ndim),fate(nmax),
     1  q_force_new(nvibw,nmono,nmax),psi_new(nmax),drift(nmax*ndim),
     1  g_top(nmax),old_psips(ndim,nmax),psi_sq_new(nmax),
     1  psi_sq_old(nmax),g_bottom(nmax),w(nmax),acceptance_ratio(nmax),
     1  acceptance(nmax),psi_old2(nmax),roh(nmax),roh2(nmax),
     1  q_f_old(ndim,nmax),q_f_new(ndim,nmax),t_old(nmax),angle(nmax),
     1  v_old(nmax),r_old(nvibw,nmono,nmax),drdx_new(nvib,ndim,nmax),
     1  dr2dx2_new(nvib,ndim,nmax),q_force_old(nvibw,nmono,nmax),
     1  first(nvib,nmax),second(nvib,nmax),first_new(nvib,nmax),
     1  second_new(nvib,nmax)
        character (len=1024) energy,wavefunction,wf
C
        read(*,*) dt
        read(*,*) endtime
        read(*,*) dseed
        read(*,*) wavefunction
        read(*,*) energy
        open (unit=12,file=energy,status='unknown')
        open (unit=20,file=wavefunction,status='unknown')
C
C
C       initialization
        sig_stretch = 1.d0
        sig_bend = 1.d0
        pi = dacos(-1.d0)
        call initial_roh_wave(sig_stretch)
        call initial_rod_wave(sig_stretch)
        h_weight = 1.00782503223d0
        o_weight = 15.99491461957d0
        re = 1.81005d0
        c1 = 1655.93078643302/219474.6
        c2 = 3887.5378889312051/219474.6
        rea = (104.508)*(pi/180.d0)
        g = (1/(h_weight*(re**2)))+(1/(h_weight*(re**2)))+
     1  (1/(o_weight*(re**2)))+(1/(o_weight*(re**2)))-
     1  ((2*dcos(rea))/(o_weight*re*re))
        g = (1/g)*1822.88852962d0
        alpha =g*c1
        call setup(dseed,n0)
        pi = dacos(-1.d0)
        n = n0
        time = 0.
c        dseed = 5.d9
        nborn = 0
        ndie = 0
c        dt = 10
        nstart = n
        call calc_r(n,psips,nmono,r)
        call calc_psi(n,r,sig_bend,nmono,psi_old)
        call calc_grad(n,r,sig_bend,psi_old,nmono,q_force_old,first,
     1  second)
c       calculate b matrix to get derivative in cartesian coordinates
        do i = 1,n
            call calc_b_mat(psips(:,i),q_force_old(:,:,i),q_f_old(:,i),
     1      nmono,drdx(:,:,i),dr2dx2(:,:,i))
        enddo
        veff = vavg(n0,first,second)
        do i = 1,n
            t_old(i) = t(i)
            v_old(i) = v(i)
        enddo
        vbar = veff
        alp = 0.5/dt
        write(20,*) 10
        do while (time.lt.endtime)
c            write(*,*) 'timestep = ',time
            ntot = n
            time = time + dt
            vtot = 0.
            sig = sqrt(dt)
            call gauss(0.d0,sig,step,ndim*nmax,dseed)
            naccept = 0
            ip = 0
c           Calculate drift term and move the walker
            do i = 1,n
                do j = 1,ndim
                    ip = ip + 1
                    step(ip) = step(ip)/sqrt(weight(j))
                    old_psips(j,i) = psips(j,i)
                    drift(ip) = (q_f_old(j,i)*dt)/(2*weight(j))
                    step(ip) = step(ip) + drift(ip)
                    psips(j,i) = psips(j,i) + step(ip)
                enddo
            enddo
            call ggubs(dseed,nmax,acceptance)
C           Calculate new Psi_t term to determine if walker can make the
C           move.
            do i = 1,n
                psi_new(i) = 1.d0
                g_top(i) = 1.d0
                g_bottom(i) = 1.d0
            enddo
            call calc_r(n,psips,nmono,r)
            call calc_psi(n,r,sig_bend,nmono,psi_new)
            call calc_grad(n,r,sig_bend,psi_new,nmono,q_force_new,
     1      first_new,second_new)
            do i = 1,n
                call calc_b_mat(psips(:,i),q_force_new(:,:,i),
     1          q_f_new(:,i),nmono,drdx_new(:,:,i),dr2dx2_new(:,:,i))
            enddo
            do i = 1,n
                psi_sq_old(i) = psi_old(i)**2
                psi_sq_new(i) = psi_new(i)**2
            enddo
C           Calculating the green's function value to determine if
C           movement will be made
            do i = 1,n
                do j = 1,ndim
                    g_top(i) = g_top(i)*exp(-((old_psips(j,i)-psips(j,i)
     1              -(q_f_new(j,i)*dt)/(2*weight(j)))**2)/((2*dt)/
     1              weight(j))) 
                    g_bottom(i) = g_bottom(i)*exp(-((psips(j,i)-
     1              old_psips(j,i)-(q_f_old(j,i)*dt)/(2*weight(j)))**2)/
     1              ((2*dt)/weight(j)))
                enddo
            enddo
            do i = 1,n
               w(i)=(psi_sq_new(i)*g_top(i))/(psi_sq_old(i)*g_bottom(i))
                acceptance_ratio(i) = min(w(i),1.d0)
C               If movement is accepted then copy everything the quantum
C               force, and value of psi along with the walker, if not,
C               then change the walker back to its previous position.
                if (acceptance_ratio(i).gt.acceptance(i)) then
                    naccept = naccept + 1
                    do j = 1,nmono
                        do k = 1,nvibw
                            q_force_old(k,j,i) = q_force_new(k,j,i)
                        enddo
                    enddo
                    do j = 1,nvib
                        first(j,i) = first_new(j,i)
                        second(j,i) = second_new(j,i)
                    enddo
                    drdx(:,:,i) = drdx_new(:,:,i)
                    dr2dx2(:,:,i) = dr2dx2_new(:,:,i)
                    do j = 1,ndim
                        q_f_old(j,i) = q_f_new(j,i)
                    enddo
                    psi_sq_old(i) = psi_sq_new(i)
                    psi_old(i) = psi_new(i)
                    t_old(i) = t(i)
                    v_old(i) = v(i)
                else 
                    do j = 1,ndim
                        psips(j,i) = old_psips(j,i)
                    enddo
                    do j = 1,nvib
                        first_new(j,i) = first(j,i)
                        second_new(j,i) = second(j,i)
                    enddo
                    t_old(i) = t_old(i)
                    v_old(i) = v_old(i)
                endif
            enddo
C           Calculate the new dt effective to use for branching.
            neff = n
            dt_eff = (dfloat(naccept)/dfloat(neff))*dt
            nborn = 0
            ndie = 0
            call ggubs(dseed,nmax,fate)
            n = ntot
            call calc_t(n,psips,psi_old,sig_bend,nmono,first,second,t)
            do i = 1,n
                do j = 1,ndim
                    psips(j,i) = psips(j,i)*0.52917721067d0
                enddo
            enddo
            do i = 1,n
                call calcpot(nmono,v(i),psips(:,i))
                v(i) = v(i)/627.509474d0
            enddo
            do i = 1,n
                do j = 1,ndim
                    psips(j,i) = psips(j,i)/0.52917721067d0
                enddo
            enddo
            do i = 1,n
                etot = t(i) + v(i)
                etot_old = t_old(i) + v_old(i)
                dv = etot - veff
                if (dv.lt.0) then
                    Pb = exp(-dv*dt) -1.
                    nb = int(Pb)
                    vtot = vtot + etot
                    if (fate(i).le.Pb-float(nb)) then
                        nb = nb + 1
                        nborn = nborn + 1
                    endif
                    do k = 1,nb
                        ntot = ntot + 1
                        psi_sq_old(ntot) = psi_sq_old(i)
                        psi_old(ntot) = psi_old(i)
                        t_old(ntot) = t_old(i)
                        t(ntot) = t(i)
                        v_old(ntot) = v_old(i)
                        v(ntot) = v(i)
                        drdx(:,:,ntot) = drdx(:,:,i)
                        dr2dx2(:,:,ntot) = dr2dx2(:,:,i)
                        do j = 1,ndim
                            psips(j,ntot) = psips(j,i)
                            q_f_old(j,ntot) = q_f_old(j,i)
                        enddo
                        do j = 1,nvib
                            first(j,ntot) = first(j,i)
                            second(j,ntot) = second(j,i)
                        enddo
                        do j = 1,nmono
                            do l = 1,nvibw
                                q_force_old(l,j,ntot)=q_force_old(l,j,i)
                            enddo
                        enddo
                        vtot = vtot + etot
                    enddo
                else
                    if (dv.gt.0) then
                        Pd = 1. - exp(-dv*dt)
                        if (fate(i).ge.Pd) then
                            vtot = vtot + etot
                        else
                            psips(1,i) = -999999.
                            ndie = ndie + 1
                        endif
                    endif
                endif
            enddo
            nt = n
            birth_rate = dfloat(nborn)/dfloat(nt)
            death_rate = dfloat(ndie)/dfloat(nt)
            n = 0
            do i = 1,ntot
                if (psips(1,i).gt.-99998.) then
                    n = n + 1
                    do j = 1,ndim 
                        psips(j,n) = psips(j,i)
                        q_f_old(j,n) = q_f_old(j,i)
                    enddo
                    do j = 1,nmono
                        do k = 1,nvibw
                            q_force_old(k,j,n) = q_force_old(k,j,i)
                        enddo
                    enddo
                    do j = 1,nvib
                        first(j,n) = first(j,i)
                        second(j,n) = second(j,i)
                    enddo
                    drdx(:,:,n) = drdx(:,:,i)
                    dr2dx2(:,:,n) = dr2dx2(:,:,i)
                    psi_sq_old(n) = psi_sq_old(i)
                    psi_old(n) = psi_old(i)
                    t_old(n) = t_old(i)
                    t(n) = t(i)
                    v_old(n) = v_old(i)
                    v(n) = v(i)
                endif
            enddo
            if ((time.gt.25000).and.mod(time,2500.).eq.0) then
c            if (time.eq.200000) then
                write(20,*) n,n0,time
                do i = 1,n
                    write(20,*) (psips(j,i),j=1,ndim)
                enddo
            endif
            veff_old = veff
            ntot = n
            vtot = vtot /dfloat(n)
            vbar = vtot
            veff = vbar + alp*log(dfloat(n0)/dfloat(n))
            write(12,*)time,veff*219474.6,n,dfloat(naccept)/dfloat(neff)
c            print*, time,veff*219474.6,n,dfloat(naccept)/dfloat(neff)
        enddo
        stop 
        end
C
C
        subroutine setup(dseed,n0)
        implicit real*8(a-h,o-z)
        parameter (ndim=27)
        parameter (nmax=25000)
        parameter (natoms=9)
        parameter (nsymm=96)
        common/pts/psips(ndim,nmax),sig_bend,t(nmax),v(nmax)
        common/kin/weight(ndim),d(ndim)
        common/coef/s2pi
        common/pcoef/alpha,c2,re,rea,c1
        dimension x(3,natoms,nsymm),y(ndim,nsymm)
C
C       initialization
        call calc_structures(dseed,n0,psips)
        pi = dacos(-1.d0)
C       m represents weight
        do j = 1,ndim
            if (j.lt.4) then
                weight(j) = 15.99491461957d0
            else if ((j.gt.9).and.(j.lt.13)) then
                weight(j) = 15.99491461957d0
            else if ((j.gt.18).and.(j.lt.22)) then
                weight(j) = 15.99491461957d0
            else
                weight(j) = 1.00782503223d0
            endif
            weight(j) = weight(j)*1822.88852962d0
        enddo
C
C       calculate the constants
        pi = dacos(-1.d0)
        s2pi = 1./sqrt(2.*pi)
        return
        end
C
C    
        function vavg(n,first,second)
        implicit real*8 (a-h,o-z)
        parameter (ndim=27)
        parameter (nvib=9)
        parameter (nvibw=3)
        parameter (nmax=25000)
        parameter (nparam=500)
        parameter (nmono=3)
        common/wav/psi_old(nmax),r(nvibw,nmono,nmax)
        common/pts/psips(ndim,nmax),sig_bend,t(nmax),v(nmax)
        common /wa1/roh1(nparam),wave_oh(nparam),feed_oh(nparam)
        common /wa2/rod1(nparam),wave_od(nparam),feed_od(nparam)
        common/pcoef/alpha,c2,re,rea,c1
        common/der/drdx(nvib,ndim,nmax),dr2dx2(nvib,ndim,nmax)
        dimension eloc(nmax),first(nvib,nmax),second(nvib,nmax)
C
        v = 0.
        call calc_r(n,psips,nmono,r)
        call calc_psi(n,r,sig_bend,nmono,psi_old)
        call calc_t(n,psips,psi_old,sig_bend,nmono,first,second,t)
        do i = 1,n
            do j = 1,ndim
                psips(j,i) = psips(j,i)*0.52917721067d0
            enddo
        enddo
        do i = 1,n
            call calcpot(nmono,v(i),psips(:,i))
            v(i) = v(i)/627.509474d0
        enddo
        do i = 1,n
            do j = 1,ndim
                psips(j,i) = psips(j,i)/0.52917721067d0
            enddo
        enddo
        vavg = 0.d0
        do i = 1,n
            eloc(i) = t(i)+v(i)
            vavg = vavg + eloc(i)
        enddo
        vavg = vavg/dfloat(n)
        print *, vavg*219474.6
        return
        end
C
C
