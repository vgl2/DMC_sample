        subroutine calc_structures(dseed,nfinal,psips)
        implicit real*8(a-h,o-z)
        parameter (nmax=25000)
        parameter (ndim=27)
        parameter (nvib=21)
        parameter (nsymm=96)
        parameter (natoms=9)
        dimension freq(nvib,nsymm),evec(ndim,nvib,nsymm),
     1  coord(ndim,nmax,nsymm),x(ndim,nsymm),ngood(nsymm),
     1  psips(ndim,nmax),coord2(3,natoms,nmax)
        call calc_nm(x,freq,evec)
        nfinal = 5000
        do i = 1,nsymm
            call mc_sample(dseed,x(:,i),nsymm,nfinal,freq(:,i),
     1      evec(:,:,i),ngood(i),coord(:,:,i))
        enddo
        ip = 0
        do i = 1,nsymm
            do j = 1,ngood(i)
                ip = ip + 1
                do k = 1,ndim
                    psips(k,ip) = coord(k,j,i)
                enddo
            enddo
        enddo
        do while (ip.lt.nfinal)
            do i = 1,nsymm
                do j = 1,ngood(i)
                    ip = ip + 1
                    do k = 1,ndim
                        psips(k,ip) = coord(k,j,i)
                    enddo
                    if (ip.eq.nfinal) then
                        exit
                    endif
                enddo
                if (ip.eq.nfinal) then
                    exit
                endif
            enddo
        enddo
        return
        end subroutine
