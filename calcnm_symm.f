        subroutine calc_nm(y,freq,evec)
        implicit real*8(a-h,o-z)
        parameter (nsymm=96)
        parameter (natoms=9)
        parameter (ndim=27)
        parameter (nvib=(3*natoms)-6)
        dimension x(3,natoms,nsymm),y(ndim,nsymm),freq(nvib,nsymm),
     1  evec(ndim,nvib,nsymm)
        call get_symm(x)
        do i = 1,nsymm
            l = 0
            do j = 1,natoms
                do k = 1,3
                    l = l + 1
                    y(l,i) = x(k,j,i)/0.52917721067d0
                enddo
            enddo
        enddo
        do i = 1,nsymm
            call calc_harmonic_frequencies(y(:,i),ndim,freq(:,i),
     1      evec(:,:,i))
        enddo
        end subroutine
