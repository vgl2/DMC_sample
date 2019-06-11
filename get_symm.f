        subroutine get_symm(x)
        implicit real*8(a-h,o-z)
        parameter (natoms=9)
        parameter (natomw=3)
        parameter (nperm=8)
        parameter (nstruc=12)
        parameter (nsymm=nperm*nstruc)
        parameter (nmono=3)
        dimension coord(3,natoms),x(3,natoms,nsymm),
     1  npermut(natoms,nperm),coord_water(3,natomw,nmono),
     1  coord_iso(3,natomw,nmono,nstruc),nwater(nmono,nstruc),
     1  coord_work(3,natoms,nstruc)
        open(unit=7,file='eq_x-trimer-uud.dat',status='old')
        open(unit=8,file='h_combo.dat',status='old')
        open(unit=9,file='water_permut.dat',status='old')
        read(7,*) ((coord(j,i),j=1,3),i=1,natoms)
        do i = 1,nperm
            read(8,*) (npermut(j,i),j=1,natoms)
        enddo
        do i = 1,nstruc
            read(9,*) (nwater(j,i),j=1,nmono)
        enddo
        ip = 0
        do i = 1,nmono
            do j = 1,natomw
                ip = ip + 1
                do k = 1,3
                    coord_water(k,j,i) = coord(k,ip)
                enddo
            enddo
        enddo
        do i = 1,nstruc
            do j = 1,nmono
                do k = 1,natomw
                    do l = 1,3
                        coord_iso(l,k,j,i) = 
     1                  coord_water(l,k,abs(nwater(j,i)))
                    enddo
                enddo
            enddo
        enddo
        do i = 1,nstruc
            do j = 1,nmono
                if (nwater(j,i).lt.0) then
                    do k = 1,natomw
                        do l = 1,3
                            coord_iso(l,k,j,i) = -coord_iso(l,k,j,i)
                        enddo
                    enddo
                else
                    continue
                endif
            enddo
        enddo
        do i = 1,nstruc
            ip = 0
            do j = 1,nmono
                do k = 1,natomw
                    ip = ip + 1
                    do l = 1,3
                        coord_work(l,ip,i) = coord_iso(l,k,j,i)
                    enddo
                enddo
            enddo
        enddo
        ip = 0
        do i = 1,nstruc
            do j = 1,nperm
                ip = ip + 1
                do k = 1,natoms
                    do l = 1,3
                        x(l,k,ip) =coord_work(l,abs(npermut(k,j)),i)
                    enddo
                enddo
            enddo
        enddo
        ip = 0
        do i = 1,nstruc
            do j = 1,nperm
                ip = ip + 1
                do k = 1,natoms
                    if (npermut(k,j).lt.0) then
                        do l = 1,3
                            x(l,k,ip) = -x(l,k,ip)
                        enddo
                    else
                        continue
                    endif
                enddo
            enddo
        enddo
        return
        end subroutine
                
