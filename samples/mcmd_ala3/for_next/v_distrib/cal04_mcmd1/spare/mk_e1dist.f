c******************************************************
        implicit real*8 (a-h,o-z)
        parameter (nwin=20000)

        real*8    e(nwin),pdf(nwin)
        integer*4 itmp(nwin)
c******************************************************
c input pdf, including 0.00000 in pdf
        i=0
 900    continue
        i=i+1
        read(50,*,end=800) e(i),pdf(i)
        goto 900
 800    continue
        ndat=i-1

c        do i=1,ndat
c        write(60,600) e(i),pdf(i)
c        enddo

        ii=0
        do i=1,ndat
        if(pdf(i).eq.0d0.or.i.eq.1) then
        ii=ii+1
        itmp(ii)=i
cc        write(6,*) ii,itmp(ii)
        endif
        enddo
        ndat0=ii

        if(ndat0.eq.1) then
        write(*,*) 'All data are accepted.'
        do i=1,ndat
        write(60,600) e(i),pdf(i)
        enddo
        stop
        endif

        imax=-100000
        ii=0
        do i=2,ndat0
        if(itmp(i)-itmp(i-1).ge.imax) then
        imax=itmp(i)-itmp(i-1)
        ii=i-1
        endif
        enddo


        do i=itmp(ii)+1,itmp(ii+1)-1
        write(60,600) e(i),pdf(i)
        enddo

 600    format(e15.7,2x,e15.7)

c******************************************************
        stop
        end
