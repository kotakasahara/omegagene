c
c*******************************************************************
        implicit real*8 (a-h,o-z)

        parameter (ndat=100000)
        dimension xx(ndat),yy(ndat)

        real*16 zz(ndat),eps
c*******************************************************************
c  Set general parameter(s).
        rgas=8.31451d0/4.184d0/1000.0d0
c*******************************************************************
c  Input T.

        read(5,*) temp
c************
c  Input ln[n(E)] for the whole range.

        icou=0

800     continue

        read(20,*,end=900) aaa,bbb
        icou=icou+1

        xx(icou)=aaa
        yy(icou)=bbb

        goto 800
c****
900     continue

         print*,' T(K) &  N of inpit data = ',temp,icou
c************
c  Generate ln[P(E,T)] & P(E,T).
c       xx <-- E
c       yy <-- ln[P(E,T)]
c       zz <-- P(E,T)

        rt=rgas * temp
        do ii=1,icou
          yy(ii)=yy(ii) - xx(ii)/rt
        enddo

        ymax=-100000.0
        do ii=1,icou
          if(yy(ii).gt.ymax) ymax=yy(ii)
        enddo
        do ii=1,icou
          yy(ii)=yy(ii)-ymax
        enddo

        do ii=1,icou
          zz(ii)=exp(yy(ii))
        enddo
c****
c  Output.

        eps=1.0e-30
        do ii=1,icou
          if(zz(ii).lt. eps ) zz(ii)=0.0d0
          write(40,133) xx(ii),yy(ii),zz(ii)
cc        print*,'  ',xx(ii),yy(ii),zz(ii)
        enddo

133     format(f12.3,2x,e20.8,2x,e20.5)
c*******************************************************************
        stop
        end
