c
c  The discontinuous ln[n(E)] --> continuous ln[n(E)].
c
c*******************************************************************
        implicit real*8 (a-h,o-z)
 
        dimension erange(2,50)

        parameter (ndat=100000)
        dimension xx(ndat,50),yy(ndat,50)

        dimension numd(50)
        dimension amove(50)
c*******************************************************************
c  Set general parameter(s).
        rgas=8.31451d0/4.184d0/1000.0d0
c*******************************************************************
c  Input N of virtual states.

        read(10,*) ninp

        print*,'  N of input files (i.e., N of virtual states) ',ninp
        print*,' '
c**********
c  Range input.
c    Only the range is input. The coefficients are not input.

        print*,'  range input: '
        do ii = 1,ninp
          idev = 200 + ii
          read(idev,*) idum

          do kk = 1,idum+1
            read(idev,*) 
          enddo

          read(idev,*)
          read(idev,*)

          read(idev,*) erange(1,ii),erange(2,ii)
          print*,'      range: ',ii,erange(1,ii),erange(2,ii)
        enddo
        print*,' '
c**********
c  Input ln[P].
c    The ln[P]s of v-states have discontinuity.

        print*,'  Input pdfs for v-states. '

        do kk=1,ninp
          iinp= 100 + kk
          icou=0
700       continue

          read(iinp,*,end=900) x,qdum,y
          icou=icou+1
          xx(icou,kk)=x 
          yy(icou,kk)=y 

          goto 700

900       continue
          numd(kk)=icou

          print*,' '
          print*,'  N of data points = ',kk,numd(kk)
        enddo
c**********
c  The discontinuity is reset. --> continuous funtion.

        amove(1)=0.0d0

        do kk=2,ninp
          av1=0.0d0
          av2=0.0d0
          icou=0

          do mm=1,numd(kk)
            do nn=1,numd(kk-1)
              dx=xx(nn,kk-1)-xx(mm,kk)
              dx=abs(dx)
              if(dx.lt.0.00001d0) then
                icou=icou+1
                av1=av1+yy(nn,kk-1)
                av2=av2+yy(mm,kk)
              endif
            enddo
          enddo

          av1 = av1/real(icou)
          av2 = av2/real(icou)
          amove(kk) = av1-av2

          print*,'  v-state = ',kk
          print*,'     icou = ',icou
          print*,'     av1 = ',av1/real(icou)
          print*,'     av2 = ',av2/real(icou)
          print*,'     amove = ',amove(kk)

          do mm=1,numd(kk)
            yy(mm,kk)=yy(mm,kk)+amove(kk)
          enddo
        enddo
c****
c  Shift the connected ln[P].

        yymax=-10000000.0
        do kk=1,ninp
          do mm=1,numd(kk)
            if(yy(mm,kk).gt.yymax) yymax=yy(mm,kk)
          enddo
        enddo

        do kk=1,ninp
          do mm=1,numd(kk)
            yy(mm,kk)=yy(mm,kk)-yymax
          enddo
        enddo
c****
        do kk=1,ninp
          do mm=1,numd(kk)
            write(500,428) xx(mm,kk),yy(mm,kk)
            print*, xx(mm,kk),yy(mm,kk)
          enddo
        enddo

428     format(e16.8,2x,e20.10)
c****
c  For V-AUS.
cc      tt=300.0d0
cc      rt=rgas*tt
cc      do kk=1,ninp
cc        iot=300 + kk

cc        do mm=1,numd(kk)
cc          zz=-rt*yy(mm,kk)
cc          write(iot,261) xx(mm,kk),yy(mm,kk),zz
cc        enddo
cc      enddo
cc261     format(f8.3,2x,e12.4,2x,e12.4)
c**********
c  Smoothing of functions.

c  Not yet made.

cc      xmin=xx(1,1)
cc      xmax=xx(numd(ninp),ninp)

cc      do ii=1,1000000

cc      enddo
c*******************************************************************
        stop
        end
