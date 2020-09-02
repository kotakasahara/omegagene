
        parameter (nstmx=20)
        parameter (nbinmx=1000)

        dimension e(nbinmx,nstmx),pdf(nbinmx,nstmx)
        dimension w_dn(nstmx),w_up(nstmx)
        dimension erng(2,nstmx)
        dimension ndat(nstmx)

        dimension erng_n(2,nstmx)

        dimension ratio(nstmx) 
        dimension e_s(nbinmx,nstmx),pdf_s(nbinmx,nstmx)
        dimension ndat_s(nstmx) 
c*********************************************
        read(5,*) nfil
        print*,'#  N of files = ',nfil
        print*,'# '

c  Once I used a scaling factor.
c  But I do not use this quantity any more.
cc      read(11,*) sfac
cc      print*,'# scaling factor for pdf = ',sfac
c********************
c  Input pdf.
        do kk=1,nfil
          print*,'#  v-st No. = ',kk
          idev = 100 + kk
          ndat(kk)=0

          do ii=1,1000000
            read(idev,*,end=900) e(ii,kk),pdf(ii,kk)
            ndat(kk)=ndat(kk)+1
cc          pdf(ii,kk)=pdf(ii,kk)*sfac
            print*,'     ',e(ii,kk),pdf(ii,kk)
          enddo
900       continue

          print*,'#       N of points input above = ',ndat(kk)
          print*,'# '
        enddo
c********************
c  Get range and weight for each v-state.
c  Output the v-state ranges.

        print*,'#  range defined for each state.'
        print*,'#     with trancating fringes. '
        print*,'# '

        read(501,*)
        read(501,*)
        read(501,*)

        do kk=1,nfil
          read(501,*) 
          read(501,*) erng(1,kk),erng(2,kk)
          read(501,*) w_dn(kk),w_up(kk)
          print*,'#  '
          print*,'#  range = ',erng(1,kk),erng(2,kk)
          print*,'#  weight = ', w_dn(kk),w_up(kk)

          if(w_dn(kk).gt.1.0d0 .and. w_up(kk).gt.1.0) then
            print*,' ??????????????????????????'
            print*,' ? strange: weight > 1.0. ?'
            print*,' ??????????????????????????'
            stop
          endif

          if(w_dn(kk).lt.0.0d0 .and. w_up(kk).lt.0.0) then
            print*,' ??????????????????????????'
            print*,' ? strange: weight < 0.0. ?'
            print*,' ??????????????????????????'
            stop
          endif

          idotr=300+kk
          write(idotr,201) erng(1,kk),erng(2,kk)
        enddo
201     format(f12.4,2x,f12.4)
c**********
c  Calculate the entire range & output.

        rmin=10000000.0
        rmax=-10000000.0

        do kk=1,nfil
          if(erng(1,kk).lt.rmin) rmin=erng(1,kk)
          if(erng(2,kk).gt.rmax) rmax=erng(2,kk)
        enddo

        print*,'# '
        write(90,201) rmin,rmax
        print*,'#  entire range: ',rmin,rmax
        print*,'# '
c********************
c << Important >>  Do not use this part. I skip here.

c  I do not use this part, because the v-pdfs were shifted so that 
c  the v-pdfs overlaps well in v_distrib.

c  Modify the pdf's when the inter-v-state transition is not 1.0. 
c  Calculate recovering factors: ratio(*,*).                           
        print*,'#  ratio of distributions among v-states (log scale). '
        ratio(nfil)=0.0

c  Skip.
        goto 834

        print*,'# ratio(pdf/pdf) ',nfil,ratio(nfil)
        do kk=nfil,2,-1
          ratio(kk-1)= log(w_up(kk-1)/w_dn(kk))
          print*,'# log[pdf/pdf] ',kk-1,ratio(kk-1)
        enddo

834     continue
c**********
c  Output pdf for each v-state.
c    Set the level of the highest v-state pdf to 1.
c    Note that pdf is given in log scale.

        print*,'# '
        print*,'#  Output pdf for each v-state '
        print*,'#  with trancating data outside the v-state range. '
        print*,'# '

        rmult=0.0

        do kk=nfil,1,-1
          idot=200+kk
          rmult=rmult + ratio(kk)

          ic=0
          do ii=1,ndat(kk)
            ee=e(ii,kk)
            if( ee.lt.erng(1,kk) . or. ee.gt.erng(2,kk) ) goto 888
            ic=ic+1
            e_s(ic,kk)=e(ii,kk)
            pdf_s(ic,kk)=pdf(ii,kk) - rmult

            write(idot,612) e_s(ic,kk),pdf_s(ic,kk)

888         continue
          enddo

          ndat_s(kk)=ic

          print*,'#     N of data remaining = ',kk,ndat_s(kk)
        enddo
612     format(e15.7,2x,e15.7)
c**********
c  Output block ranges, which are concatenated v-state ranges.

        print*,'#'

        do kk=1,nfil
          i1=kk-1
          i2=kk
          i3=kk+1

          if(i1.eq.0) then
            i1=1
            i2=2
            i3=3
            erng_n(1,kk)=erng(1,i1)
            erng_n(2,kk)=erng(2,i3)
          endif

          if(i1.ne.0 .and. i3.ne.nfil+1) then
            erng_n(1,kk)=erng(1,i1)
            erng_n(2,kk)=erng(2,i3)
          endif

          if(i3.eq.nfil+1) then
            i1=nfil-2
            i2=nfil-1
            i3=nfil
            erng_n(1,kk)=erng(1,i1)
            erng_n(2,kk)=erng(2,i3)
          endif

          print*,'# concatenated range No = ',i1,i2,i3

          idot=400+kk
          write(idot,201) erng_n(1,kk),erng_n(2,kk)
        enddo
c*********************************************
         stop
         end
