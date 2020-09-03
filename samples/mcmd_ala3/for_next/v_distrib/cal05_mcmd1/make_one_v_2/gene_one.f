c
        dimension e(10000,10),p(10000,10)
        dimension er(2,10)
        dimension icou(10000),ic(10)
        dimension pmean(10000)
c***********************************************************************
        read(5,*) de
        read(5,*) nvst
        do kk=1,nvst
          read(5,*) er(1,kk),er(2,kk)
        enddo
  
        print*,'  Parameter input: '
        print*,' DE = ',de
        print*,' overlap range = ',over_low,over_u
        print*,' N of v states = ',nvst
        do kk=1,nvst
          print*,'  E range = ',kk,er(1,kk),er(2,kk)
        enddo
        print*,' '
c********
        do ii=1,nvst
          ic(ii)=0
          idev=100+ii

          do jj=1,10000
            read(idev,*,end=901) e(jj,ii),p(jj,ii)
            p(jj,ii)=exp(p(jj,ii))
            ic(ii)=ic(ii)+1
          enddo

901       continue
        enddo
c**
        emin=1000000.0
        do jj=1,ic(1)
          if(e(jj,1) .lt.emin) emin=e(jj,1)
        enddo

        emax=-1000000.0
        do jj=1,ic(nvst)
          if(e(jj,nvst) .gt. emax) emax=e(jj,nvst)
        enddo

        print*,'  Emin & Emax (detected) = ',emin,emax
        print*,'  N od data in reach pdf:'
        do ii=1,nvst
          print*,'  N = ',ii,ic(ii)
        enddo
c********
        ndat=0
        do ii=1,1000000
          ee=(ii-1)*de+emin
          if(ee.gt.emax) goto 990
          ndat=ndat+1
        enddo

990     continue
        print*,' N o points = ',ndat
c**
        do ii=1,ndat
          icou(ii)=0
          pmean(ii)=0.0
        enddo
c********
        do ii=1,ndat
          ee=(ii-1)*de+emin

          do kk=1,nvst

            if(kk.ne.1 .and. kk.ne.nvst) then
              do n=1,ic(kk)
                if(e(n,kk).lt.er(1,kk)) goto 333
                if(e(n,kk).gt.er(2,kk)) goto 333
                if(e(n,kk).eq.ee) then
                  icou(ii)=icou(ii)+1
                  pmean(ii)=pmean(ii)+p(n,kk)
                endif
333             continue
              enddo
            endif

            if(kk.eq.1) then
              do n=1,ic(kk)
                if(e(n,kk).gt.er(2,kk)) goto 444
                if(e(n,kk).eq.ee) then
                  icou(ii)=icou(ii)+1
                  pmean(ii)=pmean(ii)+p(n,kk)
                endif
444             continue
              enddo
            endif

            if(kk.eq.nvst) then
              do n=1,ic(kk)
                if(e(n,kk).lt.er(1,kk)) goto 555
                if(e(n,kk).eq.ee) then
                  icou(ii)=icou(ii)+1
                  pmean(ii)=pmean(ii)+p(n,kk)
                endif
555             continue
              enddo
            endif
          enddo
        enddo
c********
c  Output.

        do ii=1,ndat
          ee=(ii-1)*de+emin
          gg=pmean(ii)/real(icou(ii))
          gg=log(gg)
          write(21,221) ee,gg
221       format(e15.4,2x,e15.5)
          print*,ee,gg
        enddo

c***********************************************************************
        stop
        end
