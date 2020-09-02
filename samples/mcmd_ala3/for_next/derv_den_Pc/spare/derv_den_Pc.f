c**********************************************************
        implicit real*8 (a-h,o-z)

        parameter (nmax=100000)
        dimension ene(nmax),pdf(nmax)
        dimension den(nmax),dervpc(nmax)
c**********************************************************
c  Input the polynomial parameters of Pc(E,T)].

        print*,'#  (1) Input ln[Pc(E,T)] & calculate ln[n(E)] '
        print*,'#      using temperature input. '
        print*,'#  (2) Calc. numerical1 derivatives of ln[pc(E,T)] '
        print*,'#      and output them.'
        print*,'# '
        print*,'#   NOTE: Only well-sampled energy range is output. '
        print*,'# '
c*********
        read(5,*) temp
c*********
        rt=temp*8.31451d0/4.184d0/1000.0d0

        print*,'# /////////////////////'
        print*,'# // temperature OK? //'
        print*,'# /////////////////////'

        print*,'#  T = ',temp
        print*,'#  RT (kcal / mol) = ',rt
        print*,'# ************  '
c********
c  Input ln[Pc(E,T)].
        ic=0

800     continue

        read(10,*,end=900) e,p
        ic=ic+1
        ene(ic)=e
        pdf(ic)=p

        goto 800
900     continue
        print*,'#  N of points = ',ic

        dd=ene(ic)-ene(1)
        if(dd.le.0.0) then
          print*,' ???????????????????????????????'
          print*,' ? Detected energy descending. ?'
          print*,' ???????????????????????????????'
          stop
        endif
c********
c  Detect the maximum of prob.
        pdfmx=-1000000.0
        nkeep=0
        do kk=1,ic
          if(pdf(kk).gt.pdfmx) then
            pdfmx=pdf(kk)
            etop=ene(kk)
            nkeep=kk
          endif
        enddo

c  Get the energy range for operation.
        pdfget=pdfmx-2.0d0

        do kk=1,nkeep
          if(pdf(kk).gt.pdfget) then
            elow=ene(kk)
            goto 822
          endif
        enddo
822     continue

        do kk=ic,nkeep,-1
          if(pdf(kk).gt.pdfget) then
            eup=ene(kk)
            goto 833
          endif
        enddo
833     continue

        write(6,203) pdfmx,etop
203     format('#  pdf(max) & Etop = ',f12.3,2x,f10.2)
        write(6,204) elow,eup
204     format('#     energy range: elow & eup = ',f12.3," to ",f12.3)
c********
c  calculate ln[n(E)].
        do kk=1,ic
          den(kk)=ene(kk)/rt + pdf(kk)
        enddo
c***
c  Numerical deribatives of ln[n(E)].
        dervpc(1)=(den(2)-den(1))/(ene(2)-ene(1))
        dervpc(ic)=(den(ic)-den(ic-1))/(ene(ic)-ene(ic-1))

        do kk=2,ic-1
          dervpc(kk)=(den(kk+1)-den(kk-1))/(ene(kk+1)-ene(kk-1))
        enddo
c***
        if(elow.le.ene(1)) then
          print*,' ??????????????????????????????????'
          print*,' ? Please reset: elow is too low. ?'
          print*,' ??????????????????????????????????'
          stop
        endif
        if(eup.ge.ene(ic)) then
          print*,' ???????????????????????????????????'
          print*,' ? Please reset: eup is too large. ?'
          print*,' ???????????????????????????????????'
          stop
        endif
c***
c  Output the derivatives in the range. 
        icc = 0

        do kk=1,ic
          if(ene(kk).ge.elow .and. ene(kk).le.eup) then
            icc = icc + 1
            print*,ene(kk),dervpc(kk),den(kk)
            write(20,200) ene(kk),dervpc(kk)
          endif
        enddo
200     format(e16.8,2x,e16.9)

        print*,'# N of data for output = ',icc
c**********************************************************
        stop
        end
