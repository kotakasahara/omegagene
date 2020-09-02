c
c
c**********************************************************
        implicit real*8 (a-h,o-z)
        dimension coef1(20),coef2(20)
        dimension coef_mix(20)
c**********************************************************
        print*,' #######################'
        print*,' # Mix two polynomials #'
        print*,' #######################'

        do ii=1,20
          coef1(ii)=0.0d0
          coef2(ii)=0.0d0
        enddo
c********
        read(5,*) weight

        print*,' '
        print*,'    weight = ',weight
        print*,'        poly = poly(previous) + weight * deviation '
        print*,' '

        read(40,*) elow2,eup2
        print*,'  range: ',elow2,eup2
c********
        read(70,*) iyesno
        if(iyesno .eq.0) then
          print*,'  The first md. '
          goto 833
        endif
        if(iyesno .eq.1) then
          print*,'  Not the first md. '
        endif
c********
        print*,'#  (1) Input polynomial having been used for md.'

        read(50,*) iord1

        write(6,200) iord1
200     format('#  order of fitting parameters = ',i2)

        if(iord1.gt.12) then
          print*,' ???????????????????????????????'
          print*,' ? iord1 is too large.         ?'
          print*,' ? Chenge subroutine genedens. ?'
          print*,' ???????????????????????????????'
          stop
        endif
c**
        do ii=1,iord1+1
          read(50,*) coef1(ii)
          write(6,144) ii,coef1(ii)
144       format("#  coef = ",i2,2x,e22.15)
        enddo
c**
c  grad_low & grad_up is not used.
        read(50,*) grad_low
        read(50,*) grad_up

        read(50,*) d1,d2
        write(6,149) d1,d2
149     format("#  Elow & Eup (not used) = ",e22.15,3x,e22.15)
c  temp is not used.
        read(50,*) temp
c*********
833     continue

        print*,' '
        print*,'#  (2) Input polynomials for d ln(Pmc) / dE.'
        read(51,*) iord2

        write(6,201) iord2
201     format('#  order of fitting parameters = ',i2)

        if(iord2.gt.12) then
          print*,' ???????????????????????????????'
          print*,' ? iord2 is too large.         ?'
          print*,' ? Chenge subroutine genedens. ?'
          print*,' ???????????????????????????????'
          stop
        endif
c**
        do ii=1,iord2+1
          read(51,*) coef2(ii)
          write(6,145) ii,coef2(ii)
145       format("#  coef = ",i2,2x,e22.15)
        enddo

c  grad_low & grad_up is not used.
          read(51,*) grad_low
          read(51,*) grad_up

          read(51,*) d1,d2
          write(6,169) d1,d2
169       format("#  Elow & Eup (not used) = ",e22.15,3x,e22.15)
c  temp is not used.
        read(51,*) temp
c****
c  The fitting range should be the same.
        if(elow1.ne.elow2) then
          print*,' ?????????????????????????????????'
          print*,' ? Lower fitting range mismatch.? '
          print*,' ?????????????????????????????????'
cc          stop
        endif
        if(eup1.ne.eup2) then
          print*,' ?????????????????????????????????'
          print*,' ? upper fitting range mismatch. ?'
          print*,' ?????????????????????????????????'
cc          stop
        endif
c*********
c  Mix.
        print*,' '
        print*,'#  (3) generated polynomials for dln{Pmc(E)]/dE.'

        do ii=1,20
          coef_mix(ii)=coef1(ii) + coef2(ii)*weight
        enddo

        iord=max(iord1,iord2)
        write(60,101) iord
        write(6,201) iord
101     format(i2)

        do ii=1,iord+1
          write(60,102) coef_mix(ii)
          write(6,144) ii,coef_mix(ii)
102       format(e22.15)
        enddo
        write(60,102) grad_low
        write(60,102) grad_up
        write(60,103) elow2,eup2
103     format(f16.5,1x,f16.5)
        write(60,104) temp
104     format(f6.1)

        write(6,169) elow2,eup2
c**********************************************************
        stop
        end
