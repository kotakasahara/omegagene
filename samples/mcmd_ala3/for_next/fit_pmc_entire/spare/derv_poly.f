
        implicit real*8(a-h,o-z)
        dimension coef(20),ff(20)
c******************************************************
c  Input polynomial parameters.

        read(50,*) iord
c**
c  Used.
        print*,' ****************************************'
        do ii=1,iord+1
          read(50,*) coef(ii)
          write(6,144) ii,coef(ii)
144       format("#  coef (input) = ",i2,2x,e22.15)
        enddo
        print*,' '
c**
c  Not used.
        read(50,*) grad_low
        read(50,*) grad_up
c**
c  Not used.
        read(50,*) elow1,eup1
c       write(6,149) elow1,eup1
c**
c  Not used.
        read(50,*) temp
c****************
c  Input range used.
        read(60,*) elow1,eup1
        write(6,149) elow1,eup1
149     format("#  Elow & Eup = ",e12.4,3x,e12.4)
c****************
c  Input the squeezing factor.

        read(5,*) fac_squeeze
c****************
c Calc. derivatives in poly. form and output.

        print*,' '
        print*,' Below is output: '

        write(51,101) iord-1
101     format(i2)

        do kk=2,iord+1
          jj=kk-1
          cc = coef(kk)*real(jj)
          cc=fac_squeeze*cc
          write(51,102) cc
          write(6,102) cc
102       format(e22.15)

          ff(jj)=cc
        enddo

        write(51,102) grad_low
        write(51,102) grad_up

        write(51,103) elow1,eup1
103     format(f12.4,x,f12.4)
        write(51,104) temp
104     format(f6.1)
c********
        de=0.01d0
        do ii=1,1000000
          ee=elow1+(ii-1)*de + ecut_low
          if(ee.gt.eup1-ecut_up) goto 933

          pp=0.0d0
          do jj=1,iord
            pp=pp + ff(jj)*ee**(jj-1)
          enddo

          write(52,133) ee,pp
133       format(e17.10,2x,e22.15)
        enddo
933     continue
c******************************************************
        stop
        end
