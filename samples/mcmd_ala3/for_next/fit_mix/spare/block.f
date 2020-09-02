c
        implicit real*8 (a-h,o-z)
c
        dimension iord(50)
        dimension coef(20,50)
        dimension erange(2,50)

        parameter (nstmx=50)
        dimension blk_rng(2,nstmx)
c*******************************************************************
c  Input coefficient files for v-states.

        read(10,*) ninp

        print*,'  N of input files (i.e., N of virtual states) ',ninp
        print*,' '
        print*,'  Input polynomials: '
        print*,' '

        do ii = 1,ninp
          print*,'  file No. = ',ii
          idev = 100 + ii
          read(idev,*) iord(ii)
          print*,'     poly. ord = ',iord(ii)

          do kk = 1,iord(ii)+1
            read(idev,*) coef(kk,ii)
            print*,'        coef = ',coef(kk,ii)
          enddo

          read(idev,*)
          read(idev,*)

          read(idev,*) erange(1,ii),erange(2,ii)
          print*,'   v-state range: ',erange(1,ii),erange(2,ii)
        enddo

c  Set bin size.
        read(14,*) ebin

        print*,' '
        print*,'  Please check the bin size: ',ebin
        print*,' '
c*****************************
c  Get the whole parameter range.

        emax= 0.0d0
        emin= 10000.0d0

        do ii=1,ninp
          if(erange(1,ii) .lt. emin ) emin=erange(1,ii) 
          if(erange(2,ii) .gt. emax ) emax=erange(2,ii) 
        enddo

        print*,' '
        print*,' Emin & Emax for the whole = ',emin,emax
        print*,' '
c*****************************
c  Calculate point data for dln[P]/dE and output.
c    dln[P]/dE --> pgene
c    ln[P] --> pgene_i  (not used)

        print*,'  ==========  '
        print*,'    Output point data for ln[P(L)]/dE & P(L).'
        print*,' '

        iout=80

        do ii=1,ninp

          do kk=1,1000000
            ene=emin + (kk-1)*ebin

            if(ene.ge.erange(1,ii) .and. ene.le.erange(2,ii)) then
              pgene= 0.0d0
              pgene_i= 0.0d0
c  dln[P]/dE.
              do mm=1,iord(ii)+1
                pgene=pgene+coef(mm,ii)*ene**(mm-1)
              enddo
c  ln[P] (not used)
              do mm=1,iord(ii)+1
                pgene_i=pgene_i+(1.0d0/mm)*coef(mm,ii)*ene**mm
              enddo

              write(80,*) ene,pgene
cc            write(iout,*) ene,pgene,pgene_i
            endif
          enddo
        enddo
c*****************************
c Output the ranges of blocks. 

        do kk=1,ninp

          i1=kk-1
          i2=kk
          i3=kk+1

          if(i1.eq.0) then
            i1=1
            i2=2
            i3=3
            blk_rng(1,kk)=erange(1,i1)
            blk_rng(2,kk)=erange(2,i3)
          endif

          if(i1.ne.0 .and. i3.ne.ninp+1) then
            blk_rng(1,kk)=erange(1,i1)
            blk_rng(2,kk)=erange(2,i3)
          endif

          if(i3.eq.ninp+1) then
            i1=ninp-2
            i2=ninp-1
            i3=ninp
            blk_rng(1,kk)=erange(1,i1)
            blk_rng(2,kk)=erange(2,i3)
          endif

          print*,'# concatenated range N = ',i1,i2,i3

          idot=400+kk

          write(idot,201) blk_rng(1,kk),blk_rng(2,kk)
        enddo

201     format(f16.5,2x,f16.5)
c*******************************************************************
        stop
        end
