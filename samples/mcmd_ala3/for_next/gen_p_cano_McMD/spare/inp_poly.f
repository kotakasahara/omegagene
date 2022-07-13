c
c  This pgm is avairable for V-McMD.
c
        implicit real*8 (a-h,o-z)
c
        dimension iord(50)
        dimension coef(20,50)
        dimension erange(2,50)
c*******************************************************************
c  Input bin size.
        read(11,*) ebin
        print*,'  '
        print*,'  Bin size for output data = ',ebin
c**********
c  Input polynomial coefficient files.

        read(10,*) ninp

        print*,'  N of input files (i.e., N of virtual states) ',ninp
        print*,' '
        print*,'  Input polynomials below: '
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
          print*,'      range: ',erange(1,ii),erange(2,ii)
        enddo
c*****************************
c  Get the whole parameter range.

        emax= -1.0e+10
        emin=  1.0e+10

        do ii=1,ninp
          if(erange(1,ii) .lt. emin ) emin=erange(1,ii) 
          if(erange(2,ii) .gt. emax ) emax=erange(2,ii) 
        enddo

        print*,' '
        print*,' Emin & Emax for the whole = ',emin,emax
        print*,' '
c*****************************
c  Calculate dln[P]/dE & ln[P] from polynomials.
c    dln[P]/dE --> pgene
c    ln[P] --> pgene_i
c
c   The function ln[P] involves discontinuity usually.

        print*,'  **********'
        print*,'    Output point data for ln[P(L)]/dE & P(L).'
        print*,' '

        do ii=1,ninp

          iout=1000 + ii

          do kk=1,1000000
            ene=emin + (kk-1)*ebin

            if(ene.ge.erange(1,ii) .and. ene.le.erange(2,ii)) then
              pgene= 0.0d0
              pgene_i= 0.0d0

              do mm=1,iord(ii)+1
                pgene=pgene+coef(mm,ii)*ene**(mm-1)
              enddo

              do mm=1,iord(ii)+1
                pgene_i=pgene_i+(1.0d0/mm)*coef(mm,ii)*ene**mm
              enddo

              write(iout,*) ene,pgene,pgene_i
              print*,'  E & func = ',ene,pgene,pgene_i
            endif
          enddo
        enddo
c****************
c  At the end, generate an input file for the next task.

        eset=emin-1000.0
cc      eset=emin+1000.0
        adum=999.0
        bdum=0.0
        k_order=11

        write(30,300) eset,ebin,adum,bdum,k_order,bdum
300     format(f10.1,2x,f8.3,2x,f8.0,2x,f6.1,2x,i3,2x,f6.1)

        write(30,305) emin,emax,emin,emax
305     format(4(f10.1,2x))
        write(30,302) emin,emax
302     format(2(f10.1,2x))
c*******************************************************************
        stop
        end
