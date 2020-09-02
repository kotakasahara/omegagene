c******************************************************
        include "COMDAT"
c******************************************************
        print*,' '
        print*,' ********************************************'
        print*,' * Energy distrib. for McMD-vertial system. *'
        print*,' ********************************************'
        print*,'    If pdf(E) = zero at a point of E, set ln[pdf] = 0.'
        print*,' '
        print*,'  NOTE 1: [Important] N. of vertial states = ',nvert
        print*,'          This value is hard coded in this pgm.'
        print*,'  NOTE 2:  Normalization is dome in the whole range'
        print*,'         not in each range of vertial quantity.'
        print*,' '
        print*,'  NOTE 3: This pgm resets the level of the dustrib. '
        print*,' '
c********************************************
c  Parameter input field.

        call para_inp
c********************************************
c  Loop over input files.

cc      do mm=1,1
        do mm=1,nfil
c************
c  (f1) Input a file for v-state profile.

          ipp=0
          idevinv=1000+mm
          open(unit=idevinv,file=filnamev(mm),status="old")

          print*,' '
          write(6,428) mm,filname(mm),filnamev(mm)
428       format(i4,2x,"Input files = ",a15,2x,a15)
          write(6,277) nst(mm),nen(mm)
277       format('          Used trajectory interval = ',i9,2x,i9)

844       continue
          read(idevinv,*,end=996) iii,jjj
          ipp=ipp+1
          ipos(ipp)=iii
          ivert(ipp)=jjj

c  Error check.
          if(jjj.gt.mvert) call err_msg(1)
          goto 844

          close(idevinv)

996       continue

          print*,'     N of data for vertial system = ',ipp
          
          if(ipp.eq.1) then
            print*,' '
            print*,'   NOTE: A complete intaval is not input. '
            print*,'         Convert: Last step --> last of interval.'
            print*,' '
          endif

c  Error check.
          if(ipp.gt.nvchg) call err_msg(2)
c********
c  (f2) Input a file for the react. coord. profile.
c       kkone_use: N of used snapshots in the current input file.

          idevin=1000+mm
          kkone=0
          kkone_use=0
          open(unit=idevin,file=filname(mm),status="old")

888       continue
c  Input a snapshot.
          read(idevin,101,end=900) a1
101       format(e15.7)

c  Count snapshots in the file whether used or not.
          kkone=kkone+1

c  Pick up the data in an allowed range.
          if(kkone.gt.nst(mm) .and. kkone.le.nen(mm)) then

c  Assign the vertial No. to the snapshot input.
            if(ipp.ne.1) then
              do kk=1,ipp-1
                if(kkone.ge.ipos(kk) .and. kkone.lt.ipos(kk+1)) then
                  kvert=ivert(kk)
                  goto 444
                endif
              enddo
            endif

c  For the trajectory tail.
c    Note that this part works only when the simulation was incomplete.
c    Get the interval for virtual states.
            mmm=ipos(ipp)-ipos(ipp-1)

            if(kkone.ge.ipos(ipp)) then
              if(kkone.lt.ipos(ipp) + mmm) then
                kvert=ivert(ipp)
                goto 444
              else
c    Stop reading the current file.
                goto 900
              endif
            endif

c  Exit when the v-state No is not assigned. 
            call err_msg(3)

444         continue
c****
c  Up to here, v-state No. was assigned to the snapshot.
            icou(kvert)=icou(kvert)+1
c****
c  Sum up for pdf.
c    The index ibin is universal for all v-pdf.

            call sum_pdf(mm,a1,kvert,ibin)

c****
c  Procedure for c. mat. calc.
c  Here, the v-state No. & the bin posotion has been determined.

            if(kkone_use.eq.0) call err_msg(4)

            call sum_count(kvert,ivlocal,kp_vst,ibin)
c****
          endif

c  Go to the next in the input file.
          goto 888

900       continue

          close(idevin)
c************
        enddo

c  End of input files.
c****************
c  Output the statisctics to the monitor file.
        call out_v_info
c********************************************
c  Determine pdf.
        call cal_pdf
c********
c  Output the original pdf.
        call out_pdf_orig

c  Shift the pdf & output.

        call out_pdf_shft

c****
c  Check the c. mat. and modify it slightly.

        call check_cmat

c  Calc. state vector & t. mat.

        call t_mat

        print*,' '
        print*,'  End of all. '
        print*,' '
c****
c________
cc      isel=5
cc      do ii=1,nwin
cc        if(npdf(ii,isel).ne.0) then
cc          print*,'  npdf = ',ii,npdf(ii,isel)
cc        endif
cc      enddo
c******************************************************
        stop
        end
