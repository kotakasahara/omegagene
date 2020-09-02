
        subroutine t_mat

c  Calc. state vector & t. mat. from c. mat.

c******************************************************
        include "COMDAT"
c******************************************************
        print*,' '
        print*,' ********************************************'
        print*,' * Calculate t. mat. from c. mat. & output. *'
        print*,' ********************************************'
        print*,' '
c********
c  Calculate state vector for each v-state.
c    Note that ncmat involves corrections done in check_cmat.
c    Thus, npdf involves corrections.

        do mm=1,mvert
          do ii=1,nwin
            npdf(ii,mm)=0

            do jj=1,nwin
              npdf(ii,mm)=npdf(ii,mm)+ncmat(jj,ii,mm)
            enddo
          enddo
        enddo

c  Check.
        do mm=1,mvert
          do ii=1,nwin

            if(npdf(ii,mm).lt.0) then
              print*,' ???????????????????????????????'
              print*,' ? Negative prob. is detected. ?'
              print*,' ???????????????????????????????'
              stop
            endif

            do jj=1,nwin
              if(ncmat(jj,ii,mm).lt.0) then
                print*,'?????????????????????????????'
                print*,' ? Negative count detected. ?'
                print*,'?????????????????????????????'
                stop
              endif
            enddo
          enddo
        enddo
        
c********
        write(16,142) de
142     format(f12.5,"    : bin size ")
c********
c  Get range for each v-state.
        do mm=1,mvert
c****
          if(npdf(1,mm).ne.0) then
            print*,' ???????????????????????'
            print*,' ? Make emin1 smaller. ?'
            print*,' ???????????????????????'
            stop
          endif
          if(npdf(nwin,mm).ne.0) then
            print*,' ?????????????????????'
            print*,' ? Make nwim larger. ?'
            print*,' ?????????????????????'
            stop
          endif

          kstart_bin(mm)=0
          kend_bin(mm)=0

          do ii=2,nwin
            if(npdf(ii,mm).ne.0) then
              kstart_bin(mm)=ii
              goto 643
            endif
          enddo
643       continue

          do ii=nwin-1,1,-1
            if(npdf(ii,mm).ne.0) then
              kend_bin(mm)=ii
              goto 644
            endif
          enddo

          print*,' ??????????????????????????????????'
          print*,' ? Strange. range is not defined. ?'
          print*,' ??????????????????????????????????'
          stop
644       continue

          ndifbin=kend_bin(mm) - kstart_bin(mm) + 1

          write(16,253) mm,kstart_bin(mm),kend_bin(mm),ndifbin

          write(6,253) mm,kstart_bin(mm),kend_bin(mm),ndifbin

253       format(10x,i2,2x,i6,2x,i6,2x,i4,6x,
     *         ' : v-st, bin(start), bin(end), N(bin)')
c****
c  Set the reac. coordi. range.

          iccc=0
          do kk=kstart_bin(mm),kend_bin(mm)
            iccc=iccc+1
            e1=emin1 + (kk-1)*de
            e2=emin1 + kk*de
            emid=emin1 + (kk-1)*de + de/2.0

            write(16,177) iccc,kk,emid
          enddo
177       format(i4,2x,i5,2x,e15.7)
c****
        enddo
c********
c  Check.
        do mm=1,mvert
          do ii=1,nwin
          do jj=1,nwin
            if(npdf(ii,mm).eq.0) then
              if(ncmat(jj,ii,mm).ne.0) then
                print*,' ??????????????????????????????????'
                print*,' ? Strange. ncmat should be zero. ?'
                print*,' ??????????????????????????????????'
                stop
              endif
            endif
          enddo
          enddo
        enddo
c********
c  Calc. t. mat.

        do mm=1,mvert
          k1=kstart_bin(mm)
          k2=kend_bin(mm)

cc        do ii=1,nwin
cc        do jj=1,nwin
          do ii=k1,k2
          do jj=k1,k2
            if(npdf(ii,mm).ne.0) then
              tmat(jj,ii,mm) =
     *          real(ncmat(jj,ii,mm)) / real(npdf(ii,mm))
            endif

c____
cc          if(mm.eq.1 .and. jj.eq.352 .and. ii.eq.354) then
cc            print*,'  debug: '
cc            print*,'     ',jj,ii,tmat(jj,ii,mm)
cc          endif
c____
          enddo
          enddo
        enddo

c  Output t. mat.

        do mm=1,mvert
          k1=kstart_bin(mm)
          k2=kend_bin(mm)
          nbin=k2-k1+1

          write(21,130) mm,k1,k2,nbin
130       format(i4,2x,i4,2x,i4,2x,i4,"      : info. ")

          do ii=k1,k2
          do jj=k1,k2
            write(21,129) jj,ii,tmat(jj,ii,mm)
          enddo
          enddo
        enddo
129     format(i4,2x,i4,2x,e22.12)
c******************************************************
        return
        end
