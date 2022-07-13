
        subroutine range_v

c  Get range for each v-state.
c******************************************************
        include "COMDAT"
c******************************************************
        print*,' '
        print*,'  Get v-state ranges from current denominators & output'
        print*,'  The output is used for the next procedure.'
        print*,' '
c****************
        write(16,142) de,mvert
142     format(f12.5,2x,i3,"    : bin size ")
c****************
        do mm=1,mvert
c****
          if(cpdf(1,mm).ne.0.0) then
            print*,' ???????????????????????'
            print*,' ? Make emin1 smaller. ?'
            print*,' ???????????????????????'
            stop
          endif
          if(cpdf(nwin,mm).ne.0.0) then
            print*,' ?????????????????????'
            print*,' ? Make nwim larger. ?'
            print*,' ?????????????????????'
            stop
          endif

          kstart_bin(mm)=0
          kend_bin(mm)=0

          do ii=1,nwin
            if(cpdf(ii,mm).ne.0.0) then
              kstart_bin(mm)=ii
              goto 643
            endif
          enddo

          print*,' ????????????????????????????????????'
          print*,' ? Lower boundary detection failed. ?'
          print*,' ????????????????????????????????????'
          stop

643       continue

          do ii=nwin,1,-1
            if(cpdf(ii,mm).ne.0.0) then
              kend_bin(mm)=ii
              goto 644
            endif
          enddo

          print*,' ????????????????????????????????????'
          print*,' ? Upper boundary detection failed. ?'
          print*,' ????????????????????????????????????'
          stop
644       continue

          ndifbin=kend_bin(mm) - kstart_bin(mm) + 1

          write(16,253) mm,kstart_bin(mm),kend_bin(mm),ndifbin
          write(6,253) mm,kstart_bin(mm),kend_bin(mm),ndifbin

253       format(10x,i2,2x,i6,2x,i6,2x,i4,6x,
     *         ' : v-st, bin(start), bin(end), N(bin)')
c****
c  Set the reac. coordi. range & output.

          iccc=0
          do kk=kstart_bin(mm),kend_bin(mm)
            iccc=iccc+1
            e1=emin1 + (kk-1)*de
            e2=emin1 + kk*de
            emid=emin1 + (kk-1)*de + de/2.0

            write(16,177) iccc,kk,emid,cpdf(kk,mm)
          enddo
177       format(i4,2x,i5,2x,e15.7,4x,e16.8)
c****
        enddo
c****************
c  Check: Outside the range should be zero.

        do mm=1,mvert
        do ii=1,nwin
        do jj=1,nwin
          ichkz=0
          if(cmat(jj,ii,mm).ne.0.0d0) then
            if(ii.lt.kstart_bin(mm)) ichkz=1
            if(ii.gt.kend_bin(mm)) ichkz=1

            if(jj.lt.kstart_bin(mm)) ichkz=1
            if(jj.gt.kend_bin(mm)) ichkz=1
          endif

          if(ichkz.eq.1) then
            print*,' ?????????????????????????????'
            print*,' ? Strange non-zeroposition. ?'
            print*,' ?????????????????????????????'
            print*,'  Info.: ',jj,ii,mm,cmat(jj,ii,mm)
            stop
          endif
        enddo
        enddo
        enddo
c******************************************************
        return
        end
