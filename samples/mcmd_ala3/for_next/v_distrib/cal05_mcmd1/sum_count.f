
        subroutine sum_count(kvert,ivlocal,kp_vst,ibin)
c******************************************************
        include "COMDAT"
c******************************************************
        ican=0
c********
c  First snapshot: Just memorize.
        if(kkone_use.eq.1) then
          ivlocal=1
          kp_vst=kvert
          nwork(1) = ibin
c_____________
cc        write(6,155) kkone_use,kvert
155       format("     v-change (1): ",i8,2x,i2)
c_____________
          goto 851
        endif
c****
c  When v-state changed, reset quantities.
        if(kvert.ne.kp_vst) then
          iold=ivlocal
          ivlocal=1
          kp_vst=kvert
          nwork(1) = ibin
c_____________
cc        write(6,156) kkone_use,kvert,iold
156       format("     v-change (2): ",i8,2x,i2,2x,";old ivlocal=",i6)
c_____________
          goto 851
        endif
c****
c  Add information.
        ivlocal=ivlocal+1

        if(ivlocal.le.intvc) then
          nwork(ivlocal) = ibin
          goto 851
        endif
c**
c  Counting.
        if(ivlocal.eq.intvc+1) then
          nwork(intvc+1) = ibin
          ican=1
        endif
        if(ivlocal.gt.intvc+1) then
          do kk=2,intvc+1
            nwork(kk-1)=nwork(kk)
          enddo
          nwork(intvc+1) = ibin
          ican=1
        endif

c  Accumulate count matrix.
        if(ican.eq.1) then
          ist=nwork(1)
          ien=nwork(intvc+1)
          ncmat(ien,ist,kp_vst)=ncmat(ien,ist,kp_vst)+1
        endif
c********
851     continue
c******************************************************
        return
        end
