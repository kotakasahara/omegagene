
        subroutine out_v_info
c******************************************************
        include "COMDAT"
c******************************************************
c  Output the statisctics to the monitor file.

        print*,' '
        print*,'  N of input data used = ',kktot

        do kk=1,mvert
          write(6,251) kk,icou(kk)
        enddo
251     format('      Ndetect (vert) = ',i2,2x,i10)
c******************************************************
        return
        end
