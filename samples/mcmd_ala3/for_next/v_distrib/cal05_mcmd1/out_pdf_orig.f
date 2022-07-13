
        subroutine out_pdf_orig
c******************************************************
        include "COMDAT"
c******************************************************
c  Output the original pdf.

        do mm=1,mvert
          if(mm.eq.1) idevout=101
          if(mm.eq.2) idevout=102
          if(mm.eq.3) idevout=103
          if(mm.eq.4) idevout=104
          if(mm.eq.5) idevout=105
          if(mm.eq.6) idevout=106
          if(mm.eq.7) idevout=107
          if(mm.eq.8) idevout=108
          if(mm.eq.9) idevout=109
          if(mm.eq.10) idevout=110
          if(mm.eq.11) idevout=111
          if(mm.eq.12) idevout=112
          if(mm.eq.13) idevout=113
          if(mm.eq.14) idevout=114
          if(mm.eq.15) idevout=115

c  pdf is transformed to log scale.

          do ii=1,nene
            e1=emin1+(ii-1)*de+de/2.0
cc          aa=pdf(ii,mm)/tot
cc          if(aa.ne.0.0) then
cc            pdf(ii,mm)=log(aa)
cc          endif

            write(idevout,612) e1,pdf(ii,mm)
          enddo
        enddo

612     format(e15.7,2x,e15.7)
c******************************************************
        return
        end
