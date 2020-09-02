
        subroutine cal_pdf
c******************************************************
        include "COMDAT"
c******************************************************
c  Determine pdf.

        do mm=1,mvert
        do ii=1,nene
cc        e1=emin1+(ii-1)*de+de/2.0
          aa=pdf(ii,mm)/tot
          if(aa.ne.0.0) then
            pdf(ii,mm)=log(aa)
          endif
        enddo
        enddo
c******************************************************
        return
        end
