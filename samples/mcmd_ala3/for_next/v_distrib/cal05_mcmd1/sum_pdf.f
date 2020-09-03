        subroutine sum_pdf(mm,a1,kvert,ibin)
c******************************************************
        include "COMDAT"
c******************************************************
c  Sum up for pdf.
c    The index ibin is universal for all v-pdf.

        e1=a1-emin1
        ibin=int(e1/de)+1
        if(a1.lt.emin1) goto 750
        if(ibin.gt.nene) goto 750

        kktot=kktot+1
        tot=tot+wgt(mm)
        kkone_use=kkone_use+1

        pdf(ibin,kvert)=pdf(ibin,kvert) + wgt(mm)

750     continue
c******************************************************
        return
        end
