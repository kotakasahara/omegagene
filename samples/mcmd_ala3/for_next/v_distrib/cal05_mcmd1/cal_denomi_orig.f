
        subroutine cal_denomi_orig

c  Calculate denominator for the current c. mat.
c******************************************************
        include "COMDAT"
c******************************************************
        do mm=1,mvert
          do ii=1,nwin
          cpdf_orig(ii,mm)=0.0

          do jj=1,nwin
            cpdf_orig(ii,mm)=cpdf_orig(ii,mm)+cmat_orig(jj,ii,mm)
          enddo
        enddo
        enddo
c******************************************************
        return
        end
