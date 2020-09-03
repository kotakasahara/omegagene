
        subroutine kp_cmat_orig

c  Keep the original c. mat.
c******************************************************
        include "COMDAT"
c******************************************************
        do kk=1,mvert
        do jj=1,nwin
        do ii=1,nwin
          cmat_orig(ii,jj,kk) = ncmat(ii,jj,kk)
        enddo
        enddo
        enddo
c******************************************************
        return
        end
