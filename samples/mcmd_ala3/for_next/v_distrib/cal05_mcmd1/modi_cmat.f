
        subroutine modi_cmat
c  Symmetrize c. mat. 
c******************************************************
        include "COMDAT"
c******************************************************
        print*,' 1) Simply take average of count matrix as:'
        print*,'    (c(i,j) + c(j,i))/2 --> c(i,j) = c(j,i)'
        print*,'  '
c********
c  Simply take average.

        do kk=1,mvert
        do jj=1,nwin
        do ii=1,nwin
          cmat(ii,jj,kk) = ncmat(ii,jj,kk) + ncmat(jj,ii,kk)
          cmat(ii,jj,kk) = 0.5d0 * cmat(ii,jj,kk)
        enddo
        enddo
        enddo
c****
c  Check.

        print*,' Check start. '
        do kk=1,mvert
        do jj=1,nwin
        do ii=1,nwin
          if(cmat(ii,jj,kk).lt.0.0) then
            print*,' ????????????????????????????????????'
            print*,' ? Negative c mat ele was detected. ?'
            print*,' ????????????????????????????????????'
            print*,' ii,jj,mm,cmat(ii,jj,kk)=',ii,jj,kk,cmat(ii,jj,kk)

            print*,' ii,jj,kk,ncmat(ii,jj,kk)=',ii,jj,kk,ncmat(ii,jj,kk)
            print*,' jj,ii,kk,ncmat(jj,ii,kk)=',jj,ii,kk,ncmat(jj,ii,kk)

            stop
          endif

          if(cmat(ii,jj,kk).ne.cmat(jj,ii,kk)) then
            print*,' ????????????????????????????????????'
            print*,' ? Detected: c. mat. not symmetric. ?'
            print*,' ????????????????????????????????????'
            stop
          endif
        enddo
        enddo
        enddo
        print*,' Check end. '
c******************************************************
        return
        end
