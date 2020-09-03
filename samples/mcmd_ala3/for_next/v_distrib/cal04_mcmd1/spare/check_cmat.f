
        subroutine check_cmat
c  Check the c. mat & modify slightly.
c******************************************************
        include "COMDAT"
c******************************************************
        icd=0

        do kk=1,mvert
        do jj=1,nwin
        do ii=1,nwin
          if(ncmat(ii,jj,kk) .eq. 0 .and. ncmat(jj,ii,kk) .ne. 0 ) then
            icd=icd+1
cc          write(6,100) ii,jj,ncmat(ii,jj,kk),ncmat(jj,ii,kk)
            ncmat(ii,jj,kk)=1
          endif
          if(ncmat(ii,jj,kk) .ne. 0 .and. ncmat(jj,ii,kk) .eq. 0 ) then
            icd=icd+1
cc          write(6,100) ii,jj,ncmat(ii,jj,kk),ncmat(jj,ii,kk)
            ncmat(jj,ii,kk)=1
          endif
        enddo
        enddo
        enddo
100     format("     ncmat (ii,jj,kk) & (jj,ii,kk):",
     *          i4,2x,i4,4x,i6,2x,i6)

        print*,' '
        print*,'  N. of imvalance data = ',icd
        print*,' '
c********
c  Check again.
        icd=0

        do kk=1,mvert
        do jj=1,nwin
        do ii=1,nwin
          if(ncmat(ii,jj,kk) .eq. 0 .and. ncmat(jj,ii,kk) .ne. 0 ) then
            icd=icd+1
            ncmat(ii,jj,kk)=1
          endif
          if(ncmat(ii,jj,kk) .ne. 0 .and. ncmat(jj,ii,kk) .eq. 0 ) then
            icd=icd+1
            ncmat(jj,ii,kk)=1
          endif
        enddo
        enddo
        enddo

        if(icd.ne.0) then
          print*,' ?????????????????????????????????????????????'
          print*,' ? Something is strange. Stop in check_cmat. ?'
          print*,' ?????????????????????????????????????????????'
          stop
        endif
c******************************************************
        return
        end
