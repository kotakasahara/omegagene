
        subroutine t_mat

c  Calc. state vector & t. mat. from c. mat.
c******************************************************
        include "COMDAT"
c******************************************************
        print*,' '
        print*,' ********************************************'
        print*,' * Calculate t. mat. from c. mat. & output. *'
        print*,' ********************************************'
        print*,' '
c****************
c  1) Modify c. mat. by taking simple average.
c  This modified c. mat. may be used as the initial c. mat.
c  for the maximum entropy procedure.
        call modi_cmat
c****
c  2) Do maximum entropy procedure.
c  If you do not use the maximum entropy procedure,
c  then mask the below call.
        call maximum_entropy
c****************
c  Calculate denominator using the current c. mat.
        call cal_denomi
c  Get range for each v-state & output.
        call range_v
c****************
c  c. mat. --> t. mat.

        do mm=1,mvert
          k1=kstart_bin(mm)
          k2=kend_bin(mm)

          do ii=k1,k2
            do jj=k1,k2
              tmat(jj,ii,mm) = cmat(jj,ii,mm) / cpdf(ii,mm)
            enddo
          enddo
        enddo

c  Check invariant.
        do mm=1,mvert
          k1=kstart_bin(mm)
          k2=kend_bin(mm)

          do ii=k1,k2
            asum=0.0d0

            do jj=k1,k2
              asum=asum+tmat(jj,ii,mm)
            enddo

            dsum=1.0d0-asum
            eps=1.0e-15
            if(dsum.gt.eps) then
              print*,' ????????????????????????'
              print*,' ? Invariant is failed. ?'
              print*,' ????????????????????????'
              print*,'  D(sum) = ',ii,dsum
              stop
            endif
          enddo
        enddo
c********
c  Output t. mat.

        write(21,163) mvert
163     format(i4)

        do mm=1,mvert
          k1=kstart_bin(mm)
          k2=kend_bin(mm)
          nbin=k2-k1+1

          write(21,130) mm,k1,k2,nbin
130       format(i4,2x,i4,2x,i4,2x,i4,"      : info. ")

          do ii=k1,k2
          do jj=k1,k2
            write(21,129) jj,ii,tmat(jj,ii,mm)
          enddo
          enddo
        enddo
129     format(i4,2x,i4,2x,e23.16)
c******************************************************
        return
        end
