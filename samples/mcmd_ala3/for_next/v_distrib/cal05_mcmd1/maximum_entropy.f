
        subroutine maximum_entropy

c  Do maximum entropy procedure for c. mat.
c  The entropy procedure method I use is slightly different from
c  the Noe or Pance methods.
c  
c
c
c******************************************************
        include "COMDAT"
c******************************************************
        print*,' '
        print*,' * Do maximum entropy procedure for c. mat. *'

        read(27,*) eps
        print*,'     convergence factor = ',eps
        print*,' '
c****************
c  Preparation.

c  Keep the original c. mat., which may be asymmetric.
        call kp_cmat_orig
c  Calculate denominator using the original c. mat.
        call cal_denomi_orig
c****************
        nloop=0
887     continue
        nloop=nloop+1
c****
c  Calculate the denominator for the current c. mat.
c  The current c. mat. may be from subroutine modi_cmat or may not.
        call cal_denomi
c****
c  Keep c. mat. for convergence check.
        do mm=1,mvert
          do ii=1,nwin
          do jj=1,nwin
            cmat_chk(jj,ii,mm)=cmat(jj,ii,mm)
          enddo
          enddo
        enddo
c****
c  Upate c. mat.

c  Diagonal.
        do mm=1,mvert
          do ii=1,nwin
            if(cpdf(ii,mm).ne.0.0d0) then
              cmat(ii,ii,mm)=cmat_orig(ii,ii,mm)*
     *                          cpdf_orig(ii,mm)/cpdf(ii,mm)
            endif
          enddo
        enddo

c  Off-diagonal.
        do mm=1,mvert
          do ii=1,nwin
            if(cpdf(ii,mm).eq.0.0d0) goto 881

            do jj=1,nwin
              if(cpdf(jj,mm).eq.0.0d0) goto 882
              p1=cmat(jj,ii,mm)
              p2=cmat(ii,jj,mm)

              if(p1.eq.0.0d0 .and. p2.eq.0.0d0) goto 883

              qq=cmat_orig(jj,ii,mm)+cmat_orig(ii,jj,mm)
              rr= cpdf_orig(ii,mm)/cpdf(ii,mm) +
     *                  cpdf_orig(jj,mm)/cpdf(jj,mm)
              rr=1.0d0/rr
              cmat(jj,ii,mm)=qq*rr
              cmat(ii,jj,mm)=qq*rr

883           continue
882           continue
            enddo

881         continue
          enddo
        enddo

c  Convergence check.

        dsum=0.0d0
        esum=0.0d0
        iccc=0
        do mm=1,mvert
          do ii=1,nwin
          do jj=1,nwin
            if(cmat_chk(jj,ii,mm).ne.0.0d0) then
              iccc=iccc+1
              ddd=cmat_chk(jj,ii,mm)-cmat(jj,ii,mm)
              dsum=dsum+ddd**2
              esum=esum+cmat_chk(jj,ii,mm)**2
            endif
          enddo
          enddo
        enddo

        fff=sqrt(dsum)/sqrt(esum)

        print*,'  count = ',nloop
        print*,'    change of c. mat. = ',fff

        if(fff.gt.eps) goto 887
cc      if(fff.gt.1.0e-7) goto 887
c******************************************************
        return
        end
