
c**************************************
        read(11,*) kkk
c****
        num=0
        num1=0

400     continue
        read(10,*,end=999) istp,iv
cc        write(6,300) istp,iv
300     format('  step No. & v state = ',i10,i2)

        num=num+1
        if(iv.eq.1) num1=num1+1

        goto 400
c****
999     continue
        n1=10000
        n2=990000000

        read(5,*) na
        read(5,*) nb
        read(5,*) nc
        read(5,*) nd
cc        na=64
cc        nb=16
cc        nc=4 
cc        nd=1 

        aa = real(num1)/real(num)
        if(aa.gt.0.95) then
          write(20,201) kkk,n1,n2,na,aa
          write(6,101) kkk,aa
101       format(i4,'  ',f10.3,'  very deep trap ')
        else
          if(aa.gt.0.90) then
            write(20,201) kkk,n1,n2,nb,aa
            write(6,102) kkk,aa
102         format(i4,'  ',f10.3,'  deep trap ')
          else
            if(aa.gt.0.80) then
              write(20,201) kkk,n1,n2,nc,aa
            else
              write(20,201) kkk,n1,n2,nd,aa
              write(6,103) kkk,aa
103           format(i4,'  ',f10.3,'  ')
            endif
          endif
        endif

201     format(i4,2x,i8,2x,i9,2x,i4,2x,";",f9.4)
c**************************************
        stop
        end
