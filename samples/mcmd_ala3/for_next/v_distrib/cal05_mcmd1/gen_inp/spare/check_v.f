
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

        aa = real(num1)/real(num)
        if(aa.gt.0.99) then

          n3=16

          write(20,201) kkk,n1,n2,n3
          write(6,101) aa
101       format('  ',f10.3,'  very deep trap ')
        else
          if(aa.gt.0.9) then

            n3=4

            write(20,201) kkk,n1,n2,n3
            write(6,102) aa
102         format('  ',f10.3,'  deep trap ')
          else

            n3=1

            write(20,201) kkk,n1,n2,n3
          endif
        endif

201     format(i4,2x,i8,2x,i9,2x,i4)
c**************************************
        stop
        end
