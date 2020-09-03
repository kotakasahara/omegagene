c
c  This program supposes that the simulation
c  temperature is 300 K.
c********************************
       implicit real*8 (a-h,o-z)
c********************************
        do ii=1,150
cc          c1=0.629d0 + (ii-1)*0.005d0
          c1 = (ii-1)*0.01d0 + 0.5d0

          g = 503.217d0
          temp = g / c1
          write(6,100) ii,c1,temp
100       format(" c1 & T = ",i3,2x,f10.4,2x,f10.4)
        enddo
c********************************
        stop
        end
