      program fitting2noene
c========*=========*=========*=========*=========*=========*=========*==

      implicit		double precision (a-h,o-z), integer (i-n) 
      integer		maxene
      parameter		(maxene=10000)
      integer		i, j, ie, iemin, iemax, iee,
     &			numenef, numene, numenew, numenewle, numenewde,
     &			iprob(maxene)
      double precision	eex(maxene), probly(maxene),
     &			ee(maxene), probl(maxene),
     &			eef(maxene), problf(maxene),
     &			eefle(maxene), problfle(maxene),
     &			eedew(maxene), dew(maxene),
     &			emin, bin, temp, elw, eup,
     &			eetmp, felw, feup, welw, weup, celw, ceup,
     &			dprobl,
     &			loop,
     &			alphalw, alphaup, betalw, betaup

      integer		maxdeg, mdeg, ndeg, ierr
      parameter		(maxdeg=20)
      double precision	eps, c,
ccccc
ccccc     &			w(maxene), r(maxdeg+1),
ccccc
     &			w(maxene), r(maxene),
     &			a(3*(maxene+maxdeg+1)),
     &			tc(maxdeg+1)

      character work*3

      read(5,*) emin, bin, temp, eps, mdeg, c
      read(5,*) felw, feup, welw, weup
      read(5,*) celw, ceup

c  Modifyed by Higo.
883   continue
      read(5,332) work
332   format(a3)
      if(work.eq.'###') then
        goto 883
      else
        backspace(5)
        goto 885
      endif
885   continue
c  Modification end.

      iee = 0
      do 100 i = 1, maxene
        read(5,*,end=200) ee(i), probl(i)
        iee = iee + 1
100   continue
200   continue
      numene = iee

      iee = 0
      do 300 i = 1, numene
        if ( (ee(i) .ge. felw) .and. (ee(i) .le. feup) ) then
          iee = iee + 1
          eex(iee)    = ee(i)
          probly(iee) = probl(i)
        else if ( ee(i) .gt. feup .or. ee(i).lt. felw) then

          print*,' Skip: energy out of range: ',i,ee(i)
cc          goto 400

        endif
300   continue
400   continue
      numenef = iee
ccccc
      w(1)=-1.0D0
ccccc

      rt = 1.98717D-3*temp

c========*=========*=========*=========*=========*=========*=========*==
      call dpolft(numenef, eex, probly, w, mdeg, ndeg,
     &                                   eps, r, ierr, a)
      call dpcoef(ndeg, c, tc, a)
      if (ierr .ne. 1) then
        write(6, '("ierr = ", i2)') ierr
      else
cc        write(20, '("ierr = ", i2)') ierr
cc        write(20, '("ndeg = ", i2)') ndeg
cc        write(20, '("eps  = ", f12.7)') eps
cc        write(20, '("c    = ", e15.7)') c
cc        write(20, '(e14.7)') (tc(i), i=1,ndeg+1)
cc        write(20, '()')

        write(20, '(i2)') ndeg

        write(20, '(e22.15)') (tc(i), i=1,ndeg+1)
cc        write(20, '()')
      endif       

c========*=========*=========*=========*=========*=========*=========*==
      iee = 0
      do 4000 loop = welw, weup, bin
        iee = iee + 1
        eef(iee) = dble(loop)
        problf(iee) = tc(1)
        do 3000 j = 2, ndeg+1
          problf(iee) = problf(iee) + tc(j)*((eef(iee)-c)**(j-1))
3000    continue
4000  continue
      numenew = iee

c========*=========*=========*=========*=========*=========*=========*==
      celw = dnint(celw/bin)*bin
      ceup = dnint(ceup/bin)*bin

      alphalw = tc(2)
      alphaup = tc(2)
      do 5000 i = 2, ndeg
        alphalw = alphalw + dble(i)*tc(i+1)*((celw-c)**(i-1))
        alphaup = alphaup + dble(i)*tc(i+1)*((ceup-c)**(i-1))
5000  continue
      betalw  = tc(1)
      betaup  = tc(1)
      do 5500 i = 2, ndeg+1
        betalw  = betalw  + tc(i)*((celw-c)**(i-1))
        betaup  = betaup  + tc(i)*((ceup-c)**(i-1))
5500  continue
cc      write(20, '(e14.7)') alphalw
cc      write(20, '(e14.7)') alphaup
cc      write(20, '()')
      write(20, '(e22.15)') alphalw
      write(20, '(e22.15)') alphaup
cc      write(20, '()')
cc      write(20, '(e14.7)') betalw
cc      write(20, '(e14.7)') betaup
cc      write(20, '()')
cc      write(20, '(e22.15)') betalw
cc      write(20, '(e22.15)') betaup

      write(20,'(f9.1,x,f9.1)') celw, ceup
      write(20,'(f6.1)') temp

c========*=========*=========*=========*=========*=========*=========*==
      iee = 0
      do 6500 loop = welw, weup, bin
        iee = iee + 1
        eefle(iee) = dble(loop)
        if (eefle(iee) .le. celw) then
          problfle(iee) = alphalw*(eefle(iee) - celw) + betalw
        else if (eefle(iee) .ge. ceup) then
          problfle(iee) = alphaup*(eefle(iee) - ceup) + betaup
        else
          problfle(iee) = tc(1)
          do 6000 j = 2, ndeg+1
            problfle(iee) = problfle(iee)
     &                    + tc(j)*((eefle(iee)-c)**(j-1))
6000      continue
        endif
6500  continue
      numenewle = iee

c========*=========*=========*=========*=========*=========*=========*==
      iee = 0
      do 8500 loop = welw, weup, bin
        iee = iee + 1
        eedew(iee) = dble(loop)
        if (eedew(iee) .le. celw) then
          dprobl = alphalw
        else if (eedew(iee) .ge. ceup) then
          dprobl = alphaup
        else
          dprobl = dble(1)*tc(1+1)
          do 8000 j = 2, ndeg
            dprobl = dprobl + dble(j)*tc(j+1)*((eedew(iee)-c)**(j-1))
8000      continue
        endif
        dew(iee) = 1.0D0 + rt*dprobl
8500  continue
      numenewde = iee

c========*=========*=========*=========*=========*=========*=========*==
      write(11, '(e15.7, e15.7)') (ee(i), probl(i), i=1,numene)
      write(12, '(e15.7, e15.7)')
     &                  (eef(i), problf(i), i=1,numenew)
      write(13, '(e15.7, e15.7)')
     &                  (eefle(i), problfle(i), i=1,numenewle)

      write(16, '(e15.7, e15.7)')
     &                  (eedew(i), dew(i), i=1,numenewde)

      stop
      end
