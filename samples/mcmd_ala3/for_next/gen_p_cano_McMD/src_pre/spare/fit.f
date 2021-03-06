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

      read(5,*) emin, bin, temp, eps, mdeg, c
      read(5,*) felw, feup, welw, weup
      read(5,*) celw, ceup

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
ccccc
ccccc          w(iee)      = 1.0D0
ccccc          w(iee)      = 1.0D0/(probly(iee)**2)
ccccc
        else if ( ee(i) .gt. feup ) then
          goto 400
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
        write(20, '("ierr = ", i2)') ierr
        write(20, '("ndeg = ", i2)') ndeg
        write(20, '("eps  = ", f12.7)') eps
        write(20, '("c    = ", e15.7)') c
        write(20, '(e14.7)') (tc(i), i=1,ndeg+1)
        write(20, '()')
        write(20, '(e22.15)') (tc(i), i=1,ndeg+1)
        write(20, '()')
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
      write(20, '(e14.7)') alphalw
      write(20, '(e14.7)') alphaup
      write(20, '()')
      write(20, '(e22.15)') alphalw
      write(20, '(e22.15)') alphaup
      write(20, '()')
      write(20, '(e14.7)') betalw
      write(20, '(e14.7)') betaup
      write(20, '()')
      write(20, '(e22.15)') betalw
      write(20, '(e22.15)') betaup

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
