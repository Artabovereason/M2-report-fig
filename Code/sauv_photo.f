      implicit real*8(a-h,o-z)
      common/coul/zc,pf
      common/hydr/nh,lh
!
      eau=27.2116d0
      a0=0.529d-10
      cc=137.03604d0
      pi=4.d0*datan(1.d0)
      a02 = 0.280028D0
!
      open(unit=15,file='energy.txt')
      read(15,*),ei
      close(15)
c      ei=2.102d0!5.102d0
!      ei=5.102d0
!
      ei=ei/eau
      nh=3
      lh=0
      eh=-1.d0/2.d0/nh/nh
      zc=1.d0
*
* ei incoming energy of the photon (eV)
*
      ef=ei+eh
      pf=dsqrt(2.d0*ef)
      etaf=-zc/pf
*
      if(lh.eq.0)then
      l=lh+1
      call rad(l,xr,xn)
      write(6,*)'l xn xr=',l,xn,xr
      x1=dsqrt((2*l+1.d0))*d3j30(l,1,lh)
      x2=(x1*xr)**2.d0
      xt=x2
      endif
!
      if(lh.ne.0)then
      xt=0.d0
      do 1 l=lh-1,lh+1,2
      call rad(l,xr,xn)
      write(6,*)'l xn xr=',l,xn,xr
      x1=dsqrt((2*l+1.d0))*d3j30(l,1,lh)
      x2=(x1*xr)**2.d0
      xt=xt+x2
1     continue
      endif
c
      xsec=8.d0*pi*ei/3.d0/cc/pf
      xsec=xsec*xt
!
      xsec_Kr=64.d0*pi/3.d0/dsqrt(3.d0)/cc/8.d0/ei**3.d0/nh**5.d0
c
c    xsec_Kr : Kramer approximation to the cross section, sigma propto n^{-5}
c
!
c      write(6,*)'xsec(pi a_0^2)=',xsec/pi

      open(unit=16,file='output.txt')
      write(16,*),xsec/pi
      close(16)


      write(6,*)'xsec_Kr(pi a_0^2)=',xsec_Kr/pi
!
      stop
      end
!
      complex*16 function clgam(x)
***********************************************************************
*                                                                     *
*            COMPLEX LOG-GAMMA(Z)                                     *
*          STERLING'S APPROXIMATION                                   *
*          Z COMPLEX                                                  *
***********************************************************************
      implicit real*8 (a-h,o-z)
      complex*16  x,y,z,pz,pcz,zlg,gam,sub
      dimension q(11)
      data q
     1 / 8.3333333333333d-02,-2.7777777777778d-03, 7.9365079365079d-04,
     2  -5.9523809523810d-04, 8.4175084175084d-04,-1.9175269175269d-03,
     3   6.4102564102564d-03,-2.9550653594771d-02, 1.7964437236883d-01,
     4  -1.3924322169059d0    , 13.4028640441684d0/
      data hlp /9.18938533204673d-1/
      data pi/3.14159265358979d0/
      ineg = 0
      z = x
*        if re(z) < 0  gam(z) = pi/sin(pi*z)*gam(1-z)

      if(dreal(z).lt.0.d0) then
      ineg = 1
      pz = pi*z
      pcz = pi/cdsin(pz)
      sub = cdlog(pcz)
      z = 1.d0 - z
      endif

      gam = dcmplx(0.d0,0.d0)

*       if re(z) < 6  gam(z) = gam(z+k)/(z+k-1)*...(z)

 100   if(dreal(z).lt.6.0d0) then

*  special treatment of log for re(z) = 0

      if(dreal(z).eq.0.d0) then
      yy = dimag(z)
      if(yy.gt.0.d0)  then
      zlg = dcmplx(dlog(yy),0.5d0*pi)
      elseif(yy.lt.0.d0) then
      zlg = dcmplx(dlog(-yy),-0.5d0*pi)
      else
      stop 'clgam : z = 0 or negative integer'
      endif
      else
      zlg = cdlog(z)
      endif

      gam = gam - zlg
      z = z + 1.0d0

      go to 100
      endif

      gam = gam + (z-.5d0)*cdlog(z) - z + hlp

      y = z
      do 200 i = 1,11
      y=y/(z*z)
      gam = gam + q(i)*y
 200  continue

      clgam = gam

      if(ineg.ne.0)
     *   clgam = sub - gam

      return
      end
!
! -------------------------------------------
!
      Subroutine coulfg ( rho,eta,minl,maxl,fc,fcp,gc,gcp,accur, iret )
      implicit real*8(a-h,o-z)
c
c     computes regular and irregular coulomb wavefunctions.
c
c       this subroutine returns the regular and irregular coulomb
c     wavefunctions and their derivatives for specific values of
c     eta and rho and for a range of l's.  the range of l's (orbital
c     angular momentum) need not begin at l = 0 and, in general,
c     it will be less time consuming to request only those l's
c     that are actually desired.  the conventions of abramowitz
c     and stegun (handbook of mathematical functions) are used.
c
c       this subroutine is an adaptation of the manchester subroutine
c     of the same name (but slightly different argument list) that
c     was published in computer physics communications 8, 377 (1974).
c     by barnett, feng, steed and goldfarb.  the continued fractions
c     introduced in that article are used in the range rho >
c     rho(turn) but completely different procedures are used
c     for rho < rho(turn).  the new procedures allow results
c     for rho < rho(turn) of accuracy comparable to those for
c     rho > rho(turn) in comparable time.
c
c       one of the continued fractions is not convergent for very
c     small rho and thus this routine will find only f and f' for
c     rho < .005.  in addition the routine becomes rather slow for
c     rho < .1.  for very large rho ( >> 1000 ), another continued
c     fraction converges very slowly and the subroutine will fail
c     to converge in the maximum allowed number of iterations for
c     rho > 9500 with the exception that larger rho are possible
c     if eta is near rho/2.
c
c       aside from the above limitations, the subroutine has a very
c     large range of eta,  rho, and l for which it will return
c     reliable values in reasonable times.  it has been tested for
c     -2000 < eta < 5000,  1e-6 < rho < 10000;  0 < l < 1000 .
c     the tests consisted of comparasons to existing routines to
c     locate coding errors and of comparasons to quadruple precision
c     results to determine the numerical stability of the algorithms.
c     however, no verification of the correctness of g or g' for
c     eta < 0 has been made.
c
c
c     arguments (all floating point arguments are double precision) -
c
c     rho - the value of the radial coordinate at which f, f', g, and
c           g' are desired.  0 < rho < 9500 is required.  if
c           rho < .005, only f and f' will be found.
c
c     eta - the value of the charge parameter to be used.  the routine
c           has been tested for -2000 < eta < 5000.  eta = 0 will
c           result in the computation of the spherical bessel
c           functions times rho.
c
c     minl - the minimum value of l for which the functions are desired.
c     maxl - the maximum value of l.
c
c     fc, fcp, gc, and gcp - these one dimensional arrays will be
c           set to the computed values of f, f', g, and g' respectively.
c           each array must be of length at least maxl+1 and will
c           be set as
c             array(l+1) = function(l)   minl =< l =< maxl .
c           depending on eta and rho, the first minl elements of
c           each array may be used for intermediate computations and
c           their values are not predictable upon exit from coulfg.
c
c     accur - the desired accuracy of convergence of the continued
c          fractions and other series.  in general this will be
c          the relative accuracy of the final values except near
c          a zero of one of the functions.  accur should be in
c          the range
c              1e-6 < accur < 1e-15 ,
c          if it is not, one of the above two values will be used.
c
c     iret - this integer is set to a return code to indicate that
c          coulfg was successful or the reason for failure:
c        0 - all o.k.; all functions found.
c        1 - some o.k., those for higher l's have under/over-flowed
c        2 - rho < .005; only f, f' were computed.
c        3 - rho < .005; only f, f' were computed and some of them
c                        underflowed.
c        4 - all will under/over-flow; no results defined.
c        5 - failed to avoid divide by zero in f'/f loop
c        6 - nonconvergence of r
c        7 - nonconvergence of p + iq
c        8 - nonconvergence of maclauren series
c        9 - nonconvergence of the taylor series about rho = 2*eta
c     iret's of 5 to 9 should be reported to steve pieper ( x4523 )
c
c       the following table gives sample execution times on the
c     /75 and /195.  times are for accur = 1d-14 and for
c     minl = maxl = 0.  in general execution time is weakly dependant
c     on accur and minl.  if maxl > minl, the times required for
c     each additional l are approximatley 0.01 millisec on the /195
c     and 0.15 millisec on the /75.  times for eta < 0 are comparable to
c     to those for !eta!.  the final column of the table gives
c     the number of bits lost due to round off and truncation
c     errors and indicates the maximum precision possible on a
c     given machine.  the /360 and /370 have 56 bits of precision
c     (16.8 decimal places) and each 3.3 bits lost represents
c     one decimal place lost.
c
c      eta    rho     time in milliseconds    bits lost
c                       /195      /75
c
c        0.      .01     .15       1.2            3
c        0.     1.       .13       1.2            3
c        0.   100.      1.0       11.            11
c        0.  1000.      7.6       81.
c        1.      .1     5.4       81.             9
c        1.     1.      1.0       12.5            8
c        1.    10.       .33       3.7            3
c        1.   100.      1.1       11.
c        1.  1000.      8.1       81.
c       10.     1.      1.5       20.             6
c       10.    10.       .8        8.3            6
c      100.    75.      1.9       21.            10
c      100.   100.      2.0       21.            10
c      100.  1000.      7.3       74.            13
c     1000.  5000.     29.       385.            16
c
c
c     may 23, 1976 - revised version by s. pieper.
c     nov 19, 1976 - fix choice of root for rho << eta, eta > 10
c
c
      logical frstsw
      dimension fc(1),fcp(1),gc(1),gcp(1),
     1   iveryb(4)
c
c     vrybig is the result of an overflow
c     big is representative of max number, small is its inverse
c     smalln is  log(small)  (must be accurate)
c     precis is about  100*(machine precision)
c     precln is  > -!log(machine precision)!
c     prert3 is like the cuberoot of the machine precision
c     data concerning vax machine
c
c      data prert3 / 2.1544d-5 /,  precis / 1.d-12 /,precln /32.2362d0/
c      data pi / 3.1415926535897932d0/,
c     2   big/1.7976931348623158D+308/,
c     3   small/2.2250738585072014d-308/,
c     1   vrybig/1.D+309/smalln/-7.08396418532264079d+2/
      data prert3 / 2.d-5 /,  precis / 1.d-13 /,precln /32.2d0/
      data pi / 3.1415926535897932d0/,
     2   big/1.D+300/,
     3   small/1.d-300/,
     1   vrybig/1.D+301/smalln/-6.90775527898213705d+2/
c
c     coulomb wavefunctions calculated at rho = rho by the
c     continued-fraction method of steed   minl,maxl are actual l-values
c     see barnett feng steed and goldfarb computer physics common 1974
c
c
c     here we limit accuracy to reasonable values for the machine
c
      acc  = accur
      acc = dmax1( acc, precis )
      acc = dmin1( acc, prert3 )
C      vrybig=big*big
c
      lmax = maxl
      lmin = minl
      lmin1= lmin + 1
      xll1 = lmin*lmin1
      eta2 = eta*eta
c
c     determin which region we are in
c
c     for rho < .45, q of p+iq is poorly determined so we don't use it.
c     except that for large negative eta the maclauren series also
c     has problems
c
      if ( rho .gt. .45d0 )  go to 20
      if ( eta .ge. 0 )  go to 10
      if ( -eta*rho .gt. 7 )  go to 20
c
c     for rho < .005, we only return f and f' since the p+iq recursion
c     is very slowly convergent  ( for very small eta it is possible
c     to go to smaller rho ( rho > 1e+4*eta )  but we ignor that here.
c
 10   if ( rho .gt. .005d0 )  go to 60
      igoto = 5
      go to 70
c
 20   turn = eta + sqrt(eta2 + xll1)
      igoto = 1
      if ( rho .ge. turn-1.e-4 )  go to 100
c
c     we are inside the turning point for minl, can we get outside
c     of it by reducing minl. (this is always possible for
c     eta < 0).
c
      if ( rho .lt. eta+abs(eta) )  go to 60
c
c     yes, do so - this is the same as  rho > rho(turn)  except
c     we generate some extra f(l), g(l) for  l < minl
c
      lmin = .5*( sqrt(1+4*((rho-eta)**2-eta2)) - 1 )
      lmin1 = lmin + 1
      go to 80
c
c     must use a different method to suppliment the bad i q
c     value.  always start with lmin = 0 for simplicity.
c
c     note only eta > 0 gets to here ( except when rho < .45 )
c
 60   igoto = 2
 70   lmin = 0
      lmin1 = 1
      if ( eta .lt. 10  .or.  rho .le. eta )  go to 80
      igoto = 3
c
 80   xll1 = lmin*lmin1
c
c     here we compute  f'/f  for l = maxl
c     we then recurse down to lmin to generate the unnormalized f's
c     this section is used for all rho.
c
 100  pl   = lmax + 1
      rhouse = rho
 105  plsave = pl
 110  frstsw = .true.
c     continued fraction for  r = fp(maxl)/f(maxl)
      r  = eta/pl + pl/rhouse
      dq  = (eta*rhouse)*2.0 + 6*pl**2
      dr = 12*pl + 6
      del = 0.0
      d   = 0.0
      f   = 1.0
      x   = (pl*pl - pl + (eta*rhouse))*(2.0*pl - 1.0)
      ai = rhouse*pl**2
      di = (2*pl+1)*rhouse
c
c     loop and converge on r
c
      do 139  i = 1, 100000
         h = (ai + rhouse*eta2)*(rhouse - ai)
         x   = x + dq
         d = d*h + x
c
c     if we pass near a zero of the divisor, start over at
c     larger lmax
c
         if ( abs(d) .gt. prert3*abs(dr) )  go to 130
         pl = pl + 1
         if ( pl .lt. plsave+10 )  go to 110
         iret = 5
         return
c
 130     d = 1/d
         dq = dq + dr
         dr = dr + 12
         ai = ai + di
         di = di + 2*rhouse
         del =  del*(d*x - 1.0)
         if (frstsw) del = -rhouse*(pl*pl + eta2)*(pl + 1.0)*d/pl
         frstsw = .false.
         r  = r + del
         if(d.lt.0.0) f = -f
         if ( abs(del) .lt. abs(r*acc) )  go to 140
 139  continue
      iret = 6
      return
c
c     r has converged;  did we increase lmax
c
 140  if ( pl .eq. plsave )  go to 160
c
c     recurse down on r to lmax
c     here the only part of f that is of interest is the sign
c
      pl = pl-1
 150     d = eta/pl + pl/rhouse
         f = (r+d)*f
         r = d - (1+eta2/pl**2)/(r+d)
         pl = pl - 1
         if ( pl .gt. plsave )  go to 150
c
c     now have r(lmax, rho) or if igoto=4, r(lmin, 2*eta)
c
 160  if ( igoto .eq. 4 )  go to 210
      fc (lmax+1) = f
      fcp(lmax+1) = f*r
      if( lmax.eq.lmin) go to 200
c     downward recursion to lmin for f and fp, arrays gc,gcp are storage
      l  = lmax
      pl = lmax
      ar = 1/rho
      do 189 lp  = lmin1,lmax
         gc (l+1) = eta/pl + pl*ar
         gcp(l+1) = sqrt( (eta/pl)**2 + 1 )
         fc (l)   = (gc(l+1)*fc(l+1) + fcp(l+1))/gcp(l+1)
         fcp(l)   =  gc(l+1)*fc(l)   - gcp(l+1)*fc(l+1)
         pl = pl - 1
         l  = l - 1
c
c     if we are getting near an overflow, renormalize everything down
c
         if ( abs(fc(l+1)) .lt. big )  go to 189
         do 179  ll = l, lmax
            fc(ll+1) = small*fc(ll+1)
            fcp(ll+1) = small*fcp(ll+1)
 179     continue
 189  continue
      f  = fc (lmin1)
      r = fcp(lmin1)/f
c
c     here we find
c        p + iq  =  (g'+if')/(g+if)
c     this section is used in all cases except when
c        15 < eta < rho < 2*eta
c
 200  if ( igoto .eq. 3 )  go to 500
      if ( igoto .eq. 5 )  go to 400
c
c     now obtain p + i.q for lmin from continued fraction (32)
c     real arithmetic to facilitate conversion to ibm using real*8
 210  p  = 0.0
      q  = rhouse - eta
      pl = 0.0
      ar = -(eta2 + xll1)
      ai =   eta
      br = q + q
      bi = 2.0
      wi = eta + eta
      dr =   br/(br*br + bi*bi)
      di =  -bi/(br*br + bi*bi)
      dp = -(ar*di + ai*dr)
      dq =  (ar*dr - ai*di)
c
c     loop and converge on p + iq
c
 230     p  =  p + dp
         q  =  q + dq
         pl = pl + 2.0
         ar = ar + pl
         ai = ai + wi
         bi = bi + 2.0
         d  = ar*dr - ai*di + br
         di = ai*dr + ar*di + bi
         t  = 1.0/(d*d + di*di)
         dr =  t*d
         di = -t*di
         h  = br*dr - bi*di - 1.0
         x  = bi*dr + br*di
         t  = dp*h  - dq*x
         dq = dp*x  + dq*h
         dp = t
         if(pl.gt.46000.) go to 920
         if(abs(dp)+abs(dq).ge.(abs(p)+abs(q))*acc) go to 230
      p  = p/rhouse
      q  = q/rhouse
c
c     we now have  r  and  p+iq,  is this enough
c
      if ( igoto .eq. 2 )  go to 400
c
c     solve for fp,g,gp and normalise f  at l=lmin
c
c     since this is for  rho > rho(turn), f and g are reasonable
c     numbers
c
      x = (r-p)/q
      fmag = sqrt( 1/(q*(1+x**2)) )
      w = fmag/abs(f)
      f = w*f
      g = f*x
      gp = r*g - 1/f
      if ( igoto .eq. 4 )  go to 600
      go to 800
c
c     here   rho < eta  or  rho < 2*eta < 20  or  rho < .45
c     we use the maclauren series to get  f( l=0, eta, rho )
c
c     first compute  rho*c(l=0, eta)
c
 400  c = 2*pi*eta
      if ( abs(c) .gt. .5 )  go to 410
c
c     use maclaurin expansion of  x / (exp(x)-1)
c
      x = 0
      t = 1
      ar = 1
      br = c
      ai = 1
      c = 1
 405     ai = ai + 1
         ar = ar*br/ai
         c = c + ar
         if ( abs(ar) .ge. acc*c )  go to 405
      c = 1/c
      go to 430
c
c     here eta is not tiny.
c
 410  if ( eta .gt. 0 )  go to 420
      c = -c
      x = 0
      t = 1
      go to 425
 420  x = -smalln - pi*eta
      t = small
 425  if ( c .lt. precln )  c = c / (1-exp(-c))
 430  c = rho*sqrt(c)
      b1 = 1
      b2 = eta*rho
      sum = b1 + b2
      ai = 6
      di = 6
      do 449  i = 1, 10000
         b3 = ( (2*eta*rho)*b2 - (rho**2)*b1 ) / ai
         ai = ai + di
         di = di + 2
         sum = sum + b3
         stop = abs(b1) + abs(b2) + abs(b3)
         b1 = b2
         b2 = b3
         if ( abs(sum) .lt. big )  go to 445
         x = x - smalln
         sum = sum*small
         b1 = b1*small
         b2 = b2*small
 445     if ( stop .lt. acc*abs(sum) )  go to 450
 449  continue
      iret = 8
      return
c
 450  sum = ( c*exp(x)*sum ) * t
c
c     did it underflow
c
      if ( sum .eq. 0 )  go to 900
c
c     we now have  f (=sum),  r,  and p  ( p only if rho > .005 )
c     use the wronskian as the 4th condition
c
      w = sum/f
      f = sum
      if ( igoto .eq. 5 )  go to 850
      x = (r-p)*f
      if ( abs(x) .gt. prert3 )  go to 480
c
c     here f**3 and f**4 terms are less than machine precision
c
      g = 1/x
      gp = p*g
      go to 800
c
c     here we must include f**3, f**4; we must also worry about
c     which sign of the root to use.
c     the positive root applies for g > f; else the negative root
c     we use q in determining which is correct
c
 480  b1 = .5/x
      b2 = b1 * sqrt(1-4*(x*f)**2)
      g = b1 + b2
c     g > f in all of region 2 for eta > 0
      if ( eta .ge. 0 )  go to 490
      sum = 1/q - f**2
      gp = b1 - b2
      if ( abs(g**2-sum) .gt. abs(gp**2-sum) )  g = gp
 490  gp = p*g - x*f/g
      go to 800
c
c     eta > 15  and  eta < rho < 2*eta
c
c     we find g and g' for lmin, rho=2*eta using the above method
c     consisting of r, p+iq, and w.
c
 500  rhouse = eta+eta
      pl = lmin+1
      igoto = 4
      go to 105
c
c     now we have  g, g'  at the turning point, go in using taylor
c
c
 600  del = rhouse - rho
      b1 = g
      b2  = -del*gp
      b3 = 0
      g = b1+b2
      accr = acc/2
      delinv = -1/del
      dfactr = 3*delinv
      x = del/rhouse
      ai = x+x
      di = ai+ai
      ar = 6
      dr = 6
      do 639  i = 1, 10000
         s = ( ai*b3 + (x*del**2)*b1 ) / ar
         ar = ar + dr
         dr = dr + 2
         ai = ai + di
         di = di + 2*x
         g = g + s
         gp = gp + dfactr*s
         if ( g .ge. vrybig )  go to 900
         dfactr = dfactr + delinv
         b1 = b2
         b2 = b3
         b3 = s
         if ( s .lt. accr*g )  go to 650
 639  continue
      iret = 9
      return
c
c     here we have  r = f'/f,  g,  g'
c     use wronskian as the 4th condition
c
 650  f = fc(lmin1)
      r = fcp(lmin1)/f
      sum = 1/(r*g-gp)
      w = sum/f
      f = sum
c
c     we now have  f, r = f'/f, g, g'  at lmin
c
c     upward recursion from gc(lmin) and gcp(lmin),stored values are rho
c     renormalise fc,fcp for each l-value
c
 800  gc (lmin1) = g
      gcp(lmin1) = gp
      fc(lmin1) = f
      fcp(lmin1) = r*f
      iret = 0
      if(lmax.eq.lmin)  return
      do  829  l = lmin1,lmax
         t        = gc(l+1)
         gc (l+1) = (gc(l)*gc (l+1) - gcp(l))/gcp(l+1)
         gcp(l+1) =  gc(l)*gcp(l+1) - gc(l+1)*t
         fc (l+1) = w*fc (l+1)
 829     fcp(l+1) = w*fcp(l+1)
 840  if ( abs(fc(lmax+1))+abs(fcp(lmax+1)) .eq. 0 )  iret = iret+1
      return
c
c     rho < .005;  we cannot find p or q and so return only f, f'.
c
 850  fc(lmin1) = f
      fcp(lmin1) = r*f
      iret = 2
      if ( lmax .eq. lmin )  return
      do 859  l = lmin1, lmax
         fc(l+1) = w*fc(l+1)
         fcp(l+1) = w*fcp(l+1)
 859  continue
      go to 840
c
c     f and g are out of the machine exponent range for lmin.
c     it will be even worse for  l > lmin  so give up and return
c
 900  iret = 4
      return
c
c     p + iq failed to converge
 920  iret = 7
      return
c
      end
c
c -------------------------------------------------------
c
      function rint (f,na,nb,nq,h)
      implicit doubleprecision(a-h,o-z)
c
c  this program calculates the integral of the function f from point na
c  to point nb using a nq points quadrature ( nq is any integer between
c  1 and 14 ).  h is the grid size.
c                                      written by c. c. j. roothaan
c
      dimension c(105),c1(25),c2(80),d(14),f(nb)
      equivalence (c1(1),c(1)),(c2(1),c(26))
      data c1/1.d0,2.d0,1.d0,23.d0,28.d0,9.d0,25.d0,20.d0,31.d0,8.d0,
     &1413.d0,1586.d0,1104.d0,1902.d0,475.d0,1456.d0,1333.d0,1746.d0,
     &944.d0,1982.d0,459.d0,119585.d0,130936.d0,89437.d0,177984.d0/
      data c2/54851.d0,176648.d0,36799.d0,122175.d0,111080.d0,156451.d0,
     &46912.d0,220509.d0,29336.d0,185153.d0,35584.d0,7200319.d0,
     &7783754.d0,5095890.d0,12489922.d0,-1020160.d0,16263486.d0,
     &261166.d0,11532470.d0,2082753.d0,7305728.d0,6767167.d0,9516362.d0,
     &1053138.d0,18554050.d0,-7084288.d0,20306238.d0,-1471442.d0,
     &11965622.d0,2034625.d0,952327935.d0,1021256716.d0,636547389.d0,
     &1942518504.d0,-1065220914.d0,3897945600.d0,-2145575886.d0,
     &3373884696.d0,-454944189.d0,1637546484.d0,262747265.d0,
     &963053825.d0,896771060.d0,1299041091.d0,-196805736.d0,
     &3609224754.d0,-3398609664.d0,6231334350.d0,-3812282136.d0,
     &4207237821.d0,-732728564.d0,1693103359.d0,257696640.d0,
     &5206230892907.d0,5551687979302.d0,3283609164916.d0,
     &12465244770050.d0,-13155015007785.d0,39022895874876.d0,
     &-41078125154304.d0,53315213499588.d0,-32865015189975.d0,
     &28323664941310.d0,-5605325192308.d0,9535909891802.d0,
     &1382741929621.d0,5252701747968.d0,4920175305323.d0,
     &7268021504806.d0,-3009613761932.d0,28198302087170.d0,
     &-41474518178601.d0,76782233435964.d0,-78837462715392.d0,
     &81634716670404.d0,-48598072507095.d0,34616887868158.d0,
     &-7321658717812.d0,9821965479386.d0,1360737653653.d0/
      data d/2.d0,2.d0,24.d0,24.d0,1440.d0,1440.d0,120960.d0,120960.d0,
     &7257600.d0,7257600.d0,958003200.d0,958003200.d0,5230697472000.d0,
     &5230697472000.d0/
      a=0.0
      l=na
      m=nb
      i=nq*(nq+1)/2
      do 1 j=1,nq
      a=a+c(i)*(f(l)+f(m))
      l=l+1
      m=m-1
    1 i=i-1
      a=a/d(nq)
      do 2 n=l,m
    2 a=a+f(n)
      rint=a*h
      return
      end
c
c ------------------------------------------------------
c
      doubleprecision function d3j30(l1,l2,l3)
      implicit real*8(a-b,d-h,o-z)
      implicit complex*16 (c)
      ip = (l1+l2+l3)/2
      c1 = dcmplx(dfloat(l1+l2-l3+1),0.d0)
      c2 = dcmplx(dfloat(l2+l3-l1+1),0.d0)
      c3 = dcmplx(dfloat(l3+l1-l2+1),0.d0)
      c4 = dcmplx(dfloat(l1+l2+l3+2),0.d0)
      c5 = dcmplx(dfloat(ip-l1+1),0.d0)
      c6 = dcmplx(dfloat(ip-l2+1),0.d0)
      c7 = dcmplx(dfloat(ip-l3+1),0.d0)
      c8 = dcmplx(dfloat(ip+1))
      cldel = clgam(c1)+clgam(c2)+clgam(c3)-clgam(c4)
      cl3j = 0.5d0*cldel+clgam(c8)-clgam(c5)-clgam(c6)-clgam(c7)
      d3j30 = ((-1.d0)**ip)*dreal(cdexp(cl3j))
      return
      end
!
!**********************************************************************
!
      function ifact(n)
**********************************************************************
*     FONCTION  FACTORIELLE  N : ENTIER                              *
**********************************************************************
      ifact = 1
      if (n.lt.2) return
      do 21 i = 1,n
      ifact = i*ifact
 21   continue
      return
      end
!
!**********************************************************************
!
!
      subroutine interp (pi,qi,ri,mi,ni,pf,qf,rf,mf,nf,n)
c
c  this program uses aitken's method to do interpolations
c
c  pi,qi : input functions
c  ri    : original radial grid where pi & qi are defined
c  mi,ni : first & last points where pi & qi are defined (mi < ni)
c  pf,qf : output interpolated functions
c  rf    : new radial grid where pf and qf are defined
c  mf,nf : first & last points where pf & qf are defined (mf <= nf)
c  n     : order of the interpolation method (n >= 2)
c
c  note that ri & rf must be monotonically increasing
c
c           written by k. t. cheng        version : 01/31/83
c _____________________________________________________________________
c
      implicit doubleprecision(a-h,o-z)
      dimension pi(1),qi(1),ri(1),pf(1),qf(1),rf(1)
      dimension px(20),qx(20),rx(20),py(20),qy(20)
c
      np=min0(n,ni-mi+1)
      nq=(np+1)/2+1
      nr=ni-np
      ii=mi
      do 70 i=mf,nf
      do 10 j=ii,ni
      if(ri(j)-rf(i)) 10,15,20
 10   continue
      ii=ni
      ia=nr
      go to 30
 15   ii=j
      pf(i)=pi(j)
      qf(i)=qi(j)
      go to 70
 20   ii=j
      ia=max0(ii-nq,mi-1)
      ia=min0(ia,nr)
c
c  aitken's np point method of interpolation
c
 30   do 40 k=1,np
      px(k)=pi(ia+k)
      qx(k)=qi(ia+k)
 40   rx(k)=ri(ia+k)
      do 60 j=2,np
      do 50 k=j,np
      dr=(rf(i)-rx(j-1))/(rx(k)-rx(j-1))
      py(k)=(1.d0-dr)*px(j-1)+dr*px(k)
 50   qy(k)=(1.d0-dr)*qx(j-1)+dr*qx(k)
      do 60 k=j,np
      px(k)=py(k)
 60   qx(k)=qy(k)
      pf(i)=px(np)
      qf(i)=qx(np)
 70   continue
      return
      end
!
! ------------------------------------------
!
      DOUBLE PRECISION FUNCTION P1S(Zn,R)
c
c Radial part of the 1s Bound wave function for the
c Coulomb potential. Hydrogen-like ions with charge
c Zn-1.
c
      IMPLICIT REAL*8 (A-H,O-Z)
      Zr = Zn*R
      P1S = 2.D0*Zr*DSQRT(Zn)*DEXP(-Zr)
      RETURN
      END
!
! -------------------------------------------
!
      DOUBLE PRECISION FUNCTION P2S(Zn,R)
c
c Radial part of the 2s Bound wave function for the
c Coulomb potential. Hydrogen-like ions with charge
c Zn-1.
c P(r)=R(r)*r
c
      IMPLICIT REAL*8 (A-H,O-Z)
      Zr = Zn*R
      P2S = DSQRT(Zn/2.d0)*Zr*(1.d0-0.5d0*Zr)*DEXP(-0.5d0*Zr)
      RETURN
      END
!
! -------------------------------------------
!
      DOUBLE PRECISION FUNCTION P2P(Zn,R)
c
c Radial part of the 2p Bound wave function for the
c Coulomb potential. Hydrogen-like ions with charge
c Zn-1.
c P(r)=R(r)*r
c
      IMPLICIT REAL*8 (A-H,O-Z)
      Zr = Zn*R
      P2P = DSQRT(Zn/2.d0)*Zr**2*DEXP(-0.5d0*Zr)/dsqrt(3.d0)/2.d0
      RETURN
      END
!
! -------------------------------------------
!
      DOUBLE PRECISION FUNCTION P3S(Zn,R)
c
c Radial part of the 3s Bound wave function for the
c Coulomb potential. Hydrogen-like ions with charge
c Zn-1.
c P(r)=R(r)*r
c
      IMPLICIT REAL*8 (A-H,O-Z)
      Zr = Zn*R
      P3S=2.d0*dsqrt(Zn/3.d0)**3*(1.d0-2.d0*Zr/3.d0+2.d0*Zr**2/27.d0)*
     $dexp(-Zr/3.d0)*R
      RETURN
      END
!
! -------------------------------------------
!
      DOUBLE PRECISION FUNCTION P3P(Zn,R)
c
c Radial part of the 3p bound wave function for the
c Coulomb potential. Hydrogen-like ions with charge
c Zn-1.
c P(r)=R(r)*r
c
      IMPLICIT REAL*8 (A-H,O-Z)
      rho=2.d0*Zn*R/3.d0
      R3P = (1.d0/9.d0/dsqrt(6.d0))*rho*(4.d0-rho)
      R3P = R3P*(Zn**1.5d0)*dexp(-0.5d0*rho)
      P3P = R3P*R
      RETURN
      END
!
! -------------------------------------------
!
      DOUBLE PRECISION FUNCTION P3D(Zn,R)
c
c Radial part of the 3d bound wave function for the
c Coulomb potential. Hydrogen-like ions with charge
c Zn-1.
c P(r)=R(r)*r
c
      IMPLICIT REAL*8 (A-H,O-Z)
      Zr = Zn*R
      xx = 4.d0/DSQRT(10.d0)/27.d0/3.d0/dsqrt(3.d0)
      P3D = xx*Zn*dsqrt(Zn)*(Zr**2)*r*DEXP(-Zr/3.d0)
      RETURN
      END
!
! -------------------------------------------
!
      DOUBLE PRECISION FUNCTION P4S(Zn,R)
c
c Radial part of the 4s Bound wave function for the
c Coulomb potential. Hydrogen-like ions with charge
c Zn-1.
c P(r)=R(r)*r
c
      IMPLICIT REAL*8 (A-H,O-Z)
      rho=2.d0*Zn*R/4.d0
      rho2=rho*rho
      rho3=rho2*rho
      P4S = R*(1.d0/96.d0)*(24.d0-36.d0*rho+12.d0*rho2-rho3)*
     &      dexp(-rho/2.d0)*Zn**(1.5d0)
      RETURN
      END
!
! -------------------------------------------
!
      DOUBLE PRECISION FUNCTION P4P(Zn,R)
c
c Radial part of the 4p Bound wave function for the
c Coulomb potential. Hydrogen-like ions with charge
c Zn-1.
c P(r)=R(r)*r
c
      IMPLICIT REAL*8 (A-H,O-Z)
      rho=2.d0*Zn*R/4.d0
      rho2=rho*rho
      P4P = R*(1.d0/32.d0/dsqrt(15.d0))*rho*(20.d0-10.d0*rho+rho2)*
     &      dexp(-rho/2.d0)*Zn**(1.5d0)
      RETURN
      END
!
! -------------------------------------------
!
      DOUBLE PRECISION FUNCTION P4D(Zn,R)
c
c Radial part of the 4d Bound wave function for the
c Coulomb potential. Hydrogen-like ions with charge
c Zn-1.
c P(r)=R(r)*r
c
      IMPLICIT REAL*8 (A-H,O-Z)
      rho=2.d0*Zn*R/4.d0
      rho2=rho*rho
      P4D = R*(1.d0/96.d0/dsqrt(5.d0))*rho2*(6.d0-rho)*
     &      dexp(-rho/2.d0)*Zn**(1.5d0)
      RETURN
      END
!
! -------------------------------------------
!
      DOUBLE PRECISION FUNCTION P4F(Zn,R)
c
c Radial part of the 4f Bound wave function for the
c Coulomb potential. Hydrogen-like ions with charge
c Zn-1.
c P(r)=R(r)*r
c
      IMPLICIT REAL*8 (A-H,O-Z)
      rho=2.d0*Zn*R/4.d0
      rho2=rho*rho
      rho3=rho2*rho
      P4F = R*(1.d0/96.d0/dsqrt(35.d0))*rho3*
     &      dexp(-rho/2.d0)*Zn**(1.5d0)
      RETURN
      END
!
! -------------------------------------------
!
      DOUBLE PRECISION FUNCTION P5S(Zn,R)
c
c Radial part of the 5s Bound wave function for the
c Coulomb potential. Hydrogen-like ions with charge
c Zn-1.
c P(r)=R(r)*r
c
      IMPLICIT REAL*8 (A-H,O-Z)
      rho=2.d0*Zn*R/5.d0
      rho2=rho*rho
      rho3=rho2*rho
      rho4=rho3*rho
      P5S = R*(1.d0/300.d0/dsqrt(5.d0))*(120.d0-240.d0*rho+120.d0*rho2
     &         -20.d0*rho3+rho4)*
     &         dexp(-rho/2.d0)*Zn**(1.5d0)
      RETURN
      END
c
c ----------------------------------------------
c
      subroutine rad(l,xr,xn)
      implicit real*8 (a-h,o-z)
      parameter(imax=500,lmax=50)
      dimension FC(LMAX+2),FCP(LMAX+2),GC(LMAX+2),GCP(LMAX+2)
      dimension r(imax),rp(imax),w(imax),wn(imax)
      common/coul/zc,pf
      common/hydr/nh,lh
      etaf=-zc/pf
      accur=1.d-15
      h=0.022d0
      r0=0.005d0
      r(1)=0.d0
      rp(1)=r0
      do 1 i=1,imax
      rp(i)=dexp((i-1)*h)
      r(i)=r0*(rp(i)-1.d0)
      rp(i)=r0*rp(i)
1     continue
      do 2 i=1,imax
      rhof=pf*r(i)
      call coulfg (rhof,etaf,l,l,fc,fcp,gc,gcp,accur,iret)
      PP=P1S(1.d0,r(i))
      w(i)=r(i)*fc(l+1)*PP*rp(i)
      wn(i)=PP*PP*rp(i)
2     continue
      xr=rint(w,1,imax,7,h)
      xn=rint(wn,1,imax,7,h)
      return
      end
