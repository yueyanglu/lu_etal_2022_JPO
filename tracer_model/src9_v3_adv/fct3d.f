      subroutine fct3d(iord,fld,fco,fc)
      use mod_xc                ! HYCOM communication interface
      use mod_za
      use mod_tracer
c
c --- fully 3-d version of advfct.f
c     
c     Main input vars:
c       iord    - 1: donnor cell; 2: complete with antidiffusive fluxes
c       fco/fc  - layer thicknesses at previous and new time step
c       u/vflx  - isopycnal mass fluxes [m3] (times delta1)
c       diaflx  - mass flux through layer interfaces [m], + downward
c                 diaflx(k) is the flux across the bottom of layer k
c       fld     - 'tracer', transported mixing ratio, like sal or temp
c
c       scal   - grid cell size
c       scali  - inverse of scale
c 
      implicit none
c      include 'dimensions.h'
c      include 'geopar.h'
c
      real, dimension(1-nbdy:ii+nbdy,1-nbdy:jj+nbdy,kdm) ::
     $     fld,fco,fc,vertfx,vertdv
      real, dimension(1-nbdy:ii+nbdy,1-nbdy:jj+nbdy)     ::
     $     fmx,fmn,flp,fln,flx,fly,uan,van,flxdiv
      real, dimension(1-nbdy:jj+nbdy)               :: clipj, vlumj
      real a(kdm),b(kdm),c(kdm),athird,dx,fcdx,yl,yr
      real onemu,q,clip,vlume,amount,bfore,after,slab,dslab,thkchg,
     .     fluxdv,epsil,bforej(1-nbdy:jj+nbdy),afterj(1-nbdy:jj+nbdy)
      integer iord,ip1,im1,jp1,jm1,kp,jaa,margin,i1,i2,j1,j2,
     $     ia,ib,ja,jb,l
      logical wrap,recovr
      data recovr/.false./
      data (athird = 1./3.)
c
c --- if iord=1, scheme reduces to simple donor cell scheme.
      parameter (epsil = 1.e-11, onemu = 1.e-6)
c
c 103  format (i3,4f14.1)
c
c --- get vertical flux by summing over upstream slab of thickness -w-
c
      margin = nbdy - 1
c      print*,'Tracer 1:',fld(itest,jtest,1)
!$OMP PARALLEL DO PRIVATE(jb,amount,slab,dslab,kp,a,b,c,dx,fcdx,yl,yr)
      do 26 j = 1-margin, jj+margin
      do 26 l = 1, isp(j)
c
c --- fill massless cells with data from layer above or below
      do 17 k = kk-1, 1, -1
      do 17 i = ifp(j,l), ilp(j,l)
 17   fld(i,j,k)=(fld(i,j,k)*fco(i,j,k)+fld(i,j,k+1)*onemu)
     .                 /(           fco(i,j,k)+             onemu)
      do 18 k=2,kk
      do 18 i=ifp(j,l),ilp(j,l)
 18   fld(i,j,k)=(fld(i,j,k)*fco(i,j,k)+fld(i,j,k-1)*onemu)
     .                 /(           fco(i,j,k)+             onemu)
c
      do 26 i=ifp(j,l),ilp(j,l)
c
c --- fit 0th, 1st, or 2nd deg. polynomial to fld in each cell
      a(1 )=fld(i,j,1 )
      b(1 )=0.
      c(1 )=0.
      a(kk)=fld(i,j,kk)
      b(kk)=0.
      c(kk)=0.
      do 16 k=2,kk-1
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c --- piecewise constant method:
ccc      a(k)=fld(i,j,k)
ccc      b(k)=0.
ccc      c(k)=0.
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c --- piecewise linear method:
c --- fit linear function a+bx to fld in each cell (-.5 < x < +.5)
ccc      a(k)=fld(i,j,k)
ccc      b(k)=0.
ccc      if (fld(i,j,k).le.min(fld(i,j,k-1),fld(i,j,k+1)) .or.
ccc     .    fld(i,j,k).ge.max(fld(i,j,k-1),fld(i,j,k+1))) then
ccc        b(k)=0.
ccc      else if ((fld(i,j,k+1)-fld(i,j,k-1))*(fld(i,j,k-1)+fld(i,j,k+1)
ccc     .  -2.*fld(i,j,k)).gt.0.) then
ccc        b(k)=fld(i,j,k)-fld(i,j,k-1)
ccc      else
ccc        b(k)=fld(i,j,k+1)-fld(i,j,k)
ccc      end if
ccc      c(k)=0.
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c --- piecewise parabolic method:
c --- fit parabola a+bx+cx^2 to fld in each cell (-.5 < x < +.5)
      yl=.5*(fld(i,j,k-1)+fld(i,j,k))
      yr=.5*(fld(i,j,k+1)+fld(i,j,k))
      a(k)=1.5*fld(i,j,k)-.25*(yl+yr)
      b(k)=yr-yl
      c(k)=6.*(.5*(yl+yr)-fld(i,j,k))
      if (abs(yr-yl) .lt. 6.*abs(.5*(yl+yr)-fld(i,j,k))) then
c --- apex of parabola occurs inside interval [-.5,+.5], implying an over-
c --- or undershoot situation. change curve to prevent over/undershoots.
         if (abs(yr-yl) .gt. 2.*abs(.5*(yl+yr)-fld(i,j,k))) then
c     --- put apex of parabola on edge of interval [-.5,+.5]
            if ((yr-yl)*(.5*(yl+yr)-fld(i,j,k)) .gt. 0.) then
c --- apex at x=-.5
               a(k)=.25*(3.*fld(i,j,k)+yl)
               c(k)=3.*(fld(i,j,k)-yl)
               b(k)=c(k)
            else
c     --- apex at x=+.5
               a(k)=.25*(3.*fld(i,j,k)+yr)
               c(k)=3.*(fld(i,j,k)-yr)
               b(k)=-c(k)
            end if
         else                   !  -1/6 < x < +1/6
c --- moving apex won't help. replace parabola by constant.
            a(k)=fld(i,j,k)
            b(k)=0.
            c(k)=0.
         end if
      end if
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 16   continue
c
      do 23 k=1,kk-1
      slab=onemu
      if (diaflx(i,j,k).lt.0.) then ! + downward ! interface moves down
         amount=slab*fld(i,j,k+1)
         kp=k
 24      kp=kp+1
         if (slab.ge.-diaflx(i,j,k)) goto 23
         if (fco(i,j,kp).gt.0.) then
            dslab=min(slab+fco(i,j,kp),-diaflx(i,j,k))
     .           -min(slab            ,-diaflx(i,j,k))
            dx=dslab/fco(i,j,kp)
            fcdx=a(kp)
     .           +b(kp)*.5*(dx-1.) !  not needed in pcm
     .           +c(kp)*(.25-dx*(.5-dx*athird)) !  not needed in pcm,plm
            amount=amount+fcdx*dslab
            slab=slab+dslab
         end if
         if (kp.lt.kk) go to 24
      else if (diaflx(i,j,k).gt.0.) then ! interface moves up
         amount=slab*fld(i,j,k)
         kp=k+1
 25      kp=kp-1
         if (slab.ge.diaflx(i,j,k)) goto 23
         if (fco(i,j,kp).gt.0.) then
            dslab=min(slab+fco(i,j,kp), diaflx(i,j,k))
     .           -min(slab            , diaflx(i,j,k))
            dx=dslab/fco(i,j,kp)
            fcdx=a(kp)
     .           +b(kp)*.5*(1.-dx) !  not needed in pcm
     .           +c(kp)*(.25-dx*(.5-dx*athird)) !  not needed in pcm,plm
            amount=amount+fcdx*dslab
            slab=slab+dslab
         end if
         if (kp.gt.2) go to 25
      end if
 23   vertfx(i,j,k)=diaflx(i,j,k)*amount/slab
c     
      vertfx(i,j,kk)=0.         !  don't allow flux through bottom
      vertdv(i,j,1)=vertfx(i,j,1)
      do 26 k=2,kk
      vertdv(i,j,k)=vertfx(i,j,k)-vertfx(i,j,k-1)
c
 26   continue
c
      bfore=0.
      after=0.
c
      do 4 k=1,kk
c
!$OMP PARALLEL DO SHARED(k)
      do 14 j=1,jj
      bforej(j)=0.
      do 14 l=1,isp(j)
c      do 14 i=ifp(j,l),ilp(j,l)
      do 14 i=1,ii
 14   if (ip(i,j).ne.0) bforej(j)=bforej(j)+fld(i,j,k)*
     $        fco(i,j,k)*scp2(i,j)
c
c --- compute low-order antidiffusive (high- minus low-order) fluxes
c
      margin = nbdy - 2
c
!$OMP PARALLEL DO PRIVATE(ja,jaa,jb,q,ia,ib) SHARED(k)
      do 11 j=1-margin,jj+margin
c      ja=mod(j-2+jj,jj)+1
      ja=j-1
c      jb=mod(j     ,jj)+1
      jb=j+1
c      jaa=mod(j-3+jj,jj)+1
      jaa=j-2
      margin=nbdy-2
      do 2 l=1,isu(j)
      i1=max(1-margin,ifu(j,l))
      i2=min(ii+margin,ilu(j,l))
c      do 2 i=ifu(j,l),ilu(j,l)
      do 2 i=i1,i2
      if (uflx(i,j,k).ge.0.) then
        q=fld(i-1,j,k)
      else
        q=fld(i  ,j,k)
      end if
      flx(i,j)=uflx(i,j,k)*q
      q=fld(i,j,k)+fld(i-1,j,k) !  2nd order
      if (ip(i+1,j)+iu(i-1,j).eq.2)
     .     q=1.125*q-.125*(fld(i+1,j,k)+fld(i-2,j,k)) !  4th order
 2    uan(i,j)=.5*q*uflx(i,j,k)-flx(i,j)
c
      do 3 l=1,isv(j)
      i1=max(1-margin,ifv(j,l))
      i2=min(ii+margin,ilv(j,l))
      do 3 i=i1,i2
c      do 3 i=ifv(j,l),ilv(j,l)
      if (vflx(i,j,k).ge.0.) then
        q=fld(i,ja ,k)
      else
        q=fld(i,j  ,k)
      end if
      fly(i,j)=vflx(i,j,k)*q
      q=fld(i,ja ,k)+fld(i,j,k) !  2nd order
      if (ip(i,jb )+iv(i,ja).eq.2)
     .     q=1.125*q-.125*(fld(i,jb ,k)+fld(i,j-2,k)) !  4th order
 3    van(i,j)=.5*q*vflx(i,j,k)-fly(i,j)
c
      margin=nbdy-2
      do 11 l=1,isp(j)
      i1=max(1-margin,ifp(j,l))
      i2=min(ii+margin,ilp(j,l))
c      do 11 i=ifp(j,l),ilp(j,l)
      do 11 i=i1,i2
c      ia=max( 1,i-1)
      ia=i-1
      if (ip(ia,j).eq.0) ia=i
c      ib=min(ii,i+1)
      ib=i+1
      if (ip(ib,j).eq.0) ib=i
c      ja=mod(j-2+jj,jj)+1
      ja=j-1
      if (ip(i,ja).eq.0) ja=j
c      jb=mod(j     ,jj)+1
      jb=j+1
      if (ip(i,jb).eq.0) jb=j
      fmx(i,j)=max(fld(i,j,k),
     .   fld(ia,j,k),fld(ib,j,k),fld(i,ja,k),fld(i,jb,k))
      fmn(i,j)=min(fld(i,j,k),
     .   fld(ia,j,k),fld(ib,j,k),fld(i,ja,k),fld(i,jb,k))
      if (k.lt.kk) then
        if (diaflx(i,j,k  ).lt.0.) then
          fmx(i,j)=max(fmx(i,j),vertfx(i,j,k  )/diaflx(i,j,k  ))
          fmn(i,j)=min(fmn(i,j),vertfx(i,j,k  )/diaflx(i,j,k  ))
        end if
      end if
      if (k.gt.1) then
        if (diaflx(i,j,k-1).gt.0.) then
          fmx(i,j)=max(fmx(i,j),vertfx(i,j,k-1)/diaflx(i,j,k-1))
          fmn(i,j)=min(fmn(i,j),vertfx(i,j,k-1)/diaflx(i,j,k-1))
        end if
      end if
 11   continue
      margin=nbdy-2
c      print*,'ii+margin 1:',ii+margin
!$OMP PARALLEL DO
      do 22 j=1-margin,jj+margin
      do 22 l=1,isp(j)
      i1=max(1-margin,ifp(j,l))
      i2=min(ii+margin,ilp(j,l))
c      flx(ifp(j,l)  ,j)=0.
c      flx(ilp(j,l)+1,j)=0.
c      uan(ifp(j,l)  ,j)=0.
c      uan(ilp(j,l)+1,j)=0.
      flx(i1  ,j)=0.
      flx(i2+1,j)=0.
      uan(i1  ,j)=0.
      uan(i2+1,j)=0.
  22  continue
c$OMP END PARALLEL DO
c
!$OMP PARALLEL DO PRIVATE(j)
      do 33 i=1-margin,ii+margin
c     wrap=jfv(i,1).eq.1	! true if j=1 and j=jj are both water points
      do 33 l=1,jsp(i)
c      j1=max(1-margin,jfp(i,l))
c      j2=min(jj+margin,jlp(i,l))
      if (j.gt.1 .or. .not.wrap) then
         fly(i,j)=0.
         van(i,j)=0.
      end if
c      j=mod(jlp(i,l),jj)+1
      j=min(jj+nbdy,jlp(i,l)+1)
      if (j.gt.1 .or. .not.wrap) then
         fly(i,j)=0.
         van(i,j)=0.
      end if
   33 continue
c$OMP END PARALLEL DO
c
cdiag i=itest
cdiag j=jtest
cdiag write (lp,101) 'advem(1)',i,j,k,fld(i-1,j,k),uflx(i,j,k),
cdiag. fld(i,j-1,k),vflx(i,j,k),fld(i,j,k),vflx(i,j+1,k),fld(i,j+1,k),
cdiag.  uflx(i+1,j,k),fld(i+1,j,k)
 101  format(a,2i5,i3,f18.3/1pe39.2/0pf19.3,1pe11.2,0pf9.3,
     .1pe11.2,0pf9.3/1pe39.2/0pf39.3)
c
      margin=nbdy-2
!$OMP PARALLEL DO PRIVATE(ib, jb,q,amount) SHARED(k)
      do 61 j=1-margin,jj+margin
c      jb=mod(j     ,jj)+1
      jb=j+1
      if (recovr) vlumj(j)=0.
      if (recovr) clipj(j)=0.
      do 61 l=1,isp(j)
      i1=max(1-margin,ifp(j,l))
      i2=min(ii+margin,ilp(j,l))
      do 61 i=i1,i2
c      do 61 i=ifp(j,l),ilp(j,l)
      flxdiv(i,j)=(flx(i+1,j)-flx(i,j)+fly(i,jb )-fly(i,j))*scp2i(i,j)
c
cdiag if (i.eq.itest .and. j.eq.jtest)
cdiag. write (lp,'(2i5,i3,a,4f10.5,1pe9.2)') i,j,k,'  fc,fco,divs:',
cdiag.  fc(i,j,k),fco(i,j,k),flxdiv(i,j),vertdv(i,j,k),
cdiag.  fc(i,j,k)-fco(i,j,k)+flxdiv(i,j)+vertdv(i,j,k)
c
      q=fld(i,j,k)*fco(i,j,k)-flxdiv(i,j)-vertdv(i,j,k)
      amount=max(fmn(i,j)*fc(i,j,k),min(q,fmx(i,j)*fc(i,j,k)))
      if (recovr) then
        vlumj(j)=vlumj(j)+scp2(i,j)*fc(i,j,k)
        clipj(j)=clipj(j)+(q-amount)*scp2(i,j)
      end if
 61   fld(i,j,k)=(fld(i,j,k)*onemu+amount)/(onemu+fc(i,j,k))
c$OMP END PARALLEL DO
c
cdiag call findmx(ip,fld(1,1,k),idm,ii1,jj,'fld after 61')
c
      if (iord.le.1) go to 100
c
c --- at each grid point, determine the ratio of the largest permissible
c --- pos. (neg.) change in -fld- to the sum of all incoming (outgoing) fluxes
c
c      print*,'ii+margin 2:',ii+margin
!$OMP PARALLEL DO PRIVATE(jb) SHARED(k)
      do 12 j=1-margin,jj+margin
c      jb=mod(j     ,jj)+1
      jb=j+1
      do 12 l=1,isp(j)
      i1=max(1-margin,ifp(j,l))
      i2=min(ii+margin,ilp(j,l))
      do 12 i=i1,i2
c      do 12 i=ifp(j,l),ilp(j,l)
      flp(i,j)=(fmx(i,j)-fld(i,j,k))*fc(i,j,k)
     ./((max(0.,uan(i,j))-min(0.,uan(i+1,j))
     .  +max(0.,van(i,j))-min(0.,van(i,jb ))+epsil)*scp2i(i,j))
c
 12   fln(i,j)=(fmn(i,j)-fld(i,j,k))*fc(i,j,k)
     ./((min(0.,uan(i,j))-max(0.,uan(i+1,j))
     .  +min(0.,van(i,j))-max(0.,van(i,jb ))-epsil)*scp2i(i,j))
c$OMP END PARALLEL DO
c
c---- limit antidiffusive fluxes
c
      print*,'ii+margin 3:',ii+margin
!$OMP PARALLEL DO PRIVATE(ja,clip)
      do 8 j=1-margin,jj+margin
c      ja=mod(j-2+jj,jj)+1
      ja=j-1
c
      do 7 l=1,isu(j)
      i1=max(1-margin,ifu(j,l))
      i2=min(ii+margin,ilu(j,l))
      do 7 i=i1,i2
c      do 7 i=ifu(j,l),ilu(j,l)
      if (uan(i,j).ge.0.) then
        clip=min(1.,flp(i,j),fln(i-1,j))
      else
        clip=min(1.,fln(i,j),flp(i-1,j))
      end if
 7    flx(i,j)=uan(i,j)*clip
c
      do 8 l=1,isv(j)
      i1=max(1-margin,ifv(j,l))
      i2=min(ii+margin,ilv(j,l))
      do 8 i=i1,i2
c      do 8 i=ifv(j,l),ilv(j,l)
      if (van(i,j).ge.0.) then
        clip=min(1.,flp(i,j),fln(i,ja ))
      else
        clip=min(1.,fln(i,j),flp(i,ja ))
      end if
 8    fly(i,j)=van(i,j)*clip
c$OMP END PARALLEL DO
c
cdiag i=itest
cdiag j=jtest
cdiag write (lp,101) 'advem(2)',i,j,k,fld(i-1,j,k),uflx(i,j,k),
cdiag. fld(i,j-1,k),vflx(i,j,k),fld(i,j,k),vflx(i,j+1,k),fld(i,j+1,k),
cdiag.  uflx(i+1,j,k),fld(i+1,j,k)
c
      margin=nbdy-2
!$OMP PARALLEL DO PRIVATE(jb,amount,q) SHARED(k)
      do 62 j=1-margin,jj+margin
      jb=j+1
c      jb=mod(j     ,jj)+1
      do 62 l=1,isp(j)
      i1=max(1-margin,ifp(j,l))
      i2=min(ii+margin,ilp(j,l))
      do 62 i=i1,i2
c      do 62 i=ifp(j,l),ilp(j,l)
      flxdiv(i,j)=(flx(i+1,j)-flx(i,j)+fly(i,jb )-fly(i,j))*scp2i(i,j)
      q=fld(i,j,k)*fc(i,j,k)-flxdiv(i,j)
      amount=max(fmn(i,j)*fc(i,j,k),min(q,fmx(i,j)*fc(i,j,k)))
      if (recovr) clipj(j)=clipj(j)+(q-amount)*scp2(i,j)
 62   fld(i,j,k)=(fld(i,j,k)*onemu+amount)/(onemu+fc(i,j,k))
c$OMP END PARALLEL DO
c
cdiag call findmx(ip,fld(1,1,k),idm,ii1,jj,'fld after 62')
c
 100  continue
c
c --- recover 'clipped' amount and return to field layer by layer
c
      if (recovr) then
        vlume=0.
        clip=0.
c
        do 19 j=1-margin,jj+margin
        vlume=vlume+vlumj(j)
 19     clip=clip+clipj(j)
c
        if (vlume.ne.0.) then
          clip=clip/vlume
          write (lp,'(i2,a,1pe11.3)') k,'  fld drift in fct3d',-clip
!$OMP PARALLEL DO SHARED(k)
          do 13 j=1-margin,jj+margin
          do 13 l=1,isp(j)
          i1=max(1-margin,ifp(j,l))
          i2=min(ii+margin,ilp(j,l))
          do 13 i=i1,i2
c          do 13 i=ifp(j,l),ilp(j,l)
 13       fld(i,j,k)=fld(i,j,k)+clip
c$OMP END PARALLEL DO
        end if
      end if
      print*,'Tracer 7:',fld(itest,jtest,1)
c
!$OMP PARALLEL DO SHARED(k)
      do 15 j=1,jj
      afterj(j)=0.
      do 15 l=1,isp(j)
      i1=max(1-margin,ifp(j,l))
      i2=min(ii+margin,ilp(j,l))
c      do 15 i=ifp(j,l),ilp(j,l)
      do 15 i=1,ii
 15   if (ip(i,j).ne.0) afterj(j)=afterj(j)+fld(i,j,k)*
     $        fc(i,j,k)*scp2(i,j)
c$OMP END PARALLEL DO
      do j=1,jj
        bfore=bfore+bforej(j)
        after=after+afterj(j)
      end do
c
 4    continue
c
c      if (bfore.ne.0.)
c     . write (lp,'(a,1p,3e14.6,e11.1)') 'fct3d conservation:',
c     .  bfore,after,after-bfore,(after-bfore)/bfore
      q=1.
      if (after.ne.0.) q=bfore/after
c      write (lp,'(a,f11.6)') 'fct3d: multiply fld field by',q
c
!$OMP PARALLEL DO
      do 20 j=1-margin,jj+margin
      do 20 k=1,kk
      do 20 l=1,isp(j)
      i1=max(1-margin,ifp(j,l))
      i2=min(ii+margin,ilp(j,l))
      do 20 i=i1,i2
c      do 20 i=ifp(j,l),ilp(j,l)
 20   fld(i,j,k)=fld(i,j,k)*q
c$OMP END PARALLEL DO
      print*,'Tracer 8:',bfore,after,fld(itest,jtest,1)
c
      return
      end
c
c
c> Revision history:
c>
c> Feb. 2005 - added plm,ppm options to vertical flux calc'n (so far: pcm)
c> Oct. 2017 - rewritten using buffer points
