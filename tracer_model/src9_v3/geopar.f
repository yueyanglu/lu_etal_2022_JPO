      subroutine geopar
      use mod_tracer  ! HYCOM post-processing tracer array interface
      use mod_xc      ! HYCOM communication interface
      use mod_za      ! HYCOM array I/O interface
c
      implicit none
c
c --- set up model parameters related to geography
c
      real      hmina,hminb,hmaxa,hmaxb
      integer   l,ia,ja,margin
      character preambl(5)*79,cline*80
      real epsil,valland, depmax
      real, allocatable :: drec(:)
      epsil=1.e-12
c
c --- read basin depth array
c
      write (lp,'(2a)') ' reading bathymetry file from ',
     &                  'regional.depth.[ab]'
      call xcsync(flush_lp)
      open (unit=9,file='regional.depth.b',status='old')
      read (     9,'(a79)')  preambl
      read (     9,'(a)')    cline
      i = index(cline,'=')
      read (cline(i+1:),*)   hminb,hmaxb
      close(unit=9)
      write (lp,'(/(1x,a))') preambl,cline
c
      call zaiopf('regional.depth.a','old', 9)
      call zaiord(depths,ip,.false., hmina,hmaxa, 9)     
      call zaiocl(9)
c
      if     (abs(hmina-hminb).gt.abs(hminb)*1.e-4 .or.
     &        abs(hmaxa-hmaxb).gt.abs(hmaxb)*1.e-4     ) then
        write(lp,'(/ a / a,1p3e14.6 / a,1p3e14.6 /)')
     &    'error - .a and .b files not consistent:',
     &    '.a,.b min = ',hmina,hminb,hmina-hminb,
     &    '.a,.b max = ',hmaxa,hmaxb,hmaxa-hmaxb
c        stop '(geopar)'
      endif
c
c --- start out with masks as land everywhere
c
c!$OMP PARALLEL DO PRIVATE(j,i)
c!$OMP&         SCHEDULE(STATIC,jblk)
      do j=1-nbdy,jj+nbdy
        do i=1-nbdy,ii+nbdy
          ip(i,j)=0
          iu(i,j)=0
          iv(i,j)=0
        enddo
      enddo
c
c --- mass points are defined where water depth is greater than zero
c --- but for GLBa0.08 there are large numbers over land to be masked later.
c!$OMP PARALLEL DO PRIVATE(j,i)
c!$OMP&         SCHEDULE(STATIC,jblk)
      do j=1-nbdy,jj+nbdy
         do i=1-nbdy,ii+nbdy
            if (depths(i,j).gt.0.0) then
               ip(i,j)=1
            endif
         enddo
      enddo
c
c --- u,v points are located halfway between any 2 adjoining mass points
c!$OMP PARALLEL DO PRIVATE(j,i)
c!$OMP&         SCHEDULE(STATIC,jblk)
      do j=1-nbdy,jj+nbdy
        do i=1-nbdy,ii+nbdy
           ia=max(1-nbdy,i-1)
           if (ip(ia,j).gt.0.and.ip(i,j).gt.0) then
              iu(i,j)=1
           endif
           ja=max(1-nbdy,j-1)
           if (ip(i,ja).gt.0.and.ip(i,j).gt.0) then
              iv(i,j)=1
           endif
        enddo
      enddo
c
      margin = nbdy-1
c --- determine sea-only i-1, i+1, j-1, and j+1 indexes
c!$OMP PARALLEL DO PRIVATE(j,i)
c!$OMP&         SCHEDULE(STATIC,jblk)
      do j=1-margin,jj+margin
         do i=1-margin,ii+margin
c ---     W,E,S,N if sea, otherwise center point
            if     (ip(i-1,j).ne.0) then
               ipim1(i,j) = i-1
            else
               ipim1(i,j) = i
            endif
            if     (ip(i+1,j).ne.0) then
               ipip1(i,j) = i+1
            else
               ipip1(i,j) = i
            endif
            if     (ip(i,j-1).ne.0) then
               ipjm1(i,j) = j-1
            else
               ipjm1(i,j) = j
            endif
            if     (ip(i,j+1).ne.0) then
               ipjp1(i,j) = j+1
            else
               ipjp1(i,j) = j
            endif
         enddo
      enddo
c
c --- read grid location,spacing,coriolis arrays
c
      write (lp,'(2a)') ' reading grid file from ',
     &                  'regional.grid.[ab]'
      call xcsync(flush_lp)
      open (unit=9,file='regional.grid.b',status='old')
      read (9,*) i
      read (9,*) j
      if     (i.ne.idm .or. j.ne.jdm) then
        write(lp,'(/ a /)')
     &    'error - wrong array size in grid file'
        call xcstop('(geopar)')
               stop '(geopar)'
      endif
      read (9,*) mapflg
      write(lp,*)'mapflg', mapflg
c     read (9,'(a)') cline
c     write (lp,'(a)') cline(1:len_trim(cline))
c
      call zaiopf('regional.grid.a','old', 9)
c
      do k= 1,12
        read (9,'(a)') cline
        i = index(cline,'=')
        read (cline(i+1:),*) hminb,hmaxb
        write (lp,'(a)') cline(1:len_trim(cline))
        call xcsync(flush_lp)
c
        if     (k.eq.1) then
           call zaiord(plon, ip,.true., hmina,hmaxa, 9)
        elseif (k.eq.2) then
           call zaiord(plat, ip,.true., hmina,hmaxa, 9)
          do i= 1,2
             read (9,'(a)') cline
             call zaiosk(9)
          enddo
        elseif (k.eq.3) then
           call zaiord(ulon, iu,.true., hmina,hmaxa, 9)
        elseif (k.eq.4) then
           call zaiord(ulat, iu,.true., hmina,hmaxa, 9)
        elseif (k.eq.5) then
           call zaiord(vlon, iv,.true., hmina,hmaxa, 9)
        elseif (k.eq.6) then
           call zaiord(vlat, iv,.true., hmina,hmaxa, 9)
           read (9,'(a)') cline
           call zaiosk(9)
        elseif (k.eq.7) then
           call zaiord(scpx, ip,.true., hmina,hmaxa, 9)
        elseif (k.eq.8) then
           call zaiord(scpy, ip,.true., hmina,hmaxa, 9)
           do i= 1,2
              read (9,'(a)') cline
              call zaiosk(9)
          enddo
        elseif (k.eq.9) then
           call zaiord(scux, iu,.true., hmina,hmaxa, 9)
        elseif (k.eq.10) then
           call zaiord(scuy, iu,.true., hmina,hmaxa, 9)
        elseif (k.eq.11) then
           call zaiord(scvx, iv,.true., hmina,hmaxa, 9)
        else
           call zaiord(scvy, iv,.true., hmina,hmaxa, 9)
        endif
c
        if     (abs(hmina-hminb).gt.abs(hminb)*1.e-4 .or.
     &          abs(hmaxa-hmaxb).gt.abs(hmaxb)*1.e-4     ) then
c          write(lp,'(/ a / a,1p3e14.6 / a,1p3e14.6 /)')
c     &      'error - .a and .b files not consistent:',
c     &      '.a,.b min = ',hmina,hminb,hmina-hminb,
c     &      '.a,.b max = ',hmaxa,hmaxb,hmaxa-hmaxb
c          call xcstop('(geopar)')
c                 stop '(geopar)'
        endif
      enddo
c!$OMP PARALLEL DO PRIVATE(j,i)
c!$OMP&         SCHEDULE(STATIC,jblk)
      do j= 1-nbdy,jj+nbdy
         do i= 1-nbdy,ii+nbdy
            scu2(i,j)=scux(i,j)*scuy(i,j)
            scv2(i,j)=scvx(i,j)*scvy(i,j)
            scp2(i,j)=scpx(i,j)*scpy(i,j)
            scp2i(i,j)=1/max(scp2(i,j),epsil)
         enddo
      enddo
      write(6,*)'itest, jtest, scux,scvy, scp2,scp2i',itest,jtest
      write(6,'(10e13.5)')
     .     scux(itest,jtest),
     .     scvy(itest,jtest),scp2(itest,jtest),scp2i(itest,jtest)
c --- determine if this is a global periodic grid
c
      call flush(6)
c --- is the domain periodic in longitude?
      lperiod=.false.
      write(6,*)'lperiod', lperiod
      call flush(6)
       
c --- is the domain periodic in latitude?
      larctic=.false.
      write(6,*)'larctic', larctic
      call flush(6)
      call zaiocl(9)
c
c!$OMP PARALLEL DO PRIVATE(j,i)
c!$OMP&         SCHEDULE(STATIC,jblk)
      do j= 1-nbdy,jj+nbdy
        do i= 1-nbdy,ii+nbdy
           if     (depths(i,j).gt.1.e10) depths(i,j) = 0.0
        enddo
      enddo
c recompute ip, iu,iv
      do j=1-nbdy,jj+nbdy
         do i=1-nbdy,ii+nbdy
            ip(i,j)=0
            iu(i,j)=0
            iv(i,j)=0
            if (depths(i,j).gt.0.)then
               ip(i,j)=1
            endif
         enddo
      enddo
c --- u,v points are located halfway between any 2 adjoining mass points
c!$OMP PARALLEL DO PRIVATE(j,i)
c!$OMP&         SCHEDULE(STATIC,jblk)
      do j=1-nbdy,jj+nbdy
        do i=1-nbdy,ii+nbdy
           ia=max(1-nbdy,i-1)
           if (ip(ia,j).gt.0.and.ip(i,j).gt.0) then
              iu(i,j)=1
           endif
           ja=max(1-nbdy,j-1)
           if (j.gt.1.and.ip(i,ja).gt.0.and.ip(i,j).gt.0) then
              iv(i,j)=1
           endif
        enddo
      enddo
c     write(6,*)'iu',5,100,iu(5,100)
c     write(6,*)'iu',5,500,iu(5,500)
c     write(6,*)'iv',5,100,iv(5,100)
c     write(6,*)'iv',5,500,iv(5,500)
c --- determine loop indices for mass and velocity points
      call indxi(ip,ifp,ilp,isp)
c     write(6,*)'in geopar, isp ifp ilp'
c     do j-1-nbdy,jj 
c       write(6,*)j,isp(j),(ifp(j,l),ilp(j,l),l=1,isp(j))
c     enddo
      call indxj(ip,jfp,jlp,jsp)
c     write(6,*)'in geopar, jsp jfp jlp'
c     do i=1,ii
c       write(6,*)i,jsp(j),(jfp(i,l),jlp(i,l),l=1,jsp(i))
c     enddo
      call indxi(iu,ifu,ilu,isu)
      call indxj(iu,jfu,jlu,jsu)
      call indxi(iv,ifv,ilv,isv)
      call indxj(iv,jfv,jlv,jsv)
c
c!$OMP PARALLEL DO PRIVATE(j,i)
c!$OMP&         SCHEDULE(STATIC,jblk)
      do j= 1-nbdy,jj+nbdy
        do i= 1-nbdy,ii+nbdy
           ia=max(1-nbdy,i-1)
           depthu(i,j)=min(depths(i,j),depths(ia,j))
           ja=max(1-nbdy,j-1)
           depthv(i,j)=min(depths(i,j),depths(i,ja))
        enddo
      enddo
c print first column
c     write(6,*)'in geopar'
      print *, 'onem', onem
      print *, 'ip,iu,iv, depths,depthu,depthv i=1,3 j=1000'
      do i=1,3
      print *, ip(i,1000),iu(i,1000),iv(i,1000) 
      print *, depths(i,1000),depthu(i,1000),depthv(i,1000)
      enddo
c print last column
      print *, 'ip,iu,iv, depths,depthu,depthv i=ii-2,ii,j=1000'
      do i=ii-2,ii
      print *, ip(i,1000),iu(i,1000),iv(i,1000) 
      print *, depths(i,1000),depthu(i,1000),depthv(i,1000)
      enddo

c print if,il for j=1000
      j=1000
      print *,'ifp,ilp,ifu,ilu,ifv,ilv j=1000'
      print *, ifp(1000,1:isp(j))
      print *, ilp(1000,1:isp(j))
      print *, ifu(1000,1:isu(j))
      print *, ilu(1000,1:isu(j))
      print *, ifv(1000,1:isv(j))
      print *, ilv(1000,1:isv(j))
      print *, 'ip,iu,iv, depths,depthu,depthv,itest,jtest',itest,jtest
      i=itest
      j=jtest
      print *, ip(i,j),iu(i,j),iv(i,j) 
      print *, depths(i,j),depthu(i,j),depthv(i,j)

      print *, 'end geopar'
c
      return
      end

      subroutine indxi(ipt,if,il,is)
c      use mod_xc  ! HYCOM communication interface
      use mod_tracer  ! HYCOM post-processing tracer array interface
      implicit none
c
      integer, dimension (1-nbdy:ii+nbdy,1-nbdy:jj+nbdy) :: ipt
      integer, dimension (1-nbdy:jj+nbdy,ms) :: if,il
      integer, dimension (1-nbdy:jj+nbdy) :: is
c
c --- input array ipt contains 1 at grid point locations, 0 elsewhere
c --- output is arrays if, il, is  where
c --- if(j,k) gives row index of first point in column j for k-th section
c --- il(j,k) gives row index of last point
c --- is(j) gives number of sections in column j (maximum: ms)
c
      integer last,i0,j0
c
c!$OMP PARALLEL DO PRIVATE(j,i)
c!$OMP&         SCHEDULE(STATIC,jblk)
      do j=1-nbdy,jj+nbdy
         is(j) = 0
         do k=1,ms
            if(j,k) = 0
            il(j,k) = 0
         end do
c
         k=1
         last = ipt(1-nbdy,j)
         if     (last .eq. 1) then
            if(j,k) = 1-nbdy
         endif
         do i=2-nbdy,ii+nbdy
            if      (last .eq. 1 .and. ipt(i,j) .eq. 0) then
               il(j,k) = i-1
               k = k+1
            elseif (last .eq. 0 .and. ipt(i,j) .eq. 1) then
               if     (k .gt. ms) then
                  write(lp,'(a,i5)')  'indxi problem on proc ',mnproc
                  write(lp,'(a,2i5)')
     &   ' error in indxi -- ms too small at i,j =',i0+i,j0+j
c                  call xchalt('(indxi)')
                  stop '(indxi)'
               endif
               if(j,k) = i
            endif
            last = ipt(i,j)
         enddo
         if     (last .eq. 1) then
            il(j,k) = ii+nbdy
            is(j) = k
         else
            is(j) = k-1
         endif
      enddo
      call xcsync(no_flush)
      return
      end

c
      subroutine indxj(jpt,jf,jl,js)
      use mod_xc  ! HYCOM communication interface
      use mod_tracer  ! HYCOM post-processing tracer array interface
      implicit none
c
      integer, dimension (1-nbdy:ii+nbdy,1-nbdy:jj+nbdy) :: jpt
      integer, dimension (1-nbdy:ii+nbdy,ms) :: jf,jl
      integer, dimension (1-nbdy:ii+nbdy) :: js
c
c --- input array jpt contains 1 at grid point locations, 0 elsewhere
c --- output is arrays jf, jl, js  where
c --- jf(i,k) gives column index of first point in row i for k-th section
c --- jl(i,k) gives column index of last point
c --- js(i) gives number of sections in row i (maximum: ms)
c
      integer last,i0,j0
c
c!$OMP PARALLEL DO PRIVATE(j,i)
c!$OMP&         SCHEDULE(STATIC,jblk)
      do i=1-nbdy,ii+nbdy
        js(i) = 0
        do k=1,ms
          jf(i,k) = 0
          jl(i,k) = 0
        end do
c
        k=1
        last = jpt(i,1-nbdy)
        if     (last .eq. 1) then
          jf(i,k) = 1-nbdy
        endif
        do j=2-nbdy,jj+nbdy
           if      (last .eq. 1 .and. jpt(i,j) .eq. 0) then
              jl(i,k) = j-1
              k = k+1
           elseif (last .eq. 0 .and. jpt(i,j) .eq. 1) then
              if     (k .gt. ms) then
                 write(lp,'(a,i5)')  'indxj problem on proc ',mnproc
                 write(lp,'(a,2i5)')
     &             ' error in indxj -- ms too small at i,j =',i0+i,j0+j
c                 call xchalt('(indxj)')
                 stop '(indxj)'
              endif
              jf(i,k) = j
           endif
           last = jpt(i,j)
        enddo
        if     (last .eq. 1) then
           jl(i,k) = jj+nbdy
           js(i) = k
        else
           js(i) = k-1
        endif
      enddo
      call xcsync(no_flush)
      return
      end
