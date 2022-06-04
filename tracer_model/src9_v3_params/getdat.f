c=====================================================================
c          SUBTOUTINE
c=====================================================================
      subroutine forday(dtime,yrflag, iyear,iday,ihour)
      implicit none
c
      double precision dtime
      integer          yrflag, iyear,iday,ihour
c
c --- converts model day to "calendar" date (year,julian-day,hour).
c
      double precision dtim1,day
      integer          iyr,nleap
c
      if     (yrflag.eq.0) then
c ---   360 days per model year, starting Jan 16
        iyear =  int((dtime+15.001d0)/360.d0) + 1
        iday  =  mod( dtime+15.001d0 ,360.d0) + 1
        ihour = (mod( dtime+15.001d0 ,360.d0) + 1.d0 - iday)*24.d0
c
      elseif (yrflag.eq.1) then
c ---   366 days per model year, starting Jan 16
        iyear =  int((dtime+15.001d0)/366.d0) + 1
        iday  =  mod( dtime+15.001d0 ,366.d0) + 1
        ihour = (mod( dtime+15.001d0 ,366.d0) + 1.d0 - iday)*24.d0
c
      elseif (yrflag.eq.2) then
c ---   366 days per model year, starting Jan 01
        iyear =  int((dtime+ 0.001d0)/366.d0) + 1
        iday  =  mod( dtime+ 0.001d0 ,366.d0) + 1
        ihour = (mod( dtime+ 0.001d0 ,366.d0) + 1.d0 - iday)*24.d0
c
      elseif (yrflag.eq.3) then
c ---   model day is calendar days since 01/01/1901
        iyr   = (dtime-1.d0)/365.25d0
        nleap = iyr/4
        dtim1 = 365.d0*iyr + nleap + 1.d0
        day   = dtime - dtim1 + 1.d0
        if     (dtim1.gt.dtime) then
          iyr = iyr - 1
        elseif (day.ge.367.d0) then
          iyr = iyr + 1
        elseif (day.ge.366.d0 .and. mod(iyr,4).ne.3) then
          iyr = iyr + 1
        endif
        nleap = iyr/4
        dtim1 = 365.d0*iyr + nleap + 1.d0
c
        iyear =  1901 + iyr
        iday  =  dtime +0.01d0 - dtim1 + 1
        ihour = (dtime - dtim1 + 1.01d0 - iday)*24.d0
c
      else
        stop '(forday)'
      endif
      return
      end
c=====================================================================
c          SUBTOUTINE
c=====================================================================
      subroutine getdat_name(flnm, dtime,yrflag,iyear)
      use mod_za     ! HYCOM array I/O interface
c
      implicit none
c
      character*120    flnm
      double precision dtime
      integer          yrflag
c
c --- modify the archive filename for the specified model day.
c --- flnm is assumed to be of the form uvdp_offline_????_???.b.
c
      integer   iyear,iday,ihour,l
c
      call forday(dtime,yrflag, iyear,iday,ihour)
      l = len_trim(flnm)
      write(flnm(l-9:l-6),'(i4.4)') iyear
      write(flnm(l-4:l-2),'(i3.3)') iday
c      write(flnm(l-3:l-2),'(i2.2)') ihour
      return
      end
c=====================================================================
c          SUBTOUTINE
c=====================================================================
      subroutine getdat_archo(flnm, flg)
c     READ forcing fields: u/vflx, dpm and dpi ('dpio')
c     flg: '0' read full forcing, '1' read CS foc in 'mfoc_path'
      use mod_tracer  ! HYCOM post-processing tracer archive array interface
      use mod_za      ! HYCOM array I/O interface
c
      implicit none
c
      character*120 flnm
      character*120 ffoc_path, mfoc_path, dir_flnm ! file name
     &      , mldfnm, mldfnm_ful
      double precision time
      integer          flg,i1,i2,j1,j2
     &                 ,l2

c
c --- read model fields from an archive file specific for teh offline model
c --- HYCOM 2.0 array I/O archive file.
c 
c --- from the files for the offline model, read :
c ---      instantaneous dp
c ---      mean velocity and dp to create momentum fluxes. 
c ---      mixed layer depth and other 2d fields if needed
c
      character cline*80
      character preambl(5)*79
      character cvarin*6
      real      hminb,hmaxb,thet
      integer   idmtst,irec,iversn,ios,jdmtst,l,layer,ni,nstep,recsh,nf
      integer n2drec, nrecl
      real, allocatable :: drec(:)
c
      data ni/14/
c
      l = len_trim(flnm)
      if (timflg.gt.0) write(flnm(l-9:l-6),'(a4)') 'clim'
c
      i1=1 -nbdy
      i2=ii+nbdy
      j1=1 -nbdy
      j2=jj+nbdy
c
      n2drec = (((i2-i1+1)*(j2-j1+1)+4095)/4096)*4096
      allocate (drec(n2drec))
      inquire(iolength=nrecl) drec
      deallocate (drec)
c
c     Full name of (mean, if needed) forcing files
c
      ffoc_path = '/scratch/projects/ome/hra_expt/UVDP/'
      mfoc_path = '/scratch/projects/ome/hra_expt/UVDP_CS/'
      if (flg.eq.0) then  
        dir_flnm = trim(ffoc_path)//flnm
      elseif (flg.eq.1) then
        dir_flnm = trim(mfoc_path)//flnm
      endif
c
c     Open forcing archives 
c
      write(lp,*) 'READ data from ',dir_flnm !flnm(1:l) 
      open(ni+1000, file=dir_flnm, access='direct',
     $     recl=nrecl, form='unformatted', status='old')
c
c --- Read instantaneous dp ('dpi')
c
      do k=1,kk
         read(ni+1000,rec=k)((dpio(i,j,k),i=i1,i2),j=j1,j2)
         if(k.eq.1.or.k.eq.25) then
            write(lp,*) itest,jtest,k,'dpio=',dpio(itest,jtest,k)
         endif
      enddo
c
c --- Read daily mean fileds
c
      nf=ni+1000
      if (timflg.eq.3) nf=nf+1
      do k=1,kk
         recsh=kk+nforf*(k-1)
         read(ni+1000,rec=recsh+1) (( uflx(i,j,k),i=i1,i2),j=j1,j2)
         read(ni+1000,rec=recsh+2) (( vflx(i,j,k),i=i1,i2),j=j1,j2)
         read(ni+1000,rec=recsh+3) ((  dpm(i,j,k),i=i1,i2),j=j1,j2)
      enddo                     !k
      close(ni+1000)
C
C --- READ mixed layer depth [m, MLD] to homogenize tracer above
C
      ! set the name (time) of MLD file according to uvdp file
      mldfnm = 'mld_offline.day_hr.a'
      l2 = len_trim(mldfnm)
      write(mldfnm(l2-7:l2-5),'(a3)') flnm(l-7:l-5) ! iday in tracers.f
      write(mldfnm(l2-3:l2-2),'(a2)') flnm(l-3:l-2) ! ihr
      mldfnm_ful = trim('/projects2/rsmas/ikamenkovich/yxl1496/')
     &         //'GSH_OUTPUT/MLD/'//trim(mldfnm)
      if(mxlflg.eq.1) then
         write(6,*) 'READ MLD from ',mldfnm_ful
         open(98,file=mldfnm_ful,access='direct',
     $        recl=nrecl,form='unformatted',status='old')
         read(98,rec=1) ((dpmixl(i,j),i=i1,i2),j=j1,j2) 
         close(98)
         write(6,*) itest,jtest,'mld=',dpmixl(itest,jtest)
      endif
      return
c
c --- unexpected end of file
 6    continue
        write (lp,*) '***** 2 unexpected end of archive file *****'
        write (lp,*) flnm(1:l-2)//'.b'
        call flush(lp)
        stop '(e-o-f)'
      end
c=====================================================================
c          SUBTOUTINE
c=====================================================================
      subroutine getfld(field, iunit, hminb,hmaxb)
      use mod_za ! HYCOM array I/O interface
c
      implicit none
c
      integer iunit
      real    field(idm,jdm), hminb,hmaxb
c
c --- read a single array
c
      integer i, j
      integer, dimension(idm,jdm) :: mask
      real    hmina,hmaxa
c
      call zaiord(field,mask,.false., hmina,hmaxa, iunit)
c
c      if     (abs(hmina-hminb).gt.abs(hminb)*1.e-4 .or.
c     &        abs(hmaxa-hmaxb).gt.abs(hmaxb)*1.e-4     ) then
c        write(lp,'(/ a / a,1p3e14.6 / a,1p3e14.6 /)')
c     &    'error - .a and .b files not consistent:',
c     &    '.a,.b min = ',hmina,hminb,hmina-hminb,
c     &    '.a,.b max = ',hmaxa,hmaxb,hmaxa-hmaxb
c        call flush(lp)
c        stop
c      endif
c
      do j= 1,jdm
        do i= 1,idm
          if     (field(i,j).gt.2.0**99) then
            field(i,j) = 0.0
          endif
        enddo
      enddo
      return
      end
