      program make_mld_HRA
!
!     The program reads and saves MLD for the use of offline model 
!
      implicit none
      integer idm,jdm, kdm, nrecl, i,j,k, index, nday,nyear,
     $     n
      integer npad
      character*3 dayn
      character*4 yearn
      character*2 hourn
      character*100 grd_path,atl_path,sav_path
      integer ndsh,days,daye
      character*200 flnm
      logical  ifexist
c
c     FORCING FIELDS ARE 1573 BY 1073
c
      parameter(idm=1573,jdm=1073, kdm=30)
      real, dimension(idm,jdm):: depth,  var
c
c     extra points are needed for tracer advection
c
      real, allocatable :: pad(:) 
c
      npad=4096-mod(idm*jdm,4096)
      allocate (pad(npad))
      inquire(iolength=nrecl)var,pad ! what's this for?????
c
      grd_path = '/nethome/ikamenkovich/hycom/GS_HR/'
      atl_path = '/projects2/rsmas/ikamenkovich/Atlantic_HR/'
      sav_path = '/projects2/rsmas/ikamenkovich/yxl1496/GSH_OUTPUT/MLD/'
      write(6,*) 'nrecl = ', nrecl
c=====================================================================
c     Read topography (fid = 31) and apply land mask
c=====================================================================
      open(31,file=trim(grd_path)//'regional.depth.a', access='direct',
     $     recl=nrecl, form='unformatted', status='old',
     $     convert='big_endian')
c
      read(31,rec=1) depth
      close(31)
c=====================================================================
c     read vel & thck data at t & t+dt, fid = 12 & 13
c=====================================================================
      nyear = 8  
      ndsh = 133
      days = 2 + ndsh ! HYCOM used 366-day year, we want 365 days
c      days=233+ndsh
      daye = 366 + ndsh
      write(yearn,'(i4.4)') nyear ! '0008'

c     loop over each time step
      do nday=2*days,2*daye
c
c     ================ open data at t+dt (13)
c
         if (nday+1.eq.367*2) then
            nyear = 9
            write(yearn,'(i4.4)') nyear
         endif

         write(dayn,'(i3.3)') (nday+1)/2 - 366*(nyear-8)
         write(hourn,'(i2.2)') 12*(nday+1-2*((nday+1)/2))
         flnm = trim(atl_path)//'DATA/016_archv.'//yearn//
     $        '_'//dayn//'_'//hourn//'.a'
         ! check exist
         INQUIRE(FILE=flnm, EXIST=ifexist)
         if (ifexist) then
            write(6,*) 'READ from: ', flnm
            open(13,file=flnm,access='direct',recl=nrecl,
     $           convert='big_endian',form='unformatted',status='old')
         else
            write(6,*) 'NOT exist, skip: ', flnm    
            write(6,*) ' ----------- '    
            cycle 
         endif

c     ============== SAVE file named as t+dt (relative time) (fid = 21)
         write(dayn,'(i3.3)') (nday+1)/2 - ndsh -1 
         open(21,file=trim(sav_path)//'mld_offline.'//dayn//'_'//
     $        hourn//'.a',access='direct',recl=nrecl,
     $        convert='big_endian',form='unformatted',status='unknown')

c     ========== Process and save MLD
         read(13,rec=6) var ! mld at t+dt
         do i = 1,idm
            do j = 1,jdm    
               if (depth(i,j).gt.1.e10) then
                  var(i,j) = 0.0
               endif 
               var(i,j) = var(i,j) / 9806.0 ! in [m] 
            enddo
         enddo
         write(6,'(a,f20.10,f20.10)') 'mld', var(400,500), var(800,800)
         write(21,rec=1)  ( (var(i,j), i=1,idm), j=1,jdm )  ! write into .a file
         ! same with "write(...) var"
C 
         close(13)
         close(21)
         write(6,*) 'MLD saved: '//trim(sav_path)//'mld_offline.'
     $   //dayn//'_'//hourn//'.a'
         write(6,*) ' ----------- '
      enddo                     ! nday
c
      stop
      end



