      program make_uvdp_HRA
!
!     The program reads the variables needed for the offline model
!     and saves them in the forcing files (uflx, vflx, dpm).
!
      implicit none
      integer idm,jdm, kdm, nrecl, i,j,k, index, nday,nyear,ia,ja, 
     $     n, ip1, ip2, jp1, jp2, idum, jdum, nbdy, lyear
      integer npad
      real, parameter :: onem=9806.0
      real, parameter :: lflg=1.0e10
      real fact
      character*1 yr_dir
      character*3 dayn
      character*4 yearn
      character*2 hourn
      character*100 grd_path,atl_path,sav_path
      integer ndsh,days,daye
c
c     FORCING FIELDS ARE 1573 BY 1073
c
      parameter(idum=1573,jdum=1073, kdm=30)
      real, dimension(idum,jdum):: depth, depthu, depthv, 
     $     var1, var2
      real, dimension(idum,jdum,kdm) :: dp1, dp2, dpu1, dpv1, dpu2, dpv2
c
c     extra points are needed for tracer advection
c
      real, allocatable :: pad(:) 
c
      npad=4096-mod(idum*jdum,4096)
      allocate (pad(npad))
      inquire(iolength=nrecl)var1,pad ! what's this for?????
c
      grd_path = '/nethome/ikamenkovich/hycom/GS_HR/'
      atl_path = '/projects2/rsmas/ikamenkovich/Atlantic_HR/'
      sav_path = '/projects2/rsmas/ikamenkovich/yxl1496/GSH_OUTPUT/'

c=====================================================================
c     Read topography (fid = 31) and apply land mask
c=====================================================================
      open(31,file=trim(grd_path)//'regional.depth.a', access='direct',
     $     recl=nrecl, form='unformatted', status='old',
     $     convert='big_endian')

      read(31,rec=1) depth
c     mask land
      do i=1,idum
         do j=1,jdum
            if (depth(i,j).gt.1.e10) then
               depth(i,j)=0.0
            endif
         enddo
      enddo
      close(31)
c=====================================================================
c     read data at current and next time (t & t+dt), fid = 12 & 13
c=====================================================================
      nyear = 8  
      ndsh = 133
      days = 1 + ndsh ! HYCOM used 366-day year, we want 365 days
c      days=233+ndsh
      daye = 2 + ndsh
      write(yearn,'(i4.4)') nyear ! '0008'

c     loop over each time step
      do nday=2*days,2*daye
c
c     ================ open data file at t 
c
         print*,'DAY',nday
         write(dayn,'(i3.3)') nday/2 - 366*(nyear-8)
         write(hourn,'(i2.2)') 12*(nday-2*(nday/2))
         open(12, file=trim(atl_path)//'DATA/016_archv.'//yearn//
     $        '_'//dayn//'_'//hourn//'.a', access='direct', recl=nrecl,
     $        convert='big_endian', form='unformatted', status='old')
         print*,'Data from year ',yearn,' day ',dayn,' hour ',hourn,
     $        ' to'

         if (nday+1.eq.367*2) then
            nyear = 9
            write(yearn,'(i4.4)') nyear
         endif
c
c     ================ open data at t+dt
c
         write(dayn,'(i3.3)') (nday+1)/2 - 366*(nyear-8)
         write(hourn,'(i2.2)') 12*(nday+1-2*((nday+1)/2))
         print*,'Data from year ',yearn,' day ',dayn,' hour ',hourn
         open(13,file=trim(atl_path)//'DATA/016_archv.'//yearn//
     $        '_'//dayn//'_'//hourn//'.a',access='direct',recl=nrecl,
     $        convert='big_endian',form='unformatted',status='old')

c     ============== UVDP file named as t+dt (relative time) (fid = 21)
         write(dayn,'(i3.3)') (nday+1)/2 - ndsh !-1 
         open(21,file=trim(sav_path)//'UVDP/uvdp_offline.'//dayn//'_'//
     $        hourn//'.a',access='direct',recl=nrecl,
     $        convert='big_endian',form='unformatted',status='unknown')

c     ========== Save the mean layer thck ('dpm') averaged over t & t+dt 
c                to UVDP file
c     Note that the instant dp at t+dt is also saved
         do k = 1,kdm
            read(12,rec=11+(k-1)*5+3) var1
            read(13,rec=11+(k-1)*5+3) var2
            do i=1,idum
               do j=1,jdum
                  fact=depth(i,j)
                  if (fact.eq.0.0) then
                     var1(i,j)=0.0
                     var2(i,j)=0.0
                  else
                     var1(i,j)=var1(i,j)/onem
                     var2(i,j)=var2(i,j)/onem
                  endif
                  dp1(i,j,k)=var1(i,j)
                  dp2(i,j,k)=var2(i,j)
               enddo
            enddo
            write(21,rec=k) var2
            write(21,rec=kdm+3*(k-1)+3) (var1+var2)/2.0
         enddo
c
c     ========== Calc layer thickness at UV-grid (dpu1,dpv1;dpu2,dpv2) 
c         for the conversion of velocities to mass fluxes
c
         do j=1,jdum
            do i=1,idum
               ia=max(1,i-1)
               depthu(i,j)=min(depth(i,j),depth(ia,j))
               ja=max(1,j-1)
               depthv(i,j)=min(depth(i,j),depth(i,ja))
            enddo
         enddo  
         ! t
         if (nday.eq.2*days) then
            call calc_dpuv(idum,jdum,dp1,dpu1,dpv1,depthu,depthv)
         else
            dpu1=dpu2
            dpv1=dpv2
         endif
         ! t+dt
         call calc_dpuv(idum,jdum,dp2,dpu2,dpv2,depthu,depthv)
         print*,'DPU and DPV are calculated'
c
c     ================ Calc u*dp & v*dp (not dpm) at t & t+dt and save
c      the time-mean value to UVDP file
c      if layer vanishes, Hycom sets vel(k)=vel(k-1). Set those values
c      to zero here.
c
         do k=1,kdm
c           do u      
            read(12,rec=11+(k-1)*5+1) var1
            read(13,rec=11+(k-1)*5+1) var2
            do i=1,idum
               do j=1,jdum
                  fact=depth(max(i-1,1),j)*depth(i,j)
                  if (fact.eq.0.0) var1(i,j)=0.0
                  if (fact.eq.0.0) var2(i,j)=0.0
                  if (var1(i,j).gt.lflg) var1(i,j)=0.0
                  if (var2(i,j).gt.lflg) var2(i,j)=0.0
               enddo
            enddo
            do i=1,idum
               do j=1,jdum
                  var1(i,j)=var1(i,j)*dpu1(i,j,k)
                  var2(i,j)=var2(i,j)*dpu2(i,j,k)
               enddo
            enddo
            write(21,rec=kdm+3*(k-1)+1) (var1+var2)/2.0
c
c           do v 
c     
            read(12,rec=11+(k-1)*5+2) var1
            read(13,rec=11+(k-1)*5+2) var2
            do i=1,idum
               do j=1,jdum
                  fact=depth(i,max(j-1,1))*depth(i,j)
                  if (fact.eq.0.0) var1(i,j)=0.0
                  if (fact.eq.0.0) var2(i,j)=0.0
                  if (var1(i,j).gt.lflg) var1(i,j)=0.0
                  if (var2(i,j).gt.lflg) var2(i,j)=0.0
               enddo
            enddo
            do i=1,idum
               do j=1,jdum
                  var1(i,j)=var1(i,j)*dpv1(i,j,k)
                  var2(i,j)=var2(i,j)*dpv2(i,j,k)
               enddo
            enddo
            write(21,rec=kdm+3*(k-1)+2) (var1+var2)/2.0
         enddo                  ! k
         close(12)
         close(13)
         close(21)
         write(6,*) 'UVDP saved: '//trim(sav_path)//'UVDP/uvdp_offline.'
     $   //dayn//'_'//hourn//'.a'
      enddo                     ! nday
c
      stop
      end

c=====================================================================
c       subroutine 'calc_dpuv': calc dp at UV-grid
c=====================================================================

      subroutine calc_dpuv(idum,jdum,dpm,dpu,dpv,depthu,depthv)
      implicit none
      integer, parameter :: kdm=30
      integer :: idum, jdum, i,j,k 
      real :: duk,dukm1,dvkm1,dvk,depthi(0:1,0:1,0:kdm)
      real, dimension(idum,jdum) :: depthu,depthv
      real, dimension(idum,jdum,kdm) :: dpm, dpu, dpv
c
      do j = 2,jdum
         do i = 2,idum
            duk           = 0.0
            dvk           = 0.0
            depthi(:,:,0) = 0.0
            do k= 1,kdm
               depthi(1,1,k) =
     &              depthi(1,1,k-1) + dpm(i  ,j,k)
               depthi(0,1,k) =
     &              depthi(0,1,k-1) + dpm(i-1,j,k)
               depthi(1,0,k) =
     &              depthi(1,0,k-1) + dpm(i,j-1,k)
               dukm1 = duk
               duk   = min( depthu(i,j),
     &              0.5*(depthi(0,1,k) + depthi(1,1,k)) )
               dvkm1 = dvk
               dvk   = min( depthv(i,j),
     &              0.5*(depthi(1,0,k) + depthi(1,1,k)) )
               dpu(i,j,k) = max( 0.0, duk-dukm1 )
               dpv(i,j,k) = max( 0.0, dvk-dvkm1 )
            enddo               !k
         enddo                  !i
      enddo                     !j
c     
      return
      end
