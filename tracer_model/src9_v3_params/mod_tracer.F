      module mod_tracer
      use mod_xc    ! HYCOM communication interface
      use mod_za
c
      implicit none
      public
      character*120 flnm_ai,flnm_ao,flnmi_ai,flnmi_ao,flnm_o
      real, save, allocatable, dimension (:,:,:,:) :: tracer,dp,ktensor
      real, save, allocatable, dimension (:,:,:) :: 
     &     p,dpm,dpu,dpv,dpii,dpio,hordiv,uflx,vflx,diaflx,difs,ghats
     &     , tracer_prev !------------- added (used in tracers.f)------------!
     &     , dpcons
      real, save, allocatable, dimension (:,:) :: dpmixl,
     &     depths,depthu,depthv,
     &     scpx,scpy,scux,scuy,scvx,scvy,scu2,scv2,scp2,scp2i,
     &     plon,plat,ulon,ulat,vlon,vlat, 
     &     tracfor,util1,util2
      real, save :: archfq,archin,onecm,tencm,thbase,
     &     deltfl,fldepm,tbvar,tdecri,dtturb,uturb0,
     &     dudx,dvdy,radian,pi,diagfq, vertdiv, delt1, flxdiv0
      integer, save, allocatable, dimension (:,:)   ::
     &     ip,iu,iv,ifp,ilp,jfp,jlp,ifu,ilu,jfu,jlu,ifv,ilv,jfv,jlv,
     $     ipim1, ipip1, ipjm1, ipjp1, klist
      integer, save, allocatable, dimension (:)  :: 
     &     isp,jsp,isu,jsu,isv,jsv
      integer, save :: ktr,kdm,ii,jj,kk,iar,jar,narch,ntracr,difflg,
     &     tracin,trcflg,advflg,itest,jtest,n,m,nintrp,mxlflg,timflg,
     $     ntim,ghtflg,forflg,mapflg,lperiod,larctic,jblk,nforf
      integer i,j,k
      real,    parameter   :: huge=1.0e20
      real,    parameter   :: onem=9806.0
      real,    parameter   :: thref=1.0e-3
      real,    parameter   :: g=9.806
      real,    parameter   :: temdf2=0.02    ! horiz diffus (orig 0.02): [m/s]
      real,    parameter   :: onezm=9806.e-20
      real,    parameter   :: onemm = 0.001*onem

      contains

c=====================================================================
c          SUBTOUTINE: linspace
c=====================================================================
      subroutine linspace(x, xi, xf, len)
      implicit none
      real, dimension(:), intent(out) :: x(:)
      real, intent(in)  :: xi, xf
      real              :: dx
      integer           :: len, i
      
      if (len <= 0) return
      if (len == 1) then
         x(1) = xi
         return
      endif
      ! delta x
      dx = (xf - xi) / (len - 1)
      x(1:len) = [(xi + (i - 1) * dx, i = 1,len)]
      end subroutine linspace

c=====================================================================
c          SUBTOUTINE: allocate
c=====================================================================
      subroutine tracer_alloc
      implicit none
c
      allocate( tracer(1-nbdy:ii+nbdy,1-nbdy:jj+nbdy,kk,ntracr))
      allocate(   uflx(1-nbdy:ii+nbdy,1-nbdy:jj+nbdy,kk  ) )
      allocate(   vflx(1-nbdy:ii+nbdy,1-nbdy:jj+nbdy,kk  ) )
      allocate( diaflx(1-nbdy:ii+nbdy,1-nbdy:jj+nbdy,kk  ) )
      allocate( hordiv(1-nbdy:ii+nbdy,1-nbdy:jj+nbdy,kk  ) )
      allocate(   difs(1-nbdy:ii+nbdy,1-nbdy:jj+nbdy,kk  ) )
C       
C     ----------- ADDED FOR 'fct3d' OUTPUT -------------
C    
      allocate(tracer_prev(1-nbdy:ii+nbdy,1-nbdy:jj+nbdy,kk  ) )
      allocate(dpcons(1-nbdy:ii+nbdy,1-nbdy:jj+nbdy,kk  ))
C       
      if (ghtflg.eq.1) then
         allocate(ghats(1-nbdy:ii+nbdy,1-nbdy:jj+nbdy,kk  ) )
      endif
      if (difflg.eq.1) then
c     ISOTROPIC (SCALAR) diffusivity
         allocate(ktensor(1-nbdy:ii+nbdy,1-nbdy:jj+nbdy,kk,1) )
      elseif (difflg.eq.2) then
c     ANISOTROPIC (2D TENSOR) diffusivity
         allocate(ktensor(1-nbdy:ii+nbdy,1-nbdy:jj+nbdy,kk,4) )
      elseif (difflg.eq.3) then
c     ISOTROPIC diffusivity & advective velocity
         allocate(ktensor(1-nbdy:ii+nbdy,1-nbdy:jj+nbdy,kk,3) )
      endif
c      allocate(      p(1-nbdy:ii+nbdy,1-nbdy:jj+nbdy,kk+1) )
      allocate(    dpm(1-nbdy:ii+nbdy,1-nbdy:jj+nbdy,kk  ) )
      allocate(     dp(1-nbdy:ii+nbdy,1-nbdy:jj+nbdy,kk,2) )
      allocate(   dpii(1-nbdy:ii+nbdy,1-nbdy:jj+nbdy,kk  ) )
      allocate(   dpio(1-nbdy:ii+nbdy,1-nbdy:jj+nbdy,kk  ) )
      allocate( dpmixl(1-nbdy:ii+nbdy,1-nbdy:jj+nbdy) )
      allocate( depths(1-nbdy:ii+nbdy,1-nbdy:jj+nbdy) )
      allocate( depthu(1-nbdy:ii+nbdy,1-nbdy:jj+nbdy) )
      allocate( depthv(1-nbdy:ii+nbdy,1-nbdy:jj+nbdy) )
      allocate(   plon(1-nbdy:ii+nbdy,1-nbdy:jj+nbdy) )
      allocate(   plat(1-nbdy:ii+nbdy,1-nbdy:jj+nbdy) )
      allocate(   ulon(1-nbdy:ii+nbdy,1-nbdy:jj+nbdy) )
      allocate(   ulat(1-nbdy:ii+nbdy,1-nbdy:jj+nbdy) )
      allocate(   vlon(1-nbdy:ii+nbdy,1-nbdy:jj+nbdy) )
      allocate(   vlat(1-nbdy:ii+nbdy,1-nbdy:jj+nbdy) )
      allocate(   scpx(1-nbdy:ii+nbdy,1-nbdy:jj+nbdy) )
      allocate(   scpy(1-nbdy:ii+nbdy,1-nbdy:jj+nbdy) )
      allocate(   scux(1-nbdy:ii+nbdy,1-nbdy:jj+nbdy) )
      allocate(   scuy(1-nbdy:ii+nbdy,1-nbdy:jj+nbdy) )
      allocate(   scvx(1-nbdy:ii+nbdy,1-nbdy:jj+nbdy) )
      allocate(   scvy(1-nbdy:ii+nbdy,1-nbdy:jj+nbdy) )
      allocate(   scu2(1-nbdy:ii+nbdy,1-nbdy:jj+nbdy) )
      allocate(   scv2(1-nbdy:ii+nbdy,1-nbdy:jj+nbdy) )
      allocate(   scp2(1-nbdy:ii+nbdy,1-nbdy:jj+nbdy) )
      allocate(   scp2i(1-nbdy:ii+nbdy,1-nbdy:jj+nbdy))
      allocate(   util1(1-nbdy:ii+nbdy,1-nbdy:jj+nbdy))
      allocate(   util2(1-nbdy:ii+nbdy,1-nbdy:jj+nbdy))
C       
C     ----------- ADDED FOR 'fct3d' OUTPUT (used in tracers.f) -------------
C 
      if(forflg) then
      allocate( tracfor(1-nbdy:ii+nbdy,1-nbdy:jj+nbdy) )
      endif
c
      allocate(  klist(1-nbdy:ii+nbdy,1-nbdy:jj+nbdy) )
      allocate(     ip(1-nbdy:ii+nbdy,1-nbdy:jj+nbdy) )
      allocate(     iu(1-nbdy:ii+nbdy,1-nbdy:jj+nbdy) )
      allocate(     iv(1-nbdy:ii+nbdy,1-nbdy:jj+nbdy) )
      allocate( ipim1(1-nbdy+1:ii+nbdy-1,1-nbdy+1:jj+nbdy-1) )
      allocate( ipip1(1-nbdy+1:ii+nbdy-1,1-nbdy+1:jj+nbdy-1) )
      allocate( ipjm1(1-nbdy+1:ii+nbdy-1,1-nbdy+1:jj+nbdy-1) )
      allocate( ipjp1(1-nbdy+1:ii+nbdy-1,1-nbdy+1:jj+nbdy-1) )
      allocate(    ifp(1-nbdy:jj+nbdy,ms) )
      allocate(    ilp(1-nbdy:jj+nbdy,ms) )
      allocate(    ifu(1-nbdy:jj+nbdy,ms) )
      allocate(    ilu(1-nbdy:jj+nbdy,ms) )
      allocate(    ifv(1-nbdy:jj+nbdy,ms) )
      allocate(    ilv(1-nbdy:jj+nbdy,ms) )
      allocate(    jfp(1-nbdy:ii+nbdy,ms) )
      allocate(    jlp(1-nbdy:ii+nbdy,ms) )
      allocate(    jfu(1-nbdy:ii+nbdy,ms) )
      allocate(    jlu(1-nbdy:ii+nbdy,ms) )
      allocate(    jfv(1-nbdy:ii+nbdy,ms) )
      allocate(    jlv(1-nbdy:ii+nbdy,ms) )
c
      allocate(    isp(1-nbdy:jj+nbdy)    )
      allocate(    isu(1-nbdy:jj+nbdy)    )
      allocate(    isv(1-nbdy:jj+nbdy)    )
      allocate(    jsp(1-nbdy:ii+nbdy)    )
      allocate(    jsu(1-nbdy:ii+nbdy)    )
      allocate(    jsv(1-nbdy:ii+nbdy)    )

      end subroutine tracer_alloc

c=====================================================================
c          SUBTOUTINE: tracer initialization
c=====================================================================
      subroutine tracer_init(timein)
      real timein,timerst,thet,hminb,hmaxb,aktr
      integer nti,layer,mon,ntr
      character cline*80 

c     params for analytic form of intial tracer
      real, parameter   :: pi = 3.14159265358979, c = 1.
      real, parameter   :: iwav1 = 1., iwav2 = 2.
      real, dimension (1-nbdy:ii+nbdy) :: x_lr
      real, dimension (1-nbdy:jj+nbdy) :: y_lr

      data nti/21/

      if(tracin.eq.0) then
c --- initialize tracer manually
         do ktr=1,ntracr
            write(lp,*) '-- INITIALIZE tracer-', ktr
!$OMP   PARALLEL DO PRIVATE(i,j,k)
            do k=1,kk
               do j=1-nbdy,jj+nbdy
                  do i=1-nbdy,ii+nbdy                     
                     if (ktr.eq.1) then
                        !  meridional-propagation wave No.1
                        call linspace(x_lr, 0., 0., ii+2*nbdy)
                        call linspace(y_lr, 0., 1., jj+2*nbdy)
                        tracer(i,j,k,ktr) = c*cos(2*pi*iwav1*x_lr(i) + 
     &                  2*pi*iwav1*y_lr(j)) + 2.*c
                     elseif (ktr.eq.2) then
                        !  meridional wave No.2
                        call linspace(x_lr, 0., 0., ii+2*nbdy)
                        call linspace(y_lr, 0., 1., jj+2*nbdy)
                        tracer(i,j,k,ktr) = c*cos(2*pi*iwav2*x_lr(i) + 
     &                  2*pi*iwav2*y_lr(j)) + 2.*c
                     elseif (ktr.eq.3) then
                        !  zonal wave No.1
                        call linspace(x_lr, 0., 1., ii+2*nbdy)
                        call linspace(y_lr, 0., 0., jj+2*nbdy)
                        tracer(i,j,k,ktr) = c*cos(2*pi*iwav1*x_lr(i) + 
     &                  2*pi*iwav1*y_lr(j)) + 2.*c
                     elseif (ktr.eq.4) then
                        !  zonal wave No.2
                        call linspace(x_lr, 0., 1., ii+2*nbdy)
                        call linspace(y_lr, 0., 0., jj+2*nbdy)
                        tracer(i,j,k,ktr) = c*cos(2*pi*iwav2*x_lr(i) + 
     &                  2*pi*iwav2*y_lr(j)) + 2.*c
                     elseif (ktr.eq.5) then
                        ! diagnoal (lowL-upR) wave No.1
                        call linspace(x_lr, 0., 1., ii+2*nbdy)
                        call linspace(y_lr, 0., 1., jj+2*nbdy)
                        tracer(i,j,k,ktr) = c*cos(2*pi*iwav1*x_lr(i) + 
     &                  2*pi*iwav1*y_lr(j)) + 2.*c
                     elseif (ktr.eq.6) then
                        ! diagnoal (lowL-upR) wave No.2
                        call linspace(x_lr, 0., 1., ii+2*nbdy)
                        call linspace(y_lr, 0., 1., jj+2*nbdy)
                        tracer(i,j,k,ktr) = c*cos(2*pi*iwav2*x_lr(i) + 
     &                  2*pi*iwav2*y_lr(j)) + 2.*c
                     elseif (ktr.eq.7) then
                        ! diagnoal (upL-lowR) wave No.1
                        call linspace(x_lr, 1., 0., ii+2*nbdy)
                        call linspace(y_lr, 0., 1., jj+2*nbdy)
                        tracer(i,j,k,ktr) = c*cos(2*pi*iwav1*x_lr(i) + 
     &                  2*pi*iwav1*y_lr(j)) + 2.*c
                     elseif (ktr.eq.8) then
                        ! diagnoal (upL-lowR) wave No.2
                        call linspace(x_lr, 1., 0., ii+2*nbdy)
                        call linspace(y_lr, 0., 1., jj+2*nbdy)
                        tracer(i,j,k,ktr) = c*cos(2*pi*iwav2*x_lr(i) + 
     &                  2*pi*iwav2*y_lr(j)) + 2.*c
                     elseif (ktr.eq.9) then
                        ! two waves summed
                        call linspace(x_lr, 0., 1., ii+2*nbdy)
                        call linspace(y_lr, 0., 1., jj+2*nbdy)
                        tracer(i,j,k,ktr) = .5 * c * ( cos(2*pi*iwav2*
     &                  x_lr(i)) + sin(2*pi*iwav2*y_lr(j)) ) + 2.*c
                     elseif (ktr.eq.10) then
                        ! 1/4 of a paraboloid
                        call linspace(x_lr, -1., -0.6, ii+2*nbdy)
                        call linspace(y_lr, -1., -0.6, jj+2*nbdy)
                        tracer(i,j,k,ktr) = c * ( x_lr(i)**2. / 9. + 
     &                  y_lr(j)**2. / 16. ) * 25.                    
                     endif
                  enddo
               enddo
            enddo
!$OMP END PARALLEL DO
         enddo
      elseif(tracin.eq.1) then
c --- read 'restart_tracr_in.[a,b]' 
         open (unit=nti,file='restart_tracer_in.b',form='formatted',
     &        status='old',action='read')
         call zaiopf('restart_tracer_in.a','old', nti)
         do ktr=1,ntracr
            read(nti,'(a)') cline
            i = index(cline,'=')
            read (cline(i+1:),*)  timerst,layer,hminb,hmaxb
            if(timein.ne.timerst)then
               write(6,*)'error, time NOT EQUAL time of restart file:'
               write(6,*) timein, timerst
               call flush(lp)
               stop
            endif
            do k=1,kk   
c     call getfld(tracer(1,1,k,1,ktr),nti,hminb,hmaxb)
               ntr=(ktr-1)*kk+k
               read(nti+1000,rec=ntr) ((tracer(i,j,k,ktr),i=1-nbdy,
     $              ii+nbdy),j=1-nbdy,jj+nbdy)
               write(lp,'("input  ",a," into ",a,i3)')
     $              cline(1:8),'layer ',k
               call flush(lp)
            enddo
         enddo
         close(nti)
         call zaiocl(nti)
      elseif(tracin.eq.2) then
c --- read 'relax_tracer.[a,b]' 
         open (unit=nti,file='relax.tracer.b',form='formatted',
     &        status='old',action='read')
         call zaiopf('relax.tracer.a','old', nti)
         do k=1,5
            read(nti,'(a)') cline
         enddo
         do ktr=1,ntracr
            do k=1,kk
               read(nti,'(a)') cline
               i = index(cline,'=')
               read (cline(i+1:),*)  mon,layer,aktr,hminb,hmaxb
               call getfld(tracer(1,1,k,ktr),nti,hminb,hmaxb)
               write(lp,'("input  ",a,i3,f5.2," into ",a,i3)') 
     .              cline(1:8),mon,aktr,'lay',k
               call flush(lp)
               call tracer_buffer(tracer(1,1,k,ktr))
            enddo
         enddo
      endif                     !types of initialization
c
      return
      end subroutine tracer_init

c=====================================================================
c          SUBTOUTINE: tracer initialization for SKILL TESTING 
c=====================================================================
      subroutine tracer_init_ctest(timein)
      real timein,timerst,thet,hminb,hmaxb,aktr
      integer nti,layer,mon,ntr
      character cline*80 

c     params for analytic form of intial tracer
      real, parameter   :: pi = 3.14159265358979, c = 1.
      real, parameter   :: iwav1 = 1., iwav2 = 2.
      real, dimension (1-nbdy:ii+nbdy) :: x_lr
      real, dimension (1-nbdy:jj+nbdy) :: y_lr
C       
C     param for the circular patch tracer
C     Note (x0,y0) here corresponds to (x0+nbdy,y0+nbdy) in [1:IDM,1:JDM]
      real, parameter   :: r2 = 75.**2., x0 = 898., y0 = 398.
C     conc at center, edge and outside
      real, parameter   :: c_cen = 3., c_edg = 1., c_out = 0.
      real :: dist2
C 
      data nti/21/
C 
C     New set of tracers used here
C 
      if(tracin.eq.0) then
c --- initialize tracer manually
         do ktr = 1,ntracr
            write(lp,*) '-- INITIALIZE the TEST tracer-', ktr
!$OMP   PARALLEL DO PRIVATE(i,j,k)
            do k=1,kk
               do j=1-nbdy,jj+nbdy
                  do i=1-nbdy,ii+nbdy                     
                     if (ktr.eq.0) then
                        ! TEST-1. a tracer blob
                        ! dist2 from any point to circ center
                        dist2 = (i - x0)**2. + (j - y0)**2.
                        ! outside or winthin/on edge
                        if (dist2.gt.r2) then
                           tracer(i,j,k,ktr) = 0.
                        elseif (dist2.le.r2) then
                           tracer(i,j,k,ktr) = c_cen + 
     &                     (c_edg - c_cen) * dist2 / r2;
                        endif
                        ! initialize tracer only in certain layer(s)
                        if (k.lt.20 .OR. k.gt.20) then
                           !tracer(i,j,k,ktr) = 0.0
                        endif                         

                     elseif (ktr.eq.0) then
                        ! TEST-2. CONSTANT meridional gradient 
                        call linspace(x_lr, 0., 0., ii+2*nbdy)
                        call linspace(y_lr, 4., 0., jj+2*nbdy)
                        tracer(i,j,k,ktr) = y_lr(j) + c

                     elseif (ktr.eq.0) then
                        ! TEST-3. A tracer used in KL 
                        !    meridional-propagation wave No.1
                        call linspace(x_lr, 0., 0., ii+2*nbdy)
                        call linspace(y_lr, 0., 1., jj+2*nbdy)
                        tracer(i,j,k,ktr) = c*cos(2*pi*iwav1*x_lr(i) + 
     &                  2*pi*iwav1*y_lr(j)) + 2.*c

                     elseif (ktr.eq.1) then
                        ! TEST-4. A tracer used in KL 
                        !    meridional-propagation wave No.2
                        call linspace(x_lr, 0., 0., ii+2*nbdy)
                        call linspace(y_lr, 0., 1., jj+2*nbdy)
                        tracer(i,j,k,ktr) = c*cos(2*pi*iwav2*x_lr(i) + 
     &                  2*pi*iwav2*y_lr(j)) + 2.*c          
                     endif
                  enddo
               enddo
            enddo
!$OMP END PARALLEL DO
         enddo
      else 
         stop 'error'
      endif
      end subroutine tracer_init_ctest

c====================================================================
c          SUBTOUTINE: calc total amount of c-mass of snapshot
c=====================================================================
      subroutine tracer_mass(fld,klay,cmass)
      implicit none
      real,dimension(1-nbdy:ii+nbdy,1-nbdy:jj+nbdy,kk),intent(in) :: fld
      integer,intent(in) :: klay
      real :: cmassj(1-nbdy:jj+nbdy)
      real,intent(out) :: cmass
C
      cmass = 0.
C       
      if (klay.gt.0) then
         do j = 1, jj
            cmassj(j) = 0.
            do i = 1, ii
               if (ip(i,j).ne.0) cmassj(j) = cmassj(j) + fld(i,j,klay)*
     &              dp(i,j,klay,2)*scp2(i,j)
            enddo 
            cmass = cmass + cmassj(j)          
         enddo     

      elseif (klay.eq.0) then
         do k = 1, kk
!$OMP PARALLEL DO SHARED(k)
            do j = 1, jj
               cmassj(j) = 0.
               do i = 1, ii
                  if (ip(i,j).ne.0) cmassj(j) = cmassj(j) + fld(i,j,k)*
     &               dp(i,j,k,2)*scp2(i,j)
               enddo               
            enddo                  
!$OMP END PARALLEL DO
            do j = 1, jj
               cmass = cmass + cmassj(j)
            enddo               
         enddo                     !k

      endif

      return
      end subroutine tracer_mass

c====================================================================
c          SUBTOUTINE: surface BC
c=====================================================================
      subroutine tracer_sbc(str,time_tot)
      implicit none
      integer,parameter :: its=1-nbdy,ite=1475+nbdy,
     $     jts=1181-nbdy,jte=1950+nbdy
      real, parameter :: tts=0.0,tte=365.0
      real str(1-nbdy:ii+nbdy,1-nbdy:jj+nbdy),dpm0,dpmmi,threfi,tracsfl
      double precision time_tot
c         
      if (trcflg.eq.1) then
c    
c     Employ a modified idealized age tracer (in the limited domain) 
c     
!$OMP PARALLEL DO PRIVATE(i,j)
         do j=jts,jte
            do i=its,ite
               if (ip(i,j).ne.0) str(i,j) = time_tot
            enddo
         enddo
!$OMP END PARALLEL DO
      elseif (trcflg.eq.2) then
c
c     Employ a pulse tracer during a limited time and in the limited domain
c                                                   
!$OMP PARALLEL DO PRIVATE(i,j)
         do i=1-nbdy,ii+nbdy
            do j=1-nbdy,jj+nbdy
               if (ip(i,j).ne.0) str(i,j)=0.0
            enddo
         enddo
!$OMP END PARALLEL DO
         if (time_tot.ge.tts .and. time_tot.le.tte) then
!$OMP PARALLEL DO PRIVATE(i,j)
            do i=its,ite
               do j=jts,jte
                  if (ip(i,j).ne.0) str(i,j)=1.0
               enddo
            enddo
!$OMP END PARALLEL DO
         endif
      elseif (trcflg.gt.10) then
c
c     Prescribe or calculate fluxes
c
!$OMP PARALLEL DO PRIVATE(dpm0,dpmmi,threfi,tracsfl)   
         do j=1-nbdy,jj+nbdy
            do i=1-nbdy,ii+nbdy
               if (trcflg.eq.11) then
c
c     Surface flux is prescribed
c
                  tracsfl = 0.0
               endif
               dpm0  =max(onemm,dpm(i,j,1))
               dpmmi =1.0/dpm0
               threfi=1.0/thref
               str(i,j) = str(i,j) + tracsfl*delt1*g*dpmmi*threfi   
            enddo
         enddo
!$OMP END PARALLEL DO
      endif
c
      return
      end subroutine tracer_sbc

c=====================================================================
c          SUBTOUTINE: lateral BC
c=====================================================================
      subroutine tracer_lbc(fld,time_tot,klay)
      implicit none
      integer i,j,klay
      integer, parameter :: jss=1-nbdy,jse=5,iss=1-nbdy,ise=1475+nbdy,
     $     jps=jss,jpe=jse,ips=iss,ipe=ise, kps=1,kpe=30
      real, parameter :: tps=0.0,tpe=365.0
      real fld(1-nbdy:ii+nbdy,1-nbdy:jj+nbdy)
      double precision time_tot
c         
      if (trcflg.eq.3) then
c
c     Employ a pulse tracer during a limited time and in the limited domain. 
c     Set the tracer to zero on the entire vertical surface
c 
!$OMP PARALLEL DO PRIVATE(i,j)
         do i=iss,ise
            do j=jss,jse
               fld(i,j) = 0.0
            enddo
         enddo
!$OMP END PARALLEL DO 
c
c     Introduce a pulse in a segment of the lateral surface and during 
c     a limited time
c 
         if (klay.ge.kps .and. klay.le.kpe) then
            if (time_tot.ge.tps .and. time_tot.le.tpe) then
!$OMP PARALLEL DO PRIVATE(i,j)
               do i=ips,ipe
                  do j=jps,jpe
                     if (ip(i,j).ne.0) fld(i,j) = 1.0
                  enddo
               enddo
!$OMP END PARALLEL DO
            endif
         endif
      endif
c
      return
      end subroutine tracer_lbc
c
c=====================================================================
c         Spatially smooth the tracer field (same as 'coarsen')
c===================================================================== 
      subroutine tracer_smooth(fld, nfili, nfilj, klay)
c
c --- Smooth the 2D field using a boxcar filter
c
      implicit none
      integer :: nfili, nfilj
      integer :: is, ie, js, je, ic, jc, i, j
      real, dimension(1-nbdy:ii+nbdy,1-nbdy:jj+nbdy) :: fld, flds, arsm
      real, dimension(1-nbdy:jj+nbdy) :: fsum_oj, fsum_nj
      real :: fsum_o, fsum_n, fldij, areaij
      real, parameter :: onemu = 1.e-6
      integer,intent(in) :: klay
c
c     Land or vanishing layer points are masked by 0's
c
C !OMP PARALLEL DO PRIVATE(i,j)
C       do j = 1-nbdy, jj+nbdy
C          do i = 1-nbdy, ii+nbdy
C             fld(i,j) = fld(i,j) !* scp2(i,j) ! * max(dpm(i,j,klay),onemu) 
C          enddo
C       enddo
C !OMP END PARALLEL DO

c
c     Calculate sums over j-direction at each point
c
!OMP PARALLEL DO PRIVATE(i,j,jc,js,je)
      do j = 1-nbdy, jj+nbdy
         fsum_oj(j) = 0.0
         js = max(1-nbdy, j-nfilj)
         je = min(jj+nbdy,j+nfilj)
         do i = 1-nbdy, ii+nbdy 
            ! tracermass before smooth 
            fsum_oj(j) = fsum_oj(j) + fld(i,j) * scp2(i,j)
            !
            flds(i,j) = 0.0
            arsm(i,j) = 0.0
            ! sum over current metidional segment, ignoring vanishing 
            ! or land points in the sum. But note these points can have 
            ! non-zero sums
             do jc = js, je
                if (dpm(i,jc,klay).gt.onemu .AND. fld(i,jc).gt.0.0) then
                   flds(i,j) = flds(i,j) + fld(i,jc) * scp2(i,jc)
                   arsm(i,j) = arsm(i,j) + scp2(i,jc) 
                endif
             enddo
         enddo
      enddo
!OMP END PARALLEL DO
c
c     Sum j-sums in the i-direction and divide by area
c
!OMP PARALLEL DO PRIVATE(i,j,ic,is,ie)
      do j = 1-nbdy, jj+nbdy
         fsum_nj(j) = 0.0
         do i = 1-nbdy, ii+nbdy
            is = max(1-nbdy, i-nfili)
            ie = min(ii+nbdy,i+nfili)
            ! sum the y-sum (flds, arsm) within the box. The sum at the
            ! vanishing layer or land points are 0's 
            fldij = 0.0
            areaij = 0.0
            ! 
            if  (dpm(i,j,klay).gt.onemu .AND. fld(i,j).gt.0.0) then
               do ic = is, ie
                  fldij = fldij + flds(ic,j)
                  areaij = areaij + arsm(ic,j)
               enddo
            endif
            ! smoothed the fld. The special points are unchanged.
            if (areaij .gt. 0.0) then
                  fld(i,j) = fldij / areaij
C             elseif (areaij .eq. 0.0) then
C                   fld(i,j) = fld(i,j)     
            endif
            ! tracermass after smooth
            fsum_nj(j) = fsum_nj(j) + fld(i,j) * scp2(i,j)
         enddo
      enddo
!OMP END PARALLEL DO
c
c     Increment due to the smooth
c
      fsum_o = 0.0
      fsum_n = 0.0
!OMP PARALLEL DO PRIVATE(j)
      do j = 1-nbdy, jj+nbdy
         fsum_o = fsum_o + fsum_oj(j)
         fsum_n = fsum_n + fsum_nj(j)
      enddo
!OMP END PARALLEL DO
      if (klay .eq. 20) then
           write(*,'(a,e12.3, a,e12.3, a,e12.3, a,i2, a,i3)'),
     &    ' SMOOTH:',fsum_o,',',fsum_n,',',fsum_n/fsum_o-1.0
     &    , ' ,Z', klay, ' ,ni=',nfili
      endif
c
      return
      end subroutine tracer_smooth

c=====================================================================
c          SUBTOUTINE: Laplacian operator
c=====================================================================
      subroutine tracer_diff(fld, fld0, klay, flgin)
c
c     Compute Laplacian operator (times dx) of the input variable
c
      implicit none
      integer l,margin,klay
      integer :: flgin
      real, parameter :: aspmax = 2.0
      real, parameter :: epsil = 1.e-11
      real, parameter :: onemu = 1.e-12 !very small layer thickness (m)
      real harmon,a,b,kxx,kyy
      real, dimension(1-nbdy:ii+nbdy,1-nbdy:jj+nbdy) :: fld,dtdx,dtdy,
     &     xflux,yflux,flds,fld0, dp1, dp2, fmx, fmn, lhs
      real :: aspux,aspvy,factor
      real :: ia, ib, ja, jb, q, q_cor
     &    ,lhs_sm, cm1,cm2
      integer :: constrainMethod = 0 ! 0, 1 or 2
C       
      harmon(a,b) = 2.*a*b / (a+b)
C 
c
      if (flgin.eq.0) then
          constrainMethod = 0
      else
          constrainMethod = 1
      endif
c
      if (klay.eq.1) then
         write(*,*) 'If constrain on DIFF (',flgin,'):',constrainMethod
      endif
C                
!$OMP PARALLEL DO PRIVATE(i,j)
      do j = 1-nbdy, jj+nbdy
         do i = 1-nbdy, ii+nbdy
            dtdx(i,j) = 0.0
            dtdy(i,j) = 0.0
            xflux(i,j) = 0.0
            yflux(i,j) = 0.0
            lhs(i,j) = 0.0     
         enddo
      enddo
!$OMP END PARALLEL DO
c
      margin = nbdy - 1
!$OMP PARALLEL DO PRIVATE(i,j)
      do j = 1-nbdy, jj+nbdy
         do i = 1-nbdy, ii+nbdy
            flds(i,j) = fld0(i,j)
         enddo
      enddo
!$OMP END PARALLEL DO
c
c --- Smooth tracer if use Ktensor outside
c
      if (flgin.ne.0) then
         call tracer_smooth(flds,50,50,klay)
      endif
c
c --- Reconstruct (-) c-flx using tracer gradient & diffusivity [c*m2/s]
c     
!$OMP   PARALLEL DO PRIVATE(i,j,factor, aspux, aspvy)                  
c!$OMP&           SCHEDULE(STATIC,jblk)
      do j = 1-margin, jj+margin
         do i = 1-margin, ii+margin
c
c ---       - c-uflx ([c*m2/s], equivalent to '-u*h*c')
c
            if (iu(i,j).ne.0) then
               ! (harmonic) mean layer thck [m]
               factor = harmon( max(dpm(i-1,j,klay),onemu),
     &            max(dpm(i,j,klay),onemu) )
               ! c-grad [c] = [m*c/m]
               dtdx(i,j) = factor * ( flds(i,j)-flds(i-1,j) ) / 
     &            max(scux(i,j),epsil)
c
c              Horizontal diffusivity
c
               if (flgin.eq.0) then
c     CONSTANT ISOTROPIC diffusivity [m2/s, temdf2*len], consistent with HYCOM
c       Largest grid spacing (within limits) used in all diffusion
c       coefficients: min( max(sc?x,sc?y), min(sc?x,sc?y)*aspmax )
c       divided by scux/vy
                  kxx = temdf2 * min( max(scux(i,j),scuy(i,j)),
     &              min(scux(i,j),scuy(i,j))*aspmax )
               elseif (flgin.eq.1) then
c     ISOTROPIC (SCALAR, scale-dependent) diffusivity [m2/s]
                  kxx = ktensor(i,j,klay,1)
               elseif (flgin.eq.2) then
c     ANISOTROPIC (2D TENSOR) diffusivity [m2/s]
                  kxx = ktensor(i,j,klay,1)
               endif
               xflux(i,j) = kxx * dtdx(i,j)
            endif               
c
c ---       c-vflx
c
            if (iv(i,j).ne.0) then
               factor = harmon( max(dpm(i,j-1,klay),onemu), 
     &            max(dpm(i,j,klay),onemu) )
               dtdy(i,j) = factor * ( flds(i,j)-flds(i,j-1) ) / 
     &            max(scvy(i,j),epsil)
               if (flgin.eq.0) then
                  kyy = temdf2 * min( max(scvx(i,j),scvy(i,j)),
     &              min(scvx(i,j),scvy(i,j))*aspmax )
               elseif (flgin.eq.1) then
                  kyy = ktensor(i,j,klay,1)
               elseif (flgin.eq.2) then
                  kyy = ktensor(i,j,klay,4)
               endif
               yflux(i,j) = kyy * dtdy(i,j)
            endif               
         enddo                  !i
      enddo                     !j                                     
!$OMP   END PARALLEL DO 
c
c     FOR ANISOTROPIC diffusivity tensor: ADD off-diagonal tensor terms
c     Recall the eqn: 
c           ucflx = - (kxx * dtdx + kxy * dtdy), kxx/xy on u-grid
c           vcflx = - (kyx * dtdx + kyy * dtdy), kyx/yy on v-grid
c
      if (flgin.eq.2) then       
!$OMP   PARALLEL DO PRIVATE(i,j,factor)
c!$OMP&           SCHEDULE(STATIC,jblk)
         do j = 1, jj      
            do i = 1, ii
               if (iu(i,j).ne.0) then
                  xflux(i,j) = xflux(i,j) + ktensor(i,j,klay,2)*0.25*
     &              (dtdy(i,j)+dtdy(i-1,j)+dtdy(i,j+1)+dtdy(i-1,j+1))
               endif
               if (iv(i,j).ne.0) then
                  yflux(i,j) = yflux(i,j) + ktensor(i,j,klay,3)*0.25*
     &              (dtdx(i,j)+dtdx(i,j-1)+dtdx(i+1,j)+dtdx(i+1,j-1))
               endif
            enddo                  !i
         enddo                     !j
!$OMP   END PARALLEL DO
      endif
c
c     Calculate the divergence of recons c-flx, and update tracer !
c     div is in [c], dt*d_cflx*len/vol
c
      if (constrainMethod.eq.1) then
!$OMP PARALLEL DO PRIVATE(i,j,ja,jb,ia,ib)
         do j = 1, jj
            ja = j - 1
            jb = j + 1
            do i = 1, ii
               if (ip(i,j).ne.0) then
                  !CONSTRAIN TRACER, same in 'fct3d'
                  ia = i - 1
                  if (ip(ia,j).eq.0) ia = i
                  ib = i + 1
                  if (ip(ib,j).eq.0) ib = i
                  ja = j - 1
                  if (ip(i,ja).eq.0) ja = j
                  jb = j + 1
                  if (ip(i,jb).eq.0) jb = j
                  fmx(i,j) = max( fld(i,j),
     &                 fld(ia,j),fld(ib,j),fld(i,ja),fld(i,jb) )
                  fmn(i,j) = min( fld(i,j),
     &                 fld(ia,j),fld(ib,j),fld(i,ja),fld(i,jb) )
                  endif
            enddo
          enddo
!$OMP END PARALLEL DO
      endif
C 
C     Layer thickness at current and next time steps
C     dp(i,j,klay,1/2) or both are max(dpm(i,j,klay),onemu)
C 
!$OMP PARALLEL DO PRIVATE(i,j)
      do j = 1-nbdy, jj+nbdy
         do i = 1-nbdy, ii+nbdy
            dp1(i,j) = dpm(i,j,klay)
            dp2(i,j) = dpm(i,j,klay)
         enddo
      enddo
!$OMP END PARALLEL DO
C 
      cm1 = 0.
      do j = 1, jj
         do i = 1, ii
              ! [c*m3]
              cm1 = cm1 + fld(i,j)*dp1(i,j)
         enddo
      enddo
C 
!$OMP   PARALLEL DO PRIVATE(i,j,q,q_cor)
      do j = 1, jj    
C          ja = j - 1
C          jb = j + 1
         do i = 1, ii
            if (ip(i,j).ne.0) then
               ! (ch)_2 = (ch)_1 - div_Fc = (ch)_1 + div_(K*del_<c>)
               lhs(i,j) = delt1 / scp2(i,j) * (
     &                xflux(i+1,j)*scuy(i+1,j) - xflux(i,j)*scuy(i,j) +
     &                yflux(i,j+1)*scvx(i,j+1) - yflux(i,j)*scvx(i,j) ) 
               ! [c*m] 
               q = fld(i,j)*dp1(i,j) + lhs(i,j)
C                
C              CONSTRAIN, if needed
C                
               if (constrainMethod.eq.0) then
                  q_cor = q
               elseif (constrainMethod.eq.1) then
C                   ia = i - 1
C                   if (ip(ia,j).eq.0) ia = i
C                   ib = i + 1
C                   if (ip(ib,j).eq.0) ib = i
C                   ja = j - 1
C                   if (ip(i,ja).eq.0) ja = j
C                   jb = j + 1
C                   if (ip(i,jb).eq.0) jb = j
C                   fmx2 = max( fld(i,j),
C      &                 fld(ia,j),fld(ib,j),fld(i,ja),fld(i,jb) )
C                   fmn2 = min( fld(i,j),
C      &                 fld(ia,j),fld(ib,j),fld(i,ja),fld(i,jb) )
                  q_cor = max( fmn(i,j)*dp2(i,j),
     &                         min(q, fmx(i,j)*dp2(i,j)) )
               endif
C                
c              Updated tracer conc [c], essentially q/h
               fld(i,j) = (fld(i,j)*onemu + q_cor) / (onemu + dp2(i,j))
C                
C              Eliminate very small values far from tracer patch
               if (fld(i,j) .lt. 1.e-7)  fld(i,j) = 0.0
            endif               !ip
         enddo                  !i
      enddo                     !j
!$OMP   END PARALLEL DO

c  domain-integrated lhs (flxdiv)
      lhs_sm = 0.  
      cm2 = 0.  
      do j = 1, jj
         do i = 1, ii
              ! [c*m3]  
              lhs_sm = lhs_sm + lhs(i,j)
              cm2 = cm2 + fld(i,j)*dp2(i,j)
         enddo 
      enddo  
c      if (klay.eq.20)  write(*,*) 'flxdiv_SM: ', lhs_sm, cm1, cm2
c     &  , (cm2-cm1)/cm1,' Z:', klay 

      return
      end subroutine tracer_diff
c
c=====================================================================
c          Apply the parameterized diffusion & advection 
c=====================================================================
      subroutine tracer_diffadv(fld, fld0, klay)
c 
c     Compute Laplacian operator & adv of input field
c       DIV_F (or ADV) = - K * del_del<c> + lambda * del<c>
c       LHS     -  [m/s*c], uh*c*L/L2 (or uh*c/L)
c       K       -  [m2/s]
c       lambda  -  [m/s], p-grid
c       del<c>  -  [c], multiplied by <h> or h
c       del2<c> -  [c/m], cgrad*L/L2
c
C     Update the tracer ('fld') using the gradient of tracer before 
C     advection ('fld0')
C 
      implicit none
      integer :: margin, klay, counts
      real, parameter :: aspmax = 2.0
      real, parameter :: epsil = 1.e-11
      real, parameter :: onemu = 1.e-12 !very small layer thickness (m)
      real :: harmon, a, b, kxx, kyy
      real, dimension(1-nbdy:ii+nbdy,1-nbdy:jj+nbdy) :: fld,dtdx,dtdy,
     &     flds, fld0, fmx, fmn, dp1, dp2, kdel2c, lmdcgrad, lhs
      real :: factor, ia, ib, ja, jb, q, q_cor, cm1, cm2, lhs_sm
      integer :: constrainMethod ! 0, 1 or 2
C       
      harmon(a,b) = 2.*a*b / (a+b)
C 
      constrainMethod = 1
      if (klay.eq.1) then
         write(lp,*) 'Constrain on KL: ', constrainMethod
      endif
C       
!OMP PARALLEL DO PRIVATE(i,j)
      do j = 1-nbdy, jj+nbdy
         do i = 1-nbdy, ii+nbdy
            dtdx(i,j) = 0.0
            dtdy(i,j) = 0.0
            lhs(i,j) = 0.0
         enddo
      enddo
!OMP END PARALLEL DO
c
      margin = nbdy - 1
!OMP PARALLEL DO PRIVATE(i,j)
      do j = 1-nbdy, jj+nbdy
         do i = 1-nbdy, ii+nbdy
            flds(i,j) = fld0(i,j)
         enddo
      enddo
!OMP END PARALLEL DO
c
c --- Smooth tracer if use Ktensor/advection outside
c
      if (difflg.ne.0) then
         call tracer_smooth(flds,50,50,klay)
      endif
c
c --- Calc tracer gradient [c], h*delta_c/delta_L
c     
!$OMP   PARALLEL DO PRIVATE(i,j,factor)                  
c!$OMP&           SCHEDULE(STATIC,jblk)
      do j = 1-margin, jj+margin
         do i = 1-margin, ii+margin
c
c ---       dc/dx on u-grid 
c
            if (iu(i,j).ne.0) then
               ! (harmonic) mean layer thck [m]
               factor = harmon( max(dpm(i-1,j,klay),onemu),
     &            max(dpm(i,j,klay),onemu) )
               dtdx(i,j) = factor * ( flds(i,j) - flds(i-1,j) ) / 
     &            max(scux(i,j),epsil)
            endif
c
c ---       dc/dy on v-grid 
c
            if (iv(i,j).ne.0) then
               factor = harmon( max(dpm(i,j-1,klay),onemu), 
     &            max(dpm(i,j,klay),onemu) )
               dtdy(i,j) = factor * ( flds(i,j) - flds(i,j-1) ) / 
     &            max(scvy(i,j),epsil)
            endif               
         enddo                  !i
      enddo                     !j                                     
!$OMP   END PARALLEL DO 
C 
C     Layer thickness at current and next time steps
C     dp(i,j,klay,1/2) or both are max(dpm(i,j,klay),onemu)
C 
!$OMP PARALLEL DO PRIVATE(i,j)
      do j = 1-nbdy, jj+nbdy
         do i = 1-nbdy, ii+nbdy
            dp1(i,j) = dpm(i,j,klay)
            dp2(i,j) = dpm(i,j,klay)
         enddo
      enddo
!$OMP END PARALLEL DO
      cm1 = 0.
      do j = 1, jj
         do i = 1, ii
              ! [c*m3]
              cm1 = cm1 + fld(i,j)*dp1(i,j)
         enddo
      enddo
c
c     Calculate K*DIV_cgrad - lmd*cgrad [c*m/s], and update tracer
c     When updating, compute dt*LHS/h ([c]), consistent with FLUX method
c      I.e. f = f + K * del_del<c> - lambda * del<c>
c      
      if (constrainMethod.eq.1) then
!$OMP PARALLEL DO PRIVATE(i,j,ja,jb,ia,ib)
         do j = 1, jj
            ja = j - 1
            jb = j + 1
            do i = 1, ii
               if (ip(i,j).ne.0) then
                  !CONSTRAIN TRACER, same in 'fct3d'
                  ia = i - 1
                  if (ip(ia,j).eq.0) ia = i
                  ib = i + 1
                  if (ip(ib,j).eq.0) ib = i
                  ja = j - 1
                  if (ip(i,ja).eq.0) ja = j
                  jb = j + 1
                  if (ip(i,jb).eq.0) jb = j
                  fmx(i,j) = max( fld(i,j),
     &                 fld(ia,j),fld(ib,j),fld(i,ja),fld(i,jb) )
                  fmn(i,j) = min( fld(i,j),
     &                 fld(ia,j),fld(ib,j),fld(i,ja),fld(i,jb) )
                  endif
            enddo
          enddo
!$OMP END PARALLEL DO
      endif
C 
!$OMP   PARALLEL DO PRIVATE(i,j,q,q_cor)
c!$OMP&           SCHEDULE(STATIC,jblk)
      do j = 1, jj    
         do i = 1, ii
            if (ip(i,j).ne.0) then
               ! K * DIV_del<c> [m2/s/m2 * cgrad*m] = [m/s*c]
               kdel2c(i,j) = ktensor(i,j,klay,1) / scp2(i,j) * (
     &            dtdx(i+1,j)*scuy(i+1,j) - dtdx(i,j)*scuy(i,j) +
     &            dtdy(i,j+1)*scvx(i,j+1) - dtdy(i,j)*scvx(i,j) ) 
               ! lmd*del<c> [m/s*c]
               lmdcgrad(i,j) = ktensor(i,j,klay,2) * 
     &            0.5*(dtdx(i,j) + dtdx(i+1,j))
     &            + ktensor(i,j,klay,3) * 0.5*(dtdy(i,j) + dtdy(i,j+1))
               ! add the two [m*c]
               lhs(i,j) = delt1 * (kdel2c(i,j) - lmdcgrad(i,j))
            endif               !ip
         enddo                  !i
      enddo                     !j
!$OMP   END PARALLEL DO
C 
c  domain-integrated lhs (represented eddy forcing)
      lhs_sm = 0.
      counts = 0
      do j = 1, jj
         do i = 1, ii
              ! [c*m3]
              lhs_sm = lhs_sm + lhs(i,j)
              if (lhs(i,j).ne.0.) counts = counts + 1
         enddo
      enddo
      lhs_sm = lhs_sm / float(counts)
      if (klay.eq.15)  write(*,*) 'Domain mean(LHS m*c)=', lhs_sm, 
     &   'counts=',counts,' lhs(550,550)=', lhs(550,550),'Z:', klay
C 
C     Update tracer by corrected eddy forcing
C 
      do j = 1, jj    
         do i = 1, ii
            if (ip(i,j).ne.0) then
C                
C              Calc the updated tracer conc [c]
C              With domain mean forcing SUBTRACTED!!!
C 
               if (lhs(i,j).ne.0.) then
                  lhs(i,j) = lhs(i,j) - lhs_sm
               endif
               q = fld(i,j)*dp1(i,j) + lhs(i,j)
C                
c              Choices of CONSTRAIN
c
               if (constrainMethod.eq.0) then
                  q_cor = q
               elseif (constrainMethod.eq.1) then
                  q_cor = max( fmn(i,j)*dp2(i,j),
     &                         min(q, fmx(i,j)*dp2(i,j)) )
               endif
C                
c              Updated tracer conc [c], essentially q/h
               fld(i,j) = (fld(i,j)*onemu + q_cor) / (onemu + dp2(i,j))
C                
C              Eliminate very small values far from tracer patch
               if (fld(i,j) .lt. 1.e-7)  fld(i,j) = 0.0
            endif               !ip
         enddo                  !i
      enddo                     !j


c      if (klay.eq.20)  write(*,*) 'flxdiv_SM: ', flxdiv_sm, cm1, cm2
c     &  , (cm2-cm1)/cm1,' Z:', klay
C 
      return
      end subroutine tracer_diffadv

c=====================================================================
c          SUBTOUTINE: k-profile vertical diffusion
c=====================================================================
      subroutine tracer_kpp(i,j)
c --- hycom version 1.0
      implicit none
c
c -------------------------------------------------------------
c --- k-profile vertical diffusion, single j-row (part B)
c --- vertical coordinate is z negative below the ocean surface
c -------------------------------------------------------------
c
c --- local 1-d arrays for matrix solution
      real tr1do(kk+1,ntracr),tr1dn(kk+1,ntracr),dpmm(kk),pij(kk+1),
     &     difftr(kk+1),ghat(kk+1),zm(kk+1),hm(kk),dzb(kk),
     $     hwide(kk)
c
c --- tridiagonal matrix solution arrays
      real tri(kk,0:1)      ! dt/dz/dz factors in trid. matrix
      real tcu(kk),         ! upper coeff for (k-1) on k line of trid.matrix
     &     tcc(kk),         ! central ...     (k  ) ..
     &     tcl(kk),         ! lower .....     (k-1) ..
     &     rhs(kk)          ! right-hand-side terms
      real    ghatflux, qonem
      parameter(qonem=1.0/onem)
      integer k,ka,ktr,nlayer,n,i,j
c
c      include 'stmt_fns.h'
c --- calculate z at vertical grid levels - this array is the z values in m
c --- at the mid-depth of each micom layer except for index klist+1, where it
c --- is the z value of the bottom
c
      nlayer=klist(i,j)
c
c --- locate lowest substantial mass-containing layer.
      do k=1,kk
         dpmm( k)  =max(0.001,dpm(i,j,k))
c     qdpmm(k)  =1.0/dpmm(k)
c     pij(k) = p(i,j,k)
      enddo
c --- calculate layer thickness in meters
      do k=1,kk
         if (k.eq.1) then
            hwide(k)=dpmm(k)
            zm(k)=-.5*hwide(k)  ! zgrid(i,j,k)
         else if (k.lt.nlayer) then
            hwide(k)=dpmm(k)
            zm(k)=zm(k-1)-.5*(hwide(k-1)+hwide(k)) !zgrid(i,j,k)
         else if (k.eq.nlayer) then
            hwide(k)=dpmm(k)
            zm(k)=zm(k-1)-.5*(hwide(k-1)+hwide(k)) ! zgrid(i,j,k)
            zm(k+1)=zm(k)-.5*hwide(k) ! zgrid(i,j,k+1)
         else
            hwide(k)=0.
         endif
      enddo
c
      do k=1,nlayer
         difftr(k+1)=difs(i,j,k) ! difs is saved as k+1
         if (ghtflg.eq.1) then
            ghat(k+1)= ghats(i,j,k)
         else
            ghat(k+1)=0.0
         endif
         do ktr= 1,ntracr
            tr1do(k,ktr)=tracer(i,j,k,ktr)
         enddo
         hm(k)=max(0.001,dpm(i,j,k))
c     zm(k)=zgrid(i,j,k)
      enddo                     !k
c
      k=nlayer+1
      ka=min(k,kk)
      difftr(k)=0.0
      ghat(k)=0.0
      do ktr= 1,ntracr
         tr1do(k,ktr)=tracer(i,j,ka,ktr)
      enddo
      zm(k)=zm(k-1)-0.001
c
c --- compute factors for coefficients of tridiagonal matrix elements.
c     tri(k=1:NZ,0) : dt/hwide(k)/ dzb(k-1)=z(k-1)-z(k)=dzabove)
c     tri(k=1:NZ,1) : dt/hwide(k)/(dzb(k  )=z(k)-z(k+1)=dzbelow)
c
      do k=1,nlayer
        dzb(k)=zm(k)-zm(k+1)
      enddo
c
      tri(1,1)=delt1/(hm(1)*dzb(1))
      tri(1,0)=0.
      do k=2,nlayer
        tri(k,1)=delt1/(hm(k)*dzb(k))
        tri(k,0)=delt1/(hm(k)*dzb(k-1))
      enddo
c
c --- solve the diffusion equation
c --- salflx, sswflx and surflx are positive into the ocean
c
c --- t-like tracer solution
c      do ktr= 1,ntracr
c        if     (trcflg(ktr).eq.2) then
c          ghatflux=-(surflx(i,j)-sswflx(i,j))*thref/spcifh
c          call tridrhs(hm,
c     &                 tr1do(1,ktr),difft,ghat,ghatflux,tri,nlayer,rhs)
c          call tridmat(tcu,tcc,tcl,nlayer,hm,rhs,
c     &                 tr1do(1,ktr),tr1dn(1,ktr),difft, i,j)
c        endif
c      enddo
c
c --- standard tracer solution
      if     (ntracr.gt.0) then
         call tridcof(difftr,tri,nlayer,tcu,tcc,tcl)
      endif
      do ktr = 1,ntracr
         ghatflux = 0.
         call tridrhs(hm,delt1,
     $        tr1do(1,ktr),difftr,ghat,ghatflux,tri,nlayer,rhs)
         call tridmat(tcu,tcc,tcl,nlayer,hm,rhs,
     &        tr1do(1,ktr),tr1dn(1,ktr),difftr, i,j)
      enddo
      do k = 1,nlayer
         do ktr = 1,ntracr
            tracer(i,j,k,ktr) = tr1dn(k,ktr)
         enddo                  !ktr
      enddo                     !k
c
      return
      end subroutine tracer_kpp
c

c=====================================================================
c          SUBTOUTINE: homogenization within the mixed layer
c=====================================================================
      subroutine tracer_mlhom(i,j)
      implicit none
      real sumtr, sumdp
      real p1d(kk+1)
      integer kmixl, i, j

      p1d(1) = 0
      kmixl = 0
C       
C     p1d(k+1) = dpm_z1 + dpm_z1 + ... + dpm_zk
C 
      do k = 1, kk
         p1d(k+1) = p1d(k) + dp(i,j,k,2)
      enddo
C       
C     kmixl is the # of layer whose bottom is below MLD
C 
      do k = 1, kk
         if (p1d(k+1).le.dpmixl(i,j)) then
            kmixl = k
         endif
      enddo      
C 
C     Homogenize by SUM_c(z)/SUM_h(z) for current cell (i,j). Note this 
C     is conservative for c-mass. 
C 
      sumtr = 0.
      sumdp = 0.
      do k = 1, kmixl
        sumtr = sumtr + tracer(i,j,k,ktr)*dp(i,j,k,2)
        sumdp = sumdp + dp(i,j,k,2)
      enddo
      do k = 1, kmixl
        ! > 1m
        if (sumdp.gt.1.0) then
           tracer(i,j,k,ktr) = sumtr/sumdp
        endif
      enddo                  !k
      return
      end subroutine tracer_mlhom  
      
c=====================================================================
c          SUBTOUTINE: output tracer field
c=====================================================================
      subroutine tracer_output(tim,flnm)
      double precision tim
      character*120 flnm 
      real hminb,hmaxb,xmin,xmax
c      real, dimension(idm,jdm) :: dumm2d
      integer nto,ntr
      character cline*80
      data nto/22/
      real*4, parameter   :: spval=2.0**100
      real cmass_out
C       
c     write output tracer restart file
      write(6,*)trim(flnm)//'.b'
      call flush(lp)
      open (unit=nto,file=trim(flnm)//'.b',form='formatted',
     &      status='new',action='write')
      write(6,*)trim(flnm)//'.a'
      call flush(lp)
      call zaiopf(trim(flnm)//'.a','new', nto)
c     write header
      do ktr=1,ntracr
         do k=1,kk
            ntr=(ktr-1)*kk+k
            write(nto+1000,rec=ntr) ((tracer(i,j,k,ktr),
     $           i=1-nbdy,ii+nbdy),j=1-nbdy,jj+nbdy)
            xmin=spval
            xmax=-spval
            do i=1-nbdy,ii+nbdy
               do j=1-nbdy,jj+nbdy
                  if (ip(i,j).ne.0) then
                     xmin=min(tracer(i,j,k,ktr),xmin)
                     xmax=max(tracer(i,j,k,ktr),xmax)
                  endif
               enddo
            enddo
            write(nto,117)'tracer  ',tim,k,xmin,xmax
         enddo
c        calc c mass
         call tracer_mass(tracer(1-nbdy,1-nbdy,1,ktr),0,cmass_out) 
         write(nto,*)'C-',ktr, 'cmass all layers:', cmass_out
      enddo
 117  format (a8,' =',f11.3,i3,1p2e16.7)
c
      call flush(nto)
      close(unit=nto)
      call zaiocl(nto)
c
      return
      end subroutine tracer_output

c=====================================================================
c          SUBTOUTINE: buffer points 
c=====================================================================
      subroutine tracer_buffer(fld)
c
      real fld(1-nbdy:ii+nbdy,1-nbdy:jj+nbdy)
c
c --- DEFINE BUFFER POINTS, using linear extrapolation (which may lead to
C     instabilities) OR by extending the border points
c
!$OMP PARALLEL DO PRIVATE(i,j)
      do j=1-nbdy,jj+nbdy
         do i=1,nbdy
c            fld( 1-i,j)=float(i+1)*fld( 1,j)-float(i)*fld(   2,j)
c            fld(ii+i,j)=float(i+1)*fld(ii,j)-float(i)*fld(ii-1,j)
            fld( 1-i,j)=fld(1,j)
            fld(ii+i,j)=fld(ii,j)
         enddo
      enddo
!$OMP END PARALLEL DO
!$OMP PARALLEL DO PRIVATE(i,j)
      do i=1-nbdy,ii+nbdy
         do j=1,nbdy
c            fld(i, 1-j)=float(j+1)*fld(i, 1)-float(j)*fld(i,   2)
c            fld(i,jj+j)=float(j+1)*fld(i,jj)-float(j)*fld(i,jj-1)
            fld(i, 1-j)=fld(i,1)
            fld(i,jj+j)=fld(i,jj)
         enddo
      enddo
!$OMP END PARALLEL DO
c
      return
      end subroutine tracer_buffer

c=====================================================================
c          SUBTOUTINE: volume flux (u*dp*len*dt)
c=====================================================================
      subroutine tracer_uvflx(q0,uflflg)
c
c --- CALCULATE VOLUME FLUXES FROM VELOCITIES AND DP AT P-POINTS 
C     (IF NEEDED) AND MULTIPLY BY THE GEOMETRICAL COEFFICIENTS
c
      integer uflflg
      real depthi(0:1,0:1,0:kk),duk,dukm1,dvk,dvkm1,q0
      if (uflflg.eq.1) then
c
c --- Calculate DPU, DPV and mass fluxes
c
         print*,'Calculate volume fluxes!'
!$OMP PARALLEL DO PRIVATE(i,j,k,depthi, duk, dvk, dukm1, dvkm1)
         do j= 2-nbdy,jj+nbdy
            do i= 2-nbdy,ii+nbdy
               duk           = 0.0
               dvk           = 0.0
               depthi(:,:,0) = 0.0
               do k= 1,kk
                  depthi(1,1,k) =
     &                 depthi(1,1,k-1) + dpm(i  ,j,k)
                  depthi(0,1,k) =
     &                 depthi(0,1,k-1) + dpm(i-1,j,k)
                  depthi(1,0,k) =
     &                 depthi(1,0,k-1) + dpm(i,j-1,k)
                  dukm1 = duk
                  duk   = min( depthu(i,j),
     &                 0.5*(depthi(0,1,k) + depthi(1,1,k)) )
                  dvkm1 = dvk
                  dvk   = min( depthv(i,j),
     &                 0.5*(depthi(1,0,k) + depthi(1,1,k)) )
                  dpu(i,j,k) = max( 0.0, duk-dukm1 )
                  dpv(i,j,k) = max( 0.0, dvk-dvkm1 )
               enddo            !k
            enddo               !i
         enddo                  !j
!$OMP END PARALLEL DO
c         call pr_5x5(dpu(1,1,1),ii,jj,itest,jtest,0,1.,'dpu')
c         call pr_5x5(dpv(1,1,1),ii,jj,itest,jtest,0,1.,'dpv')
         print*,'DPU=',dpu(640,1010,:)
         print*,'U=',uflx(640,1010,:)
!$OMP PARALLEL DO PRIVATE(i,j,k)
         do k=1,kk
            do j=2-nbdy,jj+nbdy
               do i=2-nbdy,ii+nbdy
                  if(iu(i,j).eq.1) then
                     uflx(i,j,k)=uflx(i,j,k)*dpu(i,j,k)
                  else
                     uflx(i,j,k)=0.
                  endif
                  if(iv(i,j).eq.1) then
                     vflx(i,j,k)=vflx(i,j,k)*dpv(i,j,k)
                  else
                     vflx(i,j,k)=0.
                  endif
               enddo
            enddo
         enddo
!$OMP END PARALLEL DO
      endif                     ! calculation of the volume flux
c
!$OMP PARALLEL DO PRIVATE(i,j,k)
      do k=1,kk
         do j=1-nbdy,jj+nbdy
            do i=1-nbdy,ii+nbdy
               if(iu(i,j).eq.1) then
                  uflx(i,j,k)=uflx(i,j,k)*scuy(i,j)*q0
               else
                  uflx(i,j,k)=0.
               endif
               if(iv(i,j).eq.1) then
                  vflx(i,j,k)=vflx(i,j,k)*scvx(i,j)*q0
               else
                  vflx(i,j,k)=0.
               endif
            enddo
         enddo
      enddo
c
      return
      end subroutine tracer_uvflx       

c=====================================================================
c          SUBTOUTINE: 5 x 5 point cluster
c=====================================================================
      subroutine pr_5x5(array,ii,jj,iz,jz,offset,scale,what)
c
c --- (same as prt9x9 except that array dimensions are passed as arguments)
c
c --- write 5 x 5 point cluster of 'array' values centered on (iz,jz).
c --- the printed numbers actually represent (array(i,j) + offset) * scale
c
      implicit none
      integer,intent(IN) :: ii,jj,iz,jz
      real   ,intent(IN) :: array(ii,jj),scale,offset
      character what*(*),text*12
      integer i,j,jwrap
 100  format(a12,5i13)
 101  format(i10,3x,5f13.5)
c
      text='            '
      text(1:min(12,len_trim(what)))=what(1:min(12,len_trim(what)))
      write (*,100) text,(i,i= iz-2,iz+2)
      write (*,101) (j,(scale*(array(i,j)+offset),
     . i=iz-2,iz+2),j=jz+2,jz-2,-1)
      return
      end subroutine pr_5x5


c=====================================================================
c          SUBTOUTINE: 3-d version of advfct.f
c=====================================================================
      subroutine fct3d(iord,fld,fco,fc)
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
c      hordivE  - 
c
c     Some inmportant vars in the adv:
c       bforej/afterj - (jdm), amount of tracermass in each lat band in 
c                       current layer at previous/new time steps
c       bfore/after   - amount of tracermass in the whole domain at
c                       previous/new time steps
c
      implicit none
c
      real, dimension(1-nbdy:ii+nbdy,1-nbdy:jj+nbdy,kdm)::
     $     fld,fco,fc,vertfx,vertdv
      real, dimension(1-nbdy:ii+nbdy,1-nbdy:jj+nbdy) ::
     $     flp,fln,flx,fly,uan,van,flxdiv,fmx,fmn
      real, dimension(1-nbdy:jj+nbdy) :: clipj, vlumj
      real a(kdm),b(kdm),c(kdm),athird,dx,fcdx,yl,yr
      real onemu,q,clip,vlume,amount,bfore,after,slab,dslab,thkchg,
     .     fluxdv,epsil,bforej(1-nbdy:jj+nbdy),afterj(1-nbdy:jj+nbdy)
      integer iord,ip1,im1,jp1,jm1,kp,jaa,margin,i1,i2,j1,j2,
     $     ia,ib,ja,jb,l
      logical wrap,recovr
      data recovr/.false./
      data (athird = 1./3.)
      real,    parameter :: ft14 = 7.0/6.0,  !4th centered inner coeff
     &                      ft24 =-1.0/6.0   !4th centered outer coeff
C   
c
c --- if iord=1, scheme reduces to simple donor cell scheme. usually it
c     is set to 2 in our model
      parameter (epsil = 1.e-11, onemu = 1.e-6)

c ---------------------------------------------------------------------
c --- get vertical c-flux by summing over upstream slab of thickness 
c ---------------------------------------------------------------------
      margin = nbdy - 1
!$OMP PARALLEL DO PRIVATE(jb,amount,slab,dslab,kp,a,b,c,dx,fcdx,yl,yr)
      DO j = 1-margin, jj+margin
         DO i = 1-margin, ii+margin
c
c ---       fill massless cells with data from layer above or below
c
            if (ip(i,j).ne.0) then
               do k = kk-1, 1, -1
                  fld(i,j,k)=(fld(i,j,k)*fco(i,j,k)+fld(i,j,k+1)*onemu)
     .                 /(           fco(i,j,k)+             onemu)
               enddo
               do k = 2, kk
                  fld(i,j,k)=(fld(i,j,k)*fco(i,j,k)+fld(i,j,k-1)*onemu)
     .                 /(           fco(i,j,k)+             onemu)
               enddo
c
c ---          Parabolic method using 3 points: fit (vertically)
c              a+bx+cx^2 to fld in each cell (-.5 < x < +.5)
c
               a(1 ) = fld(i,j,1 )
               b(1 ) = 0.
               c(1 ) = 0.
               a(kk) = fld(i,j,kk)
               b(kk) = 0.
               c(kk) = 0.
c
               do k = 2, kk-1
                  yl = .5 * ( fld(i,j,k-1) + fld(i,j,k) )
                  yr = .5 * ( fld(i,j,k+1) + fld(i,j,k) )
                  a(k) = 1.5 * fld(i,j,k) - .25 * (yl + yr)
                  b(k) = yr - yl
                  c(k) = 6. * ( .5 * (yl + yr) - fld(i,j,k) )
                  if (abs(yr-yl).lt.6.*abs(.5*(yl+yr)-fld(i,j,k))) then
c --- apex of parabola occurs inside interval [-.5,+.5], implying an over-
c     or undershoot situation. change curve to prevent over/undershoots.
                     if (abs(yr-yl).gt.2.*abs(.5*(yl+yr)-fld(i,j,k))) 
     &               then
c     --- put apex of parabola on edge of interval [-.5,+.5]
                        if ((yr-yl)*(.5*(yl+yr)-fld(i,j,k)).gt.0.) then
c     --- apex at x=-.5
                           a(k) = .25*(3.*fld(i,j,k)+yl)
                           c(k) = 3.*(fld(i,j,k)-yl)
                           b(k) = c(k)
                        else
c     --- apex at x=+.5
                           a(k) = .25*(3.*fld(i,j,k)+yr)
                           c(k) = 3.*(fld(i,j,k)-yr)
                           b(k) = -c(k)
                        end if
                     else       !  -1/6 < x < +1/6
c     --- moving apex won't help. replace parabola by constant.
                        a(k) = fld(i,j,k)
                        b(k) = 0.
                        c(k) = 0.
                     endif
                  endif
               enddo ! k = 2, kk-1  End the fiting of parabola
c
c              Calc diapycnal tracer flux ('vertfx' [m*c], diaflx * c)  
c              through layer interface, considering direction of 'diaflx'
c
               do k = 1, kk-1
                  slab = onemu
c                 water from k+1 to k at the interface (botm) of layer k
                  if (diaflx(i,j,k).lt.0.) then 
                     amount = slab * fld(i,j,k+1)
                     kp = k
 24                  kp = kp + 1
                     ! if diaflx is too large, 2nd or 3rd ajacent layers
                     ! are needed !!
                     if (slab .ge. -diaflx(i,j,k)) goto 23
                     if (fco(i,j,kp) .gt. 0.) then
                        dslab = min(slab+fco(i,j,kp), -diaflx(i,j,k))
     .                         -min(slab            , -diaflx(i,j,k))
                        dx = dslab / fco(i,j,kp) ! no dim
                        ! ??? fcdx  [c]
                        fcdx = a(kp)
     .                        +b(kp)*.5*(dx-1.) !  not needed in pcm
     .                        +c(kp)*(.25-dx*(.5-dx*athird)) !  not needed in pcm,plm
                        amount = amount + fcdx * dslab ! [c]
                        slab = slab + dslab
                     endif
                     if (kp.lt.kk) go to 24
c                 water from k to k+1 at the interface (botm) of layer k
                  else if (diaflx(i,j,k).gt.0.) then
                     amount = slab * fld(i,j,k)
                     kp = k + 1
 25                  kp = kp - 1
                     if (slab.ge.diaflx(i,j,k)) goto 23
                     if (fco(i,j,kp).gt.0.) then
                        dslab = min(slab+fco(i,j,kp), diaflx(i,j,k))
     .                         -min(slab            , diaflx(i,j,k))
                        dx = dslab / fco(i,j,kp)
                        fcdx = a(kp)
     .                        +b(kp)*.5*(1.-dx) !  not needed in pcm
     .                        +c(kp)*(.25-dx*(.5-dx*athird)) !  not needed in pcm,plm
                        amount = amount + fcdx*dslab
                        slab = slab + dslab
                     endif
                     if (kp.gt.2) go to 25
                  endif ! if (diaflx(i,j,k).lt.0.) 
 23               vertfx(i,j,k) = diaflx(i,j,k) * amount / slab
               enddo
c
c              vertical div (deeper vertfx - shallower vertfx) [m*c]
c
               vertfx(i,j,kk) = 0. ! No flux through bottom
               vertdv(i,j,1) = vertfx(i,j,1)
               do k = 2, kk
                  vertdv(i,j,k) = vertfx(i,j,k) - vertfx(i,j,k-1)
               enddo            
            endif               ! endif ip(i,j) ~= 0
         ENDDO                  ! i
      ENDDO                     ! j
!$OMP END PARALLEL DO

c ---------------------------------------------------------------------
c ---  ADV
c ---------------------------------------------------------------------
      bfore = 0.
      after = 0.
      margin = nbdy - 2
c
      DO k = 1, kk
!$OMP PARALLEL DO SHARED(k)
c
c ---    Tracer mass per lat band in layer k at previous time
c
         do j = 1, jj
            bforej(j) = 0.
            do i = 1, ii
 14            if (ip(i,j).ne.0) bforej(j) = bforej(j) + 
     $             fld(i,j,k) * fco(i,j,k) * scp2(i,j)
            enddo               !i
         enddo                  !j
!$OMP END PARALLEL DO
c
c ---    Compute antidiffusive tracer fluxes (high- minus low-order)
c ---    at x, y and z directions [m3*c]
c   -    flx/fly - low-order c-flux (uflx*c) [m3*c], used for 'flxdiv'
c   -    q       - 4th order tracer [c]
c   -    uan/van - high-order flux minus low-order flux
c   -    fmx/fmn - max/min values among tracers in 5 cells and the 
c                  adjacent inward vertfx/diaflx (both positive)
c
!$OMP PARALLEL DO PRIVATE(i,j,ja,jaa,jb,q,ia,ib,i1,i2) SHARED(k)
         do j = 1-margin, jj+margin
            ja = j - 1
            jb = j + 1
            jaa = j - 2
            do i = 1-margin, ii+margin
               ! c-flux in x dir
               if (iu(i,j).ne.0) then
                  if (uflx(i,j,k).ge.0.) then
                     q = fld(i-1,j,k)
                  else
                     q = fld(i  ,j,k)
                  endif
                  flx(i,j) = uflx(i,j,k) * q
                  q = fld(i,j,k) + fld(i-1,j,k) !  2nd order
                  if (ip(i+1,j)+iu(i-1,j).eq.2) then !  4th order
                     q = 1.125*q - .125*(fld(i+1,j,k) + fld(i-2,j,k)) 
                  endif
                  uan(i,j) = .5*q*uflx(i,j,k) - flx(i,j)
               endif    
               ! c-flux in y dir
               if (iv(i,j).ne.0) then
                  if (vflx(i,j,k).ge.0.) then
                     q = fld(i,j-1,k)
                  else
                     q = fld(i,j  ,k)
                  endif
                  fly(i,j) = vflx(i,j,k) * q
                  q = fld(i,j-1,k) + fld(i,j,k) !  2nd order
                  if (ip(i,j+1 )+iv(i,j-1).eq.2) then !  4th order
                     q = 1.125*q - .125*(fld(i,j+1,k) + fld(i,j-2,k)) 
                  endif
                  van(i,j) = .5*q*vflx(i,j,k) - fly(i,j)
               endif           
               ! c-flux in z dir
               if (ip(i,j).ne.0) then
                  ia = i - 1
                  if (ip(ia,j).eq.0) ia = i
                  ib = i + 1
                  if (ip(ib,j).eq.0) ib = i
                  ja = j - 1
                  if (ip(i,ja).eq.0) ja = j
                  jb = j + 1
                  if (ip(i,jb).eq.0) jb = j
                  fmx(i,j) = max(fld(i,j,k),
     .                 fld(ia,j,k),fld(ib,j,k),fld(i,ja,k),fld(i,jb,k))
                  fmn(i,j) = min(fld(i,j,k),
     .                 fld(ia,j,k),fld(ib,j,k),fld(i,ja,k),fld(i,jb,k))
C                   
C     fmx and fmn used in diffusion
C 
                  if (k.lt.kk) then
                     if (diaflx(i,j,k  ).lt.0.) then
                        fmx(i,j) = max( fmx(i,j), 
     &                             vertfx(i,j,k  )/diaflx(i,j,k  ) )
                        fmn(i,j) = min( fmn(i,j),
     &                             vertfx(i,j,k  )/diaflx(i,j,k  ) )
                     endif
                  endif
                  if (k.gt.1) then
                     if (diaflx(i,j,k-1).gt.0.) then
                        fmx(i,j) = max( fmx(i,j),
     &                             vertfx(i,j,k-1)/diaflx(i,j,k-1) )
                        fmn(i,j) = min( fmn(i,j),
     &                             vertfx(i,j,k-1)/diaflx(i,j,k-1) )
                     endif
                  endif
               endif            !ip
            enddo               !i
c
c ---       Boundaries for fluxes, u and v
c
            do l = 1, isp(j)
               i1 = max(1-margin, ifp(j,l))
               i2 = min(ii+margin,ilp(j,l))
               flx(i1  ,j) = 0.
               flx(i2+1,j) = 0.
               uan(i1  ,j) = 0.
               uan(i2+1,j) = 0.
            enddo               !l
         enddo                  !j
!$OMP END PARALLEL DO
!$OMP PARALLEL DO PRIVATE(j)
         do i = 1-margin, ii+margin
            wrap = jfv(i,1).eq.1 ! true if j=1 and j=jj are both water points
            do l = 1, jsp(i)
               j1 = max(1-margin, jfp(i,l)  )
               j2 = min(jj+margin,jlp(i,l)+1)
               if (j1.gt.1-margin .or. .not.wrap) then
                  fly(i,j1) = 0.
                  van(i,j1) = 0.
               endif
               if (j2.gt.1-margin .or. .not.wrap) then
                  fly(i,j2) = 0.
                  van(i,j2) = 0.
               endif
            enddo               !l
         enddo                  !i
!$OMP END PARALLEL DO
c
c ---    UPDATE the tracer in current layer with low-order c-flux
c
!$OMP PARALLEL DO PRIVATE(ib, jb,q,amount) SHARED(k)
         do j = 1-margin, jj+margin
            if (recovr) vlumj(j) = 0. ! recovr = false
            if (recovr) clipj(j) = 0.
            do i = 1-margin, ii+margin
               if (ip(i,j).ne.0) then
                  flxdiv(i,j)=(flx(i+1,j)-flx(i,j)+fly(i,j+1)-fly(i,j))
     $                 *scp2i(i,j)
                  ! updated tracer mass calc from c-eqn
                  ! (ch)_2 = (ch)_1 - div_Fc - div_wc
                  q = fld(i,j,k)*fco(i,j,k)-flxdiv(i,j)-vertdv(i,j,k) 
                  ! THIS IS TO KEEP TRACER WITHIN REASONABLE VALUES !
                  ! q(q<fmn)=fmn, q(q>fmx)=fmx, q(fmn<q<fmx)=q
                  amount = max( fmn(i,j)*fc(i,j,k), 
     $                          min(q, fmx(i,j)*fc(i,j,k)) )
                  if (recovr) then
                     vlumj(j)=vlumj(j)+scp2(i,j)*fc(i,j,k)
                     clipj(j)=clipj(j)+(q-amount)*scp2(i,j)
                  endif
                  ! updated tracer conc, essentially q/h, corrected to 
                  ! avoid negative values
                  fld(i,j,k)=(fld(i,j,k)*onemu+amount)/(onemu+fc(i,j,k))
               endif            !ip
            enddo               !i
         enddo                  !j
!$OMP END PARALLEL DO
c
         if (iord.le.1) go to 100
c
c --- at each grid point, determine the ratio of the largest permissible
c --- pos. (neg.) change in -fld- to the sum of all incoming (outgoing) fluxes
c
!$OMP PARALLEL DO SHARED(k)
         do j = 1-margin, jj+margin
ccc      jb=mod(j     ,jj)+1
            do i = 1-margin, ii+margin
               if (ip(i,j).ne.0) then
ccc      do 12 i=ifp(j,l),ilp(j,l)
                  flp(i,j) = (fmx(i,j)-fld(i,j,k))*fc(i,j,k)
     .                 /((max(0.,uan(i,j))-min(0.,uan(i+1,j))
     .                 +max(0.,van(i,j))-min(0.,van(i,j+1))+epsil)
     $                 *scp2i(i,j))
c
 12               fln(i,j) = (fmn(i,j)-fld(i,j,k))*fc(i,j,k)
     .                 /((min(0.,uan(i,j))-max(0.,uan(i+1,j))
     .                 +min(0.,van(i,j))-max(0.,van(i,j+1))-epsil)
     $                 *scp2i(i,j))
               endif            !ip
            enddo               !i
         enddo                  !j
!$OMP END PARALLEL DO
c
c----    Limit antidiffusive fluxes ("clipping") by multiplying 'clip'
c
!$OMP PARALLEL DO SHARED(k)
         do j = 1-margin, jj+margin
            do i = 1-margin, ii+margin
               if (iu(i,j).ne.0) then
ccc      do 7 i=ifu(j,l),ilu(j,l)
                  if (uan(i,j).ge.0.) then
                     clip = min(1.,flp(i,j),fln(i-1,j))
                  else
                     clip = min(1.,fln(i,j),flp(i-1,j))
                  endif
 7                flx(i,j) = uan(i,j)*clip
               endif           !iu
               if (iv(i,j).ne.0) then
ccc      do 8 i=ifv(j,l),ilv(j,l)
                  if (van(i,j).ge.0.) then
                     clip = min(1.,flp(i,j),fln(i,j-1))
                  else
                     clip = min(1.,fln(i,j),flp(i,j-1))
                  endif
 8                fly(i,j) = van(i,j)*clip
               endif           !iv
            enddo               !i
         enddo                  !j
!$OMP END PARALLEL DO
c
c ---    UPDATE (for the 2nd time) with 'clipped' low-order flux
c
!$OMP PARALLEL DO PRIVATE(amount,q) SHARED(k)
         do j = 1-margin, jj+margin
            do i = 1-margin, ii+margin
               if (ip(i,j).ne.0) then
ccc      do 62 i=ifp(j,l),ilp(j,l)
                  flxdiv(i,j)=(flx(i+1,j)-flx(i,j)+fly(i,j+1)-fly(i,j))
     $                 *scp2i(i,j)
                  q = fld(i,j,k)*fc(i,j,k) - flxdiv(i,j)
                amount=max(fmn(i,j)*fc(i,j,k),min(q,fmx(i,j)*fc(i,j,k)))
                  if (recovr) clipj(j)=clipj(j)+(q-amount)*scp2(i,j)
 62               fld(i,j,k)=(fld(i,j,k)*onemu+amount)/(onemu+fc(i,j,k))
               endif           !ip
            enddo               !i
         enddo                  !j
!$OMP END PARALLEL DO
c
cdiag call findmx(ip,fld(1,1,k),idm,ii1,jj,'fld after 62')
c
 100     continue
c
c ---    Recover (if needed) 'clipped' amount and return to field layer by layer
c
         if (recovr) then
            vlume=0.
            clip=0.
c
            do j = 1-margin, jj+margin
               vlume = vlume + vlumj(j)
 19            clip = clip + clipj(j)
            enddo
c     
            if (vlume.ne.0.) then
               clip = clip / vlume
            write (lp,'(i2,a,1pe11.3)') k,'  fld drift in fct3d',-clip
!$OMP PARALLEL DO SHARED(k)
               do j=1-margin,jj+margin
                  do i=1-margin,ii+margin
                     if (ip(i,j).ne.0) then
ccc     do 13 i=ifp(j,l),ilp(j,l)
 13                     fld(i,j,k)=fld(i,j,k)+clip
                     endif     !ip
                  enddo         !i
               enddo            !j
!$OMP END PARALLEL DO
            endif
         endif                 !recover
c
c ---    
c
!$OMP PARALLEL DO SHARED(k)
         do j = 1, jj
            afterj(j) = 0.
            do i = 1, ii
 15            if (ip(i,j).ne.0) afterj(j) = afterj(j) + 
     $              fld(i,j,k) * fc(i,j,k) * scp2(i,j)
            enddo               
         enddo                 
!$OMP END PARALLEL DO
         do j = 1, jj
            bfore = bfore + bforej(j)
            after = after + afterj(j)
         enddo       

      ENDDO               ! END of k-LOOP
c
c      if (bfore.ne.0.)
c     . write (lp,'(a,1p,3e14.6,e11.1)') 'fct3d conservation:',
c     .  bfore,after,after-bfore,(after-bfore)/bfore
      if (after.ne.0.) q = bfore / after
      write (lp,*) 'fct3d:',bfore,after,bfore/after,(after-bfore)/bfore
C       write (lp,'(a,f11.6)') 'fct3d: multiply fld field by',q
c
!$OMP PARALLEL DO
      do j = 1-margin, jj+margin
         do k = 1, kk
            do i = 1-margin, ii+margin
               if (ip(i,j).ne.0) then
 20               fld(i,j,k) = fld(i,j,k) * q
               endif            !ip
            enddo               !i
         enddo                  !k
      enddo                     !j
!$OMP END PARALLEL DO
c
c      print*,'Tracer after advection:',bfore,after,fld(itest,jtest,1)
c
      return
      end subroutine fct3d
c===================================================================
      end module mod_tracer
