      program tracers
      use mod_tracer ! HYCOM offline tracer 
      use mod_za     ! HYCOM array I/O interface
c   
      implicit none
c
c --- Offline tracer code.
c     Author: Z. Garraffo
c====================================================
c
c     This is version 3 that has buffer points, uses
c     vertical diffusion and can handle both velocities
c     and mass fluxes. 
c     
c     Author: I. Kamenkovich
c====================================================
c     
      character*200    cwd,flnm,out_path,dir_flnm_o
      character*200    k_dir, cStr, charK, kflnm, kflnm_ful
      logical          dir_exist
      integer          id, nvar
      integer          iexpt,yrflag,ka,l,iaux,iw,nstep
      integer          time_i,time_o,time_s,arch_l,iyear,iday,ihr
      double precision time_ai,time_ao,time_sp,arch_len,
     $     timetot,time1,time2
      real             time,delt, q1,q2,q0,t1,t2
      integer, parameter :: uflflg = 0   ! 'tracer_uvflx'
      integer, parameter :: forcflg = 0  ! 'getdat_archo'
      integer, parameter :: mforcflg = 1  ! 'getdat_archo'
      integer, parameter :: whatTracer = 2 ! 1: orig 10 c; 2: test
      integer :: corrflg
c
      real amin,amax
      integer iord, margin
C       
      real :: cm1, cm2
c
c --- Read the size of model grid, idm & jdm   
c     'xcspmd' in 'mod_xc.F'; 'zaiost' in 'mod_za.F'
c
      call xcspmd
      call zaiost
c
      lp = 6
c
c --- Set the output path
      write(lp,*)
      call getcwd(cwd) ! current dir
      out_path = trim(cwd)//'/out*/'
      
C     check if the dir exists
      do id = 1, 9
         l = len_trim(out_path)
         write(out_path(l-1:l-1),'(i1.1)') id
         ! 'directory' works with ifort, 'file' with gfort
         inquire(directory = out_path, exist = dir_exist)
         if ( dir_exist ) then
            write(lp,*) 'Dir exists ... ', out_path
         else
            exit
         endif
      enddo
      call system('mkdir ' // trim(out_path))
      write(lp,*) 'OUTPUT to: ', out_path
      write(lp,*)
c=====================================================================
c                      READ IN PARAMETERS 
c=====================================================================
c --- 'yrflag' = days in year flag (0=360, 1=366, 2=366J1, 3=actual, 4=365)
c --- 'time_s' = start of archive (6939 for 12/31/1919) or 0 (no real dates)
c --- 'arch_l' = length of archive (forcing) in days (for recycling in multi-year runs) 
c --- 'time_i' = initial day for the tracer integration
c --- 'time_o' = final day for the tracer integration
c --- 'archfq' = intervals between archive input (hours) 12.0
c --- 'archin' = time step for the tracer (hours) 1.0
c --- 'ii    ' = longitudinal domain size (excluding buffer points)
c --- 'jj    ' = latitudinal  domain size (excluding buffer points)
c --- 'kk    ' = actual number of top layers used
c --- 'idm   ' = longitudinal array size
c --- 'jdm   ' = latitudinal  array size
c --- 'kdm   ' = total number of layers
c
      call blkini(yrflag, 'yrflag')
      call blkini(time_s, 'time_s')
      call blkini(arch_l, 'arch_l')
      call blkini(time_i, 'time_i')
      call blkini(time_o, 'time_o')
      call blkinr(archfq,
     &     'archfq','("blkinr: ",a6," =",f11.4," hours")')
      call blkinr(archin,
     &           'archin','("blkinr: ",a6," =",f11.4," hours")')
c
c --- Times for the whole integration     
c
      time_sp = time_s
      arch_len = arch_l
      time_ai = time_i - 1.d0
      time_ao = time_o
      write(lp,*) 'First day in archive:',time_sp,' Archive length:',
     $     arch_len,' Integrating from ',time_ai,' to ',time_ao
c
c --- Units of time steps of archive & integration converted to [days]
c
      archfq = archfq / 24.
      archin = archin / 24.
      if (archin.gt.archfq) then
         write(lp,*)'error, archin can not be greater than archfq'
         stop
      endif
c     # of time steps of interp within one archive interval (integer)
      nintrp = nint(archfq / archin)
      if (archfq/archin.ne.float(nintrp)) then
         stop 'archfq must be an integer multiple of archin'
      endif
c     time steps in [seconds]
      delt = archfq * 86400.
      delt1 = archin * 86400.
      write(lp,*) 'delt as archfq in seconds', delt
c
c --- Other params set outside
c
      call blkini(ii, 'ii    ')
      call blkini(jj, 'jj    ')
      if (ii.ne.idm-2*nbdy .or. jj.ne.jdm-2*nbdy) then
         write(lp,*)
         write(lp,*) 'ii and jj (',ii,jj,') are inconsistent with'//
     $        'idm, jdm and nbdy (',idm,jdm,nbdy,')'
         write(lp,*)
         call flush(lp)
         stop 
      endif
      call blkini(kdm, 'kdm   ')
      call blkini(kk,  'kk    ')
      jblk = (jj + 2*nbdy)
c --- 'thbase'   = reference sigma
c --- 'nforf'    = number of daily-mean forcing fields
c --- 'ntracr'   = number of tracers
c --- 'tracin'   = 0: initialze to 0 
c                  1: from restart_tracr_in
c                  2: from relax.tracer.a,b
c --- 'mxlflg'   = 0: no extra mixing in mixed layer; 
c                  1: instantaneous homog of mixed layer.
c --- 'ghtflg'   = 0: no non-local mixing
c                  1: nonlocal mixing (ghat variable)
c --- 'trcflg'   = 0: no surface conditions
c                  1: modified ideal age tracer (= time at the surface)
c                  2: surface "pulse" (= unity in a limited domain/time)
c                  3: lateral "pulse" (= unity in a limited domain/time)
c --- 'difflg'   = 0: scalar constant diffusivity
c                  1: scalar space-dependent diffusivity
c                  2: space-dependent diffusivity tensor
c --- 'forflg'   = 0 no forcing (internal sources/sinks)
c --- 'timflg'   = 0: regular time forcing (unfiltered)
c                  1: smoothed archives ("mean" instead of year)
c                  2: annual mean archives (const fields)
c                  3: combination of smoothed and annual-mean archives
c --- 'diagfq'   = intervals between tracer diagnostic output (hours)
c --- 'iord  '   = 1: simple donor cell
c                  2: complete with antidiffusive fluxes
c --- 'itest ' = i for printout
c --- 'jtest ' = j for printout
      call blkinr(thbase,
     &     'thbase','("blkinr: ",a6," =",f11.4," days")')
      call blkini(nforf , 'nforf ')
      call blkini(ntracr, 'ntracr')
      call blkini(tracin, 'tracin')
      call blkini(mxlflg, 'mxlflg')
      call blkini(ghtflg, 'ghtflg')
      call blkini(trcflg, 'trcflg')
      call blkini(difflg, 'difflg')
      call blkini(forflg, 'forflg')
      call blkini(timflg, 'timflg')
      call blkinr(diagfq,
     &     'diagfq','("blkinr: ",a6," =",f11.4," hours")')
c     output interval in [days]
      diagfq = diagfq / 24.
      call blkini(iord , 'iord  ')
      call blkini(itest, 'itest ')
      call blkini(jtest, 'jtest ')
      write(lp,*) '-----------------------------------'
c=====================================================================
c                      Initial time step
c=====================================================================
c
c --- array allocation
c
      call tracer_alloc
c
c --- read and derive grids, topography and land masks
c
      call geopar
      write(lp,*)
c
c --- first and last model days from the input archives.
c
      write(lp,*) '-----------------------------------'
      write(lp,*) 'archfq [days] = ',archfq
c
c --- # of archive records needed
c
      narch = nint((time_ao - time_ai) / archfq) 
c
c --- INITIALIZE the tracer or read a restart file
c     'time' is only used here
c     
      time = time_ai + time_sp
      if (whatTracer .eq. 1) then
         call tracer_init(time)
      elseif (whatTracer .eq. 2) then
         call tracer_init_ctest(time)
      endif
c
      flnm = 'uvdp_offline.day_hr.a'
c
c --- each time step needs 2 instantaneous archives for delta dp and
c --- one set of time-mean fields mean
c 
c --- first time step, we read the dpii at the beginning of the segment
      time2 = mod(time_ai-1.d0, arch_len) + 1.d0 + time_sp
      if (time2.eq.0.0) time2 = arch_len
      if (yrflag.eq.4) then
         iyear =  int((time2 - 1.d0)/365.d0) + 1
         iday  =  mod( time2 - 1.d0 ,365.d0) + 1
         ihr   =  int((time2 - int(time2)) * 24.d0)
         l = len_trim(flnm)
         write(flnm(l-7:l-5),'(i3.3)') iday
         write(flnm(l-3:l-2),'(i2.2)') ihr
      else
         call getdat_name(flnm,time2,yrflag,iyear)
      endif
c
c --- READ forcing flds (u/vflx, dpm, dpio). flnm show in 'getdat_archo'
c
      ! read full forcing
      call getdat_archo(flnm, mforcflg) 
      call flush(lp)
!$OMP PARALLEL DO PRIVATE(i,j,k)
      do k = 1, kk
         do j = 1-nbdy, jj+nbdy
            do i = 1-nbdy, ii+nbdy
               dpii(i,j,k) = dpio(i,j,k)

               dpcons(i,j,k) = dpio(i,j,k) 
            enddo
         enddo
      enddo
!$OMP END PARALLEL DO
      write(lp,*) '-----------------------------------'
c=====================================================================
c                      MAIN LOOP of integration
c=====================================================================
      do iar = 1, narch
         if (iar.eq.1 .or. timflg.ne.2) then
c
c ---       Read dp and all mean variables from the same forcing file
c           time2 is the day of the tracers run, e.g. 1.5, 2.0...
c  
            time2 = mod(time_ai+iar*archfq-1.d0,arch_len)+1.d0+time_sp
            if (yrflag.eq.4) then
               iyear =  int((time2 - 1.d0)/365.d0) + 1
               iday  =  mod( time2 - 1.d0 ,365.d0) + 1
               ihr   =  int((time2 - int(time2)) * 24.d0)
               l = len_trim(flnm)
               write(flnm(l-7:l-5),'(i3.3)') iday
               write(flnm(l-3:l-2),'(i2.2)') ihr
            else
               call getdat_name(flnm,time2,yrflag,iyear)
            endif
c           READ full forc (u/vflx, dpm, dpio). flnm will display
            write(lp,*) '---' 
            call getdat_archo(flnm, mforcflg)
            write(lp,'(a,i4,f10.2,a,a)') 'FIELDS at iar,time2,flnm = ',
     &           iar,time2,' ',flnm(1:len_trim(flnm))
            call flush(lp)
            !  isopycnal mass flux [m3], u*h multipllied by delt1*L
            q0 = delt1
            call tracer_uvflx(q0,uflflg)            
c
c set lateral u*h*dt to zero
!$OMP PARALLEL DO PRIVATE(i,j,k)
               do k = 1, kk
                  do j = 1-nbdy, jj+nbdy
                     do i = 1-nbdy, ii+nbdy
C                         uflx(i,j,k) = 0.0
C                         vflx(i,j,k) = 0.0
C                         dpio(i,j,k) = dpcons(i,j,k)
C                         dpm(i,j,k) = dpcons(i,j,k)
                     enddo
                  enddo
               enddo
!$OMP END PARALLEL DO

c
c ---       Calc MEAN & EDDY horizontal mass divergence [m]
c
            margin = nbdy - 1
!$OMP PARALLEL DO PRIVATE(i,j,k)
            do k = 1, kk
               do j = 1-margin, jj+margin
                  do i = 1-margin, ii+margin
                     if(ip(i,j).eq.1) then
                        hordiv(i,j,k) = uflx(i+1,j,k) - uflx(i,j,k)
     &                       + vflx(i,j+1,k) - vflx(i,j,k)
                        hordiv(i,j,k) = hordiv(i,j,k) * scp2i(i,j)
                     endif
                  enddo
               enddo
            enddo
!$OMP END PARALLEL DO
c     
c ---       Calc MEAN vertical mass div [m] and fluxes through layer 
c           interfaces (diaflx, [m]) in the  interval 'delt1'
c           'diaflx' is positive downward!!
c
!$OMP PARALLEL DO PRIVATE(i,j,k,ka,vertdiv) 
            do j = 1-margin, jj+margin
               do i = 1-margin, ii+margin
                  if(ip(i,j).ne.0) then
                     diaflx(i,j,1) = 0
                     do k = 1, kk
                        ka = max(1,k-1)
                        vertdiv = - hordiv(i,j,k) + ( dpii(i,j,k) 
     $                       - dpio(i,j,k) ) / float(nintrp)
                        ! set 'LHS = 0' to turn off dia-flx of tracer
                        diaflx(i,j,k) = diaflx(i,j,ka) + vertdiv !0
                     enddo ! k
                  endif
               enddo
            enddo
!$OMP END PARALLEL DO
c            write(lp,*) diaflx(itest,jtest,10)     
         endif       ! if (iar.eq.1 .or.timflg.ne.2); load new data
c
c ---   READ diffusivity tensor in all layers at current time (if needed)
c
         if (difflg.ge.1) then
            kflnm = 'D***H**_GSH.dat'
            l = len_trim(kflnm)
            write(kflnm(l-13:l-11),'(i3.3)') iday
            write(kflnm(l-9:l-8),'(i2.2)') ihr
C 
            if (difflg.eq.1) then !
               nvar = 1
               cStr = '0103'
               k_dir=trim('/scratch/projects/ome/hra_expt/')
     &         //'K_iso_div/sm101_uecs_uece_noextraSM_99_01_C'
     &         //trim(cStr)//'_noextraK'
            elseif (difflg.eq.2) then ! tensor
               nvar = 4
               cStr = '0103' 
               k_dir=trim('/projects2/rsmas/ikamenkovich/yxl1496/')
     &         //'GSH_OUTPUT/K_KL/'         
     &         //'KisoA_div/sm101_uecs_uece_noextraSM_99_01_C'
     &         //trim(cStr)//'_extraK'
            elseif (difflg.eq.3) then ! k-Chi
               nvar = 3
               cStr = '0103050709' ! 0103050709
               k_dir=trim('/projects2/rsmas/ikamenkovich/yxl1496/')
     &         //'GSH_OUTPUT/K_KL/'
     &         //'KL_adv_iso/sm101_uecs_uece_noextraSM_99_01_C'
     &         //trim(cStr)//'_extraK'
            endif
C             
            kflnm_ful=trim(k_dir)//'/K_'//trim(kflnm)
            write(lp,*) 'READ K from: ', kflnm_ful
            open(unit=99, file=kflnm_ful, form='unformatted',
     &         status='old', action='read', access='stream', 
     &         convert='little_endian')
            read(99) ( ( ( (ktensor(i,j,k,id),j=1-nbdy,jj+nbdy), 
     &         i=1-nbdy,ii+nbdy ), k=1,kk ), id = 1,nvar )
            close(99)
         endif
      
c        ---------------------------------------------------------
c                      interp to update 'tracer' 
c        ---------------------------------------------------------
         write(lp,*) 'Start the main loop ...'
c ---    START TIME STEPPING WITHIN AN ARCHIVE SEGMENT (time-interp)
         do nstep = 1, nintrp
c
c ---       Keep fluxes constant during the entire segment (.5 day) and 
c           linearly interpolate dpi to each subsegment
c
            q1 = float(nstep-1) / float(nintrp)
            q2 = float(nstep)   / float(nintrp)
!$OMP PARALLEL DO PRIVATE(i,j,k)
            do j = 1-nbdy, jj+nbdy
               do i = 1-nbdy, ii+nbdy
                  if (ip(i,j).ne.0) then
                     do k = 1, kk
                        dp(i,j,k,1)=dpii(i,j,k)*(1.-q1)+dpio(i,j,k)*q1
                        dp(i,j,k,2)=dpii(i,j,k)*(1.-q2)+dpio(i,j,k)*q2
                     enddo
                  endif
               enddo
            enddo
!$OMP END PARALLEL DO
c
c ---       SOLVE TRACER EQUATION
c
            do ktr = 1, ntracr
               write(lp,'(a,i2)') '--C-',ktr
C --- Record the 3D tracer before being advected
!$OMP PARALLEL DO PRIVATE(i,j,k)
               do k = 1, kk
                  do j = 1-nbdy, jj+nbdy
                     do i = 1-nbdy, ii+nbdy
                        tracer_prev(i,j,k) = tracer(i,j,k,ktr)
                     enddo
                  enddo
               enddo
!$OMP END PARALLEL DO

c --- 3D Rainer advection
               call fct3d(iord,tracer(1-nbdy,1-nbdy,1,ktr),
     &              dp(1-nbdy,1-nbdy,1,1),dp(1-nbdy,1-nbdy,1,2))
C 
               call tracer_mass(tracer(1-nbdy,1-nbdy,1,ktr),0,cm1)

c --- Diffusion (parameterization of lateral eddies)
               do k = 1, kk
                  if (difflg.lt.3) then ! 0 1 2
                     call tracer_diff   (tracer(1-nbdy,1-nbdy,k,ktr),
     &                  tracer_prev(1-nbdy,1-nbdy,k), k, difflg)
                  elseif (difflg.ge.3) then ! 3 ...
                     call tracer_diffadv(tracer(1-nbdy,1-nbdy,k,ktr),
     &                  tracer_prev(1-nbdy,1-nbdy,k), k)
                  endif
               enddo
C                
               call tracer_mass(tracer(1-nbdy,1-nbdy,1,ktr),0,cm2)
               if (difflg.ge.1) then ! 1 kiso(x,y); 2 Ktens; 3 k-Chi
                  corrflg = 1      
               else
                  corrflg = 0      
               endif
               write(lp,*) 'DIFF delta cm (bf, af, bf/af, af/bf-1):',
     &         cm1, cm2, cm1/cm2, cm2/cm1-1.0

C     Correct the delta c-mass due to diffusion?
               if (corrflg.eq.1)  then
                  write(lp,*) 'CORR by delta cm due to DIFF!!'
!$OMP PARALLEL DO PRIVATE(i,j,k)
                 do k = 1, kk
                   do j = 1-nbdy, jj+nbdy
                    do i = 1-nbdy, ii+nbdy
                     if (ip(i,j).ne.0) then
                           tracer(i,j,k,ktr) = tracer(i,j,k,ktr)*cm1/cm2
                     endif
                    enddo
                   enddo
                 enddo
!$OMP END PARALLEL DO
               else
                  write(lp,*) 'NO corr.'
               endif ! corrflg 

c     (timetot is the time since the start of this run)
               timetot = time_ai + iar*archfq + nstep*archin
               do k = 1, kk
c --- calculate buffer values
                  call tracer_buffer(tracer(1-nbdy,1-nbdy,k,ktr))

c --- SMALL-SCALE DIFFUSION, if not done above
                  if (difflg.ne.0)  call tracer_diff(tracer(1-nbdy,
     &               1-nbdy,k,ktr),tracer(1-nbdy,1-nbdy,k,ktr),k,0)

c --- 3-point square smoothing to suppress grid-scale noise
c                  call tracer_smooth(tracer(1-nbdy,1-nbdy,k,ktr),1,1,k)
c --- LATERAL BOUNDARY CONDITIONS
                  call tracer_lbc(tracer(1-nbdy,1-nbdy,k,ktr),timetot,k)
               enddo            ! kk

c --- SURFACE BOUNDARY CONDITIONS
               call tracer_sbc(tracer(1-nbdy,1-nbdy,1,ktr),timetot)
c --- VERTICAL DIFFUSION
!$OMP PARALLEL DO PRIVATE(i,j)
!$OMP&         SCHEDULE(STATIC,jblk)
               do j = 1, jj
                  do i = 1, ii
                     if (ip(i,j).ne.0) then
c                        call tracer_kpp(i,j)
                        if(mxlflg.eq.1) call tracer_mlhom(i,j)
                     endif      !ip
                  enddo         !i
               enddo            !j
!$OMP END PARALLEL DO
c --- re-calculate buffer values
               do k = 1, kk
                  call tracer_buffer(tracer(1-nbdy,1-nbdy,k,ktr))
               enddo
            enddo               ! ktr
         enddo                  ! nintrp
c --- update 'dpio'
!$OMP PARALLEL DO PRIVATE(i,j,k)
         do k = 1, kdm
            do j = 1-nbdy, jj+nbdy
               do i = 1-nbdy, ii+nbdy
                  dpii(i,j,k) = dpio(i,j,k)
               enddo
            enddo
         enddo
!$OMP END PARALLEL DO

c     -----------------------------------------------------------------
         t1 = iar
         t2 = diagfq
c     --- save diagnostic file 
          write(lp,*) 'iar = ', t1, 'diagfq = ', t2
         if (mod(t1,2*t2).eq.0.0) then
c  --- set the name of the file; time1 is the total time since the start
            time1 = time_ai + iar*archfq
            l = len_trim(flnm) 
c           flnm_o='/projects2/rsmas/ikamenkovich/HYCOM_OUTPUT/'//
c     $           'Atl_offline/tracer_diag_year_day'
            flnm_o = 'tracer_diag_year_day_hr' 
            iyear  =  int((time1 - 1.d0)/365.d0) + 1
            iday   =  mod( time1 - 1.d0 ,365.d0) + 1
            l = len_trim(flnm_o)
            write(lp,*)  'time to write output', iyear, iday, ihr
            write(flnm_o(l-10:l-7),'(i4.4)') iyear
            write(flnm_o(l-5:l-3),'(i3.3)') iday
            write(flnm_o(l-1:l),'(i2.2)') ihr
            write(lp,'(a,i4,f10.2,a)') 'diagnostic file at ',
     &           iar, time2,' '// flnm_o(1:len_trim(flnm_o))
c           filename with dir
            dir_flnm_o = trim(out_path)//flnm_o
            call tracer_output(time1,dir_flnm_o)
         endif
      enddo                     ! iar
c
      stop
      end
