!     ifort -g -O3 -convert big_endian -openmp -o calc_mean_uvdp_HRA.f
      program calc_mean_uvdp_HRA
!     
!     The program reads the variables needed for the offline model
!     and saves them in the forcing files.

!     The forcing fields are calc by time-averaging and spatial filtering
!     'make_uvdp_HRA.f' only performs time-averaging.
!
      implicit none
      character*4 var_name
c
      !var_name = 'dpii'
      !call process_var(var_name)
      var_name='uflx' 
      call process_var(var_name)
      var_name='vflx'
      call process_var(var_name)
      !var_name='dpmm'
      !call process_var(var_name)
      stop
      end
c
c=====================================================================
c          SUBTOUTINE: process_var
c=====================================================================
      subroutine process_var(var_name)
      implicit none
      integer nday,i,j,k,days,daye
      integer nrec,npad,nrecl
      character*3 dayn
      character*2 hourn
      integer, parameter :: isz = 1573, jsz = 1073, ksz = 30, nvar = 3
      integer, parameter :: nfilx = 50, nfily = 50
      real, dimension(isz,jsz) :: var, sc2, scux, scuy, scvx,
     $     scvy, scpx, scpy
      real, allocatable :: pad(:)
      character*4 var_name
      character*100 foc_path, grd_path, sav_path
      integer lp
c
      npad=4096-mod(isz*jsz,4096)
      allocate (pad(npad))
      inquire(iolength=nrecl)var,pad
c
      grd_path = '/nethome/yxl1496/hycom/GS_HR/'
c     grd_path = '/scratch/projects/ome/hra_expt/'
      foc_path = '/scratch/projects/ome/hra_expt/UVDP/'
      sav_path = '/scratch/projects/ome/hra_expt/UVDP_CS3/'
c
      lp = 6
      days = 001
      daye = 003
c
c     ================ Read grid (fid=31)
c
      open(31,file=trim(grd_path)//'regional.grid.a', access='direct',
     $     recl=nrecl, form='unformatted', status='old',
     $     convert='big_endian')
      write(lp,*)  'READING '//trim(grd_path)//'regional.grid.a'
c
      read(31,rec=10) scpx
      read(31,rec=11) scpy
      read(31,rec=14) scux
      read(31,rec=15) scuy
      read(31,rec=16) scvx
      read(31,rec=17) scvy
      close(31)
c
c     ================== Read UVDP (fid=21)
c
      print*,'PROCESSING ',var_name 
c
c     ---- loop over time
c
      DO nday = 2*days-1, 2*daye

        write(dayn,'(i3.3)') (nday+1)/2 
        write(hourn,'(i2.2)') 12*(nday+1-2*((nday+1)/2))
c
c       ---- open file to be READ
        write(lp,*) 'READING '//trim(foc_path)//'uvdp_offline.'//dayn//
     $  '_'//hourn//'.a'
        open(21, file=trim(foc_path)//'uvdp_offline.'//dayn//'_'//
     $  hourn//'.a', access='direct', recl=nrecl, convert='big_endian',
     $  form='unformatted',status='unknown')
c
c       ---- open file to be SAVED
        write(lp,*) 'TO SAVE '//trim(sav_path)//'uvdp_offline.'//dayn//
     $    '_'//hourn//'.a' 
        open(12,file=trim(sav_path)//'uvdp_offline.'//dayn//'_'//
     $  hourn//'.a', access='direct', recl=nrecl, convert='big_endian',
     $  form='unformatted', status='unknown')
c
c       ---- layer by layer
c
!OMP PARALLEL DO PRIVATE(k)
        do k = 1,ksz
          print*,'K level ', k
c
c         ---- Read 'dpm', which is needed 
c
          if (var_name.eq.'dpii') nrec = k
          if (var_name.eq.'uflx') nrec = ksz + nvar*(k-1) + 1
          if (var_name.eq.'vflx') nrec = ksz + nvar*(k-1) + 2
          if (var_name.eq.'dpmm') nrec = ksz + nvar*(k-1) + 3
          read(21,rec=nrec) var
c
c         ---- Spatially running average and save (fid=12)
c
          if (var_name.eq.'uflx') then
            sc2 = scux*scuy
          elseif (var_name.eq.'vflx') then
            sc2 = scvx*scvy
          else
            sc2 = scpx*scpy
          endif
          call coarsen(var, nfilx, nfily, sc2, k)  
          write(12,rec=nrec) ((var(i,j),i=1,isz),j=1,jsz)
        enddo   ! end of k-loop
!OMP END PARALLEL DO
         close(12)
         close(21)
      ENDDO       ! end of nday-loop
      return
      end
c=====================================================================
c          SUBTOUTINE: spatially coarsen the field
c===================================================================== 
      subroutine coarsen(fld, nfili, nfilj, sc2, klay)
c
c --- Smooth the 2D field using a boxcar filter
c
      implicit none
      integer :: nfili, nfilj
      integer :: is, ie, js, je, ic, jc, i, j
      integer, parameter :: isz = 1573, jsz = 1073 
      real, dimension(isz,jsz) :: fld, flds, arsm, sc2
      real, dimension(jsz) :: fsum_oj, fsum_nj
      real :: fsum_o, fsum_n, fldij, areaij
C 
      real, parameter :: onemu = 1.e-6
      integer, intent(in) :: klay

!OMP PARALLEL DO PRIVATE(i,j)
         do i = 1, isz
            do j = 1, jsz
              ! sc2(i,j) = 1
              ! fld(i,j) = abs(fld(i,j)) 
            enddo
         enddo
!OMP END PARALLEL DO

c
c     Calculate sums over j-direction at each point
c
!OMP PARALLEL DO PRIVATE(i,j,jc,js,je)
      do j = 1, jsz
         fsum_oj(j) = 0.0
         js = max(1,  j-nfilj)
         je = min(jsz,j+nfilj)
         do i = 1, isz
            ! tracermass before smooth 
            fsum_oj(j) = fsum_oj(j) + fld(i,j) * sc2(i,j)
            !
            flds(i,j) = 0.0
            arsm(i,j) = 0.0
            ! sum over current meridional segment, ignoring vanishing 
            ! or land points in the sum. But note these points can have 
            ! non-zero sums
             do jc = js, je
                if (abs(fld(i,jc)) .gt. onemu) then
                   flds(i,j) = flds(i,j) + fld(i,jc) * sc2(i,jc)
                   arsm(i,j) = arsm(i,j) + sc2(i,jc) 
                endif
             enddo
         enddo
      enddo
!OMP END PARALLEL DO
c
c     Sum j-sums in the i-direction and divide by area
c
!OMP PARALLEL DO PRIVATE(i,j,ic,is,ie)
      do j = 1, jsz
         fsum_nj(j) = 0.0
         do i = 1, isz
            is = max(1,  i-nfili)
            ie = min(isz,i+nfili)
            ! sum the y-sum (flds, arsm) within the box. The sum at the
            ! vanishing/land points are 0's 
            fldij = 0.0
            areaij = 0.0
            if (abs(fld(i,j)) .gt. onemu) then
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
            fsum_nj(j) = fsum_nj(j) + fld(i,j) * sc2(i,j)
         enddo
      enddo
!OMP END PARALLEL DO
c
c     Correct the smoothed field, if needed
c
      fsum_o = 0.0
      fsum_n = 0.0
!OMP PARALLEL DO PRIVATE(j)
      do j = 1, jsz
         fsum_o = fsum_o + fsum_oj(j)
         fsum_n = fsum_n + fsum_nj(j)
      enddo
!OMP END PARALLEL DO
      if (klay .eq. klay) then
           write(*,'(a,e12.3,a,i2,a,i3,a,i3)'),'Delta Cmass for SM: ' 
     &    ,fsum_n/fsum_o-1.0, ' Z', klay, ' ni=',nfili, ' nj=',nfilj
!OMP PARALLEL DO PRIVATE(i,j)
         do i = 1, isz
            do j = 1, jsz
c               fld(i,j) = fld(i,j) * fsum_o/fsum_n
            enddo
         enddo
!OMP END PARALLEL DO
      endif
c
      return
      end
