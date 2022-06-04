      program test
      implicit none
      character*120    k_dir, cStr, kflnm, kflnm_ful
      integer :: lp, iday, ihr, l, ii, jj, kk, nbdy
      integer :: i, j, k, ivar
      real, save, allocatable, dimension (:,:,:,:) :: ktensor
C 
      lp = 6
      ii = 1569
      jj = 1069
      kk = 30
      nbdy = 2
C 
      iday = 22
      ihr = 00
C       
      allocate( ktensor(1-nbdy:ii+nbdy,1-nbdy:jj+nbdy,kk,3) )
C 
c --- dir for params 
      k_dir = '/scratch/projects/ome/hra_expt/KL_adv_iso/sm101'
      cStr = '01020304' ! tracers used to calc K
      write(lp,*)
      write(lp,*) 'K will be read from ... ', k_dir, cStr
c
      kflnm = 'D***H**_GSH.dat'
      l = len_trim(kflnm)
      write(kflnm(l-13:l-11),'(i3.3)') iday
      write(kflnm(l-9:l-8),'(i2.2)') ihr


      kflnm_ful = trim(k_dir)//'/K_C'//trim(cStr)//'_'//trim(kflnm)
       open(unit=99, file=kflnm_ful, form='unformatted',
     &        status='old', action='read', access='stream', 
     &        convert='little_endian')
      read(99) ( ( ( (ktensor(i,j,k,ivar),j=1-nbdy,jj+nbdy), 
     &   i=1-nbdy,ii+nbdy ), k=1,kk ), ivar = 1,3 )
      close(99)

      write(lp,*) ktensor(500,500,20,1:3)

      end
