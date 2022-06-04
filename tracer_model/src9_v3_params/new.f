c=====================================================================
c          Apply the parameterized diffusion & advection 
c=====================================================================
      subroutine tracer_diffadv(fld, klay)
c
c     Compute Laplacian operator & adv of input field
c       DIV_F (or ADV) = - K * del_del<c> + lambda * del<c>
c       LHS     -  [m/s*c], uh*c*L/L2 (or uh*c/L)
c       K       -  [m2/s]
c       lambda  -  [m/s], p-grid
c       del<c>  -  [c], multiplied by <h> or h
c       del2<c> -  [c/m], cgrad*L/L2
c
      implicit none
      integer margin, klay
      real, parameter :: aspmax = 2.0
      real, parameter :: epsil = 1.e-11
      real, parameter :: onemu = 1.e-12 !very small layer thickness (m)
      real harmon, a, b, kxx, kyy
      real, dimension(1-nbdy:ii+nbdy,1-nbdy:jj+nbdy) :: fld,dtdx,dtdy,
     $     xflux,yflux,flds
      real factor, kdel2c, lmdcgrad
c
      harmon(a,b) = 2.*a*b / (a+b)
c     
!OMP PARALLEL DO PRIVATE(i,j)
      do j = 1-nbdy, jj+nbdy
         do i = 1-nbdy, ii+nbdy
            dtdx(i,j) = 0.0
            dtdy(i,j) = 0.0
            xflux(i,j) = 0.0
            yflux(i,j) = 0.0
         enddo
      enddo
!OMP END PARALLEL DO
c
      margin = nbdy - 1
!OMP PARALLEL DO PRIVATE(i,j)
      do j = 1-nbdy, jj+nbdy
         do i = 1-nbdy, ii+nbdy
            flds(i,j) = fld(i,j)
         enddo
      enddo
!OMP END PARALLEL DO
c
c --- Spatial-mean tracer used in diffusion & advection operators! 
c
      call tracer_smooth(flds,klay,50,50)
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
c
c     Calculate K*DIV_cgrad - lmd*cgrad [c*m/s], and update tracer
c     When updating, compute dt*LHS/h ([c]), consistent with FLUX method
c      I.e. f = f + K * del_del<c> - lambda * del<c>
c      
!$OMP   PARALLEL DO PRIVATE(i,j,factor)
c!$OMP&           SCHEDULE(STATIC,jblk)
      do j = 1, jj      
         do i = 1, ii
            if (ip(i,j).ne.0) then
              ! K * DIV_cgrad [m2/s * cgrad*m/m2] 
              factor = Ktensor(i,j,klay,1) / scp2(i,j)
              kdel2c = factor * ( 
     &          dtdx(i+1,j)*scuy(i+1,j) - dtdx(i,j)*scuy(i,j) +
     &          dtdy(i,j+1)*scvx(i,j+1) - dtdy(i,j)*scvx(i,j) ) 
              ! lmd*cgrad [m/s*c]
              lmdcgrad = lmdu(i,j,klay) * 0.5*(dtdx(i,j) + dtdx(i+1,j))
     &          + lmdv(i,j,klay) * 0.5*(dtdy(i,j) + dtdy(i,j+1))
              ! update tracer [c]
              factor = delt1 / max(dpm(i,j,klay),onemu)
              fld(i,j) = fld(i,j) + factor * (kdel2c - lmdcgrad)
            endif               !ip
         enddo                  !i
      enddo                     !j
!$OMP   END PARALLEL DO
      return
      end subroutine tracer_diff

