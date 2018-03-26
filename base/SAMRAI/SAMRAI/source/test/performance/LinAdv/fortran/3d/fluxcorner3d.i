include(FORTDIR/m4fluxjt.i)dnl

      subroutine onethirdstates3d(dt,dx,
     &  ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2,
     &  advecspeed,uval,
     &  flux0,flux1,flux2,
     &  st3_0,st3_1,st3_2)
c***********************************************************************
include(FORTDIR/const.i)dnl
c***********************************************************************     
c input arrays:
      integer ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2
      REAL dt 
c variables in 1d axis indexed
c
      REAL dx(0:NDIM-1)
c variables in 2d cell indexed         
      REAL
     &     advecspeed(0:NDIM-1),
     &     uval(CELL3d(ifirst,ilast,CELLG)),
     &     flux0(FACE3d0(ifirst,ilast,FLUXG)),
     &     flux1(FACE3d1(ifirst,ilast,FLUXG)),
     &     flux2(FACE3d2(ifirst,ilast,FLUXG)),
     &     st3_0(CELL3d(ifirst,ilast,CELLG)),
     &     st3_1(CELL3d(ifirst,ilast,CELLG)),
     &     st3_2(CELL3d(ifirst,ilast,CELLG))
c
c***********************************************************************     
c
      integer ic0,ic1,ic2
      REAL trnsvers
     
c     ******************************************************************
c     * complete tracing at cell edges
c     ******************************************************************
c         
st_third(0,1,2,`ic1,ic2')dnl
c
st_third(1,2,0,`ic2,ic0')dnl
c
st_third(2,0,1,`ic0,ic1')dnl
c
      return
      end   
c
c***********************************************************************
c***********************************************************************
c***********************************************************************
      subroutine fluxthird3d(dt,dx,
     &  ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2,
     &  advecspeed,
     &  uval,
     &  st3_0,st3_1,st3_2,
     &  flux01,flux12,flux20,
     &  flux02,flux10,flux21)
     
c***********************************************************************
include(FORTDIR/const.i)dnl
c***********************************************************************
c***********************************************************************
c input arrays:
      integer ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2
      REAL dt 
      REAL 
     &     dx(0:NDIM-1)
c variables in 2d cell indexed         
      REAL
     &     advecspeed(0:NDIM-1),
     &     uval(CELL3d(ifirst,ilast,CELLG))
c variables in 2d side indexed         
      REAL
     &     flux01(FACE3d0(ifirst,ilast,FLUXG)),
     &     flux10(FACE3d1(ifirst,ilast,FLUXG)),
     &     flux20(FACE3d2(ifirst,ilast,FLUXG)),
     &     flux02(FACE3d0(ifirst,ilast,FLUXG)),
     &     flux12(FACE3d1(ifirst,ilast,FLUXG)),
     &     flux21(FACE3d2(ifirst,ilast,FLUXG)),
     &       st3_0(CELL3d(ifirst,ilast,CELLG)),
     &       st3_1(CELL3d(ifirst,ilast,CELLG)),
     &       st3_2(CELL3d(ifirst,ilast,CELLG))
c
c***********************************************************************     
c
      integer ic0,ic1,ic2
      REAL   riemst
c
c***********************************************************************
c solve riemann problems for conservative flux
c  arguments: ( axis for RP, other axis, extra cells-direction)
c***********************************************************************
c      

f_third(0,1,2,`ic1,ic2',`ic0-1,ic1,ic2')dnl
c
f_third(0,2,1,`ic1,ic2',`ic0-1,ic1,ic2')dnl
c
f_third(1,0,2,`ic2,ic0',`ic0,ic1-1,ic2')dnl
c
f_third(1,2,0,`ic2,ic0',`ic0,ic1-1,ic2')dnl
c
f_third(2,0,1,`ic0,ic1',`ic0,ic1,ic2-1')dnl
c
f_third(2,1,0,`ic0,ic1',`ic0,ic1,ic2-1')dnl
c
c      call flush(6)     
      return
      end 
c***********************************************************************
c***********************************************************************
c***********************************************************************
      subroutine fluxcorrecjt3d(dt,dx,
     &  ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2,
     &  advecspeed,uval,
     &  flux01,flux12,flux20,
     &  flux02,flux10,flux21,
     &  tracelft0,tracelft1,tracelft2,
     &  tracergt0,tracergt1,tracergt2)
c***********************************************************************
include(FORTDIR/const.i)dnl
c***********************************************************************
c***********************************************************************     
c input arrays:
      integer ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2
      REAL dt 
c variables in 1d axis indexed
c
      REAL 
     &     dx(0:3-1)
c variables in 2d cell indexed         
      REAL
     &     advecspeed(0:NDIM-1),
     &     uval(CELL3d(ifirst,ilast,CELLG)),
     &     flux01(FACE3d0(ifirst,ilast,FLUXG)),
     &     flux10(FACE3d1(ifirst,ilast,FLUXG)),
     &     flux20(FACE3d2(ifirst,ilast,FLUXG)),
     &     flux02(FACE3d0(ifirst,ilast,FLUXG)),
     &     flux12(FACE3d1(ifirst,ilast,FLUXG)),
     &     flux21(FACE3d2(ifirst,ilast,FLUXG)),
     &     tracelft0(FACE3d0(ifirst,ilast,FACEG)),
     &     tracergt0(FACE3d0(ifirst,ilast,FACEG)),
     &     tracelft1(FACE3d1(ifirst,ilast,FACEG)),
     &     tracergt1(FACE3d1(ifirst,ilast,FACEG)),
     &     tracelft2(FACE3d2(ifirst,ilast,FACEG)),
     &     tracergt2(FACE3d2(ifirst,ilast,FACEG))
c
c***********************************************************************     
c
      integer ic0,ic1,ic2
      REAL trnsvers
c     REAL ttvlft,ttvrgt     
     
c     ******************************************************************
c     * complete tracing at cell edges
c     ******************************************************************
c         
correc_fluxjt(2,0,1,`ic1,ic2',`ic2,ic0')dnl
c
correc_fluxjt(1,2,0,`ic0,ic1',`ic1,ic2')dnl
c
correc_fluxjt(0,1,2,`ic2,ic0',`ic0,ic1')dnl
c
c      call flush(6)     
      return
      end   
c
c***********************************************************************
c***********************************************************************
c***********************************************************************
