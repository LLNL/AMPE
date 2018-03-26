define(NDIM,3)dnl
define(REAL,`double precision')dnl
define(PARAMDIR,SAMRAI_FORTDIR/../examples/LinAdv/fortran)dnl
include(SAMRAI_FORTDIR/pdat_m4arrdim3d.i)dnl
include(FORTDIR/m4bdry.i)dnl
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c bcceglobal is called before riemnv to set boundary values for left and
c   right states in riemann problem.
c
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c******************************************************************
      subroutine getbdry3d(bdry_type,
     &  ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2,
     &  ibeg0,iend0,ibeg1,iend1,ibeg2,iend2,
     &  ngc0,ngc1,ngc2,
     &  bside,
     &  dx,time,
     &  uval,
     &  bdry_states,bdry_data,bdry_case)
c***********************************************************************
      implicit none
include(PARAMDIR/probparams.i)dnl

      integer ngc0,ngc1,ngc2,bside,bdry_type
      integer ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2
      integer ibeg0,iend0,ibeg1,iend1,ibeg2,iend2
      integer bdry_case
      REAL bdry_states(2*NDIM)
      REAL bdry_data(2*NDIM)
      REAL time
c
      REAL dx(0:NDIM-1)
c
      REAL
     &  uval(CELL3dVECG(ifirst,ilast,ngc))

      REAL slope,uval0
      integer ic1,ic0,ic2,ict1,ict0,ict2
      integer sn
      integer ibdebug(0:NDIM-1),iedebug(0:NDIM-1)

c      write(6,*) "IN getbdry"
c      write(6,*) "bdry_type ",bdry_type,", bside ",bside
c      write(6,*) "ibeg ",ibeg0,ibeg1,ibeg2
c      write(6,*) "iend ",iend0,iend1,iend2
c      write(6,*) "ifirst ",ifirst0,ifirst1,ifirst2
c      write(6,*) "ilast ",ilast0,ilast1,ilast2
c      write(6,*) "bdry_case ",bdry_case

      ibdebug(0) = ibeg0
      ibdebug(1) = ibeg1
      ibdebug(2) = ibeg2
      iedebug(0) = iend0
      iedebug(1) = iend1
      iedebug(2) = iend2

      if (bdry_type.eq.1) then
c
c     bside      index     bside      index      bside      index
c       0     (-1, 0, 0)     2     ( 0,-1, 0)      4     ( 0, 0,-1)
c       1     ( 1, 0, 0)     3     ( 0, 1, 0)      5     ( 0, 0, 1)
c***********************************************************************
        sn = 1+bside
        if (bside.eq.0) then
c***********************************************************************
c x0 boundary
c***********************************************************************
c           ibdebug(0) =ifirst0-ngc0
c           iedebug(0) =ifirst0-1
do_bdry_face(0,1,2,`ifirst0-ngc0,ifirst0-1',`ifirst0',`ifirst0+1')dnl
        else if (bside.eq.1) then
c
c***********************************************************************
c x1 boundary
c***********************************************************************
c           ibdebug(0) =ilast0+1
c           iedebug(0) =ilast0+ngc0
do_bdry_face(0,1,2,`ilast0+1,ilast0+ngc0',`ilast0',`ilast0-1')dnl
        else if (bside.eq.2) then
c
c***********************************************************************
c y0 boundary
c***********************************************************************
c           ibdebug(1) =ifirst1-ngc1
c           iedebug(1) =ifirst1-1
do_bdry_face(1,0,2,`ifirst1-ngc1,ifirst1-1',`ifirst1',`ifirst1+1')dnl
        else if (bside.eq.3) then
c
c***********************************************************************
c y1 boundary
c***********************************************************************
c           ibdebug(1) =ilast1+1
c           iedebug(1) =ilast1+ngc1
do_bdry_face(1,0,2,`ilast1+1,ilast1+ngc1',`ilast1',`ilast1-1')dnl
        else if (bside.eq.4) then
c
c***********************************************************************
c z0 boundary
c***********************************************************************
c           ibdebug(2) =ifirst2-ngc2
c           iedebug(2) =ifirst2-1
do_bdry_face(2,0,1,`ifirst2-ngc2,ifirst2-1',`ifirst2',`ifirst2+1')dnl
        else if (bside.eq.5) then
c
c***********************************************************************
c z1 boundary
c***********************************************************************
c           ibdebug(2) =ilast2+1
c           iedebug(2) =ilast2+ngc2
do_bdry_face(2,0,1,`ilast2+1,ilast2+ngc2',`ilast2',`ilast2-1')dnl
        endif
c***********************************************************************

      else if (bdry_type.eq.2) then
c***********************************************************************
c        write(6,*) "doing edges"
c***********************************************************************
c
c     bside     index      bside     index      bside     index
c       0     (0,-1,-1)      4     (-1,0,-1)      8     (-1,-1,0)
c       1     (0, 1,-1)      5     (-1,0, 1)      9     ( 1,-1,0)
c       2     (0,-1, 1)      6     ( 1,0,-1)     10     (-1, 1,0)
c       3     (0, 1, 1)      7     ( 1,0, 1)     11     ( 1, 1,0)
c***********************************************************************
        sn = 1+bside
        if (bside.eq.0) then
c           ibdebug(1) =ifirst1-ngc1
c           iedebug(1) =ifirst1-1
c           ibdebug(2) =ifirst2-ngc1
c           iedebug(2) =ifirst2-1
do_bdry_edge(0,1,2,`ifirst1-ngc1,ifirst1-1',`ifirst2-ngc2,ifirst2-1',`ifirst1',`ifirst1+1')dnl
        else if (bside.eq.1) then
c           ibdebug(1) =ilast1+1
c           iedebug(1) =ilast1+ngc1
c           ibdebug(2) =ifirst2-ngc2
c           iedebug(2) =ifirst2-1
do_bdry_edge(0,1,2,`ilast1+1,ilast1+ngc1',`ifirst2-ngc2,ifirst2-1',`ilast1',`ilast1-1')dnl
        else if (bside.eq.2) then
c           ibdebug(1) =ifirst1-ngc1
c           iedebug(1) =ifirst1-1
c           ibdebug(2) =ilast2+1
c           iedebug(2) =ilast2+ngc2
do_bdry_edge(0,1,2,`ifirst1-ngc1,ifirst1-1',`ilast2+1,ilast2+ngc2',`ifirst1',`ifirst1+1')dnl
        else if (bside.eq.3) then
c           ibdebug(1) =ilast1+1
c           iedebug(1) =ilast1+ngc1
c           ibdebug(2) =ilast2+1
c           iedebug(2) =ilast2+ngc2
do_bdry_edge(0,1,2,`ilast1+1,ilast1+ngc1',`ilast2+1,ilast2+ngc2',`ilast1',`ilast1-1')dnl

        else if (bside.eq.4) then
c           ibdebug(0) =ifirst0-ngc0
c           iedebug(0) =ifirst0-1
c           ibdebug(2) =ifirst2-ngc2
c           iedebug(2) =ifirst2-1
do_bdry_edge(1,2,0,`ifirst2-ngc2,ifirst2-1',`ifirst0-ngc0,ifirst0-1',`ifirst2',`ifirst2+1')dnl
        else if (bside.eq.5) then
c           ibdebug(0) =ifirst0-ngc0
c           iedebug(0) =ifirst0-1
c           ibdebug(2) =ilast2+1
c           iedebug(2) =ilast2+ngc2
do_bdry_edge(1,2,0,`ilast2+1,ilast2+ngc2',`ifirst0-ngc0,ifirst0-1',`ilast2',`ilast2-1')dnl
        else if (bside.eq.6) then
c           ibdebug(0) =ilast0+1
c           iedebug(0) =ilast0+ngc0
c           ibdebug(2) =ifirst2-ngc2
c           iedebug(2) =ifirst2-1
do_bdry_edge(1,2,0,`ifirst2-ngc2,ifirst2-1',`ilast0+1,ilast0+ngc0',`ifirst2',`ifirst2+1')dnl
        else if (bside.eq.7) then
c           ibdebug(0) =ilast0+1
c           iedebug(0) =ilast0+ngc0
c           ibdebug(2) =ilast2+1
c           iedebug(2) =ilast2+ngc2
do_bdry_edge(1,2,0,`ilast2+1,ilast2+ngc2',`ilast0+1,ilast0+ngc0',`ilast2',`ilast2-1')dnl

        else if (bside.eq.8) then
c           ibdebug(0) =ifirst0-ngc0
c           iedebug(0) =ifirst0-1
c           ibdebug(1) =ifirst1-ngc1
c           iedebug(1) =ifirst1-1
do_bdry_edge(2,0,1,`ifirst0-ngc0,ifirst0-1',`ifirst1-ngc1,ifirst1-1',`ifirst0',`ifirst0+1')dnl
        else if (bside.eq.9) then
c           ibdebug(0) =ilast0+1
c           iedebug(0) =ilast0+ngc0
c           ibdebug(1) =ifirst1-ngc1
c           iedebug(1) =ifirst1-1
do_bdry_edge(2,0,1,`ilast0+1,ilast0+ngc0',`ifirst1-ngc1,ifirst1-1',`ilast0',`ilast0-1')dnl
        else if (bside.eq.10) then
c           ibdebug(0) =ifirst0-ngc0
c           iedebug(0) =ifirst0-1
c           ibdebug(1) =ilast1+1
c           iedebug(1) =ilast1+ngc1
do_bdry_edge(2,0,1,`ifirst0-ngc0,ifirst0-1',`ilast1+1,ilast1+ngc1',`ifirst0',`ifirst0+1')dnl
        else if (bside.eq.11) then
c           ibdebug(0) =ilast0+1
c           iedebug(0) =ilast0+ngc0
c           ibdebug(1) =ilast1+1
c           iedebug(1) =ilast1+ngc1
do_bdry_edge(2,0,1,`ilast0+1,ilast0+ngc0',`ilast1+1,ilast1+ngc1',`ilast0',`ilast0-1')dnl
        endif
c
c***********************************************************************

      else if (bdry_type.eq.NDIM) then
c***********************************************************************
c        write(6,*) "IN getbdrynode"
c***********************************************************************
c     bside     index       bside     index
c       0     (-1,-1,-1)      4     (-1,-1, 1)
c       1     ( 1,-1,-1)      5     ( 1,-1, 1)
c       2     (-1, 1,-1)      6     (-1, 1, 1)
c       3     ( 1, 1,-1)      7     ( 1, 1, 1)
c***********************************************************************
        sn = 1+bside

        if (bside.eq.0) then
c           ibdebug(0) =ifirst0-ngc0
c           iedebug(0) =ifirst0-1
c           ibdebug(1) =ifirst1-ngc1
c           iedebug(1) =ifirst1-1
c           ibdebug(2) =ifirst2-ngc2
c           iedebug(2) =ifirst2-1
do_bdry_node(`ifirst0-ngc0,ifirst0-1',`ifirst1-ngc1,ifirst1-1',`ifirst2-ngc2,ifirst2-1',`ifirst0',`ifirst0+1')dnl
        else if (bside.eq.1) then
c           ibdebug(0) =ilast0+1
c           iedebug(0) =ilast0+ngc0
c           ibdebug(1) =ifirst1-ngc1
c           iedebug(1) =ifirst1-1
c           ibdebug(2) =ifirst2-ngc2
c           iedebug(2) =ifirst2-1
do_bdry_node(`ilast0+1,ilast0+ngc0',`ifirst1-ngc1,ifirst1-1',`ifirst2-ngc2,ifirst2-1',`ilast0',`ilast0-1')dnl
        else if (bside.eq.2) then
c           ibdebug(0) =ifirst0-ngc0
c           iedebug(0) =ifirst0-1
c           ibdebug(1) =ilast1+1
c           iedebug(1) =ilast1+ngc1
c           ibdebug(2) =ifirst2-ngc2
c           iedebug(2) =ifirst2-1
do_bdry_node(`ifirst0-ngc0,ifirst0-1',`ilast1+1,ilast1+ngc1',`ifirst2-ngc2,ifirst2-1',`ifirst0',`ifirst0+1')dnl
        else if (bside.eq.3) then
c           ibdebug(0) =ilast0+1
c           iedebug(0) =ilast0+ngc0
c           ibdebug(1) =ilast1+1
c           iedebug(1) =ilast1+ngc1
c           ibdebug(2) =ifirst2-ngc2
c           iedebug(2) =ifirst2-1
do_bdry_node(`ilast0+1,ilast0+ngc0',`ilast1+1,ilast1+ngc1',`ifirst2-ngc2,ifirst2-1',`ilast0',`ilast0-1')dnl
        else if (bside.eq.4) then
c           ibdebug(0) =ifirst0-ngc0
c           iedebug(0) =ifirst0-1
c           ibdebug(1) =ifirst1-ngc1
c           iedebug(1) =ifirst1-1
c           ibdebug(2) =ilast2+1
c           iedebug(2) =ilast2+ngc2
do_bdry_node(`ifirst0-ngc0,ifirst0-1',`ifirst1-ngc1,ifirst1-1',`ilast2+1,ilast2+ngc2',`ifirst0',`ifirst0+1')dnl
        else if (bside.eq.5) then
c           ibdebug(0) =ilast0+1
c           iedebug(0) =ilast0+ngc0
c           ibdebug(1) =ifirst1-ngc1
c           iedebug(1) =ifirst1-1
c           ibdebug(2) =ilast2+1
c           iedebug(2) =ilast2+ngc2
do_bdry_node(`ilast0+1,ilast0+ngc0',`ifirst1-ngc1,ifirst1-1',`ilast2+1,ilast2+ngc2',`ilast0',`ilast0-1')dnl
        else if (bside.eq.6) then
c           ibdebug(0) =ifirst0-ngc0
c           iedebug(0) =ifirst0-1
c           ibdebug(1) =ilast1+1
c           iedebug(1) =ilast1+ngc1
c           ibdebug(2) =ilast2+1
c           iedebug(2) =ilast2+ngc2
do_bdry_node(`ifirst0-ngc0,ifirst0-1',`ilast1+1,ilast1+ngc1',`ilast2+1,ilast2+ngc2',`ifirst0',`ifirst0+1')dnl
        else if (bside.eq.7) then
c           ibdebug(0) =ilast0+1
c           iedebug(0) =ilast0+ngc0
c           ibdebug(1) =ilast1+1
c           iedebug(1) =ilast1+ngc1
c           ibdebug(2) =ilast2+1
c           iedebug(2) =ilast2+ngc2
do_bdry_node(`ilast0+1,ilast0+ngc0',`ilast1+1,ilast1+ngc1',`ilast2+1,ilast2+ngc2',`ilast0',`ilast0-1')dnl
        endif
      endif
c
c***********************************************************************
      return
      end
