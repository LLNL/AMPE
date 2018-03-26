define(NDIM,3)dnl
define(REAL,`double precision')dnl
include(SAMRAI_FORTDIR/pdat_m4arrdim3d.i)dnl

      subroutine tagcells3d(
     &  ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2,
     &  gcw0,gcw1,gcw2,
     &  tags,
     &  var, 
     &  refine_tag_val,
     &  tolerance,
     &  nequ)
c***********************************************************************
c***********************************************************************     
      implicit none
c input arrays:
      integer ifirst0,ilast0,ifirst1,ilast1,ifirst2,ilast2
      integer gcw0,gcw1,gcw2
      integer refine_tag_val
      integer nequ
      REAL    tolerance(0:nequ-1)
c
c variables in 3d cell indexed         
      integer 
     &     tags(CELL3d(ifirst,ilast,0))
      REAL
     &     var(CELL3dVECG(ifirst,ilast,gcw),0:nequ-1)
c
      integer ic0,ic1,ic2,ineq
c
c***********************************************************************     
      do ic2=ifirst2,ilast2
        do ic1=ifirst1,ilast1
          do ic0=ifirst0,ilast0
            do ineq=0,nequ-1
              tags(ic0,ic1,ic2) = 0
              if (var(ic0,ic1,ic2,ineq) .gt. tolerance(ineq)) 
     &            tags(ic0,ic1,ic2) = refine_tag_val
            end do
          end do
        end do
      end do

      return
      end
