define(trace_init,`dnl
      do  ic$3=ifirst$3-FACEG,ilast$3+FACEG
         do  ic$2=ifirst$2-FACEG,ilast$2+FACEG
           ie$1=ifirst$1-FACEG
           tracergt$1(ie$1,ic$2,ic$3)=uval($5)
           tracelft$1(ie$1,ic$2,ic$3)=uval($5)

           do  ie$1=ifirst$1+1-FACEG,ilast$1+FACEG
             tracelft$1(ie$1,ic$2,ic$3)=uval($4)
             tracergt$1(ie$1,ic$2,ic$3)=uval($5)
           enddo

           ie$1=ilast$1+FACEG+1
           tracelft$1(ie$1,ic$2,ic$3)=uval($4)
           tracergt$1(ie$1,ic$2,ic$3)=uval($4)
         enddo
      enddo

 
')dnl
define(trace_call,`dnl
        do ic$3=ifirst$3-2,ilast$3+2
          do ic$2=ifirst$2-2,ilast$2+2
              do ic$1=ifirst$1-FACEG,ilast$1+FACEG
                ttraclft(ic$1) = tracelft(ic$1,ic$2,ic$3)
                ttracrgt(ic$1) = tracergt(ic$1,ic$2,ic$3)
              enddo
   
            call trace(dt,ifirst$1,ilast$1,mc,
     &        dx,idir,advecspeed,igdnv,
     &        ttraclft,ttracrgt,
     &        ttcelslp,ttedgslp)
            do ic$1=ifirst$1-FACEG,ilast$1+FACEG
                tracelft(ic$1,ic$2,ic$3) = ttraclft(ic$1)
                tracergt(ic$1,ic$2,ic$3) = ttracrgt(ic$1)
            enddo
          enddo
        enddo
')dnl
