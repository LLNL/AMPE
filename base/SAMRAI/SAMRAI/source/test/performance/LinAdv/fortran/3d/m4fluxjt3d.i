define(st_third,`dnl
      do ic$3=ifirst$3-1,ilast$3+1
         do ic$2=ifirst$2-1,ilast$2+1
           do ic$1=ifirst$1,ilast$1
             trnsvers=
     &           (flux$1(ic$1+1,$4)-flux$1(ic$1,$4))/(3*dx($1))
c
             st3_$1(ic0,ic1,ic2)=uval(ic0,ic1,ic2) -trnsvers
           enddo
         enddo
      enddo
')dnl
define(f_third,`dnl
c     write(6,*) "checking onedr sol in riemann solve  "
c     write(6,*) "         dt= ",dt
      do ic$3=ifirst$3-1,ilast$3+1
         do ic$2=ifirst$2,ilast$2
           do ic$1=ifirst$1,ilast$1+1
 
             if (advecspeed($1).ge.zero) then
               riemst = st3_$2($5)
             else
               riemst = st3_$2(ic0,ic1,ic2)
             endif
             flux$1$2(ic$1,$4)= dt*riemst*advecspeed($1)
           enddo
         enddo
      enddo
')dnl
define(correc_fluxjt,`dnl
c   correct the $1-direction with $2$3-fluxes
      do ic$3=ifirst$3,ilast$3
        do ic$2=ifirst$2,ilast$2
           ic$1=ifirst$1-1
           trnsvers=0.5*(
     &        (flux$2$3(ic$2+1,$4)-flux$2$3(ic$2,$4))/dx($2)+
     &        (flux$3$2(ic$3+1,$5)-flux$3$2(ic$3,$5))/dx($3))
  
           tracelft$1(ic$1+1,ic$2,ic$3)=tracelft$1(ic$1+1,ic$2,ic$3)
     &                                         - trnsvers
           do ic$1=ifirst$1,ilast$1
             trnsvers=0.5*(
     &          (flux$2$3(ic$2+1,$4)-flux$2$3(ic$2,$4))/dx($2)+
     &          (flux$3$2(ic$3+1,$5)-flux$3$2(ic$3,$5))/dx($3))
  
             tracelft$1(ic$1+1,ic$2,ic$3)=tracelft$1(ic$1+1,ic$2,ic$3)
     &                                           - trnsvers
               tracergt$1(ic$1  ,ic$2,ic$3)=tracergt$1(ic$1  ,ic$2,ic$3)
     &                                           - trnsvers
           enddo
           ic$1=ilast$1+1
           trnsvers=0.5*(
     &        (flux$2$3(ic$2+1,$4)-flux$2$3(ic$2,$4))/dx($2)+
     &        (flux$3$2(ic$3+1,$5)-flux$3$2(ic$3,$5))/dx($3))
  
           tracergt$1(ic$1  ,ic$2,ic$3)=tracergt$1(ic$1  ,ic$2,ic$3)
     &                                        - trnsvers
         enddo
      enddo
')dnl
