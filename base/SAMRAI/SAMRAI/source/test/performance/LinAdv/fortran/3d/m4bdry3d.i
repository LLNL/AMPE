define(cell_even,`dnl
                  uval(ic0,ic1,ic2)=uval($1)
')dnl
define(lin_boundary,`dnl
                  ict$1 = $2
                  uval0 = uval(ict0,ict1,ict2)
                  ict$1 = $3
                  slope=(uval0-uval(ict0,ict1,ict2))/($2-$3)
                  uval(ic0,ic1,ic2)=uval0 + (ic$1-$2)*slope
')dnl
define(fixed_state,`dnl
                  uval(ic0,ic1,ic2)=bdry_states(1+(sn-1))
')dnl
define(symmetry,`dnl
                  ict$1 = 3*$2-$3-ic$1
                  uval(ic0,ic1,ic2)=uval(ict0,ict1,ict2)
')dnl
define(do_bdry_face,`dnl
c          write(6,*) "ibdebug ",ibdebug
c          write(6,*) "iedebug ",iedebug
           if ((bdry_case.eq.FLOW)) then
             do ic$2=ibeg$2,iend$2
               ict$2 = ic$2
               do ic$3=ibeg$3,iend$3
                 ict$3 = ic$3
                 do ic$1=$4
                   ict$1 = $5
cell_even(`ict0,ict1,ict2')dnl
                 enddo
               enddo
             enddo
           else if (bdry_case.eq.LINEAR) then
             do ic$2=ibeg$2,iend$2
                ict$2 = ic$2
                do ic$3=ibeg$3,iend$3
                   ict$3 = ic$3
                   do ic$1=$4
lin_boundary(`$1',`$5',`$6')dnl
                   enddo
                enddo
             enddo
           else if (bdry_case.eq.FIXED) then
             do ic$2=ibeg$2,iend$2
                ict$2 = ic$2
                do ic$3=ibeg$3,iend$3
                   ict$3 = ic$3
                   do ic$1=$4
                      ict$1 = $5
fixed_state(`ict0,ict1,ict2')dnl
                   enddo
                enddo
             enddo
           else if (bdry_case.eq.SYMMETRIC) then
             do ic$2=ibeg$2,iend$2
                ict$2 = ic$2
                do ic$3=ibeg$3,iend$3
                   ict$3 = ic$3
                   do ic$1=$4
symmetry(`$1',`$5',`$6')dnl
                   enddo
                enddo
             enddo
          endif
')dnl
c
define(do_bdry_node,`dnl
c          write(6,*) "ibdebug ",ibdebug
c          write(6,*) "iedebug ",iedebug
           if ((bdry_case.eq.FLOW)) then
             do ic0=$1
                ict0 = $4
                do ic1=$2
                   ict1 = ic1
                   do ic2=$3
                      ict2 = ic2
cell_even(`ict0,ict1,ict2')dnl
                   enddo
                enddo
             enddo
           else if (bdry_case.eq.LINEAR) then
             do ic0=$1
                ict0 = $4
                do ic1=$2
                   ict1 = ic1
                   do ic2=$3
                      ict2 = ic2
lin_boundary(`0',`$4',`$5')dnl
                   enddo
                enddo
             enddo
           else if (bdry_case.eq.FIXED) then
             do ic0=$1
                ict0 = $4
                do ic1=$2
                   ict1 = ic1
                   do ic2=$3
                      ict2 = ic2
fixed_state(`ict0,ict1,ict2')dnl
                   enddo
                enddo
             enddo
           else if (bdry_case.eq.SYMMETRIC) then
             do ic0=$1
                ict0 = $4
                do ic1=$2
                   ict1 = ic1
                   do ic2=$3
                      ict2 = ic2
symmetry(`0',`$4',`$5')dnl
                   enddo
                enddo
             enddo
           endif
')dnl
c
define(do_bdry_edge,`dnl
c          write(6,*) "ibdebug ",ibdebug
c          write(6,*) "iedebug ",iedebug
          if ((bdry_case.eq.FLOW)) then
            do ic$2=$4
              ict$2 = ic$2
              do ic$3=$5
                ict$3 = ic$3
                do ic$1=ibeg$1,iend$1
                   ict$1 = ic$1
                   ict$2 = $6
cell_even(`ict0,ict1,ict2')dnl
                enddo
              enddo
            enddo
          else if (bdry_case.eq.LINEAR) then
            do ic$2=$4
              ict$2 = ic$2
              do ic$3=$5
                ict$3 = ic$3
                do ic$1=ibeg$1,iend$1
                   ict$1 = ic$1
lin_boundary($2,`$6',`$7')dnl
                enddo
              enddo
            enddo
          else if (bdry_case.eq.FIXED) then
            do ic$2=$4
              ict$2 = ic$2
              do ic$3=$5
                ict$3 = ic$3
                do ic$1=ibeg$1,iend$1
                   ict$1 = ic$1
                   ict$2 = $6
fixed_state(`ict0,ict1,ict2')dnl
                enddo
              enddo
            enddo
          else if (bdry_case.eq.SYMMETRIC) then
            do ic$2=$4
              ict$2 = ic$2
              do ic$3=$5
                ict$3 = ic$3
                do ic$1=ibeg$1,iend$1
                   ict$1 = ic$1
symmetry($2,`$6',`$7')dnl
                enddo
              enddo
            enddo
          endif
')dnl
c
