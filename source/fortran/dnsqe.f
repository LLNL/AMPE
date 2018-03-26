c Copyright (c) 2018, Lawrence Livermore National Security, LLC.
c Produced at the Lawrence Livermore National Laboratory
c Written by M.R. Dorr, J.-L. Fattebert and M.E. Wickett
c LLNL-CODE-747500
c All rights reserved.
c This file is part of AMPE. 
c For details, see https://github.com/LLNL/AMPE
c Please also read AMPE/LICENSE.
c Redistribution and use in source and binary forms, with or without 
c modification, are permitted provided that the following conditions are met:
c - Redistributions of source code must retain the above copyright notice,
c   this list of conditions and the disclaimer below.
c - Redistributions in binary form must reproduce the above copyright notice,
c   this list of conditions and the disclaimer (as noted below) in the
c   documentation and/or other materials provided with the distribution.
c - Neither the name of the LLNS/LLNL nor the names of its contributors may be
c   used to endorse or promote products derived from this software without
c   specific prior written permission.
c
c THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
c AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
c IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
c ARE DISCLAIMED. IN NO EVENT SHALL LAWRENCE LIVERMORE NATIONAL SECURITY,
c LLC, THE U.S. DEPARTMENT OF ENERGY OR CONTRIBUTORS BE LIABLE FOR ANY
c DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
c DAMAGES  (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS
c OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
c HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
c STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING
c IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
c POSSIBILITY OF SUCH DAMAGE.
c 
CC+------------------------------------------------------------------+CC
CC+------------------------------------------------------------------+CC
CC+                                                                  +CC
CC+              Lawrence Livermore National Laboratory              +CC
CC+                                                                  +CC
CC+        Livermore Computing   Mathematical Software Library       +CC
CC+         SLATEC Common Mathematical Library -- Version 4.1        +CC
CC+                                                                  +CC
CC+------------------------------------------------------------------+CC
CC+                                                                  +CC
CC+  The SLATEC Common Mathematical Library (CML) is issued by the   +CC
CC+  following:                                                      +CC
CC+                                                                  +CC
CC+       Air Force Weapons Laboratory, Albuquerque                  +CC
CC+       Lawrence Livermore National Laboratory, Livermore          +CC
CC+       Los Alamos National Scientific Laboratory, Los Alamos      +CC
CC+       National Institute for Standards and Technology,           +CC
CC+                 Washington, D.C.                                 +CC
CC+       Oak Ridge National Laboratory, Oak Ridge                   +CC
CC+       Sandia National Laboratories, Albuquerque                  +CC
CC+       Sandia National Laboratories, Livermore                    +CC
CC+                                                                  +CC
CC+  SLATEC source codes are distributed by LC                       +CC
CC+  exclusively for use in support of LLNL programs.                +CC
CC+                                                                  +CC
CC+  IMPORTANT:  *** READ THIS BEFORE COMPILING THIS ROUTINE ***     +CC
CC+  ---------                                                       +CC
CC+   1. Note above restrictions on the use of the SLATEC CML.       +CC
CC+   2. Constants in the machine-dependent modules I1MACH, R1MACH,  +CC
CC+      and D1MACH have been set for use on a Sun or IBM-PC.        +CC
CC+      If these are used by this routine, be sure to change to     +CC
CC+      constants appropriate for the target machine.               +CC
CC+                                                                  +CC
CC+  +------------------------------------------------------------+  +CC
CC+  +                        N O T I C E                         +  +CC
CC+  +------------------------------------------------------------+  +CC
CC+  +  This report was prepared as an account of work sponsored  +  +CC
CC+  +  by the United States government.  Neither the United      +  +CC
CC+  +  States government nor any of their employees, nor any of  +  +CC
CC+  +  their contractors, subcontractors, or their employees,    +  +CC
CC+  +  makes any warranty, expressed or implied, or assumes any  +  +CC
CC+  +  legal liability or responsibility for the accuracy,       +  +CC
CC+  +  completeness or usefulness of any information, apparatus, +  +CC
CC+  +  product or process disclosed, or represents that its use  +  +CC
CC+  +  would not infringe privately-owned rights.                +  +CC
CC+  +------------------------------------------------------------+  +CC
CC+                                                                  +CC
CC+  Please report any suspected errors in this routine to the LC    +CC
CC+  Client Services Hotline, (925)422-4531.                         +CC
CC+                                                                  +CC
CC+------------------------------------------------------------------+CC
CC+------------------------------------------------------------------+CC
*DECK DNSQE
      SUBROUTINE DNSQE (FCN, JAC, IOPT, N, X, FVEC, TOL, NPRINT, INFO,
     +   WA, LWA)
C***BEGIN PROLOGUE  DNSQE
C***PURPOSE  An easy-to-use code to find a zero of a system of N
C            nonlinear functions in N variables by a modification of
C            the Powell hybrid method.
C***LIBRARY   SLATEC
C***CATEGORY  F2A
C***TYPE      DOUBLE PRECISION (SNSQE-S, DNSQE-D)
C***KEYWORDS  EASY-TO-USE, NONLINEAR SQUARE SYSTEM,
C             POWELL HYBRID METHOD, ZEROS
C***AUTHOR  Hiebert, K. L. (SNLA)
C***DESCRIPTION
C
C 1. Purpose.
C
C       The purpose of DNSQE is to find a zero of a system of N
C       nonlinear functions in N variables by a modification of the
C       Powell hybrid method.  This is done by using the more general
C       nonlinear equation solver DNSQ.  The user must provide a
C       subroutine which calculates the functions.  The user has the
C       option of either to provide a subroutine which calculates the
C       Jacobian or to let the code calculate it by a forward-difference
C       approximation.  This code is the combination of the MINPACK
C       codes (Argonne) HYBRD1 and HYBRJ1.
C
C 2. Subroutine and Type Statements.
C
C       SUBROUTINE DNSQE(FCN,JAC,IOPT,N,X,FVEC,TOL,NPRINT,INFO,
C      *                  WA,LWA)
C       INTEGER IOPT,N,NPRINT,INFO,LWA
C       DOUBLE PRECISION TOL
C       DOUBLE PRECISION X(N),FVEC(N),WA(LWA)
C       EXTERNAL FCN,JAC
C
C 3. Parameters.
C
C       Parameters designated as input parameters must be specified on
C       entry to DNSQE and are not changed on exit, while parameters
C       designated as output parameters need not be specified on entry
C       and are set to appropriate values on exit from DNSQE.
C
C       FCN is the name of the user-supplied subroutine which calculates
C         the functions.  FCN must be declared in an external statement
C         in the user calling program, and should be written as follows.
C
C         SUBROUTINE FCN(N,X,FVEC,IFLAG)
C         INTEGER N,IFLAG
C         DOUBLE PRECISION X(N),FVEC(N)
C         ----------
C         Calculate the functions at X and
C         return this vector in FVEC.
C         ----------
C         RETURN
C         END
C
C         The value of IFLAG should not be changed by FCN unless the
C         user wants to terminate execution of DNSQE.  In this case set
C         IFLAG to a negative integer.
C
C       JAC is the name of the user-supplied subroutine which calculates
C         the Jacobian.  If IOPT=1, then JAC must be declared in an
C         external statement in the user calling program, and should be
C         written as follows.
C
C         SUBROUTINE JAC(N,X,FVEC,FJAC,LDFJAC,IFLAG)
C         INTEGER N,LDFJAC,IFLAG
C         DOUBLE PRECISION X(N),FVEC(N),FJAC(LDFJAC,N)
C         ----------
C         Calculate the Jacobian at X and return this
C         matrix in FJAC.  FVEC contains the function
C         values at X and should not be altered.
C         ----------
C         RETURN
C         END
C
C         The value of IFLAG should not be changed by JAC unless the
C         user wants to terminate execution of DNSQE. In this case set
C         IFLAG to a negative integer.
C
C         If IOPT=2, JAC can be ignored (treat it as a dummy argument).
C
C       IOPT is an input variable which specifies how the Jacobian will
C         be calculated.  If IOPT=1, then the user must supply the
C         Jacobian through the subroutine JAC.  If IOPT=2, then the
C         code will approximate the Jacobian by forward-differencing.
C
C       N is a positive integer input variable set to the number of
C         functions and variables.
C
C       X is an array of length N.  On input X must contain an initial
C         estimate of the solution vector.  On output X contains the
C         final estimate of the solution vector.
C
C       FVEC is an output array of length N which contains the functions
C         evaluated at the output X.
C
C       TOL is a nonnegative input variable.  Termination occurs when
C         the algorithm estimates that the relative error between X and
C         the solution is at most TOL.  Section 4 contains more details
C         about TOL.
C
C       NPRINT is an integer input variable that enables controlled
C         printing of iterates if it is positive.  In this case, FCN is
C         called with IFLAG = 0 at the beginning of the first iteration
C         and every NPRINT iterations thereafter and immediately prior
C         to return, with X and FVEC available for printing. Appropriate
C         print statements must be added to FCN(see example).  If NPRINT
C         is not positive, no special calls of FCN with IFLAG = 0 are
C         made.
C
C       INFO is an integer output variable.  If the user has terminated
C         execution, INFO is set to the (negative) value of IFLAG.  See
C         description of FCN and JAC. Otherwise, INFO is set as follows.
C
C         INFO = 0  Improper input parameters.
C
C         INFO = 1  Algorithm estimates that the relative error between
C                   X and the solution is at most TOL.
C
C         INFO = 2  Number of calls to FCN has reached or exceeded
C                   100*(N+1) for IOPT=1 or 200*(N+1) for IOPT=2.
C
C         INFO = 3  TOL is too small.  No further improvement in the
C                   approximate solution X is possible.
C
C         INFO = 4  Iteration is not making good progress.
C
C         Sections 4 and 5 contain more details about INFO.
C
C       WA is a work array of length LWA.
C
C       LWA is a positive integer input variable not less than
C         (3*N**2+13*N))/2.
C
C 4. Successful Completion.
C
C       The accuracy of DNSQE is controlled by the convergence parameter
C       TOL.  This parameter is used in a test which makes a comparison
C       between the approximation X and a solution XSOL.  DNSQE
C       terminates when the test is satisfied.  If TOL is less than the
C       machine precision (as defined by the  function D1MACH(4)), then
C       DNSQE only attempts to satisfy the test defined by the machine
C       precision.  Further progress is not usually possible.  Unless
C       high precision solutions are required, the recommended value
C       for TOL is the square root of the machine precision.
C
C       The test assumes that the functions are reasonably well behaved,
C       and, if the Jacobian is supplied by the user, that the functions
C       and the Jacobian are coded consistently. If these conditions are
C       not satisfied, then DNSQE may incorrectly indicate convergence.
C       The coding of the Jacobian can be checked by the subroutine
C       DCKDER.  If the Jacobian is coded correctly or IOPT=2, then
C       the validity of the answer can be checked, for example, by
C       rerunning DNSQE with a tighter tolerance.
C
C       Convergence Test.  If DENORM(Z) denotes the Euclidean norm of a
C         vector Z, then this test attempts to guarantee that
C
C               DENORM(X-XSOL) .LE. TOL*DENORM(XSOL).
C
C         If this condition is satisfied with TOL = 10**(-K), then the
C         larger components of X have K significant decimal digits and
C         INFO is set to 1.  There is a danger that the smaller
C         components of X may have large relative errors, but the fast
C         rate of convergence of DNSQE usually avoids this possibility.
C
C 5. Unsuccessful Completion.
C
C       Unsuccessful termination of DNSQE can be due to improper input
C       parameters, arithmetic interrupts, an excessive number of
C       function evaluations, errors in the functions, or lack of good
C       progress.
C
C       Improper Input Parameters.  INFO is set to 0 if IOPT .LT. 1, or
C         IOPT .GT. 2, or N .LE. 0, or TOL .LT. 0.E0, or
C         LWA .LT. (3*N**2+13*N)/2.
C
C       Arithmetic Interrupts.  If these interrupts occur in the FCN
C         subroutine during an early stage of the computation, they may
C         be caused by an unacceptable choice of X by DNSQE.  In this
C         case, it may be possible to remedy the situation by not
C         evaluating the functions here, but instead setting the
C         components of FVEC to numbers that exceed those in the initial
C         FVEC.
C
C       Excessive Number of Function Evaluations.  If the number of
C         calls to FCN reaches 100*(N+1) for IOPT=1 or 200*(N+1) for
C         IOPT=2, then this indicates that the routine is converging
C         very slowly as measured by the progress of FVEC, and INFO is
C         set to 2.  This situation should be unusual because, as
C         indicated below, lack of good progress is usually diagnosed
C         earlier by DNSQE, causing termination with INFO = 4.
C
C       Errors In the Functions.  When IOPT=2, the choice of step length
C         in the forward-difference approximation to the Jacobian
C         assumes that the relative errors in the functions are of the
C         order of the machine precision.  If this is not the case,
C         DNSQE may fail (usually with INFO = 4).  The user should
C         then either use DNSQ and set the step length or use IOPT=1
C         and supply the Jacobian.
C
C       Lack of Good Progress.  DNSQE searches for a zero of the system
C         by minimizing the sum of the squares of the functions.  In so
C         doing, it can become trapped in a region where the minimum
C         does not correspond to a zero of the system and, in this
C         situation, the iteration eventually fails to make good
C         progress.  In particular, this will happen if the system does
C         not have a zero.  If the system has a zero, rerunning DNSQE
C         from a different starting point may be helpful.
C
C 6. Characteristics of The Algorithm.
C
C       DNSQE is a modification of the Powell Hybrid method.  Two of
C       its main characteristics involve the choice of the correction as
C       a convex combination of the Newton and scaled gradient
C       directions, and the updating of the Jacobian by the rank-1
C       method of Broyden.  The choice of the correction guarantees
C       (under reasonable conditions) global convergence for starting
C       points far from the solution and a fast rate of convergence.
C       The Jacobian is calculated at the starting point by either the
C       user-supplied subroutine or a forward-difference approximation,
C       but it is not recalculated until the rank-1 method fails to
C       produce satisfactory progress.
C
C       Timing.  The time required by DNSQE to solve a given problem
C         depends on N, the behavior of the functions, the accuracy
C         requested, and the starting point.  The number of arithmetic
C         operations needed by DNSQE is about 11.5*(N**2) to process
C         each evaluation of the functions (call to FCN) and 1.3*(N**3)
C         to process each evaluation of the Jacobian (call to JAC,
C         if IOPT = 1).  Unless FCN and JAC can be evaluated quickly,
C         the timing of DNSQE will be strongly influenced by the time
C         spent in FCN and JAC.
C
C       Storage.  DNSQE requires (3*N**2 + 17*N)/2 single precision
C         storage locations, in addition to the storage required by the
C         program.  There are no internally declared storage arrays.
C
C *Long Description:
C
C 7. Example.
C
C       The problem is to determine the values of X(1), X(2), ..., X(9),
C       which solve the system of tridiagonal equations
C
C       (3-2*X(1))*X(1)           -2*X(2)                   = -1
C               -X(I-1) + (3-2*X(I))*X(I)         -2*X(I+1) = -1, I=2-8
C                                   -X(8) + (3-2*X(9))*X(9) = -1
C
C       **********
C
C       PROGRAM TEST
C C
C C     DRIVER FOR DNSQE EXAMPLE.
C C
C       INTEGER J,N,IOPT,NPRINT,INFO,LWA,NWRITE
C       DOUBLE PRECISION TOL,FNORM
C       DOUBLE PRECISION X(9),FVEC(9),WA(180)
C       DOUBLE PRECISION DENORM,D1MACH
C       EXTERNAL FCN
C       DATA NWRITE /6/
C C
C       IOPT = 2
C       N = 9
C C
C C     THE FOLLOWING STARTING VALUES PROVIDE A ROUGH SOLUTION.
C C
C       DO 10 J = 1, 9
C          X(J) = -1.E0
C    10    CONTINUE
C
C       LWA = 180
C       NPRINT = 0
C C
C C     SET TOL TO THE SQUARE ROOT OF THE MACHINE PRECISION.
C C     UNLESS HIGH PRECISION SOLUTIONS ARE REQUIRED,
C C     THIS IS THE RECOMMENDED SETTING.
C C
C       TOL = SQRT(D1MACH(4))
C C
C       CALL DNSQE(FCN,JAC,IOPT,N,X,FVEC,TOL,NPRINT,INFO,WA,LWA)
C       FNORM = DENORM(N,FVEC)
C       WRITE (NWRITE,1000) FNORM,INFO,(X(J),J=1,N)
C       STOP
C  1000 FORMAT (5X,' FINAL L2 NORM OF THE RESIDUALS',E15.7 //
C      *        5X,' EXIT PARAMETER',16X,I10 //
C      *        5X,' FINAL APPROXIMATE SOLUTION' // (5X,3E15.7))
C       END
C       SUBROUTINE FCN(N,X,FVEC,IFLAG)
C       INTEGER N,IFLAG
C       DOUBLE PRECISION X(N),FVEC(N)
C       INTEGER K
C       DOUBLE PRECISION ONE,TEMP,TEMP1,TEMP2,THREE,TWO,ZERO
C       DATA ZERO,ONE,TWO,THREE /0.E0,1.E0,2.E0,3.E0/
C C
C       DO 10 K = 1, N
C          TEMP = (THREE - TWO*X(K))*X(K)
C          TEMP1 = ZERO
C          IF (K .NE. 1) TEMP1 = X(K-1)
C          TEMP2 = ZERO
C          IF (K .NE. N) TEMP2 = X(K+1)
C          FVEC(K) = TEMP - TEMP1 - TWO*TEMP2 + ONE
C    10    CONTINUE
C       RETURN
C       END
C
C       RESULTS OBTAINED WITH DIFFERENT COMPILERS OR MACHINES
C       MAY BE SLIGHTLY DIFFERENT.
C
C       FINAL L2 NORM OF THE RESIDUALS  0.1192636E-07
C
C       EXIT PARAMETER                         1
C
C       FINAL APPROXIMATE SOLUTION
C
C       -0.5706545E+00 -0.6816283E+00 -0.7017325E+00
C       -0.7042129E+00 -0.7013690E+00 -0.6918656E+00
C       -0.6657920E+00 -0.5960342E+00 -0.4164121E+00
C
C***REFERENCES  M. J. D. Powell, A hybrid method for nonlinear equa-
C                 tions. In Numerical Methods for Nonlinear Algebraic
C                 Equations, P. Rabinowitz, Editor.  Gordon and Breach,
C                 1988.
C***ROUTINES CALLED  DNSQ, XERMSG
C***REVISION HISTORY  (YYMMDD)
C   800301  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890831  Modified array declarations.  (WRB)
C   890831  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  DNSQE
      INTEGER INDEX, INFO, IOPT, J, LR, LWA, MAXFEV, ML, MODE, MU, N,
     1     NFEV, NJEV, NPRINT
      DOUBLE PRECISION EPSFCN, FACTOR, FVEC(*), ONE, TOL, WA(*),
     1     X(*), XTOL, ZERO
      EXTERNAL FCN, JAC
      SAVE FACTOR, ONE, ZERO
      DATA FACTOR,ONE,ZERO /1.0D2,1.0D0,0.0D0/
C     BEGIN BLOCK PERMITTING ...EXITS TO 20
C***FIRST EXECUTABLE STATEMENT  DNSQE
         INFO = 0
C
C        CHECK THE INPUT PARAMETERS FOR ERRORS.
C
C     ...EXIT
         IF (IOPT .LT. 1 .OR. IOPT .GT. 2 .OR. N .LE. 0
     1       .OR. TOL .LT. ZERO .OR. LWA .LT. (3*N**2 + 13*N)/2)
     2      GO TO 20
C
C        CALL DNSQ.
C
         MAXFEV = 100*(N + 1)
         IF (IOPT .EQ. 2) MAXFEV = 2*MAXFEV
         XTOL = TOL
         ML = N - 1
         MU = N - 1
         EPSFCN = ZERO
         MODE = 2
         DO 10 J = 1, N
            WA(J) = ONE
   10    CONTINUE
         LR = (N*(N + 1))/2
         INDEX = 6*N + LR
         CALL DNSQ(FCN,JAC,IOPT,N,X,FVEC,WA(INDEX+1),N,XTOL,MAXFEV,ML,
     1             MU,EPSFCN,WA(1),MODE,FACTOR,NPRINT,INFO,NFEV,NJEV,
     2             WA(6*N+1),LR,WA(N+1),WA(2*N+1),WA(3*N+1),WA(4*N+1),
     3             WA(5*N+1))
         IF (INFO .EQ. 5) INFO = 4
   20 CONTINUE
      IF (INFO .EQ. 0) CALL XERMSG ('SLATEC', 'DNSQE',
     +   'INVALID INPUT PARAMETER.', 2, 1)
      RETURN
C
C     LAST CARD OF SUBROUTINE DNSQE.
C
      END
*DECK D1MACH
      DOUBLE PRECISION FUNCTION D1MACH (I)
C***BEGIN PROLOGUE  D1MACH
C***SUBSIDIARY
C***PURPOSE  Return floating point machine dependent constants.
C***LIBRARY   SLATEC
C***CATEGORY  R1
C***TYPE      DOUBLE PRECISION (R1MACH-S, D1MACH-D)
C***KEYWORDS  MACHINE CONSTANTS
C***AUTHOR  Fox, P. A., (Bell Labs)
C           Hall, A. D., (Bell Labs)
C           Schryer, N. L., (Bell Labs)
C***DESCRIPTION
C
C   D1MACH can be used to obtain machine-dependent parameters for the
C   local machine environment.  It is a function subprogram with one
C   (input) argument, and can be referenced as follows:
C
C        D = D1MACH(I)
C
C   where I=1,...,5.  The (output) value of D above is determined by
C   the (input) value of I.  The results for various values of I are
C   discussed below.
C
C   D1MACH( 1) = B**(EMIN-1), the smallest positive magnitude.
C   D1MACH( 2) = B**EMAX*(1 - B**(-T)), the largest magnitude.
C   D1MACH( 3) = B**(-T), the smallest relative spacing.
C   D1MACH( 4) = B**(1-T), the largest relative spacing.
C   D1MACH( 5) = LOG10(B)
C
C   Assume double precision numbers are represented in the T-digit,
C   base-B form
C
C              sign (B**E)*( (X(1)/B) + ... + (X(T)/B**T) )
C
C   where 0 .LE. X(I) .LT. B for I=1,...,T, 0 .LT. X(1), and
C   EMIN .LE. E .LE. EMAX.
C
C   The values of B, T, EMIN and EMAX are provided in I1MACH as
C   follows:
C   I1MACH(10) = B, the base.
C   I1MACH(14) = T, the number of base-B digits.
C   I1MACH(15) = EMIN, the smallest exponent E.
C   I1MACH(16) = EMAX, the largest exponent E.
C
C   To alter this function for a particular environment, the desired
C   set of DATA statements should be activated by removing the C from
C   column 1.  Also, the values of D1MACH(1) - D1MACH(4) should be
C   checked for consistency with the local operating system.
C
C***REFERENCES  P. A. Fox, A. D. Hall and N. L. Schryer, Framework for
C                 a portable library, ACM Transactions on Mathematical
C                 Software 4, 2 (June 1978), pp. 177-188.
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   750101  DATE WRITTEN
C   890213  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900315  Calls to XERROR changed to calls to XERMSG.  (THJ)
C   900618  Added DEC RISC constants.  (WRB)
C   900723  Added IBM RS 6000 constants.  (WRB)
C   900911  Added SUN 386i constants.  (WRB)
C   910710  Added HP 730 constants.  (SMR)
C   911114  Added Convex IEEE constants.  (WRB)
C   920121  Added SUN -r8 compiler option constants.  (WRB)
C   920229  Added Touchstone Delta i860 constants.  (WRB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C   920625  Added CONVEX -p8/-pd8 compiler option constants. (BKS, WRB)
C   930201  Added DEC Alpha and SGI constants.  (RWC and WRB)
C   011206  Activated IBM-PC block, to give decimal values. (ACH)
C   020227  Corrected DATA statements in IBM-PC block. (ACH)
C   020313  Removed the argument checking and XERMSG call. (ACH and TS)
C***END PROLOGUE  D1MACH
C
      INTEGER SMALL(4)
      INTEGER LARGE(4)
      INTEGER RIGHT(4)
      INTEGER DIVER(4)
      INTEGER LOG10(4)
C
      DOUBLE PRECISION DMACH(5)
      SAVE DMACH
C
      EQUIVALENCE (DMACH(1),SMALL(1))
      EQUIVALENCE (DMACH(2),LARGE(1))
      EQUIVALENCE (DMACH(3),RIGHT(1))
      EQUIVALENCE (DMACH(4),DIVER(1))
      EQUIVALENCE (DMACH(5),LOG10(1))
C
C     MACHINE CONSTANTS FOR THE AMIGA
C     ABSOFT FORTRAN COMPILER USING THE 68020/68881 COMPILER OPTION
C
C     DATA SMALL(1), SMALL(2) / Z'00100000', Z'00000000' /
C     DATA LARGE(1), LARGE(2) / Z'7FEFFFFF', Z'FFFFFFFF' /
C     DATA RIGHT(1), RIGHT(2) / Z'3CA00000', Z'00000000' /
C     DATA DIVER(1), DIVER(2) / Z'3CB00000', Z'00000000' /
C     DATA LOG10(1), LOG10(2) / Z'3FD34413', Z'509F79FF' /
C
C     MACHINE CONSTANTS FOR THE AMIGA
C     ABSOFT FORTRAN COMPILER USING SOFTWARE FLOATING POINT
C
C     DATA SMALL(1), SMALL(2) / Z'00100000', Z'00000000' /
C     DATA LARGE(1), LARGE(2) / Z'7FDFFFFF', Z'FFFFFFFF' /
C     DATA RIGHT(1), RIGHT(2) / Z'3CA00000', Z'00000000' /
C     DATA DIVER(1), DIVER(2) / Z'3CB00000', Z'00000000' /
C     DATA LOG10(1), LOG10(2) / Z'3FD34413', Z'509F79FF' /
C
C     MACHINE CONSTANTS FOR THE APOLLO
C
C     DATA SMALL(1), SMALL(2) / 16#00100000, 16#00000000 /
C     DATA LARGE(1), LARGE(2) / 16#7FFFFFFF, 16#FFFFFFFF /
C     DATA RIGHT(1), RIGHT(2) / 16#3CA00000, 16#00000000 /
C     DATA DIVER(1), DIVER(2) / 16#3CB00000, 16#00000000 /
C     DATA LOG10(1), LOG10(2) / 16#3FD34413, 16#509F79FF /
C
C     MACHINE CONSTANTS FOR THE BURROUGHS 1700 SYSTEM
C
C     DATA SMALL(1) / ZC00800000 /
C     DATA SMALL(2) / Z000000000 /
C     DATA LARGE(1) / ZDFFFFFFFF /
C     DATA LARGE(2) / ZFFFFFFFFF /
C     DATA RIGHT(1) / ZCC5800000 /
C     DATA RIGHT(2) / Z000000000 /
C     DATA DIVER(1) / ZCC6800000 /
C     DATA DIVER(2) / Z000000000 /
C     DATA LOG10(1) / ZD00E730E7 /
C     DATA LOG10(2) / ZC77800DC0 /
C
C     MACHINE CONSTANTS FOR THE BURROUGHS 5700 SYSTEM
C
C     DATA SMALL(1) / O1771000000000000 /
C     DATA SMALL(2) / O0000000000000000 /
C     DATA LARGE(1) / O0777777777777777 /
C     DATA LARGE(2) / O0007777777777777 /
C     DATA RIGHT(1) / O1461000000000000 /
C     DATA RIGHT(2) / O0000000000000000 /
C     DATA DIVER(1) / O1451000000000000 /
C     DATA DIVER(2) / O0000000000000000 /
C     DATA LOG10(1) / O1157163034761674 /
C     DATA LOG10(2) / O0006677466732724 /
C
C     MACHINE CONSTANTS FOR THE BURROUGHS 6700/7700 SYSTEMS
C
C     DATA SMALL(1) / O1771000000000000 /
C     DATA SMALL(2) / O7770000000000000 /
C     DATA LARGE(1) / O0777777777777777 /
C     DATA LARGE(2) / O7777777777777777 /
C     DATA RIGHT(1) / O1461000000000000 /
C     DATA RIGHT(2) / O0000000000000000 /
C     DATA DIVER(1) / O1451000000000000 /
C     DATA DIVER(2) / O0000000000000000 /
C     DATA LOG10(1) / O1157163034761674 /
C     DATA LOG10(2) / O0006677466732724 /
C
C     MACHINE CONSTANTS FOR THE CDC 170/180 SERIES USING NOS/VE
C
C     DATA SMALL(1) / Z"3001800000000000" /
C     DATA SMALL(2) / Z"3001000000000000" /
C     DATA LARGE(1) / Z"4FFEFFFFFFFFFFFE" /
C     DATA LARGE(2) / Z"4FFE000000000000" /
C     DATA RIGHT(1) / Z"3FD2800000000000" /
C     DATA RIGHT(2) / Z"3FD2000000000000" /
C     DATA DIVER(1) / Z"3FD3800000000000" /
C     DATA DIVER(2) / Z"3FD3000000000000" /
C     DATA LOG10(1) / Z"3FFF9A209A84FBCF" /
C     DATA LOG10(2) / Z"3FFFF7988F8959AC" /
C
C     MACHINE CONSTANTS FOR THE CDC 6000/7000 SERIES
C
C     DATA SMALL(1) / 00564000000000000000B /
C     DATA SMALL(2) / 00000000000000000000B /
C     DATA LARGE(1) / 37757777777777777777B /
C     DATA LARGE(2) / 37157777777777777777B /
C     DATA RIGHT(1) / 15624000000000000000B /
C     DATA RIGHT(2) / 00000000000000000000B /
C     DATA DIVER(1) / 15634000000000000000B /
C     DATA DIVER(2) / 00000000000000000000B /
C     DATA LOG10(1) / 17164642023241175717B /
C     DATA LOG10(2) / 16367571421742254654B /
C
C     MACHINE CONSTANTS FOR THE CELERITY C1260
C
C     DATA SMALL(1), SMALL(2) / Z'00100000', Z'00000000' /
C     DATA LARGE(1), LARGE(2) / Z'7FEFFFFF', Z'FFFFFFFF' /
C     DATA RIGHT(1), RIGHT(2) / Z'3CA00000', Z'00000000' /
C     DATA DIVER(1), DIVER(2) / Z'3CB00000', Z'00000000' /
C     DATA LOG10(1), LOG10(2) / Z'3FD34413', Z'509F79FF' /
C
C     MACHINE CONSTANTS FOR THE CONVEX
C     USING THE -fn OR -pd8 COMPILER OPTION
C
C     DATA DMACH(1) / Z'0010000000000000' /
C     DATA DMACH(2) / Z'7FFFFFFFFFFFFFFF' /
C     DATA DMACH(3) / Z'3CC0000000000000' /
C     DATA DMACH(4) / Z'3CD0000000000000' /
C     DATA DMACH(5) / Z'3FF34413509F79FF' /
C
C     MACHINE CONSTANTS FOR THE CONVEX
C     USING THE -fi COMPILER OPTION
C
C     DATA DMACH(1) / Z'0010000000000000' /
C     DATA DMACH(2) / Z'7FEFFFFFFFFFFFFF' /
C     DATA DMACH(3) / Z'3CA0000000000000' /
C     DATA DMACH(4) / Z'3CB0000000000000' /
C     DATA DMACH(5) / Z'3FD34413509F79FF' /
C
C     MACHINE CONSTANTS FOR THE CONVEX
C     USING THE -p8 COMPILER OPTION
C
C     DATA DMACH(1) / Z'00010000000000000000000000000000' /
C     DATA DMACH(2) / Z'7FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF' /
C     DATA DMACH(3) / Z'3F900000000000000000000000000000' /
C     DATA DMACH(4) / Z'3F910000000000000000000000000000' /
C     DATA DMACH(5) / Z'3FFF34413509F79FEF311F12B35816F9' /
C
C     MACHINE CONSTANTS FOR THE CRAY
C
C     DATA SMALL(1) / 201354000000000000000B /
C     DATA SMALL(2) / 000000000000000000000B /
C     DATA LARGE(1) / 577767777777777777777B /
C     DATA LARGE(2) / 000007777777777777774B /
C     DATA RIGHT(1) / 376434000000000000000B /
C     DATA RIGHT(2) / 000000000000000000000B /
C     DATA DIVER(1) / 376444000000000000000B /
C     DATA DIVER(2) / 000000000000000000000B /
C     DATA LOG10(1) / 377774642023241175717B /
C     DATA LOG10(2) / 000007571421742254654B /
C
C     MACHINE CONSTANTS FOR THE DATA GENERAL ECLIPSE S/200
C     NOTE - IT MAY BE APPROPRIATE TO INCLUDE THE FOLLOWING CARD -
C     STATIC DMACH(5)
C
C     DATA SMALL /    20K, 3*0 /
C     DATA LARGE / 77777K, 3*177777K /
C     DATA RIGHT / 31420K, 3*0 /
C     DATA DIVER / 32020K, 3*0 /
C     DATA LOG10 / 40423K, 42023K, 50237K, 74776K /
C
C     MACHINE CONSTANTS FOR THE DEC ALPHA
C     USING G_FLOAT
C
C     DATA DMACH(1) / '0000000000000010'X /
C     DATA DMACH(2) / 'FFFFFFFFFFFF7FFF'X /
C     DATA DMACH(3) / '0000000000003CC0'X /
C     DATA DMACH(4) / '0000000000003CD0'X /
C     DATA DMACH(5) / '79FF509F44133FF3'X /
C
C     MACHINE CONSTANTS FOR THE DEC ALPHA
C     USING IEEE_FORMAT
C
C     DATA DMACH(1) / '0010000000000000'X /
C     DATA DMACH(2) / '7FEFFFFFFFFFFFFF'X /
C     DATA DMACH(3) / '3CA0000000000000'X /
C     DATA DMACH(4) / '3CB0000000000000'X /
C     DATA DMACH(5) / '3FD34413509F79FF'X /
C
C     MACHINE CONSTANTS FOR THE DEC RISC
C
C     DATA SMALL(1), SMALL(2) / Z'00000000', Z'00100000'/
C     DATA LARGE(1), LARGE(2) / Z'FFFFFFFF', Z'7FEFFFFF'/
C     DATA RIGHT(1), RIGHT(2) / Z'00000000', Z'3CA00000'/
C     DATA DIVER(1), DIVER(2) / Z'00000000', Z'3CB00000'/
C     DATA LOG10(1), LOG10(2) / Z'509F79FF', Z'3FD34413'/
C
C     MACHINE CONSTANTS FOR THE DEC VAX
C     USING D_FLOATING
C     (EXPRESSED IN INTEGER AND HEXADECIMAL)
C     THE HEX FORMAT BELOW MAY NOT BE SUITABLE FOR UNIX SYSTEMS
C     THE INTEGER FORMAT SHOULD BE OK FOR UNIX SYSTEMS
C
C     DATA SMALL(1), SMALL(2) /        128,           0 /
C     DATA LARGE(1), LARGE(2) /     -32769,          -1 /
C     DATA RIGHT(1), RIGHT(2) /       9344,           0 /
C     DATA DIVER(1), DIVER(2) /       9472,           0 /
C     DATA LOG10(1), LOG10(2) /  546979738,  -805796613 /
C
C     DATA SMALL(1), SMALL(2) / Z00000080, Z00000000 /
C     DATA LARGE(1), LARGE(2) / ZFFFF7FFF, ZFFFFFFFF /
C     DATA RIGHT(1), RIGHT(2) / Z00002480, Z00000000 /
C     DATA DIVER(1), DIVER(2) / Z00002500, Z00000000 /
C     DATA LOG10(1), LOG10(2) / Z209A3F9A, ZCFF884FB /
C
C     MACHINE CONSTANTS FOR THE DEC VAX
C     USING G_FLOATING
C     (EXPRESSED IN INTEGER AND HEXADECIMAL)
C     THE HEX FORMAT BELOW MAY NOT BE SUITABLE FOR UNIX SYSTEMS
C     THE INTEGER FORMAT SHOULD BE OK FOR UNIX SYSTEMS
C
C     DATA SMALL(1), SMALL(2) /         16,           0 /
C     DATA LARGE(1), LARGE(2) /     -32769,          -1 /
C     DATA RIGHT(1), RIGHT(2) /      15552,           0 /
C     DATA DIVER(1), DIVER(2) /      15568,           0 /
C     DATA LOG10(1), LOG10(2) /  1142112243, 2046775455 /
C
C     DATA SMALL(1), SMALL(2) / Z00000010, Z00000000 /
C     DATA LARGE(1), LARGE(2) / ZFFFF7FFF, ZFFFFFFFF /
C     DATA RIGHT(1), RIGHT(2) / Z00003CC0, Z00000000 /
C     DATA DIVER(1), DIVER(2) / Z00003CD0, Z00000000 /
C     DATA LOG10(1), LOG10(2) / Z44133FF3, Z79FF509F /
C
C     MACHINE CONSTANTS FOR THE ELXSI 6400
C     (ASSUMING REAL*8 IS THE DEFAULT DOUBLE PRECISION)
C
C     DATA SMALL(1), SMALL(2) / '00100000'X,'00000000'X /
C     DATA LARGE(1), LARGE(2) / '7FEFFFFF'X,'FFFFFFFF'X /
C     DATA RIGHT(1), RIGHT(2) / '3CB00000'X,'00000000'X /
C     DATA DIVER(1), DIVER(2) / '3CC00000'X,'00000000'X /
C     DATA LOG10(1), LOG10(2) / '3FD34413'X,'509F79FF'X /
C
C     MACHINE CONSTANTS FOR THE HARRIS 220
C
C     DATA SMALL(1), SMALL(2) / '20000000, '00000201 /
C     DATA LARGE(1), LARGE(2) / '37777777, '37777577 /
C     DATA RIGHT(1), RIGHT(2) / '20000000, '00000333 /
C     DATA DIVER(1), DIVER(2) / '20000000, '00000334 /
C     DATA LOG10(1), LOG10(2) / '23210115, '10237777 /
C
C     MACHINE CONSTANTS FOR THE HONEYWELL 600/6000 SERIES
C
C     DATA SMALL(1), SMALL(2) / O402400000000, O000000000000 /
C     DATA LARGE(1), LARGE(2) / O376777777777, O777777777777 /
C     DATA RIGHT(1), RIGHT(2) / O604400000000, O000000000000 /
C     DATA DIVER(1), DIVER(2) / O606400000000, O000000000000 /
C     DATA LOG10(1), LOG10(2) / O776464202324, O117571775714 /
C
C     MACHINE CONSTANTS FOR THE HP 730
C
C     DATA DMACH(1) / Z'0010000000000000' /
C     DATA DMACH(2) / Z'7FEFFFFFFFFFFFFF' /
C     DATA DMACH(3) / Z'3CA0000000000000' /
C     DATA DMACH(4) / Z'3CB0000000000000' /
C     DATA DMACH(5) / Z'3FD34413509F79FF' /
C
C     MACHINE CONSTANTS FOR THE HP 2100
C     THREE WORD DOUBLE PRECISION OPTION WITH FTN4
C
C     DATA SMALL(1), SMALL(2), SMALL(3) / 40000B,       0,       1 /
C     DATA LARGE(1), LARGE(2), LARGE(3) / 77777B, 177777B, 177776B /
C     DATA RIGHT(1), RIGHT(2), RIGHT(3) / 40000B,       0,    265B /
C     DATA DIVER(1), DIVER(2), DIVER(3) / 40000B,       0,    276B /
C     DATA LOG10(1), LOG10(2), LOG10(3) / 46420B,  46502B,  77777B /
C
C     MACHINE CONSTANTS FOR THE HP 2100
C     FOUR WORD DOUBLE PRECISION OPTION WITH FTN4
C
C     DATA SMALL(1), SMALL(2) /  40000B,       0 /
C     DATA SMALL(3), SMALL(4) /       0,       1 /
C     DATA LARGE(1), LARGE(2) /  77777B, 177777B /
C     DATA LARGE(3), LARGE(4) / 177777B, 177776B /
C     DATA RIGHT(1), RIGHT(2) /  40000B,       0 /
C     DATA RIGHT(3), RIGHT(4) /       0,    225B /
C     DATA DIVER(1), DIVER(2) /  40000B,       0 /
C     DATA DIVER(3), DIVER(4) /       0,    227B /
C     DATA LOG10(1), LOG10(2) /  46420B,  46502B /
C     DATA LOG10(3), LOG10(4) /  76747B, 176377B /
C
C     MACHINE CONSTANTS FOR THE HP 9000
C
C     DATA SMALL(1), SMALL(2) / 00040000000B, 00000000000B /
C     DATA LARGE(1), LARGE(2) / 17737777777B, 37777777777B /
C     DATA RIGHT(1), RIGHT(2) / 07454000000B, 00000000000B /
C     DATA DIVER(1), DIVER(2) / 07460000000B, 00000000000B /
C     DATA LOG10(1), LOG10(2) / 07764642023B, 12047674777B /
C
C     MACHINE CONSTANTS FOR THE IBM 360/370 SERIES,
C     THE XEROX SIGMA 5/7/9, THE SEL SYSTEMS 85/86, AND
C     THE PERKIN ELMER (INTERDATA) 7/32.
C
C     DATA SMALL(1), SMALL(2) / Z00100000, Z00000000 /
C     DATA LARGE(1), LARGE(2) / Z7FFFFFFF, ZFFFFFFFF /
C     DATA RIGHT(1), RIGHT(2) / Z33100000, Z00000000 /
C     DATA DIVER(1), DIVER(2) / Z34100000, Z00000000 /
C     DATA LOG10(1), LOG10(2) / Z41134413, Z509F79FF /
C
C     MACHINE CONSTANTS FOR THE IBM PC
C     ASSUMES THAT ALL ARITHMETIC IS DONE IN DOUBLE PRECISION
C     ON 8088, I.E., NOT IN 80 BIT FORM FOR THE 8087.
C
      DATA DMACH(1) / 2.23D-308  /
      DATA DMACH(2) / 1.79D+308  /
      DATA DMACH(3) / 1.11D-16   /
      DATA DMACH(4) / 2.22D-16   /
      DATA DMACH(5) / 0.301029995663981195D0 /
C
C     MACHINE CONSTANTS FOR THE IBM RS 6000
C
C     DATA DMACH(1) / Z'0010000000000000' /
C     DATA DMACH(2) / Z'7FEFFFFFFFFFFFFF' /
C     DATA DMACH(3) / Z'3CA0000000000000' /
C     DATA DMACH(4) / Z'3CB0000000000000' /
C     DATA DMACH(5) / Z'3FD34413509F79FF' /
C
C     MACHINE CONSTANTS FOR THE INTEL i860
C
C     DATA DMACH(1) / Z'0010000000000000' /
C     DATA DMACH(2) / Z'7FEFFFFFFFFFFFFF' /
C     DATA DMACH(3) / Z'3CA0000000000000' /
C     DATA DMACH(4) / Z'3CB0000000000000' /
C     DATA DMACH(5) / Z'3FD34413509F79FF' /
C
C     MACHINE CONSTANTS FOR THE PDP-10 (KA PROCESSOR)
C
C     DATA SMALL(1), SMALL(2) / "033400000000, "000000000000 /
C     DATA LARGE(1), LARGE(2) / "377777777777, "344777777777 /
C     DATA RIGHT(1), RIGHT(2) / "113400000000, "000000000000 /
C     DATA DIVER(1), DIVER(2) / "114400000000, "000000000000 /
C     DATA LOG10(1), LOG10(2) / "177464202324, "144117571776 /
C
C     MACHINE CONSTANTS FOR THE PDP-10 (KI PROCESSOR)
C
C     DATA SMALL(1), SMALL(2) / "000400000000, "000000000000 /
C     DATA LARGE(1), LARGE(2) / "377777777777, "377777777777 /
C     DATA RIGHT(1), RIGHT(2) / "103400000000, "000000000000 /
C     DATA DIVER(1), DIVER(2) / "104400000000, "000000000000 /
C     DATA LOG10(1), LOG10(2) / "177464202324, "476747767461 /
C
C     MACHINE CONSTANTS FOR PDP-11 FORTRAN SUPPORTING
C     32-BIT INTEGERS (EXPRESSED IN INTEGER AND OCTAL).
C
C     DATA SMALL(1), SMALL(2) /    8388608,           0 /
C     DATA LARGE(1), LARGE(2) / 2147483647,          -1 /
C     DATA RIGHT(1), RIGHT(2) /  612368384,           0 /
C     DATA DIVER(1), DIVER(2) /  620756992,           0 /
C     DATA LOG10(1), LOG10(2) / 1067065498, -2063872008 /
C
C     DATA SMALL(1), SMALL(2) / O00040000000, O00000000000 /
C     DATA LARGE(1), LARGE(2) / O17777777777, O37777777777 /
C     DATA RIGHT(1), RIGHT(2) / O04440000000, O00000000000 /
C     DATA DIVER(1), DIVER(2) / O04500000000, O00000000000 /
C     DATA LOG10(1), LOG10(2) / O07746420232, O20476747770 /
C
C     MACHINE CONSTANTS FOR PDP-11 FORTRAN SUPPORTING
C     16-BIT INTEGERS (EXPRESSED IN INTEGER AND OCTAL).
C
C     DATA SMALL(1), SMALL(2) /    128,      0 /
C     DATA SMALL(3), SMALL(4) /      0,      0 /
C     DATA LARGE(1), LARGE(2) /  32767,     -1 /
C     DATA LARGE(3), LARGE(4) /     -1,     -1 /
C     DATA RIGHT(1), RIGHT(2) /   9344,      0 /
C     DATA RIGHT(3), RIGHT(4) /      0,      0 /
C     DATA DIVER(1), DIVER(2) /   9472,      0 /
C     DATA DIVER(3), DIVER(4) /      0,      0 /
C     DATA LOG10(1), LOG10(2) /  16282,   8346 /
C     DATA LOG10(3), LOG10(4) / -31493, -12296 /
C
C     DATA SMALL(1), SMALL(2) / O000200, O000000 /
C     DATA SMALL(3), SMALL(4) / O000000, O000000 /
C     DATA LARGE(1), LARGE(2) / O077777, O177777 /
C     DATA LARGE(3), LARGE(4) / O177777, O177777 /
C     DATA RIGHT(1), RIGHT(2) / O022200, O000000 /
C     DATA RIGHT(3), RIGHT(4) / O000000, O000000 /
C     DATA DIVER(1), DIVER(2) / O022400, O000000 /
C     DATA DIVER(3), DIVER(4) / O000000, O000000 /
C     DATA LOG10(1), LOG10(2) / O037632, O020232 /
C     DATA LOG10(3), LOG10(4) / O102373, O147770 /
C
C     MACHINE CONSTANTS FOR THE SILICON GRAPHICS
C
C     DATA SMALL(1), SMALL(2) / Z'00100000', Z'00000000' /
C     DATA LARGE(1), LARGE(2) / Z'7FEFFFFF', Z'FFFFFFFF' /
C     DATA RIGHT(1), RIGHT(2) / Z'3CA00000', Z'00000000' /
C     DATA DIVER(1), DIVER(2) / Z'3CB00000', Z'00000000' /
C     DATA LOG10(1), LOG10(2) / Z'3FD34413', Z'509F79FF' /
C
C     MACHINE CONSTANTS FOR THE SUN
C
C     DATA DMACH(1) / Z'0010000000000000' /
C     DATA DMACH(2) / Z'7FEFFFFFFFFFFFFF' /
C     DATA DMACH(3) / Z'3CA0000000000000' /
C     DATA DMACH(4) / Z'3CB0000000000000' /
C     DATA DMACH(5) / Z'3FD34413509F79FF' /
C
C     MACHINE CONSTANTS FOR THE SUN
C     USING THE -r8 COMPILER OPTION
C
C     DATA DMACH(1) / Z'00010000000000000000000000000000' /
C     DATA DMACH(2) / Z'7FFEFFFFFFFFFFFFFFFFFFFFFFFFFFFF' /
C     DATA DMACH(3) / Z'3F8E0000000000000000000000000000' /
C     DATA DMACH(4) / Z'3F8F0000000000000000000000000000' /
C     DATA DMACH(5) / Z'3FFD34413509F79FEF311F12B35816F9' /
C
C     MACHINE CONSTANTS FOR THE SUN 386i
C
C     DATA SMALL(1), SMALL(2) / Z'FFFFFFFD', Z'000FFFFF' /
C     DATA LARGE(1), LARGE(2) / Z'FFFFFFB0', Z'7FEFFFFF' /
C     DATA RIGHT(1), RIGHT(2) / Z'000000B0', Z'3CA00000' /
C     DATA DIVER(1), DIVER(2) / Z'FFFFFFCB', Z'3CAFFFFF'
C     DATA LOG10(1), LOG10(2) / Z'509F79E9', Z'3FD34413' /
C
C     MACHINE CONSTANTS FOR THE UNIVAC 1100 SERIES FTN COMPILER
C
C     DATA SMALL(1), SMALL(2) / O000040000000, O000000000000 /
C     DATA LARGE(1), LARGE(2) / O377777777777, O777777777777 /
C     DATA RIGHT(1), RIGHT(2) / O170540000000, O000000000000 /
C     DATA DIVER(1), DIVER(2) / O170640000000, O000000000000 /
C     DATA LOG10(1), LOG10(2) / O177746420232, O411757177572 /
C
C***FIRST EXECUTABLE STATEMENT  D1MACH
      D1MACH = DMACH(I)
      RETURN
C
      END
*DECK D1MPYQ
      SUBROUTINE D1MPYQ (M, N, A, LDA, V, W)
C***BEGIN PROLOGUE  D1MPYQ
C***SUBSIDIARY
C***PURPOSE  Subsidiary to DNSQ and DNSQE
C***LIBRARY   SLATEC
C***TYPE      DOUBLE PRECISION (R1MPYQ-S, D1MPYQ-D)
C***AUTHOR  (UNKNOWN)
C***DESCRIPTION
C
C     Given an M by N matrix A, this subroutine computes A*Q where
C     Q is the product of 2*(N - 1) transformations
C
C           GV(N-1)*...*GV(1)*GW(1)*...*GW(N-1)
C
C     and GV(I), GW(I) are Givens rotations in the (I,N) plane which
C     eliminate elements in the I-th and N-th planes, respectively.
C     Q itself is not given, rather the information to recover the
C     GV, GW rotations is supplied.
C
C     The SUBROUTINE statement is
C
C       SUBROUTINE D1MPYQ(M,N,A,LDA,V,W)
C
C     where
C
C       M is a positive integer input variable set to the number
C         of rows of A.
C
C       N IS a positive integer input variable set to the number
C         of columns of A.
C
C       A is an M by N array. On input A must contain the matrix
C         to be postmultiplied by the orthogonal matrix Q
C         described above. On output A*Q has replaced A.
C
C       LDA is a positive integer input variable not less than M
C         which specifies the leading dimension of the array A.
C
C       V is an input array of length N. V(I) must contain the
C         information necessary to recover the Givens rotation GV(I)
C         described above.
C
C       W is an input array of length N. W(I) must contain the
C         information necessary to recover the Givens rotation GW(I)
C         described above.
C
C***SEE ALSO  DNSQ, DNSQE
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   800301  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890831  Modified array declarations.  (WRB)
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900326  Removed duplicate information from DESCRIPTION section.
C           (WRB)
C   900328  Added TYPE section.  (WRB)
C***END PROLOGUE  D1MPYQ
      INTEGER I, J, LDA, M, N, NM1, NMJ
      DOUBLE PRECISION A(LDA,*), COS, ONE, SIN, TEMP, V(*), W(*)
      SAVE ONE
      DATA ONE /1.0D0/
C
C     APPLY THE FIRST SET OF GIVENS ROTATIONS TO A.
C
C***FIRST EXECUTABLE STATEMENT  D1MPYQ
      NM1 = N - 1
      IF (NM1 .LT. 1) GO TO 50
      DO 20 NMJ = 1, NM1
         J = N - NMJ
         IF (ABS(V(J)) .GT. ONE) COS = ONE/V(J)
         IF (ABS(V(J)) .GT. ONE) SIN = SQRT(ONE-COS**2)
         IF (ABS(V(J)) .LE. ONE) SIN = V(J)
         IF (ABS(V(J)) .LE. ONE) COS = SQRT(ONE-SIN**2)
         DO 10 I = 1, M
            TEMP = COS*A(I,J) - SIN*A(I,N)
            A(I,N) = SIN*A(I,J) + COS*A(I,N)
            A(I,J) = TEMP
   10       CONTINUE
   20    CONTINUE
C
C     APPLY THE SECOND SET OF GIVENS ROTATIONS TO A.
C
      DO 40 J = 1, NM1
         IF (ABS(W(J)) .GT. ONE) COS = ONE/W(J)
         IF (ABS(W(J)) .GT. ONE) SIN = SQRT(ONE-COS**2)
         IF (ABS(W(J)) .LE. ONE) SIN = W(J)
         IF (ABS(W(J)) .LE. ONE) COS = SQRT(ONE-SIN**2)
         DO 30 I = 1, M
            TEMP = COS*A(I,J) + SIN*A(I,N)
            A(I,N) = -SIN*A(I,J) + COS*A(I,N)
            A(I,J) = TEMP
   30       CONTINUE
   40    CONTINUE
   50 CONTINUE
      RETURN
C
C     LAST CARD OF SUBROUTINE D1MPYQ.
C
      END
*DECK D1UPDT
      SUBROUTINE D1UPDT (M, N, S, LS, U, V, W, SING)
C***BEGIN PROLOGUE  D1UPDT
C***SUBSIDIARY
C***PURPOSE  Subsidiary to DNSQ and DNSQE
C***LIBRARY   SLATEC
C***TYPE      DOUBLE PRECISION (R1UPDT-S, D1UPDT-D)
C***AUTHOR  (UNKNOWN)
C***DESCRIPTION
C
C     Given an M by N lower trapezoidal matrix S, an M-vector U,
C     and an N-vector V, the problem is to determine an
C     orthogonal matrix Q such that
C
C                   t
C           (S + U*V )*Q
C
C     is again lower trapezoidal.
C
C     This subroutine determines Q as the product of 2*(N - 1)
C     transformations
C
C           GV(N-1)*...*GV(1)*GW(1)*...*GW(N-1)
C
C     where GV(I), GW(I) are Givens rotations in the (I,N) plane
C     which eliminate elements in the I-th and N-th planes,
C     respectively. Q itself is not accumulated, rather the
C     information to recover the GV, GW rotations is returned.
C
C     The SUBROUTINE statement is
C
C       SUBROUTINE D1UPDT(M,N,S,LS,U,V,W,SING)
C
C     where
C
C       M is a positive integer input variable set to the number
C         of rows of S.
C
C       N is a positive integer input variable set to the number
C         of columns of S. N must not exceed M.
C
C       S is an array of length LS. On input S must contain the lower
C         trapezoidal matrix S stored by columns. On output S contains
C         the lower trapezoidal matrix produced as described above.
C
C       LS is a positive integer input variable not less than
C         (N*(2*M-N+1))/2.
C
C       U is an input array of length M which must contain the
C         vector U.
C
C       V is an array of length N. On input V must contain the vector
C         V. On output V(I) contains the information necessary to
C         recover the Givens rotation GV(I) described above.
C
C       W is an output array of length M. W(I) contains information
C         necessary to recover the Givens rotation GW(I) described
C         above.
C
C       SING is a LOGICAL output variable. SING is set TRUE if any
C         of the diagonal elements of the output S are zero. Otherwise
C         SING is set FALSE.
C
C***SEE ALSO  DNSQ, DNSQE
C***ROUTINES CALLED  D1MACH
C***REVISION HISTORY  (YYMMDD)
C   800301  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890831  Modified array declarations.  (WRB)
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900326  Removed duplicate information from DESCRIPTION section.
C           (WRB)
C   900328  Added TYPE section.  (WRB)
C***END PROLOGUE  D1UPDT
      DOUBLE PRECISION D1MACH
      INTEGER I, J, JJ, L, LS, M, N, NM1, NMJ
      DOUBLE PRECISION COS, COTAN, GIANT, ONE, P25, P5, S(*),
     1     SIN, TAN, TAU, TEMP, U(*), V(*), W(*), ZERO
      LOGICAL SING
      SAVE ONE, P5, P25, ZERO
      DATA ONE,P5,P25,ZERO /1.0D0,5.0D-1,2.5D-1,0.0D0/
C
C     GIANT IS THE LARGEST MAGNITUDE.
C
C***FIRST EXECUTABLE STATEMENT  D1UPDT
      GIANT = D1MACH(2)
C
C     INITIALIZE THE DIAGONAL ELEMENT POINTER.
C
      JJ = (N*(2*M - N + 1))/2 - (M - N)
C
C     MOVE THE NONTRIVIAL PART OF THE LAST COLUMN OF S INTO W.
C
      L = JJ
      DO 10 I = N, M
         W(I) = S(L)
         L = L + 1
   10    CONTINUE
C
C     ROTATE THE VECTOR V INTO A MULTIPLE OF THE N-TH UNIT VECTOR
C     IN SUCH A WAY THAT A SPIKE IS INTRODUCED INTO W.
C
      NM1 = N - 1
      IF (NM1 .LT. 1) GO TO 70
      DO 60 NMJ = 1, NM1
         J = N - NMJ
         JJ = JJ - (M - J + 1)
         W(J) = ZERO
         IF (V(J) .EQ. ZERO) GO TO 50
C
C        DETERMINE A GIVENS ROTATION WHICH ELIMINATES THE
C        J-TH ELEMENT OF V.
C
         IF (ABS(V(N)) .GE. ABS(V(J))) GO TO 20
            COTAN = V(N)/V(J)
            SIN = P5/SQRT(P25+P25*COTAN**2)
            COS = SIN*COTAN
            TAU = ONE
            IF (ABS(COS)*GIANT .GT. ONE) TAU = ONE/COS
            GO TO 30
   20    CONTINUE
            TAN = V(J)/V(N)
            COS = P5/SQRT(P25+P25*TAN**2)
            SIN = COS*TAN
            TAU = SIN
   30    CONTINUE
C
C        APPLY THE TRANSFORMATION TO V AND STORE THE INFORMATION
C        NECESSARY TO RECOVER THE GIVENS ROTATION.
C
         V(N) = SIN*V(J) + COS*V(N)
         V(J) = TAU
C
C        APPLY THE TRANSFORMATION TO S AND EXTEND THE SPIKE IN W.
C
         L = JJ
         DO 40 I = J, M
            TEMP = COS*S(L) - SIN*W(I)
            W(I) = SIN*S(L) + COS*W(I)
            S(L) = TEMP
            L = L + 1
   40       CONTINUE
   50    CONTINUE
   60    CONTINUE
   70 CONTINUE
C
C     ADD THE SPIKE FROM THE RANK 1 UPDATE TO W.
C
      DO 80 I = 1, M
         W(I) = W(I) + V(N)*U(I)
   80    CONTINUE
C
C     ELIMINATE THE SPIKE.
C
      SING = .FALSE.
      IF (NM1 .LT. 1) GO TO 140
      DO 130 J = 1, NM1
         IF (W(J) .EQ. ZERO) GO TO 120
C
C        DETERMINE A GIVENS ROTATION WHICH ELIMINATES THE
C        J-TH ELEMENT OF THE SPIKE.
C
         IF (ABS(S(JJ)) .GE. ABS(W(J))) GO TO 90
            COTAN = S(JJ)/W(J)
            SIN = P5/SQRT(P25+P25*COTAN**2)
            COS = SIN*COTAN
            TAU = ONE
            IF (ABS(COS)*GIANT .GT. ONE) TAU = ONE/COS
            GO TO 100
   90    CONTINUE
            TAN = W(J)/S(JJ)
            COS = P5/SQRT(P25+P25*TAN**2)
            SIN = COS*TAN
            TAU = SIN
  100    CONTINUE
C
C        APPLY THE TRANSFORMATION TO S AND REDUCE THE SPIKE IN W.
C
         L = JJ
         DO 110 I = J, M
            TEMP = COS*S(L) + SIN*W(I)
            W(I) = -SIN*S(L) + COS*W(I)
            S(L) = TEMP
            L = L + 1
  110       CONTINUE
C
C        STORE THE INFORMATION NECESSARY TO RECOVER THE
C        GIVENS ROTATION.
C
         W(J) = TAU
  120    CONTINUE
C
C        TEST FOR ZERO DIAGONAL ELEMENTS IN THE OUTPUT S.
C
         IF (S(JJ) .EQ. ZERO) SING = .TRUE.
         JJ = JJ + (M - J + 1)
  130    CONTINUE
  140 CONTINUE
C
C     MOVE W BACK INTO THE LAST COLUMN OF THE OUTPUT S.
C
      L = JJ
      DO 150 I = N, M
         S(L) = W(I)
         L = L + 1
  150    CONTINUE
      IF (S(JJ) .EQ. ZERO) SING = .TRUE.
      RETURN
C
C     LAST CARD OF SUBROUTINE D1UPDT.
C
      END
*DECK DDOGLG
      SUBROUTINE DDOGLG (N, R, LR, DIAG, QTB, DELTA, X, WA1, WA2)
C***BEGIN PROLOGUE  DDOGLG
C***SUBSIDIARY
C***PURPOSE  Subsidiary to DNSQ and DNSQE
C***LIBRARY   SLATEC
C***TYPE      DOUBLE PRECISION (DOGLEG-S, DDOGLG-D)
C***AUTHOR  (UNKNOWN)
C***DESCRIPTION
C
C     Given an M by N matrix A, an N by N nonsingular diagonal
C     matrix D, an M-vector B, and a positive number DELTA, the
C     problem is to determine the convex combination X of the
C     Gauss-Newton and scaled gradient directions that minimizes
C     (A*X - B) in the least squares sense, subject to the
C     restriction that the Euclidean norm of D*X be at most DELTA.
C
C     This subroutine completes the solution of the problem
C     if it is provided with the necessary information from the
C     QR factorization of A. That is, if A = Q*R, where Q has
C     orthogonal columns and R is an upper triangular matrix,
C     then DDOGLG expects the full upper triangle of R and
C     the first N components of (Q transpose)*B.
C
C     The subroutine statement is
C
C       SUBROUTINE DDOGLG(N,R,LR,DIAG,QTB,DELTA,X,WA1,WA2)
C
C     where
C
C       N is a positive integer input variable set to the order of R.
C
C       R is an input array of length LR which must contain the upper
C         triangular matrix R stored by rows.
C
C       LR is a positive integer input variable not less than
C         (N*(N+1))/2.
C
C       DIAG is an input array of length N which must contain the
C         diagonal elements of the matrix D.
C
C       QTB is an input array of length N which must contain the first
C         N elements of the vector (Q transpose)*B.
C
C       DELTA is a positive input variable which specifies an upper
C         bound on the Euclidean norm of D*X.
C
C       X is an output array of length N which contains the desired
C         convex combination of the Gauss-Newton direction and the
C         scaled gradient direction.
C
C       WA1 and WA2 are work arrays of length N.
C
C***SEE ALSO  DNSQ, DNSQE
C***ROUTINES CALLED  D1MACH, DENORM
C***REVISION HISTORY  (YYMMDD)
C   800301  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890831  Modified array declarations.  (WRB)
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900326  Removed duplicate information from DESCRIPTION section.
C           (WRB)
C   900328  Added TYPE section.  (WRB)
C***END PROLOGUE  DDOGLG
      DOUBLE PRECISION D1MACH,DENORM
      INTEGER I, J, JJ, JP1, K, L, LR, N
      DOUBLE PRECISION ALPHA, BNORM, DELTA, DIAG(*), EPSMCH, GNORM,
     1     ONE, QNORM, QTB(*), R(*), SGNORM, SUM, TEMP, WA1(*),
     2     WA2(*), X(*), ZERO
      SAVE ONE, ZERO
      DATA ONE,ZERO /1.0D0,0.0D0/
C
C     EPSMCH IS THE MACHINE PRECISION.
C
C***FIRST EXECUTABLE STATEMENT  DDOGLG
      EPSMCH = D1MACH(4)
C
C     FIRST, CALCULATE THE GAUSS-NEWTON DIRECTION.
C
      JJ = (N*(N + 1))/2 + 1
      DO 50 K = 1, N
         J = N - K + 1
         JP1 = J + 1
         JJ = JJ - K
         L = JJ + 1
         SUM = ZERO
         IF (N .LT. JP1) GO TO 20
         DO 10 I = JP1, N
            SUM = SUM + R(L)*X(I)
            L = L + 1
   10       CONTINUE
   20    CONTINUE
         TEMP = R(JJ)
         IF (TEMP .NE. ZERO) GO TO 40
         L = J
         DO 30 I = 1, J
            TEMP = MAX(TEMP,ABS(R(L)))
            L = L + N - I
   30       CONTINUE
         TEMP = EPSMCH*TEMP
         IF (TEMP .EQ. ZERO) TEMP = EPSMCH
   40    CONTINUE
         X(J) = (QTB(J) - SUM)/TEMP
   50    CONTINUE
C
C     TEST WHETHER THE GAUSS-NEWTON DIRECTION IS ACCEPTABLE.
C
      DO 60 J = 1, N
         WA1(J) = ZERO
         WA2(J) = DIAG(J)*X(J)
   60    CONTINUE
      QNORM = DENORM(N,WA2)
      IF (QNORM .LE. DELTA) GO TO 140
C
C     THE GAUSS-NEWTON DIRECTION IS NOT ACCEPTABLE.
C     NEXT, CALCULATE THE SCALED GRADIENT DIRECTION.
C
      L = 1
      DO 80 J = 1, N
         TEMP = QTB(J)
         DO 70 I = J, N
            WA1(I) = WA1(I) + R(L)*TEMP
            L = L + 1
   70       CONTINUE
         WA1(J) = WA1(J)/DIAG(J)
   80    CONTINUE
C
C     CALCULATE THE NORM OF THE SCALED GRADIENT AND TEST FOR
C     THE SPECIAL CASE IN WHICH THE SCALED GRADIENT IS ZERO.
C
      GNORM = DENORM(N,WA1)
      SGNORM = ZERO
      ALPHA = DELTA/QNORM
      IF (GNORM .EQ. ZERO) GO TO 120
C
C     CALCULATE THE POINT ALONG THE SCALED GRADIENT
C     AT WHICH THE QUADRATIC IS MINIMIZED.
C
      DO 90 J = 1, N
         WA1(J) = (WA1(J)/GNORM)/DIAG(J)
   90    CONTINUE
      L = 1
      DO 110 J = 1, N
         SUM = ZERO
         DO 100 I = J, N
            SUM = SUM + R(L)*WA1(I)
            L = L + 1
  100       CONTINUE
         WA2(J) = SUM
  110    CONTINUE
      TEMP = DENORM(N,WA2)
      SGNORM = (GNORM/TEMP)/TEMP
C
C     TEST WHETHER THE SCALED GRADIENT DIRECTION IS ACCEPTABLE.
C
      ALPHA = ZERO
      IF (SGNORM .GE. DELTA) GO TO 120
C
C     THE SCALED GRADIENT DIRECTION IS NOT ACCEPTABLE.
C     FINALLY, CALCULATE THE POINT ALONG THE DOGLEG
C     AT WHICH THE QUADRATIC IS MINIMIZED.
C
      BNORM = DENORM(N,QTB)
      TEMP = (BNORM/GNORM)*(BNORM/QNORM)*(SGNORM/DELTA)
      TEMP = TEMP - (DELTA/QNORM)*(SGNORM/DELTA)**2
     1       + SQRT((TEMP-(DELTA/QNORM))**2
     2               +(ONE-(DELTA/QNORM)**2)*(ONE-(SGNORM/DELTA)**2))
      ALPHA = ((DELTA/QNORM)*(ONE - (SGNORM/DELTA)**2))/TEMP
  120 CONTINUE
C
C     FORM APPROPRIATE CONVEX COMBINATION OF THE GAUSS-NEWTON
C     DIRECTION AND THE SCALED GRADIENT DIRECTION.
C
      TEMP = (ONE - ALPHA)*MIN(SGNORM,DELTA)
      DO 130 J = 1, N
         X(J) = TEMP*WA1(J) + ALPHA*X(J)
  130    CONTINUE
  140 CONTINUE
      RETURN
C
C     LAST CARD OF SUBROUTINE DDOGLG.
C
      END
*DECK DENORM
      DOUBLE PRECISION FUNCTION DENORM (N, X)
C***BEGIN PROLOGUE  DENORM
C***SUBSIDIARY
C***PURPOSE  Subsidiary to DNSQ and DNSQE
C***LIBRARY   SLATEC
C***TYPE      DOUBLE PRECISION (ENORM-S, DENORM-D)
C***AUTHOR  (UNKNOWN)
C***DESCRIPTION
C
C     Given an N-vector X, this function calculates the
C     Euclidean norm of X.
C
C     The Euclidean norm is computed by accumulating the sum of
C     squares in three different sums. The sums of squares for the
C     small and large components are scaled so that no overflows
C     occur. Non-destructive underflows are permitted. Underflows
C     and overflows do not occur in the computation of the unscaled
C     sum of squares for the intermediate components.
C     The definitions of small, intermediate and large components
C     depend on two constants, RDWARF and RGIANT. The main
C     restrictions on these constants are that RDWARF**2 not
C     underflow and RGIANT**2 not overflow. The constants
C     given here are suitable for every known computer.
C
C     The function statement is
C
C       DOUBLE PRECISION FUNCTION DENORM(N,X)
C
C     where
C
C       N is a positive integer input variable.
C
C       X is an input array of length N.
C
C***SEE ALSO  DNSQ, DNSQE
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   800301  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890831  Modified array declarations.  (WRB)
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900326  Removed duplicate information from DESCRIPTION section.
C           (WRB)
C   900328  Added TYPE section.  (WRB)
C***END PROLOGUE  DENORM
      INTEGER I, N
      DOUBLE PRECISION AGIANT, FLOATN, ONE, RDWARF, RGIANT, S1, S2, S3,
     1     X(*), X1MAX, X3MAX, XABS, ZERO
      SAVE ONE, ZERO, RDWARF, RGIANT
      DATA ONE,ZERO,RDWARF,RGIANT /1.0D0,0.0D0,3.834D-20,1.304D19/
C***FIRST EXECUTABLE STATEMENT  DENORM
      S1 = ZERO
      S2 = ZERO
      S3 = ZERO
      X1MAX = ZERO
      X3MAX = ZERO
      FLOATN = N
      AGIANT = RGIANT/FLOATN
      DO 90 I = 1, N
         XABS = ABS(X(I))
         IF (XABS .GT. RDWARF .AND. XABS .LT. AGIANT) GO TO 70
            IF (XABS .LE. RDWARF) GO TO 30
C
C              SUM FOR LARGE COMPONENTS.
C
               IF (XABS .LE. X1MAX) GO TO 10
                  S1 = ONE + S1*(X1MAX/XABS)**2
                  X1MAX = XABS
                  GO TO 20
   10          CONTINUE
                  S1 = S1 + (XABS/X1MAX)**2
   20          CONTINUE
               GO TO 60
   30       CONTINUE
C
C              SUM FOR SMALL COMPONENTS.
C
               IF (XABS .LE. X3MAX) GO TO 40
                  S3 = ONE + S3*(X3MAX/XABS)**2
                  X3MAX = XABS
                  GO TO 50
   40          CONTINUE
                  IF (XABS .NE. ZERO) S3 = S3 + (XABS/X3MAX)**2
   50          CONTINUE
   60       CONTINUE
            GO TO 80
   70    CONTINUE
C
C           SUM FOR INTERMEDIATE COMPONENTS.
C
            S2 = S2 + XABS**2
   80    CONTINUE
   90    CONTINUE
C
C     CALCULATION OF NORM.
C
      IF (S1 .EQ. ZERO) GO TO 100
         DENORM = X1MAX*SQRT(S1+(S2/X1MAX)/X1MAX)
         GO TO 130
  100 CONTINUE
         IF (S2 .EQ. ZERO) GO TO 110
            IF (S2 .GE. X3MAX)
     1         DENORM = SQRT(S2*(ONE+(X3MAX/S2)*(X3MAX*S3)))
            IF (S2 .LT. X3MAX)
     1         DENORM = SQRT(X3MAX*((S2/X3MAX)+(X3MAX*S3)))
            GO TO 120
  110    CONTINUE
            DENORM = X3MAX*SQRT(S3)
  120    CONTINUE
  130 CONTINUE
      RETURN
C
C     LAST CARD OF FUNCTION DENORM.
C
      END
*DECK DFDJC1
      SUBROUTINE DFDJC1 (FCN, N, X, FVEC, FJAC, LDFJAC, IFLAG, ML, MU,
     +   EPSFCN, WA1, WA2)
C***BEGIN PROLOGUE  DFDJC1
C***SUBSIDIARY
C***PURPOSE  Subsidiary to DNSQ and DNSQE
C***LIBRARY   SLATEC
C***TYPE      DOUBLE PRECISION (FDJAC1-S, DFDJC1-D)
C***AUTHOR  (UNKNOWN)
C***DESCRIPTION
C
C     This subroutine computes a forward-difference approximation
C     to the N by N Jacobian matrix associated with a specified
C     problem of N functions in N variables. If the Jacobian has
C     a banded form, then function evaluations are saved by only
C     approximating the nonzero terms.
C
C     The subroutine statement is
C
C       SUBROUTINE DFDJC1(FCN,N,X,FVEC,FJAC,LDFJAC,IFLAG,ML,MU,EPSFCN,
C                         WA1,WA2)
C
C     where
C
C       FCN is the name of the user-supplied subroutine which
C         calculates the functions. FCN must be declared
C         in an EXTERNAL statement in the user calling
C         program, and should be written as follows.
C
C         SUBROUTINE FCN(N,X,FVEC,IFLAG)
C         INTEGER N,IFLAG
C         DOUBLE PRECISION X(N),FVEC(N)
C         ----------
C         Calculate the functions at X and
C         return this vector in FVEC.
C         ----------
C         RETURN
C
C         The value of IFLAG should not be changed by FCN unless
C         the user wants to terminate execution of DFDJC1.
C         In this case set IFLAG to a negative integer.
C
C       N is a positive integer input variable set to the number
C         of functions and variables.
C
C       X is an input array of length N.
C
C       FVEC is an input array of length N which must contain the
C         functions evaluated at X.
C
C       FJAC is an output N by N array which contains the
C         approximation to the Jacobian matrix evaluated at X.
C
C       LDFJAC is a positive integer input variable not less than N
C         which specifies the leading dimension of the array FJAC.
C
C       IFLAG is an integer variable which can be used to terminate
C         the execution of DFDJC1. See description of FCN.
C
C       ML is a nonnegative integer input variable which specifies
C         the number of subdiagonals within the band of the
C         Jacobian matrix. If the Jacobian is not banded, set
C         ML to at least N - 1.
C
C       EPSFCN is an input variable used in determining a suitable
C         step length for the forward-difference approximation. This
C         approximation assumes that the relative errors in the
C         functions are of the order of EPSFCN. If EPSFCN is less
C         than the machine precision, it is assumed that the relative
C         errors in the functions are of the order of the machine
C         precision.
C
C       MU is a nonnegative integer input variable which specifies
C         the number of superdiagonals within the band of the
C         Jacobian matrix. If the Jacobian is not banded, set
C         MU to at least N - 1.
C
C       WA1 and WA2 are work arrays of length N. If ML + MU + 1 is at
C         least N, then the Jacobian is considered dense, and WA2 is
C         not referenced.
C
C***SEE ALSO  DNSQ, DNSQE
C***ROUTINES CALLED  D1MACH
C***REVISION HISTORY  (YYMMDD)
C   800301  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890831  Modified array declarations.  (WRB)
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900326  Removed duplicate information from DESCRIPTION section.
C           (WRB)
C   900328  Added TYPE section.  (WRB)
C***END PROLOGUE  DFDJC1
      DOUBLE PRECISION D1MACH
      INTEGER I, IFLAG, J, K, LDFJAC, ML, MSUM, MU, N
      DOUBLE PRECISION EPS, EPSFCN, EPSMCH, FJAC(LDFJAC,*),
     1     FVEC(*), H, TEMP, WA1(*), WA2(*), X(*), ZERO
      SAVE ZERO
      DATA ZERO /0.0D0/
C
C     EPSMCH IS THE MACHINE PRECISION.
C
C***FIRST EXECUTABLE STATEMENT  DFDJC1
      EPSMCH = D1MACH(4)
C
      EPS = SQRT(MAX(EPSFCN,EPSMCH))
      MSUM = ML + MU + 1
      IF (MSUM .LT. N) GO TO 40
C
C        COMPUTATION OF DENSE APPROXIMATE JACOBIAN.
C
         DO 20 J = 1, N
            TEMP = X(J)
            H = EPS*ABS(TEMP)
            IF (H .EQ. ZERO) H = EPS
            X(J) = TEMP + H
            CALL FCN(N,X,WA1,IFLAG)
            IF (IFLAG .LT. 0) GO TO 30
            X(J) = TEMP
            DO 10 I = 1, N
               FJAC(I,J) = (WA1(I) - FVEC(I))/H
   10          CONTINUE
   20       CONTINUE
   30    CONTINUE
         GO TO 110
   40 CONTINUE
C
C        COMPUTATION OF BANDED APPROXIMATE JACOBIAN.
C
         DO 90 K = 1, MSUM
            DO 60 J = K, N, MSUM
               WA2(J) = X(J)
               H = EPS*ABS(WA2(J))
               IF (H .EQ. ZERO) H = EPS
               X(J) = WA2(J) + H
   60          CONTINUE
            CALL FCN(N,X,WA1,IFLAG)
            IF (IFLAG .LT. 0) GO TO 100
            DO 80 J = K, N, MSUM
               X(J) = WA2(J)
               H = EPS*ABS(WA2(J))
               IF (H .EQ. ZERO) H = EPS
               DO 70 I = 1, N
                  FJAC(I,J) = ZERO
                  IF (I .GE. J - MU .AND. I .LE. J + ML)
     1               FJAC(I,J) = (WA1(I) - FVEC(I))/H
   70             CONTINUE
   80          CONTINUE
   90       CONTINUE
  100    CONTINUE
  110 CONTINUE
      RETURN
C
C     LAST CARD OF SUBROUTINE DFDJC1.
C
      END
*DECK DNSQ
      SUBROUTINE DNSQ (FCN, JAC, IOPT, N, X, FVEC, FJAC, LDFJAC, XTOL,
     +   MAXFEV, ML, MU, EPSFCN, DIAG, MODE, FACTOR, NPRINT, INFO, NFEV,
     +   NJEV, R, LR, QTF, WA1, WA2, WA3, WA4)
C***BEGIN PROLOGUE  DNSQ
C***PURPOSE  Find a zero of a system of a N nonlinear functions in N
C            variables by a modification of the Powell hybrid method.
C***LIBRARY   SLATEC
C***CATEGORY  F2A
C***TYPE      DOUBLE PRECISION (SNSQ-S, DNSQ-D)
C***KEYWORDS  NONLINEAR SQUARE SYSTEM, POWELL HYBRID METHOD, ZEROS
C***AUTHOR  Hiebert, K. L. (SNLA)
C***DESCRIPTION
C
C 1. Purpose.
C
C       The purpose of DNSQ is to find a zero of a system of N nonlinear
C       functions in N variables by a modification of the Powell
C       hybrid method.  The user must provide a subroutine which
C       calculates the functions.  The user has the option of either to
C       provide a subroutine which calculates the Jacobian or to let the
C       code calculate it by a forward-difference approximation.
C       This code is the combination of the MINPACK codes (Argonne)
C       HYBRD and HYBRDJ.
C
C 2. Subroutine and Type Statements.
C
C       SUBROUTINE DNSQ(FCN,JAC,IOPT,N,X,FVEC,FJAC,LDFJAC,XTOL,MAXFEV,
C      *                 ML,MU,EPSFCN,DIAG,MODE,FACTOR,NPRINT,INFO,NFEV,
C      *                 NJEV,R,LR,QTF,WA1,WA2,WA3,WA4)
C       INTEGER IOPT,N,MAXFEV,ML,MU,MODE,NPRINT,INFO,NFEV,LDFJAC,NJEV,LR
C       DOUBLE PRECISION XTOL,EPSFCN,FACTOR
C       DOUBLE PRECISION
C       X(N),FVEC(N),DIAG(N),FJAC(LDFJAC,N),R(LR),QTF(N),
C      *     WA1(N),WA2(N),WA3(N),WA4(N)
C       EXTERNAL FCN,JAC
C
C 3. Parameters.
C
C       Parameters designated as input parameters must be specified on
C       entry to DNSQ and are not changed on exit, while parameters
C       designated as output parameters need not be specified on entry
C       and are set to appropriate values on exit from DNSQ.
C
C       FCN is the name of the user-supplied subroutine which calculates
C         the functions.  FCN must be declared in an EXTERNAL statement
C         in the user calling program, and should be written as follows.
C
C         SUBROUTINE FCN(N,X,FVEC,IFLAG)
C         INTEGER N,IFLAG
C         DOUBLE PRECISION X(N),FVEC(N)
C         ----------
C         CALCULATE THE FUNCTIONS AT X AND
C         RETURN THIS VECTOR IN FVEC.
C         ----------
C         RETURN
C         END
C
C         The value of IFLAG should not be changed by FCN unless the
C         user wants to terminate execution of DNSQ.  In this case set
C         IFLAG to a negative integer.
C
C       JAC is the name of the user-supplied subroutine which calculates
C         the Jacobian.  If IOPT=1, then JAC must be declared in an
C         EXTERNAL statement in the user calling program, and should be
C         written as follows.
C
C         SUBROUTINE JAC(N,X,FVEC,FJAC,LDFJAC,IFLAG)
C         INTEGER N,LDFJAC,IFLAG
C         DOUBLE PRECISION X(N),FVEC(N),FJAC(LDFJAC,N)
C         ----------
C         Calculate the Jacobian at X and return this
C         matrix in FJAC.  FVEC contains the function
C         values at X and should not be altered.
C         ----------
C         RETURN
C         END
C
C         The value of IFLAG should not be changed by JAC unless the
C         user wants to terminate execution of DNSQ.  In this case set
C         IFLAG to a negative integer.
C
C         If IOPT=2, JAC can be ignored (treat it as a dummy argument).
C
C       IOPT is an input variable which specifies how the Jacobian will
C         be calculated.  If IOPT=1, then the user must supply the
C         Jacobian through the subroutine JAC.  If IOPT=2, then the
C         code will approximate the Jacobian by forward-differencing.
C
C       N is a positive integer input variable set to the number of
C         functions and variables.
C
C       X is an array of length N.  On input X must contain an initial
C         estimate of the solution vector.  On output X contains the
C         final estimate of the solution vector.
C
C       FVEC is an output array of length N which contains the functions
C         evaluated at the output X.
C
C       FJAC is an output N by N array which contains the orthogonal
C         matrix Q produced by the QR factorization of the final
C         approximate Jacobian.
C
C       LDFJAC is a positive integer input variable not less than N
C         which specifies the leading dimension of the array FJAC.
C
C       XTOL is a nonnegative input variable.  Termination occurs when
C         the relative error between two consecutive iterates is at most
C         XTOL.  Therefore, XTOL measures the relative error desired in
C         the approximate solution.  Section 4 contains more details
C         about XTOL.
C
C       MAXFEV is a positive integer input variable.  Termination occurs
C         when the number of calls to FCN is at least MAXFEV by the end
C         of an iteration.
C
C       ML is a nonnegative integer input variable which specifies the
C         number of subdiagonals within the band of the Jacobian matrix.
C         If the Jacobian is not banded or IOPT=1, set ML to at
C         least N - 1.
C
C       MU is a nonnegative integer input variable which specifies the
C         number of superdiagonals within the band of the Jacobian
C         matrix.  If the Jacobian is not banded or IOPT=1, set MU to at
C         least N - 1.
C
C       EPSFCN is an input variable used in determining a suitable step
C         for the forward-difference approximation.  This approximation
C         assumes that the relative errors in the functions are of the
C         order of EPSFCN.  If EPSFCN is less than the machine
C         precision, it is assumed that the relative errors in the
C         functions are of the order of the machine precision.  If
C         IOPT=1, then EPSFCN can be ignored (treat it as a dummy
C         argument).
C
C       DIAG is an array of length N.  If MODE = 1 (see below), DIAG is
C         internally set.  If MODE = 2, DIAG must contain positive
C         entries that serve as implicit (multiplicative) scale factors
C         for the variables.
C
C       MODE is an integer input variable.  If MODE = 1, the variables
C         will be scaled internally.  If MODE = 2, the scaling is
C         specified by the input DIAG.  Other values of MODE are
C         equivalent to MODE = 1.
C
C       FACTOR is a positive input variable used in determining the
C         initial step bound.  This bound is set to the product of
C         FACTOR and the Euclidean norm of DIAG*X if nonzero, or else to
C         FACTOR itself.  In most cases FACTOR should lie in the
C         interval (.1,100.).  100. is a generally recommended value.
C
C       NPRINT is an integer input variable that enables controlled
C         printing of iterates if it is positive.  In this case, FCN is
C         called with IFLAG = 0 at the beginning of the first iteration
C         and every NPRINT iterations thereafter and immediately prior
C         to return, with X and FVEC available for printing. appropriate
C         print statements must be added to FCN(see example).  If NPRINT
C         is not positive, no special calls of FCN with IFLAG = 0 are
C         made.
C
C       INFO is an integer output variable.  If the user has terminated
C         execution, INFO is set to the (negative) value of IFLAG.  See
C         description of FCN and JAC. Otherwise, INFO is set as follows.
C
C         INFO = 0  Improper input parameters.
C
C         INFO = 1  Relative error between two consecutive iterates is
C                   at most XTOL.
C
C         INFO = 2  Number of calls to FCN has reached or exceeded
C                   MAXFEV.
C
C         INFO = 3  XTOL is too small.  No further improvement in the
C                   approximate solution X is possible.
C
C         INFO = 4  Iteration is not making good progress, as measured
C                   by the improvement from the last five Jacobian
C                   evaluations.
C
C         INFO = 5  Iteration is not making good progress, as measured
C                   by the improvement from the last ten iterations.
C
C         Sections 4 and 5 contain more details about INFO.
C
C       NFEV is an integer output variable set to the number of calls to
C         FCN.
C
C       NJEV is an integer output variable set to the number of calls to
C         JAC. (If IOPT=2, then NJEV is set to zero.)
C
C       R is an output array of length LR which contains the upper
C         triangular matrix produced by the QR factorization of the
C         final approximate Jacobian, stored rowwise.
C
C       LR is a positive integer input variable not less than
C         (N*(N+1))/2.
C
C       QTF is an output array of length N which contains the vector
C         (Q transpose)*FVEC.
C
C       WA1, WA2, WA3, and WA4 are work arrays of length N.
C
C
C 4. Successful completion.
C
C       The accuracy of DNSQ is controlled by the convergence parameter
C       XTOL.  This parameter is used in a test which makes a comparison
C       between the approximation X and a solution XSOL.  DNSQ
C       terminates when the test is satisfied.  If the convergence
C       parameter is less than the machine precision (as defined by the
C       function D1MACH(4)), then DNSQ only attempts to satisfy the test
C       defined by the machine precision.  Further progress is not
C       usually possible.
C
C       The test assumes that the functions are reasonably well behaved,
C       and, if the Jacobian is supplied by the user, that the functions
C       and the Jacobian are coded consistently.  If these conditions
C       are not satisfied, then DNSQ may incorrectly indicate
C       convergence.  The coding of the Jacobian can be checked by the
C       subroutine DCKDER. If the Jacobian is coded correctly or IOPT=2,
C       then the validity of the answer can be checked, for example, by
C       rerunning DNSQ with a tighter tolerance.
C
C       Convergence Test.  If DENORM(Z) denotes the Euclidean norm of a
C         vector Z and D is the diagonal matrix whose entries are
C         defined by the array DIAG, then this test attempts to
C         guarantee that
C
C               DENORM(D*(X-XSOL)) .LE. XTOL*DENORM(D*XSOL).
C
C         If this condition is satisfied with XTOL = 10**(-K), then the
C         larger components of D*X have K significant decimal digits and
C         INFO is set to 1.  There is a danger that the smaller
C         components of D*X may have large relative errors, but the fast
C         rate of convergence of DNSQ usually avoids this possibility.
C         Unless high precision solutions are required, the recommended
C         value for XTOL is the square root of the machine precision.
C
C
C 5. Unsuccessful Completion.
C
C       Unsuccessful termination of DNSQ can be due to improper input
C       parameters, arithmetic interrupts, an excessive number of
C       function evaluations, or lack of good progress.
C
C       Improper Input Parameters.  INFO is set to 0 if IOPT .LT .1,
C         or IOPT .GT. 2, or N .LE. 0, or LDFJAC .LT. N, or
C         XTOL .LT. 0.E0, or MAXFEV .LE. 0, or ML .LT. 0, or MU .LT. 0,
C         or FACTOR .LE. 0.E0, or LR .LT. (N*(N+1))/2.
C
C       Arithmetic Interrupts.  If these interrupts occur in the FCN
C         subroutine during an early stage of the computation, they may
C         be caused by an unacceptable choice of X by DNSQ.  In this
C         case, it may be possible to remedy the situation by rerunning
C         DNSQ with a smaller value of FACTOR.
C
C       Excessive Number of Function Evaluations.  A reasonable value
C         for MAXFEV is 100*(N+1) for IOPT=1 and 200*(N+1) for IOPT=2.
C         If the number of calls to FCN reaches MAXFEV, then this
C         indicates that the routine is converging very slowly as
C         measured by the progress of FVEC, and INFO is set to 2. This
C         situation should be unusual because, as indicated below, lack
C         of good progress is usually diagnosed earlier by DNSQ,
C         causing termination with info = 4 or INFO = 5.
C
C       Lack of Good Progress.  DNSQ searches for a zero of the system
C         by minimizing the sum of the squares of the functions.  In so
C         doing, it can become trapped in a region where the minimum
C         does not correspond to a zero of the system and, in this
C         situation, the iteration eventually fails to make good
C         progress.  In particular, this will happen if the system does
C         not have a zero.  If the system has a zero, rerunning DNSQ
C         from a different starting point may be helpful.
C
C
C 6. Characteristics of The Algorithm.
C
C       DNSQ is a modification of the Powell Hybrid method.  Two of its
C       main characteristics involve the choice of the correction as a
C       convex combination of the Newton and scaled gradient directions,
C       and the updating of the Jacobian by the rank-1 method of
C       Broyden.  The choice of the correction guarantees (under
C       reasonable conditions) global convergence for starting points
C       far from the solution and a fast rate of convergence.  The
C       Jacobian is calculated at the starting point by either the
C       user-supplied subroutine or a forward-difference approximation,
C       but it is not recalculated until the rank-1 method fails to
C       produce satisfactory progress.
C
C       Timing.  The time required by DNSQ to solve a given problem
C         depends on N, the behavior of the functions, the accuracy
C         requested, and the starting point.  The number of arithmetic
C         operations needed by DNSQ is about 11.5*(N**2) to process
C         each evaluation of the functions (call to FCN) and 1.3*(N**3)
C         to process each evaluation of the Jacobian (call to JAC,
C         if IOPT = 1).  Unless FCN and JAC can be evaluated quickly,
C         the timing of DNSQ will be strongly influenced by the time
C         spent in FCN and JAC.
C
C       Storage.  DNSQ requires (3*N**2 + 17*N)/2 single precision
C         storage locations, in addition to the storage required by the
C         program.  There are no internally declared storage arrays.
C
C *Long Description:
C
C 7. Example.
C
C       The problem is to determine the values of X(1), X(2), ..., X(9),
C       which solve the system of tridiagonal equations
C
C       (3-2*X(1))*X(1)           -2*X(2)                   = -1
C               -X(I-1) + (3-2*X(I))*X(I)         -2*X(I+1) = -1, I=2-8
C                                   -X(8) + (3-2*X(9))*X(9) = -1
C C     **********
C
C       PROGRAM TEST
C C
C C     Driver for DNSQ example.
C C
C       INTEGER J,IOPT,N,MAXFEV,ML,MU,MODE,NPRINT,INFO,NFEV,LDFJAC,LR,
C      *        NWRITE
C       DOUBLE PRECISION XTOL,EPSFCN,FACTOR,FNORM
C       DOUBLE PRECISION X(9),FVEC(9),DIAG(9),FJAC(9,9),R(45),QTF(9),
C      *     WA1(9),WA2(9),WA3(9),WA4(9)
C       DOUBLE PRECISION DENORM,D1MACH
C       EXTERNAL FCN
C       DATA NWRITE /6/
C C
C       IOPT = 2
C       N = 9
C C
C C     THE FOLLOWING STARTING VALUES PROVIDE A ROUGH SOLUTION.
C C
C       DO 10 J = 1, 9
C          X(J) = -1.E0
C    10    CONTINUE
C C
C       LDFJAC = 9
C       LR = 45
C C
C C     SET XTOL TO THE SQUARE ROOT OF THE MACHINE PRECISION.
C C     UNLESS HIGH PRECISION SOLUTIONS ARE REQUIRED,
C C     THIS IS THE RECOMMENDED SETTING.
C C
C       XTOL = SQRT(D1MACH(4))
C C
C       MAXFEV = 2000
C       ML = 1
C       MU = 1
C       EPSFCN = 0.E0
C       MODE = 2
C       DO 20 J = 1, 9
C          DIAG(J) = 1.E0
C    20    CONTINUE
C       FACTOR = 1.E2
C       NPRINT = 0
C C
C       CALL DNSQ(FCN,JAC,IOPT,N,X,FVEC,FJAC,LDFJAC,XTOL,MAXFEV,ML,MU,
C      *           EPSFCN,DIAG,MODE,FACTOR,NPRINT,INFO,NFEV,NJEV,
C      *           R,LR,QTF,WA1,WA2,WA3,WA4)
C       FNORM = DENORM(N,FVEC)
C       WRITE (NWRITE,1000) FNORM,NFEV,INFO,(X(J),J=1,N)
C       STOP
C  1000 FORMAT (5X,' FINAL L2 NORM OF THE RESIDUALS',E15.7 //
C      *        5X,' NUMBER OF FUNCTION EVALUATIONS',I10 //
C      *        5X,' EXIT PARAMETER',16X,I10 //
C      *        5X,' FINAL APPROXIMATE SOLUTION' // (5X,3E15.7))
C       END
C       SUBROUTINE FCN(N,X,FVEC,IFLAG)
C       INTEGER N,IFLAG
C       DOUBLE PRECISION X(N),FVEC(N)
C       INTEGER K
C       DOUBLE PRECISION ONE,TEMP,TEMP1,TEMP2,THREE,TWO,ZERO
C       DATA ZERO,ONE,TWO,THREE /0.E0,1.E0,2.E0,3.E0/
C C
C       IF (IFLAG .NE. 0) GO TO 5
C C
C C     INSERT PRINT STATEMENTS HERE WHEN NPRINT IS POSITIVE.
C C
C       RETURN
C     5 CONTINUE
C       DO 10 K = 1, N
C          TEMP = (THREE - TWO*X(K))*X(K)
C          TEMP1 = ZERO
C          IF (K .NE. 1) TEMP1 = X(K-1)
C          TEMP2 = ZERO
C          IF (K .NE. N) TEMP2 = X(K+1)
C          FVEC(K) = TEMP - TEMP1 - TWO*TEMP2 + ONE
C    10    CONTINUE
C       RETURN
C       END
C
C       Results obtained with different compilers or machines
C       may be slightly different.
C
C       Final L2 norm of the residuals  0.1192636E-07
C
C       Number of function evaluations        14
C
C       Exit parameter                         1
C
C       Final approximate solution
C
C       -0.5706545E+00 -0.6816283E+00 -0.7017325E+00
C       -0.7042129E+00 -0.7013690E+00 -0.6918656E+00
C       -0.6657920E+00 -0.5960342E+00 -0.4164121E+00
C
C***REFERENCES  M. J. D. Powell, A hybrid method for nonlinear equa-
C                 tions. In Numerical Methods for Nonlinear Algebraic
C                 Equations, P. Rabinowitz, Editor.  Gordon and Breach,
C                 1988.
C***ROUTINES CALLED  D1MACH, D1MPYQ, D1UPDT, DDOGLG, DENORM, DFDJC1,
C                    DQFORM, DQRFAC, XERMSG
C***REVISION HISTORY  (YYMMDD)
C   800301  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890831  Modified array declarations.  (WRB)
C   890831  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  DNSQ
      DOUBLE PRECISION D1MACH,DENORM
      INTEGER I, IFLAG, INFO, IOPT, ITER, IWA(1), J, JM1, L, LDFJAC,
     1     LR, MAXFEV, ML, MODE, MU, N, NCFAIL, NCSUC, NFEV, NJEV,
     2     NPRINT, NSLOW1, NSLOW2
      DOUBLE PRECISION ACTRED, DELTA, DIAG(*), EPSFCN, EPSMCH, FACTOR,
     1     FJAC(LDFJAC,*), FNORM, FNORM1, FVEC(*), ONE, P0001, P001,
     2     P1, P5, PNORM, PRERED, QTF(*), R(*), RATIO, SUM, TEMP,
     3     WA1(*), WA2(*), WA3(*), WA4(*), X(*), XNORM, XTOL, ZERO
      EXTERNAL FCN
      LOGICAL JEVAL,SING
      SAVE ONE, P1, P5, P001, P0001, ZERO
      DATA ONE,P1,P5,P001,P0001,ZERO
     1     /1.0D0,1.0D-1,5.0D-1,1.0D-3,1.0D-4,0.0D0/
C
C     BEGIN BLOCK PERMITTING ...EXITS TO 320
C***FIRST EXECUTABLE STATEMENT  DNSQ
         EPSMCH = D1MACH(4)
C
         INFO = 0
         IFLAG = 0
         NFEV = 0
         NJEV = 0
C
C        CHECK THE INPUT PARAMETERS FOR ERRORS.
C
C     ...EXIT
         IF (IOPT .LT. 1 .OR. IOPT .GT. 2 .OR. N .LE. 0
     1       .OR. XTOL .LT. ZERO .OR. MAXFEV .LE. 0 .OR. ML .LT. 0
     2       .OR. MU .LT. 0 .OR. FACTOR .LE. ZERO .OR. LDFJAC .LT. N
     3       .OR. LR .LT. (N*(N + 1))/2) GO TO 320
         IF (MODE .NE. 2) GO TO 20
            DO 10 J = 1, N
C     .........EXIT
               IF (DIAG(J) .LE. ZERO) GO TO 320
   10       CONTINUE
   20    CONTINUE
C
C        EVALUATE THE FUNCTION AT THE STARTING POINT
C        AND CALCULATE ITS NORM.
C
         IFLAG = 1
         CALL FCN(N,X,FVEC,IFLAG)
         NFEV = 1
C     ...EXIT
         IF (IFLAG .LT. 0) GO TO 320
         FNORM = DENORM(N,FVEC)
C
C        INITIALIZE ITERATION COUNTER AND MONITORS.
C
         ITER = 1
         NCSUC = 0
         NCFAIL = 0
         NSLOW1 = 0
         NSLOW2 = 0
C
C        BEGINNING OF THE OUTER LOOP.
C
   30    CONTINUE
C           BEGIN BLOCK PERMITTING ...EXITS TO 90
               JEVAL = .TRUE.
C
C              CALCULATE THE JACOBIAN MATRIX.
C
               IF (IOPT .EQ. 2) GO TO 40
C
C                 USER SUPPLIES JACOBIAN
C
                  CALL JAC(N,X,FVEC,FJAC,LDFJAC,IFLAG)
                  NJEV = NJEV + 1
               GO TO 50
   40          CONTINUE
C
C                 CODE APPROXIMATES THE JACOBIAN
C
                  IFLAG = 2
                  CALL DFDJC1(FCN,N,X,FVEC,FJAC,LDFJAC,IFLAG,ML,MU,
     1                        EPSFCN,WA1,WA2)
                  NFEV = NFEV + MIN(ML+MU+1,N)
   50          CONTINUE
C
C     .........EXIT
               IF (IFLAG .LT. 0) GO TO 320
C
C              COMPUTE THE QR FACTORIZATION OF THE JACOBIAN.
C
               CALL DQRFAC(N,N,FJAC,LDFJAC,.FALSE.,IWA,1,WA1,WA2,WA3)
C
C              ON THE FIRST ITERATION AND IF MODE IS 1, SCALE ACCORDING
C              TO THE NORMS OF THE COLUMNS OF THE INITIAL JACOBIAN.
C
C           ...EXIT
               IF (ITER .NE. 1) GO TO 90
               IF (MODE .EQ. 2) GO TO 70
                  DO 60 J = 1, N
                     DIAG(J) = WA2(J)
                     IF (WA2(J) .EQ. ZERO) DIAG(J) = ONE
   60             CONTINUE
   70          CONTINUE
C
C              ON THE FIRST ITERATION, CALCULATE THE NORM OF THE SCALED
C              X AND INITIALIZE THE STEP BOUND DELTA.
C
               DO 80 J = 1, N
                  WA3(J) = DIAG(J)*X(J)
   80          CONTINUE
               XNORM = DENORM(N,WA3)
               DELTA = FACTOR*XNORM
               IF (DELTA .EQ. ZERO) DELTA = FACTOR
   90       CONTINUE
C
C           FORM (Q TRANSPOSE)*FVEC AND STORE IN QTF.
C
            DO 100 I = 1, N
               QTF(I) = FVEC(I)
  100       CONTINUE
            DO 140 J = 1, N
               IF (FJAC(J,J) .EQ. ZERO) GO TO 130
                  SUM = ZERO
                  DO 110 I = J, N
                     SUM = SUM + FJAC(I,J)*QTF(I)
  110             CONTINUE
                  TEMP = -SUM/FJAC(J,J)
                  DO 120 I = J, N
                     QTF(I) = QTF(I) + FJAC(I,J)*TEMP
  120             CONTINUE
  130          CONTINUE
  140       CONTINUE
C
C           COPY THE TRIANGULAR FACTOR OF THE QR FACTORIZATION INTO R.
C
            SING = .FALSE.
            DO 170 J = 1, N
               L = J
               JM1 = J - 1
               IF (JM1 .LT. 1) GO TO 160
               DO 150 I = 1, JM1
                  R(L) = FJAC(I,J)
                  L = L + N - I
  150          CONTINUE
  160          CONTINUE
               R(L) = WA1(J)
               IF (WA1(J) .EQ. ZERO) SING = .TRUE.
  170       CONTINUE
C
C           ACCUMULATE THE ORTHOGONAL FACTOR IN FJAC.
C
            CALL DQFORM(N,N,FJAC,LDFJAC,WA1)
C
C           RESCALE IF NECESSARY.
C
            IF (MODE .EQ. 2) GO TO 190
               DO 180 J = 1, N
                  DIAG(J) = MAX(DIAG(J),WA2(J))
  180          CONTINUE
  190       CONTINUE
C
C           BEGINNING OF THE INNER LOOP.
C
  200       CONTINUE
C
C              IF REQUESTED, CALL FCN TO ENABLE PRINTING OF ITERATES.
C
               IF (NPRINT .LE. 0) GO TO 210
                  IFLAG = 0
                  IF (MOD(ITER-1,NPRINT) .EQ. 0)
     1               CALL FCN(N,X,FVEC,IFLAG)
C     ............EXIT
                  IF (IFLAG .LT. 0) GO TO 320
  210          CONTINUE
C
C              DETERMINE THE DIRECTION P.
C
               CALL DDOGLG(N,R,LR,DIAG,QTF,DELTA,WA1,WA2,WA3)
C
C              STORE THE DIRECTION P AND X + P. CALCULATE THE NORM OF P.
C
               DO 220 J = 1, N
                  WA1(J) = -WA1(J)
                  WA2(J) = X(J) + WA1(J)
                  WA3(J) = DIAG(J)*WA1(J)
  220          CONTINUE
               PNORM = DENORM(N,WA3)
C
C              ON THE FIRST ITERATION, ADJUST THE INITIAL STEP BOUND.
C
               IF (ITER .EQ. 1) DELTA = MIN(DELTA,PNORM)
C
C              EVALUATE THE FUNCTION AT X + P AND CALCULATE ITS NORM.
C
               IFLAG = 1
               CALL FCN(N,WA2,WA4,IFLAG)
               NFEV = NFEV + 1
C     .........EXIT
               IF (IFLAG .LT. 0) GO TO 320
               FNORM1 = DENORM(N,WA4)
C
C              COMPUTE THE SCALED ACTUAL REDUCTION.
C
               ACTRED = -ONE
               IF (FNORM1 .LT. FNORM) ACTRED = ONE - (FNORM1/FNORM)**2
C
C              COMPUTE THE SCALED PREDICTED REDUCTION.
C
               L = 1
               DO 240 I = 1, N
                  SUM = ZERO
                  DO 230 J = I, N
                     SUM = SUM + R(L)*WA1(J)
                     L = L + 1
  230             CONTINUE
                  WA3(I) = QTF(I) + SUM
  240          CONTINUE
               TEMP = DENORM(N,WA3)
               PRERED = ZERO
               IF (TEMP .LT. FNORM) PRERED = ONE - (TEMP/FNORM)**2
C
C              COMPUTE THE RATIO OF THE ACTUAL TO THE PREDICTED
C              REDUCTION.
C
               RATIO = ZERO
               IF (PRERED .GT. ZERO) RATIO = ACTRED/PRERED
C
C              UPDATE THE STEP BOUND.
C
               IF (RATIO .GE. P1) GO TO 250
                  NCSUC = 0
                  NCFAIL = NCFAIL + 1
                  DELTA = P5*DELTA
               GO TO 260
  250          CONTINUE
                  NCFAIL = 0
                  NCSUC = NCSUC + 1
                  IF (RATIO .GE. P5 .OR. NCSUC .GT. 1)
     1               DELTA = MAX(DELTA,PNORM/P5)
                  IF (ABS(RATIO-ONE) .LE. P1) DELTA = PNORM/P5
  260          CONTINUE
C
C              TEST FOR SUCCESSFUL ITERATION.
C
               IF (RATIO .LT. P0001) GO TO 280
C
C                 SUCCESSFUL ITERATION. UPDATE X, FVEC, AND THEIR NORMS.
C
                  DO 270 J = 1, N
                     X(J) = WA2(J)
                     WA2(J) = DIAG(J)*X(J)
                     FVEC(J) = WA4(J)
  270             CONTINUE
                  XNORM = DENORM(N,WA2)
                  FNORM = FNORM1
                  ITER = ITER + 1
  280          CONTINUE
C
C              DETERMINE THE PROGRESS OF THE ITERATION.
C
               NSLOW1 = NSLOW1 + 1
               IF (ACTRED .GE. P001) NSLOW1 = 0
               IF (JEVAL) NSLOW2 = NSLOW2 + 1
               IF (ACTRED .GE. P1) NSLOW2 = 0
C
C              TEST FOR CONVERGENCE.
C
               IF (DELTA .LE. XTOL*XNORM .OR. FNORM .EQ. ZERO) INFO = 1
C     .........EXIT
               IF (INFO .NE. 0) GO TO 320
C
C              TESTS FOR TERMINATION AND STRINGENT TOLERANCES.
C
               IF (NFEV .GE. MAXFEV) INFO = 2
               IF (P1*MAX(P1*DELTA,PNORM) .LE. EPSMCH*XNORM) INFO = 3
               IF (NSLOW2 .EQ. 5) INFO = 4
               IF (NSLOW1 .EQ. 10) INFO = 5
C     .........EXIT
               IF (INFO .NE. 0) GO TO 320
C
C              CRITERION FOR RECALCULATING JACOBIAN
C
C           ...EXIT
               IF (NCFAIL .EQ. 2) GO TO 310
C
C              CALCULATE THE RANK ONE MODIFICATION TO THE JACOBIAN
C              AND UPDATE QTF IF NECESSARY.
C
               DO 300 J = 1, N
                  SUM = ZERO
                  DO 290 I = 1, N
                     SUM = SUM + FJAC(I,J)*WA4(I)
  290             CONTINUE
                  WA2(J) = (SUM - WA3(J))/PNORM
                  WA1(J) = DIAG(J)*((DIAG(J)*WA1(J))/PNORM)
                  IF (RATIO .GE. P0001) QTF(J) = SUM
  300          CONTINUE
C
C              COMPUTE THE QR FACTORIZATION OF THE UPDATED JACOBIAN.
C
               CALL D1UPDT(N,N,R,LR,WA1,WA2,WA3,SING)
               CALL D1MPYQ(N,N,FJAC,LDFJAC,WA2,WA3)
               CALL D1MPYQ(1,N,QTF,1,WA2,WA3)
C
C              END OF THE INNER LOOP.
C
               JEVAL = .FALSE.
            GO TO 200
  310       CONTINUE
C
C           END OF THE OUTER LOOP.
C
         GO TO 30
  320 CONTINUE
C
C     TERMINATION, EITHER NORMAL OR USER IMPOSED.
C
      IF (IFLAG .LT. 0) INFO = IFLAG
      IFLAG = 0
      IF (NPRINT .GT. 0) CALL FCN(N,X,FVEC,IFLAG)
      IF (INFO .LT. 0) CALL XERMSG ('SLATEC', 'DNSQ',
     +   'EXECUTION TERMINATED BECAUSE USER SET IFLAG NEGATIVE.', 1, 1)
      IF (INFO .EQ. 0) CALL XERMSG ('SLATEC', 'DNSQ',
     +   'INVALID INPUT PARAMETER.', 2, 1)
      IF (INFO .EQ. 2) CALL XERMSG ('SLATEC', 'DNSQ',
     +   'TOO MANY FUNCTION EVALUATIONS.', 9, 1)
      IF (INFO .EQ. 3) CALL XERMSG ('SLATEC', 'DNSQ',
     +   'XTOL TOO SMALL. NO FURTHER IMPROVEMENT POSSIBLE.', 3, 1)
      IF (INFO .GT. 4) CALL XERMSG ('SLATEC', 'DNSQ',
     +   'ITERATION NOT MAKING GOOD PROGRESS.', 1, 1)
      RETURN
C
C     LAST CARD OF SUBROUTINE DNSQ.
C
      END
*DECK DQFORM
      SUBROUTINE DQFORM (M, N, Q, LDQ, WA)
C***BEGIN PROLOGUE  DQFORM
C***SUBSIDIARY
C***PURPOSE  Subsidiary to DNSQ and DNSQE
C***LIBRARY   SLATEC
C***TYPE      DOUBLE PRECISION (QFORM-S, DQFORM-D)
C***AUTHOR  (UNKNOWN)
C***DESCRIPTION
C
C     This subroutine proceeds from the computed QR factorization of
C     an M by N matrix A to accumulate the M by M orthogonal matrix
C     Q from its factored form.
C
C     The subroutine statement is
C
C       SUBROUTINE DQFORM(M,N,Q,LDQ,WA)
C
C     where
C
C       M is a positive integer input variable set to the number
C         of rows of A and the order of Q.
C
C       N is a positive integer input variable set to the number
C         of columns of A.
C
C       Q is an M by M array. On input the full lower trapezoid in
C         the first MIN(M,N) columns of Q contains the factored form.
C         On output Q has been accumulated into a square matrix.
C
C       LDQ is a positive integer input variable not less than M
C         which specifies the leading dimension of the array Q.
C
C       WA is a work array of length M.
C
C***SEE ALSO  DNSQ, DNSQE
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   800301  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890831  Modified array declarations.  (WRB)
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900326  Removed duplicate information from DESCRIPTION section.
C           (WRB)
C   900328  Added TYPE section.  (WRB)
C***END PROLOGUE  DQFORM
      INTEGER I, J, JM1, K, L, LDQ, M, MINMN, N, NP1
      DOUBLE PRECISION ONE, Q(LDQ,*), SUM, TEMP, WA(*), ZERO
      SAVE ONE, ZERO
      DATA ONE,ZERO /1.0D0,0.0D0/
C
C     ZERO OUT UPPER TRIANGLE OF Q IN THE FIRST MIN(M,N) COLUMNS.
C
C***FIRST EXECUTABLE STATEMENT  DQFORM
      MINMN = MIN(M,N)
      IF (MINMN .LT. 2) GO TO 30
      DO 20 J = 2, MINMN
         JM1 = J - 1
         DO 10 I = 1, JM1
            Q(I,J) = ZERO
   10       CONTINUE
   20    CONTINUE
   30 CONTINUE
C
C     INITIALIZE REMAINING COLUMNS TO THOSE OF THE IDENTITY MATRIX.
C
      NP1 = N + 1
      IF (M .LT. NP1) GO TO 60
      DO 50 J = NP1, M
         DO 40 I = 1, M
            Q(I,J) = ZERO
   40       CONTINUE
         Q(J,J) = ONE
   50    CONTINUE
   60 CONTINUE
C
C     ACCUMULATE Q FROM ITS FACTORED FORM.
C
      DO 120 L = 1, MINMN
         K = MINMN - L + 1
         DO 70 I = K, M
            WA(I) = Q(I,K)
            Q(I,K) = ZERO
   70       CONTINUE
         Q(K,K) = ONE
         IF (WA(K) .EQ. ZERO) GO TO 110
         DO 100 J = K, M
            SUM = ZERO
            DO 80 I = K, M
               SUM = SUM + Q(I,J)*WA(I)
   80          CONTINUE
            TEMP = SUM/WA(K)
            DO 90 I = K, M
               Q(I,J) = Q(I,J) - TEMP*WA(I)
   90          CONTINUE
  100       CONTINUE
  110    CONTINUE
  120    CONTINUE
      RETURN
C
C     LAST CARD OF SUBROUTINE DQFORM.
C
      END
*DECK DQRFAC
      SUBROUTINE DQRFAC (M, N, A, LDA, PIVOT, IPVT, LIPVT, SIGMA,
     +   ACNORM, WA)
C***BEGIN PROLOGUE  DQRFAC
C***SUBSIDIARY
C***PURPOSE  Subsidiary to DNLS1, DNLS1E, DNSQ and DNSQE
C***LIBRARY   SLATEC
C***TYPE      DOUBLE PRECISION (QRFAC-S, DQRFAC-D)
C***AUTHOR  (UNKNOWN)
C***DESCRIPTION
C
C   **** Double Precision version of QRFAC ****
C
C     This subroutine uses Householder transformations with column
C     pivoting (optional) to compute a QR factorization of the
C     M by N matrix A. That is, DQRFAC determines an orthogonal
C     matrix Q, a permutation matrix P, and an upper trapezoidal
C     matrix R with diagonal elements of nonincreasing magnitude,
C     such that A*P = Q*R. The Householder transformation for
C     column K, K = 1,2,...,MIN(M,N), is of the form
C
C                           T
C           I - (1/U(K))*U*U
C
C     where U has zeros in the first K-1 positions. The form of
C     this transformation and the method of pivoting first
C     appeared in the corresponding LINPACK subroutine.
C
C     The subroutine statement is
C
C       SUBROUTINE DQRFAC(M,N,A,LDA,PIVOT,IPVT,LIPVT,SIGMA,ACNORM,WA)
C
C     where
C
C       M is a positive integer input variable set to the number
C         of rows of A.
C
C       N is a positive integer input variable set to the number
C         of columns of A.
C
C       A is an M by N array. On input A contains the matrix for
C         which the QR factorization is to be computed. On output
C         the strict upper trapezoidal part of A contains the strict
C         upper trapezoidal part of R, and the lower trapezoidal
C         part of A contains a factored form of Q (the non-trivial
C         elements of the U vectors described above).
C
C       LDA is a positive integer input variable not less than M
C         which specifies the leading dimension of the array A.
C
C       PIVOT is a logical input variable. If pivot is set .TRUE.,
C         then column pivoting is enforced. If pivot is set .FALSE.,
C         then no column pivoting is done.
C
C       IPVT is an integer output array of length LIPVT. IPVT
C         defines the permutation matrix P such that A*P = Q*R.
C         Column J of P is column IPVT(J) of the identity matrix.
C         If pivot is .FALSE., IPVT is not referenced.
C
C       LIPVT is a positive integer input variable. If PIVOT is
C             .FALSE., then LIPVT may be as small as 1. If PIVOT is
C             .TRUE., then LIPVT must be at least N.
C
C       SIGMA is an output array of length N which contains the
C         diagonal elements of R.
C
C       ACNORM is an output array of length N which contains the
C         norms of the corresponding columns of the input matrix A.
C         If this information is not needed, then ACNORM can coincide
C         with SIGMA.
C
C       WA is a work array of length N. If pivot is .FALSE., then WA
C         can coincide with SIGMA.
C
C***SEE ALSO  DNLS1, DNLS1E, DNSQ, DNSQE
C***ROUTINES CALLED  D1MACH, DENORM
C***REVISION HISTORY  (YYMMDD)
C   800301  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890831  Modified array declarations.  (WRB)
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900326  Removed duplicate information from DESCRIPTION section.
C           (WRB)
C   900328  Added TYPE section.  (WRB)
C***END PROLOGUE  DQRFAC
      INTEGER M,N,LDA,LIPVT
      INTEGER IPVT(*)
      LOGICAL PIVOT
      SAVE ONE, P05, ZERO
      DOUBLE PRECISION A(LDA,*),SIGMA(*),ACNORM(*),WA(*)
      INTEGER I,J,JP1,K,KMAX,MINMN
      DOUBLE PRECISION AJNORM,EPSMCH,ONE,P05,SUM,TEMP,ZERO
      DOUBLE PRECISION D1MACH,DENORM
      DATA ONE,P05,ZERO /1.0D0,5.0D-2,0.0D0/
C***FIRST EXECUTABLE STATEMENT  DQRFAC
      EPSMCH = D1MACH(4)
C
C     COMPUTE THE INITIAL COLUMN NORMS AND INITIALIZE SEVERAL ARRAYS.
C
      DO 10 J = 1, N
         ACNORM(J) = DENORM(M,A(1,J))
         SIGMA(J) = ACNORM(J)
         WA(J) = SIGMA(J)
         IF (PIVOT) IPVT(J) = J
   10    CONTINUE
C
C     REDUCE A TO R WITH HOUSEHOLDER TRANSFORMATIONS.
C
      MINMN = MIN(M,N)
      DO 110 J = 1, MINMN
         IF (.NOT.PIVOT) GO TO 40
C
C        BRING THE COLUMN OF LARGEST NORM INTO THE PIVOT POSITION.
C
         KMAX = J
         DO 20 K = J, N
            IF (SIGMA(K) .GT. SIGMA(KMAX)) KMAX = K
   20       CONTINUE
         IF (KMAX .EQ. J) GO TO 40
         DO 30 I = 1, M
            TEMP = A(I,J)
            A(I,J) = A(I,KMAX)
            A(I,KMAX) = TEMP
   30       CONTINUE
         SIGMA(KMAX) = SIGMA(J)
         WA(KMAX) = WA(J)
         K = IPVT(J)
         IPVT(J) = IPVT(KMAX)
         IPVT(KMAX) = K
   40    CONTINUE
C
C        COMPUTE THE HOUSEHOLDER TRANSFORMATION TO REDUCE THE
C        J-TH COLUMN OF A TO A MULTIPLE OF THE J-TH UNIT VECTOR.
C
         AJNORM = DENORM(M-J+1,A(J,J))
         IF (AJNORM .EQ. ZERO) GO TO 100
         IF (A(J,J) .LT. ZERO) AJNORM = -AJNORM
         DO 50 I = J, M
            A(I,J) = A(I,J)/AJNORM
   50       CONTINUE
         A(J,J) = A(J,J) + ONE
C
C        APPLY THE TRANSFORMATION TO THE REMAINING COLUMNS
C        AND UPDATE THE NORMS.
C
         JP1 = J + 1
         IF (N .LT. JP1) GO TO 100
         DO 90 K = JP1, N
            SUM = ZERO
            DO 60 I = J, M
               SUM = SUM + A(I,J)*A(I,K)
   60          CONTINUE
            TEMP = SUM/A(J,J)
            DO 70 I = J, M
               A(I,K) = A(I,K) - TEMP*A(I,J)
   70          CONTINUE
            IF (.NOT.PIVOT .OR. SIGMA(K) .EQ. ZERO) GO TO 80
            TEMP = A(J,K)/SIGMA(K)
            SIGMA(K) = SIGMA(K)*SQRT(MAX(ZERO,ONE-TEMP**2))
            IF (P05*(SIGMA(K)/WA(K))**2 .GT. EPSMCH) GO TO 80
            SIGMA(K) = DENORM(M-J,A(JP1,K))
            WA(K) = SIGMA(K)
   80       CONTINUE
   90       CONTINUE
  100    CONTINUE
         SIGMA(J) = -AJNORM
  110    CONTINUE
      RETURN
C
C     LAST CARD OF SUBROUTINE DQRFAC.
C
      END
*DECK XERMSG
      SUBROUTINE XERMSG (LIBRAR, SUBROU, MESSG, NERR, LEVEL)
C-----------------------------------------------------------------------
C Subroutines XERMSG, XSETF, XSETUN, and the function routine IXSAV, as
C given here, constitute a simplified version of the SLATEC error 
C handling package.  Written by A. C. Hindmarsh, 18 November 1992.
C
C All arguments are input arguments.
C LIBRAR = Library name (character array).  Prefixed to message.
C SUBROU = Routine name (character array).  Prefixed to message.
C MESSG  = The message (character array).
C NERR   = Integer error number.  Prefixed to message.
C LEVEL  = The error level..
C          0 or 1 means recoverable (control returns to caller).
C          2 means fatal (run is aborted--see note below).
C
C Note..  This routine has been simplified in the following ways..
C 1. A single prefix line is printed with NERR, SUBROU, and LIBRAR.
C 2. The message in MESSG is printed, unaltered, on lines of up to 72
C    characters each using a format of (A).
C 3. If LEVEL = 2, control passes to the statement   STOP
C    to abort the run.  This statement may be machine-dependent.
C
C For a different default logical unit number, change the data
C statement in function routine IXSAV.
C For a different run-abort command, change the statement following
C statement 100 at the end.
C-----------------------------------------------------------------------
C Subroutines called by XERMSG.. None
C Function routines called by XERMSG.. IXSAV
C Intrinsic function used by XERMSG.. LEN
C-----------------------------------------------------------------------
      CHARACTER*(*) LIBRAR, SUBROU, MESSG
      INTEGER NERR, LEVEL
      INTEGER I1, I2, IL, IXSAV, LENMSG, LLEN, LUNIT, MESFLG, NLINES
      PARAMETER (LLEN = 72)
C
C Get message print flag and logical unit number. ----------------------
      LUNIT = IXSAV (1, 0, .FALSE.)
      MESFLG = IXSAV (2, 0, .FALSE.)
      IF (MESFLG .EQ. 0) GO TO 100
C Write NERR, SUBROU, and LIBRAR. --------------------------------------
      I1 = LEN(SUBROU)
      I2 = LEN(LIBRAR)
      WRITE (LUNIT, 10) NERR, SUBROU(1:I1), LIBRAR(1:I2)
  10  FORMAT(/,'***Error number ',I6,' from ',A,' in library ',A,'***')
C Write the message. ---------------------------------------------------
      LENMSG = LEN(MESSG)
      NLINES = ( (LENMSG - 1)/LLEN ) + 1
      DO 20 IL = 1,NLINES
        I1 = 1 + (IL - 1)*LLEN
        I2 = MIN(IL*LLEN,LENMSG)
        WRITE (LUNIT,'(A)') MESSG(I1:I2)
  20    CONTINUE
C Abort the run if LEVEL = 2. ------------------------------------------
 100  IF (LEVEL .NE. 2) RETURN
      STOP
C----------------------- End of Subroutine XERMSG ----------------------
      END
*DECK XSETUN
      SUBROUTINE XSETUN (LUN)
C-----------------------------------------------------------------------
C This routine resets the logical unit number for messages.
C
C Subroutines called by XSETUN.. None
C Function routine called by XSETUN.. IXSAV
C-----------------------------------------------------------------------
      INTEGER LUN, JUNK, IXSAV
C
      IF (LUN .GT. 0) JUNK = IXSAV (1,LUN,.TRUE.)
      RETURN
C----------------------- End of Subroutine XSETUN ----------------------
      END
*DECK XSETF
      SUBROUTINE XSETF (MFLAG)
C-----------------------------------------------------------------------
C This routine resets the print control flag MFLAG.
C
C Subroutines called by XSETF.. None
C Function routine called by XSETF.. IXSAV
C-----------------------------------------------------------------------
      INTEGER MFLAG, JUNK, IXSAV
C
      IF (MFLAG .EQ. 0 .OR. MFLAG .EQ. 1) JUNK = IXSAV (2,MFLAG,.TRUE.)
      RETURN
C----------------------- End of Subroutine XSETF -----------------------
      END
*DECK IXSAV
      INTEGER FUNCTION IXSAV (IPAR, IVALUE, ISET)
      LOGICAL ISET
      INTEGER IPAR, IVALUE
C-----------------------------------------------------------------------
C IXSAV saves and recalls one of two error message parameters:
C   LUNIT, the logical unit number to which messages are printed, and
C   MESFLG, the message print flag.
C This is a modification of the SLATEC library routine J4SAVE.
C
C Saved local variables..
C  LUNIT  = Logical unit number for messages.
C           The default is 6 (machine-dependent).
C  MESFLG = Print control flag..
C           1 means print all messages (the default).
C           0 means no printing.
C
C On input..
C   IPAR   = Parameter indicator (1 for LUNIT, 2 for MESFLG).
C   IVALUE = The value to be set for the parameter, if ISET = .TRUE.
C   ISET   = Logical flag to indicate whether to read or write.
C            If ISET = .TRUE., the parameter will be given
C            the value IVALUE.  If ISET = .FALSE., the parameter
C            will be unchanged, and IVALUE is a dummy argument.
C
C On return..
C   IXSAV = The (old) value of the parameter.
C
C Subroutines/functions called by IXSAV.. None
C-----------------------------------------------------------------------
      INTEGER LUNIT, MESFLG
C-----------------------------------------------------------------------
C The following Fortran-77 declaration is to cause the values of the
C listed (local) variables to be saved between calls to this routine.
C-----------------------------------------------------------------------
      SAVE LUNIT, MESFLG
      DATA LUNIT/6/, MESFLG/1/
C
      IF (IPAR .EQ. 1) THEN
        IXSAV = LUNIT
        IF (ISET) LUNIT = IVALUE
        ENDIF
C
      IF (IPAR .EQ. 2) THEN
        IXSAV = MESFLG
        IF (ISET) MESFLG = IVALUE
        ENDIF
C
      RETURN
C----------------------- End of Function IXSAV -------------------------
      END
