*DECK D1MACH
      DOUBLE PRECISION FUNCTION D1MACH (IDUMMY)
C***BEGIN PROLOGUE  D1MACH
C***PURPOSE  Compute the unit roundoff of the machine.
C***CATEGORY  R1
C***TYPE      DOUBLE PRECISION (R1MACH-S, D1MACH-D)
C***KEYWORDS  MACHINE CONSTANTS
C***AUTHOR  Hindmarsh, Alan C., (LLNL)
C***DESCRIPTION
C *Usage:
C        DOUBLE PRECISION  A, D1MACH
C        A = D1MACH(idummy)  [The argument is ignored.]
C
C *Function Return Values:
C     A : the unit roundoff of the machine.
C
C *Description:
C     The unit roundoff is defined as the smallest positive machine
C     number u such that  1.0 + u .ne. 1.0.  This is computed by D1MACH
C     in a machine-independent manner.
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  DUMSUM
C***REVISION HISTORY  (YYYYMMDD)
C   19930216  DATE WRITTEN
C   19930818  Added SLATEC-format prologue.  (FNF)
C   20030707  Added DUMSUM to force normal storage of COMP.  (ACH)
C***END PROLOGUE  D1MACH
C
      INTEGER IDUMMY
      DOUBLE PRECISION U, COMP
C***FIRST EXECUTABLE STATEMENT  D1MACH
      U = 1.0D0
 10   U = U*0.5D0
      CALL DUMSUM(1.0D0, U, COMP)
      IF (COMP .NE. 1.0D0) GO TO 10
      D1MACH = U*2.0D0
      RETURN
C----------------------- End of Function D1MACH ------------------------
      END