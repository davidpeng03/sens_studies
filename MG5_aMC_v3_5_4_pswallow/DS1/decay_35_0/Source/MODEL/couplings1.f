ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      written by the UFO converter
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      SUBROUTINE COUP1()

      IMPLICIT NONE
      INCLUDE 'model_functions.inc'

      DOUBLE PRECISION PI, ZERO
      PARAMETER  (PI=3.141592653589793D0)
      PARAMETER  (ZERO=0D0)
      INCLUDE 'input.inc'
      INCLUDE 'coupl.inc'
      GC_70 = -((MDL_COMPLEXI*MDL_SH*MDL_YB)/MDL_SQRT__2)
      GC_72 = -((MDL_COMPLEXI*MDL_SH*MDL_YC)/MDL_SQRT__2)
      GC_74 = -((MDL_COMPLEXI*MDL_SH*MDL_YE)/MDL_SQRT__2)
      GC_76 = -((MDL_COMPLEXI*MDL_SH*MDL_YM)/MDL_SQRT__2)
      GC_80 = -((MDL_COMPLEXI*MDL_SH*MDL_YTAU)/MDL_SQRT__2)
      END
