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
      GC_67 = -(MDL_CH__EXP__3*MDL_COMPLEXI*MDL_KAP*MDL_V)+2.000000D
     $ +00*MDL_CH*MDL_COMPLEXI*MDL_KAP*MDL_SH__EXP__2*MDL_V-6.000000D
     $ +00*MDL_CH*MDL_COMPLEXI*MDL_LAM*MDL_SH__EXP__2*MDL_V-2.000000D
     $ +00*MDL_CH__EXP__2*MDL_COMPLEXI*MDL_KAP*MDL_SH*MDL_XI+6.000000D
     $ +00*MDL_CH__EXP__2*MDL_COMPLEXI*MDL_RHO*MDL_SH*MDL_XI
     $ +MDL_COMPLEXI*MDL_KAP*MDL_SH__EXP__3*MDL_XI
      GC_71 = -((MDL_CH*MDL_COMPLEXI*MDL_YC)/MDL_SQRT__2)
      END
