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
      GC_23 = (MDL_CKM1X1*MDL_EE*MDL_COMPLEXI)/(MDL_SW*MDL_SQRT__2)
      GC_24 = (MDL_CKM1X2*MDL_EE*MDL_COMPLEXI)/(MDL_SW*MDL_SQRT__2)
      GC_26 = (MDL_CKM2X1*MDL_EE*MDL_COMPLEXI)/(MDL_SW*MDL_SQRT__2)
      GC_27 = (MDL_CKM2X2*MDL_EE*MDL_COMPLEXI)/(MDL_SW*MDL_SQRT__2)
      GC_49 = (MDL_EE*MDL_COMPLEXI*MDL_VEN1)/(MDL_SW*MDL_SQRT__2)
      GC_179 = (MDL_EE*MDL_COMPLEXI*MDL_CONJG__CKM1X1)/(MDL_SW
     $ *MDL_SQRT__2)
      GC_182 = (MDL_EE*MDL_COMPLEXI*MDL_CONJG__CKM1X2)/(MDL_SW
     $ *MDL_SQRT__2)
      GC_188 = (MDL_EE*MDL_COMPLEXI*MDL_CONJG__CKM2X1)/(MDL_SW
     $ *MDL_SQRT__2)
      GC_191 = (MDL_EE*MDL_COMPLEXI*MDL_CONJG__CKM2X2)/(MDL_SW
     $ *MDL_SQRT__2)
      END
