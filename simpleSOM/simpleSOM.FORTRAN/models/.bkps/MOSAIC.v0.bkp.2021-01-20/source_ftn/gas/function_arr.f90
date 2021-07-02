      function ARR(AA,BB)
      use module_data_mosaic_main

      implicit none

      real(r8) :: ARR
      real(r8) :: AA, BB

      ARR = AA*exp(BB/te)

      return
      end function ARR
