      subroutine PhotoConstants_Fixed
      use module_data_mosaic_main
      use module_data_mosaic_gas

      implicit none

      rk_photo(1)   = 0.0	! no2 + hv    --> no + o(3P)
      rk_photo(2)   = 0.0	! no3 + hv    --> .89no2 + .89o(3P) + .11no
      rk_photo(3)   = 0.0	! hono + hv   --> oh + no
      rk_photo(4)   = 0.0	! hno3 + hv   --> oh + no2
      rk_photo(5)   = 0.0	! hno4 + hv   --> ho2 + no2
      rk_photo(6)   = 0.0	! n2o5 + hv   --> no2 + no3
      rk_photo(7)   = 0.0	! o3 + hv     --> o(3P)
      rk_photo(8)   = 0.0	! o3 + hv     --> o(1D)
      rk_photo(9)   = 0.0	! h2o2 + hv   --> 2oh
      rk_photo(10)  = 0.0	! hcho + hv   --> 2ho2 + co
      rk_photo(11)  = 0.0	! hcho + hv   --> co
      rk_photo(12)  = 0.0	! ch3ooh + hv --> hcho + ho2 + oh
      rk_photo(13)  = 0.0	! ethooh + hv --> ald2 + ho2 + oh
      rk_photo(14)  = 0.0	! ald2 + hv   -->
      rk_photo(15)  = 0.0	! aone + hv   -->
      rk_photo(16)  = 0.0	! mgly + hv   -->
      rk_photo(17)  = 0.0	! open + hv   -->
      rk_photo(18)  = 0.0	! rooh + hv   -->
      rk_photo(19)  = 0.0	! onit + hv   -->
      rk_photo(20)  = 0.0	! isoprd + hv -->

      return
      end subroutine PhotoConstants_Fixed
