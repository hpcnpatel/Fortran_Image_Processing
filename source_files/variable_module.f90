MODULE variable_module

USE CONSTANTS_MODULE

!** USER defined data types are declared in this module **!

      TYPE textent
        INTEGER                   ::a
        INTEGER                   ::e
      END TYPE textent

      TYPE tCt
        INTEGER                    ::pixel_x    ! No. of pixels in x-direction
        INTEGER                    ::pixel_y    ! No. of pixels in y-direction
        INTEGER                    ::slices     ! No. of pixels in z-direction
        REAL(Kind=rk4)             ::dx        ! Pixels spacing in x-direction
        REAL(Kind=rk4)             ::dy        ! Pixels spacing in y-direction  
        REAL(Kind=rk4)             ::dz        ! Pixels spacing in z-direction
        TYPE(textent)              ::x,y,z      ! Pixel frame in x-, y- & z-dir.
!        Character(Len=256)        ::props
        REAL(kind=rk4)             ::offset_x,offset_y,offset_z! Offset of coordinates x,y,z
        endtype
      
      TYPE vector
        REAL(kind=rk8)             ::x,y,z
      ENDTYPE vector

END MODULE variable_module

