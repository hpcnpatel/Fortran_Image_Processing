module constants_module
      integer,parameter  :: rk=8
      integer,parameter  :: dk=8  ! not less than rk
      ! PI from http://www.pibel.de/
real(kind=dk),parameter    :: PI_long=&
& 3.1415926535897932384626433832795028841971693993751058209749445923078164062862089986280348253421170679821480865132_dk
end module constants_module
