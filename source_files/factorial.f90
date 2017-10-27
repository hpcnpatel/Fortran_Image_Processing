Module factorial
contains

Function fact(n)
IMPLICIT NONE

integer ::i,n
REAL::fact
fact=1
do i = 1, n
fact = fact * i
end do
End function fact 
END module factorial
