SUBROUTINE shift
vector Y(250) 
vector Y1(250)
DO i = 72,235
 Y1(i) = Y(i-4)
ENDDO
RETURN
end
