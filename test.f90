PROGRAM array
  IMPLICIT NONE
  integer :: i,j
    do i=1,9
      do j=10,15
        if(j==13) exit
        print*, i,j
      end do
    end do
END PROGRAM array
