! ---
! simulation
! ---
module tools_mod
  implicit none
  public :: stop_mc3d, wrong_input, taste, number2string
contains

  ! ################################################################################################
  ! stop: general
  ! ---
  subroutine stop_mc3d()
    print *, "mc3d stopped."
    stop
  end subroutine stop_mc3d


  ! ################################################################################################
  ! stop because of wrong input
  ! ---
  subroutine wrong_input()
    print *, "Wrong input: Chosen value outside parameter range."
    print *, "mc3d stopped."
    stop
  end subroutine wrong_input


  ! ################################################################################################
  ! input any character before continuing
  ! ---
  subroutine taste()
    implicit none
    character(len=2) :: hchar1
    ! ---
    print *,"<taste>"
    read *, hchar1
  end subroutine taste


  ! ################################################################################################
  ! convert number (000 ... 999) to string ('000' ... '999')
  ! ---
  subroutine number2string(number,string)
    integer, intent(in)           :: number
    character(len=3), intent(out) :: string

    character(len=1), dimension(0:9) :: snumber
    integer, dimension(1:3) :: counter
    integer :: counter_sum
    ! ---
    if (number<0 .or. number>999) then
       print *, "<!> error in <sr number2string>:"
       print *, "    value argument 'number' exceeds definition range (0-999)"
       print *, "    number = ", number
       stop
    end if

    snumber(0) = '0'
    snumber(1) = '1' 
    snumber(2) = '2' 
    snumber(3) = '3' 
    snumber(4) = '4' 
    snumber(5) = '5' 
    snumber(6) = '6' 
    snumber(7) = '7' 
    snumber(8) = '8' 
    snumber(9) = '9'

    counter(:) = 0

    do
       counter_sum = counter(1) + 10*counter(2) + 100*counter(3)
       if (counter_sum==number) then
          exit
       else
          ! increase counter
          counter(1) = counter(1) + 1
          if (counter(1)>9) then
             counter(1) = 0
             counter(2) = counter(2) + 1
             if (counter(2)>9) then
                counter(2) = 0
                counter(3) = counter(3) + 1
             end if
          end if
          cycle
       end if
    end do

    ! create string
    string = snumber(counter(3)) // snumber(counter(2)) // snumber(counter(1))    
  end subroutine number2string

end module tools_mod

