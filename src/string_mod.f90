module string_mod
  implicit none
  public :: read_int, read_real, check_int, check_real, real2strg, strg2real,  &
            int2strg, strg2int, num2strg4, manual

contains

  ! ################################################################################################
  ! read integer - allow for "?" => call manual entry
  ! ---
  subroutine read_int(int,manpage)
    implicit none
    integer, intent(out) :: int
    integer, intent(in)  :: manpage
    character(len=100)   :: hch_in
    ! ---
    print *,"<'?' - explanation>"
    do
       read *, hch_in
       if ( check_int(hch_in(1:len_trim(hch_in))) ) then
          call strg2int(hch_in,int)
          exit
       else
          if (hch_in(1:len_trim(hch_in)) == "?")  then
          call manual(manpage)
             cycle
          else if (hch_in(1:len_trim(hch_in)) == "x")  then
             ! clean up & stop the mc3d
             print *, "mc3d stopped."
             stop
          else
             print *,"Wrong input: Type mismatch. Try again."
             cycle
          end if
       end if
    end do
  end subroutine read_int


  ! ################################################################################################
  ! read real - allow for "?" => call manual entry
  ! ---
  subroutine read_real(real_x,manpage)
    use datatype

    real(kind=r2), intent(out) :: real_x
    integer, intent(in)  :: manpage
    character(len=100)   :: hch_in
    ! ---
    print *,"<'?' - explanation>"
    do
       read *, hch_in
       if ( check_real(hch_in(1:len_trim(hch_in))) ) then
          call strg2real(hch_in,real_x)
          exit
       else
          if (hch_in(1:len_trim(hch_in)) == "?")  then
          call manual(manpage)
             cycle
          else if (hch_in(1:len_trim(hch_in)) == "x")  then
             ! clean up & stop the mc3d
             print *, "mc3d stopped."
             stop
          else
             print *,"Wrong input: Type mismatch. Try again."
             cycle
          end if
       end if
    end do
  end subroutine read_real


  ! ################################################################################################
  ! check: - ist der string "strg" der laenge "len_strg" eine integer-zahl oder nicht?
  ! ---
  function check_int(strg) result(check_int_result)

    character(len=*), intent(in) :: strg
    logical :: check_int_result
    
    logical :: rote_karte
    integer :: ptr_ss, ichar, ichar1
    character(len=1), dimension(1:10) :: ss
    ! ---
    check_int_result = .true.

    ss(01) = "0"
    ss(02) = "1"
    ss(03) = "2"
    ss(04) = "3"
    ss(05) = "4"
    ss(06) = "5"
    ss(07) = "6"
    ss(08) = "7"
    ss(09) = "8"
    ss(10) = "9"

    ichar1 = 1

    ! vorzeichen getrennt behandeln (da es nicht zwingend vorhanden sein muss)
    if (strg(ichar1:ichar1)=="+" .or. strg(1:1)=="-") then
       ichar1 = ichar1 + 1
    end if

    rote_karte = .false.
    do ichar=ichar1,len_trim(strg)       
       ptr_ss = 1
       if (.not. rote_karte) then
          do
             if ( strg(ichar:ichar) == ss(ptr_ss) ) then
                exit
             else
                ptr_ss = ptr_ss + 1
                if (ptr_ss>10) then
                   rote_karte = .true.
                   exit
                end if
             end if
          end do
       end if
    end do

    if (rote_karte) then
       check_int_result = .false.
    end if
  end function check_int


  ! ---
  ! check: - ist der string "strg" der laenge "len_strg" eine (fliesskomma-)zahl oder nicht?
  !                erlaubt: -1.23e+23, aber auch 1e+12, 1.2, 1234, .1123, .e-23(=0)
  !          nicht erlaubt: 1e21
  ! ---
  function check_real(strg) result(check_real_result)

    character(len=*), intent(in) :: strg
    logical :: check_real_result
    
    logical :: gelbe_karte, rote_karte, found
    integer :: ptr_step, ptr_ss, ichar, max_ptr_step, ichar1
    integer, dimension(1:4,1:2) :: bss
    logical, dimension(1:17) :: lss, lssx
    character(len=1), dimension(1:17) :: ss
    ! ---
    check_real_result = .true.

    ss(01) = "0"; lss(01) = .false.
    ss(02) = "1"; lss(02) = .false.
    ss(03) = "2"; lss(03) = .false.
    ss(04) = "3"; lss(04) = .false.
    ss(05) = "4"; lss(05) = .false.
    ss(06) = "5"; lss(06) = .false.
    ss(07) = "6"; lss(07) = .false.
    ss(08) = "7"; lss(08) = .false.
    ss(09) = "8"; lss(09) = .false.
    ss(10) = "9"; lss(10) = .false.
    ss(11) = "."; lss(11) = .true. 
    ss(12) = "e"; lss(12) = .true. 
    ss(13) = "d"; lss(13) = .true. 
    ss(14) = "E"; lss(14) = .true. 
    ss(15) = "D"; lss(15) = .true. 
    ss(16) = "+"; lss(16) = .true. 
    ss(17) = "-"; lss(17) = .true. 

    lssx(:) = .false.

    bss(1,1) =  1
    bss(1,2) = 11

    bss(2,1) = 12
    bss(2,2) = 15

    bss(3,1) = 16
    bss(3,2) = 17

    bss(4,1) =  1
    bss(4,2) = 10

    max_ptr_step = 4

    ptr_step    = 1
    gelbe_karte = .false.
    rote_karte  = .false.
    
    ichar1 = 1

    ! vorzeichen getrennt behandeln (da es nicht zwingend vorhanden sein muss)
    if (strg(ichar1:ichar1)=="+" .or. strg(1:1)=="-") then
       ichar1 = ichar1 + 1
    end if

    do ichar=ichar1,len_trim(strg)

       if (.not. rote_karte) then
          ptr_ss = bss(ptr_step,1)
          do
             do
                if ( strg(ichar:ichar) == ss(ptr_ss) ) then
                   found       = .true.
                   gelbe_karte = .false.

                   ! war dieses zeichen bereits vorhanden und darf maximal einmal vorkommen?
                   if (lss(ptr_ss)) then
                      if (lssx(ptr_ss)) then
                         rote_karte = .true.
                      else
                         lssx(ptr_ss) = .true.
                      end if
                   end if

                   exit
                else
                   ptr_ss = ptr_ss + 1
                   if (ptr_ss>bss(ptr_step,2)) then
                      found = .false.
                      exit
                   end if
                end if
             end do

             if (found) then
                exit
             else
                ! character konnte im aktuellen set erlaubter character nicht gefunden werden
                ! -> check: im naechsten vorhanden?
                if (gelbe_karte .or. ptr_step==max_ptr_step) then
                   rote_karte  = .true.
                   exit
                else
                   ptr_step    = ptr_step + 1
                   gelbe_karte = .true.
                   ptr_ss      = bss(ptr_step,1)                   
                end if
             end if
          end do
       end if
    end do

    if (rote_karte) then
       check_real_result = .false.
    end if
  end function check_real


  ! ################################################################################################
  ! convert real to string
  ! ---
  subroutine real2strg(rnumber,string)
   ! use var_global

    real, intent(in)                :: rnumber
    character(len=100), intent(out) :: string
    ! ---
    ! new version by J. Sauter:
    ! Do the type casting without reference to real files
    ! in order to avoid conflicts between multiple simultaneous running
    ! instances of mc3d
    write(Unit=string,FMT=*) rnumber

  end subroutine real2strg


  ! ################################################################################################
  ! convert string to real
  ! ---
  subroutine strg2real(string,rnumber)
    use datatype
    !use var_global

    character(len=100), intent(in) :: string
    real(kind=r2), intent(out)     :: rnumber
    ! ---
    ! new version by J. Sauter:
    ! Do the type casting without reference to real files
    ! in order to avoid conflicts between multiple simultaneous running
    ! instances of mc3d
    Read(Unit=string(1:len_trim(string)),FMT=*) rnumber

  end subroutine strg2real

  
  ! ################################################################################################
  ! convert integer to string
  ! ---
  subroutine int2strg(inumber,string)

    integer, intent(in)              :: inumber
    character(len=100), intent(out) :: string
    ! ---
    ! new version by J. Sauter:
    ! Do the type casting without reference to real files
    ! in order to avoid conflicts between multiple simultaneous running
    ! instances of mc3d
    Read(Unit=inumber, FMT=*) string
    
  end subroutine int2strg


  ! ################################################################################################
  ! convert string to integer
  ! ---
  subroutine strg2int(string,inumber)
    !use var_global

    character(len=100), intent(in) :: string
    integer, intent(out)           :: inumber
    ! ---
    ! new version by J. Sauter:
    ! Do the type casting without reference to real files
    ! in order to avoid conflicts between multiple simultaneous running
    ! instances of mc3d
    Read(Unit=string(1:len_trim(string)),FMT=*) inumber

  end subroutine strg2int


  ! ################################################################################################
  ! convert number (0000 ... 9999) to string ('0000' ... '9999')
  ! ---
  subroutine num2strg4(number,string)
    integer, intent(in)           :: number
    character(len=4), intent(out) :: string

    character(len=1), dimension(0:9) :: snumber
    integer, dimension(1:4) :: counter
    integer :: counter_sum
    ! ---
    if (number<0 .or. number>9999) then
       print *, "<!> error in <sr number2string4>:"
       print *, "    value argument 'number' exceeds definition range (0-9999)"
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
       counter_sum = counter(1) + 10*counter(2) + 100*counter(3) + 1000*counter(4)
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
                if (counter(3)>9) then
                   counter(3) = 0
                   counter(4) = counter(4) + 1
                end if
             end if
          end if
          cycle
       end if
    end do

    ! create string
    string = snumber(counter(4)) // snumber(counter(3)) // snumber(counter(2)) // snumber(counter(1))    
  end subroutine num2strg4


  ! ################################################################################################
  ! read manual entry
  ! ---
  subroutine manual(manpage)
    use var_global
    
    integer, intent(in) :: manpage
    logical :: hl1, hl2
    character(len=4) :: hch1
    character(len=100) :: hch2
    ! ---
    call num2strg4(manpage,hch1)
    
    open(unit=1,file=path_misc//"manpages.v4.txt",action="read",status="unknown",form="formatted")

    hl1 = .false. ! hl1=.true.: entry found
    hl2 = .false. ! hl2=.true.: reading of entry in progress
    do
       read(unit=1,fmt=*) hch2
       if (hch2(1:8) == "****"//hch1) then
          hl1 = .true.
          exit
       else
          if (hch2(1:7) == "****end") then
             exit
          else
             cycle
          end if
       end if
    end do
    
    if (hl1) then
       print *,"---"
       do
          read(unit=1,fmt=*) hch2
          if (hch2(1:4) /= "****") then
             print *, "> ", hch2(1:len_trim(hch2))
             cycle
          else
             exit
          end if
       end do
       print *,"---"
    else
       print *, "manual entry not available"
    end if
    close(unit=1)
    print *,">> your choice:"

  end subroutine manual
  
end module string_mod
