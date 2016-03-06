!------------------------------------------------------------------------------!
! This file is part of Mol3D.
!
!    Mol3D is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    Mol3D is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with Mol3D.  If not, see <http://www.gnu.org/licenses/>.
!
!    Diese Datei ist Teil von Mol3D.
!
!    Mol3D ist Freie Software: Sie können es unter den Bedingungen
!    der GNU General Public License, wie von der Free Software Foundation,
!    Version 3 der Lizenz oder (nach Ihrer Wahl) jeder späteren
!    veröffentlichten Version, weiterverbreiten und/oder modifizieren.
!
!    Mol3D wird in der Hoffnung, dass es nützlich sein wird, aber
!    OHNE JEDE GEWÄHRLEISTUNG, bereitgestellt; sogar ohne die implizite
!    Gewährleistung der MARKTFÄHIGKEIT oder EIGNUNG FÜR EINEN BESTIMMTEN ZWECK.
!    Siehe die GNU General Public License für weitere Details.
!
!    Sie sollten eine Kopie der GNU General Public License zusammen mit diesem
!    Programm erhalten haben. Wenn nicht, siehe <http://www.gnu.org/licenses/>.
!------------------------------------------------------------------------------!
! This module defines the Random number generator. Define your prefered
! one or use one offered here. 
!
MODULE randgen_type
    ! This is the generic type and calls the actual random number generator
    USE datatype
    USE common_type, &
       GetType_common => GetType, GetName_common => GetName, &
       Initialized_common => Initialized
    IMPLICIT NONE
    !--------------------------------------------------------------------------!
    PRIVATE
    ! 
    !--------------------------------------------------------------------------!
    ! Period parameters
!~     INTEGER, PARAMETER :: n = 624, n1 = n+1, m = 397, mata = -1727483681
    !                                    constant vector a
!~     INTEGER, PARAMETER :: umask = -2147483647 -1
    !                                    most significant w-r bits
!~     INTEGER, PARAMETER :: lmask =  2147483647
    !                                    least significant r bits
    ! Tempering parameters
!~     INTEGER, PARAMETER  :: tmaskb= -1658038656, tmaskc= -272236544
    REAL(r2), PARAMETER :: norm_r2 = 2.328306e-10
    TYPE Randgen_TYP
        TYPE(Common_TYP)                   :: gentype
        REAL(kind=r2)                      :: rndx
        ! generator array for the state vector
        INTEGER, DIMENSION(1:200) :: state_array
        
!~         INTEGER,
!~         !                     the array for the state vector
!~         INTEGER                        :: mt(0:n-1), mti = n1
!~         !                     mti==N+1 means mt[N] is not initialized
    END TYPE Randgen_TYP
    SAVE

    !--------------------------------------------------------------------------!
    PUBLIC :: &
        ! types
        Randgen_TYP, &
        ! methods
        InitRandgen, &
        CloseRandgen, &
        GetGenID, &
        GetGenName, &
        RandgenInitialized, &
        GetNewRandomNumber
    !--------------------------------------------------------------------------!
CONTAINS

    SUBROUTINE InitRandgen(this, seed, gen_name)
        IMPLICIT NONE
        !----------------------------------------------------------------------!
        TYPE(Randgen_TYP)                    :: this
        INTEGER                              :: seed
        INTEGER                              :: N, i
        CHARACTER(LEN=*)                     :: gen_name
        !----------------------------------------------------------------------!
        INTENT(IN)                           :: seed, gen_name
        INTENT(INOUT)                        :: this

        !----------------------------------------------------------------------!
        IF (RandgenInitialized(this)) THEN
            PRINT *, 'ERROR, Generator already initialized!'
            STOP 
        END IF
        
        !  gen_name = name of the desired random number generator,
        !             possible values at the moment:
        ! 'MT' , 'TAUS88', 'COMPILER', 'KISS99'


        SELECT CASE(gen_name)
!~             CASE('MT')
!~                 CALL init_MT(this, seed*5)
!~                 i = 1
            CASE('TAUS88')
                N = 3
                i = 2
            CASE('LFSR113')
                N = 4
                i = 3
            CASE('KISS99')
                N = 4
                i = 4
            CASE('CONG')
                N = 1
                i = 5
            CASE('COMPILER')
            !    CALL RANDOM_SEED(N)
                i = 6
            CASE DEFAULT
                PRINT *, 'ERROR, chosen random number generator not found.'
                STOP
        END SELECT
        N = 6
        CALL InitCommon(this%gentype, i, gen_name)
        
        ! Mol3D uses CONG to generate the seeds, CONG is initialized with a
        ! seed based on the actual thread, thus Mol3D calculates a few numbers
        ! with CONG to ensure the quality of the numbers
        this%state_array(1) = seed
        DO i = 1, 2000
            this%state_array(1) = CONG(this%state_array(1))
        END DO
        ! Now fill the state array

        DO i = 2, 200
            this%state_array(i) = CONG(this%state_array(i-1))
        END DO

        IF  (GetGenName(this) == 'COMPILER') THEN
            CALL RANDOM_SEED(N)
            CALL RANDOM_SEED(PUT=this%state_array(1:N))
        END IF
    END SUBROUTINE InitRandgen

    SUBROUTINE CloseRandgen(this)
        IMPLICIT NONE
        !----------------------------------------------------------------------!
        TYPE(Randgen_TYP), INTENT(INOUT) :: this
        !----------------------------------------------------------------------!
        CALL CloseCommon(this%gentype)
    END SUBROUTINE CloseRandgen

    SUBROUTINE GetNewRandomNumber(this, randno)
        IMPLICIT NONE
        !----------------------------------------------------------------------!
        TYPE(Randgen_TYP), INTENT(INOUT) :: this
        REAL(kind=r2),INTENT(OUT)        :: randno
        !----------------------------------------------------------------------!
        SELECT CASE(GetGenID(this))
            !CASE(1)
            !    randno = getNumber_MT(this)
            CASE(2)
                randno = getNumber_TAUS88(this)
            CASE(3)
                randno = getNumber_lfsr113(this)
            CASE(4)
                randno = getNumber_KISS99(this)
            CASE(5)
                randno = getNumber_CONG(this)
            CASE(6)
                CALL RANDOM_NUMBER(randno)
        CASE DEFAULT
            PRINT *, 'ERROR, chosen random number generator not initialized.'
            STOP
        END SELECT

    END SUBROUTINE GetNewRandomNumber

    PURE FUNCTION GetGenID(this) RESULT(ut)
        IMPLICIT NONE
        !----------------------------------------------------------------------!
        TYPE(Randgen_TYP), INTENT(IN) :: this
        INTEGER :: ut
        !----------------------------------------------------------------------!
        ut = GetType_common(this%gentype)
    END FUNCTION GetGenID


    PURE FUNCTION GetGenName(this) RESULT(un)
        IMPLICIT NONE
        !----------------------------------------------------------------------!
        TYPE(Randgen_TYP), INTENT(IN) :: this
        CHARACTER(LEN=32) :: un
        !----------------------------------------------------------------------!
        un = GetName_common(this%gentype)
    END FUNCTION GetGenName

    PURE FUNCTION RandgenInitialized(this) RESULT(i)
        IMPLICIT NONE
        !----------------------------------------------------------------------!
        TYPE(Randgen_TYP), INTENT(IN) :: this
        LOGICAL :: i

        i = Initialized_common(this%gentype)
    END FUNCTION RandgenInitialized

    !--------------------------------------------------------------------------!
    !------------ start of the individual random number generator -------------!
    !--------------------------------------------------------------------------!
    !
    ! CONG - very simple linear congruential generator
    PURE FUNCTION CONG(s1)
        IMPLICIT NONE
        !----------------------------------------------------------------------!
        INTEGER                       :: CONG
        INTEGER, INTENT(IN)           :: s1
        !----------------------------------------------------------------------!
        CONG = 69069 * s1 + 1327217885
    END FUNCTION CONG
    
    FUNCTION getNumber_CONG(this) result(random_numb)
        IMPLICIT NONE
        !----------------------------------------------------------------------!
        TYPE(Randgen_TYP)                    :: this
        REAL(r2)                             :: random_numb
        !----------------------------------------------------------------------!
        this%state_array(1) = CONG(this%state_array(1))
        random_numb = this%state_array(1) * norm_r2 + 0.5_r2
    END FUNCTION getNumber_CONG
    
    !--------------------------------------------------------------------------!

    FUNCTION getNumber_KISS99(this)  RESULT(random_numb)
        IMPLICIT NONE
        !----------------------------------------------------------------------!
        TYPE(Randgen_TYP)                    :: this
        REAL(r2)                             :: random_numb
        INTEGER                              :: t
        !----------------------------------------------------------------------!
        ! equals CONG generator
        this%state_array(1) = CONG(this%state_array(1))

        ! 3 Xorshift generators
        this%state_array(2) = IEOR(this%state_array(2), ISHFT(this%state_array(2), 13))
        this%state_array(2) = IEOR(this%state_array(2), ISHFT(this%state_array(2),-17))
        this%state_array(2) = IEOR(this%state_array(2), ISHFT(this%state_array(2), 5))

        ! 2 multiply-with-carry generators
        this%state_array(3) = 30903*IAND(this%state_array(3), 65535)+ISHFT(this%state_array(3), -16)
        this%state_array(4) = 18000*IAND(this%state_array(4), 65535)+ISHFT(this%state_array(4), -16)
        t = (ISHFT(this%state_array(4), 16) + this%state_array(3))

        ! combination
        random_numb = (IEOR(t, this%state_array(1)) + this%state_array(2)) * norm_r2 + 0.5_r2

    END FUNCTION getNumber_KISS99

    FUNCTION getNumber_MRG32k3a(this) RESULT(random_numb)
        ! Generates a random number between 0 and 1.
        ! Reference: 'Good Parameters and implementations for
        !            CMR random number generators''
        ! L'Ecuyer, P. (1998)
        ! Not sure if it's implemented correctly -> not used at the moment
        IMPLICIT NONE
        !----------------------------------------------------------------------!
        TYPE(Randgen_TYP)                    :: this
        REAL(kind=r2)                        :: random_numb
        REAL(kind=r2)                        :: p1, p2
        INTEGER                              :: k
        !----------------------------------------------------------------------!
        REAL(kind=r2), PARAMETER                   :: m1   = 4294967087.0
        REAL(kind=r2), PARAMETER                   :: m2   = 4294944443.0
        REAL(kind=r2), PARAMETER                   :: a12  =    1403580.0
        REAL(kind=r2), PARAMETER                   :: a13n =     810728.0
        REAL(kind=r2), PARAMETER                   :: a21  =     527612.0
        REAL(kind=r2), PARAMETER                   :: a23n =    1370589.0
        !----------------------------------------------------------------------!
        !
        ! Component 1
        p1 = a12 * this%state_array(2) - a13n * this%state_array(1)
        k = p1 / m1
        p1 = p1 - k * m1
        IF (p1 < 0.0) p1 = p1 + m1
        this%state_array(1) = this%state_array(2)
        this%state_array(2) = this%state_array(3)
        this%state_array(3) = p1

        ! Component 2
        p2 = a21 * this%state_array(6) - a23n * this%state_array(4)
        k = p2 / m2
        p2 = p2 - k * m2
        IF (p2 < 0.0) p2 = p2 + m2
        this%state_array(4) = this%state_array(5)
        this%state_array(5) = this%state_array(6)
        this%state_array(6) = p2

        ! Combination

        IF (p1 <= p2) THEN
            random_numb = ((p1 - p2 + m1) * norm_r2)
        ELSE
            random_numb = ((p1 - p2) * norm_r2)
        END IF

    END FUNCTION getNumber_MRG32k3a

    FUNCTION getNumber_TAUS88(this) RESULT(random_numb)
        ! Generates a random number between 0 and 1.
        ! Translated from C function in:
        ! Reference:
        ! L'Ecuyer, P. (1996) `Maximally equidistributed combined Tausworthe
        ! generators', Math. of Comput., 65, 203-213.

        ! The cycle length is claimed to be about 2^(88) or about 3E+26.
        ! Actually - (2^31 - 1).(2^29 - 1).(2^28 - 1).

        IMPLICIT NONE
        !----------------------------------------------------------------------!
        TYPE(Randgen_TYP), INTENT(INOUT) :: this
        INTEGER                          :: b
        REAL(kind=r2)                    :: random_numb
        !----------------------------------------------------------------------!

        ! N.B. ISHFT(i,j) is a bitwise (non-circular) shift operation;
        !      to the left if j > 0, otherwise to the right.
        
        b  = ISHFT( IEOR( ISHFT(this%state_array(1), 13), this%state_array(1)), -19)
        this%state_array(1) = IEOR( ISHFT( IAND(this%state_array(1),-2), 12), b)
        b  = ISHFT( IEOR( ISHFT(this%state_array(2), 2), this%state_array(2)), -25)
        this%state_array(2) = IEOR( ISHFT( IAND(this%state_array(2),-8), 4), b)
        b  = ISHFT( IEOR( ISHFT(this%state_array(3), 3), this%state_array(3)), -11)
        this%state_array(3) = IEOR( ISHFT( IAND(this%state_array(3),-16), 17), b)
        
        random_numb = IEOR( IEOR(this%state_array(1),this%state_array(2)), this%state_array(3)) *   &
                        norm_r2+ 0.5_r2
        this%rndx = random_numb
    END FUNCTION getNumber_TAUS88

    !--------------------------------------------------------------------------!
    !
    ! Lin_Feedback_Shift_Reg
    !
    ! L'Ecuyer's 1999 random number generator.
    ! Fortran version by Alan.Miller @ vic.cmis.csiro.au
    ! N.B. This version is compatible with Lahey's ELF90
    ! This version requires that the default integer type is of 32-bits
    ! http://www.ozemail.com.au/~milleraj
    ! http://users.bigpond.net.au/amiller/
    ! Latest revision - 12 January 2001

    FUNCTION getNumber_lfsr113(this) RESULT(random_numb)
        ! Generates a random number between 0 and 1.
        ! Translated from C function
        ! Reference:
        ! L'Ecuyer, P. (1999) `Tables of maximally equidistributed combined LFSR
        ! generators', Math. of Comput., 68, 261-269.

        ! The cycle length is claimed to be about 2^(113) or about 10^(34).
        ! Actually - (2^31 - 1).(2^29 - 1).(2^28 - 1).(2^25 - 1)
        IMPLICIT NONE
        !----------------------------------------------------------------------!
        TYPE(Randgen_TYP), INTENT(INOUT) :: this
        REAL(r2)  :: random_numb
        INTEGER   :: b
        !----------------------------------------------------------------------!
        ! N.B. ISHFT(i,j) is a bitwise (non-circular) shift operation;
        !      to the left if j > 0, otherwise to the right.

        b  = ISHFT( IEOR( ISHFT(this%state_array(1), 6), this%state_array(1)), -13)
        this%state_array(1) = IEOR( ISHFT( IAND(this%state_array(1),-2), 18), b)
        b  = ISHFT( IEOR( ISHFT(this%state_array(2), 2), this%state_array(2)), -27)
        this%state_array(2) = IEOR( ISHFT( IAND(this%state_array(2),-8), 2), b)
        b  = ISHFT( IEOR( ISHFT(this%state_array(3), 13), this%state_array(3)), -21)
        this%state_array(3) = IEOR( ISHFT( IAND(this%state_array(3),-16), 7), b)
        b  = ISHFT( IEOR( ISHFT(this%state_array(4), 3), this%state_array(4)), -12)
        this%state_array(4) = IEOR( ISHFT( IAND(this%state_array(4),-128), 13), b)

        ! The constant norm_r2 below is the reciprocal of (2^32 - 1)
        random_numb = IEOR( IEOR( IEOR(this%state_array(1), this%state_array(2)), &
                           this%state_array(3)), this%state_array(4)) *        &
                           norm_r2 + 0.5_r2

    END FUNCTION getNumber_lfsr113

    !--------------------------------------------------------------------------!
    !
    ! Mersenne Twister (MT)
    
    !***********************************************************************
    ! A Fortran-program for MT19937: Real number version
    !***********************************************************************
    ! Code converted using TO_F90 by Alan Miller
    ! Date: 1999-11-26  Time: 17:09:23
    ! Latest revision - 5 February 2002
    ! A new seed initialization routine has been added based upon the new
    ! C version dated 26 January 2002.
    ! This version assumes that integer overflows do NOT cause crashes.
    ! This version is compatible with Lahey's ELF90 compiler,
    ! and should be compatible with most full Fortran 90 or 95 compilers.
    ! Notice the strange way in which umask is specified for ELF90.
     
    ! genrand() generates one pseudorandom real number (double) which is
    ! uniformly distributed on [0,1]-interval, for each call.
    ! sgenrand(seed) set initial values to the working area of 624 words.
    ! Before genrand(), sgenrand(seed) must be called once.  (seed is any 32-bit
    ! integer except for 0).
    ! Integer generator is obtained by modifying two lines.
    !   Coded by Takuji Nishimura, considering the suggestions by
    ! Topher Cooper and Marc Rieffel in July-Aug. 1997.

    ! This library is free software; you can redistribute it and/or modify it
    ! under the terms of the GNU Library General Public License as published by
    ! the Free Software Foundation; either version 2 of the License, or (at your
    ! option) any later version.   This library is distributed in the hope that
    ! it will be useful, but WITHOUT ANY WARRANTY; without even the implied
    ! warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
    ! See the GNU Library General Public License for more details.
    ! You should have received a copy of the GNU Library General Public License
    ! along with this library; if not, write to the Free Foundation, Inc.,
    ! 59 Temple Place, Suite 330, Boston, MA 02111-1307  USA

    ! Copyright (C) 1997 Makoto Matsumoto and Takuji Nishimura.
    ! When you use this, send an email to: matumoto@math.keio.ac.jp
    ! with an appropriate reference to your work.

    !***********************************************************************
    ! Fortran translation by Hiroshi Takano.  Jan. 13, 1999.

    !   genrand()      -> double precision function grnd()
    !   sgenrand(seed) -> subroutine sgrnd(seed)
    !                     integer seed

    ! This program uses the following standard intrinsics.
    !   ishft(i,n): If n > 0, shifts bits in i by n positions to left.
    !               If n < 0, shifts bits in i by n positions to right.
    !   iand (i,j): Performs logical AND on corresponding bits of i and j.
    !   ior  (i,j): Performs inclusive OR on corresponding bits of i and j.
    !   ieor (i,j): Performs exclusive OR on corresponding bits of i and j.

    !***********************************************************************
!~     SUBROUTINE init_MT(this, seed)
!~         !----------------------------------------------------------------------!
!~         TYPE(Randgen_TYP), INTENT(INOUT) :: this
!~         INTEGER, INTENT(IN)              :: seed
!~         INTEGER                          :: latest
!~         INTEGER                          :: mti
!~         !----------------------------------------------------------------------!
!~         ! This initialization is based upon the multiplier given on p.106 of the
!~         ! 3rd edition of Knuth, The Art of Computer Programming Vol. 2.
!~         ! This version assumes that integer overflow does NOT cause a crash.

!~         this%mt(0) = seed
!~         latest = seed
!~         DO mti = 1, n-1
!~             this%mti = mti
!~             latest = IEOR( latest, ISHFT( latest, -30 ) )
!~             latest = latest * 1812433253 + this%mti
!~             this%mt(this%mti) = latest
!~         END DO

!~     END SUBROUTINE init_MT

!~     FUNCTION getNumber_MT(this) RESULT(fn_val)
!~         !----------------------------------------------------------------------!
!~         TYPE(Randgen_TYP), INTENT(INOUT) :: this
!~         INTEGER                          :: mag01(0:1) = (/ 0, mata /)
!~         INTEGER                          :: kk, y
!~         REAL (kind=r2)                   :: fn_val
!~         !----------------------------------------------------------------------!

!~         ! These statement functions have been replaced with separate functions
!~         ! tshftu(y) = ISHFT(y,-11)
!~         ! tshfts(y) = ISHFT(y,7)
!~         ! tshftt(y) = ISHFT(y,15)
!~         ! tshftl(y) = ISHFT(y,-18)

!~         IF(this%mti >= n) THEN
!~         !                       generate N words at one time
!~             IF(this%mti == n+1) THEN
!~         !                            if sgrnd() has not been called,
!~                 CALL init_MT(this, 4357)
!~         !                              a default initial seed is used
!~             END IF
          
!~             DO  kk = 0, n-m-1
!~                 y = IOR(IAND(this%mt(kk),umask), IAND(this%mt(kk+1),lmask))
!~                 this%mt(kk) = IEOR(IEOR(this%mt(kk+m), ISHFT(y,-1)),mag01(IAND(y,1)))
!~             END DO
!~             DO  kk = n-m, n-2
!~                 y = IOR(IAND(this%mt(kk),umask), IAND(this%mt(kk+1),lmask))
!~                 this%mt(kk) = IEOR(IEOR(this%mt(kk+(m-n)), ISHFT(y,-1)),mag01(IAND(y,1)))
!~             END DO
!~             y = IOR(IAND(this%mt(n-1),umask), IAND(this%mt(0),lmask))
!~             this%mt(n-1) = IEOR(IEOR(this%mt(m-1), ISHFT(y,-1)),mag01(IAND(y,1)))
!~             this%mti = 0
!~         END IF

!~         y = this%mt(this%mti)
!~         this%mti = this%mti + 1
!~         y = IEOR(y, tshftu(y))
!~         y = IEOR(y, IAND(tshfts(y),tmaskb))
!~         y = IEOR(y, IAND(tshftt(y),tmaskc))
!~         y = IEOR(y, tshftl(y))

!~         IF(y < 0) THEN
!~           fn_val = (DBLE(y) + 2.0D0**32) / (2.0D0**32 - 1.0D0)
!~         ELSE
!~           fn_val = DBLE(y) / (2.0D0**32 - 1.0D0)
!~         END IF

!~         this%rndx = fn_val
!~     END FUNCTION getNumber_MT

!~     FUNCTION tshftu(y) RESULT(fn_val)
!~         INTEGER, INTENT(IN) :: y
!~         INTEGER             :: fn_val

!~         fn_val = ISHFT(y,-11)
!~         RETURN
!~     END FUNCTION tshftu

!~     FUNCTION tshfts(y) RESULT(fn_val)
!~         INTEGER, INTENT(IN) :: y
!~         INTEGER             :: fn_val

!~         fn_val = ISHFT(y,7)
!~         RETURN
!~     END FUNCTION tshfts


!~     FUNCTION tshftt(y) RESULT(fn_val)
!~         INTEGER, INTENT(IN) :: y
!~         INTEGER             :: fn_val

!~         fn_val = ISHFT(y,15)
!~         RETURN
!~     END FUNCTION tshftt

!~     FUNCTION tshftl(y) RESULT(fn_val)
!~         INTEGER, INTENT(IN) :: y
!~         INTEGER             :: fn_val

!~         fn_val = ISHFT(y,-18)
!~         RETURN
!~     END FUNCTION tshftl

END MODULE randgen_type

!~ PROGRAM main
! A simple program to test the random number generators
!
!~     USE randgen_type
!~     USE datatype
!~     IMPLICIT NONE

!~     INTEGER, PARAMETER :: buckets = 13
!~     INTEGER            :: seed
!~     INTEGER(8)         :: j, i
!~     INTEGER(8)         :: count, no, error_counter
!~     INTEGER(8), DIMENSION(buckets) :: bucket
!~     REAL (r2)          :: temp, big, small, average, sumsq, stdev, chi
!~     REAL               :: t1, t0
!~     TYPE(Randgen_TYP)  :: rand_nr


!~     no = 2000000000
!~     bucket(:) = 0
!~     DO i = 4, 4
!~         temp = 0.0_r2
!~         big = 0.5_r2
!~         small = 0.5_r2
!~         count = 0
!~         average = 0.0_r2
!~         sumsq = 0.0_r2
!~         bucket(:) = 0
!~         error_counter = 0
!~         seed = 1 ! 2 , 3 ,4 ...
        
!~         IF (i == 1) THEN
!~             PRINT *, 'KISS99'
!~             CALL InitRandgen(rand_nr, seed, 'KISS99')
!~         ELSEIF (i == 2) THEN
!~             PRINT *, 'TAUS88'
!~             CALL InitRandgen(rand_nr, seed, 'TAUS88')
!~         ELSEIF (i == 3) THEN
!~             PRINT *, 'LFSR113'
!~             CALL InitRandgen(rand_nr, seed, 'LFSR113')
!~         ELSEIF (i == 4) THEN
!~             PRINT *, 'CONG'
!~             CALL InitRandgen(rand_nr, seed, 'CONG')
!~         END IF
!~         CALL cpu_time(t0)
!~         DO  j=1, no
!~             CALL GetNewRandomNumber(rand_nr, temp)
!~             !IF (temp < 1.0e-10_r2 .or. temp == 1.0_r2- 1.0e-10_r2) THEN
!~             !    error_counter = error_counter + 1
!~             !END IF
!~             IF (temp > big) THEN
!~                 big = temp
!~             ELSE IF (temp < small) THEN
!~                 small = temp
!~             END IF
!~             CALL update(temp, count, average, sumsq, bucket, buckets)
!~         END DO
!~         CALL cpu_time(t1)
!~         stdev = SQRT( sumsq / (count - 1) )
!~         chi = 0.0
!~         print *, ''
!~         DO j = 1, buckets
!~             chi = chi + (bucket(j)-REAL(no, kind=r2)/buckets)**2/(REAL(no, kind=r2)/buckets)
!~         END DO
!~         WRITE(*, *) ' Smallest = ', small, '  Largest = ', big
!~         WRITE(*, *) ' Average = ', average, '  Std. devn. = ', stdev
!~         WRITE(*, *) ' Std. devn. should be about 1/sqrt(12) = 0.288675'
!~         WRITE(*, *) ' Chi^2 = ', chi
!~         !WRITE(*, *) ' Error values found: ', error_counter
!~         WRITE(*, *) ' Calculations took = ', t1-t0 ,'s'
!~         WRITE(*, *) ''
!~         CALL CloseRandgen(rand_nr)
!~     END DO


!~     CONTAINS

!~         SUBROUTINE update(x, n, avge, sumsq, bucket, bn)
!~         REAL (r2), INTENT(IN)      :: x
!~         INTEGER(8), INTENT(IN OUT)    :: n
!~         INTEGER, INTENT(IN)    :: bn
!~         REAL (r2), INTENT(IN OUT)  :: avge, sumsq
!~         INTEGER(8), DIMENSION(bn), INTENT(IN OUT) :: bucket

!~         REAL (r2)  :: dev

!~         n = n + 1
!~         dev = x - avge
!~         avge = avge + dev / n
!~         sumsq = sumsq + dev*(x - avge)
!~         bucket(int(x*bn)+ 1) = bucket(int(x*bn) + 1) + 1

!~         END SUBROUTINE update
        
!~         FUNCTION RAN0(IDUM)
!~               IMPLICIT NONE
!~               INTEGER, INTENT(IN OUT) :: IDUM
!~               INTEGER, PARAMETER :: IA=16807, IM=2147483647, IQ=127773, IR=2836, MASK=123459876
!~               REAL(r2), PARAMETER :: AM=1.0_r2/IM
!~               REAL(r2) :: RAN0
!~               INTEGER  :: k

!~               k = IDUM/IQ
!~               IDUM=IA*(IDUM-k*IQ)-IR*k
!~               IF (idum .lt. 0) idum=idum+IM
!~               RAN0 = AM*IDUM
              

!~         END FUNCTION RAN0

!~ END PROGRAM main
