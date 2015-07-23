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
    INTEGER, PARAMETER       :: MNew=714025, IA=1366, IC=150889
    ! Period parameters
    INTEGER, PARAMETER :: n = 624, n1 = n+1, m = 397, mata = -1727483681
    !                                    constant vector a
    INTEGER, PARAMETER :: umask = -2147483647 -1
    !                                    most significant w-r bits
    INTEGER, PARAMETER :: lmask =  2147483647
    !                                    least significant r bits
    ! Tempering parameters
    INTEGER, PARAMETER :: tmaskb= -1658038656, tmaskc= -272236544

    
    TYPE Randgen_TYP
        TYPE(Common_TYP)               :: gentype
        REAL(kind=r2)                  :: rndx
        INTEGER                        :: s1, s2, s3
        !                     the array for the state vector
        INTEGER                        :: mt(0:n-1), mti = n1
        !                     mti==N+1 means mt[N] is not initialized
    END TYPE Randgen_TYP
    SAVE

    !--------------------------------------------------------------------------!
    PUBLIC :: &
        ! types
        Randgen_TYP, &
        ! methods
        InitRandgen, &
        CloseRandgen, &
        GetGenSeed, &
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
        INTEGER, DIMENSION(:), ALLOCATABLE   :: seed_array
        CHARACTER(LEN=*)                     :: gen_name
        !----------------------------------------------------------------------!
        INTENT(IN)                           :: seed, gen_name
        INTENT(INOUT)                        :: this
        !----------------------------------------------------------------------!
        CALL InitCommon(this%gentype, seed, gen_name)
        !  gen_name = name of the desired random number generator,
        !             possible values at the moment:
        ! 'MT' , 'TAUS88', 'COMPILER'

        SELECT CASE(GetGenName(this))

        CASE('MT')
            CALL init_MT(this, seed*16)
        CASE('TAUS88')
            CALL init_TAUS88(this, seed*11, seed*117, seed*16 )
        CASE('COMPILER')
            CALL RANDOM_SEED(N)
            ALLOCATE(seed_array(1:N))
            DO i = 1, N
                seed_array(i) = INT((seed * MNew) + i * IA + IC)
            END DO
            CALL RANDOM_SEED(PUT=seed_array)
            DEALLOCATE(seed_array)

        CASE DEFAULT
            PRINT *, 'ERROR, chosen random number generator not found.'
            STOP
        END SELECT

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

        SELECT CASE(GetGenName(this))

        CASE('MT')
            randno = getNumber_MT(this)
        CASE('TAUS88')
            randno = getNumber_TAUS88(this)
        CASE('COMPILER')
            CALL RANDOM_NUMBER(randno)
        CASE DEFAULT
            PRINT *, 'ERROR, chosen random number generator not initialized.'
            STOP
        END SELECT

    END SUBROUTINE GetNewRandomNumber

    PURE FUNCTION GetGenSeed(this) RESULT(ut)
        IMPLICIT NONE
        !----------------------------------------------------------------------!
        TYPE(Randgen_TYP), INTENT(IN) :: this
        INTEGER :: ut
        !----------------------------------------------------------------------!
        ut = GetType_common(this%gentype)
    END FUNCTION GetGenSeed


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
        !----------------------------------------------------------------------!
        i = Initialized_common(this%gentype)
    END FUNCTION RandgenInitialized

    SUBROUTINE init_TAUS88(this, i1, i2, i3)
        IMPLICIT NONE
        !----------------------------------------------------------------------!
        TYPE(Randgen_TYP), INTENT(INOUT) :: this
        INTEGER, INTENT(IN)              :: i1, i2, i3
        !----------------------------------------------------------------------!

        this%s1 = i1
        this%s2 = i2
        this%s3 = i3
        IF (IAND(this%s1,-2) == 0) this%s1 = i1 - 1023
        IF (IAND(this%s2,-8) == 0) this%s2 = i2 - 1023
        IF (IAND(this%s3,-16) == 0) this%s3 = i3 - 1023

    END SUBROUTINE init_TAUS88

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

        b  = ISHFT( IEOR( ISHFT(this%s1,13), this%s1), -19)
        this%s1 = IEOR( ISHFT( IAND(this%s1,-2), 12), b)
        b  = ISHFT( IEOR( ISHFT(this%s2,2), this%s2), -25)
        this%s2 = IEOR( ISHFT( IAND(this%s2,-8), 4), b)
        b  = ISHFT( IEOR( ISHFT(this%s3,3), this%s3), -11)
        this%s3 = IEOR( ISHFT( IAND(this%s3,-16), 17), b)
        random_numb = IEOR( IEOR(this%s1,this%s2), this%s3) *                  &
                           2.3283064365E-10_r2 + 0.5_r2
        this%rndx = random_numb
    END FUNCTION getNumber_TAUS88
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
     
    !   genrand() generates one pseudorandom real number (double) which is
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
    SUBROUTINE init_MT(this, seed)
        !----------------------------------------------------------------------!
        TYPE(Randgen_TYP), INTENT(INOUT) :: this
        INTEGER, INTENT(IN)              :: seed
        INTEGER                          :: latest
        INTEGER                          :: mti
        !----------------------------------------------------------------------!
        ! This initialization is based upon the multiplier given on p.106 of the
        ! 3rd edition of Knuth, The Art of Computer Programming Vol. 2.
        ! This version assumes that integer overflow does NOT cause a crash.

        this%mt(0) = seed
        latest = seed
        DO mti = 1, n-1
            this%mti = mti
            latest = IEOR( latest, ISHFT( latest, -30 ) )
            latest = latest * 1812433253 + this%mti
            this%mt(this%mti) = latest
        END DO

    END SUBROUTINE init_MT
    !***********************************************************************

    FUNCTION getNumber_MT(this) RESULT(fn_val)
        !----------------------------------------------------------------------!
        TYPE(Randgen_TYP), INTENT(INOUT) :: this
        INTEGER                          :: mag01(0:1) = (/ 0, mata /)
        INTEGER                          :: kk, y
        REAL (kind=r2)                   :: fn_val
        !----------------------------------------------------------------------!

        ! These statement functions have been replaced with separate functions
        ! tshftu(y) = ISHFT(y,-11)
        ! tshfts(y) = ISHFT(y,7)
        ! tshftt(y) = ISHFT(y,15)
        ! tshftl(y) = ISHFT(y,-18)

        IF(this%mti >= n) THEN
        !                       generate N words at one time
            IF(this%mti == n+1) THEN
        !                            if sgrnd() has not been called,
                CALL init_MT(this, 4357)
        !                              a default initial seed is used
            END IF
          
            DO  kk = 0, n-m-1
                y = IOR(IAND(this%mt(kk),umask), IAND(this%mt(kk+1),lmask))
                this%mt(kk) = IEOR(IEOR(this%mt(kk+m), ISHFT(y,-1)),mag01(IAND(y,1)))
            END DO
            DO  kk = n-m, n-2
                y = IOR(IAND(this%mt(kk),umask), IAND(this%mt(kk+1),lmask))
                this%mt(kk) = IEOR(IEOR(this%mt(kk+(m-n)), ISHFT(y,-1)),mag01(IAND(y,1)))
            END DO
            y = IOR(IAND(this%mt(n-1),umask), IAND(this%mt(0),lmask))
            this%mt(n-1) = IEOR(IEOR(this%mt(m-1), ISHFT(y,-1)),mag01(IAND(y,1)))
            this%mti = 0
        END IF

        y = this%mt(this%mti)
        this%mti = this%mti + 1
        y = IEOR(y, tshftu(y))
        y = IEOR(y, IAND(tshfts(y),tmaskb))
        y = IEOR(y, IAND(tshftt(y),tmaskc))
        y = IEOR(y, tshftl(y))

        IF(y < 0) THEN
          fn_val = (DBLE(y) + 2.0D0**32) / (2.0D0**32 - 1.0D0)
        ELSE
          fn_val = DBLE(y) / (2.0D0**32 - 1.0D0)
        END IF

        this%rndx = fn_val
    END FUNCTION getNumber_MT


    FUNCTION tshftu(y) RESULT(fn_val)
        INTEGER, INTENT(IN) :: y
        INTEGER             :: fn_val

        fn_val = ISHFT(y,-11)
        RETURN
    END FUNCTION tshftu


    FUNCTION tshfts(y) RESULT(fn_val)
        INTEGER, INTENT(IN) :: y
        INTEGER             :: fn_val

        fn_val = ISHFT(y,7)
        RETURN
    END FUNCTION tshfts


    FUNCTION tshftt(y) RESULT(fn_val)
        INTEGER, INTENT(IN) :: y
        INTEGER             :: fn_val

        fn_val = ISHFT(y,15)
        RETURN
    END FUNCTION tshftt

    FUNCTION tshftl(y) RESULT(fn_val)
        INTEGER, INTENT(IN) :: y
        INTEGER             :: fn_val

        fn_val = ISHFT(y,-18)
        RETURN
    END FUNCTION tshftl

END MODULE randgen_type
