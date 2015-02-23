!------------------------------------------------------------------------------!
! def Gas_TYP
! inspired by fosite by T. Illenseer 2011
!------------------------------------------------------------------------------!
MODULE gas_type
  
    USE datatype
    USE var_global
    USE common_type, &
        GetType_common => GetType, GetName_common => GetName, &
        Initialized_common => Initialized
    IMPLICIT NONE

    !--------------------------------------------------------------------------!
    PRIVATE
    ! 
    !--------------------------------------------------------------------------!
    TYPE Gas_TYP
        TYPE(Common_TYP) :: mtype                       ! -----------------    !
        !----------------------------------------------------------------------!
        REAL(kind=r1)                           :: dust_mol_ratio
        REAL(kind=r1)                           :: mol_abund
        REAL(kind=r1)                           :: vel_max
        
        REAL(kind=r1)                           :: mol_weight
        
        INTEGER                                 :: egy_lvl
        INTEGER                                 :: trans_lvl
        INTEGER                                 :: i_vel_chan
        INTEGER                                 :: n_tr
        INTEGER                                 :: col_partner
        INTEGER                                 :: col_trans
        INTEGER                                 :: col_temps
        
        REAL(kind=r1), DIMENSION(:), ALLOCATABLE       :: energylevel
        REAL(kind=r1), DIMENSION(:), ALLOCATABLE       :: g_level
        REAL(kind=r1), DIMENSION(:,:), ALLOCATABLE     :: col_alltemps
        REAL(kind=r1), DIMENSION(:,:,:), ALLOCATABLE   :: col_colmatrix

        
        INTEGER, DIMENSION(:,:), ALLOCATABLE           :: col_upper
        INTEGER, DIMENSION(:,:), ALLOCATABLE           :: col_lower
        INTEGER, DIMENSION(:), ALLOCATABLE          :: tr_cat
        INTEGER, DIMENSION(:), ALLOCATABLE          :: col_id
        INTEGER, DIMENSION(:), ALLOCATABLE          :: trans_upper
        INTEGER, DIMENSION(:), ALLOCATABLE          :: trans_lower
        REAL(kind=r1), DIMENSION(:), ALLOCATABLE    :: trans_einstA
        REAL(kind=r1), DIMENSION(:), ALLOCATABLE    :: trans_einstB_l
        REAL(kind=r1), DIMENSION(:), ALLOCATABLE    :: trans_einstB_u
        REAL(kind=r2), DIMENSION(:), ALLOCATABLE    :: trans_freq
        REAL(kind=r1), DIMENSION(:), ALLOCATABLE    :: trans_inneregy
        
        REAL(kind=r1), DIMENSION(:), ALLOCATABLE    :: velo_channel

    END TYPE Gas_TYP
    SAVE
    !--------------------------------------------------------------------------!
    
    PUBLIC :: &
        ! types
        Gas_TYP, &
        ! methods
        InitGas, &
        CloseGas, &
        GetGasType, &
        GetName, &
        GasInitialized
    !--------------------------------------------------------------------------!
CONTAINS

    SUBROUTINE InitGas(this, ut, un, moldustratio,                             &
                       nvratio, n_tr, tr_cat, i_vel_chan, vel_max)
        IMPLICIT NONE
        !----------------------------------------------------------------------!
        TYPE(Gas_TYP)                      :: this
        
        INTEGER                            :: ut, vch, n_tr
        INTEGER                            :: i_vel_chan
        INTEGER, DIMENSION(1:n_tr)         :: tr_cat
        CHARACTER(LEN=*)                   :: un
        
        REAL(kind=r2)                      :: moldustratio
        REAL(kind=r2)                      :: nvratio
        REAL(kind=r2)                      :: abundance
        REAL(kind=r2)                      :: vel_max

        !----------------------------------------------------------------------!
        INTENT(IN)             :: ut,un, moldustratio, nvratio, n_tr, tr_cat, &
                                  i_vel_chan, vel_max
        INTENT(INOUT)          :: this
        !----------------------------------------------------------------------!
        CALL InitCommon(this%mtype,ut,un)
    
        this%i_vel_chan     = i_vel_chan
        this%vel_max        = vel_max
        
        ALLOCATE(   this%velo_channel(-this%i_vel_chan:this%i_vel_chan),& 
                    this%tr_cat(1:n_tr) )
        this%n_tr = n_tr
        this%tr_cat = tr_cat
        
        DO vch = -this%i_vel_chan,this%i_vel_chan  
            this%velo_channel(vch) = this%vel_max*real(vch,kind=r2)            &
                                     /REAL(this%i_vel_chan,kind=r2)
        END DO
        
        
        CALL READ_mol_file(this)
        
        this%dust_mol_ratio = moldustratio
        abundance           = nvratio
        this%mol_abund      = abundance
        
    END SUBROUTINE InitGas

    SUBROUTINE READ_mol_file(this)
    
        IMPLICIT NONE
        !----------------------------------------------------------------------!
        TYPE(Gas_TYP), INTENT(INOUT) :: this
        !----------------------------------------------------------------------!
        CHARACTER(len=32)            :: mol_name        
        CHARACTER(len=32)            :: hch1
        
        INTEGER                      :: i, k, waste_i
        !----------------------------------------------------------------------!
          
        mol_name = GetName(this)
        
        OPEN(unit=1, file=path_mol//TRIM(mol_name), &
            action="read", status="unknown", form="formatted")
            
        read(unit=1,fmt=*) hch1                      !CH Molecule
        read(unit=1,fmt=*) hch1                      !Moleculename
        read(unit=1,fmt=*) hch1                      !CH Molecular weight
        read(unit=1,fmt=*) this%mol_weight           !Molecular weight number
        read(unit=1,fmt=*) hch1                      !CH Number of energy Levels
        read(unit=1,fmt=*) this%egy_lvl              !Number of energy Levels
        
        ALLOCATE(   this%energylevel( 1:this%egy_lvl), &
                    this%g_level( 1:this%egy_lvl) )
        
        READ(unit=1,fmt=*) hch1                      !CH LEVELS
        
        DO i = 1, this%egy_lvl
            read(unit=1,fmt=*) waste_i, this%energylevel(i), this%g_level(i)
        END DO
        
        READ(unit=1,fmt=*) hch1              !CH Number of radiative Transitions
        READ(unit=1,fmt=*) this%trans_lvl    !Number of radiative Transitions
        
        ALLOCATE(  this%trans_upper( 1:this%trans_lvl), &
                    this%trans_lower( 1:this%trans_lvl), &
                    this%trans_freq( 1:this%trans_lvl), &
                    this%trans_einstA( 1:this%trans_lvl), &
                    this%trans_einstB_l( 1:this%trans_lvl), &
                    this%trans_einstB_u( 1:this%trans_lvl), &
                    this%trans_inneregy( 1:this%trans_lvl) )
        
        READ(unit=1,fmt=*) hch1              !CH TRANSITIONS
        DO i = 1, this%trans_lvl
            
            read(unit=1,fmt=*)  waste_i, this%trans_upper(i),                  &
                                this%trans_lower(i), this%trans_einstA(i),     &
                                this%trans_freq(i), this%trans_inneregy(i)
                                
            this%trans_freq(i) =  this%trans_freq(i)*1e9_r2  
            
            
            ! calculate einstein B _u(pper), _l(ower) levels
            this%trans_einstB_u(i) = this%trans_einstA(i) *                    &
                                     (con_c/(this%trans_freq(i)))**2 /         &
                                     (2.0_r2*con_h*(this%trans_freq(i)))
                                     
            this%trans_einstB_l(i) = this%g_level(this%trans_upper(i)) /       &
                                     this%g_level(this%trans_lower(i)) *       &
                                     this%trans_einstB_u(i)
        END DO

        ! Read all collision data from file, 
        !
        READ(unit=1,fmt=*) hch1                  !CH number of collision partner
        READ(unit=1,fmt=*) this%col_partner
        ALLOCATE (this%col_id(1:this%col_partner))
        this%col_id = 0
        DO k = 1, this%col_partner
        
            READ(unit=1,fmt=*) hch1                   !CH collision between
            READ(unit=1,fmt=*) this%col_id(k)         !nr cols id

            READ(unit=1,fmt=*) hch1                   !CH number col transitions
            READ(unit=1,fmt=*) waste_i
            READ(unit=1,fmt=*) hch1                   !CH number temperatures
            READ(unit=1,fmt=*) this%col_temps
            IF (k ==1) THEN
                this%col_trans = waste_i
                ALLOCATE(this%col_alltemps(1:this%col_temps,1:6))
                ALLOCATE(this%col_colmatrix(1:6,1:this%col_trans,1:this%col_temps) )
                ALLOCATE(this%col_upper(1:this%col_trans,1:6))
                ALLOCATE(this%col_lower(1:this%col_trans,1:6))
                this%col_alltemps  = 0.0_r2
                this%col_colmatrix = 0.0_r2
                this%col_upper     = 0
                this%col_lower     = 0
            END IF
            ! check if all col partner have the same number of transitions 
            ! if False, stop the code. This should be done in future!
            IF (this%col_trans /= waste_i) THEN
                PRINT *, 'ERROR: The number of collistional & 
                          &transitions is not &
                          &the same for all col partners.'
                PRINT *, 'Solution: Choose an other molecule like CO'
                STOP
            END IF
            READ(unit=1,fmt=*) hch1     !CH all temps
            READ(unit=1,fmt=*) this%col_alltemps(:,this%col_id(k))
            READ(unit=1,fmt=*) hch1     !CH TRANS+ UP+ LOW+ COLLRATES(cm^3 s^-1)
            DO i = 1, this%col_trans
                READ(unit=1,fmt=*) waste_i,                                    &
                                   this%col_upper(i,this%col_id(k)),           &
                                   this%col_lower(i,this%col_id(k)),           &
                                   this%col_colmatrix(this%col_id(k),i,:)
                this%col_colmatrix(this%col_id(k),i,:) =                       &
                                   this%col_colmatrix(this%col_id(k),i,:) *    &
                                   1.0e-6_r2
            END DO
        END DO !end of col partner
        CLOSE(unit=1)
    END SUBROUTINE READ_mol_file



    SUBROUTINE CloseGas(this)
        IMPLICIT NONE
        !----------------------------------------------------------------------!
        TYPE(Gas_TYP), INTENT(INOUT) :: this
        !----------------------------------------------------------------------!
        CALL CloseCommon(this%mtype)
        DEALLOCATE( this%energylevel, &
                    this%trans_einstA, &
                    this%trans_einstB_l, &
                    this%trans_einstB_u, &
                    this%trans_upper, &
                    this%trans_lower, &
                    this%trans_freq, &
                    this%trans_inneregy, &
                    this%velo_channel, &
                    this%g_level, &
                    this%col_id, &
                    this%tr_cat)
        
    END SUBROUTINE CloseGas


    PURE FUNCTION GetGasType(this) RESULT(ut)
        IMPLICIT NONE
        !----------------------------------------------------------------------!
        TYPE(Gas_TYP), INTENT(IN) :: this
        INTEGER                   :: ut
        !----------------------------------------------------------------------!
        ut = GetType_common(this%mtype)
    END FUNCTION GetGasType


    PURE FUNCTION GetName(this) RESULT(un)
        IMPLICIT NONE
        !----------------------------------------------------------------------!
        TYPE(Gas_TYP), INTENT(IN) :: this
        CHARACTER(LEN=32)         :: un
        !----------------------------------------------------------------------!
        un = GetName_common(this%mtype)
    END FUNCTION GetName
    

    PURE FUNCTION GasInitialized(this) RESULT(i)
        IMPLICIT NONE
          !--------------------------------------------------------------------!
          TYPE(Gas_TYP), INTENT(IN) :: this
          LOGICAL                   :: i
          !--------------------------------------------------------------------!
          i = Initialized_common(this%mtype)
    END FUNCTION GasInitialized
    

End Module gas_type
