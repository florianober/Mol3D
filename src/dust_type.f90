!------------------------------------------------------------------------------!
! def Dust_TYP
! inspired by fosite by T. Illenseer 2011
!------------------------------------------------------------------------------!
MODULE Dust_type
  
    USE datatype
    USE gas_type
    USE basic_type
    USE model_type
    USE photon_type
    USE var_global
    USE math_mod
    USE common_type, &
        GetType_common => GetType, GetName_common => GetName, &
        Initialized_common => Initialized
    IMPLICIT NONE

    !--------------------------------------------------------------------------!
    PRIVATE
    ! 
    !--------------------------------------------------------------------------!
    TYPE Dust_TYP
        TYPE(Common_TYP) :: dttype                       ! -----------------   !
        !----------------------------------------------------------------------!
        INTEGER                                            :: n_dust
        INTEGER                                            :: n_lam
        INTEGER                                            :: n_lam_map
        INTEGER                                            :: aniiso      
        INTEGER                                            :: n_scatt_th      
        INTEGER                                            :: nrndpt
        INTEGER, DIMENSION(:),ALLOCATABLE                  :: num_lam_map
        INTEGER, DIMENSION(:), ALLOCATABLE                 :: cont_map
        INTEGER, DIMENSION(:,:,:),ALLOCATABLE              :: SCAANG
        
        REAL(kind=r2)                                      :: sizexp
        REAL(kind=r2)                                      :: D_ANG
        REAL(kind=r2), DIMENSION(:),ALLOCATABLE            :: den_dust
        REAL(kind=r2), DIMENSION(:),ALLOCATABLE            :: r_dust
        REAL(kind=r2), DIMENSION(:),ALLOCATABLE            :: sidi
        REAL(kind=r2), DIMENSION(:,:),ALLOCATABLE          :: sidi_par
        
        REAL(kind=r2), DIMENSION(:),ALLOCATABLE            :: lam
        REAL(kind=r2), DIMENSION(:),ALLOCATABLE            :: d_lam
        REAL(kind=r2), DIMENSION(:),ALLOCATABLE            :: nu
        REAL(kind=r2), DIMENSION(:),ALLOCATABLE            :: d_nu
        REAL(kind=r2), DIMENSION(:,:),ALLOCATABLE          :: sizepar
        REAL(kind=r2), DIMENSION(:,:),ALLOCATABLE          :: C_sca
        REAL(kind=r2), DIMENSION(:,:),ALLOCATABLE          :: C_ext
        REAL(kind=r2), DIMENSION(:,:),ALLOCATABLE          :: C_abs
        REAL(kind=r2), DIMENSION(:,:),ALLOCATABLE          :: albedo
        REAL(kind=r2), DIMENSION(:,:),ALLOCATABLE          :: QB
        REAL(kind=r2), DIMENSION(:,:),ALLOCATABLE          :: QdB_dT
        REAL(kind=r2), DIMENSION(:,:,:),ALLOCATABLE        :: QdB_dT_l
        REAL(kind=r2), DIMENSION(:,:,:),ALLOCATABLE        :: QdB_dT_l_cdf
        REAL(kind=r2), DIMENSION(:,:),ALLOCATABLE          :: planck_tab
        REAL(kind=r2), DIMENSION(:,:,:,:,:),ALLOCATABLE    :: SME
        REAL(kind=r2), DIMENSION(:),ALLOCATABLE            :: sctth
        REAL(kind=r1), DIMENSION(:),ALLOCATABLE            :: tem_tab
        
        CHARACTER(len=8), DIMENSION(:),ALLOCATABLE         :: dust_cat
        CHARACTER(len=12), DIMENSION(:,:),ALLOCATABLE      :: dust_single

    END TYPE Dust_TYP
    SAVE
    !--------------------------------------------------------------------------!
    
    PUBLIC :: &
        ! types
        Dust_TYP, &
        ! methods
        InitDust, &
        CloseDust, &
        GetDustType, &
        GetDustName, &
        DustInitialized, &
        dust_select
    !--------------------------------------------------------------------------!
CONTAINS

    SUBROUTINE InitDust(this,basics , gas, model, ut,un, ndust,dendust,        &
                        sizexp, aniso, dust_cat, n_scatt_th, nrndpt)
        
        USE math_mod
        IMPLICIT NONE
        !----------------------------------------------------------------------!
        TYPE(Basic_TYP)  , INTENT(IN)         :: basics
        TYPE(Dust_TYP)                        :: this
        TYPE(Gas_TYP)                         :: gas
        TYPE(Model_TYP), INTENT(IN)           :: model
        
        INTEGER                               :: ut
        INTEGER                               :: hi1
        INTEGER                               :: tr
        INTEGER                               :: ptr_1
        INTEGER                               :: ptr_2
        INTEGER                               :: ndust
        INTEGER                               :: i_dust
        INTEGER                               :: i_lam
        INTEGER                               :: i_lam_map

        INTEGER                               :: aniso        
        INTEGER                               :: n_scatt_th   
        INTEGER                               :: nrndpt
        INTEGER                               :: i_scatt_th   

        REAL(kind=r2)                         :: sizexp
        REAL(kind=r2)                         :: hd1, hd2, hd_dth, hd_th
        REAL(kind=r2)                         :: N_AN
        REAL(kind=r2)                         :: REF_MED
        REAL(kind=r2), DIMENSION(0:180)       :: phg
        REAL(kind=r2), DIMENSION(-1:180)      :: Fphg
        REAL(kind=r2), DIMENSION(ndust)       :: dendust

        REAL(kind=r2), DIMENSION(:,:), ALLOCATABLE    :: gscatt_all
        REAL(kind=r2), DIMENSION(:,:), ALLOCATABLE    :: REF_RE
        REAL(kind=r2), DIMENSION(:,:), ALLOCATABLE    :: REF_IM
        REAL(kind=r2), DIMENSION(:,:), ALLOCATABLE    :: C_back
        
        CHARACTER(len=8), DIMENSION(ndust)    :: dust_cat
        CHARACTER(LEN=*)                      :: un
    
        !----------------------------------------------------------------------!
        INTENT(IN)          :: ut,un, ndust, dendust, sizexp, aniso, dust_cat, &
                               n_scatt_th, nrndpt
        INTENT(INOUT)       :: this, gas
        !----------------------------------------------------------------------!
        CALL InitCommon(this%dttype,ut,un)
        
        this%n_dust = ndust
        ALLOCATE( this%dust_cat( 1:this%n_dust) )
        this%dust_cat = dust_cat       
        open(unit=1, file=path_dust_cat//this%dust_cat(1)//".dct", &
            action="read", status="unknown", form="formatted")
            read(unit=1,fmt=*) this%n_lam
        close(unit=1)
        
        this%n_lam_map          = this%n_lam
        this%aniiso             = aniso
        this%sizexp             = sizexp
        this%n_scatt_th         = n_scatt_th
        this%nrndpt             = nrndpt
        
        ALLOCATE( this%den_dust( 1:this%n_dust), &
            this%r_dust( 1:this%n_dust), &
            this%sidi( 1:this%n_dust), &
            this%sidi_par( 1:this%n_dust,1:3), &
            this%lam( 1:this%n_lam ), & 
            this%d_lam( 1:this%n_lam ), & 
            this%nu( 1:this%n_lam ), & 
            this%d_nu( 1:this%n_lam ), & 
            this%num_lam_map( 1:this%n_lam_map ), & 
            this%dust_single( 1:this%n_dust, 1:this%n_lam ), &
            this%sizepar(     1:this%n_dust, 1:this%n_lam ), &
            this%C_sca(       1:this%n_dust, 1:this%n_lam ), &
            this%C_ext(       1:this%n_dust, 1:this%n_lam ), &
            this%C_abs(       1:this%n_dust, 1:this%n_lam ), &
            this%albedo(      1:this%n_dust, 1:this%n_lam ), &
            this%QdB_dT_l(    1:this%n_lam,0:basics%n_tem, 1:this%n_dust ),    &
            this%QdB_dT_l_cdf( 1:this%n_lam,0:basics%n_tem, 1:this%n_dust ),   &
            this%QB(       -1:basics%n_tem,1:this%n_dust ),                    &
            this%QdB_dT(        0:basics%n_tem,1:this%n_dust ),                &
            this%planck_tab(  1:this%n_lam,0:basics%n_tem ),                   &

            this%tem_tab(      0:basics%n_tem ),                               &
            this%SME(1:4, 1:4, 1:this%n_dust, this%n_lam, this%n_scatt_th),    &
            this%SCAANG(1:this%n_dust,1:this%n_lam,0:this%nrndpt),             &
            this%cont_map(1:gas%n_tr),                                         &
            this%sctth(1:n_scatt_th) )

        
        this%den_dust = dendust
        
        this%r_dust(:)          = 0.0_r2
        this%sidi(:)            = 0.0_r2
        this%sidi_par(:,:)      = 0.0_r2
        this%lam(:)             = 0.0_r2
        this%d_lam(:)           = 0.0_r2
        this%dust_single( :,: ) = ''
        this%sizepar(     :,: ) = 0.0_r2
        this%C_sca(       :,: ) = 0.0_r2
        this%C_ext(       :,: ) = 0.0_r2
        this%C_abs(       :,: ) = 0.0_r2
        this%albedo(      :,: ) = 0.0_r2
        this%SME(    :,:,:,:,:) = 0.0_r2
        this%SCAANG(     :,:,:) = 0
        this%sctth(:)           = 0.0_r2
        this%QB(:,:)            = 0.0_r2
        this%QdB_dT(:,:)        = 0.0_r2
        this%QdB_dT_l(:,:,:)    = 0.0_r2
        this%QdB_dT_l_cdf(:,:,:)    = 0.0_r2
        this%tem_tab(:)         = 0.0_r2

        ALLOCATE(gscatt_all(  1:this%n_dust, 1:this%n_lam ),                   &
                 REF_RE(      1:this%n_dust, 1:this%n_lam ),                   &
                 REF_IM(      1:this%n_dust, 1:this%n_lam ),                   &
                 C_back(      1:this%n_dust, 1:this%n_lam ) )
                 
        gscatt_all(:, :) = 0.0_r2
        REF_RE(:, :)     = 0.0_r2
        REF_IM(:, :)     = 0.0_r2
        C_back(:, :)     = 0.0_r2

        ! 1. read dust name catalogues
        do i_dust = 1, this%n_dust
            open(unit=2, file=path_dust_cat//this%dust_cat(i_dust)//".dct", &
            action="read", status="unknown", form="formatted")
            read(unit=2,fmt=*)
            
            do i_lam=1, this%n_lam
                read(unit=2,fmt=*) this%dust_single(i_dust,i_lam)
            end do
            close(unit=2)
        end do

        ! 2. read individual files
        !    NOTE: this routine is a mainly unmodified version from MC3D.
        !          it is working for 1 dust grain, more are feasible in
        !          princible, but not tested.
        do i_dust=1, this%n_dust
!~             if (this%n_dust>1) then
!~                 print *,"   > type of dust [=f(chemistry,size)] :", i_dust, " / ", this%n_dust
!~             end if

            do i_lam=1, this%n_lam
                ! ---
                ! 2.1. load data from file
                open(unit=1, file=path_dust_single//this%dust_single(i_dust,i_lam)//".lot", &
                action="read", status="unknown", form="formatted")
          
                ! handling of anisotropy: aniso == 1: anisotropic       scattering
                !                         aniso == 2: isotropic         scattering
                !                         aniso == 3: henyey-greenstein scattering
                ! the desired handling of anisotropy is choosen in the init file
                ! TbD: a test to check if henyey-greenstein scattering is working

                read(unit=1,fmt=*) hi1      ! scatt. matrix / henyey-greenstein (obsolete)

                !if (this%aniiso==1) then    ! before: only anisotropy chosen
                !    this%aniiso = hi1       ! now specify handling of anisotropy (aniiso == 1 or aniiso == 3)
                !end if

                ! wavelength/ frequency
                ! Tbd: this should be generalized. We should use a frequency
                !      table and interpolate for every frequency in the MC
                !      process. Otherwise we introduce errors due to the
                !      energy integration and unit conversion.
                !      
                hd1 = this%lam(i_lam)
                read(unit=1,fmt=*)  this%lam(i_lam)
                this%lam(i_lam) = 1.0e-6_r2 * this%lam(i_lam)
                this%nu(i_lam) = con_c/this%lam(i_lam)
                ! set wavelength/frequency intervals

                IF (i_lam > 1 ) THEN
                    this%d_lam(i_lam) = (this%lam(i_lam)-this%lam(i_lam-1))
                    this%d_nu(i_lam) = this%nu(i_lam-1) - this%nu(i_lam)
                END IF

                ! verify: wavelength correct?
                if ((this%n_dust>1) .and. (i_dust>1)) then
                    if (((hd1-this%lam(i_lam))/hd1) > 1.0e-3) then
                        print *,"!!! Warning: subroutine ld_dust():"
                        print *,"    Wavelengths of different kinds of dust do not agree with each other"
                        stop
                    end if
                end if

                ! grain size distribution
                ! smallest radius [micron]
                read(unit=1,fmt=*)  this%sidi_par(i_dust, 1)
                ! largest radius  [micron]
                read(unit=1,fmt=*)  this%sidi_par(i_dust, 2)
                ! size distribution exp. (<0)
                read(unit=1,fmt=*)  this%sidi_par(i_dust, 3)
                ! [micron] -> [m]
                this%sidi_par(i_dust,1:2) = this%sidi_par(i_dust,1:2) * 1.0e-6_r2

                ! HG scattering g
                read(unit=1,fmt=*)  gscatt_all(i_dust, i_lam)
          
                ! environment
                read(unit=1,fmt=*)  REF_MED

                ! distribution of scattering angles
                read(unit=1,fmt=*)  this%D_ANG
                read(unit=1,fmt=*)  N_AN

                ! dust grain: refractive index
                read(unit=1,fmt=*)  REF_RE( i_dust,i_lam )
                read(unit=1,fmt=*)  REF_IM( i_dust,i_lam )

                ! radius & size parameter
                ! r_dust(): used to calculate the dust mass

                read(unit=1,fmt=*)  this%r_dust(i_dust)
                
                if (this%sidi_par(i_dust,3)==0.0_r2) then
                    hd2 = this%sidi_par(i_dust,1)
                else
                    hd2 = (1.0_r2+this%sidi_par(i_dust,3)) / (4.0_r2+this%sidi_par(i_dust,3)) &
                        * (this%sidi_par(i_dust,2)**(4.0_r2+this%sidi_par(i_dust,3)) - &
                        this%sidi_par(i_dust,1)**(4.0_r2+this%sidi_par(i_dust,3))) &
                        / (this%sidi_par(i_dust,2)**(1.0_r2+this%sidi_par(i_dust,3)) - &
                        this%sidi_par(i_dust,1)**(1.0_r2+this%sidi_par(i_dust,3)))
                    hd2 = hd2**(1.0_r2/3.0_r2)
                end if
                this%r_dust(i_dust) = hd2                                 ! [m]

                if (i_lam /= 1) then
                    if (this%r_dust(i_dust) /= hd2) then
                        print *,"subroutine ld_dust(): wrong grain size" 
                        print *,"              i_lam : ", i_lam
                        print *,"              i_dust: ", i_dust
                        print *,"            required: ", hd2
                        print *,"           available: ", this%r_dust(i_dust)
                        stop
                    end if
                end if
                read(unit=1,fmt=*)  this%sizepar(i_dust, i_lam)

                ! scattering cross sections
                read(unit=1,fmt=*)  this%C_sca(  i_dust, i_lam)
                read(unit=1,fmt=*)  this%C_ext(  i_dust, i_lam)
                read(unit=1,fmt=*)  C_back( i_dust, i_lam)

                this%C_abs( i_dust,i_lam) = this%C_ext(i_dust,i_lam) - this%C_sca(i_dust,i_lam)
                this%albedo(i_dust,i_lam) = this%C_sca(i_dust,i_lam) / this%C_ext(i_dust,i_lam)

                do i_scatt_th=1, n_scatt_th
                    read(unit=1,fmt=*)  this%SME(1,1, i_dust,i_lam, i_scatt_th)
                    read(unit=1,fmt=*)  this%SME(1,2, i_dust,i_lam, i_scatt_th)
                    read(unit=1,fmt=*)  this%SME(3,3, i_dust,i_lam, i_scatt_th)
                    read(unit=1,fmt=*)  this%SME(3,4, i_dust,i_lam, i_scatt_th)

                    this%SME(2,2, i_dust,i_lam, i_scatt_th) =  this%SME(1,1, i_dust,i_lam, i_scatt_th)
                    this%SME(2,1, i_dust,i_lam, i_scatt_th) =  this%SME(1,2, i_dust,i_lam, i_scatt_th)
                    this%SME(4,3, i_dust,i_lam, i_scatt_th) = -this%SME(3,4, i_dust,i_lam, i_scatt_th)
                    this%SME(4,4, i_dust,i_lam, i_scatt_th) =  this%SME(3,3, i_dust,i_lam, i_scatt_th)
                end do
          
                close(unit=1)
             
                ! ---
                ! 2.2. prepare scattering distribution function
                select case(this%aniiso)

                case (1) ! --- anisotropic scattering: Mie theory ---
                    if (N_AN/=181) then
                        print *,"!!! Warning: subroutine ld_dust()"
                        print *,"             mie-scattering prepared for 181 scattering angles only"
                        stop
                    end if

                    phg(:)  = 0.0_r2
                    Fphg(:) = 0.0_r2
                    hd_dth  = grad2rad(0.5_r2)

                    do i_scatt_th = 0, 180
                        if (i_scatt_th/=0 .and. i_scatt_th/=180) then
                            hd_th = grad2rad(real(i_scatt_th,kind=r2))
                            phg(i_scatt_th) = 2.0_r2*PI * real(this%SME(1,1, i_dust,i_lam, i_scatt_th+1),kind=r2) * &
                            (cos(hd_th-hd_dth) - cos(hd_th+hd_dth))
                        else
                            phg(i_scatt_th) = 2.0_r2*PI * real(this%SME(1,1, i_dust,i_lam, i_scatt_th+1),kind=r2) * &
                            (1.0_r2 - cos(2.0_r2*hd_dth))
                        end if
                
                        Fphg(i_scatt_th) = Fphg(i_scatt_th-1) + phg(i_scatt_th)
                    end do
                    Fphg(:) = Fphg(:) / Fphg(180)
             
                    ptr_1 = 0
                    ptr_2 = 0
                    do i_scatt_th = 0, 180
                        ptr_1 = int( Fphg(i_scatt_th-1) * real(nrndpt,kind=r2) )
                
                        if (i_scatt_th<180) then
                            ptr_2 = int( Fphg(i_scatt_th) * real(nrndpt,kind=r2) )
                            if (ptr_2>nrndpt) then
                                ptr_2 = nrndpt
                            end if
                        else
                            ptr_2 = nrndpt
                        end if
                
                        this%SCAANG(i_dust,i_lam, ptr_1:ptr_2) = i_scatt_th +1
                    end do

                case (2) ! --- isotropic scattering  ---
                    ! SCAANG: not defined

                case (3) ! --- anisotropic scattering: HENYEY GREENSTEIN approximation
                    ! assumption: streuwinkel = 0, 1, 2, ..., 180 degree
                    ! resulting <g> agrees with the g scattering parameter used as input parameter
                    if (N_AN/=181) then
                        print *,"!!! Warning: subroutine ld_dust()"
                        print *,"             mie-scattering prepared for 181 scattering angles only"
                        stop
                    end if

                    phg(:)  = 0.0_r2
                    Fphg(:) = 0.0_r2
                 
                    do i_scatt_th = 0, 180
                        hd_th = grad2rad(real(i_scatt_th,kind=r2))
                        phg(i_scatt_th) = (1.0_r2 - gscatt_all(i_dust,i_lam)**2) / &
                            ( 1.0_r2   &
                             + gscatt_all(i_dust,i_lam)**2   &
                             - (2.0_r2*gscatt_all(i_dust,i_lam) * cos(hd_th))   &
                             )**1.5
                        Fphg(i_scatt_th) = Fphg(i_scatt_th-1) + phg(i_scatt_th)*sin(hd_th)
                    end do
                    Fphg(:) = Fphg(:) / Fphg(180)
                 
                    ptr_1 = 0
                    ptr_2 = 0
                    do i_scatt_th = 0, 180
                        ptr_1 = int( Fphg(i_scatt_th-1) * real(nrndpt,kind=r2) )
                    
                        if (i_scatt_th<180) then
                            ptr_2 = int( Fphg(i_scatt_th) * real(nrndpt,kind=r2) )
                            if (ptr_2>nrndpt) then
                                ptr_2 = nrndpt
                            end if
                        else
                            ptr_2 = nrndpt
                        end if
                        this%SCAANG(i_dust,i_lam, ptr_1:ptr_2) = i_scatt_th
                    end do
                case default
                    print *,"<!> subroutine ld_dust: select case(aniiso): wrong input parameter. stopped."
                    stop
                end select
            end do ! end do i_lam
        end do   !end do i_dust

        ! ---
        ! 2.3. pre-calculation of some parameters

        ! scattering angles (D_ANG = delta.theta, theta >= 0°, [sctth(:)] = rad
        do i_scatt_th=1, this%n_scatt_th
            this%sctth(i_scatt_th) = real(i_scatt_th-1, kind=r2) * this%D_ANG
        end do
    
        ! particle size distribution <-> particle number distribution
        ! sidi ... rel. number of dust grains with radius r_dust(i_dust)
        this%sidi(:) = ( this%r_dust(:)/this%r_dust(this%n_dust) ) ** this%sizexp
        this%sidi(:) = this%sidi(:) / sum( this%sidi(:) )

        do i_lam_map=1, this%n_lam_map
            this%num_lam_map(i_lam_map) = i_lam_map
        end do
        
        ! select the suited dust wavelength for each tr(ansition) for the selected line transition
        DO tr = 1, gas%n_tr
          IF (con_c/gas%trans_freq(gas%tr_cat(tr)) .GT. MAXVAL(this%lam)) THEN
               this%cont_map(tr) = MAXLOC(this%lam, 1) !position of the maximum wavelength in the array lam
               print *, 'WARNING: the central wavelength of the gas is higher than the maximum &
                                &  wavelength provided by the choosen dust-catalogue'
            ELSEIF (con_c/gas%trans_freq(gas%tr_cat(tr)) .LT. MINVAL(this%lam)) THEN
                 this%cont_map(tr) = MINLOC(this%lam, 1) !position of the minimum wavelength in the array lam
                 print *, 'WARNING: the central wavelength of the gas is lower than the minimum &
                                   & wavelength provided by the choosen dust-catalogue'
            ELSE
            hi1 = binary_search(con_c/gas%trans_freq(gas%tr_cat(tr)),this%lam(:))
            
            IF (this%lam(hi1+1) .gt. 2*con_c/gas%trans_freq(gas%tr_cat(tr)) - this%lam(hi1) ) THEN
                this%cont_map(tr) = hi1
            ELSE
                this%cont_map(tr) = hi1+1
            END IF
          END IF
        END DO
        CALL prepare_QB(this,basics)
    
        DEALLOCATE(gscatt_all, REF_RE, REF_IM, C_back)
    END SUBROUTINE InitDust

    
    SUBROUTINE prepare_QB(this,basics)

        IMPLICIT NONE
        
        !----------------------------------------------------------------------!
        TYPE(Dust_TYP)                                 :: this
        TYPE(Basic_TYP),INTENT(IN)                     :: basics
        !----------------------------------------------------------------------!
        INTEGER                                        :: i_lam, i_dust, i_tem
        REAL(kind=r2), DIMENSION(:), allocatable       :: QBx
        REAL(kind=r2)                                  :: hd1
        !-----------------------------------------------------------------------
        
        ! ---
        print *,"tabulating: - right side of energy equation"
        print *,"            - planck function               [this may take a while]"

        allocate( QBx(1:this%n_lam) )

        this%QB(-1,:) = -1.0_r2
        
        do i_dust=1, this%n_dust
            do i_tem=0, basics%n_tem
                ! temperature (tabulate)
                ! tbd: logarithmic distribution
                this%tem_tab(i_tem) = real(i_tem, kind=r2) * basics%d_tem  +  basics%t_dust_min
                ! tabulating planck function
                this%planck_tab(:,i_tem) = planck( this%tem_tab(i_tem), this%lam(:))

                ! individual contributions C_abs * B(T)
                QBx(:) =  this%C_abs(i_dust,:) * this%planck_tab(:,i_tem)

                ! integral [W]
                this%QB(i_tem, i_dust) = integ1(this%lam(:), QBx(:), 1, this%n_lam)

                ! individual contributions C_abs * dB(T)/dT
                ! (needed for Bjorkman & Wood temperature/wavelength correction)
                this%QdB_dT_l(:,i_tem,i_dust) =  this%C_abs(i_dust,:) * dB_dT_l(this%tem_tab(i_tem), this%lam(:))
                ! integral [W]
                this%QdB_dT(i_tem,i_dust) = integ1(this%lam(:),this%QdB_dT_l(:,i_tem,i_dust) , 1, this%n_lam)

                ! normalize -> cdf ->  used to get the new wavelength in B & W approach
                hd1 = 0.0_r2
                IF (this%QdB_dT(i_tem,i_dust) .lt. 1.0e-301_r2) THEN
                    this%QdB_dT_l_cdf(:,i_tem,i_dust) = 0.0_r2
                ELSE
                    DO i_lam=2, this%n_lam
                        hd1 = hd1 +                                            &
                            integ1(this%lam(:),this%QdB_dT_l(:,i_tem,i_dust)/  &
                            this%QdB_dT(i_tem,i_dust),i_lam-1,i_lam)
                        this%QdB_dT_l_cdf(i_lam,i_tem,i_dust) =  hd1
                    END DO
                END IF
           end do
        end do
        deallocate(QBx)

    END SUBROUTINE prepare_QB

    SUBROUTINE CloseDust(this)
        IMPLICIT NONE
        !----------------------------------------------------------------------!
        TYPE(Dust_TYP), INTENT(INOUT) :: this
        !----------------------------------------------------------------------!
        CALL CloseCommon(this%dttype)
        DEALLOCATE( this%dust_cat, &
                    this%den_dust, &
                    this%r_dust, &
                    this%sidi, &
                    this%sidi_par, &
                    this%lam, & 
                    this%d_lam, & 
                    this%nu, & 
                    this%d_nu, & 
                    this%dust_single, &
                    this%cont_map, &
                    this%sizepar, &
                    this%C_sca, &
                    this%C_ext, &
                    this%C_abs, &
                    this%albedo, &
                    this%planck_tab, &
                    this%QB, &
                    this%QdB_dT, &
                    this%QdB_dT_l, &
                    this%QdB_dT_l_cdf, &
                    this%SME, &
                    this%SCAANG, &
                    this%tem_tab, &
                    this%sctth)
        
    END SUBROUTINE CloseDust

    ! ##########################################################################
    ! select this species in current cell
    ! ---
    
    FUNCTION dust_select(grid, this, rand_nr, photon) RESULT(i_dust_action)
        USE randgen_type
        USE grid_type
        IMPLICIT NONE
        !----------------------------------------------------------------------!
        TYPE(Dust_TYP),INTENT(IN)                        :: this
        
        TYPE(Randgen_TYP),INTENT(INOUT)                  :: rand_nr
        TYPE(Grid_TYP),INTENT(IN)                        :: grid
        
        TYPE(PHOTON_TYP),INTENT(IN)                      :: photon
        !----------------------------------------------------------------------!
        INTEGER                                          :: i_dust_action
        REAL(kind=r2),DIMENSION(1:this%n_dust)           :: prob_action
        REAL(kind=r2)                                    :: rndx
        real(kind=r2) :: hd1, hd2
        !----------------------------------------------------------------------!
        
        ! ---
        IF (this%n_dust==1) THEN
           i_dust_action = 1
        ELSE
           ! approach: probability of interaction of the radiation with a this 
           !           grain species is in direct proportion to the
           !           a) the number density of these this grains in the cell
           !      and  b) the extinction cross section of that species
           prob_action(:) = grid%grd_dust_density(photon%nr_cell,:) *          &
                                 this%C_ext(:,photon%nr_lam)

           CALL RAN2(rand_nr, rndx)    
           hd1 = rndx * sum( prob_action(:) )

           i_dust_action = 1
           hd2           = prob_action(i_dust_action)  +  hd1 * epsilon(1.0_r2)
           DO
              IF (hd2 >= hd1 ) THEN
                 EXIT
              ELSE
                 i_dust_action = i_dust_action + 1
                 hd2           = hd2           + prob_action(i_dust_action)
              END IF
           END DO
        END IF
    END FUNCTION dust_select

    PURE FUNCTION GetDustType(this) RESULT(ut)
        IMPLICIT NONE
        !----------------------------------------------------------------------!
        TYPE(Dust_TYP), INTENT(IN) :: this
        INTEGER :: ut
        !----------------------------------------------------------------------!
        ut = GetType_common(this%dttype)
    END FUNCTION GetDustType


    PURE FUNCTION GetDustName(this) RESULT(un)
        IMPLICIT NONE
        !----------------------------------------------------------------------!
        TYPE(Dust_TYP), INTENT(IN) :: this
        CHARACTER(LEN=32) :: un
        !----------------------------------------------------------------------!
        un = GetName_common(this%dttype)
    END FUNCTION GetDustName

    PURE FUNCTION DustInitialized(this) RESULT(i)
        IMPLICIT NONE
          !--------------------------------------------------------------------!
          TYPE(Dust_TYP), INTENT(IN) :: this
          LOGICAL :: i
          !--------------------------------------------------------------------!
          i = Initialized_common(this%dttype)
    END FUNCTION DustInitialized


End Module Dust_type
