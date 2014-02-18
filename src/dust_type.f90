!----------------------------------------------------------------------------!
! def Dust_TYP
! inspired by fosite by T. Illenseer 2011
!----------------------------------------------------------------------------!
MODULE Dust_type
  
    USE datatype
    USE gas_type
    USE basic_type
    USE model_type
    USE var_globalnew
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
        TYPE(Common_TYP) :: dttype                       ! -----------------    !
        !-----------------------------------------------------------------------!
        INTEGER                                            :: n_dust
        INTEGER                                            :: n_lam
        INTEGER                                            :: n_lam_map
        INTEGER                                            :: aniiso      
        INTEGER                                            :: n_scatt_th      
        INTEGER                                            :: nrndpt
        INTEGER, DIMENSION(:),POINTER                    :: num_lam_map
        INTEGER, DIMENSION(:), POINTER                   :: cont_map
        INTEGER, DIMENSION(:,:,:),POINTER                :: SCAANG
        
        REAL(kind=r2)                                      :: sizexp      
        REAL(kind=r2)                                      :: D_ANG      
        REAL(kind=r2), DIMENSION(:),POINTER                :: den_dust
        REAL(kind=r2), DIMENSION(:),POINTER                :: r_dust
        REAL(kind=r2), DIMENSION(:),POINTER                :: sidi
        REAL(kind=r2), DIMENSION(:,:),POINTER              :: sidi_par
        
        REAL(kind=r2), DIMENSION(:),POINTER                :: lam
        REAL(kind=r2), DIMENSION(:),POINTER                :: d_lam
        REAL(kind=r2), DIMENSION(:,:),POINTER              :: gscatt_all
        REAL(kind=r2), DIMENSION(:,:),POINTER              :: REF_RE
        REAL(kind=r2), DIMENSION(:,:),POINTER              :: REF_IM
        REAL(kind=r2), DIMENSION(:,:),POINTER              :: sizepar
        REAL(kind=r2), DIMENSION(:,:),POINTER              :: C_sca
        REAL(kind=r2), DIMENSION(:,:),POINTER              :: C_ext
        REAL(kind=r2), DIMENSION(:,:),POINTER              :: C_back
        REAL(kind=r2), DIMENSION(:,:),POINTER              :: C_abs
        REAL(kind=r2), DIMENSION(:,:),POINTER              :: Q_abs
        REAL(kind=r2), DIMENSION(:,:),POINTER              :: albedo
        REAL(kind=r2), DIMENSION(:,:),POINTER              :: QB
        REAL(kind=r2), DIMENSION(:,:),POINTER              :: planck_tab
        REAL(kind=r2), DIMENSION(:,:,:,:,:),POINTER        :: SME
        REAL(kind=r2), DIMENSION(:),POINTER                :: sctth
        REAL(kind=r2), DIMENSION(:),POINTER                :: i_star_emi
        REAL(kind=r2), DIMENSION(:),POINTER                :: c_in_star
        REAL(kind=r2), DIMENSION(:),POINTER                :: c_in_dust
        REAL(kind=r1), DIMENSION(:),POINTER                :: tem_tab
        
        CHARACTER(len=8), DIMENSION(:),POINTER             :: dust_cat
        CHARACTER(len=12), DIMENSION(:,:),POINTER          :: dust_single
        
        LOGICAL, DIMENSION(:),POINTER            :: doRT
        
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
        DustInitialized
    !--------------------------------------------------------------------------!
CONTAINS

    SUBROUTINE InitDust(this,basics , gas, model, ut,un,ndust,dendust, sizexp, aniso, dust_cat,n_scatt_th, nrndpt)
        
        USE math_mod
        IMPLICIT NONE
        !------------------------------------------------------------------------!
        TYPE(Basic_TYP)  , INTENT(IN)        :: basics
        TYPE(Dust_TYP)                        :: this
        TYPE(Gas_TYP)                         :: gas
        TYPE(Model_TYP), INTENT(IN)          :: model
        
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
        REAL(kind=r2), dimension(ndust)       :: dendust
        
        CHARACTER(len=8), DIMENSION(ndust)    :: dust_cat
        CHARACTER(LEN=*)                       :: un
    
        !------------------------------------------------------------------------!
        INTENT(IN)          :: ut,un, ndust, dendust, sizexp, aniso, dust_cat, n_scatt_th, nrndpt
        INTENT(INOUT)       :: this, gas
        !------------------------------------------------------------------------!
        CALL InitCommon(this%dttype,ut,un)
        
        this%n_dust = ndust
        ALLOCATE( this%dust_cat( 1:this%n_dust) )
        this%dust_cat = dust_cat       
        open(unit=1, file=path_dust_cat//this%dust_cat(1)//".dct", &
            action="read", status="unknown", form="formatted")
            read(unit=1,fmt=*) this%n_lam
        close(unit=1)
        
        this%n_lam_map          = this%n_lam
        this%aniiso              = aniso
        this%sizexp             = sizexp
        this%n_scatt_th         = n_scatt_th
        this%nrndpt             = nrndpt
        
        ALLOCATE( this%den_dust( 1:this%n_dust), &
            this%r_dust( 1:this%n_dust), &
            this%sidi( 1:this%n_dust), &
            this%sidi_par( 1:this%n_dust,1:3), &
            this%lam( 1:this%n_lam ), & 
            this%d_lam( 1:this%n_lam ), & 
            this%num_lam_map( 1:this%n_lam_map ), & 
            this%dust_single( 1:this%n_dust, 1:this%n_lam ), &
            this%gscatt_all(  1:this%n_dust, 1:this%n_lam ), &
            this%REF_RE(      1:this%n_dust, 1:this%n_lam ), &
            this%REF_IM(      1:this%n_dust, 1:this%n_lam ), &
            this%sizepar(     1:this%n_dust, 1:this%n_lam ), &
            this%C_sca(       1:this%n_dust, 1:this%n_lam ), &
            this%C_ext(       1:this%n_dust, 1:this%n_lam ), &
            this%C_back(      1:this%n_dust, 1:this%n_lam ), &
            this%C_abs(       1:this%n_dust, 1:this%n_lam ), &
            this%Q_abs(       1:this%n_dust, 1:this%n_lam ), &
            this%albedo(      1:this%n_dust, 1:this%n_lam ), &
            this%QB(      1:this%n_dust, -1:basics%n_tem ), &
            this%planck_tab(  0:basics%n_tem , 1:this%n_lam), &
            this%doRT(       1:this%n_lam ), &
            this%tem_tab(      0:basics%n_tem ), &
            this%i_star_emi( 1:this%n_lam ), &
            this%c_in_star(  1:this%n_lam ), &
            this%c_in_dust(  1:this%n_lam ), &
            this%SME(1:4, 1:4, 1:this%n_dust, this%n_lam, this%n_scatt_th), &
            this%SCAANG(1:this%n_dust,1:this%n_lam,0:this%nrndpt), &
            this%cont_map(1:gas%n_tr), &
            this%sctth(1:n_scatt_th) ) 
        
        
        this%doRT(:)            = .false.
        this%den_dust = dendust
        
        this%r_dust(:)          = 0.0_r2
        this%sidi(:)            = 0.0_r2
        this%sidi_par(:,:)      = 0.0_r2

        
        this%lam(:)             = 0.0_r2
        this%dust_single( :,: ) = ''
        this%gscatt_all(  :,: ) = 0.0_r2
        this%REF_RE(      :,: ) = 0.0_r2
        this%REF_IM(      :,: ) = 0.0_r2
        this%sizepar(     :,: ) = 0.0_r2
        this%C_sca(       :,: ) = 0.0_r2
        this%C_ext(       :,: ) = 0.0_r2
        this%C_back(      :,: ) = 0.0_r2
        this%C_abs(       :,: ) = 0.0_r2
        this%Q_abs(       :,: ) = 0.0_r2
        this%albedo(      :,: ) = 0.0_r2
        this%SME(    :,:,:,:,:) = 0.0_r2
        this%SCAANG(     :,:,:) = 0
        this%sctth(:)           = 0.0_r2
        this%QB(:,:)            = 0.0_r2
        this%tem_tab(:)         = 0.0_r2
        
        this%i_star_emi(:)      = 0.0_r2
        this%c_in_star(:)       = 0.0_r2
        this%c_in_dust(:)       = 0.0_r2
        
        
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
        do i_dust=1, this%n_dust
            if (this%n_dust>1) then
                print *,"   > type of dust [=f(chemistry,size)] :", i_dust, " / ", this%n_dust
            end if
       
            do i_lam=1, this%n_lam
                ! ---
                ! 2.1. load data from file
                open(unit=1, file=path_dust_single//this%dust_single(i_dust,i_lam)//".lot", &
                action="read", status="unknown", form="formatted")
          
                ! handling of anisotropy: aniso == 1:  * anisotropic       scattering *
                !                         aniso == 2:  * isotropic         scattering *
                !                         aniso == 3:  * henyey-greenstein scattering *
                read(unit=1,fmt=*) hi1      ! scatt. matrix / henyey-greenstein
                if (this%aniiso==1) then     ! * before: only anisotropy chosen
                    this%aniiso = hi1        ! * now specify handling of anisotropy (aniiso == 1 or aniiso == 3)
                end if
          
                ! wavelength
                hd1 = this%lam(i_lam)
                read(unit=1,fmt=*)  this%lam(i_lam)
                this%lam(i_lam) = 1.0e-06_r2 * this%lam(i_lam)
          
                ! verify: wavelength correct?
                if ((this%n_dust>1) .and. (i_dust>1)) then
                    if (((hd1-this%lam(i_lam))/hd1) > 1.0e-3) then
                        print *,"!!! Warning: subroutine ld_dust():"
                        print *,"    Wavelengths of different kinds of dust do not agree with each other"
                        !call stop_mc3d()  !TBD!~!
                    end if
                end if
          
                ! grain size distribution
                read(unit=1,fmt=*)  this%sidi_par(i_dust,1)               ! smallest radius [micron]
                read(unit=1,fmt=*)  this%sidi_par(i_dust,2)               ! largest radius  [micron]
                read(unit=1,fmt=*)  this%sidi_par(i_dust,3)               ! size distribution exp. (<0)
                this%sidi_par(i_dust,1:2) = this%sidi_par(i_dust,1:2) * 1.0e-6_r2 ! [micron] -> [m]

                ! HG scattering g
                read(unit=1,fmt=*)  this%gscatt_all(i_dust,i_lam)
          
                ! environment
                read(unit=1,fmt=*)  REF_MED

                ! distribution of scattering angles
                read(unit=1,fmt=*)  this%D_ANG
                read(unit=1,fmt=*)  N_AN

                ! dust grain: refractive index
                read(unit=1,fmt=*)  this%REF_RE( i_dust,i_lam )
                read(unit=1,fmt=*)  this%REF_IM( i_dust,i_lam )

                ! radius & size parameter
                ! r_dust(): used to calculate Q_abs from C_abs and thus in the process of temperature
                !           calculation. however, it effectively cancels out from the equation
                !           for the temperature calculation in sr immediate_temp. therefore,
                !           the value of r_dust() has *no* impact on anything 
                !           (i.e., it should be removed entirely from the file and mc3d source code)
                !           note: in mc3d.v3, dust parameter calculation routines r_dust=maximum grain size
                !                 of a grain size distribution (it should be the surface averaged
                !                 radius of the distribution). however -again-, it has no impact
                !                 on the calculations, whats'o'ever
                read(unit=1,fmt=*)  this%r_dust(i_dust)

                this%r_dust(i_dust) = this%r_dust(i_dust) * 1.0e-6_r2  ! [micron] -> [m]
                if (i_lam == 1) then
                    hd2 = this%r_dust(i_dust)
                else
                    if (this%r_dust(i_dust) /= hd2) then
                        print *,"subroutine ld_dust(): wrong grain size" 
                        print *,"              i_lam : ", i_lam
                        print *,"              i_dust: ", i_dust
                        print *,"            required: ", hd2
                        print *,"           available: ", this%r_dust(i_dust)
                        !call stop_mc3d() TBD!~!
                    end if
                end if
                read(unit=1,fmt=*)  this%sizepar(i_dust,i_lam)

                ! scattering cross sections
                read(unit=1,fmt=*)  this%C_sca(  i_dust,i_lam)
                read(unit=1,fmt=*)  this%C_ext(  i_dust,i_lam)
                read(unit=1,fmt=*)  this%C_back( i_dust,i_lam)

                this%C_abs( i_dust,i_lam) = this%C_ext(i_dust,i_lam) - this%C_sca(i_dust,i_lam)
                this%Q_abs( i_dust,i_lam) = this%C_abs(i_dust,i_lam) / (PI * this%r_dust(i_dust)**2)
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
                        !call stop_mc3d() !TBD!~!
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


!          case (2) ! --- isotropic scattering  ---
!             ! SCAANG: not defined


!          case (3) ! --- anisotropic scattering: HENYEY GREENSTEIN approximation ---
!             ! assumption: streuwinkel = 0, 1, 2, ..., 180 degree
!             ! resulting <g> agrees with the g scattering parameter used as input parameter
!             if (N_AN/=181) then
!                print *,"!!! Warning: subroutine ld_dust()"
!                print *,"             mie-scattering prepared for 181 scattering angles only"
!                call stop_mc3d()
!             end if

!             phg(:)  = 0.0_r2
!             Fphg(:) = 0.0_r2
             
!             do i_scatt_th = 0, 180
!                hd_th = grad2rad(real(i_scatt_th,kind=r2))
!                phg(i_scatt_th) = (1.0_r2 - gscatt_all(i_dust,i_lam)**2) / &
!                     ( 1.0_r2   &
!                     + gscatt_all(i_dust,i_lam)**2   &
!                     - (2.0_r2*gscatt_all(i_dust,i_lam) * cos(hd_th))   &
!                     )**1.5
!                Fphg(i_scatt_th) = Fphg(i_scatt_th-1) + phg(i_scatt_th)*sin(hd_th)
!             end do
!             Fphg(:) = Fphg(:) / Fphg(180)
             
!             ptr_1 = 0
!             ptr_2 = 0
!             do i_scatt_th = 0, 180
!                ptr_1 = int( Fphg(i_scatt_th-1) * real(nrndpt,kind=r2) )
                
!                if (i_scatt_th<180) then
!                   ptr_2 = int( Fphg(i_scatt_th) * real(nrndpt,kind=r2) )
!                   if (ptr_2>nrndpt) then
!                      ptr_2 = nrndpt
!                   end if
!                else
!                   ptr_2 = nrndpt
!                end if
                
!                SCAANG(i_dust,i_lam, ptr_1:ptr_2) = i_scatt_th
!             end do
             
                case default
                    print *,"<!> subroutine ld_dust: select case(aniiso): wrong input parameter. stopped."
                    !call stop_mc3d()
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
    
        ! some output
        if (this%n_dust>1) then
            do i_dust=1, this%n_dust
                print *,"   - type of dust ", i_dust, " -- radius : ", &
                real(this%r_dust(i_dust)*1.0e+6_r2), " micron"
                print *,"     rel. number of particles [%] : ", &
                real(this%sidi(i_dust)*100.0_r2)
            end do
        end if
        
        
        do i_lam_map=1, this%n_lam_map
            this%num_lam_map(i_lam_map) = i_lam_map
        end do
        
        ! select the suited dust wavelength for each tr(ansition) for the selected line transition
        DO tr = 1, gas%n_tr
            hi1 = binary_search(con_c/gas%trans_freq(gas%tr_cat(tr)),this%lam(:))
            
            IF (this%lam(hi1+1) .gt. 2*con_c/gas%trans_freq(gas%tr_cat(tr)) - this%lam(hi1) ) THEN
                this%cont_map(tr) = hi1
            ELSE
                this%cont_map(tr) = hi1+1
            END IF
        END DO
!~         print *, this%lam(this%cont_map(1)) *1e6
!~         print *, this%cont_map(1)
!~         stop
        !print *, gas%tr_cat(1)
        !print *, con_c/gas%trans_freq(gas%tr_cat(1))*1e6
!~         stop
        DO i_lam = 1, this%n_lam
            this%i_star_emi(i_lam) = PI * 4.0_r2 * PI * model%r_star**2 * planck(model%t_star,this%lam(i_lam))  ! [W/m]
            this%c_in_star( i_lam) = this%i_star_emi(i_lam) / model%n_star_emi
        END DO
        
!~         !$ this%c_in_star(:) =  this%c_in_star(:) * basics%num_core
        CALL acc_select(this)
        CALL prepare_QB(this,basics)
    
        
    END SUBROUTINE InitDust
    
    
    SUBROUTINE acc_select(this) !routine by mc3d v4 

        
        IMPLICIT NONE
        
        !------------------------------------------------------------------------!
        TYPE(Dust_TYP)                        :: this
!~         TYPE(Gas_TYP)                         :: gas
        !------------------------------------------------------------------------!
        INTEGER                             :: i_lam, hi3
        REAL(kind=r2)                       :: hd1, hd2, hd3
        LOGICAL                             :: hl1
        
        !------------------------------------------------------------------------!
        ! ---
        ! doRT = .true. -> do RT for this wavelength
        this%doRT(:) = .false.! [default]
           
        ! ---
        ! mark wavelengths
        ! 1. wavelengths for which flux > minimal flux
        hd1 = acc_select_level * maxval(this%lam(:)*this%i_star_emi(:))    !To do set acc_select_level, this%i_star_emi
        do i_lam=1,this%n_lam
            if (this%lam(i_lam)*this%i_star_emi(i_lam) > hd1) then
                this%doRT(i_lam) = .true.
            end if
        end do
           
        ! 2. wavelength to the left and the right beside the found wavelengths
        hl1 = .false.
        do i_lam=1, this%n_lam
            if (this%doRT(i_lam)) then
                hl1 = .true.
            else
                if (hl1) then
                    this%doRT(i_lam) = .true.
                end if
                if (i_lam+1 <= this%n_lam) then
                    if (this%doRT(i_lam+1)) then
                       this%doRT(i_lam) = .true.
                    end if
                end if
                hl1 = .false.
            endif
        end do
           
        ! hd2 ... neglected amount of energy
        ! hd3 ... total amount of energy (including hd2)
        hd3 = integ1( this%lam(:), this%i_star_emi(:), 1, this%n_lam )
        do i_lam=1, this%n_lam
            if (.not. this%doRT(i_lam)) then
                this%i_star_emi(i_lam) = 0.0_r2
            end if
        end do
        hd2 = hd3 - integ1(this%lam(:), this%i_star_emi(:), 1, this%n_lam)
           
        ! count number of (re-emitting cells)*wavelengths
        hi3 = 0
        do i_lam=1,this%n_lam
            if (this%doRT(i_lam)) then
                hi3 = hi3 + 1
            end if
        end do
           
        print *,"radiative transfer only for fraction of wavelengths [%]: ", &
                nint(100.0 *real(hi3)/real(this%n_lam))
        print *,"                      => Neglected amount of energy [%]: ", &
                real(hd2/hd3*100.0_r2)
    
    END SUBROUTINE 
    
    SUBROUTINE prepare_QB(this,basics)

        IMPLICIT NONE
        
        !------------------------------------------------------------------------!
        TYPE(Dust_TYP)                                 :: this
        TYPE(Basic_TYP),INTENT(IN)                    :: basics
        !------------------------------------------------------------------------!
        INTEGER                                        :: i_lam, i_dust, i_tem
        REAL(kind=r2), DIMENSION(:), allocatable    :: QBx
        !------------------------------------------------------------------------
        
        ! ---
        print *,"tabulating: - right side of energy equation"
        print *,"            - planck function               [this may take a while]"

        allocate( QBx(1:this%n_lam) )

        this%QB(:,-1) = -1.0_r2
        
        do i_dust=1, this%n_dust
           do i_tem=0, basics%n_tem
              ! temperature (tabulate)
              ! tbd: logarithmic distribution
              this%tem_tab(i_tem) = real(i_tem, kind=r2) * basics%d_tem  +  basics%t_dust_min
              
              do i_lam=1, this%n_lam
                ! tabulating planck function
                this%planck_tab(i_tem,i_lam) = planck( this%tem_tab(i_tem), this%lam(i_lam))
                 
                ! individual contributions Q_abs * B(T)
                QBx(i_lam) =  this%Q_abs(i_dust,i_lam) * this%planck_tab(i_tem,i_lam)
              enddo

              ! integral [W * m^-2]
              this%QB(i_dust,i_tem) = integ1(this%lam(:), QBx(:), 1, this%n_lam)
           end do
        end do
        
        deallocate(QBx)
        
        ! save QB
!~         print *,"saving array <QB> to ", path_results//pronam//".QB"
!~         open(unit=1, file=path_results//pronam//".QB", &
!~              action="write", status="unknown", form="formatted")
!~         do i_dust=1, n_dust
!~            do i_tem=1, n_tem
!~               write(unit=1,fmt=*) QB(i_dust,i_tem)
!~            end do
!~         end do
!~         close(unit=1)

        ! set wavelength intervals
        this%d_lam(1) = this%lam(2)-this%lam(1)
        DO i_lam=2, this%n_lam-1
           this%d_lam(i_lam) = 0.5_r2*(this%lam(i_lam)-this%lam(i_lam-1)) + &
                               0.5_r2*(this%lam(i_lam+1)-this%lam(i_lam)) 
        END DO
        this%d_lam(this%n_lam) = this%lam(this%n_lam) - this%lam(this%n_lam-1)
        
    END SUBROUTINE prepare_QB


    SUBROUTINE CloseDust(this)
        IMPLICIT NONE
        !------------------------------------------------------------------------!
        TYPE(Dust_TYP), INTENT(INOUT) :: this
        !------------------------------------------------------------------------!
        CALL CloseCommon(this%dttype)
        DEALLOCATE( this%dust_cat, &
                    this%den_dust, &
                    this%r_dust, &
                    this%sidi, &
                    this%sidi_par, &
                    this%lam, & 
                    this%d_lam, & 
                    this%dust_single, &
                    this%cont_map, &
                    this%gscatt_all, &
                    this%REF_RE, &
                    this%REF_IM, &
                    this%sizepar, &
                    this%C_sca, &
                    this%C_ext, &
                    this%C_back, &
                    this%C_abs, &
                    this%Q_abs, &
                    this%albedo, &
                    this%planck_tab, &
                    this%QB, &
                    this%doRT, &
                    this%i_star_emi, &
                    this%c_in_star, &
                    this%c_in_dust, &
                    this%SME, &
                    this%SCAANG, &
                    this%tem_tab, &
                    this%sctth)
        
    END SUBROUTINE CloseDust


    PURE FUNCTION GetDustType(this) RESULT(ut)
        IMPLICIT NONE
        !------------------------------------------------------------------------!
        TYPE(Dust_TYP), INTENT(IN) :: this
        INTEGER :: ut
        !------------------------------------------------------------------------!
        ut = GetType_common(this%dttype)
    END FUNCTION GetDustType


    PURE FUNCTION GetDustName(this) RESULT(un)
        IMPLICIT NONE
        !------------------------------------------------------------------------!
        TYPE(Dust_TYP), INTENT(IN) :: this
        CHARACTER(LEN=32) :: un
        !------------------------------------------------------------------------!
        un = GetName_common(this%dttype)
    END FUNCTION GetDustName

    PURE FUNCTION DustInitialized(this) RESULT(i)
        IMPLICIT NONE
          !------------------------------------------------------------------------!
          TYPE(Dust_TYP), INTENT(IN) :: this
          LOGICAL :: i
          !------------------------------------------------------------------------!
          i = Initialized_common(this%dttype)
    END FUNCTION DustInitialized
    


End Module Dust_type
