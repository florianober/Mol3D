module var_globalnew
  use datatype

  implicit none
  PUBLIC
 
  ! ---
  ! parameters
  ! ---
  ! Msun      ... solar mass [kg]
  ! con_gamma ... grav. constant          [m^3 / (kg * s^2)]
  ! con_K     ... Boltzmann constant      [J/K]
  ! con_sigma ... stefan/boltzmann const. [W / (m^2 * K^4)] 
  ! con_mu    ... atomic mass unit        [kg]
  ! con_ c2   ... con_h*con_c/con_k
  ! real(kind=r3), public, parameter :: &
  !     PI        = 3.14159265358979324_r3
  
  
    real(kind=r2),parameter :: &
        PI        = 3.14159265358979324_r2, &
        con_h     = 6.626176e-34_r2, &
        con_k     = 1.380662e-23_r2, &
        con_Na    = 6.02214129e23_r2, &
        con_c     = 299792458.0_r2, &
        con_c2    = 1.43878631423230881E-2_r2, &
        con_AU    = 1.496e+11_r2, &
        con_pc    = 3.0856776e+16_r2,&
        con_gamma = 6.67384e-11_r2, &
        con_sigma = 5.6704e-8_r2, &
        con_mu    = 1.660538e-27_r2, &
        L_sun     = 3.85e+26_r2, &
        M_sun     = 1.9891e+30_r2, &
        R_sun     = 0.6960e+9, &
        SBK       = 5.67e-8_r2, &
        w2erg     = 1.0e+07_r2
    REAL(kind=r2),PARAMETER,DIMENSION(1:5)    :: &
        col_p_weight = (/2.01588_r2,0.0_r2,0.0_r2,0.0_r2,0.0_r2/)
        
    REAL(kind=r2),PARAMETER,DIMENSION(1:6)       ::  &
        RK_c  = (/0.0_r2, 1.0_r2/4.0_r2, 3.0_r2/8.0_r2, 12.0_r2/13.0_r2, 1.0_r2, 0.5_r2/), &
        RK_b1 = (/16.0_r2/135.0_r2, 0.0_r2, 6656.0_r2/12825.0_r2, 28561.0_r2/56430.0_r2, &
                -9.0_r2/50.0_r2, 2.0_r2/55.0_r2/), &
        RK_b2 = (/25.0_r2/216.0_r2, 0.0_r2, 1408.0_r2/2565.0_r2, 2197.0_r2/4104.0_r2, &
                -1.0_r2/5.0_r2, 0.0_r2/)
    
        ! Runge Kutta constants
    REAL(kind=r2),PARAMETER,DIMENSION(1:6,1:6)       ::  &   
        RK_a  = RESHAPE( (/ &
                0.0_r2, 0.0_r2, 0.0_r2, 0.0_r2, 0.0_r2, 0.0_r2, & 
                1.0_r2/4.0_r2, 0.0_r2, 0.0_r2, 0.0_r2, 0.0_r2, 0.0_r2, &
                3.0_r2/32.0_r2, 9.0_r2/32.0_r2, 0.0_r2, 0.0_r2, 0.0_r2, 0.0_r2, &
                1932.0_r2/2197.0_r2, -7200.0_r2/2197.0_r2, 7296.0_r2/2197.0_r2, 0.0_r2, 0.0_r2, 0.0_r2, &
                439.0_r2/216.0_r2, -8.0_r2, 3680.0_r2/513.0_r2, -845.0_r2/4104.0_r2, 0.0_r2, 0.0_r2, &
                -8.0_r2/27.0_r2, 2.0_r2, -3544.0_r2/2565.0_r2, 1859.0_r2/4104.0_r2, -11.0_r2/40.0_r2, &
                0.0_r2/), (/6,6/))
!~         RK_a  = RESHAPE( (/ &
!~                 0.0_r2, 1.0_r2/4.0_r2, 3.0_r2/32.0_r2, 1932.0_r2/2197.0_r2, 439.0_r2/216.0_r2, -8.0_r2/27.0_r2, &
!~                  
!~                 0.0_r2, 0.0_r2, 0.0_r2, 0.0_r2, 0.0_r2, &
!~                 , 9.0_r2/32.0_r2, 0.0_r2, 0.0_r2, 0.0_r2, 0.0_r2, &
!~                  -7200.0_r2/2197.0_r2, 7296.0_r2/2197.0_r2, 0.0_r2, 0.0_r2, 0.0_r2, &
!~                  -8.0_r2, 3680.0_r2/513.0_r2, -845.0_r2/4104.0_r2, 0.0_r2, 0.0_r2, &
!~                  2.0_r2, -3544.0_r2/2565.0_r2, 1859.0_r2/4104.0_r2, -11.0_r2/40.0_r2, &
!~                 0.0_r2/), (/6,6/)))

  !integer, public, parameter :: &
  !     unit_tmp = 4           !????

  character(len=75), public, parameter :: &
       hrule = "---------------------------------------------------------------------------"
  character(len= 5), public, parameter :: path_misc        = "misc/"
  character(len=11), public, parameter :: path_dust_cat    = "input/dust/"
  character(len=18), public, parameter :: path_dust_single = "input/dust/single/"
  character(len=21), public, parameter :: path_dust_nk     = "input/dust/tables_nk/"
  character(len=21), public, parameter :: path_dust_wv     = "input/dust/tables_wv/"
  character(len=15), public, parameter :: path_dust_tmp    = "input/dust/tmp/"
  character(len=12), public, parameter :: path_model       = "input/model/"
  character(len=13), public, parameter :: path_extern      = "input/extern/"
  character(len=10), public, parameter :: path_mol         = "input/mol/"
  
  character(len=8), public, parameter :: mol3d_version    = "Feb 2014"
  


  ! ---
  ! strings
  ! ---
  !###########3character(len=8), public :: pronam_temp
  !character(len=12), public :: dust_single_tmp
  !character(len=12), public, allocatable, dimension(:,:) :: dust_single

  ! ---
  ! integer
  ! ---
  integer, public :: pre_ana, i_lam_show
  !integer, public :: nr_lam
  !integer, public :: lot_th
  integer, public :: kill_photon_count, n_wrong_temp, n_interact_max
  !integer, public, allocatable, dimension(:,:,:) :: SCAANG
  !integer, public :: n_den_par

  ! ---
  ! real(kind=r2)
  ! ---
  !real(kind=r2), public :: i_min
  real(kind=r2), public :: acc_select_level
  !real(kind=r2), public :: SINPHI, COSPHI, SINTHE, COSTHE, SIN2PH, COS2PH
  !real(kind=r2), public :: t_dust_min, t_dust_max
  !real(kind=r2), public :: zoom  

  !real(kind=r2), dimension(1:3), public :: pos_xyz, pos_xyz_li
  !real(kind=r2), dimension(1:3,1:3), public :: D

  !real(kind=r2), allocatable, dimension(:), public :: tem_tab, d_lam, diff_planck

  !real(kind=r2), allocatable, dimension(:), public :: sin_th_map, cos_th_map, sin_ph_map, cos_ph_map
  !real(kind=r2), allocatable, dimension(:), public :: den_par

  !real(kind=r2), allocatable, dimension(:,:), public :: planck_tab


!  real(kind=r2), allocatable, dimension(:,:,:), public :: i_star_abs 

!  real(kind=r2), allocatable, dimension(:,:,:),     public :: stokes_sed
!  real(kind=r2), allocatable, dimension(:,:,:,:,:), public :: stokes_map

!  real(kind=r2), allocatable, dimension(:,:,:,:,:), public :: SME

  ! ---
  ! logical
  ! ---
  logical, public :: show_error
  !logical, public :: inside, show_error, scatt_mat


end module var_globalnew
