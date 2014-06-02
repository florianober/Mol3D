! ---
! mathematical routines
! ---
module math_mod

    USE datatype
    USE var_globalnew

    implicit none
    public :: planck, planckhz, integ1, rad2grad, grad2rad, atan3, atanx, norm,cy2ca, sp2ca, ca2sp, ipol1, ipol2, &
       arcdis, dasico, gauss_arr, solvestateq, binary_search, get_expo,get_opa
    
contains

  ! ################################################################################################
  ! random number generator
!  ! ---
!  subroutine RAN2()
!    use var_globalnew
    
!    implicit none
    
!    integer :: J
!     ---
!    RM=1.0/M
!    if  (IDUM<0 .or. IFF==0)  then
!       IFF  = 1
!       IDUM = modulo(IC-IDUM,M)
!       do J=1,97
!          IDUM = modulo(IA*IDUM+IC,M)
!          IR(J)= IDUM
!       end do
!       IDUM = modulo(IA*IDUM+IC,M)
!       IY   = IDUM
!    endif
    
!    J = 1+(97*IY)/M
!    if (J > 97 .or. J < 1) then
!      print *, "subroutine RAN2(): failed."
!      stop
!    end if
!    IY   = IR( J)
!    rndx = IY*RM
!    IDUM = modulo(IA*IDUM+IC,M)
!    IR(J)= IDUM    
!  end subroutine RAN2


  ! ################################################################################################
  ! kirchhoff-planck function [B.lambda.(T)]
  ! (unsoeld, s. 111, [4.61])
  ! ...
  ! tem_in      ... temperature   [K]
  ! lam_in      ... wavelength    [m]
  ! planck      ... kirchhoff-planck function [W * m^-2 * m^-1]
  ! ---
    PURE FUNCTION planck(tem_in, lam_in) result(planck_result)
  
        !------------------------------------------------------------------------!
        real(kind=r1), intent(in)      :: tem_in
        real(kind=r2), intent(in)      :: lam_in
        real(kind=r2)                   :: planck_result
        real(kind=r2)                   :: const
        !------------------------------------------------------------------------!
        
        if ( tem_in .gt. 1.0e-30 ) then
            
            const = 2.0*con_h*con_c**2.0/(lam_in**5.0)
            
            !Use Wien approximation for h*c/lam>>k*T:
            if (con_h*con_c/lam_in .gt. 100.0*con_k*tem_in) then
                planck_result = const * exp(-con_c2/(tem_in*lam_in))
            
!~ !            elseif ( con_h*con_c/lam_in .lt. 0.001_r2*con_k*tem_in) then
!~ !                planck_result = 2.0_r2*con_c*con_k*tem_in/lam_in**4.0_r2
!~ !                print *,'here'
            else
                planck_result = const / (exp(con_c2/(tem_in*lam_in))-1.0)
            endif
            
        else
            planck_result = 0.0_r2
        endif
        
    END FUNCTION planck

    ELEMENTAL FUNCTION get_expo(var_in) result(var_out)
        !------------------------------------------------------------------------!
        real(kind=r2), intent(in)      :: var_in
        real(kind=r2)                  :: var_out
        !------------------------------------------------------------------------!
        IF (var_in .lt. -11.0) THEN
            var_out = 0.0_r2
        ELSE
            var_out = exp(var_in)
        END IF
    END FUNCTION 

    ELEMENTAL FUNCTION get_opa(var1,var2,var3) result(var_out)
        !------------------------------------------------------------------------!
        real(kind=r2), intent(in)      :: var1,var2,var3
        real(kind=r2)                  :: var_out
        !------------------------------------------------------------------------!
        
        IF (var2 .lt. 1.0e-100_r2) THEN
            var_out = var3
        ELSE
            var_out = var1*var2+var3
        END IF
    END FUNCTION 
    
! ################################################################################################
  ! kirchhoff-planck function [B.lambda.(T)]
  ! (unsoeld, s. 111, [4.61])
  ! ...
  ! tem_in      ... temperature    [K]
  ! freq_in      ... wavelength    [Hz]
  ! planck      ... kirchhoff-planck function [W * m^-2 * Hz^-1]
  ! ---
  
    ELEMENTAL FUNCTION planckhz(tem_in, freq_in) result(planck_result)
  
        !------------------------------------------------------------------------!
        real(kind=r1), intent(in)       :: tem_in
        real(kind=r1), intent(in)       :: freq_in
        real(kind=r2)                   :: planck_result
        real(kind=r2)                   :: const
        real(kind=r2)                   :: hnu_tk
        !------------------------------------------------------------------------!
        
        if ( tem_in .gt. 1.0e-30 ) then
            
!~             const = 2.0*con_h*freq_in**3/(con_c**2)
            const = 2.0*con_h*freq_in * (freq_in/con_c) *(freq_in/con_c)
            !Use Wien approximation for h*c/lam>>k*T:
            hnu_tk = con_h*freq_in/(tem_in*con_k)
!~             print *,const, hnu_tk
            if (hnu_tk .gt. 100.0) then
                planck_result = const * exp(-hnu_tk)
!~                 print *,'lala'
            
            else
                planck_result = const / (exp(hnu_tk)-1.0)
!~                 print *, 'ich'
            endif
        else
            planck_result = 0.0
        endif
        !print *, const * exp(-hnu_tk) ,con_c /freq_in
    END FUNCTION planckhz

  ! ################################################################################################
  ! "linear" integration
  ! ---
  function integ1( x, y, xlow, xup ) result(integ1_result)
    use datatype

    real(kind=r2), dimension(:), intent(in) :: x, y
    integer, intent(in)                     :: xlow, xup
    real(kind=r2)                           :: integ1_result
    integer :: i1
    ! ---
    if (xlow /= xup) then
       integ1_result = 0.0_r2
       do i1=xlow+1, xup
          integ1_result = integ1_result &
               + (x(i1)-x(i1-1)) * y(i1-1) &
               + (x(i1)-x(i1-1)) * (y(i1)-y(i1-1)) / 2.0_r2
       enddo
    else
       integ1_result = 0.0_r2
    end if
  end function integ1


!!$  ! ################################################################################################
!!$  ! "linear" integration
!!$  ! ---
!!$  function integ1( x, y, xlow, xup ) result(integ1_result)
!!$    use datatype
!!$
!!$    real(kind=r2), dimension(:), intent(in) :: x, y
!!$    integer, intent(in)                     :: xlow, xup
!!$    real(kind=r2)                           :: integ1_result
!!$    integer :: i1
!!$    ! ---
!!$    if (xlow /= xup) then
!!$       integ1_result = 0.0_r2
!!$
!!$       do i1=xlow+1, xup
!!$          integ1_result = integ1_result + (x(i1)-x(i1-1)) * (y(i1)-y(i1-1))
!!$       enddo
!!$       integ1_result = integ1_result / 2.0_r2
!!$
!!$       do i1=xlow+1, xup
!!$          integ1_result = integ1_result + (x(i1)-x(i1-1)) * y(i1-1)
!!$       enddo
!!$    else
!!$       integ1_result = 0.0_r2
!!$    end if
!!$  end function integ1


  ! ################################################################################################
  ! rad -> grad
  ! ---
  function rad2grad(rad) result(rad2grad_result)
    use datatype
    use var_globalnew

    implicit none
    real(kind=r2), intent(in) :: rad
    real(kind=r2)             :: rad2grad_result
    ! ---
    rad2grad_result = rad * 180.0_r2/PI
  end function rad2grad


  ! ################################################################################################
  ! grad -> rad
  ! ---
  function grad2rad(grad) result(grad2rad_result)
    use datatype
    use var_globalnew

    implicit none
    real(kind=r2), intent(in) :: grad
    real(kind=r2)             :: grad2rad_result
    ! ---
    grad2rad_result = grad * PI/180.0_r2
  end function grad2rad


  ! ################################################################################################
  ! atan3(y,x)
  ! - return angle in range 0...2pi [unit: radian]
  !   (y,x) = ( 0, 1) :   0°
  !           ( 1, 0) :  90°
  !           ( 1,-1) : 135°
  !           (-1, 0) : 270°
  ! ---
  PURE function atan3(y,x) result(atan3_result)
    use datatype
    use var_globalnew

    implicit none
    real(kind=r2), intent(in) :: x,y
    real(kind=r2)             :: atan3_result
    real(kind=r2)             :: hd1
    ! ---
    hd1 = atan2(y,x)
    if (hd1 > 0.0_r2) then
       atan3_result = hd1
    else
       atan3_result = hd1 + 2.0_r2*PI
    end if
  end function atan3
  

  ! ################################################################################################
  ! x & y => arcustanges(y/x)=w  [rad]                -3.14 <= w <= 3.14
  !  ---
  function atanx( yy, xx ) result(atanx_result)
    use datatype

    real(kind=r2), intent(in) :: xx, yy
    real(kind=r2) :: atanx_result
    real(kind=r2) :: rr
    ! ---
    rr = sqrt( xx**2 + yy**2 )
    if (rr==0.0_r2) then
       atanx_result = 0.0_r2   ! any value
    else
       atanx_result = sign( acos(xx/rr), yy/rr )
    end if
  end function atanx


    ! ################################################################################################
    ! norm(vector)
    ! ---
    PURE FUNCTION norm(vec) result(norm_result)

        implicit none
        real(kind=r2), dimension(3), intent(in) :: vec
        real(kind=r2) :: norm_result
        ! ---
        norm_result = sqrt( dot_product( vec, vec ) )
        !norm_result = sqrt( vec(1)**2+vec(2)**2+vec(3)**2)
    END function norm
  

  ! ################################################################################################
  ! in : spherical coordinates    spco(r, theta, phi)
  ! out: cartesian coordinates    caco(x, y, z)
  ! --
  ! (r=R,theta= 0°,phi=     0°)  =>  (x= R,y=0,z=0)
  ! (r=R,theta=90°,phi=     0°)  =>  (x= 0,y=0,z=R)
  ! (r=R,theta= 0°,phi=+/-180°)  =>  (x=-R,y=0,z=0)
  ! ---
    PURE FUNCTION sp2ca( spco) RESULT( caco )
    
        IMPLICIT NONE
        
        real(kind=r2), dimension(1:3), intent(in)   :: spco
        real(kind=r2), dimension(1:3)               :: caco 
        ! ---
        caco(1) = spco(1) * cos( spco(2) ) * cos( spco(3) )
        caco(2) = spco(1) * cos( spco(2) ) * sin( spco(3) )  
        caco(3) = spco(1) * sin( spco(2) )
        
    END FUNCTION sp2ca
    
    
    PURE FUNCTION cy2ca( cyco) RESULT( caco )
    
        IMPLICIT NONE
        
        real(kind=r2), dimension(1:3), intent(in)  :: cyco
        real(kind=r2), dimension(1:3)               :: caco 
        ! ---
        caco(1) = cyco(1) * cos( cyco(2) ) 
        caco(2) = cyco(1) * sin( cyco(2) )
        caco(3) = cyco(3)
        
    END FUNCTION cy2ca    
    

    FUNCTION ca2sp( caco ) RESULT (spco)
  
        IMPLICIT NONE
        
        real(kind=r2), dimension(1:3), intent(in)  :: caco
        real(kind=r2), dimension(1:3)              :: spco
        
        spco(1)     = norm(caco(:))
        spco(2)     = atan2( caco(3), sqrt(caco(1)**2 + caco(2)**2))
        spco(3)     = atan3( caco(2), caco(1) )
        
    
    END FUNCTION ca2sp

  ! ################################################################################################
  ! linear interpolation, typ 1: very general
  ! ---
  !  # 1. feldindex != 1
  !  # werte muessen nicht sortiert sein
  !  # sowohl x- als auch y-werte muessen vom typ "r2" sein
  !  - xi_i   : x-wert
  !    xbas_i : vektor, der mindestens 2 x-werte enthaelt
  !    ybas_i : vektor, der die zugehoerigen y-werte enthaelt
  !  - xbas_i od. ybas_i muessen keinerlei ordnung aufweisen 
  !  - min(xbas_i(*)) <= xi_i <= max(xbas_i(*))
  ! ---
  function ipol1( xi_i, xbas_i, ybas_i ) result(ipol1_result)
    use datatype
    
    real(kind=r2), intent(in) :: xi_i
    real(kind=r2), dimension(:), intent(in)  :: xbas_i, ybas_i
    real(kind=r2) :: ipol1_result
    integer :: nx, i1, icount, ib1, ib2
    integer,dimension(1:1) :: sxbas, sybas
    logical :: f1
    real(kind=r2) :: xi, hd1
    real(kind=r2), allocatable, dimension(:) :: xbas, ybas
    ! ---
    sxbas = ubound(xbas_i)
    sybas = ubound(ybas_i)

     !---
     !some tests
    if ( sxbas(1) /= sybas(1) ) then
       print *, "error: ipol1: dimensions of x- and y- achse not equal"
       stop
    endif
    
    if ( sxbas(1) == 1 ) then
       print *, "error: ipol1: dim(xbas) >! 1; here: dim(xbas)=1"
       stop
    endif

    nx = sxbas(1)  ! number of indezes resp. pairs
    allocate (xbas(1:nx), ybas(1:nx))
    xbas(:) = xbas_i(:)
    ybas(:) = ybas_i(:)
    xi      = xi_i

    ! ---
    ! 1. sort xbas/ybas: xbas(i) < xbas(i+1)
    do
       f1 = .false. !  default: no change needed
       do i1=2,nx
          if ( xbas(i1) < xbas(i1-1) ) then
             hd1        = xbas(i1) 
             xbas(i1)   = xbas(i1-1) 
             xbas(i1-1) = hd1
 
             hd1        = ybas(i1)
             ybas(i1)   = ybas(i1-1)
             ybas(i1-1) = hd1
             
             f1 = .true.
          endif
       end do

       if (f1) then
          cycle
       else
          exit
       end if
    end do

    ! ---
    ! test:  x_min <= xi <= x_max
    if ( (xi < xbas(1)) .or. (xi > xbas(nx)) ) then
       print *, "error: ipol1: xi outside definition range"
       stop
    endif

    ! ---
    ! find xbas.min <= xi <= xbas.max
    icount = 0

    do
       icount = icount+1
       if (xi > xbas(icount)) then
          cycle
       else
          exit
       end if
    end do

    ib1 = icount-1 ! left  value
    ib2 = icount   ! right value
    if (ib1 < 1) then 
       ib1=1
    end if

    ! ---
    ! interpolation
    if (ib1 == ib2) then
       ipol1_result = ybas(ib1)
    else
       ipol1_result = &
            (ybas(ib2)-ybas(ib1))*(xi-xbas(ib1))/(xbas(ib2)-xbas(ib1)) &
            + ybas(ib1)
    end if

    deallocate (xbas, ybas)
  end function ipol1


   !################################################################################################
   !linear interpolation, typ 2: as simple as possible...
   !---
  ELEMENTAL function ipol2( x_min, x_max, y_min, y_max, x_ipol ) result(ipol2_result)
    use datatype
    
    real(kind=r2), intent(in) :: x_min, x_max, y_min, y_max, x_ipol
    real(kind=r2) :: ipol2_result
    ! ---
    if (x_min == x_max) then
       ipol2_result = y_min
    else
       ipol2_result = (y_max-y_min)*(x_ipol-x_min)/(x_max-x_min) + y_min
    end if
  end function ipol2


  ! ################################################################################################
  ! angular distance of 2 points with the spherical coordinates
  ! (a1,d1) and (a2,d2) [radian]                           Meeus, p. 118
  ! theta ... angle to z axis
  ! phi   ... angle to x axis
  ! ---
  function arcdis( &
       p1_cos_th, p1_sin_th, p1_ph, &
       p2_cos_th, p2_sin_th, p2_ph &
       ) result(arcdis_result)
    use datatype
  
    implicit none
    real(kind=r2), intent(in) :: &
         p1_cos_th, p1_sin_th, p1_ph, &
         p2_cos_th, p2_sin_th, p2_ph
    real(kind=r2)             :: arcdis_result
    ! ---
    arcdis_result = acos( p1_cos_th*p2_cos_th + p1_sin_th*p2_sin_th*cos(p1_ph-p2_ph) )
  end function arcdis
  

  ! ################################################################################################
  ! DCOS(winkel) [DCOSW] & DSIN(winkel) [DSINW] => W; [W]=rad
  !                                                   -3.14 <= w <= 3.14
  ! ---
  function dasico( dsinw, dcosw ) result(dasico_result)
    use datatype

    implicit none
    real(kind=r2), intent(in) :: dsinw, dcosw
    real(kind=r2)             :: dasico_result
    ! ---
    dasico_result = sign( acos(dcosw), dsinw )
  end function dasico


  ! ################################################################################################
  ! conversion: W/m/ x sr => W/m/sr !TBD!-!
  ! ---
!~     PURE FUNCTION cnv_Wmsr( lum_x, map_x ) result(cnv_Wmsr_result)
!~         use datatype
!~         use var_globalnew
!~         
!~         IMPLICIT NONE
!~         
!~         !--------------------------------------------------------------------------!
!~         real(kind=r2), intent(in) :: lum_x
!~         integer,       intent(in) :: map_x
!~         !--------------------------------------------------------------------------!        
!~         real(kind=r2)             :: cnv_Wmsr_result
!~         !--------------------------------------------------------------------------!
!~         ! ---
!~         ! 1. normalization to W/m/sr
!~         !if ((simu_type==2 .or. simu_type==3) .and. ree_type==2) then
!~             ! raytracing: default: W/m/sr => no conversion needed
!~         cnv_Wmsr_result = lum_x
!~         !else
!~         !    ! MC RT
!~         !    if (n_ph==1 .and. project_2D) then
!~         !        ! rotational symmetry used
!~         !        cnv_Wmsr_result = lum_x / ( PIx2  * ( cos(th_map(map_x)-al_map(map_x)) - &
!~         !        cos(th_map(map_x)+al_map(map_x)) ) )
!~         !    else
!~         !        ! rotational symmetry NOT used
!~         !        cnv_Wmsr_result = lum_x / ( PIx2 * (1.0_r2 - cos(al_map(map_x))) )
!~         !    end if
!~         !end if
!~     END FUNCTION cnv_Wmsr


  ! ################################################################################################
  ! conversion: W/m => Jy
  ! ---
!~     function cnv_lum2Jy( lum_x, lam_x, map_x, distance ) result(cnv_lum2Jy_result)
!~         use datatype
!~         use var_globalnew
!~         
!~         IMPLICIT NONE
!~         
!~         !--------------------------------------------------------------------------!
!~         real(kind=r2), intent(in) :: lum_x
!~         real(kind=r2), intent(in) :: lam_x
!~         real(kind=r2), intent(in) :: distance
!~         
!~         integer,       intent(in) :: map_x
!~         !--------------------------------------------------------------------------!
!~         real(kind=r2)             :: hd_ny
!~         real(kind=r2)             :: hd_flux
!~         real(kind=r2)             :: cnv_lum2Jy_result
!~         !--------------------------------------------------------------------------!      
!~         
!~         ! ---
!~         ! 1. normalization to W/m/sr
!~         hd_flux = cnv_Wmsr( lum_x, map_x )
!~ 
!~         ! 2. frequency
!~         hd_ny = con_c  / lam_x
!~ 
!~         ! 3. [L_ny] = W * Hz^-1 * sr^-1
!~         hd_flux = hd_flux  *  con_c/(hd_ny**2)
!~         
!~         ! 4. if we know the distance to the observer we know how much energy is flowing
!~         !    through a square meter - at his location - each second and Hz. Thus, we arrive at Jansky.
!~         hd_flux = hd_flux  *  1.0e+26_r2 / ((distance*con_pc)**2)
!~ 
!~         cnv_lum2Jy_result = hd_flux
!~     end function cnv_lum2Jy
    
    ELEMENTAL FUNCTION gauss_arr(lscale, a, velo) RESULT (gaussian_result)
    
        !USE gas_type
    
        IMPLICIT NONE
        !--------------------------------------------------------------------------!

        REAL(kind=r2), INTENT(IN)                               :: lscale
        REAL(kind=r2), INTENT(IN)                               :: a
        REAL(kind=r2), INTENT(IN)                               :: velo
        
        REAL(kind=r2)                                            :: gaussian_result
        
        !--------------------------------------------------------------------------!

        gaussian_result = lscale/a * exp( -(velo/a)**2 )
      

    END FUNCTION gauss_arr
    
    PURE SUBROUTINE solvestateq(A,lenA,c,x)
        !
        ! copyright Gisela Engeln-Muellges: Numerik-Algorithmen, P.147
        !
        ! solve the linear equation Ax = c
        IMPLICIT NONE
        !--------------------------------------------------------------------------!
        INTEGER, INTENT(IN)                                       :: lenA   
        INTEGER, DIMENSION(1:lenA)                                :: p
        
        REAL(kind=r1), DIMENSION(1:lenA,1:lenA) , INTENT(INOUT)   :: A
        REAL(kind=r1), DIMENSION(1:lenA) , INTENT(INOUT)          :: c
        REAL(kind=r1), DIMENSION(1:lenA) , INTENT(INOUT)          :: x
     

        
        REAL(kind=r2), DIMENSION(1:lenA)                          :: r_vec
        REAL(kind=r2), DIMENSION(1:lenA)                          :: help

        !REAL(kind=r2)                                              :: t
        
        INTEGER                                                    :: i,j,k, i_0, t
        
        !--------------------------------------------------------------------------!
        DO i = 1 , lenA
            P(i) = i
        END DO
        DO j = 1,lenA-1
        
            ! find i_0 (pivot element)
            i_0 = int(maxloc(abs(A(j:lenA,j)),1)+(j-1))
            
            ! set pivot vector
            t = p(i_0)
            p(i_0) = p(j)
            p(j) = t
            
            ! change lines in A
            help = A(i_0,:)
            A(i_0,:) = A(j,:)
            A(j,:) = help
            !calculate values
            DO i =j+1, lenA
                A(i,j) = A(i,j)/A(j,j)
                DO k = j+1,lenA
                    A(i,k) = A(i,k) - A(j,k)*A(i,j)
                END DO
                
            END DO
        END DO
        ! just renaming here, can be commented out
!~         DO i = 1, lenA
!~             DO j = 1, lenA
!~                 IF (j .lt. i) THEN
!~                     L(i,j) = A(i,j)
!~                     R(i,j) = 0.0_r2
!~                     
!~                 ELSEIF (i == j) THEN
!~                     L(i,j) = 1.0_r2
!~                     R(i,j) = A(i,j)
!~                 ELSE
!~                     L(i,j) = 0.0_r2
!~                     R(i,j) = A(i,j)
!~                 END IF
!~             
!~             END DO
!~         END DO 
        
        ! forward solving Lr = a
        r_vec(1) = c(p(1))
        DO i = 2, lenA
            r_vec(i) = c(p(i)) - sum(A(i,1:i-1)*r_vec(1:i-1))
        END DO
        x(lenA) = r_vec(lenA)/A(lenA,lenA)
        
        ! backward solving Rx = r
        DO i = lenA-1,1,-1
            x(i) = 1.0_r2/A(i,i)*(r_vec(i)-sum(A(i,i+1:lenA)*x(i+1:lenA)))
        END DO
    END SUBROUTINE solvestateq
    
    SUBROUTINE vecmat(simu_var)
        USE simu_type
        IMPLICIT NONE
        
        TYPE(Simu_TYP), INTENT(INOUT)                     :: simu_var
        
        real(kind=r2), dimension(1:3,1:3)                 :: D_help
        REAL(kind=r2)                                      :: R,L,P
        
        ! ---
        R = -simu_var%SINTHE * simu_var%SINPHI
        L =  simu_var%SINTHE * simu_var%COSPHI
        P =  simu_var%COSTHE

        simu_var%dir_xyz(2) = -( simu_var%D(1,1) * R  +  simu_var%D(1,2) * L  +  simu_var%D(1,3) * P )
        simu_var%dir_xyz(1) =    simu_var%D(2,1) * R  +  simu_var%D(2,2) * L  +  simu_var%D(2,3) * P
        simu_var%dir_xyz(3) =    simu_var%D(3,1) * R  +  simu_var%D(3,2) * L  +  simu_var%D(3,3) * P

        D_help(1,1) =   simu_var%COSPHI
        D_help(2,1) =   simu_var%SINPHI
        D_help(3,1) =   0.0_r2
        D_help(1,2) = - simu_var%SINPHI*simu_var%COSTHE
        D_help(2,2) =   simu_var%COSPHI*simu_var%COSTHE
        D_help(3,2) = - simu_var%SINTHE
        D_help(1,3) = - simu_var%SINPHI*simu_var%SINTHE
        D_help(2,3) =   simu_var%COSPHI*simu_var%SINTHE
        D_help(3,3) =   simu_var%COSTHE

        simu_var%D = matmul(simu_var%D,D_help)
    end subroutine vecmat
    
    PURE FUNCTION binary_search(var_in, array_in) RESULT (int_result)
    
        IMPLICIT NONE
        !------------------------------------------------------------------------!
        REAL(kind=r2), DIMENSION(:), INTENT(IN)       :: array_in
        REAL(kind=r2), INTENT(IN)                     :: var_in
        
        INTEGER                                        :: int_result
        INTEGER                                        :: L, R, MID
        !------------------------------------------------------------------------!
        int_result = 0
        L = 1
        R = SIZE(array_in)
        IF ( var_in .gt. array_in(R)) THEN
            int_result = R
        ELSE
            DO WHILE (L .le. R )
                MID = INT((L+R)/2)
                IF (var_in .ge. array_in(MID) .and. var_in .le. array_in(MID+1) ) THEN
                    !print *, var_in,  array_in(MID) ,array_in(MID+1), MID
                    int_result = MID
                    EXIT
                ELSE
                    IF ( var_in < array_in(MID) ) THEN
                        R = MID-1
                    
                    ELSE IF ( var_in > array_in(MID) ) THEN
                        L = MID+1
                    END IF
                END IF
            END DO
        END IF
        
    END FUNCTION binary_search
    
    
end module math_mod
