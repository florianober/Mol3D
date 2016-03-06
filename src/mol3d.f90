!------------------------------------------------------------------------------!
! Mol3D: (SELF-CONSISTENT) N-LTE LINE AND DUST CONTINUUM RADIATIVE TRANSFER CODE
!
! MC3D: RADIATIVE TRANSFER CODE (Mol3D consists of some parts)
!       [by S. Wolf - wolf@astrophysik.uni-kiel.de]
!
! Feedback to:    Florian Ober - fober@astrophysik.uni-kiel.de 
!
!
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
!
! ---
PROGRAM Mol3D
    USE var_global
    USE datatype

    USE basic_type
    USE randgen_type
    USE fluxes_type
    USE model_type
    USE grid_type
    USE gas_type
    USE dust_type
    USE source_type

    USE initiate, ONLY       : inimol, cleanup
    USE fileio, ONLY         : vis_plane, save_ch_map
    USE grd_mod, ONLY        : make_grid
    USE simulation_mod, ONLY : run_simu
    USE math_mod

    IMPLICIT NONE

    !--------------------------------------------------------------------------!
    TYPE(Basic_TYP)        :: basics
    TYPE(Fluxes_TYP)       :: fluxes
    TYPE(Grid_TYP)         :: grid
    TYPE(Model_TYP)        :: model
    TYPE(Dust_TYP)         :: dust
    TYPE(Gas_TYP)          :: gas
    TYPE(SOURCES)          :: sources_in
    REAL                   :: t0,t1, second
    INTEGER                :: hour, minute

    !$ double precision omp_get_wtime
    !$ double precision wt0,wt1   

    !--------------------------------------------------------------------------!
    ! ---
    CALL cpu_time(t0)
    !$ wt0=omp_get_wtime()

    mol3d_version = TRIM(GITHASH)

    
    PRINT '(A)', " ======================================================"
    PRINT '(A)', " Mol3D                                         Welcome!"
    PRINT '(A)', ""
    PRINT '(A)', " git commit: "
    PRINT '(A)', " "//TRIM(mol3d_version)
    PRINT '(A)', " ======================================================"
    
    !--------------------------------------------------------------------------!
    ! 1. initiate code 
    !--------------------------------------------------------------------------!
        CALL inimol(basics,fluxes, grid, model, dust, gas, sources_in)
    !--------------------------------------------------------------------------!
    ! 2. all initialisation done, now set up the grid!
    !--------------------------------------------------------------------------!
        CALL make_grid(basics, grid, model, dust, gas)
    !--------------------------------------------------------------------------!
    ! 3. run simulation
    !--------------------------------------------------------------------------!
        CALL run_simu(basics, fluxes, grid, model, dust, gas, sources_in)
    !--------------------------------------------------------------------------!
    ! 4. clean up
    !--------------------------------------------------------------------------!
        CALL cleanup(basics, fluxes, grid, model, dust, gas, sources_in)
    !--------------------------------------------------------------------------!
    ! 5. time measurement
    !--------------------------------------------------------------------------!
        CALL cpu_time(t1)
        hour   = int((t1-t0)/REAL(3600))
        minute = int((t1-t0)/REAL(60)) - hour*60
        second = t1-t0 - hour*3600 - minute*60
        !$ wt1=omp_get_wtime()
        PRINT *,'|'
        PRINT '(3a,i3,a,i2,a,f5.2,a,f10.2,a)', " simulation ",                 &
                    Getproname(basics), " took",                               &
                    hour, 'h ', minute, 'm ', second, 's    (', t1-t0 ,'sec)'
        !$ write (*,'(a,1pg12.4)') ' omp_get_wtime:', wt1-wt0
    !--------------------------------------------------------------------------!
    
    PRINT *,''
    PRINT *, 'Bye bye'
    PRINT *

END PROGRAM Mol3D
