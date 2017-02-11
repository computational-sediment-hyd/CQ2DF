! Copyright (C) 2017 computational-sediment-hyd https://github.com/computational-sediment-hyd
! Released under the MIT license
! http://opensource.org/licenses/mit-license.php

!> @mainpage @file Cross-sectionQuasi2DFlow.f90
!! @brief Cross section Quasi 2D Flow
!! @detail F03 format
!! @detail Vesion 1.0.0
!! @date  2017.02.01
!! @date  Last Update 2017.02.01
!! @author 

!> @brief class
!! @detail 
module classModel
    implicit none
    PRIVATE
    TYPE, PUBLIC :: model
        INTEGER, PUBLIC :: ny, nz
        DOUBLE PRECISION, PRIVATE :: ib, ks, kappa, dy, dz
    CONTAINS
        PROCEDURE, PUBLIC :: initilize
        PROCEDURE, PUBLIC :: calSectionProfile
        PROCEDURE, PRIVATE :: getNut
        PROCEDURE, PRIVATE :: getInitialU
        PROCEDURE, PRIVATE :: calSectU
    END TYPE
contains

!> @brief Constructor
!! @detail 
!! @param 
!! @return 
subroutine initilize(self, ib_, ks_, kappa_, ny_, nz_, dy_, dz_)
    class(model) :: self
    INTEGER, INTENT(IN) :: ny_, nz_
    DOUBLE PRECISION, INTENT(IN) :: ib_, ks_, kappa_, dy_, dz_
    self%ib = ib_
    self%ks = ks_
    self%kappa = kappa_
    self%ny = ny_
    self%nz = nz_
    self%dy = dy_
    self%dz = dz_
end subroutine

!> @brief model
!! @detail 
!! @param 
!! @return 
subroutine calSectionProfile(self, cellcnd, un)
    class(model) :: self
    INTEGER, INTENT(IN) :: cellcnd(self%ny, self%nz)
    DOUBLE PRECISION, INTENT(OUT) :: un(self%ny, self%nz)
    DOUBLE PRECISION :: unnew(self%ny, self%nz)
    DOUBLE PRECISION :: nut(self%ny, self%nz)
    DOUBLE PRECISION :: err
    INTEGER :: n,j,k

!set ν_t
    nut = self%getNut(cellcnd)
!set initial condition
    un = self%getInitialU(cellcnd)
!main iteration
    print *, 'Calculate...'
    do n = 1, 100000
        CALL self%calSectU(un, cellcnd, nut)
        err = maxval( abs( unnew(:,:) - un(:,:) ) )
        if(err > 0.000001)then
            unnew(:,:) = un(:,:)
        else
            print *, 'iteration = ' ,n,'error = ' ,err
            exit
        end if
    end do

end subroutine

!> @brief calculate eddy viscosity
!! @detail 
!! @param 
!! @return 
function getNut(self, cellcnd) RESULT(nut)
    class(model) :: self
    INTEGER, INTENT(IN) :: cellcnd(self%ny, self%nz)
    DOUBLE PRECISION :: nut(self%ny, self%nz)
    INTEGER :: j, k , kk
    DOUBLE PRECISION :: uster, depth, zc

    DO j = 1, self%ny
        DO k = 1, self%nz
            IF(mod(cellcnd(j, k), 10) == 1)THEN ! fuild=1 , not fluid = 0
                depth = (self%nz - k + 1) * self%dz
                uster = SQRT(9.8 * depth * self%ib)
                DO kk = k, self%nz
                    zc = (kk-k)*self%dz + 0.5*self%dz
                    nut(j, kk) = self%kappa*uster*depth*zc/depth*(1-zc/depth)
                END DO
                EXIT
            ELSE
                nut(j, k) = 0.0
            END IF
        END DO
    END DO

end function

!> @brief calculate initial condition
!! @detail log-low
!! @param 
!! @return 
function getInitialU(self, cellcnd) RESULT(un)
    class(model) :: self
    INTEGER, INTENT(IN) :: cellcnd(self%ny, self%nz)
    DOUBLE PRECISION :: un(self%ny, self%nz)
    INTEGER :: j, k, kk
    DOUBLE PRECISION :: uster, depth, zc

    DO j = 1, self%ny
        DO k = 1, self%nz
            IF(mod(cellcnd(j, k), 10) == 1)THEN ! fuild=1 , not fluid = 0
                depth = (self%nz - k + 1) * self%dz
                uster = SQRT(9.8 * depth * self%ib)
                DO kk = k, self%nz
                    zc = (kk-k)*self%dz + 0.5*self%dz
                    un(j, k) = uster/self%kappa*LOG(30.0*zc/self%ks)
                END DO
                EXIT
            ELSE
                un(j, k) = -9999.0
            END IF
        END DO
    END DO

end function


!> @brief calculate velosity distribution
!! @detail reverse y and z in python
!! @param 
!! @return 
subroutine calSectU(self, un, cellcnd, nut)
    class(model) :: self
    INTEGER, INTENT(IN) :: cellcnd(self%ny, self%nz)
    DOUBLE PRECISION, INTENT(IN) :: nut(self%ny, self%nz)
    DOUBLE PRECISION, INTENT(INOUT)  :: un(self%ny, self%nz)
    INTEGER :: k, j
    DOUBLE PRECISION :: termYp,termYpc,termYm,termYmc,termZp,termZpc,termZm,termZmc
    DOUBLE PRECISION :: cfbed, cfwall, dy, dz
    DOUBLE PRECISION :: tmp
    LOGICAL :: isLBank, isRBank, isBottom

    dy = self%dy
    dz = self%dz
    cfbed  = (LOG(30.0*0.5*dz/self%ks)/self%kappa)**(-2.0)
    cfwall = (LOG(30.0*0.5*dy/self%ks)/self%kappa)**(-2.0)

    DO k = 1, self%nz
        DO j = 1, self%ny
            IF(mod(cellcnd(j, k), 10) == 0) CYCLE ! fuild=1 , not fluid = 0
            isLBank  =.FALSE. ; IF( mod( cellcnd(j, k)/1000, 10 ) == 1 ) isLBank  = .TRUE.
            isBottom =.FALSE. ; IF( mod( cellcnd(j, k)/100 , 10 ) == 1 ) isBottom = .TRUE.
            isRBank  =.FALSE. ; IF( mod( cellcnd(j, k)/10  , 10 ) == 1 ) isRBank  = .TRUE.

            IF(isLBank)THEN
                termYm  = 0.0
                termYmc = -cfwall*un(j, k)/dy
            ELSE
                tmp = -0.5*(nut(j, k) + nut(j-1, k))/dy**2.0
                termYm  = -tmp*un(j-1, k)
                termYmc =  tmp
            END IF

            IF(isRBank)THEN
                termYp  = 0.0
                termYpc = -cfwall*un(j, k)/dy
            ELSE
                tmp = 0.5*(nut(j, k) + nut(j+1, k))/dy**2.0
                termYp  =  tmp*un(j+1, k)
                termYpc = -tmp
            END IF

            IF(isBottom)THEN
                termZm  = 0.0
                termZmc = -cfbed*un(j, k)/dz
            ELSE
                tmp = -0.5*(nut(j, k) + nut(j, k-1))/dz**2.0
                termZm  = -tmp*un(j, k-1)
                termZmc =  tmp
            END IF

            IF(k == self%nz)THEN !water surface
                termZp  = 0.0
                termZpc = 0.0
            ELSE
                tmp = 0.5*(nut(j, k) + nut(j, k+1))/dz**2.0
                termZp  =  tmp*un(j, k+1)
                termZpc = -tmp
            END IF

            un(j, k) = (9.8*self%ib + termYp + termYm + termZp + termZm) &
                      /(- termYpc - termYmc - termZpc - termZmc)
        END DO
    END DO
end subroutine

end module classModel


!> @brief common module
!! @detail
module commonModel
    use classModel
    CLASS(model), ALLOCATABLE :: m
end module

!> @brief Constructor
!! @detail reverse y and z in python
!! @param 
!! @return 
subroutine initilize(ib, ks, kappa, ny, nz, dy, dz)
    !below 1line only ifort
    !DEC$ ATTRIBUTES DLLEXPORT :: initilize
    use commonModel
    implicit none
    INTEGER :: ny, nz
    DOUBLE PRECISION :: ib, ks, kappa, dy, dz

    ALLOCATE(m)
    CALL m%initilize(ib, ks, kappa, ny, nz, dy, dz)
end subroutine

!> @brief calculation
!! @detail reverse y and z in python
!! @param 
!! @return 
subroutine calSectionProfile(cellcnd, un)
    !below 1line only ifort
    !DEC$ ATTRIBUTES DLLEXPORT :: calSectionProfile
    use commonModel
    implicit none
    INTEGER, INTENT(IN) :: cellcnd(m%ny, m%nz)
    DOUBLE PRECISION, INTENT(OUT) :: un(m%ny, m%nz)

    CALL m%calSectionProfile(cellcnd, un)
end subroutine

!> @brief Destructor
!! @detail
!! @param 
!! @return 
subroutine finalize()
    !below 1line only ifort
    !DEC$ ATTRIBUTES DLLEXPORT :: finalize
    use commonModel
    implicit none
    DEALLOCATE(m)
end subroutine
