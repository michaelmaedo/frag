!
! Copyright (C) 2016 - 2017 Michael A. Maedo
!
module parameters_mod

    use precision_mod, only: DOUBLE, QUAD

    save
    real(DOUBLE), parameter :: PI    =  3.141592653589793238462643383279503_QUAD
    real(DOUBLE), parameter :: ZERO  =  0.000000000000000000000000000000000_QUAD
    real(DOUBLE), parameter :: ONE   =  1.000000000000000000000000000000000_QUAD
    real(DOUBLE), parameter :: TWO   =  2.000000000000000000000000000000000_QUAD
    real(DOUBLE), parameter :: THREE =  3.000000000000000000000000000000000_QUAD
    real(DOUBLE), parameter :: FOUR  =  4.000000000000000000000000000000000_QUAD
    real(DOUBLE), parameter :: FIVE  =  5.000000000000000000000000000000000_QUAD
    real(DOUBLE), parameter :: SIX   =  6.000000000000000000000000000000000_QUAD
    real(DOUBLE), parameter :: SEVEN =  7.000000000000000000000000000000000_QUAD
    real(DOUBLE), parameter :: EIGHT =  8.000000000000000000000000000000000_QUAD
    real(DOUBLE), parameter :: NINE  =  9.000000000000000000000000000000000_QUAD
    real(DOUBLE), parameter :: TEN   = 10.000000000000000000000000000000000_QUAD

    real(DOUBLE), parameter :: ONE_HALF    = ONE/TWO
    real(DOUBLE), parameter :: ONE_THIRD   = ONE/THREE
    real(DOUBLE), parameter :: ONE_QUARTER = ONE/FOUR
    real(DOUBLE), parameter :: ONE_FIFTH   = ONE/FIVE
    real(DOUBLE), parameter :: ONE_SIXTH   = ONE/SIX
    real(DOUBLE), parameter :: ONE_TENTH   = ONE/TEN

end module parameters_mod
