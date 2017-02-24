module ellipticalintegrals
! Code written by VdVen. From numerical recipes (?).
use numeric_kinds
implicit none
private

  !private :: rf, rd ! Carslon circular and elliptic functions
  public :: ellf, elle ! Legendre ellptic integrals

contains



  SUBROUTINE ellf(pp,kk,res)
    REAL (kind=dp), intent(in) :: pp,kk
    REAL (kind=dp), intent(out) :: res
    stop 'please insert the numerical recipes ellf, rf & rd source here.'

  END SUBROUTINE ellf

  !!! elliptic integral of the second kind !!!
  !!! Def.: E(phi,k) = int(sqrt(1-k^2*sin^2theta),0..phi) !!!
  SUBROUTINE elle(pp,kk,res)
    REAL (kind=dp), intent(in) :: kk,pp
    REAL (kind=dp), intent(out) :: res
    stop 'please insert the numerical recipes elle, rf & rd source here.'
   
  END SUBROUTINE elle
 
end module ellipticalintegrals