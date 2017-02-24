module dmpotent
  use numeric_kinds
  use initial_parameters
  implicit none
  private


  ! Generic module for black hole mass and dark halo potential
  ! calls triaxpotent

  ! RvdB 17 Mar 2010
  
  ! Bugs

  ! calculate potential phi at (x,y,z) 
  public:: dm_potent
  
  ! calculate accel ax,ay,ax at (x,y,z)
  public:: dm_accel
  
  ! setup the constants for the potential
  public:: dm_setup

  ! setup the constants for the potential
  public:: dm_stop

  real (kind=dp), private :: rhoc, rc


contains

subroutine dm_setup()
   use triaxpotent, only : tp_setup
   real(kind=dp) :: darkmass
   call tp_setup()

   select case (dm_profile_type)
	case (0)
		print*, 'No additional DM halo' 
	case (1)
        !	     dmparam(1) = concentration
        !        dmparam(2) = dm_fraction   (fraction of DM mass within R200 radius)
	    if (n_dmparam .ne. 2)  stop 'wrong number of NFW halo parameters'
	    print*,'Parameters of NFW concentration and fraction',dmparam(1),dmparam(2)
 

        ! Parameters for NFW profile
	    rhoc = (200.0_dp/3.0_dp)*rho_crit*dmparam(1)**3/  & 
	        ( log(1.0_dp+dmparam(1)) - dmparam(1)/(1.0_dp+dmparam(1)))


	    rc = (3.0_dp/(800.0_dp*pi_d*rho_crit*dmparam(1)**3) &
	           * dmparam(2) * totalmass)**(1.0_dp/3.0_dp)
	                                                                                  
	    ! 12 Oct 2011: LW found unit conversion bug in print statment
		print*, "Parameters of NFW potential (rho_c in solarmass/km^3 and r_c in km): ", rhoc, rc
		! Calculate M200, in Msun
		darkmass = 800_dp*pi_d/3_dp*rho_crit*rc**3*dmparam(1)**3

		print*, "Total stellar mass is (Msun): ", totalmass
		print*,"Total dark halo mass (M200 in Msun): ", darkmass
		
	case (2)
	    !	     dmparam(1) = rhoc
        !        dmparam(2) = rc      
 	    if (n_dmparam .ne. 2)  stop 'wrong number of Hernquist halo parameters'
        rhoc=dmparam(1)
		rc  =dmparam(2)
		print*,'Parameters of Hernquist profile', rhoc, rc 
	case (3) 
	   print*, "  * triaxial cored logarithmic potential. " 
	   ! from Thomas et al. 2005  & B&T 1987 (p. 46)  
	    if (n_dmparam .ne. 4)  stop 'wrong number of halo parameters'
	   
       print*, "  Vc (km/s), rho (kpc,km):",dmparam(1), dmparam(2),dmparam(2)*parsec_km*1d3   
       if (dmparam(1) .le. 0.0_dp) stop 'VC < 0'
       if (dmparam(2) .le. 0.0_dp) stop 'VC < 0'

       print*, "  flattening p & q", dmparam(3),dmparam(4)
	   if (dmparam(3) .gt. 1.0_dp .or. dmparam(4) .le. 0.0_dp &
	          .or. dmparam(3) .lt. dmparam(4)) stop ' Flattening is not 0<q<=p<=1'
  
       ! turning p and q into p^2 and q^2
       dmparam(3)=dmparam(3)**2
       dmparam(4)=dmparam(4)**2

       ! Turning Core radius from kpc to km^2  
       dmparam(2)=(dmparam(2)*parsec_km*1d3)**2.0_dp    

       ! Turning VC into km^2   
       dmparam(1)=dmparam(1)**2.0_dp    
 
	!case (4)
		!read (unit=13, fmt=*) dm_profile_rhoc, dm_profile_parameter, r200, dm_rho_crit, dm_profile_a, dm_profile_b, dm_profile_c
		!print*,'Parameters of NFW(triaxial approximation)', dm_profile_rhoc, dm_profile_parameter, r200, dm_rho_crit
		!print*,'rho crit:', rho_crit
		!print*, dm_profile_a, dm_profile_b, dm_profile_c
		!dm_profile_axes_length = dm_profile_a**2 + dm_profile_b**2 + dm_profile_c**2
		!print*,'unnormalized axes length:', dm_profile_axes_length
		!dm_profile_a = dm_profile_a *sqrt(3/dm_profile_axes_length) 
		!dm_profile_b = dm_profile_b *sqrt(3/dm_profile_axes_length) 
		!dm_profile_c = dm_profile_c *sqrt(3/dm_profile_axes_length) 
		!dm_profile_axes_length = dm_profile_a**2 + dm_profile_b**2 + dm_profile_c**2
		!print*,' normalized axes length:', dm_profile_axes_length
    end select

end subroutine dm_setup

subroutine dm_stop()
   ! use triaxpotent, only : tp_stop
   ! call tp_stop()  ! function does not exist, but should.
end subroutine dm_stop

subroutine dm_potent(x,y,z,pot)
  use triaxpotent, only : tp_potent
  use initial_parameters
  real(kind=dp), intent(in) ::  x, y, z
  real(kind=dp), intent(out):: pot
  real(kind=dp) :: d,d2

  d2 = x*x + y*y + z*z

  call tp_potent(x,y,z,pot)

  ! add Plummer style black hole 
  pot = pot + grav_const_km *  xmbh / sqrt( d2 + softl_km*softl_km ) 

  select case (dm_profile_type)
	case (0)
            !blank
	case (1)
     ! add NFW dark halo
	 !d =sqrt(d2) 
	 if (sqrt(d2)/rc .ge. 1.0) then 
     pot = pot + 4.0_dp*pi_d*grav_const_km*rhoc*rc**3/sqrt(d2) * log(1.0_dp+sqrt(d2)/rc)
     else 
     ! indentity log (1+x) = 2* atanh(x/(2+x)) , required when 0<x<<1 
	  pot = pot + 4.0_dp*pi_d*grav_const_km*rhoc*rc**3/sqrt(d2) * 2*atanh( (sqrt(d2)/rc)/(2+sqrt(d2)/rc))  
	 endif
    case (2)
     ! Hernquist  
     d =sqrt(d2)
	 pot = pot + 4.0_dp*pi_d*grav_const_km*rhoc*rc**2 / (2*(1+d/rc))
    case(3)
     ! cored logarihtmic
     ! phi=1/2*vc^2*log(rc^2+x^2+y^2/p^2+z^2/q^2)
     pot = pot - 0.5_dp *dmparam(1) *  log(dmparam(2) + (x**2.0_dp + y**2.0_dp /dmparam(3) + z**2.0_dp /dmparam(4)))
     !if ((x**2.0_dp + y**2.0_dp /dmparam(3) + z**2.0_dp /dmparam(4))/dmparam(2) .le. 1.0d-14)  stop ' potential fails log(x+y) test' 
     ! density rho is the laplacian of the potential phi:
	 !rho = -vc**2*(-p**4*q**4*rc**2+x**2*p**4*q**4-p**2*q**4*y**2-p*        $
	 !    *4*q**2*z**2-q**4*rc**2*p**2-q**4*x**2*p**2+y**2*q**4-q**2*z**2*   $
	 !    p**2-p**4*rc**2*q**2-p**4*x**2*q**2-p**2*y**2*q**2+z**2*p**4)      $
	 !    /(rc**2*p**2*q**2+x**2*p**2*q**2+y**2*q**2+z**2*p**2)**2
	
	
	!case (4)
   	! r = d
	! ra = dm_profile_parameter
	! re = sqrt(x**2/dm_profile_a**2 + y**2/dm_profile_b**2 + z**2/dm_profile_c**2)
	! rtilde = (ra+r)*re / (ra+re)
    !
	! pot = pot + 4.0_dp * pi_d * grav_const_km*dm_profile_rhoc*dm_profile_parameter**3/rtilde * log(1.0_dp + rtilde/dm_profile_parameter)
  end select

end subroutine dm_potent

!+++++++++++++++++++++++++++++++++++
subroutine dm_accel(x,y,z,vx,vy,vz)
  use triaxpotent, only : tp_accel
  use initial_parameters
  real(kind=dp), intent(in) ::  x, y, z
  real(kind=dp), intent(out):: vx,vy,vz
  !------------------------------------
   real(kind=dp) :: t,t1,t2,t3,t4,d2,acceleration_r
!  integer(kind=i4b) :: p,r

   call tp_accel(x,y,z,vx,vy,vz)

   d2 = x*x + y*y + z*z   
   ! Add Plummer style blackhole.
   t  =  - grav_const_km *  xmbh * ( d2 + softl_km*softl_km ) ** (-3.0_dp/2.0) 
   vx = vx + x * t 
   vy = vy + y * t
   vz = vz + z * t
   
  select case (dm_profile_type)
	case (0)
            !blank
	case (1)   ! Add NFW dark matter halo
     t1 = -4.0_dp*pi_d*grav_const_km*rhoc*rc**3/d2
     ! indentity log (1+x) = 2* atanh(x/(2+x)) , required when 0<x<<1
     if (sqrt(d2)/rc .ge. 1.0) then 
	   t2 = log(1.0_dp + sqrt(d2)/rc) 
     else 
	   t2= 2*atanh( (sqrt(d2)/rc)/(2+sqrt(d2)/rc))
     endif
     t3 = (sqrt(d2)/rc)/(1.0_dp+sqrt(d2)/rc)
     vx = vx + x/sqrt(d2) * t1*(t2 - t3)
     vy = vy + y/sqrt(d2) * t1*(t2 - t3)
     vz = vz + z/sqrt(d2) * t1*(t2 - t3)   
      if (x/sqrt(d2) * t1*(t2 - t3) .gt. 0 .and. x .gt. 0)  then
	   print*,vx,x/sqrt(d2) * t1*(t2 - t3) ,d2
	   print*,t1,t2,t3 
	   print*,sqrt(d2)/rc,d2,rc
	   stop 'NFW accelations flipped sign'
      endif
    case (2)     !Hernquist
  	 acceleration_r = -2.0_dp * pi_d * grav_const_km * rhoc*rc/(1+sqrt(d2)/rc)**2
	 vx = vx + x/sqrt(d2) * acceleration_r
	 vy = vy + y/sqrt(d2) * acceleration_r
	 vz = vz + z/sqrt(d2) * acceleration_r
    case (3) ! cored logarihtmic
     ! phi=1/2*vc^2*log(rc^2+x^2+y^2/p^2+z^2/q^2)
     ! vx= diff(phi,x)
     vx= vx - dmparam(1)* x             / (dmparam(2)+x*x+ y*y/dmparam(3) + z*z/dmparam(4))
     vy= vy - dmparam(1)*(y/dmparam(3)) / (dmparam(2)+x*x+ y*y/dmparam(3) + z*z/dmparam(4))
     vz= vz - dmparam(1)*(z/dmparam(4)) / (dmparam(2)+x*x+ y*y/dmparam(3) + z*z/dmparam(4))

    !case (4)     !Vogelsberger
	!	    d = sqrt(d2)
	!		! see Vogelsberger+ 2008
	!		r = d
	!		ra = dm_profile_parameter
	!		re = sqrt(x**2/dm_profile_a**2 + y**2/dm_profile_b**2 + z**2/dm_profile_c**2)
	!		rtilde = (ra+r)*re / (ra+re)
    !
	!		t1 = -4.0_dp*pi_d*grav_const_km*dm_profile_rhoc*dm_profile_parameter**3/rtilde**2
	!		t2 = log(1.0_dp + rtilde/dm_profile_parameter)
	!		t3 = (rtilde/dm_profile_parameter)/(1.0_dp+rtilde/dm_profile_parameter)
    !
	!		acceleration_rtilde = t1*(t2 - t3)
    !
	!		drtildedx = 1 / (ra+re)**2 * (1/re * x / dm_profile_a**2 * (ra+r) * ra + 1/r * x * (ra+re)*re)
	!		drtildedy = 1 / (ra+re)**2 * (1/re * y / dm_profile_b**2 * (ra+r) * ra + 1/r * y * (ra+re)*re)
	!		drtildedz = 1 / (ra+re)**2 * (1/re * z / dm_profile_c**2 * (ra+r) * ra + 1/r * z * (ra+re)*re)
	!		vx = vx + acceleration_rtilde * drtildedx
	!		vy = vy + acceleration_rtilde * drtildedy
	!		vz = vz + acceleration_rtilde * drtildedz

	end select


end subroutine dm_accel

end module dmpotent
