! This file contains a set of subroutines used by the main code rancat_jswml()
! that the user will need to customize
!
!
! The routines contained have the following purpose:
! 1) kcorr(z) is the k-correction (assumed to the same for all galaxies, ie an average)
!    in the band in which the catalogue is selected.
! 2) ecorr(z,u) is the corresponding e-correction and u is the luminosity evolution
!    parameter that the code and fit in the maximum likelihood procedure
! 3) decorru_du(z,u) the differential of the above function wrt to u.
! 4) The number density evolution function P(z,a). This should be a continuous
!    and differentiable function.
!
!    Notes:
!       If you don't want to fit for the evolution parameters u and a then you
!    can do this by giving zero width to the Gaussian priors. Ie set
!    prior%u=0.0 and prior%a=0 in the main calling program. In this case
!    u=0 a=0 are kept fixed. As you are free to define the evolution functions
!    ecorr(z,u) and P(z,a) you don't have to make u=a=0 correspond to
!    ecorr(z,0)=0 and P(z,0)=1. Thus if you want to have fixed and none zero
!    amount of evolution you should define the functions ecorr() and P() to
!    return this amount of evolution when u=a=0. Similarly if you want to fit
!    for the evolution but with Gaussian priors centred on a default amount
!    of evolution rather than zero evolution you should set up ecorr() and P()
!    in the same way.
!
! 5) The module galaxy_datatype is used to define the fortran90 datatype
!    that carries the information about each property carried by each 
!    galaxy in the genuine and random catalogue. An advantage of this
!    method of generating a random catalogue is that all the original
!    galaxy properties can be carried by the cloned galaxies in the random
!    catalogue. However to do this you must specify what the properties
!    are and read their values from the input catalogue.
!
! 6) The next complication is that some properties such as apparent magnitudes
!    are distance/redshift dependent. Thus when the cloned galaxy is placed
!    at a different redshift to that of the original these properties have to 
!    be adjusted for the change in redshift. Some properties such as an absolute 
!    magnitude (assuming you don't want to change it to be conistent with the fitted
!    e-correction) or a rest frame velocity dispersion or equivalent width
!    need no correcting and will automatically be propagated to the random
!    catalogue. In other cases you need to insert code to make this modification
!    in the routine user_transformations below.
!
!

!-----------------------------------------------------------------------------
! 1)
!
! The relation between apparent, mag, and absolute magnitude, absmag is
!
! mag = absmag + 5log10(dlum/10pc) + kcorr(z)  + ecorr(z)
!
! where kpluse(z) = -2.5 log (1+z) - 2.5 log (L_nu_0(1+z)^z=0/L_nu_0)
!                                  - 2.5 log (L_nu_0(1+z)^z=z/L_nu_0(1+z)^z=0)
!      Hogg et al (2002) (on astro-ph only)                    
!
!real function kcorr(z)
!  implicit none
!  real, intent(in) :: z
!  kcorr= 1.39*z**2+0.87*z  !this is a fit for SDSS r-band selected galaxies
!  return
!end function kcorr

!!! Gama version
real function kcorr(z)
  implicit none
  real, intent(in) :: z
  real :: k
  
  ! Hardcoded kcorrvals
  real, dimension(5) :: kcorrvals = [0.20848, 1.0226, 0.52366, 3.5902, 2.3843]
  integer :: i

  k = 0.0
  do i = 0, 4
     k = k + kcorrvals(i+1) * (z - 0.2)**i
  end do
  kcorr = k
  return
end function kcorr


!
! 2)
real function ecorr(z,u)
  implicit none
  real, intent(in) :: z,u
  ecorr=(u-1.75)*z !set so that the default, u=0, evolution is ecorr=-1.75*z for gama
end function ecorr
!
!
! 3)
! Differential of the above function with respect to the
! the evolution parameter u
real function decorr_du(z,u)
  implicit none
  real, intent(in) :: z,u
  decorr_du=z +0.0*u  !written in this bizzare way as
                      !some compilers complain if u is not referenced.
  return
end function decorr_du
!
! 4)
! -----------------------------------------------
!   Phi^* is assumed to evolve as Phi^*(z)=P(z) Phi^*(0)
!   with the parameterized functional form of P(z) given below
real function P(z,a)
  implicit none
  real :: z,a
  P=exp((a+0.1)*z) !set so that the default, a=0, amount of 
!                   !density evolution is P=exp(0.18*z)
  return
end function P
!
! 5)
module galaxy_datatype
  implicit none
  ! Derived data type for storing all the properties for one galaxy in the
  ! genuine or random catalogue
  type galaxy
     real :: z,mag !redshift and apparent magnitude in the selection band
                   !these must be present and their values set before rancat_jswml
                   !is called.
     real :: absmag,dm,zmax !corresponding absolute magnitude, distance modulus and
             !zmax. These should not be set by the user as rancat_jswml will set them
             !to be consistent with the specified cosmology, survey limits and
             !(fitted) e-correction. The user might want to store the output values.
     real :: dabsmag_du,pvmax,pvmaxeff,weight,pv,v !Also required by rancat_jswml
!            but probably not of interest to the user.
     integer :: id  !as an example of an extra property we have galaxy id
!    real :: Xmag ! Insert here other galaxy properties you want to propagate to the
!            random catalogue. Note that if they are redshift/distance dependent 
!            properties you will also need to add code to user_transformations (below)
!             to correct them for the change in redshift/distance.    
  end type galaxy
end module galaxy_datatype
!
! 6)
!
!  This code is required to modify any additional user extra supplied galaxy properties
!  that are redshift/distance dependent.
!
subroutine user_transformations(clone,original,u)
  use galaxy_datatype
  implicit none
  real,intent(in) :: u  ! luminosity evolution parameter
  type(galaxy),intent(in) :: original ! all properties of the original galaxy
  type(galaxy),intent(inout) :: clone ! properties of the clone with all but
!                                      ! the extra user defined properties already
!                                      ! updated to the new assigned redshift
!
! No code is required if the galaxy properties carried by the original galaxy
! are only the selection band apparent magnitude, redshift and additional
! redshift/distance independent properties.
!
! An example of the sort of transformation you might want to supply 
! if what you need to do is modify an apparent magnitude in some extra "X"band.
!
!  clone%Xmag=original%Mag-original%dm-kcorrX(original%z)-ecorr(original%z) &
!        & +clone%dm-kcorrX(clone%z)-ecorr(clone%z)
!
!  I.E. You adjust the apparent magnitude for the change in distance modulus (dm)
!  and change in the k- and e-corrections. The distance modulus at each redshift
!  is already stored for you, but you need to supply the k-correction kcorrX(z)
!  for this band. Also in the above we assume the e-correction is the same in
!  this band as in the selection band.
!
  clone%z=clone%z +u*0.0*original%z !dummy null code just to prevent some compiler complaints
end subroutine user_transformations
