module datatypes
  implicit none
  ! Derived type for carrying the model parameters and the associated posterior
  ! probability.
  ! The components are the density evolution parameter, a, the luminosity
  ! evolution parameter, u, and the Lagrangian multiplier parameter, mu
  ! and finally their associated posterior probability.
  type model_parameters
     real :: a,u,mu,p_post
  end type model_parameters
  ! Derived type carrying the parameters that specify the survey selection
  type survey_spec
     real :: zmin,zmax,magfaint,solid_angle
  end type survey_spec
  ! Derived type carrying the parameters that specify the priors
  ! include a logical variable to switch between Gaussain priors
  ! an an alternative B-spline
  type prior_parameters
     real :: u,a,fourPiJ3
     logical :: spline
  end type prior_parameters
  ! Derived type carrying redshift binned quantities that we need to pass internally
  ! between routines
  type binned_z
     real :: z,n,n_ran,pvol,vol,delta
  end type binned_z
  ! Derived type carrying magnitude binned quantities that we need to pass internally
  ! between routines
  type binned_mag
     real :: n,dn_du
  end type binned_mag
 end module datatypes


module likelihood_and_priors
contains
!Likelihood function of the model given the data
!
! ln(Likelihood) = Sum_gals ln V_pdelta_p + ln lf_i
!               where p and i are the redshift and luminosity bins in
!               which each galaxy falls.
! Equivalent to
! ln(Likelihood) = Sum_p nofz_p ln V_pdelta_p + Sum_i nmag_i ln lf_i
!               where p and i are the redshift and luminosity bins in
!               which each galaxy falls.
!
!NB Actually log(Likelihood) and neglecting overall constants
! including the final term in the likelihood expression 
real function log_likelihood(zbin,delta,nmag,lf,nzbin,nlf)
  use datatypes
  implicit none
  type(binned_z),intent(in) :: zbin(:) 
  integer,intent(in) :: nlf,nzbin
  real,intent(in) :: delta(nzbin),nmag(nlf),lf(nlf)
  !  The masked sum avoids taking log(0) for terms whose multiplier, nofz() or nmag(), 
  !  are zero anyway
  log_likelihood=sum(zbin(:)%n*log(zbin(:)%pvol*delta(:)),mask=zbin(:)%n.gt.0) + &
       &         sum(nmag(:)*log(lf(:)),mask=nmag(:).gt.0) 
end function log_likelihood


! Prior probabilities on the evolution parameters par%a and par%u
! and on the binned overdensity delta(:)
! 
real function logp_prior(par,prior,zbin,delta,nzbin)
  use datatypes
  implicit none
  type(binned_z),intent(in) :: zbin(:) 
  type(model_parameters),intent(in) :: par
  type(prior_parameters),intent(in) :: prior
  integer,intent(in) :: nzbin
  real,intent(in) :: delta(nzbin)
  real :: sigmasq(nzbin)
  logp_prior= 0.0  !initialize before incrementing by the various 
!                  !contributions to the logp_prior
  if (prior%spline) then
     ! B-spline constraints on both u and a if used.
     if (prior%a.gt.0.0) logp_prior= logp_prior+lnBspline(par%a/prior%a)
     if (prior%u.gt.0.0) logp_prior= logp_prior+lnBspline(par%u/prior%u)
  else   
     ! Gaussian constraints on both u and a if used.
     if (prior%a.gt.0.0) logp_prior= logp_prior-0.5*(par%a/prior%a)**2
     if (prior%u.gt.0.0) logp_prior= logp_prior-0.5*(par%u/prior%u)**2 
 end if    
! Gaussian constraint on the overdensity in each redshift bin based
! on Peebles cluster model
  if (prior%fourPiJ3.ge.0.0) then
     where (zbin(:)%n_ran.gt.0.0)
        sigmasq(:)=(1.0+zbin(:)%n_ran*prior%fourPiJ3/zbin(:)%vol)/zbin(:)%n_ran  
     elsewhere   
        sigmasq(:)=1.0e+30
     end where
     logp_prior=logp_prior-sum(0.5*(delta(:)-1.0)**2/sigmasq(:))
  end if   
end function logp_prior

! ln(Bspline) where Bspline is spline approximation to a Gaussian
!
real function lnBspline(x)
  real,intent(in)  :: x
  real,parameter :: inf=1.0e+30
  real :: u
  u=abs(x/sqrt(3.0))
  if (u.le.1.0) then
     lnBspline=log(1.0-1.5*u**2+0.75*u**3)
  else if (u.lt.2.0) then
     lnBspline=log(0.25*(2.0-u)**3)
  else  !Bspline would be zero for u>2 but we want to avoid infinities
     lnBspline=-inf
  end if   
end function lnBspline
!
! The derivative of lnBspline wrt x
real function dlnBspline_dx(x)
  real,intent(in)  :: x
  real,parameter :: inf=1.0e+30
  real :: u,dlnBspline_du
  u=abs(x/sqrt(3.0))
  !first calculate the derivative wrt u
  if (u.le.1.0) then
     dlnBspline_du= (-3.0*u+2.25*u**2)/(1.0-1.5*u**2+0.75*u**3)
  else if (u.lt.2.0) then
     dlnBspline_du=-3.0/(2.0-u)
  else  !Bspline would be zero for u>2 but we want to avoid infinities
     dlnBspline_du=-inf
  end if  
  ! now multiply by du/dx whose sign depends on the sign of x
  if (x.gt.0) then
     dlnBspline_dx= dlnBspline_du/sqrt(3.0)
  else   
     dlnBspline_dx=-dlnBspline_du/sqrt(3.0)
  end if   
end function dlnBspline_dx

end module likelihood_and_priors

module rancat_subs
contains

! The constraint equation whose root must be found to
! determine the density evolution parameter a
!
real function acon(a,a_ref,sigmaA,mu,nzbin,zbin,spline)
  use datatypes
  use likelihood_and_priors
  implicit none
  type(binned_z) :: zbin(:)
  real,intent(in) :: a_ref,a,sigmaA,mu
  integer,intent(in) :: nzbin
  logical,intent(in) ::spline
  integer :: iz
  real, parameter :: epsilon=1.0e-03
  real:: dlnPda,P,dP,da,nofz_scl,nbartothat,nbartothat_ref
!
! total galaxy number in the random catalogue with the reference value of a
! and the value one would get if all redshifts bins increased in proportion
! to the fractional change in P(zbin).
! As the actual number is kept fixed we can use this as a normalization constraint
! when estimating the change in any given bin.
  nbartothat_ref=0.0
  nbartothat=0.0
  do iz=1,nzbin
       nbartothat_ref=zbin(iz)%n_ran+nbartothat_ref
       nbartothat=zbin(iz)%n_ran*P(zbin(iz)%z,a)/P(zbin(iz)%z,a_ref)+nbartothat
  end do   
!
  if (sigmaA.gt.0.0) then
     if (spline) then
        acon=dlnBspline_dx(a/sigmaA)/sigmaA !contribution from spline prior
     else   
        acon=-a/sigmaA**2   !contribution from Gaussian prior on a
     end if
  end if

  do iz=1,nzbin
     ! Simple numerical differentiation of P(z) function wrt a
     if (abs(a).gt.epsilon) then
        da=epsilon*a
     else
        da=epsilon
     end if
     dP=P(zbin(iz)%z,a+0.5*da)-P(zbin(iz)%z,a-0.5*da)
     if (dP.eq.0.0) then
        dlnPda= 0.0
     else   
        dlnPda= dP/(P(zbin(iz)%z,a)*da)
     end if
!    Assumes nofz_ran(iz) is proportional to P(z,a), but subject to a normalization
!    constraint
     nofz_scl=(P(zbin(iz)%z,a)/P(zbin(iz)%z,a_ref))*(nbartothat_ref/nbartothat)
     acon = acon+ dlnPda * (zbin(iz)%n-zbin(iz)%n_ran*nofz_scl*(zbin(iz)%delta+mu))
  end do
end function acon

! The constraint equation whose root must be found to
! determine the luminosity evolution parameter u
!
real function ucon(u,sigmaU,mu,dmag,faint_maglimit,nzbin,nlf,lf,magbin,mag,zbin,zref,spline)
  use datatypes
  use cosmology
  use likelihood_and_priors
  implicit none
  integer,intent(in):: nzbin,nlf
  real,intent(in):: u,sigmaU,mu,dmag,faint_maglimit,zref
  real,intent(in) :: lf(nlf),magbin(nlf)
  logical,intent(in) :: spline
  type(binned_z),intent(in) :: zbin(:)
  type(binned_mag),intent(in) :: mag(:)
  integer :: iz,imag
  real :: absmag,magoffset,zplus1,rc,dlum,pvol,h,lf_int
  real*8 :: ucon1,ucon2
  real,external :: kpluse,dkpluse_du
!
  if (sigmaU.gt.0) then
     if (spline) then
        ucon=dlnBspline_dx(u/sigmaU)/sigmaU !contribution from spline prior
     else
        ucon=-u/sigmaU**2  !term from the Gaussian prior
     end if
  end if
  !
  !Add the term from the rate of change of the denominator in the likelihood 
  !expression combined with the Lagrange constraint
  ucon1=0.0
  do iz=1,nzbin
     ! Find the absolute magnitude corresponding to the faint apparent magnitude
     ! limit at this redshift
     zplus1=1.0+zbin(iz)%z
     absmag=0.0 !set so that that the value returned for magoffset=mag-absmag
     call search_distancemag(rc,zplus1,dlum,pvol,magoffset,2,absmag,kpluse,u,zref)  
     absmag=  faint_maglimit-magoffset !absolute corresponding to app. mag limit
     !interpolate between the nearest bins of the LF to estimate its value absnag
     imag= int(1.0+(absmag-magbin(1))/dmag) !find corresponding magnitude bin
     if (imag.ge.1 .and. imag.le.nlf) then
        h=(absmag-magbin(imag))/dmag
        if (imag.lt.nlf) then
           lf_int=lf(imag)*(1.0-h)+lf(imag+1)*h  !interpolate
        else
           lf_int=lf(imag)
        end if
        ucon1=ucon1+zbin(iz)%pvol*(zbin(iz)%delta+mu)*lf_int*dkpluse_du(zbin(iz)%z,u)/dmag
     end if
  end do
  !
  !Now add the term from the rate of change of the LF term
  ucon2=0.0
  do imag=1,nlf
     if (mag(imag)%n.gt.0.0) then 
        ucon2=ucon2+mag(imag)%dn_du*log(lf(imag)/dmag)
     end if
  end do
  ucon=ucon1+ucon2+ucon  !combine terms
end function ucon

!Posterior probabilty on all the model parameters
! evolution parameters: u and a
! binned overdensity  : delta(:)
! binned lum. func.   : lf(:)
!
!NB Actually log(Probability) and neglecting overall constants 
real function logp_post(par,prior,zbin,delta,nzbin,nmag,lf,nlf)
  use likelihood_and_priors
  use datatypes
  implicit none
  integer,intent(in) :: nzbin,nlf
  type(model_parameters),intent(in) :: par  
  type(prior_parameters),intent(in) :: prior
  type(binned_z),intent(in) :: zbin(:) 
  real,intent(in) :: delta(nzbin)
  real,intent(in) :: lf(nlf),nmag(nlf)
  logp_post=log_likelihood(zbin,delta,nmag,lf,nzbin,nlf)+ &
       & logp_prior(par,prior,zbin,delta,nzbin)
end function logp_post


! Tabulate the P-evolution weighted volume of each redshift bin
!
subroutine tabulate_binpvolumes(survey,nzbin,zbin,pvolbin,volbin,zref)  
  use datatypes
  use cosmology
  implicit none
  integer,intent(in) :: nzbin
  real,intent(out) :: volbin(nzbin),pvolbin(nzbin)
  type(survey_spec),intent(in) :: survey
  real,intent(in) :: zbin(nzbin),zref
  real :: zlo,zhi,rc,zplus1,dlum,pvol,mag,u_dummy=0.0
  integer :: i,isel
  real,external :: kpluse
  do i=1,nzbin
     if (i.gt.1) then  !lower redshift of bin edge
        zlo=0.5*(zbin(i-1)+zbin(i))
     else
        zlo=survey%zmin
     end if   
     if (i.lt.nzbin) then !upper redshift of bin edge
        zhi=0.5*(zbin(i+1)+zbin(i))
     else
        zhi=survey%zmax
     end if   
     isel=2   !search on redshift
     zplus1=1.0+zhi
     call search_distancemag(rc,zplus1,dlum,pvol,mag,isel,0.0,kpluse,u_dummy,zref)  
     pvolbin(i)=survey%solid_angle*pvol !cumulative volume to upper z (P-weighted)
     volbin(i)=(survey%solid_angle/3.0)*rc**3 !cumulative volume to upper z
     zplus1=1.0+zlo
     call search_distancemag(rc,zplus1,dlum,pvol,mag,isel,0.0,kpluse,u_dummy,zref)   
     pvolbin(i)=pvolbin(i)-survey%solid_angle*pvol !subtract off volume inside lower z
     volbin(i)=volbin(i)-(survey%solid_angle/3.0)*rc**3 !subtract off volume inside lower z
  end do
end subroutine tabulate_binpvolumes

!
! For each of the ncat galaxies in the genuine catalogue compute their absolute, %absmag,
! magnitude (consistent with the kpluse correction specified by the paramerter u),
! distance modulus %dm and their %zmax 
! (maximum redshift at which they would still satisfy the selection bounds)
!
subroutine  derived_gal_props(u,zref,survey,ncat,cat)
  use datatypes
  use galaxy_datatype
  use cosmology
  implicit none
  type(galaxy) :: cat(ncat)
  type(survey_spec),intent(in) :: survey
  integer,intent(in) :: ncat
  real,intent(in) :: u,zref
  integer :: i
  real :: zplus1,rc,dlum,pvol,mag,pvolmin,absmag_dummy
  real,parameter :: epsilon=1.0e-05
  real, external :: kpluse
! Compute pvolmin corresponding to the catalogue zsurvey_min as this has to
! be subtracted off the P-weighted volumes computed from z=0
  zplus1=1.0+survey%zmin !2 indicates search on zplus1
  call search_distancemag(rc,zplus1,dlum,pvolmin,mag,2,absmag_dummy,kpluse,u,zref)  
!
  do i=1,ncat
     !tabulate the distance modulus of each galaxy
     zplus1=1.0+cat(i)%z !2 indicates search on zplus1
     call search_distancemag(rc,zplus1,dlum,pvol,mag,2,cat(i)%absmag,kpluse,u,zref) 
     cat(i)%pv=survey%solid_angle*(pvol-pvolmin)
     cat(i)%dm=25.0+5.0*log10(max(dlum,1.0e-20))
     !compute absmag consistent with kpluse assumption and current value of u
     !this is done by using that (cat(i)%absmag-mag) equals the difference between
     !the actual apparent mag, cat(i)%mag, and the required absolute magnitude.
     cat(i)%absmag=cat(i)%absmag-mag+cat(i)%mag
     !Find the redshift at which this galaxy would fall out of the survey 
     call search_distancemag(rc,zplus1,dlum,pvol,survey%magfaint,5,cat(i)%absmag,kpluse,u,zref)  !5 indicates search on mag 
     if  (zplus1.ge.1.0+survey%zmax) then
        zplus1=1.0+survey%zmax-epsilon
        call search_distancemag(rc,zplus1,dlum,pvol,mag,2,cat(i)%absmag,kpluse,u,zref)  !2 indicates search on zplus1
     end if   
     cat(i)%zmax=zplus1-1.0
     cat(i)%pvmax=survey%solid_angle*(pvol-pvolmin)
  end do
end subroutine derived_gal_props

! Find the value of the Lagrange multiplier parameter, mu, such that weights
! pvmax/(pvmaxeff+mu.pvmax) have a mean of unity.
!
subroutine normalize_weights(cat,ncat,mu)
  use datatypes
  use galaxy_datatype
  implicit none
  integer,intent(in) :: ncat
  type(galaxy),intent(in) :: cat(ncat)
  real,intent(inout) :: mu
  integer,parameter :: niimax=30
  integer :: i,ii
  real :: muprime
  real*8 :: r,dr_dmu,ratio,mudble,tol=1.0e-08 !tolerance for convergence
  muprime=0.0d+00 !initial guess
  mudble =1.0d+30 !mudble=dble(mu) temporarily set so that we enter the loop
  ii=0
  do while (abs(mudble-muprime).gt.tol)
     mudble=muprime
     ! Accumulate the sum of the weights and its derivative wrt mu
     r=0.0
     dr_dmu=0.0
     do i=1,ncat
        ratio=cat(i)%pvmaxeff/cat(i)%pvmax
        r=r+1.0d+00/(ratio+mudble)
        dr_dmu=dr_dmu-1.0d+00/(ratio+mu)**2
     end do
     muprime=mudble+(dble(ncat)-r)/dr_dmu  !Newton-Raphson estimate 
     ii=ii+1
     if (ii.eq.niimax) then
        write(*,*) 'Error: mu not converged after ',ii,' iterations. mu=',muprime
        stop
     end if   
  end do
  mu=muprime
end subroutine normalize_weights

! Generate a random catalogue by cloning galaxies from the genuine catalogue
! and distributing them uniformly within the (P-weighted) Vmax and using the
! current values of the evolution parameters a and u
!
subroutine generate_rancat(iseed,nmult,survey,mu,u,zref,ncat,cat,nrancat,rancat,nran)
  use datatypes
  use galaxy_datatype
  use cosmology
  implicit none
  integer,intent(in) :: ncat,nmult,nrancat
  integer,intent(out) :: nran
  type(survey_spec),intent(in) :: survey
  real,intent(in) ::u,mu,zref
  integer,intent(inout) :: iseed
  type(galaxy),intent(in) :: cat(ncat)
  type(galaxy),intent(out) :: rancat(nrancat)
  real :: rep_target
  integer :: i,ir,nrep,irep
  real :: zplus1,dlum,mag,rc,pvol,pvolmax,pvolmin,absmag_dummy,weight
  real,external :: ran3,kpluse
!
! Compute pvolmin corresponding to the catalogue survey%zmin as this has to
! be subtracted off the P-weighted volumes computed from z=0
  zplus1=1.0+survey%zmin !2 indicates search on zplus1
  absmag_dummy=1.0
  call search_distancemag(rc,zplus1,dlum,pvolmin,mag,2,absmag_dummy,kpluse,u,zref)  
!
  ir=0
  do i=1,ncat
     weight=cat(i)%pvmax/(cat(i)%pvmaxeff+mu*cat(i)%pvmax)
     rep_target=real(nmult)*weight
     ! round to one of the nearest integers such that the mean
     ! value would equal the target value
     if (ran3(iseed).lt.rep_target-int(rep_target)) then
        nrep=int(rep_target)+1
     else
        nrep=int(rep_target)
     end if
     ! Generate a redshift for each of these replicas by selecting
     ! randomly within Vmax. 
     zplus1=1.0+cat(i)%zmax  !2 indicates search on zplus1
     call search_distancemag(rc,zplus1,dlum,pvolmax,mag,2,cat(i)%absmag,kpluse,u,zref)  
     do irep=1,nrep
        ! Copy this galaxy to a new slot in the random catalogue then adjust
        ! values to correspond to the assigned redshift 
        ir=ir+1
        if (ir.gt.nrancat) then
           write(0,*) 'ERROR:',i, ' of ', ncat,' . replica ',irep,' of ',nrep,' total:',ir
           stop 'generate_rancat(): rancat(nrancat) array too small'
        end if   
        rancat(ir)=cat(i) !identical clone
        pvol=pvolmin+(pvolmax-pvolmin)*ran3(iseed) !generate uniform in pvol
        if (pvol.eq.0) then
           write(0,*) 'WARNING: pvol=0.0 generated. Consider setting survey%zmax>0'
           pvol=1.0e-10*pvolmax
           write(0,*) 'reset to pvol=',pvol
        end if
        ! search for corresponding redshift etc 
        call search_distancemag(rc,zplus1,dlum,pvol,mag,4,rancat(ir)%absmag,kpluse,u,zref) 
        ! assign new redshift and replace/modify other redshift dependent quantities
        rancat(ir)%z=zplus1-1.0 
        rancat(ir)%pv= survey%solid_angle*(pvol-pvolmin)
        rancat(ir)%mag=mag
        rancat(ir)%dm=25.0+5.0*log10(max(dlum,1.0e-20))
        !modify addtional properties supplied by the user for the change in redshift
        call user_transformations(rancat(ir),cat(i),u)
     end do
  end do
  nran=ir
end subroutine generate_rancat


!Estimate overdensity, delta, in each redshift bin
!
subroutine estimate_delta(nzbin,zbin,fourPiJ3,delta,delta_raw,sigma)
  use datatypes
  implicit none
  integer,intent(in) :: nzbin
  type(binned_z) :: zbin(:)
  real,intent(in) :: fourPiJ3
  real,intent(out) :: delta(nzbin),sigma(nzbin),delta_raw(nzbin)
  real :: f(nzbin)
  integer :: i
  where (zbin(:)%n_ran.gt.0)
     delta_raw(:)=zbin(:)%n/zbin(:)%n_ran !raw estimate before applying LSS constraint
  elsewhere   
     delta_raw(:)=1.0
  endwhere   
  !Now solve the quadratic 
  ! 0 = - delta_raw + delta + f delta(delta-1)
  ! for delta where delta_raw is the raw estimate above and f is  
  ! 1/(sigma^2 nofz_ran)
  where (zbin(:)%n_ran.gt.0.0)
     sigma(:)=sqrt( (1.0+zbin(:)%n_ran*fourPiJ3/zbin(:)%vol)/zbin(:)%n_ran  )
     f(:)=1.0/(sigma(:)**2*zbin(:)%n_ran)
  elsewhere   
     sigma(:)=1.0e+30
     f(:)=1.0
  end where
  do i=1,nzbin
     if (f(i).gt.5.0e-04) then
        delta(i)=(f(i)-1.0+sqrt((1.0-f(i))**2+4.0*f(i)*delta_raw(i)))/(2.0*f(i))
     else !for smaller f use the Taylor expansion approximation
        delta(i)=delta_raw(i)+f(i)*(delta_raw(i)-delta_raw(i)**2)        
     end if
  end do
end subroutine estimate_delta

! Make a direct estimate of the N(z) expected for the random catalogue
! This involves redistributing the pvmax/(pvmaxeff+mu.pvmax) weight of each
! genuine galaxy uniformly over the accessible P-weighted Vmax.
!
subroutine estimate_nofz(ncat,cat,mu,zbin,dzbin)
  use datatypes
  use galaxy_datatype
  implicit none
  integer,intent(in) :: ncat
  real,intent(in) :: mu,dzbin
  type(galaxy),intent(in) :: cat(ncat)
  type(binned_z),intent(inout) :: zbin(:)
  real :: weight,pvsum
  integer :: i,iz
!
  zbin(:)%n_ran=0.0 !initialize redshift histogram
  do i=1,ncat !loop over the galaxies in the genuine catalogue
!    The contribution they make to each redshift bin of the random catalogue is
!    pvolbin/pvmax  * pvmax/(pvmaxeff+mu.pvmax). Cancelling the two qvmax factors
!    we can write this contribution as pvolbin(iz)*weight where weight is given by
     weight=1.0/(cat(i)%pvmaxeff+mu*cat(i)%pvmax)
     iz=1 ! assign contributions from all bins fully below zmax
     pvsum=0.0
     do while (cat(i)%zmax.gt.zbin(iz)%z+0.5*dzbin)
        zbin(iz)%n_ran=zbin(iz)%n_ran+zbin(iz)%pvol*weight
        pvsum=pvsum+zbin(iz)%pvol
        iz=iz+1
     end do
     ! then assign contribution from the final fractional bin
     ! In pvsum we've accumlated the P-weighted volume from all the fully contributing bins
     ! hence (cat(i)%pvmax-pvsum) is the P-weighted volume of the final fractionally occupied 
     ! bin which has the index iz at which we exited the above loop
     zbin(iz)%n_ran=zbin(iz)%n_ran+(cat(i)%pvmax-pvsum)*weight
  end do
end subroutine estimate_nofz

!Make 1/Vmax,eff estimate of the luminosity function, computing the
! Vmax,eff values in the process and also its derivative wrt to the
! luminosity evolution parameter u
!
subroutine estimate_lf(ncat,cat,nzbin,zbin,dzbin,delta,nlf,magbin,mag,lf,u)
  use datatypes
  use galaxy_datatype
  use histograms
  implicit none
  type(binned_z) :: zbin(:)
  type(binned_mag) :: mag(:)
  integer,intent(in) :: nzbin,nlf,ncat
  real,intent(in) :: delta(nzbin),magbin(nlf),dzbin,u
  real,intent(out) :: lf(nlf)
  type(galaxy),intent(inout) :: cat(ncat)
  real,parameter :: TOL=5.0e-03
  logical,parameter :: ngp=.true.
  integer :: i,iz
  real*8 :: pvsum,pvsumeff
  real,external :: P,kpluse,dkpluse_du
!
  do i=1,ncat
     pvsumeff=0.0  !accumulate the delta weighted and P-weighted volume bin by bin
     pvsum=0.0
     iz=1 ! accumulate the contribution from all bins fully below zmax
     do while (cat(i)%zmax.gt.zbin(iz)%z+0.5*dzbin)
        pvsumeff=pvsumeff+delta(iz)*zbin(iz)%pvol
        pvsum=pvsum+zbin(iz)%pvol
        iz=iz+1
     end do
     ! then add on the contribution from the final fractional bin
     cat(i)%pvmaxeff=pvsumeff +delta(iz)*(cat(i)%pvmax-pvsum)
     if ((pvsum-cat(i)%pvmax)/pvsum.gt.TOL) write(0,*) 'ERROR: significantly negative effective volume',&
   &cat(i)%pvmaxeff,cat(i)%pvmax,pvsum,cat(i)%zmax,zbin(iz-1)%z+0.5*dzbin
  end do
! Tabulate for each galaxy the rate of change of its absolute magnitude
! wrt to the change in the evolution parameter u.
  do i=1,ncat
     cat(i)%dabsmag_du= -dkpluse_du(cat(i)%z,u)
  end do   
  !Count the number of galaxies in each magnitude bin
  cat(:)%weight=1.0
  call wcichist_deriv(ncat,cat(:)%absmag,cat(:)%dabsmag_du,cat(:)%weight,nlf,magbin,mag(:)%n,mag(:)%dn_du)
  if (ngp) call whist(ncat,cat(:)%absmag,cat(:)%weight,nlf,magbin,mag(:)%n)
  !Estimate the LF using weights equal to 1/Vmax,eff
  !This can either be done with CIC (cloud-in-cell) weighting or NGP (nearest grid point).
  !In either case we use the same assignment for the galaxy count mag(:)%n above, although
  !its derivtaive mag(:)%dn_du is only well defined for the CIC weighting and so is always used.
  cat(:)%weight=1.0/cat(:)%pvmaxeff
  if (ngp) then
     call whist(ncat,cat(:)%absmag,cat(:)%weight,nlf,magbin,lf) 
  else
     call wcichist(ncat,cat(:)%absmag,cat(:)%weight,nlf,magbin,lf) 
  end if
end subroutine estimate_lf


! Use Brent's method to take a single step towards finding the
! maximum of the function yy(x), within the bounds xmin to xmax 
! where we supply vectors xx(n) and yy(n) of the previous i evaluations 
! of yy(x) and the routine returns as its result the next value of xx to try 
! 
! Adapted from Numerical Recipes routine brent()
real function brent_step(xmin,xmax,n,i,xx,yy,converged)
  implicit none
  logical,intent(out) :: converged
  integer,intent(in) :: n,i
  real,intent(in) :: xx(n),yy(n),xmin,xmax
  real,parameter :: CGOLD=0.3819660,ZEPS=1.0e-05,TOL=2.0e-03
  real,save :: a,b,d,e,etemp,fu,fv,fw,fx,p,q,r,tol1,tol2,u,v,w,x,xm
  logical :: parabolic
  if (i.lt.1) then
     stop 'brent_step(): must have i=3 on first call'
  else if (i.eq.1) then
     ! initialise on first call
     if (xmin.ge.xmax) stop 'brent_step(): must have xmax>xmin'
     a=xmin
     b=xmax
     v=xx(1)
     w=v
     x=v
     e=0.0
     fx=yy(1)
     fv=fx
     fw=fx
  else if (i.gt.n) then
     write(*,*) 'maximum number of iterations reached'
     brent_step=xx(n) !return previous value
     converged=.true.
     return
  else
     fu=yy(i)  !new function evaluation at previously determined point
     if(fu.ge.fx) then  !if new value sucessfully higher than previous best,
        if(u.ge.x) then !which is at point x, move in one of the bounds [a,b]
           a=x          !to the point x, keeping a<u<b.
        else            !x will then be replaced by u
           b=x
        end if
        v=w            !update v and fv to be the previous value of w,fw
        fv=fw
        w=x            !update w and fw to be the second highest point
        fw=fx
        x=u            !update x and fx to the new highest point
        fx=fu
     else  !if new value fu at u is not higher than previous higest x,fx
        if(u.lt.x) then ! then move bound in to u while keeping a<x<b
           a=u
        else
           b=u
        end if
        if(fu.ge.fw .or. w.eq.x) then
           ! here we don't update x,fx as that is still the highest point
           ! but we update w,fw (second highest point) if fu>fw or if w=x
           v=w   !update v and fv to be the previous value of w,fw
           fv=fw
           w=u   !update w and fw to be the new 2nd highest point
           fw=fu
        else if(fu.ge.fv .or. v.eq.x .or. v.eq.w) then
           ! Here fx and fw, the 1st and 2nd highest points remain unchanged,
           ! but we update the 3rd highest point fv if fu<fv or v=x or v=w
           v=u
           fv=fu
        end if
     end if
  end if
!
  xm=0.5*(a+b)  !midpoint of the two bounds
  tol1=tol*abs(x)+ZEPS !tolerance for convergence
  tol2=2.*tol1         !twice tolerance
  if(abs(x-xm).le.(tol2-0.5*(b-a))) then
     converged=.true.
     brent_step=xx(i)
     return
  end if
  if(abs(e).gt.tol1) then !construct trial parabolic fit
     parabolic=.true.
     r=(x-w)*(fx-fv)
     q=(x-v)*(fx-fw)
     p=(x-v)*q-(x-w)*r
     q=2.*(q-r)
     if(q.gt.0.) p=-p
     q=abs(q)
     etemp=e
     e=d
  else
     parabolic=.false.   
  end if
!  
  if (.not.parabolic .or. & !test if parabolic fit produces unacceptable result
       & abs(e).le.tol1.or.abs(p).ge.abs(.5*q*etemp) .or. &
       & p.le.q*(a-x).or.p.ge.q*(b-x)) then
     if(x.ge.xm) then ! take golden section step
        e=a-x
     else
        e=b-x
     end if
     d=CGOLD*e
  else   
     d=p/q   !accept the parabolic fit
     u=x+d
!    If the result is going to be within tol2 of either bound
!    then instead make the shift tol1 with sign xm-x
     if(u-a.lt.tol2 .or. b-u.lt.tol2) d=sign(tol1,xm-x)
  end if
  if(abs(d).ge.tol1) then !use this to get the update x position u
     u=x+d
  else
     u=x+sign(tol1,d)
  endif
  brent_step=u !return u the trial new x value
  converged=.false.
  return
end function brent_step



end module rancat_subs




!-----------------------------------------------------------------------------
!  The kpluse(z)=(k(z)+e(z))-(k(zref)+e(zref))
!  I.E. the combined k+e correction defined to be zero at the reference redshift zref
!
! The relation between apparent, mag, and absolute magnitude, absmag is
!
! mag = absmag + 5log10(dlum/10pc) + kpluse(z)
!
! where kpluse(z) = -2.5 log (1+z) - 2.5 log (L_nu_0(1+z)^z=0/L_nu_0)
!                                  - 2.5 log (L_nu_0(1+z)^z=z/L_nu_0(1+z)^z=0)
!      Hogg et al (2002) (on astro-ph only)                    
!
real function kpluse(z,u,zref)
  implicit none
  real, intent(in) :: z,u,zref
  real :: ecorr,kcorr
  kpluse= kcorr(z)+ecorr(z,u) -( kcorr(zref)+ecorr(zref,u) )
  return
end function kpluse
!
! Differential of the above function with respect to
! the luminosity evolution parameter u
real function dkpluse_du(z,u)
  implicit none
  real, intent(in) :: z,u
  real :: decorr_du,ecorr
  integer, save :: ifirst=1
  real :: zt,ut,eps=1.0e-04,err
!
! On the first call only do a quick check to catch cases where the
! user has failed to supply consistent ecorr and decorr_du
  if (ifirst.eq.1) then
     do zt=0.0,2.0,0.1
        do ut=-2.0,2.0,0.1
           err=abs((ecorr(zt,ut+eps)-ecorr(zt,ut))/eps-decorr_du(zt,ut))
           if (err.gt.0.01) then
              write(0,*) 'Warning: supplied decorr_du not consistent with numerical derivative'
           end if   
        end do   
     end do   
     ifirst=0
  end if
!
! On subsequent calls we just use the user supplied function.
  dkpluse_du=decorr_du(z,u)
end function dkpluse_du
