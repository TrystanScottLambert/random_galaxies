module jswml
contains
! The Joint StepWise Maximum Likelihood algorithm described in Cole (2011)
!
! This code can actually use one of two algorithms to search for the same maximum likelihood
! soloution.
!    search_for_max=.true.
!       Here the method is as described in the paper. The luminosity evolution parameter
!       u is held fixed and the iterative procedure used to find the ML solution. The
!       posterior probablity is recorded and  then u is varied and the procedure repeated.
!       To determine the value of u one searches over u for the maximum probability. This
!       search is preformed using Brent's method rather than a grid search.         
!
!    search_for_max=.false. 
!       Here the code makes use of the extra constrain equation (ucon in the code) that
!       arises from differentiating the likelihood (strictly speaking the posterior 
!       probability) with respect to the luminosity evolution parameter u. A slightly
!       more complex iterative procedure is used to find the simultaneous solution of this
!       equation along with the original set. In some ways this is more elegant but its
!       downside is one has to explicitly bin and evaluate the LF function rendering the
!       result somewhat senstive the bin size and binning algorithm.  
!
!
subroutine rancat_jswml(ncat,cat,nmult,iseed,nran,rancat,nlf,magbin,lf,nzbin,zbins,delta, &
     & zref,survey,omega0,lambda0,prior,par_it,nitermax,search_for_max,iv)
  use datatypes
  use galaxy_datatype
  use rancat_subs
  use distancemagpvol_table
  use zbrent
  use histograms
  use cosmology
  implicit none
  integer,intent(in) :: iv !verbosity
  real,intent(in) :: zref
  type(survey_spec),intent(in) :: survey 
  type(prior_parameters),intent(in) :: prior
  integer, intent(in) :: nmult,nlf,nzbin,ncat
  logical,intent(in) :: search_for_max
  integer, intent(inout) :: nitermax,iseed
  real,intent(inout) :: magbin(nlf)
  type(galaxy), intent(inout) :: cat(ncat)
  real,intent(in) :: omega0,lambda0
  type(galaxy),allocatable,intent(out) :: rancat(:)
  type(model_parameters),intent(out),allocatable :: par_it(:)
  type(model_parameters) :: par,par_last,par_prev
  real,intent(out) :: lf(nlf)
  integer,intent(out) :: nran
  real,intent(out) :: delta(nzbin),zbins(nzbin)
  ! Various tolerance parameters some of which should only be adjusted with care!
  real, parameter :: tolM=1.0e-04 !used to check that magnitudes bins are uniformly spaced
  real, parameter :: tolL=1.0e-06,tolD=1.0e-4 !allowable error in the LF and Delta(z)
  real, parameter :: tolE=1.0e-03 !allowable error in the evolution parameters
  real, parameter :: tolR=1.0e-04 !allowable in evolution parameter when finding the
  !                               !root of ucon/acon. tolR should be less than tolE
  integer,parameter :: n_srch=1000!array size used in search over u (only search_for_max=.true.)
  type(binned_mag) :: mag(nlf)
  type(binned_z) :: zbin(nzbin)
  real :: dmag,dzbin
  real :: sigma(nzbin),delta_raw(nzbin)
  real :: lf_last(nlf)
  integer :: iter,nrancat,i,ibranch,jbranch,nitermin,i_srch
  real :: u_srch(n_srch),logpmax_srch(n_srch)
  real :: ntot,ntot_ran,ntot_model,ucons,acons
  logical :: converged,convergedPE,convergedu
!
! Part 1)   Various preparations.
!
! Set up the redshift bins uniformly between the specified min and max
  dzbin=(survey%zmax-survey%zmin)/real(nzbin)
  do i=1,nzbin
     zbin(i)%z=survey%zmin+(real(i)-0.5)*dzbin
  end do
  zbins(:)=zbin(:)%z  !copy to the redshift bins to the output array
!
! Check magnitude bins are uniform and record the spacing
  dmag=(magbin(nlf)-magbin(1))/real(nlf-1)
  do i=1,nlf
     if (abs(magbin(i)-(magbin(1)+real(i-1)*dmag)).gt.tolM*dmag) then
        write(0,*) 'Warning: Magnitude bins reset to have uniform spacing'
        magbin(i)=magbin(1)+real(i-1)*dmag
        write(0,*) 'bin ',i,' magnitude ',magbin(i)
     end if
  end do 
!  
! Tabulate the distance redshift relation including the volume weighted by the
! P(z,a)-evolution model of density evolution
  if (iv.ge.3) write(0,*) 'tabulating distance redshift relation'
  par%a=0.0 !initially assume the default density evolution
  call tabulate_distance(omega0,lambda0,par%a,zref)
!
  if (iv.ge.3) write(0,*) 'Tabulating P-evolution weighted volume of each redshift bin'
  call tabulate_binpvolumes(survey,nzbin,zbin(:)%z,zbin(:)%pvol,zbin(:)%vol,zref)
!
! Compute the nofz histogram
  if (iv.ge.3) write(0,*) 'Computing N(z)'
  cat(:)%weight=1.0 !for an unweighted histogram of nofz
  call wcichist(ncat,cat(:)%z,cat(:)%weight,nzbin,zbin(:)%z,zbin(:)%n)
!
! Allocate the array for output catalogue
  nrancat=int(1+real(ncat*nmult)*1.1)
  allocate(rancat(nrancat),par_it(nitermax))
!
!  Part 2:  Start iterations
!
! Set the minimum number of iterrations. If nitermax=1 then output will
! be one iteration and the standard Vmax result
  par%u=0.0  !first guess luminosity evolution parameter u
  par%a=0.0  !first guess density evolution parameter a
  par%mu=0.0 !first guess normalization/Lagrange multiplier parameter mu
  par_prev%u=1.0e+30  !when iterations have started these
  par_prev%a=1.0e+30  !three variables will contain the 
  par_prev%mu=1.0e+30 !parameter values from the previous iteration
  par_last=par !this will store the parameter values from several iterations
               !back to the point where we last switched branches (see Note BRANCH)
  lf(:)=0.0 !initialize the LF
  lf_last(:)=0.0 !this will store LF from previous iteration
  delta(:)=1.0 !initialize overdensity in each redshift bin
!
  convergedPE=.false. !changes to true when the evolution parameters (u and a) converge
  converged=.false.   !subsequently changes to true when lf(:) and delta(:) have converged
!                     !to high accuracy
  convergedu=.false. !changes to true when the search over u (search_for_max=.true. case) completes
  ibranch=0 !used later to indicate if last iteration was on branch 1 or 2 (see Note BRANCH)
  u_srch(1)=par%u !used in search over u when search_for_max=.true.
  i_srch=0
  iter=1
do                 !outer loop used when search_for_max=.true. to search for maximum 
  i_srch=i_srch+1  !posterior probability value of the luminosity evolution parameter u
  if (i_srch.gt.n_srch) stop 'ERROR: rancat_jswml() n_srch too small'
! loop that iterates to find the converged solution to the set of constraint equations
  nitermin=min(2+iter,nitermax) !unless Vmax is what you intend then do at least 3 iterations ????
  do while(.not.converged .and. iter.le.nitermax)
     !After a minimum of nitermin iterations check for convergence 
     !First check that the evolution parameters have converged
     if (iter.gt.nitermin .and. abs(par%a-par_last%a).le.tolE &
          & .and. abs(par%u-par_last%u).le.tolE ) then
        if (.not. convergedPE) then
           if(iv.ge.1) write(0,*) 'Density and luminosity evolution parameters converged to accuracy: ',&
                &' err(a)=',abs(par%a-par_last%a),' err(u)=',abs(par%u-par_last%u)
        end if
        convergedPE=.true.
     end if   
     !Then continue iterating just the LF and density bins until they have also 
     !accurately converged
     if (convergedPE .and. maxval(abs(delta(:)-zbin(:)%delta)).lt.tolD .and. &
            maxval(abs(lf(:)-lf_last(:))).lt.tolL )  then
           converged=.true.
           if (iv.ge.1)then
              write(0,*) 'Density bins converged to accuracy: ',maxval(abs(delta(:)-zbin(:)%delta))
              write(0,*) 'LF bins converged to accuracy: ',maxval(abs(lf(:)-lf_last(:)))
           end if
     end if
     !Retabulate to update pvol(z) for new value of par%a before the call to derived props
     call tabulate_distance(omega0,lambda0,par%a,zref)
     call tabulate_binpvolumes(survey,nzbin,zbin(:)%z,zbin(:)%pvol,zbin(:)%vol,zref)
     if (iv.ge.3) write(0,*) 'retabulated distance relation and P-weighted volbins'
     !
     !Recompute and update the derived galaxy properties for the new evolution parameter
     if (iv.ge.3) write(0,*) 'Tabulating derived galaxy properties for each catalogue galaxy'
     call derived_gal_props(par%u,zref,survey,ncat,cat)
     !
     !Make 1/Vmax,eff estimate of the luminosity function
     if (iv.ge.3) write(0,*) 'Estimating Vmaxeff and LF for iteration',iter
     lf_last(:)=lf(:)
     call estimate_lf(ncat,cat,nzbin,zbin,dzbin,delta,nlf,magbin,mag(:),lf,par%u)
     !
     ! Solve for mu by normalizing the average weight, <pvmax/(pvmaxeff+mu.pvmax)>, to unity
     if (iv.ge.3) write(0,*) 'Normalizing weights by solving for mu'
     call normalize_weights(cat,ncat,par%mu)
     !Record a history of evolution of various parameters and diagnostics for each iteration
     par%p_post=logp_post(par,prior,zbin,delta,nzbin,mag(:)%n,lf,nlf)
     par_it(iter)=par
     if (iv.ge.2) write(0,*) 'iter=',iter,' Log(P_post)=',par_it(iter)%p_post
     !
     !Estimate the redshift distribution that the random catalogue would have
     if (iv.ge.3) write(0,*) 'Estimating N(z) of random catalogue'
     call estimate_nofz(ncat,cat,par%mu,zbin,dzbin)
     ntot=sum(zbin(:)%n)
     ntot_ran=sum(zbin(:)%n_ran)
     if(iv.ge.3) write(0,*)'Renormalization factor for rancat:',ntot/ntot_ran,'(not applied)'
     !
     !Estimate overdensity, delta, in each redshift bin
     if (iv.ge.3) write(0,*) 'Estimating delta()'
     ! In the zbin(:)%delta we keep the old value of delta(:) prior to the
     ! following update. This is used to check convergence and as the input
     ! value when solving the acon and ucon constraint equations on this iteration. 
     zbin(:)%delta=delta(:) 
     call estimate_delta(nzbin,zbin(:),prior%fourPiJ3,delta,delta_raw,sigma)
     ntot_model=sum(zbin(:)%n*delta(:))
     if(iv.ge.3) write(0,*) 'ntot=',ntot,'   ntot_model=',ntot_model
     !Save catalogues for selected and the final converged iteration
     if (converged .or. iter.eq.nitermax) then
        if (iter.eq.nitermax) then
           write(0,*) 'ERROR: Failed  to reach convergence after ',iter, 'iterations'
           write(0,*) 'Exiting and saving result as if converged!!!'
           converged=.true.
        end if
!
!       Generate the random catalogue using the converged parameters        
        if (.not.search_for_max) then
           if(iv.ge.2) write(0,*) 'Generating random catalogue using pvmax/(pvmaxeff+mu.pvmax) weights'
           call generate_rancat(iseed,nmult,survey,par%mu,par%u,zref,ncat,cat,nrancat,rancat,nran)
        end if   
     end if
     !
     ! Note BRANCH:
     !      If we are fitting for both luminosity and density evolution (parameters u and a)
     !      then we take one of two branches, either updating a or updating u.
     !      Trying to update both while keeping delta() and lf() fixed leads to none convergence.
     !      Here the scheme is to update the parameter whose constraint equation,
     !      acon=0 or ucon=0, has the largest residual, except if doing so has not
     !      produced a significant change in the corresponding parameter.
     !
     if (.not.convergedPE) then
        !Solve for density or luminosity evolution parameter par%a or par%u
        !Branch 1 will update a 
        !Branch 2 will update u 
        !Here we make the decisions as to which branch to take:
        !
        !First we determine the amount by which each constraint is violated 
        ucons=ucon(par%u,prior%u,par%mu,dmag,survey%magfaint,nzbin,nlf,&
     &             lf,magbin,mag(:),zbin(:),zref,prior%spline)
        acons=acon(par%a,par%a,prior%a,par%mu,nzbin,zbin(:),prior%spline)
        if (prior%a.le.0.0 .and. (prior%u.le.0.0 .or. search_for_max)) then
           jbranch=0 !update neither a nor u as both are fixed.
        else if (prior%u.le.0.0 .or. search_for_max) then
           jbranch=1 !update a as u is fixed
        else if (prior%a.le.0.0) then
           jbranch=2 !update u as a is fixed
        else if (abs(acons).gt.abs(ucons) ) then !The acon equation has the biggest residual
           if( abs(par%a-par_prev%a).gt.tolE .or. ibranch.ne.1 ) then
              jbranch=1  !update a as acon residual greater 
           else
              jbranch=2  !update u instead as we tried branch 1 for the same reason
                         !on the last iteration and it didn't alter a by more than the tolerance
           end if
        else !The ucon equation has the biggest residual
           if ((abs(par%u-par_prev%u).gt.tolE .or. ibranch.ne.2) ) then
              jbranch=2 !update u as ucon residual greater 
           else
              jbranch=1 !update a instead as we tried branch 1 for the same reason
              !on the last iteration and it didn't alter u by more than the tolerance
           end if
        end if
        if (jbranch.eq.1) then
           if(iv.ge.2) write(0,*) 'Applying density evolution constraint (Branch 1)'
           par_prev=par !save parameter values before updating 
           if (ibranch.ne.1) then !we have just switched to this branch
              par_last%a=par%a  !and so we record the starting value of a
           end if
           par%a=zbrent_acon(-4.0,4.0,tolR,par%a,prior%a,par%mu,nzbin,zbin,prior%spline)
           ibranch=1
        else if (jbranch.eq.2) then
           if (iv.ge.2) write(0,*) 'Applying luminosity evolution constraint (Branch 2)'
           par_prev=par !save parameter values before updating 
           if (ibranch.ne.2) then !we have just switched to this branch
              par_last%u=par%u  !and so we record the starting value of par%u
           end if
           par%u=zbrent_ucon(-3.46,3.46,tolR,prior%u,par%mu,dmag,survey%magfaint,nzbin &
     &     ,nlf,lf,magbin,mag(:),zbin(:),zref,prior%spline)
           ibranch=2
        end if
!       If we are always only iterating one branch as either we have a delta-function
!       prior on the other parameter or we using the search method and keeping u fixed
!       during each iteration then we need to set par_last=par_prev for the parameter
!       being iterated
        if (search_for_max .or. prior%u.le.0.0) par_last=par_prev
        if (prior%a.le.0.0) par_last=par_prev
!
        if (iv.ge.1) then
           ! Report by how much the constraints are violated after this iteration
           ucons=ucon(par%u,prior%u,par%mu,dmag,survey%magfaint,nzbin,nlf,&
                &             lf,magbin,mag(:),zbin(:),zref,prior%spline)
           acons=acon(par%a,par%a,prior%a,par%mu,nzbin,zbin,prior%spline)

           write(0,'(a,i3,2(a,f7.3,a,f10.3))') 'iter=',iter,' a=',par%a,' acon=',acons,' u=',par%u,' ucon=',ucons
        end if
           
     end if
     !
     !next iteration
     iter=iter+1
  end do
  ! We only use this outer do loop if searching over u for the maximum probability solution
  if (.not.search_for_max .and. prior%u.gt.0.0)   exit 
  ! Record the maximum posterior probability acheived at the end of the above
  ! iterations over the number density evolution parameter, a, and the LF and radial
  ! density values.
  logpmax_srch(i_srch)=logp_post(par,prior,zbin,delta,nzbin,mag(:)%n,lf,nlf)
  ! Update u based on Brent's method of searching for the maximum 
  par%u=brent_step(-3.0*prior%u,3.0*prior%u,n_srch,i_srch,u_srch,logpmax_srch,convergedu)
  u_srch(i_srch+1)=par%u !record updated u parameter for next step of the search
  write(*,*) '*** updated u in outer loop:',u_srch(i_srch),'  to', u_srch(i_srch+1)
  if (convergedu) exit !exit this loop if search over u has converged
  ! Repeat the inner loop for this new value of u, until it has again converged
  converged=.false.
  convergedPE=.false.
end do
  if (iv.ge.0 .and. converged) write(0,*) 'rancat_jswml() exited sucessfully with converged results'
  if (iv.ge.0 .and. search_for_max) write(0,*) 'Searching over u for maximum took ',i_srch,' steps'
  if(converged) nitermax=iter-1 !return in nitermax the actual number of iterations taken
!       Generate the random catalogue using the converged parameters        
  if (search_for_max) then
     if(iv.ge.2) write(0,*) 'Generating random catalogue using pvmax/(pvmaxeff+mu.pvmax) weights'
     call generate_rancat(iseed,nmult,survey,par%mu,par%u,zref,ncat,cat,nrancat,rancat,nran)
  end if
end subroutine rancat_jswml
end module jswml
