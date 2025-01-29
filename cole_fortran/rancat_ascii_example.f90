! ASCII version of rancat_example. 
!   I.E. inputs and outputs ASCII data files rather than hdf5
!
!
! This example program reads an input catalogue with corresponding selection parameters, 
! finds the JSWML Luminosity Function (LF) and generates a corresponding random catalogue. 
! The random catalogue, the LF and various auxillary parameters, together
! with a copy of the original catalogue are written to a set of ascii files.
!
! To adapt it to your own purposes you will have to modify the code to read your own
! catalogue and supply the other parameters described below.
!
! In addition you will probably need to modify the code in the file user_routines.f90
! to supply appropriate k-corrections as described in README
!
  program make_rancat
  use datatypes   !module defining derived datatypes  (rancat_subs.f90)   
  use galaxy_datatype !module defining user specified galaxy datatype (user_routines.f90)
  use rancat_subs !module containing slave subroutines (rancat_subs.f90)   
  use jswml       !module containing main JSWML subroutine (rancat_swml.f90)
  implicit none
  real, parameter :: omega0=0.3,lambda0=0.7 !cosmological parameters
  type(prior_parameters) :: prior!composite variable carrying the prior parameters (see Note 1)
  type(survey_spec) :: survey !composite variable carrying the survey parameters (see Note 2)
  integer, parameter :: nmult=10 !factor by which you want random catalogue to be larger 
                                 !than the original input catalogue
  integer, parameter :: nzbin=40 !number of redshift bins used between the min and max 
                                 !redshift of the catalogue.  
  integer, parameter :: nlf=25  !number of absolute magnitude bins in the LF
  real, parameter :: dmag=0.5  !in this example the mag bins are set up with this spacing
  integer :: ncat,nitermax !input catalogue size and maximum number of iterations
  integer :: iseed=-5875 !random number seed used in generating the random catalogue
  logical,parameter :: search_for_max=.false. !set to .true. for the method descirbed in the paper
  integer :: verbosity ![-1 to 3] controls the amount info written to stnderr during execution
  real :: zref !reference redshift for k+e=0
  real :: magbin(nlf),lf(nlf) !magnitude bins and LF array
  real :: zbin(nzbin),delta(nzbin) !arrays in which the overdensity in redshift bins returned
  type(galaxy),allocatable :: cat(:),rancat(:) !genuine and random catalogue arrays
  type(model_parameters),allocatable :: par_it(:) !stores parameter iteration history
  character*100 :: catfile,rancatfile !input and output file names
  character*1 :: hash  !used to read ascii data file
  character*40 :: label
  integer :: nran,i
!
! Note 1: Set the priors
  prior%fourPiJ3=5000.0 !The galaxy clustering parameter used in the density fluctuation prior
  prior%u=0.01 !sigma of the Gaussian prior on the luminosity evolution parameter 
               !Setting to zero will mean u is kept fixed at its input value
  prior%a=0.01 !sigma of the Gaussian prior on the density evolution parameter
               !Setting to zero will mean a is kept fixed at its input value
  prior%spline=.true. !If true use a B-spline approximation to a Gaussian rather than
               !a pure Gaussian. The B-spline has zero probability beyond +/- 2sqrt(3) sigma.
  nitermax=400 !maximum number of iterations. (default should be something like 400)
               !Setting it to 1 will produce the 1/Vmax LF and corresponding random
               !catalogue.
               !On exit nitermax will equal the actual number of iterations taken to 
               !reach convergence
  zref=0.0 !reference redshift at which to define the absolute magnitudes, i.e. the
           !redshift at which k+e=0 by definition.
! Set absolute magnitude bins for the luminosity function. They should have uniform spacing.
! The code will complain and reset them if they are not uniformly spaced.
  do i=1,nlf
     magbin(i)=-25.0+real(i)*dmag
  end do
!  
!
! Read your data catalogue into the correct elements of the composite array cat(ncat).
! See user_routines.f90 for an explanation of its elements. The value of ncat
! should be the number of galaxies in your catalogue and array cat() should have
! dimension ncat. 
!    Note 2:  Also read the corresponding selection limits and solid angle of the survey
!    into the elements of the composite variable survey
  catfile="TestData/catalogue.txt"
  write(*,*) 'reading input catalogue from ',trim(catfile)
  open(10,file=catfile,form='formatted')
  read(10,*) hash,label !reads and ignores labels in the ASCII file
  read(10,*) hash,label,survey%zmin,label,survey%zmax
  write(*,'(2(a,f7.4))') 'zmin=',survey%zmin,' zmax=',survey%zmax
  read(10,*) hash,label,label,survey%magfaint
  write(*,'(a,f7.4)') 'mag_faintlimit=',survey%magfaint
  read(10,*) hash,label,label,survey%solid_angle,label,label
  write(*,'(a,f7.4,a)') 'solid angle= ',survey%solid_angle,'  square degrees'
  read(10,*) hash,label,label,ncat
  allocate(cat(ncat))
  read(10,*) hash,label,label
  do i=1,ncat
     read(10,*)  cat(i)%mag,cat(i)%z,cat(i)%id
  end do
  close(10) 
  write(*,*) 'catalogue read. ncat=',ncat,' maglimit=',survey%magfaint
  write(*,*)
!
! Iteratively estimate the luminosity function (LF) and corresponding random galaxy catalogue
! using the method of Cole 2011, which is based on a Joint StepWise Maximum Likelihood
! estimate of the LF and run of overdensity with redshift.
!
! This is the call the main subroutine
  verbosity=2 !default setting regading how much info is written to stnderr
  write(*,*) 'Calling rancat_swml() to iterate to the ML solution'
  call rancat_jswml(ncat,cat,nmult,iseed,nran,rancat,nlf,magbin,lf,nzbin,zbin,delta, &
     & zref,survey,omega0,lambda0,prior,par_it,nitermax,search_for_max,verbosity)
!
! On sucessful exit:
!  i) The random catalogue is contained in the array rancat(nran) with nran being the
!    number of galaxies in the catalogue, which should be close to ncat*nmult
!
!  ii) The LF is contained in lf(nlf) with corresponding magnitude bins magbin(nlf)
!
!  iii) nitermax will equal the number of iterations it took to reach convergence
!
!  iv) The composite array par_it(nitermax) contains the iteration history of
!      each of the model parameters (a,u and mu) along with (to within an arbitary
!      additive constant) the natural log of the posterior probability.
!
! As an example we write out all these results below
!
  rancatfile="TestData/rancat.txt"
  write(*,*) 'Saving results in ', trim(rancatfile)
  open(10,file=rancatfile,form='formatted')
  write(10,'(a)') '# Parameters:'
  write(10,'(2(a,f7.4))') '# Omega0= ',omega0,' Lambda0= ',lambda0
  write(10,'(2(a,f7.4))') '# zmin= ',survey%zmin,' zmax= ',survey%zmax 
  write(10,'(a,f7.4)') '# faint maglimit= ',survey%magfaint
  write(10,'(a,f7.4,a)') '# solid angle= ',survey%solid_angle,'  square degrees'
  write(10,'(a,i10)') '# Catalogue: nran= ',nran
  write(10,'(a)') '# mag    z         absmag     dm       v          vmax'
  do i=1,nran
     write(10,'(f7.4,1x,f7.4,1x,f10.4,1x,f8.4,1x,f11.3,1x,f11.3,1x,i7)') &
     rancat(i)%mag,rancat(i)%z,rancat(i)%absmag,rancat(i)%dm,rancat(i)%pv,rancat(i)%pvmax,rancat(i)%id
  end do   
  close(10)
!
  write(*,*) 'Saving iteration history'
  open(10,file='TestData/iterations.txt',form='formatted')
  write(10,'(a,i10)')  '# nitermax= ',nitermax
  write(10,'(a)') '# u              a              mu             P_post'
  do i=1,nitermax
     write(10,*) par_it(i)%u,par_it(i)%a,par_it(i)%mu,par_it(i)%p_post
  end do   
  close(10)
!
  write(*,*) 'Saving Luminosity Function'
  open(10,file='TestData/LF.txt',form='formatted')
  write(10,'(a,i10)')  '# nmag= ',nlf
  write(10,'(a)') '# mag            LF'
  do i=1,nlf
     write(10,*) magbin(i),lf(i)
  end do
  close(10)
!
  write(*,*) 'Saving Delta(z)'
  open(10,file='TestData/delta.txt',form='formatted')
  write(10,'(a,i10)')  '# nzbin= ',nzbin
  write(10,'(a)') '# z            delta'
  do i=1,nlf
     write(10,*) zbin(i),delta(i)
  end do
  close(10)
!
  write(*,*) 'Saving copy of input catalogue but with extra parameters, eg absmag'  
  open(10,file='TestData/catalogue_out.txt',form='formatted')
  write(10,'(a)') '# Selection criteria:'
  write(10,'(2(a,f7.4))') '# zmin= ',survey%zmin,' zmax= ',survey%zmax
  write(10,'(a,f7.4)') '# faint maglimit= ',survey%magfaint
  write(10,'(a,f7.4,a)') '# solid angle= ',survey%solid_angle,'  square degrees'
  write(10,'(a,i10)') '# Catalogue: ncat=',ncat
  write(10,'(a)') '# mag    z         absmag     dm       v          vmax     id'
  do i=1,ncat
     write(10,'(f7.4,1x,f7.4,1x,f10.4,1x,f8.4,1x,f11.3,1x,f11.3,1x,i7)') &
     & cat(i)%mag,cat(i)%z,cat(i)%absmag,cat(i)%dm,cat(i)%pv,cat(i)%pvmax,cat(i)%id
  end do
  close(10)
!
end program make_rancat

