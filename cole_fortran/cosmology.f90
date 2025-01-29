!   Module containing the table of redshifts, comoving and luminosity 
!   distances and the evolution weighted volume integral int P(z) dV/dz dz.
!   This module is used by both tabulate_distance() and
!   search_distance(). tabulate_distance() loads up the table and
!   search_distance() searches and interpolates.
!
module distancemagpvol_table
  implicit none
  integer, parameter:: NT=1000
  real, save:: rcomov(NT),z(NT),dl(NT),pv(NT),drc=0
end module distancemagpvol_table


module cosmology
contains

! Tabulate redshift, z, comoving distance, r_c, luminosity distance, d_L and
! pvol= int P(z) dV/dz dz
! for a chosen combination of omega0, lambda0 and volume evolution parameter
! pvol
!     
subroutine  tabulate_distance(omega0,lambda0,Pparam,zref)
  use distancemagpvol_table  !contains arrays in which tables are stored
  implicit none
  real,intent(in) :: omega0,lambda0,Pparam,zref
  integer :: it,NNT,jt
  real, parameter:: EPS=1.0e-04
  real, allocatable :: zplus1t(:),rct(:)
  real :: zplus1,x,om,atanhom,rchar,dz,intp,intt,h 
  real :: rcuml
  real :: P
!
!     Check that the omega0, lambda0, Pparam combination is one that this
!     program can deal with and
!     in the case of the open universe or one with a cosmological
!     constant compute some variables or a table of values
!     that will be useful later on
  if (omega0.lt.0.0) then    !flags non-relativistic assumption
  else if (abs(omega0-1.0).le.EPS .and. (lambda0).le.EPS) then !omega=1
  else if (abs(omega0-1.0).gt.EPS .and. (lambda0).le.EPS) then !open
     om=sqrt(1.0-omega0)  !useful numbers in open model
     atanhom=0.5*log((1.0+om)/(1.0-om) )
  else if (abs(omega0+lambda0-1.0).lt.EPS ) then !flat with lambda
     !        First tabulate rcomov with redshift by integrating
     NNT=10000
     dz=2.0/real(NNT)
     allocate(zplus1t(NNT),rct(NNT))
     rcuml=0.0
     rct(1)=0.0
     zplus1t(1)=1.0
     intp=0.0
     do it=2,NNT
        zplus1t(it)=1.0+real(it-1)*dz
        intt=1.0/sqrt(omega0*zplus1t(it)**3 + 1.0 - omega0)
        rcuml=rcuml+0.5*(intt+intp)*dz
        rct(it)=rcuml*3000.0
        intp=intt
     end do
  else
     write(0,*) 'ERROR: Omega0=',omega0,' Lambda0=',lambda0
     stop 'Not programmed for this cosmology'
  end if


!     Now build up the tables using constant steps in comoving
!     distance
!
  drc=3000.0/real(NT)
  do it=1,NT   !make table uniformly spaced in comoving distance.
     rcomov(it)=real(it-1)*drc  !comoving distance
     if (omega0.lt.0.0) then          !non-relavistic approximation
        zplus1=1.0 +rcomov(it)/3000.0
        z(it)= rcomov(it)/3000.0
        dl(it)=rcomov(it)
     else if (abs(omega0-1.0).le.EPS .and. (lambda0).le.EPS) then !omega=1
        zplus1=1.0/(1.0-rcomov(it)/6000.0)**2
        z(it)=zplus1-1.0
        dl(it)=rcomov(it)*zplus1
     else if (abs(omega0-1.0).gt.EPS .and. (lambda0).le.EPS) then !open
        x=tanh(atanhom-om*rcomov(it)/6000.0)
        zplus1=((om/x)**2-om**2)/omega0
        z(it)=zplus1-1.0
        rchar=(6000.0/(zplus1*omega0**2)) * &
             &           (omega0*z(it)+(omega0-2.0)*(sqrt(1+omega0*z(it))-1.0))
        dl(it)=rchar*zplus1
     else if (abs(omega0+lambda0-1.0).lt.EPS ) then !flat
        !look up redshift corresponding to rcomov(it) from temporary table
        jt=1
        do while (rcomov(it).ge.rct(jt)) 
           jt=jt+1
           if (jt.gt.NNT) stop &
           &   'tabulate_dist(): rct() not tabulated to sufficient redshift'
        end do
        h=(rcomov(it)-rct(jt-1))/(rct(jt)-rct(jt-1))
        zplus1=zplus1t(jt-1)*(1.0-h) + zplus1t(jt)*h
        z(it)=zplus1-1.0
        dl(it)=rcomov(it)*zplus1
        if (it.eq.1) then
           pv(1)=0.0
        else
          !look up redshift corresponding to rcomov(it)-0.5drc from temporary table
           jt=1
           do while (rcomov(it)-0.5*drc.ge.rct(jt)) 
              jt=jt+1
              if (jt.gt.NNT) stop &
              &  'tabulate_dist(): rct() not tabulated to sufficient redshift'
           end do
           h=(rcomov(it)-0.5*drc-rct(jt-1))/(rct(jt)-rct(jt-1))
           zplus1=zplus1t(jt-1)*(1.0-h) + zplus1t(jt)*h
           pv(it)=pv(it-1) + &
           & (P(zplus1-1.0,Pparam)/P(zref,Pparam))*(rcomov(it)**3-rcomov(it-1)**3)/3.0
        end if
     end if
  end do
!
!    Deallocate temporary arrays in the case that they were used
if (omega0.lt.0.0) then    !flags non-relativistic assumption
else if (abs(omega0-1.0).le.EPS .and. (lambda0).le.EPS) then !omega=1
else if (abs(omega0-1.0).gt.EPS .and. (lambda0).le.EPS) then !open
else if (abs(omega0+lambda0-1.0).lt.EPS ) then
   deallocate(zplus1t,rct)      
end if
return
end subroutine tabulate_distance

!-----------------------------------------------------------------------------
! Given either rcomov, z, dlum, Pvol or apparent magnitude mag
! look up the corresponding values of the other quantities. 
! The integer isel specifies which argument is set on input.
!
! The parameter absmag gives the absolue magnitude of the object
! and the externally provided kpluse function the k+e correction
! as a function of redshift, z.
!
! The relation between apparent, mag, and absolute magnitude, absmag is
!
! mag = absmag + 5log10(dlum/10pc) + kpluse(z)
!
! where kpluse(z) = -2.5 log (1+z) - 2.5 log (L_nu_0(1+z)^z=0/L_nu_0)
!                                  - 2.5 log (L_nu_0(1+z)^z=z/L_nu_0(1+z)^z=0)
!      Hogg et al (2002) (on astro-ph only)              
!
!  The quantity pvol is the evolution weighted comoving volumes
!   \int P(z) dV/dz dz where P(z) is an externally provided function
!  parameterizing the assumed evolution of \Phi^*.      

subroutine  search_distancemag(rc,zplus1,dlum,pvol,mag,isel,absmag,kpluse,eparam,zref)
  use distancemagpvol_table
  implicit none
  real,intent(inout) :: pvol,rc,zplus1,dlum
  real,intent(in) :: eparam,zref
  integer,intent(in) :: isel
  integer, save:: itlo=1,ithi=1
  integer:: it
  real :: h,rit,zt,magt,magtlo,magthi,absmag,mag
  real,parameter :: epsilon=1.0e-10
  real, external :: kpluse
!  First test drc, the table spacing, has been set 
!  which indicates tabulate_distance() has been called
  if (drc.eq.0) stop 'ERROR: tabulate_distance must be called before search_dist'
!  Locate the two elements in the table spanning the required
!  value and then linearly interpolate.
  if (isel.eq.1) then  !rc set so lookup zplus1, and dlum 
!      This is the fast/easy case as the table is equally spaced in rc 
     rit=1.0+rc/drc
     it = int(rit)
     h= rit-real(it) 
     if (it.ge.NT) stop 'search_dist() r beyond tabulated range'
     zplus1=1.0+ z(it)*(1.0-h)    +    z(it+1)*h
     dlum=       dl(it)*(1.0-h)  +    dl(it+1)*h
     pvol=       pv(it)*(1.0-h)  +    pv(it+1)*h
     mag = absmag + 5.0*log10(max(dlum,epsilon))+25.0 + kpluse(zplus1-1.0,eparam,zref)
  else if (isel.eq.2) then !zplus1 set 
!    Search for corresponding redshift
     zt=zplus1-1.0
     if (z(ithi).lt.zt .or. z(itlo).gt.zt ) then!short cut if zt close to
        itlo=1                                  !last call
        ithi=NT                                 !otherwise do binary search
        do while (ithi-itlo.gt.1) 
           it=(ithi+itlo)/2
           if(z(it).gt.zt)then
              ithi=it
           else
              itlo=it
           endif
        end do
     end if
     h=(zt-z(itlo))/(z(ithi)-z(itlo))
     rc=     rcomov(itlo)*(1.0-h)  +rcomov(ithi)*h
     dlum=       dl(itlo)*(1.0-h)  +    dl(ithi)*h
     pvol=       pv(itlo)*(1.0-h)  +    pv(ithi)*h
     mag = absmag + 5.0*log10(max(dlum,epsilon))+25.0 + kpluse(zplus1-1.0,eparam,zref)
!
  else if (isel.eq.3) then !luminosity distance set 
!    Search for corresponding luminosity distance
     if (dl(ithi).lt.dlum .or. dl(itlo).gt.dlum ) then
        itlo=1                                  
        ithi=NT                                 
        do while (ithi-itlo.gt.1) 
           it=(ithi+itlo)/2
           if(dl(it).gt.dlum)then
              ithi=it
           else
              itlo=it
           endif
        end do
     end if
     h=(dlum-dl(itlo))/(dl(ithi)-dl(itlo))
     rc=     rcomov(itlo)*(1.0-h)  +rcomov(ithi)*h
     zplus1=1.0+ z(itlo)*(1.0-h)    +    z(ithi)*h
     pvol=       pv(itlo)*(1.0-h)  +    pv(ithi)*h
     mag = absmag + 5.0*log10(dlum)+25.0 + kpluse(zplus1-1.0,eparam,zref)
     
  else if (isel.eq.4) then !pvol set
 !   Search for corresponding pvol
     if (pv(ithi).lt.pvol .or. pv(itlo).gt.pvol ) then
        itlo=1                                  
        ithi=NT                                 
        do while (ithi-itlo.gt.1) 
           it=(ithi+itlo)/2
           if(pv(it).gt.pvol)then
              ithi=it
           else
              itlo=it
           endif
        end do
     end if
     h=(pvol-pv(itlo))/(pv(ithi)-pv(itlo))
     rc=     rcomov(itlo)*(1.0-h)  +rcomov(ithi)*h
     zplus1=1.0+ z(itlo)*(1.0-h)    +    z(ithi)*h
     pvol=       pv(itlo)*(1.0-h)  +    pv(ithi)*h
     dlum=       dl(itlo)*(1.0-h)  +    dl(ithi)*h
     mag = absmag + 5.0*log10(dlum)+25.0 + kpluse(zplus1-1.0,eparam,zref)
!     
  else if (isel.eq.5) then !apparent magnitude set 
!    Binary search for corresponding apparent magnitude
     itlo=1                                  
     ithi=NT                                 
     do while (ithi-itlo.gt.1) 
        it=(ithi+itlo)/2
        magt= absmag + 5.0*log10(dl(it))+25.0 + kpluse(z(it),eparam,zref)
        if(magt.gt.mag)then
           ithi=it
        else
           itlo=it
        endif
     end do
     if (itlo.gt.1) then  !won't have log(0) problems
        magtlo=absmag + 5.0*log10(dl(itlo))+25.0 + kpluse(z(itlo),eparam,zref)
        magthi=absmag + 5.0*log10(dl(ithi))+25.0 + kpluse(z(ithi),eparam,zref)
        h=(mag-magtlo)/(magthi-magtlo)
        rc=     rcomov(itlo)*(1.0-h)  +rcomov(ithi)*h
        zplus1=1.0+ z(itlo)*(1.0-h)    +    z(ithi)*h
        pvol=       pv(itlo)*(1.0-h)  +    pv(ithi)*h
        dlum=       dl(itlo)*(1.0-h)  +    dl(ithi)*h
     else !if we end up here it is because the sought after redshift is between
          !0 and the second bin x(ithi). To avoid taking the log(zero) we take the
          !alternative interpolation strategy which is that the distances rc and dlum
          !are proporitonal to redshift and volume, pvol, tot he cubre.
        !first compute the redshift (stored here as h)
        h=z(ithi)/dl(ithi)*10.0**(mag-absmag-kpluse(z(ithi),eparam,zref)-25.0)/5.0
        rc= rcomov(ithi)*h/z(ithi)  !compute the returned values
        zplus1=1.0+h
        pvol=pv(ithi)*(h/z(ithi))**3
        dlum= dl(ithi)*h/z(ithi)
     end if   
  end if
  return
end subroutine search_distancemag


end module cosmology
