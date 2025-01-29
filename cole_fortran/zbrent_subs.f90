module zbrent
contains
! Adapted from Numerical Recipes zbrent() routine.
! Finds the root, u, of the equation ucon(u|fixed params)=0
! The function ucon is defined in the module rancat_subs
! The root must be bounded by x1<u<x2 and is determined to
! an accuracy tol.
!
real function zbrent_ucon(x1,x2,tol,sigmaU,mu,dmag,faint_maglimit,nzbin &
  &     ,nlf,lf,magbin,mag,zbin,zref,spline)
  use datatypes
  use rancat_subs
  implicit none
  integer,intent(in) ::  nlf,nzbin
  type(binned_z),intent(in) ::  zbin(:)
  type(binned_mag),intent(in) ::  mag(:)
  real,intent(in)::  sigmaU,mu,dmag,x1,x2,tol,faint_maglimit,zref
  real,intent(in) :: lf(nlf),magbin(nlf)
  logical,intent(in) :: spline
  integer,parameter :: ITMAX=100
  real,parameter :: EPS=3.0e-08
  integer :: iter
  real :: a,b,c,d,e,fa,fb,fc,p,q,r,s,tol1,xm
  a=x1
  b=x2 
  fa=ucon(a,sigmaU,mu,dmag,faint_maglimit,nzbin,nlf,lf,magbin,mag(:),zbin(:),zref,spline)
  fb=ucon(b,sigmaU,mu,dmag,faint_maglimit,nzbin,nlf,lf,magbin,mag(:),zbin(:),zref,spline)
  if((fa.gt.0..and.fb.gt.0.).or.(fa.lt.0..and.fb.lt.0.))pause &
       & 'root must be bracketed for zbrent_ucon'
  c=b
  fc=fb
  do  iter=1,ITMAX
     if((fb.gt.0..and.fc.gt.0.).or.(fb.lt.0..and.fc.lt.0.))then
        c=a
        fc=fa
        d=b-a
        e=d
     endif
     if(abs(fc).lt.abs(fb)) then
        a=b
        b=c
        c=a
        fa=fb
        fb=fc
        fc=fa
     endif
     tol1=2.*EPS*abs(b)+0.5*tol
     xm=.5*(c-b)
     if(abs(xm).le.tol1 .or. fb.eq.0.)then
        zbrent_ucon=b
        return
     endif
     if(abs(e).ge.tol1 .and. abs(fa).gt.abs(fb)) then
        s=fb/fa
        if(a.eq.c) then
           p=2.*xm*s
           q=1.-s
        else
           q=fa/fc
           r=fb/fc
           p=s*(2.*xm*q*(q-r)-(b-a)*(r-1.))
           q=(q-1.)*(r-1.)*(s-1.)
        endif
        if(p.gt.0.) q=-q
        p=abs(p)
        if(2.*p .lt. min(3.*xm*q-abs(tol1*q),abs(e*q))) then
           e=d
           d=p/q
        else
           d=xm
           e=d
        endif
     else
        d=xm
        e=d
     endif
     a=b
     fa=fb
     if(abs(d) .gt. tol1) then
        b=b+d
     else
        b=b+sign(tol1,xm)
     endif
     fb=ucon(b,sigmaU,mu,dmag,faint_maglimit,nzbin,nlf,lf,magbin,mag,zbin,zref,spline)
  end do
  pause 'zbrent_con exceeding maximum iterations'
  zbrent_ucon=b
  return
end function zbrent_ucon

! Adapted from Numerical Recipes zbrent() routine.
! Finds the root, a, of the equation acon(a|fixed params)=0
! The function acon is defined in the module rancat_subs
! The root must be bounded by x1<a<x2 and is determined to
! an accuracy tol.
!
real function zbrent_acon(x1,x2,tol,qparam_ref,sigmaA,mu,nzbin,zbin,spline)
  use datatypes
  use rancat_subs
  implicit none
  integer,parameter :: ITMAX=100
  type(binned_z),intent(in) :: zbin(:)
  real,intent(in) :: x1,x2,tol
  logical,intent(in) :: spline
  real,parameter ::EPS=3.0e-08
  integer:: iter
  real:: a,b,c,d,e,fa,fb,fc,p,q,r,s,tol1,xm
  real,intent(in) :: qparam_ref,sigmaA,mu
  integer,intent(in) :: nzbin
  a=x1
  b=x2
  fa=acon(a,qparam_ref,sigmaA,mu,nzbin,zbin,spline)
  fb=acon(b,qparam_ref,sigmaA,mu,nzbin,zbin,spline)
  if((fa.gt.0..and.fb.gt.0.).or.(fa.lt.0..and.fb.lt.0.))pause &
       & 'root must be bracketed for zbrent'
  c=b
  fc=fb
  do iter=1,ITMAX
     if((fb.gt.0..and.fc.gt.0.).or.(fb.lt.0..and.fc.lt.0.))then
        c=a
        fc=fa
        d=b-a
        e=d
     endif
     if(abs(fc).lt.abs(fb)) then
        a=b
        b=c
        c=a
        fa=fb
        fb=fc
        fc=fa
     endif
     tol1=2.*EPS*abs(b)+0.5*tol
     xm=.5*(c-b)
     if(abs(xm).le.tol1 .or. fb.eq.0.)then
        zbrent_acon=b
        return
     endif
     if(abs(e).ge.tol1 .and. abs(fa).gt.abs(fb)) then
        s=fb/fa
        if(a.eq.c) then
           p=2.*xm*s
           q=1.-s
        else
           q=fa/fc
           r=fb/fc
           p=s*(2.*xm*q*(q-r)-(b-a)*(r-1.))
           q=(q-1.)*(r-1.)*(s-1.)
        endif
        if(p.gt.0.) q=-q
        p=abs(p)
        if(2.*p .lt. min(3.*xm*q-abs(tol1*q),abs(e*q))) then
           e=d
           d=p/q
        else
           d=xm
           e=d
        endif
     else
        d=xm
        e=d
     endif
     a=b
     fa=fb
     if(abs(d) .gt. tol1) then
        b=b+d
     else
        b=b+sign(tol1,xm)
     endif
     fb=acon(b,qparam_ref,sigmaA,mu,nzbin,zbin,spline)
  end do
  pause 'zbrent_acon exceeding maximum iterations'
  zbrent_acon=b
  return
END function zbrent_acon
end module zbrent
