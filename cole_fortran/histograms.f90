module histograms
contains
! Construct a weighted histogram, assuming uniformly spaced bins 
!
! end bins extend to +/- infinity
!
subroutine whist(ndata,values,weights,nbin,bins,hist)
  implicit none
  integer,intent(in) :: ndata,nbin 
  real,intent(in) :: values(ndata),weights(ndata),bins(nbin)
  real,intent(out) :: hist(nbin)
  real :: inv_dbin
  integer :: i,ibin
! Compute inverse of bin spacing
  inv_dbin=(real(nbin)-1.0)/(bins(nbin)-bins(1))
! Initialize histogram  
  hist(:)=0.0
! Accumulate weighted histogram
  do i=1,ndata
     ibin=int(1.5+(values(i)-bins(1))*inv_dbin)
!    If outside binned range then stuff into first or last bin
     if (ibin.le.0) ibin=1
     if (ibin.gt.nbin) ibin=nbin
     hist(ibin)=hist(ibin)+weights(i)
  end do   
end subroutine whist




! Construct a weighted histogram using cloud-in-cell assignment and  
! assuming uniformly spaced bins 
!
! end bins extend to +/- infinity
!
subroutine wcichist(ndata,values,weights,nbin,bins,hist)
  implicit none
  integer,intent(in) :: ndata,nbin
  real,intent(in) :: values(ndata),weights(ndata),bins(nbin)
  real,intent(out) :: hist(nbin)
  integer :: i,ibin,ileft,iright
  real :: wleft,inv_dbin,ri
! Compute inverse of bin spacing
  inv_dbin=(real(nbin)-1.0)/(bins(nbin)-bins(1))
! Initialize histogram  
  hist(:)=0.0
! Accumulate weighted histogram
  do i=1,ndata
     ri=1.0+(values(i)-bins(1))*inv_dbin
     ibin=int(ri+0.5)  !index of nearest bin 
     if (ri.le.real(ibin)) then !index of bin to left
        ileft=ibin-1
     else
        ileft =ibin
     end if
     iright=ileft+1 !index of bin to right
     wleft=real(iright)-ri !overlap with bin to left
!    If outside binned range then stuff into first or last bin
     if (ileft.le.0) then
        ileft=1
        iright=1
     end if
     if (iright.gt.nbin) then 
        iright=nbin
        ileft=nbin
     end if
!    split the weight between these left and right bins using
!    CIC assignment
     hist(ileft)=hist(ileft)+weights(i)*wleft
     hist(iright)=hist(iright)+weights(i)*(1.0-wleft)
  end do   
end subroutine wcichist



! Construct a weighted histogram using cloud-in-cell assignment and  
! assuming uniformly spaced bins. At the same time construct the
! rate of change of this histogram wrt x given the rate of change wrt x
! of each data value.
!
! dhist_dx is an estimate of the rate of change of hist(values) in the limit
! of infinitesimal changes in dx, resulting in infinitesimal changes 
! dvalues_dx*dx in the values.
!
! Thus for dvalues_dx which are all small compared to the bin
! width then to a good approximation^**
!   dhist_dx = hist(values+dvalues_dx)-hist(values)  (Eq 1).
! If the dvalues_dx are scaled up by some factor, alpha, such that they
! are large then the dhist_dx will just increase in proportion to alpha
! rather than be given by the equation (Eq 1). This is the desired behaviour
! for a differential.
!
! ** For small but finite displacements, dvalue_dx, (Eq 1) will be satisfied 
! for all data values which when displaced by dvalue_dx still only 
! contribute to the same two bins of the histogram. There is trap in the code
! to catch the rare case of a data value whose weight is initially 100% in
! a single bin so that no matter which way it is displaced the correct two
! bin weights are adjusted. However the (Eq 1) will not be satisfied for 
! data points which say just straddle a bin (eg 1% in the lefthand bin
! and 99% in the righthand bin) for which the displacement is large enough
! to change which bin boundary they straddle. 
!
!
! end bins extend to +/- infinity
!
subroutine wcichist_deriv(ndata,values,dvalues_dx,weights,nbin,bins,hist,dhist_dx)
  implicit none
  integer,intent(in) :: ndata,nbin
  real,intent(in) :: values(ndata),weights(ndata),dvalues_dx(ndata),bins(nbin)
  real,intent(out) :: hist(nbin),dhist_dx(nbin)
  integer :: i,ibin,ileft,iright
  real :: wleft,inv_dbin,ri,dri_dx,dwleft_dx
! Compute inverse of bin spacing
  inv_dbin=(real(nbin)-1.0)/(bins(nbin)-bins(1))
! Initialize histograms  
  hist(:)=0.0
  dhist_dx(:)=0.0
! Accumulate weighted histograms
  do i=1,ndata
     ri=1.0+(values(i)-bins(1))*inv_dbin
     dri_dx=dvalues_dx(i)*inv_dbin
     dwleft_dx=-dri_dx !movement to the right decrease the weight given to the left bin
     ibin=int(ri+0.5)  !index of nearest bin 
     if (ri.le.real(ibin)) then !index of bin to left
        ileft=ibin-1
     else
        ileft =ibin
     end if
     iright=ileft+1 !index of bin to right
     wleft=real(iright)-ri !overlap with bin to left
!    If outside binned range then stuff into first or last bin
     if (ileft.le.0) then
        ileft=1
        iright=1
     end if
     if (iright.gt.nbin) then 
        iright=nbin
        ileft=nbin
     end if
!
!    Split the weight between these left and right bins using
!    CIC assignment
     hist(ileft)=hist(ileft)+weights(i)*wleft
     hist(iright)=hist(iright)+weights(i)*(1.0-wleft)
!
!    We then record the change in weight given to the bins if the value
!    of this data point is shifted to the right by dri_dx
!    The following if clause traps the case where 100% of the weight 
!    was given to the righthand bin and so the weights that change are
!    those of the righthand bin and the bin to the right of the righthand bin.
     if (wleft.eq.0.0 .and. dwleft_dx.lt.0) then !can't decrease further
        ileft=ileft+1
        iright=iright+1
        if (iright.gt.nbin) iright=nbin
     end if   
     dhist_dx(ileft)=dhist_dx(ileft)+weights(i)*dwleft_dx
     dhist_dx(iright)=dhist_dx(iright)+weights(i)*(-dwleft_dx)
  end do   
end subroutine wcichist_deriv
end module histograms
