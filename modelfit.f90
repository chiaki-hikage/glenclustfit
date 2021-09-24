
  integer :: npow
  integer, parameter :: npowlim=100
  double precision :: powlg_r(npowlim),powlg_obs(npowlim,nlrg),powlg_cov(npowlim,nlrg,nlrg),powlg_twohalo(npowlim),powlg_cov_inv(npowlim,nlrg,nlrg)
  double precision, parameter :: delv=200,m_fid=1.4d14,rvir_fid=(3*m_fid/(4*pi*delv*rhom))**(1./3.)
  double precision, parameter :: dchi=200
  double precision, parameter :: munit=1e14

contains

  subroutine read_powlg_obsdata()

    implicit none

    integer :: i,j,j2,n
    real :: xt,yt,et,covt,skip
    double precision :: xp(256),yp(256),yp2(256)
    character :: inf*100

    powlg_cov=0.
    do j=1,nlrg
       i=0
       inf=infname(1:len_trim(infname))//clrg(j)
       open(1,file=inf,status='old')
       do n=1,npowlim
          read(1,*,end=9) xt,yt,et
          xt=xt*1e-3  ! change unit from kpc/h to Mpc/h
          if (xt>rmin.and.xt<rmax) then
             i=i+1
             powlg_r(i)=xt
             powlg_obs(i,j)=yt
             powlg_cov(i,j,j)=et*et
          endif
       enddo
9      close(1)
       npow=i
    enddo

    if (ccov=='y') then
       do j=1,nlrg
          do j2=1,j-1
             inf=infname(1:len_trim(infname))//clrg(j)(1:len_trim(clrg(j)))//'_vs_'//clrg(j2)(1:len_trim(clrg(j2)))
             write(*,*) inf
             open(1,file=inf,status='old')
             i=0
             do n=1,npowlim
                read(1,*,end=8) xt,covt
                xt=xt*1e-3
                if (xt>rmin.and.xt<rmax) then
                   i=i+1
                   powlg_cov(i,j,j2)=covt*sqrt(powlg_cov(i,j,j)*powlg_cov(i,j2,j2))
                   powlg_cov(i,j2,j)=powlg_cov(i,j,j2)
                endif
             enddo
8            close(1)
             if (i/=npow) then
                write(*,*) 'strange',i,npow
                stop
             endif
          enddo
       enddo
    endif

    do i=1,npow
       call calcinvarray(powlg_cov(i,1:nlrg,1:nlrg),powlg_cov_inv(i,1:nlrg,1:nlrg),nlrg)
    enddo

    j=0
    open(1,file=inf2h,status='old')
    do n=1,128
       read(1,*) xt,yt
       if (xt>rmin.and.xt<rmax) then
          j=j+1
          xp(j)=xt
          yp(j)=yt
       endif
    enddo
    close(1)
    call spline(xp,yp,j,1d30,1d30,yp2)

    do i=1,npow
       call splint(xp,yp,yp2,j,powlg_r(i),powlg_twohalo(i))
    enddo

!!! when 2-halo terms are neglected,
!    powlg_twohalo=0.

  end subroutine read_powlg_obsdata
    
  subroutine masspro_singlemass(m,c,soff,fsub,y,qcen,ycomp)

    implicit none

    integer :: i,j
    double precision, intent(in) :: m,soff,fsub,c,qcen  ! m in unit of 10^14M_sun/h
    double precision, intent(out) :: y(1:npow),ycomp(1:npow,1:3)
    double precision :: rhalo,rs,mnfwc,ycen,yoff,xoff,fac,x
    integer, parameter :: nmax=256
    double precision :: logrmin,logrmax,dlogr
    double precision, dimension(nmax) :: xc,yc,yc2,ys,ys2

! factor
    fac=1e-12 ! for changinge the unit of scale from Mpc^2 to pc^-2

! halo size
    rhalo=(3*m*munit/(4*pi*delv*rhom))**(1./3.)
    
! halo core size
    rs=rhalo/c

! integral factor
    mnfwc=log(1+c)-c/(1+c)

! smoothing scale normalized by rs
    xoff=soff/rs

    logrmin=-4.d0
    logrmax=4.d0
    dlogr=(logrmax-logrmin)/nmax
    do i=1,nmax
       xc(i)=10**(logrmin+(i-0.5)*dlogr)
       yc(i)=nfw2(xc(i))
       ys(i)=yc(i)
    enddo

    call fftlog_smooth(xc,yc,nmax,2.d0,-4.d0,4.d0,0.d0)
    call spline(xc,yc,nmax,1d30,1d30,yc2)
    call fftlog_smooth(xc,ys,nmax,2.d0,-4.d0,4.d0,xoff)
    call spline(xc,ys,nmax,1d30,1d30,ys2)
    
    do i=1,npow
       x=powlg_r(i)/rs
       call splint(xc,yc,yc2,nmax,x,ycen)
       call splint(xc,ys,ys2,nmax,x,yoff)
       ycomp(i,1)=fac*m*munit*qcen*ycen/mnfwc/(2*pi*rs*rs)
       ycomp(i,2)=fac*m*munit*(1-qcen)*yoff/mnfwc/(2*pi*rs*rs)
       ycomp(i,3)=fac*m*fsub*munit/pi/powlg_r(i)**2
       y(i)=ycomp(i,1)+ycomp(i,2)+ycomp(i,3)
    enddo

  end subroutine masspro_singlemass

  subroutine masspro_mdist(camp,soff,fsub,y,qcen,ycomp,ispec)

!!! take into account the variation of halo mass

    implicit none

    integer :: i,j,im
    integer, parameter :: immax=30
    integer, intent(in) :: ispec
    double precision, intent(in) :: soff,fsub,camp,qcen 
    double precision, intent(out) :: y(1:npow),ycomp(1:npow,1:3)
    double precision :: rhalo,rs,mnfwc,ycen,yoff,xoff,fac,x,amp,damp,m,fprob,c,msub
    integer, parameter :: nmax=256
    double precision :: logrmin,logrmax,dlogr
    double precision, dimension(nmax) :: xc,yc,yc2,ys,ys2,xcin,ycin,ysin

    logrmin=-4.d0
    logrmax=4.d0
    dlogr=(logrmax-logrmin)/nmax
    do i=1,nmax
       xcin(i)=10**(logrmin+(i-0.5)*dlogr)
       ycin(i)=nfw2(xcin(i))
       ysin(i)=ycin(i)
    enddo

    xc=xcin
    yc=ycin
    call fftlog_smooth(xc,yc,nmax,2.d0,-4.d0,4.d0,0.d0)
    call spline(xc,yc,nmax,1d30,1d30,yc2)

! factor
    fac=1e-12 ! for changinge the unit of scale from Mpc^2 to pc^-2

    if (clrg(ispec)=='highpcen') then
       open(1,file='data/mdist_highpcen_cen.dat',status='old')
    else
       open(1,file='data/mdist_RMCGneBCG_cen.dat',status='old')
    endif

    ycomp=0.d0
    do im=1,immax
       
       read(1,*) m,c,fprob
       if (fprob<1e-5) cycle
       c=camp*c

! halo size
       rhalo=(3*m/(4*pi*delv*rhom))**(1./3.)

! halo core size
       rs=rhalo/c

! integral factor
       mnfwc=log(1+c)-c/(1+c)

! smoothing scale normalized by rs
       xoff=(soff*rhalo)/rs

       xc=xcin
       ys=ysin

       call fftlog_smooth(xc,ys,nmax,2.d0,-4.d0,4.d0,xoff)
       call spline(xc,ys,nmax,1d30,1d30,ys2)
    
       damp=fprob*fac*m/mnfwc/(2*pi*rs*rs)
       amp=amp+damp
       do i=1,npow
          x=powlg_r(i)/rs
          call splint(xc,yc,yc2,nmax,x,ycen)
          call splint(xc,ys,ys2,nmax,x,yoff)
          ycomp(i,1)=ycomp(i,1)+damp*qcen*ycen
          ycomp(i,2)=ycomp(i,2)+damp*(1-qcen)*yoff
          ycomp(i,3)=ycomp(i,3)+fprob*fac*msub/(pi*powlg_r(i)**2)
       enddo
    enddo
    close(1)

    do i=1,npow
       y(i)=ycomp(i,1)+ycomp(i,2)+ycomp(i,3)
    enddo

  end subroutine masspro_mdist

  subroutine crosscorr_single(amp,rs,soff,y,qcen,ycomp)

    implicit none

    integer :: i
    double precision, intent(in) :: amp,rs,soff,qcen  ! m in unit of 10^14M_sun/h
    double precision, intent(out) :: y(1:npow),ycomp(1:npow,1:3)
    double precision :: c,mnfwc,xoff,x,ycen,yoff
    integer, parameter :: nmax=256
    double precision :: logrmin,logrmax,dlogr
    double precision, dimension(nmax) :: xc,yc,yc2

! halo/subhalo core size
    c=rvir_fid/rs

! integral factor
    mnfwc=log(1+c)-c/(1+c)

! smoothing scale normalized by rs
    xoff=soff/rs

    logrmin=-4.d0
    logrmax=4.d0
    dlogr=(logrmax-logrmin)/nmax
    do i=1,nmax
       xc(i)=10**(logrmin+(i-0.5)*dlogr)
       yc(i)=nfw2(xc(i))
    enddo

    call fftlog_smooth(xc,yc,nmax,0.d0,-4.d0,4.d0,xoff)
    call spline(xc,yc,nmax,1d30,1d30,yc2)
    
    do i=1,npow
       x=powlg_r(i)/rs
       call splint(xc,yc,yc2,nmax,x,yoff)
       ycen=nfw2(x)

       ycen=3*ycen
       yoff=3*yoff

       ycomp(i,1)=amp*qcen*ycen
       ycomp(i,2)=amp*(1-qcen)*yoff
       ycomp(i,3)=0.d0
       y(i)=ycomp(i,1)+ycomp(i,2)+ycomp(i,3)
    enddo

  end subroutine crosscorr_single

  subroutine crosscorr_mdist(ampin,camp,soff,y,qcen,ycomp,ispec)

    implicit none

    integer :: i,im
    integer, intent(in) :: ispec
    integer, parameter :: immax=30
    double precision, intent(in) :: ampin,camp,soff,qcen  ! m in unit of 10^14M_sun/h
    double precision, intent(out) :: y(1:npow),ycomp(1:npow,1:3)
    double precision :: m,c,rs,mnfwc,xoff,x,ycen,yoff,rhalo,fprob,amp,damp
    double precision :: coff,mnfwc_off,rs_off
    integer, parameter :: nmax=256
    double precision, parameter :: delc=1.676
    double precision :: logrmin,logrmax,dlogr
    double precision, dimension(nmax) :: xc,yc,yc2,xcin,ycin,ycoff
    double precision :: rave,fsum

    logrmin=-4.d0
    logrmax=4.d0
    dlogr=(logrmax-logrmin)/nmax
    do i=1,nmax
       xcin(i)=10**(logrmin+(i-0.5)*dlogr)
       ycin(i)=nfw2(xcin(i))
    enddo

    if (clrg(ispec)=='highpcen') then
       open(1,file='data/mdist_highpcen_cen.dat',status='old')
    else
       open(1,file='data/mdist_RMCGneBCG_cen.dat',status='old')
    endif

    ycomp=0.d0
    amp=0.d0

    rave=0.d0
    fsum=0.d0

    do im=1,immax

       read(1,*) m,c,fprob
       if (fprob<1e-5) cycle
       c=camp*c

! halo size
       rhalo=(3*m/(4*pi*delv*rhom))**(1./3.)

! halo core size
       rs=rhalo/c

! integral factor
       mnfwc=log(1+c)-c/(1+c)

! smoothing scale normalized by rs

       xc=xcin
       yc=ycin

       if (offmodel=='nfw') then
          coff=soff*c
          if (coff<1.) coff=1.
!
          mnfwc_off=log(1+coff)-coff/(1+coff)
          rs_off=rhalo/coff
          xoff=rs_off/rs
!
          do i=1,nmax
             ycoff(i)=nfw2(xc(i)/xoff)*(coff**2/(8*pi*pi*mnfwc))
!             ycoff(i)=nfw2_truncate(xc(i)/xoff,coff)*(coff/(4*pi*pi*mnfwc))
          enddo
          call fftlog_smooth_nfw(xc,yc,ycoff,nmax,0.d0,-4.d0,4.d0)
       else
          xoff=(soff*rhalo)/rs
          call fftlog_smooth(xc,yc,nmax,0.d0,-4.d0,4.d0,xoff)
       endif
       call spline(xc,yc,nmax,1d30,1d30,yc2)

       rave=rave+soff*rhalo*fprob
       fsum=fsum+fprob

       damp=fprob*m/mnfwc/(6*pi*rs*rs)/rhom
       amp=amp+damp

       do i=1,npow
          x=powlg_r(i)/rs
          call splint(xc,yc,yc2,nmax,x,yoff)
          ycen=nfw2(x)

          ycen=3*ycen
          yoff=3*yoff

          ycomp(i,1)=ycomp(i,1)+damp*qcen*ycen
          ycomp(i,2)=ycomp(i,2)+damp*(1-qcen)*yoff
          ycomp(i,3)=ycomp(i,3)+0.d0
       enddo
    enddo
    close(1)

    rave=rave/fsum
!    write(*,*) rave

!    open(12,file='check_corr_nfw.dat',status='unknown')
!    do i=1,npow
!       write(12,*) powlg_r(i),ycomp(i,1)/qcen,ycomp(i,2)/(1-qcen)
!    enddo
!    close(12)
!    stop
!    write(*,*) amp,amp*ampin

    ycomp=ycomp/dchi*ampin

    do i=1,npow
       y(i)=ycomp(i,1)+ycomp(i,2)+ycomp(i,3)
    enddo

  end subroutine crosscorr_mdist

  subroutine crosscorr_mave(ampin,camp,soff,y,qcen,ycomp,ispec)

    implicit none

    integer :: i,im
    integer, intent(in) :: ispec
    integer, parameter :: immax=30
    double precision, intent(in) :: ampin,camp,soff,qcen  ! m in unit of 10^14M_sun/h
    double precision, intent(out) :: y(1:npow),ycomp(1:npow,1:3)
    double precision :: m,c,rs,mnfwc,xoff,x,ycen,yoff,rhalo,fprob,amp,damp,mt,ct
    integer, parameter :: nmax=256
    double precision, parameter :: delc=1.676
    double precision :: logrmin,logrmax,dlogr
    double precision, dimension(nmax) :: xc,yc,yc2,xcin,ycin

    logrmin=-4.d0
    logrmax=4.d0
    dlogr=(logrmax-logrmin)/nmax
    do i=1,nmax
       xcin(i)=10**(logrmin+(i-0.5)*dlogr)
       ycin(i)=nfw2(xcin(i))
    enddo

    if (clrg(ispec)=='highpcen') then
       open(1,file='data/mdist_highpcen_cen.dat',status='old')
    else
       open(1,file='data/mdist_RMCGneBCG_cen.dat',status='old')
    endif

    m=0.d0
    c=0.d0
    do im=1,immax
       read(1,*) mt,ct,fprob
       m=m+mt*fprob
       c=c+ct*fprob
    enddo
    close(1)

    ycomp=0.d0
    amp=0.d0

! halo size
    rhalo=(3*m/(4*pi*delv*rhom))**(1./3.)

! halo/subhalo core size
!       c=camp*gzmean**0.54*5.9*(delc/sigm(m))**(-0.35)
!    c=camp*7.85*(m/2e12)**-0.081*(1+z)**-0.71
       
! halo core size
    rs=rhalo/c

! integral factor
    mnfwc=log(1+c)-c/(1+c)

! smoothing scale normalized by rs
    xoff=(soff*rhalo)/rs

    xc=xcin
    yc=ycin

    call fftlog_smooth(xc,yc,nmax,0.d0,-4.d0,4.d0,xoff)
    call spline(xc,yc,nmax,1d30,1d30,yc2)
    
    damp=m/mnfwc/(6*pi*rs*rs)/rhom
    amp=amp+damp

    do i=1,npow
       x=powlg_r(i)/rs
       call splint(xc,yc,yc2,nmax,x,yoff)
       ycen=nfw2(x)

       ycen=3*ycen
       yoff=3*yoff

       ycomp(i,1)=ycomp(i,1)+damp*qcen*ycen
       ycomp(i,2)=ycomp(i,2)+damp*(1-qcen)*yoff
       ycomp(i,3)=ycomp(i,3)+0.d0
    enddo

    ycomp=ycomp/dchi*ampin

    do i=1,npow
       y(i)=ycomp(i,1)+ycomp(i,2)+ycomp(i,3)
    enddo

  end subroutine crosscorr_mave

  subroutine masspro_genNFW(m,c,soff,fsub,y,qcen,ycomp,alpha,beta)

    implicit none

    integer :: i,j
    double precision, intent(in) :: m,soff,fsub,c,qcen,alpha,beta  ! m in unit of 10^14M_sun/h
    double precision, intent(out) :: y(1:npow),ycomp(1:npow,1:3)
    double precision :: rhalo,rs,mnfwc,ycen,yoff,xoff,fac,x
    integer, parameter :: nmax=256
    double precision :: logrmin,logrmax,dlogr
    double precision, dimension(nmax) :: xc,yc,yc2,ys,ys2

! factor
    fac=1e-12 ! for changinge the unit of scale from Mpc^2 to pc^-2

! halo size
    rhalo=(3*m*munit/(4*pi*delv*rhom))**(1./3.)
    
! halo core size
    rs=rhalo/c

! integral factor for generalized NFW profile
    call mnfwc_genNFW(alpha,beta,c,mnfwc)

! smoothing scale normalized by rs
    xoff=soff/rs

    logrmin=-4.d0
    logrmax=4.d0
    dlogr=(logrmax-logrmin)/nmax
    do i=1,nmax
       xc(i)=10**(logrmin+(i-0.5)*dlogr)
! projected generalized NFW profile
       yc(i)=projected_genNFW(alpha,beta,xc(i))
       ys(i)=yc(i)
    enddo

    call fftlog_smooth(xc,yc,nmax,2.d0,-4.d0,4.d0,0.d0)
    call spline(xc,yc,nmax,1d30,1d30,yc2)
    call fftlog_smooth(xc,ys,nmax,2.d0,-4.d0,4.d0,xoff)
    call spline(xc,ys,nmax,1d30,1d30,ys2)
    
    do i=1,npow
       x=powlg_r(i)/rs
       call splint(xc,yc,yc2,nmax,x,ycen)
       call splint(xc,ys,ys2,nmax,x,yoff)
       ycomp(i,1)=fac*m*munit*qcen*ycen/mnfwc/(2*pi*rs*rs)
       ycomp(i,2)=fac*m*munit*(1-qcen)*yoff/mnfwc/(2*pi*rs*rs)
       ycomp(i,3)=fac*m*fsub*munit/pi/powlg_r(i)**2
       y(i)=ycomp(i,1)+ycomp(i,2)+ycomp(i,3)
    enddo

  end subroutine masspro_genNFW

  subroutine crosscorr_genNFW(amp,rs,soff,y,qcen,ycomp,alpha,beta)

    implicit none

    integer :: i
    double precision, intent(in) :: amp,rs,soff,qcen,alpha,beta
    double precision, intent(out) :: y(1:npow),ycomp(1:npow,1:3)
    double precision :: c,mnfwc,xoff,x,ycen,yoff
    integer, parameter :: nmax=256
    double precision :: logrmin,logrmax,dlogr
    double precision, dimension(nmax) :: xc,yc,yc2

! halo/subhalo core size
    c=rvir_fid/rs

! integral factor for generalized NFW
    call mnfwc_genNFW(alpha,beta,c,mnfwc)

! smoothing scale normalized by rs
    xoff=soff/rs

    logrmin=-4.d0
    logrmax=4.d0
    dlogr=(logrmax-logrmin)/nmax
    do i=1,nmax
       xc(i)=10**(logrmin+(i-0.5)*dlogr)
! projected generalized NFW profile
       yc(i)=projected_genNFW(alpha,beta,xc(i))
    enddo

    call fftlog_smooth(xc,yc,nmax,0.d0,-4.d0,4.d0,xoff)
    call spline(xc,yc,nmax,1d30,1d30,yc2)
    
    do i=1,npow
       x=powlg_r(i)/rs
       call splint(xc,yc,yc2,nmax,x,yoff)
! projected generalized NFW profile
       ycen=projected_genNFW(alpha,beta,x)
       ycomp(i,1)=amp*qcen*ycen/(2*pi*rs*rs)
       ycomp(i,2)=amp*(1-qcen)*yoff/(2*pi*rs*rs)
       ycomp(i,3)=0.d0
       y(i)=ycomp(i,1)+ycomp(i,2)+ycomp(i,3)
    enddo

  end subroutine crosscorr_genNFW

  subroutine mnfwc_genNFW(a,b,c,y)

    implicit none
    integer :: n
    integer, parameter :: nmax=100
    double precision, intent(in) :: a,b,c
    double precision :: x,y,dlnx
    double precision, parameter :: xmin=1d-2

    y=0.d0
    dlnx=log(c/xmin)/nmax
    do n=1,nmax
       x=xmin*exp(dlnx*(n-0.5))
       y=y+x**(3-b)/(1+x)**(a-b)*dlnx
    enddo
    
  end subroutine mnfwc_genNFW

  double precision function projected_genNFW(a,b,t)

    implicit none
    integer :: n
    integer, parameter :: nmax=100
    double precision, intent(in) :: a,b,t
    double precision :: r,z,dlnz,sigma,zmin=1.d-6,zmax=1.d6

    dlnz=log(zmax/zmin)/nmax
    sigma=0.d0
    do n=1,nmax
       z=zmin*exp(dlnz*(n-0.5))
       r=sqrt(t*t+z*z)
       sigma=sigma+z*dlnz/(r**b*(1+r)**(a-b))
    enddo
    projected_genNFW=sigma

  end function projected_genNFW

  double precision function nfw2int(x)

    implicit none
    double precision, intent(in) :: x

    if (abs(x-1)<1.e-8) then
       nfw2int=2*(1-log(2.))
    elseif (x<1) then
       nfw2int=2/x/x*(2/sqrt(1-x*x)*atanh(sqrt((1-x)/(1+x)))+log(x/2))
    else
       nfw2int=2/x/x*(2/sqrt(x*x-1)*atan(sqrt((x-1)/(x+1)))+log(x/2))
    endif
    nfw2int=nfw2int-nfw2(x)
    
  end function nfw2int

  double precision function nfw2int_truncate(x,c)

    implicit none
    double precision, intent(in) :: x,c
    double precision :: mnfwc

    if (abs(x-1)<1.e-8) then
       nfw2int_truncate=1/(3*(1+c))*((11*c+10)*sqrt(c*c-1)/(1+c)-6*c)+2*log((1+c)/(c+sqrt(c*c-1)))
    elseif (x<1) then
       nfw2int_truncate=1/(x*x*(1+c))*((2-x*x)*sqrt(c*c-x*x)/(1-x*x)-2*c)+2/x**2*log(x*(1+c)/(c+sqrt(c*c-x*x))) &
            +(2-3*x*x)/(x*x*(1-x*x)**1.5)*acosh((x*x+c)/(x*(1+c)))
    elseif (x>c) then
       mnfwc=log(1+c)-c/(1+c)
       nfw2int_truncate=2*mnfwc/x**2
    else
       nfw2int_truncate=1/(x*x*(1+c))*((2-x*x)*sqrt(c*c-x*x)/(1-x*x)-2*c)+2/x**2*log(x*(1+c)/(c+sqrt(c*c-x*x))) &
            -(2-3*x*x)/(x*x*(x*x-1)**1.5)*acos((x*x+c)/(x*(1+c)))
    endif
    
  end function nfw2int_truncate


  double precision function nfw2_truncate(x,c)

    implicit none
    double precision, intent(in) :: x,c
    
    if (abs(x-1)<1.e-8) then
       nfw2_truncate=sqrt(c*c-1)/(3*(1+c))*(1+1/(c+1))
    elseif (x<1) then
       nfw2_truncate=-sqrt(c*c-x*x)/(1-x*x)/(1+c)+(1-x*x)**(-1.5)*acosh((x*x+c)/(x*(1+c)))
    elseif (x<c) then
       nfw2_truncate=-sqrt(c*c-x*x)/(1-x*x)/(1+c)-(x*x-1)**(-1.5)*acos((x*x+c)/(x*(1+c)))
    else
       nfw2_truncate=0.
    endif
    
  end function nfw2_truncate

  double precision function nfw2(x)

    implicit none
    double precision, intent(in) :: x
    
    if (abs(x-1)<1.e-8) then
       nfw2=1./3.
    elseif (x<1) then
       nfw2=1/(x*x-1)*(1-2/sqrt(1-x*x)*atanh(sqrt((1-x)/(1+x))))
    else
       nfw2=1/(x*x-1)*(1-2/sqrt(x*x-1)*atan(sqrt((x-1)/(x+1))))
    endif
    
  end function nfw2


  double precision function acosh(x)

    double precision, intent(in) :: x
    acosh=log(x+sqrt(x*x-1))

  end function acosh

SUBROUTINE spline(x,y,n,yp1,ypn,y2)
  IMPLICIT NONE
  INTEGER, INTENT(in) :: n
  DOUBLE PRECISION, INTENT(in) :: yp1,ypn
  DOUBLE PRECISION, DIMENSION(n), INTENT(in) :: x,y
  DOUBLE PRECISION, DIMENSION(n), INTENT(out) :: y2
  INTEGER, parameter :: NMAX=1500
  INTEGER :: i,k
  DOUBLE PRECISION :: p,qn,sig,un
  DOUBLE PRECISION, DIMENSION(NMAX) :: u
  if (yp1>.99e30) then
     y2(1)=0.
     u(1)=0.
  else
     y2(1)=-0.5
     u(1)=(3./(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
  endif
  do i=2,n-1
     sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
     p=sig*y2(i-1)+2.
     y2(i)=(sig-1.)/p
     u(i)=(6.*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1))/(x(i)-x(i-1))) &
          /(x(i+1)-x(i-1))-sig*u(i-1))/p
  enddo
  if (ypn>.99e30) then
     qn=0.
     un=0.
  else
     qn=0.5
     un=(3./(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
  endif
  y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.)
  do k=n-1,1,-1
     y2(k)=y2(k)*y2(k+1)+u(k)
  enddo
  
END SUBROUTINE spline

! *********************************************************************
SUBROUTINE splint(xa,ya,y2a,n,x,y)
! *********************************************************************
  IMPLICIT NONE
  INTEGER, INTENT(in) :: n
  DOUBLE PRECISION, INTENT(in) :: x
  DOUBLE PRECISION, INTENT(in), DIMENSION(n) :: xa, y2a, ya
  DOUBLE PRECISION, INTENT(out) :: y
  INTEGER :: k,khi,klo
  DOUBLE PRECISION :: a,b,h

  klo=1
  khi=n
  do
     k=(khi+klo)/2
     if(xa(k)>x)then
        khi=k
     else
        klo=k
     endif
     if (khi-klo<=1) exit
  end do
  h=xa(khi)-xa(klo)
  if (h==0.) then
     write(*,*) 'bad xa input in splint'
     stop
  endif
  a=(xa(khi)-x)/h
  b=(x-xa(klo))/h
  y=a*ya(klo)+b*ya(khi)+((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**2)/6.

END SUBROUTINE splint

subroutine calcinvarray(covorg,incov,nmax)

!     read the data of cov(nmax,nmax)
!     and output their inversed array incov(nmax,nmax)
!      format of array is  
!          1, 1, array(1,1)
!          2, 1, array(2,1)
!          3, 1, array(3,1)
!               ...
!          1, 2, array(1,2)
!          2, 2, array(2,2)
!               ...
!          nmax-1,nmax,array(nmax-1,nmax)
!          nmax,nmax,array(nmax,nmax)
!
!      data type of array is real*8
!
  implicit none

  integer, intent(in) ::  nmax
  double precision, intent(in) :: covorg(nmax,nmax)
  double precision, intent(out) :: incov(nmax,nmax)
  double precision :: cov(nmax,nmax)
  integer ::  k, nx, ny, indx(nmax)
  double precision :: d, col(nmax)

  cov=covorg
  call ludcmp(cov,nmax,nmax,indx,d)
  do ny=1,nmax
     col=0.d0
     col(ny)=1.d0
     call lubksb(cov,nmax,nmax,indx,col)
     do nx=1,nmax
        incov(nx,ny)=col(nx)
     enddo
  enddo

end subroutine calcinvarray

! ------------------------------------------------

SUBROUTINE ludcmp(a,n,np,indx,d)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: n, np
  INTEGER, INTENT(OUT) :: indx(n)
  INTEGER, PARAMETER :: NMAX=1300
  DOUBLE PRECISION, INTENT(OUT) :: d
  DOUBLE PRECISION, INTENT(INOUT) :: a(np,np)
  DOUBLE PRECISION, PARAMETER :: TINY=1.0d-20
  INTEGER :: i,imax,j,k
  DOUBLE PRECISION :: aamax,dum,sum,vv(NMAX)

  d=1.
  do i=1,n
     aamax=0.d0
     do j=1,n
        if (abs(a(i,j))>aamax) aamax=abs(a(i,j))
     enddo
     if (aamax==0.) then
        write(*,*) 'singular matrix in ludcmp'
        stop
     endif
     vv(i)=1./aamax
  enddo
  do j=1,n
     do i=1,j-1
        sum=a(i,j)
        do k=1,i-1
           sum=sum-a(i,k)*a(k,j)
        enddo
        a(i,j)=sum
     enddo
     aamax=0.
     do i=j,n
        sum=a(i,j)
        do k=1,j-1
           sum=sum-a(i,k)*a(k,j)
        enddo
        a(i,j)=sum
        dum=vv(i)*abs(sum)
        if (dum.ge.aamax) then
           imax=i
           aamax=dum
        endif
     enddo
     if (j/=imax)then
        do k=1,n
           dum=a(imax,k)
           a(imax,k)=a(j,k)
           a(j,k)=dum
        enddo
        d=-d
        vv(imax)=vv(j)
     endif
     indx(j)=imax
     if (a(j,j)==0.) a(j,j)=TINY
     if (j/=n) then
        dum=1./a(j,j)
        do i=j+1,n
           a(i,j)=a(i,j)*dum
        enddo
     endif
  enddo

END SUBROUTINE ludcmp

! -----------------------------------------

SUBROUTINE lubksb(a,n,np,indx,b)
  INTEGER, INTENT(IN) :: n,np
  INTEGER, INTENT(IN) :: indx(n)
  DOUBLE PRECISION, INTENT(IN) :: a(np,np)
  DOUBLE PRECISION, INTENT(INOUT) :: b(n)
  INTEGER :: i,ii,j,ll
  DOUBLE PRECISION :: sum

  ii=0
  do i=1,n
     ll=indx(i)
     sum=b(ll)
     b(ll)=b(i)
     if (ii/=0) then
        do j=ii,i-1
           sum=sum-a(i,j)*b(j)
        enddo
     else if (sum/=0.) then
        ii=i
     endif
     b(i)=sum
  enddo
  do i=n,1,-1
     sum=b(i)
     do j=i+1,n
        sum=sum-a(i,j)*b(j)
     enddo
     b(i)=sum/a(i,i)
  enddo

END SUBROUTINE lubksb


subroutine fftlog_smooth(x,y,n,mu,logrmin,logrmax,rad)

  implicit none

!!!  n must be smaller than NMAX
  integer :: n
  integer, parameter :: NMAX=4096
  integer :: lnblnk
  character :: inf*100,outf*100
  integer :: dir,i,kropt
  logical :: ok
  double precision :: x(n),y(n),rad,yorg
  double precision :: dlnr,dlogr,k,kr,logkc,logrc,logrmax,logrmin,nc,q,r,rk,mu
  double precision :: wsave(2*NMAX+3*(NMAX/2)+19)
  real :: criden,angdis
  character :: type*1,comp*1

  q=0.d0
  kr=1.d0
  kropt=2

! forward transform

  dir=1
  
  logrc=(logrmin+logrmax)/2.d0

  nc=dble(n+1)/2.d0
  dlogr=(logrmax-logrmin)/n
  dlnr=dlogr*log(10.d0)
  
  do i=1,n
     y(i)=y(i)*x(i)
  enddo

  call fhti(n,0.d0,q,dlnr,kr,kropt,wsave,ok)
  if (.not.ok) then
     write(*,*) 'strange'
     stop
  endif
  logkc=log10(kr)-logrc
  call fht(n,y,dir,wsave)

! add Gaussian window functions
  do i=1,n
     x(i)=10.d0**(logkc+(i-nc)*dlogr)
     y(i)=y(i)*exp(-0.5*(rad*x(i))**2)
  enddo

!     backward transform

  q=0.d0
  kr=1.d0
  kropt=2

  dir=-1

  call fhti(n,mu,q,dlnr,kr,kropt,wsave,ok)
  if (.not.ok) then
     write(*,*) 'strange'
     stop
  endif
  logrc=log10(kr)-logkc
  call fht(n,y,dir,wsave)

  do i=1,n
     x(i)=10.d0**(logrc+(i-nc)*dlogr)
     y(i)=y(i)/x(i)
  enddo

end subroutine fftlog_smooth

subroutine fftlog_smooth_nfw(x,y,yoff,n,mu,logrmin,logrmax)
  
  implicit none

!!!  n must be smaller than NMAX
  integer :: n
  integer, parameter :: NMAX=4096
  integer :: lnblnk
  character :: inf*128,outf*128
  integer :: dir,i,kropt
  logical :: ok
  double precision :: x(n),y(n),yoff(n),filter
  double precision :: dlnr,dlogr,k,kr,logkc,logrc,logrmax,logrmin,nc,q,r,rk,mu,qoff
  double precision :: wsave(2*NMAX+3*(NMAX/2)+19),wsave_off(2*NMAX+3*(NMAX/2)+19)
  real :: criden,angdis
  character :: type*1,comp*1

  q=0.d0
  kr=1.d0
  kropt=2

!  forward transform

  dir=1
  
  logrc=(logrmin+logrmax)/2.d0
  
  nc=dble(n+1)/2.d0
  dlogr=(logrmax-logrmin)/n
  dlnr=dlogr*log(10.d0)
  
  do i=1,n
     y(i)=y(i)*x(i)
     yoff(i)=yoff(i)*x(i)
  enddo

  call fhti(n,0.d0,q,dlnr,kr,kropt,wsave,ok)
  call fhti(n,0.d0,qoff,dlnr,kr,kropt,wsave_off,ok)
  if (.not.ok) then
     write(*,*) 'strange'
     stop
  endif
  logkc=log10(kr)-logrc
  call fht(n,y,dir,wsave)
  call fht(n,yoff,dir,wsave_off)

  do i=1,n
     x(i)=10.d0**(logkc+(i-nc)*dlogr)
     if (i<n/4) then
        filter=1
     else
        filter=(yoff(i)/x(i))/(yoff(n/4)/x(n/4))
     endif
     y(i)=y(i)*filter
  enddo

!     backward transform

  q=0.d0
  kr=1.d0
  kropt=2

  dir=-1

  call fhti(n,mu,q,dlnr,kr,kropt,wsave,ok)
  if (.not.ok) then
     write(*,*) 'strange'
     stop
  endif
  logrc=log10(kr)-logkc
  call fht(n,y,dir,wsave)
  
  do i=1,n
     x(i)=10.d0**(logrc+(i-nc)*dlogr)
     y(i)=y(i)/x(i)
  enddo

end subroutine fftlog_smooth_nfw

end MODULE Readmodel
