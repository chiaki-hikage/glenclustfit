module CalcLike
 use Random
 use settings
 use ParamDef
 use Readmodel

 implicit none

contains

  function GenericLikelihoodFunction(Params) 

    implicit none
    integer :: i,j,n,j2
    type(ParamSet)  Params 
    real :: GenericLikelihoodFunction
    double precision :: m,c,amp,rs,soff(nlrg),fsub(nlrg),qcen(nlrg),amp2h,beta,alpha
    double precision :: powlg_theo(npow,nlrg),powlg_comp(npow,3,nlrg)

! constraints on the sum of off-centering fraction of RMCG and BCG less than unity
    if (Params%P(10)+Params%P(11)>1) then
       GenericLikelihoodFunction=1000
       goto 99
    endif

    do n=1,nlrg
       if (measurement=='dsigma') then
          m=Params%P(1)
          c=Params%P(2)
       else if (measurement=='angcor') then
          amp=Params%P(1)
          c=Params%P(2)
       endif
       fsub(n)=Params%P(n+5)
       soff(n)=Params%P(n+2)
       qcen(n)=Params%P(n+8)
       if (genNFW=='n') then
          if (measurement=='dsigma') then
!             call masspro_singlemass(m,c,soff(n),fsub(n),powlg_theo(1:npow,n),qcen(n),powlg_comp(1:npow,1:3,n))
             call masspro_mdist(c,soff(n),fsub(n),powlg_theo(1:npow,n),qcen(n),powlg_comp(1:npow,1:3,n),n)
          else if (measurement=='angcor') then
!             call crosscorr_single(amp,rs,soff(n),powlg_theo(1:npow,n),qcen(n),powlg_comp(1:npow,1:3,n))
             call crosscorr_mdist(amp,c,soff(n),powlg_theo(1:npow,n),qcen(n),powlg_comp(1:npow,1:3,n),n)
          endif
       else
          alpha=Params%P(13)
          beta=Params%P(14)
          if (measurement=='dsigma') then
             call masspro_genNFW(m,c,soff(n),fsub(n),powlg_theo(1:npow,n),qcen(n),powlg_comp(1:npow,1:3,n),alpha,beta)
          else if (measurement=='angcor') then
             call crosscorr_genNFW(amp,rs,soff(n),powlg_theo(1:npow,n),qcen(n),powlg_comp(1:npow,1:3,n),alpha,beta)
          endif
       endif
    enddo

    amp2h=Params%P(12)

    GenericLikelihoodFunction=0.            

    do i=1,npow
       do j=1,nlrg
          powlg_theo(i,j) = powlg_theo(i,j) + amp2h*powlg_twohalo(i)
       enddo
    enddo

    do i=1,npow
       do j=1,nlrg
          do j2=1,nlrg
             GenericLikelihoodFunction = GenericLikelihoodFunction + &
                  (powlg_obs(i,j)-powlg_theo(i,j))*(powlg_obs(i,j2)-powlg_theo(i,j2))*powlg_cov_inv(i,j,j2)
          enddo
       enddo
    enddo

99 end function GenericLikelihoodFunction


  function GetLogLike(Params) !Get -Ln(Likelihood)
    type(ParamSet)  Params 
    real GetLogLike
    real dum(1,1)
 
    if (any(Params%P > Scales%PMax) .or. any(Params%P < Scales%PMin)) then
       GetLogLike = logZero
        return
    end if

    GetLogLike = GenericLikelihoodFunction(Params) 

   end function GetLogLike
   

end module CalcLike
