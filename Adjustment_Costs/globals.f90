!include "toolbox.f90"

module globals

    use toolbox

    implicit none

    !~~~~~~~~~~~~~~~~~!
    !                 !
    !   AGGREGATES    !
    !                 !
    !~~~~~~~~~~~~~~~~~!

    real*8               :: pstar, Ra  ! endogenous prices
    real*8               :: Agg_DI, Agg_Gov, Agg_Loan, &
                            Agg_Sec, Agg_Dep, Agg_Whole, Agg_Corp, Agg_Def2, Agg_Liq
    real*8, dimension(2) :: Agg_Def, Agg_Div, Agg_Def_stationary

    real*8, parameter    :: beta_pre = 1d0/1.0025d0  ! household discount factor
    real*8, parameter    :: Rd       = 1d0/beta_pre   ! equilibrium deposit rate

    !~~~~~~~~~~~~~~~~~~~!
    !                   !
    !   BANK PROBLEM    !
    !                   !
    !~~~~~~~~~~~~~~~~~~~!

    ! secondary security investors
    real*8, parameter                                  :: Abar = 1d0

    ! corporate tax
    real*8, parameter                                  :: tax = .33d0

    ! loan cost functions
    real*8, parameter                                  :: theta      = .024d0
    real*8, parameter                                  :: nu         = 2d0                                                   ! curvature

    real*8, dimension(20)                               :: Rl   ! bank return process
    real*8, dimension(size(Rl))                        :: prob_l

    real*8, parameter                                  :: l_mu = 1.01d0   ! median loan returns : 1-x = net median return (.04), then log

    ! old deposit grid
    real*8, parameter            :: bank_da = 0d0
    real*8                       :: bank_db
    integer, parameter           :: bank_dlen = 30
    real*8, dimension(bank_dlen) :: bank_dgrid
    real*8, parameter            :: dgrowth = .1d0

    ! networth grid
    real*8                                             :: nstar
    real*8, parameter                                  :: bank_na = 0d0
    real*8                                             :: bank_nb
    integer, parameter                                 :: bank_nlen = 30
    real*8, dimension(bank_nlen)                       :: bank_ngrid
    real*8, parameter                                  :: growth    =  .1d0 !0.001d0 !.01d0

    real*8, parameter                                  :: hair      = -.6995d0  !  29% of wholesale funding is secured with 105% collateral ==> 1.05*.29 = .3045 thus hair = -.6955

    ! funding shock process
    real*8, dimension(9), parameter                    :: delta     =  (/ 0d0 , .1d0, .25d0, .37d0, .51d0, .635d0, .765d0, &
                                                                            .92d0, 1d0  /)
    real*8, dimension(9), parameter                    :: prob_d    = (/ .5d0,  .4d0, .06d0, .02d0, .01d0, .004d0, .002d0, &
                                                                                .0035d0, .0005d0 /)

                                                        !(/ .9975d0,  0d0, 0d0, 0d0, 0d0, 0d0, 0d0, &
                                                        !                0d0, .0025d0 /)





    ! risk weights
    real*8, parameter                                  :: omega_l  = 1d0
    real*8, parameter                                  :: omega_s  = 1d0

    ! liquidity ratio violating fraction
    real*8, parameter                                  :: liq_b = 0d0

    !Capital and liquidity requirements
    real*8                                             :: ebar
    real*8                                             :: phi_lr
    real*8, parameter :: lr_hair = .083d0 !.15d0

    ! no liquidation price
    real*8, parameter :: no_liq_price = .999d0

    ! bank liquidation value
    real*8, parameter                                  :: liq    = 0.65d0  !.85d0

    !Individual policy functions
    real*8, dimension(bank_nlen,bank_dlen)             :: v, vnew, lpol, spol, dpol, div_pol, apol, epol
    real*8, dimension(bank_nlen,bank_dlen,size(delta)) :: stilde
    real*8, dimension(bank_nlen,bank_dlen)             :: F_stationary
    real*8, dimension(bank_nlen,bank_dlen)             :: default_prob, default_liq_prob

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
    !                            !
    !   Structural Parameters    !
    !                            !
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~!

    real*8                         :: kappa   ! deposit adjustment cost
    real*8, parameter              :: phi = 0d0            ! equity issuance cost
    real*8                         :: l_sigma2       ! loan return volatility
    real*8                         :: gam            ! bank discount factor
    real*8                         :: alpha          ! outside security investor elasticity
    real*8                         :: i_s             ! return on securities


    real*8                         :: beta           ! total bank discount factor
    real*8                         :: omega_outside  ! outside investor endowment

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
    !                                       !
    !   Value Function Solver Parameters    !
    !                                       !
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!

    integer :: grid_len
    real*8, allocatable, dimension(:) :: lgrid, divgrid, cgrid, sgrid, dgrid, agrid

    real*8, parameter :: la   = .000001d0
    real*8, parameter :: aa = la
    real*8, parameter :: sa = la
    real*8, parameter :: da = la

    ! upper bound of search grids
    real*8, parameter :: lb = .01d0
    real*8, parameter :: ab = lb*.5d0
    real*8, parameter :: sb = lb*.5d0
    real*8, parameter :: db = lb*.5d0

    real*8, parameter :: diva = 0d0
    real*8, parameter :: divb = .01d0

    ! policy function search length
    real*8                         :: delta_l

    real*8                         :: delta_div
    real*8                         :: delta_a
    real*8                         :: delta_d
    real*8                         :: delta_s

    ! discounting for looking at lower values
    real*8, parameter                         :: back_look = .95d0

    ! penalty for breaking constraint
    real*8, parameter                         :: penalty = -10000000d0

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
    !                                !
    !   General Calibration Stuff    !
    !                                !
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!

    ! loan-to-security, deposit-to-wholesale, leverage ratio, insolvency default rate, liquidity ratio, deposit share 1, deposit share 2, deposit share 3
    real*8, dimension(6) :: data_moms = (/ 3.35d0, 3.20d0, 9.61d0, 25.5d0, 73.3d0, 45.3d0 /)

    !~~~~~~~~~~~~~~~~!
    !                !
    !   MPI Stuff    !
    !                !
    !~~~~~~~~~~~~~~~~!
    integer :: nprocs,my_id,ierr


contains

    ! return index values of sorted array
    function index_sort(x,xsize)

        implicit none

        real*8, intent(in), dimension(xsize) :: x
        integer, intent(in) :: xsize

        integer, dimension(xsize) :: index_sort
        integer :: ii
        integer, dimension(1) :: idx
        real*8, dimension(xsize) :: x_temp

        x_temp = x

        do ii=1,xsize

            idx = maxloc( x_temp )

            index_sort(ii) = idx(1)

            x_temp(idx(1)) = -1000000000d0

        enddo

    end function

    !Origination cost
    function g(x)

        implicit none

        real*8, intent(in) :: x
        real*8 :: g

        g = theta*x**nu/nu

    end function

    function phid(oldd,newd)

        implicit none

        real*8, intent(in) :: oldd, newd
        real*8             :: phid

        phid = -kappa*( newd - oldd)**2d0

    end function

    ! dividend adjustment cost
    function phi_div(x)

        implicit none

        real*8, intent(in) :: x
        real*8 :: phi_div

        if (x <= 0) then
            phi_div = abs(phi*x) !**2d0
        else
            phi_div = 0d0
        endif

    end function

    function Aggregate_Corp_Tax()

        implicit none

        real*8  :: Aggregate_Corp_Tax
        real*8  :: temp_tax
        integer :: ii,jj,kk,ll,mm,nn
        real*8 :: excess_cash, ctilde, stilde, net_int, stock, net_temp

        temp_tax = 0d0

        ! for each bank state
        do ii=1,bank_nlen
            do jj=1,bank_dlen

                ! each funding shock
                do ll=1,size(delta)

                    ! if no liquidity default
                    if ( pstar*spol(ii,jj) > delta(ll)*apol(ii,jj) ) then

                        ! compute liquidations
                        stilde = (delta(ll)*apol(ii,jj)-ctilde)/pstar

                        ! for each loan return
                        do mm=1,size(Rl)

                            net_int = (Rl(mm)-1d0)*lpol(ii,jj) + i_s*(spol(ii,jj)-stilde) &
                                             - (Rd-1d0)*dpol(ii,jj) - (Ra-1d0)*(1d0-delta(ll))*apol(ii,jj)

                            stock = lpol(ii,jj) + (spol(ii,jj)-stilde) -&
                                    dpol(ii,jj) - (1d0-delta(ll))*apol(ii,jj)

                            if ( net_int > 0d0 ) then
                                net_temp = (1d0-tax)*net_int + stock
                            else
                                net_temp = net_int + stock
                            endif

                            ! if no default
                            if ( net_temp > nstar ) then

                                ! record tax income
                                temp_tax = temp_tax + F_stationary(ii,jj)*prob_d(ll)*prob_l(mm)*&
                                                                     tax*net_int

                            endif
                        enddo
                    endif
                enddo
            enddo
        enddo

        Aggregate_Corp_Tax = temp_tax

    end function

    function Aggregate_Gov()

        implicit none

        real*8 :: Aggregate_Gov
        real*8 :: temp_gov
        integer :: ii, jj, kk

        do ii=1,bank_nlen
            do jj=1,bank_dlen

                temp_gov = temp_gov + F_stationary(ii,jj)*spol(ii,jj)

            enddo
        enddo

        Aggregate_Gov = i_s*temp_gov

    end function

    function bank_shares()

        implicit none

        real*8, dimension(4) :: bank_shares  ! loan, securities, deposits, wholesale

        real*8 :: total_assets, total_liabs, total_loan, total_sec, total_debt, total_whole, total_net
        integer :: ii, jj, kk

        ! portfolio shares
        bank_shares = 0d0

        total_assets = 0d0
        total_liabs  = 0d0
        total_loan   = 0d0
        total_sec    = 0d0
        total_debt   = 0d0
        total_whole  = 0d0
        total_net    = 0d0

        do jj=1,bank_dlen
            do kk=1,bank_nlen

                total_assets = total_assets + F_stationary(kk,jj)*( lpol(kk,jj) + &
                spol(kk,jj) )
                total_liabs  = total_liabs  + F_stationary(kk,jj)*( dpol(kk,jj) + &
                apol(kk,jj) + bank_ngrid(kk) )

                total_loan = total_loan + F_stationary(kk,jj)*lpol(kk,jj)
                total_sec  = total_sec  + F_stationary(kk,jj)*spol(kk,jj)

                total_debt  = total_debt  + F_stationary(kk,jj)*dpol(kk,jj)
                total_whole = total_whole + F_stationary(kk,jj)*apol(kk,jj)
                total_net   = total_net   + F_stationary(kk,jj)*bank_ngrid(kk)
            enddo
        enddo

        bank_shares(1) = total_loan/total_assets
        bank_shares(2) = total_sec/total_assets
        bank_shares(3) = total_debt/total_liabs
        bank_shares(4) = total_whole/total_liabs

    end function

    function Aggregate_BS()

        implicit none

        real*8, dimension(4) :: Aggregate_BS  !(loans,security,deposits,wholesale)

        integer :: ii, jj, kk
        real*8 :: temp_loan, temp_sec, temp_dep, temp_whole

        temp_loan = 0d0
        temp_sec  = 0d0

        do ii=1,bank_nlen
            do kk=1,bank_dlen

                temp_loan = temp_loan + F_stationary(ii,kk)*lpol(ii,kk)
                temp_sec  = temp_sec  + F_stationary(ii,kk)*spol(ii,kk)

                temp_dep    = temp_dep    + F_stationary(ii,kk)*dpol(ii,kk)
                temp_whole  = temp_whole  + F_stationary(ii,kk)*apol(ii,kk)

            enddo
        enddo

        Aggregate_BS(1) = temp_loan
        Aggregate_BS(2) = temp_sec
        Aggregate_BS(3) = temp_dep
        Aggregate_BS(4) = temp_whole

    end function

    function Aggregate_Div()

        implicit none

        real*8, dimension(2) :: Aggregate_Div

        ! banking sector
        integer :: ii, jj, kk, ll, mm, nn
        real*8 :: temp_div

        ! money market sector
        real*8 :: agg_a, obj1, obj2, obj3, obj4, temp_net, rtilde_a
        real*8, dimension(bank_nlen,bank_dlen) :: bshare
        real*8 :: excess_cash, net_int, stock

        ! compute aggregate dividends from banking sector
        temp_div = 0d0

        do ii = 1,bank_nlen
            do kk=1,bank_dlen

                temp_div = temp_div + F_stationary(ii,kk)*div_pol(ii,kk)

            enddo
        enddo

        Aggregate_Div(1) = temp_div

        ! compute aggregate dividends from money market sector
        agg_a = 0d0
        do ii=1,bank_nlen
            do kk=1,bank_dlen
                agg_a = agg_a + F_stationary(ii,kk)*apol(ii,kk)
            enddo
        enddo

        do ii=1,bank_nlen
            do kk=1,bank_dlen
                bshare(ii,kk) = F_stationary(ii,kk)*apol(ii,kk)/agg_a
            enddo
        enddo

        obj1 = 0d0
        obj2 = 0d0
        obj3 = 0d0
        obj4 = 0d0

        do ii=1,bank_nlen
            do kk=1,bank_dlen

                ! liquidity shocks
                do ll=1,size(delta)

                    if ( delta(ll)*apol(ii,kk) > pstar*spol(ii,kk) ) then  ! liquidity default cash flow
                        obj1 = obj1 + bshare(ii,kk)*prob_d(ll)
                    else

                        obj2 = obj2 + bshare(ii,kk)*prob_d(ll)*delta(ll)                       ! non-default liquidity shock cash flow

                        do mm=1,size(Rl)

                            net_int = (Rl(mm)-1d0)*lpol(ii,kk) + i_s*(spol(ii,kk)-stilde(ii,kk,ll)) &
                                             - (Rd-1d0)*dpol(ii,kk) - (Ra-1d0)*(1d0-delta(ll))*apol(ii,kk)

                            stock = lpol(ii,kk) + (spol(ii,kk)-stilde(ii,kk,ll)) -&
                                    dpol(ii,kk) - (1d0-delta(ll))*apol(ii,kk)

                            if ( net_int > 0d0 ) then
                                temp_net = (1d0-tax)*net_int + stock
                            else
                                temp_net = net_int + stock
                            endif

                            if ( temp_net < nstar ) then
                                obj3 = obj3 + bshare(ii,kk)*prob_d(ll)*(1d0-delta(ll))*prob_l(mm) ! insolvency default cash flow
                            else
                                obj4 = obj4 + bshare(ii,kk)*prob_d(ll)*(1d0-delta(ll))*prob_l(mm)  ! repayment cash flow
                            endif
                        enddo
                    endif
                enddo
            enddo
        enddo

        ! gross return on wholesale lending
        rtilde_a = (1d0+hair)*obj1 + obj2 + (1d0+hair)*obj3 + Ra*obj4

        Aggregate_Div(2) = (rtilde_a-1d0)*agg_a

    end function

    function Aggregate_Dep()

        implicit none

        real*8 :: Aggregate_Dep

        integer :: ii, jj, kk
        real*8 :: temp_dep

        temp_dep = 0d0

        do ii = 1,bank_nlen
            do kk=1,bank_dlen

                temp_dep = temp_dep + F_stationary(ii,kk)*dpol(ii,kk)

            enddo
        enddo

        Aggregate_Dep = temp_dep

    end function

    function Aggregate_DI()

        implicit none

        real*8 :: Aggregate_DI

        integer :: ii, jj, kk, ll, mm, nn
        integer, dimension(1) :: min_idx
        real*8 :: r_thresh
        real*8 :: def_diff, temp_DI, excess_cash, net_temp, liq_assets, stilde, ctilde
        real*8 :: net_int, stock
        real*8 :: net_position

        ! initialize
        temp_DI = 0d0

        ! for each (n_b,dbar)
        do ii = 1,bank_nlen
            do kk=1,bank_dlen

                ! for each funding shock
                do ll=1,size(delta)

                    ! if liquidity default
                    if ( pstar*spol(ii,kk) < delta(ll)*apol(ii,kk) ) then

                        !
                        !    In default a fraction (1-liq) of assets are lost.  Wholesale lenders immediately seize collateral (1+h)a
                        !
                        !    Notice that liquidity default is "worse" in some sense then insolvency default, as the assets ( l + s ) were not allowed to mature
                        !

                        liq_assets = liq*( lpol(ii,kk) + pstar*spol(ii,kk)  -&
                                                                    (1d0+hair)*apol(ii,kk) )

                        !
                        !   after colalteral seizures, bank still has debt repayments (1) deposits R^d*d and
                        !         (2) remaining wholesale debt   max( Ra*a - (1+h)a, 0 ) ; i.e. it repays more if collateral seizures
                        !           where not enough to cover full repayment
                        !
                        !
                        !   net_position = - ( assets - deposit obligations )   so if debt >> liquidated assets, net_position > 0
                        !
                        !   all that matters for the deposit insurance formula is whether or not deposits have been paid in full
                        !

                        net_position = -( liq_assets - Rd*dpol(ii,kk) )  !- maxval((/apol(ii,jj,kk)*(Ra -(1d0+hair)) , 0d0/))

                        ! if negative position, no cost of deposit insurance, so set to 0
                        if ( net_position < 0d0 ) then
                            net_position = 0d0
                        endif

                        temp_DI = temp_DI + F_stationary(ii,kk)*prob_d(ll)*net_position

                    else

                        ! compute liquidations
                        stilde = (delta(ll)*apol(ii,kk))/pstar

                        ! for each loan return
                        do mm=1,size(Rl)

                            net_int = (Rl(mm)-1d0)*lpol(ii,kk) + i_s*(spol(ii,kk)-stilde) &
                                             - (Rd-1d0)*dpol(ii,kk) - (Ra-1d0)*(1d0-delta(ll))*apol(ii,kk)

                            stock = lpol(ii,kk) + (spol(ii,kk)-stilde) -&
                                    dpol(ii,kk) - (1d0-delta(ll))*apol(ii,kk)

                            if ( net_int > 0d0 ) then
                                net_temp = (1d0-tax)*net_int + stock
                            else
                                net_temp = net_int + stock
                            endif

                            ! if insolvency default
                            if ( net_temp <= nstar ) then

                                ! in insolvency default, assets (l+s) have interest maturity
                                liq_assets = liq*( Rl(mm)*lpol(ii,kk) + (1d0+i_s)*(spol(ii,kk)-stilde) -&
                                                (1d0+hair)*apol(ii,kk)*(1d0-delta(ll)) )

                                net_position = -( liq_assets - Rd*dpol(ii,kk) ) !-(1d0-delta(ll))*maxval( (/ apol(ii,jj,kk)*(Ra -(1d0+hair)) , 0d0 /) ) )

                                if ( net_position < 0d0 ) then
                                    net_position = 0d0
                                endif

                                ! compute cost of deposit insurance, net of assets
                                temp_DI = temp_DI + F_stationary(ii,kk)*prob_d(ll)*prob_l(mm)*net_position

                            endif

                        enddo

                    endif

                enddo
            enddo
        enddo

        Aggregate_DI = temp_DI

    end function

    function Aggregate_Liq()

        implicit none

        real*8 :: Aggregate_Liq

        integer :: ii, jj, kk, ll, mm, nn
        integer, dimension(1) :: min_idx
        real*8 :: r_thresh
        real*8 :: def_diff, excess_cash, net_temp, liq_assets, stilde, ctilde
        real*8 :: temp_liq

        ! initialize
        temp_liq = 0d0

        ! for each (n_b,theta,dbar)
        do ii = 1,bank_nlen
            do jj =1,bank_dlen

                ! for each funding shock
                do ll=1,size(delta)

                    ! if no liquidity default
                    if ( pstar*spol(ii,jj) >= delta(ll)*apol(ii,jj) ) then

                        ! compute liquidations
                        stilde = delta(ll)*apol(ii,jj)/pstar

                        temp_liq = temp_liq + stilde*prob_d(ll)*F_stationary(ii,jj)

                    else   ! liquidity default

                        temp_liq = temp_liq + spol(ii,jj)*prob_d(ll)*F_stationary(ii,jj)

                    endif

                enddo

            enddo
        enddo

        Aggregate_Liq = temp_liq

    end function

    function weighted_def()

        implicit none

        real*8, dimension(2) :: weighted_def

        integer :: ii, jj, kk, ll, mm, nn
        integer, dimension(1) :: min_idx
        real*8 :: r_thresh
        real*8 :: def_diff, agg_a, excess_cash, net_temp, liq_assets, stilde, wholesale_weight

        ! compute aggregate wholesale funding
        agg_a = 0d0

        do ii=1,bank_nlen
            do jj=1,bank_dlen

                agg_a = agg_a + F_stationary(ii,jj)*apol(ii,jj)

            enddo
        enddo

        ! compute weighted default rates
        do ii=1,bank_nlen
            do jj=1,bank_dlen

                ! compute wholesale weight
                wholesale_weight = F_stationary(ii,jj)*apol(ii,jj)/agg_a

                ! liquidity default
                weighted_def(1) = wholesale_weight*default_liq_prob(ii,jj)

                ! insolvency default
                weighted_def(2) = wholesale_weight*(1d0-default_liq_prob(ii,jj))*default_prob(ii,jj)

            enddo
        enddo

    end function

    function Aggregate_Def()

        implicit none

        ! input/output variables
        real*8, dimension(2) :: Aggregate_Def

        ! local variables
        integer :: ii,jj,kk
        real*8, dimension(2)  :: temp_def

        temp_def = 0d0

        do ii=1,bank_nlen
            do jj=1,bank_dlen

                temp_def(1) = temp_def(1) + F_stationary(ii,jj)*default_prob(ii,jj)
                temp_def(2) = temp_def(2) + F_stationary(ii,jj)*default_liq_prob(ii,jj)

            enddo
        enddo

        Aggregate_Def = temp_def

    end function


    function model_moments()

        implicit none

        integer :: ii, jj, kk
        real*8 :: total_liabs, total_dep, cap, total_eq, total_rwa, liq_num, liq_den
        real*8 :: agg_loan, agg_sec, agg_dep, agg_whole

        real*8, dimension(size(data_moms)) :: model_moments

        !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
        !    loan-to-security ratio    !
        !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
        agg_loan = 0d0
        agg_sec  = 0d0

        do ii=1,bank_nlen
            do jj=1,bank_dlen

                agg_loan = agg_loan + F_stationary(ii,jj)*lpol(ii,jj)

                agg_sec  = agg_sec  + F_stationary(ii,jj)*spol(ii,jj)

            enddo
        enddo

        model_moments(1) = agg_loan/agg_sec

        !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
        !    deposit-to-wholesale ratio    !
        !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
        agg_dep = 0d0
        agg_whole  = 0d0

        do ii=1,bank_nlen
            do jj=1,bank_dlen

                agg_dep    = agg_dep + F_stationary(ii,jj)*dpol(ii,jj)

                agg_whole  = agg_whole  + F_stationary(ii,jj)*apol(ii,jj)

            enddo
        enddo

        model_moments(2) = agg_dep/agg_whole


        !~~~~~~~~~~~~~~~~~~~~~~~~~~~!
        !    risk-weighted ratio    !
        !~~~~~~~~~~~~~~~~~~~~~~~~~~~!
        total_eq  = 0d0
        total_rwa = 0d0

        do ii=1,bank_nlen
            do jj=1,bank_dlen

                total_eq  = total_eq  + F_stationary(ii,jj)*( lpol(ii,jj) +&
                            spol(ii,jj) - (dpol(ii,jj)+apol(ii,jj)))

                total_rwa = total_rwa + F_stationary(ii,jj)*lpol(ii,jj)

            enddo
        enddo

        model_moments(3) = 100d0*total_eq/total_rwa

        !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
        !    insolvency default rate ratio    !
        !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
        Agg_Def = Aggregate_Def()

        model_moments(4) = 100d0*100d0*Agg_Def(1)

        !~~~~~~~~~~~~~~~~~~~~~~~!
        !    liquidity ratio    !
        !~~~~~~~~~~~~~~~~~~~~~~~!
        liq_num = 0d0
        liq_den = 0d0

        do ii=1,bank_nlen
            do jj=1,bank_dlen

                if (apol(ii,jj) >0d0) then

                    liq_num = liq_num + F_stationary(ii,jj)*( (1d0-lr_hair)*spol(ii,jj))

                    liq_den = liq_den + F_stationary(ii,jj)*( apol(ii,jj) )
                endif

            enddo
        enddo

        model_moments(5) = 100d0*liq_num/liq_den

        !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
        !    compute deposits shares for bank types 1,2,3    !
        !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!

        total_liabs  = 0d0
        total_dep   = 0d0

        do jj=1,bank_dlen
            do kk=1,bank_nlen

                total_liabs  = total_liabs  + F_stationary(kk,jj)*( dpol(kk,jj) + &
                apol(kk,jj) + bank_ngrid(kk) )

                total_dep  = total_dep  + F_stationary(kk,jj)*dpol(kk,jj)
            enddo
        enddo

        model_moments(6) = 100d0*total_dep/total_liabs


    end function

    function fire_elastic(s)

        implicit none

        real*8, intent(in) :: s
        real*8 :: fire_elastic

        fire_elastic = (alpha-1d0)*s/(s+omega_outside)

    end function



end module
