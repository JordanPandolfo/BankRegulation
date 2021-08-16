PROGRAM compute_eq

    use globals

    include 'mpif.h'

    integer :: tt, iterator

    call MPI_init(ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD, my_id, ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierr)

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
    !                                    !
    !   Set External Parameter Values    !
    !                                    !
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
    ebar   = .04d0
    phi_lr = 0.001d0
    nstar = 0d0

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
    !                                                 !
    !   Set Internally Calibrated Parameter Values    !
    !                                                 !
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
    ! new calibration solutions with i_s = 1.0035

    if (my_id==0)then
        call tic()
    endif

    kappa      = 0d0
    gam        = 0.96d0
    pstar      = 0.97d0
    l_sigma2   = 0.000416d0
    i_s        = .0039d0

    beta = gam*beta_pre

    alpha = 0.0002126d0
    omega_outside = (no_liq_price/alpha)**(1d0/(alpha-1d0))

    ! relevant state space and VFI grid search parameters
    !
    !   - determine ub of networth state space, and ub of deposit state space's relation to it
    !
    !   - determine width of grid search for policy functions
    !
    !   - determine # of state space grid points + # of grid search grid points
    !
    !   - determine max VFI iterations and VFI convergence criterion
    !
    bank_nb   = .01d0       ! upper bound of state space
    bank_db   = bank_nb*(1d0-ebar)/ebar  ! upper bound of state space

    delta_l   = 20d0*( bank_nb )/bank_nlen   ! local search window for policy function grid search

    delta_div = delta_l*.5d0
    delta_a   = delta_l*1d0
    delta_d   = delta_l*1d0
    delta_s   = delta_l*1d0

    grid_len = 5

    ! MPI BARRIER
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
    !                                   !
    !   Solve Stationary Equilibrium    !
    !                                   !
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!

    allocate( lgrid(grid_len) )
    allocate( divgrid(grid_len) )
    allocate( sgrid(grid_len) )
    allocate( dgrid(grid_len) )
    allocate( agrid(grid_len) )


    if (my_id==0) then
        write(*,*)
        write(*,*)
        write(*,*)
        write(*,*) '          Solving Stationary Equilibrium          '
        write(*,*)
        write(*,*)
        write(*,*)
    endif

    call initialize()

    call solve_ge()
    Agg_Def_stationary = Aggregate_Def()

    ! MPI BARRIER
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)

    ! output
    call output()

    ! MPI BARRIER
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)


    if (my_id==0) then
        call toc()
    endif

    ! MPI BARRIER
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)

contains

    !Initialize grids and constraint sets
    subroutine initialize()

        implicit none

        integer :: ii,jj,kk,iterator
        real*8, dimension(bank_nlen) :: temp

        ! inidividual bank networth grid
        call grid_Cons_Grow(bank_ngrid,bank_na,bank_nb,growth)
        call grid_Cons_Grow(bank_dgrid,bank_da,bank_db,dgrowth)

        ! initialize value function
        do ii = 1,bank_nlen
            do kk=1,bank_dlen
                v(ii,kk) =  bank_ngrid(ii)
            enddo
        enddo

        do ii=1,bank_nlen
            temp(ii) = ii
        enddo

        call normal_discrete( Rl, prob_l, l_mu, l_sigma2)

        if (my_id == 0) then
            write(*,*)
            write(*,*) 'Mean return:', sum(Rl*prob_l)
            write(*,*) 'Return vol:', ( sum(prob_l*( Rl - sum(Rl*prob_l) )**2d0) )**(0.5d0)
            write(*,*) 'Rl:', Rl
            write(*,*) 'Prob Rl:', prob_l
            write(*,*)
        endif

    end subroutine


    subroutine solve_ge()

        implicit none

        real*8  ::  Ra_implied, pstar_implied
        real*8 :: ra_a, ra_b, p_a, p_b
        integer :: jj

        ra_a = 1.0033996884d0
        ra_b = 1.00339968841d0

        p_a = 0.9710773393d0
        p_b = 0.97107733931d0

        ! solve GE
        do jj=1,1

            ! midpoint
            Ra    =    (Ra_a + Ra_b)/2d0
            pstar =    (p_a + p_b)/2d0

            ! solve model
            call vfi_bank()
            call balance_pols()
            call get_dist()

            Ra_implied = implied_ra_rate()
            Agg_Liq = Aggregate_Liq()
            pstar_implied = alpha*(Agg_Liq+omega_outside)**(alpha-1d0)

            if (my_id==0) then
                write(*,*)
                write(*,*) 'Agg Liq:', Agg_Liq
                write(*,*) 'omega:', omega_outside
                write(*,*)
                write(*,*) 'Security liquidations:',Agg_Liq
                write(*,*) 'pstar guess:',pstar
                write(*,*) 'pstar implied:', pstar_implied
                write(*,*)
                write(*,*) 'Ra guess:',Ra
                write(*,*) 'Ra implied:', Ra_implied
                write(*,*)
                write(*,*) 'Avg Error (basis points):', 100d0*100d0*(abs(Ra-Ra_implied)/Ra + abs(pstar-pstar_implied)/pstar)/2d0
                write(*,*)
            endif


            if ( Ra_implied > Ra) then
                Ra_a = Ra
            else
                Ra_b = Ra
            endif

            if ( pstar_implied > pstar) then
                p_a = pstar
            else
                p_b = pstar
            endif

            if ( 100d0*100d0*(abs(Ra-Ra_implied)/Ra + abs(pstar-pstar_implied)/pstar)/2d0 <1d0 ) then
                exit
            endif

        enddo


    end subroutine


    subroutine balance_pols()

        implicit none

        ! compute balancing period policy functions
        integer :: ii,jj,kk,ll
        real*8, dimension(bank_nlen,bank_dlen,size(delta)) :: ctilde, stilde
        real*8 :: excess_cash

        do ii=1,bank_nlen
            do kk=1,bank_dlen

                do ll=1,size(delta)

                    if (delta(ll)*apol(ii,kk) > pstar*spol(ii,kk)) then ! liquidity default

                        stilde(ii,kk,ll) = spol(ii,kk)

                    else

                        stilde(ii,kk,ll) = delta(ll)*apol(ii,kk)/pstar

                    endif

                enddo
            enddo
        enddo

    end subroutine

    ! VFI iterator for stationary equilibrium, calls bank_solve_stationary()
    subroutine vfi_bank()

        implicit none

        ! input/output variables

        ! VFI variables/parameters
        real*8, parameter :: tol = 1e-30
        real*8            :: error
        integer           :: iterator, tt
        integer, parameter :: it_max = 100

        ! Solving bank problem
        integer :: ii, jj, kk, ll, mm, nn, oo, pp, qq,rr,ss,zz
        real*8 :: v_temp, RHS, net_temp, e_temp
        real*8, dimension(5) :: v_curr_max, l_curr_max, d_curr_max, div_curr_max, s_curr_max, a_curr_max
        real*8 :: a_bound, l_bound, s_bound, d_bound, liq_temp, implied_a, &
                  implied_s, d_ubound, implied_d, implied_l, implied_div
        real*8 :: s_tilde  ! liquidation values
        integer, dimension(1) :: v_idx
        real*8 :: net_int, stock, pi_bound
        real*8 :: temp_pen

        ! linear interpolation variables
        integer :: ixl, ixr, iyl, iyr
        real*8  :: phix, phiy

        ! spline interpolation variables
        real*8, dimension( bank_dlen+2, bank_nlen+2 ) :: coeff

        !---------------!
        !   mpi stuff   !
        !---------------!
        real*8, allocatable, dimension(:) :: ltemp, dtemp, stemp, atemp, divtemp, vtemp, ctemp
        real*8, allocatable, dimension(:) :: ltemp_all, dtemp_all, stemp_all, atemp_all, divtemp_all, vtemp_all, ctemp_all
        integer, dimension(2) :: grid_idx   ! 2 states (nb,do)
        integer, allocatable, dimension(:,:,:) :: idx_full_grid

        integer :: total_grid, proc_len

        ! total size of state space grid
        total_grid = bank_nlen*bank_dlen

        if ( mod(total_grid,nprocs) == 0) then  ! perfect fit
            proc_len = total_grid/nprocs
        else
            proc_len = total_grid/nprocs + 1     ! the last proccessor will have some extra/empty slots
        endif

        ! construct temporary policy functions
        allocate( ltemp(proc_len) )
        allocate( dtemp(proc_len) )
        allocate( stemp(proc_len) )
        allocate( atemp(proc_len) )
        allocate( divtemp(proc_len) )
        allocate( vtemp(proc_len) )

        allocate( ltemp_all(proc_len*nprocs) )
        allocate( dtemp_all(proc_len*nprocs) )
        allocate( stemp_all(proc_len*nprocs) )
        allocate( atemp_all(proc_len*nprocs) )
        allocate( divtemp_all(proc_len*nprocs) )
        allocate( vtemp_all(proc_len*nprocs) )

        allocate( idx_full_grid(nprocs,proc_len,2) )

        call grid_f(nprocs,proc_len,idx_full_grid)

        ! MPI BARRIER
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)

        ! outter loop for VFI
        do iterator = 1,it_max
            ! spline interpolator over grid
            !do jj=1,size(theta)
            !    do kk=1,dbar_size
            !        call spline_interp(v(:,jj,kk), coeff(jj,kk,:))
            !    enddo
            !enddo

            !---------------------------!
            !   solve value function    !
            !---------------------------!
            do ii=1,proc_len  ! for each point in parallel grid; state space = (bank ngrid, do)

                ! determine indices for bank state space
                grid_idx = idx_full_grid(my_id+1,ii,:)

                ! initialize current problem
                v_curr_max   = 0d0
                l_curr_max   = 0d0
                d_curr_max   = 0d0
                div_curr_max = 0d0
                s_curr_max   = 0d0
                a_curr_max   = 0d0

                ! construct grids for control variables
                if ( ( grid_idx(1) == 1 ).or.( grid_idx(2) == 1 ).or.( iterator == 1 ) ) then  ! if n_b(1) or first iteratior of value function, use default grids

                    call grid_Cons_Equi( lgrid, la,lb)
                    call grid_Cons_Equi( sgrid, sa, sb )
                    call grid_Cons_Equi( dgrid, da,db)
                    call grid_Cons_Equi( agrid, aa,ab)
                    call grid_Cons_Equi( divgrid,diva, divb )

                else

                    d_bound = maxval( (/ 0.0000001d0, dpol(grid_idx(1)-1,grid_idx(2)-1) - back_look*delta_d /) )
                    a_bound = maxval( (/ 0.0000001d0, apol(grid_idx(1)-1,grid_idx(2)-1) - back_look*delta_a /) )
                    l_bound = maxval( (/ 0.0000001d0, lpol(grid_idx(1)-1,grid_idx(2)-1) - back_look*delta_l /) )
                    pi_bound = maxval( (/ 0d0, div_pol(grid_idx(1)-1,grid_idx(2)-1) - delta_div /) )
                    s_bound = maxval( (/ 0.0000001d0, spol(grid_idx(1)-1,grid_idx(2)) - back_look*delta_s /) )

                    call grid_Cons_Equi( lgrid, l_bound, lpol(grid_idx(1)-1,grid_idx(2)-1) + delta_l )
                    call grid_Cons_Equi( dgrid, d_bound, dpol(grid_idx(1)-1,grid_idx(2)-1) + delta_d )
                    call grid_Cons_Equi( agrid, a_bound, apol(grid_idx(1)-1,grid_idx(2)-1) + delta_a )
                    call grid_Cons_Equi( divgrid, pi_bound , div_pol(grid_idx(1)-1,grid_idx(2)-1) + delta_div )
                    call grid_Cons_Equi( sgrid, s_bound, spol(grid_idx(1)-1,grid_idx(2)) + delta_s )

                endif

                !
                !
                !    Case 1: binding capital requirement
                !
                !    implication:   I) d = ( l + s )*(1-ebar) - a
                !
                !                  II) div = n + a + d + Phi(do,d) - l- theta*l^2/2 - s
                !
                !    search over: (l,s,a) --> d --> div
                !
                !

                ! solve problem
                do ll=1,grid_len ! for each loan
                    do ss=1,grid_len   ! for each security
                        do zz=1,grid_len   ! for each wholesale

                            ! implied deposits
                            implied_d = (lgrid(ll) + sgrid(ss) )*(1d0-ebar) - agrid(zz)

                            ! implied dividends
                            implied_div =  bank_ngrid(grid_idx(1)) + agrid(zz)  + implied_d +&
                                           phid(bank_dgrid(grid_idx(2)),implied_d) -lgrid(ll) -&
                                           g(lgrid(ll)) - sgrid(ss)

                            ! liquidity requirement
                            liq_temp = ( (1d0-lr_hair)*sgrid(ss) )/agrid(zz)

                            if (implied_d < 0d0) then
                                v_temp = penalty
                            elseif (implied_div < 0d0) then
                                v_temp = penalty
                            elseif (liq_temp < phi_lr) then
                                v_temp = penalty
                            elseif ( sgrid(ss) < (1d0+hair)*agrid(zz)) then  ! collateral constraint
                                v_temp = penalty
                            else               ! if all constraints satisfied

                                ! initialize RHS of bellman equation
                                RHS = 0d0

                                ! for each possible funding shock
                                do qq=1,size(delta)

                                    ! check liquidity default
                                    if ( pstar*sgrid(ss) >= delta(qq)*agrid(zz)) then

                                        s_tilde = ( delta(qq)*agrid(zz) )/pstar

                                        ! no funding risk
                                        do rr=1,size(Rl)        ! for each next period loan shock

                                            ! compute networth
                                            net_int = (Rl(rr)-1d0)*lgrid(ll) + i_s*(sgrid(ss)-s_tilde) -&
                                                      (Rd-1d0)*implied_d -(Ra-1d0)*(1d0-delta(qq))*agrid(zz)

                                            stock = lgrid(ll) + (sgrid(ss)-s_tilde) - implied_d - (1d0-delta(qq))*agrid(zz)

                                            if ( net_int > 0d0) then
                                                net_temp = (1d0-tax)*net_int + stock
                                            else
                                                net_temp = net_int + stock
                                            endif

                                            if (net_temp >0d0) then

                                                ! linear interpolate for value
                                                call linint_Grow( net_temp, bank_na, bank_nb,&
                                                                    growth,bank_nlen-1, ixl, ixr, phix)

                                                call linint_Grow( implied_d, bank_da, bank_db,&
                                                                    dgrowth,bank_dlen-1, iyl, iyr, phiy)

                                                if (phix <= phiy) then

                                                    v_temp = phix*v(ixl+1,iyl+1) +&
                                                             (phiy-phix)*v(ixr+1,iyl+1) +&
                                                             (1d0-phiy)*v(ixr+1,iyr+1)

                                                else

                                                    v_temp = phiy*v(ixl+1,iyl+1) +&
                                                             (phix-phiy)*v(ixl+1,iyr+1) +&
                                                             (1d0-phix)*v(ixr+1,iyr+1)

                                                endif

                                                RHS = RHS + prob_d(qq)*prob_l(rr)*v_temp
                                            endif

                                        enddo
                                    endif
                                enddo

                                v_temp = implied_div + beta*RHS

                            endif

                            ! evaluate value of policy function combination
                            if (v_temp > v_curr_max(1)) then! if beats current, create new argmax and max value

                                v_curr_max(1)   = v_temp

                                l_curr_max(1)   = lgrid(ll)
                                s_curr_max(1)   = sgrid(ss)
                                d_curr_max(1)   = implied_d
                                a_curr_max(1)   = agrid(zz)
                                div_curr_max(1) = implied_div

                            endif
                        enddo
                    enddo
                enddo  ! end of policy function loop

                !
                !
                !    Case 2: binding liquidity requirement
                !
                !    implication:   I) a = (1-h_s)*s/phi_lr
                !
                !                  II) div = n + a + d + Phi(do,d) - l- theta*l^2/2 - s
                !
                !    search over: (l,s,d) --> a --> div
                !
                !

                ! solve problem
                do ll=1,grid_len ! for each loan
                    do oo=1,grid_len   ! for each deposit
                        do nn=1,grid_len  ! for each security

                            ! implies wholesale funding
                            implied_a = (1d0-lr_hair)*sgrid(nn)/phi_lr

                            ! implied dividend
                            implied_div =  bank_ngrid(grid_idx(1)) + implied_a  + dgrid(oo) +&
                                           phid(bank_dgrid(grid_idx(2)),dgrid(oo)) -lgrid(ll) -&
                                           g(lgrid(ll)) - sgrid(nn)

                            ! implied capital ratio
                            e_temp = (lgrid(ll) + sgrid(nn) - (dgrid(oo)+implied_a))/&
                                    ( lgrid(ll) + sgrid(nn) )

                            if (e_temp > 1d0) then ! capital requirement (above)
                                v_temp = penalty
                            elseif (e_temp < ebar) then ! capital requirement (below)
                                v_temp = penalty
                            elseif (implied_a < 0d0) then   ! liquidity requirement
                                v_temp = penalty
                            elseif (implied_div < 0d0) then     ! negative wholesale funding
                                v_temp = penalty
                            elseif ( sgrid(nn) < (1d0+hair)*implied_a) then  ! collateral constraint
                                v_temp = penalty
                            else               ! if all constraints satisfied

                                ! initialize RHS of bellman equation
                                RHS = 0d0

                                ! for each possible funding shock
                                do qq=1,size(delta)

                                    ! check liquidity default
                                    if ( pstar*sgrid(nn) >= delta(qq)*implied_a) then

                                        s_tilde = ( delta(qq)*implied_a )/pstar

                                        do rr=1,size(Rl)        ! for each next period loan shock

                                            ! compute networth
                                            net_int = (Rl(rr)-1d0)*lgrid(ll) + i_s*(sgrid(nn)-s_tilde) - &
                                                      (Rd-1d0)*dgrid(oo) - &
                                                      (Ra-1d0)*(1d0-delta(qq))*implied_a

                                            stock = lgrid(ll) + (sgrid(nn)-s_tilde) - &
                                                    dgrid(oo) - (1d0-delta(qq))*implied_a

                                            if ( net_int > 0d0) then
                                                net_temp = (1d0-tax)*net_int + stock
                                            else
                                                net_temp = net_int + stock
                                            endif

                                            if (net_temp>0d0) then
                                                ! linear interpolate for value
                                                call linint_Grow( net_temp, bank_na, bank_nb,&
                                                                    growth,bank_nlen-1, ixl, ixr, phix)

                                                call linint_Grow( dgrid(oo), bank_da, bank_db,&
                                                                    dgrowth,bank_dlen-1, iyl, iyr, phiy)

                                                if (phix <= phiy) then

                                                    v_temp = phix*v(ixl+1,iyl+1) +&
                                                             (phiy-phix)*v(ixr+1,iyl+1) +&
                                                             (1d0-phiy)*v(ixr+1,iyr+1)

                                                else

                                                    v_temp = phiy*v(ixl+1,iyl+1) +&
                                                             (phix-phiy)*v(ixl+1,iyr+1) +&
                                                             (1d0-phix)*v(ixr+1,iyr+1)

                                                endif

                                                RHS = RHS + prob_d(qq)*prob_l(rr)*v_temp
                                            endif

                                        enddo
                                    endif

                                enddo

                                v_temp = implied_div + beta*RHS

                            endif

                            ! evaluate value of policy function combination

                            if (v_temp > v_curr_max(2)) then! if beats current, create new argmax and max value

                                v_curr_max(2) = v_temp

                                l_curr_max(2)   = lgrid(ll)
                                s_curr_max(2)   = sgrid(nn)
                                d_curr_max(2)   = dgrid(oo)
                                a_curr_max(2)   = implied_a
                                div_curr_max(2) = implied_div

                            endif

                        enddo
                    enddo
                enddo ! end of policy function loop

                !
                !
                !    Case 3: binding capital & liquidity requirement
                !
                !    implication:   I) s = phi_lr*a/(1-h_s)
                !
                !                  II) l = (a+d)/(1-e) - s
                !
                !                 III) div = n + a + d + Phi(do,d) - l- theta*l^2/2 - s
                !
                !
                !    search over: (a,d)   --> s --> l --> div
                !
                !

                ! solve problem
                do ll=1,grid_len ! for each wholesale
                    do oo=1,grid_len   ! for deposit

                        ! implied securities
                        implied_s =  phi_lr*agrid(ll)/(1d0-lr_hair)

                        ! implied loans
                        implied_l = (agrid(ll) + dgrid(oo))/(1d0-ebar) - implied_s

                        ! implied dividends
                        implied_div =  bank_ngrid(grid_idx(1)) + agrid(ll)  + dgrid(oo) +&
                                           phid(bank_dgrid(grid_idx(2)),dgrid(oo)) -implied_l -&
                                           g(implied_l) - implied_s

                        if (implied_s < 0d0) then   ! liquidity requirement
                            v_temp = penalty
                        elseif (implied_div < 0d0) then
                            v_temp = penalty
                        elseif (implied_l < 0d0) then     ! negative wholesale funding
                            v_temp = penalty
                        elseif ( implied_s < (1d0+hair)*agrid(ll)) then  ! collateral constraint
                            v_temp = penalty
                        else               ! if all constraints satisfied

                            ! initialize RHS of bellman equation
                            RHS = 0d0

                            ! for each possible funding shock
                            do qq=1,size(delta)

                                ! check liquidity default
                                if ( pstar*implied_s >= delta(qq)*agrid(ll)) then

                                    s_tilde = ( delta(qq)*agrid(ll) )/pstar

                                    do rr=1,size(Rl)        ! for each next period loan shock

                                        ! compute networth
                                        net_int = (Rl(rr)-1d0)*implied_l + i_s*(implied_s-s_tilde) - &
                                                  (Rd-1d0)*dgrid(oo) - &
                                                  (Ra-1d0)*(1d0-delta(qq))*agrid(ll)

                                        stock = implied_l + (implied_s-s_tilde) - &
                                                dgrid(oo) - (1d0-delta(qq))*agrid(ll)

                                        if ( net_int > 0d0) then
                                            net_temp = (1d0-tax)*net_int + stock
                                        else
                                            net_temp = net_int + stock
                                        endif

                                        if (net_temp>0d0) then

                                            ! linear interpolate for value
                                            call linint_Grow( net_temp, bank_na, bank_nb,&
                                                                growth,bank_nlen-1, ixl, ixr, phix)

                                            call linint_Grow( dgrid(oo), bank_da, bank_db,&
                                                                dgrowth,bank_dlen-1, iyl, iyr, phiy)

                                            if (phix <= phiy) then

                                                v_temp = phix*v(ixl+1,iyl+1) +&
                                                         (phiy-phix)*v(ixr+1,iyl+1) +&
                                                         (1d0-phiy)*v(ixr+1,iyr+1)

                                            else

                                                v_temp = phiy*v(ixl+1,iyl+1) +&
                                                         (phix-phiy)*v(ixl+1,iyr+1) +&
                                                         (1d0-phix)*v(ixr+1,iyr+1)

                                            endif

                                            RHS = RHS + prob_d(qq)*prob_l(rr)*v_temp
                                        endif

                                    enddo
                                endif

                            enddo

                            v_temp = implied_div + beta*RHS

                        endif

                        ! evaluate value of policy function combination

                        if (v_temp > v_curr_max(3)) then! if beats current, create new argmax and max value

                            v_curr_max(3) = v_temp

                            l_curr_max(3)   = implied_l
                            d_curr_max(3)   = dgrid(oo)
                            a_curr_max(3)   = agrid(ll)
                            div_curr_max(3) = implied_div

                        endif

                    enddo
                enddo  ! end of policy function loop

                !
                !    Case 4: no assumed binding constraints
                !

                ! solve problem
                do ll=1,grid_len ! for each loan
                    do mm=1,grid_len ! for each deposit
                        do oo=1,grid_len   ! for each dividend
                            do zz=1,grid_len  ! for each wholesale

                                ! implied securities
                                implied_s = bank_ngrid(grid_idx(1)) + agrid(zz) + dgrid(mm) +&
                                            phid(bank_dgrid(grid_idx(2)),dgrid(mm)) - divgrid(oo) -&
                                            lgrid(ll) - g(lgrid(ll))

                                ! capital ratio
                                e_temp = (lgrid(ll) + implied_s - dgrid(mm) - agrid(zz))/( lgrid(ll) + implied_s  )

                                ! liquidity ratio
                                liq_temp = (1d0-lr_hair)*implied_s/agrid(zz)

                                if (e_temp > 1d0) then ! capital requirement (above)
                                    v_temp = penalty
                                elseif (e_temp < ebar) then ! capital requirement (below)
                                    v_temp = penalty
                                elseif (liq_temp < phi_lr) then
                                    v_temp = penalty
                                elseif ( implied_s < (1d0+hair)*agrid(zz)) then  ! collateral constraint
                                    v_temp = penalty
                                elseif (implied_s < 0d0) then
                                    v_temp = penalty
                                else               ! if all constraints satisfied

                                    ! initialize RHS of bellman equation
                                    RHS = 0d0
                                    ! for each possible funding shock
                                    do qq=1,size(delta)

                                        ! check liquidity default
                                        if ( pstar*implied_s >= delta(qq)*agrid(zz)) then

                                            s_tilde = ( delta(qq)*agrid(zz) )/pstar

                                            ! no funding risk
                                            do rr=1,size(Rl)        ! for each next period loan shock

                                                ! compute networth
                                                net_int = (Rl(rr)-1d0)*lgrid(ll) + i_s*(implied_s-s_tilde) -&
                                                          (Rd-1d0)*dgrid(mm) - (Ra-1d0)*(1d0-delta(qq))*agrid(zz)

                                                stock = lgrid(ll) + (implied_s-s_tilde) - dgrid(mm) - (1d0-delta(qq))*agrid(zz)

                                                if ( net_int > 0d0) then
                                                    net_temp = (1d0-tax)*net_int + stock
                                                else
                                                    net_temp = net_int + stock
                                                endif

                                                if (net_temp>0d0) then

                                                    ! linear interpolate for value
                                                    call linint_Grow( net_temp, bank_na, bank_nb,&
                                                                        growth,bank_nlen-1, ixl, ixr, phix)

                                                    call linint_Grow( dgrid(mm), bank_da, bank_db,&
                                                                        dgrowth,bank_dlen-1, iyl, iyr, phiy)

                                                    if (phix <= phiy) then

                                                        v_temp = phix*v(ixl+1,iyl+1) +&
                                                                 (phiy-phix)*v(ixr+1,iyl+1) +&
                                                                 (1d0-phiy)*v(ixr+1,iyr+1)

                                                    else

                                                        v_temp = phiy*v(ixl+1,iyl+1) +&
                                                                 (phix-phiy)*v(ixl+1,iyr+1) +&
                                                                 (1d0-phix)*v(ixr+1,iyr+1)

                                                    endif

                                                    RHS = RHS + prob_d(qq)*prob_l(rr)*v_temp
                                                endif

                                            enddo
                                        endif
                                    enddo

                                    v_temp = divgrid(oo) + beta*RHS

                                endif

                                ! evaluate value of policy function combination

                                if (v_temp > v_curr_max(4)) then! if beats current, create new argmax and max value

                                    v_curr_max(4) = v_temp

                                    l_curr_max(4)   = lgrid(ll)
                                    s_curr_max(4)   = implied_s
                                    d_curr_max(4)   = dgrid(mm)
                                    a_curr_max(4)   = agrid(zz)
                                    div_curr_max(4) = divgrid(oo)

                                endif

                            enddo
                        enddo
                    enddo
                enddo  ! end of policy function loop

                !
                !    Case 5: d = old value
                !

                ! solve problem
                do ll=1,grid_len ! for each loan
                    do oo=1,grid_len   ! for each dividend
                        do zz=1,grid_len  ! for each wholesale

                            ! implied securities
                            implied_s = bank_ngrid(grid_idx(1)) + agrid(zz) + bank_dgrid(grid_idx(2)) +&
                                        phid(bank_dgrid(grid_idx(2)),bank_dgrid(grid_idx(2))) - divgrid(oo) -&
                                        lgrid(ll) - g(lgrid(ll))

                            ! capital ratio
                            e_temp = (lgrid(ll) + implied_s - bank_dgrid(grid_idx(2)) - agrid(zz))/( lgrid(ll) + implied_s  )

                            ! liquidity ratio
                            liq_temp = (1d0-lr_hair)*implied_s/agrid(zz)

                            if (e_temp > 1d0) then ! capital requirement (above)
                                v_temp = penalty
                            elseif (e_temp < ebar) then ! capital requirement (below)
                                v_temp = penalty
                            elseif (liq_temp < phi_lr) then
                                v_temp = penalty
                            elseif ( implied_s < (1d0+hair)*agrid(zz)) then  ! collateral constraint
                                v_temp = penalty
                            elseif (implied_s < 0d0) then
                                v_temp = penalty
                            else               ! if all constraints satisfied

                                ! initialize RHS of bellman equation
                                RHS = 0d0
                                ! for each possible funding shock
                                do qq=1,size(delta)

                                    ! check liquidity default
                                    if ( pstar*implied_s >= delta(qq)*agrid(zz)) then

                                        s_tilde = ( delta(qq)*agrid(zz) )/pstar

                                        ! no funding risk
                                        do rr=1,size(Rl)        ! for each next period loan shock

                                            ! compute networth
                                            net_int = (Rl(rr)-1d0)*lgrid(ll) + i_s*(implied_s-s_tilde) -&
                                                      (Rd-1d0)*bank_dgrid(grid_idx(2)) - (Ra-1d0)*(1d0-delta(qq))*agrid(zz)

                                            stock = lgrid(ll) + (implied_s-s_tilde) - bank_dgrid(grid_idx(2)) -&
                                                    (1d0-delta(qq))*agrid(zz)

                                            if ( net_int > 0d0) then
                                                net_temp = (1d0-tax)*net_int + stock
                                            else
                                                net_temp = net_int + stock
                                            endif

                                            if (net_temp>0d0) then

                                                ! linear interpolate for value
                                                call linint_Grow( net_temp, bank_na, bank_nb,&
                                                                    growth,bank_nlen-1, ixl, ixr, phix)

                                                call linint_Grow( bank_dgrid(grid_idx(2)), bank_da, bank_db,&
                                                                    dgrowth,bank_dlen-1, iyl, iyr, phiy)

                                                if (phix <= phiy) then

                                                    v_temp = phix*v(ixl+1,iyl+1) +&
                                                             (phiy-phix)*v(ixr+1,iyl+1) +&
                                                             (1d0-phiy)*v(ixr+1,iyr+1)

                                                else

                                                    v_temp = phiy*v(ixl+1,iyl+1) +&
                                                             (phix-phiy)*v(ixl+1,iyr+1) +&
                                                             (1d0-phix)*v(ixr+1,iyr+1)

                                                endif

                                                RHS = RHS + prob_d(qq)*prob_l(rr)*v_temp
                                            endif

                                        enddo
                                    endif
                                enddo

                                v_temp = divgrid(oo) + beta*RHS

                            endif

                            ! evaluate value of policy function combination

                            if (v_temp > v_curr_max(5)) then! if beats current, create new argmax and max value

                                v_curr_max(5) = v_temp

                                l_curr_max(5)   = lgrid(ll)
                                s_curr_max(5)   = implied_s
                                d_curr_max(5)   = bank_dgrid(grid_idx(2))
                                a_curr_max(5)   = agrid(zz)
                                div_curr_max(5) = divgrid(oo)

                            endif

                        enddo
                    enddo
                enddo  ! end of policy function loop

                !
                ! find max values of all the cases considered
                !
                v_idx = maxloc( v_curr_max )

                ! record argmax and max value
                vtemp(ii)    = v_curr_max( v_idx(1) )

                ltemp(ii)    = l_curr_max( v_idx(1) )
                stemp(ii)    = s_curr_max( v_idx(1) )

                dtemp(ii)    = d_curr_max( v_idx(1) )
                atemp(ii)    = a_curr_max( v_idx(1) )

                divtemp(ii) = div_curr_max( v_idx(1) )

            enddo   ! state space loop

            ! MPI GATHER
            call MPI_GATHER(vtemp,proc_len,MPI_REAL8,  vtemp_all,proc_len,   MPI_REAL8,0,MPI_COMM_WORLD,ierr)
            call MPI_GATHER(ltemp,proc_len,MPI_REAL8,  ltemp_all,proc_len,   MPI_REAL8,0,MPI_COMM_WORLD,ierr)
            call MPI_GATHER(dtemp,proc_len,MPI_REAL8,  dtemp_all,proc_len,   MPI_REAL8,0,MPI_COMM_WORLD,ierr)
            call MPI_GATHER(stemp,proc_len,MPI_REAL8,  stemp_all,proc_len,   MPI_REAL8,0,MPI_COMM_WORLD,ierr)
            call MPI_GATHER(atemp,proc_len,MPI_REAL8,  atemp_all,proc_len,   MPI_REAL8,0,MPI_COMM_WORLD,ierr)
            call MPI_GATHER(divtemp,proc_len,MPI_REAL8,divtemp_all,proc_len, MPI_REAL8,0,MPI_COMM_WORLD,ierr)

            ! MPI BARRIER
            call MPI_BARRIER(MPI_COMM_WORLD,ierr)

            if ( my_id == 0 ) then  ! if root

                ! construct dimension-based policy functions
                vnew = transform_f(vtemp_all)
                lpol = transform_f(ltemp_all)
                dpol = transform_f(dtemp_all)
                spol = transform_f(stemp_all)
                apol = transform_f(atemp_all)
                div_pol = transform_f(divtemp_all)

                ! compute error
                error = sum( (v-vnew)**2d0 )/(bank_nlen*bank_dlen)

                write(*,*) 'Root node, VFI iteration:', iterator,'error:',error

                ! update value function
                v = vnew

            endif

            ! broadcast value function error
            call MPI_BCAST( error,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr  )

            ! broadcast policy functions
            call MPI_BCAST(vtemp_all,  proc_len*nprocs,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
            call MPI_BCAST(ltemp_all,  proc_len*nprocs,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
            call MPI_BCAST(dtemp_all,  proc_len*nprocs,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
            call MPI_BCAST(stemp_all,  proc_len*nprocs,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
            call MPI_BCAST(atemp_all,  proc_len*nprocs,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
            call MPI_BCAST(divtemp_all,proc_len*nprocs,MPI_REAL8,0,MPI_COMM_WORLD,ierr)

            if ( my_id /= 0) then

                ! construct dimension-based policy functions
                v = transform_f(vtemp_all)
                lpol = transform_f(ltemp_all)
                dpol = transform_f(dtemp_all)
                spol = transform_f(stemp_all)
                apol = transform_f(atemp_all)
                div_pol = transform_f(divtemp_all)

            endif

            ! MPI BARRIER
            call MPI_BARRIER(MPI_COMM_WORLD,ierr)

            if (error < tol) then
                exit
            endif

            if (my_id ==0) then
                if (iterator == it_max) then
                    write(*,*) 'Used maxinum number of iterations'
                endif
            endif


        enddo  !vfi loop

    end subroutine

    subroutine get_dist()

        implicit none

        real*8, dimension(bank_nlen,bank_dlen) :: F0, F1, G_ent
        real*8, parameter :: tolerance = 1e-40
        real*8 :: error_temp
        integer :: ii, jj, kk, ll, mm, nn, qq
        integer, parameter :: it_max = 1000
        real*8 :: net_temp
        integer, dimension(1) :: min_idx
        real*8 :: s_tilde, mstar, excess_cash, mstar_agg
        real*8 :: net_int, stock

        integer :: il, ir, iyl, iyr
        real*8  :: phi, phiy

!        default_liq_prob = 0d0
        ! compute default probabilities
        do ii=1,bank_nlen
            do jj=1,bank_dlen

                if ( bank_ngrid(ii) <= nstar) then
                    default_prob(ii,jj)     = 0d0
                    default_liq_prob(ii,jj) = 0d0
                else

                    default_prob(ii,jj)     = 0d0
                    default_liq_prob(ii,jj) = 0d0

                    ! for each liquidity shock
                    do ll=1,size(delta)

                        ! if can't meet funding shock, liquidity default
                        if ( ( pstar*spol(ii,jj)) < delta(ll)*apol(ii,jj) ) then

                            default_liq_prob(ii,jj) = default_liq_prob(ii,jj) + prob_d(ll)

                        else

                            s_tilde = ( delta(ll)*apol(ii,jj) )/pstar

                            ! for each return shock
                            do mm=1,size(Rl)

                                net_int = (Rl(mm)-1d0)*lpol(ii,jj) + i_s*(spol(ii,jj)-s_tilde) &
                                                 - (Rd-1d0)*dpol(ii,jj) - (Ra-1d0)*(1d0-delta(ll))*apol(ii,jj)

                                stock = lpol(ii,jj) + (spol(ii,jj)-s_tilde) -&
                                        dpol(ii,jj) - (1d0-delta(ll))*apol(ii,jj)

                                if ( net_int > 0d0 ) then
                                    net_temp = (1d0-tax)*net_int + stock
                                else
                                    net_temp = net_int + stock
                                endif

                                ! compute default probability
                                if (net_temp <= nstar ) then
                                    default_prob(ii,jj) = default_prob(ii,jj) + &
                                                    prob_l(mm)*prob_d(ll)
                                endif

                            enddo


                        endif

                    enddo
                endif
            enddo
        enddo

        !Initialize problem
        F0 = 1d0/(bank_nlen*bank_dlen)
        error_temp = tolerance + 1d0

        do qq = 1,it_max        !for each iteration in updating process

            F1 = 0d0

            do ii=1,bank_nlen
                do jj=1,bank_dlen

                    ! for each delta shock
                    do ll=1,size(delta)

                        ! if no default, compute liquidations
                        if ( ( pstar*spol(ii,jj)) >= delta(ll)*apol(ii,jj) ) then

                            s_tilde = ( delta(ll)*apol(ii,jj) )/pstar

                            ! for each return shock
                            do mm=1,size(Rl)

                                net_int = (Rl(mm)-1d0)*lpol(ii,jj) + i_s*(spol(ii,jj)-s_tilde) &
                                                 - (Rd-1d0)*dpol(ii,jj) - (Ra-1d0)*(1d0-delta(ll))*apol(ii,jj)

                                stock = lpol(ii,jj) + (spol(ii,jj)-s_tilde) -&
                                        dpol(ii,jj) - (1d0-delta(ll))*apol(ii,jj)

                                if ( net_int > 0d0 ) then
                                    net_temp = (1d0-tax)*net_int + stock
                                else
                                    net_temp = net_int + stock
                                endif

                                ! if no default, record new (networth,deposit level)
                                if (net_temp > nstar ) then

                                    ! find nearest grid point
                                    if ( net_temp > bank_nb ) then
                                        !min_idx(1) = bank_nlen
                                        il = bank_nlen-1
                                        ir = bank_nlen
                                        phi = 0d0
                                    else
                                        !min_idx= minloc( abs( bank_ngrid(jj,:) - net_temp ) )
                                        call linint_Grow(net_temp, bank_na, bank_nb, growth,&
                                                            bank_nlen+1, il, ir, phi)

                                    endif

                                    if ( dpol(ii,jj) > bank_db ) then
                                        !min_idx(1) = bank_nlen
                                        iyl = bank_dlen-1
                                        iyr = bank_dlen
                                        phiy = 0d0
                                    else
                                        !min_idx= minloc( abs( bank_ngrid(jj,:) - net_temp ) )
                                        call linint_Grow(dpol(ii,jj), bank_da, bank_db, dgrowth,&
                                                            bank_dlen+1, iyl, iyr, phiy)

                                    endif

                                    F1(il,iyl) = F1(il,iyl) + phi*phiy*prob_l(mm)*prob_d(ll)*F0(ii,jj)
                                    F1(il,iyr) = F1(il,iyr) + phi*(1d0-phiy)*prob_l(mm)*prob_d(ll)*F0(ii,jj)
                                    F1(ir,iyl) = F1(ir,iyl) + (1d0-phi)*phiy*prob_l(mm)*prob_d(ll)*F0(ii,jj)
                                    F1(ir,iyr) = F1(ir,iyr) + (1d0-phi)*(1d0-phiy)*prob_l(mm)*prob_d(ll)*F0(ii,jj)
                                endif
                            enddo



                        endif


                    enddo

                enddo
            enddo

            ! create entrant distribution (with replacement)
            do ii = 1,bank_nlen
                do jj=1,bank_dlen
                    G_ent(ii,jj) = (default_liq_prob(ii,jj) + default_prob(ii,jj))*F0(ii,jj)
                enddo
            enddo

            mstar_agg = 0d0
            if (sum(G_ent) == 0d0) then
                continue
            else
                G_ent = G_ent/sum(G_ent)  ! normalize to 1
            endif

            !solve for mass of entrants
            mstar = 1d0 - sum(F1)

            mstar_agg = mstar_agg + mstar

            !next-period distribution, with entrants
            F1 = F1 + mstar*G_ent

            error_temp = sum(abs(F1-F0))/(bank_nlen*bank_dlen )

            if (error_temp < tolerance) then
                exit
            else
                F0 = F1
            endif

            if ( (qq == it_max).AND.(my_id==0) ) then
                write(*,*) 'Maximum number of distribution iterations!'
            endif

        enddo

        if (my_id == 0) then
            write(*,*) 'Distribution iteration:', qq
            write(*,*) 'Distribution error:', error_temp
        endif

        F_stationary  = F1
        Agg_Def2 = mstar_agg

    end subroutine

    function implied_ra_rate()

        implicit none

        ! need to compute 4 objects
        !       (1) weighted liquidity default
        !       (2) (1-weighted liquidity default)*liquidity shock
        !       (3) (1-weighted liquidity default)*(1-liquidity shock)*weighted insolvency default
        !       (4) (1-weighted liquidity default)*(1-liquidity shock)*(1-weighted insolvency default)
        !
        !       then Ra = ( (1/beta) - (1+h)(1) - E[(2)] - (1+h)E[3] )/E[(4)]
        !

        real*8 :: implied_ra_rate

        integer :: ii,jj,kk,ll, mm, nn
        real*8, dimension(bank_nlen,bank_dlen) :: bshare
        real*8 :: agg_a
        real*8 :: obj1, obj2, obj3, obj4
        real*8 :: temp_net
        real*8 :: excess_cash
        real*8 :: net_int, stock

        agg_a = 0d0
        ! compute aggregate wholesale lending
        do ii=1,bank_nlen
            do jj=1,bank_dlen
                agg_a = agg_a + F_stationary(ii,jj)*apol(ii,jj)
            enddo
        enddo

        ! compute bank wholesalefunding shares
        do ii=1,bank_nlen
            do jj=1,bank_dlen
                bshare(ii,jj) = F_stationary(ii,jj)*apol(ii,jj)/agg_a
            enddo
        enddo

        ! initialize
        obj1 = 0d0
        obj2 = 0d0
        obj3 = 0d0
        obj4 = 0d0

        ! for each bank type today (n_b,theta,dbar)
        do ii=1,bank_nlen
            do jj=1,bank_dlen

                ! liquidity shocks
                do ll=1,size(delta)

                    if ( delta(ll)*apol(ii,jj) >  pstar*spol(ii,jj) ) then

                        ! computing object (1)
                        obj1 = obj1 + bshare(ii,jj)*prob_d(ll)

                    else

                        ! computing object (2)
                        obj2 = obj2 + bshare(ii,jj)*prob_d(ll)*delta(ll)

                        ! compute possible next-period outcomes
                        do mm=1,size(Rl)

                            net_int = (Rl(mm)-1d0)*lpol(ii,jj) + i_s*(spol(ii,jj)-stilde(ii,jj,ll)) &
                                             - (Rd-1d0)*dpol(ii,jj) - (Ra-1d0)*(1d0-delta(ll))*apol(ii,jj)

                            stock = lpol(ii,jj) + (spol(ii,jj)-stilde(ii,jj,ll)) -&
                                    dpol(ii,jj) - (1d0-delta(ll))*apol(ii,jj)

                            if ( net_int > 0d0 ) then
                                temp_net = (1d0-tax)*net_int + stock
                            else
                                temp_net = net_int + stock
                            endif

                            ! computing objects (3) and (4)
                            if ( temp_net < nstar  ) then  ! insolvency default
                                obj3 = obj3 + bshare(ii,jj)*prob_d(ll)*(1d0-delta(ll))*prob_l(mm)
                            else
                                obj4 = obj4 + bshare(ii,jj)*prob_d(ll)*(1d0-delta(ll))*prob_l(mm)
                            endif

                        enddo


                    endif

                enddo

            enddo
        enddo

!        if (my_id==0) then
!            write(*,*) 'DEF^LIQ=', obj1
!            write(*,*) '(1-DEF^LIQ)*DELTA=',obj2
!            write(*,*) '(1-DEF^LIQ)*DEF*IN*(1-DELTA) = ', obj3
!            write(*,*) '(1-DEF^LIQ)*(1-DEF^IN)*(1-DELTA)=',obj4
!        endif

        ! Ra = ( (1/beta) - (1+h)(1) - E[(2)] - (1+h)E[3] )/E[(4)]

        implied_ra_rate = ( (1d0/beta_pre) - (1d0+hair)*obj1 - obj2 - (1d0+hair)*obj3 )/obj4

    end function


    subroutine grid_f(nprocs,proc_len,idx_grid)

        integer, intent(in) :: nprocs  ! number of processors
        integer, intent(in) :: proc_len  ! length of each processor grid

        integer, dimension(nprocs,proc_len,2), intent(out)  :: idx_grid  ! for each job/id, for each job grid point --> model state space: (nb,theta,dbar)

        ! local variables
        integer :: ii,jj,kk, iterator, proc_idx

        !
        !   Inputs
        !       nprocs: # of processors in parallel code
        !       proc_len: length of each grid per processor
        !
        !   Outputs
        !       idx_grid(i,j,:):  state space index grid in serial dimension, for processor i at parallel grid point j
        !
        !

        ! initialize
        iterator = 1
        proc_idx = 1

        idx_grid = 1

        do ii=1,bank_nlen                    ! for each level of bank networth
            do jj=1,bank_dlen            ! each deposit capacity constraint

                ! record serial grid index
                idx_grid(proc_idx,iterator,:) = (/ ii,jj /)

                if (iterator == proc_len) then  ! if iteration is last of processor length
                    iterator = 1                    ! reset iteration
                    proc_idx = proc_idx+1           ! update the processor id
                else
                    iterator = iterator + 1
                endif

            enddo
        enddo

    end subroutine

    function transform_f(total_array)

        real*8, dimension(:), intent(in) :: total_array

        real*8, dimension(bank_nlen,bank_dlen) :: transform_f

        ! local variables
        integer :: ii,jj,kk, iterator

        !
        !   Inputs
        !       total_array: parallel version of policy function
        !
        !   Outputs
        !       transform_f: serial version of policy function # of dimensions = # of states
        !
        !

        iterator = 1

        do ii=1,bank_nlen                   ! for each level of bank networth
            do jj=1,bank_dlen           ! for each deposit capacity constraint

                transform_f(ii,jj) = total_array(iterator)

                iterator = iterator + 1

            enddo
        enddo

    end function

    subroutine output()

        implicit none

        real*8 :: total_assets, total_liabs, total_loan, total_sec, total_debt, total_whole, total_net,total_eq,&
                  total_rwa, elasticity, roe1_num, roe1_den, roe2_num, roe2_den, s_tilde, net_int, stock, net_temp, &
                  hold, weighted_spread_num, weighted_spread_den

        real*8, dimension(2) :: roe_size_corr, Agg_Def_Marg, port_shares
        real*8, dimension(3) :: debt_shares
        real*8, dimension(4) :: al_sheet, holder

        real*8, dimension(size(data_moms)) :: mod_moms
        real*8, dimension(6) :: size_corr
        real*8, dimension(7) :: size_corr_intra
        real*8, dimension(8,bank_nlen*bank_dlen) :: corr_objs  ! distribution weight, size, risk-weighted equity, leverage, liquidity, ins default, liq default, equity value
        real*8, dimension(9,bank_nlen*bank_dlen) :: corr_objs_intra
        real*8, dimension(2, bank_nlen*bank_dlen ) :: roe_stack
        real*8, dimension(bank_nlen,bank_dlen) :: roe_pol1, roe_pol2, lev_pol, avg_net

        integer :: ii, jj, kk, ll, mm, iterator

        if (my_id==0) then

            !~~~~~~~~~~~!
            !           !
            !   Plots   !
            !           !
            !~~~~~~~~~~~!

            ! deposit policy functions with different dbar
            do ii=1,bank_dlen
                call plot(bank_ngrid,dpol(:,ii),marker=2)
            enddo
            call execplot(title='Deposit Policy Functions',filename='deposit',filetype='png',output='deposit')

            ! loan policy functions with different dbar
            do ii=1,bank_dlen
                call plot(bank_ngrid,lpol(:,ii),marker=2)
            enddo
            call execplot(title='Loan Policy Functions',filename='loan',filetype='png',output='loan')

            ! wholesale policy functions with different dbar
            do ii=1,bank_dlen
                call plot(bank_ngrid,apol(:,ii),marker=2)
            enddo
            call execplot(title='Wholesale Policy Functions',filename='wholesale',filetype='png',output='wholesale')

            ! dividend policy functions with different dbar
            do ii=1,bank_dlen
                call plot(bank_ngrid,div_pol(:,ii),marker=2)
            enddo
            call execplot(title='Dividend Policy Functions',filename='dividend',filetype='png',output='dividend')

            ! stationary distribution
            do ii=1,bank_dlen
                call plot_hist(bank_ngrid,F_stationary(:,ii))
            enddo
            call execplot(title='Stationary Distribution of Networth' ,filename='dist',filetype='png',output='dist')

            ! leverage equity ratio policy function
            do ii=1,bank_nlen
                do jj=1,bank_dlen

                    lev_pol(ii,jj) = ( lpol(ii,jj) + spol(ii,jj) - dpol(ii,jj) - apol(ii,jj) )/(lpol(ii,jj) + spol(ii,jj))

                enddo
            enddo

            do ii=1,bank_dlen
                call plot(bank_ngrid,lev_pol(:,ii),marker=2)
            enddo
            call execplot(title='Leverage Ratio Policy Function',filename='lev',filetype='png',output='lev')

            ! average networth law of motion
            do ii=1,bank_nlen
                do jj=1,bank_dlen
                    avg_net(ii,jj) = l_mu*lpol(ii,jj) + (1d0+i_s)*spol(ii,jj) -&
                                        Rd*dpol(ii,jj) - Ra*apol(ii,jj)
                enddo
            enddo

            do ii = 1,bank_dlen
                call plot(bank_ngrid,avg_net(:,ii))
            enddo
            call plot(bank_ngrid,bank_ngrid,linewidth=4d0)
            call execplot(title='Networth Law of Motion',filename='net_law',filetype='png',output='net_law')

            ! deposit law of motion
            do ii = 1,bank_nlen
                call plot(bank_dgrid,dpol(ii,:))
            enddo
            call plot(bank_dgrid,bank_dgrid,linewidth=4d0)
            call execplot(title='Deposits Law of Motion',filename='dep_law',filetype='png',output='dep_law')

            !~~~~~~~~~~~~~~~~~~~~~~~!
            !                       !
            !   TARGETED MOMENTS    !
            !                       !
            !~~~~~~~~~~~~~~~~~~~~~~~!

            mod_moms = model_moments()

            write(*,*)
            write(*,*) 'TARGETED MOMENTS'
            write(*,*)
            write(*,*) '                               Model                   Data'
            write(*,*)
            write(*,*) 'loan-to-security ratio:    ', mod_moms(1),'        ',data_moms(1)
            write(*,*)
            write(*,*) 'deposit-to-wholesale ratio:', mod_moms(2),'        ',data_moms(2)
            write(*,*)
            write(*,*) 'risk-weighted equity ratio:', mod_moms(3),'        ',data_moms(3)
            write(*,*)
            write(*,*) 'insolvency default rate:   ', mod_moms(4),'        ',data_moms(4)
            write(*,*)
            write(*,*) 'liquidity ratio:           ', mod_moms(5),'        ',data_moms(5)
            write(*,*)
            write(*,*) 'deposit shares  :          ', mod_moms(6),'        ',data_moms(6)
            write(*,*)


            !~~~~~~~~~~~~~~~~~~~~~~~~~!
            !                         !
            !   UNTARGETED MOMENTS    !
            !                         !
            !~~~~~~~~~~~~~~~~~~~~~~~~~!

            al_sheet = bank_shares()

            write(*,*) 'UNTARGETED MOMENTS'
            write(*,*)
            write(*,*) 'BANK BALANCE SHEET'
            write(*,*) '                 Model                 Data'
            write(*,*)
            write(*,*) 'loans:     ', 100d0*al_sheet(1)               ,'        ',70d0
            write(*,*)
            write(*,*) 'securities:', 100d0*al_sheet(2)               ,'        ',25d0
            write(*,*)
            write(*,*) 'deposits:  ', 100d0*al_sheet(3)               ,'        ',75d0
            write(*,*)
            write(*,*) 'wholesale: ', 100d0*al_sheet(4)               ,'        ',20d0
            write(*,*)
            write(*,*) 'equity:    ', 100d0*(1d0 - al_sheet(3)-al_sheet(4)),'        ',5d0
            write(*,*)

            ! compute leverage ratios and size correlations with (i) equity ratios and (ii) liquidity ratios
            total_eq  = 0d0
            total_rwa = 0d0

            do ii=1,bank_nlen
                do jj=1,bank_dlen


                    total_eq  = total_eq  + F_stationary(ii,jj)*( lpol(ii,jj)  +&
                                spol(ii,jj) - (dpol(ii,jj)+apol(ii,jj)))

                    total_rwa = total_rwa + F_stationary(ii,jj)*( lpol(ii,jj) + &
                                         spol(ii,jj) )

                enddo
            enddo

            write(*,*)
            write(*,*) '                 Model                 Data'
            write(*,*)
            write(*,*) 'Leverage ratio:', 100d0*total_eq/total_rwa,'       ',7.25d0
            write(*,*)

            ! compute spreads (1) Ra-Rd (2) i_s - Rd (3) e[Rl] - i_s
            weighted_spread_num = 0d0
            weighted_spread_den = 0d0

            do ii=1,bank_nlen
                do ll=1,size(Rl)
                    do kk=1,bank_dlen


                        if ( lpol(ii,kk) .ne. 0d0 ) then

                            weighted_spread_num = weighted_spread_num + F_stationary(ii,kk)*&
                                                Rl(ll)*lpol(ii,kk)

                            weighted_spread_den = weighted_spread_den + F_stationary(ii,kk)*&
                                                (lpol(ii,kk) + g(lpol(ii,kk)) )



                        endif

                    enddo
                enddo
            enddo

            write(*,*)
            write(*,*) 'Equilibrium Spreads (Annualized)'
            write(*,*)
            write(*,*) 'ra-rd:', 100d0*4d0*(Ra-Rd)
            write(*,*) 'is-rd:', 100d0*4d0*(1d0+i_s-Rd)
            write(*,*) 'E[rl]:', 100d0*4d0*(weighted_spread_num/weighted_spread_den-1d0)
            write(*,*) 'E[rl]-is:', 100d0*4d0*(weighted_spread_num/weighted_spread_den - (1d0+i_s) )
            write(*,*)

            ! size correlations with (1) risk-weighted equity (2) leverage ratio and (3) liquidity ratio

            ! stack correlation objects
            corr_objs = 0d0
            iterator = 1

            do ii=1,bank_nlen
                do jj=1,bank_dlen

                    ! distribution weights
                    corr_objs(1,iterator) =  F_stationary(ii,jj)

                    ! size measure
                    corr_objs(2,iterator) =  lpol(ii,jj)

                    ! risk-weighted equity ratios
                    if ( lpol(ii,jj) .ne. 0d0 ) then
                        corr_objs(3,iterator) =  (lpol(ii,jj) + spol(ii,jj) - &
                                                        (dpol(ii,jj)+apol(ii,jj)))/lpol(ii,jj)
                    else
                        corr_objs(3,iterator) = ebar
                    endif

                    ! leverage ratio
                    if ( ( lpol(ii,jj) + spol(ii,jj) ) .ne. 0d0 ) then
                        corr_objs(4,iterator) =  (lpol(ii,jj) + spol(ii,jj) - &
                                    (dpol(ii,jj)+apol(ii,jj)))/( lpol(ii,jj) + spol(ii,jj) )
                    else
                        corr_objs(4,iterator) =  ebar
                    endif

                    ! liquidity ratio
                    if ( apol(ii,jj) .ne. 0d0) then
                        corr_objs(5,iterator) =  ( (1d0-lr_hair)*spol(ii,jj))/( apol(ii,jj) )
                    else
                        corr_objs(5,iterator) =  phi_lr
                    endif

                    iterator = iterator + 1

                enddo
            enddo

            size_corr(1) = weighted_corr( corr_objs(1,:),corr_objs(2,:),corr_objs(3,:) )  ! size and risk-weighted equity
            size_corr(2) = weighted_corr( corr_objs(1,:),corr_objs(2,:),corr_objs(4,:) )  ! size and leverage
            size_corr(3) = weighted_corr( corr_objs(1,:),corr_objs(2,:),corr_objs(5,:) )  ! size and liquidity ratio

            write(*,*)
            write(*,*) 'Size correlation with risk-weighted equity ratios:', size_corr(1),'       ',0d0
            write(*,*) 'Size correlation with leverage ratios:', size_corr(2),'       ',0d0
            write(*,*) 'Size correlation with liquidity ratios:', size_corr(3),'       ',0d0
            write(*,*)

            !~~~~~~~~~~~~~~~~~~~~~~~~~!
            !                         !
            !   UNOBSERVED MOMENTS    !
            !                         !
            !~~~~~~~~~~~~~~~~~~~~~~~~~!

            write(*,*) 'UNOBSERVED MOMENTS'
            write(*,*)

            ! liquidity default
            Agg_Def = Aggregate_Def()

            ! firesale elasticity
            elasticity = fire_elastic(Agg_Liq)


            write(*,*) 'Liquidity Default Rate:', 100d0*100d0*Agg_Def(2)
            write(*,*)
            write(*,*)
            write(*,*) 'Elasticity of Firesold Securities:', elasticity
            write(*,*)

            ! size correlations with (1) insolvency default (2) liquidity default and (3) equity value

            ! stack correlation objects
            iterator = 1

            do ii=1,bank_nlen
                do jj=1,bank_dlen

                    ! insolvency default
                    corr_objs(6,iterator) = default_prob(ii,jj)

                    ! liquidity default
                    corr_objs(7,iterator) = default_liq_prob(ii,jj)

                    iterator = iterator + 1

                enddo
            enddo

            size_corr(4) = weighted_corr( corr_objs(1,:),corr_objs(2,:),corr_objs(6,:) )  ! size and insolvency default
            size_corr(5) = weighted_corr( corr_objs(1,:),corr_objs(2,:),corr_objs(7,:) )  ! size and liquidity default

            write(*,*)
            write(*,*) 'Size correlation with insolvency default:', size_corr(4)
            write(*,*) 'Size correlation with liquidity default:', size_corr(5)
            write(*,*)


            ! return on equity measures
            !
            !   ROE1: E[n']/n
            !
            !   ROE2: E[n']/( l + s + c - a - d)
            !
            roe1_num = 0d0
            roe1_den = 0d0
            roe2_num = 0d0
            roe2_den = 0d0
            roe_pol1 = 0d0
            roe_pol2 = 0d0

            do ii=1,bank_nlen
                do jj=1,bank_dlen

                    ! for each liquidity shock
                    do ll=1,size(delta)

                        if ( pstar*spol(ii,jj)  >= delta(ll)*apol(ii,jj)) then

                            s_tilde = ( delta(ll)*apol(ii,jj) )/pstar

                            ! for each loan shock
                            do mm=1,size(Rl)

                                ! compute net interest
                                net_int = (Rl(mm)-1d0)*lpol(ii,jj) + i_s*(spol(ii,jj)-s_tilde) - (Rd-1d0)*dpol(ii,jj) &
                                                        - (Ra-1d0)*(1d0-delta(ll))*apol(ii,jj)

                                stock = lpol(ii,jj) + (spol(ii,jj)-s_tilde) - dpol(ii,jj) &
                                                        - (1d0-delta(ll))*apol(ii,jj)

                                if ( net_int > 0d0) then
                                    net_temp = (1d0-tax)*net_int + stock
                                else
                                    net_temp = net_int + stock
                                endif

                                if (net_temp >= 0d0) then

                                    ! record return on equity
                                    roe1_num = roe1_num + prob_l(mm)*prob_d(ll)*F_stationary(ii,jj)*&
                                                net_temp
                                    roe1_den = roe1_den + prob_l(mm)*prob_d(ll)*F_stationary(ii,jj)*&
                                                ( bank_ngrid(ii) - div_pol(ii,jj) )


                                    if ( ( bank_ngrid(ii) - div_pol(ii,jj) ) .ne. 0d0 ) then

                                        roe_pol1(ii,jj) = roe_pol1(ii,jj) + prob_l(mm)*prob_d(ll)*(net_temp/&
                                                            ( bank_ngrid(ii) - div_pol(ii,jj) ) -1d0)
                                    endif

                                    roe2_num = roe2_num + prob_l(mm)*prob_d(ll)*F_stationary(ii,jj)*&
                                                net_temp
                                    roe2_den = roe2_den + prob_l(mm)*prob_d(ll)*F_stationary(ii,jj)*&
                                        ( lpol(ii,jj) + spol(ii,jj) - dpol(ii,jj) - apol(ii,jj) )

                                    if ( ( lpol(ii,jj) + spol(ii,jj) - &
                                            dpol(ii,jj) - apol(ii,jj) ) .ne. 0d0 ) then

                                        roe_pol2(ii,jj) = roe_pol2(ii,jj) + prob_l(mm)*prob_d(ll)*&
                                            (net_temp/( lpol(ii,jj) + spol(ii,jj) - &
                                            dpol(ii,jj) - apol(ii,jj) ) -1d0)

                                    endif



                                endif

                            enddo
                        endif
                    enddo
                enddo
            enddo

            ! return on equity correlation with (1) size
            ! stack correlation objects
            iterator = 1

            do ii=1,bank_nlen
                do jj=1,bank_dlen

                    ! insolvency default
                    roe_stack(1,iterator) = roe_pol1(ii,jj)

                    ! liquidity default
                    roe_stack(2,iterator) = roe_pol2(ii,jj)

                    iterator = iterator + 1

                enddo
            enddo

            roe_size_corr(1) = weighted_corr( corr_objs(1,:),corr_objs(2,:),roe_stack(1,:) )  ! size and return on equity
            roe_size_corr(2) = weighted_corr( corr_objs(1,:),corr_objs(2,:),roe_stack(2,:) )  ! size and return on equity

            write(*,*)
            write(*,*) 'Return on Equity (Measure 1):', 100d0*( roe1_num/roe1_den -1d0)
            write(*,*) 'Size correlation:', roe_size_corr(1)
            write(*,*)
            write(*,*) 'Return on Equity (Measure 2):', 100d0*( roe2_num/roe2_den- 1d0)
            write(*,*) 'Size correlation:', roe_size_corr(2)
            write(*,*)


            write(*,*) 'other cross-section correlations'
            write(*,*)

            ! risk-weighted and leverage (3 + 4)
            hold = weighted_corr( corr_objs(1,:),corr_objs(3,:),corr_objs(4,:))
            write(*,*) 'risk-weighted equity and leverage:', hold

            ! risk-weighted and liquidity ( 3 + 5)
            hold = weighted_corr( corr_objs(1,:),corr_objs(3,:),corr_objs(5,:))
            write(*,*) 'risk-weighted equity and liquidity:', hold

            ! risk-weighted and insolvency default (3 + 6)
            hold = weighted_corr( corr_objs(1,:),corr_objs(3,:),corr_objs(6,:))
            write(*,*) 'risk-weighted equity and insolvency default:', hold

            ! risk-weighted and liquidity default (3 + 7)
            hold = weighted_corr( corr_objs(1,:),corr_objs(3,:),corr_objs(7,:))
            write(*,*) 'risk-weighted equity and liquidity default:', hold

            ! risk-weighted and roe (3 + roe 2)
            hold = weighted_corr( corr_objs(1,:),corr_objs(3,:),roe_stack(2,:))
            write(*,*) 'risk-weighted equity and roe:', hold

            ! leverage and liquidity ratio ( 4 + 5)
            hold = weighted_corr( corr_objs(1,:),corr_objs(4,:),corr_objs(5,:))
            write(*,*) 'leverage and liquidity:', hold

            ! leverage and insolvency default ( 4 + 6)
            hold = weighted_corr( corr_objs(1,:),corr_objs(4,:),corr_objs(6,:))
            write(*,*) 'leverage and insolvency default:', hold

            ! leverage and liquidity defa (4 + 7)
            hold = weighted_corr( corr_objs(1,:),corr_objs(4,:),corr_objs(7,:))
            write(*,*) 'leverage and liquidity default:', hold

            ! leverage and roe (4 + roe 2)
            hold = weighted_corr( corr_objs(1,:),corr_objs(4,:),roe_stack(2,:))
            write(*,*) 'leverage and roe:', hold

            ! liquidity ratio and insolvency default (5 + 6)
            hold = weighted_corr( corr_objs(1,:),corr_objs(5,:),corr_objs(6,:))
            write(*,*) 'liq ratio and ins default:', hold

            ! liquidity ratio and liq default (5 + 7)
            hold = weighted_corr( corr_objs(1,:),corr_objs(5,:),corr_objs(7,:))
            write(*,*) 'liq ratio and liq default:', hold

            ! liquidity ratio and roe (5 + roe 2)
            hold = weighted_corr( corr_objs(1,:),corr_objs(5,:),roe_stack(2,:))
            write(*,*) 'liq ratio and roe:', hold

            ! insolvency default and liquidity default ( 6 + 7)
            hold = weighted_corr( corr_objs(1,:),corr_objs(6,:),corr_objs(7,:))
            write(*,*) 'insolvency default and liquidity default:', hold

            ! insolvency default and roe ( 6 + roe 2)
            hold = weighted_corr( corr_objs(1,:),corr_objs(6,:),roe_stack(2,:))
            write(*,*) 'insolvency default and roe:', hold

            ! liquidity default and roe ( 7 + roe 2)
            hold = weighted_corr( corr_objs(1,:),corr_objs(7,:),roe_stack(2,:))
            write(*,*) 'liquidity default and roe:', hold


            ! compute intra-bank type correlation measures for size and roe
            corr_objs_intra = 0d0
            iterator = 1

            do ii=1,bank_nlen
                do kk=1,bank_dlen

                    ! distribution weights
                    corr_objs_intra(1,iterator) =  F_stationary(ii,kk)

                    ! size measure
                    corr_objs_intra(2,iterator) =  lpol(ii,kk)

                    ! risk-weighted equity ratios
                    if ( lpol(ii,kk) .ne. 0d0 ) then
                        corr_objs_intra(3,iterator) =  (lpol(ii,kk) + spol(ii,kk) - &
                                                        (dpol(ii,kk)+apol(ii,kk)))/lpol(ii,kk)
                    else
                        corr_objs_intra(3,iterator) = ebar
                    endif

                    ! leverage ratio
                    if ( ( lpol(ii,kk) + spol(ii,kk) ) .ne. 0d0 ) then
                        corr_objs_intra(4,iterator) =  (lpol(ii,kk) + spol(ii,kk) - &
                                    (dpol(ii,kk)+apol(ii,kk)))/( lpol(ii,kk) + spol(ii,kk) )
                    else
                        corr_objs_intra(4,iterator) =  ebar
                    endif

                    ! liquidity ratio
                    if ( apol(ii,kk) .ne. 0d0) then
                        corr_objs_intra(5,iterator) =  ( (1d0-lr_hair)*spol(ii,kk))/( apol(ii,kk) )
                    else
                        corr_objs_intra(5,iterator) =  phi_lr
                    endif

                    ! insolvency default
                    corr_objs_intra(6,iterator) = default_prob(ii,kk)

                    ! liquidity default
                    corr_objs_intra(7,iterator) = default_liq_prob(ii,kk)

                    ! roe 1
                    corr_objs_intra(8,iterator) = roe_pol1(ii,kk)

                    ! liquidity default
                    corr_objs_intra(9,iterator) = roe_pol2(ii,kk)

                    iterator = iterator + 1

                enddo
            enddo

            ! record
            size_corr_intra(1) = weighted_corr( corr_objs_intra(1,:), corr_objs_intra(2,:), corr_objs_intra(3,:) ) ! size and risk weighted equity
            size_corr_intra(2) = weighted_corr( corr_objs_intra(1,:), corr_objs_intra(2,:), corr_objs_intra(4,:) ) ! leverage
            size_corr_intra(3) = weighted_corr( corr_objs_intra(1,:), corr_objs_intra(2,:), corr_objs_intra(5,:) ) ! liqudity
            size_corr_intra(4) = weighted_corr( corr_objs_intra(1,:), corr_objs_intra(2,:), corr_objs_intra(6,:) ) ! insolv def
            size_corr_intra(5) = weighted_corr( corr_objs_intra(1,:), corr_objs_intra(2,:), corr_objs_intra(7,:) ) ! liq default
            size_corr_intra(6) = weighted_corr( corr_objs_intra(1,:), corr_objs_intra(2,:), corr_objs_intra(8,:) ) ! roe 1
            size_corr_intra(7) = weighted_corr( corr_objs_intra(1,:), corr_objs_intra(2,:), corr_objs_intra(9,:) ) ! roe 2

            write(*,*)
            write(*,*) 'Bank Size Correlations'
            write(*,*) 'risk-weighted equity ratio:', size_corr_intra(1)
            write(*,*) 'leverage ratio:', size_corr_intra(2)
            write(*,*) 'liquidity ratio:', size_corr_intra(3)
            write(*,*) 'insolvency default:', size_corr_intra(4)
            write(*,*) 'liquidity default:', size_corr_intra(5)
            write(*,*) 'ROE 1:', size_corr_intra(6)
            write(*,*) 'ROE 2:', size_corr_intra(7)

            hold = weighted_corr( corr_objs_intra(1,:), corr_objs_intra(3,:), corr_objs_intra(4,:) )
            write(*,*) 'risk weighted and leverage:',hold
            hold = weighted_corr( corr_objs_intra(1,:), corr_objs_intra(3,:), corr_objs_intra(5,:) )
            write(*,*) 'risk weighted and liquidity:',hold
            hold = weighted_corr( corr_objs_intra(1,:), corr_objs_intra(3,:), corr_objs_intra(6,:) )
            write(*,*) 'risk weighted and ins default:',hold
            hold = weighted_corr( corr_objs_intra(1,:), corr_objs_intra(3,:), corr_objs_intra(7,:) )
            write(*,*) 'risk weighted and liq default:',hold
            hold = weighted_corr( corr_objs_intra(1,:), corr_objs_intra(3,:), corr_objs_intra(8,:) )
            write(*,*) 'risk weighted and roe1:',hold
            hold = weighted_corr( corr_objs_intra(1,:), corr_objs_intra(3,:), corr_objs_intra(9,:) )
            write(*,*) 'risk weighted and roe2:',hold

            hold = weighted_corr( corr_objs_intra(1,:), corr_objs_intra(4,:), corr_objs_intra(5,:) )
            write(*,*) 'leverage and liquidity:',hold
            hold = weighted_corr( corr_objs_intra(1,:), corr_objs_intra(4,:), corr_objs_intra(6,:) )
            write(*,*) 'leverage and ins default:',hold
            hold = weighted_corr( corr_objs_intra(1,:), corr_objs_intra(4,:), corr_objs_intra(7,:) )
            write(*,*) 'leverage and liq default:',hold
            hold = weighted_corr( corr_objs_intra(1,:), corr_objs_intra(4,:), corr_objs_intra(8,:) )
            write(*,*) 'leverage and roe1:',hold
            hold = weighted_corr( corr_objs_intra(1,:), corr_objs_intra(4,:), corr_objs_intra(9,:) )
            write(*,*) 'leverage and roe2:',hold

            hold = weighted_corr( corr_objs_intra(1,:), corr_objs_intra(5,:), corr_objs_intra(6,:) )
            write(*,*) 'liquidity and ins default:',hold
            hold = weighted_corr( corr_objs_intra(1,:), corr_objs_intra(5,:), corr_objs_intra(7,:) )
            write(*,*) 'liquidity and liq default:',hold
            hold = weighted_corr( corr_objs_intra(1,:), corr_objs_intra(5,:), corr_objs_intra(8,:) )
            write(*,*) 'liquidity and roe1:',hold
            hold = weighted_corr( corr_objs_intra(1,:), corr_objs_intra(5,:), corr_objs_intra(9,:) )
            write(*,*) 'liquidity and roe2:',hold

            hold = weighted_corr( corr_objs_intra(1,:), corr_objs_intra(6,:), corr_objs_intra(7,:) )
            write(*,*) 'ins default and liq default:',hold
            hold = weighted_corr( corr_objs_intra(1,:), corr_objs_intra(6,:), corr_objs_intra(8,:) )
            write(*,*) 'ins default and roe1:',hold
            hold = weighted_corr( corr_objs_intra(1,:), corr_objs_intra(6,:), corr_objs_intra(9,:) )
            write(*,*) 'ins default and roe2:',hold

            hold = weighted_corr( corr_objs_intra(1,:), corr_objs_intra(7,:), corr_objs_intra(8,:) )
            write(*,*) 'liq default and roe1:',hold
            hold = weighted_corr( corr_objs_intra(1,:), corr_objs_intra(7,:), corr_objs_intra(9,:) )
            write(*,*) 'liq default and roe2:',hold

            hold = weighted_corr( corr_objs_intra(1,:), corr_objs_intra(8,:), corr_objs_intra(9,:) )
            write(*,*) 'roe1 and roe2:',hold


            ! bank type portfolio shares
            port_shares = 0d0

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

            port_shares(1) = total_loan/total_assets
            port_shares(2) = total_sec/total_assets

            debt_shares(1) = total_debt/total_liabs
            debt_shares(2) = total_whole/total_liabs
            debt_shares(3) = total_net/total_liabs

            write(*,*)
            write(*,*) 'BANK-LEVEL STATISTICS'
            write(*,*)
            write(*,*) 'Loan portfolio share (%):', 100d0*port_shares(1)
            write(*,*) 'Security portfolio share (%):', 100d0*port_shares(2)
            write(*,*) 'Deposit liability share (%):', 100d0*debt_shares(1)
            write(*,*) 'Wholesale liability share (%):', 100d0*debt_shares(2)
            write(*,*) 'Networth liabiloty share (%):', 100d0*debt_shares(3)
            write(*,*)

            ! record aggregates for comparison to Dodd-Frank action
            Agg_DI  = Aggregate_DI()
            Agg_Div = Aggregate_Div()
            Agg_Gov = Aggregate_Gov()
            Agg_Dep = Aggregate_Dep()
            Agg_Liq = Aggregate_Liq()
            Agg_Corp = Aggregate_Corp_Tax()

            holder    = Aggregate_BS()
            Agg_Loan  = holder(1)
            Agg_Sec   = holder(2)
            Agg_Dep   = holder(3)
            Agg_Whole = holder(4)


            write(*,*) 'Aggregates for Comparison to Dodd-Frank Outcomes'
            write(*,*)
            ! agg lending
            write(*,*) 'Aggregate Lending:', Agg_Loan

            ! agg balance sheet (l + s + c)
            write(*,*) 'Aggregate Bank Balance Sheet Size:', Agg_Loan + Agg_Sec
            write(*,*) 'Aggregate Balance Sheet Liquidity (%):', 100d0*(Agg_Sec )/( Agg_Loan + Agg_Sec )
            write(*,*) 'Aggregate Wholesale Funding:', Agg_Whole

            ! agg consumption
            write(*,*) 'Household Consumption:', (Rd-1d0)*Agg_Dep + Agg_Div(1) + Agg_Div(2) - Agg_Gov - Agg_DI + Agg_Corp
            write(*,*)

        endif

    end subroutine


    subroutine output_elasticities()

        implicit none

        integer :: ii,jj,kk
        real*8 :: ebar_new, lcr_new

        ! compute 1% increase in dodd-frank cap, liq reqs
        ebar_new = 1.01d0*.06d0
        lcr_new = 1.01*100d0

        write(*,*)
        write(*,*) 'Output Elasticities for Current Dodd-Frank Regulations'
        write(*,*)

    end subroutine

    function weighted_corr(f,x,y)

        real*8, dimension(:), intent(in) :: x,y,f
        real*8 :: weighted_corr

        real*8 :: mux, muy, sx, sy, sxy

        ! compute means
        mux = sum( f*x )
        muy = sum( f*y )

        ! compute vols
        sx = ( sum( f*(x-mux)**2d0 ) )**(0.5d0)
        sy = ( sum( f*(y-muy)**2d0 ) )**(0.5d0)

        ! compute covariance
        sxy = sum( f*(x-mux)*(y-muy) )

        ! compute correlation
        weighted_corr = sxy/(sx*sy)

    end function

    subroutine f_order(fvals,xvals)

        implicit none

        real*8, intent(inout), dimension(3)   :: fvals
        real*8, intent(inout), dimension(3,2) :: xvals

        integer :: ii
        integer, dimension(1) :: idx_hold
        real*8, dimension(3)  :: fval_hold
        real*8, dimension(3,2)  ::   xval_hold

        do ii=1,3

            idx_hold = minloc(fvals)

            fval_hold(ii)   = fvals(idx_hold(1))
            xval_hold(ii,:) = xvals(idx_hold(1),:)

            fvals(idx_hold(1)) = 100000000000d0

        enddo

        fvals = fval_hold
        xvals = xval_hold

    end subroutine


END PROGRAM compute_eq
