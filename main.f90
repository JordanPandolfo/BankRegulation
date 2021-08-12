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

    dbar_mu    = 0.0336d0
    gam        = 0.96d0
    pstar      = 0.97d0
    l_sigma2   = 0.000416d0
    i_s        = .0039d0

    beta = gam*beta_pre

    alpha = 0.0002126d0
    omega_outside = (no_liq_price/alpha)**(1d0/(alpha-1d0))

    bank_nb   = .005d0 !dbar_mu*1.45188876d0*( (l_mu + 2.5d0*l_sigma2**(0.5d0)) - Rd*(1d0-.05d0)  )


    delta_l   = 20d0*( bank_nb )/bank_nlen

    delta_div = delta_l*.1d0
    delta_a   = delta_l*1d0
    delta_d   = delta_l*1d0
    delta_s   = delta_l*1d0
    delta_c   = delta_l*.25d0

    ! MPI BARRIER
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
    !                                   !
    !   Solve Stationary Equilibrium    !
    !                                   !
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!

    grid_len = 5
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

        ! initialize value function
        do ii = 1,bank_nlen
            do kk=1,dbar_size
                v(ii,kk) =  bank_ngrid(ii)
            enddo
        enddo

        do ii=1,bank_nlen
            temp(ii) = ii
        enddo

        ! create identity matrix id_mat
        id_mat = 0d0
        do ii=1,bank_nlen*dbar_size
            do jj=1,bank_nlen*dbar_size

                if (ii==jj) then
                    id_mat(ii,jj) = 1d0
                endif

            enddo
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

        !call normal_discrete(Rl,prob_l,l_mu,l_sigma2)

        ! set deposit capacity constraints and probabilities
        dbar = dbar_mu*( 1d0 + (/ -0.45188876,-0.27113325,-0.09037775,0.09037775,0.27113325,0.45188876 /) )

        dbar_prob(1,:) = (/ 3.6355770e-01,4.0356710e-01,1.9755536e-01,3.3369986e-02,1.9130819e-03,3.6768955e-05 /)
        dbar_prob(2,:) = (/ 1.4184205e-01,3.6066970e-01,3.5845480e-01,1.2376471e-01,1.4674066e-02,5.9467011e-04 /)
        dbar_prob(3,:) = (/ 3.6311337e-02,2.0043165e-01,4.0441560e-01,2.8389830e-01,6.9047132e-02,5.8959822e-03 /)
        dbar_prob(4,:) = (/ 5.8959822e-03,6.9047132e-02,2.8389830e-01,4.0441560e-01,2.0043165e-01,3.6311337e-02 /)
        dbar_prob(5,:) = (/ 5.9467011e-04,1.4674066e-02,1.2376471e-01,3.5845480e-01,3.6066970e-01,1.4184205e-01 /)
        dbar_prob(6,:) = (/ 3.6768955e-05,1.9130819e-03,3.3369986e-02,1.9755536e-01,4.0356710e-01,3.6355770e-01 /)


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
        real*8, dimension(bank_nlen,dbar_size,size(delta)) :: ctilde, stilde
        real*8 :: excess_cash

        do ii=1,bank_nlen
            do kk=1,dbar_size

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
        integer, parameter :: it_max = 300

        ! Solving bank problem
        integer :: ii, jj, kk, ll, mm, nn, oo, pp, qq,rr,ss
        real*8 :: v_temp, RHS, net_temp, e_temp
        real*8, dimension(6) :: v_curr_max, l_curr_max, d_curr_max, div_curr_max, s_curr_max, a_curr_max
        real*8 :: a_bound, l_bound, s_bound, d_bound, liq_temp, implied_a, &
                  implied_s, implied_c, d_ubound, implied_d, implied_l, implied_div
        real*8 :: s_tilde  ! liquidation values
        integer, dimension(1) :: v_idx
        real*8 :: net_int, stock, pi_bound
        real*8 :: temp_pen

        ! linear interpolation variables
        integer :: ixl, ixr
        real*8  :: varphi

        ! spline interpolation variables
        real*8, dimension( dbar_size, bank_nlen+2 ) :: coeff

        !---------------!
        !   mpi stuff   !
        !---------------!
        real*8, allocatable, dimension(:) :: ltemp, dtemp, stemp, atemp, divtemp, vtemp, ctemp
        real*8, allocatable, dimension(:) :: ltemp_all, dtemp_all, stemp_all, atemp_all, divtemp_all, vtemp_all, ctemp_all
        integer, dimension(2) :: grid_idx   ! 2 states (nb,dbar)
        integer, allocatable, dimension(:,:,:) :: idx_full_grid

        integer :: total_grid, proc_len

        ! total size of state space grid
        total_grid = bank_nlen*dbar_size

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
            do ii=1,proc_len  ! for each point in parallel grid; state space = (bank ngrid, theta, dbar)

                ! determine indices for bank state space
                grid_idx = idx_full_grid(my_id+1,ii,:)
                !
                !
                !    Case 1: binding capital requirement
                !
                !    implication:   I) d = ( l + s )*(1-ebar)
                !
                !                  II) s = ( n - theta*l^2/2 - div - ebar*l )/ebar
                !
                !                 III) a = 0
                !
                !    search over: (l,div)
                !
                !

                ! construct grids for control variables
                if ( ( grid_idx(1) == 1 ).or.( iterator == 1 ) ) then  ! if n_b(1) or first iteratior of value function, use default grids

                    call grid_Cons_Equi( lgrid, la,lb   )
                    call grid_Cons_Equi( divgrid,diva, divb )

                else

                    l_bound = maxval( (/ 0.0000001d0, lpol(grid_idx(1)-1,grid_idx(2)) - back_look*delta_l /) )
                    pi_bound = maxval( (/ 0d0, div_pol(grid_idx(1)-1,grid_idx(2)) - delta_div /) )

                    call grid_Cons_Equi( lgrid, l_bound, lpol(grid_idx(1)-1,grid_idx(2)) + delta_l )
                    call grid_Cons_Equi( divgrid, pi_bound , div_pol(grid_idx(1)-1,grid_idx(2)) + delta_div )

                endif

                ! initialize current problem
                v_curr_max(1)   = 0d0
                l_curr_max(1)   = 0d0
                d_curr_max(1)   = 0d0
                div_curr_max(1) = 0d0
                s_curr_max(1)   = 0d0
                a_curr_max(1)   = 0d0

                ! solve problem
                do ll=1,grid_len ! for each loan
                    do oo=1,grid_len   ! for each dividend

                        ! implied securities
                        implied_s =  ( bank_ngrid(grid_idx(1)) - divgrid(oo) - &
                                       ebar*lgrid(ll) - g(lgrid(ll)) )/ebar

                        implied_d = (lgrid(ll) + implied_s )*(1d0-ebar)

                        ! implied wholesale funding
                        implied_a = 0d0

                        ! liquidity requirement naturally satisfied

                        if (implied_d < 0d0) then
                            v_temp = penalty
                        elseif (implied_d > dbar(grid_idx(2))) then   ! liquidity requirement
                            v_temp = penalty
                        elseif (implied_s < 0d0) then     ! negative securities
                            v_temp = penalty
                        else               ! if all constraints satisfied

                            ! initialize RHS of bellman equation
                            RHS = 0d0

                            ! no funding risk
                            do rr=1,size(Rl)        ! for each next period loan shock
                                do ss=1,dbar_size  ! and each deposit funding shock

                                    ! compute networth
                                    net_int = (Rl(rr)-1d0)*lgrid(ll) + i_s*implied_s - (Rd-1d0)*implied_d

                                    stock = lgrid(ll) + implied_s - implied_d

                                    if ( net_int > 0d0) then
                                        net_temp = (1d0-tax)*net_int + stock
                                    else
                                        net_temp = net_int + stock
                                    endif

                                    if (net_temp >0d0) then

                                        ! linear interpolate for value
                                        call linint_Grow( net_temp, bank_na, bank_nb,&
                                                            growth,bank_nlen-1, ixl, ixr, varphi)

                                        v_temp = ( varphi*v(ixl+1,ss) + &
                                                    (1d0-varphi)*v(ixr+1,ss) )

                                        RHS = RHS + prob_l(rr)*dbar_prob(grid_idx(2),ss)*v_temp
                                    endif

                                enddo
                            enddo

                            v_temp = divgrid(oo) + beta*RHS

                        endif

                        ! evaluate value of policy function combination
                        if (v_temp > v_curr_max(1)) then! if beats current, create new argmax and max value

                            v_curr_max(1) = v_temp

                            l_curr_max(1)   = lgrid(ll)
                            s_curr_max(1)   = implied_s
                            d_curr_max(1)   = implied_d
                            a_curr_max(1)   = 0d0
                            div_curr_max(1) = divgrid(oo)

                        endif

                    enddo
                enddo  ! end of policy function loop

                !
                !
                !    Case 2: binding deposit capacity constraint
                !
                !    implication:   I) d = dbar
                !
                !                  II) a = s + l + theta*l**2/2 + div - n - dbar
                !
                !    search over: (l,s,div) --> a
                !
                !
                ! construct grids for control variables
                if ( ( grid_idx(1) == 1 ).or.( iterator == 1 ) ) then  ! if n_b(1) or first iteratior of value function, use default grids

                    call grid_Cons_Equi( lgrid, la,lb)
                    call grid_Cons_Equi( sgrid, sa, sb)
                    call grid_Cons_Equi( divgrid,diva, divb )

                else

                    l_bound = maxval( (/ 0.0000001d0, lpol(grid_idx(1)-1,grid_idx(2)) - back_look*delta_l /) )
                    s_bound = maxval( (/ 0.0000001d0, spol(grid_idx(1)-1,grid_idx(2)) - back_look*delta_s /) )
                    pi_bound = maxval( (/ 0d0, div_pol(grid_idx(1)-1,grid_idx(2)) - delta_div /) )

                    call grid_Cons_Equi( lgrid, l_bound, lpol(grid_idx(1)-1,grid_idx(2)) + delta_l )
                    call grid_Cons_Equi( sgrid, s_bound, spol(grid_idx(1)-1,grid_idx(2)) + delta_s )
                    call grid_Cons_Equi( divgrid, pi_bound , div_pol(grid_idx(1)-1,grid_idx(2)) + delta_div )

                endif

                ! initialize current problem
                v_curr_max(2)   = 0d0
                l_curr_max(2)   = 0d0
                d_curr_max(2)   = 0d0
                div_curr_max(2) = 0d0
                s_curr_max(2)   = 0d0
                a_curr_max(2)   = 0d0

                ! solve problem
                do ll=1,grid_len ! for each loan
                    do oo=1,grid_len   ! for each dividend
                        do nn=1,grid_len  ! for each security

                            ! implies wholesale funding
                            implied_a = sgrid(nn) + lgrid(ll) + g(lgrid(ll)) + &
                                        divgrid(oo) - &
                                        bank_ngrid(grid_idx(1)) - dbar(grid_idx(2))

                            ! implied capital ratio
                            e_temp = (lgrid(ll) + sgrid(nn) - (dbar(grid_idx(2))+implied_a))/&
                                    ( lgrid(ll) + sgrid(nn) )

                            ! implies a liquidity ratio
                            liq_temp = ( (1d0-lr_hair)*sgrid(nn) )/( implied_a ) !- &
                                                !(1d0-lr_hair)*(1d0+hair)/lr_run

                            ! deposit capacity constraint already satisfied

                            if (e_temp > 1d0) then ! capital requirement (above)
                                v_temp = penalty
                            elseif (e_temp < ebar) then ! capital requirement (below)
                                v_temp = penalty
                            elseif (liq_temp < phi_lr) then   ! liquidity requirement
                                v_temp = penalty
                            elseif (implied_a < 0d0) then     ! negative wholesale funding
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
                                            do ss=1,dbar_size  ! and each deposit funding shock

                                                ! compute networth
                                                net_int = (Rl(rr)-1d0)*lgrid(ll) + i_s*(sgrid(nn)-s_tilde) - &
                                                          (Rd-1d0)*dbar(grid_idx(2)) - &
                                                          (Ra-1d0)*(1d0-delta(qq))*implied_a

                                                stock = lgrid(ll) + (sgrid(nn)-s_tilde) - &
                                                        dbar(grid_idx(2)) - (1d0-delta(qq))*implied_a

                                                if ( net_int > 0d0) then
                                                    net_temp = (1d0-tax)*net_int + stock
                                                else
                                                    net_temp = net_int + stock
                                                endif

                                                if (net_temp>0d0) then
                                                    ! linear interpolate for value
                                                    call linint_Grow( net_temp, bank_na, bank_nb,&
                                                                        growth,bank_nlen-1, ixl, ixr, varphi)

                                                    v_temp = ( varphi*v(ixl+1,ss) + &
                                                                (1d0-varphi)*v(ixr+1,ss) )

                                                    RHS = RHS + prob_d(qq)*prob_l(rr)*&
                                                                dbar_prob(grid_idx(2),ss)*v_temp
                                                endif

                                            enddo
                                        enddo
                                    endif

                                enddo

                                v_temp = divgrid(oo) + beta*RHS

                            endif

                            ! evaluate value of policy function combination

                            if (v_temp > v_curr_max(2)) then! if beats current, create new argmax and max value

                                v_curr_max(2) = v_temp

                                l_curr_max(2)   = lgrid(ll)
                                s_curr_max(2)   = sgrid(nn)
                                d_curr_max(2)   = dbar(grid_idx(2))
                                a_curr_max(2)   = implied_a
                                div_curr_max(2) = divgrid(oo)

                            endif

                        enddo
                    enddo
                enddo ! end of policy function loop

                !
                !
                !    Case 3: binding capital & deposit capacity constraint
                !
                !    implication:   I) d = dbar
                !
                !                  II) a = ( l + s )*(1-ebar) - dbar
                !
                !                 III) s = ( n - theta*l^2/2 - div - ebar*l )/ebar
                !
                !
                !    search over: (l,div)   --> s --> a
                !
                !
                ! construct grids for control variables
                if ( ( grid_idx(1) == 1 ).or.( iterator == 1 ) ) then  ! if n_b(1) or first iteratior of value function, use default grids

                    call grid_Cons_Equi( lgrid, la,lb)
                    call grid_Cons_Equi( divgrid,diva, divb )

                else

                    l_bound = maxval( (/ 0.0000001d0, lpol(grid_idx(1)-1,grid_idx(2)) - back_look*delta_l /) )
                    pi_bound = maxval( (/ 0d0, div_pol(grid_idx(1)-1,grid_idx(2)) - delta_div /) )

                    call grid_Cons_Equi( lgrid, l_bound, lpol(grid_idx(1)-1,grid_idx(2)) + delta_l )
                    call grid_Cons_Equi( divgrid, pi_bound , div_pol(grid_idx(1)-1,grid_idx(2)) + delta_div )

                endif

                ! initialize current problem
                v_curr_max(3)   = 0d0
                l_curr_max(3)   = 0d0
                d_curr_max(3)   = 0d0
                div_curr_max(3) = 0d0
                s_curr_max(3)   = 0d0
                a_curr_max(3)   = 0d0

                ! solve problem
                do ll=1,grid_len ! for each loan
                    do oo=1,grid_len   ! for each dividend

                        ! implied securities
                        implied_s =  ( bank_ngrid(grid_idx(1)) - divgrid(oo) - &
                                        ebar*lgrid(ll) - g(lgrid(ll)) )/( ebar )

                        ! implied wholesale funding
                        implied_a = (lgrid(ll) + implied_s )*(1d0-ebar) - dbar(grid_idx(2))

                        ! implied liquidity ratio
                        liq_temp = ( (1d0-lr_hair)*implied_s )/( implied_a )

                        if (liq_temp < phi_lr) then   ! liquidity requirement
                            v_temp = penalty
                        elseif (implied_a < 0d0) then     ! negative wholesale funding
                            v_temp = penalty
                        elseif (implied_s < 0d0) then     ! negative securities
                            v_temp = penalty
                        elseif ( implied_s < (1d0+hair)*implied_a) then  ! collateral constraint
                            v_temp = penalty
                        else               ! if all constraints satisfied

                            ! initialize RHS of bellman equation
                            RHS = 0d0

                            ! for each possible funding shock
                            do qq=1,size(delta)

                                ! check liquidity default
                                if ( pstar*implied_s >= delta(qq)*implied_a) then

                                    s_tilde = ( delta(qq)*implied_a )/pstar

                                    do rr=1,size(Rl)        ! for each next period loan shock
                                        do ss=1,dbar_size  ! and each deposit funding shock

                                            ! compute networth
                                            net_int = (Rl(rr)-1d0)*lgrid(ll) + i_s*(implied_s-s_tilde) - &
                                                      (Rd-1d0)*dbar(grid_idx(2)) - &
                                                      (Ra-1d0)*(1d0-delta(qq))*implied_a

                                            stock = lgrid(ll) + (implied_s-s_tilde) - &
                                                    dbar(grid_idx(2)) - (1d0-delta(qq))*implied_a

                                            if ( net_int > 0d0) then
                                                net_temp = (1d0-tax)*net_int + stock
                                            else
                                                net_temp = net_int + stock
                                            endif

                                            if (net_temp>0d0) then

                                                ! linear interpolate for value
                                                call linint_Grow( net_temp, bank_na, bank_nb,&
                                                                    growth,bank_nlen-1, ixl, ixr, varphi)

                                                v_temp = ( varphi*v(ixl+1,ss) + &
                                                            (1d0-varphi)*v(ixr+1,ss) )

                                                RHS = RHS + prob_d(qq)*prob_l(rr)*&
                                                            dbar_prob(grid_idx(2),ss)*v_temp
                                            endif

                                        enddo
                                    enddo
                                endif

                            enddo

                            v_temp = divgrid(oo) + beta*RHS

                        endif

                        ! evaluate value of policy function combination

                        if (v_temp > v_curr_max(3)) then! if beats current, create new argmax and max value

                            v_curr_max(3) = v_temp

                            l_curr_max(3)   = lgrid(ll)
                            s_curr_max(3)   = implied_s
                            d_curr_max(3)   = dbar(grid_idx(2))
                            a_curr_max(3)   = implied_a
                            div_curr_max(3) = divgrid(oo)

                        endif

                    enddo
                enddo  ! end of policy function loop

                !
                !
                !    Case 4: binding liquidity requirement & deposit capacity constraint
                !
                !    implication:   I) d = dbar
                !
                !                  II) a = (1-h^s)*s/phi_lr
                !
                !                 III) s = ( n + dbar - l - theta*l^2/2 - div )/( 1 - (1-hs)/phi_lr )
                !
                !
                !    search over: (l,div)   --> s --> a
                !
                !

                ! construct grids for control variables
                if ( ( grid_idx(1) == 1 ).or.( iterator == 1 ) ) then  ! if n_b(1) or first iteratior of value function, use default grids

                    call grid_Cons_Equi( lgrid, la,lb)
                    call grid_Cons_Equi( divgrid,diva, divb )

                else

                    l_bound = maxval( (/ 0.0000001d0, lpol(grid_idx(1)-1,grid_idx(2)) - back_look*delta_l /) )
                    pi_bound = maxval( (/ 0d0, div_pol(grid_idx(1)-1,grid_idx(2)) - delta_div /) )

                    call grid_Cons_Equi( lgrid, l_bound, lpol(grid_idx(1)-1,grid_idx(2)) + delta_l )
                    call grid_Cons_Equi( divgrid, pi_bound , div_pol(grid_idx(1)-1,grid_idx(2)) + delta_div )

                endif

                ! initialize current problem
                v_curr_max(4)   = 0d0
                l_curr_max(4)   = 0d0
                d_curr_max(4)   = 0d0
                div_curr_max(4) = 0d0
                s_curr_max(4)   = 0d0
                a_curr_max(4)   = 0d0

                ! solve problem
                do ll=1,grid_len ! for each loan
                    do oo=1,grid_len   ! for each dividend

                        ! implied securities
                        implied_s = (bank_ngrid(grid_idx(1)) + dbar(grid_idx(2)) - lgrid(ll) - g(lgrid(ll)) - divgrid(oo))/&
                                    ( 1d0 - (1d0-lr_hair)/phi_lr )

                        ! implied wholesale funding
                        implied_a = (1d0-lr_hair)*implied_s/phi_lr

                        ! implied capital ratio
                        e_temp = (lgrid(ll) + implied_s - (dbar(grid_idx(2))+implied_a))/&
                                ( lgrid(ll) + implied_s )

                        ! deposit capacity constraint already satisfied
                        if (e_temp > 1d0) then ! capital requirement (above)
                            v_temp = penalty
                        elseif (e_temp < ebar) then ! capital requirement (below)
                            v_temp = penalty
                        elseif (implied_a < 0d0) then     ! negative wholesale funding
                            v_temp = penalty
                        elseif (implied_s < 0d0) then     ! negative securities
                            v_temp = penalty
                        elseif ( implied_s < (1d0+hair)*implied_a) then  ! collateral constraint
                            v_temp = penalty
                        else               ! if all constraints satisfied

                            ! initialize RHS of bellman equation
                            RHS = 0d0

                            ! for each possible funding shock
                            do qq=1,size(delta)

                                ! check liquidity default
                                if ( pstar*implied_s >= delta(qq)*implied_a) then

                                    s_tilde = ( delta(qq)*implied_a )/pstar

                                    do rr=1,size(Rl)        ! for each next period loan shock
                                        do ss=1,dbar_size  ! and each deposit funding shock

                                            ! compute networth
                                            net_int = (Rl(rr)-1d0)*lgrid(ll) + i_s*(implied_s-s_tilde) - &
                                                      (Rd-1d0)*dbar(grid_idx(2)) - &
                                                      (Ra-1d0)*(1d0-delta(qq))*implied_a

                                            stock = lgrid(ll) + (implied_s-s_tilde) - &
                                                    dbar(grid_idx(2)) - (1d0-delta(qq))*implied_a

                                            if ( net_int > 0d0) then
                                                net_temp = (1d0-tax)*net_int + stock
                                            else
                                                net_temp = net_int + stock
                                            endif

                                            if (net_temp>0d0) then

                                                ! linear interpolate for value
                                                call linint_Grow( net_temp, bank_na, bank_nb,&
                                                                    growth,bank_nlen-1, ixl, ixr, varphi)

                                                v_temp = ( varphi*v(ixl+1,ss) + &
                                                            (1d0-varphi)*v(ixr+1,ss) )

                                                RHS = RHS + prob_d(qq)*prob_l(rr)*&
                                                            dbar_prob(grid_idx(2),ss)*v_temp
                                            endif

                                        enddo
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
                            d_curr_max(4)   = dbar(grid_idx(2))
                            a_curr_max(4)   = implied_a
                            div_curr_max(4) = divgrid(oo)

                        endif

                    enddo
                enddo  ! end of policy function loop


                !
                !
                !    Case 5: binding liquidity & capital requirement & deposit capacity constraint
                !
                !    implication:   I) d = dbar
                !
                !                  II) s = ( phi_lr*a )/(1-h^s)
                !
                !                 III) l = ( dbar + a*( 1 - phi_lr*(1-ebar)/(1-hs)) )/(1-ebar)
                !
                !
                !                  IV) div = n + a + dbar - l(a) - theta*l(a)^2/2 - s(a)
                !
                !
                !    search over: (a)   --> s--> l --> div

                ! construct grids for control variables
                if ( ( grid_idx(1) == 1 ).or.( iterator == 1 ) ) then  ! if n_b(1) or first iteratior of value function, use default grids

                    call grid_Cons_Equi( agrid, aa, ab)

                else

                    a_bound = maxval( (/ 0.0000001d0, apol(grid_idx(1)-1,grid_idx(2)) - back_look*delta_a /) )

                    call grid_Cons_Equi( agrid, a_bound, apol(grid_idx(1)-1,grid_idx(2)) + delta_a )

                endif

                ! initialize current problem
                v_curr_max(5)   = 0d0
                l_curr_max(5)   = 0d0
                d_curr_max(5)   = 0d0
                div_curr_max(5) = 0d0
                s_curr_max(5)   = 0d0
                a_curr_max(5)   = 0d0

                ! solve problem
                do oo=1,grid_len   ! for each wholesale funding

                    !implied securities
                    implied_s = ( phi_lr*agrid(oo) )/(1d0-lr_hair)

                    implied_l = (dbar(grid_idx(2)) + agrid(oo)*( 1d0 - phi_lr*(1d0-ebar)/(1d0-lr_hair) ) )/(1d0-ebar)

                    implied_div = ( bank_ngrid(grid_idx(1)) + agrid(oo) + dbar(grid_idx(2)) ) - &
                                    implied_l - g(implied_l) - implied_s

                    if (implied_s < 0d0) then     ! negative securities
                        v_temp = penalty
                    elseif (implied_div < 0d0) then     ! negative securities
                        v_temp = penalty
                    elseif (implied_l < 0d0) then     ! negative cash
                        v_temp = penalty
                    elseif ( implied_s < (1d0+hair)*agrid(oo)) then  ! collateral constraint
                        v_temp = penalty
                    else               ! if all constraints satisfied

                        ! initialize RHS of bellman equation
                        RHS = 0d0

                        ! for each possible funding shock
                        do qq=1,size(delta)

                            ! check liquidity default
                            if ( pstar*implied_s >= delta(qq)*agrid(oo)) then

                                s_tilde = ( delta(qq)*agrid(oo) )/pstar

                                do rr=1,size(Rl)        ! for each next period loan shock
                                    do ss=1,dbar_size  ! and each deposit funding shock

                                        ! compute networth
                                        net_int = (Rl(rr)-1d0)*implied_l + i_s*(implied_s-s_tilde) - &
                                                (Rd-1d0)*dbar(grid_idx(2)) - (Ra-1d0)*(1d0-delta(qq))*agrid(oo)

                                        stock = implied_l + (implied_s-s_tilde) - &
                                                dbar(grid_idx(2)) - (1d0-delta(qq))*agrid(oo)

                                        if ( net_int > 0d0) then
                                            net_temp = (1d0-tax)*net_int + stock
                                        else
                                            net_temp = net_int + stock
                                        endif

                                        if (net_temp>0d0) then

                                            ! linear interpolate for value
                                            call linint_Grow( net_temp, bank_na, bank_nb,&
                                                                growth,bank_nlen-1, ixl, ixr, varphi)

                                            v_temp = ( varphi*v(ixl+1,ss) + &
                                                        (1d0-varphi)*v(ixr+1,ss) )

                                            RHS = RHS + prob_d(qq)*prob_l(rr)*&
                                                        dbar_prob(grid_idx(2),ss)*v_temp
                                        endif


                                    enddo
                                enddo
                            endif

                        enddo

                        v_temp = implied_div + beta*RHS

                    endif

                    ! evaluate value of policy function combination

                    if (v_temp > v_curr_max(5)) then! if beats current, create new argmax and max value

                        v_curr_max(5) = v_temp

                        l_curr_max(5)   = implied_l
                        s_curr_max(5)   = implied_s
                        d_curr_max(5)   = dbar(grid_idx(2))
                        a_curr_max(5)   = agrid(oo)
                        div_curr_max(5) = implied_div

                    endif

                enddo  ! end of policy function loop


                !
                !    Case 6: no assumed binding constraints
                !

                ! construct grids for control variables
                if ( ( grid_idx(1) == 1 ).or.( iterator == 1 ) ) then  ! if n_b(1) or first iteratior of value function, use default grids

                    call grid_Cons_Equi( lgrid, la,lb)

                    if ( da >= dbar(grid_idx(2))) then
                        call grid_Cons_Equi( dgrid, .99d0*dbar(grid_idx(2)), dbar(grid_idx(2)) )
                    else
                        call grid_Cons_Equi( dgrid, da, dbar(grid_idx(2)))
                    endif

                    call grid_Cons_Equi( divgrid,diva, divb )

                else

                    d_bound = maxval( (/ 0.0000001d0, dpol(grid_idx(1)-1,grid_idx(2)) - back_look*delta_d /) )
                    d_ubound = minval( (/ dpol(grid_idx(1)-1,grid_idx(2)) + delta_d, dbar(grid_idx(2)) /) )
                    l_bound = maxval( (/ 0.0000001d0, lpol(grid_idx(1)-1,grid_idx(2)) - back_look*delta_l /) )
                    pi_bound = maxval( (/ 0d0, div_pol(grid_idx(1)-1,grid_idx(2)) - delta_div /) )

                    call grid_Cons_Equi( lgrid, l_bound, lpol(grid_idx(1)-1,grid_idx(2)) + delta_l )

                    if (d_bound >= d_ubound) then
                        call grid_Cons_Equi( dgrid, .99d0*d_bound, d_bound )
                    else
                        call grid_Cons_Equi( dgrid, d_bound, d_ubound )
                    endif

                    call grid_Cons_Equi( divgrid, pi_bound , div_pol(grid_idx(1)-1,grid_idx(2)) + delta_div )


                endif

                ! initialize current problem
                v_curr_max(6)   = 0d0
                l_curr_max(6)   = 0d0
                d_curr_max(6)   = 0d0
                div_curr_max(6) = 0d0
                s_curr_max(6)   = 0d0
                a_curr_max(6)   = 0d0

                ! solve problem
                do ll=1,grid_len ! for each loan
                    do mm=1,grid_len ! for each deposit
                        do oo=1,grid_len   ! for each dividend

                            implied_a = 0d0

                            ! implied securities
                            implied_s = bank_ngrid(grid_idx(1)) + dgrid(mm) - divgrid(oo) -&
                                        lgrid(ll) - g(lgrid(ll))

                            ! implies a capital ratio
                            e_temp = (lgrid(ll) + implied_s - dgrid(mm))/( lgrid(ll) + implied_s  )

                            ! liquidity requirement naturally satisfied

                            if (e_temp > 1d0) then ! capital requirement (above)
                                v_temp = penalty
                            elseif (e_temp < ebar) then ! capital requirement (below)
                                v_temp = penalty
                            elseif (implied_s < 0d0) then
                                v_temp = penalty
                            else               ! if all constraints satisfied

                                ! initialize RHS of bellman equation
                                RHS = 0d0

                                ! no funding risk
                                do rr=1,size(Rl)        ! for each next period loan shock
                                    do ss=1,dbar_size  ! and each deposit funding shock

                                        ! compute networth
                                        net_int = (Rl(rr)-1d0)*lgrid(ll) + i_s*implied_s - (Rd-1d0)*dgrid(mm)

                                        stock = lgrid(ll) + implied_s - dgrid(mm)

                                        if ( net_int > 0d0) then
                                            net_temp = (1d0-tax)*net_int + stock
                                        else
                                            net_temp = net_int + stock
                                        endif

                                        if (net_temp>0d0) then

                                            ! linear interpolate for value
                                            call linint_Grow( net_temp, bank_na, bank_nb,&
                                                                growth,bank_nlen-1, ixl, ixr, varphi)

                                            v_temp = ( varphi*v(ixl+1,ss) + &
                                                        (1d0-varphi)*v(ixr+1,ss) )

                                            RHS = RHS + prob_l(rr)*dbar_prob(grid_idx(2),ss)*v_temp
                                        endif

                                    enddo
                                enddo

                                v_temp = divgrid(oo) + beta*RHS

                            endif

                            ! evaluate value of policy function combination

                            if (v_temp > v_curr_max(6)) then! if beats current, create new argmax and max value

                                v_curr_max(6) = v_temp

                                l_curr_max(6)   = lgrid(ll)
                                s_curr_max(6)   = implied_s
                                d_curr_max(6)   = dgrid(mm)
                                a_curr_max(6)   = 0d0
                                div_curr_max(6) = divgrid(oo)

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
                error = sum( (v-vnew)**2d0 )/(bank_nlen*dbar_size)

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

        real*8, dimension(bank_nlen,dbar_size) :: F0, F1, G_ent
        real*8, parameter :: tolerance = 1e-40
        real*8 :: error_temp
        integer :: ii, jj, kk, ll, mm, nn, qq
        integer, parameter :: it_max = 1000
        real*8 :: net_temp
        integer, dimension(1) :: min_idx
        real*8 :: s_tilde, mstar, excess_cash, mstar_agg
        real*8 :: net_int, stock

        integer :: il, ir
        real*8  :: phi

!        default_liq_prob = 0d0
        ! compute default probabilities
        do ii=1,bank_nlen
            do jj=1,dbar_size

                if ( bank_ngrid(ii) <= nstar(jj)) then
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
                                ! for each capacity constraint shock
                                do nn=1,dbar_size

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
                                    if (net_temp <= nstar(jj) ) then
                                        default_prob(ii,jj) = default_prob(ii,jj) + &
                                                        prob_l(mm)*dbar_prob(jj,nn)*prob_d(ll)
                                    endif

                                enddo
                            enddo


                        endif

                    enddo
                endif
            enddo
        enddo

        !Initialize problem
        F0 = 1d0/(bank_nlen*dbar_size)
        error_temp = tolerance + 1d0

        do qq = 1,it_max        !for each iteration in updating process

            F1 = 0d0

            do ii=1,bank_nlen
                do jj=1,dbar_size

                    ! for each delta shock
                    do ll=1,size(delta)

                        ! if no default, compute liquidations
                        if ( ( pstar*spol(ii,jj)) >= delta(ll)*apol(ii,jj) ) then

                            s_tilde = ( delta(ll)*apol(ii,jj) )/pstar

                            ! for each return shock
                            do mm=1,size(Rl)
                                ! for each capacity constraint
                                do nn=1,dbar_size


                                    net_int = (Rl(mm)-1d0)*lpol(ii,jj) + i_s*(spol(ii,jj)-s_tilde) &
                                                     - (Rd-1d0)*dpol(ii,jj) - (Ra-1d0)*(1d0-delta(ll))*apol(ii,jj)

                                    stock = lpol(ii,jj) + (spol(ii,jj)-s_tilde) -&
                                            dpol(ii,jj) - (1d0-delta(ll))*apol(ii,jj)

                                    if ( net_int > 0d0 ) then
                                        net_temp = (1d0-tax)*net_int + stock
                                    else
                                        net_temp = net_int + stock
                                    endif

                                    ! if no default, record new (networth,theta,capacity constraint)
                                    if (net_temp > nstar(jj) ) then

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

                                        F1(il,nn) = F1(il,nn) + phi*prob_l(mm)*dbar_prob(jj,nn)*prob_d(ll)*&
                                                                                        F0(ii,jj)
                                        F1(ir,nn) = F1(ir,nn) + (1d0-phi)*prob_l(mm)*dbar_prob(jj,nn)*&
                                                                                        prob_d(ll)*F0(ii,jj)
                                    endif
                                enddo
                            enddo



                        endif


                    enddo

                enddo
            enddo

            ! create entrant distribution (with replacement)
            do ii = 1,bank_nlen
                do jj=1,dbar_size
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

            error_temp = sum(abs(F1-F0))/(bank_nlen*dbar_size )

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
        real*8, dimension(bank_nlen,dbar_size) :: bshare
        real*8 :: agg_a
        real*8 :: obj1, obj2, obj3, obj4
        real*8 :: temp_net
        real*8 :: excess_cash
        real*8 :: net_int, stock

        agg_a = 0d0
        ! compute aggregate wholesale lending
        do ii=1,bank_nlen
            do jj=1,dbar_size
                agg_a = agg_a + F_stationary(ii,jj)*apol(ii,jj)
            enddo
        enddo

        ! compute bank wholesalefunding shares
        do ii=1,bank_nlen
            do jj=1,dbar_size
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
            do jj=1,dbar_size

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
                            do nn=1,dbar_size

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
                                if ( temp_net < nstar(jj)  ) then  ! insolvency default
                                    obj3 = obj3 + bshare(ii,jj)*prob_d(ll)*(1d0-delta(ll))*prob_l(mm)*dbar_prob(jj,nn)
                                else
                                    obj4 = obj4 + bshare(ii,jj)*prob_d(ll)*(1d0-delta(ll))*prob_l(mm)*dbar_prob(jj,nn)
                                endif


                            enddo
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
            do jj=1,dbar_size            ! each deposit capacity constraint

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

        real*8, dimension(bank_nlen,dbar_size) :: transform_f

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
            do jj=1,dbar_size           ! for each deposit capacity constraint

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
        real*8, dimension(8,bank_nlen*dbar_size) :: corr_objs  ! distribution weight, size, risk-weighted equity, leverage, liquidity, ins default, liq default, equity value
        real*8, dimension(9,bank_nlen*dbar_size) :: corr_objs_intra
        real*8, dimension(2, bank_nlen*dbar_size ) :: roe_stack
        real*8, dimension(bank_nlen,dbar_size) :: roe_pol1, roe_pol2, lev_pol

        integer :: ii, jj, kk, ll, mm, iterator

        if (my_id==0) then

            !~~~~~~~~~~~!
            !           !
            !   Plots   !
            !           !
            !~~~~~~~~~~~!
            ! bank policy functions
            call plot(bank_ngrid,lpol(:,3),legend='loan',marker=2)
            call plot(bank_ngrid,spol(:,3),legend='securities')
            call plot(bank_ngrid,dpol(:,3),legend='deposit')
            call plot(bank_ngrid,apol(:,3),legend='wholesale')
            call plot(bank_ngrid,div_pol(:,3),legend='dividends',marker=2)
            call execplot(title='Bank Policy Functions',filename='policies',filetype='png',output='policies')

            ! deposit policy functions with different dbar
            call plot(bank_ngrid,dpol(:,1),legend='dbar = 1',marker=2)
            call plot(bank_ngrid,dpol(:,3),legend='dbar = 3',marker=2)
            call plot(bank_ngrid,dpol(:,6),legend='dbar = 6',marker=2)
            call execplot(title='Deposit Policy Functions',filename='deposit',filetype='png',output='deposit')

            ! loan policy functions with different dbar
            call plot(bank_ngrid,lpol(:,1),legend='dbar = 1',marker=2)
            call plot(bank_ngrid,lpol(:,3),legend='dbar = 3',marker=2)
            call plot(bank_ngrid,lpol(:,6),legend='dbar = 6',marker=2)
            call execplot(title='Loan Policy Functions',filename='loan',filetype='png',output='loan')

            ! wholesale policy functions with different dbar
            call plot(bank_ngrid,apol(:,1),legend='dbar = 1',marker=2)
            call plot(bank_ngrid,apol(:,3),legend='dbar = 3',marker=2)
            call plot(bank_ngrid,apol(:,6),legend='dbar = 6',marker=2)
            call execplot(title='Wholesale Policy Functions',filename='wholesale',filetype='png',output='wholesale')

            ! dividend policy functions with different dbar
            call plot(bank_ngrid,div_pol(:,1),legend='dbar = 1',marker=2)
            call plot(bank_ngrid,div_pol(:,3),legend='dbar = 3',marker=2)
            call plot(bank_ngrid,div_pol(:,6),legend='dbar = 6',marker=2)
            call execplot(title='Dividend Policy Functions',filename='dividend',filetype='png',output='dividend')

            ! stationary distribution
            call plot_hist(bank_ngrid,F_stationary(:,1)        ,legend='dbar=1')
            call plot_hist(bank_ngrid,F_stationary(:,6)        ,legend='dbar=6')
            call execplot(title='Stationary Distribution of Networth' ,filename='dist',filetype='png',output='dist')

            ! leverage equity ratio policy function
            do ii=1,bank_nlen
                do jj=1,dbar_size

                    lev_pol(ii,jj) = ( lpol(ii,jj) + spol(ii,jj) - dpol(ii,jj) - apol(ii,jj) )/(lpol(ii,jj) + spol(ii,jj))

                enddo
            enddo

            call plot(bank_ngrid,lev_pol(:,1),legend='dbar = 1',marker=2)
            call plot(bank_ngrid,lev_pol(:,3),legend='dbar = 3',marker=2)
            call plot(bank_ngrid,lev_pol(:,6),legend='dbar = 6',marker=2)
            call execplot(title='Leverage Ratio Policy Function',filename='lev',filetype='png',output='lev')

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
                do jj=1,dbar_size


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
                    do kk=1,dbar_size


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
                do jj=1,dbar_size


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
                do jj=1,dbar_size

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
                do jj=1,dbar_size

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
                do jj=1,dbar_size

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
                do kk=1,dbar_size

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

            do jj=1,dbar_size
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
