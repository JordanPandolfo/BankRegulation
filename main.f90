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

    dbar_mu(1) = 0.0122553249d0
    dbar_mu(2) = 0.0335536442d0
    dbar_mu(3) = 0.0106243001d0
    gam        = 0.961326537d0
    pstar      = 0.971990977d0
    l_sigma2   = 0.000415561944d0
    i_s        = .00390620392d0

    beta = gam*beta_pre

    alpha = 0.00021261271199274906
    omega_outside = (no_liq_price/alpha)**(1d0/(alpha-1d0))

    bank_nb   = (/ .005d0, .005d0, .005d0 /)!(/ dbar_mu(1)*1.45188876d0*( (l_mu + 2.5d0*l_sigma2**(0.5d0)) - Rd*(1d0-.05d0)  ), &
                !   dbar_mu(2)*1.41676404d0*( (l_mu + 2.5d0*l_sigma2**(0.5d0)) - Rd*(1d0-.05d0)  ), &
                !   dbar_mu(3)*1.22691827d0*( (l_mu + 2.5d0*l_sigma2**(0.5d0)) - Rd*(1d0-.05d0)  ) &
                !/)

    delta_l   = 20d0*( bank_nb(1) )/bank_nlen ! .002d0 !.01d0

    delta_div = delta_l*.1d0
    delta_a   = delta_l*1d0
    delta_d   = delta_l*1d0
    delta_s   = delta_l*1d0
    delta_c   = delta_l*.25d0

    ! MPI BARRIER
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
    !                                                  !
    !   Initialize Objects for Transitional Dynamics   !
    !                                                  !
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!

    ! set shock process
    dshock = 0d0
    ishock = 0d0

    dshock(2) = 0.05d0
    !dshock(3) = .2d0
    !dshock(4) = .1d0

    !ishock(2) = .01d0
    !ishock(3) = .0075d0
    !ishock(4) = .005d0
    !ishock(4) = .0025d0

    ! intialize equilibrium objects
    pguess_prime = 0d0
    rguess_prime = 0d0
    vtran        = 0d0
    ltran        = 0d0
    stran        = 0d0
    ctran        = 0d0
    dtran        = 0d0
    div_tran     = 0d0
    atran        = 0d0
    Ftran        = 0d0
    indeftran    = 0d0
    liqdeftran   = 0d0

    tran_error_temp = 10000d0
    iterator = 1

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
    allocate( cgrid(grid_len) )
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

    ! set guess for sequence of prices (just set it to stationary eqm values)
    pguess =  pstar

    rguess =  Ra

    ! populate t=1 and t=T transition objects
    pguess_prime(1)       = pstar
    rguess_prime(1)       = Ra
    vtran(1,:,:,:)        = v
    ltran(1,:,:,:)        = lpol
    stran(1,:,:,:)        = spol
    ctran(1,:,:,:)        = cpol
    dtran(1,:,:,:)        = dpol
    div_tran(1,:,:,:)     = div_pol
    atran(1,:,:,:)        = apol
    Ftran(1,:,:,:)        = F_stationary
    indeftran(1,:,:,:)    = default_prob
    liqdeftran(1,:,:,:)   = default_liq_prob

    pguess_prime(T)       = pstar
    rguess_prime(T)       = Ra
    vtran(T,:,:,:)        = v
    ltran(T,:,:,:)        = lpol
    stran(T,:,:,:)        = spol
    ctran(T,:,:,:)        = cpol
    dtran(T,:,:,:)        = dpol
    div_tran(T,:,:,:)     = div_pol
    atran(T,:,:,:)        = apol
    Ftran(T,:,:,:)        = F_stationary
    indeftran(T,:,:,:)    = default_prob
    liqdeftran(T,:,:,:)   = default_liq_prob

    ! MPI BARRIER
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)

    !~~~~~~~~~~~~~~~~~~~~~~~!
    !                       !
    !   Solve Transition    !
    !                       !
    !~~~~~~~~~~~~~~~~~~~~~~~!
    deallocate( lgrid )
    deallocate( divgrid )
    deallocate( cgrid )
    deallocate( sgrid )
    deallocate( dgrid )
    deallocate( agrid )

    grid_len = 5
    allocate( lgrid(grid_len) )
    allocate( divgrid(grid_len) )
    allocate( cgrid(grid_len) )
    allocate( sgrid(grid_len) )
    allocate( dgrid(grid_len) )
    allocate( agrid(grid_len) )


    ! while error criterion is not met
    !do while ( (tran_error_temp > .1d0).and.(iterator<20) )

    ! record prices
    if (my_id==0) then
        open(19,file='ra_guess.csv',access = 'append')
        write(19,*) rguess
        flush(19)
        close(19)

        open(20,file='p_guess.csv',access = 'append')
        write(20,*) pguess
        flush(20)
        close(20)
    endif

    ! for each time period t = T-1 to 2
    do tt = T-1,2,-1

        if (my_id==0) then
            write(*,*) 'transition period:',tt
        endif

        ! solve value function, policy functions at time t
            ! receive value function v_{tt+1}, agg shocks s_{tt+1} and prices p_{tt} --> ( v_{tt}, y_{tt} )

        call bank_update( vtran(tt+1,:,:,:), (/ dshock(tt+1),ishock(tt+1) /), (/ pguess(tt), rguess(tt) /), &
                          vtran(tt,:,:,:), ltran(tt,:,:,:)   , stran(tt,:,:,:), ctran(tt,:,:,:), &
                          dtran(tt,:,:,:), div_tran(tt,:,:,:), atran(tt,:,:,:) )

!        vtran(tt,:,:,:) = vtran(1,:,:,:)
!        ltran(tt,:,:,:) = ltran(1,:,:,:)
!        stran(tt,:,:,:) = stran(1,:,:,:)
!        ctran(tt,:,:,:) = ctran(1,:,:,:)
!        dtran(tt,:,:,:) = dtran(1,:,:,:)
!        div_tran(tt,:,:,:) = div_tran(1,:,:,:)
!        atran(tt,:,:,:) = atran(1,:,:,:)

        ! MPI BARRIER
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)

    enddo

    ! compute joint distribution for periods 2,...,T-1
    do tt=1,T-2

        if (my_id==0) then
            write(*,*) 'distribution period:',tt
        endif

        ! generate joint distribution
        call get_dist_tran(Ftran(tt,:,:,:),ltran(tt,:,:,:),stran(tt,:,:,:),ctran(tt,:,:,:),dtran(tt,:,:,:),&
                                div_tran(tt,:,:,:),atran(tt,:,:,:),&
                             (/ dshock(tt+1),ishock(tt+1) /), (/ pguess(tt), rguess(tt) /), &
                             Ftran(tt+1,:,:,:)  )

        !Ftran(tt+1,:,:,:) = Ftran(1,:,:,:)


        ! MPI BARRIER
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)

    enddo

        if (my_id==0) then

            call plot(bank_ngrid(1,:),Ftran(1,:,3,3),legend='t=1')
            call plot(bank_ngrid(1,:),Ftran(2,:,3,3),legend='t=2')
            call plot(bank_ngrid(1,:),Ftran(3,:,3,3),legend='t=3')
            call plot(bank_ngrid(1,:),Ftran(4,:,3,3),legend='t=4')
            call plot(bank_ngrid(1,:),Ftran(5,:,3,3),legend='t=5')
            call plot(bank_ngrid(1,:),Ftran(6,:,3,3),legend='t=6')

            call execplot(title='dist',filename='dist_beg',filetype='png',output='dist_beg')

            call plot(bank_ngrid(1,:),Ftran(T-3,:,3,3),legend='t=T-3')
            call plot(bank_ngrid(1,:),Ftran(T-2,:,3,3),legend='t=T-2')
            call plot(bank_ngrid(1,:),Ftran(T-1,:,3,3),legend='t=T-1')
            call plot(bank_ngrid(1,:),Ftran(T,:,3,3),legend='t=T')

            call execplot(title='dist',filename='dist_end',filetype='png',output='dist_end')

            call plot(bank_ngrid(1,:),ltran(1,:,3,3),legend='t=1')
            call plot(bank_ngrid(1,:),ltran(2,:,3,3),legend='t=2')
            call plot(bank_ngrid(1,:),ltran(3,:,3,3),legend='t=3')
            call plot(bank_ngrid(1,:),ltran(T-2,:,3,3),legend='t=T-2')
            call plot(bank_ngrid(1,:),ltran(T-1,:,3,3),legend='t=T-1')
            call plot(bank_ngrid(1,:),ltran(T,:,3,3),legend='t=T')
            call execplot(title='loans',filename='all_loan',filetype='png',output='all_loan')


            call plot(bank_ngrid(1,:),ltran(T,:,3,3),legend='loan')
            call plot(bank_ngrid(1,:),stran(T,:,3,3),legend='securities')
            call plot(bank_ngrid(1,:),ctran(T,:,3,3),legend='cash')
            call plot(bank_ngrid(1,:),dtran(T,:,3,3),legend='deposit')
            call plot(bank_ngrid(1,:),atran(T,:,3,3),legend='wholesale')
            call plot(bank_ngrid(1,:),div_tran(T,:,3,3),legend='div')
            call execplot(title='last_loans',filename='last_per',filetype='png',output='last_per')

            call plot(bank_ngrid(1,:),ltran(T-1,:,3,3),legend='loan')
            call plot(bank_ngrid(1,:),stran(T-1,:,3,3),legend='securities')
            call plot(bank_ngrid(1,:),ctran(T-1,:,3,3),legend='cash')
            call plot(bank_ngrid(1,:),dtran(T-1,:,3,3),legend='deposit')
            call plot(bank_ngrid(1,:),atran(T-1,:,3,3),legend='wholesale')
            call plot(bank_ngrid(1,:),div_tran(T-1,:,3,3),legend='div')
            call execplot(title='second_last_loans',filename='second_last_per',filetype='png',output='second_last_per')

        endif

    ! compute implied prices for periods 2,...,T-1
    do tt=2,T-1

        ! compute implied rate Ra and liquidation price p^*
        rguess_prime(tt) = implied_ra_rate_tran(Ftran(tt,:,:,:),ltran(tt,:,:,:),stran(tt,:,:,:),ctran(tt,:,:,:),&
                                                                dtran(tt,:,:,:),div_tran(tt,:,:,:), atran(tt,:,:,:), &
                                                            (/ dshock(tt+1),ishock(tt+1) /),(/ pguess(tt), rguess(tt) /) )

        Agg_Liq = Aggregate_Liq_tran(Ftran(tt,:,:,:),stran(tt,:,:,:),ctran(tt,:,:,:),atran(tt,:,:,:),&
        dshock(tt+1), pguess(tt) )

        pguess_prime(tt) = alpha*(Agg_Liq(2)+omega_outside)**(alpha-1d0)


    enddo

    ! compute price error (avg % deviation from guess, in basis points)
    tran_error_temp = 100d0*100d0*( sum(abs(rguess - rguess_prime)/rguess_prime)/T + &
                                    sum(abs(pguess - pguess_prime)/pguess_prime)/T      )/2d0


    if (my_id==0) then
        open(17,file='price_iter.csv',access = 'append')
        write(17,*) iterator
        flush(17)
        close(17)

        open(18,file='price_error.csv',access = 'append')
        write(18,*) tran_error_temp
        flush(18)
        close(18)
    endif

    if (my_id==0) then
        write(*,*)
        write(*,*) 'PRICE GUESS ITERATION:', iterator
        write(*,*)
        write(*,*) 'Avg Error (Basis Points):', tran_error_temp
        write(*,*)
    endif

    ! update price guess with midpoint rule
    rguess = .85d0*rguess + (1d0-.85d0)*rguess_prime
    pguess = .85d0*pguess + (1d0-.85d0)*pguess_prime

    iterator = iterator + 1

    !enddo

    ! ouput
    call output_tran()

    ! MPI BARRIER
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)


    if (my_id==0) then
        call toc()
    endif

    ! MPI BARRIER
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)

contains

    subroutine bank_update( vp, shock, price, vrec, lrec, srec, crec, drec, div_rec, arec )

        implicit none

        ! input/output variables
        real*8, dimension(bank_nlen,size(theta),dbar_size), intent(in)  :: vp
        real*8, dimension(2), intent(in)                                :: shock, price
        real*8, dimension(bank_nlen,size(theta),dbar_size), intent(inout) :: vrec, lrec, srec, crec, drec, div_rec, arec

        ! Solving bank problem
        integer :: ii, jj, kk, ll, mm, nn, oo, pp, qq,rr,ss
        real*8 :: v_temp, RHS, net_temp, e_temp
        real*8, dimension(6) :: v_curr_max, l_curr_max, d_curr_max, div_curr_max, c_curr_max, s_curr_max, a_curr_max
        real*8 :: a_bound, l_bound, s_bound, d_bound,c_bound, liq_temp, implied_a, &
                  implied_s, implied_c, d_ubound, implied_d, implied_l, implied_div
        real*8 :: c_tilde, s_tilde, excess_cash  ! liquidation values
        integer, dimension(1) :: v_idx
        real*8 :: net_int, stock, pi_bound
        real*8 :: temp_pen
        real*8            :: error


        ! linear interpolation variables
        integer :: ixl, ixr
        real*8  :: varphi

        ! spline interpolation variables
        real*8, dimension( size(theta), dbar_size, bank_nlen+2 ) :: coeff

        !---------------!
        !   mpi stuff   !
        !---------------!
        real*8, allocatable, dimension(:) :: ltemp, dtemp, stemp, atemp, divtemp, vtemp, ctemp
        real*8, allocatable, dimension(:) :: ltemp_all, dtemp_all, stemp_all, atemp_all, divtemp_all, vtemp_all, ctemp_all
        integer, dimension(3) :: grid_idx   ! 3 states (nb,theta,dbar)
        integer, allocatable, dimension(:,:,:) :: idx_full_grid

        integer :: total_grid, proc_len

        ! total size of state space grid
        total_grid = bank_nlen*size(theta)*dbar_size

        if ( mod(total_grid,nprocs) == 0) then  ! perfect fit
            proc_len = total_grid/nprocs
        else
            proc_len = total_grid/nprocs + 1     ! the last proccessor will have some extra/empty slots
        endif

        ! construct temporary policy functions
        allocate( ltemp(proc_len) )
        allocate( dtemp(proc_len) )
        allocate( stemp(proc_len) )
        allocate( ctemp(proc_len) )
        allocate( atemp(proc_len) )
        allocate( divtemp(proc_len) )
        allocate( vtemp(proc_len) )

        allocate( ltemp_all(proc_len*nprocs) )
        allocate( dtemp_all(proc_len*nprocs) )
        allocate( stemp_all(proc_len*nprocs) )
        allocate( ctemp_all(proc_len*nprocs) )
        allocate( atemp_all(proc_len*nprocs) )
        allocate( divtemp_all(proc_len*nprocs) )
        allocate( vtemp_all(proc_len*nprocs) )

        allocate( idx_full_grid(nprocs,proc_len,3) )

        call grid_f(nprocs,proc_len,idx_full_grid)

        ! MPI BARRIER
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)

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
            !                 III) a = 0, c = 0
            !
            !    search over: (l,div)
            !
            !

            ! construct grids for control variables
            if ( ( grid_idx(1) == 1 ) ) then  ! if n_b(1) or first iteratior of value function, use default grids

                call grid_Cons_Equi( lgrid, la(grid_idx(2)),lb(grid_idx(2))   )
                call grid_Cons_Equi( divgrid,diva(grid_idx(2)), divb(grid_idx(2)) )

            else

                l_bound = maxval( (/ 0.0000001d0, lpol(grid_idx(1)-1,grid_idx(2),grid_idx(3)) - back_look*delta_l /) )
                pi_bound = maxval( (/ 0d0, div_pol(grid_idx(1)-1,grid_idx(2),grid_idx(3)) - delta_div /) )

                call grid_Cons_Equi( lgrid, l_bound, lpol(grid_idx(1)-1,grid_idx(2),grid_idx(3)) + delta_l )
                call grid_Cons_Equi( divgrid, pi_bound , div_pol(grid_idx(1)-1,grid_idx(2),grid_idx(3)) + delta_div )

            endif

            ! initialize current problem
            v_curr_max(1)   = 0d0
            l_curr_max(1)   = 0d0
            d_curr_max(1)   = 0d0
            div_curr_max(1) = 0d0
            s_curr_max(1)   = 0d0
            c_curr_max(1)   = 0d0
            a_curr_max(1)   = 0d0

            ! solve problem
            do ll=1,grid_len ! for each loan
                do oo=1,grid_len   ! for each dividend

                    ! implied securities
                    implied_s =  ( bank_ngrid(grid_idx(2),grid_idx(1)) - divgrid(oo) - &
                                   ebar*lgrid(ll) - g(lgrid(ll),theta(grid_idx(2))) )/ebar

                    implied_d = (lgrid(ll) + implied_s )*(1d0-ebar)

                    ! implied wholesale funding
                    implied_a = 0d0
                    implied_c = 0d0

                    ! liquidity requirement naturally satisfied

                    if (implied_d < 0d0) then
                        v_temp = penalty
                    elseif (implied_d > dbar(grid_idx(2),grid_idx(3))) then   ! liquidity requirement
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
                                net_int = (Rl(rr)-1d0-shock(2))*lgrid(ll) + i_s*implied_s - (Rd-1d0)*implied_d

                                stock = lgrid(ll) + implied_s - implied_d

                                if ( net_int > 0d0) then
                                    net_temp = (1d0-tax)*net_int + stock
                                else
                                    net_temp = net_int + stock
                                endif

                                if (net_temp >0d0) then

                                    ! linear interpolate for value
                                    call linint_Grow( net_temp, bank_na(grid_idx(2)), bank_nb(grid_idx(2)),&
                                                        growth,bank_nlen-1, ixl, ixr, varphi)

                                    v_temp = ( varphi*vp(ixl+1,grid_idx(2),ss) + &
                                                (1d0-varphi)*vp(ixr+1,grid_idx(2),ss) )

                                    RHS = RHS + prob_l(rr)*dbar_prob(grid_idx(2),grid_idx(3),ss)*v_temp
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
                        c_curr_max(1)   = 0d0
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
            !                  II) a = c + s + l + theta*l**2/2 + div - n - dbar
            !
            !    search over: (l,c,s,div) --> a
            !
            !
            ! construct grids for control variables
            if ( ( grid_idx(1) == 1 ) ) then  ! if n_b(1) or first iteratior of value function, use default grids

                call grid_Cons_Equi( lgrid, la(grid_idx(2)),lb(grid_idx(2)))
                call grid_Cons_Equi( sgrid, sa(grid_idx(2)), sb(grid_idx(2)))
                call grid_Cons_Equi( cgrid, ca(grid_idx(2)), cb(grid_idx(2)))
                call grid_Cons_Equi( divgrid,diva(grid_idx(2)), divb(grid_idx(2)) )

            else

                l_bound = maxval( (/ 0.0000001d0, lpol(grid_idx(1)-1,grid_idx(2),grid_idx(3)) - back_look*delta_l /) )
                s_bound = maxval( (/ 0.0000001d0, spol(grid_idx(1)-1,grid_idx(2),grid_idx(3)) - back_look*delta_s /) )
                c_bound = maxval( (/ 0.0000001d0, cpol(grid_idx(1)-1,grid_idx(2),grid_idx(3)) - back_look*delta_c /) )
                pi_bound = maxval( (/ 0d0, div_pol(grid_idx(1)-1,grid_idx(2),grid_idx(3)) - delta_div /) )

                call grid_Cons_Equi( lgrid, l_bound, lpol(grid_idx(1)-1,grid_idx(2),grid_idx(3)) + delta_l )
                call grid_Cons_Equi( sgrid, s_bound, spol(grid_idx(1)-1,grid_idx(2),grid_idx(3)) + delta_s )
                call grid_Cons_Equi( cgrid, c_bound, cpol(grid_idx(1)-1,grid_idx(2),grid_idx(3)) + delta_c )
                call grid_Cons_Equi( divgrid, pi_bound , div_pol(grid_idx(1)-1,grid_idx(2),grid_idx(3)) + delta_div )

            endif

            ! initialize current problem
            v_curr_max(2)   = 0d0
            l_curr_max(2)   = 0d0
            d_curr_max(2)   = 0d0
            div_curr_max(2) = 0d0
            s_curr_max(2)   = 0d0
            c_curr_max(2)   = 0d0
            a_curr_max(2)   = 0d0

            ! solve problem
            do ll=1,grid_len ! for each loan
                do mm=1,grid_len ! for each cash
                    do nn=1,grid_len  ! for each security
                        do oo=1,grid_len   ! for each dividend

                            ! implies wholesale funding
                            implied_a = cgrid(mm) + sgrid(nn) + lgrid(ll) + g(lgrid(ll),theta(grid_idx(2))) + &
                                        divgrid(oo) - &
                                        bank_ngrid(grid_idx(2),grid_idx(1)) - dbar(grid_idx(2),grid_idx(3))

                            ! implied capital ratio
                            e_temp = (lgrid(ll) + cgrid(mm) + sgrid(nn) - (dbar(grid_idx(2),grid_idx(3))+implied_a))/&
                                    ( lgrid(ll) + sgrid(nn) + cgrid(mm) )

                            ! implies a liquidity ratio
                            liq_temp = ( cgrid(mm) + (1d0-lr_hair)*sgrid(nn) )/( lr_run*implied_a ) !- &
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

                                    ! compute cash liquidations
                                    excess_cash = cgrid(mm)

                                    ! check liquidity default
                                    if ( excess_cash + price(1)*sgrid(nn) >= (delta(qq)+shock(1))*implied_a) then

                                        c_tilde =  minval( (/ (delta(qq)+shock(1))*implied_a, excess_cash /) )

                                        if ( c_tilde == excess_cash ) then  ! if not enough cash to settle
                                            ! have to liquidate securities
                                            s_tilde = ( (delta(qq)+shock(1))*implied_a - excess_cash )/price(1)
                                        else
                                            ! if enough cash, security liquidations = 0
                                            s_tilde = 0d0
                                        endif

                                        do rr=1,size(Rl)        ! for each next period loan shock
                                            do ss=1,dbar_size  ! and each deposit funding shock

                                                ! compute networth
                                                net_int = (Rl(rr)-shock(2)-1d0)*lgrid(ll) + i_s*(sgrid(nn)-s_tilde) - &
                                                          (Rd-1d0)*dbar(grid_idx(2),grid_idx(3)) - &
                                                          (price(2)-1d0)*(1d0-delta(qq)-shock(1))*implied_a

                                                stock = lgrid(ll) + (sgrid(nn)-s_tilde) + (cgrid(mm)-c_tilde) - &
                                                        dbar(grid_idx(2),grid_idx(3)) - (1d0-delta(qq)-shock(1))*implied_a

                                                if ( net_int > 0d0) then
                                                    net_temp = (1d0-tax)*net_int + stock
                                                else
                                                    net_temp = net_int + stock
                                                endif

                                                if (net_temp>0d0) then
                                                    ! linear interpolate for value
                                                    call linint_Grow( net_temp, bank_na(grid_idx(2)), bank_nb(grid_idx(2)),&
                                                                        growth,bank_nlen-1, ixl, ixr, varphi)

                                                    v_temp = ( varphi*vp(ixl+1,grid_idx(2),ss) + &
                                                                (1d0-varphi)*vp(ixr+1,grid_idx(2),ss) )

                                                    RHS = RHS + prob_d(qq)*prob_l(rr)*&
                                                                dbar_prob(grid_idx(2),grid_idx(3),ss)*v_temp
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
                                c_curr_max(2)   = cgrid(mm)
                                d_curr_max(2)   = dbar(grid_idx(2),grid_idx(3))
                                a_curr_max(2)   = implied_a
                                div_curr_max(2) = divgrid(oo)

                            endif


                        enddo
                    enddo
                enddo
            enddo ! end of policy function loop

            !
            !
            !    Case 3: binding capital & deposit capacity constraint
            !
            !    implication:   I) d = dbar
            !
            !                  II) a = ( l + s + c )*(1-ebar) - dbar
            !
            !                 III) c = ( n - theta*l^2/2 - div - ebar*s - ebar*l )/ebar
            !
            !
            !    search over: (l,s,div)   --> c --> a
            !
            !
            ! construct grids for control variables
            if ( ( grid_idx(1) == 1 ) ) then  ! if n_b(1) or first iteratior of value function, use default grids

                call grid_Cons_Equi( lgrid, la(grid_idx(2)),lb(grid_idx(2)))
                call grid_Cons_Equi( divgrid,diva(grid_idx(2)), divb(grid_idx(2)) )
                call grid_Cons_Equi( sgrid, sa(grid_idx(2)), sb(grid_idx(2)))

            else

                s_bound = maxval( (/ 0.0000001d0, spol(grid_idx(1)-1,grid_idx(2),grid_idx(3)) - back_look*delta_s /) )
                l_bound = maxval( (/ 0.0000001d0, lpol(grid_idx(1)-1,grid_idx(2),grid_idx(3)) - back_look*delta_l /) )
                pi_bound = maxval( (/ 0d0, div_pol(grid_idx(1)-1,grid_idx(2),grid_idx(3)) - delta_div /) )

                call grid_Cons_Equi( lgrid, l_bound, lpol(grid_idx(1)-1,grid_idx(2),grid_idx(3)) + delta_l )
                call grid_Cons_Equi( sgrid, s_bound, spol(grid_idx(1)-1,grid_idx(2),grid_idx(3)) + delta_s )
                call grid_Cons_Equi( divgrid, pi_bound , div_pol(grid_idx(1)-1,grid_idx(2),grid_idx(3)) + delta_div )

            endif

            ! initialize current problem
            v_curr_max(3)   = 0d0
            l_curr_max(3)   = 0d0
            d_curr_max(3)   = 0d0
            div_curr_max(3) = 0d0
            s_curr_max(3)   = 0d0
            c_curr_max(3)   = 0d0
            a_curr_max(3)   = 0d0

            ! solve problem
            do ll=1,grid_len ! for each loan
                do mm=1,grid_len ! for each security
                    do oo=1,grid_len   ! for each dividend

                        ! implied securities [ n - Phi(div) - ebar*wl*l - theta*l^2/2 - ebar*wc*p_c*c ]/(  ebar*ws*p_s )
                        implied_c =  ( bank_ngrid(grid_idx(2),grid_idx(1)) - divgrid(oo) - &
                                        ebar*lgrid(ll) - g(lgrid(ll),theta(grid_idx(2))) - &
                                        ebar*sgrid(mm) )/( ebar )

                        ! implied wholesale funding a = l*(1-ebar*wl) + p_s*s*(1-ebar*ws) + p_c*c*(1-ebar*wc) - d
                        implied_a = (lgrid(ll) + sgrid(mm) + implied_c)*(1d0-ebar) - dbar(grid_idx(2),grid_idx(3))

                        ! implied liquidity ratio
                        liq_temp = ( implied_c + (1d0-lr_hair)*sgrid(mm) )/( lr_run*implied_a ) !- &
                                                !(1d0-lr_hair)*(1d0+hair)/lr_run

                        if (liq_temp < phi_lr) then   ! liquidity requirement
                            v_temp = penalty
                        elseif (implied_a < 0d0) then     ! negative wholesale funding
                            v_temp = penalty
                        elseif (implied_c < 0d0) then     ! negative securities
                            v_temp = penalty
                        elseif ( sgrid(mm) < (1d0+hair)*implied_a) then  ! collateral constraint
                            v_temp = penalty
                        else               ! if all constraints satisfied

                            ! initialize RHS of bellman equation
                            RHS = 0d0

                            ! for each possible funding shock
                            do qq=1,size(delta)

                                ! compute cash liquidations
                                excess_cash = implied_c

                                ! check liquidity default
                                if ( excess_cash + price(1)*sgrid(mm) >= (delta(qq)+shock(1))*implied_a) then

                                    c_tilde =  minval( (/ (delta(qq)+shock(1))*implied_a, excess_cash /) )

                                    if ( c_tilde == excess_cash ) then  ! if not enough cash to settle
                                        ! have to liquidate securities
                                        s_tilde = ( (delta(qq)+shock(1))*implied_a - excess_cash )/price(1)
                                    else
                                        ! if enough cash, security liquidations = 0
                                        s_tilde = 0d0
                                    endif

                                    do rr=1,size(Rl)        ! for each next period loan shock
                                        do ss=1,dbar_size  ! and each deposit funding shock

                                            ! compute networth
                                            net_int = (Rl(rr)-shock(2)-1d0)*lgrid(ll) + i_s*(sgrid(mm)-s_tilde) - &
                                                      (Rd-1d0)*dbar(grid_idx(2),grid_idx(3)) - &
                                                      (price(2)-1d0)*(1d0-delta(qq)-shock(1))*implied_a

                                            stock = lgrid(ll) + (sgrid(mm)-s_tilde) + (implied_c-c_tilde) - &
                                                    dbar(grid_idx(2),grid_idx(3)) - (1d0-delta(qq)-shock(1))*implied_a

                                            if ( net_int > 0d0) then
                                                net_temp = (1d0-tax)*net_int + stock
                                            else
                                                net_temp = net_int + stock
                                            endif

                                            if (net_temp>0d0) then

                                                ! linear interpolate for value
                                                call linint_Grow( net_temp, bank_na(grid_idx(2)), bank_nb(grid_idx(2)),&
                                                                    growth,bank_nlen-1, ixl, ixr, varphi)

                                                v_temp = ( varphi*vp(ixl+1,grid_idx(2),ss) + &
                                                            (1d0-varphi)*vp(ixr+1,grid_idx(2),ss) )

                                                RHS = RHS + prob_d(qq)*prob_l(rr)*&
                                                            dbar_prob(grid_idx(2),grid_idx(3),ss)*v_temp
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
                            s_curr_max(3)   = sgrid(mm)
                            c_curr_max(3)   = implied_c
                            d_curr_max(3)   = dbar(grid_idx(2),grid_idx(3))
                            a_curr_max(3)   = implied_a
                            div_curr_max(3) = divgrid(oo)

                        endif


                    enddo
                enddo
            enddo  ! end of policy function loop

            !
            !
            !    Case 4: binding liquidity requirement & deposit capacity constraint
            !
            !    implication:   I) d = dbar
            !
            !                  II) c = phi_lr*delta^a*a - (1-h^s*s)
            !
            !                 III) a = ( l + theta*l^2/2 + h^s*s + div - n - dbar )/( 1 - phi_lr*delta^a )
            !
            !
            !    search over: (l,c,div)   --> a --> c
            !
            !

            ! construct grids for control variables
            if ( ( grid_idx(1) == 1 ) ) then  ! if n_b(1) or first iteratior of value function, use default grids

                call grid_Cons_Equi( lgrid, la(grid_idx(2)),lb(grid_idx(2)))
                call grid_Cons_Equi( divgrid,diva(grid_idx(2)), divb(grid_idx(2)) )
                call grid_Cons_Equi( sgrid, sa(grid_idx(2)), sb(grid_idx(2)))

            else

                s_bound = maxval( (/ 0.0000001d0, spol(grid_idx(1)-1,grid_idx(2),grid_idx(3)) - back_look*delta_s /) )
                l_bound = maxval( (/ 0.0000001d0, lpol(grid_idx(1)-1,grid_idx(2),grid_idx(3)) - back_look*delta_l /) )
                pi_bound = maxval( (/ 0d0, div_pol(grid_idx(1)-1,grid_idx(2),grid_idx(3)) - delta_div /) )

                call grid_Cons_Equi( lgrid, l_bound, lpol(grid_idx(1)-1,grid_idx(2),grid_idx(3)) + delta_l )
                call grid_Cons_Equi( sgrid, s_bound, spol(grid_idx(1)-1,grid_idx(2),grid_idx(3)) + delta_s )
                call grid_Cons_Equi( divgrid, pi_bound , div_pol(grid_idx(1)-1,grid_idx(2),grid_idx(3)) + delta_div )

            endif

            ! initialize current problem
            v_curr_max(4)   = 0d0
            l_curr_max(4)   = 0d0
            d_curr_max(4)   = 0d0
            div_curr_max(4) = 0d0
            s_curr_max(4)   = 0d0
            c_curr_max(4)   = 0d0
            a_curr_max(4)   = 0d0

            ! solve problem
            do ll=1,grid_len ! for each loan
                do mm=1,grid_len ! for each security
                    do oo=1,grid_len   ! for each dividend


                        ! implied wholesale funding
!                            implied_a = ( lgrid(ll) + g(lgrid(ll),theta(grid_idx(2))) + lr_hair*sgrid(mm) + divgrid(oo) - &
!                                            bank_ngrid(grid_idx(2),grid_idx(1)) - dbar(grid_idx(2),grid_idx(3)) )/&
!                                                ( 1d0 - (phi_lr+(1d0-lr_hair)*(1d0+hair)/lr_run)*lr_run )

                        implied_a = ( lgrid(ll) + g(lgrid(ll),theta(grid_idx(2))) + lr_hair*sgrid(mm) + divgrid(oo) - &
                                        bank_ngrid(grid_idx(2),grid_idx(1)) - dbar(grid_idx(2),grid_idx(3)) )/&
                                            ( 1d0 - phi_lr*lr_run )

                        !implied_c = (phi_lr+(1d0-lr_hair)*(1d0+hair)/lr_run)*lr_run*implied_a - (1d0-lr_hair)*sgrid(mm)
                        implied_c = phi_lr*lr_run*implied_a - (1d0-lr_hair)*sgrid(mm)

                        ! implied capital ratio
                        e_temp = (lgrid(ll) + implied_c + sgrid(mm) - (dbar(grid_idx(2),grid_idx(3))+implied_a))/&
                                ( lgrid(ll) + sgrid(mm) + implied_c )

                        ! deposit capacity constraint already satisfied
                        if (e_temp > 1d0) then ! capital requirement (above)
                            v_temp = penalty
                        elseif (e_temp < ebar) then ! capital requirement (below)
                            v_temp = penalty
                        elseif (implied_c < 0d0) then     ! negative securities funding
                            v_temp = penalty
                        elseif (implied_a < 0d0) then     ! negative wholesale funding
                            v_temp = penalty
                        elseif ( sgrid(mm) < (1d0+hair)*implied_a) then  ! collateral constraint
                            v_temp = penalty
                        else               ! if all constraints satisfied

                            ! initialize RHS of bellman equation
                            RHS = 0d0

                            ! for each possible funding shock
                            do qq=1,size(delta)

                                ! compute cash liquidations
                                excess_cash = implied_c

                                ! check liquidity default
                                if ( excess_cash + price(1)*sgrid(mm) >= (delta(qq)+shock(1))*implied_a) then

                                    c_tilde =  minval( (/ (delta(qq)+shock(1))*implied_a, excess_cash /) )

                                    if ( c_tilde == excess_cash ) then  ! if not enough cash to settle
                                        ! have to liquidate securities
                                        s_tilde = ( (delta(qq)+shock(1))*implied_a - excess_cash )/price(1)
                                    else
                                        ! if enough cash, security liquidations = 0
                                        s_tilde = 0d0
                                    endif

                                    do rr=1,size(Rl)        ! for each next period loan shock
                                        do ss=1,dbar_size  ! and each deposit funding shock

                                            ! compute networth
                                            net_int = (Rl(rr)-shock(2)-1d0)*lgrid(ll) + i_s*(sgrid(mm)-s_tilde) - &
                                                      (Rd-1d0)*dbar(grid_idx(2),grid_idx(3)) - &
                                                      (price(2)-1d0)*(1d0-delta(qq)-shock(1))*implied_a

                                            stock = lgrid(ll) + (sgrid(mm)-s_tilde) + (implied_c-c_tilde) - &
                                                    dbar(grid_idx(2),grid_idx(3)) - (1d0-delta(qq)-shock(1))*implied_a

                                            if ( net_int > 0d0) then
                                                net_temp = (1d0-tax)*net_int + stock
                                            else
                                                net_temp = net_int + stock
                                            endif

                                            if (net_temp>0d0) then

                                                ! linear interpolate for value
                                                call linint_Grow( net_temp, bank_na(grid_idx(2)), bank_nb(grid_idx(2)),&
                                                                    growth,bank_nlen-1, ixl, ixr, varphi)

                                                v_temp = ( varphi*vp(ixl+1,grid_idx(2),ss) + &
                                                            (1d0-varphi)*vp(ixr+1,grid_idx(2),ss) )

                                                RHS = RHS + prob_d(qq)*prob_l(rr)*&
                                                            dbar_prob(grid_idx(2),grid_idx(3),ss)*v_temp
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
                            s_curr_max(4)   = sgrid(mm)
                            c_curr_max(4)   = implied_c
                            d_curr_max(4)   = dbar(grid_idx(2),grid_idx(3))
                            a_curr_max(4)   = implied_a
                            div_curr_max(4) = divgrid(oo)

                        endif


                    enddo
                enddo
            enddo  ! end of policy function loop


            !
            !
            !    Case 5: binding liquidity & capital requirement & deposit capacity constraint
            !
            !    implication:   I) d = dbar
            !
            !                  II) s = ( phi_lr*delta^a*a - c )/(1-h^s)
            !
            !                 III) l = ( a + dbar - s*(1-ebar) - c*(1-ebar) )/(1-ebar)
            !                        = ( a + dbar - (( phi_lr*delta^a*a - c )/(1-h^s))*(1-ebar) - c*(1-ebar) )/(1-ebar)
            !
            !                  IV) div = n + a + dbar - l - theta*l^2/2 - c - s
            !                          = n + a + dbar - l - theta*l^2/2 - c - ( phi_lr*delta^a*a - c )/(1-h^s)
            !
            !
            !    search over: (a,c)   --> s--> div --> l

            ! construct grids for control variables
            if ( ( grid_idx(1) == 1 ) ) then  ! if n_b(1) or first iteratior of value function, use default grids

                call grid_Cons_Equi( agrid, aa(grid_idx(2)), ab(grid_idx(2)))
                call grid_Cons_Equi( cgrid, ca(grid_idx(2)), cb(grid_idx(2)) )

            else

                a_bound = maxval( (/ 0.0000001d0, apol(grid_idx(1)-1,grid_idx(2),grid_idx(3)) - back_look*delta_a /) )
                c_bound = maxval( (/ 0.0000001d0, cpol(grid_idx(1)-1,grid_idx(2),grid_idx(3)) - back_look*delta_c  /) )

                call grid_Cons_Equi( agrid, a_bound, apol(grid_idx(1)-1,grid_idx(2),grid_idx(3)) + delta_a )
                call grid_Cons_Equi( cgrid, c_bound, cpol(grid_idx(1)-1,grid_idx(2),grid_idx(3)) + delta_c )

            endif

            ! initialize current problem
            v_curr_max(5)   = 0d0
            l_curr_max(5)   = 0d0
            d_curr_max(5)   = 0d0
            div_curr_max(5) = 0d0
            s_curr_max(5)   = 0d0
            c_curr_max(5)   = 0d0
            a_curr_max(5)   = 0d0

            ! solve problem
            do ll=1,grid_len ! for cash
                do oo=1,grid_len   ! for each wholesale funding

                    !implied_s = ( (phi_lr+(1d0-lr_hair)*(1d0+hair)/lr_run)*lr_run*agrid(oo) - cgrid(ll) )/(1d0-lr_hair)
                    implied_s = ( phi_lr*lr_run*agrid(oo) - cgrid(ll) )/(1d0-lr_hair)

                    implied_l = ( agrid(oo) + dbar(grid_idx(2),grid_idx(3)) - implied_s*(1d0-ebar) - cgrid(ll)*(1d0-ebar) )/&
                                                (1d0-ebar)

                    implied_div = ( bank_ngrid(grid_idx(2),grid_idx(1)) + agrid(oo) + dbar(grid_idx(2),grid_idx(3)) ) - &
                                    implied_l - g(implied_l,theta(grid_idx(2))) - cgrid(ll) - implied_s

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

                            ! compute cash liquidations
                            excess_cash = cgrid(ll)

                            ! check liquidity default
                            if ( excess_cash + price(1)*implied_s >= (delta(qq)+shock(1))*agrid(oo)) then

                                c_tilde =  minval( (/ (delta(qq)+shock(1))*agrid(oo), excess_cash /) )

                                if ( c_tilde == excess_cash ) then  ! if not enough cash to settle
                                    ! have to liquidate securities
                                    s_tilde = ( (delta(qq)+shock(1))*agrid(oo) - excess_cash )/price(1)
                                else
                                    ! if enough cash, security liquidations = 0
                                    s_tilde = 0d0
                                endif

                                do rr=1,size(Rl)        ! for each next period loan shock
                                    do ss=1,dbar_size  ! and each deposit funding shock

                                        ! compute networth
                                        net_int = (Rl(rr)-shock(2)-1d0)*implied_l + i_s*(implied_s-s_tilde) - &
                                                (Rd-1d0)*dbar(grid_idx(2),grid_idx(3)) - &
                                                (price(2)-1d0)*(1d0-delta(qq)-shock(1))*agrid(oo)

                                        stock = implied_l + (implied_s-s_tilde) + (cgrid(ll)-c_tilde) - &
                                                dbar(grid_idx(2),grid_idx(3)) - (1d0-delta(qq)-shock(1))*agrid(oo)

                                        if ( net_int > 0d0) then
                                            net_temp = (1d0-tax)*net_int + stock
                                        else
                                            net_temp = net_int + stock
                                        endif

                                        if (net_temp>0d0) then

                                            ! linear interpolate for value
                                            call linint_Grow( net_temp, bank_na(grid_idx(2)), bank_nb(grid_idx(2)),&
                                                                growth,bank_nlen-1, ixl, ixr, varphi)

                                            v_temp = ( varphi*vp(ixl+1,grid_idx(2),ss) + &
                                                        (1d0-varphi)*vp(ixr+1,grid_idx(2),ss) )

                                            RHS = RHS + prob_d(qq)*prob_l(rr)*&
                                                        dbar_prob(grid_idx(2),grid_idx(3),ss)*v_temp
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
                        c_curr_max(5)   = cgrid(ll)
                        d_curr_max(5)   = dbar(grid_idx(2),grid_idx(3))
                        a_curr_max(5)   = agrid(oo)
                        div_curr_max(5) = implied_div

                    endif

                enddo
            enddo  ! end of policy function loop


            !
            !    Case 6: no assumed binding constraints
            !

            ! construct grids for control variables
            if ( ( grid_idx(1) == 1 ) ) then  ! if n_b(1) or first iteratior of value function, use default grids

                call grid_Cons_Equi( lgrid, la(grid_idx(2)),lb(grid_idx(2)))

                if ( da(grid_idx(2)) >= dbar(grid_idx(2),grid_idx(3))) then
                    call grid_Cons_Equi( dgrid, .99d0*dbar(grid_idx(2),grid_idx(3)), dbar(grid_idx(2),grid_idx(3)) )
                else
                    call grid_Cons_Equi( dgrid, da(grid_idx(2)), dbar(grid_idx(2),grid_idx(3)))
                endif

                call grid_Cons_Equi( divgrid,diva(grid_idx(2)), divb(grid_idx(2)) )

            else

                d_bound = maxval( (/ 0.0000001d0, dpol(grid_idx(1)-1,grid_idx(2),grid_idx(3)) - back_look*delta_d /) )
                d_ubound = minval( (/ dpol(grid_idx(1)-1,grid_idx(2),grid_idx(3)) + delta_d, dbar(grid_idx(2),grid_idx(3)) /) )
                l_bound = maxval( (/ 0.0000001d0, lpol(grid_idx(1)-1,grid_idx(2),grid_idx(3)) - back_look*delta_l /) )
                pi_bound = maxval( (/ 0d0, div_pol(grid_idx(1)-1,grid_idx(2),grid_idx(3)) - delta_div /) )

                call grid_Cons_Equi( lgrid, l_bound, lpol(grid_idx(1)-1,grid_idx(2),grid_idx(3)) + delta_l )

                if (d_bound >= d_ubound) then
                    call grid_Cons_Equi( dgrid, .99d0*d_bound, d_bound )
                else
                    call grid_Cons_Equi( dgrid, d_bound, d_ubound )
                endif

                call grid_Cons_Equi( divgrid, pi_bound , div_pol(grid_idx(1)-1,grid_idx(2),grid_idx(3)) + delta_div )


            endif

            ! initialize current problem
            v_curr_max(6)   = 0d0
            l_curr_max(6)   = 0d0
            d_curr_max(6)   = 0d0
            div_curr_max(6) = 0d0
            s_curr_max(6)   = 0d0
            c_curr_max(6)   = 0d0
            a_curr_max(6)   = 0d0

            ! solve problem
            do ll=1,grid_len ! for each loan
                do mm=1,grid_len ! for each deposit
                    do oo=1,grid_len   ! for each dividend

                        implied_a = 0d0
                        implied_c = 0d0

                        ! implied securities
                        implied_s = bank_ngrid(grid_idx(2),grid_idx(1)) + dgrid(mm) - divgrid(oo) -&
                                    lgrid(ll) - g(lgrid(ll),theta(grid_idx(2)))

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
                                    net_int = (Rl(rr)-shock(2)-1d0)*lgrid(ll) + i_s*implied_s - (Rd-1d0)*dgrid(mm)

                                    stock = lgrid(ll) + implied_s - dgrid(mm)

                                    if ( net_int > 0d0) then
                                        net_temp = (1d0-tax)*net_int + stock
                                    else
                                        net_temp = net_int + stock
                                    endif

                                    if (net_temp>0d0) then

                                        ! linear interpolate for value
                                        call linint_Grow( net_temp, bank_na(grid_idx(2)), bank_nb(grid_idx(2)),&
                                                            growth,bank_nlen-1, ixl, ixr, varphi)

                                        v_temp = ( varphi*vp(ixl+1,grid_idx(2),ss) + &
                                                    (1d0-varphi)*vp(ixr+1,grid_idx(2),ss) )

                                        RHS = RHS + prob_l(rr)*dbar_prob(grid_idx(2),grid_idx(3),ss)*v_temp
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
                            c_curr_max(6)   = 0d0
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
            ctemp(ii)    = c_curr_max( v_idx(1) )

            dtemp(ii)    = d_curr_max( v_idx(1) )
            atemp(ii)    = a_curr_max( v_idx(1) )

            divtemp(ii) = div_curr_max( v_idx(1) )

        enddo   ! state space loop

        ! MPI GATHER
        call MPI_GATHER(vtemp,proc_len,MPI_REAL8,  vtemp_all,proc_len,   MPI_REAL8,0,MPI_COMM_WORLD,ierr)
        call MPI_GATHER(ltemp,proc_len,MPI_REAL8,  ltemp_all,proc_len,   MPI_REAL8,0,MPI_COMM_WORLD,ierr)
        call MPI_GATHER(dtemp,proc_len,MPI_REAL8,  dtemp_all,proc_len,   MPI_REAL8,0,MPI_COMM_WORLD,ierr)
        call MPI_GATHER(stemp,proc_len,MPI_REAL8,  stemp_all,proc_len,   MPI_REAL8,0,MPI_COMM_WORLD,ierr)
        call MPI_GATHER(ctemp,proc_len,MPI_REAL8,  ctemp_all,proc_len,   MPI_REAL8,0,MPI_COMM_WORLD,ierr)
        call MPI_GATHER(atemp,proc_len,MPI_REAL8,  atemp_all,proc_len,   MPI_REAL8,0,MPI_COMM_WORLD,ierr)
        call MPI_GATHER(divtemp,proc_len,MPI_REAL8,divtemp_all,proc_len, MPI_REAL8,0,MPI_COMM_WORLD,ierr)

        ! MPI BARRIER
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)

        if ( my_id == 0 ) then  ! if root

            ! construct dimension-based policy functions
            vrec    = transform_f(vtemp_all)
            lrec    = transform_f(ltemp_all)
            srec    = transform_f(stemp_all)
            crec    = transform_f(ctemp_all)
            drec    = transform_f(dtemp_all)
            div_rec = transform_f(divtemp_all)
            arec    = transform_f(atemp_all)

        endif

        ! broadcast value function error
        call MPI_BCAST( error,1,MPI_REAL8,0,MPI_COMM_WORLD,ierr  )

        ! broadcast policy functions
        call MPI_BCAST(vtemp_all,  proc_len*nprocs,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(ltemp_all,  proc_len*nprocs,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(dtemp_all,  proc_len*nprocs,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(stemp_all,  proc_len*nprocs,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(atemp_all,  proc_len*nprocs,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(ctemp_all,  proc_len*nprocs,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(divtemp_all,proc_len*nprocs,MPI_REAL8,0,MPI_COMM_WORLD,ierr)

        if ( my_id /= 0) then

            ! construct dimension-based policy functions
            vrec = transform_f(vtemp_all)
            lrec = transform_f(ltemp_all)
            drec = transform_f(dtemp_all)
            srec = transform_f(stemp_all)
            arec = transform_f(atemp_all)
            crec = transform_f(ctemp_all)
            div_rec = transform_f(divtemp_all)

        endif

        ! MPI BARRIER
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)

    end subroutine


    !Initialize grids and constraint sets
    subroutine initialize()

        implicit none

        integer :: ii,jj,kk,iterator
        real*8, dimension(bank_nlen) :: temp
        real*8, dimension(bank_nlen*size(theta)*dbar_size) :: holder

        ! inidividual bank networth grid
        do ii=1,size(theta)
            call grid_Cons_Grow(bank_ngrid(ii,:),bank_na(ii),bank_nb(ii),growth)
        enddo

        ! initialize value function
        do ii = 1,bank_nlen
            do jj=1,size(theta)
                do kk=1,dbar_size
                    v(ii,jj,kk) =  bank_ngrid(jj,ii)
                enddo
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

!        if (my_id==0) then
!            call plot( temp, bank_ngrid(1,:))
!            call plot( temp, bank_ngrid(2,:))
!            call execplot(title='grid',filename='grids',filetype='png',output='grids')
!        endif

        ! discretize idiosyncratic bank  processes
!        call log_norm( Rl, prob_l, l_mu, l_sigma2)
!        Rl = 2d0 - Rl
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
        dbar(1,:) = dbar_mu(1)*( 1d0 + (/ -0.45188876,-0.27113325,-0.09037775,0.09037775,0.27113325,0.45188876 /) )
        dbar(2,:) = dbar_mu(2)*( 1d0 + (/ -0.41676404,-0.25005843,-0.08335281,0.08335281,0.25005843,0.41676404 /) )
        dbar(3,:) = dbar_mu(3)*( 1d0 +(/  -0.22691827,-0.13615096,-0.04538365,0.04538365,0.13615096,0.22691827 /) )

        dbar_prob(1,1,:) = (/ 3.2277057e-01,3.8918364e-01,2.3078261e-01,5.2564357e-02,4.5494089e-03,1.4941799e-04 /)
        dbar_prob(1,2,:) = (/ 1.3758505e-01,3.3365003e-01,3.5692069e-01,1.4719323e-01,2.3233293e-02,1.4177126e-03 /)
        dbar_prob(1,3,:) = (/ 4.2489315e-02,1.9840338e-01,3.8296720e-01,2.8515319e-01,8.1694866e-02,9.2920457e-03 /)
        dbar_prob(1,4,:) = (/ 9.2920457e-03,8.1694866e-02,2.8515319e-01,3.8296720e-01,1.9840338e-01,4.2489315e-02 /)
        dbar_prob(1,5,:) = (/ 1.4177126e-03,2.3233293e-02,1.4719323e-01,3.5692069e-01,3.3365003e-01,1.3758505e-01 /)
        dbar_prob(1,6,:) = (/ 1.4941799e-04,4.5494089e-03,5.2564357e-02,2.3078261e-01,3.8918364e-01,3.2277057e-01 /)

        dbar_prob(2,1,:) = (/ 3.6355770e-01,4.0356710e-01,1.9755536e-01,3.3369986e-02,1.9130819e-03,3.6768955e-05 /)
        dbar_prob(2,2,:) = (/ 1.4184205e-01,3.6066970e-01,3.5845480e-01,1.2376471e-01,1.4674066e-02,5.9467011e-04 /)
        dbar_prob(2,3,:) = (/ 3.6311337e-02,2.0043165e-01,4.0441560e-01,2.8389830e-01,6.9047132e-02,5.8959822e-03 /)
        dbar_prob(2,4,:) = (/ 5.8959822e-03,6.9047132e-02,2.8389830e-01,4.0441560e-01,2.0043165e-01,3.6311337e-02 /)
        dbar_prob(2,5,:) = (/ 5.9467011e-04,1.4674066e-02,1.2376471e-01,3.5845480e-01,3.6066970e-01,1.4184205e-01 /)
        dbar_prob(2,6,:) = (/ 3.6768955e-05,1.9130819e-03,3.3369986e-02,1.9755536e-01,4.0356710e-01,3.6355770e-01 /)

        dbar_prob(3,1,:) = (/ 3.1194917e-01,3.8448179e-01,2.3914483e-01,5.8608838e-02,5.6057289e-03,2.0963264e-04 /)
        dbar_prob(3,2,:) = (/ 1.3614670e-01,3.2668569e-01,3.5607137e-01,1.5339403e-01,2.5953103e-02,1.7491070e-03 /)
        dbar_prob(3,3,:) = (/ 4.4072976e-02,1.9765215e-01,3.7760875e-01,2.8526807e-01,8.5029134e-02,1.0368924e-02 /)
        dbar_prob(3,4,:) = (/ 1.0368924e-02,8.5029134e-02,2.8526807e-01,3.7760875e-01,1.9765215e-01,4.4072976e-02 /)
        dbar_prob(3,5,:) = (/ 1.7491070e-03,2.5953103e-02,1.5339403e-01,3.5607137e-01,3.2668569e-01,1.3614670e-01 /)
        dbar_prob(3,6,:) = (/ 2.0963264e-04,5.6057289e-03,5.8608838e-02,2.3914483e-01,3.8448179e-01,3.1194917e-01 /)

    end subroutine

!    subroutine solve_ge()
!
!        implicit none
!
!        real*8, parameter :: alpha_nm = 1d0
!        real*8, parameter :: gamma_nm = 2d0
!        real*8, parameter :: rho_nm = .5d0
!        real*8, parameter :: sigma_nm = .5d0
!        integer, parameter :: max_iter = 10
!        real*8, parameter :: error_tol = 1d0
!
!        integer :: ii, iterator
!        integer, dimension(1) :: idx_hold
!
!        real*8, dimension(3,2) :: xvals_nm
!        real*8, dimension(3) :: fvals_nm
!        real*8 :: temp_error, fr, fe, fc
!        real*8, dimension(2) :: xo, xr, xe, xc
!        real*8 :: pstar_implied, Ra_implied
!
!        if (my_id==0) then
!            open(12,file='nm_output.csv')
!        endif
!
!
!        ! implements Nelder-Mead method to search for GE prices
!        !
!        !   - 2 equilibrium prices, thus initial simplex x(3,2) for 3 guesses and corresponding
!        !           values f(3)
!        !
!        !   - prices: Ra, pstar
!        !
!        !
!        !
!
!        ! set initial points to evaluate
!        xvals_nm(1,:) = (/ 1.0035159602d0, 0.97237587d0 /)
!        xvals_nm(2,:) = (/ 1.0037d0, 0.97237587d0 /)
!        xvals_nm(3,:) = (/ 1.0035159602d0, .975d0 /)
!
!
!        if (my_id ==0) then
!            write(12,*) 'Initial Evaluation'
!            flush(12)
!        endif
!
!        ! evaluate initial points
!        do ii=1,size(fvals_nm)
!
!            ! set parameters
!            Ra    = xvals_nm(ii,1)
!            pstar = xvals_nm(ii,2)
!
!            ! solve model
!            call vfi_bank()
!            call balance_pols()
!            call get_dist()
!            Ra_implied = implied_ra_rate()
!            Agg_Liq = Aggregate_Liq()
!            pstar_implied = alpha*(Agg_Liq(2)+omega_outside)**(alpha-1d0)
!
!            ! evaluate point
!            fvals_nm(ii) = 100d0*100d0*( abs(Ra-Ra_implied)/Ra + abs(pstar-pstar_implied)/pstar )/2d0
!
!            if (my_id ==0) then
!                write(12,*) 'Initial Guess:', ii
!                write(12,*) 'Ra:', Ra
!                write(12,*) 'pstar:', pstar
!                write(12,*) 'Ra implied:', Ra_implied
!                write(12,*) 'pstar implied:', pstar_implied
!                write(12,*) 'value:', fvals_nm(ii)
!                flush(12)
!            endif
!
!        enddo
!
!        ! initialize error
!        temp_error = 10000d0
!        iterator = 0
!
!        do while ( (temp_error > error_tol) .and.  (iterator < max_iter) )
!
!            iterator = iterator + 1
!
!            if (my_id ==0) then
!                write(12,*) 'Main Iteration:', iterator
!                flush(12)
!            endif
!
!            ! sort points/values by smallest to largest
!            call f_order(fvals_nm,xvals_nm)
!
!            ! compute centroid
!            xo = ( xvals_nm(1,:) + xvals_nm(2,:) + xvals_nm(3,:) )/size(fvals_nm)
!
!            ! reflect
!            xr = xo + alpha_nm*( xo - xvals_nm(size(fvals_nm),:) )
!
!            if ( xr(2) > no_liq_price ) then
!                xr(2) = no_liq_price
!            endif
!
!            if ( xr(1) < Rd ) then
!                xr(1) = Rd + .0001d0
!            endif
!
!            Ra = xr(1)
!            pstar = xr(2)
!
!            ! solve model
!            call vfi_bank()
!            call balance_pols()
!            call get_dist()
!            Ra_implied = implied_ra_rate()
!            Agg_Liq = Aggregate_Liq()
!            pstar_implied = alpha*(Agg_Liq(2)+omega_outside)**(alpha-1d0)
!
!            ! evaluate point
!            fr =  100d0*100d0*( abs(Ra-Ra_implied)/Ra + abs(pstar-pstar_implied)/pstar )/2d0
!
!            if (my_id ==0) then
!                write(12,*) 'Reflected Point:', xr
!                write(12,*) 'Ra implied:', Ra_implied
!                write(12,*) 'pstar implied:', pstar_implied
!                write(12,*) 'value:', fr
!                flush(12)
!            endif
!
!            ! if reflect greater than 1 and less than n
!            if ((fr.ge.fvals_nm(1)).AND.(fr < fvals_nm(size(fvals_nm)-1))) then
!
!                ! replace n+1 with reflect
!                xvals_nm(size(fvals_nm),:) = xr
!                fvals_nm(size(fvals_nm))   = fr
!
!                temp_error = minval( fvals_nm )
!
!            ! else if reflect less than 1
!            else if (fr < fvals_nm(1) ) then
!
!                !expand
!                xe = xo + gamma_nm*(xr - xo)
!
!                if ( xe(2) > no_liq_price ) then
!                    xe(2) = no_liq_price
!                endif
!
!                if ( xe(1) < Rd ) then
!                    xe(1) = Rd + .0001d0
!                endif
!
!                Ra    = xe(1)
!                pstar = xe(2)
!
!                ! solve model
!                call vfi_bank()
!                call balance_pols()
!                call get_dist()
!                Ra_implied = implied_ra_rate()
!                Agg_Liq = Aggregate_Liq()
!                pstar_implied = alpha*(Agg_Liq(2)+omega_outside)**(alpha-1d0)
!
!                ! evaluate point
!                fe = 100d0*100d0*( abs(Ra-Ra_implied)/Ra + abs(pstar-pstar_implied)/pstar )/2d0
!
!                if (my_id ==0) then
!                    write(12,*) 'Expanded Point:', xe
!                    write(12,*) 'Ra implied:', Ra_implied
!                    write(12,*) 'pstar implied:', pstar_implied
!                    write(12,*) 'value:', fe
!                    flush(12)
!                endif
!
!
!                ! if expand is better than reflect
!                if (fe < fr) then
!                    ! replace n+1 with expand, go to 1
!                    xvals_nm(size(fvals_nm),:) = xe
!                    fvals_nm(size(fvals_nm)) = fe
!
!                    temp_error = minval( fvals_nm )
!
!                else
!                    ! replace n+1 with reflect, go to 1
!                    xvals_nm(size(fvals_nm),:) = xr
!                    fvals_nm(size(fvals_nm))   = fr
!
!                    temp_error = minval( fvals_nm )
!
!                endif
!
!            !else (case when reflect greater than xn)
!            else
!                ! contract
!                xc = xo + rho_nm*( xvals_nm(size(fvals_nm),:) - xo )
!
!                if ( xc(2) > no_liq_price ) then
!                    xc(2) = no_liq_price
!                endif
!
!                if ( xc(1) < Rd ) then
!                    xc(1) = Rd + .0001d0
!                endif
!
!                Ra    = xc(1)
!                pstar = xc(2)
!
!                ! solve model
!                call vfi_bank()
!                call balance_pols()
!                call get_dist()
!                Ra_implied = implied_ra_rate()
!                Agg_Liq = Aggregate_Liq()
!                pstar_implied = alpha*(Agg_Liq(2)+omega_outside)**(alpha-1d0)
!
!                ! evaluate point
!                fc =  100d0*100d0*( abs(Ra-Ra_implied)/Ra + abs(pstar-pstar_implied)/pstar )/2d0
!
!                if (my_id ==0) then
!                    write(12,*) 'Contracted Point:', xc
!                    write(12,*) 'Ra implied:', Ra_implied
!                    write(12,*) 'pstar implied:', pstar_implied
!                    write(12,*) 'value:', fc
!                    flush(12)
!                endif
!
!                ! if contract less than n+1
!                if ( fc < fvals_nm(size(fvals_nm)) ) then
!                    ! replace n+1 with contract, go to 1
!                    xvals_nm(size(fvals_nm),:) = xc
!                    fvals_nm(size(fvals_nm))   = fc
!
!                    temp_error = minval( fvals_nm )
!
!                else
!                    ! shrink all points, except 1, go to 1
!                    do ii=2,size(fvals_nm)
!                        xvals_nm(ii,:) = xvals_nm(1,:) + sigma_nm*( xvals_nm(ii,:) - xvals_nm(1,:) )
!
!                        if ( xvals_nm(ii,2) > no_liq_price ) then
!                            xvals_nm(ii,2) = no_liq_price
!                        endif
!
!                        if ( xvals_nm(ii,1) < Rd ) then
!                            xvals_nm(ii,1) = Rd + .0001d0
!                        endif
!
!                        Ra    = xvals_nm(ii,1)
!                        pstar = xvals_nm(ii,2)
!
!                        ! solve model
!                        call vfi_bank()
!                        call balance_pols()
!                        call get_dist()
!                        Ra_implied = implied_ra_rate()
!                        Agg_Liq = Aggregate_Liq()
!                        pstar_implied = alpha*(Agg_Liq(2)+omega_outside)**(alpha-1d0)
!
!                        ! evaluate point
!                        fvals_nm(ii) = 100d0*100d0*( abs(Ra-Ra_implied)/Ra + abs(pstar-pstar_implied)/pstar )/2d0
!
!                        if (my_id ==0) then
!                            write(12,*) 'Shrink Point:', xvals_nm(ii,:)
!                            write(12,*) 'Ra implied:', Ra_implied
!                            write(12,*) 'pstar implied:', pstar_implied
!                            write(12,*) 'value:', fvals_nm(ii)
!                            flush(12)
!                        endif
!
!
!                    enddo
!
!                    temp_error = minval( fvals_nm )
!
!                endif
!
!            endif
!
!        enddo
!
!        ! locate min
!        idx_hold = minloc(fvals_nm)
!
!        ! set prices
!        Ra    = xvals_nm( idx_hold(1), 1)
!        pstar = xvals_nm( idx_hold(1), 2)
!
!        if (my_id ==0 ) then
!
!            write(12,*) 'minimizing point!'
!            write(12,*) 'Ra:', Ra
!            write(12,*) 'pstar:',pstar
!            write(12,*) 'value:', fvals_nm(idx_hold(1))
!            flush(12)
!
!        endif
!
!        ! solve economy with correct eq'm prices
!        call vfi_bank()
!        call balance_pols()
!        call get_dist()
!
!
!    end subroutine


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
            pstar_implied = alpha*(Agg_Liq(2)+omega_outside)**(alpha-1d0)

            if (my_id==0) then
                write(*,*)
                write(*,*) 'Agg Liq:', Agg_Liq(2)
                write(*,*) 'omega:', omega_outside
                write(*,*) 'Agg Liq:', Agg_Liq(2)
                write(*,*)
                write(*,*) 'Security liquidations:',Agg_Liq(2)
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
        real*8, dimension(bank_nlen,size(theta),dbar_size,size(delta)) :: ctilde, stilde
        real*8 :: excess_cash

        do ii=1,bank_nlen
            do jj=1,size(theta)
                do kk=1,dbar_size

                    excess_cash = cpol(ii,jj,kk)

                    do ll=1,size(delta)

                        if (delta(ll)*apol(ii,jj,kk) > excess_cash + pstar*spol(ii,jj,kk)) then ! liquidity default

                            ctilde(ii,jj,kk,ll) = excess_cash
                            stilde(ii,jj,kk,ll) = spol(ii,jj,kk)

                        else

                            ctilde(ii,jj,kk,ll) = minval( (/ delta(ll)*apol(ii,jj,kk), excess_cash /) )
                            stilde(ii,jj,kk,ll) = (delta(ll)*apol(ii,jj,kk) - ctilde(ii,jj,kk,ll))/pstar

                        endif

                    enddo
                enddo
            enddo
        enddo

    end subroutine

    ! VFI iterator for stationary equilibrium, calls bank_solve_stationary()
    subroutine vfi_bank()

        implicit none

        ! input/output variables

        ! VFI variables/parameters
        real*8, parameter :: tol = 1e-14
        real*8            :: error
        integer           :: iterator, tt
        integer, parameter :: it_max = 75

        ! Solving bank problem
        integer :: ii, jj, kk, ll, mm, nn, oo, pp, qq,rr,ss
        real*8 :: v_temp, RHS, net_temp, e_temp
        real*8, dimension(6) :: v_curr_max, l_curr_max, d_curr_max, div_curr_max, c_curr_max, s_curr_max, a_curr_max
        real*8 :: a_bound, l_bound, s_bound, d_bound,c_bound, liq_temp, implied_a, &
                  implied_s, implied_c, d_ubound, implied_d, implied_l, implied_div
        real*8 :: c_tilde, s_tilde, excess_cash  ! liquidation values
        integer, dimension(1) :: v_idx
        real*8 :: net_int, stock, pi_bound
        real*8 :: temp_pen

        ! linear interpolation variables
        integer :: ixl, ixr
        real*8  :: varphi

        ! spline interpolation variables
        real*8, dimension( size(theta), dbar_size, bank_nlen+2 ) :: coeff

        !---------------!
        !   mpi stuff   !
        !---------------!
        real*8, allocatable, dimension(:) :: ltemp, dtemp, stemp, atemp, divtemp, vtemp, ctemp
        real*8, allocatable, dimension(:) :: ltemp_all, dtemp_all, stemp_all, atemp_all, divtemp_all, vtemp_all, ctemp_all
        integer, dimension(3) :: grid_idx   ! 3 states (nb,theta,dbar)
        integer, allocatable, dimension(:,:,:) :: idx_full_grid

        integer :: total_grid, proc_len

        ! total size of state space grid
        total_grid = bank_nlen*size(theta)*dbar_size

        if ( mod(total_grid,nprocs) == 0) then  ! perfect fit
            proc_len = total_grid/nprocs
        else
            proc_len = total_grid/nprocs + 1     ! the last proccessor will have some extra/empty slots
        endif

        ! construct temporary policy functions
        allocate( ltemp(proc_len) )
        allocate( dtemp(proc_len) )
        allocate( stemp(proc_len) )
        allocate( ctemp(proc_len) )
        allocate( atemp(proc_len) )
        allocate( divtemp(proc_len) )
        allocate( vtemp(proc_len) )

        allocate( ltemp_all(proc_len*nprocs) )
        allocate( dtemp_all(proc_len*nprocs) )
        allocate( stemp_all(proc_len*nprocs) )
        allocate( ctemp_all(proc_len*nprocs) )
        allocate( atemp_all(proc_len*nprocs) )
        allocate( divtemp_all(proc_len*nprocs) )
        allocate( vtemp_all(proc_len*nprocs) )

        allocate( idx_full_grid(nprocs,proc_len,3) )

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
                !                 III) a = 0, c = 0
                !
                !    search over: (l,div)
                !
                !

                ! construct grids for control variables
                if ( ( grid_idx(1) == 1 ).or.( iterator == 1 ) ) then  ! if n_b(1) or first iteratior of value function, use default grids

                    call grid_Cons_Equi( lgrid, la(grid_idx(2)),lb(grid_idx(2))   )
                    call grid_Cons_Equi( divgrid,diva(grid_idx(2)), divb(grid_idx(2)) )

                else

                    l_bound = maxval( (/ 0.0000001d0, lpol(grid_idx(1)-1,grid_idx(2),grid_idx(3)) - back_look*delta_l /) )
                    pi_bound = maxval( (/ 0d0, div_pol(grid_idx(1)-1,grid_idx(2),grid_idx(3)) - delta_div /) )

                    call grid_Cons_Equi( lgrid, l_bound, lpol(grid_idx(1)-1,grid_idx(2),grid_idx(3)) + delta_l )
                    call grid_Cons_Equi( divgrid, pi_bound , div_pol(grid_idx(1)-1,grid_idx(2),grid_idx(3)) + delta_div )

                endif

                ! initialize current problem
                v_curr_max(1)   = 0d0
                l_curr_max(1)   = 0d0
                d_curr_max(1)   = 0d0
                div_curr_max(1) = 0d0
                s_curr_max(1)   = 0d0
                c_curr_max(1)   = 0d0
                a_curr_max(1)   = 0d0

                ! solve problem
                do ll=1,grid_len ! for each loan
                    do oo=1,grid_len   ! for each dividend

                        ! implied securities
                        implied_s =  ( bank_ngrid(grid_idx(2),grid_idx(1)) - divgrid(oo) - &
                                       ebar*lgrid(ll) - g(lgrid(ll),theta(grid_idx(2))) )/ebar

                        implied_d = (lgrid(ll) + implied_s )*(1d0-ebar)

                        ! implied wholesale funding
                        implied_a = 0d0
                        implied_c = 0d0

                        ! liquidity requirement naturally satisfied

                        if (implied_d < 0d0) then
                            v_temp = penalty
                        elseif (implied_d > dbar(grid_idx(2),grid_idx(3))) then   ! liquidity requirement
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
                                        call linint_Grow( net_temp, bank_na(grid_idx(2)), bank_nb(grid_idx(2)),&
                                                            growth,bank_nlen-1, ixl, ixr, varphi)

                                        v_temp = ( varphi*v(ixl+1,grid_idx(2),ss) + &
                                                    (1d0-varphi)*v(ixr+1,grid_idx(2),ss) )

                                        RHS = RHS + prob_l(rr)*dbar_prob(grid_idx(2),grid_idx(3),ss)*v_temp
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
                            c_curr_max(1)   = 0d0
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
                !                  II) a = c + s + l + theta*l**2/2 + div - n - dbar
                !
                !    search over: (l,c,s,div) --> a
                !
                !
                ! construct grids for control variables
                if ( ( grid_idx(1) == 1 ).or.( iterator == 1 ) ) then  ! if n_b(1) or first iteratior of value function, use default grids

                    call grid_Cons_Equi( lgrid, la(grid_idx(2)),lb(grid_idx(2)))
                    call grid_Cons_Equi( sgrid, sa(grid_idx(2)), sb(grid_idx(2)))
                    call grid_Cons_Equi( cgrid, ca(grid_idx(2)), cb(grid_idx(2)))
                    call grid_Cons_Equi( divgrid,diva(grid_idx(2)), divb(grid_idx(2)) )

                else

                    l_bound = maxval( (/ 0.0000001d0, lpol(grid_idx(1)-1,grid_idx(2),grid_idx(3)) - back_look*delta_l /) )
                    s_bound = maxval( (/ 0.0000001d0, spol(grid_idx(1)-1,grid_idx(2),grid_idx(3)) - back_look*delta_s /) )
                    c_bound = maxval( (/ 0.0000001d0, cpol(grid_idx(1)-1,grid_idx(2),grid_idx(3)) - back_look*delta_c /) )
                    pi_bound = maxval( (/ 0d0, div_pol(grid_idx(1)-1,grid_idx(2),grid_idx(3)) - delta_div /) )

                    call grid_Cons_Equi( lgrid, l_bound, lpol(grid_idx(1)-1,grid_idx(2),grid_idx(3)) + delta_l )
                    call grid_Cons_Equi( sgrid, s_bound, spol(grid_idx(1)-1,grid_idx(2),grid_idx(3)) + delta_s )
                    call grid_Cons_Equi( cgrid, c_bound, cpol(grid_idx(1)-1,grid_idx(2),grid_idx(3)) + delta_c )
                    call grid_Cons_Equi( divgrid, pi_bound , div_pol(grid_idx(1)-1,grid_idx(2),grid_idx(3)) + delta_div )

                endif

                ! initialize current problem
                v_curr_max(2)   = 0d0
                l_curr_max(2)   = 0d0
                d_curr_max(2)   = 0d0
                div_curr_max(2) = 0d0
                s_curr_max(2)   = 0d0
                c_curr_max(2)   = 0d0
                a_curr_max(2)   = 0d0

                ! solve problem
                do ll=1,grid_len ! for each loan
                    do mm=1,grid_len ! for each cash
                        do nn=1,grid_len  ! for each security
                            do oo=1,grid_len   ! for each dividend

                                ! implies wholesale funding
                                implied_a = cgrid(mm) + sgrid(nn) + lgrid(ll) + g(lgrid(ll),theta(grid_idx(2))) + &
                                            divgrid(oo) - &
                                            bank_ngrid(grid_idx(2),grid_idx(1)) - dbar(grid_idx(2),grid_idx(3))

                                ! implied capital ratio
                                e_temp = (lgrid(ll) + cgrid(mm) + sgrid(nn) - (dbar(grid_idx(2),grid_idx(3))+implied_a))/&
                                        ( lgrid(ll) + sgrid(nn) + cgrid(mm) )

                                ! implies a liquidity ratio
                                liq_temp = ( cgrid(mm) + (1d0-lr_hair)*sgrid(nn) )/( lr_run*implied_a ) !- &
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

                                        ! compute cash liquidations
                                        excess_cash = cgrid(mm)

                                        ! check liquidity default
                                        if ( excess_cash + pstar*sgrid(nn) >= delta(qq)*implied_a) then

                                            c_tilde =  minval( (/ delta(qq)*implied_a, excess_cash /) )

                                            if ( c_tilde == excess_cash ) then  ! if not enough cash to settle
                                                ! have to liquidate securities
                                                s_tilde = ( delta(qq)*implied_a - excess_cash )/pstar
                                            else
                                                ! if enough cash, security liquidations = 0
                                                s_tilde = 0d0
                                            endif

                                            do rr=1,size(Rl)        ! for each next period loan shock
                                                do ss=1,dbar_size  ! and each deposit funding shock

                                                    ! compute networth
                                                    net_int = (Rl(rr)-1d0)*lgrid(ll) + i_s*(sgrid(nn)-s_tilde) - &
                                                              (Rd-1d0)*dbar(grid_idx(2),grid_idx(3)) - &
                                                              (Ra-1d0)*(1d0-delta(qq))*implied_a

                                                    stock = lgrid(ll) + (sgrid(nn)-s_tilde) + (cgrid(mm)-c_tilde) - &
                                                            dbar(grid_idx(2),grid_idx(3)) - (1d0-delta(qq))*implied_a

                                                    if ( net_int > 0d0) then
                                                        net_temp = (1d0-tax)*net_int + stock
                                                    else
                                                        net_temp = net_int + stock
                                                    endif

                                                    if (net_temp>0d0) then
                                                        ! linear interpolate for value
                                                        call linint_Grow( net_temp, bank_na(grid_idx(2)), bank_nb(grid_idx(2)),&
                                                                            growth,bank_nlen-1, ixl, ixr, varphi)

                                                        v_temp = ( varphi*v(ixl+1,grid_idx(2),ss) + &
                                                                    (1d0-varphi)*v(ixr+1,grid_idx(2),ss) )

                                                        RHS = RHS + prob_d(qq)*prob_l(rr)*&
                                                                    dbar_prob(grid_idx(2),grid_idx(3),ss)*v_temp
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
                                    c_curr_max(2)   = cgrid(mm)
                                    d_curr_max(2)   = dbar(grid_idx(2),grid_idx(3))
                                    a_curr_max(2)   = implied_a
                                    div_curr_max(2) = divgrid(oo)

                                endif


                            enddo
                        enddo
                    enddo
                enddo ! end of policy function loop

                !
                !
                !    Case 3: binding capital & deposit capacity constraint
                !
                !    implication:   I) d = dbar
                !
                !                  II) a = ( l + s + c )*(1-ebar) - dbar
                !
                !                 III) c = ( n - theta*l^2/2 - div - ebar*s - ebar*l )/ebar
                !
                !
                !    search over: (l,s,div)   --> c --> a
                !
                !
                ! construct grids for control variables
                if ( ( grid_idx(1) == 1 ).or.( iterator == 1 ) ) then  ! if n_b(1) or first iteratior of value function, use default grids

                    call grid_Cons_Equi( lgrid, la(grid_idx(2)),lb(grid_idx(2)))
                    call grid_Cons_Equi( divgrid,diva(grid_idx(2)), divb(grid_idx(2)) )
                    call grid_Cons_Equi( sgrid, sa(grid_idx(2)), sb(grid_idx(2)))

                else

                    s_bound = maxval( (/ 0.0000001d0, spol(grid_idx(1)-1,grid_idx(2),grid_idx(3)) - back_look*delta_s /) )
                    l_bound = maxval( (/ 0.0000001d0, lpol(grid_idx(1)-1,grid_idx(2),grid_idx(3)) - back_look*delta_l /) )
                    pi_bound = maxval( (/ 0d0, div_pol(grid_idx(1)-1,grid_idx(2),grid_idx(3)) - delta_div /) )

                    call grid_Cons_Equi( lgrid, l_bound, lpol(grid_idx(1)-1,grid_idx(2),grid_idx(3)) + delta_l )
                    call grid_Cons_Equi( sgrid, s_bound, spol(grid_idx(1)-1,grid_idx(2),grid_idx(3)) + delta_s )
                    call grid_Cons_Equi( divgrid, pi_bound , div_pol(grid_idx(1)-1,grid_idx(2),grid_idx(3)) + delta_div )

                endif

                ! initialize current problem
                v_curr_max(3)   = 0d0
                l_curr_max(3)   = 0d0
                d_curr_max(3)   = 0d0
                div_curr_max(3) = 0d0
                s_curr_max(3)   = 0d0
                c_curr_max(3)   = 0d0
                a_curr_max(3)   = 0d0

                ! solve problem
                do ll=1,grid_len ! for each loan
                    do mm=1,grid_len ! for each security
                        do oo=1,grid_len   ! for each dividend

                            ! implied securities [ n - Phi(div) - ebar*wl*l - theta*l^2/2 - ebar*wc*p_c*c ]/(  ebar*ws*p_s )
                            implied_c =  ( bank_ngrid(grid_idx(2),grid_idx(1)) - divgrid(oo) - &
                                            ebar*lgrid(ll) - g(lgrid(ll),theta(grid_idx(2))) - &
                                            ebar*sgrid(mm) )/( ebar )

                            ! implied wholesale funding a = l*(1-ebar*wl) + p_s*s*(1-ebar*ws) + p_c*c*(1-ebar*wc) - d
                            implied_a = (lgrid(ll) + sgrid(mm) + implied_c)*(1d0-ebar) - dbar(grid_idx(2),grid_idx(3))

                            ! implied liquidity ratio
                            liq_temp = ( implied_c + (1d0-lr_hair)*sgrid(mm) )/( lr_run*implied_a ) !- &
                                                    !(1d0-lr_hair)*(1d0+hair)/lr_run

                            if (liq_temp < phi_lr) then   ! liquidity requirement
                                v_temp = penalty
                            elseif (implied_a < 0d0) then     ! negative wholesale funding
                                v_temp = penalty
                            elseif (implied_c < 0d0) then     ! negative securities
                                v_temp = penalty
                            elseif ( sgrid(mm) < (1d0+hair)*implied_a) then  ! collateral constraint
                                v_temp = penalty
                            else               ! if all constraints satisfied

                                ! initialize RHS of bellman equation
                                RHS = 0d0

                                ! for each possible funding shock
                                do qq=1,size(delta)

                                    ! compute cash liquidations
                                    excess_cash = implied_c

                                    ! check liquidity default
                                    if ( excess_cash + pstar*sgrid(mm) >= delta(qq)*implied_a) then

                                        c_tilde =  minval( (/ delta(qq)*implied_a, excess_cash /) )

                                        if ( c_tilde == excess_cash ) then  ! if not enough cash to settle
                                            ! have to liquidate securities
                                            s_tilde = ( delta(qq)*implied_a - excess_cash )/pstar
                                        else
                                            ! if enough cash, security liquidations = 0
                                            s_tilde = 0d0
                                        endif

                                        do rr=1,size(Rl)        ! for each next period loan shock
                                            do ss=1,dbar_size  ! and each deposit funding shock

                                                ! compute networth
                                                net_int = (Rl(rr)-1d0)*lgrid(ll) + i_s*(sgrid(mm)-s_tilde) - &
                                                          (Rd-1d0)*dbar(grid_idx(2),grid_idx(3)) - &
                                                          (Ra-1d0)*(1d0-delta(qq))*implied_a

                                                stock = lgrid(ll) + (sgrid(mm)-s_tilde) + (implied_c-c_tilde) - &
                                                        dbar(grid_idx(2),grid_idx(3)) - (1d0-delta(qq))*implied_a

                                                if ( net_int > 0d0) then
                                                    net_temp = (1d0-tax)*net_int + stock
                                                else
                                                    net_temp = net_int + stock
                                                endif

                                                if (net_temp>0d0) then

                                                    ! linear interpolate for value
                                                    call linint_Grow( net_temp, bank_na(grid_idx(2)), bank_nb(grid_idx(2)),&
                                                                        growth,bank_nlen-1, ixl, ixr, varphi)

                                                    v_temp = ( varphi*v(ixl+1,grid_idx(2),ss) + &
                                                                (1d0-varphi)*v(ixr+1,grid_idx(2),ss) )

                                                    RHS = RHS + prob_d(qq)*prob_l(rr)*&
                                                                dbar_prob(grid_idx(2),grid_idx(3),ss)*v_temp
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
                                s_curr_max(3)   = sgrid(mm)
                                c_curr_max(3)   = implied_c
                                d_curr_max(3)   = dbar(grid_idx(2),grid_idx(3))
                                a_curr_max(3)   = implied_a
                                div_curr_max(3) = divgrid(oo)

                            endif


                        enddo
                    enddo
                enddo  ! end of policy function loop

                !
                !
                !    Case 4: binding liquidity requirement & deposit capacity constraint
                !
                !    implication:   I) d = dbar
                !
                !                  II) c = phi_lr*delta^a*a - (1-h^s*s)
                !
                !                 III) a = ( l + theta*l^2/2 + h^s*s + div - n - dbar )/( 1 - phi_lr*delta^a )
                !
                !
                !    search over: (l,c,div)   --> a --> c
                !
                !

                ! construct grids for control variables
                if ( ( grid_idx(1) == 1 ).or.( iterator == 1 ) ) then  ! if n_b(1) or first iteratior of value function, use default grids

                    call grid_Cons_Equi( lgrid, la(grid_idx(2)),lb(grid_idx(2)))
                    call grid_Cons_Equi( divgrid,diva(grid_idx(2)), divb(grid_idx(2)) )
                    call grid_Cons_Equi( sgrid, sa(grid_idx(2)), sb(grid_idx(2)))

                else

                    s_bound = maxval( (/ 0.0000001d0, spol(grid_idx(1)-1,grid_idx(2),grid_idx(3)) - back_look*delta_s /) )
                    l_bound = maxval( (/ 0.0000001d0, lpol(grid_idx(1)-1,grid_idx(2),grid_idx(3)) - back_look*delta_l /) )
                    pi_bound = maxval( (/ 0d0, div_pol(grid_idx(1)-1,grid_idx(2),grid_idx(3)) - delta_div /) )

                    call grid_Cons_Equi( lgrid, l_bound, lpol(grid_idx(1)-1,grid_idx(2),grid_idx(3)) + delta_l )
                    call grid_Cons_Equi( sgrid, s_bound, spol(grid_idx(1)-1,grid_idx(2),grid_idx(3)) + delta_s )
                    call grid_Cons_Equi( divgrid, pi_bound , div_pol(grid_idx(1)-1,grid_idx(2),grid_idx(3)) + delta_div )

                endif

                ! initialize current problem
                v_curr_max(4)   = 0d0
                l_curr_max(4)   = 0d0
                d_curr_max(4)   = 0d0
                div_curr_max(4) = 0d0
                s_curr_max(4)   = 0d0
                c_curr_max(4)   = 0d0
                a_curr_max(4)   = 0d0

                ! solve problem
                do ll=1,grid_len ! for each loan
                    do mm=1,grid_len ! for each security
                        do oo=1,grid_len   ! for each dividend


                            ! implied wholesale funding
!                            implied_a = ( lgrid(ll) + g(lgrid(ll),theta(grid_idx(2))) + lr_hair*sgrid(mm) + divgrid(oo) - &
!                                            bank_ngrid(grid_idx(2),grid_idx(1)) - dbar(grid_idx(2),grid_idx(3)) )/&
!                                                ( 1d0 - (phi_lr+(1d0-lr_hair)*(1d0+hair)/lr_run)*lr_run )

                            implied_a = ( lgrid(ll) + g(lgrid(ll),theta(grid_idx(2))) + lr_hair*sgrid(mm) + divgrid(oo) - &
                                            bank_ngrid(grid_idx(2),grid_idx(1)) - dbar(grid_idx(2),grid_idx(3)) )/&
                                                ( 1d0 - phi_lr*lr_run )

                            !implied_c = (phi_lr+(1d0-lr_hair)*(1d0+hair)/lr_run)*lr_run*implied_a - (1d0-lr_hair)*sgrid(mm)
                            implied_c = phi_lr*lr_run*implied_a - (1d0-lr_hair)*sgrid(mm)

                            ! implied capital ratio
                            e_temp = (lgrid(ll) + implied_c + sgrid(mm) - (dbar(grid_idx(2),grid_idx(3))+implied_a))/&
                                    ( lgrid(ll) + sgrid(mm) + implied_c )

                            ! deposit capacity constraint already satisfied
                            if (e_temp > 1d0) then ! capital requirement (above)
                                v_temp = penalty
                            elseif (e_temp < ebar) then ! capital requirement (below)
                                v_temp = penalty
                            elseif (implied_c < 0d0) then     ! negative securities funding
                                v_temp = penalty
                            elseif (implied_a < 0d0) then     ! negative wholesale funding
                                v_temp = penalty
                            elseif ( sgrid(mm) < (1d0+hair)*implied_a) then  ! collateral constraint
                                v_temp = penalty
                            else               ! if all constraints satisfied

                                ! initialize RHS of bellman equation
                                RHS = 0d0

                                ! for each possible funding shock
                                do qq=1,size(delta)

                                    ! compute cash liquidations
                                    excess_cash = implied_c

                                    ! check liquidity default
                                    if ( excess_cash + pstar*sgrid(mm) >= delta(qq)*implied_a) then

                                        c_tilde =  minval( (/ delta(qq)*implied_a, excess_cash /) )

                                        if ( c_tilde == excess_cash ) then  ! if not enough cash to settle
                                            ! have to liquidate securities
                                            s_tilde = ( delta(qq)*implied_a - excess_cash )/pstar
                                        else
                                            ! if enough cash, security liquidations = 0
                                            s_tilde = 0d0
                                        endif

                                        do rr=1,size(Rl)        ! for each next period loan shock
                                            do ss=1,dbar_size  ! and each deposit funding shock

                                                ! compute networth
                                                net_int = (Rl(rr)-1d0)*lgrid(ll) + i_s*(sgrid(mm)-s_tilde) - &
                                                          (Rd-1d0)*dbar(grid_idx(2),grid_idx(3)) - &
                                                          (Ra-1d0)*(1d0-delta(qq))*implied_a

                                                stock = lgrid(ll) + (sgrid(mm)-s_tilde) + (implied_c-c_tilde) - &
                                                        dbar(grid_idx(2),grid_idx(3)) - (1d0-delta(qq))*implied_a

                                                if ( net_int > 0d0) then
                                                    net_temp = (1d0-tax)*net_int + stock
                                                else
                                                    net_temp = net_int + stock
                                                endif

                                                if (net_temp>0d0) then

                                                    ! linear interpolate for value
                                                    call linint_Grow( net_temp, bank_na(grid_idx(2)), bank_nb(grid_idx(2)),&
                                                                        growth,bank_nlen-1, ixl, ixr, varphi)

                                                    v_temp = ( varphi*v(ixl+1,grid_idx(2),ss) + &
                                                                (1d0-varphi)*v(ixr+1,grid_idx(2),ss) )

                                                    RHS = RHS + prob_d(qq)*prob_l(rr)*&
                                                                dbar_prob(grid_idx(2),grid_idx(3),ss)*v_temp
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
                                s_curr_max(4)   = sgrid(mm)
                                c_curr_max(4)   = implied_c
                                d_curr_max(4)   = dbar(grid_idx(2),grid_idx(3))
                                a_curr_max(4)   = implied_a
                                div_curr_max(4) = divgrid(oo)

                            endif


                        enddo
                    enddo
                enddo  ! end of policy function loop


                !
                !
                !    Case 5: binding liquidity & capital requirement & deposit capacity constraint
                !
                !    implication:   I) d = dbar
                !
                !                  II) s = ( phi_lr*delta^a*a - c )/(1-h^s)
                !
                !                 III) l = ( a + dbar - s*(1-ebar) - c*(1-ebar) )/(1-ebar)
                !                        = ( a + dbar - (( phi_lr*delta^a*a - c )/(1-h^s))*(1-ebar) - c*(1-ebar) )/(1-ebar)
                !
                !                  IV) div = n + a + dbar - l - theta*l^2/2 - c - s
                !                          = n + a + dbar - l - theta*l^2/2 - c - ( phi_lr*delta^a*a - c )/(1-h^s)
                !
                !
                !    search over: (a,c)   --> s--> div --> l

                ! construct grids for control variables
                if ( ( grid_idx(1) == 1 ).or.( iterator == 1 ) ) then  ! if n_b(1) or first iteratior of value function, use default grids

                    call grid_Cons_Equi( agrid, aa(grid_idx(2)), ab(grid_idx(2)))
                    call grid_Cons_Equi( cgrid, ca(grid_idx(2)), cb(grid_idx(2)) )

                else

                    a_bound = maxval( (/ 0.0000001d0, apol(grid_idx(1)-1,grid_idx(2),grid_idx(3)) - back_look*delta_a /) )
                    c_bound = maxval( (/ 0.0000001d0, cpol(grid_idx(1)-1,grid_idx(2),grid_idx(3)) - back_look*delta_c  /) )

                    call grid_Cons_Equi( agrid, a_bound, apol(grid_idx(1)-1,grid_idx(2),grid_idx(3)) + delta_a )
                    call grid_Cons_Equi( cgrid, c_bound, cpol(grid_idx(1)-1,grid_idx(2),grid_idx(3)) + delta_c )

                endif

                ! initialize current problem
                v_curr_max(5)   = 0d0
                l_curr_max(5)   = 0d0
                d_curr_max(5)   = 0d0
                div_curr_max(5) = 0d0
                s_curr_max(5)   = 0d0
                c_curr_max(5)   = 0d0
                a_curr_max(5)   = 0d0

                ! solve problem
                do ll=1,grid_len ! for cash
                    do oo=1,grid_len   ! for each wholesale funding

                        !implied_s = ( (phi_lr+(1d0-lr_hair)*(1d0+hair)/lr_run)*lr_run*agrid(oo) - cgrid(ll) )/(1d0-lr_hair)
                        implied_s = ( phi_lr*lr_run*agrid(oo) - cgrid(ll) )/(1d0-lr_hair)

                        implied_l = ( agrid(oo) + dbar(grid_idx(2),grid_idx(3)) - implied_s*(1d0-ebar) - cgrid(ll)*(1d0-ebar) )/&
                                                    (1d0-ebar)

                        implied_div = ( bank_ngrid(grid_idx(2),grid_idx(1)) + agrid(oo) + dbar(grid_idx(2),grid_idx(3)) ) - &
                                        implied_l - g(implied_l,theta(grid_idx(2))) - cgrid(ll) - implied_s

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

                                ! compute cash liquidations
                                excess_cash = cgrid(ll)

                                ! check liquidity default
                                if ( excess_cash + pstar*implied_s >= delta(qq)*agrid(oo)) then

                                    c_tilde =  minval( (/ delta(qq)*agrid(oo), excess_cash /) )

                                    if ( c_tilde == excess_cash ) then  ! if not enough cash to settle
                                        ! have to liquidate securities
                                        s_tilde = ( delta(qq)*agrid(oo) - excess_cash )/pstar
                                    else
                                        ! if enough cash, security liquidations = 0
                                        s_tilde = 0d0
                                    endif

                                    do rr=1,size(Rl)        ! for each next period loan shock
                                        do ss=1,dbar_size  ! and each deposit funding shock

                                            ! compute networth
                                            net_int = (Rl(rr)-1d0)*implied_l + i_s*(implied_s-s_tilde) - &
                                                    (Rd-1d0)*dbar(grid_idx(2),grid_idx(3)) - (Ra-1d0)*(1d0-delta(qq))*agrid(oo)

                                            stock = implied_l + (implied_s-s_tilde) + (cgrid(ll)-c_tilde) - &
                                                    dbar(grid_idx(2),grid_idx(3)) - (1d0-delta(qq))*agrid(oo)

                                            if ( net_int > 0d0) then
                                                net_temp = (1d0-tax)*net_int + stock
                                            else
                                                net_temp = net_int + stock
                                            endif

                                            if (net_temp>0d0) then

                                                ! linear interpolate for value
                                                call linint_Grow( net_temp, bank_na(grid_idx(2)), bank_nb(grid_idx(2)),&
                                                                    growth,bank_nlen-1, ixl, ixr, varphi)

                                                v_temp = ( varphi*v(ixl+1,grid_idx(2),ss) + &
                                                            (1d0-varphi)*v(ixr+1,grid_idx(2),ss) )

                                                RHS = RHS + prob_d(qq)*prob_l(rr)*&
                                                            dbar_prob(grid_idx(2),grid_idx(3),ss)*v_temp
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
                            c_curr_max(5)   = cgrid(ll)
                            d_curr_max(5)   = dbar(grid_idx(2),grid_idx(3))
                            a_curr_max(5)   = agrid(oo)
                            div_curr_max(5) = implied_div

                        endif

                    enddo
                enddo  ! end of policy function loop


                !
                !    Case 6: no assumed binding constraints
                !

                ! construct grids for control variables
                if ( ( grid_idx(1) == 1 ).or.( iterator == 1 ) ) then  ! if n_b(1) or first iteratior of value function, use default grids

                    call grid_Cons_Equi( lgrid, la(grid_idx(2)),lb(grid_idx(2)))

                    if ( da(grid_idx(2)) >= dbar(grid_idx(2),grid_idx(3))) then
                        call grid_Cons_Equi( dgrid, .99d0*dbar(grid_idx(2),grid_idx(3)), dbar(grid_idx(2),grid_idx(3)) )
                    else
                        call grid_Cons_Equi( dgrid, da(grid_idx(2)), dbar(grid_idx(2),grid_idx(3)))
                    endif

                    call grid_Cons_Equi( divgrid,diva(grid_idx(2)), divb(grid_idx(2)) )

                else

                    d_bound = maxval( (/ 0.0000001d0, dpol(grid_idx(1)-1,grid_idx(2),grid_idx(3)) - back_look*delta_d /) )
                    d_ubound = minval( (/ dpol(grid_idx(1)-1,grid_idx(2),grid_idx(3)) + delta_d, dbar(grid_idx(2),grid_idx(3)) /) )
                    l_bound = maxval( (/ 0.0000001d0, lpol(grid_idx(1)-1,grid_idx(2),grid_idx(3)) - back_look*delta_l /) )
                    pi_bound = maxval( (/ 0d0, div_pol(grid_idx(1)-1,grid_idx(2),grid_idx(3)) - delta_div /) )

                    call grid_Cons_Equi( lgrid, l_bound, lpol(grid_idx(1)-1,grid_idx(2),grid_idx(3)) + delta_l )

                    if (d_bound >= d_ubound) then
                        call grid_Cons_Equi( dgrid, .99d0*d_bound, d_bound )
                    else
                        call grid_Cons_Equi( dgrid, d_bound, d_ubound )
                    endif

                    call grid_Cons_Equi( divgrid, pi_bound , div_pol(grid_idx(1)-1,grid_idx(2),grid_idx(3)) + delta_div )


                endif

                ! initialize current problem
                v_curr_max(6)   = 0d0
                l_curr_max(6)   = 0d0
                d_curr_max(6)   = 0d0
                div_curr_max(6) = 0d0
                s_curr_max(6)   = 0d0
                c_curr_max(6)   = 0d0
                a_curr_max(6)   = 0d0

                ! solve problem
                do ll=1,grid_len ! for each loan
                    do mm=1,grid_len ! for each deposit
                        do oo=1,grid_len   ! for each dividend

                            implied_a = 0d0
                            implied_c = 0d0

                            ! implied securities
                            implied_s = bank_ngrid(grid_idx(2),grid_idx(1)) + dgrid(mm) - divgrid(oo) -&
                                        lgrid(ll) - g(lgrid(ll),theta(grid_idx(2)))

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
                                            call linint_Grow( net_temp, bank_na(grid_idx(2)), bank_nb(grid_idx(2)),&
                                                                growth,bank_nlen-1, ixl, ixr, varphi)

                                            v_temp = ( varphi*v(ixl+1,grid_idx(2),ss) + &
                                                        (1d0-varphi)*v(ixr+1,grid_idx(2),ss) )

                                            RHS = RHS + prob_l(rr)*dbar_prob(grid_idx(2),grid_idx(3),ss)*v_temp
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
                                c_curr_max(6)   = 0d0
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
                ctemp(ii)    = c_curr_max( v_idx(1) )

                dtemp(ii)    = d_curr_max( v_idx(1) )
                atemp(ii)    = a_curr_max( v_idx(1) )

                divtemp(ii) = div_curr_max( v_idx(1) )

            enddo   ! state space loop

            ! MPI GATHER
            call MPI_GATHER(vtemp,proc_len,MPI_REAL8,  vtemp_all,proc_len,   MPI_REAL8,0,MPI_COMM_WORLD,ierr)
            call MPI_GATHER(ltemp,proc_len,MPI_REAL8,  ltemp_all,proc_len,   MPI_REAL8,0,MPI_COMM_WORLD,ierr)
            call MPI_GATHER(dtemp,proc_len,MPI_REAL8,  dtemp_all,proc_len,   MPI_REAL8,0,MPI_COMM_WORLD,ierr)
            call MPI_GATHER(stemp,proc_len,MPI_REAL8,  stemp_all,proc_len,   MPI_REAL8,0,MPI_COMM_WORLD,ierr)
            call MPI_GATHER(ctemp,proc_len,MPI_REAL8,  ctemp_all,proc_len,   MPI_REAL8,0,MPI_COMM_WORLD,ierr)
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
                cpol = transform_f(ctemp_all)
                div_pol = transform_f(divtemp_all)

                ! compute error
                error = sum( (v-vnew)**2d0 )/(bank_nlen*size(theta)*dbar_size)

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
            call MPI_BCAST(ctemp_all,  proc_len*nprocs,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
            call MPI_BCAST(divtemp_all,proc_len*nprocs,MPI_REAL8,0,MPI_COMM_WORLD,ierr)

            if ( my_id /= 0) then

                ! construct dimension-based policy functions
                v = transform_f(vtemp_all)
                lpol = transform_f(ltemp_all)
                dpol = transform_f(dtemp_all)
                spol = transform_f(stemp_all)
                apol = transform_f(atemp_all)
                cpol = transform_f(ctemp_all)
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



    subroutine get_dist_tran(Fcurr,lrec,srec,crec,drec,div_rec,arec,shock,price,Fp)

        implicit none

        ! input/output variables
        real*8, dimension(bank_nlen,size(theta),dbar_size), intent(in) :: Fcurr, lrec, srec, crec, drec, div_rec, arec
        real*8, dimension(2), intent(in) :: shock, price  ! (pstar, Ra) (dshock, ishock)
        real*8, dimension(bank_nlen,size(theta),dbar_size), intent(inout) :: Fp

        ! other variables
        real*8, dimension(bank_nlen,size(theta),dbar_size) :: G_ent
        integer :: ii, jj, kk, ll, mm, nn, qq
        real*8 :: net_temp
        integer, dimension(1) :: min_idx
        real*8 :: c_tilde, s_tilde, mstar, excess_cash, mstar_agg
        real*8 :: net_int, stock
        real*8, dimension(bank_nlen,size(theta),dbar_size) :: Ftemp

        integer :: il, ir
        real*8  :: phi


        ! compute default probabilities
        do ii=1,bank_nlen
            do jj=1,size(theta)
                do kk=1,dbar_size

                    if ( bank_ngrid(jj,ii) <= nstar(jj,kk)) then
                        default_prob(ii,jj,kk)     = 0d0
                        default_liq_prob(ii,jj,kk) = 0d0
                    else

                        default_prob(ii,jj,kk)     = 0d0
                        default_liq_prob(ii,jj,kk) = 0d0

                        excess_cash = crec(ii,jj,kk)

                        ! for each liquidity shock
                        do ll=1,size(delta)

                            ! if can't meet funding shock, liquidity default
                            if ( (excess_cash + price(1)*srec(ii,jj,kk)) < (delta(ll)+shock(1))*arec(ii,jj,kk) ) then

                                default_liq_prob(ii,jj,kk) = default_liq_prob(ii,jj,kk) + prob_d(ll)

                            else

                                ! compute cash and security liquidations
                                c_tilde =  minval( (/ (delta(ll)+shock(1))*arec(ii,jj,kk), excess_cash /) )

                                if ( c_tilde == excess_cash ) then  ! if not enough cash to settle
                                    ! have to liquidate securities
                                    s_tilde = ( (delta(ll)+shock(1))*arec(ii,jj,kk) - c_tilde )/price(1)
                                else
                                    ! if enough cash, security liquidations = 0
                                    s_tilde = 0d0
                                endif

                                ! for each return shock
                                do mm=1,size(Rl)
                                    ! for each capacity constraint shock
                                    do nn=1,dbar_size

                                        net_int = (Rl(mm)-shock(2)-1d0)*lrec(ii,jj,kk) + i_s*(srec(ii,jj,kk)-s_tilde) &
                                                         - (Rd-1d0)*drec(ii,jj,kk) - &
                                                         (price(2)-1d0)*(1d0-delta(ll)-shock(1))*arec(ii,jj,kk)

                                        stock = lrec(ii,jj,kk) + (srec(ii,jj,kk)-s_tilde) +&
                                                (crec(ii,jj,kk)-c_tilde) - drec(ii,jj,kk) - &
                                                (1d0-delta(ll)-shock(1))*arec(ii,jj,kk)

                                        if ( net_int > 0d0 ) then
                                            net_temp = (1d0-tax)*net_int + stock
                                        else
                                            net_temp = net_int + stock
                                        endif

                                        ! compute default probability
                                        if (net_temp <= nstar(jj,kk) ) then
                                            default_prob(ii,jj,kk) = default_prob(ii,jj,kk) + &
                                                            prob_l(mm)*dbar_prob(jj,kk,nn)*prob_d(ll)
                                        endif

                                    enddo
                                enddo


                            endif

                        enddo
                    endif
                enddo
            enddo
        enddo

        ! intialize
        Ftemp = 0d0
        G_ent = 0d0
        ! update next-period distribution
        do ii=1,bank_nlen
            do jj=1,size(theta)
                do kk=1,dbar_size

                    excess_cash = crec(ii,jj,kk)

                    ! for each delta shock
                    do ll=1,size(delta)

                        ! if no default, compute liquidations
                        if ( (excess_cash + price(1)*srec(ii,jj,kk)) >= (delta(ll)+shock(1))*arec(ii,jj,kk) ) then

                            c_tilde =  minval( (/ (delta(ll)+shock(1))*arec(ii,jj,kk), excess_cash /) )
                            if ( c_tilde == excess_cash ) then  ! if not enough cash to settle
                                s_tilde = ( (delta(ll)+shock(1))*arec(ii,jj,kk) - c_tilde )/price(1)
                            else
                                s_tilde = 0d0
                            endif

                            ! for each return shock
                            do mm=1,size(Rl)
                                ! for each capacity constraint
                                do nn=1,dbar_size


                                    net_int = (Rl(mm)-shock(2)-1d0)*lrec(ii,jj,kk) + i_s*(srec(ii,jj,kk)-s_tilde) &
                                                     - (Rd-1d0)*drec(ii,jj,kk) -  &
                                                     (price(2)-1d0)*(1d0-delta(ll)-shock(1))*arec(ii,jj,kk)

                                    stock = lrec(ii,jj,kk) + (srec(ii,jj,kk)-s_tilde) +&
                                            (crec(ii,jj,kk)-c_tilde) - drec(ii,jj,kk) - &
                                             (1d0-delta(ll)-shock(1))*arec(ii,jj,kk)

                                    if ( net_int > 0d0 ) then
                                        net_temp = (1d0-tax)*net_int + stock
                                    else
                                        net_temp = net_int + stock
                                    endif

                                    ! if no default, record new (networth,theta,capacity constraint)
                                    if (net_temp > nstar(jj,kk) ) then

                                        ! find nearest grid point
                                        if ( net_temp > bank_nb(jj) ) then
                                            !min_idx(1) = bank_nlen
                                            il = bank_nlen-1
                                            ir = bank_nlen
                                            phi = 0d0
                                        else
                                            !min_idx= minloc( abs( bank_ngrid(jj,:) - net_temp ) )
                                            call linint_Grow(net_temp, bank_na(jj), bank_nb(jj), growth,&
                                                                bank_nlen+1, il, ir, phi)

                                        endif

                                        Ftemp(il,jj,nn) = Ftemp(il,jj,nn) + phi*&
                                                       prob_l(mm)*dbar_prob(jj,kk,nn)*prob_d(ll)*Fcurr(ii,jj,kk)
                                        Ftemp(ir,jj,nn) = Ftemp(ir,jj,nn) + (1d0-phi)*&
                                                       prob_l(mm)*dbar_prob(jj,kk,nn)*prob_d(ll)*Fcurr(ii,jj,kk)

                                    endif
                                enddo
                            enddo



                        endif


                    enddo

                enddo
            enddo
        enddo

        ! create entrant distribution (with replacement)
        do ii = 1,bank_nlen
            do jj=1,size(theta)
                do kk=1,dbar_size
                    G_ent(ii,jj,kk) = (default_liq_prob(ii,jj,kk) + default_prob(ii,jj,kk))*Fcurr(ii,jj,kk)
                enddo
            enddo
        enddo

        mstar_agg = 0d0
        do jj=1,size(theta)
            if (sum(G_ent(:,jj,:)) == 0d0) then
                continue
            else
                G_ent(:,jj,:) = G_ent(:,jj,:)/sum(G_ent(:,jj,:))  ! normalize to 1
            endif

            !solve for mass of entrants
            mstar = 1d0 - sum(Ftemp(:,jj,:))

            mstar_agg = mstar_agg + prob_theta(jj)*mstar

            !next-period distribution, with entrants
            Ftemp(:,jj,:) = Ftemp(:,jj,:) + mstar*G_ent(:,jj,:)
        enddo

        Fp = Ftemp

    end subroutine



    subroutine get_dist()

        implicit none

        real*8, dimension(bank_nlen,size(theta),dbar_size) :: F0, F1, G_ent
        real*8, parameter :: tolerance = 1e-40
        real*8 :: error_temp
        integer :: ii, jj, kk, ll, mm, nn, qq
        integer, parameter :: it_max = 1000
        real*8 :: net_temp
        integer, dimension(1) :: min_idx
        real*8 :: c_tilde, s_tilde, mstar, excess_cash, mstar_agg
        real*8 :: net_int, stock

        integer :: il, ir
        real*8  :: phi

!        default_liq_prob = 0d0
        ! compute default probabilities
        do ii=1,bank_nlen
            do jj=1,size(theta)
                do kk=1,dbar_size

                    if ( bank_ngrid(jj,ii) <= nstar(jj,kk)) then
                        default_prob(ii,jj,kk)     = 0d0
                        default_liq_prob(ii,jj,kk) = 0d0
                    else

                        default_prob(ii,jj,kk)     = 0d0
                        default_liq_prob(ii,jj,kk) = 0d0

                        excess_cash = cpol(ii,jj,kk)

                        ! for each liquidity shock
                        do ll=1,size(delta)

                            ! if can't meet funding shock, liquidity default
                            if ( (excess_cash + pstar*spol(ii,jj,kk)) < delta(ll)*apol(ii,jj,kk) ) then

                                default_liq_prob(ii,jj,kk) = default_liq_prob(ii,jj,kk) + prob_d(ll)

                            else

                                ! compute cash and security liquidations
                                c_tilde =  minval( (/ delta(ll)*apol(ii,jj,kk), excess_cash /) )

                                if ( c_tilde == excess_cash ) then  ! if not enough cash to settle
                                    ! have to liquidate securities
                                    s_tilde = ( delta(ll)*apol(ii,jj,kk) - c_tilde )/pstar
                                else
                                    ! if enough cash, security liquidations = 0
                                    s_tilde = 0d0
                                endif

                                ! for each return shock
                                do mm=1,size(Rl)
                                    ! for each capacity constraint shock
                                    do nn=1,dbar_size

                                        net_int = (Rl(mm)-1d0)*lpol(ii,jj,kk) + i_s*(spol(ii,jj,kk)-s_tilde) &
                                                         - (Rd-1d0)*dpol(ii,jj,kk) - (Ra-1d0)*(1d0-delta(ll))*apol(ii,jj,kk)

                                        stock = lpol(ii,jj,kk) + (spol(ii,jj,kk)-s_tilde) +&
                                                (cpol(ii,jj,kk)-c_tilde) - dpol(ii,jj,kk) - (1d0-delta(ll))*apol(ii,jj,kk)

                                        if ( net_int > 0d0 ) then
                                            net_temp = (1d0-tax)*net_int + stock
                                        else
                                            net_temp = net_int + stock
                                        endif

                                        ! compute default probability
                                        if (net_temp <= nstar(jj,kk) ) then
                                            default_prob(ii,jj,kk) = default_prob(ii,jj,kk) + &
                                                            prob_l(mm)*dbar_prob(jj,kk,nn)*prob_d(ll)
                                        endif

                                    enddo
                                enddo


                            endif

                        enddo
                    endif
                enddo
            enddo
        enddo

        !Initialize problem
        F0 = 1d0/(bank_nlen*dbar_size)   ! still solving for marginal distributions (i.e. sum(F0(:,i,:)) = 1 )
        error_temp = tolerance + 1d0

        do qq = 1,it_max        !for each iteration in updating process

            F1 = 0d0

            do ii=1,bank_nlen
                do jj=1,size(theta)
                    do kk=1,dbar_size

                        excess_cash = cpol(ii,jj,kk)

                        ! for each delta shock
                        do ll=1,size(delta)

                            ! if no default, compute liquidations
                            if ( (excess_cash + pstar*spol(ii,jj,kk)) >= delta(ll)*apol(ii,jj,kk) ) then

                                c_tilde =  minval( (/ delta(ll)*apol(ii,jj,kk), excess_cash /) )
                                if ( c_tilde == excess_cash ) then  ! if not enough cash to settle
                                    s_tilde = ( delta(ll)*apol(ii,jj,kk) - c_tilde )/pstar
                                else
                                    s_tilde = 0d0
                                endif

                                ! for each return shock
                                do mm=1,size(Rl)
                                    ! for each capacity constraint
                                    do nn=1,dbar_size


                                        net_int = (Rl(mm)-1d0)*lpol(ii,jj,kk) + i_s*(spol(ii,jj,kk)-s_tilde) &
                                                         - (Rd-1d0)*dpol(ii,jj,kk) - (Ra-1d0)*(1d0-delta(ll))*apol(ii,jj,kk)

                                        stock = lpol(ii,jj,kk) + (spol(ii,jj,kk)-s_tilde) +&
                                                (cpol(ii,jj,kk)-c_tilde) - dpol(ii,jj,kk) - (1d0-delta(ll))*apol(ii,jj,kk)

                                        if ( net_int > 0d0 ) then
                                            net_temp = (1d0-tax)*net_int + stock
                                        else
                                            net_temp = net_int + stock
                                        endif

                                        ! if no default, record new (networth,theta,capacity constraint)
                                        if (net_temp > nstar(jj,kk) ) then

                                            ! find nearest grid point
                                            if ( net_temp > bank_nb(jj) ) then
                                                !min_idx(1) = bank_nlen
                                                il = bank_nlen-1
                                                ir = bank_nlen
                                                phi = 0d0
                                            else
                                                !min_idx= minloc( abs( bank_ngrid(jj,:) - net_temp ) )
                                                call linint_Grow(net_temp, bank_na(jj), bank_nb(jj), growth,&
                                                                    bank_nlen+1, il, ir, phi)

                                            endif
                                            !F1(min_idx(1),jj,nn) = F1(min_idx(1),jj,nn) + &
                                            !                prob_l(mm)*dbar_prob(jj,kk,nn)*prob_d(ll)*F0(ii,jj,kk)

                                            F1(il,jj,nn) = F1(il,jj,nn) + phi*prob_l(mm)*dbar_prob(jj,kk,nn)*prob_d(ll)*&
                                                                                            F0(ii,jj,kk)
                                            F1(ir,jj,nn) = F1(ir,jj,nn) + (1d0-phi)*prob_l(mm)*dbar_prob(jj,kk,nn)*&
                                                                                            prob_d(ll)*F0(ii,jj,kk)
                                        endif
                                    enddo
                                enddo



                            endif


                        enddo

                    enddo
                enddo
            enddo

            ! create entrant distribution (with replacement)
            do ii = 1,bank_nlen
                do jj=1,size(theta)
                    do kk=1,dbar_size
                        G_ent(ii,jj,kk) = (default_liq_prob(ii,jj,kk) + default_prob(ii,jj,kk))*F0(ii,jj,kk)
                    enddo
                enddo
            enddo

            mstar_agg = 0d0
            do jj=1,size(theta)
                if (sum(G_ent(:,jj,:)) == 0d0) then
                    continue
                else
                    G_ent(:,jj,:) = G_ent(:,jj,:)/sum(G_ent(:,jj,:))  ! normalize to 1
                endif

                !solve for mass of entrants
                mstar = 1d0 - sum(F1(:,jj,:))

                mstar_agg = mstar_agg + prob_theta(jj)*mstar

                !next-period distribution, with entrants
                F1(:,jj,:) = F1(:,jj,:) + mstar*G_ent(:,jj,:)
            enddo

            error_temp = sum(abs(F1-F0))/(bank_nlen*size(theta)*dbar_size )

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
        real*8, dimension(bank_nlen,size(theta),dbar_size) :: bshare
        real*8 :: agg_a
        real*8 :: obj1, obj2, obj3, obj4
        real*8 :: temp_net
        real*8 :: excess_cash
        real*8 :: net_int, stock

        agg_a = 0d0
        ! compute aggregate wholesale lending
        do ii=1,bank_nlen
            do jj=1,size(theta)
                do kk=1,dbar_size

                    agg_a = agg_a + prob_theta(jj)*F_stationary(ii,jj,kk)*apol(ii,jj,kk)

                enddo
            enddo
        enddo

        ! compute bank wholesalefunding shares
        do ii=1,bank_nlen
            do jj=1,size(theta)
                do kk=1,dbar_size

                    bshare(ii,jj,kk) = prob_theta(jj)*F_stationary(ii,jj,kk)*apol(ii,jj,kk)/agg_a

                enddo
            enddo
        enddo

        ! initialize
        obj1 = 0d0
        obj2 = 0d0
        obj3 = 0d0
        obj4 = 0d0

        ! for each bank type today (n_b,theta,dbar)
        do ii=1,bank_nlen
            do jj=1,size(theta)
                do kk=1,dbar_size

                    excess_cash = cpol(ii,jj,kk)

                    ! liquidity shocks
                    do ll=1,size(delta)

                        if ( delta(ll)*apol(ii,jj,kk) > excess_cash + pstar*spol(ii,jj,kk) ) then

                            ! computing object (1)
                            obj1 = obj1 + bshare(ii,jj,kk)*prob_d(ll)

                        else

                            ! computing object (2)
                            obj2 = obj2 + bshare(ii,jj,kk)*prob_d(ll)*delta(ll)

                            ! compute possible next-period outcomes
                            do mm=1,size(Rl)
                                do nn=1,dbar_size

                                    net_int = (Rl(mm)-1d0)*lpol(ii,jj,kk) + i_s*(spol(ii,jj,kk)-stilde(ii,jj,kk,ll)) &
                                                     - (Rd-1d0)*dpol(ii,jj,kk) - (Ra-1d0)*(1d0-delta(ll))*apol(ii,jj,kk)

                                    stock = lpol(ii,jj,kk) + (spol(ii,jj,kk)-stilde(ii,jj,kk,ll)) +&
                                            (cpol(ii,jj,kk)-ctilde(ii,jj,kk,ll)) - dpol(ii,jj,kk) - (1d0-delta(ll))*apol(ii,jj,kk)

                                    if ( net_int > 0d0 ) then
                                        temp_net = (1d0-tax)*net_int + stock
                                    else
                                        temp_net = net_int + stock
                                    endif

                                    ! computing objects (3) and (4)
                                    if ( temp_net < nstar(jj,kk)  ) then  ! insolvency default
                                        obj3 = obj3 + bshare(ii,jj,kk)*prob_d(ll)*(1d0-delta(ll))*prob_l(mm)*dbar_prob(jj,kk,nn)
                                    else
                                        obj4 = obj4 + bshare(ii,jj,kk)*prob_d(ll)*(1d0-delta(ll))*prob_l(mm)*dbar_prob(jj,kk,nn)
                                    endif


                                enddo
                            enddo


                        endif

                    enddo

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

    function implied_ra_rate_tran(Fcurr,lrec,srec,crec,drec,div_rec,arec,shock,price)

        implicit none

        ! input/output variables
        real*8, dimension(bank_nlen,size(theta),dbar_size), intent(in) :: Fcurr, lrec, srec, crec, drec, div_rec, arec
        real*8, dimension(2), intent(in) :: shock, price  ! (pstar, Ra) (dshock, ishock)
        real*8 :: implied_ra_rate_tran

        ! need to compute 4 objects
        !       (1) weighted liquidity default
        !       (2) (1-weighted liquidity default)*liquidity shock
        !       (3) (1-weighted liquidity default)*(1-liquidity shock)*weighted insolvency default
        !       (4) (1-weighted liquidity default)*(1-liquidity shock)*(1-weighted insolvency default)
        !
        !       then Ra = ( (1/beta) - (1+h)(1) - E[(2)] - (1+h)E[3] )/E[(4)]
        !

        integer :: ii,jj,kk,ll, mm, nn
        real*8, dimension(bank_nlen,size(theta),dbar_size) :: bshare
        real*8 :: agg_a
        real*8 :: obj1, obj2, obj3, obj4
        real*8 :: temp_net
        real*8 :: excess_cash
        real*8 :: net_int, stock

        agg_a = 0d0
        ! compute aggregate wholesale lending
        do ii=1,bank_nlen
            do jj=1,size(theta)
                do kk=1,dbar_size

                    agg_a = agg_a + prob_theta(jj)*Fcurr(ii,jj,kk)*arec(ii,jj,kk)

                enddo
            enddo
        enddo

        ! compute bank wholesalefunding shares
        do ii=1,bank_nlen
            do jj=1,size(theta)
                do kk=1,dbar_size

                    bshare(ii,jj,kk) = prob_theta(jj)*Fcurr(ii,jj,kk)*arec(ii,jj,kk)/agg_a

                enddo
            enddo
        enddo

        ! initialize
        obj1 = 0d0
        obj2 = 0d0
        obj3 = 0d0
        obj4 = 0d0

        ! for each bank type today (n_b,theta,dbar)
        do ii=1,bank_nlen
            do jj=1,size(theta)
                do kk=1,dbar_size

                    excess_cash = crec(ii,jj,kk)

                    ! liquidity shocks
                    do ll=1,size(delta)

                        if ( (delta(ll)+shock(1))*arec(ii,jj,kk) > excess_cash + price(1)*srec(ii,jj,kk) ) then

                            ! computing object (1)
                            obj1 = obj1 + bshare(ii,jj,kk)*prob_d(ll)

                        else

                            ! computing object (2)
                            obj2 = obj2 + bshare(ii,jj,kk)*prob_d(ll)*(delta(ll)+shock(1))

                            ! compute possible next-period outcomes
                            do mm=1,size(Rl)
                                do nn=1,dbar_size

                                    net_int = (Rl(mm)-shock(2)-1d0)*lrec(ii,jj,kk) + i_s*(srec(ii,jj,kk)-stilde(ii,jj,kk,ll)) &
                                                - (Rd-1d0)*drec(ii,jj,kk) - (price(2)-1d0)*(1d0-delta(ll)-shock(1))*arec(ii,jj,kk)

                                    stock = lrec(ii,jj,kk) + (srec(ii,jj,kk)-stilde(ii,jj,kk,ll)) +&
                                            (crec(ii,jj,kk)-ctilde(ii,jj,kk,ll)) - drec(ii,jj,kk) - &
                                            (1d0-delta(ll)-shock(1))*arec(ii,jj,kk)

                                    if ( net_int > 0d0 ) then
                                        temp_net = (1d0-tax)*net_int + stock
                                    else
                                        temp_net = net_int + stock
                                    endif

                                    ! computing objects (3) and (4)
                                    if ( temp_net < nstar(jj,kk)  ) then  ! insolvency default
                                        obj3 = obj3 + bshare(ii,jj,kk)*prob_d(ll)*(1d0-delta(ll)-shock(1))*prob_l(mm)*&
                                        dbar_prob(jj,kk,nn)
                                    else
                                        obj4 = obj4 + bshare(ii,jj,kk)*prob_d(ll)*(1d0-delta(ll)-shock(1))*prob_l(mm)*&
                                        dbar_prob(jj,kk,nn)
                                    endif


                                enddo
                            enddo


                        endif

                    enddo

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

        implied_ra_rate_tran = ( (1d0/beta_pre) - (1d0+hair)*obj1 - obj2 - (1d0+hair)*obj3 )/obj4

    end function


    function Aggregate_Liq_tran(Fcurr,srec,crec,arec,shock,price)

        implicit none

        ! input/output variables
        real*8, dimension(bank_nlen,size(theta),dbar_size), intent(in) :: Fcurr, srec, crec, arec
        real*8, intent(in) :: shock, price  !  (dshock)  (pstar)
        real*8, dimension(2) :: Aggregate_Liq_tran

        integer :: ii, jj, kk, ll, mm, nn
        integer, dimension(1) :: min_idx
        real*8 :: r_thresh
        real*8 :: def_diff, excess_cash, net_temp, liq_assets, stilde, ctilde
        real*8, dimension(2) :: temp_liq

        ! initialize
        temp_liq = 0d0

        ! for each (n_b,theta,dbar)
        do ii = 1,bank_nlen
            do jj =1,size(theta)
                do kk=1,dbar_size

                    ! for each funding shock
                    do ll=1,size(delta)

                        excess_cash  = crec(ii,jj,kk)

                        ! if no liquidity default
                        if ( excess_cash + price*srec(ii,jj,kk) >= (delta(ll)+shock)*arec(ii,jj,kk) ) then

                            ! compute liquidations
                            ctilde = minval( (/ (delta(ll)+shock)*arec(ii,jj,kk), excess_cash /) )

                            if ( ctilde == excess_cash ) then
                                stilde = ( (delta(ll)+shock)*arec(ii,jj,kk)-ctilde)/price
                            else
                                stilde = 0d0
                            endif

                            temp_liq(1) = temp_liq(1) + ctilde*prob_theta(jj)*prob_d(ll)*Fcurr(ii,jj,kk)
                            temp_liq(2) = temp_liq(2) + stilde*prob_theta(jj)*prob_d(ll)*Fcurr(ii,jj,kk)

                        else   ! liquidity default

                            temp_liq(1) = temp_liq(1) + excess_cash*prob_theta(jj)*prob_d(ll)*Fcurr(ii,jj,kk)
                            temp_liq(2) = temp_liq(2) + srec(ii,jj,kk)*prob_theta(jj)*prob_d(ll)*Fcurr(ii,jj,kk)

                        endif

                    enddo

                enddo
            enddo
        enddo

        Aggregate_Liq_tran = temp_liq

    end function


    subroutine grid_f(nprocs,proc_len,idx_grid)

        integer, intent(in) :: nprocs  ! number of processors
        integer, intent(in) :: proc_len  ! length of each processor grid

        integer, dimension(nprocs,proc_len,3), intent(out)  :: idx_grid  ! for each job/id, for each job grid point --> model state space: (nb,theta,dbar)

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
            do jj=1,size(theta)              ! for each theta lending type
                do kk=1,dbar_size            ! each deposit capacity constraint

                    ! record serial grid index
                    idx_grid(proc_idx,iterator,:) = (/ ii,jj,kk /)

                    if (iterator == proc_len) then  ! if iteration is last of processor length
                        iterator = 1                    ! reset iteration
                        proc_idx = proc_idx+1           ! update the processor id
                    else
                        iterator = iterator + 1
                    endif

                enddo
            enddo
        enddo

    end subroutine

    function transform_f(total_array)

        real*8, dimension(:), intent(in) :: total_array

        real*8, dimension(bank_nlen,size(theta),dbar_size) :: transform_f

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
            do jj=1,size(theta)             ! for each theta lending type
                do kk=1,dbar_size           ! for each deposit capacity constraint

                    transform_f(ii,jj,kk) = total_array(iterator)

                    iterator = iterator + 1

                enddo
            enddo
        enddo

    end function

    subroutine output()

        implicit none

        real*8, dimension(size(theta),3) :: port_shares
        real*8, dimension(size(theta),3) :: debt_shares
        real*8 :: total_assets, total_liabs, total_loan, total_cash, total_sec, total_debt, total_whole, total_net

        real*8, dimension(size(data_moms)) :: mod_moms
        real*8, dimension(5) :: al_sheet

        integer :: ii, jj, kk, ll, mm
        real*8 :: total_eq, total_rwa, elasticity

        real*8, dimension(6) :: size_corr
        real*8, dimension(8,bank_nlen*dbar_size*size(theta)) :: corr_objs  ! distribution weight, size, risk-weighted equity, leverage, liquidity, ins default, liq default, equity value

        real*8, dimension(7) :: size_corr_intra
        real*8, dimension(9,bank_nlen*dbar_size) :: corr_objs_intra

        real*8, dimension( 2, bank_nlen*dbar_size*size(theta) ) :: roe_stack
        real*8, dimension(2) :: roe_size_corr
        integer :: iterator
        real*8 :: roe1_num, roe1_den, roe2_num, roe2_den
        real*8 :: c_tilde, s_tilde, net_int, stock, net_temp
        real*8, dimension(bank_nlen,size(theta),dbar_size) :: roe_pol1, roe_pol2
        real*8 :: weighted_spread_num, weighted_spread_den
        real*8, dimension(5) :: holder
        real*8 :: hold
        real*8, dimension(size(theta),2) :: Agg_Def_Marg



        if (my_id==0) then

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
            write(*,*) 'deposit shares 1:          ', mod_moms(6),'        ',data_moms(6)
            write(*,*)
            write(*,*) 'deposit shares 2:          ', mod_moms(7),'        ',data_moms(7)
            write(*,*)
            write(*,*) 'deposit shares 3:          ', mod_moms(8),'        ',data_moms(8)
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
            write(*,*) 'cash:      ', 100d0*al_sheet(3)               ,'        ',5d0
            write(*,*)
            write(*,*) 'deposits:  ', 100d0*al_sheet(4)               ,'        ',75d0
            write(*,*)
            write(*,*) 'wholesale: ', 100d0*al_sheet(5)               ,'        ',20d0
            write(*,*)
            write(*,*) 'equity:    ', 100d0*(1d0 - al_sheet(4)-al_sheet(5)),'        ',5d0
            write(*,*)

            ! compute leverage ratios and size correlations with (i) equity ratios and (ii) liquidity ratios
            total_eq  = 0d0
            total_rwa = 0d0

            do ii=1,bank_nlen
                do jj=1,size(theta)
                    do kk=1,dbar_size

                        total_eq  = total_eq  + prob_theta(jj)*F_stationary(ii,jj,kk)*( lpol(ii,jj,kk) + cpol(ii,jj,kk) +&
                                    spol(ii,jj,kk) - (dpol(ii,jj,kk)+apol(ii,jj,kk)))

                        total_rwa = total_rwa + prob_theta(jj)*F_stationary(ii,jj,kk)*( lpol(ii,jj,kk) + &
                                            cpol(ii,jj,kk) + spol(ii,jj,kk) )

                    enddo
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
                do jj=1,size(theta)
                    do kk=1,dbar_size
                        do ll=1,size(Rl)

                            if ( lpol(ii,jj,kk) .ne. 0d0 ) then

                                weighted_spread_num = weighted_spread_num + prob_theta(jj)*F_stationary(ii,jj,kk)*&
                                                    Rl(ll)*lpol(ii,jj,kk)

                                weighted_spread_den = weighted_spread_den + prob_theta(jj)*F_stationary(ii,jj,kk)*&
                                                    (lpol(ii,jj,kk) + g(lpol(ii,jj,kk),theta(jj)) )



                            endif

                        enddo
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

            Agg_Def_Marg = Aggregate_Def_Marg()

            write(*,*)
            write(*,*) 'MARGINAL DEFAULT STATISICS'
            write(*,*)
            write(*,*) 'Bank 1 Insolvency Default Rate (bps):', 100d0*100d0*Agg_Def_Marg(1,1)
            write(*,*) 'Bank 1 Liquidity Default Rate (bps):', 100d0*100d0*Agg_Def_Marg(1,2)
            write(*,*)
            write(*,*) 'Bank 2 Insolvency Default Rate (bps):', 100d0*100d0*Agg_Def_Marg(2,1)
            write(*,*) 'Bank 2 Liquidity Default Rate (bps):', 100d0*100d0*Agg_Def_Marg(2,2)
            write(*,*)
            write(*,*) 'Bank 3 Insolvency Default Rate (bps):', 100d0*100d0*Agg_Def_Marg(3,1)
            write(*,*) 'Bank 3 Liquidity Default Rate (bps):', 100d0*100d0*Agg_Def_Marg(3,2)
            write(*,*)



            ! size correlations with (1) risk-weighted equity (2) leverage ratio and (3) liquidity ratio

            ! stack correlation objects
            corr_objs = 0d0
            iterator = 1

            do ii=1,bank_nlen
                do jj=1,size(theta)
                    do kk=1,dbar_size

                        ! distribution weights
                        corr_objs(1,iterator) =  prob_theta(jj)*F_stationary(ii,jj,kk)

                        ! size measure
                        corr_objs(2,iterator) =  lpol(ii,jj,kk)

                        ! risk-weighted equity ratios
                        if ( lpol(ii,jj,kk) .ne. 0d0 ) then
                            corr_objs(3,iterator) =  (lpol(ii,jj,kk) + cpol(ii,jj,kk) + spol(ii,jj,kk) - &
                                                            (dpol(ii,jj,kk)+apol(ii,jj,kk)))/lpol(ii,jj,kk)
                        else
                            corr_objs(3,iterator) = ebar
                        endif

                        ! leverage ratio
                        if ( ( lpol(ii,jj,kk) + cpol(ii,jj,kk) + spol(ii,jj,kk) ) .ne. 0d0 ) then
                            corr_objs(4,iterator) =  (lpol(ii,jj,kk) + cpol(ii,jj,kk) + spol(ii,jj,kk) - &
                                        (dpol(ii,jj,kk)+apol(ii,jj,kk)))/( lpol(ii,jj,kk) + cpol(ii,jj,kk) + spol(ii,jj,kk) )
                        else
                            corr_objs(4,iterator) =  ebar
                        endif

                        ! liquidity ratio
                        if ( lr_run*apol(ii,jj,kk) .ne. 0d0) then
                            corr_objs(5,iterator) =  (cpol(ii,jj,kk) + (1d0-lr_hair)*spol(ii,jj,kk))/( lr_run*apol(ii,jj,kk) )
                        else
                            corr_objs(5,iterator) =  phi_lr
                        endif

                        iterator = iterator + 1

                    enddo
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
            elasticity = fire_elastic(Agg_Liq(2))


            write(*,*) 'Liquidity Default Rate:', 100d0*100d0*Agg_Def(2)
            write(*,*)
            write(*,*)
            write(*,*) 'Elasticity of Firesold Securities:', elasticity
            write(*,*)

            ! size correlations with (1) insolvency default (2) liquidity default and (3) equity value

            ! stack correlation objects
            iterator = 1

            do ii=1,bank_nlen
                do jj=1,size(theta)
                    do kk=1,dbar_size

                        ! insolvency default
                        corr_objs(6,iterator) = default_prob(ii,jj,kk)

                        ! liquidity default
                        corr_objs(7,iterator) = default_liq_prob(ii,jj,kk)

                        iterator = iterator + 1

                    enddo
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
                do jj=1,size(theta)
                    do kk=1,dbar_size

                        ! for each liquidity shock
                        do ll=1,size(delta)

                            if ( pstar*spol(ii,jj,kk) + cpol(ii,jj,kk) >= delta(ll)*apol(ii,jj,kk)) then

                                c_tilde =  minval( (/ delta(ll)*apol(ii,jj,kk), cpol(ii,jj,kk) /) )

                                if ( c_tilde == cpol(ii,jj,kk) ) then  ! if not enough cash to settle
                                    ! have to liquidate securities
                                    s_tilde = ( delta(ll)*apol(ii,jj,kk) - cpol(ii,jj,kk) )/pstar
                                else
                                    s_tilde = 0d0
                                endif

                                ! for each loan shock
                                do mm=1,size(Rl)

                                    ! compute net interest
                                    net_int = (Rl(mm)-1d0)*lpol(ii,jj,kk) + i_s*(spol(ii,jj,kk)-s_tilde) - (Rd-1d0)*dpol(ii,jj,kk) &
                                                            - (Ra-1d0)*(1d0-delta(ll))*apol(ii,jj,kk)

                                    stock = lpol(ii,jj,kk) + (spol(ii,jj,kk)-s_tilde) + (cpol(ii,jj,kk)-c_tilde) - dpol(ii,jj,kk) &
                                                            - (1d0-delta(ll))*apol(ii,jj,kk)

                                    if ( net_int > 0d0) then
                                        net_temp = (1d0-tax)*net_int + stock
                                    else
                                        net_temp = net_int + stock
                                    endif

                                    if (net_temp >= 0d0) then

                                        ! record return on equity
                                        roe1_num = roe1_num + prob_l(mm)*prob_d(ll)*prob_theta(jj)*F_stationary(ii,jj,kk)*&
                                                    net_temp
                                        roe1_den = roe1_den + prob_l(mm)*prob_d(ll)*prob_theta(jj)*F_stationary(ii,jj,kk)*&
                                                    ( bank_ngrid(jj,ii) - div_pol(ii,jj,kk) )


                                        if ( ( bank_ngrid(jj,ii) - div_pol(ii,jj,kk) ) .ne. 0d0 ) then

                                            roe_pol1(ii,jj,kk) = roe_pol1(ii,jj,kk) + prob_l(mm)*prob_d(ll)*(net_temp/&
                                                                ( bank_ngrid(jj,ii) - div_pol(ii,jj,kk) ) -1d0)
                                        endif

                                        roe2_num = roe2_num + prob_l(mm)*prob_d(ll)*prob_theta(jj)*F_stationary(ii,jj,kk)*&
                                                    net_temp
                                        roe2_den = roe2_den + prob_l(mm)*prob_d(ll)*prob_theta(jj)*F_stationary(ii,jj,kk)*&
                                            ( lpol(ii,jj,kk) + cpol(ii,jj,kk) + spol(ii,jj,kk) - dpol(ii,jj,kk) - apol(ii,jj,kk) )

                                        if ( ( lpol(ii,jj,kk) + cpol(ii,jj,kk) + spol(ii,jj,kk) - &
                                                dpol(ii,jj,kk) - apol(ii,jj,kk) ) .ne. 0d0 ) then

                                            roe_pol2(ii,jj,kk) = roe_pol2(ii,jj,kk) + prob_l(mm)*prob_d(ll)*&
                                                (net_temp/( lpol(ii,jj,kk) + cpol(ii,jj,kk) + spol(ii,jj,kk) - &
                                                dpol(ii,jj,kk) - apol(ii,jj,kk) ) -1d0)

                                        endif



                                    endif

                                enddo
                            endif
                        enddo
                    enddo
                enddo
            enddo

            ! return on equity correlation with (1) size
            ! stack correlation objects
            iterator = 1

            do ii=1,bank_nlen
                do jj=1,size(theta)
                    do kk=1,dbar_size

                        ! insolvency default
                        roe_stack(1,iterator) = roe_pol1(ii,jj,kk)

                        ! liquidity default
                        roe_stack(2,iterator) = roe_pol2(ii,jj,kk)

                        iterator = iterator + 1

                    enddo
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
            do jj=1,size(theta)

                corr_objs_intra = 0d0
                iterator = 1

                do ii=1,bank_nlen
                    do kk=1,dbar_size

                        ! distribution weights
                        corr_objs_intra(1,iterator) =  F_stationary(ii,jj,kk)

                        ! size measure
                        corr_objs_intra(2,iterator) =  lpol(ii,jj,kk)

                        ! risk-weighted equity ratios
                        if ( lpol(ii,jj,kk) .ne. 0d0 ) then
                            corr_objs_intra(3,iterator) =  (lpol(ii,jj,kk) + cpol(ii,jj,kk) + spol(ii,jj,kk) - &
                                                            (dpol(ii,jj,kk)+apol(ii,jj,kk)))/lpol(ii,jj,kk)
                        else
                            corr_objs_intra(3,iterator) = ebar
                        endif

                        ! leverage ratio
                        if ( ( lpol(ii,jj,kk) + cpol(ii,jj,kk) + spol(ii,jj,kk) ) .ne. 0d0 ) then
                            corr_objs_intra(4,iterator) =  (lpol(ii,jj,kk) + cpol(ii,jj,kk) + spol(ii,jj,kk) - &
                                        (dpol(ii,jj,kk)+apol(ii,jj,kk)))/( lpol(ii,jj,kk) + cpol(ii,jj,kk) + spol(ii,jj,kk) )
                        else
                            corr_objs_intra(4,iterator) =  ebar
                        endif

                        ! liquidity ratio
                        if ( lr_run*apol(ii,jj,kk) .ne. 0d0) then
                            corr_objs_intra(5,iterator) =  (cpol(ii,jj,kk) + (1d0-lr_hair)*spol(ii,jj,kk))/( lr_run*apol(ii,jj,kk) )
                        else
                            corr_objs_intra(5,iterator) =  phi_lr
                        endif

                        ! insolvency default
                        corr_objs_intra(6,iterator) = default_prob(ii,jj,kk)

                        ! liquidity default
                        corr_objs_intra(7,iterator) = default_liq_prob(ii,jj,kk)

                        ! roe 1
                        corr_objs_intra(8,iterator) = roe_pol1(ii,jj,kk)

                        ! liquidity default
                        corr_objs_intra(9,iterator) = roe_pol2(ii,jj,kk)

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
                write(*,*) 'Bank Size Correlations, Type:', jj
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


            enddo

            ! compute intra-bank type correlation measures for net income and roe
            do jj=1,size(theta)

                corr_objs_intra = 0d0
                iterator = 1

                do ii=1,bank_nlen
                    do kk=1,dbar_size

                        ! distribution weights
                        corr_objs_intra(1,iterator) =  F_stationary(ii,jj,kk)

                        ! size measure
                        corr_objs_intra(2,iterator) =  bank_ngrid(jj,ii)

                        ! risk-weighted equity ratios
                        if ( lpol(ii,jj,kk) .ne. 0d0 ) then
                            corr_objs_intra(3,iterator) =  (lpol(ii,jj,kk) + cpol(ii,jj,kk) + spol(ii,jj,kk) - &
                                                            (dpol(ii,jj,kk)+apol(ii,jj,kk)))/lpol(ii,jj,kk)
                        else
                            corr_objs_intra(3,iterator) = ebar
                        endif

                        ! leverage ratio
                        if ( ( lpol(ii,jj,kk) + cpol(ii,jj,kk) + spol(ii,jj,kk) ) .ne. 0d0 ) then
                            corr_objs_intra(4,iterator) =  (lpol(ii,jj,kk) + cpol(ii,jj,kk) + spol(ii,jj,kk) - &
                                        (dpol(ii,jj,kk)+apol(ii,jj,kk)))/( lpol(ii,jj,kk) + cpol(ii,jj,kk) + spol(ii,jj,kk) )
                        else
                            corr_objs_intra(4,iterator) =  ebar
                        endif

                        ! liquidity ratio
                        if ( lr_run*apol(ii,jj,kk) .ne. 0d0) then
                            corr_objs_intra(5,iterator) =  (cpol(ii,jj,kk) + (1d0-lr_hair)*spol(ii,jj,kk))/( lr_run*apol(ii,jj,kk) )
                        else
                            corr_objs_intra(5,iterator) =  phi_lr
                        endif

                        ! insolvency default
                        corr_objs_intra(6,iterator) = default_prob(ii,jj,kk)

                        ! liquidity default
                        corr_objs_intra(7,iterator) = default_liq_prob(ii,jj,kk)

                        ! roe 1
                        corr_objs_intra(8,iterator) = roe_pol1(ii,jj,kk)

                        ! liquidity default
                        corr_objs_intra(9,iterator) = roe_pol2(ii,jj,kk)

                        iterator = iterator + 1

                    enddo
                enddo

                ! record
                size_corr_intra(1) = weighted_corr( corr_objs_intra(1,:), corr_objs_intra(2,:), corr_objs_intra(3,:) )
                size_corr_intra(2) = weighted_corr( corr_objs_intra(1,:), corr_objs_intra(2,:), corr_objs_intra(4,:) )
                size_corr_intra(3) = weighted_corr( corr_objs_intra(1,:), corr_objs_intra(2,:), corr_objs_intra(5,:) )
                size_corr_intra(4) = weighted_corr( corr_objs_intra(1,:), corr_objs_intra(2,:), corr_objs_intra(6,:) )
                size_corr_intra(5) = weighted_corr( corr_objs_intra(1,:), corr_objs_intra(2,:), corr_objs_intra(7,:) )
                size_corr_intra(6) = weighted_corr( corr_objs_intra(1,:), corr_objs_intra(2,:), corr_objs_intra(8,:) )
                size_corr_intra(7) = weighted_corr( corr_objs_intra(1,:), corr_objs_intra(2,:), corr_objs_intra(9,:) )

                write(*,*)
                write(*,*) 'Bank Net Income Correlations, Type:', jj
                write(*,*) 'risk-weighted equity ratio:', size_corr_intra(1)
                write(*,*) 'leverage ratio:', size_corr_intra(2)
                write(*,*) 'liquidity ratio:', size_corr_intra(3)
                write(*,*) 'insolvency default:', size_corr_intra(4)
                write(*,*) 'liquidity default:', size_corr_intra(5)
                write(*,*) 'ROE 1:', size_corr_intra(6)
                write(*,*) 'ROE 2:', size_corr_intra(7)


            enddo

            ! bank type portfolio shares
            port_shares = 0d0

            do ii =1,size(theta)

                total_assets = 0d0
                total_liabs  = 0d0
                total_loan   = 0d0
                total_cash   = 0d0
                total_sec    = 0d0
                total_debt   = 0d0
                total_whole  = 0d0
                total_net    = 0d0

                do jj=1,dbar_size
                    do kk=1,bank_nlen

                        total_assets = total_assets + F_stationary(kk,ii,jj)*( lpol(kk,ii,jj) + &
                        pc*cpol(kk,ii,jj) + ps*spol(kk,ii,jj) )
                        total_liabs  = total_liabs  + F_stationary(kk,ii,jj)*( dpol(kk,ii,jj) + &
                        apol(kk,ii,jj) + bank_ngrid(ii,kk) )

                        total_loan = total_loan + F_stationary(kk,ii,jj)*lpol(kk,ii,jj)
                        total_cash = total_cash + F_stationary(kk,ii,jj)*pc*cpol(kk,ii,jj)
                        total_sec  = total_sec  + F_stationary(kk,ii,jj)*ps*spol(kk,ii,jj)

                        total_debt  = total_debt  + F_stationary(kk,ii,jj)*dpol(kk,ii,jj)
                        total_whole = total_whole + F_stationary(kk,ii,jj)*apol(kk,ii,jj)
                        total_net   = total_net   + F_stationary(kk,ii,jj)*bank_ngrid(ii,kk)
                    enddo
                enddo

                port_shares(ii,1) = total_loan/total_assets
                port_shares(ii,2) = total_cash/total_assets
                port_shares(ii,3) = total_sec/total_assets

                debt_shares(ii,1) = total_debt/total_liabs
                debt_shares(ii,2) = total_whole/total_liabs
                debt_shares(ii,3) = total_net/total_liabs
            enddo

            write(*,*)
            write(*,*) 'BANK-LEVEL STATISTICS'
            write(*,*)
            write(*,*) 'THETA = 1 BANKS'
            write(*,*) 'Loan portfolio share (%):', 100d0*port_shares(1,1)
            write(*,*) 'Cash portfolio share (%):',100d0*port_shares(1,2)
            write(*,*) 'Security portfolio share (%):', 100d0*port_shares(1,3)
            write(*,*) 'Deposit liability share (%):', 100d0*debt_shares(1,1)
            write(*,*) 'Wholesale liability share (%):', 100d0*debt_shares(1,2)
            write(*,*) 'Networth liabiloty share (%):', 100d0*debt_shares(1,3)
            write(*,*)
            write(*,*) 'THETA = 2 BANKS'
            write(*,*) 'Loan portfolio share (%):', 100d0*port_shares(2,1)
            write(*,*) 'Cash portfolio share (%):',100d0*port_shares(2,2)
            write(*,*) 'Security portfolio share (%):', 100d0*port_shares(2,3)
            write(*,*) 'Deposit liability share (%):', 100d0*debt_shares(2,1)
            write(*,*) 'Wholesale liability share (%):', 100d0*debt_shares(2,2)
            write(*,*) 'Networth liabiloty share (%):', 100d0*debt_shares(2,3)
            write(*,*)
            write(*,*)
            write(*,*) 'THETA = 3 BANKS'
            write(*,*) 'Loan portfolio share (%):', 100d0*port_shares(3,1)
            write(*,*) 'Cash portfolio share (%):',100d0*port_shares(3,2)
            write(*,*) 'Security portfolio share (%):', 100d0*port_shares(3,3)
            write(*,*) 'Deposit liability share (%):', 100d0*debt_shares(3,1)
            write(*,*) 'Wholesale liability share (%):', 100d0*debt_shares(3,2)
            write(*,*) 'Networth liabiloty share (%):', 100d0*debt_shares(3,3)
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
            Agg_Cash  = holder(2)
            Agg_Sec   = holder(3)
            Agg_Dep   = holder(4)
            Agg_Whole = holder(5)


            write(*,*) 'Aggregates for Comparison to Dodd-Frank Outcomes'
            write(*,*)
            ! agg lending
            write(*,*) 'Aggregate Lending:', Agg_Loan

            ! agg balance sheet (l + s + c)
            write(*,*) 'Aggregate Bank Balance Sheet Size:', Agg_Loan + Agg_Sec + Agg_Cash
            write(*,*) 'Aggregate Balance Sheet Liquidity (%):', 100d0*(Agg_Sec + Agg_Cash)/( Agg_Loan + Agg_Sec + Agg_Cash )
            write(*,*) 'Aggregate Wholesale Funding:', Agg_Whole

            ! agg consumption
            write(*,*) 'Household Consumption:', (Rd-1d0)*Agg_Dep + Agg_Div(1) + Agg_Div(2) - Agg_Gov - Agg_DI + Agg_Corp
            write(*,*)

        endif

    end subroutine

    subroutine output_tran()

        implicit none


        integer :: ii,jj,kk,ll,mm,nn,tt

        real*8, dimension(T) :: time, ins_def_time, liq_def_time, agg_loan_time, agg_whole_time, rwcr_time, lev_time,&
                                lr_time, agg_liq_time

        real*8, dimension(2)  :: temp_def
        integer, dimension(1) :: min_idx

        real*8 :: temp_l, temp_a, net_temp, c_tilde, s_tilde, mstar, excess_cash, mstar_agg, net_int, stock, temp_liq, &
                  ctilde, stilde, total_eq, total_rwa, liq_num, liq_den

        integer :: il, ir
        real*8  :: phi

        do ii=1,T
            time(ii) = ii
        enddo


        !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
        !                                                       !
        !   insolvency and liquidity default rates over time    !
        !                                                       !
        !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!

        ins_def_time(1) = 100d0*Agg_Def_stationary(1)
        liq_def_time(1) = 100d0*Agg_Def_stationary(2)

        do tt=1,T-1 ! for each time period

            ! for each bank type
            do ii=1,bank_nlen
                do jj=1,size(theta)
                    do kk=1,dbar_size

                        if ( bank_ngrid(jj,ii) <= nstar(jj,kk)) then
                            default_prob(ii,jj,kk)     = 0d0
                            default_liq_prob(ii,jj,kk) = 0d0
                        else

                            default_prob(ii,jj,kk)     = 0d0
                            default_liq_prob(ii,jj,kk) = 0d0

                            excess_cash = ctran(tt,ii,jj,kk)

                            ! for each liquidity shock
                            do ll=1,size(delta)

                                ! if can't meet funding shock, liquidity default
                                if ( (excess_cash + pguess(tt)*stran(tt,ii,jj,kk)) < &
                                                            (delta(ll)+dshock(tt+1))*atran(tt,ii,jj,kk) ) then

                                    default_liq_prob(ii,jj,kk) = default_liq_prob(ii,jj,kk) + prob_d(ll)

                                else

                                    ! compute cash and security liquidations
                                    c_tilde =  minval( (/ (delta(ll)+dshock(tt+1))*atran(tt,ii,jj,kk), excess_cash /) )

                                    if ( c_tilde == excess_cash ) then  ! if not enough cash to settle
                                        ! have to liquidate securities
                                        s_tilde = ( (delta(ll)+dshock(tt+1))*atran(tt,ii,jj,kk) - c_tilde )/pguess(tt)
                                    else
                                        ! if enough cash, security liquidations = 0
                                        s_tilde = 0d0
                                    endif

                                    ! for each return shock
                                    do mm=1,size(Rl)
                                        ! for each capacity constraint shock
                                        do nn=1,dbar_size

                                            net_int = (Rl(mm)-ishock(tt+1)-1d0)*ltran(tt,ii,jj,kk) + &
                                                            i_s*(stran(tt,ii,jj,kk)-s_tilde) - (Rd-1d0)*dtran(tt,ii,jj,kk) - &
                                                             (rguess(tt)-1d0)*(1d0-delta(ll)-dshock(tt+1))*atran(tt,ii,jj,kk)

                                            stock = ltran(tt,ii,jj,kk) + (stran(tt,ii,jj,kk)-s_tilde) +&
                                                    (ctran(tt,ii,jj,kk)-c_tilde) - dtran(tt,ii,jj,kk) - &
                                                    (1d0-delta(ll)-dshock(tt+1))*atran(tt,ii,jj,kk)

                                            if ( net_int > 0d0 ) then
                                                net_temp = (1d0-tax)*net_int + stock
                                            else
                                                net_temp = net_int + stock
                                            endif

                                            ! compute default probability
                                            if (net_temp <= nstar(jj,kk) ) then
                                                default_prob(ii,jj,kk) = default_prob(ii,jj,kk) + &
                                                                prob_l(mm)*dbar_prob(jj,kk,nn)*prob_d(ll)
                                            endif

                                        enddo
                                    enddo


                                endif

                            enddo
                        endif
                    enddo
                enddo
            enddo

            ! compute total default rates
            temp_def = 0d0

            do ii=1,bank_nlen
                do jj=1,size(theta)
                    do kk=1,dbar_size

                        temp_def(1) = temp_def(1) + prob_theta(jj)*Ftran(tt,ii,jj,kk)*default_prob(ii,jj,kk)
                        temp_def(2) = temp_def(2) + prob_theta(jj)*Ftran(tt,ii,jj,kk)*default_liq_prob(ii,jj,kk)

                    enddo
                enddo
            enddo

            ins_def_time(tt+1) = 100d0*temp_def(1)
            liq_def_time(tt+1) = 100d0*temp_def(2)

        enddo


        !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
        !                                              !
        !   Aggregate Lending and Wholesale Lending    !
        !                                              !
        !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
        do tt=1,T
            temp_l = 0d0
            temp_a = 0d0

            do ii=1,bank_nlen
                do jj=1,size(theta)
                    do kk=1,dbar_size

                        temp_l = temp_l + prob_theta(jj)*Ftran(tt,ii,jj,kk)*ltran(tt,ii,jj,kk)
                        temp_a = temp_a + prob_theta(jj)*Ftran(tt,ii,jj,kk)*atran(tt,ii,jj,kk)

                    enddo
                enddo
            enddo

            agg_loan_time(tt)  = temp_l
            agg_whole_time(tt) = temp_a
        enddo

        ! normalize (% deviations from steady state)
        agg_loan_time  = 100d0*(agg_loan_time -agg_loan_time(1) )/agg_loan_time(1)
        agg_whole_time = 100d0*(agg_whole_time-agg_whole_time(1))/agg_whole_time(1)


        !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
        !                                    !
        !   Secondary Market Liquidations    !
        !                                    !
        !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!

        ! funding sohcks delta' occur in the intra-period Settlement Stage.
        !    I assume shocks at date t+1 hit during the Settlement Stage of period t

        ! fill in stationary outcome at date t=1
        temp_liq = 0d0

        do ii = 1,bank_nlen
            do jj =1,size(theta)
                do kk=1,dbar_size

                    ! for each funding shock
                    do ll=1,size(delta)

                        excess_cash  = ctran(1,ii,jj,kk)

                        ! if no liquidity default
                        if ( excess_cash + pguess(1)*stran(1,ii,jj,kk) >= (delta(ll)*atran(1,ii,jj,kk) ) ) then

                            ! compute liquidations
                            ctilde = minval( (/ delta(ll)*atran(1,ii,jj,kk), excess_cash /) )

                            if ( ctilde == excess_cash ) then
                                stilde = (delta(ll)*atran(1,ii,jj,kk)-ctilde)/pguess(1)
                            else
                                stilde = 0d0
                            endif

                            temp_liq = temp_liq + stilde*prob_theta(jj)*prob_d(ll)*Ftran(1,ii,jj,kk)

                        else   ! liquidity default

                            temp_liq = temp_liq + stran(1,ii,jj,kk)*prob_theta(jj)*prob_d(ll)*Ftran(1,ii,jj,kk)

                        endif

                    enddo

                enddo
            enddo
        enddo
        agg_liq_time(1) = temp_liq


        ! fill in outcomes at date t=2,...,T
        do tt=1,T-1
            ! initialize
            temp_liq = 0d0

            ! for each (n_b,theta,dbar)
            do ii = 1,bank_nlen
                do jj =1,size(theta)
                    do kk=1,dbar_size

                        ! for each funding shock
                        do ll=1,size(delta)

                            excess_cash  = ctran(tt,ii,jj,kk)

                            ! if no liquidity default
                            if ( excess_cash + pguess(tt)*stran(tt,ii,jj,kk) >= &
                                            (delta(ll)+dshock(tt+1))*atran(tt,ii,jj,kk) ) then

                                ! compute liquidations
                                ctilde = minval( (/ (delta(ll)+dshock(tt+1))*atran(tt,ii,jj,kk), excess_cash /) )

                                if ( ctilde == excess_cash ) then
                                    stilde = ((delta(ll)+dshock(tt+1))*atran(tt,ii,jj,kk)-ctilde)/pguess(tt)
                                else
                                    stilde = 0d0
                                endif

                                temp_liq = temp_liq + stilde*prob_theta(jj)*prob_d(ll)*Ftran(tt,ii,jj,kk)

                            else   ! liquidity default

                                temp_liq = temp_liq + stran(tt,ii,jj,kk)*prob_theta(jj)*prob_d(ll)*Ftran(tt,ii,jj,kk)

                            endif

                        enddo

                    enddo
                enddo
            enddo

            agg_liq_time(tt+1) = temp_liq

        enddo

        ! normalize
        agg_liq_time  = 100d0*(agg_liq_time -agg_liq_time(1) )/agg_liq_time(1)


        !~~~~~~~~~~~~~~~~~~~~~!
        !                     !
        !   Leverage Ratio    !
        !                     !
        !~~~~~~~~~~~~~~~~~~~~~!
        do tt=1,T

            total_eq  = 0d0
            total_rwa = 0d0

            do ii=1,bank_nlen
                do jj=1,size(theta)
                    do kk=1,dbar_size

                        total_eq  = total_eq  + prob_theta(jj)*Ftran(tt,ii,jj,kk)*( ltran(tt,ii,jj,kk) + ctran(tt,ii,jj,kk) +&
                                    stran(tt,ii,jj,kk) - (dtran(tt,ii,jj,kk)+atran(tt,ii,jj,kk)))

                        total_rwa = total_rwa + prob_theta(jj)*Ftran(tt,ii,jj,kk)*(ltran(tt,ii,jj,kk)+&
                                                stran(tt,ii,jj,kk)+ctran(tt,ii,jj,kk))

                    enddo
                enddo
            enddo

            lev_time(tt) = 100d0*total_eq/total_rwa

        enddo
        !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
        !                                 !
        !   Risk-weighted Equity Ratio    !
        !                                 !
        !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
        do tt=1,T

            total_eq  = 0d0
            total_rwa = 0d0

            do ii=1,bank_nlen
                do jj=1,size(theta)
                    do kk=1,dbar_size

                        total_eq  = total_eq  + prob_theta(jj)*Ftran(tt,ii,jj,kk)*( ltran(tt,ii,jj,kk) + ctran(tt,ii,jj,kk) +&
                                    stran(tt,ii,jj,kk) - (dtran(tt,ii,jj,kk)+atran(tt,ii,jj,kk)))

                        total_rwa = total_rwa + prob_theta(jj)*Ftran(tt,ii,jj,kk)*ltran(tt,ii,jj,kk)

                    enddo
                enddo
            enddo

            rwcr_time(tt) = 100d0*total_eq/total_rwa

        enddo

        !~~~~~~~~~~~~~~~~~~~~~~!
        !                      !
        !   Liquidity Ratio    !
        !                      !
        !~~~~~~~~~~~~~~~~~~~~~~!
        do tt=1,T

            liq_num = 0d0
            liq_den = 0d0

            do ii=1,bank_nlen
                do jj=1,size(theta)
                    do kk=1,dbar_size

                        if (atran(tt,ii,jj,kk) >0d0) then

                            liq_num = liq_num + prob_theta(jj)*Ftran(tt,ii,jj,kk)*(ctran(tt,ii,jj,kk) +&
                                                (1d0-lr_hair)*stran(tt,ii,jj,kk))

                            liq_den = liq_den + prob_theta(jj)*Ftran(tt,ii,jj,kk)*( lr_run*atran(tt,ii,jj,kk) )
                        endif

                    enddo
                enddo
            enddo

            lr_time(tt) = 100d0*liq_num/liq_den

        enddo




        !~~~~~~~~~~~~~~~~~~~!
        !                   !
        !   Graph Output    !
        !                   !
        !~~~~~~~~~~~~~~~~~~~!

        if (my_id ==0 ) then

            call plot(time,4d0*ins_def_time,legend='insolvency default',markersize=3d0)
            call execplot(title='Insolvency Default Rates (Annualized %)',filename='ins_def',filetype='png',output='ins_def')

            call plot(time,4d0*liq_def_time,legend='liquidity default',markersize=3d0)
            call execplot(title='Liquidity Default Rates (Annualized %)',filename='liq_def',filetype='png',output='liq_def')

            call plot(time,dshock,legend='liquidity shock',markersize=3d0)
            call plot(time,ishock,legend='loan shock',markersize=3d0)
            call execplot(title='Aggregate Shocks',filename='shocks',filetype='png',output='shocks')

            call plot(time,rguess,markersize=3d0)
            call execplot(title='Wholesale Funding Rate',filename='ra',filetype='png',output='ra')

            call plot(time,pguess,markersize=3d0)
            call execplot(title='Firesale Price',filename='pstar',filetype='png',output='pstar')

            call plot(time,agg_loan_time,markersize=3d0)
            call execplot(title='Aggregate Lending (% Deviation from Steady State)',filename='loan',filetype='png',output='loan')

            call plot(time,agg_whole_time,markersize=3d0)
            call execplot(title='Aggregate Wholesale Lending (% Deviation from Steady State)',filename='whole',&
                                    filetype='png',output='whole')

            call plot(time,agg_liq_time,markersize=3d0)
            call execplot(title='Aggregate Secondary Market Liquidations (% Deviation from Steady State)',filename='liq',&
                                filetype='png',output='liq')

            call plot(time,rwcr_time,markersize=3d0)
            call execplot(title='Risk-Weighted Equity Ratio',filename='rwcr',filetype='png',output='rwcr')

            call plot(time,lev_time,markersize=3d0)
            call execplot(title='Leverage Ratio',filename='lev',filetype='png',output='lev')

            call plot(time,lr_time,markersize=3d0)
            call execplot(title='Liquidity Ratio Ratio',filename='lr',filetype='png',output='lr')

            ! write/print these quantities to csv
            open(40,file='loan_tran.csv',access = 'append')
            open(41,file='whole_tran.csv',access = 'append')
            open(42,file='secondary_tran.csv',access = 'append')

            open(43,file='ishock.csv',access = 'append')
            open(44,file='dshock.csv',access = 'append')

            open(45,file='rwcr_tran.csv',access = 'append')
            open(46,file='lev_tran.csv',access = 'append')
            open(47,file='liq_tran.csv',access = 'append')

            open(48,file='ra_tran.csv',access = 'append')
            open(49,file='pstar_tran.csv',access = 'append')

            open(50,file='ins_def_tran.csv',access = 'append')
            open(51,file='liq_def_tran.csv',access = 'append')


            write(40,*) agg_loan_time
            write(41,*) agg_whole_time
            write(42,*) agg_liq_time

            write(43,*) ishock
            write(44,*) dshock

            write(45,*) rwcr_time
            write(46,*) lev_time
            write(47,*) lr_time

            write(48,*) rguess
            write(49,*) pguess

            write(50,*) 4d0*ins_def_time
            write(51,*) 4d0*liq_def_time


            close(40)
            close(41)
            close(42)
            close(43)
            close(44)
            close(45)
            close(46)
            close(47)
            close(48)
            close(49)
            close(50)
            close(51)

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
