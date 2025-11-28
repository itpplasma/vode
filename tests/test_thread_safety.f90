!> @file test_thread_safety.f90
!> @brief Thread-safety test for VODE ODE solver with OpenMP
!>
!> This test verifies that the VODE solver produces consistent results when
!> called concurrently from multiple OpenMP threads. Each thread independently
!> solves the same simple exponential decay ODE:
!>
!>   dy/dt = -y
!>
!> With initial conditions y = [1.0, 2.0] and integrating from t=0 to t=1,
!> the exact solution is y(1) = exp(-1) ≈ 0.367879...
!>
!> If the solver is thread-safe, all threads should produce identical results
!> within numerical tolerance. Without proper threadprivate declarations,
!> threads would corrupt each other's solver state leading to different
!> results or runtime errors.

program test_thread_safety
    use omp_lib
    use dvode_f90_m
    implicit none

    integer, parameter :: NUM_THREADS = 8
    integer, parameter :: NUM_EQUATIONS = 2
    real(kind=8), parameter :: EXPECTED_RESULT = 0.36787944117144232d0  ! exp(-1)
    real(kind=8), parameter :: TOLERANCE = 1.0d-8

    real(kind=8) :: results(NUM_THREADS)
    real(kind=8) :: max_deviation
    integer :: i, num_failures
    logical :: test_passed

    ! Set number of threads for the test
    call omp_set_num_threads(NUM_THREADS)

    write(*,'(A)') "VODE Thread-Safety Test"
    write(*,'(A)') "========================"
    write(*,'(A,I0,A)') "Running with ", NUM_THREADS, " OpenMP threads"
    write(*,'(A)') ""

    ! Initialize results to invalid values to detect failures
    results = -1.0d0

    ! Run VODE in parallel from multiple threads
    !$omp parallel
    block
        integer :: thread_id, istate, itask
        real(kind=8) :: y(NUM_EQUATIONS), t, tout
        real(kind=8) :: atol(NUM_EQUATIONS), rtol
        type(vode_opts) :: opts

        thread_id = omp_get_thread_num() + 1

        ! Set up the ODE problem: dy/dt = -y (exponential decay)
        y = [1.0d0, 2.0d0]
        t = 0.0d0
        tout = 1.0d0

        ! Error tolerances
        rtol = 1.0d-9
        atol = 1.0d-10

        ! Integration control
        itask = 1   ! Normal computation
        istate = 1  ! First call to solver

        ! Initialize solver options
        opts = set_normal_opts(abserr_vector=atol, relerr=rtol)

        ! Solve the ODE
        call dvode_f90(rhs_exponential_decay, NUM_EQUATIONS, y, t, tout, &
                       itask, istate, opts)

        ! Store result for comparison
        results(thread_id) = y(1)
    end block
    !$omp end parallel

    ! Analyze results
    write(*,'(A)') "Results from each thread:"
    do i = 1, NUM_THREADS
        write(*,'(A,I0,A,ES23.16)') "  Thread ", i, ": y(1) = ", results(i)
    end do
    write(*,'(A)') ""

    ! Check that all results match the expected value
    num_failures = 0
    max_deviation = 0.0d0

    do i = 1, NUM_THREADS
        max_deviation = max(max_deviation, abs(results(i) - EXPECTED_RESULT))
        if (abs(results(i) - EXPECTED_RESULT) > TOLERANCE) then
            num_failures = num_failures + 1
            write(*,'(A,I0,A,ES10.3,A)') "  WARNING: Thread ", i, &
                " result deviates by ", abs(results(i) - EXPECTED_RESULT), &
                " from expected value"
        end if
    end do

    ! Check that all threads got the same result
    do i = 2, NUM_THREADS
        if (abs(results(i) - results(1)) > 1.0d-14) then
            write(*,'(A,I0,A)') "  ERROR: Thread ", i, &
                " produced different result than thread 1"
            num_failures = num_failures + 1
        end if
    end do

    write(*,'(A,ES23.16)') "Expected value (exp(-1)):  ", EXPECTED_RESULT
    write(*,'(A,ES10.3)')  "Maximum deviation:         ", max_deviation
    write(*,'(A)') ""

    ! Report test outcome
    test_passed = (num_failures == 0)

    if (test_passed) then
        write(*,'(A)') "TEST PASSED: All threads produced consistent, correct results"
        write(*,'(A)') "             Thread-safety of VODE solver verified."
    else
        write(*,'(A,I0,A)') "TEST FAILED: ", num_failures, " check(s) failed"
        write(*,'(A)') "             VODE solver may not be thread-safe!"
        error stop 1
    end if

contains

    !> Right-hand side function for exponential decay ODE
    !> dy/dt = -y
    subroutine rhs_exponential_decay(neq, t, y, ydot)
        integer, intent(in) :: neq
        real(kind=8), intent(in) :: t
        real(kind=8), intent(in) :: y(neq)
        real(kind=8), intent(out) :: ydot(neq)

        ! Suppress unused argument warning
        if (.false.) print *, t

        ! Simple exponential decay: dy/dt = -y
        ydot = -y
    end subroutine rhs_exponential_decay

end program test_thread_safety
