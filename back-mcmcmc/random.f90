module random
    implicit none
    integer, parameter :: MAXDEPTH = 2**20
    integer, parameter, private :: P = 250*2**2 !magic number
    integer, parameter, private :: Q = 103*2**2 !magic number
    integer, parameter, private:: MASK = MAXDEPTH - 1
    integer, parameter, private :: IMAX = 2147483647 ! == 2**31 - 1
    real(8), parameter, private :: DMAX = dble(IMAX) + 1d0
    integer, private :: rndstep
    integer, private :: DRN(0:MAXDEPTH-1)
    integer :: ir(MAXDEPTH)
    real(8) :: ur(MAXDEPTH)
    
    private :: coin, initR, setDRN

contains

    integer function coin(mrand)
        integer mrand
        
        mrand = mrand * 16807
        if(mrand < 0)then
            coin = 0
        else
            coin = 1
        end if
        
    end function coin
    
    integer function initR(mrand)
        integer mrand, R, i
        
        R = coin(mrand)
        do i = 1, 31
            R = ishft(R, 1)
            R = R + coin(mrand)
        end do
        initR = R
        
    end function initR
    
    subroutine setDRN(iseed)
        integer, intent(in) :: iseed
        integer :: mrand, k
        
        mrand = iseed
        do k = 0, MAXDEPTH - 1
            DRN(k) = initR(mrand)
        end do
    
    end subroutine setDRN
    
    subroutine rfr(nsamp)
        integer, intent(in) :: nsamp
        integer :: i, j, k, n
        
        do n = rndstep + 1, rndstep + nsamp
            i = IAND(n, MASK)
            j = IAND(n - P, MASK)
            k = IAND(n - Q, MASK)
            DRN(i) = IEOR(DRN(j), DRN(k))
        end do
        
        do i = 1, nsamp
            n = IAND(i + rndstep, MASK)
            ir(i) = IAND(DRN(n), IMAX)
            ur(i) = ir(i) / DMAX
        end do
        
        rndstep = IAND(rndstep + nsamp, MASK)
    
    end subroutine rfr
    
    subroutine warmDRN(iseed)
        integer iseed, i
        
        call setDRN(iseed)
        rndstep = MAXDEPTH - 1
        do i = 1, 2
            call rfr(MAXDEPTH)
        end do
    
    end subroutine warmDRN
    
end module random
