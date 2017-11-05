program main
  use random
  implicit none
  integer, parameter :: L=16, N=L*L
  double precision, parameter :: XX=0.8d0
  integer :: target_spin(0:N-1), spin0(0:N-1), spin1(0:N-1), o(0:N)
  double precision :: hist(0:N), f(0:N), k
  integer :: seed, mc, E, ki

  read(*, *) seed
  call warmDRN(seed)
  call init_spin

  f(:) = 0.d0
  spin0(:) = 0
  o(:) = 0
  E = sum(target_spin(:))
  ki = 1
  k = 1.d0
  write(*, "('init done')")
  do
     hist(:) = 0.d0
     do
        do mc = 1, N*N
           call one_mc
        end do
        if (is_flat() == 1) then
           exit
        end if
     end do
     if (k < 1.d-8) then
        exit
     end if
     k = k * 0.5d0
     ki = ki + 1
  end do
  call output

contains

  subroutine output
    integer :: i
    character(99) :: name
    write(name, "('data/L', i4.4, '_seed', i4.4, '.spin')") L, seed
    open(10, file=trim(name), action="write")
    do i = 0, L - 1
       write(10, "(9999i1)") target_spin(i*L:(i+1)*L)
    end do
    close(10)

    write(name, "('data/L', i4.4, '_seed', i4.4, '.hist')") L, seed
    open(10, file=trim(name), action="write")
    do i = 0, N
       write(10, *) i, f(i)
    end do
    close(10)
  end subroutine output

  subroutine init_spin
    integer :: i

    call rfr(N)
    do i = 0, N-1
       if (ur(i) < 0.5d0) then
          target_spin(i) = 0
       else
          target_spin(i) = 1
       end if
    end do
  end subroutine init_spin

  subroutine one_mc
    integer :: mmc, i, E2
    double precision :: df

    call rfr(N*2)
    do mmc = 0, N - 1
       i = mod(ir(mmc*2), N)
       spin0(i) = 1 - spin0(i)
       call calc_E(E2)
       df = f(E) - f(E2)
       if (ur(mmc*2+1) < exp(df)) then
          E = E2
       else
          spin0(i) = 1 - spin0(i)
       end if
       hist(E) = hist(E) + 1.d0
       f(E) = f(E) + k
       o(E) = 1
    end do
  end subroutine one_mc

  subroutine calc_E(E2)
    integer, intent(out) :: E2
    integer :: c, x, y
    do y = 0, L - 1
       do x = 0, L - 1
          c = spin0(mod(y-1+L, L)*L + mod(x-1+L, L))
          c = c + spin0(mod(y-1+L, L)*L + x)
          c = c + spin0(mod(y-1+L, L)*L + mod(x+1, L))
          c = c + spin0(y*L + mod(x-1+L, L))
          c = c + spin0(y*L + x)
          c = c + spin0(y*L + mod(x+1, L))
          c = c + spin0(mod(y+1, L)*L + mod(x-1+L, L))
          c = c + spin0(mod(y+1, L)*L + x)
          c = c + spin0(mod(y+1, L)*L + mod(x+1, L))
          if (c == 3 .or. c == 4) then
             spin1(y*L + x) = 1
          else
             spin1(y*L + x) = 0
          end if
       end do
    end do
    E2 = sum(abs(spin1(:) - target_spin(:)))
  end subroutine calc_E

  integer function is_flat()
    integer :: i, hist_n
    double precision :: min_hist=-1.d0, min_f=-1.d0

    hist_n = 0
    do i = 0, N
       if (o(i) == 1) then
          hist_n = hist_n + 1
          if (min_hist < 0.d0 .or. min_hist > hist(i)) min_hist = hist(i)
          if (min_f < 0.d0 .or. min_f > f(i)) min_f = f(i)
       end if
    end do
    f(:) = f(:) - min_f
    if (min_hist > sum(hist) * XX / hist_n) then
       is_flat = 1
    else
       is_flat = 0
    end if
    write(*, "(i3, f8.4, i9)") ki, min_hist / sum(hist) * hist_n, hist_n
  end function is_flat

end program main
