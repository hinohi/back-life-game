program main
  use random
  implicit none
  integer, parameter :: L=8, N=L*L
  double precision, parameter :: XX=0.8d0
  integer :: target_spin(0:N-1), spin0(0:N-1), spin1(0:N-1), o(0:N)
  double precision :: hist(0:N), f(0:N), k
  integer :: seed, p, mc, E, ki

  read(*, *) seed, p
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
        do mc = 1, N*N*10
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
    double precision :: min_F
    character(99) :: name
    write(name, "('data/L', i4.4, '_seed', i4.4, '_p', i4.4, '.spin')") L, seed, p
    open(10, file=trim(name), action="write")
    do i = 0, L - 1
       write(10, "(9999i1)") target_spin(i*L:(i+1)*L)
    end do
    close(10)

    write(name, "('data/L', i4.4, '_seed', i4.4, '_p', i4.4, '.hist')") L, seed, p
    open(10, file=trim(name), action="write")
    min_f = -1.d0
    do i = 0, N
       if (o(i) == 1) then
          if (min_f < 0.d0 .or. min_f > f(i)) min_F = f(i)
       end if
    end do
    do i = 0, N
       if (o(i) == 1) then
          write(10, *) i, f(i) - min_f
       end if
    end do
    close(10)
  end subroutine output

  subroutine init_spin
    integer :: i, done

    target_spin(:) = 0
    done = 0
    do
       if (done >= p) exit
       call rfr(1)
       i = mod(ir(1), N)
       if (target_spin(i) == 0) then
          target_spin(i) = 1
          done = done + 1
       end if
    end do
  end subroutine init_spin

  subroutine one_mc
    integer :: mmc, i, E2
    double precision :: df

    call rfr(N*2)
    do mmc = 0, N - 1
       i = mod(ir(mmc*2+1), N)
       spin0(i) = 1 - spin0(i)
       call calc_E(E2)
       df = f(E) - f(E2)
       if (df > 0.d0 .or. ur(mmc*2+2) < exp(df)) then
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
    double precision :: min_hist, min_f

    min_hist = -1.d0
    min_f = -1.d0
    hist_n = 0
    do i = 0, N
       if (o(i) == 1) then
          hist_n = hist_n + 1
          if (min_hist < 0.d0 .or. min_hist > hist(i)) min_hist = hist(i)
          if (min_f < 0.d0 .or. min_f > f(i)) min_f = f(i)
       end if
    end do
    if (min_hist > sum(hist) * XX / hist_n) then
       is_flat = 1
    else
       is_flat = 0
    end if
    write(*, "(i3, f8.4, f8.4)") ki, min_hist / sum(hist) * hist_n, hist_n*1.d0/(N+1)
  end function is_flat

end program main
