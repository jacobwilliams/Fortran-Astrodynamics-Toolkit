!*****************************************************************************************
!>
!  Test program for the string routines.

    program string_test

    use string_module

    implicit none

    character(len=:),allocatable :: s1

    s1 = 'Hello, World!'
    call uppercase(s1)
    if (s1 /= 'HELLO, WORLD!') error stop 'uppercase test failed: '//s1

    s1 = 'Hello, World!'
    call lowercase(s1)
    if (s1 /= 'hello, world!') error stop 'lowercase test failed: '//s1

    s1 = replace_char('banana', 'a', 'o')
    if (s1 /= 'bonono') error stop 'replace char test failed: '//s1

    s1 = reverse('stressed')
    if (s1 /= 'desserts') error stop 'reverse test failed: '//s1

    s1 = lchop('abc efg', 'abc')
    if (s1 /= ' efg' .and. len(s1)/=4) error stop 'lchop test failed: '//s1

    s1 = rchop('abc efg', 'efg')
    if (s1 /= 'abc' .and. s1(4:4) /= ' ') error stop 'rchop test failed: '//s1

    s1 = strip('   hello world   ')
    if (s1 /= 'hello world' .and. len_trim(s1)/=11) error stop 'strip test failed: '//s1

    s1 = strip('xxxxhelloxworldxxxx', 'x')
    if (s1 /= 'hello world' .and. len_trim(s1)/=11) error stop 'strip test failed: '//s1

    s1 = strip('abchelloxworldcccccab', 'bac')
    if (s1 /= 'hello world' .and. len_trim(s1)/=11) error stop 'strip test failed: '//s1

    s1 = strip('abchelloxworldcccccab ', 'bac')
    if (s1 /= 'hello worldcccccab ' .and. len(s1)/=19) error stop 'strip test failed: '//s1

    s1 = strip('abcaaabbbcccaaa', 'bac')
    if (s1 /= '' .and. len(s1)/=0) error stop 'strip test failed: '//s1

    end program string_test
!*****************************************************************************************