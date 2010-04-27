program quadpack_prb
!
!******************************************************************************
!
!! QUADPACK_PRB runs the QUADPACK tests.
!
  implicit none
!
  call timestamp ( )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'QUADPACK_PRB'
  write ( *, '(a)' ) '  Tests for QUADPACK.'

  call test01
  call test02
  call test03
  call test04
  call test05
  call test06
  call test07
  call test08
  call test09
  call test10

  call test11
  call test12
  call test13
  call test14
  call test15

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'QUADPACK_PRB'
  write ( *, '(a)' ) '  Normal end of execution.'

  stop
end
subroutine test01
!
!******************************************************************************
!
!! TEST01 tests QAG.
!
!  QAG is an adaptive automatic integrator using a gauss-kronrod rule.
!
!  integrate cos(100*sin(x)) from 0 to pi.
!  exact answer is pi*j0(100), or roughly 0.06278740.
!
!  KEY chooses the order of the integration rule, from 1 to 6.
!
  implicit none
!
  real, parameter :: a = 0.0E+00
  real abserr
  real b
  real, parameter :: epsabs = 0.0E+00
  real, parameter :: epsrel = 0.001E+00
  real, external :: f02
  integer ier
  integer, parameter :: key = 6
  integer neval
  real pi
  real result
  real, parameter :: true = 0.06278740E+00
!
  b = pi()

  call qag ( f02, a, b, epsabs, epsrel, key, result, abserr, neval, ier )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST01'
  write ( *, '(a)' ) '  Test QAG'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Integrand is COS(100*SIN(X))'
  write ( *, '(a,g14.6)' ) '  Integral left endpoint A =    ', a
  write ( *, '(a,g14.6)' ) '  Integral right endpoint B =   ', b
  write ( *, '(a,g14.6)' ) '  Exact integral is             ', true
  write ( *, '(a,g14.6)' ) '  Estimated integral is         ', result
  write ( *, '(a,g14.6)' ) '  Estimated integral error =    ', abserr
  write ( *, '(a,g14.6)' ) '  Exact integral error =        ', true - result
  write ( *, '(a,i4)' ) '  Number of function evaluations, NEVAL = ', neval
  write ( *, '(a,i4)' ) '  Error return code IER = ', ier

  return
end
subroutine test02
!
!******************************************************************************
!
!! TEST02 tests QAGI.
!
!  QAGI is an adaptive quadrature routine for infinite intervals.
!
!  integrate log(x)/(1+100*x*x) from 0 to infinity.
!  exact answer is -pi*log(10)/20 = -0.3616892
!
!  give the type of infinity
!
!  inf=1 means a to infinity
!     -1      -infinity to a
!      2      -infinity to infinity
!
  implicit none
!
  real, parameter :: a = 0.0E+00
  real abserr
  real, parameter :: epsabs = 0.0E+00
  real, parameter :: epsrel = 0.001E+00
  real, external :: f05
  integer ier
  integer, parameter :: inf = 1
  integer neval
  real pi
  real result
  real true
!
  call qagi ( f05, a, inf, epsabs, epsrel, result, abserr, neval, ier )

  true=- pi() * log ( 10.0E+00 ) / 20.0E+00

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST02'
  write ( *, '(a)' ) '  Test QAGI'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Integrand is log(x)/(1+100*x*x)'
  write ( *, '(a,g14.6)' ) '  Integral left endpoint A =    ', a
  write ( *, '(a,g14.6)' ) '  Integral right endpoint B =    Infinity'
  write ( *, '(a,g14.6)' ) '  Exact integral is             ', true
  write ( *, '(a,g14.6)' ) '  Estimated integral is         ', result
  write ( *, '(a,g14.6)' ) '  Estimated integral error =    ', abserr
  write ( *, '(a,g14.6)' ) '  Exact integral error =        ', true - result
  write ( *, '(a,i4)' ) '  Number of function evaluations, NEVAL = ', neval
  write ( *, '(a,i4)' ) '  Error return code IER = ', ier

  return
end
subroutine test03
!
!******************************************************************************
!
!! TEST03 tests QAGP.
!
!
!  QAGP is an adaptive integrator that can handle singularities
!  of the integrand at user specified points,
!
!  Integrate 
!
!    x**3 * log(abs( (x*x-1)*(x*x-2) )) 
!
!  from 0 to 3.
!
!  Exact answer is 61*log(2)+77*log(7)/4 - 27
!
  implicit none
!
  integer, parameter :: npts = 2
  integer, parameter :: npts2 = 2 * npts
!
  real, parameter :: a = 0.0E+00
  real abserr
  real, parameter :: b = 3.0E+00
  real, parameter :: epsabs = 0.0E+00
  real, parameter :: epsrel = 0.001E+00
  real, external :: f04
  integer ier
  integer neval
  real points(npts2)
  real result
  real true
!
!  Singularity points:
!
  points(1) = 1.0E+00
  points(2) = sqrt ( 2.0E+00 )

  call qagp ( f04, a, b, npts2, points, epsabs, epsrel, result, abserr, &
    neval, ier )

  true = 61.0E+00 * log ( 2.0E+00 ) &
    + 77.0E+00 * log ( 7.0E+00 ) / 4.0E+00 - 27.0E+00

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST03'
  write ( *, '(a)' ) '  Test QAGP'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Integrand is x**3 * log(abs((x*x-1)*(x*x-2)))'
  write ( *, '(a,g14.6)' ) '  Integral left endpoint A =    ', a
  write ( *, '(a,g14.6)' ) '  Integral right endpoint B =   ', b
  write ( *, '(a,g14.6)' ) '  Exact integral is             ', true
  write ( *, '(a,g14.6)' ) '  Estimated integral is         ', result
  write ( *, '(a,g14.6)' ) '  Estimated integral error =    ', abserr
  write ( *, '(a,g14.6)' ) '  Exact integral error =        ', true - result
  write ( *, '(a,i4)' ) '  Number of function evaluations, NEVAL = ', neval
  write ( *, '(a,i4)' ) '  Error return code IER = ', ier

  return
end
subroutine test04
!
!******************************************************************************
!
!! TEST04 tests QAGS.
!
!  QAGS is an adaptive integrator for endpoint singularities.
!
!  integrate log(x)/sqrt(x) from 0 to 1.
!  exact answer is -4.
!
  implicit none
!
  real, parameter :: a = 0.0E+00
  real abserr
  real, parameter :: b = 1.0E+00
  real, parameter :: epsabs = 0.0E+00
  real, parameter :: epsrel = 0.001E+00
  real, external :: f03
  integer ier
  integer neval
  real result
  real, parameter :: true = -4.0E+00
!
  call qags ( f03, a, b, epsabs, epsrel, result, abserr, neval, ier )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST04'
  write ( *, '(a)' ) '  Test QAGS'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Integrand is LOG(X)/SQRT(X)'
  write ( *, '(a,g14.6)' ) '  Integral left endpoint A =    ', a
  write ( *, '(a,g14.6)' ) '  Integral right endpoint B =   ', b
  write ( *, '(a,g14.6)' ) '  Exact integral is             ', true
  write ( *, '(a,g14.6)' ) '  Estimated integral is         ', result
  write ( *, '(a,g14.6)' ) '  Estimated integral error =    ', abserr
  write ( *, '(a,g14.6)' ) '  Exact integral error =        ', true - result
  write ( *, '(a,i4)' ) '  Number of function evaluations, NEVAL = ', neval
  write ( *, '(a,i4)' ) '  Error return code IER = ', ier

  return
end
subroutine test05
!
!******************************************************************************
!
!! TEST05 tests QAWC.
!
!  QAWC is an adaptive integrator for finding the Cauchy
!  principal value of the integral of f(x)*w(x) over (a,b)
!  where w(x)=1/(x-c), c between a and b.
!
!  Integrate 1/(x*(5*x*x*x+6)) from -1 to 5
!
!  Exact answer is log(125/631) / 18 = -0.08994401
!
  implicit none
!
  real, parameter :: a = -1.0E+00
  real abserr
  real, parameter :: b = 5.0E+00
  real, parameter :: c = 0.0E+00
  real, parameter :: epsabs = 0.0E+00
  real, parameter :: epsrel = 0.001E+00
  real, external :: f09
  integer ier
  integer neval
  real result
  real true
!
  call qawc ( f09, a, b, c, epsabs, epsrel, result, abserr, neval, ier )

  true = log ( 125.0E+00 / 631.0E+00 ) / 18.0E+00

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST05'
  write ( *, '(a)' ) '  Test QAWC'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Integrand is 1/(x*(5*x**3+6)'
  write ( *, '(a,g14.6)' ) '  Integral left endpoint A =    ', a
  write ( *, '(a,g14.6)' ) '  Integral right endpoint B =   ', b
  write ( *, '(a,g14.6)' ) '  Point of singularity c =      ', c
  write ( *, '(a,g14.6)' ) '  Exact integral is             ', true
  write ( *, '(a,g14.6)' ) '  Estimated integral is         ', result
  write ( *, '(a,g14.6)' ) '  Estimated integral error =    ', abserr
  write ( *, '(a,g14.6)' ) '  Exact integral error =        ', true - result
  write ( *, '(a,i4)' ) '  Number of function evaluations, NEVAL = ', neval
  write ( *, '(a,i4)' ) '  Error return code IER = ', ier

  return
end
subroutine test06
!
!******************************************************************************
!
!! TEST06 tests QAWF.
!
!
!  QAWF handles fourier integration of f(x)*w(x) from
!  a to infinity, with w(x)=cos(omega*x) or sine(omega*x)
!
!  integrate cos(pi*x/2) /sqrt(x) from 0 to infinity.
!  exact answer is 1.0
!
  implicit none
!
  real, parameter :: a = 0.0E+00
  real abserr
  real, parameter :: epsabs = 0.001E+00
  real, external :: f07
  integer ier
  integer, parameter :: integr = 1
  integer neval
  real omega
  real pi
  real result
  real, parameter :: true = 1.0e+00
!
!  set argument of sine or cosine
!  set integr=1 for cosine, 2 for sine
!
  omega = 0.5E+00 * pi()

  call qawf ( f07, a, omega, integr, epsabs, result, abserr, neval, ier )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST06'
  write ( *, '(a)' ) '  Test QAWF'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Integrand is cos(pi*x/2)/sqrt(x)'
  write ( *, '(a,g14.6)' ) '  Integral left endpoint A =    ', a
  write ( *, '(a,g14.6)' ) '  Exact integral is             ', true
  write ( *, '(a,g14.6)' ) '  Estimated integral is         ', result
  write ( *, '(a,g14.6)' ) '  Estimated integral error =    ', abserr
  write ( *, '(a,g14.6)' ) '  Exact integral error =        ', true - result
  write ( *, '(a,i4)' ) '  Number of function evaluations, NEVAL = ', neval
  write ( *, '(a,i4)' ) '  Error return code IER = ', ier

  return
end
subroutine test07
!
!******************************************************************************
!
!! TEST07 tests QAWO.
!
!
!  QAWO integrates functions of the form f(x)*sin(omega*x)
!  or f(x)*cos(omega*x)
!
!  integrate log(x)*sin(10*pi*x)dx from 0 to 1
!  exact answer is -(gamma+log(10*pi)-ci(10*pi))/(10*pi)
!  =-0.1281316
!
!  gamma is euler's constant.
!  ci is the cosine integral ci(x)=integral(x to infinity) cos(v)/v dv.
!
!  specify integr=1 for cosine integral, 2 for sine integral
!
  implicit none
!
  real, parameter :: a = 0.0E+00
  real abserr
  real, parameter :: b = 1.0E+00
  real, parameter :: ci = -0.000842E+00
  real, parameter :: epsabs = 0.0E+00
  real, parameter :: epsrel = 0.001E+00
  real, external :: f06
  real, parameter :: gamma = 0.5772156649E+00
  integer ier
  integer, parameter :: integr = 2
  integer neval
  real omega
  real pi
  real result
  real true
!
!  set argument of the sine or cosine function
!
  omega = 10.0E+00 * pi()

  call qawo ( f06, a, b, omega, integr, epsabs, epsrel, result, abserr, &
    neval, ier )
!
!  I can't find an evaluation of ci(10*pi)
!  the following is an estimate
!
  true = - ( gamma + log ( 10.0E+00 * pi() ) - ci ) / ( 10.0E+00 * pi() )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST07'
  write ( *, '(a)' ) '  Test QAWO'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Integrand is log(x)*sin(10*pi*x)'
  write ( *, '(a,g14.6)' ) '  Integral left endpoint A =    ', a
  write ( *, '(a,g14.6)' ) '  Integral right endpoint B =   ', b
  write ( *, '(a,g14.6)' ) '  Exact integral is             ', true
  write ( *, '(a,g14.6)' ) '  Estimated integral is         ', result
  write ( *, '(a,g14.6)' ) '  Estimated integral error =    ', abserr
  write ( *, '(a,g14.6)' ) '  Exact integral error =        ', true - result
  write ( *, '(a,i4)' ) '  Number of function evaluations, NEVAL = ', neval
  write ( *, '(a,i4)' ) '  Error return code IER = ', ier

  return
end
subroutine test08
!
!******************************************************************************
!
!! TEST08 tests QAWS.
!
!
!  QAWS is an adaptive integrator for integrands with
!  algebraic or logarithmic singularities at the endpoints.
!
!  integrate log(x)/(1+(log(x))**2)**2 from 0 to 1
!  exact answer is -(ci(1.0)*sin(1.0)+(0.5*pi-si(1.0))*cos(1.0))/pi
!  (apparently this 'exact' answer is incorrect.)
!  (the numerical value seems correct, only by comparison with
!  the results of the program)
!  =-0.1892752
!
!  ci is the cosine integral ci(x)=integral(x to infinity) cos(v)/v dv
!  si is the sine integral si(x)=integral(0 to x) sin(v)/v dv
!  ci(1.0)=0.33740
!  si(1.0)=0.94608
!
  implicit none
!
  real, parameter :: a = 0.0E+00
  real abserr
  real, parameter :: alfa = 0.0E+00
  real, parameter :: b = 1.0E+00
  real, parameter :: beta = 0.0E+00
  real, parameter :: ci = 0.33740E+00
  real, parameter :: epsabs = 0.0E+00
  real, parameter :: epsrel = 0.001E+00
  real, external :: f08
  integer ier
  integer integr
  integer neval
  real pi
  real result
  real, parameter :: si = 0.94608E+00
  real true
!
!  describe integrand weight function
!  integr=2 means weight function is
!  (x-a)**alfa * (b-x)**beta * log(x-a)
!
  integr = 2

  call qaws ( f08, a, b, alfa, beta, integr, epsabs, epsrel, result, &
    abserr, neval, ier )
!
!  The following formula doesn't seem to come out correctly.
!
  true = - ( ci * sin ( 1.0E+00 ) + ( 0.5E+00 * pi() - si ) &
    * cos ( 1.0E+00 ) ) / pi()

  true = -0.1892752E+00

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST08'
  write ( *, '(a)' ) '  Test QAWS'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Integrand is log(x)/(1+(log(x))**2)**2'
  write ( *, '(a,g14.6)' ) '  Integral left endpoint A =    ', a
  write ( *, '(a,g14.6)' ) '  Integral right endpoint B =   ', b
  write ( *, '(a,g14.6)' ) '  Exact integral is             ', true
  write ( *, '(a,g14.6)' ) '  Estimated integral is         ', result
  write ( *, '(a,g14.6)' ) '  Estimated integral error =    ', abserr
  write ( *, '(a,g14.6)' ) '  Exact integral error =        ', true - result
  write ( *, '(a,i4)' ) '  Number of function evaluations, NEVAL = ', neval
  write ( *, '(a,i4)' ) '  Error return code IER = ', ier

  return
end
subroutine test09
!
!******************************************************************************
!
!! TEST09 tests QK15.
!
!
!  QK15 is a Gauss-Kronrod quadrature rule.
!
  implicit none
!
  real, parameter :: a = 0.0E+00
  real abserr
  real, parameter :: b = 1.0E+00
  real, external :: f01
  real resabs
  real resasc
  real result
  real, parameter :: true = - 4.0E+00 / 9.0E+00
!
  call qk15 ( f01, a, b, result, abserr, resabs, resasc )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST09'
  write ( *, '(a)' ) '  Test QK15'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Integrand is SQRT(X)*LOG(X)'
  write ( *, '(a,g14.6)' ) '  Integral left endpoint A =    ', a
  write ( *, '(a,g14.6)' ) '  Integral right endpoint B =   ', b
  write ( *, '(a,g14.6)' ) '  Exact integral is             ', true
  write ( *, '(a,g14.6)' ) '  Estimated integral is         ', result
  write ( *, '(a,g14.6)' ) '  Estimated integral error =    ', abserr
  write ( *, '(a,g14.6)' ) '  Exact integral error =        ', true - result
  write ( *, '(a,g14.6)' ) '  RESABS =                      ', resabs
  write ( *, '(a,g14.6)' ) '  RESASC =                      ', resasc

  return
end
subroutine test10
!
!******************************************************************************
!
!! TEST10 tests QK21.
!
!  QK21 is a Gauss-Kronrod quadrature rule.
!
  implicit none
!
  real, parameter :: a = 0.0E+00
  real abserr
  real, parameter :: b = 1.0E+00
  real, external :: f01
  real resabs
  real resasc
  real result
  real, parameter :: true = - 4.0E+00 / 9.0E+00
!
  call qk21 ( f01, a, b, result, abserr, resabs, resasc )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST10'
  write ( *, '(a)' ) '  Test QK21'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Integrand is SQRT(X)*LOG(X)'
  write ( *, '(a,g14.6)' ) '  Integral left endpoint A =    ', a
  write ( *, '(a,g14.6)' ) '  Integral right endpoint B =   ', b
  write ( *, '(a,g14.6)' ) '  Exact integral is             ', true
  write ( *, '(a,g14.6)' ) '  Estimated integral is         ', result
  write ( *, '(a,g14.6)' ) '  Estimated integral error =    ', abserr
  write ( *, '(a,g14.6)' ) '  Exact integral error =        ', true - result
  write ( *, '(a,g14.6)' ) '  RESABS =                      ', resabs
  write ( *, '(a,g14.6)' ) '  RESASC =                      ', resasc

  return
end
subroutine test11
!
!******************************************************************************
!
!! TEST11 tests QK31.
!
!  QK31 is a Gauss-Kronrod quadrature rule.
!
  implicit none
!
  real, parameter :: a = 0.0E+00
  real abserr
  real, parameter :: b = 1.0E+00
  real, external :: f01
  real resabs
  real resasc
  real result
  real, parameter :: true = - 4.0E+00 / 9.0E+00
!
  call qk31 ( f01, a, b, result, abserr, resabs, resasc )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST11'
  write ( *, '(a)' ) '  Test QK31'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Integrand is SQRT(X)*LOG(X)'
  write ( *, '(a,g14.6)' ) '  Integral left endpoint A =    ', a
  write ( *, '(a,g14.6)' ) '  Integral right endpoint B =   ', b
  write ( *, '(a,g14.6)' ) '  Exact integral is             ', true
  write ( *, '(a,g14.6)' ) '  Estimated integral is         ', result
  write ( *, '(a,g14.6)' ) '  Estimated integral error =    ', abserr
  write ( *, '(a,g14.6)' ) '  Exact integral error =        ', true - result
  write ( *, '(a,g14.6)' ) '  RESABS =                      ', resabs
  write ( *, '(a,g14.6)' ) '  RESASC =                      ', resasc

  return
end
subroutine test12
!
!******************************************************************************
!
!! TEST12 tests QK41.
!
!  QK41 is a Gauss-Kronrod quadrature rule.
!
  implicit none
!
  real, parameter :: a = 0.0E+00
  real abserr
  real, parameter :: b = 1.0E+00
  real, external :: f01
  real resabs
  real resasc
  real result
  real, parameter :: true = - 4.0E+00 / 9.0E+00
!
  call qk41 ( f01, a, b, result, abserr, resabs, resasc )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST12'
  write ( *, '(a)' ) '  Test QK41'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Integrand is SQRT(X)*LOG(X)'
  write ( *, '(a,g14.6)' ) '  Integral left endpoint A =    ', a
  write ( *, '(a,g14.6)' ) '  Integral right endpoint B =   ', b
  write ( *, '(a,g14.6)' ) '  Exact integral is             ', true
  write ( *, '(a,g14.6)' ) '  Estimated integral is         ', result
  write ( *, '(a,g14.6)' ) '  Estimated integral error =    ', abserr
  write ( *, '(a,g14.6)' ) '  Exact integral error =        ', true - result
  write ( *, '(a,g14.6)' ) '  RESABS =                      ', resabs
  write ( *, '(a,g14.6)' ) '  RESASC =                      ', resasc

  return
end
subroutine test13
!
!******************************************************************************
!
!! TEST13 tests QK51.
!
!  QK51 is a Gauss-Kronrod quadrature rule.
!
  implicit none
!
  real, parameter :: a = 0.0E+00
  real abserr
  real, parameter :: b = 1.0E+00
  real, external :: f01
  real resabs
  real resasc
  real result
  real, parameter :: true = - 4.0E+00 / 9.0E+00
!
  call qk51 ( f01, a, b, result, abserr, resabs, resasc )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST13'
  write ( *, '(a)' ) '  Test QK51'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Integrand is SQRT(X)*LOG(X)'
  write ( *, '(a,g14.6)' ) '  Integral left endpoint A =    ', a
  write ( *, '(a,g14.6)' ) '  Integral right endpoint B =   ', b
  write ( *, '(a,g14.6)' ) '  Exact integral is             ', true
  write ( *, '(a,g14.6)' ) '  Estimated integral is         ', result
  write ( *, '(a,g14.6)' ) '  Estimated integral error =    ', abserr
  write ( *, '(a,g14.6)' ) '  Exact integral error =        ', true - result
  write ( *, '(a,g14.6)' ) '  RESABS =                      ', resabs
  write ( *, '(a,g14.6)' ) '  RESASC =                      ', resasc

  return
end
subroutine test14
!
!******************************************************************************
!
!! TEST14 tests QK61.
!
!  QK61 is a Gauss-Kronrod quadrature rule.
!
  implicit none
!
  real, parameter :: a = 0.0E+00
  real abserr
  real, parameter :: b = 1.0E+00
  real, external :: f01
  real resabs
  real resasc
  real result
  real, parameter :: true = - 4.0E+00 / 9.0E+00
!
  call qk61 ( f01, a, b, result, abserr, resabs, resasc )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST14'
  write ( *, '(a)' ) '  Test QK61'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Integrand is SQRT(X)*LOG(X)'
  write ( *, '(a,g14.6)' ) '  Integral left endpoint A =    ', a
  write ( *, '(a,g14.6)' ) '  Integral right endpoint B =   ', b
  write ( *, '(a,g14.6)' ) '  Exact integral is             ', true
  write ( *, '(a,g14.6)' ) '  Estimated integral is         ', result
  write ( *, '(a,g14.6)' ) '  Estimated integral error =    ', abserr
  write ( *, '(a,g14.6)' ) '  Exact integral error =        ', true - result
  write ( *, '(a,g14.6)' ) '  RESABS =                      ', resabs
  write ( *, '(a,g14.6)' ) '  RESASC =                      ', resasc

  return
end
subroutine test15
!
!******************************************************************************
!
!! TEST15 tests QNG.
!
!  QNG is a nonadaptive automatic integrator using a Gauss-Kronrod or 
!  Patterson rule.
!
  implicit none
!
  real, parameter :: a = 0.0e+00
  real abserr
  real, parameter :: b = 1.0e+00
  real, parameter :: epsabs = 0.0E+00
  real, parameter :: epsrel = 0.001E+00
  real, external :: f01
  integer ier
  integer neval
  real result
  real, parameter :: true = - 4.0E+00 / 9.0E+00
!
  call qng ( f01, a, b, epsabs, epsrel, result, abserr, neval, ier )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST15'
  write ( *, '(a)' ) '  Test QNG'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '  Integrand is SQRT(X)*LOG(X)'
  write ( *, '(a,g14.6)' ) '  Integral left endpoint A =    ', a
  write ( *, '(a,g14.6)' ) '  Integral right endpoint B =   ', b
  write ( *, '(a,g14.6)' ) '  Exact integral is             ', true
  write ( *, '(a,g14.6)' ) '  Estimated integral is         ', result
  write ( *, '(a,g14.6)' ) '  Estimated integral error =    ', abserr
  write ( *, '(a,g14.6)' ) '  Exact integral error =        ', true - result
  write ( *, '(a,i4)' ) '  Number of function evaluations, NEVAL = ', neval
  write ( *, '(a,i4)' ) '  Error return code IER = ', ier

  return
end
function f01 ( x )
!
!******************************************************************************
!
!! F01 is the integrand function SQRT(X) * LOG(X).
!
  implicit none
!
  real f01
  real x
!
  if ( x <= 0.0E+00 )then
    f01 = 0.0E+00
  else
    f01 = sqrt ( x ) * log ( x )
  end if

  return
end
function f02 ( x )
!
!******************************************************************************
!
!! F02 is the integrand function COS(100*SIN(X)).
!
  implicit none
!
  real f02
  real x

  f02 = cos ( 100.0E+00 * sin ( x ) )

  return
end
function f03 ( x )
!
!******************************************************************************
!
!! F03 is the integrand function LOG(X)/SQRT(X).
!
  implicit none
!
  real f03
  real x

  if ( x <= 0.0E+00 ) then
    f03 = 0.0E+00
  else
    f03 = log ( x ) / sqrt ( x )
  end if

  return
end
function f04 ( x )
!
!******************************************************************************
!
!! F04 is the integrand function X^3 LOG((X^2-1)*(X^2-2))
!
  implicit none
!
  real f04
  real x

  f04 = x**3 * log ( abs ( ( x**2 - 1.0E+00 ) * ( x**2 - 2.0E+00 ) ) )

  return
end
function f05 ( x )
!
!******************************************************************************
!
!! F05 is the integrand function LOG(X)/(1+100X^2).
!
  implicit none
!
  real f05
  real x

  f05 = log ( x ) / ( 1.0E+00 + 100.0E+00 * x**2 )

  return
end
function f06 ( x )
!
!******************************************************************************
!
!! F06 is the integrand function LOG(X).
!
  implicit none
!
  real f06
  real x

  if ( x <= 0.0E+00 ) then
    f06 = 0.0E+00
  else
    f06 = log ( x )
  end if

  return
end
function f07 ( x )
!
!******************************************************************************
!
!! F07 is the integrand function 1/SQRT(X).
!
  implicit none
!
  real f07
  real x
!
  if ( x <= 0.0E+00 ) then
    f07 = 0.0E+00
  else
    f07 = 1.0E+00 / sqrt ( x )
  end if

  return
end
function f08 ( x )
!
!******************************************************************************
!
!! F08 is the integrand function 1/(1+LOG(X)**2)**2
!
  implicit none
!
  real f08
  real x
!
  if ( x > 0.0E+00 ) then
    f08 = 1.0e+00 / ( 1.0E+00 + log ( x )**2 )**2
  else
    f08 = 0.0e+00
  end if

  return
end
function f09 ( x )
!
!******************************************************************************
!
!! F09 is the integrand function 1 / ( 5 X^3 + 6 ).
!
  implicit none
!
  real f09
  real x
!
  f09 = 1.0E+00 / ( 5.0E+00 * x**3 + 6.0E+00 )

  return
end
