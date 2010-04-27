module error_handler
  implicit none
  
  private
  public :: fileError, subError
  
  contains

  !##################################################
  ! error handler for file opening errors -> kill program
  subroutine fileError(fname,error,callsub,flag)

    ! fname: file which threw error trying to open
    ! error: integer error returned from open() call
    ! callsub: string identifying which subroutine/ what location error happened
    ! flag:  =0 -> write problem;  <>0 -> read problem
    character(len=128), intent(in) :: fname, callsub
    integer, intent(in) :: error, flag
    
    if (flag == 0) then
       print *, trim(callsub),': error',error,'opening file ',trim(fname),' for writing'
    else
       print *, trim(callsub),': error',error,'opening file ',trim(fname),' for reading'
    end if
    stop 'quitting due to file error'

  end subroutine fileError

  !##################################################
  ! error handler for external subroutine errors (bessel functions, matrix decomposition)
  subroutine subError(sub,error,callsub)

    ! sub: name of subroutine which passed back error code
    ! error: integer error returned by subroutine
    ! callsub: string identifying which subroutine/ what location error happened
    character(len=128), intent(in) :: sub, callsub
    integer, intent(in) :: error
    
    print *, trim(callsub),': error',error,' returned from ',trim(sub)
    stop 'quitting due to subroutine error'
    
  end subroutine subError
end module error_handler
