module my_diff_xs
use iso_c_binding
implicit none

! Fortran version of GENIE's KinePhaseSpace_t enumerated type
enum, bind (c) !:: KinePhaseSpace_t
  ! Fortran is case-insensitive, so "QQ2" should be interpreted as "q2"
  enumerator :: kPSNull = 0
  enumerator :: kPSfE, kPSxfE, kPSlogxfE, kPSxfEy, kPSlogxfEy, kPSyfE
  enumerator :: kPSlogyfE, kPSyfEx, kPSlogyfEx, kPSxyfE, kPSlogxlogyfE
  enumerator :: kPSQ2fE, kPSQD2fE, kPSlogQ2fE, kPSQ2fEW, kPSlogQ2fEW
  enumerator :: kPSQQ2fE, kPSQQ2fEW, kPSWfE, kPSWfEQ2, kPSWfEQQ2, kPSWQ2fE
  enumerator :: kPSWQD2fE, kPSW2Q2fE, kPSWlogQ2fE, kPSW2logQ2fE, kPSWQQ2fE
  enumerator :: kPSW2QQ2fE, kPSxytfE, kPSQ2yfE, kPSlogQ2logyfE, kPSTlctl
  enumerator :: kPSElOlOpifE, kPSElOlTpifE, kPSTkTlctl, kPSTnctnBnctl
end enum

contains

subroutine compute_my_diff_xsec(interaction, kps, xsec) &
  bind (c, name="compute_my_diff_xsec")

  type (c_ptr), intent (in), value :: interaction
  integer ( kind(kPSNull) ), intent (in), value :: kps
  real (c_double), intent (out) :: xsec

  write (*,*) "Fortran side: kps == ", kps
  if (kps .eq. kPSWQ2fE) then
    write (*,*) "kps == kPSWQ2fE"
  else
    write (*,*) "kps != kPSWQ2fE"
  end if

  xsec = 5.1
end subroutine

end module
