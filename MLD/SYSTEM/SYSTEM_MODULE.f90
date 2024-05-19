module SYSTEM_MODULE
  implicit none 
  integer,parameter :: SP = kind(1.0)
  integer,parameter :: DP = selected_real_kind(2*precision(1.0_SP))
  integer,parameter :: QP = selected_real_kind(2*precision(1.0_DP))
  !integer,parameter :: MYPRE = selected_real_kind(15, 307)
  integer,parameter :: MYPRE = DP
  !
  integer,parameter :: IP1 = selected_int_kind(1)
  integer,parameter :: IP2 = selected_int_kind(2)
  integer,parameter :: IP4 = selected_int_kind(4)
  integer,parameter :: IP8 = selected_int_kind(8)
  !
  integer,parameter :: MYIP = IP8
end module SYSTEM_MODULE

