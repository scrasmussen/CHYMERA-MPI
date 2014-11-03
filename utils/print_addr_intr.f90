module print_addr_intr

  interface
     subroutine print_addr(addr) bind(C, name="print_addr")
       use :: iso_c_binding, only : C_PTR
       type(C_PTR), intent(in), value :: addr
     end subroutine print_addr
  end interface

end module print_addr_intr


subroutine test_print_addr
  use :: iso_c_binding, only : C_LOC
  use :: print_addr_intr, only : print_addr
  real :: A(10)

  call print_addr(C_LOC(A))

end subroutine test_print_addr

