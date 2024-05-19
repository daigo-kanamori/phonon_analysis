module IBRAV_MODULE
 use SYSTEM_MODULE
 contains
  subroutine ibrav_LatticeVector(Lattice_Vector,celldm,ibrav)
    real(MYPRE) :: Lattice_Vector(3,3)
    real(MYPRE),intent(in) ::celldm(6)
    integer(MYIP),intent(in) ::           ibrav 
    !
    real(MYPRE) ::a(3),b(3),c(3)
    print *, "IBRAV TRANSFORM START"
    if(ibrav == 2) then 
       a = [0.50_MYPRE,0.0_MYPRE,0.50_MYPRE]
       b = [0.0_MYPRE,0.50_MYPRE,0.50_MYPRE]
       c = [0.50_MYPRE,0.50_MYPRE,0.0_MYPRE]
       !
       a(:) = a(:) * celldm(1)  
       b(:) = b(:) * celldm(1)  
       c(:) = c(:) * celldm(1)  
       !
      !a(:) = a(:) / 1.0d-10
      !b(:) = b(:) / 1.0d-10
      !c(:) = c(:) / 1.0d-10
      Lattice_Vector(:,1) = a(:)
      Lattice_Vector(:,2) = b(:)
      Lattice_Vector(:,3) = c(:)
    endif
    if(ibrav == 4) then 
      !
      a = [1.00_MYPRE,0.0_MYPRE,0.0_MYPRE]
      b = [-0.500_MYPRE,sqrt(3.00_MYPRE)/2.00_MYPRE,0.0_MYPRE]
      c = [0.0_MYPRE,0.0_MYPRE,celldm(3)]
      !
      a(:) = a(:) * celldm(1)  
      b(:) = b(:) * celldm(1)  
      c(:) = c(:) * celldm(1)  

      !a(:) = a(:) / 1.0d-10
      !b(:) = b(:) / 1.0d-10
      !c(:) = c(:) / 1.0d-10
      Lattice_Vector(:,1) = a(:)
      Lattice_Vector(:,2) = b(:)
      Lattice_Vector(:,3) = c(:)
    endif 
    print *, "IBRAV TRANSFORM DONE"
    !
    !
    return 
 end subroutine 
end module IBRAV_MODULE
