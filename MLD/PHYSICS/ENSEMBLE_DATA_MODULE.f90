module ENSEMBLE_DATA_MODULE
  use SYSTEM_MODULE
  use FFTW_MODULE
  use VECTOR_MODULE
  use FITTING_MODULE
  use ENSEMBLE_MODULE
  use KAHAN_MODULE
  implicit none
  !
  !=========================================================================================!
  !                   IMPLEMENT DATA STRUCTURE                                              !
  !=========================================================================================!
  type Ensemble_data_Varialbles
    type(Vector_d),allocatable :: Displacement(:)
    type(Vector_d),allocatable :: Velocity(:)
    type(Vector_d),allocatable :: Force(:)
  end type 
  !
  type Ensemble_data_Information
    type(Vector_d) :: Mass
  end type Ensemble_data_Information
  !
  type Ensemble_data
    type(Ensemble_data_Varialbles) :: state
    type(Ensemble_data_Information) :: info
  end type
  !
  type GrandEnsemble_data
    type(Ensemble_data),allocatable :: System(:)
    integer(MYIP),allocatable :: Ensemble_data_list(:,:,:)
    integer(MYIP) :: Ensemble_data_Number
    integer(MYIP)  :: time
    real(MYPRE)   :: step
  end type
  !
  type Cell_energy_space
    type(KahanSum),allocatable :: energy_time(:)
    integer(MYIP) :: start,fin
  end type Cell_energy_space
  !
  type SED_space_0
    complex(MYPRE),allocatable :: xin(:)
    complex(MYPRE),allocatable :: xout(:)
    type(KahanSum),allocatable :: sed(:)
    integer(MYIP) :: start,fin
  end type SED_space_0
  !
  type SED_space_1
    complex(MYPRE),allocatable :: xin(:,:)
    complex(MYPRE),allocatable :: xout(:,:)
    type(KahanSum),allocatable :: sed(:,:)
    integer(MYIP) :: start,fin
    integer(MYIP) :: xlength
  end type SED_space_1
  !
  !
  contains 
  !

  !=============================================================================================================!
  !   GrandEnsemble -> GrandEnsemble_data
  !=============================================================================================================!
  !
  subroutine Copy_GrandEnsemble_data_info(sys_data,sys,fin)
    type(GrandEnsemble_data) :: sys_data
    type(GrandEnsemble),intent(in) :: sys
    integer(MYIP) ,intent(in) :: fin

    integer(MYIP) :: il,jl,kl
    integer(MYIP) :: ir,jr,kr
    integer(MYIP) :: i
    integer(MYIP) :: ide
    integer(MYIP) :: Nx(2),Ny(2),Nz(2)
    integer(MYIP) :: nfree
    real(MYPRE) :: time_step

    sys_data%Ensemble_data_Number = sys%Ensemble_Number


    Nx = [lbound(sys%Ensemble_list,1),ubound(sys%Ensemble_list,1)]
    Ny = [lbound(sys%Ensemble_list,2),ubound(sys%Ensemble_list,2)]
    Nz = [lbound(sys%Ensemble_list,3),ubound(sys%Ensemble_list,3)]

    sys_data%step = sys%step


    call allocate_GrandEnsemble_data(sys_data,Nx,Ny,Nz,fin)

    do kl = Nz(1),Nz(2)
    do jl = Ny(1),Ny(2)
    do il = Nx(1),Nx(2)
      ide = sys%Ensemble_list(il,jl,kl)
      sys_data%Ensemble_data_list(il,jl,kl) = ide
      if( allocated(sys%System(ide)%info%Mass%vector)) then 
        nfree = size(sys%System(ide)%info%Mass%vector)
        call allocate_Ensemble_data(sys_data%System(ide),nfree)
        sys_data%System(ide)%info%Mass%vector(:) = sys%System(ide)%info%Mass%vector(:)
      end if
    end do 
    end do 
    end do 
    !
    return 
  end subroutine Copy_GrandEnsemble_data_info
  !
  subroutine Copy_GrandEnsemble_data_state(sys_data,sys,now)
    type(GrandEnsemble_data) :: sys_data
    type(GrandEnsemble),intent(in) :: sys
    integer(MYIP),intent(in) :: now

    integer(MYIP) :: il,jl,kl
    integer(MYIP) :: ide
    integer(MYIP) :: i
    integer(MYIP) :: nfree
    integer(MYIP) :: Nx(2),Ny(2),Nz(2)

    Nx = [lbound(sys_data%Ensemble_data_list,1),ubound(sys_data%Ensemble_data_list,1)]
    Ny = [lbound(sys_data%Ensemble_data_list,2),ubound(sys_data%Ensemble_data_list,2)]
    Nz = [lbound(sys_data%Ensemble_data_list,3),ubound(sys_data%Ensemble_data_list,3)]
    !$omp parallel 
    !$omp do private(jl,il,ide,nfree,i)
    do kl = Nz(1),Nz(2)
    do jl = Ny(1),Ny(2)
    do il = Nx(1),Nx(2)
      ide = sys_data%Ensemble_data_list(il,jl,kl)
      nfree = size(sys_data%System(ide)%info%Mass%vector)
      do i = 1,nfree
        sys_data%System(ide)%state%Displacement(now)%vector(i) = sys%System(ide)%state%Displacement(1)%vector(i)%sum
        sys_data%System(ide)%state%Velocity(now)%vector(i) = sys%System(ide)%state%Velocity(1)%vector(i)%sum
        sys_data%System(ide)%state%Force(now)%vector(i) = sys%System(ide)%state%Force(1)%vector(i)%sum
      end do 
    end do
    end do
    end do
    !$omp end do
    !$omp end parallel

    return 
  end subroutine Copy_GrandEnsemble_data_state
  !
  !=============================================================================================================!
  !                               WRITE DATA (STATE) IN TIME YOU SPECIFY                                        !
  !=============================================================================================================!
  !
  subroutine Write_Ensemble_data_state(ensem,fd,it)
    type(Ensemble_data),intent(in) :: ensem
    integer(MYIP),intent(in) :: fd,it
    integer(MYIP) :: freenumber
    integer(MYIP)  :: i
    !
    freenumber = size(ensem%info%Mass%vector)
    !
    write(fd,*) (ensem%state%Displacement(it)%vector(i),i=1,freenumber)
    write(fd,*) (ensem%state%Velocity(it)%vector(i),i=1,freenumber)
    write(fd,*) (ensem%state%Force(it)%vector(i),i=1,freenumber)
    !
    return 
    !
  end subroutine 
  !
  subroutine Write_GrandEnsemble_data_state(sys,fd,it)
    type(GrandEnsemble_data),intent(in) :: sys
    integer(MYIP),intent(in) :: fd,it
    integer(MYIP) :: Nx(2),Ny(2),Nz(2)
    integer(MYIP),allocatable :: xpos(:),ypos(:),zpos(:)
    integer(MYIP) :: il,jl,kl
    integer(MYIP) :: ide,id
    !
    allocate(xpos(sys%Ensemble_data_Number))  
    allocate(ypos(sys%Ensemble_data_Number))  
    allocate(zpos(sys%Ensemble_data_Number))  
    !
    xpos(:) = 0
    ypos(:) = 0
    zpos(:) = 0
    !
    Nx = [lbound(sys%Ensemble_data_list,1),ubound(sys%Ensemble_data_list,1)]
    Ny = [lbound(sys%Ensemble_data_list,2),ubound(sys%Ensemble_data_list,2)]
    Nz = [lbound(sys%Ensemble_data_list,3),ubound(sys%Ensemble_data_list,3)]
    !
    do kl = Nz(1),Nz(2)
    do jl = Ny(1),Ny(2)
    do il = Nx(1),Nx(2)
      ide = sys%Ensemble_data_list(il,jl,kl)
      xpos(ide) = il
      xpos(ide) = jl
      zpos(ide) = kl
    end do 
    end do 
    end do 
    !
    write(fd,*) sys%Ensemble_data_Number
    !
    do id = 1,sys%Ensemble_data_Number
      write(fd,*) xpos(id),ypos(id),zpos(id)
      call Write_Ensemble_data_state(sys%System(id),fd,it)
    end do 
  end subroutine  
  !
  !=============================================================================================================!
  !                                 READ PRE DATA AND ALLOCATE DATA STURUCTURE                                  !
  !=============================================================================================================!
  !
  !
  subroutine Read_GrandEnsemble_data_info(sys_data,fd)
    type(GrandEnsemble_data) :: sys_data
    integer(MYIP),intent(in) :: fd
    integer(MYIP) :: Nx(2) ,Ny(2), Nz(2)
    integer(MYIP) :: ir,jr,kr
    integer(MYIP) :: time_step
    integer(MYIP) :: kind_cell
    integer(MYIP) :: ide
    integer(MYIP) :: i
    !
    print *, "START READ PRE DATA "
    read(fd,*) time_step
    read(fd,*) Nx(1),Nx(2)
    read(fd,*) Ny(1),Ny(2)
    read(fd,*) Nz(1),Nz(2)
    read(fd,*) kind_cell
    !
    sys_data%Ensemble_data_Number = kind_cell
    !
    call allocate_GrandEnsemble_data(sys_data,Nx,Ny,Nz,time_step)
    !
    !
    do i = 1,kind_cell
      read(fd,*)ir,jr,kr,ide
      sys_data%Ensemble_data_list(ir,jr,kr) = ide
      call Read_Ensemble_data_info(sys_data%System(ide),fd)
    end do 
    !
    !
    print *, "ALLOCATE GRAND ENSEMBLE DATA DONE"
    !
    return 
  end subroutine Read_GrandEnsemble_data_info
  !
  subroutine Read_Ensemble_data_info(ensem,fd)
    type(Ensemble_data) :: ensem
    integer(MYIP),intent(in) ::  fd
    integer(MYIP) :: freenumber
    integer(MYIP) :: i
    !
    read(fd,*)freenumber 
    call allocate_Ensemble_data(ensem,freenumber)
    !
    read(fd,*)(ensem%info%Mass%vector(i),i=1,freenumber)
    !
    !
    return 
  end subroutine Read_Ensemble_data_info
  !
  subroutine allocate_Ensemble_data(ensem,freenumber)
    type(Ensemble_data) :: ensem
    integer(MYIP),intent(in) :: freenumber
    integer(MYIP) :: it
    !
    !
    do it = lbound(ensem%state%Displacement,1),ubound(ensem%state%Displacement,1)
      call allocate_Vector_d(ensem%state%Displacement(it),1,freenumber)
    end do 
    do it = lbound(ensem%state%Velocity,1),ubound(ensem%state%Velocity,1)
      call allocate_Vector_d(ensem%state%Velocity(it),1,freenumber)
    end do
    do it = lbound(ensem%state%Force,1),ubound(ensem%state%Force,1)
      call allocate_Vector_d(ensem%state%Force(it),1,freenumber)
    end do 
    !
    call allocate_Vector_d(ensem%info%Mass,1,freenumber)
    !
    return 
  end subroutine allocate_Ensemble_data
  !
  !
  subroutine allocate_GrandEnsemble_data(sys_data,Nx,Ny,Nz,time)
    type(GrandEnsemble_data) :: sys_data
    integer(MYIP),intent(in) :: Nx(2),Ny(2),Nz(2) 
    integer(MYIP),intent(in) :: time
    integer(MYIP) :: cell_number
    integer(MYIP) :: ide
    !
    allocate(sys_data%Ensemble_data_list(Nx(1):Nx(2),Ny(1):Ny(2),Nz(1):Nz(2)))
    !
    cell_number = sys_data%Ensemble_data_Number
    allocate(sys_data%System(cell_number))
    !
    do ide = 1,cell_number
      allocate(sys_data%System(ide)%state%Displacement(0:time))
      allocate(sys_data%System(ide)%state%Velocity(0:time))
      allocate(sys_data%System(ide)%state%Force(0:time))
    end do 
    !
    sys_data%time = time
    !
    return 
  end subroutine allocate_GrandEnsemble_data
  !
  !
  subroutine Read_GrandEnsemble_data(sys_data,info_fd,data_fd)
    type(GrandEnsemble_data) :: sys_data
    integer(MYIP),intent(in) :: data_fd,info_fd
    !
    print *, "START SET GRAND ENSEMBLE DATA"
    call Read_GrandEnsemble_data_info(sys_data,info_fd)
    call Read_GrandEnsemble_data_state(sys_data,data_fd)
    !
    print *, "DONE"
    return 
  end subroutine  Read_GrandEnsemble_data
  !
  !==================================================================================================!
  !                          READ   ENSEMBLE STATE                                     !  
  !==================================================================================================!
  !
  !
  subroutine Read_Ensemble_data_state(ensem,fd,it)
    type(Ensemble_data) :: ensem
    integer(MYIP),intent(in) :: fd,it
    integer(MYIP) :: freenumber
    integer(MYIP) :: i
    !
    freenumber = size(ensem%info%Mass%vector)
    read(fd,*)(ensem%state%Displacement(it)%vector(i),i=1,freenumber)    
    read(fd,*)(ensem%state%Velocity(it)%vector(i),i=1,freenumber)    
    read(fd,*)(ensem%state%Force(it)%vector(i),i=1,freenumber)    
    !
    return 
  end subroutine Read_Ensemble_data_state
  !
  subroutine Read_GrandEnsemble_data_state(sys,fd)
    type(GrandEnsemble_data) :: sys
    integer(MYIP) , intent(in) :: fd
    integer(MYIP) :: it
    integer(MYIP) :: ide,i
    integer(MYIP) :: time_now
    !
    print *, "START GRAND ENSEMBLE DATA STATE"
    !
    !
    do it = 0,sys%time
      read(fd,*)time_now
      do i = 1,sys%Ensemble_data_Number
        read(fd,*)ide
        call Read_Ensemble_data_state(sys%System(ide),fd,time_now)
      end do 
    end do 
    !
    print *, "READ GRAND ENSEMBLE DATA DONE"
  end subroutine Read_GrandEnsemble_data_state
  !
  !==================================================================================!
  !                            FREE SED SPACE                                        !
  !==================================================================================!
  !
  subroutine Free_SED_space_0(sed_space)
    type(SED_space_0) :: sed_space
    !
    if(allocated(sed_space%xin))then 
      deallocate(sed_space%xin)
      deallocate(sed_space%xout)
      deallocate(sed_space%sed)
    end if
    sed_space%start = 0
    sed_space%fin = 0
    !
    return 
  end subroutine Free_SED_space_0 
  !
  subroutine Free_SED_space_1(sed_space) 
    type(SED_space_1) :: sed_space
    !
    if(allocated(sed_space%xin))then 
      deallocate(sed_space%xin)
      deallocate(sed_space%xout)
      deallocate(sed_space%sed)
    end if
    sed_space%start = 0
    sed_space%fin = 0
    !
    return 
  end subroutine Free_SED_space_1
  !
  !
  !
  !===================================================================================!
  !                 ANALYZE KINETIC ENERGY                                            !                  
  !===================================================================================!
  !
  subroutine Calculate_Kinetic_Energy_GrandEnsemble_data(sys_data,kegy,now)
    type(GrandEnsemble_data),intent(in) :: sys_data
    real(MYPRE) :: kegy
    integer(MYIP),intent(in) :: now 

    integer(MYIP) :: ide
    type(KahanSum) :: kaha

    call kahan_zero(kaha)

    kegy = 0_MYPRE

    do ide = 1,sys_data%Ensemble_data_Number
      call Calculate_Kinetic_Energy_Ensemble_data(sys_data%System(ide),kaha,now)  
    end do

    kegy = kaha%sum
    
    return 
  end subroutine Calculate_Kinetic_Energy_GrandEnsemble_data  !

  subroutine Calculate_Kinetic_Energy_Ensemble_data(ensem_data,kegy,now)
    type(Ensemble_data),intent(in) :: ensem_data
    type(KahanSum) :: kegy
    integer(MYIP) :: now

    real(MYPRE) :: mass,velo
    real(MYPRE) :: tmp
    integer(MYIP) :: nfree
    integer(MYIP) :: freenumber

    freenumber = size(ensem_data%info%Mass%vector)

    do nfree = 1,freenumber
      mass = ensem_data%info%Mass%vector(nfree)
      velo = ensem_data%state%Velocity(now)%vector(nfree)
      tmp = 0.5_MYPRE * mass * velo * velo
      call add_kahan(kegy,tmp)
    end do 
    return 
  end subroutine Calculate_Kinetic_Energy_Ensemble_data

  subroutine Calculate_Potensial_Energy_GrandEnsemble_data(sys_data,pegy,now)
    type(GrandEnsemble_data),intent(in) :: sys_data
    real(MYPRE) :: pegy
    integer(MYIP),intent(in) :: now 

    integer(MYIP) :: ide
    type(KahanSum) :: kaha

    call kahan_zero(kaha)

    pegy = 0_MYPRE

    do ide = 1,sys_data%Ensemble_data_Number
      call Calculate_Potensial_Energy_Ensemble_data(sys_data%System(ide),kaha,now)  
    end do

    pegy = kaha%sum
    
    return 
  end subroutine Calculate_Potensial_Energy_GrandEnsemble_data  !
  !
  subroutine Calculate_Potensial_Energy_Ensemble_data(ensem_data,pegy,now)
    type(Ensemble_data),intent(in) :: ensem_data
    type(KahanSum) :: pegy
    integer(MYIP) :: now

    real(MYPRE) :: frc,disp
    real(MYPRE) :: tmp
    integer(MYIP) :: nfree
    integer(MYIP) :: freenumber

    freenumber = size(ensem_data%info%Mass%vector)

    do nfree = 1,freenumber
      frc = ensem_data%state%Force(now)%vector(nfree)
      disp = ensem_data%state%Displacement(now)%vector(nfree)
      tmp = -0.5_MYPRE * frc * disp
      call add_kahan(pegy,tmp)
    end do 
    return 
  end subroutine Calculate_Potensial_Energy_Ensemble_data
  !===================================================================================!
  !                 ANALYZE KINETIC ENERGY                                            !                  
  !===================================================================================!
  !
  subroutine Calculate_kinetic_energy_freq(sys,xout,start,fin) 
    type(GrandEnsemble_data),intent(in) :: sys
    integer(MYIP),intent(in) :: start,fin
    complex(MYPRE),allocatable ::  xout(:)
    complex(MYPRE),allocatable ::  xin(:)
    real(MYPRE) :: tmp
    type(KahanSum) :: kaha
    integer(MYIP) :: it,i,ide
    integer(MYIP) :: freenumber
    !
    if(allocated(xout)) then 
      deallocate(xout)
    end if
    !
    allocate(xout(0:(fin-start)))
    allocate(xin(0:(fin-start)))
    !
    do it = start,fin
      call kahan_zero(kaha)
      do ide = 1,sys%Ensemble_data_Number
        freenumber = size(sys%System(ide)%info%Mass%vector)
        do i = 1,freenumber
          tmp = 0.5 * sys%System(ide)%info%Mass%vector(i) &
              * sys%System(ide)%state%Velocity(it)%vector(i) &
              * sys%System(ide)%state%Velocity(it)%vector(i) 
          call add_Kahan(kaha,tmp)
          end do 
        end do 
       xin(it-start) = CMPLX(kaha%sum,0d0,kind=kind(0d0))
    end do 
    !
    call fftw_1d(xin,xout,FFTW_FORWARD)
    !
    deallocate(xin)
    !
    return 
  end subroutine Calculate_kinetic_energy_freq
  !
  !===================================================================================!
  !                         ANALYZE TIME CONSTANTS                                    !
  !===================================================================================!
  !
  subroutine allocate_Cell_energy_time(ene_t,start,fin)
    type(Cell_energy_space) :: ene_t
    integer(MYIP),intent(in) :: start,fin
    integer(MYIP) :: it
    !
    allocate(ene_t%energy_time(0:(fin-start)))
    ene_t%fin = fin
    ene_t%start = start
    do it = 0,fin-start
      call kahan_zero(ene_t%energy_time(it))
    end do 
    !
    return 
  end subroutine allocate_Cell_energy_time
  !
  !
  subroutine Set_Cell_energy_time(sys,ene_t,ide,freenumber)
    type(GrandEnsemble_data),intent(in) :: sys
    type(Cell_energy_space) :: ene_t
    integer(MYIP),intent(in) :: ide,freenumber
    integer(MYIP) :: it
    real(MYPRE) :: ke,pe,mass,disp,velo,frc
    !
    !ide = sys%Ensemble_data_list(Lx,Ly,Lz)
    mass = sys%System(ide)%info%Mass%vector(freenumber)
    !
    do it = ene_t%start,ene_t%fin
      disp = sys%System(ide)%state%Displacement(it)%vector(freenumber)
      velo = sys%System(ide)%state%Velocity(it)%vector(freenumber)
      frc = sys%System(ide)%state%Force(it)%vector(freenumber)
      ke = 0.5d0 * mass * velo * velo
      pe = -0.5d0 * frc * disp
      !
      call add_kahan(ene_t%energy_time(it-ene_t%start),ke)
      call add_kahan(ene_t%energy_time(it-ene_t%start),pe)
      !
    end do 
    !
    return 
  end subroutine Set_Cell_energy_time 
  !
  !
  subroutine Calculate_Time_Constants(ene_t,U,tau)
    type(Cell_energy_space),intent(in) :: ene_t
    real(MYPRE) :: U,tau
    integer(MYIP) :: it
    real(MYPRE) :: a,b
    real(MYPRE),allocatable :: x(:),y(:)
    !
    allocate(x(0:(ene_t%fin-ene_t%start)))
    allocate(y(0:(ene_t%fin-ene_t%start)))
    !
    do it = 0,ene_t%fin-ene_t%start
      x(it) = it
      y(it) = ene_t%energy_time(it)%sum
    end do 
    !
    !ax + b;
    call LinearFitting(x,y,a,b)
    !
    print *, a,b
    U = b
    tau = abs(a/(U*0.5d0))
    tau = 1d0/tau
    !
    return 
  end subroutine Calculate_Time_Constants 
  !
  !
  !=========================================================================================!
  !                           SED IMPLEMENTAION                                             !
  !=========================================================================================!
  !
  subroutine allocate_SED_0(sed_space,start,fin)
    type(SED_space_0) :: sed_space
    integer(MYIP),intent(in) :: start,fin
    integer(MYIP) :: it
    integer(MYIP) :: N 
    !
    N = fin - start
    !
    sed_space%fin = fin
    sed_space%start = start
    allocate(sed_space%xin(0:N))
    allocate(sed_space%sed(0:N))
    !
    sed_space%xin(:) = CMPLX(0d0,0d0,kind=kind(0d0))
    !
    !
    do it = 0,N
      call kahan_zero(sed_space%sed(it))
    end do 
    !
    return 
  end subroutine allocate_SED_0
  !
  subroutine Set_xin_SED_0(ensem,sed_space,freenumber)
    type(Ensemble_data),intent(in) :: ensem
    type(SED_space_0) :: sed_space
    integer(MYIP),intent(in) ::  freenumber
    integer(MYIP) :: it
    !
    do it = sed_space%start,sed_space%fin
      sed_space%xin(it-sed_space%start) = ensem%state%Velocity(it)%vector(freenumber)
    end do 
    !
    return
  end subroutine Set_xin_SED_0
  !
  subroutine Calculate_SED_0(sed_space) 
    type(SED_space_0) :: sed_space
    real(MYPRE) :: tmp
    integer(MYIP) :: it
    !
    !
    call fftw_1d(sed_space%xin,sed_space%xout,FFTW_FORWARD)
    !
    do it = 0,(sed_space%fin-sed_space%start)
      tmp = abs(sed_space%xout(it))
      tmp = tmp * tmp
      call add_kahan(sed_space%sed(it),tmp)
    end do 
    !
    return 
  end subroutine Calculate_SED_0
  !
  !
  subroutine Calculate_SED_0_Kinetic_Energy(sed_space,mass) 
    type(SED_space_0) :: sed_space
    real(MYPRE),intent(in) :: mass
    real(MYPRE) :: tmp
    integer(MYIP) :: it
    !
    !
    call fftw_1d(sed_space%xin,sed_space%xout,FFTW_FORWARD)
    !
    do it = 0,(sed_space%fin-sed_space%start)
      tmp = abs(sed_space%xout(it))
      tmp = mass * tmp * tmp
      call add_kahan(sed_space%sed(it),tmp)
    end do 
    !
    return 
  end subroutine Calculate_SED_0_Kinetic_Energy
  !
  subroutine Calculate_SED_0_Ensemble_energy(ensem,sed_space)
    type(SED_space_0)              :: sed_space
    type(Ensemble_data),intent(in) :: ensem
    integer(MYIP) :: freenumber,i
    real(MYPRE) :: mass
    !
    freenumber = size(ensem%info%Mass%vector)
    do i = 1,freenumber
      call Set_xin_SED_0(ensem,sed_space,i)
      mass = ensem%info%Mass%vector(i)
      call Calculate_SED_0_Kinetic_Energy(sed_space,mass)
    end do 
    !
    return 
  end subroutine Calculate_SED_0_Ensemble_energy
  !
  subroutine Write_Result_SED_0(sed_space,fd)
    type(SED_space_0),intent(in) :: sed_space
    integer(MYIP) , intent(in) :: fd
    integer(MYIP)  :: wlen
    integer(MYIP) :: it
    !
    wlen = size(sed_space%sed,1)
    do it = 0,ubound(sed_space%sed,1)
      write(fd,fmt='(I8, A)',advance='no') it,","
      write(fd,fmt='(E20.7)')(sed_space%sed(it)%sum)/(DBLE(wlen)*DBLE(wlen))
    end do 
    !
    return 
  end subroutine 
  !
  !
  subroutine allocate_SED_0_group(sed0_group,start,fin)
    type(SED_space_0),allocatable :: sed0_group(:)
    integer(MYIP) ,intent(in) :: start,fin
    integer(MYIP) :: i
    !
    do i = lbound(sed0_group,1),ubound(sed0_group,1)
      call allocate_SED_0(sed0_group(i),start,fin)
    end do 
    !
    return 
  end subroutine allocate_SED_0_group
  !
  subroutine Write_Result_SED_0_group(sed0_group,fd)
    implicit none
    type(SED_space_0),allocatable,intent(in) :: sed0_group(:)
    integer(MYIP),intent(in) :: fd
    integer(MYIP) it,i
    integer(MYIP) :: N,M
    integer(MYIP) :: wlen
    
    N = ubound(sed0_group,1)
    M = lbound(sed0_group,1)
    !
    do it = 0,ubound(sed0_group(0)%sed,1)
      wlen = size(sed0_group(0)%sed,1)
      write(fd,fmt='(I8, A)',advance='no') it,","
      do i = M,N-1
        write(fd,fmt='(E20.7,A)',advance='no')sed0_group(i)%sed(it)%sum/(DBLE(wlen)*DBLE(wlen)),","
      end do 
      !
      write(fd,fmt='(E20.7)')(sed0_group(i)%sed(it)%sum/(DBLE(wlen)*DBLE(wlen)))
      !
    end do 
    !
    return 
  end subroutine Write_Result_SED_0_group 
  !
  !
  !====================================================================================================
  !                                         SED 1                                                     =
  !====================================================================================================
  !
  subroutine allocate_SED_1(sed_space,xlength,start,fin)
    integer(MYIP),intent(in) :: xlength,start,fin
    type(SED_space_1) ::  sed_space
    integer(MYIP) :: it
    integer(MYIP) :: il
    !
    !
    allocate(sed_space%xin(0:(xlength-1),0:(fin-start)))
    allocate(sed_space%sed(0:(xlength-1),0:(fin-start)))
    !
    sed_space%xin(:,:) = CMPLX(0d0,0d0,kind=kind(0d0))
    !
    do it = 0,fin-start
    do il = 0,xlength-1
      call kahan_zero(sed_space%sed(il,it))
    end do 
    end do 
    !
    sed_space%start = start
    sed_space%fin = fin
    sed_space%xlength = xlength
    !
    return 
  end subroutine allocate_SED_1
  !
  subroutine Set_xin_SED_1(sed_space,ensem,pos,freenumber)
    type(SED_space_1) :: sed_space
    type(Ensemble_data),intent(in) :: ensem
    integer(MYIP),intent(in) :: pos,freenumber
    integer(MYIP) :: it
    !
    !
    do it = sed_space%start,sed_space%fin
      sed_space%xin(pos,it-sed_space%start) = ensem%state%Velocity(it)%vector(freenumber)
    end do 
    !
    return 
  end subroutine  Set_xin_SED_1
  !
  !
  subroutine Calculate_SED_1(sed_space)
    type(SED_space_1) :: sed_space
    real(MYPRE) :: tmp
    integer(MYIP) :: it, i
    !
    call fftw_2d(sed_space%xin,sed_space%xout,FFTW_FORWARD)
    !
    do it =0,ubound(sed_space%sed,2) 
    do i =0,ubound(sed_space%sed,1) 
      tmp = abs(sed_space%xout(i,it))
      tmp = tmp * tmp
      call add_kahan(sed_space%sed(i,it),tmp) 
    end do 
    end do 
    !
    return 
  end subroutine Calculate_SED_1
  !
  subroutine Calculate_SED_1_Kinetic_Energy(sed_space,mass)
    type(SED_space_1) :: sed_space
    real(MYPRE),intent(in) :: mass
    real(MYPRE) :: tmp
    integer(MYIP) :: it, i
    !
    call fftw_2d(sed_space%xin,sed_space%xout,FFTW_FORWARD)
    !
    do it =0,ubound(sed_space%sed,2) 
    do i =0,ubound(sed_space%sed,1) 
      tmp = abs(sed_space%xout(i,it))
      tmp = 0.5d0 * mass * tmp * tmp
      call add_kahan(sed_space%sed(i,it),tmp) 
    end do 
    end do 
    !
    return
  end subroutine Calculate_SED_1_Kinetic_Energy
  !
  subroutine Write_Result_SED_1(sed_space,fd)
    type(SED_space_1),intent(in) :: sed_space
    integer(MYIP) , intent(in) :: fd
    integer(MYIP) :: xlen,wlen
    integer(MYIP) :: it,ix
    !
    xlen = ubound(sed_space%sed,1)+1
    wlen = ubound(sed_space%sed,2)+1
    print *, "xlen",xlen
    print *, "xlen",wlen
    !
    do it = 0,ubound(sed_space%sed,2)
        write(fd,fmt='(I8,A)',advance='no')it ,", "
      do ix = 0,ubound(sed_space%sed,1)-1
        write(fd,fmt='(E20.7,A)',advance='no')(sed_space%sed(ix,it)%sum)/(DBLE(xlen*wlen)*DBLE(xlen*wlen)) ,", "
      end do 
      ix = ubound(sed_space%sed,1)
      write(fd,fmt='(E20.7)')(sed_space%sed(ix,it)%sum)/(DBLE(xlen*wlen)*DBLE(xlen*wlen))   
    end do 
    !
  end subroutine 
  !
  !
  !
  !=========================================================================================!
  !                              NMP MODULE IMPLEMENT                                       !
  !=========================================================================================!
  !
  ! YET
end module ENSEMBLE_DATA_MODULE
