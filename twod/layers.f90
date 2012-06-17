
module layers
  use global_com,only:pid,eps0d,rmu0d,c1,dp,wod
  use global_geom,only:rho_max,z_min,z_max,estimated_edge_av
  implicit none

  ! multilayered media parameters for IC application, the media is isotropic and mur=1
  real(kind=dp),parameter::euler=0.577215664901533d0
  logical::is_multilayer
  integer::nlayers
  ! number of layers
  real(kind=dp),dimension(:),allocatable::h_of_layer,zlow_of_layer
  real(kind=dp)::green_x_max_pos = -1.d0, green_x_min_pos = 1.d0
  ! h_of_layer: height of each layer                                                                 
  ! zlow_of_layer: z coordinate of the interface of each layer
  real(kind=dp),dimension(:),allocatable::eps_t
  ! eps_t: transverse relative permitivity of each layer                                         
  real(kind=dp)::eps_t_max,k0,k_prop
  real(kind=dp),dimension(:),allocatable::k_prop2
  complex(kind=dp),dimension(:),allocatable::kz_wave
  ! maximum eps
  real(kind=dp),dimension(:),allocatable::Z_0
  real(kind=dp),dimension(:),allocatable::GammaL_mn,GammaR_mn       

  ! GammaL_mn, GammaR_mn, refer to (23) in the document  
  ! Characteristic impdance of each section
  real(kind=dp)::height(5)
  real(kind=dp)::ejkh(5)

  ! Layered geometry parameter
  integer::zlayer_min,zlayer_max,green_index,green_mode
  integer::layer_me,layer_ne ! store the layer number of patches
  integer::nlayers_eff ! number of layers containing conductors
  logical,allocatable::exist_cond(:) ! store whether each layer contains conductors
  integer,allocatable::layers_eff(:),map_layer(:)

  ! Sommerfeld integral parameters
  real(kind=dp)::h1,h2
  ! contour parameter
  real(kind=dp)::threshold
  ! loop stop criterion
  real(kind=dp),dimension(10)::qp,wght
  complex(kind=dp)::ft,fh
  complex(kind=dp)::kernel,kernel_t,kernel_h
  
  ! position parameters
  ! TODO: remove rs,ro globals
  real(kind=dp),dimension(2)::rs,ro
  ! rs: coordinate of source point
  ! ro: coordinate of observation point 
  integer::layer_s,layer_o
  real(kind=dp)::delta_x
  ! rho: transverse distance between source and observation points

  ! Transmission line Green's function
  real(kind=dp),dimension(4)::C_s
  ! refer to (20) in the document
  ! C_s=C_vi
  ! C_vi={1,1,1,1}
  complex(kind=dp)::TLGF,TLGF_t,TLGF_h

  ! Singularity subtraction
  real(kind=dp)::d_sub
  real(kind=dp),dimension(4)::R_n_sub
  complex(kind=dp),dimension(4)::ejkr_ns_sub
  complex(kind=dp)::TauL_vv_prod_sub,TauR_vv_prod_sub,TauL_vv_prod_coef,TauR_vv_prod_coef
  ! TauL_vv_prod/right: product of the TauL_vv or TauR_vv of all the sections between the source section and observation section, refer to (65) and (67) in the document
  complex(kind=dp)::TLGF_sub,TLGF_sub_t,TLGF_sub_h

  ! Interpolation
  real(kind=dp),allocatable::drho(:,:),dz(:)
  real(kind=dp),allocatable::drho_ij(:,:,:),dz_ij(:,:)
  integer,allocatable::num_rho(:,:),num_z(:)
  integer::layer_max,layer_min,layer_num
  type samelayer
     complex(kind=dp),pointer::Gf_grid_array_t(:,:),Gf_grid_array_h(:,:) ! for the same layer
  end type samelayer
  type(samelayer),allocatable::gf_table_same(:)
  type difflayer
     complex(kind=dp),pointer::Gf_grid_array(:,:,:) ! for the different layer
  end type difflayer
  type(difflayer),allocatable::gf_table_diff(:,:)
    
  complex(kind=dp),allocatable::green_array(:,:)
  real(kind=dp),allocatable::src_obs_array(:,:)

  save
  contains
    subroutine inlayers
      use global_com,only:real_mem,complex_mem
      implicit none
      integer::i
      real(kind=dp)::mem_est
      ! read multiplayered media parameters
      read(66,*) is_multilayer
      if (.not. is_multilayer) then  
         nlayers=1
      else
         read(66,*) nlayers
      end if
         
      if (nlayers==1) then
         !print*,'The layered medium is only free space'
      else if (nlayers<1) then
         !print*,'ERROR::NLAYERS MUST BE A POSITIVE INTEGER GREAT THAN 1!'
         stop
      else
         !print*,'The number of layered medium is: ',nlayers
      end if

      mem_est=(4.d0*nlayers+1)*real_mem+(5.d0*nlayers)*complex_mem
      !print*,'The layered medium requires memory (MB): ',mem_est/1024.d0/1024.d0
      ! Essential array
      allocate(h_of_layer(1:nlayers),zlow_of_layer(1:nlayers+1),eps_t(1:nlayers),&
           Z_0(1:nlayers),GammaL_mn(1:nlayers),GammaR_mn(1:nlayers),&
           kz_wave(1:nlayers),k_prop2(1:nlayers))
      
      if (.not. is_multilayer) then           
         eps_t(1)=1.d0; h_of_layer(1)=-1.d0
         zlow_of_layer(1)=-1.d0; zlow_of_layer(2)=-1.d0
      else
         read(66,*) zlow_of_layer(2)
         ! In order to handle the infinity case, we use the following notations.
         ! h_of_layer(i)==-1.0 means the height is infinite
         ! In order to handle PEC case, we use the following notation: eps_t=-1 (Inf), Z_0=0
         do i=1,nlayers
            read(66,*) eps_t(i),h_of_layer(i)
            if (i==1) then
               if (h_of_layer(i)==-1.d0) then           
                  zlow_of_layer(i)=-1.d0 ! the height of the bottom layer is infinite (half space)
               else
                  zlow_of_layer(i)=zlow_of_layer(2)  ! the height of the bottom layer is 0 (PEC or PMC)
               end if
            else if (i==nlayers) then
               if (h_of_layer(i)==-1.d0) then           
                  zlow_of_layer(i+1)=-1.d0 ! the height of the top layer is infinite
               else
                  zlow_of_layer(i+1)=zlow_of_layer(i-1)+h_of_layer(i-1)  ! the height of the top layer is 0 (PEC or PMC)
               end if
            else
               ! when layer_s==nlayers, zlow_of_layer(nlayers+1) should be handled carefully
               zlow_of_layer(i+1)=zlow_of_layer(i)+h_of_layer(i) 
            end if
         end do
         read(66,*) threshold     ! Sommerfeld integral stop criterion
      end if

      eps_t_max=maxval(eps_t(:)) 
      C_s(1:4)=(/1.d0,1.d0,1.d0,1.d0/)
      return
    end subroutine inlayers


    subroutine init_layers
      implicit none

      integer::i,j

      ! Initialization
      Z_0(:)=0.d0
      ! Calculate the propogation constant and characteristic impedance of each section
      ! PEC should be handled individually
      do j=1,nlayers
         if (eps_t(j)/=-1.d0) then ! j-th layer is dielectric
            Z_0(j)=1.d0/eps_t(j) ! characteristic impedance
            k_prop2(j)=(wod*dsqrt(eps0d*rmu0d*eps_t(j)))**2
         else ! PEC
            Z_0(j)=0.d0 ! characteristic impedance
            k_prop2(j)=0.d0
         end if
      end do

      k0=wod*dsqrt(eps0d*rmu0d)
      h1=5.d-3*k0
      h2=1.2_dp*k0*dsqrt(eps_t_max)

      if (.not. is_multilayer) then
         GammaL_mn(1)=0.d0  ! GammaL and GammaL_mn are equal to 0 in the first section
         GammaR_mn(nlayers)=0.d0 ! GammaR and GammaR_mn are equal to 0 in the last section
      else
         GammaL_mn(1)=0.d0  ! GammaL and GammaL_mn are equal to 0 in the first section
         GammaR_mn(nlayers)=0.d0 ! GammaR and GammaR_mn are equal to 0 in the last section
         
         GammaL_mn(2)=(Z_0(1)-Z_0(2))/(Z_0(1)+Z_0(2))
         GammaR_mn(nlayers-1)=(Z_0(nlayers)-Z_0(nlayers-1))/(Z_0(nlayers)+Z_0(nlayers-1))
         do i=2,nlayers-1
            GammaL_mn(i+1)=(Z_0(i)-Z_0(i+1))/(Z_0(i)+Z_0(i+1))
            GammaR_mn(nlayers-i)=(Z_0(nlayers-i+1)-Z_0(nlayers-i))/(Z_0(nlayers-i+1)+Z_0(nlayers-i)) 
         end do
      end if

      return 
    end subroutine init_layers

    subroutine find_rho_z
      use global_com,only:max_INF
      use global_geom,only:nsuinf,sunod, nsuedgn
      implicit none

      integer::ne,me
      real(kind=dp)::r1(2),r2(2),dist
      
      ! find effective layer
      !allocate(exist_cond(1:nlayers))
      !exist_cond(:)=.false.
      !do ne=1,nsuinf(1)
      !   r1=sunod(:,ne)
      !   call find_layer(r1(1:2),layer_ne)
      !   exist_cond(layer_ne)=.true.
      !end do

      nlayers_eff=0
      do ne=1,nlayers
         if (exist_cond(ne)) then
            nlayers_eff=nlayers_eff+1
         end if
      end do

      allocate(layers_eff(1:nlayers_eff))
      allocate(map_layer(1:nlayers))
      allocate(z_min(1:nlayers_eff),z_max(1:nlayers_eff))
      allocate(rho_max(1:nlayers_eff,1:nlayers_eff))  
      nlayers_eff=0
      map_layer(:)=0
      z_max(:)=-1.d0
      z_min(:)=-1.d0
      rho_max(:,:)=-1.d0
      do ne=1,nlayers
         if (exist_cond(ne)) then
            nlayers_eff=nlayers_eff+1
            layers_eff(nlayers_eff)=ne
            map_layer(ne)=nlayers_eff
            if (ne==nlayers) then
               z_max(map_layer(ne))=1.5d0*zlow_of_layer(ne) ! This will not happen for Altan's cases
            else
               z_max(map_layer(ne))=zlow_of_layer(ne+1)-1.d-3*h_of_layer(ne) ! need to be validated
            end if
            z_min(map_layer(ne))=zlow_of_layer(ne)
         end if
      end do

      !print*,'The number of layers containing conductors is: ',nlayers_eff
      !print*,'They are: ',layers_eff(1:nlayers_eff)

      ! find the range of layer the structure occupies and boundary of the structure in rho and z direction
      ! This is O(N^2/P) operation
      ! find z_max and z_min 

      do ne=1,nlayers
         if (exist_cond(ne)) then
            do me=1,nlayers
               if (exist_cond(me)) then
                  rho_max(map_layer(ne),map_layer(me))=green_x_max_pos-green_x_min_pos
               end if
            end do
         end if
      end do
    
      return
    end subroutine find_rho_z

    subroutine find_layer(qp,layer_qp)
      implicit none

      real(kind=dp),intent(in)::qp(2)
      integer,intent(out)::layer_qp
      integer::i

      if (is_multilayer) then
         if (qp(2)<zlow_of_layer(2)) then
            layer_qp=1
         elseif (qp(2)>=zlow_of_layer(nlayers)) then
            layer_qp=nlayers
         else
            i=2
            do  ! if the qp is on the interface, it belongs to the upper layer
               if (qp(2)<zlow_of_layer(i+1)) then
                  layer_qp=i
                  exit
               end if
               i=i+1
            end do
         end if
      else
         layer_qp=1
      end if
      
      return
    end subroutine find_layer

    subroutine init_interpolation
      implicit none

      integer::i,j,ne,me

      ! The interpolation spacing is determined by the average edge length

      ! find drho and drho is needed for any structure
      ! interpolated interval is [0,rho_max+4*dz]
      allocate(num_rho(1:nlayers_eff,1:nlayers_eff),num_z(1:nlayers_eff))
      allocate(drho(1:nlayers_eff,1:nlayers_eff),dz(1:nlayers_eff))
      drho(:,:)=2.d0*estimated_edge_av
      !print*,"drho is: ",drho(1,1)

      ! +6 should be changed to +3 or +1 to save more
      do i=1,nlayers_eff
         do j=1,nlayers_eff
            num_rho(j,i)=rho_max(j,i)/drho(j,i)+4
            if (num_rho(j,i)<6) then
               num_rho(j,i)=6
            end if
            !print*,"For layers",layers_eff(j),layers_eff(i),"max_rho is: ",rho_max(j,i)
            !print*,"For layers",layers_eff(j),layers_eff(i),"num_rho is: ",num_rho(j,i)
         end do
      end do

      ! find dz and interpolated interval [z_min,z_max]
      do i=1,nlayers_eff
         if (z_max(i)==z_min(i)) then
            num_z(i)=1
         else
            dz(i)=2.d0*min(z_max(i)-z_min(i),estimated_edge_av)
            num_z(i)=(z_max(i)-z_min(i))/dz(i)
            if (num_z(i)<6) then
               num_z(i)=6
            end if
            dz(i)=(z_max(i)-z_min(i))/num_z(i)
         end if
         !print*,"For layer",layers_eff(i),"min z is: ",z_min(i)
         !print*,"For layer",layers_eff(i),"max_z is: ",z_max(i)
         !print*,"For layer",layers_eff(i),"num_z is: ",num_z(i)
      end do

      allocate(gf_table_same(1:nlayers_eff))
      allocate(gf_table_diff(1:nlayers_eff,1:nlayers_eff))
      do i=1,nlayers_eff
         allocate(gf_table_same(i)%Gf_grid_array_t(0:num_rho(i,i),0:num_z(i)))
         allocate(gf_table_same(i)%Gf_grid_array_h(0:num_rho(i,i),0:2*num_z(i)+1))
         gf_table_same(i)%Gf_grid_array_t(0:num_rho(i,i),0:num_z(i))=&
              cmplx(0.d0,0.d0,dp)
         gf_table_same(i)%Gf_grid_array_h(0:num_rho(i,i),0:2*num_z(i)+1)=&
              cmplx(0.d0,0.d0,dp)
         do j=1,nlayers_eff
            if (i/=j) then
               allocate(gf_table_diff(j,i)%Gf_grid_array(0:num_rho(j,i),0:num_z(j),0:num_z(i)))
               gf_table_diff(j,i)%Gf_grid_array(0:num_rho(j,i),0:num_z(j),0:num_z(i))=&
                    cmplx(0.d0,0.d0,dp)
            end if
         end do
      end do

      allocate(drho_ij(1:6,1:nlayers_eff,1:nlayers_eff),dz_ij(1:6,1:nlayers_eff))
      drho_ij(1:6,1:nlayers_eff,1:nlayers_eff)=1.d0
      dz_ij(1:6,1:nlayers_eff)=1.d0

      do ne=1,nlayers_eff
         do me=1,nlayers_eff
            do i=1,6
               do j=1,6
                  if (i/=j) then
                     drho_ij(i,me,ne)=drho_ij(i,me,ne)/drho(me,ne)/(i-j)
                  end if
               end do
            end do
         end do
      end do

      do ne=1,nlayers_eff
         do i=1,6
            do j=1,6
               if (i/=j) then
                  dz_ij(i,ne)=dz_ij(i,ne)/dz(ne)/(i-j)
               end if
            end do
         end do
      end do
      return
    end subroutine init_interpolation

    subroutine fill_Green_stored_array
      use global_com,only:real_mem
      use global_geom,only:sunod
      implicit none

      complex(kind=dp)::Gf_t,Gf,Gf_tmp,Gf_tmp2,Gf_tmp3,Gf_tmp4,Gf_tmp5,Gf_sub
      real(kind=dp),dimension(2)::rs,ro
      integer::irho,iz,jz
      integer::counter,i,j
      real(kind=dp)::r1(2),r2(2),rho
      real(kind=dp)::mem_est
      
      
      !print*,'Fill self term Gf table'
      
      do i=1,nlayers_eff
         ! Toeplitz and Hankel [0,h]
         rs(:)=(/green_x_min_pos,z_min(i)/)
         call find_layer(rs(:),layer_s)      
      
         counter=0
         do iz=1,num_z(i)+1
            do irho=1,num_rho(i,i)+1
               counter=counter+1
               ro(:)=rs(:)+(/(irho-1)*drho(i,i),(iz-1)*dz(i)/)
               call find_layer(ro(:),layer_o)
               delta_x=ro(1)-rs(1)
               call find_height(rs,ro)
            
               call fill_Layered_Green(rs,ro,Gf,Gf_tmp,Gf_tmp2,Gf_tmp3,Gf_tmp4,Gf_tmp5)
               gf_table_same(i)%Gf_grid_array_t(irho-1,iz-1)=Gf_tmp2
               gf_table_same(i)%Gf_grid_array_h(irho-1,iz-1)=Gf_tmp3

               if (green_mode==2) then
                  if (modulo(counter-1,50)==0) then
                     !print*,'Layer',layers_eff(i),'TH',counter,'of',(num_z(i)+1)*(num_rho(i,i)+1)
                  end if
               end if
            end do
         end do
         
         ! Hankel [h,2h]
         rs(:)=(/green_x_min_pos,z_max(i)/)
         call find_layer(rs(:),layer_s)  
      
         counter=0
         do iz=1,num_z(i)+1
            do irho=1,num_rho(i,i)+1
               counter=counter+1
               ro(:)=(/rs(1),z_min(i)/)+(/(irho-1)*drho(i,i),(iz-1)*dz(i)/)
               call find_layer(ro(:),layer_o)
               delta_x=ro(1)-rs(1)
               call find_height(rs,ro)
               
               call fill_Layered_Green(rs,ro,Gf,Gf_tmp,Gf_tmp2,Gf_tmp3,Gf_tmp4,Gf_tmp5)
               gf_table_same(i)%Gf_grid_array_h(irho-1,iz+num_z(i))=Gf_tmp3

               if (green_mode==2) then
                  if (modulo(counter-1,50)==0) then
                     !print*,'Layer',layers_eff(i),'H',counter,'of',(num_z(i)+1)*(num_rho(i,i)+1)
                  end if
               end if
            end do
         end do
      end do

      if (nlayers_eff/=1) then
         !print*,'Fill mutual term Gf table'
      
         ! Fill Gf layer_s<layer_o
         do i=1,nlayers_eff ! src
            r1(:)=(/green_x_min_pos,z_min(i)/)
            do j=i+1,nlayers_eff ! obs
               r2(:)=(/green_x_min_pos,z_min(j)/)
   
               counter=0
               do iz=1,num_z(i)+1 ! soruce
                  do jz=1,num_z(j)+1 ! observer
                     do irho=1,num_rho(j,i)+1
                        counter=counter+1
                        rs(:)=r1(:)+(/0.d0,(iz-1)*dz(i)/)
                        call find_layer(rs(:),layer_s)
                        ro=r2(:)+(/(irho-1)*drho(j,i),(jz-1)*dz(j)/)
                        call find_layer(ro(:),layer_o)
                        delta_x=ro(1)-rs(1)
                        call find_height(rs,ro)

                        call fill_Layered_Green(rs,ro,Gf,Gf_tmp,Gf_tmp2,Gf_tmp3,Gf_tmp4,Gf_tmp5)
                        gf_table_diff(j,i)%Gf_grid_array(irho-1,jz-1,iz-1)=Gf_tmp ! rho,z,z'
                        gf_table_diff(i,j)%Gf_grid_array(irho-1,iz-1,jz-1)=Gf

                        if (green_mode==2) then
                           if (modulo(counter-1,50)==0) then 
                              !print*,'Layer',layers_eff(j),layers_eff(i),counter,'of',(num_z(i)+1)*(num_z(j)+1)*(num_rho(j,i)+1)
                           end if
                        end if
                     end do
                  end do
               end do
            end do
         end do
         
         ! Fill Gf layer_s>layer_o using reciprocity
         do j=1,nlayers_eff ! obs
            r2(:)=(/green_x_min_pos,z_min(j)/)
            !r2(:)=sunod(:,1); r2(2)=z_min(j)
            do i=j+1,nlayers_eff ! src
               !r1(:)=sunod(:,1); r1(2)=z_min(i)
               r1(:)=(/green_x_min_pos,z_min(i)/)
               counter=0
               do jz=1,num_z(j)+1 ! observer
                  do iz=1,num_z(i)+1 ! source
                     do irho=1,num_rho(i,j)+1
                        counter=counter+1
                        rs(:)=r1(:)+(/0.d0,(iz-1)*dz(i)/)
                        call find_layer(rs(:),layer_s)
                        ro=r2(:)+(/(irho-1)*drho(i,j),(jz-1)*dz(j)/)
                        call find_layer(ro(:),layer_o)
                        delta_x=ro(1)-rs(1)
                        call find_height(rs,ro)
                        
                        ! direct calculation
!                        call fill_Layered_Green(Gf,Gf_tmp,Gf_tmp2,Gf_tmp3,Gf_tmp4,Gf_tmp5)
!                        gf_table_diff(j,i)%Gf_grid_array(irho-1,jz-1,iz-1)=Gf_tmp
                        
                        green_index = green_index + 1
                        call find_subtraction(rs,ro,Gf_sub)
                        gf_table_diff(j,i)%Gf_grid_array(irho-1,jz-1,iz-1)=&
                             gf_table_diff(j,i)%Gf_grid_array(irho-1,jz-1,iz-1)-Gf_sub
                        
                        if (green_mode==2) then
                           if (modulo(counter-1,50)==0) then
                              !print*,'Layer',layers_eff(j),layers_eff(i),counter,'of',(num_z(i)+1)*(num_z(j)+1)*(num_rho(i,j)+1)
                           end if
                        end if
                     end do
                  end do
               end do
            end do
         end do
      end if
      return
    end subroutine fill_Green_stored_array

    subroutine fill_Layered_Green(src, obs, Gf,Gf_nsigu,Gf_t_nsigu,Gf_h_nsigu,Gf_t,Gf_h)
      implicit none
      real(kind=dp),dimension(2),intent(in)::src,obs
      complex(kind=dp),intent(out)::Gf,Gf_t,Gf_h,Gf_nsigu,Gf_t_nsigu,Gf_h_nsigu
      ! fi (i=1,...,4) are the different integrand
      complex(kind=dp),external::f1_2d,f2_2d,f3_2d,f4_2d
      ! Define numericial integral parameters
      real(kind=dp)::a,b,eps(2),ratio,step,start
      integer::i,j,tmp,ier
      integer::qrule,limit,num_int
      complex(kind=dp)::k_rho,Gf_tmp,Gf_tmp_t,Gf_tmp_h,num,num_t,num_h
      ! extraction part
      real(kind=dp)::dist1,dist2,rho_array(0:9),z_array(0:9),dint
      complex(kind=dp)::h02(0:9)
      complex(kind=dp)::Gf_sub,Gf_sub_t,Gf_sub_h
      complex(kind=dp)::num_sta,num_sta_t,num_sta_h

      ! green_mode == 0 is to count how many entries in green_array
      ! green_mode == 2 is to calculate and store in green_array
      ! green_mode == 1 is for lookup from green_array
      if (green_mode == 0) then
         green_index = green_index + 1
         return
      else if (green_mode == 1) then
         Gf = green_array(1,green_index)
         Gf_nsigu=green_array(2,green_index)
         Gf_t_nsigu=green_array(3,green_index)
         Gf_h_nsigu=green_array(4,green_index)
         Gf_t=green_array(5,green_index)
         Gf_h=green_array(6,green_index)
         !(/Gf,Gf_nsigu,Gf_t_nsigu,Gf_h_nsigu,Gf_t,Gf_h/) = green_array(:,green_index)
         green_index = green_index + 1
         return
      else if (green_mode == 3) then
         src_obs_array(:,green_index) = (/src(1),src(2),obs(1),obs(2),0.d0,0.d0/)
         !print *, 'GGG', green_index, src_obs_array(:,green_index)
         !if (green_index > 250) then
         !   call exit
         !end if
         green_index = green_index + 1
         return
      else if (green_mode == 4) then
         Gf = cmplx(src_obs_array(5,green_index), src_obs_array(6,green_index))
         green_index = green_index + 1
         return
      else if (green_mode == 2) then
         ! compute green function below
      else
         print *, "bad green mode:", green_mode
         return
      end if
      
      call find_R_n_sub
      delta_x=obs(1)-src(1)

      Gf_tmp=cmplx(0.d0,0.d0,dp)
      Gf_tmp_t=cmplx(0.d0,0.d0,dp)
      Gf_tmp_h=cmplx(0.d0,0.d0,dp)
      Gf_nsigu=cmplx(0.d0,0.d0,dp)

      ! eps(1) for first 3 integrals; eps(2) for the last integral
      eps(1)=1.d-3
      eps(2)=1.d-6

      rs = src
      ro = obs

      qrule=7; num_int=1; dint=1.d0/num_int
      call oneD_quadrature(qrule,qp,wght)
      
      ! Contour integral along C1
      a=0.d0;b=h1
      call quadrature(f1_2d,a,b,qrule,num_sta,num_sta_t,num_sta_h)
      call quadrature_1D(f1_2d,a,b,qrule,eps(1),limit,num_sta,num,num_t,num_h,ier)    

      Gf_tmp=Gf_tmp+num
      Gf_tmp_t=Gf_tmp_t+num_t
      Gf_tmp_h=Gf_tmp_h+num_h

      ! Contour integral along C2
      do j=1,num_int
         a=(j-1)*dint*h2;b=j*dint*h2
         call quadrature(f2_2d,a,b,qrule,num_sta,num_sta_t,num_sta_h)
         call quadrature_1D(f2_2d,a,b,qrule,eps(1),limit,num_sta,num,num_t,num_h,ier)    

         Gf_tmp=Gf_tmp+num
         Gf_tmp_t=Gf_tmp_t+num_t
         Gf_tmp_h=Gf_tmp_h+num_h
      end do

      ! Contour integral along C3
      a=h1;b=0.d0
      call quadrature(f3_2d,a,b,qrule,num_sta,num_sta_t,num_sta_h)
      call quadrature_1D(f3_2d,a,b,qrule,eps(1),limit,num_sta,num,num_t,num_h,ier)    

      Gf_tmp=Gf_tmp+num
      Gf_tmp_t=Gf_tmp_t+num_t
      Gf_tmp_h=Gf_tmp_h+num_h

      ! As there is no pole on real krho axis, the contour is chosen the same as real krho axis
      if (dabs(delta_x)<1.d-5) then
         step=5.d4
      else if (dabs(delta_x)>=1.d-5 .and. dabs(delta_x)<5.d-3) then
         step=5.d3
      else if (dabs(delta_x)>=5.d-3 .and. dabs(delta_x)<5.d-2) then
         step=2.d3
      else
         step=5.d2
      end if

      ! Remainder Contour integral along on real axis
      j=0
      do
         a=h2+j*step;b=h2+(j+1)*step
         call quadrature(f4_2d,a,b,qrule,num_sta,num_sta_t,num_sta_h)
         call quadrature_1D(f4_2d,a,b,qrule,eps(2),limit,num_sta,num,num_t,num_h,ier)
         
         Gf_tmp=Gf_tmp+num
         Gf_tmp_t=Gf_tmp_t+num_t
         Gf_tmp_h=Gf_tmp_h+num_h
         
!         ratio=cdabs(num/Gf_tmp)
!         print*,'ratio',j,ratio
         if (cdabs(num)<=cdabs(Gf_tmp)*threshold) then
            exit
         else
            j=j+1
         end if
      end do

      if (isnan(cdabs(Gf_tmp))) then
         !print*,'Gf is NAN'
         !print*,'src',src
         !print*,'obs',obs
         stop
      end if
         
      k_prop=k0*sqrt(eps_t(layer_s))
      ! All the ri here are actually 1/ri
      ! As the triangle is sub-divided, no exact singularity exists. But the integral might be inaccurate. 
      if (layer_s==layer_o) then
         z_array(0)=obs(2)-src(2)
         z_array(1)=2.d0*zlow_of_layer(layer_s+1)-src(2)-obs(2)
         z_array(2)=src(2)+obs(2)-2.d0*zlow_of_layer(layer_s)
         z_array(3)=2.d0*h_of_layer(layer_s)-src(2)+obs(2)
         z_array(4)=2.d0*h_of_layer(layer_s)+src(2)-obs(2)
         ! The singular term for EFIE is handled seperatly
         rho_array(0:4)=dsqrt(delta_x**2+z_array(0:4)**2)

         h02(0:4)=(1.d0-c1*2.d0/pid*(dlog(0.5d0*k_prop*rho_array(0:4))+euler))

         if (zlow_of_layer(layer_s+1)==-1.d0) then
            h02(1)=cmplx(0.d0,0.d0,dp);h02(3:4)=cmplx(0.d0,0.d0,dp)
         else if (zlow_of_layer(layer_s)==-1.d0) then
            h02(2)=cmplx(0.d0,0.d0,dp);h02(3:4)=cmplx(0.d0,0.d0,dp)
         end if
         
         Gf_sub_t=h02(0)&
              +GammaR_mn(layer_s)*GammaL_mn(layer_s)*h02(3)&
              +GammaR_mn(layer_s)*GammaL_mn(layer_s)*h02(4)
            
         if (rho_array(1)==0.d0 .or. rho_array(2)==0.d0) then
            Gf_sub_h=cmplx(0.d0,0.d0,dp)
         else
            Gf_sub_h=GammaR_mn(layer_s)*h02(1)&
                 +GammaL_mn(layer_s)*h02(2)
         end if
         Gf_sub_t=0.5d0*Z_0(layer_s)*Gf_sub_t/(2.d0*c1)*pid ! pid will be divided later
         Gf_sub_h=0.5d0*Z_0(layer_s)*Gf_sub_h/(2.d0*c1)*pid
         Gf_sub=Gf_sub_t+Gf_sub_h
         
      else if (layer_s<layer_o) then
         ! zlow_of_layer(layer_s+1) and zlow_of_layer(layer_o) cannot be infinity
         z_array(0)=zlow_of_layer(layer_s+1)-src(2)+obs(2)-zlow_of_layer(layer_o)+d_sub
         z_array(1)=zlow_of_layer(layer_s+1)-src(2)+obs(2)-zlow_of_layer(layer_o)+d_sub
         z_array(2)=h_of_layer(layer_s)-zlow_of_layer(layer_s)+src(2)+obs(2)-zlow_of_layer(layer_o)+d_sub
         z_array(3)=2.d0*h_of_layer(layer_s)+zlow_of_layer(layer_s+1)-src(2)+obs(2)-zlow_of_layer(layer_o)+d_sub
         z_array(4)=2.d0*h_of_layer(layer_s)-zlow_of_layer(layer_s+1)+src(2)+obs(2)-zlow_of_layer(layer_o)+d_sub
         z_array(5)=zlow_of_layer(layer_s+1)-src(2)+2.d0*zlow_of_layer(layer_o+1)-obs(2)-zlow_of_layer(layer_o)+d_sub
         z_array(6)=2.d0*zlow_of_layer(layer_s+1)-zlow_of_layer(layer_s+1)-src(2)+2.d0&
              *zlow_of_layer(layer_o+1)-obs(2)-zlow_of_layer(layer_o)+d_sub
         z_array(7)=zlow_of_layer(layer_s+1)-2.d0*zlow_of_layer(layer_s)+src(2)+2.d0&
              *zlow_of_layer(layer_o+1)-obs(2)-zlow_of_layer(layer_o)+d_sub
         z_array(8)=2.d0*h_of_layer(layer_s)+zlow_of_layer(layer_s+1)-src(2)+2.d0&
              *zlow_of_layer(layer_o+1)-obs(2)-zlow_of_layer(layer_o)+d_sub
         z_array(9)=2.d0*h_of_layer(layer_s)-zlow_of_layer(layer_s+1)+src(2)+2.d0&
              *zlow_of_layer(layer_o+1)-obs(2)-zlow_of_layer(layer_o)+d_sub
         
         rho_array(0:9)=dsqrt(delta_x**2+z_array(0:9)**2)
         h02(0:9)=(1.d0-c1*2.d0/pid*(dlog(0.5d0*k_prop*rho_array(0:9))+euler))
         ! Whenever zlow_of_layer(layer_s+1/layer_s)=-1, h_of_layer(layer_s)==-1.d0
         if (zlow_of_layer(layer_o+1)==-1.d0) then
            h02(5:9)=cmplx(0.d0,0.d0,dp)
         end if
         
         if (zlow_of_layer(layer_s)==-1.d0) then
            h02(2:4)=cmplx(0.d0,0.d0,dp);h02(7:9)=cmplx(0.d0,0.d0,dp)
         end if
         
         Gf_sub=h02(0)&
              +GammaR_mn(layer_s)*h02(1)&
              +GammaL_mn(layer_s)*h02(2)&
              +GammaR_mn(layer_s)*GammaL_mn(layer_s)*h02(3)&
              +GammaR_mn(layer_s)*GammaL_mn(layer_s)*h02(4)& ! first 5 terms
              
              +GammaR_mn(layer_o)*h02(5)&
              +GammaR_mn(layer_s)*GammaR_mn(layer_o)*h02(6)&
              +GammaL_mn(layer_s)*GammaR_mn(layer_o)*h02(7)&
              +GammaR_mn(layer_s)*GammaL_mn(layer_s)*GammaR_mn(layer_o)*h02(8)&
              +GammaR_mn(layer_s)*GammaL_mn(layer_s)*GammaR_mn(layer_o)*h02(9)
         
         Gf_sub=0.5d0*Z_0(layer_s)*Gf_sub*TauR_vv_prod_coef/(2.d0*c1)*pid
      else
         ! zlow_of_layer(layer_s) and zlow_of_layer(layer_o+1) cannot be infinity
         z_array(0)=dabs(zlow_of_layer(layer_s)-src(2))-obs(2)+zlow_of_layer(layer_o+1)+d_sub
         z_array(1)=2.d0*zlow_of_layer(layer_s+1)-zlow_of_layer(layer_s)-src(2)-obs(2)+zlow_of_layer(layer_o+1)+d_sub
         z_array(2)=zlow_of_layer(layer_s)-2.d0*zlow_of_layer(layer_s)+src(2)-obs(2)+zlow_of_layer(layer_o+1)+d_sub
         z_array(3)=2.d0*h_of_layer(layer_s)+zlow_of_layer(layer_s)-src(2)-obs(2)+zlow_of_layer(layer_o+1)+d_sub
         z_array(4)=2.d0*h_of_layer(layer_s)-zlow_of_layer(layer_s)+src(2)-obs(2)+zlow_of_layer(layer_o+1)+d_sub
         z_array(5)=dabs(zlow_of_layer(layer_s)-src(2))+obs(2)+zlow_of_layer(layer_o+1)-2.d0*zlow_of_layer(layer_o)+d_sub
         z_array(6)=2.d0*zlow_of_layer(layer_s+1)-zlow_of_layer(layer_s)-src(2)+obs(2)&
              +zlow_of_layer(layer_o+1)-2.d0*zlow_of_layer(layer_o)+d_sub
         z_array(7)=zlow_of_layer(layer_s)+src(2)-2.d0*zlow_of_layer(layer_s)+obs(2)&
              +zlow_of_layer(layer_o+1)-2.d0*zlow_of_layer(layer_o)+d_sub
         z_array(8)=2.d0*h_of_layer(layer_s)+zlow_of_layer(layer_s)-src(2)+obs(2)&
              +zlow_of_layer(layer_o+1)-2.d0*zlow_of_layer(layer_o)+d_sub
         z_array(9)=2.d0*h_of_layer(layer_s)-zlow_of_layer(layer_s)+src(2)+obs(2)&
              +zlow_of_layer(layer_o+1)-2.d0*zlow_of_layer(layer_o)+d_sub
         
         rho_array(0:9)=dsqrt(delta_x**2+z_array(0:9)**2)
         h02(0:9)=(1.d0-c1*2.d0/pid*(dlog(0.5d0*k_prop*rho_array(0:9))+euler))

         if (zlow_of_layer(layer_o)==-1.d0) then
            h02(5:9)=cmplx(0.d0,0.d0,dp)
         end if
         
         if (zlow_of_layer(layer_s+1)==-1.d0) then
            h02(1)=cmplx(0.d0,0.d0,dp);h02(3:4)=cmplx(0.d0,0.d0,dp);&
                 h02(6)=cmplx(0.d0,0.d0,dp);h02(8:9)=cmplx(0.d0,0.d0,dp)
         end if
         
         Gf_sub=h02(0)&
              +GammaR_mn(layer_s)*h02(1)&
              +GammaL_mn(layer_s)*h02(2)&
              +GammaR_mn(layer_s)*GammaL_mn(layer_s)*h02(3)&
              +GammaR_mn(layer_s)*GammaL_mn(layer_s)*h02(4)& ! first 5 terms
              
              +GammaL_mn(layer_o)*h02(5)&
              +GammaR_mn(layer_s)*GammaL_mn(layer_o)*h02(6)&
              +GammaL_mn(layer_s)*GammaL_mn(layer_o)*h02(7)&
              +GammaR_mn(layer_s)*GammaL_mn(layer_s)*GammaL_mn(layer_o)*h02(8)&
              +GammaR_mn(layer_s)*GammaL_mn(layer_s)*GammaL_mn(layer_o)*h02(9)
         
         Gf_sub=0.5d0*Z_0(layer_s)*Gf_sub*TauL_vv_prod_coef/(2.d0*c1)*pid
      end if

      ! Gf: total Green's function; Gf_nsigu: total Green's function after extracting direct and potential terms
      ! Gf_nsigu=Gf_t+Gf_h
      ! Gf_t: Toeplitz part of Green's function after extracting direct and potential terms
      ! Gf_h: Hankel part of Green's function after extracting potential terms
      ! For GA_zx and GA_xz, there is no direct or potential terms. Therefore, Gf_nsigu=Gf
      ! gf_rule=6,7 are for calculating scattering field. No singularity is extracted. Gf_nsigu=Gf
      ! Gf_t and Gf_h are meaningless for 6 and 7 cases.
      
!      Gf=Gf_tmp/pid
      Gf=(Gf_tmp+Gf_sub)/pid

      Gf_nsigu=Gf_tmp/pid
      Gf_t_nsigu=Gf_tmp_t/pid
      Gf_h_nsigu=Gf_tmp_h/pid
      Gf_t=(Gf_tmp_t+Gf_sub_t)/pid
      Gf_h=(Gf_tmp_h+Gf_sub_h)/pid
      
      ! temp temp temp
      !Gf_nsigu=cmplx(0.d0,0.d0,dp)
      !Gf_t_nsigu=cmplx(0.d0,0.d0,dp)
      !Gf_h_nsigu=cmplx(0.d0,0.d0,dp)
      !Gf_h=cmplx(0.d0,0.d0,dp)

!      !print*,Gf!,Gf_sub/pid
      
      green_array(:,green_index) = (/Gf,Gf_nsigu,Gf_t_nsigu,Gf_h_nsigu,Gf_t,Gf_h/)
      !write(990,*) green_index,src,obs,Gf
      green_index = green_index + 1
      return
    end subroutine fill_Layered_Green

    subroutine find_height(src,obs)
      implicit none  
      real(kind=dp),dimension(2),intent(in)::src,obs

      if (is_multilayer) then
         if (layer_s==layer_o) then
            height(1)=2.d0*zlow_of_layer(layer_s+1)-(obs(2)+src(2))
            height(2)=(obs(2)+src(2))-2.d0*zlow_of_layer(layer_s)
            height(3)=2.d0*h_of_layer(layer_s)+(obs(2)-src(2))
            height(4)=2.d0*h_of_layer(layer_s)-(obs(2)-src(2))
            height(5)=dabs(obs(2)-src(2))
         else if (layer_s>layer_o) then
            height(1)=2.d0*zlow_of_layer(layer_s+1)-(zlow_of_layer(layer_s)+src(2))
            height(2)=(zlow_of_layer(layer_s)+src(2))-2.d0*zlow_of_layer(layer_s)
            height(3)=2.d0*h_of_layer(layer_s)+(zlow_of_layer(layer_s)-src(2))
            height(4)=2.d0*h_of_layer(layer_s)-(zlow_of_layer(layer_s)-src(2))
            height(5)=dabs(zlow_of_layer(layer_s)-src(2))
         else
            height(1)=2.d0*zlow_of_layer(layer_s+1)-(zlow_of_layer(layer_s+1)+src(2))
            height(2)=(zlow_of_layer(layer_s+1)+src(2))-2.d0*zlow_of_layer(layer_s)
            height(3)=2.d0*h_of_layer(layer_s)+(zlow_of_layer(layer_s+1)-src(2))
            height(4)=2.d0*h_of_layer(layer_s)-(zlow_of_layer(layer_s+1)-src(2))
            height(5)=dabs(zlow_of_layer(layer_s+1)-src(2))
         end if
         ejkh(1:5)=dexp(height(1:5))
      else
         height(1:4)=-1.d0
         ejkh(1:4)=0.d0
         ejkh(5)=dexp(height(5))
      end if
      return
    end subroutine find_height
    
    subroutine find_R_n_sub
      implicit none
      
      R_n_sub(1)=GammaR_mn(layer_s)
      R_n_sub(2)=GammaL_mn(layer_s)
      R_n_sub(3)=GammaL_mn(layer_s)*GammaR_mn(layer_s)
      R_n_sub(4)=GammaL_mn(layer_s)*GammaR_mn(layer_s)
      return
    end subroutine find_R_n_sub

    subroutine find_subtraction(src,obs,Gf_sub)
      implicit none
      real(kind=dp),dimension(2),intent(in)::src,obs

      complex(kind=dp),intent(out)::Gf_sub

      integer::i
      real(kind=dp)::rho_array(0:9),z_array(0:9),d_sub,TauL_vv_prod_coef
      complex(kind=dp)::h02(0:9)

      if (layer_s>layer_o) then
         k_prop=k0*sqrt(eps_t(layer_s))
         
         d_sub=0.d0; TauL_vv_prod_coef=1.d0
         do i=layer_o+1,layer_s-1
            if (h_of_layer(i)/=-1.d0) then
               TauL_vv_prod_coef=TauL_vv_prod_coef*(1.d0+GammaL_mn(i))
               d_sub=d_sub+h_of_layer(i)
            end if
         end do

         ! zlow_of_layer(layer_s) and zlow_of_layer(layer_o+1) cannot be infinity
         z_array(0)=dabs(zlow_of_layer(layer_s)-src(2))-obs(2)+zlow_of_layer(layer_o+1)+d_sub
         z_array(1)=2.d0*zlow_of_layer(layer_s+1)-zlow_of_layer(layer_s)-src(2)-obs(2)+zlow_of_layer(layer_o+1)+d_sub
         z_array(2)=zlow_of_layer(layer_s)-2.d0*zlow_of_layer(layer_s)+src(2)-obs(2)+zlow_of_layer(layer_o+1)+d_sub
         z_array(3)=2.d0*h_of_layer(layer_s)+zlow_of_layer(layer_s)-src(2)-obs(2)+zlow_of_layer(layer_o+1)+d_sub
         z_array(4)=2.d0*h_of_layer(layer_s)-zlow_of_layer(layer_s)+src(2)-obs(2)+zlow_of_layer(layer_o+1)+d_sub
         z_array(5)=dabs(zlow_of_layer(layer_s)-src(2))+obs(2)+zlow_of_layer(layer_o+1)-2.d0*zlow_of_layer(layer_o)+d_sub
         z_array(6)=2.d0*zlow_of_layer(layer_s+1)-zlow_of_layer(layer_s)-src(2)+obs(2)&
              +zlow_of_layer(layer_o+1)-2.d0*zlow_of_layer(layer_o)+d_sub
         z_array(7)=zlow_of_layer(layer_s)+src(2)-2.d0*zlow_of_layer(layer_s)+obs(2)&
              +zlow_of_layer(layer_o+1)-2.d0*zlow_of_layer(layer_o)+d_sub
         z_array(8)=2.d0*h_of_layer(layer_s)+zlow_of_layer(layer_s)-src(2)+obs(2)&
              +zlow_of_layer(layer_o+1)-2.d0*zlow_of_layer(layer_o)+d_sub
         z_array(9)=2.d0*h_of_layer(layer_s)-zlow_of_layer(layer_s)+src(2)+obs(2)&
              +zlow_of_layer(layer_o+1)-2.d0*zlow_of_layer(layer_o)+d_sub
         
         rho_array(0:9)=dsqrt(delta_x**2+z_array(0:9)**2)
         h02(0:9)=(1.d0-c1*2.d0/pid*(dlog(0.5d0*k_prop*rho_array(0:9))+euler))

         if (zlow_of_layer(layer_o)==-1.d0) then
            h02(5:9)=cmplx(0.d0,0.d0,dp)
         end if
         
         if (zlow_of_layer(layer_s+1)==-1.d0) then
            h02(1)=cmplx(0.d0,0.d0,dp);h02(3:4)=cmplx(0.d0,0.d0,dp);&
                 h02(6)=cmplx(0.d0,0.d0,dp);h02(8:9)=cmplx(0.d0,0.d0,dp)
         end if
         
         Gf_sub=h02(0)&
              +GammaR_mn(layer_s)*h02(1)&
              +GammaL_mn(layer_s)*h02(2)&
              +GammaR_mn(layer_s)*GammaL_mn(layer_s)*h02(3)&
              +GammaR_mn(layer_s)*GammaL_mn(layer_s)*h02(4)& ! first 5 terms
              
              +GammaL_mn(layer_o)*h02(5)&
              +GammaR_mn(layer_s)*GammaL_mn(layer_o)*h02(6)&
              +GammaL_mn(layer_s)*GammaL_mn(layer_o)*h02(7)&
              +GammaR_mn(layer_s)*GammaL_mn(layer_s)*GammaL_mn(layer_o)*h02(8)&
              +GammaR_mn(layer_s)*GammaL_mn(layer_s)*GammaL_mn(layer_o)*h02(9)
         
         Gf_sub=0.5d0*Z_0(layer_s)*Gf_sub*TauL_vv_prod_coef/(2.d0*c1)
      else
         print*,'This subroutine should not be called'
         stop
      end if

      return
    end subroutine find_subtraction

    subroutine fill_TL_Green(section_s,section_o,k_rho)
      ! Calculate the transmission line Green's functions(TLGF)
      ! inputs: rs,ro,section_s,section_o,k_rho
      ! inputs: nlayers,eps_t_eff(nlayers),eps_z(nlayers),mu_t(nlayers),mu_z(nlayers),sigma(nlayers),h_of_layer(nlayers),zlow_of_layer(nlayers)
      ! output: TLGF
      
      implicit none  

      integer,intent(in)::section_s,section_o
      ! section_s: section number of the source point
      ! section_o: section number of observation point
      ! gf_num: specify the component of GF need to be calculate
      complex(kind=dp),intent(in)::k_rho
      ! k_rho: transverse wave number
      complex(kind=dp)::kz_s,kz_o
      complex(kind=dp)::D_n                                                                                    
      ! D_n: refer to (20) in the document                                                                                  
      complex(kind=dp),dimension(4)::R_n
      ! refer to (20) in the document
      complex(kind=dp),dimension(4)::ejkr_ns
      ! r_ns is r_n in (19) in the document
      ! r_ns={2*z_(n+1)-(z_s+z_o),(z_s+z_o)-2*z_n,2*d_n+(z_o-z_s),2*d_n-(z_o-z_s)}
      complex(kind=dp),dimension(nlayers)::GammaL,GammaR                                                                       !    |                                 |
      ! GammaL/right: voltage reflection coefficient seen at the terminals of section looking into the left or right           !    |GammaL                     GammaR| 
      ! GammaL(n,:)=[GammaL_h,GammaL_e], refer to Fig. 5 in the document or (50) in the document                               !    |<---                         --->|
      ! GammaR(n,:)=[GammaR_h,GammaR_e], refer to Fig. 5 in the document or (51) in the document                               !    |   |                         |   |
      complex(kind=dp),dimension(nlayers)::TauL_vv,TauR_vv
      ! TauL_vv: In one section, the voltage generated by the current source located at the right terminal of the section, refer to (62) in the document
      ! TauR_vv: In one section, the voltage generated by the current source located at the left terminal of the section, refer to (62) in the document
      ! TauL_vv(n,:)=[TauL_vv_h,TauL_vv_e]
      ! TauR_vv(n,:)=[TauR_vv_h,TauR_vv_e] 
      complex(kind=dp)::TauL_vv_prod,TauR_vv_prod 
      ! TauL_vv_prod/right: product of the TauL_vv or TauR_vv of all the sections between the source section and observation section, refer to (65) and (67) in the document
      integer::i,j,k
      ! i,j,k are loop variable
      
      ! Calculate the wave number, propogation constant and characteristic impedance of each section
      ! PEC and PMC should be handled individually
      do i=1,nlayers
         if (eps_t(i)/=-1.d0) then ! i-th layer is dielectric
            kz_wave(i)=cdsqrt(k_prop2(i)-k_rho**2)
            ! Re(kz_wave)>=0.d0: wave is propogating along +z direction
            ! Im(kz_wave)<=0.d0: wave will decay to zero when it propogates to infinity 
            if (real(kz_wave(i),8)<0.d0 .or. dimag(kz_wave(i))>0.d0) then
               kz_wave(i)=-kz_wave(i)
            end if
         end if
      end do

      GammaL(1)=cmplx(0.d0,0.d0,dp)
      GammaR(nlayers)=cmplx(0.d0,0.d0,dp)
      TauL_vv(1)=cmplx(1.d0,0.d0,dp)
      TauR_vv(nlayers)=cmplx(1.d0,0.d0,dp)
      
    do i=2,nlayers
       if (h_of_layer(i-1)==-1.d0) then
          GammaL(i)=GammaL_mn(i)
       else
          GammaL(i)=(GammaL_mn(i)+GammaL(i-1)*cdexp(-2.d0*c1*kz_wave(i-1)*h_of_layer(i-1)))&
               /(1.d0+GammaL_mn(i)*GammaL(i-1)*cdexp(-2.d0*c1*kz_wave(i-1)*h_of_layer(i-1)))
       end if

       if (h_of_layer(i)==-1.d0) then
          TauL_vv(i)=0.d0
       else
          TauL_vv(i)=(1.d0+GammaL(i))*cdexp(-c1*kz_wave(i)*h_of_layer(i))/(1.d0+GammaL(i)*cdexp(-2.d0*c1*kz_wave(i)*h_of_layer(i))) 
       end if


       if (h_of_layer(nlayers-i+2)==-1.d0) then
          GammaR(nlayers-i+1)=GammaR_mn(nlayers-i+1)
       else
          GammaR(nlayers-i+1)=(GammaR_mn(nlayers-i+1)+GammaR(nlayers-i+2)&
               *cdexp(-2.d0*c1*kz_wave(nlayers-i+2)*h_of_layer(nlayers-i+2)))&
               /(1.d0+GammaR_mn(nlayers-i+1)*GammaR(nlayers-i+2)*cdexp(-2.d0&
               *c1*kz_wave(nlayers-i+2)*h_of_layer(nlayers-i+2)))
       end if

       if (h_of_layer(nlayers-i+1)==-1.d0) then
          TauR_vv(nlayers-i+1)=0.d0
       else
          TauR_vv(nlayers-i+1)=(1.d0+GammaR(nlayers-i+1))&
               *cdexp(-c1*kz_wave(nlayers-i+1)*h_of_layer(nlayers-i+1))& 
               /(1.d0+GammaR(nlayers-i+1)*cdexp(-2.d0*c1*kz_wave(nlayers-i+1)*h_of_layer(nlayers-i+1))) 
       end if
    end do

    if (eps_t(nlayers)==-1.d0) then ! PEC
       TauL_vv(nlayers)=1.d0
    else
       TauL_vv(nlayers)=0.d0
    end if
    
    if (eps_t(1)==-1.d0) then ! PEC
       TauR_vv(1)=1.d0
    else
       TauR_vv(1)=0.d0
    end if
    
    kz_s=kz_wave(section_s)
    kz_o=kz_wave(section_o)
    ! refer to (56) in the document
    if (h_of_layer(section_s)==-1.d0) then
       D_n=1.d0
    else
       D_n=1.d0-GammaL(section_s)*GammaR(section_s)*cdexp(-2.d0*c1*kz_s*h_of_layer(section_s))
    end if
    
    R_n(1)=GammaR(section_s)
    R_n(2)=GammaL(section_s)
    R_n(3)=GammaL(section_s)*GammaR(section_s)
    R_n(4)=GammaL(section_s)*GammaR(section_s)   
    
    TLGF=cmplx(0.d0,0.d0,dp)
    TLGF_sub=cmplx(0.d0,0.d0,dp)
    TLGF_t=cmplx(0.d0,0.d0,dp)
    TLGF_sub_t=cmplx(0.d0,0.d0,dp)
    TLGF_h=cmplx(0.d0,0.d0,dp)
    TLGF_sub_h=cmplx(0.d0,0.d0,dp)

    ! refer to (57) in the document
    ejkr_ns(1)=ejkh(1)**(-c1*kz_s)
    ejkr_ns(2)=ejkh(2)**(-c1*kz_s)
    ejkr_ns(3)=ejkh(3)**(-c1*kz_s)
    ejkr_ns(4)=ejkh(4)**(-c1*kz_s)
    
    if (section_s==1) then
       if (h_of_layer(section_s)==-1.d0) then
          ejkr_ns(2)=cmplx(0.d0,0.d0,dp)
          ejkr_ns(3)=cmplx(0.d0,0.d0,dp)
          ejkr_ns(4)=cmplx(0.d0,0.d0,dp)
       end if
    elseif(section_s==nlayers) then
       if (h_of_layer(section_s)==-1.d0) then
          ejkr_ns(1)=cmplx(0.d0,0.d0,dp)
          ejkr_ns(3)=cmplx(0.d0,0.d0,dp)
          ejkr_ns(4)=cmplx(0.d0,0.d0,dp)
       end if
    end if

    d_sub=0
    if (section_s==section_o) then ! source and observation points are in the same layer 
       ! Calculate the second term in the bracket of (55) and (58)-(60) in the document
       do k=1,2 ! Hankel
          TLGF_h=TLGF_h+C_s(k)*R_n(k)*ejkr_ns(k)
          TLGF_sub_h=TLGF_sub_h+C_s(k)*R_n_sub(k)*ejkr_ns(k)
       end do

       do k=3,4 ! Toeplitz
          TLGF_t=TLGF_t+C_s(k)*R_n(k)*ejkr_ns(k)
          TLGF_sub_t=TLGF_sub_t+C_s(k)*R_n_sub(k)*ejkr_ns(k)
       end do

       TLGF_h=TLGF_h/D_n
       TLGF_t=TLGF_t/D_n
       
       TLGF_h=0.5d0*Z_0(section_s)*TLGF_h/(c1*kz_s)
       TLGF_t=0.5d0*Z_0(section_s)*TLGF_t/(c1*kz_s) ! without direct term
!       TLGF_t=0.5d0*Z_0(section_s)*(ejkh(5)**(-c1*kz_s)+TLGF_t)/(c1*kz_s) 
       ! with direct term for testing singularity subtraction
       TLGF=TLGF_t+TLGF_h

       ! singularity subtraction
       TLGF_sub_t=0.5d0*Z_0(section_s)*TLGF_sub_t/(c1*kz_s)
       TLGF_sub_h=0.5d0*Z_0(section_s)*TLGF_sub_h/(c1*kz_s)
       TLGF_sub=TLGF_sub_t+TLGF_sub_h

    else if (section_s>section_o) then ! source section is above observation section        
       TauL_vv_prod=cmplx(1.d0,0.d0,dp)
       TauL_vv_prod_sub=cmplx(1.d0,0.d0,dp)
       TauL_vv_prod_coef=cmplx(1.d0,0.d0,dp)
       if (section_o+1<=section_s-1) then     
          do i=section_o+1,section_s-1
             TauL_vv_prod=TauL_vv_prod*TauL_vv(i)
             if (h_of_layer(i)/=-1.d0) then
                TauL_vv_prod_sub=TauL_vv_prod_sub*(1.d0+GammaL_mn(i))*cdexp(-c1*kz_s*h_of_layer(i))
                TauL_vv_prod_coef=TauL_vv_prod_coef*(1.d0+GammaL_mn(i))
                d_sub=d_sub+h_of_layer(i)
             end if
          end do
       end if

       do k=1,4
          TLGF=TLGF+C_s(k)*R_n(k)*ejkr_ns(k)
          TLGF_sub=TLGF_sub+C_s(k)*R_n_sub(k)*ejkr_ns(k)
       end do

       TLGF=TLGF/D_n
       TLGF=0.5d0*Z_0(section_s)*(ejkh(5)**(-c1*kz_s)+TLGF)
       TLGF_sub=0.5d0*Z_0(section_s)*(ejkh(5)**(-c1*kz_s)+TLGF_sub)

       if (h_of_layer(section_o)==-1.d0) then
          TLGF=TLGF*TauL_vv_prod&
               *cdexp(-c1*kz_o*(zlow_of_layer(section_o+1)-ro(2))) 
          TLGF_sub=TLGF_sub*TauL_vv_prod_sub&
               *cdexp(-c1*kz_s*(zlow_of_layer(section_o+1)-ro(2))) 
       else
          TLGF=TLGF*TauL_vv_prod/(1.d0+GammaL(section_o)*cdexp(-2.d0*c1*kz_o*h_of_layer(section_o)))&
               *(1.d0+GammaL(section_o)*cdexp(-2.d0*c1*kz_o*(ro(2)-zlow_of_layer(section_o))))&
               *cdexp(-c1*kz_o*(zlow_of_layer(section_o+1)-ro(2))) 
          TLGF_sub=TLGF_sub*TauL_vv_prod_sub&
               *(1.d0+GammaL_mn(section_o)*cdexp(-2.d0*c1*kz_s*(ro(2)-zlow_of_layer(section_o))))*&
               cdexp(-c1*kz_s*(zlow_of_layer(section_o+1)-ro(2))) 
       end if
       TLGF=TLGF/(c1*kz_s)
       TLGF_sub=TLGF_sub/(c1*kz_s)
       ! Singularity subtraction       
    else if (section_s<section_o) then ! source section is below observation section 
       TauR_vv_prod=cmplx(1.d0,0.d0,dp)
       TauR_vv_prod_sub=cmplx(1.d0,0.d0,dp)       
       TauR_vv_prod_coef=cmplx(1.d0,0.d0,dp)       
       if (section_s+1<=section_o-1) then     
          do i=section_s+1,section_o-1
             TauR_vv_prod=TauR_vv_prod*TauR_vv(i)
             if (h_of_layer(i)/=-1.d0) then
                TauR_vv_prod_sub=TauR_vv_prod_sub*(1.d0+GammaR_mn(i))*cdexp(-c1*kz_s*h_of_layer(i))
                TauR_vv_prod_coef=TauR_vv_prod_coef*(1.d0+GammaR_mn(i))
                d_sub=d_sub+h_of_layer(i)
             end if
          end do
       end if

       do k=1,4
          TLGF=TLGF+C_s(k)*R_n(k)*ejkr_ns(k)
          TLGF_sub=TLGF_sub+C_s(k)*R_n_sub(k)*ejkr_ns(k)
       end do

       TLGF=TLGF/D_n
       TLGF=0.5d0*Z_0(section_s)*(ejkh(5)**(-c1*kz_s)+TLGF)
       TLGF_sub=0.5d0*Z_0(section_s)*(ejkh(5)**(-c1*kz_s)+TLGF_sub)

       if (h_of_layer(section_o)==-1.d0) then
          TLGF=TLGF*TauR_vv_prod&
               *cdexp(-c1*kz_o*(ro(2)-zlow_of_layer(section_o))) 
          TLGF_sub=TLGF_sub*TauR_vv_prod_sub&
               *cdexp(-c1*kz_s*(ro(2)-zlow_of_layer(section_o))) 
       else
          TLGF=TLGF*TauR_vv_prod/(1.d0+GammaR(section_o)*cdexp(-2.d0*c1*kz_o*h_of_layer(section_o)))&
               *(1.d0+GammaR(section_o)*cdexp(-2.d0*c1*kz_o*(zlow_of_layer(section_o+1)-ro(2))))&
               *cdexp(-c1*kz_o*(ro(2)-zlow_of_layer(section_o))) 
          TLGF_sub=TLGF_sub*TauR_vv_prod_sub&
               *(1.d0+GammaR_mn(section_o)*cdexp(-2.d0*c1*kz_s*(zlow_of_layer(section_o+1)-ro(2))))&
               *cdexp(-c1*kz_s*(ro(2)-zlow_of_layer(section_o))) 
       end if
       TLGF=TLGF/(c1*kz_s)
       TLGF_sub=TLGF_sub/(c1*kz_s)
    end if
    return
    end subroutine fill_TL_Green
    
    recursive subroutine quadrature_1D(f,a,b,qrule,eps,limit,num_sta,num,num_t,num_h,ier)
      implicit none
      integer,intent(in)::qrule,limit
      real(kind=dp),intent(in)::a,b,eps
      complex(kind=dp),intent(in)::num_sta
      complex(kind=dp),intent(out)::num,num_t,num_h
      integer,intent(out)::ier
      real(kind=dp)::x,diff
      complex(kind=dp)::num_tmp1,num_tmp2,num_tmp3,num_tmp1_t,num_tmp2_t,num_tmp3_t,&
           num_tmp1_h,num_tmp2_h,num_tmp3_h
      complex(kind=dp),external::f
      
      num=cmplx(0.d0,0.d0,dp)
      num_t=cmplx(0.d0,0.d0,dp)
      num_h=cmplx(0.d0,0.d0,dp)

      diff=b-a
      call quadrature(f,a,a+0.5d0*diff,qrule,num_tmp2,num_tmp2_t,num_tmp2_h)
      call quadrature(f,a+0.5d0*diff,b,qrule,num_tmp3,num_tmp3_t,num_tmp3_h)

      if (cdabs(num_sta-num_tmp2-num_tmp3)>eps*cdabs(num_tmp2+num_tmp3)) then
!      if (cdabs(num_sta-num_tmp2-num_tmp3)>eps*cdabs(num_tmp2+num_tmp3) .and. cdabs(num_sta-num_tmp2-num_tmp3)>1.d-16) then
         call quadrature_1D(f,a,a+0.5d0*diff,qrule,eps,limit,num_tmp2,num_tmp1,num_tmp1_t,num_tmp1_h,ier)
         num=num+num_tmp1
         num_t=num_t+num_tmp1_t
         num_h=num_h+num_tmp1_h
         call quadrature_1D(f,a+0.5d0*diff,b,qrule,eps,limit,num_tmp3,num_tmp2,num_tmp2_t,num_tmp2_h,ier)
         num=num+num_tmp2
         num_t=num_t+num_tmp2_t
         num_h=num_h+num_tmp2_h
      else
         num=num+num_tmp2+num_tmp3
         num_t=num_t+num_tmp2_t+num_tmp3_t
         num_h=num_h+num_tmp2_h+num_tmp3_h
      end if
      
      ier=0
      return
    end subroutine quadrature_1D
    
    subroutine quadrature(f,a,b,qrule,num,num_t,num_h)
      use global_com,only:dp
      implicit none
      integer,intent(in)::qrule
      real(kind=dp),intent(in)::a,b
      complex(kind=dp),intent(out)::num,num_t,num_h
      real(kind=dp)::x,wght_tmp
      complex(kind=dp)::fval,fval_t,fval_h
      integer::i
      complex(kind=dp),external::f
      
      num=0.d0
      num_t=0.d0
      num_h=0.d0
      
      do i=1,qrule
         x=(b-a)*qp(i)+a
         wght_tmp=wght(i)
         fval=f(x)
         ! only for m==n
         num_t=num_t+ft*wght_tmp
         num_h=num_h+fh*wght_tmp
         ! for m/=n and m==n
         num=num+fval*wght_tmp
!         print*,'num',i,fval
      end do
      ! only for m==n
      num_t=num_t*(b-a)
      num_h=num_h*(b-a)
!      num=num_t+num_h
      ! for m/=n and m==n
      num=num*(b-a)
     ! stop
      return
    end subroutine quadrature
    
    subroutine oneD_quadrature(qrule,qp,wght)  
      use global_com,only:dp
      implicit none
      integer,intent(in)::qrule
      real(kind=dp),intent(out),dimension(1:qrule)::qp,wght
      
      select case (qrule)
      case(1)
         qp(1)=0.5d0
         wght(1)=1.0d0
      case(2)
         qp(1:2)=(/0.211324865405187d0,0.788675134594813d0/)
         wght(1:2)=(/0.5d0,0.5d0/)
      case(3)
         qp(1:3)=(/0.5d0,0.112701665379258d0,0.887298334620742d0/)
         wght(1:3) = (/0.444444444444444d0,0.277777777777776d0,&
              0.277777777777776d0/)
      case(4)
         qp(1:4)= (/0.330009478207572d0,6.943184420297372d-2,&
              0.669990521792428d0,0.930568155797026d0/)
         wght(1:4)=(/0.326072577431273d0,0.173927422568727d0,&
              0.326072577431273d0,0.173927422568727d0/)
      case(5)
         qp(1:5)=(/0.5d0,4.691007703066802d-2,0.230765344947158d0,&
              0.769234655052841d0,0.953089922969332d0/)
         wght(1:5)=(/0.284444444444444d0,0.118463442528095d0,&
              0.239314335249683d0,0.239314335249683d0,&
              0.118463442528095d0/)
      case(6)
         qp(1:6)=(/0.619309593041598d0,3.376524289842397d-2,&
              0.169395306766868d0,0.380690406958402d0,&
              0.830604693233132d0,0.966234757101576d0/)
         wght(1:6)=(/0.233956967286346d0,8.566224618958525d-2,&
              0.180380786524069d0,0.233956967286346d0,&
              0.180380786524069d0,8.566224618958525d-2/)
      case(7)
         qp(1:7)=(/0.5d0,2.544604382862076d-2,0.129234407200303d0,&
              0.297077424311301d0,0.702922575688699d0,&
              0.870765592799697d0,0.974553956171379d0/)
         wght(1:7)=(/0.208979591836735d0,6.474248308443483d-2,&
              0.139852695744638d0,0.190915025252560d0,&
              0.190915025252560d0,0.139852695744638d0,&
              6.474248308443483d-2/)
      case(8)
         qp(1:8)=(/0.408282678752175d0,1.985507175123186d-2,&
              0.101666761293187d0,0.237233795041836d0,&
              0.591717321247825d0,0.762766204958164d0,&
              0.898333238706813d0,0.980144928248768d0/)
         wght(1:8)=(/0.181341891689181d0,5.061426814518839d-2,&
              0.111190517226687d0,0.156853322938944d0,&
              0.181341891689181d0,0.156853322938944d0,&
              0.111190517226687d0,5.061426814518839d-2/)
      case(9)
         qp(1:9)=(/0.5d0,1.591988024618696d-2,8.198444633668206d-2,&
              0.193314283649705d0,0.337873288298096d0,&
              0.662126711701904d0,0.806685716350295d0,&
              0.918015553663318d0,0.984080119753813d0/)
         wght(1:9)=(/0.165119677500630d0,4.063719418078731d-2,&
              9.032408034742879d-2,0.130305348201468d0,&
              0.156173538520001d0,0.156173538520001d0,&
              0.130305348201468d0,9.032408034742879d-2,&
              4.063719418078731d-2/)
      case(10)
         qp(1:10)=(/0.574437169490816d0,1.304673574141413d-2,&
              6.746831665550773d-2,0.160295215850488d0,&
              0.283302302935376d0,0.425562830509184d0,&
              0.716697697064624d0,0.839704784149512d0,&
              0.932531683344492d0,0.986953264258586d0/)
         wght(1:10)=(/0.147762112357376d0,3.333567215434186d-2,&
              7.472567457529027d-2,0.109543181257991d0,&
              0.134633359654998d0,0.147762112357376d0,&
              0.134633359654998d0,0.109543181257991d0,&
              7.472567457529027d-2,3.333567215434186d-2/)
      end select
      return

    end subroutine oneD_quadrature

    subroutine Gf_interpolation_2d(src,obs,Gf_intpl_t,Gf_intpl_h)
      implicit none
      
      real(kind=dp),intent(in)::src(2),obs(2)
      complex(kind=dp),intent(out)::Gf_intpl_t,Gf_intpl_h
      
      real(kind=dp)::l_coef_rho(6),l_coef_zt(6),l_coef_zh(6)
      real(kind=dp)::rho_tmp,z_t_tmp,z_h_tmp
      integer::i,j,k,rho_index(6),zt_index(6),zh_index(6)
      integer::index1,index2,index3 ! rho,z-z',z+z'
      integer::ns,no
      complex(kind=dp)::Gf_intpl_tmp_t(6),Gf_intpl_tmp_h(6)
      complex(kind=dp),allocatable::Gf_table_t(:,:),Gf_table_h(:,:)
      integer::h_sta

      ns=map_layer(layer_s)
      no=map_layer(layer_o)
      allocate(Gf_table_t(0:num_rho(ns,no),0:num_z(ns)))
      allocate(Gf_table_h(0:num_rho(ns,no),0:2*num_z(ns)+1))
      Gf_table_t(:,:)=gf_table_same(ns)%Gf_grid_array_t(:,:)
      Gf_table_h(:,:)=gf_table_same(ns)%Gf_grid_array_h(:,:)

      Gf_intpl_t=cmplx(0.d0,0.d0,dp)
      Gf_intpl_h=cmplx(0.d0,0.d0,dp)

      Gf_intpl_tmp_t(:)=cmplx(0.d0,0.d0,dp)
      Gf_intpl_tmp_h(:)=cmplx(0.d0,0.d0,dp)

      rho_index(:)=0
      zt_index(:)=0
      zh_index(:)=0

      rho_tmp=dabs(obs(1)-src(1))
      z_t_tmp=abs(obs(2)-src(2))
      z_h_tmp=obs(2)+src(2)
         
      index1=dint(rho_tmp/drho(ns,no))

      ! To find the approximated Green's functions becomes 2 2D array interpolation problem
      ! Toeplitz part
      if (index1<2) then
         rho_index(1:6)=(/0,1,2,3,4,5/)
      else if (index1>=2) then
         ! central
         rho_index(1:6)=(/index1-2,index1-1,index1,index1+1,index1+2,index1+3/)
      end if
  
      if (num_z(ns)/=1) then
         index2=dint(z_t_tmp/dz(ns))
         if (z_h_tmp>z_min(ns)+z_max(ns)) then
            index3=dint((z_h_tmp-z_max(ns)-z_min(ns))/dz(ns)) ! h2 array
            z_h_tmp=z_h_tmp-z_max(ns)-z_min(ns) ! set z_h_tmp starting at 0
            h_sta=num_z(ns)+1
         else if (z_h_tmp<z_min(ns)+z_max(ns)) then
            index3=dint((z_h_tmp-2*z_min(ns))/dz(ns)) ! h1 array
            z_h_tmp=z_h_tmp-2*z_min(ns) ! set z_h_tmp starting at 0
            h_sta=0
         else
            if (src(2)==z_min(ns)) then
               index3=num_z(ns) ! h1 array
               z_h_tmp=z_h_tmp-2*z_min(ns)
               h_sta=0
            else
               index3=0 ! h2 array
               z_h_tmp=0.d0
               h_sta=num_z(ns)+1
            end if
         end if

         if (index2<2) then
            ! forward
            zt_index(1:6)=(/0,1,2,3,4,5/)
         else if (index2>=num_z(ns)-2) then
            ! backward
            zt_index(1:6)=(/num_z(ns)-5,num_z(ns)-4,num_z(ns)-3,&
                 num_z(ns)-2,num_z(ns)-1,num_z(ns)/)
         else
            ! central
            zt_index(1:6)=(/index2-2,index2-1,index2,index2+1,index2+2,index2+3/)
         end if
         
         if (index3<2) then
            ! forward
            zh_index(1:6)=(/0,1,2,3,4,5/)
         else if (index3>=num_z(ns)-2) then
            ! backward
            zh_index(1:6)=(/num_z(ns)-5,num_z(ns)-4,num_z(ns)-3,&
                 num_z(ns)-2,num_z(ns)-1,num_z(ns)/)
         else
            ! central
            zh_index(1:6)=(/index3-2,index3-1,index3,index3+1,index3+2,index3+3/)
         end if

         l_coef_rho(:)=1.d0
         do i=1,6
            do j=1,6
               if (j/=i) then
                 ! l_coef_rho(i)=l_coef_rho(i)*(rho_tmp-(rho_index(j)*drho))/(rho_index(i)-rho_index(j))/drho
                  l_coef_rho(i)=l_coef_rho(i)*(rho_tmp-(rho_index(j)*drho(ns,no)))
               end if
            end do
            l_coef_rho(i)=l_coef_rho(i)*drho_ij(i,ns,no)
         end do

         l_coef_zt(:)=1.d0
         l_coef_zh(:)=1.d0

         do k=1,6
            do j=1,6
               Gf_intpl_tmp_t(k)=Gf_intpl_tmp_t(k)+Gf_table_t(rho_index(j),zt_index(k))*l_coef_rho(j)
               Gf_intpl_tmp_h(k)=Gf_intpl_tmp_h(k)+Gf_table_h(rho_index(j),h_sta+zh_index(k))*l_coef_rho(j)
            end do

            do j=1,6
               if (j/=k) then
                  l_coef_zt(k)=l_coef_zt(k)*(z_t_tmp-(zt_index(j)*dz(ns)))
                  l_coef_zh(k)=l_coef_zh(k)*(z_h_tmp-(zh_index(j)*dz(ns)))
               end if
            end do
            
            l_coef_zt(k)=l_coef_zt(k)*dz_ij(k,ns)
            l_coef_zh(k)=l_coef_zh(k)*dz_ij(k,ns)

            Gf_intpl_t=Gf_intpl_t+Gf_intpl_tmp_t(k)*l_coef_zt(k)
            Gf_intpl_h=Gf_intpl_h+Gf_intpl_tmp_h(k)*l_coef_zh(k)
         end do
      else
         index2=0;index3=0
         l_coef_rho(:)=1.d0
         do i=1,6
            do j=1,6
               if (j/=i) then
                  l_coef_rho(i)=l_coef_rho(i)*(rho_tmp-(rho_index(j)*drho(ns,no)))&
                       /(rho_index(i)-rho_index(j))/drho(ns,no)
               end if
            end do
            Gf_intpl_t=Gf_intpl_t+Gf_table_t(rho_index(i),index2)*l_coef_rho(i)
            Gf_intpl_h=Gf_intpl_h+Gf_table_h(rho_index(i),index3)*l_coef_rho(i)
         end do
      end if
      return
    end subroutine Gf_interpolation_2d

    subroutine Gf_interpolation_3d(src,obs,Gf_intpl)
      implicit none
      
      real(kind=dp),intent(in)::src(2),obs(2)
      complex(kind=dp),intent(out)::Gf_intpl
      
      real(kind=dp)::l_coef_rho(6),l_coef_zs(6),l_coef_zo(6)
      real(kind=dp)::rho_tmp,zs_tmp,zo_tmp
      integer::i,j,k,rho_index(6),zs_index(6),zo_index(6)
      integer::index1,index2,index3 ! rho,zs,zo
      integer::ns,no
      complex(kind=dp)::Gf_intpl_tmp(6),Gf_intpl_tmp_2d(6,6)
      complex(kind=dp),allocatable::gf_table(:,:,:)
      integer::h_sta

      ns=map_layer(layer_s)
      no=map_layer(layer_o)
      allocate(Gf_table(0:num_rho(no,ns),0:num_z(no),0:num_z(ns)))
      Gf_table(0:num_rho(no,ns),0:num_z(no),0:num_z(ns))=&
           gf_table_diff(no,ns)%Gf_grid_array(0:num_rho(no,ns),0:num_z(no),0:num_z(ns))

      Gf_intpl=cmplx(0.d0,0.d0,dp)
      Gf_intpl_tmp(:)=cmplx(0.d0,0.d0,dp)
      Gf_intpl_tmp_2d(:,:)=cmplx(0.d0,0.d0,dp)

      rho_index(:)=0
      zs_index(:)=0
      zo_index(:)=0

      rho_tmp=dabs(obs(1)-src(1))
      zs_tmp=src(2)-z_min(ns)
      zo_tmp=obs(2)-z_min(no)
         
      index1=dint(rho_tmp/drho(no,ns))

      ! To find the approximated Green's functions becomes 2 2D array interpolation problem
      ! Toeplitz part
      if (index1<2) then
         rho_index(1:6)=(/0,1,2,3,4,5/)
      else if (index1>=2) then
         ! central
         rho_index(1:6)=(/index1-2,index1-1,index1,index1+1,index1+2,index1+3/)
      end if
  

      index2=dint(zs_tmp/dz(ns))
      if (index2<2) then
         ! forward
         zs_index(1:6)=(/0,1,2,3,4,5/)
      else if (index2>=num_z(ns)-2) then
         ! backward
         zs_index(1:6)=(/num_z(ns)-5,num_z(ns)-4,num_z(ns)-3,&
              num_z(ns)-2,num_z(ns)-1,num_z(ns)/)
      else
         ! central
         zs_index(1:6)=(/index2-2,index2-1,index2,index2+1,index2+2,index2+3/)
      end if
      
      index3=dint(zo_tmp/dz(no))
      if (index3<2) then
         ! forward
         zo_index(1:6)=(/0,1,2,3,4,5/)
      else if (index3>=num_z(no)-2) then
         ! backward
         zo_index(1:6)=(/num_z(no)-5,num_z(no)-4,num_z(no)-3,&
              num_z(no)-2,num_z(no)-1,num_z(no)/)
      else
         ! central
         zo_index(1:6)=(/index3-2,index3-1,index3,index3+1,index3+2,index3+3/)
      end if
      
      l_coef_rho(:)=1.d0
      l_coef_zs(:)=1.d0
      l_coef_zo(:)=1.d0
      do i=1,6
         do j=1,6
            if (j/=i) then
               ! l_coef_rho(i)=l_coef_rho(i)*(rho_tmp-(rho_index(j)*drho))/(rho_index(i)-rho_index(j))/drho
               l_coef_rho(i)=l_coef_rho(i)*(rho_tmp-(rho_index(j)*drho(no,ns)))
               l_coef_zs(i)=l_coef_zs(i)*(zs_tmp-(zs_index(j)*dz(ns)))
               l_coef_zo(i)=l_coef_zo(i)*(zo_tmp-(zo_index(j)*dz(no)))
            end if
         end do
         l_coef_rho(i)=l_coef_rho(i)*drho_ij(i,no,ns)
         l_coef_zs(i)=l_coef_zs(i)*dz_ij(i,ns)
         l_coef_zo(i)=l_coef_zo(i)*dz_ij(i,no)
      end do
      

      ! 3D to 2D
      do k=1,6
         do j=1,6
            do i=1,6
               Gf_intpl_tmp_2d(j,k)=Gf_intpl_tmp_2d(j,k)+Gf_table(rho_index(i),zo_index(j),zs_index(k))*l_coef_rho(i)
            end do
         end do
      end do

      ! 2D to 1D
      do j=1,6
         do i=1,6
            Gf_intpl_tmp(j)=Gf_intpl_tmp(j)+Gf_intpl_tmp_2d(i,j)*l_coef_zo(i)
         end do
      end do

      ! 1D to point
      do i=1,6
         Gf_intpl=Gf_intpl+Gf_intpl_tmp(i)*l_coef_zs(i)
      end do

      return
    end subroutine Gf_interpolation_3d
  end module layers
  
  function f1_2d(x)
    use layers
    implicit none     
    
    real(kind=dp),intent(in)::x 
    complex(kind=dp)::k_rho,f1_2d
    
    k_rho=c1*x
    call fill_TL_Green(layer_s,layer_o,k_rho)

    kernel_t=TLGF_t
    kernel_h=TLGF_h
    kernel=TLGF

!    ft=kernel_t*cdcos(k_rho*delta_x)*c1
!    fh=kernel_h*cdcos(k_rho*delta_x)*c1
!    f1_2d=kernel*cdcos(k_rho*delta_x)*c1

    ft=(kernel_t-TLGF_sub_t)*cdcos(k_rho*delta_x)*c1
    fh=(kernel_h-TLGF_sub_h)*cdcos(k_rho*delta_x)*c1
    f1_2d=(kernel-TLGF_sub)*cdcos(k_rho*delta_x)*c1
    return
  end function f1_2d

  function f2_2d(x)
    use layers
    implicit none     
    
    real(kind=dp),intent(in)::x 
    complex(kind=dp)::k_rho,f2_2d
    
    k_rho=x+c1*h1
    call fill_TL_Green(layer_s,layer_o,k_rho)

    kernel_t=TLGF_t
    kernel_h=TLGF_h
    kernel=TLGF

!    ft=kernel_t*cdcos(k_rho*delta_x)
!    fh=kernel_h*cdcos(k_rho*delta_x)
!    f2_2d=kernel*cdcos(k_rho*delta_x)

    ft=(kernel_t-TLGF_sub_t)*cdcos(k_rho*delta_x)
    fh=(kernel_h-TLGF_sub_h)*cdcos(k_rho*delta_x)
    f2_2d=(kernel-TLGF_sub)*cdcos(k_rho*delta_x)
    return
  end function f2_2d

  function f3_2d(x)
    use layers
    implicit none     
    
    real(kind=dp),intent(in)::x 
    complex(kind=dp)::k_rho,f3_2d
    
    k_rho=h2+c1*x
    call fill_TL_Green(layer_s,layer_o,k_rho)

    kernel_t=TLGF_t
    kernel_h=TLGF_h
    kernel=TLGF

!    ft=kernel_t*cdcos(k_rho*delta_x)*c1
!    fh=kernel_h*cdcos(k_rho*delta_x)*c1
!    f3_2d=kernel*cdcos(k_rho*delta_x)*c1

    ft=(kernel_t-TLGF_sub_t)*cdcos(k_rho*delta_x)*c1
    fh=(kernel_h-TLGF_sub_h)*cdcos(k_rho*delta_x)*c1
    f3_2d=(kernel-TLGF_sub)*cdcos(k_rho*delta_x)*c1
    return
  end function f3_2d

  function f4_2d(x)
    use layers
    implicit none     
    
    real(kind=dp),intent(in)::x 
    complex(kind=dp)::k_rho,f4_2d
    
    k_rho=x
    call fill_TL_Green(layer_s,layer_o,k_rho)

    kernel_t=TLGF_t
    kernel_h=TLGF_h
    kernel=TLGF

!    ft=kernel_t*cdcos(k_rho*delta_x)
!    fh=kernel_h*cdcos(k_rho*delta_x)
!    f4_2d=kernel*cdcos(k_rho*delta_x)

    ft=(kernel_t-TLGF_sub_t)*cdcos(k_rho*delta_x)
    fh=(kernel_h-TLGF_sub_h)*cdcos(k_rho*delta_x)
    f4_2d=(kernel-TLGF_sub)*cdcos(k_rho*delta_x)
    return
  end function f4_2d
