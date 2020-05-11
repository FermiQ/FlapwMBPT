      subroutine output_resp_mbpt_static(p00,p0z,pz0,pw_00,pw_zz,vc,pw,
     &                                   key)
c     key = 0 Non-interacting    
c     key = 1 Vertex-Corrected      
	use atom_mod
	use manager_mod
	use models_mod
	use parallel_mod
	use solid_mod
	use units_mod
      use vertex_mod
      implicit none
      integer, intent(in) :: key
      real*8, intent(in) :: pw_00(nc_line),pw_zz(nc_line)
	complex*16, intent(in) :: p00(n_pbtot,n_pbtot,nc_line),
     & 	                      p0z(n_pbtot,n_pbtot,nc_line),
     &                          pz0(n_pbtot,n_pbtot,nc_line),
     &                          vc(n_pbtot,n_pbtot,nc_line),
     &                          pw(n_pbtot,nc_line)
      character*2 :: ch
	integer :: k,i_len,i,iq,n,ka
	real*8 :: fc,dpwz,qq,pi2a
	real*8, allocatable :: epw(:),epw1(:),dpw_00(:),dpw_zz(:)
	complex*16, allocatable :: epb(:,:),ddr(:,:),r0z(:,:),pv1(:,:)
	pi2a=(pi+pi)/par
	if(key==0) ch='_0'
	if(key==1) ch='_1'
	allocate(epb(n_pbtot,n_pbtot))
	allocate(ddr(n_pbtot,n_pbtot))
	allocate(r0z(n_pbtot,n_pbtot))
	allocate(pv1(n_pbtot,n_pbtot))
	allocate(epw(nc_line))
	allocate(epw1(nc_line))
	allocate(dpw_00(nc_line))
	allocate(dpw_zz(nc_line))
	epw=0.d0
	epw1=0.d0
	dpw_00=0.d0
	dpw_zz=0.d0
	do iq=1,nc_line
	  k=kline_in_npnt_c(iq)
        ka=k_a_from_c(k)
	  n=n_pbmt+nplwgw(ka)
c -------- Prepare (1-PV)^(-1) --------------------------------------
	  call zgemm('n','n',n,n,n,(1.d0,0.d0),p00(1,1,iq),n_pbtot,
     &	         vc(1,1,iq),n_pbtot,(0.d0,0.d0),pv1,n_pbtot)
        pv1=-pv1
        do i=1,n
          pv1(i,i)=(1.d0,0.d0)+pv1(i,i)
        enddo
        call invers_z(n,pv1,n_pbtot)
	  if(chi_cmp(1:1)/='0') then
c  ------- Dielectric function  ------------------------------
	    call zgemm('n','n',n,n,n,(1.d0,0.d0),vc(1,1,iq),n_pbtot,
     &	           p00(1,1,iq),n_pbtot,(0.d0,0.d0),epb,n_pbtot)
          epb=-epb
          do i=1,n
            epb(i,i)=(1.d0,0.d0)+epb(i,i)
          enddo
          call pw_from_pb_line(n,epb,epw(iq),2,pw(1,iq))
c  ------- Inverse Dielectric function  ------------------------------
          call invers_z(n,epb,n_pbtot)
          call pw_from_pb_line(n,epb,epw1(iq),2,pw(1,iq))
c --------- Charge-Charge Response  ----------------------------------
	    call zgemm('n','n',n,n,n,(1.d0,0.d0),pv1,n_pbtot,
     &	             p00(1,1,iq),n_pbtot,(0.d0,0.d0),ddr,n_pbtot)
          ddr=ddr-p00(:,:,iq)
          call pw_from_pb_line(n,ddr,dpw_00(iq),3,pw(1,iq))
          dpw_00(iq)=dpw_00(iq)+pw_00(iq)
        endif
	  if(chi_cmp(3:3)/='0') then
c --------- Spin-Spin (ZZ) Response  ----------------------------------
          dpw_zz(iq)=pw_zz(iq)
          if(nspin==2) then
	      call zgemm('n','n',n,n,n,(1.d0,0.d0),pv1,n_pbtot,
     &	               p0z(1,1,iq),n_pbtot,(0.d0,0.d0),r0z,n_pbtot)
	      call zgemm('n','n',n,n,n,(1.d0,0.d0),vc(1,1,iq),n_pbtot,
     &	               r0z,n_pbtot,(0.d0,0.d0),ddr,n_pbtot)
	      call zgemm('n','n',n,n,n,(1.d0,0.d0),pz0(1,1,iq),n_pbtot,
     &	               ddr,n_pbtot,(0.d0,0.d0),r0z,n_pbtot)
            call pw_from_pb_line(n,r0z,dpwz,3,pw(1,iq))
            dpw_zz(iq)=dpw_zz(iq)+dpwz
          endif
        endif
      enddo  !! over iq
	deallocate(epb,ddr,r0z,pv1)
	if(maswrk) then
	  fc=amega
	  i_len=len_trim(allfile)
	  if(ubi==' gw') then
	    open(3,file=allfile(1:i_len)//'_Responses_MBPT'//ch//'.gw')
	  else if(ubi==' hf') then
	    open(3,file=allfile(1:i_len)//'_Responses_MBPT'//ch//'.hf')
	  else if(ubi==' qp') then
	    open(3,file=allfile(1:i_len)//'_Responses_MBPT'//ch//'.qp')
	  else if(ubi=='psi') then
	    open(3,file=allfile(1:i_len)//'_Responses_MBPT'//ch//'.psi')
	  else if(ubi=='bsp') then
	    open(3,file=allfile(1:i_len)//'_Responses_MBPT'//ch//'.bsp')
	  endif
	  write(3,*)'Q     E    E^-1    Xi_00(1/Ry/u.c.) Xi_ZZ(1/Ry/u.c.)'
	  do iq=1,nc_line
          k=kline_in_npnt_c(iq)
          qq=sqrt(dot_product(pnt_c(:,k),pnt_c(:,k)))*pi2a
 	    write(3,'(f11.5,4(1x,e12.5))')qq,epw(iq),epw1(iq),
     &	                                dpw_00(iq)*fc,dpw_zz(iq)*fc
	  enddo
	  close(3)
	endif
	deallocate(epw,epw1,dpw_00,dpw_zz)
      end
