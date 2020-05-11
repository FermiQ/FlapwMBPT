      subroutine sigma_g_lambda_heg_0(iq,nom,ind_nu,lambda_dyn,g_om,
     &                                w3_nu_d,aw,bw)
      use atom_mod
      use heg_mod
      use manager_mod
      use models_mod
      use parallel_mod
	use solid_mod
	use units_mod
	use vertex_mod
      implicit none
      integer, intent(in) :: nom,ind_nu,iq
      real*8, intent(in) :: w3_nu_d
      complex*16, intent(in) :: lambda_dyn(nrr_div_c,nom),
     &                          g_om(nstar_c,nom,2)
	complex*16, intent(inout) :: aw(nq_sph_c,0:n_tau),
     &	                         bw(nq_sph_c,0:n_tau)
      integer :: i_omega,k,kq,kq0,in,k0
      real*8 :: fc,v(3)
      complex*16, allocatable :: tmp(:,:),tmp1(:,:)
      allocate(tmp1(nq_sph_c,nom))
      allocate(tmp(nq_sph_c,0:n_tau))
      do in=1,2
        tmp1=(0.d0,0.d0)
	  do k0=1,nq_sph_c
          k=q_sph_in_cube_c(k0)
	    v=q_mesh_heg_c(:,k)+q_irr_c(:,iq)
	    if(in==1) v=-v
	    call zone1_number(v,rb0_c,nr_div_c,kq)
	    kq0=q_cube_in_sph_c(kq)
	    if(kq0==0) cycle
	    kq0=i_kref_c(kq0)
	    if(in==1) then
	      v=-q_mesh_heg_c(:,k)
	      call zone1_number(v,rb0_c,nr_div_c,kq)
	    endif
	    do i_omega=1,nom
	      tmp1(k0,i_omega)=g_om(kq0,i_omega,3-in)
     &                      *lambda_dyn(kq,i_omega)
          enddo
        enddo  !! over k
        call from_omega_nu_to_tau_nu_aa(nq_sph_c,ind_nu,tmp1,tmp,in)
c ------------------------------------------------------------------
        fc=w3_nu_d*wgt_c(iq)/amega
	  do k=1,nq_sph_c
	    if(in==1) aw(k,:)=aw(k,:)+tmp(k,:)*fc
	    if(in==2) bw(k,:)=bw(k,:)+tmp(k,:)*fc
        enddo  !! over k
      enddo   !! over in
      deallocate(tmp1,tmp)
      end