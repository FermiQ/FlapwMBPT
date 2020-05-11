      subroutine self_energy_1_dyn_b(aw,sig)
      use atom_mod
      use heg_mod
      use manager_mod
      use models_mod
      use parallel_mod
	  use solid_mod
	  use units_mod
	  use vertex_mod
      implicit none
      complex*16, intent(inout) :: aw(nra_red,nra_red,nqdiv_c,0:n_tau)
	  complex*16, intent(out) :: sig(nra_red,nra_red,nqdiv_c,0:n_tau)
      integer :: i_tau,k,j,i,k1,i1,j1,ka,i0,j0,gbs_number,ka0,na,k1a
      real*8 :: v(3),ga(3)
      do k=1,nqdiv_c
        ka=k_a_from_c(k)
        ka0=i_kref(ka)
        na=nfun_red+nbask_red(ka0)
	    v=-pnt_c(:,k)
	    call zone1_number(v,rb0_c,ndiv_c,k1)
	    k1=index_k1_c(k1)     !! for -K here
        k1a=k_a_from_c(k1)
        ga=pnt_c(:,k1)+pnt_c(:,k)
	    do i_tau=0,n_tau
          do j=1,na
            if(j<=nfun_red) then
              j1=j
            else
	          j0=indgb_red(j-nfun_red,ka)  !! G in GBS-list
              v=-gbs(:,j0)-ga
              j1=gbs_number(v)
              j1=nfun_red+iplf_bk_red(j1,k1a)
            endif
            do i=1,na
              if(i<=nfun_red) then
                i1=i
              else
	          i0=indgb_red(i-nfun_red,ka)  !! G in GBS-list
                v=-gbs(:,i0)-ga
                i1=gbs_number(v)
                i1=nfun_red+iplf_bk_red(i1,k1a)
              endif
	          sig(i,j,k,i_tau)=-(aw(i,j,k,i_tau)
     &	                       +conjg(aw(i1,j1,k1,i_tau)))
            enddo
	      enddo
        enddo
      enddo   !! over k
      end
