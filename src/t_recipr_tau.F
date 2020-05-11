      subroutine t_recipr_tau(qt_pw,v2_q,iq)
	  use atom_mod
	  use manager_mod
	  use parallel_mod
	  use solid_mod
	  use units_mod
	  use vertex_mod
      implicit none
      integer, intent(in) :: iq
      complex*16, intent(in) :: v2_q(n_pbtot_red,n_pbtot_red,nqdiv_c)
      complex*16, intent(inout) :: qt_pw(n_pbtot_red,n_pbtot_red,
     &                                   nqdiv_c,2,2,0:n_tau/2)
      integer :: iq1,iq2,iq10,np,iq20,np2,i_tau,ka,ka2,it,il
      real*8 :: v(3)
      complex*16, allocatable :: tmp(:,:)
      allocate(tmp(n_pbtot_red,n_pbtot_red))
      do iq1=1,nqdiv_c
	    iq10=i_kref_c(iq1)
        ka=k_a_from_c(iq10)
	    np=n_pbmt_red+nplwgw_red(ka)
        v=pnt_c(:,iq1)-pnt_c(:,iq)
        call zone1_number(v,rb0_c,ndiv_c,iq2)  ! iq2 - index in BZ_1
        iq2=index_k1_c(iq2)  ! Now iq2 - index in BZ_0
	    iq20=i_kref_c(iq2)
        ka2=k_a_from_c(iq20)
	    np2=n_pbmt_red+nplwgw_red(ka2)
        do i_tau=0,n_tau/2
          do it=1,2
            do il=1,2
	          call zgemm('n','n',np,np2,np,(1.d0,0.d0),v2_q(1,1,iq1),
     &	                 n_pbtot_red,qt_pw(1,1,iq1,il,it,i_tau),
     &                   n_pbtot_red,(0.d0,0.d0),tmp,n_pbtot_red)
              call zgemm('n','n',np,np2,np2,(1.d0,0.d0),tmp,n_pbtot_red,
     &                   v2_q(1,1,iq2),n_pbtot_red,(0.d0,0.d0),
     &                   qt_pw(1,1,iq1,il,it,i_tau),n_pbtot_red)
            enddo   !! over il
          enddo   !! over it
        enddo  !! over i_tau
      enddo  !! over iq1
      deallocate(tmp)
      end
