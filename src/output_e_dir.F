      subroutine output_e_dir(p_dir,w_dir,e_dir,e1_dir)
      use atom_mod
      use manager_mod
      use models_mod
      use parallel_mod
      use solid_mod
      use units_mod
      use vertex_mod
      implicit none
      real*8, intent(in) :: p_dir(n_phead_dir,ndim3_nu),
     &                      w_dir(n_phead_dir,ndim3_nu),
     &                      e_dir(n_phead_dir,ndim3_nu),
     &                      e1_dir(n_phead_dir,ndim3_nu)
      integer :: ind_nu,i_nu
      real*8, allocatable :: f(:,:)
      if(me_k_b/=0) return
      allocate(f(0:n_nu,n_phead_dir))
      f=0.d0
      do ind_nu=1,ndim3_nu
        i_nu=me_t*ndim3_nu+ind_nu-1
        f(i_nu,:)=p_dir(:,ind_nu)
      enddo
      if(nproc_t/=1) call dgop(f,(n_nu+1)*n_phead_dir,'  +',comm_t)
      call output_boson_dir(f,'_____P')
      f=0.d0
      do ind_nu=1,ndim3_nu
        i_nu=me_t*ndim3_nu+ind_nu-1
        f(i_nu,:)=w_dir(:,ind_nu)
      enddo
      if(nproc_t/=1) call dgop(f,(n_nu+1)*n_phead_dir,'  +',comm_t)
      call output_boson_dir(f,'_____W')
      f=0.d0
      do ind_nu=1,ndim3_nu
        i_nu=me_t*ndim3_nu+ind_nu-1
        f(i_nu,:)=e_dir(:,ind_nu)
      enddo
      if(nproc_t/=1) call dgop(f,(n_nu+1)*n_phead_dir,'  +',comm_t)
      call output_boson_dir(f,'_____E')
      f=0.d0
      do ind_nu=1,ndim3_nu
        i_nu=me_t*ndim3_nu+ind_nu-1
        f(i_nu,:)=e1_dir(:,ind_nu)
      enddo
      if(nproc_t/=1) call dgop(f,(n_nu+1)*n_phead_dir,'  +',comm_t)
      call output_boson_dir(f,'____E1')
      deallocate(f)
      end
