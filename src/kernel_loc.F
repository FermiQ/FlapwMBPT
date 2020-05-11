      subroutine kernel_loc_1_1(n,n2,th0,nu,v4,wexa,wasy,wspl_geom,
     &                          wspl_asy)
      use atom_mod
      use etot_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use units_mod
      use vertex_mod
      implicit none
      integer, intent(in) :: n,n2
      real*8, intent(in) :: nu,v4(n,n,n,n,nspin,nspin),
     &                      wexa(0:n_nu_exa,n2,nspin,nspin),
     &                      wspl_geom(0:n_nu_geom+2,4,n2,nspin,nspin),
     &                      wspl_asy(n_nu_asy+1,4,n2,nspin,nspin),
     &                      wasy(n2,nspin,nspin)
      complex*16, intent(out) :: th0(2*n,2*n,2*n,2*n)
      integer :: isp1,isp2,isp3,isp4,is1,is2,is3,is4,i10,i20,i30,i40,nn
      real*8, allocatable :: w4(:,:,:,:)
      nn=n*n
      allocate(w4(n,n,n,n))
      th0=(0.d0,0.d0)
      do isp4=1,2
        is4=min(isp4,nspin)
        do isp3=1,2
          is3=min(isp4,nspin)
          do isp2=1,2
            is2=min(isp4,nspin)
            do isp1=1,2
              is1=min(isp4,nspin)
              if(isp1/=isp3) cycle
              if(isp2/=isp4) cycle
              i40=n*(isp4-1)
              i30=n*(isp3-1)
              i20=n*(isp2-1)
              i10=n*(isp1-1)
              th0(i10+1:i10+n,i20+1:i20+n,i30+1:i30+n,i40+1:i40+n)=
     &          th0(i10+1:i10+n,i20+1:i20+n,i30+1:i30+n,i40+1:i40+n)
     &         -v4(:,:,:,:,is1,is2)
              call interp_w_nu_spl(wexa(0,1,is1,is2),
     &                             wspl_geom(0,1,1,is1,is2),
     &                             wspl_asy(1,1,1,is1,is2),
     &                             wasy(1,is1,is2),w4,nu,n2)
              th0(i10+1:i10+n,i20+1:i20+n,i30+1:i30+n,i40+1:i40+n)=
     &          th0(i10+1:i10+n,i20+1:i20+n,i30+1:i30+n,i40+1:i40+n)
     &         -w4(:,:,:,:)
            enddo
          enddo
        enddo
      enddo
      deallocate(w4)
      end
      
      
      
      subroutine kernel_loc_1_2(n,n2,th0,o1,o2,nu,v4,w4,wexa,wasy,
     &                          wspl_geom,wspl_asy,gt,
     &                          trn,gexa,gspl_geom,gspl_asy,gasy)
      use atom_mod
      use etot_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use units_mod
      use vertex_mod
      implicit none
      integer, intent(in) :: n,n2
      real*8, intent(in) :: o1,o2,nu,v4(n,n,n,n,nspin,nspin),
     &                      wexa(0:n_nu_exa,n2,nspin,nspin),
     &                      wspl_geom(0:n_nu_geom+2,4,n2,nspin,nspin),
     &                      wspl_asy(n_nu_asy+1,4,n2,nspin,nspin),
     &                      wasy(n2,nspin,nspin),
     &                      gt(n,n,0:n_tau,nspin),trn(0:n_nu),
     &                      gexa(0:n_omega_exa,2,n,n,nspin),
     &                      gspl_geom(0:n_omega_geom+2,4,2,n,n,nspin),
     &                      gspl_asy(n_omega_asy+1,4,2,n,n,nspin),
     &                      gasy(2,n,n,nspin),
     &                      w4(n,n,n,n,0:n_nu,nspin,nspin)
      complex*16, intent(out) :: th0(2*n,2*n,2*n,2*n)
      integer :: isp1,isp2,isp3,isp4,is1,is2,is3,is4,i10,i20,i30,i40,nn,
     &           nud,it,i,j,k,l,ind_tau,i_tau,i1_tau,ind_nu,i_nu
      real*8 :: pi2b,a,b
      complex*16 :: cc
      real*8, allocatable :: ww(:,:,:,:),vv(:,:,:,:),gg(:,:,:,:,:,:)
      complex*16, allocatable :: wgg(:,:,:,:),g1(:,:),g2(:,:),
     &                           wg(:,:,:,:),wggg(:,:,:,:)
      pi2b=(pi+pi)/betta_t
      nn=n*n
      allocate(ww(n,n,n,n))
      allocate(vv(n,n,n,n))
      th0=(0.d0,0.d0)
c ------------- Static part ---------------------------------------
      allocate(gg(n,n,n,n,ndim3_tau,2))
      nud=idnint((o1-o2)/pi2b)
      cc=(0.d0,1.d0)
      if(nud<0) cc=(0.d0,-1.d0)
      nud=iabs(nud)
      do isp4=1,2
        is4=min(isp4,nspin)
        do isp3=1,2
          is3=min(isp4,nspin)
          do isp2=1,2
            is2=min(isp4,nspin)
            do isp1=1,2
              is1=min(isp4,nspin)
              if(isp1/=isp2) cycle
              if(isp3/=isp4) cycle
              i40=n*(isp4-1)
              i30=n*(isp3-1)
              i20=n*(isp2-1)
              i10=n*(isp1-1)
              do l=1,n
                do k=1,n
                  do j=1,n
                    do i=1,n
                      vv(i,j,k,l)=v4(i,k,j,l,is1,is3)
                    enddo
                  enddo
                enddo
              enddo
              do ind_tau=1,ndim3_tau
                i_tau=me_t*ndim3_tau+ind_tau-1
                i1_tau=n_tau-i_tau
                do l=1,n
                  do k=1,n
                    do j=1,n
                      do i=1,n
                        a=gt(i,j,i_tau,is1)*gt(k,l,i1_tau,is3)
                        b=gt(i,j,i1_tau,is1)*gt(k,l,i_tau,is3)
                        gg(i,l,j,k,ind_tau,1)=a+b
                        gg(i,l,j,k,ind_tau,2)=a-b
                      enddo
                    enddo
                  enddo
                enddo
                do it=1,2
	            call dgemm('n','n',nn,nn,nn,1.d0,vv,nn,
     &                       gg(1,1,1,1,ind_tau,it),nn,0.d0,ww,nn)
	            call dgemm('n','t',nn,nn,nn,1.d0,ww,nn,vv,nn,0.d0,
     &                       gg(1,1,1,1,ind_tau,it),nn)
                enddo
              enddo   !! over ind_tau
              do l=1,n
                do k=1,n
                  do j=1,n
                    do i=1,n
                      a=dot_product(gg(i,k,j,l,:,1),
     &                              ab_nu_from_tau(:,nud,1))
                      b=dot_product(gg(i,k,j,l,:,2),
     &                              ab_nu_from_tau(:,nud,2))
                      th0(i10+i,i20+j,i30+k,i40+l)=
     &                  th0(i10+i,i20+j,i30+k,i40+l)+a+cc*b
                    enddo
                  enddo
                enddo
              enddo
            enddo
          enddo
        enddo
      enddo
      deallocate(gg)
	if(nproc_t/=1) call dgop(th0,32*n**4,'  +',comm_t)
c ------------- Semidynamic part ---------------------------------------
      allocate(wggg(n,n,n,n))
      allocate(wgg(n,n,n,n))
      allocate(wg(n,n,n,n))
      allocate(g1(n,n))
      allocate(g2(n,n))
      do isp4=1,2
        is4=min(isp4,nspin)
        do isp3=1,2
          is3=min(isp4,nspin)
          do isp2=1,2
            is2=min(isp4,nspin)
            do isp1=1,2
              is1=min(isp4,nspin)
              if(isp1/=isp2) cycle
              if(isp3/=isp4) cycle
              i40=n*(isp4-1)
              i30=n*(isp3-1)
              i20=n*(isp2-1)
              i10=n*(isp1-1)
              do l=1,n
                do k=1,n
                  do j=1,n
                    do i=1,n
                      vv(i,j,k,l)=v4(i,k,j,l,is1,is3)
                    enddo
                  enddo
                enddo
              enddo
              wggg=(0.d0,0.d0)
              do ind_nu=1,ndim3_nu
                i_nu=me_t*ndim3_nu+ind_nu-1
                call interp_g_spl(gexa(0,1,1,1,is1),
     &                            gspl_geom(0,1,1,1,1,is1),
     &                            gspl_asy(1,1,1,1,1,is1),
     &                            gasy(1,1,1,is1),g1,o1+w_nu(i_nu),nn,1)
                call interp_g_spl(gexa(0,1,1,1,is3),
     &                            gspl_geom(0,1,1,1,1,is3),
     &                            gspl_asy(1,1,1,1,1,is3),
     &                            gasy(1,1,1,is3),g2,o2+w_nu(i_nu),nn,1)
	          call dgemm('n','n',2*n,nn*n,n,1.d0,g1,2*n,
     &                     w4(1,1,1,1,i_nu,is1,is3),n,0.d0,wg,2*n)
	          call zgemm('n','t',nn*n,n,n,(1.d0,0.d0),wg,nn*n,
     &                     g2,n,(0.d0,0.d0),wgg,nn*n)
                if(i_nu/=0) then
                  call interp_g_spl(gexa(0,1,1,1,is1),
     &                              gspl_geom(0,1,1,1,1,is1),
     &                              gspl_asy(1,1,1,1,1,is1),
     &                              gasy(1,1,1,is1),g1,o1-w_nu(i_nu),
     &                              nn,1)
                  call interp_g_spl(gexa(0,1,1,1,is3),
     &                              gspl_geom(0,1,1,1,1,is3),
     &                              gspl_asy(1,1,1,1,1,is3),
     &                              gasy(1,1,1,is3),g2,o2-w_nu(i_nu),
     &                              nn,1)
	            call dgemm('n','n',2*n,nn*n,n,1.d0,g1,2*n,
     &                       w4(1,1,1,1,i_nu,is1,is3),n,0.d0,wg,2*n)
	            call zgemm('n','t',nn*n,n,n,(1.d0,0.d0),wg,nn*n,
     &                       g2,n,(1.d0,0.d0),wgg,nn*n)
                endif
                do l=1,n
                  do k=1,n
                    do j=1,n
                      do i=1,n
                        wggg(j,k,i,l)=wggg(j,k,i,l)+trn(i_nu)
     &                                             *wgg(i,j,k,l)
                      enddo
                    enddo
                  enddo
                enddo
              enddo   !! over ind_nu
	        if(nproc_t/=1) call dgop(wggg,2*n**4,'  +',comm_t)
	        call dgemm('n','t',2*nn,nn,nn,1.d0,wggg,2*nn,
     &                   vv,nn,0.d0,wg,2*nn)
              do l=1,n
                do k=1,n
                  do j=1,n
                    do i=1,n
                      th0(i10+i,i20+j,i30+k,i40+l)=
     &                  th0(i10+i,i20+j,i30+k,i40+l)-wg(i,k,j,l)
                    enddo
                  enddo
                enddo
              enddo
              wggg=(0.d0,0.d0)
              do ind_nu=1,ndim3_nu
                i_nu=me_t*ndim3_nu+ind_nu-1
                call interp_g_spl(gexa(0,1,1,1,is1),
     &                            gspl_geom(0,1,1,1,1,is1),
     &                            gspl_asy(1,1,1,1,1,is1),
     &                            gasy(1,1,1,is1),g1,o1-nu+w_nu(i_nu),
     &                            nn,1)
                call interp_g_spl(gexa(0,1,1,1,is3),
     &                            gspl_geom(0,1,1,1,1,is3),
     &                            gspl_asy(1,1,1,1,1,is3),
     &                            gasy(1,1,1,is3),g2,o2-nu+w_nu(i_nu),
     &                            nn,1)
	          call dgemm('n','n',2*n,nn*n,n,1.d0,g1,2*n,
     &                     w4(1,1,1,1,i_nu,is1,is3),n,0.d0,wg,2*n)
	          call zgemm('n','n',nn*n,n,n,(1.d0,0.d0),wg,nn*n,
     &                     g2,n,(0.d0,0.d0),wgg,nn*n)
                if(i_nu/=0) then
                  call interp_g_spl(gexa(0,1,1,1,is1),
     &                              gspl_geom(0,1,1,1,1,is1),
     &                              gspl_asy(1,1,1,1,1,is1),
     &                              gasy(1,1,1,is1),g1,o1-nu-w_nu(i_nu),
     &                              nn,1)
                  call interp_g_spl(gexa(0,1,1,1,is3),
     &                              gspl_geom(0,1,1,1,1,is3),
     &                              gspl_asy(1,1,1,1,1,is3),
     &                              gasy(1,1,1,is3),g2,o2-nu-w_nu(i_nu),
     &                              nn,1)
	            call dgemm('n','n',2*n,nn*n,n,1.d0,g1,2*n,
     &                       w4(1,1,1,1,i_nu,is1,is3),n,0.d0,wg,2*n)
	            call zgemm('n','t',nn*n,n,n,(1.d0,0.d0),wg,nn*n,
     &                       g2,n,(1.d0,0.d0),wgg,nn*n)
                endif
                do l=1,n
                  do k=1,n
                    do j=1,n
                      do i=1,n
                        wggg(j,k,i,l)=wggg(j,k,i,l)+trn(i_nu)
     &                                             *wgg(i,j,k,l)
                      enddo
                    enddo
                  enddo
                enddo
              enddo   !! over ind_nu
	        if(nproc_t/=1) call dgop(wggg,2*n**4,'  +',comm_t)
	        call dgemm('n','t',2*nn,nn,nn,1.d0,wggg,2*nn,
     &                   vv,nn,0.d0,wg,2*nn)
              do l=1,n
                do k=1,n
                  do j=1,n
                    do i=1,n
                      th0(i10+i,i20+j,i30+k,i40+l)=
     &                  th0(i10+i,i20+j,i30+k,i40+l)-wg(j,l,i,k)
                    enddo
                  enddo
                enddo
              enddo           
            enddo
          enddo
        enddo
      enddo
      deallocate(wgg,g1,g2,wg,wggg)
c ------------- Dynamic part ---------------------------------------
      allocate(wggg(n,n,n,n))
      allocate(wgg(n,n,n,n))
      allocate(wg(n,n,n,n))
      allocate(g1(n,n))
      allocate(g2(n,n))
      do isp4=1,2
        is4=min(isp4,nspin)
        do isp3=1,2
          is3=min(isp4,nspin)
          do isp2=1,2
            is2=min(isp4,nspin)
            do isp1=1,2
              is1=min(isp4,nspin)
              if(isp1/=isp2) cycle
              if(isp3/=isp4) cycle
              i40=n*(isp4-1)
              i30=n*(isp3-1)
              i20=n*(isp2-1)
              i10=n*(isp1-1)
              wggg=(0.d0,0.d0)
              do ind_nu=1,ndim3_nu
                i_nu=me_t*ndim3_nu+ind_nu-1
                call interp_g_spl(gexa(0,1,1,1,is1),
     &                            gspl_geom(0,1,1,1,1,is1),
     &                            gspl_asy(1,1,1,1,1,is1),
     &                            gasy(1,1,1,is1),g1,o1+w_nu(i_nu),
     &                            nn,1)
                call interp_g_spl(gexa(0,1,1,1,is3),
     &                            gspl_geom(0,1,1,1,1,is3),
     &                            gspl_asy(1,1,1,1,1,is3),
     &                            gasy(1,1,1,is3),g2,o2+w_nu(i_nu),
     &                            nn,1)
                call interp_w_nu_spl(wexa(0,1,is1,is3),
     &                               wspl_geom(0,1,1,is1,is3),
     &                               wspl_asy(1,1,1,is1,is3),
     &                               wasy(1,is1,is3),ww,-w_nu(i_nu)-nu,
     &                               n2)
	          call dgemm('n','n',2*n,nn*n,n,1.d0,g1,2*n,ww,n,0.d0,wg,
     &                     2*n)
	          call zgemm('n','n',nn*n,n,n,(1.d0,0.d0),wg,nn*n,
     &                     g2,n,(0.d0,0.d0),wgg,nn*n)
                if(i_nu/=0) then
                  call interp_g_spl(gexa(0,1,1,1,is1),
     &                              gspl_geom(0,1,1,1,1,is1),
     &                              gspl_asy(1,1,1,1,1,is1),
     &                              gasy(1,1,1,is1),g1,o1-w_nu(i_nu),
     &                              nn,1)
                  call interp_g_spl(gexa(0,1,1,1,is3),
     &                              gspl_geom(0,1,1,1,1,is3),
     &                              gspl_asy(1,1,1,1,1,is3),
     &                              gasy(1,1,1,is3),g2,o2-w_nu(i_nu),
     &                              nn,1)
                  call interp_w_nu_spl(wexa(0,1,is1,is3),
     &                                 wspl_geom(0,1,1,is1,is3),
     &                                 wspl_asy(1,1,1,is1,is3),
     &                                 wasy(1,is1,is3),ww,w_nu(i_nu)-nu,
     &                                 n2)
	            call dgemm('n','n',2*n,nn*n,n,1.d0,g1,2*n,ww,n,0.d0,
     &                       wg,2*n)
	            call zgemm('n','t',nn*n,n,n,(1.d0,0.d0),wg,nn*n,
     &                       g2,n,(1.d0,0.d0),wgg,nn*n)
                endif
                do l=1,n
                  do k=1,n
                    do j=1,n
                      do i=1,n
                        wg(j,k,i,l)=wgg(i,j,k,l)
                        ww(i,k,j,l)=w4(i,j,k,l,i_nu,is1,is3)
                      enddo
                    enddo
                  enddo
                enddo
	          call dgemm('n','t',2*nn,nn,nn,1.d0,wg,2*nn,ww,nn,0.d0,
     &                     wgg,2*nn)
                do l=1,n
                  do k=1,n
                    do j=1,n
                      do i=1,n
                        wggg(k,i,l,j)=wggg(k,i,l,j)+trn(i_nu)
     &                                             *wgg(i,j,k,l)
                      enddo
                    enddo
                  enddo
                enddo
              enddo   !! over ind_nu
	        if(nproc_t/=1) call dgop(wggg,2*n**4,'  +',comm_t)
              do l=1,n
                do k=1,n
                  do j=1,n
                    do i=1,n
                      th0(i10+i,i20+j,i30+k,i40+l)=
     &                  th0(i10+i,i20+j,i30+k,i40+l)-wggg(i,j,k,l)
                    enddo
                  enddo
                enddo
              enddo
            enddo
          enddo
        enddo
      enddo
      deallocate(wgg,g1,g2,wg,wggg)
      deallocate(ww,vv)
      end
      
      
      
      subroutine kernel_loc_2_2(n,n2,th0,o1,o2,nu,v4,w4,wexa,wasy,
     &                          wspl_geom,wspl_asy,gt,
     &                          trn,gexa,gspl_geom,gspl_asy,gasy)
      use atom_mod
      use etot_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use units_mod
      use vertex_mod
      implicit none
      integer, intent(in) :: n,n2
      real*8, intent(in) :: o1,o2,nu,v4(n,n,n,n,nspin,nspin),
     &                      wexa(0:n_nu_exa,n2,nspin,nspin),
     &                      wspl_geom(0:n_nu_geom+2,4,n2,nspin,nspin),
     &                      wspl_asy(n_nu_asy+1,4,n2,nspin,nspin),
     &                      wasy(n2,nspin,nspin),
     &                      gt(n,n,0:n_tau,nspin),trn(0:n_nu),
     &                      gexa(0:n_omega_exa,2,n,n,nspin),
     &                      gspl_geom(0:n_omega_geom+2,4,2,n,n,nspin),
     &                      gspl_asy(n_omega_asy+1,4,2,n,n,nspin),
     &                      gasy(2,n,n,nspin),
     &                      w4(n,n,n,n,0:n_nu,nspin,nspin)
      complex*16, intent(out) :: th0(2*n,2*n,2*n,2*n)
      integer :: isp1,isp2,isp3,isp4,is1,is2,is3,is4,i10,i20,i30,i40,nn,
     &           nud,it,i,j,k,l,ind_tau,i_tau,i1_tau,ind_nu,i_nu
      real*8 :: pi2b,a,b
      complex*16 :: cc
      real*8, allocatable :: ww(:,:,:,:),vv(:,:,:,:),gg(:,:,:,:,:,:)
      complex*16, allocatable :: wgg(:,:,:,:),g1(:,:),g2(:,:),
     &                           wg(:,:,:,:),wggg(:,:,:,:)
      pi2b=(pi+pi)/betta_t
      nn=n*n
      allocate(ww(n,n,n,n))
      allocate(vv(n,n,n,n))
      th0=(0.d0,0.d0)
c ------------- Static part ---------------------------------------
      allocate(gg(n,n,n,n,ndim3_tau,2))
      nud=idnint((o1+o2-nu)/pi2b)
      cc=(0.d0,1.d0)
      if(nud<0) cc=(0.d0,-1.d0)
      nud=iabs(nud)
      do isp4=1,2
        is4=min(isp4,nspin)
        do isp3=1,2
          is3=min(isp4,nspin)
          do isp2=1,2
            is2=min(isp4,nspin)
            do isp1=1,2
              is1=min(isp4,nspin)
              if(isp1/=isp2) cycle
              if(isp3/=isp4) cycle
              i40=n*(isp4-1)
              i30=n*(isp3-1)
              i20=n*(isp2-1)
              i10=n*(isp1-1)
              do l=1,n
                do k=1,n
                  do j=1,n
                    do i=1,n
                      vv(i,j,k,l)=v4(i,k,j,l,is1,is3)
                    enddo
                  enddo
                enddo
              enddo
              do ind_tau=1,ndim3_tau
                i_tau=me_t*ndim3_tau+ind_tau-1
                i1_tau=n_tau-i_tau
                do l=1,n
                  do k=1,n
                    do j=1,n
                      do i=1,n
                        a=gt(i,j,i_tau,is1)*gt(k,l,i_tau,is3)
                        b=gt(i,j,i1_tau,is1)*gt(k,l,i1_tau,is3)
                        gg(i,k,j,l,ind_tau,1)=a+b
                        gg(i,k,j,l,ind_tau,2)=a-b
                      enddo
                    enddo
                  enddo
                enddo
                do it=1,2
	            call dgemm('n','n',nn,nn,nn,1.d0,vv,nn,
     &                       gg(1,1,1,1,ind_tau,it),nn,0.d0,ww,nn)
	            call dgemm('n','t',nn,nn,nn,1.d0,ww,nn,vv,nn,0.d0,
     &                       gg(1,1,1,1,ind_tau,it),nn)
                enddo
              enddo   !! over ind_tau
              do l=1,n
                do k=1,n
                  do j=1,n
                    do i=1,n
                      a=dot_product(gg(i,l,j,k,:,1),
     &                              ab_nu_from_tau(:,nud,1))
                      b=dot_product(gg(i,l,j,k,:,2),
     &                              ab_nu_from_tau(:,nud,2))
                      th0(i10+i,i20+j,i30+k,i40+l)=
     &                  th0(i10+i,i20+j,i30+k,i40+l)-a-cc*b
                    enddo
                  enddo
                enddo
              enddo
            enddo
          enddo
        enddo
      enddo
      deallocate(gg)
	if(nproc_t/=1) call dgop(th0,32*n**4,'  +',comm_t)
c ------------- Semidynamic part ---------------------------------------
      allocate(wggg(n,n,n,n))
      allocate(wgg(n,n,n,n))
      allocate(wg(n,n,n,n))
      allocate(g1(n,n))
      allocate(g2(n,n))
      do isp4=1,2
        is4=min(isp4,nspin)
        do isp3=1,2
          is3=min(isp4,nspin)
          do isp2=1,2
            is2=min(isp4,nspin)
            do isp1=1,2
              is1=min(isp4,nspin)
              if(isp1/=isp2) cycle
              if(isp3/=isp4) cycle
              i40=n*(isp4-1)
              i30=n*(isp3-1)
              i20=n*(isp2-1)
              i10=n*(isp1-1)
              do l=1,n
                do k=1,n
                  do j=1,n
                    do i=1,n
                      vv(i,j,k,l)=v4(i,k,j,l,is1,is3)
                    enddo
                  enddo
                enddo
              enddo
              wggg=(0.d0,0.d0)
              do ind_nu=1,ndim3_nu
                i_nu=me_t*ndim3_nu+ind_nu-1
                call interp_g_spl(gexa(0,1,1,1,is1),
     &                            gspl_geom(0,1,1,1,1,is1),
     &                            gspl_asy(1,1,1,1,1,is1),
     &                            gasy(1,1,1,is1),g1,o1+w_nu(i_nu),
     &                            nn,1)
                call interp_g_spl(gexa(0,1,1,1,is3),
     &                            gspl_geom(0,1,1,1,1,is3),
     &                            gspl_asy(1,1,1,1,1,is3),
     &                            gasy(1,1,1,is3),g2,o2-nu-w_nu(i_nu),
     &                            nn,1)
	          call dgemm('n','n',2*n,nn*n,n,1.d0,g1,2*n,
     &                     w4(1,1,1,1,i_nu,is1,is3),n,0.d0,wg,2*n)
	          call zgemm('n','t',nn*n,n,n,(1.d0,0.d0),wg,nn*n,
     &                     g2,n,(0.d0,0.d0),wgg,nn*n)
                if(i_nu/=0) then
                  call interp_g_spl(gexa(0,1,1,1,is1),
     &                              gspl_geom(0,1,1,1,1,is1),
     &                              gspl_asy(1,1,1,1,1,is1),
     &                              gasy(1,1,1,is1),g1,o1-w_nu(i_nu),
     &                              nn,1)
                  call interp_g_spl(gexa(0,1,1,1,is3),
     &                              gspl_geom(0,1,1,1,1,is3),
     &                              gspl_asy(1,1,1,1,1,is3),
     &                              gasy(1,1,1,is3),g2,o2-nu+w_nu(i_nu),
     &                              nn,1)
	            call dgemm('n','n',2*n,nn*n,n,1.d0,g1,2*n,
     &                       w4(1,1,1,1,i_nu,is1,is3),n,0.d0,wg,2*n)
	            call zgemm('n','t',nn*n,n,n,(1.d0,0.d0),wg,nn*n,
     &                       g2,n,(1.d0,0.d0),wgg,nn*n)
                endif
                do l=1,n
                  do k=1,n
                    do j=1,n
                      do i=1,n
                        wggg(j,k,i,l)=wggg(j,k,i,l)+trn(i_nu)
     &                                             *wgg(i,j,k,l)
                      enddo
                    enddo
                  enddo
                enddo
              enddo   !! over ind_nu
	        if(nproc_t/=1) call dgop(wggg,2*n**4,'  +',comm_t)
	        call dgemm('n','t',2*nn,nn,nn,1.d0,wggg,2*nn,
     &                   vv,nn,0.d0,wg,2*nn)
              do l=1,n
                do k=1,n
                  do j=1,n
                    do i=1,n
                      th0(i10+i,i20+j,i30+k,i40+l)=
     &                  th0(i10+i,i20+j,i30+k,i40+l)-wg(i,l,j,k)
                    enddo
                  enddo
                enddo
              enddo
              wggg=(0.d0,0.d0)
              do ind_nu=1,ndim3_nu
                i_nu=me_t*ndim3_nu+ind_nu-1
                call interp_g_spl(gexa(0,1,1,1,is1),
     &                            gspl_geom(0,1,1,1,1,is1),
     &                            gspl_asy(1,1,1,1,1,is1),
     &                            gasy(1,1,1,is1),g1,o2-w_nu(i_nu),
     &                            nn,1)
                call interp_g_spl(gexa(0,1,1,1,is3),
     &                            gspl_geom(0,1,1,1,1,is3),
     &                            gspl_asy(1,1,1,1,1,is3),
     &                            gasy(1,1,1,is3),g2,o1-nu+w_nu(i_nu),
     &                            nn,1)
	          call dgemm('n','n',2*n,nn*n,n,1.d0,g1,2*n,
     &                     w4(1,1,1,1,i_nu,is1,is3),n,0.d0,wg,2*n)
	          call zgemm('n','n',nn*n,n,n,(1.d0,0.d0),wg,nn*n,
     &                     g2,n,(0.d0,0.d0),wgg,nn*n)
                if(i_nu/=0) then
                  call interp_g_spl(gexa(0,1,1,1,is1),
     &                              gspl_geom(0,1,1,1,1,is1),
     &                              gspl_asy(1,1,1,1,1,is1),
     &                              gasy(1,1,1,is1),g1,o2+w_nu(i_nu),
     &                              nn,1)
                  call interp_g_spl(gexa(0,1,1,1,is3),
     &                              gspl_geom(0,1,1,1,1,is3),
     &                              gspl_asy(1,1,1,1,1,is3),
     &                              gasy(1,1,1,is3),g2,o1-nu-w_nu(i_nu),
     &                              nn,1)
	            call dgemm('n','n',2*n,nn*n,n,1.d0,g1,2*n,
     &                       w4(1,1,1,1,i_nu,is1,is3),n,0.d0,wg,2*n)
	            call zgemm('n','t',nn*n,n,n,(1.d0,0.d0),wg,nn*n,
     &                       g2,n,(1.d0,0.d0),wgg,nn*n)
                endif
                do l=1,n
                  do k=1,n
                    do j=1,n
                      do i=1,n
                        wggg(j,k,i,l)=wggg(j,k,i,l)+trn(i_nu)
     &                                             *wgg(i,j,k,l)
                      enddo
                    enddo
                  enddo
                enddo
              enddo   !! over ind_nu
	        if(nproc_t/=1) call dgop(wggg,2*n**4,'  +',comm_t)
	        call dgemm('n','t',2*nn,nn,nn,1.d0,wggg,2*nn,
     &                   vv,nn,0.d0,wg,2*nn)
              do l=1,n
                do k=1,n
                  do j=1,n
                    do i=1,n
                      th0(i10+i,i20+j,i30+k,i40+l)=
     &                  th0(i10+i,i20+j,i30+k,i40+l)-wg(j,k,i,l)
                    enddo
                  enddo
                enddo
              enddo           
            enddo
          enddo
        enddo
      enddo
      deallocate(wgg,g1,g2,wg,wggg)
c ------------- Dynamic part ---------------------------------------
      allocate(wggg(n,n,n,n))
      allocate(wgg(n,n,n,n))
      allocate(wg(n,n,n,n))
      allocate(g1(n,n))
      allocate(g2(n,n))
      do isp4=1,2
        is4=min(isp4,nspin)
        do isp3=1,2
          is3=min(isp4,nspin)
          do isp2=1,2
            is2=min(isp4,nspin)
            do isp1=1,2
              is1=min(isp4,nspin)
              if(isp1/=isp2) cycle
              if(isp3/=isp4) cycle
              i40=n*(isp4-1)
              i30=n*(isp3-1)
              i20=n*(isp2-1)
              i10=n*(isp1-1)
              wggg=(0.d0,0.d0)
              do ind_nu=1,ndim3_nu
                i_nu=me_t*ndim3_nu+ind_nu-1
                call interp_g_spl(gexa(0,1,1,1,is1),
     &                            gspl_geom(0,1,1,1,1,is1),
     &                            gspl_asy(1,1,1,1,1,is1),
     &                            gasy(1,1,1,is1),g1,o1+w_nu(i_nu),
     &                            nn,1)
                call interp_g_spl(gexa(0,1,1,1,is3),
     &                            gspl_geom(0,1,1,1,1,is3),
     &                            gspl_asy(1,1,1,1,1,is3),
     &                            gasy(1,1,1,is3),g2,o2-nu-w_nu(i_nu),
     &                            nn,1)
                call interp_w_nu_spl(wexa(0,1,is1,is3),
     &                               wspl_geom(0,1,1,is1,is3),
     &                               wspl_asy(1,1,1,is1,is3),
     &                               wasy(1,is1,is3),ww,
     &                               o1-o2+w_nu(i_nu),n2)
	          call dgemm('n','n',2*n,nn*n,n,1.d0,g1,2*n,ww,n,0.d0,wg,
     &                     2*n)
	          call zgemm('n','n',nn*n,n,n,(1.d0,0.d0),wg,nn*n,
     &                     g2,n,(0.d0,0.d0),wgg,nn*n)
                if(i_nu/=0) then
                  call interp_g_spl(gexa(0,1,1,1,is1),
     &                              gspl_geom(0,1,1,1,1,is1),
     &                              gspl_asy(1,1,1,1,1,is1),
     &                              gasy(1,1,1,is1),g1,o1-w_nu(i_nu),
     &                              nn,1)
                  call interp_g_spl(gexa(0,1,1,1,is3),
     &                              gspl_geom(0,1,1,1,1,is3),
     &                              gspl_asy(1,1,1,1,1,is3),
     &                              gasy(1,1,1,is3),g2,o2-nu+w_nu(i_nu),
     &                              nn,1)
                  call interp_w_nu_spl(wexa(0,1,is1,is3),
     &                                 wspl_geom(0,1,1,is1,is3),
     &                                 wspl_asy(1,1,1,is1,is3),
     &                                 wasy(1,is1,is3),ww,
     &                                 o1-o2-w_nu(i_nu),n2)
	            call dgemm('n','n',2*n,nn*n,n,1.d0,g1,2*n,ww,n,0.d0,
     &                       wg,2*n)
	            call zgemm('n','t',nn*n,n,n,(1.d0,0.d0),wg,nn*n,
     &                       g2,n,(1.d0,0.d0),wgg,nn*n)
                endif
                do l=1,n
                  do k=1,n
                    do j=1,n
                      do i=1,n
                        wg(j,k,i,l)=wgg(i,j,k,l)
                        ww(j,l,i,k)=w4(i,j,k,l,i_nu,is1,is3)
                      enddo
                    enddo
                  enddo
                enddo
	          call dgemm('n','n',2*nn,nn,nn,1.d0,wg,2*nn,ww,nn,0.d0,
     &                     wgg,2*nn)
                do l=1,n
                  do k=1,n
                    do j=1,n
                      do i=1,n
                        wggg(k,i,j,l)=wggg(k,i,j,l)+trn(i_nu)
     &                                             *wgg(i,j,k,l)
                      enddo
                    enddo
                  enddo
                enddo
              enddo   !! over ind_nu
	        if(nproc_t/=1) call dgop(wggg,2*n**4,'  +',comm_t)
              do l=1,n
                do k=1,n
                  do j=1,n
                    do i=1,n
                      th0(i10+i,i20+j,i30+k,i40+l)=
     &                  th0(i10+i,i20+j,i30+k,i40+l)-wggg(i,j,k,l)
                    enddo
                  enddo
                enddo
              enddo
            enddo
          enddo
        enddo
      enddo
      deallocate(wgg,g1,g2,wg,wggg)
      deallocate(ww,vv)
      end
      
      
      
      subroutine kernel_loc_3_2(n,n2,th0,o1,o2,nu,v4,w4,wexa,wasy,
     &                          wspl_geom,wspl_asy,gt,
     &                          trn,gexa,gspl_geom,gspl_asy,gasy)
      use atom_mod
      use etot_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use units_mod
      use vertex_mod
      implicit none
      integer, intent(in) :: n,n2
      real*8, intent(in) :: o1,o2,nu,v4(n,n,n,n,nspin,nspin),
     &                      wexa(0:n_nu_exa,n2,nspin,nspin),
     &                      wspl_geom(0:n_nu_geom+2,4,n2,nspin,nspin),
     &                      wspl_asy(n_nu_asy+1,4,n2,nspin,nspin),
     &                      wasy(n2,nspin,nspin),
     &                      gt(n,n,0:n_tau,nspin),trn(0:n_nu),
     &                      gexa(0:n_omega_exa,2,n,n,nspin),
     &                      gspl_geom(0:n_omega_geom+2,4,2,n,n,nspin),
     &                      gspl_asy(n_omega_asy+1,4,2,n,n,nspin),
     &                      gasy(2,n,n,nspin),
     &                      w4(n,n,n,n,0:n_nu,nspin,nspin)
      complex*16, intent(out) :: th0(2*n,2*n,2*n,2*n)
      integer :: isp1,isp2,isp3,isp4,is1,is2,is3,is4,i10,i20,i30,i40,nn,
     &           nud,it,i,j,k,l,ind_tau,i_tau,i1_tau,ind_nu,i_nu
      real*8 :: pi2b,a,b
      complex*16 :: cc
      real*8, allocatable :: ww(:,:,:,:),vv(:,:,:,:),gg(:,:,:,:,:,:),
     &                       v1(:,:,:,:)
      complex*16, allocatable :: wgg(:,:,:,:),g1(:,:),g2(:,:),
     &                           wg(:,:,:,:),wggg(:,:,:,:)
      pi2b=(pi+pi)/betta_t
      nn=n*n
      allocate(ww(n,n,n,n))
      allocate(vv(n,n,n,n))
      allocate(v1(n,n,n,n))
      th0=(0.d0,0.d0)
c ------------- Static + Semi_2 part --------------------------
      allocate(gg(n,n,n,n,ndim3_tau,2))
      nud=idnint(nu/pi2b)
      cc=(0.d0,1.d0)
      if(nud<0) cc=(0.d0,-1.d0)
      nud=iabs(nud)
      do isp4=1,2
        is4=min(isp4,nspin)
        do isp3=1,2
          is3=min(isp4,nspin)
          do isp2=1,2
            is2=min(isp4,nspin)
            do isp1=1,2
              is1=min(isp4,nspin)
              if(isp1/=isp3) cycle
              if(isp2/=isp4) cycle
              i40=n*(isp4-1)
              i30=n*(isp3-1)
              i20=n*(isp2-1)
              i10=n*(isp1-1)
              call interp_w_nu_spl(wexa(0,1,is1,is2),
     &                             wspl_geom(0,1,1,is1,is2),
     &                             wspl_asy(1,1,1,is1,is2),
     &                             wasy(1,is1,is2),ww,-nu,n2)
              do l=1,n
                do k=1,n
                  do j=1,n
                    do i=1,n
                      vv(i,j,k,l)=v4(i,k,j,l,is2,is2)
                      ww(i,j,k,l)=ww(i,j,k,l)+v4(i,j,k,l,is1,is2)
                    enddo
                  enddo
                enddo
              enddo
              do ind_tau=1,ndim3_tau
                i_tau=me_t*ndim3_tau+ind_tau-1
                i1_tau=n_tau-i_tau
                do l=1,n
                  do k=1,n
                    do j=1,n
                      do i=1,n
                        a=gt(i,j,i_tau,is2)*gt(k,l,i1_tau,is2)
                        b=gt(i,j,i1_tau,is2)*gt(k,l,i_tau,is2)
                        gg(i,l,j,k,ind_tau,1)=a+b
                        gg(i,l,j,k,ind_tau,2)=a-b
                      enddo
                    enddo
                  enddo
                enddo
                do it=1,2
	            call dgemm('n','t',nn,nn,nn,1.d0,vv,nn,
     &                       gg(1,1,1,1,ind_tau,it),nn,0.d0,v1,nn)
	            call dgemm('n','t',nn,nn,nn,1.d0,ww,nn,v1,nn,0.d0,
     &                       gg(1,1,1,1,ind_tau,it),nn)
                enddo
              enddo   !! over ind_tau
              do l=1,n
                do k=1,n
                  do j=1,n
                    do i=1,n
                      a=dot_product(gg(i,k,j,l,:,1),
     &                              ab_nu_from_tau(:,nud,1))
                      b=dot_product(gg(i,k,j,l,:,2),
     &                              ab_nu_from_tau(:,nud,2))
                      th0(i10+i,i20+j,i30+k,i40+l)=
     &                  th0(i10+i,i20+j,i30+k,i40+l)-a-cc*b
                    enddo
                  enddo
                enddo
              enddo
            enddo
          enddo
        enddo
      enddo
      deallocate(gg)
	if(nproc_t/=1) call dgop(th0,32*n**4,'  +',comm_t)
c ------------- Semi_1 + Dynamic part ----------------------------
      allocate(wggg(n,n,n,n))
      allocate(wgg(n,n,n,n))
      allocate(wg(n,n,n,n))
      allocate(g1(n,n))
      allocate(g2(n,n))
      do isp4=1,2
        is4=min(isp4,nspin)
        do isp3=1,2
          is3=min(isp4,nspin)
          do isp2=1,2
            is2=min(isp4,nspin)
            do isp1=1,2
              is1=min(isp4,nspin)
              if(isp1/=isp3) cycle
              if(isp2/=isp4) cycle
              i40=n*(isp4-1)
              i30=n*(isp3-1)
              i20=n*(isp2-1)
              i10=n*(isp1-1)
              call interp_w_nu_spl(wexa(0,1,is1,is2),
     &                             wspl_geom(0,1,1,is1,is2),
     &                             wspl_asy(1,1,1,is1,is2),
     &                             wasy(1,is1,is2),ww,-nu,n2)
              do l=1,n
                do k=1,n
                  do j=1,n
                    do i=1,n
                      ww(i,j,k,l)=ww(i,j,k,l)+v4(i,j,k,l,is1,is2)
                    enddo
                  enddo
                enddo
              enddo
              wggg=(0.d0,0.d0)
              do ind_nu=1,ndim3_nu
                i_nu=me_t*ndim3_nu+ind_nu-1
                call interp_g_spl(gexa(0,1,1,1,is2),
     &                            gspl_geom(0,1,1,1,1,is2),
     &                            gspl_asy(1,1,1,1,1,is2),
     &                            gasy(1,1,1,is2),g1,o2+w_nu(i_nu),
     &                            nn,1)
                call interp_g_spl(gexa(0,1,1,1,is2),
     &                            gspl_geom(0,1,1,1,1,is2),
     &                            gspl_asy(1,1,1,1,1,is2),
     &                            gasy(1,1,1,is2),g2,o2-nu+w_nu(i_nu),
     &                            nn,1)
	          call dgemm('n','n',2*n,nn*n,n,1.d0,g1,2*n,
     &                     w4(1,1,1,1,i_nu,is2,is2),n,0.d0,wg,
     &                     2*n)
	          call zgemm('n','n',nn*n,n,n,(1.d0,0.d0),wg,nn*n,
     &                     g2,n,(0.d0,0.d0),wgg,nn*n)
                if(i_nu/=0) then
                  call interp_g_spl(gexa(0,1,1,1,is2),
     &                              gspl_geom(0,1,1,1,1,is2),
     &                              gspl_asy(1,1,1,1,1,is2),
     &                              gasy(1,1,1,is2),g1,o2-w_nu(i_nu),
     &                              nn,1)
                  call interp_g_spl(gexa(0,1,1,1,is2),
     &                              gspl_geom(0,1,1,1,1,is2),
     &                              gspl_asy(1,1,1,1,1,is2),
     &                              gasy(1,1,1,is2),g2,o2-nu-w_nu(i_nu),
     &                              nn,1)
	            call dgemm('n','n',2*n,nn*n,n,1.d0,g1,2*n,
     &                       w4(1,1,1,1,i_nu,is2,is2),n,0.d0,wg,
     &                       2*n)
	            call zgemm('n','n',nn*n,n,n,(1.d0,0.d0),wg,nn*n,
     &                       g2,n,(1.d0,0.d0),wgg,nn*n)
                endif
                do l=1,n
                  do k=1,n
                    do j=1,n
                      do i=1,n
                        wg(j,k,i,l)=wgg(i,j,k,l)
                      enddo
                    enddo
                  enddo
                enddo
	          call dgemm('n','t',2*nn,nn,nn,1.d0,wg,2*nn,ww,nn,0.d0,
     &                     wgg,2*nn)
                do l=1,n
                  do k=1,n
                    do j=1,n
                      do i=1,n
                        wggg(k,i,l,j)=wggg(k,i,l,j)+trn(i_nu)
     &                                             *wgg(i,j,k,l)
                      enddo
                    enddo
                  enddo
                enddo
              enddo   !! over ind_nu
	        if(nproc_t/=1) call dgop(wggg,2*n**4,'  +',comm_t)
              do l=1,n
                do k=1,n
                  do j=1,n
                    do i=1,n
                      th0(i10+i,i20+j,i30+k,i40+l)=
     &                  th0(i10+i,i20+j,i30+k,i40+l)+wggg(i,j,k,l)
                    enddo
                  enddo
                enddo
              enddo
            enddo
          enddo
        enddo
      enddo
      deallocate(wgg,g1,g2,wg,wggg)
      deallocate(ww,vv,v1)
      end
      
      
      
      subroutine kernel_loc_4_2(n,n2,th0,o1,o2,nu,v4,w4,wexa,wasy,
     &                          wspl_geom,wspl_asy,gt,
     &                          trn,gexa,gspl_geom,gspl_asy,gasy)
      use atom_mod
      use etot_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use units_mod
      use vertex_mod
      implicit none
      integer, intent(in) :: n,n2
      real*8, intent(in) :: o1,o2,nu,v4(n,n,n,n,nspin,nspin),
     &                      wexa(0:n_nu_exa,n2,nspin,nspin),
     &                      wspl_geom(0:n_nu_geom+2,4,n2,nspin,nspin),
     &                      wspl_asy(n_nu_asy+1,4,n2,nspin,nspin),
     &                      wasy(n2,nspin,nspin),
     &                      gt(n,n,0:n_tau,nspin),trn(0:n_nu),
     &                      gexa(0:n_omega_exa,2,n,n,nspin),
     &                      gspl_geom(0:n_omega_geom+2,4,2,n,n,nspin),
     &                      gspl_asy(n_omega_asy+1,4,2,n,n,nspin),
     &                      gasy(2,n,n,nspin),
     &                      w4(n,n,n,n,0:n_nu,nspin,nspin)
      complex*16, intent(out) :: th0(2*n,2*n,2*n,2*n)
      integer :: isp1,isp2,isp3,isp4,is1,is2,is3,is4,i10,i20,i30,i40,nn,
     &           nud,it,i,j,k,l,ind_tau,i_tau,i1_tau,ind_nu,i_nu
      real*8 :: pi2b,a,b
      complex*16 :: cc
      real*8, allocatable :: ww(:,:,:,:),vv(:,:,:,:),gg(:,:,:,:,:,:),
     &                       v1(:,:,:,:)
      complex*16, allocatable :: wgg(:,:,:,:),g1(:,:),g2(:,:),
     &                           wg(:,:,:,:),wggg(:,:,:,:)
      pi2b=(pi+pi)/betta_t
      nn=n*n
      allocate(ww(n,n,n,n))
      allocate(vv(n,n,n,n))
      allocate(v1(n,n,n,n))
      th0=(0.d0,0.d0)
c ------------- Static + Semi_2 part --------------------------
      allocate(gg(n,n,n,n,ndim3_tau,2))
      nud=idnint(nu/pi2b)
      cc=(0.d0,1.d0)
      if(nud<0) cc=(0.d0,-1.d0)
      nud=iabs(nud)
      do isp4=1,2
        is4=min(isp4,nspin)
        do isp3=1,2
          is3=min(isp4,nspin)
          do isp2=1,2
            is2=min(isp4,nspin)
            do isp1=1,2
              is1=min(isp4,nspin)
              if(isp1/=isp3) cycle
              if(isp2/=isp4) cycle
              i40=n*(isp4-1)
              i30=n*(isp3-1)
              i20=n*(isp2-1)
              i10=n*(isp1-1)
              call interp_w_nu_spl(wexa(0,1,is1,is2),
     &                             wspl_geom(0,1,1,is1,is2),
     &                             wspl_asy(1,1,1,is1,is2),
     &                             wasy(1,is1,is2),ww,nu,n2)
              do l=1,n
                do k=1,n
                  do j=1,n
                    do i=1,n
                      vv(i,j,k,l)=v4(i,k,j,l,is1,is1)
                      ww(i,j,k,l)=ww(i,j,k,l)+v4(i,j,k,l,is1,is2)
                    enddo
                  enddo
                enddo
              enddo
              do ind_tau=1,ndim3_tau
                i_tau=me_t*ndim3_tau+ind_tau-1
                i1_tau=n_tau-i_tau
                do l=1,n
                  do k=1,n
                    do j=1,n
                      do i=1,n
                        a=gt(i,j,i_tau,is2)*gt(k,l,i1_tau,is2)
                        b=gt(i,j,i1_tau,is2)*gt(k,l,i_tau,is2)
                        gg(i,l,k,j,ind_tau,1)=a+b
                        gg(i,l,k,j,ind_tau,2)=a-b
                      enddo
                    enddo
                  enddo
                enddo
                do it=1,2
	            call dgemm('n','n',nn,nn,nn,1.d0,vv,nn,
     &                       gg(1,1,1,1,ind_tau,it),nn,0.d0,v1,nn)
	            call dgemm('n','t',nn,nn,nn,1.d0,v1,nn,ww,nn,0.d0,
     &                       gg(1,1,1,1,ind_tau,it),nn)
                enddo
              enddo   !! over ind_tau
              do l=1,n
                do k=1,n
                  do j=1,n
                    do i=1,n
                      a=dot_product(gg(i,k,j,l,:,1),
     &                              ab_nu_from_tau(:,nud,1))
                      b=dot_product(gg(i,k,j,l,:,2),
     &                              ab_nu_from_tau(:,nud,2))
                      th0(i10+i,i20+j,i30+k,i40+l)=
     &                  th0(i10+i,i20+j,i30+k,i40+l)-a-cc*b
                    enddo
                  enddo
                enddo
              enddo
            enddo
          enddo
        enddo
      enddo
      deallocate(gg)
	if(nproc_t/=1) call dgop(th0,32*n**4,'  +',comm_t)
c ------------- Semi_1 + Dynamic part ----------------------------
      allocate(wggg(n,n,n,n))
      allocate(wgg(n,n,n,n))
      allocate(wg(n,n,n,n))
      allocate(g1(n,n))
      allocate(g2(n,n))
      do isp4=1,2
        is4=min(isp4,nspin)
        do isp3=1,2
          is3=min(isp4,nspin)
          do isp2=1,2
            is2=min(isp4,nspin)
            do isp1=1,2
              is1=min(isp4,nspin)
              if(isp1/=isp3) cycle
              if(isp2/=isp4) cycle
              i40=n*(isp4-1)
              i30=n*(isp3-1)
              i20=n*(isp2-1)
              i10=n*(isp1-1)
              call interp_w_nu_spl(wexa(0,1,is1,is2),
     &                             wspl_geom(0,1,1,is1,is2),
     &                             wspl_asy(1,1,1,is1,is2),
     &                             wasy(1,is1,is2),ww,nu,n2)
              do l=1,n
                do k=1,n
                  do j=1,n
                    do i=1,n
                      ww(i,j,k,l)=ww(i,j,k,l)+v4(i,j,k,l,is1,is2)
                    enddo
                  enddo
                enddo
              enddo
              wggg=(0.d0,0.d0)
              do ind_nu=1,ndim3_nu
                i_nu=me_t*ndim3_nu+ind_nu-1
                call interp_g_spl(gexa(0,1,1,1,is1),
     &                            gspl_geom(0,1,1,1,1,is1),
     &                            gspl_asy(1,1,1,1,1,is1),
     &                            gasy(1,1,1,is1),g1,o1+w_nu(i_nu),
     &                            nn,1)
                call interp_g_spl(gexa(0,1,1,1,is1),
     &                            gspl_geom(0,1,1,1,1,is1),
     &                            gspl_asy(1,1,1,1,1,is1),
     &                            gasy(1,1,1,is1),g2,o1-nu+w_nu(i_nu),
     &                            nn,1)
	          call dgemm('n','n',2*n,nn*n,n,1.d0,g1,2*n,
     &                     w4(1,1,1,1,i_nu,is1,is1),n,0.d0,wg,
     &                     2*n)
	          call zgemm('n','n',nn*n,n,n,(1.d0,0.d0),wg,nn*n,
     &                     g2,n,(0.d0,0.d0),wgg,nn*n)
                if(i_nu/=0) then
                  call interp_g_spl(gexa(0,1,1,1,is1),
     &                              gspl_geom(0,1,1,1,1,is1),
     &                              gspl_asy(1,1,1,1,1,is1),
     &                              gasy(1,1,1,is1),g1,o1-w_nu(i_nu),
     &                              nn,1)
                  call interp_g_spl(gexa(0,1,1,1,is1),
     &                              gspl_geom(0,1,1,1,1,is1),
     &                              gspl_asy(1,1,1,1,1,is1),
     &                              gasy(1,1,1,is1),g2,o1-nu-w_nu(i_nu),
     &                              nn,1)
	            call dgemm('n','n',2*n,nn*n,n,1.d0,g1,2*n,
     &                       w4(1,1,1,1,i_nu,is1,is1),n,0.d0,wg,
     &                       2*n)
	            call zgemm('n','n',nn*n,n,n,(1.d0,0.d0),wg,nn*n,
     &                       g2,n,(1.d0,0.d0),wgg,nn*n)
                endif
                do l=1,n
                  do k=1,n
                    do j=1,n
                      do i=1,n
                        wg(j,k,l,i)=wgg(i,j,k,l)
                      enddo
                    enddo
                  enddo
                enddo
	          call dgemm('n','n',2*nn,nn,nn,1.d0,wg,2*nn,ww,nn,0.d0,
     &                     wgg,2*nn)
                do l=1,n
                  do k=1,n
                    do j=1,n
                      do i=1,n
                        wggg(i,k,j,l)=wggg(i,k,j,l)+trn(i_nu)
     &                                             *wgg(i,j,k,l)
                      enddo
                    enddo
                  enddo
                enddo
              enddo   !! over ind_nu
	        if(nproc_t/=1) call dgop(wggg,2*n**4,'  +',comm_t)
              do l=1,n
                do k=1,n
                  do j=1,n
                    do i=1,n
                      th0(i10+i,i20+j,i30+k,i40+l)=
     &                  th0(i10+i,i20+j,i30+k,i40+l)+wggg(i,j,k,l)
                    enddo
                  enddo
                enddo
              enddo
            enddo
          enddo
        enddo
      enddo
      deallocate(wgg,g1,g2,wg,wggg)
      deallocate(ww,vv,v1)
      end
      
      
      
      subroutine kernel_loc_5_2(n,n2,th0,o1,o2,nu,v4,w4,wexa,wasy,
     &                          wspl_geom,wspl_asy,gt,
     &                          trn,gexa,gspl_geom,gspl_asy,gasy)
      use atom_mod
      use etot_mod
      use manager_mod
      use parallel_mod
      use solid_mod
      use units_mod
      use vertex_mod
      implicit none
      integer, intent(in) :: n,n2
      real*8, intent(in) :: o1,o2,nu,v4(n,n,n,n,nspin,nspin),
     &                      wexa(0:n_nu_exa,n2,nspin,nspin),
     &                      wspl_geom(0:n_nu_geom+2,4,n2,nspin,nspin),
     &                      wspl_asy(n_nu_asy+1,4,n2,nspin,nspin),
     &                      wasy(n2,nspin,nspin),
     &                      gt(n,n,0:n_tau,nspin),trn(0:n_nu),
     &                      gexa(0:n_omega_exa,2,n,n,nspin),
     &                      gspl_geom(0:n_omega_geom+2,4,2,n,n,nspin),
     &                      gspl_asy(n_omega_asy+1,4,2,n,n,nspin),
     &                      gasy(2,n,n,nspin),
     &                      w4(n,n,n,n,0:n_nu,nspin,nspin)
      complex*16, intent(out) :: th0(2*n,2*n,2*n,2*n)
      integer :: isp1,isp2,isp3,isp4,is1,is2,is3,is4,i10,i20,i30,i40,nn,
     &           nud,it,i,j,k,l,ind_tau,i_tau,i1_tau,ind_nu,i_nu
      real*8 :: pi2b,a,b
      complex*16 :: cc
      real*8, allocatable :: ww(:,:,:,:),vv(:,:,:,:),gg(:,:,:,:,:,:)
      complex*16, allocatable :: wgg(:,:,:,:),g1(:,:),g2(:,:),
     &                           wg(:,:,:,:),wggg(:,:,:,:)
      pi2b=(pi+pi)/betta_t
      nn=n*n
      allocate(ww(n,n,n,n))
      allocate(vv(n,n,n,n))
      th0=(0.d0,0.d0)
c ------------- Static part ---------------------------------------
      allocate(gg(n,n,n,n,ndim3_tau,2))
      nud=idnint((o1+o2-nu)/pi2b)
      cc=(0.d0,1.d0)
      if(nud<0) cc=(0.d0,-1.d0)
      nud=iabs(nud)
      do isp4=1,2
        is4=min(isp4,nspin)
        do isp3=1,2
          is3=min(isp4,nspin)
          do isp2=1,2
            is2=min(isp4,nspin)
            do isp1=1,2
              is1=min(isp4,nspin)
              if(isp1/=isp3) cycle
              if(isp2/=isp4) cycle
              i40=n*(isp4-1)
              i30=n*(isp3-1)
              i20=n*(isp2-1)
              i10=n*(isp1-1)
              do l=1,n
                do k=1,n
                  do j=1,n
                    do i=1,n
                      vv(i,j,k,l)=v4(i,k,j,l,is1,is2)
                    enddo
                  enddo
                enddo
              enddo
              do ind_tau=1,ndim3_tau
                i_tau=me_t*ndim3_tau+ind_tau-1
                i1_tau=n_tau-i_tau
                do l=1,n
                  do k=1,n
                    do j=1,n
                      do i=1,n
                        a=gt(i,j,i_tau,is1)*gt(k,l,i_tau,is2)
                        b=gt(i,j,i1_tau,is1)*gt(k,l,i1_tau,is2)
                        gg(i,k,j,l,ind_tau,1)=a+b
                        gg(i,k,j,l,ind_tau,2)=a-b
                      enddo
                    enddo
                  enddo
                enddo
                do it=1,2
	            call dgemm('n','n',nn,nn,nn,1.d0,vv,nn,
     &                       gg(1,1,1,1,ind_tau,it),nn,0.d0,ww,nn)
	            call dgemm('n','t',nn,nn,nn,1.d0,ww,nn,vv,nn,0.d0,
     &                       gg(1,1,1,1,ind_tau,it),nn)
                enddo
              enddo   !! over ind_tau
              do l=1,n
                do k=1,n
                  do j=1,n
                    do i=1,n
                      a=dot_product(gg(i,l,k,j,:,1),
     &                              ab_nu_from_tau(:,nud,1))
                      b=dot_product(gg(i,l,k,j,:,2),
     &                              ab_nu_from_tau(:,nud,2))
                      th0(i10+i,i20+j,i30+k,i40+l)=
     &                  th0(i10+i,i20+j,i30+k,i40+l)+a+cc*b
                    enddo
                  enddo
                enddo
              enddo
            enddo
          enddo
        enddo
      enddo
      deallocate(gg)
	if(nproc_t/=1) call dgop(th0,32*n**4,'  +',comm_t)
c ------------- Semidynamic part ---------------------------------------
      allocate(wggg(n,n,n,n))
      allocate(wgg(n,n,n,n))
      allocate(wg(n,n,n,n))
      allocate(g1(n,n))
      allocate(g2(n,n))
      do isp4=1,2
        is4=min(isp4,nspin)
        do isp3=1,2
          is3=min(isp4,nspin)
          do isp2=1,2
            is2=min(isp4,nspin)
            do isp1=1,2
              is1=min(isp4,nspin)
              if(isp1/=isp3) cycle
              if(isp2/=isp4) cycle
              i40=n*(isp4-1)
              i30=n*(isp3-1)
              i20=n*(isp2-1)
              i10=n*(isp1-1)
              do l=1,n
                do k=1,n
                  do j=1,n
                    do i=1,n
                      vv(i,j,k,l)=v4(i,k,j,l,is1,is2)
                    enddo
                  enddo
                enddo
              enddo
              wggg=(0.d0,0.d0)
              do ind_nu=1,ndim3_nu
                i_nu=me_t*ndim3_nu+ind_nu-1
                call interp_g_spl(gexa(0,1,1,1,is1),
     &                            gspl_geom(0,1,1,1,1,is1),
     &                            gspl_asy(1,1,1,1,1,is1),
     &                            gasy(1,1,1,is1),g1,o2+w_nu(i_nu),
     &                            nn,1)
                  call interp_g_spl(gexa(0,1,1,1,is2),
     &                              gspl_geom(0,1,1,1,1,is2),
     &                              gspl_asy(1,1,1,1,1,is2),
     &                              gasy(1,1,1,is2),g2,o1-nu-w_nu(i_nu),
     &                              nn,1)
	          call dgemm('n','n',2*n,nn*n,n,1.d0,g1,2*n,
     &                     w4(1,1,1,1,i_nu,is1,is2),n,0.d0,wg,2*n)
	          call zgemm('n','n',nn*n,n,n,(1.d0,0.d0),wg,nn*n,
     &                     g2,n,(0.d0,0.d0),wgg,nn*n)
                if(i_nu/=0) then
                  call interp_g_spl(gexa(0,1,1,1,is1),
     &                              gspl_geom(0,1,1,1,1,is1),
     &                              gspl_asy(1,1,1,1,1,is1),
     &                              gasy(1,1,1,is1),g1,o2-w_nu(i_nu),
     &                              nn,1)
                  call interp_g_spl(gexa(0,1,1,1,is2),
     &                              gspl_geom(0,1,1,1,1,is2),
     &                              gspl_asy(1,1,1,1,1,is2),
     &                              gasy(1,1,1,is2),g2,o1-nu+w_nu(i_nu),
     &                              nn,1)
	            call dgemm('n','n',2*n,nn*n,n,1.d0,g1,2*n,
     &                       w4(1,1,1,1,i_nu,is1,is2),n,0.d0,wg,2*n)
	            call zgemm('n','n',nn*n,n,n,(1.d0,0.d0),wg,nn*n,
     &                       g2,n,(1.d0,0.d0),wgg,nn*n)
                endif
                do l=1,n
                  do k=1,n
                    do j=1,n
                      do i=1,n
                        wggg(k,j,i,l)=wggg(k,j,i,l)+trn(i_nu)
     &                                             *wgg(i,j,k,l)
                      enddo
                    enddo
                  enddo
                enddo
              enddo   !! over ind_nu
	        if(nproc_t/=1) call dgop(wggg,2*n**4,'  +',comm_t)
	        call dgemm('n','t',2*nn,nn,nn,1.d0,wggg,2*nn,
     &                   vv,nn,0.d0,wg,2*nn)
              do l=1,n
                do k=1,n
                  do j=1,n
                    do i=1,n
                      th0(i10+i,i20+j,i30+k,i40+l)=
     &                  th0(i10+i,i20+j,i30+k,i40+l)+wg(j,k,i,l)
                    enddo
                  enddo
                enddo
              enddo
              wggg=(0.d0,0.d0)
              do ind_nu=1,ndim3_nu
                i_nu=me_t*ndim3_nu+ind_nu-1
                call interp_g_spl(gexa(0,1,1,1,is1),
     &                            gspl_geom(0,1,1,1,1,is1),
     &                            gspl_asy(1,1,1,1,1,is1),
     &                            gasy(1,1,1,is1),g1,o1-w_nu(i_nu),
     &                            nn,1)
                call interp_g_spl(gexa(0,1,1,1,is2),
     &                            gspl_geom(0,1,1,1,1,is2),
     &                            gspl_asy(1,1,1,1,1,is2),
     &                            gasy(1,1,1,is2),g2,o2-nu+w_nu(i_nu),
     &                            nn,1)
	          call dgemm('n','n',2*n,nn*n,n,1.d0,g1,2*n,
     &                     w4(1,1,1,1,i_nu,is1,is2),n,0.d0,wg,2*n)
	          call zgemm('n','n',nn*n,n,n,(1.d0,0.d0),wg,nn*n,
     &                     g2,n,(0.d0,0.d0),wgg,nn*n)
                if(i_nu/=0) then
                  call interp_g_spl(gexa(0,1,1,1,is1),
     &                              gspl_geom(0,1,1,1,1,is1),
     &                              gspl_asy(1,1,1,1,1,is1),
     &                              gasy(1,1,1,is1),g1,o1+w_nu(i_nu),
     &                              nn,1)
                  call interp_g_spl(gexa(0,1,1,1,is2),
     &                              gspl_geom(0,1,1,1,1,is2),
     &                              gspl_asy(1,1,1,1,1,is2),
     &                              gasy(1,1,1,is2),g2,o2-nu-w_nu(i_nu),
     &                              nn,1)
	            call dgemm('n','n',2*n,nn*n,n,1.d0,g1,2*n,
     &                       w4(1,1,1,1,i_nu,is2,is1),n,0.d0,wg,2*n)
	            call zgemm('n','n',nn*n,n,n,(1.d0,0.d0),wg,nn*n,
     &                       g2,n,(1.d0,0.d0),wgg,nn*n)
                endif
                do l=1,n
                  do k=1,n
                    do j=1,n
                      do i=1,n
                        wggg(j,k,i,l)=wggg(j,k,i,l)+trn(i_nu)
     &                                             *wgg(i,j,k,l)
                      enddo
                    enddo
                  enddo
                enddo
              enddo   !! over ind_nu
	        if(nproc_t/=1) call dgop(wggg,2*n**4,'  +',comm_t)
	        call dgemm('n','t',2*nn,nn,nn,1.d0,wggg,2*nn,
     &                   vv,nn,0.d0,wg,2*nn)
              do l=1,n
                do k=1,n
                  do j=1,n
                    do i=1,n
                      th0(i10+i,i20+j,i30+k,i40+l)=
     &                  th0(i10+i,i20+j,i30+k,i40+l)+wg(i,l,k,j)
                    enddo
                  enddo
                enddo
              enddo           
            enddo
          enddo
        enddo
      enddo
      deallocate(wgg,g1,g2,wg,wggg)
c ------------- Dynamic part ---------------------------------------
      allocate(wggg(n,n,n,n))
      allocate(wgg(n,n,n,n))
      allocate(wg(n,n,n,n))
      allocate(g1(n,n))
      allocate(g2(n,n))
      do isp4=1,2
        is4=min(isp4,nspin)
        do isp3=1,2
          is3=min(isp4,nspin)
          do isp2=1,2
            is2=min(isp4,nspin)
            do isp1=1,2
              is1=min(isp4,nspin)
              if(isp1/=isp3) cycle
              if(isp2/=isp4) cycle
              i40=n*(isp4-1)
              i30=n*(isp3-1)
              i20=n*(isp2-1)
              i10=n*(isp1-1)
              wggg=(0.d0,0.d0)
              do ind_nu=1,ndim3_nu
                i_nu=me_t*ndim3_nu+ind_nu-1
                call interp_g_spl(gexa(0,1,1,1,is1),
     &                            gspl_geom(0,1,1,1,1,is1),
     &                            gspl_asy(1,1,1,1,1,is1),
     &                            gasy(1,1,1,is1),g1,o2+w_nu(i_nu),
     &                            nn,1)
                call interp_g_spl(gexa(0,1,1,1,is2),
     &                            gspl_geom(0,1,1,1,1,is2),
     &                            gspl_asy(1,1,1,1,1,is2),
     &                            gasy(1,1,1,is2),g2,o1-nu-w_nu(i_nu),
     &                            nn,1)
                call interp_w_nu_spl(wexa(0,1,is1,is2),
     &                               wspl_geom(0,1,1,is1,is2),
     &                               wspl_asy(1,1,1,is1,is2),
     &                               wasy(1,is1,is2),ww,
     &                               o2-o1+w_nu(i_nu),n2)
	          call dgemm('n','n',2*n,nn*n,n,1.d0,g1,2*n,ww,n,0.d0,wg,
     &                     2*n)
	          call zgemm('n','n',nn*n,n,n,(1.d0,0.d0),wg,nn*n,
     &                     g2,n,(0.d0,0.d0),wgg,nn*n)
                if(i_nu/=0) then
                  call interp_g_spl(gexa(0,1,1,1,is1),
     &                              gspl_geom(0,1,1,1,1,is1),
     &                              gspl_asy(1,1,1,1,1,is1),
     &                              gasy(1,1,1,is1),g1,o2-w_nu(i_nu),
     &                              nn,1)
                  call interp_g_spl(gexa(0,1,1,1,is2),
     &                              gspl_geom(0,1,1,1,1,is2),
     &                              gspl_asy(1,1,1,1,1,is2),
     &                              gasy(1,1,1,is2),g2,o1-nu+w_nu(i_nu),
     &                              nn,1)
                  call interp_w_nu_spl(wexa(0,1,is1,is2),
     &                                 wspl_geom(0,1,1,is1,is2),
     &                                 wspl_asy(1,1,1,is1,is2),
     &                                 wasy(1,is1,is2),ww,
     &                                 o2-o1-w_nu(i_nu),n2)
	            call dgemm('n','n',2*n,nn*n,n,1.d0,g1,2*n,ww,n,0.d0,
     &                       wg,2*n)
	            call zgemm('n','n',nn*n,n,n,(1.d0,0.d0),wg,nn*n,
     &                       g2,n,(1.d0,0.d0),wgg,nn*n)
                endif
                do l=1,n
                  do k=1,n
                    do j=1,n
                      do i=1,n
                        wg(j,k,i,l)=wgg(i,j,k,l)
                        ww(j,l,k,i)=w4(i,j,k,l,i_nu,is1,is2)
                      enddo
                    enddo
                  enddo
                enddo
	          call dgemm('n','n',2*nn,nn,nn,1.d0,wg,2*nn,ww,nn,0.d0,
     &                     wgg,2*nn)
                do l=1,n
                  do k=1,n
                    do j=1,n
                      do i=1,n
                        wggg(i,k,l,j)=wggg(i,k,l,j)+trn(i_nu)
     &                                             *wgg(i,j,k,l)
                      enddo
                    enddo
                  enddo
                enddo
              enddo   !! over ind_nu
	        if(nproc_t/=1) call dgop(wggg,2*n**4,'  +',comm_t)
              do l=1,n
                do k=1,n
                  do j=1,n
                    do i=1,n
                      th0(i10+i,i20+j,i30+k,i40+l)=
     &                  th0(i10+i,i20+j,i30+k,i40+l)+wggg(i,j,k,l)
                    enddo
                  enddo
                enddo
              enddo
            enddo
          enddo
        enddo
      enddo
      deallocate(wgg,g1,g2,wg,wggg)
      deallocate(ww,vv)
      end