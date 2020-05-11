      subroutine k_mesh_coarse
c	------- Coarse Mesh in Brillouin Zone ----------------------------
	  use manager_mod
	  use parallel_mod
      use sitegen_mod
	  use solid_mod
	  use vertex_mod
      implicit none
      integer :: k,i1,i2,i3,k0,n1,n2,n3,m
      nqdiv_c=ndiv_c(1)*ndiv_c(2)*ndiv_c(3)
      allocate(k_a_from_c(nqdiv_c))
      allocate(k_c_from_a(nqdiv))
      allocate(pnt_c(3,nqdiv_c))
      k_c_from_a=0
      n1=ndiv(1)/ndiv_c(1)
      n2=ndiv(2)/ndiv_c(2)
      n3=ndiv(3)/ndiv_c(3)
      k=0
      do k0=1,nqdiv
        i1=idnint(dot_product(pnt(:,k0),rb0(:,1)))
        if(mod(i1,n1)==0) then
          i2=idnint(dot_product(pnt(:,k0),rb0(:,2)))
          if(mod(i2,n2)==0) then
            i3=idnint(dot_product(pnt(:,k0),rb0(:,3)))
            if(mod(i3,n3)==0) then
              k=k+1
              pnt_c(:,k)=pnt(:,k0)
              k_a_from_c(k)=k0
              k_c_from_a(k0)=k
              if(k0<=npnt) npnt_c=k
            endif
          endif
        endif
      enddo
	  do m = 1, 3
        rb0_c(m,1) = rbas(m,1)*ndiv_c(1)
        rb0_c(m,2) = rbas(m,2)*ndiv_c(2)
        rb0_c(m,3) = rbas(m,3)*ndiv_c(3)
        qb0_c(m,1) = gbas(m,1)/ndiv_c(1)
        qb0_c(m,2) = gbas(m,2)/ndiv_c(2)
        qb0_c(m,3) = gbas(m,3)/ndiv_c(3)
      enddo
      allocate(i_kref_c(nqdiv_c))
	  allocate(k_group_c(nqdiv_c))
	  allocate(index_k1_c(nqdiv_c))
	  allocate(ip_k_c(nqdiv_c,ngroup))
	  allocate(k_star_c(npnt_c))
	  allocate(wgt_c(npnt_c))
	  call inf_qnt_c
	  do k=1,npnt_c
	    wgt_c(k)=dfloat(k_star_c(k))/nqdiv_c
	  enddo
      end

