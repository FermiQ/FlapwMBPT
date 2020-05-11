      subroutine mesh_omega_for_omega(ist,ifin,n_ex,n_as,nn,mesh_mats,
     &                                key)
c     key = -1 asymptotic is from the left
c     key =  1 asymptotic is from the right
c ------- MESH for Omega ------------------------------------------------
	use manager_mod
	use models_mod
	use parallel_mod
	use units_mod
	use vertex_mod
      implicit none
      integer, intent(in) :: ist,ifin,key
      integer, intent(inout) :: nn
      integer, intent(out) :: n_ex,n_as,mesh_mats(0:n_omega)
	integer :: nh,ii,i,iii,n
	real*8 :: pi_bet,h1
	integer, allocatable :: mesh(:)
	if(ist>ifin) then
	  n_ex=-1
	  n_as=0
	  return
	endif
	pi_bet=pi/betta_t
	h1=pi_bet+pi_bet
	n=ifin-ist
	allocate(mesh(0:n+1))
	n_ex=min(n_omega_exa,n)
	do i=0,n_ex
	  mesh(i)=i
	enddo
	n_as=0
	if(n>mesh(n_ex)) then
        n_as=n_omega_asy
c ---------- Find asymptotic part of the mesh ----------------------
        nh=2
        do i=n_ex+1,n_ex+n_omega_asy
          mesh(i)=mesh(i-1)+nh
          if(n-mesh(i)<=nh) then
            n_as=i-n_ex
            mesh(i)=n
            exit
          endif
          nh=2*nh
        enddo
        if(mod(n_as,2)==1) then
          n_ex=n_ex-1
          n_as=n_as+1
        endif
      endif
c ------------------------------------------------------------------      
	if(key==-1) iii=n_ex+n_as
	if(key==1) iii=0
	do ii=0,n_ex+n_as   !! over intervals
	  if(key==-1) mesh_mats(ii)=ifin-mesh(iii)
	  if(key==1) mesh_mats(ii)=ist+mesh(iii)
	  iii=iii+key
	enddo
	nn=nn+n_ex+n_as+1
      deallocate(mesh)
      end