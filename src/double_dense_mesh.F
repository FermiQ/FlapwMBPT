      subroutine double_dense_mesh(key,i0,i1,nn,n,key1,mesh)
c     key = 1 Mesh is denser from the left only 
c     key = 2 Mesh is denser from the right only
c     key = 3 Mesh is denser from both sides
c     i0 and i1 should be given as bondaries
c     n is the number of points obtained
c     We work with INTEGER numbers here
c     key1 =0 - just to get NN 
c     key1 =1 - to fill in MESH array         
      implicit none
	integer, intent(in) :: key,i0,i1,key1,nn
      integer, intent(out) :: mesh(0:nn),n
	integer :: ih,n1,im,jj
	integer, allocatable :: mesh0(:)
	if(key1==1) allocate(mesh0(0:nn))
	if(key==1) then
	  ih=1
	  n=0
	  im=i0
	  if(key1==1) mesh(n)=im
	  do while(im<i1)
	    im=im+ih
	    if(im>i1) then
	      if(key1==1) mesh(n)=i1
	      exit
	    endif
	    ih=ih+ih
	    n=n+1
	    if(key1==1) mesh(n)=im
	  enddo
	else if(key==2) then
	  ih=1
	  n=0
	  im=i1
	  if(key1==1) mesh0(nn)=im
	  do while(im>i0)
	    im=im-ih
	    if(im<i0) then
	      if(key1==1) mesh0(nn-n)=i0
	      exit
	    endif
	    ih=ih+ih
	    n=n+1
	    if(key1==1) mesh0(nn-n)=im
	  enddo
	  if(key1==1) mesh(0:n)=mesh0(nn-n:nn)
	else if(key==3) then
	  jj=(i0+i1)/2
	  ih=1
	  n=0
	  im=i0
	  if(key1==1) mesh(n)=im
	  do while(im<jj)
	    im=im+ih
	    if(im>jj) then
	      if(key1==1) mesh(n)=jj
	      exit
	    endif
	    ih=ih+ih
	    n=n+1
	    if(key1==1) mesh(n)=im
	  enddo
	  ih=1
	  n1=0
	  im=i1
	  if(key1==1) mesh0(nn)=im
	  do while(im>jj)
	    im=im-ih
	    if(im<jj) then
	      if(key1==1) mesh0(nn-n1)=jj
	      exit
	    endif
	    ih=ih+ih
	    n1=n1+1
	    if(key1==1) mesh0(nn-n1)=im
	  enddo
	  if(key1==1) mesh(n:n+n1)=mesh0(nn-n1:nn)
	  n=n+n1
	endif
	if(key1==1) deallocate(mesh0)
      end