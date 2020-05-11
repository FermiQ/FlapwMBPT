      subroutine get_fi0_red
      use atom_mod
      use manager_mod
      use solid_mod
      use units_mod
      use vertex_mod
      implicit none
      integer :: ispin,isort,le1,mt1,mt,le,ir,i,mtpb,jspin,le0,le10,ie,
     &           ie1,li1,l1,li,l
      real*8 :: c2,dqdall
      real*8, allocatable :: pp(:),pp1(:),qq(:),qq1(:),work(:)
      allocate(pp(0:maxnrad),pp1(0:maxnrad),qq(0:maxnrad),
     &	     qq1(0:maxnrad),work(0:maxnrad))
      c2=clight**2
      do jspin=1,nspin
        do ispin=1,nspin
          do isort=1,nsort
            do le10=1,lfun_red(isort)
              le1=le_red(le10,isort)
              li1=li_val(le1,isort)
              if(irel.le.1) l1=li1-1
              if(irel.eq.2) l1=li1/2
              ie1=ie_val(le1,isort)
              if(correlated(ie1,l1,isort)=='A') then
                do ir=0,nrad(isort)
                  pp1(ir)=p_f(ir,ie1,li1,isort,jspin)
                  if(irel.ge.1) qq1(ir)=q_f(ir,ie1,li1,isort,jspin)
                enddo
              else if(correlated(ie1,l1,isort)=='L') then 
                mt1=ind_wf(le1,isort)
                do ir=0,nrad(isort)
                  pp1(ir)=gfun(mt1+ir,jspin)*r(ir,isort)
                  if(irel>=1) qq1(ir)=gfund(mt1+ir,jspin)*r(ir,isort)
                enddo
              endif
              do le0=1,lfun_red(isort)
                le=le_red(le0,isort)
                li=li_val(le,isort)
                if(irel.le.1) l=li-1
                if(irel.eq.2) l=li/2
                ie=ie_val(le,isort)
                if(correlated(ie,l,isort)=='A') then
                  do ir=0,nrad(isort)
                    pp(ir)=p_f(ir,ie,li,isort,ispin)
                    if(irel.ge.1) qq(ir)=q_f(ir,ie,li,isort,ispin)
                  enddo
                else if(correlated(ie1,l1,isort)=='L') then
                  mt=ind_wf(le,isort)
                  do ir=0,nrad(isort)
                    pp(ir)=gfun(mt+ir,ispin)*r(ir,isort)
                    if(irel>=1) qq(ir)=gfund(mt+ir,ispin)*r(ir,isort)
                  enddo
                endif
                do i=1,lfun_pb_red(isort)
                  mtpb=ind_pb_red(i,isort)
                  do ir=0,nrad(isort)
                    work(ir)=pp(ir)*pp1(ir)
                    if(irel>=1) work(ir)=work(ir)+qq(ir)*qq1(ir)/c2
                    work(ir)=work(ir)*fun_pb_red(ir+mtpb)*dr(ir,isort)
                  enddo
                  fi0_red(le0,le10,i,isort,ispin,jspin)=
     &              dqdall(h(isort),work,nrad(isort))
                enddo   !! over i
              enddo   !! over le0
            enddo   !! over le10
          enddo  !! over isort
        enddo  !! over ispin
      enddo  !! over jspin
      deallocate(pp,pp1,qq,qq1,work)
      end
      
      
      
      subroutine get_fi0_f_red
      use atom_mod
      use manager_mod
      use solid_mod
      use units_mod
      use vertex_mod
      implicit none
      integer :: ispin,isort,le1,mt1,mt,le,ir,i,mtpb,jspin,ie,
     &           ie1,li1,l1,li,l
      real*8 :: c2,dqdall
      real*8, allocatable :: pp(:),pp1(:),qq(:),qq1(:),work(:)
      allocate(pp(0:maxnrad),pp1(0:maxnrad),qq(0:maxnrad),
     &	     qq1(0:maxnrad),work(0:maxnrad))
      c2=clight**2
      do jspin=1,nspin
        do ispin=1,nspin
          do isort=1,nsort
            do le1=1,lfun(isort)
              li1=li_val(le1,isort)
              if(irel.le.1) l1=li1-1
              if(irel.eq.2) l1=li1/2
              ie1=ie_val(le1,isort)
              do ir=0,nrad(isort)
                pp1(ir)=p_f(ir,ie1,li1,isort,jspin)
                if(irel.ge.1) qq1(ir)=q_f(ir,ie1,li1,isort,jspin)
              enddo
              mt1=ind_wf(le1,isort)
              do ir=0,nrad(isort)
                pp1(ir)=gfun(mt1+ir,jspin)*r(ir,isort)
                if(irel>=1) qq1(ir)=gfund(mt1+ir,jspin)*r(ir,isort)
              enddo
              do le=1,lfun(isort)
                li=li_val(le,isort)
                if(irel.le.1) l=li-1
                if(irel.eq.2) l=li/2
                ie=ie_val(le,isort)
                do ir=0,nrad(isort)
                  pp(ir)=p_f(ir,ie,li,isort,ispin)
                  if(irel.ge.1) qq(ir)=q_f(ir,ie,li,isort,ispin)
                enddo
                mt=ind_wf(le,isort)
                do ir=0,nrad(isort)
                  pp(ir)=gfun(mt+ir,ispin)*r(ir,isort)
                  if(irel>=1) qq(ir)=gfund(mt+ir,ispin)*r(ir,isort)
                enddo
                do i=1,lfun_pb_red(isort)
                  mtpb=ind_pb_red(i,isort)
                  do ir=0,nrad(isort)
                    work(ir)=pp(ir)*pp1(ir)
                    if(irel>=1) work(ir)=work(ir)+qq(ir)*qq1(ir)/c2
                    work(ir)=work(ir)*fun_pb_red(ir+mtpb)*dr(ir,isort)
                  enddo
                  fi0_f_red(le,le1,i,isort,ispin,jspin)=
     &              dqdall(h(isort),work,nrad(isort))
                enddo   !! over i
              enddo   !! over le
            enddo   !! over le1
          enddo  !! over isort
        enddo  !! over ispin
      enddo  !! over jspin
      deallocate(pp,pp1,qq,qq1,work)
      end