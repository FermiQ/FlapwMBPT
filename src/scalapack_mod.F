      module scalapack_mod
      logical :: lsclp
      integer :: nblock_pw,
     &           ncol_b,
     &           ncol_pw_max,
     &           nrow_b,
     &           nrow_pw_max,
     &           sctxt_pw
      integer, allocatable :: desc_pw(:,:),
     &                        map_pbr(:,:),
     &                        ncol_pw(:),
     &                        nrow_pw(:),
     &                        refcol_pw(:,:),
     &                        refrow_pw(:,:)
      end
