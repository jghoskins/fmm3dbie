        implicit real *8 (a-h,o-z)
        complex *16 ima
        real *8, allocatable :: srcvals(:,:),srccoefs(:,:)
        integer, allocatable :: norders(:),ixyzs(:),
     1      iptype(:)
        real *8, allocatable :: uvs_pts(:,:)
        integer, allocatable :: ipatch_id(:)
        complex *16, allocatable :: srcvals_c(:,:),srccoefs_c(:,:)
        complex *16, allocatable :: cms_c(:,:)
        real *8, allocatable :: rads(:),rads_near(:)
        integer, allocatable :: row_ptr(:),col_ind(:),
     1      iquad(:)
        complex *16, allocatable :: wnear_c(:)
        integer, allocatable :: irowind(:),icolind(:)
        data ima /(0,1)/

        call prini(6,13)
        print *, 'Enter n0'
        read *, n0

        r    = 1
        dl   = 10
        rmax = 1
        call get_cyl_mem(r,dl,rmax,npatches)

        norder = 10
        npol   = (norder+1)*(norder+2)/2
        npts   = npatches*npol
        allocate(srcvals(12,npts),srccoefs(9,npts),
     1      norders(npatches),ixyzs(npatches+1),
     1      iptype(npatches))

        call get_cyl_geom(r,dl,rmax,norder,npatches,
     1      npts,norders,ixyzs,iptype,srcvals,srccoefs)

        iw = 21
        call arrwrite(iw,srcvals,12*npts)


        deallocate(norders,ixyzs)
        deallocate(iptype)

cc
cc          .   .   .   now complex

        r    = 1
        dl   = 10
        rmax = 1
        call get_cyl_mem(r,dl,rmax,npatches)

        norder = 5
        npol   = (norder+1)*(norder+2)/2
        npts   = npatches*npol
        allocate(srcvals_c(12,npts),srccoefs_c(9,npts),
     1      norders(npatches),ixyzs(npatches+1),
     1      iptype(npatches))

        call get_cyl_geom_c(r,dl,rmax,norder,npatches,
     1      npts,norders,ixyzs,iptype,srcvals_c,srccoefs_c)

        iw = 22
        call arrwrite(iw,srcvals_c,2*12*npts)


cc
cc          .   .   .   prepare for quad generation

        allocate(uvs_pts(2,npts),ipatch_id(npts))

cc          .   .   .   in surface_routs
cc          get the patch id for each point in the mesh

        call get_patch_id_uvs(npatches,norders,ixyzs,
     1      iptype,npts,ipatch_id,uvs_pts)

cc          .   .   .   in near_field_routs
cc          get the default radii for computing near field quads 

        call get_rfacs(norder,iptype,rfac,rfac0)


cc          .   .   .   in surface_routs
cc          get the centers (cmplx) and radii (real) for each patch

        allocate(cms_c(3,npatches))
        allocate(rads(npatches))
        
        ifre = 1
        call get_centroid_rads_c(npatches,norders,ixyzs,iptype,
     1      npts,srccoefs_c,cms_c,rads,ifre)

        call prin2('rads = *',rads,18)
        call prin2('cms_c = *',cms_c,18)

        allocate(rads_near(npatches))

        do i=1,npatches
        rads_near(i) = rads(i)*rfac
        enddo

cc          .   .   .   in near field routs        
cc          process geometry to get nearby targets per patch
        ifre = 0

        ndt = 12
        call findnearmem_c(cms_c,npatches,rads_near,ndt,
     1        srcvals_c,npts,nnz,ifre)


        call prinf('nnz = *',nnz,1)
        call prinf('npts= *',npts,1)

        allocate(row_ptr(npts+1))
        allocate(col_ind(nnz))

        ndt = 12
        call findnear_c(cms_c,npatches,rads_near,ndt,srcvals_c,
     1      npts,row_ptr,col_ind,ifre) 

        allocate(iquad(nnz+1))
        call get_iquad_rsc(npatches,ixyzs,npts,nnz,row_ptr,
     1    col_ind,iquad)

        nquad = iquad(nnz+1)-1
        call prinf('nquad = *',nquad,1)
        call prin2('quad ratio = *',nquad*1.0d0/(npts**2),1)

        allocate(wnear_c(nquad))
        allocate(irowind(nquad),icolind(nquad))


cc          .   .   .   helm_comb_dir

        iquadtype = 1
        eps = 1.0d-4
        ndtarg = 12
        ntarg  = npts

        call getnearquad_helm_comb_dir_c(npatches,norders,
     1   ixyzs,iptype,npts,srccoefs_c,srcvals_c,ndtarg,ntarg,
     1   srcvals_c,
     2   ipatch_id,uvs_targ,eps,zpars,iquadtype,nnz,row_ptr,col_ind,
     3   iquad,rfac0,nquad,wnear_c,ifre)

        call prinf('completed *',0,0)

        stop
        end
c
c
c
c
c
        subroutine arrwrite(iw,arr,n)
        implicit real *8 (a-h,o-z)
        dimension arr(n)
c
 1200 format(6x,d38.32)
        write(iw,1200) (arr(i),i=1,n)
c
        return
        end
c
c
c
c
c
        subroutine get_cyl_mem(r,dl,rmax,npatches)
        implicit real *8 (a-h,o-z)
        real *8, intent(in) :: r,dl,rmax
        integer, intent(out) :: npatches

        done = 1
        pi   = atan(done)*4

        nthet = ceiling(2*pi*r/rmax)
        nleng = ceiling(2*dl/rmax)

        npatches = nthet*nleng*2

        return
        end
c
c
c
c
c
        subroutine get_cyl_geom(r,dl,rmax,norder,npatches,
     1      npts,norders,ixyzs,iptype,srcvals,srccoefs)
        implicit real *8 (a-h,o-z)
        real *8, intent(in) :: r,dl,rmax
        integer, intent(in) :: npatches,norder,npts
        integer, intent(out) :: norders(npatches),ixyzs(npatches+1)
        integer, intent(out) :: iptype(npatches)
        real *8, intent(out) :: srcvals(12,npts),srccoefs(9,npts)
        real *8, pointer :: ptr1,ptr2,ptr3,ptr4
        real *8, allocatable :: triaskel(:,:,:)
        real *8, allocatable :: uvs(:,:),wts(:),umatr(:,:),vmatr(:,:)
        external xtri_cyl_eval

        done = 1
        pi   = atan(done)*4

        npols = (norder+1)*(norder+2)/2
        do i=1,npatches
            norders(i) = norder
            ixyzs(i)   = (i-1)*npols+1
            iptype(i)  = 1
        enddo
        ixyzs(npatches + 1) = npts + 1

        allocate(triaskel(3,3,npatches))

        nthet = ceiling(2*pi*r/rmax)
        nleng = ceiling(2*dl/rmax)
        nover = 0
        n0    = 0

        umin =  0
        umax =  2*pi*r
        vmin = -dl
        vmax =  dl

        call xtri_rectmesh_ani(umin,umax,vmin,vmax,nthet,nleng,
     1      nover,npatches,npatches,triaskel)

        npols = (norder+1)*(norder+2)/2
        allocate(uvs(2,npols),wts(npols),umatr(npols,npols),
     1      vmatr(npols,npols))
        call vioreanu_simplex_quad(norder,npols,uvs,umatr,vmatr,
     1      wts)
        
        call getgeominfo(npatches,xtri_cyl_eval,triaskel,
     1      r,pjunk1,pjunk2,npols,uvs,umatr,srcvals,srccoefs)

        call prin2('srcvals = *',srcvals,20)

        return
        end
c
c
c
c
c
        subroutine get_cyl_geom_c(r,dl,rmax,norder,npatches,
     1      npts,norders,ixyzs,iptype,srcvals,srccoefs)
        implicit real *8 (a-h,o-z)
        real *8, intent(in) :: r,dl,rmax
        integer, intent(in) :: npatches,norder,npts
        integer, intent(out) :: norders(npatches),ixyzs(npatches+1)
        integer, intent(out) :: iptype(npatches)
        complex *16, intent(out) :: srcvals(12,npts),srccoefs(9,npts)
        real *8, allocatable :: triaskel(:,:,:)
        real *8, allocatable :: uvs(:,:),wts(:),umatr(:,:),vmatr(:,:)
        external xtri_cyl_eval_null_c

        done = 1
        pi   = atan(done)*4

        npols = (norder+1)*(norder+2)/2
        do i=1,npatches
            norders(i) = norder
            ixyzs(i)   = (i-1)*npols+1
            iptype(i)  = 1
        enddo
        ixyzs(npatches + 1) = npts + 1

        allocate(triaskel(3,3,npatches))

        nthet = ceiling(2*pi*r/rmax)
        nleng = ceiling(2*dl/rmax)
        nover = 0
        n0    = 0

        umin =  0
        umax =  2*pi*r
        vmin = -dl
        vmax =  dl

        call xtri_rectmesh_ani(umin,umax,vmin,vmax,nthet,nleng,
     1      nover,npatches,npatches,triaskel)

        npols = (norder+1)*(norder+2)/2
        allocate(uvs(2,npols),wts(npols),umatr(npols,npols),
     1      vmatr(npols,npols))
        call vioreanu_simplex_quad(norder,npols,uvs,umatr,vmatr,
     1      wts)
        
        call getgeominfo_c(npatches,xtri_cyl_eval_null_c,triaskel,
     1      r,pjunk1,pjunk2,npols,uvs,umatr,srcvals,srccoefs)

        call prin2('srcvals = *',srcvals,20)

        return
        end
c
c
c
c
c
        subroutine xtri_cyl_eval(itri,u,v,xyz,dxyzduv,
     1      triainfo,r,p1,p2)
        implicit real *8 (a-h,o-z)
        real *8 :: xyz(3), dxyzduv(3,2),triainfo(3,3,*)

        x0 = triainfo(1,1,itri)
        y0 = triainfo(2,1,itri)
        z0 = triainfo(3,1,itri)

        x1 = triainfo(1,2,itri)
        y1 = triainfo(2,2,itri)
        z1 = triainfo(3,2,itri)

        x2 = triainfo(1,3,itri)
        y2 = triainfo(2,3,itri)
        z2 = triainfo(3,3,itri)

        x = x0 + u*(x1-x0) + v*(x2-x0)
        y = y0 + u*(y1-y0) + v*(y2-y0)
        z = z0 + u*(z1-z0) + v*(z2-z0)

        dxdu = x1-x0
        dydu = y1-y0
        dzdu = z1-z0

        dxdv = x2-x0
        dydv = y2-y0
        dzdv = z2-z0

        xyz(1) = r*cos(x/r)
        xyz(2) = r*sin(x/r)
        xyz(3) = y

        dxyzduv(1,1) = -sin(x/r)*dxdu + 0*dydu
        dxyzduv(2,1) =  cos(x/r)*dxdu + 0*dydu
        dxyzduv(3,1) =  0*dxdu + 1*dydu

        dxyzduv(1,2) = -sin(x/r)*dxdv
        dxyzduv(2,2) =  cos(x/r)*dxdv
        dxyzduv(3,2) =  dydv

        return
        end
c
c
c
c
c
        subroutine xtri_cyl_eval_null_c(itri,u,v,xyz,dxyzduv,
     1      triainfo,r,p1,p2)
        implicit real *8 (a-h,o-z)
        complex *16 :: xyz(3), dxyzduv(3,2)
        real *8 :: triainfo(3,3,*)

        x0 = triainfo(1,1,itri)
        y0 = triainfo(2,1,itri)
        z0 = triainfo(3,1,itri)

        x1 = triainfo(1,2,itri)
        y1 = triainfo(2,2,itri)
        z1 = triainfo(3,2,itri)

        x2 = triainfo(1,3,itri)
        y2 = triainfo(2,3,itri)
        z2 = triainfo(3,3,itri)

        x = x0 + u*(x1-x0) + v*(x2-x0)
        y = y0 + u*(y1-y0) + v*(y2-y0)
        z = z0 + u*(z1-z0) + v*(z2-z0)

        dxdu = x1-x0
        dydu = y1-y0
        dzdu = z1-z0

        dxdv = x2-x0
        dydv = y2-y0
        dzdv = z2-z0

        xyz(1) = r*cos(x/r)
        xyz(2) = r*sin(x/r)
        xyz(3) = y

        dxyzduv(1,1) = -sin(x/r)*dxdu + 0*dydu
        dxyzduv(2,1) =  cos(x/r)*dxdu + 0*dydu
        dxyzduv(3,1) =  0*dxdu + 1*dydu

        dxyzduv(1,2) = -sin(x/r)*dxdv
        dxyzduv(2,2) =  cos(x/r)*dxdv
        dxyzduv(3,2) =  dydv

        return
        end
