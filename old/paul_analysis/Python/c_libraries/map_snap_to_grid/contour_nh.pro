pro contour_nh, frun, snapnum, xlen, sendto, $
                        center=center, $
						use_calc_center=use_calc_center, $
                        colorbar=colorbar, $
                        crude=crude, $
                        filename=filename, $
                        fitstoo=fitstoo, $
                        loadedsnap=loadedsnap, $
                        msg=msg, $
						nodust=nodust, $
                        nolabels=nolabels, $
                        particlesonly=particlesonly, $
                        pixels=pixels, $
                        pubstyle=pubstyle, $
                        plotpts=plotpts, $
                        rotate_phi=rotate_phi, rotate_theta=rotate_theta, $
                        startid=startid, numpart=numpart, $
                        thumbnail=thumbnail, $
                        track_to_draw= track_to_draw, $
                        xthickness=xthickness, ythickness=ythickness, zthickness=zthickness, $
						showbhs=showbhs, $
                        xz=xz, yz=yz

if not keyword_set(xlen) then xlen=100.0
if not keyword_set(sendto) then sendto='x'
if not keyword_set(snapnum) then snapnum= 0

setup_plot_ps,'idl.ps'


; --------------------------------
;  Center this on something other
;  than 0,0,0
; --------------------------------
if not keyword_set(center) then center=[0.0,0.0,0.0]


set_maxden=10.0
set_dynrng=1.0e+5

; --------------------------------
;  Load variables for smoothing
; --------------------------------
    if not keyword_set(loadedsnap) then begin
      ;if (fload_snapshot(frun, snapnum)) then begin
      if (fload_snapshot_bh(frun, snapnum, nopot=1)) then begin
        print, "PROBLEM: opening file"
        return
      endif
    endif

    if keyword_set(use_calc_center) then center=fload_center_alreadycomp(1)

    Ngas= fload_npart(0)
    
    ;; setup the background grid of fake stars: 
	N_grid_side = 30L
		Ngrid = N_grid_side*N_grid_side
		d_x_grid = 2.*xlen/float(N_grid_side)
		x_grid_l = -xlen + d_x_grid*findgen(N_grid_side)
		x_grid = fltarr(Ngrid)
		y_grid = fltarr(Ngrid)
		print, center[2]
		for i = 0, N_grid_side-1 do begin
			ni = i*N_grid_side
			x_grid[0+ni:0+ni+N_grid_side-1] = 0.0*x_grid_l + center[0] + x_grid_l[i]
			y_grid[0+ni:0+ni+N_grid_side-1] = center[1] + x_grid_l
		endfor
		;; plot,x_grid,y_grid,PSYM=7
		;; and put it in the background
		z_offset_default = -200.	
			;; make sure to keep the sign negative, that's what puts them in the background!
			z_offset_min = center[2]-xlen
			if (z_offset_min LT z_offset_default) then z_offset_default = z_offset_min
 		z_grid = 0.0*x_grid + z_offset_default
	
		;; set the smoothing length to be comparable to the grid scale
		hsml = 0.0*x_grid + 2.0*d_x_grid


	; need it in radians
        if keyword_set(rotate_theta) then theta = (!PI / 180.0) * rotate_theta else theta= 0.0
        if keyword_set(rotate_phi) then phi = (!PI / 180.0) * rotate_phi else phi= 0.0



	
	proper_rotation = 0	;; if want to scale the background sources to the nhat vector
	if (proper_rotation EQ 1) then begin

	;; (need to add a rotation matrix for the background plane)
	
	;; really not so bad, b/c contour code choses a fixed axis to project onto; 
	;;   so only need to switch up for different choice of axes
	xi = x_grid
	yi = y_grid
	zi = z_grid
	if (keyword_set(xz)) then begin
		zi = y_grid + center[2] - center[1]
		yi = z_grid
	endif
	if (keyword_set(yz)) then begin
		zi = x_grid + center[2] - center[0]
		xi = z_grid
	endif
	x_grid = xi
	y_grid = yi
	z_grid = zi
	
	
		if (keyword_set(rotate_theta) OR keyword_set(rotate_phi)) then begin
		;; rezero to the center of the plane
		x = x_grid - center[0]
		y = y_grid - center[1]
		z = 0.*z_grid
		r_i = SQRT(x*x + y*y)
			nx = SIN(theta)*COS(phi)
			ny = SIN(theta)*SIN(phi)
			nz = COS(theta)

		phi_0 = atan(y/x)
			div = where((abs(y) GE abs(x)) AND ((y/x) GT 0.), n_div)
				if (n_div GT 0) then phi_0[div] = !PI/2.
			div = where((abs(y) GE abs(x)) AND ((y/x) LT 0.), n_div)
				if (n_div GT 0) then phi_0[div] = -!PI/2.
		xf = r_i * COS(phi + phi_0) * COS(theta)
		yf = r_i * SIN(phi + phi_0) * COS(theta)
		s0 = - (x*nx + y*ny + z*nz)
			s00 = 1.*(s0 GT 0.) - 1.*(s0 LT 0.)
		zf = r_i * SIN(theta) * s00

		;; redo to make sure signs aren't messed up
		xf = (x*COS(phi) - y*SIN(phi)) * COS(theta)
		yf = (x*SIN(phi) + y*COS(phi)) * COS(theta)
		zf =-(x*COS(phi) + y*SIN(phi)) * SIN(theta)

		r00 = z_offset_default
		x0  = center[0] + r00*nx
		y0  = center[1] + r00*ny
		z0  = center[2] + r00*nz
		
		x_grid = xf + x0
		y_grid = yf + y0
		z_grid = zf + z0
	endif	
	endif
	

        ; fields to pass
        Coord = fltarr(13,Ngas+Ngrid)

        ; gas
        x_gas = fload_gas_xyz('x', center=[0,0,0])
		y_gas = fload_gas_xyz('y', center=[0,0,0])
		z_gas = fload_gas_xyz('z', center=[0,0,0])
			;; rezero to the center of the plane
			x = x_gas - center[0]
			y = y_gas - center[1]
			z = z_gas - center[2]

			n_hat_zz_x = SIN(theta) * COS(phi)
			n_hat_zz_y = SIN(theta) * SIN(phi)
			n_hat_zz_z = COS(theta)
			n_hat_xx_x = SIN(theta+!PI/2.) * COS(phi)
			n_hat_xx_y = SIN(theta+!PI/2.) * SIN(phi)
			n_hat_xx_z = COS(theta+!PI/2.)
			n_hat_yy_x = SIN(theta+!PI/2.) * COS(phi+!PI/2.)
			n_hat_yy_y = SIN(theta+!PI/2.) * SIN(phi+!PI/2.)
			n_hat_yy_z = COS(theta+!PI/2.)

			x_proj = x*n_hat_xx_x + y*n_hat_xx_y + z*n_hat_xx_z
			y_proj = x*n_hat_yy_x + y*n_hat_yy_y + z*n_hat_yy_z
			z_proj = x*n_hat_zz_x + y*n_hat_zz_y + z*n_hat_zz_z

			x_gas  = x_proj + center[0] 
			y_gas  = y_proj + center[1] 
			z_gas  = z_proj + center[2] 
	;; since theta, phi set by rotating the galaxy above, they should not be set here
		theta = 0.
		rotate_theta = 0.
		phi   = 0.
		rotate_phi	 = 0.
	
        Coord(0,0:Ngas-1) = float(x_gas)
        Coord(1,0:Ngas-1) = float(y_gas)
        Coord(2,0:Ngas-1) = float(z_gas)

        Coord(3,0:Ngas-1) = fload_gas_u(1)
        Coord(4,0:Ngas-1) = fload_gas_rho(1)
        Coord(5,0:Ngas-1) = fload_gas_hsml(1)
        Coord(6,0:Ngas-1) = fload_gas_numh(1)
        Coord(7,0:Ngas-1) = fload_gas_nume(1)
        Coord(8,0:Ngas-1) = fload_gas_metallicity(1)
        Coord(9,0:Ngas-1) = fload_gas_mass(1)

        Coord(10,0:Ngas-1) = fload_gas_v('x')
        Coord(11,0:Ngas-1) = fload_gas_v('y')
        Coord(12,0:Ngas-1) = fload_gas_v('z')

        ; background grid points
        Coord(0,Ngas:Ngrid+Ngas-1) = x_grid
        Coord(1,Ngas:Ngrid+Ngas-1) = y_grid
        Coord(2,Ngas:Ngrid+Ngas-1) = z_grid
              
        Nstars = long(Ngrid)
        Nbh = 0L

        los_NH= fltarr(Ngrid)
        los_Z= fltarr(Ngrid)

        print, "PASSING: "
        print, "N_gas= ", Ngas
        print, "N_grid= ", Nstars
        print, "theta= ", theta
        print, "phi= ", phi
        help, Coord

	N_VEL_GRID = 200L
	N_LOS_NH_VEL_GRID = Ngrid * N_VEL_GRID	;; allow for each to be a grid in velocity
	los_NH = fltarr(N_LOS_NH_VEL_GRID)
	los_Z  = fltarr(N_LOS_NH_VEL_GRID)

	nh_code_exe = return_idl_c_dir(0)+'/LOS_column_background/getnh.so'
	;;nh_code_exe = return_idl_c_dir(0)+'/LOS_column_code/getnh.so'
	S = CALL_EXTERNAL(nh_code_exe, $
                'getnh', $
                Ngas, $ 
                Nstars, $
                Nbh, $ 
                theta, $
                phi, $ 
                Coord, $
                los_NH, $
				los_Z)

	; trap for really low NH values
	;   and zero metallicity (make it small instead)
	idx= where(los_NH lt 1.0e+10)
	if idx(0) ne -1 then los_NH(idx)= 1.0e+10
	idx= where(los_Z le 0.0)
	if idx(0) ne -1 then los_Z(idx)= 1.0e-5
	
	;; convert the single-vector los_NH into a matrix -- velocity x background location
	;;   (i.e. items 0-200 are velocity bins for first of N_grid, etc.
	;;
	los_NH_v = fltarr(Ngrid,N_VEL_GRID)
	los_Z_v  = fltarr(Ngrid,N_VEL_GRID)
	for i=0,Ngrid-1 do begin
		los_NH_v[i,*] = los_NH[i*N_VEL_GRID:i*N_VEL_GRID+N_VEL_GRID-1]
		los_Z_v[i,*]  =  los_Z[i*N_VEL_GRID:i*N_VEL_GRID+N_VEL_GRID-1]
	endfor
	
	
vgrid = -5000.0 + 50.0*findgen(200)
plot,vgrid,0.*vgrid,xstyle=1,xrange=[MIN(vgrid),MAX(vgrid)],ystyle=1,yrange=[17.,22.], $
	xtitle=textoidl("v_{los} [km s^{-1}]"),ytitle=textoidl("log( d N_{H} / dln v_{los}  )  [cm^-2]")
for i=0,(Ngrid-1)/10 do begin
	oplot,vgrid,alog10(los_NH_v[10*i,*]),color=250.*i/((Ngrid-1)/10)
endfor
;; oplot simplest expectation for iso-sphere mass distribution and v=r/t
v0 = 300.0
	oplot,vgrid,alog10(5.0d20/(1.0 + (vgrid/v0)^2)),linestyle=0
v0 = 500.0
	oplot,vgrid,alog10(3.0d20/(1.0 + (vgrid/v0)^2)),linestyle=0


;; crude local maxima selection
v_all_maxs  = [0.]
NH_all_maxs = [0.]
for i=0,Ngrid-1 do begin
  x=vgrid
  y=alog10(los_NH_v[i,*])
  ny=n_elements(y)
  y0=y[0:ny-3]
  y1=y[1:ny-2]
  y2=y[2:ny-1]
  ymaxs=where((y1 GT y0) AND (y1 GT y2), n_maxs)
  if (n_maxs GT 0) then begin
    ym=y1[ymaxs]
    xm=x[1:ny-2]
    xm=xm[ymaxs]
    v_all_maxs  = [v_all_maxs,  xm]
    NH_all_maxs = [NH_all_maxs, ym]
  endif
endfor
v_all_maxs=v_all_maxs[1:n_elements(v_all_maxs)-1]
NH_all_maxs=NH_all_maxs[1:n_elements(NH_all_maxs)-1]
v_all_maxs=ABS(v_all_maxs)
plot,findgen(2),/nodata,xstyle=1,xrange=[1.,4.d3],ystyle=1,yrange=[17.5,21.]
  plotsym,0,/fill,0.1
  oplot,v_all_maxs,NH_all_maxs,psym=8,color=90

BIN=100.0
h=HISTOGRAM(v_all_maxs,MIN=0.,MAX=3.d3,BINS=BIN,LOCATIONS=x)
  plot,x+BIN/2.,h,psym=10
BIN=0.1
h=HISTOGRAM(NH_all_maxs,MIN=17.5,MAX=21.,BINS=BIN,LOCATIONS=x)
  plot,x+BIN/2.,h,psym=10


;; covering factors ::

;; NH, averaged over all
f=(findgen(n_elements(NH_all_maxs))+1.0)/float(Ngrid)
plot,findgen(2),xrange=[17.5,21.5],xstyle=1,/nodata,yrange=[0.01,4.],ystyle=1,/ylog
	oplot,NH_all_maxs[REVERSE(sort(NH_all_maxs))],f,PSYM=10 
		;; not quite -- this includes multiple components per line, not the same 
		;;   exactly as probability of seeing *something*
	

;; velocity, averaged over all
f=(findgen(n_elements(NH_all_maxs))+1.0)/float(Ngrid)
plot,findgen(2),xrange=[1.,4.d3],xstyle=1,yrange=[0.01,4.],ystyle=1,/ylog
	oplot,v_all_maxs[REVERSE(sort(v_all_maxs))],f,PSYM=10

	
if (2 EQ 0) then begin
for i=0,9 do begin
	;; project the column density outputs
	ProjectThisNH = los_NH_v[*,20*i]
		print, ' column density range = ',MIN(alog10(los_NH)),' to ',MAX(alog10(los_NH))


; ----------------------------------------------
;  Call our generic contour plotting procedure
; ----------------------------------------------

LOSMAX = 1.0d24
LOSMAX = MAX(los_NH)
set_maxden=LOSMAX * ((1.0/d_x_grid)^2)
set_dynrng=LOSMAX/1.0d19
set_dynrng=LOSMAX/MIN(los_NH)


contour_makeplot, x_grid, y_grid, z_grid, ProjectThisNH, hsml, xlen, sendto, xz=xz, yz=yz, $
					filename=filename, thumbnail=thumbnail, fitstoo=fitstoo, $
					xthickness=xthickness, ythickness=ythickness, colorbar=colorbar, $
					pixels=pixels, zthickness=zthickness, $
					crude=crude, center=center, msg=msg, $
					plotpts=plotpts, particlesonly=particlesonly, $
					rotate_phi=rotate_phi, rotate_theta=rotate_theta, $
					nolabels=nolabels, pubstyle=pubstyle, $
					set_maxden=set_maxden, set_dynrng=set_dynrng, $
					track_to_draw=track_to_draw, showbhs=showbhs


endfor
endif

; -------------
;  Done
; -------------




device,/close
end


