pro test
fname_set=['9','5','6','7','8','10','11']
color_set=[30,80,150,210,250,100,110]
mdot_set=[0.025,0.01,0.005,0.0025,0.001,0.01,0.01]

;; can show in plot as range: w. intrinsic -2.+/-0.2, giving HR vs. NH+z; 
;;   to compare w. model SEDs

plot,findgen(2),/nodata,xstyle=1,ystyle=1,xrange=[10.,22.],yrange=[38.,44.]

l_sun=3.9d33
for i=0,n_elements(fname_set)-1 do begin

  fname='obslum.dat'+fname_set[i]
  readcol,fname,log_nu,log_nuLnu,/silent
	
	;;;log_nuLnu=log_nuLnu+log_nu-18.
	
   nu=10^(double(log_nu)) & cx=0.*nu
   for i_nu=0,n_elements(nu)-1 do cx(i_nu)=cross_section(nu(i_nu)) 
    nh=1.0d14 & tau=cx*nh 
    log_nuLnu=log_nuLnu-tau*alog10(exp(1.))

   ;;;log_nu=log_nu-alog10(1.+3.)

   oplot,log_nu,log_nuLnu,thick=2.,linestyle=0,color=color_set(i)
   
   
   d_log_nu=log_nu(1)-log_nu(0)
   d_ln_nu=d_log_nu*alog(10.)
   d_L=10^(double(log_nuLnu)) * d_ln_nu(0)
   L=total(d_L)

   L_0=1.0d8 * 3.3d4 * l_sun
   print, (L/L_0)
   
   rad_corr=mdot_set(i)
    rad_corr_break=0.1
    if (rad_corr lt rad_corr_break) then rad_corr=rad_corr*(rad_corr/rad_corr_break)   
   L_max=rad_corr*L_0
   
   nulnu_disk=agn_spectrum(10^log_nu,alog10(L_max/l_sun))+alog10(l_sun)


    nh=1.0d22 & tau=cx*nh 
   nulnu_disk=nulnu_disk - 0.2*(log_nu-18.)

   nulnu_disk=nulnu_disk-tau*alog10(exp(1.))
   
   oplot,log_nu,nulnu_disk,thick=2.,linestyle=2,color=color_set(i)


   nu=10^(double(log_nu)) & nuLnu=10^(double(log_nuLnu)) & Lnu=nuLnu/nu 
   h=6.6261d-27 & E_nu=h*nu & Nnu=Lnu/E_nu & dnu=nu*d_ln_nu(0) & dN=Nnu*dnu
   eV=1.6022d-12 & keV=eV*1.0d3

   z=1.

   emin=0.5*keV*(1.+z) & emax=2.0*keV*(1.+z)
	ok=(E_nu ge emin)*(E_nu le emax) & N=total(dN*ok)
	N_soft=N
   emin=2.0*keV*(1.+z) & emax=10.0*keV*(1.+z)
	ok=(E_nu ge emin)*(E_nu le emax) & N=total(dN*ok)
	N_hard=N
   
   HR=(N_hard-N_soft)/(N_hard+N_soft)
   print, HR



   nu=10^(double(log_nu)) & nuLnu=10^(double(nulnu_disk)) & Lnu=nuLnu/nu 
   h=6.6261d-27 & E_nu=h*nu & Nnu=Lnu/E_nu & dnu=nu*d_ln_nu(0) & dN=Nnu*dnu
   eV=1.6022d-12 & keV=eV*1.0d3

   emin=0.5*keV*(1.+z) & emax=2.0*keV*(1.+z)
	ok=(E_nu ge emin)*(E_nu le emax) & N=total(dN*ok)
	N_soft=N
   emin=2.0*keV*(1.+z) & emax=10.0*keV*(1.+z)
	ok=(E_nu ge emin)*(E_nu le emax) & N=total(dN*ok)
	N_hard=N
   
   HR=(N_hard-N_soft)/(N_hard+N_soft)
   print, HR
   


   nu=10^(double(log_nu)) & nuLnu=10^(double(0.*nulnu_disk+42.-tau*alog10(exp(1.)))) & Lnu=nuLnu/nu 
	oplot,alog10(nu),alog10(nuLnu),color=0
   h=6.6261d-27 & E_nu=h*nu & Nnu=Lnu/E_nu & dnu=nu*d_ln_nu(0) & dN=Nnu*dnu
   eV=1.6022d-12 & keV=eV*1.0d3

   emin=0.5*keV*(1.+z) & emax=2.0*keV*(1.+z)
	ok=(E_nu ge emin)*(E_nu le emax) & N=total(dN*ok)
	N_soft=N
   emin=2.0*keV*(1.+z) & emax=10.0*keV*(1.+z)
	ok=(E_nu ge emin)*(E_nu le emax) & N=total(dN*ok)
	N_hard=N
   
   HR=(N_hard-N_soft)/(N_hard+N_soft)
   print, HR
   
   
   gamma=1.5
    alpha=1.-gamma + 1.
    oplot,log_nu,alpha*(log_nu-17.)+42.
   
   gamma=2.0
    alpha=1.-gamma + 1.
    oplot,log_nu,alpha*(log_nu-17.)+42.

   gamma=1.0
    alpha=1.-gamma + 1.
    ;;oplot,log_nu,alpha*(log_nu-17.)+42.
   
   
   
endfor
end



pro vs_z
 
 setup_plot_ps,'hardness_vs_z_template.ps'
 
 nh_grid=double( [21.5,22.,22.5] )
  color_nh=[ 30, 150, 250 ]
 gamma_grid=double( [1.8,2.0,2.2] )
  lsty_gam=[  2,   0,   3 ]



 log_nu=double(6.+0.01*findgen((24.-6.)/0.01+1.))
   nu=10^(double(log_nu)) & cx=0.*nu
   for i_nu=0,n_elements(nu)-1 do cx(i_nu)=double(cross_section(nu(i_nu)))
   l_sun=3.9d33 & h=6.6261d-27 & E_nu=h*nu & eV=1.6022d-12 & keV=eV*1.0d3 
   d_log_nu=log_nu(1)-log_nu(0) & d_ln_nu=d_log_nu*alog(10.) & dnu=nu*d_ln_nu(0)
   dN_conv=dnu * 1./(nu*E_nu)
   nulnu_disk=agn_spectrum(nu,double(14.))+alog10(l_sun)
   gamma_assumed=double(1.7)

 z_grid=double( 0.+0.1*findgen(4./0.1+1.) )
 HR_z=fltarr(n_elements(z_grid),n_elements(gamma_grid),n_elements(nh_grid))
 for i_z=0,n_elements(z_grid)-1 do begin
	z=z_grid(i_z)
	soft=(E_nu ge 0.5*keV*(1.+z))*(E_nu le 2.0*keV*(1.+z))
	hard=(E_nu ge 2.0*keV*(1.+z))*(E_nu le 10.0*keV*(1.+z))
	
	for i_gamma=0,n_elements(gamma_grid)-1 do begin
		gamma=gamma_grid(i_gamma)
		nuLnu=nuLnu_disk - (gamma - gamma_assumed)*(log_nu-18.)
		for i_nh=0,n_elements(nh_grid)-1 do begin
			tau = (10^nh_grid(i_nh)) * cx
			nuLnu_obs = nuLnu - tau*alog10(exp(1.))
			dN = (10^nuLnu_obs) * dN_conv			
			 N_soft=total(dN*soft)
			 N_hard=total(dN*hard)
			 HR = (N_hard-N_soft)/(N_hard+N_soft)
			 HR_z[i_z,i_gamma,i_nh] = HR(0)
		endfor
	endfor
 endfor
 
 
 plot,findgen(2),/nodata,xstyle=1,ystyle=1,xrange=[0.,4.],yrange=[-1.,1.]
  x=-1.+0.1*findgen(81) & oplot,x,0.*x-0.3,thick=4.,linestyle=1,color=0
 
  for i_gamma=0,n_elements(gamma_grid)-1 do begin
  for i_nh=0,n_elements(nh_grid)-1 do begin
    oplot,z_grid,HR_z[*,i_gamma,i_nh],thick=4.,color=color_nh[i_nh],linestyle=lsty_gam[i_gamma]  
  endfor
  endfor
  
end


