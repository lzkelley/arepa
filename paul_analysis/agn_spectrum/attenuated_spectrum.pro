function attenuated_spectrum, log_nu_input, log_l_bol, log_NH
	nu_min = 9.
	nu_max = 21.
	d_log_nu = 0.1
		log_nu = nu_min + d_log_nu*findgen((nu_max-nu_min)/d_log_nu+1.)
		nu = 10^(DOUBLE(log_nu))
		cx = 0.0*log_nu
	for i=0, n_elements(nu)-1 do cx[i] = cross_section(nu[i])
	tau = cx * 10^(DOUBLE(log_NH))
	atten = EXP(-DOUBLE(tau))
	bad = where((FINITE(atten) NE 1) OR (FINITE(atten,/NAN) EQ 1), n_bad)
	if (n_bad GT 0) then atten[bad] = 1.0d-40
	
	l_bol = 10^(DOUBLE(log_l_bol))
	un_atten_spectrum = 10^(DOUBLE(intrinsic_spectrum(log_nu,log_l_bol)))
	;; renormalize it as a safety check for heavily obscured systems
	un_atten_spectrum = un_atten_spectrum * $
		(l_bol / TOTAL(DOUBLE(un_atten_spectrum)*alog(10.)*d_log_nu))	
	atten_spectrum = un_atten_spectrum * atten


	;; mock BB spectrum at temp T (arbitrary normalization radius R0)
	;; model as a smooth disk with T ~ r^-1/2
		pc = 3.086d18
	r_inner = 0.85  * pc
	r_outer = 300. * pc
	T_inner = 2000.
	
	d_log_r = 0.01
	r = 10^(alog10(r_inner) + d_log_r*findgen((alog10(r_outer/r_inner)/d_log_r+1.)))
		kB= 1.38d-16
		h = 6.626d-27
		nu = 10^(DOUBLE(log_nu))
		c = 3.0d10
	L = 0.0*nu
	for i = 0, n_elements(r)-1 do begin
		T_r = T_inner * ((r[i]/r_inner)^(-0.55))
		dA_r= 2.*!PI*alog(10.)*r[i]*r[i]*d_log_r
		Fnu = !Pi * (2.*h*(nu*nu*nu)/(c*c)) / (EXP(h*nu/(kB*T_r))-1.)
		dL = nu * (dA_r * Fnu) / (3.9d33)
		L = L + dL
	endfor
	L_IR_bump = L 	
	;; now that provides a very good model of the observed IR-bump portion 
	;;	of the spectrum -- rescale it, not some arbitrary cutoff

	;; now cut above ~1.5 microns & add all the attenuated 
	;; energy back in the IR bump
	nu_c = 3.0d8/(1.5d-6)
	l_bol_obs = TOTAL(DOUBLE(atten_spectrum)*alog(10.)*d_log_nu)
	l_bol_abs = l_bol - l_bol_obs
		if (l_bol_abs GE l_bol) then l_bol_abs = l_bol

	l_tot_ir_c = TOTAL(DOUBLE(L_IR_bump)*alog(10.)*d_log_nu)
	f_ir = L_IR_bump * (l_bol_abs/l_tot_ir_c)
	
	;;
	;; its actually much more consistent with the radiative 
	;;   transfer models and the observed SED models of canonical 
	;;   type-II qsos to assume that the absorbed energy is 
	;;   re-radiated in a single, thermal cold dust bump (T~50K, 
	;;   with a rough range ~30-70 K for the most part)
	;;
	h = 6.625d-27
	k = 1.381d-16
	c = 2.998d10
	T = 50.0	
	x = (h*(10^(log_nu)))/(k*T)	
	nu_L_nu = (15./(!PI^(4.))) * (x^(4.) / (exp(x)-1.))  
	f_ir = l_bol_abs * nu_L_nu
	;;
	;; note that even with high obscuration, this 
	;;  only changes the IR sed owing to re-radiation at 
	;;  wavelengths >~ 30 microns, so our 8-15 micron 
	;;  calculations (the hot dust bump) are basically unaffected
	;;
	
	f = atten_spectrum + f_ir

if (2 eq 0) then begin
	plot,log_nu,alog10(f),xstyle=1,xrange=[8.5,19.2],ystyle=1, $
		yrange=[log_l_bol-5.,log_l_bol],THICK=2.0
	oplot,log_nu,alog10(atten_spectrum),color=80
	oplot,log_nu,alog10(f_ir),color=250
	oplot,[-10.,200.],[14.6,14.6]
	oplot,[-10.,200.],[13.6,13.6]
	oplot,[-10.,200.],[12.6,12.6]
endif
	
	return, INTERPOL(alog10(f),log_nu,log_nu_input)
end
