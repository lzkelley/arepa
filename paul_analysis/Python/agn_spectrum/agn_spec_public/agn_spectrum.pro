;;
;; function to return the AGN spectrum
;;   given nu (in Hz), and log(L_bol/L_sun), returns 
;;   log(nu*L_nu/L_sun)  :: NOTE: L_sun here *always* has the same units, 
;;                           namely *bolometric* solar, i.e. 3.9x10^33 erg/s
;;                           (given this, it should be trivial to modify the code 
;;                            if you prefer to work directly in CGS units)
;;
;;	keywords to shortcut for various bands:
;;	  BB = B-band (4400 Angstroms), 
;;    IR = mid-IR (15 microns), 
;;    SX = soft X-rays (0.5-2 keV), 
;;    HX = hard X-rays (2-10 keV)
;;
;;
;;  keywords for different spectra:
;;
;;    HRH = Hopkins, Richards, & Hernquist 2006 luminosity-dependent template spectrum
;;	    compiled from a number of observations therein. If no keyword is set, 
;;      this is the default spectrum used.
;;
;;    RICHARDS = Richards et al. 2006 all-quasar mean SED: if you want a 
;;      non-luminosity-dependent spectrum for whatever reason, this is essentially 
;;      an update of Elvis et al. 1994 with an expanded and more homogenously-selected 
;;      quasar sample. Note this does not include an X-ray reflection component, 
;;      instead simply assuming a flat X-ray spectrum in nu*L_nu
;;
;;    SDSS is an optional keyword that can be added as a keyword with *either* 
;;      of the model spectra. It will overlay the vanden Berk et al. 2001 
;;      median SDSS SED over the relevant 
;;      portion of the spectrum :: basically does so with a sliding continuum determination
;;      that then overlays this spectrum, such that the integrated luminosity in the 
;;      entire bandpass and the continuum spectrum are conserved (i.e. you can put this 
;;      over an arbitrary continuum)
;;
;;
;;	  These model spectra are compiled in the papers above, but draw heavily from the 
;;    observations in a number of papers, including :
;;      Richards et al. 2006, Steffen et al. 2006, Hatziminaoglou et al. 2005, 
;;      Tozzi et al. 2006, Strateva et al. 2005, Telfer et al. 2002, Vanden Berk et al. 2001, 
;;      George et al. 1998, Elvis et al. 1994
;;		also, Marconi et al. 2004 (we follow a similar methodology to their derivation)
;;		and note that we generate the X-ray reflection component 
;;        following Ueda et al. 2003 in the PEXRAV code, Magdziarz & Zdziarski 1995
;;
;;
;;	Example use :: 
;;     nu_L_nu = agn_spectrum( 4.5e14, 12.0 )
;;		  returns nu*L_nu at nu=4.5x10^14 Hz for a quasar with L_bol=10^12 L_sun
;;		  since no keyword were set, the HRH06 spectrum is returned, without
;;			overlaying the SDSS template spectrum in the optical
;;	   nu_L_nu = agn_spectrum( 10^[12., 13., 14., 15., 16., 17., 18.], 12.0 )
;;		  returns nu*L_nu at nu=10^12, 10^13, 10^14, 10^15, 10^16, 10^17, 10^18 Hz, 
;;          for a quasar with L_bol=10^12 L_sun (i.e. nu can be a vector)
;;	   nu_L_nu = agn_spectrum( 4.5e14, [12.0, 13.0, 14.0, 15.0, 16.0] )
;;		  returns nu*L_nu at nu=4.5x10^14 Hz, 
;;          for quasars with L_bol=10^12, 10^13, 10^14, 10^15, 10^16 L_sun 
;;          (i.e. L can also be a vector)
;;	   nu_L_nu = agn_spectrum( 10^[13., 14., 15.], [ 13.0, 14.0 ] )
;;		  if both nu and L are vectors, it will return every combination, in an 
;;          array with [N_ELEMENTS(nu), N_ELEMENTS(L)], where the first set of 
;;          elements give the nu index, the second set the L index, for each 
;;          corresponding value of nu*L_nu at some L
;;	   nu_L_nu = agn_spectrum( 10^[13., 14., 15.], [ 13.0, 14.0 ] , /BB)
;;	   nu_L_nu = agn_spectrum( 0., [ 13.0, 14.0 ] , /BB)
;;		  regardless of what nu is set to, setting the keyword BB will override it 
;;          and return the result *just* for one value of nu, namely 
;;          that of B-band. Likewise for the similar band keywords
;;          BB, IR, SX, and HX
;;	   nu_L_nu = agn_spectrum( 10^[13., 14., 15.], [ 13.0, 14.0 ] , /RICHARDS)
;;		  acts just as usual, but returns nu*L_nu based on the template 
;;          spectrum in Richards et al. 2006, instead of that in HRH06. 
;;          This can be used in combination with any other keyword(s), except
;;          of course HRH, which explicitly chooses the HRH06 spectrum (although 
;;          that is already the default)
;;	   nu_L_nu = agn_spectrum( 10^[13., 14., 15.], [ 13.0, 14.0 ] , /SDSS)
;;	 	  acts just as usual, but returns nu*L_nu by taking the model spectrum 
;;          (either HRH06 or Richards et al. 2006, if the /RICHARDS keyword is set), 
;;          and overlaying the SDSS vanden Berk et al. 2001 spectrum on the optical 
;;          portion of the template spectrum. This can be combined with any other 
;;          keyword(s).
;;
;;
function agn_spectrum, nu_in_Hz, log_l_bol, $
	BB=BB, IR=IR, SX=SX, HX=HX, $
	HRH=HRH, RICHARDS=RICHARDS, $
	SDSS=SDSS

	exec_call='./agn_spectrum_shared.so'   
			;; location where the code "agn_spectrum.c"
			;; has been compiled as a shared library

	nu_in_Hz     = DOUBLE(nu_in_Hz)
		if (keyword_set(BB)) then nu_in_Hz = -1.0d0
		if (keyword_set(IR)) then nu_in_Hz = -2.0d0
		if (keyword_set(SX)) then nu_in_Hz = -3.0d0
		if (keyword_set(HX)) then nu_in_Hz = -4.0d0	
		N_nu  = LONG(n_elements(nu_in_Hz))
		
	;; just puts everything in the right units for the code
	;; recall, the c-code takes nu=-1.0 to be the shortcut for the B-band, 
	;;   -2.0 = mid-IR (15 micron), -3.0 = soft X-ray (0.5-2 keV), 
	;;   -4.0 = hard X-ray (2-10 keV)

	log_l_bol    = DOUBLE(log_l_bol)
		N_lum = LONG(n_elements(log_l_bol))

	spectrum_key = 0L
		if (keyword_set(HRH)) 		then spectrum_key = 0L
		if (keyword_set(RICHARDS)) 	then spectrum_key = 1L	
	sloan_key    = 0L
		if (keyword_set(SDSS))		then sloan_key = 1L

	l_band_all = DOUBLE(fltarr(N_nu,N_lum))

	;;;;
	;;;; the code now loops over the (possibly scalars, possibly vector lists) 
	;;;; values of nu and log_L given to it, and calls the c-code for each 
	;;;; pair, to determine nu*L_nu in the desired band
	;;;;
	for i=0, N_lum-1 do begin
	l_band_vec = DOUBLE(fltarr(N_nu))
	log_l_bol_pass = DOUBLE(log_l_bol[i])
	S = CALL_EXTERNAL(exec_call, $
           'main',N_nu,nu_in_Hz,log_l_bol_pass,spectrum_key,sloan_key,l_band_vec)
	l_band_all[*,i] = l_band_vec
	endfor


	;; return the nu*L_nu corresponding to each nu given to this routine
	return, l_band_all
end

