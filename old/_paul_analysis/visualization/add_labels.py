import visualization.movie_maker as mm
import numpy as np
import os

# simple function to (optionally) add redshift labels to set of dumped 
#   frames for movies -- much of this needs to be customized for the 
#   particular set of frames, so play around
#
def add_frame_labels( prefix='m12_mr_s', suffix='_t090_gas_N3c_r200',
      fileext='jpeg', frames_per_gyr=100., a_min=0.0322581, a_max=1.0, 
      myfont='Palatino-Bold', mysize='45', mycolor='white', mylocation='text 20,46 ' ):
    
    time_min = mm.cosmological_time(a_min) ## input in 'a_scale', not time (Gyr)
    time_max = mm.cosmological_time(a_max)
    time_frame_grid = np.arange(time_min, time_max, 1./np.float(frames_per_gyr))
    a_tmp = 10.**np.arange(-3.,0.001,0.001)
    t_tmp = mm.cosmological_time(a_tmp)
    a_scale_grid = np.exp(np.interp(np.log(time_frame_grid),np.log(t_tmp),np.log(a_tmp)))
    a_scale_grid[a_scale_grid < a_tmp[0]] = a_tmp[0]
    z_scale_grid=1./(1.+a_scale_grid)

    ## ok these are the times we used, each should be associated with the corresponding frame
    #for i in range(10):
    for i in range(a_scale_grid.size):
        lbl = get_time_label(a_scale_grid[i],cosmo=1) # ok this is the string
        snum = snap_ext(i,four_char=1) # string for frame number

        fname=prefix+snum+suffix
        infi=fname+'.'+fileext
        outfi=fname+'_lbl.'+fileext
        
        #myfont = 'AGaramondPro-Bold' # palatino, georgia, AGaramondPro-Bold
        #myfont = 'Georgia-Bold'
        #myfont = 'Palatino-Bold'
        #myfont = 'Cambria-Bold'
        #mysize = '45'
        #mycolor= 'white'
        #mylocation = 'text 20,46 '
        cmd = 'convert -font '+myfont+' -antialias -pointsize '+mysize+' -fill '+mycolor+\
          ' -draw '+"'"+mylocation+'"'+lbl+'"'+"'"+' -quality 95 '+infi+' '+outfi
        print cmd
        os.system(cmd)
        #!{cmd} # only works within ipython directly
        #convert -font  -antialias -pointsize 45 -fill white -draw 'text 20,46 "z=10.5"' -quality 95 $infi $outfi



def snap_ext(snum,four_char=0):
	ext='00'+str(snum);
	if (snum>=10): ext='0'+str(snum)
	if (snum>=100): ext=str(snum)
	if (four_char==1): ext='0'+ext
	if (snum>=1000): ext=str(snum)
	return ext;    

def get_time_label(time, cosmo=0, c0='w', n_sig_add=0, tunit_suffix='Gyr'):
    if(cosmo==0): 
        time_to_use = time ## absolute time in Gyr
        prefix = ''
        suffix = tunit
    else: ## 
        time_to_use = 1./time - 1. ## redshift
        prefix = 'z='
        suffix = ''
    n_sig = 2
    if time_to_use >= 10: n_sig += 0
    if time_to_use < 1.: n_sig += 1
    t_str = round_to_n( time_to_use, n_sig+n_sig_add )
    label_str = prefix+t_str+suffix
    return label_str

def round_to_n(x, n):
    ''' Utility function used to round labels to significant figures for display purposes  from: http://mail.python.org/pipermail/tutor/2004-July/030324.html'''
    if n < 1:
        raise ValueError("number of significant digits must be >= 1")

    # show everything as floats (preference; can switch using code below to showing eN instead
    format = "%." +str(n-1) +"f"
    as_string=format % x
    return as_string
