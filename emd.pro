;+
; function emd, data, shiftfactor = shiftfactor, maxsiftings = maxsiftings, epsilon = epsilon
;
; :DESCRIPTION:
;    Given a time series 'data', this function returns an empirical mode decomposition (EMD) of the signal
;    in a basis derived directly from the data by iterative searching for the local time scales naturally
;    appearing in the signal. The modes, or intrinsic mode functions (IMFs), are returned in a 2d array
;    with elements of length equal to data, and second index referring to the IMF. Usually returns between
;    3-15 modes. Shiftfactor controls sensitivity of the decomposition via limiting the standard deviation
;    between two subsequent siftings, and should be varied. For references/more information use links below.
;    http://adsabs.harvard.edu/abs/1998RSPSA.454..903E, http://adsabs.harvard.edu/abs/2008RvGeo..46.2006H
;
; :PARAMETERS:
;    data         is the series to be decomposed. Assumed to be a time series of equal step spacing. Note
;                 the data should have mean zero.
;
;    shiftfactor  is the value of the standard deviation between two subsequent siftings, controlling the
;                 sensitivity of the decomposition, and should be tested. Default = 0.2. In some cases,
;                 it may be reduced to 0.001-0.01.
;
;    maxsiftings  is the maximum number of siftings for each IMF. Has default= 4096. Suggest leave at
;                 default. Terminates sifting if shiftfactor set too low and cannot be achieved.
;
;    epsilon      is the lowest modal amplitude, has default 0.001. As with maxsiftings, sifting terminates
;                 if unable to find IMFs with the chosen shiftfactor. Suggest leave at default.
;
;  :RETURNS:
;    modes        is the two dimensional array containing the modes (IMFs).

;  :CALLING SEQUENCE:
;    Should remove mean signal i.e. IDL> data4emd = data - mean(data) then call
;     IDL> modes = emd(data4emd)
;
;  :SUBSEQUENT PLOTTING:
;    Get number of modes (IMFs). The last, or superposition of the last few modes, will be trend:
;     IDL> n_modes=(size(modes))[2].
;    Plot data against 'detrended' data i.e. minus final mode.
;     IDL> window,1,title='Detrended data' & plot,data,thick=2
;     IDL> oplot,modes[*,n_modes-1]+mean(data),color=55
;    After this should investigate each mode separately. Use commands like:
;     IDL> window,title='Mode 5 (red) against data (white) and detrended data (blue)'
;     IDL> plot,data4emd,thick=2 & oplot,data4emd-modes[*,n_modes-1],color=55 & oplot, modes[*,5],color=155
;
;    Number of extrema can be used to estimate mean period of mode (note mode is NOT necessarily sinusoidal,
;    arbitrary waveform).
;     IDL> length = float(n_elements(data4emd))     ; First get length of signal (in signal units)
;     IDL> n_ext = fltarr(n_modes-1) & for i=0,n_modes-1 do n_ext[i]=float((size(extrema(modes[*,i])))[1])
;     IDL> mode_period = 2d*length / n_ext          ; Estimate periods
;     IDL> print,mode_period ; doubling (dyadic) behaviour may indicate the random nature of these IMFs.
;    Breaks from doubling may imply real (non random) signal. See the below reference for further details.
;    http://adsabs.harvard.edu/abs/2004ISPL...11..112F, http://adsabs.harvard.edu/abs/2004RSPSA.460.1597W )
;
;  :EMD NOISE TESTING:
;    If working with realistic observational signals where quasi-periodic oscillatory phenomena are superimposed
;    on random frequency-dependant background processes (white or coloured noise), NOT ALL EMD modes will have
;    statistically significant properties. Usually, significant components are indicated among the whole set of
;    intrinsic modes by larger amplitudes (not counting the background slowly varying trend) and deviations from
;    dyadic (doubling) behaviour. For more detailed scheme of distinguishing between clear oscillations and noise
;    in EMD, see http://adsabs.harvard.edu/abs/2016A%26A...592A.153K.
;
;    This method normalises data before applying EMD + compares relative IMF amplitudes to estimate significance:
;     IDL> data4emd = data - mean(data)
;     IDL> data4emd = data4emd/sqrt(total(data4emd^2))
;     IDL> modes = emd(data4emd) & n_modes = (size(modes))[2] ; Also repeat estimating period from above
;     IDL> mode_energy = dblarr(n_modes) & for i=0,n_modes-1 do mode_energy[i]=total(modes[*,i]^2)
;     IDL> window,2,title='Logarithmic normalised mode energies'
;     IDL> plot,mode_period,mode_energy,/xlog,/ylog,psym=8,yrange=[1d-6,1],/ystyle,xtitle='Mean period',ytitle='Relative IMF energy'
;    Different noises give rise to modes lying on straight lines in log space. Look for deviation from straight lines.
;
;  :History:
;    Uploaded to CFSA library 07/08/2017
;
;  :Authors:
;    Sergey Anfinogentov (Sergey.Anfinogentov@warwick.ac.uk), Tim Duckenfield (T.Duckenfield@warwick.ac.uk)
;    & Dmitrii Kolotkov (D.Kolotkov@warwick.ac.uk)
;-
function spline_interp,x,y,x_out
  Y2 = SPL_INIT(X, Y)
  return,SPL_INTERP(X, Y, Y2, X_out)
end
pro sift,imf,local_mean
   n = n_elements(imf)
   foo = extrema(imf, max = maxpos, min = minpos )
   n_max= n_elements(maxpos)
   n_min = n_elements(minpos)
   ; Period of beginning wave
    period0 = 2 * abs ( maxpos[0] - minpos[0] )
    ; Period of end wave
    period1 = 2 * abs( maxpos[n_max-1] - minpos[n_min-1] )
    maxval = imf[maxpos]
    minval = imf[minpos]
    max0 = imf[maxpos[0]]
    min0 = imf[minpos[0]]
    max1 = imf[maxpos[n_max-1]]
    min1 = imf[minpos[n_min-1]]
    maxval = [max0,max0,imf[maxpos],max1,max1]
    minval = [min0,min0,imf[minpos],min1,min1]
    maxpos = [[-2l,-1l]*period0, maxpos, maxpos[n_max-1]+[1l,2l]*period1]
    minpos = [[-2l,-1l]*period0, minpos, minpos[n_min-1]+[1l,2l]*period1]
   
   maxenv = spline_interp( maxpos, maxval, dindgen(n))
   minenv = spline_interp( minpos, minval, dindgen(n))
;   maxenv = spline( maxpos, maxval, dindgen(n))
;   minenv = spline( minpos, minval, dindgen(n))
;   maxenv = interpol(maxval,maxpos, dindgen(n),/q)
;   minenv = interpol(minval,minpos, dindgen(n),/q)
   local_mean = ( maxenv + minenv )*0.5d
   imf -= local_mean
end

function next_imf,data,last_imf, no_residual, shiftfactor = shiftfactor, maxsiftings = maxsiftings
  n = n_elements(data)
  last_imf = 0
  residual = dblarr(n)
  imf = data
  n_sift = maxsiftings
  sd_all = dblarr(n_sift)
  sqr = total(data^2)
  previous_imf = dblarr(n)
  for i =0, n_sift -1 do begin
    previous_imf = imf
    sift,imf,local_mean
    ;if sqr_new gt sqr then break
    sd = total(local_mean^2/previous_imf^2)
    if sd lt shiftfactor then break
    if i mod 5 eq 0 then begin
;      plot,data
;      oplot, residual, color =250
;      wait,0.01
;      print,i
    endif
    ;print,i
  endfor
  ;print, 'IMF found for '+strcompress(i)+' siftings'
  return,imf

end

function emd,data, shiftfactor = shiftfactor, maxsiftings = maxsiftings, epsilon = epsilon
  if not keyword_set(shiftfactor) then shiftfactor =0.2
  if not keyword_set(maxsiftings) then maxsiftings =4096
  if not keyword_set(epsilon) then epsilon =1d-4
  x = double(data)
  n = n_elements(data)
  n_modes = 50
  result = dblarr(n,n_modes)
  data_std_dev = stddev(data)
  for i =0, n_modes -1 do begin
    imf = next_imf(x,last_imf, no_residual, shiftfactor = shiftfactor, maxsiftings = maxsiftings)  
    result[*,i] = imf
    x-= imf
    if stddev(x)/data_std_dev lt epsilon then begin
      result[*,i] += x
      result = result[*,0:i]
      break
    endif 
  endfor
  
  
  return, result
  
  
  
end

