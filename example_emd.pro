
pro example_emd
  emd_synthetic_model,t,x, x_clean, trend_clean
  dt = t[1] - t[0]
  
  x_original = x
  
  ;Caculate EMD modes
  modes = emd(x,shiftfactor=0.15,  maxsiftings=1d5)
  
  ;Calcilate trend
  trend_emd = modes[*,-1]
  
  
  ;-------PLOTTONG------------
  window,0,  xsize = 1000, ysize =700
  !p.BACKGROUND =255
  !p.color = 0
  !p.multi = [0,2,2]
  plot, t, x, ystyle = 1, title = 'Original signal'
  loadct, 39
  oplot, t, x_clean,thick =2, color =64
  loadct, 39
  oplot, t, trend_clean,thick =1, color =250;, linest = 1
  ;------- END of PLOTTING---------
  
  ;subtract trend from the signal
  x=x-trend_emd 
  
  ;Estimating noise parameters from the tetrended signal
  fit_fft = fft_alpha(x,dt,fap =0.005)

  ;convert frequency to period
  period = 1d/fit_fft.frequency
  
  ;------PLOTTING--------------
  plot,period, fit_fft.power,/xlog,/ylog,xstyle =1, ystyle =1$ 
    ,yrange = minmax(fit_fft.power)*[0.1d,10d], xrange = minmax(period)*[0.5d,2d],$
     title = 'FFT Power'
  oplot,period, fit_fft.expectation, color = 64, thick =2
  oplot,period,fit_fft.confidence_level,  color =250
  ;------END OF PLOTTING-------
  
  ;colored noise index
  alpha = fit_fft.pl_index
  
  ; Confidence intervals for coloured noise
  conf_c = emd_noise_conf(n_elements(x),alpha,nsampl =100, energy = fit_fft.color_energy, fap = 0.005d)
  
  ; Confidence intervals for the wight noies
  conf_w = emd_noise_conf(n_elements(x),0d,nsampl =100, energy = fit_fft.white_energy, fap = 0.005d)
  
  ;Calculate EMD power spectrum
  sp = emd_energy_spectrum(modes)
  
  ;Upper confidence interval for the combined noises
  conf_up = conf_c.up + conf_w.up
  
  ;Lower confidence interval for the combined noises
  conf_period = conf_c.period
  
  ;Recover "CLEAN" signal by summing up significant modes and the trend
  emd_clean = emd_recover_signal(modes, conf_period, conf_up)
  
 ;------PLOTTING-------------- 
 !p.multi = [1,2,2] 
  plot,sp.period*dt,sp.energy, /xlog, /ylog,psym = 1,$
     yrange = minmax(sp.energy)*[1d/20d,20d], title = 'EMD Ppower spectrum'
     
  
  oplot, conf_period*dt,conf_c.mean_energy + conf_w.mean_energy, color =64, thick =2 
  oplot,conf_period*dt,conf_c.up + conf_w.up, color =250
  oplot,conf_period*dt,conf_c.down + conf_w.down, color =250

  
  !p.multi = [2,2,2]
  
  plot, t,x_original, title = 'Recovered Signal',ystyle = 1
  oplot, t, emd_clean + trend_emd, color =64, thick =2
  oplot, t, trend_emd,thick =1, color =250;, linest = 1
  !p.multi = 0
  ;---PLOTTING-----------

end