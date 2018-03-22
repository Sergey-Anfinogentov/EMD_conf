;+
  ; :Description:
  ;    This procedure creates a synthetic sugnal to be used as an example
  ;
  ; :Params:
  ;    t            - time (out)
  ;    x            - synthetic signal (out)
  ;    x_clean      - synthetic signal without noise
  ;    trend_clean  - prescribed trend
  ;    
  ;  :Keywords:
  ;   N             - number of data points
  ;   dt            - cadence time
  ;   alpha         - Power law index of the coloured noise
  ;   energy_white  - Energy of the white noise relative to the  energy of the trend
  ;   energy_color  - Energy of the coloured noise relative to the  energy of the trend
  ;   energy_signal - Energy of the oscillatory component relative to the  energy of the trend
  ;   period        - Period of the oscillatory component
  ;
  ;
  ;
  ; :Authors: Dmitrii Kolotkov (D.Kolotkov@warwick.ac.uk) & Sergey Anfinogentov (anfinogentov@iszf.irk.ru) 
  ;-
pro emd_synthetic_model, t, x, x_clean, trend_clean, N=N, dt = dt, alpha = alpha, energy_white, energy_color,$
            energy_signal = energy_signal, period = period
  if not keyword_set(N) then N=400l                ;N= NUMBER OF DATAPOINTS
  if not keyword_set(dt) then dt=0.1d              ;dt=CADENCE
  if not keyword_set(alpha) then alpha=1.5d           ;COLOURED NOISE indec
  t=findgen(N)*dt            ;ARTIFICIAL TIME ARRAY
  trend=exp(-((t - 0.5d*max(t))/(max(t)*1d))^2)       ;ARTIFICIAL TREND
  
  if not keyword_set(energy_white) then energy_white = 0.1d                  ;Energy density in the    white   noise   relative to the trend energy density 
  if not keyword_set(energy_color) then energy_color = 0.1d                 ;Energy density in the    colored noise   relative to the trend energy density 
  if not keyword_set(energy_signal) then energy_signal = 0.2d                     ;Energy density of the oscillatory signal relative to the trend energy density 
  if not keyword_set(period) then period = 7d                          ;PRESCRIBED PERIOD OF SIGNAL
  
  trend_energy = stddev(trend)^2
  s =100l
  
  An0=sqrt(energy_white*trend_energy)
  An1=sqrt(energy_color*trend_energy)
  
  noise_0=An0*random_alpha(s,0d,N); white noise component
  noise_1=An1*random_alpha(s,alpha,N); colored noise component
  As=sqrt(energy_signal*trend_energy)
  

  signal=sin((2*!Pi/Period)*t)*exp(-((t - 0.65d*max(t))/(max(t)*0.25))^2)        ;PRESCRIBED SINUSOIDAL SIGNAL
  signal = As*signal/stddev(signal)

  x=noise_0+noise_1+trend+signal          ;TOTAL SIGNAL
  x_clean = trend+signal                  ;CLEAN SIGNAL + TREND
  trend_clean = trend                     ;TREND
end
