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
  ;   Num             - number of data points
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
pro emd_synthetic_model, t, x, x_clean, trend_clean, Num = Num, dt = dt, alpha = alpha, energy_white = energy_white, energy_color = energy_color,$
            energy_signal = energy_signal, period = period, non_stationary = non_stationary
  if n_elements(Num) eq 0 then Num=400l                ;Num= NUMBER OF DATAPOINTS
  if n_elements(dt) eq 0 then dt=0.1d              ;dt=CADENCE
  if n_elements(alpha) eq 0 then alpha=1.5d           ;COLOURED NOISE index
  t=findgen(Num)*dt            ;ARTIFICIAL TIME ARRAY
  
  if keyword_set(non_stationary) then begin
    trend = t^2*exp(-0.25*t)
    trend/=max(trend)
    signal=fltarr(Num)
    for i=Num/6.,Num-1 do signal[i]=exp(-((t[i]-max(t)/4.)/(1.*period))^2)*sin((2d*!Pi/((period)*(0.95^t[i])))*(t[i]-t[Num/6.]))
    ;signal=exp(-((t)/(max(t)*0.7))^2)*sin((2d*!Pi/((period)*(0.99^t)))*(t))
  endif else begin
    trend=exp(-((t - 0.5d*max(t))/(max(t)*1d))^2)       ;ARTIFICIAL TREND
    signal=sin((2*!Pi/Period)*t)*exp(-((t - 0.65d*max(t))/(max(t)*0.25))^2)        ;PRESCRIBED SINUSOIDAL SIGNAL
  endelse
 
  if n_elements(energy_white) eq 0 then energy_white = 0.1d                  ;Energy density in the    white   noise   relative to the trend energy density 
  if n_elements(energy_color) eq 0 then energy_color = 0.1d                 ;Energy density in the    colored noise   relative to the trend energy density 
  if n_elements(energy_signal) eq 0 then energy_signal = 0.2d                     ;Energy density of the oscillatory signal relative to the trend energy density 
  if n_elements(period) eq 0 then period = 7d                          ;PRESCRIBED PERIOD OF SIGNAL
  
  trend_energy = stddev(trend)^2
  s =100l
  
  An0=sqrt(energy_white*trend_energy)
  An1=sqrt(energy_color*trend_energy)
  
  noise_0=An0*random_alpha(s,0d,Num); white noise component
  noise_1=An1*random_alpha(s,alpha,Num); colored noise component
  As=sqrt(energy_signal*trend_energy)
 
  signal = As*signal/stddev(signal)

  x=noise_0+noise_1+trend+signal          ;TOTAL SIGNAL
  x_clean = trend+signal                  ;CLEAN SIGNAL + TREND
  trend_clean = trend                     ;TREND
end
