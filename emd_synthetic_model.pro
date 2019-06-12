function emd_synthetic_trend, t, maximum_time = maximum_time, grow_time = graw_time, decay_time = decay_time
  if not keyword_Set(decay_time) then decay_time = (max(t) - min(t))*0.3d
  if not keyword_Set(maximum_time) then maximum_time = (max(t) - min(t))*0.3d +min(t)
  if not keyword_Set(grow_time) then grow_time = (max(t) - min(t))*0.2d
  result = t*0d
  grow = (t le maximum_time) * exp(-((t-maximum_time)/grow_time)^2)
  decay = (t gt maximum_time) * exp(-((t-maximum_time)/decay_time))
  return, grow + decay
end

function emd_synthetic_non_stationary, t, start_period = start_period, period_change_rate = period_change_rate, start_time = start_time
  duration = max(t) - min(t)
  if not keyword_Set(start_time) then start_time = duration*0.3d +min(t)
  
  if not keyword_set(start_period) then start_period = duration * 0.2d
  if not keyword_Set(period_change_time) then period_change_time = duration*0.5
  period = start_period * exp(-(t - start_time)/period_change_time)
  omega = 2d*!dpi/period
  osc = sin(omega  * (t - start_time))
  
  decay_time = (max(t) - min(t))*0.5d
  grow_time = (max(t) - min(t))*0.05d
  
  grow  = (t le start_time) * exp(-((t-start_time)/grow_time)^2)
  decay = (t gt start_time) * exp(-((t-start_time)/decay_time))
  
  envelope = (t le start_time) * exp(-((t-start_time)/grow_time)^2) + (t gt start_time) * exp(-((t-start_time)/decay_time)^2)
  
  return, osc*envelope
end

function emd_synthetic_stationary, t, period = period, start_time = start_time
  duration = max(t) - min(t)
  if not keyword_Set(start_time) then start_time = duration*0.3d +min(t)

  if not keyword_set(period) then period = duration * 0.15d

  omega = 2d*!dpi/period
  osc = sin(omega  * (t - start_time))

  decay_time = (max(t) - min(t))*0.3d
  grow_time = (max(t) - min(t))*0.05d

  grow  = (t le start_time) * exp(-((t-start_time)/grow_time)^2)
  decay = (t gt start_time) * exp(-((t-start_time)/decay_time))

  envelope = (t le start_time) * exp(-((t-start_time)/grow_time)^2) + (t gt start_time) * exp(-((t-start_time)/decay_time)^2)

  return, osc*envelope
end

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
    signal = emd_synthetic_non_stationary(t)
  endif else begin
    signal = emd_synthetic_stationary(t)     ;PRESCRIBED SINUSOIDAL SIGNAL
  endelse
  
  trend = emd_synthetic_trend(t)
 
  if n_elements(energy_white) eq 0 then energy_white = 0.1d                  ;Energy density in the    white   noise   relative to the trend energy density 
  if n_elements(energy_color) eq 0 then energy_color = 0.1d                 ;Energy density in the    colored noise   relative to the trend energy density 
  if n_elements(energy_signal) eq 0 then energy_signal = 0.05d                     ;Energy density of the oscillatory signal relative to the trend energy density 
  if n_elements(period) eq 0 then period = 7d                          ;PRESCRIBED PERIOD OF SIGNAL
  
  trend_energy = stddev(trend)^2
  s = random_seed();100l
  
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
