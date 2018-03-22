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
  ;
  ;
  ; :Authors: Dmitrii Kolotkov (D.Kolotkov@warwick.ac.uk) & Sergey Anfinogentov (anfinogentov@iszf.irk.ru) 
  ;-
pro emd_synthetic_model, t, x, x_clean, trend_clean
  N=400l                ;N= NUMBER OF DATAPOINTS
  cad=0.1d              ;CAD=CADENCE
  alpha_pres=1.5d           ;COLOURED NOISE indec
  t=findgen(N)*cad            ;ARTIFICIAL TIME ARRAY
  trend=exp(-((t - 0.5d*max(t))/(max(t)*1d))^2)       ;ARTIFICIAL TREND
  
  en0 = 0.1d                  ;Energy density in the    white   noise   relative to the trend energy density 
  en1 = 0.1d                  ;Energy density in the    colored noise   relative to the trend energy density 
  ens = 0.2d                  ;Energy density of the oscillatory signal relative to the trend energy density 
  Per=7d                      ;PRESCRIBED PERIOD OF SIGNAL
  
  trend_energy = stddev(trend)^2
  s =100l
  
  An0=sqrt(en0*trend_energy)
  An1=sqrt(en1*trend_energy)
  
  noise_0=An0*random_alpha(s,0d,N); white noise component
  noise_1=An1*random_alpha(s,alpha_pres,N); colored noise component
  As=sqrt(ens*trend_energy)
  

  signal=sin((2*!Pi/Per)*t)*exp(-((t - 0.65d*max(t))/(max(t)*0.25))^2)        ;PRESCRIBED SINUSOIDAL SIGNAL
  signal = As*signal/stddev(signal)

  x=noise_0+noise_1+trend+signal          ;TOTAL SIGNAL
  x_clean = trend+signal                  ;CLEAN SIGNAL + TREND
  trend_clean = trend                     ;TREND
end
