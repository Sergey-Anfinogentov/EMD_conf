function emd_logspace, x1, x2, n
  compile_opt idl2
  if x1 le 0 or x2 le 0 then message,'Both limits must be positive'
  dx = (double(x2)/double(x1))^(1d/(n-1d))
  return, product([x1,replicate(dx,n-1)], /cumulative)
end

;+
  ; :Description:
  ;    Fits a linear function in log-log scale.
  ;    This function should not be called directly
  ;
  ; :Params:
  ;    period   - the EMD period (independend variable)
  ;    value    - Dependend variable (Enrgy or DoF)
  ;    period1  - Period interval to make a make a fitted dependence
  ;    period2  - Period interval to make a make a fitted dependence
  ;
  ; :Keywords:
  ;    n        - the length of the array to generate
  ;
  ; :Author: Sergey Anfinogentov (anfinogentov@iszf.irk.ru)
  ;-
function make_lin_dep_loglog, period, value, period1, period2,n = n
  If not keyword_set(period1) then period1 = min(period)
  If not keyword_set(period2) then period2 = max(period)
  if not keyword_set(n) then n = 500l
  
  par = ladfit(alog(period), alog(value))
  
  p = emd_logspace(period1, period2, n)
  
  v = exp(par[0] + par[1]*alog(p))
  
  return,{period:p, value:v}
end


;+
  ; :Description:
  ;    Computes modified Chisqr probability distribution function$
  ;    with the expectation of e_mean and DoF degrees of freedom
  ;
  ; :Params:
  ;    x      - independend variable i.e. Energy
  ;    e_mean - energy expectation, i.e. average model energy
  ;    dof    - number of the degree of freedom
  ; 
  ; :Output:
  ;   Return the probability density
  ;
  ;
  ; :Author: Sergey Anfinogentov (anfinogentov@iszf.irk.ru)
  ;-
function chisqr_emd,x,e_mean,dof
  h = 1d-5
  x1  = x - h*0.5d
  x2 = x + h*0.5d
  
  y1 = dof*(x1/e_mean)
  y2 = dof*(x2/e_mean)
  
  h_y = (y2 - y1)
  
  return,  (h_y/h)*(chisqr_pdf(y2,dof) - chisqr_pdf(y1,dof))/h_y
end

;+
  ; :Description:
  ;    A function to be used for fitting
  ;    This function should not be called directly
  ;
  ; :Params:
  ;    x - [in] Energy
  ;    a - [in] parameter set
  ;    f - [out] probability density
  ;
  ;
  ;
  ; :Author: Sergey Anfinogentov (anfinogentov@iszf.irk.ru)
  ;-
pro chisqr_fun, x, a, f
  mean_energy = abs(a[0])
  dof = abs(a[1])
  f = chisqr_emd(x,mean_energy,dof)
end

;+
  ; :Description:
  ;    Fits modified chisqr PDF in energy histograms
  ;
  ; :Params:
  ;    energy - Centers of the energy bins
  ;    hist   - Energy histogram
  ;
  ;
  ;
  ; :Author: Sergey Anfinogentov (anfinogentov@iszf.irk.ru)
  ;-
function chisqr_fit,energy, hist
  weights=replicate(1.0, N_ELEMENTS(hist))
  
  foo = max(hist,ind)
  parms = [energy[ind], 1d]
  
  h_fit=CURVEFIT( energy, hist, weigths, parms, fita=[1,1], FUNCTION_NAME= 'chisqr_fun',/noderivative, status = status)
  


  return,{mean_energy:abs(parms[0]),dof:abs(parms[1])}

end



;+
  ; :Description:
  ;    Processes synthetic EMD modes and estimates mean_energy,
  ;    mean_period and number of degrees of freedom for every mode number
  ;
  ; :Params:
  ;    p        - [in] mode periods
  ;    e        - [in] mode energies
  ;    n_mode   - [in] mode number
  ;
  ; :Outup:
  ;   returns a structure with the following fields:
  ;     - mean_energy
  ;     - DoF
  ;     - period
  ;
  ; :Author: Sergey Anfinogentov (anfinogentov@iszf.irk.ru)
  ;-
function emd_noise_fit, p, e, n_mode
  max_modes = max(n_mode)-1
  
  dof  = dblarr(max_modes) ; DoF
  energy = dblarr(max_modes) ; Mean modal energy
  period = dblarr(max_modes) ; Mean modal period
  
  for i = 1, max_modes do begin
    ind  = where(n_mode  eq  i)
    period_i = p[ind]
    energy_i = e[ind]
    

    period[i-1] = mean(period_i)
    
    ;calculating and normalising energy histogram
    h = histogram(energy_i,loc = loc,nbins =200)
    bin = loc[1]-loc[0]
    h = h/total(bin*h)
    
    ;fitting Chisqr PDF in the histogram
    par  = chisqr_fit(loc+bin*0.5, h)
    energy[i-1] = par.mean_energy
    dof[i-1] = par.dof
  endfor
  
  return,{mean_energy:energy,dof:dof, period: period}

end 

;conf_intervals**********************************

;+
; :Description:
;    Computes confidense intervals for EMD power spectrum
;
; :Params:
;    n      - [in] length of the signal in data points 
;    alpha  - [in] Coloured noise index estimated with FFT_ALPHA
;
; :Keywords:
;    energy           - [in] Energy density of the noise estimated with FFT_ALPHA
;    nsamples         - [in] Number of synthetic signals to generate (defualt 500)
;    fap              - [in] False Alarm Probability
;
; :Output:
;   Returns a structure with the following fields:
;       up            - upper confindence level as a function of period
;       down          - bottom confidence level as a function of period
;       DoF           - number of degrees of freedom as afunction of period
;       mean_energy   - average modal energy as a function of period
;       period        - period
;
; :Author: Sergey Anfinogentov (anfinogentov@iszf.irk.ru) &
;          Dmitrii Kolotkov (D.Kolotkov@warwick.ac.uk)
; 
;-
function emd_noise_conf, n, alpha, energy = energy, nsamples = nsamples, fap = fap

  if not keyword_set(energy) then energy = 1d ; energy density
  if not keyword_set(nsamples) then nsamples =500l
  if not keyword_set(fap) then fap = 0.01d
  confidence_level = 1d - fap*0.5
  
  nmoden = []
  pn=[]
  en=[]
  s = systime(1) ; random seed
  for j=0,nsamples-1 do begin
    
    ; generating random noise
    fnoise=random_alpha(s,alpha,N)
    fnoise -= mean(fnoise)
    fnoise/=stddev(fnoise)
    fnoise*=sqrt(energy)

    tmp=emd(fnoise,shi=0.3)
    ;estimating energy and period for each mode
    for i=0,(size(tmp))[2]-1 do begin
      length=float(n_elements(fnoise))
      n_ext=float((size(extrema(tmp[*,i])))[1])
      p=(emd_period_energy(tmp[*,i])).period
      e=(stddev(tmp[*,i]))^2
      pn = [pn, p]
      en = [en, e]
      nmoden = [nmoden,i+1]
    endfor
    message,'processing noise samples: '+strcompress(j)+ ' of '+strcompress(nsamples-1),/info 
  endfor

  p = pn
  e = en
  nmode = nmoden
  
  noise =emd_noise_fit(p,e,nmode)
  
  n_conf =500l
  ind = where(noise.period lt n/2)
  ind = ind[1:*]
  
  foo = make_lin_dep_loglog(noise.period[ind], noise.dof[ind], 2,n,n = n_conf)
  period = foo.period
  dof = foo.value
  foo = make_lin_dep_loglog(noise.period[ind], noise.mean_energy[ind], 2,n,n = n_conf)
  mean_energy = foo.value
  
  n = n_elements(dof)
  up = dblarr(n)
  down = dblarr(n)

  for i = 0, n-1 do begin
    down[i] = chisqr_cvf(confidence_level,dof[i])*mean_energy[i]/dof[i]
    up[i] = chisqr_cvf(1d -confidence_level,dof[i])*mean_energy[i]/dof[i]
  endfor


return,{up:up, down:down,dof:dof, mean_energy: mean_energy, period:period}

end