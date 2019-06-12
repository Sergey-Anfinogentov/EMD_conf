;+
    ; :Description:
    ;    A simple gradient descent maximisation algorithm.
    ;    Used bt FFT_ALPHA
    ;
    ; :Params:
    ;    fun_name - name of the function to maximised
    ;    parms - initial guess of the dunction parameters
    ;
    ; :Keywords:
    ;    x - independend variable values, will be passed to the fun_name function
    ;    y - onserved values of the dependend values be passed to the fun_name function
    ;
    ; :Author: Sergey Anfinogentov (anfinogentov@iszf.irk.ru)
    ;-
pro grad_max, fun_name, parms, x = x, y = y
  k = 1d-2
  max_iter =5000
  current = call_function(fun_name,parms,dp,x=x,y=y)
  for i =0, max_iter do begin
    new_parms = parms + dp*k
    new = call_function(fun_name,new_parms,dp,_extra={x:x,y:y})
    if new lt current then  k*=0.7d else begin
      rate = (new - current)/current
      current = new
      parms = new_parms
      if abs(rate) lt 1d-15 then break 
      k*= 1.1d
    endelse
  endfor
  print,i

end

function model_spectrum,x,p
  white_energy =  exp(p[0]);>0
  color_energy =  exp(p[1]);>0
  color_index  =  p[2]
  nf = n_elements(x)

  colored_noise = x ^ (-color_index)
  colored_noise = color_energy* colored_noise/mean(colored_noise)

  model = white_energy + colored_noise
  return, model
end

;+
  ; :Description:
  ;    Calculates the log-likelihood function for the following Fourier power spectrum model:
  ;    WHITE_NOISE + COLOURED_NOISE + POSSIBLE_OUTLIERS
  ;
  ; :Params:
  ;    p        - [in] free parameters = [alog(WHITE_NOISE_ENERGY), alog(COLOURED_NOISE_ENERGY), COLOURED_NOISE_INDEX]
  ;    dp       - [out] gradient of the function
  ;
  ; :Keywords:
  ;    x        - [in] independend variable i.e. frequency 
  ;    y        - [out] depended variable i.e.
  ;    model    - [out] modeled spectrum
  ;    no_grad  - [in] If the keyword is set, gradient DP is not computed
  ;
  ; :Author: Sergey Anfinogentov (anfinogentov@iszf.irk.ru)
  ;-
function noise_model_alpha_tn, p, dp, x=x, y=y, model = model, no_grad = no_grad
  white_energy =  exp(p[0]);>0
  color_energy =  exp(p[1]);>0
  color_index  =  p[2]
  nf = n_elements(x)
  
  colored_noise = x ^ (-color_index)
  colored_noise = color_energy* colored_noise/mean(colored_noise)
  
  model = white_energy + colored_noise
  max_y = max(y)
  min_y = 0d
  
  outlier_prob = 1d/nf
  outlier_density = outlier_prob/(max_y-min_y)
  
  temp = exp(-y/model)/(model) * (1d - outlier_prob) + outlier_density
  result = total(alog(temp))
  
  if not keyword_set(no_grad) then begin
    dp = dblarr(3)
    h = (abs(p[0])*1d-8)>1d-10
    p1 = p+[h,0d,0d]
    dp[0] = (noise_model_alpha_tn( p1, undefined, x=x, y=y,/no_grad) - result)/h
    h = (abs(p[1])*1d-8)>1d-10
    p1 = p+[0d,h,0d]
    dp[1] = (noise_model_alpha_tn( p1, undefined, x=x, y=y,/no_grad) - result)/h
    h = p[2]*1d-8
    p1 = p+[0d,0d,h]
    dp[2] = (noise_model_alpha_tn( p1, undefined, x=x, y=y,/no_grad) - result)/h
    ;print,[p,result,dp]
  endif
  return, result
end



;+
  ; :Description:
  ;    Estimates the parameters of 2 component (white + coloured) noise using FFT
  ;
  ; :Params:
  ;    signal - [in] signal to investigate 
  ;    dt     - [in] time cadence
  ;
  ; :Keywords:
  ;    fap    - [in] False alarm probability
  ;    
  ; :Outputs:
  ;     Returns a structure with the following fields
  ;         FREQUENCY         - FFT frequency
  ;         EXPECTATION       - Frequency dependend expectation of the FFT power
  ;         CONFIDENCE_LEVEL  - Frequency dependend confidence level based on
  ;                             the provided False alarm probability
  ;         white_energy      - Estimated energy density of the white noise component
  ;         color_energy      - Estimates energy density of the coloured noise component
  ;         pl_index          - Estimated Coloured noise index
  ;
  ; :Author: Sergey Anfinogentov (anfinogentov@iszf.irk.ru)
  ;-
function fft_alpha, signal, dt, fap = fap
  nt = N_ELEMENTS(signal)
  
  if not keyword_set(fap) then fap = 0.05d
  
  IF (nt mod 2) eq 0 THEN signal_ = [signal,0d] else signal_ = signal
  
  nt = N_ELEMENTS(signal_)
  nf = nt/2l
  k = DINDGEN(nf) +1d
  Frequency = k/(nt*dt)
  
  spectrum = FFT(signal_)
  pow = (nt/(nf))*abs(SPECTRUM[1:nf])^2
  
  
  parinfo = replicate({value:0d, fixed:0d, limited:[1d,1d], limits:[0d,0d],parname:'',TNSIDE:2},3)
  pow_alog = alog(pow)
  parinfo[0].value = alog(mean(pow))
  parinfo[0].limits = alog([min(pow),mean(pow)])+[-10d,0d]
  parinfo[0].parname = 'White noise energy'
  
  parinfo[1].value = alog(mean(pow))
  parinfo[1].limits = alog([min(pow),mean(pow)])+[-10d,0d]
  parinfo[1].parname = 'Colored noise energy'
  
  parinfo[2].value = 4d
  parinfo[2].limits = [0.5,3d]
  parinfo[2].parname = 'Power law index for coloured noise'

  
  FUNCTARGS = {x: frequency, y: pow}
   
   parms = parinfo.value
   grad_max, 'noise_model_alpha_tn', parms, x = frequency, y = pow
  
  foo = noise_model_alpha_tn(parms, x=frequency[2:*], y=pow[2:*], model = yfit)
  
  yfit = model_spectrum( frequency, parms) 
  
    threshold = 1d - (1d - fap)^(1d/nf)
    
    confidence_level = chisqr_cvf(threshold,2)*yfit/2d
    
    if parms[2] gt 3 then parms[2] =3
  
    message, 'white noise energy density:'+strcompress(exp(parms[0])*nf*(nt-1)/nt),/info;, exp(parinfo[0].limits)
    message, 'colored noise energy density:' +strcompress(exp(parms[1])*nf*(nt-1)/nt),/info;, exp(parinfo[1].limits)
    message, 'PL_index:'+strcompress((parms[2])),/info
    message, 'original signal energy density:'+strcompress(stddev(signal)^2),/info
    message, 'fitted spectrum energy density:'+strcompress(total(yfit)),/info


    
  return,{power:pow,$
     frequency: frequency,$
     expectation:yfit,$
     white_energy:exp(parms[0])*nf*(nt-1)/nt,$
     color_energy:exp(parms[1])*nf*(nt-1)/nt,$
     pl_index:parms[2],$
     confidence_level:confidence_level}
  
end