;+
; :Description:
;    This function calulate  Global Wavelet Spectrum using  fast global wavelet spectrum computation algorithm
; :Categories:
;   Pixelize Wavelet Filtration library
; :Returns:
;   Filtration result as fltarr(period)  
; :Params:
;    x: input,required,type = fltarr(nx,nt,ntime)
;      3d data to analize 
;    t1: input,optional, default = 2
;     The low limit of period interval to analyze
;    t2: input,optional, default = n/3
;     The high limit of period interval to analyze
;
; :Keywords:
;    period: out
;     named variable for period values output
;    k_mor: in, default = 6
;    Morlet wavelet parameter
;    counts: in, default = 50
;    Period counts number
;     normalize: Set this keyword to do frequency normalization in spectrum
;     nopad: Set this keyword not to use zero padding
; 
;     linpad: Set this keyword to use linear interpolation padding
;     instead of zero padding. This option can reduce edge effect 
;     and fix some problems with low frequencies
; :Examples: 
; sp=pwf_gws(x,3,100,period=period)
; plot,period,sp   
; :Uses:
;  morletfft  
; :Author: Sergey Anfinogentov
;-
function pwf_gwsfast,x,t1,t2,period=period,k_mor=k_mor,counts=counts,normalize=normalize,nopad=nopad,linpad=linpad,test=test
  if not keyword_set(k_mor) then k_mor=6.
  if n_elements(counts) eq 1 then num=counts else num=50
  nt=(size(x))(3)
  nx=(size(x))(1)
  ny=(size(x))(2)
  if not keyword_set(t1) then t1=2.
  if not keyword_set(t2) then t2=nt/3.
  foo=WV_FN_MORLET (k_mor,20.0,0)
  sen=float(t2)/foo.fourier_period; period to scale
  sst=float(t1)/foo.fourier_period
  ds=(sen/sst)^(1/float(num-1))
  s=sst
  ;res=fltarr(nt,num)
  linpad=1
  if keyword_set(nopad)then begin
    xx=x
    nnt=nt
  endif else begin
     xx=pwf_expand(x-mean(x),linpad=linpad)
     nnt=2*nt

  endelse
  if (size(xx))[0] eq 1 then begin
    xx=reform(xx,1,1,n_elements(xx))
     nt=(size(xx))(3)
    nx=(size(xx))(1)
    ny=(size(xx))(2)
  endif
  fx=fft(xx,dim=3)
  fxm=sqrt(total(total(abs(fx)^2,1),1))
  ;fxm=total(total(abs(fx),1),1)
  fx0=fx
  fx=fxm;*exp(complex(0,randomu(systime(/sec)*2.*!pi)))
  res=fltarr(num)
   period=fltarr(num)
  for i=0,num-1 do begin
    fm=morletfft(s,nnt,k_mor=k_mor)  
    period(i)=s*foo.fourier_period
    if keyword_set(normalize) then fm=fm/sqrt(s)
    ;fm=rebin(reform(fm,1,1,nnt),nx,ny,nnt)
    s=s*ds
    if keyword_set(test) then begin
      fm=rebin(reform(fm,1,1,nnt),nx,ny,nnt)
      sig=total(abs(fx0*fm)^2,3)
      back=total(abs(fx0*(1.-fm))^2,3)
      res[i]=mean(sig/back)  
      print,period[i]
      tvcon,sig/back
      stop
      continue
    endif
    ;res(*,i)=2*total(total((abs(fft(fx*fm,+1,dim=3))^2)(*,*,0:nt-1),1),1)
    res(i)=total(abs(fft(fx*fm,+1))^2)/(nx*ny)
  endfor
  return,res;;total(abs(res)^2,1)
end
