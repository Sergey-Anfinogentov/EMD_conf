function morletcor,s,n,k_mor=k_mor
if not keyword_set(k_mor) then k_mor=6.
  ;foo=WV_FN_MORLET (k_mor,float(s),float(2*n),wavelet=w,/spat)
  ;res=fastcon(w,sin(2*!pi*findgen(n)/s))
  ;xx=[complex(cos(2*!pi*findgen(n)/s),sin(2*!pi*findgen(n)/s)),complexarr(n)]
  w=morletfft(s,2*n,k_mor=k_mor)
  w=fft(w,+1)/sqrt(s)
  xx=[fltarr(n)+1.0,fltarr(n)]
  
  res=n*2*fft(fft(xx,-1)*fft(abs(w),-1),+1)

  ;res=shift(res,-n)
  ;plot,abs(res)
  res=abs(res(0:n-1))
  ;print,max(res)
  k=fltarr(20)
  k=[7.6238909,0]
  return,res/max(res)
end
