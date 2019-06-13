function  morletfft,s,n,k_mor=k_mor
  if not keyword_set(k_mor) then k_mor=6.
  foo=WV_FN_MORLET (k_mor,s,n,wavelet=w)
  return,w
end
function morletfft2,s,n,k_mor=k_mor
  t=findgen(n); fixed �������� �������, ���� �� ���� ������ �� �������
  if not keyword_set(k_mor) then k_mor=6.
  alp=1/(k_mor*sqrt(2)); ����������� � ������� �����
  return, sqrt(s)*exp(-((1-2*!dpi*s*t/n)/(2*alp))^2)*(1/alp)*sqrt(!dpi)/n
 ; return, sqrt(s)*exp(-((1-2*!dpi*s*t/n)/(2*alp))^2)*(1/alp)*sqrt(!dpi)/n
   ;���������� ������������ ������������� ��� ����������-����������� ������� �����
  ;FFT(f(x)*)=reverse(f(x))$ f(x)*--����������-����������� ������� f(x)
end