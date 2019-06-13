
;+
; :Description:
;    The functions pads 1d or 3d dataset with zero or with linear interpolated data.
;  :Categories: 
;   Pixelize Wavelet Filtration library
; :Returns: 
; The function returns global wavelet spectrum as 1d (fltarr(2*nt)) array
; :Params:
;    x: in, required, type: fltarr(nt) or fltarr(nx,ny,nt) 
;    input data 1d or 3d array
;
; :Keywords:
;    linpad: Set this keyword to use linear interpolation padding instead of zero padding
;
; :Author: Sergey Anfinogentov
;-
function pwf_expand,x,linpad=linpad,nopad=nopad
if keyword_set(nopad) then return,x
temp=size(x)
if temp(0) eq 3 then return,pwf_expand3(x,linpad=linpad)
  nt=temp(1)
  if keyword_set(linpad) then begin
    t=findgen(nt)
    b=x(nt-1)   
    k=(x(0)-b)/nt
    e=k*t+b
  endif else begin
    e=fltarr(nt)
  endelse
  return,[x,e]
end