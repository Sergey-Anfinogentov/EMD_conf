;+
; :Description:
;    Finds maximum position with subpixel resolution
;    
;
; :Params:
;    var: 1d or 2d array
;    coord: coordinate
;
;
;
; :Author: Sergey Anfinogentov (email: anfinogentov@iszf.irk.ru)
;-
Function maxi,var,coord
  si=size(var)
  if si(0) eq 1 then begin
    n=si(1)
    foo=max(var,xm)
    xst=xm-10
    xen=xm+10
    if xst lt 0 then xst=0
    if xen ge n then xen=n-1
    ivar= Congrid(var(xst:xen),(xen-xst+1)*10,cubic=-0.5,/minus_one)
    res=max(ivar,xmi)
    coord=float(xst)+xmi/10.
    return,res    
  endif
  if si(0) eq 2 then begin
    nx=si(1)
    ny=si(2)
    foo=max(var,ind)
    xm=ind mod nx
    ym=floor(ind/float(nx))
    xst=xm-10
    xen=xm+10
    yst=ym-10
    yen=ym+10
    if xst lt 0 then xst=0
    if xen ge nx then xen=nx-1
    if yst lt 0 then yst=0
    if yen ge ny then yen=ny-1
    snx=(xen-xst+1)*10
    sny=(yen-yst+1)*10
    ivar= Congrid(var(xst:xen,yst:yen),snx,sny,cubic=-0.5,/minus_one)
    res=max(ivar,ind)
    xmi=ind mod snx
    ymi=floor(ind/float(snx))
    coord=[float(xst)+xmi/10.,float(yst)+ymi/10.]
    return,res    
  endif
end