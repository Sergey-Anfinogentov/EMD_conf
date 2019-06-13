function logspace, x1, x2, n
  compile_opt idl2
  if x1 le 0 or x2 le 0 then message,'Both limits must be positive'
  dx = (double(x2)/double(x1))^(1d/(n-1d))
  return, product([x1,replicate(dx,n-1)], /cumulative)
end



function at_emd, x, show = show
 
  x_original = x
  x = x - mean(x)
    
  modes = emd(x)
  sz = size(modes)
  n = sz[2]
  n_sh = 100
  
  
 ; stop
  shiftfactors = logspace(0.01d,0.5,n_sh)
  
 chr = dblarr(n_sh) 
 result = list()
 threshold = 0.03
 var_orig = stddev(x)^2
 if keyword_set(show) then begin
   window, 0
   plot, x,/psym
 endif
 
 for k = 0, 100 do begin
  modes_list = list()
  
   for i =0, n_sh-1 do begin
     modes = emd(x, shiftfactor = shiftfactors[i])
     
     en = total(modes^2,1)
     nm = (size(modes))[2]
     chr_0 = dblarr(nm)
     for j =0, nm-1 do begin
      chr_0[j] =  total((modes[*,j] - x)^2)
     endfor
     foo = min(chr_0, mode_ind)    
      mode = modes[*,mode_ind] 
     modes_list.add, mode
     
     chr[i] = total((mode - x)^2)
     print,k, i, shiftfactors[i]
   endfor
   
   foo = min(chr, ind)
   mode = modes_list[ind]
   sh_used = shiftfactors[ind]
   
   result.add, mode
   
   if keyword_set(show) then begin
     window, k
     plot, x,/psym
     oplot, mode
   endif
   
   x -= mode
   print,'residual energy:',stddev(x)^2/var_orig, ', shift factor used: ', sh_used
   if stddev(x)^2 lt threshold*var_orig then break
 endfor
 
 ;wdel, alll
  
  nm = n_elements(result)
  nt = n_elements(result[0])
  modes = dblarr(nt, nm + 1)
  for i =0,nm-1 do begin
    modes[*,i] = result[i]
  endfor
  modes[*,nm] = x
  
   sp = emd_energy_spectrum(modes)
   
   ind = sort(sp.period)
   
   modes = modes[*,ind]
   sp = emd_energy_spectrum(modes)
   plot,sp.period,sp.energy, /xlog, /ylog,psym = 1,$
     yrange = minmax(sp.energy)*[1d/20d,20d], title = 'EMD Power spectrum'
     
 x = x_original 
return, modes

end