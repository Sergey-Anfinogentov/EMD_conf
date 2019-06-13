function emd_trend, modes
  sp = emd_energy_spectrum(modes)
  nt = (size(modes))[1]
  ind = where(sp.period gt 0.7*nt)
  result = modes[*,0]*0
  for i =0, n_elements(ind) -1 do begin
    result += modes[*,ind[i]]
  endfor
  return, result
end