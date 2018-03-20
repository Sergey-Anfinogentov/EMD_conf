
;+
  ; :Description:
  ;    Computes EMD power spectrum (modal energy densities) from EMD modes
  ;
  ; :Params:
  ;    modes - EMD modes computed (in)
  ;    dt - time cadence (in)
  ; :Output:
  ;   Returns a structure with 2 fields: period and energy
  ;
  ;
  ; :Author: Sergey Anfinogentov (anfinogentov@iszf.irk.ru)
  ;-
function emd_energy_spectrum, modes, dt
  sz = size(modes)
  if not keyword_set(dt) then dt =1d
  n = sz[2]
  length = sz[1]
  energy = dblarr(n)
  period = dblarr(n)
  for i =0, n-1 do begin
    n_ext=float((size(extrema(modes[*,i])))[1])
    period[i]=2d*length / n_ext
    energy[i] = stddev(modes[*,i])^2
  endfor
  return,{period:period*dt,energy:energy}
end