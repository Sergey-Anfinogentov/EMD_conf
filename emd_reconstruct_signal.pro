;+
; :Description:
;    Recover "clean" signal from EMD modes as a superposition
;    of significant modes and the trend (last mode)
;
; :Params:
;    modes        - EMD modes calculated with EMD function
;    conf_period  - period returned by EMD_NOISE_CONF function
;    conf_up      - upper confidentce level computed by EMD_NOISE_CONF function
;
;
;
; :Author: Sergey Anfinogentov (anfinogentov@iszf.irk.ru)
;-
function EMD_reconstruct_signal,modes, conf_period, conf_up
  sp = emd_energy_spectrum(modes)
  nt = (size(modes))[1]
  mode_period = sp.period[1:-2]
  mode_energy = sp.energy[1:-2]
  conf_mod = interpol(conf_up,conf_period,mode_period)
  ind = where(mode_energy gt conf_mod and mode_period le 0.7*nt)
  emd_clean = dblarr(n_elements(modes[*,0]))
  if ind[0] ne -1 then begin
    for i =0, n_elements(ind) -1 do begin
      emd_clean += modes[*,ind[i]+1]
    endfor
  endif
  return, emd_clean
end