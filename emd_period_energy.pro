function emd_period_energy, s
  energy = total(s^2)
  foo = pwf_gws(s,2,n_elements(s), per = per, counts = 150,/norm,/nopad)
  bar = maxi(foo, ind)
  period = interpolate(per,ind, cubic = -0.5)
  ;nex = extrema(s)
  ;period_ex = 2d*n_elements(s)/double(nex)
  return,{period:period, energy:energy};,period_ex : period_ex}
end