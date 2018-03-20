Function EXTREMA,  Data, MAXIMA=maxima, MINIMA=minima
  n = n_elements(data)
  minima = lonarr(n)
  maxima = lonarr(n)
  n_ma = 0
  n_mi = 0

  
  for i = 1l, n-2 do begin
    if (data[i] ge data[i-1]) and (data[i] gt data[i+1]) then begin
      maxima[n_ma] = i
      n_ma+=1
    endif
    if (data[i] le data[i-1]) and (data[i] lt data[i+1]) then begin
      minima[n_mi] = i
      n_mi+=1
    endif  
  endfor
 
  minima = minima[0:n_mi-1]
  maxima = maxima[0:n_ma-1]
  
  
  if (data[0] gt data[maxima[0]]) or (data[0] lt data[minima[0]]) then begin
    if data[0] gt data[1] then begin
      maxima = [0l,maxima]
      n_ma+=1
    endif
    if data[0] lt data[1] then begin
      minima = [0l,minima]
      n_mi+=1
    endif
  endif
  
  if (data[n-1] gt data[maxima[n_ma-1]]) or (data[n-1] lt data[minima[n_mi-1]]) then begin
    if data[n-1] gt data[n-2] then begin
      maxima = [maxima,n-1l]
    endif
    if data[n-1] lt data[n-2] then begin
      minima = [minima,n-1l]
    endif
  endif
  
  
  if n_mi eq 0 then foo = min(data,minima)
  if n_ma eq 0 then foo = max(data,maxima)
  result = [maxima,minima]
  return, result[sort(result)]
end