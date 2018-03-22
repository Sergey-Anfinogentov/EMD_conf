;+
   ; :Description:
   ;    Generates coloured noise
   ;
   ; :Params:
   ;    seed  - random number generator seed
   ;    alpha - Coloured noise index
   ;    n     - size of the array to generate 
   ;
   ;
   ;
   ; :Author: Sergey Anfinogentov (anfinogentov@iszf.irk.ru)
   ;-
function random_alpha,seed,alpha,n
   n2 = long(n*2l)
   r = randomn(seed,n2,/double)
   f = [dindgen(n),reverse(dindgen(n)+1d)]
   s = 1d/f^(alpha/2d)
   s[0] =0d
   s/= sqrt(total(s^2)/(n2))
   result = double(fft(s*fft(r),/inv)  )
   result = (result-mean(result))/stddev(result)
   return, result[0:n-1]

end
