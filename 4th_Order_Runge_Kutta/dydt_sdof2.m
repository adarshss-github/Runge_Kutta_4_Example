function dydt = dydt_sdof2(y,f,m,k,c)



dydt = [ y(2); (f - ( c*y(2) ) - ( k*y(1) ) )/m] ; 



end