function dydt = dydt_sdof(y,f,m,k,c)



dydt = [ y(2); (f - ( k*y(1) ) - ( c*y(2) ) )/m] ; 



end