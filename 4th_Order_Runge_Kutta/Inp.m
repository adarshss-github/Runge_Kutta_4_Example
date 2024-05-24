[t,y] = rk4(@(y,fk)dydt_sdof(y,fk,m,k,c),ft,1/200,f,2,[0;0]) ;
[tup,yup] = rk4(@(y,fk)dydt_sdof(y,fk,m,k,c),ftup,1/2000,fup,2,[0;0]) ;


[t,y] = rk4(@(y,fk)dydt_sdof(y,fk,m,k,c),ft,1/200,f,[0;0]) ; %Latest


[t1,y1] = rk4(@(y,fk)dydt_sdof(y,fk,m,k,c),ft,1/200,f,[0;0]) ; %Latest

[t2,y2] = rk4(@(y,fk)dydt_sdof2(y,fk,m,k,c),ft,1/200,f,[0;0]) ;


