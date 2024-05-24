function [t,y] = rk4_2mdof(dydt,ft,dt,f,y0,M,K,C)

%===========================================================================================================================================================
%***********************************************************************************************************************************************************
%
%                                      --------------------------------------------------
%                                                          PREAMBLE
%                                      --------------------------------------------------
%
%[t,y] = rk4(dydt,ft,dt,f,y0)
%
%Adarsh S, Ph.D. Candidate, IIT Kanpur
%
%Function description:
%-------------------
%This function implements the 4th order Runge-Kutta time integration scheme with fixed time steps
%
%Input arguments:
%---------------
%1)dydt: Function handle that takes a given state vector, x, and input vector, f, and returns the values of dy/dt 
% dydt(x,f)
%2)ft: Time vector of input forces  (row vector)
%3)dt: Sampling time of ft
%4)f: Input forces arranged in rows
%5)y0: Column vector of state vector at initial time
%----------------
%
%Outputs arguments:
%-----------------
%1)y: The evaluated states arranged in rows
%2)t: Time vector for y
%-----------------
%
%
%Ex:  [t,y] = rk4(@(y,fk)dydt_sdof(y,fk,m,k,c),ft,1/200,f,[0;0]) ; 
%
%
%
%                                      -----------------------------------------------
%                            *********|| All rights reserved; Adarsh S; November, 2017 || *********
%                                      -----------------------------------------------
%
%
%***********************************************************************************************************************************************************
%===========================================================================================================================================================

yn = length(y0) ;
t = ft ;
[mf,~] = size(f) ;

y = zeros(yn,length(t)) ;
y(:,1) = y0 ;

for i = 2:1:length(t)

dy1 = dydt(y(:,i-1),f(:,i-1),M(:,:,i),K(:,:,i),C(:,:,i))*dt ;
f2 = zeros(mf,1) ;

for j = 1:1:mf

f2(j,1) = interp1(t,f(j,:),( t(i-1) + dt/2 ));


end

dy2 = dydt(y(:,i-1)+(dy1/2),f2,M(:,:,i),K(:,:,i),C(:,:,i))*dt ;

dy3 = dydt(y(:,i-1)+(dy2/2),f2,M(:,:,i),K(:,:,i),C(:,:,i))*dt ;

dy4 = dydt(y(:,i-1)+(dy3/2),f(:,i),M(:,:,i),K(:,:,i),C(:,:,i))*dt ;

y(:,i) = y(:,i-1) + (dy1 + 2*dy2 + 2*dy3 + dy4)/6 ;

end





end