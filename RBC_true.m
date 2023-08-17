clear , clf , clc
global k

ks = 8:-0.01:1; %range of ks to solve for
Ra=7000; %start guess for k = 8 
delta_Ra = 0.001 ; %step size for derivative

tolerance = 0.001 ; %tolerance for the absolute error
Result = [ ] ; %vector for results

for k = ks
dR = 1 ; %dummy value

while abs(dR) > tolerance
f = detA(Ra);
df = (detA(Ra+delta_Ra)-detA(Ra-delta_Ra) ) / ( 2*delta_Ra) ;
dR = - f /df ;
Ra = Ra + dR;
end
Result =[Result Ra];
end
plot(Result, ks)




function ret = detA(Ra)
global k
q0= 1i*k*(-1+(Ra/k^4) ^ ( 1 / 3 ) ) ^ ( 1/2 ) ;
q1= k*(1+(Ra/k^4) ^ ( 1 / 3 ) *(1/2+1i * sqrt ( 3 ) / 2 ) ) ^ ( 1 / 2 ) ;
q2= k*(1+(Ra/k^4) ^ ( 1 / 3 ) *(1/2-1i * sqrt ( 3 ) / 2 ) ) ^ ( 1 / 2 ) ;


A = [1 1 1 1 1 1; exp(q0) exp(-q0) exp(q1) exp(-q1) exp(q2) exp(-q2); q0 -q0 q1 -q1 q2 -q2; q0*exp(q0) -q0*exp(-q0) q1*exp(q1) -q1*exp(-q1) q2*exp(q2) -q2*exp(-q2); (q0^2-k^2)^2 (q0^2-k^2)^2 (q1^2-k^2)^2 (q1^2-k^2)^2 (q2^2-k^2)^2 (q2^2-k^2)^2; (q0^2-k^2)^2*exp(q0) (q0^2-k^2)^2*exp(-q0) (q1^2-k^2)^2*exp(q1) (q1^2-k^2)^2*exp(-q1) (q2^2-k^2)^2*exp(q2) (q2^2-k^2)^2*exp(-q2)];

ret=det(A) ;
end
