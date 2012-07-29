function[r] = cirpath(t,a,b,s,r0)
% CIRPATH   Simulate Cox-Ingersoll-Ross process 
% INPUTS  : t     - observation times, n*1 vector
%           a,b,s - process parameters, positive scalars
%           r0    - starting value, cir(t(1))
% OUTPUTS : r     - realized values, n*1 vector, r(i) = cir(t(i))
% NOTES   : CIR process r(t) is defined by dr = a(b-r)dt + s*sqrt(r)*dW, 
%           where W(t) is  standard Brownian motion. This implementation
%           follows Glasserman (2004, p. 124), but calls NCX2RND to draw
%           random values from the non-central chi-squared distribution.
% EXAMPLE : r = cirpath(0:.1:1,.2,.05,.1,.04)
% AUTHOR  : Dimitri Shvorob, dimitri.shvorob@vanderbilt.edu, 10/1/07
if ~isvector(t)
   error('Input argument "t" must be a vector')
end  
if any(t < 0)
   error('Input vector "t" must contain non-negative values')
end
dt = diff(t(:));
if any(dt < 0)
   error('Input vector "t" must contain increasing values')
end
par = {'a','b','s','r0'};
for i = 1:4
    y = eval(par{i}); 
    if ~(isscalar(y) && y > 0)
        error(['Input argument "' par{i} '" must be a positive scalar']) 
    end
end
n = length(t);
r = [r0; nan*dt];
v = s^2;
d = 4*a*b/v;
e = exp(-a*dt);
c = v.*(1-e)/(4*a);
for i = 1:(n-1)
    l = r(i)*e(i)/c(i); 
    r(i+1) = c(i)*ncx2rnd(d,l);
end   

