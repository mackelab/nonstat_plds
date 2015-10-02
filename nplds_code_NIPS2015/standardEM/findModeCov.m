function [mufwrd, sigfwrd] = findModeCov(mutilde, sigtilde, C, d, y, opts)

% p = size(C,1);
% maxexp = 1e2; 

%(1) update mu
fun = @(a) C'*(y-exp(C*a+d)) - sigtilde\(a - mutilde); 
% fun = @(a) C'*(y-exp(min(C*a+d, maxexp*ones(p,1)))) - sigtilde\(a - mutilde); 
a0 = zeros(length(mutilde),1);
mufwrd = fsolve(fun, a0, opts);

%(2) update V
sigfwrd = inv(C'*diag(exp(C*mufwrd+d))*C + inv(sigtilde)); 
% sigfwrd = inv(C'*diag(exp(min(C*mufwrd+d, maxexp*ones(p,1))))*C + inv(sigtilde)); 
