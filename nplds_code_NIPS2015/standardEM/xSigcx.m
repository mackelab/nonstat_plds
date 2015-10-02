function [quad, lin] = xSigcx(x, covC)

% covC is k by k by p
% x is a column vector of length k
% this function computes the quadratic form: x'*covC_i*x for each covC_i
% also linear term: covC_i*x for each covC_i

p = size(covC,3);
k = size(x,1);
quad = zeros(p,1);
lin = zeros(k,p);

for i=1:p
    quad(i) = x'*covC(:,:,i)*x;
    lin(:,i) = covC(:,:,i)*x;
end
