function W = computeHessian(mufwrd, covC, covd, C, d, maxexp)


[quad, lin] = xSigcx(mufwrd, covC);
p = size(C,1);
g = exp(min(d + 0.5*diag(covd) + C*mufwrd + 0.5*quad, maxexp*ones(p,1)));
frstTrm = bsxfun(@times, covC, reshape(g,[1 1 p]));
scndTrm = (C'+ lin)*diag(g)*(C'+ lin)';

W = sum(frstTrm,3) + scndTrm; 
