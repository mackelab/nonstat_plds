function bs = backwardsmoothing(inpn, Yn, params, opts)

%% unpack params
A = params.A;
B = params.B;
C = params.C;
d = params.d;
x0 = params.x0;
V0 = params.V0;

AA = A'*A;

%% matrices to store mean/cov of backward msgs

T = size(Yn,2);
k = size(A,1);
p = size(C,1);
pinp = size(B,2);

if pinp>0
    AB = A'*B;
end

maxexp = 1e3; 
I = eye(k);

etatilde = zeros(k, T);
inv_phitilde = zeros(k, k, T);

inv_phistar = zeros(k, k, T);

mubwrd = zeros(k, T);
inv_sigbwrd = zeros(k, k, T);

% for t = T, inv(sigbwrd) = 0 to ensure beta(x_T) = 1
t = T;
fun = @(a) C'*(Yn(:,t) - exp(C*a+d));
% fun = @(a) C'*(Yn(:,t) - exp(min(C*a+d, maxexp*ones(p,1))));
a0 = zeros(k,1);
etatilde(:,t) = fsolve(fun, a0, opts);
inv_phitilde(:,:,t) = C'*diag(exp(C*etatilde(:,t)+d))*C;

inv_phistar(:,:,t) = I + inv_phitilde(:,:,t);

% inv_sigbwrd(:,:,t-1) = AA - A'*(inv_phistar(:,:,t)\A);

tmp = AA - A'*(inv_phistar(:,:,t)\A);
% tmp = AA - A'*(inv_phistar(:,:,t)\A);
% diagtmp = diag(tmp);
% diagtmp(diagtmp==0) = 1./maxexp;
% tmp(logical(eye(k))) = diagtmp;

inv_sigbwrd(:,:,t-1) = tmp; 

if pinp >0 
    mubwrd(:,t-1) = inv_sigbwrd(:,:,t-1)\( (A'/inv_phistar(:,:,t))*(B*inpn(:,t) + inv_phitilde(:,:,t)*etatilde(:,t)) - AB*inpn(:,t) );
else
    mubwrd(:,t-1) = inv_sigbwrd(:,:,t-1)\( (A'/inv_phistar(:,:,t))*(inv_phitilde(:,:,t)*etatilde(:,t)));
end

%%

for t = T-1:-1:2

%     fun = @(a) C'*(Yn(:,t) - exp(min(C*a + d, maxexp*ones(p,1)))) - inv_sigbwrd(:,:,t)*(a - mubwrd(:,t));
    fun = @(a) C'*(Yn(:,t) - exp(C*a + d)) - inv_sigbwrd(:,:,t)*(a - mubwrd(:,t));
    a0 = etatilde(:,t+1);
    etatilde(:,t) = fsolve(fun, a0, opts);
    inv_phitilde(:,:,t) = C'*diag(exp(C*etatilde(:,t)+d))*C + inv_sigbwrd(:,:,t);

%     inv_phitilde(:,:,t) = C'*diag(exp(min(C*etatilde(:,t)+d, maxexp*ones(p,1))))*C + inv_sigbwrd(:,:,t);
    
    inv_phistar(:,:,t) = I + inv_phitilde(:,:,t);
    
    
    inv_sigbwrd(:,:,t-1) = AA - A'*(inv_phistar(:,:,t)\A);

    % to make sure inv_sigbwrd is always divisible (to avoid getting NaN in
    % mubwrd)
    
    tmp = AA - A'*(inv_phistar(:,:,t)\A);
%     diagtmp = diag(tmp);
%     diagtmp(diagtmp < 1/maxexp) = 1/maxexp;
%     tmp(logical(eye(k))) = diagtmp; 
    
    inv_sigbwrd(:,:,t-1) = tmp;
    
    if pinp>0
        mubwrd(:,t-1) = inv_sigbwrd(:,:,t-1)\( (A'/inv_phistar(:,:,t))*(B*inpn(:,t) + inv_phitilde(:,:,t)*etatilde(:,t)) - AB*inpn(:,t) );
    else
        mubwrd(:,t-1) = inv_sigbwrd(:,:,t-1)\( (A'/inv_phistar(:,:,t))*(inv_phitilde(:,:,t)*etatilde(:,t)));
    end
    
end


%%
t = 1;
fun = @(a) C'*(Yn(:,t) - exp(C*a+d)) - inv_sigbwrd(:,:,t)*(a - mubwrd(:,t));
% fun = @(a) C'*(Yn(:,t) - exp(min(C*a+d, maxexp*ones(p,1)))) - inv_sigbwrd(:,:,t)*(a - mubwrd(:,t));
a0 = etatilde(:,t+1);
etatilde(:,t) = fsolve(fun, a0, opts);
% inv_phitilde(:,:,t) = C'*diag(exp(min(C*etatilde(:,t)+d, maxexp*ones(p,1))))*C + inv_sigbwrd(:,:,t);
inv_phitilde(:,:,t) = C'*diag(exp(C*etatilde(:,t)+d))*C + inv_sigbwrd(:,:,t);

inv_phistar(:,:,t) = I+inv_phitilde(:,:,t);

% inv_sigbwrd0 = AA - A'*(inv_phistar(:,:,t)\A);
tmp = AA - A'*(inv_phistar(:,:,t)\A);
% tmp = AA - A'*(inv_phistar(:,:,t)\A);
% diagtmp = diag(tmp);
% diagtmp(diagtmp==0) = 1./maxexp;
% tmp(logical(eye(k))) = diagtmp;

inv_sigbwrd0 = tmp;

if pinp>0
    mubwrd0 = inv_sigbwrd0\( (A'/inv_phistar(:,:,t))*(B*inpn(:,t) + inv_phitilde(:,:,t)*etatilde(:,t)) - AB*inpn(:,t) );
else
    mubwrd0 = inv_sigbwrd0\( (A'/inv_phistar(:,:,t))*(inv_phitilde(:,:,t)*etatilde(:,t)));
end

%% return these

bs.mubwrd0 = mubwrd0;
bs.inv_sigbwrd0 = inv_sigbwrd0;
bs.mubwrd = mubwrd;
bs.inv_sigbwrd = inv_sigbwrd;


