function ff = forwardfiltering(inpn, Yn, params, opts)

%% unpack params
A = params.A;
B = params.B;
C = params.C;
d = params.d;
x0 = params.x0;
V0 = params.V0;

AA = A'*A;

%% forward filtering

T = size(Yn,2);
k = length(x0);
pinp = size(B,2);

if pinp>0
    AB = A'*B;
end

% matrices for storing mean/cov of forward messages
sigstar = zeros(k, k, T);

mutilde = zeros(k, T);
sigtilde = zeros(k, k, T);

mufwrd = zeros(k, T);
sigfwrd = zeros(k, k, T);

mufwrd(:,1) = x0; 
sigfwrd(:,:,1) = V0;

I = eye(k);

for t=1:T

    if t==1
        
        sigstar0 = inv(inv(V0) + AA);
        
        sigtilde(:,:,t) = inv(I - A*sigstar0*A');
        
        if pinp>0
            mutilde(:,t) = sigtilde(:,:,t)*(B*inpn(:,t) + A*sigstar0*inv(V0)*x0 - A*sigstar0*AB*inpn(:,t));
        else
            mutilde(:,t) = sigtilde(:,:,t)*(A*sigstar0*inv(V0)*x0);
        end
        
        [mufwrd(:,t), sigfwrd(:,:,t)] = findModeCov(mutilde(:,t), sigtilde(:,:,t), C, d, Yn(:,t), opts);
        
    else
        
        
        %(1) compute sigstar
        sigstar(:,:,t-1) = inv(inv(sigfwrd(:,:,t-1)) + AA);
        
        %(2) compute mutilde and sigtilde
%         sigtilde(:,:,t) = I + A*sigfwrd(:,:,t-1)*A';
%         mutilde(:,t) = A*mufwrd(:,t-1) + B*inpn(:,t); 
        sigtilde(:,:,t) = inv(I - A*sigstar(:,:,t-1)*A');
        if pinp>0
            mutilde(:,t) = sigtilde(:,:,t)*(B*inpn(:,t) + A*(sigstar(:,:,t-1)/sigfwrd(:,:,t-1))*mufwrd(:,t-1) - A*sigstar(:,:,t-1)*AB*inpn(:,t));
        else
            mutilde(:,t) = sigtilde(:,:,t)*(A*(sigstar(:,:,t-1)/sigfwrd(:,:,t-1))*mufwrd(:,t-1));
        end
        
        [mufwrd(:,t), sigfwrd(:,:,t)] = findModeCov(mutilde(:,t), sigtilde(:,:,t), C, d, Yn(:,t), opts);
        
    end

end

%% return these:

ff.mufwrd = mufwrd;
ff.sigfwrd = sigfwrd;
ff.sigstar = sigstar;
ff.sigstar0 = sigstar0;
