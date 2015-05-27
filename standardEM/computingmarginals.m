function fromEstep = computingmarginals(inpn, Yn, ff, bs, params)

%% unpack msgs

mufwrd = ff.mufwrd;
sigfwrd = ff.sigfwrd;
sigstar = ff.sigstar;
sigstar0 = ff.sigstar0;

mubwrd0 = bs.mubwrd0;
inv_sigbwrd0 = bs.inv_sigbwrd0;
mubwrd = bs.mubwrd;
inv_sigbwrd = bs.inv_sigbwrd;

%% unpack params
A = params.A;
B = params.B;
C = params.C;
d = params.d;
x0 = params.x0;
V0 = params.V0;

%% computing marignals and cross-covariances
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[k T] = size(mufwrd);
pinp = size(inpn,1);
p = size(Yn,1);
mumarg = zeros(k, T);
inv_sigmarg = zeros(k, k, T);
crsscov = zeros(k, k, T-1);

% second deriv of log-likelihood
maxexp = 20;
% W = @(a) C'*diag(exp(min(C*a+d, maxexp*ones(p,1))))*C;
W = @(a) C'*diag(exp(C*a+d))*C;


for t=1:T
    
    tmp = inv(sigfwrd(:,:,t)) + inv_sigbwrd(:,:,t);
%     diagtmp = diag(tmp);
%     diagtmp(diagtmp==0) = 1./maxexp;
%     tmp(logical(eye(k))) = diagtmp;
    
    inv_sigmarg(:,:,t) = tmp;
    
    mumarg(:,t) = inv_sigmarg(:,:,t)\(sigfwrd(:,:,t)\mufwrd(:,t) + inv_sigbwrd(:,:,t)*mubwrd(:,t));
end

inv_sigmarg0 = inv(V0) + inv_sigbwrd0;


% to make sure inv_sigmarg0 is always dividable (to avoid getting NaN in
% mumarg0)

tmp = inv_sigmarg0;

% diagtmp = diag(tmp);
% diagtmp(diagtmp==0) = 1./maxexp;
% tmp(logical(eye(k))) = diagtmp;

inv_sigmarg0 = tmp;

mumarg0 = inv_sigmarg0\(inv_sigbwrd0*mubwrd0 + V0\x0);

%% cross-covaraince

% t=0,1
I = eye(k);
RRterm = W(mumarg(:,1)) + inv_sigbwrd(:,:,1) + I;
% crosscov0 = (sigstar0*A')/(RRterm - A*sigstar0*A');
% crosscov0 = sigstar0/((A/RRterm)\I - A'*sigstar0);
crosscov0 = (inv(sigstar0) - (A'/RRterm)*A)\(A'/RRterm); 
% crosscov0 = (sigstar0*A')/(W(mumarg(:,1)) - A*sigstar0*A');

% crossvar

for t=1:T-1
    %     t
    
    RRterm = W(mumarg(:,t+1)) + inv_sigbwrd(:,:,t+1) + I;
    crsscov(:, :, t) = inv(inv(sigstar(:,:,t)) - (A'/RRterm)*A)*(A'/RRterm);
    %     crsscov(:, :, t) = sigstar(:,:,t)/((A/RRterm)\I - A'*sigstar(:,:,t));
    %     crsscov(:, :, t) = (sigstar(:,:,t)*A')/( RRterm - A*sigstar(:,:,t)*A' );
    %     crsscov(:, :, t) = (sigstar(:,:,t)*A')/( W(mumarg(:,t+1)) - A*sigstar(:,:,t)*A' );
end


%% computing sufficient statistics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

WA = zeros(k, k, T);
SA = zeros(k, k, T);
GA = zeros(k, pinp, T);
Mtil = zeros(pinp, k, T);
WC = zeros(k, k, T);
SC = zeros(k, p, T);

expComgonehalfF = zeros(p, T);
highThrsh = 1e6;

for t = 1: T
    
    if pinp==0
        
        WC(:,:,t) =  inv(inv_sigmarg(:,:,t)) + mumarg(:,t)*mumarg(:,t)';
        SC(:,:,t) = mumarg(:,t)*Yn(:,t)';
        
        tmp = exp(C*mumarg(:,t) + 0.5*diag((C/inv_sigmarg(:,:,t))*C'));
%         tmp = exp(C*mumarg(:,t) + 0.5*diag((C/inv_sigmarg(:,:,t))*C'));
%         tmp(isinf(tmp)) = highThrsh;
%         tmp((tmp>highThrsh)) = highThrsh;
        expComgonehalfF(:,t) = tmp;
        
        if t ==1
            WA(:,:,t) = inv(inv_sigmarg0) + mumarg0*mumarg0';
            SA(:,:,t) = crosscov0 + mumarg0*mumarg(:,t)';
        else
            WA(:,:,t) = inv(inv_sigmarg(:,:,t-1)) + mumarg(:,t-1)*mumarg(:,t-1)';
            SA(:,:,t) = crsscov(:,:,t-1) + mumarg(:,t-1)*mumarg(:,t)';
        end
        
        suffstat.WA = sum(WA,3);
        suffstat.SA = sum(SA,3);
        suffstat.WC = sum(WC,3);
        suffstat.SC = sum(SC,3);
        suffstat.expComgonehalfFsum = sum(expComgonehalfF,2);
        
        
    else
        
        Mtil(:,:,t) = inpn(:,t)*mumarg(:,t)';
        WC(:,:,t) =  inv(inv_sigmarg(:,:,t)) + mumarg(:,t)*mumarg(:,t)';
        SC(:,:,t) = mumarg(:,t)*Yn(:,t)';
        
        tmp = exp(C*mumarg(:,t) + 0.5*diag((C/inv_sigmarg(:,:,t))*C'));
%         tmp(isinf(tmp)) = highThrsh;
%         tmp((tmp>highThrsh)) = highThrsh;
        expComgonehalfF(:,t) = tmp;
        
        if t ==1
            WA(:,:,t) = inv(inv_sigmarg0) + mumarg0*mumarg0';
            SA(:,:,t) = crosscov0 + mumarg0*mumarg(:,t)';
            GA(:,:,t) = mumarg0*inpn(:,t)';
        else
            WA(:,:,t) = inv(inv_sigmarg(:,:,t-1)) + mumarg(:,t-1)*mumarg(:,t-1)';
            SA(:,:,t) = crsscov(:,:,t-1) + mumarg(:,t-1)*mumarg(:,t)';
            GA(:,:,t) = mumarg(:,t-1)*inpn(:,t)';
        end
        
        suffstat.WA = sum(WA,3);
        suffstat.SA = sum(SA,3);
        suffstat.GA = sum(GA,3);
        suffstat.Mtil = sum(Mtil,3);
        suffstat.WC = sum(WC,3);
        suffstat.SC = sum(SC,3);
        suffstat.expComgonehalfFsum = sum(expComgonehalfF,2);
    end
    
end

%% return these:

fromEstep.mumarg = mumarg;
fromEstep.inv_sigmarg = inv_sigmarg;
fromEstep.mumarg0 = mumarg0;
fromEstep.inv_sigmarg0 = inv_sigmarg0;
fromEstep.crosscov0 = crosscov0;
fromEstep.crosscov = crsscov;
fromEstep.suffstat = suffstat;

