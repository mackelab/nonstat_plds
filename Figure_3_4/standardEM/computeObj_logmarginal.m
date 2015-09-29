function obj_f = computeObj_logmarginal(Yn, mumarg, inv_sigmarg, mumarg0, inv_sigmarg0, Anew, Bnew, Cnew, dnew, x0new, V0new)

T = size(mumarg, 2);


% compute likelihood term/prior/initprior
likeli_trm = zeros(T, 1);
prior_trm = zeros(T, 1);
initprior_trm = zeros(T,1);

highThrsh = 1e6;

AA = Anew'*Anew;
BB = Bnew'*Bnew; 
AB = Anew'*Bnew; 

for i=1:T
    
    %% (1) likeli_trm
    CUpsilon = Cnew/inv_sigmarg(:,:,t);
    f = diag(CUpsilon*Cnew');
    g = Cnew*mumarg(:,t) + 0.5*f + dnew;
    
    tmp = exp(g);
    tmp(isinf(tmp)) = highThrsh;
    tmp(tmp>highThrsh) = highThrsh;
    
    likeli_trm(i) = Yn(:,i)'*(Cnew*mumarg(:,i)+ dnew) - sum(tmp); 
    
    %% (2) prior_trm
    if i==1
        mumarg(:,i)'*mumarg(:,i)
    
    
    
end

