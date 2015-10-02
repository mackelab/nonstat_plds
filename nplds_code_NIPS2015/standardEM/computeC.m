function Cnew = computeC(SC, mumarg, inv_sigmarg, C, d, Yn, opts1)

[k p] = size(SC); 
% size(C) is p by k

prs0 = C(:); % vectorized C
fun = @(prs)(Lossfunc(prs, SC, mumarg, inv_sigmarg, d, Yn, p, k));
DerivCheck(fun,prs0+.1);

prs = fminunc(fun, prs0, opts1);
Cnew = reshape(prs, p, k);

% ===========================================================
function [l, dl] = Lossfunc(prs, SC, mumarg, inv_sigmarg, d, Yn, p, k)
% Computes SC - sum( exp(g)*w_t' + diag(g)*C*Upsilon_t 

C = reshape(prs, p, k); 
T = size(mumarg,2);

% that we want to maximize!
% obj = sum_t( y_t'*(C w_t + d) - 1'*sum(exp(g_t)) ); 

obj_frstTrm = zeros(T,1); 
obj_scndTrm = zeros(T,1);

% d obj / dC 
% = SC' - sum( exp(g_t)* w_t') - sum(diag(exp(g_t))*C*Upsilon_t)

der_scndTrm = zeros(p, k, T);
der_thrdTrm = zeros(p, k, T);
% 
% highThrsh = 1e6;

for t=1:T
    CUpsilon = C/inv_sigmarg(:,:,t); 
    f = diag(CUpsilon*C');
    g = C*mumarg(:,t) + 0.5*f + d;
    
%     g = zeros(p,1);
%     for s=1:p
%         cs = C(s,:)';
%         ds = d(s);
%         g(s) = cs'*mumarg(:,t) + 0.5*(cs'/inv_sigmarg(:,:,t))*cs + ds;
%     end

    tmp = exp(g); 
%     tmp(isinf(tmp)) = highThrsh; 
%     tmp(tmp>highThrsh) = highThrsh; 
    
    obj_frstTrm(t) = Yn(:,t)'*(C*mumarg(:,t)+d);
    obj_scndTrm(t) = sum(tmp);
    
    der_scndTrm(:,:,t) = tmp*mumarg(:,t)';
    der_thrdTrm(:,:,t) = diag(tmp)*CUpsilon;
end


obj = sum(obj_frstTrm - obj_scndTrm);
l = - obj; 

d_obj = SC' - sum(der_scndTrm,3) - sum(der_thrdTrm, 3); 
dl = - d_obj(:);
