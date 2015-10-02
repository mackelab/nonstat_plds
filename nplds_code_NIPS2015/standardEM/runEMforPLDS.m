function datastruct = runEMforPLDS(xyzinpn, k, pinp, p, ddotU, trueparams)

% datastruct has the following:

% mumarg: posterior mean of x
% inv_sigmarg: inverse posterior covariance of x
% crosscov0, crosscov
% estimates of A, B, C, x0, V0

%% unpack xyzinpn

inpn = xyzinpn.inp;
Yn = xyzinpn.y;

%% initial params

initparams = generate_params(k, pinp, p);

%% Iterations until convergence.

nIter = 0;
maxIter = 100;
T = size(Yn,2);

% to store mumarg, inv_sigmarg, crosscov, A, B, C, and mse values
fromEstepcell = cell(maxIter,1);
fromMstepcell = cell(maxIter,1);
msevec = zeros(maxIter, 1); 

% lam_tst from true x for one-step-ahead prediction
nsamps = 1000;
truez = [xyzinpn.z xyzinpn.z_tst];
lam_tst = exp(truez(:,2:end));

% mse for stopping criterion
mse = @(a) norm(a-lam_tst);

while (nIter<=maxIter)
    %% increase iteration number
    nIter = nIter + 1;
    
    [nIter maxIter]
    
    %% E-step
    
    fromEstep = runEstep(inpn, Yn, initparams);
    
    %% M-step
    
    fromMstep = runMstep(Yn, fromEstep, initparams, ddotU, trueparams);
    
    %% compute one-step head prediction
    
    estlam_tst = generate_one_step_ahead_prediction_tst(fromEstep, fromMstep, xyzinpn, nsamps);
    
    msevec(nIter) = mse(estlam_tst);
    
%     figure(1);
%     subplot(515); plot(msevec(1:nIter), 'o-'); set(gca, 'xlim', [0 maxIter]); pause(0.1);
    
%     %% plotting
%     
%     subplot(512); plot(1:T, xyzinpn.x', 'k', 1:T, fromEstep.mumarg', 'r'); title('true latent states');
%     
%     figure(2);
%     subplot(321); hinton(trueparams.A,['true A  ' num2str(max(max(trueparams.A)),'%.4f')],'standard');
%     subplot(322); hinton(fromMstep.A,['estimate A  ' num2str(max(max(fromMstep.A)),'%.4f')],'standard');
%     subplot(323); hinton(trueparams.B,['true B  ' num2str(max(max(trueparams.B)),'%.4f')],'standard');
%     subplot(324); hinton(fromMstep.B,['estimate B  ' num2str(max(max(fromMstep.B)),'%.4f')],'standard');
%     subplot(325); hinton(trueparams.C,['true C  ' num2str(max(max(trueparams.C)),'%.4f')],'standard');
%     subplot(326); hinton(fromMstep.C,['estimate C  ' num2str(max(max(fromMstep.C)),'%.4f')],'standard');
%    
    %% store these:
    
    fromEstepcell{nIter} = fromEstep;
    fromMstepcell{nIter} = fromMstep;
    
    %% update params
    
    initparams = fromMstep;
    
end


%% choose the best results

[~, minidx] = min(msevec);

datastruct.Estep = fromEstepcell{minidx};
datastruct.Mstep = fromMstepcell{minidx};

