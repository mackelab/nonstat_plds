function datastruct = VBEM_PLDSnonstationary(xyzinpn, r, trueparams, Model, varargin)
% everything = VBEM_PLDSnonstationary_A(xyzinpn, r, params);

% May 08, 2014
% wrote by Mijung Park

%% essential quantities

k = size(trueparams.A(:,:,1),1);
p = size(trueparams.C,1);
d = size(trueparams.B,2);

%% initial params

% initparams = generate_params_from_prior(k, p, r);
ind_train = trueparams.ind_train;
initparams = generate_params_from_prior(k, d, p, r, Model, ind_train, trueparams.tau_init);
initparams.d = trueparams.d; 
initparams.inpn = trueparams.inpn;

initparams.m_h = trueparams.m_h_init;
initparams.tau2 = trueparams.tau_init;
initparams.sig2 = trueparams.sig_init;

if strcmp(Model, 'PLDS')
  initparams.h = zeros(size(initparams.h));
  initparams.covh = repmat(zeros(k),[1,1,r]);
end
  

%% VBEM (for now r==1)

nIter = 0;
maxIter = 3; %FITPARAM Iteration number
% T = size(Yn,2);

% to store mumarg, inv_sigmarg, crosscov, A, B, C, and mse values
fromVBEstepcell = cell(maxIter,1);
fromVBMstepcell = cell(maxIter,1);
% msevec = zeros(maxIter, 1);

% alphahist = zeros(maxIter, 1);
% betahist = zeros(maxIter, 1);
% sig2hist = zeros(maxIter, 1);
% tau2hist = zeros(maxIter, 1);

%%
while (nIter<maxIter)
    
    %% update #Iter
    
    nIter = nIter + 1;
    
    
    
    fprintf(['Iteration ' num2str(nIter) '/' num2str(maxIter) ' started ' datestr(clock, 'yy-mm-dd-HH:MM:SS') '\n']);
    
    %Check if iteration has already been done
    if nargin > 4
        %varargin{1} is output folder for intermediate results
        if exist([varargin{1} filesep 'datastruct_' Model '_iter_' num2str(nIter) '.mat'], 'file')
          load([varargin{1} filesep 'datastruct_' Model '_iter_' num2str(nIter) '.mat'], 'datastruct_cur_iter');
          fromVBEstep = datastruct_cur_iter.Estep;
          fromVBMstep = datastruct_cur_iter.Mstep;
          initparams = update_params(fromVBMstep, initparams);
          fromVBEstepcell{nIter} = fromVBEstep;
          fromVBMstepcell{nIter} = fromVBMstep;
          continue;
        end
    end
    
    
    %% (1) VB E-step:
    
%     tic;
    fromVBEstep = runVBEstep(xyzinpn, initparams, r, Model);
    fprintf(['E-step done ' datestr(clock, 'yy-mm-dd-HH:MM:SS') '\n']);
%     toc;
    
%   
        
    %% (2) VB M-step
    
    fromVBMstep = runVBMstep(xyzinpn, fromVBEstep, initparams, r, ind_train,Model);
    
    fprintf(['M-step done ' datestr(clock, 'yy-mm-dd-HH:MM:SS') '\n']);

    
    %% store these:
    
    fromVBEstepcell{nIter} = fromVBEstep;
    fromVBMstepcell{nIter} = fromVBMstep;
    
    initparams = update_params(fromVBMstep, initparams);
    
    %Save output at current iter into output folder
    if nargin > 4
        %varargin{1} is output folder for intermediate results
        datastruct_cur_iter.Estep = fromVBEstep;
        datastruct_cur_iter.Mstep = fromVBMstep;
        save([varargin{1} filesep 'datastruct_' Model '_iter_' num2str(nIter) '.mat'], 'datastruct_cur_iter');
    end
    
    
end


%% output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% choose the best results

% [~, minidx] = min(msevec);
% 
% datastruct.Estep = fromVBEstepcell{minidx};
% datastruct.Mstep = fromVBMstepcell{minidx};
% datastruct.msevec = msevec;

datastruct.Estep = fromVBEstepcell;
datastruct.Mstep = fromVBMstepcell;

save([varargin{1} filesep 'datastruct_' Model '_final.mat'], 'datastruct')

  function initparams = update_params(fromVBMstep, initparams)
    %% update params (that you want to update)
    
    initparams.h = fromVBMstep.h;
    initparams.covh = fromVBMstep.covh;
    
    initparams.C = fromVBMstep.C;

    initparams.A = fromVBMstep.A;
    initparams.AA = fromVBMstep.AA;
    initparams.covA = fromVBMstep.covA;
    initparams.B = fromVBMstep.B;
    initparams.AB = fromVBMstep.AB;
    initparams.covB = fromVBMstep.covB;
    
    initparams.alpha = fromVBMstep.alpha;
    initparams.beta = fromVBMstep.beta;
    initparams.sig2 = fromVBMstep.sig2;
    initparams.tau2 = fromVBMstep.tau2;
  end

end
