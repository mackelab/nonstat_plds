function fromEstep = runEstep(inpn, Yn, params)

opts = optimset('Display','off');

%% (1) forward filtering
ff = forwardfiltering(inpn, Yn, params, opts);

%% (2) backward smoothing
bs = backwardsmoothing(inpn, Yn, params, opts);

%% (3) marginals / cross-covariance
fromEstep = computingmarginals(inpn, Yn, ff, bs, params);

