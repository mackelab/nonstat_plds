function params = generate_params(k, pinp, p)

% prior on x
params.V0 = eye(k);
params.x0 = zeros(k,1);

% latent noise covariance
params.Q = eye(k);

% A
Anrm = 0.7;
A = randn(k);
params.A = A./max(abs(eigs(A)))*Anrm;

% B
if pinp==0
    params.B = [];
else
    Bnrm = 1.2;
    B = randn(k, pinp);
    params.B = Bnrm*B/norm(B);
end

% C
Cnrm = 1; 
C = randn(p,k);
params.C = C/norm(C)*Cnrm;
% params.C = Cnrm*[eye(k)/sqrt(k); zeros(p-k, k)];

% d
params.d = -2*ones(p,1);
