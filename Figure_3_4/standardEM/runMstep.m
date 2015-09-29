function fromMstep = runMstep(Yn, fromEstep, params, ddotU, trueparams)

%% unpack sufficient statistics

suffstat = fromEstep.suffstat;

pinp = size(params.B,2);

if pinp>0
    GA = suffstat.GA;
    Mtil = suffstat.Mtil;
end

WA = suffstat.WA;
SA = suffstat.SA;
SC = suffstat.SC;

%% unpack initial parameters

A = params.A;
% B = params.B;
C = params.C;
% d = params.d;

% A = trueparams.A;
B = trueparams.B;
d = trueparams.d;
x0 = trueparams.x0;
V0 = trueparams.V0;

%% (1) update x0 and V0

% x0 = fromEstep.mumarg0;
% V0 = inv(fromEstep.inv_sigmarg0);

%% (2) update A

if pinp==0
%     A = SA'/WA;
    A = SA'/WA + 0.001*eye(size(A));
else
    A = (SA' - B*GA')/WA;
end
% to avoid numerical instatilibity
% A = A./max(abs(eigs(A))).*0.9;
maxeig = max(abs(eigs((A))));
if maxeig>1
    A = A./maxeig*0.5;
end

%% (3) update B

if pinp == 0
    B = [];
else
    B = (Mtil' - A*GA)/ddotU;
end

%% (4) update C

opts1 = optimset('Display', 'off', 'GradObj' , 'on', 'Hessian', 'off', 'maxIter', 1e3);
C = computeC(SC, fromEstep.mumarg, fromEstep.inv_sigmarg, C, d, Yn, opts1);

%% return these:

fromMstep.A = A;
fromMstep.B = B;
fromMstep.C = C;
fromMstep.d = d; 

fromMstep.x0 = x0;
fromMstep.V0 = V0;
fromMstep.Q = trueparams.Q;



