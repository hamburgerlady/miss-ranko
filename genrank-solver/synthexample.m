clear all
N = 100;
M = 100;
rk = 5;
noisesig = 1e-3;
bw0 = 20;
U0 = randn(N,rk);
V0 = randn(rk,M);
X = U0*V0+randn(N,M)*noisesig;
W = ones(N,M)-triu(ones(N,M),bw0)-tril(ones(N,M),-bw0);


% setup parameters
params.rk = 5;

params.robust = 0;
params.inlierbnd = 0.01;

params.nrinliersbnd = params.rk;
params.initnn = 10;
params.minhalfn = 4;
params.bundleiter = 2;

params.glueraniter = 5;
params.extendraniter = 5;

params.maxiter = 5000;
params.finalnormbnd = 100;
params.cutty = 0.75;
params.gksigge = 3;
params.gksize = 20;
params.maxstatic = 1;
params.finN = N;
params.finM = M;

% run solver
tic;
sol = mr_solver_rankn(X,W,params);
t1 = toc;
disp(['Final norm: ' num2str(sol.resnorm)])
disp(['Time: ' num2str(t1) 's'])
