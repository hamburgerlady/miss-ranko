load book.mat
[N,M]=size(X);

params.rk = 3;
params.inlierbnd = 0.01;
params.initnn = 1;
params.finalnormbnd = 5;
params.gksigge = 3;
params.gksize = 20;
params.minhalfn = 2;
params.cutty = 0.75;
params.nrinliersbnd = params.rk +1;
params.bundleiter = 2;
params.robust = 0;
params.maxiter = 500;
params.maxstatic = 3;
params.finN = N;
params.finM = M;
params.glueraniter = 0;
%%

tic;
sol = mr_solver_rankn(X,W,params);
t1 = toc;
disp(['Final norm: ' num2str(sol.resnorm)])
disp(['Time: ' num2str(t1) 's'])


    
    