clear all
load dinodata
% loads data matrices X and W for the dino-seq
% also loads ground thruth (projective rec) 3D points in Y0 and 
% affine projections of ground truth projective cameras in P0

[N,M]=size(X);
params.rk = 4;
params.robust = 0;
params.inlierbnd = 0.1;
params.nrinliersbnd = params.rk;
%params.ratioinlierbnd = 0.8;
%params.initnrinliersbnd = 60;
%params.initnn = 1;
params.initnn = 3;
params.minhalfn = 4;
params.bundleiter = 2;
%params.bundleiter = 3;
params.glueraniter = 5;
params.maxiter = 5000;
params.finalnormbnd = 100;
%params.reskulfact = 1;
params.cutty = 0.75;
params.gksigge = 3;
params.gksize = 20;
params.maxbundlesize = 200000;
params.maxstatic = 1;

params.finN = N;
params.finM = M;


tic;
sol = mr_solver_rank4(X,W,params);
t1 = toc;
disp(['Final norm: ' num2str(sol.resnorm)])
disp(['Time: ' num2str(t1) 's'])


clf

V = sol.V(1:3,:);
V0 = Y0(1:3,sol.indyj);

V = V-repmat(mean(V,2),1,M);
V0 = V0-repmat(mean(V0,2),1,M);


%A = (V0(:,proids(1:3))-repmat(V0(:,proids(4)),1,3))/(V(:,proids(1:3))-repmat(V(:,proids(4)),1,3));
%b = V0(:,proids(4))-A*V(:,proids(4));
%Vpr = A*V+repmat(b,1,length(V));

A = V0/V;
Vpr = A*V;
figure(1)
hold off
plot3(Vpr(1,:),Vpr(2,:),Vpr(3,:),'.')
hold on
plot3(V0(1,:),V0(2,:),V0(3,:),'g.')





