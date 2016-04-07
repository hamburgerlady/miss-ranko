function [P,Y,H,sol] = calibrate_rank4sol(sol,K)

U = sol.U;
V = sol.V;
nc = (length(U)-1)/2;
Ki = inv(K);
% Uh = inv(K)UH 
%|Uh(1,1:3)| = 1 
%|Uh(2,1:3)| = 1 
% Uh(1,1:3)Uh(2,1:3)^T = 0

[~,ordoi]=sort(sol.indyi);
[~,ordoj]=sort(sol.indyj);

U = U(ordoi,:);
V = V(:,ordoj);

u1 = zeros(nc,3);
u2 = zeros(nc,3);
for iii = 1:nc,
    tmp = Ki*[U(2*iii-1:2*iii,1:3);0 0 0];
    u1(iii,:)=tmp(1,:);
    u2(iii,:)=tmp(2,:);
end


boffa = [];


for iii = 1:nc,
    p11 = u1(iii,:)'*u1(iii,:);
    p12 = u1(iii,:)'*u2(iii,:);
    p22 = u2(iii,:)'*u2(iii,:);
    
    boffa = [boffa;p11(:)' -1;p22(:)' -1;p12(:)' 0];
end
CC = [1 0 0 0 0 0 0 0 0 0;0 1 0 1 0 0 0 0 0 0;0 0 1 0 0 0 1 0 0 0;0 0 0 0 1 0 0 0 0 0;0 0 0 0 0 1 0 1 0 0;0 0 0 0 0 0 0 0 1 0;0 0 0 0 0 0 0 0 0 1]';
boffa = boffa*CC;


[~,~,v]=svd(boffa);
ss = v(:,end);
ss = ss/ss(end);
HH = [ss(1) ss(2) ss(3);ss(2) ss(4) ss(5);ss(3) ss(5) ss(6)];
mm= min(eig(HH));
lille = 1e-8;
if mm<0,
    HH=HH+eye(size(HH))*(-mm+lille);
end
H = chol(HH);
H = [H' [0 0 0]';0 0 0 1];

P = cell(1,nc);
for iii=1:nc,
    tmp = Ki*[U(2*iii-1:2*iii,:);0 0 0 1]*H;
    tmp(3,1:3)=cross(tmp(1,1:3),tmp(2,1:3));
    [u,~,v]=svd(tmp(:,1:3));
    tmp(:,1:3)=u*v';
    P{iii}=tmp;
    
    
end
Y = H\V;

sol.indyi = sol.indyi(ordoi);
sol.indyj = sol.indyj(ordoj);

sol.U = U*H;
sol.V = Y;
sol.Wloc = sol.Wloc(ordoi,ordoj);






