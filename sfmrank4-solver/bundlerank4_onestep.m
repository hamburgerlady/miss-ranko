function [sol,res] = bundlerank4_onestep(X,sol,params)
if nargin<3,
    robust = 0;
    vvin = 0;
else
    robust = params.robust;
    vvin = params.inlierbnd^2;
end

lille = 1e-7;
N0 = size(X,1);
sol.V(end,:)=1;
id = find(sol.indyi==N0);
sol.indyi(id)=[];
sol.indyi=[N0 sol.indyi];
sol.U(id,:)=[];
sol.U = [0 0 0 1;sol.U];

U0 = sol.U;
V0 = sol.V;
[N,rk] = size(U0);
M = size(V0,2);
sol.Wloc(id,:)=[];
sol.Wloc = [ones(1,M);sol.Wloc];




indyi = sol.indyi;
indyj = sol.indyj;
res0 = U0*V0-X(indyi,indyj);
if robust
    res0 = res0(sol.Wloc==1);
else
    res0 = res0(sol.Wloc~=0);
end

resgrad = calcresgrad_rank4(sol,robust);

gids = sum(resgrad~=0)>rk;
x_gn = zeros(rk*(N-rk)+(rk-1)*M,1);
x_gn(gids) = -resgrad(:,gids)\res0;



U_gn = [zeros(rk,rk);reshape(x_gn(1:rk*(N-rk)),rk,(N-rk))'];
%V_gn = reshape(x_gn(rk*(N-rk)+1:end),rk,M);
V_gn = [reshape(x_gn(rk*(N-rk)+1:end),rk-1,M);zeros(1,M)];

sol.U = U0 + U_gn;
sol.V = V0 + V_gn;

res = sol.U*sol.V-X(indyi,indyj);
if robust
    res = res(sol.Wloc==1);
else
    res = res(sol.Wloc~=0);
end

while norm(res)>(lille + norm(res0)),
    U_gn = U_gn/10;
    V_gn = V_gn/10;
    sol.U = U0 + U_gn;
    sol.V = V0 + V_gn;

    res = sol.U*sol.V-X(indyi,indyj);
    if robust
        res = res(sol.Wloc==1);
    else
        res = res(sol.Wloc~=0);
    end
end

sol.resnorm = sqrt(sum(res.^2)+(sum(sol.Wloc(:)~=0)-numel(res))*vvin);

