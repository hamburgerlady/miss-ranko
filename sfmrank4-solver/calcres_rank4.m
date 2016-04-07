function res = calcres_rank4(X,sol,robust)
if nargin<3,
    robust = 0;
end

res = sol.U*sol.V-X(sol.indyi,sol.indyj);
if robust
    res = res(sol.Wloc==1);
else
    res = res(sol.Wloc~=0);
end