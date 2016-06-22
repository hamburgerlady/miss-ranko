function sol = updateW_rankn(X,sol,params)

res = sol.U*sol.V-X(sol.indyi,sol.indyj);
sol.Wloc(abs(res)<params.inlierbnd & sol.Wloc==-1)=1;
