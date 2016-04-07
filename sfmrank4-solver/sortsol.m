function sol = sortsol(sol)

[~,ordoi]=sort(sol.indyi);
[~,ordoj]=sort(sol.indyj);

sol.indyi = sol.indyi(ordoi);
sol.indyj = sol.indyj(ordoj);
sol.U = sol.U(ordoi,:);
sol.V = sol.V(:,ordoj);
sol.Wloc = sol.Wloc(ordoi,ordoj);
