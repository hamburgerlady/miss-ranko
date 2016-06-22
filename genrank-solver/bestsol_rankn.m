function sol = bestsol_rankn(X,goodsols,params)

bestnorm = inf;
for iii = 1:length(goodsols),
    if goodsols{iii}.resnorm<bestnorm,
        bestnorm = goodsols{iii}.resnorm;
        sol = goodsols{iii};
    end
end
for kkk = 1:params.bundleiter,
    sol=bundlerankn_onestep(X,sol,params);
    sol = updateW_rankn(X,sol,params);
end

sol = sortsol(sol);
