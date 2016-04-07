function sol = bestsol_rank4(X,goodsols,params)

bestnorm = inf;
for iii = 1:length(goodsols),
    if goodsols{iii}.resnorm<bestnorm,
        bestnorm = goodsols{iii}.resnorm;
        sol = goodsols{iii};
    end
end
for kkk = 1:params.bundleiter,
    sol=bundlerank4_onestep(X,sol,params);
    sol = updateW_rank4(X,sol,params);
end

sol = sortsol(sol);
