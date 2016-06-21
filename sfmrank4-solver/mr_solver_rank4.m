function sol = mr_solver_rank4(X,W,params)

[params.N,params.M] = size(X);
kul = 1;
iters = 0;
goodsols = {};
cursols = initsols_rank4(X,W,params);

while kul,
    iters = iters+1;
    [cursols,goodsols]=extendsols_rank4(X,W,cursols,goodsols,params);
    [cursols,goodsols]=gluesols_rank4(X,W,cursols,goodsols,params);

    
    NN = 0;
    for iii = 1:length(cursols),
        NN = max(length(cursols{iii}.indyi),NN);
    end
    disp(NN)
    
    
    if isempty(goodsols) && isempty(cursols),
        cursols = initsols_rank4(X,W,params);
    end
    kul = (iters<=params.maxiter) && ~isempty(cursols);
end

sol = bestsol_rank4(X,goodsols,params);

        
 