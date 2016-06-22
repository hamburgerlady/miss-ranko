function sol = mr_solver_rankn(X,W,params)

[params.N,params.M] = size(X);
kul = 1;
iters = 0;
goodsols = {};
cursols = initsols_rankn(X,W,params);

while kul,
    iters = iters+1;
    [cursols,goodsols]=extendsols_rankn(X,W,cursols,goodsols,params);
    [cursols,goodsols]=gluesols_rankn(X,W,cursols,goodsols,params);

    
    NN = 0;
    MM = 0;
    for iii = 1:length(cursols),
        NN = max(length(cursols{iii}.indyi),NN);
        MM = max(length(cursols{iii}.indyj),MM);
        
    end
    disp([iters NN MM length(cursols) length(goodsols)])
    
    
    if isempty(goodsols) && isempty(cursols),
        cursols = initsols_rankn(X,W,params);
    end
    kul = (iters<=params.maxiter) && ~isempty(cursols);
end

sol = bestsol_rankn(X,goodsols,params);

        
 