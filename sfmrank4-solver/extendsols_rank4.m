function [cursols,goodsols]=extendsols_rank4(X,W,cursols,goodsols,params)

badids = [];
goodids = [];
nn = length(cursols);
for iii = 1:nn,
    [cursols{iii},ok1] = extendj_rank4_ransac(X,W,cursols{iii},params);
    [cursols{iii},ok2] = extendi_rank4_ransac(X,W,cursols{iii},params);
    for kkk = 1:params.bundleiter,
        cursols{iii}=bundlerank4_onestep(X,cursols{iii},params);
        cursols{iii} = updateW_rank4(X,cursols{iii},params);
    end
    
    if ~ok1 && ~ok2,
        cursols{iii}.static = cursols{iii}.static + 1;
    end
    if cursols{iii}.static>params.maxstatic,
        badids = [badids iii];
    end
    if cursols{iii}.resnorm>params.finalnormbnd,
        badids = [badids iii];
    else
        
        maxN = length(cursols{iii}.indyi);
        maxM = length(cursols{iii}.indyj);
        if (maxN==params.finN) && (maxM==params.finM),
            goodids = [goodids iii];
            badids = [badids iii];
        end
    end
    
end

goodsols = [goodsols cursols(goodids)];
cursols(badids)=[];
