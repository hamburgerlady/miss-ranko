function [sol,ok] = gluei_rank4_ransac(X,W,sol1,sol2,params)

rk = params.rk;
inlierbnd = params.inlierbnd;
nrinliersbnd = params.nrinliersbnd;
raniter = params.glueraniter;
robust = params.robust;
indyi1 = sol1.indyi;
indyj1 = sol1.indyj;
indyi2 = sol2.indyi;
indyj2 = sol2.indyj;
U1 = sol1.U;
V1 = sol1.V;
U2 = sol2.U;
V2 = sol2.V;
ok = 0;

[indyi,idsi1] = intersect(indyi1,indyi2);
[~,idsj1] = intersect(indyj1,indyj2);

N = size(X,1);

if length(indyi)>=rk,
    for kkk = 1:raniter,
        ppi = randperm(length(indyi));
        ppi = ppi(1:rk);
        if ~any(indyi(ppi)==N),
            ppi(end)=find(indyi==N);
        end
        
        ww = W(indyi(ppi),indyj2);
        gids = find(sum(ww)==rk);
        if length(gids)>=rk,
            ppj = randperm(length(gids));
            ppj = ppj(1:rk);
            xx = X(indyi(ppi),indyj2(gids(ppj)));
            Hsol = (U1(idsi1(ppi),:)\xx)/V2(:,gids(ppj));
            
            
            sol.U = [U1;U2/Hsol];
            sol.U(idsi1,:)=[];
            sol.V = [V1 Hsol*V2];
            sol.V(:,idsj1)=[];
            sol.indyi = [indyi1 indyi2];
            sol.indyi(idsi1)=[];
            sol.indyj = [indyj1 indyj2];
            sol.indyj(idsj1)=[];
            sol.ground(2)=mean(sol.indyj);
            sol.ground(1)=mean(setdiff(sol.indyi,N));
            sol.static = 0;
            sol.Wloc = W(sol.indyi,sol.indyj);
            res = abs(sol.U*sol.V-X(sol.indyi,sol.indyj));
            res(sol.Wloc==0)=inf;
            
            if all(sum(res<inlierbnd,1)>=nrinliersbnd) && all(sum(res<inlierbnd,2)>=nrinliersbnd),
                ok = 1;
            end
            if ~robust,
                ok = 1;
            end
        else
            ok = 0;
        end
        
        if ok,
            break;
        end
    end
    
    
else
    ok = 0;
end

if ~ok
    sol = -1;
end

