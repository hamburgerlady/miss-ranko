function [sol,ok] = extendi_rank4_ransac(X,W,sol,params)


inlierbnd = params.inlierbnd;
robust = params.robust;

ok = 0;
rk = params.rk;
U = sol.U;
V = sol.V;

indyi = sol.indyi;
indyj = sol.indyj;


if robust,
    nrinliersbnd = params.nrinliersbnd;
    raniter = params.extendraniter;

    Wj = W(:,indyj);
    extindyi = find(sum(Wj,2)>=nrinliersbnd);
    extindyi = setdiff(extindyi,indyi);
    extindyi = extindyi(:);
    
    
    if ~isempty(extindyi),
        badidsi=[];
        ne = length(extindyi);
        extU = zeros(ne,rk);
        for iii=1:ne,
            
            
            idj = find(Wj(extindyi(iii),:)==1);
            ranok = 0;
            for kkk = 1:raniter,
                idjminset = idj(randperm(length(idj)));
                idjminset = idjminset(1:rk);
                extU(iii,:) = X(extindyi(iii),indyj(idjminset))/V(:,idjminset);
                res = [U; extU(iii,:)]*V-X([indyi extindyi(iii)],indyj);
                resi = res(end,:);
                resi = resi(W(extindyi(iii),indyj)==1);
                if  sum(abs(resi)<inlierbnd)>=nrinliersbnd
                    ranok = 1;
                    idok = abs(extU(iii,:)*V-X(extindyi(iii),indyj))<=inlierbnd & W(extindyi(iii),indyj);
                    extU(iii,:) = X(extindyi(iii),indyj(idok))/V(:,idok);
                    break;
                end
            end
            if ~ranok
                badidsi = [badidsi iii];
            else
                ok = 1;
            end
        end
        extU(badidsi,:)=[];
        extindyi(badidsi)=[];
        sol.indyi = [indyi extindyi'];
        sol.U = [U;extU];
        sol.Wloc = W(sol.indyi,indyj);
        sol.Wloc(abs(sol.U*V-X(sol.indyi,indyj))>inlierbnd & sol.Wloc==1)=-1;
    else
        ok = 0;
    end
else
    Wj = W(:,indyj);
    extindyi = find(sum(Wj,2)>=rk);
    extindyi = setdiff(extindyi,indyi);
    extindyi = extindyi(:);
    
    
    if ~isempty(extindyi),
        ok = 1;
        ne = length(extindyi);
        extU = zeros(ne,rk);
        for iii=1:ne,
            
            
            idj = find(Wj(extindyi(iii),:)==1);
            extU(iii,:) = X(extindyi(iii),indyj(idj))/V(:,idj);
        end
        sol.indyi = [indyi extindyi'];
        sol.U = [U;extU];
        sol.Wloc = W(sol.indyi,indyj);
        sol.Wloc(abs(sol.U*V-X(sol.indyi,indyj))>inlierbnd & sol.Wloc==1)=-1;
    else
        ok = 0;
    end
end

