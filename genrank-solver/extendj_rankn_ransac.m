function [sol,ok] = extendj_rankn_ransac(X,W,sol,params)



inlierbnd = params.inlierbnd;
robust = params.robust;

ok = 0;
rk = params.rk;
U = sol.U;
V = sol.V;

%keyboard

indyi = sol.indyi;
indyj = sol.indyj;

if robust,
    nrinliersbnd = params.nrinliersbnd;
    raniter = params.extendraniter;

    Wi = W(indyi,:);
    extindyj = find(sum(Wi)>=nrinliersbnd);
    extindyj = setdiff(extindyj,indyj);
    
    if ~isempty(extindyj),
        badidsj=[];
        ne = length(extindyj);
        extV = zeros(rk,ne);
        for jjj=1:ne,
            idi = find(Wi(:,extindyj(jjj))==1);
            ranok = 0;
            for kkk = 1:raniter,
                idiminset = idi(randperm(length(idi)));
                idiminset = idiminset(1:rk);
                
                
                extV(:,jjj) = U(idiminset,:)\X(indyi(idiminset),extindyj(jjj));
                res = U*[V extV(:,jjj)]-X(indyi,[indyj  extindyj(jjj)]);
                resj = res(:,end);
                resj = resj(W(indyi,extindyj(jjj))==1);
                if  sum(abs(resj)<inlierbnd)>=nrinliersbnd
                    ranok = 1;
                     idok = abs(U*extV(:,jjj)-X(indyi,extindyj(jjj)))<=inlierbnd & W(indyi,extindyj(jjj));
                     extV(:,jjj) = U(idok,:)\X(indyi(idok),extindyj(jjj));
                    break;
                end
            end
            if ~ranok
                badidsj = [badidsj jjj];
            else
                ok = 1;
            end
        end
        extV(:,badidsj)=[];
        extindyj(badidsj)=[];
        
        sol.indyj = [indyj extindyj];
        sol.V = [V extV];
        sol.Wloc = W(indyi,sol.indyj);
        sol.Wloc(abs(U*sol.V-X(indyi,sol.indyj))>inlierbnd & sol.Wloc==1)=-1;
    else
        ok = 0;
    end
    
    
else
    Wi = W(indyi,:);
    extindyj = find(sum(Wi)>=rk);
    extindyj = setdiff(extindyj,indyj);
    
    if ~isempty(extindyj),
        ok=1;
        ne = length(extindyj);
        extV = zeros(rk,ne);
        for jjj=1:ne,
            idi = find(Wi(:,extindyj(jjj))==1);
            extV(:,jjj) = U(idi,:)\X(indyi(idi),extindyj(jjj));
        end
        
        sol.indyj = [indyj extindyj];
        sol.V = [V extV];
        sol.Wloc = W(indyi,sol.indyj);
        sol.Wloc(abs(U*sol.V-X(indyi,sol.indyj))>inlierbnd & sol.Wloc==1)=-1;
    else
        ok = 0;
    end
end