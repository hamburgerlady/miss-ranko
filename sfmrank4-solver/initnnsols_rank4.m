function cursols = initnnsols_rank4(X,W,params)

gksize = params.gksize;
gksigge = params.gksigge;
minhalfn = params.minhalfn;
cutty = params.cutty;
inlierbnd = params.inlierbnd;
nrinliersbnd = params.nrinliersbnd;
rk = params.rk;
N = params.N;
M = params.M;
initnn = params.initnn;

[xx,yy]=meshgrid(-gksize:gksize);
gk = exp(-(xx.^2+yy.^2)/2/gksigge^2);
W2 = conv2(W,gk,'same');
W2(1:minhalfn,:)=0;
W2(:,1:minhalfn)=0;
W2(end-minhalfn:end,:)=0;
W2(:,end-minhalfn:end)=0;
W2 = W2/max(W2(:));


W2c_ind = find(W2>cutty);
nn = length(W2c_ind);
oksols = 0;
cursols = cell(1,initnn);

while oksols<initnn
    [idi,idj]=ind2sub([N,M],W2c_ind(ceil(rand*nn)));
    indyi = (idi-minhalfn):(idi+minhalfn-1);
    indyj = (idj-minhalfn):(idj+minhalfn);
    I = [W(indyi,indyj);ones(1,2*minhalfn+1)];
    XI = [X(indyi,indyj);ones(1,2*minhalfn+1)];
    indyi = [indyi N];
    
    Iout = getrandomminimalh1config(I,rk);
    [U,V,ok] = solve_r4h1(XI,Iout);
    if ok,
        
        inliers = abs(U*V-XI)<inlierbnd;
        if all(sum(inliers,1)>=nrinliersbnd) && all(sum(inliers,2)>=nrinliersbnd)
            [ii,jj]=find(~inliers);
            Wtmp = ones(N,M);
            Wtmp(sub2ind([N,M],indyi(ii),indyj(jj)))=-1 ;
            Wtmp(W==0)=0;
            
            lr = U(end,:);
            S = [null(lr) lr'/(lr*lr')];
            U = U*S;
            V = S\V;
            
            
            
            oksols = oksols+1;
            sol.U = U;
            sol.V = V;
            sol.indyi = indyi;
            sol.indyj = indyj;
            sol.ground = [mean(indyi(1:end-1)) mean(indyj)];
            sol.Wloc = Wtmp(indyi,indyj);
            sol.static = 0;
            cursols{oksols}=sol;
            
            
            
            
        end
        
        
    end
    
    
    
end



