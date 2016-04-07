function [U,V,ok] = solve_r4h1(X,I)

mini = 1e-11;
r = 4;
I0=I;
ids=[];
[N,M]=size(I);
U = zeros(N,r);
V = zeros(r,M);

idi = 1:N;
idj = 1:M;
ok = 1;
while length(I)>r,
    rs = sum(I);
    tmpids=find(rs==r);
    if ~isempty(tmpids),
        ids=[ids [2;idj(tmpids(1))]];
        I(:,tmpids(1))=[];
        idj(tmpids(1))=[];
    else
        cs = sum(I,2);
        tmpids=find(cs==r);
        if ~isempty(tmpids),
            ids=[ids [1;idi(tmpids(1))]];
            I(tmpids(1),:)=[];
            idi(tmpids(1))=[];
        else
            ok = length(I)==r & sum(I(:))==r^2;
            I = [];
        end
    end   
end

if ok,
    ok = length(I)==r & sum(I(:))==r^2;
end

if ok,
    idi0 = idi;
    idj0 = idj;
    [u,s,v]=svd(X(idi,idj));
    U(idi,:)=u(:,1:r)*s(1:r,1:r);
    V(:,idj)=v(:,1:r)';
    for iii = length(ids):-1:1,
        if ids(1,iii)==1,
           tmpids = find(I0(ids(2,iii),idj));
           if min(svd(V(:,idj(tmpids))))>mini,
                U(ids(2,iii),:)=X(ids(2,iii),idj(tmpids))/(V(:,idj(tmpids)));
                idi = sort([idi ids(2,iii)]);
           else
               ok = 0;
               break;
           end
               
        else
           tmpids = find(I0(idi,ids(2,iii)));
           if min(svd(U(idi(tmpids),:)))>mini,
                V(:,ids(2,iii))=(U(idi(tmpids),:))\X(idi(tmpids),ids(2,iii));
                idj = sort([idj ids(2,iii)]);
           else
               ok = 0;
               break;
           end
        end
        
    end
    %keyboard
    
end    
    
