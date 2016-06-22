function resgrad = calcresgrad_rankn(sol,robust)

if nargin<2,
    robust = 0;
end

U = sol.U;
V = sol.V;
N = size(U,1);
M = size(V,2);
rk = size(U,2);

if robust,
    [idsi,idsj] = find(sol.Wloc==1);
else
    [idsi,idsj] = find(sol.Wloc~=0);
end
nn = length(idsi);


spdata = zeros(2*nn*rk,3);
count = 0;

for id=1:nn,
    ii = idsi(id);
    jj = idsj(id);
    for kk = 1:rk,
        if ii>rk
            count = count + 1;
            spdata(count,:)=[id rk*(ii-rk-1)+kk V(kk,jj)];
            %resgrad(id,rk*(ii-rk-1)+kk)=V(kk,jj);
        end
        count = count + 1;
        spdata(count,:)=[id rk*(N-rk)+rk*(jj-1)+kk U(ii,kk)];
            %resgrad(id,rk*(N-rk)+rk*(jj-1)+kk)=U(ii,kk);
       
    end
end

resgrad = sparse(spdata(1:count,1),spdata(1:count,2),spdata(1:count,3),nn,(N-rk)*rk+M*rk);


%res = sol.U*sol.V-X(sol.indyi,sol.indyj);
%res = res(sol.Wloc==1);

