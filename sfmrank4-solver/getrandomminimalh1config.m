function Iout = getrandomminimalh1config(I,r)

[N,M]=size(I);
Iout = 0*I;
nrmeas = r*N+r*M-r*r;
if sum(I(:))>=nrmeas,
    idde=find(I);
    rp = randperm(length(idde));
    Iout(idde(rp(1:nrmeas)))=1;    
end
