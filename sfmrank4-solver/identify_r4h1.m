function [ok,ids] = identify_r4h1(I)

ids=[];
r = 4;
ok = 1;
while length(I)>r,
    rs = sum(I);
    tmpids=find(rs==r);
    if ~isempty(tmpids),
        ids=[ids [tmpids(1);0]];
        I(:,tmpids(1))=[];
    else
        cs = sum(I');
        tmpids=find(cs==r);
        if ~isempty(tmpids),
            ids=[ids [0;tmpids(1)]];
            I(tmpids(1),:)=[];
        else
            ok = length(I)==r & sum(I(:))==r^2;
            I = [];
        end
    end   
end

if ok,
    ok = length(I)==r & sum(I(:))==r^2;
end
