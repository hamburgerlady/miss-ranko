R = 4;
I0 = ones(R);
nn = 8;
rh1sets = cell(nn-R+1);
rh1sets{1,1}={I0};
for jjj = (R+1):nn,
    rh1sets{1,jjj-R+1}={[rh1sets{1,jjj-R}{1} ones(R,1)]};
end

for iii = (R+1):nn,
    rh1sets{iii-R+1,1}={[rh1sets{iii-R,1}{1}; ones(1,R)]};
end





for iii=(R+1):nn,
    for jjj=(R+1):nn,
        disp([iii jjj])
    newsets = {};

    
    oldsets = rh1sets{iii-R,jjj-R+1};
    indy = nchoosek(1:jjj,R);
    for kkk = 1:length(oldsets),
        for lll = 1:size(indy,1),
            I = [oldsets{kkk}; zeros(1,jjj)];
            I(end,indy(lll,:))=1;
            allnok = 1;
            for mmm = 1:length(newsets),
                if isSameIndexSet(I,newsets{mmm}),
                    allnok = 0;
                    break;
                end
            end
            if allnok 
                newsets=[newsets {I}];
            end
        end
    end
    
    oldsets = rh1sets{iii-R+1,jjj-R};
    indy = nchoosek(1:iii,R);
    for kkk = 1:length(oldsets),
        for lll = 1:size(indy,1),
            I = [oldsets{kkk} zeros(iii,1)];
            I(indy(lll,:),end)=1;
            allnok = 1;
            for mmm = 1:length(newsets),
                if isSameIndexSet(I,newsets{mmm}),
                    allnok = 0;
                    break;
                end
            end
            if allnok 
                newsets=[newsets {I}];
            end
        end
    end
    
    
    
    
    rh1sets{iii-R+1,jjj-R+1}=newsets;
    end
end

