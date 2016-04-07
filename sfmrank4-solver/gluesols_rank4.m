function [cursols,goodsols] = gluesols_rank4(X,W,cursols,goodsols,params)


%disp(['size before glue: ' num2str(length(cursols))])

N = params.N;
M = params.M;
raniter = params.glueraniter;
nn = length(cursols);
ok_indyi = zeros(1,N);
ok_indyj = zeros(1,M);
okind_indyi = zeros(nn,N);
okind_indyj = zeros(nn,M);


for iii = 1:nn,
    ok_indyi(cursols{iii}.indyi)=ok_indyi(cursols{iii}.indyi)+1;
    ok_indyj(cursols{iii}.indyj)=ok_indyj(cursols{iii}.indyj)+1;
    okind_indyi(iii,cursols{iii}.indyi)=1;
    okind_indyj(iii,cursols{iii}.indyj)=1;
end

gok_indyi = find(ok_indyi>1);
gok_indyj = find(ok_indyj>1);

badids = [];
gluesols = {};
ok = 0;
for kkk = 1:raniter,
    if ~isempty(gok_indyi),
        idis = find(okind_indyi(:,gok_indyi(randperm(length(gok_indyi),1)))==1);
        if length(idis)>=2,
            idis = idis(randperm(length(idis),2));
            [sol,ok] = gluei_rank4_ransac(X,W,cursols{idis(1)},cursols{idis(2)},params);
            if ok
                %disp('glue!');
                 for kkk = 1:params.bundleiter,
                    sol=bundlerank4_onestep(X,sol,params);
                    sol = updateW_rank4(X,sol,params);
                end
  
                maxN = length(sol.indyi);
                maxM = length(sol.indyj);
                if (maxN==params.finN) && (maxM==params.finM),
                    goodsols = [goodsols {sol}];
                else
                    gluesols = [gluesols {sol}];
                end
                badids = [badids idis];
            end
        end
    end
    if ok,
        break;
    end
    if ~isempty(gok_indyj),
        idis = find(okind_indyj(:,gok_indyj(randperm(length(gok_indyj),1)))==1);
        if length(idis)>=2,
            idis = idis(randperm(length(idis),2));
            [sol,ok] = gluei_rank4_ransac(X,W,cursols{idis(1)},cursols{idis(2)},params);
            if ok
                %disp('glue!');
                 for kkk = 1:params.bundleiter,
                    sol=bundlerank4_onestep(X,sol,params);
                    sol = updateW_rank4(X,sol,params);
                end

                maxN = length(sol.indyi);
                maxM = length(sol.indyj);
                if (maxN==params.finN) && (maxM==params.finM),
                    goodsols = [goodsols {sol}];
                else
                    gluesols = [gluesols {sol}];
                end
                badids = [badids idis];
            end
        end
    end
    if ok,
        break;
    end
end
cursols = [cursols gluesols];
cursols(badids)=[];
%disp(['size after glue: ' num2str(length(cursols))])


