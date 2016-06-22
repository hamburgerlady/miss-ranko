function cursols = initsols_rankn(X,W,params)


cursols = initnnsols_rankn(X,W,params);
for iii = 1:length(cursols),
    for kkk = 1:2*params.bundleiter,
        cursols{iii}=bundlerankn_onestep(X,cursols{iii},params);
        cursols{iii} = updateW_rankn(X,cursols{iii},params);
    end
end
