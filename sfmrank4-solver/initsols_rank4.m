function cursols = initsols_rank4(X,W,params)


cursols = initnnsols_rank4(X,W,params);
for iii = 1:length(cursols),
    for kkk = 1:2*params.bundleiter,
        cursols{iii}=bundlerank4_onestep(X,cursols{iii},params);
        cursols{iii} = updateW_rank4(X,cursols{iii},params);
    end
end
