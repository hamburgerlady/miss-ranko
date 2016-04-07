function plotsolsboxes(sols,W,cc)
if nargin<3,
    cc = 'b';
end

N=size(W,1);
nn = length(sols);
[ii,jj]=find(W==1);
plot(jj,ii,'.b');
for iii=1:nn,
    x1 = min(sols{iii}.indyj);x2 = max(sols{iii}.indyj);
    y1 = min(sols{iii}.indyi);y2 = max(setdiff(sols{iii}.indyi,N));
       [idi,idj]=find(W(sols{iii}.indyi,sols{iii}.indyj)==1);
    plot(sols{iii}.indyj(idj),sols{iii}.indyi(idi),['.g']);
    
    [idi,idj]=find(sols{iii}.Wloc==-1);
    plot(sols{iii}.indyj(idj),sols{iii}.indyi(idi),['.r']);
    
    ll=rectangle('position',[x1 y1 x2-x1 y2-y1]);
    set(ll,'EdgeColor',cc);
    set(ll,'LineWidth',2);
    

    
end

