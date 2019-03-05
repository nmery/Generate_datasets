function [ddh_comp]=compos_rosebel(ddh,lc)
%lc: length of composites

n=size(ddh,1);
cl=0;cd=0;
ifin=2;
icomp=1;
comp=0;
lcomp=0;
xdeb=ddh(1,1:3);
ddh_comp=ones(n,4)*-1;

while ifin<n
    cl=sqrt(sum((ddh(ifin,1:3)-ddh(ifin-1,1:3)).^2));
    lcomp=lcomp+cl;
    if lcomp<lc
        comp=comp+cl*ddh(ifin-1,end);
    elseif lcomp>5*lc % on vient sûrement de changer de forage;
        lcomp=0;
        comp=0;
        xdeb=ddh(ifin,1:3);
    else   % on a franchi lc
        comp=comp+(lc-(lcomp-cl))*ddh(ifin-1,end);
        xfin=xdeb+(lc/lcomp)*(ddh(ifin,1:3)-xdeb);
        xcomp=(xdeb+xfin)/2;
        ddh_comp(icomp,:)=[xcomp comp/lc];
        icomp=icomp+1;
        xdeb=xfin;
        lcomp=lcomp-lc;
        comp=lcomp*ddh(ifin,end);
    end
    ifin=ifin+1;
end
id=ddh_comp(:,end)>=0;
ddh_comp=ddh_comp(id,:);

    
    