function [Cap,Ptot]=ZFDPCapacity_MultipleAntenna(H,rxant,Power)
r=length(rxant);
rowidx=cumsum([0 rxant]);
d=[];
for n=1:r
    Hn=H(rowidx(n)+1:rowidx(n)+rxant(n),:);
    if n==1
        Hnbar=[];
        d=svd(Hn);
    else
        Hnbar=H(rowidx(1)+1:rowidx(n-1)+rxant(n-1),:);
        Vn=null(Hnbar);
        d=[d;svd(Hn*Vn)];
    end
end
d=abs(d).^2;
[wline,Ptot]=wfill(1./d,Power,1e-7);
Cap=sum(max(log2(wline*d),0));
