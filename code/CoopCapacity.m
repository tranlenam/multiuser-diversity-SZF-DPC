function [Cap,egmode]=CoopCapacity(H,Power)
egmode=svd(H);
u=abs(egmode).^2;
wline=wfill(1./u,Power,1e-7);
Cap=sum(max(log2(wline*u),0));
