function [wline,ptot]=wfill(vec, pcon, tol)
%   Inc, 1991.

N=length(vec);

%first step of water filling
wline=min(vec)+pcon/N; %initial waterline
ptot=sum(max(wline-vec,0)); %total power for current waterline
 
%gradual water filling
while abs(pcon-ptot)>tol
    wline=wline+(pcon-ptot)/N;
    ptot=sum(max(wline-vec,0));
end

