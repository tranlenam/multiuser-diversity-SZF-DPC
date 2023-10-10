function [users,Cap]=C_algorithm(H,rxant,Power)
N=size(H,2);
r=length(rxant);
rowidx=cumsum([0 rxant]);

S=[];
userselected=[];
usersremain=(1:length(rxant));
% Initialization, choose the user that has miximum capacity
Cmax=0;
for k=1:r
    C=CoopCapacity(H(rowidx(k)+1:rowidx(k)+rxant(k),:),Power);
    if (C>Cmax)
        Cmax=C;
        selecteduser=k;
    end
end

userselected=[userselected selecteduser];   
usersremain(selecteduser)=[];    
Hselected=H(rowidx(selecteduser)+1:rowidx(selecteduser)+rxant(selecteduser),:);
[G,V]=grams(Hselected);
eigens=svd(Hselected);



n=2; 
% LOOP
while (n<=ceil(N/mean(rxant)));
    selecteduser=0;
    for k=1:length(usersremain)
        Hk=H(rowidx(usersremain(k))+1:rowidx(usersremain(k))+rxant(usersremain(k)),:);
        Hktilde=Hk-Hk*V'*V; % project Hk to the null space of V
        temp=[eigens;svd(Hktilde)];
        u=abs(temp).^2;
        wline=wfill(1./u,Power,1e-7);
        C=sum(max(log2(wline*u),0));
        
        if (C>Cmax)
            Cmax=C;
            selecteduser=k;
            temp1=Hk;
            temp2=Hktilde;
            eigen_temp=temp;
        end
    end
    if (selecteduser>0)
        [G,W]=grams(temp2);
        V=[V;W];
        userselected=[userselected usersremain(selecteduser)];
        usersremain(selecteduser)=[];
        Hselected=[Hselected temp1];
        eigens=eigen_temp;
    end
    n=n+1;
end
users=userselected;
Cap=Cmax;