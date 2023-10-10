function [users,Cap]=E_algorithm(H,rxant,Power)
t=size(H,2);
r=length(rxant);
rowidx=cumsum([0 rxant]);

S=[];
userselected=[];
usersremain=(1:length(rxant));
% Initialization, choose the user that has miximum capacity
%% The first step, choose the user whose set of single values majorizes the remaining ones 
Vmax=0;
for k=1:r
    temp=svd(H(rowidx(k)+1:rowidx(k)+rxant(k),:));
    temp=abs(temp).^2;
    if (compvectors(temp,Vmax))
        Vmax=temp;
        selecteduser=k;
    end
end
%%
userselected=[userselected selecteduser];   
usersremain(selecteduser)=[];    
Hselected=H(rowidx(selecteduser)+1:rowidx(selecteduser)+rxant(selecteduser),:);
[G,V]=grams(Hselected);
eigens=svd(Hselected);



n=2; 
%% LOOP
while (n<=ceil(t/mean(rxant)));
    selecteduser=0;
    for k=1:length(usersremain)
        Hk=H(rowidx(usersremain(k))+1:rowidx(usersremain(k))+rxant(usersremain(k)),:);
        Hktilde=Hk-Hk*V'*V; % project Hk to the null space of V
        temp=[eigens;svd(Hktilde)];
        temp=sort(temp,'descend');
        u=abs(temp).^2;
        
        if (compvectors(u,Vmax))
            Vmax=u;
            selecteduser=k;
            temp1=Hk;
            temp2=Hktilde;
            eigens_temp=temp;
        end
    end
    if (selecteduser>0)
        [G,W]=grams(temp2);
        V=[V;W];
        userselected=[userselected usersremain(selecteduser)];
        usersremain(selecteduser)=[];
        Hselected=[Hselected;temp1];
        eigens=eigens_temp;
    end
    n=n+1;
end
%% Compute the capacity
users=userselected;
wline=wfill(1./Vmax,Power,1e-7);
Cap=sum(max(log2(wline*Vmax),0));
