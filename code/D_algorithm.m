function [users,Cap]=D_algorithm(H,rxant,Power)
t=size(H,2);
r=length(rxant);
rowidx=cumsum([0 rxant]);

S=[];
userselected=[];
usersremain=(1:length(rxant));
%% Initialization, choose the user that has miximum capacity
H_diag_max=0;
for k=1:r
    Hk=H(rowidx(k)+1:rowidx(k)+rxant(k),:);
    temp=diag(Hk*Hk');
    temp=sort(temp,'descend');
    if (compvectors(temp,H_diag_max))
        H_diag_max=temp;
        selecteduser=k;
    end
end

userselected=[userselected selecteduser];   
usersremain(selecteduser)=[];    
Hselected=H(rowidx(selecteduser)+1:rowidx(selecteduser)+rxant(k),:);
[G,V]=grams(Hselected);
SelectedDiagonalElements=H_diag_max;


n=2; 
%% MAIN LOOP
while (n<=ceil(t/mean(rxant)));
    selecteduser=0;
    H_diag_max_temp=0;
    for k=1:length(usersremain)
        Hk=H(rowidx(usersremain(k))+1:rowidx(usersremain(k))+rxant(usersremain(k)),:);
        Hktilde=Hk-Hk*V'*V; % project Hk to the null space of V
        temp=[H_diag_max;diag(Hktilde*Hktilde')];
        temp=sort(temp,'descend');

        if (compvectors(temp,H_diag_max_temp))
            H_diag_max_temp=temp;
            selecteduser=k;
            temp1=Hk;
            temp2=Hktilde;
        end
    end
    if (selecteduser>0)
        [G,W]=grams(temp2);
        V=[V;W];
        userselected=[userselected usersremain(selecteduser)];
        usersremain(selecteduser)=[];
        Hselected=[Hselected ;temp1];
        H_diag_max=H_diag_max_temp;
    else
        break;
    end
    n=n+1;
end
rxantnew=rxant(userselected);
% Hnew=[];
% for m=1:length(userselected)
%     Hnew=[Hnew;Hselected(2*(m-1)+1:2*m,:)];
% end

[Cap,Total_Power]=ZFDPCapacity_MultipleAntenna(Hselected,rxantnew,Power);
users=userselected;
