clear
clc
M = 2; % number of antennas per user
K = 4:4:40; % number of users
N = 8; % number of Tx antennas

PowerdB=20; % Tx power  
Power=10.^(PowerdB./10);

nChannels=10000;


%% SUCCESSIVE ZERO-FORCING DIRTY PAPER CODING WITHOUT ANTENNA SELECTION
SumRate_C_Algorithm=zeros(length(K),length(Power));
SumRate_E_Algorithm=zeros(length(K),length(Power));
SumRate_D_Algorithm=zeros(length(K),length(Power));

%% COOPERATIVE ASSUMPTION
%CapacityCoop=zeros(length(K),length(Power));

%% ROUND-ROBIN APPROACH
%CapacityRR=zeros(length(K),length(Power));
%% MAIN
for iUser=1:length(K)
    for iChannel=1:nChannels
        H=(randn(M*K(iUser),N)+1i*randn(M*K(iUser),N))/sqrt(2);
        RxArray=M*ones(1,K(iUser));

        [SelectedUsers_C_Algorithm,temp]=C_algorithm(H,RxArray,Power);
        SumRate_C_Algorithm(iUser)=SumRate_C_Algorithm(iUser)+temp;

        [SelectedUsers_E_Algorithm,temp]=E_algorithm(H,RxArray,Power);
        SumRate_E_Algorithm(iUser)=SumRate_E_Algorithm(iUser)+temp;

        [SelectedUsers_D_Algorithm,temp]=D_algorithm(H,RxArray,Power);
        SumRate_D_Algorithm(iUser)=SumRate_D_Algorithm(iUser)+temp;

    end
end

SumRate_C_Algorithm=SumRate_C_Algorithm/nChannels;
SumRate_E_Algorithm=SumRate_E_Algorithm/nChannels;
SumRate_D_Algorithm=SumRate_D_Algorithm/nChannels;


plot(K,SumRate_C_Algorithm)
hold on 
plot(K,SumRate_E_Algorithm)
plot(K,SumRate_D_Algorithm)
xlabel('Number of users')
ylabel('Avarage Sum Rate (b/s/Hz)')
legend('C-algorith','E-algorithm', 'D-algorithm','Location','southeast')
saveas(gcf,'../results/Ergodic sum rate.png')