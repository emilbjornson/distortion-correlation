%This Matlab script can be used to generate Figure 11 in the article:
%
%Emil Bjornson, Luca Sanguinetti, Jakob Hoydis, "Hardware Distortion
%Correlation Has Negligible Impact on UL Massive MIMO Spectral Efficiency,"
%IEEE Transactions on Communications, To appear.
%
%Download article: https://arxiv.org/abs/1811.02007
%
%This is version 1.0 (Last edited: 2018-10-18)
%
%License: This code is licensed under the GPLv2 license. If you in any way
%use this code for research that results in publications, please cite our
%original article listed above.


close all;
clear;

%Range of number of antennas
Mrange = round(logspace(1,3,10));

%Number of UEs
K = 1;

%Number of channel realizations
nbrOfRealizations = 1000;

%UE distortion parameter (when there is non-ideal hardware)
kappaOriginal = 0.99;

%BS distortion parameter with non-linearities (when there is non-ideal hardware)
alpha = 1/3;

%Back-off factor
b_off = db2pow(7);

%Signal-to-noise ratio
p = 1;


%Prepare to save simulation results
SE_MMSE = zeros(length(Mrange),nbrOfRealizations,4);
SE_MR = zeros(length(Mrange),nbrOfRealizations,4);


%% Go through all number of antennas
for m = 1:length(Mrange)
    
    %Write out the progress at every iterations
    disp(['Antenna iteration ' num2str(m) ' out of ' num2str(length(Mrange))]);
    
    %Extract number of antennas
    M = Mrange(m);
    
    %Generate channel realization
    H = (randn(M,K,nbrOfRealizations) + 1i*randn(M,K,nbrOfRealizations))/sqrt(2);
    
    %Create an identity matrix
    I_M = eye(M);
    
    %Go through all channel realizations
    for n = 1:nbrOfRealizations
        
        %Write out the progress at every 100 channel realizations
        if mod(n,100) == 0
            disp(['Channel ' num2str(n) ' realization of ' num2str(nbrOfRealizations)]);
        end
        
        %Consider four different cases of hardware impairments
        %l=1: Non-ideal BS and UE hardware
        %l=2: Non-ideal BS hardware, ideal UE hardware
        %l=3: Non-ideal UE hardware, ideal BS hardware
        %l=4: Ideal BS and UE hardware
        for l = 1:4
            
            if (l == 1) || (l == 2)
                a = alpha/(p*K*b_off);
            else
                a = 0;
            end
            
            if (l == 1) || (l == 3)
                kappa = kappaOriginal;
            else
                kappa = 1;
            end
            
            
            %Extract one channel realization
            Hreal = H(:,:,n);
            
            %Compute correlation matrix of received signal
            Cuu = p*(Hreal*Hreal');
            
            %Compute C_{eta eta} using (23)
            Cee = 2*a^2*Cuu.*Cuu.*conj(Cuu);
            
            %Compute effective channel using (21)
            D = eye(M)-2*a*diag(diag(Cuu));
            DH = D*Hreal;
            
            
            
            %% MMSE combining with correlated distortion
            V = p*(p*(DH*DH')+Cee+I_M)\DH;
            
            %Compute terms in numerator and denominator of the effective SINR
            channelproducts = abs(DH'*V).^2;
            numerators = kappa*p*diag(channelproducts)';
            denominators = p*sum(channelproducts,1) - numerators + sum(abs(V).^2,1) + real(sum((Cee*V).*conj(V),1));
            
            %Compute achievable SE
            SE_MMSE(m,n,l) = log2(1+numerators./denominators);
            
            
            %% MR combining with correlated distortion
            V = DH;
            
            %Compute terms in numerator and denominator of the effective SINR
            channelproducts = abs(DH'*V).^2;
            numerators = kappa*p*diag(channelproducts)';
            denominators = p*sum(channelproducts,1) - numerators + sum(abs(V).^2,1) + real(sum((Cee*V).*conj(V),1));
            
            %Compute achievable SE
            SE_MR(m,n,l) = log2(1+numerators./denominators);
            
            
        end
        
    end
    
end

%Compute SE per UE by averaging over channel realizations
avgSE_MMSE = mean(SE_MMSE,2);
avgSE_MR = mean(SE_MR,2);


%% Plot Figure 11
figure;
hold on; box on;
plot(Mrange,avgSE_MMSE(:,1,4),'ro-','LineWidth',1);
plot(Mrange,avgSE_MMSE(:,1,2),'ks--','LineWidth',1);
plot(Mrange,avgSE_MMSE(:,1,3),'bd-','LineWidth',1);
plot(Mrange,avgSE_MMSE(:,1,1),'r-','LineWidth',1);
xlabel('Number of antennas ($M$)','Interpreter','Latex');
ylabel('SE per UE [bit/s/Hz]','Interpreter','Latex');
legend({'Ideal','Non-ideal BS hardware','Non-ideal UE hardware','Non-ideal BS/UE hardware'},'Location','NorthWest','Interpreter','Latex');
set(gca,'XScale','log');


figure;
hold on; box on;
plot(Mrange,avgSE_MR(:,1,4),'ro-','LineWidth',1);
plot(Mrange,avgSE_MR(:,1,2),'ks--','LineWidth',1);
plot(Mrange,avgSE_MR(:,1,3),'bd-','LineWidth',1);
plot(Mrange,avgSE_MR(:,1,1),'r-','LineWidth',1);
xlabel('Number of antennas ($M$)','Interpreter','Latex');
ylabel('SE per UE [bit/s/Hz]','Interpreter','Latex');
legend({'Ideal','Non-ideal BS hardware','Non-ideal UE hardware','Non-ideal BS/UE hardware'},'Location','NorthWest','Interpreter','Latex');
set(gca,'XScale','log');
