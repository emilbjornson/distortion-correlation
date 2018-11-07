%This Matlab script can be used to generate Figure 5 in the article:
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

%UE distortion parameter
kappa = 0.99;

%BS distortion parameter with non-linearities
alpha = 1/3;

%Back-off factor
b_off = db2pow(7);

%Signal-to-noise ratio
SNR = 1;

%Number of antennas
M = 200;

%Range of number of UEs
Krange = [5 10 15];

%Number of channel realizations
nbrOfRealizations = 1000;


%Prepare to save simulation results
eigenvalues_corr = zeros(M,length(Krange),nbrOfRealizations);



%% Go through all number of UEs
for k = 1:length(Krange)
    
    %Output simulation progress
    disp([num2str(k) ' user set out of ' num2str(length(Krange))]);
    
    %Extract number of UEs
    K = Krange(k);
    
    %Generate channel realizations
    H = (randn(M,K,nbrOfRealizations) + 1i*randn(M,K,nbrOfRealizations))/sqrt(2);
    
    
    for n = 1:nbrOfRealizations
        
        %Extract channel realization
        Hreal = H(:,:,n);
        
        %Compute correlation matrix of signal and distortion
        a = alpha/(SNR*K*b_off);
        Cuu = SNR*(Hreal*Hreal');
        Cee = 2*a^2*Cuu.*Cuu.*conj(Cuu);
        
        %Compute eigenvalues
        eigenvalues_corr(:,k,n) = sort(eig(Cee),'descend')/trace(Cee);
        
    end
    
end



%% Plot Figure 5
figure;
hold on; box on;

plot(1:M,mean(eigenvalues_corr(:,1,:),3),'k--','LineWidth',1);
plot(1:M,mean(eigenvalues_corr(:,2,:),3),'r-','LineWidth',1);
plot(1:M,mean(eigenvalues_corr(:,3,:),3),'b-.','LineWidth',1);

set(gca,'YScale','log');
ylim([1e-4 1e-1]);
xlabel('Eigenvalue index (decaying order)','Interpreter','Latex');
ylabel('Eigenvalue','Interpreter','Latex');
legend({'$K=5$','$K=10$','$K=15$'},'Interpreter','Latex');
