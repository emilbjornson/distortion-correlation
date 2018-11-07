%This Matlab script can be used to generate Figure 4 in the article:
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


%Maximum number of antennas
Mmax = 200;

%Number of channel realizations
nbrOfRealizations = 1000;


%UE distortion parameter
kappa = 0.99;

%BS distortion parameter with non-linearities
alpha = 1/3;

%Back-off factor
b_off = db2pow(7);

%Signal-to-noise ratio
SNR = 1;


%Generate channel realizations
h = (randn(Mmax,nbrOfRealizations)+1i*randn(Mmax,nbrOfRealizations))/sqrt(2);


%Prepare to save simulation results
signalpower_DAMR = zeros(nbrOfRealizations,Mmax);
distortionpower_DAMR = zeros(nbrOfRealizations,Mmax);
signalpower_cancelation = zeros(nbrOfRealizations,Mmax);


%% Go through all channel realizations
for n = 1:nbrOfRealizations
    
    %Display simulation progress
    disp(['Iteration ' num2str(n) ' out of ' num2str(nbrOfRealizations)]);
    
    %Compute C_{uu} and C_{eta eta}
    Cuu = SNR*h(:,n)*h(:,n)';
    Cee = 2*(alpha/(SNR*b_off))^2*Cuu.*Cuu.*conj(Cuu);
    
    
    %Go through all number of natennas
    for m = 1:Mmax
        
        %Compute DA-MR with normalization
        D = eye(m)-2*alpha/(SNR*b_off)*diag(diag(Cuu(1:m,1:m)));
        DAMR = D*h(1:m,n)/norm(D*h(1:m,n));
        
        %Desired signal power with DA-MR
        signalpower_DAMR(n,m) = kappa*real(DAMR'*D*Cuu(1:m,1:m)*D'*DAMR);
        
        %BS distortion power with DA-MR
        distortionpower_DAMR(n,m) = real(DAMR'*Cee(1:m,1:m)*DAMR);
        
        %Compute DA-ZF
        basis = orth(Cee(1:m,1:m));
        DAZF = (eye(m)-basis/(basis'*basis)*basis')*DAMR;
        
        %Set the normalization for DA-ZF
        if norm(DAZF)==0 %If M=1
            DAZF = zeros(size(DAZF));
        else %If M>1
            DAZF = DAZF/norm(DAZF);
        end
        
        %Desired signal power with DA-ZF
        signalpower_cancelation(n,m) = kappa*real(DAZF'*D*Cuu(1:m,1:m)*D'*DAZF);
        
    end
    
end


%% Plot Figure 4
figure;
hold on; box on;
plot(1:Mmax,mean(signalpower_DAMR,1),'r','LineWidth',1);
plot(1:Mmax,mean(signalpower_cancelation,1),'b--','LineWidth',1);
plot(1:Mmax,mean(distortionpower_DAMR,1),'k-.','LineWidth',1);
xlabel('Number of antennas ($M$)','Interpreter','Latex');
ylabel('Signal/distortion power over noise','Interpreter','Latex');
legend({'DA-MR: Desired signal','DA-ZF: Desired signal','DA-MR: BS distortion'},'Interpreter','Latex','Location','NorthWest');
ylim([0 120]);
