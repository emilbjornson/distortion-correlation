%This Matlab script can be used to generate Figure 8 in the article:
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

%Set maximum ADC resolution
ADCresolution = 6;

%Load the quantization levels for non-uniform quantization of a standard
%Gaussian random variable, obtained using the Lloyd algorithm
load quantizationLevels;

%UE distortion parameter
kappa = 0.99;

%BS distortion parameter with non-linearities
alpha = 1/3;

%Back-off factor
b_off = db2pow(7);

%Number of antennas
M = 100;

%Range of number of UEs
Krange = 1:50;

%Number of channel realizations
nbrOfRealizations = 400;

%Number of signal transmissions, used to compute Cee by Monte-Carlo methods
signalTransmissions = 10000;

%Signal-to-noise ratio
SNR = 1;

%Create an identity matrix
I_M = eye(M);


%Prepare to save simulation results
sumSE_corr = zeros(length(Krange),2);
sumSE_uncorr = zeros(length(Krange),2);


%% Go through all channel realizations
for n = 1:nbrOfRealizations
    
    %Write out the progress at every iterations
    disp(['Progress: ' num2str(n) ' out of ' num2str(nbrOfRealizations) ' realizations']);
    
    for k = 1:length(Krange)
        
        %Extract number of UEs
        K = Krange(k);
        
        %Scaling factor that is used to make the real and imaginary part of
        %the received signal standard Gaussian random variables
        scaling = (SNR*K)/2*eye(M);
        
        %Generate channel realization
        H = (randn(M,K)+1i*randn(M,K))/sqrt(2);
        
        %Generate the signals to be transmitted in the Monte-Carlo transmission
        S = (randn(K,signalTransmissions)+1i*randn(K,signalTransmissions))/sqrt(2);
        
        %Compute %Compute C_{uu} for the given channel realization
        Cuu = SNR*(H*H');
        
        %Compute the noise-free received signal
        U = sqrt(SNR)*H*S;
        
        %Normalize the received signal so that the real and imaginary parts
        %have unit variance, not for the given channel realization but on the
        %average
        U_normalized = sqrtm(scaling)\U;
        
        %Apply non-linear distortion
        U_nonlinear = U_normalized - alpha/(2*b_off)*abs(U_normalized).^2.*U_normalized;
        
        
        
        %Compute the vector with diagonal elements of D using expression
        %given between (7) and (8)
        dvec = zeros(M,1);
        
        
        U_quant = zeros(M,signalTransmissions);
        
        for l = 1:signalTransmissions
            
            %Quantization of real and imaginary parts of received signal
            real_index = quantiz(real(U_nonlinear(:,l)),partitions{ADCresolution});
            imag_index = quantiz(imag(U_nonlinear(:,l)),partitions{ADCresolution});
            
            %Reconstruct the quantized signal using the quantization levels
            U_quant(:,l) = sqrtm(scaling)*(codebooks{ADCresolution}(real_index+1)+1i*codebooks{ADCresolution}(imag_index+1));
            
            %Approximate the elements of D using a sample correlation matrix
            dvec = dvec + (U_quant(:,l).*conj(U(:,l)))/signalTransmissions;
            
        end
        
        %Compute estimate of D matrix, including normalization
        dvec = dvec./diag(Cuu);
        D = diag(real(dvec));
        
        %Compute C_{eta eta} for the given channel realization using
        %Monte-Carlo methods
        Cee = zeros(M,M);
        
        for l = 1:signalTransmissions
            
            %Approximate C_{eta eta} by a sample correlation matrix
            Cee = Cee + ((U_quant(:,l)-D*U(:,l))*(U_quant(:,l)-D*U(:,l))')/signalTransmissions;
            
        end
        
        %Create the approximate diagonal version of C_{eta eta}
        Ceediag = diag(diag(Cee));
        
        %Product of D and H
        DH = D*H;
        
        
        %% DA-MMSE combining with correlated distortion
        V = SNR*(SNR*(DH*DH')+Cee+I_M)\DH;
        
        %Compute terms in numerator and denominator of the effective SINR
        channelproducts = abs(DH'*V).^2;
        numerators = kappa*SNR*diag(channelproducts)';
        denominators = SNR*sum(channelproducts,1) - numerators + sum(abs(V).^2,1) + real(sum((Cee*V).*conj(V),1));
        
        %Compute achievable SE
        sumSE_corr(k,1) = sumSE_corr(k,1) + sum(log2(1+numerators./denominators))/nbrOfRealizations;
        
        
        %% DA-MR combining with correlated distortion
        V = DH;
        
        %Compute terms in numerator and denominator of the effective SINR
        channelproducts = abs(DH'*V).^2;
        numerators = kappa*SNR*diag(channelproducts)';
        denominators = SNR*sum(channelproducts,1) - numerators + sum(abs(V).^2,1) + real(sum((Cee*V).*conj(V),1));
        
        %Compute achievable SE
        sumSE_corr(k,2) = sumSE_corr(k,2) + sum(log2(1+numerators./denominators))/nbrOfRealizations;
        
        
        
        %% DA-MMSE combining with uncorrelated distortion
        V = SNR*(SNR*(DH*DH')+Ceediag+I_M)\DH;
        
        %Compute terms in numerator and denominator of the effective SINR
        channelproducts = abs(DH'*V).^2;
        numerators = kappa*SNR*diag(channelproducts)';
        denominators = SNR*sum(channelproducts,1) - numerators + sum(abs(V).^2,1) + real(sum((Ceediag*V).*conj(V),1));
        
        %Compute achievable SE
        sumSE_uncorr(k,1) = sumSE_uncorr(k,1) + sum(log2(1+numerators./denominators))/nbrOfRealizations;
        
        
        
        %% DA-MR combining with uncorrelated distortion
        V = DH;
        
        %Compute terms in numerator and denominator of the effective SINR
        channelproducts = abs(DH'*V).^2;
        numerators = kappa*SNR*diag(channelproducts)';
        denominators = SNR*sum(channelproducts,1) - numerators + sum(abs(V).^2,1) + real(sum((Ceediag*V).*conj(V),1));
        
        %Compute achievable SE
        sumSE_uncorr(k,2) = sumSE_uncorr(k,2) + sum(log2(1+numerators./denominators))/nbrOfRealizations;
        
        
    end
    
end


%Compute average SE per UE
avgSE_MMSE_uncorr = sumSE_uncorr(:,1)./Krange';
avgSE_MMSE_corr = sumSE_corr(:,1)./Krange';
avgSE_MR_uncorr = sumSE_uncorr(:,2)./Krange';
avgSE_MR_corr = sumSE_corr(:,2)./Krange';


%% Plot Figure 8
figure;
hold on; box on;
plot(Krange(1),avgSE_MMSE_uncorr(1),'bo--','LineWidth',1);
plot(Krange(1),avgSE_MMSE_corr(1),'ro-','LineWidth',1);
plot(Krange(1),avgSE_MR_uncorr(1),'b--','LineWidth',1);
plot(Krange(1),avgSE_MR_corr(1),'k-','LineWidth',1);
xlabel('Number of UEs ($K$)','Interpreter','Latex');
ylabel('SE per UE [bit/s/Hz]','Interpreter','Latex');
legend({'DA-MMSE, uncorr','DA-MMSE, corr','DA-MR, uncorr','DA-MR, corr'},'Location','SouthWest','Interpreter','Latex');


y = [avgSE_MMSE_uncorr(1:5)'; avgSE_MMSE_corr(1:5)'-avgSE_MMSE_uncorr(1:5)']';
h = area(Krange(1:5), y, 0);
set(h(1), 'FaceColor', 'none');
set(h, 'LineStyle', 'none');
h(2).FaceColor = [0.85 0.85 0.85];


y = [avgSE_MR_uncorr(1:5)'; avgSE_MR_corr(1:5)'-avgSE_MR_uncorr(1:5)']';
h = area(Krange(1:5), y, 0);
set(h(1), 'FaceColor', 'none');
set(h, 'LineStyle', 'none');
h(2).FaceColor = [0.85 0.85 0.85];

plot(Krange,avgSE_MMSE_uncorr,'b--','LineWidth',1);
plot(Krange,avgSE_MMSE_corr,'r-','LineWidth',1);
plot(Krange,avgSE_MR_uncorr,'b--','LineWidth',1);
plot(Krange,avgSE_MR_corr,'k-','LineWidth',1);

plot(Krange([1 5:5:50]),avgSE_MMSE_uncorr([1 5:5:50]),'bo','LineWidth',1);
plot(Krange([1 5:5:50]),avgSE_MMSE_corr([1 5:5:50]),'ro','LineWidth',1);
