%This Matlab script can be used to generate Figure 12 in the article:
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

%Number of UEs
K = 5;

%Number of channel realizations
nbrOfRealizations = 400;

%Number of signal transmissions, used to compute Cee by Monte-Carlo methods
signalTransmissions = 10000;

%Range of signal-to-noise ratios (in dB)
SNRdBrange = -10:5:10;
SNRrange = db2pow(SNRdBrange);

%Create an identity matrix
I_M = eye(M);

%Create a pilot matrix based on DFT
pilotMatrix = fft(eye(K));


%Prepare to save simulation results
signal_perfect = zeros(length(SNRdBrange),1);
signal = zeros(length(SNRdBrange),1);

interference_perfect = zeros(length(SNRdBrange),1);
interference = zeros(length(SNRdBrange),1);

BSdistortion_corr_perfect = zeros(length(SNRdBrange),1);
BSdistortion_corr = zeros(length(SNRdBrange),1);

BSdistortion_uncorr_perfect = zeros(length(SNRdBrange),1);
BSdistortion_uncorr = zeros(length(SNRdBrange),1);



%% Go through all channel realizations
for n = 1:nbrOfRealizations
    
    %Write out the progress at every iterations
    disp(['Progress: ' num2str(n) ' out of ' num2str(nbrOfRealizations) ' realizations']);
    
    for k = 1:length(SNRdBrange)
        
        %Extract the SNR
        SNR = db2pow(SNRdBrange(k));
        
        %Scaling factor that is used to make the real and imaginary part of
        %the received signal standard Gaussian random variables
        scaling = (SNR*K)/2*eye(M);
        scalingPilot = (SNR*K+1)/2*eye(M);
        
        %Generate channel realization
        H = (randn(M,K)+1i*randn(M,K))/sqrt(2);
        
        %Generate noise realization and UE distortion in pilot transmission
        N = (randn(M,K)+1i*randn(M,K))/sqrt(2);
        UE_dist = (randn(K,K)+1i*randn(K,K))/sqrt(2);
        
        %Generate the signals to be transmitted in the Monte-Carlo transmission
        S = (randn(K,signalTransmissions)+1i*randn(K,signalTransmissions))/sqrt(2);
        
        %Compute %Compute C_{uu} for the given channel realization
        Cuu = SNR*(H*H');
        
        %Compute the noise-free received signal from dummy data
        U = sqrt(SNR)*H*S;
        
        %Compute the noise-free received signal from pilots
        U_pilot = sqrt(SNR)*H*(sqrt(kappa)*pilotMatrix+sqrt(1-kappa)*UE_dist)+N;
        
        
        %Normalize the received signal so that the real and imaginary parts
        %have unit variance, not for the given channel realization but on the
        %average
        U_normalized = sqrtm(scaling)\U;
        U_pilot_normalized = sqrtm(scalingPilot)\U_pilot;
        
        %Apply non-linear distortion
        U_nonlinear = U_normalized - alpha/(2*b_off)*abs(U_normalized).^2.*U_normalized;
        U_pilot_nonlinear = U_pilot_normalized - alpha/(2*b_off)*abs(U_pilot_normalized).^2.*U_pilot_normalized;
        
        
        %Compute the vector with diagonal elements of D using expression
        %given between (7) and (8)
        dvec = zeros(M,1);
        
        
        U_quant = zeros(M,signalTransmissions);
        U_pilot_quant = zeros(M,K);
        
        for l = 1:signalTransmissions
            
            %Quantization of real and imaginary parts of received signal
            real_index = quantiz(real(U_nonlinear(:,l)),partitions{ADCresolution});
            imag_index = quantiz(imag(U_nonlinear(:,l)),partitions{ADCresolution});
            
            
            %Reconstruct the quantized signal using the quantization levels
            U_quant(:,l) = sqrtm(scaling)*(codebooks{ADCresolution}(real_index+1)+1i*codebooks{ADCresolution}(imag_index+1));
            
            
            if l<=K
                
                %Quantization of real and imaginary parts of received pilot signal
                real_pilot_index = quantiz(real(U_pilot_nonlinear(:,l)),partitions{ADCresolution});
                imag_pilot_index = quantiz(imag(U_pilot_nonlinear(:,l)),partitions{ADCresolution});
                U_pilot_quant(:,l) = sqrtm(scalingPilot)*(codebooks{ADCresolution}(real_pilot_index+1)+1i*codebooks{ADCresolution}(imag_pilot_index+1));
                
            end
            
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
        
        
        %Compute LS channel estimates
        H_est = U_pilot_quant * pilotMatrix' / sqrt(K*SNR);
        DH_Est = D*H_est;
        
        
        %% DA-MR combining with perfect CSI
        V = DH./repmat(sqrt(sum(abs(DH).^2,1)),[M 1]);
        
        %Signal term
        signal_perfect(k) = signal_perfect(k) + SNR*mean(abs(diag(DH'*V)).^2)/nbrOfRealizations;
        
        %Interference term
        interference_perfect(k) = interference_perfect(k) + SNR*mean(sum(abs(DH'*V).^2,1)'-abs(diag(DH'*V)).^2)/nbrOfRealizations;

        %BS distortion
        BSdistortion_corr_perfect(k) = BSdistortion_corr_perfect(k) + mean(real(sum((Cee*V).*conj(V),1)))/nbrOfRealizations;
        BSdistortion_uncorr_perfect(k) = BSdistortion_uncorr_perfect(k) + mean(real(sum((Ceediag*V).*conj(V),1)))/nbrOfRealizations;
        
        
        %% DA-MR combining with imperfect CSI
        V = DH_Est./repmat(sqrt(sum(abs(DH_Est).^2,1)),[M 1]);
        
        %Signal term
        signal(k) = signal(k) + SNR*mean(abs(diag(DH'*V)).^2)/nbrOfRealizations;
        
        %Interference term
        interference(k) = interference(k) + SNR*mean(sum(abs(DH'*V).^2,1)'-abs(diag(DH'*V)).^2)/nbrOfRealizations;
        
        %BS distortion
        BSdistortion_corr(k) = BSdistortion_corr(k) + mean(real(sum((Cee*V).*conj(V),1)))/nbrOfRealizations;
        BSdistortion_uncorr(k) = BSdistortion_uncorr(k) + mean(real(sum((Ceediag*V).*conj(V),1)))/nbrOfRealizations;
        
        
    end
    
end


%Compute the SINRs
SINR_corr_perfect = kappa*signal_perfect./(BSdistortion_corr_perfect+interference_perfect+(1-kappa)*signal_perfect+1);
SINR_uncorr_perfect = kappa*signal_perfect./(BSdistortion_uncorr_perfect+interference_perfect+(1-kappa)*signal_perfect+1);
SINR_corr = kappa*signal./(BSdistortion_corr+interference+(1-kappa)*signal+1);
SINR_uncorr = kappa*signal./(BSdistortion_uncorr+interference+(1-kappa)*signal+1);



%% Plot simulation results
figure;
hold on; box on;

plot(SNRdBrange,10*log10(SINR_uncorr_perfect),'r-','LineWidth',1);
plot(SNRdBrange,10*log10(SINR_corr_perfect),'b--','LineWidth',1);
plot(SNRdBrange,10*log10(SINR_uncorr),'ro-','LineWidth',1);
plot(SNRdBrange,10*log10(SINR_corr),'bo--','LineWidth',1);

xlabel('SNR ($p/\sigma^2$)','Interpreter','Latex');
ylabel('SINR [dB]','Interpreter','Latex');
legend({'Perfect CSI, uncorr','Perfect CSI, corr','Imperfect CSI, uncorr','Imperfect CSI, corr'},'Interpreter','Latex','Location','SouthEast');
