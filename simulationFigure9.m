%This Matlab script can be used to generate Figure 9 in the article:
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

%Maximum difference between SNRs
pathlossDifference = 20;

%Number of user drops
nbrofScenarios = 1000;

%Number of channel realizations
nbrOfRealizationsOriginal = 200;

%Number of signal transmissions, used to compute Cee by Monte-Carlo methods
signalTransmissions = 10000;

%Number of channel models to be compared
nbrOfChannelModels = 3;

%Compute an identity matrix
I_M = eye(M);


%Prepare to save simulation results
SE_MMSE_corr = zeros(K,nbrofScenarios,nbrOfChannelModels);
SE_MMSE_uncorr = zeros(K,nbrofScenarios,nbrOfChannelModels);
SE_MR_corr = zeros(K,nbrofScenarios,nbrOfChannelModels);
SE_MR_uncorr = zeros(K,nbrofScenarios,nbrOfChannelModels);


%% Go through all user drops
for s = 1:nbrofScenarios
    
    %Write out the progress at every iterations
    disp(['Progress: ' num2str(s) ' out of ' num2str(nbrofScenarios) ' scenarios']);
    
    %Assign random signal-to-noise ratios of all users
    D_SNR = diag(db2pow(rand(K,1)*pathlossDifference-pathlossDifference/2));
    
    
    
    %Go through all channel models
    for model = 1:nbrOfChannelModels
        
        
        %Determine how many small-scale fading channel realizations are
        %needed, depending on if the channel model is fading or not
        if model == 1 || model == 2 || model == 3
            
            nbrOfRealizations = nbrOfRealizationsOriginal;
            
        elseif model == 4
            
            nbrOfRealizations = 1;
        end
        
        
        %Compute large-scale fading properties for the users
        if model == 2
            
            %Random user angles from -60 degrees to +60 degrees
            thetaValues = (2*pi/3)*rand(K,1)-pi/3;
            
            %Generate channel correlation matrices using the local
            %scattering model with 10 degree angular standard deviation
            R = zeros(M,M,K);
            
            for  k = 1:K
                
                R(:,:,k) = functionRlocalscatteringApprox(M,thetaValues(k),10,0.5);
                
            end
            
        elseif model == 3
            
            %Random user angles from -60 degrees to +60 degrees
            thetaValues = (2*pi/3)*rand(K,1)-pi/3;
            
        end
        
        
        %Go through all channel realizations
        for n = 1:nbrOfRealizations
            
            %Write out the progress at every iterations
            %disp(['Progress: ' num2str(n) ' out of ' num2str(nbrOfRealizations) ' realizations']);
            
            
            %Scaling factor that is used to make the real and imaginary part of
            %the received signal standard Gaussian random variables
            scaling = trace(D_SNR)/2*eye(M);
            
            %Generate channel realization
            if model == 1
                
                H = (randn(M,K)+1i*randn(M,K))/sqrt(2);
                
            elseif model == 2
                
                H = zeros(M,K);
                
                for k = 1:K
                    
                    H(:,k) = sqrtm(R(:,:,k))*(randn(M,1)+1i*randn(M,1))/sqrt(2);
                    
                end
                
            elseif model == 3
                
                H =  exp( repmat((0:M-1)',[1 K]) .* repmat(-2i*pi*0.5*sin(thetaValues'),[M 1]) );
                
            end
            
            %Generate the signals to be transmitted in the Monte-Carlo transmission
            S = (randn(K,signalTransmissions)+1i*randn(K,signalTransmissions))/sqrt(2);
            
            %Compute %Compute C_{uu} for the given channel realization
            Cuu = (H*D_SNR*H');
            
            %Compute the noise-free received signal
            U = H*sqrtm(D_SNR)*S;
            
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
            DH = D*H*sqrtm(D_SNR);
            
            
            %% DA-MMSE combining with correlated distortion
            V = ((DH*DH')+Cee+I_M)\DH;
            
            channelproducts = abs(DH'*V).^2;
            
            numerators = kappa*diag(channelproducts)';
            denominators = sum(channelproducts,1) - numerators + sum(abs(V).^2,1) + real(sum((Cee*V).*conj(V),1));
            
            SE_MMSE_corr(:,s,model) = SE_MMSE_corr(:,s,model) + (log2(1+numerators./denominators))'/nbrOfRealizations;
            
            
            %% DA-MR combining with correlated distortion
            V = DH;
            
            channelproducts = abs(DH'*V).^2;
            
            numerators = kappa*diag(channelproducts)';
            denominators = sum(channelproducts,1) - numerators + sum(abs(V).^2,1) + real(sum((Cee*V).*conj(V),1));
            
            SE_MR_corr(:,s,model) = SE_MR_corr(:,s,model) + (log2(1+numerators./denominators))'/nbrOfRealizations;
            
            
            %% DA-MMSE combining with uncorrelated distortion
            V = ((DH*DH')+Ceediag+I_M)\DH;
            
            channelproducts = abs(DH'*V).^2;
            
            numerators = kappa*diag(channelproducts)';
            denominators = sum(channelproducts,1) - numerators + sum(abs(V).^2,1) + real(sum((Ceediag*V).*conj(V),1));
            
            SE_MMSE_uncorr(:,s,model) = SE_MMSE_uncorr(:,s,model) + (log2(1+numerators./denominators))'/nbrOfRealizations;
            
            
            %% DA-MR combining with uncorrelated distortion
            V = DH;
            
            channelproducts = abs(DH'*V).^2;
            
            numerators = kappa*diag(channelproducts)';
            denominators = sum(channelproducts,1) - numerators + sum(abs(V).^2,1) + real(sum((Ceediag*V).*conj(V),1));
            
            SE_MR_uncorr(:,s,model) = SE_MR_uncorr(:,s,model) + (log2(1+numerators./denominators))'/nbrOfRealizations;
            
            
        end
        
    end
    
end


%Create vertical axis to CDF curves
yaxis = linspace(0,1,nbrofScenarios*K);


%Compute horizontal axis to CDF curves
SE_MMSE_uncorr_sorted = sort(reshape(SE_MMSE_uncorr,[nbrofScenarios*K nbrOfChannelModels]),1,'ascend');
SE_MMSE_corr_sorted = sort(reshape(SE_MMSE_corr,[nbrofScenarios*K nbrOfChannelModels]),1,'ascend');
SE_MR_uncorr_sorted = sort(reshape(SE_MR_uncorr,[nbrofScenarios*K nbrOfChannelModels]),1,'ascend');
SE_MR_corr_sorted = sort(reshape(SE_MR_corr,[nbrofScenarios*K nbrOfChannelModels]),1,'ascend');

markerlocation = round(linspace(1,nbrofScenarios*K,10));



%% Plot Figure 9a
figure;
hold on; box on;

plot(SE_MMSE_uncorr_sorted(1,1),yaxis(1),'bo--','LineWidth',1);
plot(SE_MMSE_corr_sorted(1,1),yaxis(1),'ro-','LineWidth',1);

plot(SE_MR_uncorr_sorted(:,1),yaxis,'b--','LineWidth',1);
plot(SE_MR_corr_sorted(:,1),yaxis,'k-','LineWidth',1);

plot(SE_MMSE_uncorr_sorted(:,1),yaxis,'b--','LineWidth',1);
plot(SE_MMSE_corr_sorted(:,1),yaxis,'r-','LineWidth',1);

plot(SE_MMSE_uncorr_sorted(markerlocation,1),yaxis(markerlocation),'bo','LineWidth',1);
plot(SE_MMSE_corr_sorted(markerlocation,1),yaxis(markerlocation),'ro','LineWidth',1);

xlabel('SE per UE [bit/s/Hz]','Interpreter','Latex');
ylabel('CDF','Interpreter','Latex');
legend({'DA-MMSE, uncorr','DA-MMSE, corr','DA-MR, uncorr','DA-MR, corr'},'Location','NorthWest','Interpreter','Latex');
xlim([0 7]);


% Plot Figure 9b
figure;
hold on; box on;

plot(SE_MMSE_uncorr_sorted(1,2),yaxis(1),'bo--','LineWidth',1);
plot(SE_MMSE_corr_sorted(1,2),yaxis(1),'ro-','LineWidth',1);

plot(SE_MR_uncorr_sorted(:,2),yaxis,'b--','LineWidth',1);
plot(SE_MR_corr_sorted(:,2),yaxis,'k-','LineWidth',1);

plot(SE_MMSE_uncorr_sorted(:,2),yaxis,'b--','LineWidth',1);
plot(SE_MMSE_corr_sorted(:,2),yaxis,'r-','LineWidth',1);

plot(SE_MMSE_uncorr_sorted(markerlocation,2),yaxis(markerlocation),'bo','LineWidth',1);
plot(SE_MMSE_corr_sorted(markerlocation,2),yaxis(markerlocation),'ro','LineWidth',1);

xlabel('SE per UE [bit/s/Hz]','Interpreter','Latex');
ylabel('CDF','Interpreter','Latex');
legend({'DA-MMSE, uncorr','DA-MMSE, corr','DA-MR, uncorr','DA-MR, corr'},'Location','NorthWest','Interpreter','Latex');
xlim([0 7]);


% Plot Figure 9c
figure;
hold on; box on;

plot(SE_MMSE_uncorr_sorted(1,3),yaxis(1),'bo--','LineWidth',1);
plot(SE_MMSE_corr_sorted(1,3),yaxis(1),'ro-','LineWidth',1);

plot(SE_MR_uncorr_sorted(:,3),yaxis,'b--','LineWidth',1);
plot(SE_MR_corr_sorted(:,3),yaxis,'k-','LineWidth',1);

plot(SE_MMSE_uncorr_sorted(:,3),yaxis,'b--','LineWidth',1);
plot(SE_MMSE_corr_sorted(:,3),yaxis,'r-','LineWidth',1);

plot(SE_MMSE_uncorr_sorted(markerlocation,3),yaxis(markerlocation),'bo','LineWidth',1);
plot(SE_MMSE_corr_sorted(markerlocation,3),yaxis(markerlocation),'ro','LineWidth',1);

xlabel('SE per UE [bit/s/Hz]','Interpreter','Latex');
ylabel('CDF','Interpreter','Latex');
legend({'DA-MMSE, uncorr','DA-MMSE, corr','DA-MR, uncorr','DA-MR, corr'},'Location','NorthWest','Interpreter','Latex');
xlim([0 7]);
