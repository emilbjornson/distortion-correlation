%This Matlab script can be used to generate Figure 3 in the article:
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
Krange = 1:50;


%Compute BS distortion with correlation, according to Eq. (30)
BSdistortion_corr = SNR*2*alpha^2*(Krange+6+9./Krange+4./Krange.^2 + 2*M*(Krange+1)./Krange.^2 )/b_off^2;

%Compute BS distortion without correlation, according to to Eq. (31)
BSdistortion_uncorr = SNR*2*alpha^2*(Krange+6+11./Krange+6./Krange.^2 )/b_off^2;

%Compute UE distortion according to Eq. (33)
UEdistortion = SNR*(1-kappa)*((M+1)-4*alpha*(M*Krange+Krange+M+3)./Krange/b_off+4*alpha^2*(M*Krange.^2+8*Krange+11+2*M*Krange+Krange.^2+M)./Krange.^2/b_off^2);



%% Plot Figure 3
figure;
hold on; box on;
plot(Krange,10*log10(UEdistortion),'b--','LineWidth',1);
plot(Krange,10*log10(BSdistortion_corr),'r','LineWidth',1);
plot(Krange,10*log10(BSdistortion_uncorr),'k-.','LineWidth',1);
xlabel('Number of UEs ($K$)','Interpreter','Latex');
ylabel('Distortion over noise [dB]','Interpreter','Latex');
legend({'UE distortion','BS distortion, correlated','BS distortion, uncorrelated'},'Interpreter','Latex');

%Create a shaded area between the curves for K between 1 and 10
y = [10*log10(BSdistortion_corr(1:10)); 10*(log10(BSdistortion_uncorr(1:10))-log10(BSdistortion_corr(1:10)))]';
h = area(Krange(1:10), y, -15);
set(h(1), 'FaceColor', 'none');
set(h, 'LineStyle', 'none');
h(2).FaceColor = [0.85 0.85 0.85];

%Plot the curves again to put them on top of the shaded area
plot(Krange,10*log10(UEdistortion),'b--','LineWidth',1);
plot(Krange,10*log10(BSdistortion_corr),'r','LineWidth',1);
plot(Krange,10*log10(BSdistortion_uncorr),'k-.','LineWidth',1);
