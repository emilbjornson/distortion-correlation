%This Matlab script can be used to generate Figure 2 in the article:
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

%Define range of input amplitudes
u = 0:0.01:1;

%Different values of the non-linearity parameter
alpha1 = 0.1340;
alpha2 = 1/3;

%Compute output of the amplifiers based on (19)
y1 = u - alpha1*u.^3;
y2 = u - alpha2*u.^3;


%% Plot Figure 2
figure;
hold on; box on;
plot(u,u,'r','LineWidth',1);
plot(u,y1,'b--','LineWidth',1);
plot(u,y2,'k-.','LineWidth',1);
xlabel('Normalized input $|u_m|$','Interpreter','Latex');
ylabel('Normalized output $|g_m(u_m)|$','Interpreter','Latex');
legend({'Linear amplifier','Third-order non-linear, $\alpha=0.134$','Third-order non-linear, $\alpha=1/3$'},'Location','NorthWest','Interpreter','Latex');
