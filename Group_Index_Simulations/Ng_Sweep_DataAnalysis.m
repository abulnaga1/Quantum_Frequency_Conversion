clear all;
load('2020_10_20_AlGaAs_x0p176_Ngsweep.mat');

delta_ng = ng_946 - ng_1550;

%%
figure()
hold on;
title('x=0.176');
ylabel('n_g^9^4^6 - n_g^1^5^5^0');
xlabel('Device Area (um^2)');

for i = 1:length(AR)
    plot(Area*1e12,delta_ng(i,:),'Linewidth',2);
    legend_string(i) = strcat("AR = ", num2str(AR(i)));
end
legend(legend_string);
set(gca,'FontSize',14)
set(gca,'Box','on');
% ylim([0, 0.1]);