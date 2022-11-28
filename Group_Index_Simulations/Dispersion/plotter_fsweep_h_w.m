clear all;
c = physconst('lightspeed');
load('2020_10_23_filenames.mat');

%Create empty matrices to store all data
f_mat = zeros(length(fnames),freqPnt);
f_vg_mat = f_mat;
f_D_mat = f_mat;

neff_mat = f_mat;
vg_mat = f_mat;
D_mat = f_mat;

for k = 1:length(fnames)
    datafile = fullfile("data",fnames(k));
    data = importdata(strcat(datafile,".mat"));
    
    f_mat(k,:) = data.f;
    f_vg_mat(k,:) = data.f_vg;
    f_D_mat(k,:) = data.f_D;
    
    neff_mat(k,:) = data.neff;
    vg_mat(k,:) = data.vg;
    D_mat(k,:) = data.D;
end

%%
%Plot neff for the different widths
n = 1;
for i = 1:length(w_mat)
    figure()
    hold on;
    ylabel('neff');
    xlabel('\lambda (nm)');
    title(strcat("Width = ", num2str(w_mat(i).*1e9), "nm"));
    for j = 1:length(h_mat)
        plot(1e9*c./f_mat(n,:), neff_mat(n,:), 'Linewidth', 2);
        legend_string_w(j) = strcat("H=",num2str(h_mat(j).*1e9), "nm");
        n = n+1;
    end
    xlim([850 1650]);
    xline(946,'k--');
    xline(1550,'k--');
    legend(legend_string_w, "Location", "Northeast")
    set(gca,'FontSize',14)
    set(gca,'Box','on');
    saveas(gcf,strcat("neff_w_",num2str(w_mat(i)*1e9),"_h_Sweep",".png"));
end

%%
%Plot Vg, sweeping area for different widths, sweeping over height
n = 1;
for i = 1:length(w_mat)
    figure()
    hold on;
    ylabel('n_g');
    xlabel('\lambda (nm)');
    title(strcat("Width = ", num2str(w_mat(i).*1e9), "nm"));
    for j = 1:length(h_mat)
        plot(1e9*c./f_vg_mat(n,:), c./vg_mat(n,:), 'Linewidth', 2);
        legend_string_w(j) = strcat("H=",num2str(h_mat(j).*1e9), "nm");
        n = n+1;
    end
%     ylim([0.25 0.325]);
    xlim([850 1650]);
    set(gca,'FontSize',14);
    set(gca,'Box','on');
    xline(946,'k--');
    xline(1550,'k--');
    legend(legend_string_w, "Location", "Northeast")
    saveas(gcf,strcat("vg_w_",num2str(w_mat(i)*1e9),"_h_Sweep",".png"));
end

%%
%Plot D, sweeping area for fixed aspect ratios
n = 1;
for i = 1:length(w_mat)
    figure()
    hold on;
    ylabel('D [ps / (nm km)]');
    xlabel('\lambda (nm)');
    title(strcat("Width = ", num2str(w_mat(i).*1e9), "nm"));
    for j = 1:length(h_mat)
        plot(1e9*c./f_D_mat(n,:), D_mat(n,:)*1e6, 'Linewidth', 2);
        legend_string_w(j) = strcat("H=",num2str(h_mat(j).*1e9), "nm");
        n = n+1;
    end
    yline(0, 'k--');
    xlim([850 1650]);
    ylim([-4e3 2e3]);
    set(gca,'FontSize',14);
    set(gca,'Box','on');
    xline(946,'k--');
    xline(1550,'k--');
    legend(legend_string_w, "Location", "Northwest")
    saveas(gcf,strcat("D_w_",num2str(w_mat(i)*1e9),"_h_Sweep",".png"));
end