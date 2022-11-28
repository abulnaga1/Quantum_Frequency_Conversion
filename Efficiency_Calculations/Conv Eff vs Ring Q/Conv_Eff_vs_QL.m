clear all; clc; close all;

incl_3order = 1; % choose 1 if the 3rd order is included, i.e. omega_m = omega_o + D_1*mu + 1/2*D_2*mu^2 + 1/6*D_3*mu^3

L = 2*pi*(25e-6); % cavity round-trip length
D1_946 = 3.3334e12;     % FSR at 946nm
D1_1550 = 3.3342e12;    % FSR at 1550nm
D2_946 = -7.6336e8;
D2_1550 = -2.775e8;
D3_946 = -6.0565e6;
D3_1550 = 8.2346e7;

tR_946 = 2*pi/D1_946; % cavity round-trip time at 946nm
tR_1550 = 2*pi/D1_1550; % cavity round-trip time at 1550nm

h_bar = 1.0546e-34; % Planck's constant
c_const = physconst('lightspeed'); % speed of light

lambda_p1 = 946.6e-9;
lambda_p2 = 1547.8e-9;
cavityomega_p1 = 2*pi*c_const/lambda_p1; % cavity resonance near pump 1 at 946nm
cavityomega_p2 = 2*pi*c_const/lambda_p2; % cavity resonance near pump 2 at 1550nm
P_s = 0.1e-9; % power of the signal
E_in = [sqrt(P_s);0;0;0]; % input driving field for frequency downconversion

gamma_p1 = 287.697; % nonlinear coefficient for pump 1 at 946nm (NO no/ng PREFACTOR)
gamma_p2 = 152.5; % nonlinear coefficient for pump 2 at 1550nm (NO no/ng PREFACTOR)

mu = 1; % signal separation from pump

% Qc_array = (1:0.25:10)*1e4;
Qc_array = linspace(0.5e-4,20e4,200);
P_tot_array = [10 30 50]*1e-3;

i_plus_flux = zeros(length(Qc_array),length(P_tot_array));
i_minus_flux = zeros(length(Qc_array),length(P_tot_array));
lambda_i_plus = zeros(length(Qc_array),length(P_tot_array));
lambda_i_minus = zeros(length(Qc_array),length(P_tot_array));

for m = 1:length(P_tot_array)
    P_p1 = (1/2)*P_tot_array(m); % power of pump 1 at 946nm
    P_p2 = (1/2)*P_tot_array(m); % power of pump 2 at 1550nm

    for j=1:length(Qc_array)
        Qc_946 = Qc_array(j);%Coupling Q at 946nm
        Qc_1550 = Qc_946;    %Coupling Q at 1550nm
        QL_946 = Qc_946/2;   % loaded Q at 946nm
        QL_1550 = Qc_1550/2; % loaded Q at 1550nm

        alpha_p1 = cavityomega_p1*tR_946/(2*QL_946); % total loss rate for pump 1 at 946nm
        alpha_p2 = cavityomega_p2*tR_1550/(2*QL_1550); % total loss rate for pump 2 at 1550nm
        theta_p1 = cavityomega_p1*tR_946/Qc_946; %waveguide power coupling rate for pump 1 at 946nm
        theta_p2 = cavityomega_p2*tR_1550/Qc_1550; %waveguide power coupling rate for pump 2 at 1550nm
        theta_array = [theta_p1 theta_p1 theta_p2 theta_p2];
        theta_matrix = diag(theta_array);

        detuning_array = [-5:0.1:5]'*alpha_p1;

        % for frequency downconversion
        cavityomega_s = cavityomega_p1 + D1_946*mu + (1/2)*D2_946*(mu^2) + incl_3order*(1/6)*D3_946*(mu^3); % cavity resonance near signal
        cavityomega_b = cavityomega_p1 - D1_946*mu + (1/2)*D2_946*(mu^2) - incl_3order*(1/6)*D3_946*(mu^3); % cavity resonance near auxiliary tone
        cavityomega_c = cavityomega_p2 + D1_1550*abs(mu) + (1/2)*D2_1550*(abs(mu)^2) + incl_3order*(1/6)*D3_1550*(abs(mu)^3); % cavity resonance near idler plus
        cavityomega_d = cavityomega_p2 - D1_1550*abs(mu) + (1/2)*D2_1550*(abs(mu)^2) - incl_3order*(1/6)*D3_1550*(abs(mu)^3); % cavity resonance near idler minus

        deltaphi_p1 = 0; % effective detuning for pump 1 at 946nm
        deltaphi_p2 = 0;%-0.1*alpha_p2; % effective detuning for pump 2 at 1550nm

        E_p1 = 1i*sqrt(theta_p1*P_p1)/(alpha_p1+1i*deltaphi_p1); % intracavity electric field for pump 1 at 946nm
        E_p2 = 1i*sqrt(theta_p2*P_p2)/(alpha_p2+1i*deltaphi_p2); % intracavity electric field for pump 2 at 1550nm

        omega_p1 = cavityomega_p1 - (deltaphi_p1+gamma_p1*L*(E_p1*conj(E_p1)+2*E_p2*conj(E_p2)))/tR_946; % actual frequency for pump 1 at 946nm 
        omega_p2 = cavityomega_p2 - (deltaphi_p2+gamma_p2*L*(2*E_p1*conj(E_p1)+E_p2*conj(E_p2)))/tR_1550; % actual frequency for pump 2 at 1550nm

        for k=1:length(detuning_array)
            deltaphi_s = detuning_array(k); % effective detuning for signal
    %         deltaphi_b = (cavityomega_s+cavityomega_b-2*cavityomega_p1)*tR_946...
    %             - deltaphi_s + 2*deltaphi_p1 - 2*gamma_p1*L*E_p1*conj(E_p1); % effective detuning for auxiliary tone
    %         deltaphi_c = ((cavityomega_c-cavityomega_p2) - (cavityomega_s-cavityomega_p1))*tR_1550...
    %             + deltaphi_s*tR_1550/tR_946 - (deltaphi_p1*tR_1550/tR_946-deltaphi_p2)...
    %             + L*(gamma_p1*E_p1*conj(E_p1)*tR_1550/tR_946-gamma_p2*E_p2*conj(E_p2)); % effective detuning for idler plus 
    %         deltaphi_d = ((cavityomega_d-cavityomega_p2) + (cavityomega_s-cavityomega_p1))*tR_1550...
    %             - deltaphi_s*tR_1550/tR_946 + (deltaphi_p1*tR_1550/tR_946+deltaphi_p2)...
    %             - L*(gamma_p1*E_p1*conj(E_p1)+gamma_p2*E_p2*conj(E_p2)); % effective detuning for idler minus

            omega_s = cavityomega_s - (deltaphi_s+2*gamma_p1*L*(E_p1*conj(E_p1)+E_p2*conj(E_p2)))/tR_946; % actual frequency for signal
            omega_b = 2*omega_p1 - omega_s; % actual frequency for auxiliary tone
            omega_i_plus = omega_p2 + abs(omega_p1 - omega_s); % actual frequency for idler plus
            omega_i_minus = omega_p2 - abs(omega_p1 - omega_s); % actual frequency for idle minus

            deltaphi_b = (cavityomega_b - omega_b)*tR_946 - 2*gamma_p1*L*(E_p1*conj(E_p1)+E_p2*conj(E_p2));
            deltaphi_c = (cavityomega_c - omega_i_plus)*tR_1550 - 2*gamma_p2*L*(E_p1*conj(E_p1)+E_p2*conj(E_p2));
            deltaphi_d = (cavityomega_d - omega_i_minus)*tR_1550 - 2*gamma_p2*L*(E_p1*conj(E_p1)+E_p2*conj(E_p2));

            flux_matrix = (1/h_bar)./[omega_s;omega_b;omega_i_plus;omega_i_minus];
            in_flux = P_s/(h_bar*omega_s);

            M_matrix = zeros(4,4);
            M_matrix(1,:) = [-alpha_p1-1i*deltaphi_s,...
                1i*gamma_p1*L*E_p1*E_p1,...
                1i*2*gamma_p1*L*E_p1*conj(E_p2),...
                1i*2*gamma_p1*L*E_p1*E_p2];
            M_matrix(2,:) = [-1i*gamma_p1*L*conj(E_p1)*conj(E_p1),...
                -alpha_p1+1i*deltaphi_b,...
                -1i*2*gamma_p1*L*conj(E_p1)*conj(E_p2),...
                -1i*2*gamma_p1*L*conj(E_p1)*E_p2];
            M_matrix(3,:) = [1i*2*gamma_p2*L*conj(E_p1)*E_p2,...
                1i*2*gamma_p2*L*E_p1*E_p2,...
                -alpha_p2-1i*deltaphi_c,...
                1i*gamma_p2*L*E_p2*E_p2];
            M_matrix(4,:) = [-1i*2*gamma_p2*L*conj(E_p1)*conj(E_p2),...
                -1i*2*gamma_p2*L*E_p1*conj(E_p2),...
                -1i*gamma_p2*L*conj(E_p2)*conj(E_p2),...
                -alpha_p2+1i*deltaphi_d];

            E_r = M_matrix\(-1i*sqrt(theta_matrix)*E_in);

            E_WG = E_in + 1i*sqrt(theta_matrix)*E_r;
            P_WG = E_WG.*conj(E_WG);

            out_flux = P_WG.*flux_matrix;

            if (out_flux(3)/in_flux) > i_plus_flux(j,m)
                i_plus_flux(j,m) = out_flux(3)/in_flux;
                lambda_i_plus(j,m) = 2*pi*c_const/omega_i_plus;
            end
            if (out_flux(4)/in_flux) > i_minus_flux(j,m)
                i_minus_flux(j,m) = out_flux(4)/in_flux;
                lambda_i_minus(j,m) = 2*pi*c_const/omega_i_minus;
            end
        end
    end
end
colors = ["#EDB120", "#7E2F8E", "#0072BD"];
figure;
for m = 1:length(P_tot_array)
    legend_string(m) = strcat(num2str(P_tot_array(m)*1e3), "mW");
    semilogx(Qc_array/2,10*log10(i_plus_flux(:,m)),'LineWidth',2, 'Color', colors(m))
    hold on;
end
set(gca,'FontSize',16)
xlim([1e4 1e5])
ylim([-38 0])
set(gca,'YTick',[-30, -20, -10, 0]);
% x = xticks;
% xticks([10000 20000 30000 40000 50000])
xlabel('Q_L') 
% ylabel('Conversion efficiency (dB)')
xline(3.75e4,'k--');
legend(legend_string,'Location','Southeast');

