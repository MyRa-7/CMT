1;
clc; clear; close all;

Tc = 485;
f = 0.28;
bet = 1.008;

K_ferro = (266625/259)*(f-1)*log(1+bet)/(293*f - 790);
K_para = -1*(79875/74)*(f)*log(1+bet)/(293*f - 790);

temp=Tc; % we may supply different T's here
syms Tau;

dTau = (Tc);
% Tau = (temp/Tc);

coeff_ph = c_ph_coeff((temp/Tc),K_ferro,K_para);
c_ph = coeff_ph(1)*Tau^3 + coeff_ph(2)*Tau^9 + coeff_ph(3)*Tau^15 + coeff_ph(4)*Tau^(-5) + coeff_ph(5)*Tau^(-15) + coeff_ph(6)*Tau^(-25);
s_ph = int(c_ph/Tau);
h_ph = int(c_ph)*dTau;
g_ph = (h_ph - Tau*s_ph);

G_at_T = double(subs(g_ph,Tau,(temp/Tc)));
fprintf('The value of gibbs energy magnetic contribution : %d J/mol.\n',G_at_T);

function eff = c_ph_coeff(hom_temp,K_ferro,K_para)
    % array is in this order : T^3, T^9, T^15, T^(-5), T^(-15), T^(-25)
    R = 8.314; % in SI Units
    
    if (hom_temp>1)
        eff=K_para*R*2*[0; 0; 0; 1; (1/3); (1/5)];
    elseif (hom_temp>0)
        eff=K_ferro*R*2*[1; (1/3); (1/5); 0; 0; 0];
    else
        eff=zeros(6,1); % just a fallback
    end

    % normalising because Tc is a constant and must not disturb integration
%     eff(1)=eff(1)/Tc^3; eff(2)=eff(2)/Tc^9; eff(3)=eff(3)/Tc^15;
%     eff(4)=eff(4)*Tc^5; eff(5)=eff(5)*Tc^15; eff(6)=eff(6)*Tc^25;
end