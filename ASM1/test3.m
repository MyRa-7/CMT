1;
clc; clear; close all;  

% material specific constants
Tc = 485;
f = 0.28;
bet = 1.008;

K_ferro = (266625/259)*(f-1)*log(1+bet)/(293*f - 790);
K_para = -1*(79875/74)*(f)*log(1+bet)/(293*f - 790);

temp=T_array(i); % we may supply different T's here
syms T;

coeff_ch = gibbs_ch_coeff(temp);
g_ch = coeff_ch(1) + coeff_ch(2)*T + coeff_ch(3)*T^2 + coeff_ch(4)*T^4 + coeff_ch(5)*T*log(T) + coeff_ch(6)/T;
s_ch = -1*diff(g_ch);
h_ch = g_ch + T*s_ch;
c_ch = diff(h_ch);
    
coeff_ph = c_ph_coeff(temp,Tc,K_ferro,K_para);
c_ph = coeff_ph(1)*T^3 + coeff_ph(2)*T^9 + coeff_ph(3)*T^15 + coeff_ph(4)*T^(-5) + coeff_ph(5)*T^(-15) + coeff_ph(6)*T^(-25);

