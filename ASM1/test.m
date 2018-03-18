1;
clc; clear; close all;  


temp=120; % we may supply different T's here
coeff = gibbs_ch_coeff(temp);
syms T;
g_ch = coeff(1) + coeff(2)*T + coeff(3)*T^2 + coeff(4)*T^4 + coeff(5)*T*log(T) + coeff(6)/T;
s_ch = -1*diff(g_ch);
h = g_ch + T*s_ch;
c_ch = diff(h);

S_at_T = double(subs(s_ch,T,temp));
H_at_T = double(subs(h,T,temp));
C_at_T = double(subs(c_ch,T,temp));

fprintf('Note to self : Check UNITS!\n');
fprintf('At T = %d :\n',temp);
fprintf('The value of entropy : %d J/K.\n',S_at_T);
fprintf('The value of enthalpy : %d J/mol.\n',H_at_T);
fprintf('The value of heat capacity : %d J/K-mol.\n\n',C_at_T);


function f = gibbs_ch_coeff(temp)
    % array is in this order : T^0, T, T^2, T^4, T*ln(T), T^(-1)     

    if (temp>2000)
        f=zeros(6,1); % just a fallback
    elseif (temp>163)
        f=[-10195.860754; 690.949887637; -0.0007; 0; -118.47637; 590527];
    elseif (temp>43)
        f=[11622.647246; -59.537709263; 0.27565; 0; 15.74232; 0];
    elseif (temp>0)
        f=[11369.937746; -5.641259263; 0; -8.333*10^(-6); 0; 0];
    else
        f=zeros(6,1); % just a fallback
    end
end