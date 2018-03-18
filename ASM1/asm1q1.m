% a symbolic maths attempt
1;
clc; clear; close all;

% material specific constants
Tc = 485;
f = 0.28;
bet = 1.008;

K_ferro = (266625/259)*(f-1)*log(1+bet)/(293*f - 790);
K_para = -1*(79875/74)*(f)*log(1+bet)/(293*f - 790);

fprintf('Part (a)\n\n');
T_array = [200, 400, 500];
calculate_properties(T_array,Tc,K_ferro,K_para,1);

fprintf('\nPart (b)\n\n');
T_array=5:5:1500;
C = calculate_properties(T_array,Tc,K_ferro,K_para,2);

fprintf('Part (c)\n\n');
T_array=Tc;
calculate_properties(T_array,Tc,K_ferro,K_para,3);

function C = calculate_properties(T_array,Tc,K_ferro,K_para,mode)

if (mode==2)
    C=zeros(size(T_array));
    C_mag=zeros(size(T_array));
end
maxT=-1; % fake lower limit, can be made as low as necessary.
% times = 0;
for i=1:size(T_array,2)
    temp=T_array(i); % we may supply different T's here
    syms T;
    % SER values from a handbook
%     H_ser_Fe=4489; 
%     H_ser_C=1054;
%     H_ser=3*H_ser_Fe+H_ser_C;
H_ser=0; % we are writing wrt SER
    % dTau = (Tc);
    
    if (temp>maxT) % assumes the T_array is in ascending order
        [coeff_ch, coeff_ph, para, maxT] = coeff(temp,Tc,K_ferro,K_para);
        
        g_ch = coeff_ch(1) + coeff_ch(2)*T + coeff_ch(3)*T^2 + ...
            coeff_ch(4)*T^4 + coeff_ch(5)*T*log(T) + coeff_ch(6)/T;
        s_ch = -1*diff(g_ch);
        h_ch = g_ch + T*s_ch + H_ser;
        c_ch = diff(h_ch);
 
        c_ph = coeff_ph(1)*T^3 + coeff_ph(2)*T^9 + coeff_ph(3)*T^15 + ...
            coeff_ph(4)*T^(-5) + coeff_ph(5)*T^(-15) + coeff_ph(6)*T^(-25);
        para_part = para(1)*T^3 + para(2)*T^9 + para(3)*T^15 + ...
            para(4)*T^(-5) + para(5)*T^(-15) + para(6)*T^(-25);
        
        if(temp>Tc)
            s_ph = int(c_ph/T);
        elseif (temp<=Tc)
            s_ph_para=int(para_part/T,T,Tc, inf);
            s_ph_ferro=int(c_ph/T,T,0,Tc);
            s_ph_total=s_ph_para+s_ph_ferro;
            s_ph=(int(c_ph/T)-s_ph_total);  % need to check signs    
        end
        
        if(temp>Tc)
            h_ph = int(c_ph);
        elseif (temp<=Tc)
            h_ph_para=int(para_part,T,inf,Tc);
            h_ph_ferro=int(c_ph,T,Tc,0);
            h_ph_total=h_ph_para+h_ph_ferro;
            h_ph=int(c_ph)+h_ph_total;  % need to check signs
        end
        
        g_ph = h_ph - T*s_ph;
%         
%         times=times+1;
%         disp(num2str(times));
        c = c_ph + c_ch; 
        h = h_ph + h_ch;
        s = s_ph + s_ch; 
        
        c_func=matlabFunction(c);
        h_func=matlabFunction(h);
        s_func=matlabFunction(s);
        g_ph_func=matlabFunction(g_ph);
        c_ph_func=matlabFunction(c_ph);
        
    end
   
    if (mode==1) 
        H_at_T = h_func(temp);
        C_at_T = c_func(temp);
        S_at_T = s_func(temp);

        %fprintf('Note to self : Check UNITS!\n');
        fprintf('At T = %d :\n',temp);
        fprintf('The value of enthalpy : %d J/mfu.\n',H_at_T);
        fprintf('The value of entropy : %d J/K.\n',S_at_T);
        fprintf('The value of heat capacity : %d J/K-mfu.\n\n',C_at_T);   
    elseif (mode==2)
       
        C(i)= c_func(temp);
        C_mag(i)=c_ph_func(temp);
        
    elseif (mode==3)
        
        G_mo_at_T=g_ph_func(temp);
%         G_ch_at_T=g_ch_func(temp);
        fprintf('The value of gibbs energy magnetic contribution : %d J/mol.\n',G_mo_at_T);
%         fprintf('The value of gibbs energy lattice contribution : %d J/mol.\n',G_ch_at_T);
    end
end

if(mode==2)
%     figure(1);
%     plot(T_array,C,'-o','linewidth',1.5);
%     figure(2);
%     plot(T_array,C_mag,"-o","linewidth",1.5);
%     
    C = [transpose(5:5:1500) transpose(C) ];    
    figure(1); hold on;
    plot(C(1:8,1),C(1:8,2),'-','linewidth',3);
    plot(C(8:33,1),C(8:33,2),'-','linewidth',3);
    plot(C(33:96,1),C(33:96,2),'-','linewidth',3);
    plot(C(98:300,1),C(98:300,2),'-','linewidth',3);
    title('Plot of C_p vs T (T=0 to T=1500 K)');
    legend('T=0 to T=43 K','T=44 to T=163 K','T=164 to T=485 K','T=486 to T=1500 K',"location","southeast");
    xlabel('Temperature (T) in Kelvin');
    ylabel('C_p in J/K-mfu');
    axis('square');
end

end

% Geting Functions for a range
% The gibbs function is a function of time, 
% the c function is a function of normalised temperature
function [g_eff_ch, c_eff_ph, para, maxT] = coeff(temp,Tc,K_ferro,K_para)
[g_eff_ch, maxT_ch] = gibbs_ch_coeff(temp);
[c_eff_ph, para, maxT_ph] = c_ph_coeff(temp,Tc,K_ferro,K_para);
maxT = min(maxT_ch,maxT_ph);
end

% Definition of coefficients of Gibbs Energy's Lattice Contribution
function [eff, maxT] = gibbs_ch_coeff(temp)
% array is in this order : T^0, T, T^2, T^4, T*ln(T), T^(-1)

if (temp>2000)
    eff=zeros(6,1); % just a fallback
    maxT=1e5; % fake max limit. Can me made as high as necessary.
elseif (temp>163)
    eff=[-10195.860754; 690.949887637; -0.0007; 0; -118.47637; 590527];
    maxT=2000;
elseif (temp>43)
    eff=[11622.647246; -59.537709263; -0.27565; 0; 15.74232; 0];
    maxT=163;
elseif (temp>0)
    eff=[11369.937746; -5.641259263; 0; -8.333*10^(-6); 0; 0];
    maxT=43;
else
    eff=zeros(6,1); % just a fallback
    maxT=0;
end
end

% Definition of coefficients of Heat Capacity's Physical Contribution
function [eff,eff_para, maxT] = c_ph_coeff(temp,Tc,K_ferro,K_para)
% array is in this order : T^3, T^9, T^15, T^(-5), T^(-15), T^(-25)

R = 8.314; % in SI Units
if (temp>Tc)
    eff=K_para*R*2*[0; 0; 0; 1; (1/3); (1/5)];
    eff_para=zeros(6,1); % we don't have the para contribution
    maxT=1e5; % fake max limit. Can be made as high as necessary.
elseif (temp>0)
    eff=K_ferro*R*2*[1; (1/3); (1/5); 0; 0; 0];
    eff_para=K_para*R*2*[0; 0; 0; 1; (1/3); (1/5)];
    maxT=Tc;
else
    eff=zeros(6,1); % just a fallback
    maxT=0;
end

% % normalising because Tc is a constant and must not disturb integration
eff(1)=eff(1)/Tc^3; eff(2)=eff(2)/Tc^9; eff(3)=eff(3)/Tc^15;
eff(4)=eff(4)*Tc^5; eff(5)=eff(5)*Tc^15; eff(6)=eff(6)*Tc^25;

eff_para(1)=eff_para(1)/Tc^3; eff_para(2)=eff_para(2)/Tc^9; eff_para(3)=eff_para(3)/Tc^15;
eff_para(4)=eff_para(4)*Tc^5; eff_para(5)=eff_para(5)*Tc^15; eff_para(6)=eff_para(6)*Tc^25;
end