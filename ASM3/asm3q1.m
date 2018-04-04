1;

% constants (no independent composition variable)
R=8.314;
T=1200;
x_Ni=0.6; 
x_Al=1-x_Ni;
a1=0.5;
a2=0.5;

% to find (one independent site fraction)
y_Al_1=0.1;
y_Al_2=(x_Al-a1*y_Al_1)/a2;
y_Ni_1=(1-y_Al_1);
y_Ni_2=(x_Ni-a1*y_Ni_1)/a2;

% since the number of independent site fractions is not the same as the 
% number of independent mole fractions, we must use G minimisation to get
% the site fractions.

G_SER_Al=(-1)*11278.4+188.684*T-31.7482*T*log(T)-1.231e+028*T^(-9);
G_SER_Ni=(-1)*5179.16+117.854*T-22.096*T*log(T)-0.0048407*T^2;
G_Al_Al=10083-4.813*T+G_SER_Al;
G_Ni_Ni=8715.08-3.556*T+G_SER_Ni;
G_Al_Ni=(-1)*56500-10.7*T+1.4975*T*log(T)+(0.5)*(G_SER_Al+G_SER_Ni);

G_ref=y_Al_1*y_Al_2*G_Al_Al+y_Ni_1*y_Ni_2*G_Ni_Ni+(y_Al_1*y_Ni_2+y_Ni_1*y_Al_2)*G_Al_Ni;

G_conf=R*T*(0.5)*(y_Al_1*log(y_Al_1)*y_Al_2*log(y_Al_2)+y_Ni_1*log(y_Ni_1)+y_Ni_2*log(y_Ni_2));

L0_AlNi_Al=(-1)*14225-5.625*T;
L1_AlNi_Al=0;

L0_AlNi_Ni=(-1)*22050;
L1_AlNi_Ni=1115;

L0_Al_AlNi=L0_AlNi_Al;
L1_Al_AlNi=L1_AlNi_Al;

L0_Ni_AlNi=L0_AlNi_Ni;
L1_Ni_AlNi=L1_AlNi_Ni;

L_AlNi_Al=L0_AlNi_Al+L1_AlNi_Al*(y_Al_1-y_Ni_1);
L_AlNi_Ni=L0_AlNi_Ni+L1_AlNi_Ni*(y_Al_1-y_Ni_1);

L_Al_AlNi=L0_Al_AlNi+L1_Al_AlNi*(y_Al_2-y_Ni_2);
L_Ni_AlNi=L0_Ni_AlNi+L1_Ni_AlNi*(y_Al_2-y_Ni_2);

G_xs=y_Al_1*y_Ni_1*(y_Al_2*L_AlNi_Al+y_Ni_2*L_AlNi_Ni)+y_Al_2*y_Ni_2*(y_Al_1*L_Al_AlNi+y_Ni_1*L_Ni_AlNi);

G=G_ref+G_conf+G_xs;