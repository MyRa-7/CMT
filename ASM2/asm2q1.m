% Using Wang's Formula to get Muggianu, Colinet and Kohler Extrapolations
1;
clc; clear; close all;

c=4; % quaternary

n=100;
z=60/n;
% order [Co; Cu; Fe; Ni]
x=zeros(n,4);
for jj=1:n
    x(jj,1)=0.1;
    x(jj,2)=0.3;
    x(jj,3)=(60-jj*z)/100;
    x(jj,4)=(jj*z)/100;
end

fprintf('Part (a)\n\n');
h_e_muggianu=schemeswitch('m',c,x);
h_e_kohler=schemeswitch('k',c,x);
h_e_colinet=schemeswitch('c',c,x);

figure(1);
hold on;
plot(x(:,4),h_e_muggianu,'g-o',"linewidth",1.5);
plot(x(:,4),h_e_kohler,'b-+',"linewidth",1.5);
plot(x(:,4),h_e_colinet,'r-x',"linewidth",1.5);
title("Plot of \Delta_{mix} H_m (J-mol^{-1})","FontSize",18);
xlabel("Mole Fraction of Ni","FontSize",14);
ylabel("\Delta_{mix} H_m (J-mol^{-1})","FontSize",14);
axis('square');
set(gca,'FontSize',16);
legend("Muggianu","Kohler","Colinet","location","NorthEast");

fprintf('Part (b)\n');
x=0.25*ones(1,c);
h_e_muggianu=schemeswitch('m',c,x);
h_e_kohler=schemeswitch('k',c,x);
h_e_colinet=schemeswitch('c',c,x);
fprintf('The value using Muggianu Scheme : %d J/mol.\n',h_e_muggianu);
fprintf('The value using Kohler Scheme : %d J/mol.\n',h_e_kohler);
fprintf('The value using Colinet Scheme : %d J/mol.\n',h_e_colinet);

function h_E = schemeswitch(scheme,c,x)
% getting scheme related values
switch scheme
    case 'm' % Muggianu
        t_ij=1;
        beta_ij=1;
    case 'k' % Kohler
        t_ij=1;
        beta_ij=1;
    case 'c' % Colinet
        t_ij=2;
        beta_ij=0.5; % since it's the same for both
end
h_E=zeros(size(x,1),1);
% The general schema :
for ii=1:size(x,1)
    for i=1:c-1
        for j=i+1:c
            for k=1:t_ij
                % lambda also depends on scheme
                % % muggianu
                switch scheme
                    case 'm' % Muggianu
                        lambda_ij=0;
                        lambda_ji=0;
                    case 'k' % Kohler
                        lambda_ij=(x(ii,i)-x(ii,j))/(x(ii,i)+x(ii,j));
                        lambda_ji=(x(ii,j)-x(ii,i))/(x(ii,i)+x(ii,j));
                    case 'c' % Colinet
                        if(t_ij==1)
                            lambda_ij=1;
                            lambda_ji=1;
                        else
                            lambda_ij=-1;
                            lambda_ji=-1;
                        end
                end
                
                [p,q,r]=parameters(x(ii,i),x(ii,j),lambda_ij,lambda_ji);
                
                if(i==1)
                    if(j==2)
                        gE=fCoCu(p,q);
                    elseif(j==3)
                        gE=fCoFe(p,q);
                    else
                        gE=fCoNi(p,q);
                    end
                elseif(i==2)
                    if(j==3)
                        gE=fCuFe(p,q);
                    else
                        gE=fCuNi(p,q);
                    end
                else
                    gE=fFeNi(p,q);
                end
                h_E(ii)=h_E(ii)+beta_ij*r*gE;
            end
        end
    end
end
end

function [x_ij,x_ji,f_ij]=parameters(x_1,x_2,lambda_ij,lambda_ji)
x_ij=0.5*((1+x_1-x_2)+lambda_ij*(1-x_1-x_2));
x_ji=0.5*((1+x_2-x_1)+lambda_ji*(1-x_1-x_2));
f_ij=(x_1*x_2)/(x_ij*x_ji);
end

function h = fCoCu(x_1,x_2)
h = x_1*x_2*(39332 - 1356*(x_1-x_2) + 7953*(x_1-x_2)^2 - 1119*(x_1-x_2)^3);
end

function h = fCoFe(x_1,x_2)
h = x_1*x_2*(-9312 - 1752*(x_1-x_2));
end

function h = fCoNi(x_1,x_2)
h = x_1*x_2*(1331);
end

function h = fCuFe(x_1,x_2)
h = x_1*x_2*(35626 - 1530*(x_1-x_2) + 12714*(x_1-x_2)^2 + 1177*(x_1-x_2)^3);
end

function h = fCuNi(x_1,x_2)
h = x_1*x_2*(12049 - 1862*(x_1-x_2));
end

function h = fFeNi(x_1,x_2)
h = x_1*x_2*(-18379 + 9228*(x_1-x_2));
end