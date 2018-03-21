% Using Wang's Formula to get Muggianu, Colinet and Kohler Approximations
1;
c=4; % quaternary
% x=[0.25 0.25 0.25 0.25]; % order [Co; Cu; Fe; Ni]
% format = ^0 ^1 ^2 ^3
% coeff=zeros(6,4);
coeff=[ 39332 -1356 7953 -1119;
      -9312 -1752 0 0;
      1331 0 0 0;
      35626 -1530 12714 -1177;
      12049 -1862 0 0;
      -18379 -9228 0 0];

n=100;
z=60/n;
x=zeros(n,4);
for jj=1:n
    x(jj,1)=0.1;
    x(jj,2)=0.3;
    x(jj,3)=(60-jj*z)/100;
    x(jj,4)=(jj*z)/100;
end
h_E=zeros(size(x,1),1);

% The general schema :
for ii=1:size(x,1)
    for i=1:c-1
        for j=i+1:c
            for k=1:t_ij
                % disp(num2str(i));
                % disp(num2str(j));
                
                % scheme related values
                % % muggianu
                % t_ij=1;
                % beta_ij=1;
                % lambda_ij=0;
                
                % % Kohler
                % t_ij=1;
                % beta_ij=1;
                % lambda_ij=(x(ii,i)-x(ii,j))/(x(ii,i)+x(ii,j));
                
                % % Colinet
                t_ij=2;
                beta_ij=0.5;
                if(t_ij==1)
                    lambda_ij=1;
                else
                    lambda_ij=-1;
                end
                
                % disp(num2str(lambda_ij)); 
                [p,q,r]=parameters(x(ii,i),x(ii,j),lambda_ij);
                
                gE = p*q*( - 1356*(x_1-x_2) + 7953*(x_1-x_2)^2 - 1119*(x_1-x_2)^3);
                disp(num2str([p,q]));
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

function [x_ij,x_ji,f_ij]=parameters(x_1,x_2,lambda_ij)
x_ij=0.5*((1+x_1-x_2)+lambda_ij*(1-x_1-x_2));
x_ji=0.5*((1+x_2-x_1)+lambda_ij*(1-x_1-x_2));
f_ij=(x_1*x_2)/(x_ij*x_ji);
end