% Using Wang's Formula to get Muggianu, Colinet and Kohler Approximations
1;
c=4; % quaternary
x=[0.25 0.25 0.25 0.25]; % order [Co; Cu; Fe; Ni]
% scheme related values
% % muggianu
% t_ij=1;
% beta_ij=1;

% % Kohler
% t_ij=1;
% beta_ij=1;

% % Colinet
t_ij=2;
beta_ij=0.5;


n=100;
z=60/n;
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
                %           disp(num2str(i));
                %           disp(num2str(j));
                
                % % muggianu
                % lambda_ij=0;
                % % Kohler
                % lambda_ij=(x(ii,i)-x(ii,j))/(x(ii,i)+x(ii,j));
                % disp(num2str(lambda_ij));
                % % Colinet
                if(t_ij==1)
                    lambda_ij=1;
                else
                    lambda_ij=-1;
                end
                
                [p,q,r]=parameters(x(ii,i),x(ii,j),lambda_ij);
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
h = x_1*x_2*(35626 - 1530*(x_1-x_2) + 12714*(x_1-x_2)^2 - 1177*(x_1-x_2)^3);
end

function h = fCuNi(x_1,x_2)
h = x_1*x_2*(12049 - 1862*(x_1-x_2));
end

function h = fFeNi(x_1,x_2)
h = x_1*x_2*(-18379 - 9228*(x_1-x_2));
end