clear all 

%% enter data to fit 

%%%% 1st mice
%xdata(:,1) is time and xdata(:,2) is the corresponding tumor number
fitData{1}.xdata(:,1) = [0;3;6;7;9;11;14;16;18;23;25;28;29;38;45;49;52;53];
fitData{1}.xdata(:,2) = [347.98;146.79;125.02;154.83;317.93;140.31;...
          44.98;62.5;108.58;87.81;222.43;116.64;339.21;418.97;...
          835.3;1353.21;2176;1633.51]/(1.2882*10^(-6)); 

%ydata(:,1) is time and ydata(:,2) is the corresponding T cell number
fitData{1}.ydata(:,1) = [0;3;7;10;14;18]; 
fitData{1}.ydata(:,2) = [3259100;8019100;11118000;3378700;3923800;5528400];

%% select fixed parameters 
 global a b theta_1 rho delta_E delta_R 

 a = 0.2253;
 b = 1/(4000/(1.2882*10^(-6)));

 theta_1 = 11;
 rho = 0.0016;
 delta_E = 0.0330;
 delta_R = 0.1;

global  beta_1  m  n theta_2 alpha_1 beta_2 alpha_2 gamma 
global  beta_3 alpha_3 delta_X R_0 X_0	

beta_1 = 0.10077; %*1.1
m = 6.1562E-10; %*0.9
n =4.1965E-07;  %*1.1;
theta_2 =1.4985E-07;
alpha_1 = 2.4625E+07;
beta_2 =0.10249; %*0.9
alpha_2 = 1.0132E+07;
gamma = 1.1002E+03;
beta_3 = 0.23858;
alpha_3 = 2.5556E+07;
delta_X =0.085333;
R_0 =0.017916;
X_0 = 0.50019; 

%% Simulations
%%%------------------ mice 1--------------------%%%
figure;
hold on
for i = 1 
   
    opts = odeset('AbsTol',1e-5);
    tumplottime = fitData{i}.xdata(:,1);
    tcelplottime = fitData{i}.ydata(:,1);

    CTint= [fitData{i}.xdata(1,2);fitData{i}.ydata(1,2)];
    % Note that total  Tcell = E + R + X
    EstmCTint = [CTint(1); CTint(2)*(1 - R_0- X_0);...
                CTint(2)*R_0; CTint(2)*X_0];


    [t,modFit] = ode45(@(t,y)fCE_final(t,y), tumplottime, EstmCTint,opts);

    tcellEndtime = size(tcelplottime);
    %%% effector T cell 
    eftcellEsm = modFit(:,2);
    eftcell = eftcellEsm(1:tcellEndtime);
    %%% Treg cell 
    tregcellEsm = modFit(:,3);
    tergcell = tregcellEsm(1:tcellEndtime);
    %%% exhausted T cell 
    extcellEsm = modFit(:,4);
    extcell = extcellEsm(1:tcellEndtime);
    %%% total T cell 
    tcellSum = eftcell + tergcell + extcell;
   

    subplot(1,2,1)
    hold on
    plot(tumplottime,modFit(:,1),'-r','LineWidth',3)
    plot(fitData{i}.xdata(:,1),fitData{i}.xdata(:,2),'or','MarkerSize',8,'MarkerFaceColor','r')
    xlabel('Time in Days','FontSize',18)
    ylabel('Cell Count','FontSize',18)
    legend({'C','data C'},'Location','best')
    set(gca,'FontSize',18)
    hold off
    
    subplot(1,2,2)
    hold on
    plot(tcelplottime,tcellSum,':b','LineWidth',3)
    plot(tcelplottime,eftcell,'LineWidth',3)
    plot(tcelplottime,tergcell,'LineWidth',3)
    plot(tcelplottime,extcell,'LineWidth',3)
    plot(fitData{i}.ydata(:,1),fitData{i}.ydata(:,2),'ob','MarkerSize',8,'MarkerFaceColor','b')
    legend({'Total T ','E', 'R','X','data T'},'Location','best')
    xlabel('Time in Days','FontSize',18)
    ylabel('Cell Count','FontSize',18)
    ylim([5e3 7e7])
    set(gca,'YScale', 'log')
    set(gca,'FontSize',18)
    hold off
   
end


sgt = sgtitle('Mouse 1: if the value of a decreased by 10%','Color','black');
sgt.FontSize = 20;

hold off


%%%----------------------------END------------------------------%%%


function value = fCE_final( t, y)  

global a b theta_1 delta_E delta_R rho 
global n beta_1 m theta_2 alpha_1 beta_2 alpha_2 gamma 
global beta_3 alpha_3 delta_X

   
C = y(1); % cancer cell
E = y(2); % effector CD8+ T cell
R = y(3); % treg  T cell
X = y(4); % exhausted CD8+ T cell



value(1,1) =  a*C*(1-b*C) - n*C*E*( 1/(1+theta_1*R/E) ); 

value(2,1) = -delta_E*E + beta_1*E*(C/(alpha_1+C))*(1/(1+theta_2*R)) - m*E*C - rho*E;
          
value(3,1) = -delta_R*R + beta_2*R*(C/(alpha_2+C))*(E/(gamma+E)) + rho*E;

value(4,1) = -delta_X*X + beta_3*X*(C/(alpha_3+C))*(1/(1+theta_2*R)) + m*E*C;


end    
