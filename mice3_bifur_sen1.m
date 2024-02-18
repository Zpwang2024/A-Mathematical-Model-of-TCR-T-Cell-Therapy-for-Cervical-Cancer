clear all 

%% enter data to fit 

%%%% 3rd mice
% tumor cells
fitData{3}.xdata(:,1) = [0;3;6;7;9;11;14;16;18;23;25;28;29;38;45;49;52;53];
                        % time points
fitData{3}.xdata(:,2) = [58.77;267.57;151.1;93.23;258.57;283.89;...
          109.91;131.07;130.54;154.87;347.94;201.07;118.38;237.28;...
          665.8;733.67;1744.19;1529.78]/(1.2882*10^(-6)); % cell numbers

% total T cells
fitData{3}.ydata(:,1) = [0;3;7;18]; % time points
fitData{3}.ydata(:,2) = [628320;2678600;5744000;8555900];% cell numbers

%% select fixed parameters 
global a b theta_1 delta_E delta_R rho 
a = 0.2253;
b = 1/(4000/(1.2882*10^(-6)));
theta_1 = 11;
delta_E = 0.0330;
delta_R = 0.1;
rho = 0.0016; 

global  beta_1  m  n theta_2 alpha_1 beta_2 alpha_2 gamma 
global    beta_3 alpha_3 delta_X R_0 X_0																	
beta_1 =0.12577;
n = 6.0612E-07;
m = 2.5E-11; % which is smaller than m = 2.8743E-11;   
theta_2 =1.5319E-07;
alpha_1 =2.3122E+07;
beta_2 =0.11031;
alpha_2 =2.4036E+07;
gamma =6.3815E+04;
beta_3 =0.47822;
alpha_3 =1.4425E+07;
delta_X =0.070368;
R_0 =0.0100;
X_0 =0.5294;

%% long run

%%%------------------ mice 3--------------------%%%

figure;
hold on
for i = 3 % the 3rd mice
   
   
    opts = odeset('AbsTol',1e-7);

    tumplottime_org = fitData{i}.xdata(:,1);
    tcelplottime_org = fitData{i}.ydata(:,1);
   % long run
    tumplottime = [tumplottime_org(1):tumplottime_org(end)+17];
    tcelplottime = [tcelplottime_org(1):tcelplottime_org(end)+17];

    C_0 = 4.5622e+07;
    T_01 = 4.4E+5;
    T_02 =  4.7E+5;

    EstmCTint1 = [C_0; T_01*(1 - R_0- X_0);T_01*R_0; T_01*X_0];
    [t,modFit1] = ode45(@(t,y)fCE_final(t,y), tumplottime, EstmCTint1,opts);

    EstmCTint2 = [C_0; T_02*(1 - R_0- X_0);T_02*R_0; T_02*X_0];
    [t,modFit2] = ode45(@(t,y)fCE_final(t,y), tumplottime, EstmCTint2,opts);

   
    hold on
    plot(tumplottime,modFit1(:,1),'-r','LineWidth',3)
    plot(tumplottime,modFit2(:,1),'-b','LineWidth',3)
    hold off
    
    
    xlabel('Time [Days]','FontSize',18)
    ylabel('Tumor Cell Count C(t)','FontSize',18)
    legend({'T(0) = 4.4 \cdot 10^5','T(0) = 4.7 \cdot 10^5'},'Location','best')
    set(gca,'YScale', 'log')
    set(gca,'FontSize',18)


end

sgt = sgtitle('When m = 2.5000 \cdot 10^{-11} in mouse 3','Color','black');
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
