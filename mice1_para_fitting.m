clear all 

addpath( './DRAM_Code/'); 
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
%% define parameters to estimate 
% { ParamName, starting value, uniform prior bounds }     
params = {   
    %%%%%%%%% first eq
    {'n',     3*0.1^7,      3*0.1^9,     10*0.1^7}
    %%%%%%%%% second eq
    {'beta_1',     0.12,          0.1,       0.3/2}  
    {'alpha_1',       2e7,    1e7,   10e7}
    {'theta_2', 2*0.1^7,   0.1^8,    1*0.1^4 }
    {'m',     3*0.1^10,   3*0.1^11,  3*0.1^7}
    %%%%%%%%% third eq
    {'beta_2',    0.12,          0.1,        0.3/2}
    {'alpha_2',    2e7,       1e7,         10e7}
    {'gamma',    4e4,         1e3,        1e5}
    %%%%%%%%%%% fourth eq
    {'beta_3',     0.2,           0.1,       0.3}
    {'alpha_3', 2e7,         1e7,      1e7*10}
    {'delta_X',    0.0330 ,  0.0330,  0.15}
     %%%%%%%%%%%% initial procentage
    {'R_0', 0.05,         0.01,      0.1}
    {'X_0', 0.6,         0.5,      1}
        };
    
%% This is the error function 
model.ssfun = @LLHfunc; 

options.nsimu = 1*10^4; % number of samples for DRAM

%% Run mcmc 
[results, chain, s2chain, ss2chain] = mcmcrun(model,fitData,params,options);

%% Plot results 

%%% plot to check chain convergence 
% figure; mcmcplot(chain,[],results,'chainpanel');

%%% plot to check pair of chain relation  
% %%The 'pairs' options makes pairwise scatterplots of the columns of the chain.
% figure; mcmcplot(chain,[],results,'pairs');

%%%% plot this to check pair of chain distribution 
figure; mcmcplot(chain,[],results,'denspanel');

%% Compare data to calibrated model
ind = find(ss2chain == min(ss2chain)); % find the index where we get the min error
ind = ind(1);
params = chain(ind,:); %These are your fitted parameter values

if( norm( params - mean(chain) )/norm( mean(chain) ) > 0.1 ); 
    disp( 'Fitted parameter is off from the mean of chain - check parameter distribution' ); 
end 


%%%------------------ mice 1--------------------%%%
figure;
hold on
for i = 1 
   
    opts = odeset('AbsTol',1e-5);

    tumplottime = fitData{i}.xdata(:,1);
    tcelplottime = fitData{i}.ydata(:,1);

    CTint= [fitData{i}.xdata(1,2);fitData{i}.ydata(1,2)];
    % Note that we need to estimate effector T cell (E), Tregs (R), and exhausted T
    % cell (X). That is, Tcell = E + R + X
    EstmCTint = [CTint(1); CTint(2)*(1 - params(12)- params(13));...
                CTint(2)*params(12); CTint(2)*params(13)];


    [t,modFit] = ode45(@(t,y)fCE(t,y,params), tumplottime, EstmCTint,opts);

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
    legend({'Total T','E', 'R','X','data T'},'Location','best')
    xlabel('Time in Days','FontSize',18)
    ylabel('Cell Count','FontSize',18)
    ylim([5e3 7e7])
    set(gca,'YScale', 'log')
    set(gca,'FontSize',18)
    hold off
   

end

sgt = sgtitle('Mouse 1','Color','black');
sgt.FontSize = 20;

hold off


%% Functions of errors

%%% Choose 'Sum of squares' 

function SS = LLHfunc(params,fitData)
SS = 0;
for i = 1 
    
    opts = odeset('AbsTol',1e-1);
   
    CTint= [fitData{i}.xdata(1,2);fitData{i}.ydata(1,2)];
    % Note that we need to estimate effector T cell (E), Tregs (R), and exhausted T
    % cell (X). That is, Tcell = E + R + X
    estmCTint = [CTint(1);CTint(2)*(1 - params(12)-params(13));...
                 CTint(2)*params(12); CTint(2)*params(13)];
    
    tspan= 0:1:53;
    [t,tumTcel] = ode45(@(t,y)fCE(t,y,params), tspan, estmCTint,opts);
    %tumTcel: the first column is tumor cells; the second, third and fourth 
    % ones are the T cells 

    % Remove the other time points in estimated tumor number
    Lia = ~ismember(tspan,fitData{i}.xdata(:,1));
    tumEvl = tumTcel(:,1);
    tumEvlfixtim = tumEvl(~Lia);
    
    % Remove the other time points in estimated T cell
    Lia = ~ismember(tspan,fitData{i}.ydata(:,1));
    % effector CD8+ T cell at our fixed time points
    efTcellEvl = tumTcel(:,2);
    efTcellEvlfixtim = efTcellEvl(~Lia);
    % Treg cell at our fixed time points
    TegcellEvl = tumTcel(:,3);
    TregcellEvlfixtim = TegcellEvl(~Lia);
    % exhausted CD8+ T cell at our fixed time points
    excellEvl = tumTcel(:,4);
    excellEvlfixtim = excellEvl(~Lia);
    % SUM 
    TcellEvlfixtim = efTcellEvlfixtim + TregcellEvlfixtim + excellEvlfixtim;
    

    %%% tumor vulome error
    relerr =  fitData{i}.xdata(:,2); 
    if( size( tumEvlfixtim,1 ) == size( fitData{i}.xdata(:,2),1 ) )
        SS = SS + sum(((tumEvlfixtim(:,1) - fitData{i}.xdata(:,2)).^2)./(relerr.^2));
    else
        SS = SS + 10^10;
    end

    %%% T cell error
    relerr =  fitData{i}.ydata(:,2); 
    if( size( TcellEvlfixtim,1 ) == size( fitData{i}.ydata(:,2),1 ) )
        SS = SS + sum(((TcellEvlfixtim(:,1) - fitData{i}.ydata(:,2)).^2)./(relerr.^2));
    else
        SS = SS + 10^10;
    end
    
   
    
end

end



function value = fCE( t, y, param)  

 global a b theta_1 rho delta_E delta_R 


C = y(1); % cancer cell
E = y(2); % effector CD8+ T cell
R = y(3); % treg  T cell
X = y(4); % exhausted CD8+ T cell


%%%%%%%%%%% first eq: cancer eq
n = param(1);
%%%%%%%%%%% second eq: effector CD8+ T cell
beta_1 = param(2);
alpha_1 = param(3);
theta_2 = param(4);
m = param(5);
%%%%%%%%%%% third eq: Treg cell
beta_2 = param(6);
alpha_2 = param(7);
gamma = param(8);
%%%%%%%%%%% fourth eq: exhausted CD8+ T cell
beta_3 = param(9);
alpha_3 = param(10);
delta_X = param(11);
%%%%%%%%%%%%


value(1,1) =  a*C*(1-b*C) - n*C*E*( 1/(1+theta_1*R/E) ); 

value(2,1) = -delta_E*E + beta_1*E*(C/(alpha_1+C))*(1/(1+theta_2*R)) ...
             - m*E*C - rho*E;
          
value(3,1) = -delta_R*R + beta_2*R*(C/(alpha_2+C))*(E/(gamma +E)) + rho*E;

value(4,1) = -delta_X*X + beta_3*X*(C/(alpha_3+C))*(1/(1+theta_2*R)) + m*E*C;


end    
