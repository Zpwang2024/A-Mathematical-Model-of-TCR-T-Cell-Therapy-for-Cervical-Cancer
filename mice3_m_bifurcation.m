%%Initialize the environment

close all;
clc;
clear model data params options
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

% orignial estimated parameter values
 ao = 0.2253;
 bo = 4000/(1.2882*10^(-6));
 theta_1o = 11;
 delta_Eo = 0.0330;
 delta_Ro = 0.1;
 rhoo = 0.0016; 
 beta_1o = 0.12577;
% mo = 5.2261E-10; % here bifurcation for m value
 no = 6.0612E-07;
 theta_2o = 1.5319E-07;
 alpha_1o = 2.3122E+07;
 beta_2o =  0.11031;
 alpha_2o = 2.4036E+07;
 gammao = 6.3815E+04;

global a b
global delta_E delta_R theta_1 rho 
global  beta_1  m  n theta_2 alpha_1 beta_2 alpha_2 gamma 

% order-ofmagnitude concentration scale for 
% the C, E and R cell populations
Co = 1e6;
Eo = 1e6;
Ro = 1e6;
% scaled parameters
a = ao/(no*Co);
b = bo/Co;
beta_1 = beta_1o/(no*Co);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
m = 3e-4; % actually, here m = mo/no
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n = 1;
theta_2 = theta_2o*Ro;
alpha_1 = alpha_1o/Co;
beta_2 = beta_2o/(no*Co);
alpha_2 = alpha_2o/Co;
gamma = gammao/Eo;
theta_1 = theta_1o;
delta_E = delta_Eo/(no*Co);
delta_R = delta_Ro/(no*Co);
rho = rhoo/(no*Co);

%% Calculate Cstar, Estar and Rstar and Analysis the systerm   

%%%% the trivial point (0,0,0)
y00 = [1e-8,1e-8,1e-8]
options = optimoptions('fsolve','Display','iter',...
       'MaxFunctionEvaluations',1e5,'FunctionTolerance',1e-22,...
       'MaxIterations',1e5,'OptimalityTolerance',1e-22,...
       'StepTolerance',1e-22);
[y0, fval0] = fsolve(@fCER,y00,options);
disp('The critical point y0 :'); 
disp(vpa([y0(1),y0(2),y0(3)],5))
% check y0 satisify the conditions
disp('check y0 satisify the conditions:'); 
disp(vpa([fval0(1),fval0(2),fval0(3)],5))

%%%% the capacity point (b,0,0) 
y01 = [b,1e-8,1e-8];
[y1, fval1] = fsolve(@fCER,y01,options);
disp('The critical point y1 :'); 
disp(vpa([y1(1),y1(2),y1(3)],5))
% check y1 satisify the conditions
disp('check y1 satisify the conditions:'); 
disp(vpa([fval1(1),fval1(2),fval1(3)],5))

%%%% the critical piont (Cstar, Estar, Rstar)
y02 = [8,0.4,0.01];
[y2, fval2] = fsolve(@fCER,y02,options);
disp('The critical point y2 :'); disp(vpa([y2(1),y2(2),y2(3)],4))
disp('check y2 satisify the conditions:'); 
disp(vpa([fval2(1),fval2(2),fval2(3)],5))

%%%% the critical piont (Cstar, Estar, Rstar)
y03 = [250,3,3];
[y3, fval3] = fsolve(@fCER,y03,options);
disp('The critical point y3 :'); 
disp(vpa([y3(1),y3(2),y3(3)],4))
disp('check y3 satisify the conditions:'); 
disp(vpa([fval3(1),fval3(2),fval3(3)],5))

% store critical pionts
yvalue = zeros(4,3);
yvalue(1,:) = y0;
yvalue(2,:) = y1;
yvalue(3,:) = y2;
yvalue(4,:) = y3;

%% Find the eigen-value and eigen-vector at the equibrium

%%%%% fixed parameters
global a b
global delta_E delta_R theta_1 rho 
global m beta_1 n alpha_1 theta_2 beta_2 gamma 


syms C E R 

sys1 =  a*C*(1-C/b) - n*C*E*( 1/(1+theta_1*R/E) ); 

sys2 = -delta_E*E + beta_1*E*(C/(alpha_1+C))*(1/(1+theta_2*R)) - m*E*C - rho*E;
          
sys3 = -delta_R*R + beta_2*R*(C/(alpha_2+C))*(E/(gamma+E)) + rho*E;

% computes the Jacobian matrix of symbolic function f 
% with respect to [C,E,R]
fCERjac = jacobian([sys1 sys2 sys3],[C,E,R]);

for i=2:1:4 % note the first one (0,0,0) is always unstable.

disp('The critical point is :'); 
disp( vpa(yvalue(i,:), 5) )

fCERjac_subNum = subs(fCERjac,[C,E,R],[yvalue(i,1),yvalue(i,2),yvalue(i,3)] );
disp('The linearization of the system at this point:');
disp(vpa(fCERjac_subNum, 5));

fCERlambda = eig(fCERjac_subNum);
disp('Eigenvalues associated with the linearized system at this point:'); 
disp(vpa(fCERlambda, 5));

end


% Lastly, restore the old value of digits for further calculations.
digits(digitsOld)

%%%----------------------------END------------------------------%%%

function fvalue = fCER(y)

global a b 
global n beta_1 m theta_2 alpha_1 beta_2 alpha_2 gamma 
global delta_E delta_R theta_1 rho
 

fvalue(1) = a*y(1)*(1-y(1)/b) - ...
            n*y(1)*y(2)*( 1/(1+ theta_1*( y(3)/y(2) ) ) );
fvalue(2) = -delta_E*y(2) + beta_1*y(2)*(y(1)/(alpha_1+y(1)))*...
            (1/(1+theta_2*y(3))) - m*y(2)*y(1) - rho*y(2);
fvalue(3) = -delta_R*y(3) + beta_2*y(3)*(y(1)/(alpha_2+y(1)))*...
            (y(2)/(gamma+y(2))) + rho*y(2);

end
