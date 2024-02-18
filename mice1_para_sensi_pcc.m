function test_prcc() 

N = 10000; % # or samples 
rDim = 19; % # of parameters 
LHSample= lhsdesign(N,rDim) ; 

paramFit = set_param(); 

errorsd = 0.1; 
ub = paramFit + paramFit*errorsd;
lb = paramFit - paramFit*errorsd;

for n = 1:rDim 
    LHSample(:,n) = LHSample(:,n) *(ub(n)-lb(n)) + lb(n); 
end

for n = 1:N 
   yy(n,:) = run_fCERX( LHSample(n,:) ); 
end 


[rho,pval] = partialcorr([LHSample,yy(:,1)], 'Type', 'Pearson');  % linear cc 
figure; bar( rho(1:end-1,end) )
set(gca,'XTick', 1:rDim, 'XTicklabel',{ 'a' 'b' 'n' '\theta_1'...
         '\delta_E' '\beta_1' '\alpha_1' '\theta_2' 'm' '\lambda'...
         '\delta_R' '\beta_2' '\alpha_2' '\gamma' ...
         '\delta_X' '\beta_3' '\alpha_3' 'R(0)' 'X(0)'}); 
ylabel('PCC'); 
title('Mouse 1: Parameter Change by 10 Percent');
set(gca,'FontSize',14)
    

end 

function param = set_param() 

%%% mice 1
a = 0.2253;
b = 1/(4000/(1.2882*10^(-6)));

beta_1 = 0.10077;
m = 6.1562E-10;
n =4.1965E-07;
theta_2 =1.4985E-07;
alpha_1 = 2.4625E+07
beta_R =0.10249;
alpha_R = 1.0132E+07;
gamma = 1.1002E+03
beta_X =0.11361*2.1
alpha_X = 2.5556E+07;
delta_X =0.085333

theta_1 = 11;
delta_E = 0.0330;
delta_R = 0.1;
rho = 0.0016; 

R_0 =0.017916;
X_0 = 0.50019;

param = [a b n theta_1...
         delta_E beta_1 alpha_1 theta_2 m rho...
         delta_R beta_R alpha_R gamma...
         delta_X beta_X alpha_X...
         R_0 X_0]

end 



function result = run_fCERX( param ) 

CTint(1) = 347.98/(1.2882*10^(-6)); % initial cancer cells
CTint(2) = 3259100; % initial total T cells
 
EstmCTint = [CTint(1); CTint(2)*(1 - param(18)- param(19));...
                CTint(2)*param(18); CTint(2)*param(19)];

opts = odeset('AbsTol',1e-5);
[time, result] = ode45(@(t,y)fCERX(t,y,param), [0:50], EstmCTint,opts);
result = result(end,:); 

end 



function value = fCERX(t, y, param)  

   
C = y(1); % cancer cell
E = y(2); % effector CD8+ T cell
R = y(3); % treg  T cell
X = y(4); % exhausted CD8+ T cell


value(1,1) =  param(1)*C*(1-param(2)*C) - param(3)*C*E*( 1/(1+param(4)*R/E) ); 

value(2,1) = -param(5)*E + param(6)*E*(C/(param(7)+C))*(1/(1+param(8)*R)) - param(9)*E*C - param(10)*E;
          
value(3,1) = -param(11)*R + param(12)*R*(C/(param(13)+C))*(E/(param(14)+E)) + param(10)*E;

value(4,1) = -param(15)*X + param(16)*X*(C/(param(17)+C))*(1/(1+param(8)*R)) + param(9)*E*C;


end    


