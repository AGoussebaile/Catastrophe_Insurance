%FINDS THE OPTIMAL CONTRACT WITH DEFAULT
close all;
clear all;
clc;

pkg load statistics

global w w_i eta eta_i n p L




disp('')
choix = 3 %1 for variation on c, 2 for variation on rho and 3 for variation on p and L
disp('')




%Parameters

w = 50*10^3; %wealth of the insured
w_i = 100*10^3; %wealth of the insurer
eta = 2*10^(-5); %risk aversion of insured
eta_i = 1*10^(-5); %risk aversion of insurer

%n = 4; %number of agents over number of investors
c = @(x) 3.*x;
rho = 0.08; %correlation across insureds
p = 0.004; %risk proability of insureds
L = 10*10^3; %risk loss of insureds
pL_value = 100;

c_list = [0.5,1.5,3,6,10];
rho_list = [0.01,0.04,0.08,0.15,0.25];
p_list = [0.001,0.002,0.004,0.008,0.015];



epsilon_list = [];
shortfall_list = [];
Lambda_list = [];
Irat_list = [];
u_opt_list = [];
Wrat_list = [];







%Begin the 5-rounds analysis


for k = 1 : 5



disp('')
disp('')
disp('START')
disp('')

if (choix == 1)
c = @(x) c_list(k).*x
endif
if (choix == 2)
rho = rho_list(k)
endif
if (choix == 3)
p = p_list(k)
endif

L = pL_value/p

disp('')



%Utility functions

function u = u(z)
global eta;
u = -1./eta .* (exp(-eta.*z)) ;
end;

function u_prime = u_prime(z)
global eta;
u_prime = exp(-eta.*z);
end;

function v = v(z)
global w_i eta_i n;
v =z ;
end;

function v_prime = v_prime(z)
global w_i eta_i n;
v_prime = 1 ;
end;



%Distribution function for fraction of insureds affected

a = p*(1/rho-1);
b = (1-p)*(1/rho-1);

f =@(x)  betapdf(x,a,b);
F =@(x)  betacdf(x,a,b);
F_inverse =@(x) betainv(x,a,b)

p_bis = quadgk(@(x) x.*f(x),0,1);
rho_bis = quadgk(@(x) x.^2.*f(x)./(p.*(1-p)),0,1) + quadgk(@(x) -2.*x.*p.*f(x)./(p.*(1-p)),0,1) + p^2/(p*(1-p));

x_moy =  quadgk(@(x) x.*f(x),0,1);
sigma_2 = a*b/((a+b)^2*(a+b+1));
sigmaL_2 = sigma_2*L^2;
pL_val = p*L;



%Solve the problem for different policies (ie, different defaults or shortfalls)

nb_round = 47
epsilon_list(:,k) = [0.00001,0.00005,0.0001,0.0002,0.0004,0.0007,0.0010,0.0015,0.0020,0.0025,0.0030,0.0035,0.004,0.005,0.006,0.007,0.008,0.009,0.010,0.011,0.012,0.013,0.014,0.015,0.016,0.017,0.018,0.019,0.020,0.021,0.022,0.023,0.024,0.025,0.026,0.028,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.10,0.12,0.14,0.16];
shortfall_list(:,k) = [1:nb_round];
Lambda_list(:,k) = [1:nb_round];
Irat_list(:,k) = [1:nb_round];
u_opt_list(:,k) = [1:nb_round];
Wrat_list(:,k) = [1:nb_round];


for i = 1 : nb_round
  
disp('')
disp('ROUND BEGIN')
disp('')

disp('k')
disp(k)
disp('')

disp('ROUND')
disp(i)
disp('')

disp('OUT OF')
disp(nb_round)
disp('')

epsilon = epsilon_list(i);
F_inv = F_inverse(1-epsilon);

%FIND THE SUPPPLY FUNCTION P_STAR (insurer's participation constraint)
fun_a = @(I,P) P- I.*(quadgk(@(x) x.*f(x),0,F_inv) + F_inv*epsilon + quadgk(@(x) c(x).*(1-F(x)),P./I,F_inv)  ); %Expected utility of insurer
P_star0 = @(I) fsolve(@(P) fun_a(I,P),p*L); %Solves fun d to give P(I)
P_star = @(I) arrayfun(P_star0,I);

%FIND THE OPTIMAL INSURANCE LEVEL I (Solves the agent program, taking into account the participation constraint of the insurer)
fun_g = @(I) quadgk(@(x) u(w-L+I-P_star(I)).*x.*f(x),0,F_inv);
fun_i = @(I) quadgk(@(x) u(w-L+I.*F_inv./x-P_star(I)).*x.*f(x),F_inv,1);
EXP = @(I) (1-p).*u(w-P_star(I)) + fun_g(I) + fun_i(I);
EXP1 = @(I) arrayfun(EXP,I);
OBJ = @(I) - EXP1(I);

%SOLUTION
disp('SOLUTION OF THE CONTRACT')
tic
I_guess = L;
options = optimset('TolFun',1e-03);
I_opt = fminsearch(OBJ,I_guess,options);
toc
I_opt = min(I_opt,L);
P_opt = P_star(I_opt);
pL_default = quadgk(@(x) L.*x.*f(x),0,F_inv) + quadgk(@(x) L.*F_inv.*f(x),F_inv,1);


%Welfare with no insurance
u_no_ins = (1-p).*u(w) + p.*u(w-L);

%Welfare with free insurance
u_free = u(w - p*L);

%Welfare with optimal default
u_opt = (1-p).*u(w-P_opt) + quadgk(@(x) u(w-L+I_opt-P_opt).*x.*f(x),0,F_inv) + quadgk(@(x) u(w-L+I_opt.*F_inv./x-P_opt).*x.*f(x),F_inv,1); %Welfare at optimal default contract allocation

%Welfare gain
u_gain = (u_opt - u_no_ins)./(u_free - u_no_ins);

%Results
epsilon_list(i,k) = epsilon;
shortfall_list(i,k) = I_opt.*quadgk(@(x) F_inverse(x),1-epsilon,1) - I_opt.*quadgk(@(x) F_inv*f(x),F_inv,1);
Lambda_list(i,k) = (P_star(L)-pL_default)/pL_default;
Irat_list(i,k) = I_opt/L;
u_opt_list(i,k) = u_opt;
Wrat_list(i,k) = u_gain;
%end

disp('')
disp('RESULTS')
printf('I_opt : %4.2f, P_opt : %4.2f, pL_default : %4.2f', I_opt, P_opt, pL_default);
disp('')
printf('Load-fullcov : %4.2f, Indem-rate : %4.2f, u_opt : %4.2f', Lambda_list(i,k), Irat_list(i,k), u_opt_list(i,k));
disp('')

disp('')
disp('ROUND END')
disp('')

end



disp('')
disp('')
disp('RESULTS SUMMARY')
disp('')
disp('')
printf('w : %4.2f, w_i : %4.2f, eta : %4.2E, eta_i : %4.2E', w, w_i, eta, eta_i);
disp('')
printf('n : %4.2f, rho_bis : %4.2f, p_bis : %4.2E , L : %4.2f, pL : %4.2f', n, rho_bis, p_bis, L, pL_val);
disp('')
printf('a : %4.2f, b : %4.2f, x_moy : %4.2E, sigmaL_2 : %4.2E', a, b, x_moy, sigmaL_2);
disp('')
disp('')



end






%Save the data in csv file


if (choix == 1)
delete '1nVar.csv'
fichier = '1nVar.csv'; %Pour variation de n
endif
if (choix == 2)
delete '2corVar.csv'
fichier = '2corVar.csv'; %Pour variation de rho
endif
if (choix == 3)
delete '3proVar.csv'
fichier = '3proVar.csv'; %Pour variation de p et L
endif


col_names = [];
for k = 1 : 5
nb = 1+5*(k-1)
col_names(nb) = 1+5*(k-1);
nb = 2+5*(k-1)
col_names(nb) = 2+5*(k-1);
nb = 3+5*(k-1)
col_names(nb) = 3+5*(k-1);
nb = 4+5*(k-1)
col_names(nb) = 4+5*(k-1);
nb = 5+5*(k-1)
col_names(nb) = 5+5*(k-1);
end
dlmwrite(fichier,col_names,'-append')


for i = 1 : nb_round

%para_list = [w,w_i,eta,eta_i,n,rho_bis,p_bis,L,pL_val,a,b,x_moy,sigmaL_2];
%dlmwrite(fichier,para_list,'-append','precision',6)

list = [];
for k = 1 : 5
nb = 1+5*(k-1)
list(nb) = epsilon_list(i,k);
nb = 2+5*(k-1)
list(nb) = shortfall_list(i,k);
nb = 3+5*(k-1)
list(nb) = Lambda_list(i,k);
nb = 4+5*(k-1)
list(nb) = Irat_list(i,k);
nb = 5+5*(k-1)
list(nb) = Wrat_list(i,k);
end
dlmwrite(fichier,list,'-append','precision',5)


end





##%Plot some graphs
##
##
##for k = 1 : 5
##
##
##%Graphs with respect to default
##
##if (choix == 1)
##figure(11)
##endif
##if (choix == 2)
##figure(12)
##endif
##if (choix == 3)
##figure(13)
##endif
##plot(epsilon_list(:,k)*100,Lambda_list(:,k)*100, 'LineWidth', 2,'Color',[0 0 0]+(k-1)./5);
##hold on
##xlabel ("Default rate (%)");
##ylabel ("Loading factor (%)");
##axis([0 2.5 0 60]);
##if (choix == 1)
##figure(14)
##endif
##if (choix == 2)
##figure(15)
##endif
##if (choix == 3)
##figure(16)
##endif
##plot(epsilon_list(:,k)*100,Irat_list(:,k)*100, 'LineWidth', 2,'Color',[0 0 0]+(k-1)./5);
##hold on
##xlabel ("Default rate (%)");
##ylabel ("Coverage rate (%)");
##axis([0 2.5 0 100]);
##if (choix == 1)
##figure(17)
##endif
##if (choix == 2)
##figure(18)
##endif
##if (choix == 3)
##figure(19)
##endif
##plot(epsilon_list(:,k)*100,Wrat_list(:,k)*100, 'LineWidth', 2,'Color',[0 0 0]+(k-1)./5);
##hold on
##axis([0 2.5 0 100]);
##xlabel("Default rate (%)");
##ylabel("Relative welfare gain (%)");
##
##
##%Graphs with respect to shortfall
##
##if (choix == 1)
##figure(21)
##endif
##if (choix == 2)
##figure(22)
##endif
##if (choix == 3)
##figure(23)
##endif
##plot(shortfall_list(:,k),Lambda_list(:,k)*100, 'LineWidth', 2,'Color',[0 0 0]+(k-1)./5);
##hold on
##xlabel ("Shortfall");
##ylabel ("Loading factor (%)");
##axis([0 30 0 60]);
##if (choix == 1)
##figure(24)
##endif
##if (choix == 2)
##figure(25)
##endif
##if (choix == 3)
##figure(26)
##endif
##plot(shortfall_list(:,k),Irat_list(:,k)*100, 'LineWidth', 2,'Color',[0 0 0]+(k-1)./5);
##hold on
##xlabel ("Shortfall");
##ylabel ("Coverage rate (%)");
##axis([0 30 0 100]);
##if (choix == 1)
##figure(27)
##endif
##if (choix == 2)
##figure(28)
##endif
##if (choix == 3)
##figure(29)
##endif
##plot(shortfall_list(:,k),Wrat_list(:,k)*100, 'LineWidth', 2,'Color',[0 0 0]+(k-1)./5);
##hold on
##axis([0 30 0 100]);
##xlabel("Shortfall");
##ylabel("Relative welfare gain (%)");
##
##
##end