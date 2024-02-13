%DRAW DISTRIBUTIONS

clear all;
clc;

pkg load statistics

disp('')
list_choix = [1,2,3] %1 for variation on n, 2 for variation on rho and 3 for variation on p and L
disp('')


%Data to save

z_1=[0.001:0.001:0.999];
z_2 = [];
z_tot = [z_1]

%File to save

delete 'distributions.csv'
fichier = 'distributions.csv'

col_names = [];
col_names(1) = 1;
for j = 1 : 3
for k = 1 : 5
nb = 1 + k + 5*(j-1);
col_names(nb) = 1 + k + 5*(j-1);
end
end
dlmwrite(fichier,col_names,'-append')




for j = 1 : 3
  

choix = list_choix(j)


%Parameters

w = 50*10^3; %wealth of the insured
w_i = 100*10^3; %wealth of the insurer
eta = 2*10^(-5); %risk aversion of insured
eta_i = 1*10^(-5); %risk aversion of insurer

c = @(x) 3.*x;
rho = 0.08; %correlation across insureds
p = 0.004; %risk proability of insureds
L = 10*10^3; %risk loss of insureds
pL_value = 100;

n_list = [0.5,1.5,3,6,10];
rho_list = [0.01,0.04,0.08,0.15,0.25];
p_list = [0.001,0.002,0.004,0.008,0.015];


for k = 1 : 5

if (choix == 1)
n = n_list(k);
endif
if (choix == 2)
rho = rho_list(k);
endif
if (choix == 3)
p = p_list(k);
endif

%Distribution function for fraction of insureds affected

a = p*(1/rho-1)
b = (1-p)*(1/rho-1)

f =@(x)  betapdf(x,a,b);
F =@(x)  betacdf(x,a,b);


%Compute distribution function

nb = 1 + k + 5*(j-1);
z_2=f(z_1); 
z_tot = [z_tot; z_2]


##if (choix == 1)
##figure(01)
##semilogx(z_1*100,z_2, 'LineWidth', 2,'Color',[0 0 0]+0.4);
##endif 
##if (choix == 2)
##figure(02)
##semilogx(z_1*100,z_2, 'LineWidth', 2,'Color',[0 0 0]+(k-1)./5);
##hold on
##endif 
##if (choix == 3)
##figure(03)
##semilogx(z_1*100,z_2, 'LineWidth', 2,'Color',[0 0 0]+(k-1)./5);
##hold on
##endif 

end


end




for i = 1 : 999
  
z = transpose(z_tot(:,i));

dlmwrite(fichier,z,'-append','precision',5);

end











