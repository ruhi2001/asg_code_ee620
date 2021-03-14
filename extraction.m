%for nmos
clear; close all
d1 = csvread('C_valS.csv');
d2 = csvread('V_valS.csv');
Vs = d2(1,1:end-1);
Cs = d1(1,2:end);
%%
q = 1.6e-19;             %charge
eps_0 = 8.85e-12; 
Vt = 26e-3;
k_si = 12; 
eps_si = k_si*eps_0;
k_ox = 4;
eps_ox = k_ox*eps_0;
ni = 1.5e10*1e6;
Eg = 1.1*q;

%%
figure(1);
plot(Vs,Cs);

cox = Cs(1);
cmin = Cs(end);
tox = eps_ox/cox;
cdep = 1/((1/cmin) - (1/cox));
wdep = eps_si/cdep;

Na = 1.5e17*1e6; %initial guess
for i = 1:10
    phi_b = Vt*log(Na/ni);
    Na = 4*eps_si*phi_b/(q*wdep^2);
end

phi_b = Vt*log(Na/ni);
Vfb = -Eg/2/q - phi_b;

wmg = (2*eps_si*phi_b/q/Na)^0.5;
cmg = 1/((1/cox) + wmg/eps_si);

vmg_dep = Vfb + (2*eps_si*phi_b*q*Na)^0.5/cox + phi_b;
%vmg_m = interp(Cs,Vs,cmg);
%dvmg = vmg_m - vmg_dep;
%Qox = -dvmg*cox;



