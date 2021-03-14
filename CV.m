%%Parameters

%constants
q = 1.6e-19;             %charge
eps_0 = 8.85e-12; 
VkBT = 26e-3;            %KBT/q at 300k

%semiconductor
Nsub = -1e18*1e6;        %negative for NMOS and positive for PMOS
k_si = 12;               %Dielectric Constant
ni = 1.5e10*1e6;         %Intriensic carrier concentration
Eg = 1.1*q;              %Bandgap
L = 250e-9;                %substrate thickness
eps_si = k_si*eps_0;     %Permitivity
chi_si = 4.05*q;

if Nsub < 0
    psub = abs(Nsub); nsub = ni^2/abs(Nsub);  
elseif Nsub > 0
    psub = ni^2/Nsub; nsub = Nsub;
else
    psub = ni; nsub = ni;
end

%Oxide
tox = 3e-9;            %Thickness
k_ox = 4;               %Dielectric constant
Eg_ox = 9*q;            %Bandgap
eps_ox = k_ox*eps_0;    %Permittivity
Ec_off = 3.1*q;         %Conduction band offset
C_ox = eps_ox/tox;

%Flatband/Metal related
tm = 3e-9;   %Thickness
phi_m = chi_si/q;
phi_b = -sign(Nsub)*VkBT*log(abs(Nsub)/ni); %(Ei - Ef)/q of semiconductor
phi_s = chi_si/q -sign(Nsub)*Eg/(2*q) + phi_b;
Vfb = phi_m - phi_s;            %n+/nmos, p+/pmos

%solver
dx = 1e-9;    %meshing size
NRmax = 10;   %Maximum NR iterations
tol = 1e-6;    %Error tolerance in potential

%% mesh
xox = -tox:dx:0;
xsi = dx:dx:L;
x = [xox,xsi];
N = length(x);    %number of points
iox = (x <= 0);
isi = (x > 0);
i0 = find(x==0);
 
%% Newton Raphson - Matrices

ki = k_ox*iox + k_si*isi;
kip = 0.5*(ki + circshift(ki,-1));
kim = 0.5*(ki + circshift(ki,1));
kipm = kip + kim;
A = (diag(kipm,0) - diag(kip(1:end-1),1) - diag(kim(2:end),-1));
A(1,:) = 0; A(1,1) = 1;
A(N,:) = 0; A(N,N) = 1;

b = zeros(N,1);
psi = zeros(N,1);

%% Input

 Vs = (-3:0.01:3);%different Vg's
 Vwig = 10e-3;
 Vright = 0;
 %sweep_type = 1;0 = LFCV,1 = HFCV,2 = Fast sweep
 
 
%% Newton-Raphson - Iterations
for sweep_type = 1
    Caps = zeros(length(Vs),1)';
    psi_ss = zeros(length(Vs),1)';
    
    %voltage sweep
    for j=1:length(Vs)
        Vg = Vs(j);
        %inputs at each voltage
        V1 = Vg -Vfb;
        V2 = V1 + Vwig;
        Qsub = [];
        
        %bias point
        Vleft = V1;
        for i = 1:NRmax
            %figure(2);plot(x,psi);hold on
            %generate b
            psi_si = psi(isi);
            psi_s = psi(i0);
            p = psub*exp(-psi_si/VkBT);
            n = nsub*exp(psi_si/VkBT);
            delp = -1/VkBT*p;
            deln = 1/VkBT*n;
            
            if sweep_type == 2 
                minc = 0;
                delminc = 0;
    
                if sign(Nsub) == 1
                    majc = n;
                    delmajc = deln;
                elseif sign(Nsub) == -1
                    majc = p;
                    delmajc = delp;
                end
                
            else 
                if sign(Nsub) == 1
                    majc = n;
                    delmajc = deln;
                    minc = p;
                    delminc = delp;
                elseif sign(Nsub) == -1
                    majc = p;
                    delmajc = delp;
                    minc = n;
                    delminc = deln;
                end
            end
           rho_ox = zeros(size(xox)); 
           rho_si = q*(Nsub + sign(Nsub)*(minc - majc));
           rho = [rho_ox'; rho_si];
           b = rho*dx^2/eps_0;
           b(1) = Vleft; b(N) = Vright;
           f = A*psi - b;
           
           %jacobian
           delrho_ox = zeros(size(xox));
           delrho_si = q*(sign(Nsub)*(delminc - delmajc));
           
           delrho = [delrho_ox'; delrho_si];
           delb = delrho*dx^2/eps_0;
           delb(1) = 0; delb(N) = 0;
           J = A - diag(delb);
           % J(i0,i0 + 1) = J(i0,i0 + 1) - ditdq*dx^2/eps_0;
           dV = -J\f;
           if max(abs(dV)) < tol
               break;
               
           end
           psi = psi + dV;
        end
        
        Q1 = sum(rho_si)*dx;
        minc1 = minc;
        
        %wiggle point
        Vleft = V2;
        for i = 1:NRmax
            %figure(2);plot(x,psi);hold on
            %generate b
            psi_si = psi(isi);
            psi_s = psi(i0);
            
            p = psub*exp(-psi_si/VkBT);
            n = nsub*exp(psi_si/VkBT);
            delp = -1/VkBT*p;
            deln = 1/VkBT*n;
            
            if sweep_type == 2 
                minc = 0;
                delminc = 0;
                if sign(Nsub) == 1
                    majc = n;
                    delmajc = deln;
                elseif sign(Nsub) == -1
                    majc = p;
                    delmajc = delp;
                end
            elseif sweep_type == 1
                minc = minc1;
                delminc = 0;
                if sign(Nsub) == 1
                    majc = n;
                    delmajc = deln;
                elseif sign(Nsub) == -1
                    majc = p;
                    delmajc = delp;
                end
             elseif sweep_type == 0 
                if sign(Nsub) == 1
                    majc = n;
                    delmajc = deln;
                    minc = p;
                    delminc = delp;
                elseif sign(Nsub) == -1
                    majc = p;
                    delmajc = delp;
                    minc = n;
                    delminc = deln;
                end
            end
           rho_ox = zeros(size(xox)); 
           rho_si = q*(Nsub + sign(Nsub)*(minc - majc));
           rho = [rho_ox'; rho_si];
           b = rho*dx^2/eps_0;
           b(1) = Vleft; b(N) = Vright;
           f = A*psi - b;
           
           %jacobian
           delrho_ox = zeros(size(xox));
           delrho_si = q*(sign(Nsub)*(delminc - delmajc));
           
           delrho = [delrho_ox'; delrho_si];
           delb = delrho*dx^2/eps_0;
           delb(1) = 0; delb(N) = 0;
           J = A - diag(delb);
           % J(i0,i0 + 1) = J(i0,i0 + 1) - ditdq*dx^2/eps_0;
           dV = -J\f;
           if max(abs(dV)) < tol
               break;
               
           end
           psi = psi + dV;
        end
        
        Q2 = sum(rho_si)*dx ;
        psi_ss(j) = psi_s;
        Cap = 1/Vwig*(Q1 - Q2);
        Caps(j) = Cap;
    end  
    figure(3);
    plot(Vs(Vs>-2.96),Caps(Vs>-2.96));
    
    xlabel('Vg (V)');
    ylabel('Capacitance (F/m2)');
    grid on
    hold on
   
end
csvwrite('C_valS.csv',Caps)
csvwrite('V_valS.csv',Vs)

