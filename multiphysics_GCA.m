function multiphysics_GCA
% This code is used to calculate the laser-induced crystallization on phase change materials
% It combines laser heating, GCA, effective medium theory and Fresnel equations
% Sample structure: 50nm SiN + 30nm GST
% More information can be found in this paper:https://arxiv.org/abs/2107.02035
% Welcome to citing this paper if it is useful for you.

%% %%%%% Initiallization %%%%%
disp('More information can be found in this paper:https://arxiv.org/abs/2107.02035.');
disp('Please cite this paper if you fully or partially use it in your paper.');
%% %%%Thermaldynamics parameters%%%
%GST
k_GST_a = 0.2;%thermal conductivity of a-GST, W/(mK)
k_GST_c = 0.5;%thermal conductivity of c-GST, W/(mK)
Cp_GST = 1.25e6;%density X specific heat, J/(m^3 K)
alpha_GST_a = k_GST_a/Cp_GST;
alpha_GST_c = k_GST_c/Cp_GST;
%SiN
k_SiN = 20;%thermal conductivity of SiN, W/(mK)
Cp_SiN = 700*3100;%density X specific heat, J/(m^3 K)
alpha_SiN = k_SiN/Cp_SiN;

%% %%%Geometry and time parameters%%%
Lx = 20e-6;%x length
Ly = 20e-6;%y length
Lz_GST = 30e-9;%GST thickness
Lz_SiN = 50e-9;%SiN thickness
%nonuniform meshing
h = 15e-9;%minimum space step
%x/y 500nm@0-6um 200nm@6.2-8um 100nm@8.1-8.5um 50nm@8.55-9.1um 30nm@9.13-10.9um 
     %50nm@10.95-11.5um 100nm@11.6-12um 200nm@12.2-14um 500nm@14.5-20um
%z   [0 15nm 30nm 50nm 80nm]
x = [0:0.5:6 6.2:0.2:8 8.1:0.1:8.5 8.55:0.05:9.1 9.13:0.03:10.9 10.95:0.05:11.5 11.6:0.1:12 12.2:0.2:14 14.5:0.5:20]*1e-6;
y = x;
z = [0 0.015 0.03 0.05 0.08]*1e-6;
%space interval
dx=diff(x);dy=diff(y);dz=diff(z);
dx=[dx dx(end)];dy=[dy dy(end)];dz=[dz dz(end)];
[dX,dY,dZ]=meshgrid(dx,dy,dz);
%time interval
FN = 1/6;%Stable condition: Fourier number<=1/6
dt = h^2*FN/max([alpha_SiN alpha_GST_a alpha_GST_c]);
dt = min(dt,4e-12);
%time rotation strategy
dt1 = 2*(1+sqrt(2)/2)*dt; dt2 = 2*(1-sqrt(2)/2)*dt;
dt_select = 1; %which step is selected

%% %%%Laser parameters%%%
w0 = 0.6e-6;%laser radius = 0.6 um;
Pin = 5e-3;%laser power
pulseWidth = 300e-9;%laser pulse width
edgeWidth = 8e-9;%laser edge width
[X,Y] = ndgrid(x,y);
I = 2*Pin/pi/w0^2*exp(-2*((X-Lx/2).^2+(Y-Ly/2).^2)/w0^2);%incident laser intensity
t_final = 2*pulseWidth; %end time

%% %%%optical dielectric parameters%%% 
wl = 660e-9;%center wavelength
%GST
N_GST_c = 4.61 - 4.01i;%complex refractive index of c-GST
N_GST_a = 4.36 - 1.79i;%complex refractive index of a-GST
epsilon_a = N_GST_a^2;%dielectric constant of c-GST
epsilon_c = N_GST_c^2;%dielectric constant of c-GST
%SiN
N_SiN = 2.0193 - 0i;%complex refractive index of SiN

%% %%%initial thermal conditions%%%
T0 = 298.15;%room temperature = 25¡æ
T = T0*ones(length(x),length(y),length(z));%temperature array
Temp = T;

%% %%%data to save%%%
t_save = 1e-10; %time interval to save data
Num_save = ceil(t_final/t_save)+1;%total length of saved data
t_heat = zeros(1,Num_save)*nan; t_heat(1) = 0; %heat calculation time
Tmax_heat = zeros(1,Num_save)*nan; Tmax_heat(1) = T0;%peak temperature array
Xf_save = zeros(1,Num_save)*nan;%crystalized fraction
t_Xf =zeros(1,Num_save)*nan;%Xf calculation time
Reflectivity = zeros(1,Num_save)*nan;%Reflectivity
Transmissivity = zeros(1,Num_save)*nan;%Transmissivity

%% %%%GCA parameters%%%
kB_eV = 8.617343e-5;%1st Boltzmann constant, eV/K 
kB_J = 1.38e-23;%2nd Boltzmann constant, J/K 
r_atom = 0.1365e-9;%atomic radius
Rh=0.1365e-9;%hydrodynamic radius
Ea = 2.45;%activation energy, eV 
sigma = 0.06;%interfacial energy,J/m^2 
vm = 2.9e-28;%volume of GST monomer, m^3 
d_GST = 2*(3/4/pi*vm)^(1/3);%dimeter of a GST molecule, m
lambdaj = 2.99e-10;%jump distance, m 
Tm=900;%melting temperature, K 
DHf=625e6;%entropy of fusion, J/m^3 
eta_inf = 0.012;%viscosity at infinite temperature
Tg = 472;%glass transition temperature, K
m = 140;%fragility
Tglass = 534;%transition temerature between supercooled state and glass state
theta = 130*pi/180;%wetting angle
f = (2-3*cos(theta)+cos(theta)^3)/4;%wetting factor
nij = 8;%number of neighbors in square lattice

%%%%calculation area of GCA%%%%
N = 400;
D = 6*d_GST;% space interval
h_GCA = D;
x_GCA_min = Lx/2-N*h_GCA/2;
x_GCA_max = Lx/2+N*h_GCA/2;
x_GCA = x_GCA_min:h_GCA:x_GCA_max-h_GCA;
y_GCA = x_GCA;
[X_GCA,Y_GCA]=ndgrid(x_GCA,y_GCA);
F_interp = griddedInterpolant(X,Y,T(:,:,1));%temperatur interpolation
T_GCA = F_interp({x_GCA,y_GCA});%temperature profile for GCA

%%%%piecewise viscosity%%%%
tmp = log10(eta_inf)+(12-log10(eta_inf))*Tg./Tglass.*exp((m/(12-log10(eta_inf))-1)*(Tg./Tglass-1));
eta_Tglass =10.^tmp;
eta0 = eta_Tglass*exp(-Ea/kB_eV/Tglass);
eta1 = eta0*exp(Ea/kB_eV./T_GCA);%T < Tglass
tmp = log10(eta_inf)+(12-log10(eta_inf))*Tg./T_GCA.*exp((m/(12-log10(eta_inf))-1)*(Tg./T_GCA-1));
eta2 = 10.^tmp;%T >= Tglass
eta1(T_GCA>=Tglass) = 0;eta2(T_GCA<Tglass) = 0;
eta = eta1+eta2;%for all T

%%%%nucleation%%%%
Dg = vm*DHf*(Tm-T_GCA)/Tm*2.*T_GCA./(Tm+T_GCA);%bulk Gibbs energy difference
nc = f*32*pi/3*vm^2*sigma^3./Dg.^3;%critical number
DGc = f*16*pi/3*vm^2*sigma^3./Dg.^2;%critical energy barrier
gamma = kB_J*T_GCA/3/pi/lambdaj^3./eta;%jump frequency
Z = sqrt(Dg/6/pi/kB_J./T_GCA./nc);%Zeldovich factor

%%%%Probabilities%%%%
Pgr = 4*r_atom*kB_J*T_GCA/3/pi/Rh/lambdaj.^2/D./eta.*(1);%growth
Pdi = 4*r_atom*kB_J*T_GCA/3/pi/Rh/lambdaj.^2/D./eta.*(exp(-Dg/kB_J./T_GCA));%dissociation
Pnu = 4*gamma.*nc.^(2/3).*Z.*exp(-DGc/kB_J./T_GCA).*nc*(D/d_GST)^2;%nucleation
Pnu(T_GCA>=Tm) = 0;%nucleation probability = 0 for T>=Tm

%%%%number of neighbors%%%%
N_am = nij*ones(N,N);%number of amorphous neighbors
N_cr = zeros(N,N);%number of crystalline neighbors
N_phi = zeros(N,N);%%number of crystalline neighbors belongs to same grain

%step 1 of GCA
t_total_GCA = 0;
r = zeros(N,N);%rij of as-deposited amorphous state
phi = pi*rand(N,N);%random distribution of phi_ij

%step 2 of GCA
%all sites have the same rate coefficients
Cnu = ones(N,N)*Pnu*8;
Cgr = zeros(N,N);
Cdi = zeros(N,N);

%%%%Calculate initial crystalized fraction%%%%
mask = r*0;
mask(sqrt((Y_GCA-Ly/2).^2+(X_GCA-Lx/2).^2)<=sqrt(2*log(2))*w0/2)=1; %only the sites within laser FWHM
N_mask = sum(mask(:));
Xf = sum(sum(r.*mask))/N_mask;
Xf_save(1) = Xf;
t_Xf(1) = t_total_GCA ;

%% %%%Effective Medium Theory%%%
right = Xf*(epsilon_c-1)/(epsilon_c+2)+(1-Xf)*(epsilon_a-1)/(epsilon_a+2);
epsilon_eff = (1+2*right)./(1-right);
N_eff = sqrt(epsilon_eff);
n_eff = real(N_eff);k_eff = -imag(N_eff);
alpha_eff = 4*pi*k_eff/wl;%new absorption coefficient

%% %%%Fresnel Equations%%%
%incident angle = 0
psi_SiN = 2*pi/wl*N_SiN*Lz_SiN;
M_SiN = [cos(psi_SiN) (1i/N_SiN)*sin(psi_SiN);1i*N_SiN*sin(psi_SiN) cos(psi_SiN)];
psi_GST = 2*pi/wl*N_eff*Lz_GST;
M_GST = [cos(psi_GST) (1i/N_eff)*sin(psi_GST);1i*N_eff*sin(psi_GST) cos(psi_GST)];
M = M_GST*M_SiN;
m11 = M(1,1);m12 = M(1,2);m21 = M(2,1);m22 = M(2,2);
refl = (m11-m22+m12-m21)/(m11+m22+m12+m21);
trans = 2/(m11+m22+m12+m21);
Refl = abs(refl).^2;
Trans = abs(trans)^2;
Reflectivity(1) = Refl;Transmissivity(1) = Trans;%save

%% display
load('mymap.mat','mymap');
mapflag = 1;%select colormap
figure,
hImage1 = imagesc((x_GCA-Lx/2)*1e6,(y_GCA-Lx/2)*1e6,r);
colormap([0 0 0.52]);axis equal;
axis([-N*h_GCA/2 N*h_GCA/2 -N*h_GCA/2 N*h_GCA/2]*1e6);
xlabel('X direction (um)');ylabel('Y direction (um)');
title('Time = 0 ns');
set(gca,'FontSize',12);
drawnow;

%% %%%main loop%%%
%%
t_total_heat = 0;%total calculation time
dtau =0 ; %GCA time step
count = 1;
t_show = 1e-9;%show results every 1ns

while(1)
    %time rotation strategy
    dt = dt_select*dt1+abs(dt_select-1)*dt2;
    dt_select = mod(dt_select+1,2);
    t_total_heat= t_total_heat+dt;
   %% heat calculation
    Q = I*(1-Refl)*wfm(t_total_heat);%laser heat source
    
    %GST top layer
    i = 2:length(x)-1;  j = 2:length(y)-1;  k = 1;
    Temp(i,j,k)=T(i,j,k)+alpha_GST_a*dt*2*((dX(i-1,j,k).*(T(i+1,j,k)-T(i,j,k))-dX(i,j,k).*(T(i,j,k)-T(i-1,j,k)))./(dX(i,j,k).*dX(i-1,j,k))./(dX(i,j,k)+dX(i-1,j,k))+...
        (dY(i,j-1,k).*(T(i,j+1,k)-T(i,j,k))-dY(i,j,k).*(T(i,j,k)-T(i,j-1,k)))./(dY(i,j,k).*dY(i,j-1,k))./(dY(i,j,k)+dY(i,j-1,k))+...
        (dZ(i,j,k).*(T(i,j,k+1)-T(i,j,k))-dZ(i,j,k).*(T(i,j,k)-T(i,j,k+1)))./(dZ(i,j,k).*dZ(i,j,k))./(dZ(i,j,k)+dZ(i,j,k)))+dt/Cp_GST*Q(i,j)*2/h*(exp(-alpha_eff*z(k))-exp(-alpha_eff*(z(k)+h/2)));
    
    %GST inner layer
    i = 2:length(x)-1;  j = 2:length(y)-1;  
    for k = 2:find(z>=Lz_GST,1)-1
        Temp(i,j,k)=T(i,j,k)+alpha_GST_a*dt*2*((dX(i-1,j,k).*(T(i+1,j,k)-T(i,j,k))-dX(i,j,k).*(T(i,j,k)-T(i-1,j,k)))./(dX(i,j,k).*dX(i-1,j,k))./(dX(i,j,k)+dX(i-1,j,k))+...
            (dY(i,j-1,k).*(T(i,j+1,k)-T(i,j,k))-dY(i,j,k).*(T(i,j,k)-T(i,j-1,k)))./(dY(i,j,k).*dY(i,j-1,k))./(dY(i,j,k)+dY(i,j-1,k))+...
            (dZ(i,j,k-1).*(T(i,j,k+1)-T(i,j,k))-dZ(i,j,k).*(T(i,j,k)-T(i,j,k-1)))./(dZ(i,j,k).*dZ(i,j,k-1))./(dZ(i,j,k)+dZ(i,j,k-1)))+dt/Cp_GST*Q(i,j)/h*(exp(-alpha_eff*(z(k)-h/2))-exp(-alpha_eff*(z(k)+h/2)));
    end
    
    %interface between SiN and GST
    i = 2:length(x)-1;  j = 2:length(y)-1;  k = find(z>=Lz_GST,1);
    Temp(i,j,k)=(k_SiN*dZ(i,j,k-1).*T(i,j,k+1)+k_GST_a*dZ(i,j,k).*T(i,j,k-1))./(dZ(i,j,k)*k_GST_a+dZ(i,j,k-1)*k_SiN)+dt/Cp_GST*Q(i,j)*2/(dz(k)+dz(k-1))*(exp(-alpha_eff*(z(k)-h/2))-exp(-alpha_eff*(z(k))));
     
    %SiN inner layer
    i = 2:length(x)-1;  j = 2:length(y)-1;  k = (find(z>=Lz_GST,1)+1):(length(z)-1);    
    Temp(i,j,k)=T(i,j,k)+alpha_SiN*dt*2*((dX(i-1,j,k).*(T(i+1,j,k)-T(i,j,k))-dX(i,j,k).*(T(i,j,k)-T(i-1,j,k)))./(dX(i,j,k).*dX(i-1,j,k))./(dX(i,j,k)+dX(i-1,j,k))+...
        (dY(i,j-1,k).*(T(i,j+1,k)-T(i,j,k))-dY(i,j,k).*(T(i,j,k)-T(i,j-1,k)))./(dY(i,j,k).*dY(i,j-1,k))./(dY(i,j,k)+dY(i,j-1,k))+...
        (dZ(i,j,k-1).*(T(i,j,k+1)-T(i,j,k))-dZ(i,j,k).*(T(i,j,k)-T(i,j,k-1)))./(dZ(i,j,k).*dZ(i,j,k-1))./(dZ(i,j,k)+dZ(i,j,k-1)));
    
    %SiN bottom layer
    i = 2:length(x)-1;  j = 2:length(y)-1;  k = length(z);    
    Temp(i,j,k)=T(i,j,k)+alpha_SiN*dt*2*((dX(i-1,j,k).*(T(i+1,j,k)-T(i,j,k))-dX(i,j,k).*(T(i,j,k)-T(i-1,j,k)))./(dX(i,j,k).*dX(i-1,j,k))./(dX(i,j,k)+dX(i-1,j,k))+...
        (dY(i,j-1,k).*(T(i,j+1,k)-T(i,j,k))-dY(i,j,k).*(T(i,j,k)-T(i,j-1,k)))./(dY(i,j,k).*dY(i,j-1,k))./(dY(i,j,k)+dY(i,j-1,k))+...
        (dZ(i,j,k-1).*(T(i,j,k-1)-T(i,j,k))-dZ(i,j,k).*(T(i,j,k)-T(i,j,k-1)))./(dZ(i,j,k-1).*dZ(i,j,k))./(dZ(i,j,k-1)+dZ(i,j,k)));
          
    %update temperature
    T=Temp;
    fprintf('Time is %5.2f ns\n',t_total_heat*1e9);
   %% GCA 
    %%%%re-sample temperature%%%%
    F_interp.Values = T(:,:,1);
    T_GCA = F_interp({x_GCA,y_GCA});
    
    %%%%piecewise viscosity%%%%
    eta1 = eta0*exp(Ea/kB_eV./T_GCA);%T < Tglass
    tmp = log10(eta_inf)+(12-log10(eta_inf))*Tg./T_GCA.*exp((m/(12-log10(eta_inf))-1)*(Tg./T_GCA-1));
    eta2 = 10.^tmp;%T >= Tglass
    eta1(T_GCA>=Tglass) = 0;eta2(T_GCA<Tglass) = 0;
    eta = eta1+eta2;

    %%%%%nucleation%%%%
    Dg = vm*DHf*(Tm-T_GCA)/Tm*2.*T_GCA./(Tm+T_GCA);%bulk Gibbs energy difference
    nc = f*32*pi/3*vm^2*sigma^3./Dg.^3;%critical number
    DGc = f*16*pi/3*vm^2*sigma^3./Dg.^2;%critical energy barrier
    gamma = kB_J*T_GCA/3/pi/lambdaj^3./eta;%jump frequency
    Z = sqrt(Dg/6/pi/kB_J./T_GCA./nc);%Zeldovich factor
    
    %%%%Probability%%%%
    Pgr = 4*r_atom*kB_J*T_GCA/3/pi/Rh/lambdaj.^2/D./eta.*(1);
    Pdi = 4*r_atom*kB_J*T_GCA/3/pi/Rh/lambdaj.^2/D./eta.*(exp(-Dg/kB_J./T_GCA));
    Pnu = 4*gamma.*nc.^(2/3).*Z.*exp(-DGc/kB_J./T_GCA).*nc*(D/d_GST)^2;
    Pnu(T_GCA>=Tm) = 0;%nucleation probability = 0 for T>=Tm

    %%%%rates C%%%%
    Cnu = N_am.*Pnu; Cnu(r==1) = 0;
    Cgr = N_cr.*Pgr; Cgr(r==1) = 0;
    Cdi = (nij-N_phi).*Pdi; Cdi(r==0) = 0;
    
    %call GCA 
    GCA;
    
    %%%%crystallized fraction%%%%
    Xf = sum(sum(r.*mask))/N_mask;
    
    %%%%display crystal microstructures%%%%
    if t_total_heat>=t_show
        t_show = t_show+1e-9;
        c=phi; c(r==0)=0;
        set(hImage1,'cdata',c);title(sprintf('Time = %.2f ns',t_total_heat*1e9));
        if Xf>0 && mapflag
            colormap(mymap);
            mapflag = 0;
        elseif Xf == 0 && mapflag ==0
            colormap([0 0 0.52])
            mapflag = 1;
        end
        drawnow;
    end
   %% %%%Effective Medium Theory%%%
    right = Xf*(epsilon_c-1)/(epsilon_c+2)+(1-Xf)*(epsilon_a-1)/(epsilon_a+2);
    epsilon_eff = (1+2*right)./(1-right);
    N_eff = sqrt(epsilon_eff);
    n_eff = real(N_eff);k_eff = -imag(N_eff);
    alpha_eff = 4*pi*k_eff/wl;
    
   %% %%%Fresnel Equations%%%
    psi_GST = 2*pi/wl*N_eff*Lz_GST;
    M_GST = [cos(psi_GST) (1i/N_eff)*sin(psi_GST);1i*N_eff*sin(psi_GST) cos(psi_GST)];
    M = M_GST*M_SiN;
    m11 = M(1,1);m12 = M(1,2);m21 = M(2,1);m22 = M(2,2);
    refl = (m11-m22+m12-m21)/(m11+m22+m12+m21);
    trans = 2/(m11+m22+m12+m21);
    Refl = abs(refl).^2;
    Trans = abs(trans)^2;     
    
   %% %%%save data to array%%%
    if t_total_heat>=t_save
        count = count+1;
        t_save = t_save + 1e-10;
        Tmax_heat(count) = max(T(:));
        t_heat(count) = t_total_heat;
        Xf_save(count) = Xf;
        t_Xf(count) = t_total_GCA ;
        Reflectivity(count) = Refl;
        Transmissivity(count) = Trans;
    end
    %% %%%if finished%%%
    if t_total_heat>=t_final
        break;
    end
end
%% %%%save data to file after finish%%%
save(sprintf('Pin = %f mW.mat',Pin*1e3),'t_heat','Tmax_heat','t_Xf','Xf_save','Reflectivity','Transmissivity');
disp('More information can be found in this paper:https://arxiv.org/abs/2107.02035.');
disp('Please cite this paper if you fully or partially use it in your paper.');

%% waveform function
% t: time; p: intensity ration [0 1]
% isosceles trapezoid
function p = wfm(t)
    if t<edgeWidth
        p = 1/edgeWidth*t;
    elseif t<pulseWidth+edgeWidth
        p = 1;
    elseif t<pulseWidth+2*edgeWidth
        p = 1-(t-pulseWidth-edgeWidth)/edgeWidth;
    else
        p = 0;
    end
end
   
%% GCA function
function GCA
while (1)
    %step 3 compute total probabilities
    a0 = sum(Cnu(:)+Cgr(:)+Cdi(:));%sum of all rates
    
    %step 4 determine the increment time
    n1 = rand(); n2 = rand();%two random numbers
    dtau = 1/a0*log(1/n1);%GCA time step
    t_total_GCA = t_total_GCA+dtau;
    if t_total_GCA>t_total_heat %if GCA time > heat time, quit
        t_total_GCA = t_total_GCA-dtau;
        break;
    end
    
    % step 5 identify the event and site (i,j)
    s=0;%sum
    a=0;%which operation
    for ii = 1:N
        for jj = 1:N
            s=s+Cnu(ii,jj);
            if s>=n2*a0%if nucleation
                a=1;
                break;
            end
            s=s+Cgr(ii,jj);
            if s>=n2*a0%if growth
                a=2;
                break;
            end
            s=s+Cdi(ii,jj);
            if s>=n2*a0%if dissociation
                a=3;
                break;
            end
        end
        if a ~= 0
            break;
        end
    end
    
    % step 6 update r & phi
    switch a
        case 1%nucleation
            r(ii,jj) = 1;
        case 2%growth
            [kk,ll] = grow(ii,jj);
            r(ii,jj) = 1;
            phi(ii,jj) = phi(kk,ll);
        case 3%dissociation
            r(ii,jj) = 0;
            phi(ii,jj) = pi*rand();
    end
    
    %step 7 update Cnu, Cgr, Cdi
    %only update those sites that have changed and their neighbors
    updateC(ii,jj);
    
    %step 8 return to step 3 
    
end
end
%% nested functions for GCA
    %number of amorphous neighbors, n_am
    function y = n_am(i,j)
        %periodic boundary condition
        y=nij-sum(sum(r([mod(i-2+N,N) mod(i-1+N,N) mod(i+N,N)]+1,[mod(j-2+N,N) mod(j-1+N,N) mod(j+N,N)]+1)))+r(i,j);
        N_am(i,j) = y;
    end
    %number of crystalline neighbors at orientation phi (for dissociation)
    function y = n_phi(i,j)
        t1 = phi([mod(i-2+N,N) mod(i-1+N,N) mod(i+N,N)]+1,[mod(j-2+N,N) mod(j-1+N,N) mod(j+N,N)]+1);
        y = length(find(t1 == t1(2,2)))-1;
        N_phi(i,j) = y;
    end
    %number of all crystalline neighbors (for growth)
    function y = n_gr(i,j)
        y=sum(sum(r([mod(i-2+N,N) mod(i-1+N,N) mod(i+N,N)]+1,[mod(j-2+N,N) mod(j-1+N,N) mod(j+N,N)]+1)))-r(i,j);
        N_cr(i,j) = y;
    end       
    %detemine which grain the site will be grow onto
    function [k,l] = grow(i,j)
        %check the number of crystalline neighbors and randomly grow onto one neighbors
        n_cr = n_gr(i,j);
        loc = ceil(n_cr*rand());
        s1 = 0;
        for k = [mod(i-2+N,N) mod(i-1+N,N) mod(i+N,N)]+1
            for l = [mod(j-2+N,N) mod(j-1+N,N) mod(j+N,N)]+1
                if k == i && l == j
                    continue;
                end
                if r(k,l) == 1
                    s1 = s1+1;
                    if s1 == loc
                        return
                    end
                end
            end
        end        
    end    
    %update rates Cnu, Cgr, Cdi of one site and its neighbors
    function updateC(i,j)
        for row = [mod(i-2+N,N) mod(i-1+N,N) mod(i+N,N)]+1
            for col = [mod(j-2+N,N) mod(j-1+N,N) mod(j+N,N)]+1
                n_am(row,col);n_gr(row,col);n_phi(row,col);
                if r(row,col) == 0
                    Cnu(row,col) = n_am(row,col)*Pnu(row,col);
                    if n_gr(row,col)==0
                        Cgr(row,col) = 0;
                    else
                        Cgr(row,col) = n_gr(row,col)*Pgr(row,col);
                    end
                    Cdi(row,col) = 0;
                else
                    Cnu(row,col) = 0;
                    Cgr(row,col) = 0;
                    Cdi(row,col) = (nij-n_phi(row,col))*Pdi(row,col);
                end
            end
        end        
    end
end
