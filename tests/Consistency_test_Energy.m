clear all
clc
%% Test correct function of the code by comparing with some previous results
addpath(genpath('/home/allard/Code/New_LOV3D/lov3d/'))
set(0,'defaulttextInterpreter','latex') 
mfile_name          = mfilename('fullpath');
[pathstr,name,ext]  = fileparts(mfile_name);
cd(pathstr);
cd('..')
addpath(genpath(pwd))

% NOTHING WILL BE SAVED WHEN RUNNING THIS UNIT TEST!
% Location where any data will be saved, needs to end with /
save_location = 'files/output/test_internal_consistency/';
if ~isfolder(save_location)
    mkdir(save_location)
end
Startall = tic;

%% NOTICE!
% This test can take a long time to complete and should therefore be used
% with caution!

%% Loop stuff 
Nbase_array = [5,10,20,50,100,200,500,1000];
N_iteration = length(Nbase_array);
k2_array = zeros(1,N_iteration);
e00_array = zeros(1,N_iteration);
e_01_Uni = zeros(1,N_iteration);
e_02_Uni = zeros(1,N_iteration);
e_01 = zeros(1,N_iteration);
e_02 = zeros(1,N_iteration);
Nr_array = zeros(1,N_iteration);

for ii=1:N_iteration

disp(['Starting iteration: ' num2str(ii)])

%% Starting-pattern
% Maximum number of degrees (related to number of gridpoints)
l_max = 100;

% Lat-lon grid creation
delta_lon=180/(2*(2*l_max-1)); %nominal as delta_lon
delta_lat=180/(2*(2*l_max-1)); 
lon0=-180+delta_lon/2; 
lonM=180-delta_lon/2;
lat0=-90+delta_lat/2; 
latM=90-delta_lon/2;

lon=lon0:delta_lon:lonM;
lat=lat0:delta_lat:latM; 
colat=90-lat;

colat=colat*pi/180;
lon=lon*pi/180;
P_lm=legendre(2,cos(colat)');
gv=zeros(length(colat),length(lon));
phase = 0;

for m=1:3 %loop in order
     % (-1) added for consistency with definition in my paper
     % % (-1)^m added bc I am using P_l^m and not P_{l,m}
    if m == 1
        fac = -(33/7);
    elseif m == 3
        fac = (9/14);
    else 
        fac = 0;
    end

    v1=(-1)^(m-1)*fac*squeeze(P_lm(m,:));

    if m == 1
        for i=1:length(v1)
            gv(i,:) = gv(i,:) + v1(i);
        end
    else
        for i=1:length(v1)
            gv(i,:) = gv(i,:) + v1(i)*( cos(phase)*cos((m-1)*lon) + sin(phase)*sin((m-1)*lon) );
        end
    end
end
z_val = (21/5 + 0.5*gv)/(21/5);

lon=lon*180/pi; 
lon0=min(lon);
lonM=max(lon);
lat0=min(lat);
latM=max(lat);
[lon_g,lat_g] = meshgrid(lon0:delta_lon:lonM,lat0:delta_lat:latM);

Psi_latlon.lon=lon_g;
Psi_latlon.lat=lat_g;
Psi_latlon.z=z_val; 
Psi_latlon.lmax = 2*l_max-1;

% Heat dissipation
Q_mean = 2.3; % W m^-2
Q_prod_latlon = Psi_latlon.z * Q_mean;
Q_diff_latlon = Q_prod_latlon - Q_mean;

% Melt fraction
c = 0.01; % D2-0 model
Phi_mean = 0.1;
Phi_diff_latlon = c*Q_diff_latlon;

% Viscosity
B_eta = 20;
eta_start_latlon = exp(-B_eta*Phi_diff_latlon);

% Shear modulus
B_mu = 67/15;
mu_start_latlon = (1+B_mu*Phi_mean)./(1+B_mu*(Phi_diff_latlon + Phi_mean));

% Convert lat_lon to SPH
mu_start.lon = Psi_latlon.lon;
mu_start.lat = Psi_latlon.lat;
mu_start.z = mu_start_latlon;
mu_start.lmax = Psi_latlon.lmax;

eta_start.lon = Psi_latlon.lon;
eta_start.lat = Psi_latlon.lat;
eta_start.z = eta_start_latlon;
eta_start.lmax = Psi_latlon.lmax;

%% Io non-dimensional
% Values will be non-dimensionalised based on the outer layer such that all
% radii are fractions. This approached might need to be changed when a
% dedicated crust is added as that will likely not be modeled.
G=6.67e-11;

% All values are taken from Steinke et al. (2020a)
% non-dimensional parameters for Io (see Section 2.2 and Table 1)
omega0=4.1086E-05; %Io's orbital frequency 
T=2*pi/omega0; %Io's orbital period
mu_eff=5.2; % mu_mean/(gs*rho_mean*R) gs=gravity acceleration at the surface
R_Io = 1821.6; % Radius of Io in km
R_CMB = 965; % Radius of Io's core

% Interior_Model goes from core (1) to surface (n) with indices
% Interior_Model values for the core
Interior_Model(1).R0 = R_CMB; % CMB in km 
Interior_Model(1).rho0 = 5150; % average density in kg m^-3
Interior_Model(1).rho0_2 = 5150; % rho_2 in kg m^-3 (used for icy moons)

% Interior_Model values for mantle/asthenosphere/crust layers
% Number of layers on top of the core should be equal to Numerics.Nlayers
% If the layer is elastic eta should be NaN
Interior_Model(2).R0 = 1591.6; % outer radius of this layer
Interior_Model(2).rho0 = 3244; % density in kg m^-3
Interior_Model(2).Ks0 = 200e12; % bulk modulus in Pa
Interior_Model(2).mu0 = 6e10; % shear modulus in Pa
Interior_Model(2).eta0 = 1e20; % viscosity in Pa s

Interior_Model(3).R0 = 1791.6; % outer radius of this layer
Interior_Model(3).rho0 = 3244; % density in kg m^-3
Interior_Model(3).Ks0 = 200e12; % bulk modulus in Pa
Interior_Model(3).mu0 = 7.8e5; % shear modulus in Pa
Interior_Model(3).eta0 = 1e11; % viscosity in Pa s
Interior_Model(3).mu_latlon = mu_start; % Full spectrum (struct) as input
Interior_Model(3).eta_latlon = eta_start; % Full spectrum (struct) as input

Interior_Model(4).R0 = R_Io; % outer radius of this layer
Interior_Model(4).rho0 = 3244; % density in kg m^-3
Interior_Model(4).Ks0 = 200e12; % bulk modulus in Pa
Interior_Model(4).mu0 = 6.5e10; % shear modulus in Pa
Interior_Model(4).eta0 = 1e23; % viscosity in Pa s

% Set the rho difference such that the core can also represent an ocean
Interior_Model(1).Delta_rho0 = Interior_Model(1).rho0_2 - Interior_Model(2).rho0; 

%% NUMERICS
Numerics.Nlayers = 4; % number of concentric layers. Including the core!
Numerics.method = 'combination'; % method of setting the radial points per layer: combination, variable, fixed and manual
Numerics.Nrbase = Nbase_array(ii); % depending on the method this will determine the number of points per layer
Numerics.perturbation_order = 2; %maximum order to which couplings are considered
Numerics.solution_cutoff = 12; % maximum degree of solution, not used if perturbation order is specified
Numerics.load_couplings = 1; % 0=no loading, 1=loading of specific file, 2=searches for big enough file
Numerics.Nenergy = 12; % maximum degree to which energy dissipation is expanded 
Numerics.rheology_cutoff = 2; % maximum difference (in log) up to which rheology is still used 
Numerics.minimum_rheology_value = -13; % minimum value for the exponent of a rheology term
Numerics.parallel_sol = 0; % Use a parfor-loop to call get_Love, either 0 or 1
Numerics.parallel_gen = 0; % Calculate potential coupling files and the propagation inside get_solution using parfor-loops, either 0 or 1
Numerics.coupling_file_location = 'files/couplings/'; % MUST END WITH A '/' , Location where coupling files are saved or should be saved 

[Numerics, Interior_Model] = set_boundary_indices(Numerics, Interior_Model,'verbose');
Nr_array(ii) = Numerics.Nr;

%% Small but necessary additions/checks of the code
% Make sure that Numerics.Nlayers is equal to length of Interior_Model
if ~(Numerics.Nlayers == length(Interior_Model))
    error('Error. \nSize of Interior_Model struct must be equal to Numerics.Nlayers')
end

%% FORCING 
% Forcing corresponding to a synchronous moon in an eccentric orbit (see Appendix D)
Forcing(1).Td=2*pi/omega0;
Forcing(1).n=2; 
Forcing(1).m=0; 
Forcing(1).F=3/4*sqrt(1/5); 
Forcing(2).Td=2*pi/omega0;
Forcing(2).n=2; 
Forcing(2).m=-2; 
Forcing(2).F=-7/8*sqrt(6/5);
Forcing(3).Td=2*pi/omega0;
Forcing(3).n=2; 
Forcing(3).m=2; 
Forcing(3).F=1/8*sqrt(6/5);

%% Convert the real lateral changes into complex spherical harmonics
% Calculates the rheology that will be used. First computes the complex 
% spherical harmonics and then performs a fourier transform to obtain the 
% coefficients that will be used throughout the computations 
Interior_Model_Uni = Interior_Model;
for ilayer=2:Numerics.Nlayers
    if isfield(Interior_Model_Uni(ilayer), "mu_latlon")
        Interior_Model_Uni = rmfield(Interior_Model_Uni, "mu_latlon");
    end
    if isfield(Interior_Model_Uni(ilayer), "eta_latlon")
        Interior_Model_Uni = rmfield(Interior_Model_Uni, "eta_latlon");
    end
    if isfield(Interior_Model_Uni(ilayer), "k_latlon")
        Interior_Model_Uni = rmfield(Interior_Model_Uni, "k_latlon");
    end
end

Interior_Model = get_rheology(Interior_Model,Numerics,Forcing);

% prepares a spherically-symmetric model for comparison   
Interior_Model_Uni = get_rheology(Interior_Model_Uni,Numerics,Forcing);

%% GET LOVE NUMBERS
%obtains the Love number spectra

for i=1:length(Forcing)
    [Love_Spectra_Uni(i),y_Uni(i)]=get_Love(Interior_Model_Uni,Forcing(i),Numerics, ...                                            'verbose');
    % lateral variations
    [Love_Spectra(i),y(i)]=get_Love(Interior_Model,Forcing(i),Numerics,'verbose');
end

%% GET ENERGY SPECTRA, can be commentted if elastic
[Energy_Spectra]=get_energy(y,Numerics,Forcing,Interior_Model,'verbose',1,'calc_E_contribution',0);
[Energy_Spectra_Uni]=get_energy(y_Uni,Numerics,Forcing,Interior_Model_Uni,'verbose',1,'calc_E_contribution',0);

E_k_Uni=0;
E_k=0;
for i=1:length(Forcing)
    n_f=Forcing(i).n;
    m_f=Forcing(i).m;
    for j=1:length(Forcing)
        ind1=find(Love_Spectra_Uni(j).n==n_f & Love_Spectra_Uni(j).m==m_f);
        ind2=find(Love_Spectra(j).n==n_f & Love_Spectra(j).m==m_f);
        if isempty(ind1)==0
            E_k_Uni=E_k_Uni-Forcing(i).F*Forcing(j).F*imag(Love_Spectra_Uni(j).k(ind1));
        end
        if isempty(ind2)==0
            E_k=E_k-Forcing(i).F*Forcing(j).F*imag(Love_Spectra(j).k(ind2));
        end
    end
end
% Uni
e_01_Uni(ii)=2*pi*10/(4*pi*Interior_Model_Uni(end).Gg)*E_k_Uni;
e_02_Uni(ii)=Energy_Spectra_Uni.energy_integral(1);
e_01(ii)=2*pi*10/(4*pi*Interior_Model(end).Gg)*E_k;
e_02(ii)=Energy_Spectra.energy_integral(1);

for i=1:length(Love_Spectra(1).n)
    if Love_Spectra(1).n(i)==Forcing(1).n && Love_Spectra(1).m(i)==Forcing(1).m
        k2 = y(1).y(end,8,i)-1;
    end
end

k2_array(ii) = k2;
e00_array(ii) = Energy_Spectra.energy_integral(1);

clear Interior_Model Interior_Model_Uni Numerics Love_Spectra Love_Spectra_Uni
clear Forcing Energy_Spectra Energy_Spectra_Uni eta_start mu_start

disp(["Time spent on all calculations" num2str(toc(Startall)) " s"])

end
%%
set(0,'defaulttextInterpreter','latex') 
fig = figure(1);
% fontsize(fig, 24, 'points')
ratio_Uni_variable = real(abs(e_01_Uni - e_02_Uni)) ./ e_01_Uni * 100;
ratio_variable = real(abs(e_01 - e_02)) ./ e_01 * 100;
loglog(Nbase_array,ratio_Uni_variable)
hold on
loglog(Nbase_array,ratio_variable,'--')
hold off
grid on
legend('Spherically symmetric', 'Lateral variations','Location','best')
xlabel('Total number of radial points')
ylabel('$ |E_1 - E_2|/ E_1 \cdot 100 $ [\%]')
title('Diff. between calculating energy using integrated stress or love numbers')

