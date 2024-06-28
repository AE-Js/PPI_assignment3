clear all
clc

set(0,'defaulttextInterpreter','latex') 
mfile_name          = mfilename('fullpath');
if contains(mfile_name,'LiveEditorEvaluationHelper')
    mfile_name=matlab.desktop.editor.getActiveFilename;
end
[pathstr,name,ext]  = fileparts(mfile_name);
cd(pathstr);
cd('..')
addpath(genpath(pwd))

% Location where any data will be saved, needs to end with /
save_location = 'files/output/multi_layer_example/';
if ~isfolder(save_location)
    mkdir(save_location)
end
Startall = tic;

%% Io parameters
G = 6.67430E-11;
r_l = 1821.6; % Lithosphere radius in [km] 
r_a = 1791.6; % Asthenosphere radius in [km] 
r_m = 1741.6; % Mantle radius in [km] 
r_c = 980; % Core radius in [km] 

rho_l = 2700; % Lithosphere density [kg/m³]
rho_a = 3000; % Asthenosphere density [kg/m³]
rho_m = 3200; % Mantle density [kg/m³]
rho_c = 5150; % Core density [kg/m³]

eta_l = 1e23; % Lithosphere shear viscosity [Pa s]
eta_a = 1e23; % Asthenosphere shear viscosity [Pa s]
eta_m = 1e23; % Mantle shear viscosity [Pa s]
eta_c = 0.0; % Core shear viscosity [Pa s]

mu_l = 6.5e10; % Lithosphere shear modulus [Pa s]
mu_a = 6e10; % Asthenosphere shear modulus [Pa s]
mu_m = 6e10; % Mantle shear modulus [Pa s]
mu_c = 0.0; % Core shear modulus [Pa s]

K_l = 200e9; % Lithosphere bulk modulus [Pa]
K_a = 200e9; % Asthenosphere bulk modulus [Pa]
K_m = 200e9; % Mantle bulk modulus [Pa]
K_c = 200e9; % Core bulk modulus [Pa]

M_io = 8.931938e22; % mass in kg
omega_io = 4.1086E-05; % orbital frequency in rad/s
T_io = 2 * pi / omega_io; % orbital period in s

% % % % % % % % % % % % % %
% Layered Elastic Model % %
% % % % % % % % % % % % % %

% Core
Interior_Model_A(1).R0 = r_c;
Interior_Model_A(1).rho0 = rho_c; 
Interior_Model_A(1).Ks0 = K_c;
Interior_Model_A(1).mu0 = mu_c;


% Mantle
Interior_Model_A(2).R0 = r_m;
Interior_Model_A(2).rho0 = rho_m;
Interior_Model_A(2).Ks0 = K_m;
Interior_Model_A(2).mu0 = mu_m;


% Asthenosphere
Interior_Model_A(3).R0 = r_a;
Interior_Model_A(3).rho0 = rho_a;
Interior_Model_A(3).Ks0 = K_a;
Interior_Model_A(3).mu0 = mu_a;


% Lithosphere
Interior_Model_A(4).R0 = r_l;
Interior_Model_A(4).rho0 = rho_l;
Interior_Model_A(4).Ks0 = K_l;
Interior_Model_A(4).mu0 = mu_l;


% % % % % % % % % % % % % % % % % % % %
% Layered Viscoelastic Model  % % % % % 
% % % % % % % % % % % % % % % % % % % %

% Core
Interior_Model_B(1).R0 = r_c;
Interior_Model_B(1).rho0 = rho_c; 
Interior_Model_B(1).Ks0 = K_c;
Interior_Model_B(1).mu0 = mu_c;
Interior_Model_B(1).eta0 = eta_c;

% Mantle
Interior_Model_B(2).R0 = r_m;
Interior_Model_B(2).rho0 = rho_m;
Interior_Model_B(2).Ks0 = K_m;
Interior_Model_B(2).mu0 = mu_m;
Interior_Model_B(2).eta0 = eta_m;

% Asthenosphere
Interior_Model_B(3).R0 = r_a;
Interior_Model_B(3).rho0 = rho_a;
Interior_Model_B(3).Ks0 = K_a;
Interior_Model_B(3).mu0 = mu_a;
Interior_Model_B(3).eta0 = eta_a;

% Lithosphere
Interior_Model_B(4).R0 = r_l;
Interior_Model_B(4).rho0 = rho_l;
Interior_Model_B(4).Ks0 = K_l;
Interior_Model_B(4).mu0 = mu_l;
Interior_Model_B(4).eta0 = eta_l;


% % % % % % % % % % % % % % % % % % % %
% Layered Ideal Fluid  Model  % % % % % 
% % % % % % % % % % % % % % % % % % % %

% Core
Interior_Model_C(1).R0 = r_c;
Interior_Model_C(1).rho0 = rho_c; 
Interior_Model_C(1).Ks0 = K_c;
Interior_Model_C(1).mu0 = mu_c;
Interior_Model_C(1).eta0 = eta_c;

% Mantle
Interior_Model_C(2).R0 = r_m;
Interior_Model_C(2).rho0 = rho_m;
Interior_Model_C(2).Ks0 = K_m;
Interior_Model_C(2).mu0 = mu_m;
Interior_Model_C(2).eta0 = eta_m;

% Asthenosphere
Interior_Model_C(3).R0 = r_a;
Interior_Model_C(3).rho0 = rho_a;
Interior_Model_C(3).Ks0 = K_a;
Interior_Model_C(3).mu0 = mu_a;
Interior_Model_C(3).eta0 = eta_a;

% Lithosphere
Interior_Model_C(4).R0 = r_l;
Interior_Model_C(4).rho0 = rho_l;
Interior_Model_C(4).Ks0 = K_l;
Interior_Model_C(4).mu0 = mu_l;
Interior_Model_C(4).eta0 = eta_l;

% % % % % % % % % % % % % %
% Numerics  % % % % % % % %
% % % % % % % % % % % % % %

%radial discretization
Numerics.Nlayers = 4; % number of concentric layers. Including the core!
Numerics.method = 'variable'; % method of setting the radial points per layer, here fixed number of layers
Numerics.Nrbase = 2000; % depending on the method this will determine the number of points per layer
%code parallelization
Numerics.parallel_sol = 0; % Use a parfor-loop to call get_Love, either 0 or 1
Numerics.parallel_gen = 0; % Calculate potential coupling files and the propagation inside get_solution using parfor-loops, either 0 or 1
% lateral variations
Numerics.perturbation_order = 2; %maximum order to which couplings are considered
Numerics.solution_cutoff = 12; % maximum degree of solution, not used if perturbation order is specified
Numerics.load_couplings = 1; % 0=no loading, 1=loading of specific file, 2=searches for big enough file
Numerics.Nenergy = 12; % maximum degree to which energy dissipation is expanded 
Numerics.rheology_cutoff = 2; % maximum order of difference (so in log) up to which rheology is still used 

[Numerics, Interior_Model_A] = set_boundary_indices(Numerics, Interior_Model_A,'verbose');
[Numerics, Interior_Model_B] = set_boundary_indices(Numerics, Interior_Model_B,'verbose');
[Numerics, Interior_Model_C] = set_boundary_indices(Numerics, Interior_Model_C,'verbose');

%% Small but necessary additions/checks of the code
% Make sure that Numerics.Nlayers is equal to length of Interior_Model
if ~(Numerics.Nlayers == length(Interior_Model_A))
    error('Error. \nSize of Interior_Model struct must be equal to Numerics.Nlayers')
end

if ~(Numerics.Nlayers == length(Interior_Model_B))
    error('Error. \nSize of Interior_Model struct must be equal to Numerics.Nlayers')
end

if ~(Numerics.Nlayers == length(Interior_Model_C))
    error('Error. \nSize of Interior_Model struct must be equal to Numerics.Nlayers')
end

% % % % % % % % % % % % % %
% Forcing  % % % % % % % %
% % % % % % % % % % % % % %

% Set-up forcing
Forcing_A(1).Td=2*pi/omega_io;
Forcing_A(1).n=2; 
Forcing_A(1).m=0; 
Forcing_A(1).F=1;

Forcing_B(1).Td=2*pi/omega_io;
Forcing_B(1).n=2; 
Forcing_B(1).m=0; 
Forcing_B(1).F=1;

Forcing_C(1).Td=1e12;
Forcing_C(1).n=2; 
Forcing_C(1).m=0; 
Forcing_C(1).F=1;

% % % % % % % % % % % % % %
% Tidal response  % % % % %
% % % % % % % % % % % % % %

% write interior model in the right format for the code
Interior_Model_A= get_rheology(Interior_Model_A,Numerics,Forcing_A);
[Love_Spectra_A, y_A]=get_Love(Interior_Model_A,Forcing_A,Numerics,'verbose');

Interior_Model_B= get_rheology(Interior_Model_B,Numerics,Forcing_B);
[Love_Spectra_B, y_B]=get_Love(Interior_Model_B,Forcing_B,Numerics,'verbose');

Interior_Model_C= get_rheology(Interior_Model_C,Numerics,Forcing_C);
[Love_Spectra_C, y_C]=get_Love(Interior_Model_C,Forcing_C,Numerics,'verbose');


disp(['k_2 elastic model ' num2str(Love_Spectra_A.k)])
disp(['h_2 elastic model ' num2str(Love_Spectra_A.h)])

disp(['k_2 viscoelastic model ' num2str(Love_Spectra_B.k)])
disp(['h_2 viscoelastic model ' num2str(Love_Spectra_B.h)])

disp(['k_2 ideal fluid model ' num2str(Love_Spectra_C.k)])
disp(['h_2 ideal fluid model ' num2str(Love_Spectra_C.h)])