clearvars
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

% Io parameters
G = 6.67430E-11;
R_io = 1821.6; % Radius in km !!!!!!!!!!!!!
M_io = 8.931938e22; % mass in kg
rho_avg_io = 3528; % avg density in kg/m3
omega_io = 4.1086E-05; % orbital frequency in rad/s
T_io = 2 * pi / omega_io; % orbital period in s

% % Interior Model --- IN DIMENSIONAL UNITS
eta_avg = 1e19; % avg shear viscosity in Pa.s
mu_avg = 60e9; % avg shear modulus in Pa
rho_avg = 3263; % avg density in kg/m3

% Core boundary layer (1)
Interior_Model(1).R0 = 1e-3 * R_io; % approx 0
Interior_Model(1).rho0 = rho_avg;

% Single body layer (2)
Interior_Model(2).R0 = R_io;
Interior_Model(2).rho0 = rho_avg;
Interior_Model(2).mu0 = mu_avg;
Interior_Model(2).eta0 = eta_avg;

% LATERAL VARIATIONS
% Spherically-symmetric model assumed
% compute effective shear modulus
mu_eff = mu_avg / ( rho_avg ^ 2 * ( R_io * 1e3 ) ^ 2 * 4 / 3 * pi * G );

% forcing
Forcing(1).Td = T_io;
Forcing(1).n = 2; 
Forcing(1).m = 0; 
Forcing(1).F = 1;

% radial discretization
Numerics.Nlayers = length(Interior_Model); % number of concentric layers. Including the core!
Numerics.method = 'variable'; % method of setting the radial points per layer
Numerics.Nrbase = 200; % depending on the method this will determine the number of points per layer
% code parallelization
Numerics.parallel_sol = 0; % Use a parfor-loop to call get_Love, either 0 or 1
Numerics.parallel_gen = 0; % Calculate potential coupling files and the propagation inside get_solution using parfor-loops, either 0 or 1
% lateral variations
Numerics.perturbation_order = 2; %maximum order to which couplings are considered
Numerics.solution_cutoff = 12; % maximum degree of solution, not used if perturbation order is specified
Numerics.load_couplings = 1; % 0=no loading, 1=loading of specific file, 2=searches for big enough file
Numerics.Nenergy = 12; % maximum degree to which energy dissipation is expanded 
Numerics.rheology_cutoff = 2; % maximum order of difference (so in log) up to which rheology is still used 
[Numerics, Interior_Model] = set_boundary_indices(Numerics, Interior_Model,'verbose');
%% Compute the tidal response 
% The tidal response is obtained 

% write interior model in the right format for the code
Interior_Model = get_rheology(Interior_Model, Numerics, Forcing);
[Love_Spectra, y] = get_Love(Interior_Model, Forcing, Numerics,'verbose');

% obtain the Fourier-transformed effective shear modulus
mu_eff_hat = mu_eff * Interior_Model(2).muC;
% compute Love numbers using analytical expression
n = Forcing.n; 
mu_n = ( 2 * n ^ 2 + 4 * n + 3 ) / n * mu_eff_hat;
k2_analytic = 1 / ( 1 + mu_n ) * 3 / 2 / ( n - 1 );
h2_analytic = 1 / ( 1 + mu_n ) * ( 2 * n + 1 ) / 2 / ( n - 1 );
% display results
disp(['k Love number analytical expression: ' num2str(k2_analytic)])
disp(['k Love number LOV3D: ' num2str(Love_Spectra.k)])
disp(['Normalized difference: ' num2str((Love_Spectra.k-k2_analytic)/k2_analytic*100)  '%'])
disp(['h Love number analytical expression: ' num2str(h2_analytic)])
disp(['h Love number LOV3D: ' num2str(Love_Spectra.h)])
disp(['Normalized difference: ' num2str((Love_Spectra.h-h2_analytic)/h2_analytic*100)  '%'])