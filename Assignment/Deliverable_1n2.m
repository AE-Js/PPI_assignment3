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
R_io = 1821.6; % Radius in [km] !!!!!!!!!!!!! WHY MUST THE BE IN KM? THE REST IS ALL BASE SI UNITS. THIS CODE IS STUPID
M_io = 8.931938e22; % mass in kg
rho_avg_io = 3528; % avg density in kg/m3
omega_io = 4.1086E-05; % orbital frequency in rad/s
T_io = 2 * pi / omega_io; % orbital period in s

eta_avg = 1e19; % avg shear viscosity in Pa.s
mu_avg = 60e9; % avg shear modulus in Pa
rho_avg = 3263; % avg density in kg/m3
K_avg = 200e9; % avg bulk modulus in Pa

% % % % % % % % % % % % % %
% General stuff           %
% % % % % % % % % % % % % %

% radial discretization
Numerics.Nlayers = 2; % number of concentric layers. Including the core!
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

mu_eff = mu_avg / ( rho_avg ^ 2 * ( R_io * 1e3 ) ^ 2 * 4 / 3 * pi * G );

% % % % % % % % % % % % % %
% Elastic model : eta = 0 %
% % % % % % % % % % % % % %

% Core boundary layer (1)
Elastic_model(1).R0 = 1e-3 * R_io; % approx 0
Elastic_model(1).rho0 = rho_avg;

% Single body layer (2)
Elastic_model(2).R0 = R_io;
Elastic_model(2).rho0 = rho_avg;
Elastic_model(2).mu0 = mu_avg;
Elastic_model(2).Ks0 = K_avg;

% forcing
Forcing_e(1).Td = T_io;
Forcing_e(1).n = 2; 
Forcing_e(1).m = 0; 
Forcing_e(1).F = 1;

[Numerics, Elastic_model] = set_boundary_indices(Numerics, Elastic_model,'empty');
Elastic_model = get_rheology(Elastic_model, Numerics, Forcing_e);
[Love_Spectra_e, y_e] = get_Love(Elastic_model, Forcing_e, Numerics,'empty');

% obtain the Fourier-transformed effective shear modulus
mu_eff_hat_e = mu_eff * Elastic_model(2).muC;
% compute Love numbers using analytical expression
n_e = Forcing_e.n; 
mu_n_e = ( 2 * n_e ^ 2 + 4 * n_e + 3 ) / n_e * mu_eff_hat_e;
k2_analytic_e = 1 / ( 1 + mu_n_e ) * 3 / 2 / ( n_e - 1 );
h2_analytic_e = 1 / ( 1 + mu_n_e ) * ( 2 * n_e + 1 ) / 2 / ( n_e - 1 );

% display results
disp(['ELASTIC CASE:'])
disp(['k analytical: ' num2str(k2_analytic_e)])
disp(['k LOV3D: ' num2str(Love_Spectra_e.k)])
disp(['Normalized difference: ' num2str((Love_Spectra_e.k - k2_analytic_e) / k2_analytic_e * 100)  '%'])
disp(['h analytical: ' num2str(h2_analytic_e)])
disp(['h LOV3D: ' num2str(Love_Spectra_e.h)])
disp(['Normalized difference: ' num2str((Love_Spectra_e.h - h2_analytic_e) / h2_analytic_e * 100)  '%'])
disp([' % % % % % % % % % % % % % '])

% % % % % % % % % % % % % %
% Viscoelastic model      %
% % % % % % % % % % % % % %

% Core boundary layer (1)
Visco_elas_model(1).R0 = 1e-3 * R_io; % approx 0
Visco_elas_model(1).rho0 = rho_avg;

% Single body layer (2)
Visco_elas_model(2).R0 = R_io;
Visco_elas_model(2).rho0 = rho_avg;
Visco_elas_model(2).mu0 = mu_avg;
Visco_elas_model(2).eta0 = eta_avg;
Visco_elas_model(2).Ks0 = K_avg;

% forcing
Forcing_ve(1).Td = T_io;
Forcing_ve(1).n = 2; 
Forcing_ve(1).m = 0; 
Forcing_ve(1).F = 1;

[Numerics, Visco_elas_model] = set_boundary_indices(Numerics, Visco_elas_model,'empty');
Visco_elas_model = get_rheology(Visco_elas_model, Numerics, Forcing_ve);
[Love_Spectra_ve, y_ve] = get_Love(Visco_elas_model, Forcing_ve, Numerics,'empty');

% obtain the Fourier-transformed effective shear modulus
mu_eff_hat_ve = mu_eff * Visco_elas_model(2).muC;
% compute Love numbers using analytical expression
n_ve = Forcing_ve.n; 
mu_n_ve = ( 2 * n_ve ^ 2 + 4 * n_ve + 3 ) / n_ve * mu_eff_hat_ve;
k2_analytic_ve = 1 / ( 1 + mu_n_ve ) * 3 / 2 / ( n_ve - 1 );
h2_analytic_ve = 1 / ( 1 + mu_n_ve ) * ( 2 * n_ve + 1 ) / 2 / ( n_ve - 1 );

% display results
disp(['VISCOELASTIC CASE:'])
disp(['k analytical: ' num2str(k2_analytic_ve)])
disp(['k LOV3D: ' num2str(Love_Spectra_ve.k)])
disp(['Normalized difference: ' num2str((Love_Spectra_ve.k - k2_analytic_ve) / k2_analytic_ve * 100)  '%'])
disp(['h analytical: ' num2str(h2_analytic_ve)])
disp(['h LOV3D: ' num2str(Love_Spectra_ve.h)])
disp(['Normalized difference: ' num2str((Love_Spectra_ve.h - h2_analytic_ve) / h2_analytic_ve * 100)  '%'])
disp([' % % % % % % % % % % % % % '])

% % % % % % % % % % % % % %
% Ideal fluid model       %
% % % % % % % % % % % % % %

% Core boundary layer (1)
Ideal_fluid_model(1).R0 = 1e-3 * R_io; % approx 0
Ideal_fluid_model(1).rho0 = rho_avg;

% Single body layer (2)
Ideal_fluid_model(2).R0 = R_io;
Ideal_fluid_model(2).rho0 = rho_avg;
Ideal_fluid_model(2).mu0 = mu_avg;
Ideal_fluid_model(2).eta0 = eta_avg;
Ideal_fluid_model(2).Ks0 = K_avg;

% forcing
Forcing_if(1).Td = 1E12;
Forcing_if(1).n = 2; 
Forcing_if(1).m = 0; 
Forcing_if(1).F = 1;

[Numerics, Ideal_fluid_model] = set_boundary_indices(Numerics, Ideal_fluid_model,'empty');
Ideal_fluid_model = get_rheology(Ideal_fluid_model, Numerics, Forcing_if);
[Love_Spectra_if, y_if] = get_Love(Ideal_fluid_model, Forcing_if, Numerics,'empty');

% obtain the Fourier-transformed effective shear modulus
mu_eff_hat_if = mu_eff * Ideal_fluid_model(2).muC;
% compute Love numbers using analytical expression
n_if = Forcing_if.n; 
mu_n_if = ( 2 * n_if ^ 2 + 4 * n_if + 3 ) / n_if * mu_eff_hat_if;
k2_analytic_if = 1 / ( 1 + mu_n_if ) * 3 / 2 / ( n_if - 1 );
h2_analytic_if = 1 / ( 1 + mu_n_if ) * ( 2 * n_if + 1 ) / 2 / ( n_if - 1 );

% display results
disp(['IDEAL FLUID CASE:'])
disp(['k analytical: ' num2str(k2_analytic_if)])
disp(['k LOV3D: ' num2str(Love_Spectra_if.k)])
disp(['Normalized difference: ' num2str((Love_Spectra_if.k - k2_analytic_if) / k2_analytic_if * 100)  '%'])
disp(['h analytical: ' num2str(h2_analytic_if)])
disp(['h LOV3D: ' num2str(Love_Spectra_if.h)])
disp(['Normalized difference: ' num2str((Love_Spectra_if.h - h2_analytic_if) / h2_analytic_if * 100)  '%'])
disp([' % % % % % % % % % % % % % '])
