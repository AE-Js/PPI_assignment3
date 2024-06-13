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

% % % % % % % % % % % % % %
% Viscoelastic model      %
% % % % % % % % % % % % % %

% Core boundary layer (1)
Visco_elas_model(1).R0 = 1e-3 * R_io; % approx 0
Visco_elas_model(1).rho0 = rho_avg;

% Single body layer (2)
Visco_elas_model(2).R0 = R_io;
Visco_elas_model(2).rho0 = rho_avg;
Visco_elas_model(2).Ks0 = K_avg;

% forcing
Forcing_ve(1).Td = T_io;
Forcing_ve(1).n = 2; 
Forcing_ve(1).m = 0; 
Forcing_ve(1).F = 1;

% Define range for shear modulus and viscosity
mu_values = linspace(0.1 * mu_avg, 10 * mu_avg, 10);
eta_values = linspace(1e-5 * eta_avg, 1e5 * eta_avg, 10);

% Initialize matrices for storing results
k2_array = zeros(length(mu_values), length(eta_values));

% Compute Love numbers for each combination of mu and eta
for i = 1:length(mu_values)
    for j = 1:length(eta_values)
        mu = mu_values(i);
        eta = eta_values(j);

        Visco_elas_model(2).mu0 = mu;
        Visco_elas_model(2).eta0 = eta;
        
        [Numerics, Visco_elas_model] = set_boundary_indices(Numerics, Visco_elas_model,'empty');
        Visco_elas_model = get_rheology(Visco_elas_model, Numerics, Forcing_ve);
        [Love_Spectra_ve, y_ve] = get_Love(Visco_elas_model, Forcing_ve, Numerics,'empty');

        k2_ve = real(Love_Spectra_ve.k);

        % Store the result
        k2_array(i, j) = k2_ve;
    end
end

% Plot the results
figure;
% pcolor(log10(mu_values / mu_avg), log10(eta_values / eta_avg), k2_analytic');
[C, h] = contour(mu_values / mu_avg, eta_values / eta_avg, k2_array');
clabel(C, h, 'FontSize', 12, 'Color', 'k');
colorbar;
set(gca, 'YScale', 'log')
xlabel('$$\mu / \mu_0$$', 'Interpreter', 'latex');
ylabel('$$\eta / \eta_0$$', 'Interpreter', 'latex');
title('$$k_2$$ Love Numbers', 'Interpreter', 'latex');