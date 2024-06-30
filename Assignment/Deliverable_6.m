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

% Location where any data will be saved, needs to end with /
save_location = 'files/output/multi_layer_example/';
if ~isfolder(save_location)
    mkdir(save_location)
end
Startall = tic;

% Define parameters
K_avg = 200e9; % Average bulk modulus [Pa]
eta_avg = 1e23; % Average shear viscosity [Pa s]

% % % % % % % % % % % % % %
% Variations  % % % % % % %
% % % % % % % % % % % % % %

% Set up individual element variation arrays
num = 25;
eta_scales = logspace(-5, 5, num);
K_scales = logspace(-5, 5, num);
mu_scales = logspace(-5, 5, num);

k_eta_variation = zeros(1, num);
h_eta_variation = zeros(1, num);
k_K_i_variation = zeros(1, num);
h_K_i_variation = zeros(1, num);
k_mu_variation = zeros(1, num);
h_mu_variation = zeros(1, num);

k_eta_mu_variation = zeros(num);
h_eta_mu_variation = zeros(num);

% Perform individual parameter variations
for i = 1:num
    % Bulk modulus variation
    [k_i, h_i] = tidal_response(K_scales(i), 1.0, 1.0);
    k_K_i_variation(i) = k_i;
    h_K_i_variation(i) = h_i;

    % Shear viscosity variation
    [k_i, h_i] = tidal_response(1.0, eta_scales(i), 1.0);
    k_eta_variation(i) = k_i;
    h_eta_variation(i) = h_i;

    % Shear modulus variation
    [k_i, h_i] = tidal_response(1.0, 1.0, mu_scales(i));
    k_mu_variation(i) = k_i;
    h_mu_variation(i) = h_i;
end

% Perform variation of two elements simultaneously
for i = 1:num
    for j = 1:num
        [k_i, h_i] = tidal_response(1.0, eta_scales(j), mu_scales(i));
        k_eta_mu_variation(i, j) = k_i;
        h_eta_mu_variation(i, j) = h_i;
    end

end

% % % % % % % % % % % % % %
% Plotting  % % % % % % % %
% % % % % % % % % % % % % %

% Set smaller figure size
figure_width = 600; % Width in pixels
figure_height = 200; % Height in pixels

% Figure 1: Bulk Modulus Variation
figure(1);
set(gcf, 'Position', [100, 100, figure_width, figure_height]);

% Subplot 1: Real and Imaginary parts of k variation
subplot(1, 2, 1);
yyaxis left;
semilogx(K_scales, real(k_K_i_variation), 'b');
ylabel('$\Re(k)$', 'Interpreter', 'latex');
yyaxis right;
semilogx(K_scales, imag(k_K_i_variation), 'r');
ylabel('$\Im(k)$', 'Interpreter', 'latex');
xlabel('$K/K_0$', 'Interpreter', 'latex');
title('$k$ Variation', 'Interpreter', 'latex');
legend({'$\Re(k)$', '$\Im(k)$'}, 'Interpreter', 'latex');
grid on;

% Subplot 2: Real and Imaginary parts of h variation
subplot(1, 2, 2);
yyaxis left;
semilogx(K_scales, real(h_K_i_variation), 'b');
ylabel('$\Re(h)$', 'Interpreter', 'latex');
yyaxis right;
semilogx(K_scales, imag(h_K_i_variation), 'r');
ylabel('$\Im(h)$', 'Interpreter', 'latex');
xlabel('$K/K_0$', 'Interpreter', 'latex');
title('$h$ Variation', 'Interpreter', 'latex');
legend({'$\Re(h)$', '$\Im(h)$'}, 'Interpreter', 'latex');
grid on;

% Save figure
exportgraphics(gcf,'Figures/K_variation.png','Resolution',300);

% Figure 2: Shear Viscosity Variation
figure(2);
set(gcf, 'Position', [200, 200, figure_width, figure_height]);

% Subplot 1: Real and Imaginary parts of k variation
subplot(1, 2, 1);
yyaxis left;
semilogx(eta_scales, real(k_eta_variation), 'b');
ylabel('$\Re(k)$', 'Interpreter', 'latex');
yyaxis right;
semilogx(eta_scales, imag(k_eta_variation), 'r');
ylabel('$\Im(k)$', 'Interpreter', 'latex');
xlabel('$\eta/\eta_0$', 'Interpreter', 'latex');
title('$k$ Variation', 'Interpreter', 'latex');
legend({'$\Re(k)$', '$\Im(k)$'}, 'Interpreter', 'latex');
grid on;

% Subplot 2: Real and Imaginary parts of h variation
subplot(1, 2, 2);
yyaxis left;
semilogx(eta_scales, real(h_eta_variation), 'b');
ylabel('$\Re(h)$', 'Interpreter', 'latex');
yyaxis right;
semilogx(eta_scales, imag(h_eta_variation), 'r');
ylabel('$\Im(h)$', 'Interpreter', 'latex');
xlabel('$\eta/\eta_0$', 'Interpreter', 'latex');
title('$h$ Variation', 'Interpreter', 'latex');
legend({'$\Re(h)$', '$\Im(h)$'}, 'Interpreter', 'latex');
grid on;

% Save figure
exportgraphics(gcf,'Figures/eta_variation.png','Resolution',300);

% Figure 3: Shear Modulus Variation
figure(3);
set(gcf, 'Position', [300, 300, figure_width, figure_height]);

% Subplot 1: Real and Imaginary parts of k variation
subplot(1, 2, 1);
yyaxis left;
semilogx(mu_scales, real(k_mu_variation), 'b');
ylabel('$\Re(k)$', 'Interpreter', 'latex');
yyaxis right;
semilogx(mu_scales, imag(k_mu_variation), 'r');
ylabel('$\Im(k)$', 'Interpreter', 'latex');
xlabel('$\mu/\mu_0$', 'Interpreter', 'latex');
title('$k$ Variation', 'Interpreter', 'latex');
legend({'$\Re(k)$', '$\Im(k)$'}, 'Interpreter', 'latex');
grid on;

% Subplot 2: Real and Imaginary parts of h variation
subplot(1, 2, 2);
yyaxis left;
semilogx(mu_scales, real(h_mu_variation), 'b');
ylabel('$\Re(h)$', 'Interpreter', 'latex');
yyaxis right;
semilogx(mu_scales, imag(h_mu_variation), 'r');
ylabel('$\Im(h)$', 'Interpreter', 'latex');
xlabel('$\mu/\mu_0$', 'Interpreter', 'latex');
title('$h$ Variation', 'Interpreter', 'latex');
legend({'$\Re(h)$', '$\Im(h)$'}, 'Interpreter', 'latex');
grid on;

% Save figure
exportgraphics(gcf,'Figures/mu_variation.png','Resolution',300);

% Figure 1: Real part of k variation
figure(4);
set(gcf, 'Position', [100, 100, figure_width, figure_height]);

% Subplot 1: Real part of k variation
subplot(1, 2, 1);
pcolor(mu_scales, eta_scales, real(k_eta_mu_variation'));
shading interp
colorbar;
set(gca, 'YScale', 'log');
set(gca, 'XScale', 'log');
title('Real $k$ Variation', 'Interpreter', 'latex');
xlabel('$\mu/\mu_0$', 'Interpreter', 'latex');
ylabel('$\eta/\eta_0$', 'Interpreter', 'latex');

% Subplot 2: Imaginary part of k variation
subplot(1, 2, 2);
pcolor(mu_scales, eta_scales, imag(k_eta_mu_variation'));
shading interp
colorbar;
set(gca, 'YScale', 'log');
set(gca, 'XScale', 'log'); 
title('Imaginary $k$ Variation', 'Interpreter', 'latex');
xlabel('$\mu/\mu_0$', 'Interpreter', 'latex');
ylabel('$\eta/\eta_0$', 'Interpreter', 'latex');

% Save figure
exportgraphics(gcf,'Figures/k_variation_matrix.png','Resolution',300);

% Figure 2: Real part of h variation
figure(5);
set(gcf, 'Position', [200, 200, figure_width, figure_height]);

% Subplot 1: Real part of h variation
subplot(1, 2, 1);
pcolor(mu_scales, eta_scales, real(h_eta_mu_variation'));
shading interp
colorbar;
set(gca, 'YScale', 'log');
set(gca, 'XScale', 'log'); 
title('Real $h$ Variation', 'Interpreter', 'latex');
xlabel('$\mu/\mu_0$', 'Interpreter', 'latex');
ylabel('$\eta/\eta_0$', 'Interpreter', 'latex');

% Subplot 2: Imaginary part of h variation
subplot(1, 2, 2);
pcolor(mu_scales, eta_scales, imag(h_eta_mu_variation'));
shading interp
colorbar;
set(gca, 'YScale', 'log');
set(gca, 'XScale', 'log'); 
title('Imaginary $h$ Variation', 'Interpreter', 'latex');
xlabel('$\mu/\mu_0$', 'Interpreter', 'latex');
ylabel('$\eta/\eta_0$', 'Interpreter', 'latex');

% Save figure
exportgraphics(gcf,'Figures/h_variation_matrix.png','Resolution',300);

% define function for iterations
function [k_i, h_i] = tidal_response(scale_K, scale_eta, scale_mu)

    %% Define Io parameters

    r_cr = 955 + 836.6 + 30; % Crust radius in [km] 
    r_m = 955 + 836.6; % Mantle radius in [km] 
    r_c = 955; % Core radius in [km] 

    rho_cr = 3000; % Crust density [kg/m³]
    rho_m = 3263; % Mantle density [kg/m³]
    rho_c = 5165; % Core density [kg/m³]


    eta_cr = 1e23; % Crust shear viscosity [Pa s]
    eta_m = 1e19; % Mantle shear viscosity [Pa s]
    eta_c = 0.0; % Core shear viscosity [Pa s]

    mu_cr = 6.5e10; % Crust shear modulus [Pa s]
    mu_m = 6e10; % Mantle shear modulus [Pa s]
    mu_c = 0.0; % Core shear modulus [Pa s]

    K_cr = 200e9; % Crust bulk modulus [Pa]
    K_m = 200e9; % Mantle bulk modulus [Pa]
    K_c = 200e9; % Core bulk modulus [Pa]

    omega_io = 4.1086E-05; % orbital frequency in rad/s

    % Set up layers

    % Core
    Interior_Model(1).R0 = r_c;
    Interior_Model(1).rho0 = rho_c;
    Interior_Model(1).Ks0 = scale_K * K_c;
    Interior_Model(1).mu0 = scale_mu * mu_c;
    Interior_Model(1).eta0 = scale_eta * eta_c;

    % Mantle
    Interior_Model(2).R0 = r_m;
    Interior_Model(2).rho0 = rho_m;
    Interior_Model(2).Ks0 = scale_K * K_m;
    Interior_Model(2).mu0 = scale_mu * mu_m;
    Interior_Model(2).eta0 = scale_eta * eta_m;

    % Crust
    Interior_Model(3).R0 = r_cr;
    Interior_Model(3).rho0 = rho_cr;
    Interior_Model(3).Ks0 = scale_K * K_cr;
    Interior_Model(3).mu0 = scale_mu * mu_cr;
    Interior_Model(3).eta0 = scale_eta * eta_cr;

    %radial discretization
    Numerics.Nlayers = 3; % number of concentric layers. Including the core!
    Numerics.method = 'variable'; % method of setting the radial points per layer, here fixed number of layers
    Numerics.Nrbase = 200; % depending on the method this will determine the number of points per layer
    %code parallelization
    Numerics.parallel_sol = 0; % Use a parfor-loop to call get_Love, either 0 or 1
    Numerics.parallel_gen = 0; % Calculate potential coupling files and the propagation inside get_solution using parfor-loops, either 0 or 1
    % lateral variations
    Numerics.perturbation_order = 2; %maximum order to which couplings are considered
    Numerics.solution_cutoff = 12; % maximum degree of solution, not used if perturbation order is specified
    Numerics.load_couplings = 1; % 0=no loading, 1=loading of specific file, 2=searches for big enough file
    Numerics.Nenergy = 12; % maximum degree to which energy dissipation is expanded 
    Numerics.rheology_cutoff = 2; % maximum order of difference (so in log) up to which rheology is still used 

    [Numerics, Interior_Model] = set_boundary_indices(Numerics, Interior_Model,'verbose');
        
    % Check that Numerics.Nlayers is equal to length of Interior_Model

    if ~(Numerics.Nlayers == length(Interior_Model))
        error('Error. \nSize of Interior_Model struct must be equal to Numerics.Nlayers')
    end

    % forcing

    Forcing(1).Td=2*pi/omega_io;
    Forcing(1).n=2; 
    Forcing(1).m=0; 
    Forcing(1).F=1;

    % Tidal response

    Interior_Model= get_rheology(Interior_Model,Numerics,Forcing);
    [Love_Spectra, ~]=get_Love(Interior_Model,Forcing,Numerics,'verbose');

    k_i = Love_Spectra.k;
    h_i = Love_Spectra.h;

end