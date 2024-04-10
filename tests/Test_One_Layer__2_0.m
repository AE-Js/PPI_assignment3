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
save_location = 'files/output/test_one_layer_2_0/';
if ~isfolder(save_location)
    mkdir(save_location)
end
Startall = tic;

%% Io non-dimensional
% Values will be non-dimensionalised based on the outer layer such that all
% radii are fractions. 
% All values are taken from Steinke et al. (2020a)
% non-dimensional parameters for Io (see Section 2.2 and Table 1)
omega0=4.1086E-05; %Io's orbital frequency 
T=2*pi/omega0; %Io's orbital period
R_Io = 1822; % Radius of Io in km
R_CMB = 965.66; % Radius of Io's core

% Interior_Model goes from core (1) to surface (n) with indices
% Interior_Model values for the core
Interior_Model(1).R0 = R_CMB; % CMB in km 
Interior_Model(1).rho0 = 5157.96; % average density in kg m^-3
Interior_Model(1).rho0_2 = 5157.96; % rho_2 in kg m^-3 (used for icy moons)

% Interior_Model values for mantle/asthenosphere/crust layers
% Number of layers on top of the core should be equal to Numerics.Nlayers
% If the layer is elastic eta should be NaN
Interior_Model(2).R0 = R_Io; % outer radius of this layer
Interior_Model(2).rho0 = 3244; % density in kg m^-3
Interior_Model(2).Ks0 = 200e9; % bulk modulus in Pa
Interior_Model(2).mu0 = 60e9; % shear modulus in Pa
Interior_Model(2).eta0 = 4.942e15; % viscosity in Pa s
Interior_Model(2).eta_variable_p2p(:,1) = 2; % degree of lateral variations
Interior_Model(2).eta_variable_p2p(:,2) = 0; % order of lateral variations
Interior_Model(2).eta_variable_p2p(:,3) = 50; %peak-to-peak variations of viscosity (in %)

% Set the rho difference such that the core can also represent an ocean
Interior_Model(1).Delta_rho0 = Interior_Model(1).rho0_2 - Interior_Model(2).rho0; 

%% NUMERICS
Numerics.Nlayers = 2; % number of concentric layers. Including the core!
Numerics.method = 'fixed'; % method of setting the radial points per layer
Numerics.Nrbase = 500; % depending on the method this will determine the number of points per layer
Numerics.perturbation_order = 2; %maximum order to which couplings are considered
Numerics.solution_cutoff = 12; % maximum degree of solution, not used if perturbation order is specified
Numerics.load_couplings = 1; % 0=no loading, 1=loading of specific file, 2=searches for big enough file
Numerics.Nenergy = 12; % maximum degree to which energy dissipation is expanded 
Numerics.rheology_cutoff = 2; % maximum order of difference (so in log) up to which rheology is still used 
Numerics.parallel_sol = 0; % Calculate the solution using a parfor-loop either 0 or 1
Numerics.parallel_gen = 0; % Calculate potential coupling files using parfor-loops either 0 or 1
Numerics.coupling_file_location = 'files/couplings/'; % MUST END WITH A '/' , Location where coupling files are saved or should be saved 

[Numerics, Interior_Model] = set_boundary_indices(Numerics, Interior_Model,'verbose');

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
% Build a uniform model based on the values of the lateral varying one
Interior_Model_Uni = Interior_Model;
for ilayer=2:Numerics.Nlayers
    Interior_Model_Uni(ilayer).eta_variable_p2p(:,1) = 0; % degree of lateral variations
    Interior_Model_Uni(ilayer).eta_variable_p2p(:,2) = 0; % order of lateral variations
    Interior_Model_Uni(ilayer).eta_variable_p2p(:,3) = 0; % peak-to-peak variations of viscosity (in %) 
end

% Calculates the rheology that will be used. First computes the complex 
% spherical harmonics and then performs a fourier transform to obtain the 
% coefficients that will be used throughout the computations 
Interior_Model = get_rheology(Interior_Model,Numerics,Forcing,'calculate_G',false,'include_sine',0);

% prepares a spherically-symmetric model for comparison   
Interior_Model_Uni = get_rheology(Interior_Model_Uni,Numerics,Forcing,'calculate_G',false,'include_sine',0);

%% GET LOVE NUMBERS
%obtains the Love number spectra
if Numerics.parallel_sol==1
    parfor i=1:length(Forcing)
        % spherically-symmetric model 
        [Love_Spectra_Uni(i),y_Uni(i)]=get_Love(Interior_Model_Uni,Forcing(i),Numerics,'verbose');
        % lateral variations
        [Love_Spectra(i),y(i)]=get_Love(Interior_Model,Forcing(i),Numerics,'verbose');
    end
else
    for i=1:length(Forcing)
        % spherically-symmetric model 
        [Love_Spectra_Uni(i),y_Uni(i)]=get_Love(Interior_Model_Uni,Forcing(i),Numerics,'verbose');
        % lateral variations
        [Love_Spectra(i),y(i)]=get_Love(Interior_Model,Forcing(i),Numerics,'verbose');
    end
end
%% GET ENERGY SPECTRA, can be commentted if elastic
[Energy_Spectra]=get_energy(y,Numerics,Forcing,Interior_Model,'verbose',1,'calc_E_contribution',1);
[Energy_Spectra_Uni]=get_energy(y_Uni,Numerics,Forcing,Interior_Model_Uni,'verbose',1,'calc_E_contribution',1);

disp(["Time spent on all calculations" num2str(toc(Startall)) " s"])

%% Check consistency with previous results/reference files, Non-Uniform
old_e_c_matrix2 = load("data/tests/nr2_mr0_rheoper2_per2_nE12_mup0_kp0_etap50/Energy_contribution_matrix__2_0__4_0__forc__2_0__2_-2__2_2__N__12__2.mat");
new_e_c_matrix2 = Energy_Spectra.energy_contribution;

C = abs(old_e_c_matrix2.energy_contribution_s - new_e_c_matrix2);

Max_diff_energy_contribution = max(C(:));
disp(['Maximum difference in the energy contribution matrix: ' num2str(Max_diff_energy_contribution)])

old_e_s2 = load("data/tests/nr2_mr0_rheoper2_per2_nE12_mup0_kp0_etap50/Energy_s__2_0__4_0__forc__2_0__2_-2__2_2__N__12__2.mat");
new_e_s2 = Energy_Spectra.energy_integral;

C2 = abs(old_e_s2.energy_s - new_e_s2);

Max_diff_energy = max(C2(:));
disp(['Maximum difference in the energy matrix: ' num2str(Max_diff_energy)])

Mean_diff_energy = mean(C2(:),"omitnan");
disp(['Mean difference in the energy matrix: ' num2str(Mean_diff_energy)])

C2_frac_diff = abs(C2./new_e_s2);
Max_frac_diff_energy = max(C2_frac_diff(:));
disp(['Maximum fractional difference in the energy matrix: ' num2str(Max_frac_diff_energy)])

Mean_frac_diff_energy = mean(C2_frac_diff(:),"omitnan");
disp(['Mean fractional difference in the energy matrix: ' num2str(Mean_frac_diff_energy)])

y_ref = load("data/tests/nr2_mr0_rheoper2_per2_nE12_mup0_kp0_etap50/y__2_0__4_0_per2__forc__2_0_per2.mat");
y_test = y(1).y; 

y_diff = abs(y_test - y_ref.y_sol);
y_diff_max = max(y_diff(:));
disp(['Maximum difference in the solution matrix: ' num2str(y_diff_max)])

y_diff_mean = mean(y_diff(:),"omitnan");
disp(['Mean difference in the solution matrix: ' num2str(y_diff_mean)])

k2_ref = y_ref.y_sol(end,8,2);
k2_test = y_test(end,8,2);
k2_im_diff = -(imag(k2_ref) - imag(k2_test))/imag(k2_test) * 100;
k2_im_test_arr = imag(y_test(end,8,:));
k2_im_ref_arr = imag(y_ref.y_sol(end,8,:));
k2_real_test_arr = real(y_test(end,8,:));
k2_real_ref_arr = real(y_ref.y_sol(end,8,:));

k2_im_diff_arr = (imag(y_ref.y_sol(end,8,:)) - imag(y_test(end,8,:)))./imag(y_test(end,8,:));
k2_real_diff_arr = (real(y_ref.y_sol(end,8,:)) - real(y_test(end,8,:)))./real(y_test(end,8,:));
disp(['Reference k2_1: ' num2str(y_ref.y_sol(end,8,1)) ',    New k2_1: ' num2str(y_test(end,8,1))])

disp(['k2 reference: ' num2str(k2_ref) ',    k2 test: ' num2str(k2_test)])
disp(['Percentage difference in k2: ' num2str(k2_im_diff)])


%% Uniform solution
old_e_c_matrix_uni = load("data/tests/nr2_mr0_rheoper2_per2_nE12_mup0_kp0_etap50/Energy_contribution_matrix__0_0__forc__2_0__2_-2__2_2__N__12__2.mat");
new_e_c_matrix_uni = Energy_Spectra_Uni.energy_contribution;

C_uni = abs(old_e_c_matrix_uni.energy_contribution_s - new_e_c_matrix_uni);

Max_diff_energy_contribution_uni = max(C_uni(:));
A_uni_mean = mean(C_uni(:),"omitnan");
C_uni_diff = abs(C_uni./new_e_c_matrix_uni);
A_uni_frac_max = max(C_uni_diff(:));
A_uni_frac_mean = mean(C_uni_diff(:),"omitnan");

disp(' ')
disp('Differences between UNIFORM bodies: ')
disp(['UNIFORM: Maximum difference in the energy contribution matrix: ' num2str(Max_diff_energy_contribution_uni)])
disp(['UNIFORM: Mean difference in the energy contribution matrix: ' num2str(A_uni_mean)])
disp(['UNIFORM: Maximum fractional difference in the energy contribution matrix: ' num2str(A_uni_frac_max)])
disp(['UNIFORM: Mean fractional difference in the energy contribution matrix: ' num2str(A_uni_frac_mean)])

old_e_s_uni = load("data/tests/nr2_mr0_rheoper2_per2_nE12_mup0_kp0_etap50/Energy_s__0_0__forc__2_0__2_-2__2_2__N__12__2.mat");
new_e_s_uni = Energy_Spectra_Uni.energy_integral;
C1_uni = abs(old_e_s_uni.energy_s - new_e_s_uni);

A1_uni = max(C1_uni(:));
A1_uni_mean = mean(C1_uni(:),"omitnan");
C1_uni_diff = abs(C1_uni./new_e_s_uni);
A1_uni_frac_max = max(C1_uni_diff(:));
A1_uni_frac_mean = mean(C1_uni_diff(:),"omitnan");

Max_diff_energy_Uni = max(C1_uni(:));
disp(['UNIFORM: Maximum difference in the energy matrix: ' num2str(Max_diff_energy_Uni)])

Mean_diff_energy_uni = mean(C1_uni(:),"omitnan");
disp(['UNIFORM: Mean difference in the energy matrix: ' num2str(Mean_diff_energy_uni)])

C1_frac_diff_uni = abs(C1_uni./new_e_s_uni);
Max_frac_diff_energy_uni = max(C1_frac_diff_uni(:));
disp(['UNIFORM: Maximum fractional difference in the energy matrix: ' num2str(Max_frac_diff_energy_uni)])

Mean_frac_diff_energy_uni = mean(C1_frac_diff_uni(:),"omitnan");
disp(['UNIFORM: Mean fractional difference in the energy matrix: ' num2str(Mean_frac_diff_energy_uni)])

