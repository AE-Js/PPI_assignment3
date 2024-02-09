%% SCRIPT USED TO TEST GET_LOVE
%close all
clear all
clc
addpath(genpath('/home/allard/Code/New_LOV3D/lov3d/'))
set(0,'defaulttextInterpreter','latex') 
mfile_name          = mfilename('fullpath');
[pathstr,name,ext] = fileparts(mfile_name);
cd(pathstr);
cd('..')
addpath(genpath(pwd))

% Location where any data will be saved, needs to end with /
save_location = 'files/output/multi_layer_example/';
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
mu_eff=5.2; % mu_mean/(gs*rho_mean*R) gs=gravity acceleration at the surface
R_Io = 1822; % Radius of Io in km
R_CMB = 965.66;%965; % Radius of Io's core

% Interior_Model goes from core (1) to surface (n) with indices
% Interior_Model values for the core
Interior_Model(1).R0 = R_CMB; % CMB in km 
Interior_Model(1).rho0 = 5157.96;%5150; % average density in kg m^-3
Interior_Model(1).rho0_2 = 5157.96;%5150; % rho_2 in kg m^-3 (used for icy moons)

% Interior_Model values for mantle/asthenosphere/crust layers
% Number of layers on top of the core should be equal to Numerics.Nlayers
% If the layer is elastic eta should be NaN
Interior_Model(2).R0 = 1200; % outer radius of this layer
Interior_Model(2).rho0 = 3244; % density in kg m^-3
Interior_Model(2).Ks0 = 200e9; % bulk modulus in Pa
Interior_Model(2).mu0 = 60e9; % shear modulus in Pa
Interior_Model(2).eta0 = 4.942e15; % viscosity in Pa s
Interior_Model(2).nR = 2; % degree of lateral variations
Interior_Model(2).mR = 0; % order of lateral variations
Interior_Model(2).variable_mu_p = 0; %peak-to-peak variations of shear modulus (in %)
Interior_Model(2).variable_eta_p = 50; %peak-to-peak variations of viscosity (in %)
Interior_Model(2).variable_K_p = 0; %peak-to-peak variations of bulk modulus (in %)

Interior_Model(3).R0 = 1500; % outer radius of this layer
Interior_Model(3).rho0 = 3244; % density in kg m^-3
Interior_Model(3).Ks0 = 200e9; % bulk modulus in Pa
Interior_Model(3).mu0 = 60e9; % shear modulus in Pa
Interior_Model(3).eta0 = 4.942e15; % viscosity in Pa s
Interior_Model(3).nR = 2; % degree of lateral variations
Interior_Model(3).mR = 0; % order of lateral variations
Interior_Model(3).variable_mu_p = 0; %peak-to-peak variations of shear modulus (in %)
Interior_Model(3).variable_eta_p = 50; %peak-to-peak variations of viscosity (in %)
Interior_Model(3).variable_K_p = 0; %peak-to-peak variations of bulk modulus (in %)

Interior_Model(4).R0 = R_Io; % outer radius of this layer
Interior_Model(4).rho0 = 3244; % density in kg m^-3
Interior_Model(4).Ks0 = 200e9; % bulk modulus in Pa
Interior_Model(4).mu0 = 60e9; % shear modulus in Pa
Interior_Model(4).eta0 = 4.942e15; % viscosity in Pa s
Interior_Model(4).nR = 2; % degree of lateral variations
Interior_Model(4).mR = 0; % order of lateral variations
Interior_Model(4).variable_mu_p = 0; %peak-to-peak variations of shear modulus (in %)
Interior_Model(4).variable_eta_p = 50; %peak-to-peak variations of viscosity (in %)
Interior_Model(4).variable_K_p = 0; %peak-to-peak variations of bulk modulus (in %)

% Set the rho difference such that the core can also represent an ocean
Interior_Model(1).Delta_rho0 = Interior_Model(1).rho0_2 - Interior_Model(2).rho0;

%% NUMERICS
Numerics.Nlayers = 4; % number of concentric layers. Including the core!
Numerics.method = 'combination'; % method of setting the radial points per layer: combination, variable, fixed and manual
Numerics.Nrbase = 500; % depending on the method this will determine the number of points per layer
Numerics.perturbation_order = 2; %maximum order to which couplings are considered
Numerics.solution_cutoff = 12; % maximum degree of solution, not used if perturbation order is specified
Numerics.load_couplings = 1; % 0=no loading, 1=loading of specific file, 2=searches for big enough file
Numerics.Nenergy = 12; % maximum degree to which energy dissipation is expanded 
Numerics.rheology_cutoff = 2; % maximum order of difference (so in log) up to which rheology is still used 
Numerics.parallel_sol = 0; % Calculate the solution using a parfor-loop either 0 or 1
Numerics.parallel_gen = 1; % Calculate potential coupling files using parfor-loops either 0 or 1
Numerics.coupling_file_location = 'files/couplings/'; % MUST END WITH A '/' , Location where coupling files are saved or should be saved 

[Numerics, Interior_Model] = set_boundary_indices(Numerics, Interior_Model,'verbose');

%% Small but necessary additions/checks of the code
% Make sure that Numerics.Nlayers is equal to length of Interior_Model
if ~(Numerics.Nlayers == length(Interior_Model))
    error('Error. \nSize of Interior_Model struct must be equal to Numerics.Nlayers')
end

for ilayer=2:Numerics.Nlayers
    if ~(size(Interior_Model(ilayer).variable_mu_p) == size(Interior_Model(ilayer).nR) & ...
            size(Interior_Model(ilayer).variable_mu_p) == size(Interior_Model(ilayer).mR))
        error(['Error. \nSize of variable_mu_p in layer #' num2str(ilayer) ' is not equal to the number of modes'])
    end
    if ~(size(Interior_Model(ilayer).variable_eta_p) == size(Interior_Model(ilayer).nR) & ...
            size(Interior_Model(ilayer).variable_eta_p) == size(Interior_Model(ilayer).mR))
        error(['Error. \nSize of variable_eta_p in layer #' num2str(ilayer) ' is not equal to the number of modes'])
    end
    if ~(size(Interior_Model(ilayer).variable_K_p) == size(Interior_Model(ilayer).nR) & ...
            size(Interior_Model(ilayer).variable_K_p) == size(Interior_Model(ilayer).mR))
        error(['Error. \nSize of variable_K_p in layer #' num2str(ilayer) ' is not equal to the number of modes'])
    end
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
Interior_Model_Uni=Interior_Model;
for ilayer=2:Numerics.Nlayers
    Interior_Model_Uni(ilayer).variable_mu_p=0; 
    Interior_Model_Uni(ilayer).variable_eta_p=0; 
    Interior_Model_Uni(ilayer).variable_K_p=0; 
end
% Calculates the rheology that will be used. First computes the complex 
% spherical harmonics and then performs a fourier transform to obtain the 
% coefficients that will be used throughout the computations 
Interior_Model = get_rheology(Interior_Model,Numerics,Forcing);

% prepares a spherically-symmetric model for comparison   
Interior_Model_Uni = get_rheology(Interior_Model_Uni,Numerics,Forcing);

%% GET LOVE NUMBERS
%obtains the Love number spectra
if Numerics.parallel_sol==1
    parfor i=1:length(Forcing)
        % spherically-symmetric model 
        [Love_Spectra_Uni(i),y_Uni(i)]=get_Love(Interior_Model_Uni,Forcing(i),Numerics, ...
                                                'verbose','out_file','T1','save_location',save_location);
        % lateral variations
        if i==1 %plot rheology maps
            [Love_Spectra(i),y(i)]=get_Love(Interior_Model,Forcing(i),Numerics, ...
                                            'plot_rheology_spectra',1,'verbose','out_file','T1','save_location',save_location);
        else
            [Love_Spectra(i),y(i)]=get_Love(Interior_Model,Forcing(i),Numerics, ...
                                            'plot_rheology_spectra',0,'verbose','out_file','T1','save_location',save_location);
        end
    end
else
    for i=1:length(Forcing)
        % spherically-symmetric model 
        [Love_Spectra_Uni(i),y_Uni(i)]=get_Love(Interior_Model_Uni,Forcing(i),Numerics, ...
                                                'verbose','out_file','T1','save_location',save_location);
        % lateral variations
        if i==1 %plot rheology maps
            [Love_Spectra(i),y(i)]=get_Love(Interior_Model,Forcing(i),Numerics, ...
                                            'plot_rheology_spectra',1,'verbose','out_file','T1','save_location',save_location);
        else
            [Love_Spectra(i),y(i)]=get_Love(Interior_Model,Forcing(i),Numerics, ...
                                            'plot_rheology_spectra',0,'verbose','out_file','T1','save_location',save_location);
        end
    end
end
%% GET ENERGY SPECTRA, can be commentted if elastic
[Energy_Spectra]=get_energy(y,Numerics,Forcing,Interior_Model,'verbose',1,'out_file','T1', ...
                    'save_location',save_location,'calc_E_contribution',1);
[Energy_Spectra_Uni]=get_energy(y_Uni,Numerics,Forcing,Interior_Model_Uni,'verbose',1,'out_file','T1', ...
                        'save_location',save_location,'calc_E_contribution',1);
Energy_Spectra_Norm=Energy_Spectra_Uni; 
Energy_Spectra_Norm.energy_integral_v=Energy_Spectra_Norm.energy_integral_v/Energy_Spectra_Uni.energy_integral_v(1); 
Energy_Spectra_Diff=Energy_Spectra;
Energy_Spectra_Diff.energy_integral_v=(Energy_Spectra.energy_integral_v-Energy_Spectra_Uni.energy_integral_v)/Energy_Spectra_Uni.energy_integral_v(1);

disp(["Time spent on all calculations" num2str(toc(Startall)) " s"])
%% PLOT ENERGY SPECTRA, can be commented if elastic 
% uniform, total
plot_energy_map(Energy_Spectra_Norm,'type','total','projection','mollweide', ...
    'title','Uniform Normalized','save_plot',[save_location 'Figure 3']);
% variable, total
plot_energy_map(Energy_Spectra_Norm,'type','total','projection','mollweide', ...
    'title','Variable','save_plot',[save_location 'Figure 4']);
% difference + cut
plot_energy_map(Energy_Spectra_Diff,'type','difference','projection','cut', ...
    'title','Difference','limits',[-0.4 0.4],'label','$\Delta\dot{e}/\dot{e}_0^u$','save_plot',[save_location 'Figure 5']);
%% GET MAPS AND PLOT
% (2,0)
[y_LatLon] = get_map(y(1),Interior_Model(end));
[y_LatLonUni] = get_map(y_Uni(1),Interior_Model_Uni(end));
yDif=y_LatLon; 
yDif.y=y_LatLon.y-y_LatLonUni.y;
%plot_map(y_LatLonUni,Interior_Model_UniU,'field_name','all','radial_point',floor(Numerics.Nr*0.5),'plot_title','Uniform Model (2,0)') %floor(Numerics.Nr*0.1)
%plot_map(yDif,Interior_ModelU,'field_name','all','radial_point',floor(Numerics.Nr*0.5),'plot_title','Difference (2,0)') %floor(Numerics.Nr*0.1)
plot_map(y_LatLon,Interior_Model(end),'field_name','all','radial_point',floor(Numerics.Nr*0.5), ...
    'plot_title','Variable Model','save_plot',[save_location 'Figure 6']);
plot_map(y_LatLon,Interior_Model(end),'field_name','test_strain_stress','radial_point',floor(Numerics.Nr*0.5), ...
    'plot_title','(2,0)','save_plot',[save_location 'Figure 7']);%,'save_plot',plot_name1)
plot_map(y_LatLonUni,Interior_Model_Uni(end),'field_name','test_strain_stress','radial_point',floor(Numerics.Nr*0.5), ...
    'plot_title','(2,0) Uniform','save_plot',[save_location 'Figure 8']);
%% (2,2)
[y_LatLon] = get_map(y(2),Interior_Model(end));
[y_LatLonUni] = get_map(y_Uni(2),Interior_Model_Uni(end));
yDif=y_LatLon; 
yDif.y=y_LatLon.y-y_LatLonUni.y;
%plot_map(y_LatLonUni,Interior_Model_UniU,'field_name','all','radial_point',floor(Numerics.Nr*0.5),'plot_title','Uniform Model (2,0)') %floor(Numerics.Nr*0.1)
%plot_map(yDif,Interior_ModelU,'field_name','all','radial_point',floor(Numerics.Nr*0.5),'plot_title','Difference (2,2)') %floor(Numerics.Nr*0.1)
%plot_map(y_LatLon,Interior_ModelU,'field_name','all','radial_point',floor(Numerics.Nr*0.5),'plot_title','Variable Model')
plot_map(y_LatLon,Interior_Model(end),'field_name','test_strain_stress','radial_point',floor(Numerics.Nr*0.5), ...
    'plot_title','(2,2)','save_plot',[save_location 'Figure 9']);

