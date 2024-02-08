clear all;
close all;
clc

Numerics.Nrheo_max = 4;
Forcing.n = 2;
Forcing.m = 2;
Numerics.perturbation_order = 1;
Numerics.parallel_sol = 0;
Numerics.parallel_gen = 1;

coupling_file_name=['Files_Coupling/L__Nrheomax__' num2str(Numerics.Nrheo_max) '__forc__' ...
        num2str(Forcing.n) '_' num2str(Forcing.m) '_per' num2str(Numerics.perturbation_order) '.mat'];

% Couplings = get_couplings_all(Numerics.perturbation_order,Numerics.Nrheo_max,Forcing,'verbose');
verbose=1;
Nrheo_max = Numerics.Nrheo_max;
perturbation_order = Numerics.perturbation_order;

%% compute couplings 
disp('compute couplings')
Couplings=get_couplings_all(perturbation_order,Nrheo_max,Forcing,Numerics,'verbose');

save(coupling_file_name,'Couplings','-v7.3')