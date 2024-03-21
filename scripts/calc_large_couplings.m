clear all;
close all;
clc

%% Set relevant structs
% Set the number of rheology modes, the specific forcing and 
% the perturbation order and are taken into account. 
Numerics.Nrheo_max = 20;
Forcing.n = 2;
Forcing.m = 2;
Numerics.perturbation_order = 1;
Numerics.parallel_sol = 0;
Numerics.parallel_gen = 1;

coupling_file_name=['files/couplings/L_struct__Nrheomax__' num2str(Numerics.Nrheo_max) '__forc__' ...
        num2str(Forcing.n) '_' num2str(Forcing.m) '_per' num2str(Numerics.perturbation_order) '.mat'];

verbose=1;
Nrheo_max = Numerics.Nrheo_max;
perturbation_order = Numerics.perturbation_order;

%% compute couplings 
disp('compute couplings')
Couplings = get_couplings_all(perturbation_order,Nrheo_max,Forcing,Numerics,'verbose');

%% RECOMMENDED. Save file using '-v7.3' With compression.
save(coupling_file_name,'-struct','Couplings','-v7.3')

%% Save file using '-v7.3' No compression.
% save(coupling_file_name,'-struct','Couplings','-v7.3','-nocompression')