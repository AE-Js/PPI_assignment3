clear all;
close all;
clc

%% Set relevant structs
% Set the number of rheology modes, the forcing, the perturbation order and
% the number of energy degrees that are taken into account. 
Numerics.Nrheo_max = 10;
Forcing(1).n = 2;
Forcing(1).m = 0;
Forcing(2).n = 2;
Forcing(2).m = -2;
Forcing(3).n = 2;
Forcing(3).m = 2;
Numerics.perturbation_order = 2;
Numerics.Nenergy = 12;
Numerics.parallel_sol = 0;
Numerics.parallel_gen = 1;

% Generate the forcing term string
str_forc=[];
for i=1:length(Forcing)
    if i==length(Forcing)
        str_forc=[str_forc num2str(Forcing(i).n) '_' num2str(Forcing(i).m)];
    else
        str_forc=[str_forc num2str(Forcing(i).n) '_' num2str(Forcing(i).m) '__'];
    end
end

% Generate coupling file name
coupling_file_name=['files/couplings/E_struct__Nrheomax__' num2str(Numerics.Nrheo_max) '__forc__' ...
    str_forc '__N__' num2str(Numerics.Nenergy) '__per' num2str(Numerics.perturbation_order) '.mat'];
verbose=1;

%% compute couplings 
disp('compute couplings')
Couplings = get_energy_couplings_all(Numerics.perturbation_order,Numerics.Nrheo_max,Numerics.Nenergy,Forcing,Numerics,'verbose');

%% RECOMMENDED. Save file using '-v7.3' with compression.
save(coupling_file_name,'-struct','Couplings','-v7.3')

%% 	'-v7.3' No compression
% save(coupling_file_name,'-struct','Couplings','-v7.3','-nocompression')