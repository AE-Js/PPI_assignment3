%% GET_ENERGY
% AUTHOR: M. Rovira-Navarro 
% USE: obtain the energy spectra 
%% INPUT    
% Numerics
    %Numerics.Nlayers: number of layers, simply length(Interior_Model)
    %Numerics.method: used for the radial discretization in set_boundary_indices. Four methods are possible 
    %   - variable: This method sets the number of radial points in each layer
    %       equal to Numerics.Nrbase. The variable name thus comes from the fact
    %       that the stepsize is variable. 
    %   - fixed: This method sets the total number of points equal to Nrbase.
    %       Same as for the variable method, the name fixed thus means that the
    %       stepsize is kept fixed. WARNING: This method changes the physical
    %       boundary locations to ensure that they occur at a integer number of
    %       points.
    %    - combination: This method is a combination of the fixed and variable
    %       method. Numerics.Nrbase serves as the base number of points per layer
    %       onto which an additional number of points is added based on the
    %       physical size of the layer. The total number of points will depend on
    %       the number of layers but is roughly given by (Nlayers+1)*Nrbase
    %   - manual: This method allows the user to set the number of points per 
    %       layer manually. In order for it to work an array needs to be
    %       provided as a varargin after putting 'manual' as a varargin as well. 
    %       The syntax of the array needs to be: [0, points in layer #1, points
    %       in layer #2, etc]
    %Numerics.Nrbase: used for the radial discretization  (see above) 
    %Numerics.Nr: number of radial points (computed inside set_boundary_indices)
    %Numerics.perturbation_order: maximum order of the perturbation. Default 2
    %Numerics.rheology_cutoff: determines which terms of the rheology are included (only relevant for viscoelastic). terms with log10(mu_n^m)-log10(mu_n^m(leading))>=-Numerics.rheology_cutoff are included. Default 0 (only leading terms)
    %Numerics.Nenergy: maximum degree to which energy dissipation is expanded. Defaulft 8. 
    %Numerics.load_couplings: 
         % (0) compute coupling coefficients
         % (1) load coupling coefficients from file that contains ALL coupling up to rheology variations of a degree higher than those specied, if such a file does not exist compute it   (default) 
         % (2) load coupling coefficintes from a file that contains specifically the coupling coefficients for the rheology specied, if it does not exisit, compute it. 
    % Numerics.parallel_sol  Calculate the solution using a parfor-loop either 0 or 1
    % Numerics.parallel_gen  Calculate potential coupling files using parfor-loops either 0 or 1
% ------------------
% y: vector containing radial functions
    % y.n: degree
    % y.m: order
    % y.y(radial_point,X,mode) 
        % y(radial_point,1,mode): r radial position
        % y(radial_point,2,mode): U radial displacement
        % y(radial_point,3,mode): V tangential displacement
        % y(radial_point,4,mode): R normal stress
        % y(radial_point,5,mode): S tangential stress
        % y(radial_point,6,mode): \phi gravitational potential
        % y(radial_point,7,mode): \dot\phi potential stress
        % y(radial_point,8,mode): W toroidal displacement
        % y(radial_point,9,mode): T toroidal stress        
        % y(radial_point,10,mode): u_{n,n-1}
        % y(radial_point,11,mode): u_{n,n}
        % y(radial_point,12,mode): u_{n,n+1}
        % y(radial_point,13,mode): \sigma_{n,n,0}       
        % y(radial_point,14,mode): \sigma_{n,n-2,2}
        % y(radial_point,15,mode): \sigma_{n,n-1,2}
        % y(radial_point,16,mode): \sigma_{n,n,2}
        % y(radial_point,17,mode): \sigma_{n,n+1,2}
        % y(radial_point,18,mode): \sigma_{n,n+2,2}      
        % y(radial_point,19,mode): \epsilon_{n,n,0}
        % y(radial_point,20,mode): \epsilon_{n,n-2,2}
        % y(radial_point,21,mode): \epsilon_{n,n-1,2}
        % y(radial_point,22,mode): \epsilon_{n,n,2}
        % y(radial_point,23,mode): \epsilon{n,n+1,2}
        % y(radial_point,24,mode): \epsilon_{n,n+2,2} 
% ------------------
% Forcing: Vector containing the forcing information
    % Forcing.Td: forcing period
    % Forcing.n: degree of the forcing 
    % Forcing.m: order of the forcing 
    % Forcing.F: amplitude of the component 
% ------------------
% Interior_Model: Structure containing the interior model information
    % Mean properties: 
        %Interior_Model(ilayer).R0: radius [km]
        %Interior_Model(ilayer).rho0: density [kg.m^{-3}]       
        %Interior_Model(ilayer).mu0: shear modulus [Pa]
        %Interior_Model(ilayer).Ks0: bulk modulus, if not given, the model is assumed incompressible, [Pa]
        %Interior_Model(ilayer).eta0: viscosity, if not give, the body is assumed elastic [Pa.s]
        %Interior_Model(ilayer).MaxTime: normalized Maxwell time, it can be given instead of the viscosity, else it is computed [-]
        %Interior_Model(ilayer).ocean: (0) no ocean, (1) ocean
        %Interior_Model(ilayer).Gg0: gravitational constant (should be equal for all layers) [Nm^2kg^{-2}]
    % Lateral variations 
    % Lateral variations can be provided in three different formats
    %%% (1) in complex spherical harmonics: 
        %Interio_Model(ilayer).mu_variable: shear modulus variations
            %mu_variable(:,1): degree of variation 
            %mu_variable(:,2): order of variation 
            %mu_variable(:,3): amplitude of the variation (mu_l^m/mu^0_0)
        %Interio_Model(ilayer).K_variable: bulk modulus variations 
            %K_variable(:,1): degree of variation 
            %K_variable(:,2): order of variation 
            %K_variable(:,3):  amplitude of bulk modulus variations (K_l^m/K^0_0)
        %Interio_Model(ilayer).eta_variable: viscosity
            %eta_variable(:,1): degree of variation 
            %eta_variable(:,2): order of variation 
            %eta_variable(:,3):  amplitude of viscosity variations (eta_l^m/eta^0_0)
    %%% (2) in peak to peak variation amplitude wrt to the mean value (in percent)
        %Interio_Model(ilayer).mu_variable_p2p: shear modulus variations
            %mu_variable_p2p(:,1): degree of variation 
            %mu_variable_p2p(:,2): order of variation 
            %mu_variable_p2p:,3): amplitude of the variation (mu_l^m/mu^0_0)
        %Interio_Model(ilayer).K_variable_p2p: bulk modulus variations 
            %K_variable_p2p(:,1): degree of variation 
            %K_variable_p2p(:,2): order of variation 
            %K_variable_p2p(:,3):  amplitude of bulk modulus variations (K_l^m/K^0_0)
        %Interio_Model(ilayer).eta_variable: viscosity
            %eta_variable_p2p(:,1): degree of variation 
            %eta_variable_p2p(:,2): order of variation 
            %eta_variable_p2p(:,3):  amplitude of viscosity variations (eta_l^m/eta^0_0)  
    %%% (3) as a map, in which case get_rheology will transform it to
            % a spherical harmonics expansion. The maps can be generated by generated SPH_LatLon in SPH_Tools and 
            % for efficient transformation should be evaluated in the lat lon grid specified there, 
            % which depends on the maximum degree of the expansion. 
            % IMPORTANT NOTICE: The map(s) should already be scaled by their relevant average value 
            % (Interior_Model(ilayer).eta0, Interior_Model(ilayer).mu0 or Interior_Model(ilayer).Ks0) as: 
                    % Interior_Model(ilayer).mu_latlon/Interior_Model(ilayer).mu0
                    % Interior_Model(ilayer).eta_latlon/Interior_Model(ilayer).eta0
                    % Interior_Model(ilayer).K_latlon/Interior_Model(ilayer).Ks0
            % This ensures that the maps that you input have a mean of 1, which is what
            % the code expects, and that the actual fields have an average equal
            % to the prescribed average (Interior_Model(ilayer).XXX0)
            % Each of the input fields has the components 
                    %.lon:  longitude [IN DEGREES]
                    %.lat:  latitude [IN DEGREES]
                    %.z:    pxq grid where p stands for latitude and q for
                    %.lmax: maximum harmonic degree
            % for both (1) and (2) the expansion in spherical harmonics should result in 
            % a real field (there should be both m>0 and m<0 components)
            % If this is not the case the and only >0 components are provided, the <0 component is computed
    % Assigned in get_rheology 
    %%% Mean properties: 
        % Interior_Model(ilayer).R: normalized radius (normalized with Interior_Model(end).R0)
        % Interior_Model(ilayer).rho: normalized density  (normalized with Interior_Model(end).rho0)     
        % Interior_Model(ilayer).mu: normalized shear modulus of the layer (normalized with Interior_Model(end).mu0)
        % Interior_Model(ilayer).Ks: normalized bulk modulus of layer (normalized with Interior_Model(end).mu0)
        % Interior_Model(ilayer).Ks: normalized viscosity of layer (normalized with Interior_Model(end).mu0*Forcing.T)
        % Interior_Model(ilayer).Gg: normalized gravitational constant (can also be given by user in model given in non-dimensional units)
        % Interior_Model(ilayer).gs: normalized gravity at Interior_Model(ilayer).R0
        % Interior_Model(ilayer).gs0: gravity at Interior_Model(ilayer).R0
        %Interior_Model(ilayer).muC: normalized complex shear modulus 
        % Interior_Model(ilayer).mu00R: normalized real component of the shear modulus 
         % Interior_Model(ilayer).rho0_av: averaged denseity at Interior_Model(ilayer).R0
        % Interior_Model(ilayer).rho_av: normalized averaged denseity at Interior_Model(ilayer).R0
     %%% Lateral variations:
        % Interior_Model(ilayer).rheology_variable: rheology variable 
        % rheology_variable(:,1): degree of variation 
        % rheology_variable(:,2): order of variation 
        % rheology_variable(:,3): amplitude of bulk modulus variations (K_l^m/K^0_0)  
        % rheology_variable(:,4): amplitude of complex shear modulus variations (\hat mu_l^m/mu^0_0(N))
%% OUPUT
    %Energy_Spectra: 
        %Energy_Spectra.n: degrees with non-zero energy
        %Energy_Spectra.m: orders with non-zero energy 
        %Energy_Spectra.n_v: degrees from 0 to Numerics.Nenergy 
        %Energy_Spectra.n_v: orders from 0 to Numerics.Nenergy 
        %Energy_Spectra.energy(radial_point,mode): radial profile of energy spectra
        %Energy_Spectra.energy_integral(mode): radially integrated energy for all non-zero degrees an orders (n,m)
        %Energy_Spectra.energy_integral_v(mode): radially integrated energy for all degrees an orders (n_v,m_v)
%% FUNCTION ------------------------------------------------------------------   
function [Energy_Spectra]=get_energy(y,Numerics,Forcing,Interior_Model,varargin)
%% (0) OPTIONAL INPUTS
verbose = 0;
out_file = 0; 
calc_E_contributions = 0;
save_energy_vec = 0;
if isfield(Numerics,'load_couplings')==0
    Numerics.load_couplings=1; 
end
if isfield(Numerics,'Nenergy')==0
    Numerics.Nenergy=8; 
end

for k = 1:length(varargin)
    if strcmpi(varargin{k},'verbose')
        verbose=varargin{k+1}; 
        varargin{k+1}=[]; 
        varargin{k}=[];
    end
    if strcmpi(varargin{k},'out_file')
        out_file=varargin{k+1}; 
        varargin{k+1}=[]; 
        varargin{k}=[];
    end
    if strcmpi(varargin{k},'save_location')
        save_location=varargin{k+1}; 
        varargin{k+1}=[]; 
        varargin{k}=[];
    end
    if strcmpi(varargin{k},'calc_E_contribution')
        calc_E_contributions=varargin{k+1}; 
        varargin{k+1}=[]; 
        varargin{k}=[];
    end
    if strcmpi(varargin{k},'save_energy_vec')
        save_energy_vec=varargin{k+1}; 
        varargin{k+1}=[]; 
        varargin{k}=[];
    end
end
%% (1) COMBINE SOLUTIONS & SET IN CORRECT FORMAT
% set solution in correct format 
Nf = length(y); %number of forcings
n_sol = [];
m_sol = [];
for i=1:Nf
    n_sol = [n_sol y(i).n];
    m_sol = [m_sol y(i).m];
end
n_aux = unique(sort(abs(n_sol)));
m_sol_t = [];
n_sol_t = [];
k = 1;
for i=1:length(n_aux)
    index = find(n_aux(i)==n_sol);
    m_aux = unique(abs(m_sol(index)));
    m_aux = sort(unique([m_aux -m_aux]));
    orders_l = length(m_aux);
    m_sol_t(k:k+orders_l-1) = m_aux;
    n_sol_t(k:k+orders_l-1) = n_aux(i);
    k = k + orders_l;
end

% Build the total solution 
y_total = zeros(Numerics.Nr+1,24,length(n_sol_t));
for i=1:length(n_sol_t)
    for j=1:Nf
        ind = find(n_sol_t(i)==y(j).n & m_sol_t(i)==y(j).m);
        y_total(:,1,i) = y(1).y(:,1,1); 
        if isempty(ind)==0
            y_total(:,2:end,i) = y_total(:,2:end,i) + Forcing(j).F*y(j).y(:,2:end,ind); 
        end
    end
end
%% strain and stress tensor and the - component 
stressP = zeros(Numerics.Nr+1,6,length(n_sol_t));
strainP = zeros(Numerics.Nr+1,6,length(n_sol_t));
stressN = zeros(Numerics.Nr+1,6,length(n_sol_t));
strainN = zeros(Numerics.Nr+1,6,length(n_sol_t));

% Fill the Stress and Strain tensors
for i=1:length(n_sol_t)
    index_auxN = find(n_sol_t==n_sol_t(i) & m_sol_t==-m_sol_t(i));
    stressP(:,:,i) = y_total(:,[14 15 16 17 18 13],i);
    strainP(:,:,i) = y_total(:,[20 21 22 23 24 19],i);
    stressN(:,:,i) = conj(y_total(:,[14 15 16 17 18 13],index_auxN)); 
    strainN(:,:,i) = conj(y_total(:,[20 21 22 23 24 19],index_auxN));
end

% Retrieve the radial points
r = y(1).y(:,1,1);

%% (2) GET COUPLINGS ENERGY SPECTRA
% obtain name of the coupling file
str_forc = [];

% Check whether there are any uniform layers
uniformlayers = [];
for ilayer=2:Numerics.Nlayers
    uniformlayers = [uniformlayers Interior_Model(ilayer).uniform];
end
uniformlayers = uniformlayers == 0;

% Generate the forcing term string
for i=1:length(Forcing)
    if i==length(Forcing)
        str_forc = [str_forc num2str(Forcing(i).n) '_' num2str(Forcing(i).m)];
    else
        str_forc = [str_forc num2str(Forcing(i).n) '_' num2str(Forcing(i).m) '__'];
    end
end

% Obtain the highest degree of rheology across the layers that is used
Numerics.Nrheo_max = 0;
for ilayer=2:Numerics.Nlayers
    prevmax = max(Interior_Model(ilayer).rheology_variable(:,1));
    Numerics.Nrheo_max = max([Numerics.Nrheo_max prevmax]);
end

% If there is only one non-uniform layer obtain couplings from large file
if any(uniformlayers)   
    % Decide whether to load from a file containing all the modes for this specific rheology and settings 
    % or coupling coefficients from a large generic file

    if Numerics.load_couplings == 1 % Load from specific file
        % Remove rheology terms with 0 amplitude 
        % Fill an array with all the used degree and orders (contains duplicates)
        for ilayer=2:Numerics.Nlayers
            non_zero_rheo = find(abs(Interior_Model(ilayer).rheology_variable(:,4)) > 0);
            Interior_Model(ilayer).rheology_variable = Interior_Model(ilayer).rheology_variable(non_zero_rheo,:);
            if ilayer == 2
                rheo_degree_orders = Interior_Model(ilayer).rheology_variable(:,1:2);
            else
                rheo_degree_orders = [rheo_degree_orders ; Interior_Model(ilayer).rheology_variable(:,1:2)];
            end
        end
    
        % Remove duplicates from rheo_degree_orders
        rheo_degree_orders = unique(rheo_degree_orders,'rows');
        
        % Fill character array with all unique rheology modes
        str_rheo = [];
        for i=1:size(rheo_degree_orders,1)
            str_rheo = [str_rheo num2str(rheo_degree_orders(i,1)) '_' num2str(rheo_degree_orders(i,2)) '__'];
        end

        % Generate coupling file name
        coupling_file_name = [Numerics.coupling_file_location 'E_rheo__' str_rheo 'forc__' str_forc '__N__' ... 
                              num2str(Numerics.Nenergy) '__per' num2str(Numerics.perturbation_order) '.mat'];
        
        % Check if the file has already been computed
        if isfile(coupling_file_name)==1 % load couplings
            if verbose==1
                disp(['Coupling Loaded from: ' coupling_file_name])
                tic
            end
            
            % Load couplings file
            Couplings = load(coupling_file_name);

            if verbose==1
                disp(['Loading coupling file took: ' num2str(toc) ' seconds'])
            end
        else % compute couplings
            if verbose==1
                tic
                disp(['File ' coupling_file_name ' not found. Computing all coupling coefficients, this might take some time..'])
                Couplings = get_energy_couplings(n_sol_t,m_sol_t,Numerics,'verbose');
                disp(['Time Spent: ' num2str(toc) 's'])
                disp(['File stored in: ' coupling_file_name])
            else
                Couplings = get_energy_couplings(n_sol_t,m_sol_t,Numerics,'verbose');
            end
            
            % Save couplings for all terms
            save(coupling_file_name,'-struct','Couplings','-v7.3')
        end

    elseif Numerics.load_couplings == 2 % Load from bigger/generic file
        % Generate coupling file name
        coupling_file_name = [Numerics.coupling_file_location 'E_struct__Nrheomax__' num2str(Numerics.Nrheo_max) '__forc__' ...
            str_forc '__N__' num2str(Numerics.Nenergy) '__per' num2str(Numerics.perturbation_order) '.mat'];
        
        % look if there is a file that contains the couplings
        coupling_file_name_search = [Numerics.coupling_file_location 'E_struct__Nrheomax__*__forc__' str_forc '__N__*__per*.mat'];
        possible_couplings_files = dir(coupling_file_name_search);
        file_found = 0;
        i = 1; 
        potential_coupling_files = [];
        rheo_max_file_list = [];
    
        % Select the smallest possible file that contains all the couplings
        while i<=length(possible_couplings_files)
            per_start = strfind(possible_couplings_files(i).name,'_per')+4;
            per_end = strfind(possible_couplings_files(i).name,'.mat')-1;
            perturbation_order_file = str2num(possible_couplings_files(i).name(per_start:per_end));
            rheo_start = strfind(possible_couplings_files(i).name,'Nrheomax__')+10;
            rheo_end = strfind(possible_couplings_files(i).name,'__forc')-1;
            Nrheomax_file = str2num(possible_couplings_files(i).name(rheo_start:rheo_end));
            Nenergy_start = strfind(possible_couplings_files(i).name,'__N__')+5;
            Nenergy_end = strfind(possible_couplings_files(i).name,'__per')-1;
            Nenergy_file = str2num(possible_couplings_files(i).name(Nenergy_start:Nenergy_end));
    
            % Check whether the file fullfills the criteria and add it to the list
            if Nrheomax_file>=Numerics.Nrheo_max && perturbation_order_file>=Numerics.perturbation_order && Nenergy_file >= Numerics.Nenergy
                interim_name = [Numerics.coupling_file_location possible_couplings_files(i).name];
                potential_coupling_files = [potential_coupling_files convertCharsToStrings(interim_name)];
                rheo_max_file_list = [rheo_max_file_list Nrheomax_file];
            end
            i=i+1;
        end
        
        % If suitable files have been found select the smallest one
        if ~isempty(potential_coupling_files)
            [~,ind] = min(rheo_max_file_list);
            coupling_file_name = potential_coupling_files(ind);
            file_found=1; 
        end
    
        % Coupling file exist 
        if file_found==1 
            if verbose==1
                disp([' Energy Couplings Loaded from: ' convertStringsToChars(coupling_file_name)])
                tic
            end
    
            % Retrieve the coupling coefficients that are required
            Couplings = retrieve_energy_couplings_from_file(n_sol_t,m_sol_t,Numerics.Nenergy,coupling_file_name);
    
            if verbose==1
                disp(['Loading coupling file took: ' num2str(toc) ' seconds'])
            end
        else % Compute couplings 
            if verbose==1
                tic
                disp(['File ' coupling_file_name ' not found. Computing all coupling coefficients, this might take some time..'])
                Couplings = get_energy_couplings_all(Forcing,Numerics,'verbose');
                disp(['Time Spent: ' num2str(toc) 's'])
                disp(['File stored in: ' coupling_file_name])
            else
                Couplings = get_energy_couplings_all(Forcing,Numerics);
            end
    
            % Save the file
            save(coupling_file_name,'-struct','Couplings','-v7.3')
            
            % Retrieve the coupling coefficients that are required
            Couplings = retrieve_energy_couplings(n_sol_t,m_sol_t,Numerics.Nenergy,Couplings);
        end
    end

    % Extract relevant information from struct
    EC = Couplings.EC;
    n_en = Couplings.n_en;
    m_en = Couplings.m_en;

    % Double-check whether the solution modes used in the coupling matrix are the
    % same as in the actual solution
    if any(~(n_sol_t == Couplings.n_s)) && any(~(m_sol_t == Couplings.m_s))
        error('Error. \nModes used in the coupling matrix are not the same as in the current solution')
    end

    % Delete the potentially large struct from memory
    clear Couplings

else
    if verbose==1
        disp('Only uniform layers so coupling matrix is simple')
    end

    % Set coupling file name
    if isfield(Numerics,'coupling_file_location')==0
        Numerics.coupling_file_location='data/couplings/';
    end
    coupling_file_name = [Numerics.coupling_file_location 'E__rheo__0_0__forc__' str_forc '__N__' ...
                          num2str(Numerics.Nenergy) '__per' num2str(Numerics.perturbation_order) '.mat'];
    % Check whether the file exists
    if isfile(coupling_file_name)==1
        % coupling file exist and does not need to be computed
        if verbose==1
            disp([' Energy couplings loaded from: ' coupling_file_name])
        end
        Couplings = load(coupling_file_name);
        
        % Extract relevant information from struct
        EC = Couplings.EC;
        n_en = Couplings.n_en;
        m_en = Couplings.m_en;

    else
        if verbose==1
            disp(' Computing energy couplings...')
            tic
            Couplings = get_energy_couplings(n_sol_t,m_sol_t,Numerics,'verbose');
        else
            Couplings = get_energy_couplings(n_sol_t,m_sol_t,Numerics);
        end
         
        % Save the file
        save(coupling_file_name,'-struct','Couplings','-v7.3')

        % Extract relevant information from struct
        EC = Couplings.EC;
        n_en = Couplings.n_en;
        m_en = Couplings.m_en;

        if verbose==1
            disp(['Time Spent: ' num2str(toc) 's'])
            disp(['File stored in: ' coupling_file_name])
        end
    end
end

if out_file~=0
    out_file_name = [save_location out_file 'E__Nrheomax__' num2str(Numerics.Nrheo_max) '__forc__' ...
                    str_forc '__N__' num2str(Numerics.Nenergy) '__per' num2str(Numerics.perturbation_order) '.txt'];
    FID=fopen(out_file_name,'w');
    out_file_name_E_contribution = [save_location 'Energy_contribution_matrix__Nrheomax__' num2str(Numerics.Nrheo_max) ...
                                   '__forc__' str_forc '__N__' num2str(Numerics.Nenergy) '__per' ...
                                   num2str(Numerics.perturbation_order) '.mat'];
    out_file_name_E_s = [save_location 'Energy_s__Nrheomax__' num2str(Numerics.Nrheo_max) '__forc__' str_forc ...
                        '__N__' num2str(Numerics.Nenergy) '__per' num2str(Numerics.perturbation_order) '.mat'];
end

%% (3) COMPUTE ENERGY SPECTRA
if verbose==1
    disp('Computing energy spectra...')
    tic
end
% [energy,energy_s] = get_energy_matrix_mex(r,n_en,n_sol_t,m_sol_t,EC,stressN,stressP,strainN,strainP);

energy = zeros(length(r),length(n_en));
energy_s = zeros(1,length(n_en)); 

% variable to store all the contributions 
if calc_E_contributions == 1
    energy_contribution = zeros(length(r),length(n_en),length(n_sol_t),length(n_sol_t));
    energy_contribution_s = zeros(1,length(n_en),length(n_sol_t),length(n_sol_t));
end

% Radial difference between all radial points
radialdiff = ((r(1:end-1)+r(2:end))/2).^2 .* ((r(2:end)-r(1:end-1)));

% Find the indices where the coupling term is non-zero
[nonz_sol_mode1,nonz_sol_mode2,nonz_en,nonz_n2a,nonz_n2b] = ind2sub(size(EC),find(EC));
sol_eq = zeros(length(r),length(nonz_en));

% Calculate the energy for every combination of modes based on the non-zero
% energy couplings
for i=1:length(nonz_en)
    k = nonz_en(i);
    i1 = nonz_sol_mode1(i);
    i2 = nonz_sol_mode2(i);
    i3 = nonz_n2a(i);
    i4 = nonz_n2b(i);
    n2a = [n_sol_t(i1)-2:1:n_sol_t(i1)+2 n_sol_t(i1)];
    n2b = [n_sol_t(i2)-2:1:n_sol_t(i2)+2 n_sol_t(i2)];
    term1 = 1i*2*pi*(-1)^(n2a(i3)+n_sol_t(i1)-m_sol_t(i1))*stressN(:,i3,i1).*strainP(:,i4,i2)*EC(i1,i2,k,i3,i4);
    term2 = 1i*2*pi*(-1)^(n2b(i4)+n_sol_t(i2)-m_sol_t(i2))*stressP(:,i3,i1).*strainN(:,i4,i2)*EC(i1,i2,k,i3,i4);
    sol_eq(:,i) = term1-term2;
end

% Delete some large arrays from memory
clear nonz_n2a nonz_n2b

% Fill the energy and (potentially) the energy_contribution arrays
for n=1:length(n_en)
    ind_emodes = find(nonz_en == n);
    emodes_mask = nonz_en == n;
    if ~isempty(ind_emodes)
        energy(:,n) = sum(sol_eq(:,ind_emodes),2);
        if calc_E_contributions == 1
            for i1=1:length(n_sol_t)
                for i2=1:length(n_sol_t)
                    ind_n_sol_mode1 = nonz_sol_mode1 == i1;
                    ind_n_sol_mode2 = nonz_sol_mode2 == i2;
                    ind = emodes_mask & ind_n_sol_mode1 & ind_n_sol_mode2;
                    if any(ind)
                        energy_contribution(:,n,i1,i2) = sum(sol_eq(:,ind),2);
                    end
                end
            end
        end
    end

    % Radial integration 
    energy_s(n) = sum(radialdiff.*(energy(2:end,n) + energy(1:end-1,n))/2, 1);
    % energy_s(n) = trapz(r, energy(:,n) .* r.^2);
    if calc_E_contributions == 1
        energy_contribution_s(1,n,:,:) = sum(radialdiff.*(energy_contribution(2:end,n,:,:)+energy_contribution(1:end-1,n,:,:))/2, 1);
    end
end

if verbose==1
    disp(['Time Spent: ' num2str(toc) 's'])
end

%% (4) SET IN OUTPUT FORMAT 
Energy_Spectra.n = n_en;
Energy_Spectra.m = m_en;
Energy_Spectra.energy_integral = energy_s;
Energy_Spectra.energy = energy;
if calc_E_contributions == 1
    Energy_Spectra.energy_contribution = energy_contribution_s;
end
% get also the complete spectra
i=1; 
for n=0:Numerics.Nenergy
    for m=-n:n
        Energy_Spectra.n_v(i)=n; 
        Energy_Spectra.m_v(i)=m; 
        ind=find(n==Energy_Spectra.n & m==Energy_Spectra.m);
        if isempty(ind)==0
            Energy_Spectra.energy_integral_v(i)=Energy_Spectra.energy_integral(ind);
            Energy_Spectra.energy_v(:,i) = Energy_Spectra.energy(:,ind);
        else
            Energy_Spectra.energy_integral_v(i)=0;
            Energy_Spectra.energy_v(:,i) = zeros(length(r),1);
        end
        i=i+1; 
    end
end
%% PRINT SOME INFORMATION 
if out_file~=0
    fprintf(FID, '%s\n','------------- INTERIOR MODEL -----------');
    fprintf(FID, '%s\n', 'AVERAGE PROPERTIES (dimensional)');
    fprintf(FID, '%s\n', ['Layer#     R[m]    rho[kg.m^{-3}]    mu[Pa]    K[Pa]    eta[Pa.s]' ]);
    fprintf(FID, '%s\n', ['1    ' num2str(Interior_Model(1).R0,'%10.5e') '    ' num2str(Interior_Model(1).rho0,'%10.5e') '    0'  '    -'   ]);
    for ilayer=2:Numerics.Nlayers
        fprintf(FID, '%s\n', [num2str(ilayer) '    ' num2str(Interior_Model(ilayer).R0,'%10.5e') '    ' num2str(Interior_Model(ilayer).rho0,'%10.5e') '    ' num2str(Interior_Model(ilayer).mu0,'%10.5e')  '    ' num2str(Interior_Model(ilayer).Ks0,'%10.5e') ' ' num2str(Interior_Model(ilayer).eta0,'%10.5e')   ]);
    end
    fprintf(FID, '%s\n', 'AVERAGE PROPERTIES (non-dimensional)');
    fprintf(FID, '%s\n', ['Layer#     R[-]    rho[-]    mu[-]    K [-]     eta[-]' ]);
    fprintf(FID, '%s\n', ['1    ' num2str(Interior_Model(1).R) '    ' num2str(Interior_Model(1).rho) '    0'  '    -'   ]);
    for ilayer=2:Numerics.Nlayers
        fprintf(FID, '%s\n', [num2str(ilayer) '    ' num2str(Interior_Model(ilayer).R) '    ' num2str(Interior_Model(ilayer).rho) '    ' num2str(Interior_Model(ilayer).mu)  '    ' num2str(Interior_Model(ilayer).Ks) ' ' num2str(Interior_Model(ilayer).eta)   ]);
    end
    fprintf(FID, '%s\n', 'RHEOLOGY  VARIATIONS');
    fprintf(FID, '%s\n', ['Rheology cutoff =        ' num2str(Numerics.rheology_cutoff)]);
    fprintf(FID, '%s\n', 'Shear Modulus');
    for ilayer=2:Numerics.Nlayers
        if isfield(Interior_Model,'mu_variable')==1
            fprintf(FID, '%s\n', 'Layer#      (n,m)           amplitude[mu_n^m/mu_0^0]');
            for i=1:size(Interior_Model(ilayer).mu_variable,1)
                fprintf(FID, '%s\n', [num2str(ilayer) '    (' num2str(Interior_Model(ilayer).mu_variable(i,1)) ',' num2str(Interior_Model(ilayer).mu_variable(i,2)) ')   '  num2str(Interior_Model(ilayer).mu_variable(i,3),'%10.5e') ]);
            end
        else
            fprintf(FID, '%s\n', 'None');
        end 
    end
    
    fprintf(FID, '%s\n', 'Bulk Modulus');
    for ilayer=2:Numerics.Nlayers
        if isfield(Interior_Model,'K_variable')==1
            fprintf(FID, '%s\n', 'Layer#      (n,m)           amplitude[K_n^m/K_0^0]');
            for i=1:size(Interior_Model(ilayer).K_variable,1)
                fprintf(FID, '%s\n', [num2str(ilayer) '    (' num2str(Interior_Model(ilayer).K_variable(i,1)) ',' num2str(Interior_Model(ilayer).K_variable(i,2)) ')   '  num2str(Interior_Model(ilayer).K_variable(i,3),'%10.5e') ]);
            end
        else
            fprintf(FID, '%s\n', 'None')
        end
    end

    fprintf(FID, '%s\n', 'Viscosity');
    for ilayer=2:Numerics.Nlayers
        if isfield(Interior_Model,'eta_variable')==1
            fprintf(FID, '%s\n', 'Layer#      (n,m)          amplitude[eta_n^m/K_0^0]');
            for i=1:size(Interior_Model(ilayer).eta_variable,1)
                fprintf(FID, '%s\n', [num2str(ilayer) '    (' num2str(Interior_Model(ilayer).eta_variable(i,1)) ',' num2str(Interior_Model(ilayer).eta_variable(i,2)) ')   '  num2str(Interior_Model(ilayer).eta_variable(i,3),'%10.5e') ]);
            end
        else
            fprintf(FID, '%s\n', 'None');
        end
    end
    
    fprintf(FID, '%s\n', 'Complex Shear Modulus');
    for ilayer=2:Numerics.Nlayers
        if isfield(Interior_Model,'rheology_variable')==1
            fprintf(FID, '%s\n', 'Layer#      (n,m)         amplitude[\hat\mu_n^m/K_0^0]');
            for i=1:size(Interior_Model(ilayer).rheology_variable,1)
                fprintf(FID, '%s\n', [num2str(ilayer) '    (' num2str(Interior_Model(ilayer).rheology_variable(i,1)) ',' num2str(Interior_Model(ilayer).rheology_variable(i,2)) ')    '  num2str(Interior_Model(ilayer).rheology_variable(i,4),'%10.5e') ]);
            end
        else
            fprintf(FID, '%s\n', 'None');
        end
    end
    fprintf(FID, '%s\n', ' ');
    fprintf(FID, '%s\n', '----------- FORCING ----------');
    fprintf(FID, '%s\n', 'TIDAL POTENTIAL');
    fprintf(FID, '%s\n', 'Period [s]');
    fprintf(FID, '%s\n', num2str(Forcing(1).Td,'%10.5e'));
    fprintf(FID, '%s\n', '(n,m)      Amp');
    for i=1:length(Forcing)
        fprintf(FID, '%s\n', ['(' num2str(Forcing(i).n) ',' num2str(Forcing(i).m) ')       '    num2str(Forcing(i).F)]);
    end 
    
    fprintf(FID, '%s\n', '----------- RESPONSE----------');
    for i=1:length(Forcing)
        fprintf(FID, '%s\n', ['Forcing  ' '(' num2str(Forcing(i).n) ',' num2str(Forcing(i).m) ')']);
        fprintf(FID, '%s\n','k2 Love numbers');
        fprintf(FID, '%s\n','(n,m)           k_n^m');
        for j=1:length(y(i).n(:))
            if y(i).n(j)==Forcing(i).n && y(i).m(j)==Forcing(i).m
                k2=y(i).y(end,8,j)-1;
            else
                k2=y(i).y(end,8,j);
            end
            fprintf(FID, '%s\n',['(' num2str(y(i).n(j)) ',' num2str(y(i).m(j)) ')           '  num2str(k2,'%10.5e')   ]);
        end
    end
    fprintf(FID, '%s\n', '----------- ENERGY SPECTRA ----------');
    fprintf(FID, '%s\n','(n,m)           \dot{e}_n^m');
    for i=1:length(Energy_Spectra.n)
        fprintf(FID, '%s\n',['(' num2str(Energy_Spectra.n(i)) ',' num2str(Energy_Spectra.m(i)) ')           '  num2str(Energy_Spectra.energy_integral(i),'%10.5e')   ]);
    end
    if calc_E_contributions == 1
        save(out_file_name_E_contribution,'energy_contribution_s');
    end
    if save_energy_vec == 1
        save(out_file_name_E_s,'energy_s');
    end
end
end