%% GET_LOVE
% AUTHOR: M. Rovira-Navarro 
% USE: obtain the tidal love number spectra for an interior model with lateral variations
%% INPUT 
% Interior_Model: Structure containing the interior model information
        %Interior_Model(ilayer).R0: radius
            %R0(1) core radius (solid+liquid core)
            %R0(2) surface radius
        %Interior_Model(ilayer).rho0: layers density        
            %rho0(1) density of interior layer (solid+liquid core)
            %rho0(2) density of outermost solid layer
        %Interior_Model(ilayer).Delta_rho0: Density difference between the liquid core and overlying solid layer. If not provided it is computed assuming that the two innermost layers have rho0(1).
        %Interior_Model(ilayer).mu0: shear modulus of the the outermost layer
        %Interior_Model(ilayer).Ks0: bulk modulus of the outermost layer
        % Interio_Model.mu_variable: shear modulus variations
            %mu_variable(:,1): degree of variation 
            %mu_variable(:,2): order of variation 
            %mu_variable(:,3): amplitude of the variation (mu_l^m/mu^0_0)
        % Interio_Model.K_variable: bulk modulus variations 
            %K_variable(:,1): degree of variation 
            %K_variable(:,2): order of variation 
            %K_variable(:,3):  amplitude of bulk modulus variations (K_l^m/K^0_0)
        % Interio_Model.eta_variable: viscosity
            %eta_variable(:,1): degree of variation 
            %eta_variable(:,2): order of variation 
            %eta_variable(:,3):  amplitude of viscosity variations (eta_l^m/eta^0_0)    
        % Interior_Model(ilayer).rheology_variable: rheology variable (assigned inside the code)
            %rheology_variable(:,1): degree of variation 
            %rheology_variable(:,2): order of variation 
            %rheology_variable(:,3):  amplitude of bulk modulus variations (K_l^m/K^0_0)  
            %rheology_variable(:,4):  amplitude of complex shear modulus variations (mu_l^m/mu^0_0) 
        % Interior_Model(ilayer).muC: complex shear modulus (assigned inside the code)   

    % Forcing: Structure containing forcing information
        % Forcing.Td: forcing period
        % Forcing.n: degree of the forcing 
        % Forcing.m: order of the forcing 
        % Forcing.F: amplitude of the component 

    % Numerics: Structure containing the numerical information
        %Numerics.Nr: number of radial points
        %Numerics.perturbation_order: maximum order of the perturbation. Default 2
        %Numerics.solution_cutoff: cut off degree (if specified instead of perturbation_order)
        %Numerics.rheology_cutoff: determines which terms of the rheology are included (only relevant for viscoelastic).
            % terms with log10(mu_n^m)-log10(mu_n^m(leading))>=-Numerics.rheology_cutoff are included 
            % Default 0 (only leading terms)
        %Numerics.load_couplings: 
                % (0) compute coupling coefficients
                % (1) load coupling coefficients from file that contains
                % exactly the same coupling coefficintes 
                % (2) load coupling coefficintes from a file.

    % optional variables 
        %plot_rheology_spectra: plot rheology spectra (1)
        %verbose: print information in screen (1)
        %out_file: print output in file named out_file
%% OUTPUT 
    % Love_Spectra: Love number spectra
        % Love_Spectra.nf: degree of the forcing 
        % Love_Spectra.mf: order of the forcing 
        % Love_Spectra.n: degree of the solution
        % Love_Spectra.m: order of the solution
        % Love_Spectra.k: gravity Love numbers
        % Love_Spectra.h: radial displacement Love numbers
    % y_rad: radial functions
        % y_rad.nf: degree of the forcing
        % y_rad.mf: order of the forcing 
        % y_rad.n: degree of the solution
        % y_rad.m: order of the solution 
        % y_rad.y(radial_point,X,mode) 
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
   % varargout{1}: Interior_Model updated with variables non-dimensionalized         
%% FUNCTION ------------------------------------------------------------------   
function [Love_Spectra,y_rad,varargout] = get_Love(Interior_Model,Forcing,Numerics,varargin)
verbose = 0; 
out_file = 0; 
save_solution_vec = 0;
%% (0) Optional Inputs Defaults -----------------------------------------------------
if isfield(Numerics,'Nr')==0
    Numerics.Nr=100;
end
if isfield(Numerics,'rheology_cutoff')==0
    Numerics.rheology_cutoff=2;
end
if isfield(Numerics,'load_couplings')==0
    Numerics.load_couplings=1;
end
if isfield(Numerics,'perturbation_order')==0
    Numerics.perturbation_orders=2;
end
if isfield(Numerics,'solution_cutoff')==0
    Numerics.perturbation_orders=0;
    Numerics.solution_cutoff=8;
end

for k = 1:length(varargin)
    if strcmpi(varargin{k},'verbose')
        verbose=1; 
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
    if strcmpi(varargin{k},'save_solution')
        save_solution_vec=varargin{k+1}; 
        varargin{k+1}=[]; 
        varargin{k}=[];
    end
end

%% (0.1) File for printing -------------------------------------------------------------------
% Obtain the highest degree of rheology across the layers that is used
Numerics.Nrheo_max = 0;
for ilayer=2:Numerics.Nlayers
    prevmax = max(Interior_Model(ilayer).rheology_variable(:,1));
    Numerics.Nrheo_max = max([Numerics.Nrheo_max prevmax]);
end

if out_file~=0
    out_file_name=[save_location out_file '__' num2str(Numerics.Nrheo_max) '_per' num2str(Numerics.rheology_cutoff) ...
        '__forc__' num2str(Forcing.n) '_' num2str(Forcing.m) '_per' num2str(Numerics.perturbation_order) '.txt'];
    FID=fopen(out_file_name,'w');
    out_file_name_y=[save_location 'y__' num2str(Numerics.Nrheo_max) '_per' num2str(Numerics.rheology_cutoff) ...
        '__forc__' num2str(Forcing.n) '_' num2str(Forcing.m) '_per' num2str(Numerics.perturbation_order) '.mat'];
end

%% (2) OBTAIN COUPLINGS --------------------------------------------------------------
if verbose==1
    disp('Obtaining Couplings....')
    tic
end

% Check whether there are any uniform layers
uniformlayers = [];
for ilayer=2:Numerics.Nlayers
    uniformlayers = [uniformlayers Interior_Model(ilayer).uniform];
end
uniformlayers = uniformlayers == 0;

% If there is only one non-uniform layer obtain couplings
if any(uniformlayers) 
    % Generate coupling file name
    coupling_file_name=[Numerics.coupling_file_location 'L_struct__Nrheomax__' num2str(Numerics.Nrheo_max) '__forc__' ...
        num2str(Forcing.n) '_' num2str(Forcing.m) '_per' num2str(Numerics.perturbation_order) '.mat'];

    % Look if there is a file that contains the couplings
    coupling_file_name_search=[Numerics.coupling_file_location 'L_struct__Nrheomax__*__forc__' num2str(Forcing.n) ...
                              '_' num2str(Forcing.m) '_per*.mat'];
    possible_couplings_files=dir(coupling_file_name_search);
    file_found=0;
    i=1; 
    potential_coupling_files = [];
    rheo_max_file_list = [];

    % Select the smallest possible file that contains all the couplings
    while i<=length(possible_couplings_files)
        per_start=strfind(possible_couplings_files(i).name,'_per')+4;
        per_end=strfind(possible_couplings_files(i).name,'.mat')-1;
        perturbation_order_file=str2num(possible_couplings_files(i).name(per_start:per_end));
        rheo_start=strfind(possible_couplings_files(i).name,'Nrheomax__')+10;
        rheo_end=strfind(possible_couplings_files(i).name,'__forc')-1;
        Nrheomax_file=str2num(possible_couplings_files(i).name(rheo_start:rheo_end));

        % Check whether the file fullfills the criteria and add it to the list
        if Nrheomax_file>=Numerics.Nrheo_max && perturbation_order_file>=Numerics.perturbation_order 
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
    
    % remove terms with 0 amplitude 
    % fill an array with all the used degree and orders (contains duplicates)
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
    
    % Load the coupling file
    if file_found==1 %coupling file exist 
        if verbose==1
            disp([' Couplings Loaded from: ' convertStringsToChars(coupling_file_name)])
        end

        % Retrieve the coupling coefficients that are required
        tic
        Couplings = retrieve_couplings_from_file(Numerics.perturbation_order,rheo_degree_orders,Numerics.Nrheo_max,Forcing,coupling_file_name);
        disp(['Loading coupling file took: ' num2str(toc) ' seconds'])
    else %compute couplings 
        if verbose==1
            tic
            disp(['File ' coupling_file_name 'not found. Computing all coupling coefficients, this might take some time..'])
            Couplings = get_couplings_all(Numerics.perturbation_order,Numerics.Nrheo_max,Forcing,Numerics,'verbose');
            disp(['Time Spent: ' num2str(toc) 's'])
            disp([' File stored in: ' coupling_file_name])
        else
            Couplings = get_couplings_all(Numerics.perturbation_order,Numerics.Nrheo_max,Forcing,Numerics);
        end

        % Save couplings for all terms
        save(coupling_file_name,'-struct','Couplings','-v7.3')

        % Retrieve the coupling coefficients that are required
        Couplings = retrieve_couplings(Numerics.perturbation_order,rheo_degree_orders,Numerics.Nrheo_max,Forcing,Couplings);
    end
    
else
    if verbose==1
        disp('Only uniform layers so coupling matrix is zero')
    end
    Couplings.n_s=Forcing.n;
    Couplings.m_s=Forcing.m;
    Couplings.order=0;
    Couplings.Coup=zeros(1,1,27,1); 
end

Nsol=length(Couplings.n_s);

%% PRINT MODEL INFORMATION IN SCREEN 
if verbose==1
    disp('------------- INTERIOR MODEL -----------')
    disp('AVERAGE PROPERTIES (dimensional)')
    disp(['Layer#     R[m]    rho[kg.m^{-3}]    mu[Pa]    K[Pa]    eta[Pa.s]' ])
    disp(['1    ' num2str(Interior_Model(1).R0,'%10.5e') '    ' num2str(Interior_Model(1).rho0,'%10.5e') '    0'  '    -'   ])
    for ilayer=2:Numerics.Nlayers
        disp([num2str(ilayer) '    ' num2str(Interior_Model(ilayer).R0,'%10.5e') '    ' num2str(Interior_Model(ilayer).rho0,'%10.5e') '    ' num2str(Interior_Model(ilayer).mu0,'%10.5e')  '    ' num2str(Interior_Model(ilayer).Ks0,'%10.5e') ' ' num2str(Interior_Model(ilayer).eta0,'%10.5e')   ])
    end
    disp('AVERAGE PROPERTIES (non-dimensional)')
    disp(['Layer#     R[-]    rho[-]    mu[-]    K [-]     eta[-]' ])
    disp(['1    ' num2str(Interior_Model(1).R) '    ' num2str(Interior_Model(1).rho) '    0'  '    -'   ])
    for ilayer=2:Numerics.Nlayers
        disp([num2str(ilayer) '    ' num2str(Interior_Model(ilayer).R) '    ' num2str(Interior_Model(ilayer).rho) '    ' num2str(Interior_Model(ilayer).mu)  '    ' num2str(Interior_Model(ilayer).Ks) ' ' num2str(Interior_Model(ilayer).eta)   ])
    end
    disp('RHEOLOGY  VARIATIONS')
    disp('Shear Modulus')
    for ilayer=2:Numerics.Nlayers
        if isfield(Interior_Model(ilayer),'mu_variable')==1
            disp('Layer#      (n,m)           amplitude[mu_n^m/mu_0^0]')
            for i=1:size(Interior_Model(ilayer).mu_variable,1)
                disp([num2str(ilayer) '    (' num2str(Interior_Model(ilayer).mu_variable(i,1)) ',' num2str(Interior_Model(ilayer).mu_variable(i,2)) ')   '  num2str(Interior_Model(ilayer).mu_variable(i,3),'%10.5e') ])
            end
        else
            disp('None')
        end 
    end
    
    disp('Bulk Modulus')
    for ilayer=2:Numerics.Nlayers
        if isfield(Interior_Model,'K_variable')==1
            disp('Layer#      (n,m)           amplitude[K_n^m/K_0^0]')
            for i=1:size(Interior_Model(ilayer).K_variable,1)
                disp([num2str(ilayer) '    (' num2str(Interior_Model(ilayer).K_variable(i,1)) ',' num2str(Interior_Model(ilayer).K_variable(i,2)) ')   '  num2str(Interior_Model(ilayer).K_variable(i,3),'%10.5e') ])
            end
        else
            disp('None')
        end
    end

    disp('Viscosity')
    for ilayer=2:Numerics.Nlayers
        if isfield(Interior_Model,'eta_variable')==1
            disp('Layer#      (n,m)          amplitude[eta_n^m/eta_0^0]')
            for i=1:size(Interior_Model(ilayer).eta_variable,1)
                disp([num2str(ilayer) '    (' num2str(Interior_Model(ilayer).eta_variable(i,1)) ',' num2str(Interior_Model(ilayer).eta_variable(i,2)) ')   '  num2str(Interior_Model(ilayer).eta_variable(i,3),'%10.5e') ])
            end
        else
            disp('None')
        end
    end
    
    disp('Complex Shear Modulus')
    for ilayer=2:Numerics.Nlayers
        if isfield(Interior_Model,'rheology_variable')==1
            disp('Layer#      (n,m)         amplitude[\hat\mu_n^m/mu_0^0]')
            for i=1:size(Interior_Model(ilayer).rheology_variable,1)
                disp([num2str(ilayer) '    (' num2str(Interior_Model(ilayer).rheology_variable(i,1)) ',' num2str(Interior_Model(ilayer).rheology_variable(i,2)) ')    '  num2str(Interior_Model(ilayer).rheology_variable(i,4),'%10.5e') ])
            end
        else
            disp('None')
        end
    end
    disp(' ')
    disp('----------- FORCING ----------')
    disp('TIDAL POTENTIAL')
    disp('Period [s]')
    disp(num2str(Forcing.Td,'%10.5e'))
    disp('(n,m)')
    disp(['(' num2str(Forcing.n) ',' num2str(Forcing.m) ')'])
    disp(' ')
    disp('-----------  RESPONSE SPECTRUM ----------- ')
    for j=0:Numerics.perturbation_order
    index_order=find(Couplings.order==j);
    str_h=[];
    for k=1:length(index_order)
        str_h=[str_h '(' num2str(Couplings.n_s(index_order(k))) ',' num2str(Couplings.m_s(index_order(k))) '),   '];
    end
    disp([num2str(j) 'th Order Modes'])
    disp([num2str(length(index_order)) ' modes'])
    disp(str_h)
    end
    disp('----------- NUMERICAL INFORMATION ----------')
    disp(['Number of Modes ' num2str(Nsol)])
    disp(['Radial Points '  num2str(Numerics.Nr)])
end
%% PRINT MODEL INFORMATION IN FILE
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
    fprintf(FID, '%s\n', num2str(Forcing.Td,'%10.5e'));
    fprintf(FID, '%s\n', '(n,m)');
    fprintf(FID, '%s\n', ['(' num2str(Forcing.n) ',' num2str(Forcing.m) ')']);
    fprintf(FID, '%s\n', ' ');
    fprintf(FID, '%s\n', '-----------  RESPONSE SPECTRUM ----------- ');
    for j=0:Numerics.perturbation_order
    index_order=find(Couplings.order==j);
    str_h=[];
    for k=1:length(index_order)
        str_h=[str_h '(' num2str(Couplings.n_s(index_order(k))) ',' num2str(Couplings.m_s(index_order(k))) '),   '];
    end
    fprintf(FID, '%s\n',[num2str(j) 'th Order Modes']);
    fprintf(FID, '%s\n',[num2str(length(index_order)) ' modes']);
    fprintf(FID, '%s\n',str_h);
    end
    fprintf(FID, '%s\n', '----------- NUMERICAL INFORMATION ----------');
    fprintf(FID, '%s\n', ['Number of Modes ' num2str(Nsol)]);
    fprintf(FID, '%s\n', ['Radial Points '  num2str(Numerics.Nr)]);
end
%% (3) PROPAGATE AND OBTAIN SOLUTION ------------------------------------------------------------
if verbose==1
    tStart = tic;
end
if Numerics.parallel_sol == 0 && Numerics.parallel_gen == 1
    y_sol = get_solution_parallel(Interior_Model,Forcing,Numerics,Couplings,verbose,out_file);
else
    y_sol = get_solution(Interior_Model,Forcing,Numerics,Couplings,verbose,out_file);
end
if verbose==1
    disp(['Time Spent: ' num2str(toc(tStart)) 's'])
end
%% (4) PRINT SOME INFORMATION 
if verbose==1
    Nmodes=length(Couplings.n_s);
    disp('#####################################')
    disp('---------- RESPONSE ----------------')
    disp('k2 Love numbers')
    disp('(n,m)           k_n^m')
    for i=1:Nmodes
        if Couplings.n_s(i)==Forcing.n && Couplings.m_s(i)==Forcing.m
            k2=y_sol(end,8,i)-1;
        else
            k2=y_sol(end,8,i);
        end
        disp(['(' num2str(Couplings.n_s(i)) ',' num2str(Couplings.m_s(i)) ')           '  num2str(k2,'%10.5e')   ])
    end
end
if out_file~=0
    Nmodes=length(Couplings.n_s);
    fprintf(FID, '%s\n','#####################################');
    fprintf(FID, '%s\n','---------- RESPONSE ----------------');
    fprintf(FID, '%s\n','k2 Love numbers');
    fprintf(FID, '%s\n','(n,m)           k_n^m');
    for i=1:Nmodes
        if Couplings.n_s(i)==Forcing.n && Couplings.m_s(i)==Forcing.m
            k2=y_sol(end,8,i)-1;
        else
            k2=y_sol(end,8,i);
        end
        fprintf(FID, '%s\n',['(' num2str(Couplings.n_s(i)) ',' num2str(Couplings.m_s(i)) ')           '  num2str(k2,'%10.5e')   ]);
    end
    if save_solution_vec == 1
        save(out_file_name_y,"y_sol");
    end
end
%% (5) Rearrange solution 
Nmodes=length(Couplings.n_s);
y_rad.y=y_sol; 
Love_Spectra.nf=Forcing.n;
Love_Spectra.mf=Forcing.m;
y_rad.nf=Forcing.n;
y_rad.mf=Forcing.m;
for i=1:Nmodes
    Love_Spectra.n(i)=Couplings.n_s(i);
    Love_Spectra.m(i)=Couplings.m_s(i);
    Love_Spectra.order(i)=Couplings.order(i);
    if Couplings.n_s(i)==Forcing.n && Couplings.m_s(i)==Forcing.m
        Love_Spectra.k(i)=y_sol(end,8,i)-1;
    else
        Love_Spectra.k(i)=y_sol(end,8,i);
    end
    Love_Spectra.h(i)=-Interior_Model(end).gs*y_sol(end,2,i);
    Love_Spectra.l(i)=-Interior_Model(end).gs*y_sol(end,3,i);
    y_rad.n(i)=Couplings.n_s(i);
    y_rad.m(i)=Couplings.m_s(i);
end

end
