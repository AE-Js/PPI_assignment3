%% GET_ENERGY
% AUTHOR: M. Rovira-Navarro 
% USE: obtain the energy spectra 
%% INPUT       
    % Numerics: Structure containing the numerical information
        %Numerics.Nr: number of radial points
        %Numerics.Nenergy: up to which degree the energy dissipation spectra is computed. Default 8
        %Numerics.perturbation_order: maximum order of the perturbation. Default 2
        %Numerics.solution_cutoff: cut off degree (if specified instead of perturbation_order)
        %Numerics.rheology_cutoff: up to which order the rheology is expanded. Default 0 (only leading terms)
        %Numerics.load_couplings: load coupling coefficients (1) or not (0). Default 1.
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
    %Forcing: Vector containing the forcing information
        % Forcing.Td: forcing period
        % Forcing.n: degree of the forcing 
        % Forcing.m: order of the forcing 
        % Forcing.F: amplitude of the component 
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
Nf=length(y); %number of forcings
n_sol=[];
m_sol=[];
for i=1:Nf
    n_sol=[n_sol y(i).n];
    m_sol=[m_sol y(i).m];
end
n_aux=unique(sort(abs(n_sol)));
m_sol_t=[];
n_sol_t=[];
k=1;
for i=1:length(n_aux)
    index=find(n_aux(i)==n_sol);
    m_aux=unique(abs(m_sol(index)));
    m_aux=sort(unique([m_aux -m_aux]));
    orders_l=length(m_aux);
    m_sol_t(k:k+orders_l-1)=m_aux;
    n_sol_t(k:k+orders_l-1)=n_aux(i);
    k=k+orders_l;
end

% Build the total solution 
y_total=zeros(Numerics.Nr+1,24,length(n_sol_t));
for i=1:length(n_sol_t)
    for j=1:Nf
        ind=find(n_sol_t(i)==y(j).n & m_sol_t(i)==y(j).m);
        y_total(:,1,i)=y(1).y(:,1,1); 
        if isempty(ind)==0
            y_total(:,2:end,i)=y_total(:,2:end,i) + Forcing(j).F*y(j).y(:,2:end,ind); 
        end
    end
end
%% strain and stress tensor and the - component 
stressP=zeros(Numerics.Nr+1,6,length(n_sol_t));
strainP=zeros(Numerics.Nr+1,6,length(n_sol_t));
stressN=zeros(Numerics.Nr+1,6,length(n_sol_t));
strainN=zeros(Numerics.Nr+1,6,length(n_sol_t));

% Fill the Stress and Strain tensors
for i=1:length(n_sol_t)
    index_auxN=find(n_sol_t==n_sol_t(i) & m_sol_t==-m_sol_t(i));
    stressP(:,:,i)=y_total(:,[14 15 16 17 18 13],i);
    strainP(:,:,i)=y_total(:,[20 21 22 23 24 19],i);
    stressN(:,:,i)=conj(y_total(:,[14 15 16 17 18 13],index_auxN)); 
    strainN(:,:,i)=conj(y_total(:,[20 21 22 23 24 19],index_auxN));
end

% Retrieve the radial points
r=y(1).y(:,1,1);

%% (2) GET COUPLINGS ENERGY SPECTRA
% obtain name of the coupling file
str_forc=[];
% Check whether there are any uniform layers
uniformlayers = [];
for ilayer=2:Numerics.Nlayers
    uniformlayers = [uniformlayers Interior_Model(ilayer).uniform];
end
uniformlayers = uniformlayers == 0;

% Generate the forcing term string
for i=1:length(Forcing)
    if i==length(Forcing)
        str_forc=[str_forc num2str(Forcing(i).n) '_' num2str(Forcing(i).m)];
    else
        str_forc=[str_forc num2str(Forcing(i).n) '_' num2str(Forcing(i).m) '__'];
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
    % Generate coupling file name
    coupling_file_name=[Numerics.coupling_file_location 'E_struct__Nrheomax__' num2str(Numerics.Nrheo_max) '__forc__' ...
        str_forc '__N__' num2str(Numerics.Nenergy) '__per' num2str(Numerics.perturbation_order) '.mat'];
    
    % look if there is a file that contains the couplings
    coupling_file_name_search=[Numerics.coupling_file_location 'E_struct__Nrheomax__*__forc__' str_forc '__N__*__per*.mat'];
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
        Nenergy_start=strfind(possible_couplings_files(i).name,'__N__')+5;
        Nenergy_end=strfind(possible_couplings_files(i).name,'__per')-1;
        Nenergy_file=str2num(possible_couplings_files(i).name(Nenergy_start:Nenergy_end));

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
    else %compute couplings 
        if verbose==1
            tic
            disp(['File ' coupling_file_name 'not found. Computing all coupling coefficients, this might take some time..'])
            Couplings = get_energy_couplings_all(Numerics.perturbation_order,Numerics.Nrheo_max,Numerics.Nenergy,Forcing,Numerics,'verbose');
            disp(['Time Spent: ' num2str(toc) 's'])
            disp([' File stored in: ' coupling_file_name])
        else
            Couplings = get_energy_couplings_all(Numerics.perturbation_order,Numerics.Nrheo_max,Numerics.Nenergy,Forcing,Numerics);
        end

        % Save the file
        save(coupling_file_name,'-struct','Couplings','-v7.3')
        
        % Retrieve the coupling coefficients that are required
        Couplings = retrieve_energy_couplings(n_sol_t,m_sol_t,Numerics.Nenergy,Couplings);
    end
    
    % Extract relevant information from struct
    EC = Couplings.EC;
    n_en = Couplings.n_en;
    m_en = Couplings.m_en;

    % Dubble-check whether the solution modes used in the coupling matrix are the
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
    coupling_file_name=[Numerics.coupling_file_location 'E__rheo__0_0__forc__' str_forc '__N__' ...
                       num2str(Numerics.Nenergy) '__per' num2str(Numerics.perturbation_order) '.mat'];
    
    % Check whether the file exists
    if isfile(coupling_file_name)==1
        % coupling file exist and does not need to be computed
        if verbose==1
            disp([' Energy couplings loaded from: ' coupling_file_name])
        end
        load(coupling_file_name);
    else
        if verbose==1
            disp([' Computing energy couplings...'])
            tic
        end
        [EC,n_en,m_en] = couplings_energy(n_sol_t,m_sol_t,Numerics.Nenergy);
        save(coupling_file_name,'EC','n_en','m_en');
        if verbose==1
            disp(['Time Spent: ' num2str(toc) 's'])
            disp([' File stored in: ' coupling_file_name])
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
energy=zeros(length(r),length(n_en));
energy_s=zeros(1,length(n_en)); 

% variable to store all the contributions 
if calc_E_contributions == 1
    energy_contribution=zeros(length(r),length(n_en),length(n_sol_t),length(n_sol_t));
    energy_contribution_s=zeros(1,length(n_en),length(n_sol_t),length(n_sol_t));
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
Energy_Spectra.n=n_en;
Energy_Spectra.m=m_en;
Energy_Spectra.energy_integral=energy_s;
Energy_Spectra.energy=energy;
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