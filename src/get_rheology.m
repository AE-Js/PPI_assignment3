%% Description 
% Computes the complex rheology of the body based on either percentage
% variations of the average rheology or full, real rheology fields. Also
% contains the option to run using legacy values in order to compare
% against previous results. 
% 
% -  It first normalizes the input values. 
% -  Then (potentially) converts lateral rheology variations based on
%    percentages into SPH.
% -  It then converts the SPH, either from percentage based variations or
%    real rheology fields, into a lat-lon field.
% -  Using the lat-lon it computes the complex rheology in the Fourier
%    domain and filters the relevant modes based on relative magnitude. 
%
%%
function [Interior_Model] = get_rheology(Interior_Model,Numerics,Forcing,varargin)
G=6.67430E-11;
plot_rheology_spectra=0; 
calculate_G = true;
verbose=0; 
out_file=0; 

% Determines whether sin components are taken into account when expanding 
% the rheology into the fourier domain. Default is 1.
sine_component_factor = 1; % 0 or 1 
%% (0) Optional Inputs Defaults -----------------------------------------------------
if isfield(Numerics,'Nr')==0
    Numerics.Nr=100;
end
if isfield(Numerics,'rheology_cutoff')==0
    Numerics.rheology_cutoff=2;
end
if isfield(Numerics,'minimum_rheology_value')==0
    % If a full rheology field is given as input the highest rheology mode must be large
    % enough to not crash due to numerical noise imitating rheology modes.
    Numerics.minimum_rheology_value=-14;
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
    if strcmpi(varargin{k},'plot_rheology_spectra')
        plot_rheology_spectra=varargin{k+1}; 
        varargin{k+1}=[]; 
        varargin{k}=[];
    end
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
    if strcmpi(varargin{k},'calculate_G')
        calculate_G=varargin{k+1}; 
        varargin{k+1}=[]; 
        varargin{k}=[];
    end
    if strcmpi(varargin{k},'include_sine')
        sine_component_factor=varargin{k+1}; 
        varargin{k+1}=[]; 
        varargin{k}=[];
    end
end

%% (0) Nondimensional ------------------------------------------------------------------
% Non-dimensionalize the parameters based on the values in the outer layer

% If the normalized gravitational constant is already provided, use that
if isfield(Interior_Model(end),"Gg")==1 
    Gg = Interior_Model(end).Gg;
else
    Gg = G*(Interior_Model(end).R0*1000)^2*Interior_Model(end).rho0^2/Interior_Model(end).mu0;
end

% If calculate_G=false a legacy version of the code is used
if calculate_G && isfield(Interior_Model(end),"Gg")==0
    Gg = G*(Interior_Model(end).R0*1000)^2*Interior_Model(end).rho0^2/Interior_Model(end).mu0;
else
    % Test values for Gg, Legacy values (only used when comparing against old data)
    r_ratio = 0.53;
    rho_r = 1.59;
    mu_eff = 5.2;
    rho_av = rho_r*r_ratio^3+(1-r_ratio^3);
    Gg = 3/(4*pi)/mu_eff/rho_av^2;
    % G=0.0388;
    % Gg = 0.0388;
end

for ilayer=1:Numerics.Nlayers
    if ilayer==1
        Interior_Model(ilayer).R = Interior_Model(ilayer).R0/Interior_Model(end).R0;
        Interior_Model(ilayer).rho = Interior_Model(ilayer).rho0/Interior_Model(end).rho0;
        if isfield(Interior_Model(1),"Delta_rho0")==0
            Interior_Model(ilayer).Delta_rho0= (Interior_Model(1).rho0-Interior_Model(2).rho0);
            Interior_Model(ilayer).Delta_rho = (Interior_Model(1).rho0-Interior_Model(2).rho0)/Interior_Model(end).rho0;
            % Else the Delta_rho0 is provided by the user 
        end
        Interior_Model(ilayer).Delta_rho = Interior_Model(ilayer).Delta_rho0/Interior_Model(end).rho0;
        Interior_Model(ilayer).Gg0 = G; % G
        Interior_Model(ilayer).Gg  = Gg; % Non-dimensionalised G  

        % Build the (normalized) average density of the body
        rho_av0 = 4/3*pi*Interior_Model(ilayer).R0^3*Interior_Model(ilayer).rho0;
        rho_av = Interior_Model(ilayer).rho*Interior_Model(ilayer).R^3;
    else
        Interior_Model(ilayer).R = Interior_Model(ilayer).R0/Interior_Model(end).R0;
        Interior_Model(ilayer).rho = Interior_Model(ilayer).rho0/Interior_Model(end).rho0;
        Interior_Model(ilayer).Ks = Interior_Model(ilayer).Ks0/Interior_Model(end).mu0;
        Interior_Model(ilayer).mu = Interior_Model(ilayer).mu0/Interior_Model(end).mu0;
        if isfield(Interior_Model(1),"MaxTime")==1
            Interior_Model(ilayer).eta=Interior_Model(ilayer).MaxTime/(2*pi);
            Interior_Model(ilayer).eta0=Interior_Model(ilayer).MaxTime/(2*pi);
        else
            Interior_Model(ilayer).eta = Interior_Model(ilayer).eta0/(Interior_Model(end).mu0*Forcing(1).Td);
            % 1\omega is used to non-dimensionalize time
            Interior_Model(ilayer).MaxTime = 2*pi*Interior_Model(ilayer).eta0/Interior_Model(ilayer).mu0/Forcing(1).Td; 
        end
        Interior_Model(ilayer).Gg0 = G; % G
        Interior_Model(ilayer).Gg  = Gg; % Non-dimensionalised G

        % Density difference with layer below taken such that it should be positive
        Interior_Model(ilayer).Delta_rho = Interior_Model(ilayer-1).rho-Interior_Model(ilayer).rho; 

        % Add mass of the layer to obtain the (normalized) average density
        rho_av = rho_av + Interior_Model(ilayer).rho*(Interior_Model(ilayer).R^3-Interior_Model(ilayer-1).R^3);
        rho_av0 = rho_av0 + 4/3*pi*Interior_Model(ilayer).rho0*(Interior_Model(ilayer).R0^3-Interior_Model(ilayer-1).R0^3);
    end
end

% Compute average density, required to get h and l Love numbers
rho_av0 = rho_av0/(4/3*pi*Interior_Model(end).R0^3);
for ilayer=1:Numerics.Nlayers
    Interior_Model(ilayer).rho_av0 = rho_av0;
    Interior_Model(ilayer).rho_av = rho_av;
    Interior_Model(ilayer).gs0 = 4/3*pi*Interior_Model(end).Gg0*rho_av0*Interior_Model(end).R0*1000; % Correct for radius in km
    Interior_Model(ilayer).gs = 4/3*pi*Interior_Model(end).Gg*rho_av*Interior_Model(end).R;
end

%% convert lateral percentage variations to SPH
% this part converts the variations given in real spherical harmonics and percentage into complex spherical harmonics
l_max=30; % Higher is better but slower
for ilayer=2:Numerics.Nlayers
    % Retrieve relevant values from Interior_Model
    nR = Interior_Model(ilayer).nR;
    mR = Interior_Model(ilayer).mR;
    if length(nR) > 1
        variable_mu_p = Interior_Model(ilayer).variable_mu_p;
        variable_eta_p = Interior_Model(ilayer).variable_eta_p;
        variable_K_p = Interior_Model(ilayer).variable_K_p;
    
        Ynm_stokes.clm=zeros(2*l_max,2*l_max);
        Ynm_stokes.slm=zeros(2*l_max,2*l_max);
        Ynm_stokes.lmax=2*l_max-1;
        k=1;
        for i=1:length(nR)
            nR_fixed = nR(i);
            mR_fixed = mR(i);
            Ynm_stokes.clm(nR_fixed+1,mR_fixed+1)=1;
            [Ynm_zlonlat] = SPH_LatLon(Ynm_stokes);
            Ynm_z=Ynm_zlonlat.z;
            Delta=max(Ynm_z(:))-min(Ynm_z(:));
            if mR_fixed==0
                % shear modulus 
                Interior_Model(ilayer).mu_variable(k,1)=nR_fixed;
                Interior_Model(ilayer).mu_variable(k,2)=mR_fixed;
                Interior_Model(ilayer).mu_variable(k,3)=variable_mu_p(i)/100/Delta;
                % bulk modulus
                Interior_Model(ilayer).K_variable(k,1)=nR_fixed;
                Interior_Model(ilayer).K_variable(k,2)=mR_fixed;
                Interior_Model(ilayer).K_variable(k,3)=variable_K_p(i)/100/Delta;
                %viscosity 
                Interior_Model(ilayer).eta_variable(k,1)=nR_fixed;
                Interior_Model(ilayer).eta_variable(k,2)=mR_fixed;
                Interior_Model(ilayer).eta_variable(k,3)=variable_eta_p(i)/100/Delta;
                k = k+1;
            elseif mR_fixed>0
                % shear modulus 
                Interior_Model(ilayer).mu_variable(k,1)=nR_fixed;
                Interior_Model(ilayer).mu_variable(k+1,1)=nR_fixed;
                Interior_Model(ilayer).mu_variable(k,2)=mR_fixed;
                Interior_Model(ilayer).mu_variable(k+1,2)=-mR_fixed;    
                Interior_Model(ilayer).mu_variable(k,3)=sqrt(2)/2*variable_mu_p(i)/100/Delta;
                Interior_Model(ilayer).mu_variable(k+1,3)=(-1)^mR_fixed*sqrt(2)/2*variable_mu_p(i)/100/Delta;
                %bulk modulus
                Interior_Model(ilayer).K_variable(k,1)=nR_fixed;
                Interior_Model(ilayer).K_variable(k+1,1)=nR_fixed;
                Interior_Model(ilayer).K_variable(k,2)=mR_fixed;
                Interior_Model(ilayer).K_variable(k+1,2)=-mR_fixed;    
                Interior_Model(ilayer).K_variable(k,3)=sqrt(2)/2*variable_K_p(i)/100/Delta;
                Interior_Model(ilayer).K_variable(k+1,3)=(-1)^mR_fixed*sqrt(2)/2*variable_K_p(i)/100/Delta;
                % viscosity 
                Interior_Model(ilayer).eta_variable(k,1)=nR_fixed;
                Interior_Model(ilayer).eta_variable(k+1,1)=nR_fixed;
                Interior_Model(ilayer).eta_variable(k,2)=mR_fixed;
                Interior_Model(ilayer).eta_variable(k+1,2)=-mR_fixed;    
                Interior_Model(ilayer).eta_variable(k,3)=sqrt(2)/2*variable_eta_p(i)/100/Delta;
                Interior_Model(ilayer).eta_variable(k+1,3)=(-1)^mR_fixed*sqrt(2)/2*variable_eta_p(i)/100/Delta;
                k = k+2;
            else
                % shear modulus 
                Interior_Model(ilayer).mu_variable(k,1)=nR_fixed;
                Interior_Model(ilayer).mu_variable(k+1,1)=nR_fixed;
                Interior_Model(ilayer).mu_variable(k,2)=-mR_fixed;
                Interior_Model(ilayer).mu_variable(k+1,2)=mR_fixed;    
                Interior_Model(ilayer).mu_variable(k,3)=-1i*(-1)^mR_fixed*sqrt(2)/2*variable_mu_p(i)/100/Delta;
                Interior_Model(ilayer).mu_variable(k+1,3)=1i*sqrt(2)/2*variable_mu_p(i)/100/Delta;
                % bulk modulus
                Interior_Model(ilayer).K_variable(k,1)=nR_fixed;
                Interior_Model(ilayer).K_variable(k+1,1)=nR_fixed;
                Interior_Model(ilayer).K_variable(k,2)=-mR_fixed;
                Interior_Model(ilayer).K_variable(k+1,2)=mR_fixed;    
                Interior_Model(ilayer).K_variable(k,3)=-1i*(-1)^mR_fixed*sqrt(2)/2*variable_K_p(i)/100/Delta;
                Interior_Model(ilayer).K_variable(k+1,3)=1i*sqrt(2)/2*variable_K_p(i)/100/Delta;
                %viscosity
                Interior_Model(ilayer).eta_variable(k,1)=nR_fixed;
                Interior_Model(ilayer).eta_variable(k+1,1)=nR_fixed;
                Interior_Model(ilayer).eta_variable(k,2)=-mR_fixed;
                Interior_Model(ilayer).eta_variable(k+1,2)=mR_fixed;    
                Interior_Model(ilayer).eta_variable(k,3)=-1i*(-1)^mR_fixed*sqrt(2)/2*variable_eta_p(i)/100/Delta;
                Interior_Model(ilayer).eta_variable(k+1,3)=1i*sqrt(2)/2*variable_eta_p(i)/100/Delta;
                k = k+2;
            end
        end
    else
        variable_mu_p = Interior_Model(ilayer).variable_mu_p;
        variable_eta_p = Interior_Model(ilayer).variable_eta_p;
        variable_K_p = Interior_Model(ilayer).variable_K_p;
    
        Ynm_stokes.clm=zeros(2*l_max,2*l_max);
        Ynm_stokes.slm=zeros(2*l_max,2*l_max);
        Ynm_stokes.lmax=2*l_max-1;
        Ynm_stokes.clm(nR+1,mR+1)=1;
        [Ynm_zlonlat] = SPH_LatLon(Ynm_stokes);
        Ynm_z=Ynm_zlonlat.z;
        Delta=max(Ynm_z(:))-min(Ynm_z(:));
    
        if mR==0
            % shear modulus 
            Interior_Model(ilayer).mu_variable(1,1)=nR;
            Interior_Model(ilayer).mu_variable(1,2)=mR;
            Interior_Model(ilayer).mu_variable(1,3)=variable_mu_p/100/Delta;
            % bulk modulus
            Interior_Model(ilayer).K_variable(1,1)=nR;
            Interior_Model(ilayer).K_variable(1,2)=mR;
            Interior_Model(ilayer).K_variable(1,3)=variable_K_p/100/Delta;
            %viscosity 
            Interior_Model(ilayer).eta_variable(1,1)=nR;
            Interior_Model(ilayer).eta_variable(1,2)=mR;
            Interior_Model(ilayer).eta_variable(1,3)=variable_eta_p/100/Delta;
        elseif mR>0
            % shear modulus 
            Interior_Model(ilayer).mu_variable(1,1)=nR;
            Interior_Model(ilayer).mu_variable(2,1)=nR;
            Interior_Model(ilayer).mu_variable(1,2)=mR;
            Interior_Model(ilayer).mu_variable(2,2)=-mR;    
            Interior_Model(ilayer).mu_variable(1,3)=sqrt(2)/2*variable_mu_p/100/Delta;
            Interior_Model(ilayer).mu_variable(2,3)=(-1)^mR*sqrt(2)/2*variable_mu_p/100/Delta;
            %bulk modulus
            Interior_Model(ilayer).K_variable(1,1)=nR;
            Interior_Model(ilayer).K_variable(2,1)=nR;
            Interior_Model(ilayer).K_variable(1,2)=mR;
            Interior_Model(ilayer).K_variable(2,2)=-mR;    
            Interior_Model(ilayer).K_variable(1,3)=sqrt(2)/2*variable_K_p/100/Delta;
            Interior_Model(ilayer).K_variable(2,3)=(-1)^mR*sqrt(2)/2*variable_K_p/100/Delta;
            % viscosity 
            Interior_Model(ilayer).eta_variable(1,1)=nR;
            Interior_Model(ilayer).eta_variable(2,1)=nR;
            Interior_Model(ilayer).eta_variable(1,2)=mR;
            Interior_Model(ilayer).eta_variable(2,2)=-mR;    
            Interior_Model(ilayer).eta_variable(1,3)=sqrt(2)/2*variable_eta_p/100/Delta;
            Interior_Model(ilayer).eta_variable(2,3)=(-1)^mR*sqrt(2)/2*variable_eta_p/100/Delta;
        else
            % shear modulus 
            Interior_Model(ilayer).mu_variable(1,1)=nR;
            Interior_Model(ilayer).mu_variable(2,1)=nR;
            Interior_Model(ilayer).mu_variable(1,2)=-mR;
            Interior_Model(ilayer).mu_variable(2,2)=mR;    
            Interior_Model(ilayer).mu_variable(1,3)=-1i*(-1)^mR*sqrt(2)/2*variable_mu_p/100/Delta;
            Interior_Model(ilayer).mu_variable(2,3)=1i*sqrt(2)/2*variable_mu_p/100/Delta;
            % bulk modulus
            Interior_Model(ilayer).K_variable(1,1)=nR;
            Interior_Model(ilayer).K_variable(2,1)=nR;
            Interior_Model(ilayer).K_variable(1,2)=-mR;
            Interior_Model(ilayer).K_variable(2,2)=mR;    
            Interior_Model(ilayer).K_variable(1,3)=-1i*(-1)^mR*sqrt(2)/2*variable_K_p/100/Delta;
            Interior_Model(ilayer).K_variable(2,3)=1i*sqrt(2)/2*variable_K_p/100/Delta;
            %viscosity
            Interior_Model(ilayer).eta_variable(1,1)=nR;
            Interior_Model(ilayer).eta_variable(2,1)=nR;
            Interior_Model(ilayer).eta_variable(1,2)=-mR;
            Interior_Model(ilayer).eta_variable(2,2)=mR;    
            Interior_Model(ilayer).eta_variable(1,3)=-1i*(-1)^mR*sqrt(2)/2*variable_eta_p/100/Delta;
            Interior_Model(ilayer).eta_variable(2,3)=1i*sqrt(2)/2*variable_eta_p/100/Delta;
        end
    end
end

%% Check for uniformity and elasticity
for ilayer=2:Numerics.Nlayers
    Interior_Model(ilayer).elastic=0;
    Interior_Model(ilayer).uniform=1;
    if isnan(Interior_Model(ilayer).eta0) || isnan(Interior_Model(ilayer).eta)
        Interior_Model(ilayer).elastic=1; 
    end
    if isfield(Interior_Model(ilayer),'eta_variable')==1 
        if max(abs(Interior_Model(ilayer).eta_variable(:,3)))>0
            Interior_Model(ilayer).uniform=0;
        end
    end
    if isfield(Interior_Model(ilayer),"eta_latlon")
        if ~isempty(Interior_Model(ilayer).eta_latlon)
            Interior_Model(ilayer).uniform=0;
        end
    end
    if isfield(Interior_Model(ilayer),'mu_variable')==1
        if max(abs(Interior_Model(ilayer).mu_variable(:,3)))>0
            Interior_Model(ilayer).uniform=0;
        end
    end
    if isfield(Interior_Model(ilayer),"mu_latlon")
        if ~isempty(Interior_Model(ilayer).mu_latlon)
            Interior_Model(ilayer).uniform=0;
        end
    end
    if isfield(Interior_Model(ilayer),'K_variable')==1
        if max(abs(Interior_Model(ilayer).K_variable(:,3)))>0
            Interior_Model(ilayer).uniform=0;
        end
    end
    if isfield(Interior_Model(ilayer),"k_latlon")
        if ~isempty(Interior_Model(ilayer).k_latlon)
            Interior_Model(ilayer).uniform=0;
        end
    end
end

%% (1) GET RHEOLOGY ------------------------------------------------------------------
l_max = 25; %the higher the better the transformation but also slower 

% If there is a direct insertion of a lat-lon grid, that lmax value is used
if isfield(Interior_Model(ilayer),"mu_latlon")
    if ~isempty(Interior_Model(ilayer).mu_latlon)
        l_max = Interior_Model(ilayer).mu_latlon.lmax;
    end
end
if isfield(Interior_Model(ilayer),"eta_latlon")
    if ~isempty(Interior_Model(ilayer).eta_latlon)
        l_max = Interior_Model(ilayer).eta_latlon.lmax;
    end
end
if isfield(Interior_Model(ilayer),"k_latlon")
    if ~isempty(Interior_Model(ilayer).k_latlon)
        l_max = Interior_Model(ilayer).k_latlon.lmax;
    end
end

% Construct the viscosity and shear modulus in the fourier domain
for ilayer=2:Numerics.Nlayers
    if Interior_Model(ilayer).elastic==0 % viscoelastic model
        if Interior_Model(ilayer).uniform==0
            % build mu
            mu_stokes.clm=zeros(2*l_max,2*l_max);
            mu_stokes.slm=zeros(2*l_max,2*l_max);
            mu_stokes.lmax=2*l_max-1;
            for i=1:size(Interior_Model(ilayer).mu_variable,1)
                % Convert to real spherical harmonics
                if Interior_Model(ilayer).mu_variable(i,2)>0
                    nR=Interior_Model(ilayer).mu_variable(i,1);
                    mR=Interior_Model(ilayer).mu_variable(i,2);
                    amplitude=2/sqrt(2)*Interior_Model(ilayer).mu_variable(i,3);
                    mu_stokes.clm(nR+1,mR+1)=amplitude;
                elseif Interior_Model(ilayer).mu_variable(i,2)<0               
                else
                    nR=Interior_Model(ilayer).mu_variable(i,1);
                    mR=Interior_Model(ilayer).mu_variable(i,2);
                    amplitude=Interior_Model(ilayer).mu_variable(i,3);
                    mu_stokes.clm(nR+1,mR+1)=amplitude;
                end     
            end

            % Check whether direct insertion of a lat-lon grid is used
            if isfield(Interior_Model(ilayer),"mu_latlon")
                % Check if the element is empty
                if ~isempty(Interior_Model(ilayer).mu_latlon)
                    [mu_zlonlat] = Interior_Model(ilayer).mu_latlon;
                else
                    [mu_zlonlat] = SPH_LatLon(mu_stokes);
                end
            else
                [mu_zlonlat] = SPH_LatLon(mu_stokes);
            end

            % build eta 
            eta_stokes.clm=zeros(2*l_max,2*l_max);
            eta_stokes.slm=zeros(2*l_max,2*l_max);
            eta_stokes.lmax=2*l_max-1;
            for i=1:size(Interior_Model(ilayer).eta_variable,1)
                % Convert to real spherical harmonics
                if Interior_Model(ilayer).eta_variable(i,2)>0
                    nR=Interior_Model(ilayer).eta_variable(i,1);
                    mR=Interior_Model(ilayer).eta_variable(i,2);
                    amplitude=2/sqrt(2)*Interior_Model(ilayer).eta_variable(i,3);
                    eta_stokes.clm(nR+1,mR+1)=amplitude;
                elseif Interior_Model(ilayer).eta_variable(i,2)<0               
                else
                    nR=Interior_Model(ilayer).eta_variable(i,1);
                    mR=Interior_Model(ilayer).eta_variable(i,2);
                    amplitude=Interior_Model(ilayer).eta_variable(i,3);
                    eta_stokes.clm(nR+1,mR+1)=amplitude;
                end
            end

            % Check whether direct insertion of a lat-lon grid is used
            if isfield(Interior_Model(ilayer),"eta_latlon")
                % Check if the element is empty
                if ~isempty(Interior_Model(ilayer).eta_latlon)
                    [eta_zlonlat] = Interior_Model(ilayer).eta_latlon;
                else
                    [eta_zlonlat] = SPH_LatLon(eta_stokes);
                end
            else
                [eta_zlonlat] = SPH_LatLon(eta_stokes);
            end

            % Maxwell time
            MaxTime_zlonlat.lon=mu_zlonlat.lon;
            MaxTime_zlonlat.lat=mu_zlonlat.lat;
            MaxTime_zlonlat.lmax=mu_zlonlat.lmax;
            
            % When the viscosity and/or shear modulus fields are directly
            % given, the procedure is slightly different as it is assumed
            % that what is given is the entire field not just the deviation
            % from the mean.

            % Check whether both mu and eta are given as field or just one
            if isfield(Interior_Model(ilayer),"eta_latlon") && isfield(Interior_Model(ilayer),"mu_latlon")
                % Check whether one of the two (or both) are empty fields
                if ~isempty(Interior_Model(ilayer).eta_latlon) && ~isempty(Interior_Model(ilayer).mu_latlon)
                    MaxTime_zlonlat.z = eta_zlonlat.z./mu_zlonlat.z;
                elseif ~isempty(Interior_Model(ilayer).eta_latlon)
                    MaxTime_zlonlat.z = eta_zlonlat.z./(1+mu_zlonlat.z);
                elseif ~isempty(Interior_Model(ilayer).mu_latlon)
                    MaxTime_zlonlat.z = (1+eta_zlonlat.z)./mu_zlonlat.z;
                else
                    MaxTime_zlonlat.z = (1+eta_zlonlat.z)./(1+mu_zlonlat.z);
                end
            elseif isfield(Interior_Model(ilayer),"eta_latlon")
                if ~isempty(Interior_Model(ilayer).eta_latlon)
                    MaxTime_zlonlat.z = eta_zlonlat.z./(1+mu_zlonlat.z);
                else
                    MaxTime_zlonlat.z = (1+eta_zlonlat.z)./(1+mu_zlonlat.z);
                end
            elseif isfield(Interior_Model(ilayer),"mu_latlon")
                if ~isempty(Interior_Model(ilayer).mu_latlon)
                    MaxTime_zlonlat.z = (1+eta_zlonlat.z)./mu_zlonlat.z;
                else
                    MaxTime_zlonlat.z = (1+eta_zlonlat.z)./(1+mu_zlonlat.z);
                end
            else
                MaxTime_zlonlat.z = (1+eta_zlonlat.z)./(1+mu_zlonlat.z);
            end

            % Complex Shear modulus
            Cmu_zlonlat.lon=mu_zlonlat.lon;
            Cmu_zlonlat.lat=mu_zlonlat.lat;
            Cmu_zlonlat.lmax=mu_zlonlat.lmax;
            if isfield(Interior_Model(ilayer),"mu_latlon")
                % Check if the element is empty
                if ~isempty(Interior_Model(ilayer).mu_latlon)
                    Cmu_zlonlat.z = Interior_Model(ilayer).mu*(mu_zlonlat.z)./ ...
                            (1-1i./MaxTime_zlonlat.z/Interior_Model(ilayer).MaxTime);
                else
                    Cmu_zlonlat.z = Interior_Model(ilayer).mu*(1+mu_zlonlat.z)./ ...
                            (1-1i./MaxTime_zlonlat.z/Interior_Model(ilayer).MaxTime);
                end
            else
                Cmu_zlonlat.z = Interior_Model(ilayer).mu*(1+mu_zlonlat.z)./ ...
                            (1-1i./MaxTime_zlonlat.z/Interior_Model(ilayer).MaxTime);
            end

            % find spherical harmonics coefficients
            % real part 
            muR_zlonlat.lon=mu_zlonlat.lon;
            muR_zlonlat.lat=mu_zlonlat.lat;
            muR_zlonlat.lmax=mu_zlonlat.lmax;
            muR_zlonlat.z=real(Cmu_zlonlat.z);
            [muR_SPH]=LatLon_SPH(muR_zlonlat);

            %imaginary part
            muI_zlonlat.lon=mu_zlonlat.lon;
            muI_zlonlat.lat=mu_zlonlat.lat;
            muI_zlonlat.lmax=mu_zlonlat.lmax;
            muI_zlonlat.z=imag(Cmu_zlonlat.z);
            [muI_SPH]=LatLon_SPH(muI_zlonlat);

            % set components in the right format for the code
            vec_n=0:1:2*l_max-1;
            vec_m=-2*l_max+1:1:2*l_max-1;
            [m_g, n_g]=meshgrid(vec_m,vec_n);
            muR_aux=zeros(size(m_g));
            muI_aux=zeros(size(m_g));
            for i=1:length(vec_n) %degree 
                for j=1:vec_n(i)+1 %order
                    if vec_n(j)==0
                        indexN=find(vec_n(j)==m_g & vec_n(i)==n_g);
                        muR_aux(indexN)=muR_aux(indexN)+muR_SPH.clm(i,j);
                        muI_aux(indexN)=muI_aux(indexN)+muI_SPH.clm(i,j);
                    else
                        indexN1=find(vec_n(j)==m_g & vec_n(i)==n_g); %positive m
                        indexN2=find(-vec_n(j)==m_g & vec_n(i)==n_g); %negative m
                        % REAL mu
                        %positive m
                        muR_aux(indexN1)=muR_aux(indexN1)+1/sqrt(2)*muR_SPH.clm(i,j);
                        muR_aux(indexN1)=muR_aux(indexN1)-(-1)^vec_n(j)*1i/sqrt(2)*muR_SPH.slm(i,j)*sine_component_factor;
                        %negative m
                        muR_aux(indexN2)=muR_aux(indexN2)+(-1)^vec_n(j)/sqrt(2)*muR_SPH.clm(i,j);
                        muR_aux(indexN2)=muR_aux(indexN2)+1i/sqrt(2)*muR_SPH.slm(i,j)*sine_component_factor;
                        % IMAGINARY mu 
                        muI_aux(indexN1)=muI_aux(indexN1)+1/sqrt(2)*muI_SPH.clm(i,j);
                        muI_aux(indexN1)=muI_aux(indexN1)-(-1)^vec_n(j)*1i/sqrt(2)*muI_SPH.slm(i,j)*sine_component_factor;
                        %negative m
                        muI_aux(indexN2)=muI_aux(indexN2)+(-1)^vec_n(j)/sqrt(2)*muI_SPH.clm(i,j);
                        muI_aux(indexN2)=muI_aux(indexN2)+1i/sqrt(2)*muI_SPH.slm(i,j)*sine_component_factor;
                    end
                end
            end

            % 0,0 component
            index0=find(n_g==0 & m_g==0);
            mu00R=muR_aux(index0);
            mu00I=muI_aux(index0);
            mu00=mu00R+1i*mu00I; % Will be normalized later
            Interior_Model(ilayer).muC=mu00;

            % Find components that are bigger than some quantity, this is an approximation 
            muI_aux2=log10(abs(muI_aux/muI_SPH.clm(1,1)));
            muR_aux2=log10(abs(muR_aux/muR_SPH.clm(1,1)));
            muI_aux2(index0)=log10(0);
            muR_aux2(index0)=log10(0); 
            muI_max=max(muI_aux2(:));
            muR_max=max(muR_aux2(:));
            mu_max=max([muR_max muI_max]);
            non_zero_indexes=find((muI_aux2-mu_max>=-Numerics.rheology_cutoff | muR_aux2-mu_max>=-Numerics.rheology_cutoff) & ...
                                   mu_max > Numerics.minimum_rheology_value);
            k=1;
            % Check whether any rheology modes are big enough to be usable
            if ~isempty(non_zero_indexes)
                for i=1:length(non_zero_indexes)
                    if n_g(non_zero_indexes(i))>0
                        Interior_Model(ilayer).rheology_variable(k,1)=n_g(non_zero_indexes(i)); 
                        Interior_Model(ilayer).rheology_variable(k,2)=m_g(non_zero_indexes(i)); 
                        Interior_Model(ilayer).rheology_variable(k,4)=(muR_aux(non_zero_indexes(i))+1i*muI_aux(non_zero_indexes(i)));
                        k=k+1;
                    end
                end
                Interior_Model(ilayer).mu00R=mu00R;
    
                % PLOT THE RHEOLOGY STRUCTURE TO CHECK
                % There is an internal check whether it will plot or not.
                plot_rheology_map(plot_rheology_spectra,mu_zlonlat,eta_zlonlat, ...
                                  MaxTime_zlonlat,Interior_Model(ilayer),Cmu_zlonlat, ...
                                  m_g,n_g,muI_aux,muI_SPH,muR_aux,muR_SPH, ...
                                  mu00I,non_zero_indexes)

            % If there are no modes big enough the layer becomes uniform
            else
                Interior_Model(ilayer).uniform=1;
                muC=Interior_Model(ilayer).mu*(1-1i/Interior_Model(ilayer).MaxTime)^(-1);
                mu00R=real(muC);
                mu00I=imag(muC);
                mu00=mu00R+1i*mu00I;
                Interior_Model(ilayer).muC=mu00;
                Interior_Model(ilayer).mu00R=mu00R;
                Interior_Model(ilayer).rheology_variable=[0 0 0 0];
            end
            
        % Uniform model
        else
            muC=Interior_Model(ilayer).mu*(1-1i/Interior_Model(ilayer).MaxTime)^(-1);
            mu00R=real(muC);
            mu00I=imag(muC);
            mu00=mu00R+1i*mu00I;
            Interior_Model(ilayer).muC=mu00;
            Interior_Model(ilayer).mu00R=mu00R;
            Interior_Model(ilayer).rheology_variable=[0 0 0 0];
        end
    else %elastic
        if Interior_Model(ilayer).uniform==1
            Interior_Model(ilayer).rheology_variable=[0 0 0 0];
            Interior_Model(ilayer).muC=Interior_Model(ilayer).mu;
            Interior_Model(ilayer).mu00R=Interior_Model(ilayer).mu;
        else
            Interior_Model(ilayer).rheology_variable(:,[1 2])=Interior_Model(ilayer).mu_variable(:,[1 2]);
            Interior_Model(ilayer).rheology_variable(:,3)=0;
            Interior_Model(ilayer).rheology_variable(:,4)=Interior_Model(ilayer).mu_variable(:,3);
            Interior_Model(ilayer).muC=Interior_Model(ilayer).mu;
            Interior_Model(ilayer).mu00R=Interior_Model(ilayer).mu;
        end
    end
end


%% Re-normalize if asked for
if calculate_G
    for ilayer=1:Numerics.Nlayers
        if ilayer > 1
            % Calculate Lambda
            Interior_Model(ilayer).lambda=Interior_Model(ilayer).Ks - 2/3*Interior_Model(ilayer).muC;

            % Redundancy
            if Interior_Model(ilayer).uniform==1
                Interior_Model(ilayer).rheology_variable=[0 0 0 0];
            end
            
            % Fill an array with all the modes that are active/used
            if ilayer == 2
                rheo_degree_orders = Interior_Model(ilayer).rheology_variable(:,1:2);
            else
                rheo_degree_orders = [rheo_degree_orders ; Interior_Model(ilayer).rheology_variable(:,1:2)];
            end
        end
    end
else % Legacy normalization only used for comparisons to older data
    for ilayer=1:Numerics.Nlayers
        % Re-normalize Gg only using the values of the last layer
        Interior_Model(ilayer).Gg = Interior_Model(ilayer).Gg/Interior_Model(end).mu00R;
        if ilayer > 1
            % Re-normalize the 00 element of the rheology
            Interior_Model(ilayer).muC = Interior_Model(ilayer).muC/Interior_Model(end).mu00R;
            
            % Calculate Lambda
            Interior_Model(ilayer).lambda=Interior_Model(ilayer).Ks/Interior_Model(end).mu00R - 2/3*Interior_Model(ilayer).muC;
            if Interior_Model(ilayer).uniform==1
                Interior_Model(ilayer).rheology_variable=[0 0 0 0];
            end
            
            % Re-normalize the rheology_variable
            Interior_Model(ilayer).rheology_variable(:,4) = Interior_Model(ilayer).rheology_variable(:,4)/...
                                                                Interior_Model(end).mu00R;
            
            % Fill an array with all the modes that are active/used
            if ilayer == 2
                rheo_degree_orders = Interior_Model(ilayer).rheology_variable(:,1:2);
            else
                rheo_degree_orders = [rheo_degree_orders ; Interior_Model(ilayer).rheology_variable(:,1:2)];
            end
        end
    end
end

% Next we make sure that the size of the rheology_variable is the same for all
% layers. Modes that are not part of the original response are set to zero

% Remove duplicates from rheo_degree_orders
rheo_degree_orders = unique(rheo_degree_orders,'rows');

for ilayer=2:Numerics.Nlayers
    if (Interior_Model(ilayer).uniform==0) & (Interior_Model(ilayer).elastic==0)
        % Matrix to be filled, will replace Interior_Model.rheology_variable
        temp_rheology_var = zeros(length(rheo_degree_orders(:,1)),4);
    
        % Fill with all the active/used degrees and orders
        temp_rheology_var(:,1:2) = rheo_degree_orders;
    
        % Loop over all active modes and test whether a layer has a value for 
        % that mode. Fill the temp_rheology_var with the amplitudes
        for i=1:length(temp_rheology_var(:,1))
            n_rheo_tot = temp_rheology_var(i,1);
            m_rheo_tot = temp_rheology_var(i,2);
            ind_mode = find(Interior_Model(ilayer).rheology_variable(:,1)==n_rheo_tot & ...
                           Interior_Model(ilayer).rheology_variable(:,2)==m_rheo_tot);
            if ~isempty(ind_mode)
                temp_rheology_var(i,4) = Interior_Model(ilayer).rheology_variable(ind_mode,4);
            end
        end
    
        % Reset the rheology_variable such that all layers are the same size
        Interior_Model(ilayer).rheology_variable = temp_rheology_var;
    end
end

end
