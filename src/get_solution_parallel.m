%% OBTAIN_SOLUTION
%AUTHOR: M. Rovira-Navarro - A. Veenstra
%USE: Obtain solution using propagator method
%MODIFICATION: Computes the propagation matrix using parallel computing

%% INPUT 
    % Interior_Model: Structure containing the interior model information
        %Interior_Model.R0: radius
            %R0(1) radius core
            %R0(2) surface radius
        %Interior_Model.rho0: layers density        
            %rho0(1) density of liquid core 
            %rho0(2) density of rocky part 
        %Interior_Model.mu0: shear modulus of the mantle
        %Interior_Model.eta0: shear modulus of the mantle
        %Interior_Model.Ks0: bulk modulus of the mantle
        % Interio_Model.mu_variable: shear modulus variations
            %mu_variable(:,1): degree of variation 
            %mu_variable(:,2): order of variation 
            %mu_variable(:,3): amplitude of the variation (mu_l^m/mu^0_0)
        % Interio_Model.K_variable: first bulk modulus 
            %K_variable(:,1): degree of variation 
            %K_variable(:,2): order of variation 
            %K_variable(:,3):  amplitude of bulk modulus variations (K_l^m/K^0_0)
        % Interio_Model.eta_variable: first bulk modulus 
            %eta_variable(:,1): degree of variation 
            %eta_variable(:,2): order of variation 
            %eta_variable(:,3):  amplitude of viscosity variations (eta_l^m/eta^0_0)    
        % Interior_Model.rheology_variable: rheology variable (assigned inside the code)
            %rheology_variable(:,1): degree of variation 
            %rheology_variable(:,2): order of variation 
            %rheology_variable(:,3):  amplitude of bulk modulus variations (K_l^m/eta^0_0)  
            %rheology_variable(:,4):  amplitude of complex shear modulus variations (mu_l^m/mu^0_0) 
        % Interior_Model.muC: complex shear modulus (assigned inside the code)
        
    % Forcing: Structure containing forcing information
        %Forcing.Td: forcing period
        %Forcing.n: degree of the forcing 
        %Forcing.m: order of the forcing 
        
    % Numerics: Structure containing the numerical information
        %Numerics.Nr: number of radial points
        %Numerics.perturbation_order: maximum order of the perturbation. Default 2
        %Numerics.solution_cutoff: cut off degree (if specified instead of perturbation_order)
        %Numerics.rheology_cutoff: cut off degree for rheology. Default 30
        %Numerics.load_couplings: load coupling coefficients (1) or not (0). Default 1. 
     
    % Couplings
        % Couplings.n_s: degrees for which solution needs to be solved
        % Couplings.m_s: order for which solution needs to be solved
        % Coup: Matrix containing the coupling coefficients
            %%Coupl_re(ic,ia,:,ir)
                %n_s(ic) is the degree of equation 
                %m_s(ic) is the order of the equation 
                %n_s(ia) is the degree of contribution
                %m_s(ia) is the order of contribution
                %Coupl_re(ic,ia,27,ir): 1 if any coefficient is not-zero 
                %ir indicates the rheology coupling considered
    % verbose (1) print information 
    % out_file_name: name output file (0) if no printing
%% OUTPUT 
    % y_sol: solution matrix
        % y_sol(radial_point,1,mode)): r
        % y_sol(radial_point,2,mode)): U_n^m, radial displacement
        % y_sol(radial_point,3,mode): V_n^m, tangential displacement
        % y_sol(radial_point,4,mode): W_n^m, toroidal displacement
        % y_sol(radial_point,5,mode): R_n^m, radial stress
        % y_sol(radial_point,6,mode): S_n^m, tangential stress
        % y_sol(radial_point,7,mode): W_n^m, toroidal stress
        % y_sol(radial_point,8,mode): \phi_n^m, gravitational potential
        % y_sol(radial_point,9,mode): \dot\phi_n^m, gradient grav potential
        % y_sol(radial_point,10,mode): u_{n,n-1}
        % y_sol(radial_point,11,mode): u_{n,n}
        % y_sol(radial_point,12,mode): u_{n,n+1}
        % y_sol(radial_point,13,mode): \sigma_{n,n,0}       
        % y_sol(radial_point,14,mode): \sigma_{n,n-2,2}
        % y_sol(radial_point,15,mode): \sigma_{n,n-1,2}
        % y_sol(radial_point,16,mode): \sigma_{n,n,2}
        % y_sol(radial_point,17,mode): \sigma_{n,n+1,2}
        % y_sol(radial_point,18,mode): \sigma_{n,n+2,2}      
        % y_sol(radial_point,19,mode): \epsilon_{n,n,0}
        % y_sol(radial_point,20,mode): \epsilon_{n,n-2,2}
        % y_sol(radial_point,21,mode): \epsilon_{n,n-1,2}
        % y_sol(radial_point,22,mode): \epsilon_{n,n,2}
        % y_sol(radial_point,23,mode): \epsilon{n,n+1,2}
        % y_sol(radial_point,24,mode): \epsilon_{n,n+2,2}
%% FUNCTION
function [y_sol] = get_solution_parallel(Interior_Model,Forcing,Numerics,Couplings,verbose,out_file_name)
%% OBTAIN MATRICES THAT ONLY NEED TO BE COMPUTED ONCE
    % \sigma=A1\dot u+A2/r u 
    % \epsilon=A14\dot u+A15/r u
    % U=A3u 
    % A13\Sigma=A4\sigma
    % A13\dot\Sigma=A5*\sigma/r+A6\dot{U}+g/r*A71*U+dg*A72*U+A81*\Phi+A82*\Phi/r
    % A9\dot\Phi=A100\Phi+A101/r*\Phi+A102/r*\Phi+A11/r*U+A12\dot{U}
    % with:
    % \sigma=[\sigma_{n,n,0}, \sigma_{n,n-2,2}, \sigma_{n,n-1,2}, \sigma_{n,n,2}, \sigma_{n,n+1,2}, \sigma_{n,n+2,0}]
    % \epsilon=[\epsilon_{n,n,0}, \epsilon_{n,n-2,2}, \epsilon_{n,n-1,2}, \epsilon_{n,n,2}, \epsilon_{n,n+1,2}, \epsilon_{n,n+2,0}]
    % \Sigma=[R_{n,m} S_{n,m} W_{n,m}]
    % u=[u_{n-1} u_n u_{n+1}]
    % y=[U_n^m, V_n^m, W_n^m, R_n^m, S_n^m, W_n^m, \phi_n^m, \dot\phi_n^m]

% Extract Coupling degrees and modes
n_s = Couplings.n_s;
m_s = Couplings.m_s;

ind = find(n_s==0, 1);
if isempty(ind)==1
    deg0=0;
else
    deg0=1; 
end

% Initialise the solution matrices
Nmodes = length(n_s);
y = zeros(8*Nmodes,8*Nmodes,Numerics.Nr+1);
AProp_arr = zeros(8*Nmodes,8*Nmodes,Numerics.Nr+1);
AProp_arr_aux = zeros(3*Nmodes,8*Nmodes,Numerics.Nr+1);
A1_arr=zeros(6*Nmodes,3*Nmodes,Numerics.Nlayers-1);
A2_arr=zeros(6*Nmodes,3*Nmodes,Numerics.Nlayers-1);
one_arr = zeros(8*Nmodes,8*Nmodes,1);
for i=1:8*Nmodes
    y(i,i,1)=1; %set Cs. 
    one_arr(i,i) = 1; % Identity matrix
end

% Set the indice of the first boundary layer if applicable
blayer_ind = 0; % Set default value
iblayer = 1; % First indices of Numerics.BCindices
k_corr = 1; % Indices correction %0
if ~isempty(Numerics.BCindices)
    % Add 1 to make sure that the Boundary is at the previous node
    blayer_ind = Numerics.BCindices(iblayer) + 1;
end

% Start counter of layer number
ilayer = 2; % The core does not contribute to tidal heating

% Calculate some necessary values
Delta_r = (Interior_Model(ilayer).R - Interior_Model(1).R)/(Numerics.Nrlayer(ilayer)); % Changes with layer
Rc = Interior_Model(1).R; % Radius of the core (is fixed)
Rin = Interior_Model(1).R; % Radius of the inner sphere (can change)
rhoC = Interior_Model(1).rho; % Density of the core
Gg = Interior_Model(end).Gg; 
Mc = 4/3*pi*rhoC*Rc^3; % Mass of the core
Min = 4/3*pi*rhoC*Rc^3; % Mass of the inner sphere (can change)
gc = Gg*Mc/Rc^2; % g at the CMB
rhoK = Interior_Model(ilayer).rho; % Density of the current layer

% Initialise the radial distance vector
r = zeros(1,Numerics.Nr+1);
r(1) = Rc;

% Prime the continuity correction matrix
cont_condition = zeros(8*Nmodes,8*Nmodes,1);

% Constants used in RK integration (see Cash-Karp method e.g.,https://doi.org/10.1063/1.4823060)
AA2=1/5; AA3=3/10; AA4=3/5; AA5=1; AA6=7/8;
B21=1/5;
B31=3/40; B32=9/40;
B41=3/10; B42=-9/10; B43=6/5;
B51=-11/54; B52=5/2; B53=-70/27; B54=35/27;
B61=1631/55296; B62=175/512; B63=575/13824; B64=44275/110592; B65=253/4096; 
AC1=37/378; AC2=0; AC3=250/621; AC4=125/594; AC5=0; AC6=512/1771; 

% Fill matrices that do not depend on the rheology
[A14,A15]=get_A14A15(Couplings); 
A3=get_A3(Couplings); 
A4=get_A4(Couplings);
A5=get_A5(Couplings);
A3_inv=inv(A3);
if deg0==1
    A3_inv([1,2],:)=0;
end

% Initialize matrices needed for paralellization
Delta_r_arr = zeros(1,Numerics.Nr+1);
ilayer_arr = zeros(1,Numerics.Nr+1);
Min_arr = zeros(1,Numerics.Nr+1);
Rin_arr = zeros(1,Numerics.Nr+1);
rhoK_arr = zeros(1,Numerics.Nr+1);
rK_arr = zeros(6,Numerics.Nr+1);
gK_arr = zeros(6,Numerics.Nr+1);
dgK_arr = zeros(6,Numerics.Nr+1);

% Fill first indice of said matrices
Delta_r_arr(1) = Delta_r;
Min_arr(1) = 4/3*pi*rhoC*Rc^3;
Rin_arr(1) = Interior_Model(1).R;
rhoK_arr(1) = Interior_Model(ilayer).rho;
ilayer_arr(1) = ilayer;
rK_arr(:,1) = Rin;
gK_arr(:,1) = Gg*(Min)/Rin^2;
dgK_arr(:,1) = Gg*(2*(4/3*pi*rhoK*Rin^3 - Min)/Rin^3 + 4/3*pi*rhoK);

% Compute matrices that depend on rheology/layer
[A2,A1] = get_A1A2(Interior_Model(ilayer),Couplings); 
[A13, A6, A71, A72, A81, A82, A9, A100, A101, A102, A11, A12] = get_others(Couplings,Interior_Model(ilayer));

% Caculate propagation matrix for auxilliary properties at index 1
Aprop = get_Aprop(Rc,gK_arr(1,1),dgK_arr(1,1),Nmodes,A1,A2,A3_inv,A4,A5,A6,A71,A72,A81,A82,A9,A100,A101,A102,A11,A12,A13,deg0,Gg);
AProp_arr_aux(:,:,1) = Aprop(1:3*Nmodes,:);

%% PRECOMPUTE THE PROPAGATION MATRIX
if verbose==1
    tic
end

% Pre-compute some arrays to allow parallelization later
for k=2:Numerics.Nr+1
    % Check whether the integration is at a new layer
    if k == blayer_ind 
        % Go to next layer
        ilayer = ilayer + 1;
        
        % If there is another layer reprime blayer_ind
        if iblayer < length(Numerics.BCindices)
            iblayer = iblayer + 1;
            blayer_ind = Numerics.BCindices(iblayer) + 1;
        else
            % Update regardless to keep consistency in last layer
            iblayer = iblayer + 1;
        end

        % Update density 
        rhoK = Interior_Model(ilayer).rho;

        % Update delta r
        Delta_r = (Interior_Model(ilayer).R - Interior_Model(ilayer-1).R)/(Numerics.Nrlayer(ilayer));

        % Update radial points vector before new Delta_r and Rin
        r(k) = Rin + (k-k_corr)*Delta_r;

        % Update mass and radius of the encapsulated sphere
        Min = Min + 4/3*pi*Interior_Model(ilayer-1).rho*(Interior_Model(ilayer-1).R^3 - Interior_Model(ilayer-2).R^3);
        Rin = Interior_Model(ilayer-1).R;

        % Update indices correction, needed to 
        k_corr = Numerics.BCindices(iblayer-1);

    else
        % Update radial points vector
        r(k) = Rin + (k-k_corr)*Delta_r;
    end
    
    % Save these values in an array for later
    Min_arr(k) = Min;
    Rin_arr(k) = Rin;
    rhoK_arr(k) = rhoK;
    ilayer_arr(k) = ilayer;
    Delta_r_arr(k) = Delta_r;
    
    %%%%%%%%%%% Step 1
    rK = Rin + ((k-k_corr)-1)*Delta_r;
    rK_arr(1,k) = rK;
    gK_arr(1,k) = Gg*(Min + 4/3*pi*rhoK*(rK^3 - Rin^3))/rK^2;
    dgK_arr(1,k) = Gg*(2*(4/3*pi*rhoK*Rin^3 - Min)/rK^3 + 4/3*pi*rhoK);

    %%%%%%%%%%% Step 2
    rK = Rin + ((k-k_corr)-1)*Delta_r + AA2*Delta_r; 
    rK_arr(2,k) = rK; 
    gK_arr(2,k) = Gg*(Min + 4/3*pi*rhoK*(rK^3 - Rin^3))/rK^2;
    dgK_arr(2,k) = Gg*(2*(4/3*pi*rhoK*Rin^3 - Min)/rK^3 + 4/3*pi*rhoK);

    %%%%%%%%%%% Step 3
    rK = Rin + ((k-k_corr)-1)*Delta_r + AA3*Delta_r; 
    rK_arr(3,k) = rK; 
    gK_arr(3,k) = Gg*(Min + 4/3*pi*rhoK*(rK^3 - Rin^3))/rK^2;
    dgK_arr(3,k) = Gg*(2*(4/3*pi*rhoK*Rin^3 - Min)/rK^3 + 4/3*pi*rhoK);

    %%%%%%%%%%% Step 4
    rK = Rin + ((k-k_corr)-1)*Delta_r + AA4*Delta_r; 
    rK_arr(4,k) = rK; 
    gK_arr(4,k) = Gg*(Min + 4/3*pi*rhoK*(rK^3 - Rin^3))/rK^2;
    dgK_arr(4,k) = Gg*(2*(4/3*pi*rhoK*Rin^3 - Min)/rK^3 + 4/3*pi*rhoK);

    %%%%%%%%%%% Step 5
    rK = Rin + ((k-k_corr)-1)*Delta_r + AA5*Delta_r;
    rK_arr(5,k) = rK; 
    gK_arr(5,k) = Gg*(Min + 4/3*pi*rhoK*(rK^3 - Rin^3))/rK^2;
    dgK_arr(5,k) = Gg*(2*(4/3*pi*rhoK*Rin^3 - Min)/rK^3 + 4/3*pi*rhoK);

    %%%%%%%%%%% Step 6
    rK = Rin + ((k-k_corr)-1)*Delta_r + AA6*Delta_r; 
    rK_arr(6,k) = rK; 
    gK_arr(6,k) = Gg*(Min + 4/3*pi*rhoK*(rK^3 - Rin^3))/rK^2;
    dgK_arr(6,k) = Gg*(2*(4/3*pi*rhoK*Rin^3 - Min)/rK^3 + 4/3*pi*rhoK);
end

% Calculate the propagation matrix 
for ii=2:Numerics.Nlayers
    if ii > 2
        % Recalculate the parts of the propagation matrix that depend on the rheology
        [A2,A1] = get_A1A2(Interior_Model(ii),Couplings); 
        [A13, A6, A71, A72, A81, A82, A9, A100, A101, A102, A11, A12] = get_others(Couplings,Interior_Model(ii));
    end
    inds = find(ilayer_arr == ii);
    start_ind = inds(1);
    end_ind = inds(end);
    if ii == 2
        start_ind = start_ind+1;
    end
    A1_arr(:,:,ii-1) = A1;
    A2_arr(:,:,ii-1) = A2;

    parfor k = start_ind:end_ind
        %%%%%%%%%%% Step 1
        Aprop = get_Aprop(rK_arr(1,k),gK_arr(1,k),dgK_arr(1,k),Nmodes,A1,A2,A3_inv,A4,A5,A6,A71,A72,A81,A82,A9,A100,A101,A102,A11,A12,A13,deg0,Gg);
        k1 = Delta_r_arr(k)*Aprop*one_arr;
        %%%%%%%%%%% Step 2
        Aprop = get_Aprop(rK_arr(2,k),gK_arr(2,k),dgK_arr(2,k),Nmodes,A1,A2,A3_inv,A4,A5,A6,A71,A72,A81,A82,A9,A100,A101,A102,A11,A12,A13,deg0,Gg);
        k2 = Delta_r_arr(k)*Aprop*(one_arr + B21*k1);
        %%%%%%%%%%% Step 3
        Aprop = get_Aprop(rK_arr(3,k),gK_arr(3,k),dgK_arr(3,k),Nmodes,A1,A2,A3_inv,A4,A5,A6,A71,A72,A81,A82,A9,A100,A101,A102,A11,A12,A13,deg0,Gg);
        k3 = Delta_r_arr(k)*Aprop*(one_arr + B31*k1 + B32*k2);
        %%%%%%%%%%% Step 4
        Aprop = get_Aprop(rK_arr(4,k),gK_arr(4,k),dgK_arr(4,k),Nmodes,A1,A2,A3_inv,A4,A5,A6,A71,A72,A81,A82,A9,A100,A101,A102,A11,A12,A13,deg0,Gg);
        k4 = Delta_r_arr(k)*Aprop*(one_arr + B41*k1 + B42*k2 + B43*k3);
        %%%%%%%%%%% Step 5
        Aprop = get_Aprop(rK_arr(5,k),gK_arr(5,k),dgK_arr(5,k),Nmodes,A1,A2,A3_inv,A4,A5,A6,A71,A72,A81,A82,A9,A100,A101,A102,A11,A12,A13,deg0,Gg);
        k5 = Delta_r_arr(k)*Aprop*(one_arr + B51*k1 + B52*k2 + B53*k3 + B54*k4);
        %%%%%%%%%%% Step 6
        Aprop = get_Aprop(rK_arr(6,k),gK_arr(6,k),dgK_arr(6,k),Nmodes,A1,A2,A3_inv,A4,A5,A6,A71,A72,A81,A82,A9,A100,A101,A102,A11,A12,A13,deg0,Gg);
        k6 = Delta_r_arr(k)*Aprop*(one_arr + B61*k1 + B62*k2 + B63*k3 + B64*k4 + B65*k5);
    
        % Get next step
        AProp_arr(:,:,k) = AC1*k1 + AC3*k3 + AC4*k4 + AC6*k6;
    
        % Compute Aprop for auxilliary coefficients
        % Retrieve values
        rK = r(k); 
        Min = Min_arr(k);
        Rin = Rin_arr(k);
        rhoK = rhoK_arr(k);
    
        % Calculate gravity acceleration necessary
        gK = Gg*(Min + 4/3*pi*rhoK*(rK^3 - Rin^3))/rK^2;
        dgK = Gg*(2*(4/3*pi*rhoK*Rin^3 - Min)/rK^3 + 4/3*pi*rhoK);
    
        % Caculate propagation matrix
        Aprop = get_Aprop(rK,gK,dgK,Nmodes,A1,A2,A3_inv,A4,A5,A6,A71,A72,A81,A82,A9,A100,A101,A102,A11,A12,A13,deg0,Gg);
        AProp_arr_aux(:,:,k) = Aprop(1:3*Nmodes,:);
    end
end

% Clear some big matrices from memory
clear Couplings rK_arr gK_arr dgK_arr Delta_r_arr

if verbose==1
    disp(['Time Spent computing the propagation matrix ' num2str(toc) ' s'])
end

%% PROPAGATE THE SOLUTION
if verbose == 1
    tic
end

%%%%%%%%%%% Reset all auxiliary variables %%%%%%%%%%%%
% Set the indice of the first boundary layer if applicable
blayer_ind = 0; % Set default value
iblayer = 1; % First indices of Numerics.BCindices
if ~isempty(Numerics.BCindices)
    % Add 1 to make sure that the Boundary is at the previous node
    blayer_ind = Numerics.BCindices(iblayer) + 1;
end

% Start counter of layer number
ilayer = 2; % The core does not contribute to tidal heating

% Start the radial integration using the propagation matrix
for k=2:Numerics.Nr+1
    % Check whether the integration is at a new layer
    if k == blayer_ind 
        % Go to next layer
        ilayer = ilayer + 1;

        % Obtain density difference with previous layer (is positive)
        rho_diff = Interior_Model(ilayer).Delta_rho;
        
        % If there is another layer reprime blayer_ind
        if iblayer < length(Numerics.BCindices)
            iblayer = iblayer + 1;
            blayer_ind = Numerics.BCindices(iblayer) + 1;
        else
            % Update regardless to keep consistency in last layer
            iblayer = iblayer + 1;
        end

        % Apply the effect of the density discontinuity
        cont_condition(8*(0:(Nmodes-1))+8,:,1) = ...
            4*pi*Gg*rho_diff*y(8*(0:(Nmodes-1))+1,:,k-1);
        y_old = y(:,:,k-1) + cont_condition;
    else
        % Update the solution at the previous node
        y_old = y(:,:,k-1);
    end

    % get next step
    y(:,:,k) = y_old + AProp_arr(:,:,k)*y_old;
end

% Clear a big matrix from memory
clear AProp_arr y_old

% Reset some necessary values
rhoK = Interior_Model(end).rho; % Density of the end layer

if verbose==1
    disp(['Time Spent Integrating radially ' num2str(toc) ' s'])
end

%% REARRANGE SOLUTION VECTOR 
if verbose==1
    tic
end

index_arr = zeros(8*Nmodes,1);
for i=1:Nmodes
    index1 = 3*(i-1)+(1:3);
    index2 = 3*Nmodes+3*(i-1)+(1:3);
    index3 = 6*Nmodes+2*(i-1)+(1:2);
    index_arr(8*(i-1)+(1:8)) = [index1 index2 index3];
end

% Rearrange the matrix
y = y(index_arr,index_arr,:);

if verbose==1
    disp(['Time Spent rearranging the solution vector ' num2str(toc) ' s'])
end

%% ASSEMBLE MATRIX FOR INVERSION
if verbose==1
    tic
end
B=zeros(8*Nmodes,8*Nmodes);
B2=zeros(8*Nmodes,1);
if verbose==1
    disp(' ###########################################################################')
    disp(' Propagating Solution ...')
end
for i=1:Nmodes
    n=n_s(i);
    m=m_s(i);

    % Written such that the code handles icy moons correctly
    rho2 = Interior_Model(1).Delta_rho + Interior_Model(2).rho; 

    for j=1:8*Nmodes
    %Core-Mantle BC--------------
    % BC1 radial stress
    if n==-3
    B(8*(i-1)+1,j)=B(8*(i-1)+1,j)+y(8*(i-1)+4,j,1); %R=0
    else
    B(8*(i-1)+1,j)=B(8*(i-1)+1,j)+y(8*(i-1)+1,j,1); %U
    B(8*(i-1)+1,j)=B(8*(i-1)+1,j)-1/gc/rho2*y(8*(i-1)+4,j,1); %R 
    B(8*(i-1)+1,j)=B(8*(i-1)+1,j)+1/gc*y(8*(i-1)+7,j,1); %\Phi
    end
    % BC2 no tangential stress
    B(8*(i-1)+2,j)=B(8*(i-1)+2,j)+y(8*(i-1)+5,j,1);  % S 
    % BC3 no toroidal stress 
    if n==1
        B(8*(i-1)+3,j)=B(8*(i-1)+3,j)+y(8*(i-1)+3,j,1);   %W
    else
        B(8*(i-1)+3,j)=B(8*(i-1)+3,j)+y(8*(i-1)+6,j,1);   %T 
    end
    % BC4 potential stress
    if n==0
        B(8*(i-1)+4,j)=B(8*(i-1)+4,j)+y(8*(i-1)+8,j,end); %k=0
    else
        fac=4*pi*Gg/gc/rho2*(Interior_Model(2).rho - rho2);
        B(8*(i-1)+4,j)=B(8*(i-1)+4,j)-fac*y(8*(i-1)+4,j,1); %R
        B(8*(i-1)+4,j)=B(8*(i-1)+4,j)+(n/Rc + fac*rho2)*y(8*(i-1)+7,j,1); %\Phi
        B(8*(i-1)+4,j)= B(8*(i-1)+4,j)-y(8*(i-1)+8,j,1);  %\dot\Phi
    end
    % Surface BC------------------- 
    % BC5 no radial stress
    B(8*(i-1)+5,j)=B(8*(i-1)+5,j)+y(8*(i-1)+4,j,end); %R=0
    % BC6 no tangential stress
    if n==0
        B(8*(i-1)+6,j)=B(8*(i-1)+6,j)+y(8*(i-1)+2,j,end); %V=0
    else
        B(8*(i-1)+6,j)=B(8*(i-1)+6,j)+y(8*(i-1)+5,j,end); %S=0
    end
    % BC7 no toroidal stress 
    if n==0
        B(8*(i-1)+7,j)=B(8*(i-1)+7,j)+y(8*(i-1)+3,j,end); %W=0 
    else
        B(8*(i-1)+7,j)=B(8*(i-1)+7,j)+y(8*(i-1)+6,j,end); %T=0  
    end
    % BC8 potential stress 
    if n<2
        B(8*(i-1)+8,j)=B(8*(i-1)+8,j)+y(8*(i-1)+7,j,end); %k=0
    else
        B(8*(i-1)+8,j)=B(8*(i-1)+8,j)+4*pi*Gg*rhoK*y(8*(i-1)+1,j,end); %U
        B(8*(i-1)+8,j)=B(8*(i-1)+8,j)+(n+1)*y(8*(i-1)+7,j,end); %\Phi
        B(8*(i-1)+8,j)=B(8*(i-1)+8,j)+y(8*(i-1)+8,j,end); %\dot\Phi
    end
    end
    if n==Forcing.n && m==Forcing.m
        B2(8*(i-1)+8,1)=2*n+1;
    end
end

if verbose==1
    disp(['Time Spent assembling the boundary conditions matrix ' num2str(toc) ' s'])
end

%% OBTAIN COEFFICIENTS
if verbose==1
    tic
end
C = B\B2;
if verbose==1
    disp(['Time Spent calculating the coefficients ' num2str(toc) ' s'])
end

% Clear some big matrices from memory
clear B B2

%% ASSEMBLE SOLUTION
if verbose==1
    tic
end
% THIS STEP CAN BE VECTORIZED
y_sol = zeros(Numerics.Nr+1,8*Nmodes);
for i=1:Numerics.Nr+1
    y_sol(i,:)=transpose(y(:,:,i)*C);
end
if verbose==1
    disp(['Time Spent assembling the solution ' num2str(toc) ' s'])
end

% Clear a big matrix from memory
clear y

%% Check BC
if verbose==1
    if Nmodes > 60
        disp('Skipping display of boundary conditions because Nmodes > 60')
    else
        for i=1:Nmodes
            n=n_s(i);
            mode=8*(i-1); 
            disp(['Mode (n,m) (' num2str(n_s(i)) ',' num2str(m_s(i)) ')'])
            disp(['CMB:' num2str(y_sol(1,mode+1)+y_sol(1,mode+7)/gc-y_sol(1,mode+4)/rho2/gc) ', ' num2str(y_sol(1,mode+5)) ',' num2str(y_sol(1,mode+6)) ', ' num2str(-y_sol(1,mode+8)+(n/Rc+fac*rho2)*y_sol(1,mode+7)-fac*y_sol(1,mode+4))])
            disp(['Surface: ' num2str(y_sol(end,mode+4)) ', ' num2str(y_sol(end,mode+5)) ', ' num2str(y_sol(end,mode+6)) ', ' num2str(y_sol(end,mode+8)+(n+1)*y_sol(end,mode+7)+4*pi*Gg*rhoK*y_sol(end,mode+1))])
        end
    end
end

%% COMPUTE AUXILIARY VARIABLES
if verbose==1
    tic
end

ind1 = zeros(1,3*Nmodes);% i_U
ind2 = zeros(1,3*Nmodes);
ind3 = zeros(1,2*Nmodes);
for i=1:Nmodes
    ind1(3*(i-1)+(1:3)) = 8*(i-1)+(1:3);
    ind2(3*(i-1)+(1:3)) = 8*(i-1)+(4:6);
    ind3(2*(i-1)+(1:2)) = 8*(i-1)+(7:8);
end

U = y_sol(:,ind1);
sigma = zeros(Numerics.Nr+1,6*Nmodes);
epsilon = zeros(Numerics.Nr+1,6*Nmodes);

% Compute u 
u = transpose(A3_inv*transpose(U));

% Compue u_dot (numerically)
u_dot = zeros(Numerics.Nr+1,3*Nmodes);
u_aux = zeros(size(u));
u_aux(1:end-1,:) = u(1:end-1,:);
y_sol_temp = y_sol(:,[ind1 ind2 ind3]);

for ii=2:Numerics.Nlayers
    inds = find(ilayer_arr == ii);
    start_ind = inds(1);
    end_ind = inds(end);
    if ii == Numerics.Nlayers
        end_ind = end_ind - 1;
    end
    A1 = A1_arr(:,:,ii-1);
    A2 = A2_arr(:,:,ii-1);
    
    parfor k = start_ind:end_ind
        % Retrieve values
        rK = r(k); 
        
        % Compute values
        x_dot = AProp_arr_aux(:,:,k)*transpose(y_sol_temp(k,:));
        u_dot(k,:) = transpose(A3_inv*x_dot); 
        sigma(k,:) = transpose(A1*transpose(u_dot(k,:)) + A2*transpose(u_aux(k,:)) / rK);
    end
end
epsilon(:,:)=transpose(A14*transpose(u_dot)+A15*transpose(u_aux)./r);

% Clear some big matrices from memory
clear AProp_arr_aux y_sol_temp

%compue u_dot (numerically)
% u_dot=diff(u)/Delta_r; 
% r_av=(r(1:end-1)'+r(2:end)')/2;
% u_av=(u(1:end-1,:)+u(2:end,:))/2;
% sigma(2:end,:)=transpose(A1*transpose(u_dot)+A2*transpose(u_av)./r_av');
% epsilon(2:end,:)=transpose(A14*transpose(u_dot)+A15*transpose(u_av)./r_av');
if verbose==1
    disp(['Time Spent calculating the auxilliary variables ' num2str(toc) ' s'])
end

%% REARRANGE
if verbose==1
    tic
end

y_sol2 = y_sol;
y_sol = zeros(Numerics.Nr+1,24,Nmodes);

for i=1:Nmodes
    % r
    y_sol(:,1,i)=r;
    %U,V,W,R,S,T,\phi,\dot\phi
    y_sol(:,2:9,i)=y_sol2(:,8*(i-1)+(1:8));
    %u 
    y_sol(:,10:12,i)=u(:,3*(i-1)+(1:3));
    %\sigma
    y_sol(:,13:18,i)=sigma(:,6*(i-1)+(1:6));
    %\epsilon
    y_sol(:,19:24,i)=epsilon(:,6*(i-1)+(1:6));
end

if verbose==1
    disp(['Time Spent rearranging the solution vector ' num2str(toc) ' s'])
end

end

%% FUNCTION USED IN ROUTINE
%% get_A1A2
% \sigma=A1\dot u+A2/r u 
%\sigma=[\sigma_{n,n,0}, \sigma_{n,n-2,2}, \sigma_{n,n-1,2}, \sigma_{n,n,2}, \sigma_{n,n+1,2}, \sigma_{n,n+2,2},]
% u=[u_{n-1} u_n u_{n+1}]
% y=[U_n^m, V_n^m, W_n^m, R_n^m, S_n^m, W_n^m, \phi_n^m, \dot\phi_n^m]
function [A1,A2]=get_A1A2(Interior_Model,Couplings)
N=length(Couplings.n_s);
A1=zeros(6*N,3*N);
A2=zeros(6*N,3*N);
mu=Interior_Model.muC; 
lambda=Interior_Model.lambda;
variable_rheology=Interior_Model.rheology_variable; 
Nreo=size(variable_rheology,1);
Coup=Couplings.Coup; 
n_s=Couplings.n_s;
m_s=Couplings.m_s;
r=1; 
for i=1:N
    n=Couplings.n_s(i);
    ieq=i;
    if n>0
        %DIAGONAL TERMS-----------------------------   
        %\sigma_{n,n,0}------------------------------
        %u
        A1(6*(i-1)+1,3*(i-1)+1)=(3*lambda+2*mu)/sqrt(3)/sqrt(2*n+1)*sqrt(n)*(n-1);
        % A1(6*(i-1)+1,3*(i-1)+2)=0;
        A1(6*(i-1)+1,3*(i-1)+3)=(3*lambda+2*mu)/sqrt(3)/sqrt(2*n+1)*sqrt(n+1)*(n+2);  
        %\dot u
        A2(6*(i-1)+1,3*(i-1)+1)=-(3*lambda+2*mu)/sqrt(3)/sqrt(2*n+1)*sqrt(n);
        % A2(6*(i-1)+1,3*(i-1)+2)=0;
        A2(6*(i-1)+1,3*(i-1)+3)=+(3*lambda+2*mu)/sqrt(3)/sqrt(2*n+1)*sqrt(n+1);    
        %\sigma_{n,n-2,2}----------------------------
        %u
        A1(6*(i-1)+2,3*(i-1)+1)=2*mu*sqrt((n-1)/(2*n-1))*n;
        % A1(6*(i-1)+2,3*(i-1)+2)=0;
        % A1(6*(i-1)+2,3*(i-1)+3)=0;   
        %\dot u
        A2(6*(i-1)+2,3*(i-1)+1)=2*mu*sqrt((n-1)/(2*n-1));
        % A2(6*(i-1)+2,3*(i-1)+2)=0;
        % A2(6*(i-1)+2,3*(i-1)+3)=0;    
        %\sigma_{n,n-1,2}-----------------------------
        %u
        % A1(6*(i-1)+3,3*(i-1)+1)=0;
        A1(6*(i-1)+3,3*(i-1)+2)=2*mu/sqrt(2)*sqrt((n-1)/(2*n+1))*(n+1);
        % A1(6*(i-1)+3,3*(i-1)+3)=0;    
        %\dot u
        % A2(6*(i-1)+3,3*(i-1)+1)=0;
        A2(6*(i-1)+3,3*(i-1)+2)=2*mu/sqrt(2)*sqrt((n-1)/(2*n+1));
        % A2(6*(i-1)+3,3*(i-1)+3)=0;  
        %\sigma_{n,n,2}-------------------------------
        %u
        A1(6*(i-1)+4,3*(i-1)+1)=2*mu*sqrt((2*n+3)*(2*n+2)/12/(2*n-1)/(2*n+1))*(n-1);
        % A1(6*(i-1)+4,3*(i-1)+2)=0;
        A1(6*(i-1)+4,3*(i-1)+3)=2*mu*sqrt(n*(2*n-1)*(n+1)/3/(2*n+3)/(2*n+2)/(2*n+1))*(n+2);    
        %\dot u
        A2(6*(i-1)+4,3*(i-1)+1)=-2*mu*sqrt((2*n+3)*(2*n+2)/12/(2*n-1)/(2*n+1));
        % A2(6*(i-1)+4,3*(i-1)+2)=0;
        A2(6*(i-1)+4,3*(i-1)+3)=2*mu*sqrt(n*(2*n-1)*(n+1)/3/(2*n+3)/(2*n+2)/(2*n+1));   
        %\sigma_{n,n+1,2}-----------------------------
        %u
        % A1(6*(i-1)+5,3*(i-1)+1)=0;
        A1(6*(i-1)+5,3*(i-1)+2)=2*mu/sqrt(2)*sqrt((n+2)/(2*n+1))*n;
        % A1(6*(i-1)+5,3*(i-1)+3)=0;    
        %\dot u
        % A2(6*(i-1)+5,3*(i-1)+1)=0;
        A2(6*(i-1)+5,3*(i-1)+2)=-2*mu/sqrt(2)*sqrt((n+2)/(2*n+1));
        % A2(6*(i-1)+5,3*(i-1)+3)=0;
        %\sigma_{n,n+2,2}-----------------------------
        %u
        % A1(6*(i-1)+6,3*(i-1)+1)=0;
        % A1(6*(i-1)+6,3*(i-1)+2)=0;
        A1(6*(i-1)+6,3*(i-1)+3)=2*mu*sqrt((n+2)/(2*n+3))*(n+1);    
        %\dot u
        % A2(6*(i-1)+6,3*(i-1)+1)=0;
        % A2(6*(i-1)+6,3*(i-1)+2)=0;
        A2(6*(i-1)+6,3*(i-1)+3)=-2*mu*sqrt((n+2)/(2*n+3));   
    else
        %u
        A1(6*(i-1)+1,3*(i-1)+3)=(3*lambda+2*mu)/sqrt(3)/sqrt(2*n+1)*sqrt(n+1)*(n+2);  
        %\dot u
        A2(6*(i-1)+1,3*(i-1)+3)=+(3*lambda+2*mu)/sqrt(3)/sqrt(2*n+1)*sqrt(n+1);
        %\sigma_{n,n+2,2}-----------------------------
        %u
        A1(6*(i-1)+6,3*(i-1)+3)=2*mu*sqrt((n+2)/(2*n+3))*(n+1);    
        %\dot u
        A2(6*(i-1)+6,3*(i-1)+3)=-2*mu*sqrt((n+2)/(2*n+3));
    end
    % Coupling terms--------------------------------
    for ireo=1:Nreo % Nreo loop over rheology 
        ialpha=find(Coup(ieq,:,27,ireo)>0);%n_s(ialpha) m_s(ialpha) are the degrees and orders that appear in the equation
        K_nm=variable_rheology(ireo,3);
        mu_nm=variable_rheology(ireo,4);
        for j=1:length(ialpha)
            ia=ialpha(j); %index of coefficient of coupling 
            na=n_s(ia);
            ma=m_s(ia);
            if na>0
            %\sigma_{n,n,0}------------------------------
            %n,0--------
            Cp=Coup(ieq,ia,1,ireo); 
            %u
            A1(6*(i-1)+1,3*(ia-1)+1)=A1(6*(i-1)+1,3*(ia-1)+1)+...
                                    K_nm/sqrt(3)/sqrt(2*na+1)*sqrt(na)*(na-1)*Cp;
            % A1(6*(i-1)+1,3*(ia-1)+2)=A1(6*(i-1)+1,3*(ia-1)+2)+...
            %                         0;
            A1(6*(i-1)+1,3*(ia-1)+3)=A1(6*(i-1)+1,3*(ia-1)+3)+...
                                    K_nm/sqrt(3)/sqrt(2*na+1)*sqrt(na+1)*(na+2)*Cp;
            %\dot u
            A2(6*(i-1)+1,3*(ia-1)+1)=A2(6*(i-1)+1,3*(ia-1)+1)+...
                                    -K_nm/sqrt(3)/sqrt(2*na+1)*sqrt(na)*Cp;
            % A2(6*(i-1)+1,3*(ia-1)+2)=A2(6*(i-1)+1,3*(ia-1)+2)+...
            %                         0;
            A2(6*(i-1)+1,3*(ia-1)+3)=A2(6*(i-1)+1,3*(ia-1)+3)+...
                                    +K_nm/sqrt(3)/sqrt(2*na+1)*sqrt(na+1)*Cp;

            %\sigma_{n,n-2,2}----------------------------
            comp=2;% couplings 7-11
            %n--------
            Cp=Coup(ieq,ia,7,ireo); 
            %u
            A1(6*(i-1)+comp,3*(ia-1)+1)=A1(6*(i-1)+comp,3*(ia-1)+1)+...
                                    2*mu_nm*sqrt((2*na+3)*(2*na+2)/12/(2*na-1)/(2*na+1))*(na-1)/r*Cp;
            % A1(6*(i-1)+comp,3*(ia-1)+2)=A1(6*(i-1)+comp,3*(ia-1)+2)+...
            %                         0;
            A1(6*(i-1)+comp,3*(ia-1)+3)=A1(6*(i-1)+comp,3*(ia-1)+3)+...
                                    2*mu_nm*sqrt(na*(2*na-1)*(na+1)/3/(2*na+3)/(2*na+2)/(2*na+1))*(na+2)/r*Cp;
            %\dot u
            A2(6*(i-1)+comp,3*(ia-1)+1)=A2(6*(i-1)+comp,3*(ia-1)+1)+...
                                    -2*mu_nm*sqrt((2*na+3)*(2*na+2)/12/(2*na-1)/(2*na+1))*Cp;
            % A2(6*(i-1)+comp,3*(ia-1)+2)=A2(6*(i-1)+comp,3*(ia-1)+2)+...
            %                         0;
            A2(6*(i-1)+comp,3*(ia-1)+3)=A2(6*(i-1)+comp,3*(ia-1)+3)+...
                                    2*mu_nm*sqrt(na*(2*na-1)*(na+1)/3/(2*na+3)/(2*na+2)/(2*na+1))*Cp;
            %n-2--------
            Cp=Coup(ieq,ia,8,ireo);
            %u
            A1(6*(i-1)+comp,3*(ia-1)+1)=A1(6*(i-1)+comp,3*(ia-1)+1)+...
                                    2*mu_nm*sqrt((na-1)/(2*na-1))*na/r*Cp;
            % A1(6*(i-1)+comp,3*(ia-1)+2)=A1(6*(i-1)+comp,3*(ia-1)+2)+...
            %                         0;
            % A1(6*(i-1)+comp,3*(ia-1)+3)=A1(6*(i-1)+comp,3*(ia-1)+3)+...
            %                         0;
            %\dot u
            A2(6*(i-1)+comp,3*(ia-1)+1)=A2(6*(i-1)+comp,3*(ia-1)+1)+...
                                    2*mu_nm*sqrt((na-1)/(2*na-1))*Cp;
            % A2(6*(i-1)+comp,3*(ia-1)+2)=A2(6*(i-1)+comp,3*(ia-1)+2)+...
            %                         0;
            % A2(6*(i-1)+comp,3*(ia-1)+3)=A2(6*(i-1)+comp,3*(ia-1)+3)+...
            %                         0;
            %n+2--------
            Cp=Coup(ieq,ia,9,ireo); 
            %u
            % A1(6*(i-1)+comp,3*(ia-1)+1)=A1(6*(i-1)+comp,3*(ia-1)+1)+...
            %                         0;
            % A1(6*(i-1)+comp,3*(ia-1)+2)=A1(6*(i-1)+comp,3*(ia-1)+2)+...
            %                         0;
            A1(6*(i-1)+comp,3*(ia-1)+3)=A1(6*(i-1)+comp,3*(ia-1)+3)+...
                                    2*mu_nm*sqrt((na+2)/(2*na+3))*(na+1)/r*Cp;
            %\dot u
            % A2(6*(i-1)+comp,3*(ia-1)+1)=A2(6*(i-1)+comp,3*(ia-1)+1)+...
            %                         0;
            % A2(6*(i-1)+comp,3*(ia-1)+2)=A2(6*(i-1)+comp,3*(ia-1)+2)+...
            %                         0;
            A2(6*(i-1)+comp,3*(ia-1)+3)=A2(6*(i-1)+comp,3*(ia-1)+3)+...
                                    -2*mu_nm*sqrt((na+2)/(2*na+3))*Cp;
            %n-1,0--------
            Cp=Coup(ieq,ia,10,ireo); 
            %u
            % A1(6*(i-1)+comp,3*(ia-1)+1)=A1(6*(i-1)+comp,3*(ia-1)+1)+...
            %                         0;
            A1(6*(i-1)+comp,3*(ia-1)+2)=A1(6*(i-1)+comp,3*(ia-1)+2)+...
                                    2*mu_nm/sqrt(2)*sqrt((na-1)/(2*na+1))*(na+1)/r*Cp;
            % A1(6*(i-1)+comp,3*(ia-1)+3)=A1(6*(i-1)+comp,3*(ia-1)+3)+...
            %                         0;
            %\dot u
            % A2(6*(i-1)+comp,3*(ia-1)+1)=A2(6*(i-1)+comp,3*(ia-1)+1)+...
            %                         0;
            A2(6*(i-1)+comp,3*(ia-1)+2)=A2(6*(i-1)+comp,3*(ia-1)+2)+...
                                    2*mu_nm/sqrt(2)*sqrt((na-1)/(2*na+1))*Cp;
            % A2(6*(i-1)+comp,3*(ia-1)+3)=A2(6*(i-1)+comp,3*(ia-1)+3)+...
            %                         0;
            %n+1,0--------
            Cp=Coup(ieq,ia,11,ireo);
            %u
            % A1(6*(i-1)+comp,3*(ia-1)+1)=A1(6*(i-1)+comp,3*(ia-1)+1)+...
            %                         0;
            A1(6*(i-1)+comp,3*(ia-1)+2)=A1(6*(i-1)+comp,3*(ia-1)+2)+...
                                    2*mu_nm/sqrt(2)*sqrt((na+2)/(2*na+1))*na/r*Cp;
            % A1(6*(i-1)+comp,3*(ia-1)+3)=A1(6*(i-1)+comp,3*(ia-1)+3)+...
            %                         0;
            %\dot u
            % A2(6*(i-1)+comp,3*(ia-1)+1)=A2(6*(i-1)+comp,3*(ia-1)+1)+...
            %                         0;
            A2(6*(i-1)+comp,3*(ia-1)+2)=A2(6*(i-1)+comp,3*(ia-1)+2)+...
                                    -2*mu_nm/sqrt(2)*sqrt((na+2)/(2*na+1))*Cp;
            % A2(6*(i-1)+comp,3*(ia-1)+3)=A2(6*(i-1)+comp,3*(ia-1)+3)+...
            %                         0;
            %\sigma_{n,n-1,2}-----------------------------
            comp=3; %couplings 17-21
            %n--------
            Cp=Coup(ieq,ia,17,ireo); 
            %u
            A1(6*(i-1)+comp,3*(ia-1)+1)=A1(6*(i-1)+comp,3*(ia-1)+1)+...
                                    2*mu_nm*sqrt((2*na+3)*(2*na+2)/12/(2*na-1)/(2*na+1))*(na-1)/r*Cp;
            % A1(6*(i-1)+comp,3*(ia-1)+2)=A1(6*(i-1)+comp,3*(ia-1)+2)+...
            %                         0;
            A1(6*(i-1)+comp,3*(ia-1)+3)=A1(6*(i-1)+comp,3*(ia-1)+3)+...
                                    2*mu_nm*sqrt(na*(2*na-1)*(na+1)/3/(2*na+3)/(2*na+2)/(2*na+1))*(na+2)/r*Cp;
            %\dot u
            A2(6*(i-1)+comp,3*(ia-1)+1)=A2(6*(i-1)+comp,3*(ia-1)+1)+...
                                    -2*mu_nm*sqrt((2*na+3)*(2*na+2)/12/(2*na-1)/(2*na+1))*Cp;
            % A2(6*(i-1)+comp,3*(ia-1)+2)=A2(6*(i-1)+comp,3*(ia-1)+2)+...
            %                         0;
            A2(6*(i-1)+comp,3*(ia-1)+3)=A2(6*(i-1)+comp,3*(ia-1)+3)+...
                                    2*mu_nm*sqrt(na*(2*na-1)*(na+1)/3/(2*na+3)/(2*na+2)/(2*na+1))*Cp;
            %n-2--------
            Cp=Coup(ieq,ia,18,ireo);
            %u
            A1(6*(i-1)+comp,3*(ia-1)+1)=A1(6*(i-1)+comp,3*(ia-1)+1)+...
                                    2*mu_nm*sqrt((na-1)/(2*na-1))*na/r*Cp;
            % A1(6*(i-1)+comp,3*(ia-1)+2)=A1(6*(i-1)+comp,3*(ia-1)+2)+...
            %                         0;
            % A1(6*(i-1)+comp,3*(ia-1)+3)=A1(6*(i-1)+comp,3*(ia-1)+3)+...
            %                         0;
            %\dot u
            A2(6*(i-1)+comp,3*(ia-1)+1)=A2(6*(i-1)+comp,3*(ia-1)+1)+...
                                    2*mu_nm*sqrt((na-1)/(2*na-1))*Cp;
            % A2(6*(i-1)+comp,3*(ia-1)+2)=A2(6*(i-1)+comp,3*(ia-1)+2)+...
            %                         0;
            % A2(6*(i-1)+comp,3*(ia-1)+3)=A2(6*(i-1)+comp,3*(ia-1)+3)+...
            %                         0;
            %n+2--------
            Cp=Coup(ieq,ia,19,ireo); 
            %u
            % A1(6*(i-1)+comp,3*(ia-1)+1)=A1(6*(i-1)+comp,3*(ia-1)+1)+...
            %                         0;
            % A1(6*(i-1)+comp,3*(ia-1)+2)=A1(6*(i-1)+comp,3*(ia-1)+2)+...
            %                         0;
            A1(6*(i-1)+comp,3*(ia-1)+3)=A1(6*(i-1)+comp,3*(ia-1)+3)+...
                                    2*mu_nm*sqrt((na+2)/(2*na+3))*(na+1)/r*Cp;
            %\dot u
            % A2(6*(i-1)+comp,3*(ia-1)+1)=A2(6*(i-1)+comp,3*(ia-1)+1)+...
            %                         0;
            % A2(6*(i-1)+comp,3*(ia-1)+2)=A2(6*(i-1)+comp,3*(ia-1)+2)+...
            %                         0;
            A2(6*(i-1)+comp,3*(ia-1)+3)=A2(6*(i-1)+comp,3*(ia-1)+3)+...
                                    -2*mu_nm*sqrt((na+2)/(2*na+3))*Cp;
            %n-1,0--------
            Cp=Coup(ieq,ia,20,ireo); 
            %u
            % A1(6*(i-1)+comp,3*(ia-1)+1)=A1(6*(i-1)+comp,3*(ia-1)+1)+...
            %                         0;
            A1(6*(i-1)+comp,3*(ia-1)+2)=A1(6*(i-1)+comp,3*(ia-1)+2)+...
                                    2*mu_nm/sqrt(2)*sqrt((na-1)/(2*na+1))*(na+1)/r*Cp;
            % A1(6*(i-1)+comp,3*(ia-1)+3)=A1(6*(i-1)+comp,3*(ia-1)+3)+...
            %                         0;
            %\dot u
            % A2(6*(i-1)+comp,3*(ia-1)+1)=A2(6*(i-1)+comp,3*(ia-1)+1)+...
            %                         0;
            A2(6*(i-1)+comp,3*(ia-1)+2)=A2(6*(i-1)+comp,3*(ia-1)+2)+...
                                    2*mu_nm/sqrt(2)*sqrt((na-1)/(2*na+1))*Cp;
            % A2(6*(i-1)+comp,3*(ia-1)+3)=A2(6*(i-1)+comp,3*(ia-1)+3)+...
            %                         0;
            %n+1,0--------
            Cp=Coup(ieq,ia,21,ireo);
            %u
            % A1(6*(i-1)+comp,3*(ia-1)+1)=A1(6*(i-1)+comp,3*(ia-1)+1)+...
            %                         0;
            A1(6*(i-1)+comp,3*(ia-1)+2)=A1(6*(i-1)+comp,3*(ia-1)+2)+...
                                    2*mu_nm/sqrt(2)*sqrt((na+2)/(2*na+1))*na/r*Cp;
            % A1(6*(i-1)+comp,3*(ia-1)+3)=A1(6*(i-1)+comp,3*(ia-1)+3)+...
            %                         0;
            %\dot u
            % A2(6*(i-1)+comp,3*(ia-1)+1)=A2(6*(i-1)+comp,3*(ia-1)+1)+...
            %                         0;
            A2(6*(i-1)+comp,3*(ia-1)+2)=A2(6*(i-1)+comp,3*(ia-1)+2)+...
                                    -2*mu_nm/sqrt(2)*sqrt((na+2)/(2*na+1))*Cp;
            % A2(6*(i-1)+comp,3*(ia-1)+3)=A2(6*(i-1)+comp,3*(ia-1)+3)+...
            %                         0;

            %\sigma_{n,n,2}-------------------------------
            comp=4; %couplings 2-6
            %n--------
            %n--------
            Cp=Coup(ieq,ia,2,ireo); 
            %u
            A1(6*(i-1)+comp,3*(ia-1)+1)=A1(6*(i-1)+comp,3*(ia-1)+1)+...
                                    2*mu_nm*sqrt((2*na+3)*(2*na+2)/12/(2*na-1)/(2*na+1))*(na-1)/r*Cp;
            % A1(6*(i-1)+comp,3*(ia-1)+2)=A1(6*(i-1)+comp,3*(ia-1)+2)+...
            %                         0;
            A1(6*(i-1)+comp,3*(ia-1)+3)=A1(6*(i-1)+comp,3*(ia-1)+3)+...
                                    2*mu_nm*sqrt(na*(2*na-1)*(na+1)/3/(2*na+3)/(2*na+2)/(2*na+1))*(na+2)/r*Cp;
            %\dot u
            A2(6*(i-1)+comp,3*(ia-1)+1)=A2(6*(i-1)+comp,3*(ia-1)+1)+...
                                    -2*mu_nm*sqrt((2*na+3)*(2*na+2)/12/(2*na-1)/(2*na+1))*Cp;
            % A2(6*(i-1)+comp,3*(ia-1)+2)=A2(6*(i-1)+comp,3*(ia-1)+2)+...
            %                         0;
            A2(6*(i-1)+comp,3*(ia-1)+3)=A2(6*(i-1)+comp,3*(ia-1)+3)+...
                                    2*mu_nm*sqrt(na*(2*na-1)*(na+1)/3/(2*na+3)/(2*na+2)/(2*na+1))*Cp;
            %n-2--------
            Cp=Coup(ieq,ia,3,ireo);
            %u
            A1(6*(i-1)+comp,3*(ia-1)+1)=A1(6*(i-1)+comp,3*(ia-1)+1)+...
                                    2*mu_nm*sqrt((na-1)/(2*na-1))*na/r*Cp;
            % A1(6*(i-1)+comp,3*(ia-1)+2)=A1(6*(i-1)+comp,3*(ia-1)+2)+...
            %                         0;
            % A1(6*(i-1)+comp,3*(ia-1)+3)=A1(6*(i-1)+comp,3*(ia-1)+3)+...
            %                         0;
            %\dot u
            A2(6*(i-1)+comp,3*(ia-1)+1)=A2(6*(i-1)+comp,3*(ia-1)+1)+...
                                    2*mu_nm*sqrt((na-1)/(2*na-1))*Cp;
            % A2(6*(i-1)+comp,3*(ia-1)+2)=A2(6*(i-1)+comp,3*(ia-1)+2)+...
            %                         0;
            % A2(6*(i-1)+comp,3*(ia-1)+3)=A2(6*(i-1)+comp,3*(ia-1)+3)+...
            %                         0;
            %n+2--------
            Cp=Coup(ieq,ia,4,ireo); 
            %u
            % A1(6*(i-1)+comp,3*(ia-1)+1)=A1(6*(i-1)+comp,3*(ia-1)+1)+...
            %                         0;
            % A1(6*(i-1)+comp,3*(ia-1)+2)=A1(6*(i-1)+comp,3*(ia-1)+2)+...
            %                         0;
            A1(6*(i-1)+comp,3*(ia-1)+3)=A1(6*(i-1)+comp,3*(ia-1)+3)+...
                                    2*mu_nm*sqrt((na+2)/(2*na+3))*(na+1)/r*Cp;
            %\dot u
            % A2(6*(i-1)+comp,3*(ia-1)+1)=A2(6*(i-1)+comp,3*(ia-1)+1)+...
            %                         0;
            % A2(6*(i-1)+comp,3*(ia-1)+2)=A2(6*(i-1)+comp,3*(ia-1)+2)+...
            %                         0;
            A2(6*(i-1)+comp,3*(ia-1)+3)=A2(6*(i-1)+comp,3*(ia-1)+3)+...
                                    -2*mu_nm*sqrt((na+2)/(2*na+3))*Cp;
            %n-1,0--------
            Cp=Coup(ieq,ia,5,ireo); 
            %u
            % A1(6*(i-1)+comp,3*(ia-1)+1)=A1(6*(i-1)+comp,3*(ia-1)+1)+...
            %                         0;
            A1(6*(i-1)+comp,3*(ia-1)+2)=A1(6*(i-1)+comp,3*(ia-1)+2)+...
                                    2*mu_nm/sqrt(2)*sqrt((na-1)/(2*na+1))*(na+1)/r*Cp;
            % A1(6*(i-1)+comp,3*(ia-1)+3)=A1(6*(i-1)+comp,3*(ia-1)+3)+...
            %                         0;
            %\dot u
            % A2(6*(i-1)+comp,3*(ia-1)+1)=A2(6*(i-1)+comp,3*(ia-1)+1)+...
            %                         0;
            A2(6*(i-1)+comp,3*(ia-1)+2)=A2(6*(i-1)+comp,3*(ia-1)+2)+...
                                    2*mu_nm/sqrt(2)*sqrt((na-1)/(2*na+1))*Cp;
            % A2(6*(i-1)+comp,3*(ia-1)+3)=A2(6*(i-1)+comp,3*(ia-1)+3)+...
            %                         0;
            %n+1,0--------
            Cp=Coup(ieq,ia,6,ireo);
            %u
            % A1(6*(i-1)+comp,3*(ia-1)+1)=A1(6*(i-1)+comp,3*(ia-1)+1)+...
            %                         0;
            A1(6*(i-1)+comp,3*(ia-1)+2)=A1(6*(i-1)+comp,3*(ia-1)+2)+...
                                    2*mu_nm/sqrt(2)*sqrt((na+2)/(2*na+1))*na/r*Cp;
            % A1(6*(i-1)+comp,3*(ia-1)+3)=A1(6*(i-1)+comp,3*(ia-1)+3)+...
            %                         0;
            %\dot u
            % A2(6*(i-1)+comp,3*(ia-1)+1)=A2(6*(i-1)+comp,3*(ia-1)+1)+...
            %                         0;
            A2(6*(i-1)+comp,3*(ia-1)+2)=A2(6*(i-1)+comp,3*(ia-1)+2)+...
                                    -2*mu_nm/sqrt(2)*sqrt((na+2)/(2*na+1))*Cp;
            % A2(6*(i-1)+comp,3*(ia-1)+3)=A2(6*(i-1)+comp,3*(ia-1)+3)+...
            %                         0;
 

            %\sigma_{n,n+1,2}-----------------------------
            comp=5; %couplings 22-26
            %n--------
            %n--------
            Cp=Coup(ieq,ia,22,ireo); 
            %u
            A1(6*(i-1)+comp,3*(ia-1)+1)=A1(6*(i-1)+comp,3*(ia-1)+1)+...
                                    2*mu_nm*sqrt((2*na+3)*(2*na+2)/12/(2*na-1)/(2*na+1))*(na-1)/r*Cp;
            % A1(6*(i-1)+comp,3*(ia-1)+2)=A1(6*(i-1)+comp,3*(ia-1)+2)+...
            %                         0;
            A1(6*(i-1)+comp,3*(ia-1)+3)=A1(6*(i-1)+comp,3*(ia-1)+3)+...
                                    2*mu_nm*sqrt(na*(2*na-1)*(na+1)/3/(2*na+3)/(2*na+2)/(2*na+1))*(na+2)/r*Cp;
            %\dot u
            A2(6*(i-1)+comp,3*(ia-1)+1)=A2(6*(i-1)+comp,3*(ia-1)+1)+...
                                    -2*mu_nm*sqrt((2*na+3)*(2*na+2)/12/(2*na-1)/(2*na+1))*Cp;
            % A2(6*(i-1)+comp,3*(ia-1)+2)=A2(6*(i-1)+comp,3*(ia-1)+2)+...
            %                         0;
            A2(6*(i-1)+comp,3*(ia-1)+3)=A2(6*(i-1)+comp,3*(ia-1)+3)+...
                                    2*mu_nm*sqrt(na*(2*na-1)*(na+1)/3/(2*na+3)/(2*na+2)/(2*na+1))*Cp;
            %n-2--------
            Cp=Coup(ieq,ia,23,ireo);
            %u
            A1(6*(i-1)+comp,3*(ia-1)+1)=A1(6*(i-1)+comp,3*(ia-1)+1)+...
                                    2*mu_nm*sqrt((na-1)/(2*na-1))*na/r*Cp;
            % A1(6*(i-1)+comp,3*(ia-1)+2)=A1(6*(i-1)+comp,3*(ia-1)+2)+...
            %                         0;
            % A1(6*(i-1)+comp,3*(ia-1)+3)=A1(6*(i-1)+comp,3*(ia-1)+3)+...
            %                         0;
            %\dot u
            A2(6*(i-1)+comp,3*(ia-1)+1)=A2(6*(i-1)+comp,3*(ia-1)+1)+...
                                    2*mu_nm*sqrt((na-1)/(2*na-1))*Cp;
            % A2(6*(i-1)+comp,3*(ia-1)+2)=A2(6*(i-1)+comp,3*(ia-1)+2)+...
            %                         0;
            % A2(6*(i-1)+comp,3*(ia-1)+3)=A2(6*(i-1)+comp,3*(ia-1)+3)+...
            %                         0;
            %n+2--------
            Cp=Coup(ieq,ia,24,ireo); 
            %u
            % A1(6*(i-1)+comp,3*(ia-1)+1)=A1(6*(i-1)+comp,3*(ia-1)+1)+...
            %                         0;
            % A1(6*(i-1)+comp,3*(ia-1)+2)=A1(6*(i-1)+comp,3*(ia-1)+2)+...
            %                         0;
            A1(6*(i-1)+comp,3*(ia-1)+3)=A1(6*(i-1)+comp,3*(ia-1)+3)+...
                                    2*mu_nm*sqrt((na+2)/(2*na+3))*(na+1)/r*Cp;
            %\dot u
            % A2(6*(i-1)+comp,3*(ia-1)+1)=A2(6*(i-1)+comp,3*(ia-1)+1)+...
            %                         0;
            % A2(6*(i-1)+comp,3*(ia-1)+2)=A2(6*(i-1)+comp,3*(ia-1)+2)+...
            %                         0;
            A2(6*(i-1)+comp,3*(ia-1)+3)=A2(6*(i-1)+comp,3*(ia-1)+3)+...
                                    -2*mu_nm*sqrt((na+2)/(2*na+3))*Cp;
            %n-1,0--------
            Cp=Coup(ieq,ia,25,ireo); 
            %u
            % A1(6*(i-1)+comp,3*(ia-1)+1)=A1(6*(i-1)+comp,3*(ia-1)+1)+...
            %                         0;
            A1(6*(i-1)+comp,3*(ia-1)+2)=A1(6*(i-1)+comp,3*(ia-1)+2)+...
                                    2*mu_nm/sqrt(2)*sqrt((na-1)/(2*na+1))*(na+1)/r*Cp;
            % A1(6*(i-1)+comp,3*(ia-1)+3)=A1(6*(i-1)+comp,3*(ia-1)+3)+...
            %                         0;
            %\dot u
            % A2(6*(i-1)+comp,3*(ia-1)+1)=A2(6*(i-1)+comp,3*(ia-1)+1)+...
            %                         0;
            A2(6*(i-1)+comp,3*(ia-1)+2)=A2(6*(i-1)+comp,3*(ia-1)+2)+...
                                    2*mu_nm/sqrt(2)*sqrt((na-1)/(2*na+1))*Cp;
            % A2(6*(i-1)+comp,3*(ia-1)+3)=A2(6*(i-1)+comp,3*(ia-1)+3)+...
            %                         0;
            %n+1,0--------
            Cp=Coup(ieq,ia,26,ireo);
            %u
            % A1(6*(i-1)+comp,3*(ia-1)+1)=A1(6*(i-1)+comp,3*(ia-1)+1)+...
            %                         0;
            A1(6*(i-1)+comp,3*(ia-1)+2)=A1(6*(i-1)+comp,3*(ia-1)+2)+...
                                    2*mu_nm/sqrt(2)*sqrt((na+2)/(2*na+1))*na/r*Cp;
            % A1(6*(i-1)+comp,3*(ia-1)+3)=A1(6*(i-1)+comp,3*(ia-1)+3)+...
            %                         0;
            %\dot u
            % A2(6*(i-1)+comp,3*(ia-1)+1)=A2(6*(i-1)+comp,3*(ia-1)+1)+...
            %                         0;
            A2(6*(i-1)+comp,3*(ia-1)+2)=A2(6*(i-1)+comp,3*(ia-1)+2)+...
                                    -2*mu_nm/sqrt(2)*sqrt((na+2)/(2*na+1))*Cp;
            % A2(6*(i-1)+comp,3*(ia-1)+3)=A2(6*(i-1)+comp,3*(ia-1)+3)+...
            %                         0;


            %\sigma_{n,n+2,2}-----------------------------
            comp=6; %couplings 12-16
            %n--------
            %n--------
            Cp=Coup(ieq,ia,12,ireo); 
            %u
            A1(6*(i-1)+comp,3*(ia-1)+1)=A1(6*(i-1)+comp,3*(ia-1)+1)+...
                                    2*mu_nm*sqrt((2*na+3)*(2*na+2)/12/(2*na-1)/(2*na+1))*(na-1)/r*Cp;
            % A1(6*(i-1)+comp,3*(ia-1)+2)=A1(6*(i-1)+comp,3*(ia-1)+2)+...
            %                         0;
            A1(6*(i-1)+comp,3*(ia-1)+3)=A1(6*(i-1)+comp,3*(ia-1)+3)+...
                                    2*mu_nm*sqrt(na*(2*na-1)*(na+1)/3/(2*na+3)/(2*na+2)/(2*na+1))*(na+2)/r*Cp;
            %\dot u
            A2(6*(i-1)+comp,3*(ia-1)+1)=A2(6*(i-1)+comp,3*(ia-1)+1)+...
                                    -2*mu_nm*sqrt((2*na+3)*(2*na+2)/12/(2*na-1)/(2*na+1))*Cp;
            % A2(6*(i-1)+comp,3*(ia-1)+2)=A2(6*(i-1)+comp,3*(ia-1)+2)+...
            %                         0;
            A2(6*(i-1)+comp,3*(ia-1)+3)=A2(6*(i-1)+comp,3*(ia-1)+3)+...
                                    2*mu_nm*sqrt(na*(2*na-1)*(na+1)/3/(2*na+3)/(2*na+2)/(2*na+1))*Cp;
            %n-2--------
            Cp=Coup(ieq,ia,13,ireo);
            %u
            A1(6*(i-1)+comp,3*(ia-1)+1)=A1(6*(i-1)+comp,3*(ia-1)+1)+...
                                    2*mu_nm*sqrt((na-1)/(2*na-1))*na/r*Cp;
            % A1(6*(i-1)+comp,3*(ia-1)+2)=A1(6*(i-1)+comp,3*(ia-1)+2)+...
            %                         0;
            % A1(6*(i-1)+comp,3*(ia-1)+3)=A1(6*(i-1)+comp,3*(ia-1)+3)+...
            %                         0;
            %\dot u
            A2(6*(i-1)+comp,3*(ia-1)+1)=A2(6*(i-1)+comp,3*(ia-1)+1)+...
                                    2*mu_nm*sqrt((na-1)/(2*na-1))*Cp;
            % A2(6*(i-1)+comp,3*(ia-1)+2)=A2(6*(i-1)+comp,3*(ia-1)+2)+...
            %                         0;
            % A2(6*(i-1)+comp,3*(ia-1)+3)=A2(6*(i-1)+comp,3*(ia-1)+3)+...
            %                         0;
            %n+2--------
            Cp=Coup(ieq,ia,14,ireo); 
            %u
            % A1(6*(i-1)+comp,3*(ia-1)+1)=A1(6*(i-1)+comp,3*(ia-1)+1)+...
            %                         0;
            % A1(6*(i-1)+comp,3*(ia-1)+2)=A1(6*(i-1)+comp,3*(ia-1)+2)+...
            %                         0;
            A1(6*(i-1)+comp,3*(ia-1)+3)=A1(6*(i-1)+comp,3*(ia-1)+3)+...
                                    2*mu_nm*sqrt((na+2)/(2*na+3))*(na+1)/r*Cp;
            %\dot u
            % A2(6*(i-1)+comp,3*(ia-1)+1)=A2(6*(i-1)+comp,3*(ia-1)+1)+...
            %                         0;
            % A2(6*(i-1)+comp,3*(ia-1)+2)=A2(6*(i-1)+comp,3*(ia-1)+2)+...
            %                         0;
            A2(6*(i-1)+comp,3*(ia-1)+3)=A2(6*(i-1)+comp,3*(ia-1)+3)+...
                                    -2*mu_nm*sqrt((na+2)/(2*na+3))*Cp;
            %n-1,0--------
            Cp=Coup(ieq,ia,15,ireo); 
            %u
            % A1(6*(i-1)+comp,3*(ia-1)+1)=A1(6*(i-1)+comp,3*(ia-1)+1)+...
            %                         0;
            A1(6*(i-1)+comp,3*(ia-1)+2)=A1(6*(i-1)+comp,3*(ia-1)+2)+...
                                    2*mu_nm/sqrt(2)*sqrt((na-1)/(2*na+1))*(na+1)/r*Cp;
            % A1(6*(i-1)+comp,3*(ia-1)+3)=A1(6*(i-1)+comp,3*(ia-1)+3)+...
            %                         0;
            %\dot u
            % A2(6*(i-1)+comp,3*(ia-1)+1)=A2(6*(i-1)+comp,3*(ia-1)+1)+...
            %                         0;
            A2(6*(i-1)+comp,3*(ia-1)+2)=A2(6*(i-1)+comp,3*(ia-1)+2)+...
                                    2*mu_nm/sqrt(2)*sqrt((na-1)/(2*na+1))*Cp;
            % A2(6*(i-1)+comp,3*(ia-1)+3)=A2(6*(i-1)+comp,3*(ia-1)+3)+...
            %                         0;
            %n+1,0--------
            Cp=Coup(ieq,ia,16,ireo);
            %u
            % A1(6*(i-1)+comp,3*(ia-1)+1)=A1(6*(i-1)+comp,3*(ia-1)+1)+...
            %                         0;
            A1(6*(i-1)+comp,3*(ia-1)+2)=A1(6*(i-1)+comp,3*(ia-1)+2)+...
                                    2*mu_nm/sqrt(2)*sqrt((na+2)/(2*na+1))*na/r*Cp;
            % A1(6*(i-1)+comp,3*(ia-1)+3)=A1(6*(i-1)+comp,3*(ia-1)+3)+...
            %                         0;
            %\dot u
            % A2(6*(i-1)+comp,3*(ia-1)+1)=A2(6*(i-1)+comp,3*(ia-1)+1)+...
            %                         0;
            A2(6*(i-1)+comp,3*(ia-1)+2)=A2(6*(i-1)+comp,3*(ia-1)+2)+...
                                    -2*mu_nm/sqrt(2)*sqrt((na+2)/(2*na+1))*Cp;
            % A2(6*(i-1)+comp,3*(ia-1)+3)=A2(6*(i-1)+comp,3*(ia-1)+3)+...
            %                         0;
            else
            %\sigma_{n,n,-2}------------------------------
            comp=2; 
            %n+2--------
            Cp=Coup(ieq,ia,9,ireo); 
            %u
            A1(6*(i-1)+comp,3*(ia-1)+3)=A1(6*(i-1)+comp,3*(ia-1)+3)+...
                                    2*mu_nm*sqrt((na+2)/(2*na+3))*(na+1)/r*Cp;
            %\dot u
            A2(6*(i-1)+comp,3*(ia-1)+3)=A2(6*(i-1)+comp,3*(ia-1)+3)+...
                                    -2*mu_nm*sqrt((na+2)/(2*na+3))*Cp;
            %\sigma_{n,n-1,2}-----------------------------
            comp=3; %couplings 17-21
            Cp=Coup(ieq,ia,19,ireo); 
            %u
            % A1(6*(i-1)+comp,3*(ia-1)+1)=A1(6*(i-1)+comp,3*(ia-1)+1)+...
            %                         0;
            % A1(6*(i-1)+comp,3*(ia-1)+2)=A1(6*(i-1)+comp,3*(ia-1)+2)+...
            %                         0;
            A1(6*(i-1)+comp,3*(ia-1)+3)=A1(6*(i-1)+comp,3*(ia-1)+3)+...
                                    2*mu_nm*sqrt((na+2)/(2*na+3))*(na+1)/r*Cp;
            %\dot u
            % A2(6*(i-1)+comp,3*(ia-1)+1)=A2(6*(i-1)+comp,3*(ia-1)+1)+...
            %                         0;
            % A2(6*(i-1)+comp,3*(ia-1)+2)=A2(6*(i-1)+comp,3*(ia-1)+2)+...
            %                         0;
            A2(6*(i-1)+comp,3*(ia-1)+3)=A2(6*(i-1)+comp,3*(ia-1)+3)+...
                                    -2*mu_nm*sqrt((na+2)/(2*na+3))*Cp;           
            %\sigma_{n,n,2}-------------------------------
            comp=4; %couplings 2-6
            %n+2--------
            Cp=Coup(ieq,ia,4,ireo); 
            %u
            % A1(6*(i-1)+comp,3*(ia-1)+1)=A1(6*(i-1)+comp,3*(ia-1)+1)+...
            %                         0;
            % A1(6*(i-1)+comp,3*(ia-1)+2)=A1(6*(i-1)+comp,3*(ia-1)+2)+...
            %                         0;
            A1(6*(i-1)+comp,3*(ia-1)+3)=A1(6*(i-1)+comp,3*(ia-1)+3)+...
                                    2*mu_nm*sqrt((na+2)/(2*na+3))*(na+1)/r*Cp;
            %\dot u
            % A2(6*(i-1)+comp,3*(ia-1)+1)=A2(6*(i-1)+comp,3*(ia-1)+1)+...
            %                         0;
            % A2(6*(i-1)+comp,3*(ia-1)+2)=A2(6*(i-1)+comp,3*(ia-1)+2)+...
            %                         0;
            A2(6*(i-1)+comp,3*(ia-1)+3)=A2(6*(i-1)+comp,3*(ia-1)+3)+...
                                    -2*mu_nm*sqrt((na+2)/(2*na+3))*Cp;            
            %\sigma_{n,n+1,2}-----------------------------
            comp=5; %couplings 22-26            
            %n+2--------
            Cp=Coup(ieq,ia,24,ireo); 
            %u
            % A1(6*(i-1)+comp,3*(ia-1)+1)=A1(6*(i-1)+comp,3*(ia-1)+1)+...
            %                         0;
            % A1(6*(i-1)+comp,3*(ia-1)+2)=A1(6*(i-1)+comp,3*(ia-1)+2)+...
            %                         0;
            A1(6*(i-1)+comp,3*(ia-1)+3)=A1(6*(i-1)+comp,3*(ia-1)+3)+...
                                    2*mu_nm*sqrt((na+2)/(2*na+3))*(na+1)/r*Cp;
            %\dot u
            % A2(6*(i-1)+comp,3*(ia-1)+1)=A2(6*(i-1)+comp,3*(ia-1)+1)+...
            %                         0;
            % A2(6*(i-1)+comp,3*(ia-1)+2)=A2(6*(i-1)+comp,3*(ia-1)+2)+...
            %                         0;
            A2(6*(i-1)+comp,3*(ia-1)+3)=A2(6*(i-1)+comp,3*(ia-1)+3)+...
                                    -2*mu_nm*sqrt((na+2)/(2*na+3))*Cp;            
            %\sigma_{n,n+2,2}-----------------------------
            comp=6; %couplings 12-16
            %n+2--------
            Cp=Coup(ieq,ia,14,ireo); 
            %u
            % A1(6*(i-1)+comp,3*(ia-1)+1)=A1(6*(i-1)+comp,3*(ia-1)+1)+...
            %                         0;
            % A1(6*(i-1)+comp,3*(ia-1)+2)=A1(6*(i-1)+comp,3*(ia-1)+2)+...
            %                         0;
            A1(6*(i-1)+comp,3*(ia-1)+3)=A1(6*(i-1)+comp,3*(ia-1)+3)+...
                                    2*mu_nm*sqrt((na+2)/(2*na+3))*(na+1)/r*Cp;
            %\dot u
            % A2(6*(i-1)+comp,3*(ia-1)+1)=A2(6*(i-1)+comp,3*(ia-1)+1)+...
            %                         0;
            % A2(6*(i-1)+comp,3*(ia-1)+2)=A2(6*(i-1)+comp,3*(ia-1)+2)+...
            %                         0;
            A2(6*(i-1)+comp,3*(ia-1)+3)=A2(6*(i-1)+comp,3*(ia-1)+3)+...
                                    -2*mu_nm*sqrt((na+2)/(2*na+3))*Cp;
            end
        end
    end
end
end
%% 
% \epsilon=A14\dot u+A15/r u
function [A14,A15]=get_A14A15(Couplings)
N=length(Couplings.n_s);
A14=zeros(6*N,3*N);
A15=zeros(6*N,3*N);
for i=1:N
    n=Couplings.n_s(i);
    if n>0
    %\sigma_{n,n,0}------------------------------
    %u dot
    A14(6*(i-1)+1,3*(i-1)+1)=-1/sqrt(3)*1/sqrt(2*n+1)*sqrt(n);
    % A14(6*(i-1)+1,3*(i-1)+2)=0;
    A14(6*(i-1)+1,3*(i-1)+3)=1/sqrt(3)*1/sqrt(2*n+1)*sqrt(n+1);  
    %u
    A15(6*(i-1)+1,3*(i-1)+1)=1/sqrt(3)*1/sqrt(2*n+1)*(n-1)*sqrt(n);
    % A15(6*(i-1)+1,3*(i-1)+2)=0;
    A15(6*(i-1)+1,3*(i-1)+3)=1/sqrt(3)*1/sqrt(2*n+1)*sqrt(n+1)*(n+2);
    %\sigma_{n,n-2,2}------------------------------
    %u dot
    A14(6*(i-1)+2,3*(i-1)+1)=sqrt((n-1)/(2*n-1));
    % A14(6*(i-1)+2,3*(i-1)+2)=0;
    % A14(6*(i-1)+2,3*(i-1)+3)=0;  
    %u
    A15(6*(i-1)+2,3*(i-1)+1)=sqrt((n-1)/(2*n-1))*n;
    % A15(6*(i-1)+2,3*(i-1)+2)=0;
    % A15(6*(i-1)+2,3*(i-1)+3)=0;
    %\sigma_{n,n-1,2}------------------------------
    %u dot
    % A14(6*(i-1)+3,3*(i-1)+1)=0;
    A14(6*(i-1)+3,3*(i-1)+2)=1/sqrt(2)*sqrt((n-1)/(2*n+1));
    % A14(6*(i-1)+3,3*(i-1)+3)=0;  
    %u
    % A15(6*(i-1)+3,3*(i-1)+1)=0;
    A15(6*(i-1)+3,3*(i-1)+2)=1/sqrt(2)*sqrt((n-1)/(2*n+1))*(n+1);
    % A15(6*(i-1)+3,3*(i-1)+3)=0;
    %\sigma_{n,n,2}------------------------------
    %u dot
    A14(6*(i-1)+4,3*(i-1)+1)=-sqrt((2*n+3)*(2*n+2)/12/(2*n-1)/(2*n+1));
    % A14(6*(i-1)+4,3*(i-1)+2)=0;
    A14(6*(i-1)+4,3*(i-1)+3)=sqrt(n*(2*n-1)*(n+1)/3/(2*n+3)/(2*n+2)/(2*n+1));  
    %u
    A15(6*(i-1)+4,3*(i-1)+1)=sqrt((2*n+3)*(2*n+2)/12/(2*n-1)/(2*n+1))*(n-1);
    % A15(6*(i-1)+4,3*(i-1)+2)=0;
    A15(6*(i-1)+4,3*(i-1)+3)=sqrt(n*(2*n-1)*(n+1)/3/(2*n+3)/(2*n+2)/(2*n+1))*(n+2);
    %\sigma_{n,n+1,2}------------------------------
    %u dot
    % A14(6*(i-1)+5,3*(i-1)+1)=0;
    A14(6*(i-1)+5,3*(i-1)+2)=-1/sqrt(2)*sqrt((n+2)/(2*n+1));
    % A14(6*(i-1)+5,3*(i-1)+3)=0;  
    %u
    % A15(6*(i-1)+5,3*(i-1)+1)=0;
    A15(6*(i-1)+5,3*(i-1)+2)=1/sqrt(2)*sqrt((n+2)/(2*n+1))*n;
    % A15(6*(i-1)+5,3*(i-1)+3)=0;
    %\sigma_{n,n+2,2}------------------------------
    %u dot
    % A14(6*(i-1)+6,3*(i-1)+1)=0;
    % A14(6*(i-1)+6,3*(i-1)+2)=0;
    A14(6*(i-1)+6,3*(i-1)+3)=-sqrt((n+2)/(2*n+3));  
    %u
    % A15(6*(i-1)+6,3*(i-1)+1)=0;
    % A15(6*(i-1)+6,3*(i-1)+2)=0;
    A15(6*(i-1)+6,3*(i-1)+3)=sqrt((n+2)/(2*n+3))*(n+1);
    else
    %\sigma_{n,n,0}------------------------------
    %u dot
    A14(6*(i-1)+1,3*(i-1)+3)=1/sqrt(3)*1/sqrt(2*n+1)*sqrt(n+1);  
    %u
    A15(6*(i-1)+1,3*(i-1)+3)=1/sqrt(3)*1/sqrt(2*n+1)*sqrt(n+1)*(n+2);
    %\sigma_{n,n+2,2}------------------------------
    %u dot
    A14(6*(i-1)+6,3*(i-1)+3)=-sqrt((n+2)/(2*n+3));  
    %u
    A15(6*(i-1)+6,3*(i-1)+3)=sqrt((n+2)/(2*n+3))*(n+1);    
    end
end
end

%% get_A3
% U=A3u with
% U=[U_n^m, V_n^m, W_n^m]
function [A3]=get_A3(Couplings)
N=length(Couplings.n_s);
A3=zeros(3*N,3*N);
for i=1:N
    n=Couplings.n_s(i);
    if n>0
    %U_n^m
    A3(3*(i-1)+1,3*(i-1)+1)=sqrt(n)/sqrt(2*n+1); %u_{n-1}
    % A3(3*(i-1)+1,3*(i-1)+2)=0; %u_{n}
    A3(3*(i-1)+1,3*(i-1)+3)=-sqrt(n+1)/sqrt(2*n+1); %u_{n+1}
    %V_n^m
    A3(3*(i-1)+2,3*(i-1)+1)=1/sqrt(2*n+1)/sqrt(n); %u_{n-1}
    % A3(3*(i-1)+2,3*(i-1)+2)=0; %u_{n}
    A3(3*(i-1)+2,3*(i-1)+3)=1/sqrt(2*n+1)/sqrt(n+1); %u_{n+1}
    %W_n^m
    % A3(3*(i-1)+3,3*(i-1)+1)=0; %u_{n-1}
    A3(3*(i-1)+3,3*(i-1)+2)=1i/sqrt(n*(n+1)); %u_{n}
    % A3(3*(i-1)+3,3*(i-1)+3)=0; %u_{n+1}  
    else
    %U_n^m    
    A3(3*(i-1)+1,3*(i-1)+3)=-sqrt(n+1)/sqrt(2*n+1); %u_{n+1} 
    A3(3*(i-1)+3,3*(i-1)+2)=1; %u_{n}
    A3(3*(i-1)+2,3*(i-1)+1)=1;
    end
end
end
%% get_A4
% \Sigma=A4\sigma
% \Sigma=[R_n^m, S_n^m, W_n^m]
%\sigma=[\sigma_{n,n,0} \sigma_{n,n-2,2}, \sigma_{n,n-1,2}, \sigma_{n,n,2}, \sigma_{n,n+1,2}, \sigma_{n,n+2,0}]
function [A4]=get_A4(Couplings)
N=length(Couplings.n_s);
A4=zeros(3*N,6*N);
for i=1:N
    n=Couplings.n_s(i);
    if n>0
    %R_n^m
    A4(3*(i-1)+1,6*(i-1)+1)=-1/sqrt(3); %\sigma_{n,n,0}
    A4(3*(i-1)+1,6*(i-1)+2)=sqrt(n*(n-1)/(2*n-1)/(2*n+1)); %\sigma_{n,n-2,2}
    % A4(3*(i-1)+1,6*(i-1)+3)=0; %\sigma_{n,n-1,2}
    A4(3*(i-1)+1,6*(i-1)+4)=-sqrt(n)/(2*n+1)*(sqrt((2*n+3)*(2*n+2)/12/(2*n-1))+sqrt((2*n-1)*(n+1)^2/3/(2*n+2)/(2*n+3))); %\sigma_{n,n,2}
    % A4(3*(i-1)+1,6*(i-1)+5)=0; %\sigma_{n,n+1,2}
    A4(3*(i-1)+1,6*(i-1)+6)=sqrt((n+1)*(n+2)/(2*n+1)/(2*n+3)); %\sigma_{n,n+2,2}   
    %S_n^m
    % A4(3*(i-1)+2,6*(i-1)+1)=0; %\sigma_{n,n,0}
    A4(3*(i-1)+2,6*(i-1)+2)=sqrt((n-1)/(2*n-1)/n/(2*n+1)); %\sigma_{n,n-2,2}
    % A4(3*(i-1)+2,6*(i-1)+3)=0; %\sigma_{n,n-1,2}
    A4(3*(i-1)+2,6*(i-1)+4)=(-sqrt((2*n+3)*(2*n+2)/(12*n)/(2*n-1)/(2*n+1)^2)+sqrt(n*(2*n-1)/3/(2*n+1)^2/(2*n+2)/(2*n+3))); %\sigma_{n,n,2}
    % A4(3*(i-1)+2,6*(i-1)+5)=0; %\sigma_{n,n+1,2}
    A4(3*(i-1)+2,6*(i-1)+6)=-sqrt((n+2)/(2*n+3)/(n+1)/(2*n+1)); %\sigma_{n,n+2,2} 
    %T_n^m
    % A4(3*(i-1)+3,6*(i-1)+1)=0; %\sigma_{n,n,0}
    % A4(3*(i-1)+3,6*(i-1)+2)=0; %\sigma_{n,n-2,2}
    A4(3*(i-1)+3,6*(i-1)+3)=1i/sqrt(2*n*(n+1)*(2*n+1))*sqrt(n-1); %\sigma_{n,n-1,2}
    % A4(3*(i-1)+3,6*(i-1)+4)=0; %\sigma_{n,n,2}
    A4(3*(i-1)+3,6*(i-1)+5)=-1i/sqrt(2*n*(n+1)*(2*n+1))*sqrt(n+2); %\sigma_{n,n+1,2}
    % A4(3*(i-1)+3,6*(i-1)+6)=0; %\sigma_{n,n+2,2} 
    else
    %R_n^m
    A4(3*(i-1)+1,6*(i-1)+1)=-1/sqrt(3); %\sigma_{n,n,0}
    A4(3*(i-1)+1,6*(i-1)+6)=sqrt((n+1)*(n+2)/(2*n+1)/(2*n+3)); %\sigma_{n,n+2,2}   
    %S_n^m, undefined     
    end
end
end
%% get_A5
% \dot\Sigma=A5*\sigma/r+A6
% with:
%\sigma=[\sigma_{n,n,0}, \sigma_{n,n-2,2}, \sigma_{n,n-1,2}, \sigma_{n,n,2}, \sigma_{n,n+1,2}, \sigma_{n,n+2,0},]
%\Sigma=[R_{n,m} S_{n,m} W_{n,m}]
function [A5]=get_A5(Couplings)
N=length(Couplings.n_s);
A5=zeros(3*N,6*N);
for i=1:N
    n=Couplings.n_s(i);
    if n>0
    %R_n^m
    % A5(3*(i-1)+1,6*(i-1)+1)=0; %\sigma_{n,n,0}
    A5(3*(i-1)+1,6*(i-1)+2)=-(n-2)*sqrt(n*(n-1)/(2*n-1)/(2*n+1)); %\sigma_{n,n-2,2}
    % A5(3*(i-1)+1,6*(i-1)+3)=0; %\sigma_{n,n-1,2}
    A5(3*(i-1)+1,6*(i-1)+4)=1/(2*n+1)*(...
        -(n+1)*sqrt((2*n+3)*(2*n+2)*n/12/(2*n-1))+...
        +n*sqrt(n*(n+1)^2*(2*n-1)/3/(2*n+2)/(2*n+3))); %\sigma_{n,n,2}
    % A5(3*(i-1)+1,6*(i-1)+5)=0; %\sigma_{n,n+1,2}
    A5(3*(i-1)+1,6*(i-1)+6)=(n+3)*sqrt((n+1)*(n+2)/(2*n+3)/(2*n+1)); %\sigma_{n,n+2,2}   
    %S_n^m
    A5(3*(i-1)+2,6*(i-1)+1)=-1/sqrt(3); %\sigma_{n,n,0}
    A5(3*(i-1)+2,6*(i-1)+2)=-(n-2)*sqrt((n-1)/(2*n-1)/(2*n+1)/n); %\sigma_{n,n-2,2}
    % A5(3*(i-1)+2,6*(i-1)+3)=0; %\sigma_{n,n-1,2}
    A5(3*(i-1)+2,6*(i-1)+4)=-1/(2*n+1)*(...
        (n+1)*sqrt((2*n+3)*(2*n+2)/(2*n-1)/12/n)+...
        n*sqrt(n*(2*n-1)/3/(2*n+2)/(2*n+3))); %\sigma_{n,n,2}
    % A5(3*(i-1)+2,6*(i-1)+5)=0; %\sigma_{n,n+1,2}
    A5(3*(i-1)+2,6*(i-1)+6)=-(n+3)*sqrt((n+2)/(2*n+3)/(2*n+1)/(n+1)); %\sigma_{n,n+2,2} 
    %T_n^m
    % A5(3*(i-1)+3,6*(i-1)+1)=0; %\sigma_{n,n,0}
    % A5(3*(i-1)+3,6*(i-1)+2)=0; %\sigma_{n,n-2,2}
    A5(3*(i-1)+3,6*(i-1)+3)=-(n-1)*sqrt(n-1)/sqrt(2*n*(n+1)*(2*n+1))*1i; %\sigma_{n,n-1,2}
    % A5(3*(i-1)+3,6*(i-1)+4)=0; %\sigma_{n,n,2}
    A5(3*(i-1)+3,6*(i-1)+5)=-(n+2)*sqrt(n+2)/sqrt(2*n*(n+1)*(2*n+1))*1i; %\sigma_{n,n+1,2}
    % A5(3*(i-1)+3,6*(i-1)+6)=0; %\sigma_{n,n+2,2}
    else
    %R_n^m
    A5(3*(i-1)+1,6*(i-1)+6)=(n+3)*sqrt((n+1)*(n+2)/(2*n+3)/(2*n+1)); %\sigma_{n,n+2,2}   
    %S_n^m, undefined    
    end
end
A5=-A5;
end
%% 
%A13\dot\Sigma=A5*\sigma/r+A6\dot{U}+g/r*A71*U+dg*A72*U+A81*\Phi+A82*\Phi/r
%A9\dot\Phi=A100\Phi+A101/r*\Phi+A102/r^2*\Phi+A11/r*U+A12\dot{U}
function [A, A6, A71, A72, A81, A82, A9, A100, A101, A102, A11, A12]=get_others(Couplings,Interior_Model)
N=length(Couplings.n_s);
A6=zeros(3*N,3*N);
A=zeros(3*N,3*N);
A71=zeros(3*N,3*N);
A72=zeros(3*N,3*N);
A81=zeros(3*N,2*N);
A82=zeros(3*N,2*N);
A9=zeros(2*N,2*N);
A100=zeros(2*N,2*N);
A101=zeros(2*N,2*N);
A102=zeros(2*N,2*N);
A11=zeros(2*N,3*N);
A12=zeros(2*N,3*N);
G=Interior_Model.Gg;
rho=Interior_Model.rho;
for i=1:N
    n=Couplings.n_s(i);
    % A-----------------------------
    A(3*(i-1)+1,3*(i-1)+1)=1; 
    A(3*(i-1)+2,3*(i-1)+2)=1;
    A(3*(i-1)+3,3*(i-1)+3)=1;
    % A7-----------------------------
    A71(3*(i-1)+1,3*(i-1)+1)=-2*rho;
    A71(3*(i-1)+1,3*(i-1)+2)=+rho*n*(n+1);
    A72(3*(i-1)+1,3*(i-1)+1)=+rho; 
    A71(3*(i-1)+2,3*(i-1)+1)=+rho;
    if n>0
    % A8-----------------------------
        A81(3*(i-1)+1,2*(i-1)+2)=+rho;
        A82(3*(i-1)+2,2*(i-1)+1)=+rho;
        % A9-----------------------------
        A9(2*(i-1)+1,2*(i-1)+1)=1;
        A9(2*(i-1)+2,2*(i-1)+2)=1;    
        % A10-----------------------------
        A100(2*(i-1)+1,2*(i-1)+2)=1;
        A101(2*(i-1)+2,2*(i-1)+2)=-2;
        A102(2*(i-1)+2,2*(i-1)+1)=+n*(n+1);
        % A11-----------------------------
        A11(2*(i-1)+2,3*(i-1)+1)=-2*4*pi*G*rho;
        A11(2*(i-1)+2,3*(i-1)+2)=+4*pi*G*rho*n*(n+1);
        % A12-----------------------------
        A12(2*(i-1)+2,3*(i-1)+1)=-4*pi*G*rho;
    else %from Longman 1963 
        % A9-----------------------------
        A9(2*(i-1)+1,2*(i-1)+1)=1;
        A9(2*(i-1)+2,2*(i-1)+2)=1; 
        %A12(2*(i-1)+1,3*(i-1)+1)=-4*pi*G*rho; 
        %A11(2*(i-1)+1,3*(i-1)+1)=-rho; 
        A100(2*(i-1)+1,2*(i-1)+1)=1;
        A100(2*(i-1)+2,2*(i-1)+2)=1;
    end
end

end
%% 
% function Aprop=get_Aprop(rK,gK,dgK,Nmodes,A1,A2,A3_inv,A4,A5,A6,A71,A72,A81,A82,A9,A100,A101,A102,A11,A12,A13,deg_0,Gg)
% %% (1) Assemble matrix
%     % rheology equation 
%     Adotx1=sparse(3*Nmodes,8*Nmodes);
%     Ax1=sparse(3*Nmodes,8*Nmodes);
%     Adotx1(:,1:3*Nmodes)=A4*A1*A3_inv; %\dot{U}
%     Ax1(:,1:3*Nmodes)=-A4*A2*A3_inv/rK;% U
%     Ax1(:,3*Nmodes+(1:3*Nmodes))=A13;% \Sigma
%     % for uniform model this looks OK. 
%     % momentum equations 
%     Adotx2=sparse(3*Nmodes,8*Nmodes);
%     Ax2=sparse(3*Nmodes,8*Nmodes);
%     Adotx2(:,1:3*Nmodes)=-A5*A1*A3_inv/rK+A6;% \dot{U}
%     Adotx2(:,3*Nmodes+(1:3*Nmodes))=A13;% \dot{\Sigma{U}}
%     Ax2(:,1:3*Nmodes)=A5*A2*A3_inv/rK^2+gK/rK*A71+dgK*A72;% U
%     Ax2(:,6*Nmodes+(1:2*Nmodes))=A81+A82/rK;% \Phi
%     % for uniform model, this looks OK 
%     % Poisson equation 
%     Adotx3=sparse(2*Nmodes,8*Nmodes);
%     Ax3=sparse(2*Nmodes,8*Nmodes);
%     Adotx3(:,1:3*Nmodes)=-A12; %\dot{U}
%     Adotx3(:,6*Nmodes+(1:2*Nmodes))=A9; %\dot{\Phi}
%     Ax3(:,1:3*Nmodes)=A11/rK; % U
%     Ax3(:,6*Nmodes+(1:2*Nmodes))=A100+A101/rK+A102/rK^2; % \Phi
%     % combine matrices
%     Adotx=[Adotx1; Adotx2; Adotx3];
%     Ax=[Ax1; Ax2; Ax3];
%     if deg_0==1
%         Ax([2,2+3*Nmodes],:)=0; 
%         Ax([3,3+3*Nmodes],:)=0;
%         Adotx([2,2+3*Nmodes],:)=0;
%         Adotx([3,3+3*Nmodes],:)=0;
%         Adotx(2+3*Nmodes,2+3*Nmodes)=1; 
%         Adotx(3+3*Nmodes,3+3*Nmodes)=1; 
%         Adotx(2,2)=1; 
%         Adotx(3,3)=1;
%         Ax(2,2)=1; 
%         Ax(3,3)=1; 
%         Ax(2+3*Nmodes,2+3*Nmodes)=1; 
%         Ax(3+3*Nmodes,3+3*Nmodes)=1;
%         %test
%         Ax(6*Nmodes+1,:)=0;
%         Ax(6*Nmodes+1,1)=-4*pi*Gg;
%     end
%     %% (2) Propgate solution 
%     Aprop=Adotx\Ax;
% end

function Aprop=get_Aprop(rK,gK,dgK,Nmodes,A1,A2,A3_inv,A4,A5,A6,A71,A72,A81,A82,A9,A100,A101,A102,A11,A12,A13,deg_0,Gg)
%% (1) Assemble matrix
    % rheology equation 
    Adotx1=zeros(3*Nmodes,8*Nmodes);
    Ax1=zeros(3*Nmodes,8*Nmodes);
    Adotx1(:,1:3*Nmodes)=A4*A1*A3_inv; %\dot{U}
    Ax1(:,1:3*Nmodes)=-A4*A2*A3_inv/rK;% U
    Ax1(:,3*Nmodes+(1:3*Nmodes))=A13;% \Sigma
    % for uniform model this looks OK. 
    % momentum equations 
    Adotx2=zeros(3*Nmodes,8*Nmodes);
    Ax2=zeros(3*Nmodes,8*Nmodes);
    Adotx2(:,1:3*Nmodes)=-A5*A1*A3_inv/rK+A6;% \dot{U}
    Adotx2(:,3*Nmodes+(1:3*Nmodes))=A13;% \dot{\Sigma{U}}
    Ax2(:,1:3*Nmodes)=A5*A2*A3_inv/rK^2+gK/rK*A71+dgK*A72;% U
    Ax2(:,6*Nmodes+(1:2*Nmodes))=A81+A82/rK;% \Phi
    % for uniform model, this looks OK 
    % Poisson equation 
    Adotx3=zeros(2*Nmodes,8*Nmodes);
    Ax3=zeros(2*Nmodes,8*Nmodes);
    Adotx3(:,1:3*Nmodes)=-A12; %\dot{U}
    Adotx3(:,6*Nmodes+(1:2*Nmodes))=A9; %\dot{\Phi}
    Ax3(:,1:3*Nmodes)=A11/rK; % U
    Ax3(:,6*Nmodes+(1:2*Nmodes))=A100+A101/rK+A102/rK^2; % \Phi
    % combine matrices
    Adotx=[Adotx1; Adotx2; Adotx3];
    Ax=[Ax1; Ax2; Ax3];
    if deg_0==1
        Ax([2,2+3*Nmodes],:)=0; 
        Ax([3,3+3*Nmodes],:)=0;
        Adotx([2,2+3*Nmodes],:)=0;
        Adotx([3,3+3*Nmodes],:)=0;
        Adotx(2+3*Nmodes,2+3*Nmodes)=1; 
        Adotx(3+3*Nmodes,3+3*Nmodes)=1; 
        Adotx(2,2)=1; 
        Adotx(3,3)=1;
        Ax(2,2)=1; 
        Ax(3,3)=1; 
        Ax(2+3*Nmodes,2+3*Nmodes)=1; 
        Ax(3+3*Nmodes,3+3*Nmodes)=1;
        %test
        Ax(6*Nmodes+1,:)=0;
        Ax(6*Nmodes+1,1)=-4*pi*Gg;
    end
    %% (2) Propgate solution 
    Aprop=Adotx\Ax;
end