%% OBTAIN_SOLUTION
%AUTHOR: M. Rovira-Navarro 
%USE: Obtain solution using propagator method

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
function [y_sol] = get_solution(Interior_Model,Forcing,Numerics,Couplings,verbose,out_file_name)
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

ind=find(Couplings.n_s==0);
if isempty(ind)==1
    deg0=0;
else
    deg0=1; 
end
% check if there is a subsurface ocean
ocean_layer=0;
ocean_flag=0; 
if isfield(Interior_Model,'ocean')
    for i=1:length(Interior_Model)
        if Interior_Model(i).ocean==1
            ocean_layer=i; 
            ocean_start=Numerics.BCindices(ocean_layer-2);
            ocean_end=Numerics.BCindices(ocean_layer-1);
        end
    end
end
% Initialise the solution matrix
Nmodes=length(Couplings.n_s);
y1=zeros(8*Nmodes,8*Nmodes,Numerics.Nr+1); %from core to ocean 
y2=zeros(8*Nmodes,8*Nmodes,Numerics.Nr+1); %from ocean to outer shell 
y3=zeros(8*Nmodes,8*Nmodes,Numerics.Nr+1); %outer shell
y_old=zeros(8*Nmodes,8*Nmodes,1);
for i=1:8*Nmodes
    y1(i,i,1)=1; %set Cs.
    if ocean_layer>0
        y2(i,i,ocean_start)=1; %set Cs to 1
        y3(i,i,ocean_end)=1; %set Cs to 1
    end
end


% Set the indice of the first boundary layer if applicable
blayer_ind = 0; % Set default value
iblayer = 1; % First indices of Numerics.BCindices
k_corr = 1; % Indices correction
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

% Fill matrices that do depend on the rheology using the rheology of the first layer
[A2,A1]=get_A1A2(Interior_Model(ilayer),Couplings); 
[A13, A6, A71, A72, A81, A82, A9, A100, A101, A102, A11, A12]=get_others(Couplings,Interior_Model(ilayer));

if verbose==1
    disp(' ###########################################################################')
    disp(' Integrating radially ...')
end
%% INTEGRATE RADIALLY
tic
% Start the radial integration
for k=2:Numerics.Nr+1    
    % Check whether the integration is at a new layer & change to the y solution needed depending if we are below or above the ocean
    if k==2
        y=y1;
    end
    if k == blayer_ind 
        % Go to next layer
        ilayer = ilayer + 1;
        % Update density 
        rhoK = Interior_Model(ilayer).rho;
        % Update mass and radius of the encapsulated sphere
        Min = Min + 4/3*pi*Interior_Model(ilayer-1).rho*(Interior_Model(ilayer-1).R^3 - Interior_Model(ilayer-2).R^3);
        Rin = Interior_Model(ilayer-1).R;
        % Update delta r
        Delta_r = (Interior_Model(ilayer).R - Interior_Model(ilayer-1).R)/(Numerics.Nrlayer(ilayer));
         % If there is another layer reprime blayer_ind
        if iblayer < length(Numerics.BCindices)
            iblayer = iblayer + 1;
            blayer_ind = Numerics.BCindices(iblayer) + 1;
        else
            % Update regardless to keep consistency in last layer
            iblayer = iblayer + 1;
        end 
        % Update indices correction, needed to 
        k_corr = Numerics.BCindices(iblayer-1);
        r(k) = Rin + (k-k_corr)*Delta_r;
        % Obtain density difference with previous layer (is positive)
        rho_diff = Interior_Model(ilayer).Delta_rho;
        if ilayer==ocean_layer % we are at the ocean layer
            y1=y;%store propagation of inner layer 
            y=y2;%get ready to propagate ocean layer
            y_old=y(:,:,k-1); %set y_old 
            % Recalculate propagation matrix
            [A13, A6, A71, A72, A81, A82, A9, A100, A101, A102, A11, A12]=get_others(Couplings,Interior_Model(ilayer));
            % flag that we are in an ocean layer 
            ocean_flag=1; 
        elseif ilayer==ocean_layer+1 %we are in the first layer of the shell 
            ocean_flag=0; 
            y2=y;% store ocean propagation
            y=y3;%get ready to propagate ice layer
            y_old=y(:,:,k-1); %set y_old 
            % Recalculate propagation matrix 
            [A2,A1]=get_A1A2(Interior_Model(ilayer),Couplings); 
            [A13, A6, A71, A72, A81, A82, A9, A100, A101, A102, A11, A12]=get_others(Couplings,Interior_Model(ilayer));
        else %just another layer 
            % Recalculate propagation matrix
             ocean_flag=0;
            [A2,A1]=get_A1A2(Interior_Model(ilayer),Couplings); 
            [A13, A6, A71, A72, A81, A82, A9, A100, A101, A102, A11, A12]=get_others(Couplings,Interior_Model(ilayer));
             % Apply the effect of the density discontinuity
            cont_condition(8*(0:(Nmodes-1))+8,:,1) = ...
            4*pi*Gg*rho_diff*y(8*(0:(Nmodes-1))+1,:,k-1);
            y_old=y(:,:,k-1) + cont_condition;
        end       
    else
        % Update the solution at the previous node
        y_old = y(:,:,k-1);
        % Update radial points vector
        r(k) = Rin + (k-k_corr)*Delta_r;
    end
    %%%%%%%%%%% Step 1
    rK = Rin + ((k-k_corr)-1)*Delta_r;
    gK = Gg*(Min + 4/3*pi*rhoK*(rK^3 - Rin^3))/rK^2;
    dgK = Gg*(2*(4/3*pi*rhoK*Rin^3 - Min)/rK^3 + 4/3*pi*rhoK);
    Aprop = get_Aprop(rK,gK,dgK,Nmodes,A1,A2,A3_inv,A4,A5,A6,A71,A72,A81,A82,A9,A100,A101,A102,A11,A12,A13,deg0,Gg,ocean_flag);
    k1 = Delta_r*Aprop*y_old;
    %%%%%%%%%%% Step 2
    rK = Rin + ((k-k_corr)-1)*Delta_r + AA2*Delta_r; 
    gK = Gg*(Min + 4/3*pi*rhoK*(rK^3 - Rin^3))/rK^2;
    dgK = Gg*(2*(4/3*pi*rhoK*Rin^3 - Min)/rK^3 + 4/3*pi*rhoK);
    Aprop = get_Aprop(rK,gK,dgK,Nmodes,A1,A2,A3_inv,A4,A5,A6,A71,A72,A81,A82,A9,A100,A101,A102,A11,A12,A13,deg0,Gg,ocean_flag);
    k2 = Delta_r*Aprop*(y_old + B21*k1);
    %%%%%%%%%%% Step 3
    rK = Rin + ((k-k_corr)-1)*Delta_r + AA3*Delta_r; 
    gK = Gg*(Min + 4/3*pi*rhoK*(rK^3 - Rin^3))/rK^2;
    dgK = Gg*(2*(4/3*pi*rhoK*Rin^3 - Min)/rK^3 + 4/3*pi*rhoK);
    Aprop = get_Aprop(rK,gK,dgK,Nmodes,A1,A2,A3_inv,A4,A5,A6,A71,A72,A81,A82,A9,A100,A101,A102,A11,A12,A13,deg0,Gg,ocean_flag);
    k3 = Delta_r*Aprop*(y_old + B31*k1 + B32*k2);
    %%%%%%%%%%% Step 4
    rK = Rin + ((k-k_corr)-1)*Delta_r + AA4*Delta_r; 
    gK = Gg*(Min + 4/3*pi*rhoK*(rK^3 - Rin^3))/rK^2;
    dgK = Gg*(2*(4/3*pi*rhoK*Rin^3 - Min)/rK^3 + 4/3*pi*rhoK);
    Aprop = get_Aprop(rK,gK,dgK,Nmodes,A1,A2,A3_inv,A4,A5,A6,A71,A72,A81,A82,A9,A100,A101,A102,A11,A12,A13,deg0,Gg,ocean_flag);
    k4 = Delta_r*Aprop*(y_old + B41*k1 + B42*k2 + B43*k3);
    %%%%%%%%%%% Step 5
    rK = Rin + ((k-k_corr)-1)*Delta_r + AA5*Delta_r; 
    gK = Gg*(Min + 4/3*pi*rhoK*(rK^3 - Rin^3))/rK^2;
    dgK = Gg*(2*(4/3*pi*rhoK*Rin^3 - Min)/rK^3 + 4/3*pi*rhoK);
    Aprop = get_Aprop(rK,gK,dgK,Nmodes,A1,A2,A3_inv,A4,A5,A6,A71,A72,A81,A82,A9,A100,A101,A102,A11,A12,A13,deg0,Gg,ocean_flag);
    k5 = Delta_r*Aprop*(y_old + B51*k1 + B52*k2 + B53*k3 + B54*k4);
    %%%%%%%%%%% Step 6
    rK = Rin + ((k-k_corr)-1)*Delta_r + AA6*Delta_r; 
    gK = Gg*(Min + 4/3*pi*rhoK*(rK^3 - Rin^3))/rK^2;
    dgK = Gg*(2*(4/3*pi*rhoK*Rin^3 - Min)/rK^3 + 4/3*pi*rhoK);
    Aprop = get_Aprop(rK,gK,dgK,Nmodes,A1,A2,A3_inv,A4,A5,A6,A71,A72,A81,A82,A9,A100,A101,A102,A11,A12,A13,deg0,Gg,ocean_flag);
    k6 = Delta_r*Aprop*(y_old + B61*k1 + B62*k2 + B63*k3 + B64*k4 + B65*k5);
    %get next step
    y(:,:,k) = y_old + AC1*k1 + AC2*k2 + AC3*k3 + AC4*k4 + AC5*k5 + AC6*k6;
end
if ocean_layer>0
    y3=y; 
else
    y1=y; 
end
if verbose == 1
    disp(['Time spent ingtegrating radially ' num2str(toc) ' s'])
end

%% REARRANGE SOLUTION VECTOR 
tic
y1_re=zeros(size(y));
y1_re2=zeros(size(y));
y2_re=zeros(size(y));
y2_re2=zeros(size(y));
y3_re=zeros(size(y));
y3_re2=zeros(size(y));
for i=1:Nmodes
    index1=3*(i-1)+(1:3);
    index2=3*Nmodes+3*(i-1)+(1:3);
    index3=6*Nmodes+2*(i-1)+(1:2);
    y1_re(8*(i-1)+(1:8),:,:)=y1([index1 index2 index3],:,:);
    y2_re(8*(i-1)+(1:8),:,:)=y2([index1 index2 index3],:,:);
    y3_re(8*(i-1)+(1:8),:,:)=y3([index1 index2 index3],:,:);
end
for i=1:Nmodes
    index1=3*(i-1)+(1:3);
    index2=3*Nmodes+3*(i-1)+(1:3);
    index3=6*Nmodes+2*(i-1)+(1:2);
    y1_re2(:,8*(i-1)+(1:8),:)=y1_re(:,[index1 index2 index3],:); 
    y2_re2(:,8*(i-1)+(1:8),:)=y2_re(:,[index1 index2 index3],:); 
    y3_re2(:,8*(i-1)+(1:8),:)=y3_re(:,[index1 index2 index3],:); 
end
if verbose == 1
    disp(['Time spent rearranging the solution vector' num2str(toc) ' s'])
end
%% ASSEMBLE MATRIX FOR INVERSION
if ocean_layer>0 %there is an ocean
    % there are 3*8*Nmodes integration constants
    B=zeros(24*Nmodes,24*Nmodes);
    B2=zeros(24*Nmodes,1);
    y1=y1_re2;
    y2=y2_re2;
    y3=y3_re2;
else %there is not an ocean 
    B=zeros(8*Nmodes,8*Nmodes); %there are 8*Nmodes integration constants
    B2=zeros(8*Nmodes,1);
    y1=y1_re2;
    y3=y1_re2;
end
tic
if verbose==1
    disp(' ###########################################################################')
    disp(' Applying BC to get the integration constants ...')
end

if ocean_layer==0 % the model does not have an ocean layer 
for i=1:Nmodes
    n=Couplings.n_s(i);
    m=Couplings.m_s(i);
    % Written such that the code handles icy moons correctly
    rho2 = Interior_Model(1).Delta_rho + Interior_Model(2).rho; 
    for j=1:8*Nmodes
        %Core-Mantle BC -----------------------
        % BC1 radial stress
        if n==-3
        B(8*(i-1)+1,j)=B(8*(i-1)+1,j)+y1(8*(i-1)+4,j,1); %R=0
        else
        B(8*(i-1)+1,j)=B(8*(i-1)+1,j)+y1(8*(i-1)+1,j,1); %U
        B(8*(i-1)+1,j)=B(8*(i-1)+1,j)-1/gc/rho2*y1(8*(i-1)+4,j,1); %R 
        B(8*(i-1)+1,j)=B(8*(i-1)+1,j)+1/gc*y1(8*(i-1)+7,j,1); %\Phi
        end
        % BC2 no tangential stress
        B(8*(i-1)+2,j)=B(8*(i-1)+2,j)+y1(8*(i-1)+5,j,1);  % S 
        % BC3 no toroidal stress 
        if n==1
            B(8*(i-1)+3,j)=B(8*(i-1)+3,j)+y1(8*(i-1)+3,j,1);   %W
        else
            B(8*(i-1)+3,j)=B(8*(i-1)+3,j)+y1(8*(i-1)+6,j,1);   %T 
        end
        % BC4 potential stress
        if n==0
            B(8*(i-1)+4,j)=B(8*(i-1)+4,j)+y1(8*(i-1)+8,j,end); %k=0
        else
            fac=4*pi*Gg/gc/rho2*(Interior_Model(2).rho - rho2);
            B(8*(i-1)+4,j)=B(8*(i-1)+4,j)-fac*y1(8*(i-1)+4,j,1); %R
            B(8*(i-1)+4,j)=B(8*(i-1)+4,j)+(n/Rc + fac*rho2)*y1(8*(i-1)+7,j,1); %\Phi
            B(8*(i-1)+4,j)= B(8*(i-1)+4,j)-y1(8*(i-1)+8,j,1);  %\dot\Phi
        end
        % Surface BC ----------------------- 
        % BC5 no radial stress
        B(8*(i-1)+5,j)=B(8*(i-1)+5,j)+y3(8*(i-1)+4,j,end); %R=0
        % BC6 no tangential stress
        if n==0
            B(8*(i-1)+6,j)=B(8*(i-1)+6,j)+y3(8*(i-1)+2,j,end); %V=0
        else
            B(8*(i-1)+6,j)=B(8*(i-1)+6,j)+y3(8*(i-1)+5,j,end); %S=0
        end
        % BC7 no toroidal stress 
        if n==0
            B(8*(i-1)+7,j)=B(8*(i-1)+7,j)+y3(8*(i-1)+3,j,end); %W=0 
        else
            B(8*(i-1)+7,j)=B(8*(i-1)+7,j)+y3(8*(i-1)+6,j,end); %T=0  
        end
        % BC8 potential stress 
        if n<2
            B(8*(i-1)+8,j)=B(8*(i-1)+8,j)+y3(8*(i-1)+7,j,end); %k=0
        else
            B(8*(i-1)+8,j)=B(8*(i-1)+8,j)+4*pi*Gg*rhoK*y3(8*(i-1)+1,j,end); %U
            B(8*(i-1)+8,j)=B(8*(i-1)+8,j)+(n+1)*y3(8*(i-1)+7,j,end); %\Phi
            B(8*(i-1)+8,j)=B(8*(i-1)+8,j)+y3(8*(i-1)+8,j,end); %\dot\Phi
        end
    end
    if n==Forcing.n && m==Forcing.m
        B2(8*(i-1)+8,1)=2*n+1;
    end
end
else %----------------------------------------- There is a subsurface ocean!
 for i=1:Nmodes
    n=Couplings.n_s(i);
    m=Couplings.m_s(i);
    % Written such that the code handles icy moons correctly
    rho2 = Interior_Model(1).Delta_rho + Interior_Model(2).rho; 
    for j=1:8*Nmodes
        %Core-Mantle BC -----------------------
        % BC1 radial stress
        if n==-3
        B(8*(i-1)+1,j)=B(8*(i-1)+1,j)+y1(8*(i-1)+4,j,1); %R=0
        else
        B(8*(i-1)+1,j)=B(8*(i-1)+1,j)+y1(8*(i-1)+1,j,1); %U
        B(8*(i-1)+1,j)=B(8*(i-1)+1,j)-1/gc/rho2*y1(8*(i-1)+4,j,1); %R 
        B(8*(i-1)+1,j)=B(8*(i-1)+1,j)+1/gc*y1(8*(i-1)+7,j,1); %\Phi
        end
        % BC2 no tangential stress
        B(8*(i-1)+2,j)=B(8*(i-1)+2,j)+y1(8*(i-1)+5,j,1);  % S 
        % BC3 no toroidal stress 
        if n==1
            B(8*(i-1)+3,j)=B(8*(i-1)+3,j)+y1(8*(i-1)+3,j,1);   %W
        else
            B(8*(i-1)+3,j)=B(8*(i-1)+3,j)+y1(8*(i-1)+6,j,1);   %T 
        end
        % BC4 potential stress
        if n==0
            B(8*(i-1)+4,j)=B(8*(i-1)+4,j)+y1(8*(i-1)+8,j,end); %k=0
        else
            fac=4*pi*Gg/gc/rho2*(Interior_Model(2).rho - rho2);
            B(8*(i-1)+4,j)=B(8*(i-1)+4,j)-fac*y1(8*(i-1)+4,j,1); %R
            B(8*(i-1)+4,j)=B(8*(i-1)+4,j)+(n/Rc + fac*rho2)*y1(8*(i-1)+7,j,1); %\Phi
            B(8*(i-1)+4,j)= B(8*(i-1)+4,j)-y1(8*(i-1)+8,j,1);  %\dot\Phi
        end
        % Surface BC ----------------------- 
        % BC5 no radial stress
        B(8*(i-1)+5,16*Nmodes+j)= B(8*(i-1)+5,16*Nmodes+j)+y3(8*(i-1)+4,j,end); %R=0
        % BC6 no tangential stress
        if n==0
            B(8*(i-1)+6,16*Nmodes+j)=B(8*(i-1)+6,16*Nmodes+j)+y3(8*(i-1)+2,j,end); %V=0
        else
            B(8*(i-1)+6,16*Nmodes+j)=B(8*(i-1)+6,16*Nmodes+j)+y3(8*(i-1)+5,j,end); %S=0
        end
        % BC7 no toroidal stress 
        if n==0
            B(8*(i-1)+7,16*Nmodes+j)= B(8*(i-1)+7,16*Nmodes+j)+y3(8*(i-1)+3,j,end); %W=0 
        else
            B(8*(i-1)+7,16*Nmodes+j)=B(8*(i-1)+7,16*Nmodes+j)+y3(8*(i-1)+6,j,end); %T=0  
        end
        % BC8 potential stress 
        if n<2
            B(8*(i-1)+8,16*Nmodes+j)=y3(8*(i-1)+7,j,end); %k=0
        else
            B(8*(i-1)+8,16*Nmodes+j)=B(8*(i-1)+8,16*Nmodes+j)+4*pi*Gg*rhoK*y3(8*(i-1)+1,j,end); %U
            B(8*(i-1)+8,16*Nmodes+j)=B(8*(i-1)+8,16*Nmodes+j)+(n+1)*y3(8*(i-1)+7,j,end); %\Phi
            B(8*(i-1)+8,16*Nmodes+j)=B(8*(i-1)+8,16*Nmodes+j)+y3(8*(i-1)+8,j,end); %\dot\Phi
        end
        %BC in case there is an ocean 
        % Mantle-ocean BC --------------
        gO=Interior_Model(ocean_layer-1).gs; 
        gI=Interior_Model(ocean_layer).gs;
        rhoO=Interior_Model(ocean_layer).rho; 
        % BC9 radial stress
        B(8*(i-1)+9,j)=y1(8*(i-1)+1,j,ocean_start)-1/gO/rhoO*y1(8*(i-1)+4,j,ocean_start)+1/gO*y1(8*(i-1)+7,j,ocean_start); %      
        % BC10 no tangential stress
        B(8*(i-1)+10,j)=y1(8*(i-1)+5,j,ocean_start); % S=0
        % BC11 no toroidal stress
        B(8*(i-1)+11,j)=y1(8*(i-1)+6,j,ocean_start); % T=0
        % BC12 continuity potential 
        B(8*(i-1)+12,j)=y1(8*(i-1)+7,j,ocean_start); % phi below ocean
        B(8*(i-1)+12,8*Nmodes+j)=-y2(8*(i-1)+7,j,ocean_start); % phi ocean
        % BC13 potential stres
        B(8*(i-1)+13,j)=y1(8*(i-1)+8,j,ocean_start)...
                        +4*pi*Gg*(Interior_Model(ocean_layer-1).rho-Interior_Model(ocean_layer).rho)*y1(8*(i-1)+1,j,ocean_start); %dphi below ocean
        B(8*(i-1)+13,8*Nmodes+j)=-y2(8*(i-1)+8,j,ocean_start); %dphi at ocean
        % additional BC bc some things are not defined in the ocean 
        B(8*(i-1)+14,8*Nmodes+j)=y2(8*(i-1)+1,j,ocean_start);
        B(8*(i-1)+15,8*Nmodes+j)=y2(8*(i-1)+2,j,ocean_start);
        B(8*(i-1)+16,8*Nmodes+j)=y2(8*(i-1)+3,j,ocean_start);
        B(8*(i-1)+17,8*Nmodes+j)=y2(8*(i-1)+4,j,ocean_start);
        B(8*(i-1)+18,8*Nmodes+j)=y2(8*(i-1)+5,j,ocean_start);
        B(8*(i-1)+19,8*Nmodes+j)=y2(8*(i-1)+6,j,ocean_start);
        % Ocean-shell BC -----------------
        % BC20 radial stress
        B(8*(i-1)+20,16*Nmodes+j)=y3(8*(i-1)+1,j,ocean_end)-1/gI/rhoO*y3(8*(i-1)+4,j,ocean_end)+1/gI*y3(8*(i-1)+7,j,ocean_end); % U+\phi/g-R/\rho g       
        % BC21 no tangential stress
        B(8*(i-1)+21,16*Nmodes+j)=y3(8*(i-1)+5,j,ocean_end); % S=0
        % BC22 no toroidal stress
        B(8*(i-1)+22,16*Nmodes+j)=y3(8*(i-1)+6,j,ocean_end); % T=0
        % BC23 continuity potential 
        B(8*(i-1)+23,16*Nmodes+j)=y3(8*(i-1)+7,j,ocean_end); % phi above ocean
        B(8*(i-1)+23,8*Nmodes+j)=-y2(8*(i-1)+7,j,ocean_end); % phi ocean
        % BC24 potential stres
        B(8*(i-1)+24,16*Nmodes+j)=y3(8*(i-1)+8,j,ocean_end)...
                                +4*pi*Gg*(Interior_Model(ocean_layer+1).rho-Interior_Model(ocean_layer).rho)*y3(8*(i-1)+1,j,ocean_end); % U at the boundary 
        B(8*(i-1)+24,8*Nmodes+j)=-y2(8*(i-1)+8,j,ocean_end); %dphi at ocean
    end
    if n==Forcing.n && m==Forcing.m
        B2(8*(i-1)+8,1)=2*n+1;
    end
end
end

if verbose == 1
    disp(['Time spent assembling BC matrix' num2str(toc) ' s'])
end

%% OBTAIN COEFFICIENTS
tic
C=B\B2;
if verbose == 1
    disp(['Time spent obtaining integration constants' num2str(toc) ' s'])
end

%% ASSEMBLE SOLUTION
y_sol = zeros(Numerics.Nr+1,8*Nmodes);
if ocean_layer==0
    for i=1:Numerics.Nr+1
        y_sol(i,:)=transpose(y1(:,:,i)*C);
    end
else % there is an ocean and the three solutions need to be combined. 
    y2(:,:,ocean_start)=0; % point at ocean_start and at ocean_end belongs to layer below and above the ocean
    y2(:,:,ocean_end)=0; 
    for i=1:Numerics.Nr+1
        y_sol(i,:)=transpose(y1(:,:,i)*C(1:8*Nmodes,1)+y2(:,:,i)*C(8*Nmodes+1:16*Nmodes,1)+y3(:,:,i)*C(16*Nmodes+1:24*Nmodes,1));
    end
end

%% Check BC
if verbose==1
    for i=1:Nmodes
        n=Couplings.n_s(i);
        mode=8*(i-1); 
        disp(['Mode (n,m) (' num2str(Couplings.n_s(i)) ',' num2str(Couplings.m_s(i)) ')'])
        disp(['CMB:' num2str(y_sol(1,mode+1)+y_sol(1,mode+7)/gc-y_sol(1,mode+4)/rho2/gc) ', ' num2str(y_sol(1,mode+5)) ',' num2str(y_sol(1,mode+6)) ', ' num2str(-y_sol(1,mode+8)+(n/Rc+fac*rho2)*y_sol(1,mode+7)-fac*y_sol(1,mode+4))])
        disp(['Surface: ' num2str(y_sol(end,mode+4)) ', ' num2str(y_sol(end,mode+5)) ', ' num2str(y_sol(end,mode+6)) ', ' num2str(y_sol(end,mode+8)+(n+1)*y_sol(end,mode+7)+4*pi*Gg*rhoK*y_sol(end,mode+1))])
    end
end

%% COMPUTE AUXILIARY VARIABLES
tic
i_U=[];
for i=1:Nmodes
    i_U=[i_U 8*(i-1)+(1:3)];
end
U=y_sol(:,i_U);
sigma=zeros(Numerics.Nr+1,6*Nmodes);
epsilon=zeros(Numerics.Nr+1,6*Nmodes);
% compute u 
u=transpose(A3_inv*transpose(U));
%compute u_dot (numerically)
index1=[];
index2=[];
index3=[];
for i=1:Nmodes
    index1=[index1 8*(i-1)+(1:3)];
    index2=[index2 8*(i-1)+(4:6)];
    index3=[index3 8*(i-1)+(7:8)];
end
u_dot=zeros(Numerics.Nr+1,3*Nmodes);
u_aux=zeros(size(u));
u_aux(1:end-1,:)=u(1:end-1,:);

%%%%%%%%%%% Reset all auxiliary variables %%%%%%%%%%%%
% Set the indice of the first boundary layer if applicable
blayer_ind = 0; % Set default value
iblayer = 1; % First indices of Numerics.BCindices
if ~isempty(Numerics.BCindices)
    % Add 1 to make sure that the Boundary is at the previous node
    blayer_ind = Numerics.BCindices(iblayer) + 1;
end

% Start counter of layer number
ilayer = 2; 
% Reset some necessary values
Rin = Interior_Model(1).R; % Radius of the inner sphere (can change)
Min = 4/3*pi*rhoC*Rc^3; % Mass of the inner sphere (can change)
rhoK = Interior_Model(ilayer).rho; % Density of the current layer
% Fill in matrices forthe first layer
[A2,A1]=get_A1A2(Interior_Model(ilayer),Couplings); 
[A13, A6, A71, A72, A81, A82, A9, A100, A101, A102, A11, A12]=get_others(Couplings,Interior_Model(ilayer));
for k=1:Numerics.Nr
    % Check whether the integration is at a new layer
    if k == blayer_ind 
        % Go to next layer
        ilayer = ilayer + 1;
        if ilayer==ocean_layer
            ocean_flag=1;
            [A13, A6, A71, A72, A81, A82, A9, A100, A101, A102, A11, A12]=get_others(Couplings,Interior_Model(ilayer));
        else
            ocean_flag=0;
            [A2,A1]=get_A1A2(Interior_Model(ilayer),Couplings); 
            [A13, A6, A71, A72, A81, A82, A9, A100, A101, A102, A11, A12]=get_others(Couplings,Interior_Model(ilayer));
        end
        % If there is another layer reprime blayer_ind
        if iblayer < length(Numerics.BCindices)
            iblayer = iblayer + 1;
            blayer_ind = Numerics.BCindices(iblayer) + 1;
        end
        % Update density 
        rhoK = Interior_Model(ilayer).rho;
        % Update mass and radius of the encapsulated sphere
        Min = Min + 4/3*pi*Interior_Model(ilayer-1).rho*(Interior_Model(ilayer-1).R^3 - Interior_Model(ilayer-2).R^3);
        Rin = Interior_Model(ilayer-1).R;
    end
    rK = r(k); 
    gK = Gg*(Min + 4/3*pi*rhoK*(rK^3 - Rin^3))/rK^2;
    dgK = Gg*(2*(4/3*pi*rhoK*Rin^3 - Min)/rK^3 + 4/3*pi*rhoK);
    Aprop = get_Aprop(rK,gK,dgK,Nmodes,A1,A2,A3_inv,A4,A5,A6,A71,A72,A81,A82,A9,A100,A101,A102,A11,A12,A13,deg0,Gg,ocean_flag);
    x_dot = Aprop(1:3*Nmodes,:)*transpose(y_sol(k,[index1 index2 index3])); 
    u_dot(k,:) = transpose(A3_inv*x_dot); 
    sigma(k,:) = transpose(A1*transpose(u_dot(k,:)) + A2*transpose(u_aux(k,:)) / r(k));
end
epsilon(:,:)=transpose(A14*transpose(u_dot)+A15*transpose(u_aux)./r);
%compue u_dot (numerically)
% u_dot=diff(u)/Delta_r; 
% r_av=(r(1:end-1)'+r(2:end)')/2;
% u_av=(u(1:end-1,:)+u(2:end,:))/2;
% sigma(2:end,:)=transpose(A1*transpose(u_dot)+A2*transpose(u_av)./r_av');
% epsilon(2:end,:)=transpose(A14*transpose(u_dot)+A15*transpose(u_av)./r_av');
if verbose == 1
    disp(['Time Spent auxilliary variables' num2str(toc) ' s'])
end

%% REARRANGE
tic
y_sol2=y_sol;
y_sol=zeros(Numerics.Nr+1,24,Nmodes);
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
if verbose == 1
    disp(['Time Spent rearranging the solurtion' num2str(toc) ' s'])
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
        A1(6*(i-1)+1,3*(i-1)+2)=0;
        A1(6*(i-1)+1,3*(i-1)+3)=(3*lambda+2*mu)/sqrt(3)/sqrt(2*n+1)*sqrt(n+1)*(n+2);  
        %\dot u
        A2(6*(i-1)+1,3*(i-1)+1)=-(3*lambda+2*mu)/sqrt(3)/sqrt(2*n+1)*sqrt(n);
        A2(6*(i-1)+1,3*(i-1)+2)=0;
        A2(6*(i-1)+1,3*(i-1)+3)=+(3*lambda+2*mu)/sqrt(3)/sqrt(2*n+1)*sqrt(n+1);    
        %\sigma_{n,n-2,2}----------------------------
        %u
        A1(6*(i-1)+2,3*(i-1)+1)=2*mu*sqrt((n-1)/(2*n-1))*n;
        A1(6*(i-1)+2,3*(i-1)+2)=0;
        A1(6*(i-1)+2,3*(i-1)+3)=0;   
        %\dot u
        A2(6*(i-1)+2,3*(i-1)+1)=2*mu*sqrt((n-1)/(2*n-1));
        A2(6*(i-1)+2,3*(i-1)+2)=0;
        A2(6*(i-1)+2,3*(i-1)+3)=0;    
        %\sigma_{n,n-1,2}-----------------------------
        %u
        A1(6*(i-1)+3,3*(i-1)+1)=0;
        A1(6*(i-1)+3,3*(i-1)+2)=2*mu/sqrt(2)*sqrt((n-1)/(2*n+1))*(n+1);
        A1(6*(i-1)+3,3*(i-1)+3)=0;    
        %\dot u
        A2(6*(i-1)+3,3*(i-1)+1)=0;
        A2(6*(i-1)+3,3*(i-1)+2)=2*mu/sqrt(2)*sqrt((n-1)/(2*n+1));
        A2(6*(i-1)+3,3*(i-1)+3)=0;  
        %\sigma_{n,n,2}-------------------------------
        %u
        A1(6*(i-1)+4,3*(i-1)+1)=2*mu*sqrt((2*n+3)*(2*n+2)/12/(2*n-1)/(2*n+1))*(n-1);
        A1(6*(i-1)+4,3*(i-1)+2)=0;
        A1(6*(i-1)+4,3*(i-1)+3)=2*mu*sqrt(n*(2*n-1)*(n+1)/3/(2*n+3)/(2*n+2)/(2*n+1))*(n+2);    
        %\dot u
        A2(6*(i-1)+4,3*(i-1)+1)=-2*mu*sqrt((2*n+3)*(2*n+2)/12/(2*n-1)/(2*n+1));
        A2(6*(i-1)+4,3*(i-1)+2)=0;
        A2(6*(i-1)+4,3*(i-1)+3)=2*mu*sqrt(n*(2*n-1)*(n+1)/3/(2*n+3)/(2*n+2)/(2*n+1));   
        %\sigma_{n,n+1,2}-----------------------------
        %u
        A1(6*(i-1)+5,3*(i-1)+1)=0;
        A1(6*(i-1)+5,3*(i-1)+2)=2*mu/sqrt(2)*sqrt((n+2)/(2*n+1))*n;
        A1(6*(i-1)+5,3*(i-1)+3)=0;    
        %\dot u
        A2(6*(i-1)+5,3*(i-1)+1)=0;
        A2(6*(i-1)+5,3*(i-1)+2)=-2*mu/sqrt(2)*sqrt((n+2)/(2*n+1));
        A2(6*(i-1)+5,3*(i-1)+3)=0;
        %\sigma_{n,n+2,2}-----------------------------
        %u
        A1(6*(i-1)+6,3*(i-1)+1)=0;
        A1(6*(i-1)+6,3*(i-1)+2)=0;
        A1(6*(i-1)+6,3*(i-1)+3)=2*mu*sqrt((n+2)/(2*n+3))*(n+1);    
        %\dot u
        A2(6*(i-1)+6,3*(i-1)+1)=0;
        A2(6*(i-1)+6,3*(i-1)+2)=0;
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
    % find couplings that appear in these equations
%     icou = zeros(Nreo,1);
%     for ireo=1:Nreo
%         icou{ireo}=find(Coup(ieq,:,27,ireo)>0);
%     end
    for ireo=1:Nreo % Nreo loop over rheology 
        ialpha=find(Coup(ieq,:,27,ireo)>0);%n_s(ialpha) m_s(ialpha) are the degrees and orders that appear in the equation
        K_nm=variable_rheology(ireo,3);
        mu_nm=variable_rheology(ireo,4);
%         ialpha=icou{ireo}; %n_s(ialpha) m_s(ialpha) are the degrees and orders that appear in the equation
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
            A1(6*(i-1)+1,3*(ia-1)+2)=A1(6*(i-1)+1,3*(ia-1)+2)+...
                                    0;
            A1(6*(i-1)+1,3*(ia-1)+3)=A1(6*(i-1)+1,3*(ia-1)+3)+...
                                    K_nm/sqrt(3)/sqrt(2*na+1)*sqrt(na+1)*(na+2)*Cp;
            %\dot u
            A2(6*(i-1)+1,3*(ia-1)+1)=A2(6*(i-1)+1,3*(ia-1)+1)+...
                                    -K_nm/sqrt(3)/sqrt(2*na+1)*sqrt(na)*Cp;
            A2(6*(i-1)+1,3*(ia-1)+2)=A2(6*(i-1)+1,3*(ia-1)+2)+...
                                    0;
            A2(6*(i-1)+1,3*(ia-1)+3)=A2(6*(i-1)+1,3*(ia-1)+3)+...
                                    +K_nm/sqrt(3)/sqrt(2*na+1)*sqrt(na+1)*Cp;

            %\sigma_{n,n-2,2}----------------------------
            comp=2;% couplings 7-11
            %n--------
            Cp=Coup(ieq,ia,7,ireo); 
            %u
            A1(6*(i-1)+comp,3*(ia-1)+1)=A1(6*(i-1)+comp,3*(ia-1)+1)+...
                                    2*mu_nm*sqrt((2*na+3)*(2*na+2)/12/(2*na-1)/(2*na+1))*(na-1)/r*Cp;
            A1(6*(i-1)+comp,3*(ia-1)+2)=A1(6*(i-1)+comp,3*(ia-1)+2)+...
                                    0;
            A1(6*(i-1)+comp,3*(ia-1)+3)=A1(6*(i-1)+comp,3*(ia-1)+3)+...
                                    2*mu_nm*sqrt(na*(2*na-1)*(na+1)/3/(2*na+3)/(2*na+2)/(2*na+1))*(na+2)/r*Cp;
            %\dot u
            A2(6*(i-1)+comp,3*(ia-1)+1)=A2(6*(i-1)+comp,3*(ia-1)+1)+...
                                    -2*mu_nm*sqrt((2*na+3)*(2*na+2)/12/(2*na-1)/(2*na+1))*Cp;
            A2(6*(i-1)+comp,3*(ia-1)+2)=A2(6*(i-1)+comp,3*(ia-1)+2)+...
                                    0;
            A2(6*(i-1)+comp,3*(ia-1)+3)=A2(6*(i-1)+comp,3*(ia-1)+3)+...
                                    2*mu_nm*sqrt(na*(2*na-1)*(na+1)/3/(2*na+3)/(2*na+2)/(2*na+1))*Cp;
            %n-2--------
            Cp=Coup(ieq,ia,8,ireo);
            %u
            A1(6*(i-1)+comp,3*(ia-1)+1)=A1(6*(i-1)+comp,3*(ia-1)+1)+...
                                    2*mu_nm*sqrt((na-1)/(2*na-1))*na/r*Cp;
            A1(6*(i-1)+comp,3*(ia-1)+2)=A1(6*(i-1)+comp,3*(ia-1)+2)+...
                                    0;
            A1(6*(i-1)+comp,3*(ia-1)+3)=A1(6*(i-1)+comp,3*(ia-1)+3)+...
                                    0;
            %\dot u
            A2(6*(i-1)+comp,3*(ia-1)+1)=A2(6*(i-1)+comp,3*(ia-1)+1)+...
                                    2*mu_nm*sqrt((na-1)/(2*na-1))*Cp;
            A2(6*(i-1)+comp,3*(ia-1)+2)=A2(6*(i-1)+comp,3*(ia-1)+2)+...
                                    0;
            A2(6*(i-1)+comp,3*(ia-1)+3)=A2(6*(i-1)+comp,3*(ia-1)+3)+...
                                    0;
            %n+2--------
            Cp=Coup(ieq,ia,9,ireo); 
            %u
            A1(6*(i-1)+comp,3*(ia-1)+1)=A1(6*(i-1)+comp,3*(ia-1)+1)+...
                                    0;
            A1(6*(i-1)+comp,3*(ia-1)+2)=A1(6*(i-1)+comp,3*(ia-1)+2)+...
                                    0;
            A1(6*(i-1)+comp,3*(ia-1)+3)=A1(6*(i-1)+comp,3*(ia-1)+3)+...
                                    2*mu_nm*sqrt((na+2)/(2*na+3))*(na+1)/r*Cp;
            %\dot u
            A2(6*(i-1)+comp,3*(ia-1)+1)=A2(6*(i-1)+comp,3*(ia-1)+1)+...
                                    0;
            A2(6*(i-1)+comp,3*(ia-1)+2)=A2(6*(i-1)+comp,3*(ia-1)+2)+...
                                    0;
            A2(6*(i-1)+comp,3*(ia-1)+3)=A2(6*(i-1)+comp,3*(ia-1)+3)+...
                                    -2*mu_nm*sqrt((na+2)/(2*na+3))*Cp;
            %n-1,0--------
            Cp=Coup(ieq,ia,10,ireo); 
            %u
            A1(6*(i-1)+comp,3*(ia-1)+1)=A1(6*(i-1)+comp,3*(ia-1)+1)+...
                                    0;
            A1(6*(i-1)+comp,3*(ia-1)+2)=A1(6*(i-1)+comp,3*(ia-1)+2)+...
                                    2*mu_nm/sqrt(2)*sqrt((na-1)/(2*na+1))*(na+1)/r*Cp;
            A1(6*(i-1)+comp,3*(ia-1)+3)=A1(6*(i-1)+comp,3*(ia-1)+3)+...
                                    0;
            %\dot u
            A2(6*(i-1)+comp,3*(ia-1)+1)=A2(6*(i-1)+comp,3*(ia-1)+1)+...
                                    0;
            A2(6*(i-1)+comp,3*(ia-1)+2)=A2(6*(i-1)+comp,3*(ia-1)+2)+...
                                    2*mu_nm/sqrt(2)*sqrt((na-1)/(2*na+1))*Cp;
            A2(6*(i-1)+comp,3*(ia-1)+3)=A2(6*(i-1)+comp,3*(ia-1)+3)+...
                                    0;
            %n+1,0--------
            Cp=Coup(ieq,ia,11,ireo);
            %u
            A1(6*(i-1)+comp,3*(ia-1)+1)=A1(6*(i-1)+comp,3*(ia-1)+1)+...
                                    0;
            A1(6*(i-1)+comp,3*(ia-1)+2)=A1(6*(i-1)+comp,3*(ia-1)+2)+...
                                    2*mu_nm/sqrt(2)*sqrt((na+2)/(2*na+1))*na/r*Cp;
            A1(6*(i-1)+comp,3*(ia-1)+3)=A1(6*(i-1)+comp,3*(ia-1)+3)+...
                                    0;
            %\dot u
            A2(6*(i-1)+comp,3*(ia-1)+1)=A2(6*(i-1)+comp,3*(ia-1)+1)+...
                                    0;
            A2(6*(i-1)+comp,3*(ia-1)+2)=A2(6*(i-1)+comp,3*(ia-1)+2)+...
                                    -2*mu_nm/sqrt(2)*sqrt((na+2)/(2*na+1))*Cp;
            A2(6*(i-1)+comp,3*(ia-1)+3)=A2(6*(i-1)+comp,3*(ia-1)+3)+...
                                    0;
            %\sigma_{n,n-1,2}-----------------------------
            comp=3; %couplings 17-21
            %n--------
            Cp=Coup(ieq,ia,17,ireo); 
            %u
            A1(6*(i-1)+comp,3*(ia-1)+1)=A1(6*(i-1)+comp,3*(ia-1)+1)+...
                                    2*mu_nm*sqrt((2*na+3)*(2*na+2)/12/(2*na-1)/(2*na+1))*(na-1)/r*Cp;
            A1(6*(i-1)+comp,3*(ia-1)+2)=A1(6*(i-1)+comp,3*(ia-1)+2)+...
                                    0;
            A1(6*(i-1)+comp,3*(ia-1)+3)=A1(6*(i-1)+comp,3*(ia-1)+3)+...
                                    2*mu_nm*sqrt(na*(2*na-1)*(na+1)/3/(2*na+3)/(2*na+2)/(2*na+1))*(na+2)/r*Cp;
            %\dot u
            A2(6*(i-1)+comp,3*(ia-1)+1)=A2(6*(i-1)+comp,3*(ia-1)+1)+...
                                    -2*mu_nm*sqrt((2*na+3)*(2*na+2)/12/(2*na-1)/(2*na+1))*Cp;
            A2(6*(i-1)+comp,3*(ia-1)+2)=A2(6*(i-1)+comp,3*(ia-1)+2)+...
                                    0;
            A2(6*(i-1)+comp,3*(ia-1)+3)=A2(6*(i-1)+comp,3*(ia-1)+3)+...
                                    2*mu_nm*sqrt(na*(2*na-1)*(na+1)/3/(2*na+3)/(2*na+2)/(2*na+1))*Cp;
            %n-2--------
            Cp=Coup(ieq,ia,18,ireo);
            %u
            A1(6*(i-1)+comp,3*(ia-1)+1)=A1(6*(i-1)+comp,3*(ia-1)+1)+...
                                    2*mu_nm*sqrt((na-1)/(2*na-1))*na/r*Cp;
            A1(6*(i-1)+comp,3*(ia-1)+2)=A1(6*(i-1)+comp,3*(ia-1)+2)+...
                                    0;
            A1(6*(i-1)+comp,3*(ia-1)+3)=A1(6*(i-1)+comp,3*(ia-1)+3)+...
                                    0;
            %\dot u
            A2(6*(i-1)+comp,3*(ia-1)+1)=A2(6*(i-1)+comp,3*(ia-1)+1)+...
                                    2*mu_nm*sqrt((na-1)/(2*na-1))*Cp;
            A2(6*(i-1)+comp,3*(ia-1)+2)=A2(6*(i-1)+comp,3*(ia-1)+2)+...
                                    0;
            A2(6*(i-1)+comp,3*(ia-1)+3)=A2(6*(i-1)+comp,3*(ia-1)+3)+...
                                    0;
            %n+2--------
            Cp=Coup(ieq,ia,19,ireo); 
            %u
            A1(6*(i-1)+comp,3*(ia-1)+1)=A1(6*(i-1)+comp,3*(ia-1)+1)+...
                                    0;
            A1(6*(i-1)+comp,3*(ia-1)+2)=A1(6*(i-1)+comp,3*(ia-1)+2)+...
                                    0;
            A1(6*(i-1)+comp,3*(ia-1)+3)=A1(6*(i-1)+comp,3*(ia-1)+3)+...
                                    2*mu_nm*sqrt((na+2)/(2*na+3))*(na+1)/r*Cp;
            %\dot u
            A2(6*(i-1)+comp,3*(ia-1)+1)=A2(6*(i-1)+comp,3*(ia-1)+1)+...
                                    0;
            A2(6*(i-1)+comp,3*(ia-1)+2)=A2(6*(i-1)+comp,3*(ia-1)+2)+...
                                    0;
            A2(6*(i-1)+comp,3*(ia-1)+3)=A2(6*(i-1)+comp,3*(ia-1)+3)+...
                                    -2*mu_nm*sqrt((na+2)/(2*na+3))*Cp;
            %n-1,0--------
            Cp=Coup(ieq,ia,20,ireo); 
            %u
            A1(6*(i-1)+comp,3*(ia-1)+1)=A1(6*(i-1)+comp,3*(ia-1)+1)+...
                                    0;
            A1(6*(i-1)+comp,3*(ia-1)+2)=A1(6*(i-1)+comp,3*(ia-1)+2)+...
                                    2*mu_nm/sqrt(2)*sqrt((na-1)/(2*na+1))*(na+1)/r*Cp;
            A1(6*(i-1)+comp,3*(ia-1)+3)=A1(6*(i-1)+comp,3*(ia-1)+3)+...
                                    0;
            %\dot u
            A2(6*(i-1)+comp,3*(ia-1)+1)=A2(6*(i-1)+comp,3*(ia-1)+1)+...
                                    0;
            A2(6*(i-1)+comp,3*(ia-1)+2)=A2(6*(i-1)+comp,3*(ia-1)+2)+...
                                    2*mu_nm/sqrt(2)*sqrt((na-1)/(2*na+1))*Cp;
            A2(6*(i-1)+comp,3*(ia-1)+3)=A2(6*(i-1)+comp,3*(ia-1)+3)+...
                                    0;
            %n+1,0--------
            Cp=Coup(ieq,ia,21,ireo);
            %u
            A1(6*(i-1)+comp,3*(ia-1)+1)=A1(6*(i-1)+comp,3*(ia-1)+1)+...
                                    0;
            A1(6*(i-1)+comp,3*(ia-1)+2)=A1(6*(i-1)+comp,3*(ia-1)+2)+...
                                    2*mu_nm/sqrt(2)*sqrt((na+2)/(2*na+1))*na/r*Cp;
            A1(6*(i-1)+comp,3*(ia-1)+3)=A1(6*(i-1)+comp,3*(ia-1)+3)+...
                                    0;
            %\dot u
            A2(6*(i-1)+comp,3*(ia-1)+1)=A2(6*(i-1)+comp,3*(ia-1)+1)+...
                                    0;
            A2(6*(i-1)+comp,3*(ia-1)+2)=A2(6*(i-1)+comp,3*(ia-1)+2)+...
                                    -2*mu_nm/sqrt(2)*sqrt((na+2)/(2*na+1))*Cp;
            A2(6*(i-1)+comp,3*(ia-1)+3)=A2(6*(i-1)+comp,3*(ia-1)+3)+...
                                    0;

            %\sigma_{n,n,2}-------------------------------
            comp=4; %couplings 2-6
            %n--------
            %n--------
            Cp=Coup(ieq,ia,2,ireo); 
            %u
            A1(6*(i-1)+comp,3*(ia-1)+1)=A1(6*(i-1)+comp,3*(ia-1)+1)+...
                                    2*mu_nm*sqrt((2*na+3)*(2*na+2)/12/(2*na-1)/(2*na+1))*(na-1)/r*Cp;
            A1(6*(i-1)+comp,3*(ia-1)+2)=A1(6*(i-1)+comp,3*(ia-1)+2)+...
                                    0;
            A1(6*(i-1)+comp,3*(ia-1)+3)=A1(6*(i-1)+comp,3*(ia-1)+3)+...
                                    2*mu_nm*sqrt(na*(2*na-1)*(na+1)/3/(2*na+3)/(2*na+2)/(2*na+1))*(na+2)/r*Cp;
            %\dot u
            A2(6*(i-1)+comp,3*(ia-1)+1)=A2(6*(i-1)+comp,3*(ia-1)+1)+...
                                    -2*mu_nm*sqrt((2*na+3)*(2*na+2)/12/(2*na-1)/(2*na+1))*Cp;
            A2(6*(i-1)+comp,3*(ia-1)+2)=A2(6*(i-1)+comp,3*(ia-1)+2)+...
                                    0;
            A2(6*(i-1)+comp,3*(ia-1)+3)=A2(6*(i-1)+comp,3*(ia-1)+3)+...
                                    2*mu_nm*sqrt(na*(2*na-1)*(na+1)/3/(2*na+3)/(2*na+2)/(2*na+1))*Cp;
            %n-2--------
            Cp=Coup(ieq,ia,3,ireo);
            %u
            A1(6*(i-1)+comp,3*(ia-1)+1)=A1(6*(i-1)+comp,3*(ia-1)+1)+...
                                    2*mu_nm*sqrt((na-1)/(2*na-1))*na/r*Cp;
            A1(6*(i-1)+comp,3*(ia-1)+2)=A1(6*(i-1)+comp,3*(ia-1)+2)+...
                                    0;
            A1(6*(i-1)+comp,3*(ia-1)+3)=A1(6*(i-1)+comp,3*(ia-1)+3)+...
                                    0;
            %\dot u
            A2(6*(i-1)+comp,3*(ia-1)+1)=A2(6*(i-1)+comp,3*(ia-1)+1)+...
                                    2*mu_nm*sqrt((na-1)/(2*na-1))*Cp;
            A2(6*(i-1)+comp,3*(ia-1)+2)=A2(6*(i-1)+comp,3*(ia-1)+2)+...
                                    0;
            A2(6*(i-1)+comp,3*(ia-1)+3)=A2(6*(i-1)+comp,3*(ia-1)+3)+...
                                    0;
            %n+2--------
            Cp=Coup(ieq,ia,4,ireo); 
            %u
            A1(6*(i-1)+comp,3*(ia-1)+1)=A1(6*(i-1)+comp,3*(ia-1)+1)+...
                                    0;
            A1(6*(i-1)+comp,3*(ia-1)+2)=A1(6*(i-1)+comp,3*(ia-1)+2)+...
                                    0;
            A1(6*(i-1)+comp,3*(ia-1)+3)=A1(6*(i-1)+comp,3*(ia-1)+3)+...
                                    2*mu_nm*sqrt((na+2)/(2*na+3))*(na+1)/r*Cp;
            %\dot u
            A2(6*(i-1)+comp,3*(ia-1)+1)=A2(6*(i-1)+comp,3*(ia-1)+1)+...
                                    0;
            A2(6*(i-1)+comp,3*(ia-1)+2)=A2(6*(i-1)+comp,3*(ia-1)+2)+...
                                    0;
            A2(6*(i-1)+comp,3*(ia-1)+3)=A2(6*(i-1)+comp,3*(ia-1)+3)+...
                                    -2*mu_nm*sqrt((na+2)/(2*na+3))*Cp;
            %n-1,0--------
            Cp=Coup(ieq,ia,5,ireo); 
            %u
            A1(6*(i-1)+comp,3*(ia-1)+1)=A1(6*(i-1)+comp,3*(ia-1)+1)+...
                                    0;
            A1(6*(i-1)+comp,3*(ia-1)+2)=A1(6*(i-1)+comp,3*(ia-1)+2)+...
                                    2*mu_nm/sqrt(2)*sqrt((na-1)/(2*na+1))*(na+1)/r*Cp;
            A1(6*(i-1)+comp,3*(ia-1)+3)=A1(6*(i-1)+comp,3*(ia-1)+3)+...
                                    0;
            %\dot u
            A2(6*(i-1)+comp,3*(ia-1)+1)=A2(6*(i-1)+comp,3*(ia-1)+1)+...
                                    0;
            A2(6*(i-1)+comp,3*(ia-1)+2)=A2(6*(i-1)+comp,3*(ia-1)+2)+...
                                    2*mu_nm/sqrt(2)*sqrt((na-1)/(2*na+1))*Cp;
            A2(6*(i-1)+comp,3*(ia-1)+3)=A2(6*(i-1)+comp,3*(ia-1)+3)+...
                                    0;
            %n+1,0--------
            Cp=Coup(ieq,ia,6,ireo);
            %u
            A1(6*(i-1)+comp,3*(ia-1)+1)=A1(6*(i-1)+comp,3*(ia-1)+1)+...
                                    0;
            A1(6*(i-1)+comp,3*(ia-1)+2)=A1(6*(i-1)+comp,3*(ia-1)+2)+...
                                    2*mu_nm/sqrt(2)*sqrt((na+2)/(2*na+1))*na/r*Cp;
            A1(6*(i-1)+comp,3*(ia-1)+3)=A1(6*(i-1)+comp,3*(ia-1)+3)+...
                                    0;
            %\dot u
            A2(6*(i-1)+comp,3*(ia-1)+1)=A2(6*(i-1)+comp,3*(ia-1)+1)+...
                                    0;
            A2(6*(i-1)+comp,3*(ia-1)+2)=A2(6*(i-1)+comp,3*(ia-1)+2)+...
                                    -2*mu_nm/sqrt(2)*sqrt((na+2)/(2*na+1))*Cp;
            A2(6*(i-1)+comp,3*(ia-1)+3)=A2(6*(i-1)+comp,3*(ia-1)+3)+...
                                    0;
 

            %\sigma_{n,n+1,2}-----------------------------
            comp=5; %couplings 22-26
            %n--------
            %n--------
            Cp=Coup(ieq,ia,22,ireo); 
            %u
            A1(6*(i-1)+comp,3*(ia-1)+1)=A1(6*(i-1)+comp,3*(ia-1)+1)+...
                                    2*mu_nm*sqrt((2*na+3)*(2*na+2)/12/(2*na-1)/(2*na+1))*(na-1)/r*Cp;
            A1(6*(i-1)+comp,3*(ia-1)+2)=A1(6*(i-1)+comp,3*(ia-1)+2)+...
                                    0;
            A1(6*(i-1)+comp,3*(ia-1)+3)=A1(6*(i-1)+comp,3*(ia-1)+3)+...
                                    2*mu_nm*sqrt(na*(2*na-1)*(na+1)/3/(2*na+3)/(2*na+2)/(2*na+1))*(na+2)/r*Cp;
            %\dot u
            A2(6*(i-1)+comp,3*(ia-1)+1)=A2(6*(i-1)+comp,3*(ia-1)+1)+...
                                    -2*mu_nm*sqrt((2*na+3)*(2*na+2)/12/(2*na-1)/(2*na+1))*Cp;
            A2(6*(i-1)+comp,3*(ia-1)+2)=A2(6*(i-1)+comp,3*(ia-1)+2)+...
                                    0;
            A2(6*(i-1)+comp,3*(ia-1)+3)=A2(6*(i-1)+comp,3*(ia-1)+3)+...
                                    2*mu_nm*sqrt(na*(2*na-1)*(na+1)/3/(2*na+3)/(2*na+2)/(2*na+1))*Cp;
            %n-2--------
            Cp=Coup(ieq,ia,23,ireo);
            %u
            A1(6*(i-1)+comp,3*(ia-1)+1)=A1(6*(i-1)+comp,3*(ia-1)+1)+...
                                    2*mu_nm*sqrt((na-1)/(2*na-1))*na/r*Cp;
            A1(6*(i-1)+comp,3*(ia-1)+2)=A1(6*(i-1)+comp,3*(ia-1)+2)+...
                                    0;
            A1(6*(i-1)+comp,3*(ia-1)+3)=A1(6*(i-1)+comp,3*(ia-1)+3)+...
                                    0;
            %\dot u
            A2(6*(i-1)+comp,3*(ia-1)+1)=A2(6*(i-1)+comp,3*(ia-1)+1)+...
                                    2*mu_nm*sqrt((na-1)/(2*na-1))*Cp;
            A2(6*(i-1)+comp,3*(ia-1)+2)=A2(6*(i-1)+comp,3*(ia-1)+2)+...
                                    0;
            A2(6*(i-1)+comp,3*(ia-1)+3)=A2(6*(i-1)+comp,3*(ia-1)+3)+...
                                    0;
            %n+2--------
            Cp=Coup(ieq,ia,24,ireo); 
            %u
            A1(6*(i-1)+comp,3*(ia-1)+1)=A1(6*(i-1)+comp,3*(ia-1)+1)+...
                                    0;
            A1(6*(i-1)+comp,3*(ia-1)+2)=A1(6*(i-1)+comp,3*(ia-1)+2)+...
                                    0;
            A1(6*(i-1)+comp,3*(ia-1)+3)=A1(6*(i-1)+comp,3*(ia-1)+3)+...
                                    2*mu_nm*sqrt((na+2)/(2*na+3))*(na+1)/r*Cp;
            %\dot u
            A2(6*(i-1)+comp,3*(ia-1)+1)=A2(6*(i-1)+comp,3*(ia-1)+1)+...
                                    0;
            A2(6*(i-1)+comp,3*(ia-1)+2)=A2(6*(i-1)+comp,3*(ia-1)+2)+...
                                    0;
            A2(6*(i-1)+comp,3*(ia-1)+3)=A2(6*(i-1)+comp,3*(ia-1)+3)+...
                                    -2*mu_nm*sqrt((na+2)/(2*na+3))*Cp;
            %n-1,0--------
            Cp=Coup(ieq,ia,25,ireo); 
            %u
            A1(6*(i-1)+comp,3*(ia-1)+1)=A1(6*(i-1)+comp,3*(ia-1)+1)+...
                                    0;
            A1(6*(i-1)+comp,3*(ia-1)+2)=A1(6*(i-1)+comp,3*(ia-1)+2)+...
                                    2*mu_nm/sqrt(2)*sqrt((na-1)/(2*na+1))*(na+1)/r*Cp;
            A1(6*(i-1)+comp,3*(ia-1)+3)=A1(6*(i-1)+comp,3*(ia-1)+3)+...
                                    0;
            %\dot u
            A2(6*(i-1)+comp,3*(ia-1)+1)=A2(6*(i-1)+comp,3*(ia-1)+1)+...
                                    0;
            A2(6*(i-1)+comp,3*(ia-1)+2)=A2(6*(i-1)+comp,3*(ia-1)+2)+...
                                    2*mu_nm/sqrt(2)*sqrt((na-1)/(2*na+1))*Cp;
            A2(6*(i-1)+comp,3*(ia-1)+3)=A2(6*(i-1)+comp,3*(ia-1)+3)+...
                                    0;
            %n+1,0--------
            Cp=Coup(ieq,ia,26,ireo);
            %u
            A1(6*(i-1)+comp,3*(ia-1)+1)=A1(6*(i-1)+comp,3*(ia-1)+1)+...
                                    0;
            A1(6*(i-1)+comp,3*(ia-1)+2)=A1(6*(i-1)+comp,3*(ia-1)+2)+...
                                    2*mu_nm/sqrt(2)*sqrt((na+2)/(2*na+1))*na/r*Cp;
            A1(6*(i-1)+comp,3*(ia-1)+3)=A1(6*(i-1)+comp,3*(ia-1)+3)+...
                                    0;
            %\dot u
            A2(6*(i-1)+comp,3*(ia-1)+1)=A2(6*(i-1)+comp,3*(ia-1)+1)+...
                                    0;
            A2(6*(i-1)+comp,3*(ia-1)+2)=A2(6*(i-1)+comp,3*(ia-1)+2)+...
                                    -2*mu_nm/sqrt(2)*sqrt((na+2)/(2*na+1))*Cp;
            A2(6*(i-1)+comp,3*(ia-1)+3)=A2(6*(i-1)+comp,3*(ia-1)+3)+...
                                    0;


            %\sigma_{n,n+2,2}-----------------------------
            comp=6; %couplings 12-16
            %n--------
            %n--------
            Cp=Coup(ieq,ia,12,ireo); 
            %u
            A1(6*(i-1)+comp,3*(ia-1)+1)=A1(6*(i-1)+comp,3*(ia-1)+1)+...
                                    2*mu_nm*sqrt((2*na+3)*(2*na+2)/12/(2*na-1)/(2*na+1))*(na-1)/r*Cp;
            A1(6*(i-1)+comp,3*(ia-1)+2)=A1(6*(i-1)+comp,3*(ia-1)+2)+...
                                    0;
            A1(6*(i-1)+comp,3*(ia-1)+3)=A1(6*(i-1)+comp,3*(ia-1)+3)+...
                                    2*mu_nm*sqrt(na*(2*na-1)*(na+1)/3/(2*na+3)/(2*na+2)/(2*na+1))*(na+2)/r*Cp;
            %\dot u
            A2(6*(i-1)+comp,3*(ia-1)+1)=A2(6*(i-1)+comp,3*(ia-1)+1)+...
                                    -2*mu_nm*sqrt((2*na+3)*(2*na+2)/12/(2*na-1)/(2*na+1))*Cp;
            A2(6*(i-1)+comp,3*(ia-1)+2)=A2(6*(i-1)+comp,3*(ia-1)+2)+...
                                    0;
            A2(6*(i-1)+comp,3*(ia-1)+3)=A2(6*(i-1)+comp,3*(ia-1)+3)+...
                                    2*mu_nm*sqrt(na*(2*na-1)*(na+1)/3/(2*na+3)/(2*na+2)/(2*na+1))*Cp;
            %n-2--------
            Cp=Coup(ieq,ia,13,ireo);
            %u
            A1(6*(i-1)+comp,3*(ia-1)+1)=A1(6*(i-1)+comp,3*(ia-1)+1)+...
                                    2*mu_nm*sqrt((na-1)/(2*na-1))*na/r*Cp;
            A1(6*(i-1)+comp,3*(ia-1)+2)=A1(6*(i-1)+comp,3*(ia-1)+2)+...
                                    0;
            A1(6*(i-1)+comp,3*(ia-1)+3)=A1(6*(i-1)+comp,3*(ia-1)+3)+...
                                    0;
            %\dot u
            A2(6*(i-1)+comp,3*(ia-1)+1)=A2(6*(i-1)+comp,3*(ia-1)+1)+...
                                    2*mu_nm*sqrt((na-1)/(2*na-1))*Cp;
            A2(6*(i-1)+comp,3*(ia-1)+2)=A2(6*(i-1)+comp,3*(ia-1)+2)+...
                                    0;
            A2(6*(i-1)+comp,3*(ia-1)+3)=A2(6*(i-1)+comp,3*(ia-1)+3)+...
                                    0;
            %n+2--------
            Cp=Coup(ieq,ia,14,ireo); 
            %u
            A1(6*(i-1)+comp,3*(ia-1)+1)=A1(6*(i-1)+comp,3*(ia-1)+1)+...
                                    0;
            A1(6*(i-1)+comp,3*(ia-1)+2)=A1(6*(i-1)+comp,3*(ia-1)+2)+...
                                    0;
            A1(6*(i-1)+comp,3*(ia-1)+3)=A1(6*(i-1)+comp,3*(ia-1)+3)+...
                                    2*mu_nm*sqrt((na+2)/(2*na+3))*(na+1)/r*Cp;
            %\dot u
            A2(6*(i-1)+comp,3*(ia-1)+1)=A2(6*(i-1)+comp,3*(ia-1)+1)+...
                                    0;
            A2(6*(i-1)+comp,3*(ia-1)+2)=A2(6*(i-1)+comp,3*(ia-1)+2)+...
                                    0;
            A2(6*(i-1)+comp,3*(ia-1)+3)=A2(6*(i-1)+comp,3*(ia-1)+3)+...
                                    -2*mu_nm*sqrt((na+2)/(2*na+3))*Cp;
            %n-1,0--------
            Cp=Coup(ieq,ia,15,ireo); 
            %u
            A1(6*(i-1)+comp,3*(ia-1)+1)=A1(6*(i-1)+comp,3*(ia-1)+1)+...
                                    0;
            A1(6*(i-1)+comp,3*(ia-1)+2)=A1(6*(i-1)+comp,3*(ia-1)+2)+...
                                    2*mu_nm/sqrt(2)*sqrt((na-1)/(2*na+1))*(na+1)/r*Cp;
            A1(6*(i-1)+comp,3*(ia-1)+3)=A1(6*(i-1)+comp,3*(ia-1)+3)+...
                                    0;
            %\dot u
            A2(6*(i-1)+comp,3*(ia-1)+1)=A2(6*(i-1)+comp,3*(ia-1)+1)+...
                                    0;
            A2(6*(i-1)+comp,3*(ia-1)+2)=A2(6*(i-1)+comp,3*(ia-1)+2)+...
                                    2*mu_nm/sqrt(2)*sqrt((na-1)/(2*na+1))*Cp;
            A2(6*(i-1)+comp,3*(ia-1)+3)=A2(6*(i-1)+comp,3*(ia-1)+3)+...
                                    0;
            %n+1,0--------
            Cp=Coup(ieq,ia,16,ireo);
            %u
            A1(6*(i-1)+comp,3*(ia-1)+1)=A1(6*(i-1)+comp,3*(ia-1)+1)+...
                                    0;
            A1(6*(i-1)+comp,3*(ia-1)+2)=A1(6*(i-1)+comp,3*(ia-1)+2)+...
                                    2*mu_nm/sqrt(2)*sqrt((na+2)/(2*na+1))*na/r*Cp;
            A1(6*(i-1)+comp,3*(ia-1)+3)=A1(6*(i-1)+comp,3*(ia-1)+3)+...
                                    0;
            %\dot u
            A2(6*(i-1)+comp,3*(ia-1)+1)=A2(6*(i-1)+comp,3*(ia-1)+1)+...
                                    0;
            A2(6*(i-1)+comp,3*(ia-1)+2)=A2(6*(i-1)+comp,3*(ia-1)+2)+...
                                    -2*mu_nm/sqrt(2)*sqrt((na+2)/(2*na+1))*Cp;
            A2(6*(i-1)+comp,3*(ia-1)+3)=A2(6*(i-1)+comp,3*(ia-1)+3)+...
                                    0;
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
            A1(6*(i-1)+comp,3*(ia-1)+1)=A1(6*(i-1)+comp,3*(ia-1)+1)+...
                                    0;
            A1(6*(i-1)+comp,3*(ia-1)+2)=A1(6*(i-1)+comp,3*(ia-1)+2)+...
                                    0;
            A1(6*(i-1)+comp,3*(ia-1)+3)=A1(6*(i-1)+comp,3*(ia-1)+3)+...
                                    2*mu_nm*sqrt((na+2)/(2*na+3))*(na+1)/r*Cp;
            %\dot u
            A2(6*(i-1)+comp,3*(ia-1)+1)=A2(6*(i-1)+comp,3*(ia-1)+1)+...
                                    0;
            A2(6*(i-1)+comp,3*(ia-1)+2)=A2(6*(i-1)+comp,3*(ia-1)+2)+...
                                    0;
            A2(6*(i-1)+comp,3*(ia-1)+3)=A2(6*(i-1)+comp,3*(ia-1)+3)+...
                                    -2*mu_nm*sqrt((na+2)/(2*na+3))*Cp;           
            %\sigma_{n,n,2}-------------------------------
            comp=4; %couplings 2-6
            %n+2--------
            Cp=Coup(ieq,ia,4,ireo); 
            %u
            A1(6*(i-1)+comp,3*(ia-1)+1)=A1(6*(i-1)+comp,3*(ia-1)+1)+...
                                    0;
            A1(6*(i-1)+comp,3*(ia-1)+2)=A1(6*(i-1)+comp,3*(ia-1)+2)+...
                                    0;
            A1(6*(i-1)+comp,3*(ia-1)+3)=A1(6*(i-1)+comp,3*(ia-1)+3)+...
                                    2*mu_nm*sqrt((na+2)/(2*na+3))*(na+1)/r*Cp;
            %\dot u
            A2(6*(i-1)+comp,3*(ia-1)+1)=A2(6*(i-1)+comp,3*(ia-1)+1)+...
                                    0;
            A2(6*(i-1)+comp,3*(ia-1)+2)=A2(6*(i-1)+comp,3*(ia-1)+2)+...
                                    0;
            A2(6*(i-1)+comp,3*(ia-1)+3)=A2(6*(i-1)+comp,3*(ia-1)+3)+...
                                    -2*mu_nm*sqrt((na+2)/(2*na+3))*Cp;            
            %\sigma_{n,n+1,2}-----------------------------
            comp=5; %couplings 22-26            
            %n+2--------
            Cp=Coup(ieq,ia,24,ireo); 
            %u
            A1(6*(i-1)+comp,3*(ia-1)+1)=A1(6*(i-1)+comp,3*(ia-1)+1)+...
                                    0;
            A1(6*(i-1)+comp,3*(ia-1)+2)=A1(6*(i-1)+comp,3*(ia-1)+2)+...
                                    0;
            A1(6*(i-1)+comp,3*(ia-1)+3)=A1(6*(i-1)+comp,3*(ia-1)+3)+...
                                    2*mu_nm*sqrt((na+2)/(2*na+3))*(na+1)/r*Cp;
            %\dot u
            A2(6*(i-1)+comp,3*(ia-1)+1)=A2(6*(i-1)+comp,3*(ia-1)+1)+...
                                    0;
            A2(6*(i-1)+comp,3*(ia-1)+2)=A2(6*(i-1)+comp,3*(ia-1)+2)+...
                                    0;
            A2(6*(i-1)+comp,3*(ia-1)+3)=A2(6*(i-1)+comp,3*(ia-1)+3)+...
                                    -2*mu_nm*sqrt((na+2)/(2*na+3))*Cp;            
            %\sigma_{n,n+2,2}-----------------------------
            comp=6; %couplings 12-16
            %n+2--------
            Cp=Coup(ieq,ia,14,ireo); 
            %u
            A1(6*(i-1)+comp,3*(ia-1)+1)=A1(6*(i-1)+comp,3*(ia-1)+1)+...
                                    0;
            A1(6*(i-1)+comp,3*(ia-1)+2)=A1(6*(i-1)+comp,3*(ia-1)+2)+...
                                    0;
            A1(6*(i-1)+comp,3*(ia-1)+3)=A1(6*(i-1)+comp,3*(ia-1)+3)+...
                                    2*mu_nm*sqrt((na+2)/(2*na+3))*(na+1)/r*Cp;
            %\dot u
            A2(6*(i-1)+comp,3*(ia-1)+1)=A2(6*(i-1)+comp,3*(ia-1)+1)+...
                                    0;
            A2(6*(i-1)+comp,3*(ia-1)+2)=A2(6*(i-1)+comp,3*(ia-1)+2)+...
                                    0;
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
    A14(6*(i-1)+1,3*(i-1)+2)=0;
    A14(6*(i-1)+1,3*(i-1)+3)=1/sqrt(3)*1/sqrt(2*n+1)*sqrt(n+1);  
    %u
    A15(6*(i-1)+1,3*(i-1)+1)=1/sqrt(3)*1/sqrt(2*n+1)*(n-1)*sqrt(n);
    A15(6*(i-1)+1,3*(i-1)+2)=0;
    A15(6*(i-1)+1,3*(i-1)+3)=1/sqrt(3)*1/sqrt(2*n+1)*sqrt(n+1)*(n+2);
    %\sigma_{n,n-2,2}------------------------------
    %u dot
    A14(6*(i-1)+2,3*(i-1)+1)=sqrt((n-1)/(2*n-1));
    A14(6*(i-1)+2,3*(i-1)+2)=0;
    A14(6*(i-1)+2,3*(i-1)+3)=0;  
    %u
    A15(6*(i-1)+2,3*(i-1)+1)=sqrt((n-1)/(2*n-1))*n;
    A15(6*(i-1)+2,3*(i-1)+2)=0;
    A15(6*(i-1)+2,3*(i-1)+3)=0;
    %\sigma_{n,n-1,2}------------------------------
    %u dot
    A14(6*(i-1)+3,3*(i-1)+1)=0;
    A14(6*(i-1)+3,3*(i-1)+2)=1/sqrt(2)*sqrt((n-1)/(2*n+1));
    A14(6*(i-1)+3,3*(i-1)+3)=0;  
    %u
    A15(6*(i-1)+3,3*(i-1)+1)=0;
    A15(6*(i-1)+3,3*(i-1)+2)=1/sqrt(2)*sqrt((n-1)/(2*n+1))*(n+1);
    A15(6*(i-1)+3,3*(i-1)+3)=0;
    %\sigma_{n,n,2}------------------------------
    %u dot
    A14(6*(i-1)+4,3*(i-1)+1)=-sqrt((2*n+3)*(2*n+2)/12/(2*n-1)/(2*n+1));
    A14(6*(i-1)+4,3*(i-1)+2)=0;
    A14(6*(i-1)+4,3*(i-1)+3)=sqrt(n*(2*n-1)*(n+1)/3/(2*n+3)/(2*n+2)/(2*n+1));  
    %u
    A15(6*(i-1)+4,3*(i-1)+1)=sqrt((2*n+3)*(2*n+2)/12/(2*n-1)/(2*n+1))*(n-1);
    A15(6*(i-1)+4,3*(i-1)+2)=0;
    A15(6*(i-1)+4,3*(i-1)+3)=sqrt(n*(2*n-1)*(n+1)/3/(2*n+3)/(2*n+2)/(2*n+1))*(n+2);
    %\sigma_{n,n+1,2}------------------------------
    %u dot
    A14(6*(i-1)+5,3*(i-1)+1)=0;
    A14(6*(i-1)+5,3*(i-1)+2)=-1/sqrt(2)*sqrt((n+2)/(2*n+1));
    A14(6*(i-1)+5,3*(i-1)+3)=0;  
    %u
    A15(6*(i-1)+5,3*(i-1)+1)=0;
    A15(6*(i-1)+5,3*(i-1)+2)=1/sqrt(2)*sqrt((n+2)/(2*n+1))*n;
    A15(6*(i-1)+5,3*(i-1)+3)=0;
    %\sigma_{n,n+2,2}------------------------------
    %u dot
    A14(6*(i-1)+6,3*(i-1)+1)=0;
    A14(6*(i-1)+6,3*(i-1)+2)=0;
    A14(6*(i-1)+6,3*(i-1)+3)=-sqrt((n+2)/(2*n+3));  
    %u
    A15(6*(i-1)+6,3*(i-1)+1)=0;
    A15(6*(i-1)+6,3*(i-1)+2)=0;
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
    A3(3*(i-1)+1,3*(i-1)+2)=0; %u_{n}
    A3(3*(i-1)+1,3*(i-1)+3)=-sqrt(n+1)/sqrt(2*n+1); %u_{n+1}
    %V_n^m
    A3(3*(i-1)+2,3*(i-1)+1)=1/sqrt(2*n+1)/sqrt(n); %u_{n-1}
    A3(3*(i-1)+2,3*(i-1)+2)=0; %u_{n}
    A3(3*(i-1)+2,3*(i-1)+3)=1/sqrt(2*n+1)/sqrt(n+1); %u_{n+1}
    %W_n^m
    A3(3*(i-1)+3,3*(i-1)+1)=0; %u_{n-1}
    A3(3*(i-1)+3,3*(i-1)+2)=1i/sqrt(n*(n+1)); %u_{n}
    A3(3*(i-1)+3,3*(i-1)+3)=0; %u_{n+1}  
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
    A4(3*(i-1)+1,6*(i-1)+3)=0; %\sigma_{n,n-1,2}
    A4(3*(i-1)+1,6*(i-1)+4)=-sqrt(n)/(2*n+1)*(sqrt((2*n+3)*(2*n+2)/12/(2*n-1))+sqrt((2*n-1)*(n+1)^2/3/(2*n+2)/(2*n+3))); %\sigma_{n,n,2}
    A4(3*(i-1)+1,6*(i-1)+5)=0; %\sigma_{n,n+1,2}
    A4(3*(i-1)+1,6*(i-1)+6)=sqrt((n+1)*(n+2)/(2*n+1)/(2*n+3)); %\sigma_{n,n+2,2}   
    %S_n^m
    A4(3*(i-1)+2,6*(i-1)+1)=0; %\sigma_{n,n,0}
    A4(3*(i-1)+2,6*(i-1)+2)=sqrt((n-1)/(2*n-1)/n/(2*n+1)); %\sigma_{n,n-2,2}
    A4(3*(i-1)+2,6*(i-1)+3)=0; %\sigma_{n,n-1,2}
    A4(3*(i-1)+2,6*(i-1)+4)=(-sqrt((2*n+3)*(2*n+2)/(12*n)/(2*n-1)/(2*n+1)^2)+sqrt(n*(2*n-1)/3/(2*n+1)^2/(2*n+2)/(2*n+3))); %\sigma_{n,n,2}
    A4(3*(i-1)+2,6*(i-1)+5)=0; %\sigma_{n,n+1,2}
    A4(3*(i-1)+2,6*(i-1)+6)=-sqrt((n+2)/(2*n+3)/(n+1)/(2*n+1)); %\sigma_{n,n+2,2} 
    %T_n^m
    A4(3*(i-1)+3,6*(i-1)+1)=0; %\sigma_{n,n,0}
    A4(3*(i-1)+3,6*(i-1)+2)=0; %\sigma_{n,n-2,2}
    A4(3*(i-1)+3,6*(i-1)+3)=1i/sqrt(2*n*(n+1)*(2*n+1))*sqrt(n-1); %\sigma_{n,n-1,2}
    A4(3*(i-1)+3,6*(i-1)+4)=0; %\sigma_{n,n,2}
    A4(3*(i-1)+3,6*(i-1)+5)=-1i/sqrt(2*n*(n+1)*(2*n+1))*sqrt(n+2); %\sigma_{n,n+1,2}
    A4(3*(i-1)+3,6*(i-1)+6)=0; %\sigma_{n,n+2,2} 
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
    A5(3*(i-1)+1,6*(i-1)+1)=0; %\sigma_{n,n,0}
    A5(3*(i-1)+1,6*(i-1)+2)=-(n-2)*sqrt(n*(n-1)/(2*n-1)/(2*n+1)); %\sigma_{n,n-2,2}
    A5(3*(i-1)+1,6*(i-1)+3)=0; %\sigma_{n,n-1,2}
    A5(3*(i-1)+1,6*(i-1)+4)=1/(2*n+1)*(...
        -(n+1)*sqrt((2*n+3)*(2*n+2)*n/12/(2*n-1))+...
        +n*sqrt(n*(n+1)^2*(2*n-1)/3/(2*n+2)/(2*n+3))); %\sigma_{n,n,2}
    A5(3*(i-1)+1,6*(i-1)+5)=0; %\sigma_{n,n+1,2}
    A5(3*(i-1)+1,6*(i-1)+6)=(n+3)*sqrt((n+1)*(n+2)/(2*n+3)/(2*n+1)); %\sigma_{n,n+2,2}   
    %S_n^m
    A5(3*(i-1)+2,6*(i-1)+1)=-1/sqrt(3); %\sigma_{n,n,0}
    A5(3*(i-1)+2,6*(i-1)+2)=-(n-2)*sqrt((n-1)/(2*n-1)/(2*n+1)/n); %\sigma_{n,n-2,2}
    A5(3*(i-1)+2,6*(i-1)+3)=0; %\sigma_{n,n-1,2}
    A5(3*(i-1)+2,6*(i-1)+4)=-1/(2*n+1)*(...
        (n+1)*sqrt((2*n+3)*(2*n+2)/(2*n-1)/12/n)+...
        n*sqrt(n*(2*n-1)/3/(2*n+2)/(2*n+3))); %\sigma_{n,n,2}
    A5(3*(i-1)+2,6*(i-1)+5)=0; %\sigma_{n,n+1,2}
    A5(3*(i-1)+2,6*(i-1)+6)=-(n+3)*sqrt((n+2)/(2*n+3)/(2*n+1)/(n+1)); %\sigma_{n,n+2,2} 
    %T_n^m
    A5(3*(i-1)+3,6*(i-1)+1)=0; %\sigma_{n,n,0}
    A5(3*(i-1)+3,6*(i-1)+2)=0; %\sigma_{n,n-2,2}
    A5(3*(i-1)+3,6*(i-1)+3)=-(n-1)*sqrt(n-1)/sqrt(2*n*(n+1)*(2*n+1))*1i; %\sigma_{n,n-1,2}
    A5(3*(i-1)+3,6*(i-1)+4)=0; %\sigma_{n,n,2}
    A5(3*(i-1)+3,6*(i-1)+5)=-(n+2)*sqrt(n+2)/sqrt(2*n*(n+1)*(2*n+1))*1i; %\sigma_{n,n+1,2}
    A5(3*(i-1)+3,6*(i-1)+6)=0; %\sigma_{n,n+2,2}
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
function Aprop=get_Aprop(rK,gK,dgK,Nmodes,A1,A2,A3_inv,A4,A5,A6,A71,A72,A81,A82,A9,A100,A101,A102,A11,A12,A13,deg_0,Gg,ocean_flag)
if ocean_flag==0
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
else
    % Laplace equation
    Aprop=zeros(8*Nmodes,8*Nmodes); 
    Aprop(6*Nmodes+1:end,6*Nmodes+1:end)=A100+A101/rK+A102/rK^2; 
end
end