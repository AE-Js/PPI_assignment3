%% get_couplings
%% INPUT
% Variations (array): Contains all the unique rheology modes
% Numerics (struct):
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
    %Numerics.rheology_cutoff: determines which terms of the rheology are included (only relevant for viscoelastic).
    %                          Terms with log10(mu_n^m)-log10(mu_n^m(leading))>=-Numerics.rheology_cutoff are included. 
    %                          Default 0 (only leading terms)
    %Numerics.Nenergy: maximum degree to which energy dissipation is expanded. Defaulft 8. 
    %Numerics.load_couplings: 
         % (1) load coupling coefficintes from a file that contains specifically the coupling coefficients for the rheology specied, 
         %     if it does not exisit, compute it. 
         % (2) load coupling coefficients from file that contains ALL couplings up to rheology variations of a certain degree, 
         %     if such a file does not exist compute it (default) 
    % Numerics.parallel_sol:  Use a parfor-loop to call get_Love, either 0 or 1
    % Numerics.parallel_gen:  Calculate potential coupling files and the propagation inside get_solution using parfor-loops, either 0 or 1
    % Numerics.Nrheo_max:     Highest degree of rheology across the layers that is used 
% ------------------
% Forcing (struct): Contains the forcing information
    % Forcing.Td: forcing period
    % Forcing.n: degree of the forcing 
    % Forcing.m: order of the forcing 
    % Forcing.F: amplitude of the component 
%% OUTPUT
% Couplings 
    % Couplings.n_s: degree of modes that intervene in the solution
    % Couplings.m_s: order of modes that intervene in the solution
    % Couplings.order: order of the coupling 
    % Couplings.Coup: matrix containing coupling coefficients 
        % Coup(Nsol,Nsol,27,Nreo)
            % Coup(i,j,1,k)  (T^ma_{na,na,0}Y_nb^mb,T^mc_{nc,nc,0})        
            % Coup(i,j,2,k)  (T^ma_{na,na,2}Y_nb^mb,T^mc_{nc,nc,2})
            % Coup(i,j,3,k)  (T^ma_{na,na-2,2}Y_nb^mb,T^mc_{nc,nc,2})
            % Coup(i,j,4,k)  (T^ma_{na,na+2,2}Y_nb^mb,T^mc_{nc,nc,2})
            % Coup(i,j,5,k)  (T^ma_{na,na-1,2}Y_nb^mb,T^mc_{nc,nc,2})
            % Coup(i,j,6,k)  (T^ma_{na,na+1,2}Y_nb^mb,T^mc_{nc,nc,2})        
            % Coup(i,j,7,k)   (T^ma_{na,na,2}Y_nb^mb,T^mc_{nc,nc-2,2})
            % Coup(i,j,8,k)   (T^ma_{na,na-2,2}Y_nb^mb,T^mc_{nc,nc-2,2})
            % Coup(i,j,9,k)   (T^ma_{na,na+2,2}Y_nb^mb,T^mc_{nc,nc-2,2})
            % Coup(i,j,10,k)  (T^ma_{na,na-1,2}Y_nb^mb,T^mc_{nc,nc-2,2})
            % Coup(i,j,11,k)  (T^ma_{na,na+1,2}Y_nb^mb,T^mc_{nc,nc-2,2})        
            % Coup(i,j,12,k)  (T^ma_{na,na,2}Y_nb^mb,T^mc_{nc,nc+2,2})
            % Coup(i,j,13,k)  (T^ma_{na,na-2,2}Y_nb^mb,T^mc_{nc,nc+2,2})
            % Coup(i,j,14,k)  (T^ma_{na,na+2,2}Y_nb^mb,T^mc_{nc,nc+2,2})
            % Coup(i,j,15,k)  (T^ma_{na,na-1,2}Y_nb^mb,T^mc_{nc,nc+2,2})
            % Coup(i,j,16,k)  (T^ma_{na,na+1,2}Y_nb^mb,T^mc_{nc,nc+2,2})        
            % Coup(i,j,17,k)  (T^ma_{na,na,2}Y_nb^mb,T^mc_{nc,nc-1,2})
            % Coup(i,j,18,k)  (T^ma_{na,na-2,2}Y_nb^mb,T^mc_{nc,nc-1,2})
            % Coup(i,j,19,k)  (T^ma_{na,na+2,2}Y_nb^mb,T^mc_{nc,nc-1,2})
            % Coup(i,j,20,k)  (T^ma_{na,na-1,2}Y_nb^mb,T^mc_{nc,nc-1,2})
            % Coup(i,j,21,k)  (T^ma_{na,na+1,2}Y_nb^mb,T^mc_{nc,nc-1,2})        
            % Coup(i,j,22,k)  (T^ma_{na,na,2}Y_nb^mb,T^mc_{nc,nc+1,2})
            % Coup(i,j,23,k)  (T^ma_{na,na-2,2}Y_nb^mb,T^mc_{nc,nc+1,2})
            % Coup(i,j,24,k)  (T^ma_{na,na+2,2}Y_nb^mb,T^mc_{nc,nc+1,2})
            % Coup(i,j,25,k)  (T^ma_{na,na-1,2}Y_nb^mb,T^mc_{nc,nc+1,2})
            % Coup(i,j,26,k)  (T^ma_{na,na+1,2}Y_nb^mb,T^mc_{nc,nc+1,2})
            % Coup(i,j,27,k) 0 if all coupling coefficents are 0

%% FUNCTION STARTS HERE 
function [Couplings] = get_couplings(variations,Forcing,Numerics,varargin) 
%% OPTIONAL INPUTS
verbose=0;
for k = 1:length(varargin)
    if strcmpi(varargin{k},'verbose')
        verbose=1;
        varargin{k}=[];
    end
end
%% GET MODES THAT ARE COUPLED USING SELECTION RULES
% Retrieve perturbation order from Numerics
max_order = Numerics.perturbation_order;

% Number of rheology modes
Nreo = size(variations,1);

% Retrieve all activated modes, based on rheology, forcing, and perturbation order
modes4 = get_active_modes(max_order,variations,Forcing);

% Set first entries for the Couplings struct
Couplings.n_s = modes4(:,1);%single(modes4(:,1));
Couplings.m_s = modes4(:,2);%single(modes4(:,2));
Couplings.order = modes4(:,3);%single(modes4(:,3));
Nsol = length(Couplings.n_s);

% GET COUPLING COEFFICIENTS 
if verbose==1
     disp("Getting in the loop to compute coupling coefficients...")
     disp([ num2str(Nsol) ' Modes'])
     disp([ num2str(Nreo) ' Rheology'])
     disp([ num2str((Nreo*Nsol^2)) ' Couplings to be computed'])
     tic
end
Cp=zeros(Nsol,Nsol,27,Nreo);
nb_arr = variations(:,1);
mb_arr = variations(:,2);
na_arr = Couplings.n_s;
ma_arr = Couplings.m_s;

% Check whether it is possible to use parfor
if Numerics.parallel_gen == 1 && Numerics.parallel_sol == 0
    parfor ireo=1:Nreo 
        nb=nb_arr(ireo);
        mb=mb_arr(ireo);
        for imod1=1:Nsol
            n = na_arr(imod1);
            m = ma_arr(imod1);
            for imod2=1:Nsol
                na = na_arr(imod2);
                ma = ma_arr(imod2);
    
                % check if the coefficient is non-zero
                % abs(n-nb):n+nb 
                if na>=abs(n-nb) && na<=abs(n+nb) && ma==m-mb 
                    Cp(imod1,imod2,:,ireo) = coupling_coefficients(n,m,na,ma,nb,mb);              
                end
    
                if verbose==1
                    if mod(imod2,25) == 0
                        disp([ num2str(Nsol) ' Modes'])
                        disp([ num2str(Nreo) ' Rheology'])
                        disp([ num2str((Nreo*Nsol^2)) ' Couplings to be computed'])
                        disp(['Rheology mode: ' num2str(ireo) ', Progress: ' num2str(((imod1-1)*Nsol + imod2)/(Nsol^2)*100) ' % Couplings computed'])
                    end
                end
            end
        end
    end
else
    jj = 1;
    for ireo=1:Nreo
        nb = nb_arr(ireo);
        mb = mb_arr(ireo);
        for imod1=1:Nsol
            n = na_arr(imod1);
            m = ma_arr(imod1);
            for imod2=1:Nsol
                na = na_arr(imod2);
                ma = ma_arr(imod2);
    
                % check if the coefficient is non-zero
                % abs(n-nb):n+nb 
                if na>=abs(n-nb) && na<=abs(n+nb) && ma==m-mb 
                    Cp(imod1,imod2,:,ireo) = coupling_coefficients(n,m,na,ma,nb,mb);              
                end
                jj = jj+1;
                if verbose==1
                    if mod(imod2,25) == 0
                        disp([ num2str(Nsol) ' Modes'])
                        disp([ num2str(Nreo) ' Rheology'])
                        disp([ num2str((Nreo*Nsol^2)) ' Couplings to be computed'])
                        disp(['Progress: ' num2str((jj)/(Nreo*Nsol^2)*100) ' % Couplings computed'])
                    end
                end
            end
        end
    end
end

Couplings.Coup = Cp;%single(Cp);

end

%% FUNCTION TO GET THE COUPLINGS 
function Cp = coupling_coefficients(n,m,na,ma,nb,mb)
% Define Cp array
Cp = zeros(1,1,27);

% Create the nc1 and na1 arrays
nc1v=[n-1 n n+1];
nc1 = repmat(nc1v,1,78)';
na1v=[na-1 na na+1];
na1 = repmat(repelem(na1v,3),1,26)';

% Generate the lc and la arrays
lc=repelem([0; zeros(25,1)+2],9);
la=repelem([0; zeros(25,1)+2],9);

% Generate the nc2 and na2 arrays
nc2_base=[n, n-2, n+2, n-1, n+1];
nc2 = repelem([n repelem(nc2_base,5)]',9);
na2_base=[na, na-2, na+2, na-1, na+1];
na2 = repelem([na repmat(na2_base,1,5)]',9);

% Generate the fixed arrays
nc_vec = zeros(234,1) + n;
mc_vec = zeros(234,1) + m;
na_vec = zeros(234,1) + na;
ma_vec = zeros(234,1) + ma;
nb_vec = zeros(234,1) + nb;
mb_vec = zeros(234,1) + mb; 

% Calculate coupling coefficients
Cp(1,1,1:end-1) = couplingT(nc_vec,nc2,lc,mc_vec,na_vec,na2,la,ma_vec,nb_vec,mb_vec,nc1,na1);
                                      
if max(abs(squeeze(Cp(:))))>0
    Cp(1,1,27)=1;
else
    Cp(1,1,27)=0;
end
end
%%
function C = couplingT(nc,nc2,lc,mc,na,na2,la,ma,nb,mb,nc1,na1)
one_arr = zeros(234,1)+1;

Lama = sqrt((2*la+1).*(2*na1+1));
Lamc = sqrt((2*lc+1).*(2*nc1+1));
Caux = (-1).^(na+na2+la+nc+nc2+lc).*Lama.*Lamc.* ...
       Wigner6j(one_arr, la, one_arr, na, na1, na2).* ...
       Wigner6j(one_arr, lc, one_arr, nc, nc1, nc2).* ...
       couplingY(nc,mc,nc1,nc2,na,ma,na1,na2,nb,mb);

C_int = reshape(Caux,9,[]);
C = sum(C_int);
end
%%
function [C] = couplingY(nc,mc,nc1,nc2,na,ma,na1,na2,nb,mb)
one_arr = zeros(234,1)+1;
zero_arr = zeros(234,1);

Lam = sqrt((2*nc2+1).*(2*nc1+1).*(2*nc+1).*...
        (2*na2+1).*(2*na1+1).*(2*na+1).*...
        (2*nb+1));   
C = (-1).^(na+na1+nc+nc1+mc).*Lam.*...
       Wigner6j(na, na1, one_arr, nc1, nc, nb).*...
       Wigner6j(na1, na2, one_arr, nc2, nc1, nb).*...
       Wigner3j(na2, nc2, nb, zero_arr, zero_arr, zero_arr).*...
       Wigner3j(na,nc,nb,ma,-mc,mb);
end
