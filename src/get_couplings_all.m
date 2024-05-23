%% Functions used to compute all coupling coefficintes
%% INPUT
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
    % Couplings: Structure containing all the Couplings, see get_couplings
    % for more details about the structure. 
    
%% FUNCTION
function [Couplings] = get_couplings_all(Forcing,Numerics,varargin)
verbose=0;
for k = 1:length(varargin)
    if strcmpi(varargin{k},'verbose')
        verbose=1; 
        varargin{k}=[];
    end
end
%% initialize rheology structure
Nrheo_max = Numerics.Nrheo_max;
mu_v = zeros((Nrheo_max+1)*(Nrheo_max+1)-1,2);
jj=1;
for l=1:Nrheo_max
    for m=0:l
        if m==0
            mu_v(jj,1) = l; 
            mu_v(jj,2) = m;
            jj=jj+1; 
        else
            mu_v(jj,1) = l; 
            mu_v(jj,2) = m; 
            mu_v(jj+1,1) = l;
            mu_v(jj+1,2) = -m;
            jj=jj+2; 
        end
    end
end

%% compute couplings 
if verbose==1
    Couplings = get_couplings(mu_v,Forcing,Numerics,'verbose');
else
    Couplings = get_couplings(mu_v,Forcing,Numerics);
end
Couplings.n_r = mu_v(:,1);%single(mu_v(:,1));
Couplings.m_r = mu_v(:,2);%single(mu_v(:,2));
end