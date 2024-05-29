%% set_boundary_indices
%% SHORT DESCRIPTION
% Small helper function that selects the radial boundary locations based 
% on the selected method. It sets the number of points in each layer as well
% as the indices at which the boundary's are located. 
%% INPUT
% Numerics: struct that contains the necessary data and that will be altered
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
    %Numerics.Nr: total number of radial points (computed inside set_boundary_indices)
    %Numerics.Nrlayer(layer): number of radial points per layer (computed inside set_boundary_indices)
    %Numerics.BCindices(layer): indec of the layer point at the boundaries
    %Numerics.perturbation_order: maximum order of the perturbation. Default 2
    %Numerics.rheology_cutoff: determines which terms of the rheology are included (only relevant for viscoelastic). terms with log10(mu_n^m)-log10(mu_n^m(leading))>=-Numerics.rheology_cutoff are included. Default 0 (only leading terms)
    %Numerics.Nenergy: maximum degree to which energy dissipation is expanded. Defaulft 8. 
    %Numerics.load_couplings: 
         % (0) compute coupling coefficients
         % (1) load coupling coefficients from file that contains ALL coupling up to rheology variations of a degree higher than those specied, if such a file does not exist compute it   (default) 
         % (2) load coupling coefficintes from a file that contains specifically the coupling coefficients for the rheology specied, if it does not exisit, compute it. 
    % Numerics.parallel_sol  Calculate the solution using a parfor-loop either 0 or 1
    % Numerics.parallel_gen  Calculate potential coupling files using parfor-loops either 0 or 1

% Interior_Model: struct that contains all the physical data, in case of the 'fixed' method the struct will also be altered.

%% OUTPUT
% Numerics: Changed struct, same as the input struct
% Interior_model: Changed struct (only in case of the 'fixed' method)
%% FUNCTION ------------------------------------------------------------------
function [Numerics, Interior_Model] = set_boundary_indices(Numerics,Interior_Model,varargin)
verbose = 0;
Nrlayer_manual = [];
for k = 1:length(varargin)
    if strcmpi(varargin{k},'verbose')
        verbose=1; 
        varargin{k}=[];
    end
    if strcmpi(varargin{k},'manual')
        Nrlayer_manual=varargin{k+1}; 
        varargin{k+1}=[]; 
        varargin{k}=[];
    end
end

% First check if there is more than 1 active layer
if Numerics.Nlayers > 2
    % Fill Nrlayer, Nr and BCindices based on the selected method
    if Numerics.method == "combination"
        frac_r = [];
        for ilayer=2:Numerics.Nlayers
            frac_r = [frac_r (Interior_Model(ilayer).R0-Interior_Model(ilayer-1).R0)/ ...
                        (Interior_Model(end).R0-Interior_Model(1).R0)];
        end
        Numerics.Nrlayer = [0 floor(frac_r*Numerics.Nrbase)+Numerics.Nrbase];
        Numerics.Nr = sum(Numerics.Nrlayer); % Total number of radial points

        % Set boundary indices
        for ilayer=2:Numerics.Nlayers
            if ilayer == 2
                ind_boundary = 1+Numerics.Nrlayer(ilayer);
            else
                next_boundary = 1+sum(Numerics.Nrlayer(1:(ilayer)));
                ind_boundary = [ind_boundary next_boundary];
            end
        end
        Numerics.BCindices = ind_boundary;

    elseif Numerics.method == "variable" %updated by Marc 
        Numerics.Nrlayer(1)=0; 
        % Set boundary indices
        for ilayer=2:Numerics.Nlayers
            if ilayer == 2
                Numerics.Nrlayer(ilayer)=Numerics.Nrbase;
                Numerics.BCindices(1)= 1+Numerics.Nrlayer(ilayer);
            else
                Numerics.Nrlayer(ilayer)=Numerics.Nrbase;
                Numerics.BCindices(ilayer-1)=1+sum(Numerics.Nrlayer(1:(ilayer)));
            end
        end
        Numerics.Nr = sum(Numerics.Nrlayer); % Total number of radial points, 
    elseif Numerics.method == "fixed"
        delta_r = (Interior_Model(end).R0-Interior_Model(1).R0)/Numerics.Nrbase;
        for ilayer=2:Numerics.Nlayers
            if ilayer == 2
                Numerics.Nrlayer(ilayer)=floor((Interior_Model(ilayer).R0-Interior_Model(1).R0)/delta_r);
                Numerics.BCindices(1)=Numerics.Nrlayer(ilayer)+1; 
                Interior_Model(ilayer).R0=Interior_Model(1).R0 + Numerics.Nrlayer(ilayer)*delta_r;
            elseif ilayer==Numerics.Nlayers
                Numerics.Nrlayer(ilayer)=Numerics.Nrbase-sum(Numerics.Nrlayer(1:(ilayer-1)));
                Numerics.BCindices(ilayer-1)=1+sum(Numerics.Nrlayer(1:(ilayer)));
                Interior_Model(ilayer).R0 = Interior_Model(ilayer-1).R0+Numerics.Nrlayer(ilayer)*delta_r;
            else
                Numerics.Nrlayer(ilayer)=floor((Interior_Model(ilayer).R0-Interior_Model(ilayer-1).R0)/delta_r);
                Numerics.BCindices(ilayer-1)=1+sum(Numerics.Nrlayer(1:(ilayer)));  
                Interior_Model(ilayer).R0 = Interior_Model(ilayer-1).R0+Numerics.Nrlayer(ilayer)*delta_r;               
            end
        end
        Numerics.Nr = sum(Numerics.Nrlayer);
        if ~(Numerics.Nrbase == Numerics.Nr)
            error('Error. \nNumerics.Nr is not equal to Numerics.Nrbase')
        end
        %Numerics.BCindices = ind_boundary;
    
    % Provide a way to manual set the number of points in each layer
    elseif Numerics.method == "manual"
        if isempty(Nrlayer_manual)
            error(['Error. \nAttempted to manually set Nrlayer but was not provided with an array. ' ...
                   'Put manual followed by the array in the varargin to fix this error'])
        end
        Numerics.Nrlayer = Nrlayer_manual;
        Numerics.Nr = sum(Numerics.Nrlayer);

        % Set BCindices based on the provided Nrlayer array
        for ilayer=2:Numerics.Nlayers-1
            if ilayer == 2
                ind_boundary = Numerics.Nrlayer(ilayer);
            else
                next_boundary = sum(Numerics.Nrlayer(1:(ilayer)));
                ind_boundary = [ind_boundary next_boundary];
            end
        end
        Numerics.BCindices = ind_boundary;

    else
        error('Error. \nNot a valid selection of Numerics.method')
    end

% If the model only contains a single layer set Numerics.Nr equal to
% Numerics.Nrbase and Numerics.BCindices is set to empty
else
    Numerics.BCindices = [];
    Numerics.Nrlayer = [0, Numerics.Nrbase];
    Numerics.Nr = Numerics.Nrbase;
end

% Print results if verbose
if verbose == 1
    disp('Inside set_boundary_indices: ')
    disp(['Nrlayer: ' num2str(Numerics.Nrlayer)])
    disp(['BCindices: ' num2str(Numerics.BCindices)])
    disp(['Nr: ' num2str(Numerics.Nr)])
    disp(' ')
end

end