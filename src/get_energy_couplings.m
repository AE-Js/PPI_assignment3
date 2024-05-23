%% Function used to compute the couplings for energy dissipation 
%% INPUT 
% n_sol: degrees of modes involved in the solution \alpha
% ------------------
% m_sol: orders involved in the solution \alpha
% ------------------
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
% varargin: can be 'verbose' if output of text is desired
%% OUTPUT
% ECouplings: Struct containing couplingmatrix and list of modes
    % EC(n_mode1,n_mode2,n_mode3,n2_mode1, n2_mode2):
    %   with n_mode1: mode 1 of the solution with n_s(n_mode1), m_s(n_mode1)
    %   with n_mode2: mode 2 of the solution with n_s(n_mode1), m_s(n_mode1)
    %   with n_mode3: mode of the energy spectra with n_en(n_mode1), m_en(n_mode1)
    %   with n2_mode1: n2 of mode 1
    %   with n2_mode2: n2 of mode 2
    % n_s: degrees of the solution
    % m_s: orders of the solution
    % n_en: degrees of the energy spectrum 
    % m_en: orders of the energy spectrum 

%% FUNCTION
function [ECouplings] = get_energy_couplings(n_sol,m_sol,Numerics,varargin)
verbose = 0;
for k = 1:length(varargin)
    if strcmpi(varargin{k},'verbose')
        verbose = 1; 
        varargin{k}=[];
    end
end

%% Compute n_en 
k = 1;
for i=0:Numerics.Nenergy
    for j=-i:i
        n_en(k) = i; 
        m_en(k) = j; 
        k = k+1; 
    end   
end

%% Initialize matrices
EC = zeros(length(n_sol),length(n_sol),length(n_en),6,6);
la = [2 2 2 2 2 0];
lb = [2 2 2 2 2 0];
la_vec = repelem(la,54)'; % 54 comes from 324/6, 324 comes from 3x3x6x6
lb_vec = repmat(repelem(lb,9),1,6)'; % 9 comes from number of n1
N_en = length(n_en);
N_sol = length(n_sol);

%% Fill coupling coefficient matrix
if verbose == 1
    tic
end

% Check whether a parallel method can be used
if Numerics.parallel_gen == 1
    parfor i1=1:N_sol %loop modes 1
        for i2=1:N_sol
            %(1) loop over modes that result from coupling
            % Selection rule (2). Triangular inequality.
            n_aux = [abs(n_sol(i1)-n_sol(i2)) n_sol(i1)+n_sol(i2)]; % Changed to allow for lower Nmax
    
            % Selection rule (7)
            m_aux = m_sol(i1) + m_sol(i2); 
            ind_resulting_modes = find(n_en>=n_aux(1) & n_en<=n_aux(end) & m_en==m_aux); 
            
            k = 1;
            %(2) loop over modes that result from coupling 
            for i3=1:N_en
                if k > length(ind_resulting_modes)
                    break
                end
                if i3 == ind_resulting_modes(k)
                    ind_res_mode = ind_resulting_modes(k);
        
                    % Allocate values for na and nb
                    na_val = n_sol(i1);
                    nb_val = n_sol(i2);
        
                    % Fill complete nc vector (resulting energy modes)
                    nc = zeros(324,1) + n_en(ind_res_mode);
        
                    % Generate the n2 values based on na_val and make complete vector
                    na2 = [na_val-2:na_val+2 na_val];
                    nb2 = [nb_val-2:nb_val+2 nb_val];
                    na2_vec = repelem(na2,54)';
                    nb2_vec = repmat(repelem(nb2,9),1,6)';
                    
                    % Apply selection rule 3 and 4 to reduce the number of calculations
                    non_zero_elems = find( mod(na2_vec+nb2_vec+nc,2)==0 & nb2_vec>=abs(na2_vec-nc) & nb2_vec<=abs(na2_vec+nc) );
                    
                    % Find row and column indices of non zero elements for later
                    i4_vec = ceil(non_zero_elems/54);
                    i5_vec = mod(ceil((non_zero_elems)/9),6);
                    i5_vec(i5_vec==0) = 6; % Fix modulo "mistake"
                    i4_vec = i4_vec(1:9:end);
                    i5_vec = i5_vec(1:9:end);
                    
                    % Fill fixed value arrays using the number of non zero elements
                    n_non_zero = length(non_zero_elems);
                    na = zeros(n_non_zero,1) + n_sol(i1);
                    ma = zeros(n_non_zero,1) + m_sol(i1);
                    nb = zeros(n_non_zero,1) + n_sol(i2);
                    mb = zeros(n_non_zero,1) + m_sol(i2); 
                    mc = zeros(n_non_zero,1) + m_en(ind_res_mode);
        
                    % Remove zero elements from relevant vectors
                    nc = nc(non_zero_elems);
                    na2_vec = na2_vec(non_zero_elems);
                    nb2_vec = nb2_vec(non_zero_elems);
                    la_vec2 = la_vec(non_zero_elems);
                    lb_vec2 = lb_vec(non_zero_elems);
        
                    % Generate n1 values (previously inside another function)
                    na1v = [na_val-1 na_val na_val+1];
                    nb1v = [nb_val-1 nb_val nb_val+1];
                    
                    % Convert n1 array to proper size by repeating elements/vectors
                    repna1 = floor(n_non_zero/9);
                    repnb1 = floor(n_non_zero/3);
                    na1 = repmat(repelem(na1v,3),1,repna1)';
                    nb1 = repmat(nb1v,1,repnb1)';
                    
                    % Calculate coupling coefficients
                    EC_vec = couplings_coefficient(na,na2_vec,la_vec2,ma,nb,nb2_vec,lb_vec2,mb,nc,mc,na1,nb1);
        
                    % Fill a (6,6) n2 matrix with the obtained values
                    store_mat = zeros(36,1);
                    ind = sub2ind([6 6],i4_vec,i5_vec);
                    store_mat(ind) = EC_vec;
        
                    % Put the values in the EC matrix
                    EC(i1,i2,i3,:,:) = reshape(store_mat,[6,6]);
                    k = k+1;
                end
            end
    
            % Print progress update
            if verbose == 1
                if mod(i2,25) == 0
                    disp([ num2str(N_sol) ' Solution modes'])
                    disp(['Solution mode 1: ' num2str(i1) ', Progress: ' num2str((i2)/N_sol * 100) ' % Iterations done'])
                end
            end
        end
    end
else
    for i1=1:N_sol %loop modes 1
        for i2=1:N_sol
            %(1) loop over modes that result from coupling
            % Selection rule (2). Triangular inequality.
            n_aux = [abs(n_sol(i1)-n_sol(i2)) n_sol(i1)+n_sol(i2)]; % Changed to allow for lower Nmax
    
            % Selection rule (7)
            m_aux = m_sol(i1)+m_sol(i2); 
            ind_resulting_modes = find(n_en>=n_aux(1) & n_en<=n_aux(end) & m_en==m_aux); 
            
            %(2) loop over modes that result from coupling 
            for i3=1:length(ind_resulting_modes)
                % Allocate values for na and nb
                na_val = n_sol(i1);
                nb_val = n_sol(i2);
    
                % Fill complete nc vector (resulting energy modes)
                nc = zeros(324,1)+n_en(ind_resulting_modes(i3));
    
                % Generate the n2 values based on na_val and make complete vector
                na2 = [na_val-2:na_val+2 na_val];
                nb2 = [nb_val-2:nb_val+2 nb_val];
                na2_vec = repelem(na2,54)';
                nb2_vec = repmat(repelem(nb2,9),1,6)';
                
                % Apply selection rule 3 and 4 to reduce the number of calculations
                non_zero_elems = find( mod(na2_vec+nb2_vec+nc,2)==0 & nb2_vec>=abs(na2_vec-nc) & nb2_vec<=abs(na2_vec+nc) );
                
                % Find row and column indices of non zero elements for later
                i4_vec = ceil(non_zero_elems/54);
                i5_vec = mod(ceil((non_zero_elems)/9),6);
                i5_vec(i5_vec==0) = 6; % Fix modulo "mistake"
                i4_vec = i4_vec(1:9:end);
                i5_vec = i5_vec(1:9:end);
                
                % Fill fixed value arrays using the number of non zero elements
                n_non_zero = length(non_zero_elems);
                na = zeros(n_non_zero,1) + n_sol(i1);
                ma = zeros(n_non_zero,1) + m_sol(i1);
                nb = zeros(n_non_zero,1) + n_sol(i2);
                mb = zeros(n_non_zero,1) + m_sol(i2); 
                mc = zeros(n_non_zero,1) + m_en(ind_resulting_modes(i3));
    
                % Remove zero elements from relevant vectors
                nc = nc(non_zero_elems);
                na2_vec = na2_vec(non_zero_elems);
                nb2_vec = nb2_vec(non_zero_elems);
                la_vec2 = la_vec(non_zero_elems);
                lb_vec2 = lb_vec(non_zero_elems);
    
                % Generate n1 values (previously inside another function)
                na1v = [na_val-1 na_val na_val+1];
                nb1v = [nb_val-1 nb_val nb_val+1];
                
                % Convert n1 array to proper size by repeating elements/vectors
                repna1 = floor(n_non_zero/9);
                repnb1 = floor(n_non_zero/3);
                na1 = repmat(repelem(na1v,3),1,repna1)';
                nb1 = repmat(nb1v,1,repnb1)';
                
                % Calculate coupling coefficients
                EC_vec = couplings_coefficient(na,na2_vec,la_vec2,ma,nb,nb2_vec,lb_vec2,mb,nc,mc,na1,nb1);
    
                % Fill a (6,6) n2 matrix with the obtained values
                store_mat = zeros(36,1);
                ind = sub2ind([6 6],i4_vec,i5_vec);
                store_mat(ind) = EC_vec;
    
                % Put the values in the EC matrix
                EC(i1,i2,ind_resulting_modes(i3),:,:) = reshape(store_mat,[6,6]);
            end
    
            % Print progress update
            if verbose ==1
                if mod(i2,25) == 0
                    disp([ num2str(N_sol) ' Solution modes'])
                    disp(['Progress: ' num2str(((i1-1)*N_sol + i2)/N_sol^2 * 100) ' % Iterations done'])
                end
            end
        end
    end
end
if verbose == 1
    disp(["Time spent on all calculating energy modes: " num2str(toc) " s"])
end

%% find the non-zero coefficients 
k = 1;
non_zero_index = []; 
n_en2 = [];
m_en2 = [];
for i=1:length(n_en)
    EC_deg = abs(EC(:,:,i,:,:));
    max_c = max(EC_deg(:));
    if max_c>0
        n_en2(k) = n_en(i);
        m_en2(k) = m_en(i); 
        non_zero_index(k) = i;
        k = k+1; 
    end
end

ECouplings.EC = EC(:,:,non_zero_index,:,:);%single(EC(:,:,non_zero_index,:,:)); 
ECouplings.n_s = n_sol;%single(n_sol);
ECouplings.m_s = m_sol;%single(m_sol);
ECouplings.n_en = n_en2;%single(n_en2);
ECouplings.m_en = m_en2;%single(m_en2);

end

%% explicit computation of the coupling coefficients
function [C] = couplings_coefficient(na,na2,la,ma,nb,nb2,lb,mb,nc,mc,na1,nb1)

one_arr = zeros(324,1)+1;
zero_arr = zeros(324,1);

Lama = sqrt((2*la+1).*(2*na1+1));
Lamb = sqrt((2*lb+1).*(2*nb1+1));
CC = (-1).^(mc+nb+nb2).* ...
     sqrt((2*na2+1).*(2*na1+1).*(2*na+1)).* ...
     sqrt((2*nb2+1).*(2*nb1+1).*(2*nb+1)).*sqrt(2*nc+1).* ...
     Wigner3j(na2, nb2, nc, zero_arr, zero_arr, zero_arr).* ...
     Wigner3j(na,nb,nc,ma,mb,-mc).* ...
     Wigner9j(na,na1,one_arr,nc,na2,nb2,nb,one_arr,nb1);
Caux = (-1).^(na+na2+la+nb+nb2+lb).*Lama.*Lamb.*Wigner6j(one_arr, la, one_arr, na, na1, na2).*Wigner6j(one_arr, lb, one_arr, nb, nb1, nb2).*CC;
C_int = reshape(Caux,9,[]);
C = sum(C_int);
end
