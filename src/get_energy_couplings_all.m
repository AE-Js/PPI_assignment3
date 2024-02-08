%% Functions used to compute all coupling coefficintes
%% INPUT
    % perturbation_order: maximum perturbation order considered
    % Nrheo_max: maximum degree to which rheology is expanded 
    % Forcing: Forcing for which the coupling coefficients are computed
%% OUTPUT
    % Couplings: Structure containing all the Couplings, see get_couplings
    % for more details about the structure. 
%% Function
function [ECouplings]=get_energy_couplings_all(perturbation_order,Nrheo_max,N_en_max,Forcing,Numerics,varargin)
verbose=0;
for k = 1:length(varargin)
    if strcmpi(varargin{k},'verbose')
        verbose=1; 
        varargin{k}=[];
    end
end
%% compute n_en 
k=1;
for i=0:N_en_max
    for j=-i:i
        n_en(k)=i; 
        m_en(k)=j; 
        k=k+1; 
    end   
end

%% compute n_sol and m_sol for a given rheology and perturbation order
% Order of perturbations in the solution
max_order = perturbation_order;

% Forcing modes
Nf=length(Forcing);
for ii=1:Nf
    Forcing_internal.n(ii) = Forcing(ii).n;
    Forcing_internal.m(ii) = Forcing(ii).m;
end

% initialize rheology structure
mu_v=zeros((Nrheo_max+1)*(Nrheo_max+1)-1,2);
jj=1;
for l=1:Nrheo_max
    for m=0:l
        if m==0
            mu_v(jj,1)=l; 
            mu_v(jj,2)=m;
            jj=jj+1; 
        else
            mu_v(jj,1)=l; 
            mu_v(jj,2)=m; 
            mu_v(jj+1,1)=l;
            mu_v(jj+1,2)=-m;
            jj=jj+2; 
        end
    end
end

% Rename some things
variations = mu_v;

% Loop over the three forcings
for iforce=1:length(Forcing_internal.n)
    modes4 = get_active_modes(max_order,variations,Forcing(iforce));
    
    % Store the degree,order and order of perturbation for each forcing
    Coupstore(iforce).n_s = modes4(:,1);
    Coupstore(iforce).m_s = modes4(:,2);
    Coupstore(iforce).order = modes4(:,3);
end

% Retrieve all the activated modes from all forcings, remove duplicates
Nf=length(Forcing_internal.n); 
n_sol=[];
m_sol=[];
for i=1:Nf
    n_sol=[n_sol Coupstore(i).n_s'];
    m_sol=[m_sol Coupstore(i).m_s'];
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
n_sol=n_sol_t; 
m_sol=m_sol_t;  
N_sol = length(n_sol);

%% initialize matrix
tic
EC=zeros(length(n_sol),length(n_sol),length(n_en),6,6);
la=[2 2 2 2 2 0];
lb=[2 2 2 2 2 0];
la_vec = repelem(la,54)'; % 54 comes from 324/6, 324 comes from 3x3x6x6
lb_vec = repmat(repelem(lb,9),1,6)'; % 9 comes from number of n1
N_en = length(n_en);

% Check whether a parallel method can be used
if Numerics.parallel_gen == 1
    parfor i1=1:N_sol %loop modes 1
        for i2=1:N_sol
            %(1) loop over modes that result from coupling
            % Selection rule (2). Triangular inequality.
            n_aux=[abs(n_sol(i1)-n_sol(i2)) n_sol(i1)+n_sol(i2)]; % Changed to allow for lower Nmax
    
            % Selection rule (7)
            m_aux=m_sol(i1)+m_sol(i2); 
            ind_resulting_modes=find(n_en>=n_aux(1) & n_en<=n_aux(end) & m_en==m_aux); 
            
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
                    nc=zeros(324,1)+n_en(ind_res_mode);
        
                    % Generate the n2 values based on na_val and make complete vector
                    na2=[na_val-2:na_val+2 na_val];
                    nb2=[nb_val-2:nb_val+2 nb_val];
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
                    na=zeros(n_non_zero,1)+n_sol(i1);
                    ma=zeros(n_non_zero,1)+m_sol(i1);
                    nb=zeros(n_non_zero,1)+n_sol(i2);
                    mb=zeros(n_non_zero,1)+m_sol(i2); 
                    mc=zeros(n_non_zero,1)+m_en(ind_res_mode);
        
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
            if verbose ==1
                if mod(i2,25) == 0
                    % clc
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
            n_aux=[abs(n_sol(i1)-n_sol(i2)) n_sol(i1)+n_sol(i2)]; % Changed to allow for lower Nmax
    
            % Selection rule (7)
            m_aux=m_sol(i1)+m_sol(i2); 
            ind_resulting_modes=find(n_en>=n_aux(1) & n_en<=n_aux(end) & m_en==m_aux); 
            
            %(2) loop over modes that result from coupling 
            for i3=1:length(ind_resulting_modes)
                % Allocate values for na and nb
                na_val = n_sol(i1);
                nb_val = n_sol(i2);
    
                % Fill complete nc vector (resulting energy modes)
                nc=zeros(324,1)+n_en(ind_resulting_modes(i3));
    
                % Generate the n2 values based on na_val and make complete vector
                na2=[na_val-2:na_val+2 na_val];
                nb2=[nb_val-2:nb_val+2 nb_val];
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
                na=zeros(n_non_zero,1)+n_sol(i1);
                ma=zeros(n_non_zero,1)+m_sol(i1);
                nb=zeros(n_non_zero,1)+n_sol(i2);
                mb=zeros(n_non_zero,1)+m_sol(i2); 
                mc=zeros(n_non_zero,1)+m_en(ind_resulting_modes(i3));
    
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
%                     clc
                    disp([ num2str(N_sol) ' Solution modes'])
                    disp(['Progress: ' num2str(((i1-1)*N_sol + i2)/N_sol^2 * 100) ' % Iterations done'])
                end
            end
        end
    end
end
disp(["Time spent on all calculating energy modes: " num2str(toc) " s"])

%% find the non-zero coefficients 
k=1;
non_zero_index=[]; 
n_en2=[];
m_en2=[];
for i=1:length(n_en)
    EC_deg=abs(EC(:,:,i,:,:));
    max_c=max(EC_deg(:));
    if max_c>0
        n_en2(k)=n_en(i);
        m_en2(k)=m_en(i); 
        non_zero_index(k)=i;
        k=k+1; 
    end
end

ECouplings.EC = EC(:,:,non_zero_index,:,:);%single(EC(:,:,non_zero_index,:,:)); 
ECouplings.n_s=n_sol;%single(n_sol);
ECouplings.m_s=m_sol;%single(m_sol);
ECouplings.n_en=n_en2;%single(n_en2);
ECouplings.m_en=m_en2;%single(m_en2);

end

%% explicit computation of the coupling coefficients
function [C] = couplings_coefficient(na,na2,la,ma,nb,nb2,lb,mb,nc,mc,na1,nb1)

one_arr = zeros(324,1)+1;
zero_arr = zeros(324,1);

Lama=sqrt((2*la+1).*(2*na1+1));
Lamb=sqrt((2*lb+1).*(2*nb1+1));
CC=(-1).^(mc+nb+nb2).* ...
sqrt((2*na2+1).*(2*na1+1).*(2*na+1)).* ...
sqrt((2*nb2+1).*(2*nb1+1).*(2*nb+1)).*sqrt(2*nc+1).* ...
Wigner3j(na2, nb2, nc, zero_arr, zero_arr, zero_arr).* ...
Wigner3j(na,nb,nc,ma,mb,-mc).* ...
Wigner9j(na,na1,one_arr,nc,na2,nb2,nb,one_arr,nb1);
Caux=(-1).^(na+na2+la+nb+nb2+lb).*Lama.*Lamb.*Wigner6j(one_arr, la, one_arr, na, na1, na2).*Wigner6j(one_arr, lb, one_arr, nb, nb1, nb2).*CC;
C_int = reshape(Caux,9,[]);
C = sum(C_int);
end
