%% Functions used to compute all coupling coefficintes
%% INPUT
    % perturbation_order: maximum perturbation order considered
    % Nrheo_max: maximum degree to which rheology is expanded 
    % Forcing: Forcing for which the coupling coefficients are computed
%% OUTPUT
    % Couplings: Structure containing all the Couplings, see get_couplings
    % for more details about the structure. 
function [Couplings]=get_couplings_all(perturbation_order,Nrheo_max,Forcing,Numerics,varargin)
verbose=0;
for k = 1:length(varargin)
    if strcmpi(varargin{k},'verbose')
        verbose=1; 
        varargin{k}=[];
    end
end
%% initialize rheology structure
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

%% compute couplings 
if verbose==1
    Couplings = get_couplings(perturbation_order,mu_v,Nrheo_max,Forcing,Numerics,'verbose');
else
    Couplings = get_couplings(perturbation_order,mu_v,Nrheo_max,Forcing,Numerics);
end
Couplings.n_r = mu_v(:,1);%single(mu_v(:,1));
Couplings.m_r = mu_v(:,2);%single(mu_v(:,2));
end