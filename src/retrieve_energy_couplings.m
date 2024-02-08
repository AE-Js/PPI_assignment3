%% Description
% Extracts the relevant couplings for the given solution modes and requested energy modes
% This function retrieves the relevant energy couplings when the
% coupling coefficient file has just been calculated instead of loaded from
% a pre-computed file.
%%
function [Couplings] = retrieve_energy_couplings(n_sol,m_sol,Nen,Couplings)
% Define the indexes vectors
sol_indexes=zeros(length(n_sol),1);

% Find the Solution modes indexes
for i=1:length(n_sol)
    sol_indexes(i)=find(n_sol(i)==Couplings.n_s & m_sol(i)==Couplings.m_s);
end

% find the maximum possible Energy modes indexes. 
% (This already filters out some energy modes that have a zero coupling term)
Nen_indexes = find(Couplings.n_en <= Nen);

% Find all the energy modes for which there are non-zero coefficients 
k=1;
non_zero_index=[]; 
n_en=[];
m_en=[];
for i=1:length(Nen_indexes)
    EC_deg=abs(Couplings.EC(sol_indexes,sol_indexes,i,:,:));
    max_c=max(EC_deg(:));
    if max_c>0
        n_en(k)=Couplings.n_en(i);
        m_en(k)=Couplings.m_en(i); 
        non_zero_index(k)=i;
        k=k+1; 
    end
end

% Test for good selection
if any(~(n_en == Couplings.n_en(non_zero_index))) && any(~(m_en == Couplings.m_en(non_zero_index)))
    error('Error. \nInside retrieve energy_couplings.m \nProblems with consistency')
end

% Select the relevant coefficients
Couplings.n_s=Couplings.n_s(sol_indexes);
Couplings.m_s=Couplings.m_s(sol_indexes);
Couplings.n_en=Couplings.n_en(non_zero_index);
Couplings.m_en=Couplings.m_en(non_zero_index);
Couplings.EC=Couplings.EC(sol_indexes,sol_indexes,non_zero_index,:,:);

end

