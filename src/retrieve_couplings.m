%% Description
% Extracts the relevant couplings for the given rheology, forcing modes, 
% and solution modes that are based on the perturbation depth. 
% This function retrieves the relevant couplings when the
% coupling coefficient file has just been calculated instead of loaded from
% a pre-computed file.

%% FUNCTION
function [Couplings] = retrieve_couplings(max_order,variations,Forcing,Couplings)
% for the given rheology compute the modes that are excited using the
% selection rules
modes4 = get_active_modes(max_order,variations,Forcing);

n_s = modes4(:,1);
m_s = modes4(:,2);
order = modes4(:,3);

modes_indexes = length(n_s); 
rheo_indexes = size(variations,1);

% find the rheology indexes
for i=1:size(variations,1)
    rheo_indexes(i) = find(variations(i,1)==Couplings.n_r & variations(i,2)==Couplings.m_r);
end

% find the modes indexes. 
for i=1:length(n_s)
    modes_indexes(i) = find(n_s(i)==Couplings.n_s & m_s(i)==Couplings.m_s);
end

% select coefficients
Couplings.n_s = Couplings.n_s(modes_indexes);
Couplings.m_s = Couplings.m_s(modes_indexes);
Couplings.order = order; 
Couplings.n_r = variations(:,1);
Couplings.m_r = variations(:,2);
Couplings.Coup = Couplings.Coup(modes_indexes,modes_indexes,:,rheo_indexes);

end

