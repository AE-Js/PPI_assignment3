%% Description 
% Computes all the activated modes based on the rheology and forcing modes 
% as well the perturbation depth/order that is requested. While computing
% all the modes it already filters duplicates in order for the computation
% to go faster.
%%
function [active_modes] = get_active_modes(max_order,variations,Forcing)
% for the given rheology compute the modes that are excited using the
% selection rules
rheo = unique((variations(:,[1 2])),'rows');

% Generate modes matrix (is oversized) for first order
modes=zeros( size(rheo,1)*6, 4);
order=0; 

% Start with the forcing 
mode0(1)=Forcing.n;
mode0(2)=Forcing.m;
mode0(3)=1;
mode0(4)=0; 
modes(1,:)=mode0;
index_old=1;
order=order+1;

% The modes array has the following structure:
% 1: degree of mode, 2: order of mode, 3: spheroidal or toroidal, 4: the order

% Find all the modes using the next_coupling function
while order<max_order+1
    index_new = index_old(end); 
    for k=1:size(rheo,1)
        modeR(1)=rheo(k,1);
        modeR(2)=rheo(k,2);
        for j=1:length(index_old)
            % Calculate the activated modes based on selection rules
            modes_aux = next_coupling(modes(index_old(j),:),modeR,order);
            size_val = size(modes_aux,1);

            % Put activated modes into array
            modes((index_new+1):(index_new+size_val),:) = modes_aux;
            index_new = index_new + size_val;
        end
    end
    
    % Trim the excess zeros of the modes array
    ind = find(modes(:,3));
    modes = modes(ind,:);

    % Retrieve all the unique active modes
    modes2 = unique(modes(:,[4, 1:3]),'rows');
    modes3 = unique(modes(:,[1:2]),'rows');

%     % Remove artefacts from the zeros array - NOT NEEDED ANYMORE
%     modes2(modes2(:,end) == 0,:) = [];
%     if find(modes2(:,2) == 0 & modes(:,3) == 0)
%         
%     else
%         ind = find(modes3(:,1) == 0 & modes3(:,2) == 0);
%         modes3(ind,:) = [];
%     end

    % Find potential other forcing mode
    ind_f = find(modes2(2:end,2)==Forcing.n & modes2(2:end,3)==Forcing.m);
    if isempty(ind_f)==0 && modes2(1,1)==0
        modes2(1,1) = modes2(ind_f(1)+1,1);
    end

    % Combine active modes into one array
    modes4 = zeros(size(modes3,1),4);
    for j=1:length(modes3)
        ind = find(modes2(:,2)==modes3(j,1) & modes2(:,3)==modes3(j,2));
        
        % In case there are two identical modes select based on order
        [val,ind_2] = min(modes2(ind,1));
        pot_arr = modes2(ind,end);

        % Populate the intermediate active modes array
        modes4(j,:) = [modes3(j,:), pot_arr(ind_2), val];
    end

    % Create new matrix big enough to contain all modes from the next order
    modes = zeros( size(rheo,1)*length(modes2)*20, 4);

    % Copy over the activated modes
    modes(1:size(modes4,1),:) = modes4;
    
    % Go to next order
    index_old = find(modes4(:,4) == order);
    order=order+1; 
end

% Trim the excess zeros of the modes array
ind = find(modes(:,3));
modes = modes(ind,:);

% Retrieve all the active modes
modes2=unique(modes(:,[4, 1:3]),'rows');
modes3=unique(modes(:,[1:2]),'rows');
ind_f = find(modes2(2:end,2)==Forcing.n & modes2(2:end,3)==Forcing.m);
if isempty(ind_f)==0 && modes2(1,1)==0
        modes2(1,1) = modes2(ind_f(1)+1,1);
end
%ind_f=find(modes2(2:end,2)==Forcing.n & modes2(2:end,3)==Forcing.m);
% if isempty(ind_f)==0
%     modes2(ind_f(1)+1,1)=modes2(ind_f(1)+1,1);
% end

active_modes = zeros(length(modes3),3);
for j=1:length(modes3)
    ind = find(modes2(:,2)==modes3(j,1) & modes2(:,3)==modes3(j,2));
    active_modes(j,:) = [modes3(j,:), min(modes2(ind,1))];
end

end