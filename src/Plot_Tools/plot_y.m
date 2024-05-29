%% plot_y
%AUTHOR: M. Rovira-Navarro 
%USE: Plot the y functions

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
        
% y: radial functions
        % y.nf: degree of the forcing
        % y.mf: order of the forcing 
        % y.n: degree of the solution
        % y.m: order of the solution 
        % y.y(radial_point,X,mode) 
            % y(radial_point,1,mode): r radial position
            % y(radial_point,2,mode): U radial displacement
            % y(radial_point,3,mode): V tangential displacement
            % y(radial_point,4,mode): W toroidal displacement 
            % y(radial_point,5,mode): R radial stress
            % y(radial_point,6,mode): S tangential stress
            % y(radial_point,7,mode): T toroidal stress
            % y(radial_point,8,mode): \phi toroidal potential
            % y(radial_point,9,mode): \partial_r\phi potentia derivative       
            % y(radial_point,10,mode): u_{n,n-1}
            % y(radial_point,11,mode): u_{n,n}
            % y(radial_point,12,mode): u_{n,n+1}
            % y(radial_point,13,mode): \sigma_{n,n,0}       
            % y(radial_point,14,mode): \sigma_{n,n-2,2}
            % y(radial_point,15,mode): \sigma_{n,n-1,2}
            % y(radial_point,16,mode): \sigma_{n,n,2}
            % y(radial_point,17,mode): \sigma_{n,n+1,2}
            % y(radial_point,18,mode): \sigma_{n,n+2,2}      
            % y(radial_point,19,mode): \epsilon_{n,n,0}
            % y(radial_point,20,mode): \epsilon_{n,n-2,2}
            % y(radial_point,21,mode): \epsilon_{n,n-1,2}
            % y(radial_point,22,mode): \epsilon_{n,n,2}
            % y(radial_point,23,mode): \epsilon{n,n+1,2}
            % y(radial_point,24,mode): \epsilon_{n,n+2,2}
%% Optional inputs
            % uniform solution,  in case a model with lateral variations is provided, the solution to the spherically-symmetric case is provided 
%% OUTPUT 
% y-functions plots
    % The first plot contains the solution for the spherically-symmetric case
    % The second plot contains the solution including lateral variations


%%

function [] = plot_y(Interior_Model,y,varargin)
uniform_solution=0; 
for k = 1:length(varargin)
    if strcmpi(varargin{k},'uniform solution')
        uniform_solution=1; 
        y_uni=varargin{k+1};
    end
end

nmodes=length(y.n); 
for  k=2:length(Interior_Model)
    if isempty(Interior_Model(k).mu)==0
        mu(k-1)=Interior_Model(k).mu;
    else
        mu(k-1)=0;
    end
   if isempty(Interior_Model(k).Ks)==0
        K(k-1)=Interior_Model(k).Ks;
    else
        K(k-1)=0;
    end
    if isempty(Interior_Model(k).eta)==0
        eta(k-1)=Interior_Model(k).eta;
    else
        eta(k-1)=0;
    end
end
y_limits(1)=max(mu); 
y_limits(2)=max(K); 
y_limits(3)=max(eta);
if isnan(y_limits(3))
    y_limits(3)=1; 
end
y_title={'$\mu_0$','$K_0$','$\eta_0$'...
    '$U_n^m$', '$V_n^m$','$W_n^m$'...
        '$R_n^m$', '$S_n^m$','$T_n^m$'...
        '$\phi_n^m$', '$\partial_r \phi_n^m$'};
fig=figure;
set(fig, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 0.5, 0.9]);
color_layers=customcolormap([0 1],[240, 216, 242; 216, 242, 231]/255,length(Interior_Model));
color_modes=cmocean('phase',nmodes+2); 
color_modes=color_modes(2:end-1,:);
p = uipanel('Parent',fig,'BorderType','none'); 
p.Title = 'Spherically-symmetric'; 
p.TitlePosition = 'centertop'; 
p.FontSize = 18;
p.FontWeight = 'bold';
if uniform_solution==1
    y_p=y_uni.y;
else
    y_p=y.y; 
end
    for j=1:8
            if j==1
            % plot interior model 
            for k=2:length(Interior_Model)
                Rdown=Interior_Model(k-1).R;
                Rup=Interior_Model(k).R;
                %mu 
                subplot(4,3,1,'Parent',p)
                ydown=0;
                yup=abs(y_limits(1));
                pgon = polyshape([Rdown Rdown Rup Rup],[ydown yup yup ydown]);
                if Interior_Model(k).ocean==1
                    color_face=[104, 207, 227]/255;
                else
                    color_face=color_layers(k,:);
                end
                plot(pgon,'FaceColor',color_face,'FaceAlpha',0.5)
                hold on
                plot([Rdown Rup],[mu(k-1) mu(k-1)],'LineWidth',2,'Color','k')
                hold on              
                %K
                subplot(4,3,2,'Parent',p)
                ydown=0;
                yup=abs(y_limits(2));
                pgon = polyshape([Rdown Rdown Rup Rup],[ydown yup yup ydown]);
                if Interior_Model(k).ocean==1
                    color_face=[104, 207, 227]/255;
                else
                    color_face=color_layers(k,:);
                end
                plot(pgon,'FaceColor',color_face,'FaceAlpha',0.5)
                hold on
                plot([Rdown Rup],[K(k-1) K(k-1)],'LineWidth',2,'Color','k')
                hold on
                %eta
                subplot(4,3,3,'Parent',p)
                ydown=0;
                yup=abs(y_limits(3));
                pgon = polyshape([Rdown Rdown Rup Rup],[ydown yup yup ydown]);
                if Interior_Model(k).ocean==1
                    color_face=[104, 207, 227]/255;
                else
                    color_face=color_layers(k,:);
                end
                plot(pgon,'FaceColor',color_face,'FaceAlpha',0.5)
                hold on
                plot([Rdown Rup],[eta(k-1) eta(k-1)],'LineWidth',2,'Color','k')
                hold on  
            end
            end
            % get y limits
      
            y_limits(3+j)=abs(max(abs(y_p(:,1+j,:)),[],"all"));
            if y_limits(3+j)==0
                y_limits(3+j)=1e-4; 
            end
            subplot(4,3,3+j,'Parent',p)
            for k=2:length(Interior_Model)
                hold on
                Rdown=Interior_Model(k-1).R;
                Rup=Interior_Model(k).R;
                ydown=-abs(y_limits(3+j));
                yup=abs(y_limits(3+j));
                pgon = polyshape([Rdown Rdown Rup Rup],[ydown yup yup ydown]);
                if Interior_Model(k).ocean==1
                    color_face=[104, 207, 227]/255;
                else
                    color_face=color_layers(k,:);
                end
                plot(pgon,'FaceColor',color_face,'FaceAlpha',0.5)
                hold on
            end
        subplot(4,3,3+j,'Parent',p)
        plot(y_p(:,1,1),real(y_p(:,1+j,1)),'LineWidth',2,'LineStyle','-','color','k')
        hold on
        plot(y_p(:,1,1),imag(y_p(:,1+j,1)),'LineWidth',2,'LineStyle','--','color','k')
        hold on
    end


for j=1:11
    subplot(4,3,j,'Parent',p)
    title(y_title{j},'Interpreter','latex')
    set(gca,'Box','on');
    set(gca,'fontsize', 18); 
    xlim([Interior_Model(1).R Interior_Model(end).R])
    if j<4
        ylim([0 abs(y_limits(j))])
    else
        ylim([-abs(y_limits(j)) abs(y_limits(j))])
    end
end



% perturbation solution, plot the difference with the uniform solution
if uniform_solution==1
% build interior model to plot 
ind_F=find(y.n==y.nf &  y.m==y.mf);
y_p=[];
y.y(:,2:end,ind_F)=y.y(:,2:end,ind_F)-y_uni.y(:,2:end,:);
y_p(:,1,:)=y.y(:,1,:);
y_p(:,5:12,:)=y.y(:,2:9,:);
for i=1:nmodes
    n_mode=y.n(i); 
    m_mode=y.m(i); 
    for j=2:length(Interior_Model)
        indxmu=find(Interior_Model(j).mu_variable(:,1)==n_mode & Interior_Model(j).mu_variable(:,2)==m_mode);
        indxeta=find(Interior_Model(j).eta_variable(:,1)==n_mode & Interior_Model(j).eta_variable(:,2)==m_mode);
        indxK=find(Interior_Model(j).K_variable(:,1)==n_mode & Interior_Model(j).K_variable(:,2)==m_mode);
        indxRd=find(abs(y_p(:,1,1)-Interior_Model(j-1).R)<1e-10);
        indxRu=find(abs(y_p(:,1,1)-Interior_Model(j).R)<1e-10);
        if isempty(indxmu)==0
            y_p(indxRd:indxRu,2,i)=Interior_Model(j).mu_variable(indxmu,3);
        end
        if isempty(indxeta)==0
            y_p(indxRd:indxRu,3,i)=Interior_Model(j).eta_variable(indxeta,3);
        end
        if isempty(indxK)==0
            y_p(indxRd:indxRu,4,i)=Interior_Model(j).K_variable(indxK,3);
        end
    end
end
% 
y_title={'$\mu_n^m$', '$\eta_n^m$','$K_n^m$'...
        '$U_n^m$', '$V_n^m$','$W_n^m$'...
        '$R_n^m$', '$S_n^m$','$T_n^m$'...
        '$\phi_n^m$', '$\partial_r \phi_n^m$'};
fig=figure;
set(fig, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 0.5, 0.9]);
color_layers=customcolormap([0 1],[240, 216, 242; 216, 242, 231]/255,length(Interior_Model));
color_modes=cmocean('phase',nmodes+2); 
color_modes=color_modes(2:end-1,:);
tcl = tiledlayout(4,3);
title(tcl,'Perturbation arising from lateral variations','FontSize',18)
for i=1:nmodes
    leg{i}=[ '$(n,m)=(' num2str(y.n(i)) ',' num2str(y.m(i)) ')$'];
    for j=1:11
        if i==1 %plot layers
            % get y limits
            y_limits(j)=abs(max(abs(y_p(:,1+j,:)),[],"all"));
            if y_limits(j)==0
                y_limits(j)=1e-4; 
            end
            % plot interior model 
            nexttile(j);
            for k=2:length(Interior_Model)
                hold on
                Rdown=Interior_Model(k-1).R;
                Rup=Interior_Model(k).R;
                ydown=-abs(y_limits(j));
                yup=abs(y_limits(j));
                pgon = polyshape([Rdown Rdown Rup Rup],[ydown yup yup ydown]);
                if Interior_Model(k).ocean==1
                    color_face=[104, 207, 227]/255;
                else
                    color_face=color_layers(k,:);
                end
                plot(pgon,'FaceColor',color_face,'FaceAlpha',0.5)
                hold on
            end
        end
        nexttile(j);
        if j==1
            leg_p(i)=plot(y_p(:,1,i),real(y_p(:,1+j,i)),'LineWidth',2,'LineStyle','-','color',color_modes(i,:));
        else
            plot(y_p(:,1,i),real(y_p(:,1+j,i)),'LineWidth',2,'LineStyle','-','color',color_modes(i,:));
        end
        hold on
        plot(y_p(:,1,i),imag(y_p(:,1+j,i)),'LineWidth',2,'LineStyle','--','color',color_modes(i,:));
        hold on
    end
end

for j=1:11
    nexttile(j);
    title(y_title{j},'Interpreter','latex')
    set(gca,'Box','on');
    set(gca,'fontsize', 18); 
    xlim([Interior_Model(1).R Interior_Model(end).R])
    ylim([-abs(y_limits(j)) abs(y_limits(j))])
end
hL = legend(leg_p,leg,'interpreter','latex','fontsize',16); 
hL.Layout.Tile = 12; %
end

end