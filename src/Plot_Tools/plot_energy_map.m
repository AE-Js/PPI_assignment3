%% PLOT_ENERGY_MAP
% function used to plot map of energy dissipatoion
%% INPUT
%Energy_Spectra: 
        %Energy_Spectra.n: degrees with non-zero energy
        %Energy_Spectra.m: orders with non-zero energy 
        %Energy_Spectra.n_v: degrees from 0 to Numerics.Nenergy 
        %Energy_Spectra.n_v: orders from 0 to Numerics.Nenergy 
        %Energy_Spectra.energy(radial_point,mode): radial profile of energy spectra
        %Energy_Spectra.energy_integral(mode): radially integrated energy for all non-zero degrees an orders (n,m)
        %Energy_Spectra.energy_integral_v(mode): radially integrated energy for all degrees an orders (n_v,m_v)

% optional variables
    % 'type': 'total' or 'difference'
    % 'save_plot': name of the file where plot is saved 
    % 'limits': specify the limits of the colorscale 
    % 'label': label of the colorbar
    % 'projection': type of projection used for the plot, 'mollweide' is the default and accepts
        % mollweide, 
        % flat 
        % cut: flat plus includes integral along longitude and latitude
%% OUTPUT 

%%
function [varargout] = plot_energy_map(Energy_Spectra,varargin)
%% TOTAL TIDAL DISSIPATION
n_v=Energy_Spectra.n_v;
m_v=Energy_Spectra.m_v;
energy_s=Energy_Spectra.energy_integral_v;
limits=[];
cut=0; 
projection='mollweide';
label_z=[];
no_colorbar=0; 
no_grid=0; 
type='total';
title_plot='';
save_plot=0;
for k = 1:length(varargin)
    if strcmpi(varargin{k},'limits')
        limits=varargin{k+1}; 
        varargin{k+1}=[]; 
        varargin{k}=[];
    end
    if strcmpi(varargin{k},'save_plot')
        save_plot = 1;
        save_name=varargin{k+1};  
        varargin{k}=[];
    end
    if strcmpi(varargin{k},'label')
        label_z=varargin{k+1}; 
        varargin{k+1}=[]; 
        varargin{k}=[];
    end
    if strcmpi(varargin{k},'type')
        type=varargin{k+1}; 
        varargin{k+1}=[]; 
        varargin{k}=[];
    end
    if strcmpi(varargin{k},'no_colorbar')
        no_colorbar=1; 
        varargin{k}=[];
    end
    if strcmpi(varargin{k},'no_grid')
        no_grid=1; 
        varargin{k}=[];
    end
    if strcmpi(varargin{k},'title')
        title_plot=varargin{k+1}; 
        varargin{k+1}=[]; 
        varargin{k}=[];
    end
    if  strcmpi(varargin{k},'projection')
        projection=varargin{k+1};
        varargin{k+1}=[]; 
        varargin{k}=[];
    end
end
%% Get Map
th = linspace(0,pi,500);    % inclination
phi = linspace(-pi,pi,500); % azimuth
[th,phi] = meshgrid(th,phi); 
EnergyV=zeros(size(th)); 
for i=1:length(n_v)
    if abs(energy_s(i))>0
        EnergyV=EnergyV+energy_s(i)*harmonicY(n_v(i),m_v(i),th,phi);
    end
end
aux2=real(EnergyV);
% get colobar limits
if strcmpi(type,'total')
    max_E=max(aux2(:));
    min_E=min(aux2(:)); 
else
    max_E=max(aux2(:));
    min_E=min(aux2(:));  
    max_E2=max([abs(max_E), abs(min_E)]);
end
%% Integrate energy along longitude and latitude 
th_cut=th(1,:); 
phi_cut=phi(:,1); 
Delta_th=th_cut(2)-th_cut(1);
Delta_phi=phi_cut(2)-phi_cut(1);
Energy_Lat=zeros(1,length(th_cut));
Energy_Lon=zeros(1,length(th_cut));
 for i=1:length(th_cut)
     % keep latitude constant, integrate longitude
     Energy_Lat(i)=sum(real(EnergyV(:,i))*Delta_phi);
     % keep longitude constant, integrate latitude
     Energy_Lon(i)=sum(real(EnergyV(i,:)).*sin(th_cut)*Delta_th);
 end
varargout{1}=Energy_Lat; 
varargout{2}=Energy_Lon;
sum(Energy_Lat.*sin(th_cut)*Delta_th)/(4*pi);
sum(Energy_Lon*Delta_phi)/(4*pi);
%% Make plot
if strcmpi(projection,'mollweide') 
    fig=figure;
    set(fig, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 0.6, 0.95]);
    m_proj('mollweide','clongitude',0);
    set(0,'defaulttextInterpreter','none') 
    obj=m_pcolor(phi*180/pi,90-th*180/pi,real(EnergyV)); 
    set(0,'defaulttextInterpreter','latex') 
    shading interp;
    if strcmpi(type,'total')
        colormap(cmocean('solar',100))
    else
        colormap(cmocean('balance',100))
    end
    set(0,'defaulttextInterpreter','none')
    if no_grid==0
        m_grid('xaxislocation','middle','xtick',[-90:90:180],'ytick',[-90:90:90],'fontsize',20,'gridcolor',[0.8 0.8 0.8],'linewidth',3,'linestyle','--','gridcolor','k');
    else
        m_grid('xaxislocation','middle','xtick',[-90:90:180],'ytick',[-90:90:90],'gridcolor',[0.8 0.8 0.8],'linewidth',3,'linestyle','--','gridcolor','k','fontsize',1);
    end
    if no_colorbar==0
    cb=colorbar('southoutside');
    set(0,'defaulttextInterpreter','latex')
    cb.Label.Interpreter = 'latex';
    if isempty(label_z)==1
        if strcmpi(type,'total')
            cb.Label.String='$\dot{e}$ [-]';
        else
            cb.Label.String='$\Delta\dot{e}$ [-]';
        end
    else
        cb.Label.String=label_z; 
        cb.FontSize=50;
    end
    cb.FontName = 'CMU Serif';
    end
    set(gca,'fontsize', 50);
    set(gcf,'color','w');
     if isempty(limits)==0
        caxis([limits(1) limits(2)])
     else
         if strcmpi(type,'total')
            caxis([0 max_E])
        else
           caxis([-max_E2 max_E2])
        end
     end
     if isempty(title_plot)==0
        title(title_plot,'interpreter','latex')
    end
elseif strcmpi(projection,'flat') 
    fig=figure;
    set(fig, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 0.6, 0.95]);
    pcolor(phi*180/pi,90-th*180/pi,real(EnergyV)); 
    shading interp; 
    set(0,'defaulttextInterpreter','latex') 
    shading interp;
    if strcmpi(type,'total')
        colormap(cmocean('solar',100))
    else
        colormap(cmocean('balance',100))
    end
    set(0,'defaulttextInterpreter','none')
    if no_grid==0
        m_grid('xaxislocation','middle','xtick',[-90:90:180],'ytick',[-90:90:90],'fontsize',20,'gridcolor',[0.8 0.8 0.8],'linewidth',3,'linestyle','--','gridcolor','k');
    else
        m_grid('xaxislocation','middle','xtick',[-90:90:180],'ytick',[-90:90:90],'gridcolor',[0.8 0.8 0.8],'linewidth',3,'linestyle','--','gridcolor','k','fontsize',1);
    end
    if no_colorbar==0
    cb=colorbar('southoutside');
    set(0,'defaulttextInterpreter','latex')
    cb.Label.Interpreter = 'latex';
    if isempty(label_z)==1
        if strcmpi(type,'total')
            cb.Label.String='$\dot{e}$ [-]';
        else
            cb.Label.String='$\Delta\dot{e}$ [-]';
        end
    else
        cb.Label.String=label_z; 
        cb.FontSize=50;
    end
    cb.FontName = 'CMU Serif';
    end
    set(gca,'fontsize', 50);
    set(gcf,'color','w');
     if isempty(limits)==0
        caxis([limits(1) limits(2)])
     else
         if strcmpi(type,'total')
            caxis([0 max_E])
        else
           caxis([-max_E2 max_E2])
        end
     end
     if isempty(title_plot)==0
        title(title_plot,'interpreter','latex')
    end
elseif strcmpi(projection,'cut') 
    fig=figure;
    if strcmpi(type,'units')
        set(fig, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 1]); 
    else
        set(fig, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 0.75, 1]);   
    end
    % MAP 
    if no_colorbar==1
        pos1 = [0.1 0.15 0.7 0.53];
    else
        if strcmpi(type,'units')
            pos1 = [0.1 0.13 0.8 0.85];
        else
            pos1 = [0.1 0.08 0.7 0.6];
        end
    end
    pos2 = [0.8 0.15 0.15 0.53];
    pos3=[0.1 0.68 0.7 0.2];
    subplot('Position',pos1)
    set(gca,'fontsize', 20);
    pcolor(phi*180/pi,90-th*180/pi,real(EnergyV)); 
    shading interp;  
    if strcmpi(type,'total')
        colormap(cmocean('solar',100))
    elseif strcmpi(type,'units')
        colormap(cmocean('thermal',50))
    else
        colormap(cmocean('balance',100))
    end
    set(0,'defaulttextInterpreter','none')
    if no_colorbar==0
    cb=colorbar('southoutside');
    set(0,'defaulttextInterpreter','latex')
    if isempty(label_z)==1
        cb.Label.String='$\Delta\dot{e}$' ;
    else
        cb.Label.String=label_z; 
        cb.FontSize=30;
    end
    cb.Label.Interpreter = 'latex';
    cb.FontName = 'CMU Serif';
    end
    set(gcf,'color','w');
    if isempty(limits)==0
        caxis([limits(1) limits(2)])
     else
         if strcmpi(type,'total')
            caxis([0 max_E])
        else
           caxis([-max_E2 max_E2])
        end
     end    
    xticks(-180:45:180)
    xlim([-180 180])
    yticks(-90:45:89)
    ylim([-90 90])
    if ~strcmpi(type,'units')
        xticklabels('')
    end
    set(gca,'fontsize', 30);
    set(gca,'TickLabelInterpreter','latex')
    ylabel('Latitude [deg]','interpreter','latex' )
    grid minor
    set(gca, 'layer', 'top');
    text(20,-75,title_plot,'Interpreter','latex','FontSize',40,'color',[0.5 0.5 0.5])
    if ~strcmpi(type,'units')
        % Latitude averaged 
        subplot('Position',pos2)
        box on 
        plot(Energy_Lat/(2*pi),90-th_cut*180/pi,'LineWidth',3,'color','k');
        set(gca,'TickLabelInterpreter','latex')
        yticks(-90:45:90)
        ylim([-90 90])
        ax = gca;
        ax.YAxis.FontSize =30;
        ax.XAxis.FontSize =30;
        xlabel('$\frac{1}{2\pi\dot{e}_0^u}\int\Delta\dot{e}d\phi$','interpreter','latex','FontSize',25)
        yticklabels('')
        grid minor
        %Longitude averaged
        subplot('Position',pos3)
        box on
        plot(phi_cut*180/pi,Energy_Lon/2,'LineWidth',3,'color','k');
        xticks(-180:45:180)
        xlim([-180 180]) 
        ax = gca;
        ax.YAxis.FontSize =30;
        ax.XAxis.FontSize =30;
        set(gca, 'XAxisLocation', 'top')
        xlabel('Longitude [deg]','interpreter','latex','FontSize',30)
        yyaxis right
        ylabel('$\frac{1}{2\dot{e}_0^u}\int\Delta\dot{e}\sin\theta d\theta$','color','k','interpreter','latex','FontSize',25)
        set(gca,'YTick',[])
        set(gca,'ycolor','k') 
        grid minor
        set(gca,'TickLabelInterpreter','latex')
        set(gcf,'color','w');
    else
        xticks(-180:45:180)
        xlabel('Longitude [deg]','interpreter','latex','FontSize',30)
    end
else
    error('Incorrect entry for projection')
end
set(gcf,'color','w');
%% Save Plot
if (save_plot)==1
    %export_fig(fig,save_name,'-pdf','-opengl')
    saveas(fig,save_name)%,'-png','-opengl','-r400')
end
    
end