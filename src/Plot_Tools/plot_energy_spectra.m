function [varargout] = plot_energy_spectra(n_var,m_var,n_v,m_v,energy_s,type,save_name,varargin)
%% Convert into phase and amplitude
energy_spectra_auxR=zeros(size(energy_s));
energy_spectra_amplitude=[];
energy_spectra_phase=[];
spectra_l=size(energy_s,1);
multiple=0; 

for k = 1:length(varargin)
    if strcmpi(varargin{k},'multiple')
        multiple=varargin{k+1}; 
        varargin{k+1}=[]; 
        varargin{k}=[];
    end
end

for i=1:size(energy_s,1)
    k=1; 
    for n_aux=0:max(n_v) %loop over degrees
        for m_aux=0:n_aux %loop over oders
            if m_aux==0
                index=find(m_aux==m_v & n_aux==n_v);
                energy_spectra_auxR(i,index)=energy_s(i,index);
                energy_spectra_amplitude(i,k)=abs(energy_spectra_auxR(i,index));
                energy_spectra_phase(i,k)=0; 
            else
                indexP=find(m_aux==m_v & n_aux==n_v);
                indexN=find(-m_aux==m_v & n_aux==n_v);
                %energy_spectra_auxR(i,indexP)=1/sqrt(2)*((-1)^m_aux*energy_s(i,indexP)+energy_s(i,indexN));
                %energy_spectra_auxR(i,indexN)=1/sqrt(2)*(1i*(-1)^m_aux*energy_s(i,indexP)-1i*energy_s(i,indexN));                
                energy_spectra_auxR(i,indexP)=sqrt(2)*real(energy_s(i,indexP));
                energy_spectra_auxR(i,indexN)=-sqrt(2)*imag(energy_s(i,indexP));                
                energy_spectra_amplitude(i,k)=sqrt(energy_spectra_auxR(i,indexP)^2+energy_spectra_auxR(i,indexN)^2);
                energy_spectra_phase(i,k)=atan2d(real(energy_spectra_auxR(i,indexN)),real(energy_spectra_auxR(i,indexP)));
            end
            n_v2(k)=n_aux;
            m_v2(k)=m_aux;
            k=k+1;
        end
    end
end

%% Make plot
energy_uniform_plot=log10(energy_spectra_amplitude);
energy_uniform_plot(end+1,:)=0;
energy_uniform_plot(:,end+1)=0;
aux=abs(energy_uniform_plot);
if multiple==0
    fig=figure;
    set(fig, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.45]);
else
    fig=figure;
    set(fig, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 0.6, 1]);
end
set(fig,'defaultAxesColorOrder',[0 0 0; 0 0 0]);
h=pcolor(real((energy_uniform_plot)))
set(h, 'EdgeColor', 'white');
colormap(cmocean('thermal'))
cb=colorbar
if strcmpi(type,'total') 
    cb.Label.String='$\sqrt{{\dot\mathcal{E}_n^m}^2+{\dot\mathcal{E}_n^{-m}}^2}$';
else
    cb.Label.String='$\frac{1}{\dot{e}_0^u}\sqrt{{\Delta\dot\mathcal{E}_n^m}^2+{\Delta\dot\mathcal{E}_n^{-m}}^2}$';
end
cb.Label.Interpreter = 'latex';
cb.Location='northoutside';
cb.Ticks =[-10:1:0]; 
cb.Label.Interpreter = 'latex';
cb.TickLabels ={'10^{-10}','','10^{-8}','','10^{-6}','','10^{-4}','','10^{-2}','','1',}; 
cb.Label.Interpreter = 'latex';
cb.FontName = 'CMU Serif';
caxis([-10 0]);
i_L=[];
j_l=[];
k=1;
kk=1; 
for i=1:length(n_v2)
        if m_v2(i)==0
           j_l=[j_l k-0.5];
           if n_v2(i)==0
               labels_P{kk}=[''];
           else
                labels_P{kk}=['(' num2str(n_v2(i)) ',' num2str(m_v2(i)) ')'];
           end
          kk=kk+1; 
        elseif m_v2(i)==n_v2(i)
          i_L=[i_L k];
        end
        k=k+1;
end
hold on
plot([2 2],[0 2+spectra_l],'LineWidth',2,'color','w');
hold on
plot([1 1],[0 2+spectra_l],'LineWidth',2,'color','w');
for i=1:length(i_L)
    plot([i_L(i)+1 i_L(i)+1],[0 1+spectra_l],'LineWidth',2,'color','w');
    hold on
end

y_l=[];
for i=1:length(n_var)
    y_l=[y_l i+0.5];
    label_y{i}=['(' num2str(n_var(i)) ',' num2str(m_var(i)) ')'];
    plot([0 length(n_v2)+1],[i i],'LineWidth',2,'color','w');
    hold on 
end

xticks(j_l+1)
xticklabels(labels_P)
yticks(y_l)
yticklabels(label_y)
axis equal
% set-labels y-axis
label_variations_location=[];
label_variations=[];
label_vis_location=[];
ylim([1 1+spectra_l])
if multiple==4
     plot([0 length(n_v2)+1],[1 1],'LineWidth',6,'color','w');
     hold on
     plot([0 length(n_v2)+1],[26 26],'LineWidth',6,'color','w');
     hold on 
    for j=1:4
        plot([0 length(n_v2)+1],[1+5*j 1+5*j],'LineWidth',6,'color','w');
        hold on
        yticks_loc(j)=3.5+5*(j-1);
    end
    labels_y2={'$Y_2^0\textrm{e}^{\textrm{i}\omega t}+c.c$','$\frac{\sqrt{2}}{4}(Y_2^2+Y_2^{-2})\textrm{e}^{\textrm{i}\omega t}+c.c$',...
            '$Y_2^2\textrm{e}^{\textrm{i}\omega t}+c.c$','ecc. synch. moon'};
    yl = ylim;
    yyaxis right
    ylim(yl)
    yticks(yticks_loc)
    yticklabels(labels_y2)
end

set(gca,'fontsize', 30);
set(gcf,'color','w');
set(gca,'TickLabelInterpreter','latex')
ax = gca;
ax.YAxis.FontSize =16;
ax.XAxis.FontSize =16;

if isempty(save_name)==0
        export_fig(fig,[save_name '_s'],'-pdf','-opengl')
end 

%% PLOT THE PHASE 
if strcmpi(type,'total') 
else
energy_spectra_phase(end+1,:)=0;
energy_spectra_phase(:,end+1)=0;
if multiple==0
    fig=figure;
    set(fig, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.45]);
else
    fig=figure;
    set(fig, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 0.6, 1]);
end
set(fig,'defaultAxesColorOrder',[0 0 0; 0 0 0]);
h=pcolor(energy_spectra_phase)
set(h, 'EdgeColor', 'white');
colormap_p_aux=cmocean('phase',1001);
colormap_p=[1 1 1; colormap_p_aux(1:500,:); [0 0 0]; colormap_p_aux(502:end,:) ;1 1 1];
colormap(colormap_p)
cb=colorbar;
cb.Label.String='$\varphi_n^m$ [deg]';
cb.Label.Interpreter = 'latex';
cb.Location='northoutside';
cb.Label.Interpreter = 'latex';
cb.FontName = 'CMU Serif';
caxis([-180 180]);

i_L=[];
j_l=[];
k=1;
kk=1; 
for i=1:length(n_v2)
        if m_v2(i)==0
           j_l=[j_l k-0.5];
           if n_v2(i)==0
               labels_P{kk}=[''];
           else
                labels_P{kk}=['(' num2str(n_v2(i)) ',' num2str(m_v2(i)) ')'];
           end
          kk=kk+1; 
        elseif m_v2(i)==n_v2(i)
          i_L=[i_L k];
        end
        k=k+1;
end
hold on
plot([2 2],[0 2+spectra_l],'LineWidth',2,'color','w');
hold on
plot([1 1],[0 2+spectra_l],'LineWidth',2,'color','w');
for i=1:length(i_L)
    plot([i_L(i)+1 i_L(i)+1],[0 1+spectra_l],'LineWidth',2,'color','w');
    hold on
end

y_l=[];
for i=1:length(n_var)
    y_l=[y_l i+0.5];
    label_y{i}=['(' num2str(n_var(i)) ',' num2str(m_var(i)) ')'];
    plot([0 length(n_v2)+1],[i i],'LineWidth',2,'color','w');
    hold on 
end
xticks(j_l+1)
xticklabels(labels_P)
yticks(y_l)
yticklabels(label_y)
axis equal
% set-labels y-axis
label_variations_location=[];
label_variations=[];
label_vis_location=[];
ylim([1 1+spectra_l]) 
if multiple==4
     plot([0 length(n_v2)+1],[1 1],'LineWidth',6,'color','w');
     hold on
     plot([0 length(n_v2)+1],[26 26],'LineWidth',6,'color','w');
     hold on 
    for j=1:4
        plot([0 length(n_v2)+1],[1+5*j 1+5*j],'LineWidth',6,'color','w');
        hold on
        yticks_loc(j)=3.5+5*(j-1);
    end
    labels_y2={'$Y_2^0\textrm{e}^{\textrm{i}\omega t}+c.c$','$\frac{\sqrt{2}}{4}(Y_2^2+Y_2^{-2})\textrm{e}^{\textrm{i}\omega t}+c.c$',...
            '$Y_2^2\textrm{e}^{\textrm{i}\omega t}+c.c$','ecc. synch. moon'};
    yl = ylim;
    yyaxis right
    ylim(yl)
    yticks(yticks_loc)
    yticklabels(labels_y2)
end

set(gca,'fontsize', 30);
set(gcf,'color','w');
set(gca,'TickLabelInterpreter','latex')
ax = gca;
ax.YAxis.FontSize =16;
ax.XAxis.FontSize =16;

if isempty(save_name)==0
        export_fig(fig,[save_name '_phase'],'-pdf','-opengl')
end 


end

varargout{2}=energy_spectra_phase; 
varargout{1}=energy_spectra_amplitude;

end