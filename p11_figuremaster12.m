% this script is used to create figures % figure master 8 and loaddata11
% for 1st submission

% there are 11 datasets in total


%%%%% LOAD in the all gridded datase %%%%%%%%%%%%%%%%

source_names = {'EM1','EM2','EM3','EM4','EM5','EO','LUNA','RF'};

vData.Vall = nan(360,720,11);

vData.Vall(:,:,1) = vData.Vn1; %EM1
vData.Vall(:,:,2) = vData.Vn2; %EM2
vData.Vall(:,:,3) = vData.Vn3; %EM3
vData.Vall(:,:,4) = vData.Vn4; %EM4
vData.Vall(:,:,5) = vData.Venv1; %EM5
vData.Vall(:,:,6) = vData.Venv2; %EO
vData.Vall(:,:,7) = vData.Vne; %LUNA
vData.Vall(:,:,8) = vData.annVc; %RF


% remove non-vegetated land surface
veg_land = LC05D_renum > 0;

for i = 1:8
    tmp = squeeze(vData.Vall(:,:,i));
    tmp(~veg_land) = NaN;
    vData.Vall(:,:,i) = tmp;
end

% fLNR from all Vcmax
vData.Fall(:,:,:) = vData.Vall./vData.Na./47.34./6.5 * 100;



%%%%%%%%% Figure 1, Varations of Vcmax25, ability of models, Vcmax25 and fLNR map
f1=figure('Name','histo of estimates','Units', 'centimeters','Color','white', 'Position', [2, 2, 15, 12], ...
    'OuterPosition', [2, 2, 15, 12]);

lat=-89.75:0.5:89.75;
lon=-179.75:0.5:179.75;

% sub_panel figure
SpacingVert = 0.14;
SpacingHoriz = 0.1;
MR = 0.07;
ML = 0.14;
MarginTop = 0.03;
MarginBottom = 0.12;

%%%%%% Fig.1a. variations of Vcmax25
cd('/Volumes/RemiLBNL/project6/code4Remi/functions'); 
% use a third party function https://www.mathworks.com/matlabcentral/fileexchange/3696-subaxis-subplot
subaxis(2,3,1,'SpacingVert',SpacingVert,'SpacingHoriz',SpacingHoriz,'MR',MR,'ML',ML,'MarginTop',MarginTop,'MarginBottom',MarginBottom); 

tmp_v = [];

for i = 1: size(vData.clnALL,1)
    
    lat_n = round((89.75 - vData.clnALL(i,1))./0.5); 
    lon_n = round((vData.clnALL(i,2) + 179.75)./0.5); 
    
    tmp_v(i,1) = squeeze(vData.Vall(lat_n,lon_n,8));
    
end

    ydata_f = tmp_v;
    xdata_f = vData.clnALL(:,4);
    
    xmin = nanmin(xdata_f);
    xmax = nanmax(xdata_f);
    
    ymin = nanmin(ydata_f);
    ymax = nanmax(ydata_f);

    xdata_all=[ones(size(ydata_f)) xdata_f];
    
    % type II regression, use a thrid-party function gmregress
    cd('/Users/xzluo/Dropbox (Personal)/Code/common_functions');
    [b,bint,~] = gmregress(xdata_f,ydata_f, 0.05); 

    mdl = fitlm(xdata_f,ydata_f);
    [ypred,yci] = predict(mdl);
    ydata_low=yci(:,1);
    ydata_up=yci(:,2);


    % scatter plot of obs. and est. gray gradient to show sample densities
    fg_z = [];
    inter_n = 100;
    fg_x=xmin:((xmax-xmin)/inter_n):xmax;
    fg_y=ymin:((ymax-ymin)/inter_n):ymax;
    for c_x=1:inter_n
        for c_y=1:inter_n
            c_ind=xdata_f>fg_x(c_x) & xdata_f<fg_x(c_x+1)...
                & ydata_f>fg_y(c_y) & ydata_f<fg_y(c_y+1);

            fg_z(c_y,c_x)=sum(c_ind);
        end
    end

    fg_z(fg_z==0)=NaN;

    h=pcolor(fg_x(1:100),fg_y(1:100),fg_z);%if interval too large, then the color cannot show
    set(h, 'EdgeColor', 'none');

    color_s=flip(gray(10));
    colormap(gca,color_s(3:10,:)); 
    
    hold on;
    % fitted line
    pxx=plot(xdata_f, xdata_f.*b(2)+b(1),'--','color','r');

    [r1,r2]=corrcoef(xdata_f,ydata_f,'rows','complete');
    sig=round(r2(2)*100)./100; 

    r = round(r1(1,2)^2.*100)./100;
    b = round(b.*100)./100;

    if sig < 0.01
        text(xmin+0.02*(xmax-xmin),18,strcat('r^2=',num2str(r),' (p<0.01) '),'Color','r','FontSize',8);
    else
        text(xmin+0.02*(xmax-xmin),18,strcat('r^2=',num2str(r),'(p=',num2str(sig),')'),'Color','r','FontSize',8);
    end

    text(xmin+0.02*(xmax-xmin),5,strcat('y =',num2str(b(2)),'x+',num2str(b(1))),'Color','r','FontSize',8);

    ylim([0 120]);
    xlim([0 120]);

    ylabel({'RF Vc_{max}^2^5 (\mumol m^-^2s^-^1)'},'FontSize',8);
    xlabel({'Obs. Vc_{max}^2^5 (\mumol m^-^2s^-^1)'},'FontSize',8);


    set(gca, 'box', 'off');


%%%%%% Fig. 1c Explanary ability of every models
cd('/Volumes/RemiLBNL/project6/code4Remi/functions');
subaxis(2,3,4,'SpacingVert',SpacingVert,'SpacingHoriz',SpacingHoriz,'MR',MR,'ML',ML,'MarginTop',MarginTop,'MarginBottom',MarginBottom); 

bar_res = [];
tmp_v = [];

for j = 1:7
    
    bar_res(j) = corr(reshape(vData.Vall(:,:,j),[],1),reshape(vData.Vall(:,:,8),[],1),'rows','complete');
    
end


b = bar(bar_res,'FaceColor',[166 206 227]/255,'EdgeColor','none');

set(gca, 'box', 'off');

ylim([-0.2 0.5]);

x_names = source_names(1:7);

ylabel({'Correlation (r) between';'models and RF Vc_{max}^2^5'},'FontSize',8);
h = gca;
h.XTick = 1:length(x_names);
h.XTickLabel = x_names;
h.XTickLabelRotation = 60;
h.TickLabelInterpreter = 'none';
set(gca,'XTickLabel',x_names,'Fontsize',8)


%%%%% Fig. 1d map of fLNR
cd('/Volumes/RemiLBNL/project6/code4Remi/functions');
subaxis(2,3,[5,6],'SpacingVert',0,'SpacingHoriz',SpacingHoriz,'MR',0.05,'ML',0.01,'MarginTop',MarginTop,'MarginBottom',0); 

vcmax2 = squeeze(vData.Fall(:,:,8));    
mapdata=flip(vcmax2);
    
cd('/Volumes/RemiLBNL/project6/code4Remi/functions');
mapdata= smooth2a(mapdata,3,3); 
mapdata(180,359)=-50;
mapdata(180,360)=50;

cd('/Volumes/RemiLBNL/project6/code4Remi/functions/m_map');
m_proj('Robinson','lon',[-180 180],'lat',[-60 90]); % map projections, range of lon and lat of your data
m_coast('linewidth',0.5,'color',[0.2 0.2 0.2]); % coast line settings
hold on;
m_pcolor(lon,lat,mapdata); % draw your map here
shading INTERP;   % can be flat or INTERP

colorscheme = [84,48,5;...
            140,81,10;...
            191,129,45;...
            223,194,125;...
            246,232,195;...
            199,234,229;...
            128,205,193;...
            53,151,143;...
            1,102,94;...
            0,60,48];
colorscheme = colorscheme./255;
% colorscheme = jet(12);       

colormap(gca,colorscheme);

caxis([0 50]); 
    
h=colorbar;%
set(h,'YTick',0:10:50);
ylabel(h,{'fLNR (%)'},'FontSize',8);
cd('/Volumes/RemiLBNL/project6/code4Remi/functions');cbfreeze(h);

cd('/Volumes/RemiLBNL/project6/code4Remi/functions/m_map');
m_grid('box','off','tickdir','in','yticklabels',[],'xticklabels',[],'fontsize',7);


hold off;
    


%%%%% Fig. 1b map of Vcmax
cd('/Volumes/RemiLBNL/project6/code4Remi/functions');
subaxis(2,3,[2,3],'SpacingVert',0,'SpacingHoriz',SpacingHoriz,'MR',0.05,'ML',0.01,'MarginTop',MarginTop,'MarginBottom',0); 

vcmax2 = squeeze(vData.Vall(:,:,8));
mapdata=flip(vcmax2);
    
cd('/Volumes/RemiLBNL/project6/code4Remi/functions');
mapdata= smooth2a(mapdata,3,3); 
mapdata(180,359)=-50;
mapdata(180,360)=50;

cd('/Volumes/RemiLBNL/project6/code4Remi/functions/m_map');
m_proj('Robinson','lon',[-180 180],'lat',[-60 90]); % map projections, range of lon and lat of your data
m_coast('linewidth',0.5,'color',[0.2 0.2 0.2]); % coast line settings
hold on;
m_pcolor(lon,lat,mapdata); % draw your map here
shading INTERP;   % can be flat or INTERP

colormap(gca,colorscheme);

caxis([20 120]); 
    
h=colorbar;%
set(h,'YTick',20:20:120);
ylabel(h,{'Vc_{max}^2^5 (\mumol m^-^2s^-^1)'},'FontSize',8);
cd('/Volumes/RemiLBNL/project6/code4Remi/functions');cbfreeze(h);

cd('/Volumes/RemiLBNL/project6/code4Remi/functions/m_map');
m_grid('box','off','tickdir','in','yticklabels',[],'xticklabels',[],'fontsize',7);

hold off;
    
saveas(gcf,'/Volumes/RemiLBNL/project11_Vcmax_val/Figures/Figure_1.pdf');

set(f1,'PaperPositionMode','auto');
print('/Volumes/RemiLBNL/project11_Vcmax_val/Figures/Figure_1','-djpeg','-r600');

%%%%%%%%% End: Figure 1 %%%%%%%%%%%%%%%








%%%%%%%%%%%% Figure 2, Response of fLNR to PFT, leaf traits, climate and soil %%%%%%%%%%%%%%%%%%

f2=figure('Name','histo of estimates','Units', 'centimeters','Color','white', 'Position', [2, 2, 15, 15], ...
    'OuterPosition', [2, 2, 15, 15]);


% PFT fLNR
cd('/Volumes/RemiLBNL/project6/code4Remi/functions');
subaxis(2,3,1:3,'SpacingVert',0.1,'SpacingHoriz',0.12,'MR',0.05,'ML',0.1,'MarginTop',.02,'MarginBottom',.05); 

% organise the data
bar_pft=nan(10000,8,8);
 
 for i=1:8 % for every PFT
     
     for mm = 1:8 % for every models
         
         if mm == 1 % plot the fLNR first
             tmp = reshape(vData.fLNR,[],1);
             ind = reshape(EULC_2D,[],1) == i;
             tmp2 = tmp(ind);
             
             bar_pft(1:length(tmp2),i,mm)=tmp2;
         
         else
             tmp = reshape(vData.Fall(:,:,mm-1),[],1);
             ind = reshape(EULC_2D,[],1) == i;
             tmp2 = tmp(ind);
             
             bar_pft(1:length(tmp2),i,mm)=tmp2;
             
         end
         
     end
 
 end
 
bar_pft(bar_pft<1) = NaN;
 
% color gradient
color_s = {[0,0,0],...
        [178,24,43]/255,...
        [214,96,77]/255,...
        [244,165,130]/255,...
        [253,219,199]/255,...
        [140,81,10]/255,... 
        [67,147,195]/255,...
        [84,39,136]/255};

ld_text={'RF','EM1','EM2','EM3','EM4','EM5','EO','LUNA'};
 
% use a third-party function iosr.statistics
cd('/Users/xzluo/Dropbox (Personal)/Code/common_functions/IoSR-Surrey-MatlabToolbox-4bff1bb/');
bp=iosr.statistics.boxPlot(bar_pft,'theme','colorlines','lineColor',color_s...
,'lineWidth',0.6,'medianColor',color_s,'method','R-5','meanColor',color_s,'boxWidth',0.05,'limit',[10,90]);
bp.showMean=true;
bp.showOutliers=false;
bp.meanSize=4;
bp.xSeparator =true;
bp.symbolColor =[0 0 0];
bp.outlierSize =0.5;

bp.boxWidth = 0.1;

bp.showScatter = false;
bp.scatterSize = 3;
bp.scatterColor = [1 0.4 0.4];
bp.scatterMarker = '.';
bp.scatterAlpha = 0.6;
bp.showLegend = true; % had to add a squeeze function in boxPlot function L1974
bp.groupLabels = ld_text;
bp.handles.legend.Location = 'eastoutside';


ylabel({'fLNR(%)'}, 'FontSize',9.5);
ylim([0 50]);
text(0.5,49,'a','FontWeight','bold');

ld_text={'CRO','DBF','EBF','ENF','MF','GRA','SH','WET'};
set(gca,'XTickLabel',ld_text,'FontSize',10,'box', 'off','FontSize',9.5);


%%%% partial response of fLNR to leaf traits

filename = '/Volumes/RemiLBNL/project11_Vcmax_val/data/fann_pca_gam_dec.csv';
fann_gam = readtable(filename);

filename = '/Volumes/RemiLBNL/project11_Vcmax_val/data/fenv1_pca_gam_dec.csv';
fenv1_gam = readtable(filename);

filename = '/Volumes/RemiLBNL/project11_Vcmax_val/data/fenv2_pca_gam_dec.csv';
fenv2_gam = readtable(filename);

filename = '/Volumes/RemiLBNL/project11_Vcmax_val/data/fne_pca_gam_dec.csv';
fne_gam = readtable(filename);

sz = 10;

color_s = [0,0,0;...
140,81,10;...
33,102,172;...
84,39,136]/255;

for i = 1 : 3 % for leaf traits, climate and soil
    
    % each panel
    cd('/Volumes/RemiLBNL/project6/code4Remi/functions');
    subaxis(2,3,i+3,'SpacingVert',0.1,'SpacingHoriz',0.05,'MR',0.08,'ML',0.1,'MarginTop',.1,'MarginBottom',.15); 
    
    % define ydata
    switch i 
        case 1 % leaf traits
            ind_y = [2,12];
            x_title = 'leaf traits PC1';
            
        case 2 % climate
            ind_y = [5,15];
            x_title = 'climate PC1';
            
        case 3 % soil
            ind_y = [8,18];
            x_title = 'soil PC1';
    end
    
    
    % define xdata
    for j = 1 : 4 % for RF, EM5, EO and LUNA.
        
        switch j
            case 1 % RF
                tmp_xdata = horzcat(gam_pcas,reshape(EULC_2D,[],1),reshape(vData.annVc,[],1));%
                xdata = tmp_xdata(:,(i-1)*3 + 1);
                
                ydata = table2array(fann_gam(:,ind_y(1)));
                ydata_std = table2array(fann_gam(:,ind_y(2)));
                
            case 2 % EM5
                tmp_xdata = horzcat(gam_pcas,reshape(EULC_2D,[],1),reshape(vData.Venv1,[],1));%
                xdata = tmp_xdata(:,(i-1)*3 + 1);
                
                ydata = table2array(fenv1_gam(:,ind_y(1)));
                ydata_std = table2array(fenv1_gam(:,ind_y(2)));
                   
            case 3 % EO
                tmp_xdata = horzcat(gam_pcas,reshape(EULC_2D,[],1),reshape(vData.Venv2,[],1));%
                xdata = tmp_xdata(:,(i-1)*3 + 1);
                
                ydata = table2array(fenv2_gam(:,ind_y(1)));
                ydata_std = table2array(fenv2_gam(:,ind_y(2)));
                
            case 4 % LUNA
                tmp_xdata = horzcat(gam_pcas,reshape(EULC_2D,[],1),reshape(vData.Vne,[],1));%
                xdata = tmp_xdata(:,(i-1)*3 + 1);
                
                ydata = table2array(fne_gam(:,ind_y(1)));
                ydata_std = table2array(fne_gam(:,ind_y(2)));
        end
        
        xdata(any(isnan(tmp_xdata), 2), :) = [];

        
        [xdata2,sort_ind] = sort(xdata);
        hold on;
        
        cd('/Volumes/RemiLBNL/project6/code4Remi/functions');
        
        H=shadedErrorBar(xdata2,ydata(sort_ind),sqrt(100)*ydata_std(sort_ind),{'-'},0); 

        set(H.patch,'FaceColor',color_s(j,:),'EdgeColor',color_s(j,:),'FaceAlpha',0.2,'EdgeAlpha',0.2)
        set(H.mainLine,'color',color_s(j,:));
        set(H.edge,'color','none');
        
        pp(j) = plot(xdata2,ydata(sort_ind),'-','LineWidth',1.5,'Color',color_s(j,:));


        hold on;
        
 
    end
    
    if i == 1
        ylabel('\DeltafLNR (%)','FontSize',9);
    else
        set(gca,'yticklabels',[]);
    end
    
    xlabel(x_title,'FontSize',9);
    ylim([-20 20]);
    
    text(min(xdata)+0.05*(max(xdata)-min(xdata)),18, char(97+i),'FontWeight','bold');

    
end

l2=legend(pp,{'RF','EM5','EO','LUNA'},'FontSize',8,'box','off','Orientation','horizontal','Location',[0.25,-0.01,0.5,0.1]);  



saveas(gcf,'/Volumes/RemiLBNL/project11_Vcmax_val/Figures/Figure_2.pdf');

set(f2,'PaperPositionMode','auto');
print('/Volumes/RemiLBNL/project11_Vcmax_val/Figures/Figure_2','-djpeg','-r600');  



%%%%%%%%%%%% End: Figure 2, Response of fLNR to PFT, leaf traits, climate and soil %%%%%%%%%%%%%%%%%%






%%%%%%%%%%%% Figure 3, Response of fLNR to leaf traits%%%%%%%%%%%%%%%%%%

f3=figure('Name','fLNRvsClimate','Units', 'centimeters','Color','white', 'Position', [2, 2, 18, 20], ...
    'OuterPosition', [2, 2, 18, 20]);

% sub_panel figure
SpacingVert = 0.01;
SpacingHoriz = 0.01;
MR = 0.05;
ML = 0.05;
MarginTop = 0.05;
MarginBottom = 0.05;


%%%%% the first part, dominate factor

filename = '/Volumes/RemiLBNL/project11_Vcmax_val/data/fann_pca_gam_dec.csv';
fann_gam = readtable(filename);

tmp1 = abs(fann_gam.s_leaf1_ + fann_gam.s_leaf2_);
tmp2 = abs(fann_gam.s_clim1_ + fann_gam.s_clim2_);
tmp3 = abs(fann_gam.s_soil1_ + fann_gam.s_soil2_);

tmp_xdata = horzcat(gam_pcas,reshape(EULC_2D,[],1),reshape(vData.annVc,[],1));%

ind = any(isnan(tmp_xdata), 2);

contri = nan(360*720,3);

contri(~ind,2) = tmp1./(tmp1 + tmp2 + tmp3).*1; % leaf traits
contri(~ind,1) = tmp2./(tmp1 + tmp2 + tmp3).*1.5; % climate
contri(~ind,3) = tmp3./(tmp1 + tmp2 + tmp3).*1.5; % soil

%%% smooth the contribution map over space, cannot smooth later as the ind has no physical meaning

for i = 1:3
    
    contri_2D = reshape(contri(:,i),360,720);
    cd('/Volumes/RemiLBNL/project6/code4Remi/functions');
    contri_2D= smooth2a(contri_2D,10,10); 
    contri_2D(180,359)=-50;
    contri_2D(180,360)=50;
    
    contri(:,i) = reshape(contri_2D,[],1);
    
end


%%%% triangle map
cd('/Volumes/RemiLBNL/project6/code4Remi/functions');
subaxis(3,9,[1:6,10:15],'SpacingVert',SpacingVert,'SpacingHoriz',SpacingHoriz + 0.03,'MR',MR,'ML',ML,'MarginTop',MarginTop,'MarginBottom',MarginBottom); 


color_scheme = [];

A = [0:0.1:1; 0:0.1:1; 0:0.1:1;]';

[ca, cb, cc] = ndgrid(A);
comb = [ca(:), cb(:), cc(:)];
comb2 = unique(comb, 'rows');

% one decimal places
contri = round(contri.*10)./10;


comb3 = comb2(2:end,:);

for i = 1:length(contri)
    
    if isnan(contri(i,1))
        Locb(i) = NaN;
    else
        [Lia, Locb_tmp] = ismembertol(contri(i,:),comb2,'ByRows',true);
        Locb(i) = Locb_tmp;
    end
end

Locb(ind) = NaN;


vcmax2 = reshape(Locb,360,720);    
mapdata=flip(vcmax2);
    
cd('/Volumes/RemiLBNL/project6/code4Remi/functions');
mapdata(180,359)=-50;
mapdata(180,360)=50;

cd('/Volumes/RemiLBNL/project6/code4Remi/functions/m_map');
m_proj('Robinson','lon',[-180 180],'lat',[-60 90]); % map projections, range of lon and lat of your data
m_coast('linewidth',0.5,'color',[0.2 0.2 0.2]); % coast line settings
hold on;
m_pcolor(lon,lat,mapdata); % draw your map here
m_grid('box','off','tickdir','in','fontsize',7);

shading INTERP;   % can be flat or INTERP
  
colorscheme = comb2;
colormap(gca,colorscheme);
caxis([1,1331]);


cd('/Volumes/RemiLBNL/project6/code4Remi/functions/m_map');
m_grid('box','off','tickdir','in','fontsize',7);



%%%%%% rgb map to show the contribution from each factor %%%%%%%%% 
axis;

xy=[0 0; 1 0; 0.5 sqrt(3)/2] - [2.7,0.7]; % add the second compoenent to move the position of the triangle
col=[1 0 0; 0 1 0; 0 0 1];
pp = patch('Vertices',xy, 'Faces',[1:size(xy,1)], 'EdgeColor','none','FaceVertexCData', col,'FaceColor','interp');
axis off
text(xy(1,1), xy(1,2)-0.05,'climate','HorizontalAlignment', 'center','FontSize',9);
text(xy(2,1), xy(2,2)-0.05,'leaf traits','HorizontalAlignment', 'center','FontSize',9);
text(xy(3,1), xy(3,2)+0.05,'soil','HorizontalAlignment', 'center','FontSize',9);


hold off;



%%%%%% the second part, key variables for the global Vcmax

% fit the multi-variate linear regression

ydata_raw = nan(size(reshape(vData.fLNR,[],1)));


%%%% get the index for PCAs that have values
tmp_xdata = horzcat(gam_pcas,reshape(EULC_2D,[],1),reshape(vData.annVc,[],1));%
nan_ind = any(isnan(tmp_xdata),2);
ydata_pft = nan(size(ydata_raw));
ydata_pft(~nan_ind) = fann_gam.PFT;


contri_all = [];
contri_all_std = [];
emp_coef = [];
emp_const = [];

for j = 1:3 % for each type of factor: leaf traits, climate, nutrient
    
    
    switch j
        case 1
            
            xdata_raw = horzcat(reshape(vData.Pa,[],1),reshape(vData.LMA,[],1));
            ydata_raw(~nan_ind) = fann_gam.s_leaf1_ + fann_gam.s_leaf2_ + fann_gam.s_leaf3_;

        case 2
           
            xdata_raw = grid_predictors(:,6:10);
            ydata_raw(~nan_ind) = fann_gam.s_clim1_+ fann_gam.s_clim2_ + fann_gam.s_clim3_;

            
        case 3
         
            xdata_raw = grid_predictors(:,[13:20]);
            ydata_raw(~nan_ind) = fann_gam.s_soil1_+ fann_gam.s_soil2_ + fann_gam.s_soil3_;
    
    end


    ydata_raw = ydata_raw; 
    xdata =xdata_raw(~nan_ind,:);
    ydata =ydata_raw(~nan_ind);

    xdata(xdata == -Inf) = NaN;

    mdl = fitlm(xdata,ydata); 
    ci = coefCI(mdl); 
    tmp = table2array(mdl.Coefficients);
     
    contri_val = xdata.*tmp(2:end,1)';
    contri_all = horzcat(contri_all,contri_val);
    
    % constant intercept for empirical models
    emp_const = horzcat(emp_const,tmp(1,1)');
    emp_coef = horzcat(emp_coef,tmp(2:end,1)');
    
end

emp_const = horzcat(emp_const,nanmean(fann_gam.PFT));

cd('/Volumes/RemiLBNL/project6/code4Remi/functions');
subaxis(3,9,[7:9,16:18],'SpacingVert',SpacingVert+0.12,'SpacingHoriz',SpacingHoriz + 0.1,'MR',MR,'ML',ML+0.05,'MarginTop',MarginTop+0.15,'MarginBottom',MarginBottom+0.13); 


groupC2 = [27,120,55;27,120,55;...
    178,24,43;178,24,43;178,24,43;178,24,43;178,24,43;...
    33,102,172;33,102,172;33,102,172;33,102,172;33,102,172;33,102,172;33,102,172;33,102,172]/255;
    
hold on;

x = 1:15;

for i = 1:15
    boxchart(x(i)*ones(size(contri_all(:,i))), contri_all(:,i), 'BoxFaceColor', groupC2(i,:),...
        'BoxFaceAlpha', 0.3, 'WhiskerLineColor',groupC2(i,:), 'MarkerStyle' ,'none')
end
          
hline = refline([0 0]);
hline.Color = 'k';

set(gca, 'box', 'off');

ylim([-27,10]);

factor_names = {'LPC';'LMA';...
    'Tair';'PP';'PAR';'VPD';'SWC';...
    'soilC';'soilN';'CN';'PH';'sand';'silt';'bulkD';'CEC'};

ylabel('\DeltafLNR (%)','FontSize',9);

h = gca;
h.XTick = 1:length(factor_names);
h.XTickLabel = factor_names;
h.TickLabelInterpreter = 'none';


camroll(-90);




%%%% The thrid part, key variables for each PFT

% fit the multi-variate linear regression
%initiate the y data
ydata_raw = nan(size(reshape(vData.fLNR,[],1)));
ydata_raw_pft = nan(size(reshape(vData.fLNR,[],1)));

%%%% get the index for PCAs that have values
tmp_xdata = horzcat(gam_pcas,reshape(EULC_2D,[],1),reshape(vData.annVc,[],1));%
nan_ind = any(isnan(tmp_xdata),2);

LC1D = reshape(EULC_2D,[],1);
ld_text={'CRO','DBF','EBF','ENF','MF','GRA','SH','WET'};


%%%%

contri_all = [];
contri_all = nan(8,10500,15);
contri_all_std = [];

emp_coef = [];
emp_const = [];

for j = 1:3 % for each type of environment factor
    
    % ydata is the total response in fLNR quantified by PCA, xdata is the individual variables
    switch j
        case 1
            
            xdata_raw = horzcat(reshape(vData.Pa,[],1),reshape(vData.LMA,[],1));
            ydata_raw(~nan_ind) = fann_gam.s_leaf1_ + fann_gam.s_leaf2_ + fann_gam.s_leaf3_;

        case 2
           
            xdata_raw = grid_predictors(:,6:10);
            ydata_raw(~nan_ind) = fann_gam.s_clim1_+ fann_gam.s_clim2_ + fann_gam.s_clim3_;
           
        case 3
        
            xdata_raw = grid_predictors(:,[13:20]);
            ydata_raw(~nan_ind) = fann_gam.s_soil1_+ fann_gam.s_soil2_ + fann_gam.s_soil3_;
            
    end

    %xdata_std = nanstd(xdata_raw)

    for ll = 1:8 % for every PFT

        contri_val = [];

        LC_ind = LC1D == ll & nan_ind < 1;
        
        disp(sum(LC_ind));

        xdata =xdata_raw(LC_ind,:);
        ydata =ydata_raw(LC_ind);

        xdata(xdata == -Inf) = NaN;

        mdl = fitlm(xdata,ydata);
        ci = coefCI(mdl);
        tmp = table2array(mdl.Coefficients); % coefficients for each x variables

        contri_val = xdata.*tmp(2:end,1)';
      

        if j == 1
            contri_all(ll,1:size(contri_val,1),1:2) = contri_val;
        elseif j == 2
            contri_all(ll,1:size(contri_val,1),3:7) = contri_val;
        else
            contri_all(ll,1:size(contri_val,1),8:15) = contri_val;
        end
         
        contri_val= nanstd(xdata,1,1).*tmp(2:end,1)';
        
         
        % constant intercept for empirical models
        if j == 1
            emp_coef(ll,1:2) = tmp(2:end,1)';
            emp_const(ll,j) = tmp(1,1);
        elseif j == 2
            emp_coef(ll,3:7) = tmp(2:end,1)';
            emp_const(ll,j) = tmp(1,1);
        else
            emp_coef(ll,8:15) = tmp(2:end,1)';
            emp_const(ll,j) = tmp(1,1);
        end
         
        ydata_raw_pft(~nan_ind) = fann_gam.PFT;
        emp_const(ll,4) = nanmean(ydata_raw_pft(LC_ind));

    end
    
  
end

% plot response
for ll = 1:8

    cd('/Volumes/RemiLBNL/project6/code4Remi/functions');
    if ll == 1
        subaxis(3,9,[20],'SpacingVert',SpacingVert,'SpacingHoriz',SpacingHoriz-0.01,'MR',MR,'ML',ML-0.02,'MarginTop',MarginTop,'MarginBottom',MarginBottom); 

    else
        subaxis(3,9,ll+19,'SpacingVert',SpacingVert,'SpacingHoriz',SpacingHoriz,'MR',MR,'ML',ML,'MarginTop',MarginTop,'MarginBottom',MarginBottom); 

    end

    contri_pft = squeeze(contri_all(ll,:,:));

    groupC2 = [27,120,55;27,120,55;...
    178,24,43;178,24,43;178,24,43;178,24,43;178,24,43;...
    33,102,172;33,102,172;33,102,172;33,102,172;33,102,172;33,102,172;33,102,172;33,102,172]/255;
    
    hold on;

    x = 1:15;

    for i = 1:15
        boxchart(x(i)*ones(size(contri_pft(:,i))), contri_pft(:,i), 'BoxFaceColor', groupC2(i,:),...
            'BoxFaceAlpha', 0.3, 'WhiskerLineColor',groupC2(i,:), 'MarkerStyle' ,'none')
    end

    hline = refline([0 0]);
    hline.Color = 'k';
    set(gca, 'box', 'off');

    ylim([-25, 10]);

    factor_names = {'LPC';'LMA';...
        'Tair';'PP';'PAR';'VPD';'SWC';...
        'soilC';'soilN';'CN';'PH';'sand';'silt';'bulkD';'CEC'};
    
    ylabel(ld_text{ll},'FontSize',9);
    h = gca;
    h.XTick = 1:length(factor_names);
    h.TickLabelInterpreter = 'none';
    if ll==1
        h.XTickLabel = factor_names;
    else
        set(gca,'xtick',[]);
        h.XAxis.Visible = 'off';
    end

    camroll(-90);
end



saveas(gcf,'/Volumes/RemiLBNL/project11_Vcmax_val/Figures/Figure_3.pdf');
    

set(f3,'PaperPositionMode','auto');
print('/Volumes/RemiLBNL/project11_Vcmax_val/Figures/Figure_3','-djpeg','-r600');  


%%%%%%%%%%%% End: Figure 3, Response of fLNR to leaf traits%%%%%%%%%%%%%%%%%%








%%%%%%%%%%%% Figure S1, Spatial distribution of Vcmax25 %%%%%%%%%%%%%%%%%

fs1=figure('Name','Gridded Vcmax','Units', 'centimeters','Color','white', 'Position', [2, 2, 20, 10], ...
    'OuterPosition', [2, 2, 20, 10]);

lat=-89.75:0.5:89.75;
lon=-179.75:0.5:179.75;

colorscheme = [84,48,5;...
            140,81,10;...
            191,129,45;...
            223,194,125;...
            246,232,195;...
            199,234,229;...
            128,205,193;...
            53,151,143;...
            1,102,94;...
            0,60,48];
colorscheme = colorscheme./255;

for i = 1: 8
    
    vcmax2 = squeeze(vData.Vall(:,:,i));
    
    mapdata=flip(vcmax2);
    
    
    cd('/Volumes/RemiLBNL/project6/code4Remi/functions');
    axis1=subaxis(2,4,i,'SpacingVert',-0.0,'SpacingHoriz',-0.04,'MR',0.05,'ML',0.05,'MarginTop',.02,'MarginBottom',.1); 
       
            
    cd('/Volumes/RemiLBNL/project6/code4Remi/functions');
    mapdata= smooth2a(mapdata,3,3); 
    mapdata(180,359)=-50;
    mapdata(180,360)=50;


    cd('/Volumes/RemiLBNL/project6/code4Remi/functions/m_map');
    m_proj('Miller','lon',[-180 180],'lat',[-60 90]); % map projections, range of lon and lat of your data
    m_coast('linewidth',0.5,'color',[0.2 0.2 0.2]); % coast line settings
    hold on;
    m_pcolor(lon,lat,mapdata); % draw your map here
    shading INTERP;   % can be flat or INTERP

    colormap(gca,colorscheme);


    caxis([0 150]); 
    
    h=colorbar;%
    set(h,'YTick',0:30:150);
    ylabel(h,{'Vc_{max}^2^5 (\mumol m^-^2s^-^1)'},'FontSize',8);
    cd('/Volumes/RemiLBNL/project6/code4Remi/functions');cbfreeze(h);

    if mod(i,4) > 0
        h.Visible = 'off';
    end

    cd('/Volumes/RemiLBNL/project6/code4Remi/functions/m_map');
    
    if i > 5     
        m_grid('box','off','tickdir','in','yticklabels',[],'fontsize',5);
    elseif i == 1 
        m_grid('box','off','tickdir','in','xticklabels',[],'fontsize',5);
    elseif i == 5
        m_grid('box','off','tickdir','in','fontsize',5);
    else
        m_grid('box','off','tickdir','in','yticklabels',[],'xticklabels',[],'fontsize',5);
    end
    

    m_text(-170,85,source_names(i),'fontsize',8);

    hold off;
    
    
end


set(fs1,'PaperPositionMode','auto');
print('/Volumes/RemiLBNL/project11_Vcmax_val/Figures/Figure_s1','-djpeg','-r600');

%%%%%%%%%%%% End: Figure S1, Spatial distribution of Vcmax25 %%%%%%%%%%%%%%%%%






%%%%%%%%%%%% Figure S2, Spatial distribution of fLNR %%%%%%%%%%%%%%%%%

fs2=figure('Name','Gridded fLNR','Units', 'centimeters','Color','white', 'Position', [2, 2, 20, 10], ...
    'OuterPosition', [2, 2, 20, 10]);

lat=-89.75:0.5:89.75;
lon=-179.75:0.5:179.75;

for i = 1: 8
    
    vcmax2 = squeeze(vData.Fall(:,:,i));
    
    mapdata=flip(vcmax2);
    
    cd('/Volumes/RemiLBNL/project6/code4Remi/functions');
    axis1=subaxis(2,4,i,'SpacingVert',-0.0,'SpacingHoriz',-0.04,'MR',0.05,'ML',0.05,'MarginTop',.02,'MarginBottom',.1); 
       
    mapdata= smooth2a(mapdata,3,3); 
    mapdata(180,359)=-50;
    mapdata(180,360)=50;

    cd('/Volumes/RemiLBNL/project6/code4Remi/functions/m_map');
    m_proj('Miller','lon',[-180 180],'lat',[-60 90]); % map projections, range of lon and lat of your data
    m_coast('linewidth',0.5,'color',[0.2 0.2 0.2]); % coast line settings
    hold on;
    m_pcolor(lon,lat,mapdata); % draw your map here
    shading INTERP;   % can be flat or INTERP

    colormap(gca,colorscheme);


    caxis([0 50]); 
    
    h=colorbar;%
    set(h,'YTick',0:10:50);
    ylabel(h,{'fLNR (%)'},'FontSize',8);
    cd('/Volumes/RemiLBNL/project6/code4Remi/functions');cbfreeze(h);

    if mod(i,4) > 0
        h.Visible = 'off';
    end

    cd('/Volumes/RemiLBNL/project6/code4Remi/functions/m_map');
    
    if i > 5     
        m_grid('box','off','tickdir','in','yticklabels',[],'fontsize',5);
    elseif i == 1 
        m_grid('box','off','tickdir','in','xticklabels',[],'fontsize',5);
    elseif i == 5
        m_grid('box','off','tickdir','in','fontsize',5);
    else
        m_grid('box','off','tickdir','in','yticklabels',[],'xticklabels',[],'fontsize',5);
    end
    

    m_text(-170,85,source_names(i),'fontsize',8);

    hold off;
    
    
end



set(fs2,'PaperPositionMode','auto');
print('/Volumes/RemiLBNL/project11_Vcmax_val/Figures/Figure_s2','-djpeg','-r600');

%%%%%%%%%%%% End: Figure S2, Spatial distribution of fLNR %%%%%%%%%%%%%%%%%













%%%%%%%%%%%% Figure S5, importance of factors %%%%%%%%%%%%%%%%%%%%%


fs5=figure('Name','Importance of factors','Units', 'centimeters','Color','white', 'Position', [2, 2, 20, 15], ...
    'OuterPosition', [2, 2, 20, 15]);


dataset_names = {'TRY + GlobV + LUNA','TRY', 'GlobV', 'LUNA'};

for j = 1:4
    
    %%%% prepare X and Y
    switch j
        case 1
            cln_all_rec = vData.gfALL;
            site_target = cln_all_rec(:,4);
        case 2
            cln_all_rec = vData.gfTRY;
            site_target = cln_all_rec(:,4);
        case 3
            cln_all_rec = vData.gfGlobV;
            site_target = cln_all_rec(:,4);
        case 4
            cln_all_rec = vData.gfLUNA;
            site_target = cln_all_rec(:,4);
    end
    
    site_predictors = [];
    
    for i = 1 : length(cln_all_rec)
    
        lat_n = round((89.75-cln_all_rec(i,1))./0.5);
        lon_n = round((cln_all_rec(i,2)+179.75)./0.5);


        site_predictors(i,:) = ...
            [vData.Na(lat_n,lon_n),vData.Pa(lat_n,lon_n),vData.SLA(lat_n,lon_n),chl2(lat_n,lon_n),...% leaf traits: N, P, SLA, Chl
        EULC_2D(lat_n,lon_n),...% canopy form: PFT
        env.temp2D(lat_n,lon_n),env.pre2D(lat_n,lon_n), env.par2D(lat_n,lon_n),env.vpd2D(lat_n,lon_n),...
        env.swc2D(lat_n,lon_n),env.alpha2D(lat_n,lon_n), Kp_map(lat_n,lon_n), ...% climate: Tair, PP, PAR, VPD, Koppen climate region
        vData.orgC(lat_n,lon_n),vData.totN(lat_n,lon_n),vData.CN(lat_n,lon_n),vData.PH(lat_n,lon_n),...% soil: orgC, totN, CN, PH
        vData.sand(lat_n,lon_n),vData.silt(lat_n,lon_n),vData.bulkD(lat_n,lon_n),vData.waterC(lat_n,lon_n)...% soil: sand, silt, bulkD, waterC
        ];
    
    end
    
    X = site_predictors; 
    Y = site_target;
    
    %%%%%%%% end prepare X and Y
    
    
    
    t = templateTree('NumVariablesToSample','all',...
    'PredictorSelection','interaction-curvature','Surrogate','on');

    %'interaction-curvature'

    rng(1); % For reproducibility
    Mdl = fitrensemble(X,Y,'Method','Bag','NumLearningCycles',200, ...
        'Learners',t);

    yHat = oobPredict(Mdl);
    R2 = corr(Mdl.Y,yHat)^2;

    impOOB = oobPermutedPredictorImportance(Mdl);
    
    cd('/Volumes/RemiLBNL/project6/code4Remi/functions');
    axis1=subaxis(2,2,j,'SpacingVert',0.05,'SpacingHoriz',0.1,'MR',0.05,'ML',0.1,'MarginTop',.02,'MarginBottom',.15); 
    
    b = bar(impOOB./max(impOOB),'FaceColor',[166 206 227]/255,'EdgeColor','none');
    
    b.FaceColor = 'flat';
    for i = 1:5 % leaf traits
        b.CData(i,:) = [206 227 166]/255;
    end

    for i = 6:12
        b.CData(i,:) = [227 166 206]/255;
    end

    for i = 13:20
        b.CData(i,:) = [166 206 227]/255;
    end
    
    set(gca, 'box', 'off');


    
    h = gca;
    %h.XTickLabel = Mdl.PredictorNames;
    h.XTick = 1:length(predictor_names);
    h.FontSize = 8;
    
    if j > 2
        h.XTickLabel = predictor_names;
        h.FontSize = 8;
        h.XTickLabelRotation = 60;
        h.TickLabelInterpreter = 'none';
        xlabel('Predictor variable','FontSize',10);
    else
        h.XTickLabel ={''};
    end
    
    ylabel({'Normalised importance'; 'of predictors'},'FontSize',10);
    
    
    ylim([0 1]);
    
    text(0.2,0.95,strcat(char(96+j),{' '},dataset_names{j}),'FontWeight','bold');
    
end


set(fs5,'PaperPositionMode','auto');
print('/Volumes/RemiLBNL/project11_Vcmax_val/Figures/Figure_s5','-djpeg','-r600');

%%%%%%%%%%%% End: Figure S5, importance of factors %%%%%%%%%%%%%%%%%%%%%












%%%%%%%%%%%% Figure s4, Univariate response of fLNR to all variables%%%%%%%%%%%%%%%%%%


fs4=figure('Name','Univariate response','Units', 'centimeters','Color','white', 'Position', [2, 2, 15, 18], ...
    'OuterPosition', [2, 2, 15, 18]);
%
xdata_raw = ...
    horzcat(reshape(vData.Na,[],1),reshape(vData.Pa,[],1),reshape(vData.LMA,[],1),...% leaf traits: N, P, SLA, Chl
    reshape(env.temp2D,[],1),reshape(env.pre2D,[],1), reshape(env.par2D,[],1),reshape(env.vpd2D,[],1),...
    reshape(env.swc2D,[],1),reshape(env.alpha2D,[],1), ...% climate: Tair, PP, PAR, VPD
    reshape(vData.orgC,[],1),reshape(vData.totN,[],1),reshape(vData.CN,[],1),reshape(vData.PH,[],1),...% soil: orgC, totN, CN, PH
    reshape(vData.sand,[],1),reshape(vData.silt,[],1),reshape(vData.bulkD,[],1),reshape(vData.waterC,[],1)...% soil: sand, silt, bulkD, waterC
    );

variable_names = {'LNC';'LPC';'LMA';...
    'Tair';'PP';'PAR';'VPD';'SWC';'alpha';...
    'soilC';'soilN';'CN';'PH';'sand';'silt';'bulkD';'CEC'};

ydata_raw = reshape(vData.fLNR,[],1);

tmp1 = reshape(EULC_2D,[],1);
L_ind = tmp1 > 0;
xdata = xdata_raw(L_ind,:);
ydata = ydata_raw(L_ind);

for j = 1:size(xdata,2) % for each leaf nutrient indicator 
    
    
    ydata_f = log(ydata);
    xdata_f = xdata(:,j);
    
    xmin = nanmin(xdata_f);
    xmax = nanmax(xdata_f);
    
    ymin = nanmin(ydata_f);
    ymax = nanmax(ydata_f);

    xdata_all=[ones(size(ydata_f)) xdata_f];


    % each panel
    cd('/Volumes/RemiLBNL/project6/code4Remi/functions');
    axis1=subaxis(5,4,j,'SpacingVert',0.08,'SpacingHoriz',0.06,'MR',0.08,'ML',0.1,'MarginTop',.02,'MarginBottom',.1); 

    %regression
    [b,bint,~,~,~] = regress(ydata_f,xdata_all,0.05); 

    mdl = fitlm(xdata_f,ydata_f);
    [ypred,yci] = predict(mdl);
    ydata_low=yci(:,1);
    ydata_up=yci(:,2);


    % gradient map, that can replace scatter
    fg_z = [];
    inter_n = 100;
    fg_x=xmin:((xmax-xmin)/inter_n):xmax;
    fg_y=ymin:((ymax-ymin)/inter_n):ymax;
    for c_x=1:inter_n
        for c_y=1:inter_n
            c_ind=xdata_f>fg_x(c_x) & xdata_f<fg_x(c_x+1)...
                & ydata_f>fg_y(c_y) & ydata_f<fg_y(c_y+1);

            fg_z(c_y,c_x)=sum(c_ind);
        end
    end

    fg_z(fg_z==0)=NaN;

    h=pcolor(fg_x(1:100),fg_y(1:100),fg_z);%if interval too large, then the color cannot show
    set(h, 'EdgeColor', 'none');

    color_s=flip(gray(10));
    colormap(gca,color_s(3:10,:)); 
    hold on;



    % fitted line
    pxx=plot(xdata_f, xdata_f.*b(2)+b(1),'--','color','r');

    hold on;

    [r1,r2]=corrcoef(xdata_f,ydata_f,'rows','complete');
    sig=round(r2(2)*100)./100; 


    r = round(r1(1,2).*100)./100;
    b = round(b.*100)./100;

    hold on;

    %text(1.2*nanmedian(xdata),4,strcat('r=',num2str(r)),'Color','r');
    if sig < 0.01
        text(xmin+0.02*(xmax-xmin),4.5,strcat('r=',num2str(r),'(p<0.01) '),'Color','r','FontSize',8);
    else
        text(xmin+0.02*(xmax-xmin),4.5,strcat('r=',num2str(r),'(p=',num2str(sig),')'),'Color','r','FontSize',8);
    end

    text(xmin+0.02*(xmax-xmin),3.8,strcat('y =',num2str(b(2)),'x+',num2str(b(1))),'Color','r','FontSize',8);

    ylim([1 5]);
    
    h = gca;
    h.YTick = [log(10),log(50),log(100)];
    h.YTickLabel = {'10','50','100'};
    
    if mod(j,4) == 1
        ylabel('fLNR(%)','FontSize',9);
    end
    xlabel(variable_names(j),'FontSize',9);


    
end

set(fs4,'PaperPositionMode','auto');
print('/Volumes/RemiLBNL/project11_Vcmax_val/Figures/Figure_s4','-djpeg','-r600');  


%%%%%%%%%%%% End: Figure s4, Response of fLNR to many variables %%%%%%%%%%%%%%%%%%






%%%%%% Figure s3: Uncertainty of fLNR %%%%%%%%%%%%%%%%%%%%%%%%

f9=figure('Name','unc of fLNR','Units', 'centimeters','Color','white', 'Position', [2, 2, 15, 14], ...
    'OuterPosition', [2, 2, 15, 14]);

lat=-89.75:0.5:89.75;
lon=-179.75:0.5:179.75;


% sub_panel figure
SpacingVert = 0.00;
SpacingHoriz = 0.05;
MR = 0.05;
ML = 0.01;
MarginTop = 0.01;
MarginBottom = 0.02;


for i = 1:5
    
    switch i
        case 1
            vcmax2 = vData.fLNR_std;    
        case 5
            vcmax2 = vData.fLNR_std_a25;    
        case 4
            vcmax2 = vData.fLNR_std_fNR;    
        case 3
            vcmax2 = vData.fLNR_std_LNC;   
        case 2
            vcmax2 = vData.fLNR_std_VC;   
    end
    
    disp(nanmean(reshape(vcmax2,[],1)));
    disp(nanstd(reshape(vcmax2,[],1)));
    
    %%%%% map of fLNR
    cd('/Volumes/RemiLBNL/project6/code4Remi/functions');
    subaxis(3,2,i,'SpacingVert',SpacingVert,'SpacingHoriz',SpacingHoriz,'MR',MR,'ML',ML,'MarginTop',MarginTop,'MarginBottom',MarginBottom); 

    
    mapdata=flip(vcmax2);

    cd('/Volumes/RemiLBNL/project6/code4Remi/functions');
    mapdata= smooth2a(mapdata,3,3); 
    mapdata(180,359)=-50;
    mapdata(180,360)=50;

    cd('/Volumes/RemiLBNL/project6/code4Remi/functions/m_map');
    m_proj('Miller','lon',[-180 180],'lat',[-60 90]); % map projections, range of lon and lat of your data
    m_coast('linewidth',0.5,'color',[0.2 0.2 0.2]); % coast line settings
    hold on;
    m_pcolor(lon,lat,mapdata); % draw your map here
    shading INTERP;   % can be flat or INTERP

    colorscheme = [84,48,5;...
                140,81,10;...
                191,129,45;...
                223,194,125;...
                246,232,195;...
                199,234,229;...
                128,205,193;...
                53,151,143;...
                1,102,94;...
                0,60,48];
    colorscheme = colorscheme./255;  

    colormap(gca,colorscheme);

    caxis([0 10]); 

    h=colorbar;%
    set(h,'YTick',0:2:10);
    ylabel(h,{'unc. of fLNR (%)'},'FontSize',8);
    cd('/Volumes/RemiLBNL/project6/code4Remi/functions');cbfreeze(h);

    cd('/Volumes/RemiLBNL/project6/code4Remi/functions/m_map');
    m_grid('box','off','tickdir','in','yticklabels',[],'xticklabels',[],'fontsize',7);


   m_text(-170,85,char(96+i),'fontsize',10,'FontWeight','bold');
    
    hold off;
end


set(fs5,'PaperPositionMode','auto');
print('/Volumes/RemiLBNL/project11_Vcmax_val/Figures/Figure_s5','-djpeg','-r600');



%%%%%% End: Figure s5. Uncertainty of fLNR %%%%%%%%%%%%%%%%%%%%%%%%






%%%%%%% Figure s6: distribution of obs %%%%%%%%%%% 
% plot the figure showing the distribution of TRY and GlobV dataset

load('/Volumes/RemiLBNL/project11_Vcmax_val/data/fData.mat');
 
fs6=figure('Name','obs. distribution','Units', 'centimeters','Color','white', 'Position', [2, 2, 15, 10], ...
    'OuterPosition', [2, 2, 15, 10]);


for i = 1 : 2
    
    switch i
        case 1      
            target_rec = vData.gfALL; 
            
        case 2
            target_rec = fData.gfALL;
    end
    
    [target_u, ia, ic] = unique(target_rec(:,[1,2]), 'rows');
    disp(length(ia));
    
    for uu = 1:length(target_u)
        
        Lia = ismember(target_rec(:,[1,2]),target_u(uu,1:2),'rows');
        target_u(uu,3) = sum(Lia);
        
    end
    
    
    
    cd('/Volumes/RemiLBNL/project6/code4Remi/functions');
    axis1=subaxis(2,2,i,'SpacingVert',0.1,'SpacingHoriz',0.03,'MR',0.05,'ML',0.05,'MarginTop',.05,'MarginBottom',.05); 
    
    cd('/Volumes/RemiLBNL/project6/code4Remi/functions/m_map');
    m_proj('Miller','lon',[-180 180],'lat',[-60 90]); % map projections, range of lon and lat of your data
    m_coast('linewidth',0.5,'color',[0.2 0.2 0.2]); % coast line settings
    
    
    target_size = max(target_u(:,3)/5,0.5);
    target_size = min(target_u(:,3)/5,8);

    
    for uu = 1:length(target_u)
           h1=m_line(target_u(uu,2), target_u(uu,1),'marker','o','color',[1,0.2,0.2],'linewi',0.5,...
          'linest','none','markersize',target_size(uu,1),'markerfacecolor','w');
      
           hold on;
    end

    
    m_text(-170,85,char(96+i),'fontsize',10,'FontWeight','bold');
    m_grid('box','off','tickdir','in','yticklabels',[],'xticklabels',[],'fontsize',7);

    
end

for i = 1:2
    
    switch i
        case 1
            vcmax2 = vData.fLNR;
        case 2
            vcmax2 = fData.annVc;
    end
    
    cd('/Volumes/RemiLBNL/project6/code4Remi/functions');
    subaxis(2,2,2+i,'SpacingVert',0.1,'SpacingHoriz',0.05,'MR',0.05,'ML',0.1,'MarginTop',.05,'MarginBottom',.05);
    
    mapdata=flip(vcmax2);

    cd('/Volumes/RemiLBNL/project6/code4Remi/functions');
    mapdata= smooth2a(mapdata,3,3); 
    mapdata(180,359)=-50;
    mapdata(180,360)=50;

    cd('/Volumes/RemiLBNL/project6/code4Remi/functions/m_map');
    m_proj('Miller','lon',[-180 180],'lat',[-60 90]); % map projections, range of lon and lat of your data
    m_coast('linewidth',0.5,'color',[0.2 0.2 0.2]); % coast line settings
    hold on;
    m_pcolor(lon,lat,mapdata); % draw your map here
    shading INTERP;   % can be flat or INTERP

    colorscheme = [84,48,5;...
                140,81,10;...
                191,129,45;...
                223,194,125;...
                246,232,195;...
                199,234,229;...
                128,205,193;...
                53,151,143;...
                1,102,94;...
                0,60,48];
    colorscheme = colorscheme./255;   

    colormap(gca,colorscheme);

    caxis([0 50]); 
    
    h=colorbar;%
    set(h,'YTick',0:10:50);
    ylabel(h,{'fLNR (%)'},'FontSize',8);
    cd('/Volumes/RemiLBNL/project6/code4Remi/functions');cbfreeze(h);

    cd('/Volumes/RemiLBNL/project6/code4Remi/functions/m_map');
    m_grid('box','off','tickdir','in','yticklabels',[],'xticklabels',[],'fontsize',7);


    m_text(-170,85,char(96+2+i),'fontsize',10,'FontWeight','bold');
    
    hold off;
end

    

set(fs6,'PaperPositionMode','auto');
print('/Volumes/RemiLBNL/project11_Vcmax_val/Figures/Figure_s6','-djpeg','-r600');
%%%%%%% End: Figure s6: distribution of obs %%%%%%%%%%% 






%%%%%%% Fig s7. RF response curve %%%%%%%%
fs7=figure('Name','obs. distribution','Units', 'centimeters','Color','white', 'Position', [2, 2, 15, 12], ...
    'OuterPosition', [2, 2, 15, 12]);

for i = 1:4
    
    cd('/Volumes/RemiLBNL/project6/code4Remi/functions');
    axis1=subaxis(2,2,i,'SpacingVert',0.13,'SpacingHoriz',0.1,'MR',0.05,'ML',0.1,'MarginTop',.05,'MarginBottom',.12); 
    
    [ydata,xdata,~] = partialDependence(final_RF,i);

    switch i
        case 1 %chl
            plot(xdata,ydata,'-','color','k');
            xlabel('Chl (\mug/cm^2)','Fontsize',9);
        case 2 % PFT
            
            xdata2 = 1:8;
            ydata2 = [ydata([1,18,34,51,67,84,100]);ydata(100)];
            scatter(xdata2,ydata2,20,'.','k');
            ld_text={'CRO','DBF','EBF','ENF','MF','GRA','SH','WET'};
            %set(gca,'XTickLabel',ld_text,'FontSize',9,'box', 'off');
            h = gca;
            h.XTick = 1:length(ld_text);
            h.XTickLabel = ld_text;
            h.XTickLabelRotation = 60;
            h.TickLabelInterpreter = 'none';
            set(gca,'XTickLabel',ld_text,'Fontsize',8);
            
            
        case 3 % PP
            plot(xdata,ydata,'-','color','k');
            xlabel('Precipitation (mm/yr)','Fontsize',9);
        case 4 % soil pH
            plot(xdata,ydata,'-','color','k');
            xlabel('soil pH','Fontsize',9);
    end
    
    set(gca,'box','off');
    ylim([45 75]);
    ylabel('Vc_{max}^2^5 (\mumol m^-^2s^-^1)','Fontsize',9);
end

set(fs7,'PaperPositionMode','auto');
print('/Volumes/RemiLBNL/project11_Vcmax_val/Figures/Figure_s7','-djpeg','-r600');
%%%%%%%% End Fig s7. RF response curve %%%%%%%





%%%%%%%% Fig s8. uncertainty of Vcmax %%%%%%%

fs8=figure('Name','Gridded Vcmax NS_Chl','Units', 'centimeters','Color','white', 'Position', [2, 2, 20, 10], ...
    'OuterPosition', [2, 2, 15, 6]);

lat=-89.75:0.5:89.75;
lon=-179.75:0.5:179.75;

for i = 1: 2
    
   switch i
        case 1
            vcmax2 = squeeze(vData.Vall(:,:,8)); % as the vData.annVC has been changed to fLNR run
        case 2
            vcmax2 = vData.annVC_std;
    end
    
    cd('/Volumes/RemiLBNL/project6/code4Remi/functions');
    subaxis(1,2,i,'SpacingVert',0.1,'SpacingHoriz',0.05,'MR',0.05,'ML',0.05,'MarginTop',.05,'MarginBottom',.05);
    
    mapdata=flip(vcmax2);

    cd('/Volumes/RemiLBNL/project6/code4Remi/functions');
    mapdata= smooth2a(mapdata,3,3); 
    mapdata(180,359)=-50;
    mapdata(180,360)=50;

    cd('/Volumes/RemiLBNL/project6/code4Remi/functions/m_map');
    m_proj('Miller','lon',[-180 180],'lat',[-60 90]); % map projections, range of lon and lat of your data
    m_coast('linewidth',0.5,'color',[0.2 0.2 0.2]); % coast line settings
    hold on;
    m_pcolor(lon,lat,mapdata); % draw your map here
    shading INTERP;   % can be flat or INTERP

    colorscheme = [84,48,5;...
                140,81,10;...
                191,129,45;...
                223,194,125;...
                246,232,195;...
                199,234,229;...
                128,205,193;...
                53,151,143;...
                1,102,94;...
                0,60,48];
    colorscheme = colorscheme./255;     

    colormap(gca,colorscheme);
    caxis([20 120]); 
    h=colorbar;%
    set(h,'YTick',20:20:120);

    
    switch i
        case 1
            caxis([20 120]); 
            h=colorbar;%
            set(h,'YTick',20:20:120);
            ylabel(h,{'Vc_{max}^2^5 (\mumol m^-^2s^-^1)'},'FontSize',8);
        case 2
            caxis([10 60]); 
            h=colorbar;%
            set(h,'YTick',10:10:60);
            ylabel(h,{'Unc. Vc_{max}^2^5 (\mumol m^-^2s^-^1)'},'FontSize',8);
    end
    cd('/Volumes/RemiLBNL/project6/code4Remi/functions');cbfreeze(h);

    cd('/Volumes/RemiLBNL/project6/code4Remi/functions/m_map');
    m_grid('box','off','tickdir','in','yticklabels',[],'xticklabels',[],'fontsize',7);


    m_text(-170,85,char(96+i),'fontsize',10,'FontWeight','bold');
    
    hold off;
    
    
end



set(fs8,'PaperPositionMode','auto');
print('/Volumes/RemiLBNL/project11_Vcmax_val/Figures/Figure_s8','-djpeg','-r600');


%%%%%%%% end: Fig. S8. uncertainty of Vcmax %%%%%%%




%%%%%%%% Fig S9. different LNC product and correlations %%%%%%%

fs9=figure('Name','LNCs','Units', 'centimeters','Color','white', 'Position', [2, 2, 20, 10], ...
    'OuterPosition', [2, 2, 15, 15]);

lat=-89.75:0.5:89.75;
lon=-179.75:0.5:179.75;

for i = 1: 3
    
   switch i
        case 1
            vcmax2 = vData.Na; % EB17
            title = 'EB17';
        case 2
            vcmax2 =AMM.Na; % AMM18
            title = 'AMM18';
        case 3
            vcmax2 = CB.Na; % CB20
            title = 'CB20';
    end
    
    cd('/Volumes/RemiLBNL/project6/code4Remi/functions');
    subaxis(3,2,(i-1)*2 + 1,'SpacingVert',0.05,'SpacingHoriz',0.01,'MR',0.05,'ML',0.01,'MarginTop',.02,'MarginBottom',.05);
    
    mapdata=flip(vcmax2);

    cd('/Volumes/RemiLBNL/project6/code4Remi/functions');
    mapdata= smooth2a(mapdata,3,3); 
    mapdata(180,359)=-50;
    mapdata(180,360)=50;

    cd('/Volumes/RemiLBNL/project6/code4Remi/functions/m_map');
    m_proj('Miller','lon',[-180 180],'lat',[-60 90]); % map projections, range of lon and lat of your data
    m_coast('linewidth',0.5,'color',[0.2 0.2 0.2]); % coast line settings
    hold on;
    m_pcolor(lon,lat,mapdata); % draw your map here
    shading INTERP;   % can be flat or INTERP

    colorscheme = [84,48,5;...
                140,81,10;...
                191,129,45;...
                223,194,125;...
                246,232,195;...
                199,234,229;...
                128,205,193;...
                53,151,143;...
                1,102,94;...
                0,60,48];
    colorscheme = colorscheme./255;
    % colorscheme = jet(12);       

    colormap(gca,colorscheme);
    caxis([0 3]); 
    h=colorbar;%
    set(h,'YTick',0:0.5:3);

    ylabel(h,{'LNC (g m^-^2)'},'FontSize',9);
    cd('/Volumes/RemiLBNL/project6/code4Remi/functions');cbfreeze(h);

    cd('/Volumes/RemiLBNL/project6/code4Remi/functions/m_map');
    m_grid('box','off','tickdir','in','yticklabels',[],'xticklabels',[],'fontsize',7);


    m_text(-170,85,strcat(char(96+i),{' '},title),'fontsize',10,'FontWeight','bold');
    
    hold off;
    
    
end


for j = 1:3 % for each leaf nutrient indicator 
    
    switch j 
        case 1
            xdata_f = reshape(vData.Na,[],1);
            ydata_f = reshape(AMM.Na,[],1);
            title = 'EB17 vs AMM18';
        case 2
            xdata_f = reshape(AMM.Na,[],1);
            ydata_f = reshape(CB.Na,[],1);
            title = 'AMM18 vs CB20';
        case 3
            xdata_f = reshape(CB.Na,[],1);
            ydata_f = reshape(vData.Na,[],1);
            title = 'CB20 vs EB17';
    end
    
    xmin = nanmin(xdata_f);
    xmax = nanmax(xdata_f);
    
    ymin = nanmin(ydata_f);
    ymax = nanmax(ydata_f);

    xdata_all=[ones(size(ydata_f)) xdata_f];


    % each panel
    cd('/Volumes/RemiLBNL/project6/code4Remi/functions');
    axis1=subaxis(3,2,j*2,'SpacingVert',0.1,'SpacingHoriz',0.1,'MR',0.03,'ML',0.1,'MarginTop',.05,'MarginBottom',.1); 

    %regression
    [b,bint,~,~,~] = regress(ydata_f,xdata_all,0.05); 

    mdl = fitlm(xdata_f,ydata_f);
    [ypred,yci] = predict(mdl);
    ydata_low=yci(:,1);
    ydata_up=yci(:,2);


    % gradient map, that can replace scatter
    fg_z = [];
    inter_n = 100;
    fg_x=xmin:((xmax-xmin)/inter_n):xmax;
    fg_y=ymin:((ymax-ymin)/inter_n):ymax;
    for c_x=1:inter_n
        for c_y=1:inter_n
            c_ind=xdata_f>fg_x(c_x) & xdata_f<fg_x(c_x+1)...
                & ydata_f>fg_y(c_y) & ydata_f<fg_y(c_y+1);

            fg_z(c_y,c_x)=sum(c_ind);
        end
    end

    fg_z(fg_z==0)=NaN;

    h=pcolor(fg_x(1:100),fg_y(1:100),fg_z);%if interval too large, then the color cannot show
    set(h, 'EdgeColor', 'none');

    color_s=flip(gray(10));
    colormap(gca,color_s(3:10,:)); 
    hold on;



    % fitted line
    pxx=plot(xdata_f, xdata_f.*b(2)+b(1),'--','color','r');

    hold on;

    [r1,r2]=corrcoef(xdata_f,ydata_f,'rows','complete');
    sig=round(r2(2)*100)./100; 


    r = round(r1(1,2).*100)./100;
    b = round(b.*100)./100;

    hold on;

    if sig < 0.01
        text(xmin+0.02*(xmax-xmin),2.8,strcat('r=',num2str(r),'(p<0.01) '),'Color','r','FontSize',8);
    else
        text(xmin+0.02*(xmax-xmin),2.8,strcat('r=',num2str(r),'(p=',num2str(sig),')'),'Color','r','FontSize',8);
    end

    text(xmin+0.02*(xmax-xmin),2.4,strcat('y =',num2str(b(2)),'x+',num2str(b(1))),'Color','r','FontSize',8);

    ylim([0.5 3.5]);
    xlim([0.5 3.5]);
    
    h = gca;
    
    ylabel('LNC (g m^-^2)','FontSize',9);
    
    if j == 3
    xlabel('LNC (g m^-^2)','FontSize',9);
    end

    text(0.55,3.25,strcat(char(99+j),{' '},title),'fontsize',10,'FontWeight','bold');
    
end

set(fs9,'PaperPositionMode','auto');
print('/Volumes/RemiLBNL/project11_Vcmax_val/Figures/Figure_s9','-djpeg','-r600');


%%%%%%%% end: Fig s9. different LNC product and correlations %%%%%%%



%%%%%%%% Fig S12. sample vs globe distribution %%%%%%%

fs12=figure('Name','Sample vs globe distribution','Units', 'centimeters','Color','white', 'Position', [2, 2, 20, 10], ...
    'OuterPosition', [2, 2, 15, 15]);


color_s = [255, 109, 182;  0, 0, 0];
color_s = color_s./255;
count=1;


sample_names = {'LNC';'LPC';'LMA';'Chl';...
    'Tair';'PP';'PAR';'VPD';'SWC';'alpha';...
    'soilC';'soilN';'CN';'PH';'sand';'silt';'bulkD';'CEC'};

sample_data = site_predictors(:,[1:4,6:11,13:end]);
global_data = grid_predictors(:,[1:4,6:11,13:end]);


for i = 1 : size(sample_data,2) % for each variables
    
    
    cd('/Volumes/RemiLBNL/project6/code4Remi/functions');
    axis1=subaxis(5,4,i,'SpacingVert',0.11,'SpacingHoriz',0.06,'MR',0.08,'ML',0.1,'MarginTop',.02,'MarginBottom',.1); 

    
    for j = 1:2
        
        switch j % for samples and global map
            case 1
                xdata = sample_data(:,i);
            case 2
                xdata = global_data(:,i);
        end
        
            ind = ~isnan(xdata);
            xdata2 = xdata(ind);
            
            if ~isempty(xdata2) %sum(~isnan(xdata))>0
                % plot histogram
                 h1=histfit(xdata2,[],'kernel');
                 h1(1).Visible='off';
                 step=h1(2).XData(2)-h1(2).XData(1);
                 h1(2).YData=h1(2).YData./sum(step.*h1(2).YData)./100*1000;
                 h1(2).Color=color_s(j,:);
                 h1(2).LineWidth=1;
                 %h1(2).LineStyle='--';
                 sh(count)=h1(2);
                 count=count+1;

                 hold on;
            end

            
    end
    
    xlabel(sample_names{i},'FontSize',10);
    
    if i == 9
    ylabel('probability (10^-^3)','FontSize',10);
    end
         
end


legend([sh(1),sh(2)],{'sample','globe'},'FontSize',8,'box','off','Orientation','vertical','Location',[0.7,0.08,0.25,0.1]);


set(fs12,'PaperPositionMode','auto');
print('/Volumes/RemiLBNL/project11_Vcmax_val/Figures/Figure_s12','-djpeg','-r600');


%%%%%%%% end: Fig s12. sample vs globe distribution %%%%%%%




%%%%%%%% Fig s11. cross validation + spatial cross validation %%%%%%%

fs11=figure('Name','Cross Validation','Units', 'centimeters','Color','white', 'Position', [2, 2, 15, 8], ...
    'OuterPosition', [2, 2, 15, 8]);

cd('/Volumes/RemiLBNL/project6/code4Remi/functions');
axis1=subaxis(1,2,1,'SpacingVert',0.06,'SpacingHoriz',0.15,'MR',0.1,'ML',0.1,'MarginTop',.05,'MarginBottom',.1); 

b = bar(nanmean(r2_test_all),'FaceColor',[166 206 227]/255,'EdgeColor','none');

b.FaceColor = 'flat';

b.CData(1,:) = [206 227 166]/255;
b.CData(2,:) = [227 166 206]/255;
b.CData(3,:) = [166 206 227]/255;


hold on;

er = errorbar(nanmean(r2_test_all),nanstd(r2_test_all));    
er.Color = [0 0 0];                            
er.LineStyle = 'none';  
er.CapSize = 10;

hold on;
scatter(ones(30,1)+0.2,r2_test_all(:,1),10,'k','filled');
scatter(ones(30,1)*2+0.2,r2_test_all(:,2),10,'k','filled');
scatter(ones(30,1)*3+0.2,r2_test_all(:,3),10,'k','filled');

set(gca, 'box', 'off');

ylim([0,1]);

factor_names = {'C.V.';...
    'S.C.V.';'All'};

ylabel('R^2','FontSize',9);

h = gca;
h.XTick = 1:length(factor_names);
h.XTickLabel = factor_names;
h.TickLabelInterpreter = 'none';
h.FontSize = 8;

text(0.3,0.9,char(96+1),'fontsize',10,'FontWeight','bold');

%RMSE
cd('/Volumes/RemiLBNL/project6/code4Remi/functions');
axis1=subaxis(1,2,2,'SpacingVert',0.06,'SpacingHoriz',0.15,'MR',0.1,'ML',0.1,'MarginTop',.05,'MarginBottom',.1); 

b = bar(nanmean(rmse_test_all),'FaceColor',[166 206 227]/255,'EdgeColor','none');

b.FaceColor = 'flat';

b.CData(1,:) = [206 227 166]/255;
b.CData(2,:) = [227 166 206]/255;
b.CData(3,:) = [166 206 227]/255;


hold on;

er = errorbar(nanmean(rmse_test_all),nanstd(rmse_test_all));    
er.Color = [0 0 0];                            
er.LineStyle = 'none';  
er.CapSize = 10;

hold on;
scatter(ones(30,1)+0.2,rmse_test_all(:,1),10,'k','filled');
scatter(ones(30,1)*2+0.2,rmse_test_all(:,2),10,'k','filled');
scatter(ones(30,1)*3+0.2,rmse_test_all(:,3),10,'k','filled');

set(gca, 'box', 'off');

ylim([0,60]);

factor_names = {'C.V.';...
    'S.C.V.';'All'};

ylabel('RMSE (\mumol m^-^2s^-^1)','FontSize',9);

h = gca;
h.XTick = 1:length(factor_names);
h.XTickLabel = factor_names;
h.TickLabelInterpreter = 'none';
h.FontSize = 8;

text(0.3,55,char(97+1),'fontsize',10,'FontWeight','bold');

set(fs11,'PaperPositionMode','auto');
print('/Volumes/RemiLBNL/project11_Vcmax_val/Figures/Figure_s11','-djpeg','-r600');






%%%%%%%%%% Fig. S10. using different LNC to get fLNR %%%%%%%%%%%%%%%%%%%%%%%%%

fr4=figure('Name','Gridded Vcmax NS_Chl','Units', 'centimeters','Color','white', 'Position', [2, 2, 20, 15], ...
    'OuterPosition', [2, 2, 20, 18]);

lat=-89.75:0.5:89.75;
lon=-179.75:0.5:179.75;

fLNR1 = vData.annVc./vData.Na./alpha25./fNR * 100; % EB17
fLNR2 = vData.annVc./AMM.Na./alpha25./fNR * 100;% AMM18
fLNR3 = vData.annVc./CB.Na./alpha25./fNR * 100; % CB20

for i = 1: 3
    
   switch i
        case 1
            vcmax2 = fLNR1; % EB17
            title = 'EB17';
        case 2
            vcmax2 = fLNR2;% AMM18
            title = 'AMM18';
        case 3
            vcmax2 = fLNR3; % CB20
            title = 'CB20';
    end
    
    cd('/Volumes/RemiLBNL/project6/code4Remi/functions');
    subaxis(3,3,(i-1)*3 + 1,'SpacingVert',0.05,'SpacingHoriz',0.01,'MR',0.05,'ML',0.01,'MarginTop',.02,'MarginBottom',.05);
    
    mapdata=flip(vcmax2);

    cd('/Volumes/RemiLBNL/project6/code4Remi/functions');
    mapdata= smooth2a(mapdata,3,3); 
    mapdata(180,359)=-50;
    mapdata(180,360)=50;

    cd('/Volumes/RemiLBNL/project6/code4Remi/functions/m_map');
    m_proj('Miller','lon',[-180 180],'lat',[-60 90]); % map projections, range of lon and lat of your data
    m_coast('linewidth',0.5,'color',[0.2 0.2 0.2]); % coast line settings
    hold on;
    m_pcolor(lon,lat,mapdata); % draw your map here
    shading INTERP;   % can be flat or INTERP

    colorscheme = [84,48,5;...
                140,81,10;...
                191,129,45;...
                223,194,125;...
                246,232,195;...
                199,234,229;...
                128,205,193;...
                53,151,143;...
                1,102,94;...
                0,60,48];
    colorscheme = colorscheme./255;
    % colorscheme = jet(12);       

    colormap(gca,colorscheme);
    caxis([0 50]); 
    h=colorbar;%
    set(h,'YTick',0:10:50);

    ylabel(h,{'fLNR (%)'},'FontSize',9);
    cd('/Volumes/RemiLBNL/project6/code4Remi/functions');cbfreeze(h);

    cd('/Volumes/RemiLBNL/project6/code4Remi/functions/m_map');
    %m_grid('box','off','tickdir','in','fontsize',7);
    m_grid('box','off','tickdir','in','yticklabels',[],'xticklabels',[],'fontsize',7);


    m_text(-170,85,strcat(char(96+i),{' '},title),'fontsize',10,'FontWeight','bold');
    
    hold off;
    
    
end


for j = 1:3 % for each leaf nutrient indicator 
    
    switch j 
        case 1
            xdata_f = reshape(fLNR1,[],1);
            ydata_f = reshape(fLNR2,[],1);
            title = 'EB17 vs AMM18';
        case 2
            xdata_f = reshape(fLNR2,[],1);
            ydata_f = reshape(fLNR3,[],1);
            title = 'AMM18 vs CB20';
        case 3
            xdata_f = reshape(fLNR3,[],1);
            ydata_f = reshape(fLNR1,[],1);
            title = 'CB20 vs EB17';
    end
    
    xmin = nanmin(xdata_f);
    xmax = nanmax(xdata_f);
    
    ymin = nanmin(ydata_f);
    ymax = nanmax(ydata_f);

    xdata_all=[ones(size(ydata_f)) xdata_f];


    % each panel
    cd('/Volumes/RemiLBNL/project6/code4Remi/functions');
    axis1=subaxis(3,3,(j-1)*3 + 2,'SpacingVert',0.1,'SpacingHoriz',0.1,'MR',0.03,'ML',0.1,'MarginTop',.05,'MarginBottom',.1); 

    %regression
    [b,bint,~,~,~] = regress(ydata_f,xdata_all,0.05); 

    mdl = fitlm(xdata_f,ydata_f);
    [ypred,yci] = predict(mdl);
    ydata_low=yci(:,1);
    ydata_up=yci(:,2);


    % gradient map, that can replace scatter
    fg_z = [];
    inter_n = 100;
    fg_x=xmin:((xmax-xmin)/inter_n):xmax;
    fg_y=ymin:((ymax-ymin)/inter_n):ymax;
    for c_x=1:inter_n
        for c_y=1:inter_n
            c_ind=xdata_f>fg_x(c_x) & xdata_f<fg_x(c_x+1)...
                & ydata_f>fg_y(c_y) & ydata_f<fg_y(c_y+1);

            fg_z(c_y,c_x)=sum(c_ind);
        end
    end

    fg_z(fg_z==0)=NaN;

    h=pcolor(fg_x(1:100),fg_y(1:100),fg_z);%if interval too large, then the color cannot show
    set(h, 'EdgeColor', 'none');

    color_s=flip(gray(10));
    colormap(gca,color_s(3:10,:)); 
    hold on;



    % fitted line
    pxx=plot(xdata_f, xdata_f.*b(2)+b(1),'--','color','r');

    hold on;

    [r1,r2]=corrcoef(xdata_f,ydata_f,'rows','complete');
    sig=round(r2(2)*100)./100; 


    r = round(r1(1,2).*100)./100;
    b = round(b.*100)./100;

    hold on;

    %text(1.2*nanmedian(xdata),4,strcat('r=',num2str(r)),'Color','r');
    if sig < 0.01
        text(xmin+0.02*(xmax-xmin),35,strcat('r=',num2str(r),'(p<0.01) '),'Color','r','FontSize',8);
    else
        text(xmin+0.02*(xmax-xmin),35,strcat('r=',num2str(r),'(p=',num2str(sig),')'),'Color','r','FontSize',8);
    end

    text(xmin+0.02*(xmax-xmin),30,strcat('y =',num2str(b(2)),'x+',num2str(b(1))),'Color','r','FontSize',8);

    ylim([0 50]);
    xlim([0 50]);
    
    h = gca;
    
    ylabel('fLNR (%)','FontSize',9);
    
    if j == 3
    xlabel('fLNR (%)','FontSize',9);
    end

    text(0.55,45,strcat(char(99+j),{' '},title),'fontsize',10,'FontWeight','bold');
    
end



% % fit the multi-variate linear regression

for nn = 1:3
    
    switch nn
        case 1
            filename = '/Volumes/RemiLBNL/project11_Vcmax_val/data/fann_pca_gam_dec.csv';
            fann_gam = readtable(filename);
            %%%% get the index for PCAs that have values
            tmp_xdata = horzcat(gam_pcas,reshape(EULC_2D,[],1),reshape(vData.annVc,[],1),reshape(fLNR1,[],1));%
            title = 'EB17';

        case 2
            filename = '/Volumes/RemiLBNL/project11_Vcmax_val/data/fAMM_pca_gam_altN1.csv';
            fann_gam = readtable(filename);
            %%%% get the index for PCAs that have values
            tmp_xdata = horzcat(gam_pcas,reshape(EULC_2D,[],1),reshape(vData.annVc,[],1),reshape(fLNR2,[],1));%
            title = 'AMM18';

        case 3
             filename = '/Volumes/RemiLBNL/project11_Vcmax_val/data/fCB_pca_gam_altN1.csv';
             fann_gam = readtable(filename);
             %%%% get the index for PCAs that have values
            tmp_xdata = horzcat(gam_pcas,reshape(EULC_2D,[],1),reshape(vData.annVc,[],1),reshape(fLNR3,[],1));%
            title = 'CB20';

    end
    

    ydata_raw = nan(size(reshape(vData.fLNR,[],1)));

    
    nan_ind = any(isnan(tmp_xdata),2);
    ydata_pft = nan(size(ydata_raw));
    ydata_pft(~nan_ind) = fann_gam.PFT;

    %%%%

    contri_all = [];
    contri_all_std = [];
    emp_coef = [];
    emp_const = [];

    for j = 1:3 % for each leaf nutrient indicator 


        switch j
            case 1
    %         
                xdata_raw = horzcat(reshape(vData.Pa,[],1),reshape(vData.LMA,[],1));

                ydata_raw(~nan_ind) = fann_gam.s_leaf1_ + fann_gam.s_leaf2_ + fann_gam.s_leaf3_;

            case 2

                xdata_raw = grid_predictors(:,6:10);

                ydata_raw(~nan_ind) = fann_gam.s_clim1_+ fann_gam.s_clim2_ + fann_gam.s_clim3_;


            case 3
    %         
                xdata_raw = grid_predictors(:,[13:20]); % 20: 27

                ydata_raw(~nan_ind) = fann_gam.s_soil1_+ fann_gam.s_soil2_ + fann_gam.s_soil3_;


        end



        ydata_raw = ydata_raw; %nanmean(ydata_raw); %ydata_pft;


        xdata =xdata_raw(~nan_ind,:);
        ydata =ydata_raw(~nan_ind);

        xdata(xdata == -Inf) = NaN;

        mdl = fitlm(xdata,ydata);
        ci = coefCI(mdl);
        tmp = table2array(mdl.Coefficients);
        contri_val= xdata.*tmp(2:end,1)';
        contri_all = horzcat(contri_all,contri_val);

        contri_val= nanstd(xdata,1,1).*tmp(2:end,1)';
        contri_all_std = horzcat(contri_all_std,contri_val);


        % constant intercept for empirical models
        emp_const = horzcat(emp_const,tmp(1,1)');
        emp_coef = horzcat(emp_coef,tmp(2:end,1)');

    end

    emp_const = horzcat(emp_const,nanmean(fann_gam.PFT));

    cd('/Volumes/RemiLBNL/project6/code4Remi/functions');
    axis1=subaxis(3,3,(nn-1)*3 + 3,'SpacingVert',0.06,'SpacingHoriz',0.15,'MR',0.1,'ML',0.1,'MarginTop',.05,'MarginBottom',.05); 

    groupC2 = [27,120,55;27,120,55;...
    178,24,43;178,24,43;178,24,43;178,24,43;178,24,43;...
    33,102,172;33,102,172;33,102,172;33,102,172;33,102,172;33,102,172;33,102,172;33,102,172]/255;
    
     hold on;

     x = 1:15;

    for i = 1:15
        boxchart(x(i)*ones(size(contri_all(:,i))), contri_all(:,i), 'BoxFaceColor', groupC2(i,:),...
            'BoxFaceAlpha', 0.3, 'WhiskerLineColor',groupC2(i,:), 'MarkerStyle' ,'none')
    end

    hline = refline([0 0]);
    hline.Color = 'k';

    set(gca, 'box', 'off');

    ylim([-22,10]);

    factor_names = {'LPC';'LMA';...
        'Tair';'PP';'PAR';'VPD';'SWC';...
        'soilC';'soilN';'CN';'PH';'sand';'silt';'bulkD';'CEC'};

    if nn == 1
    ylabel('\DeltafLNR (%)','FontSize',9);
    end

    h = gca;
    h.XTick = 1:length(factor_names);
    h.XTickLabel = factor_names;
    h.TickLabelInterpreter = 'none';
    h.FontSize = 8;

    text(0.6,-20.5,strcat(char(102+nn),{' '},title),'fontsize',10,'FontWeight','bold');

    camroll(-90);
    
    
    
end


set(fs10,'PaperPositionMode','auto');
print('/Volumes/RemiLBNL/project11_Vcmax_val/Figures/Figure_s10','-djpeg','-r600');



%%%%%% Fig. r1 Compare N deposition map and climate contribution %%%%%%%%%%%%%%%

fr1=figure('Name','Nitrogen deposition and climate influence','Units', 'centimeters','Color','white', 'Position', [2, 2, 17, 8], ...
    'OuterPosition', [2, 2, 17, 8]);

% nitrogen deposition map
filename = '/Users/xzluo/Dropbox (Personal)/Data/1860_1993_2050_NITROGEN_830/data/N-deposition1993.tif';

[A,R] = readgeoraster(filename);
A(A<100) = NaN; % check nan values

cd('/Volumes/RemiLBNL/project6/code4Remi/functions');
subaxis(1,2,1,'SpacingVert',0.05,'SpacingHoriz',0.1,'MR',0.05,'ML',0.05,'MarginTop',.02,'MarginBottom',.05);

mapdata=flip(A)./100;

lat=-88.125:3.75:88.125;
lon=-177.5:5:177.5;
    
cd('/Volumes/RemiLBNL/project6/code4Remi/functions/m_map');
m_proj('Miller','lon',[-180 180],'lat',[-60 90]); % map projections, range of lon and lat of your data
m_coast('linewidth',0.5,'color',[0.2 0.2 0.2]); % coast line settings
hold on;
m_pcolor(lon,lat,mapdata); % draw your map here
shading INTERP;   % can be flat or INTERP

colorscheme = [84,48,5;...
            140,81,10;...
            191,129,45;...
            223,194,125;...
            246,232,195;...
            199,234,229;...
            128,205,193;...
            53,151,143;...
            1,102,94;...
            0,60,48];
colorscheme = colorscheme./255;
% colorscheme = jet(12);       

colormap(gca,colorscheme);
caxis([0 40]); 
h=colorbar;%
set(h,'YTick',0:5:40);
ylabel(h,{'N deposition'; '(0.1 g N/m2/year)'},'FontSize',9);

cd('/Volumes/RemiLBNL/project6/code4Remi/functions/m_map');
m_grid('box','off','tickdir','in','fontsize',7);

m_text(-170,85,'a','fontsize',10, 'FontWeight','bold');


% climate contribution map
filename = '/Volumes/RemiLBNL/project11_Vcmax_val/data/fann_pca_gam_dec.csv';
fann_gam = readtable(filename);

tmp1 = abs(fann_gam.s_leaf1_ + fann_gam.s_leaf2_);
tmp2 = abs(fann_gam.s_clim1_ + fann_gam.s_clim2_);
tmp3 = abs(fann_gam.s_soil1_ + fann_gam.s_soil2_);

tmp_xdata = horzcat(gam_pcas,reshape(EULC_2D,[],1),reshape(vData.annVc,[],1));%

ind = any(isnan(tmp_xdata), 2);

contri = nan(360*720,3);


contri(~ind,2) = tmp1./(tmp1 + tmp2 + tmp3).*1; % leaf traits
contri(~ind,1) = tmp2./(tmp1 + tmp2 + tmp3).*1; % climate
contri(~ind,3) = tmp3./(tmp1 + tmp2 + tmp3).*1; % soil

%%% smooth the contribution map over space, cannot smooth later as the ind has no physical meaning


contri_2D = reshape(contri(:,1),360,720);
cd('/Volumes/RemiLBNL/project6/code4Remi/functions');
contri_2D= smooth2a(contri_2D,10,10); 
contri_2D(180,359)=-50;
contri_2D(180,360)=50;


lat=-89.75:0.5:89.75;
lon=-179.75:0.5:179.75;

cd('/Volumes/RemiLBNL/project6/code4Remi/functions');
subaxis(1,2,2,'SpacingVert',0.05,'SpacingHoriz',0.1,'MR',0.05,'ML',0.05,'MarginTop',.02,'MarginBottom',.05);

mapdata=flip(contri_2D).*100;
    
cd('/Volumes/RemiLBNL/project6/code4Remi/functions/m_map');
m_proj('Miller','lon',[-180 180],'lat',[-60 90]); % map projections, range of lon and lat of your data
m_coast('linewidth',0.5,'color',[0.2 0.2 0.2]); % coast line settings
hold on;
m_pcolor(lon,lat,mapdata); % draw your map here
shading INTERP;   % can be flat or INTERP

colorscheme = [84,48,5;...
            140,81,10;...
            191,129,45;...
            223,194,125;...
            246,232,195;...
            199,234,229;...
            128,205,193;...
            53,151,143;...
            1,102,94;...
            0,60,48];
colorscheme = colorscheme./255;
% colorscheme = jet(12);       

colormap(gca,colorscheme);
caxis([0 100]); 
h=colorbar;%
set(h,'YTick',0:10:100);

ylabel(h,{'Climate contribution'; 'to fLNR (%)'},'FontSize',9);

cd('/Volumes/RemiLBNL/project6/code4Remi/functions/m_map');
m_grid('box','off','tickdir','in','fontsize',7);

m_text(-170,85,'b','fontsize',10,'FontWeight','bold');

set(fr1,'PaperPositionMode','auto');
print('/Volumes/RemiLBNL/project11_Vcmax_val/Figures/Figure_r1','-djpeg','-r600');

%%%%%% Compare N deposition map and climate contribution %%%%%%%%%%%%%%%



%%%%%% Fig. R2 Overlap EBF and leaf phoporhus content %%%%%%%%%%%%%%%

fr2=figure('Name','EBF LNP content','Units', 'centimeters','Color','white', 'Position', [2, 2, 12, 12], ...
    'OuterPosition', [2, 2, 15, 15]);

ind = EULC_2D == 3;
global_LNP = vData.Pa;

global_LNP(~ind) = NaN;

tropical_LNP = global_LNP(120:239,:);

lat=-29.75:0.5:29.75;
lon=-179.75:0.5:179.75;

mapdata=flip(tropical_LNP);
    
cd('/Volumes/RemiLBNL/project6/code4Remi/functions/m_map');
m_proj('Miller','lon',[-180 180],'lat',[-30 30]); % map projections, range of lon and lat of your data
m_coast('linewidth',0.5,'color',[0.2 0.2 0.2]); % coast line settings
hold on;
m_pcolor(lon,lat,mapdata); % draw your map here
shading INTERP;   % can be flat or INTERP

colorscheme = [84,48,5;...
            140,81,10;...
            %191,129,45;...
            223,194,125;...
            %246,232,195;...
            %199,234,229;...
            128,205,193;...
            53,151,143;...
            %1,102,94;...
            0,60,48];
colorscheme = colorscheme./255;

colormap(gca,colorscheme);
caxis([0.02 0.12]); 
h=colorbar;%
set(h,'YTick',0.02:0.02:0.12);

ylabel(h,{'Leaf phosphorus content'; '(g/m2)'},'FontSize',9);

cd('/Volumes/RemiLBNL/project6/code4Remi/functions/m_map');
m_grid('box','on','tickdir','in','fontsize',7);

%m_text(-170,85,'(b)','fontsize',10);

set(fr2,'PaperPositionMode','auto');
print('/Volumes/RemiLBNL/project11_Vcmax_val/Figures/Figure_r2','-djpeg','-r600');

%%%%%% Fig. R2 Overlap EBF and leaf phoporhus content %%%%%%%%%%%%%%%




%%%% output data into a netCDF file
output_file = '/Volumes/RemiLBNL/project11_Vcmax_val/fLNR_Vcmax25_globe7.nc';


nccreate(output_file,'lat','Dimensions',{'lat' 360});
nccreate(output_file,'lon','Dimensions',{'lon' 720});

ncwrite(output_file,'lat',flip(-89.75:0.5:89.75));
ncwrite(output_file,'lon',-179.75:0.5:179.75);


% output fLNR
% fLNR_RF, fLNR_un, fLNR_RF_AMM18, fLNR_RF_CB20

nccreate(output_file,'fLNR_RF','Dimensions', {'lat', 360, 'lon', 720}, 'Format', 'netcdf4_classic');
nccreate(output_file,'fLNR_un','Dimensions', {'lat', 360, 'lon', 720}, 'Format', 'netcdf4_classic');
nccreate(output_file,'fLNR_RF_AMM18','Dimensions', {'lat', 360, 'lon', 720}, 'Format', 'netcdf4_classic');
nccreate(output_file,'fLNR_RF_CB20','Dimensions', {'lat', 360, 'lon', 720}, 'Format', 'netcdf4_classic');
nccreate(output_file,'fLNR_EM1','Dimensions', {'lat', 360, 'lon', 720}, 'Format', 'netcdf4_classic');
nccreate(output_file,'fLNR_EM2','Dimensions', {'lat', 360, 'lon', 720}, 'Format', 'netcdf4_classic');
nccreate(output_file,'fLNR_EM3','Dimensions', {'lat', 360, 'lon', 720}, 'Format', 'netcdf4_classic');
nccreate(output_file,'fLNR_EM4','Dimensions', {'lat', 360, 'lon', 720}, 'Format', 'netcdf4_classic');
nccreate(output_file,'fLNR_EM5','Dimensions', {'lat', 360, 'lon', 720}, 'Format', 'netcdf4_classic');
nccreate(output_file,'fLNR_EO','Dimensions', {'lat', 360, 'lon', 720}, 'Format', 'netcdf4_classic');
nccreate(output_file,'fLNR_LUNA','Dimensions', {'lat', 360, 'lon', 720}, 'Format', 'netcdf4_classic');


vData.fLNR_std(vData.fLNR_std > 50) = NaN;
ncwrite(output_file,'fLNR_RF',vData.fLNR);
ncwrite(output_file,'fLNR_un',vData.fLNR_std);
ncwrite(output_file,'fLNR_RF_AMM18',fLNR2);
ncwrite(output_file,'fLNR_RF_CB20',fLNR3);
ncwrite(output_file,'fLNR_EM1',squeeze(vData.Fall(:,:,1)));
ncwrite(output_file,'fLNR_EM2',squeeze(vData.Fall(:,:,2)));
ncwrite(output_file,'fLNR_EM3',squeeze(vData.Fall(:,:,3)));
ncwrite(output_file,'fLNR_EM4',squeeze(vData.Fall(:,:,4)));
ncwrite(output_file,'fLNR_EM5',squeeze(vData.Fall(:,:,5)));
ncwrite(output_file,'fLNR_EO',squeeze(vData.Fall(:,:,6)));
ncwrite(output_file,'fLNR_LUNA',squeeze(vData.Fall(:,:,7)));


% output Vcmax25
% Vcmax25_RF, Vcmax25_un, Vcmax25_EM1, EM2, EM3, EM4, EM5, EO, LUNA
nccreate(output_file,'Vcmax25_RF','Dimensions', {'lat', 360, 'lon', 720}, 'Format', 'netcdf4_classic');
nccreate(output_file,'Vcmax25_un','Dimensions', {'lat', 360, 'lon', 720}, 'Format', 'netcdf4_classic');
nccreate(output_file,'Vcmax25_EM1','Dimensions', {'lat', 360, 'lon', 720}, 'Format', 'netcdf4_classic');
nccreate(output_file,'Vcmax25_EM2','Dimensions', {'lat', 360, 'lon', 720}, 'Format', 'netcdf4_classic');
nccreate(output_file,'Vcmax25_EM3','Dimensions', {'lat', 360, 'lon', 720}, 'Format', 'netcdf4_classic');
nccreate(output_file,'Vcmax25_EM4','Dimensions', {'lat', 360, 'lon', 720}, 'Format', 'netcdf4_classic');
nccreate(output_file,'Vcmax25_EM5','Dimensions', {'lat', 360, 'lon', 720}, 'Format', 'netcdf4_classic');
nccreate(output_file,'Vcmax25_EO','Dimensions', {'lat', 360, 'lon', 720}, 'Format', 'netcdf4_classic');
nccreate(output_file,'Vcmax25_LUNA','Dimensions', {'lat', 360, 'lon', 720}, 'Format', 'netcdf4_classic');


ncwrite(output_file,'Vcmax25_RF',vData.annVc);
ncwrite(output_file,'Vcmax25_un',vData.annVC_std);
ncwrite(output_file,'Vcmax25_EM1',vData.Vn1);
ncwrite(output_file,'Vcmax25_EM2',vData.Vn2);
ncwrite(output_file,'Vcmax25_EM3',vData.Vn3);
ncwrite(output_file,'Vcmax25_EM4',vData.Vn4);
ncwrite(output_file,'Vcmax25_EM5',vData.Venv1);
ncwrite(output_file,'Vcmax25_EO',vData.Venv2);
ncwrite(output_file,'Vcmax25_LUNA',vData.Vne);

% write attributes
ncwriteatt(output_file,'/','Vcmax25 (maximum carboxylation rate at 25C degree)','(ug m-2 s-1)');
ncwriteatt(output_file,'/','fLNR (fraction leaf nitrogen allocated to RuBisCO','(%)');
ncwriteatt(output_file,'/','spatial resolution','0.5 degree globe grid');
ncwriteatt(output_file,'/','distributed by','X.Luo');
ncwriteatt(output_file,'/','reference','title: Global variation in the fraction of leaf nitrogen allocated to photosynthesis');

