% Semi-Empirical model for PM2.5 with optimization for the weighting
% coefficients for mass balance of elements in PM2.5
% by Rong Wang et al Oct 16 2019
% this code shows data for Beijing, Shanghai and Guangzhou as examples

tic
clear;

% %% define variables
global M Emi X
% parameters for machline learning (ML)
myoptions = optimset('Display','off','FunValCheck','on','algorithm','sqp','MaxFunEvals',2000);

%% load data
cityid=1; % 1 Beijing; 2 Shanghai;3 Guangzhou
cityname={'beijing'; 'shanghai'; 'guangzhou'};
load(strcat(cityname{cityid},'.dat'),'-mat');
% loading variable "X": 1 day ID (1 for January 1 2013); 2 month (1-12); 3 year (2013-2018);
% 4-9 surface (2 m height) meteorological data: 
% 4 air temperature, 5 precipitation, 6 wind speed, 
% 7 zonal and 8 meridional wind, 9 relative humidity
% 10-15 surface pollution data:  
% 10 carbon monoxide (CO), 11 ozone (O3), 12 sulfur dioxide (SO2), 13 nitrogen oxides (NOx)
% 14  particulate matter with a diameter of 2.5 - 10 um (PM10-2.5)
% 15 particulate matter with a diameter of less than 2.5 um (PM2.5)

    %% Perform PLS regression (To represent the impact of atmospheric circulation and gas-particle partitioning on PM2.5, 
    % our study accounts for the autocorrelation of independent variables. 
    % While the meteorological variables are correlated to each other, SO2 and NOx are both correlated with CO-informed emission (Fig S7). 
    % Recently, Brown and Caldeira explored the application of partial linear square regression (PLSR) to 
    % reduce the systematic error due to collinearity of independent variables when applying MLR or PCR (Brown and Caldeira, 2017).)
    
    M=normalize(X(:,4:9)); % surface (2 m height) meteorological data: 
    O=normalize(X(:,11)); % ozone (O3)
    Y4=log(X(:,10)); % carbon monoxide (CO)
    [b4,bint4,r4,rint4,stats4]=regress(Y4,[ones(d2,1) M]); % PLS for CO
    Emi=r4; % emission proxy
    fc=X(:,14)./(X(:,15)+X(:,14)); % observed fraction of particles that have a diameter of 2.5-10 micrometer in PM10 (f10)
    fcmax=max(fc,[],1); % upper limit for the fraction < 1
    gamamax=1/fcmax-0.0001;  % upper limit for natural PM2.5 fraction < 1
    % Main program for city-specific coefficients
    % alpha, bea and gama are found by machine learning so that variance in A is best explained by emission and meteorology
    % In brief, we test a hypothesis that, in the transfer of gases (G) to particle (P), the abundance of ([P]+q[G]) 
    % is independent on chemistry when q weights for mass balance. In an occlusive condition,
    % the transfer of sulfur dioxide (SO2) to sulfate (SO42-) informs q=1.5. For ambient PM2.5,
    % such q is determined for SO2 and nitrogen oxides (NOx) that the variance in ([P]+q?[G]) 
    % is best explained by emission and meteorology.
    [aopt, fval] = fmincon(@residu,[1.5 2 0],[],[],[],[],[0 0 0.01],[10 10 gamamax], [], myoptions); 
    
    %% Plot the dependence of rou_A^2 agains alpha and beta
    if cityid==1 || cityid==2
        ms=zeros(41,26);
        for i=0:40
            for j=0:25
                ms(i+1,j+1)=-residu([i/4 j/5 aopt(3)]);
            end
        end
    else
        ms=zeros(28,28);
        for i=-2:25
            for j=-2:25
                ms(i+3,j+3)=-residu([i/5 j/5 aopt(3)]);
            end
        end
    end
    mmax=max(max(ms,[],2),[],1);
    mmin=min(min(ms,[],2),[],1);
    mmax1=mmax-(mmax-mmin)/128;
    idx1=find(ms==mmax);
    idx2=find(ms>mmax1);
    ms(idx2)=mmax1;
    ms(idx1)=mmax;
    
    surf(ms);
    colormap(cmap)
    caxis([mmin mmax]);
    nx=size(ms,2); ny=size(ms,1);
    
    if cityid==1
        axis([1 26 1 41 0.77 0.85])
        set(gca,'xtick',1:5:26); set(gca,'ytick',1:8:41);
        set(gca,'ztick',0.77:0.02:0.85);
        xlabel('Belta 0 to 5','Fontname', 'Arial','FontSize',10);
        ylabel('Alpha 0 to 10','Fontname', 'Arial','FontSize',10);
    elseif cityid==2
        axis([1 26 1 41 0.71 0.83])
        set(gca,'xtick',1:5:26); set(gca,'ytick',1:8:41);
        set(gca,'ztick',0.71:0.03:0.83);
        xlabel('Belta 0 to 5','Fontname', 'Arial','FontSize',10);
        ylabel('Alpha 0 to 10','Fontname', 'Arial','FontSize',10);
    elseif cityid==3
        caxis([0.46 mmax]);
        axis([1 28 1 28 0.1 0.6])
        set(gca,'xtick',3:5:28); set(gca,'ytick',3:5:28);
        set(gca,'ztick',0.1:0.1:0.6);
        xlabel('Belta -0.6 to 5','Fontname', 'Arial','FontSize',10);
        ylabel('Alpha -0.6 to 10','Fontname', 'Arial','FontSize',10);
    end
    zlabel('R2_A','Fontname', 'Arial','FontSize',10);
    title(cityname{cityid},'Fontname', 'Arial','FontSize',10)





