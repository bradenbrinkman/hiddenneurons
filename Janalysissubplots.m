% This script reproduces Figure 5 in the main paper.

% Data to produce different plots is contained in the folder Jdata. To
% produce specific figures from the paper change the parameters below to
% match parameters given in Figure captions and text.

% preamble stuff

p = 0.2;
beta = 0.3;
weighttype = 'normal';
mu = -1.0;
mustr = '-1.0';
J0 = 1.0;
J0str = '1.0';

N = 1e3;
Nhidvec = [990 890 790 690 590 490 390 290 190 90 1];

J0vec = [0.25 0.5 0.75 1.0];
J0strvec = {'0.25','0.5','0.75','1.0'};

numtrials = 100;

Jstats = struct('Jrravg',[],'Jrrvar',[],'Jrrskew',[],'Jrrkur',[],'Jeffavg',[],'Jeffvar',[],'Jeffskew',[],'Jeffkur',[]);
markers = {'o-', 's-', '*-', '^-'};
rgbcolors = {[111/255 0 118/255],[31/255 138/255 206/255],[181/255 0 0],[72/255 72/255 72/255]};

seedrange = 1:10;

figure;

%%

% subplot 1 -- Erdos-Reyni mixed synapses

subplot(141)

scalingtype = 'balanced';
signtype = 'signed';
networktype = 'ER';

for J0ind=1:length(J0vec)
    J0 = J0vec(J0ind);
    J0strout = J0strvec{J0ind};
    if strcmp(signtype,'signed')
        varJeff_theory(J0ind,:) = 1 + exp(-2.0)*J0^2*Nhidvec/N + 3*exp(-4.0)*J0^4*Nhidvec.^2/N^2;
    else
        %varJeff_theory(J0ind,:) = 1 + exp(-2.0)*J0^2*Nhidvec/N; %may not be correct!
    end
    scalingtype = 'balanced';
    for seed = seedrange;
        
        filenameinJrrvar = ['Jdata/Jrrvar_' weighttype signtype scalingtype networktype '_N' num2str(N) '_p' num2str(p) '_mu' mustr '_J0' J0strout '_seed' num2str(seed) '.txt'];
        filenameinJeffvar = ['Jdata/Jeffvar_' weighttype signtype scalingtype networktype '_N' num2str(N) '_p' num2str(p) '_mu' mustr '_J0' J0strout '_seed' num2str(seed) '.txt'];
        
        
        Jstats(seed).Jrrvar = load(filenameinJrrvar,'Delimiter',' ');
        Jstats(seed).Jeffvar = load(filenameinJeffvar,'Delimiter',' ');
        
        
    end
    
    Janalysis_avg = struct('Jrravg',[],'Jrrvar',[],'Jrrskew',[],'Jrrkur',[],'Jeffavg',[],'Jeffvar',[],'Jeffskew',[],'Jeffkur',[]);
    
    for seed = seedrange
        Janalysis_avg(seed).Jrrvar = mean(Jstats(seed).Jrrvar,2);
        Janalysis_avg(seed).Jeffvar = mean(Jstats(seed).Jeffvar,2);
    end
    
    
    scalingtype == 'balanced';
        
        hold on;
        
        Jeffavg = [Janalysis_avg(:).Jeffvar];
        Jeffavg = mean(Jeffavg,2);
        Jefferr = [Janalysis_avg(:).Jeffvar];
        Jefferr = sqrt(var(Jefferr,0,2)/length(seedrange));
        
        Jrravg = [Janalysis_avg(:).Jrrvar];
        Jrravg = mean(Jrravg,2);
        Jrrerr = [Janalysis_avg(:).Jrrvar];
        Jrrerr = sqrt(var(Jrrerr,0,2)/length(seedrange));

        temp1 = [Janalysis_avg(:).Jeffvar];
        temp2 = [Janalysis_avg(:).Jrrvar];
        
        JeffJrrcov = zeros(size(Nhidvec'));
        Jeffvar = zeros(size(Nhidvec'));
        Jrrvar= zeros(size(Nhidvec'));
        for j=1:length(temp1(:,1))
            temp3 = cov(temp1(j,:),temp2(j,:));
            JeffJrrcov(j) = temp3(1,2);
            Jeffvar(j) = temp3(1,1);
            Jrrvar(j) = temp3(2,2);
        end
        
        Jratio = Jeffavg./Jrravg;
        Jratio_err = Jratio.*sqrt(Jeffvar./(Jeffavg.^2) + Jrrvar./(Jrravg.^2)-2*JeffJrrcov./Jeffavg./Jrravg)/sqrt(length(seedrange));
        
        p1(J0ind) = errorbar((N-Nhidvec)/N,sqrt(Jratio-1),0.5*Jratio_err./sqrt(Jratio-1),markers{J0ind},'MarkerSize',12,'LineWidth',5,'Color',rgbcolors{1});
        
%         xlim([0 1]);
%         ylim([0.95-1 0.5]);
%         set(gca,'FontSize',18);
%  %       xlabel('$N_{\rm rec}/N$','Interpreter','LaTeX','FontSize',20)
%  %       ylabel('$\sigma[\mathcal J_{r,r\prime}^{\rm eff}-\mathcal J_{r,r\prime}]/\sigma[\mathcal J_{r,r\prime}]$','Interpreter','LaTeX','FontSize',20)
%         set(gca,'FontSize',28);
%         axis square

end

scalingtype = 'classical';
J0strout = '1.0';
J0 = 1.0;



    for seed = seedrange;
        
        filenameinJrrvar = ['Jdata/Jrrvar_' weighttype signtype scalingtype networktype '_N' num2str(N) '_p' num2str(p) '_mu' mustr '_J0' J0strout '_seed' num2str(seed) '.txt'];
        filenameinJeffvar = ['Jdata/Jeffvar_' weighttype signtype scalingtype networktype '_N' num2str(N) '_p' num2str(p) '_mu' mustr '_J0' J0strout '_seed' num2str(seed) '.txt'];
        
        
        Jstats(seed).Jrrvar = load(filenameinJrrvar,'Delimiter',' ');
        Jstats(seed).Jeffvar = load(filenameinJeffvar,'Delimiter',' ');
        
        
    end
    
    Janalysis_avg = struct('Jrravg',[],'Jrrvar',[],'Jrrskew',[],'Jrrkur',[],'Jeffavg',[],'Jeffvar',[],'Jeffskew',[],'Jeffkur',[]);
    
    for seed = seedrange
        Janalysis_avg(seed).Jrrvar = mean(Jstats(seed).Jrrvar,2);
        Janalysis_avg(seed).Jeffvar = mean(Jstats(seed).Jeffvar,2);
    end

        Jeffavg = [Janalysis_avg(:).Jeffvar];
        Jeffavg = mean(Jeffavg,2);
        Jefferr = [Janalysis_avg(:).Jeffvar];
        Jefferr = sqrt(var(Jefferr,0,2)/length(seedrange));
        
        Jrravg = [Janalysis_avg(:).Jrrvar];
        Jrravg = mean(Jrravg,2);
        Jrrerr = [Janalysis_avg(:).Jrrvar];
        Jrrerr = sqrt(var(Jrrerr,0,2)/length(seedrange));
        
                temp1 = [Janalysis_avg(:).Jeffvar];
        temp2 = [Janalysis_avg(:).Jrrvar];
        
        JeffJrrcov = zeros(size(Nhidvec'));
        Jeffvar = zeros(size(Nhidvec'));
        Jrrvar= zeros(size(Nhidvec'));
        for j=1:length(temp1(:,1))
            temp3 = cov(temp1(j,:),temp2(j,:));
            JeffJrrcov(j) = temp3(1,2);
            Jeffvar(j) = temp3(1,1);
            Jrrvar(j) = temp3(2,2);
        end
        
        Jratio = Jeffavg./Jrravg;
        Jratio_err = Jratio.*sqrt(Jeffvar./(Jeffavg.^2) + Jrrvar./(Jrravg.^2)-2*JeffJrrcov./Jeffavg./Jrravg)/sqrt(length(seedrange));


        p1w = errorbar((N-Nhidvec)/N,sqrt(Jratio-1),0.5*Jratio_err./sqrt(Jratio-1),'p-','MarkerSize',12,'LineWidth',5,'Color',rgbcolors{4})
        
        
        xlim([0 1]);
        ylim([0.95-1 0.55]);
        set(gca,'FontSize',18);
        xlabel('$N_{\rm rec}/N$','Interpreter','LaTeX','FontSize',18)
        ylabel('$\sigma[\mathcal J_{r,r\prime}^{\rm eff}-\mathcal J_{r,r\prime}]/\sigma[\mathcal J_{r,r\prime}]$','Interpreter','LaTeX','FontSize',18)
        title('ER + mixed synapses','Interpreter','LaTeX','FontSize',16);
        set(gca,'FontSize',18);
        axis square
        box on;

for J0ind=1:length(J0vec)
p1t = plot((N-Nhidvec)/N,sqrt(varJeff_theory(J0ind,:)-1),'k--','LineWidth',5)
end
varJeff_theory_weak = 1 + exp(-2.0)*J0^2*Nhidvec/N/(p*N) + 3*exp(-4.0)*J0^4*Nhidvec.^2/N^2/(p*N)^2;
plot((N-Nhidvec)/N,sqrt(varJeff_theory_weak-1),'k--','LineWidth',5)

sh = subplot(144);
shp = get(sh,'position');
axis off
lg = legend(sh,[p1(:);p1w;p1t],'$J_0 = 0.25$ (Strong)','$J_0 = 0.5$ (Strong)','$J_0 = 0.75$ (Strong)','$J_0 = 1.0$ (Strong)','$J_0 = 1.0$ (Weak)','Length-3 path calc.','Location','Northeast');
%lg = legend('$J_0 = 0.25$ (Strong)','$J_0 = 0.5$ (Strong)','$J_0 = 0.75$ (Strong)','$J_0 = 1.0$ (Strong)','$J_0 = 1.0$ (Weak)','Length-3 path calc.','Location','Northeast');
legend boxoff
set(lg,'position',shp);
set(findobj(lg,'type','text'),'FontSize',14,'Interpreter','LaTeX');




%%

% subplot 2 -- Erdos-Reyni Dale's Law

subplot(142)

scalingtype = 'balanced';
signtype = 'DL';
networktype = 'ER';

for J0ind=1:length(J0vec)
    J0 = J0vec(J0ind);
    J0strout = J0strvec{J0ind};
    if strcmp(signtype,'signed')
        varJeff_theory(J0ind,:) = 1 + exp(-2.0)*J0^2*Nhidvec/N + 3*exp(-4.0)*J0^4*Nhidvec.^2/N^2;
    else
        %varJeff_theory(J0ind,:) = 1 + exp(-2.0)*J0^2*Nhidvec/N; %may not be correct!
    end
    scalingtype = 'balanced';
    for seed = seedrange;
        
        filenameinJrrvar = ['Jdata/Jrrvar_' weighttype signtype scalingtype networktype '_N' num2str(N) '_p' num2str(p) '_mu' mustr '_J0' J0strout '_seed' num2str(seed) '.txt'];
        filenameinJeffvar = ['Jdata/Jeffvar_' weighttype signtype scalingtype networktype '_N' num2str(N) '_p' num2str(p) '_mu' mustr '_J0' J0strout '_seed' num2str(seed) '.txt'];
        
        
        Jstats(seed).Jrrvar = load(filenameinJrrvar,'Delimiter',' ');
        Jstats(seed).Jeffvar = load(filenameinJeffvar,'Delimiter',' ');
        
        
    end
    
    Janalysis_avg = struct('Jrravg',[],'Jrrvar',[],'Jrrskew',[],'Jrrkur',[],'Jeffavg',[],'Jeffvar',[],'Jeffskew',[],'Jeffkur',[]);
    
    for seed = seedrange
        Janalysis_avg(seed).Jrrvar = mean(Jstats(seed).Jrrvar,2);
        Janalysis_avg(seed).Jeffvar = mean(Jstats(seed).Jeffvar,2);
    end
    
    
    scalingtype == 'balanced';
        
        hold on;
        
        Jeffavg = [Janalysis_avg(:).Jeffvar];
        Jeffavg = mean(Jeffavg,2);
        Jefferr = [Janalysis_avg(:).Jeffvar];
        Jefferr = sqrt(var(Jefferr,0,2)/length(seedrange));
        
        Jrravg = [Janalysis_avg(:).Jrrvar];
        Jrravg = mean(Jrravg,2);
        Jrrerr = [Janalysis_avg(:).Jrrvar];
        Jrrerr = sqrt(var(Jrrerr,0,2)/length(seedrange));

        temp1 = [Janalysis_avg(:).Jeffvar];
        temp2 = [Janalysis_avg(:).Jrrvar];
        
        JeffJrrcov = zeros(size(Nhidvec'));
        Jeffvar = zeros(size(Nhidvec'));
        Jrrvar= zeros(size(Nhidvec'));
        for j=1:length(temp1(:,1))
            temp3 = cov(temp1(j,:),temp2(j,:));
            JeffJrrcov(j) = temp3(1,2);
            Jeffvar(j) = temp3(1,1);
            Jrrvar(j) = temp3(2,2);
        end
        
        Jratio = Jeffavg./Jrravg;
        Jratio_err = Jratio.*sqrt(Jeffvar./(Jeffavg.^2) + Jrrvar./(Jrravg.^2)-2*JeffJrrcov./Jeffavg./Jrravg)/sqrt(length(seedrange));
        
        p2 = errorbar((N-Nhidvec)/N,sqrt(Jratio-1),0.5*Jratio_err./sqrt(Jratio-1),markers{J0ind},'MarkerSize',12,'LineWidth',5,'Color',rgbcolors{1});
        
%         xlim([0 1]);
%         ylim([0.95-1 0.35]);
%         set(gca,'FontSize',18);
%  %       xlabel('$N_{\rm rec}/N$','Interpreter','LaTeX','FontSize',18)
%  %       ylabel('$\sigma[\mathcal J_{r,r\prime}^{\rm eff}-\mathcal J_{r,r\prime}]/\sigma[\mathcal J_{r,r\prime}]$','Interpreter','LaTeX','FontSize',28)
%         set(gca,'FontSize',28);
%         axis square

end

scalingtype = 'classical';
J0strout = '1.0';
J0 = 1.0;
    for seed = seedrange;
        
        filenameinJrrvar = ['Jdata/Jrrvar_' weighttype signtype scalingtype networktype '_N' num2str(N) '_p' num2str(p) '_mu' mustr '_J0' J0strout '_seed' num2str(seed) '.txt'];
        filenameinJeffvar = ['Jdata/Jeffvar_' weighttype signtype scalingtype networktype '_N' num2str(N) '_p' num2str(p) '_mu' mustr '_J0' J0strout '_seed' num2str(seed) '.txt'];
        
        
        Jstats(seed).Jrrvar = load(filenameinJrrvar,'Delimiter',' ');
        Jstats(seed).Jeffvar = load(filenameinJeffvar,'Delimiter',' ');
        
        
    end
    
    Janalysis_avg = struct('Jrravg',[],'Jrrvar',[],'Jrrskew',[],'Jrrkur',[],'Jeffavg',[],'Jeffvar',[],'Jeffskew',[],'Jeffkur',[]);
    
    for seed = seedrange
        Janalysis_avg(seed).Jrrvar = mean(Jstats(seed).Jrrvar,2);
        Janalysis_avg(seed).Jeffvar = mean(Jstats(seed).Jeffvar,2);
    end

        Jeffavg = [Janalysis_avg(:).Jeffvar];
        Jeffavg = mean(Jeffavg,2);
        Jefferr = [Janalysis_avg(:).Jeffvar];
        Jefferr = sqrt(var(Jefferr,0,2)/length(seedrange));
        
        Jrravg = [Janalysis_avg(:).Jrrvar];
        Jrravg = mean(Jrravg,2);
        Jrrerr = [Janalysis_avg(:).Jrrvar];
        Jrrerr = sqrt(var(Jrrerr,0,2)/length(seedrange));
        
                temp1 = [Janalysis_avg(:).Jeffvar];
        temp2 = [Janalysis_avg(:).Jrrvar];
        
        JeffJrrcov = zeros(size(Nhidvec'));
        Jeffvar = zeros(size(Nhidvec'));
        Jrrvar= zeros(size(Nhidvec'));
        for j=1:length(temp1(:,1))
            temp3 = cov(temp1(j,:),temp2(j,:));
            JeffJrrcov(j) = temp3(1,2);
            Jeffvar(j) = temp3(1,1);
            Jrrvar(j) = temp3(2,2);
        end
        
        Jratio = Jeffavg./Jrravg;
        Jratio_err = Jratio.*sqrt(Jeffvar./(Jeffavg.^2) + Jrrvar./(Jrravg.^2)-2*JeffJrrcov./Jeffavg./Jrravg)/sqrt(length(seedrange));


        errorbar((N-Nhidvec)/N,sqrt(Jratio-1),0.5*Jratio_err./sqrt(Jratio-1),'p-','MarkerSize',12,'LineWidth',5,'Color',rgbcolors{4});
        
        
        xlim([0 1]);
        ylim([0.95-1 0.55]);
        set(gca,'FontSize',18);
        xlabel('$N_{\rm rec}/N$','Interpreter','LaTeX','FontSize',18)
        title('ER + Dale''s law','Interpreter','LaTeX','FontSize',16);
 %       ylabel('$\sigma[\mathcal J_{r,r\prime}^{\rm eff}-\mathcal J_{r,r\prime}]/\sigma[\mathcal J_{r,r\prime}]$','Interpreter','LaTeX','FontSize',28)
        set(gca,'FontSize',18);
        axis square
        box on;

%%

% subplot 3 -- Watts-Strogatz mixed synapses

subplot(143)

scalingtype = 'balanced';
signtype = 'signed';
networktype = 'WS';

for J0ind=1:length(J0vec)
    J0 = J0vec(J0ind);
    J0strout = J0strvec{J0ind};
    varJeff_theory(J0ind,:) = 1 + exp(-2.0)*J0^2*Nhidvec/N + 1.5*exp(-4.0)*J0^4*Nhidvec.^2/N^2;
    scalingtype = 'balanced';
    for seed = seedrange;
        
        filenameinJrrvar = ['Jdata/Jrrvar_' weighttype signtype scalingtype networktype '_N' num2str(N) '_p' num2str(p) '_beta' num2str(beta) '_mu' mustr '_J0' J0strout '_seed' num2str(seed) '.txt'];
        filenameinJeffvar = ['Jdata/Jeffvar_' weighttype signtype scalingtype networktype '_N' num2str(N) '_p' num2str(p) '_beta' num2str(beta) '_mu' mustr '_J0' J0strout '_seed' num2str(seed) '.txt'];
        
        
        Jstats(seed).Jrrvar = load(filenameinJrrvar,'Delimiter',' ');
        Jstats(seed).Jeffvar = load(filenameinJeffvar,'Delimiter',' ');
        
        
    end
    
    Janalysis_avg = struct('Jrravg',[],'Jrrvar',[],'Jrrskew',[],'Jrrkur',[],'Jeffavg',[],'Jeffvar',[],'Jeffskew',[],'Jeffkur',[]);
    
    for seed = seedrange
        Janalysis_avg(seed).Jrrvar = mean(Jstats(seed).Jrrvar,2);
        Janalysis_avg(seed).Jeffvar = mean(Jstats(seed).Jeffvar,2);
    end
    
    
    if(strcmp(scalingtype,'balanced'))
        
        hold on;
        
        Jeffavg = [Janalysis_avg(:).Jeffvar];
        Jeffavg = mean(Jeffavg,2);
        Jefferr = [Janalysis_avg(:).Jeffvar];
        Jefferr = sqrt(var(Jefferr,0,2)/length(seedrange));
        
        Jrravg = [Janalysis_avg(:).Jrrvar];
        Jrravg = mean(Jrravg,2);
        Jrrerr = [Janalysis_avg(:).Jrrvar];
        Jrrerr = sqrt(var(Jrrerr,0,2)/length(seedrange));

                temp1 = [Janalysis_avg(:).Jeffvar];
        temp2 = [Janalysis_avg(:).Jrrvar];
        
        JeffJrrcov = zeros(size(Nhidvec'));
        Jeffvar = zeros(size(Nhidvec'));
        Jrrvar= zeros(size(Nhidvec'));
        for j=1:length(temp1(:,1))
            temp3 = cov(temp1(j,:),temp2(j,:));
            JeffJrrcov(j) = temp3(1,2);
            Jeffvar(j) = temp3(1,1);
            Jrrvar(j) = temp3(2,2);
        end
        
        Jratio = Jeffavg./Jrravg;
        Jratio_err = Jratio.*sqrt(Jeffvar./(Jeffavg.^2) + Jrrvar./(Jrravg.^2)-2*JeffJrrcov./Jeffavg./Jrravg)/sqrt(length(seedrange));
        
        p3 = errorbar((N-Nhidvec)/N,sqrt(Jratio-1),0.5*Jratio_err./sqrt(Jratio-1),markers{J0ind},'MarkerSize',12,'LineWidth',5,'Color',rgbcolors{1});
        
%         xlim([0 1]);
%         ylim([0.95-1 0.6]);
%         set(gca,'FontSize',18);
%  %       xlabel('$N_{\rm rec}/N$','Interpreter','LaTeX','FontSize',18)
% %        ylabel('$\sigma[\mathcal J_{r,r\prime}^{\rm eff}-\mathcal J_{r,r\prime}]/\sigma[\mathcal J_{r,r\prime}]$','Interpreter','LaTeX','FontSize',28)
%         set(gca,'FontSize',28);
%         axis square
     end
end

% Plot classical scaling J0 = 1.0 case

scalingtype = 'classical';
Jstats = struct('Jrravg',[],'Jrrvar',[],'Jrrskew',[],'Jrrkur',[],'Jeffavg',[],'Jeffvar',[],'Jeffskew',[],'Jeffkur',[]);
J0strout = '1.0';
J0 = 1.0;
    for seed = seedrange;
        
        filenameinJrrvar = ['Jdata/Jrrvar_' weighttype signtype scalingtype networktype '_N' num2str(N) '_p' num2str(p) '_beta' num2str(beta) '_mu' mustr '_J0' J0strout '_seed' num2str(seed) '.txt'];
        filenameinJeffvar = ['Jdata/Jeffvar_' weighttype signtype scalingtype networktype '_N' num2str(N) '_p' num2str(p) '_beta' num2str(beta) '_mu' mustr '_J0' J0strout '_seed' num2str(seed) '.txt'];
        
        
        Jstats(seed).Jrrvar = load(filenameinJrrvar,'Delimiter',' ');
        Jstats(seed).Jeffvar = load(filenameinJeffvar,'Delimiter',' ');
        
        
    end
    
    Janalysis_avg = struct('Jrravg',[],'Jrrvar',[],'Jrrskew',[],'Jrrkur',[],'Jeffavg',[],'Jeffvar',[],'Jeffskew',[],'Jeffkur',[]);
    
    for seed = seedrange
        Janalysis_avg(seed).Jrrvar = mean(Jstats(seed).Jrrvar,2);
        Janalysis_avg(seed).Jeffvar = mean(Jstats(seed).Jeffvar,2);
    end

        Jeffavg = [Janalysis_avg(:).Jeffvar];
        Jeffavg = mean(Jeffavg,2);
        Jefferr = [Janalysis_avg(:).Jeffvar];
        Jefferr = sqrt(var(Jefferr,0,2)/length(seedrange));
        
        Jrravg = [Janalysis_avg(:).Jrrvar];
        Jrravg = mean(Jrravg,2);
        Jrrerr = [Janalysis_avg(:).Jrrvar];
        Jrrerr = sqrt(var(Jrrerr,0,2)/length(seedrange));

                temp1 = [Janalysis_avg(:).Jeffvar];
        temp2 = [Janalysis_avg(:).Jrrvar];
        
        JeffJrrcov = zeros(size(Nhidvec'));
        Jeffvar = zeros(size(Nhidvec'));
        Jrrvar= zeros(size(Nhidvec'));
        for j=1:length(temp1(:,1))
            temp3 = cov(temp1(j,:),temp2(j,:));
            JeffJrrcov(j) = temp3(1,2);
            Jeffvar(j) = temp3(1,1);
            Jrrvar(j) = temp3(2,2);
        end
        
        Jratio = Jeffavg./Jrravg;
        Jratio_err = Jratio.*sqrt(Jeffvar./(Jeffavg.^2) + Jrrvar./(Jrravg.^2)-2*JeffJrrcov./Jeffavg./Jrravg)/sqrt(length(seedrange));
        
       

        errorbar((N-Nhidvec)/N,sqrt(Jratio-1),0.5*Jratio_err./sqrt(Jratio-1),'p-','MarkerSize',12,'LineWidth',5,'Color',rgbcolors{4})
        
        
        xlim([0 1]);
        ylim([0.95-1 0.55]);
        set(gca,'FontSize',18);
        xlabel('$N_{\rm rec}/N$','Interpreter','LaTeX','FontSize',18)
 %       ylabel('$\sigma[\mathcal J_{r,r\prime}^{\rm eff}-\mathcal J_{r,r\prime}]/\sigma[\mathcal J_{r,r\prime}]$','Interpreter','LaTeX','FontSize',28)
        set(gca,'FontSize',18);
        axis square
        box on;

title('WS + mixed synapses','Interpreter','LaTeX','FontSize',16);

