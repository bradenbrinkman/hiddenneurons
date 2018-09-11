% This script produces either figure 7 or 8 in the supplementary info (SI),
% showing the typical strength of the effective synaptic connections as a
% function of the fraction of observed neurons for a sigmoidal
% nonlinearity.

% The parameters (and axis labels) as set in this file produce Figure 8. To
% produce Figure 7, change signtype to 'signed' (see below). Axis labels
% must be adjusted by hand.

% Data is retrieved from the folder Jdata; this script belongs in the
% parent folder of Jdata.

p = 0.2;
beta = 0.3;
weighttype = 'normal';
%weighttype = 'lognormal';
scalingtype = 'balanced';
%scalingtype = 'classical';
signtype = 'DL';
%signtype = 'signed';
networktype = 'ER';
mu = -1.0;
mustr = '-1.0';
J0 = 1.0;
J0str = '1.0';
amp = 2.0;

rgbcolors = {[111/255 0 118/255],[31/255 138/255 206/255],[181/255 0 0],[72/255 72/255 72/255]};

N = 1e3;
Nhidvec = [990 890 790 690 590 490 390 290 190 90 1];

J0vec = [0.25 0.5 0.75 1.0];
J0strvec = {'0.25','0.5','0.75','1.0'};

numtrials = 100;

if(strcmp(scalingtype,'balanced'))
    scale = 1/sqrt(p*N);
elseif(strcmp(scalingtype,'classical'))
    scale = 1/(p*N);
else
    'Error! scale not assigned!'
end

%%
seedrange = 1:10;
Jstats = struct('Jrravg',[],'Jrrvar',[],'Jrrskew',[],'Jrrkur',[],'Jeffavg',[],'Jeffvar',[],'Jeffskew',[],'Jeffkur',[]);
markers = {'o-', 's-', '*-', '^-'};
figure;
for J0ind=1:length(J0vec)
    J0 = J0vec(J0ind);
    J0strout = J0strvec{J0ind};
    if strcmp(signtype,'signed')
        varJeff_theory(J0ind,:) = 1 + exp(-2.0)*J0^2*Nhidvec/N + 1.5*exp(-4.0)*J0^4*Nhidvec.^2/N^2;
    end
    scalingtype = 'balanced';
    for seed = seedrange
        
        filenameinJrrvar = ['Jdata/Jrrvar_sigmoidnonlin_' weighttype signtype scalingtype networktype '_N' num2str(N) '_p' num2str(p) '_amp' num2str(amp) '_mu' mustr '_J0' J0strout '_seed' num2str(seed) '.txt'];
        filenameinJeffvar = ['Jdata/Jeffvar_sigmoidnonlin_' weighttype signtype scalingtype networktype '_N' num2str(N) '_p' num2str(p) '_amp' num2str(amp) '_mu' mustr '_J0' J0strout '_seed' num2str(seed) '.txt'];
        
        
        Jstats(seed).Jrrvar = load(filenameinJrrvar,'Delimiter',' ');
        Jstats(seed).Jeffvar = load(filenameinJeffvar,'Delimiter',' ');
        
        
    end
    
    Janalysis_avg = struct('Jrravg',[],'Jrrvar',[],'Jrrskew',[],'Jrrkur',[],'Jeffavg',[],'Jeffvar',[],'Jeffskew',[],'Jeffkur',[]);
    
    for seed = seedrange
        Janalysis_avg(seed).Jrrvar = mean(Jstats(seed).Jrrvar,2);
        Janalysis_avg(seed).Jeffvar = mean(Jstats(seed).Jeffvar,2);
    end
    
    
    scalingtype = 'balanced';
        
        hold on;
        
        Jeffavg = [Janalysis_avg(:).Jeffvar];
        Jeffavg = mean(Jeffavg,2);
        Jefferr = [Janalysis_avg(:).Jeffvar];
        Jefferr = sqrt(var(Jefferr,0,2)/length(seedrange));
        
        Jrravg = [Janalysis_avg(:).Jrrvar];
        Jrravg = mean(Jrravg,2);
        Jrrerr = [Janalysis_avg(:).Jrrvar];
        Jrrerr = sqrt(var(Jrrerr,0,2)/length(seedrange));
        
%         Jratio = Jeffavg./Jrravg;
%         Jratio_err = Jratio.*sqrt((Jefferr./Jeffavg).^2 + (Jrrerr./Jrravg).^2);

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
        
        errorbar((N-Nhidvec)/N,sqrt(Jratio-1),0.5*Jratio_err./sqrt(Jratio-1),markers{J0ind},'MarkerSize',20,'LineWidth',5,'Color',rgbcolors{1})
        
        xlim([0 1]);
        ylim([0.95-1 0.35]);
        set(gca,'FontSize',18);
        xlabel('$N_{\rm rec}/N$','Interpreter','LaTeX','FontSize',28)
        ylabel('$\sigma[\mathcal J_{r,r\prime}^{\rm eff}-\mathcal J_{r,r\prime}]/\sigma[\mathcal J_{r,r\prime}]$','Interpreter','LaTeX','FontSize',28)
        set(gca,'FontSize',28);
%         lg = legend('Strong $1/\sqrt{N}$ scaling','Weak $1/N$ scaling','Location','Northeast');
%         legend boxoff
%         set(findobj(lg,'type','text'),'FontSize',28,'Interpreter','LaTeX');
        axis square

end

%%

% Plot classical scaling J0 = 1.0 case

scalingtype = 'classical';
Jstats = struct('Jrravg',[],'Jrrvar',[],'Jrrskew',[],'Jrrkur',[],'Jeffavg',[],'Jeffvar',[],'Jeffskew',[],'Jeffkur',[]);
J0strout = '1.0';
J0 = 1.0;
    for seed = seedrange
        
        filenameinJrrvar = ['Jdata/Jrrvar_sigmoidnonlin_' weighttype signtype scalingtype networktype '_N' num2str(N) '_p' num2str(p) '_amp' num2str(amp) '_mu' mustr '_J0' J0strout '_seed' num2str(seed) '.txt'];
        filenameinJeffvar = ['Jdata/Jeffvar_sigmoidnonlin_' weighttype signtype scalingtype networktype '_N' num2str(N) '_p' num2str(p) '_amp' num2str(amp) '_mu' mustr '_J0' J0strout '_seed' num2str(seed) '.txt'];
        
        
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
        
%         Jratio = Jeffavg./Jrravg;
%         Jratio_err = Jratio.*(Jefferr./Jeffavg + Jrrerr./Jrravg);

%         Jeffavg = [Janalysis_avg(:).Jeffvar];
%         Jrravg = [Janalysis_avg(:).Jrrvar];
%         Jratio_dist = Jeffavg./Jrravg;
%         Jratio = mean(Jratio_dist,2);
%         Jratio_err = var(Jratio_dist,0,2);


        errorbar((N-Nhidvec)/N,sqrt(Jratio-1),0.5*Jratio_err./sqrt(Jratio-1),'p-','MarkerSize',20,'LineWidth',5,'Color',rgbcolors{4})
        
        
        xlim([0 1]);
        ylim([0.95-1 0.6]);
        set(gca,'FontSize',18);
        xlabel('$N_{\rm rec}/N$','Interpreter','LaTeX','FontSize',28)
        ylabel('$\sigma[\mathcal J_{r,r\prime}^{\rm eff}-\mathcal J_{r,r\prime}]/\sigma[\mathcal J_{r,r\prime}]$','Interpreter','LaTeX','FontSize',28)
        title('ER + Dale''s law, sigmoidal nonlinearity','Interpreter','LaTeX','FontSize',24)
        set(gca,'FontSize',28);
%         lg = legend('Strong $1/\sqrt{N}$ scaling','Weak $1/N$ scaling','Location','Northeast');
%         legend boxoff
%         set(findobj(lg,'type','text'),'FontSize',28,'Interpreter','LaTeX');
        axis square

% for J0ind=1:length(J0vec)
% plot((N-Nhidvec)/N,sqrt(varJeff_theory(J0ind,:)-1),'k--','LineWidth',5)
% end
% varJeff_theory_weak = 1 + exp(-2.0)*J0^2*Nhidvec/N/(p*N) + 1.5*exp(-4.0)*J0^4*Nhidvec.^2/N^2/(p*N)^2;
% plot((N-Nhidvec)/N,sqrt(varJeff_theory_weak-1),'k--','LineWidth',5)

lg = legend('$J_0 = 0.25$ (Strong)','$J_0 = 0.5$ (Strong)','$J_0 = 0.75$ (Strong)','$J_0 = 1.0$ (Strong)','$J_0 = 1.0$ (Weak)','Length-3 path calc.','Location','Northeast');
legend boxoff
set(findobj(lg,'type','text'),'FontSize',22,'Interpreter','LaTeX');
