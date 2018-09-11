% This file produces Figures 10, 11, & 12 in the Supplementary Information (SI). 

% The first figure produced demonstrates the agreement between mean field 
% predictions of the firing rates with empirical estimates from full 
% simulations. These data use the entire network of neurons. This forms
% Figure 10 in the SI.

% The second figure produced compares empirically estimated firing rates
% from full simulations with predictions based on: top row) mean field treatment
% of the hidden neurons only (i.e., recorded neurons removed from network)
% and bottom row) prediction of firing rates accounting for first order correction
% from activity of recorded neurons. For the variable choice Nhid = 900
% (set below) this plot corresponds to Figure 11. If Nhid = 500, this plot
% corresponds to Figure 12.

%See SI for full explanations. 

% Data files necessary to produce plots are included in this folder

p = 0.2;
beta = 0.3;
weighttype = 'normal';
%weighttype = 'lognormal';
scalingtype = 'balanced';
%scalingtype = 'classical';
%signtype = 'DL';
signtype = 'signed';
networktype = 'ER';
mu = -1.0;
mustr = '-1.0';
J0 = 1.0;
J0str = '1.0';
seed = 1;

rgbcolors = {[111/255 0 118/255],[31/255 138/255 206/255],[181/255 0 0],[72/255 72/255 72/255]};

N = 1e3;

Nhid = 900; %change to 500 to reproduce Figure S3.

J0vec = [0.25 0.5 0.75 1.0];
J0strvec = {'0.25','0.5','0.75','1.0'};

%%

Avgratedata = struct('approx',[],'emp',[],'MFThid',[],'MFTfull',[]);


for J0ind=1:length(J0vec)
    J0 = J0vec(J0ind);
    J0strout = J0strvec{J0ind};
        
        filenameinAvgapprox = ['SFigs/Avgrates_approx_hiddenlinearrates_' weighttype signtype scalingtype networktype '_N' num2str(N) '_Nhid' num2str(Nhid) '_p' num2str(p) '_mu' mustr '_J0' J0strout '_seed' num2str(seed) '.txt'];
        filenameinAvgempirical = ['SFigs/Avgrates_empirical_hiddenlinearrates_' weighttype signtype scalingtype networktype '_N' num2str(N) '_Nhid' num2str(Nhid) '_p' num2str(p) '_mu' mustr '_J0' J0strout '_seed' num2str(seed) '.txt'];
        filenameinAvgMFThid = ['SFigs/Avgrates_MFT_hiddenlinearrates_' weighttype signtype scalingtype networktype '_N' num2str(N) '_Nhid' num2str(Nhid) '_p' num2str(p) '_mu' mustr '_J0' J0strout '_seed' num2str(seed) '.txt'];
        filenameinAvgMFTfull = ['SFigs/Avgrates_MFTfull_hiddenlinearrates_' weighttype signtype scalingtype networktype '_N' num2str(N) '_Nhid' num2str(Nhid) '_p' num2str(p) '_mu' mustr '_J0' J0strout '_seed' num2str(seed) '.txt'];

        
        Avgratedata(J0ind).approx = load(filenameinAvgapprox,'Delimiter',' ');
        Avgratedata(J0ind).emp = load(filenameinAvgempirical,'Delimiter',' ');
        Avgratedata(J0ind).MFThid = load(filenameinAvgMFThid,'Delimiter',' ');
        Avgratedata(J0ind).MFTfull = load(filenameinAvgMFTfull,'Delimiter',' ');
        
end
        
        %Figure to establish validity of mean field theory
        figure;
        subplot(141)
        plot(0:0.1:1.5,0:0.1:1.5,'k--','Linewidth',2)
        hold on;
        scatter(Avgratedata(1).emp,Avgratedata(1).MFTfull(1:Nhid),'o') %J0 = 0.25
        xlim([0 1.5]);
        ylim([0 1.5]);
        set(gca,'FontSize',14);
        xlabel('$\langle \dot{n}_h\rangle^{\rm emp}$','Interpreter','LaTeX','FontSize',16)
        ylabel('$\langle \dot{n}_h\rangle^{\rm full~MFT}$','Interpreter','LaTeX','FontSize',16)
        title('$J_0 = 0.25$','Interpreter','LaTeX','FontSize',16)
        set(gca,'FontSize',16);
        axis square
        
        subplot(142) %J0 = 0.5
        plot(0:0.1:1.5,0:0.1:1.5,'k--','Linewidth',2)
        hold on;
        scatter(Avgratedata(2).emp,Avgratedata(2).MFTfull(1:Nhid),'o') %J0 = 0.25
        xlim([0 1.5]);
        ylim([0 1.5]);
        set(gca,'FontSize',14);
        xlabel('$\langle \dot{n}_h\rangle^{\rm emp}$','Interpreter','LaTeX','FontSize',16)
        ylabel('$\langle \dot{n}_h\rangle^{\rm full~MFT}$','Interpreter','LaTeX','FontSize',16)
        title('$J_0 = 0.5$','Interpreter','LaTeX','FontSize',16)
        set(gca,'FontSize',16);
        axis square
        
        subplot(143) %J0 = 0.75
        plot(0:0.1:1.5,0:0.1:1.5,'k--','Linewidth',2)
        hold on;
        scatter(Avgratedata(3).emp,Avgratedata(3).MFTfull(1:Nhid),'o') %J0 = 0.25
        xlim([0 1.5]);
        ylim([0 1.5]);
        set(gca,'FontSize',14);
        xlabel('$\langle \dot{n}_h\rangle^{\rm emp}$','Interpreter','LaTeX','FontSize',16)
        ylabel('$\langle \dot{n}_h\rangle^{\rm full~MFT}$','Interpreter','LaTeX','FontSize',16)
        title('$J_0 = 0.75$','Interpreter','LaTeX','FontSize',16)
        set(gca,'FontSize',16);
        axis square
        
        subplot(144) %J0 = 1.0
        plot(0:0.1:1.5,0:0.1:1.5,'k--','Linewidth',2)
        hold on;
        scatter(Avgratedata(4).emp,Avgratedata(4).MFTfull(1:Nhid),'o') %J0 = 0.25
        xlim([0 1.5]);
        ylim([0 1.5]);
        set(gca,'FontSize',14);
        xlabel('$\langle \dot{n}_h\rangle^{\rm emp}$','Interpreter','LaTeX','FontSize',16)
        ylabel('$\langle \dot{n}_h\rangle^{\rm full~MFT}$','Interpreter','LaTeX','FontSize',16)
        title('$J_0 = 1.0$','Interpreter','LaTeX','FontSize',16)
        set(gca,'FontSize',16);
        axis square
        
        %% Testing linear approximation of hidden unit firing rates
        
        figure;
        subplot(241) %J0 = 0.25, subMFT vs emp
        plot(0:0.1:1.5,0:0.1:1.5,'k--','Linewidth',2)
        hold on;
        scatter(Avgratedata(1).emp,Avgratedata(1).MFThid,'o') %J0 = 0.25
        xlim([0 1.5]);
        ylim([0 1.5]);
        set(gca,'FontSize',14);
        xlabel('$\langle \dot{n}_h\rangle^{\rm emp}$','Interpreter','LaTeX','FontSize',16)
        ylabel('$\nu_h$','Interpreter','LaTeX','FontSize',16)
        title('$J_0 = 0.25$','Interpreter','LaTeX','FontSize',16)
        set(gca,'FontSize',16);
        axis square
        
        subplot(245) %J0 = 0.25, subMFT vs emp
        plot(0:0.1:1.5,0:0.1:1.5,'k--','Linewidth',2)
        hold on;
        scatter(Avgratedata(1).emp,Avgratedata(1).approx,'o') %J0 = 0.25
        xlim([0 1.5]);
        ylim([0 1.5]);
        set(gca,'FontSize',14);
        xlabel('$\langle \dot{n}_h\rangle^{\rm emp}$','Interpreter','LaTeX','FontSize',16)
        ylabel('$\langle \dot{n}_h\rangle^{\rm approx}$','Interpreter','LaTeX','FontSize',16)
        title('$J_0 = 0.25$','Interpreter','LaTeX','FontSize',16)
        set(gca,'FontSize',16);
        axis square
        
        subplot(242) %J0 = 0.5
        plot(0:0.1:1.5,0:0.1:1.5,'k--','Linewidth',2)
        hold on;
        scatter(Avgratedata(2).emp,Avgratedata(2).MFThid,'o') %J0 = 0.25
        xlim([0 1.5]);
        ylim([0 1.5]);
        set(gca,'FontSize',14);
        xlabel('$\langle \dot{n}_h\rangle^{\rm emp}$','Interpreter','LaTeX','FontSize',16)
        ylabel('$\nu_h$','Interpreter','LaTeX','FontSize',16)
        title('$J_0 = 0.5$','Interpreter','LaTeX','FontSize',16)
        set(gca,'FontSize',16);
        axis square
        
        subplot(246) %J0 = 0.5
        plot(0:0.1:1.5,0:0.1:1.5,'k--','Linewidth',2)
        hold on;
        scatter(Avgratedata(2).emp,Avgratedata(2).approx,'o') %J0 = 0.25
        xlim([0 1.5]);
        ylim([0 1.5]);
        set(gca,'FontSize',14);
        xlabel('$\langle \dot{n}_h\rangle^{\rm emp}$','Interpreter','LaTeX','FontSize',16)
        ylabel('$\langle \dot{n}_h\rangle^{\rm approx}$','Interpreter','LaTeX','FontSize',16)
        title('$J_0 = 0.5$','Interpreter','LaTeX','FontSize',16)
        set(gca,'FontSize',16);
        axis square
        
        
        subplot(243) %J0 = 0.75
        plot(0:0.1:1.5,0:0.1:1.5,'k--','Linewidth',2)
        hold on;
        scatter(Avgratedata(3).emp,Avgratedata(3).MFThid,'o') %J0 = 0.25
        xlim([0 1.5]);
        ylim([0 1.5]);
        set(gca,'FontSize',14);
        xlabel('$\langle \dot{n}_h\rangle^{\rm emp}$','Interpreter','LaTeX','FontSize',16)
        ylabel('$\nu_h$','Interpreter','LaTeX','FontSize',16)
        title('$J_0 = 0.75$','Interpreter','LaTeX','FontSize',16)
        set(gca,'FontSize',16);
        axis square
        
        subplot(247) %J0 = 0.75
        plot(0:0.1:1.5,0:0.1:1.5,'k--','Linewidth',2)
        hold on;
        scatter(Avgratedata(3).emp,Avgratedata(3).approx,'o') %J0 = 0.25
        xlim([0 1.5]);
        ylim([0 1.5]);
        set(gca,'FontSize',14);
        xlabel('$\langle \dot{n}_h\rangle^{\rm emp}$','Interpreter','LaTeX','FontSize',16)
        ylabel('$\langle \dot{n}_h\rangle^{\rm approx}$','Interpreter','LaTeX','FontSize',16)
        title('$J_0 = 0.75$','Interpreter','LaTeX','FontSize',16)
        set(gca,'FontSize',16);
        axis square
        
        subplot(244) %J0 = 1.0
        plot(0:0.1:1.5,0:0.1:1.5,'k--','Linewidth',2)
        hold on;
        scatter(Avgratedata(4).emp,Avgratedata(4).MFThid,'o') %J0 = 0.25
        xlim([0 1.5]);
        ylim([0 1.5]);
        set(gca,'FontSize',14);
        xlabel('$\langle \dot{n}_h\rangle^{\rm emp}$','Interpreter','LaTeX','FontSize',16)
        ylabel('$\nu_h$','Interpreter','LaTeX','FontSize',16)
        title('$J_0 = 1.0$','Interpreter','LaTeX','FontSize',16)
        set(gca,'FontSize',16);
        axis square
        
        subplot(248) %J0 = 1.0       
        plot(0:0.1:1.5,0:0.1:1.5,'k--','Linewidth',2)
        hold on;
        scatter(Avgratedata(4).emp,Avgratedata(4).approx,'o')
        xlim([0 1.5]);
        ylim([0 1.5]);
        set(gca,'FontSize',14);
        xlabel('$\langle \dot{n}_h\rangle^{\rm emp}$','Interpreter','LaTeX','FontSize',16)
        ylabel('$\langle \dot{n}_h\rangle^{\rm approx}$','Interpreter','LaTeX','FontSize',16)
        title('$J_0 = 1.0$','Interpreter','LaTeX','FontSize',16)
        set(gca,'FontSize',16);
        axis square
        
