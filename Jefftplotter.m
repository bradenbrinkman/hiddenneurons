% This script reproduces Figure 6 from the main text. Necessary data files
% are included in the same folder as this script.

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
seed = 5;
t = 0:1:199;
tau = 10;
gfun = t.*exp(-t/tau)/tau^2;

N = 1e3;
%Nhidvec = [990 890 790 690 590 490 390 290 190 90 1];
Nhid = 997;

J0vec = [1.0];
J0strvec = {'1.0'};

%%

Jdata = struct('Jeff',[]);
rgbcolors = {[0.25 0.25 0.25],[111/255 0 118/255]};


for J0ind=1:length(J0vec)
    J0 = J0vec(J0ind);
    J0strout = J0strvec{J0ind};
    
    filenameinWrr = ['Fig6/Wrr_' weighttype signtype scalingtype networktype '_N' num2str(N) '_Nhid' num2str(Nhid) '_p' num2str(p) '_mu' mustr '_J0' J0strout '_seed' num2str(seed) '.txt'];

    Wrr = load(filenameinWrr,'Delimiter',' ');
    
    for ii=1:3
        for jj=1:3
        
        filenameinJeff = ['Fig6/Jefft_r1' num2str(ii-1) '_r2' num2str(jj-1) '_' weighttype signtype scalingtype networktype '_N' num2str(N) '_Nhid' num2str(Nhid) '_p' num2str(p) '_mu' mustr '_J0' J0strout '_seed' num2str(seed) '.txt'];

        
        Jdata(ii,jj).Jeff = load(filenameinJeff,'Delimiter',' ');
        end
    end
        
end
        
        %Plot grid of Jeffs vs. Jtrues
        figure;
        subplot(331)
        plot(t,Wrr(1,1)*gfun,':','Linewidth',5,'Color',rgbcolors{1})
        hold on;
        plot(t,Jdata(1,1).Jeff,'-','LineWidth',5,'Color',rgbcolors{2}) %J0 = 0.25
        xlim([0 200]);
        ylim([-2e-4 2e-4]);
        text(0.4,0.25,'$1 \leftarrow 1$','Units','Normalized','Interpreter','LaTeX','FontSize',16);
     %   ylabel('$J^{\rm eff}(t)$','Interpreter','LaTeX','FontSize',16)
    %    xlabel('$t$','Interpreter','LaTeX','FontSize',16)
 %       title('$J_0 = 0.25$','Interpreter','LaTeX','FontSize',16)
        set(gca,'FontSize',16);
        axis square
        
        subplot(332)
        plot(t,Wrr(1,2)*gfun,':','Linewidth',5,'Color',rgbcolors{1})
        hold on;
        plot(t,Jdata(1,2).Jeff,'-','LineWidth',5,'Color',rgbcolors{2}) %J0 = 0.25
        xlim([0 200]);
        ylim([-4e-3 0e-3]);
        text(0.4,0.25,'$1 \leftarrow 2$','Units','Normalized','Interpreter','LaTeX','FontSize',16);
     %   ylabel('$J^{\rm eff}(t)$','Interpreter','LaTeX','FontSize',16)
    %    xlabel('$t$','Interpreter','LaTeX','FontSize',16)
 %       title('$J_0 = 0.25$','Interpreter','LaTeX','FontSize',16)
        set(gca,'FontSize',16);
        axis square
        
        subplot(333)
        plot(t,Wrr(1,3)*gfun,':','Linewidth',5,'Color',rgbcolors{1})
        hold on;
        plot(t,Jdata(1,3).Jeff,'-','LineWidth',5,'Color',rgbcolors{2}) %J0 = 0.25
        xlim([0 200]);
        ylim([-5e-4 5e-4]);
        text(0.4,0.25,'$1 \leftarrow 3$','Units','Normalized','Interpreter','LaTeX','FontSize',16);
      %  ylabel('$J^{\rm eff}(t)$','Interpreter','LaTeX','FontSize',16)
     %   xlabel('$t$','Interpreter','LaTeX','FontSize',16)
 %       title('$J_0 = 0.25$','Interpreter','LaTeX','FontSize',16)
        set(gca,'FontSize',16);
        axis square
        
        subplot(334)
        plot(t,Wrr(2,1)*gfun,':','Linewidth',5,'Color',rgbcolors{1})
        hold on;
        plot(t,Jdata(2,1).Jeff,'-','LineWidth',5,'Color',rgbcolors{2}) %J0 = 0.25
        xlim([0 200]);
        ylim([-2e-4 2e-4]);
        text(0.4,0.25,'$2 \leftarrow 1$','Units','Normalized','Interpreter','LaTeX','FontSize',16);
     %   ylabel('$J^{\rm eff}(t)$','Interpreter','LaTeX','FontSize',16)
    %    xlabel('$t$','Interpreter','LaTeX','FontSize',16)
 %       title('$J_0 = 0.25$','Interpreter','LaTeX','FontSize',16)
        set(gca,'FontSize',16);
        axis square
        
        subplot(335)
        plot(t,Wrr(2,2)*gfun,':','Linewidth',5,'Color',rgbcolors{1})
        hold on;
        plot(t,Jdata(2,2).Jeff,'-','LineWidth',5,'Color',rgbcolors{2}) %J0 = 0.25
        xlim([0 200]);
        ylim([-5e-4 5e-4]);
        text(0.4,0.25,'$2 \leftarrow 2$','Units','Normalized','Interpreter','LaTeX','FontSize',16);
   %     ylabel('$J^{\rm eff}(t)$','Interpreter','LaTeX','FontSize',16)
  %      xlabel('$t$','Interpreter','LaTeX','FontSize',16)
 %       title('$J_0 = 0.25$','Interpreter','LaTeX','FontSize',16)
        set(gca,'FontSize',16);
        axis square
        
        subplot(336)
        plot(t,Wrr(2,3)*gfun,':','Linewidth',5,'Color',rgbcolors{1})
        hold on;
        plot(t,Jdata(2,3).Jeff,'-','LineWidth',5,'Color',rgbcolors{2}) %J0 = 0.25
        xlim([0 200]);
        ylim([-2e-3 2e-3]);
        text(0.4,0.25,'$2 \leftarrow 3$','Units','Normalized','Interpreter','LaTeX','FontSize',16);
     %   ylabel('$J^{\rm eff}(t)$','Interpreter','LaTeX','FontSize',16)
    %    xlabel('$t$','Interpreter','LaTeX','FontSize',16)
 %       title('$J_0 = 0.25$','Interpreter','LaTeX','FontSize',16)
        set(gca,'FontSize',16);
        axis square
        
        subplot(337)
        plot(t,Wrr(3,1)*gfun,':','Linewidth',5,'Color',rgbcolors{1})
        hold on;
        plot(t,Jdata(3,1).Jeff,'-','LineWidth',5,'Color',rgbcolors{2}) %J0 = 0.25
        xlim([0 200]);
        ylim([-2e-3 2e-3]);
        set(gca,'FontSize',14);
        ylabel('$J^{\rm eff}(t)$','Interpreter','LaTeX','FontSize',16)
        xlabel('Time from spike (ms)','Interpreter','LaTeX','FontSize',16)
 %       title('$J_0 = 0.25$','Interpreter','LaTeX','FontSize',16)
        text(0.4,0.25,'$3 \leftarrow 1$','Units','Normalized','Interpreter','LaTeX','FontSize',16);
        set(gca,'FontSize',16);
        axis square
        
        subplot(338)
        plot(t,Wrr(3,2)*gfun,':','Linewidth',5,'Color',rgbcolors{1})
        hold on;
        plot(t,Jdata(3,2).Jeff,'-','LineWidth',5,'Color',rgbcolors{2}) %J0 = 0.25
        xlim([0 200]);
        ylim([-1e-3 1e-3]);
        text(0.4,0.25,'$3 \leftarrow 2$','Units','Normalized','Interpreter','LaTeX','FontSize',16);
    %    ylabel('$J^{\rm eff}(t)$','Interpreter','LaTeX','FontSize',16)
   %     xlabel('$t$','Interpreter','LaTeX','FontSize',16)
 %       title('$J_0 = 0.25$','Interpreter','LaTeX','FontSize',16)
        set(gca,'FontSize',16);
        axis square
        
        subplot(339)
        plot(t,Wrr(3,3)*gfun,':','Linewidth',5,'Color',rgbcolors{1})
        hold on;
        plot(t,Jdata(3,3).Jeff,'-','LineWidth',5,'Color',rgbcolors{2}) %J0 = 0.25
        xlim([0 200]);
        ylim([-5e-4 5e-4]);
        text(0.4,0.25,'$3 \leftarrow 3$','Units','Normalized','Interpreter','LaTeX','FontSize',16);
    %    ylabel('$J^{\rm eff}(t)$','Interpreter','LaTeX','FontSize',16)
   %     xlabel('$t$','Interpreter','LaTeX','FontSize',16)
 %       title('$J_0 = 0.25$','Interpreter','LaTeX','FontSize',16)
        set(gca,'FontSize',16);
        axis square
