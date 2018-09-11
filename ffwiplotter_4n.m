% This script plots the theoretically predicted effective interactions
% between the simple 4 neuron circuit shown in Figure 4. 

%The first figure produced shows the effective interaction in purple, 
% with the direct excitatory contribution shown in red and the indirect inhibitory
% contribution shown in blue. This figure is not used in the main paper.

% The second figure produced by this file plots contributions to the
% effective interaction out as separate contributions from different paths
% through the circuit. The first plot is the effective interaction. The
% fifth plot is the contribution of both infinite sums shown in Figure 4B
% combined into a single plot. This version is not used in the main paper.

% The third figure produced by this file is the same as above, only the
% infinite sums shown in Figure 4B are separated into two terms, and so
% there are six plots in total. These are the plots used in the final
% version of Figure 4.

% The final figure produced by this file demonstrates that the analytical
% expression for the total sum of indirect contributions (dJt_4n) is indeed
% equal to the sum of the contributions of the individual paths (variables
% dJt1 to dJt4).

% Plot feedforward inhibitory effective & true interactions

t=0:0.01:10;

% Set parameters

% Synaptic strengths (e.g., J21 is the synaptic strength from neuron 1 TO
% 2)
J21 = 1.0;
J23 = -3.0; 
J31 = 1.0;
J41 = 1.0;
J = 0.9; % magnitude
J34 = J; %this is the magnitude; assumed to be negative in eqn below.
J43 = J; %this is the magnitude; assumed to be negative in eqn below.



a21 = 1.0; % neuron 1 to neuron 2 interaction timescale
a = 1.294; % timescale for all interactions other than the 1 -> 2 interaction.

J21t = J21*a21^2*t.*exp(-a21*t); % direct interaction from 1 to 2.


% Analytically calculated expressions for the indirect contribution to
% the effective interaction (split into two terms).

dJt_4na = J23*J31/(4*J34^(3/4)*J43^(3/4))*a*(-2*exp(-a*t).*(sin(a*(J34*J43)^(1/4)*t) - sinh(a*(J34*J43)^(1/4)*t)));

dJt_4nb = a*J23*J41*exp(-a*t).*(2*a*J34^(1/4)*J43^(1/4)*t - sin(a*(J34*J43)^(1/4)*t) - sinh(a*(J34*J43)^(1/4)*t))/(2*J34^(1/4)*J43^(5/4));

dJt_4n = dJt_4na + dJt_4nb; % sum of above two terms.

% Indirect contributions to the effective interactions split into
% contributions from each path shown in Figure 4.
dJt1 = 1/6*a^4*J23*J31*exp(-a*t).*t.^3;
dJt2 = -(1/120)*a^6*J*J23*J41*exp(-a*t).*t.^5;
dJt3 = (a*J23*J31*(cosh(a*t) - sinh(a*t)).*(-2*a^3*J^(3/2).*t.^3 - 6*sin(a*sqrt(J)*t) + 6*sinh(a*sqrt(J)*t)))/(12*J^(3/2));
dJt4 = (a*exp(-a*t)*J23*J41.*(a*sqrt(J).*t.*(120 + a^4*J^2.*t.^4) - 60*sin(a*sqrt(J)*t) - 60*sinh(a*sqrt(J)*t)))/(120*J^(3/2));

Jefft_4n = J21t + dJt_4n;

rgbcolors = {[111/255 0 118/255],[0 0 1],[181/255 0 0]};

figure;
plot(t,Jefft_4n,'Linewidth',5,'Color',rgbcolors{1})
hold on;
plot(t,J21t,':','Linewidth',5,'Color',rgbcolors{3})
plot(t,dJt_4n,'--','Linewidth',5,'Color',rgbcolors{2})
 set(gca,'FontSize',14);
%ylabel('$J^{\rm eff}(t)$','Interpreter','LaTeX','FontSize',16)
xlabel('time from spike $\lambda_0 \tau$ (normalized)','Interpreter','LaTeX','FontSize',16)
set(gca,'FontSize',16);
lg = legend('$J^{\rm eff}_{21}(\tau)$','$J_{21}(\tau)$','$J^{\rm eff}_{21}(\tau)-J_{21}(\tau)$','location','Northeast');
set(lg,'Interpreter','LaTeX','FontSize',14);
legend boxoff;
axis square

%%

figure;

subplot(151)
plot(t,Jefft_4n,'Linewidth',5,'Color',rgbcolors{1})
hold on;
plot(t,zeros(size(t)),'k--','Linewidth',3)
set(gca,'FontSize',14,'YTick',-0.8:0.4:0.4);
ylim([-0.8 0.401]);
%xlabel('time from spike (normalized)','Interpreter','LaTeX','FontSize',16)
%ylabel('membrane response (normalized)','Interpreter','LaTeX','FontSize',16)
axis square

subplot(152)
plot(t,J21t,'Linewidth',5,'Color',rgbcolors{3})
hold on;
plot(t,zeros(size(t)),'k--','Linewidth',3)
plot(t,Jefft_4n,':','Linewidth',3,'Color',rgbcolors{1})
ylim([-0.8 0.401]);
set(gca,'FontSize',14,'YTick',-0.8:0.4:0.4);
axis square

subplot(153)
plot(t,dJt1,'Linewidth',5,'Color',rgbcolors{2})
hold on;
plot(t,zeros(size(t)),'k--','Linewidth',3)
plot(t,Jefft_4n,':','Linewidth',3,'Color',rgbcolors{1})
ylim([-0.8 0.401]);
set(gca,'FontSize',14,'YTick',-0.8:0.4:0.4);
axis square

subplot(154)
plot(t,dJt2,'Linewidth',5,'Color',rgbcolors{3})
hold on;
plot(t,zeros(size(t)),'k--','Linewidth',3)
plot(t,Jefft_4n,':','Linewidth',3,'Color',rgbcolors{1})
ylim([-0.8 0.401]);
set(gca,'FontSize',14,'YTick',-0.8:0.4:0.4);
axis square

subplot(155)
plot(t,dJt3+dJt4,'Linewidth',5,'Color',rgbcolors{2})
hold on;
plot(t,zeros(size(t)),'k--','Linewidth',3)
plot(t,Jefft_4n,':','Linewidth',3,'Color',rgbcolors{1})
ylim([-0.8 0.401]);
set(gca,'FontSize',14,'YTick',-0.8:0.4:0.4);
axis square

%%

figure;

subplot(161)
plot(t,Jefft_4n,'Linewidth',5,'Color',rgbcolors{1})
hold on;
plot(t,zeros(size(t)),'k--','Linewidth',3)
set(gca,'FontSize',14,'YTick',-0.8:0.4:0.8);
ylim([-1 0.8]);
%xlabel('time from spike (normalized)','Interpreter','LaTeX','FontSize',16)
%ylabel('membrane response (normalized)','Interpreter','LaTeX','FontSize',16)
axis square

subplot(162)
plot(t,J21t,'Linewidth',5,'Color',rgbcolors{3})
hold on;
plot(t,zeros(size(t)),'k--','Linewidth',3)
plot(t,Jefft_4n,':','Linewidth',3,'Color',rgbcolors{1})
ylim([-1 0.8]);
set(gca,'FontSize',14,'YTick',-0.8:0.4:0.8);
axis square

subplot(163)
plot(t,dJt1,'Linewidth',5,'Color',rgbcolors{2})
hold on;
plot(t,zeros(size(t)),'k--','Linewidth',3)
plot(t,Jefft_4n,':','Linewidth',3,'Color',rgbcolors{1})
ylim([-1 0.8]);
set(gca,'FontSize',14,'YTick',-0.8:0.4:0.8);
axis square

subplot(164)
plot(t,dJt2,'Linewidth',5,'Color',rgbcolors{3})
hold on;
plot(t,zeros(size(t)),'k--','Linewidth',3)
plot(t,Jefft_4n,':','Linewidth',3,'Color',rgbcolors{1})
ylim([-1 0.8]);
set(gca,'FontSize',14,'YTick',-0.8:0.4:0.8);
axis square

subplot(165)
plot(t,dJt3,'Linewidth',5,'Color',rgbcolors{2})
hold on;
plot(t,zeros(size(t)),'k--','Linewidth',3)
plot(t,Jefft_4n,':','Linewidth',3,'Color',rgbcolors{1})
ylim([-1 0.8]);
set(gca,'FontSize',14,'YTick',-0.8:0.4:0.8);
axis square

subplot(166)
plot(t,dJt4,'Linewidth',5,'Color',rgbcolors{3})
hold on;
plot(t,zeros(size(t)),'k--','Linewidth',3)
plot(t,Jefft_4n,':','Linewidth',3,'Color',rgbcolors{1})
ylim([-1 0.8]);
set(gca,'FontSize',14,'YTick',-0.8:0.4:0.8);
axis square

%%

% Plots to demonstrate that the contributions really do cancel to give
% the effective interaction

figure;
plot(t,dJt_4n,'b')
hold on;
plot(t,dJt1,'r')
plot(t,dJt1+dJt2,'g')
plot(t,dJt1+dJt2+dJt3,'m')
plot(t,dJt1+dJt2+dJt3,'k')
plot(t,dJt1+dJt2+dJt3+dJt4,'k--','LineWidth',5)
