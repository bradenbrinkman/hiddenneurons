% This script plots the theoretically predicted effective interactions
% between the simple 3 neuron circuit shown in Figure 3. The first figure
% produced shows the effective interaction in purple, with the direct 
% excitatory contribution shown in red and the indirect inhibitory
% contribution shown in blue. This figure is not used in the main paper.
% The second figure produced by this file plots the aforementioned curves
% on separate plots that are shown in Figure 3. Full details are given in
% the Methods section.


% Plot feedforward inhibitory effective & true interactions

% Set times
t=0:0.1:10;

% Set parameters

% Synaptic strengths (e.g., J21 is the synaptic strength from neuron 1 TO
% 2)
J21 = 1.0;
J23 = -2.0; 
J31 = 2.0;
J33 = -0.9;

% Time scale of the temporal component of the synaptic filters (a's are
% inter-neuron coupling timescales labeled alpha in text, the b's are the 
% self-coupling timescales labeled beta in text. Note that there are
% restrictions on the allowed values of the timescales: a31 is assumed to
% not equal a21 or a23, nor can it be equal to b33*(J33-1)
a21 = 1.0;
a23 = 1.0;
a31 = 1.8;
b33 = 1.0;

% Define temporal components. True interactions are alpha functions.
% Components of the effective interactions were computed analytically using Mathematica.

J21t = J21*a21^2*t.*exp(-a21*t);

% The following expressions for the components explicitly assume the
% restrictions on timescales noted above.

dJta = (b33*J33*exp(-b33*(1-J33)*t))/((a23 - b33*(1 - J33))^2*(a31 - b33*(1 - J33))^2);
dJtb = (exp(-a31*t)*(-2*a31^2 - a31*b33*(-4 + J33) + b33*(2*b33*(-1 + J33) - a23*J33)))/((a23 - a31)^3*(a31 + b33*(-1 + J33))^2);
dJtc = (exp(-a23*t)*(-2*a23^2 - a23*b33*(-4 + J33) + b33*(2*b33*(-1 + J33) - a31*J33)))/((a31 - a23)^3*(a23 + b33*(-1 + J33))^2);
dJtd = ((a23 - b33)*exp(-a23*t).*t)/((a23 - a31)^2*(a23 + b33*(-1 + J33)));
dJte = ((a31 - b33)*exp(-a31*t).*t)/((a31 - a23)^2*(a31 + b33*(-1 + J33)));

dJt = J23*J31*(a23^2)*(a31^2)*(dJta+dJtb+dJtc+dJtd+dJte); 

Jefft = J21t + dJt;

rgbcolors = {[111/255 0 118/255],[0 0 1],[1 0 0]};

figure;
plot(t,Jefft,'Linewidth',5,'Color',rgbcolors{1})
hold on;
plot(t,J21t,':','Linewidth',5,'Color',rgbcolors{3})
plot(t,dJt,'--','Linewidth',5,'Color',rgbcolors{2})
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

subplot(131)
plot(t,Jefft,'Linewidth',5,'Color',rgbcolors{1})
hold on;
plot(t,zeros(size(t)),'k--','Linewidth',3)
set(gca,'FontSize',14,'YTick',-0.8:0.4:0.4);
ylim([-0.8 0.401]);
%xlabel('time from spike (normalized)','Interpreter','LaTeX','FontSize',16)
%ylabel('membrane response (normalized)','Interpreter','LaTeX','FontSize',16)
axis square

subplot(132)
plot(t,J21t,'Linewidth',5,'Color',rgbcolors{3})
hold on;
plot(t,zeros(size(t)),'k--','Linewidth',3)
plot(t,Jefft,':','Linewidth',3,'Color',rgbcolors{1})
ylim([-0.8 0.401]);
set(gca,'FontSize',14,'YTick',-0.8:0.4:0.4);
axis square

subplot(133)
plot(t,dJt,'Linewidth',5,'Color',rgbcolors{2})
hold on;
plot(t,zeros(size(t)),'k--','Linewidth',3)
plot(t,Jefft,':','Linewidth',3,'Color',rgbcolors{1})
ylim([-0.8 0.401]);
set(gca,'FontSize',14,'YTick',-0.8:0.4:0.4);
axis square
