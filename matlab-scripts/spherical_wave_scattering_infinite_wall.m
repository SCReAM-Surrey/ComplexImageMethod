
%% infinite wall on xy plane - source on wall
close all; clear all;

% number of points to evaluate on
Nx = 50;
Ny = 100;

%angle of receiver from image source
theta0 = linspace(0,pi/2-0.001,Ny);

%wall admittance
beta = [0,0.1,1,10]*1i;

% constant term, wave number times source distance from listener
kr = linspace(0,50,Nx);

%2D mesh grid of different kr and angles
[KR, THETA0] = meshgrid(kr, theta0);

gamma0 = cos(THETA0);    


for i = 1:length(beta)

    %source strength at image position
    rho = 1i*KR.*((gamma0 - beta(i)).^2)./(2*(1+gamma0.*beta(i)));
    R0 = (gamma0 - beta(i))./(gamma0 + beta(i));
    Q = R0 + ((1-R0).*(rho.^(0.25)).*exp(rho/2).*whittakerW(-0.25,0.25,rho));
    
    % image strength should be smaller
    assert(all(all(abs(Q) <= 1.0)),['Image strength for admittance = ', num2str(beta),' too large']);
  
    %pressure at receiver in decibels
    p = (exp(1j*KR)./KR).*(1+Q);
    

    %look at the directional pressure at certain values of kr
    fig = figure('Units','inches', 'Position',[0 0 2.2 3.5],'PaperPositionMode','auto');

    idx = [2,15,50];
    for j = 1:length(idx)
        polar_p = p(:,idx(j));
        subplot(length(idx),1,j);polarplot(linspace(-pi/2,pi/2,2*Ny), abs([flipud(polar_p); polar_p]));
        thetalim([-90,90]);
        title(['kr = ', num2str(round(kr(idx(j)),3))]); rlim([0, max(abs(polar_p))])
        set(gca, 'FontUnits','points', 'FontWeight','normal', 'FontSize',8, 'FontName','Times');


    end
    sgtitle(['admittance = ', num2str(beta(i))]);
    print(fig,['../figures/S-on-wall-beta=',num2str(beta(i)),'.eps'], '-depsc');

end



%% infinite wall on xy plane, source above wall
close all;

kh = 5;    %h is the distance of source above wall     

% number of points to evaluate on
Nx = 50;
Ny = 100;

%angle of receiver from imagec source;
theta0 = linspace(0,pi/2-0.01,Ny);

%wall admittance
beta = [0,0.1,1,10]*1i;

% constant term, wave number times source distance from listener
kr2 = linspace(0,50,Nx);

%2D mesh grid of different kr and angles
[KR2, THETA0] = meshgrid(kr2, theta0);

%from triangle cosine law
gamma0 = cos(THETA0);
KR1 = sqrt(4*kh^2 + KR2.^2 - 4*kh.*KR2.*gamma0);

for i = 1:length(beta)

    %source strength at image position
    beta(i) = beta(i)*(1+exp(-1j*kr2(i)));
    rho = 1i*KR2.*((gamma0 - beta(i)).^2)./(2*(1+gamma0.*beta(i)));
    R0 = (gamma0 - beta(i))./(gamma0 + beta(i));
    Q = R0 + ((1-R0).*(rho.^(0.25)).*exp(rho/2).*whittakerW(-0.25,0.25,rho));

    frac = KR2./KR1;
    p = exp(1i*KR1)./(KR1).*(1 + frac.*Q.*exp(1i*KR1.*(frac-1)));

    %look at the directional pressure at certain values of kr
    fig = figure('Units','inches', 'Position',[0 0 2.5 4.2],'PaperPositionMode','auto');

    
    idx = [2,15,50];
    for j = 1:length(idx)
        polar_p = p(:,idx(j));
        subplot(length(idx),1,j);
        polarplot(linspace(-pi/2,pi/2,2*Ny), abs([flipud(polar_p); polar_p]));
        thetalim([-90,90]);
        title(['kr = ', num2str(round(kr2(idx(j)),3))], 'FontSize',8, 'FontName','Times');
        set(gca, 'FontUnits','points', 'FontWeight','normal', 'FontSize',8, 'FontName','Times');


    end
    sgtitle(['admittance = ', num2str(beta(i))], 'FontSize',8, 'FontName','Times');
    print(['../figures/S-above-wall-beta=',num2str(round(beta(i),2)),'.eps'], '-depsc');

end

%% Comments
% Theta0 is the angle between the image source and the receiver, in case of
% a finite wall, theta0 = pi/2 is impossible. When theta0 = 0, the receiver
% is directly above the image source and will have maximum pressure.
%
% When the source is on the wall and if the wall is rigid (0 admittance), 
% then an exact spherical wave will be reflected with the same strength as 
% the incoming wave. Therefore, the pressure magnitude will be same 
% regardless of the angle of the receiver (circular polar plot).



