%% finite wall on xy plane, source above wall
close all; clear all;

% number of points to evaluate on
Nx = 50;
Ny = 100;

L = 5;     % Length of wall in m
h = 0.5;     % height of source from wall
k = linspace(0,50,Nx);    % wavenumber = 2pi/lambda

%source position
sourcePos = [L/2,L/2,h];
%image position 
imagePos = [sourcePos(1), sourcePos(2), -h];

minDist = h;
r2 = minDist+2;    % distance from image source to listener


% limits of incoming angle
thetaMax = acos(h/getDistance(sourcePos, [L,0,0]));
thetaMin = -acos(h/getDistance(sourcePos, [0,0,0]));


%angle of receiver from image source - only points above plane
theta0 = linspace(-pi/2,pi/2,Ny);


%% draw the setup
figure(1);

% plot source and image source
plot3(sourcePos(1), sourcePos(2), sourcePos(3),'bo');hold on;
text(sourcePos(1), sourcePos(2), h+0.1, 'S');
plot3(imagePos(1), imagePos(2), imagePos(3), 'ko');hold on;
text(imagePos(1), imagePos(2), -h-0.1, 'Q');
line([sourcePos(1),0], [sourcePos(2),0],[sourcePos(3),0]);hold on;
text(sourcePos(1)-0.25, sourcePos(2)-0.25, h-0.2, '\theta_{min}');
line([sourcePos(1),L], [sourcePos(2),0],[sourcePos(3),0]);hold on;
text(sourcePos(1)+0.25, sourcePos(2)+0.25, h-0.2, '\theta_{max}');
line([sourcePos(1),  sourcePos(1)], [sourcePos(2), sourcePos(2)], [h,-h], 'LineStyle','--'); hold on;

% plot reflecting wall
[x,y] = meshgrid(linspace(0,L,100));
z = zeros(100,100);
plot3(x,y,z,'k');hold on;

% plot receiver positions
xr = imagePos(1)+(r2.*cos(theta0+pi/2));
yr = imagePos(2)+(r2.*sin(theta0+pi/2));
zr = r2.*ones(1, Ny);
plot3(xr,yr,zr,'rx'); hold off;

xlim([0,L]); ylim([0,L]); zlim([-h-0.2, r2]);

%%

% create 2D mesh over wave numbers and receiver angles
[k, theta0] = meshgrid(k, theta0);


%Limits of integration 
etaMax = theta0 - repmat(thetaMin, [Ny,Nx]);
etaMin = theta0 - repmat(thetaMax, [Ny,Nx]);

%cos eta = 1+it
tMin = -1i*(1-cos(etaMin));
tMax = -1i*(1-cos(etaMax));

%from triangle cosine law
gamma0 = cos(theta0);
kr1 = sqrt(4*(k*h).^2 + (k*r2).^2 - 4*k.^2*h*r2.*gamma0);

%wall admittance
beta = [0,0.1,1,10]*1i;
idx = [2,10,50];
% beta = 1 + 1i*k(1,idx); 

for i = 1:length(beta)
        R0 = (gamma0 - beta(i))./(gamma0 + beta(i));

        rho = (gamma0+beta(i))./sqrt(2*(1+gamma0.*beta(i)));

        % for infinite wall scattering
        rho = 1i*k*r2.*(rho.^2);
        Q = R0 + ((1-R0).*sqrt(rho).*exp(rho).*whittakerW(-0.25,0.25,rho));

        frac = (k.*r2)./kr1;
        p = exp(1i*kr1)./(kr1).*(1 + frac.*Q.*exp(1i*kr1.*(frac-1)));

        figure(2*i+1);    
        for j = 1:length(idx)
            polar_p = p(:,idx(j));
            subplot(length(idx),1,j);polarplot(theta0(:,1), abs(polar_p));
            thetalim([-90,90]);
            title(['kr = ', num2str(k(1,idx(j))*r2)]);
        end
        sgtitle(['beta = ', num2str(beta(i))]);
        print(['../figures/infinite-wall/S-beta=',num2str(beta(i)),'.eps'], '-depsc');

end


%% helper functions

function [dist] = getDistance(point1, point2)
    dist = sqrt((point1(1) - point2(1))^2 + (point1(2) - point2(2))^2 + (point1(3) - point2(3))^2);
end


function [output] = evaluateIntegral(c,B,t)
    output = exp(-1i.*B.*c).*sqrt(c.*(t - 1i*B).*igamma(0.5,c.*(t-1i*B)))./...
        (c.*sqrt(B + 1i*t));
end

