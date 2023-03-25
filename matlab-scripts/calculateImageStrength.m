function [Q] = calculateImageStrength(kr2,theta0,beta, thetaMin, thetaMax)
%% 
% Calculate intensity of source at image position behind a reflecting wall
% Q - output image source strength
% kr2 - 2D grid of wave numbers times distance from IS to receiver
% theta0 - 2D grid of receiver angles from image source
% beta - wall admittance
% thetaMin, thetaMax - lower limits of integration; 
% angles made between source and wall edges.

%size of 2D grid
[Ny, Nx] = size(kr2);

%Limits of integration 
etaMax = theta0 - repmat(thetaMin, [Ny,Nx]);
etaMin = theta0 - repmat(thetaMax, [Ny,Nx]);

%cos eta = 1+it
tMin = -1i*(1-cos(etaMin));
tMax = -1i*(1-cos(etaMax));


gamma0 = cos(theta0);
R0 = (gamma0 - beta)./(gamma0 + beta);
rho = (gamma0+beta)./sqrt(2*(1+gamma0.*beta));

% infinite wall scattering
rho = 1i*kr2.*(rho.^2);
Q = R0 + (1-R0).*((rho.^(0.25)).*exp(rho/2).*whittakerW(-0.25,0.25,rho));

% % finite wall scattering
% integral = evaluateIntegral(kr2,rho,tMax) - evaluateIntegral(kr2, rho, tMin);
% Q = R0 + (1-R0).*(1 - rho.*exp(-1i*kr2.*rho).*integral);


end