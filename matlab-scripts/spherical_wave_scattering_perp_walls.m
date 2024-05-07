%% finite walls on xy plane and yz planes
close all; clear all;

% number of points to evaluate on
Nx = 50;
Ny = 100;

Nwalls = 2;
L = [5,5];        % Length of walls in m
h = [2.5,2];      % height of source from walls
k = linspace(0.01,10,Nx);    % wavenumber = 2pi/lambda


%source position
xs = h(2); zs=h(1) ; ys = 3;
sourcePos = [xs;ys;zs];
%first order image position 
imagePos = [[xs; ys; -zs], [-xs; ys; zs]];
%second order image position
secondImagePos = [-xs; ys; -zs];


% limits of incoming angle
thetaMax = [acos(h(1)/getDistance(sourcePos, [L(1);0;0]));acos(h(2)/getDistance(sourcePos,[0;0;L(2)]))];
thetaMin = -[acos(h(1)/getDistance(sourcePos, [0;0;0]));acos(h(2)/getDistance(sourcePos, [0;0;0]))];


% receivers are in an arc around xz plane
r = 3;    
theta = linspace(0,pi/2,Ny);
rx = r.*cos(theta);
ry = L(1)/2*ones(1,Ny);
rz = r.*sin(theta);
receiverPos = [rx;ry;rz];

%% draw the setup

fig = figure('Units','inches', 'Position',[0 0 3.29 2.2],'PaperPositionMode','auto');


% plot source and image source
plot3(sourcePos(1), sourcePos(2), sourcePos(3),'bo');hold on;
text(sourcePos(1)+0.2, sourcePos(2),  sourcePos(3)+0.2, 'S', 'FontSize',8);

plot3(imagePos(1,:), imagePos(2,:), imagePos(3,:), 'ko');hold on;
text(imagePos(1,:)-0.2, imagePos(2,:), imagePos(3,:)-0.6, 'Q', 'FontSize',8);

plot3(secondImagePos(1), secondImagePos(2), secondImagePos(3), 'ko');hold on;
text(secondImagePos(1)-0.25, secondImagePos(2)-0.1, secondImagePos(3)-0.7, 'Q_s', 'FontSize',8); hold on;

line([sourcePos(1),0], [sourcePos(2),0],[sourcePos(3),0]);hold on;
text(sourcePos(1)-0.25, sourcePos(2)-0.25, sourcePos(3)-0.2, '\theta_{min}');

line([sourcePos(1),L(1)], [sourcePos(2),0],[sourcePos(3),0]);hold on;
text(sourcePos(1)+0.2, sourcePos(2), sourcePos(3)-0.2, '\theta_{max}');

line([sourcePos(1),0], [sourcePos(2),0],[sourcePos(3),L(2)]);hold on;
text(sourcePos(1)+0.25, sourcePos(2)+0.25, sourcePos(3)-0.2, '\theta_{max}');

line([sourcePos(1),  sourcePos(1)], [sourcePos(2), sourcePos(2)], [h(1),-h(1)], 'LineStyle','--'); hold on;
line([h(2),-h(2)], [sourcePos(2), sourcePos(2)], [sourcePos(3), sourcePos(3)], 'LineStyle','--'); hold on;


% plot reflecting wall
[x,y] = meshgrid(linspace(0,L(1),100));
z = zeros(100,100);
plot3(x,y,z,'k');hold on;

[z,y] = meshgrid(linspace(0, L(2),100));
x = zeros(100,100);
plot3(x,y,z,'k'); hold off;

% plot receiver positions
plot3(rx,ry,rz,'rx'); hold off;
view(-10,10);

xlim([-h(2)-0.2,L(1)]); zlim([-h(1)-0.2, L(2)]);  ylim([0,L(1)]); 
set(gca, 'Visible', 'off');
print(['../figures/perp_walls_setup.eps'], '-depsc');



% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% contribution due to the individual image sources on each wall,
% plus direct contribution of source to receiver
pI = zeros(Ny, Nx);
Q = zeros(Ny, Nx, Nwalls);
beta =[1.5, 10]*1i;

% distance between source and receiver
r1 = getDistance(sourcePos, receiverPos);
R1 = repmat(r1, [Nx,1]).';


for n = 1:Nwalls

    %angle and distance from image source to receiver
    r2 = getDistance(imagePos(:,n), receiverPos);
    if n == 1
        theta0 = acos((h(n) + rz)./r2);
    else
        theta0 = acos((h(n) + rx)./r2);
    end
    R2 = repmat(r2, [Nx,1]).';

    % create 2D mesh over wave numbers and receiver angles
    [K, theta0] = meshgrid(k, theta0);
    kr1 = K.*R1;
    kr2 = K.*R2;

    % strength of primary image source
    Q(:,:,n) = calculateImageStrength(kr2, theta0 ,beta(n), thetaMin(n), thetaMax(n));

    frac = kr2./kr1;
    pI = pI + exp(1i*kr1)./R1 + (Q(:,:,n).*exp(1i*kr2)./R2);
    

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% contribution due to the secondary image source, 
% models interaction between walls)


% limits of incoming angle
alphaMin = [acos((h(1)+L(2))/getDistance(imagePos(:,1), [0;0;L(2)]));...
    acos((h(2)+L(1))/getDistance(imagePos(:,2), [L(1);0;0]))];
alphaMax = [acos(h(1)/getDistance(imagePos(:,1), [0;0;0]));...
    acos(h(2)/getDistance(imagePos(:,2), [0;0;0]))];

pS = zeros(Ny,Nx);

for n = 1:Nwalls

    otherWall = setdiff(1:Nwalls,n);

    %distance between secondary image source and receiver
    rs = getDistance(secondImagePos, receiverPos);
    if n == 1
        %contribution on xy wall due to xz wall
        alpha =  acos((h(n) + rz)./rs);
    else
        %contribution on xz wall due to xy wall
        alpha = acos((h(n) + rx)./rs);
    end

    Rs = repmat(rs, [Nx,1]).';

    % create 2D mesh over wave numbers and receiver angles
    [K, alpha] = meshgrid(k, alpha);
    kr = K.*Rs;

    % strength of secondary image source
    % incoming source is on the other wall
    Qs = calculateImageStrength(kr, alpha ,beta(n), alphaMin(otherWall), alphaMax(otherWall));
    
    % the secondary image strength is scaled by the strength of the image 
    % source on other wall, since that acts as the virtual source now
    pS = pS + (exp(1i*kr)./Rs).*(Qs.*Q(:,:,otherWall));

end

%% plot

close all;
generatePressurePlots(pI,k,theta, 'title', 'Pressure due to 1st image source');
generatePressurePlots(pS,k,theta, 'title', 'Pressure due to 2nd image source');
generatePressurePlots(pS+pI,k,theta, 'title', 'Total pressure');


%% helper functions


function [dist] = getDistance(point1, point2)
    dist = sqrt((point1(1,:) - point2(1,:)).^2 + (point1(2,:) - point2(2,:)).^2 ...
        + (point1(3,:) - point2(3,:)).^2);
end



