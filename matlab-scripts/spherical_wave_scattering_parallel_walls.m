%% finite and parallel walls on xy planes
close all; clear all;

% number of points to evaluate on
Nx = 50;
Ny = 100;

Nwalls = 2;
L = [5,5];        % Length of walls in m
H = 6;            % distance between walls 
h0 = 3;           % height of source from floor
h = [h0,H-h0];    % height of source from walls
ref_order = 3;    % maximum number of reflections
wallLimits = zeros(Nwalls, 2, 3);       %limits of walls (used for calculating limitng angles)
wallLimits(1,:,:) = [[0,0,0];[L(1),0,0]];
wallLimits(2,:,:) = [[0,0,H];[L(2),0,H]];

%source position
xs = 2; zs=h0 ; ys = L(2);
sourcePos = [xs,ys,zs];
%image source positions for each wall
imagePos = zeros(Nwalls,ref_order,3);
%limiting angles
thetaMin = zeros(Nwalls, ref_order);
thetaMax = zeros(Nwalls, ref_order);


% receivers are in an arc around xy plane
r = 2.5;    
theta = linspace(0,pi,Ny);
rx = r.*cos(theta).'  + L(1)/2;
rz = r.*sin(theta).';
ry = r.*ones(1,Ny).';
receiverPos = [rx,ry,rz];


for j = 1:ref_order

    if j == 1
        imagePos(1,j,:) = [xs,ys,-h0];
        imagePos(2,j,:) = [xs,ys,2*H-h0];

        thetaMax(:,j) = [acos(h0/getDistance(sourcePos,wallLimits(1,1,:))); ...
            acos((H-h0)/getDistance(sourcePos, wallLimits(2,1,:)))];
        thetaMin(:,j) = -[acos(h0/getDistance(sourcePos,wallLimits(1,end,:))); ...
            acos((H-h0)/getDistance(sourcePos, wallLimits(2,end,:)))];

 
    else
        imagePos(1,j,:) = [xs,ys,-imagePos(2,j-1,3)];
        imagePos(2,j,:) = [xs,ys,-imagePos(1,j-1,3)+2*H];

        thetaMax(:,j) = [acos(imagePos(2,j-1,3)/getDistance(imagePos(2,j-1,:),wallLimits(1,1,:))); ...
            acos((-imagePos(1,j-1,3)+H)/getDistance(imagePos(1,j-1,:), wallLimits(2,1,:)))];
        thetaMin(:,j) = -[acos(imagePos(2,j-1,3)/getDistance(imagePos(2,j-1,:),wallLimits(1,2,:))); ...
            acos((-imagePos(1,j-1,3)+H)/getDistance(imagePos(1,j-1,:), wallLimits(2,2,:)))];
    end
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% draw the setup

figure(1); 


% plot source and image source
plot3(sourcePos(1), sourcePos(2), sourcePos(3),'bo');hold on;
text(sourcePos(1)+0.2, sourcePos(2),  sourcePos(3)+0.2, 'S', 'FontSize',8);

plot3(imagePos(1,:,1), imagePos(1,:,2), imagePos(1,:,3), 'ko');hold on;
text(imagePos(1,:,1)-0.5, imagePos(1,:,2), imagePos(1,:,3)-0.2, 'Q_f', 'FontSize',8);

plot3(imagePos(2,:,1), imagePos(2,:,2), imagePos(2,:,3), 'ko');hold on;
text(imagePos(2,:,1)-0.5, imagePos(2,:,2)-0.1, imagePos(2,:,3)-0.2, 'Q_c', 'FontSize',8); hold on;



% plot reflecting wall
[x,y] = meshgrid(linspace(0,L(1),100));
z = zeros(100,100);
plot3(x,y,z,'k');hold on;

[x,y] = meshgrid(linspace(0, L(2),100));
z = H.*ones(100,100);
plot3(x,y,z,'k'); hold on;

% plot receiver positions
plot3(rx,ry,rz,'rx'); hold off;
view(10,5);

xlim([0,L(1)]); zlim([-1.5*H, L(2)+1.8*H]);  ylim([0,L(1)]); 
set(gca, 'Visible', 'off');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% pressure at receivers

k = linspace(0.01,10,Nx);               % wave numbers = 2pi/lambda
pR = zeros(Ny, Nx);                     % reflected (scattered) pressure
Q = zeros(Ny,Nx,ref_order,Nwalls);      % strength of each image source
pR_approx = zeros(Ny,Nx);               % plane wave approximation
beta =[1.5i, 10+1i];                    % wall admittance


% distance between source and receiver
r1 = getDistance(sourcePos, receiverPos);
R1 = repmat(r1, [1,Nx]);


for n = 1:Nwalls
    otherWall = setdiff(1:Nwalls,n);
    for j = 1:ref_order

        %angle and distance from image source to receiver
        r2 = getDistance(imagePos(n,j,:), receiverPos);
        if n == 1
            theta0 = acos((abs(imagePos(n,j,3)) + rz)./r2);
        else
            theta0 = acos((imagePos(n,j,3) - H + rz)./r2);
        end

        R2 = repmat(r2, [1,Nx]);

        % create 2D mesh over wave numbers and receiver angles
        [K, theta0] = meshgrid(k, theta0);
     
        kr2 = K.*R2;

        % strength of each image-source
        Q(:,:,j,n) = calculateImageStrength(kr2, theta0 ,beta(n), thetaMin(n,j), thetaMax(n,j));
        
        % the strength needs to be multiplied with the strength of original
        % source
        if j > 1
            Q(:,:,j,n) = Q(:,:,j-1,otherWall).*Q(:,:,j,n);
        end
        pR = pR + (Q(:,:,j,n).*exp(1i*kr2)./R2);

        % plane wave approximation
        if j == 1
            pR_approx = pR_approx + exp(1i*kr2).*((beta(n) - cos(theta0))./(beta(n) + cos(theta0)));
        end
    end


end

% add incident pressure to get total pressure
kr1 = K.*R1;
pI = exp(1i.*kr1)./R1;


%% plot

close all;

generatePressurePlots(pI,k,theta, 'title', 'Incident pressure');
generatePressurePlots(pR,k,theta, 'title', 'Scattered Pressure');
generatePressurePlots(pR+pI,k,theta, 'title', 'Total pressure');



%% helper functions

function [dist] = getDistance(point1, point2)
    dist = sqrt((point1(:,1) - point2(:,1)).^2 + (point1(:,2) - point2(:,2)).^2 ...
        + (point1(:,3) - point2(:,3)).^2);
end