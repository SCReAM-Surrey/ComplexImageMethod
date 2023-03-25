close all; clear all;

addpath('../sdn-matlab/.');
addpath('../utils/.');

% create a shoebox room with specified dimensions and get the node
% positions for a first order SDN

L = 5;
fs = 44100;
c = 343;
nWalls = 6;

room = Room();
room.shape = Cuboid(L,L,L);

source = Source();
source.position = Position(0.3,0.5,0.9);
source.signal = Signal([1, zeros(1,999)], fs);

microphone = Microphone();
microphone.position = Position(0.4,0.1,0.4);

% sim = Simulation();
% sim.room = room;
% sim.source = source;
% sim.microphone = microphone;
% sim.NSamples = 100;
% 
% sim.frameLength = 1;
% 
% output = sim.run();

% to keep track of walls - specified as planes ax+by+cx+d=0
Walls = zeros(6,4);
Walls(1,:) = [0,-L^2,0,0];      %floor
Walls(2,:) = [0,L^2,0,-L^3];    %ceiling
Walls(3,:) = [L^2,0,0,0];       %left
Walls(4,:) = [L^2,0,0,-L^3];    %right
Walls(5,:) = [0,0,-L^2,0];      %back
Walls(6,:) = [0, 0, L^2, -L^3]; %front


% to keep track of nodes on each wall
NodePos = zeros(6,3);
NodePos(1,:) = [0.3833,0,0.4833];
NodePos(2,:) = [0.3479, 5, 0.6606];
NodePos(3,:) = [0, 0.3286, 0.6857];
NodePos(4,:) = [5, 0.2978, 0.6473];
NodePos(5,:) = [0.3692, 0.2231, 0];
NodePos(6,:) = [0.3471, 0.3115,5];


% specify the scattering matrix at these frequencies
Nfreq = 512;
freqs = linspace(0, fs/2, Nfreq);
waveNum = 2*pi*freqs/c;
S = zeros(nWalls, nWalls-1, nWalls-1, Nfreq);

%wall admittance filters 
order = 2;
wallImpFreq = zeros(nWalls, Nfreq);
for i = 1:nWalls
    [b,a] = butter(order, 0.3,"low");
    wallImpFreq(i,:) = freqz(b,a,freqs,fs);
end


%% form scattering matrix
for n = 1:nWalls
    thisWall = n;
    thisNode = NodePos(n,:);
    incCount = 1;
    beta = wallImpFreq(n,:);

    %incoming nodes
    for m = 1:nWalls
        if n == m
            continue;
        else
            %image of node behind wall
            nodeImage = pointReflectionFromPlane(NodePos(m,:),Walls(thisWall,:));
            % perpendicular distance from source node to wall
            h = getDistance(NodePos(m,:), nodeImage)/2;
            %distance from source node to wall node
            d = getDistance(NodePos(m,:), thisNode);
            %angle between source node and wall node
            theta = acos(h/d);
            % incident pressure on wall (integrated over all theta and phi)
            pinc = 1i*exp(1i*waveNum.*h) - (waveNum*h).*igamma(0,-1i*waveNum.*h);

            outCount = 1;

            %outgoing nodes
            for k = 1:nWalls
                if k == m
                    continue;
                else
                     intersectionAtWall = lineIntersectsPlane(nodeImage, NodePos(k,:), Walls(thisWall,:));
                     %angle between node image and receiver
                     theta0 = getAngle(nodeImage - NodePos(m,:), nodeImage - NodePos(k,:));
                     gamma0 = cos(theta0);

                     R0 = (gamma0 - beta)./(gamma0 + beta);
                     % distance of receiver node from image source
                     r = getDistance(intersectionAtWall, nodeImage) + getDistance(intersectionAtWall, NodePos(k,:));
                     
                     % Q =  strength of virtual source at image position
%                      % if beta varies with the angle of incidence as beta = beta0 cos(theta)
%                      Q = (1-beta)./(1+beta);
%                      %otherwise (from eq 13)
                     rho = (1i*waveNum*r).*(((gamma0 + beta).^2)./(2*(1+gamma0.*beta)));
                     Q = R0 + ((1-R0).*sqrt(rho).*exp(rho).*whittakerW(-0.25,0.25,rho));
                     % reflected pressure
                     pref = Q.*exp(1i.*waveNum*r)./(waveNum*r);

                     S(n, incCount, outCount, :) = pref./pinc;
                     outCount = outCount+1;
                end
            end
            incCount = incCount+1;
        end
    end

    figure(n);
    Swall = reshape(S(n,:,:,:), [nWalls-1, nWalls-1, Nfreq]);

    % magnitude response should be flat
    plotImpulseResponseMatrix(freqs, 20*log10(abs(Swall)), 'xlim',[0,20000],...
    'ylim', [-10, max(max(max(20*log10(abs(Swall)))))], 'xlabel', 'Freq(Hz)', 'ylabel', 'Magnitude(dB)', 'stemFlag',0, 'save',0);
    print(['../figures/S_node=',num2str(n),'.eps'], '-depsc');

    
    % should be delta functions scaled by h/r and delayed by (r-h)/c
%     Swallfilt = real(ifft(fftshift(Swall,3),[],3));
%     plotImpulseResponseMatrix(1:Nfreq, Swallfilt(:,:,1:Nfreq), 'stemFlag',0);
end


%% NOTE 

% this just introduces 1/r scaling and r/c delay - that is equivalent to having 
% a unity scattering matrix and delay lines taking care of the rest


%% helper functions

function [image] = pointReflectionFromPlane(point, plane)
    
    a= plane(1); b = plane(2); c = plane(3); d = plane(4);
    x = point(1); y = point(2); z = point(3);
    
    % ﻿ # equation of line from (x1,y1,z1) to where it intersects plane
    %  # (x - x1) / a = (y - y1) / b = (z - z1) / c = k
    %  # replace x = ak + x1 etc in ax + by + cz + d = 0 to find k
        
        k = -(a * x + b * y + c * z + d) / (a ^ 2 + b ^ 2 + c ^ 2);
        
%         # where line from point intersects plane, is the midpoint between (x1,y1,z1)
%         # and its reflection
        
        xR = 2 * (a * k + x)  - x;
        yR = 2 * (b * k + y)  - y;
        zR = 2 * (c * k + z)  - z;

        image = [xR, yR, zR];

end


function [intersect] = lineIntersectsPlane(pointA, pointB, plane)
    
    intersect = zeros(1,3);
    xA = pointA(1); 
    yA = pointA(2); 
    zA = pointA(3);   
    xB = pointB(1); 
    yB = pointB(2); 
    zB = pointB(3);

    l = xB - xA;
    m = yB - yA;
    n = zB - zA;
    
    a = plane(1); b = plane(2); c = plane(3); d = plane(4);
    k = -(a*xA + b*yA + c*zA + d)/(a*l + b*m + c*n);
    
    intersect(1) = k*l + xA;
    intersect(2) = k*m + yA;
    intersect(3) = k*n + zA;
end



function [dist] = getDistance(point1, point2)
    dist = sqrt((point1(1) - point2(1))^2 + (point1(2) - point2(2))^2 + (point1(3) - point2(3))^2);
end

function angle = getAngle(vec1, vec2)
    %these are row vectors
    angle = acos((vec1*vec2.')/(norm(vec1)*norm(vec2)));
end


