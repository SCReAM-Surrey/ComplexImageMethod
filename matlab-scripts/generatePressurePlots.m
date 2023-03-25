function [h1,h2] = generatePressurePlots(p, k, theta, varargin)
% polar plot (presure vs receiver angle) 
% and pressure vs wave number vs receiver angle plots
% p = complex pressire
% k - all wave numbers
% theta - all receiver angles
% optional input - plot title (this really helps)

pr = inputParser;
pr.KeepUnmatched = true;
addParameter(pr,'title',[]);
parse(pr,varargin{:});
tt = pr.Results.title;

% polar plot
idx = [2,5,10,50];
h1 = figure;
for j = 1:length(idx)
    polar_p = p(:,idx(j));
    subplot(length(idx),1,j);polarplot(theta, abs(polar_p));
    thetalim([0,90]);
    title(['k = ', num2str(k(idx(j)))]);
end
sgtitle(tt);

% pressure distribution

h2 = figure;
[K, THETA] = meshgrid(k, theta);
surf(K, THETA/pi, 20*log10(abs(p)));
set(get(gca,'ylabel')); %where angle is in degrees
ylabel('Angle of receiver from origin (norm. radians)'); xlabel('Wave number'); 
axis tight;
view(0,90);
colorbar;
title(tt);

end