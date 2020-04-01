% Colors:
blue =    [0.0000, 0.4470, 0.7410];
rust =    [0.8500, 0.3250, 0.0980];
gold =    [0.9290, 0.6940, 0.1250];
purple =  [0.4940, 0.1840, 0.5560];
green =   [0.4660, 0.6740, 0.1880];
ciel =    [0.3010, 0.7450, 0.9330];
cherry =  [0.6350, 0.0780, 0.1840];
red =     [1.0, 0.0, 0.0];
black =   [0.0, 0.0, 0.0];
myTextFont = 30;
infile = {'skin3min-ENDPOS-PH4.txt', 'skin3min-ENDPOS-PH7.txt', ...
    'skin3min-ENDPOS-PH10.txt'};

for iLoop = 1:length(infile)
    data = load(char(infile(iLoop)), '-ascii');
    x = data(:,1); % column 1
    y = data(:,2); % column 2
    z = data(:,3); % column 3
    myColor = data(:,4); % column 4
    
    figure
    set(gca, 'FontSize', myTextFont);
    histogram(myColor);
    %myTitle = sprintf('photon distribution for %s', char(infile(iLoop)));
    %title(myTitle);
    xlabel('Final wavenumber (cm-1)', 'FontSize', myTextFont);
    xticklabels({'N/C','1082','1430','1584','1702','other'})
    ylabel('Number of photons', 'FontSize', myTextFont);
end

function [h1, h2, h3] = plotCylinderWithCaps(r,cnt,height,nSides,color)
[X,Y,Z] = cylinder(r,nSides);
capColor = [0., 0., 0.]; % black
X = X + cnt(1); 
Y = Y + cnt(2); 
Z = Z* height-1*height/2; % center at 0
%h1 = surf(X,Z,Y,'facecolor',color,'LineStyle','none');
h1 = mesh(X,Z,Y,'edgecolor', 'k');
%h1 = surf(X,Z,Y,'LineStyle','none');
h2 = fill3(X(1,:),Z(1,:),Y(1,:),capColor); % Note: z<->y are flipped to obtain desired orientation
h3 = fill3(X(2,:),Z(2,:),Y(2,:),capColor); 
xlabel('x [cm]');
ylabel('y [cm]');
zlabel('z [cm]');
set(gca, 'ZDir','reverse')
hold on;
end  