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

for iLoop = 1:1 %length(infile)
    data = load(char(infile(iLoop)), '-ascii');
    x = data(:,1); % column 1
    y = data(:,2); % column 2
    z = data(:,3); % column 3
    myColor = data(:,4); % column 4

    figure
    grid on

    for i = 1:200 % takes a long time length(x)
        switch myColor(i)
            case 0
                plotColor = black; % no change in wavelength
            case 1
                plotColor = red; % change to 1082 cm-1
            case 2
                plotColor = green; % change to 1430 cm-1
            case 3
                plotColor = blue; % change to 1584 cm-1
            case 4
                plotColor = purple; % change to 1702 cm-1
            case 5
                plotColor = rust; % change to none of the above
        end
        plot3(x(i),y(i),z(i),'.','Color',plotColor);
        hold on;
        pause(0.001);
    end;

    x = 0.;
    y = 0.9;
    z = 1.0;
    deltaY = 0.4;
    text(x, y, z, 'no change', 'Color', black, 'FontSize', myTextFont);
    text(x, y, z, '_________', 'Color', black, 'FontSize', myTextFont);
    y = y - deltaY;
    text(x, y, z, '1082 cm-1', 'Color', red, 'FontSize', myTextFont);
    text(x, y, z, '_________', 'Color', red, 'FontSize', myTextFont);
    y = y - deltaY;
    text(x, y, z, '1430 cm-1', 'Color', green, 'FontSize', myTextFont);
    text(x, y, z, '_________', 'Color', green, 'FontSize', myTextFont);
    y = y - deltaY;
    text(x, y, z, '1584 cm-1', 'Color', blue, 'FontSize', myTextFont);
    text(x, y, z, '_________', 'Color', blue, 'FontSize', myTextFont);
    y = y - deltaY;
    text(x, y, z, '1702 cm-1', 'Color', purple, 'FontSize', myTextFont);
    text(x, y, z, '_________', 'Color', purple, 'FontSize', myTextFont);
    y = y - deltaY;
    text(x, y, z, 'other values', 'Color', rust, 'FontSize', myTextFont);
    text(x, y, z, '_____________', 'Color', rust, 'FontSize', myTextFont);
    y = y - deltaY;

    % draw the SERS-active hydrogel
    r = [0.05];         % radii of each cyl
    cnt = [0,0.2];      % [x,y] center of each cyl
    height = 1;         % height of each cyl
    color = red;        % color of each cyl
    nSides = 100;       % number of "sides" of the cyl
    plotCylinderWithCaps(r,cnt(:),height,nSides,color(:));
    myTitle = sprintf('final position of photons for %s', char(infile(iLoop)));
    title(myTitle);
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