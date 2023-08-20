stem = "raman10min";

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

% Now look at the maxZ
filename = stem + "_STEPS.txt";
data = load(char(filename), '-ascii');
maxZ = data(:,3); % column 3

% filename = stem + "-ENDPOS-PH4.txt";
filename = stem + "-ENDPOS-PH7.txt";
% filename = stem + "-ENDPOS-PH10.txt";

infile = filename;

% Initialize bins to zero
wLallbins  = zeros(10,10);
wLnoChangebins  = zeros(10,10);
wL1072bins  = zeros(10,10);
wL1430bins  = zeros(10,10);
wL1582bins  = zeros(10,10);
wL1702bins  = zeros(10,10);
wLotherbins = zeros(10,10);

global depths0105
global depths0106

global depths0205
global depths0206

global depths0305
global depths0306

global depths0405
global depths0406

global depths0505
global depths0506

global depths0605
global depths0606

global depths0705
global depths0706

global depths0805
global depths0806

global depths0905
global depths0906

global depths1005
global depths1006

close all

for iLoop = 1:length(infile)
    data = load(char(infile(iLoop)), '-ascii');
    x = data(:,1); % column 1
    y = data(:,2); % column 2
    z = data(:,3); % column 3
    myColor = data(:,4); % column 4
    
    % draw the SERS-active hydrogel
    r = [0.05];         % radii of each cyl
    cnt = [0,0.2];      % [x,y] center of each cyl
    height = 1;         % height of each cyl
    color = gold;        % color of each cyl
    nSides = 100;       % number of "sides" of the cyl
    plotCylinderWithCaps(r,cnt(:),height,nSides,color(:));
%     myTitle = sprintf('final position of photons for %s', char(infile(iLoop)));
%     title(myTitle);
    xlabel('x [cm]','FontSize', myTextFont);
    ylabel('y [cm]','FontSize', myTextFont);
    zlabel('z [cm]','FontSize', myTextFont);
    set(gca, 'ZDir','reverse', 'FontSize', myTextFont);
    
    % NEW 5/24/2022
    for i = 1:length(x)
        % use value of x and y to determine which bin this photon falls in
        [row, col] = determineBin(x(i),y(i));
        %fprintf('final position (%f,%f) goes into bin(%d,%d)\n',...
        %    x(i),y(i),row,col);
        
        % puts photons of all resultant wavelengths in the same bin
        if row > 0 && col > 0
            wLallbins(row,col) = wLallbins(row,col) + 1;
            switch myColor(i)
                case 0 % no change in wavelength
                    wLnoChangebins(row,col) = wLnoChangebins(row,col) + 1;            
                case 1 % change to 1072 cm-1 
                    wL1072bins(row,col) = wL1072bins(row,col) + 1;
                case 2 % change to 1430 cm-1
                    wL1430bins(row,col) = wL1430bins(row,col) + 1;            
                case 3 % change to 1582 cm-1    
                    wL1582bins(row,col) = wL1582bins(row,col) + 1;        
                case 4 % change to 1702 cm-1
                    wL1702bins(row,col) = wL1702bins(row,col) + 1;              
                case 5 % change to none of the above
                    wLotherbins(row,col) = wLotherbins(row,col) + 1; 
            end
            saveThisDepth(myColor(i), row, col, maxZ(i));
        end
    end
    
    fprintf('Photons of all wavelengths\n');
    for j = 1:10
        for i = 1:10
            fprintf('%d ', wLallbins(i,j));
        end
        fprintf('\n');
    end
    figure
    heatmap(wLallbins', 'Colormap', parula); % transpose to flip column major order to row
%    title('Photons of all wavelengths, vertical laser, pH7');
    title('Photons of all wavelengths, 45 degree laser, pH7');
    xlabel('x (mm)');
    ylabel('y (mm)');
    fprintf('Photons reflected\n');
    for j = 1:10
        for i = 1:10
            fprintf('%d ', wLnoChangebins(i,j));
        end
        fprintf('\n');
    end
    figure
    grid on
    heatmap(wLnoChangebins', 'Colormap', parula); % transpose to flip column major order to row
%    title('Photons reflected, vertical laser, pH7');
    title('Photons reflected, 45 degree laser, pH7');
    xlabel('x (mm)');
    ylabel('y (mm)');
    fprintf('Photons with 1072 cm-1 wavenumber\n');
    for j = 1:10
        for i = 1:10
            fprintf('%d ', wL1072bins(i,j));
        end
        fprintf('\n');
    end
    figure
    grid on
    heatmap(wL1072bins', 'Colormap', parula); % transpose to flip column major order to row
%    title('Photons with 1072 cm^-^1 wavenumber, vertical laser, pH7');
    title('Photons with 1072 cm^-^1 wavenumber, 45 degree laser, pH7');
    xlabel('x (mm)');
    ylabel('y (mm)');    
    fprintf('Photons with 1430 cm-1 wavenumber\n');
    for j = 1:10
        for i = 1:10
            fprintf('%d ', wL1430bins(i,j));
        end
        fprintf('\n');
    end
    figure
    grid on
    heatmap(wL1430bins', 'Colormap', parula); % transpose to flip column major order to row
%    title('Photons with 1430 cm^-^1 wavenumber, vertical laser, pH7');
    title('Photons with 1430 cm^-^1 wavenumber, 45 degree laser, pH7');
    xlabel('x (mm)');
    ylabel('y (mm)');    
    fprintf('Photons with 1582 cm-1 wavenumber\n');
    for j = 1:10
        for i = 1:10
            fprintf('%d ', wL1582bins(i,j));
        end
        fprintf('\n');
    end
    figure
    grid on
    heatmap(wL1582bins', 'Colormap', parula); % transpose to flip column major order to row
%    title('Photons with 1582 cm^-^1 wavenumber, vertical laser, pH7');
    title('Photons with 1582 cm^-^1 wavenumber, 45 degree laser, pH7');
    xlabel('x (mm)');
    ylabel('y (mm)');    
    fprintf('Photons with 1702 cm-1 wavenumber\n');
    for j = 1:10
        for i = 1:10
            fprintf('%d ', wL1702bins(i,j));
        end
        fprintf('\n');
    end
    figure
    grid on
    heatmap(wL1702bins', 'Colormap', parula); % transpose to flip column major order to row
%    title('Photons with 1702 cm^-^1 wavenumber, vertical laser, pH7');
    title('Photons with 1702 cm^-^1 wavenumber, 45 degree laser, pH7');   
    xlabel('x (mm)');
    ylabel('y (mm)');    
    fprintf('Photons with other wavenumbers\n');
    for j = 1:10
        for i = 1:10
            fprintf('%d ', wLotherbins(i,j));
        end
        fprintf('\n');
    end
    figure
    grid on
    heatmap(wLotherbins', 'Colormap', parula); % transpose to flip column major order to row
%    title('Photons with other wavenumbers, vertical laser, pH7');
    title('Photons with other wavenumbers, 45 degree laser, pH7');
    xlabel('x (mm)');
    ylabel('y (mm)');
    
    figure
    histogram(depths0105)
    title('Max depth achieved for photons exiting surface at bin(1,5)');
    xlabel('Max depth achieved (cm)');
    ylabel('Number of photons in bin');
    figure
    histogram(depths0106)
    title('Max depth achieved for photons exiting surface at bin(1,6)');
    xlabel('Max depth achieved (cm)');
    ylabel('Number of photons in bin');

    figure
    histogram(depths0205)
    title('Max depth achieved for photons exiting surface at bin(2,5)');
    xlabel('Max depth achieved (cm)');
    ylabel('Number of photons in bin');
    figure
    histogram(depths0206)
    title('Max depth achieved for photons exiting surface at bin(2,6)');
    xlabel('Max depth achieved (cm)');
    ylabel('Number of photons in bin');
    
    figure
    histogram(depths0305)
    title('Max depth achieved for photons exiting surface at bin(3,5)');
    xlabel('Max depth achieved (cm)');
    ylabel('Number of photons in bin');
    figure
    histogram(depths0306)
    title('Max depth achieved for photons exiting surface at bin(3,6)');
    xlabel('Max depth achieved (cm)');
    ylabel('Number of photons in bin');
    
    figure
    histogram(depths0405)
    title('Max depth achieved for photons exiting surface at bin(4,5)');
    xlabel('Max depth achieved (cm)');
    ylabel('Number of photons in bin');
    figure
    histogram(depths0406)
    title('Max depth achieved for photons exiting surface at bin(4,6)');
    xlabel('Max depth achieved (cm)');
    ylabel('Number of photons in bin');
    
    figure
    histogram(depths0505)
    title('Max depth achieved for photons exiting surface at bin(5,5)');
    xlabel('Max depth achieved (cm)');
    ylabel('Number of photons in bin');
    figure
    histogram(depths0506)
    title('Max depth achieved for photons exiting surface at bin(5,6)');
    xlabel('Max depth achieved (cm)');
    ylabel('Number of photons in bin');

    figure
    histogram(depths0605)
    title('Max depth achieved for photons exiting surface at bin(6,5)');
    xlabel('Max depth achieved (cm)');
    ylabel('Number of photons in bin');
    figure
    histogram(depths0606)
    title('Max depth achieved for photons exiting surface at bin(6,6)');
    xlabel('Max depth achieved (cm)');
    ylabel('Number of photons in bin');
    
    figure
    histogram(depths0705)
    title('Max depth achieved for photons exiting surface at bin(7,5)');
    xlabel('Max depth achieved (cm)');
    ylabel('Number of photons in bin');
    figure
    histogram(depths0706)
    title('Max depth achieved for photons exiting surface at bin(7,6)');
    xlabel('Max depth achieved (cm)');
    ylabel('Number of photons in bin');

    figure
    histogram(depths0805)
    title('Max depth achieved for photons exiting surface at bin(8,5)');
    xlabel('Max depth achieved (cm)');
    ylabel('Number of photons in bin');
    figure
    histogram(depths0806)
    title('Max depth achieved for photons exiting surface at bin(8,6)');
    xlabel('Max depth achieved (cm)');
    ylabel('Number of photons in bin');
    
    figure
    histogram(depths0905)
    title('Max depth achieved for photons exiting surface at bin(9,5)');
    xlabel('Max depth achieved (cm)');
    ylabel('Number of photons in bin');
    figure
    histogram(depths0906)
    title('Max depth achieved for photons exiting surface at bin(9,6)');
    xlabel('Max depth achieved (cm)');
    ylabel('Number of photons in bin');
    
    figure
    histogram(depths1005)
    title('Max depth achieved for photons exiting surface at bin(10,5)');
    xlabel('Max depth achieved (cm)');
    ylabel('Number of photons in bin');
    figure
    histogram(depths1006)
    title('Max depth achieved for photons exiting surface at bin(10,6)');
    xlabel('Max depth achieved (cm)');
    ylabel('Number of photons in bin');
end


function saveThisDepth(myColor, row, col, maxZ)

    global depths0105
    global depths0106

    global depths0205
    global depths0206

    global depths0305
    global depths0306

    global depths0405
    global depths0406

    global depths0505
    global depths0506
    
    global depths0605
    global depths0606
    
    global depths0705
    global depths0706

    global depths0805
    global depths0806
    
    global depths0905
    global depths0906
    
    global depths1005
    global depths1006
    
    % append the value maxZ to an array of depths for the square at 
    % (row,col), start with all colors together
    % maxZ1072 = [maxZ1072 maxZ(i)];
    switch row
        case 1
            switch col
%                 case 1
%                     depths0101 = [depths0101 maxZ];
%                 case 2
%                     depths0102 = [depths0102 maxZ];
%                 case 3
%                     depths0103 = [depths0103 maxZ];
%                 case 4
%                     depths0104 = [depths0104 maxZ];
                case 5
                    depths0105 = [depths0105 maxZ];
                case 6
                    depths0106 = [depths0106 maxZ];
%                 case 7
%                     depths0107 = [depths0107 maxZ];
%                 case 8
%                     depths0108 = [depths0108 maxZ];
%                 case 9
%                     depths0109 = [depths0109 maxZ];
%                 case 10
%                     depths0110 = [depths0110 maxZ];
            end
        case 2
            switch col
%                 case 1
%                     depths0201 = [depths0201 maxZ];
%                 case 2
%                     depths0202 = [depths0202 maxZ];
%                 case 3
%                     depths0203 = [depths0203 maxZ];
%                 case 4
%                     depths0204 = [depths0204 maxZ];
                case 5
                    depths0205 = [depths0205 maxZ];
                case 6
                    depths0206 = [depths0206 maxZ];
%                 case 7
%                     depths0207 = [depths0207 maxZ];
%                 case 8
%                     depths0208 = [depths0208 maxZ];
%                 case 9
%                     depths0209 = [depths0209 maxZ];
%                 case 10
%                     depths0210 = [depths0210 maxZ];
            end
        case 3
            switch col
%                 case 1
%                     depths0301 = [depths0301 maxZ];
%                 case 2
%                     depths0302 = [depths0302 maxZ];
%                 case 3
%                     depths0303 = [depths0303 maxZ];
%                 case 4
%                     depths0304 = [depths0304 maxZ];
                case 5
                    depths0305 = [depths0305 maxZ];
                case 6
                    depths0306 = [depths0306 maxZ];
%                 case 7
%                     depths0307 = [depths0307 maxZ];
%                 case 8
%                     depths0308 = [depths0308 maxZ];
%                 case 9
%                     depths0309 = [depths0309 maxZ];
%                 case 10
%                     depths0310 = [depths0310 maxZ];
            end
        case 4
            switch col
%                 case 1
%                     depths0401 = [depths0401 maxZ];
%                 case 2
%                     depths0402 = [depths0402 maxZ];
%                 case 3
%                     depths0403 = [depths0403 maxZ];
%                 case 4
%                     depths0404 = [depths0404 maxZ];
                case 5
                    depths0405 = [depths0405 maxZ];
                case 6
                    depths0406 = [depths0406 maxZ];
%                 case 7
%                     depths0407 = [depths0407 maxZ];
%                 case 8
%                     depths0408 = [depths0408 maxZ];
%                 case 9
%                     depths0409 = [depths0409 maxZ];
%                 case 10
%                     depths0410 = [depths0410 maxZ];
            end
        case 5
            switch col
%                 case 1
%                     depths0501 = [depths0501 maxZ];
%                 case 2
%                     depths0502 = [depths0502 maxZ];
%                 case 3
%                     depths0503 = [depths0503 maxZ];
%                 case 4
%                     depths0504 = [depths0504 maxZ];
                case 5
                    depths0505 = [depths0505 maxZ];
                case 6
                    depths0506 = [depths0506 maxZ];
%                 case 7
%                     depths0507 = [depths0507 maxZ];
%                 case 8
%                     depths0508 = [depths0508 maxZ];
%                 case 9
%                     depths0509 = [depths0509 maxZ];
%                 case 10
%                     depths0510 = [depths0510 maxZ];
            end
        case 6
            switch col
%                 case 1
%                     depths0601 = [depths0601 maxZ];
%                 case 2
%                     depths0602 = [depths0602 maxZ];
%                 case 3
%                     depths0603 = [depths0603 maxZ];
%                 case 4
%                     depths0604 = [depths0604 maxZ];
                case 5
                    depths0605 = [depths0605 maxZ];
                case 6
                    depths0606 = [depths0606 maxZ];
%                 case 7
%                     depths0607 = [depths0607 maxZ];
%                 case 8
%                     depths0608 = [depths0608 maxZ];
%                 case 9
%                     depths0609 = [depths0609 maxZ];
%                 case 10
%                     depths0610 = [depths0610 maxZ];
            end
        case 7
            switch col
%                 case 1
%                     depths0701 = [depths0701 maxZ];
%                 case 2
%                     depths0702 = [depths0702 maxZ];
%                 case 3
%                     depths0703 = [depths0703 maxZ];
%                 case 4
%                     depths0704 = [depths0704 maxZ];
                case 5
                    depths0705 = [depths0705 maxZ];
                case 6
                    depths0706 = [depths0706 maxZ];
%                 case 7
%                     depths0707 = [depths0707 maxZ];
%                 case 8
%                     depths0708 = [depths0708 maxZ];
%                 case 9
%                     depths0709 = [depths0709 maxZ];
%                 case 10
%                     depths0710 = [depths0710 maxZ];
            end
        case 8
            switch col
%                 case 1
%                     depths0801 = [depths0801 maxZ];
%                 case 2
%                     depths0802 = [depths0802 maxZ];
%                 case 3
%                     depths0803 = [depths0803 maxZ];
%                 case 4
%                     depths0804 = [depths0804 maxZ];
                case 5
                    depths0805 = [depths0805 maxZ];
                case 6
                    depths0806 = [depths0806 maxZ];
%                 case 7
%                     depths0807 = [depths0807 maxZ];
%                 case 8
%                     depths0808 = [depths0808 maxZ];
%                 case 9
%                     depths0809 = [depths0809 maxZ];
%                 case 10
%                     depths0810 = [depths0810 maxZ];
            end
        case 9
            switch col
%                 case 1
%                     depths0901 = [depths0901 maxZ];
%                 case 2
%                     depths0902 = [depths0902 maxZ];
%                 case 3
%                     depths0903 = [depths0903 maxZ];
%                 case 4
%                     depths0904 = [depths0904 maxZ];
                case 5
                    depths0905 = [depths0905 maxZ];
                case 6
                    depths0906 = [depths0906 maxZ];
%                 case 7
%                     depths0907 = [depths0907 maxZ];
%                 case 8
%                     depths0908 = [depths0908 maxZ];
%                 case 9
%                     depths0909 = [depths0909 maxZ];
%                 case 10
%                     depths0910 = [depths0910 maxZ];
            end
        case 10
            switch col
%                 case 1
%                     depths1001 = [depths1001 maxZ];
%                 case 2
%                     depths1002 = [depths1002 maxZ];
%                 case 3
%                     depths1003 = [depths1003 maxZ];
%                 case 4
%                     depths1004 = [depths1004 maxZ];
                case 5
                    depths1005 = [depths1005 maxZ];
                case 6
                    depths1006 = [depths1006 maxZ];
%                 case 7
%                     depths1007 = [depths1007 maxZ];
%                 case 8
%                     depths1008 = [depths1008 maxZ];
%                 case 9
%                     depths1009 = [depths1009 maxZ];
%                 case 10
%                     depths1010 = [depths1010 maxZ];
            end
    end
end

% 5/24/2022 
% This is confusing b/c Matlab is column major and I wrote this as row
% major. i is the horizontal dimension, where x ranges from 0 to 1 cm in
% steps of 0.1 cm. j is the vertical dimension, where y ranges from -0.5 to 0.5
% cm in steps of 0.1 cm.
function [i,j] = determineBin(x,y)
    i = 0;
    j = 0;
    
    % Choose i based on x value
    if x < 0.1 && x >= 0.
        i = 1;
    else 
        if x < 0.2 && x >= 0.1
            i = 2;
        else
            if x < 0.3 && x >= 0.2
                i = 3;
            else
                if x < 0.4 && x >= 0.3
                    i = 4;
                else
                    if x < 0.5 && x >= 0.4
                        i = 5;
                    else
                        if x < 0.6 && x >= 0.5
                            i = 6;
                        else
                            if x < 0.7 && x >= 0.6
                                i = 7;
                            else
                                if x < 0.8 && x >= 0.7
                                    i = 8;
                                else
                                    if x < 0.9 && x >= 0.8
                                        i = 9;
                                    else
                                        if x < 1.0 && x >= 0.9
                                            i = 10;
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end

    % Choose i based on x value
    if y < -0.4 && y >= -0.5
        j = 1;
    else 
        if y < -0.3 && y >= -0.4
            j = 2;
        else
            if y < -0.2 && y >= -0.3
                j = 3;
            else
                if y < -0.1 && y >= -0.2
                    j = 4;
                else
                    if y < 0.0 && y >= -0.1
                        j = 5;
                    else
                        if y < 0.1 && y >= 0.
                            j = 6;
                        else
                            if y < 0.2 && y >= 0.1
                                j = 7;
                            else
                                if y < 0.3 && y >= 0.2
                                    j = 8;
                                else
                                    if y < 0.4 && y >= 0.3
                                        j = 9;
                                    else
                                        if y < 0.5 && y >= 0.4
                                            j = 10;
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end

function [h1, h2, h3] = plotCylinderWithCaps(r,cnt,height,nSides,color)

[X,Y,Z] = cylinder(r,nSides);
%capColor = [0., 0., 0.]; % black
capColor = [0.9290, 0.6940, 0.1250];
X = X + cnt(1); 
Y = Y + cnt(2); 
Z = Z* height-1*height/2; % center at 0
%h1 = surf(X,Z,Y,'facecolor',color,'LineStyle','none');
%h1 = mesh(X,Z,Y,'edgecolor', 'k');
h1 = mesh(X,Z,Y,'edgecolor', color);
%h1 = surf(X,Z,Y,'LineStyle','none');
h2 = fill3(X(1,:),Z(1,:),Y(1,:),capColor); % Note: z<->y are flipped to obtain desired orientation
h3 = fill3(X(2,:),Z(2,:),Y(2,:),capColor); 
hold on;
end  