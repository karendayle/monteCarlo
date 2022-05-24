close all;

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
infile = {'skin3min_STEPS.txt', 'skin3min-ENDPOS-PH4.txt', ...
    'skin3min-ENDPOS-PH7.txt', ...
    'skin3min-ENDPOS-PH10.txt'};

for iLoop = 1:length(infile)
    data = load(char(infile(iLoop)), '-ascii');
    switch iLoop
        case 1
            photon = data(:,1); % column 1
            steps = data(:,2); % column 2
            maxz = data(:,3); % column 3
        case 2
            x = data(:,1); % column 1
            y = data(:,2); % column 2
            z = data(:,3); % column 3
            pH4Color = data(:,4); % column 4
        case 3
            pH7Color = data(:,4); % column 4
        case 4
            pH10Color = data(:,4); % column 4
    end
end

% for iLoop = 1:length(infile)
%     switch iLoop
%         case 1
%             figure
%             set(gca, 'FontSize', myTextFont);
%             histogram(maxz);
%             %myTitle = sprintf('photon distribution for %s', char(infile(iLoop)));
%             %title(myTitle);
%             xlabel('Maximum depth achieved (cm)', 'FontSize', myTextFont);
%             ylabel('Number of photons', 'FontSize', myTextFont);
%         case {2,3,4}
%             figure
%             set(gca, 'FontSize', myTextFont);
%             if iLoop == 2
%                 histogram(pH4Color);
%             else
%                 if iLoop == 3
%                     histogram(pH7Color);
%                 else 
%                     histogram(pH10Color);
%                 end
%             end
%             %myTitle = sprintf('photon distribution for %s', char(infile(iLoop)));
%             %title(myTitle);
%             xlabel('Final wavenumber (cm-1)', 'FontSize', myTextFont);
%             xticklabels({'N/C','1082','1430','1584','1702','other'})
%             ylabel('Number of photons', 'FontSize', myTextFont);
%     end
% end

% Make a heat map of the max z vs final x and y
% figure
% %heatmap(x,y,maxz); complains of duplicate x values
% plot3(x,y,maxz);

% Repeat but only looking at the photons with final z = -0.01
k = 1;
for j = 1:10000 %length(x)
    if z(j)== -0.01
        % photon leaves out the top. This is the case we want
        cdata(k,1) = x(j);
        cdata(k,2) = y(j);
        %stepsNew(k) = steps(j);
        cdata(k,3) = maxz(j);
        k = k+1;
    end
end

clear x
clear y
clear maxz

% Make a heat map of the max z vs final x and y
figure
%heatmap('final x (cm)', 'final y (cm)', cdata);
heatmap(cdata);
%plot3(xNew,yNew,maxZNew,'+');
    

