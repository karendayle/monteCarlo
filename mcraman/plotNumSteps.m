infile = 'skin3min_STEPS.txt';
data = load(infile, '-ascii');
x = data(:,1); % column 1
numSteps = data(:,2); % column 2
maxZ = data(:,3);

myLabelFont = 30;

figure
set(gca, 'FontSize', myLabelFont);
%histogram(y,length(x));
plot(x,numSteps);
%title('Number of steps taken by each photon');
xlabel('Photon #', 'FontSize', myLabelFont);
ylabel('Number of steps', 'FontSize', myLabelFont);

figure
set(gca, 'FontSize', myLabelFont);
%histogram(y,length(x));
plot(x,maxZ);
%title('Maximum depth achieved by each photon');
xlabel('Photon #', 'FontSize', myLabelFont);
ylabel('Maximum depth achieved (cm)', 'FontSize', myLabelFont);