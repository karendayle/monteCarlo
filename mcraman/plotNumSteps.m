infile = '../mcraman/skinvessel_STEPS.txt';
data = load(infile, '-ascii');
x = data(:,1); % column 1
y = data(:,2); % column 2

figure

histogram(y,length(x));
myLabelFont = 30;
title('Number of steps taken by each photon');

xlabel('Photon #', 'FontSize', myLabelFont);
ylabel('Number of steps', 'FontSize', myLabelFont);