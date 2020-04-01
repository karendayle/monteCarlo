infile = 'skin3min_STEPS.txt';
data = load(infile, '-ascii');
x = data(:,1); % column 1
y = data(:,2); % column 2
myLabelFont = 30;
figure
set(gca, 'FontSize', myLabelFont);
%histogram(y,length(x));
plot(x,y);
%title('Number of steps taken by each photon');
xlabel('Photon #', 'FontSize', myLabelFont);
ylabel('Number of steps', 'FontSize', myLabelFont);