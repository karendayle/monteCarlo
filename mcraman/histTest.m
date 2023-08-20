data = 10*rand(1,100);
%x = histogram(data,'BinWidth',1);
% x = histogram(data,0.1,1);
labels = num2str(x.BinCounts');
xt = get(gca, 'XTick');
xVals = xt - ([0 diff(xt)]/2);
xVals(1) = [];
text(xVals, x.BinCounts, labels, 'HorizontalAlignment','center', 'VerticalAlignment','bottom')
% But first check your data to plot histogram there is something wrong about it:
% A = randi([90 100,100,1]); % look randi function and its inputs
% B = rand(50,1);
% C = A + B; % check A and B have the same sizes to ad