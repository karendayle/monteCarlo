% https://www.mathworks.com/matlabcentral/answers/433958-histogram-with-values-above-the-bars
% I would not depend upon the tick locations as Luna's solution does. Use the data stored in the histogram object itself instead.
% Generate some sample data and create the histogram.
data = 10*rand(1,100);
x = histogram(data,'BinWidth',1);
% Retrieve the locations of the edges and the heights of the bars.
E = x.BinEdges;
y = x.BinCounts;
% Find the X coordinates of the center of the bars. This takes the left edge of each bar and adds half the distance to the next bar edge.
%xloc = E(1:end-1)+diff(E)/2;
% Put the labels 1 unit above the top of the bar center.
%text(xloc, y+1, string(y))