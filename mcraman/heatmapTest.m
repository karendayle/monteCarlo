figure;
my_matrix = rand(3);
hm = heatmap(my_matrix, 'Colormap', parula(3), 'ColorbarVisible', 'on', 'XLabel', 'Time', 'YLabel', 'September')
ax = gca;
% ax.XData = [1 2 3]
properties(ax)
% Here we see the documented property XDisplayLabels, and changing this gives the desires result.
ax.XDisplayLabels = num2cell(0:2);
