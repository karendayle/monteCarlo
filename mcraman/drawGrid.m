figure

% vertical lines
x = 0.;
y = -0.5;
for i = 1:11
    for j = 1:11
        plot([x x], [y y+1], 'Color', 'k');
        hold on      
    end
    x = x + 0.1;
end
     
% horizontal lines
x = 0.;
y = -0.5;
for i = 1:11
    for j = 1:11
        plot([x x+1], [y y], 'Color', 'k');
        hold on     
    end
    y = y + 0.1;
end