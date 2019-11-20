function [h]=plot3d_errorbars(x, y, z, ez, varargin)

% create the standard 3d scatterplot
hold off;
h=scatter3(x, y, z,25);

hold on

% now draw the vertical errorbar for each point
for i=1:length(x)
        xV = [x(i); x(i)];
        yV = [y(i); y(i)];
       
        zMin = z(i) + ez(i);
        zMax = z(i) - ez(i);

        
        zB = [zMin, zMax];

        % draw error bars
        h=plot3(xV, yV, zB, '-k');
        
end