% Script to generate a custom color map 

%Example:
% interptype = 'cubic';
% xn = round(linspace(1,256,255));
% colors = [1 .7 .4 .4 0;0 .4 .4 .7 1;1 .7 .4 .4 0]';
% xb = round(linspace(1,256,size(colors,1)));
% cmap = [interp1(xb,colors(:,1),xn, interptype);...
%         interp1(xb,colors(:,2),xn, interptype);...
%         interp1(xb,colors(:,3),xn, interptype)];

% Blue-Gray-Red

interptype = 'cubic';
xn = round(linspace(1,256,255));
colors = [ 0 0 1;.8 .8 .8;1 0 0]';
xb = round(linspace(1,256,size(colors,1)));
cmap = [interp1(xb,colors(:,1),xn, interptype);...
        interp1(xb,colors(:,2),xn, interptype);...
        interp1(xb,colors(:,3),xn, interptype)];