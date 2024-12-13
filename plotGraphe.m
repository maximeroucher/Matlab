function []=plotGraphe(a,h)
%___________________________________________________________
% h : pas de la discrétisation spaciale
%___________________________________________________________
x=-a:h:a;
y=-a:h:a;
%
[x,y]=meshgrid(x,y);
z= J(x,y);
%
surf(x,y,z);
%
hold on
colormap hsv