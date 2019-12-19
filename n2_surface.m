function [S] = n2_surface(h)
%N2_SURFACE evaluates surface area of the interface between 3M outer sphere
%and N2 bubble inside 
r = 1.46;
S = 2*pi*r*h;
end

