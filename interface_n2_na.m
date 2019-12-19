function [area] = interface_n2_na(h)
%INTERFACE_N2_NA returns surface area of the interface between sodium and
%nitrogen in m2
r = 1.46;
area = pi*(2*r-h)*(h);
end

