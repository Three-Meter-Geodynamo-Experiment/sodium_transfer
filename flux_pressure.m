function [flux] = flux_pressure(pressure)
%FLUX_PRESSURE returns the volume flow of sodium 
% pressure in psi, flux in L/s
% assuming there is a 1.5' d pipe, L=30 m, roughness 0.1mm
% scaling made by using data from http://www.pressure-drop.com/Online-Calculator/
flux = (pressure/1.2)^(1/1.92);
% for zero roughness 
% flux = (pressure/0.9)^(1/1.77);
end

