function [h] = sodium_height(V)
%evaluates the height of sodium in 3M on the axis as a function of volume
%of sodium inside

r = 1.46;
p = [pi/3 -pi*r 0 4/3*pi*r^3-V];
r1 = roots(p);

answ = r1(r1>=0 & r1<=2*r);
h = answ(1);
end

