function [xdot] = pendulumEOM(t,x,m,g,L)

xdot = [x(2), -g*sin(x(1))/L];

end