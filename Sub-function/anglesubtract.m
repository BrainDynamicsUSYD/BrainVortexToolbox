function xydiff = anglesubtract(x, y, angleFlag)
% ANGLESUBTRACT Subtracts the matrix y from the equally sized matrix x and
% returns a value between -pi and +pi. This is valid if x and y both
% contain angular data in radians. If angleFlag is false, returns a
% non-angular subtraction.
%
% Rory Townsend, Oct 2017
% rory.townsend@sydney.edu.au

% Assume input is angular if unspecified
if nargin < 3
    angleFlag = 1;
end

if angleFlag == 1
    % METHOD 1: Modulo method
    xydiff = mod(x - y + pi, 2*pi) - pi;
    
    % METHOD 2: Trig method, works but is slower
    % xydiff2 = atan2(sin(x-y), cos(x-y));
    
else
    xydiff = x - y;
end

end