function Y = step(X)
%HEAVISIDE    Step function.
%    HEAVISIDE(X) is 0 for X < 0, 1 for X > 0, and .5 for X == 0.
%    HEAVISIDE(X) is not a function in the strict sense.
%    See also DIRAC.

%   Copyright 1993-2008 The MathWorks, Inc.
%   $Revision: 1.1.6.3.2.1 $  $Date: 2008/07/17 04:41:04 $

Y = zeros(size(X));
Y(X > 0) = 1;
eng = symengine;
if strcmp(eng.kind,'maple')
    Y(X == 0) = nan;
else
    Y(X == 0) = 1;
end