function [aX] = skew(a)
    [n, m] = size(a);
    
    assert(min(n, m) == 1 & max(n, m) == 3, ...
    "Input must be a vector of size 3x1 or 1x3")
    
    aX = [ 0   -a(3)  a(2);
          a(3)   0   -a(1);
         -a(2)  a(1)   0 ];
end