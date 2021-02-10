%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: christopher@todo.todo                               %
%                                                             %
% Monkey code to retrieve the Sx and Su matrices from a state %
% space model (A,B) for the "Reduced space" batch approach v2 %
% for solving dynamic QPs                                     %
%                                                             %
% A = System matrix                                           %                                            
% B = Gain matrix                                             %
% n = number of states                                        %
%                                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Sx, Su] = reduce_space(A, B, n)
    [An, Am] = size(A);
    [Bn, Bm] = size(B);
    
    assert(An == Am, 'A must be a square matrix.');
    assert(Bn >= Bm, 'Bn must be larger than or equal to Bm.');
    assert(n > 0, 'n must be a positive number.');
    
    Sx = cell(1, n);
    for i = 1:n
        Sx(1,i) = {A^i};
    end
    
    Su = cell(n, n);
    for i = 1:n
        for j = 1:n
            if j <= i
                Su(i, j) = {(A^(i-j))*B};
            else
                Su(i, j) = {zeros(Bn, Bm)};
            end
        end
    end
    
    Sx = cell2mat(Sx);
    Su = cell2mat(Su);
    
end
