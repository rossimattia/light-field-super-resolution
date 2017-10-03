
% =========================================================================
% =========================================================================
%
% Author:
% Mattia Rossi (rossi.mattia@gmail.com)
% Signal Processing Laboratory 4 (LTS4)
% Ecole Polytechnique Federale de Lausanne (Switzerland)
%
% =========================================================================
% =========================================================================

function Z = merge(zSub, vRes, hRes, yRes, xRes, yResSub, xResSub, yCoord, xCoord, alpha)
% MERGE stiches all the input sub light fields into one single light field.
% Each sub light field view is weighted with a Tukey window and multiple
% estimates of the same pixel are average.
%
% INPUT:
% zSub - sub light fields in vectorized form.
% vRes - the vertical angular resolution of the full light field.
% hRes - the horizontal angular resolution of the full light field.
% yRes - the vertical spatial resolution of the full light field.
% xRes - the horizontal spatial resolution of the full light field.
% yResSub - the vertical spatial resolution of a sub light field.
% xResSub - the horizontal spatial resolution of a sub light field.
% yCoord - the vertical spatial coordinates of the sub light fields.
% xCoord - the horizontal spatial coordinates of the sub light fields.
% alpha - the Tukey window parameter.
%
% OUTPUT:
% Z - the full light field with spatial resolution (yRes x xRes).

% =========================================================================

% Number of views.
M = vRes * hRes;

% Number of sub light fields.
subNum = size(zSub, 2);

% Allocate a view stack and a normalization one, as some pixels may have
% multiple estimates due to the overlapping.
stack = zeros(yRes, xRes, M);
wei = zeros(yRes, xRes, M);

% Build the merging window.
% We build a window which is one pixel larger on each side, and then we
% crop it in order to have a (yRes x xRes) window with no zeros at the borders.
win = tukeywin(yResSub + 2, alpha) * tukeywin(xResSub + 2, alpha)';
win = win(2:(end - 1), 2:(end - 1));

% Filter each sub light fields with the merging window. This will limit
% the border effects when merging the different sub light fields.
aux = reshape(zSub, yResSub, xResSub, []);
aux = bsxfun(@times, aux, win);

% Merge the sub light fields, stored as a stack of view.
counter = 0;
for k = 1:1:subNum
    
    stack( ...
        yCoord(k):(yCoord(k) + yResSub - 1), ...
        xCoord(k):(xCoord(k) + xResSub - 1), ...
        :) = stack( ...
        yCoord(k):(yCoord(k) + yResSub - 1), ...
        xCoord(k):(xCoord(k) + xResSub - 1), ...
        :) + aux(:, :, (counter + 1):(counter + M));
    
    wei( ...
        yCoord(k):(yCoord(k) + yResSub - 1), ...
        xCoord(k):(xCoord(k) + xResSub - 1), ...
        :) = ...
        bsxfun( ...
        @plus, ...
        wei(yCoord(k):(yCoord(k) + yResSub - 1), xCoord(k):(xCoord(k) + xResSub - 1), :), ...
        win);
    
    counter = counter + M;
    
end

% Perform the normalization.
stack(:, :, :) = stack ./ wei;

% Arrange the light field in a "vRes x hRes" cell array.
Z = col2lf(stack(:), vRes, hRes, yRes, xRes, 1);

end

