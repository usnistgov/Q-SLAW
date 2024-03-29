% Computes the deformation gradient tensor (Fij) and Jacobian of the
% deformation gradient tensor
%
%
%
% INPUTS :
% -------------------------------------------------------------------------
%
%   u_         = displacement field vector with rigid drift removed.
%                       Format: cell array, each containing a 2D matrix for each
%                       time point (components in x,y)
%                           unew{time}{1} = displacement in x-direction
%                           unew{time}{2} = displacement in y-direction
%                           unew{time}{3} = magnitude
%   dm         = meshgrid spacing
%   m2px       = meter to pixel conversion of original images in [x y]
%   type       = spatial differentiation kernel used for gradientN
%                options: 'fb', 'prewitt', 'sobel', 'scharr',
%                       'stencil', or 'optimal#' (# are odd numbers from 5
%                       to 19)
%                Suggested: 'optimal 9'
%
%
% OUTPUTS
% -------------------------------------------------------------------------
%   Fij     = deformation gradient tensor calculated on the mesh grid
%               with spacing dm.
%               Format: Fij{time}{2x2 deformation gradient tensor matrix}
%   J       = determinant of the deformation gradient tensor defined on the
%               mesh grid with spacing dm.
%               Format: J{time}
%
%
% NOTES
% -------------------------------------------------------------------------
% none
%
%
%%
function [Fij, J] = calculateFij_2d(varargin)

%Establish variables and inputs
[u,spacing,m2px,type] = parseInputs(varargin{:});
maxTime = length(u);
Fij = cell(maxTime,1);
J = cell(maxTime,1);

for i = 1:maxTime
    Fij{i} = funCalculateFij(u{i},spacing,m2px,type);
    J{i} = funCalculateJ(Fij{i});
end

end

function Fij = funCalculateFij(u,spacing,m2px,type)

% Calculate Displacment Gradient
Fij = cell(2,2);
for i = 1:2
    [Fij{i,1}, Fij{i,2}] = gradientN(u{i}*m2px(i),type);
end

%Calculate Deformation Gradient
for i = 1:2
    for j = 1:2
        Fij{i,j} = Fij{i,j}/(spacing*m2px(j));
    end
end

for i = 1:2, Fij{i,i} = Fij{i,i} + 1; end

end

function J = funCalculateJ(Fij)
%Calculate Jacobian of Deformation gradient
J =     Fij{1,1}.*Fij{2,2};
J = J - Fij{1,2}.*Fij{2,1};

end

function [u,spacing,m2px,type] = parseInputs(varargin)
u = varargin{1};
spacing = varargin{2};
m2px = varargin{3};

if length(m2px) < 2
    m2px = [m2px,m2px];
end

if length(varargin) < 4, type = 'optimal9';
else 
    type = varargin{4};
end

for i = 1:length(u)
    u{i} = cellfun(@double, u{i}, 'UniformOutput', 0);
end

end