% Computes the lagrangian and infinitesimal strain tensors from the deformation 
% gradient tensor. 
%
% INPUTS : 
% -------------------------------------------------------------------------
%   Fij     = deformation gradient tensor calculated on the mesh grid
%               with spacing dm.
%               Format: Fij{time}{2x2 deformation gradient tensor matrix}
%                     
%
% OUTPUTS
% -------------------------------------------------------------------------
%   Eij     = Lagrangian strain tensor calculated on the mesh grid with
%               spacing dm.
%                   Format: Eij{time}{2x2 Lagrangian strain tensor
%                   matrix}
%   eij     = Infinitesimal strain tensor calculated on the mesh grid with
%               spacing dm. 
%                   Format: eij{time}{2x2 infinitesimal strain tensor
%                   matrix}
%
%
% NOTES
% -------------------------------------------------------------------------
% none
%%

function [Eij, eij] = calculateEij_2d(Fij)

%Set up variables
maxTime = length(Fij);
Eij = cell(maxTime,1);
eij = cell(maxTime,1);

for i = 1:maxTime
  Eij{i} = funCalculateLagrangianEij(Fij{i});
  eij{i} = funCalculateEij(Fij{i});
end


end

function eij = funCalculateEij(Fij)
%Calculate Infinitesimal strain
% A Bower, "solid mechanics", pg 22
eij = cell(2,2);

eij{1,1} = Fij{1,1} - 1;
eij{2,2} = Fij{2,2} - 1;

eij{1,2} = 0.5*(Fij{1,2} + Fij{2,1});

eij{2,1} = eij{1,2};

end

function Eij = funCalculateLagrangianEij(Fij)
%Calcualte Lagrangian Strain
% A. Bower,"Solid mechanics", pg 20

Eij{1,1} = 0.5*(Fij{1,1}.*Fij{1,1} + Fij{2,1}.*Fij{2,1} - 1);
Eij{2,2} = 0.5*(Fij{1,2}.*Fij{1,2} + Fij{2,2}.*Fij{2,2} - 1);
Eij{1,2} = 0.5*(Fij{1,1}.*Fij{1,2} + Fij{2,1}.*Fij{2,2} - 1);
Eij{2,1} = Eij{1,2}; 

end



