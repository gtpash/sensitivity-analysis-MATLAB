function output = parameter_subset_selection(S, Toll)
% function [output] = parameter_subset_selection(S, Toll)
%
% Description: Run parameter subset selection from sensitivity matrix.
%
% Author(s) Paul R. Miles | November 14, 2018
%           Graham T. Pash
%           Lider S. Leon
%           Alexandros Solomou
%
% Input:
% - S: Sensitivity matrix (Nxp)
% - Toll: Algorithm stopping criteria (scalar)
%
% Output: Structure
% - 
%
% References:
% [1] Smith, R. C. (2013). Uncertainty quantification: theory, implementation,
%     and applications (Vol. 12). SIAM.
% [2] Leon, L. S., Smith, R. C., Oates, W. S., Miles, P. R. (2017, September). 
%     Identifiability and Active Subspace Analysis for a Polydomain Ferroelectric Phase Field Model. 
%     In ASME 2017 Conference on Smart Materials, Adaptive Structures and Intelligent Systems (pp. V002T03A022).

[~, p] = size(S);

i = 1; % Iterations
index_p = 1:p;
removed_index_p = zeros(1,p);

Diag_eigenvalues_sorted=zeros(1,p);% INITIALIZE THE ARRAY
if nargin < 2 || isempty(Toll)
    Toll = 1e-3;
end
% Toll=1e-3;
eigenvalues_cell = [];
eigenvector_cell = [];
while Diag_eigenvalues_sorted(1) < Toll || i==1
    
    fishM=S'*S;
    [eigenvectors,eigenvalues] = eig(fishM); 
    
    eigenvalues = eigenvalues/max(max(eigenvalues));
    
    Diag_eigenvalues = diag(eigenvalues);
    [Diag_eigenvalues_sorted,Ind] = sort(abs(Diag_eigenvalues));
    
    if Diag_eigenvalues_sorted(1) > Toll
    break;
    end
    target_eigenvector=eigenvectors(:,Ind(1));% THIS IS THE EIGEN VECTOR CORRESPONDING TO THE LOWEST MAGNITUDE EIGENVALUE
    [~,ud_P] = max(abs(target_eigenvector)); % ud_P(i) is the index of the component with the maximum magnitude OF THE TARGET EIGENVECTOR
    S(:,ud_P) = []; % REMOVE THAT COLUMN
    removed_index_p(i) = index_p(ud_P); % BASED ON THE INDEXES NUMBERING STORED IN INDEX_P, SAVE WHICH INDEX HAS BEEN REMOVED.
    index_p(ud_P) = [];
    
    
    eigenvalues_cell{i} = Diag_eigenvalues;
    eigenvector_cell{i} = eigenvectors;
    
    i = i + 1;
end

% assign output
output.eigenvalues_cell = eigenvalues_cell;
output.eigenvector_cell = eigenvector_cell;
output.removed_index_p = removed_index_p;
output.Diag_eigenvalues_sorted = Diag_eigenvalues_sorted;
output.Diag_eigenvalues = Diag_eigenvalues;
output.toleranace = Toll;
