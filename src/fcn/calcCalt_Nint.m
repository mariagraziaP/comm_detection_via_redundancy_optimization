function [CN,I_lvl,I_lvl_max] = calcCalt_Nint(COV,nth)

% ------------------------------------------------------------------------------
% COV:          covariance matrix of X
% nth:          threshold for subset samples
% CN:           neural complexity
% I_lvl:        average integration at level 1 to N
% I_lvl_max:    maximal integration at level 1 to N
% ------------------------------------------------------------------------------
% Returns the neural complexity for a system X.
% System dynamics is characterised by its covariance matrix COV.
% Computes complexity from the determinant of the covariance matrix.
% Olaf Sporns, Indiana University, 2003, 2021
% ------------------------------------------------------------------------------
N = size(COV,1);

% calculate I at level N
I_n = calcI_alt(COV);

% initialize all other levels - note: I(1) = 0
I_lvl = zeros(1,N-1);   I_lvl(1) = 0;
I_lvl_max = zeros(1,N-1);  I_lvl_max(1) = I_n/N;

% suppress warning resulting from large lvl-out-of-N
warning off MATLAB:nchoosek:LargeCoefficient

% loop over levels 2 to N-1
for lvl=2:N-1
    % how many subsets are there?
    n_subsets = nchoosek(N,lvl);
    % if below nth use all
    if (n_subsets < nth)
        subsets = nchoosek(1:N,lvl);
    end
    % if greater than nth create a subsample of size nth
    if (n_subsets >= nth)
        subsets = zeros(nth,lvl);
        for n=1:nth
            rp = randperm(N);
            subsets(n,:) = rp(1:lvl);
        end
    end
    % loop over individual subsets
    for s=1:size(subsets,1)
        COVs = COV(subsets(s,:),subsets(s,:));
        I_lvl(lvl) = I_lvl(lvl) + calcI_alt(COVs);
    end
    I_lvl(lvl) = I_lvl(lvl)./size(subsets,1);
    I_lvl_max(lvl) = (lvl/N)*I_n;
end
I_lvl(N) = I_n;
I_lvl_max(N) = I_n;

CN = sum(I_lvl_max-I_lvl);
