%
% SCRIPT: KNN_EXAMPLE
%
%   Sample code to read dataset, perform kNN, match labels, and
%   visualize kNN graph and adjacency matrix.
%
% MAT-FILES
% 
%   mnist_train.mat
%   mnist_train_svd.mat
%


%% CLEAN-UP

%clear
%close all


%% PARAMETERS

% data parameters
filename = 'train_testKNN';
varX     = 'train_X';
varL     = 'train_labels';

% kNN parameters
%kNbr = 1;


%% (BEGIN)

fprintf('\n *** Running %s ***\n\n',mfilename);


%% READ DATA

fprintf('Reading data...\n');
fprintf('   - using file %s\n', filename)

% IO data file
ioData = matfile( [filename '.mat'] );

% read variables
X = ioData.(varX);
L = ioData.(varL);

% number of points
nPoint = size( X, 1 );


%% KNN SEARCH

fprintf('Matlab knn search...'); tic;

[IDX, DIST] = knnsearch( X, X, 'k', kNbr+1 );

% drop first (self sources and targets are the same)
IDX  = IDX(:, 2:end);
DIST = DIST(:, 2:end);

fprintf('DONE in %.2f sec\n', toc);


%% MATCH LABELS

fprintf('Match labels...'); tic;

% labels of nearest neighbors
Lnn = L(IDX);

% find most frequent values in array
Mnn = mode( Lnn, 2 );

% find matches and incosistency
matches = ~(L - Mnn);

fprintf('DONE in %.2f sec\n', toc);


%% DISPLAY RESULTS

%fprintf('...displaying match results...\n');

fprintf('   - match percentage: %3.1f %%\n', ...
        nnz(matches) ./ nPoint * 100);


%% (END)

fprintf('\n *** End %s ***\n\n',mfilename);



%%------------------------------------------------------------
%
% AUTHORS
%
%   Dimitris Floros                         fcdimitr@auth.gr
%
% VERSION
%
%   0.1 - December 05, 2017
%
% CHANGELOG
%
%   0.1 (Dec 05, 2017) - Dimitris
%       * initial implementation
%
% ------------------------------------------------------------
