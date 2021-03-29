%% SCRIPT: PIPELINE_NON_LOCAL_MEANS
%
% Pipeline for non local means algorithm as described in [1].
%
% The code thus far is implemented in CPU.
%
% DEPENDENCIES
%
%   adaptNonLocalMeans
%   nonLocalMeans
%
  
%   clear all %#ok
%   close all

  %% PARAMETERS
  
  % input image
  pathImg   = 'test_nlm.mat';
  strImgVar = 'test_image';
  

  %% USEFUL FUNCTIONS

  % image normalizer
  normImg = @(I) (I - min(I(:))) ./ max(I(:) - min(I(:)));
  
  % minimum and maximum of input
  minmax = @(x) [min(x(:)) max(x(:))];
  
  
  %% (BEGIN)

  fprintf('...begin %s...\n',mfilename);  
  
  %% INPUT DATA
  
  fprintf('...loading input data...\n')
  
  ioImg = matfile( pathImg );
  I     = ioImg.(strImgVar);
  
  %% PREPROCESS
  
  fprintf(' - normalizing image...\n')
  I = normImg( I );
  
  figure('Name','Original Image');
  imagesc(I); axis image;
  colormap gray;
  
  %% NOISE
  
  fprintf(' - applying noise...\n')
  J = imnoise( I, noiseParams{:} );
  figure('Name','Noisy-Input Image');
  imagesc(J); axis image;
  colormap gray;
  
  
  %% NON LOCAL MEANS
  
  tic;
%   J   = ioImg.(strImgVar); %%%%%%
  If = nonLocalMeans( J, patchSize, filtSigma, patchSigma );
  toc
  
  
  %% PREPARE ADAPTIVE NLM INPUTS
  
  % distinct levels
  L = round( J .* (nLevel-1) );
  
  % adaptive sigma
  adFiltSigma = zeros( size( J ) );
  
  % sigma in each region (STD of each region)
  for i = 0 : (nLevel-1)
    
    adFiltSigma( L == i ) = std( J( L == i ) );
    
  end
  
  % visualize inputs adaptive NLM inputs
  figure('Name', 'Irregular search regions');
  imagesc(L); axis image;
  colormap parula;
  
  figure('Name', 'Adaptive sigma value');
  imagesc(adFiltSigma); axis image;
  colormap parula;
  
  
  %% ADAPTIVE NON LOCAL MEANS
    
  tic;
  Ia = adaptNonLocalMeans( J, L, patchSize, ...
                           adFiltSigma, patchSigma );
  toc
  
  
  %% VISUALIZE RESULT
  
  % get residuals
  res = cat(3, If-J, Ia-J);
  
  figure('Name', 'Filtered image');
  imagesc(If); axis image;
  colormap gray;
  
  figure('Name', 'Residual');
  imagesc( res(:,:,1) ); axis image;
  colormap gray;
  caxis( minmax( res ) )
  
  figure('Name', 'Adaptive filtered image');
  imagesc( Ia ); axis image;
  colormap gray;
  
  figure('Name', 'Adaptive residual');
  imagesc( res(:,:,2) ); axis image;
  colormap gray;
  caxis( minmax( res ) )
  
  %% (END)

  fprintf('...end %s...\n',mfilename);


%%------------------------------------------------------------
%
% AUTHORS
%
%   Dimitris Floros                         fcdimitr@auth.gr
%
% VERSION
%
%   0.1 - January 10, 2018
%
% CHANGELOG
%
%   0.1 (Jan 10, 2018) - Dimitris
%       * initial implementation
%
% ------------------------------------------------------------
