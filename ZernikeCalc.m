function [output1 output2] = ZernikeCalc( ZernikeList, Zdata, mask, ... 
                                          ZernikeDef, ShapeParameter, ...
                                          unitCircle, MMESorC)
%ZERNIKECALC Uses 'mask' region to fit circular (or other shape) Zernikes to surface data.
%
% VERSION:  2013-02-07   (YYYY-MM-DD)
%
% Fits circular, hexagonal, rectangular, square, elliptical, or annulus
% orthogonal Zernike polynomials to the surface data provided.  If no surface 
% data is provided then plots of these polynomial functions over the 
% region specified are produced.  Several different convesions (sign,
% normalization) can be explicitly specified.
%
% function [output1 output2] = ZernikeCalc( ZernikeList, Zdata, mask, ... 
%                                           ZernikeDef, ShapeParameter, ...
%                                           unitCircle, MMESorC)
% 
%  INPUT:
%    ZernikeList     = if size() = (1,k) then this is a Zernike j-ordering 
%                      list of numbers (j).
%                      if size() = (2,k) then this is Zernike (n, m) list.
%                      row 1 is the list of n values, row 2 is the list of
%                      m values.  m is positive or negative values to indicate
%                      whether sin() or cos() factor for the term.   See MMESorC
%                      parameter for the meaning of a minus m value.
%                      DEFAULT: [1:15]. First 15 j-ordered Zernikes.
%
%    Zdata           = If this is a single column then it contains the
%                      Zernike coefficients and we are to calculate the
%                      surface data.
%                      If this is more than a single column then it is a
%                      surface data and we are to calculate the 
%                      coefficients for the Zernike polynomials specified.
%                      The surface data must be a phase unwrapped surface data.
%                      DEFAULT: ones(15,1).
% 
%    mask            = Matrix (same size as surface data) indicating
%                      the individual pixels in the surface data that are
%                      to be used for fitting of the Zernike polynomials.
%                      '0' = don't use this pixel, '1' = use this pixel.
%                      If mask is a scalar number n, then an nxn mask matrix
%                      is created with an n diameter circle as the mask.
%                      DEFAULT: 100x100 with circluar mask.
% 
%    ZernikeDef      = One of 'FRINGE', 'ISO', 'WYANT', 'MAHAJAN', 'NOLL', 
%                      'B&W', 'STANDARD'
%                      'HEXAGON', 'HEXAGONR30', 'ELLIPSE',
%                      'RECTANGLE', 'SQUARE', 'ANNULUS'
%                      See table below for possible j-value ranges.
%                      NOTE: 'MAHAJAN' and 'NOLL' are programmed to be the same.
%                      NOTE: 'B&W' and 'STANDARD' are programmed to be the same.
%                      DEFAULT: 'FRINGE'.
% 
%    ShapeParameter  = For 'ELLIPSE' this is the aspect ratio (b).
%                      For 'RECTANGLE' this is half-width along the x-axis (a).
%                      For 'ANNULUS' this is the obscuration ratio (e).
%                      For all other shapes this is ignored.
%                      DEFAULT: 0, which is not valid for 'ELLIPSE' and 'RECTANGLE'.
%
%    unitCircle      = Gives the unit circle around the mask which the circular 
%                      Zernikes are defined on.  It is at most a 1x3 vector
%                      specifying centerRow, centerCol, radiusInPixels.  The
%                      permutations and ording allowed are as follows:
%                      []  empty so all 3 defaults are calculated.
%                      [radiusInPixels]  the default centerRow and centerCol
%                                        are calculated.
%                      [centerRow, centerCol] the default radiusInPixels is
%                                        calculated.
%                      [centerRow, centerCol, radiusInPixels]  no defaults.
%                      DEFAULT: The defaults for centerRow and centerCol are 
%                               calculated by  calculating the weighted
%                               average of row and column for the mask's 1's.
%                               The default radiusInPixels is calculated to
%                               be the largest distance from (centerRow,centerCol) 
%                               to all mask's pixels with a '1' value.
% 
%    MMESorC         = Either '-m=sin' or '-m=cos'.  Indicates, for (n, m) ordering
%                      what a minus m value represents, a sine factor or a
%                      cosine factor.
%                      DEFAULT: '-m=sin'. Ignored when j ordering is used.
% 
%  OUTPUT:
%    - When no output result is requested (ZernikeCalc(...)) this function
%      produces plots of the Zernike polynomials requested.
%      If Zdata is a single column specifying the coefficients then 
%          a plot of the sum of each of the Zernikes specified is produced.
%      If Zdata is an n x m matrix of surface data then 3 plots are
%      produced: 1) of the surface data, 2) of the fitted surface, 3) the
%      difference between the surface data and the fitted surface.
%    - When one output result is requested (results = ZernikeCalc(...)) then 1
%      of 2 possible results are returned:
%        if the input Zdata is a single column specifying the coefficients
%           to multiply the Zernike polynomials by then the result retruned
%           is the Zernike polynomial value matrices (across the mask area).
%        if the input Zdata is a matrix for which the Zernikes are to be
%           fit then the results returned is a column vector of the Zernike 
%           coefficients corrseponding to (and in the same order as) the 
%           polynomials identified by ZernikeList.
%    - If 2 output results are requested ([out1 out2] = ZernikeCalc(...)) then 
%            if the input for Zdata is a single column, giving the coefficients
%               to multiply the Zernike polynomials by, then out1 is the sum
%               of the Zernike polynomials requested, and out2 is a 3 dimensional
%               matrix of all the n Zernike polynomials requested [:,:,n].
%            if the input for Zdata is not a column vector then out1 is the 
%                3 dimensional data cube of the n Zernike polynomials requested [:,:,n]
%                and out2 are the coefficients used.
% 
%
% Examples: 
%
%   ZernikeCalc  
%   - Displays the first 15 Fringe Zernikes in 15 color plots.
%
%   ZernikeCalc([4 5 7 8], [0.4; -0.6; 1.2; 0.25])  
%   - Displays Fringe Zernikes Zj=4, Zj=5, Zj=7, Zj=8 multiplied by the 
%   scaling factors 0.4, -0.6, 1.2 and 0.25, respectively in 4 separate
%   color plots.
%
%   ZernikeCalc([4 5 7 8], [0.4; -0.6; 1.2; 0.25], [], 'standard')   
%   - Same as the last case only using standard Zernikes rather than Fringe
%   Zernikes.
%
%   ZernikeCalc([2 2; 2 0; 3 3; 3 1]', [0.4; -0.6; 1.2; 0.25], [], 'standard')   
%   - Same as last case now using Z(n,m) notation to specify which Zernikes
%   to use.
%
%   Let SD be an n x m matrix of surface data to which the specified 
%   Zernikes are to be fit.  Then
%   
%   coeff = ZernikeCalc([2 2; 2 0; 3 3; 3 1]', SD, [], 'standard')   
%    
%   returns a column consisting of the calculated fitting coefficients 
%   in coeff.  No plots are produced.
%
%   [DC, coeff] = ZernikeCalc([2 2; 2 0; 3 3; 3 1]', SD, [], 'standard')   
%    
%   returns a column consisting of the calculated fitting coefficients 
%   in coeff and a n x m x 4 data cube.  ('4' because 4 Zernikes were 
%   specified) Each DC(:, :, i) (i=1,2,3,4) is the ith specified Zernike  
%   fitted to the surface data SD across the (default) mask area.
%
%   [DC, coeff] = ZernikeCalc([4 5 7 8], SD, [], 'annulus', 0.25)
%   
%   This uses the annular Zernikes, with a central obscuration radius ratio of 0.25
%   to fit the surface data.  See Ref. 1 for details on noncircular Zernikes.
%   
%
% Circular Zernike polynomials are available in several optical design
% software packages, including Code V(R), OSLO(R), Zemax(R), etc.
%
%      Table 1 Software Conventions
%--------------------------------------------
%|INPUT PARAM |APERTURE |SOFTWARE|ORDER &   |
%|ZernikeDef  | SHAPE   |        | RANGE    |
%|------------|---------|--------|----------|			
%|'B&W',      |Circle   |CODE V  | (n, ±m), |
%|'STANDARD'  |         |        | j=1...   |
%|------------|---------|--------|----------|
%|'MAHAJAN',  |Circle   |ZEMAX   | (n, ±m), | 
%|'NOLL'      |         |        | j=1...   |
%|------------|---------|--------|----------|
%|'FRINGE'    |Circle   |CODE V, | j=1...37 |
%|            |         |ZEMAX   |          |
%|------------|---------|--------|----------|
%|'ISO'       |Circle   |        | j=0..35  |
%|------------|---------|--------|----------|
%|'WYANT'     |Circle   | OSLO   | j=0...36 |
%|------------|---------|--------|----------|
%|'HEXAGON'   |Hexagon  |        | j=1...45 |  
%|------------|---------|--------|----------|
%|'HEXAGONR30'|Heaxgon  |        | j=1...28 |
%|            |rotated  |        |          |
%|            |30 deg.  |        |          |
%|------------|---------|--------|----------|
%|'ELLIPSE'   |Ellipse  |        | j=1..15  |
%|------------|---------|--------|----------|
%|'RECTANGLE' |Rectangle|        | j=1...15 |
%|------------|---------|--------|----------|
%|'SQUARE'    |Square   |        | j=1...45 |
%|------------|---------|--------|----------|
%|'ANNULUS'   |Annulus  |        | j=1...35,|
%|            |         |        | j<>29,30,|
%|            |         |        |    31,32 |
%--------------------------------------------
% 
% Ref. 1:  Mahajan, V.N., G.-m. Dai, "Orthonormal polynomials in wavefront
%          analysis: analytical solution," J. Opt. Soc. Am. A, Vol. 24, No. 9
%          Sept. 2007.
%

%
% Updates:  2012-01-08  (YYYY-MM-DD)
%           RWG - Added default mask shapes for the different ZernikeDef
%                 input parameter values.
%
%           2012-01-08  (YYYY-MM-DD)
%           RWG - When no output requested ZernikeCalc will print all 
%                 Zernike polynomials specified.
%


%
% Code Copyrighted, 2011-2013 by Robert W. Gray and Joseph M. Howard.  All
% rights reserved.
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% For validation of the equations, uncomment the next 2 lines.
% Then, it doesn't matter what input parameters are specified.
%

%          validateZ();
%          return;

%
% end validation of equations.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Assign default values to input parameters that are empty.
  % Empty is '' for strings and [] for vectors and matrices.
  %
  
     if nargin == 0 || isempty(ZernikeList)
       ZernikeList = 1:15;    % default is first 15 fringe Zernikes.  
     end % if statement    
     
     if nargin <= 1 || isempty(Zdata)
       theSize = size(ZernikeList);  
       Zdata = ones(theSize(1,2),1);  % all coefficients are 1.
     end % if statement
         
            
     if nargin <= 3 || isempty(ZernikeDef)
       ZernikeDef = 'FRINGE';  
     end % if statement       
        
     % Convert to upper case
     ZernikeDef = upper(ZernikeDef);
     
     
     if nargin <= 4 || isempty(ShapeParameter)
       ShapeParameter = 0;  
     end % if statement
     
     
     if nargin <= 2 || isempty(mask)
       % make a default mask  
       
       defaultRows = 600;
       defaultCols = 600;
       
       theZdataSize = size(Zdata);
       if (theZdataSize(1,1) > 1) && (theZdataSize(1,2) > 1)
         defaultRows = theZdataSize(1,1);
         defaultCols = theZdataSize(1,2);
       end % if statement
       
       mask = makeDefaultMask(ZernikeDef, defaultRows, defaultCols, ShapeParameter);
       
     end % if no mask statement
     
     sm = size(mask);
     if (sm(1,1) == 1) && (sm(1,2) == 1)
       % if mask is a scalar n then create nxn matrix 
       % make a circular mask  
       defaultRows = mask(1,1);
       defaultCols = mask(1,1);
       
       mask = makeDefaultMask(ZernikeDef, defaultRows, defaultCols, ShapeParameter);
             
     end % end if statement    
     



     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     % Calculate default centerRow and centerCol values.
     % These might have been specified, but if they haven't 
     % use these values.
     
     radiusInPixels = 0;
     
     theSize = size(mask);
     numrows = theSize(1,1);
     numcols = theSize(1,2);
       
     % calculate center by weighted averages.
     sumMaskRow = sum(mask,2);
     sumMaskCol = sum(mask,1);
     
     sumMaskAll = sum(sum(mask));
     
     centerRow = sum(sumMaskRow .* ((1:numrows)')) / sumMaskAll;
     centerCol = sum(sumMaskCol .* (1:numcols)) / sumMaskAll;
     
     %
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       
     if nargin <= 5 || isempty(unitCircle)
       % The unit circle associated with the mask has not been specified.    
       % unitCircle is empty.
  
       unitCircle = [centerRow, centerCol];
       
     end % if statement 
     
     % At this point, unitCircle is not empty.
     theSize = size(unitCircle);
     nc = theSize(1,2);
     switch nc
         case 1
             % only the radiusInPixels has been specified
             radiusInPixels = unitCircle(1,1);
             
         case 2
             % the centerRow and centerCol have been specified
             % so can now calculate the radius in Pixels.
             centerRow = unitCircle(1,1);
             centerCol = unitCircle(1,2);
             
             % a matrix such that each element (r,c) has the value r-centerRow    
             rm = (((1:numrows)-centerRow)'*ones(1,numcols)).*mask;
             % a matrix such that each element (r,c) has the value c-centerCol
             cm = (ones(numrows,1)*((1:numcols)-centerCol)).*mask;
             % sqrt(rm.^2 + cm.^2) is a matrix such that (r,c) contains the distance
             % of (r,c) to the center (centerRow, centerCol).  
             radiusInPixels = max(max(sqrt(cm.^2 + rm.^2)));
           
         case 3
             % the centerRow, centerCol, radiusInPixels have been specified
             centerRow = unitCircle(1,1);
             centerCol = unitCircle(1,2);
             radiusInPixels = unitCircle(1,3);
          
         otherwise
             % error.
     end % switch statement

     if nargin <= 6 || isempty(MMESorC)
       MMESorC = '-m=sin';  
     end % if stateemnt
    
  %
  % end of section on default input assignments
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Input validation
     % Too much code is distracting, so put validation in separate function.
     validateInput(ZernikeList, Zdata, mask,  ...
                   ZernikeDef, ShapeParameter, ...
                   centerRow, centerCol, radiusInPixels,  MMESorC);
  %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % The Zernike polynomials are caluclated using (n, m) ordering
  % so need to calculate (n,m,sc) even when j ordering is specified.
  %  sc = 's' means sin() factor
  %  sc = 'c' means cos() factor
  %  sc = ' ' means m = 0 so has no sin() or cos() factor
  %
  
     theSize = size(ZernikeList);
     maxNumOfZernikeTerms = theSize(1,2);
     n = zeros(1, theSize(1,2));
     m = zeros(1, theSize(1,2));

     sc = ' ';
     if theSize(1,1) == 1
       % the ZernikeList is list of j order values.
       for k=1:theSize(1,2)
         % convert the j ordering of type ZernikeDef to (n,m,sc) of
         % same ZernikeDef.
         [n(k) m(k) sc(k)] = convertj(ZernikeList(1, k), ZernikeDef);  
       end % for statement    
     else
       % the ZernikeList is list of (n,m) pairs. 
       for k=1:theSize(1,2)
         % convert (n,m) to (n,m,sc) using MMESorC
         n(k) = ZernikeList(1, k);
         m(k) = abs(ZernikeList(2, k));
         sc(k) = ' ';
         switch MMESorC
             case '-m=sin'
                    if ZernikeList(2, k) < 0
                      sc(k) = 's';
                    end % if statement
                    if ZernikeList(2, k) > 0
                      sc(k) = 'c';  
                    end % if statement
             case '-m=cos'
                    if ZernikeList(2, k) < 0
                      sc(k) = 'c';
                    end % if statement
                    if ZernikeList(2, k) > 0
                      sc(k) = 's';  
                    end % if statement                 
         end % switch statement
         
       end % for k statement
        
     end % if j or (n,m) ordering    
  
  %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Preallocate some of the matrices.
  %
  
     % We'll need the size of the mask matrix (same as surface data).
     theSize = size(mask);
     rows = theSize(1,1);
     cols = theSize(1,2);
  
     % pre-allocate vectors and matrices.
     numOfMaskPixels = sum(sum(mask));

     zcoeff = zeros(maxNumOfZernikeTerms, 1); 
     xMatrices = zeros(rows, cols, maxNumOfZernikeTerms);
     xVectors = zeros(numOfMaskPixels, maxNumOfZernikeTerms);
     yVector = zeros(numOfMaskPixels, 1);
  
  %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  % Check if surface data is the input or if a coefficent vector is
  % the input.
  
  interferogramInput = true;
  
  theSize = size(Zdata);
  if theSize(1,2) == 1
    % There is one column in Zdata so it is a column vector of coefficents.  
    zcoeff = Zdata;
    interferogramInput = false;
  end % if statement
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %  for each pixel (matrix element) we need the distance from the 
  %  center (centerRow, centerCol).  So we make a matrix that has
  %  as element values the distance of that (row, col) element to 
  %  (centerRow, centerCol).
  %
     
    % a matrix such that each element (r,c) has the value r-centerRow    
    rm = ((1:rows)-centerRow)'*ones(1,cols);
    % a matrix such that each element (r,c) has the value c-centerCol
    cm = ones(rows,1)*((1:cols)-centerCol);
    
    % sqrt(rm.^2 + cm.^2)./radiusInPixels is a matrix such that (r,c) contains 
    % the normalized distance of (r,c) to the center (centerRow, centerCol).  
    % Then we reshape this into a column vector.
    rho = reshape(sqrt(cm.^2 + rm.^2)./radiusInPixels, rows*cols, 1);
    % atan2(rm, cm) is a matrix such that each element (r,c) contains
    % the angle (radians) from the centerRow axis to (r,c).  Then we 
    % reshape this into a vector.
    theta = reshape(atan2(rm, cm), rows*cols, 1);
    % reshape the mask into a vector.
    vmask = reshape(mask, rows*cols, 1);
   
    if interferogramInput
      % we have surface data so reshape it just like rho and theta.  
      yVector = reshape(Zdata, rows*cols, 1);
      % and remove all the elements that don't have a corresponding '1'
      % in the mask.
      yVector(vmask ==0) = [];
    end % if statement
    
    for i=1:maxNumOfZernikeTerms
      % calculate the ith Zernike polynomial value for each pixel
      % in the mask matrix (now a vector).
          
      switch ZernikeDef
         % Handle the special shapes.
         case 'HEXAGON'
            hldZ = ZHexagon(ZernikeList(1,i), rho, theta);  
         case 'HEXAGONR30'
            hldZ = ZHexagonR30(ZernikeList(1,i), rho, theta);  
         case 'ELLIPSE'
            hldZ = ZEllipse(ZernikeList(1,i), rho, theta, ShapeParameter);  
         case 'RECTANGLE'
            hldZ = ZRectangle(ZernikeList(1,i), rho, theta, ShapeParameter);  
         case 'SQUARE'
            hldZ = ZSquare(ZernikeList(1,i), rho, theta);  
         case 'ANNULUS'
            hldZ = ZAnnulus(ZernikeList(i), n(i), m(i), rho, theta, ShapeParameter);        
         otherwise
            % Otherwise, its a circle Zernike. 
            hldZ = calculateZernike(n(i), m(i), sc(i), rho, theta, ZernikeDef);
               
       end % switch ZernikeShape
  
       % reshape the column vector result into a (rows, cols) matrix.
       xMatrices(:, :, i) = reshape(hldZ,rows,cols);
       
       % remove the elements from the Zernike calculation that do not 
       % have a corresponding '1' in the mask.
       hldZ(vmask == 0) = [];
       % this is one of the Zernike polynomial results for each pixel
       % for which the mask is a '1'.
       xVectors(:, i) = hldZ; 
       
    end % for i statement
  
  %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Do the regression (least squares fit) only if the 
  % input Zdata is surface data.
  %
  
    if interferogramInput          
      % Use least squares fit to determine the coefficients.
      zcoeff = xVectors\yVector;
    end % if statement
  %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % multiply Zernike polynomial matrices by the calculated coefficients.
  % since we already have the Zernike matrices calculated, we do this here.
  %
    zm = zeros(rows, cols, maxNumOfZernikeTerms);
    for i=1:maxNumOfZernikeTerms
      % multiply the Zernike matrix by the corresponding coefficient
      % and by the mask to zero out pixels that we don't care about.
      zm(:,:,i) =  xMatrices(:, :, i) .* zcoeff(i) .* mask; 
    end % for statement    
  %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %
  % Change the output depending on input and output specified.
  %
    if nargout < 1
      % plot Zernike figures
      sizeZD = size(Zdata);
      
      if sizeZD(1,2) == 1
        % plot each of the Zernike polynomials  
        figs(zm);
        % plot the sum of the request Zernike Polynomials
 %       figs(sum(zm,3)); 
        return;
      end % if statement
      
      % plot 1) input surface data, 2) the fit results, 3) their difference
      
      toPlot(:,:,1) = Zdata;
      
      theFit = sum(zm, 3);
      toPlot(:,:,2) = theFit;
      
      theDiff = Zdata - theFit;
      toPlot(:,:,3) = theDiff;
      
      figs(toPlot, mask, 3, {'Input Data', 'Zernike Fit', 'Input Data - Zernike Fit'});

      return;
    end % if nargout < 1
    
    % At this point we know nargout >= 1
    
    if nargout < 2
      % This means nargout == 1.
      if interferogramInput
        % return calculated coefficients only
        output1 = zcoeff;
        return;
      end % if interferogramInput
        
      % Not an interferogram as Input.  Must be coefficients as input
      % so return only the interferogram.  They already know the 
      % coefficients.
        
      output1 = zm;
      return;
    end % if nargout < 2

    % At this point we know the output requested is 2 (or more).
    
    theZdataSize = size(Zdata);
    if theZdataSize(1,2) == 1
      output1 = sum(zm,3);
      output2 = zm;
      return;  
    end % if statement
    
    output1 = zm;
    output2 = zcoeff;
  
end % ZernikeCalc


function handle = figs(data,mask,numplotcols,titles,labeloption,labelprecision,handle)
%displays array data using imagesc function
%   function h = figs(data,mask,numplotcols,titles,labeloption,labelprecision)
%   INPUT:  data, 2, 3 or 4 dimensional stack of 2d arrays
%               2-d data of column vectors (i.e. non-square) are reshaped
%                   into 3-d stack of square figures
%               3-d data is subplotted with max of 6 cols and rows, with multiple
%                   figures created.
%               4-d data is plotted in a single figure with the 3rd and 4th
%                   dimension determining the subplot array size.
%           mask is optional, either one array for all, or a stack similar
%               in size to data.
%           numplotcols = set number of plot columns desired (default = 6 max)
%           titles = cell array of text titles for each plot (optional)
%           labeloption = 1, RMS data given
%                       = 2, Avg, RMS, and Peak-to-Valley given
%           labelprecision = number of digits in label (default = 4)
%           handle = handle of current figure to use
% 
%   OUTPUT:   handle = handles for each figure

%  To do:   1. add 5d capability
%           2. ensure mask is permuted along with data for 4d input

if nargin<1, handle = figure; return; end
if nargin<2, mask = []; end
if nargin<3, numplotcols = []; end
if nargin<4, titles = {}; end; if ~iscell(titles), titles = {}; end
if nargin<5, labeloption = []; end; if isempty(labeloption), labeloption = 1; end
if nargin<6, labelprecision = []; end; if isempty(labelprecision), labelprecision = 4; end
if nargin<7, handle = []; end

sizedata = size(data);

%extract data from structure if given as such
if isstruct(data)
    datastruct = data; clear data;
    for i=1:length(datastruct)
        if isfield(datastruct,'opd'), data(:,:,i) = datastruct(i).opd; end
        if isfield(datastruct,'mask'), mask(:,:,i) = datastruct(i).mask; end
        if isfield(datastruct,'psf'), data(:,:,i) = datastruct(i).psf; end
    end
end
% reshape data if 2d column or row input and square (e.g. size(col) = 100^2)
if size(data,1) ~= size(data,2) % if non-square input matrix
    sqrtval = sqrt(length(data(:,1)));
    if mod(length(data(:,1)),10000) == 0 %if data is integer lengths of 100x100
        disp('Data appears to be integer number of 100x100 matrices in a column vector.')
        disp('Reshaping data into 2-D stack.')
        num100x100 = length(data(:,1))/10000;
        if num100x100 == 1
            data = reshape(data,100,100,size(data,2)); %3d data
        else
            data = reshape(data,100,100,num100x100,size(data,2));  %4d data
        end
    elseif mod(length(data(:,1)), sqrtval) == 0 % if data square
        disp('Data appears to be square in column format, reshaping into 2d.');
        data = reshape(data,sqrtval,sqrtval,size(data,2));
    else
        sqrtval = sqrt(length(data(:,2)));
        if mod(length(data(:,1)), sqrtval) == 0 % if data square
            disp('Data appears to be square in row format, reshaping into 2d.');
            data = reshape(data',sqrtval,sqrtval,size(data',2));
        end
    end
end

% determine number of plots, and convert 5d and 4d data to 3d stack
if ndims(data)==5
    [s1,s2,s3,s4,s5] = size(data);
    data = reshape( permute(data,[1 2 5 4 3]),s1,s2,s3*s4*s5);
end
if ndims(data)==4
    [s1,s2,s3,s4] = size(data);
    numplots = s3*s4;
    data = reshape( permute(data,[1 2 4 3]),s1,s2,s3*s4);
    data4d = s4;
else
    numplots = size(data,3);
    data4d = 0;
end


%create mask data if needed
if isempty(mask), mask = ~isnan(data) & data~=0; end 

% determine subplot format: rows and cols 
if isempty(numplotcols)
    if data4d, numplotcols = data4d; %let 4d data size determine rows and cols 
    elseif numplots==1, numplotcols=1; 
    elseif numplots<5, numplotcols = 2;
    elseif numplots<7, numplotcols = 3;
    elseif numplots<9, numplotcols = 4;
    elseif numplots<11, numplotcols = 5;
    else numplotcols = 6;
    end
end
numplotrows = ceil(numplots/numplotcols);
if numplotrows>6 && data4d<1, numplotrows=6; end
if numplots>numplotrows*numplotcols
    numfigs = ceil(numplots/(numplotrows*numplotcols)); %number of figures
else
    numfigs = 1;
end

if ~iscell(titles) %array of text given
    for i=1:size(titles,1)
        t{i} = titles(i,:);
    end
    titles = t;
end
numtitles = length(titles(:));

disp(['Total number of figures     = ' int2str(numfigs)]);
disp(['Number of plots per figure  = ' int2str(numplotrows*numplotcols)]);
disp(['Number of plots rows        = ' int2str(numplotrows)]);
disp(['Number of plots cols        = ' int2str(numplotcols)]);

% plot data
for k=1:numfigs
    scrsz = get(0,'ScreenSize');
    if numfigs>1 || isempty(handle)
        handle = figure('Position',[scrsz(3)/4 scrsz(4)/4 0.6*scrsz(3) 0.6*scrsz(4)],'color',[1 1 1]);
    else figure(handle);
    end
        set(gcf,'Name',inputname(1));
    for i=1:numplotrows
        for j=1:numplotcols
            plotnum = (i-1)*numplotcols+j+(k-1)*numplotrows*numplotcols; %data to plot
            subplotnum = (i-1)*numplotcols+j; %subplot location to plot
            if plotnum<numplots+1,
                plotdata = data(:,:,plotnum)*1; % *1 converts logical to doubles
                if size(mask,3)>1, maskdata = mask(:,:,plotnum); else maskdata = mask(:,:,1); end
                subplot(numplotrows,numplotcols,subplotnum);
                plotdata(maskdata==0)=nan;
                imagesc(plotdata);
                    axis square; 
                    axis xy;
                alpha(1*maskdata);
                set(gca,'XTick',[]);
                set(gca,'YTick',[]);
                set(gca,'ZTick',[]);
                if numtitles<1
                    title(int2str(plotnum));
                elseif plotnum>numtitles %do nothing if partial title list given
                else title(titles{plotnum}); 
                end 
                put_xlabel(plotdata,maskdata,labeloption,labelprecision); %subfunction included below
            end
            if nargin>3, if length(titles)>=plotnum, title(titles(plotnum)); end; end
        end
    end
end
end % function figs

% X Label for each plot
function put_xlabel(plotdata,maskdata,labeloption,labelprecision)
if nargin<2
    maskdata = ~isnan(maskdata) & data~=0;
    datavals = plotdata(maskdata);
else
    if ~islogical(maskdata), maskdata = logical(maskdata); end
    datavals = plotdata(maskdata);
end
if ~isempty(datavals) %only analyze statistics on a real data set
    stdval = std(plotdata(maskdata));
 %   pk2val = pv(plotdata,maskdata);
    maxval = max(datavals);
    minval = min(datavals);
    meanval = mean(datavals);
    p95 = sort(datavals);
    p95 = p95(ceil(.95*length(datavals)));
else
    stdval = 0; pk2val = 0; maxval = 0; minval = 0; meanval = 0; p95 = 0;
end
if labeloption==1
    xlabel(['RMS = ' num2str(stdval,labelprecision)]);
elseif labeloption==2
    xlabel({['AVG = ' num2str(meanval,labelprecision)];['RMS = ' num2str(stdval,labelprecision)];['PV =  ' num2str(pk2val,labelprecision)]});
elseif labeloption==3
    xlabel(['RMS = ' num2str(stdval*1e9,labelprecision) ' nm']);    
end

end % function put_xlabel




function result = calculateZernike(n, m, sc, rho, theta, ZernikeDef)
%
%  Calculates the Zernike polynomial value for the given pixel location.
%
%   INPUT:
%     n           = radial order. 
%     m           = azimuthal order.
%     sc          = ' ' for m = 0,  = 's' for sin() term,  = 'c' for cos() term.
%     rho         = normalized radial distance to pixel location.
%     theta       = angle from the x-axis in radians of pixel location.
%     ZernikeDef  = One of 'MAHAJAN', 'NOLL', 'FRINGE', 'ISO', 'WYANT',
%                          'B&W', 'CIRCLE', 'HEXAGON', 'HEXAGONR30',
%                          'RECTANGLE', 'SQUARE', 'ANNULUS'
%
%   OUTPUT
%     results = The Zernike polynomial (n,m,sc) value for the given pixel (rho, theta).
%

   % calculate radial part Rnm
   Rnm = zeros(size(rho));
   for s=0:(n-m)/2
     numerator = (-1)^s * factorial(n-s);
     denominator = factorial(s)*factorial((n+m)/2-s)*factorial((n-m)/2-s);
     Rnm = Rnm + (numerator / denominator) * rho.^(n-2*s);  
   end % for s statement    
  
   % 3 cases.  sc=' ', 's', 'c'.
   theFactor = 1;
   switch sc
     case ' '
       % means m=0  
       theFactor = sqrt(n+1); 
       result = Rnm;  
     case 's'
       theFactor = sqrt(2*(n+1));  
       result = Rnm .* sin(m*theta);  
     case 'c'  
       theFactor = sqrt(2*(n+1));  
       result = Rnm .* cos(m*theta);  
   end % switch sc  
   
   switch ZernikeDef
       case {'MAHAJAN', 'NOLL', ...
             'HEXAGON', 'HEXAGONR30', 'RECTANGLE', 'SQUARE', 'ANNULUS'}
         result = theFactor * result;  
   end % switch
   

end % functon calculateZernike



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% The following functions are used to calculate (n,m,sc) from a j-order
% number.
%
%

function [n, m, sc] = convertj(j, ztype)
%CONVERTJ Convert ordering parameter j to n, m, sc parameters.
%
% function [n, m, sc] = convertj(j, ztype)
%
%   INPUT:
%     j     = Zernike polynomial order number.  j >= 0 and integer.
%             For 'FRINGE' 1 <= j <= 37.
%             For 'ISO' 0 <= j <= 35.
%     ztype = type of odering: 'FRINGE', 'ISO', 'WYANT', 'B&W', 
%                              'STANDARD', 'MAHAJAN', 'NOLL',
%                              'HEXAGON', 'HEXAGONR30', 'ELLIPSE',
%                              'RECTANGLE', 'SQUARE', 'ANNULUS'
%
%   OUTPUT:
%     n  = radial order. 
%     m  = azimuthal order
%     sc = indicate whether sine or cosine or neither factor
%          's' means sin() factor
%          'c' means cos() factor
%          ' ' means m = 0 so has no sin() or cos() factor
%

  %
  % DANGER: Validation of input is NOT performed. 
  %         This function should only be called from ZernikeCalc 
  %         which does the input validation. 
  %

  
  switch ztype
    case 'FRINGE' 
        [n,m,sc] = fringe(j);
        
    case {'ISO', 'WYANT'}
        [n,m,sc] = fringe(j+1);
        
    case {'B&W', 'STANDARD'}
        [n,m,sc] = bw(j);
 
    case {'MAHAJAN','NOLL', ...
          'HEXAGON', 'HEXAGONR30', 'ELLIPSE', ...
          'RECTANGLE', 'SQUARE', 'ANNULUS'}
        [n,m,sc] = mahajan(j);
        
  end % switch ztype

end % function convertj


function [n, m, sc] = fringe(j)
%
% Note that 'FRINGE', 'ISO', 'WYANT' use this function
% to assign (n, m) pairs to the j values.
%
% sc = 's' = sin(),  sc = 'c' = cos(),  sc = ' ' for neither.
%

  switch j
    case 1
      n = 0; m = 0; sc = ' ';
    case 2
      n = 1; m = 1; sc = 'c';  
    case 3
      n = 1; m = 1; sc = 's';
    case 4
      n = 2; m = 0; sc = ' ';
    case 5
      n = 2; m = 2; sc = 'c';
    case 6
      n = 2; m = 2; sc = 's';
    case 7
      n = 3; m = 1; sc = 'c';
    case 8
      n = 3; m = 1; sc = 's';
    case 9
      n = 4; m = 0; sc = ' ';
    case 10
      n = 3; m = 3; sc = 'c';
    case 11
      n = 3; m = 3; sc = 's';
    case 12
      n = 4; m = 2; sc = 'c';
    case 13
      n = 4; m = 2; sc = 's';
    case 14
      n = 5; m = 1; sc = 'c';
    case 15
      n = 5; m = 1; sc = 's';
    case 16
      n = 6; m = 0; sc = ' ';
    case 17
      n = 4; m = 4; sc = 'c';
    case 18
      n = 4; m = 4; sc = 's';
    case 19
      n = 5; m = 3; sc = 'c';
    case 20
      n = 5; m = 3; sc = 's';
    case 21
      n = 6; m = 2; sc = 'c';
    case 22
      n = 6; m = 2; sc = 's';
    case 23
      n = 7; m = 1; sc = 'c';
    case 24
      n = 7; m = 1; sc = 's';
    case 25
      n = 8; m = 0; sc = ' ';
    case 26
      n = 5; m = 5; sc = 'c';
    case 27
      n = 5; m = 5; sc = 's';
    case 28
      n = 6; m = 4; sc = 'c';
    case 29
      n = 6; m = 4; sc = 's';
    case 30
      n = 7; m = 3; sc = 'c';
    case 31
      n = 7; m = 3; sc = 's';
    case 32
      n = 8; m = 2; sc = 'c';
    case 33
      n = 8; m = 2; sc = 's';
    case 34
      n = 9; m = 1; sc = 'c';
    case 35
      n = 9; m = 1; sc = 's';
    case 36
      n = 10; m = 0; sc = ' ';
    case 37
      n = 12; m = 0; sc = ' ';

  end % switch j
  
end % function fringe


function [n, m, sc] = bw(j)
  % sc = 's' = sin(),  sc = 'c' = cos(),  sc = ' ' for neither.
  % sc = 1 = sin(),  sc = 2 = cos().
 
  % calculate the n value
  n1 = (-1 + sqrt(1 + 8 * j)) / 2;
  n = floor(n1);
  if n1 == n
   n = n - 1;   
  end % if statement
  
  % calculate the m value
  k = (n+1)*(n+2)/2;
  d = k - j;
  m1 = n - 2*d;
  m = abs(m1);
  
  % calculate the sc value
  sc = ' ';
  if m1 > 0 
    sc = 's';
  end % if statement
  if m1 < 0 
    sc = 'c';
  end % if statement
      
end % function bw



function [n, m, sc] = mahajan(j)
  % sc = 's' = sin(),  sc = 'c' = cos(),  sc = ' ' for neither.
  % sc = 1 = sin(),  sc = 2 = cos().
 
  % calculate the n value
  n1 = (-1 + sqrt(1 + 8 * j)) / 2;
  n = floor(n1);
  if n1 == n
   n = n - 1;   
  end % if statement
  
  % calculate the m value
  k = (n+1)*(n+2)/2;
  m = n - 2 * floor((k - j)/2);
  
  % calculate the sc value
  sc = ' ';
  if (m ~= 0) && (mod(j,2) ~= 0)
    % m ~= 0 and j odd  
    sc = 's';
  end % if statement
  if (m ~= 0) && (mod(j,2) == 0)
    % m ~= 0 and j even  
    sc = 'c';
  end % if statement
      
end % function mahajan

%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function ok = validateInput(ZernikeList, Zdata, mask, ...
                            ZernikeDef, ShapeParameter, ...
                            centerRow, centerCol, radiusInPixels, MMESorC)
%
%  Validates the input.  
%  See above function definition for definition of input. 
%
%   OUTPUT
%     ok  = true.
%

  ok = true;
  
  % Check validity of input parameters.
  
  theIOCSize = size(Zdata);
  theMaskSize = size(mask);
  
  if (theIOCSize(1,2) > 1) && (sum(theIOCSize == theMaskSize) ~= 2)
    ME = MException('VerifyData:InvalidData', ...
          'ZernikeCalc: Surface data and mask matrices are not the same size.');
    throwAsCaller(ME);       
  end % if statement    
  
  if centerRow < 1
    ME = MException('VerifyData:InvalidData', ...
          'ZernikeCalc: centerRow must be positive.');
    throwAsCaller(ME);        
  end % if statement    
  
  if centerCol < 1
    ME = MException('VerifyData:InvalidData', ...
          'ZernikeCalc: centerCol must be positive.');
    throwAsCaller(ME);        
  end % if statement
  
  if radiusInPixels < 1
    ME = MException('VerifyData:InvalidData', ...
          'ZernikeCalc: radiusInPixels must be positive.');
    throwAsCaller(ME);        
  end % if statement
  
  hlda = mask == 0;
  hldb = mask == 1;
  if sum(sum(hlda + hldb)) ~= (theMaskSize(1,1)*theMaskSize(1,2))
     ME = MException('VerifyData:InvalidData', ...
                     'ZernikeCalc: mask matrix must contain 0 or 1 only.');
     throwAsCaller(ME);    
  end % if statement    
   
  
  % Now for the fun stuff: Validating the Zernike ordering.
  
  switch ZernikeDef
      case {'FRINGE', 'ISO', 'WYANT' 'MAHAJAN', 'NOLL', 'B&W', 'STANDARD', ...
            'HEXAGON', 'HEXAGONR30', 'ELLIPSE', 'RECTANGLE', ...
            'SQUARE', 'ANNULUS'}
        % These are the valid values.  
      otherwise
        % ZernikeType is not valid
        ME = MException('VerifyData:InvalidData', ...
                     'ZernikeCalc: ZernikeDef is not valid.');
        throwAsCaller(ME); 
  end % switch statement        
      

  theSize = size(ZernikeList);
  rows = theSize(1,1);
  cols = theSize(1,2);
    
  hld1 = sum(sum(zeros(rows, cols) == (abs(ZernikeList) - floor(abs(ZernikeList)))));
  if hld1 ~= rows*cols
    % a number in ZernikeList is not an integer
    ME = MException('VerifyData:InvalidData', ...
                     'ZernikeCalc: ZernikeList must contain only integers.');
    throwAsCaller(ME);    
  end % if statement
  
  if rows == 1
    % j ordering
    if sum(abs(ZernikeList(1, :)) ~= ZernikeList(1, :)) ~= 0
      ME = MException('VerifyData:InvalidData', ...
                      'ZernikeCalc: ZernikeList j values must be positive or 0.');
      throwAsCaller(ME);            
    end % if statement    
    
    switch ZernikeDef
      case 'FRINGE'
        if sum(ZernikeList(1, :) < 1) ~= 0
          ME = MException('VerifyData:InvalidData', ...
                          'ZernikeCalc: FRINGE order value j must be greater than 0.');
          throwAsCaller(ME);             
        end % if statement 
        
        if sum(ZernikeList(1, :) > 37) ~= 0
          ME = MException('VerifyData:InvalidData', ...
                          'ZernikeCalc: FRINGE order value j must be less than 38.');
          throwAsCaller(ME);             
        end % if statement 
        
      case 'ISO' 
        if sum(ZernikeList(1, :) > 35) ~= 0
          ME = MException('VerifyData:InvalidData', ...
                          'ZernikeCalc: ISO order value j out of range.');
          throwAsCaller(ME);             
        end % if statement  
        
      case 'WYANT' 
        if sum(ZernikeList(1, :) > 36) ~= 0
          ME = MException('VerifyData:InvalidData', ...
                          'ZernikeCalc: WYANT order value j out of range.');
          throwAsCaller(ME);             
        end % if statement  
        
      case {'MAHAJAN', 'NOLL', 'B&W', 'STANDARD', ...
            'HEXAGON', 'HEXAGONR30', 'ELLIPSE', 'RECTANGLE', ...
            'SQUARE', 'ANNULUS'}
        if sum(ZernikeList(1, :) == 0) ~= 0
          ME = MException('VerifyData:InvalidData', ...
                          'ZernikeCalc: Zernike order value j can not be zero.');
          throwAsCaller(ME);             
        end % if statement  

    end % switch statement           
    
  else
    % (n, m) specification
    if sum(abs(ZernikeList(1, :)) ~= ZernikeList(1, :)) ~= 0
      ME = MException('VerifyData:InvalidData', ...
                      'ZernikeCalc: ZernikeList n values must be positive or 0.');
      throwAsCaller(ME);            
    end % if statement
    
    if sum(abs(ZernikeList(2, :)) > ZernikeList(1, :)) ~= 0
      ME = MException('VerifyData:InvalidData', ...
                      'ZernikeCalc: abs(m) values must not exceed n values.');
      throwAsCaller(ME);     
    end % if statement   
    
    switch ZernikeDef
        case 'FRINGE'
              if rows > 1
                ME = MException('VerifyData:InvalidData', ...
                          'ZernikeCalc: FRINGE can only be specified with j ordering.');
                throwAsCaller(ME);              
              end % if statement
              
        case 'ISO'
              if rows > 1
                ME = MException('VerifyData:InvalidData', ...
                          'ZernikeCalc: ISO can only be specified with j ordering.');
                throwAsCaller(ME);              
              end % if statement
              
        case 'WYANT'
              if rows > 1
                ME = MException('VerifyData:InvalidData', ...
                          'ZernikeCalc: WYANT can only be specified with j ordering.');
                throwAsCaller(ME);              
              end % if statement   
              
        case {'HEXAGON', 'HEXAGONR30', 'ELLIPSE', 'RECTANGLE', ...
              'SQUARE', 'ANNULUS'}
              if rows > 1
                ME = MException('VerifyData:InvalidData', ...
                          'ZernikeCalc: ZernikeDef value can only be specified with j ordering.');
                throwAsCaller(ME);              
              end % if statement   
              
        otherwise
          if (sum(mod(ZernikeList(1, :) - abs(ZernikeList(2, :)), 2)) ~= 0)
            ME = MException('VerifyData:InvalidData', ...
                      'ZernikeCalc: n-abs(m) must be even.');
            throwAsCaller(ME);     
          end % if statement   
    end % switch
    
    switch MMESorC
       case {'-m=sin', '-m=cos'}
       otherwise
           ME = MException('VerifyData:InvalidData', ...
                      'ZernikeCalc: Sign convention value for (n, m) not valid.');
           throwAsCaller(ME);  
    end % switch
    
  end % (n,m) order specified  

  % Validate ZernikeShape parameter
  
  switch ZernikeDef    
      case {'ELLIPSE', 'RECTANGLE','ANNULUS'}
         % Check that the ShapeParameter is valid
         if (ShapeParameter <= 0) || (ShapeParameter >= 1)
            ME = MException('VerifyData:InvalidData', ...
                      'ZernikeCalc: The ShapeParameter is out of range: (0,1).');
            throwAsCaller(ME);
         end % if statement
  end % switch ZernikeDef
   
end % function validateInput

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  This section is for implementing non-circular Zernike polynomials.
%  This immplements the tables in Mahajan's paper "Orthonomal polynomials
%  in wavefront analysis: analytical soultion," J. Opt. Soc. Am. A, 24(9)
%  pp. 2994-3016 2007
%

function result = Z(j, rho, theta)
%
% Caluclates the Mahajan Zernike polynomial value.
%

  [n, m, sc] = convertj(j, 'MAHAJAN');
  
  result = calculateZernike(n, m, sc, rho, theta, 'MAHAJAN');

end % function Z
      

function result = ZHexagon(j, rho, theta)

% Orthonormal hexagonal polynomials
  
  if (j < 1) || (j > 45)
      ME = MException('VerifyData:InvalidData', ...
                      'ZernikeCalc: Hexagon order j out of range.');
      throwAsCaller(ME);        
  end % if statement
      
  switch j
      case 1 
          result = Z(1,rho,theta);
      case 2 
          result = sqrt(6/5)*Z(2,rho,theta);
      case 3 
          result = sqrt(6/5)*Z(3,rho,theta);
      case 4 
          result = sqrt(5/43)*Z(1,rho,theta)+(2*sqrt(15/43))*Z(4,rho,theta);
      case 5 
          result = sqrt(10/7)*Z(5,rho,theta);
      case 6 
          result = sqrt(10/7)*Z(6,rho,theta);
      case 7 
          result = 16*sqrt(14/11055)*Z(3,rho,theta)+10*sqrt(35/2211)*Z(7,rho,theta);
      case 8 
          result = 16*sqrt(14/11055)*Z(2,rho,theta)+10*sqrt(35/2211)*Z(8,rho,theta);
      case 9 
          result = (2*sqrt(5)/3)*Z(9,rho,theta);
      case 10 
          result = 2*sqrt(35/103)*Z(10,rho,theta);
      case 11 
          result = (521/sqrt(1072205))*Z(1,rho,theta) ...
                     +88*sqrt(15/214441)*Z(4,rho,theta) ...
                     +14*sqrt(43/4987)*Z(11,rho,theta);
      case 12 
          result = 225*sqrt(6/492583)*Z(6,rho,theta)+42*sqrt(70/70369)*Z(12,rho,theta);
      case 13 
          result = 225*sqrt(6/492583)*Z(5,rho,theta)+42*sqrt(70/70369)*Z(13,rho,theta);
      case 14 
          result = -2525*sqrt(14/297774543)*Z(6,rho,theta) ...
                   -(1495*sqrt(70/99258181)/3)*Z(12,rho,theta) ...
                     +(sqrt(378910/18337)/3)*Z(14,rho,theta);
      case 15 
          result = 2525*sqrt(14/297774543)*Z(5,rho,theta) ... 
                   +(1495*sqrt(70/99258181)/3)*Z(13,rho,theta) ...
                    +(sqrt(378910/18337)/3)*Z(15,rho,theta);
      case 16 
          result = 30857*sqrt(2/3268147641)*Z(2,rho,theta) ...
                    +(49168/sqrt(3268147641))*Z(8,rho,theta) ...
                     +42*sqrt(1474/1478131)*Z(16,rho,theta);
      case 17 
          result = 30857*sqrt(2/3268147641)*Z(3,rho,theta) ...
                     +(49168/sqrt(3268147641))*Z(7,rho,theta) ...
                      +42*sqrt(1474/1478131)*Z(17,rho,theta);
      case 18 
          result = 386*sqrt(770/295894589)*Z(10,rho,theta) ...
                    +6*sqrt(118965/2872763)*Z(18,rho,theta);
      case 19 
          result = 6*sqrt(10/97)*Z(9,rho,theta) ...
                     +14*sqrt(5/291)*Z(19,rho,theta);
      case 20 
          result = -0.71499593*Z(2,rho,theta) ...
                    -0.72488884*Z(8,rho,theta)-0.46636441*Z(16,rho,theta) ...
                     +1.72029850*Z(20,rho,theta);
      case 21 
          result = 0.71499594*Z(3,rho,theta)+0.72488884*Z(7,rho,theta) ...
                   +0.46636441*Z(17,rho,theta)+1.72029850*Z(21,rho,theta);
      case 22 
          result = 0.58113135*Z(1,rho,theta)+0.89024136*Z(4,rho,theta) ...
                    +0.89044507*Z(11,rho,theta)+1.32320623*Z(22,rho,theta);
      case 23  
         result = 1.15667686*Z(5,rho,theta)+1.10775599*Z(13,rho,theta) ...
                    +0.43375081*Z(15,rho,theta)+1.39889072*Z(23,rho,theta);
      case 24 
          result = 1.15667686*Z(6,rho,theta)+1.10775599*Z(12,rho,theta) ...
                    -0.43375081*Z(14,rho,theta)+1.39889072*Z(24,rho,theta);
      case 25 
          result = 1.31832566*Z(5,rho,theta)+1.14465174*Z(13,rho,theta) ...
                     +1.94724032*Z(15,rho,theta)+0.67629133*Z(23,rho,theta) ...
                     +1.75496998*Z(25,rho,theta);
      case 26 
          result = -1.31832566*Z(6,rho,theta)-1.14465174*Z(12,rho,theta) ...
                    +1.94724032*Z(14,rho,theta)-0.67629133*Z(24,rho,theta) ...
                     +1.75496998*Z(26,rho,theta);
      case 27 
          result = 2*sqrt(77/93)*Z(27,rho,theta);
      case 28 
          result = -1.07362889*Z(1,rho,theta)-1.52546162*Z(4,rho,theta) ...
                     -1.28216588*Z(11,rho,theta)-0.70446308*Z(22,rho,theta) ...
                      +2.09532473*Z(28,rho,theta);
      case 29 
          result = 0.97998834*Z(3,rho,theta)+1.16162002*Z(7,rho,theta) ...
                    +1.04573775*Z(17,rho,theta)+0.40808953*Z(21,rho,theta) ...
                     +1.36410394*Z(29,rho,theta);
      case 30 
          result = 0.97998834*Z(2,rho,theta)+1.16162002*Z(8,rho,theta) ... 
                     +1.04573775*Z(16,rho,theta)-0.40808953*Z(20,rho,theta) ...
                      +1.36410394*Z(30,rho,theta);
      case 31 
          result = 3.63513758*Z(9,rho,theta)+2.92084414*Z(19,rho,theta) ...
                   +2.11189625*Z(31,rho,theta);
      case 32 
          result = 0.69734874*Z(10,rho,theta)+0.67589740*Z(18,rho,theta) ...
                     +1.22484055*Z(32,rho,theta);
      case 33 
          result = 1.56189763*Z(3,rho,theta)+1.69985309*Z(7,rho,theta) ... 
                   +1.29338869*Z(17,rho,theta)+2.57680871*Z(21,rho,theta) ...
                    +0.67653220*Z(29,rho,theta)+1.95719339*Z(33,rho,theta);
      case 34 
          result = -1.56189763*Z(2,rho,theta)-1.69985309*Z(8,rho,theta) ...
                    -1.29338869*Z(16,rho,theta)+2.57680871*Z(20,rho,theta) ...
                     -0.67653220*Z(30,rho,theta)+1.95719339*Z(34,rho,theta);
      case 35 
          result = -1.63832594*Z(3,rho,theta)-1.74759886*Z(7,rho,theta) ...
                    -1.27572528*Z(17,rho,theta)-0.77446421*Z(21,rho,theta) ...
                      -0.60947360*Z(29,rho,theta)-0.36228537*Z(33,rho,theta) ...
                       +2.24453237*Z(35,rho,theta);
      case 36 
          result = -1.63832594*Z(2,rho,theta)-1.74759886*Z(8,rho,theta) ...
                    -1.27572528*Z(16,rho,theta)+0.77446421*Z(20,rho,theta) ...
                      -0.60947360*Z(30,rho,theta)+0.36228537*Z(34,rho,theta) ...
                       +2.24453237*Z(36,rho,theta);
      case 37 
          result = 0.82154671*Z(1,rho,theta)+1.27988084*Z(4,rho,theta) ...
                    +1.32912377*Z(11,rho,theta)+1.11636637*Z(22,rho,theta) ...
                     -0.54097038*Z(28,rho,theta)+1.37406534*Z(37,rho,theta);
      case 38 
          result = 1.54526522*Z(6,rho,theta)+1.57785242*Z(12,rho,theta) ...
                     -0.89280081*Z(14,rho,theta)+1.28876176*Z(24,rho,theta) ...
                       -0.60514082*Z(26,rho,theta)+1.43097780*Z(38,rho,theta);
      case 39 
          result = 1.54526522*Z(5,rho,theta)+1.57785242*Z(13,rho,theta) ...
                    +0.89280081*Z(15,rho,theta)+1.28876176*Z(23,rho,theta) ...
                     +0.60514082*Z(25,rho,theta)+1.43097780*Z(39,rho,theta);
      case 40 
          result = -2.51783502*Z(6,rho,theta)-2.38279377*Z(12,rho,theta) ...
                     +3.42458933*Z(14,rho,theta)-1.69296616*Z(24,rho,theta) ...
                       +2.56612920*Z(26,rho,theta)-0.85703819*Z(38,rho,theta) ...
                        +1.89468756*Z(40,rho,theta);
      case 41 
          result = 2.51783502*Z(5,rho,theta)+2.38279377*Z(13,rho,theta) ...
                   +3.42458933*Z(15,rho,theta)+1.69296616*Z(23,rho,theta) ...
                    +2.56612920*Z(25,rho,theta)+0.85703819*Z(39,rho,theta) ...
                     +1.89468756*Z(41,rho,theta);
      case 42 
          result = -2.72919646*Z(1,rho,theta)-4.02313214*Z(4,rho,theta) ...
                    -3.69899239*Z(11,rho,theta)-2.49229315*Z(22,rho,theta) ...
                      +4.36717121*Z(28,rho,theta)-1.13485132*Z(37,rho,theta) ...
                        +2.52330106*Z(42,rho,theta);
      case 43 
          result = 1362*sqrt(77/20334667)*Z(27,rho,theta)+(260/3) ...
                    *sqrt(341/655957)*Z(43,rho,theta);
      case 44 
          result = -2.76678413*Z(6,rho,theta)-2.50005278*Z(12,rho,theta) ...
                    +1.48041348*Z(14,rho,theta)-1.62947374*Z(24,rho,theta) ...
                     +0.95864121*Z(26,rho,theta)-0.69034812*Z(38,rho,theta) ...
                      +0.40743941*Z(40,rho,theta)+2.56965299*Z(44,rho,theta);
      case 45 
          result = -2.76678413*Z(5,rho,theta)-2.50005278*Z(13,rho,theta) ...
                     -1.48041348*Z(15,rho,theta)-1.62947374*Z(23,rho,theta) ...
                      -0.95864121*Z(25,rho,theta)-0.69034812*Z(39,rho,theta) ...
                        -0.40743941*Z(41,rho,theta)+2.56965299*Z(45,rho,theta);

  end % switch j statement

end % function ZHexagon 


function result = ZHexagonR30(j, rho, theta)

% Orthonormal hexagonal polynomials (hexagon rotated 30 degrees)

  if (j < 1) || (j > 28)
      ME = MException('VerifyData:InvalidData', ...
                      'ZernikeCalc: Hexagon R30 order j out of range.');
      throwAsCaller(ME);        
  end % if statement
      
  switch j
      case 1
        result = Z(1,rho,theta);   
      case 2
        result = sqrt(6/5)*Z(2, rho, theta);   
      case 3
        result = sqrt(6/5)*Z(3, rho, theta);   
      case 4
        result = sqrt(5/43)*Z(1, rho, theta)+2*sqrt(15/43)*Z(4, rho, theta);   
      case 5
        result = sqrt(10/7)*Z(5,rho,theta);   
      case 6
        result = sqrt(10/7)*Z(6,rho,theta);   
      case 7
        result = 16*sqrt(14/11055)*Z(3,rho,theta)+10*sqrt(35/2211)*Z(7,rho,theta);   
      case 8
        result = 16*sqrt(14/11055)*Z(2,rho,theta)+10*sqrt(35/2211)*Z(8,rho,theta);   
      case 9
        result = 2*sqrt(35/103)*Z(9,rho,theta);   
      case 10
        result = (2*sqrt(5)/3)*Z(10,rho,theta);   
      case 11
        result = (521/sqrt(1072205))*Z(1,rho,theta)+88*sqrt(15/214441)*Z(4,rho,theta) ...
                    + 14*sqrt(43/4987)*Z(11,rho,theta);   
      case 12
        result = 225*sqrt(6/492583)*Z(6,rho,theta)+42*sqrt(70/70369)*Z(12,rho,theta);   
      case 13 
        result = 225*sqrt(6/492583)*Z(5,rho,theta)+42*sqrt(70/70369)*Z(13,rho,theta);   
      case 14
        result = 2525*sqrt(14/297774543)*Z(6,rho,theta)+(1495*sqrt(70/99258181)/3)*Z(12,rho,theta) ...
                 + (sqrt(378910/18337)/3)*Z(14,rho,theta);   
      case 15
        result = -2525*sqrt(14/297774543)*Z(5,rho,theta)-(1495*sqrt(70/99258181)/3)*Z(13,rho,theta) ...
                   + (sqrt(378910/18337)/3)*Z(15,rho,theta);   
      case 16
        result = 30857*sqrt(2/3268147641)*Z(2,rho,theta)+49168/sqrt(3268147641)*Z(8,rho,theta) ...
                       + 42*sqrt(1474/1478131)*Z(16,rho,theta);   
      case 17
        result = 30857*sqrt(2/3268147641)*Z(3,rho,theta)+(49168/sqrt(3268147641))*Z(7,rho,theta) ...
                    + 42*sqrt(1474/1478131)*Z(17,rho,theta);   
      case 18
        result = 6*sqrt(10/97)*Z(10,rho,theta)+14*sqrt(5/291)*Z(18,rho,theta);   
      case 19
        result = 386*sqrt(770/295894589)*Z(9,rho,theta)+6*sqrt(118965/2872763)*Z(19,rho,theta);   
      case 20
        result = 0.71499593*Z(2,rho,theta)+0.72488884*Z(8,rho,theta) ...
                 + 0.46636441*Z(16,rho,theta)+1.72029850*Z(20,rho,theta);   
      case 21  
        result = -0.71499593*Z(3,rho,theta)-0.72488884*Z(7,rho,theta) ...
                  - 0.46636441*Z(17,rho,theta)+1.72029850*Z(21,rho,theta);   
      case 22
        result = 0.58113135*Z(1,rho,theta)+0.89024136*Z(4,rho,theta) ...
                 + 0.89044507*Z(11,rho,theta)+1.32320623*Z(22,rho,theta);   
      case 23
        result = 1.15667686*Z(5,rho,theta)+1.10775599*Z(13,rho,theta) ...
                  - 0.43375081*Z(15,rho,theta)+1.39889072*Z(23,rho,theta);   
      case 24
        result = 1.15667686*Z(6,rho,theta)+1.10775599*Z(12,rho,theta) ... 
                   + 0.43375081*Z(14,rho,theta)+1.39889072*Z(24,rho,theta);   
      case 25
        result = -1.31832566*Z(5,rho,theta)-1.14465174*Z(13,rho,theta) ... 
                  + 1.94724032*Z(15,rho,theta)-0.67629133*Z(23,rho,theta)+1.75496998*Z(25,rho,theta);   
      case 26
        result = 1.31832566*Z(6,rho,theta)+1.14465174*Z(12,rho,theta) ... 
                  + 1.94724032*Z(14,rho,theta)+0.67629133*Z(24,rho,theta)+1.75496998*Z(26,rho,theta);   
      case 27
        result = 1.81984283*Z(27,rho,theta);   
      case 28
        result = 1.07362889*Z(1,rho,theta)+1.52546162*Z(4,rho,theta) ... 
                  + 1.28216588*Z(11,rho,theta)+0.70446308*Z(22,rho,theta)+2.09532473*Z(28,rho,theta);   
  end % switch j   
  
end % function ZHexagonR30 


function result = ZEllipse(j, rho, theta, b)

% Orthonormal elliptical polynomials

  if (j < 1) || (j > 15)
      ME = MException('VerifyData:InvalidData', ...
                      'ZernikeCalc: Elliptical order j out of range.');
      throwAsCaller(ME);        
  end % if statement
      
  
  alpha = sqrt(45-60*b^2+94*b^4-60*b^6+45*b^8);

  beta = sqrt(1575 - 4800*b^2 + 12020*b^4 - 17280*b^6 + 21066*b^8 - 17280*b^10 ...
         + 12020*b^12 - 4800*b^14 + 1575*b^16);

  gamma = sqrt(35*b^8 - 60*b^6 + 114*b^4 - 60*b^2 + 35);

  delta = sqrt(5 - 6*b^2 + 5*b^4);

   switch j
      case 1 
          result = Z(1,rho,theta); 
      case 2 
          result = Z(2,rho,theta); 
      case 3    
          result = Z(3,rho,theta)/b; 
      case 4 
          result = (sqrt(3)*(1-b^2)/sqrt(3-2*b^2+3*b^4))*Z(1,rho,theta) ...
                   +(2/sqrt(3-2*b^2+3*b^4))*Z(4,rho,theta); 
      case 5 
          result = Z(5,rho,theta)/b; 
      case 6 
          result = -(sqrt(3)*(3-4*b^2+b^4)/(2*b^2*sqrt(6-4*b^2+6*b^4)))*Z(1,rho,theta) ...
                    -3*((1-b^4)/(2*b^2*sqrt(6-4*b^2+6*b^4)))*Z(4,rho,theta) ...
                    +(sqrt(3-2*b^2+3*b^4)/(2*b^2))*Z(6,rho,theta); 
      case 7 
          result = (6*(1-b^2)/(b*sqrt(5-6*b^2+9*b^4)))*Z(3,rho,theta) ...
                    +(2*sqrt(2)/(b*sqrt(5-6*b^2+9*b^4)))*Z(7,rho,theta); 
      case 8 
          result = (2*(1-b^2)/sqrt(9-6*b^2+5*b^4))*Z(2,rho,theta) ...
                    +(2*sqrt(2)/sqrt(9-6*b^2+5*b^4))*Z(8,rho,theta); 
      case 9 
          result = -((5-8*b^2+3*b^4)/(b^3*sqrt(5-6*b^2+9*b^4)))*Z(3,rho,theta) ...
                   -((5-2*b^2-3*b^4)/(2*sqrt(2)*b^3*sqrt(5-6*b^2+9*b^4)))*Z(7,rho,theta) ...
                   +(sqrt(5-6*b^2+9*b^4)/(2*sqrt(2)*b^3))*Z(9,rho,theta); 
      case 10
          result = -((3-4*b^2+b^4)/(b^2*sqrt(9-6*b^2+5*b^4)))*Z(2,rho,theta) ...
                   -((3+2*b^2-5*b^4)/(2*sqrt(2)*b^2*sqrt(9-6*b^2+5*b^4)))*Z(8,rho,theta) ...
                   +(sqrt(9-6*b^2+5*b^4)/(2*sqrt(2)*b^2))*Z(10,rho,theta);           
      case 11 
          result = sqrt(5)*(7 - 10*b^2 + 3*b^4)*alpha^(-1)*Z(1,rho,theta) ... 
                   + 4*sqrt(15)*(1 - b^2)*alpha^(-1)*Z(4,rho,theta) ...
                   -2*sqrt(30)*(1 - b^2)*alpha^(-1)*Z(6,rho,theta) ...
                   +8*alpha^(-1)*Z(11,rho,theta);
      case 12 
         result = (b^2-1)*(sqrt(10)/4)*alpha^(-1)*gamma^(-1)*b^(-2)* ...
                   (195 - 280*b^2+278*b^4-144*b^6+15*b^8)*Z(1,rho,theta) ...
                   +(b^2-1)*(sqrt(30)/4)*b^(-2)*alpha^(-1)*gamma^(-1)* ...
                    (105-100*b^2+94*b^4-20*b^6-15*b^8)*Z(4,rho,theta) ...
                   -(b^2-1)*(sqrt(15)/2)*b^(-2)*alpha^(-1)*gamma^(-1)* ...
                    (75-80*b^2+94*b^4-40*b^6+15*b^8)*Z(6,rho,theta) ...
                   -10*sqrt(2)*b^(-2)*alpha^(-1)*gamma^(-1)*(3-2*b^2+2*b^6-3*b^8)*Z(11,rho,theta) ...
                   + alpha*b^(-2)*gamma^(-1)*Z(12,rho,theta);
               
      case 13 
          result = sqrt(15)*(1 - b^2)/(b*sqrt(5 - 6*b^2 + 5*b^4))*Z(5,rho,theta) ...
                   +2/(b*sqrt(5 - 6*b^2 + 5*b^4))*Z(13,rho,theta); 
      case 14 
            result =  (b^2 - 1)^2*(sqrt(10)/8)*(35-10*b^2-b^4)*b^(-4)*gamma^(-1)*Z(1,rho,theta) ...
                     +(b^2 - 1)^2*5*(sqrt(30)/16)*(7+2*b^2-b^4)*b^(-4)*gamma^(-1)*Z(4,rho,theta) ...
                     -(sqrt(15)/8)*(35-70*b^2+56*b^4-26*b^6+5*b^8)*b^(-4)*gamma^(-1)*Z(6,rho,theta) ...
                     +(b^2 - 1)^2*5*(sqrt(2)/16)*(7+10*b^2+7*b^4)*b^(-4)*gamma^(-1)*Z(11,rho,theta) ...
                     -(5/8)*(7-6*b^2+6*b^6-7*b^8)*b^(-4)*gamma^(-1)*Z(12,rho,theta) ...
                     +(1/8)*gamma*b^(-4)*Z(14,rho,theta);                
                 
      case 15 
            result = -(sqrt(15)/4)*(5-8*b^2+3*b^4)*b^(-3)*delta^(-1)*Z(5,rho,theta) ...
                      +(b^4-1)*(5/4)*b^(-3)*delta^(-1)*Z(13,rho,theta) ...
                      +(1/2)*b^(-3)*delta*Z(15,rho,theta);       
                   
  end % switch statement
  
end % function ZElliptical



function result = ZRectangle(j, rho, theta, a)

% Orthonormal rectangle polynomials

  if (j < 1) || (j > 15)
      ME = MException('VerifyData:InvalidData', ...
                      'ZernikeCalc: Rectangle order j out of range.');
      throwAsCaller(ME);        
  end % if statement
    

  mu = sqrt(9-36*a^2+103*a^4-134*a^6+67*a^8);

  nu = sqrt(49-196*a^2+330*a^4-268*a^6+134*a^8);

  tau = 1/(128*nu*a^4*(1-a^2)^2);

  eta = 9 - 45*a^2 + 139*a^4 - 237*a^6 + 210*a^8 - 67*a^10;

   switch j

      case 1 
           result = Z(1,rho,theta);
      case 2 
           result = (sqrt(3)/(2*a))*Z(2,rho,theta);
      case 3 
           result = (sqrt(3)/(2*sqrt(1-a^2)))*Z(3,rho,theta);
      case 4 
           result = (sqrt(5)/(4*sqrt(1-2*a^2+2*a^4)))*(Z(1,rho,theta)+sqrt(3)*Z(4,rho,theta));
      case 5 
           result = (sqrt(3/2)/(2*a*sqrt(1-a^2)))*Z(5,rho,theta);
      case 6 
           result = (sqrt(5)/(8*a^2*(1-a^2)*sqrt(1-2*a^2+2*a^4)))*((3-10*a^2+12*a^4-8*a^6)*Z(1,rho,theta) ...
                     +sqrt(3)*(1-2*a^2)*Z(4,rho,theta) + sqrt(6)*(1-2*a^2+2*a^4)*Z(6,rho,theta));
      case 7 
           result = (sqrt(21)/(4*sqrt(2)*sqrt(27-81*a^2+116*a^4-62*a^6))) ...
                      * (sqrt(2)*(1+4*a^2)*Z(3,rho,theta)+5*Z(7,rho,theta));
      case 8 
           result = (sqrt(21)/(4*sqrt(2)*a*sqrt(35-70*a^2+62*a^4))) ...
                     * (sqrt(2)*(5-4*a^2)*Z(2,rho,theta)+5*Z(8,rho,theta));
      case 9 
           result = (sqrt(5/2)*sqrt((27-54*a^2+62*a^4)/(1-a^2)) ...
                    / (16*a^2*(27-81*a^2+116*a^4-62*a^6))) ...
                    * (2*sqrt(2)*(9-36*a^2+52*a^4-60*a^6)*Z(3,rho,theta) ...
                    + (9-18*a^2-26*a^4)*Z(7,rho,theta) ...
                    + (27-54*a^2+62*a^4)*Z(9,rho,theta));
      case 10 
           result = (sqrt(5/2)/(16*a^3*(1-a^2)*sqrt(35-70*a^2+62*a^4))) ...
                    * (2*sqrt(2)*(35-112*a^2+128*a^4-60*a^6)*Z(2,rho,theta) ...
                    + (35-70*a^2+26*a^4)*Z(8,rho,theta)+(35-70*a^2+62*a^4)*Z(10,rho,theta));
      case 11 
           result = (1/(16*mu))*(8*(3+4*a^2-4*a^4)*Z(1,rho,theta)+25*sqrt(3)*Z(4,rho,theta) ...
                      + 10*sqrt(6)*(1-2*a^2)*Z(6,rho,theta) ...
	              + 21*sqrt(5)*Z(11,rho,theta));
      case 12     
              alpha = (134*a^8 - 268*a^6 + 330*a^4 - 196*a^2 + 49)^(1/2);
              beta = (67*a^8 - 134*a^6 + 103*a^4 - 36*a^2 + 9)^(1/2);

              result = ((3234*a^10 - 8085*a^8 + 8508*a^6 - 4677*a^4 + 1650*a^2 - 315)/((16*a^4 - 16*a^2)*beta*alpha))* ...                    
                        Z(1,rho,theta) + ...   
                        -(15*3^(1/2)*(14-74*a^2+205*a^4-360*a^6+335*a^8-134*a^10))/(16*a^2*(a^2 - 1)*beta*alpha)* ...  
                        Z(4,rho,theta) + ...
                        -(6^(1/2)*(2340-3975*a^2+3975*a^4+(525)/(a^2*(a^2-1)))/(64*alpha*beta)) * ...
                        Z(6,rho,theta) + ...
                        -(63*sqrt(5)*(1-4*a^2+6*a^4-4*a^6)/(16*a^2*(a^2-1)*alpha*beta))* ...
                        Z(11,rho,theta) + ...
                        -(21*sqrt(10)*beta/(64*a^2*(a^2-1)*alpha))* ...
                        Z(12,rho,theta);
                            
      case 13 
           result = (sqrt(21)/(16*sqrt(2)*a*sqrt(1-3*a^2+4*a^4-2*a^6))) ...
                     * (sqrt(3)*Z(5,rho,theta) + sqrt(5)*Z(13,rho,theta));
      case 14 
           result = tau*(6*(245-1400*a^2+3378*a^4-4452*a^6+3466*a^8-1488*a^10+496*a^12)*Z(1,rho,theta) ...
	                + 15*sqrt(3)*(49-252*a^2+522*a^4-540*a^6+270*a^8)*Z(4,rho,theta) ...
                        + 15*sqrt(6)*(49-252*a^2+534*a^4-596*a^6+360*a^8-144*a^10)*Z(6,rho,theta) ...
                        + 3*sqrt(5)*(49-196*a^2+282*a^4-172*a^6+86*a^8)*Z(11,rho,theta) ...
                        + 147*sqrt(10)*(1-4*a^2+6*a^4-4*a^6)*Z(12,rho,theta) ...
                        + 3*sqrt(10)*nu^2*Z(14,rho,theta));
      case 15 
           result = (1/(32*a^3*(1-a^2)*sqrt(1-3*a^2+4*a^4-2*a^6))) ...
                       * (3*sqrt(7/2)*(5-18*a^2+24*a^4-16*a^6)*Z(5,rho,theta) ...
                       + sqrt(105/2)*(1-2*a^2)*Z(13,rho,theta) ...
                       + sqrt(210)*(1-2*a^2+2*a^4)*Z(15,rho,theta));       
       
   end % switch statement
   
end % function ZRectangle


function result = ZSquare(j, rho, theta)

% Orthonormal square polynomials

  if (j < 1) || (j > 45)
      ME = MException('VerifyData:InvalidData', ...
                      'ZernikeCalc: Square order j out of range.');
      throwAsCaller(ME);        
  end % if statement
    
   switch j
      case 1
          result = Z(1,rho,theta);
      case 2
          result = sqrt(3/2)*Z(2,rho,theta);
      case 3
          result = sqrt(3/2)*Z(3,rho,theta);
      case 4
          result = (sqrt(5/2)/2)*Z(1,rho,theta)+(sqrt(15/2)/2)*Z(4,rho,theta);
      case 5
          result = sqrt(3/2)*Z(5,rho,theta);
      case 6
          result = (sqrt(15)/2)*Z(6,rho,theta);
      case 7
          result = (3*sqrt(21/31)/2)*Z(3,rho,theta)+(5*sqrt(21/62)/2)*Z(7,rho,theta);
      case 8
          result = (3*sqrt(21/31)/2)*Z(2,rho,theta)+(5*sqrt(21/62)/2)*Z(8,rho,theta);
      case 9
          result = -(7*sqrt(5/31)/2)*Z(3,rho,theta)-(13*sqrt(5/62)/4)*Z(7,rho,theta) ...
                      +(sqrt(155/2)/4)*Z(9,rho,theta);
      case 10
          result = (7*sqrt(5/31)/2)*Z(2,rho,theta)+(13*sqrt(5/62)/4)*Z(8,rho,theta) ...
                     +(sqrt(155/2)/4)*Z(10,rho,theta);
      case 11
          result = (8/sqrt(67))*Z(1,rho,theta)+(25*sqrt(3/67)/4)*Z(4,rho,theta) ...
                      +(21*sqrt(5/67)/4)*Z(11,rho,theta);
      case 12
          result = (45*sqrt(3)/16)*Z(6,rho,theta)+(21*sqrt(5)/16)*Z(12,rho,theta);
      case 13
          result = (3*sqrt(7)/8)*Z(5,rho,theta)+(sqrt(105)/8)*Z(13,rho,theta);
      case 14
          result = 261/(8*sqrt(134))*Z(1,rho,theta)+(345*sqrt(3/134)/16)*Z(4,rho,theta) ...
                   +(129*sqrt(5/134)/16)*Z(11,rho,theta)+(3*sqrt(335)/16)*Z(14,rho,theta);
      case 15
          result = (sqrt(105)/4)*Z(15,rho,theta);
      case 16
          result = ((41*sqrt(55/1966))/4)*Z(2,rho,theta)+((29*sqrt(55/983))/4)*Z(8,rho,theta) ...
                     +((11*sqrt(55/983))/4)*Z(10,rho,theta)+((21*sqrt(165/1966))/4)*Z(16,rho,theta);
      case 17
          result = ((41*sqrt(55/1966))/4)*Z(3,rho,theta)+((29*sqrt(55/983))/4)*Z(7,rho,theta) ...
                     -((11*sqrt(55/983))/4)*Z(9,rho,theta)+((21*sqrt(165/1966))/4)*Z(17,rho,theta); 
      case 18
          result = ((34843*sqrt(3/844397))/16)*Z(2,rho,theta)+((20761*sqrt(3/1688794))/8)*Z(8,rho,theta) ...
                     +((32077*sqrt(3/1688794))/8)*Z(10,rho,theta) ...
                     +(22323/(16*sqrt(844397)))*Z(16,rho,theta)+((21*sqrt(983/859))/8)*Z(18,rho,theta);
      case 19
          result = -((34843*sqrt(3/844397))/16)*Z(3,rho,theta)-((20761*sqrt(3/1688794))/8)*Z(7,rho,theta) ...
                     +((32077*sqrt(3/1688794))/8)*Z(9,rho,theta) ...
                      -(22323/(16*sqrt(844397)))*Z(17,rho,theta)+((21*sqrt(983/859))/8)*Z(19,rho,theta);
      case 20
          result = ((1975*sqrt(7/859))/32)*Z(2,rho,theta)+((557*sqrt(7/1718))/8)*Z(8,rho,theta) ...
                     +((377*sqrt(7/1718))/8)*Z(10,rho,theta)+((349*sqrt(21/859))/32)*Z(16,rho,theta) ...
                      +((239*sqrt(21/859))/32)*Z(18,rho,theta)+(sqrt(18039)/32)*Z(20,rho,theta); 
      case 21
          result = ((1975*sqrt(7/859))/32)*Z(3,rho,theta)+((557*sqrt(7/1718))/8)*Z(7,rho,theta) ...
                     -((377*sqrt(7/1718))/8)*Z(9,rho,theta)+((349*sqrt(21/859))/32)*Z(17,rho,theta) ... 
                    -((239*sqrt(21/859))/32)*Z(19,rho,theta)+(sqrt(18039)/32)*Z(21,rho,theta); 
      case 22
          result = ((77*sqrt(65/849))/16)*Z(1,rho,theta)+((65*sqrt(65/283))/16)*Z(4,rho,theta) ...
                    +((75*sqrt(39/283))/16)*Z(11,rho,theta)+((5*sqrt(39/566))/2)*Z(14,rho,theta) ...
                    +((11*sqrt(1365/283))/16)*Z(22,rho,theta); 
      case 23
          result = ((51*sqrt(11/7846))/2)*Z(5,rho,theta)+(7*sqrt(165/7846))*Z(13,rho,theta) ...
                    +((15*sqrt(231/7846))/2)*Z(23,rho,theta); 
      case 24
          result = ((147*sqrt(195/2698))/4)*Z(6,rho,theta)+(105*sqrt(13/2698))*Z(12,rho,theta) ...
                    +((33*sqrt(455/2698))/4)*Z(24,rho,theta); 
      case 25
          result = ((7*sqrt(165))/16)*Z(15,rho,theta)+((3*sqrt(231))/16)*Z(25,rho,theta); 
      case 26
          result = (20525/(64*sqrt(849)))*Z(1,rho,theta)+(15077/(64*sqrt(283)))*Z(4,rho,theta) ...
                     +((2565*sqrt(15/283))/64)*Z(11,rho,theta) ...
                    +((2665*sqrt(15/566))/32)*Z(14,rho,theta)+((749*sqrt(21/283))/64)*Z(22,rho,theta) ...
                     +((3*sqrt(5943/2))/32)*Z(26,rho,theta);
      case 27
          result = ((4911*sqrt(3/3923))/32)*Z(5,rho,theta)+((2429*sqrt(5/3923))/32)*Z(13,rho,theta) ...
                    +((641*sqrt(7/3923))/32)*Z(23,rho,theta)+(sqrt(27461)/32)*Z(27,rho,theta); 
      case 28
          result = ((16877*sqrt(3/2698))/32)*Z(6,rho,theta)+((8295*sqrt(5/2698))/32)*Z(12,rho,theta) ...
                    +((2247*sqrt(7/2698))/32)*Z(24,rho,theta) ...
                   +((3*sqrt(9443/2))/32)*Z(28,rho,theta);
      case 29
          result = 2.42764289*Z(3,rho,theta)+2.69721906*Z(7,rho,theta) ...
                     -1.56598065*Z(9,rho,theta)+2.12208902*Z(17,rho,theta) ...
                     -0.93135654*Z(19,rho,theta)+0.25252773*Z(21,rho,theta) ...
                     +1.59017528*Z(29,rho,theta);
      case 30
          result = 2.42764289*Z(2,rho,theta)+2.69721906*Z(8,rho,theta)+1.56598064*Z(10,rho,theta) ...
                      +2.12208902*Z(16,rho,theta)+0.93135653*Z(18,rho,theta) ...
                    +0.25252773*Z(20,rho,theta)+1.59017528*Z(30,rho,theta);
      case 31
          result = -9.10300982*Z(3,rho,theta)-8.79978208*Z(7,rho,theta)+10.69381427*Z(9,rho,theta) ...
                    -5.37383386*Z(17,rho,theta)+7.01044701*Z(19,rho,theta) ...
                   -1.26347273*Z(21,rho,theta)-1.90131757*Z(29,rho,theta)+3.07960207*Z(31,rho,theta);
      case 32
          result = 9.10300982*Z(2,rho,theta)+8.79978208*Z(8,rho,theta)+10.69381427*Z(10,rho,theta) ...
                    +5.37383385*Z(16,rho,theta)+7.01044701*Z(18,rho,theta) ...
                    +1.26347272*Z(20,rho,theta)+1.90131756*Z(30,rho,theta)+3.07960207*Z(32,rho,theta);
      case 33
          result = 21.39630883*Z(3,rho,theta)+19.76696884*Z(7,rho,theta)-12.70550260*Z(9,rho,theta) ...
                    +11.05819453*Z(17,rho,theta)-7.02178757*Z(19,rho,theta) ...
                    +15.80286172*Z(21,rho,theta)+3.29259996*Z(29,rho,theta) ...
                    -2.07602716*Z(31,rho,theta)+5.40902889*Z(33,rho,theta);
      case 34
          result = 21.39630883*Z(2,rho,theta)+19.76696884*Z(8,rho,theta)+12.70550260*Z(10,rho,theta) ...
                    +11.05819453*Z(16,rho,theta)+7.02178756*Z(18,rho,theta) ...
                    +15.80286172*Z(20,rho,theta)+3.29259996*Z(30,rho,theta) ...
                    +2.07602718*Z(32,rho,theta)+5.40902889*Z(34,rho,theta);
      case 35
          result = -16.54454463*Z(3,rho,theta)-14.89205549*Z(7,rho,theta)+22.18054997*Z(9,rho,theta) ... 
                    -7.94524850*Z(17,rho,theta)+11.85458952*Z(19,rho,theta) ...
                    -6.18963458*Z(21,rho,theta)-2.19431442*Z(29,rho,theta)+3.24324400*Z(31,rho,theta) ...
                     -1.72001172*Z(33,rho,theta)+8.16384008*Z(35,rho,theta);
      case 36
          result = 16.54454462*Z(2,rho,theta)+14.89205549*Z(8,rho,theta)+22.18054997*Z(10,rho,theta) ...
                    +7.94524849*Z(16,rho,theta)+11.85458952*Z(18,rho,theta) ...
                    +6.18963457*Z(20,rho,theta)+2.19431441*Z(30,rho,theta)+3.24324400*Z(32,rho,theta) ...
                    +1.72001172*Z(34,rho,theta)+8.16384008*Z(36,rho,theta);
      case 37
          result = 1.75238960*Z(1,rho,theta)+2.72870567*Z(4,rho,theta)+2.76530671*Z(11,rho,theta) ...
                    +1.43647360*Z(14,rho,theta)+2.12459170*Z(22,rho,theta) ...
                   +0.92450043*Z(26,rho,theta)+1.58545010*Z(37,rho,theta);
      case 38
          result = 19.24848143*Z(6,rho,theta)+16.41468913*Z(12,rho,theta)+9.76776798*Z(24,rho,theta) ...
                    +1.47438007*Z(28,rho,theta)+3.83118509*Z(38,rho,theta);
      case 39
          result = 0.46604820*Z(5,rho,theta)+0.84124290*Z(13,rho,theta)+1.00986774*Z(23,rho,theta) ...
                    -0.42520747*Z(27,rho,theta)+1.30579570*Z(39,rho,theta);
      case 40
          result = 28.18104531*Z(1,rho,theta)+38.52219208*Z(4,rho,theta)+30.18363661*Z(11,rho,theta) ...
                    +36.44278147*Z(14,rho,theta)+15.52577202*Z(22,rho,theta) ...
                    +19.21524879*Z(26,rho,theta)+4.44731721*Z(37,rho,theta)+6.00189814*Z(40,rho,theta);
      case 41
          result = 9.12899823*Z(15,rho,theta)+6.15821579*Z(25,rho,theta)+2.96653218*Z(41,rho,theta);
      case 42
          result = 85.33469748*Z(6,rho,theta)+64.01249391*Z(12,rho,theta)+30.59874671*Z(24,rho,theta) ... 
                    +34.09158819*Z(28,rho,theta)+7.75796322*Z(38,rho,theta) ...
                     +9.37150432*Z(42,rho,theta);
      case 43
          result = 14.30642479*Z(5,rho,theta)+11.17404702*Z(13,rho,theta)...
                    +5.68231935*Z(23,rho,theta)+18.15306055*Z(27,rho,theta)...
                    +1.54919583*Z(39,rho,theta)+5.90178984*Z(43,rho,theta);
      case 44
          result = 36.12567424*Z(1,rho,theta)+47.95305224*Z(4,rho,theta)+35.30691679*Z(11,rho,theta) ...
                    +56.72014548*Z(14,rho,theta)+16.36470429*Z(22,rho,theta) ...
                    +26.32636277*Z(26,rho,theta)+3.95466397*Z(37,rho,theta)+6.33853092*Z(40,rho,theta) ...
                    +12.38056785*Z(44,rho,theta);
      case 45
          result = 21.45429746*Z(15,rho,theta)+9.94633083*Z(25,rho,theta) ...
                    +2.34632890*Z(41,rho,theta)+10.39130049*Z(45,rho,theta);

   end % switch statement
   
end % function ZSquare


function result = ZAnnulus(j, n, m, rho, theta, e)

% Orthonormal annulus polynomials

  if (j < 1) || (j > 37) 
      ME = MException('VerifyData:InvalidData', ...
                      'ZernikeCalc: Annulus order j out of range.');
      throwAsCaller(ME);        
  end % if statement
  
   switch j
      case 1
          % n=0, m=0
          result = ones(size(rho));
       case {2, 3}
          % n=1, m=1 
          result = rho./sqrt(1+e^2); 
       case 4
          % n=2, m=0 
          result = (2.*rho.^2-1-e^2)./(1-e^2);
       case {5, 6}
          % n=2, m=2 
          result = rho.^2./sqrt(1+e^2+e^4);
       case {7, 8}
          % n=3, m=1 
          result = (3*(1+e^2).*rho.^3-2*(1+e^2+e^4).*rho)./...
                    ((1-e^2)*sqrt((1+e^2)*(1+4*e^2+e^4)));
       case {9, 10}
          % n=3, m=3 
          result = rho.^3./sqrt(1+e^2+e^4+e^6);
       case 11
          % n=4, m=0 
          result = (6.*rho.^4-6*(1+e^2).*rho.^2+1+4*e^2+e^4)/(1-e^2)^2;
       case {12, 13}
          % n=4, m=2 
          result = (4.*rho.^4-3*((1-e^8)/(1-e^6)).*rho.^2)./...
                    sqrt((1-e^2)^(-1)*(16*(1-e^10)-15*(1-e^8)^2/(1-e^6)));
       case {14, 15}
          % n=4, m=4 
          result = rho.^4./sqrt(1+e^2+e^4+e^6+e^8);
       case {16, 17}
          % n=5, m=1 
          result = (10*(1+4*e^2+e^4).*rho.^5-12*(1+4*e^2+4*e^4+e^6).*rho.^3 ...
                       + 3*(1+4*e^2+10*e^4+4*e^6+e^8).*rho)./...
                      ((1-e^2)^2*sqrt((1+4*e^2+e^4)*(1+9*e^2+9*e^4+e^6)));
       case {18, 19}
          % n=5, m=3 
          result = (5.*rho.^5-4*((1-e^10)/(1-e^8).*rho.^3))./...
                    sqrt(1/(1-e^2)*(25*(1-e^12)-24*(1-e^10)^2/(1-e^8)));
       case {20, 21}
          % n=5, m=5 
          result = rho.^5./sqrt(1+e^2+e^4+e^6+e^8+e^10);
       case 22
          % n=6, m=0 
          result = (20.*rho.^6-30*(1+e^2).*rho.^4+12*(1+3*e^2+e^4).*rho.^2- ...
                       (1+9*e^2+9*e^4+e^6))/(1-e^2)^3; 
       case {23, 24}
          % n=6, m=2 
          result = (15*(1+4*e^2+10*e^4+4*e^6+e^8).*rho.^6-...
                     20*(1+4*e^2+10*e^4+10*e^6+4*e^8+e^10).*rho.^4+...
                      6*(1+4*e^2+10*e^4+20*e^6+10*e^8+4*e^10+e^12).*rho.^2)./...
                      ((1-e^2)^2*sqrt((1+4*e^2+10*e^4+4*e^6+e^8)*...
                         (1+9*e^2+45*e^4+65*e^6+45*e^8+9*e^10+e^12)));
       case {25, 26}
          % n=6, m=4 
          result = (6.*rho.^6-5*((1-e^12)/(1-e^10)).*rho.^4)./...
                     sqrt((1/(1-e^2))*(36*(1-e^14)-35*(1-e^12)^2/(1-e^10))); 
       case {27, 28}
          % n=6, m=6 
          result = rho.^6/sqrt(1+e^2+e^4+e^6+e^8+e^10+e^12);
       case {29, 30}
          % n=7, m=1
            A710 = (1-e^2)^3 * sqrt(1 + 9*e^2 + 9*e^4 + e^6) * ...
                   sqrt(1 + 16*e^2 + 36*e^4 + 16*e^6 + e^8);
            a71 = 35*(1 + 9*e^2 + 9*e^4 + e^6)/A710;
            b71 = -60 * (1 + 9*e^2 + 15*e^4 + 9*e^6 + e^8)/A710;
            c71 = 30*(1 + 9*e^2 + 25*e^4 + 25*e^6 + 9*e^8 + e^10)/A710;
            d71 = -4*(1 + 9*e^2 + 45*e^4 + 65*e^6 + 45*e^8 + 9*e^10 + e^12)/A710;
            
            result = a71.*rho.^7+b71.*rho.^5+c71.*rho.^3+d71.*rho;  
       case {31, 32}
           % n=7, m=3 
            A730 = (1-e^2)^2 * sqrt(1 + 4*e^2 + 10*e^4 + 20*e^6 + 10*e^8 + 4*e^10 + e^12) * ...
                   sqrt(1 + 9*e^2 + 45*e^4 + 165*e^6 + 270*e^8 + 270*e^10 + ...
                         165*e^12 + 45*e^14 + 9*e^16 + e^18);
            a73 = 21*(1 + 4*e^2 + 10*e^4 + 20*e^6 + 10*e^8 + 4*e^10 + e^12)/A730;
            b73 = -30*(1 + 4*e^2 + 10*e^4 + 20*e^6 + 20*e^8 + 10*e^10 + 4*e^12 + e^14)/A730;
            c73 = 10*(1 + 4*e^2 + 10*e^4 + 20*e^6 + 35*e^8 + 20*e^10 + 10*e^12 + 4*e^14 + e^16)/A730;

            result = a73.*rho.^7+b73.*rho.^5+c73.*rho.^3;
       case {33, 34}   
          % n=7, m=5 
          result = (7.*rho.^7-6*((1-e^14)/(1-e^12)).*rho.^5)/...
                     sqrt((1/(1-e^2))*(49*(1-e^16)-48*(1-e^14)^2/(1-e^12)));
       case {35, 36}
          % n=7, m=7 
          result = rho.^7/sqrt(1+e^2+e^4+e^6+e^8+e^10+e^12+e^14);
       case 37
          % n=8, m=0 
          result = (70*rho.^8-140*(1+e^2)*rho.^6+30*(3+8*e^2+3*e^4)*rho.^4 ...
                     -20*(1+6*e^2+6*e^4+e^6)*rho.^2+(1+16*e^2+36*e^4+16*e^6+e^8))./(1-e^2)^4;
           
   end % switch j   
       
   % Calculate normalization factor.
   
   if m == 0
     factor = sqrt(n+1);
   else
     factor = sqrt(2*(n+1));  
   end % if statement
   
   result = result * factor;
   
   % Calculate sin(), cos() factor
   if (m ~= 0) &&  (mod(j,2) == 0)
     % j is even  
     result = result .* cos(m.*theta);
   end % if statement
   if (m ~= 0) && (mod(j,2) ~= 0)       
     result = result .* sin(m.*theta);
   end % if statement
   
end % function ZAnnulus

function mask=makeDefaultMask(maskType, defaultRows, defaultCols, ShapeParameter)
% This function makes a default mask.  Since it is a default mask, the 
% mask matrix is to be square.

  mask = zeros(defaultRows, defaultCols);
 
  % the circle into which the mask shape is to fit
  cr = (defaultRows+1)/2; 
  cc = (defaultCols+1)/2;
  defaultRadiusInPixels = (min(defaultRows, defaultCols) - 1)/2;
          
  switch maskType
      case 'HEXAGON'
         % make a Hexagon mask
         for r=1:defaultRows
           for c=1:defaultCols
              x = (c-cc); 
              y = (r-cr);
              rho = sqrt(x^2+y^2); 
              eTheta = atan2(y,x);
              
              if (eTheta >= 0) && (eTheta <= (60*pi/180))
                 beta = 120*pi/180 - eTheta;
                 R = defaultRadiusInPixels * sind(60) / sin(beta);
                 if rho < R 
                   mask(r,c) = 1;
                 end % if R statement
              end % if statement    
              
              if eTheta >= (120*pi/180)
                 beta = 120*pi/180 - (pi - eTheta);
                 R = defaultRadiusInPixels * sind(60) / sin(beta);
                 if rho < R 
                   mask(r,c) = 1;
                 end % if R statement
              end % if statement    
              
              if eTheta <= -(120*pi/180)
                 beta = 120*pi/180 - (pi - abs(eTheta));
                 R = defaultRadiusInPixels * sind(60) / sin(beta);
                 if rho < R 
                   mask(r,c) = 1;
                 end % if R statement
              end % if statement 
              
              if (eTheta <= 0) && (eTheta >= -(60*pi/180))
                 beta = 120*pi/180 - abs(eTheta);
                 R = defaultRadiusInPixels * sind(60) / sin(beta);
                 if rho < R 
                   mask(r,c) = 1;
                 end % if R statement
              end % if statement                                 
              
              if (eTheta >= (60*pi/180)) && (eTheta <= (120*pi/180))
                if abs(cr-r) <= ((sqrt(3)/2)*defaultRadiusInPixels)
                  mask(r,c) = 1;  
                end % if statement
              end % if statement
              
              if (eTheta <= (-60*pi/180)) && (eTheta >= (-120*pi/180))
                if abs(cr-r) <= ((sqrt(3)/2)*defaultRadiusInPixels)
                  mask(r,c) = 1;  
                end % if statement
              end % if statement
              
              
            end % for c statement
          end % for r statement
  
         
      case 'HEXAGONR30' 
         % make a Hexagon mask rotated 30 degrees
         for r=1:defaultRows
           for c=1:defaultCols
              x = (c-cc); 
              y = (r-cr);
              rho = sqrt(x^2+y^2); 
              eTheta = atan2(y,x);
              
              if (eTheta >= (30*pi/180)) && (eTheta <= (90*pi/180))
                 beta = 120*pi/180 - (eTheta - 30*pi/180);
                 R = defaultRadiusInPixels * sind(60) / sin(beta);
                 if rho < R 
                   mask(r,c) = 1;
                 end % if R statement
              end % if statement    
              
              if (eTheta >= (90*pi/180)) && (eTheta <= (150*pi/180))
                 beta = 120*pi/180 - (pi - eTheta - 30*pi/180);
                 R = defaultRadiusInPixels * sind(60) / sin(beta);
                 if rho < R 
                   mask(r,c) = 1;
                 end % if R statement
              end % if statement    
              
              if (eTheta <= -(90*pi/180)) && (eTheta >= -(150*pi/180))
                 beta = 120*pi/180 - (pi - abs(eTheta) - 30*pi/180);
                 R = defaultRadiusInPixels * sind(60) / sin(beta);
                 if rho < R 
                   mask(r,c) = 1;
                 end % if R statement
              end % if statement 
              
              if (eTheta <= -30*pi/180) && (eTheta >= -(90*pi/180))
                 beta = 120*pi/180 - (abs(eTheta) - 30*pi/180);
                 R = defaultRadiusInPixels * sind(60) / sin(beta);
                 if rho < R 
                   mask(r,c) = 1;
                 end % if R statement
              end % if statement                                 
              
              if (eTheta >= (-30*pi/180)) && (eTheta <= (30*pi/180))
                if abs(cc-c) <= ((sqrt(3)/2)*defaultRadiusInPixels)
                  mask(r,c) = 1;  
                end % if statement
              end % if statement
              
              if (eTheta >= (150*pi/180)) || (eTheta <= (-150*pi/180))
                if abs(cc-c) <= ((sqrt(3)/2)*defaultRadiusInPixels)
                  mask(r,c) = 1;  
                end % if statement
              end % if statement
              
              
            end % for c statement
          end % for r statement          
          
      case 'ELLIPSE'
          % an ellipse mask is needed
          a = defaultRadiusInPixels;
          b = ShapeParameter * defaultRadiusInPixels;
          for r=1:defaultRows
            for c=1:defaultCols
              rho = sqrt((cr-r)^2 + (cc-c)^2);
              eTheta = atan2((cr-r),(cc-c));
              er = a*b/sqrt((b*cos(eTheta))^2+(a*sin(eTheta))^2);
              if rho <= er
                mask(r,c) = 1;
              end % if statement  
            end % for c statement
          end % for r stateemnt
      case 'RECTANGLE'
          % a rectangular mask is needed
          halfEdgec = ShapeParameter * defaultRadiusInPixels;
          halfEdger = sqrt(1 - ShapeParameter^2) * defaultRadiusInPixels;
          for r=1:defaultRows
            for c=1:defaultCols
              if (r > abs(cr-halfEdger)) && (r < (cr+halfEdger)) && ...
                 (c > abs(cr-halfEdgec)) && (c < (cr+halfEdgec))     
                mask(r,c) = 1;
              end % if statement  
            end % for c statement
          end % for r stateemnt     
      case 'SQUARE'
          % a square mask is needed
          halfEdge = (1/sqrt(2)) * defaultRadiusInPixels;
          for r=1:defaultRows
            for c=1:defaultCols
              if (r > abs(cr-halfEdge)) && (r < (cr+halfEdge)) && ...
                 (c > abs(cr-halfEdge)) && (c < (cr+halfEdge))     
                mask(r,c) = 1;
              end % if statement  
            end % for c statement
          end % for r stateemnt
      case 'ANNULUS'
          % an annulus mask is needed
          obscuration = ShapeParameter * defaultRadiusInPixels;
          for r=1:defaultRows
            for c=1:defaultCols
              rho = sqrt((cr-r)^2 + (cc-c)^2);
              if (rho <= defaultRadiusInPixels) && (rho > obscuration)
                mask(r,c) = 1;
              end % if statement  
            end % for c statement
          end % for r stateemnt
      otherwise
          % a circle mask is needed     
          for r=1:defaultRows
            for c=1:defaultCols
              rho = sqrt((cr-r)^2 + (cc-c)^2);
              if rho <= defaultRadiusInPixels
                mask(r,c) = 1;
              end % if statement  
            end % for c statement
          end % for r stateemnt
       
  end % switch statement

end % function makeDefaultMask

function validateZ()

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
  %  Validate Circle
    
%     fprintf(1,'\n\n    CIRCLE  VALIDATION \n');
%     
%     ZernikeType = 'MAHAJAN';
%     
%     for j1=1:37
%         
%         j2 = j1 - 1;  
%         if j1 == 1
%           j2 = 2;
%         end % if statement
%       
%         [n1,m1,sc1] = convertj(j1, ZernikeType);
%         [n2,m2,sc2] = convertj(j2, ZernikeType);
%            
%         fun1 = @(rho,theta) calculateZernike(n1,m1,sc1,rho,theta,ZernikeType).*calculateZernike(n1,m1,sc1,rho,theta,ZernikeType).*rho; 
%         Q1 = quad2d(fun1, 0,1, 0,2*pi);
%         diff1 = pi - Q1;
%       
%         fun2 = @(rho,theta) calculateZernike(n1,m1,sc1,rho,theta,ZernikeType).*calculateZernike(n2,m2,sc2,rho,theta,ZernikeType).*rho;
%         Q2 = quad2d(fun2, 0,1, 0,2*pi);     
%         diff2 = 0.0 - Q2;    
%       
%         fprintf(1,'j1=%d, n1=%d, m1=%d, result = %+20.18f  Difference = %+20.18f \n',j1,n1,m1,Q1,diff1);
%         fprintf(1,'j2=%d, n2=%d, m2=%d, result = %+20.18f  Difference = %+20.18f \n\n',j2,n2,m2,Q2,diff2);
%           
%      
%     end % for statement

  %
  % end Circle validation 
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  


  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
  %  Validate Circle Fringe
    
%     fprintf(1,'\n\n    CIRCLE (FRINGE)  VALIDATION \n');
%     
%     ZernikeType = 'FRINGE';
%     
%     for j1=1:37
%         
%         j2 = j1 - 1;  
%         if j1 == 1
%           j2 = 2;
%         end % if statement
%       
%         [n1,m1,sc1] = convertj(j1, ZernikeType);
%         [n2,m2,sc2] = convertj(j2, ZernikeType);
%            
%         fun1 = @(rho,theta) calculateZernike(n1,m1,sc1,rho,theta,ZernikeType).*calculateZernike(n1,m1,sc1,rho,theta,ZernikeType).*rho; 
%         Q1 = quad2d(fun1, 0,1, 0,2*pi);
%         % The normalization is different.
%         em = 1;
%         if m1 == 0
%           em = 2;
%         end % if statement  
%         diff1 = (em*pi)/(2*n1+2) - Q1;
%       
%         fun2 = @(rho,theta) calculateZernike(n1,m1,sc1,rho,theta,ZernikeType).*calculateZernike(n2,m2,sc2,rho,theta,ZernikeType).*rho;
%         Q2 = quad2d(fun2, 0,1, 0,2*pi);     
%         diff2 = 0.0 - Q2;    
%       
%         fprintf(1,'j1=%d, n1=%d, m1=%d, result = %+20.18f  Difference = %+20.18f \n',j1,n1,m1,Q1,diff1);
%         fprintf(1,'j2=%d, n2=%d, m2=%d, result = %+20.18f  Difference = %+20.18f \n\n',j2,n2,m2,Q2,diff2);
%           
%      
%     end % for statement

  %
  % end Circle Fringe validation 
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  


  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
  %  Validate ZHexagon
    
%     fprintf(1,'\n\n    HEXAGON  VALIDATION \n');
%     for j1=1:45
%      
%         % Need another j value to compare orthogonality condition
%         j2 = j1 - 1;  
%         if j1 == 1
%           j2 = 2;
%         end % if statement
%            
%       
%         ymina = @(x) -(sqrt(3)*x + sqrt(3));
%         ymaxa = @(x)   sqrt(3)*x + sqrt(3);
%         
%         yminc = @(x) -(-sqrt(3)*x + sqrt(3));
%         ymaxc = @(x)   -sqrt(3)*x + sqrt(3);        
%   
%         fun1 = @(x,y) ZHexagon(j1,sqrt(x.*x+y.*y),atan2(y,x)).*ZHexagon(j1,sqrt(x.*x+y.*y),atan2(y,x)); 
%         Q1a = quad2d(fun1, -1,-0.5, ymina,ymaxa);
%         Q1b = quad2d(fun1, -0.5,0.5, -sqrt(3)/2,sqrt(3)/2);       
%         Q1c = quad2d(fun1, 0.5,1, yminc,ymaxc);
%         diff1 = 3*sqrt(3)/2 - (Q1a+Q1b+Q1c);
%       
%         fun2 = @(x,y) ZHexagon(j1,sqrt(x.*x+y.*y),atan2(y,x)).*ZHexagon(j2,sqrt(x.*x+y.*y),atan2(y,x));
%         Q2a = quad2d(fun2, -1,-0.5, ymina,ymaxa);
%         Q2b = quad2d(fun2, -0.5,0.5, -sqrt(3)/2,sqrt(3)/2);       
%         Q2c = quad2d(fun2, 0.5,1, yminc,ymaxc);  
%         diff2 = 0.0 - (Q2a+Q2b+Q2c);    
%       
%         fprintf(1,'j1=%d, result = %+20.18f  Difference = %+20.18f \n',j1,(Q1a+Q1b),diff1);
%         fprintf(1,'j2=%d, result = %+20.18f  Difference = %+20.18f \n\n',j2,(Q2a+Q2b),diff2);
%      
%       
%     end % for statement

  %
  % end ZHexagon validation 
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
  %  Validate ZHexagonR30
    
%     fprintf(1,'\n\n    HEXAGON ROTATED 30 DEGREES  VALIDATION \n');
%     for j1=1:28
%      
%         % Need another j value to compare orthogonality condition
%         j2 = j1 - 1;  
%         if j1 == 1
%           j2 = 2;
%         end % if statement
%            
%       
%         ymina = @(x) -(x/sqrt(3) + 1);
%         ymaxa = @(x)   x/sqrt(3) + 1;
%         
%         yminb = @(x) -(-x/sqrt(3) + 1);
%         ymaxb = @(x)   -x/sqrt(3) + 1;        
%   
%         fun1 = @(x,y) ZHexagonR30(j1,sqrt(x.*x+y.*y),atan2(y,x)).*ZHexagonR30(j1,sqrt(x.*x+y.*y),atan2(y,x)); 
%         Q1a = quad2d(fun1, -sqrt(3)/2,0, ymina,ymaxa);
%         Q1b = quad2d(fun1, 0,sqrt(3)/2, yminb,ymaxb);
%         diff1 = 3*sqrt(3)/2 - (Q1a+Q1b);
%       
%         fun2 = @(x,y) ZHexagonR30(j1,sqrt(x.*x+y.*y),atan2(y,x)).*ZHexagonR30(j2,sqrt(x.*x+y.*y),atan2(y,x));
%         Q2a = quad2d(fun2, -sqrt(3)/2,0, ymina,ymaxa);
%         Q2b = quad2d(fun2, 0,sqrt(3)/2, yminb,ymaxb);    
%         diff2 = 0.0 - (Q2a+Q2b);    
%       
%         fprintf(1,'j1=%d, result = %+20.18f  Difference = %+20.18f \n',j1,(Q1a+Q1b),diff1);
%         fprintf(1,'j2=%d, result = %+20.18f  Difference = %+20.18f \n\n',j2,(Q2a+Q2b),diff2);
%      
%       
%     end % for statement

  %
  % end ZHexagonR30 validation 
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  

  

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
  %  Validate ZEllipse
% 
%     b = 0.8;  % semi-minor axis
%     %b = 1.0;  % semi-minor axis
%      
%     fprintf(1,'\n\n    ELLIPSE  VALIDATION (b=%f)\n',b);
%     for j1=1:15
%      
%         % Need another j value to compare orthogonality condition
%         j2 = j1 - 1;  
%         if j1 == 1
%           j2 = 2;
%         end % if statement
%            
%         
%         ymax = @(x) b*sqrt(1-x.*x);
%         ymin = @(x) -b*sqrt(1-x.*x);
%         
%         fun1 = @(x,y) ZEllipse(j1,sqrt(x.*x+y.*y),atan2(y,x),b).*ZEllipse(j1,sqrt(x.*x+y.*y),atan2(y,x),b); 
%         Q1 = quad2d(fun1, -1,1, ymin,ymax);
%         diff1 = pi*b - Q1;
%       
%         fun2 = @(x,y) ZEllipse(j1,sqrt(x.*x+y.*y),atan2(y,x),b).*ZEllipse(j2,sqrt(x.*x+y.*y),atan2(y,x),b);
%         Q2 = quad2d(fun2, -1,1, ymin,ymax);     
%         diff2 = 0.0 - Q2;    
%       
%         fprintf(1,'j1=%d, result = %+20.18f  Difference = %+20.18f \n',j1,Q1,diff1);
%         fprintf(1,'j2=%d, result = %+20.18f  Difference = %+20.18f \n\n',j2,Q2,diff2);
%      
%       
%     end % for statement

  %
  % end ZEllipse validation 
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
  %  Validate ZRectangle

%     a = 0.37;  % half length of rectangle
%     %a = 1/sqrt(2);  % for a square
%     
%     fprintf(1,'\n\n    RECTANGLE  VALIDATION (a=%f)\n',a);
%     for j1=1:15
%      
%         % Need another j value to compare orthogonality condition
%         j2 = j1 - 1;  
%         if j1 == 1
%           j2 = 2;
%         end % if statement
%            
%         fun1 = @(x,y) ZRectangle(j1,sqrt(x.*x+y.*y),atan2(y,x),a).*ZRectangle(j1,sqrt(x.*x+y.*y),atan2(y,x),a); 
%         Q1 = quad2d(fun1, -a,a, -sqrt(1-a*a),sqrt(1-a*a));
%         diff1 = (2*a * 2*sqrt(1-a*a)) - Q1;
%       
%         fun2 = @(x,y) ZRectangle(j1,sqrt(x.*x+y.*y),atan2(y,x),a).*ZRectangle(j2,sqrt(x.*x+y.*y),atan2(y,x),a);
%         Q2 = quad2d(fun2, -a,a, -sqrt(1-a*a),sqrt(1-a*a));     
%         diff2 = 0.0 - Q2;    
%       
%         fprintf(1,'j1=%d, result = %+20.18f  Difference = %+20.18f \n',j1,Q1,diff1);
%         fprintf(1,'j2=%d, result = %+20.18f  Difference = %+20.18f \n\n',j2,Q2,diff2);
%      
%       
%     end % for statement

  %
  % end ZRectangle validation 
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
  %  Validate ZSquare

%     fprintf(1,'\n\n    SQUARE  VALIDATION \n');
%     for j1=1:45
%      
%         % Need another j value to compare orthogonality condition
%         j2 = j1 - 1;  
%         if j1 == 1
%           j2 = 2;
%         end % if statement
%            
%         fun1 = @(x,y) ZSquare(j1,sqrt(x.*x+y.*y),atan2(y,x)).*ZSquare(j1,sqrt(x.*x+y.*y),atan2(y,x)); 
%         Q1 = quad2d(fun1, -1/sqrt(2),1/sqrt(2), -1/sqrt(2),1/sqrt(2));
%         diff1 = (sqrt(2)*sqrt(2)) - Q1;
%       
%         fun2 = @(x,y) ZSquare(j1,sqrt(x.*x+y.*y),atan2(y,x)).*ZSquare(j2,sqrt(x.*x+y.*y),atan2(y,x));
%         Q2 = quad2d(fun2, -1/sqrt(2),1/sqrt(2), -1/sqrt(2),1/sqrt(2));     
%         diff2 = 0.0 - Q2;    
%       
%         fprintf(1,'j1=%d, result = %+20.18f  Difference = %+20.18f \n',j1,Q1,diff1);
%         fprintf(1,'j2=%d, result = %+20.18f  Difference = %+20.18f \n\n',j2,Q2,diff2);
%      
%       
%     end % for statement

  %
  % end ZSquare validation 
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
  %  Validate ZAnnulus
  
%     fprintf(1,'\n\n    ANNULUS  VALIDATION \n');
%     for j1=1:37 
%  
%         j2 = j1 - 1;  
%         if j1 == 1
%           j2 = 2;
%         end % if statement
%         
%       
%         [n1,m1,sc1] = convertj(j1, 'MAHAJAN');
%         [n2,m2,sc2] = convertj(j2, 'MAHAJAN');
%            
%         ir = 0.3;  % annulus inner radius
%         fun1 = @(rho,theta) ZAnnulus(j1,n1,m1,rho,theta,ir).*ZAnnulus(j1,n1,m1,rho,theta,ir).*rho; 
%         Q1 = quad2d(fun1, ir,1, 0,2*pi);
%         diff1 = (pi - pi*ir^2) - Q1;
%       
%         fun2 = @(rho,theta) ZAnnulus(j1,n1,m1,rho,theta,ir).*ZAnnulus(j2,n2,m2,rho,theta,ir).*rho;
%         Q2 = quad2d(fun2, ir,1, 0,2*pi);     
%         diff2 = 0.0 - Q2;    
%       
%         fprintf(1,'j1=%d, n1=%d, m1=%d, result = %+20.18f  Difference = %+20.18f \n',j1,n1,m1,Q1,diff1);
%         fprintf(1,'j2=%d, n2=%d, m2=%d, result = %+20.18f  Difference = %+20.18f \n\n',j2,n2,m2,Q2,diff2);
%      
%     end % for statement
%     
  %
  % end annulus validation 
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  
end % validateZ
          
          