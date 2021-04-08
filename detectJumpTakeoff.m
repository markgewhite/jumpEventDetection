% ************************************************************************
% Function: detectJumpTakeoff
% Purpose:  Detect the jump takeoff from the signal alone
%
% Parameters:
%       signal:      cell vector array time series (must be triaxial)
%       idxLanding:  vector array of indexes specifying actual landings
%                    to serve as a backstop;
%                    run detectJumpLanding first 
%
%       opt: options, which include:
%          .nSmoothTO:          moving average time size
%          .idxMaxDivergence    maximum tolerable difference between 
%                               idxXMax and idxDZMax in order to prefer
%                               idxDZMax; if not then prefer idxDXMax
%          .idxOffsetTO:        final fixed adjustment to takeoff index 
%
%
% Output:
%       idxTakeoff:  array of detected takeoffs as indices
%
% ************************************************************************


function idxTakeoff = detectJumpTakeoff( signal, idxLanding, opt )

if size( signal{1}, 2 )~=3
    error('Signal is not 3D.');
end

nCases = size( signal, 1 ); % number of time series
idxXMin = zeros( nCases, 1 );
idxXMax = zeros( nCases, 1 );
idxDXMax = zeros( nCases, 1 );
idxDZMax = zeros( nCases, 1 );
idxTakeoff = zeros( nCases, 1 ); % practical

for i = 1:nCases
       
    % work with the 'vertical' axis
    x = sig( 1:idxLanding(i), 1 );
    z = sig( 1:idxLanding(i), 3 );
    
    % smooth the signal using a moving average
    x = movmean( x, opt.nSmoothTO*2+1 );
    z = movmean( z, opt.nSmoothTO*2+1 );   
    
    % find the rate of change
    dx = centraldiff( x );
    dz = centraldiff( z );
    
    % find the X trough
    [ ~, idxXMin(i) ] = min( x );
    
    % now find the X max after this point
    [~, idxXMaxPt ] = max( x(idxXMin(i):end) );
    idxXMax(i) = idxXMaxPt+idxXMin(i)-1;
       
    % work backwards from XMax and find the peak in the derivative of X
    [ ~, idxDXMax(i) ] = max( dx(1:idxXMax(i)) );
    
    % work backwards from XMax and find the peak in the derivative of Z
    [ ~, idxDZMax(i) ] = max( dz(1:idxXMax(i)) );

    % determine which to use
    % first preference is idxDZMax, then idxDXMax and idxXZero
    if idxXMax(i) - idxDZMax(i) < opt.idxMaxDivergence
        idxTakeoff(i) = idxDZMax(i);
    else
        idxTakeoff(i) = idxDXMax(i);
    end
    
    % apply fixed bias offset
    idxTakeoff(i) = idxTakeoff(i) + opt.idxOffsetTO;   

end

end


