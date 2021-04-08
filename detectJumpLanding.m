% ************************************************************************
% Function: detectJumpLanding
% Purpose:  Detect the jump landing from the accelerometer signal alone
%
% Parameters:
%       signal:      cell vector array time series
%       opt: options, which include:
%          .nSmooth:            moving average time size
%          .freefallThreshold:  acceleration below which body is assumed
%                               to be in freefall (i.e. < 1 g)
%          .maxSpikeWidth:      max gap between acceleration spikes after 
%                               landing to prevent false detections
%          .freefallRange:      search range when the body is expect to be
%                               in freefall
%          .idxOffset:          final fixed adjustment to landing index 
%
% Output:
%       idxLanding:  detected landing index
%       idxImpact:   detect impact index
%
% ************************************************************************


function [ idxLanding, idxImpact ] = detectJumpLanding( signal, opt )

nCases = size( signal, 1 ); % number of time series
idxLanding = zeros( nCases, 1 ); % practical
idxImpact = zeros( nCases, 1 ); % alternative

for i = 1:nCases

    % work with the resultant (shorthand)
    r = sqrt( sum( signal{i}.^2, 2) );
    
    % smooth the signal using a moving average
    r = movmean( r, opt.nSmooth*2+1 );
    
    % detect the peak landing impact
    [ pkSpike, idxSpike ] = findpeaks( r, ...
                                'MinPeakHeight', opt.freefallThreshold );
    % validate each candidates
    % before each there must be a verifiable period of freefall
    nCandidates = length( idxSpike );
    idxCandidate = zeros( nCandidates, 1 );
    rFreefall = zeros( nCandidates, 1 );
    
    for j = 1:nCandidates
        % find the first point when the r falls below a threshold
        % working backwards from the peak impact
        if any(r( max(idxSpike(j)-opt.maxSpikeWidth,1) : idxSpike(j) ) ...
                                    < opt.freefallThreshold)
            % accept the candidate
            idxCandidate(j) = find( ...
                    r( 1:idxSpike(j)) < opt.freefallThreshold, 1, 'last' );
            % find the mean prior to this point over a set period
            rFreefall(j) = mean( r( max(idxCandidate(j)-opt.freefallRange,1) ...
                                : idxCandidate(j) ) );
        end
    end
    
    % remove candidates that did not meet criteria
    retained = idxCandidate > 0;
    if sum( retained ) > 0
        % valid candidate identified
        idxCandidate = idxCandidate( retained );
        rFreefall = rFreefall( retained );
        pkSpike = pkSpike( retained );
        idxSpike = idxSpike( retained );

        % the best candidate is the one with the highest peak-to-lull ratio
        [ ~, best ] = max( pkSpike./rFreefall );

        % record the impact spike
        idxImpact(i) = idxSpike( best );

        % use the nearest apply fixed bias offset
        idxLanding(i) = idxCandidate( best ) + opt.idxOffset;
    
    else
        % no valid candidate identified
        idxImpact(i) = -1;
        idxLanding(i) = -1;
    
    end

end

end


