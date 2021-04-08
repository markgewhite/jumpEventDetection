# jumpEventDetection
Algorithms for detecting takeoff and landing in a vertical jump from body-worn accelerometer data

The algorithms were developed for accelerometers attached to the lower back, but should work when the accelerometers are attached to the upper back.

The algorithms are intended to process a set of accelerometer signals stored as a cell array.
e.g. signal{i} refers to an where the rows indicate time and columns indicate the three accelerometer axes


**STEP 1**

[ *landingIdx*, *impactIdx* ] = **detectJumpLanding**( *signals*, *params* );

Inputs:
- *signals* : cell array of triaxial time series 
- *params* : structure specifying the parameters

Parameters:
- *params.nSmoothLD* : moving average time size
- *params.freefallThreshold* : acceleration below which body is assumed to be in freefall (i.e. < 1 g)
- *params.maxSpikeWidth* : max gap between acceleration spikes after landing to prevent false detections
- *params.freefallRange* : search range when the body is expect to be in freefall
- *params.idxOffsetLD* : final fixed adjustment to landing index 

Outputs:
- *impactIdx* : defined as the first acceleration spike upon landing (which may not be the highest)
- *landingIdx* : defined as the estimate of the true landing; it will come before the impact



**STEP 2**

*takeoffIdx* = **detectJumpTakeoff**( *signals*, *landingIdx*, *params* );

Inputs:
- *signals* : cell array of triaxial time series
- *landingIdx* : times obtained from the above algorithm (used as a backstop here)
- *params* : structure specifying the parameters

Parameters:
- *params.nSmoothTO* : moving average time size
- *params.idxMaxDivergence* : maximum tolerable difference between idxXMax and idxDZMax in order to prefer idxDZMax; if not then prefer idxDXMax
- *params.maxSpikeWidth* : max gap between acceleration spikes after landing to prevent false detections
- *params.idxOffsetTO* : final fixed adjustment to takeoff index 

Outputs:
- *takeoffIdx* : defined as the estimate of the true takeoff





