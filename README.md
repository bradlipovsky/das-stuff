# DAS Stuff - Local Earthquakes
  
We would like work in this branch to ultimately result in an automated workflow for seismic velocity inversion using local earthquake tomography from earthquakes recorded on the Whidbey Island or SeaDAS-N distributed acoustic sensing arrays. The end goal is a publication in a short format journal. 


To date we have automated the following pieces of the workflow:

- Earthquake catalog search and data retrieval 
- Read in the Sintela Onyx waveforms and do simple preprocesses (detrend, taper, filter)
- Pick and correct travel times for P-wave arrivals
- Interpolate latitude and longitude for channel locations along fiber route

We would like to continue work on this project into the future to impliment the velocity tomography once we have reliable travel times for multiple earthquakes. Future work would include development of the code and improvements to what we already have. This would include:

- Developing the auto-picker to generalize across the 2 DAS arrays (Whidbey and SeaDAS) and be robust to variable data quality (i.e. different events have variable data quality).
- Update and constrain the fiber location algorithm with more tap testing 
- Implement the tomography on our processed data


### (c) John-Morgan Manos and Parker Sprinkle


