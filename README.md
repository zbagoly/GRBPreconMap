# GRBPreconMap
This code is licenced under GNU GENERAL PUBLIC LICENSE, Version 3, 29 June 2007.

You are welcome to use these program for your own scientific purposes.  
If you do so we would appreciate it if you include a reference to the paper announcing it: 
Recovering Gamma-Ray Burst Redshift Completeness Maps via Spherical Generalized Additive Models
          Z. Bagoly and  I. I. Racz  submitted to Universe 

We use a comprehensive database of the Swift GRBs with measured positions and spectroscopic redshifts.  Here we use the Swift XRT enhanced positional information, if available for a given event, otherwise we use the Swift BAT position.

The selection results in a total of 1665 GRBs in the blue.csv file.  Of this sample, 427 GRBs possess precise spectroscopic redshifts, with 300 derived from afterglow observations (bluea.csv) and 127 from host galaxy measurements (blueh.csv).

In the hpixlb5.csv there are 12288 HEALPix grid positions.

GRBPreconMapAHSplit.R is the main CRAN R program. 
