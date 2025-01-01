# Radar_Sounder_Model
Model inputs
freq   % center frequency
tp     % 1/4 of pulse length
theta0 % incident angle
gx     % along-track beam pattern
gy     % cross-track beam pattern
sig1   % upper surface rms height
L1     % upper surface correlation length
sig2   % lower surface rms height
L2     % lower surface correlation length
h0     % orbit height
ep2r   % dielectric constant of dry regolith
ep3r   % dielectric constant of the lower layer
d      % layer depth

Radar corss sections are represented as functions of tau (time difference from the 1st pulse peak). 
The x-axis parameter is tau. It should be increased when layer depth is large.
