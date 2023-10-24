# Microrheology-analysis
The Matlab script to process the track files (xml.) generated by imageJ TrackMate, generating MSD curves, alpha value list, box and histogram plots
Inputs
%
%   FileTag: File name stem for plots and tables to save. MUST BE IN ''. 
%               ex: 'Sample1'
%   path: file path address dor target directory to save to. MUST BE IN ''. 
%               ex: 'Sample1'
%   TracksFile: Filename of the Trackmate .xml output. MUST BE IN ''. 
%               ex: 'Water1_Tracks.xml'
%   varargin: A stand-in variable for the inputparser, which allows inputs
%             to be optional and take default values when not explicitly given.
%   
%   Optional POSITIONAL Inputs: TracksFile(2-6):Filename of the Trackmate 
%             .xml output. MUST BE IN ''. Must be listed BEFORE
%             nonpositional optional inputs and after required inputs.
%               ex: 'Water2_Tracks.xml'
%
%   Optional NONpositional Inputs: (period, fps, sz, plotopt, tableopt)
%
%   period: total length of video in seconds. Defaults to 8s
%   fps: framerate of video in 1/seconds. Defaults to 30fps.
%   sz: size of beads in pixels. Defauls to 13px (per um).
%   plotopt: option to save plots generated. Script saves if plotopt = 1
%            (default)
%   tableopt:option to save tables generated. Script saves if plotopt = 1
%            (default)
%   
%
%Hailey Currie
%Fall 2023
