%% Analyze TrackMate .xml outputs to plot tracks and calculate alphas. 
%
%Inputs
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
%Xuening edited version
%
function [Alpha, msd] = MicroRheology(FileTag, path, TracksFile, varargin)

%% Input Parser
%Input parser efficiently sets defaults for inputs and checks if given
%inputs are valid/correct data type. Parameters are not positional, so the
%function must be fed the variable name and new value when called if not
%using defaults. ex: Rheology(objs_link, dt = 1/10). Unmentioned variables default.
p = inputParser;

% For doing various checks on the inputs
% valid_positive_integer = @(x) isnumeric(x) && isscalar(x) && (x > 0) && (mod(x,1)<eps); % last check enforces integer
% valid_positive_integer_array = @(x) isnumeric(x)  && sum(x > 0)==length(x) && sum((mod(x,1)<eps))==length(x); % last check enforces integer
valid_positive_number = @(x) isnumeric(x) && isscalar(x) && (x > 0);

%Optional additional track files for batch processing. These are positional
%and must come before other optional variables. NOT passed in the same way
%as "Parameters' below. Just provide filename in '' like the first
addOptional(p,'TracksFile2',[], @mustBeText); 
addOptional(p,'TracksFile3',[], @mustBeText); 
addOptional(p,'TracksFile4',[], @mustBeText); 
addOptional(p,'TracksFile5',[], @mustBeText); 
addOptional(p,'TracksFile6',[], @mustBeText); 
addOptional(p,'TracksFile7',[], @mustBeText); 
addOptional(p,'TracksFile8',[], @mustBeText); 
addOptional(p,'TracksFile9',[], @mustBeText); 
addOptional(p,'TracksFile10',[], @mustBeText); 
addOptional(p,'TracksFile11',[], @mustBeText); 
addOptional(p,'TracksFile12',[], @mustBeText); 

% Default values, and check validity
addParameter(p,'period', 8, valid_positive_number); %Duration of video (s). Default 8s
addParameter(p,'fps', 30, valid_positive_number); %Frame rate of video (1/s). Default 30fps.
addParameter(p,'sz', 13, valid_positive_number); %Spatial scale (px/um). Default 13px/um.
addParameter(p,'plotopt', 1); %Option to save plot. Default 1 (save plot)
addParameter(p,'tableopt', 1); %Option to save data tables. Default 1 (save tables)

parse(p, varargin{:});

%To check if any inputs did nothing/were not interpreted:
if ~isempty(fieldnames(p.Unmatched))
   disp('Extra inputs:')
   disp(p.Unmatched)
end

%To see which variables were set to their defaults, and what those values
%are:
if ~isempty(p.UsingDefaults)
   disp('Using defaults: ')
   disp(p.UsingDefaults)
end

%Call & rename variables from parser data structure
TracksFile2 = p.Results.TracksFile2;
TracksFile3 = p.Results.TracksFile3;
TracksFile4 = p.Results.TracksFile4;
TracksFile5 = p.Results.TracksFile5;
TracksFile6 = p.Results.TracksFile6;
TracksFile7 = p.Results.TracksFile6;
TracksFile8 = p.Results.TracksFile6;
TracksFile9 = p.Results.TracksFile6;
TracksFile10 = p.Results.TracksFile6;
TracksFile11 = p.Results.TracksFile6;
TracksFile12 = p.Results.TracksFile6;
period = p.Results.period;
fps = p.Results.fps;
sz = p.Results.sz;
plotopt = p.Results.plotopt;
tableopt = p.Results.tableopt;

%% Import Tracks Data
%Import tracks. clipZ is set to true to remove the z column preemptively
Tracks = importTrackMateTracks(TracksFile, true);

%% Import optional additional track files for batch processing
%Check if TracksFile2 has been provided (empty [] by default if not)
if 0 == isempty(TracksFile2)
   %Import tracks data if file provided
   Tracks2 = importTrackMateTracks(TracksFile2, true);
   % Concatenate tracks data to first file's tracks data
   Tracks = vertcat(Tracks,Tracks2);
end
% So on for other tracks files...
if 0 == isempty(TracksFile3)
   Tracks3 = importTrackMateTracks(TracksFile3, true);
   Tracks = vertcat(Tracks,Tracks3);
end
if 0 == isempty(TracksFile4)
   Tracks4 = importTrackMateTracks(TracksFile4, true); 
   Tracks = vertcat(Tracks,Tracks4);
end
if 0 == isempty(TracksFile5)
   Tracks5 = importTrackMateTracks(TracksFile5, true);
   Tracks = vertcat(Tracks,Tracks5);
end
if 0 == isempty(TracksFile6)
   Tracks6 = importTrackMateTracks(TracksFile6, true);
   Tracks = vertcat(Tracks,Tracks6);
end
if 0 == isempty(TracksFile7)
   Tracks7 = importTrackMateTracks(TracksFile7, true);
   Tracks = vertcat(Tracks,Tracks7);
end
if 0 == isempty(TracksFile8)
   Tracks8 = importTrackMateTracks(TracksFile8, true);
   Tracks = vertcat(Tracks,Tracks8);
end
if 0 == isempty(TracksFile9)
   Tracks9 = importTrackMateTracks(TracksFile9, true);
   Tracks = vertcat(Tracks,Tracks9);
end
if 0 == isempty(TracksFile10)
   Tracks10 = importTrackMateTracks(TracksFile10, true);
   Tracks = vertcat(Tracks,Tracks10);
end
if 0 == isempty(TracksFile11)
   Tracks11 = importTrackMateTracks(TracksFile11, true);
   Tracks = vertcat(Tracks,Tracks11);
end
if 0 == isempty(TracksFile12)
   Tracks12 = importTrackMateTracks(TracksFile12, true);
   Tracks = vertcat(Tracks,Tracks12);
end

%% Reformat Tracks data
%Convert frames to seconds and pixels to microns
beads = size(Tracks, 1);
frames = fps*period;
for i=1:beads
Tracks{i,1}(:,1)=(period/frames)*Tracks{i,1}(:,1);
Tracks{i,1}(:,2)=(1/sz)*Tracks{i,1}(:,2);
Tracks{i,1}(:,3)=(1/sz)*Tracks{i,1}(:,3);
end

%% Begin Analysis
%Generate analyzer object and feed in tracks
msd = msdanalyzer(2,'microns','seconds');
msd = msd.addAll(Tracks);
%Calculate and account for drift
msd = msd.computeDrift('velocity');
msd = msd.fitLogLogMSD(.1);


%% Screen out bad tracks
%find good r2 indices, take only MSD for tracks corresponding to good r2
r2array = msd.loglogfit.r2fit;
goodindicesr2 = find(r2array>0.7);
alphaarray = msd.loglogfit.alpha;
goodindicesalpha = find(alphaarray>0 & alphaarray<1);
goodindices = intersect(goodindicesr2, goodindicesalpha);
h=histogram(alphaarray)
%% (Inactive) Convert MSD data to array
% Commented out but saved in case it will be useful
% MSD in cell array type. Change it to a matrix [dt mean std N particlenumber] 
% % MSDarray = cell2mat(msd.msd(goodindices));
% % N = size(cell2mat(msd.msd(1,1)),1);
% % MSDarray = transpose(MSDarray);
% % goodbeads = (size(MSDarray,2))/N;
% % v = 1:goodbeads;
% % particlenum = repelem(v,N);
% % MSDTable = vertcat(MSDarray,particlenum);

%% Generate plots
% Clear any pre-existing graphics from current axes
cla reset
%plot mean MSD
% msd.plotMeanMSD(gca, false, goodindices)
% axis([0,1,0,10])
% if plotopt == 1 %optionally save plot
% plotaddress1 = append(FileTag,'MeanMSDPlot.jpeg');
%     saveas(gcf,fullfile(path, plotaddress1));
% end
% Clear any pre-existing graphics from current axes
cla reset
%plot MSD

msd.plotMSD(gca, goodindices, false)
set(gca, 'ColorOrder', colormap(gray),'XScale','log','YScale','log') % color and log-log scale set here
msd.plotMeanMSD(gca, false, goodindices)
axis([0,1,0,10])

if plotopt == 1 %optionally save plot
plotaddress2 = append(FileTag,'MSDPlot.jpeg');
saveas(gcf,fullfile(path, plotaddress2));
end
%plot alpha distribution
cla reset
alpha_values = msd.loglogfit.alpha(goodindices)
Hist, bin_edges = np.histogram(alpha_values, bins=100, density=True)
bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
hist_percentage = Hist / np.sum(Hist) * 100
plt.bar(bin_centers, hist_percentage, width=0.01)
ax.set_xlim(0, 1) 

%histogram(msd.loglogfit.alpha(goodindices),100),axis([0,1,0,20])
xlabel('Alpha')
ylabel('Frequency')
plotaddress3 = append(FileTag,'-Alpha Histogram.jpeg');
saveas(gcf,fullfile(path, plotaddress3));
%plot box plot
cla reset
boxchart(msd.loglogfit.alpha(goodindices)),ylim([0,1])
ylabel('Alpha')
plotaddress4 = append(FileTag,'-Alpha Box-whisker chart.jpeg');
saveas(gcf,fullfile(path, plotaddress4));

%% Table Data
if tableopt == 1
% First table is rows of particle x columns alpha, r^2
headers1 = {'Alpha', 'R2'};
Alpha = table(msd.loglogfit.alpha(goodindices),msd.loglogfit.r2fit(goodindices),'VariableNames', headers1);
filename1 = append(FileTag,'Alpha_R2_Table.xlsx');
fileaddress1 = fullfile(path, filename1);
writetable(Alpha,fileaddress1,'Sheet','AlphaR2')
% This table is Mean MSD, an N x 4 [ dT M STD N ] (see getMeanMSD.m)
headers3 = {'dT', 'Mean', 'std', 'N'};
MeanMSD = array2table(msd.getMeanMSD(goodindices), 'VariableNames', headers3);
filename3 = append(FileTag,'MeanMSD.Table.xlsx');
fileaddress3 = fullfile(path, filename3);
writetable(MeanMSD,fileaddress3,'Sheet','MeanMSD')
% Table code for MSD table, saved for posterity
% This table is  [dt mean std N (particle number)]
% MSDinfo = table(MSDTable);
% filename2 = append(FileTag,'MSDInfo.Table.xlsx');
% fileaddress2 = fullfile(path, filename2);
% writetable(MSDinfo,fileaddress2,'Sheet','MSDInfo')
end
