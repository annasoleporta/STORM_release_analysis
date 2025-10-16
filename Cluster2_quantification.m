function [C, ClustSize, Loc2particleMain, diam]=Cluster2_quantification(FileNameMain, FileNameRef, FileNameFid, Bandwidth, MinPts, MaxDiam, FactorMaxDist, varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% CLUSTER2: perform mean-shift clustering of (X,Y,T) coordinates of
%%%%%% REFERENCE and FIDUCIAL channel in order to roughly identify centers 
%%%%%  of valid particles and fiducial markers.
%%%%%% Then, the localizations in MAIN channel within a defined distance
%%%%%% from centers are stored.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%% REQUIRED INPUTS:
%%%%%% FileNameMain: txt file with XYT coords of main channel, in nm and frames
%%%%%% FileNameRef: txt file with XYT coords of ref channel, in nm and frames
%%%%%% FileNameFid: txt file with XYT coords of fid channel, in nm and frames
%%%%%% Bandwidth: parameter for clustering, in nm
%%%%%% MinPts: minimum number of localiz in a cluster
%%%%%% MaxDiam: maximum diameter (longest axis), in nm
%%%%%% FactorMaxDist: factor*radius sets the maximum distance between particle center and localization in main channel to be considered attached

%%%%%% OPTIONAL INPUTS:
%%%%%% Elong: max ellipse elongation allowed, default=10.0 
%%%%%% ScaleFactor: scale factor in ellipse fit, default=1.0
%%%%%% MinClustDist: minimum distance between clusters (in nm)to be non-aggregated, default=300 
%%%%%% DistFidRef: below this distance, a NC is considered a fiducial marker

%%%%%% OUTPUTS:
%%%%%% Loc2particleMain: cell array with XYT of MAIN channel for each selected cluster
%%%%%% C: XY coordinates of detected particles centers
%%%%%% ClustSizeMain: for each selected REF cluster, number of counted localizations in the MAIN channel
%%%%%% ClustSizeRef: for each selected REF cluster, number of counted localizations in the REF channel
%%%%%% diam: diameter of each selected cluster (in nm)

%%%%%%%%% beginning of function...%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Parse input data =====================================================

tic

p = inputParser; %init parser object
validNum = @(x) isnumeric(x) && (x > 0); %define valid inputs: positive num
validChar = @(x) ischar(x);

%define defaults values for optional param:
defaultElong = 10.0; %elongation in ellipse fit
defaultScaleFactor= 1.0; %scale factor in ellipse fit
defaultMinClustDist = 300.0; % min distance between clusters
defaultDistFidRef = 100.0; % distance between NC and fid

%define required and optional input parameters:
addRequired(p,'FileNameMain',validChar);
addRequired(p,'FileNameRef',validChar);
addRequired(p,'FileNameFid',validChar);
addRequired(p,'Bandwidth',validNum); 
addRequired(p,'MinPts',validNum);
addRequired(p,'MaxDiam', validNum);
addRequired(p,'FactorMaxDist', validNum);

addParameter(p,'Elong', defaultElong, validNum);
addParameter(p,'ScaleFactor', defaultScaleFactor, validNum);
addParameter(p,'MinClustDist', defaultMinClustDist, validNum);
addParameter(p,'DistFidRef', defaultDistFidRef, validNum);

%read input values:
parse(p,FileNameMain, FileNameRef, FileNameFid, Bandwidth, MinPts, MaxDiam, FactorMaxDist, varargin{:});

%assign the parsed values:
Last = Inf; % Limit the number of points processed (set to Inf to process all)
MaxParticleElongation = p.Results.Elong; % max elongation allowed in ellipse fit
EllipseFitScaleFactor = p.Results.ScaleFactor; % scale factor in ellipse fit 
MinClustDist = p.Results.MinClustDist; % min distance of closest cluster to be considered isolated
DistFidRef = p.Results.DistFidRef;


%%% Read coordinates ======================================================

disp('Reading data...');

% read XYT coords on MAIN channel in txt file:
CoordinatesMain=importdata(FileNameMain);
% %%%read txt with more complex file structure:
% % delimiterIn = ' ';
% % headerlinesIn = 1;
% % Coordinates = importdata(FileNameMain,delimiterIn,headerlinesIn);

% assign X-Y-T coordinates of MAIN channel:
XCoords488 = CoordinatesMain(:,1); %X coords in first column
YCoords488 = CoordinatesMain(:,2); %Y coords in second column
TCoords488 = CoordinatesMain(:,3); %T coords in third column

% read XYT coords on REF channel in txt file:
CoordinatesRef=importdata(FileNameRef);
% %%%read txt with more complex file structure:
% % delimiterIn = ' ';
% % headerlinesIn = 1;
% % Coordinates = importdata(FileNameFid, FileNameRef,delimiterIn,headerlinesIn);

% assign X-Y-T coordinates of REF channel:
XCoordsRef = CoordinatesRef(:,1); %X coords in first column
YCoordsRef = CoordinatesRef(:,2); %Y coords in second column
%TCoordsRef = CoordinatesRef(:,3); %T coords in third column

% read XYT coords on FID channel in txt file:
CoordinatesFid=importdata(FileNameFid);
% %%%read txt with more complex file structure:
% % delimiterIn = ' ';
% % headerlinesIn = 1;
% % Coordinates = importdata(FileNameFid,delimiterIn,headerlinesIn);

% assign X-Y-T coordinates of FID channel:
XCoordsFid = CoordinatesFid(:,1); %X coords in first column
YCoordsFid = CoordinatesFid(:,2); %Y coords in second column
%TCoordsFid = CoordinatesFid(:,3); %T coords in third column

% Plot raw coordinates data:
XCoords488 = XCoords488(1:min(Last,numel(XCoords488))); 
YCoords488 = YCoords488(1:min(Last,numel(YCoords488)));
XCoordsRef = XCoordsRef(1:min(Last,numel(XCoordsRef))); 
YCoordsRef = YCoordsRef(1:min(Last,numel(YCoordsRef)));
XCoordsFid = XCoordsFid(1:min(Last,numel(XCoordsFid))); 
YCoordsFid = YCoordsFid(1:min(Last,numel(YCoordsFid)));
figure(1);
plot(XCoords488,YCoords488,'.g'); axis equal; hold on;
plot(XCoordsRef,YCoordsRef,'.r'); axis equal; hold on;
plot(XCoordsFid,YCoordsFid,'.b'); axis equal; hold on;
title({'Magenta: valid, Cyan: elongated, Yellow: aggregated, Blue: fiducials'});
xlabel("x (nm)"); ylabel("y (nm)");


%%% Clustering in FID channel =============================================

% Clustering using mean-shift:
disp('Clustering...');
PtsFid = [XCoordsFid YCoordsFid];
[clustCentFid,~,clustMembsCellFid] = MeanShiftCluster(PtsFid.',Bandwidth);
%ClustSizeFid = cellfun(@numel, clustMembsCellFid);
NClustFid = numel(clustMembsCellFid);
disp(strcat(['Found ' num2str(NClustFid) ' clusters in FID channel']));

%%% Clustering in REF channel =============================================

% Clustering using mean-shift:
disp('Clustering...');
PtsRef = [XCoordsRef YCoordsRef];
[clustCentRef,~,clustMembsCellRef] = MeanShiftCluster(PtsRef.',Bandwidth);
%ClustSizeRef = cellfun(@numel, clustMembsCellRef);
NClustRef = numel(clustMembsCellRef);
disp(strcat(['Found ' num2str(NClustRef) ' clusters in REF channel']));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FILTERING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Filtering REF channel...');
RightClustRef=1:1:NClustRef; % Initialization cluster indexes;
massCenter=clustCentRef; %Initialization centers coords of clusters;

%%% Filter by localization number =========================================

IndexRightClustRef_loc=false(NClustRef,1); % Initialization logical array to identify valid clusters
r_loc=0; % Counter for clusters excluded due to few localizations
for i = 1:NClustRef
    A = [XCoordsRef(cell2mat(clustMembsCellRef(i)))'; YCoordsRef(cell2mat(clustMembsCellRef(i)))'];
    if(size(A,2) >= MinPts)
        IndexRightClustRef_loc(i) = true; % Mark this cluster index as valid (1) if it meets the condition
    else
        r_loc = r_loc + 1; % Increment counter for clusters excluded due to few localizations
    end
end
disp(strcat( [num2str(r_loc) ' REF clusters have been excluded because of few localizations']));

% Update filtered clusters
RightClustRef_loc = RightClustRef(IndexRightClustRef_loc);
massCenter_loc = massCenter(:, IndexRightClustRef_loc);
clustMembsCellRef_loc = clustMembsCellRef(IndexRightClustRef_loc);

%%% Filter aggregated nanoparticles =======================================

% Calculate the distance matrix between the filtered clusters
distCenters = pdist2(massCenter_loc',massCenter_loc','euclidean');

% Find indices of clusters that are too close to other clusters
[rowsAgg, ~] = find((distCenters~=0) & (distCenters<MinClustDist)); %index of clusters closer than MinClustDist and different from 0
Aggregates = rowsAgg';

% Initialization of logic index with all values = true
IndexRightClustRef_agg = true(length(RightClustRef_loc),1);

% Assign false to indexes corresponding to aggregates
IndexRightClustRef_agg(Aggregates) = false;

% Draw yellow ellipse for aggregates
for w = Aggregates
    A = [XCoordsRef(cell2mat(clustMembsCellRef_loc(w)))'; YCoordsRef(cell2mat(clustMembsCellRef_loc(w)))'];
    ellipse_t = fit_ellipse_aggregates(A(1,:), A(2,:), gcf, EllipseFitScaleFactor, MaxParticleElongation, MaxDiam, true); 
end

% Remove aggregates from filtered clusters
RightClustRef_agg = RightClustRef_loc(IndexRightClustRef_agg);
massCenter_agg = massCenter_loc(:, IndexRightClustRef_agg);
clustMembsCellRef_agg = clustMembsCellRef_loc(IndexRightClustRef_agg);
disp(strcat( [num2str(length(Aggregates)) ' REF clusters have been excluded because they are aggregates']));

%%% Filter by elongation ==================================================

IndexRightClustRef_elong=false(length(RightClustRef_agg),1); % Initialization logical array to identify valid clusters
isolatedRef = IndexRightClustRef_agg(IndexRightClustRef_agg == 1);
r_elong = 0; % Counter for clusters excluded due to elongation

for j = 1:length(RightClustRef_agg)
    A = [XCoordsRef(cell2mat(clustMembsCellRef_agg(j)))'; YCoordsRef(cell2mat(clustMembsCellRef_agg(j)))'];
   
    % Filter particle based on ellipse elongation & major axis length
    ellipse_t = fit_ellipse(A(1,:), A(2,:), gcf, EllipseFitScaleFactor, MaxParticleElongation, MaxDiam, isolatedRef(j));           
    valid = ellipse_t.valid;
    
    if valid == 1
        IndexRightClustRef_elong(j) = true; % Mark this cluster index as valid (1) if it meets the elongation condition  
    else 
        % disp(strcat(['#' num2str(j) ' cluster is elongated and has been excluded'])); 
        r_elong = r_elong + 1; % Increment counter for clusters excluded due to elongation
    end
end

disp(strcat([num2str(r_elong) ' REF clusters have been excluded due to elongation']));

% Update filtered clusters based on elongation
RightClustRef_elong = RightClustRef_agg(IndexRightClustRef_elong);
massCenter_elong = massCenter_agg(:, IndexRightClustRef_elong);
clustMembsCellRef_elong = clustMembsCellRef_agg(IndexRightClustRef_elong);

%%% Filter fiducial markers ===============================================

%%% Filter by localization number in FID channel
disp('Filtering FID channel...');
 RightClustFid=1:1:NClustFid; % Intialization cluster indexes; contains the indices of all FID clusters, from 1 to NClustFid
 massCenterFid=clustCentFid; %Initialization centers coords of clusters; contains the coordinates of the centers of the FID clusters

 IndexRightClustFid=false(NClustFid,1); % Initialization logical array to identify valid clusters
 r_fid=0; % Counter for clusters excluded due to few localizations
 for i = 1:NClustFid
    B = [XCoordsFid(cell2mat(clustMembsCellFid(i)))'; YCoordsFid(cell2mat(clustMembsCellFid(i)))'];
    if(size(B,2) >= MinPts)
        IndexRightClustFid(i) = true; % Mark this cluster index as valid (1) if it meets the condition
    else
        r_fid = r_fid + 1; % Increment counter for clusters excluded due to few localizations
    end
end
disp(strcat( [num2str(r_fid) ' FID clusters have been excluded because of few localizations']));

RightClustFid=RightClustFid(IndexRightClustFid); %indices of filtered FID cluster 
massCenterFid=massCenterFid(:,IndexRightClustFid); %centers of filtered FID clusters
clustMembsCellFid=clustMembsCellFid(IndexRightClustFid); %point indices of filtered FID clusters

% Calculate the distance matrix between the filtered clusters and filtered fiducials
distance_FidRef = pdist2(massCenter_elong', massCenterFid', 'euclidean');

% Find indices of fiducials that are too close to the filtered clusters of NCs
[rowsFid, ~] = find(distance_FidRef < DistFidRef); % Adjust threshold as needed
FidRefCoincidences = rowsFid';

% Initialization of logic index with all values = true
IndexRightClustRef_fid = true(size(RightClustRef_elong));

% Assign false to indexes corresponding to aggregates
IndexRightClustRef_fid(FidRefCoincidences) = false;

% Draw blue ellipse for fiducials
for k = FidRefCoincidences
    A = [XCoordsRef(cell2mat(clustMembsCellRef_elong(k)))'; YCoordsRef(cell2mat(clustMembsCellRef_elong(k)))'];
    % Filter particle based on ellipse elongation & major axis length
    ellipse_t = fit_ellipse_fiducial(A(1,:), A(2,:), gcf, EllipseFitScaleFactor, MaxParticleElongation, MaxDiam, true); 
end

% Remove fiducials from filtered clusters
RightClustRef_fid = RightClustRef_elong(IndexRightClustRef_fid);
massCenter_fid = massCenter_elong(:, IndexRightClustRef_fid);
clustMembsCellRef_fid = clustMembsCellRef_elong(IndexRightClustRef_fid);

disp(strcat( [num2str(length(FidRefCoincidences)) ' clusters have been excluded because they are fiducials']));

disp(strcat(['Identified ' num2str(length(RightClustRef_fid)) ' nanoparticles candidates from ' num2str(NClustRef) ' candidate clusters']));

NPMembs=clustMembsCellRef(RightClustRef_fid); %for every selected cluster which points are in it 
NPSizeRef = cellfun(@numel, NPMembs);  % for each selected cluster, counts the number of points in it

% Continue figure(1)
hold on
plot(massCenter_fid(1,:), massCenter_fid(2,:), 'xk', 'LineWidth',1, 'MarkerSize',10); 
%this is to visualize the label of the identified clusters
[~, colsmC, ~] = size(massCenter_fid);
if colsmC == 1
    txt1 = ['\leftarrow ' num2str(1)];
    text(massCenter_fid(1),massCenter_fid(2),txt1)
else
    for m = 1 : length(massCenter_fid)
        txt1 = ['\leftarrow ' num2str(m)];
        text(massCenter_fid(1,m), massCenter_fid(2,m), txt1)
    end
end

% Start figure(2)

figure(2)
hold on
for w = Aggregates
    A = [XCoordsRef(cell2mat(clustMembsCellRef_loc(w)))'; YCoordsRef(cell2mat(clustMembsCellRef_loc(w)))'];
    plot(A(:,1), A(:,2), 'xy')
end

%%% Size Check ============================================================

DataType='SolidSphere';
nclusters=length(RightClustRef_fid);
%SizeCheck=zeros(nclusters,1);
C=zeros(nclusters,2);
R=zeros(nclusters,1);
Rcheck=zeros(nclusters,1);
RCheckLow=5;
RCheckHigh=400;
discard=0;
IndexToRemove=true(nclusters,1); % logical array to remove index of cluster discarded
for i=1:nclusters
    switch DataType
        %case 'HollowSphere'
           % [C(i,1),C(i,2),R(i)]=circfit(StoreClusterCoords{i,1}(:,1),StoreClusterCoords{i,1}(:,2));
        case 'SolidSphere'
            ClusterData=[XCoordsRef(cell2mat(NPMembs(i))), YCoordsRef(cell2mat(NPMembs(i)))];
            plot(ClusterData(:,1), ClusterData(:,2), 'xr')
            Cinitial=[mean(ClusterData(:,1)) mean(ClusterData(:,2))];
            FracThreshold=0.85;
            CovMat=cov(ClusterData);
            Rinitial=1.5*mean([sqrt(CovMat(1,1)) sqrt(CovMat(2,2))]);
            
            % Now the location of the center of the smallest sphere 
            % encompassing 95% of the datapoints is determined.
            CenterSampleSizeAz=10; CenterSampleSizeRad=4;
            CenterLocation=zeros(CenterSampleSizeAz+1,2,CenterSampleSizeRad);
            TrackLocNumber=zeros(CenterSampleSizeAz+1,CenterSampleSizeRad);
            
            for n=1:CenterSampleSizeRad
                t2 = linspace(0,2*pi,CenterSampleSizeAz);
                XCurrentCircle=0.05*n*Rinitial*cos(t2)+Cinitial(1); 
                YCurrentCircle=0.05*n*Rinitial*sin(t2)+Cinitial(2);
                for o1=1:CenterSampleSizeAz
                        CenterLocation(o1,:,n)=[XCurrentCircle(o1) YCurrentCircle(o1)];
                        LocDistCheck=find(((ClusterData(:,1)-XCurrentCircle(o1)).^2+(ClusterData(:,2)-YCurrentCircle(o1)).^2) < Rinitial^2);
                        TrackLocNumber(o1,n)=length(LocDistCheck);
                end
            end
            
            RefLocCheck=find(((ClusterData(:,1)-Cinitial(1)).^2+(ClusterData(:,2)-Cinitial(2)).^2) < Rinitial^2);
            RefLocNumber=length(RefLocCheck);
            MaxLoc=max(max(TrackLocNumber));
            if MaxLoc > RefLocNumber
                [I,J]=find(TrackLocNumber==MaxLoc);
                MaxLocMinRadTemp=[I J];
                MaxLocMinRadInd=find(J==min(J));
                MaxLocMinRad=MaxLocMinRadTemp(MaxLocMinRadInd,:);
                Cfinal=mean(CenterLocation(MaxLocMinRad(:,1),:,min(J)),1);
            else
                Cfinal=Cinitial;
            end
            
            RadiusVec=linspace(0,3*Rinitial,1500);
            RadiusNumLoc=zeros(length(RadiusVec),1);
            for p=1:length(RadiusVec)
                LocDistCheckRad=find(((ClusterData(:,1)-Cfinal(1)).^2+(ClusterData(:,2)-Cfinal(2)).^2) < RadiusVec(p)^2);
                RadiusNumLoc(p)=length(LocDistCheckRad);
            end
            
            RadTrack=1;
            while RadiusNumLoc(RadTrack) < FracThreshold * size(ClusterData,1)
                RadTrack=RadTrack+1;
            end
            Rfinal=RadiusVec(RadTrack);   
            
            C(i,:)=Cfinal';
            Rcheck(i)=Rfinal;       
    end
    if Rcheck(i) >= RCheckLow && Rcheck(i) <= RCheckHigh
        R(i)=Rfinal;
    else
        disp(['Cluster #' num2str(i) ' has been excluded from the analysis due to an unrealistic size:' num2str(round(Rcheck(i))*2)])
        IndexToRemove(i) = false; % Pietro: collect the indexes of the clusters to remove
        discard=discard+1;
    end
end

disp(strcat(['Identified ' num2str(nclusters-discard) ' valid nanoparticles from ' num2str(nclusters) ' candidate nanoparticles']));
NPselect= R~=0;
R=R(NPselect);
R=round(R);
C=[C(NPselect,1), C(NPselect,2)];
NPSizeRef = NPSizeRef(IndexToRemove);  %Pietro: remove the discarded cluster
%NPMembs=NPMembs(IndexToRemove); %Pietro: remove discarded clusters, after size check

for m=1:length(R)
       t = linspace(0,2*pi,100);
    plot(R(m)*cos(t)+C(m,1),R(m)*sin(t)+C(m,2),'k','LineWidth', 1)
    hold on
%     axis image
    %axis([C(m,1)-500 C(m,1)+500 C(m,2)-500 C(m,2)+500])
end


%%% Identifying localizations in MAIN around each center ==================

disp('Creating output...');

A=[XCoords488 YCoords488 TCoords488];
figure(2)
hold on
plot(XCoords488,YCoords488,'.k'); axis equal; hold on;
xlabel("x (nm)"); ylabel("y (nm)");

[rowsmC, ~, ~] = size(C);
if rowsmC == 1
    Loc2particleMain = cell(1,1); %init: for every REF particle which MAIN-localiz are around, also 0 localizations are included    
    distanceMassCenter = pdist2(C(1,:), [XCoords488 YCoords488]);   % distance between center and all the MAIN-localiz
    Close = (distanceMassCenter <= FactorMaxDist*R); % select MAIN-localiz closer than FactorMaxDist*Radius
    ClosePoint=A(Close,:);                         
    %DistClose=distanceMassCenter(Close);
    %Y = prctile(DistClose,100);                   % option for a further selction of points within a percentage
    %Select=DistClose < Y;
    %DistAttach=DistClose; %(Select);
    Localization=ClosePoint; %(Select,:);
    Loc2particleMain{1} = Localization;
    plot(FactorMaxDist*R*cos(t)+C(1,1),FactorMaxDist*R*sin(t)+C(1,2),'g','LineWidth', 1)
    hold on
    plot(ClosePoint(:,1), ClosePoint(:,2), 'c.'); hold on;
    plot(Localization(:,1), Localization(:,2), 'g.'); hold on;
    %this is to visualize the label of the identified clusters
    txt2 = ['\leftarrow ' num2str(1)];
    plot(C(1,1), C(1,2), 'xk', 'LineWidth',1, 'MarkerSize',5); hold on;
    text(C(1,1),C(1,2),txt2)
else
    Loc2particleMain = cell(length(C),1); %init: for every REF particle which MAIN-localiz are around, also 0 localizations are included
    for k = 1:length(C)
    distanceMassCenter = pdist2(C(k,:), [XCoords488 YCoords488]);   % distance between center and all the MAIN-localiz
    Close = (distanceMassCenter <= FactorMaxDist*R(k)); % select MAIN-localiz closer than FactorMaxDist*Radius
    ClosePoint=A(Close,:);                         
    %DistClose=distanceMassCenter(Close);
    %Y = prctile(DistClose,100);                   % option for a further selction of points within a percentage
    %Select=DistClose < Y;
    %DistAttach=DistClose; %(Select);
    Localization=ClosePoint; %(Select,:);
    Loc2particleMain{k} = Localization;
    plot(FactorMaxDist*R(k)*cos(t)+C(k,1),FactorMaxDist*R(k)*sin(t)+C(k,2),'g','LineWidth', 1)
    hold on
    plot(ClosePoint(:,1), ClosePoint(:,2), 'c.'); hold on;
    plot(Localization(:,1), Localization(:,2), 'g.'); hold on;
    %this is to visualize the label of the identified clusters
    txt2 = ['\leftarrow ' num2str(k)];
    plot(C(k,1), C(k,2), 'xk', 'LineWidth',1, 'MarkerSize',5); hold on;
    text(C(k,1),C(k,2),txt2)
    end

end

title({'Black ellipse: circle fitting';'Red dots: NP, Green dots: BSA attached to NP, Black dots: free BSA'});grid on;axis equal;

% Statistics
ClustSize = cellfun(@numel, Loc2particleMain); %how many MAIN-localiz for particle
ClustSize=ClustSize/3; %IMPORTANT! Because count 1 localization as 3 (3 coordinators)

% Create figure(3)
figure(3);histogram(ClustSize,60);grid on;title('Number of protein localizations/NC');
xlabel("Number of protein localizations"); ylabel("Number of NCs");

%Save data

%save example/Cluster2_quantification_results/Loc2particleMain.txt A -ascii %XYT data of MAIN channel for each selected cluster
save example/Cluster2_quantification_results/CentersRef.txt C -ascii %XY data of centers of REF channel
clustCentFid_T = clustCentFid';
save example/Cluster2_quantification_results/CentersFid.txt clustCentFid_T -ascii %XY data of centers of FID channel


save example/Cluster2_quantification_results/ClustSizeMain.txt ClustSize -ascii %number of protein localizations of each NC
save example/Cluster2_quantification_results/ClustSizeRef.txt NPSizeRef -ascii %number of cy5 localizations of each NC

diam=R*2;
save example/Cluster2_quantification_results/DiametersRef.txt diam -ascii %diameter of NCs

%Save figures

figure(1); 
saveas(gcf, 'example/Cluster2_quantification_results/Fig1_ValidNC.fig');
saveas(gcf, 'example/Cluster2_quantification_results/Fig1_ValidNC.png');

figure(2); 
saveas(gcf, 'example/Cluster2_quantification_results/Fig2_SelectedProt.fig');
saveas(gcf, 'example/Cluster2_quantification_results/Fig2_SelectedProt.png');
 
figure(3); 
saveas(gcf, 'example/Cluster2_quantification_results/Fig3_Histogram.fig');
saveas(gcf, 'example/Cluster2_quantification_results/Fig3_Histogram.png');

clearvars
toc

end
