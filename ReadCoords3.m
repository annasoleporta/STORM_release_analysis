function [XCoords647,YCoords647,TCoords647,XCoords488,YCoords488,TCoords488,XCoordsFid,YCoordsFid,TCoordsFid]=ReadCoords3(FileName, InputType, varargin)

        %%%%%   READCOORDS: reads raw data from SMLM
        %%%%% read raw data (TXT or CSV) obtained from SMLM analysis in
        %%%%% NIS-elements (Nikon N-STORM) or ONI software, and extract
        %%%%% coordinates of interest for further processing. It generates
        %%%%% txt file of 3 columns containing XYT coordinates to use for
        %%%%% further processing.
        %%%%%   
        
        %------------------------------------------------------------------
        % INPUTS:
        % FileName: name of the file with extension, e.g. 'MyFile.txt'
        % InputType: denotes the type of file, type 'N-STORM' for Nikon
        % software TXT files, or 'ONI' for CSV files from ONI.
        % 
        % N.B.
        % Three channels can be read: ref-channel (named 647), main-channel (named 488) and
        % fiducial markers channels (named Fid). File format should be
        % checked.
        %
        % N-STORM files (TXT) are supposed to be 26-column, columns 4-5-13 
        % are read as X-Y-T.        
        % ONI files (CSV) are supposed to be 11-columns, columns 3-4-2 are
        % read as X-Y-T. 
        %
        %------------------------------------------------------------------
        % OUTPUTS:
        % X-Y-TCoords647: X, Y, T(frames) coordinates of localization in
        % the ref channel (NCs), named 647, expressed in nanometers
        %
        % X-Y-TCoords488: X, Y, T(frames) coordinates of localization in
        % the main channel (protein), named 488, expressed in nanometers
        %
        % X-Y-TCoordsFid: X, Y, T(frames) coordinates of localization in
        % the fiducial markers channel, named Fid, expressed in
        % nanometers
        % 
        % 
        %------------------------------------------------------------------
        
p = inputParser; %init parser object
validChar = @(x) ischar(x);
validNum = @(x) isreal(x);
% here add something to read the coordinate file


%define defaults values for optional param:
defaultSTORMref = '405/647'; % ref channel name for N-STORM data
defaultSTORMmain = '405/488'; % main channel name for N-STORM data
defaultSTORMfid = 'Bead Drift Correction'; % fid channel name for N-STORM data
defaultONIref = 1 ; % ref channel name for ONI
defaultONImain = 2 ; % main channel name for ONI
defaultONIfid = 0 ; % fid channel name for ONI

%define required and optional input parameters:
addRequired(p,'FileName',validChar);
addRequired(p,'InputType',validChar); 

addParameter(p,'STORMref', defaultSTORMref, validChar);
addParameter(p,'STORMmain', defaultSTORMmain, validChar);
addParameter(p,'STORMfid', defaultSTORMfid, validChar);
addParameter(p,'ONIref', defaultONIref, validNum);
addParameter(p,'ONImain', defaultONImain, validNum);
addParameter(p,'ONIfid', defaultONIfid, validNum);

%read input values:
parse(p, FileName, InputType, varargin{:});

%assign the parsed values:
FileName = p.Results.FileName; % file name
InputType = p.Results.InputType;
STORMref = p.Results.STORMref;
STORMmain = p.Results.STORMmain;
STORMfid = p.Results.STORMfid;
ONIref = p.Results.ONIref;
ONImain = p.Results.ONImain;
ONIfid = p.Results.ONIfid;

switch InputType
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%------------------ TXT file from N-STORM ---------------------%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    case 'N-STORM'

        %%%----- Read text file (check format)---
        disp('Importing N-STORM data...');
        fileID = fopen(FileName,'r');
        DataIn = textscan(fileID,'%s %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f','Delimiter','\t','HeaderLines',1); % Saving data as a cell, check the numb of columns 
        fclose(fileID);

        %%%--- Retrieve useful information (channel, corrected XY location and frame number or the localization)
        XCoords = cell2mat(DataIn(4)); %% colum 2 X non drift corrected, 4 X drift corrected
        YCoords = cell2mat(DataIn(5)); %% colum 3 Y non drift corrected, 5 Y drift corrected
        TCoords = cell2mat(DataIn(13)); %% column 13 Frame
        Channel = DataIn(1);
        clear DataIn;

        %%%--- Split channel based on first column and generate output
        Channel647=strcmp(Channel{1},STORMref);%Check in the .txt if 647 or 405/647
        XCoords647=XCoords(Channel647);
        YCoords647=YCoords(Channel647);
        TCoords647=TCoords(Channel647);
        Channel488=strcmp(Channel{1},STORMmain);%Check in the .txt if 488 or 405/488
        XCoords488=XCoords(Channel488);
        YCoords488=YCoords(Channel488);
        TCoords488=TCoords(Channel488);
        ChannelFid=strcmp(Channel{1},STORMfid);%Check in the .txt the name
        XCoordsFid=XCoords(ChannelFid);
        YCoordsFid=YCoords(ChannelFid);
        TCoordsFid=TCoords(ChannelFid);
        
        %%%--- Export a txt file with XYT coords
        nref=length(XCoords647);
        nmain=length(XCoords488);
        nfid=length(XCoordsFid);
        Ref=zeros(nref,3);
        Main=zeros(nmain,3);
        Fid=zeros(nfid,3);
        Ref(:,1)=XCoords647;
        Ref(:,2)=YCoords647;
        Ref(:,3)=TCoords647;
        Main(:,1)=XCoords488;
        Main(:,2)=YCoords488;
        Main(:,3)=TCoords488;
        Fid(:,1)=XCoordsFid;
        Fid(:,2)=YCoordsFid;
        Fid(:,3)=TCoordsFid;
        save example/ReadCoords3_results/XYTref.txt Ref -ascii
        save example/ReadCoords3_results/XYTcoordinates.txt Main -ascii
        save example/ReadCoords3_results/XYTfid.txt Fid -ascii
        
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%--------------------- CSV file from ONI ----------------------%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    case 'ONI'
        
        %%%----- Read text file (check format)---
        disp('Importing ONI data...');
        fileID = fopen(FileName,'r');
        DataIn = textscan(fileID,'%f %f %f %f %f %f %f %f %f %f %f','Delimiter',',','HeaderLines',1); % Saving data as a cell, check the numb of columns 
        fclose(fileID);
        
        %%%--- Retrieve useful information (channel, corrected XY location and frame number or the localization)
        XCoords = cell2mat(DataIn(3)); %% colum 3 X coords
        YCoords = cell2mat(DataIn(4)); %% colum 4 Y coords
        TCoords = cell2mat(DataIn(2)); %% column 2 Frame coords
        Channel = DataIn(1);
        clear DataIn;
        
        
        %%%--- Split channel based on first column and generate output
        Channel647= (Channel{1}==ONIref); % 647-channel denoted with 1
        XCoords647=XCoords(Channel647);
        YCoords647=YCoords(Channel647);
        TCoords647=TCoords(Channel647);
        Channel488= (Channel{1}==ONImain); % 488-channel denoted with 2
        XCoords488=XCoords(Channel488);
        YCoords488=YCoords(Channel488);
        TCoords488=TCoords(Channel488);
        ChannelFid= (Channel{1}==ONIfid);%C Fiducial-channels denoted with 0
        XCoordsFid=XCoords(ChannelFid);
        YCoordsFid=YCoords(ChannelFid);
        TCoordsFid=TCoords(ChannelFid);
        
        %%%--- Export a txt file with XYT coords
        nref=length(XCoords647);
        nmain=length(XCoords488);
        nfid=length(XCoordsFid);
        Ref=zeros(nref,3);
        Main=zeros(nmain,3);
        Fid=zeros(nfid,3);
        Ref(:,1)=XCoords647;
        Ref(:,2)=YCoords647;
        Ref(:,3)=TCoords647;
        Main(:,1)=XCoords488;
        Main(:,2)=YCoords488;
        Main(:,3)=TCoords488;
        Fid(:,1)=XCoordsFid;
        Fid(:,2)=YCoordsFid;
        Fid(:,3)=TCoordsFid;
        save example/ReadCoords3_results/XYTref.txt Ref -ascii
        save example/ReadCoords3_results/XYTcoordinates.txt Main -ascii
        save example/ReadCoords3_results/XYTfid.txt Fid -ascii
        

    otherwise
        disp('invalid InputType!');
end
clearvars

disp('Coordinates file generated :)');
end