% code to read and assemble the downloaded CMEMS (AVISO) along track data
% this code assumes you are on a UNIX-like system

% this assume you have downloaded the global, delayed-time, along-track, 
% unfiltered, (absolute dynamic topography data from AVISO) Sea Level 
% Anomalies from 
% CMEMS (https://resources.marine.copernicus.eu/?option=com_csw&view=details&product_id=SEALEVEL_GLO_PHY_L3_REP_OBSERVATIONS_008_062)

% data for Jason 3, 2, 1, and TOPEX/Poseidon are organized in subfolders by years/month/; individual files are zipped with extension .gz
% are organized by individual netcdf files with extension .nc
% data for Jason 2 are organized in subfolders by years; individual netcdf files with extension .nc
% the data are now distributed through COPERNICUS MARINE ENVIRONMENT MONITORING SERVICE
% http://marine.copernicus.edu 
% SEALEVEL_GLO_PHY_L3_REP_OBSERVATIONS_008_062
% the naming of files may have changed compared to what is shown here 

% please edit the directory names below to match your file system

% basedir = 'aviso/global/delayed-time/along-track/unfiltered/adt/tp/';
basedir = '/Users/gnovelli/Desktop/sla-cmems/tp/';

n = 0;

% netvar = {'ADT','cycle','latitude','longitude','time'};
netvar = {'sla_unfiltered','mdt','cycle','latitude','longitude','time'};

%
% for year = 1992:2002;
for year = 1993:2002
    for month = 1:12
        
    filelist = dir([basedir num2str(year) '/' num2str(month,'%02d')]);
    % remove directories
    q = find([filelist.isdir] == 1);
    filelist(q) = [];
        
    for k = 1:length(filelist)
        
            eval(['!cp ' basedir num2str(year) '/' num2str(month,'%02d') '/' filelist(k).name ' ./foo.nc']);
%             eval(['!gunzip -f foo.nc.gz']);
            pause(.1);
            track = ncread(['foo.nc'],'track');
            % should load ADT, cycle, latitude, longitude, time, track, all of same dimension
            q96 = find(track==96);
            if ~isempty(q96)
                n = n+1;
                for m = 1:length(netvar)
                    eval([netvar{m}  '= ncread([''foo.nc''],''' netvar{m} ''');']);
                end
%                 adt96_tp_vxxc(n).ADT = ADT(q96);
                adt96_tp_vxxc(n).ADT = sla_unfiltered(q96)+mdt(q96);
                adt96_tp_vxxc(n).cycle = cycle(q96);
                adt96_tp_vxxc(n).latitude = latitude(q96);
                adt96_tp_vxxc(n).longitude = longitude(q96);
                adt96_tp_vxxc(n).time = time(q96);
                adt96_tp_vxxc(n).track = track(q96);
                adt96_tp_vxxc(n).file = filelist(k).name;
                disp(n);
                sprintf('%d',n);
                clear sla_unfiltered mdt cycle latitude longitude time track q96
            end
            eval(['!rm foo.nc']);
    end
    end
end
%    
% jason 1

% basedir = 'aviso/global/delayed-time/along-track/unfiltered/adt/j1/';
basedir = '/Users/gnovelli/Desktop/sla-cmems/j1/';
n = 0;

% for year = 2002:2008;
for year = 2002:2008
    for month = 1:12
        
    filelist = dir([basedir num2str(year) '/' num2str(month,'%02d')]);
    % remove directories
    q = find([filelist.isdir] == 1);
    filelist(q) = [];
        
    for k = 1:length(filelist)
        
            eval(['!cp ' basedir num2str(year) '/' num2str(month,'%02d') '/' filelist(k).name ' ./foo.nc']);
%             eval(['!gunzip -f foo.nc.gz']);
            pause(.1);
            track = ncread(['foo.nc'],'track');
            % should load ADT, cycle, latitude, longitude, time, track, all of same dimension
            q96 = find(track==96);
            if ~isempty(q96)
                n = n+1;
                for m = 1:length(netvar)
                    eval([netvar{m}  '= ncread([''foo.nc''],''' netvar{m} ''');']);
                end
                adt96_j1_vxxc(n).ADT = sla_unfiltered(q96)+mdt(q96);
                adt96_j1_vxxc(n).cycle = cycle(q96);
                adt96_j1_vxxc(n).latitude = latitude(q96);
                adt96_j1_vxxc(n).longitude = longitude(q96);
                adt96_j1_vxxc(n).time = time(q96);
                adt96_j1_vxxc(n).track = track(q96);
                adt96_j1_vxxc(n).file = filelist(k).name;
                disp(n);
                sprintf('%d',n);
                clear sla_unfiltered mdt cycle latitude longitude time track q96
            end
            eval(['!rm foo.nc']);
    end
    end
end


% jason 2

% basedir = 'aviso/global/delayed-time/along-track/unfiltered/adt/j2/';
basedir = '/Users/gnovelli/Desktop/sla-cmems/j2/';
n = 0;

% for year = 2008:2016;
for year = 2008:2016
    for month = 1:12
        
    filelist = dir([basedir num2str(year) '/' num2str(month,'%02d')]);
    % remove directories
    q = find([filelist.isdir] == 1);
    filelist(q) = [];
        
    for k = 1:length(filelist)
        
            eval(['!cp ' basedir num2str(year) '/' num2str(month,'%02d') '/' filelist(k).name ' ./foo.nc']);
%             eval(['!gunzip -f foo.nc.gz']);
            pause(.1);
            track = ncread(['foo.nc'],'track');
            % should load ADT, cycle, latitude, longitude, time, track, all of same dimension
            q96 = find(track==96);
            if ~isempty(q96)
                n = n+1;
                for m = 1:length(netvar)
                    eval([netvar{m}  '= ncread([''foo.nc''],''' netvar{m} ''');']);
                end
                adt96_j2_vxxc(n).ADT =  sla_unfiltered(q96)+mdt(q96);
                adt96_j2_vxxc(n).cycle = cycle(q96);
                adt96_j2_vxxc(n).latitude = latitude(q96);
                adt96_j2_vxxc(n).longitude = longitude(q96);
                adt96_j2_vxxc(n).time = time(q96);
                adt96_j2_vxxc(n).track = track(q96);
                adt96_j2_vxxc(n).file = filelist(k).name;
                disp(n);
                sprintf('%d',n);
                clear sla_unfiltered mdt cycle latitude longitude time track q96
            end
            eval(['!rm foo.nc']);
    end
    end
    
end
%
% jason 3

% basedir = 'aviso/global/delayed-time/along-track/unfiltered/adt/j2/';
basedir = '/Users/gnovelli/Desktop/sla-cmems/j3/';
n = 0;


for year = 2016:2020
    for month = 1:12
        
    filelist = dir([basedir num2str(year) '/' num2str(month,'%02d')]);
    % remove directories
    q = find([filelist.isdir] == 1);
    filelist(q) = [];
        
    for k = 1:length(filelist)
        
            eval(['!cp ' basedir num2str(year) '/' num2str(month,'%02d') '/' filelist(k).name ' ./foo.nc']);
%             eval(['!gunzip -f foo.nc.gz']);
            pause(.1);
            track = ncread(['foo.nc'],'track');
            % should load ADT, cycle, latitude, longitude, time, track, all of same dimension
            q96 = find(track==96);
            if ~isempty(q96)
                n = n+1;
                for m = 1:length(netvar)
                    eval([netvar{m}  '= ncread([''foo.nc''],''' netvar{m} ''');']);
                end
                adt96_j3_vxxc(n).ADT =  sla_unfiltered(q96)+mdt(q96);
                adt96_j3_vxxc(n).cycle = cycle(q96);
                adt96_j3_vxxc(n).latitude = latitude(q96);
                adt96_j3_vxxc(n).longitude = longitude(q96);
                adt96_j3_vxxc(n).time = time(q96);
                adt96_j3_vxxc(n).track = track(q96);
                adt96_j3_vxxc(n).file = filelist(k).name;
                disp(n);
                sprintf('%d',n);
                clear sla_unfiltered mdt cycle latitude longitude time track q96
            end
            eval(['!rm foo.nc']);
    end
    end
    
end

save('adt96_tp.mat','adt96_tp_vxxc');
save('adt96_j1.mat','adt96_j1_vxxc');
save('adt96_j2.mat','adt96_j2_vxxc');
save('adt96_j3.mat','adt96_j3_vxxc');

