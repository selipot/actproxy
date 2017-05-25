% code to read and assemble the downloaded AVISO along track data
% this code assumes you are on a UNIX-like system

% this assume you have downloaded the global, delayed-time, along-track, unfiltered, absolute dynamic topography data from AVISO

% data for TOPEX/Poseidon are organized in subfolders by years; individual files are zipped with extension .gz

% note that in 2017 the data are now distributed through COPERNICUS MARINE ENVIRONMENT MONITORING SERVICE
% http://marine.copernicus.edu 
% the naming of files may have changed compared to what is shown here 

% please edit the directory names below to match your file system

basedir = 'aviso/global/delayed-time/along-track/unfiltered/adt/tp/';

n = 0;

netvar = {'ADT','cycle','latitude','longitude','time'};

for year = 1992:2002;
    
    filelist = dir([basedir num2str(year)]);
    % remove directories
    q = find([filelist.isdir] == 1);
    filelist(q) = [];
        
    for k = 1:length(filelist)
        
            eval(['!cp ' basedir num2str(year) '/' filelist(k).name ' ./foo.nc.gz']);
            eval(['!gunzip -f foo.nc.gz']);
            pause(.1);
            track = ncread(['foo.nc'],'track');
            % should load ADT, cycle, latitude, longitude, time, track, all of same dimension
            q96 = find(track==96);
            if ~isempty(q96)
                n = n+1;
                for m = 1:length(netvar)
                    eval([netvar{m}  '= ncread([''foo.nc''],''' netvar{m} ''');']);
                end
                adt96_tp_vxxc(n).ADT = ADT(q96);
                adt96_tp_vxxc(n).cycle = cycle(q96);
                adt96_tp_vxxc(n).latitude = latitude(q96);
                adt96_tp_vxxc(n).longitude = longitude(q96);
                adt96_tp_vxxc(n).time = time(q96);
                adt96_tp_vxxc(n).track = track(q96);
                adt96_tp_vxxc(n).file = filelist(k).name;
                disp(n);
                sprintf('%d',n);
                clear ADT cycle latitude longitude time track q96
            end
            eval(['!rm foo.nc']);
    end
end
    
% jason 1

basedir = 'aviso/global/delayed-time/along-track/unfiltered/adt/j1/';

n = 0;

for year = 2002:2008;
    
    filelist = dir([basedir num2str(year)]);
    % remove directories
    q = find([filelist.isdir] == 1);
    filelist(q) = [];
        
    for k = 1:length(filelist)
        
            eval(['!cp ' basedir num2str(year) '/' filelist(k).name ' ./foo.nc.gz']);
            eval(['!gunzip -f foo.nc.gz']);
            pause(.1);
            track = ncread(['foo.nc'],'track');
            % should load ADT, cycle, latitude, longitude, time, track, all of same dimension
            q96 = find(track==96);
            if ~isempty(q96)
                n = n+1;
                for m = 1:length(netvar)
                    eval([netvar{m}  '= ncread([''foo.nc''],''' netvar{m} ''');']);
                end
                adt96_j1_vxxc(n).ADT = ADT(q96);
                adt96_j1_vxxc(n).cycle = cycle(q96);
                adt96_j1_vxxc(n).latitude = latitude(q96);
                adt96_j1_vxxc(n).longitude = longitude(q96);
                adt96_j1_vxxc(n).time = time(q96);
                adt96_j1_vxxc(n).track = track(q96);
                adt96_j1_vxxc(n).file = filelist(k).name;
                disp(n);
                sprintf('%d',n);
                clear ADT cycle latitude longitude time track q96
            end
            eval(['!rm foo.nc']);
    end
end


% jason 2

basedir = 'aviso/global/delayed-time/along-track/unfiltered/adt/j2/';

n = 0;

for year = 2008:2016;
    
    filelist = dir([basedir num2str(year)]);
    % remove directories
    q = find([filelist.isdir] == 1);
    filelist(q) = [];
        
    for k = 1:length(filelist)
        
            eval(['!cp ' basedir num2str(year) '/' filelist(k).name ' ./foo.nc.gz']);
            eval(['!gunzip -f foo.nc.gz']);
            pause(.1);
            track = ncread(['foo.nc'],'track');
            % should load ADT, cycle, latitude, longitude, time, track, all of same dimension
            q96 = find(track==96);
            if ~isempty(q96)
                n = n+1;
                for m = 1:length(netvar)
                    eval([netvar{m}  '= ncread([''foo.nc''],''' netvar{m} ''');']);
                end
                adt96_j2_vxxc(n).ADT = ADT(q96);
                adt96_j2_vxxc(n).cycle = cycle(q96);
                adt96_j2_vxxc(n).latitude = latitude(q96);
                adt96_j2_vxxc(n).longitude = longitude(q96);
                adt96_j2_vxxc(n).time = time(q96);
                adt96_j2_vxxc(n).track = track(q96);
                adt96_j2_vxxc(n).file = filelist(k).name;
                disp(n);
                sprintf('%d',n);
                clear ADT cycle latitude longitude time track q96
            end
            eval(['!rm foo.nc']);
    end
    
end

save('adt96_tp.mat','adt96_tp_vxxc');
save('adt96_j1.mat','adt96_j1_vxxc');
save('adt96_j2.mat','adt96_j2_vxxc');

