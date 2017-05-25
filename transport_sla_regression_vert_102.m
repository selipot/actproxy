%% code to derive and update the proxy time series of Agulhas jet and boundary transports as in Beal and Elipot (2016) (doi:10.1038/nature19853)

 latM = [-33.3144,... % shoreline
        -33.5564,... % Mooring A
        -33.6536,... % Mooring B
        -33.7823,... % Mooring C  
        -34.0207,... % Mooring D
        -34.2865,... % Mooring E
        -34.5387,... % Mooring F
        -34.8219,... % Mooring G
     -34.959,... % CPIES P3
      -35.346,...% CPIES P4
      -35.734,...% CPIES P5
     ]';

lonM = [27.4426,... %
        27.5986,... % 
        27.6579,... % 
        27.7198,... % 
        27.8633,... % 
        28.0259,... % 
        28.1622,... % 
        28.3417,... % 
        28.427 ,... %
        28.661,... %
        28.897]';%

% convert the coordinates onto coordinates on plane tangent to some point
lato = 0.5*sum(latM([1 end]));
lono = 0.5*sum(lonM([1 end]));

% projection of mooring positions
[xM,yM,dM] = latlon2xy(latM,lonM,lato,lono);
dM(xM<0) = -dM(xM<0);

% let's define an axis of projection as a line going from the coast to P5
xy12 = latlon2xy(latM([1 end]),lonM([1 end]),lato,lono);
% y = a x +b
a = (imag(xy12(2))-imag(xy12(1)))/(real(xy12(2))-real(xy12(1)));
b = imag(xy12(1)) - real(xy12(1))*(imag(xy12(2))-imag(xy12(1)))/(real(xy12(2))-real(xy12(1)));

% projection of mooring positions on axis 
xM2 = (xM+a*(yM-b))./(1+a.^2);
yM2 = a*xM2+b;
dM2 = abs(xM2+1i*yM2);
dM2(xM2<0) = -dM2(xM2<0);

dM3(1:7) = dM2(2:8);% mooring A to G
dM3(8) = 0.5*(dM2(9)+dM2(10)); %P3-P4
dM3(9) = 0.5*(dM2(10)+dM2(11)); %P3-P4

% load the assembled sla data from the read_alongtrack.m code
load('adt96_tp.mat','adt96_tp_vxxc');
load('adt96_j1.mat','adt96_j1_vxxc');
load('adt96_j2.mat','adt96_j2_vxxc');

% convert int16 to double
for k = 1:length(adt96_tp_vxxc)
    adt96_tp_vxxc(k).track = double(adt96_tp_vxxc(k).track);
    adt96_tp_vxxc(k).cycle = double(adt96_tp_vxxc(k).cycle);
end
for k = 1:length(adt96_j1_vxxc)
    adt96_j1_vxxc(k).track = double(adt96_j1_vxxc(k).track);
    adt96_j1_vxxc(k).cycle = double(adt96_j1_vxxc(k).cycle);
end
for k = 1:length(adt96_j2_vxxc)
    adt96_j2_vxxc(k).track = double(adt96_j2_vxxc(k).track);
    adt96_j2_vxxc(k).cycle = double(adt96_j2_vxxc(k).cycle);
end
    
cycle_tp = cell2col({adt96_tp_vxxc.cycle});
cycle_j1 = 1000 + cell2col({adt96_j1_vxxc.cycle});
cycle_j2 = 2000 + cell2col({adt96_j2_vxxc.cycle});

ADT_tp = cell2col({adt96_tp_vxxc.ADT});
latitude_tp = cell2col({adt96_tp_vxxc.latitude});
longitude_tp = cell2col({adt96_tp_vxxc.longitude});
time_tp = datenum(1950,1,1) + cell2col({adt96_tp_vxxc.time});

ADT_j1 = cell2col({adt96_j1_vxxc.ADT});
latitude_j1 = cell2col({adt96_j1_vxxc.latitude});
longitude_j1 = cell2col({adt96_j1_vxxc.longitude});
time_j1 = datenum(1950,1,1) + cell2col({adt96_j1_vxxc.time});

ADT_j2 = cell2col({adt96_j2_vxxc.ADT});
latitude_j2 = cell2col({adt96_j2_vxxc.latitude});
longitude_j2 = cell2col({adt96_j2_vxxc.longitude});
time_j2 = datenum(1950,1,1) + cell2col({adt96_j2_vxxc.time});

ADT = [ADT_tp;ADT_j1;ADT_j2];
clear ADT_tp ADT_j1 ADT_j2;
latitude = [latitude_tp;latitude_j1;latitude_j2];
clear latitude_tp latitude_j1 latitude_j2;
longitude = [longitude_tp;longitude_j1;longitude_j2];
clear longitude_tp longitude_j1 longitude_j2;
time = [time_tp;time_j1;time_j2];
clear time_tp time_j1 time_j2;

cycle = [cycle_tp ; cycle_j1 ; cycle_j2];
ncycle = unique(cycle(~isnan(cycle)));%

for k = 1:length(ncycle)
    n = ncycle(k);
    q = find(cycle==n);
    adt96_vxxc(k).ADT = ADT(q);
    adt96_vxxc(k).latitude = latitude(q);
    adt96_vxxc(k).longitude = longitude(q);
    adt96_vxxc(k).time = time(q);
    adt96_vxxc(k).cycle = n;
end

clear cycle_*
clear adt96_tp_vxxc adt96_j1_vxxc adt96_j2_vxxc
clear ADT cycle latitude longitude time

ncycle_adt = ncycle;
clear ncycle

% projection onto a straight line going through the array
for k = 1:length(ncycle_adt)
    q = find(adt96_vxxc(k).latitude<0);
    if ~isempty(q)
        adt96_vxxc(k).qxy = q;
        [x,y,d] = latlon2xy(adt96_vxxc(k).latitude(q),adt96_vxxc(k).longitude(q),lato,lono);
        adt96_vxxc(k).xiy = x+1i*y;
        % calculate the projection of the points perpendicularly onto a straight line
        % based on minimization of distance squared and equation of straight line
        x2 = (x+a*(y-b))./(1+a.^2);
        y2 = a*x2+b;
        adt96_vxxc(k).xiy2 = x2+1i*y2;
        % distance axis from origin (0,0) in km
        xd = abs(x2+1i*y2);
        xd(x2<0) = -xd(x2<0);
        adt96_vxxc(k).xd = xd;
    end
end

% calculate for each mooring the time series of SLA and slope
% length scales over which to estimate slope
hh_a = [27.8127426375306 26.9926868661308 49.3793997130902 84.911738809326 89.5898643727928 ...
    67.0051815604118 101.715196175176 93.0415542036546 84.4151849397396];

% recalculate all the slopes
if 1
    
    f = 2*pi*phi2f(lato)/(60*60*24);
    fac = 9.81/f/10^3;% because gradient is in m per km
    
    % default is quadratic but also try linear and cubic
    p = 2;

    H = NaN*ones(length(adt96_vxxc),length(dM3),2);
    V = NaN*ones(length(adt96_vxxc),length(dM3),2);
    A = NaN*ones(length(adt96_vxxc),length(dM3),2);
    
    
    for k = 1:length(adt96_vxxc)
        %sprintf('%s',num2str(k));
        %disp(num2str(k));
        di = adt96_vxxc(k).xd;
        zi = adt96_vxxc(k).ADT(adt96_vxxc(k).qxy);
       
        % coefficient of polynomial
        if ~isempty(di)
            
            if 1
                % from ADT
                cp1 = NaN*ones(9,p-1+1);
                cp2 = NaN*ones(9,p+1);
                
                for m = 1:9
                        cp2(m,:) = LocalPolyFit(di,zi,dM3(m),p,hh_a(m));
                        cp1(m,:) = LocalPolyFit(di,zi,dM3(m),p-1,hh_a(m));
                end
                % ssh
                H(k,:,1) = cp1(:,1);
                H(k,:,2) = cp2(:,1);
                %slope
                V(k,:,1) = cp1(:,2);
                V(k,:,2) = cp2(:,2);
                % "acceleration"
                A(k,:,1) = 0;
                A(k,:,2) = cp2(:,3);
                
            end
            
        end
        
    end
    
end

% define the altimetry times
M = length(adt96_vxxc);
talt = NaN*ones(length(adt96_vxxc),1);
for k = 1:M;
    if ~isempty([adt96_vxxc(k).qxy])
        talt(k) = adt96_vxxc(k).time(adt96_vxxc(k).qxy(1));
    end
end
q = find(isnan(talt));
talt(q) = 0.5*(talt(q-1)+talt(q+1));

% regression coefficients to convert slope to transport

% slope from ADT
% linear for A, B, quadratic for others
s = cat(2,squeeze(V(:,1,1)),squeeze(V(:,2,1)),...
    squeeze(V(:,3,2)),squeeze(V(:,4,2)),...
    squeeze(V(:,5,2)),squeeze(V(:,6,2)),...
    squeeze(V(:,7,2)),squeeze(V(:,8,2)),...
    squeeze(V(:,9,2)));

% SLA from ADT
sla = cat(2,squeeze(H(:,1,1)),squeeze(H(:,2,1)),...
    squeeze(H(:,3,2)),squeeze(H(:,4,2)),...
    squeeze(H(:,5,2)),squeeze(H(:,6,2)),...
    squeeze(H(:,7,2)),squeeze(H(:,8,2)),...
    squeeze(H(:,9,2)));


% FIRST remove outliers in s and sla
s2 = s;
sla2 = sla;

if 1
    for k = [1 2 3 4 5 6]
        q = find(s(:,k)>prctile(s(:,k),99.75));
        s2(q,k) = NaN;
        sla2(q,k) = NaN;
        q = find(s(:,k)<prctile(s(:,k),.25));
        s2(q,k) = NaN;
        sla2(q,k) = NaN;
    end
end

% how many gaps per mooring
[qt2,qm2] = find(isnan(s2));% 
for k = 1:9
    g(k) = sum(qm2==k);
end

% fill data gaps, one mooring at the time
qall = find(all(~isnan([s2 sla2]).'));% times when all data are available
Y = [s2(qall,:) sla2(qall,:)];
clear xstats2 xstats
for m = 1:18
    x = Y(:,m);
    if m < 10
        index = setdiff(1:18,[m m+9]); %can't use sla at same location
    else
        index = setdiff(1:18,[m m-9]); %can't use slope at same location
    end
    y = Y(:,index);
    xstats(m) = regstats(x,y);
    q = find(xstats(m).tstat.pval(2:end)<0.01);
    xstats2{m} = regstats(x,y(:,q));
    xstats2{m}.q = index(q);
    disp([xstats(m).rsquare xstats2{m}.rsquare]);
end

s3 = s2;
sla3 = sla2;
Y = [s2 sla2];
for m = 1:9
    yr = Y(:,xstats2{m}.q);
    qnan =find(isnan(Y(:,m)) & all(~isnan(yr).').');
    xr = yr(qnan,:)*xstats2{m}.beta(2:end);
    s3(qnan,m) = xr; 
    qnan0 = find(isnan(Y(:,m)));
    disp(['filled ' num2str(length(qnan)) ' slope points out of ' num2str(length(qnan0)) ' at mooring ' num2str(m)]);
end

for m = 10:18
    yr = Y(:,xstats2{m}.q);
    qnan =find(isnan(Y(:,m)) & all(~isnan(yr).').');
    xr = yr(qnan,:)*xstats2{m}.beta(2:end);
    sla3(qnan,m-9) = xr; 
    qnan0 = find(isnan(Y(:,m)));
    disp(['filled ' num2str(length(qnan)) ' sla points out of ' num2str(length(qnan0)) ' at mooring ' num2str(m)]);
end

% next step, iterative procedure, starting from offshore
% #6 from offshore
clear xstats2 xstats
m = [7:9];
for k = 6;
    qall = find(all(~isnan([s2(:,[k m]) sla2(:,m)]).'));% times when all data are available %can't use sla at same location
    y = [s2(qall,m) sla2(qall,m)];
    x = s2(qall,k);
    xstats(k) = regstats(x,y);
    q = find(xstats(k).tstat.pval(2:end)<0.01);
    xstats2{k} = regstats(x,y(:,q));
    index = [m m+9];
    xstats2{k}.q = index(q);
    %disp(xstats2{k}.rsquare);
end

for k = 6;
    qall = find(all(~isnan([s2(:,[m]) sla2(:,[k m])]).'));% times when all data are available %can't use sla slope at same location
    y = [s2(qall,m) sla2(qall,m)];
    x = sla2(qall,k);
    xstats(k+9) = regstats(x,y);
    q = find(xstats(k+9).tstat.pval(2:end)<0.01);
    xstats2{k+9} = regstats(x,y(:,q));
    index = [m m+9];
    xstats2{k+9}.q = index(q);
    %disp(xstats2{k+9}.rsquare);
end

s4 = s3;
sla4 = sla3;
Y = [s3 sla3];
for k = 6
    yr = Y(:,xstats2{k}.q);
    qnan =find(isnan(s3(:,k)) & all(~isnan(yr).').');
    xr = yr(qnan,:)*xstats2{k}.beta(2:end);
    s4(qnan,k) = xr; 
    qnan0 = find(isnan(s3(:,k)));
    disp(['filled ' num2str(length(qnan)) ' slope points out of ' num2str(length(qnan0)) ' at mooring ' num2str(k)]);
end
for k = 6
    yr = Y(:,xstats2{k+9}.q);
    qnan =find(isnan(sla3(:,k)) & all(~isnan(yr).').');
    xr = yr(qnan,:)*xstats2{k+9}.beta(2:end);
    sla4(qnan,k) = xr; 
    qnan0 = find(isnan(sla3(:,k)));
    disp(['filled ' num2str(length(qnan)) ' sla points out of ' num2str(length(qnan0)) ' at mooring ' num2str(k)]);
end

% #5 from offshore
clear xstats2 xstats
m = [6:9];
for k = 5;
    qall = find(all(~isnan([s2(:,[k m]) sla2(:,m)]).'));% times when all data are available %can't use sla at same location
    y = [s2(qall,m) sla2(qall,m)];
    x = s2(qall,k);
    xstats(k) = regstats(x,y);
    q = find(xstats(k).tstat.pval(2:end)<0.01);
    xstats2{k} = regstats(x,y(:,q));
    index = [m m+9];
    xstats2{k}.q = index(q);
    disp(xstats2{k}.rsquare);
end
for k = 5;
    qall = find(all(~isnan([s2(:,[m]) sla2(:,[k m])]).'));% times when all data are available %can't use sla slope at same location
    y = [s2(qall,m) sla2(qall,m)];
    x = sla2(qall,k);
    xstats(k+9) = regstats(x,y);
    q = find(xstats(k+9).tstat.pval(2:end)<0.01);
    xstats2{k+9} = regstats(x,y(:,q));
    index = [m m+9];
    xstats2{k+9}.q = index(q);
    disp(xstats2{k+9}.rsquare);
end

s5 = s4;
sla5 = sla4;
Y = [s4 sla4];
for k = 5
    yr = Y(:,xstats2{k}.q);
    qnan =find(isnan(s4(:,k)) & all(~isnan(yr).').');
    xr = yr(qnan,:)*xstats2{k}.beta(2:end);
    s5(qnan,k) = xr; 
    qnan0 = find(isnan(s4(:,k)));
    disp(['filled ' num2str(length(qnan)) ' slope points out of ' num2str(length(qnan0)) ' at mooring ' num2str(k)]);
end
for k = 5
    yr = Y(:,xstats2{k+9}.q);
    qnan =find(isnan(sla4(:,k)) & all(~isnan(yr).').');
    xr = yr(qnan,:)*xstats2{k+9}.beta(2:end);
    sla5(qnan,k) = xr; 
    qnan0 = find(isnan(sla4(:,k)));
    disp(['filled ' num2str(length(qnan)) ' sla points out of ' num2str(length(qnan0)) ' at mooring ' num2str(k)]);
end

% #4 from offshore
clear xstats2 xstats
m = [5:9];
for k = 4;
    qall = find(all(~isnan([s2(:,[k m]) sla2(:,m)]).'));% times when all data are available %can't use sla at same location
    y = [s2(qall,m) sla2(qall,m)];
    x = s2(qall,k);
    xstats(k) = regstats(x,y);
    q = find(xstats(k).tstat.pval(2:end)<0.01);
    xstats2{k} = regstats(x,y(:,q));
    index = [m m+9];
    xstats2{k}.q = index(q);
    disp(xstats2{k}.rsquare);
end
for k = 4;
    qall = find(all(~isnan([s2(:,[m]) sla2(:,[k m])]).'));% times when all data are available %can't use sla slope at same location
    y = [s2(qall,m) sla2(qall,m)];
    x = sla2(qall,k);
    xstats(k+9) = regstats(x,y);
    q = find(xstats(k+9).tstat.pval(2:end)<0.01);
    xstats2{k+9} = regstats(x,y(:,q));
    index = [m m+9];
    xstats2{k+9}.q = index(q);
    disp(xstats2{k+9}.rsquare);
end

s6 = s5;
sla6 = sla5;
Y = [s5 sla5];
for k = 4
    yr = Y(:,xstats2{k}.q);
    qnan =find(isnan(s5(:,k)) & all(~isnan(yr).').');
    xr = yr(qnan,:)*xstats2{k}.beta(2:end);
    s6(qnan,k) = xr; 
    qnan0 = find(isnan(s5(:,k)));
    disp(['filled ' num2str(length(qnan)) ' slope points out of ' num2str(length(qnan0)) ' at mooring ' num2str(k)]);
end
for k = 4
    yr = Y(:,xstats2{k+9}.q);
    qnan =find(isnan(sla5(:,k)) & all(~isnan(yr).').');
    xr = yr(qnan,:)*xstats2{k+9}.beta(2:end);
    sla6(qnan,k) = xr; 
    qnan0 = find(isnan(sla5(:,k)));
    disp(['filled ' num2str(length(qnan)) ' sla points out of ' num2str(length(qnan0)) ' at mooring ' num2str(k)]);
end

% #3 from offshore
clear xstats2 xstats
m = [4:9];
for k = 3;
    qall = find(all(~isnan([s2(:,[k m]) sla2(:,m)]).'));% times when all data are available %can't use sla at same location
    y = [s2(qall,m) sla2(qall,m)];
    x = s2(qall,k);
    xstats(k) = regstats(x,y);
    q = find(xstats(k).tstat.pval(2:end)<0.01);
    xstats2{k} = regstats(x,y(:,q));
    index = [m m+9];
    xstats2{k}.q = index(q);
    disp(xstats2{k}.rsquare);
end
for k = 3;
    qall = find(all(~isnan([s2(:,[m]) sla2(:,[k m])]).'));% times when all data are available %can't use sla slope at same location
    y = [s2(qall,m) sla2(qall,m)];
    x = sla2(qall,k);
    xstats(k+9) = regstats(x,y);
    q = find(xstats(k+9).tstat.pval(2:end)<0.01);
    xstats2{k+9} = regstats(x,y(:,q));
    index = [m m+9];
    xstats2{k+9}.q = index(q);
    disp(xstats2{k+9}.rsquare);
end

s7 = s6;
sla7 = sla6;
Y = [s6 sla6];
for k = 3
    yr = Y(:,xstats2{k}.q);
    qnan =find(isnan(s6(:,k)) & all(~isnan(yr).').');
    xr = yr(qnan,:)*xstats2{k}.beta(2:end);
    s7(qnan,k) = xr; 
    qnan0 = find(isnan(s6(:,k)));
    disp(['filled ' num2str(length(qnan)) ' slope points out of ' num2str(length(qnan0)) ' at mooring ' num2str(k)]);
end
for k = 3
    yr = Y(:,xstats2{k+9}.q);
    qnan =find(isnan(sla6(:,k)) & all(~isnan(yr).').');
    xr = yr(qnan,:)*xstats2{k+9}.beta(2:end);
    sla7(qnan,k) = xr; 
    qnan0 = find(isnan(sla6(:,k)));
    disp(['filled ' num2str(length(qnan)) ' sla points out of ' num2str(length(qnan0)) ' at mooring ' num2str(k)]);
end

% #2 from offshore
clear xstats2 xstats
m = [3:9];
for k = 2;
    qall = find(all(~isnan([s2(:,[k m]) sla2(:,m)]).'));% times when all data are available %can't use sla at same location
    y = [s2(qall,m) sla2(qall,m)];
    x = s2(qall,k);
    xstats(k) = regstats(x,y);
    q = find(xstats(k).tstat.pval(2:end)<0.01);
    xstats2{k} = regstats(x,y(:,q));
    index = [m m+9];
    xstats2{k}.q = index(q);
    disp(xstats2{k}.rsquare);
end
for k = 2;
    qall = find(all(~isnan([s2(:,[m]) sla2(:,[k m])]).'));% times when all data are available %can't use sla slope at same location
    y = [s2(qall,m) sla2(qall,m)];
    x = sla2(qall,k);
    xstats(k+9) = regstats(x,y);
    q = find(xstats(k+9).tstat.pval(2:end)<0.01);
    xstats2{k+9} = regstats(x,y(:,q));
    index = [m m+9];
    xstats2{k+9}.q = index(q);
    disp(xstats2{k+9}.rsquare);
end

s8 = s7;
sla8 = sla7;
Y = [s7 sla7];
for k = 2
    yr = Y(:,xstats2{k}.q);
    qnan =find(isnan(s7(:,k)) & all(~isnan(yr).').');
    xr = yr(qnan,:)*xstats2{k}.beta(2:end);
    s8(qnan,k) = xr; 
    qnan0 = find(isnan(s7(:,k)));
    disp(['filled ' num2str(length(qnan)) ' slope points out of ' num2str(length(qnan0)) ' at mooring ' num2str(k)]);
end
for k = 2
    yr = Y(:,xstats2{k+9}.q);
    qnan =find(isnan(sla7(:,k)) & all(~isnan(yr).').');
    xr = yr(qnan,:)*xstats2{k+9}.beta(2:end);
    sla8(qnan,k) = xr; 
    qnan0 = find(isnan(sla7(:,k)));
    disp(['filled ' num2str(length(qnan)) ' sla points out of ' num2str(length(qnan0)) ' at mooring ' num2str(k)]);
end

% #1 from offshore
clear xstats2 xstats
m = [2:9];
for k = 1;
    qall = find(all(~isnan([s2(:,[k m]) sla2(:,m)]).'));% times when all data are available %can't use sla at same location
    y = [s2(qall,m) sla2(qall,m)];
    x = s2(qall,k);
    xstats(k) = regstats(x,y);
    q = find(xstats(k).tstat.pval(2:end)<0.01);
    xstats2{k} = regstats(x,y(:,q));
    index = [m m+9];
    xstats2{k}.q = index(q);
    disp(xstats2{k}.rsquare);
end
for k = 1;
    qall = find(all(~isnan([s2(:,[m]) sla2(:,[k m])]).'));% times when all data are available %can't use sla slope at same location
    y = [s2(qall,m) sla2(qall,m)];
    x = sla2(qall,k);
    xstats(k+9) = regstats(x,y);
    q = find(xstats(k+9).tstat.pval(2:end)<0.01);
    xstats2{k+9} = regstats(x,y(:,q));
    index = [m m+9];
    xstats2{k+9}.q = index(q);
    disp(xstats2{k+9}.rsquare);
end

s9 = s8;
sla9 = sla8;
Y = [s8 sla8];
for k = 1
    yr = Y(:,xstats2{k}.q);
    qnan =find(isnan(s8(:,k)) & all(~isnan(yr).').');
    xr = yr(qnan,:)*xstats2{k}.beta(2:end);
    s9(qnan,k) = xr; 
    qnan0 = find(isnan(s8(:,k)));
    disp(['filled ' num2str(length(qnan)) ' slope points out of ' num2str(length(qnan0)) ' at mooring ' num2str(k)]);
end
for k = 1
    yr = Y(:,xstats2{k+9}.q);
    qnan =find(isnan(sla8(:,k)) & all(~isnan(yr).').');
    xr = yr(qnan,:)*xstats2{k+9}.beta(2:end);
    sla9(qnan,k) = xr; 
    qnan0 = find(isnan(sla8(:,k)));
    disp(['filled ' num2str(length(qnan)) ' sla points out of ' num2str(length(qnan0)) ' at mooring ' num2str(k)]);
end

% final data

sf = s9;
slaf = sla9;

%%
% convert slope to transport per unit distance Tx1
% define regression models

% convert slope into barotropic transport per unit distance
f = 2*pi*phi2f([latM(2:8)' 0.5*(sum(latM(9:10)))  0.5*(sum(latM(10:11))) ])/(60*60*24);
dp = [320 1260 2210 3607 3701 3994 4276 4389 4645];
fac = dp.*9.81.*(1./f)/10^3;% because gradient is in m per km
       
abr{1,1}.meanx  = -241.154580248103;
abr{2,1}.meanx = -1137.04488260785;
abr{3,1}.meanx = -1929.30641823704;
abr{4,1}.meanx = -2371.91240578889;
abr{5,1}.meanx = -1536.74883245249;
abr{6,1}.meanx = -787.11394987845;
abr{7,1}.meanx = -361.99547122813;
abr{8,1}.meanx =  -6.89062783586677;
abr{9,1}.meanx =  200.107616655137;

abr{1,1}.beta = [-244.948916379702         0.784239077663288];
abr{2,1}.beta = [-772.394101118125         0.456482023519002];
abr{3,1}.beta = [ -829.542414061001         0.368644651864055];
abr{4,1}.beta = [-635.527162874393         0.351781786384768];
abr{5,1}.beta = [-376.578162776363         0.302565551971023];
abr{6,1}.beta = [  -215.801902477374         0.195831349823038];
abr{7,1}.beta = [ -122.485214791723         0.305821715011084];
abr{8,1}.beta = [44.6918498033721         0.327727412769925];
abr{9,1}.beta = [ 93.0929217731398         0.313940401385852];

abr{1,2}.beta = [-255.678646239298         0.698312720164358];
abr{2,2}.beta = [-809.418188453189         0.402015721517372];
abr{3,2}.beta = [-892.360326127802         0.277934891548132];
abr{4,2}.beta = [-728.314794840215         0.265946846932697];
abr{5,2}.beta = [-427.348765133512          0.27600380739157];
abr{6,2}.beta = [-277.780231952403         0.176757914787721];
abr{7,2}.beta = [-221.533262065019         0.241884813261591];
abr{8,2}.beta = [-139.252371099845         0.194994052520951];
abr{9,2}.beta = [-101.613264034412         0.181361019362705];


clear foo;
for k = 1:9;
    foo(k) = abr{k,1}.meanx;
end

dum = bsxfun(@plus,bsxfun(@times,fac,sf),-foo);

Tx1 = NaN*dum;
clear foo;
for k = 1:9;
    foo(k,1) = abr{k,1}.beta(2);
    foo(k,2) = abr{k,1}.beta(1);
end

for k = 1:size(dum,1)
    Tx1(k,:) = bsxfun(@plus,bsxfun(@times,dum(k,:),foo(:,1)'),foo(:,2)');
end

clear dum* foo*

% convert slope to SOUTHWARD transport per unit distance Tx1sw
% need to remove the mean from the relationship with the NET southward transport
clear foo;
for k = 1:9;
    foo(k) = abr{k,1}.meanx;
end
dum = bsxfun(@plus,bsxfun(@times,fac,sf),-foo);

Tx1sw = NaN*dum;
clear foo;
for k = 1:9;
    foo(k,1) = abr{k,2}.beta(2);
    foo(k,2) = abr{k,2}.beta(1);
end

for k = 1:size(dum,1)
    Tx1sw(k,:) = bsxfun(@plus,bsxfun(@times,dum(k,:),foo(:,1)'),foo(:,2)');
end

clear dum* foo*

%%
x = dM3'- dM2(1);
xscale = 0:1:300;

Tx = NaN*ones(size(Tx1,1),length(xscale));
Txsw = Tx;

dx = 1;
L = 56;% 
wa = (1/L)/(1/(2*dx));
[b10k,a10k] = butter(3,wa,'low');

for m = 1:size(Tx1,1)
    qan = find(~isnan(Tx1(m,:)));
    if ~isempty(qan)
        Ym = [0 Tx1(m,qan)];
        Xm = [ 0 x(qan)'];
        % piecewise cubic hermite interpolating spline; shape preserving
        foo3 = pchiptx(Xm,Ym,xscale);
        Tx(m,:) = foo3(:)';
    end
end

qpos = find(Tx1sw>0);
Tx1sw(qpos) = 0;

for m = 1:size(Tx1sw,1)
    qan = find(~isnan(Tx1sw(m,:)));
    if ~isempty(qan)
        Ym = [0 Tx1sw(m,qan)];
        Xm = [ 0 x(qan)'];
        % piecewise cubic hermite interpolating spline
        foo3 = pchiptx(Xm,Ym,xscale);
        Txsw(m,:) = foo3(:)';
    end
end

% now look for first maximum past 100 km
dTxdx = diff(Tx,1,2);
d2Txdx2 = diff(dTxdx,1,2);

x0 = cell(length(talt),2);
for k = 1:length(talt)
        foo = findzero(xscale(1:end-1)+0.5,squeeze(dTxdx(k,:)));% transport extrema
        if ~isempty(foo)
            x0{k,1} = foo;
            foo2 = interp1(xscale(2:end-1),squeeze(d2Txdx2(k,:)),foo);% figure out minima and maxima
            x0{k,2} = foo2;
        end
        clear foo*
end

X =  NaN*ones(length(talt),1);
% look for maximum further away than mid distance of E and D
for k = 1:length(talt)
        foo1 = x0{k,1};
        foo2 = x0{k,2};
        if ~isempty(foo1)
            dum = find(foo2<0&foo1>104);
            if ~isempty(dum)
                X(k,1) = foo1(min(dum));
            else
                X(k,1) = 300;
            end
        else
            X(k,1) = 300;
        end
 end


% calculate transport for each profile at each time step
% zero the positive values of southwest transport
qpos = find(Txsw>0);
foo = Txsw;
foo(qpos) = 0;
TTsw = NaN*ones(length(talt),1);
for k = 1:length(talt)
        q = find(abs(X(k) - xscale)<0.5);
        TTsw(k) = sum(foo(k,1:q))*1000;
end

TTbox = sum(Tx(:,1:220),2)*1000;

% correct for the means; in situ mean from ACT experiment for Tjet is -83726616
% mean for Tbox is -77488032
Toffset1 =  -83726616 - mean(TTsw(talt>=734245.5 & talt<=735284));
TTsw = TTsw + Toffset1;
Toffset2 =  -77488032 - mean(TTbox(talt>=734245.5 & talt<=735284));
TTbox = TTbox + Toffset2;

if 1
    figure
    subplot(2,1,1)
    hold on
    plot(talt,TTbox/10^6);
    title('Agulhas boundary transport');
    ylabel('Sv');
    datetick('x','yy');
    subplot(2,1,2)
    hold on
    plot(talt,TTsw/10^6);
    title('Agulhas jet transport');
    ylabel('Sv');
    datetick('x','yy');
end

