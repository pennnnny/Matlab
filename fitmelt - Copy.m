
function Coverage = fitmelt()
%% get data

 %% choose directory that want to grab .txt from and change to that directory
        curPath = pwd;
        folderName = uigetdir(curPath,'Choose folder with cvdata in it.');
        cd(folderName);
    
         % get all the .txt files
        files = dir('*.txt');
        cvFiles = [];
        for x = 1:length(files)
            name = lower(files(x).name);
            if ~isempty(strfind(name, '.txt')) || ...
                    ~isempty(strfind(name, '.txt')) || ...
                    ~isempty(strfind(name, '.txt'))
                cvFiles = [cvFiles files(x)];  %#ok<*AGROW>
            end
        end
        %% read all data
        cvData = [];
        for i = 1:length(cvFiles)
            cvID = fopen(cvFiles(i).name);
            cvDat = textscan(cvID, '%f%f', ' ');
            cvData = [cvData cvDat];
        end
      

%% get Initial parameter.
intparDir = uigetdir('Choose int parameter folder.');
intparName = uigetfile('Choose int parameter file.');
intparPath = [intparDir, '\', intparName];

%% Read columns of data as strings:
% For more information, see the TEXTSCAN documentation.
formatSpec = '%7s%2s%1s%26s%7s%7s%s%[^\n\r]';
fileID = fopen(intparPath,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', '', 'WhiteSpace', '',  'ReturnOnError', false);
fclose(fileID);

%% Convert the contents of columns containing numeric strings to numbers.
% Replace non-numeric strings with NaN.
    raw = repmat({''},length(dataArray{1}),length(dataArray)-1);
    for col=1:length(dataArray)-1
        raw(1:length(dataArray{col}),col) = dataArray{col};
    end
    numericData = NaN(size(dataArray{1},1),size(dataArray,2));

    for col=[4,5,6,7]
        % Converts strings in the input cell array to numbers. Replaced non-numeric
        % strings with NaN.
        rawData = dataArray{col};
        for row=1:size(rawData, 1);
            % Create a regular expression to detect and remove non-numeric prefixes and
            % suffixes.
            regexstr = '(?<prefix>.*?)(?<numbers>([-]*(\d+[\,]*)+[\.]{0,1}\d*[eEdD]{0,1}[-+]*\d*[i]{0,1})|([-]*(\d+[\,]*)*[\.]{1,1}\d+[eEdD]{0,1}[-+]*\d*[i]{0,1}))(?<suffix>.*)';
            try
                result = regexp(rawData{row}, regexstr, 'names');
                numbers = result.numbers;

                % Detected commas in non-thousand locations.
                invalidThousandsSeparator = false;
                if any(numbers==',');
                    thousandsRegExp = '^\d+?(\,\d{3})*\.{0,1}\d*$';
                    if isempty(regexp(thousandsRegExp, ',', 'once'));
                        numbers = NaN;
                        invalidThousandsSeparator = true;
                    end
                end
                % Convert numeric strings to numbers.
                if ~invalidThousandsSeparator;
                    numbers = textscan(strrep(numbers, ',', ''), '%f');
                    numericData(row, col) = numbers{1};
                    raw{row, col} = numbers{1};
                end
            catch me
            end
        end
    end

%% Split data into numeric and cell columns.
rawNumericColumns = raw(:, [4,5,6,7]);
rawCellColumns = raw(:, [1,2,3]);


%% Replace non-numeric cells with NaN
R = cellfun(@(x) ~isnumeric(x) && ~islogical(x),rawNumericColumns); % Find non-numeric cells
rawNumericColumns(R) = {NaN}; % Replace non-numeric cells

%% Allocate imported array to column variable names
interce = rawCellColumns(:, 1);
pt = rawCellColumns(:, 2);
f = rawCellColumns(:, 3);
para1 = cell2mat(rawNumericColumns(:, 1)); %initial value
para2 = cell2mat(rawNumericColumns(:, 2));
para3 = cell2mat(rawNumericColumns(:, 3));
para4 = cell2mat(rawNumericColumns(:, 4));

%% for forward CV scan 
para = para1(1:8);

%% Clear temporary variables
clearvars filename delimiter formatSpec fileID dataArray ans;

%% use cvcurveFit to find area
cvlength = 2*length(files);
fitresultA1 = [];
fitresultA2 = [];
fitresultA3 = [];
fitresultA4 = [];
fitresultA5 = [];
fitresultA6 = [];
fitresultA7 = [];
fitresultA7f = [];
fitresultA8 = [];
fitR = [];
for x = 1:2:cvlength
        datax = cvData{x};
        dataxf = datax(200:700);
        datay = cvData{x+1};
        datayf = datay(200:700);
        if x == 1
            paraf = para;
        else
            paraf = fitR;
        end
        fitresult = cvcurveFit(dataxf, datayf, para);
        [pk,Vindex] = getPeak(datayf,dataxf);
        peak = dataxf(Vindex);
        % creat parameter array
        fitresultA1 = [fitresultA1, fitresult.a1];
        fitresultA2 = [fitresultA2, fitresult.a2];
        fitresultA3 = [fitresultA3, fitresult.a3];
        fitresultA4 = [fitresultA4, fitresult.a4];
        fitresultA5 = [fitresultA5, fitresult.a5];
        fitresultA6 = [fitresultA6, fitresult.a6];
        %peak value
        fitresultA7 = [fitresultA7, peak];
        fitresultA7f = [fitresultA7f, fitresult.a7];
        fitresultA8 = [fitresultA8, fitresult.a8];    
        fitR = [fitresult.a1,fitresult.a2,fitresult.a3,fitresult.a4,fitresult.a5,fitresult.a6,fitresult.a7,fitresult.a8];
        %generate baseline
        cb = @(x) fitresult.a1+fitresult.a2*x+fitresult.a3*(exp((x.^fitresult.a4)/fitresult.a5)-1);
        dataytb = [];
        for x = 0.2:0.001:0.7
            datab = cb(x);
            dataytb = [dataytb, datab];
        end
        hold on 
        plot(dataxf(Vindex), pk, 'or')
        plot(dataxf, dataytb, 'k')
        hold off
end
paraResult = [fitresultA1; fitresultA2; fitresultA3;fitresultA4;fitresultA5;fitresultA6;fitresultA7f;fitresultA8];
paraSize =  size(paraResult);
paralength = paraSize(2);
peakArea = [];
    for x = 1:paralength
        peakarea = findArea(0.2,0.7,paraResult(:,x));
        peakArea = [peakArea peakarea];
    end
%% Convert to coverage
Coverage = [];    
for x = 1:length(peakArea)
    Cover = getcover(peakArea(x), 1.43);
    Coverage = [Coverage, Cover];
end
end 
%% CVfit function
function [fitresult, gof] = cvcurveFit(datax, datay, para)

    [xData, yData] = prepareCurveData( datax, datay);

    % Set up fittype and options.
    ft = fittype( 'a1+a2*x+a3*(exp((x^a4)/a5)-1)+(a6*exp(a8*(x-a7))/(1+exp(a8*(x-a7)))^2)', 'independent', 'x', 'dependent', 'y' );
    opts = fitoptions( ft );
    %opts.Algorithm = 'Levenberg-Marquardt';
    opts.Display = 'Off';
    opts.Lower = [-Inf -Inf -Inf -Inf -Inf -Inf -Inf -Inf];
    opts.MaxIter = 2000;
    opts.Robust = 'LAR';
    opts.StartPoint = para;
    opts.Upper = [Inf Inf Inf Inf Inf Inf Inf Inf];

    % Fit model to data.
    [fitresult, gof] = fit( xData, yData, ft, opts );

    % Plot fit with data.
    figure( 'Name', 'cvfit 1' );
    h = plot( fitresult, xData, yData );
    legend( h, 'datay vs. datax', 'cvfit 1', 'Location', 'NorthEast' );
    % Label axes
    xlabel( 'datax' );
    ylabel( 'datay' );
    grid on
end

%% function to find peak area
function area = findArea(bi,bf,para)
    a1 = para(1);
    a2 = para(2);
    a3 = para(3);
    a4 = para(4);
    a5 = para(5);
    a6 = para(6);
    a7 = para(7);
    a8 = para(8);
    cp = @(x) a1+a2*x+a3*(exp((x.^a4)/a5)-1)+(a6*exp(a8*(x-a7))/(1+exp(a8*(x-a7))).^2);
    dataxt = bi:0.001:bf;
    datayt = [];
    for x = bi:0.001:bf
        datat = cp(x);
        datayt = [datayt, datat];
    end
    Ap = trapz(dataxt, datayt);
    cb = @(x) a1+a2*x+a3*(exp((x.^a4)/a5)-1);
    dataytb = [];
    for x = bi:0.001:bf
        datab = cb(x);
        dataytb = [dataytb, datab];
    end
    Ab = trapz(dataxt, dataytb);
    area = Ap-Ab;
end
%% get Peak position
function [pk,Vindex] = getPeak(x, y)
    [pk,Vindex] = findpeaks(x,'NPeaks',1,'SortStr','descend');
    peaky = y(Vindex);
    %plot
    %hold on 
    %plot(y,x)
    %plot(y(Vindex), pk, 'or')
    %hold off
end
%% Get DNA coverage
function coverage = getcover(area,rough)
f =@(area) (area/(20*96485.3*0.0201*rough))*6.022E+23;
coverage = f(area);
end

    