% Finds all the picture files (jpeg and png) in a folder and calculate the
% angles the water drop makes with the surface. Can either a single value,
% a vector of values (that is the same length as the number of pictures
% files in the folder), or no values in which you select threshold values
% from pictures. 
function angles = findAngles()
    
    % turn off all warnings
    warning('off','all');
    % display warning
    message = sprintf(['HOWDY! How are you doing this fine day? Great!\n'...
                        'Before we get started I wanted to make sure of a few things.\n'...
                        '1. All the pictures files are .jpg or .png\n'...
                        '2. There arne''t any random .jpg or .png files in that folder'...
                        'If not I''m ready to get started! Click Yes if you are too!']);
    if ~strcmp(questdlg(message),'Yes')
        msgbox('Oh, well that''s unfortunate');
    else
        %% choose directory that want to grab .jpegs from and change to that directory
        curPath = pwd;
        folderName = uigetdir(curPath,'Choose folder with pictures in it.');
        cd(folderName);
    
        %% get all the .jpeg and png files
        files = dir();
        imageFiles = [];
        for x = 1:length(files)
            name = lower(files(x).name);
            if ~isempty(strfind(name, '.jpeg')) || ...
                    ~isempty(strfind(name, '.jpg')) || ...
                    ~isempty(strfind(name, '.png'))
                imageFiles = [imageFiles files(x)];  %#ok<*AGROW>
            end
        end
        
        %% choose what type of approximation, circle or linear
        type = questdlg('Choose what type of approximation',...
            'Approximation Type',...
            'linear', 'circle','linear');
        %% deal with threshold vector
%         if nargin == 0
%             threshold = zeros(1,length(imageFiles));
%             for x = 1:length(imageFiles)
%                 fig = figure; i = imshow(imageFiles(x).name);
%                 pos = impoint().getPosition;
%                 threshold(x) = i.CData(round(pos(2)),round(pos(1)));
%                 close(fig);
%             end
%         elseif length(varargin{1}) == 1
%             threshold = ones(1,length(imageFiles)) * varargin{1};
%         else
%             if length(varargin{1}) == length(imageFiles)
%                 threshold = varargin{1};
%             else
%                 error('Length of input threshold doens''t match the length of the number of Images');
%             end
%         end
        
        %% use getAngle to calculate all the angles
        angles = zeros(1,length(imageFiles));
        h = figure(1);
        for y = 1:length(imageFiles)
            imArr = imread(imageFiles(y).name);
            arr = imArr(:,:,1);
            [prof, threshold] = findProfile(arr);
            [angle, outVals] = getAngle(arr, prof, type);
            angles(y) = angle;
            [prof, changeLoc, xvals, upvals, downvals] = deal(outVals{1},outVals{2},outVals{3},outVals{4},outVals{5});
            %subplot(3,length(imageFiles),y);
            %imshow(arr)
            %if y == 1
            %    ylabel('Image');
            %end
            %title(imageFiles(y).name);
            subplot(2,length(imageFiles),y);%+length(imageFiles)
            imshow(arr)
            hold on
            for z = 1:length(prof)
                plot(prof(z,2),prof(z,1),'.b');
            end
            if y == 1
                ylabel('Profile')
            end
            xlabel(sprintf('Threshold: %d',threshold));
            subplot(2,length(imageFiles),y+length(imageFiles));%+2*length(imageFiles)
            imshow(arr);
            hold on;
            if strcmp(type,'circle') && min(xvals) < prof(changeLoc,2)
                plot(xvals,upvals,'b',xvals+200,downvals,'r','LineWidth',1);
            else
                plot(xvals,upvals,'b',xvals,downvals,'r','LineWidth',1);
            end
            plot(prof(changeLoc,2),prof(changeLoc,1),'y*');
            if y == 1
                ylabel('Angle Diagram');
            end
            xlabel(sprintf('%.2f%c',angle, char(176)));
            set(h,'position',get(0,'screensize'));
        end
        
        % change back to previous directory
        cd(curPath);
    end
end

% Takes in image file and a threshold for that image and gives you an angle
% and a plot showing you what it looks like
function [angle, outVals] = getAngle(arr, prof, type)
    
    %% find changing location
    yloc = round(length(prof)/2);
    found = false;
    direction = 1;
    if prof(yloc,2) >= prof(yloc+10,2) % left
        while ~found
            if prof(yloc,2) < prof(yloc+1,2)
                count = 0;
                for n = yloc+1:yloc+20
                    if prof(yloc,2) < prof(n,2)
                        count = count + 1;
                    else
                        break;
                    end
                end
                if count == 20
                    found = true;
                end
            end
            yloc = yloc + 1;
        end
    else % right
        while ~found
            direction = -1;
            if prof(yloc,2) > prof(yloc+1,2)
                count = 0;
                for n = yloc+1:yloc+20
                    if prof(yloc,2) > prof(n,2)
                        count = count + 1;
                    else
                        break;
                    end
                end
                if count == 20
                    found = true;
                end
            end
            yloc = yloc + 1;
        end
    end
    changeLoc = yloc - 1;

    %% determine angle using helper functions
    [r,~] = size(arr);
    profR = [r-prof(:,1), prof(:,2)];
    switch type
        case 'linear'
            % this is for linear approximation averaging 50 points
            [angle,upDis,downDis] = linearCalcAngle(profR,changeLoc,direction);
        case 'circle'
            % this is for tangent line circle approximation
            [angle,upDis,downDis] = circleCalcAngle(profR,changeLoc,direction);
        otherwise
            error('Please choose either "linear" or "circle"');
    end
    
    %% prepare output
    xvals = linspace(prof(changeLoc,2), prof(changeLoc,2)+200*direction,200);
    upvals = prof(changeLoc,1):upDis:prof(changeLoc,1)+199*upDis;
    downvals = prof(changeLoc,1):downDis:prof(changeLoc,1)+199*downDis;
    if isempty(downvals)
        downvals = zeros(1,length(xvals)) + prof(changeLoc,1);
    end
    outVals = {prof, changeLoc, xvals, upvals, downvals};
    
    %% plot image
    
    % imshow(imArr);
    % hold on;
    % plot(xvals,upvals,xvals,downvals,'LineWidth',1);
    % plot(prof(changeLoc,2),prof(changeLoc,1),'*')
    
end

function [prof, threshold] = findProfile(layer)
    % load image and start at top
    [r,c] = size(layer);
    p = [1, c];
    
    % get threshold
    y = round(r/2);
    avl = mean(layer(y,1:50));
    avr = mean(layer(y,end-50:end));
    threshold = round(mean([avl avr]));
    
    % find top most point
    while layer(p(1),p(2)) > threshold
        p(1) = p(1) + 1;
    end
    
    while layer(p(1),p(2)) < threshold
        p(2) = p(2) - 1;
    end
    
    % constantly only look 5 pixels previous for the profile
    prof = p;
    intact = true;
    box = 10;
    while intact
        newp = [p(1) + 1, p(2) + 50];
        while intact && layer(newp(1),newp(2)) < threshold && newp(2) ~= 1
            newp(2) = newp(2) - 1;
        end
        if newp(2) > p(2) + box || newp(2) < p(2) - box
               intact = false; 
        end
        if intact
            prof = [prof; newp]; %#ok<*AGROW>
            p = newp;
        end
    end
    
    % cut off bottom tip cuz it's usually pretty whack
     prof(end-40:end,:) = [];
     
    xvals = prof(:,1);
    yvals = prof(:,2);
    newy = spline(xvals, yvals, linspace(xvals(1),xvals(end),length(xvals)));
%     for z = 1:length(prof)
%         plot(newy(z),prof(z,1),'.r');
%     end
%     title('Splined points');
    prof(:,2) = newy;
end

function [angle,upDis,downDis] = linearCalcAngle(prof,changeLoc,direction)
    if length(prof) - changeLoc < 60
        delta = length(prof)-changeLoc;
    else
        delta = 59;
    end
    upCurve = prof(changeLoc-delta:changeLoc,:);
    downCurve = prof(changeLoc:changeLoc+delta,:);
    upCoef = polyfit(upCurve(:,2),upCurve(:,1),1);
    downCoef = polyfit(downCurve(:,2),downCurve(:,1),1);
    upDis = -abs(upCoef(1)); downDis = abs(downCoef(1));
    angle1 = acosd(dot([direction upDis],[1 0])/(norm([direction upDis])*norm([1 0])));
    angle2 = acosd(dot([direction downDis],[1 0])/(norm([direction downDis])*norm([1 0])));
    angle = mean([angle1,angle2]);
end

function [angle,upDis,downDis] = circleCalcAngle(prof,changeLoc,direction)
    % first get 3 points including changeLoc, maybe experiment with good
    % choices
    point1 = [prof(changeLoc, 2),prof(changeLoc,1)];
    point2 = [prof(changeLoc-100, 2),prof(changeLoc-100,1)];
    point3 = [prof(changeLoc-200, 2),prof(changeLoc-200,1)];
    [ucenter,uradius] = calc_circle(point1,point2,point3);
    uvec = (point1-ucenter)/uradius;
    upDis = [direction*uvec(2),-1*direction*uvec(1)];
%     point5 = [prof(changeLoc+20, 2),prof(changeLoc+20,1)];
%     point6 = [prof(changeLoc+40, 2),prof(changeLoc+40,1)];
%     [dcenter,dradius] = calc_circle(point1,point5,point6);
%     dvec = (point1-dcenter)/dradius;
%     downDis = [direction*dvec(2),-1*direction*dvec(1)];
    upDis = -abs(upDis(2)/upDis(1));
    angle = acosd(dot([direction upDis],[1 0])/(norm([direction upDis])*norm([1 0])));
    downDis = 0;
%     angle2 = acosd(dot(downDis,[1 0])/(norm(downDis)*norm([1 0])));
%     angle = mean([angle1,angle2]);
end

function [centre, radius] = calc_circle(pt1, pt2, pt3)

delta_a = pt2 - pt1;
delta_b = pt3 - pt2;

ax_is_0 = abs(delta_a(1)) <= 0.000000001;
bx_is_0 = abs(delta_b(1)) <= 0.000000001;

% check whether both lines are vertical - collinear
if (ax_is_0 && bx_is_0)
    centre = [0 0];
    radius = -1;
    return
end

% make sure delta gradients are not vertical
% re-arrange points to change deltas
if (ax_is_0)
    [centre, radius] = calc_circle(pt1, pt3, pt2);
    return
end
if (bx_is_0)
    [centre, radius] = calc_circle(pt2, pt1, pt3);
    return
end 

grad_a = delta_a(2) / delta_a(1); 
grad_b = delta_b(2) / delta_b(1); 

% check whether the given points are collinear
if (abs(grad_a-grad_b) <= 0.000000001)
    centre = [0 0]; 
    radius = -1;
    return
end

% swap grads and points if grad_a is 0
if abs(grad_a) <= 0.000000001
    tmp = grad_a;
    grad_a = grad_b;
    grad_b = tmp;
    tmp = pt1;
    pt1 = pt3;
    pt3 = tmp;
end 

% calculate centre - where the lines perpendicular to thecentre of
% segments a and b intersect.
centre(1) = ( grad_a*grad_b*(pt1(2)-pt3(2)) +...
grad_b*(pt1(1)+pt2(1)) - grad_a*(pt2(1)+pt3(1)) ) /...
(2*(grad_b-grad_a));
centre(2) = ((pt1(1)+pt2(1))/2 - centre(1)) / grad_a +...
(pt1(2)+pt2(2))/2;

% calculate radius
radius = norm(centre - pt1);
end