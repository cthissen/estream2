function [XY] = estream2(varargin)
%ESTREAM2  2D evenly spaced streamlines.
%   XY = ESTREAM2(X,Y,U,V,STARTX,STARTY) computes streamlines from vector
%   data U,V. The arrays X,Y define the coordinates for U,V and must be
%   monotonic and 2D plaid (as if produced by MESHGRID).  STARTX and
%   STARTY define the starting positions of the streamlines. A cell array
%   of vertex arrays is returned in XY.
%   
%   XY = ESTREAM2(U,V,STARTX,STARTY) assumes [X Y] = meshgrid(1:N, 1:M)
%        where[M,N]=SIZE(U).  
%
%   XY = ESTREAM2(...,OPTIONS) specifies the options used in creating the
%   streamlines. OPTIONS is specified as a one, two, or three element vector
%   containing the step size, the maximum number of points in a stream line, and
%   the separation distance between streamlines in grid units.  If OPTIONS is
%   not specified the default step size is 1/10th the minimum cell distance, the
%   default maximum number of points is 10000, and the default separation
%   distance is 5*stepsize. OPTIONS can either be [stepsize], [stepsize
%   maxpoints], or [stepsize maxpoints sepdistance].
% 
%   XY = ESTREAM2(...,OPTIONS,'plot') plots each complete sequence onto figure(1)
%   as soon as it is complete. This can save time, because it allows you to
%   evaluate your choice of sepdistance without calculating all the streamlines.
%   For example, if you make sepdistance too small, the program may require
%   several minutes to complete.
%   
%   Example:
%     load wind
%     startX = 80;
%     startY = 35;
%     stepSize = 0.01;
%     maxIter = 10000;
%     dSep = 2;
%     XY = estream2(x(:,:,5),y(:,:,5),u(:,:,5),v(:,:,5),startX,startY,[stepSize,maxIter,dSep],'plot');
%     streamline(XY);
%   See also STREAM2, STREAMLINE, STREAM3, CONEPLOT, ISOSURFACE, SMOOTH3,
%   SUBVOLUME,
%            REDUCEVOLUME.
% 
% Method overview:
% Streamlines are calculated using a modified version of the
% method presented in B. Jobard and W. Lefer, Creating Evenly-Spaced Streamlines
% of Arbitrary Density, Proc. Eighth Eurographics Workshop on Visualization in
% Scientific Computing, pp. 45-55, 1997.
%
% The first streamline is calculated from the point (STARTX,STARTY). New seed
% points are generated a distance dSep away from each streamline. Streamlines
% are stopped if they come within 0.5*dSep of extant streamlines.
%
% No toolboxes are required. Required subfunctions are included (dispstat.m and
% interparc.m).
%
% Things to note: This code uses a 4th order Runge-Kutta method to integrate the
% streamlines. The stepsize is fixed.
% 
%%% Input variables
% U and V are the x and y direction cosines of the 2d lineation
% X and Y give the location of xDir and yDir.
% STARTX and STARTY see the first streamline.
% stepSize gives the step size for the Runge Kutta method
% iterMax is the maximum number of iterations allowed for any one trajectory
% dSep is the desired separation distance between streamlnes.
% 
%%% Output variables
% output is stored as a cell array that is equivalent to matlab's stream2 output
% 
% Christopher Thissen, Yale University.
% 
% V 1.0 Dec 2015

%% parse and check inputs
narginchk(4,8);

plotCase = 'noplot'; % default is do not plot
xGrid = [];
options = [];

% check if final input is a string
if ischar(varargin{end})
    nin = nargin-1;
    plotCase = varargin{end};
else
    nin = nargin;
end

% determine if grid has been input
if nin==4 || nin==5       % stream2(u,v,sx,sy)
    xDir   = varargin{1};
    yDir   = varargin{2};
    startX = varargin{3}; % these don't make sense if 
    startY = varargin{4};
    if nin==5, 
        options = varargin{5}; 
    end
  
elseif nin==6 || nin==7 || nin==8   % stream2(x,y,u,v,sx,sy)
    xGrid  = varargin{1};
    yGrid  = varargin{2};
    xDir   = varargin{3};
    yDir   = varargin{4};
    startX = varargin{5};
    startY = varargin{6};
    if nin==7, 
        options = varargin{7}; 
    end
else
  error(message('MATLAB:estream2:WrongNumberOfInputs')); 
end

startX = startX(:); 
startY = startY(:); 



%%% check inputs
% [msg, x, y] = xyuvcheck(x,y,u,v);  
% error(msg);
if isempty(xGrid)
    [nX,nY] = size(xDir);    
    xGrid = linspace(1,nX,nX);
    yGrid = linspace(1,nY,nY);
    [xGrid,yGrid] = meshgrid(yGrid,xGrid);
%     keyboard
end
[nX,nY] = size(xGrid);
xLim(1) = min(xGrid(:));
xLim(2) = max(xGrid(:));
yLim(1) = min(yGrid(:));
yLim(2) = max(yGrid(:));
cellSize = [(xLim(2)-xLim(1))/nX; (yLim(2)-yLim(1))/nY];

% set default values for options
stepsize = 0.1*min(cellSize); % default values first
 iterMax = 10000;
    dSep = 5*stepsize;
    
% parse optional inputs (same as stream2)
optargs = [stepsize iterMax dSep];
optargs(1:numel(options)) = options;
stepsize = optargs(1);
 iterMax = optargs(2);
    dSep = optargs(3);
         

validateattributes(startX,{'numeric'},  {'real','scalar'});
validateattributes(startY,{'numeric'},  {'real','scalar'});
validateattributes(stepsize,{'numeric'},{'real','positive','scalar'});
validateattributes(iterMax,{'numeric'}, {'real','positive','integer','scalar'});
validateattributes(dSep,{'numeric'},    {'real','positive','scalar'});
validateattributes(xLim,{'numeric'},    {'real','numel',2,'increasing'});
validateattributes(yLim,{'numeric'},    {'real','numel',2,'increasing'});
validatestring(plotCase,{'plot','noplot'});

checksamesize(xGrid,yGrid);
checksamesize(xDir,yDir);

% if plotting, setup figure
switch plotCase
    case 'plot'
        hFig = figure(1); 
        hAx = axes;
        hold(hAx,'on');
        hAx.XLim = xLim;
        hAx.YLim = yLim;
end
%% setup variables

% setup output variable
XY = mat2cell([startX(:),startY(:)],1);


%%% set internal variables
iStm = 1; % initialize streamline counts
tNow = 0; % initialize time
dt = stepsize; % set parametric step distance (grid units)
xNow = startX(:)';
yNow = startY(:)';
iterCt = 0; % initialize iteration counts 
idxCand = 0; % for indexing potential streamline starting points
nCand = 1; % initialize number of candidates
dTest = 0.5*dSep; % for testing if streamlines are too close (0.5 factor recommended by Jobard and Leifer)



nSyms = 20; % number of candidate starting points from each established streamline
t = linspace(0,1,nSyms); % interparc.m requires input as t       

xAll = []; % store all streamline points here
yAll = [];

xCand = []; % candidate points to start new streamline
yCand = [];


dispstat('Cleavage Trajectory Calculation\n','timestamp','init','keepthis');
while idxCand <= nCand
    dispstat(sprintf('Iter:%i Streamline Number: %i',...
        iterCt,iStm),'timestamp');
    iterCt = iterCt+1;
    
    %% Calculate next point in streamline
    
    % set even/odd flag. odd indicates forward calculation
    oddStm = mod(iStm,2)==1;                
    if oddStm
        % odd streamlines go forward
        [xNow,yNow,tNow] = RKstep(xGrid,yGrid, xDir, yDir,xNow,yNow,tNow,dt);
    else
        % even streamlines go backward
        [xNow,yNow,tNow] = RKstep(xGrid,yGrid,-xDir,-yDir,xNow,yNow,tNow,dt);
    end

    
    %% check if particle has exited domain

    out = checkdomain(xNow,yNow,xLim,yLim);

    % check if too close to other streamlines
    iDist = min(sqrt((xAll-xNow).^2 + (yAll-yNow).^2));
    outDist = iDist < dTest;    
    
    % check if point hasn't moved
    % (simple stagnation detection)
    xPrev = XY{iStm}(end,1);
    yPrev = XY{iStm}(end,2);
    moved = sqrt((xPrev-xNow)^2 + (yPrev-yNow)^2);
    mvOut = moved < 10*eps;
    
    % check if point is too close to any other points in the same streamline
    % (simplistic loop detection)
    nBack = 2*round(dTest/stepsize);
    iDistStm = min(sqrt((XY{iStm}(1:end-nBack,1)-xNow).^2 + (XY{iStm}(1:end-nBack,2)-yNow).^2));
    outDistStm = iDistStm < dTest;
%     outDistStm=[];
    
    out = [out,outDist,mvOut,outDistStm];        

    %% if point is "out", start new streamline
    if any(out) || iterCt > iterMax    
        % streamline is finished. 
      
        % add previous streamline points to total streamline list
        %... don't do this if odd. Otherwise reverse streamlines always fail
        if oddStm
        else
            xAll = [xAll(:); XY{iStm-1}(:,1); XY{iStm}(:,1)];
            yAll = [yAll(:); XY{iStm-1}(:,2); XY{iStm}(:,2)];
        end

        % generate new candidate points from the just-established streamline
        if numel(XY{iStm}(:,1))<2 
            % streamline has only 1 point, orthogonal points not defined
        else
            % generate evenly spaced points along streamline
            [pt,dudt,~] = interparc(t,XY{iStm}(:,1),XY{iStm}(:,2));
            
            % get points orthogonal to line (both sides)            
            for iP = 1:numel(pt(:,1));
                rotAng = atan2d(dudt(iP,2),dudt(iP,1))-90; % positive is CW
                
                xTmp(1) = pt(iP,1) + cosd(rotAng)*dSep;
                xTmp(2) = pt(iP,1) - cosd(rotAng)*dSep;
                
                yTmp(1)= pt(iP,2) + sind(rotAng)*dSep;
                yTmp(2)= pt(iP,2) - sind(rotAng)*dSep;
                
                xCand = [xCand(:); xTmp(:)];
                yCand = [yCand(:); yTmp(:)];
            end
            nCand = numel(xCand); % reset number of candidates
            
        end
                
        if oddStm   
            % restart from initial point in the reverse direction
            iterCt=0;
            xNow = XY{iStm}(1,1);
            yNow = XY{iStm}(1,2);

            iStm = iStm+1;                    
            XY{iStm} = [xNow,yNow];            
            
        else

            % find next viable starting location
            while idxCand < nCand
                idxCand = idxCand+1;            
                xCandNow = xCand(idxCand);
                yCandNow = yCand(idxCand);
            
                outCand = checkdomain(xCandNow,yCandNow,xLim,yLim);                
                iDist = min((sqrt((xAll-xCandNow).^2 + (yAll-yCandNow).^2)));                                        
                
               
                if any(outCand) || iDist < dTest
                    % outside domain or too close. try next candidate
                else
                    % inside domain and not too close to extant streamlines.
                    % start new streamline
                    iterCt=0;
                    iStm = iStm+1;                    
                    xNow = xCandNow;
                    yNow = yCandNow;
                    XY{iStm} = [xNow,yNow];                    
                    break
                end
            end
        end
      
        %%% plot just-finished streamline
        switch plotCase
            case 'plot'
                plot(hAx,XY{iStm-1}(:,1) ,XY{iStm-1}(:,2) ,'-k');
                drawnow
        end

    else
        % point is okay, keep going.
        % store output
        XY{iStm}(end+1,:) = [xNow,yNow];
    end
    
    if iterCt > iterMax
        break
    end   
    
    if idxCand == nCand
        XY{iStm}(end+1,:) = [xNow,yNow];
        break
    end
        
end

% get here only when run out of possible starting locations
dispstat('Finished','timestamp','keepprev');

% keyboard
end

%% subfunctions
function [xNow,yNow,tNow] = RKstep(xGrid,yGrid,xDir,yDir,xNow,yNow,tNow,dt)
% separate out RK so we can internally handle the reverse trajectory

    %% RK Steps
    % RK Step 1
    %... k1 = f(xn,yn)
    xDirK1 = dt*interp2(xGrid,yGrid,xDir, xNow, yNow);
    yDirK1 = dt*interp2(xGrid,yGrid,yDir, xNow, yNow);

    % RK Step 2
    %... k2 = f(xn+0.5dt, yn+0.5*k1)
    xDirK2 = dt*interp2(xGrid,yGrid,xDir, xNow+0.5*xDirK1, yNow+0.5*yDirK1);
    yDirK2 = dt*interp2(xGrid,yGrid,yDir, xNow+0.5*xDirK1, yNow+0.5*yDirK1);
    % RK Step 3
    %... k3 = f(xn+0.5dt, yn+0.5k2)
    xDirK3 = dt*interp2(xGrid,yGrid,xDir, xNow+0.5*xDirK2, yNow+0.5*yDirK2);
    yDirK3 = dt*interp2(xGrid,yGrid,yDir, xNow+0.5*xDirK2, yNow+0.5*yDirK2);
    
    % RK: Step 4
    %... k4 = f(xn+h, yn+k3)
    xDirK4 = dt*interp2(xGrid,yGrid,xDir, xNow+xDirK3, yNow+yDirK3);
    yDirK4 = dt*interp2(xGrid,yGrid,yDir, xNow+xDirK3, yNow+yDirK3);

    % RK: Combine all intermediate steps
    %... x_(n+1) = x_n + (dt/6)*(k1 + 2*k2 + 2*k3 + k4)
    xNow = xNow + (1/6)*(xDirK1 + 2*xDirK2 + 2*xDirK3 + xDirK4);
    yNow = yNow + (1/6)*(yDirK1 + 2*yDirK2 + 2*yDirK3 + yDirK4);
    
    tNow = tNow + dt;

end

function [out] = checkdomain(xNow,yNow,xLim,yLim)
    %... top surface
    outTop = find(yNow > yLim(2));
    
    %... bottom surface
    outBot = find(yNow < yLim(1));
    
    %... right surface
    outRhs = find(xNow > xLim(2));
    
    %... left surface
    outLhs = find(xNow < xLim(1));
    
    %... out Nan
    outNan = find(isnan(xNow));

    
    out = [outTop(:); outBot(:); outRhs(:); outLhs(:); outNan(:)];
%     out = any(out);
        
end
    
function [] = checksamesize(data1,data2)
% check if data1 and data2 have the same dimensions

narginchk(2,2);
validateattributes(data1,{'numeric'},{});
validateattributes(data2,{'numeric'},{});

%%
nData1 = size(data1);
nData2 = size(data2);

nDims1 = numel(nData1);
nDims2 = numel(nData2);

if nDims1 ~= nDims2
    error('data do not have the same number of dimensions');
end
if sum(nData1 == nData2) ~= nDims1
    error('data have the same number of dimensions but are not the same size');
end

end


%% dispstat.m
% Copyright (c) 2013, kasim tasdemir
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are
% met:
% 
%     * Redistributions of source code must retain the above copyright
%       notice, this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above copyright
%       notice, this list of conditions and the following disclaimer in
%       the documentation and/or other materials provided with the distribution
% 
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
% POSSIBILITY OF SUCH DAMAGE.

function dispstat(TXT,varargin)
% Prints overwritable message to the command line. If you dont want to keep
% this message, call dispstat function with option 'keepthis'. If you want to
% keep the previous message, use option 'keepprev'. First argument must be
% the message.
% IMPORTANT! In the firt call, option 'init' must be used for initialization purposes.
% Options:
%     'init'      this must be called in the begining. Otherwise, it can overwrite the previous outputs on the command line.
%     'keepthis'    the message will be persistent, wont be overwritable,
%     'keepprev'  the previous message wont be overwritten. New message will start from next line,
%     'timestamp' current time hh:mm:ss will be appended to the begining of the message.
% Example:
%   clc;
%   fprintf('12345677890\n');
%   dispstat('','init')      %Initialization. Does not print anything.
%   dispstat('Time stamp will be written over this text.'); % First output
%   dispstat('is current time.','timestamp','keepthis'); % Overwrites the previous output but this output wont be overwritten.
%   dispstat(sprintf('*********\nDeveloped by %s\n*********','Kasim')); % does not overwrites the previous output
%   dispstat('','timestamp','keepprev','keepthis'); % does not overwrites the previous output
%   dispstat('this wont be overwriten','keepthis');
%   dispstat('dummy dummy dummy');
%   dispstat('final stat');
% % Output:
%     12345677890
%     15:15:34 is current time.
%     *********
%     Developed by Kasim
%     *********
%     15:15:34 
%     this wont be overwriten
%     final stat

% **********
% **** Options
keepthis = 0; % option for not overwriting
keepprev = 0;
timestamp = 0; % time stamp option
init = 0; % is it initialization step?
if ~isstr(TXT)
    return
end
persistent prevCharCnt;
if isempty(prevCharCnt)
    prevCharCnt = 0;
end
if nargin == 0
    return
elseif nargin > 1
    for i = 2:nargin
        eval([varargin{i-1} '=1;']);
    end
end
if init == 1
    prevCharCnt = 0;
    return;
end
if isempty(TXT) && timestamp == 0
    return
end
if timestamp == 1
    c = clock; % [year month day hour minute seconds]
    txtTimeStamp = sprintf('%02d:%02d:%02d ',c(4),c(5),round(c(6)));
else
    txtTimeStamp = '';
end
if keepprev == 1
    prevCharCnt = 0;
end
% *************** Make safe for fprintf, replace control charachters
TXT = strrep(TXT,'%','%%');
TXT = strrep(TXT,'\','\\');
% *************** Print
TXT = [txtTimeStamp TXT '\n'];
fprintf([repmat('\b',1, prevCharCnt) TXT]);
nof_extra = length(strfind(TXT,'%%'));
nof_extra = nof_extra + length(strfind(TXT,'\\'));
nof_extra = nof_extra + length(strfind(TXT,'\n'));
prevCharCnt = length(TXT) - nof_extra; %-1 is for \n
if keepthis == 1
    prevCharCnt = 0;
end


end


%% interparc.m
% Copyright (c) 2012, John D'Errico
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are
% met:
% 
%     * Redistributions of source code must retain the above copyright
%       notice, this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above copyright
%       notice, this list of conditions and the following disclaimer in
%       the documentation and/or other materials provided with the distribution
% 
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
% POSSIBILITY OF SUCH DAMAGE.
function [pt,dudt,fofthandle] = interparc(t,px,py,varargin)
% interparc: interpolate points along a curve in 2 or more dimensions
% usage: pt = interparc(t,px,py)    % a 2-d curve
% usage: pt = interparc(t,px,py,pz) % a 3-d curve
% usage: pt = interparc(t,px,py,pz,pw,...) % a 4-d or higher dimensional curve
% usage: pt = interparc(t,px,py,method) % a 2-d curve, method is specified
% usage: [pt,dudt,fofthandle] = interparc(t,px,py,...) % also returns derivatives, and a function handle
%
% Interpolates new points at any fractional point along
% the curve defined by a list of points in 2 or more
% dimensions. The curve may be defined by any sequence
% of non-replicated points.
%
% arguments: (input)
%  t   - vector of numbers, 0 <= t <= 1, that define
%        the fractional distance along the curve to
%        interpolate the curve at. t = 0 will generate
%        the very first point in the point list, and
%        t = 1 yields the last point in that list.
%        Similarly, t = 0.5 will yield the mid-point
%        on the curve in terms of arc length as the
%        curve is interpolated by a parametric spline.
%
%        If t is a scalar integer, at least 2, then
%        it specifies the number of equally spaced
%        points in arclength to be generated along
%        the curve.
%
%  px, py, pz, ... - vectors of length n, defining
%        points along the curve. n must be at least 2.
%        Exact Replicate points should not be present
%        in the curve, although there is no constraint
%        that the curve has replicate independent
%        variables.
%
%  method - (OPTIONAL) string flag - denotes the method
%        used to compute the points along the curve.
%
%        method may be any of 'linear', 'spline', or 'pchip',
%        or any simple contraction thereof, such as 'lin',
%        'sp', or even 'p'.
%        
%        method == 'linear' --> Uses a linear chordal
%               approximation to interpolate the curve.
%               This method is the most efficient.
%
%        method == 'pchip' --> Uses a parametric pchip
%               approximation for the interpolation
%               in arc length.
%
%        method == 'spline' --> Uses a parametric spline
%               approximation for the interpolation in
%               arc length. Generally for a smooth curve,
%               this method may be most accurate.
%
%        method = 'csape' --> if available, this tool will
%               allow a periodic spline fit for closed curves.
%               ONLY use this method if your points should
%               represent a closed curve.
%               
%               If the last point is NOT the same as the
%               first point on the curve, then the curve
%               will be forced to be periodic by this option.
%               That is, the first point will be replicated
%               onto the end.
%
%               If csape is not present in your matlab release,
%               then an error will result.
%
%        DEFAULT: 'spline'
%
%
% arguments: (output)
%  pt - Interpolated points at the specified fractional
%        distance (in arc length) along the curve.
%
%  dudt - when a second return argument is required,
%       interparc will return the parametric derivatives
%       (dx/dt, dy/dt, dz/dt, ...) as an array.
%
%  fofthandle - a function handle, taking numbers in the interval [0,1]
%       and evaluating the function at those points.
%
%       Extrapolation will not be permitted by this call.
%       Any values of t that lie outside of the interval [0,1]
%       will be clipped to the endpoints of the curve.
%
% Example:
% % Interpolate a set of unequally spaced points around
% % the perimeter of a unit circle, generating equally
% % spaced points around the perimeter.
% theta = sort(rand(15,1))*2*pi;
% theta(end+1) = theta(1);
% px = cos(theta);
% py = sin(theta);
%
% % interpolate using parametric splines
% pt = interparc(100,px,py,'spline');
%
% % Plot the result
% plot(px,py,'r*',pt(:,1),pt(:,2),'b-o')
% axis([-1.1 1.1 -1.1 1.1])
% axis equal
% grid on
% xlabel X
% ylabel Y
% title 'Points in blue are uniform in arclength around the circle'
%
%
% Example:
% % For the previous set of points, generate exactly 6
% % points around the parametric splines, verifying
% % the uniformity of the arc length interpolant.
% pt = interparc(6,px,py,'spline');
%
% % Convert back to polar form. See that the radius
% % is indeed 1, quite accurately.
% [TH,R] = cart2pol(pt(:,1),pt(:,2))
% % TH =
% %       0.86005
% %        2.1141
% %       -2.9117
% %        -1.654
% %      -0.39649
% %       0.86005
% % R =
% %             1
% %        0.9997
% %        0.9998
% %       0.99999
% %        1.0001
% %             1
%
% % Unwrap the polar angles, and difference them.
% diff(unwrap(TH))
% % ans =
% %        1.2541
% %        1.2573
% %        1.2577
% %        1.2575
% %        1.2565
%
% % Six points around the circle should be separated by
% % 2*pi/5 radians, if they were perfectly uniform. The
% % slight differences are due to the imperfect accuracy
% % of the parametric splines.
% 2*pi/5
% % ans =
% %        1.2566
%
%
% See also: arclength, spline, pchip, interp1
%
% Author: John D'Errico
% e-mail: woodchips@rochester.rr.com
% Release: 1.0
% Release date: 3/15/2010

% unpack the arguments and check for errors
if nargin < 3
  error('ARCLENGTH:insufficientarguments', ...
    'at least t, px, and py must be supplied')
end

t = t(:);
if (numel(t) == 1) && (t > 1) && (rem(t,1) == 0)
  % t specifies the number of points to be generated
  % equally spaced in arclength
  t = linspace(0,1,t)';
elseif any(t < 0) || any(t > 1)
  error('ARCLENGTH:impropert', ...
    'All elements of t must be 0 <= t <= 1')
end

% how many points will be interpolated?
nt = numel(t);

% the number of points on the curve itself
px = px(:);
py = py(:);
n = numel(px);

% are px and py both vectors of the same length?
if ~isvector(px) || ~isvector(py) || (length(py) ~= n)
  error('ARCLENGTH:improperpxorpy', ...
    'px and py must be vectors of the same length')
elseif n < 2
  error('ARCLENGTH:improperpxorpy', ...
    'px and py must be vectors of length at least 2')
end

% compose px and py into a single array. this way,
% if more dimensions are provided, the extension
% is trivial.
pxy = [px,py];
ndim = 2;

% the default method is 'linear'
method = 'spline';

% are there any other arguments?
if nargin > 3
  % there are. check the last argument. Is it a string?
  if ischar(varargin{end})
    method = varargin{end};
    varargin(end) = [];
    
    % method may be any of {'linear', 'pchip', 'spline', 'csape'.}
    % any any simple contraction thereof.
    valid = {'linear', 'pchip', 'spline', 'csape'};
    [method,errstr] = validstring(method,valid);
    if ~isempty(errstr)
      error('INTERPARC:incorrectmethod',errstr)
    end
  end
  
  % anything that remains in varargin must add
  % an additional dimension on the curve/polygon
  for i = 1:numel(varargin)
    pz = varargin{i};
    pz = pz(:);
    if numel(pz) ~= n
      error('ARCLENGTH:improperpxorpy', ...
        'pz must be of the same size as px and py')
    end
    pxy = [pxy,pz]; %#ok
  end
  
  % the final number of dimensions provided
  ndim = size(pxy,2);
end

% if csape, then make sure the first point is replicated at the end.
% also test to see if csape is available
if method(1) == 'c'
  if exist('csape','file') == 0
    error('CSAPE was requested, but you lack the necessary toolbox.')
  end
  
  p1 = pxy(1,:);
  pend = pxy(end,:);
  
  % get a tolerance on whether the first point is replicated.
  if norm(p1 - pend) > 10*eps(norm(max(abs(pxy),[],1)))
    % the two end points were not identical, so wrap the curve
    pxy(end+1,:) = p1;
    nt = nt + 1;
  end
end

% preallocate the result, pt
pt = NaN(nt,ndim);

% Compute the chordal (linear) arclength
% of each segment. This will be needed for
% any of the methods.
chordlen = sqrt(sum(diff(pxy,[],1).^2,2));

% Normalize the arclengths to a unit total
chordlen = chordlen/sum(chordlen);

% cumulative arclength
cumarc = [0;cumsum(chordlen)];

% The linear interpolant is trivial. do it as a special case
if method(1) == 'l'
  % The linear method.
  
  % which interval did each point fall in, in
  % terms of t?
  [junk,tbins] = histc(t,cumarc); %#ok
  
  % catch any problems at the ends
  tbins((tbins <= 0) | (t <= 0)) = 1;
  tbins((tbins >= n) | (t >= 1)) = n - 1;
  
  % interpolate
  s = (t - cumarc(tbins))./chordlen(tbins);
  % be nice, and allow the code to work on older releases
  % that don't have bsxfun
  pt = pxy(tbins,:) + (pxy(tbins+1,:) - pxy(tbins,:)).*repmat(s,1,ndim);
  
  % do we need to compute derivatives here?
  if nargout > 1
    dudt = (pxy(tbins+1,:) - pxy(tbins,:))./repmat(chordlen(tbins),1,ndim);
  end
  
  % do we need to create the spline as a piecewise linear function?
  if nargout > 2
    spl = cell(1,ndim);
    for i = 1:ndim
      coefs = [diff(pxy(:,i))./diff(cumarc),pxy(1:(end-1),i)];
      spl{i} = mkpp(cumarc.',coefs);
    end
    
    %create a function handle for evaluation, passing in the splines
    fofthandle = @(t) foft(t,spl);
  end
  
  % we are done at this point
  return
end

% If we drop down to here, we have either a spline
% or csape or pchip interpolant to work with.

% compute parametric splines
spl = cell(1,ndim);
spld = spl;
diffarray = [3 0 0;0 2 0;0 0 1;0 0 0];
for i = 1:ndim
  switch method
    case 'pchip'
      spl{i} = pchip(cumarc,pxy(:,i));
    case 'spline'
      spl{i} = spline(cumarc,pxy(:,i));
      nc = numel(spl{i}.coefs);
      if nc < 4
        % just pretend it has cubic segments
        spl{i}.coefs = [zeros(1,4-nc),spl{i}.coefs];
        spl{i}.order = 4;
      end
    case 'csape'
      % csape was specified, so the curve is presumed closed,
      % therefore periodic
      spl{i} = csape(cumarc,pxy(:,i),'periodic');
      nc = numel(spl{i}.coefs);
      if nc < 4
        % just pretend it has cubic segments
        spl{i}.coefs = [zeros(1,4-nc),spl{i}.coefs];
        spl{i}.order = 4;
      end
  end
  
  % and now differentiate them
  xp = spl{i};
  xp.coefs = xp.coefs*diffarray;
  xp.order = 3;
  spld{i} = xp;
end

% catch the case where there were exactly three points
% in the curve, and spline was used to generate the
% interpolant. In this case, spline creates a curve with
% only one piece, not two.
if (numel(cumarc) == 3) && (method(1) == 's')
  cumarc = spl{1}.breaks;
  n = numel(cumarc);
  chordlen = sum(chordlen);
end

% Generate the total arclength along the curve
% by integrating each segment and summing the
% results. The integration scheme does its job
% using an ode solver.

% polyarray here contains the derivative polynomials
% for each spline in a given segment
polyarray = zeros(ndim,3);
seglen = zeros(n-1,1);

% options for ode45
opts = odeset('reltol',1.e-9);
for i = 1:spl{1}.pieces
  % extract polynomials for the derivatives
  for j = 1:ndim
    polyarray(j,:) = spld{j}.coefs(i,:);
  end
  
  % integrate the arclength for the i'th segment
  % using ode45 for the integral. I could have
  % done this part with quad too, but then it
  % would not have been perfectly (numerically)
  % consistent with the next operation in this tool.
  [tout,yout] = ode45(@(t,y) segkernel(t,y),[0,chordlen(i)],0,opts); %#ok
  seglen(i) = yout(end);
end

% and normalize the segments to have unit total length
totalsplinelength = sum(seglen);
cumseglen = [0;cumsum(seglen)];

% which interval did each point fall into, in
% terms of t, but relative to the cumulative
% arc lengths along the parametric spline?
[junk,tbins] = histc(t*totalsplinelength,cumseglen); %#ok

% catch any problems at the ends
tbins((tbins <= 0) | (t <= 0)) = 1;
tbins((tbins >= n) | (t >= 1)) = n - 1;

% Do the fractional integration within each segment
% for the interpolated points. t is the parameter
% used to define the splines. It is defined in terms
% of a linear chordal arclength. This works nicely when
% a linear piecewise interpolant was used. However,
% what is asked for is an arclength interpolation
% in terms of arclength of the spline itself. Call s
% the arclength traveled along the spline.
s = totalsplinelength*t;

% the ode45 options will now include an events property
% so we can catch zero crossings.
opts = odeset('reltol',1.e-9,'events',@ode_events);

ti = t;
for i = 1:nt
  % si is the piece of arc length that we will look
  % for in this spline segment.
  si = s(i) - cumseglen(tbins(i));
  
  % extract polynomials for the derivatives
  % in the interval the point lies in
  for j = 1:ndim
    polyarray(j,:) = spld{j}.coefs(tbins(i),:);
  end
  
  % we need to integrate in t, until the integral
  % crosses the specified value of si. Because we
  % have defined totalsplinelength, the lengths will
  % be normalized at this point to a unit length.
  %
  % Start the ode solver at -si, so we will just
  % look for an event where y crosses zero.
  [tout,yout,te,ye] = ode45(@(t,y) segkernel(t,y),[0,chordlen(tbins(i))],-si,opts); %#ok
  
  % we only need that point where a zero crossing occurred
  % if no crossing was found, then we can look at each end.
  if ~isempty(te)
    ti(i) = te(1) + cumarc(tbins(i));
  else
    % a crossing must have happened at the very
    % beginning or the end, and the ode solver
    % missed it, not trapping that event.
    if abs(yout(1)) < abs(yout(end))
      % the event must have been at the start.
      ti(i) = tout(1) + cumarc(tbins(i));
    else
      % the event must have been at the end.
      ti(i) = tout(end) + cumarc(tbins(i));
    end
  end
end

% Interpolate the parametric splines at ti to get
% our interpolated value.
for L = 1:ndim
  pt(:,L) = ppval(spl{L},ti);
end

% do we need to compute first derivatives here at each point?
if nargout > 1
  dudt = zeros(nt,ndim);
  for L = 1:ndim
    dudt(:,L) = ppval(spld{L},ti);
  end
end

% create a function handle for evaluation, passing in the splines
if nargout > 2
  fofthandle = @(t) foft(t,spl);
end

% ===============================================
%  nested function for the integration kernel
% ===============================================
  function val = segkernel(t,y) %#ok
    % sqrt((dx/dt)^2 + (dy/dt)^2 + ...)
    val = zeros(size(t));
    for k = 1:ndim
      val = val + polyval(polyarray(k,:),t).^2;
    end
    val = sqrt(val);
    
  end % function segkernel

% ===============================================
%  nested function for ode45 integration events
% ===============================================
  function [value,isterminal,direction] = ode_events(t,y) %#ok
    % ode event trap, looking for zero crossings of y.
    value = y;
    isterminal = ones(size(y));
    direction = ones(size(y));
  end % function ode_events

end % mainline - interparc


% ===============================================
%       end mainline - interparc
% ===============================================
%       begin subfunctions
% ===============================================

% ===============================================
%  subfunction for evaluation at any point externally
% ===============================================
function f_t = foft(t,spl)
% tool allowing the user to evaluate the interpolant at any given point for any values t in [0,1]
pdim = numel(spl);
f_t = zeros(numel(t),pdim);

% convert t to a column vector, clipping it to [0,1] as we do.
t = max(0,min(1,t(:)));

% just loop over the splines in the cell array of splines
for i = 1:pdim
  f_t(:,i) = ppval(spl{i},t);
end
end % function foft

function [str,errorclass] = validstring(arg,valid)
% validstring: compares a string against a set of valid options
% usage: [str,errorclass] = validstring(arg,valid)
%
% If a direct hit, or any unambiguous shortening is found, that
% string is returned. Capitalization is ignored.
%
% arguments: (input)
%  arg - character string, to be tested against a list
%        of valid choices. Capitalization is ignored.
%
%  valid - cellstring array of alternative choices
%
% Arguments: (output)
%  str - string - resulting choice resolved from the
%        list of valid arguments. If no unambiguous
%        choice can be resolved, then str will be empty.
%
%  errorclass - string - A string argument that explains
%        the error. It will be one of the following
%        possibilities:
%
%        ''  --> No error. An unambiguous match for arg
%                was found among the choices.
%
%        'No match found' --> No match was found among 
%                the choices provided in valid.
%
%        'Ambiguous argument' --> At least two ambiguous
%                matches were found among those provided
%                in valid.
%        
%
% Example:
%  valid = {'off' 'on' 'The sky is falling'}
%  
%
% See also: parse_pv_pairs, strmatch, strcmpi
%
% Author: John D'Errico
% e-mail: woodchips@rochester.rr.com
% Release: 1.0
% Release date: 3/25/2010

ind = find(strncmpi(lower(arg),valid,numel(arg)));
if isempty(ind)
  % No hit found
  errorclass = 'No match found';
  str = '';
elseif (length(ind) > 1)
  % Ambiguous arg, hitting more than one of the valid options
  errorclass = 'Ambiguous argument';
  str = '';
  return
else
  errorclass = '';
  str = valid{ind};
end

end % function validstring

