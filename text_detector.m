clc; clear;

% ############ VARY THE PARAMETERS BELOW ########

% List of filenames on which the algorithm is to be performed
imageFiles = [
%     'image6.png                    '
%     'image3.png                    '
    'image4.png                    '
%     'image7.png                    '
    'image9.png                    '
];
for i=1:0
    imgNum = '';
    j = i;
    while j<1000
        imgNum = strcat(imgNum, '0');
        j = j*10;
    end;
    file = strcat('text_database\text_img', imgNum, int2str(i),'.png');
    imageFiles = [imageFiles; file];
end;

% Threshold for Canny Edge Detector
threshold = 0.3;

% If there are some pixels for which the actual stroke-width is not
% calculated, then set their width as the median width of the whole group

% Upper limit to the acceptable stroke width to median ratio
strokeWidthToMedianRatio = 3;
% Upper limit of the ratio of (the number of pixels not satisfying the above
% criterion to the total pixels in the group)
ratioNotCloseToMedian = 0.4;

% Lower limit on the height and width of the character
% Note the both height and width should be below threshold to restrict a
% character
minHeight = 0;
minWidth = 0;

% For a group having a bunch of stroke-widths, discard the character if
% AA * variance > BB * mean^CC
AA = 2;
BB = 1;
CC = 2;

% Discard the character not fulfilling the below ratios
heightToStrokeUpper = 15;
heightToStrokeLower = 1.5;  % To include thick letters
widthToStrokeUpper = 15;
widthToStrokeLower = 0.8;   % To include I

% ############# END OF PARAMETERS ########################

tic;
for imageFile = imageFiles'
    im_colored = imread(strtrim(imageFile'));

    im = double(rgb2gray(im_colored));
    % figure, imshow(uint8(im));
    [HEIGHT WIDTH] = size(im);

    % ###############  GAUSSIAN SMOOTHING  ###############
    % Commented code written at the end
    % Mean filter
    % imfiltered = conv2(im, ones(3)) / 9;
    % imshow(uint8(imfiltered));
    
    imfiltered = im;

    % #############  CANNY EDGE DETECTION  ##################
    % Sobel Mask
    % Gradient along x direction
    % +ve gradient from Black to White
    gx = conv2(imfiltered, [-1 0 1; -2 0 2; -1 0 1]);
    % Gradient along y direction
    % +ve gradient from Black to White
    gy = conv2(imfiltered, [-1 -2 -1; 0 0 0; 1 2 1]);

    % Gradient is always between 0 and 360 (non-negative)
    alpha = atand(gy./gx);
    alpha(gx<0) = alpha(gx<0) + 180;
    alpha = mod(alpha + 360, 360);

    M = edge(imfiltered, 'canny', threshold);
    stepsize = 1;

    % There are still some pixels where alpha is present in a neighbour pixel
    % from the pixel where M is present
    % Simply copy the nearest pixel
    for i=1:3
        for y = 1:HEIGHT
            for x = 1:WIDTH
                if not(M(y,x)==0) && isnan(alpha(y,x))
                    if not(isnan(alpha(y+1,x)))
                        alpha(y,x) = alpha(y+1,x);
                        continue;
                    end;
                    if not(isnan(alpha(y-1,x)))
                        alpha(y,x) = alpha(y-1,x);
                        continue;
                    end;
                    if not(isnan(alpha(y,x+1)))
                        alpha(y,x) = alpha(y,x+1);
                        continue;
                    end;
                    if not(isnan(alpha(y,x-1)))
                        alpha(y,x) = alpha(y,x-1);
                        continue;
                    end;
                end;
            end;
        end;
    end;

    % ##############  STROKE WIDTH TRANSFORM   ##############
    tic;
    SWT = ones(HEIGHT, WIDTH) ./ zeros(HEIGHT, WIDTH);

    % Temporary array for fast access
    pixelGroupArray = ones(HEIGHT, WIDTH)./zeros(HEIGHT, WIDTH);
    pixelGroups = zeros(1,1,2); %(group index, pixel index, (y,x) coordinates)
    pixelGroupCount = zeros(1,1);   %(group index, num pixel in that group)
    pixelGroupBound = zeros(1,2,2); %(group index, [(y1,x1) (y2,x2)] coordinates)

    pixelVisited = zeros(HEIGHT, WIDTH);
    num_groups = 0;

    for y_ite = 1:HEIGHT
        for x_ite = 1:WIDTH
            rays_f = zeros(1,1,2);
            rays_b = zeros(1,1,2);
            rays_pixelCount_f = zeros(1,1);
            rays_pixelCount_b = zeros(1,1);
            num_rays_f = 0;
            num_rays_b = 0;
            num_pixels_f = 0;
            num_pixels_b = 0;
            % Found a point on the edge which is not yet traversed
            if not(M(y_ite,x_ite) == 0) && pixelVisited(y_ite,x_ite) == false
                % Rays will be thrown from this point
                queue = zeros(1,2);
                queue(1,:) = [y_ite x_ite];
                front = 2;
                ymin = y_ite; ymax = y_ite; xmin = x_ite; xmax = x_ite;

                % Take all the points on the edge (by applying Floodfill) and
                % throw rays inwards and outwards
                count = 0;

                while not(front==1)
                    % Dequeue
                    count = count+1;
                    pixel = queue(1,:);
                    queue(1,:) = [];
                    front = front - 1;

                    y = pixel(1); x = pixel(2);

                    % Extending the bounding rectangle
                    if y < ymin
                        ymin = y;
                    end;
                    if y > ymax
                        ymax = y;
                    end;
                    if x < xmin
                        xmin = x;
                    end;
                    if x > xmax
                        xmax = x;
                    end;

                    % Search for White text on Black background
                    raylength = 0;
                    ray = zeros(1,2);   % 1 point ray (dynamic length)
                    firstTime = true;
                    rayDiscarded = false;
                    dir = alpha(y,x);
                    strokeWidth = 0;
                    ynext = y;
                    xnext = x;
                    
                    % Need a do-while here**
                    while firstTime == true || M(ynextInt, xnextInt)==0
                        xnext = xnext + stepsize*cosd(dir); xnextInt = int16(xnext);
                        ynext = ynext + stepsize*sind(dir); ynextInt = int16(ynext);
                        if (xnextInt >= 1 && xnextInt <= WIDTH && ynextInt >= 1 && ynextInt <= HEIGHT)
                            raylength = raylength + 1;
                            ray(raylength,:) = [ynextInt xnextInt];
                            strokeWidth = strokeWidth + stepsize;
                        else
                            rayDiscarded = true;
                            break;
                        end;
                        firstTime = false;
                    end;
                    if rayDiscarded == false
                        % [ynextInt xnextInt] contains the other side of the edge
                        % Check the gradient at that position
    %                     if abs(dir - (alpha(ynextInt, xnextInt) - 180)) <= 30
                            num_rays_f = num_rays_f + 1;
                            % Can this be vectorized? **
                            for i = 1:raylength
                                SWT(ray(i,1), ray(i,2)) = min(SWT(ray(i,1), ray(i,2)), strokeWidth);
                                rays_f(num_rays_f,i,:) = [ray(i,1) ray(i,2)];
                            end;
                            rays_pixelCount_f(num_rays_f) = raylength;
    %                     end
                    end;

                    % Search for Black text on White background
                    raylength = 0;
                    ray = zeros(1,2);   % 1 point ray (dynamic length)
                    firstTime = true;
                    rayDiscarded = false;
                    dir = alpha(y,x);
                    strokeWidth = 0;
                    ynext = y;
                    xnext = x;
                    % Need a do-while here**
                    while firstTime == true || M(ynextInt, xnextInt)==0
                        xnext = xnext - stepsize*cosd(dir); xnextInt = int16(xnext);
                        ynext = ynext - stepsize*sind(dir); ynextInt = int16(ynext);
                        if (xnextInt >= 1 && xnextInt <= WIDTH && ynextInt >= 1 && ynextInt <= HEIGHT)
                            raylength = raylength + 1;
                            ray(raylength,:) = [ynextInt xnextInt];
                            strokeWidth = strokeWidth + stepsize;
                        else
                            % Discard the ray if outside the boundary
                            rayDiscarded = true;
                            break;
                        end;
                        firstTime = false;
                    end;
                    if rayDiscarded == false
                        % [ynextInt xnextInt] contains the other side of the edge
                        % Check the gradient at that position
    %                     if abs(dir - (alpha(ynextInt, xnextInt) + 180)) <= 30
                            num_rays_b = num_rays_b + 1;
                            % Can this be vectorized? **
                            for i = 1:raylength
                                SWT(ray(i,1), ray(i,2)) = min(SWT(ray(i,1), ray(i,2)), strokeWidth);
                                rays_b(num_rays_b,i,:) = [ray(i,1) ray(i,2)];
                            end;
                            rays_pixelCount_b(num_rays_b) = raylength;
    %                     end;
                    end;

                    % Add another ray source to the queue if it is not already
                    % present and belongs to an edge
                    if x < WIDTH && pixelVisited(y,x+1) == false
                        % Right
                        if not(M(y,x+1)==0)
                            queue(front,:) = [y x+1];
                            front = front + 1;
                            pixelVisited(y,x+1) = true;
                        end;
                        % Bottom Right
                        if y < HEIGHT && pixelVisited(y+1,x+1) == false
                            if not(M(y+1,x+1)==0)
                                queue(front,:) = [y+1 x+1];
                                front = front + 1;
                                pixelVisited(y+1,x+1) = true;
                            end;
                        end;
                        % Top Right
                        if y > 1 && pixelVisited(y-1,x+1) == false
                            if not(M(y-1,x+1)==0)
                                queue(front,:) = [y-1 x+1];
                                front = front + 1;
                                pixelVisited(y-1,x+1) = true;
                            end;
                        end;
                    end;
                    if x > 1 && pixelVisited(y,x-1) == false
                        if not(M(y,x-1)==0)
                            queue(front,:) = [y x-1];
                            front = front + 1;
                            pixelVisited(y,x-1) = true;
                        end;
                        if y < HEIGHT && pixelVisited(y+1,x-1) == false
                            if not(M(y+1,x-1)==0)
                                queue(front,:) = [y+1 x-1];
                                front = front + 1;
                                pixelVisited(y+1,x-1) = true;
                            end;
                        end;
                        if y > 1 && pixelVisited(y-1,x-1) == false
                            if not(M(y-1,x-1)==0)
                                queue(front,:) = [y-1 x-1];
                                front = front + 1;
                                pixelVisited(y-1,x-1) = true;
                            end;
                        end;
                    end;
                    if y < HEIGHT && pixelVisited(y+1,x) == false
                        if not(M(y+1,x)==0)
                            queue(front,:) = [y+1 x];
                            front = front + 1;
                            pixelVisited(y+1,x) = true;
                        end;
                    end;
                    if y > 1 && pixelVisited(y-1,x) == false
                        if not(M(y-1,x)==0)
                            queue(front,:) = [y-1 x];
                            front = front + 1;
                            pixelVisited(y-1,x) = true;
                        end;
                    end;
                end;

                num_groups = num_groups + 1;
                pixelGroupBound(num_groups,:,:) = [ymin xmin; ymax xmax];
                h = ymax - ymin;
                w = xmax - xmin;
                
                % FILTERING
                A = zeros(1, sum(rays_pixelCount_f));
                k = 1;
                for i = 1:num_rays_f
                    for j = 1:rays_pixelCount_f(i)
                        A(k) = SWT(rays_f(i,j,1), rays_f(i,j,2));
                        k = k+1;
                    end;
                end;
                % Count the number of values having strokewidth near to mean
                % value
                B = A(A./median(A)>=strokeWidthToMedianRatio);
                sizes = [size(A) size(B)];
                % If there are not too many such values,
                if sizes(4)/sizes(2) < ratioNotCloseToMedian
                    % then set their stroke width as the median of the whole
                    A(A./median(A)>=strokeWidthToMedianRatio) = median(A);
                end;

                discarded_f = false;
                % If the group does not have enough elements
                if k<20 || (h<minHeight && w<minWidth)
                    discarded_f = true;
                end;
                % Discard the group of rays having too much variance
                if AA*var(double(int64(A))) > BB*mean(double(int64(A)))^CC
                    discarded_f = true;
                else
                    median_ = median(A);
                    % Discard if height (or width) to stroke ratio is not within
                    % the limits
                    if h/median_ < heightToStrokeLower || heightToStrokeUpper < h/median_ || w/median_ < widthToStrokeLower || widthToStrokeUpper < w/median_
                        discarded_f = true;
                    end;
                end;

                if discarded_f == true
                    for i = 1:num_rays_f
                        for j = 1:rays_pixelCount_f(i)
                            SWT(rays_f(i,j,1), rays_f(i,j,2)) = Inf;
                        end;
                    end;
                    rays_f(i,:,:) = [];
                    num_rays_f = num_rays_f - 1;
                end;

                if discarded_f == false
                    im_colored(ymin, xmin:xmax,1) = 255;
                    im_colored(ymin, xmin:xmax,2) = 0;
                    im_colored(ymin, xmin:xmax,3) = 0;
                    im_colored(ymax, xmin:xmax,1) = 255;
                    im_colored(ymax, xmin:xmax,2) = 0;
                    im_colored(ymax, xmin:xmax,3) = 0;
                    im_colored(ymin:ymax, xmin,1) = 255;
                    im_colored(ymin:ymax, xmin,2) = 0;
                    im_colored(ymin:ymax, xmin,3) = 0;
                    im_colored(ymin:ymax, xmax,1) = 255;
                    im_colored(ymin:ymax, xmax,2) = 0;
                    im_colored(ymin:ymax, xmax,3) = 0;
                end;


                A = zeros(1, sum(rays_pixelCount_b));
                k = 1;
                for i = 1:num_rays_b
                    for j = 1:rays_pixelCount_b(i)
                        A(k) = SWT(rays_b(i,j,1), rays_b(i,j,2));
                        k = k+1;
                    end;
                end;

                B = A(A./median(A)>=strokeWidthToMedianRatio);
                sizes = [size(A) size(B)];
                % If there are not too many such values,
                if sizes(4)/sizes(2) < ratioNotCloseToMedian
                    % then set their stroke width as the median of the whole
                    A(A./median(A)>=strokeWidthToMedianRatio) = median(A);
                end;

                discarded_b = false;
                % If the group does not have enough elements
                if k<20 || (h<minHeight && w<minWidth)
                    discarded_b = true;
                end;
                % Discard the group of rays having too much variance
                if AA*var(double(int64(A))) > BB*mean(double(int64(A)))^CC
                    discarded_b = true;
                else
                    median_ = median(A);
                    % Discard if height (or width) to stroke ratio is not within
                    % the limits
                    if h/median_ < heightToStrokeLower || heightToStrokeUpper < h/median_ || w/median_ < widthToStrokeLower || widthToStrokeUpper < w/median_
                        discarded_b = true;
                    end;
                end;

                if discarded_b == true
                    for i = 1:num_rays_b
                        for j = 1:rays_pixelCount_b(i)
                            SWT(rays_b(i,j,1), rays_b(i,j,2)) = Inf;
                        end;
                    end;
                    rays_b(i,:,:) = [];
                    num_rays_b = num_rays_b - 1;
                end;

                if discarded_b == false
                    % [ymin xmin; ymax xmax; h w; mean(A) var(A)]
                    im_colored(ymin-1, xmin:xmax,1) = 0;
                    im_colored(ymin-1, xmin:xmax,2) = 255;
                    im_colored(ymin-1, xmin:xmax,3) = 0;
                    im_colored(ymax+1, xmin:xmax,1) = 0;
                    im_colored(ymax+1, xmin:xmax,2) = 255;
                    im_colored(ymax+1, xmin:xmax,3) = 0;
                    im_colored(ymin:ymax, xmin-1,1) = 0;
                    im_colored(ymin:ymax, xmin-1,2) = 255;
                    im_colored(ymin:ymax, xmin-1,3) = 0;
                    im_colored(ymin:ymax, xmax+1,1) = 0;
                    im_colored(ymin:ymax, xmax+1,2) = 255;
                    im_colored(ymin:ymax, xmax+1,3) = 0;
                end;

            end;
        end;
    end;
    figure, imshow(uint8(im_colored));
    toc;
end;

toc;

% ############ GAUSSIAN FILTER #######################
% a = 5;
% sigma_sq = 3;
% imfiltered = zeros(HEIGHT, WIDTH);
% for y = 1:HEIGHT
%     for x = 1:WIDTH
%         left = x-a;
%         right = x+a;
%         top = y-a;
%         bottom = y+a;
%         if x-a < 1
%             left = 1; 
%         end
%         if x+a > WIDTH
%             right = WIDTH; 
%         end
%         if y-a < 1 
%             top = 1; 
%         end
%         if y+a > HEIGHT
%             bottom = HEIGHT; 
%         end
%         wsum_ = 0;
%         sum_ = 0;
%         for j = top:bottom
%             for i = left:right
%                 weight = exp( -((y-j)^2 + (x-i)^2)/2*sigma_sq );
%                 wsum_ = wsum_ + im(j,i) * weight;
%                 sum_ = sum_ + weight;
%             end;
%         end;
%         imfiltered(y,x) = wsum_/sum_;
%     end;
% end;