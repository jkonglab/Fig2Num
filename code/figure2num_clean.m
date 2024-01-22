%This code is for data extraction from a bar chart figure
%Author: Jun Kong; jkong@gsu.edu
%Mathematics and Statistics at GSU

function figure2num_clean(varargin)

imtool close all;
switch nargin
    case 0
        %input setup
        data_path = '../data/';
        filename = 'figure1_20190808.png';
        result_path = '../result/';
        
        X_ORIG = 73.783871;
        X_LIMIT = 716.261290;
        Y_ORIG = 354.977419;
        Y_LIMIT = 27.958065;
        X_LEGEND = 267.848387;
        Y_LEGEND = 220.370968;

    case 2
        data_path = varargin{1};
        filename = varargin{2};
        result_path = '../result/';
    case 3
        data_path = varargin{1};
        filename = varargin{2};
        result_path = varargin{3};
    otherwise
        error('This number of arguments is not supported')
end

NUM_X_LABEL = 22;
NUM_BIN_BTW_X_LABEL = 2;
ADDITION_BIN_IN_END = 0;
NUM_BIN = NUM_X_LABEL*(NUM_BIN_BTW_X_LABEL+1)+1 + ADDITION_BIN_IN_END;
Y_MAX_MEASURE = 140;
SQUARE_SIZE = 11;
ONE_PERC = 0.95;
colors = [0 0 0.4;...
    0 0.6 0;...
    0.9290 0.6940 0.1250;...
    0.6350 0.0780 0.1840;...
    0.8500 0.3250 0.0980;...
    0.7529 0.7529 0.7529];
legend_label = {'Mabalako', 'Mandima', 'Beni', 'Katwa and Butembo', 'Musienene', 'Other health zones'};
XTickLabel = {'30-Apr', '', '', '21-May', '', '', '11-Jun', '', '', '02-Jul', '', '',...
    '23-Jul', '', '', '13-Aug', '', '', '03-Sep', '', '', '24-Sep', '', '',...
    '15-Oct', '', '', '05-Nov', '', '', '26-Nov', '', '', '17-Dec', '', '',...
    '07-Jan', '', '', '28-Jan', '', '', '18-Feb', '', '', '11-Mar', '', '',...
    '01-Apr', '', '', '22-Apr', '', '', '13-May', '', '', '03-Jun', '', '',...
    '24-Jun', '', '', '15-Jul', '', '', '05-Aug' };


% filename = 'SITREP_report_20190901.png';
% I = imread(filename);
% 
% X_ORIG = 215.346216;
% X_LIMIT = 2106.425121;
% Y_ORIG = 1078.535427; 	 
% Y_LIMIT = 41.562802;
% X_LEGEND = 652.254428; 	 
% Y_LEGEND = 235.379227; 
% NUM_X_LABEL = 22;
% NUM_BIN_BTW_X_LABEL = 2;
% ADDITION_BIN_IN_END = 3;
% NUM_BIN = NUM_X_LABEL*(NUM_BIN_BTW_X_LABEL+1)+1 + ADDITION_BIN_IN_END;
% Y_MAX_MEASURE = 140;
% SQUARE_SIZE = 31;
% ONE_PERC = 0.95;
% colors = [192, 0, 0;...
%          0, 112, 192]/255;
% legend_label = {'Confirmed', 'Probable'};
% XTickLabel = {'30-Apr', '', '', '21-May', '', '', '11-Jun', '', '', '02-Jul', '', '',...
%               '23-Jul', '', '', '13-Aug', '', '', '03-Sep', '', '', '24-Sep', '', '',...
%               '15-Oct', '', '', '05-Nov', '', '', '26-Nov', '', '', '17-Dec', '', '',...
%               '07-Jan', '', '', '28-Jan', '', '', '18-Feb', '', '', '11-Mar', '', '',...
%               '01-Apr', '', '', '22-Apr', '', '', '13-May', '', '', '03-Jun', '', '',...
%               '24-Jun', '', '', '15-Jul', '', '', '05-Aug' };

%% load image and provide four key points from users
I = imread([data_path filename]);
if ~exist('X_ORIG', 'var') || ~exist('X_LIMIT', 'var') || ~exist('Y_ORIG', 'var') || ~exist('Y_LIMIT', 'var')
    
    fprintf('please click on origin, end point on X-axis and end point on Y-axis in a moment\n');
    pause(2);
    imshow(I, []);
    set(gcf, 'Position', get(0, 'Screensize'));
    
    [x, y] = ginput(3);
    
    [~,indx] = min( abs([diff(x(1:2)), diff(x(2:3)), diff(x([1,3]))]) );
    [~,indy] = min( abs([diff(y(1:2)), diff(y(2:3)), diff(y([1,3]))]) );
    
    ind_orig = intersect([indx, mod(indx, 3)+1], [indy, mod(indy, 3)+1]);
    ind_on_yaxis = setdiff([indx, mod(indx, 3)+1], ind_orig);
    ind_on_xaxis = setdiff([indy, mod(indy, 3)+1], ind_orig);
    
    X_ORIG = x(ind_orig); X_LIMIT = x(ind_on_xaxis);
    Y_ORIG = y(ind_orig); Y_LIMIT = y(ind_on_yaxis);
    
    fprintf(' X_ORIG: %f \t X_LIMIT: %f \n Y_ORIG: %f \t Y_LIMIT: %f\n',...
        X_ORIG, X_LIMIT, Y_ORIG, Y_LIMIT);
    
    fprintf('\n please click the bottom right corner of the legend area\n');
    pause(2);
    figure(gcf);
    
    [X_LEGEND, Y_LEGEND] = ginput(1);
    fprintf(' X_LEGEND: %f \t Y_LEGEND: %f \n',X_LEGEND, Y_LEGEND);
    close gcf
end

%% find legend squares
r = I(:,:,1); g = I(:,:,2); b = I(:,:,3);
X = [r(:) g(:) b(:)];

bw = any( X<255, 2); 
bw = reshape(bw, size(r));

bb = regionprops(bw, 'BoundingBox');
bb = cat(1, bb.BoundingBox);
TF = abs(bb(:, 3)-bb(:, 4)) <= 3 &  abs(bb(:, 3)-SQUARE_SIZE) <=2;
selected_bb = bb(TF, :);

idx = [];
for i = 1:size(selected_bb,1)
    col = ceil(selected_bb(i, 1))+1;
    row = ceil(selected_bb(i, 2))+1;
    wid = floor(selected_bb(i, 3))-1;
    hei = floor(selected_bb(i, 4))-1;
    temp = bw(row:row+hei-1, col:col+wid-1);
    
    if sum(temp(:)) / length(temp(:)) <= ONE_PERC
        idx = [idx; i];
    end
end
selected_bb(idx, :) = [];
loc = round([selected_bb(:,2)+selected_bb(:,4)/2, selected_bb(:,1)+selected_bb(:,3)/2]);

% verify color location on image
figure; imshow(I,[]); hold on;
for i = 1: size(selected_bb, 1)
    plot(loc(i,2),...
        loc(i,1),...
        'kX', 'Markersize', 11);
end

%% identify colors of interest
C = [];
for i = 1: size(loc,1)
    row = loc(i, 1);
    col = loc(i, 2);
    C = [C; [I(row, col, 1), I(row, col, 2), I(row, col, 3)]];
end

%% find individual masks
delta = 20;
se = strel('line',2,90);
masks = false( size(bw,1), size(bw,2), size(loc,1) );
for i = 1:size(loc,1)
    r = I(:,:,1); g = I(:,:,2); b = I(:,:,3);
    masks(:,:,i) = abs(double(r) - double(C(i,1)))<=delta...
        & abs(double(g) - double(C(i,2)))<=delta...
        & abs(double(b) - double(C(i,3)))<=delta;
    
    masks(:,:,i) = imdilate(masks(:,:,i), se);
end

%% find bar heights
WIDTH_BIN = (X_LIMIT - X_ORIG) / NUM_BIN;
MEASURE_HALF_WIDTH = 2;
AREA_SMALL_THR = (WIDTH_BIN/2) * 2;
AREA_LARGE_THR = (WIDTH_BIN*5) * 2;
h = flipud(fspecial('prewitt'));
data_location = round(Y_ORIG) * ones(NUM_BIN, size(loc,1));

for i = 1:size(loc,1)
    edge = imfilter(masks(:,:,i), h,'replicate'); %bottom edge
    
    %remove legend
    r1 = floor(selected_bb(i,2));
    r2 = r1 + ceil(selected_bb(i,4));
    c1 = floor(selected_bb(i,1));
    c2 = c1 + ceil(selected_bb(i,3));
    edge(r1:r2, c1:c2) = 0;
    edge(1:ceil(Y_LEGEND), 1:ceil(X_LEGEND)) = 0;
    
    %remove space left to y-axis
    edge(:, 1:ceil(X_ORIG)) = 0;
    
    %remove space below x-axis
    edge(round(Y_ORIG):end, :) = 0;
    
    %remove small edges
    L = bwlabel(edge);
    stats = regionprops(L, 'area');
    allArea = [stats.Area];
    edge_good = ismember(L, find(allArea >= AREA_SMALL_THR & allArea <= AREA_LARGE_THR) );
    L = bwlabel(edge_good);
    
    figure;imshow(I); hold on;
    
    for j = 1:NUM_BIN 
       bin_center = round(X_ORIG + j*WIDTH_BIN - WIDTH_BIN/2);
           
       tempL = L(1:round(Y_ORIG), (bin_center-MEASURE_HALF_WIDTH):(bin_center+MEASURE_HALF_WIDTH));
       tempL = unique(tempL(:));
       tempL(tempL==0) = [];
       
       if length(tempL) > 1
           temp_count = zeros(length(tempL), 1);
          for k = 1:length(tempL) 
             temp_count(k,1) = max(max( L(1:round(Y_ORIG), (bin_center-MEASURE_HALF_WIDTH):(bin_center+MEASURE_HALF_WIDTH))==tempL(k) ));
          end
          [~, max_ind] = max(temp_count);
          tempL = tempL(max_ind(1));
       end
       
       if length(tempL) > 1
            error(sprintf('FIND MORE THAN 1 LABEL in a bin\n Location:%d, Bin Center:%d', i, bin_center));
       end
       if isempty(tempL)
            data_location(j, i) = round(Y_ORIG);
            continue;
       end
       
       temp_bw = zeros(size(L));
       temp_bw(1:round(Y_ORIG), (bin_center-MEASURE_HALF_WIDTH):(bin_center+MEASURE_HALF_WIDTH)) = L(1:round(Y_ORIG), (bin_center-MEASURE_HALF_WIDTH):(bin_center+MEASURE_HALF_WIDTH)) == tempL;
       allPixelList = regionprops(temp_bw, 'PixelList');
       
       unique_y = unique(allPixelList.PixelList(:,2));
       [counts, centers] = hist(allPixelList.PixelList(:,2), unique_y(1):unique_y(end)); 
       idx = find(counts == max(counts));
       temp_y = centers(idx(end));
       
       
       if ~all( L(temp_y, (bin_center-MEASURE_HALF_WIDTH):(bin_center+MEASURE_HALF_WIDTH))==tempL) 
            fprintf('\t NOT ALL LABELs in L[%d, %d:%d] has the same label\n',...
                  temp_y, bin_center-MEASURE_HALF_WIDTH, bin_center+MEASURE_HALF_WIDTH);
            data_location(j, i) = round(Y_ORIG);
            continue;  
       end
       
       plot(bin_center, temp_y, '+', 'Markersize', 16);
       
       data_location(j, i) = temp_y;
       
    end
    
 fprintf('Data of location %d out %d have been extracted.\n', i, size(loc,1));       
end


%% convert from y-coordinates to y-values
data = (round(Y_ORIG)-data_location) / (round(Y_ORIG)-round(Y_LIMIT)) * Y_MAX_MEASURE;
data = round(data);

diff_data = data(:,1);
for c = 2:size(data,2)
    M = max(data(:, 1:c-1), [], 2);
    temp = data(:,c) - M;
    temp(temp<0) = 0;
    diff_data = [diff_data, temp];
end

if ADDITION_BIN_IN_END > 0
    for i = 1:ADDITION_BIN_IN_END
            XTickLabel = [XTickLabel, {''}];
    end
end

%% save as an excel sheet
T= table(XTickLabel', diff_data);
writetable(T,[result_path filename '.xlsx'],'Sheet',1,'Range','A1');
T= table(XTickLabel', data);
writetable(T,[result_path filename '.xlsx'],'Sheet',2,'Range','A1');

%% draw bar chart to verify results
figure; b=bar(diff_data, 'stacked', 'LineStyle', 'none', 'BarWidth', 1);
for i = 1:length(b)
    set(b(i),'FaceColor', colors(i,:));
end
set(gca,'XTick', [1:NUM_BIN]);    
set(gca,'XTickLabel', XTickLabel);
set(gca,'XTickLabelRotation',90);
xlabel('Week of illness onset', 'FontWeight','bold', 'FontSize', 20);
ylabel('Number of cases', 'FontWeight','bold', 'FontSize', 20);
legend(legend_label,'Location', 'NorthWest');
print([result_path 'reconstruction_' filename], '-dpng');


%% If certain data points in spread sheet need to be modified, run the following code
% %read from an excel sheet
% T = readtable([filename '.xlsx'],'Sheet',1);
% XTickLabel = table2cell(T(:,1));
% diff_data = table2array(T(:,2:end));
% 
% %draw bar chart to verify results
% figure; b=bar(diff_data, 'stacked', 'LineStyle', 'none', 'BarWidth', 1);
% for i = 1:length(b)
%     set(b(i),'FaceColor', colors(i,:));
% end
% set(gca,'XTick', [1:NUM_BIN]);    
% set(gca,'XTickLabel', XTickLabel);
% set(gca,'XTickLabelRotation',90);
% xlabel('Week of illness onset', 'FontWeight','bold', 'FontSize', 20);
% ylabel('Number of cases', 'FontWeight','bold', 'FontSize', 20);
% legend(legend_label,'Location', 'NorthWest');
% 
% print([result_path 'reconstruction_' filename], '-dpng');