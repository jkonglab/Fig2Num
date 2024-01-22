%This code is for data extraction from a bar chart figure on DRC Ebola
%epidemic from Situation Reports
%Author: Jun Kong; jkong@gsu.edu
%Mathematics and Statistics at GSU

function figure2num_DRC(varargin)

imtool close all;
switch nargin
    case 0
        %input setup
        data_path = '../data/';
        filename = 'SITREP_EVD_DRC_20191126-eng.png'; %no space between legend boxes
        result_path = '../result/';
        
        dates=[2018 4 30 2019 11 18];

        X_ORIG = 64.180613;
        X_LIMIT = 2118.882353;
        Y_ORIG = 878.482187;
        Y_LIMIT = 48.407622;
        X_LEGEND = 747.299917;
        Y_LEGEND = 296.004971;
    case 4
        data_path = varargin{1};
        filename = varargin{2};
        dates=varargin{3};
        result_path = varargin{4};
    otherwise
        error('This number of arguments is not supported')
end


I = imread([data_path filename]);

datenumIni = datenum(dates(1:3));
datenumEnd = datenum(dates(4:6));
dates_str = datestr(datenumIni:21:datenumEnd,19);
NUM_X_LABEL = size(dates_str,1)-1;
NUM_BIN_BTW_X_LABEL = 2;
ADDITION_BIN_IN_END = 1;
NUM_BIN = NUM_X_LABEL*(NUM_BIN_BTW_X_LABEL+1)+1 + ADDITION_BIN_IN_END;
dates_str_Complete = repmat(blanks(size(dates_str,2)), [NUM_BIN, 1]);
dates_str_Complete(1:3:end, :) = dates_str;
XTickLabel = cellstr(dates_str_Complete)';
Y_MAX_MEASURE = 125;
SQUARE_SIZE = 40;
ONE_PERC = 0.50;
legend_label = {'Beni', 'Mabalako', 'Mandima', 'Butembo', 'Kalunguta', 'Mambasa', 'Oicha',  'Non-active zones'};

%% load image and provide four key points from users
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

bw = any( X<254, 2);
bw = reshape(bw, size(r));

bw(:, round(X_LEGEND):end) = false;
bw(round(Y_LEGEND):end, :) = false;

%separate legend squares (filename = 'SITREP_EVD_DRC_20191126-eng.png';)
h = flipud(fspecial('prewitt')); 
edge = imfilter(double(rgb2gray(I)), h,'replicate'); %bottom edge
bw(edge>0) = 0;

stats = regionprops(bw, 'BoundingBox', 'Extent');
bb = cat(1, stats.BoundingBox);
TF = abs(bb(:, 3)-bb(:, 4)) <= 20 & abs(bb(:, 3)-SQUARE_SIZE) <=2; %% filename = 'SITREP_EVD_DRC_20191126-eng.png'; %no space between legend boxes

selected_bb = bb(TF, :); % col, row, width, height

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

% identify colors of interest
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
       
       plot(bin_center, temp_y, '+', 'Markersize', 26);
       
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

%% save as an excel sheet
varNames = {'date', 'Beni', 'Mabalako', 'Mandima', 'Butembo', 'Kalunguta', 'Mambasa', 'Oicha',  'Nonactive'};
s= 'T= table(XTickLabel'',';
for i = 1:length(legend_label)
   s = [s 'diff_data(:,' num2str(i) '), '];
end
s= [s  '''VariableNames'', varNames);'];
eval(s);

writetable(T,[result_path filename '.xlsx'],'Sheet',1,'Range','A1');

%% draw bar chart to verify results
figure; b=bar(diff_data, 'stacked', 'LineStyle', 'none', 'BarWidth', 1);
for i = 1:length(b)
    set(b(i),'FaceColor', double(C(i,:))./255);
end
set(gca,'YLim', [0, Y_MAX_MEASURE]);
set(gca,'XTick', [1:NUM_BIN]);    
set(gca,'XTickLabel', XTickLabel);
set(gca,'XTickLabelRotation',90);
xlabel('Week of illness onset', 'FontWeight','bold', 'FontSize', 20);
ylabel('Extracted number of cases', 'FontWeight','bold', 'FontSize', 20);
legend(legend_label,'Location', 'NorthWest');
legend boxoff
print([result_path 'reconstruction_' filename], '-dpng');

%% If certain data points in spread sheet need to be modified, run the following code
% %read from an excel sheet
% T = readtable([ result_path filename '.xlsx'],'Sheet',1);
% XTickLabel = table2cell(T(:,1));
% diff_data = table2array(T(:,2:end));
% 
% %draw bar chart to verify results
% figure; b=bar(diff_data, 'stacked', 'LineStyle', 'none', 'BarWidth', 1);
% for i = 1:length(b)
%     set(b(i),'FaceColor', double(C(i,:))./255);
% end
% set(gca,'XTick', [1:NUM_BIN]);    
% set(gca,'XTickLabel', XTickLabel);
% set(gca,'XTickLabelRotation',90);
% xlabel('Week of illness onset', 'FontWeight','bold', 'FontSize', 20);
% ylabel('Number of cases', 'FontWeight','bold', 'FontSize', 20);
% legend(legend_label,'Location', 'NorthWest');
% 
% print([result_path 'reconstruction_' filename], '-dpng');