%This code is for data extraction from a bar chart figure
%Author: Jun Kong; jkong@gsu.edu
%Mathematics and Statistics at GSU

function figure2num_group(varargin)

imtool close all;
switch nargin
    case 0
        %input setup
        data_path = '../data/';
        filename = 'nejmsr_fig4.tiff';
        result_path = '../result/';
        
        dates=[2018 8 1 2019 5 1];

        X_ORIG = 240.566667;   %x coordinate of the first xlabel
        X_LIMIT = 1480.471429; %x coordinate of the last xlabel
        Y_ORIG = 911.404762;   %y coordinate of the first xlabel
        Y_LIMIT = 181.785714;  %y coordinate of Y limit
        X_LEGEND = 1481.366667;
        Y_LEGEND = 145.976190;
    case 4
        data_path = varargin{1};
        filename = varargin{2};
        dates=varargin{3};
        result_path = varargin{4};
    otherwise
        error('This number of arguments is not supported')
end

I = imread([data_path filename]); 
I = I(:,:,1:3);

datenumIni = datenum(dates(1:3));
datenumEnd = datenum(dates(4:6));
dates_str = datestr(datenumIni:1:datenumEnd,28);
XTickLabel = cellstr(dates_str)';
XTickLabel = unique(XTickLabel, 'stable');
TF_label = ~cellfun(@isempty, XTickLabel);          
NUM_X_LABEL = sum(TF_label);
Y_MAX_MEASURE = 60;
SQUARE_SIZE = 33;
ONE_PERC = 0.95;
WIDTH_BIN = 17;
HALF_WIDTH_LABEL = 45;
legend_label = {'Armed conflict', 'Civil unrest', 'Community resistance', 'Crime', 'Hazards'};
X_Label_Space = linspace(X_ORIG, X_LIMIT, NUM_X_LABEL);


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


r = I(:,:,1); g = I(:,:,2); b = I(:,:,3);
X = [r(:) g(:) b(:)];

bw = any( X<254, 2); 
bw = reshape(bw, size(r));

%only keep foreground in the legend area
bw(:, round(X_LEGEND):end) = false;
bw(round(Y_LEGEND):end, :) = false;


bb = regionprops(bw, 'BoundingBox', 'Extent');
% square legend
bb = cat(1, bb.BoundingBox);
TF = abs(bb(:, 3)-bb(:, 4)) <= 3 &  abs(bb(:, 3)-SQUARE_SIZE) <=2; %%
%TF = abs(bb(:, 3)-SQUARE_SIZE) <=2; %%
% rectangle legend
%bb = cat(1, bb.Extent);
%TF = bb > 0.9; %%

selected_bb = bb(TF, :); % col, row, width, height


idx = [];
for i = 1:size(selected_bb,1)
    col = ceil(selected_bb(i, 1))+1;
    row = ceil(selected_bb(i, 2))+1;
    wid = floor(selected_bb(i, 3))-1;
    hei = floor(selected_bb(i, 4))-1;
    temp = bw(row:row+hei-1, col:col+wid-1);
    %if ~all(temp(:))
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

C = [];
for i = 1: size(loc,1)
    row = loc(i, 1);
    col = loc(i, 2);
    C = [C; [I(row, col, 1), I(row, col, 2), I(row, col, 3)]];
end

delta = 20;
se = strel('line',2,90); %%
masks = false( size(bw,1), size(bw,2), size(loc,1) );
for i = 1:size(loc,1)
    r = I(:,:,1); g = I(:,:,2); b = I(:,:,3);
    masks(:,:,i) = abs(double(r) - double(C(i,1)))<=delta...
        & abs(double(g) - double(C(i,2)))<=delta...
        & abs(double(b) - double(C(i,3)))<=delta;
    
    masks(:,:,i) = imdilate(masks(:,:,i), se);
    
    %r(~masks(:,:,i))=255; g(~masks(:,:,i))=255; b(~masks(:,:,i))=255;
    %maskI = cat(3, r, g, b );
    %figure; imshow(maskI,[]);
    %imwrite(maskI, sprintf('mask_%d.jpg', i));
end


AREA_SMALL_THR = (WIDTH_BIN/2) * 1.5;
AREA_LARGE_THR = (WIDTH_BIN*5) * 2;
h = flipud(fspecial('prewitt'));
data_location = round(Y_ORIG) * ones(NUM_X_LABEL, size(loc,1));

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
    edge(:, 1:ceil(X_ORIG-HALF_WIDTH_LABEL)) = 0;
    
    %remove space below x-axis
    edge(round(Y_ORIG):end, :) = 0;
    
    %remove small edges
    L = bwlabel(edge);
    stats = regionprops(L, 'area');
    allArea = [stats.Area];
    edge_good = ismember(L, find(allArea >= AREA_SMALL_THR & allArea <= AREA_LARGE_THR) );
    L = bwlabel(edge_good);
    stats = regionprops(L, 'PixelList');
    N = max(L(:));
    
    figure;imshow(I); hold on;
    for j = 1:N
        allPixelList = stats(j);
        unique_y = unique(allPixelList.PixelList(:,2));
        [counts, centers] = hist(allPixelList.PixelList(:,2), unique_y(1):unique_y(end));
        idx = find(counts == max(counts));
        temp_y = centers(idx(end));
        bin_center = mean(allPixelList.PixelList(:,1));
        [~, ind] = sort(abs(X_Label_Space - bin_center));
        if temp_y < data_location(ind(1), i)
            data_location(ind(1), i) = temp_y;
        end
        
        plot(bin_center, temp_y, '+', 'Markersize', 16);
    end
    
    fprintf('Data of location %d out %d have been extracted.\n', i, size(loc,1));
end

%convert from y-coordinates to y-values
data = (round(Y_ORIG)-data_location) / (round(Y_ORIG)-round(Y_LIMIT)) * Y_MAX_MEASURE;
data = round(data);

%save as an excel sheet
T= table(XTickLabel', data);
writetable(T,[result_path filename '.xlsx'],'Sheet',1,'Range','A1');

%draw bar chart to verify results
%figure; b=bar(diff_data, 'stacked', 'LineStyle', 'none', 'BarWidth', 1);
figure; b=bar(data, 'grouped', 'LineStyle', 'none', 'BarWidth', 1);
for i = 1:length(b)
    %set(b(i),'FaceColor', colors(i,:));
    set(b(i),'FaceColor', double(C(i,:))./255);
end
set(gca,'YLim', [0, Y_MAX_MEASURE]);
set(gca,'XTick', [1:NUM_X_LABEL]);    
set(gca,'XTickLabel', XTickLabel);
set(gca,'XTickLabelRotation',90);
xlabel('Week of illness onset', 'FontWeight','bold', 'FontSize', 20);
ylabel('Extracted number of cases', 'FontWeight','bold', 'FontSize', 20);
legend(legend_label,'Location', 'NorthWest');
legend boxoff
print([result_path 'reconstruction_' filename], '-dpng');


if ~exist([data_path filename(1:end-4) '.mat'],'var')
    return;
end

% plot ground truth
load([data_path filename(1:end-4) '.mat'], 'barss');
figure; b=bar(barss, 'stacked', 'LineStyle', 'none');
for i = 1:length(b)
    set(b(i),'FaceColor', colors(i,:));
end
set(gca,'YLim', [0, Y_MAX_MEASURE]);
set(gca,'XTick', [1:NUM_X_LABEL]);    
set(gca,'XTickLabel', XTickLabel);
set(gca,'XTickLabelRotation',90);
ylabel('True number of cases', 'FontWeight','bold', 'FontSize', 20);
legend(legend_label,'Location', 'NorthWest');
legend boxoff

% plot difference between extracted data with ground truth
r = min(size(diff_data, 1), size(barss, 1));
diff = abs(barss(1:r,:)-diff_data(1:r,:));
figure; b=bar(diff, 'stacked', 'LineStyle', 'none');
for i = 1:length(b)
    set(b(i),'FaceColor', colors(i,:));
end
%set(gca,'YLim', [0, Y_MAX_MEASURE]);
set(gca,'XTick', [1:NUM_X_LABEL]);    
set(gca,'XTickLabel', XTickLabel);
set(gca,'XTickLabelRotation',90);
xlabel('Week of illness onset', 'FontWeight','bold', 'FontSize', 20);
ylabel('Case number difference', 'FontWeight','bold', 'FontSize', 20);
legend(legend_label,'Location', 'NorthWest');
legend boxoff
print([result_path 'diff_' filename], '-dpng');

% plot difference percentage between extracted data with ground truth
zero_TF = (barss(1:r,:) == 0);
diff_perc = diff./(barss(1:r,:)+eps) * 100;
diff_perc(zero_TF) = 0;
figure; b=bar(diff_perc, 'stacked', 'LineStyle', 'none');
for i = 1:length(b)
    set(b(i),'FaceColor', colors(i,:));
end
%set(gca,'YLim', [0, Y_MAX_MEASURE]);
set(gca,'XTick', [1:NUM_X_LABEL]);    
set(gca,'XTickLabel', XTickLabel);
set(gca,'XTickLabelRotation',90);
xlabel('Week of illness onset', 'FontWeight','bold', 'FontSize', 20);
ylabel('Case number difference percentage', 'FontWeight','bold', 'FontSize', 20);
legend(legend_label,'Location', 'NorthWest');
legend boxoff


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
% set(gca,'XTick', [1:NUM_X_LABEL]);    
% set(gca,'XTickLabel', XTickLabel);
% set(gca,'XTickLabelRotation',90);
% xlabel('Week of illness onset', 'FontWeight','bold', 'FontSize', 20);
% ylabel('Number of cases', 'FontWeight','bold', 'FontSize', 20);
% legend(legend_label,'Location', 'NorthWest');
% 
% print([result_path 'reconstruction_' filename], '-dpng');