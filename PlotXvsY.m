function PlotXvsY(x,y,varargin)
    %set up default values. all default values can be changed by varargin
    %by putting them in the command line like so
    %plotXvsY(...,'color','[0,0,1]');
    if(size(x,1) > 1)
        graphType = 'line';
    else
        graphType = 'scatter';
    end
    
    plotColor = [];
    error = [];
    lineStyle = '-';
    
    for ii = 1:2:length(varargin)
        eval([varargin{ii} '= varargin{' num2str(ii+1) '};']);
    end

    if isempty(plotColor)
        plotColor = lines(size(y,2));
    end
    
    plotColor = repmat(plotColor,[ceil(size(y,2)/size(plotColor,1)) 1]);
    plotColor = plotColor(1:size(y,2),:);
    
    switch graphType
        case 'line'
            if isempty(error)
%                 set(gca, 'ColorOrder', color, 'NextPlot', 'add');
                plottedLines = plot(x,y,'lineStyle',lineStyle);
                colorCell = mat2cell(plotColor, ones(size(plotColor, 1),1), 3);
                [plottedLines.Color] = deal(colorCell{:});
            else
                PlotErrorPatch(x,y,error,plotColor);
            end
        case 'scatter'
            if isempty(error)
                scatter(x,y,50,plotColor);
            else
                hStat = ishold;
                hold on;
                for c = 1:size(x,2)
                    scatter(x(:,c),y(:,c),50,plotColor(c,:));
                    errorbar(x(:,c),y(:,c),error(:,c),'color',plotColor(c,:),'LineStyle','none');
                end
                hold off;
                if hStat, hold on; end
            end
        case 'bar'
            colormap(plotColor);
            bar(x,y);
            
            if ~isempty(error)
                set(gca, 'ColorOrder', plotColor, 'NextPlot', 'replace');
                hStat = ishold;
                hold on;
                numbars = length(x);
                groupwidth = min(0.8, numbars/(numbars+1.5));
                relBarPos = 1:size(y, 2);
                groupShift = -groupwidth/(2*numbars) + (2*(relBarPos)-1) * groupwidth / (2*numbars);
                x = repmat(x,[1 ceil(size(y,2)/size(x,2))]);
                x = x(:,1:size(y,2));
                groupShift = repmat(groupShift, [size(x, 1), 1]);
                x = x + groupShift;
                if size(error, 3) == 1
                    errorbar(x,y,error,'LineStyle','none', 'Color', 'k');
                else
                    errorbar(x,y,error(:, :, 1), error(:, :, 2),'LineStyle','none');
                end
                hold off;
            end
            
            if hStat, hold on; end
    end
    
%     ConfAxis;
end