function g=PlotErrorPatch(x,y,e,plotColor,plotMarkers)
    %x=x(:)';
    %y=y(:)';
    %e=e(:)';
    if nargin<4
        plotColor = [];
    end
    
    if nargin<5
        plotMarkers = 1;
    end
    
    if size(e, 3)==1
        yErrorTop=y+e;
        yErrorBottom=y-e;
    else
        yErrorBottom=y-e(:, :, 1);
        yErrorTop=y+e(:, :, 2);
    end
    
    numPointsToPlotErrorBars = 8;
    
    %Make the bottom run the opposite direction to plot around the eventual
    %shape of the error patch clockwise
    yErrorBottom=yErrorBottom(end:-1:1,:);
    ye=[yErrorTop;yErrorBottom];
    
    %Similarily run the x back
    xe=[x;x(end:-1:1,:)];
    xe = repmat(xe,[1 size(ye,2)/size(xe,2)]);
    x = repmat(x,[1 size(y,2)/size(x,2)]);

    % if the number of colors provided is less than the number of lines,
    % then repeat colors
    
    hStat = ishold;
    
    hold on;
    if ~all(e(:)==0)
        h=fill(xe,ye,repmat(0:size(xe,2)-1,[size(xe,1) 1]),'linestyle','none','FaceAlpha',0.25);
        
        for cc = 1:length(h)
            h(cc).FaceColor = plotColor(cc,:);
        end
        
        hAnnotation = get(h,'Annotation');

        if ~iscell(hAnnotation)
            hAnnotation = {hAnnotation};
        end

        for ii = 1:length(h)
            hLegendEntry = get(hAnnotation{ii},'LegendInformation');
            set(hLegendEntry,'IconDisplayStyle','off');
        end
    end
    
    if (size(x,1) < numPointsToPlotErrorBars) && (plotMarkers)
        if size(e, 3) == 1 
            g=errorbar(x,y,e,'marker','o');
        else
            g=errorbar(x,y,e(:, :, 1), e(:, :, 2),'marker','o');
        end
    else
        g=plot(x,y);
    end
    
    for cc = 1:length(g)
        g(cc).Color = plotColor(cc,:);
    end
    hold off;
    
    if hStat, hold on; end
end
