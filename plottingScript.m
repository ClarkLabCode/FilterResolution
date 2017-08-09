% all kernels for linescan

% all kernels for 2d scan


% plots for each ROI

% figure;
% subplot(2, 1, 1);
% hold on;
% for i = 1:size(kernels, 2)
%     plot(kernels(:, i));
% end
% title('Individual ROI kernels');
% xlabel('Time (frames)');
% ylabel('Kernel');
% 
% subplot(2, 1, 2);
% code taken from Emilio's PlotRoisOnMask function
% imagesc(meanMovie);axis off;axis tight; axis equal;
% colormap gray; hold on;
% roiMaskHere = roiMask;
% if size(roiMaskHere, 1) == 1
%     roiMaskHere = repmat(roiMaskHere, size(roiMaskHere, 2), 1);
% end
% roiMaskChoices = false(size(roiMaskHere));
% roiChoiceInds = unique(roiMaskHere(:));
% roiChoiceInds(roiChoiceInds==0) = [];
% roiMaskOutlines = cell(0, 0);
% roiInd = 1;
% for i = 1:length(roiChoiceInds)
%     roiMaskChoices(roiMaskHere==roiChoiceInds(i)) = true;
%     bndrs = bwboundaries(roiMaskChoices);
%     roiMaskOutlines(end+1, 1:length(bndrs)) = bndrs;
%     [indRows, indCols] = find(roiMaskChoices);
%     roiCenterOfMass(i, 1) = mean(indRows);
%     roiCenterOfMass(i, 2) = mean(indCols);
%     roiMaskChoices(roiMaskHere==roiChoiceInds(i)) = false;
%     roiInd = roiInd+length(bndrs);
% end
% for i = 1:size(roiMaskOutlines, 1)
%     for j = 1:size(roiMaskOutlines, 2)
%         if ~isempty(roiMaskOutlines{i, j})
%             lnOut = plot(roiMaskOutlines{i, j}(:, 2), roiMaskOutlines{i, j}(:, 1), 'LineWidth', 3);
%         end
%     end
%     text(roiCenterOfMass(i, 2), roiCenterOfMass(i, 1), num2str(i), 'HorizontalAlignment', 'center', 'Color',  [1 1 1]);
% end