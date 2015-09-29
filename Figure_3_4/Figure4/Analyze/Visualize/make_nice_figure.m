function make_nice_figure(varargin)
%MAKE_NICE_FIGURE Summary of this function goes here
%   Detailed explanation goes here

% h, tit, xlab, ylab, leg

varargin(numel(varargin)+1:6) = {''};

[h, name, tit, xlab, ylab, leg] = varargin{:};

hTitle = [];
hXLabel = [];
hYLabel = [];
hLegend = [];

if ~isempty(tit)
    hTitle  = title (h, tit);
end
if ~isempty(xlab)
    hXLabel = xlabel(h, xlab);
end
if ~isempty(ylab)
    hYLabel = ylabel(h, ylab);
end
if ~isempty(leg)
    hLegend = legend(h, leg{:});
end


set( gca                       , ...
    'FontName'   , 'Helvetica' );
set([hTitle, hXLabel, hYLabel], ...
    'FontName'   , 'AvantGarde');
set([hLegend, gca]             , ...
    'FontSize'   , 8           );
set([hXLabel, hYLabel]  , ...
    'FontSize'   , 10          );
set( hTitle                    , ...
    'FontSize'   , 12          , ...
    'FontWeight' , 'bold'      );

% set(gca, ...
%   'Box'         , 'off'     , ...
%   'TickDir'     , 'out'     , ...
%   'TickLength'  , [.02 .02] , ...
%   'XMinorTick'  , 'on'      , ...
%   'YMinorTick'  , 'on'      , ...
%   'YGrid'       , 'on'      , ...
%   'XColor'      , [.3 .3 .3], ...
%   'YColor'      , [.3 .3 .3], ...
%   'YTick'       , 0:500:2500, ...
%   'LineWidth'   , 1         );


set(gcf, 'PaperPositionMode', 'auto');
if isempty(name)
    name = 'current';
end
print_to = ['Figures\' name '.eps'];
print('-depsc2', print_to);

end

