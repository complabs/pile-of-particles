
font = { 'Times', 'Courier', 'Helvetica', 'Symbol' };

% 'ZapfDingbats', 'AvantGarde', 'ZapfChancery', 'Bookman', 'NewCenturySchlbk', 'Palatino',  ...
    
weihgt = { 'normal', 'bold', 'light', 'demi' };
angle = { 'normal', 'italic', 'oblique' };

figure;
grid off;
axis off;
for i = 1:length(font)
    for j = 1:length(weihgt)
        for k = 1:length(angle)
            y = ( ( i - 1 ) * 3 + ( k - 1 ) ) / length(font) / 3;
            x = 0.01 + ( j - 1 ) / 4;
            if j == 1 && k == 1
                str = font{i};
            else
                str = [ weihgt{j}, ' ', angle{k} ];
            end
            text( x, 1-y, str, ...
                'Units', 'normalized', 'VerticalAlign', 'top', ...
                'FontName', font{i}, 'FontWeight', weihgt{j}, 'FontAngle', angle{k} ); 
        end
    end
end

lab2_main( '@exportFigure', gcf, '-depsc2', 'test' );
% close( gcf );

