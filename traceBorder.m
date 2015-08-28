function [border,mark] = traceBorder(nacqT,acqTime,han)

i=1;
while i<=nacqT
    if i==1, col = 'r'; else col = 'y'; end
    [a,b,status] = ginput_alt(1,col);
    if isempty(status)
        break
    else
        hold on
        x(i) = a; y(i) = b;
        mark(i) = plot(x(i),y(i),'Color',col,'Marker','*','MarkerSize',6,'Parent',han);
        if i==1
            orig = get(han,'ylim');
            new = [y(1)-20 y(1)+20];
            if new(1)<orig(1), new(1) = orig(1); end
            if new(2)>orig(2), new(2) = orig(2); end
            set(han,'ylim',new);
        end
        i=i+1;
    end
end

set(han,'ylim',orig);
border = interp1(x,y,acqTime);
border = smooth(border',5);

