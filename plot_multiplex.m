function plot_multiplex(x,y,xlab,ylab,ttitle,fig)
% this function plots multiple data sets with different 
% x axis into one figure, but with just one label
% arguments: x data, y data, x label, ylabel, title, figure number

% scrsz = get(0, 'ScreenSize');
fig.h = figure(fig);
% semilogy(1,1);
% maximize(fig.h);
% set(gcf,'Position', [1 scrsz(4) scrsz(3) scrsz(4)]);
% figure(fig);
clf;
% error('jetzt')
% ax1 = gca;
% maximize(gca)
% pos1 = get(gca,'Position');
% uni1 = get(gca,'Units');
anz = size(x,1);
colors = ['g','r','b','k'];
for i_ = 1:anz
    h(i_) = axes;
%     set(gca,'Parent',fig.h);
%     set(h(i_),'Units',uni1);
%     set(h(i_),'Position',pos1);
    set(h(i_),'Units','points')
    set(gca,'Linewidth',1.5);
    pos =  get(h(i_),'Position');
    set(h(i_),'Fontsize',30);
    height = get(h(i_),'Fontsize');
    if i_ <= floor(anz/2) && i_ ~= anz
        set(h(i_),'Position',pos-[0,(ceil(height*1.5))*(floor(anz/2)-i_+0),0,0]);
        set(h(i_),'Color','none');
        set(h(i_),'Ytick',[]);
        set(h(i_),'YColor',get(gcf,'Color'));
    elseif i_ > floor(anz/2) && i_ ~= anz
        set(h(i_),'Position',pos-[0,(ceil(height*1.5))*(floor(anz/2)-i_+0),0,0]);
        set(h(i_),'Color','none');
        set(h(i_),'Ytick',[]);
        set(h(i_),'YColor',get(gcf,'Color'));    
    elseif i_ == anz
        set(h(i_),'Position',pos-[0,(ceil(height*1.5))*(floor(anz/2)-i_+0),0,0]);
    end
    set(h(i_),'XColor',colors(end-anz+i_));
    set(h(i_),'Xlim',[x(anz-i_+1,1),x(anz-i_+1,end)]);
    set(h(i_),'XTick',x(anz-i_+1,1:2:end));
    if i_ == 1
        xlabel(xlab,'Color','k');
    end
end
hold on
for i_ = 1:1
%     plot(x(1,:),y(i_,:),'Color',colors(end-i_+1),'Linewidth',2);
    plot(x(1,:),y(i_,:),'rx-','Linewidth',2.5,'Markersize',26);
end
ylabel(ylab);
% title(ttitle);
