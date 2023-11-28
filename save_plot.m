function [] = save_plot(fig,path)
    set(fig,'Visible','off');
    set(fig,'Units','inches');
    fig.Position = [100,100,5,4]; %[left, bottom, width, height] 
    screenposition = get(fig,'Position');
    set(fig,...
        'PaperPosition',[0 0 screenposition(3:4)],...
        'PaperSize',[screenposition(3:4)]);
    %print -dpdf -painters path
    print(path,'-dpdf', '-vector');

end