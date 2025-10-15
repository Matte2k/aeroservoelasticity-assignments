function [] = exportFigure(fig,plotFlag)
%EXPORT FIGURE - Funciton exporting all the figures generated

warning('off');
    
    if plotFlag.aircraft_model == true
        set(fig.aircraft,'units','centimeters','position',[0,0,7,6]);
        exportname = 'fig/aircraft.eps';
        exportgraphics(fig.aircraft,exportname,'ContentType','vector');
    end
    
    if plotFlag.original_modes == true
        numModePlotted = length(fig.original_modes);
        for i = 1:numModePlotted
            set(fig.original_modes{i},'units','centimeters','position',[0,0,20,7]);
            exportname = ['fig/original_mode',num2str(i),'.eps'];
            exportgraphics(fig.original_modes{i},exportname,'ContentType','vector');
        end
    end
    
    if plotFlag.actStiff_modes == true
        numModePlotted = length(fig.actStiff_modes);
        for i = 1:numModePlotted
            set(fig.actStiff_modes{i},'units','centimeters','position',[0,0,20,7]);
            exportname = ['fig/actStiff_mode',num2str(i),'.eps'];
            exportgraphics(fig.actStiff_modes{i},exportname,'ContentType','vector');
        end
    end
    
    if plotFlag.actDamp_modes == true
        numModePlotted = length(fig.actDamp_modes);
        for i = 1:numModePlotted
            set(fig.actDamp_modes{i},'units','centimeters','position',[0,0,20,7]);
            exportname = ['fig/actDamp_mode',num2str(i),'.eps'];
            exportgraphics(fig.actDamp_modes{i},exportname,'ContentType','vector');
        end
    end

    if plotFlag.hingeSens == true
        set(fig.hingeSens,'units','centimeters','position',[0,0,16,8]);
        exportname = 'fig/hingeSens.eps';
        exportgraphics(fig.hingeSens,exportname,'ContentType','vector');
    end

    if plotFlag.flutter_std == true
        set(fig.flutterStd.rootlocus,'units','centimeters','position',[0,0,7,6]);
        exportname = 'fig/stdFlutterRoot.eps';
        exportgraphics(fig.flutterStd.rootlocus,exportname,'ContentType','vector');

        set(fig.flutterStd.damping,'units','centimeters','position',[0,0,7,6]);
        exportname = 'fig/stdFlutterDamp.eps';
        exportgraphics(fig.flutterStd.damping,exportname,'ContentType','vector');

        set(fig.flutterStd.freq,'units','centimeters','position',[0,0,7,6]);
        exportname = 'fig/stdFlutterFreq.eps';
        exportgraphics(fig.flutterStd.freq,exportname,'ContentType','vector');
    end

    if plotFlag.flutter_cont == true
        % set(fig.flutterCont.damping,'units','centimeters','position',[0,0,7,6]);
        % exportname = 'fig/contFlutterDamp.eps';
        % exportgraphics(fig.flutterCont.damping,exportname,'ContentType','vector');
        % 
        % set(fig.flutterCont.freq,'units','centimeters','position',[0,0,7,6]);
        % exportname = 'fig/contFlutterFreq.eps';
        % exportgraphics(fig.flutterCont.freq,exportname,'ContentType','vector');
        
        set(fig.flutterCont.flutter,'units','centimeters','position',[0,0,18,8]);
        exportname = 'fig/contFlutter.eps';
        exportgraphics(fig.flutterCont.flutter,exportname,'ContentType','vector');

        % set(fig.flutterModalCont.damping,'units','centimeters','position',[0,0,7,6]);
        % exportname = 'fig/contModalFlutterDamp.eps';
        % exportgraphics(fig.flutterModalCont.damping,exportname,'ContentType','vector');
        % 
        % set(fig.flutterModalCont.freq,'units','centimeters','position',[0,0,7,6]);
        % exportname = 'fig/contModalFlutterFreq.eps';
        % exportgraphics(fig.flutterModalCont.freq,exportname,'ContentType','vector');

        set(fig.flutterModalCont.flutter,'units','centimeters','position',[0,0,18,8]);
        exportname = 'fig/contModalFlutter.eps';
        exportgraphics(fig.flutterModalCont.flutter,exportname,'ContentType','vector');
    end

    if plotFlag.flutter_sens == true
        set(fig.flutterSens,'units','centimeters','position',[0,0,15,12]);
        exportname = 'fig/flutterSens.eps';
        exportgraphics(fig.flutterSens,exportname,'ContentType','vector');

        set(fig.flutterSens,'units','centimeters','position',[0,0,15,12]);
        exportname = 'fig/flutterSens.png';
        exportgraphics(fig.flutterSens,exportname,Resolution=1500);
    end

warning('on');

end

