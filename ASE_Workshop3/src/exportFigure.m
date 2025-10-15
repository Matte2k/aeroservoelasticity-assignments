function [] = exportFigure(fig,plotFlag)
%EXPORT FIGURE - Funciton exporting all the figures generated

warning('off');
    
    if plotFlag.bodeAirframe == true
        set(fig.bodeAirframe,'units','centimeters','position',[0,0,20,8]);
        exportname = 'fig/bodeAirframe.eps';
        exportgraphics(fig.bodeAirframe,exportname,'ContentType','vector');
    end

    if plotFlag.bodeTiltrotor == true
        set(fig.bodeTiltrotor,'units','centimeters','position',[0,0,20,8]);
        exportname = 'fig/bodeTiltrotor.eps';
        exportgraphics(fig.bodeTiltrotor,exportname,'ContentType','vector');
    end
    
    if plotFlag.flutterStd == true
        set(fig.flutterStd.freq,'units','centimeters','position',[0,0,10,7]);
        exportname = 'fig/flutterStd_freq.eps';
        exportgraphics(fig.flutterStd.freq,exportname,'ContentType','vector');

        set(fig.flutterStd.damping,'units','centimeters','position',[0,0,10,7]);
        exportname = 'fig/flutterStd_damp.eps';
        exportgraphics(fig.flutterStd.damping,exportname,'ContentType','vector');
    end

    if plotFlag.flutterCont == true
        set(fig.flutterCont.vgvf,'units','centimeters','position',[0,0,20,9]);
        exportname = 'fig/flutterCont_vgvf.eps';
        exportgraphics(fig.flutterCont.vgvf,exportname,'ContentType','vector');
    end
    
    if plotFlag.freqSweep == true
        set(fig.freqSweep,'units','centimeters','position',[0,0,30,12]);
        exportname = 'fig/freqSweep.eps';
        exportgraphics(fig.freqSweep,exportname,'ContentType','vector');
    end

    if plotFlag.freqDwell == true
        set(fig.freqDwell,'units','centimeters','position',[0,0,30,12]);
        exportname = 'fig/freqDwell.eps';
        exportgraphics(fig.freqDwell,exportname,'ContentType','vector');
    end
    
    if plotFlag.logDecrement == true
        set(fig.logDecrement,'units','centimeters','position',[0,0,30,8]);
        exportname = 'fig/logDecrement.eps';
        exportgraphics(fig.logDecrement,exportname,'ContentType','vector');
    end   

warning('on');

end

