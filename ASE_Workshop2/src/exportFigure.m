function [] = exportFigure(fig,plotFlag)
%EXPORT FIGURE - Funciton exporting all the figures generated

warning('off');
    
    if plotFlag.mesh.fullModel == true
        set(fig.mesh.fullModel,'units','centimeters','position',[0,0,10,7]);
        exportname = 'fig/meshFullModel.eps';
        exportgraphics(fig.mesh.fullModel,exportname,'ContentType','vector');
    end

    if plotFlag.rootLocus.airframe == true
        set(fig.rootLocus.airframe,'units','centimeters','position',[0,0,10,7]);
        exportname = 'fig/rootLocusAirframe.eps';
        exportgraphics(fig.rootLocus.airframe,exportname,'ContentType','vector');
    end
    
    if plotFlag.rootLocus.grdRotor == true
        set(fig.rootLocus.grdRotor,'units','centimeters','position',[0,0,10,7]);
        exportname = 'fig/rootLocusGrdRotor.eps';
        exportgraphics(fig.rootLocus.grdRotor,exportname,'ContentType','vector');
    end

    if plotFlag.rootLocus.tiltrotor == true
        set(fig.rootLocus.tiltrotor,'units','centimeters','position',[0,0,10,7]);
        exportname = 'fig/rootLocusTiltrotor.eps';
        exportgraphics(fig.rootLocus.tiltrotor,exportname,'ContentType','vector');
    end
    
    if plotFlag.rootLocus.completeXV15 == true
        set(fig.rootLocus.completeXV15,'units','centimeters','position',[0,0,12,7]);
        exportname = 'fig/rootLocusCompleteXV15.eps';
        exportgraphics(fig.rootLocus.completeXV15,exportname,'ContentType','vector');
    end

    if plotFlag.rootLocus.gimbalXV15 == true
        set(fig.rootLocus.gimbalXV15,'units','centimeters','position',[0,0,12,7]);
        exportname = 'fig/rootLocusGimbalXV15.eps';
        exportgraphics(fig.rootLocus.gimbalXV15,exportname,'ContentType','vector');
    end
    
    if plotFlag.rootLocus.reducedXV15 == true
        set(fig.rootLocus.reducedXV15,'units','centimeters','position',[0,0,12,7]);
        exportname = 'fig/rooLocusReducedXV15.eps';
        exportgraphics(fig.rootLocus.reducedXV15,exportname,'ContentType','vector');
    end

    if plotFlag.modes.grdRotor == true
        numModePlotted = length(fig.grdRotor);
        for i = 1:numModePlotted
            set(fig.grdRotor{i},'units','centimeters','position',[0,0,20,7]);
            exportname = ['fig/modesGrdRotor',num2str(i),'.eps'];
            exportgraphics(fig.grdRotor{i},exportname,'ContentType','vector');
        end
    end
    
    if plotFlag.modes.tiltRotor == true
        numModePlotted = length(fig.tiltRotor);
        for i = 1:numModePlotted
            set(fig.tiltRotor{i},'units','centimeters','position',[0,0,20,7]);
            exportname = ['fig/modesTiltRotor',num2str(i),'.eps'];
            exportgraphics(fig.tiltRotor{i},exportname,'ContentType','vector');
        end
    end
    
    if plotFlag.modes.gimbalXV15 == true
        numModePlotted = length(fig.gimbalXV15);
        for i = 1:numModePlotted
            set(fig.gimbalXV15{i},'units','centimeters','position',[0,0,20,7]);
            exportname = ['fig/modesGimbalXV15',num2str(i),'.eps'];
            exportgraphics(fig.gimbalXV15{i},exportname,'ContentType','vector');
        end
    end

    if plotFlag.modes.reducedXV15 == true
        numModePlotted = length(fig.reducedXV15);
        for i = 1:numModePlotted
            set(fig.reducedXV15{i},'units','centimeters','position',[0,0,20,7]);
            exportname = ['fig/modesReducedXV15',num2str(i),'.eps'];
            exportgraphics(fig.reducedXV15{i},exportname,'ContentType','vector');
        end
    end
   

warning('on');

end

