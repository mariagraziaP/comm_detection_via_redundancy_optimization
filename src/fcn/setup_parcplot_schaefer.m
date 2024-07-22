function [surfStruct, allAnnots, currName] = setup_parcplot_schaefer(nodes)

%

addpath(genpath('C:\Users\maria\Dropbox\work\VizStuff\parc_plotter-master'))


% load surface info and save as mat

% mkdir([pwd '/data/fsaverage/mat/'])

%for fff = {'sphere','smoothwm','inflated_pre','inflated'} 
for fff = {'inflated'} 
    surfStruct = load_surfStruct('C:\Users\maria\Dropbox\work\VizStuff\parc_plotter-master\data', 'fsaverage', fff{1}) ;
    fileName = ['C:\Users\maria\Dropbox\work\VizStuff\parc_plotter-master\data\fsaverage\mat\fsaverage_', fff{1}, '.mat' ] ;
%     surfStruct = load_surfStruct([pwd '/data/'],'fsaverage',fff{1}) ;
%     fileName = [pwd '/data/fsaverage/mat/fsaverage_',fff{1},'.mat' ] ;
    save(fileName,'surfStruct','-v7.3') ;
end

% load annotations
    
% initialize a map
allAnnots = containers.Map ;

currName = ['schaefer' num2str(nodes) '-yeo17']
% currName = ['Schaefer2018_' num2str(nodes) 'Parcels_7Networks_order']
tmpAnnot = load_annotStruct('C:\Users\maria\Dropbox\work\VizStuff\parc_plotter-master\data', 'fsaverage', currName) ;
%     currName = ['schaefer' num2str(iii) '-yeo17']
% tmpAnnot = load_annotStruct([pwd '/data/'],'fsaverage',currName) ;

tmpAnnot.combo_table = [ tmpAnnot.LH.ct.table(2:end,:) ; tmpAnnot.RH.ct.table(2:end,:) ] ;
tmpAnnot.roi_ids = [ tmpAnnot.LH.ct.table(2:end,5) ; tmpAnnot.RH.ct.table(2:end,5) ] ;
tmpAnnot.combo_names = [ tmpAnnot.LH.ct.struct_names(2:end) ; tmpAnnot.RH.ct.struct_names(2:end) ] ;

tmpAnnot.LH.border = get_parc_borders(...
    tmpAnnot.LH.labs,surfStruct.LH.nbrs,0) ;
tmpAnnot.RH.border = get_parc_borders(...
    tmpAnnot.RH.labs,surfStruct.RH.nbrs,0) ;

allAnnots(currName) = tmpAnnot ;



