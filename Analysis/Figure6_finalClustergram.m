simpMatrix = matfile('simpleMatrix.mat');
simpleMatrix = simpMatrix.simpleMatrix;
cgo = clustergram(simpleMatrix, 'Cluster', 'Row');
set(0, 'ShowHiddenHandles', 'on')
allhnds = get(0, 'Children');
h = findall(allhnds, 'Tag', 'HeatMapAxes');
set(h, 'FontSize', 14)

get(cgo)