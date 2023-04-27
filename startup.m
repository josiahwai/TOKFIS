root = fileparts(mfilename('fullpath'));

addpath(genpath([root '/matlab']))
addpath(genpath([root '/tutorials']))

% select a tokamak
warning('off','MATLAB:rmpath:DirNotFound')
rmpath(genpath([root '/tokamaks']))
warning('on','MATLAB:rmpath:DirNotFound')
d = dir([root '/tokamaks']);
d(1:2) = [];
d(~[d(:).isdir]) = [];
n = length(d);
fprintf('Choose a tokamak [default, 1]: \n')
for i = 1:n
  fprintf('[%d] %s \n', i, d(i).name)
end
itok = input('');
if isempty(itok) || ~ismember(itok, 1:n)
  itok = 1; 
end
addpath([root '/tokamaks'])


% addpaths related to the specific tokamak
tokpath = [d(itok).folder '/' d(itok).name];
addpath(tokpath)
addpath([tokpath '/eq'])
addpath([tokpath '/externaldata'])
addpath([tokpath '/matlab'])
addpath([tokpath '/tutorial_codes'])

% clear vars
clear tokpath d itok n i 