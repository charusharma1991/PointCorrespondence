function end2end()
    addpath('common')
    addpath('Sphere')
    load('data/Desktop.mat')
    %N = length(data);
    err = GM_RipsSphere(data,7);
    disp(err)
end
