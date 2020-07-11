function main()
    addpath('common')
    addpath('Sphere')

%-------------spherical vs spherical table 2 results in main paper----------------------
    results_SS=zeros(3,5);
    load('data/Kamaishi.mat')
    err = GM_RipsSphere_SS(data,15);
    result_SS(1,1) = err;
    disp(err)
    load('data/Chessboard.mat')
    err = GM_RipsSphere_SS(data,28);
    result_SS(1,2) = err;
    disp(err)
    load('data/Desktop.mat')
    err = GM_RipsSphere_SS(data,7);
    result_SS(1,3) = err;
    disp(err)
    load('data/Parking.mat')
    err = GM_RipsSphere_SS(data,2);
    result_SS(1,4) = err;
    disp(err)
    load('data/Table.mat')
    err = GM_RipsSphere_SS(data,9);
    result_SS(1,5) = err;
    disp(err)
 

    Header = 'Kamaishi,Chessboard,Desktop,Parking,Table';
    %write header to file
    fid = fopen('Spherical_Spherical.csv','w'); 
    fprintf(fid,'%s\n',Header)
    fclose(fid)
    %write data to end of file
    dlmwrite('Spherical_Spherical.csv',result_SS,'-append');


%---------------------------------------------Cone and Ellipsoid results-----------------------------
%You can match warped images on Cone and Ellipsoid surfaces.
    addpath('Conical')
    addpath('Ellipsoidal')
    load('data/Desktop.mat')
    err = GM_RipsCone_SS(data,7);
    result_SS(2,3) = err;
    disp(err)
% %Ellipsoidal results
%     load('data/Desktop.mat')
%     err = GM_RipsEllipse_SS(data,7);
    load('Desktop_RipsEllipse_SS.mat')
    result_SS(3,3)=err;
    disp(err)


    
