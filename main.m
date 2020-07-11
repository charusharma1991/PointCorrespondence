function main()
    addpath('common')
    addpath('Sphere')
    addpath('Conical')
    addpath('Ellipsoidal')
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
 
    load('data/Kamaishi.mat')
    err = GM_RipsCone_SS(data,15);
    result_SS(2,1) = err;
    disp(err)
    load('data/Chessboard.mat')
    err = GM_RipsCone_SS(data,28);
    result_SS(2,2) = err;
    disp(err)
    load('data/Desktop.mat')
    err = GM_RipsCone_SS(data,7);
    result_SS(2,3) = err;
    disp(err)
    load('data/Parking.mat')
    err = GM_RipsCone_SS(data,2);
    result_SS(2,4) = err;
    disp(err)
    load('data/Table.mat')
    err = GM_RipsCone_SS(data,9);
    result_SS(2,5) = err;
    disp(err)
% %Ellipsoidal results (will take sometime to run for Chessboard and Table dataset)

%     load('data/Kamaishi.mat')
%     err = GM_RipsEllipse_SS(data,15);
    load('Kamaishi_RipsEllipse_SS.mat')
    result_SS(3,1)=err;
    disp(err)
%     load('data/Chessboard.mat')
%     err = GM_RipsEllipse_SS(data,28);
    load('Chessboard_RipsEllipse_SS.mat')
    result_SS(3,2)=err;
    disp(err)
%     load('data/Desktop.mat')
%     err = GM_RipsEllipse_SS(data,7);
    load('Desktop_RipsEllipse_SS.mat')
    result_SS(3,3)=err;
    disp(err)
%     load('data/Parking.mat')
%     err = GM_RipsEllipse_SS(data,2);
    load('Parking_RipsEllipse_SS.mat')
    result_SS(3,4)=err;
    disp(err)
%     load('data/Table.mat')
%     err = GM_RipsEllipse_SS(data,9);
    load('Table_RipsEllipse_SS.mat')
    result_SS(3,5)=err;
    disp(err)

    Header = 'Kamaishi,Chessboard,Desktop,Parking,Table';
    %write header to file
    fid = fopen('Spherical_Spherical.csv','w'); 
    fprintf(fid,'%s\n',Header)
    fclose(fid)
    %write data to end of file
    dlmwrite('Spherical_Spherical.csv',result_SS,'-append');
    
%-------------spherical vs planar table 3 results in main paper-----------------------------
    result_SP=zeros(3,3);
    load('data/Desktop.mat')
    err = GM_RipsSphere_SP(data,7,9);
    result_SP(1,1) = err;
    disp(err)
    load('data/Parking.mat')
    err = GM_RipsSphere_SP(data,2,37);
    result_SP(1,2) = err;
    disp(err)
    load('data/Table.mat')
    err = GM_RipsSphere_SP(data,9,4);
    result_SP(1,3) = err;
    disp(err)

    load('data/Desktop.mat')
    err = GM_RipsCone_SP(data,7,9);
    result_SP(2,1) = err;
    disp(err)
    load('data/Parking.mat')
    err = GM_RipsCone_SP(data,2,37);
    result_SP(2,2) = err;
    disp(err)
    load('data/Table.mat')
    err = GM_RipsCone_SP(data,9,4);
    result_SP(2,3) = err;
    disp(err)

    
% %Ellipsoidal results (will take sometime to run for Chessboard and Table dataset, so loading computed matchings)
%     load('data/Desktop.mat')
%     err = GM_RipsEllipse_SP(data,7,9);
    load('Desktop_RipsEllipse_SP.mat')
    result_SP(3,1)=err;
    disp(err)
%     load('data/Parking.mat')
%     err = GM_RipsEllipse_SP(data,2,37);
    load('Parking_RipsEllipse_SP.mat')
    result_SP(3,2)=err;
    disp(err)
%     load('data/Table.mat')
%     err = GM_RipsEllipse_SP(data,9,4);
    load('Table_RipsEllipse_SP.mat')
    result_SP(3,3)=err;
    disp(err)
    
    Header = 'Desktop,Parking,Table';
    %write header to file
    fid = fopen('Spherical_Planar.csv','w'); 
    fprintf(fid,'%s\n',Header)
    fclose(fid)
    %write data to end of file
    dlmwrite('Spherical_Planar.csv',result_SP,'-append');
end