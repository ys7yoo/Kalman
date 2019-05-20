


switch computer
    case 'MACI64'
        if exist('/opt/intel/','dir')
            USE_ICC = input('Compile with icc (intel compiler)? (1 or 0)');
        else
            USE_ICC = 0;
        end
    otherwise 
        USE_ICC = 0;
end
        
 
if USE_ICC
    % compile with intel compiler 
    mex CXX='icpc' estimParamEM.cpp kalman_smth_1d.cpp m_step_1d.cpp  CFLAGS="\$CFLAGS" LDFLAGS="\$LDFLAGS -L/opt/intel/lib  -lirc"  -D_MEX_EM -I../../VipMatrix ../../VipMatrix/VipMatrix.cpp

else
    % compile with gcc
    mex estimParamEM.cpp kalman_smth_1d.cpp m_step_1d.cpp  -D_MEX_EM -I../../VipMatrix ../../VipMatrix/VipMatrix.cpp
end

disp('compiled result')
ls -l estimParamEM.mexmaci64

    
    
return 


% % %% need to rune this one time
% % setenv('DYLD_LIBRARY_PATH',[getenv('DYLD_LIBRARY_PATH') ':/opt/intel/lib:/opt/intel/mkl/lib:/opt/intel/ipp/lib'])
% % setenv('PATH',[getenv('PATH') ':/opt/intel/lib:/opt/intel/mkl/lib:/opt/intel/ipp/lib'])





% call export NUM_OMP_THREAD=16  % before launching matlab 















%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% test routines from here 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% single thread version
mex estimParamEM.cpp kalman_smth_1d.cpp m_step_1d.cpp  -D_MEX_EM -I../../VipMatrix ../../VipMatrix/VipMatrix.cpp
ls -l estimParamEM.mexmaci64
%=> -rwxr-xr-x  1 yyoo  staff  45716 May 18 10:36 estimParamEM.mexmaci64
%% with icc 
mex CXX='icpc' estimParamEM.cpp kalman_smth_1d.cpp m_step_1d.cpp  CFLAGS="\$CFLAGS" LDFLAGS="\$LDFLAGS -L/opt/intel/lib  -lirc"  -D_MEX_EM -I../../VipMatrix ../../VipMatrix/VipMatrix.cpp
ls -l estimParamEM.mexmaci64
%=> -rwxr-xr-x  1 yyoo  staff  114904 May 18 10:34 estimParamEM.mexmaci64



% gcc (openMP) on unhygenix (sub, rep=10, ns=1000)
% timeDuration =
% 
%    12.7884
%     7.2463
%     6.6722
%     7.2397
%     8.7764
    

% icc (openMP) on unhygenix (sub, rep=10, ns=1000)
% timeDuration =
%    11.7927
%     6.9438
%     6.4087
%     7.1850
%     8.6981
    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Multi-thread versions 

%% multi thread version with gcc
mex estimParamEM.cpp kalman_smth_1d.cpp m_step_1d.cpp  CFLAGS="\$CFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp" -DRUN_MULTI_THREAD  -D_MEX_EM -I../../VipMatrix ../../VipMatrix/VipMatrix.cpp
ls -l estimParamEM.mexmaci64
% -rwxr-xr-x  1 yyoo  staff  57612 May 18 01:12 estimParamEM.mexmaci64



%% icc with dynamic link => run only from terminal
mex CXX='icpc' estimParamEM.cpp kalman_smth_1d.cpp m_step_1d.cpp  CFLAGS="\$CFLAGS -openmp -openmp-lib=compat" LDFLAGS="\$LDFLAGS -L/opt/intel/lib  -liomp5 -lirc" -DRUN_MULTI_THREAD  -D_MEX_EM -I../../VipMatrix ../../VipMatrix/VipMatrix.cpp
ls -l estimParamEM.mexmaci64
% vital xcode 4.1
% -rwxr-xr-x  1 yyoo  staff  113872 May 18 00:28 estimParamEM.mexmaci64
% %unhugenix xcode 4.2.1
% -rwxr-xr-x  1 yyoo  staff  115256 May 18 01:15 estimParamEM.mexmaci64





% gcc (openMP) on unhygenix 
% timeDuration =
% 
%    12.4089
%     7.1632
%     6.5879
%     7.1482
%     8.6578

% icc (openMP) on unhygenix 
% timeDuration =
% 11.7360
% 6.9487
% 6.4016
% 7.0353
% 8.4761
% not much difference ....






%%
maxNumCompThreads (24)    

setenv('OMP_NUM_THREADS','24')

%% icc with static link => crash 
mex CXX='icpc' estimParamEM.cpp kalman_smth_1d.cpp m_step_1d.cpp  CFLAGS="\$CFLAGS -openmp -openmp-lib=compat" LDFLAGS="\$LDFLAGS  /opt/intel/lib/libiomp5.a /opt/intel/lib/libirc.a" -DRUN_MULTI_THREAD  -D_MEX_EM -I../../VipMatrix ../../VipMatrix/VipMatrix.cpp
% a lot of warning but ...
% -rwxr-xr-x  1 yyoo  staff  887520 May 18 01:26 estimParamEM.mexmaci64



% OMP: Error #15: Initializing libiomp5.dylib, but found libiomp5.a already initialized.
% OMP: Hint: This may cause performance degradation and correctness issues. Set environment variable KMP_DUPLICATE_LIB_OK=TRUE to ignore this problem and force the program to continue anyway. Please note that the use of KMP_DUPLICATE_LIB_OK is unsupported and using it may cause undefined behavior. For more information, please see http://www.intel.com/software/products/support/.




%% below is test on unhygenix


% % compile with gcc on catalix
% % -rwxr-xr-x  1 yyoo  staff  56976 May 17 16:56 Kalman/estimParamEM.mexmaci64
% % 
% % compile with gcc MEX-files on unhygie
% % -rwxr-xr-x  1 yyoo  staff  60600 May 17 16:57 Kalman/estimParamEM.mexmaci64
% % 
% % compile with llvm-gcc MEX-files on unhygie
% % -rwxr-xr-x  1 yyoo  staff  57052 May 17 16:59 Kalman/estimParamEM.mexmaci64

% % compile with llvm-gcc (Xcode 4.3)
% % -rwxr-xr-x  1 yyoo  staff  57052 May 17 22:38 estimParamEM.mexmaci64

% % compile with gcc and link with mex    => no change...
% % -rwxr-xr-x  1 yyoo  staff  88244 May 17 17:18 estimParamEM.mexmaci64
return 



%% use icc  or gcc from terminal 
% % icpc -c -openmp -openmp-lib=compat estimParamEM.cpp -DRUN_MULTI_THREAD  -D_MEX_EM -I../../VipMatrix ../../VipMatrix/VipMatrix.cpp -I/Applications/MATLAB_R2012a.app/extern/include
% % 
% % g++ -c -fopenmp estimParamEM.cpp kalman_smth_1d.cpp m_step_1d.cpp  -DRUN_MULTI_THREAD  -D_MEX_EM -I../../VipMatrix ../../VipMatrix/VipMatrix.cpp -I/Applications/MATLAB_R2012a.app/extern/include
% % 
% % 
% % LINK icc => need more library 
% % 
% % LINK gcc
% % mex -cxx estimParamEM.o VipMatrix.o kalman_smth_1d.o m_step_1d.o LDFLAGS="\$LDFLAGS -fopenmp" -output estimParamEM
% % mex -cxx CXX='icpc' estimParamEM.o VipMatrix.o kalman_smth_1d.o m_step_1d.o LDFLAGS="\$LDFLAGS -fopenmp" -output estimParamEM
% % 
% % 
% % Error with icc
% % 
% % rix/VipMatrix.cpp -I/Applications/MATLAB_R2012a.app/extern/include
% % ../../VipMatrix/VipMatrix.h(61): error: invalid combination of type specifiers
% %   	typedef long long __int64;
% %   	                  ^
% % 
% % compilation aborted for estimParamEM.cpp (code 2)
% % 
% % ../../VipMatrix/VipMatrix.h(61): error: invalid combination of type specifiers
% %   	typedef long long __int64;
% %   	                  ^
% % 
% % compilation aborted for ../../VipMatrix/VipMatrix.cpp (code 2)=======


