
% Increase counter:

if (exist('idx', 'var'));
  idx = idx + 1;
else;
  idx = 1;
end;

% Version, title and date:

VERSION                   (idx, [1: 14])  = 'Serpent 2.1.30' ;
COMPILE_DATE              (idx, [1: 20])  = 'Feb 27 2018 16:26:15' ;
DEBUG                     (idx, 1)        = 0 ;
TITLE                     (idx, [1:  8])  = 'Untitled' ;
CONFIDENTIAL_DATA         (idx, 1)        = 0 ;
INPUT_FILE_NAME           (idx, [1:  8])  = 'slab_d2o' ;
WORKING_DIRECTORY         (idx, [1: 35])  = '/home/abrate/phytra_ANE/XS_2G/UO2_1' ;
HOSTNAME                  (idx, [1:  7])  = 'vpcen13' ;
CPU_TYPE                  (idx, [1: 41])  = 'Intel(R) Xeon(R) CPU E5-2630 v3 @ 2.40GHz' ;
CPU_MHZ                   (idx, 1)        = 4294967295.0 ;
START_DATE                (idx, [1: 24])  = 'Sat Sep 29 12:28:58 2018' ;
COMPLETE_DATE             (idx, [1: 24])  = 'Sat Sep 29 14:03:50 2018' ;

% Run parameters:

POP                       (idx, 1)        = 100000 ;
CYCLES                    (idx, 1)        = 1000 ;
SKIP                      (idx, 1)        = 500 ;
BATCH_INTERVAL            (idx, 1)        = 1 ;
SRC_NORM_MODE             (idx, 1)        = 2 ;
SEED                      (idx, 1)        = 1538216938 ;
UFS_MODE                  (idx, 1)        = 0 ;
UFS_ORDER                 (idx, 1)        = 1.00000;
NEUTRON_TRANSPORT_MODE    (idx, 1)        = 1 ;
PHOTON_TRANSPORT_MODE     (idx, 1)        = 0 ;
GROUP_CONSTANT_GENERATION (idx, 1)        = 1 ;
B1_CALCULATION            (idx, [1:  3])  = [ 0 0 0 ];
B1_BURNUP_CORRECTION      (idx, 1)        = 0 ;
IMPLICIT_REACTION_RATES   (idx, 1)        = 1 ;

% Optimization:

OPTIMIZATION_MODE         (idx, 1)        = 4 ;
RECONSTRUCT_MICROXS       (idx, 1)        = 1 ;
RECONSTRUCT_MACROXS       (idx, 1)        = 1 ;
DOUBLE_INDEXING           (idx, 1)        = 0 ;
MG_MAJORANT_MODE          (idx, 1)        = 0 ;

% Parallelization:

MPI_TASKS                 (idx, 1)        = 1 ;
OMP_THREADS               (idx, 1)        = 30 ;
MPI_REPRODUCIBILITY       (idx, 1)        = 0 ;
OMP_REPRODUCIBILITY       (idx, 1)        = 1 ;
OMP_HISTORY_PROFILE       (idx, [1:  30]) = [  1.08888E+00  9.94430E-01  1.02065E+00  1.01759E+00  1.01579E+00  9.71243E-01  9.57860E-01  1.03162E+00  1.05536E+00  9.61553E-01  9.87352E-01  9.60729E-01  9.74345E-01  1.03210E+00  9.67696E-01  1.01849E+00  9.58672E-01  1.02126E+00  9.57751E-01  1.01910E+00  9.97507E-01  1.02151E+00  1.01712E+00  1.02707E+00  9.57848E-01  1.01812E+00  9.61310E-01  1.03389E+00  9.55931E-01  9.97215E-01  ];
SHARE_BUF_ARRAY           (idx, 1)        = 0 ;
SHARE_RES2_ARRAY          (idx, 1)        = 1 ;

% File paths:

XS_DATA_FILE_PATH         (idx, [1: 46])  = '/opt/serpent/xsdata/jeff311/sss_jeff311.xsdata' ;
DECAY_DATA_FILE_PATH      (idx, [1:  3])  = 'N/A' ;
SFY_DATA_FILE_PATH        (idx, [1:  3])  = 'N/A' ;
NFY_DATA_FILE_PATH        (idx, [1:  3])  = 'N/A' ;
BRA_DATA_FILE_PATH        (idx, [1:  3])  = 'N/A' ;

% Collision and reaction sampling (neutrons/photons):

MIN_MACROXS               (idx, [1:   4]) = [  5.00000E-02 0.0E+00  0.00000E+00 0.0E+00 ];
DT_THRESH                 (idx, [1:  2])  = [  9.00000E-01  9.00000E-01 ];
ST_FRAC                   (idx, [1:   4]) = [  4.45997E-04 0.00038  0.00000E+00 0.0E+00 ];
DT_FRAC                   (idx, [1:   4]) = [  9.99554E-01 1.7E-07  0.00000E+00 0.0E+00 ];
DT_EFF                    (idx, [1:   4]) = [  6.17261E-01 2.1E-05  0.00000E+00 0.0E+00 ];
REA_SAMPLING_EFF          (idx, [1:   4]) = [  1.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
REA_SAMPLING_FAIL         (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
TOT_COL_EFF               (idx, [1:   4]) = [  6.17425E-01 2.1E-05  0.00000E+00 0.0E+00 ];
AVG_TRACKING_LOOPS        (idx, [1:   8]) = [  2.18175E+00 3.9E-05  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
AVG_TRACKS                (idx, [1:   4]) = [  2.42050E+02 0.00030  0.00000E+00 0.0E+00 ];
AVG_REAL_COL              (idx, [1:   4]) = [  2.42024E+02 0.00030  0.00000E+00 0.0E+00 ];
AVG_VIRT_COL              (idx, [1:   4]) = [  1.49968E+02 0.00036  0.00000E+00 0.0E+00 ];
AVG_SURF_CROSS            (idx, [1:   4]) = [  3.27516E-02 0.00053  0.00000E+00 0.0E+00 ];
LOST_PARTICLES            (idx, 1)        = 0 ;

% Run statistics:

CYCLE_IDX                 (idx, 1)        = 1000 ;
SOURCE_POPULATION         (idx, 1)        = 100000976 ;
MEAN_POP_SIZE             (idx, [1:  2])  = [  1.00001E+05 0.00019 ];
MEAN_POP_WGT              (idx, [1:  2])  = [  1.00001E+05 0.00019 ];
SIMULATION_COMPLETED      (idx, 1)        = 1 ;

% Running times:

TOT_CPU_TIME              (idx, 1)        =  2.70421E+03 ;
RUNNING_TIME              (idx, 1)        =  9.48730E+01 ;
INIT_TIME                 (idx, [1:  2])  = [  3.18833E-02  3.18833E-02 ];
PROCESS_TIME              (idx, [1:  2])  = [  2.66668E-04  2.66668E-04 ];
TRANSPORT_CYCLE_TIME      (idx, [1:  3])  = [  9.48409E+01  9.48409E+01  0.00000E+00 ];
MPI_OVERHEAD_TIME         (idx, [1:  2])  = [  0.00000E+00  0.00000E+00 ];
ESTIMATED_RUNNING_TIME    (idx, [1:  2])  = [  9.48715E+01  0.00000E+00 ];
CPU_USAGE                 (idx, 1)        = 28.50343 ;
TRANSPORT_CPU_USAGE       (idx, [1:   2]) = [  2.87208E+01 7.9E-05 ];
OMP_PARALLEL_FRAC         (idx, 1)        =  9.50408E-01 ;

% Memory usage:

AVAIL_MEM                 (idx, 1)        = 96872.39 ;
ALLOC_MEMSIZE             (idx, 1)        = 1036.29;
MEMSIZE                   (idx, 1)        = 752.09;
XS_MEMSIZE                (idx, 1)        = 86.48;
MAT_MEMSIZE               (idx, 1)        = 9.01;
RES_MEMSIZE               (idx, 1)        = 3.37;
MISC_MEMSIZE              (idx, 1)        = 653.23;
UNKNOWN_MEMSIZE           (idx, 1)        = 0.00;
UNUSED_MEMSIZE            (idx, 1)        = 284.20;

% Geometry parameters:

TOT_CELLS                 (idx, 1)        = 7 ;
UNION_CELLS               (idx, 1)        = 0 ;

% Neutron energy grid:

NEUTRON_ERG_TOL           (idx, 1)        =  0.00000E+00 ;
NEUTRON_ERG_NE            (idx, 1)        = 97927 ;
NEUTRON_EMIN              (idx, 1)        =  1.00000E-11 ;
NEUTRON_EMAX              (idx, 1)        =  2.00000E+01 ;

% Unresolved resonance probability table sampling:

URES_DILU_CUT             (idx, 1)        =  1.00000E-09 ;
URES_EMIN                 (idx, 1)        =  1.00000E+37 ;
URES_EMAX                 (idx, 1)        = -1.00000E+37 ;
URES_AVAIL                (idx, 1)        = 2 ;
URES_USED                 (idx, 1)        = 0 ;

% Nuclides and reaction channels:

TOT_NUCLIDES              (idx, 1)        = 5 ;
TOT_TRANSPORT_NUCLIDES    (idx, 1)        = 5 ;
TOT_DOSIMETRY_NUCLIDES    (idx, 1)        = 0 ;
TOT_DECAY_NUCLIDES        (idx, 1)        = 0 ;
TOT_PHOTON_NUCLIDES       (idx, 1)        = 0 ;
TOT_REA_CHANNELS          (idx, 1)        = 114 ;
TOT_TRANSMU_REA           (idx, 1)        = 0 ;

% Neutron physics options:

USE_DELNU                 (idx, 1)        = 1 ;
USE_URES                  (idx, 1)        = 0 ;
USE_DBRC                  (idx, 1)        = 0 ;
IMPL_CAPT                 (idx, 1)        = 0 ;
IMPL_NXN                  (idx, 1)        = 1 ;
IMPL_FISS                 (idx, 1)        = 0 ;
DOPPLER_PREPROCESSOR      (idx, 1)        = 0 ;
TMS_MODE                  (idx, 1)        = 0 ;
SAMPLE_FISS               (idx, 1)        = 1 ;
SAMPLE_CAPT               (idx, 1)        = 1 ;
SAMPLE_SCATT              (idx, 1)        = 1 ;

% Radioactivity data:

TOT_ACTIVITY              (idx, 1)        =  0.00000E+00 ;
TOT_DECAY_HEAT            (idx, 1)        =  0.00000E+00 ;
TOT_SF_RATE               (idx, 1)        =  0.00000E+00 ;
ACTINIDE_ACTIVITY         (idx, 1)        =  0.00000E+00 ;
ACTINIDE_DECAY_HEAT       (idx, 1)        =  0.00000E+00 ;
FISSION_PRODUCT_ACTIVITY  (idx, 1)        =  0.00000E+00 ;
FISSION_PRODUCT_DECAY_HEAT(idx, 1)        =  0.00000E+00 ;
INHALATION_TOXICITY       (idx, 1)        =  0.00000E+00 ;
INGESTION_TOXICITY        (idx, 1)        =  0.00000E+00 ;
ACTINIDE_INH_TOX          (idx, 1)        =  0.00000E+00 ;
ACTINIDE_ING_TOX          (idx, 1)        =  0.00000E+00 ;
FISSION_PRODUCT_INH_TOX   (idx, 1)        =  0.00000E+00 ;
FISSION_PRODUCT_ING_TOX   (idx, 1)        =  0.00000E+00 ;
SR90_ACTIVITY             (idx, 1)        =  0.00000E+00 ;
TE132_ACTIVITY            (idx, 1)        =  0.00000E+00 ;
I131_ACTIVITY             (idx, 1)        =  0.00000E+00 ;
I132_ACTIVITY             (idx, 1)        =  0.00000E+00 ;
CS134_ACTIVITY            (idx, 1)        =  0.00000E+00 ;
CS137_ACTIVITY            (idx, 1)        =  0.00000E+00 ;
PHOTON_DECAY_SOURCE       (idx, 1)        =  0.00000E+00 ;
NEUTRON_DECAY_SOURCE      (idx, 1)        =  0.00000E+00 ;
ALPHA_DECAY_SOURCE        (idx, 1)        =  0.00000E+00 ;
ELECTRON_DECAY_SOURCE     (idx, 1)        =  0.00000E+00 ;

% Normalization coefficient:

NORM_COEF                 (idx, [1:   4]) = [  9.96625E-06 8.7E-05  0.00000E+00 0.0E+00 ];

% Analog reaction rate estimators:

CONVERSION_RATIO          (idx, [1:   2]) = [  1.12808E+00 0.00022 ];
U235_FISS                 (idx, [1:   4]) = [  3.38897E-01 0.00016  8.43815E-01 6.7E-05 ];
U238_FISS                 (idx, [1:   4]) = [  6.27283E-02 0.00041  1.56185E-01 0.00036 ];
U235_CAPT                 (idx, [1:   4]) = [  8.08771E-02 0.00036  1.41112E-01 0.00034 ];
U238_CAPT                 (idx, [1:   4]) = [  4.73595E-01 0.00015  8.26313E-01 6.5E-05 ];

% Neutron balance (particles/weight):

BALA_SRC_NEUTRON_SRC     (idx, [1:  2])  = [ 0 0.00000E+00 ];
BALA_SRC_NEUTRON_FISS    (idx, [1:  2])  = [ 100000976 1.00000E+08 ];
BALA_SRC_NEUTRON_NXN     (idx, [1:  2])  = [ 0 3.51983E+05 ];
BALA_SRC_NEUTRON_VR      (idx, [1:  2])  = [ 0 0.00000E+00 ];
BALA_SRC_NEUTRON_TOT     (idx, [1:  2])  = [ 100000976 1.00352E+08 ];

BALA_LOSS_NEUTRON_CAPT    (idx, [1:  2])  = [ 57302596 5.75083E+07 ];
BALA_LOSS_NEUTRON_FISS    (idx, [1:  2])  = [ 40165004 4.02985E+07 ];
BALA_LOSS_NEUTRON_LEAK    (idx, [1:  2])  = [ 2533376 2.54514E+06 ];
BALA_LOSS_NEUTRON_CUT     (idx, [1:  2])  = [ 0 0.00000E+00 ];
BALA_LOSS_NEUTRON_ERR     (idx, [1:  2])  = [ 0 0.00000E+00 ];
BALA_LOSS_NEUTRON_TOT     (idx, [1:  2])  = [ 100000976 1.00352E+08 ];

BALA_NEUTRON_DIFF         (idx, [1:  2])  = [ 0 -1.02967E-05 ];

% Normalized total reaction rates (neutrons):

TOT_POWER                 (idx, [1:   2]) = [  1.30508E-11 8.0E-05 ];
TOT_POWDENS               (idx, [1:   2]) = [  0.00000E+00 0.0E+00 ];
TOT_GENRATE               (idx, [1:   2]) = [  1.00255E+00 7.9E-05 ];
TOT_FISSRATE              (idx, [1:   2]) = [  4.01548E-01 8.0E-05 ];
TOT_CAPTRATE              (idx, [1:   2]) = [  5.73087E-01 5.8E-05 ];
TOT_ABSRATE               (idx, [1:   2]) = [  9.74635E-01 1.5E-05 ];
TOT_SRCRATE               (idx, [1:   2]) = [  9.96625E-01 8.7E-05 ];
TOT_FLUX                  (idx, [1:   2]) = [  5.10661E+02 0.00029 ];
TOT_PHOTON_PRODRATE       (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
TOT_LEAKRATE              (idx, [1:   2]) = [  2.53654E-02 0.00059 ];
ALBEDO_LEAKRATE           (idx, [1:   2]) = [  0.00000E+00 0.0E+00 ];
TOT_LOSSRATE              (idx, [1:   2]) = [  1.00000E+00 0.0E+00 ];
TOT_CUTRATE               (idx, [1:   2]) = [  0.00000E+00 0.0E+00 ];
TOT_RR                    (idx, [1:   2]) = [  2.42265E+02 0.00030 ];
INI_FMASS                 (idx, 1)        =  0.00000E+00 ;
TOT_FMASS                 (idx, 1)        =  0.00000E+00 ;

% Six-factor formula:

SIX_FF_ETA                (idx, [1:   2]) = [  1.83338E+00 0.00010 ];
SIX_FF_F                  (idx, [1:   2]) = [  9.50966E-01 4.1E-05 ];
SIX_FF_P                  (idx, [1:   2]) = [  3.04073E-01 0.00016 ];
SIX_FF_EPSILON            (idx, [1:   2]) = [  1.94747E+00 0.00016 ];
SIX_FF_LF                 (idx, [1:   2]) = [  1.00000E+00 0.0E+00 ];
SIX_FF_LT                 (idx, [1:   2]) = [  9.74549E-01 1.5E-05 ];
SIX_FF_KINF               (idx, [1:   2]) = [  1.03242E+00 0.00013 ];
SIX_FF_KEFF               (idx, [1:   2]) = [  1.00614E+00 0.00013 ];

% Fission neutron and energy production:

NUBAR                     (idx, [1:   2]) = [  2.49670E+00 5.9E-06 ];
FISSE                     (idx, [1:   2]) = [  2.02856E+02 5.4E-07 ];

% Criticality eigenvalues:

ANA_KEFF                  (idx, [1:   6]) = [  1.00619E+00 0.00013  9.99121E-01 0.00013  7.02222E-03 0.00182 ];
IMP_KEFF                  (idx, [1:   2]) = [  1.00608E+00 8.0E-05 ];
COL_KEFF                  (idx, [1:   2]) = [  1.00595E+00 0.00011 ];
ABS_KEFF                  (idx, [1:   2]) = [  1.00608E+00 8.0E-05 ];
ABS_KINF                  (idx, [1:   2]) = [  1.03236E+00 7.8E-05 ];
GEOM_ALBEDO               (idx, [1:   6]) = [  1.00000E+00 0.0E+00  1.00000E+00 0.0E+00  1.00000E+00 0.0E+00 ];

% ALF (Average lethargy of neutrons causing fission):
% Based on E0 = 2.000000E+01 MeV

ANA_ALF                   (idx, [1:   2]) = [  1.34865E+01 8.9E-05 ];
IMP_ALF                   (idx, [1:   2]) = [  1.34871E+01 6.8E-05 ];

% EALF (Energy corresponding to average lethargy of neutrons causing fission):

ANA_EALF                  (idx, [1:   2]) = [  2.78130E-05 0.00119 ];
IMP_EALF                  (idx, [1:   2]) = [  2.77881E-05 0.00092 ];

% AFGE (Average energy of neutrons causing fission):

ANA_AFGE                  (idx, [1:   2]) = [  5.25106E-01 0.00038 ];
IMP_AFGE                  (idx, [1:   2]) = [  5.24974E-01 0.00022 ];

% Forward-weighted delayed neutron parameters:

FWD_ANA_BETA_ZERO         (idx, [1:  18]) = [  8.37638E-03 0.00117  2.05386E-04 0.00713  1.12987E-03 0.00301  6.23115E-04 0.00398  1.48908E-03 0.00268  2.62622E-03 0.00198  1.07392E-03 0.00304  8.32284E-04 0.00353  3.96514E-04 0.00490 ];
FWD_ANA_LAMBDA            (idx, [1:  18]) = [  5.38799E-01 0.00161  1.24667E-02 7.2E-09  2.82917E-02 0.0E+00  4.25244E-02 0.0E+00  1.33042E-01 0.0E+00  2.92467E-01 0.0E+00  6.66488E-01 2.7E-09  1.63478E+00 4.2E-09  3.55460E+00 2.9E-09 ];

% Beta-eff using Meulekamp's method:

ADJ_MEULEKAMP_BETA_EFF    (idx, [1:  18]) = [  7.12372E-03 0.00181  1.77695E-04 0.01153  9.83290E-04 0.00484  5.31750E-04 0.00667  1.27579E-03 0.00431  2.24693E-03 0.00315  8.90216E-04 0.00496  6.93766E-04 0.00571  3.24278E-04 0.00817 ];
ADJ_MEULEKAMP_LAMBDA      (idx, [1:  18]) = [  5.27878E-01 0.00263  1.24667E-02 7.2E-09  2.82917E-02 0.0E+00  4.25244E-02 0.0E+00  1.33042E-01 0.0E+00  2.92467E-01 0.0E+00  6.66488E-01 2.7E-09  1.63478E+00 4.2E-09  3.55460E+00 2.8E-09 ];

% Adjoint weighted time constants using Nauchi's method:

ADJ_NAUCHI_GEN_TIME       (idx, [1:   6]) = [  1.90111E-03 0.00063  1.90002E-03 0.00064  2.05477E-03 0.00672 ];
ADJ_NAUCHI_LIFETIME       (idx, [1:   6]) = [  1.91284E-03 0.00061  1.91175E-03 0.00062  2.06739E-03 0.00671 ];
ADJ_NAUCHI_BETA_EFF       (idx, [1:  18]) = [  6.98322E-03 0.00186  1.75027E-04 0.01249  9.70819E-04 0.00514  5.16308E-04 0.00674  1.24334E-03 0.00450  2.20892E-03 0.00334  8.68910E-04 0.00528  6.81497E-04 0.00616  3.18387E-04 0.00879 ];
ADJ_NAUCHI_LAMBDA         (idx, [1:  18]) = [  5.28224E-01 0.00282  1.24667E-02 7.1E-09  2.82917E-02 0.0E+00  4.25244E-02 0.0E+00  1.33042E-01 0.0E+00  2.92467E-01 0.0E+00  6.66488E-01 2.7E-09  1.63478E+00 4.2E-09  3.55460E+00 2.9E-09 ];

% Adjoint weighted time constants using IFP:

ADJ_IFP_GEN_TIME          (idx, [1:   6]) = [  2.04747E-03 0.00151  2.04597E-03 0.00151  2.25135E-03 0.01799 ];
ADJ_IFP_LIFETIME          (idx, [1:   6]) = [  2.06010E-03 0.00150  2.05860E-03 0.00150  2.26525E-03 0.01799 ];
ADJ_IFP_IMP_BETA_EFF      (idx, [1:  18]) = [  7.06305E-03 0.00648  1.78447E-04 0.03926  9.69265E-04 0.01762  5.34667E-04 0.02252  1.23453E-03 0.01525  2.25953E-03 0.01137  8.98532E-04 0.01748  6.83517E-04 0.02047  3.04563E-04 0.02884 ];
ADJ_IFP_IMP_LAMBDA        (idx, [1:  18]) = [  5.20992E-01 0.00941  1.24667E-02 0.0E+00  2.82917E-02 0.0E+00  4.25244E-02 0.0E+00  1.33042E-01 0.0E+00  2.92467E-01 0.0E+00  6.66488E-01 2.7E-09  1.63478E+00 4.4E-09  3.55460E+00 2.7E-09 ];
ADJ_IFP_ANA_BETA_EFF      (idx, [1:  18]) = [  7.06705E-03 0.00637  1.76679E-04 0.03824  9.71410E-04 0.01723  5.35958E-04 0.02184  1.24027E-03 0.01514  2.25908E-03 0.01108  8.96590E-04 0.01712  6.82861E-04 0.02004  3.04202E-04 0.02780 ];
ADJ_IFP_ANA_LAMBDA        (idx, [1:  18]) = [  5.20340E-01 0.00913  1.24667E-02 0.0E+00  2.82917E-02 0.0E+00  4.25244E-02 0.0E+00  1.33042E-01 0.0E+00  2.92467E-01 0.0E+00  6.66488E-01 2.7E-09  1.63478E+00 4.4E-09  3.55460E+00 2.8E-09 ];
ADJ_IFP_ROSSI_ALPHA       (idx, [1:   2]) = [ -3.45836E+00 0.00659 ];

% Adjoint weighted time constants using perturbation technique:

ADJ_PERT_GEN_TIME         (idx, [1:   2]) = [  1.98808E-03 0.00032 ];
ADJ_PERT_LIFETIME         (idx, [1:   2]) = [  2.00035E-03 0.00029 ];
ADJ_PERT_BETA_EFF         (idx, [1:   2]) = [  7.08691E-03 0.00123 ];
ADJ_PERT_ROSSI_ALPHA      (idx, [1:   2]) = [ -3.56512E+00 0.00128 ];

% Inverse neutron speed :

ANA_INV_SPD               (idx, [1:   2]) = [  3.21717E-06 7.4E-05 ];

% Analog slowing-down and thermal neutron lifetime (total/prompt/delayed):

ANA_SLOW_TIME             (idx, [1:   6]) = [  1.08864E-05 0.00010  1.08860E-05 0.00010  1.09300E-05 0.00112 ];
ANA_THERM_TIME            (idx, [1:   6]) = [  5.10764E-03 0.00032  5.10973E-03 0.00032  4.84644E-03 0.00373 ];
ANA_THERM_FRAC            (idx, [1:   6]) = [  3.21917E-01 0.00015  3.22064E-01 0.00016  3.05013E-01 0.00234 ];
ANA_DELAYED_EMTIME        (idx, [1:   2]) = [  1.11429E+01 0.00252 ];
ANA_MEAN_NCOL             (idx, [1:   4]) = [  2.42024E+02 0.00030  2.56192E+02 0.00042 ];

% Group constant generation:

GC_UNIVERSE_NAME          (idx, [1:  3])  = 'UO2' ;

% Micro- and macro-group structures:

MICRO_NG                  (idx, 1)        = 70 ;
MICRO_E                   (idx, [1:  71]) = [  1.00000E-11  5.00000E-09  1.00000E-08  1.50000E-08  2.00000E-08  2.50000E-08  3.00000E-08  3.50000E-08  4.20000E-08  5.00000E-08  5.80000E-08  6.70000E-08  8.00000E-08  1.00000E-07  1.40000E-07  1.80000E-07  2.20000E-07  2.50000E-07  2.80000E-07  3.00000E-07  3.20000E-07  3.50000E-07  4.00000E-07  5.00000E-07  6.25000E-07  7.80000E-07  8.50000E-07  9.10000E-07  9.50000E-07  9.72000E-07  9.96000E-07  1.02000E-06  1.04500E-06  1.07100E-06  1.09700E-06  1.12300E-06  1.15000E-06  1.30000E-06  1.50000E-06  1.85500E-06  2.10000E-06  2.60000E-06  3.30000E-06  4.00000E-06  9.87700E-06  1.59680E-05  2.77000E-05  4.80520E-05  7.55014E-05  1.48728E-04  3.67262E-04  9.06898E-04  1.42510E-03  2.23945E-03  3.51910E-03  5.50000E-03  9.11800E-03  1.50300E-02  2.47800E-02  4.08500E-02  6.74300E-02  1.11000E-01  1.83000E-01  3.02500E-01  5.00000E-01  8.21000E-01  1.35300E+00  2.23100E+00  3.67900E+00  6.06550E+00  2.00000E+01 ];

MACRO_NG                  (idx, 1)        = 2 ;
MACRO_E                   (idx, [1:   3]) = [  1.00000E+37  6.25000E-07  0.00000E+00 ];

% Micro-group spectrum:

INF_MICRO_FLX             (idx, [1: 140]) = [  3.96655E+05 0.00093  1.80013E+06 0.00040  4.22715E+06 0.00026  4.67613E+06 0.00024  5.46400E+06 0.00025  1.22109E+07 0.00021  9.46061E+06 0.00018  1.32891E+07 0.00015  1.31620E+07 0.00018  1.21791E+07 0.00018  1.11398E+07 0.00020  9.71939E+06 0.00023  7.98239E+06 0.00025  6.29860E+06 0.00026  4.86373E+06 0.00026  3.28564E+06 0.00034  2.54820E+06 0.00030  2.00151E+06 0.00044  1.52504E+06 0.00050  2.00577E+06 0.00046  1.20517E+06 0.00058  5.64016E+05 0.00074  2.75311E+05 0.00074  2.38194E+05 0.00092  1.83100E+05 0.00096  2.02951E+05 0.00094  2.36367E+05 0.00097  7.53505E+04 0.00173  1.06576E+05 0.00137  1.13756E+05 0.00158  6.50431E+04 0.00149  1.23378E+05 0.00125  8.51131E+04 0.00147  6.37734E+04 0.00136  1.00579E+04 0.00246  9.69501E+03 0.00257  9.94213E+03 0.00250  1.03419E+04 0.00292  1.02468E+04 0.00282  1.01392E+04 0.00267  1.04698E+04 0.00279  9.84368E+03 0.00295  1.84914E+04 0.00250  2.91073E+04 0.00233  3.65091E+04 0.00185  9.18133E+04 0.00130  8.57142E+04 0.00154  7.49098E+04 0.00144  3.74678E+04 0.00135  2.17130E+04 0.00163  1.43051E+04 0.00217  1.46478E+04 0.00261  2.39317E+04 0.00159  2.90290E+04 0.00152  5.63114E+04 0.00133  1.03653E+05 0.00081  2.19550E+05 0.00067  1.88807E+05 0.00059  1.60434E+05 0.00055  1.26902E+05 0.00061  1.22572E+05 0.00059  1.27371E+05 0.00057  1.10984E+05 0.00062  7.55889E+04 0.00061  6.94303E+04 0.00083  6.08536E+04 0.00064  4.92416E+04 0.00082  3.43910E+04 0.00087  1.86241E+04 0.00122  4.32998E+03 0.00151 ];

% Integral parameters:

INF_KINF                  (idx, [1:   2]) = [  1.04769E+00 0.00012 ];

% Flux spectra in infinite geometry:

INF_FLX                   (idx, [1:   4]) = [  6.57924E+01 0.00013  9.12291E-01 0.00024 ];
INF_FISS_FLX              (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Reaction cross sections:

INF_TOT                   (idx, [1:   4]) = [  4.20034E-01 1.1E-05  6.99832E-01 4.8E-05 ];
INF_CAPT                  (idx, [1:   4]) = [  7.42036E-03 8.2E-05  7.61231E-02 1.0E-04 ];
INF_ABS                   (idx, [1:   4]) = [  1.03113E-02 7.7E-05  3.07793E-01 0.00011 ];
INF_FISS                  (idx, [1:   4]) = [  2.89091E-03 0.00011  2.31670E-01 0.00011 ];
INF_NSF                   (idx, [1:   4]) = [  7.41212E-03 0.00011  5.64394E-01 0.00011 ];
INF_NUBAR                 (idx, [1:   4]) = [  2.56394E+00 9.6E-06  2.43620E+00 0.0E+00 ];
INF_KAPPA                 (idx, [1:   4]) = [  2.03507E+02 8.7E-07  2.02270E+02 0.0E+00 ];
INF_INVV                  (idx, [1:   4]) = [  9.97366E-09 0.00038  2.82354E-06 9.6E-05 ];

% Total scattering cross sections:

INF_SCATT0                (idx, [1:   4]) = [  4.09722E-01 1.2E-05  3.91978E-01 8.9E-05 ];
INF_SCATT1                (idx, [1:   4]) = [  3.69026E-02 0.00015  7.90637E-03 0.00534 ];
INF_SCATT2                (idx, [1:   4]) = [  1.36889E-02 0.00027  2.15587E-04 0.13145 ];
INF_SCATT3                (idx, [1:   4]) = [  5.57517E-03 0.00054  4.37499E-05 0.59592 ];
INF_SCATT4                (idx, [1:   4]) = [  3.20519E-03 0.00087  2.80032E-05 0.75353 ];
INF_SCATT5                (idx, [1:   4]) = [  1.54991E-03 0.00162  2.54632E-06 1.00000 ];
INF_SCATT6                (idx, [1:   4]) = [  7.84208E-04 0.00309  6.33780E-06 1.00000 ];
INF_SCATT7                (idx, [1:   4]) = [  3.21437E-04 0.00653 -3.56515E-05 0.43298 ];

% Total scattering production cross sections:

INF_SCATTP0               (idx, [1:   4]) = [  4.09759E-01 1.2E-05  3.91978E-01 8.9E-05 ];
INF_SCATTP1               (idx, [1:   4]) = [  3.69031E-02 0.00015  7.90637E-03 0.00534 ];
INF_SCATTP2               (idx, [1:   4]) = [  1.36890E-02 0.00027  2.15587E-04 0.13145 ];
INF_SCATTP3               (idx, [1:   4]) = [  5.57517E-03 0.00054  4.37499E-05 0.59592 ];
INF_SCATTP4               (idx, [1:   4]) = [  3.20518E-03 0.00087  2.80032E-05 0.75353 ];
INF_SCATTP5               (idx, [1:   4]) = [  1.54992E-03 0.00162  2.54632E-06 1.00000 ];
INF_SCATTP6               (idx, [1:   4]) = [  7.84192E-04 0.00310  6.33780E-06 1.00000 ];
INF_SCATTP7               (idx, [1:   4]) = [  3.21449E-04 0.00654 -3.56515E-05 0.43298 ];

% Diffusion parameters:

INF_TRANSPXS              (idx, [1:   4]) = [  3.44868E-01 3.2E-05  6.55478E-01 9.6E-05 ];
INF_DIFFCOEF              (idx, [1:   4]) = [  9.66554E-01 3.2E-05  5.08535E-01 9.6E-05 ];

% Reduced absoption and removal:

INF_RABSXS                (idx, [1:   4]) = [  1.02748E-02 7.8E-05  3.07793E-01 0.00011 ];
INF_REMXS                 (idx, [1:   4]) = [  1.03860E-02 0.00011  3.08408E-01 0.00016 ];

% Poison cross sections:

INF_I135_YIELD            (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_XE135_YIELD           (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_PM147_YIELD           (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_PM148_YIELD           (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_PM148M_YIELD          (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_PM149_YIELD           (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_SM149_YIELD           (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_I135_MICRO_ABS        (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_XE135_MICRO_ABS       (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_PM147_MICRO_ABS       (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_PM148_MICRO_ABS       (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_PM148M_MICRO_ABS      (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_PM149_MICRO_ABS       (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_SM149_MICRO_ABS       (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_XE135_MACRO_ABS       (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_SM149_MACRO_ABS       (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Fission spectra:

INF_CHIT                  (idx, [1:   4]) = [  1.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_CHIP                  (idx, [1:   4]) = [  1.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_CHID                  (idx, [1:   4]) = [  1.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Scattering matrixes:

INF_S0                    (idx, [1:   8]) = [  4.09648E-01 1.2E-05  7.38218E-05 0.00153  5.54046E-04 0.00417  3.91424E-01 8.9E-05 ];
INF_S1                    (idx, [1:   8]) = [  3.69230E-02 0.00015 -2.03571E-05 0.00313  5.22838E-06 0.30784  7.90114E-03 0.00535 ];
INF_S2                    (idx, [1:   8]) = [  1.36901E-02 0.00027 -1.23244E-06 0.04590 -1.93372E-05 0.06175  2.34924E-04 0.12163 ];
INF_S3                    (idx, [1:   8]) = [  5.57528E-03 0.00054 -1.09801E-07 0.42797 -1.28674E-05 0.06971  5.66173E-05 0.46005 ];
INF_S4                    (idx, [1:   8]) = [  3.20533E-03 0.00087 -1.37582E-07 0.24350 -5.86871E-06 0.14532  3.38719E-05 0.62074 ];
INF_S5                    (idx, [1:   8]) = [  1.54999E-03 0.00162 -8.23210E-08 0.41416 -4.97894E-06 0.14406  7.52526E-06 1.00000 ];
INF_S6                    (idx, [1:   8]) = [  7.84260E-04 0.00309 -5.26271E-08 0.66456 -2.04799E-06 0.29845  8.38579E-06 1.00000 ];
INF_S7                    (idx, [1:   8]) = [  3.21477E-04 0.00654 -3.96302E-08 0.76954 -2.62717E-06 0.21702 -3.30243E-05 0.46991 ];

% Scattering production matrixes:

INF_SP0                   (idx, [1:   8]) = [  4.09685E-01 1.2E-05  7.38218E-05 0.00153  5.54046E-04 0.00417  3.91424E-01 8.9E-05 ];
INF_SP1                   (idx, [1:   8]) = [  3.69235E-02 0.00015 -2.03571E-05 0.00313  5.22838E-06 0.30784  7.90114E-03 0.00535 ];
INF_SP2                   (idx, [1:   8]) = [  1.36902E-02 0.00027 -1.23244E-06 0.04590 -1.93372E-05 0.06175  2.34924E-04 0.12163 ];
INF_SP3                   (idx, [1:   8]) = [  5.57528E-03 0.00054 -1.09801E-07 0.42797 -1.28674E-05 0.06971  5.66173E-05 0.46005 ];
INF_SP4                   (idx, [1:   8]) = [  3.20532E-03 0.00087 -1.37582E-07 0.24350 -5.86871E-06 0.14532  3.38719E-05 0.62074 ];
INF_SP5                   (idx, [1:   8]) = [  1.55000E-03 0.00162 -8.23210E-08 0.41416 -4.97894E-06 0.14406  7.52526E-06 1.00000 ];
INF_SP6                   (idx, [1:   8]) = [  7.84245E-04 0.00310 -5.26271E-08 0.66456 -2.04799E-06 0.29845  8.38579E-06 1.00000 ];
INF_SP7                   (idx, [1:   8]) = [  3.21489E-04 0.00654 -3.96302E-08 0.76954 -2.62717E-06 0.21702 -3.30243E-05 0.46991 ];

% Micro-group spectrum:

B1_MICRO_FLX              (idx, [1: 140]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Integral parameters:

B1_KINF                   (idx, [1:   2]) = [  0.00000E+00 0.0E+00 ];
B1_KEFF                   (idx, [1:   2]) = [  0.00000E+00 0.0E+00 ];
B1_B2                     (idx, [1:   2]) = [  0.00000E+00 0.0E+00 ];
B1_ERR                    (idx, [1:   2]) = [  0.00000E+00 0.0E+00 ];

% Critical spectra in infinite geometry:

B1_FLX                    (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_FISS_FLX               (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Reaction cross sections:

B1_TOT                    (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_CAPT                   (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_ABS                    (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_FISS                   (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_NSF                    (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_NUBAR                  (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_KAPPA                  (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_INVV                   (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Total scattering cross sections:

B1_SCATT0                 (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATT1                 (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATT2                 (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATT3                 (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATT4                 (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATT5                 (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATT6                 (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATT7                 (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Total scattering production cross sections:

B1_SCATTP0                (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATTP1                (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATTP2                (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATTP3                (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATTP4                (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATTP5                (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATTP6                (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATTP7                (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Diffusion parameters:

B1_TRANSPXS               (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_DIFFCOEF               (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Reduced absoption and removal:

B1_RABSXS                 (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_REMXS                  (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Poison cross sections:

B1_I135_YIELD             (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_XE135_YIELD            (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_PM147_YIELD            (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_PM148_YIELD            (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_PM148M_YIELD           (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_PM149_YIELD            (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SM149_YIELD            (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_I135_MICRO_ABS         (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_XE135_MICRO_ABS        (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_PM147_MICRO_ABS        (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_PM148_MICRO_ABS        (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_PM148M_MICRO_ABS       (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_PM149_MICRO_ABS        (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SM149_MICRO_ABS        (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_XE135_MACRO_ABS        (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SM149_MACRO_ABS        (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Fission spectra:

B1_CHIT                   (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_CHIP                   (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_CHID                   (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Scattering matrixes:

B1_S0                     (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_S1                     (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_S2                     (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_S3                     (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_S4                     (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_S5                     (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_S6                     (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_S7                     (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Scattering production matrixes:

B1_SP0                    (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SP1                    (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SP2                    (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SP3                    (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SP4                    (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SP5                    (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SP6                    (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SP7                    (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Additional diffusion parameters:

CMM_TRANSPXS              (idx, [1:   4]) = [  3.63059E-01 0.00019  2.01836E-03 0.00046 ];
CMM_TRANSPXS_X            (idx, [1:   4]) = [  4.21907E-01 0.00028  9.64074E-02 0.00137 ];
CMM_TRANSPXS_Y            (idx, [1:   4]) = [  3.39319E-01 0.00022  1.35499E-03 0.00071 ];
CMM_TRANSPXS_Z            (idx, [1:   4]) = [  3.39461E-01 0.00030  1.35511E-03 0.00061 ];
CMM_DIFFCOEF              (idx, [1:   4]) = [  9.18127E-01 0.00019  1.65152E+02 0.00046 ];
CMM_DIFFCOEF_X            (idx, [1:   4]) = [  7.90066E-01 0.00028  3.45787E+00 0.00137 ];
CMM_DIFFCOEF_Y            (idx, [1:   4]) = [  9.82363E-01 0.00022  2.46010E+02 0.00071 ];
CMM_DIFFCOEF_Z            (idx, [1:   4]) = [  9.81953E-01 0.00030  2.45988E+02 0.00061 ];

% Delayed neutron parameters (Meulekamp method):

BETA_EFF                  (idx, [1:  18]) = [  7.12372E-03 0.00181  1.77695E-04 0.01153  9.83290E-04 0.00484  5.31750E-04 0.00667  1.27579E-03 0.00431  2.24693E-03 0.00315  8.90216E-04 0.00496  6.93766E-04 0.00571  3.24278E-04 0.00817 ];
LAMBDA                    (idx, [1:  18]) = [  5.27878E-01 0.00263  1.24667E-02 7.2E-09  2.82917E-02 0.0E+00  4.25244E-02 0.0E+00  1.33042E-01 0.0E+00  2.92467E-01 0.0E+00  6.66488E-01 2.7E-09  1.63478E+00 4.2E-09  3.55460E+00 2.8E-09 ];


% Increase counter:

if (exist('idx', 'var'));
  idx = idx + 1;
else;
  idx = 1;
end;

% Version, title and date:

VERSION                   (idx, [1: 14])  = 'Serpent 2.1.30' ;
COMPILE_DATE              (idx, [1: 20])  = 'Feb 27 2018 16:26:15' ;
DEBUG                     (idx, 1)        = 0 ;
TITLE                     (idx, [1:  8])  = 'Untitled' ;
CONFIDENTIAL_DATA         (idx, 1)        = 0 ;
INPUT_FILE_NAME           (idx, [1:  8])  = 'slab_d2o' ;
WORKING_DIRECTORY         (idx, [1: 35])  = '/home/abrate/phytra_ANE/XS_2G/UO2_1' ;
HOSTNAME                  (idx, [1:  7])  = 'vpcen13' ;
CPU_TYPE                  (idx, [1: 41])  = 'Intel(R) Xeon(R) CPU E5-2630 v3 @ 2.40GHz' ;
CPU_MHZ                   (idx, 1)        = 4294967295.0 ;
START_DATE                (idx, [1: 24])  = 'Sat Sep 29 12:28:58 2018' ;
COMPLETE_DATE             (idx, [1: 24])  = 'Sat Sep 29 14:03:50 2018' ;

% Run parameters:

POP                       (idx, 1)        = 100000 ;
CYCLES                    (idx, 1)        = 1000 ;
SKIP                      (idx, 1)        = 500 ;
BATCH_INTERVAL            (idx, 1)        = 1 ;
SRC_NORM_MODE             (idx, 1)        = 2 ;
SEED                      (idx, 1)        = 1538216938 ;
UFS_MODE                  (idx, 1)        = 0 ;
UFS_ORDER                 (idx, 1)        = 1.00000;
NEUTRON_TRANSPORT_MODE    (idx, 1)        = 1 ;
PHOTON_TRANSPORT_MODE     (idx, 1)        = 0 ;
GROUP_CONSTANT_GENERATION (idx, 1)        = 1 ;
B1_CALCULATION            (idx, [1:  3])  = [ 0 0 0 ];
B1_BURNUP_CORRECTION      (idx, 1)        = 0 ;
IMPLICIT_REACTION_RATES   (idx, 1)        = 1 ;

% Optimization:

OPTIMIZATION_MODE         (idx, 1)        = 4 ;
RECONSTRUCT_MICROXS       (idx, 1)        = 1 ;
RECONSTRUCT_MACROXS       (idx, 1)        = 1 ;
DOUBLE_INDEXING           (idx, 1)        = 0 ;
MG_MAJORANT_MODE          (idx, 1)        = 0 ;

% Parallelization:

MPI_TASKS                 (idx, 1)        = 1 ;
OMP_THREADS               (idx, 1)        = 30 ;
MPI_REPRODUCIBILITY       (idx, 1)        = 0 ;
OMP_REPRODUCIBILITY       (idx, 1)        = 1 ;
OMP_HISTORY_PROFILE       (idx, [1:  30]) = [  1.08888E+00  9.94430E-01  1.02065E+00  1.01759E+00  1.01579E+00  9.71243E-01  9.57860E-01  1.03162E+00  1.05536E+00  9.61553E-01  9.87352E-01  9.60729E-01  9.74345E-01  1.03210E+00  9.67696E-01  1.01849E+00  9.58672E-01  1.02126E+00  9.57751E-01  1.01910E+00  9.97507E-01  1.02151E+00  1.01712E+00  1.02707E+00  9.57848E-01  1.01812E+00  9.61310E-01  1.03389E+00  9.55931E-01  9.97215E-01  ];
SHARE_BUF_ARRAY           (idx, 1)        = 0 ;
SHARE_RES2_ARRAY          (idx, 1)        = 1 ;

% File paths:

XS_DATA_FILE_PATH         (idx, [1: 46])  = '/opt/serpent/xsdata/jeff311/sss_jeff311.xsdata' ;
DECAY_DATA_FILE_PATH      (idx, [1:  3])  = 'N/A' ;
SFY_DATA_FILE_PATH        (idx, [1:  3])  = 'N/A' ;
NFY_DATA_FILE_PATH        (idx, [1:  3])  = 'N/A' ;
BRA_DATA_FILE_PATH        (idx, [1:  3])  = 'N/A' ;

% Collision and reaction sampling (neutrons/photons):

MIN_MACROXS               (idx, [1:   4]) = [  5.00000E-02 0.0E+00  0.00000E+00 0.0E+00 ];
DT_THRESH                 (idx, [1:  2])  = [  9.00000E-01  9.00000E-01 ];
ST_FRAC                   (idx, [1:   4]) = [  4.45997E-04 0.00038  0.00000E+00 0.0E+00 ];
DT_FRAC                   (idx, [1:   4]) = [  9.99554E-01 1.7E-07  0.00000E+00 0.0E+00 ];
DT_EFF                    (idx, [1:   4]) = [  6.17261E-01 2.1E-05  0.00000E+00 0.0E+00 ];
REA_SAMPLING_EFF          (idx, [1:   4]) = [  1.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
REA_SAMPLING_FAIL         (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
TOT_COL_EFF               (idx, [1:   4]) = [  6.17425E-01 2.1E-05  0.00000E+00 0.0E+00 ];
AVG_TRACKING_LOOPS        (idx, [1:   8]) = [  2.18175E+00 3.9E-05  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
AVG_TRACKS                (idx, [1:   4]) = [  2.42050E+02 0.00030  0.00000E+00 0.0E+00 ];
AVG_REAL_COL              (idx, [1:   4]) = [  2.42024E+02 0.00030  0.00000E+00 0.0E+00 ];
AVG_VIRT_COL              (idx, [1:   4]) = [  1.49968E+02 0.00036  0.00000E+00 0.0E+00 ];
AVG_SURF_CROSS            (idx, [1:   4]) = [  3.27516E-02 0.00053  0.00000E+00 0.0E+00 ];
LOST_PARTICLES            (idx, 1)        = 0 ;

% Run statistics:

CYCLE_IDX                 (idx, 1)        = 1000 ;
SOURCE_POPULATION         (idx, 1)        = 100000976 ;
MEAN_POP_SIZE             (idx, [1:  2])  = [  1.00001E+05 0.00019 ];
MEAN_POP_WGT              (idx, [1:  2])  = [  1.00001E+05 0.00019 ];
SIMULATION_COMPLETED      (idx, 1)        = 1 ;

% Running times:

TOT_CPU_TIME              (idx, 1)        =  2.70421E+03 ;
RUNNING_TIME              (idx, 1)        =  9.48730E+01 ;
INIT_TIME                 (idx, [1:  2])  = [  3.18833E-02  3.18833E-02 ];
PROCESS_TIME              (idx, [1:  2])  = [  2.66668E-04  2.66668E-04 ];
TRANSPORT_CYCLE_TIME      (idx, [1:  3])  = [  9.48409E+01  9.48409E+01  0.00000E+00 ];
MPI_OVERHEAD_TIME         (idx, [1:  2])  = [  0.00000E+00  0.00000E+00 ];
ESTIMATED_RUNNING_TIME    (idx, [1:  2])  = [  9.48715E+01  0.00000E+00 ];
CPU_USAGE                 (idx, 1)        = 28.50343 ;
TRANSPORT_CPU_USAGE       (idx, [1:   2]) = [  2.87208E+01 7.9E-05 ];
OMP_PARALLEL_FRAC         (idx, 1)        =  9.50408E-01 ;

% Memory usage:

AVAIL_MEM                 (idx, 1)        = 96872.39 ;
ALLOC_MEMSIZE             (idx, 1)        = 1036.29;
MEMSIZE                   (idx, 1)        = 752.09;
XS_MEMSIZE                (idx, 1)        = 86.48;
MAT_MEMSIZE               (idx, 1)        = 9.01;
RES_MEMSIZE               (idx, 1)        = 3.37;
MISC_MEMSIZE              (idx, 1)        = 653.23;
UNKNOWN_MEMSIZE           (idx, 1)        = 0.00;
UNUSED_MEMSIZE            (idx, 1)        = 284.20;

% Geometry parameters:

TOT_CELLS                 (idx, 1)        = 7 ;
UNION_CELLS               (idx, 1)        = 0 ;

% Neutron energy grid:

NEUTRON_ERG_TOL           (idx, 1)        =  0.00000E+00 ;
NEUTRON_ERG_NE            (idx, 1)        = 97927 ;
NEUTRON_EMIN              (idx, 1)        =  1.00000E-11 ;
NEUTRON_EMAX              (idx, 1)        =  2.00000E+01 ;

% Unresolved resonance probability table sampling:

URES_DILU_CUT             (idx, 1)        =  1.00000E-09 ;
URES_EMIN                 (idx, 1)        =  1.00000E+37 ;
URES_EMAX                 (idx, 1)        = -1.00000E+37 ;
URES_AVAIL                (idx, 1)        = 2 ;
URES_USED                 (idx, 1)        = 0 ;

% Nuclides and reaction channels:

TOT_NUCLIDES              (idx, 1)        = 5 ;
TOT_TRANSPORT_NUCLIDES    (idx, 1)        = 5 ;
TOT_DOSIMETRY_NUCLIDES    (idx, 1)        = 0 ;
TOT_DECAY_NUCLIDES        (idx, 1)        = 0 ;
TOT_PHOTON_NUCLIDES       (idx, 1)        = 0 ;
TOT_REA_CHANNELS          (idx, 1)        = 114 ;
TOT_TRANSMU_REA           (idx, 1)        = 0 ;

% Neutron physics options:

USE_DELNU                 (idx, 1)        = 1 ;
USE_URES                  (idx, 1)        = 0 ;
USE_DBRC                  (idx, 1)        = 0 ;
IMPL_CAPT                 (idx, 1)        = 0 ;
IMPL_NXN                  (idx, 1)        = 1 ;
IMPL_FISS                 (idx, 1)        = 0 ;
DOPPLER_PREPROCESSOR      (idx, 1)        = 0 ;
TMS_MODE                  (idx, 1)        = 0 ;
SAMPLE_FISS               (idx, 1)        = 1 ;
SAMPLE_CAPT               (idx, 1)        = 1 ;
SAMPLE_SCATT              (idx, 1)        = 1 ;

% Radioactivity data:

TOT_ACTIVITY              (idx, 1)        =  0.00000E+00 ;
TOT_DECAY_HEAT            (idx, 1)        =  0.00000E+00 ;
TOT_SF_RATE               (idx, 1)        =  0.00000E+00 ;
ACTINIDE_ACTIVITY         (idx, 1)        =  0.00000E+00 ;
ACTINIDE_DECAY_HEAT       (idx, 1)        =  0.00000E+00 ;
FISSION_PRODUCT_ACTIVITY  (idx, 1)        =  0.00000E+00 ;
FISSION_PRODUCT_DECAY_HEAT(idx, 1)        =  0.00000E+00 ;
INHALATION_TOXICITY       (idx, 1)        =  0.00000E+00 ;
INGESTION_TOXICITY        (idx, 1)        =  0.00000E+00 ;
ACTINIDE_INH_TOX          (idx, 1)        =  0.00000E+00 ;
ACTINIDE_ING_TOX          (idx, 1)        =  0.00000E+00 ;
FISSION_PRODUCT_INH_TOX   (idx, 1)        =  0.00000E+00 ;
FISSION_PRODUCT_ING_TOX   (idx, 1)        =  0.00000E+00 ;
SR90_ACTIVITY             (idx, 1)        =  0.00000E+00 ;
TE132_ACTIVITY            (idx, 1)        =  0.00000E+00 ;
I131_ACTIVITY             (idx, 1)        =  0.00000E+00 ;
I132_ACTIVITY             (idx, 1)        =  0.00000E+00 ;
CS134_ACTIVITY            (idx, 1)        =  0.00000E+00 ;
CS137_ACTIVITY            (idx, 1)        =  0.00000E+00 ;
PHOTON_DECAY_SOURCE       (idx, 1)        =  0.00000E+00 ;
NEUTRON_DECAY_SOURCE      (idx, 1)        =  0.00000E+00 ;
ALPHA_DECAY_SOURCE        (idx, 1)        =  0.00000E+00 ;
ELECTRON_DECAY_SOURCE     (idx, 1)        =  0.00000E+00 ;

% Normalization coefficient:

NORM_COEF                 (idx, [1:   4]) = [  9.96625E-06 8.7E-05  0.00000E+00 0.0E+00 ];

% Analog reaction rate estimators:

CONVERSION_RATIO          (idx, [1:   2]) = [  1.12808E+00 0.00022 ];
U235_FISS                 (idx, [1:   4]) = [  3.38897E-01 0.00016  8.43815E-01 6.7E-05 ];
U238_FISS                 (idx, [1:   4]) = [  6.27283E-02 0.00041  1.56185E-01 0.00036 ];
U235_CAPT                 (idx, [1:   4]) = [  8.08771E-02 0.00036  1.41112E-01 0.00034 ];
U238_CAPT                 (idx, [1:   4]) = [  4.73595E-01 0.00015  8.26313E-01 6.5E-05 ];

% Neutron balance (particles/weight):

BALA_SRC_NEUTRON_SRC     (idx, [1:  2])  = [ 0 0.00000E+00 ];
BALA_SRC_NEUTRON_FISS    (idx, [1:  2])  = [ 100000976 1.00000E+08 ];
BALA_SRC_NEUTRON_NXN     (idx, [1:  2])  = [ 0 3.51983E+05 ];
BALA_SRC_NEUTRON_VR      (idx, [1:  2])  = [ 0 0.00000E+00 ];
BALA_SRC_NEUTRON_TOT     (idx, [1:  2])  = [ 100000976 1.00352E+08 ];

BALA_LOSS_NEUTRON_CAPT    (idx, [1:  2])  = [ 57302596 5.75083E+07 ];
BALA_LOSS_NEUTRON_FISS    (idx, [1:  2])  = [ 40165004 4.02985E+07 ];
BALA_LOSS_NEUTRON_LEAK    (idx, [1:  2])  = [ 2533376 2.54514E+06 ];
BALA_LOSS_NEUTRON_CUT     (idx, [1:  2])  = [ 0 0.00000E+00 ];
BALA_LOSS_NEUTRON_ERR     (idx, [1:  2])  = [ 0 0.00000E+00 ];
BALA_LOSS_NEUTRON_TOT     (idx, [1:  2])  = [ 100000976 1.00352E+08 ];

BALA_NEUTRON_DIFF         (idx, [1:  2])  = [ 0 -1.02967E-05 ];

% Normalized total reaction rates (neutrons):

TOT_POWER                 (idx, [1:   2]) = [  1.30508E-11 8.0E-05 ];
TOT_POWDENS               (idx, [1:   2]) = [  0.00000E+00 0.0E+00 ];
TOT_GENRATE               (idx, [1:   2]) = [  1.00255E+00 7.9E-05 ];
TOT_FISSRATE              (idx, [1:   2]) = [  4.01548E-01 8.0E-05 ];
TOT_CAPTRATE              (idx, [1:   2]) = [  5.73087E-01 5.8E-05 ];
TOT_ABSRATE               (idx, [1:   2]) = [  9.74635E-01 1.5E-05 ];
TOT_SRCRATE               (idx, [1:   2]) = [  9.96625E-01 8.7E-05 ];
TOT_FLUX                  (idx, [1:   2]) = [  5.10661E+02 0.00029 ];
TOT_PHOTON_PRODRATE       (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
TOT_LEAKRATE              (idx, [1:   2]) = [  2.53654E-02 0.00059 ];
ALBEDO_LEAKRATE           (idx, [1:   2]) = [  0.00000E+00 0.0E+00 ];
TOT_LOSSRATE              (idx, [1:   2]) = [  1.00000E+00 0.0E+00 ];
TOT_CUTRATE               (idx, [1:   2]) = [  0.00000E+00 0.0E+00 ];
TOT_RR                    (idx, [1:   2]) = [  2.42265E+02 0.00030 ];
INI_FMASS                 (idx, 1)        =  0.00000E+00 ;
TOT_FMASS                 (idx, 1)        =  0.00000E+00 ;

% Six-factor formula:

SIX_FF_ETA                (idx, [1:   2]) = [  1.83338E+00 0.00010 ];
SIX_FF_F                  (idx, [1:   2]) = [  9.50966E-01 4.1E-05 ];
SIX_FF_P                  (idx, [1:   2]) = [  3.04073E-01 0.00016 ];
SIX_FF_EPSILON            (idx, [1:   2]) = [  1.94747E+00 0.00016 ];
SIX_FF_LF                 (idx, [1:   2]) = [  1.00000E+00 0.0E+00 ];
SIX_FF_LT                 (idx, [1:   2]) = [  9.74549E-01 1.5E-05 ];
SIX_FF_KINF               (idx, [1:   2]) = [  1.03242E+00 0.00013 ];
SIX_FF_KEFF               (idx, [1:   2]) = [  1.00614E+00 0.00013 ];

% Fission neutron and energy production:

NUBAR                     (idx, [1:   2]) = [  2.49670E+00 5.9E-06 ];
FISSE                     (idx, [1:   2]) = [  2.02856E+02 5.4E-07 ];

% Criticality eigenvalues:

ANA_KEFF                  (idx, [1:   6]) = [  1.00619E+00 0.00013  9.99121E-01 0.00013  7.02222E-03 0.00182 ];
IMP_KEFF                  (idx, [1:   2]) = [  1.00608E+00 8.0E-05 ];
COL_KEFF                  (idx, [1:   2]) = [  1.00595E+00 0.00011 ];
ABS_KEFF                  (idx, [1:   2]) = [  1.00608E+00 8.0E-05 ];
ABS_KINF                  (idx, [1:   2]) = [  1.03236E+00 7.8E-05 ];
GEOM_ALBEDO               (idx, [1:   6]) = [  1.00000E+00 0.0E+00  1.00000E+00 0.0E+00  1.00000E+00 0.0E+00 ];

% ALF (Average lethargy of neutrons causing fission):
% Based on E0 = 2.000000E+01 MeV

ANA_ALF                   (idx, [1:   2]) = [  1.34865E+01 8.9E-05 ];
IMP_ALF                   (idx, [1:   2]) = [  1.34871E+01 6.8E-05 ];

% EALF (Energy corresponding to average lethargy of neutrons causing fission):

ANA_EALF                  (idx, [1:   2]) = [  2.78130E-05 0.00119 ];
IMP_EALF                  (idx, [1:   2]) = [  2.77881E-05 0.00092 ];

% AFGE (Average energy of neutrons causing fission):

ANA_AFGE                  (idx, [1:   2]) = [  5.25106E-01 0.00038 ];
IMP_AFGE                  (idx, [1:   2]) = [  5.24974E-01 0.00022 ];

% Forward-weighted delayed neutron parameters:

FWD_ANA_BETA_ZERO         (idx, [1:  18]) = [  8.37638E-03 0.00117  2.05386E-04 0.00713  1.12987E-03 0.00301  6.23115E-04 0.00398  1.48908E-03 0.00268  2.62622E-03 0.00198  1.07392E-03 0.00304  8.32284E-04 0.00353  3.96514E-04 0.00490 ];
FWD_ANA_LAMBDA            (idx, [1:  18]) = [  5.38799E-01 0.00161  1.24667E-02 7.2E-09  2.82917E-02 0.0E+00  4.25244E-02 0.0E+00  1.33042E-01 0.0E+00  2.92467E-01 0.0E+00  6.66488E-01 2.7E-09  1.63478E+00 4.2E-09  3.55460E+00 2.9E-09 ];

% Beta-eff using Meulekamp's method:

ADJ_MEULEKAMP_BETA_EFF    (idx, [1:  18]) = [  7.12372E-03 0.00181  1.77695E-04 0.01153  9.83290E-04 0.00484  5.31750E-04 0.00667  1.27579E-03 0.00431  2.24693E-03 0.00315  8.90216E-04 0.00496  6.93766E-04 0.00571  3.24278E-04 0.00817 ];
ADJ_MEULEKAMP_LAMBDA      (idx, [1:  18]) = [  5.27878E-01 0.00263  1.24667E-02 7.2E-09  2.82917E-02 0.0E+00  4.25244E-02 0.0E+00  1.33042E-01 0.0E+00  2.92467E-01 0.0E+00  6.66488E-01 2.7E-09  1.63478E+00 4.2E-09  3.55460E+00 2.8E-09 ];

% Adjoint weighted time constants using Nauchi's method:

ADJ_NAUCHI_GEN_TIME       (idx, [1:   6]) = [  1.90111E-03 0.00063  1.90002E-03 0.00064  2.05477E-03 0.00672 ];
ADJ_NAUCHI_LIFETIME       (idx, [1:   6]) = [  1.91284E-03 0.00061  1.91175E-03 0.00062  2.06739E-03 0.00671 ];
ADJ_NAUCHI_BETA_EFF       (idx, [1:  18]) = [  6.98322E-03 0.00186  1.75027E-04 0.01249  9.70819E-04 0.00514  5.16308E-04 0.00674  1.24334E-03 0.00450  2.20892E-03 0.00334  8.68910E-04 0.00528  6.81497E-04 0.00616  3.18387E-04 0.00879 ];
ADJ_NAUCHI_LAMBDA         (idx, [1:  18]) = [  5.28224E-01 0.00282  1.24667E-02 7.1E-09  2.82917E-02 0.0E+00  4.25244E-02 0.0E+00  1.33042E-01 0.0E+00  2.92467E-01 0.0E+00  6.66488E-01 2.7E-09  1.63478E+00 4.2E-09  3.55460E+00 2.9E-09 ];

% Adjoint weighted time constants using IFP:

ADJ_IFP_GEN_TIME          (idx, [1:   6]) = [  2.04747E-03 0.00151  2.04597E-03 0.00151  2.25135E-03 0.01799 ];
ADJ_IFP_LIFETIME          (idx, [1:   6]) = [  2.06010E-03 0.00150  2.05860E-03 0.00150  2.26525E-03 0.01799 ];
ADJ_IFP_IMP_BETA_EFF      (idx, [1:  18]) = [  7.06305E-03 0.00648  1.78447E-04 0.03926  9.69265E-04 0.01762  5.34667E-04 0.02252  1.23453E-03 0.01525  2.25953E-03 0.01137  8.98532E-04 0.01748  6.83517E-04 0.02047  3.04563E-04 0.02884 ];
ADJ_IFP_IMP_LAMBDA        (idx, [1:  18]) = [  5.20992E-01 0.00941  1.24667E-02 0.0E+00  2.82917E-02 0.0E+00  4.25244E-02 0.0E+00  1.33042E-01 0.0E+00  2.92467E-01 0.0E+00  6.66488E-01 2.7E-09  1.63478E+00 4.4E-09  3.55460E+00 2.7E-09 ];
ADJ_IFP_ANA_BETA_EFF      (idx, [1:  18]) = [  7.06705E-03 0.00637  1.76679E-04 0.03824  9.71410E-04 0.01723  5.35958E-04 0.02184  1.24027E-03 0.01514  2.25908E-03 0.01108  8.96590E-04 0.01712  6.82861E-04 0.02004  3.04202E-04 0.02780 ];
ADJ_IFP_ANA_LAMBDA        (idx, [1:  18]) = [  5.20340E-01 0.00913  1.24667E-02 0.0E+00  2.82917E-02 0.0E+00  4.25244E-02 0.0E+00  1.33042E-01 0.0E+00  2.92467E-01 0.0E+00  6.66488E-01 2.7E-09  1.63478E+00 4.4E-09  3.55460E+00 2.8E-09 ];
ADJ_IFP_ROSSI_ALPHA       (idx, [1:   2]) = [ -3.45836E+00 0.00659 ];

% Adjoint weighted time constants using perturbation technique:

ADJ_PERT_GEN_TIME         (idx, [1:   2]) = [  1.98808E-03 0.00032 ];
ADJ_PERT_LIFETIME         (idx, [1:   2]) = [  2.00035E-03 0.00029 ];
ADJ_PERT_BETA_EFF         (idx, [1:   2]) = [  7.08691E-03 0.00123 ];
ADJ_PERT_ROSSI_ALPHA      (idx, [1:   2]) = [ -3.56512E+00 0.00128 ];

% Inverse neutron speed :

ANA_INV_SPD               (idx, [1:   2]) = [  3.21717E-06 7.4E-05 ];

% Analog slowing-down and thermal neutron lifetime (total/prompt/delayed):

ANA_SLOW_TIME             (idx, [1:   6]) = [  1.08864E-05 0.00010  1.08860E-05 0.00010  1.09300E-05 0.00112 ];
ANA_THERM_TIME            (idx, [1:   6]) = [  5.10764E-03 0.00032  5.10973E-03 0.00032  4.84644E-03 0.00373 ];
ANA_THERM_FRAC            (idx, [1:   6]) = [  3.21917E-01 0.00015  3.22064E-01 0.00016  3.05013E-01 0.00234 ];
ANA_DELAYED_EMTIME        (idx, [1:   2]) = [  1.11429E+01 0.00252 ];
ANA_MEAN_NCOL             (idx, [1:   4]) = [  2.42024E+02 0.00030  2.56192E+02 0.00042 ];

% Group constant generation:

GC_UNIVERSE_NAME          (idx, [1:  3])  = 'D2O' ;

% Micro- and macro-group structures:

MICRO_NG                  (idx, 1)        = 70 ;
MICRO_E                   (idx, [1:  71]) = [  1.00000E-11  5.00000E-09  1.00000E-08  1.50000E-08  2.00000E-08  2.50000E-08  3.00000E-08  3.50000E-08  4.20000E-08  5.00000E-08  5.80000E-08  6.70000E-08  8.00000E-08  1.00000E-07  1.40000E-07  1.80000E-07  2.20000E-07  2.50000E-07  2.80000E-07  3.00000E-07  3.20000E-07  3.50000E-07  4.00000E-07  5.00000E-07  6.25000E-07  7.80000E-07  8.50000E-07  9.10000E-07  9.50000E-07  9.72000E-07  9.96000E-07  1.02000E-06  1.04500E-06  1.07100E-06  1.09700E-06  1.12300E-06  1.15000E-06  1.30000E-06  1.50000E-06  1.85500E-06  2.10000E-06  2.60000E-06  3.30000E-06  4.00000E-06  9.87700E-06  1.59680E-05  2.77000E-05  4.80520E-05  7.55014E-05  1.48728E-04  3.67262E-04  9.06898E-04  1.42510E-03  2.23945E-03  3.51910E-03  5.50000E-03  9.11800E-03  1.50300E-02  2.47800E-02  4.08500E-02  6.74300E-02  1.11000E-01  1.83000E-01  3.02500E-01  5.00000E-01  8.21000E-01  1.35300E+00  2.23100E+00  3.67900E+00  6.06550E+00  2.00000E+01 ];

MACRO_NG                  (idx, 1)        = 2 ;
MACRO_E                   (idx, [1:   3]) = [  1.00000E+37  6.25000E-07  0.00000E+00 ];

% Micro-group spectrum:

INF_MICRO_FLX             (idx, [1: 140]) = [  1.61359E+05 0.00161  6.82505E+05 0.00077  1.38241E+06 0.00041  1.43199E+06 0.00037  1.29664E+06 0.00036  1.94599E+06 0.00027  1.69187E+06 0.00027  2.27126E+06 0.00022  2.36481E+06 0.00021  2.44429E+06 0.00019  2.50735E+06 0.00028  2.53864E+06 0.00022  2.53316E+06 0.00028  2.50841E+06 0.00024  2.50080E+06 0.00023  2.17617E+06 0.00022  2.16800E+06 0.00021  2.13542E+06 0.00023  2.09850E+06 0.00026  4.08308E+06 0.00023  3.94011E+06 0.00018  2.85466E+06 0.00027  1.85450E+06 0.00031  2.20899E+06 0.00025  2.14965E+06 0.00025  1.84257E+06 0.00031  3.37953E+06 0.00025  7.15970E+05 0.00033  8.95841E+05 0.00034  8.08851E+05 0.00033  4.73733E+05 0.00040  8.21799E+05 0.00035  5.60112E+05 0.00041  4.82211E+05 0.00048  9.36607E+04 0.00099  9.24695E+04 0.00075  9.50981E+04 0.00077  9.78000E+04 0.00080  9.62738E+04 0.00092  9.53880E+04 0.00091  9.75910E+04 0.00078  9.20824E+04 0.00084  1.73415E+05 0.00066  2.76923E+05 0.00054  3.52576E+05 0.00055  9.26476E+05 0.00034  9.65016E+05 0.00036  1.01369E+06 0.00038  6.41945E+05 0.00038  4.59062E+05 0.00046  3.55873E+05 0.00051  4.21654E+05 0.00045  8.71068E+05 0.00028  1.54013E+06 0.00030  5.26387E+06 0.00030  1.74996E+07 0.00033  5.90262E+07 0.00033  6.71909E+07 0.00034  6.80879E+07 0.00035  6.03573E+07 0.00035  6.52596E+07 0.00034  7.58992E+07 0.00033  7.45349E+07 0.00035  5.69696E+07 0.00035  5.84157E+07 0.00034  5.82004E+07 0.00035  5.51835E+07 0.00035  4.75190E+07 0.00035  3.49888E+07 0.00035  1.38582E+07 0.00034 ];

% Integral parameters:

INF_KINF                  (idx, [1:   2]) = [  0.00000E+00 0.0E+00 ];

% Flux spectra in infinite geometry:

INF_FLX                   (idx, [1:   4]) = [  3.30884E+01 0.00018  4.10870E+02 0.00035 ];
INF_FISS_FLX              (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Reaction cross sections:

INF_TOT                   (idx, [1:   4]) = [  3.42622E-01 8.2E-06  4.93236E-01 2.9E-06 ];
INF_CAPT                  (idx, [1:   4]) = [  2.92027E-05 0.00096  3.52315E-05 9.3E-06 ];
INF_ABS                   (idx, [1:   4]) = [  2.92027E-05 0.00096  3.52315E-05 9.3E-06 ];
INF_FISS                  (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_NSF                   (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_NUBAR                 (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_KAPPA                 (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_INVV                  (idx, [1:   4]) = [  1.08969E-07 8.4E-05  3.98200E-06 9.3E-06 ];

% Total scattering cross sections:

INF_SCATT0                (idx, [1:   4]) = [  3.42593E-01 8.2E-06  4.93201E-01 2.9E-06 ];
INF_SCATT1                (idx, [1:   4]) = [  7.31567E-02 8.1E-05  5.30536E-02 3.8E-05 ];
INF_SCATT2                (idx, [1:   4]) = [  1.15769E-02 0.00036  2.10543E-04 0.00827 ];
INF_SCATT3                (idx, [1:   4]) = [  8.00734E-04 0.00454 -5.51418E-03 0.00020 ];
INF_SCATT4                (idx, [1:   4]) = [  6.59315E-05 0.05721 -1.10419E-02 0.00011 ];
INF_SCATT5                (idx, [1:   4]) = [  7.40351E-05 0.04493 -5.76834E-03 0.00018 ];
INF_SCATT6                (idx, [1:   4]) = [ -1.01582E-04 0.02481 -1.15913E-02 7.1E-05 ];
INF_SCATT7                (idx, [1:   4]) = [  5.67691E-06 0.40211 -2.89266E-03 0.00028 ];

% Total scattering production cross sections:

INF_SCATTP0               (idx, [1:   4]) = [  3.42627E-01 8.2E-06  4.93201E-01 2.9E-06 ];
INF_SCATTP1               (idx, [1:   4]) = [  7.31816E-02 8.1E-05  5.30536E-02 3.8E-05 ];
INF_SCATTP2               (idx, [1:   4]) = [  1.15914E-02 0.00036  2.10543E-04 0.00827 ];
INF_SCATTP3               (idx, [1:   4]) = [  8.07673E-04 0.00451 -5.51418E-03 0.00020 ];
INF_SCATTP4               (idx, [1:   4]) = [  6.87174E-05 0.05475 -1.10419E-02 0.00011 ];
INF_SCATTP5               (idx, [1:   4]) = [  7.50120E-05 0.04428 -5.76834E-03 0.00018 ];
INF_SCATTP6               (idx, [1:   4]) = [ -1.01258E-04 0.02490 -1.15913E-02 7.1E-05 ];
INF_SCATTP7               (idx, [1:   4]) = [  5.79902E-06 0.39420 -2.89266E-03 0.00028 ];

% Diffusion parameters:

INF_TRANSPXS              (idx, [1:   4]) = [  2.62743E-01 3.7E-05  4.25465E-01 6.1E-06 ];
INF_DIFFCOEF              (idx, [1:   4]) = [  1.26867E+00 3.7E-05  7.83457E-01 6.1E-06 ];

% Reduced absoption and removal:

INF_RABSXS                (idx, [1:   4]) = [ -4.29615E-06 0.02656  3.52315E-05 9.3E-06 ];
INF_REMXS                 (idx, [1:   4]) = [  9.68352E-03 7.8E-05  4.27776E-05 0.00073 ];

% Poison cross sections:

INF_I135_YIELD            (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_XE135_YIELD           (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_PM147_YIELD           (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_PM148_YIELD           (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_PM148M_YIELD          (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_PM149_YIELD           (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_SM149_YIELD           (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_I135_MICRO_ABS        (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_XE135_MICRO_ABS       (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_PM147_MICRO_ABS       (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_PM148_MICRO_ABS       (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_PM148M_MICRO_ABS      (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_PM149_MICRO_ABS       (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_SM149_MICRO_ABS       (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_XE135_MACRO_ABS       (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_SM149_MACRO_ABS       (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Fission spectra:

INF_CHIT                  (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_CHIP                  (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_CHID                  (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Scattering matrixes:

INF_S0                    (idx, [1:   8]) = [  3.32939E-01 8.0E-06  9.65420E-03 7.7E-05  7.53196E-06 0.00170  4.93193E-01 2.9E-06 ];
INF_S1                    (idx, [1:   8]) = [  7.39639E-02 7.8E-05 -8.07171E-04 0.00102  3.08064E-06 0.00299  5.30505E-02 3.8E-05 ];
INF_S2                    (idx, [1:   8]) = [  1.24781E-02 0.00032 -9.01218E-04 0.00069  1.01849E-06 0.00619  2.09525E-04 0.00831 ];
INF_S3                    (idx, [1:   8]) = [  1.04583E-03 0.00341 -2.45097E-04 0.00256  1.58681E-07 0.03375 -5.51434E-03 0.00020 ];
INF_S4                    (idx, [1:   8]) = [  1.73286E-04 0.02207 -1.07354E-04 0.00487 -1.69982E-07 0.02174 -1.10417E-02 0.00011 ];
INF_S5                    (idx, [1:   8]) = [  7.52821E-05 0.04300 -1.24704E-06 0.40517 -2.66125E-07 0.01466 -5.76808E-03 0.00018 ];
INF_S6                    (idx, [1:   8]) = [ -1.84820E-05 0.13144 -8.31003E-05 0.00547 -2.56468E-07 0.01395 -1.15911E-02 7.1E-05 ];
INF_S7                    (idx, [1:   8]) = [ -3.01670E-05 0.07213  3.58439E-05 0.01097 -1.86156E-07 0.01868 -2.89247E-03 0.00028 ];

% Scattering production matrixes:

INF_SP0                   (idx, [1:   8]) = [  3.32972E-01 8.0E-06  9.65420E-03 7.7E-05  7.53196E-06 0.00170  4.93193E-01 2.9E-06 ];
INF_SP1                   (idx, [1:   8]) = [  7.39887E-02 7.8E-05 -8.07171E-04 0.00102  3.08064E-06 0.00299  5.30505E-02 3.8E-05 ];
INF_SP2                   (idx, [1:   8]) = [  1.24927E-02 0.00032 -9.01218E-04 0.00069  1.01849E-06 0.00619  2.09525E-04 0.00831 ];
INF_SP3                   (idx, [1:   8]) = [  1.05277E-03 0.00339 -2.45097E-04 0.00256  1.58681E-07 0.03375 -5.51434E-03 0.00020 ];
INF_SP4                   (idx, [1:   8]) = [  1.76071E-04 0.02167 -1.07354E-04 0.00487 -1.69982E-07 0.02174 -1.10417E-02 0.00011 ];
INF_SP5                   (idx, [1:   8]) = [  7.62591E-05 0.04240 -1.24704E-06 0.40517 -2.66125E-07 0.01466 -5.76808E-03 0.00018 ];
INF_SP6                   (idx, [1:   8]) = [ -1.81578E-05 0.13383 -8.31003E-05 0.00547 -2.56468E-07 0.01395 -1.15911E-02 7.1E-05 ];
INF_SP7                   (idx, [1:   8]) = [ -3.00449E-05 0.07255  3.58439E-05 0.01097 -1.86156E-07 0.01868 -2.89247E-03 0.00028 ];

% Micro-group spectrum:

B1_MICRO_FLX              (idx, [1: 140]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Integral parameters:

B1_KINF                   (idx, [1:   2]) = [  0.00000E+00 0.0E+00 ];
B1_KEFF                   (idx, [1:   2]) = [  0.00000E+00 0.0E+00 ];
B1_B2                     (idx, [1:   2]) = [  0.00000E+00 0.0E+00 ];
B1_ERR                    (idx, [1:   2]) = [  0.00000E+00 0.0E+00 ];

% Critical spectra in infinite geometry:

B1_FLX                    (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_FISS_FLX               (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Reaction cross sections:

B1_TOT                    (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_CAPT                   (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_ABS                    (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_FISS                   (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_NSF                    (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_NUBAR                  (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_KAPPA                  (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_INVV                   (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Total scattering cross sections:

B1_SCATT0                 (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATT1                 (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATT2                 (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATT3                 (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATT4                 (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATT5                 (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATT6                 (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATT7                 (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Total scattering production cross sections:

B1_SCATTP0                (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATTP1                (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATTP2                (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATTP3                (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATTP4                (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATTP5                (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATTP6                (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATTP7                (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Diffusion parameters:

B1_TRANSPXS               (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_DIFFCOEF               (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Reduced absoption and removal:

B1_RABSXS                 (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_REMXS                  (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Poison cross sections:

B1_I135_YIELD             (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_XE135_YIELD            (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_PM147_YIELD            (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_PM148_YIELD            (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_PM148M_YIELD           (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_PM149_YIELD            (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SM149_YIELD            (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_I135_MICRO_ABS         (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_XE135_MICRO_ABS        (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_PM147_MICRO_ABS        (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_PM148_MICRO_ABS        (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_PM148M_MICRO_ABS       (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_PM149_MICRO_ABS        (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SM149_MICRO_ABS        (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_XE135_MACRO_ABS        (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SM149_MACRO_ABS        (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Fission spectra:

B1_CHIT                   (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_CHIP                   (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_CHID                   (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Scattering matrixes:

B1_S0                     (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_S1                     (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_S2                     (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_S3                     (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_S4                     (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_S5                     (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_S6                     (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_S7                     (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Scattering production matrixes:

B1_SP0                    (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SP1                    (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SP2                    (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SP3                    (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SP4                    (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SP5                    (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SP6                    (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SP7                    (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Additional diffusion parameters:

CMM_TRANSPXS              (idx, [1:   4]) = [  3.21400E-01 0.00014 -4.50212E+00 0.00034 ];
CMM_TRANSPXS_X            (idx, [1:   4]) = [  3.14757E-01 0.00017 -4.07726E+00 0.00033 ];
CMM_TRANSPXS_Y            (idx, [1:   4]) = [  3.24785E-01 0.00018 -4.74897E+00 0.00037 ];
CMM_TRANSPXS_Z            (idx, [1:   4]) = [  3.24870E-01 0.00021 -4.75020E+00 0.00041 ];
CMM_DIFFCOEF              (idx, [1:   4]) = [  1.03713E+00 0.00014 -7.40396E-02 0.00034 ];
CMM_DIFFCOEF_X            (idx, [1:   4]) = [  1.05902E+00 0.00017 -8.17547E-02 0.00033 ];
CMM_DIFFCOEF_Y            (idx, [1:   4]) = [  1.02632E+00 0.00018 -7.01911E-02 0.00037 ];
CMM_DIFFCOEF_Z            (idx, [1:   4]) = [  1.02605E+00 0.00021 -7.01731E-02 0.00041 ];

% Delayed neutron parameters (Meulekamp method):

BETA_EFF                  (idx, [1:  18]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
LAMBDA                    (idx, [1:  18]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
