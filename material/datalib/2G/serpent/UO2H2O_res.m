
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
INPUT_FILE_NAME           (idx, [1:  8])  = 'slab_h2o' ;
WORKING_DIRECTORY         (idx, [1: 33])  = '/home/abrate/SERPENT-2/phytra/h2o' ;
HOSTNAME                  (idx, [1:  7])  = 'vpcen13' ;
CPU_TYPE                  (idx, [1: 41])  = 'Intel(R) Xeon(R) CPU E5-2630 v3 @ 2.40GHz' ;
CPU_MHZ                   (idx, 1)        = 4294967295.0 ;
START_DATE                (idx, [1: 24])  = 'Fri Oct 19 02:33:47 2018' ;
COMPLETE_DATE             (idx, [1: 24])  = 'Fri Oct 19 03:04:58 2018' ;

% Run parameters:

POP                       (idx, 1)        = 100000 ;
CYCLES                    (idx, 1)        = 1000 ;
SKIP                      (idx, 1)        = 500 ;
BATCH_INTERVAL            (idx, 1)        = 1 ;
SRC_NORM_MODE             (idx, 1)        = 2 ;
SEED                      (idx, 1)        = 1539909227 ;
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
OMP_HISTORY_PROFILE       (idx, [1:  30]) = [  1.35442E+00  1.01048E+00  9.97119E-01  9.93778E-01  9.93544E-01  1.02143E+00  9.97955E-01  9.82677E-01  9.98220E-01  9.83732E-01  9.97989E-01  9.87917E-01  9.74039E-01  9.87954E-01  9.81608E-01  9.88555E-01  9.69563E-01  9.88927E-01  9.75858E-01  9.89385E-01  9.70688E-01  9.88898E-01  9.76523E-01  9.90822E-01  9.81058E-01  9.90898E-01  9.77066E-01  9.88386E-01  9.76905E-01  9.83602E-01  ];
SHARE_BUF_ARRAY           (idx, 1)        = 0 ;
SHARE_RES2_ARRAY          (idx, 1)        = 1 ;

% File paths:

XS_DATA_FILE_PATH         (idx, [1: 45])  = '/opt/serpent/xsdata/jeff33/sss2_jeff33.xsdata' ;
DECAY_DATA_FILE_PATH      (idx, [1:  3])  = 'N/A' ;
SFY_DATA_FILE_PATH        (idx, [1:  3])  = 'N/A' ;
NFY_DATA_FILE_PATH        (idx, [1:  3])  = 'N/A' ;
BRA_DATA_FILE_PATH        (idx, [1:  3])  = 'N/A' ;

% Collision and reaction sampling (neutrons/photons):

MIN_MACROXS               (idx, [1:   4]) = [  5.00000E-02 1.7E-09  0.00000E+00 0.0E+00 ];
DT_THRESH                 (idx, [1:  2])  = [  9.00000E-01  9.00000E-01 ];
ST_FRAC                   (idx, [1:   4]) = [  2.99589E-04 0.00061  0.00000E+00 0.0E+00 ];
DT_FRAC                   (idx, [1:   4]) = [  9.99700E-01 1.8E-07  0.00000E+00 0.0E+00 ];
DT_EFF                    (idx, [1:   4]) = [  6.47950E-01 8.1E-05  0.00000E+00 0.0E+00 ];
REA_SAMPLING_EFF          (idx, [1:   4]) = [  1.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
REA_SAMPLING_FAIL         (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
TOT_COL_EFF               (idx, [1:   4]) = [  6.48052E-01 8.1E-05  0.00000E+00 0.0E+00 ];
AVG_TRACKING_LOOPS        (idx, [1:   8]) = [  3.10097E+00 6.9E-05  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
AVG_TRACKS                (idx, [1:   4]) = [  7.50474E+01 0.00014  0.00000E+00 0.0E+00 ];
AVG_REAL_COL              (idx, [1:   4]) = [  7.50474E+01 0.00014  0.00000E+00 0.0E+00 ];
AVG_VIRT_COL              (idx, [1:   4]) = [  4.07569E+01 0.00016  0.00000E+00 0.0E+00 ];
AVG_SURF_CROSS            (idx, [1:   4]) = [  1.02305E-03 0.00312  0.00000E+00 0.0E+00 ];
LOST_PARTICLES            (idx, 1)        = 0 ;

% Run statistics:

CYCLE_IDX                 (idx, 1)        = 1000 ;
SOURCE_POPULATION         (idx, 1)        = 100002620 ;
MEAN_POP_SIZE             (idx, [1:  2])  = [  1.00003E+05 0.00024 ];
MEAN_POP_WGT              (idx, [1:  2])  = [  1.00003E+05 0.00024 ];
SIMULATION_COMPLETED      (idx, 1)        = 1 ;

% Running times:

TOT_CPU_TIME              (idx, 1)        =  8.15880E+02 ;
RUNNING_TIME              (idx, 1)        =  3.11838E+01 ;
INIT_TIME                 (idx, [1:  2])  = [  3.58000E-02  3.58000E-02 ];
PROCESS_TIME              (idx, [1:  2])  = [  3.50000E-04  3.50000E-04 ];
TRANSPORT_CYCLE_TIME      (idx, [1:  3])  = [  3.11477E+01  3.11477E+01  0.00000E+00 ];
MPI_OVERHEAD_TIME         (idx, [1:  2])  = [  0.00000E+00  0.00000E+00 ];
ESTIMATED_RUNNING_TIME    (idx, [1:  2])  = [  3.11822E+01  0.00000E+00 ];
CPU_USAGE                 (idx, 1)        = 26.16354 ;
TRANSPORT_CPU_USAGE       (idx, [1:   2]) = [  2.66144E+01 0.00024 ];
OMP_PARALLEL_FRAC         (idx, 1)        =  8.65600E-01 ;

% Memory usage:

AVAIL_MEM                 (idx, 1)        = 96872.24 ;
ALLOC_MEMSIZE             (idx, 1)        = 1074.40;
MEMSIZE                   (idx, 1)        = 849.31;
XS_MEMSIZE                (idx, 1)        = 174.20;
MAT_MEMSIZE               (idx, 1)        = 18.51;
RES_MEMSIZE               (idx, 1)        = 3.37;
MISC_MEMSIZE              (idx, 1)        = 653.23;
UNKNOWN_MEMSIZE           (idx, 1)        = 0.00;
UNUSED_MEMSIZE            (idx, 1)        = 225.08;

% Geometry parameters:

TOT_CELLS                 (idx, 1)        = 7 ;
UNION_CELLS               (idx, 1)        = 0 ;

% Neutron energy grid:

NEUTRON_ERG_TOL           (idx, 1)        =  0.00000E+00 ;
NEUTRON_ERG_NE            (idx, 1)        = 201574 ;
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
TOT_REA_CHANNELS          (idx, 1)        = 136 ;
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

NORM_COEF                 (idx, [1:   4]) = [  9.97140E-06 9.5E-05  0.00000E+00 0.0E+00 ];

% Analog reaction rate estimators:

CONVERSION_RATIO          (idx, [1:   2]) = [  1.96211E+00 0.00028 ];
U235_FISS                 (idx, [1:   4]) = [  1.87058E-01 0.00023  7.20834E-01 0.00012 ];
U238_FISS                 (idx, [1:   4]) = [  7.24446E-02 0.00039  2.79166E-01 0.00032 ];
U235_CAPT                 (idx, [1:   4]) = [  5.14156E-02 0.00044  6.94389E-02 0.00043 ];
U238_CAPT                 (idx, [1:   4]) = [  4.68062E-01 0.00016  6.32136E-01 0.00012 ];

% Neutron balance (particles/weight):

BALA_SRC_NEUTRON_SRC     (idx, [1:  2])  = [ 0 0.00000E+00 ];
BALA_SRC_NEUTRON_FISS    (idx, [1:  2])  = [ 100002620 1.00000E+08 ];
BALA_SRC_NEUTRON_NXN     (idx, [1:  2])  = [ 0 2.83554E+05 ];
BALA_SRC_NEUTRON_VR      (idx, [1:  2])  = [ 0 0.00000E+00 ];
BALA_SRC_NEUTRON_TOT     (idx, [1:  2])  = [ 100002620 1.00284E+08 ];

BALA_LOSS_NEUTRON_CAPT    (idx, [1:  2])  = [ 74034448 7.42569E+07 ];
BALA_LOSS_NEUTRON_FISS    (idx, [1:  2])  = [ 25966140 2.60246E+07 ];
BALA_LOSS_NEUTRON_LEAK    (idx, [1:  2])  = [ 2032 2.03244E+03 ];
BALA_LOSS_NEUTRON_CUT     (idx, [1:  2])  = [ 0 0.00000E+00 ];
BALA_LOSS_NEUTRON_ERR     (idx, [1:  2])  = [ 0 0.00000E+00 ];
BALA_LOSS_NEUTRON_TOT     (idx, [1:  2])  = [ 100002620 1.00284E+08 ];

BALA_NEUTRON_DIFF         (idx, [1:  2])  = [ 0 -6.55651E-06 ];

% Normalized total reaction rates (neutrons):

TOT_POWER                 (idx, [1:   2]) = [  8.46664E-12 1.0E-04 ];
TOT_POWDENS               (idx, [1:   2]) = [  0.00000E+00 0.0E+00 ];
TOT_GENRATE               (idx, [1:   2]) = [  6.58012E-01 9.9E-05 ];
TOT_FISSRATE              (idx, [1:   2]) = [  2.59525E-01 0.00010 ];
TOT_CAPTRATE              (idx, [1:   2]) = [  7.40455E-01 3.5E-05 ];
TOT_ABSRATE               (idx, [1:   2]) = [  9.99980E-01 4.4E-07 ];
TOT_SRCRATE               (idx, [1:   2]) = [  9.97140E-01 9.5E-05 ];
TOT_FLUX                  (idx, [1:   2]) = [  9.52194E+01 6.8E-05 ];
TOT_PHOTON_PRODRATE       (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
TOT_LEAKRATE              (idx, [1:   2]) = [  2.02638E-05 0.02147 ];
ALBEDO_LEAKRATE           (idx, [1:   2]) = [  0.00000E+00 0.0E+00 ];
TOT_LOSSRATE              (idx, [1:   2]) = [  1.00000E+00 0.0E+00 ];
TOT_CUTRATE               (idx, [1:   2]) = [  0.00000E+00 0.0E+00 ];
TOT_RR                    (idx, [1:   2]) = [  7.50241E+01 0.00012 ];
INI_FMASS                 (idx, 1)        =  0.00000E+00 ;
TOT_FMASS                 (idx, 1)        =  0.00000E+00 ;

% Six-factor formula:

SIX_FF_ETA                (idx, [1:   2]) = [  1.82454E+00 0.00020 ];
SIX_FF_F                  (idx, [1:   2]) = [  2.92036E-01 0.00029 ];
SIX_FF_P                  (idx, [1:   2]) = [  3.04374E-01 0.00023 ];
SIX_FF_EPSILON            (idx, [1:   2]) = [  4.06928E+00 0.00038 ];
SIX_FF_LF                 (idx, [1:   2]) = [  9.99984E-01 3.9E-07 ];
SIX_FF_LT                 (idx, [1:   2]) = [  9.99996E-01 2.1E-07 ];
SIX_FF_KINF               (idx, [1:   2]) = [  6.59860E-01 0.00018 ];
SIX_FF_KEFF               (idx, [1:   2]) = [  6.59847E-01 0.00018 ];

% Fission neutron and energy production:

NUBAR                     (idx, [1:   2]) = [  2.53545E+00 8.9E-06 ];
FISSE                     (idx, [1:   2]) = [  2.03620E+02 1.1E-06 ];

% Criticality eigenvalues:

ANA_KEFF                  (idx, [1:   6]) = [  6.59836E-01 0.00018  6.54904E-01 0.00018  4.94300E-03 0.00222 ];
IMP_KEFF                  (idx, [1:   2]) = [  6.59878E-01 1.0E-04 ];
COL_KEFF                  (idx, [1:   2]) = [  6.59904E-01 0.00012 ];
ABS_KEFF                  (idx, [1:   2]) = [  6.59878E-01 1.0E-04 ];
ABS_KINF                  (idx, [1:   2]) = [  6.59892E-01 1.0E-04 ];
GEOM_ALBEDO               (idx, [1:   6]) = [  1.00000E+00 0.0E+00  1.00000E+00 0.0E+00  1.00000E+00 0.0E+00 ];

% ALF (Average lethargy of neutrons causing fission):
% Based on E0 = 2.000000E+01 MeV

ANA_ALF                   (idx, [1:   2]) = [  9.01442E+00 0.00017 ];
IMP_ALF                   (idx, [1:   2]) = [  9.01496E+00 0.00015 ];

% EALF (Energy corresponding to average lethargy of neutrons causing fission):

ANA_EALF                  (idx, [1:   2]) = [  2.43572E-03 0.00154 ];
IMP_EALF                  (idx, [1:   2]) = [  2.43382E-03 0.00137 ];

% AFGE (Average energy of neutrons causing fission):

ANA_AFGE                  (idx, [1:   2]) = [  9.42531E-01 0.00033 ];
IMP_AFGE                  (idx, [1:   2]) = [  9.42338E-01 0.00019 ];

% Forward-weighted delayed neutron parameters:

FWD_ANA_BETA_ZERO         (idx, [1:  18]) = [  1.47902E-02 0.00111  3.03133E-04 0.00711  1.86524E-03 0.00281  9.75542E-04 0.00386  2.45780E-03 0.00258  4.52288E-03 0.00190  2.18201E-03 0.00265  1.59658E-03 0.00309  8.87013E-04 0.00420 ];
FWD_ANA_LAMBDA            (idx, [1:  18]) = [  6.06132E-01 0.00146  1.24667E-02 7.2E-09  2.82917E-02 0.0E+00  4.25244E-02 0.0E+00  1.33042E-01 0.0E+00  2.92467E-01 0.0E+00  6.66488E-01 2.7E-09  1.63478E+00 4.2E-09  3.55460E+00 2.9E-09 ];

% Beta-eff using Meulekamp's method:

ADJ_MEULEKAMP_BETA_EFF    (idx, [1:  18]) = [  7.77991E-03 0.00179  1.61375E-04 0.01314  9.79460E-04 0.00518  5.11242E-04 0.00722  1.28703E-03 0.00452  2.38318E-03 0.00318  1.14001E-03 0.00468  8.51471E-04 0.00537  4.66140E-04 0.00730 ];
ADJ_MEULEKAMP_LAMBDA      (idx, [1:  18]) = [  6.07980E-01 0.00266  1.24667E-02 7.2E-09  2.82917E-02 0.0E+00  4.25244E-02 0.0E+00  1.33042E-01 0.0E+00  2.92467E-01 0.0E+00  6.66488E-01 2.7E-09  1.63478E+00 4.2E-09  3.55460E+00 2.9E-09 ];

% Adjoint weighted time constants using Nauchi's method:

ADJ_NAUCHI_GEN_TIME       (idx, [1:   6]) = [  3.51055E-05 0.00082  3.50275E-05 0.00082  4.54399E-05 0.00747 ];
ADJ_NAUCHI_LIFETIME       (idx, [1:   6]) = [  2.31632E-05 0.00080  2.31118E-05 0.00080  2.99794E-05 0.00746 ];
ADJ_NAUCHI_BETA_EFF       (idx, [1:  18]) = [  7.49038E-03 0.00224  1.51808E-04 0.01568  9.45641E-04 0.00641  4.92440E-04 0.00878  1.23990E-03 0.00569  2.29799E-03 0.00393  1.09080E-03 0.00564  8.22658E-04 0.00683  4.49132E-04 0.00929 ];
ADJ_NAUCHI_LAMBDA         (idx, [1:  18]) = [  6.08233E-01 0.00338  1.24667E-02 7.1E-09  2.82917E-02 0.0E+00  4.25244E-02 0.0E+00  1.33042E-01 0.0E+00  2.92467E-01 0.0E+00  6.66488E-01 2.7E-09  1.63478E+00 4.2E-09  3.55460E+00 2.9E-09 ];

% Adjoint weighted time constants using IFP:

ADJ_IFP_GEN_TIME          (idx, [1:   6]) = [  3.37506E-05 0.00224  3.36730E-05 0.00227  4.45041E-05 0.02173 ];
ADJ_IFP_LIFETIME          (idx, [1:   6]) = [  2.22693E-05 0.00224  2.22182E-05 0.00227  2.93607E-05 0.02172 ];
ADJ_IFP_IMP_BETA_EFF      (idx, [1:  18]) = [  7.46324E-03 0.00767  1.36733E-04 0.05298  9.32062E-04 0.02188  4.95025E-04 0.02967  1.25424E-03 0.01937  2.33365E-03 0.01428  1.06990E-03 0.02103  7.84238E-04 0.02428  4.57398E-04 0.03181 ];
ADJ_IFP_IMP_LAMBDA        (idx, [1:  18]) = [  6.05248E-01 0.01167  1.24667E-02 0.0E+00  2.82917E-02 0.0E+00  4.25244E-02 0.0E+00  1.33042E-01 0.0E+00  2.92467E-01 0.0E+00  6.66488E-01 2.4E-09  1.63478E+00 4.7E-09  3.55460E+00 2.7E-09 ];
ADJ_IFP_ANA_BETA_EFF      (idx, [1:  18]) = [  7.45573E-03 0.00756  1.36449E-04 0.05246  9.27685E-04 0.02158  4.92577E-04 0.02915  1.25622E-03 0.01919  2.33353E-03 0.01411  1.06648E-03 0.02057  7.90128E-04 0.02349  4.52665E-04 0.03149 ];
ADJ_IFP_ANA_LAMBDA        (idx, [1:  18]) = [  6.04938E-01 0.01149  1.24667E-02 0.0E+00  2.82917E-02 0.0E+00  4.25244E-02 0.0E+00  1.33042E-01 0.0E+00  2.92467E-01 0.0E+00  6.66488E-01 2.5E-09  1.63478E+00 4.7E-09  3.55460E+00 2.7E-09 ];
ADJ_IFP_ROSSI_ALPHA       (idx, [1:   2]) = [ -2.22823E+02 0.00806 ];

% Adjoint weighted time constants using perturbation technique:

ADJ_PERT_GEN_TIME         (idx, [1:   2]) = [  3.45296E-05 0.00053 ];
ADJ_PERT_LIFETIME         (idx, [1:   2]) = [  2.27832E-05 0.00051 ];
ADJ_PERT_BETA_EFF         (idx, [1:   2]) = [  7.47329E-03 0.00149 ];
ADJ_PERT_ROSSI_ALPHA      (idx, [1:   2]) = [ -2.16498E+02 0.00159 ];

% Inverse neutron speed :

ANA_INV_SPD               (idx, [1:   2]) = [  4.79418E-07 0.00033 ];

% Analog slowing-down and thermal neutron lifetime (total/prompt/delayed):

ANA_SLOW_TIME             (idx, [1:   6]) = [  1.65501E-06 0.00022  1.65442E-06 0.00022  1.72236E-06 0.00231 ];
ANA_THERM_TIME            (idx, [1:   6]) = [  1.47455E-04 0.00021  1.47546E-04 0.00021  1.37088E-04 0.00225 ];
ANA_THERM_FRAC            (idx, [1:   6]) = [  3.04429E-01 0.00023  3.06326E-01 0.00023  1.78275E-01 0.00219 ];
ANA_DELAYED_EMTIME        (idx, [1:   2]) = [  1.02440E+01 0.00239 ];
ANA_MEAN_NCOL             (idx, [1:   4]) = [  7.50474E+01 0.00014  4.83032E+01 0.00028 ];

% Group constant generation:

GC_UNIVERSE_NAME          (idx, [1:  3])  = 'UO2' ;

% Micro- and macro-group structures:

MICRO_NG                  (idx, 1)        = 70 ;
MICRO_E                   (idx, [1:  71]) = [  1.00000E-11  5.00000E-09  1.00000E-08  1.50000E-08  2.00000E-08  2.50000E-08  3.00000E-08  3.50000E-08  4.20000E-08  5.00000E-08  5.80000E-08  6.70000E-08  8.00000E-08  1.00000E-07  1.40000E-07  1.80000E-07  2.20000E-07  2.50000E-07  2.80000E-07  3.00000E-07  3.20000E-07  3.50000E-07  4.00000E-07  5.00000E-07  6.25000E-07  7.80000E-07  8.50000E-07  9.10000E-07  9.50000E-07  9.72000E-07  9.96000E-07  1.02000E-06  1.04500E-06  1.07100E-06  1.09700E-06  1.12300E-06  1.15000E-06  1.30000E-06  1.50000E-06  1.85500E-06  2.10000E-06  2.60000E-06  3.30000E-06  4.00000E-06  9.87700E-06  1.59680E-05  2.77000E-05  4.80520E-05  7.55014E-05  1.48728E-04  3.67262E-04  9.06898E-04  1.42510E-03  2.23945E-03  3.51910E-03  5.50000E-03  9.11800E-03  1.50300E-02  2.47800E-02  4.08500E-02  6.74300E-02  1.11000E-01  1.83000E-01  3.02500E-01  5.00000E-01  8.21000E-01  1.35300E+00  2.23100E+00  3.67900E+00  6.06550E+00  2.00000E+01 ];

MACRO_NG                  (idx, 1)        = 2 ;
MACRO_E                   (idx, [1:   3]) = [  1.00000E+37  6.25000E-07  0.00000E+00 ];

% Micro-group spectrum:

INF_MICRO_FLX             (idx, [1: 140]) = [  4.52801E+05 0.00070  2.06997E+06 0.00041  4.91248E+06 0.00028  5.36467E+06 0.00024  6.34426E+06 0.00022  1.49506E+07 0.00019  1.17486E+07 0.00019  1.66659E+07 0.00020  1.66802E+07 0.00024  1.54899E+07 0.00027  1.42271E+07 0.00029  1.23637E+07 0.00030  9.95077E+06 0.00032  7.59558E+06 0.00032  5.58995E+06 0.00039  3.54619E+06 0.00046  2.55902E+06 0.00047  1.85137E+06 0.00051  1.27646E+06 0.00060  1.39208E+06 0.00067  6.26157E+05 0.00092  2.36163E+05 0.00129  1.06692E+05 0.00137  9.24659E+04 0.00169  7.34245E+04 0.00152  8.46435E+04 0.00206  1.04196E+05 0.00154  3.45805E+04 0.00225  4.91343E+04 0.00190  5.25641E+04 0.00166  3.00380E+04 0.00220  5.69646E+04 0.00197  3.92125E+04 0.00178  2.91484E+04 0.00273  4.52212E+03 0.00413  4.34829E+03 0.00405  4.52141E+03 0.00427  4.69416E+03 0.00434  4.67286E+03 0.00409  4.59588E+03 0.00395  4.79689E+03 0.00404  4.48345E+03 0.00299  8.41668E+03 0.00268  1.32293E+04 0.00317  1.66198E+04 0.00246  4.17536E+04 0.00263  3.90021E+04 0.00214  3.36082E+04 0.00224  1.65565E+04 0.00230  9.58592E+03 0.00278  6.28192E+03 0.00335  6.24927E+03 0.00311  1.00119E+04 0.00284  1.18614E+04 0.00251  2.20198E+04 0.00174  3.75933E+04 0.00166  7.46397E+04 0.00112  6.24237E+04 0.00114  5.05628E+04 0.00109  4.02024E+04 0.00118  3.78910E+04 0.00124  3.80956E+04 0.00119  3.26204E+04 0.00127  2.16100E+04 0.00121  1.95427E+04 0.00133  1.65737E+04 0.00130  1.31173E+04 0.00150  9.15641E+03 0.00176  4.86861E+03 0.00188  1.18288E+03 0.00298 ];

% Integral parameters:

INF_KINF                  (idx, [1:   2]) = [  8.43912E-01 0.00019 ];

% Flux spectra in infinite geometry:

INF_FLX                   (idx, [1:   4]) = [  7.81576E+01 0.00022  3.06749E-01 0.00068 ];
INF_FISS_FLX              (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Reaction cross sections:

INF_TOT                   (idx, [1:   4]) = [  4.15999E-01 1.3E-05  6.75547E-01 8.9E-05 ];
INF_CAPT                  (idx, [1:   4]) = [  6.41099E-03 7.3E-05  7.15833E-02 0.00019 ];
INF_ABS                   (idx, [1:   4]) = [  8.87838E-03 6.5E-05  2.88966E-01 0.00020 ];
INF_FISS                  (idx, [1:   4]) = [  2.46740E-03 0.00011  2.17382E-01 0.00021 ];
INF_NSF                   (idx, [1:   4]) = [  6.34984E-03 0.00011  5.27237E-01 0.00021 ];
INF_NUBAR                 (idx, [1:   4]) = [  2.57350E+00 1.0E-05  2.42539E+00 4.2E-09 ];
INF_KAPPA                 (idx, [1:   4]) = [  2.04087E+02 1.2E-06  2.02270E+02 0.0E+00 ];
INF_INVV                  (idx, [1:   4]) = [  5.85231E-09 0.00035  2.63838E-06 0.00018 ];

% Total scattering cross sections:

INF_SCATT0                (idx, [1:   4]) = [  4.07121E-01 1.4E-05  3.86606E-01 0.00015 ];
INF_SCATT1                (idx, [1:   4]) = [  3.70921E-02 0.00015  7.86255E-03 0.00865 ];
INF_SCATT2                (idx, [1:   4]) = [  1.63277E-02 0.00020  1.99532E-04 0.23653 ];
INF_SCATT3                (idx, [1:   4]) = [  5.52822E-03 0.00053  6.50407E-05 0.70588 ];
INF_SCATT4                (idx, [1:   4]) = [  3.20279E-03 0.00082  3.92045E-05 1.00000 ];
INF_SCATT5                (idx, [1:   4]) = [  1.50872E-03 0.00152  6.50427E-06 1.00000 ];
INF_SCATT6                (idx, [1:   4]) = [  7.57631E-04 0.00283  1.65527E-05 1.00000 ];
INF_SCATT7                (idx, [1:   4]) = [  3.11643E-04 0.00685  7.40486E-06 1.00000 ];

% Total scattering production cross sections:

INF_SCATTP0               (idx, [1:   4]) = [  4.07157E-01 1.4E-05  3.86606E-01 0.00015 ];
INF_SCATTP1               (idx, [1:   4]) = [  3.70927E-02 0.00015  7.86255E-03 0.00865 ];
INF_SCATTP2               (idx, [1:   4]) = [  1.63278E-02 0.00020  1.99532E-04 0.23653 ];
INF_SCATTP3               (idx, [1:   4]) = [  5.52822E-03 0.00053  6.50407E-05 0.70588 ];
INF_SCATTP4               (idx, [1:   4]) = [  3.20277E-03 0.00082  3.92045E-05 1.00000 ];
INF_SCATTP5               (idx, [1:   4]) = [  1.50872E-03 0.00152  6.50427E-06 1.00000 ];
INF_SCATTP6               (idx, [1:   4]) = [  7.57608E-04 0.00283  1.65527E-05 1.00000 ];
INF_SCATTP7               (idx, [1:   4]) = [  3.11619E-04 0.00684  7.40486E-06 1.00000 ];

% Diffusion parameters:

INF_TRANSPXS              (idx, [1:   4]) = [  3.42056E-01 3.6E-05  6.33528E-01 0.00017 ];
INF_DIFFCOEF              (idx, [1:   4]) = [  9.74501E-01 3.6E-05  5.26155E-01 0.00017 ];

% Reduced absoption and removal:

INF_RABSXS                (idx, [1:   4]) = [  8.84221E-03 6.6E-05  2.88966E-01 0.00020 ];
INF_REMXS                 (idx, [1:   4]) = [  8.90666E-03 9.3E-05  2.89661E-01 0.00030 ];

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

INF_S0                    (idx, [1:   8]) = [  4.07092E-01 1.5E-05  2.82110E-05 0.00250  7.19945E-04 0.00678  3.85886E-01 0.00015 ];
INF_S1                    (idx, [1:   8]) = [  3.70999E-02 0.00015 -7.81982E-06 0.00412  1.05911E-05 0.25072  7.85196E-03 0.00870 ];
INF_S2                    (idx, [1:   8]) = [  1.63282E-02 0.00020 -4.60016E-07 0.06515 -2.40302E-05 0.07488  2.23563E-04 0.21109 ];
INF_S3                    (idx, [1:   8]) = [  5.52827E-03 0.00053 -4.77995E-08 0.43957 -1.82321E-05 0.10549  8.32728E-05 0.54879 ];
INF_S4                    (idx, [1:   8]) = [  3.20283E-03 0.00082 -4.29127E-08 0.46795 -1.07526E-05 0.13333  4.99571E-05 0.93626 ];
INF_S5                    (idx, [1:   8]) = [  1.50877E-03 0.00152 -4.66612E-08 0.46723 -8.79102E-06 0.16329  1.52953E-05 1.00000 ];
INF_S6                    (idx, [1:   8]) = [  7.57638E-04 0.00283 -7.30702E-09 1.00000 -4.74747E-06 0.27115  2.13002E-05 1.00000 ];
INF_S7                    (idx, [1:   8]) = [  3.11675E-04 0.00685 -3.18401E-08 0.55022 -3.12466E-06 0.33851  1.05295E-05 1.00000 ];

% Scattering production matrixes:

INF_SP0                   (idx, [1:   8]) = [  4.07129E-01 1.5E-05  2.82110E-05 0.00250  7.19945E-04 0.00678  3.85886E-01 0.00015 ];
INF_SP1                   (idx, [1:   8]) = [  3.71006E-02 0.00015 -7.81982E-06 0.00412  1.05911E-05 0.25072  7.85196E-03 0.00870 ];
INF_SP2                   (idx, [1:   8]) = [  1.63283E-02 0.00020 -4.60016E-07 0.06515 -2.40302E-05 0.07488  2.23563E-04 0.21109 ];
INF_SP3                   (idx, [1:   8]) = [  5.52827E-03 0.00053 -4.77995E-08 0.43957 -1.82321E-05 0.10549  8.32728E-05 0.54879 ];
INF_SP4                   (idx, [1:   8]) = [  3.20281E-03 0.00082 -4.29127E-08 0.46795 -1.07526E-05 0.13333  4.99571E-05 0.93626 ];
INF_SP5                   (idx, [1:   8]) = [  1.50876E-03 0.00151 -4.66612E-08 0.46723 -8.79102E-06 0.16329  1.52953E-05 1.00000 ];
INF_SP6                   (idx, [1:   8]) = [  7.57616E-04 0.00283 -7.30702E-09 1.00000 -4.74747E-06 0.27115  2.13002E-05 1.00000 ];
INF_SP7                   (idx, [1:   8]) = [  3.11651E-04 0.00684 -3.18401E-08 0.55022 -3.12466E-06 0.33851  1.05295E-05 1.00000 ];

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

CMM_TRANSPXS              (idx, [1:   4]) = [  4.08102E-01 0.00019  2.12974E-02 0.00072 ];
CMM_TRANSPXS_X            (idx, [1:   4]) = [  4.44436E-01 0.00024  2.88684E-02 0.00142 ];
CMM_TRANSPXS_Y            (idx, [1:   4]) = [  3.92054E-01 0.00031  1.88300E-02 0.00066 ];
CMM_TRANSPXS_Z            (idx, [1:   4]) = [  3.92099E-01 0.00025  1.88280E-02 0.00092 ];
CMM_DIFFCOEF              (idx, [1:   4]) = [  8.16790E-01 0.00019  1.56518E+01 0.00072 ];
CMM_DIFFCOEF_X            (idx, [1:   4]) = [  7.50016E-01 0.00024  1.15478E+01 0.00141 ];
CMM_DIFFCOEF_Y            (idx, [1:   4]) = [  8.50226E-01 0.00031  1.77026E+01 0.00066 ];
CMM_DIFFCOEF_Z            (idx, [1:   4]) = [  8.50128E-01 0.00025  1.77049E+01 0.00092 ];

% Delayed neutron parameters (Meulekamp method):

BETA_EFF                  (idx, [1:  18]) = [  7.77991E-03 0.00179  1.61375E-04 0.01314  9.79460E-04 0.00518  5.11242E-04 0.00722  1.28703E-03 0.00452  2.38318E-03 0.00318  1.14001E-03 0.00468  8.51471E-04 0.00537  4.66140E-04 0.00730 ];
LAMBDA                    (idx, [1:  18]) = [  6.07980E-01 0.00266  1.24667E-02 7.2E-09  2.82917E-02 0.0E+00  4.25244E-02 0.0E+00  1.33042E-01 0.0E+00  2.92467E-01 0.0E+00  6.66488E-01 2.7E-09  1.63478E+00 4.2E-09  3.55460E+00 2.9E-09 ];


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
INPUT_FILE_NAME           (idx, [1:  8])  = 'slab_h2o' ;
WORKING_DIRECTORY         (idx, [1: 33])  = '/home/abrate/SERPENT-2/phytra/h2o' ;
HOSTNAME                  (idx, [1:  7])  = 'vpcen13' ;
CPU_TYPE                  (idx, [1: 41])  = 'Intel(R) Xeon(R) CPU E5-2630 v3 @ 2.40GHz' ;
CPU_MHZ                   (idx, 1)        = 4294967295.0 ;
START_DATE                (idx, [1: 24])  = 'Fri Oct 19 02:33:47 2018' ;
COMPLETE_DATE             (idx, [1: 24])  = 'Fri Oct 19 03:04:58 2018' ;

% Run parameters:

POP                       (idx, 1)        = 100000 ;
CYCLES                    (idx, 1)        = 1000 ;
SKIP                      (idx, 1)        = 500 ;
BATCH_INTERVAL            (idx, 1)        = 1 ;
SRC_NORM_MODE             (idx, 1)        = 2 ;
SEED                      (idx, 1)        = 1539909227 ;
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
OMP_HISTORY_PROFILE       (idx, [1:  30]) = [  1.35442E+00  1.01048E+00  9.97119E-01  9.93778E-01  9.93544E-01  1.02143E+00  9.97955E-01  9.82677E-01  9.98220E-01  9.83732E-01  9.97989E-01  9.87917E-01  9.74039E-01  9.87954E-01  9.81608E-01  9.88555E-01  9.69563E-01  9.88927E-01  9.75858E-01  9.89385E-01  9.70688E-01  9.88898E-01  9.76523E-01  9.90822E-01  9.81058E-01  9.90898E-01  9.77066E-01  9.88386E-01  9.76905E-01  9.83602E-01  ];
SHARE_BUF_ARRAY           (idx, 1)        = 0 ;
SHARE_RES2_ARRAY          (idx, 1)        = 1 ;

% File paths:

XS_DATA_FILE_PATH         (idx, [1: 45])  = '/opt/serpent/xsdata/jeff33/sss2_jeff33.xsdata' ;
DECAY_DATA_FILE_PATH      (idx, [1:  3])  = 'N/A' ;
SFY_DATA_FILE_PATH        (idx, [1:  3])  = 'N/A' ;
NFY_DATA_FILE_PATH        (idx, [1:  3])  = 'N/A' ;
BRA_DATA_FILE_PATH        (idx, [1:  3])  = 'N/A' ;

% Collision and reaction sampling (neutrons/photons):

MIN_MACROXS               (idx, [1:   4]) = [  5.00000E-02 1.7E-09  0.00000E+00 0.0E+00 ];
DT_THRESH                 (idx, [1:  2])  = [  9.00000E-01  9.00000E-01 ];
ST_FRAC                   (idx, [1:   4]) = [  2.99589E-04 0.00061  0.00000E+00 0.0E+00 ];
DT_FRAC                   (idx, [1:   4]) = [  9.99700E-01 1.8E-07  0.00000E+00 0.0E+00 ];
DT_EFF                    (idx, [1:   4]) = [  6.47950E-01 8.1E-05  0.00000E+00 0.0E+00 ];
REA_SAMPLING_EFF          (idx, [1:   4]) = [  1.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
REA_SAMPLING_FAIL         (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
TOT_COL_EFF               (idx, [1:   4]) = [  6.48052E-01 8.1E-05  0.00000E+00 0.0E+00 ];
AVG_TRACKING_LOOPS        (idx, [1:   8]) = [  3.10097E+00 6.9E-05  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
AVG_TRACKS                (idx, [1:   4]) = [  7.50474E+01 0.00014  0.00000E+00 0.0E+00 ];
AVG_REAL_COL              (idx, [1:   4]) = [  7.50474E+01 0.00014  0.00000E+00 0.0E+00 ];
AVG_VIRT_COL              (idx, [1:   4]) = [  4.07569E+01 0.00016  0.00000E+00 0.0E+00 ];
AVG_SURF_CROSS            (idx, [1:   4]) = [  1.02305E-03 0.00312  0.00000E+00 0.0E+00 ];
LOST_PARTICLES            (idx, 1)        = 0 ;

% Run statistics:

CYCLE_IDX                 (idx, 1)        = 1000 ;
SOURCE_POPULATION         (idx, 1)        = 100002620 ;
MEAN_POP_SIZE             (idx, [1:  2])  = [  1.00003E+05 0.00024 ];
MEAN_POP_WGT              (idx, [1:  2])  = [  1.00003E+05 0.00024 ];
SIMULATION_COMPLETED      (idx, 1)        = 1 ;

% Running times:

TOT_CPU_TIME              (idx, 1)        =  8.15880E+02 ;
RUNNING_TIME              (idx, 1)        =  3.11838E+01 ;
INIT_TIME                 (idx, [1:  2])  = [  3.58000E-02  3.58000E-02 ];
PROCESS_TIME              (idx, [1:  2])  = [  3.50000E-04  3.50000E-04 ];
TRANSPORT_CYCLE_TIME      (idx, [1:  3])  = [  3.11477E+01  3.11477E+01  0.00000E+00 ];
MPI_OVERHEAD_TIME         (idx, [1:  2])  = [  0.00000E+00  0.00000E+00 ];
ESTIMATED_RUNNING_TIME    (idx, [1:  2])  = [  3.11822E+01  0.00000E+00 ];
CPU_USAGE                 (idx, 1)        = 26.16354 ;
TRANSPORT_CPU_USAGE       (idx, [1:   2]) = [  2.66144E+01 0.00024 ];
OMP_PARALLEL_FRAC         (idx, 1)        =  8.65600E-01 ;

% Memory usage:

AVAIL_MEM                 (idx, 1)        = 96872.24 ;
ALLOC_MEMSIZE             (idx, 1)        = 1074.40;
MEMSIZE                   (idx, 1)        = 849.31;
XS_MEMSIZE                (idx, 1)        = 174.20;
MAT_MEMSIZE               (idx, 1)        = 18.51;
RES_MEMSIZE               (idx, 1)        = 3.37;
MISC_MEMSIZE              (idx, 1)        = 653.23;
UNKNOWN_MEMSIZE           (idx, 1)        = 0.00;
UNUSED_MEMSIZE            (idx, 1)        = 225.08;

% Geometry parameters:

TOT_CELLS                 (idx, 1)        = 7 ;
UNION_CELLS               (idx, 1)        = 0 ;

% Neutron energy grid:

NEUTRON_ERG_TOL           (idx, 1)        =  0.00000E+00 ;
NEUTRON_ERG_NE            (idx, 1)        = 201574 ;
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
TOT_REA_CHANNELS          (idx, 1)        = 136 ;
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

NORM_COEF                 (idx, [1:   4]) = [  9.97140E-06 9.5E-05  0.00000E+00 0.0E+00 ];

% Analog reaction rate estimators:

CONVERSION_RATIO          (idx, [1:   2]) = [  1.96211E+00 0.00028 ];
U235_FISS                 (idx, [1:   4]) = [  1.87058E-01 0.00023  7.20834E-01 0.00012 ];
U238_FISS                 (idx, [1:   4]) = [  7.24446E-02 0.00039  2.79166E-01 0.00032 ];
U235_CAPT                 (idx, [1:   4]) = [  5.14156E-02 0.00044  6.94389E-02 0.00043 ];
U238_CAPT                 (idx, [1:   4]) = [  4.68062E-01 0.00016  6.32136E-01 0.00012 ];

% Neutron balance (particles/weight):

BALA_SRC_NEUTRON_SRC     (idx, [1:  2])  = [ 0 0.00000E+00 ];
BALA_SRC_NEUTRON_FISS    (idx, [1:  2])  = [ 100002620 1.00000E+08 ];
BALA_SRC_NEUTRON_NXN     (idx, [1:  2])  = [ 0 2.83554E+05 ];
BALA_SRC_NEUTRON_VR      (idx, [1:  2])  = [ 0 0.00000E+00 ];
BALA_SRC_NEUTRON_TOT     (idx, [1:  2])  = [ 100002620 1.00284E+08 ];

BALA_LOSS_NEUTRON_CAPT    (idx, [1:  2])  = [ 74034448 7.42569E+07 ];
BALA_LOSS_NEUTRON_FISS    (idx, [1:  2])  = [ 25966140 2.60246E+07 ];
BALA_LOSS_NEUTRON_LEAK    (idx, [1:  2])  = [ 2032 2.03244E+03 ];
BALA_LOSS_NEUTRON_CUT     (idx, [1:  2])  = [ 0 0.00000E+00 ];
BALA_LOSS_NEUTRON_ERR     (idx, [1:  2])  = [ 0 0.00000E+00 ];
BALA_LOSS_NEUTRON_TOT     (idx, [1:  2])  = [ 100002620 1.00284E+08 ];

BALA_NEUTRON_DIFF         (idx, [1:  2])  = [ 0 -6.55651E-06 ];

% Normalized total reaction rates (neutrons):

TOT_POWER                 (idx, [1:   2]) = [  8.46664E-12 1.0E-04 ];
TOT_POWDENS               (idx, [1:   2]) = [  0.00000E+00 0.0E+00 ];
TOT_GENRATE               (idx, [1:   2]) = [  6.58012E-01 9.9E-05 ];
TOT_FISSRATE              (idx, [1:   2]) = [  2.59525E-01 0.00010 ];
TOT_CAPTRATE              (idx, [1:   2]) = [  7.40455E-01 3.5E-05 ];
TOT_ABSRATE               (idx, [1:   2]) = [  9.99980E-01 4.4E-07 ];
TOT_SRCRATE               (idx, [1:   2]) = [  9.97140E-01 9.5E-05 ];
TOT_FLUX                  (idx, [1:   2]) = [  9.52194E+01 6.8E-05 ];
TOT_PHOTON_PRODRATE       (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
TOT_LEAKRATE              (idx, [1:   2]) = [  2.02638E-05 0.02147 ];
ALBEDO_LEAKRATE           (idx, [1:   2]) = [  0.00000E+00 0.0E+00 ];
TOT_LOSSRATE              (idx, [1:   2]) = [  1.00000E+00 0.0E+00 ];
TOT_CUTRATE               (idx, [1:   2]) = [  0.00000E+00 0.0E+00 ];
TOT_RR                    (idx, [1:   2]) = [  7.50241E+01 0.00012 ];
INI_FMASS                 (idx, 1)        =  0.00000E+00 ;
TOT_FMASS                 (idx, 1)        =  0.00000E+00 ;

% Six-factor formula:

SIX_FF_ETA                (idx, [1:   2]) = [  1.82454E+00 0.00020 ];
SIX_FF_F                  (idx, [1:   2]) = [  2.92036E-01 0.00029 ];
SIX_FF_P                  (idx, [1:   2]) = [  3.04374E-01 0.00023 ];
SIX_FF_EPSILON            (idx, [1:   2]) = [  4.06928E+00 0.00038 ];
SIX_FF_LF                 (idx, [1:   2]) = [  9.99984E-01 3.9E-07 ];
SIX_FF_LT                 (idx, [1:   2]) = [  9.99996E-01 2.1E-07 ];
SIX_FF_KINF               (idx, [1:   2]) = [  6.59860E-01 0.00018 ];
SIX_FF_KEFF               (idx, [1:   2]) = [  6.59847E-01 0.00018 ];

% Fission neutron and energy production:

NUBAR                     (idx, [1:   2]) = [  2.53545E+00 8.9E-06 ];
FISSE                     (idx, [1:   2]) = [  2.03620E+02 1.1E-06 ];

% Criticality eigenvalues:

ANA_KEFF                  (idx, [1:   6]) = [  6.59836E-01 0.00018  6.54904E-01 0.00018  4.94300E-03 0.00222 ];
IMP_KEFF                  (idx, [1:   2]) = [  6.59878E-01 1.0E-04 ];
COL_KEFF                  (idx, [1:   2]) = [  6.59904E-01 0.00012 ];
ABS_KEFF                  (idx, [1:   2]) = [  6.59878E-01 1.0E-04 ];
ABS_KINF                  (idx, [1:   2]) = [  6.59892E-01 1.0E-04 ];
GEOM_ALBEDO               (idx, [1:   6]) = [  1.00000E+00 0.0E+00  1.00000E+00 0.0E+00  1.00000E+00 0.0E+00 ];

% ALF (Average lethargy of neutrons causing fission):
% Based on E0 = 2.000000E+01 MeV

ANA_ALF                   (idx, [1:   2]) = [  9.01442E+00 0.00017 ];
IMP_ALF                   (idx, [1:   2]) = [  9.01496E+00 0.00015 ];

% EALF (Energy corresponding to average lethargy of neutrons causing fission):

ANA_EALF                  (idx, [1:   2]) = [  2.43572E-03 0.00154 ];
IMP_EALF                  (idx, [1:   2]) = [  2.43382E-03 0.00137 ];

% AFGE (Average energy of neutrons causing fission):

ANA_AFGE                  (idx, [1:   2]) = [  9.42531E-01 0.00033 ];
IMP_AFGE                  (idx, [1:   2]) = [  9.42338E-01 0.00019 ];

% Forward-weighted delayed neutron parameters:

FWD_ANA_BETA_ZERO         (idx, [1:  18]) = [  1.47902E-02 0.00111  3.03133E-04 0.00711  1.86524E-03 0.00281  9.75542E-04 0.00386  2.45780E-03 0.00258  4.52288E-03 0.00190  2.18201E-03 0.00265  1.59658E-03 0.00309  8.87013E-04 0.00420 ];
FWD_ANA_LAMBDA            (idx, [1:  18]) = [  6.06132E-01 0.00146  1.24667E-02 7.2E-09  2.82917E-02 0.0E+00  4.25244E-02 0.0E+00  1.33042E-01 0.0E+00  2.92467E-01 0.0E+00  6.66488E-01 2.7E-09  1.63478E+00 4.2E-09  3.55460E+00 2.9E-09 ];

% Beta-eff using Meulekamp's method:

ADJ_MEULEKAMP_BETA_EFF    (idx, [1:  18]) = [  7.77991E-03 0.00179  1.61375E-04 0.01314  9.79460E-04 0.00518  5.11242E-04 0.00722  1.28703E-03 0.00452  2.38318E-03 0.00318  1.14001E-03 0.00468  8.51471E-04 0.00537  4.66140E-04 0.00730 ];
ADJ_MEULEKAMP_LAMBDA      (idx, [1:  18]) = [  6.07980E-01 0.00266  1.24667E-02 7.2E-09  2.82917E-02 0.0E+00  4.25244E-02 0.0E+00  1.33042E-01 0.0E+00  2.92467E-01 0.0E+00  6.66488E-01 2.7E-09  1.63478E+00 4.2E-09  3.55460E+00 2.9E-09 ];

% Adjoint weighted time constants using Nauchi's method:

ADJ_NAUCHI_GEN_TIME       (idx, [1:   6]) = [  3.51055E-05 0.00082  3.50275E-05 0.00082  4.54399E-05 0.00747 ];
ADJ_NAUCHI_LIFETIME       (idx, [1:   6]) = [  2.31632E-05 0.00080  2.31118E-05 0.00080  2.99794E-05 0.00746 ];
ADJ_NAUCHI_BETA_EFF       (idx, [1:  18]) = [  7.49038E-03 0.00224  1.51808E-04 0.01568  9.45641E-04 0.00641  4.92440E-04 0.00878  1.23990E-03 0.00569  2.29799E-03 0.00393  1.09080E-03 0.00564  8.22658E-04 0.00683  4.49132E-04 0.00929 ];
ADJ_NAUCHI_LAMBDA         (idx, [1:  18]) = [  6.08233E-01 0.00338  1.24667E-02 7.1E-09  2.82917E-02 0.0E+00  4.25244E-02 0.0E+00  1.33042E-01 0.0E+00  2.92467E-01 0.0E+00  6.66488E-01 2.7E-09  1.63478E+00 4.2E-09  3.55460E+00 2.9E-09 ];

% Adjoint weighted time constants using IFP:

ADJ_IFP_GEN_TIME          (idx, [1:   6]) = [  3.37506E-05 0.00224  3.36730E-05 0.00227  4.45041E-05 0.02173 ];
ADJ_IFP_LIFETIME          (idx, [1:   6]) = [  2.22693E-05 0.00224  2.22182E-05 0.00227  2.93607E-05 0.02172 ];
ADJ_IFP_IMP_BETA_EFF      (idx, [1:  18]) = [  7.46324E-03 0.00767  1.36733E-04 0.05298  9.32062E-04 0.02188  4.95025E-04 0.02967  1.25424E-03 0.01937  2.33365E-03 0.01428  1.06990E-03 0.02103  7.84238E-04 0.02428  4.57398E-04 0.03181 ];
ADJ_IFP_IMP_LAMBDA        (idx, [1:  18]) = [  6.05248E-01 0.01167  1.24667E-02 0.0E+00  2.82917E-02 0.0E+00  4.25244E-02 0.0E+00  1.33042E-01 0.0E+00  2.92467E-01 0.0E+00  6.66488E-01 2.4E-09  1.63478E+00 4.7E-09  3.55460E+00 2.7E-09 ];
ADJ_IFP_ANA_BETA_EFF      (idx, [1:  18]) = [  7.45573E-03 0.00756  1.36449E-04 0.05246  9.27685E-04 0.02158  4.92577E-04 0.02915  1.25622E-03 0.01919  2.33353E-03 0.01411  1.06648E-03 0.02057  7.90128E-04 0.02349  4.52665E-04 0.03149 ];
ADJ_IFP_ANA_LAMBDA        (idx, [1:  18]) = [  6.04938E-01 0.01149  1.24667E-02 0.0E+00  2.82917E-02 0.0E+00  4.25244E-02 0.0E+00  1.33042E-01 0.0E+00  2.92467E-01 0.0E+00  6.66488E-01 2.5E-09  1.63478E+00 4.7E-09  3.55460E+00 2.7E-09 ];
ADJ_IFP_ROSSI_ALPHA       (idx, [1:   2]) = [ -2.22823E+02 0.00806 ];

% Adjoint weighted time constants using perturbation technique:

ADJ_PERT_GEN_TIME         (idx, [1:   2]) = [  3.45296E-05 0.00053 ];
ADJ_PERT_LIFETIME         (idx, [1:   2]) = [  2.27832E-05 0.00051 ];
ADJ_PERT_BETA_EFF         (idx, [1:   2]) = [  7.47329E-03 0.00149 ];
ADJ_PERT_ROSSI_ALPHA      (idx, [1:   2]) = [ -2.16498E+02 0.00159 ];

% Inverse neutron speed :

ANA_INV_SPD               (idx, [1:   2]) = [  4.79418E-07 0.00033 ];

% Analog slowing-down and thermal neutron lifetime (total/prompt/delayed):

ANA_SLOW_TIME             (idx, [1:   6]) = [  1.65501E-06 0.00022  1.65442E-06 0.00022  1.72236E-06 0.00231 ];
ANA_THERM_TIME            (idx, [1:   6]) = [  1.47455E-04 0.00021  1.47546E-04 0.00021  1.37088E-04 0.00225 ];
ANA_THERM_FRAC            (idx, [1:   6]) = [  3.04429E-01 0.00023  3.06326E-01 0.00023  1.78275E-01 0.00219 ];
ANA_DELAYED_EMTIME        (idx, [1:   2]) = [  1.02440E+01 0.00239 ];
ANA_MEAN_NCOL             (idx, [1:   4]) = [  7.50474E+01 0.00014  4.83032E+01 0.00028 ];

% Group constant generation:

GC_UNIVERSE_NAME          (idx, [1:  3])  = 'H2O' ;

% Micro- and macro-group structures:

MICRO_NG                  (idx, 1)        = 70 ;
MICRO_E                   (idx, [1:  71]) = [  1.00000E-11  5.00000E-09  1.00000E-08  1.50000E-08  2.00000E-08  2.50000E-08  3.00000E-08  3.50000E-08  4.20000E-08  5.00000E-08  5.80000E-08  6.70000E-08  8.00000E-08  1.00000E-07  1.40000E-07  1.80000E-07  2.20000E-07  2.50000E-07  2.80000E-07  3.00000E-07  3.20000E-07  3.50000E-07  4.00000E-07  5.00000E-07  6.25000E-07  7.80000E-07  8.50000E-07  9.10000E-07  9.50000E-07  9.72000E-07  9.96000E-07  1.02000E-06  1.04500E-06  1.07100E-06  1.09700E-06  1.12300E-06  1.15000E-06  1.30000E-06  1.50000E-06  1.85500E-06  2.10000E-06  2.60000E-06  3.30000E-06  4.00000E-06  9.87700E-06  1.59680E-05  2.77000E-05  4.80520E-05  7.55014E-05  1.48728E-04  3.67262E-04  9.06898E-04  1.42510E-03  2.23945E-03  3.51910E-03  5.50000E-03  9.11800E-03  1.50300E-02  2.47800E-02  4.08500E-02  6.74300E-02  1.11000E-01  1.83000E-01  3.02500E-01  5.00000E-01  8.21000E-01  1.35300E+00  2.23100E+00  3.67900E+00  6.06550E+00  2.00000E+01 ];

MACRO_NG                  (idx, 1)        = 2 ;
MACRO_E                   (idx, [1:   3]) = [  1.00000E+37  6.25000E-07  0.00000E+00 ];

% Micro-group spectrum:

INF_MICRO_FLX             (idx, [1: 140]) = [  8.59043E+04 0.00174  3.28275E+05 0.00093  6.57008E+05 0.00067  6.72379E+05 0.00066  6.07455E+05 0.00064  7.60239E+05 0.00060  5.42427E+05 0.00059  5.25151E+05 0.00059  4.38116E+05 0.00056  3.75124E+05 0.00058  3.34594E+05 0.00048  3.04993E+05 0.00053  2.86653E+05 0.00052  2.74452E+05 0.00046  2.68650E+05 0.00049  2.32022E+05 0.00048  2.31234E+05 0.00044  2.28363E+05 0.00054  2.25835E+05 0.00043  4.44493E+05 0.00049  4.36868E+05 0.00047  3.22406E+05 0.00052  2.12323E+05 0.00048  2.55941E+05 0.00052  2.52534E+05 0.00046  2.18654E+05 0.00044  4.09566E+05 0.00048  8.75490E+04 0.00062  1.09161E+05 0.00059  9.85233E+04 0.00070  5.75435E+04 0.00062  9.95619E+04 0.00056  6.77683E+04 0.00057  5.83871E+04 0.00065  1.13907E+04 0.00122  1.12077E+04 0.00131  1.14906E+04 0.00120  1.17429E+04 0.00112  1.16521E+04 0.00128  1.15461E+04 0.00117  1.18476E+04 0.00107  1.10924E+04 0.00121  2.08649E+04 0.00076  3.32544E+04 0.00083  4.24929E+04 0.00064  1.11476E+05 0.00052  1.15714E+05 0.00060  1.19140E+05 0.00061  7.28207E+04 0.00074  4.99925E+04 0.00079  3.67830E+04 0.00075  4.06909E+04 0.00073  7.20134E+04 0.00061  9.74202E+04 0.00065  2.36702E+05 0.00055  5.82050E+05 0.00048  1.72454E+06 0.00049  1.91955E+06 0.00050  1.78851E+06 0.00050  1.66027E+06 0.00050  1.76897E+06 0.00052  1.99912E+06 0.00053  1.98910E+06 0.00051  1.50310E+06 0.00051  1.55470E+06 0.00052  1.52506E+06 0.00050  1.43848E+06 0.00050  1.25688E+06 0.00051  8.96173E+05 0.00054  3.49214E+05 0.00052 ];

% Integral parameters:

INF_KINF                  (idx, [1:   2]) = [  0.00000E+00 0.0E+00 ];

% Flux spectra in infinite geometry:

INF_FLX                   (idx, [1:   4]) = [  5.38965E+00 0.00046  1.13659E+01 0.00046 ];
INF_FISS_FLX              (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Reaction cross sections:

INF_TOT                   (idx, [1:   4]) = [  1.01365E+00 6.6E-05  3.24133E+00 1.3E-05 ];
INF_CAPT                  (idx, [1:   4]) = [  4.71368E-04 0.00015  1.89070E-02 2.3E-05 ];
INF_ABS                   (idx, [1:   4]) = [  4.71368E-04 0.00015  1.89070E-02 2.3E-05 ];
INF_FISS                  (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_NSF                   (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_NUBAR                 (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_KAPPA                 (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_INVV                  (idx, [1:   4]) = [  8.01804E-08 0.00013  3.86690E-06 2.3E-05 ];

% Total scattering cross sections:

INF_SCATT0                (idx, [1:   4]) = [  1.01318E+00 6.6E-05  3.22243E+00 1.4E-05 ];
INF_SCATT1                (idx, [1:   4]) = [  5.94946E-01 8.3E-05  6.33504E-01 4.6E-05 ];
INF_SCATT2                (idx, [1:   4]) = [  2.27371E-01 0.00013  9.80508E-02 0.00030 ];
INF_SCATT3                (idx, [1:   4]) = [  9.29877E-03 0.00224  2.49223E-02 0.00096 ];
INF_SCATT4                (idx, [1:   4]) = [ -3.35216E-02 0.00046 -2.91791E-02 0.00059 ];
INF_SCATT5                (idx, [1:   4]) = [ -3.69104E-03 0.00369  1.34763E-02 0.00120 ];
INF_SCATT6                (idx, [1:   4]) = [  1.02404E-02 0.00120 -3.88421E-02 0.00039 ];
INF_SCATT7                (idx, [1:   4]) = [  4.09342E-04 0.02596  1.33293E-02 0.00101 ];

% Total scattering production cross sections:

INF_SCATTP0               (idx, [1:   4]) = [  1.01318E+00 6.6E-05  3.22243E+00 1.4E-05 ];
INF_SCATTP1               (idx, [1:   4]) = [  5.94946E-01 8.3E-05  6.33504E-01 4.6E-05 ];
INF_SCATTP2               (idx, [1:   4]) = [  2.27371E-01 0.00013  9.80508E-02 0.00030 ];
INF_SCATTP3               (idx, [1:   4]) = [  9.29877E-03 0.00224  2.49223E-02 0.00096 ];
INF_SCATTP4               (idx, [1:   4]) = [ -3.35216E-02 0.00046 -2.91791E-02 0.00059 ];
INF_SCATTP5               (idx, [1:   4]) = [ -3.69104E-03 0.00369  1.34763E-02 0.00120 ];
INF_SCATTP6               (idx, [1:   4]) = [  1.02404E-02 0.00120 -3.88421E-02 0.00039 ];
INF_SCATTP7               (idx, [1:   4]) = [  4.09342E-04 0.02596  1.33293E-02 0.00101 ];

% Diffusion parameters:

INF_TRANSPXS              (idx, [1:   4]) = [  2.68213E-01 0.00017  2.25923E+00 2.5E-05 ];
INF_DIFFCOEF              (idx, [1:   4]) = [  1.24280E+00 0.00017  1.47543E-01 2.5E-05 ];

% Reduced absoption and removal:

INF_RABSXS                (idx, [1:   4]) = [  4.71368E-04 0.00015  1.89070E-02 2.3E-05 ];
INF_REMXS                 (idx, [1:   4]) = [  5.66067E-02 0.00011  1.89960E-02 0.00023 ];

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

INF_S0                    (idx, [1:   8]) = [  9.57042E-01 6.5E-05  5.61369E-02 0.00011  9.15565E-05 0.00317  3.22234E+00 1.4E-05 ];
INF_S1                    (idx, [1:   8]) = [  5.78064E-01 8.2E-05  1.68822E-02 0.00031  6.36531E-05 0.00342  6.33441E-01 4.6E-05 ];
INF_S2                    (idx, [1:   8]) = [  2.32292E-01 0.00012 -4.92152E-03 0.00068  3.94979E-05 0.00473  9.80113E-02 0.00030 ];
INF_S3                    (idx, [1:   8]) = [  1.52411E-02 0.00134 -5.94232E-03 0.00062  1.93455E-05 0.00801  2.49029E-02 0.00096 ];
INF_S4                    (idx, [1:   8]) = [ -3.14594E-02 0.00048 -2.06218E-03 0.00120  5.01159E-06 0.02402 -2.91842E-02 0.00059 ];
INF_S5                    (idx, [1:   8]) = [ -3.61156E-03 0.00361 -7.94831E-05 0.03048 -2.78552E-06 0.03841  1.34791E-02 0.00120 ];
INF_S6                    (idx, [1:   8]) = [  1.08895E-02 0.00113 -6.49090E-04 0.00452 -5.66457E-06 0.01627 -3.88365E-02 0.00039 ];
INF_S7                    (idx, [1:   8]) = [  9.94813E-04 0.01057 -5.85471E-04 0.00442 -5.98188E-06 0.01468  1.33353E-02 0.00101 ];

% Scattering production matrixes:

INF_SP0                   (idx, [1:   8]) = [  9.57042E-01 6.5E-05  5.61369E-02 0.00011  9.15565E-05 0.00317  3.22234E+00 1.4E-05 ];
INF_SP1                   (idx, [1:   8]) = [  5.78064E-01 8.2E-05  1.68822E-02 0.00031  6.36531E-05 0.00342  6.33441E-01 4.6E-05 ];
INF_SP2                   (idx, [1:   8]) = [  2.32292E-01 0.00012 -4.92152E-03 0.00068  3.94979E-05 0.00473  9.80113E-02 0.00030 ];
INF_SP3                   (idx, [1:   8]) = [  1.52411E-02 0.00134 -5.94232E-03 0.00062  1.93455E-05 0.00801  2.49029E-02 0.00096 ];
INF_SP4                   (idx, [1:   8]) = [ -3.14594E-02 0.00048 -2.06218E-03 0.00120  5.01159E-06 0.02402 -2.91842E-02 0.00059 ];
INF_SP5                   (idx, [1:   8]) = [ -3.61156E-03 0.00361 -7.94831E-05 0.03048 -2.78552E-06 0.03841  1.34791E-02 0.00120 ];
INF_SP6                   (idx, [1:   8]) = [  1.08895E-02 0.00113 -6.49090E-04 0.00452 -5.66457E-06 0.01627 -3.88365E-02 0.00039 ];
INF_SP7                   (idx, [1:   8]) = [  9.94813E-04 0.01057 -5.85471E-04 0.00442 -5.98188E-06 0.01468  1.33353E-02 0.00101 ];

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

CMM_TRANSPXS              (idx, [1:   4]) = [  7.10415E-02 0.00024 -5.55342E-01 0.00090 ];
CMM_TRANSPXS_X            (idx, [1:   4]) = [  7.13497E-02 0.00029 -5.53692E-01 0.00089 ];
CMM_TRANSPXS_Y            (idx, [1:   4]) = [  7.08887E-02 0.00031 -5.56140E-01 0.00091 ];
CMM_TRANSPXS_Z            (idx, [1:   4]) = [  7.08883E-02 0.00026 -5.56202E-01 0.00093 ];
CMM_DIFFCOEF              (idx, [1:   4]) = [  4.69211E+00 0.00024 -6.00255E-01 0.00089 ];
CMM_DIFFCOEF_X            (idx, [1:   4]) = [  4.67185E+00 0.00029 -6.02042E-01 0.00089 ];
CMM_DIFFCOEF_Y            (idx, [1:   4]) = [  4.70223E+00 0.00031 -5.99393E-01 0.00091 ];
CMM_DIFFCOEF_Z            (idx, [1:   4]) = [  4.70225E+00 0.00026 -5.99329E-01 0.00093 ];

% Delayed neutron parameters (Meulekamp method):

BETA_EFF                  (idx, [1:  18]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
LAMBDA                    (idx, [1:  18]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
