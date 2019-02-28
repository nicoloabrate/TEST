
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
INPUT_FILE_NAME           (idx, [1:  7])  = 'slab_Be' ;
WORKING_DIRECTORY         (idx, [1: 32])  = '/home/abrate/SERPENT-2/phytra/Be' ;
HOSTNAME                  (idx, [1:  7])  = 'vpcen13' ;
CPU_TYPE                  (idx, [1: 41])  = 'Intel(R) Xeon(R) CPU E5-2630 v3 @ 2.40GHz' ;
CPU_MHZ                   (idx, 1)        = 4294967295.0 ;
START_DATE                (idx, [1: 24])  = 'Thu Oct 18 12:24:54 2018' ;
COMPLETE_DATE             (idx, [1: 24])  = 'Thu Oct 18 13:02:56 2018' ;

% Run parameters:

POP                       (idx, 1)        = 100000 ;
CYCLES                    (idx, 1)        = 1000 ;
SKIP                      (idx, 1)        = 500 ;
BATCH_INTERVAL            (idx, 1)        = 1 ;
SRC_NORM_MODE             (idx, 1)        = 2 ;
SEED                      (idx, 1)        = 1539858294 ;
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
OMP_HISTORY_PROFILE       (idx, [1:  30]) = [  1.07673E+00  9.95795E-01  1.01378E+00  1.01578E+00  1.03450E+00  9.98450E-01  1.01727E+00  9.77822E-01  1.01735E+00  9.81090E-01  1.00668E+00  9.99945E-01  1.02627E+00  9.68767E-01  9.53429E-01  1.03623E+00  1.01636E+00  1.02624E+00  9.57534E-01  9.98284E-01  9.59594E-01  1.03622E+00  9.68629E-01  1.01591E+00  9.57620E-01  9.66122E-01  9.58367E-01  1.02610E+00  9.63976E-01  1.02916E+00  ];
SHARE_BUF_ARRAY           (idx, 1)        = 0 ;
SHARE_RES2_ARRAY          (idx, 1)        = 1 ;

% File paths:

XS_DATA_FILE_PATH         (idx, [1: 45])  = '/opt/serpent/xsdata/jeff33/sss2_jeff33.xsdata' ;
DECAY_DATA_FILE_PATH      (idx, [1:  3])  = 'N/A' ;
SFY_DATA_FILE_PATH        (idx, [1:  3])  = 'N/A' ;
NFY_DATA_FILE_PATH        (idx, [1:  3])  = 'N/A' ;
BRA_DATA_FILE_PATH        (idx, [1:  3])  = 'N/A' ;

% Collision and reaction sampling (neutrons/photons):

MIN_MACROXS               (idx, [1:   4]) = [  5.00000E-02 3.2E-09  0.00000E+00 0.0E+00 ];
DT_THRESH                 (idx, [1:  2])  = [  9.00000E-01  9.00000E-01 ];
ST_FRAC                   (idx, [1:   4]) = [  1.98443E-03 0.00022  0.00000E+00 0.0E+00 ];
DT_FRAC                   (idx, [1:   4]) = [  9.98016E-01 4.5E-07  0.00000E+00 0.0E+00 ];
DT_EFF                    (idx, [1:   4]) = [  8.13914E-01 2.1E-05  0.00000E+00 0.0E+00 ];
REA_SAMPLING_EFF          (idx, [1:   4]) = [  1.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
REA_SAMPLING_FAIL         (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
TOT_COL_EFF               (idx, [1:   4]) = [  8.14269E-01 2.1E-05  0.00000E+00 0.0E+00 ];
AVG_TRACKING_LOOPS        (idx, [1:   8]) = [  2.37996E+00 4.4E-05  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
AVG_TRACKS                (idx, [1:   4]) = [  1.06563E+02 0.00019  0.00000E+00 0.0E+00 ];
AVG_REAL_COL              (idx, [1:   4]) = [  1.06547E+02 0.00019  0.00000E+00 0.0E+00 ];
AVG_VIRT_COL              (idx, [1:   4]) = [  2.43025E+01 0.00013  0.00000E+00 0.0E+00 ];
AVG_SURF_CROSS            (idx, [1:   4]) = [  2.53911E-02 0.00062  0.00000E+00 0.0E+00 ];
LOST_PARTICLES            (idx, 1)        = 0 ;

% Run statistics:

CYCLE_IDX                 (idx, 1)        = 1000 ;
SOURCE_POPULATION         (idx, 1)        = 100001907 ;
MEAN_POP_SIZE             (idx, [1:  2])  = [  1.00002E+05 0.00021 ];
MEAN_POP_WGT              (idx, [1:  2])  = [  1.00002E+05 0.00021 ];
SIMULATION_COMPLETED      (idx, 1)        = 1 ;

% Running times:

TOT_CPU_TIME              (idx, 1)        =  1.02221E+03 ;
RUNNING_TIME              (idx, 1)        =  3.80203E+01 ;
INIT_TIME                 (idx, [1:  2])  = [  4.71667E-02  4.71667E-02 ];
PROCESS_TIME              (idx, [1:  2])  = [  4.83333E-04  4.83333E-04 ];
TRANSPORT_CYCLE_TIME      (idx, [1:  3])  = [  3.79726E+01  3.79726E+01  0.00000E+00 ];
MPI_OVERHEAD_TIME         (idx, [1:  2])  = [  0.00000E+00  0.00000E+00 ];
ESTIMATED_RUNNING_TIME    (idx, [1:  2])  = [  3.80187E+01  0.00000E+00 ];
CPU_USAGE                 (idx, 1)        = 26.88582 ;
TRANSPORT_CPU_USAGE       (idx, [1:   2]) = [  2.73107E+01 0.00022 ];
OMP_PARALLEL_FRAC         (idx, 1)        =  8.91196E-01 ;

% Memory usage:

AVAIL_MEM                 (idx, 1)        = 96872.24 ;
ALLOC_MEMSIZE             (idx, 1)        = 1094.02;
MEMSIZE                   (idx, 1)        = 869.75;
XS_MEMSIZE                (idx, 1)        = 194.50;
MAT_MEMSIZE               (idx, 1)        = 18.64;
RES_MEMSIZE               (idx, 1)        = 3.37;
MISC_MEMSIZE              (idx, 1)        = 653.23;
UNKNOWN_MEMSIZE           (idx, 1)        = 0.00;
UNUSED_MEMSIZE            (idx, 1)        = 224.27;

% Geometry parameters:

TOT_CELLS                 (idx, 1)        = 7 ;
UNION_CELLS               (idx, 1)        = 0 ;

% Neutron energy grid:

NEUTRON_ERG_TOL           (idx, 1)        =  0.00000E+00 ;
NEUTRON_ERG_NE            (idx, 1)        = 203000 ;
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
TOT_REA_CHANNELS          (idx, 1)        = 177 ;
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

NORM_COEF                 (idx, [1:   4]) = [  9.72815E-06 9.6E-05  0.00000E+00 0.0E+00 ];

% Analog reaction rate estimators:

CONVERSION_RATIO          (idx, [1:   2]) = [  1.66982E+00 0.00024 ];
U235_FISS                 (idx, [1:   4]) = [  2.44356E-01 0.00019  7.82063E-01 9.2E-05 ];
U238_FISS                 (idx, [1:   4]) = [  6.80951E-02 0.00039  2.17937E-01 0.00033 ];
U235_CAPT                 (idx, [1:   4]) = [  6.72895E-02 0.00039  1.00411E-01 0.00038 ];
U238_CAPT                 (idx, [1:   4]) = [  5.20510E-01 0.00014  7.76716E-01 7.3E-05 ];

% Neutron balance (particles/weight):

BALA_SRC_NEUTRON_SRC     (idx, [1:  2])  = [ 0 0.00000E+00 ];
BALA_SRC_NEUTRON_FISS    (idx, [1:  2])  = [ 100001907 1.00000E+08 ];
BALA_SRC_NEUTRON_NXN     (idx, [1:  2])  = [ 0 2.79086E+06 ];
BALA_SRC_NEUTRON_VR      (idx, [1:  2])  = [ 0 0.00000E+00 ];
BALA_SRC_NEUTRON_TOT     (idx, [1:  2])  = [ 100001907 1.02791E+08 ];

BALA_LOSS_NEUTRON_CAPT    (idx, [1:  2])  = [ 67178403 6.88868E+07 ];
BALA_LOSS_NEUTRON_FISS    (idx, [1:  2])  = [ 31171345 3.21182E+07 ];
BALA_LOSS_NEUTRON_LEAK    (idx, [1:  2])  = [ 1652159 1.78581E+06 ];
BALA_LOSS_NEUTRON_CUT     (idx, [1:  2])  = [ 0 0.00000E+00 ];
BALA_LOSS_NEUTRON_ERR     (idx, [1:  2])  = [ 0 0.00000E+00 ];
BALA_LOSS_NEUTRON_TOT     (idx, [1:  2])  = [ 100001907 1.02791E+08 ];

BALA_NEUTRON_DIFF         (idx, [1:  2])  = [ 0 -5.76675E-06 ];

% Normalized total reaction rates (neutrons):

TOT_POWER                 (idx, [1:   2]) = [  1.01793E-11 8.2E-05 ];
TOT_POWDENS               (idx, [1:   2]) = [  0.00000E+00 0.0E+00 ];
TOT_GENRATE               (idx, [1:   2]) = [  7.84499E-01 8.2E-05 ];
TOT_FISSRATE              (idx, [1:   2]) = [  3.12479E-01 8.3E-05 ];
TOT_CAPTRATE              (idx, [1:   2]) = [  6.70149E-01 3.8E-05 ];
TOT_ABSRATE               (idx, [1:   2]) = [  9.82628E-01 1.4E-05 ];
TOT_SRCRATE               (idx, [1:   2]) = [  9.72815E-01 9.6E-05 ];
TOT_FLUX                  (idx, [1:   2]) = [  1.84101E+02 0.00014 ];
TOT_PHOTON_PRODRATE       (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
TOT_LEAKRATE              (idx, [1:   2]) = [  1.73724E-02 0.00078 ];
ALBEDO_LEAKRATE           (idx, [1:   2]) = [  0.00000E+00 0.0E+00 ];
TOT_LOSSRATE              (idx, [1:   2]) = [  1.00000E+00 0.0E+00 ];
TOT_CUTRATE               (idx, [1:   2]) = [  0.00000E+00 0.0E+00 ];
TOT_RR                    (idx, [1:   2]) = [  1.09069E+02 0.00018 ];
INI_FMASS                 (idx, 1)        =  0.00000E+00 ;
TOT_FMASS                 (idx, 1)        =  0.00000E+00 ;

% Six-factor formula:

SIX_FF_ETA                (idx, [1:   2]) = [  1.82531E+00 0.00016 ];
SIX_FF_F                  (idx, [1:   2]) = [  6.65608E-01 0.00016 ];
SIX_FF_P                  (idx, [1:   2]) = [  2.17175E-01 0.00024 ];
SIX_FF_EPSILON            (idx, [1:   2]) = [  3.11187E+00 0.00028 ];
SIX_FF_LF                 (idx, [1:   2]) = [  9.99963E-01 7.4E-07 ];
SIX_FF_LT                 (idx, [1:   2]) = [  9.82178E-01 1.4E-05 ];
SIX_FF_KINF               (idx, [1:   2]) = [  8.21016E-01 0.00015 ];
SIX_FF_KEFF               (idx, [1:   2]) = [  8.06354E-01 0.00015 ];

% Fission neutron and energy production:

NUBAR                     (idx, [1:   2]) = [  2.51057E+00 7.9E-06 ];
FISSE                     (idx, [1:   2]) = [  2.03324E+02 9.8E-07 ];

% Criticality eigenvalues:

ANA_KEFF                  (idx, [1:   6]) = [  8.06361E-01 0.00015  8.00760E-01 0.00015  5.59335E-03 0.00211 ];
IMP_KEFF                  (idx, [1:   2]) = [  8.06390E-01 8.5E-05 ];
COL_KEFF                  (idx, [1:   2]) = [  8.06428E-01 0.00012 ];
ABS_KEFF                  (idx, [1:   2]) = [  8.06390E-01 8.5E-05 ];
ABS_KINF                  (idx, [1:   2]) = [  8.21052E-01 8.3E-05 ];
GEOM_ALBEDO               (idx, [1:   6]) = [  1.00000E+00 0.0E+00  1.00000E+00 0.0E+00  1.00000E+00 0.0E+00 ];

% ALF (Average lethargy of neutrons causing fission):
% Based on E0 = 2.000000E+01 MeV

ANA_ALF                   (idx, [1:   2]) = [  1.06539E+01 0.00013 ];
IMP_ALF                   (idx, [1:   2]) = [  1.06552E+01 0.00011 ];

% EALF (Energy corresponding to average lethargy of neutrons causing fission):

ANA_EALF                  (idx, [1:   2]) = [  4.72626E-04 0.00142 ];
IMP_EALF                  (idx, [1:   2]) = [  4.71881E-04 0.00121 ];

% AFGE (Average energy of neutrons causing fission):

ANA_AFGE                  (idx, [1:   2]) = [  7.30708E-01 0.00036 ];
IMP_AFGE                  (idx, [1:   2]) = [  7.30553E-01 0.00022 ];

% Forward-weighted delayed neutron parameters:

FWD_ANA_BETA_ZERO         (idx, [1:  18]) = [  1.17846E-02 0.00114  2.72029E-04 0.00710  1.54476E-03 0.00297  8.39573E-04 0.00411  2.02367E-03 0.00272  3.65053E-03 0.00204  1.61538E-03 0.00284  1.21189E-03 0.00345  6.26767E-04 0.00448 ];
FWD_ANA_LAMBDA            (idx, [1:  18]) = [  5.69057E-01 0.00158  1.24667E-02 7.2E-09  2.82917E-02 0.0E+00  4.25244E-02 0.0E+00  1.33042E-01 0.0E+00  2.92467E-01 0.0E+00  6.66488E-01 2.7E-09  1.63478E+00 4.2E-09  3.55460E+00 3.0E-09 ];

% Beta-eff using Meulekamp's method:

ADJ_MEULEKAMP_BETA_EFF    (idx, [1:  18]) = [  7.16556E-03 0.00177  1.57627E-04 0.01159  9.53016E-04 0.00498  5.03394E-04 0.00685  1.22501E-03 0.00447  2.22433E-03 0.00319  9.77719E-04 0.00477  7.44218E-04 0.00560  3.80245E-04 0.00769 ];
ADJ_MEULEKAMP_LAMBDA      (idx, [1:  18]) = [  5.70082E-01 0.00263  1.24667E-02 7.2E-09  2.82917E-02 0.0E+00  4.25244E-02 0.0E+00  1.33042E-01 0.0E+00  2.92467E-01 0.0E+00  6.66488E-01 2.7E-09  1.63478E+00 4.2E-09  3.55460E+00 2.9E-09 ];

% Adjoint weighted time constants using Nauchi's method:

ADJ_NAUCHI_GEN_TIME       (idx, [1:   6]) = [  3.77377E-04 0.00073  3.77240E-04 0.00074  3.96394E-04 0.00786 ];
ADJ_NAUCHI_LIFETIME       (idx, [1:   6]) = [  3.04295E-04 0.00072  3.04185E-04 0.00072  3.19645E-04 0.00786 ];
ADJ_NAUCHI_BETA_EFF       (idx, [1:  18]) = [  6.94294E-03 0.00212  1.54986E-04 0.01376  9.23598E-04 0.00574  4.85386E-04 0.00802  1.19133E-03 0.00506  2.15075E-03 0.00381  9.47914E-04 0.00563  7.15905E-04 0.00662  3.73073E-04 0.00884 ];
ADJ_NAUCHI_LAMBDA         (idx, [1:  18]) = [  5.71100E-01 0.00309  1.24667E-02 7.1E-09  2.82917E-02 0.0E+00  4.25244E-02 0.0E+00  1.33042E-01 0.0E+00  2.92467E-01 0.0E+00  6.66488E-01 2.7E-09  1.63478E+00 4.2E-09  3.55460E+00 2.9E-09 ];

% Adjoint weighted time constants using IFP:

ADJ_IFP_GEN_TIME          (idx, [1:   6]) = [  4.24621E-04 0.00176  4.24325E-04 0.00177  4.66854E-04 0.01958 ];
ADJ_IFP_LIFETIME          (idx, [1:   6]) = [  3.42385E-04 0.00175  3.42147E-04 0.00176  3.76471E-04 0.01960 ];
ADJ_IFP_IMP_BETA_EFF      (idx, [1:  18]) = [  7.00216E-03 0.00756  1.58637E-04 0.05008  9.58060E-04 0.02055  4.70082E-04 0.02863  1.21467E-03 0.01856  2.20366E-03 0.01348  9.50740E-04 0.01991  7.14426E-04 0.02303  3.31888E-04 0.03320 ];
ADJ_IFP_IMP_LAMBDA        (idx, [1:  18]) = [  5.47820E-01 0.01090  1.24667E-02 0.0E+00  2.82917E-02 0.0E+00  4.25244E-02 0.0E+00  1.33042E-01 0.0E+00  2.92467E-01 0.0E+00  6.66488E-01 2.3E-09  1.63478E+00 4.5E-09  3.55460E+00 2.6E-09 ];
ADJ_IFP_ANA_BETA_EFF      (idx, [1:  18]) = [  6.99668E-03 0.00739  1.56935E-04 0.04907  9.54247E-04 0.02006  4.69135E-04 0.02793  1.21494E-03 0.01825  2.20132E-03 0.01311  9.52459E-04 0.01952  7.15454E-04 0.02231  3.32186E-04 0.03258 ];
ADJ_IFP_ANA_LAMBDA        (idx, [1:  18]) = [  5.48448E-01 0.01063  1.24667E-02 0.0E+00  2.82917E-02 0.0E+00  4.25244E-02 0.0E+00  1.33042E-01 0.0E+00  2.92467E-01 0.0E+00  6.66488E-01 2.3E-09  1.63478E+00 4.5E-09  3.55460E+00 2.8E-09 ];
ADJ_IFP_ROSSI_ALPHA       (idx, [1:   2]) = [ -1.65502E+01 0.00769 ];

% Adjoint weighted time constants using perturbation technique:

ADJ_PERT_GEN_TIME         (idx, [1:   2]) = [  4.04844E-04 0.00041 ];
ADJ_PERT_LIFETIME         (idx, [1:   2]) = [  3.26443E-04 0.00039 ];
ADJ_PERT_BETA_EFF         (idx, [1:   2]) = [  7.00569E-03 0.00128 ];
ADJ_PERT_ROSSI_ALPHA      (idx, [1:   2]) = [ -1.73067E+01 0.00130 ];

% Inverse neutron speed :

ANA_INV_SPD               (idx, [1:   2]) = [  1.60556E-06 0.00021 ];

% Analog slowing-down and thermal neutron lifetime (total/prompt/delayed):

ANA_SLOW_TIME             (idx, [1:   6]) = [  1.20639E-05 7.2E-05  1.20635E-05 7.3E-05  1.21154E-05 0.00082 ];
ANA_THERM_TIME            (idx, [1:   6]) = [  1.29914E-03 0.00028  1.29980E-03 0.00028  1.21354E-03 0.00332 ];
ANA_THERM_FRAC            (idx, [1:   6]) = [  2.31208E-01 0.00023  2.32175E-01 0.00023  1.50302E-01 0.00272 ];
ANA_DELAYED_EMTIME        (idx, [1:   2]) = [  1.07279E+01 0.00241 ];
ANA_MEAN_NCOL             (idx, [1:   4]) = [  1.06547E+02 0.00019  1.02508E+02 0.00031 ];

% Group constant generation:

GC_UNIVERSE_NAME          (idx, [1:  1])  = 'UO2' ;

% Micro- and macro-group structures:

MICRO_NG                  (idx, 1)        = 70 ;
MICRO_E                   (idx, [1:  71]) = [  1.00000E-11  5.00000E-09  1.00000E-08  1.50000E-08  2.00000E-08  2.50000E-08  3.00000E-08  3.50000E-08  4.20000E-08  5.00000E-08  5.80000E-08  6.70000E-08  8.00000E-08  1.00000E-07  1.40000E-07  1.80000E-07  2.20000E-07  2.50000E-07  2.80000E-07  3.00000E-07  3.20000E-07  3.50000E-07  4.00000E-07  5.00000E-07  6.25000E-07  7.80000E-07  8.50000E-07  9.10000E-07  9.50000E-07  9.72000E-07  9.96000E-07  1.02000E-06  1.04500E-06  1.07100E-06  1.09700E-06  1.12300E-06  1.15000E-06  1.30000E-06  1.50000E-06  1.85500E-06  2.10000E-06  2.60000E-06  3.30000E-06  4.00000E-06  9.87700E-06  1.59680E-05  2.77000E-05  4.80520E-05  7.55014E-05  1.48728E-04  3.67262E-04  9.06898E-04  1.42510E-03  2.23945E-03  3.51910E-03  5.50000E-03  9.11800E-03  1.50300E-02  2.47800E-02  4.08500E-02  6.74300E-02  1.11000E-01  1.83000E-01  3.02500E-01  5.00000E-01  8.21000E-01  1.35300E+00  2.23100E+00  3.67900E+00  6.06550E+00  2.00000E+01 ];

MACRO_NG                  (idx, 1)        = 2 ;
MACRO_E                   (idx, [1:   3]) = [  1.00000E+37  6.25000E-07  0.00000E+00 ];

% Micro-group spectrum:

INF_MICRO_FLX             (idx, [1: 140]) = [  4.27803E+05 0.00075  1.95159E+06 0.00038  4.64727E+06 0.00029  5.34537E+06 0.00024  6.29962E+06 0.00020  1.44415E+07 0.00021  1.16258E+07 0.00015  1.64775E+07 0.00016  1.65927E+07 0.00016  1.55813E+07 0.00018  1.44972E+07 0.00019  1.27906E+07 0.00018  1.04961E+07 0.00018  8.22948E+06 0.00022  6.28040E+06 0.00028  4.16653E+06 0.00030  3.16216E+06 0.00032  2.43257E+06 0.00041  1.80317E+06 0.00045  2.23325E+06 0.00054  1.20608E+06 0.00065  5.12937E+05 0.00080  2.33085E+05 0.00118  1.93749E+05 0.00117  1.42327E+05 0.00115  1.50689E+05 0.00148  1.70752E+05 0.00134  5.24724E+04 0.00205  7.44307E+04 0.00224  7.99253E+04 0.00177  4.57326E+04 0.00198  8.72395E+04 0.00182  6.02354E+04 0.00193  4.48185E+04 0.00196  6.92005E+03 0.00305  6.71987E+03 0.00354  6.93043E+03 0.00317  7.24576E+03 0.00323  7.20155E+03 0.00360  7.14971E+03 0.00372  7.27410E+03 0.00349  6.90660E+03 0.00327  1.29296E+04 0.00310  2.02983E+04 0.00269  2.52316E+04 0.00255  6.26377E+04 0.00182  5.73681E+04 0.00178  4.87953E+04 0.00164  2.40355E+04 0.00234  1.38227E+04 0.00210  9.02350E+03 0.00327  9.04569E+03 0.00297  1.45265E+04 0.00218  1.71244E+04 0.00250  3.10634E+04 0.00171  5.25216E+04 0.00147  1.03911E+05 0.00080  8.75979E+04 0.00094  7.40268E+04 0.00088  6.01475E+04 0.00084  5.69088E+04 0.00102  5.92058E+04 0.00084  5.30388E+04 0.00085  3.73647E+04 0.00113  3.04980E+04 0.00111  2.80099E+04 0.00099  2.46412E+04 0.00113  1.69503E+04 0.00130  9.26736E+03 0.00133  6.83439E+03 0.00147 ];

% Integral parameters:

INF_KINF                  (idx, [1:   2]) = [  8.70761E-01 0.00013 ];

% Flux spectra in infinite geometry:

INF_FLX                   (idx, [1:   4]) = [  7.91452E+01 0.00013  4.50281E-01 0.00043 ];
INF_FISS_FLX              (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Reaction cross sections:

INF_TOT                   (idx, [1:   4]) = [  4.19690E-01 1.1E-05  6.93798E-01 9.9E-05 ];
INF_CAPT                  (idx, [1:   4]) = [  7.03695E-03 8.4E-05  7.58453E-02 0.00021 ];
INF_ABS                   (idx, [1:   4]) = [  9.67147E-03 8.2E-05  3.06750E-01 0.00022 ];
INF_FISS                  (idx, [1:   4]) = [  2.63452E-03 0.00011  2.30904E-01 0.00022 ];
INF_NSF                   (idx, [1:   4]) = [  6.72603E-03 0.00010  5.60032E-01 0.00022 ];
INF_NUBAR                 (idx, [1:   4]) = [  2.55304E+00 1.0E-05  2.42539E+00 0.0E+00 ];
INF_KAPPA                 (idx, [1:   4]) = [  2.03850E+02 1.2E-06  2.02270E+02 0.0E+00 ];
INF_INVV                  (idx, [1:   4]) = [  7.57898E-09 0.00042  2.78133E-06 0.00020 ];

% Total scattering cross sections:

INF_SCATT0                (idx, [1:   4]) = [  4.10018E-01 1.2E-05  3.87072E-01 0.00015 ];
INF_SCATT1                (idx, [1:   4]) = [  3.54860E-02 0.00013  7.79194E-03 0.00734 ];
INF_SCATT2                (idx, [1:   4]) = [  1.54527E-02 0.00022  3.04173E-04 0.16058 ];
INF_SCATT3                (idx, [1:   4]) = [  5.16246E-03 0.00049  4.04896E-05 0.89869 ];
INF_SCATT4                (idx, [1:   4]) = [  2.97284E-03 0.00086  4.99554E-05 0.68322 ];
INF_SCATT5                (idx, [1:   4]) = [  1.38246E-03 0.00159  1.23899E-05 1.00000 ];
INF_SCATT6                (idx, [1:   4]) = [  6.95515E-04 0.00281 -2.43276E-05 1.00000 ];
INF_SCATT7                (idx, [1:   4]) = [  2.84457E-04 0.00652 -8.99320E-06 1.00000 ];

% Total scattering production cross sections:

INF_SCATTP0               (idx, [1:   4]) = [  4.10051E-01 1.2E-05  3.87072E-01 0.00015 ];
INF_SCATTP1               (idx, [1:   4]) = [  3.54865E-02 0.00013  7.79194E-03 0.00734 ];
INF_SCATTP2               (idx, [1:   4]) = [  1.54528E-02 0.00022  3.04173E-04 0.16058 ];
INF_SCATTP3               (idx, [1:   4]) = [  5.16247E-03 0.00049  4.04896E-05 0.89869 ];
INF_SCATTP4               (idx, [1:   4]) = [  2.97283E-03 0.00086  4.99554E-05 0.68322 ];
INF_SCATTP5               (idx, [1:   4]) = [  1.38247E-03 0.00159  1.23899E-05 1.00000 ];
INF_SCATTP6               (idx, [1:   4]) = [  6.95500E-04 0.00280 -2.43276E-05 1.00000 ];
INF_SCATTP7               (idx, [1:   4]) = [  2.84444E-04 0.00652 -8.99320E-06 1.00000 ];

% Diffusion parameters:

INF_TRANSPXS              (idx, [1:   4]) = [  3.48031E-01 2.8E-05  6.43825E-01 0.00016 ];
INF_DIFFCOEF              (idx, [1:   4]) = [  9.57770E-01 2.8E-05  5.17740E-01 0.00016 ];

% Reduced absoption and removal:

INF_RABSXS                (idx, [1:   4]) = [  9.63857E-03 8.3E-05  3.06750E-01 0.00022 ];
INF_REMXS                 (idx, [1:   4]) = [  9.71171E-03 7.4E-05  3.07431E-01 0.00030 ];

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

INF_CHIT                  (idx, [1:   4]) = [  1.00000E+00 9.2E-09  9.32927E-09 1.00000 ];
INF_CHIP                  (idx, [1:   4]) = [  1.00000E+00 9.6E-09  9.41375E-09 1.00000 ];
INF_CHID                  (idx, [1:   4]) = [  1.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Scattering matrixes:

INF_S0                    (idx, [1:   8]) = [  4.09978E-01 1.2E-05  4.04464E-05 0.00175  7.05990E-04 0.00689  3.86366E-01 0.00015 ];
INF_S1                    (idx, [1:   8]) = [  3.54972E-02 0.00013 -1.12223E-05 0.00366  9.37947E-06 0.21838  7.78256E-03 0.00739 ];
INF_S2                    (idx, [1:   8]) = [  1.54533E-02 0.00022 -6.59956E-07 0.04858 -2.28228E-05 0.08989  3.26995E-04 0.14965 ];
INF_S3                    (idx, [1:   8]) = [  5.16254E-03 0.00049 -8.38279E-08 0.33613 -1.79404E-05 0.09078  5.84300E-05 0.61592 ];
INF_S4                    (idx, [1:   8]) = [  2.97288E-03 0.00086 -4.21496E-08 0.62999 -1.01382E-05 0.13827  6.00936E-05 0.56515 ];
INF_S5                    (idx, [1:   8]) = [  1.38254E-03 0.00159 -7.92023E-08 0.24540 -5.86215E-06 0.21949  1.82521E-05 1.00000 ];
INF_S6                    (idx, [1:   8]) = [  6.95498E-04 0.00281  1.67481E-08 1.00000 -3.84861E-06 0.22075 -2.04790E-05 1.00000 ];
INF_S7                    (idx, [1:   8]) = [  2.84472E-04 0.00652 -1.48212E-08 1.00000 -1.24204E-06 0.67474 -7.75115E-06 1.00000 ];

% Scattering production matrixes:

INF_SP0                   (idx, [1:   8]) = [  4.10011E-01 1.2E-05  4.04464E-05 0.00175  7.05990E-04 0.00689  3.86366E-01 0.00015 ];
INF_SP1                   (idx, [1:   8]) = [  3.54977E-02 0.00013 -1.12223E-05 0.00366  9.37947E-06 0.21838  7.78256E-03 0.00739 ];
INF_SP2                   (idx, [1:   8]) = [  1.54535E-02 0.00022 -6.59956E-07 0.04858 -2.28228E-05 0.08989  3.26995E-04 0.14965 ];
INF_SP3                   (idx, [1:   8]) = [  5.16255E-03 0.00049 -8.38279E-08 0.33613 -1.79404E-05 0.09078  5.84300E-05 0.61592 ];
INF_SP4                   (idx, [1:   8]) = [  2.97287E-03 0.00086 -4.21496E-08 0.62999 -1.01382E-05 0.13827  6.00936E-05 0.56515 ];
INF_SP5                   (idx, [1:   8]) = [  1.38255E-03 0.00159 -7.92023E-08 0.24540 -5.86215E-06 0.21949  1.82521E-05 1.00000 ];
INF_SP6                   (idx, [1:   8]) = [  6.95483E-04 0.00281  1.67481E-08 1.00000 -3.84861E-06 0.22075 -2.04790E-05 1.00000 ];
INF_SP7                   (idx, [1:   8]) = [  2.84459E-04 0.00651 -1.48212E-08 1.00000 -1.24204E-06 0.67474 -7.75115E-06 1.00000 ];

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

CMM_TRANSPXS              (idx, [1:   4]) = [  3.78181E-01 0.00016  7.64168E-03 0.00062 ];
CMM_TRANSPXS_X            (idx, [1:   4]) = [  4.32573E-01 0.00022  5.58566E-02 0.00132 ];
CMM_TRANSPXS_Y            (idx, [1:   4]) = [  3.55841E-01 0.00017  5.34251E-03 0.00069 ];
CMM_TRANSPXS_Z            (idx, [1:   4]) = [  3.55781E-01 0.00023  5.33334E-03 0.00074 ];
CMM_DIFFCOEF              (idx, [1:   4]) = [  8.81414E-01 0.00016  4.36213E+01 0.00062 ];
CMM_DIFFCOEF_X            (idx, [1:   4]) = [  7.70585E-01 0.00022  5.96817E+00 0.00132 ];
CMM_DIFFCOEF_Y            (idx, [1:   4]) = [  9.36748E-01 0.00017  6.23940E+01 0.00069 ];
CMM_DIFFCOEF_Z            (idx, [1:   4]) = [  9.36908E-01 0.00023  6.25016E+01 0.00074 ];

% Delayed neutron parameters (Meulekamp method):

BETA_EFF                  (idx, [1:  18]) = [  7.16556E-03 0.00177  1.57627E-04 0.01159  9.53016E-04 0.00498  5.03394E-04 0.00685  1.22501E-03 0.00447  2.22433E-03 0.00319  9.77719E-04 0.00477  7.44218E-04 0.00560  3.80245E-04 0.00769 ];
LAMBDA                    (idx, [1:  18]) = [  5.70082E-01 0.00263  1.24667E-02 7.2E-09  2.82917E-02 0.0E+00  4.25244E-02 0.0E+00  1.33042E-01 0.0E+00  2.92467E-01 0.0E+00  6.66488E-01 2.7E-09  1.63478E+00 4.2E-09  3.55460E+00 2.9E-09 ];


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
INPUT_FILE_NAME           (idx, [1:  7])  = 'slab_Be' ;
WORKING_DIRECTORY         (idx, [1: 32])  = '/home/abrate/SERPENT-2/phytra/Be' ;
HOSTNAME                  (idx, [1:  7])  = 'vpcen13' ;
CPU_TYPE                  (idx, [1: 41])  = 'Intel(R) Xeon(R) CPU E5-2630 v3 @ 2.40GHz' ;
CPU_MHZ                   (idx, 1)        = 4294967295.0 ;
START_DATE                (idx, [1: 24])  = 'Thu Oct 18 12:24:54 2018' ;
COMPLETE_DATE             (idx, [1: 24])  = 'Thu Oct 18 13:02:56 2018' ;

% Run parameters:

POP                       (idx, 1)        = 100000 ;
CYCLES                    (idx, 1)        = 1000 ;
SKIP                      (idx, 1)        = 500 ;
BATCH_INTERVAL            (idx, 1)        = 1 ;
SRC_NORM_MODE             (idx, 1)        = 2 ;
SEED                      (idx, 1)        = 1539858294 ;
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
OMP_HISTORY_PROFILE       (idx, [1:  30]) = [  1.07673E+00  9.95795E-01  1.01378E+00  1.01578E+00  1.03450E+00  9.98450E-01  1.01727E+00  9.77822E-01  1.01735E+00  9.81090E-01  1.00668E+00  9.99945E-01  1.02627E+00  9.68767E-01  9.53429E-01  1.03623E+00  1.01636E+00  1.02624E+00  9.57534E-01  9.98284E-01  9.59594E-01  1.03622E+00  9.68629E-01  1.01591E+00  9.57620E-01  9.66122E-01  9.58367E-01  1.02610E+00  9.63976E-01  1.02916E+00  ];
SHARE_BUF_ARRAY           (idx, 1)        = 0 ;
SHARE_RES2_ARRAY          (idx, 1)        = 1 ;

% File paths:

XS_DATA_FILE_PATH         (idx, [1: 45])  = '/opt/serpent/xsdata/jeff33/sss2_jeff33.xsdata' ;
DECAY_DATA_FILE_PATH      (idx, [1:  3])  = 'N/A' ;
SFY_DATA_FILE_PATH        (idx, [1:  3])  = 'N/A' ;
NFY_DATA_FILE_PATH        (idx, [1:  3])  = 'N/A' ;
BRA_DATA_FILE_PATH        (idx, [1:  3])  = 'N/A' ;

% Collision and reaction sampling (neutrons/photons):

MIN_MACROXS               (idx, [1:   4]) = [  5.00000E-02 3.2E-09  0.00000E+00 0.0E+00 ];
DT_THRESH                 (idx, [1:  2])  = [  9.00000E-01  9.00000E-01 ];
ST_FRAC                   (idx, [1:   4]) = [  1.98443E-03 0.00022  0.00000E+00 0.0E+00 ];
DT_FRAC                   (idx, [1:   4]) = [  9.98016E-01 4.5E-07  0.00000E+00 0.0E+00 ];
DT_EFF                    (idx, [1:   4]) = [  8.13914E-01 2.1E-05  0.00000E+00 0.0E+00 ];
REA_SAMPLING_EFF          (idx, [1:   4]) = [  1.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
REA_SAMPLING_FAIL         (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
TOT_COL_EFF               (idx, [1:   4]) = [  8.14269E-01 2.1E-05  0.00000E+00 0.0E+00 ];
AVG_TRACKING_LOOPS        (idx, [1:   8]) = [  2.37996E+00 4.4E-05  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
AVG_TRACKS                (idx, [1:   4]) = [  1.06563E+02 0.00019  0.00000E+00 0.0E+00 ];
AVG_REAL_COL              (idx, [1:   4]) = [  1.06547E+02 0.00019  0.00000E+00 0.0E+00 ];
AVG_VIRT_COL              (idx, [1:   4]) = [  2.43025E+01 0.00013  0.00000E+00 0.0E+00 ];
AVG_SURF_CROSS            (idx, [1:   4]) = [  2.53911E-02 0.00062  0.00000E+00 0.0E+00 ];
LOST_PARTICLES            (idx, 1)        = 0 ;

% Run statistics:

CYCLE_IDX                 (idx, 1)        = 1000 ;
SOURCE_POPULATION         (idx, 1)        = 100001907 ;
MEAN_POP_SIZE             (idx, [1:  2])  = [  1.00002E+05 0.00021 ];
MEAN_POP_WGT              (idx, [1:  2])  = [  1.00002E+05 0.00021 ];
SIMULATION_COMPLETED      (idx, 1)        = 1 ;

% Running times:

TOT_CPU_TIME              (idx, 1)        =  1.02221E+03 ;
RUNNING_TIME              (idx, 1)        =  3.80203E+01 ;
INIT_TIME                 (idx, [1:  2])  = [  4.71667E-02  4.71667E-02 ];
PROCESS_TIME              (idx, [1:  2])  = [  4.83333E-04  4.83333E-04 ];
TRANSPORT_CYCLE_TIME      (idx, [1:  3])  = [  3.79726E+01  3.79726E+01  0.00000E+00 ];
MPI_OVERHEAD_TIME         (idx, [1:  2])  = [  0.00000E+00  0.00000E+00 ];
ESTIMATED_RUNNING_TIME    (idx, [1:  2])  = [  3.80187E+01  0.00000E+00 ];
CPU_USAGE                 (idx, 1)        = 26.88581 ;
TRANSPORT_CPU_USAGE       (idx, [1:   2]) = [  2.73107E+01 0.00022 ];
OMP_PARALLEL_FRAC         (idx, 1)        =  8.91196E-01 ;

% Memory usage:

AVAIL_MEM                 (idx, 1)        = 96872.24 ;
ALLOC_MEMSIZE             (idx, 1)        = 1094.02;
MEMSIZE                   (idx, 1)        = 869.75;
XS_MEMSIZE                (idx, 1)        = 194.50;
MAT_MEMSIZE               (idx, 1)        = 18.64;
RES_MEMSIZE               (idx, 1)        = 3.37;
MISC_MEMSIZE              (idx, 1)        = 653.23;
UNKNOWN_MEMSIZE           (idx, 1)        = 0.00;
UNUSED_MEMSIZE            (idx, 1)        = 224.27;

% Geometry parameters:

TOT_CELLS                 (idx, 1)        = 7 ;
UNION_CELLS               (idx, 1)        = 0 ;

% Neutron energy grid:

NEUTRON_ERG_TOL           (idx, 1)        =  0.00000E+00 ;
NEUTRON_ERG_NE            (idx, 1)        = 203000 ;
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
TOT_REA_CHANNELS          (idx, 1)        = 177 ;
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

NORM_COEF                 (idx, [1:   4]) = [  9.72815E-06 9.6E-05  0.00000E+00 0.0E+00 ];

% Analog reaction rate estimators:

CONVERSION_RATIO          (idx, [1:   2]) = [  1.66982E+00 0.00024 ];
U235_FISS                 (idx, [1:   4]) = [  2.44356E-01 0.00019  7.82063E-01 9.2E-05 ];
U238_FISS                 (idx, [1:   4]) = [  6.80951E-02 0.00039  2.17937E-01 0.00033 ];
U235_CAPT                 (idx, [1:   4]) = [  6.72895E-02 0.00039  1.00411E-01 0.00038 ];
U238_CAPT                 (idx, [1:   4]) = [  5.20510E-01 0.00014  7.76716E-01 7.3E-05 ];

% Neutron balance (particles/weight):

BALA_SRC_NEUTRON_SRC     (idx, [1:  2])  = [ 0 0.00000E+00 ];
BALA_SRC_NEUTRON_FISS    (idx, [1:  2])  = [ 100001907 1.00000E+08 ];
BALA_SRC_NEUTRON_NXN     (idx, [1:  2])  = [ 0 2.79086E+06 ];
BALA_SRC_NEUTRON_VR      (idx, [1:  2])  = [ 0 0.00000E+00 ];
BALA_SRC_NEUTRON_TOT     (idx, [1:  2])  = [ 100001907 1.02791E+08 ];

BALA_LOSS_NEUTRON_CAPT    (idx, [1:  2])  = [ 67178403 6.88868E+07 ];
BALA_LOSS_NEUTRON_FISS    (idx, [1:  2])  = [ 31171345 3.21182E+07 ];
BALA_LOSS_NEUTRON_LEAK    (idx, [1:  2])  = [ 1652159 1.78581E+06 ];
BALA_LOSS_NEUTRON_CUT     (idx, [1:  2])  = [ 0 0.00000E+00 ];
BALA_LOSS_NEUTRON_ERR     (idx, [1:  2])  = [ 0 0.00000E+00 ];
BALA_LOSS_NEUTRON_TOT     (idx, [1:  2])  = [ 100001907 1.02791E+08 ];

BALA_NEUTRON_DIFF         (idx, [1:  2])  = [ 0 -5.76675E-06 ];

% Normalized total reaction rates (neutrons):

TOT_POWER                 (idx, [1:   2]) = [  1.01793E-11 8.2E-05 ];
TOT_POWDENS               (idx, [1:   2]) = [  0.00000E+00 0.0E+00 ];
TOT_GENRATE               (idx, [1:   2]) = [  7.84499E-01 8.2E-05 ];
TOT_FISSRATE              (idx, [1:   2]) = [  3.12479E-01 8.3E-05 ];
TOT_CAPTRATE              (idx, [1:   2]) = [  6.70149E-01 3.8E-05 ];
TOT_ABSRATE               (idx, [1:   2]) = [  9.82628E-01 1.4E-05 ];
TOT_SRCRATE               (idx, [1:   2]) = [  9.72815E-01 9.6E-05 ];
TOT_FLUX                  (idx, [1:   2]) = [  1.84101E+02 0.00014 ];
TOT_PHOTON_PRODRATE       (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
TOT_LEAKRATE              (idx, [1:   2]) = [  1.73724E-02 0.00078 ];
ALBEDO_LEAKRATE           (idx, [1:   2]) = [  0.00000E+00 0.0E+00 ];
TOT_LOSSRATE              (idx, [1:   2]) = [  1.00000E+00 0.0E+00 ];
TOT_CUTRATE               (idx, [1:   2]) = [  0.00000E+00 0.0E+00 ];
TOT_RR                    (idx, [1:   2]) = [  1.09069E+02 0.00018 ];
INI_FMASS                 (idx, 1)        =  0.00000E+00 ;
TOT_FMASS                 (idx, 1)        =  0.00000E+00 ;

% Six-factor formula:

SIX_FF_ETA                (idx, [1:   2]) = [  1.82531E+00 0.00016 ];
SIX_FF_F                  (idx, [1:   2]) = [  6.65608E-01 0.00016 ];
SIX_FF_P                  (idx, [1:   2]) = [  2.17175E-01 0.00024 ];
SIX_FF_EPSILON            (idx, [1:   2]) = [  3.11187E+00 0.00028 ];
SIX_FF_LF                 (idx, [1:   2]) = [  9.99963E-01 7.4E-07 ];
SIX_FF_LT                 (idx, [1:   2]) = [  9.82178E-01 1.4E-05 ];
SIX_FF_KINF               (idx, [1:   2]) = [  8.21016E-01 0.00015 ];
SIX_FF_KEFF               (idx, [1:   2]) = [  8.06354E-01 0.00015 ];

% Fission neutron and energy production:

NUBAR                     (idx, [1:   2]) = [  2.51057E+00 7.9E-06 ];
FISSE                     (idx, [1:   2]) = [  2.03324E+02 9.8E-07 ];

% Criticality eigenvalues:

ANA_KEFF                  (idx, [1:   6]) = [  8.06361E-01 0.00015  8.00760E-01 0.00015  5.59335E-03 0.00211 ];
IMP_KEFF                  (idx, [1:   2]) = [  8.06390E-01 8.5E-05 ];
COL_KEFF                  (idx, [1:   2]) = [  8.06428E-01 0.00012 ];
ABS_KEFF                  (idx, [1:   2]) = [  8.06390E-01 8.5E-05 ];
ABS_KINF                  (idx, [1:   2]) = [  8.21052E-01 8.3E-05 ];
GEOM_ALBEDO               (idx, [1:   6]) = [  1.00000E+00 0.0E+00  1.00000E+00 0.0E+00  1.00000E+00 0.0E+00 ];

% ALF (Average lethargy of neutrons causing fission):
% Based on E0 = 2.000000E+01 MeV

ANA_ALF                   (idx, [1:   2]) = [  1.06539E+01 0.00013 ];
IMP_ALF                   (idx, [1:   2]) = [  1.06552E+01 0.00011 ];

% EALF (Energy corresponding to average lethargy of neutrons causing fission):

ANA_EALF                  (idx, [1:   2]) = [  4.72626E-04 0.00142 ];
IMP_EALF                  (idx, [1:   2]) = [  4.71881E-04 0.00121 ];

% AFGE (Average energy of neutrons causing fission):

ANA_AFGE                  (idx, [1:   2]) = [  7.30708E-01 0.00036 ];
IMP_AFGE                  (idx, [1:   2]) = [  7.30553E-01 0.00022 ];

% Forward-weighted delayed neutron parameters:

FWD_ANA_BETA_ZERO         (idx, [1:  18]) = [  1.17846E-02 0.00114  2.72029E-04 0.00710  1.54476E-03 0.00297  8.39573E-04 0.00411  2.02367E-03 0.00272  3.65053E-03 0.00204  1.61538E-03 0.00284  1.21189E-03 0.00345  6.26767E-04 0.00448 ];
FWD_ANA_LAMBDA            (idx, [1:  18]) = [  5.69057E-01 0.00158  1.24667E-02 7.2E-09  2.82917E-02 0.0E+00  4.25244E-02 0.0E+00  1.33042E-01 0.0E+00  2.92467E-01 0.0E+00  6.66488E-01 2.7E-09  1.63478E+00 4.2E-09  3.55460E+00 3.0E-09 ];

% Beta-eff using Meulekamp's method:

ADJ_MEULEKAMP_BETA_EFF    (idx, [1:  18]) = [  7.16556E-03 0.00177  1.57627E-04 0.01159  9.53016E-04 0.00498  5.03394E-04 0.00685  1.22501E-03 0.00447  2.22433E-03 0.00319  9.77719E-04 0.00477  7.44218E-04 0.00560  3.80245E-04 0.00769 ];
ADJ_MEULEKAMP_LAMBDA      (idx, [1:  18]) = [  5.70082E-01 0.00263  1.24667E-02 7.2E-09  2.82917E-02 0.0E+00  4.25244E-02 0.0E+00  1.33042E-01 0.0E+00  2.92467E-01 0.0E+00  6.66488E-01 2.7E-09  1.63478E+00 4.2E-09  3.55460E+00 2.9E-09 ];

% Adjoint weighted time constants using Nauchi's method:

ADJ_NAUCHI_GEN_TIME       (idx, [1:   6]) = [  3.77377E-04 0.00073  3.77240E-04 0.00074  3.96394E-04 0.00786 ];
ADJ_NAUCHI_LIFETIME       (idx, [1:   6]) = [  3.04295E-04 0.00072  3.04185E-04 0.00072  3.19645E-04 0.00786 ];
ADJ_NAUCHI_BETA_EFF       (idx, [1:  18]) = [  6.94294E-03 0.00212  1.54986E-04 0.01376  9.23598E-04 0.00574  4.85386E-04 0.00802  1.19133E-03 0.00506  2.15075E-03 0.00381  9.47914E-04 0.00563  7.15905E-04 0.00662  3.73073E-04 0.00884 ];
ADJ_NAUCHI_LAMBDA         (idx, [1:  18]) = [  5.71100E-01 0.00309  1.24667E-02 7.1E-09  2.82917E-02 0.0E+00  4.25244E-02 0.0E+00  1.33042E-01 0.0E+00  2.92467E-01 0.0E+00  6.66488E-01 2.7E-09  1.63478E+00 4.2E-09  3.55460E+00 2.9E-09 ];

% Adjoint weighted time constants using IFP:

ADJ_IFP_GEN_TIME          (idx, [1:   6]) = [  4.24621E-04 0.00176  4.24325E-04 0.00177  4.66854E-04 0.01958 ];
ADJ_IFP_LIFETIME          (idx, [1:   6]) = [  3.42385E-04 0.00175  3.42147E-04 0.00176  3.76471E-04 0.01960 ];
ADJ_IFP_IMP_BETA_EFF      (idx, [1:  18]) = [  7.00216E-03 0.00756  1.58637E-04 0.05008  9.58060E-04 0.02055  4.70082E-04 0.02863  1.21467E-03 0.01856  2.20366E-03 0.01348  9.50740E-04 0.01991  7.14426E-04 0.02303  3.31888E-04 0.03320 ];
ADJ_IFP_IMP_LAMBDA        (idx, [1:  18]) = [  5.47820E-01 0.01090  1.24667E-02 0.0E+00  2.82917E-02 0.0E+00  4.25244E-02 0.0E+00  1.33042E-01 0.0E+00  2.92467E-01 0.0E+00  6.66488E-01 2.3E-09  1.63478E+00 4.5E-09  3.55460E+00 2.6E-09 ];
ADJ_IFP_ANA_BETA_EFF      (idx, [1:  18]) = [  6.99668E-03 0.00739  1.56935E-04 0.04907  9.54247E-04 0.02006  4.69135E-04 0.02793  1.21494E-03 0.01825  2.20132E-03 0.01311  9.52459E-04 0.01952  7.15454E-04 0.02231  3.32186E-04 0.03258 ];
ADJ_IFP_ANA_LAMBDA        (idx, [1:  18]) = [  5.48448E-01 0.01063  1.24667E-02 0.0E+00  2.82917E-02 0.0E+00  4.25244E-02 0.0E+00  1.33042E-01 0.0E+00  2.92467E-01 0.0E+00  6.66488E-01 2.3E-09  1.63478E+00 4.5E-09  3.55460E+00 2.8E-09 ];
ADJ_IFP_ROSSI_ALPHA       (idx, [1:   2]) = [ -1.65502E+01 0.00769 ];

% Adjoint weighted time constants using perturbation technique:

ADJ_PERT_GEN_TIME         (idx, [1:   2]) = [  4.04844E-04 0.00041 ];
ADJ_PERT_LIFETIME         (idx, [1:   2]) = [  3.26443E-04 0.00039 ];
ADJ_PERT_BETA_EFF         (idx, [1:   2]) = [  7.00569E-03 0.00128 ];
ADJ_PERT_ROSSI_ALPHA      (idx, [1:   2]) = [ -1.73067E+01 0.00130 ];

% Inverse neutron speed :

ANA_INV_SPD               (idx, [1:   2]) = [  1.60556E-06 0.00021 ];

% Analog slowing-down and thermal neutron lifetime (total/prompt/delayed):

ANA_SLOW_TIME             (idx, [1:   6]) = [  1.20639E-05 7.2E-05  1.20635E-05 7.3E-05  1.21154E-05 0.00082 ];
ANA_THERM_TIME            (idx, [1:   6]) = [  1.29914E-03 0.00028  1.29980E-03 0.00028  1.21354E-03 0.00332 ];
ANA_THERM_FRAC            (idx, [1:   6]) = [  2.31208E-01 0.00023  2.32175E-01 0.00023  1.50302E-01 0.00272 ];
ANA_DELAYED_EMTIME        (idx, [1:   2]) = [  1.07279E+01 0.00241 ];
ANA_MEAN_NCOL             (idx, [1:   4]) = [  1.06547E+02 0.00019  1.02508E+02 0.00031 ];

% Group constant generation:

GC_UNIVERSE_NAME          (idx, [1:  1])  = 'Be' ;

% Micro- and macro-group structures:

MICRO_NG                  (idx, 1)        = 70 ;
MICRO_E                   (idx, [1:  71]) = [  1.00000E-11  5.00000E-09  1.00000E-08  1.50000E-08  2.00000E-08  2.50000E-08  3.00000E-08  3.50000E-08  4.20000E-08  5.00000E-08  5.80000E-08  6.70000E-08  8.00000E-08  1.00000E-07  1.40000E-07  1.80000E-07  2.20000E-07  2.50000E-07  2.80000E-07  3.00000E-07  3.20000E-07  3.50000E-07  4.00000E-07  5.00000E-07  6.25000E-07  7.80000E-07  8.50000E-07  9.10000E-07  9.50000E-07  9.72000E-07  9.96000E-07  1.02000E-06  1.04500E-06  1.07100E-06  1.09700E-06  1.12300E-06  1.15000E-06  1.30000E-06  1.50000E-06  1.85500E-06  2.10000E-06  2.60000E-06  3.30000E-06  4.00000E-06  9.87700E-06  1.59680E-05  2.77000E-05  4.80520E-05  7.55014E-05  1.48728E-04  3.67262E-04  9.06898E-04  1.42510E-03  2.23945E-03  3.51910E-03  5.50000E-03  9.11800E-03  1.50300E-02  2.47800E-02  4.08500E-02  6.74300E-02  1.11000E-01  1.83000E-01  3.02500E-01  5.00000E-01  8.21000E-01  1.35300E+00  2.23100E+00  3.67900E+00  6.06550E+00  2.00000E+01 ];

MACRO_NG                  (idx, 1)        = 2 ;
MACRO_E                   (idx, [1:   3]) = [  1.00000E+37  6.25000E-07  0.00000E+00 ];

% Micro-group spectrum:

INF_MICRO_FLX             (idx, [1: 140]) = [  8.90276E+04 0.00194  3.97134E+05 0.00061  7.94015E+05 0.00061  1.77088E+06 0.00053  1.83039E+06 0.00046  2.08058E+06 0.00038  2.45928E+06 0.00036  2.49363E+06 0.00033  2.41330E+06 0.00032  2.38463E+06 0.00033  2.35444E+06 0.00033  2.28208E+06 0.00031  2.22229E+06 0.00032  2.16228E+06 0.00031  2.14910E+06 0.00032  1.86904E+06 0.00029  1.85998E+06 0.00032  1.82727E+06 0.00030  1.79130E+06 0.00032  3.46896E+06 0.00032  3.32588E+06 0.00033  2.40031E+06 0.00034  1.55561E+06 0.00033  1.85018E+06 0.00035  1.79945E+06 0.00033  1.53453E+06 0.00036  2.78055E+06 0.00036  5.82258E+05 0.00042  7.20013E+05 0.00040  6.44634E+05 0.00041  3.74613E+05 0.00049  6.44629E+05 0.00038  4.36264E+05 0.00040  3.74673E+05 0.00047  7.25628E+04 0.00077  7.16105E+04 0.00075  7.33681E+04 0.00081  7.53084E+04 0.00072  7.41994E+04 0.00081  7.30246E+04 0.00082  7.49018E+04 0.00078  7.02830E+04 0.00071  1.32347E+05 0.00055  2.10249E+05 0.00056  2.65830E+05 0.00046  6.91896E+05 0.00041  7.10597E+05 0.00045  7.30555E+05 0.00039  4.52076E+05 0.00045  3.13853E+05 0.00051  2.33636E+05 0.00054  2.61997E+05 0.00062  4.73069E+05 0.00053  6.52268E+05 0.00051  1.56685E+06 0.00043  3.98815E+06 0.00048  1.17179E+07 0.00048  1.28516E+07 0.00049  1.26148E+07 0.00046  1.12749E+07 0.00049  1.19590E+07 0.00054  1.39555E+07 0.00056  1.35818E+07 0.00049  1.05703E+07 0.00055  1.02318E+07 0.00059  1.09881E+07 0.00051  9.51622E+06 0.00054  8.47022E+06 0.00054  5.94948E+06 0.00056  2.18184E+06 0.00065 ];

% Integral parameters:

INF_KINF                  (idx, [1:   2]) = [  0.00000E+00 0.0E+00 ];

% Flux spectra in infinite geometry:

INF_FLX                   (idx, [1:   4]) = [  2.89941E+01 0.00031  7.55131E+01 0.00046 ];
INF_FISS_FLX              (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Reaction cross sections:

INF_TOT                   (idx, [1:   4]) = [  6.75717E-01 1.7E-05  7.40919E-01 8.8E-06 ];
INF_CAPT                  (idx, [1:   4]) = [  3.33409E-04 0.00029  9.18986E-04 3.1E-05 ];
INF_ABS                   (idx, [1:   4]) = [  3.33409E-04 0.00029  9.18986E-04 3.1E-05 ];
INF_FISS                  (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_NSF                   (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_NUBAR                 (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_KAPPA                 (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_INVV                  (idx, [1:   4]) = [  9.70940E-08 8.6E-05  3.85271E-06 3.1E-05 ];

% Total scattering cross sections:

INF_SCATT0                (idx, [1:   4]) = [  6.75383E-01 1.7E-05  7.40001E-01 8.8E-06 ];
INF_SCATT1                (idx, [1:   4]) = [  5.62488E-02 0.00016 -2.75215E-02 0.00034 ];
INF_SCATT2                (idx, [1:   4]) = [  4.36948E-03 0.00160 -4.92700E-02 0.00013 ];
INF_SCATT3                (idx, [1:   4]) = [  1.01770E-03 0.00609 -3.50297E-02 0.00019 ];
INF_SCATT4                (idx, [1:   4]) = [ -3.03560E-04 0.01535 -2.22133E-02 0.00028 ];
INF_SCATT5                (idx, [1:   4]) = [  1.76705E-04 0.02634 -1.17675E-02 0.00044 ];
INF_SCATT6                (idx, [1:   4]) = [ -5.74554E-04 0.00689 -4.75108E-03 0.00109 ];
INF_SCATT7                (idx, [1:   4]) = [  1.05524E-04 0.03487  2.55237E-03 0.00217 ];

% Total scattering production cross sections:

INF_SCATTP0               (idx, [1:   4]) = [  6.76230E-01 1.6E-05  7.40001E-01 8.8E-06 ];
INF_SCATTP1               (idx, [1:   4]) = [  5.64531E-02 0.00016 -2.75215E-02 0.00034 ];
INF_SCATTP2               (idx, [1:   4]) = [  4.41151E-03 0.00159 -4.92700E-02 0.00013 ];
INF_SCATTP3               (idx, [1:   4]) = [  1.02750E-03 0.00610 -3.50297E-02 0.00019 ];
INF_SCATTP4               (idx, [1:   4]) = [ -3.01073E-04 0.01554 -2.22133E-02 0.00028 ];
INF_SCATTP5               (idx, [1:   4]) = [  1.78666E-04 0.02617 -1.17675E-02 0.00044 ];
INF_SCATTP6               (idx, [1:   4]) = [ -5.73190E-04 0.00697 -4.75108E-03 0.00109 ];
INF_SCATTP7               (idx, [1:   4]) = [  1.06593E-04 0.03455  2.55237E-03 0.00217 ];

% Diffusion parameters:

INF_TRANSPXS              (idx, [1:   4]) = [  5.47663E-01 5.4E-05  6.87852E-01 4.3E-05 ];
INF_DIFFCOEF              (idx, [1:   4]) = [  6.08647E-01 5.4E-05  4.84600E-01 4.3E-05 ];

% Reduced absoption and removal:

INF_RABSXS                (idx, [1:   4]) = [ -5.13160E-04 0.00096  9.18986E-04 3.1E-05 ];
INF_REMXS                 (idx, [1:   4]) = [  8.14172E-03 0.00010  9.77799E-04 0.00043 ];

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

INF_S0                    (idx, [1:   8]) = [  6.67575E-01 1.7E-05  7.80855E-03 9.5E-05  5.89378E-05 0.00171  7.39942E-01 8.8E-06 ];
INF_S1                    (idx, [1:   8]) = [  5.83658E-02 0.00016 -2.11705E-03 0.00041  1.20509E-05 0.00392 -2.75336E-02 0.00034 ];
INF_S2                    (idx, [1:   8]) = [  4.63289E-03 0.00148 -2.63409E-04 0.00323 -5.05494E-06 0.00711 -4.92649E-02 0.00013 ];
INF_S3                    (idx, [1:   8]) = [  1.04135E-03 0.00592 -2.36572E-05 0.02989 -4.67599E-06 0.00719 -3.50250E-02 0.00019 ];
INF_S4                    (idx, [1:   8]) = [ -2.60719E-04 0.01848 -4.28408E-05 0.01579 -3.05477E-06 0.00948 -2.22103E-02 0.00028 ];
INF_S5                    (idx, [1:   8]) = [  1.40163E-04 0.03302  3.65422E-05 0.01809 -3.99338E-07 0.06252 -1.17671E-02 0.00044 ];
INF_S6                    (idx, [1:   8]) = [ -5.29467E-04 0.00753 -4.50861E-05 0.01170 -1.26364E-06 0.01908 -4.74981E-03 0.00109 ];
INF_S7                    (idx, [1:   8]) = [  5.22207E-05 0.06947  5.33035E-05 0.00771  1.58251E-07 0.10578  2.55221E-03 0.00217 ];

% Scattering production matrixes:

INF_SP0                   (idx, [1:   8]) = [  6.68421E-01 1.6E-05  7.80855E-03 9.5E-05  5.89378E-05 0.00171  7.39942E-01 8.8E-06 ];
INF_SP1                   (idx, [1:   8]) = [  5.85701E-02 0.00016 -2.11705E-03 0.00041  1.20509E-05 0.00392 -2.75336E-02 0.00034 ];
INF_SP2                   (idx, [1:   8]) = [  4.67492E-03 0.00147 -2.63409E-04 0.00323 -5.05494E-06 0.00711 -4.92649E-02 0.00013 ];
INF_SP3                   (idx, [1:   8]) = [  1.05115E-03 0.00593 -2.36572E-05 0.02989 -4.67599E-06 0.00719 -3.50250E-02 0.00019 ];
INF_SP4                   (idx, [1:   8]) = [ -2.58232E-04 0.01872 -4.28409E-05 0.01579 -3.05477E-06 0.00948 -2.22103E-02 0.00028 ];
INF_SP5                   (idx, [1:   8]) = [  1.42124E-04 0.03271  3.65423E-05 0.01809 -3.99338E-07 0.06252 -1.17671E-02 0.00044 ];
INF_SP6                   (idx, [1:   8]) = [ -5.28104E-04 0.00762 -4.50862E-05 0.01170 -1.26364E-06 0.01908 -4.74981E-03 0.00109 ];
INF_SP7                   (idx, [1:   8]) = [  5.32890E-05 0.06820  5.33035E-05 0.00771  1.58251E-07 0.10578  2.55221E-03 0.00217 ];

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

CMM_TRANSPXS              (idx, [1:   4]) = [  4.40113E-01 0.00013 -1.97279E+00 0.00055 ];
CMM_TRANSPXS_X            (idx, [1:   4]) = [  4.49403E-01 0.00018 -1.80311E+00 0.00049 ];
CMM_TRANSPXS_Y            (idx, [1:   4]) = [  4.35664E-01 0.00017 -2.07014E+00 0.00067 ];
CMM_TRANSPXS_Z            (idx, [1:   4]) = [  4.35558E-01 0.00018 -2.07026E+00 0.00061 ];
CMM_DIFFCOEF              (idx, [1:   4]) = [  7.57382E-01 0.00013 -1.68968E-01 0.00055 ];
CMM_DIFFCOEF_X            (idx, [1:   4]) = [  7.41726E-01 0.00018 -1.84867E-01 0.00049 ];
CMM_DIFFCOEF_Y            (idx, [1:   4]) = [  7.65117E-01 0.00017 -1.61023E-01 0.00067 ];
CMM_DIFFCOEF_Z            (idx, [1:   4]) = [  7.65303E-01 0.00018 -1.61014E-01 0.00061 ];

% Delayed neutron parameters (Meulekamp method):

BETA_EFF                  (idx, [1:  18]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
LAMBDA                    (idx, [1:  18]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
