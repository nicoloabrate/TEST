
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
INPUT_FILE_NAME           (idx, [1: 11])  = 'reactor_hfp' ;
WORKING_DIRECTORY         (idx, [1: 35])  = '/home/abrate/SERPENT-2/c5g7/c5g7_2g' ;
HOSTNAME                  (idx, [1:  7])  = 'vpcen13' ;
CPU_TYPE                  (idx, [1: 41])  = 'Intel(R) Xeon(R) CPU E5-2630 v3 @ 2.40GHz' ;
CPU_MHZ                   (idx, 1)        = 4294967295.0 ;
START_DATE                (idx, [1: 24])  = 'Thu Oct 18 14:24:16 2018' ;
COMPLETE_DATE             (idx, [1: 24])  = 'Thu Oct 18 21:32:41 2018' ;

% Run parameters:

POP                       (idx, 1)        = 1000000 ;
CYCLES                    (idx, 1)        = 1000 ;
SKIP                      (idx, 1)        = 550 ;
BATCH_INTERVAL            (idx, 1)        = 1 ;
SRC_NORM_MODE             (idx, 1)        = 2 ;
SEED                      (idx, 1)        = 1539865456 ;
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
OMP_HISTORY_PROFILE       (idx, [1:  30]) = [  1.05998E+00  9.95045E-01  9.85876E-01  1.01079E+00  9.86691E-01  1.00524E+00  9.91255E-01  9.91210E-01  1.00613E+00  1.00907E+00  1.00228E+00  1.00326E+00  1.00056E+00  1.02414E+00  9.87812E-01  9.90574E-01  1.01188E+00  9.89691E-01  1.00034E+00  9.95425E-01  9.96001E-01  9.91779E-01  9.92047E-01  9.96249E-01  9.90974E-01  1.00485E+00  9.91236E-01  9.86254E-01  9.95371E-01  1.00800E+00  ];
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
ST_FRAC                   (idx, [1:   4]) = [  3.97956E-02 5.3E-05  0.00000E+00 0.0E+00 ];
DT_FRAC                   (idx, [1:   4]) = [  9.60204E-01 2.2E-06  0.00000E+00 0.0E+00 ];
DT_EFF                    (idx, [1:   4]) = [  6.74317E-01 6.7E-06  0.00000E+00 0.0E+00 ];
REA_SAMPLING_EFF          (idx, [1:   4]) = [  1.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
REA_SAMPLING_FAIL         (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
TOT_COL_EFF               (idx, [1:   4]) = [  6.75845E-01 6.7E-06  0.00000E+00 0.0E+00 ];
AVG_TRACKING_LOOPS        (idx, [1:   8]) = [  3.12822E+00 2.2E-05  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
AVG_TRACKS                (idx, [1:   4]) = [  2.84241E+01 2.4E-05  0.00000E+00 0.0E+00 ];
AVG_REAL_COL              (idx, [1:   4]) = [  2.83913E+01 2.4E-05  0.00000E+00 0.0E+00 ];
AVG_VIRT_COL              (idx, [1:   4]) = [  1.36173E+01 3.6E-05  0.00000E+00 0.0E+00 ];
AVG_SURF_CROSS            (idx, [1:   4]) = [  1.56776E+00 5.6E-05  0.00000E+00 0.0E+00 ];
LOST_PARTICLES            (idx, 1)        = 0 ;

% Run statistics:

CYCLE_IDX                 (idx, 1)        = 1000 ;
SOURCE_POPULATION         (idx, 1)        = 1000000700 ;
MEAN_POP_SIZE             (idx, [1:  2])  = [  1.00000E+06 0.00005 ];
MEAN_POP_WGT              (idx, [1:  2])  = [  1.00000E+06 0.00005 ];
SIMULATION_COMPLETED      (idx, 1)        = 1 ;

% Running times:

TOT_CPU_TIME              (idx, 1)        =  9.77436E+03 ;
RUNNING_TIME              (idx, 1)        =  4.28418E+02 ;
INIT_TIME                 (idx, [1:  2])  = [  3.70833E-01  3.70833E-01 ];
PROCESS_TIME              (idx, [1:  2])  = [  8.25000E-03  8.23333E-03 ];
TRANSPORT_CYCLE_TIME      (idx, [1:  3])  = [  4.28039E+02  4.28039E+02  0.00000E+00 ];
MPI_OVERHEAD_TIME         (idx, [1:  2])  = [  0.00000E+00  0.00000E+00 ];
ESTIMATED_RUNNING_TIME    (idx, [1:  2])  = [  4.28413E+02  0.00000E+00 ];
CPU_USAGE                 (idx, 1)        = 22.81501 ;
TRANSPORT_CPU_USAGE       (idx, [1:   2]) = [  2.30476E+01 0.00138 ];
OMP_PARALLEL_FRAC         (idx, 1)        =  7.52072E-01 ;

% Memory usage:

AVAIL_MEM                 (idx, 1)        = 96872.24 ;
ALLOC_MEMSIZE             (idx, 1)        = 10422.49;
MEMSIZE                   (idx, 1)        = 10198.59;
XS_MEMSIZE                (idx, 1)        = 3343.50;
MAT_MEMSIZE               (idx, 1)        = 325.59;
RES_MEMSIZE               (idx, 1)        = 5.37;
MISC_MEMSIZE              (idx, 1)        = 6524.13;
UNKNOWN_MEMSIZE           (idx, 1)        = 0.00;
UNUSED_MEMSIZE            (idx, 1)        = 223.90;

% Geometry parameters:

TOT_CELLS                 (idx, 1)        = 3 ;
UNION_CELLS               (idx, 1)        = 0 ;

% Neutron energy grid:

NEUTRON_ERG_TOL           (idx, 1)        =  0.00000E+00 ;
NEUTRON_ERG_NE            (idx, 1)        = 887926 ;
NEUTRON_EMIN              (idx, 1)        =  1.00000E-11 ;
NEUTRON_EMAX              (idx, 1)        =  2.00000E+01 ;

% Unresolved resonance probability table sampling:

URES_DILU_CUT             (idx, 1)        =  1.00000E-09 ;
URES_EMIN                 (idx, 1)        =  1.00000E+37 ;
URES_EMAX                 (idx, 1)        = -1.00000E+37 ;
URES_AVAIL                (idx, 1)        = 27 ;
URES_USED                 (idx, 1)        = 0 ;

% Nuclides and reaction channels:

TOT_NUCLIDES              (idx, 1)        = 52 ;
TOT_TRANSPORT_NUCLIDES    (idx, 1)        = 52 ;
TOT_DOSIMETRY_NUCLIDES    (idx, 1)        = 0 ;
TOT_DECAY_NUCLIDES        (idx, 1)        = 0 ;
TOT_PHOTON_NUCLIDES       (idx, 1)        = 0 ;
TOT_REA_CHANNELS          (idx, 1)        = 1853 ;
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

NORM_COEF                 (idx, [1:   4]) = [  2.98886E+09 3.2E-05  0.00000E+00 0.0E+00 ];

% Analog reaction rate estimators:

CONVERSION_RATIO          (idx, [1:   2]) = [  4.73627E-01 7.4E-05 ];
U235_FISS                 (idx, [1:   4]) = [  7.47552E+14 6.0E-05  5.57643E-01 4.9E-05 ];
U238_FISS                 (idx, [1:   4]) = [  7.81159E+13 0.00019  5.82712E-02 0.00018 ];
PU239_FISS                (idx, [1:   4]) = [  4.22658E+14 8.5E-05  3.15286E-01 7.7E-05 ];
PU240_FISS                (idx, [1:   4]) = [  3.58108E+12 0.00091  2.67133E-03 0.00091 ];
PU241_FISS                (idx, [1:   4]) = [  8.71458E+13 0.00019  6.50071E-02 0.00018 ];
U235_CAPT                 (idx, [1:   4]) = [  1.67626E+14 0.00013  1.07819E-01 0.00013 ];
U238_CAPT                 (idx, [1:   4]) = [  6.02469E+14 7.7E-05  3.87513E-01 5.4E-05 ];
PU239_CAPT                (idx, [1:   4]) = [  2.21212E+14 0.00011  1.42286E-01 0.00010 ];
PU240_CAPT                (idx, [1:   4]) = [  1.87703E+14 0.00013  1.20732E-01 0.00012 ];
PU241_CAPT                (idx, [1:   4]) = [  2.96373E+13 0.00031  1.90630E-02 0.00031 ];

% Neutron balance (particles/weight):

BALA_SRC_NEUTRON_SRC     (idx, [1:  2])  = [ 0 0.00000E+00 ];
BALA_SRC_NEUTRON_FISS    (idx, [1:  2])  = [ 1000000700 1.00000E+09 ];
BALA_SRC_NEUTRON_NXN     (idx, [1:  2])  = [ 0 1.46559E+06 ];
BALA_SRC_NEUTRON_VR      (idx, [1:  2])  = [ 0 0.00000E+00 ];
BALA_SRC_NEUTRON_TOT     (idx, [1:  2])  = [ 1000000700 1.00147E+09 ];

BALA_LOSS_NEUTRON_CAPT    (idx, [1:  2])  = [ 519374634 5.20166E+08 ];
BALA_LOSS_NEUTRON_FISS    (idx, [1:  2])  = [ 447869141 4.48518E+08 ];
BALA_LOSS_NEUTRON_LEAK    (idx, [1:  2])  = [ 32756925 3.27816E+07 ];
BALA_LOSS_NEUTRON_CUT     (idx, [1:  2])  = [ 0 0.00000E+00 ];
BALA_LOSS_NEUTRON_ERR     (idx, [1:  2])  = [ 0 0.00000E+00 ];
BALA_LOSS_NEUTRON_TOT     (idx, [1:  2])  = [ 1000000700 1.00147E+09 ];

BALA_NEUTRON_DIFF         (idx, [1:  2])  = [ 0 -8.16584E-05 ];

% Normalized total reaction rates (neutrons):

TOT_POWER                 (idx, [1:   2]) = [  4.40000E+04 0.0E+00 ];
TOT_POWDENS               (idx, [1:   2]) = [  3.87002E-02 4.1E-09 ];
TOT_GENRATE               (idx, [1:   2]) = [  3.51593E+15 3.4E-06 ];
TOT_FISSRATE              (idx, [1:   2]) = [  1.34045E+15 7.0E-07 ];
TOT_CAPTRATE              (idx, [1:   2]) = [  1.55462E+15 3.3E-05 ];
TOT_ABSRATE               (idx, [1:   2]) = [  2.89507E+15 1.8E-05 ];
TOT_SRCRATE               (idx, [1:   2]) = [  2.98886E+15 3.2E-05 ];
TOT_FLUX                  (idx, [1:   2]) = [  1.34975E+17 4.2E-05 ];
TOT_PHOTON_PRODRATE       (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
TOT_LEAKRATE              (idx, [1:   2]) = [  9.79798E+13 0.00020 ];
ALBEDO_LEAKRATE           (idx, [1:   2]) = [  0.00000E+00 0.0E+00 ];
TOT_LOSSRATE              (idx, [1:   2]) = [  2.99305E+15 1.9E-05 ];
TOT_CUTRATE               (idx, [1:   2]) = [  0.00000E+00 0.0E+00 ];
TOT_RR                    (idx, [1:   2]) = [  8.49754E+16 3.5E-05 ];
INI_FMASS                 (idx, 1)        =  1.13695E+00 ;
TOT_FMASS                 (idx, 1)        =  1.13695E+00 ;

% Six-factor formula:

SIX_FF_ETA                (idx, [1:   2]) = [  1.79192E+00 3.0E-05 ];
SIX_FF_F                  (idx, [1:   2]) = [  8.81927E-01 1.6E-05 ];
SIX_FF_P                  (idx, [1:   2]) = [  5.64442E-01 2.7E-05 ];
SIX_FF_EPSILON            (idx, [1:   2]) = [  1.36355E+00 2.8E-05 ];
SIX_FF_LF                 (idx, [1:   2]) = [  9.67231E-01 6.3E-06 ];
SIX_FF_LT                 (idx, [1:   2]) = [  9.99986E-01 1.2E-07 ];
SIX_FF_KINF               (idx, [1:   2]) = [  1.21630E+00 3.5E-05 ];
SIX_FF_KEFF               (idx, [1:   2]) = [  1.17643E+00 3.5E-05 ];

% Fission neutron and energy production:

NUBAR                     (idx, [1:   2]) = [  2.62295E+00 4.1E-06 ];
FISSE                     (idx, [1:   2]) = [  2.04876E+02 7.0E-07 ];

% Criticality eigenvalues:

ANA_KEFF                  (idx, [1:   6]) = [  1.17642E+00 3.6E-05  1.16979E+00 3.5E-05  6.63645E-03 0.00061 ];
IMP_KEFF                  (idx, [1:   2]) = [  1.17644E+00 1.9E-05 ];
COL_KEFF                  (idx, [1:   2]) = [  1.17634E+00 3.1E-05 ];
ABS_KEFF                  (idx, [1:   2]) = [  1.17644E+00 1.9E-05 ];
ABS_KINF                  (idx, [1:   2]) = [  1.21631E+00 1.7E-05 ];
GEOM_ALBEDO               (idx, [1:   6]) = [  1.00000E+00 0.0E+00  1.00000E+00 0.0E+00  1.00000E+00 0.0E+00 ];

% ALF (Average lethargy of neutrons causing fission):
% Based on E0 = 2.000000E+01 MeV

ANA_ALF                   (idx, [1:   2]) = [  1.70670E+01 1.4E-05 ];
IMP_ALF                   (idx, [1:   2]) = [  1.70671E+01 7.0E-06 ];

% EALF (Energy corresponding to average lethargy of neutrons causing fission):

ANA_EALF                  (idx, [1:   2]) = [  7.74329E-07 0.00023 ];
IMP_EALF                  (idx, [1:   2]) = [  7.74243E-07 0.00012 ];

% AFGE (Average energy of neutrons causing fission):

ANA_AFGE                  (idx, [1:   2]) = [  2.19933E-01 0.00019 ];
IMP_AFGE                  (idx, [1:   2]) = [  2.19959E-01 8.0E-05 ];

% Forward-weighted delayed neutron parameters:

FWD_ANA_BETA_ZERO         (idx, [1:  18]) = [  4.82130E-03 0.00042  1.30725E-04 0.00257  7.56929E-04 0.00104  3.84487E-04 0.00152  8.71609E-04 0.00098  1.52212E-03 0.00076  5.37696E-04 0.00126  4.47156E-04 0.00138  1.70578E-04 0.00228 ];
FWD_ANA_LAMBDA            (idx, [1:  18]) = [  4.76271E-01 0.00068  1.24667E-02 7.2E-09  2.82917E-02 0.0E+00  4.25244E-02 0.0E+00  1.33042E-01 0.0E+00  2.92467E-01 0.0E+00  6.66488E-01 2.7E-09  1.63478E+00 4.2E-09  3.55460E+00 2.9E-09 ];

% Beta-eff using Meulekamp's method:

ADJ_MEULEKAMP_BETA_EFF    (idx, [1:  18]) = [  5.75544E-03 0.00059  1.58974E-04 0.00363  8.90692E-04 0.00143  4.65287E-04 0.00214  1.04627E-03 0.00140  1.82129E-03 0.00109  6.39182E-04 0.00181  5.31401E-04 0.00196  2.02338E-04 0.00320 ];
ADJ_MEULEKAMP_LAMBDA      (idx, [1:  18]) = [  4.74834E-01 0.00096  1.24667E-02 7.2E-09  2.82917E-02 0.0E+00  4.25244E-02 0.0E+00  1.33042E-01 0.0E+00  2.92467E-01 0.0E+00  6.66488E-01 2.7E-09  1.63478E+00 4.2E-09  3.55460E+00 2.8E-09 ];

% Adjoint weighted time constants using Nauchi's method:

ADJ_NAUCHI_GEN_TIME       (idx, [1:   6]) = [  1.37554E-05 8.8E-05  1.37431E-05 8.8E-05  1.58697E-05 0.00092 ];
ADJ_NAUCHI_LIFETIME       (idx, [1:   6]) = [  1.61821E-05 7.8E-05  1.61677E-05 7.9E-05  1.86694E-05 0.00092 ];
ADJ_NAUCHI_BETA_EFF       (idx, [1:  18]) = [  5.64128E-03 0.00061  1.54190E-04 0.00385  8.77323E-04 0.00153  4.54172E-04 0.00229  1.02330E-03 0.00144  1.78177E-03 0.00110  6.28404E-04 0.00194  5.22085E-04 0.00208  2.00030E-04 0.00344 ];
ADJ_NAUCHI_LAMBDA         (idx, [1:  18]) = [  4.76249E-01 0.00102  1.24667E-02 7.2E-09  2.82917E-02 0.0E+00  4.25244E-02 0.0E+00  1.33042E-01 0.0E+00  2.92467E-01 0.0E+00  6.66488E-01 2.7E-09  1.63478E+00 4.2E-09  3.55460E+00 2.9E-09 ];

% Adjoint weighted time constants using IFP:

ADJ_IFP_GEN_TIME          (idx, [1:   6]) = [  1.42662E-05 0.00020  1.42537E-05 0.00020  1.63573E-05 0.00224 ];
ADJ_IFP_LIFETIME          (idx, [1:   6]) = [  1.67831E-05 0.00019  1.67683E-05 0.00019  1.92430E-05 0.00224 ];
ADJ_IFP_IMP_BETA_EFF      (idx, [1:  18]) = [  5.93767E-03 0.00196  1.64262E-04 0.01215  9.10746E-04 0.00496  4.89669E-04 0.00702  1.08477E-03 0.00439  1.88124E-03 0.00338  6.56079E-04 0.00590  5.45767E-04 0.00655  2.05135E-04 0.01031 ];
ADJ_IFP_IMP_LAMBDA        (idx, [1:  18]) = [  4.71972E-01 0.00300  1.24667E-02 7.2E-09  2.82917E-02 0.0E+00  4.25244E-02 0.0E+00  1.33042E-01 0.0E+00  2.92467E-01 0.0E+00  6.66488E-01 2.7E-09  1.63478E+00 4.2E-09  3.55460E+00 2.9E-09 ];
ADJ_IFP_ANA_BETA_EFF      (idx, [1:  18]) = [  5.93395E-03 0.00191  1.63590E-04 0.01180  9.11693E-04 0.00490  4.88490E-04 0.00685  1.08573E-03 0.00431  1.87904E-03 0.00327  6.55468E-04 0.00580  5.45385E-04 0.00633  2.04551E-04 0.01004 ];
ADJ_IFP_ANA_LAMBDA        (idx, [1:  18]) = [  4.71630E-01 0.00292  1.24667E-02 7.2E-09  2.82917E-02 0.0E+00  4.25244E-02 0.0E+00  1.33042E-01 0.0E+00  2.92467E-01 0.0E+00  6.66488E-01 2.7E-09  1.63478E+00 4.2E-09  3.55460E+00 2.9E-09 ];
ADJ_IFP_ROSSI_ALPHA       (idx, [1:   2]) = [ -4.16587E+02 0.00197 ];

% Adjoint weighted time constants using perturbation technique:

ADJ_PERT_GEN_TIME         (idx, [1:   2]) = [  1.40127E-05 5.4E-05 ];
ADJ_PERT_LIFETIME         (idx, [1:   2]) = [  1.64848E-05 4.0E-05 ];
ADJ_PERT_BETA_EFF         (idx, [1:   2]) = [  5.84573E-03 0.00040 ];
ADJ_PERT_ROSSI_ALPHA      (idx, [1:   2]) = [ -4.17175E+02 0.00041 ];

% Inverse neutron speed :

ANA_INV_SPD               (idx, [1:   2]) = [  2.67832E-07 5.4E-05 ];

% Analog slowing-down and thermal neutron lifetime (total/prompt/delayed):

ANA_SLOW_TIME             (idx, [1:   6]) = [  2.71638E-06 4.1E-05  2.71631E-06 4.1E-05  2.72612E-06 0.00050 ];
ANA_THERM_TIME            (idx, [1:   6]) = [  1.83718E-05 5.0E-05  1.83643E-05 5.0E-05  1.95828E-05 0.00061 ];
ANA_THERM_FRAC            (idx, [1:   6]) = [  5.46201E-01 2.8E-05  5.45460E-01 2.8E-05  6.99248E-01 0.00068 ];
ANA_DELAYED_EMTIME        (idx, [1:   2]) = [  1.22823E+01 0.00095 ];
ANA_MEAN_NCOL             (idx, [1:   4]) = [  2.83913E+01 2.4E-05  2.90621E+01 2.6E-05 ];

% Group constant generation:

GC_UNIVERSE_NAME          (idx, [1:  4])  = 'C5G7' ;

% Micro- and macro-group structures:

MICRO_NG                  (idx, 1)        = 70 ;
MICRO_E                   (idx, [1:  71]) = [  1.00000E-11  5.00000E-09  1.00000E-08  1.50000E-08  2.00000E-08  2.50000E-08  3.00000E-08  3.50000E-08  4.20000E-08  5.00000E-08  5.80000E-08  6.70000E-08  8.00000E-08  1.00000E-07  1.40000E-07  1.80000E-07  2.20000E-07  2.50000E-07  2.80000E-07  3.00000E-07  3.20000E-07  3.50000E-07  4.00000E-07  5.00000E-07  6.25000E-07  7.80000E-07  8.50000E-07  9.10000E-07  9.50000E-07  9.72000E-07  9.96000E-07  1.02000E-06  1.04500E-06  1.07100E-06  1.09700E-06  1.12300E-06  1.15000E-06  1.30000E-06  1.50000E-06  1.85500E-06  2.10000E-06  2.60000E-06  3.30000E-06  4.00000E-06  9.87700E-06  1.59680E-05  2.77000E-05  4.80520E-05  7.55014E-05  1.48728E-04  3.67262E-04  9.06898E-04  1.42510E-03  2.23945E-03  3.51910E-03  5.50000E-03  9.11800E-03  1.50300E-02  2.47800E-02  4.08500E-02  6.74300E-02  1.11000E-01  1.83000E-01  3.02500E-01  5.00000E-01  8.21000E-01  1.35300E+00  2.23100E+00  3.67900E+00  6.06550E+00  2.00000E+01 ];

MACRO_NG                  (idx, 1)        = 2 ;
MACRO_E                   (idx, [1:   3]) = [  1.00000E+37  6.25000E-07  0.00000E+00 ];

% Micro-group spectrum:

INF_MICRO_FLX             (idx, [1: 140]) = [  6.84061E+06 0.00023  2.69217E+07 0.00011  5.60136E+07 6.5E-05  6.38667E+07 5.9E-05  6.50627E+07 7.4E-05  7.76104E+07 8.4E-05  5.78963E+07 9.3E-05  5.02982E+07 9.0E-05  4.02456E+07 0.00011  3.04419E+07 9.5E-05  2.57043E+07 9.6E-05  2.17951E+07 7.6E-05  2.16680E+07 9.7E-05  1.85316E+07 9.1E-05  1.67936E+07 6.1E-05  1.51960E+07 9.6E-05  1.55088E+07 9.0E-05  1.58988E+07 9.4E-05  1.52994E+07 9.5E-05  2.94150E+07 7.5E-05  2.74326E+07 6.1E-05  2.01105E+07 7.3E-05  1.27923E+07 9.9E-05  1.47807E+07 9.2E-05  1.37882E+07 8.3E-05  1.18459E+07 0.00010  1.96117E+07 9.0E-05  4.34537E+06 0.00012  5.20701E+06 0.00010  4.71887E+06 0.00011  2.72796E+06 0.00014  4.67887E+06 0.00011  3.07132E+06 0.00015  2.38972E+06 0.00015  4.11759E+05 0.00030  3.88913E+05 0.00031  3.88357E+05 0.00031  3.94743E+05 0.00038  3.91230E+05 0.00026  3.92472E+05 0.00034  4.11682E+05 0.00026  3.95494E+05 0.00026  7.63817E+05 0.00022  1.26024E+06 0.00021  1.64216E+06 0.00014  4.33726E+06 0.00012  4.33525E+06 0.00012  4.13995E+06 0.00013  2.25395E+06 0.00017  1.42435E+06 0.00017  9.94903E+05 0.00016  1.05446E+06 0.00021  1.74081E+06 0.00015  2.04509E+06 0.00016  3.52796E+06 0.00011  5.00538E+06 9.4E-05  8.00436E+06 9.6E-05  6.03962E+06 0.00010  4.83739E+06 0.00010  3.91806E+06 0.00012  3.80212E+06 0.00012  4.00248E+06 0.00011  3.61573E+06 0.00011  2.58700E+06 0.00015  2.52820E+06 0.00012  2.31517E+06 0.00013  2.06462E+06 0.00015  1.70932E+06 0.00012  1.13710E+06 0.00015  4.15048E+05 0.00029 ];

% Integral parameters:

INF_KINF                  (idx, [1:   2]) = [  1.21621E+00 3.0E-05 ];

% Flux spectra in infinite geometry:

INF_FLX                   (idx, [1:   4]) = [  1.23991E+17 5.4E-05  1.09838E+16 4.1E-05 ];
INF_FISS_FLX              (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Reaction cross sections:

INF_TOT                   (idx, [1:   4]) = [  5.48269E-01 1.3E-05  1.54728E+00 1.7E-05 ];
INF_CAPT                  (idx, [1:   4]) = [  7.39945E-03 4.2E-05  5.80082E-02 5.7E-05 ];
INF_ABS                   (idx, [1:   4]) = [  1.01894E-02 4.0E-05  1.48553E-01 4.4E-05 ];
INF_FISS                  (idx, [1:   4]) = [  2.78992E-03 5.4E-05  9.05446E-02 4.0E-05 ];
INF_NSF                   (idx, [1:   4]) = [  7.56024E-03 5.4E-05  2.34757E-01 4.5E-05 ];
INF_NUBAR                 (idx, [1:   4]) = [  2.70984E+00 7.4E-06  2.59272E+00 8.4E-06 ];
INF_KAPPA                 (idx, [1:   4]) = [  2.05880E+02 1.2E-06  2.04527E+02 1.4E-06 ];
INF_INVV                  (idx, [1:   4]) = [  4.97658E-08 4.7E-05  2.72948E-06 1.9E-05 ];

% Total scattering cross sections:

INF_SCATT0                (idx, [1:   4]) = [  5.38078E-01 1.3E-05  1.39872E+00 1.9E-05 ];
INF_SCATT1                (idx, [1:   4]) = [  2.09891E-01 3.8E-05  3.31481E-01 4.4E-05 ];
INF_SCATT2                (idx, [1:   4]) = [  8.56698E-02 3.5E-05  7.31107E-02 0.00010 ];
INF_SCATT3                (idx, [1:   4]) = [  6.82012E-03 0.00021  1.99137E-02 0.00032 ];
INF_SCATT4                (idx, [1:   4]) = [ -8.00313E-03 0.00019 -8.02894E-03 0.00071 ];
INF_SCATT5                (idx, [1:   4]) = [  3.56332E-04 0.00356  3.75573E-03 0.00126 ];
INF_SCATT6                (idx, [1:   4]) = [  4.23861E-03 0.00028 -1.31923E-02 0.00039 ];
INF_SCATT7                (idx, [1:   4]) = [  4.02820E-04 0.00185  1.77922E-03 0.00264 ];

% Total scattering production cross sections:

INF_SCATTP0               (idx, [1:   4]) = [  5.38114E-01 1.3E-05  1.39872E+00 1.9E-05 ];
INF_SCATTP1               (idx, [1:   4]) = [  2.09891E-01 3.8E-05  3.31481E-01 4.4E-05 ];
INF_SCATTP2               (idx, [1:   4]) = [  8.56699E-02 3.5E-05  7.31107E-02 0.00010 ];
INF_SCATTP3               (idx, [1:   4]) = [  6.82013E-03 0.00021  1.99137E-02 0.00032 ];
INF_SCATTP4               (idx, [1:   4]) = [ -8.00314E-03 0.00019 -8.02894E-03 0.00071 ];
INF_SCATTP5               (idx, [1:   4]) = [  3.56318E-04 0.00357  3.75573E-03 0.00126 ];
INF_SCATTP6               (idx, [1:   4]) = [  4.23860E-03 0.00028 -1.31923E-02 0.00039 ];
INF_SCATTP7               (idx, [1:   4]) = [  4.02823E-04 0.00185  1.77922E-03 0.00264 ];

% Diffusion parameters:

INF_TRANSPXS              (idx, [1:   4]) = [  2.46126E-01 2.3E-05  1.04737E+00 2.2E-05 ];
INF_DIFFCOEF              (idx, [1:   4]) = [  1.35432E+00 2.3E-05  3.18257E-01 2.2E-05 ];

% Reduced absoption and removal:

INF_RABSXS                (idx, [1:   4]) = [  1.01540E-02 4.1E-05  1.48553E-01 4.4E-05 ];
INF_REMXS                 (idx, [1:   4]) = [  2.34235E-02 4.3E-05  1.49378E-01 4.6E-05 ];

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

INF_CHIT                  (idx, [1:   4]) = [  1.00000E+00 0.0E+00  9.99046E-10 1.00000 ];
INF_CHIP                  (idx, [1:   4]) = [  1.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_CHID                  (idx, [1:   4]) = [  1.00000E+00 1.8E-07  1.76900E-07 1.00000 ];

% Scattering matrixes:

INF_S0                    (idx, [1:   8]) = [  5.24845E-01 1.3E-05  1.32331E-02 5.1E-05  8.20127E-04 0.00058  1.39790E+00 1.9E-05 ];
INF_S1                    (idx, [1:   8]) = [  2.06139E-01 3.8E-05  3.75181E-03 8.2E-05  3.16698E-04 0.00105  3.31164E-01 4.4E-05 ];
INF_S2                    (idx, [1:   8]) = [  8.69367E-02 3.5E-05 -1.26689E-03 0.00016  1.86957E-04 0.00143  7.29238E-02 0.00010 ];
INF_S3                    (idx, [1:   8]) = [  8.21137E-03 0.00017 -1.39125E-03 0.00012  9.04429E-05 0.00216  1.98233E-02 0.00033 ];
INF_S4                    (idx, [1:   8]) = [ -7.58737E-03 0.00020 -4.15761E-04 0.00037  2.15557E-05 0.00737 -8.05049E-03 0.00070 ];
INF_S5                    (idx, [1:   8]) = [  3.35079E-04 0.00377  2.12526E-05 0.00644 -1.58589E-05 0.01019  3.77159E-03 0.00125 ];
INF_S6                    (idx, [1:   8]) = [  4.39362E-03 0.00027 -1.55008E-04 0.00081 -2.98868E-05 0.00469 -1.31624E-02 0.00039 ];
INF_S7                    (idx, [1:   8]) = [  5.50483E-04 0.00128 -1.47662E-04 0.00092 -3.17351E-05 0.00398  1.81095E-03 0.00260 ];

% Scattering production matrixes:

INF_SP0                   (idx, [1:   8]) = [  5.24880E-01 1.3E-05  1.32331E-02 5.1E-05  8.20127E-04 0.00058  1.39790E+00 1.9E-05 ];
INF_SP1                   (idx, [1:   8]) = [  2.06139E-01 3.8E-05  3.75181E-03 8.2E-05  3.16698E-04 0.00105  3.31164E-01 4.4E-05 ];
INF_SP2                   (idx, [1:   8]) = [  8.69368E-02 3.5E-05 -1.26689E-03 0.00016  1.86957E-04 0.00143  7.29238E-02 0.00010 ];
INF_SP3                   (idx, [1:   8]) = [  8.21138E-03 0.00017 -1.39125E-03 0.00012  9.04429E-05 0.00216  1.98233E-02 0.00033 ];
INF_SP4                   (idx, [1:   8]) = [ -7.58738E-03 0.00020 -4.15761E-04 0.00037  2.15557E-05 0.00737 -8.05049E-03 0.00070 ];
INF_SP5                   (idx, [1:   8]) = [  3.35065E-04 0.00378  2.12526E-05 0.00644 -1.58589E-05 0.01019  3.77159E-03 0.00125 ];
INF_SP6                   (idx, [1:   8]) = [  4.39361E-03 0.00027 -1.55008E-04 0.00081 -2.98868E-05 0.00469 -1.31624E-02 0.00039 ];
INF_SP7                   (idx, [1:   8]) = [  5.50485E-04 0.00128 -1.47662E-04 0.00092 -3.17351E-05 0.00398  1.81095E-03 0.00260 ];

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

CMM_TRANSPXS              (idx, [1:   4]) = [  3.05442E-01 5.6E-05  9.12860E-01 0.00051 ];
CMM_TRANSPXS_X            (idx, [1:   4]) = [  3.22378E-01 7.9E-05  9.21348E-01 0.00054 ];
CMM_TRANSPXS_Y            (idx, [1:   4]) = [  3.22338E-01 6.5E-05  9.21476E-01 0.00062 ];
CMM_TRANSPXS_Z            (idx, [1:   4]) = [  2.76429E-01 8.5E-05  8.96235E-01 0.00060 ];
CMM_DIFFCOEF              (idx, [1:   4]) = [  1.09132E+00 5.6E-05  3.65157E-01 0.00051 ];
CMM_DIFFCOEF_X            (idx, [1:   4]) = [  1.03398E+00 7.9E-05  3.61794E-01 0.00054 ];
CMM_DIFFCOEF_Y            (idx, [1:   4]) = [  1.03411E+00 6.5E-05  3.61745E-01 0.00062 ];
CMM_DIFFCOEF_Z            (idx, [1:   4]) = [  1.20586E+00 8.5E-05  3.71933E-01 0.00060 ];

% Delayed neutron parameters (Meulekamp method):

BETA_EFF                  (idx, [1:  18]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
LAMBDA                    (idx, [1:  18]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
