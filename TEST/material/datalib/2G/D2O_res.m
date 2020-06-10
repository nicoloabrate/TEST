
% Increase counter:

if (exist('idx', 'var'));
  idx = idx + 1;
else;
  idx = 1;
end;

% Version, title and date:

VERSION                   (idx, [1: 14])  = 'Serpent 2.1.30' ;
COMPILE_DATE              (idx, [1: 20])  = 'May 22 2018 17:57:40' ;
DEBUG                     (idx, 1)        = 0 ;
TITLE                     (idx, [1:  8])  = 'Untitled' ;
CONFIDENTIAL_DATA         (idx, 1)        = 0 ;
INPUT_FILE_NAME           (idx, [1:  8])  = 'slab_d2o' ;
WORKING_DIRECTORY         (idx, [1: 47])  = '/home/caron/shared_memory/script_serpent/PHYTRA' ;
HOSTNAME                  (idx, [1:  4])  = 'case' ;
CPU_TYPE                  (idx, [1: 39])  = 'Intel(R) Core(TM) i7-4770 CPU @ 3.40GHz' ;
CPU_MHZ                   (idx, 1)        = 34.0 ;
START_DATE                (idx, [1: 24])  = 'Thu May 24 17:21:22 2018' ;
COMPLETE_DATE             (idx, [1: 24])  = 'Thu May 24 19:34:47 2018' ;

% Run parameters:

POP                       (idx, 1)        = 100000 ;
CYCLES                    (idx, 1)        = 1000 ;
SKIP                      (idx, 1)        = 500 ;
BATCH_INTERVAL            (idx, 1)        = 1 ;
SRC_NORM_MODE             (idx, 1)        = 2 ;
SEED                      (idx, 1)        = 1527175282 ;
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
OMP_THREADS               (idx, 1)        = 4 ;
MPI_REPRODUCIBILITY       (idx, 1)        = 0 ;
OMP_REPRODUCIBILITY       (idx, 1)        = 1 ;
OMP_HISTORY_PROFILE       (idx, [1:   4]) = [  9.87767E-01  1.01297E+00  1.00591E+00  9.93348E-01  ];
SHARE_BUF_ARRAY           (idx, 1)        = 0 ;
SHARE_RES2_ARRAY          (idx, 1)        = 1 ;

% File paths:

XS_DATA_FILE_PATH         (idx, [1: 46])  = '/opt/serpent/xsdata/jeff311/sss_jeff311.xsdata' ;
DECAY_DATA_FILE_PATH      (idx, [1:  3])  = 'N/A' ;
SFY_DATA_FILE_PATH        (idx, [1:  3])  = 'N/A' ;
NFY_DATA_FILE_PATH        (idx, [1:  3])  = 'N/A' ;
BRA_DATA_FILE_PATH        (idx, [1:  3])  = 'N/A' ;

% Collision and reaction sampling (neutrons/photons):

MIN_MACROXS               (idx, [1:   4]) = [  5.00000E-02 1.9E-09  0.00000E+00 0.0E+00 ];
DT_THRESH                 (idx, [1:  2])  = [  9.00000E-01  9.00000E-01 ];
ST_FRAC                   (idx, [1:   4]) = [  1.11150E-03 0.00030  0.00000E+00 0.0E+00 ];
DT_FRAC                   (idx, [1:   4]) = [  9.98888E-01 3.3E-07  0.00000E+00 0.0E+00 ];
DT_EFF                    (idx, [1:   4]) = [  6.22161E-01 3.0E-05  0.00000E+00 0.0E+00 ];
REA_SAMPLING_EFF          (idx, [1:   4]) = [  1.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
REA_SAMPLING_FAIL         (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
TOT_COL_EFF               (idx, [1:   4]) = [  6.22563E-01 3.0E-05  0.00000E+00 0.0E+00 ];
AVG_TRACKING_LOOPS        (idx, [1:   8]) = [  2.36595E+00 5.1E-05  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
AVG_TRACKS                (idx, [1:   4]) = [  8.93460E+01 0.00019  0.00000E+00 0.0E+00 ];
AVG_REAL_COL              (idx, [1:   4]) = [  8.92797E+01 0.00019  0.00000E+00 0.0E+00 ];
AVG_VIRT_COL              (idx, [1:   4]) = [  5.41278E+01 0.00026  0.00000E+00 0.0E+00 ];
AVG_SURF_CROSS            (idx, [1:   4]) = [  7.26528E-02 0.00037  0.00000E+00 0.0E+00 ];
LOST_PARTICLES            (idx, 1)        = 0 ;

% Run statistics:

CYCLE_IDX                 (idx, 1)        = 1000 ;
SOURCE_POPULATION         (idx, 1)        = 100001134 ;
MEAN_POP_SIZE             (idx, [1:  2])  = [  1.00001E+05 0.00017 ];
MEAN_POP_WGT              (idx, [1:  2])  = [  1.00001E+05 0.00017 ];
SIMULATION_COMPLETED      (idx, 1)        = 1 ;

% Running times:

TOT_CPU_TIME              (idx, 1)        =  5.25935E+02 ;
RUNNING_TIME              (idx, 1)        =  1.33415E+02 ;
INIT_TIME                 (idx, [1:  2])  = [  2.64000E-02  2.64000E-02 ];
PROCESS_TIME              (idx, [1:  2])  = [  8.33350E-05  8.33350E-05 ];
TRANSPORT_CYCLE_TIME      (idx, [1:  3])  = [  1.33388E+02  1.33388E+02  0.00000E+00 ];
MPI_OVERHEAD_TIME         (idx, [1:  2])  = [  0.00000E+00  0.00000E+00 ];
ESTIMATED_RUNNING_TIME    (idx, [1:  2])  = [  1.33414E+02  0.00000E+00 ];
CPU_USAGE                 (idx, 1)        = 3.94211 ;
TRANSPORT_CPU_USAGE       (idx, [1:   2]) = [  3.94995E+00 9.1E-05 ];
OMP_PARALLEL_FRAC         (idx, 1)        =  9.84206E-01 ;

% Memory usage:

AVAIL_MEM                 (idx, 1)        = 31809.05 ;
ALLOC_MEMSIZE             (idx, 1)        = 776.29;
MEMSIZE                   (idx, 1)        = 707.11;
XS_MEMSIZE                (idx, 1)        = 43.47;
MAT_MEMSIZE               (idx, 1)        = 8.98;
RES_MEMSIZE               (idx, 1)        = 2.15;
MISC_MEMSIZE              (idx, 1)        = 652.51;
UNKNOWN_MEMSIZE           (idx, 1)        = 0.00;
UNUSED_MEMSIZE            (idx, 1)        = 69.18;

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

NORM_COEF                 (idx, [1:   4]) = [  9.95413E-06 8.4E-05  0.00000E+00 0.0E+00 ];

% Analog reaction rate estimators:

CONVERSION_RATIO          (idx, [1:   2]) = [  1.07340E+00 0.00021 ];
U235_FISS                 (idx, [1:   4]) = [  3.32383E-01 0.00016  8.08442E-01 7.6E-05 ];
U238_FISS                 (idx, [1:   4]) = [  7.87581E-02 0.00037  1.91558E-01 0.00032 ];
U235_CAPT                 (idx, [1:   4]) = [  7.76317E-02 0.00036  1.48541E-01 0.00033 ];
U238_CAPT                 (idx, [1:   4]) = [  4.40228E-01 0.00014  8.42337E-01 6.0E-05 ];

% Neutron balance (particles/weight):

BALA_SRC_NEUTRON_SRC     (idx, [1:  2])  = [ 0 0.00000E+00 ];
BALA_SRC_NEUTRON_FISS    (idx, [1:  2])  = [ 100001134 1.00000E+08 ];
BALA_SRC_NEUTRON_NXN     (idx, [1:  2])  = [ 0 4.72037E+05 ];
BALA_SRC_NEUTRON_VR      (idx, [1:  2])  = [ 0 0.00000E+00 ];
BALA_SRC_NEUTRON_TOT     (idx, [1:  2])  = [ 100001134 1.00472E+08 ];

BALA_LOSS_NEUTRON_CAPT    (idx, [1:  2])  = [ 52240768 5.25036E+07 ];
BALA_LOSS_NEUTRON_FISS    (idx, [1:  2])  = [ 41133201 4.13036E+07 ];
BALA_LOSS_NEUTRON_LEAK    (idx, [1:  2])  = [ 6627165 6.66489E+06 ];
BALA_LOSS_NEUTRON_CUT     (idx, [1:  2])  = [ 0 0.00000E+00 ];
BALA_LOSS_NEUTRON_ERR     (idx, [1:  2])  = [ 0 0.00000E+00 ];
BALA_LOSS_NEUTRON_TOT     (idx, [1:  2])  = [ 100001134 1.00472E+08 ];

BALA_NEUTRON_DIFF         (idx, [1:  2])  = [ 0 3.72529E-06 ];

% Normalized total reaction rates (neutrons):

TOT_POWER                 (idx, [1:   2]) = [  1.33702E-11 6.9E-05 ];
TOT_POWDENS               (idx, [1:   2]) = [  0.00000E+00 0.0E+00 ];
TOT_GENRATE               (idx, [1:   2]) = [  1.03475E+00 6.9E-05 ];
TOT_FISSRATE              (idx, [1:   2]) = [  4.11107E-01 6.9E-05 ];
TOT_CAPTRATE              (idx, [1:   2]) = [  5.22550E-01 6.3E-05 ];
TOT_ABSRATE               (idx, [1:   2]) = [  9.33657E-01 2.8E-05 ];
TOT_SRCRATE               (idx, [1:   2]) = [  9.95413E-01 8.4E-05 ];
TOT_FLUX                  (idx, [1:   2]) = [  2.30324E+02 0.00014 ];
TOT_PHOTON_PRODRATE       (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
TOT_LEAKRATE              (idx, [1:   2]) = [  6.63429E-02 0.00040 ];
ALBEDO_LEAKRATE           (idx, [1:   2]) = [  0.00000E+00 0.0E+00 ];
TOT_LOSSRATE              (idx, [1:   2]) = [  1.00000E+00 0.0E+00 ];
TOT_CUTRATE               (idx, [1:   2]) = [  0.00000E+00 0.0E+00 ];
TOT_RR                    (idx, [1:   2]) = [  8.93341E+01 0.00018 ];
INI_FMASS                 (idx, 1)        =  0.00000E+00 ;
TOT_FMASS                 (idx, 1)        =  0.00000E+00 ;

% Six-factor formula:

SIX_FF_ETA                (idx, [1:   2]) = [  1.89104E+00 0.00012 ];
SIX_FF_F                  (idx, [1:   2]) = [  9.80878E-01 3.1E-05 ];
SIX_FF_P                  (idx, [1:   2]) = [  2.20999E-01 0.00021 ];
SIX_FF_EPSILON            (idx, [1:   2]) = [  2.71730E+00 0.00022 ];
SIX_FF_LF                 (idx, [1:   2]) = [  9.99704E-01 1.7E-06 ];
SIX_FF_LT                 (idx, [1:   2]) = [  9.33627E-01 2.8E-05 ];
SIX_FF_KINF               (idx, [1:   2]) = [  1.11384E+00 0.00012 ];
SIX_FF_KEFF               (idx, [1:   2]) = [  1.03961E+00 0.00012 ];

% Fission neutron and energy production:

NUBAR                     (idx, [1:   2]) = [  2.51698E+00 7.1E-06 ];
FISSE                     (idx, [1:   2]) = [  2.02989E+02 5.9E-07 ];

% Criticality eigenvalues:

ANA_KEFF                  (idx, [1:   6]) = [  1.03963E+00 0.00012  1.03224E+00 0.00012  7.36933E-03 0.00176 ];
IMP_KEFF                  (idx, [1:   2]) = [  1.03964E+00 7.0E-05 ];
COL_KEFF                  (idx, [1:   2]) = [  1.03952E+00 9.6E-05 ];
ABS_KEFF                  (idx, [1:   2]) = [  1.03964E+00 7.0E-05 ];
ABS_KINF                  (idx, [1:   2]) = [  1.11389E+00 6.8E-05 ];
GEOM_ALBEDO               (idx, [1:   6]) = [  1.00000E+00 0.0E+00  1.00000E+00 0.0E+00  1.00000E+00 0.0E+00 ];

% ALF (Average lethargy of neutrons causing fission):
% Based on E0 = 2.000000E+01 MeV

ANA_ALF                   (idx, [1:   2]) = [  1.09306E+01 0.00012 ];
IMP_ALF                   (idx, [1:   2]) = [  1.09305E+01 9.0E-05 ];

% EALF (Energy corresponding to average lethargy of neutrons causing fission):

ANA_EALF                  (idx, [1:   2]) = [  3.58329E-04 0.00128 ];
IMP_EALF                  (idx, [1:   2]) = [  3.58243E-04 0.00099 ];

% AFGE (Average energy of neutrons causing fission):

ANA_AFGE                  (idx, [1:   2]) = [  7.07842E-01 0.00032 ];
IMP_AFGE                  (idx, [1:   2]) = [  7.07711E-01 0.00019 ];

% Forward-weighted delayed neutron parameters:

FWD_ANA_BETA_ZERO         (idx, [1:  18]) = [  8.44270E-03 0.00111  1.98916E-04 0.00683  1.11990E-03 0.00289  6.10968E-04 0.00405  1.46853E-03 0.00263  2.61659E-03 0.00190  1.13107E-03 0.00288  8.58977E-04 0.00336  4.37737E-04 0.00464 ];
FWD_ANA_LAMBDA            (idx, [1:  18]) = [  5.60803E-01 0.00158  1.24667E-02 7.2E-09  2.82917E-02 0.0E+00  4.25244E-02 0.0E+00  1.33042E-01 0.0E+00  2.92467E-01 0.0E+00  6.66488E-01 2.7E-09  1.63478E+00 4.2E-09  3.55460E+00 2.9E-09 ];

% Beta-eff using Meulekamp's method:

ADJ_MEULEKAMP_BETA_EFF    (idx, [1:  18]) = [  7.31477E-03 0.00165  1.68574E-04 0.01045  9.94792E-04 0.00453  5.19082E-04 0.00639  1.27531E-03 0.00398  2.27709E-03 0.00287  9.68238E-04 0.00443  7.35500E-04 0.00515  3.76186E-04 0.00722 ];
ADJ_MEULEKAMP_LAMBDA      (idx, [1:  18]) = [  5.56784E-01 0.00250  1.24667E-02 7.2E-09  2.82917E-02 0.0E+00  4.25244E-02 0.0E+00  1.33042E-01 0.0E+00  2.92467E-01 0.0E+00  6.66488E-01 2.7E-09  1.63478E+00 4.2E-09  3.55460E+00 3.0E-09 ];

% Adjoint weighted time constants using Nauchi's method:

ADJ_NAUCHI_GEN_TIME       (idx, [1:   6]) = [  4.62729E-04 0.00052  4.62264E-04 0.00052  5.27231E-04 0.00563 ];
ADJ_NAUCHI_LIFETIME       (idx, [1:   6]) = [  4.81063E-04 0.00052  4.80580E-04 0.00052  5.48124E-04 0.00562 ];
ADJ_NAUCHI_BETA_EFF       (idx, [1:  18]) = [  7.09057E-03 0.00181  1.64489E-04 0.01236  9.64144E-04 0.00501  5.02669E-04 0.00673  1.23407E-03 0.00463  2.21230E-03 0.00313  9.38182E-04 0.00526  7.11330E-04 0.00590  3.63388E-04 0.00820 ];
ADJ_NAUCHI_LAMBDA         (idx, [1:  18]) = [  5.56040E-01 0.00288  1.24667E-02 7.2E-09  2.82917E-02 0.0E+00  4.25244E-02 0.0E+00  1.33042E-01 0.0E+00  2.92467E-01 0.0E+00  6.66488E-01 2.7E-09  1.63478E+00 4.2E-09  3.55460E+00 2.8E-09 ];

% Adjoint weighted time constants using IFP:

ADJ_IFP_GEN_TIME          (idx, [1:   6]) = [  4.84149E-04 0.00129  4.83653E-04 0.00129  5.54303E-04 0.01423 ];
ADJ_IFP_LIFETIME          (idx, [1:   6]) = [  5.03333E-04 0.00129  5.02818E-04 0.00129  5.76273E-04 0.01423 ];
ADJ_IFP_IMP_BETA_EFF      (idx, [1:  18]) = [  7.13133E-03 0.00594  1.59640E-04 0.04160  9.77148E-04 0.01546  5.00680E-04 0.02235  1.25540E-03 0.01443  2.22579E-03 0.01085  9.38960E-04 0.01640  7.08185E-04 0.01823  3.65528E-04 0.02588 ];
ADJ_IFP_IMP_LAMBDA        (idx, [1:  18]) = [  5.54951E-01 0.00868  1.24667E-02 0.0E+00  2.82917E-02 0.0E+00  4.25244E-02 0.0E+00  1.33042E-01 0.0E+00  2.92467E-01 0.0E+00  6.66488E-01 2.7E-09  1.63478E+00 4.3E-09  3.55460E+00 3.0E-09 ];
ADJ_IFP_ANA_BETA_EFF      (idx, [1:  18]) = [  7.13522E-03 0.00582  1.60426E-04 0.04056  9.78173E-04 0.01518  4.96763E-04 0.02204  1.25811E-03 0.01421  2.22878E-03 0.01063  9.40140E-04 0.01603  7.08907E-04 0.01798  3.63928E-04 0.02525 ];
ADJ_IFP_ANA_LAMBDA        (idx, [1:  18]) = [  5.54217E-01 0.00850  1.24667E-02 0.0E+00  2.82917E-02 0.0E+00  4.25244E-02 0.0E+00  1.33042E-01 0.0E+00  2.92467E-01 0.0E+00  6.66488E-01 2.7E-09  1.63478E+00 4.3E-09  3.55460E+00 2.8E-09 ];
ADJ_IFP_ROSSI_ALPHA       (idx, [1:   2]) = [ -1.47712E+01 0.00610 ];

% Adjoint weighted time constants using perturbation technique:

ADJ_PERT_GEN_TIME         (idx, [1:   2]) = [  4.76062E-04 0.00029 ];
ADJ_PERT_LIFETIME         (idx, [1:   2]) = [  4.94923E-04 0.00026 ];
ADJ_PERT_BETA_EFF         (idx, [1:   2]) = [  7.17129E-03 0.00113 ];
ADJ_PERT_ROSSI_ALPHA      (idx, [1:   2]) = [ -1.50651E+01 0.00117 ];

% Inverse neutron speed :

ANA_INV_SPD               (idx, [1:   2]) = [  1.93990E-06 0.00016 ];

% Analog slowing-down and thermal neutron lifetime (total/prompt/delayed):

ANA_SLOW_TIME             (idx, [1:   6]) = [  1.07069E-05 0.00011  1.07067E-05 0.00011  1.07283E-05 0.00122 ];
ANA_THERM_TIME            (idx, [1:   6]) = [  1.63246E-03 0.00024  1.63278E-03 0.00024  1.59449E-03 0.00267 ];
ANA_THERM_FRAC            (idx, [1:   6]) = [  2.72744E-01 0.00018  2.72793E-01 0.00018  2.67304E-01 0.00243 ];
ANA_DELAYED_EMTIME        (idx, [1:   2]) = [  1.09133E+01 0.00239 ];
ANA_MEAN_NCOL             (idx, [1:   4]) = [  8.92797E+01 0.00019  8.82382E+01 0.00029 ];

% Group constant generation:

GC_UNIVERSE_NAME          (idx, [1:  3])  = 'UO2' ;

% Micro- and macro-group structures:

MICRO_NG                  (idx, 1)        = 70 ;
MICRO_E                   (idx, [1:  71]) = [  1.00000E-11  5.00000E-09  1.00000E-08  1.50000E-08  2.00000E-08  2.50000E-08  3.00000E-08  3.50000E-08  4.20000E-08  5.00000E-08  5.80000E-08  6.70000E-08  8.00000E-08  1.00000E-07  1.40000E-07  1.80000E-07  2.20000E-07  2.50000E-07  2.80000E-07  3.00000E-07  3.20000E-07  3.50000E-07  4.00000E-07  5.00000E-07  6.25000E-07  7.80000E-07  8.50000E-07  9.10000E-07  9.50000E-07  9.72000E-07  9.96000E-07  1.02000E-06  1.04500E-06  1.07100E-06  1.09700E-06  1.12300E-06  1.15000E-06  1.30000E-06  1.50000E-06  1.85500E-06  2.10000E-06  2.60000E-06  3.30000E-06  4.00000E-06  9.87700E-06  1.59680E-05  2.77000E-05  4.80520E-05  7.55014E-05  1.48728E-04  3.67262E-04  9.06898E-04  1.42510E-03  2.23945E-03  3.51910E-03  5.50000E-03  9.11800E-03  1.50300E-02  2.47800E-02  4.08500E-02  6.74300E-02  1.11000E-01  1.83000E-01  3.02500E-01  5.00000E-01  8.21000E-01  1.35300E+00  2.23100E+00  3.67900E+00  6.06550E+00  2.00000E+01 ];

MACRO_NG                  (idx, 1)        = 2 ;
MACRO_E                   (idx, [1:   3]) = [  1.00000E+37  6.25000E-07  0.00000E+00 ];

% Micro-group spectrum:

INF_MICRO_FLX             (idx, [1: 140]) = [  5.12903E+05 0.00077  2.15587E+06 0.00040  4.14984E+06 0.00029  5.49518E+06 0.00030  9.59953E+06 0.00028  2.09049E+07 0.00031  2.82690E+07 0.00028  2.73772E+07 0.00026  2.25262E+07 0.00029  1.72669E+07 0.00029  1.29182E+07 0.00034  8.05041E+06 0.00036  4.24789E+06 0.00057  2.44728E+06 0.00078  1.64414E+06 0.00068  1.08827E+06 0.00082  9.19887E+05 0.00102  8.25432E+05 0.00116  6.56261E+05 0.00126  9.97733E+05 0.00099  7.06215E+05 0.00077  3.77987E+05 0.00114  1.94681E+05 0.00125  1.99887E+05 0.00152  1.66712E+05 0.00151  1.90281E+05 0.00148  2.26390E+05 0.00107  8.41049E+04 0.00203  1.03740E+05 0.00222  1.29916E+05 0.00198  5.12281E+04 0.00248  1.14243E+05 0.00227  7.17605E+04 0.00285  3.83089E+04 0.00298  4.20346E+03 0.00437  3.97593E+03 0.00457  4.25475E+03 0.00450  4.68115E+03 0.00411  4.99040E+03 0.00442  5.24538E+03 0.00415  5.72031E+03 0.00471  5.63653E+03 0.00468  1.09904E+04 0.00424  1.82321E+04 0.00341  2.33466E+04 0.00287  5.83683E+04 0.00174  5.17107E+04 0.00203  4.17441E+04 0.00220  1.97346E+04 0.00252  1.10919E+04 0.00280  7.26540E+03 0.00355  7.49076E+03 0.00303  1.26165E+04 0.00210  1.60562E+04 0.00201  3.16711E+04 0.00155  5.71401E+04 0.00126  1.15550E+05 0.00079  9.62176E+04 0.00081  8.03872E+04 0.00067  6.25136E+04 0.00073  5.97990E+04 0.00102  6.14925E+04 0.00079  5.30182E+04 0.00072  3.57587E+04 0.00094  3.24928E+04 0.00089  2.82025E+04 0.00086  2.24821E+04 0.00115  1.54125E+04 0.00098  8.15979E+03 0.00142  1.79475E+03 0.00192 ];

% Integral parameters:

INF_KINF                  (idx, [1:   2]) = [  1.11833E+00 0.00012 ];

% Flux spectra in infinite geometry:

INF_FLX                   (idx, [1:   4]) = [  8.70280E+01 0.00017  4.62769E-01 0.00032 ];
INF_FISS_FLX              (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Reaction cross sections:

INF_TOT                   (idx, [1:   4]) = [  2.68659E-01 2.0E-05  6.81228E-01 9.6E-05 ];
INF_CAPT                  (idx, [1:   4]) = [  5.43183E-03 9.8E-05  9.74174E-02 0.00014 ];
INF_ABS                   (idx, [1:   4]) = [  8.35978E-03 9.5E-05  4.35158E-01 0.00015 ];
INF_FISS                  (idx, [1:   4]) = [  2.92795E-03 0.00012  3.37741E-01 0.00015 ];
INF_NSF                   (idx, [1:   4]) = [  7.51466E-03 0.00012  8.22804E-01 0.00015 ];
INF_NUBAR                 (idx, [1:   4]) = [  2.56652E+00 9.3E-06  2.43620E+00 0.0E+00 ];
INF_KAPPA                 (idx, [1:   4]) = [  2.03430E+02 7.8E-07  2.02270E+02 4.2E-09 ];
INF_INVV                  (idx, [1:   4]) = [  6.00193E-09 0.00046  2.73688E-06 0.00013 ];

% Total scattering cross sections:

INF_SCATT0                (idx, [1:   4]) = [  2.60298E-01 2.2E-05  2.46028E-01 0.00020 ];
INF_SCATT1                (idx, [1:   4]) = [  4.10270E-02 0.00010  7.49803E-04 0.05249 ];
INF_SCATT2                (idx, [1:   4]) = [  1.27072E-02 0.00024 -1.81332E-05 1.00000 ];
INF_SCATT3                (idx, [1:   4]) = [  5.94472E-03 0.00033 -8.77279E-06 1.00000 ];
INF_SCATT4                (idx, [1:   4]) = [  3.22449E-03 0.00066  2.20651E-05 0.94989 ];
INF_SCATT5                (idx, [1:   4]) = [  1.47865E-03 0.00129  7.59642E-06 1.00000 ];
INF_SCATT6                (idx, [1:   4]) = [  7.54602E-04 0.00171 -1.27011E-05 1.00000 ];
INF_SCATT7                (idx, [1:   4]) = [  3.20895E-04 0.00433  1.49691E-05 1.00000 ];

% Total scattering production cross sections:

INF_SCATTP0               (idx, [1:   4]) = [  2.60341E-01 2.2E-05  2.46028E-01 0.00020 ];
INF_SCATTP1               (idx, [1:   4]) = [  4.10277E-02 0.00010  7.49803E-04 0.05249 ];
INF_SCATTP2               (idx, [1:   4]) = [  1.27074E-02 0.00024 -1.81332E-05 1.00000 ];
INF_SCATTP3               (idx, [1:   4]) = [  5.94474E-03 0.00033 -8.77279E-06 1.00000 ];
INF_SCATTP4               (idx, [1:   4]) = [  3.22453E-03 0.00066  2.20651E-05 0.94989 ];
INF_SCATTP5               (idx, [1:   4]) = [  1.47866E-03 0.00129  7.59642E-06 1.00000 ];
INF_SCATTP6               (idx, [1:   4]) = [  7.54595E-04 0.00171 -1.27011E-05 1.00000 ];
INF_SCATTP7               (idx, [1:   4]) = [  3.20895E-04 0.00433  1.49691E-05 1.00000 ];

% Diffusion parameters:

INF_TRANSPXS              (idx, [1:   4]) = [  2.03163E-01 4.0E-05  6.07966E-01 0.00015 ];
INF_DIFFCOEF              (idx, [1:   4]) = [  1.64072E+00 4.0E-05  5.48277E-01 0.00015 ];

% Reduced absoption and removal:

INF_RABSXS                (idx, [1:   4]) = [  8.31690E-03 9.6E-05  4.35158E-01 0.00015 ];
INF_REMXS                 (idx, [1:   4]) = [  8.36559E-03 0.00013  4.35622E-01 0.00017 ];

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

INF_CHIT                  (idx, [1:   4]) = [  1.00000E+00 9.8E-09  9.92803E-09 1.00000 ];
INF_CHIP                  (idx, [1:   4]) = [  1.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_CHID                  (idx, [1:   4]) = [  9.99999E-01 1.1E-06  1.14226E-06 1.00000 ];

% Scattering matrixes:

INF_S0                    (idx, [1:   8]) = [  2.60293E-01 2.2E-05  5.01586E-06 0.00436  4.22090E-04 0.00664  2.45606E-01 0.00020 ];
INF_S1                    (idx, [1:   8]) = [  4.10283E-02 0.00010 -1.24683E-06 0.01091 -5.69498E-05 0.02549  8.06753E-04 0.04868 ];
INF_S2                    (idx, [1:   8]) = [  1.27073E-02 0.00024 -7.67013E-08 0.15818 -1.79031E-05 0.07165 -2.30107E-07 1.00000 ];
INF_S3                    (idx, [1:   8]) = [  5.94478E-03 0.00033 -5.21196E-08 0.18701 -4.94815E-06 0.24254 -3.82464E-06 1.00000 ];
INF_S4                    (idx, [1:   8]) = [  3.22450E-03 0.00066 -9.84009E-09 0.78180 -2.85968E-06 0.32152  2.49248E-05 0.84588 ];
INF_S5                    (idx, [1:   8]) = [  1.47866E-03 0.00129 -7.56944E-09 0.88777 -3.33505E-06 0.27585  1.09315E-05 1.00000 ];
INF_S6                    (idx, [1:   8]) = [  7.54607E-04 0.00172 -5.79097E-09 1.00000 -1.33458E-06 0.63165 -1.13665E-05 1.00000 ];
INF_S7                    (idx, [1:   8]) = [  3.20892E-04 0.00433  2.86907E-09 1.00000  3.52480E-07 1.00000  1.46166E-05 1.00000 ];

% Scattering production matrixes:

INF_SP0                   (idx, [1:   8]) = [  2.60336E-01 2.2E-05  5.01586E-06 0.00436  4.22090E-04 0.00664  2.45606E-01 0.00020 ];
INF_SP1                   (idx, [1:   8]) = [  4.10289E-02 0.00010 -1.24683E-06 0.01091 -5.69498E-05 0.02549  8.06753E-04 0.04868 ];
INF_SP2                   (idx, [1:   8]) = [  1.27075E-02 0.00024 -7.67013E-08 0.15818 -1.79031E-05 0.07165 -2.30107E-07 1.00000 ];
INF_SP3                   (idx, [1:   8]) = [  5.94479E-03 0.00033 -5.21196E-08 0.18701 -4.94815E-06 0.24254 -3.82464E-06 1.00000 ];
INF_SP4                   (idx, [1:   8]) = [  3.22454E-03 0.00066 -9.84009E-09 0.78180 -2.85968E-06 0.32152  2.49248E-05 0.84588 ];
INF_SP5                   (idx, [1:   8]) = [  1.47867E-03 0.00129 -7.56944E-09 0.88777 -3.33505E-06 0.27585  1.09315E-05 1.00000 ];
INF_SP6                   (idx, [1:   8]) = [  7.54600E-04 0.00172 -5.79097E-09 1.00000 -1.33458E-06 0.63165 -1.13665E-05 1.00000 ];
INF_SP7                   (idx, [1:   8]) = [  3.20892E-04 0.00433  2.86907E-09 1.00000  3.52480E-07 1.00000  1.46166E-05 1.00000 ];

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

CMM_TRANSPXS              (idx, [1:   4]) = [  2.39525E-01 0.00019  2.46574E-03 0.00032 ];
CMM_TRANSPXS_X            (idx, [1:   4]) = [  2.85282E-01 0.00029  1.22569E-02 0.00094 ];
CMM_TRANSPXS_Y            (idx, [1:   4]) = [  2.21770E-01 0.00023  1.76163E-03 0.00046 ];
CMM_TRANSPXS_Z            (idx, [1:   4]) = [  2.21715E-01 0.00027  1.76236E-03 0.00048 ];
CMM_DIFFCOEF              (idx, [1:   4]) = [  1.39165E+00 0.00019  1.35187E+02 0.00032 ];
CMM_DIFFCOEF_X            (idx, [1:   4]) = [  1.16844E+00 0.00029  2.71968E+01 0.00094 ];
CMM_DIFFCOEF_Y            (idx, [1:   4]) = [  1.50306E+00 0.00023  1.89221E+02 0.00046 ];
CMM_DIFFCOEF_Z            (idx, [1:   4]) = [  1.50343E+00 0.00027  1.89143E+02 0.00049 ];

% Delayed neutron parameters (Meulekamp method):

BETA_EFF                  (idx, [1:  18]) = [  7.31477E-03 0.00165  1.68574E-04 0.01045  9.94792E-04 0.00453  5.19082E-04 0.00639  1.27531E-03 0.00398  2.27709E-03 0.00287  9.68238E-04 0.00443  7.35500E-04 0.00515  3.76186E-04 0.00722 ];
LAMBDA                    (idx, [1:  18]) = [  5.56784E-01 0.00250  1.24667E-02 7.2E-09  2.82917E-02 0.0E+00  4.25244E-02 0.0E+00  1.33042E-01 0.0E+00  2.92467E-01 0.0E+00  6.66488E-01 2.7E-09  1.63478E+00 4.2E-09  3.55460E+00 3.0E-09 ];


% Increase counter:

if (exist('idx', 'var'));
  idx = idx + 1;
else;
  idx = 1;
end;

% Version, title and date:

VERSION                   (idx, [1: 14])  = 'Serpent 2.1.30' ;
COMPILE_DATE              (idx, [1: 20])  = 'May 22 2018 17:57:40' ;
DEBUG                     (idx, 1)        = 0 ;
TITLE                     (idx, [1:  8])  = 'Untitled' ;
CONFIDENTIAL_DATA         (idx, 1)        = 0 ;
INPUT_FILE_NAME           (idx, [1:  8])  = 'slab_d2o' ;
WORKING_DIRECTORY         (idx, [1: 47])  = '/home/caron/shared_memory/script_serpent/PHYTRA' ;
HOSTNAME                  (idx, [1:  4])  = 'case' ;
CPU_TYPE                  (idx, [1: 39])  = 'Intel(R) Core(TM) i7-4770 CPU @ 3.40GHz' ;
CPU_MHZ                   (idx, 1)        = 34.0 ;
START_DATE                (idx, [1: 24])  = 'Thu May 24 17:21:22 2018' ;
COMPLETE_DATE             (idx, [1: 24])  = 'Thu May 24 19:34:47 2018' ;

% Run parameters:

POP                       (idx, 1)        = 100000 ;
CYCLES                    (idx, 1)        = 1000 ;
SKIP                      (idx, 1)        = 500 ;
BATCH_INTERVAL            (idx, 1)        = 1 ;
SRC_NORM_MODE             (idx, 1)        = 2 ;
SEED                      (idx, 1)        = 1527175282 ;
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
OMP_THREADS               (idx, 1)        = 4 ;
MPI_REPRODUCIBILITY       (idx, 1)        = 0 ;
OMP_REPRODUCIBILITY       (idx, 1)        = 1 ;
OMP_HISTORY_PROFILE       (idx, [1:   4]) = [  9.87767E-01  1.01297E+00  1.00591E+00  9.93348E-01  ];
SHARE_BUF_ARRAY           (idx, 1)        = 0 ;
SHARE_RES2_ARRAY          (idx, 1)        = 1 ;

% File paths:

XS_DATA_FILE_PATH         (idx, [1: 46])  = '/opt/serpent/xsdata/jeff311/sss_jeff311.xsdata' ;
DECAY_DATA_FILE_PATH      (idx, [1:  3])  = 'N/A' ;
SFY_DATA_FILE_PATH        (idx, [1:  3])  = 'N/A' ;
NFY_DATA_FILE_PATH        (idx, [1:  3])  = 'N/A' ;
BRA_DATA_FILE_PATH        (idx, [1:  3])  = 'N/A' ;

% Collision and reaction sampling (neutrons/photons):

MIN_MACROXS               (idx, [1:   4]) = [  5.00000E-02 1.9E-09  0.00000E+00 0.0E+00 ];
DT_THRESH                 (idx, [1:  2])  = [  9.00000E-01  9.00000E-01 ];
ST_FRAC                   (idx, [1:   4]) = [  1.11150E-03 0.00030  0.00000E+00 0.0E+00 ];
DT_FRAC                   (idx, [1:   4]) = [  9.98888E-01 3.3E-07  0.00000E+00 0.0E+00 ];
DT_EFF                    (idx, [1:   4]) = [  6.22161E-01 3.0E-05  0.00000E+00 0.0E+00 ];
REA_SAMPLING_EFF          (idx, [1:   4]) = [  1.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
REA_SAMPLING_FAIL         (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
TOT_COL_EFF               (idx, [1:   4]) = [  6.22563E-01 3.0E-05  0.00000E+00 0.0E+00 ];
AVG_TRACKING_LOOPS        (idx, [1:   8]) = [  2.36595E+00 5.1E-05  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
AVG_TRACKS                (idx, [1:   4]) = [  8.93460E+01 0.00019  0.00000E+00 0.0E+00 ];
AVG_REAL_COL              (idx, [1:   4]) = [  8.92797E+01 0.00019  0.00000E+00 0.0E+00 ];
AVG_VIRT_COL              (idx, [1:   4]) = [  5.41278E+01 0.00026  0.00000E+00 0.0E+00 ];
AVG_SURF_CROSS            (idx, [1:   4]) = [  7.26528E-02 0.00037  0.00000E+00 0.0E+00 ];
LOST_PARTICLES            (idx, 1)        = 0 ;

% Run statistics:

CYCLE_IDX                 (idx, 1)        = 1000 ;
SOURCE_POPULATION         (idx, 1)        = 100001134 ;
MEAN_POP_SIZE             (idx, [1:  2])  = [  1.00001E+05 0.00017 ];
MEAN_POP_WGT              (idx, [1:  2])  = [  1.00001E+05 0.00017 ];
SIMULATION_COMPLETED      (idx, 1)        = 1 ;

% Running times:

TOT_CPU_TIME              (idx, 1)        =  5.25935E+02 ;
RUNNING_TIME              (idx, 1)        =  1.33415E+02 ;
INIT_TIME                 (idx, [1:  2])  = [  2.64000E-02  2.64000E-02 ];
PROCESS_TIME              (idx, [1:  2])  = [  8.33350E-05  8.33350E-05 ];
TRANSPORT_CYCLE_TIME      (idx, [1:  3])  = [  1.33388E+02  1.33388E+02  0.00000E+00 ];
MPI_OVERHEAD_TIME         (idx, [1:  2])  = [  0.00000E+00  0.00000E+00 ];
ESTIMATED_RUNNING_TIME    (idx, [1:  2])  = [  1.33414E+02  0.00000E+00 ];
CPU_USAGE                 (idx, 1)        = 3.94211 ;
TRANSPORT_CPU_USAGE       (idx, [1:   2]) = [  3.94995E+00 9.1E-05 ];
OMP_PARALLEL_FRAC         (idx, 1)        =  9.84206E-01 ;

% Memory usage:

AVAIL_MEM                 (idx, 1)        = 31809.05 ;
ALLOC_MEMSIZE             (idx, 1)        = 776.29;
MEMSIZE                   (idx, 1)        = 707.11;
XS_MEMSIZE                (idx, 1)        = 43.47;
MAT_MEMSIZE               (idx, 1)        = 8.98;
RES_MEMSIZE               (idx, 1)        = 2.15;
MISC_MEMSIZE              (idx, 1)        = 652.51;
UNKNOWN_MEMSIZE           (idx, 1)        = 0.00;
UNUSED_MEMSIZE            (idx, 1)        = 69.18;

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

NORM_COEF                 (idx, [1:   4]) = [  9.95413E-06 8.4E-05  0.00000E+00 0.0E+00 ];

% Analog reaction rate estimators:

CONVERSION_RATIO          (idx, [1:   2]) = [  1.07340E+00 0.00021 ];
U235_FISS                 (idx, [1:   4]) = [  3.32383E-01 0.00016  8.08442E-01 7.6E-05 ];
U238_FISS                 (idx, [1:   4]) = [  7.87581E-02 0.00037  1.91558E-01 0.00032 ];
U235_CAPT                 (idx, [1:   4]) = [  7.76317E-02 0.00036  1.48541E-01 0.00033 ];
U238_CAPT                 (idx, [1:   4]) = [  4.40228E-01 0.00014  8.42337E-01 6.0E-05 ];

% Neutron balance (particles/weight):

BALA_SRC_NEUTRON_SRC     (idx, [1:  2])  = [ 0 0.00000E+00 ];
BALA_SRC_NEUTRON_FISS    (idx, [1:  2])  = [ 100001134 1.00000E+08 ];
BALA_SRC_NEUTRON_NXN     (idx, [1:  2])  = [ 0 4.72037E+05 ];
BALA_SRC_NEUTRON_VR      (idx, [1:  2])  = [ 0 0.00000E+00 ];
BALA_SRC_NEUTRON_TOT     (idx, [1:  2])  = [ 100001134 1.00472E+08 ];

BALA_LOSS_NEUTRON_CAPT    (idx, [1:  2])  = [ 52240768 5.25036E+07 ];
BALA_LOSS_NEUTRON_FISS    (idx, [1:  2])  = [ 41133201 4.13036E+07 ];
BALA_LOSS_NEUTRON_LEAK    (idx, [1:  2])  = [ 6627165 6.66489E+06 ];
BALA_LOSS_NEUTRON_CUT     (idx, [1:  2])  = [ 0 0.00000E+00 ];
BALA_LOSS_NEUTRON_ERR     (idx, [1:  2])  = [ 0 0.00000E+00 ];
BALA_LOSS_NEUTRON_TOT     (idx, [1:  2])  = [ 100001134 1.00472E+08 ];

BALA_NEUTRON_DIFF         (idx, [1:  2])  = [ 0 3.72529E-06 ];

% Normalized total reaction rates (neutrons):

TOT_POWER                 (idx, [1:   2]) = [  1.33702E-11 6.9E-05 ];
TOT_POWDENS               (idx, [1:   2]) = [  0.00000E+00 0.0E+00 ];
TOT_GENRATE               (idx, [1:   2]) = [  1.03475E+00 6.9E-05 ];
TOT_FISSRATE              (idx, [1:   2]) = [  4.11107E-01 6.9E-05 ];
TOT_CAPTRATE              (idx, [1:   2]) = [  5.22550E-01 6.3E-05 ];
TOT_ABSRATE               (idx, [1:   2]) = [  9.33657E-01 2.8E-05 ];
TOT_SRCRATE               (idx, [1:   2]) = [  9.95413E-01 8.4E-05 ];
TOT_FLUX                  (idx, [1:   2]) = [  2.30324E+02 0.00014 ];
TOT_PHOTON_PRODRATE       (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
TOT_LEAKRATE              (idx, [1:   2]) = [  6.63429E-02 0.00040 ];
ALBEDO_LEAKRATE           (idx, [1:   2]) = [  0.00000E+00 0.0E+00 ];
TOT_LOSSRATE              (idx, [1:   2]) = [  1.00000E+00 0.0E+00 ];
TOT_CUTRATE               (idx, [1:   2]) = [  0.00000E+00 0.0E+00 ];
TOT_RR                    (idx, [1:   2]) = [  8.93341E+01 0.00018 ];
INI_FMASS                 (idx, 1)        =  0.00000E+00 ;
TOT_FMASS                 (idx, 1)        =  0.00000E+00 ;

% Six-factor formula:

SIX_FF_ETA                (idx, [1:   2]) = [  1.89104E+00 0.00012 ];
SIX_FF_F                  (idx, [1:   2]) = [  9.80878E-01 3.1E-05 ];
SIX_FF_P                  (idx, [1:   2]) = [  2.20999E-01 0.00021 ];
SIX_FF_EPSILON            (idx, [1:   2]) = [  2.71730E+00 0.00022 ];
SIX_FF_LF                 (idx, [1:   2]) = [  9.99704E-01 1.7E-06 ];
SIX_FF_LT                 (idx, [1:   2]) = [  9.33627E-01 2.8E-05 ];
SIX_FF_KINF               (idx, [1:   2]) = [  1.11384E+00 0.00012 ];
SIX_FF_KEFF               (idx, [1:   2]) = [  1.03961E+00 0.00012 ];

% Fission neutron and energy production:

NUBAR                     (idx, [1:   2]) = [  2.51698E+00 7.1E-06 ];
FISSE                     (idx, [1:   2]) = [  2.02989E+02 5.9E-07 ];

% Criticality eigenvalues:

ANA_KEFF                  (idx, [1:   6]) = [  1.03963E+00 0.00012  1.03224E+00 0.00012  7.36933E-03 0.00176 ];
IMP_KEFF                  (idx, [1:   2]) = [  1.03964E+00 7.0E-05 ];
COL_KEFF                  (idx, [1:   2]) = [  1.03952E+00 9.6E-05 ];
ABS_KEFF                  (idx, [1:   2]) = [  1.03964E+00 7.0E-05 ];
ABS_KINF                  (idx, [1:   2]) = [  1.11389E+00 6.8E-05 ];
GEOM_ALBEDO               (idx, [1:   6]) = [  1.00000E+00 0.0E+00  1.00000E+00 0.0E+00  1.00000E+00 0.0E+00 ];

% ALF (Average lethargy of neutrons causing fission):
% Based on E0 = 2.000000E+01 MeV

ANA_ALF                   (idx, [1:   2]) = [  1.09306E+01 0.00012 ];
IMP_ALF                   (idx, [1:   2]) = [  1.09305E+01 9.0E-05 ];

% EALF (Energy corresponding to average lethargy of neutrons causing fission):

ANA_EALF                  (idx, [1:   2]) = [  3.58329E-04 0.00128 ];
IMP_EALF                  (idx, [1:   2]) = [  3.58243E-04 0.00099 ];

% AFGE (Average energy of neutrons causing fission):

ANA_AFGE                  (idx, [1:   2]) = [  7.07842E-01 0.00032 ];
IMP_AFGE                  (idx, [1:   2]) = [  7.07711E-01 0.00019 ];

% Forward-weighted delayed neutron parameters:

FWD_ANA_BETA_ZERO         (idx, [1:  18]) = [  8.44270E-03 0.00111  1.98916E-04 0.00683  1.11990E-03 0.00289  6.10968E-04 0.00405  1.46853E-03 0.00263  2.61659E-03 0.00190  1.13107E-03 0.00288  8.58977E-04 0.00336  4.37737E-04 0.00464 ];
FWD_ANA_LAMBDA            (idx, [1:  18]) = [  5.60803E-01 0.00158  1.24667E-02 7.2E-09  2.82917E-02 0.0E+00  4.25244E-02 0.0E+00  1.33042E-01 0.0E+00  2.92467E-01 0.0E+00  6.66488E-01 2.7E-09  1.63478E+00 4.2E-09  3.55460E+00 2.9E-09 ];

% Beta-eff using Meulekamp's method:

ADJ_MEULEKAMP_BETA_EFF    (idx, [1:  18]) = [  7.31477E-03 0.00165  1.68574E-04 0.01045  9.94792E-04 0.00453  5.19082E-04 0.00639  1.27531E-03 0.00398  2.27709E-03 0.00287  9.68238E-04 0.00443  7.35500E-04 0.00515  3.76186E-04 0.00722 ];
ADJ_MEULEKAMP_LAMBDA      (idx, [1:  18]) = [  5.56784E-01 0.00250  1.24667E-02 7.2E-09  2.82917E-02 0.0E+00  4.25244E-02 0.0E+00  1.33042E-01 0.0E+00  2.92467E-01 0.0E+00  6.66488E-01 2.7E-09  1.63478E+00 4.2E-09  3.55460E+00 3.0E-09 ];

% Adjoint weighted time constants using Nauchi's method:

ADJ_NAUCHI_GEN_TIME       (idx, [1:   6]) = [  4.62729E-04 0.00052  4.62264E-04 0.00052  5.27231E-04 0.00563 ];
ADJ_NAUCHI_LIFETIME       (idx, [1:   6]) = [  4.81063E-04 0.00052  4.80580E-04 0.00052  5.48124E-04 0.00562 ];
ADJ_NAUCHI_BETA_EFF       (idx, [1:  18]) = [  7.09057E-03 0.00181  1.64489E-04 0.01236  9.64144E-04 0.00501  5.02669E-04 0.00673  1.23407E-03 0.00463  2.21230E-03 0.00313  9.38182E-04 0.00526  7.11330E-04 0.00590  3.63388E-04 0.00820 ];
ADJ_NAUCHI_LAMBDA         (idx, [1:  18]) = [  5.56040E-01 0.00288  1.24667E-02 7.2E-09  2.82917E-02 0.0E+00  4.25244E-02 0.0E+00  1.33042E-01 0.0E+00  2.92467E-01 0.0E+00  6.66488E-01 2.7E-09  1.63478E+00 4.2E-09  3.55460E+00 2.8E-09 ];

% Adjoint weighted time constants using IFP:

ADJ_IFP_GEN_TIME          (idx, [1:   6]) = [  4.84149E-04 0.00129  4.83653E-04 0.00129  5.54303E-04 0.01423 ];
ADJ_IFP_LIFETIME          (idx, [1:   6]) = [  5.03333E-04 0.00129  5.02818E-04 0.00129  5.76273E-04 0.01423 ];
ADJ_IFP_IMP_BETA_EFF      (idx, [1:  18]) = [  7.13133E-03 0.00594  1.59640E-04 0.04160  9.77148E-04 0.01546  5.00680E-04 0.02235  1.25540E-03 0.01443  2.22579E-03 0.01085  9.38960E-04 0.01640  7.08185E-04 0.01823  3.65528E-04 0.02588 ];
ADJ_IFP_IMP_LAMBDA        (idx, [1:  18]) = [  5.54951E-01 0.00868  1.24667E-02 0.0E+00  2.82917E-02 0.0E+00  4.25244E-02 0.0E+00  1.33042E-01 0.0E+00  2.92467E-01 0.0E+00  6.66488E-01 2.7E-09  1.63478E+00 4.3E-09  3.55460E+00 3.0E-09 ];
ADJ_IFP_ANA_BETA_EFF      (idx, [1:  18]) = [  7.13522E-03 0.00582  1.60426E-04 0.04056  9.78173E-04 0.01518  4.96763E-04 0.02204  1.25811E-03 0.01421  2.22878E-03 0.01063  9.40140E-04 0.01603  7.08907E-04 0.01798  3.63928E-04 0.02525 ];
ADJ_IFP_ANA_LAMBDA        (idx, [1:  18]) = [  5.54217E-01 0.00850  1.24667E-02 0.0E+00  2.82917E-02 0.0E+00  4.25244E-02 0.0E+00  1.33042E-01 0.0E+00  2.92467E-01 0.0E+00  6.66488E-01 2.7E-09  1.63478E+00 4.3E-09  3.55460E+00 2.8E-09 ];
ADJ_IFP_ROSSI_ALPHA       (idx, [1:   2]) = [ -1.47712E+01 0.00610 ];

% Adjoint weighted time constants using perturbation technique:

ADJ_PERT_GEN_TIME         (idx, [1:   2]) = [  4.76062E-04 0.00029 ];
ADJ_PERT_LIFETIME         (idx, [1:   2]) = [  4.94923E-04 0.00026 ];
ADJ_PERT_BETA_EFF         (idx, [1:   2]) = [  7.17129E-03 0.00113 ];
ADJ_PERT_ROSSI_ALPHA      (idx, [1:   2]) = [ -1.50651E+01 0.00117 ];

% Inverse neutron speed :

ANA_INV_SPD               (idx, [1:   2]) = [  1.93990E-06 0.00016 ];

% Analog slowing-down and thermal neutron lifetime (total/prompt/delayed):

ANA_SLOW_TIME             (idx, [1:   6]) = [  1.07069E-05 0.00011  1.07067E-05 0.00011  1.07283E-05 0.00122 ];
ANA_THERM_TIME            (idx, [1:   6]) = [  1.63246E-03 0.00024  1.63278E-03 0.00024  1.59449E-03 0.00267 ];
ANA_THERM_FRAC            (idx, [1:   6]) = [  2.72744E-01 0.00018  2.72793E-01 0.00018  2.67304E-01 0.00243 ];
ANA_DELAYED_EMTIME        (idx, [1:   2]) = [  1.09133E+01 0.00239 ];
ANA_MEAN_NCOL             (idx, [1:   4]) = [  8.92797E+01 0.00019  8.82382E+01 0.00029 ];

% Group constant generation:

GC_UNIVERSE_NAME          (idx, [1:  3])  = 'D2O' ;

% Micro- and macro-group structures:

MICRO_NG                  (idx, 1)        = 70 ;
MICRO_E                   (idx, [1:  71]) = [  1.00000E-11  5.00000E-09  1.00000E-08  1.50000E-08  2.00000E-08  2.50000E-08  3.00000E-08  3.50000E-08  4.20000E-08  5.00000E-08  5.80000E-08  6.70000E-08  8.00000E-08  1.00000E-07  1.40000E-07  1.80000E-07  2.20000E-07  2.50000E-07  2.80000E-07  3.00000E-07  3.20000E-07  3.50000E-07  4.00000E-07  5.00000E-07  6.25000E-07  7.80000E-07  8.50000E-07  9.10000E-07  9.50000E-07  9.72000E-07  9.96000E-07  1.02000E-06  1.04500E-06  1.07100E-06  1.09700E-06  1.12300E-06  1.15000E-06  1.30000E-06  1.50000E-06  1.85500E-06  2.10000E-06  2.60000E-06  3.30000E-06  4.00000E-06  9.87700E-06  1.59680E-05  2.77000E-05  4.80520E-05  7.55014E-05  1.48728E-04  3.67262E-04  9.06898E-04  1.42510E-03  2.23945E-03  3.51910E-03  5.50000E-03  9.11800E-03  1.50300E-02  2.47800E-02  4.08500E-02  6.74300E-02  1.11000E-01  1.83000E-01  3.02500E-01  5.00000E-01  8.21000E-01  1.35300E+00  2.23100E+00  3.67900E+00  6.06550E+00  2.00000E+01 ];

MACRO_NG                  (idx, 1)        = 2 ;
MACRO_E                   (idx, [1:   3]) = [  1.00000E+37  6.25000E-07  0.00000E+00 ];

% Micro-group spectrum:

INF_MICRO_FLX             (idx, [1: 140]) = [  1.44415E+05 0.00159  5.82431E+05 0.00071  1.09231E+06 0.00037  1.17272E+06 0.00037  1.18922E+06 0.00044  1.78398E+06 0.00032  1.85556E+06 0.00028  2.31280E+06 0.00029  2.35627E+06 0.00029  2.40547E+06 0.00027  2.43899E+06 0.00026  2.41463E+06 0.00028  2.35679E+06 0.00029  2.29740E+06 0.00028  2.26463E+06 0.00032  1.95259E+06 0.00032  1.93189E+06 0.00033  1.89140E+06 0.00027  1.84939E+06 0.00034  3.57515E+06 0.00029  3.42905E+06 0.00032  2.47372E+06 0.00032  1.60325E+06 0.00031  1.91027E+06 0.00030  1.85987E+06 0.00030  1.59067E+06 0.00035  2.91965E+06 0.00028  6.19113E+05 0.00032  7.72746E+05 0.00039  6.98332E+05 0.00043  4.08167E+05 0.00049  7.07235E+05 0.00036  4.81103E+05 0.00046  4.14334E+05 0.00048  8.02741E+04 0.00081  7.91489E+04 0.00056  8.14510E+04 0.00105  8.36953E+04 0.00092  8.25878E+04 0.00090  8.15220E+04 0.00100  8.36136E+04 0.00067  7.89131E+04 0.00098  1.48390E+05 0.00068  2.37381E+05 0.00057  3.02068E+05 0.00054  7.93463E+05 0.00042  8.25676E+05 0.00038  8.65584E+05 0.00050  5.44532E+05 0.00060  3.81965E+05 0.00062  2.86004E+05 0.00049  3.20773E+05 0.00046  5.84555E+05 0.00043  8.25287E+05 0.00047  2.09434E+06 0.00042  5.58121E+06 0.00041  1.69080E+07 0.00041  1.85435E+07 0.00040  1.85378E+07 0.00041  1.63326E+07 0.00040  1.75693E+07 0.00039  2.03608E+07 0.00038  1.99491E+07 0.00038  1.52144E+07 0.00040  1.55703E+07 0.00040  1.55024E+07 0.00039  1.46833E+07 0.00040  1.26283E+07 0.00039  9.30193E+06 0.00040  3.68585E+06 0.00042 ];

% Integral parameters:

INF_KINF                  (idx, [1:   2]) = [  0.00000E+00 0.0E+00 ];

% Flux spectra in infinite geometry:

INF_FLX                   (idx, [1:   4]) = [  2.98067E+01 0.00022  1.13028E+02 0.00037 ];
INF_FISS_FLX              (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Reaction cross sections:

INF_TOT                   (idx, [1:   4]) = [  3.43412E-01 6.8E-06  4.90167E-01 6.6E-06 ];
INF_CAPT                  (idx, [1:   4]) = [  2.82226E-05 0.00081  3.45933E-05 2.2E-05 ];
INF_ABS                   (idx, [1:   4]) = [  2.82226E-05 0.00081  3.45933E-05 2.2E-05 ];
INF_FISS                  (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_NSF                   (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_NUBAR                 (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_KAPPA                 (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_INVV                  (idx, [1:   4]) = [  1.04269E-07 9.4E-05  3.90984E-06 2.2E-05 ];

% Total scattering cross sections:

INF_SCATT0                (idx, [1:   4]) = [  3.43383E-01 6.9E-06  4.90132E-01 6.6E-06 ];
INF_SCATT1                (idx, [1:   4]) = [  7.27902E-02 6.9E-05  5.36373E-02 6.2E-05 ];
INF_SCATT2                (idx, [1:   4]) = [  1.13905E-02 0.00044  5.71271E-04 0.00479 ];
INF_SCATT3                (idx, [1:   4]) = [  6.86652E-04 0.00512 -5.23732E-03 0.00047 ];
INF_SCATT4                (idx, [1:   4]) = [  4.63264E-05 0.07301 -1.07636E-02 0.00018 ];
INF_SCATT5                (idx, [1:   4]) = [  7.45568E-05 0.04378 -5.58359E-03 0.00034 ];
INF_SCATT6                (idx, [1:   4]) = [ -9.26642E-05 0.03371 -1.13582E-02 0.00019 ];
INF_SCATT7                (idx, [1:   4]) = [  3.19688E-06 1.00000 -2.85768E-03 0.00061 ];

% Total scattering production cross sections:

INF_SCATTP0               (idx, [1:   4]) = [  3.43416E-01 6.8E-06  4.90132E-01 6.6E-06 ];
INF_SCATTP1               (idx, [1:   4]) = [  7.28141E-02 6.9E-05  5.36373E-02 6.2E-05 ];
INF_SCATTP2               (idx, [1:   4]) = [  1.14045E-02 0.00044  5.71271E-04 0.00479 ];
INF_SCATTP3               (idx, [1:   4]) = [  6.93316E-04 0.00507 -5.23732E-03 0.00047 ];
INF_SCATTP4               (idx, [1:   4]) = [  4.90357E-05 0.06875 -1.07636E-02 0.00018 ];
INF_SCATTP5               (idx, [1:   4]) = [  7.55660E-05 0.04315 -5.58359E-03 0.00034 ];
INF_SCATTP6               (idx, [1:   4]) = [ -9.22612E-05 0.03383 -1.13582E-02 0.00019 ];
INF_SCATTP7               (idx, [1:   4]) = [  3.38834E-06 1.00000 -2.85768E-03 0.00061 ];

% Diffusion parameters:

INF_TRANSPXS              (idx, [1:   4]) = [  2.64203E-01 3.2E-05  4.21275E-01 1.2E-05 ];
INF_DIFFCOEF              (idx, [1:   4]) = [  1.26166E+00 3.2E-05  7.91250E-01 1.2E-05 ];

% Reduced absoption and removal:

INF_RABSXS                (idx, [1:   4]) = [ -4.19792E-06 0.02550  3.45933E-05 2.2E-05 ];
INF_REMXS                 (idx, [1:   4]) = [  9.21335E-03 8.5E-05  5.81574E-05 0.00110 ];

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

INF_S0                    (idx, [1:   8]) = [  3.34198E-01 6.9E-06  9.18512E-03 8.6E-05  2.34214E-05 0.00196  4.90109E-01 6.6E-06 ];
INF_S1                    (idx, [1:   8]) = [  7.35619E-02 6.7E-05 -7.71754E-04 0.00125  9.61895E-06 0.00248  5.36277E-02 6.2E-05 ];
INF_S2                    (idx, [1:   8]) = [  1.22468E-02 0.00041 -8.56253E-04 0.00087  3.20263E-06 0.00691  5.68068E-04 0.00481 ];
INF_S3                    (idx, [1:   8]) = [  9.19183E-04 0.00392 -2.32530E-04 0.00266  4.94820E-07 0.03919 -5.23781E-03 0.00047 ];
INF_S4                    (idx, [1:   8]) = [  1.48020E-04 0.02208 -1.01694E-04 0.00553 -5.15802E-07 0.03115 -1.07631E-02 0.00018 ];
INF_S5                    (idx, [1:   8]) = [  7.56592E-05 0.04163 -1.10236E-06 0.39425 -8.25484E-07 0.01482 -5.58276E-03 0.00034 ];
INF_S6                    (idx, [1:   8]) = [ -1.31443E-05 0.23868 -7.95199E-05 0.00591 -7.95670E-07 0.01599 -1.13574E-02 0.00019 ];
INF_S7                    (idx, [1:   8]) = [ -3.08043E-05 0.11299  3.40012E-05 0.01288 -5.89112E-07 0.02154 -2.85710E-03 0.00061 ];

% Scattering production matrixes:

INF_SP0                   (idx, [1:   8]) = [  3.34231E-01 6.8E-06  9.18512E-03 8.6E-05  2.34214E-05 0.00196  4.90109E-01 6.6E-06 ];
INF_SP1                   (idx, [1:   8]) = [  7.35859E-02 6.7E-05 -7.71754E-04 0.00125  9.61895E-06 0.00248  5.36277E-02 6.2E-05 ];
INF_SP2                   (idx, [1:   8]) = [  1.22608E-02 0.00041 -8.56253E-04 0.00087  3.20263E-06 0.00691  5.68068E-04 0.00481 ];
INF_SP3                   (idx, [1:   8]) = [  9.25846E-04 0.00389 -2.32530E-04 0.00266  4.94820E-07 0.03919 -5.23781E-03 0.00047 ];
INF_SP4                   (idx, [1:   8]) = [  1.50730E-04 0.02159 -1.01694E-04 0.00553 -5.15802E-07 0.03115 -1.07631E-02 0.00018 ];
INF_SP5                   (idx, [1:   8]) = [  7.66684E-05 0.04102 -1.10236E-06 0.39425 -8.25484E-07 0.01482 -5.58276E-03 0.00034 ];
INF_SP6                   (idx, [1:   8]) = [ -1.27413E-05 0.24610 -7.95199E-05 0.00591 -7.95670E-07 0.01599 -1.13574E-02 0.00019 ];
INF_SP7                   (idx, [1:   8]) = [ -3.06128E-05 0.11373  3.40012E-05 0.01288 -5.89112E-07 0.02154 -2.85710E-03 0.00061 ];

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

CMM_TRANSPXS              (idx, [1:   4]) = [  2.02586E-01 0.00018 -7.83338E-01 0.00028 ];
CMM_TRANSPXS_X            (idx, [1:   4]) = [  2.05516E-01 0.00026 -7.89027E-01 0.00032 ];
CMM_TRANSPXS_Y            (idx, [1:   4]) = [  2.01154E-01 0.00021 -7.80533E-01 0.00030 ];
CMM_TRANSPXS_Z            (idx, [1:   4]) = [  2.01150E-01 0.00023 -7.80517E-01 0.00032 ];
CMM_DIFFCOEF              (idx, [1:   4]) = [  1.64540E+00 0.00018 -4.25531E-01 0.00028 ];
CMM_DIFFCOEF_X            (idx, [1:   4]) = [  1.62194E+00 0.00026 -4.22464E-01 0.00032 ];
CMM_DIFFCOEF_Y            (idx, [1:   4]) = [  1.65711E+00 0.00021 -4.27060E-01 0.00030 ];
CMM_DIFFCOEF_Z            (idx, [1:   4]) = [  1.65715E+00 0.00023 -4.27070E-01 0.00032 ];

% Delayed neutron parameters (Meulekamp method):

BETA_EFF                  (idx, [1:  18]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
LAMBDA                    (idx, [1:  18]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
