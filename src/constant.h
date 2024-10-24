
// maximum length of a protein sequence
#define MAX_SEQ_LEN 65536

// maximum length of comments
#define MAX_COM_LEN 65536

// maximum number of index lines per DPU
#define MAX_NB_INDEX_LINE_DPU 4194304

// max number of db sequences per DPU (16 K)
#define MAX_NB_SEQ_DPU  16384

// max size of the db per dpu (8 MB)
#define MAX_SIZE_DB_DPU 8388608

// size of offset index 
#define SIZE_OFFSET_INDEX 65536

// maximum length of a query sequence on DPU
#define MAX_QUERY_LEN_DPU 1024

// max distance to merge 2 hits
#define MAX_HIT_DIST 20

// maximum number of hits per DPU
#define MAX_HIT_DPU 65536

// maximum size of a display comment
#define MAX_SIZE_DISPLAY_COM 30

// default e-value
#define DEFAULT_EVALUE 1e-6

// xdrop value
#define XDROP 16

// size of the seed
#define SEED_SIZE 4

// seed mask
#define SEED_MASK 0xFFFF

// seed neighborhood size
#define SIZE_NEIGHBOR 8

// max number of ranks used on UPMEM device
#define MAX_NB_RANKS 40

// number of DPU per rank
#define NB_DPU_PER_RANK 64

// =======================================
// Karling Altschul statistical parameters
// =======================================
//
// https://www.ncbi.nlm.nih.gov/CBBresearch/Spouge/html_ncbi/html/blast/index.html
//
// WARNING : only for nogap alignment

// K et LAMBDA
#define BLOSUM45_LAMBDA 0.229091
#define BLOSUM45_K      0.092395
#define BLOSUM50_LAMBDA 0.231827
#define BLOSUM50_K      0.111519
#define BLOSUM62_LAMBDA 0.317606
#define BLOSUM62_K      0.133741
#define BLOSUM80_LAMBDA 0.342969
#define BLOSUM80_K      0.176899
#define BLOSUM90_LAMBDA 0.334565
#define BLOSUM90_K      0.189948
#define PAM30_LAMBDA    0.340025
#define PAM30_K         0.282799
#define PAM70_LAMBDA    0.334503
#define PAM70_K         0.228671
#define PAM250_LAMBDA   0.225153
#define PAM250_K        0.086737

// default subsitution martrix
#define DEFAULT_SUB_MATRIX "blosum62"
#define DEFAULT_LAMBDA BLOSUM62_LAMBDA
#define DEFAULT_K BLOSUM62_K
