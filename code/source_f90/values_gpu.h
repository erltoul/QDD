// values for the repartitions on the GPU
// it seems that the 1 can not be changed ..why ?
// GRIDX has to be odd..why ?
#define BLOCK1 384
// for older GPUs 192 is better but performance is reduced
//#define BLOCK1 192
#define BLOCK2 1
#define BLOCK3 1
 #define GRIDX 257
// for older GPUs 65 is better but performance is reduced
//#define GRIDX 65
#define GRIDY 1
#define GRIDZ 1
#define BLOCK1B 4
#define WARPSIZE 32
#define BATCH 256 //The number of batched ffts 
// for older GPUs 64 is better but performance is reduced
//#define BATCH 64 //The number of batched ffts 


