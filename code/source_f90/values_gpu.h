// values for the repartitions on the GPU
// it seems that the 1 can not be changed ..why ?
// GRIDX has to be odd..why ?
#define BLOCK1 384
#define BLOCK2 1
#define BLOCK3 1
#define GRIDX 257
#define GRIDY 1
#define GRIDZ 1
#define BLOCK1B 4
#define WARPSIZE 64
#define BATCH 256 //The number of batched ffts 


