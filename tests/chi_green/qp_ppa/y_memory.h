 implicit none
#define YAMBO_ALLOC(x,s) allocate(x s)
#define YAMBO_FREE(x) if(allocated(x)) deallocate(x)
#define YAMBO_ALLOC_GPU(x,s) allocate(x s)
#define YAMBO_ALLOC_GPU_MOLD(x,y) allocate(x, mold=y)
#define YAMBO_FREE_GPU(x) if(allocated(x)) deallocate(x)
