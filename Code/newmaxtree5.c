#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <math.h>
#include <pthread.h>
#include <sys/types.h>
#include <sys/times.h>
#include <unistd.h>
#include <assert.h>
#include <FreeImage.h>


/*
Floating-point Max-tree algorithm

Michael Wilkinson

*/
#define BOTTOM (-1)

#define PI 3.14159265358979323846
#define false 0
#define true  1
#define MAXTHREADS 128
int MULFACTOR = 1;     

#define CONNECTIVITY  6
#define MIN(a,b)  ((a<=b) ? (a) : (b))
#define MAX(a,b)  ((a>=b) ? (a) : (b))
#define LWB(self) (size2D*(((self)*depth)/nthreads))
#define UPB(self) (size2D*(((self+1)*depth)/nthreads))
typedef short bool;
typedef unsigned char ubyte;

int ParCount=0;
int nthreads;
pthread_t threadID[MAXTHREADS];

int width, height, depth, size;  /* precondition: width <= size/nthreads */
int size2D;
double lambda;
clock_t start,built, merged,finish;
struct tms tstruct;


typedef unsigned char Pixel;

Pixel *gval=NULL, *out=NULL;

typedef struct MaxNode  
{ 
  int parent;
  int Area;
  bool filtered; /* indicates whether or not the filtered value is OK */
  Pixel gval;
  Pixel outval;
} MaxNode;

#define bottom (-1)


MaxNode *node; 

/************************** safe malloc and calloc ************************************/

pthread_mutex_t mutex = PTHREAD_MUTEX_INITIALIZER;

void *SafeMalloc(int n) {
  void *ptr;
  pthread_mutex_lock(&mutex);
  ptr = malloc(n);
  if (ptr==NULL) {
    fprintf (stderr, "Error: out of memory. "
	     "Could not allocate %d bytes\n", n);
  }
  pthread_mutex_unlock(&mutex);
  return ptr;
}

void *SafeCalloc(int nmemb, int size) {
  void *ptr;
  pthread_mutex_lock(&mutex);
  ptr = calloc(nmemb, size);
  if (ptr==NULL) {
    fprintf (stderr, "Error: out of memory. Could not "
	     "allocate %d bytes\n", nmemb*size);
  }
  pthread_mutex_unlock(&mutex);
  return ptr;
}

/************************** semaphore ************************************/

pthread_mutex_t samut[MAXTHREADS];
pthread_cond_t  sacv[MAXTHREADS];
int             saval[MAXTHREADS];

void Psa(int p) {
  pthread_mutex_lock(&samut[p]);
  while (saval[p] <= 0)
    pthread_cond_wait(&sacv[p], &samut[p]);
  saval[p]--;
  pthread_mutex_unlock(&samut[p]);
}

void Vsa(int p) {
  pthread_mutex_lock(&samut[p]);
  saval[p]++;
  pthread_mutex_unlock(&samut[p]);
  pthread_cond_broadcast(&sacv[p]);
}

/************************** barrier ************************************/

pthread_mutex_t barriermutex = PTHREAD_MUTEX_INITIALIZER;
pthread_cond_t barriercv = PTHREAD_COND_INITIALIZER;
int barcnt = 0;

void Barrier(int self) {
  pthread_mutex_lock(&barriermutex);
  barcnt++;
  if (barcnt == nthreads) {
    barcnt = 0;  /* for reuse of routine */
    pthread_cond_broadcast (&barriercv);  
  } else {
    pthread_cond_wait (&barriercv, &barriermutex);
  }
  pthread_mutex_unlock(&barriermutex);
}

/************************** level root ************************************/

ulong levroot(ulong x) {
  int r=x, y ,gv=gval[x];
  if (r==bottom) return bottom;
  while ((node[r].parent!=bottom) && (gv==gval[node[r].parent]))
      r = node[r].parent;

  while (x!=r){
    y=node[x].parent;
    node[x].parent=r;
    x=y;
  }
  return r;
}

ulong Par(ulong x) {
  ParCount++;
  return levroot(node[x].parent);
}

void levrootfix(ulong lwb, ulong upb) {
  ulong z, u, x;

  for (x=lwb; x<upb;x++) {
    u = levroot(x);
    if (x!=u) node[x].parent=u;
    else {
      z = Par(x); 
      node[x].parent=z;
    }
  }
}

typedef struct{
  int  curpos, maxsize;
  int *stack;  
} pStack;


pStack *pStackCreate(long maxsize, int *stackArray){
  pStack *newStack = (pStack *) malloc(sizeof(pStack));
  
  newStack->stack = stackArray;
  newStack->maxsize= maxsize;
  newStack->curpos = newStack->maxsize;

  return newStack;
}

void pStackDelete(pStack *oldstack){
  free(oldstack);
}

#define IsEmptyStack(stack)        ((stack->curpos)==stack->maxsize)
#define IsFullStack(stack)         ( stack->curpos == 0 )
#define IsEmpty(queue)             ((queue->size)==0)
#define IsFull(queue)              ( queue->size == queue->maxsize )
#define pStackTop(stack)       ( stack->stack[stack->curpos + 1])
#define pStackPop(stack)       (stack->stack[++stack->curpos])
#define pStackPush(stack,elem)      (stack->stack[stack->curpos--]=elem)

typedef struct{
  int size,maxsize;
  int *queue;
} pQueue;

pQueue *pQueueCreate(long maxsize){
  pQueue *newQueue = (pQueue *) malloc(sizeof(pQueue));
  newQueue->size = 0;
  newQueue->queue = (int *)malloc((maxsize+1)*sizeof(int));
  newQueue->maxsize=maxsize;
  return newQueue;
}


#define pQueueFront(queue)       (queue->queue[1])

void pQueueDelete(pQueue *oldqueue){
  free(oldqueue->queue);
  free(oldqueue);

}

int pQueuePop(pQueue *queue, Pixel *priority){
  int outval = queue->queue[1];
  int current = 1,
    moved;
  Pixel curval;

  moved = queue->queue[queue->size];
  queue->size--;
  curval = priority[moved];

  while ( ((current*2<=queue->size) &&
	   (curval< priority[queue->queue[current*2]]))
	  ||
	  ((current*2+1<=queue->size) &&
	   (curval< priority[queue->queue[current*2+1]]))
	  ){
    if ((current*2+1<=queue->size) && 
	(priority[queue->queue[current*2]]< 
	 priority[queue->queue[current*2+1]])){
      queue->queue[current]= queue->queue[current*2+1];
      current+=current+1;
    } else {
      queue->queue[current]= queue->queue[current*2];
      current+=current;
    }
  }
  queue->queue[current]=moved;
  return outval;
}

void pQueuePush(pQueue *queue, Pixel *priority, int pixpos){
  long current;
  Pixel curval = priority[pixpos];
  queue->size++;
  current=queue->size;
  
  while ((current/2 !=0) && (priority[queue->queue[current/2]]<curval)){
    queue->queue[current]= queue->queue[current/2];
    current=current/2;
  }

  queue->queue[current]=pixpos;
}

int GetNeighbors(int p, int x, int y, int z, 
		 int *neighbors, int width, int height, int depth,
		 int size2D, int lwb, int upb)
{

   int n=0;


   if (x<(width-1))       neighbors[n++] = p+1;
   if (y>0)               neighbors[n++] = p-width;
   if (x>0)               neighbors[n++] = p-1;
   if (y<(height-1))      neighbors[n++] = p+width;
   if (depth>1) {
      if (z>lwb)          neighbors[n++] = p-size2D;
      if (z<(upb-1))           neighbors[n++] = p+size2D;
   }
   return(n);
} /* GetNeighbors */


void Flood(int self, pQueue *queue, pStack *stack, 
	   Pixel *gval, int width, int height, int depth,
	   MaxNode *node){
  int neighbors[CONNECTIVITY];
  int size = width*height*depth;
  int size2D = width*height;
  int lwb, upb;

  int xm,i, x,y,z, nextpix, p, q, numneighbors, oldtop, oldpix;

  lwb = LWB(self); upb = UPB(self);
  xm = lwb;
  for (x=xm; x<upb; x++) { 
    node[x].parent = bottom;
    node[x].filtered = false;
    if (gval[xm]>gval[x]) xm = x;
  }



  node[xm].parent = xm;
  pStackPush(stack,xm);
  node[xm].Area = 1;
   
  nextpix = xm;
  /*  printf("stack max size = %d, queue max size =%d\n", stack->maxsize, queue->maxsize);
   */
  do{ 
    p = nextpix;
    assert( ( p == pStackTop(stack) ) || ( p == pQueueFront(queue) ) );
    oldpix = p;
    x = p % width; 
    y = (p % size2D) / width;
    z = p / size2D;

 
    

    numneighbors = GetNeighbors(p, x, y, z, neighbors, width, height,
				depth, size2D, lwb/size2D, upb/size2D);    

    for (i=0; i<numneighbors; i++) {
      q = neighbors[i];
      if (node[q].parent == BOTTOM){
	node[q].parent=q;
	node[q].Area=1;
        if (gval[q]>gval[p]){
	  pStackPush(stack,q);
	  nextpix = q;

          assert (nextpix == pStackTop(stack) );
          assert (nextpix != oldpix );
	  break;
	}
	pQueuePush(queue,gval,q);
      }
    }

    if (nextpix==p){      /* No break occurred */

      if (p != pStackTop(stack)){      /* p is in queue and 
				          processing is finished  */

	assert(gval[p] == gval[pStackTop(stack)]);

	p =  pQueuePop(queue,gval);     /* remove from queue */

	node[p].parent = pStackTop(stack);
	node[pStackTop(stack)].Area++;

	if (!IsEmpty(queue)){

	  nextpix = pQueueFront(queue);

	  if (gval[nextpix] < gval[p]){   
	    /* moving down, but first process top of stack */
	    
	    nextpix = pStackTop(stack);
	  }
	} else {
	    nextpix = pStackTop(stack);
	}
   
	assert( ( nextpix != oldpix ) &&
	        (( nextpix == pStackTop(stack) ) || 
		 ( nextpix == pQueueFront(queue) ) ));
        
	
      } else {                          /* p is top of stack */

	if (!IsEmpty(queue)){
	  
	  nextpix = pQueueFront(queue);    /* candidate next pixel */

	  assert( oldpix != nextpix );

	  if (gval[nextpix] < gval[p]){ /* moving down, top of stack done */

	    oldtop = pStackPop(stack);

	    p = pStackTop(stack);

	    assert( oldpix != p );
 
	    if (gval[nextpix] > gval[p]){

	      nextpix = pQueuePop(queue, gval); /* move nextpix from 
					     queue to stack     */
 
	      pStackPush(stack, nextpix);

              p = nextpix; 

	      assert( oldpix != nextpix );
        
	    } else if (gval[nextpix] < gval[p]){

              nextpix = p;              /* flood from current stack top */

	    }

	    node[oldtop].parent = p;           /* == stack top */
	    node[p].Area += node[oldtop].Area; /* update area */
           
	  }
	  assert( oldpix != nextpix );
	} else {

	    oldtop = pStackPop(stack);
            if (!IsEmptyStack(stack)){
	      p = pStackTop(stack);
	      node[oldtop].parent = p;           /* == stack top */
	      node[p].Area += node[oldtop].Area; /* update area */
	      nextpix = p;
	      assert( oldpix != nextpix );
	    }
	}
	
      }

    }

  } while (!IsEmpty(queue) || (!IsEmptyStack(stack)));

    
  node[xm].parent=bottom;
}


void MaxTreeAreaFilter(int self, double lambda)
{
   ulong v, u, w, parent, lwb, upb;
   Pixel val;
 
   lwb = LWB(self); upb = UPB(self);
   for (v=lwb; v<upb; v++) {
      if (!node[v].filtered) {
         w = v;
         parent = node[w].parent;
         while ((parent != bottom) && (!node[w].filtered) && 
            ((gval[w] == gval[parent]) || (node[w].Area < lambda))) {
                  w = parent;
                  parent = node[w].parent;
         }
         if (node[w].filtered) val = node[w].outval;
         else if (node[w].Area >= lambda) val = gval[w]; 
         else val = 0; /* criterion cannot be satisfied */
         u = v;
         while (u!=w) {
           if ((lwb<=u) && (u<upb)){ 
	     node[u].outval = val;
	     node[u].filtered=true;
	   }

           u = node[u].parent;
           }
         if ((lwb<=w) && (w<upb)){ 
	   node[w].outval = val;
	   node[w].filtered=true;
	 }
       }
      out[v]=node[v].outval;
   }
} /* MaxTreeAreaFilter */


/** Generic image loader
@param lpszPathName Pointer to the full file name
@param flag Optional load flag constant
@return Returns the loaded dib if successful, returns NULL otherwise
*/
FIBITMAP* GenericLoader(const char* lpszPathName, int flag) {
  FREE_IMAGE_FORMAT fif = FIF_UNKNOWN;
//  check the file signature and deduce its format
// (the second argument is currently not used by FreeImage)
  fif = FreeImage_GetFileType(lpszPathName, 0);
  if(fif == FIF_UNKNOWN) {
    // no signature ?
    // try to guess the file format from the file extension
    fif = FreeImage_GetFIFFromFilename(lpszPathName);
  }
  // check that the plugin has reading capabilities ...
  if((fif != FIF_UNKNOWN) && FreeImage_FIFSupportsReading(fif)) {
    // ok, let's load the file
    FIBITMAP *dib = FreeImage_Load(fif, lpszPathName, flag);
    // unless a bad file format, we are done !
    return dib;
  }
  return NULL;
}

Pixel *ReadTIFF(char *fnm, int *width, int *height){
   FIBITMAP *dib = GenericLoader(fnm,0);
   unsigned long  bitsperpixel;
   Pixel *im;
   unsigned int x,y,i,imsize;
   if (dib == NULL) return NULL;
     
   bitsperpixel =  FreeImage_GetBPP(dib);
   *height = FreeImage_GetHeight(dib), 
   *width = FreeImage_GetWidth(dib);
   imsize = (*width)*(*height);
   im = calloc((size_t)imsize, sizeof(Pixel));
   assert(im!=NULL);
       switch(bitsperpixel) {
         case 8:
           i=0;
           for(y = 0; y < *height; y++) {
               BYTE *bits = (BYTE *)FreeImage_GetScanLine(dib, y);
               for(x = 0; x < *width; x++,i++) {
                 im[i] = bits[x];
               }
             }
           
           FreeImage_Unload(dib);
           return im;
         case 16:
	   i=0;
           for(y = 0; y < *height; y++) {
             unsigned short *bits = (unsigned short *)FreeImage_GetScanLine(dib, y);
               for(x = 0; x < *width; x++,i++) {
                 im[i] = bits[x];
               }
           }
           FreeImage_Unload(dib);
           return im;  
         default : 
           FreeImage_Unload(dib);

           return NULL; 
       }
 
}

void WriteTIFF( char *fname, Pixel *im, long width, long height, int bitspp){
  FIBITMAP *outmap;
  long i,j,y,x; 
  FREE_IMAGE_FORMAT fif = FreeImage_GetFIFFromFilename(fname);
  
  if (bitspp == 8){
    ubyte *imagebuf;
    RGBQUAD *pal;
    outmap = FreeImage_AllocateT(FIT_BITMAP,width,height,bitspp,0xFF,0xFF,0xFF);
    pal = FreeImage_GetPalette(outmap);
    for (i = 0; i < 256; i++) {
      pal[i].rgbRed = i;
      pal[i].rgbGreen = i;
      pal[i].rgbBlue = i;
    }
    i = 0;
    for (y=0; y<height; y++){      
      imagebuf = FreeImage_GetScanLine(outmap,y);
      for (x=0;x<width;x++,i++)
	 imagebuf[x]=im[i];
	
    }

  } else {
    unsigned short *imagebuf;
    outmap = FreeImage_AllocateT(FIT_UINT16,width,height,16,0xFFFF,0xFFFF,0xFFFF);
    i = 0;
    for (y=0; y<height; y++){      
      imagebuf = (unsigned short *)FreeImage_GetScanLine(outmap,y);
      for (x=0;x<width;x++,i++)
	 imagebuf[x]=im[i];
	
    }
  }
  FreeImage_Save(fif,outmap,fname,0); 
  FreeImage_Unload(outmap);

}

                                                    
Pixel *ReadPGM (char *fnm)
{
  FILE *f;
  int i;

  gval = ReadTIFF(fnm,&width,&depth);
  height=1;
  size = width*depth;
  size2D = width;

  srand(42);
  if (MULFACTOR>0){
    for (i=0; i<size;i++)
      gval[i] = gval[i]*MULFACTOR + (MULFACTOR-1)*(float)rand()/(float)RAND_MAX;
  } else {
    for (i=0; i<size;i++)
      gval[i] = gval[i]/(-MULFACTOR);
  }

  return gval;

}

Pixel FindMax(Pixel *im){
  Pixel max = im[0];
  int i;
  for (i = 0; i<size;i++){
     max = (im[i]>max) ? im[i] : max;
  }
  printf("max = %d\n", max);
  return max;
}

void WritePGM(char *fname, Pixel *gval)
{
  int i;
  if (MULFACTOR>0)
    for (i=0; i<size;i++)
      gval[i] = gval[i]/MULFACTOR;
  if (FindMax(gval)>255)
    WriteTIFF(fname,gval, width, depth, 16);
  else
    WriteTIFF(fname, gval, width, depth, 8);

  
} /* WritePGM */


void Connect(ulong x, ulong y) {

  ulong area = 0, area1 = 0;
  ulong h, z;

  x = levroot(x);
  y = levroot(y);
  if (gval[x] < gval[y]) {
    h=x; x=y; y=h;
  }
  while ((x!=y) && (y!=bottom)) {
    z = Par(x);
    if ((z!=bottom) && (gval[z]>=gval[y])) {
      node[x].Area += area;
      x = z;
    } else {
      area1 = node[x].Area + area;
      area = node[x].Area ;
      node[x].Area =area1;
      node[x].parent = y;
      x = y;
      y = z;
    }
  }
  if (y==bottom) {
    while(x!=bottom) {
      node[x].Area += area;
      x = Par(x);
    }
  }
}

void Fuse2(int self, int i) /* fuse regions [LWB(self), UPB(self+i-1)) and  [LWB(self+i), UPB(self+2i-1)) vertically */
{
  ulong lwb, upb;
  ulong p, xm, x, y;
   
  lwb = LWB(self);
  upb = UPB(self+2*i-1); 
  if (upb>size) upb=size;

  /* get the horizontal boundary */
  xm = LWB(self+i); 

  x = xm;
  p = x % width;
  if ((p>0) && (x-1>=lwb)) {
     y = x-1;
     Connect(x, y);
  }  

  for (p=0; p<width && x<upb; p++){
      if (x>=lwb+width) {
         y= x-width;
         Connect(x, y);
      }
      x++;
  }

  if (depth>1) {
    x = xm;
    for (p=0; p<size2D && x<upb; p++){
      if (x>=lwb+size2D) {
	y= x-size2D;
	Connect(x, y);
      }
      x++;
    }
  }

  /* levrootfix(lwb, upb); */
}


void Fuse(int self, int i) /* fuse regions [LWB(self), UPB(self+i-1)) and  [LWB(self+i), UPB(self+2i-1)) vertically */
{
  ulong p, q, x, y, z;

  Pixel prevmin, curmin, nextmin;

  /* get the horizontal boundary */ 

  p  = LWB(self+i);
  q = p - size2D;

  x = p % width;
  y = (p /width) % height;
  z = p /size2D;
  
  /*  printf("Region %d merger with %d: (%d,%d,%d)\n",self, self+i, x,y,z);*/

  for ( y = 0 ; y < height ; y++ ){
    Connect(p, q);
    prevmin = MIN(gval[p],gval[q]);
    p++;
    q++;
    curmin = MIN(gval[p],gval[q]);
    for ( x = 1 ; x < width-1 ; x++, p++, q++ ){
      nextmin = MIN(gval[p+1],gval[q+1]);
      if ((curmin > prevmin) && (curmin>=nextmin)){
	Connect(p, q);
      }
      prevmin = curmin;    
      curmin = nextmin;
    }  
    p++;
    q++;
    Connect(p, q);

  }

  /* levrootfix(lwb, upb); */
}


/**************** Concurrent construction and filter of Maxtree  ***********/

typedef struct { 
  int self;
  pQueue *thisqueue;
  pStack *thisstack;
} ThreadData;

ThreadData *MakeThreadData(int numthreads){
  ThreadData *data = malloc(numthreads *sizeof(ThreadData));
  int i;

  for (i=0; i<numthreads; i++){
    data[i].self=i;
    data[i].thisqueue= pQueueCreate( UPB(i)-LWB(i));
    data[i].thisstack= pStackCreate( UPB(i)-LWB(i), data[i].thisqueue->queue);
    printf("Thread %d, LWB = %d, UPB = %d\n",i,LWB(i),UPB(i));
  }    
  return(data);
}

void FreeThreadData(ThreadData *data, int numthreads)
{
  int i;
  for (i=0;i<numthreads; i++){
    pQueueDelete(data[i].thisqueue);
    pStackDelete(data[i].thisstack);
  }
  free(data);
}


void *ccaf(void *arg) {
  ThreadData *thdata = (ThreadData *) arg;
  int self = thdata->self, q, i;
  ulong x, area=0;


  Flood( self, thdata->thisqueue, thdata->thisstack, 
	 gval, width, height, depth, node);


  i = 1;
  q = self;
  if (self==0)
    built = times(&tstruct);

  while ((self+i<nthreads) && (q%2 == 0)) {
    Psa(self+i);  /* wait to glue with righthand neighbor */
    Fuse(self, i);
    i = 2*i;
    q = q/2;
  }
  if (self != 0) {
    Vsa(self);  /* signal lefthand neighbor */
  }
  Barrier(self);
  if (self==0)
    merged = times(&tstruct);
  MaxTreeAreaFilter(self, lambda);

  return NULL;
}


void BuildTreeAndFilter(ThreadData *thdata, int nthreads) {
  int thread;
  for (thread=0; thread<nthreads; ++thread) {
    pthread_create(threadID+thread, NULL, ccaf, (void *) (thdata + thread));
  }
  for (thread=0; thread<nthreads; ++thread) {
    pthread_join(threadID[thread], NULL);
  }
}

int main (int argc, char *argv[]) {
  
   char *imgfname, *outfname = "out.tif";
   ThreadData *thdata;
   int r;
   ulong i;
   long tickspersec = sysconf(_SC_CLK_TCK);  
   float musec;
   FreeImage_Initialise(0);

   if (argc<4)
   {
      printf("Usage: %s <nthreads> <input image> <lambda>  [output image] [MULFACTOR]\n", argv[0]);
      exit(0);
   }

   nthreads = MIN(atoi(argv[1]), MAXTHREADS);
   imgfname = argv[2];

   lambda = atof(argv[3]);
   if (argc>4)  outfname = argv[4];
   if (argc>5)  MULFACTOR = atoi(argv[5]);
 
   
   if (!ReadPGM(imgfname)) 
     return(-1);
   size2D = width;
   
   printf("Filtering image '%s' using attribute area with lambda=%f\n", imgfname, lambda);
   printf("Image: Width=%d Height=%d Depth=%d\n", width, height, depth);
   printf ("nthreads: %d\n", nthreads);
   
   node = calloc((size_t)size, sizeof(MaxNode));
   if (node==NULL) {
      fprintf(stderr, "out of memory! \n");
      free(gval);
      return(-1);
   }
 
   out =  malloc(size*sizeof(Pixel));
   if (out==NULL) {
     fprintf(stderr, "Can't create output image! \n");
     free(node);
     free(gval);
     return(-1);
   }
 
   thdata = MakeThreadData(nthreads); 
   printf("Data read, start filtering.\n");
   start = times(&tstruct);
   
   BuildTreeAndFilter(thdata,nthreads);
   
   finish = times(&tstruct);
   musec = (float)(built - start)/((float)tickspersec);

   printf("build time: %f s\n",musec);
   musec = (float)(merged - built)/((float)tickspersec);

   printf("merge time: %f s\n",musec);
   musec = (float)(finish - merged)/((float)tickspersec);

   printf("filtering time: %f s\n",musec);
   musec = (float)(finish - start)/((float)tickspersec);

   printf("wall-clock time: %f s\n",musec);
   

   printf("Parcount = %d\n",ParCount); 

   WritePGM(outfname,out);
   free(out);  
   if (r)  printf("Filtered image written to '%s'\n", outfname);
   
   FreeThreadData(thdata,nthreads);
   free(node);
   free(gval);
  FreeImage_DeInitialise();
   return(0);
} /* main */
