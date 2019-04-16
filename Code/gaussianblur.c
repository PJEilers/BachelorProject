#include <assert.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <sys/times.h>
#include <unistd.h>
#include <FreeImage.h>
#include <math.h>

long ImageWidth;
long ImageHeight;
long eWidth;
long eHeight;
long ImageSize;
long NumPlanes;
int ksize;

typedef unsigned short greyval;
typedef unsigned char ubyte;
greyval N;

#define NUMLEVELS     65536

greyval **ORI;    /* Denotes the original ... in the unproceesed input image */
greyval **new;    /* Denotes the processed input image */
greyval **eImage;



/** Generic image loader
@param lpszPathName Pointer to the full file name
@param flag Optional load flag constant
@return Returns the loaded dib if successful, returns NULL otherwise
*/
FIBITMAP* GenericLoader(const char* lpszPathName, int flag) {
  FREE_IMAGE_FORMAT fif = FIF_UNKNOWN;
/*  check the file signature and deduce its format */
/* (the second argument is currently not used by FreeImage) */
  fif = FreeImage_GetFileType(lpszPathName, 0);
  if(fif == FIF_UNKNOWN) {
    /* no signature ? */
    /* try to guess the file format from the file extension */
    fif = FreeImage_GetFIFFromFilename(lpszPathName);
  }
  /* check that the plugin has reading capabilities ... */
  if((fif != FIF_UNKNOWN) && FreeImage_FIFSupportsReading(fif)) {
    /* ok, let's load the file */
    FIBITMAP *dib = FreeImage_Load(fif, lpszPathName, flag);
    /* unless a bad file format, we are done ! */
    return dib;
  }
  return NULL;
}

void WriteTIFF(char *fname, long width, long height, int bitspp, long numplanes, greyval **im){
  FIBITMAP *outmap;
  long i,j,y,x,plane; 
  FREE_IMAGE_FORMAT fif = FreeImage_GetFIFFromFilename(fname);

  if (bitspp == 8){
    ubyte *imagebuf;
    RGBQUAD *pal;
    outmap = FreeImage_AllocateT(FIT_BITMAP,width,height,numplanes*bitspp,0xFF,0xFF,0xFF);
    if (numplanes==1){
      pal = FreeImage_GetPalette(outmap);
      for (i = 0; i < 256; i++) {
	pal[i].rgbRed = i;
	pal[i].rgbGreen = i;
	pal[i].rgbBlue = i;
      }
    }
    i = 0;
    for (j=height-1; j>=0; j--){      
      imagebuf = FreeImage_GetScanLine(outmap,j);
      for (x=0;x<width;x++,i++)
        for (plane=0;plane<numplanes;plane++)
	 imagebuf[numplanes*x+plane]=im[plane][i];	
    }

  } else {
    greyval *imagebuf;
    if (numplanes==1)
      outmap = FreeImage_AllocateT(FIT_UINT16,width,height,16,0xFFFF,0xFFFF,0xFFFF);
    else   
      outmap = FreeImage_AllocateT(FIT_RGB16,width,height,48,0xFFFF,0xFFFF,0xFFFF);

    i = 0;
    for (j=height-1; j>=0; j--){      
      imagebuf = (greyval *)FreeImage_GetScanLine(outmap,j);
      for (x=0;x<width;x++,i++)
        for (plane=0;plane<numplanes;plane++)
	 imagebuf[numplanes*x+plane]=im[plane][i];	
	
    }
  }

  FreeImage_Save(fif,outmap,fname,0); 
  FreeImage_Unload(outmap);

}


greyval **ReadTIFF(char *fnm, long *width, long *height, long *numplanes){
   FIBITMAP *dib = GenericLoader(fnm,0);
   unsigned long  bitsperpixel;
   greyval **im;
   unsigned int x,y,i,j,imsize;
   if (dib == NULL) return NULL;
     
   bitsperpixel =  FreeImage_GetBPP(dib);
   *numplanes=1;
   if ((bitsperpixel==24) ||(bitsperpixel==48))
     *numplanes=3;
   *height = FreeImage_GetHeight(dib), 
   *width = FreeImage_GetWidth(dib);
   
   imsize = (*width)*(*height);
   im = calloc((size_t)*numplanes, sizeof(greyval *));
   assert(im!=NULL);
   for (i=0;i<*numplanes;i++){
     im[i] = calloc((size_t)imsize, sizeof(greyval));
     assert(im[i]!=NULL);
   }
   switch(bitsperpixel) {
   case 8 :
     i=0;
     for(y = 0; y <*height; y++) {
       BYTE *bits = (BYTE *)FreeImage_GetScanLine(dib, *height - y -1);
       /*printf("y=%d\n",y);*/
       for(x = 0; x < *width; x++,i++) {
	 im[0][i] = bits[x];
       }
     }
     
     FreeImage_Unload(dib);
     return im;
   case 16 :
     i=0;
     for(y = 0; y < *height; y++) {
       unsigned short *bits = (unsigned short *)FreeImage_GetScanLine(dib,*height - y -1);
       for(x = 0; x < *width; x++,i++) {
	 im[0][i] = bits[x];
       }
     }
     FreeImage_Unload(dib);
     return im; 
   case 24 : 
     i=0;
     for(y = 0; y <*height; y++) {
       BYTE *bits = (BYTE *)FreeImage_GetScanLine(dib, *height - y -1);
       /*printf("y=%d\n",y);*/
       for(x = 0; x < *width; x++,i++) {
         for (j=0;j<*numplanes;j++)
	   im[j][i] = bits[3*x+j];
       }
     }
     
     FreeImage_Unload(dib);
     return im;
   case 48 : 
     i=0;
     for(y = 0; y <*height; y++) {
       unsigned short *bits = (unsigned short *)FreeImage_GetScanLine(dib, *height - y -1);
       /*printf("y=%d\n",y); */
       for(x = 0; x < *width; x++,i++) {
         for (j=0;j<*numplanes;j++)
	   im[j][i] = bits[3*x+j];
       }
     }
     
     FreeImage_Unload(dib);
     return im;
     
   default : 
     FreeImage_Unload(dib);
     fprintf(stderr, "unsupported format\n"); exit(-1);
     return NULL; 
   }
 
}

greyval findMax(greyval *buf, long width, long height){
  greyval max = buf[0];
  long i,size = width*height;
  for (i = 1; i<size;i++){
    max = (buf[i]>max) ? buf[i] : max;
  }
  return max;
}

/* Getting standard deviation */

double getStdDev () {;
  return 0.3*((ksize-1)*0.5 - 1) + 0.8;
}

/* Sum of array */

double sumArray(double *arr, int size) {
  int i;
  double sum = 0;
  for (i = 0; i < size; i++){
    sum+= arr[i];
  }
  return sum;
}

/** Generating kernel used for gaussian blurring
 * @param sigma, the standard deviation
 * @returns the kernel
 */


double *generateKernel(double sigma) {
  int k = ksize;
  int kh = ksize/2;
  double *kernel = malloc(k*k*sizeof(double));
  printf("half k:%d\n", kh);
  

  int i = 0, j = 0;
  
  for(i = 0; i < k; i++) {
    for (j = 0; j < k; j++) {
      double alpha = (double) -(((i-kh)*(i-kh)) + ((j-kh)*(j-kh)));
      kernel[i*k+j] = exp(alpha/(2*sigma*sigma));
    }
  }
  
    double sum = sumArray(kernel, k*k);
  
  for(i = 0; i < k*k; i ++) {
    kernel[i] /= sum;
  }
 /* for(i = 0; i < k; i++) {
    for(j = 0; j < k; j++) {
      printf("%lf ", kernel[i*k + j]);
    }
    printf("\n");
  } */
  
  return kernel;
}

/* Extending image corners and edges for easier gaussian blurring */


void extendImage() {
  eImage = calloc((size_t)NumPlanes, sizeof(greyval *));
  assert(eImage != NULL);
  long imsize = (eWidth) * (eHeight);
  int i,j;
  for (i=0;i<NumPlanes;i++){
    eImage[i] = calloc((size_t)imsize, sizeof(greyval));
    assert(eImage[i] != NULL);
  }
  
  //Top left
  
  for(i = 0; i < ksize/2; i++) {
    for(j = 0; j < ksize/2; j++) {
      eImage[0][i*(eWidth) + j] = ORI[0][0];

    }
  }
  
  //Top right
  
  for(i = 0; i < ksize/2; i++) {
    for(j = eWidth - ksize/2; j < eWidth; j++) {
      eImage[0][i*(eWidth) + j] = ORI[0][ImageWidth-1];

    }
  }
  
  //Bottom left
  
  for(i = eHeight - ksize/2; i < eHeight; i++) {
    for(j = 0; j < ksize/2; j++) {
      eImage[0][i*(eWidth) + j] = ORI[0][(ImageHeight-1)*ImageWidth];

    }
  }
  
  //Bottom right
  
  for(i = eHeight - ksize/2; i < eHeight; i++) {
    for(j = eWidth - ksize/2; j < eWidth; j++) {
      eImage[0][i*(eWidth) + j] = ORI[0][ImageHeight*ImageWidth-1];
    }
  }
  
  //Horizontal Edges
  
  for(i = ksize/2; i < eWidth - ksize/2; i++) {
    for(j = 0 ; j < ksize/2;j++) {
      eImage[0][i + j*eWidth] = ORI[0][i-ksize/2];
    } 
  }  
  
  for(i = ksize/2; i < eWidth - ksize/2; i++) {
    for(j = 0 ; j < ksize/2;j++) {
      eImage[0][i + (j+eHeight-ksize/2)*eWidth] = ORI[0][i+(ImageHeight-1)*ImageWidth-ksize/2];
    } 
  }
  
  
  //Vertical edges
  
  for(i = ksize/2; i < eHeight - ksize/2; i++) {
    for(j = 0 ; j < ksize/2;j++) {
      eImage[0][i*eWidth + j] = ORI[0][ImageWidth*(i-ksize/2)];
    } 
  } 
  
  for(i = ksize/2; i < eHeight - ksize/2; i++) {
    for(j = eWidth -ksize/2 ; j < eWidth;j++) {
      eImage[0][i*eWidth + j] = ORI[0][ImageWidth*(i-ksize/2) +ImageWidth -1];
    } 
  }
  
  //Rest of image;
  
  for(i = ksize/2; i < eHeight- ksize/2; i++) {
    for(j = ksize/2; j < eWidth - ksize/2; j++) {
      eImage[0][i * eWidth + j] = ORI[0][(i-ksize/2)*ImageWidth + j-ksize/2];
    }
  }
  
}

/* Gaussian blurring of image */

void gaussianBlur() {
  int i,j,k, l;
  extendImage();
  long imsize = ImageWidth*ImageHeight;
  new = calloc((size_t)NumPlanes, sizeof(greyval *));
  assert(new!=NULL);
  for (i=0;i<NumPlanes;i++){
    new[i] = calloc((size_t)imsize, sizeof(greyval));
    assert(new[i]!=NULL);
  }
  double *kernel = generateKernel(getStdDev());
  
  /* Blurring is happeninng here */
  
  for(i = 0; i < ImageHeight; i++) {
    for(j= 0; j < ImageWidth; j++) {
      for(k = 0; k < ksize; k++) {
        for(l = 0; l < ksize; l++) {
          new[0][i*ImageWidth+j] += kernel[k*ksize + l] * eImage[0][i*eWidth + j +  l + k*eWidth];
       //   if(j == 82 && i < 1) printf("%d %lf %d %d\n", ORI[0][i*ImageWidth + j], kernel[k*ksize + l], i*eWidth + j +  l + k*eWidth, eWidth);
        }
      }
    }
  }
  free(kernel);
}

/* Original image minus low pass filter */


void subtractImages(int beta) {
  greyval **pos, **neg;
  long i, plane;
  long imsize = ImageWidth*ImageHeight;
  
  //Positive values
  
  pos = calloc((size_t)NumPlanes, sizeof(greyval *));
  assert(pos!=NULL);
  for (i=0;i<NumPlanes;i++){
    pos[i] = calloc((size_t)imsize, sizeof(greyval));
    assert(pos[i]!=NULL);
  }
  
  //Negative values
  
  neg = calloc((size_t)NumPlanes, sizeof(greyval *));
  assert(neg!=NULL);
  for (i=0;i<NumPlanes;i++){
    neg[i] = calloc((size_t)imsize, sizeof(greyval));
    assert(neg[i]!=NULL);
  }
  
  for(i = 0; i <ImageWidth*ImageHeight; i++) {
    int pixel = (int) (ORI[0][i] - new[0][i]);
    if(pixel < 0) {
      neg[0][i] = (pixel*beta*-1 > N) ? N :(greyval) pixel*-1*beta;
    } else {
      pos[0][i] = (pixel*beta > N) ? N : (greyval) pixel*beta;
    }
  }
  if ( N> 255) {
    WriteTIFF("pos.tif", ImageWidth, ImageHeight, 16, NumPlanes, pos);
    WriteTIFF("neg.tif", ImageWidth, ImageHeight, 16, NumPlanes, neg);
  } else {
    WriteTIFF("pos.tif", ImageWidth, ImageHeight, 8, NumPlanes, pos);
    WriteTIFF("neg.tif", ImageWidth, ImageHeight, 8, NumPlanes, neg);
    
  }
  
  
  for(plane=0; plane<NumPlanes; plane++){
    free(pos[plane]);
    free(neg[plane]);
  }
  
  free(pos);
  free(neg);
  
}

int main(int argc, char *argv[]) {
  
  FreeImage_Initialise(0);
  int i,plane, beta =1;
  
  if (argc < 3)  {
    printf("Usage: %s inputfile ksize [beta][\n", argv[0]);
    exit(0);
  }
  ksize = atoi(argv[2]);
  if(ksize % 2 == 0) {
    printf("ksize must be odd\n");
    exit(0);
  }
  
  if (argc  == 4) {
    beta = atoi(argv[3]);
  }
  
  long width, height;
  ORI = ReadTIFF(argv[1], &width, &height, &NumPlanes);
  N = (findMax(ORI[0], width, height) > 255 ? NUMLEVELS-1 : 255);
  ImageWidth = width;
  ImageHeight = height;
  eWidth = width + ksize -1;
  eHeight = height + ksize -1;
  gaussianBlur();
  subtractImages(beta);

 // printf("ImageHeight: %ld, ImageWidth: %ld\n", ImageHeight, ImageWidth);
  
  for(plane=0; plane<NumPlanes; plane++){
     free(ORI[plane]);
     free(eImage[plane]);
     free(new[plane]);
   }
  
  free(ORI);
  free(eImage);
  free(new);
    
  FreeImage_DeInitialise();
}