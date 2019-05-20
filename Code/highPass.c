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

typedef unsigned long greyval;
typedef unsigned short greyvalold;
typedef unsigned char ubyte;

greyvalold N;

#define NUMLEVELS     65536

greyval **ORI;    /* Denotes the original ... in the unproceesed input image */
greyval **new;    /* Denotes the output image */



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

void WriteTIFF(char *fname, long width, long height, int bitspp, long numplanes, greyvalold **im){
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
    greyvalold *imagebuf;
    if (numplanes==1)
      outmap = FreeImage_AllocateT(FIT_UINT16,width,height,16,0xFFFF,0xFFFF,0xFFFF);
    else   
      outmap = FreeImage_AllocateT(FIT_RGB16,width,height,48,0xFFFF,0xFFFF,0xFFFF);

    i = 0;
    for (j=height-1; j>=0; j--){      
      imagebuf = (greyvalold *)FreeImage_GetScanLine(outmap,j);
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

/* Find maximum of array*/

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
  

  int i = 0, j = 0;
  double alpha = (double) -(((-kh)*(-kh)) + ((-kh)*(-kh)));
  double div = exp(alpha/(2*sigma*sigma));
  for(i = 0; i < k; i++) {
    for (j = 0; j < k; j++) {
      alpha = (double) -(((i-kh)*(i-kh)) + ((j-kh)*(j-kh)));
      kernel[i*k+j] = exp(alpha/(2*sigma*sigma))/div;
    }
  }
  /*for(i = 0; i < k; i++) {
    for(j = 0; j < k; j++) {
      printf("%lf ", kernel[i*k + j]);
    }
    printf("\n");
  }*/ 
  
  return kernel;
}

/* Extending image edges for easier gaussian blurring */


greyval **extendImage(greyval **img) {
  greyval **eImage = calloc((size_t)NumPlanes, sizeof(greyval *));
  assert(eImage != NULL);
  long imsize = (eWidth) * (eHeight);
  int i,j, plane;
  for (i=0;i<NumPlanes;i++){
    eImage[i] = calloc((size_t)imsize, sizeof(greyval));
    assert(eImage[i] != NULL);
  }
  
  for (plane = 0; plane < NumPlanes; plane++) {
      
    //Horizontal Edges
    
    for(i = ksize/2; i < eWidth - ksize/2; i++) {
      for(j = 0 ; j < ksize/2;j++) {
        eImage[plane][i + j*eWidth] = img[plane][i-ksize/2];
      } 
    }  
    
    for(i = ksize/2; i < eWidth - ksize/2; i++) {
      for(j = 0 ; j < ksize/2;j++) {
        eImage[plane][i + (j+eHeight-ksize/2)*eWidth] = img[plane][i+(ImageHeight-1)*ImageWidth-ksize/2];
      } 
    }
    
    
    //Vertical edges
    
    for(i = ksize/2; i < eHeight - ksize/2; i++) {
      for(j = 0 ; j < ksize/2;j++) {
        eImage[plane][i*eWidth + j] = img[plane][ImageWidth*(i-ksize/2)];
      } 
    } 
    
    for(i = ksize/2; i < eHeight - ksize/2; i++) {
      for(j = eWidth -ksize/2 ; j < eWidth;j++) {
        eImage[plane][i*eWidth + j] = img[plane][ImageWidth*(i-ksize/2) +ImageWidth -1];
      } 
    }
    
    //Rest of image;
    
    for(i = ksize/2; i < eHeight- ksize/2; i++) {
      for(j = ksize/2; j < eWidth - ksize/2; j++) {
        eImage[plane][i * eWidth + j] = img[plane][(i-ksize/2)*ImageWidth + j-ksize/2];
      }
    }
  }
  return eImage;
}

/* Gaussian blurring of image */

void gaussianBlur() {
  int i,j,k, l, plane;
  int kh = ksize/2;
  greyval **eImage = extendImage(ORI);
  long imsize = ImageWidth*ImageHeight;
  new = calloc((size_t)NumPlanes, sizeof(greyval *));

  assert(new!=NULL);
  for (i=0;i<NumPlanes;i++){
    new[i] = calloc((size_t)imsize, sizeof(greyval));
    assert(new[i]!=NULL);
  }
  double *kernel = generateKernel(getStdDev());
  double sum = sumArray(kernel,ksize*ksize);
  //printf("sum:%lf\n", sum);
  
  /* Blurring is happening here */
  
  
  for (plane = 0; plane < NumPlanes; plane++) {
      
      
    // Horizontal
    for(i = 0; i < ImageHeight; i++) {
      for(j = 0; j < ImageWidth; j++) {    
        for(k = 0; k < ksize; k++) {
          
          new[plane][i*ImageWidth + j] += kernel[k] * eImage[plane][(i+kh)*eWidth + k+j];
          
        }
      }      
    }
    eImage = extendImage(new); //extending resulting image again
    // Vertical
     for(i = 0; i < ImageHeight; i++) {
      for(j = 0; j < ImageWidth; j++) {   
        new[plane][i*ImageWidth + j] = 0;
        for(k = 0; k < ksize; k++) {
            
          new[plane][i*ImageWidth + j] += kernel[k] * eImage[plane][i*eWidth + k*eWidth +kh+j];
          
        }
        new[plane][i*ImageWidth+j]/=sum; // Normalizing
      }
    }
  }
    
    for(plane=0; plane<NumPlanes; plane++){
     free(eImage[plane]);
   }
  free(eImage);
  free(kernel);
  }

/* Original image minus low pass filter */


void highPass(int beta) {
  greyvalold **pos, **neg;
  long i, plane;
  long imsize = ImageWidth*ImageHeight;
  
  //Positive values
  pos = calloc((size_t)NumPlanes, sizeof(greyvalold *));    
  assert(pos!=NULL);
  for (i=0;i<NumPlanes;i++){
    pos[i] = calloc((size_t)imsize, sizeof(greyvalold));
    assert(pos[i]!=NULL);
  }
  
  //Negative values
  
  neg = calloc((size_t)NumPlanes, sizeof(greyvalold *));
  assert(neg!=NULL);
  for (i=0;i<NumPlanes;i++){
    neg[i] = calloc((size_t)imsize, sizeof(greyvalold));
    assert(neg[i]!=NULL);
  }
  
  //High-pass
  for (plane = 0; plane < NumPlanes; plane++) {
    for(i = 0; i <ImageWidth*ImageHeight; i++) {
      int pixel =  (ORI[plane][i] - new[plane][i]);
      if(pixel < 0) {
        neg[plane][i] = (pixel*beta*-1 > N) ? N :(greyvalold) pixel*-1*beta;
      } else {
        pos[plane][i] = (pixel*beta > N) ? N : (greyvalold) pixel*beta;
      }
    }
  }
  
  //Writing images
  
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
    printf("Usage: %s <inputfile> <ksize> [beta]\n", argv[0]);
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
  highPass(beta);

 // printf("ImageHeight: %ld, ImageWidth: %ld\n", ImageHeight, ImageWidth);
  
  for(plane=0; plane<NumPlanes; plane++){
     free(ORI[plane]);
     free(new[plane]);
   }
  
  free(ORI);
  free(new);
    
  FreeImage_DeInitialise();
}