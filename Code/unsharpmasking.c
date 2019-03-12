
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

long ImageSize;
long NumPlanes;

typedef unsigned short greyval;
typedef unsigned char ubyte;

#define NUMLEVELS     65536

greyval **ORI;    /* Denotes the original ... in the unproceesed input image */
greyval **neg;    /* Denotes the processed input image */
greyval **pos;

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

void usm (int beta) {
  int i;
  for(i = 0; i < ImageSize; i++) {
    int pixel = (int) ((int)ORI[0][i] + pos[0][i]*beta - neg[0][i]*beta);
    if(pixel > NUMLEVELS-1) {
      ORI[0][i] = NUMLEVELS-1;
    } else if(pixel < 0) {
      ORI[0][i] = 0;
    } else { 
      ORI[0][i] = (greyval) pixel;
    }
  }
}


int main(int argc, char *argv[]) {
    
  FreeImage_Initialise(0);
  int beta = 1;
  if (argc < 5) {
    printf ("Usage: %s original positive negative beta\n", argv[0]);
  }
  
  ORI =ReadTIFF(argv[1], &ImageWidth, &ImageHeight, &NumPlanes);
  pos = ReadTIFF(argv[2], &ImageWidth, &ImageHeight, &NumPlanes);
  neg = ReadTIFF(argv[3], &ImageWidth, &ImageHeight, &NumPlanes);
  

  
  ImageSize = ImageWidth * ImageHeight;
  
  if(argc == 5) {
    beta = atoi(argv[4]);
  }
  
  usm(beta);
  
  FreeImage_DeInitialise();
  
  WriteTIFF("usm.tif", ImageWidth, ImageHeight, 16, NumPlanes, ORI);
  
  
}