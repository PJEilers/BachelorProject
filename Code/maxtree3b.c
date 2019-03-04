/* maxtree3b.c
 * October 25, 2005  Erik R. Urbach
 * Email: erik@cs.rug.nl
 * Max-tree with a single attribute parameter and an optional template
 * Attribute: I/(A^2) (default) and others (start program without arguments for
 *            a complete list of attributes available)
 *            Attributes here can use gray value information of the pixels used
 * Decision: Min, Direct, Max, Subtractive (default)
 * Input images: raw (P5) and plain (P2) PGM 8-bit gray-scale images
 * Output image: raw (P5) PGM 8-bit gray-scale image
 * Compilation: gcc -ansi -pedantic -Wall -O3 -o maxtree3b maxtree3b.c -lm
 *
 * Related papers:
 * [1] E. J. Breen and R. Jones.
 *     Attribute openings, thinnings and granulometries.
 *     Computer Vision and Image Understanding.
 *     Vol.64, No.3, Pages 377-389, 1996.
 * [2] P. Salembier and A. Oliveras and L. Garrido.
 *     Anti-extensive connected operators for image and sequence processing.
 *     IEEE Transactions on Image Processing,
 *     Vol.7, Pages 555-570, 1998.
 * [3] E. R. Urbach and M. H. F. Wilkinson.
 *     Shape-Only Granulometries and Grey-Scale Shape Filters.
 *     Proceedings of the ISMM2002,
 *     Pages 305-314, 2002.
 * [4] E. R. Urbach and J. B. T. M. Roerdink and M. H. F. Wilkinson.
 *     Connected Rotation-Invariant Size-Shape Granulometries.
 *     Proceedings of the 17th Int. Conf. Pat. Rec.,
 *     Vol.1, Pages 688-691, 2004.
 */

#include <math.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <sys/types.h>
#include <sys/times.h>
#include <unistd.h>



#define NUMLEVELS     256
#define CONNECTIVITY  4
#define PI 3.14159265358979323846



typedef short bool;
#define false 0
#define true  1

typedef unsigned char ubyte;
typedef unsigned int uint;
typedef unsigned long ulong;


#define MIN(a,b)  ((a<=b) ? (a) : (b))
#define MAX(a,b)  ((a>=b) ? (a) : (b))



typedef struct ImageGray ImageGray;
clock_t start,built, merged,finish;
struct tms tstruct;

struct ImageGray
{
   ulong Width;
   ulong Height;
   ubyte *Pixmap;
};



typedef struct HQueue
{
   ulong *Pixels;
   ulong Head;
   ulong Tail; /* First free place in queue, or -1 when the queue is full */
} HQueue;



typedef struct MaxNode MaxNode;

struct MaxNode
{
   ulong Parent;
   ulong Area;
   void *Attribute;
   ubyte Level;
   ubyte NewLevel;  /* gray level after filtering */
};



/* Status stores the information of the pixel status: the pixel can be
 * NotAnalyzed, InTheQueue or assigned to node k at level h. In this
 * last case Status(p)=k. */
#define ST_NotAnalyzed  -1
#define ST_InTheQueue   -2

typedef struct MaxTree MaxTree;

struct MaxTree
{
   long *Status;
   ulong *NumPixelsBelowLevel;
   ulong *NumNodesAtLevel; /* Number of nodes C^k_h at level h */
   MaxNode *Nodes;
   void *(*NewAuxData)(ulong, ulong, int, ulong *, ImageGray *);
   void (*AddToAuxData)(void *, ulong, ulong, int, ulong *, ImageGray *);
   void (*MergeAuxData)(void *, void *);
   void (*DeleteAuxData)(void *);
};



void MaxTreeDelete(MaxTree *mt);



/****** Typedefs and functions for area attributes ******************************/

typedef struct AreaData
{
   ulong Area;
} AreaData;

void *NewAreaData(ulong x, ulong y, int numneighbors, ulong *neighbors, ImageGray *img)
{
   AreaData *areadata;

   areadata = malloc(sizeof(AreaData));
   areadata->Area = 1;
   return(areadata);
} /* NewAreaData */

void DeleteAreaData(void *areaattr)
{
   free(areaattr);
} /* DeleteAreaData */

void AddToAreaData(void *areaattr, ulong x, ulong y, int numneighbors, ulong *neighbors, ImageGray *img)
{
   AreaData *areadata = areaattr;

   areadata->Area ++;
} /* AddToAreaData */

void MergeAreaData(void *areaattr, void *childattr)
{
   AreaData *areadata = areaattr;
   AreaData *childdata = childattr;

   areadata->Area += childdata->Area;
} /* MergeAreaData */

double AreaAttribute(void *areaattr)
{
   AreaData *areadata = areaattr;
   double area;

   area = areadata->Area;
   return(area);
} /* AreaAttribute */



/****** Typedefs and functions for minimum enclosing rectangle attributes *******/

typedef struct EnclRectData
{
   ulong MinX;
   ulong MinY;
   ulong MaxX;
   ulong MaxY;
} EnclRectData;

void *NewEnclRectData(ulong x, ulong y, int numneighbors, ulong *neighbors, ImageGray *img)
{
   EnclRectData *rectdata;

   rectdata = malloc(sizeof(EnclRectData));
   rectdata->MinX = rectdata->MaxX = x;
   rectdata->MinY = rectdata->MaxY = y;
   return(rectdata);
} /* NewEnclRectData */

void DeleteEnclRectData(void *rectattr)
{
   free(rectattr);
} /* DeleteEnclRectData */

void AddToEnclRectData(void *rectattr, ulong x, ulong y, int numneighbors, ulong *neighbors, ImageGray *img)
{
   EnclRectData *rectdata = rectattr;

   rectdata->MinX = MIN(rectdata->MinX, x);
   rectdata->MinY = MIN(rectdata->MinY, y);
   rectdata->MaxX = MAX(rectdata->MaxX, x);
   rectdata->MaxY = MAX(rectdata->MaxY, y);
} /* AddToEnclRectData */

void MergeEnclRectData(void *rectattr, void *childattr)
{
   EnclRectData *rectdata = rectattr;
   EnclRectData *childdata = childattr;

   rectdata->MinX = MIN(rectdata->MinX, childdata->MinX);
   rectdata->MinY = MIN(rectdata->MinY, childdata->MinY);
   rectdata->MaxX = MAX(rectdata->MaxX, childdata->MaxX);
   rectdata->MaxY = MAX(rectdata->MaxY, childdata->MaxY);
} /* MergeEnclRectData */

double EnclRectAreaAttribute(void *rectattr)
{
   EnclRectData *rectdata = rectattr;
   double area;

   area = (rectdata->MaxX - rectdata->MinX + 1)
        * (rectdata->MaxY - rectdata->MinY + 1);
   return(area);
} /* EnclRectAreaAttribute */

double EnclRectDiagAttribute(void *rectattr)
/* Computes the square of the length of the diagonal */
{
   EnclRectData *rectdata = rectattr;
   double minx, miny, maxx, maxy, l;

   minx = rectdata->MinX;
   miny = rectdata->MinY;
   maxx = rectdata->MaxX;
   maxy = rectdata->MaxY;
   l = (maxx-minx+1) * (maxx-minx+1)
     + (maxy-miny+1) * (maxy-miny+1);
   return(l);
} /* EnclRectDiagAttribute */



/****** Typedefs and functions for cityblock (huge) perimeter attributes *********************/

typedef struct PeriCBData
{
   ulong Area;
   ulong Perimeter;
} PeriCBData;

void *NewPeriCBData(ulong x, ulong y, int numneighbors, ulong *neighbors, ImageGray *img)
{
   PeriCBData *peridata;
   ulong lx=x, ly=y, p, peri=0, q;
   int i;
   ubyte h;

   p = ly*(img->Width) + lx;
   h = img->Pixmap[p];
   peri += CONNECTIVITY-numneighbors;
   for (i=0; i<numneighbors; i++)
   {
      q = neighbors[i];
      if (img->Pixmap[q]<h)  peri++;
      if (img->Pixmap[q]>h)  peri--;
   }
   peridata = malloc(sizeof(PeriCBData));
   peridata->Area = 1;
   peridata->Perimeter = peri;
   return(peridata);
} /* NewPeriCBData */

void DeletePeriCBData(void *periattr)
{
   free(periattr);
} /* DeletePeriCBData */

void AddToPeriCBData(void *periattr, ulong x, ulong y, int numneighbors, ulong *neighbors, ImageGray *img)
{
   PeriCBData *peridata = periattr;
   ulong lx=x, ly=y, p, peri=0, q;
   int i;
   ubyte h;

   p = ly*(img->Width) + lx;
   h = img->Pixmap[p];
   peri += CONNECTIVITY-numneighbors;
   for (i=0; i<numneighbors; i++)
   {
      q = neighbors[i];
      if (img->Pixmap[q]<h)  peri++;
      if (img->Pixmap[q]>h)  peri--;
   }
   peridata->Area ++;
   peridata->Perimeter += peri;
} /* AddToPeriCBData */

void MergePeriCBData(void *periattr, void *childattr)
{
   PeriCBData *peridata = periattr;
   PeriCBData *childdata = childattr;

   peridata->Area += childdata->Area;
   peridata->Perimeter += childdata->Perimeter;
} /* MergePeriCBData */

double PeriCBAreaAttribute(void *periattr)
{
   PeriCBData *peridata = periattr;
   double area;

   area = peridata->Area;
   return(area);
} /* PeriCBAreaAttribute */

double PeriCBPerimeterAttribute(void *periattr)
{
   PeriCBData *peridata = periattr;
   double peri;

   peri = peridata->Perimeter;
   return(peri);
} /* PeriCBPerimeterAttribute */

double PeriCBComplexityAttribute(void *periattr)
{
   PeriCBData *peridata = periattr;
   double area, peri;

   area = peridata->Area;
   peri = peridata->Perimeter;
   return(peri/area);
} /* PeriCBComplexityAttribute */

double PeriCBSimplicityAttribute(void *periattr)
{
   PeriCBData *peridata = periattr;
   double area, peri;

   area = peridata->Area;
   peri = peridata->Perimeter;
   return(area/peri);
} /* PeriCBSimplicityAttribute */

double PeriCBCompactnessAttribute(void *periattr)
{
   PeriCBData *peridata = periattr;
   double area, peri;

   area = peridata->Area;
   peri = peridata->Perimeter;
   return((peri*peri)/(4.0*PI*area));
} /* PeriCBCompactnessAttribute */



int Get8NeighValues(ImageGray *img, ulong *neighbors, ulong x, ulong y)
{
   ulong p;
   int i;

   p = y*(img->Width) + x;
   for (i=0; i<8; i++)  neighbors[i] = -1;
   if (y>0)
   {
      neighbors[1] = img->Pixmap[p-(img->Width)];
      if (x>0)  neighbors[0] = img->Pixmap[p-(img->Width)-1];
      if (x<(img->Width)-1)  neighbors[2] = img->Pixmap[p-(img->Width)+1];
   }
   if (x<(img->Width)-1)  neighbors[3] = img->Pixmap[p+1];
   if (y<(img->Height)-1)
   {
      neighbors[5] = img->Pixmap[p+(img->Width)];
      if (x>0)  neighbors[6] = img->Pixmap[p+(img->Width)-1];
      if (x<(img->Width)-1)  neighbors[4] = img->Pixmap[p+(img->Width)+1];
   }
   if (x>0)  neighbors[7] = img->Pixmap[p-1];
   return(img->Pixmap[p]);
} /* Get8NeighValues */

/****** Typedefs and functions for large (8-conn.) perimeter attributes *********************/

typedef struct PeriLargeData
{
   ulong Area;
   double Perimeter;
} PeriLargeData;

double PeriLargeCalcSide(int h, int a, int b, int c, int d, int e)
{
   double addarray[8] =
   {
      4.0, 4.0, 2.0, 2.0, 4.0, 2.0, 2.0, 2.7017787186
   };
   double l=0.0;

#define PeriLargeBit(v)  ((v)>=(h) ? 1 : 0)

   if ((b<h) && (c>h) && (e>=h))  l += sqrt(2.0) - 2.0;
   if ((b<h) && (a>h) && (d>=h))  l += sqrt(2.0) - 2.0;
   if (b>h)
   {
      if ((c<h) && (e>h))  l += sqrt(8.0) - 4.0;
      else if ((a>h) && (d>h))
      {
         if ((c>h) && (e>h))  l -= sqrt(20.0)-sqrt(8.0);
         else  l -= sqrt(8.0);
      }
      if ((a<=h) && (c<=h))  l -= 2.0;
      else if ((d<=h) && (e<=h))  l -= 2.0;
      else if ((a<=h) && (e<=h))  l -= 2.0;
      else if ((c<=h) && (d<=h))  l -= 2.0;
   }
   if (b>=h)  return(l);
   else
   {
      if (a<h)
      {
         if ((c>=h) && (e>=h))  return(l+sqrt(2.0));
         return(l+2.0);
      } else  return(l+sqrt(addarray[(PeriLargeBit(c)<<2)+(PeriLargeBit(d)<<1)+PeriLargeBit(e)]));
   }
} /* PeriLargeCalcSide */

void *NewPeriLargeData(ulong x, ulong y, int numneighbors, ulong *neighbors, ImageGray *img)
{
   PeriLargeData *peridata;
   double peri;
   ulong neigh8[8];
   ubyte h;

   h = Get8NeighValues(img, neigh8, x, y);
   peri = 0.5*PeriLargeCalcSide(h, neigh8[0], neigh8[1], neigh8[2], neigh8[7], neigh8[3]);
   peri += 0.5*PeriLargeCalcSide(h, neigh8[2], neigh8[3], neigh8[4], neigh8[1], neigh8[5]);
   peri += 0.5*PeriLargeCalcSide(h, neigh8[4], neigh8[5], neigh8[6], neigh8[3], neigh8[7]);
   peri += 0.5*PeriLargeCalcSide(h, neigh8[6], neigh8[7], neigh8[0], neigh8[5], neigh8[1]);
   peridata = malloc(sizeof(PeriLargeData));
   peridata->Area = 1;
   peridata->Perimeter = peri;
   return(peridata);
} /* NewPeriLargeData */

void DeletePeriLargeData(void *periattr)
{
   free(periattr);
} /* DeletePeriLargeData */

void AddToPeriLargeData(void *periattr, ulong x, ulong y, int numneighbors, ulong *neighbors, ImageGray *img)
{
   PeriLargeData *peridata = periattr;
   double peri;
   ulong neigh8[8];
   ubyte h;

   h = Get8NeighValues(img, neigh8, x, y);
   peri = 0.5*PeriLargeCalcSide(h, neigh8[0], neigh8[1], neigh8[2], neigh8[7], neigh8[3]);
   peri += 0.5*PeriLargeCalcSide(h, neigh8[2], neigh8[3], neigh8[4], neigh8[1], neigh8[5]);
   peri += 0.5*PeriLargeCalcSide(h, neigh8[4], neigh8[5], neigh8[6], neigh8[3], neigh8[7]);
   peri += 0.5*PeriLargeCalcSide(h, neigh8[6], neigh8[7], neigh8[0], neigh8[5], neigh8[1]);
   peridata->Area ++;
   peridata->Perimeter += peri;
} /* AddToPeriLargeData */

void MergePeriLargeData(void *periattr, void *childattr)
{
   PeriLargeData *peridata = periattr;
   PeriLargeData *childdata = childattr;

   peridata->Area += childdata->Area;
   peridata->Perimeter += childdata->Perimeter;
} /* MergePeriLargeData */

double PeriLargeAreaAttribute(void *periattr)
{
   PeriLargeData *peridata = periattr;
   double area;

   area = peridata->Area;
   return(area);
} /* PeriLargeAreaAttribute */

double PeriLargePerimeterAttribute(void *periattr)
{
   PeriLargeData *peridata = periattr;
   double peri;

   peri = peridata->Perimeter;
   return(peri);
} /* PeriLargePerimeterAttribute */

double PeriLargeCompactnessAttribute(void *periattr)
{
   PeriLargeData *peridata = periattr;
   double area, peri;

   area = peridata->Area;
   peri = peridata->Perimeter;
   return((peri*peri)/(4.0*PI*area));
} /* PeriLargeCompactnessAttribute */



/****** Typedefs and functions for small (8-conn.) perimeter attributes *********************/

typedef struct PeriSmallData
{
   ulong Area;
   double Perimeter;
} PeriSmallData;

double PeriSmallCalcSide(int h, int a, int b, int c)
{
   double addarray[8] =
   {
      0.0, 0.25, 0.0, 2.0/9.0, 0.25, 2.0/9.0, 2.0/9.0, 0.0
   };
   double l;

#define PeriSmallBit(v)  ((v)>=(h) ? 1 : 0)

   l = sqrt(addarray[(PeriSmallBit(a)<<2)+(PeriSmallBit(b)<<1)+PeriSmallBit(c)]);
   if (c>h)
   {
      if (b<h)
      {
         if (a<h)
	 {
	     return(l+1.0);
	 }
         else if (a>h)
	 {
             return(l+sqrt(2.0)*2.0/3.0);
	 }
      }
      else if (b>h)
      {
         if (a==h)
	 {
	   return(l-1.0);
	 }
         else if (a>h)
	 {
	   return(l-sqrt(2.0));
	 }
      }
   }
   return(l);
} /* PeriSmallCalcSide */

void *NewPeriSmallData(ulong x, ulong y, int numneighbors, ulong *neighbors, ImageGray *img)
{
   PeriSmallData *peridata;
   double peri;
   ulong neigh8[8];
   ubyte h;

   h = Get8NeighValues(img, neigh8, x, y);
   peri = PeriSmallCalcSide(h, neigh8[7], neigh8[0], neigh8[1]);
   peri += PeriSmallCalcSide(h, neigh8[1], neigh8[2], neigh8[3]);
   peri += PeriSmallCalcSide(h, neigh8[3], neigh8[4], neigh8[5]);
   peri += PeriSmallCalcSide(h, neigh8[5], neigh8[6], neigh8[7]);
   peridata = malloc(sizeof(PeriSmallData));
   peridata->Area = 1;
   peridata->Perimeter = peri;
   return(peridata);
} /* NewPeriSmallData */

void DeletePeriSmallData(void *periattr)
{
   free(periattr);
} /* DeletePeriSmallData */

void AddToPeriSmallData(void *periattr, ulong x, ulong y, int numneighbors, ulong *neighbors, ImageGray *img)
{
   PeriSmallData *peridata = periattr;
   double peri;
   ulong neigh8[8];
   ubyte h;

   h = Get8NeighValues(img, neigh8, x, y);
   peri = PeriSmallCalcSide(h, neigh8[7], neigh8[0], neigh8[1]);
   peri += PeriSmallCalcSide(h, neigh8[1], neigh8[2], neigh8[3]);
   peri += PeriSmallCalcSide(h, neigh8[3], neigh8[4], neigh8[5]);
   peri += PeriSmallCalcSide(h, neigh8[5], neigh8[6], neigh8[7]);
   peridata->Area ++;
   peridata->Perimeter += peri;
} /* AddToPeriSmallData */

void MergePeriSmallData(void *periattr, void *childattr)
{
   PeriSmallData *peridata = periattr;
   PeriSmallData *childdata = childattr;

   peridata->Area += childdata->Area;
   peridata->Perimeter += childdata->Perimeter;
} /* MergePeriSmallData */

double PeriSmallAreaAttribute(void *periattr)
{
   PeriSmallData *peridata = periattr;
   double area;

   area = peridata->Area;
   return(area);
} /* PeriSmallAreaAttribute */

double PeriSmallPerimeterAttribute(void *periattr)
{
   PeriSmallData *peridata = periattr;
   double peri;

   peri = peridata->Perimeter;
   return(peri);
} /* PeriSmallPerimeterAttribute */

double PeriSmallCompactnessAttribute(void *periattr)
{
   PeriSmallData *peridata = periattr;
   double area, peri;

   area = peridata->Area;
   peri = peridata->Perimeter;
   return((peri*peri)/(4.0*PI*area));
} /* PeriSmallCompactnessAttribute */



/****** Typedefs and functions for moment of inertia attributes **************************/

typedef struct InertiaData
{
   ulong Area;
   double SumX, SumY, SumX2, SumY2;
} InertiaData;

void *NewInertiaData(ulong x, ulong y, int numneighbors, ulong *neighbors, ImageGray *img)
{
   InertiaData *inertiadata;

   inertiadata = malloc(sizeof(InertiaData));
   inertiadata->Area = 1;
   inertiadata->SumX = x;
   inertiadata->SumY = y;
   inertiadata->SumX2 = x*x;
   inertiadata->SumY2 = y*y;
   return(inertiadata);
} /* NewInertiaData */

void DeleteInertiaData(void *inertiaattr)
{
   free(inertiaattr);
} /* DeleteInertiaData */

void AddToInertiaData(void *inertiaattr, ulong x, ulong y, int numneighbors, ulong *neighbors, ImageGray *img)
{
   InertiaData *inertiadata = inertiaattr;

   inertiadata->Area ++;
   inertiadata->SumX += x;
   inertiadata->SumY += y;
   inertiadata->SumX2 += x*x;
   inertiadata->SumY2 += y*y;
} /* AddToInertiaData */

void MergeInertiaData(void *inertiaattr, void *childattr)
{
   InertiaData *inertiadata = inertiaattr;
   InertiaData *childdata = childattr;

   inertiadata->Area += childdata->Area;
   inertiadata->SumX += childdata->SumX;
   inertiadata->SumY += childdata->SumY;
   inertiadata->SumX2 += childdata->SumX2;
   inertiadata->SumY2 += childdata->SumY2;
} /* MergeInertiaData */

double InertiaAttribute(void *inertiaattr)
{
   InertiaData *inertiadata = inertiaattr;
   double area, inertia;

   area = inertiadata->Area;
   inertia = inertiadata->SumX2 + inertiadata->SumY2 -
             (inertiadata->SumX * inertiadata->SumX +
              inertiadata->SumY * inertiadata->SumY) / area
             + area / 6.0;
   return(inertia);
} /* InertiaAttribute */

double InertiaDivA2Attribute(void *inertiaattr)
{
   InertiaData *inertiadata = inertiaattr;
   double inertia, area;

   area = (double)(inertiadata->Area);
   inertia = inertiadata->SumX2 + inertiadata->SumY2 -
             (inertiadata->SumX * inertiadata->SumX +
              inertiadata->SumY * inertiadata->SumY) / area
             + area / 6.0;
   return(inertia*2.0*PI/(area*area));
} /* InertiaDivA2Attribute */

double MeanXAttribute(void *inertiaattr)
{
   InertiaData *inertiadata = inertiaattr;
   double area, sumx;

   area = inertiadata->Area;
   sumx = inertiadata->SumX;
   return(sumx/area);
} /* MeanXAttribute */

double MeanYAttribute(void *inertiaattr)
{
   InertiaData *inertiadata = inertiaattr;
   double area, sumy;

   area = inertiadata->Area;
   sumy = inertiadata->SumY;
   return(sumy/area);
} /* MeanYAttribute */





/****** Typedefs and functions for jaggedness attributes *********************/

typedef struct JaggedData
{
   ulong Area;
   ulong Perimeter;
   double SumX, SumY, SumX2, SumY2;
} JaggedData;

void *NewJaggedData(ulong x, ulong y, int numneighbors, ulong *neighbors, ImageGray *img)
{
   JaggedData *jaggeddata;
   ulong lx=x, ly=y, p, peri=0, q;
   int i;
   ubyte h;

   p = ly*(img->Width) + lx;
   h = img->Pixmap[p];
   peri += CONNECTIVITY-numneighbors;
   for (i=0; i<numneighbors; i++)
   {
      q = neighbors[i];
      if (img->Pixmap[q]<h)  peri++;
      if (img->Pixmap[q]>h)  peri--;
   }
   jaggeddata = malloc(sizeof(JaggedData));
   jaggeddata->Area = 1;
   jaggeddata->Perimeter = peri;
   jaggeddata->SumX = x;
   jaggeddata->SumY = y;
   jaggeddata->SumX2 = x*x;
   jaggeddata->SumY2 = y*y;
   return(jaggeddata);
} /* NewJaggedData */

void DeleteJaggedData(void *jaggedattr)
{
   free(jaggedattr);
} /* DeleteJaggedData */

void AddToJaggedData(void *jaggedattr, ulong x, ulong y, int numneighbors, ulong *neighbors, ImageGray *img)
{
   JaggedData *jaggeddata = jaggedattr;
   ulong lx=x, ly=y, p, peri=0, q;
   int i;
   ubyte h;

   p = ly*(img->Width) + lx;
   h = img->Pixmap[p];
   peri += CONNECTIVITY-numneighbors;
   for (i=0; i<numneighbors; i++)
   {
      q = neighbors[i];
      if (img->Pixmap[q]<h)  peri++;
      if (img->Pixmap[q]>h)  peri--;
   }
   jaggeddata->Area ++;
   jaggeddata->Perimeter += peri;
   jaggeddata->SumX += x;
   jaggeddata->SumY += y;
   jaggeddata->SumX2 += x*x;
   jaggeddata->SumY2 += y*y;
} /* AddToJaggedData */

void MergeJaggedData(void *jaggedattr, void *childattr)
{
   JaggedData *jaggeddata = jaggedattr;
   JaggedData *childdata = childattr;

   jaggeddata->Area += childdata->Area;
   jaggeddata->Perimeter += childdata->Perimeter;
   jaggeddata->SumX += childdata->SumX;
   jaggeddata->SumY += childdata->SumY;
   jaggeddata->SumX2 += childdata->SumX2;
   jaggeddata->SumY2 += childdata->SumY2;
} /* MergeJaggedData */

double JaggedAttribute(void *jaggedattr)
{
   JaggedData *jaggeddata = jaggedattr;
   double peri;

   peri = jaggeddata->Perimeter;
   return(peri);
} /* JaggedPerimeterAttribute */

double JaggedCompactnessAttribute(void *jaggedattr)
{
   JaggedData *jaggeddata = jaggedattr;
   double area, peri;

   area = (double)(jaggeddata->Area);
   peri = jaggeddata->Perimeter;
   return((peri*peri)/(4.0*PI*area));
} /* JaggedCompactnessAttribute */

double JaggedInertiaDivA2Attribute(void *jaggedattr)
{
   JaggedData *jaggeddata = jaggedattr;
   double inertia, area;

   area = (double)(jaggeddata->Area);
   inertia = jaggeddata->SumX2 + jaggeddata->SumY2 -
             (jaggeddata->SumX * jaggeddata->SumX +
              jaggeddata->SumY * jaggeddata->SumY) / area
             + area / 6.0;
   return(inertia*2.0*PI/(area*area));
} /* JaggedInertiaDivA2Attribute */

double JaggednessAttribute(void *jaggedattr)
{
   JaggedData *jaggeddata = jaggedattr;
   double area, peri, inertia;

   area = (double)(jaggeddata->Area);
   peri = jaggeddata->Perimeter;
   inertia = jaggeddata->SumX2 + jaggeddata->SumY2 -
             (jaggeddata->SumX * jaggeddata->SumX +
              jaggeddata->SumY * jaggeddata->SumY) / area
             + area / 6.0;
   return(area*peri*peri/(8.0*PI*PI*inertia));
} /* JaggednessAttribute */



/****** Typedefs and functions for Entropy attributes ******************************/
/* TODO: check all attribute functions below for use of pixelsize */

typedef struct EntropyData
{
   ulong Hist[NUMLEVELS];
} EntropyData;

void *NewEntropyData(ulong x, ulong y, int numneighbors, ulong *neighbors, ImageGray *img)
{
   EntropyData *entropydata;
   ulong lx=x, ly=y, p;
   int i;

   p = ly*(img->Width) + lx;
   entropydata = malloc(sizeof(EntropyData));
   for (i=0; i<NUMLEVELS; i++)  entropydata->Hist[i] = 0;
   entropydata->Hist[img->Pixmap[p]] = 1;
   return(entropydata);
} /* NewEntropyData */

void DeleteEntropyData(void *entropyattr)
{
   free(entropyattr);
} /* DeleteEntropyData */

void AddToEntropyData(void *entropyattr, ulong x, ulong y, int numneighbors, ulong *neighbors, ImageGray *img)
{
   EntropyData *entropydata = entropyattr;
   ulong lx=x, ly=y, p;

   p = ly*(img->Width) + lx;
   entropydata->Hist[img->Pixmap[p]] ++;
} /* AddToEntropyData */

void MergeEntropyData(void *entropyattr, void *childattr)
{
   EntropyData *entropydata = entropyattr;
   EntropyData *childdata = childattr;
   int i;

   for (i=0; i<NUMLEVELS; i++)  entropydata->Hist[i] += childdata->Hist[i];
} /* MergeEntropyData */

double EntropyAttribute(void *entropyattr)
{
   EntropyData *entropydata = entropyattr;
   double p[NUMLEVELS];
   double num=0.0, entropy = 0.0;
   int i;

   for (i=0; i<NUMLEVELS; i++)  num += entropydata->Hist[i];
   for (i=0; i<NUMLEVELS; i++)  p[i] = (entropydata->Hist[i])/num;
   for (i=0; i<NUMLEVELS; i++)  entropy += p[i] * (log(p[i]+0.00001)/log(2.0));
   return(-entropy);
} /* EntropyAttribute */



/****** Typedefs and functions for lambda-max attributes *************************/

typedef struct LambdamaxData
{
   ubyte MinLevel;
   ubyte MaxLevel;
} LambdamaxData;

void *NewLambdamaxData(ulong x, ulong y, int numneighbors, ulong *neighbors, ImageGray *img)
{
   LambdamaxData *lambdadata;

   lambdadata = malloc(sizeof(LambdamaxData));
   lambdadata->MinLevel = img->Pixmap[y*(img->Width)+x];
   lambdadata->MaxLevel = img->Pixmap[y*(img->Width)+x];
   return(lambdadata);
} /* NewLambdamaxData */

void DeleteLambdamaxData(void *lambdaattr)
{
   free(lambdaattr);
} /* DeleteLambdamaxData */

void AddToLambdamaxData(void *lambdaattr, ulong x, ulong y, int numneighbors, ulong *neighbors, ImageGray *img)
{
} /* AddToLambdamaxData */

void MergeLambdamaxData(void *lambdaattr, void *childattr)
{
   LambdamaxData *lambdadata = lambdaattr;
   LambdamaxData *childdata = childattr;

   if (childdata->MaxLevel > lambdadata->MaxLevel)  lambdadata->MaxLevel = childdata->MaxLevel;
} /* MergeLambdamaxData */

double LambdamaxAttribute(void *lambdaattr)
{
   LambdamaxData *lambdadata = lambdaattr;
   double height;

   height = ((double)lambdadata->MaxLevel) - ((double)lambdadata->MinLevel);
   return(height);
} /* LambdamaxAttribute */



/****** Typedefs and functions for graylevel attributes ******************************/

typedef struct LevelData
{
   ubyte level;
} LevelData;

void *NewLevelData(ulong x, ulong y, int numneighbors, ulong *neighbors, ImageGray *img)
{
   LevelData *leveldata;

   leveldata = malloc(sizeof(LevelData));
   leveldata->level = img->Pixmap[y*(img->Width)+x];
   return(leveldata);
} /* NewLevelData */

void DeleteLevelData(void *levelattr)
{
   free(levelattr);
} /* DeleteLevelData */

void AddToLevelData(void *levelattr, ulong x, ulong y, int numneighbors, ulong *neighbors, ImageGray *img)
{
} /* AddToLevelData */

void MergeLevelData(void *levelattr, void *childattr)
{
} /* MergeLevelData */

double LevelAttribute(void *levelattr)
{
   LevelData *leveldata = levelattr;
   double level;

   level = leveldata->level;
   return(level);
} /* LevelAttribute */



/****** Image create/read/write functions ******************************/

ImageGray *ImageGrayCreate(ulong width, ulong height)
{
   ImageGray *img;

   img = malloc(sizeof(ImageGray));
   if (img==NULL)  return(NULL);
   img->Width = width;
   img->Height = height;
   img->Pixmap = malloc(width*height);
   if (img->Pixmap==NULL)
   {
      free(img);
      return(NULL);
   }
   return(img);
} /* ImageGrayCreate */



void ImageGrayDelete(ImageGray *img)
{
   free(img->Pixmap);
   free(img);
} /* ImageGrayDelete */



void ImageGrayInit(ImageGray *img, ubyte h)
{
   memset(img->Pixmap, h, (img->Width)*(img->Height));
} /* ImageGrayInit */



ImageGray *ImagePGMAsciiRead(char *fname)
{
   FILE *infile;
   ImageGray *img;
   ulong width, height, i;
   int c;

   infile = fopen(fname, "r");
   if (infile==NULL)  return(NULL);
   fscanf(infile, "P2\n");
   while ((c=fgetc(infile)) == '#')
      while ((c=fgetc(infile)) != '\n');
   ungetc(c, infile);
   fscanf(infile, "%lu %lu\n255\n", &width, &height);
   img = ImageGrayCreate(width, height);
   if (img==NULL)
   {
      fclose(infile);
      return(NULL);
   }
   for (i=0; i<width*height; i++)
   {
      fscanf(infile, "%d", &c);
      img->Pixmap[i] = c;
   }
   fclose(infile);
   return(img);
} /* ImagePGMAsciiRead */



ImageGray *ImagePGMBinRead(char *fname)
{
   FILE *infile;
   ImageGray *img;
   ulong width, height;
   int c;

   infile = fopen(fname, "rb");
   if (infile==NULL)  return(NULL);
   fscanf(infile, "P5\n");
   while ((c=fgetc(infile)) == '#')
      while ((c=fgetc(infile)) != '\n');
   ungetc(c, infile);
   fscanf(infile, "%lu %lu\n255\n", &width, &height);
   img = ImageGrayCreate(width, height);
   if (img)  fread(img->Pixmap, 1, width*height, infile);
   fclose(infile);
   return(img);
} /* ImagePGMBinRead */



ImageGray *ImagePGMRead(char *fname)
{
   FILE *infile;
   char id[4];

   infile = fopen(fname, "r");
   if (infile==NULL)  return(NULL);
   fscanf(infile, "%3s", id);
   fclose(infile);
   if (strcmp(id, "P2")==0)  return(ImagePGMAsciiRead(fname));
   else if (strcmp(id, "P5")==0)  return(ImagePGMBinRead(fname));
   else  return(NULL);
} /* ImagePGMRead */



int ImagePGMBinWrite(ImageGray *img, char *fname)
{
   FILE *outfile;

   outfile = fopen(fname, "wb");
   if (outfile==NULL)  return(-1);
   fprintf(outfile, "P5\n%ld %ld\n255\n", img->Width, img->Height);
   fwrite(img->Pixmap, 1, (size_t)((img->Width)*(img->Height)), outfile);
   fclose(outfile);
   return(0);
} /* ImagePGMBinWrite */



/****** Max-tree routines ******************************/

HQueue *HQueueCreate(ulong imgsize, ulong *numpixelsperlevel)
{
   HQueue *hq;
   int i;

   hq = calloc(NUMLEVELS, sizeof(HQueue));
   if (hq==NULL)  return(NULL);
   hq->Pixels = calloc(imgsize, sizeof(ulong));
   if (hq->Pixels==NULL)
   {
      free(hq);
      return(NULL);
   }
   hq->Head = hq->Tail = 0;
   for (i=1; i<NUMLEVELS; i++)
   {
      hq[i].Pixels = hq[i-1].Pixels + numpixelsperlevel[i-1];
      hq[i].Head = hq[i].Tail = 0;
   }
   return(hq);
} /* HQueueCreate */



void HQueueDelete(HQueue *hq)
{
   free(hq->Pixels);
   free(hq);
} /* HQueueDelete */



#define HQueueFirst(hq,h)  (hq[h].Pixels[hq[h].Head++])
#define HQueueAdd(hq,h,p)  hq[h].Pixels[hq[h].Tail++] = p
#define HQueueNotEmpty(hq,h)  (hq[h].Head != hq[h].Tail)



int GetNeighbors(ubyte *shape, ulong imgwidth, ulong imgsize, ulong p,
                 ulong *neighbors)
{
   ulong x;
   int n=0;

   x = p % imgwidth;
   if ((x<(imgwidth-1)) && (shape[p+1]))      neighbors[n++] = p+1;
   if ((p>=imgwidth) && (shape[p-imgwidth]))  neighbors[n++] = p-imgwidth;
   if ((x>0) && (shape[p-1]))                 neighbors[n++] = p-1;
   p += imgwidth;
   if ((p<imgsize) && (shape[p]))             neighbors[n++] = p;
   return(n);
} /* GetNeighbors */



int MaxTreeFlood(MaxTree *mt, HQueue *hq, ulong *numpixelsperlevel,
                 bool *nodeatlevel, ImageGray *img, ubyte *shape, int h,
                 ulong *thisarea,
                 void **thisattr)
/* Returns value >=NUMLEVELS if error */
{
   ulong neighbors[CONNECTIVITY];
   ubyte *pixmap;
   void *attr = NULL, *childattr;
   ulong imgwidth, imgsize, p, q, idx, x, y;
   ulong area = *thisarea, childarea;
   MaxNode *node;
   int numneighbors, i;
   int m;

   imgwidth = img->Width;
   imgsize = imgwidth * (img->Height);
   pixmap = img->Pixmap;
   while(HQueueNotEmpty(hq, h))
   {
      area++;
      p = HQueueFirst(hq, h);
      numneighbors = GetNeighbors(shape, imgwidth, imgsize, p, neighbors);
      x = p % imgwidth;
      y = p / imgwidth;
      if (attr)  mt->AddToAuxData(attr, x, y, numneighbors, neighbors, img);
      else
      {
         attr = mt->NewAuxData(x, y, numneighbors, neighbors, img);
         if (attr==NULL)  return(NUMLEVELS);
         if (*thisattr)  mt->MergeAuxData(attr, *thisattr);
      }
      mt->Status[p] = mt->NumNodesAtLevel[h];
      for (i=0; i<numneighbors; i++)
      {
         q = neighbors[i];
         if (mt->Status[q]==ST_NotAnalyzed)
         {
            HQueueAdd(hq, pixmap[q], q);
            mt->Status[q] = ST_InTheQueue;
            nodeatlevel[pixmap[q]] = true;
            if (pixmap[q] > pixmap[p])
            {
               m = pixmap[q];
               childarea = 0;
               childattr = NULL;
               do
               {
                  m = MaxTreeFlood(mt,hq,numpixelsperlevel,nodeatlevel,img,shape,m, &childarea, &childattr);
                  if (m>=NUMLEVELS)
                  {
                     mt->DeleteAuxData(attr);
                     return(m);
                  }
               } while (m!=h);
               area += childarea;
               mt->MergeAuxData(attr, childattr);
            }
         }
      }
   }
   mt->NumNodesAtLevel[h] = mt->NumNodesAtLevel[h]+1;
   m = h-1;
   while ((m>=0) && (nodeatlevel[m]==false))  m--;
   if (m>=0)
   {
      node = mt->Nodes + (mt->NumPixelsBelowLevel[h] + mt->NumNodesAtLevel[h]-1);
      node->Parent = mt->NumPixelsBelowLevel[m] + mt->NumNodesAtLevel[m];
   } else {
      idx = mt->NumPixelsBelowLevel[h];
      node = mt->Nodes + idx;
      node->Parent = idx;
   }
   node->Area = area;
   node->Attribute = attr;
   node->Level = h;
   nodeatlevel[h] = false;
   *thisarea = area;
   *thisattr = attr;
   return(m);
} /* MaxTreeFlood */



MaxTree *MaxTreeCreate(ImageGray *img, ImageGray *template,
                       void *(*newauxdata)(ulong, ulong, int, ulong *, ImageGray *),
                       void (*addtoauxdata)(void *, ulong, ulong, int, ulong *, ImageGray *),
                       void (*mergeauxdata)(void *, void *),
                       void (*deleteauxdata)(void *))
{
   ulong numpixelsperlevel[NUMLEVELS];
   bool nodeatlevel[NUMLEVELS];
   HQueue *hq;
   MaxTree *mt;
   ubyte *pixmap = img->Pixmap;
   void *attr = NULL;
   ulong imgsize, p, m=0, area=0;
   int l;

   /* Allocate structures */
   mt = malloc(sizeof(MaxTree));
   if (mt==NULL)  return(NULL);
   imgsize = (img->Width)*(img->Height);
   mt->Status = calloc((size_t)imgsize, sizeof(long));
   if (mt->Status==NULL)
   {
      free(mt);
      return(NULL);
   }
   mt->NumPixelsBelowLevel = calloc(NUMLEVELS, sizeof(ulong));
   if (mt->NumPixelsBelowLevel==NULL)
   {
      free(mt->Status);
      free(mt);
      return(NULL);
   }
   mt->NumNodesAtLevel = calloc(NUMLEVELS, sizeof(ulong));
   if (mt->NumNodesAtLevel==NULL)
   {
      free(mt->NumPixelsBelowLevel);
      free(mt->Status);
      free(mt);
      return(NULL);
   }
   mt->Nodes = calloc((size_t)imgsize, sizeof(MaxNode));
   if (mt->Nodes==NULL)
   {
      free(mt->NumNodesAtLevel);
      free(mt->NumPixelsBelowLevel);
      free(mt->Status);
      free(mt);
      return(NULL);
   }

   /* Initialize structures */
   for (p=0; p<imgsize; p++)  mt->Status[p] = ST_NotAnalyzed;
   bzero(nodeatlevel, NUMLEVELS*sizeof(bool));
   bzero(numpixelsperlevel, NUMLEVELS*sizeof(ulong));
   /* Following bzero is redundant, array is initialized by calloc */
   /* bzero(mt->NumNodesAtLevel, NUMLEVELS*sizeof(ulong)); */
   for (p=0; p<imgsize; p++)  numpixelsperlevel[pixmap[p]]++;
   mt->NumPixelsBelowLevel[0] = 0;
   for (l=1; l<NUMLEVELS; l++)
   {
      mt->NumPixelsBelowLevel[l] = mt->NumPixelsBelowLevel[l-1] + numpixelsperlevel[l-1];
   }
   hq = HQueueCreate(imgsize, numpixelsperlevel);
   if (hq==NULL)
   {
      free(mt->Nodes);
      free(mt->NumNodesAtLevel);
      free(mt->NumPixelsBelowLevel);
      free(mt->Status);
      free(mt);
      return(NULL);
   }

   /* Find pixel m which has the lowest intensity l in the image */
   for (p=0; p<imgsize; p++)
   {
      if (pixmap[p]<pixmap[m])  m = p;
   }
   l = pixmap[m];

   /* Add pixel m to the queue */
   nodeatlevel[l] = true;
   HQueueAdd(hq, l, m);
   mt->Status[m] = ST_InTheQueue;

   /* Build the Max-tree using a flood-fill algorithm */
   mt->NewAuxData = newauxdata;
   mt->AddToAuxData = addtoauxdata;
   mt->MergeAuxData = mergeauxdata;
   mt->DeleteAuxData = deleteauxdata;
   l = MaxTreeFlood(mt, hq, numpixelsperlevel, nodeatlevel, img,
                    template->Pixmap, l, &area, &attr);
   if (l>=NUMLEVELS)  MaxTreeDelete(mt);
   HQueueDelete(hq);
   return(mt);
} /* MaxTreeCreate */



void MaxTreeDelete(MaxTree *mt)
{
   void *attr;
   ulong i;
   int h;

   for (h=0; h<NUMLEVELS; h++)
   {
      for (i=0; i<mt->NumNodesAtLevel[h]; i++)
      {
         attr = mt->Nodes[mt->NumPixelsBelowLevel[h]+i].Attribute;
         if (attr)  mt->DeleteAuxData(attr);
      }
   }
   free(mt->Nodes);
   free(mt->NumNodesAtLevel);
   free(mt->NumPixelsBelowLevel);
   free(mt->Status);
   free(mt);
} /* MaxTreeDelete */



void MaxTreeFilterMin(MaxTree *mt, ImageGray *img, ImageGray *template,
                      ImageGray *out, double (*attribute)(void *),
                      double lambda)
{
   MaxNode *node, *parnode;
   ubyte *shape = template->Pixmap;
   ulong i, idx, parent;
   int l;

   for (l=0; l<NUMLEVELS; l++)
   {
      for (i=0; i<mt->NumNodesAtLevel[l]; i++)
      {
         idx = mt->NumPixelsBelowLevel[l] + i;
         node = &(mt->Nodes[idx]);
         parent = node->Parent;
         if (idx!=parent)
         {
            parnode = &(mt->Nodes[parent]);
            if (((*attribute)(node->Attribute) < lambda) || (parnode->Level!=parnode->NewLevel))
            {
               node->NewLevel = parnode->NewLevel;
            } else  node->NewLevel = node->Level;
         }
      }
   }
   for (i=0; i<(img->Width)*(img->Height); i++)
   {
      if (shape[i])
      {
         idx = mt->NumPixelsBelowLevel[img->Pixmap[i]] + mt->Status[i];
         out->Pixmap[i] = mt->Nodes[idx].NewLevel;
      }
   }
} /* MaxTreeFilterMin */



void MaxTreeFilterDirect(MaxTree *mt, ImageGray *img, ImageGray *template,
                         ImageGray *out, double (*attribute)(void *),
                         double lambda)
{
   MaxNode *node;
   ubyte *shape = template->Pixmap;
   ulong i, idx, parent;
   int l;

   for (l=0; l<NUMLEVELS; l++)
   {
      for (i=0; i<mt->NumNodesAtLevel[l]; i++)
      {
         idx = mt->NumPixelsBelowLevel[l] + i;
         node = &(mt->Nodes[idx]);
         parent = node->Parent;
         if (idx!=parent)
         {
            if ((*attribute)(node->Attribute) < lambda)  node->NewLevel = mt->Nodes[parent].NewLevel;
            else  node->NewLevel = node->Level;
         }
      }
   }
   for (i=0; i<(img->Width)*(img->Height); i++)
   {
      if (shape[i])
      {
         idx = mt->NumPixelsBelowLevel[img->Pixmap[i]] + mt->Status[i];
         out->Pixmap[i] = mt->Nodes[idx].NewLevel;
      }
   }
} /* MaxTreeFilterDirect */



void MaxTreeFilterMax(MaxTree *mt, ImageGray *img, ImageGray *template,
                      ImageGray *out, double (*attribute)(void *),
                      double lambda)
{
   MaxNode *node;
   ubyte *shape = template->Pixmap;
   ulong i, idx, parent;
   int l;

   for (l=0; l<NUMLEVELS; l++)
   {
      for (i=0; i<mt->NumNodesAtLevel[l]; i++)
      {
         idx = mt->NumPixelsBelowLevel[l] + i;
         node = &(mt->Nodes[idx]);
         parent = node->Parent;
         if (idx!=parent)
         {
            if ((*attribute)(node->Attribute) < lambda)  node->NewLevel = mt->Nodes[parent].NewLevel;
            else  node->NewLevel = node->Level;
         }
      }
   }
   for (l=NUMLEVELS-1; l>0; l--)
   {
      for (i=0; i<mt->NumNodesAtLevel[l]; i++)
      {
         idx = mt->NumPixelsBelowLevel[l] + i;
         node = &(mt->Nodes[idx]);
         parent = node->Parent;
         if ((idx!=parent) && (node->NewLevel==node->Level))
         {
            mt->Nodes[parent].NewLevel = mt->Nodes[parent].Level;
         }
      }
   }
   for (i=0; i<(img->Width)*(img->Height); i++)
   {
      if (shape[i])
      {
         idx = mt->NumPixelsBelowLevel[img->Pixmap[i]] + mt->Status[i];
         out->Pixmap[i] = mt->Nodes[idx].NewLevel;
      }
   }
} /* MaxTreeFilterMax */



void MaxTreeFilterSubtractive(MaxTree *mt, ImageGray *img, ImageGray *template,
                              ImageGray *out, double (*attribute)(void *),
                              double lambda)
{
   MaxNode *node, *parnode;
   ubyte *shape = template->Pixmap;
   ulong i, idx, parent;
   int l;

   for (l=0; l<NUMLEVELS; l++)
   {
      for (i=0; i<mt->NumNodesAtLevel[l]; i++)
      {
         idx = mt->NumPixelsBelowLevel[l] + i;
         node = &(mt->Nodes[idx]);
         parent = node->Parent;
         if (idx!=parent)
         {
            parnode = &(mt->Nodes[parent]);
            if ((*attribute)(node->Attribute) < lambda)  node->NewLevel = parnode->NewLevel;
            else  node->NewLevel = ((int)(node->Level)) + ((int)(parnode->NewLevel)) - ((int)(parnode->Level));
         }
      }
   }
   for (i=0; i<(img->Width)*(img->Height); i++)
   {
      if (shape[i])
      {
         idx = mt->NumPixelsBelowLevel[img->Pixmap[i]] + mt->Status[i];
         out->Pixmap[i] = mt->Nodes[idx].NewLevel;
      }
   }
} /* MaxTreeFilterSubtractive */



ImageGray *GetTemplate(char *templatefname, ImageGray *img)
{
   ImageGray *template;

   if (templatefname)
   {
      template = ImagePGMRead(templatefname);
      if (template==NULL)  return(NULL);
      if ((img->Width != template->Width) || (img->Height != template->Height))
      {
	 ImageGrayDelete(template);
         return(NULL);
      }
   } else {
      template = ImageGrayCreate(img->Width, img->Height);
      if (template)  ImageGrayInit(template, NUMLEVELS-1);
   }
   return(template);
} /* GetTemplate */


typedef struct AttribStruct AttribStruct;

struct AttribStruct
{
   char *Name;
   void *(*NewAuxData)(ulong, ulong, int, ulong *, ImageGray *);
   void (*DeleteAuxData)(void *);
   void (*AddToAuxData)(void *, ulong, ulong, int, ulong *, ImageGray *);
   void (*MergeAuxData)(void *, void *);
   double (*Attribute)(void *);
};

#define NUMATTR 19

AttribStruct Attribs[NUMATTR] =
{
  {"Area", NewAreaData, DeleteAreaData, AddToAreaData, MergeAreaData, AreaAttribute},
  {"Area of min. enclosing rectangle", NewEnclRectData, DeleteEnclRectData, AddToEnclRectData, MergeEnclRectData, EnclRectAreaAttribute},
  {"Square of diagonal of min. enclosing rectangle", NewEnclRectData, DeleteEnclRectData, AddToEnclRectData, MergeEnclRectData, EnclRectDiagAttribute},
  {"Cityblock perimeter", NewPeriCBData, DeletePeriCBData, AddToPeriCBData, MergePeriCBData, PeriCBPerimeterAttribute},
  {"Cityblock complexity (Perimeter/Area)", NewPeriCBData, DeletePeriCBData, AddToPeriCBData, MergePeriCBData, PeriCBComplexityAttribute},
  {"Cityblock simplicity (Area/Perimeter)", NewPeriCBData, DeletePeriCBData, AddToPeriCBData, MergePeriCBData, PeriCBSimplicityAttribute},
  {"Cityblock compactness (Perimeter^2/(4*PI*Area))", NewPeriCBData, DeletePeriCBData, AddToPeriCBData, MergePeriCBData, PeriCBCompactnessAttribute},
  {"Large perimeter", NewPeriLargeData, DeletePeriLargeData, AddToPeriLargeData, MergePeriLargeData, PeriLargePerimeterAttribute},
  {"Large compactness (Perimeter^2/(4*PI*Area))", NewPeriLargeData, DeletePeriLargeData, AddToPeriLargeData, MergePeriLargeData, PeriLargeCompactnessAttribute},
  {"Small perimeter", NewPeriSmallData, DeletePeriSmallData, AddToPeriSmallData, MergePeriSmallData, PeriSmallPerimeterAttribute},
  {"Small compactness (Perimeter^2/(4*PI*Area))", NewPeriSmallData, DeletePeriSmallData, AddToPeriSmallData, MergePeriSmallData, PeriSmallCompactnessAttribute},
  {"Moment of Inertia", NewInertiaData, DeleteInertiaData, AddToInertiaData, MergeInertiaData, InertiaAttribute},
  {"Elongation: (Moment of Inertia) / (area)^2", NewInertiaData, DeleteInertiaData, AddToInertiaData, MergeInertiaData, InertiaDivA2Attribute},
  {"Mean X position", NewInertiaData, DeleteInertiaData, AddToInertiaData, MergeInertiaData, MeanXAttribute},
  {"Mean Y position", NewInertiaData, DeleteInertiaData, AddToInertiaData, MergeInertiaData, MeanYAttribute},
  {"Jaggedness: Area*Perimeter^2/(8*PI^2*Inertia)", NewJaggedData, DeleteJaggedData, AddToJaggedData, MergeJaggedData, JaggednessAttribute},
  {"Entropy", NewEntropyData, DeleteEntropyData, AddToEntropyData, MergeEntropyData, EntropyAttribute},
  {"Lambda-max (Max.child gray level - current gray level)", NewLambdamaxData, DeleteLambdamaxData, AddToLambdamaxData, MergeLambdamaxData, LambdamaxAttribute},
  {"Gray level", NewLevelData, DeleteLevelData, AddToLevelData, MergeLevelData, LevelAttribute}
};


typedef struct DecisionStruct DecisionStruct;

struct DecisionStruct
{
   char *Name;
   void (*Filter)(MaxTree *, ImageGray *, ImageGray *, ImageGray *, double (*attribute)(void *), double);
};

#define NUMDECISIONS 4

DecisionStruct Decisions[NUMDECISIONS] =
{
   {"Min", MaxTreeFilterMin},
   {"Direct", MaxTreeFilterDirect},
   {"Max", MaxTreeFilterMax},
   {"Subtractive", MaxTreeFilterSubtractive},
};



int main(int argc, char *argv[])
{
   ImageGray *img, *template, *out;
   MaxTree *mt;
   char *imgfname, *templatefname = NULL, *outfname = "out.pgm";
   double lambda;
   int attrib, decision=3, r;
   long tickspersec = sysconf(_SC_CLK_TCK);  
   float musec;

   if (argc<4)
   {
      printf("Usage: %s <input image> <attrib> <lambda> [decision] [output image] [template]\n", argv[0]);
      printf("Where attrib is:\n");
      for (attrib=0; attrib<NUMATTR; attrib++)
      {
         printf("\t%d - %s\n", attrib, Attribs[attrib].Name);
      }
      printf("and decision is:\n");
      for (r=0; r<NUMDECISIONS; r++)
      {
         printf("\t%d - %s", r, Decisions[r].Name);
	 if (r==decision)  printf(" (default)");
	 printf("\n");
      }
      exit(0);
   }
   imgfname = argv[1];
   attrib = atoi(argv[2]);
   lambda = atof(argv[3]);
   if (argc>=5)  decision = atoi(argv[4]);
   if (argc>=6)  outfname = argv[5];
   if (argc>=7)  templatefname = argv[6];
   img = ImagePGMRead(imgfname);
   if (img==NULL)
   {
      fprintf(stderr, "Can't read image '%s'\n", imgfname);
      return(-1);
   }
   template = GetTemplate(templatefname, img);
   if (template==NULL)
   {
      fprintf(stderr, "Can't create template\n");
      ImageGrayDelete(img);
      return(-1);
   }
   printf("Filtering image '%s' using attribute '%s' with lambda=%f\n", imgfname, Attribs[attrib].Name, lambda);
   printf("Decision rule: %s   Template: ", Decisions[decision].Name);
   if (templatefname==NULL)  printf("<not used>\n");
   else  printf("%s\n", templatefname);
   printf("Image: Width=%ld Height=%ld\n", img->Width, img->Height);
   out = ImageGrayCreate(img->Width, img->Height);
   if (out==NULL)
   {
      fprintf(stderr, "Can't create output image\n");
      ImageGrayDelete(template);
      ImageGrayDelete(img);
      return(-1);
   }
   start = times(&tstruct);
   mt = MaxTreeCreate(img, template, Attribs[attrib].NewAuxData, Attribs[attrib].AddToAuxData, Attribs[attrib].MergeAuxData, Attribs[attrib].DeleteAuxData);
   if (mt==NULL)
   {
      fprintf(stderr, "Can't create Max-tree\n");
      ImageGrayDelete(out);
      ImageGrayDelete(template);
      ImageGrayDelete(img);
      return(-1);
   }
   Decisions[decision].Filter(mt, img, template, out, Attribs[attrib].Attribute, lambda);
   finish = times(&tstruct);
   musec = (float)(finish - start)/((float)tickspersec);
   MaxTreeDelete(mt);
   printf("wall-clock time: %f s\n",musec);
   r = ImagePGMBinWrite(out, outfname);
   if (r)  fprintf(stderr, "Error writing image '%s'\n", outfname);
   else  printf("Filtered image written to '%s'\n", outfname);
   ImageGrayDelete(out);
   ImageGrayDelete(template);
   ImageGrayDelete(img);
   return(0);
} /* main */
