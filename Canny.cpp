/* source: http://marathon.csee.usf.edu/edge/edge_detection.html */
/* URL: ftp://figment.csee.usf.edu/pub/Edge_Comparison/source_code/canny.src */

/* CENG 381 Assignment 6 solution */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "systemc.h"

#define VERBOSE 1

#define NOEDGE 255
#define POSSIBLE_EDGE 128
#define EDGE 0
#define BOOSTBLURFACTOR 90.0
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#define SIGMA 0.6
#define TLOW  0.3
#define THIGH 0.8

#define COLS 3840
#define ROWS 2160
#define SIZE COLS*ROWS
#define VIDEONAME "Keck"
#define IMG_IN    "video/" VIDEONAME "%03d.pgm"
#define IMG_OUT   VIDEONAME "%03d_edges.pgm"
#define IMG_NUM   30 /* number of images processed (1 or more) */
#define AVAIL_IMG 30 /* number of different image frames (1 or more) */
#define SET_STACK_SIZE set_stack_size(128*1024*1024);
/* upper bound for the size of the gaussian kernel
 * SIGMA must be less than 4.0
 * check for 'windowsize' below
 */
#define WINSIZE 21

typedef struct Image_s
{
	unsigned char img[SIZE];

	Image_s(void)
	{
	   for (int i=0; i<SIZE; i++)
	   {
	      img[i] = 0;
	   }
	}

	Image_s& operator=(const Image_s& copy)
	{
	   for (int i=0; i<SIZE; i++)
	   {
	      img[i] = copy.img[i];
	   }
	   return *this;
	}

	operator unsigned char*()
	{
	   return img;
	}

	unsigned char& operator[](const int index)
	{
	   return img[index];
	}
} IMAGE;

typedef struct Image_i
{
        short int img[SIZE];

        Image_i(void)
        {
           for (int i=0; i<SIZE; i++)
           {
              img[i] = 0;
           }
	}

	Image_i& operator=(const Image_i& copy)
        {
           for (int i=0; i<SIZE; i++)
           {
              img[i] = copy.img[i];
           }
           return *this;
        }

	operator short int*()
        {
           return img;
        }

	short int& operator[](const int index)
        {
           return img[index];
        }
} SIMAGE;

typedef struct kernel_i
{
        float img[WINSIZE];

        kernel_i(void)
        {
           for (int i=0; i<WINSIZE; i++)
           {
              img[i] = 0;
           }
	}

        kernel_i& operator=(const kernel_i& copy)
        {
           for (int i=0; i<WINSIZE; i++)
           {
              img[i] = copy.img[i];
           }
           return *this;
        }

	operator float*()
        {
           return img;
        }

        float& operator[](const int index)
        {
           return img[index];
        }
} FLOAT;

typedef struct kernel_s
{
        float img[SIZE];

        kernel_s(void)
        {
           for (int i=0; i<SIZE; i++)
           {
              img[i] = 0;
           }
        }

	kernel_s& operator=(const kernel_s& copy)
        {
           for (int i=0; i<SIZE; i++)
           {
              img[i] = copy.img[i];
           }
           return *this;
        }

	operator float*()
        {
           return img;
        }

	float& operator[](const int index)
        {
           return img[index];
        }
} GOAT;


SC_MODULE(Stimulus){
	IMAGE imageout;
	sc_fifo_out<IMAGE> ImgOut;
	/******************************************************************************
	* Function: read_pgm_image
	* Purpose: This function reads in an image in PGM format. The image can be
	* read in from either a file or from standard input. The image is only read
	* from standard input when infilename = NULL. Because the PGM format includes
	* the number of columns and the number of rows in the image, these are read
	* from the file. Memory to store the image is allocated OUTSIDE this function.
	* The found image size is checked against the expected rows and cols.
	* All comments in the header are discarded in the process of reading the
	* image. Upon failure, this function returns 0, upon sucess it returns 1.
	******************************************************************************/
	int read_pgm_image(const char *infilename, unsigned char *image, int rows, int cols)
	{
	   FILE *fp;
	   char buf[71];
	   int r, c;

	   /***************************************************************************
	   * Open the input image file for reading if a filename was given. If no
	   * filename was provided, set fp to read from standard input.
	   ***************************************************************************/
	   if(infilename == NULL) fp = stdin;
	   else{
	      if((fp = fopen(infilename, "r")) == NULL){
	         fprintf(stderr, "Error reading the file %s in read_pgm_image().\n",
	            infilename);
	         return(0);
	      }
	   }

	   /***************************************************************************
	   * Verify that the image is in PGM format, read in the number of columns
	   * and rows in the image and scan past all of the header information.
	   ***************************************************************************/
	   fgets(buf, 70, fp);
	   if(strncmp(buf,"P5",2) != 0){
	      fprintf(stderr, "The file %s is not in PGM format in ", infilename);
	      fprintf(stderr, "read_pgm_image().\n");
	      if(fp != stdin) fclose(fp);
	      return(0);
	   }
	   do{ fgets(buf, 70, fp); }while(buf[0] == '#');  /* skip all comment lines */
	   sscanf(buf, "%d %d", &c, &r);
	   if(c != cols || r != rows){
	      fprintf(stderr, "The file %s is not a %d by %d image in ", infilename, cols, rows);
	      fprintf(stderr, "read_pgm_image().\n");
	      if(fp != stdin) fclose(fp);
	      return(0);
	   }
	   do{ fgets(buf, 70, fp); }while(buf[0] == '#');  /* skip all comment lines */

	   /***************************************************************************
	   * Read the image from the file.
	   ***************************************************************************/
	   if((unsigned)rows != fread(image, cols, rows, fp)){
	      fprintf(stderr, "Error reading the image data in read_pgm_image().\n");
	      if(fp != stdin) fclose(fp);
	      return(0);
	   }

	   if(fp != stdin) fclose(fp);
	   return(1);
	}

	void main(void)
	{
	   int i=0, n=0;
	   char infilename[40];

	   for(i=0; i<IMG_NUM; i++)
	   {
	      n = i % AVAIL_IMG;
	      sprintf(infilename, IMG_IN, n+1);

	      /****************************************************************************
	      * Read in the image.
	      ****************************************************************************/
	      if(VERBOSE) printf("Reading the image %s.\n", infilename);
	      if(read_pgm_image(infilename, imageout, ROWS, COLS) == 0){
	         fprintf(stderr, "Error reading the input image, %s.\n", infilename);
	         exit(1);
	      }
	      ImgOut.write(imageout);
	   }
	}

	SC_CTOR(Stimulus)
	{
	   SC_THREAD(main);
	   SET_STACK_SIZE
	}
};

SC_MODULE(Monitor){
	IMAGE imagein;
	sc_fifo_in<IMAGE>  ImgIn;

	/******************************************************************************
	* Function: write_pgm_image
	* Purpose: This function writes an image in PGM format. The file is either
	* written to the file specified by outfilename or to standard output if
	* outfilename = NULL. A comment can be written to the header if coment != NULL.
	******************************************************************************/
	int write_pgm_image(const char *outfilename, unsigned char *image, int rows,
	    int cols, const char *comment, int maxval)
	{
	   FILE *fp;

	   /***************************************************************************
	   * Open the output image file for writing if a filename was given. If no
	   * filename was provided, set fp to write to standard output.
	   ***************************************************************************/
	   if(outfilename == NULL) fp = stdout;
	   else{
	      if((fp = fopen(outfilename, "w")) == NULL){
	         fprintf(stderr, "Error writing the file %s in write_pgm_image().\n",
	            outfilename);
	         return(0);
	      }
	   }

	   /***************************************************************************
	   * Write the header information to the PGM file.
	   ***************************************************************************/
	   fprintf(fp, "P5\n%d %d\n", cols, rows);
	   if(comment != NULL)
	      if(strlen(comment) <= 70) fprintf(fp, "# %s\n", comment);
	   fprintf(fp, "%d\n", maxval);

	   /***************************************************************************
	   * Write the image data to the file.
	   ***************************************************************************/
	   if((unsigned)rows != fwrite(image, cols, rows, fp)){
	      fprintf(stderr, "Error writing the image data in write_pgm_image().\n");
	      if(fp != stdout) fclose(fp);
	      return(0);
	   }

	   if(fp != stdout) fclose(fp);
	   return(1);
	}

	void main(void)
	{
	   char outfilename[128];    /* Name of the output "edge" image */
	   int i, n;

	   for(i=0; i<IMG_NUM; i++)
	   {
	      ImgIn.read(imagein);

	      /****************************************************************************
	      * Write out the edge image to a file.
	      ****************************************************************************/
	      n = i % AVAIL_IMG;
	      sprintf(outfilename, IMG_OUT, n+1);
	      if(VERBOSE) printf("Writing the edge image in the file %s.\n", outfilename);
	      if(write_pgm_image(outfilename, imagein, ROWS, COLS,"", 255) == 0){
	         fprintf(stderr, "Error writing the edge image, %s.\n", outfilename);
	         exit(1);
	      }
	   }
	   if(VERBOSE) printf("Monitor exits simulation.\n");
	      sc_stop();	// done testing, quit the simulation
	}

	SC_CTOR(Monitor)
	{
	   SC_THREAD(main);
	   SET_STACK_SIZE
	}
};

SC_MODULE(DataIn)
{
	sc_fifo_in<IMAGE> ImgIn;
	sc_fifo_out<IMAGE> ImgOut;
	IMAGE Image;

	void main()
	{
	   while(1)
	   {
	      ImgIn.read(Image);
	      ImgOut.write(Image);
	   }
	}

	SC_CTOR(DataIn)
	{
	   SC_THREAD(main);
	   SET_STACK_SIZE
	}
};

SC_MODULE(DataOut)
{
	sc_fifo_in<IMAGE> ImgIn;
	sc_fifo_out<IMAGE> ImgOut;
	IMAGE Image;

	void main()
	{
	   while(1)
	   {
	      ImgIn.read(Image);
	      ImgOut.write(Image);
	   }
	}

	SC_CTOR(DataOut)
	{
	   SC_THREAD(main);
	   SET_STACK_SIZE
	}
};

SC_MODULE(DUT)
{
	sc_fifo_in<IMAGE> ImgIn;
	sc_fifo_out<IMAGE> ImgOut;
	sc_fifo<SIMAGE> q1;
        sc_fifo<SIMAGE> q2;
	sc_fifo<SIMAGE> q3;
        sc_fifo<SIMAGE> q4;
	sc_fifo<SIMAGE> q5;
        sc_fifo<SIMAGE> q6;
	sc_fifo<SIMAGE> q7;
        sc_fifo<IMAGE> q8;
	IMAGE imageout;
	IMAGE imagein;

	SC_MODULE(mgk){

		sc_fifo_in<IMAGE> ImgIn;
                sc_fifo_out<SIMAGE> ImgOut;
		sc_fifo<FLOAT> kern1;
		sc_fifo<FLOAT> kern2;
		sc_fifo<int> c1;
		sc_fifo<int> c2;
		sc_fifo<GOAT> temporary;
		SIMAGE smoothdim;
		IMAGE imagein;

		SC_MODULE(Gauss){
			sc_fifo_out<FLOAT> k1;
			sc_fifo_out<FLOAT> k2;
			sc_fifo_out<int> center1;
			sc_fifo_out<int> center2;

			FLOAT kernel;
			int center;

			void make_gaussian_kernel()
		        {
           			int i, windowsize;
           			float x, fx, sum=0.0;

           			windowsize = 1 + 2 * ceil(2.5 * SIGMA);
           			center = (windowsize) / 2;

           			if(VERBOSE) printf("      The kernel has %d elements.\n", windowsize);

           			for(i=0;i<(windowsize);i++){
              				x = (float)(i - center);
              				fx = pow(2.71828, -0.5*x*x/(SIGMA*SIGMA)) / (SIGMA * sqrt(6.2831853));
              				kernel[i] = fx;
              				sum += fx;
           			}

           			for(i=0;i<(windowsize);i++) kernel[i] /= sum;

           			if(VERBOSE){
              				printf("The filter coefficients are:\n");
              				for(i=0;i<(windowsize);i++)
                 			printf("kernel[%d] = %f\n", i, kernel[i]);
           			}
			}

			void action(){
				while(1){
					make_gaussian_kernel();
					k1.write(kernel);
					k2.write(kernel);
					center1.write(center);
					center2.write(center);
				}


			}

			SC_CTOR(Gauss)
        		{
           			SC_THREAD(action);
           			SET_STACK_SIZE
        		}


		};

		SC_MODULE(Blurx){
			sc_fifo_in<IMAGE> ImgIn;
			sc_fifo_in<FLOAT> k1;
			sc_fifo_in<int> c1;
                        sc_fifo_out<GOAT> t;
                        IMAGE imagein;
			FLOAT kernel;
			int center;
			GOAT tempim;

			void action(){
				float dot;
				float sum;
				while(1){

					k1.read(kernel);
					c1.read(center);
					ImgIn.read(imagein);

					if(VERBOSE) printf("   Bluring the image in the X-direction.\n");
           				for(int r=0;r<ROWS;r++){
              					for(int c=0;c<COLS;c++){
                 					dot = 0.0;
                 					sum = 0.0;
                 					for(int cc=(-center);cc<=center;cc++){
                    						if(((c+cc) >= 0) && ((c+cc) < COLS)){
                       							dot += (float)imagein[r*COLS+(c+cc)] * kernel[center+cc];
                       							sum += kernel[center+cc];
                    						}
                 					}
                 					tempim[r*COLS+c] = dot/sum;
              					}
           				}

					t.write(tempim);

				}
			}

			SC_CTOR(Blurx){
				SC_THREAD(action);
				SET_STACK_SIZE
			}
		};

		SC_MODULE(Blury){

			sc_fifo_in<FLOAT> k2;
			sc_fifo_in<int> c2;
			sc_fifo_in<GOAT> tempin;
                        sc_fifo_out<SIMAGE> ImgOut;
			FLOAT kernel;
			int center;
			GOAT tempim;
			SIMAGE smoothedim;

                        void action(){
				float sum;
				float dot;
	                        while(1){

					k2.read(kernel);
					c2.read(center);
					tempin.read(tempim);

					if(VERBOSE) printf("   Bluring the image in the Y-direction.\n");
           				for(int c=0;c<COLS;c++){
              					for(int r=0;r<ROWS;r++){
                 					sum = 0.0;
                 					dot = 0.0;
                 					for(int rr=(-center);rr<=center;rr++){
                    						if(((r+rr) >= 0) && ((r+rr) < ROWS)){
                       							dot += tempim[(r+rr)*COLS+c] * kernel[center+rr];
                       							sum += kernel[center+rr];
                    						}
                 					}
                 					smoothedim[r*COLS+c] = (short int)(dot*BOOSTBLURFACTOR/sum + 0.5);
              					}
           				}

					ImgOut.write(smoothedim);

				}
			}
			SC_CTOR(Blury){
                                SC_THREAD(action);
                                SET_STACK_SIZE
                        }

		};

		Gauss gauss;
		Blurx blurx;
		Blury blury;

		void before_end_of_elaboration(){

			blurx.ImgIn.bind(ImgIn);
			gauss.k1.bind(kern1);
			gauss.center1.bind(c1);
			blurx.k1.bind(kern1);
			blurx.c1.bind(c1);
			gauss.k2.bind(kern2);
			gauss.center2.bind(c2);
			blurx.t.bind(temporary);
			blury.k2.bind(kern2);
			blury.c2.bind(c2);
			blury.tempin.bind(temporary);
			blury.ImgOut.bind(ImgOut);



		}

		SC_CTOR(mgk)
		        :kern1("kern1")
			,kern2("kern2")
			,c1("c1")
			,c2("c2")
			,temporary("temporary")
			,gauss("gauss")
			,blurx("blurx")
			,blury("blury")
	        {}



	};


		SC_MODULE(derivx_y){

			sc_fifo_in<SIMAGE> ImgIn;
                        sc_fifo_out<SIMAGE> deltax_to_mag;
			sc_fifo_out<SIMAGE> deltay_to_mag;
			sc_fifo_out<SIMAGE> deltax_to_non;
			sc_fifo_out<SIMAGE> deltay_to_non;
			SIMAGE smoothedim;
			SIMAGE delta_x;
			SIMAGE delta_y;

			void derivative_x_y(int rows, int cols)
        		{
           			int r, c, pos;

           			/****************************************************************************
           			* Compute the x-derivative. Adjust the derivative at the borders to avoid
           			* losing pixels.
           			****************************************************************************/
           			if(VERBOSE) printf("   Computing the X-direction derivative.\n");
           			for(r=0;r<rows;r++){
              				pos = r * cols;
              				delta_x[pos] = smoothedim[pos+1] - smoothedim[pos];
              				pos++;
              				for(c=1;c<(cols-1);c++,pos++){
                 				delta_x[pos] = smoothedim[pos+1] - smoothedim[pos-1];
              				}
              				delta_x[pos] = smoothedim[pos] - smoothedim[pos-1];
           			}

           			/****************************************************************************
           			* Compute the y-derivative. Adjust the derivative at the borders to avoid
           			* losing pixels.
           			****************************************************************************/
           			if(VERBOSE) printf("   Computing the Y-direction derivative.\n");
           			for(c=0;c<cols;c++){
              				pos = c;
              				delta_y[pos] = smoothedim[pos+cols] - smoothedim[pos];
              				pos += cols;
              				for(r=1;r<(rows-1);r++,pos+=cols){
                 				delta_y[pos] = smoothedim[pos+cols] - smoothedim[pos-cols];
              				}
             			 		delta_y[pos] = smoothedim[pos] - smoothedim[pos-cols];
	   			}
			}


			void action(){
				while(1){
					ImgIn.read(smoothedim);
					/****************************************************************************
           				* Compute the first derivative in the x and y directions.
           				****************************************************************************/
           				if(VERBOSE) printf("Computing the X and Y first derivatives.\n");
           				derivative_x_y(ROWS, COLS);
					deltax_to_mag.write(delta_x);
					deltay_to_mag.write(delta_y);
					deltax_to_non.write(delta_x);
					deltay_to_non.write(delta_y);
				};
			}

			SC_CTOR(derivx_y){
				SC_THREAD(action);
				SET_STACK_SIZE
			}

		};

		SC_MODULE(magx_y){

			sc_fifo_in<SIMAGE> Delta_x;
			sc_fifo_in<SIMAGE> Delta_y;
			sc_fifo_out<SIMAGE> mag1;
                        sc_fifo_out<SIMAGE> mag2;
                        SIMAGE delta_x;
                        SIMAGE delta_y;
			SIMAGE magnitude;

			void magnitude_x_y(int rows, int cols)
        		{
           			int r, c, pos, sq1, sq2;
           			for(r=0,pos=0;r<rows;r++){
              				for(c=0;c<cols;c++,pos++){
                 				sq1 = (int)delta_x[pos] * (int)delta_x[pos];
                 				sq2 = (int)delta_y[pos] * (int)delta_y[pos];
                 				magnitude[pos] = (short)(0.5 + sqrt((float)sq1 + (float)sq2));
              				}
           			}

			}


			void action(){
				while(1){
					Delta_x.read(delta_x);
					Delta_y.read(delta_y);
					if(VERBOSE) printf("Computing the magnitude of the gradient.\n");
           				magnitude_x_y(ROWS, COLS);
					mag1.write(magnitude);
					mag2.write(magnitude);
				}

			}

			SC_CTOR(magx_y){
				SC_THREAD(action);
				SET_STACK_SIZE
			}

		};


		SC_MODULE(non_max){
			sc_fifo_in<SIMAGE> pmag;
			sc_fifo_in<SIMAGE> pgraddy;
			sc_fifo_in<SIMAGE> pgraddx;
                        sc_fifo_out<IMAGE> pnms;
                        SIMAGE mag;
                        SIMAGE gradx;
			SIMAGE grady;
			IMAGE result;

			/*****************************************************************************************************/
			void non_max_supp(int nrows, int ncols)
        		{
            			int rowcount, colcount,count;
            			short *magrowptr,*magptr;
            			short *gxrowptr,*gxptr;
            			short *gyrowptr,*gyptr,z1,z2;
            			short m00; short gx=0; short gy=0;
            			float mag1, mag2; float xperp=0;float yperp=0;
            			unsigned char *resultrowptr, *resultptr;

           			/****************************************************************************
           			* Zero the edges of the result image.
           			****************************************************************************/
           			for(count=0,resultrowptr=result,resultptr=result+ncols*(nrows-1);
                			count<ncols; resultptr++,resultrowptr++,count++){
                			*resultrowptr = *resultptr = (unsigned char) 0;
            			}

            			for(count=0,resultptr=result,resultrowptr=result+ncols-1;
                			count<nrows; count++,resultptr+=ncols,resultrowptr+=ncols){
                			*resultptr = *resultrowptr = (unsigned char) 0;
            			}

           			/****************************************************************************
           			* Suppress non-maximum points.
           			****************************************************************************/
           			for(rowcount=1,magrowptr=mag+ncols+1,gxrowptr=gradx+ncols+1,
              				gyrowptr=grady+ncols+1,resultrowptr=result+ncols+1;
              				rowcount<=nrows-2;        // bug fix 3/29/17, RD
              				rowcount++,magrowptr+=ncols,gyrowptr+=ncols,gxrowptr+=ncols,
              				resultrowptr+=ncols){
              					for(colcount=1,magptr=magrowptr,gxptr=gxrowptr,gyptr=gyrowptr,
                 				resultptr=resultrowptr;colcount<=ncols-2;	// bug fix 3/29/17, RD
                 				colcount++,magptr++,gxptr++,gyptr++,resultptr++){
                 				m00 = *magptr;
                 				if(m00 == 0){
                    					*resultptr = (unsigned char) NOEDGE;
                 				}
                 				else{
                    					xperp = -(gx = *gxptr)/((float)m00);
                    					yperp = (gy = *gyptr)/((float)m00);
                 				}

                 				if(gx >= 0){
                    					if(gy >= 0){
                            					if (gx >= gy)
                            					{
                             						/* 111 */
                                					/* Left point */
                                					z1 = *(magptr - 1);
                                					z2 = *(magptr - ncols - 1);

                                					mag1 = (m00 - z1)*xperp + (z2 - z1)*yperp;

                                					/* Right point */
                                					z1 = *(magptr + 1);
                                					z2 = *(magptr + ncols + 1);

                                					mag2 = (m00 - z1)*xperp + (z2 - z1)*yperp;
                            					}
                            					else
                            					{
                             						/* 110 */
	 				                                /* Left point */
                                					z1 = *(magptr - ncols);
                                					z2 = *(magptr - ncols - 1);

                                					mag1 = (z1 - z2)*xperp + (z1 - m00)*yperp;

                                					/* Right point */
                                					z1 = *(magptr + ncols);
                                					z2 = *(magptr + ncols + 1);

                                					mag2 = (z1 - z2)*xperp + (z1 - m00)*yperp;
                            					}
                        				}
                        				else
                        				{
                            					if (gx >= -gy)
                            					{
                             						/* 101 */
                                					/* Left point */
                                					z1 = *(magptr - 1);
                                					z2 = *(magptr + ncols - 1);

                                					mag1 = (m00 - z1)*xperp + (z1 - z2)*yperp;

                                					/* Right point */
                                					z1 = *(magptr + 1);
                                					z2 = *(magptr - ncols + 1);

                                					mag2 = (m00 - z1)*xperp + (z1 - z2)*yperp;
                            					}
                            					else
                            					{
                             						/* 100 */
                                					/* Left point */
                                					z1 = *(magptr + ncols);
                                					z2 = *(magptr + ncols - 1);

                                					mag1 = (z1 - z2)*xperp + (m00 - z1)*yperp;

                                					/* Right point */
                                					z1 = *(magptr - ncols);
                                					z2 = *(magptr - ncols + 1);

                                					mag2 = (z1 - z2)*xperp  + (m00 - z1)*yperp;
                            					}
                        				}
                    				}
                    				else
                    				{
                     					if ((gy = *gyptr) >= 0)
                        				{
                            					if (-gx >= gy)
                            					{
                             						/* 011 */
                                					/* Left point */
                                					z1 = *(magptr + 1);
                                					z2 = *(magptr - ncols + 1);

                                					mag1 = (z1 - m00)*xperp + (z2 - z1)*yperp;

                                					/* Right point */
                                					z1 = *(magptr - 1);
                                					z2 = *(magptr + ncols - 1);

                                					mag2 = (z1 - m00)*xperp + (z2 - z1)*yperp;
                            					}
                            					else
                            					{
                             						/* 010 */
                                					/* Left point */
                                					z1 = *(magptr - ncols);
                                					z2 = *(magptr - ncols + 1);

                                					mag1 = (z2 - z1)*xperp + (z1 - m00)*yperp;

                                					/* Right point */
                                					z1 = *(magptr + ncols);
                                					z2 = *(magptr + ncols - 1);

                                					mag2 = (z2 - z1)*xperp + (z1 - m00)*yperp;
                            					}
                        				}
                        				else
                        				{
                            					if (-gx > -gy)
                            					{
                             						/* 001 */
                                					/* Left point */
                                					z1 = *(magptr + 1);
                                					z2 = *(magptr + ncols + 1);

                                					mag1 = (z1 - m00)*xperp + (z1 - z2)*yperp;

                                					/* Right point */
                                					z1 = *(magptr - 1);
                                					z2 = *(magptr - ncols - 1);

                                					mag2 = (z1 - m00)*xperp + (z1 - z2)*yperp;
                            					}
                            					else
                            					{
                             						/* 000 */
                                					/* Left point */
                                					z1 = *(magptr + ncols);
                                					z2 = *(magptr + ncols + 1);

                                					mag1 = (z2 - z1)*xperp + (m00 - z1)*yperp;

                                					/* Right point */
                                					z1 = *(magptr - ncols);
                                					z2 = *(magptr - ncols - 1);

                                					mag2 = (z2 - z1)*xperp + (m00 - z1)*yperp;
                            					}
                        				}
                    				}

                    				/* Now determine if the current point is a maximum point */

                    				if ((mag1 > 0.0) || (mag2 > 0.0))
                    				{
                     					*resultptr = (unsigned char) NOEDGE;
                    				}
                    				else
                   	 			{
                     					if (mag2 == 0.0)
                            					*resultptr = (unsigned char) NOEDGE;
                        				else
                            					*resultptr = (unsigned char) POSSIBLE_EDGE;
                    				}
                			}
            			}
			}



			/*****************************************************************************************************/



			void action(){
                                while(1){

					pmag.read(mag);
					pgraddy.read(grady);
					pgraddx.read(gradx);
					/****************************************************************************
           				* Perform non-maximal suppression.
           				****************************************************************************/
           				if(VERBOSE) printf("Doing the non-maximal suppression.\n");
           				non_max_supp(ROWS, COLS);
					pnms.write(result);

                       		}
			}

                        SC_CTOR(non_max){
                                SC_THREAD(action);
                                SET_STACK_SIZE
                        }

		};

		SC_MODULE(apply_h){

			sc_fifo_in<IMAGE> pnms;
                        sc_fifo_in<SIMAGE> pmag;
			sc_fifo_out<IMAGE> pedge;
                        IMAGE nms;
                        SIMAGE mag;
			IMAGE edge;

			/*******************************************************************************************/
		        void follow_edges(unsigned char *edgemapptr, short *edgemagptr, short lowval,
           			int cols)
        		{
           			short *tempmagptr;
           			unsigned char *tempmapptr;
           			int i;
           			int x[8] = {1,1,0,-1,-1,-1,0,1},
               			y[8] = {0,1,1,1,0,-1,-1,-1};

           			for(i=0;i<8;i++){
              				tempmapptr = edgemapptr - y[i]*cols + x[i];
              				tempmagptr = edgemagptr - y[i]*cols + x[i];

              				if((*tempmapptr == POSSIBLE_EDGE) && (*tempmagptr > lowval)){
                 				*tempmapptr = (unsigned char) EDGE;
                 				follow_edges(tempmapptr,tempmagptr, lowval, cols);
              				}
           			}
        		}


			void apply_hysteresis( int rows, int cols)
        		{
           			int r, c, pos, numedges, highcount, lowthreshold, highthreshold, hist[32768];
           			short int maximum_mag=0;

           			/****************************************************************************
           			* Initialize the edge map to possible edges everywhere the non-maximal
           			* suppression suggested there could be an edge except for the border. At
           			* the border we say there can not be an edge because it makes the
           			* follow_edges algorithm more efficient to not worry about tracking an
           			* edge off the side of the image.
           			****************************************************************************/
           			for(r=0,pos=0;r<rows;r++){
              				for(c=0;c<cols;c++,pos++){
                 				if(nms[pos] == POSSIBLE_EDGE) edge[pos] = POSSIBLE_EDGE;
                 				else edge[pos] = NOEDGE;
              				}
           			}

           			for(r=0,pos=0;r<rows;r++,pos+=cols){
              				edge[pos] = NOEDGE;
              				edge[pos+cols-1] = NOEDGE;
           			}
           			pos = (rows-1) * cols;
           			for(c=0;c<cols;c++,pos++){
              				edge[c] = NOEDGE;
              				edge[pos] = NOEDGE;
           			}

           			/****************************************************************************
           			* Compute the histogram of the magnitude image. Then use the histogram to
           			* compute hysteresis thresholds.
           			****************************************************************************/
           			for(r=0;r<32768;r++) hist[r] = 0;
           			for(r=0,pos=0;r<rows;r++){
              				for(c=0;c<cols;c++,pos++){
                 				if(edge[pos] == POSSIBLE_EDGE) hist[mag[pos]]++;
              				}
           			}

           			/****************************************************************************
           			* Compute the number of pixels that passed the nonmaximal suppression.
           			****************************************************************************/
           			for(r=1,numedges=0;r<32768;r++){
              				if(hist[r] != 0) maximum_mag = r;
              				numedges += hist[r];
           			}

           			highcount = (int)(numedges * THIGH + 0.5);

           			/****************************************************************************
           			* Compute the high threshold value as the (100 * thigh) percentage point
           			* in the magnitude of the gradient histogram of all the pixels that passes
           			* non-maximal suppression. Then calculate the low threshold as a fraction
           			* of the computed high threshold value. John Canny said in his paper
           			* "A Computational Approach to Edge Detection" that "The ratio of the
           			* high to low threshold in the implementation is in the range two or three
           			* to one." That means that in terms of this implementation, we should
           			* choose tlow ~= 0.5 or 0.33333.
           			****************************************************************************/
           			r = 1;
           			numedges = hist[1];
           			while((r<(maximum_mag-1)) && (numedges < highcount)){
              				r++;
              				numedges += hist[r];
           			}
           			highthreshold = r;
           			lowthreshold = (int)(highthreshold * TLOW + 0.5);

           			if(VERBOSE){
              				printf("The input low and high fractions of %f and %f computed to\n",
                 			TLOW, THIGH);
              				printf("magnitude of the gradient threshold values of: %d %d\n",
                 			lowthreshold, highthreshold);
           			}

           			/****************************************************************************
           			* This loop looks for pixels above the highthreshold to locate edges and
           			* then calls follow_edges to continue the edge.
           			****************************************************************************/
           			for(r=0,pos=0;r<rows;r++){
              				for(c=0;c<cols;c++,pos++){
                 				if((edge[pos] == POSSIBLE_EDGE) && (mag[pos] >= highthreshold)){
                    					edge[pos] = EDGE;
                    					follow_edges((edge+pos), (mag+pos), lowthreshold, cols);
                 				}
              				}
           			}

           			/****************************************************************************
           			* Set all the remaining possible edges to non-edges.
           			****************************************************************************/
           			for(r=0,pos=0;r<rows;r++){
              				for(c=0;c<cols;c++,pos++) if(edge[pos] != EDGE) edge[pos] = NOEDGE;
           			}
			}




			/*******************************************************************************************/


			void action(){
				while(1){
					pnms.read(nms);
					pmag.read(mag);
					/****************************************************************************
           				* Use hysteresis to mark the edge pixels.
           				****************************************************************************/
           				if(VERBOSE) printf("Doing hysteresis thresholding.\n");
           				apply_hysteresis(ROWS,COLS);
					pedge.write(edge);
				}
			}

			SC_CTOR(apply_h){
				SC_THREAD(action);
				SET_STACK_SIZE
			}

		};





	mgk mgk1;
	derivx_y derivx_y1;
	magx_y magx_y1;
	non_max non_max1;
	apply_h apply_h1;

	void before_end_of_elaboration(){
           mgk1.ImgIn.bind(ImgIn);
           mgk1.ImgOut.bind(q1);
           derivx_y1.ImgIn.bind(q1);
           derivx_y1.deltax_to_mag.bind(q2);
	   derivx_y1.deltay_to_mag.bind(q3);
	   derivx_y1.deltax_to_non.bind(q4);
           derivx_y1.deltay_to_non.bind(q5);
	   magx_y1.Delta_x.bind(q2);
	   magx_y1.Delta_y.bind(q3);
	   magx_y1.mag1.bind(q6);
	   magx_y1.mag2.bind(q7);
	   non_max1.pmag.bind(q4);
	   non_max1.pgraddy.bind(q5);
	   non_max1.pgraddx.bind(q7);
	   non_max1.pnms.bind(q8);
	   apply_h1.pnms.bind(q8);
	   apply_h1.pmag.bind(q6);
	   apply_h1.pedge.bind(ImgOut);
        }

	SC_CTOR(DUT)
        :q1("q1",1)
        ,q2("q2",1)
	,q3("q3",1)
        ,q4("q4",1)
	,q5("q5",1)
        ,q6("q6",1)
	,q7("q7",1)
        ,q8("q8",1)
        ,mgk1("mgk1")
        ,derivx_y1("derivx_y1")
        ,magx_y1("magx_y1")
	,non_max1("non_max1")
	,apply_h1("apply_h1")
        {}


};

SC_MODULE(Platform)
{
	sc_fifo_in<IMAGE> ImgIn;
	sc_fifo_out<IMAGE> ImgOut;
	sc_fifo<IMAGE> q1;
	sc_fifo<IMAGE> q2;

	DataIn din;
	DUT canny;
	DataOut dout;

	void before_end_of_elaboration(){
	   din.ImgIn.bind(ImgIn);
	   din.ImgOut.bind(q1);
	   canny.ImgIn.bind(q1);
	   canny.ImgOut.bind(q2);
	   dout.ImgIn.bind(q2);
	   dout.ImgOut.bind(ImgOut);
	}

	SC_CTOR(Platform)
	:q1("q1",1)
	,q2("q2",1)
	,din("din")
	,canny("canny")
	,dout("dout")
	{}
};

SC_MODULE(Top)
{
	sc_fifo<IMAGE> q1;
	sc_fifo<IMAGE> q2;
	Stimulus stimulus;
	Platform platform;
	Monitor monitor;

	void before_end_of_elaboration(){
	   stimulus.ImgOut.bind(q1);
	   platform.ImgIn.bind(q1);
	   platform.ImgOut.bind(q2);
	   monitor.ImgIn.bind(q2);
	}

	SC_CTOR(Top)
	:q1("q1",1)
	,q2("q2",1)
	,stimulus("stimulus")
	,platform("platform")
	,monitor("monitor")
	{}
};

Top top("top");

int sc_main(int argc, char* argv[])
{
	sc_start();
	return 0;
}
