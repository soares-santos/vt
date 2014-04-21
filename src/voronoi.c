#include "voronoi.h"

void set_defaults(void){
  int i;
  Verbose=0;
  Holes=0;
  XYcoord=0;
  Ebeling=1;
  Poisson=0;
  Fits=1;
  NotFind=0;
  RawFind=0;
  area_mean=0.;
  tot_area_covered=0.;
  background_density=0.;
  conf_lev=0.8;
  reject_lev=0.1;
  w_amp=0.005;
  w_pow=1.7;
  ndata_ok=0;
  FieldID=0;
  for(i=0;i<4;i++){
    ttype[i]  = (char *) malloc(FLEN_VALUE);
    ttypeH[i] = (char *) malloc(FLEN_VALUE);
  }
  ttype[0]="id";    ttypeH[0]="ramin";    frame[0]=0; 
  ttype[1]="ra";    ttypeH[1]="ramax";    frame[1]=360;
  ttype[2]="dec";   ttypeH[2]="decmin";   frame[2]=-90;
  ttype[3]="z";     ttypeH[3]="decmax";   frame[3]=90;
}

void print_init(void){
  if(Verbose){
    printf("Initial parameters and flags:\n");
    printf("name       value\n");
    printf("A          %f\n",w_amp);
    printf("g          %f\n",w_pow);
    printf("r          %f\n",reject_lev);
    if(!Ebeling) printf("s          %f\n",conf_lev);
    printf("C          %s %s %s %s\n",ttype[0],ttype[1],ttype[2],ttype[3]);
    printf("B          %s %s %s %s\n",ttypeH[0],ttypeH[1],ttypeH[2],ttypeH[3]);
    printf("F          %f %f %f %f\n",frame[0],frame[1],frame[2],frame[3]);
    if(Holes) printf("holesfile  %s\n",file_holes);
    printf("datafile   %s\n",file_name);
    printf("Verbose    %d\n",Verbose);
    printf("Holes      %d\n",Holes);
    printf("Ebeling    %d\n",Ebeling);
    printf("Poisson    %d\n",Poisson);
    printf("Fits       %d\n",Fits);
    printf("XYcoord    %d\n",XYcoord);
    printf("NotFind    %d\n",NotFind);
    printf("RawFind    %d\n",RawFind);
    printf("FieldID    %lld\n",FieldID);
  }
}

void set_TF_params(){
  static REAL TFpara[20][10]= //alpha 
    { {4.     ,4.     ,4.     ,4.     ,4.     ,4.      ,4.     ,4.     ,4.     ,4.    },
      {3.90733,3.86026,3.93163,3.91549,3.8146 ,3.96033,3.91312,3.86322,3.81397,3.93777},
      {3.94332,3.87495,3.93277,3.874  ,3.82526,3.89633,3.95259,3.78706,3.83847,3.89183},
      {3.81324,3.89603,3.93585,3.8392 ,3.85502,3.96963,3.80767,3.80239,3.81277,3.86103},
      {3.86312,3.80879,3.78786,3.86579,3.84952,3.969  ,3.86682,3.81766,3.90912,3.89701},
      {3.80722,3.85894,3.84802,3.79815,3.8298 ,3.87444,3.80618,3.89002,3.95618,3.88399},
      {3.90213,3.77684,3.84351,3.8769 ,3.82756,3.85622,3.71072,3.85908,3.92118,3.91385},
      {3.84708,3.89556,3.84301,3.81757,3.88597,3.808  ,3.83862,3.77147,3.75256,3.86434},
      {3.85582,3.86727,3.9242 ,3.81602,3.83066,3.98575,3.86409,3.87609,3.8761 ,3.86004},
      {3.90152,3.95895,3.96218,3.78374,3.90743,3.86134,3.92896,3.84033,3.84875,3.94683},
      {3.90049,3.93273,3.92551,3.82574,3.93574,3.7996 ,3.87878,3.81869,3.87702,3.81574},
      {3.9328 ,3.83791,3.8751 ,3.92334,3.90891,3.75224,3.78868,3.93861,3.93153,3.88375},
      {3.73865,3.87765,3.90657,3.89199,3.89202,3.77922,3.78771,3.82214,3.84844,3.85797},
      {3.91013,3.89401,3.91808,3.7988 ,3.9347 ,3.96778,3.85649,3.82778,3.77219,3.81238},
      {3.87384,3.85139,3.79959,3.88255,3.89718,3.95738,3.8504 ,3.89141,3.76614,3.71983},
      {3.88291,3.85755,3.80339,3.90986,3.85826,3.76541,3.88196,3.72906,3.75739,3.81976},
      {3.85453,3.74402,3.87299,3.86607,3.85779,3.90817,3.78404,3.80403,3.78011,3.871  },
      {3.86361,3.90158,3.87814,3.8911 ,3.79042,3.78792,3.90431,3.76375,3.69806,3.92489},
      {3.87031,3.92482,3.97383,3.83817,3.94345,3.74538,3.94648,3.78796,3.91634,3.86448},
      {3.86695,3.91418,3.93322,3.85489,3.858  ,3.88375,3.85168,3.79989,3.73856,3.75517}
    };
  static REAL TFparb[20][10]= //beta
    { {4.     ,4.     ,4.     ,4.     ,4.     ,4.     ,4.     ,4.     ,4.     ,4.     },
      {3.65798,3.61274,3.6826 ,3.66807,3.56888,3.71208,3.65469,3.61503,3.56672,3.6946 },
      {3.70586,3.62721,3.68964,3.62823,3.57834,3.65822,3.70262,3.54818,3.59306,3.64921},
      {3.56253,3.65117,3.69038,3.59607,3.60731,3.72108,3.56535,3.55447,3.56547,3.61036},
      {3.62441,3.56328,3.54607,3.62183,3.60498,3.73081,3.62803,3.56904,3.66496,3.64874},
      {3.55676,3.60657,3.59899,3.54745,3.58359,3.63186,3.56736,3.64525,3.69275,3.63923},
      {3.65864,3.53545,3.60548,3.62894,3.58854,3.60811,3.47476,3.60775,3.67349,3.65714},
      {3.59523,3.64289,3.604  ,3.57243,3.64275,3.55649,3.5912 ,3.52175,3.50378,3.61306},
      {3.60738,3.6248 ,3.68473,3.57191,3.58975,3.72901,3.61207,3.63019,3.62307,3.59754},
      {3.65609,3.70244,3.70517,3.53193,3.65121,3.62769,3.66595,3.59069,3.60356,3.69434},
      {3.65013,3.68771,3.67827,3.57753,3.69118,3.5552 ,3.62682,3.56678,3.62355,3.56251},
      {3.68463,3.60022,3.62635,3.6698 ,3.6591 ,3.50324,3.54812,3.68265,3.67449,3.61736},
      {3.48998,3.62914,3.66446,3.64605,3.63116,3.53819,3.53823,3.57267,3.58521,3.59715},
      {3.65556,3.65211,3.66215,3.55802,3.67779,3.7073 ,3.61068,3.57421,3.50776,3.54477},
      {3.63164,3.59817,3.54631,3.62592,3.63562,3.69417,3.58501,3.6222 ,3.49923,3.44695},
      {3.62394,3.60621,3.5469 ,3.65687,3.59589,3.51836,3.61431,3.46163,3.48937,3.55022},
      {3.6056 ,3.4964 ,3.61948,3.61035,3.60182,3.65008,3.5315 ,3.54242,3.50412,3.58632},
      {3.60428,3.63853,3.62489,3.63258,3.5377 ,3.53322,3.63575,3.50439,3.42142,3.62975},
      {3.609  ,3.67394,3.7076 ,3.58518,3.66587,3.49304,3.67251,3.52723,3.63369,3.58337},
      {3.61502,3.6589 ,3.68932,3.5966 ,3.59605,3.61797,3.5726 ,3.5304 ,3.46445,3.49211}
    };
  static REAL TFpar0[20][10]= //c0
    { {-0.047 ,-0.047 ,-0.047 ,-0.047 ,-0.047 ,-0.047 ,-0.047 ,-0.047 ,-0.047 ,-0.047  },
      {-0.0361,-0.0177,-0.0365,-0.0306,-0.0492,-0.0173,-0.0315,-0.0269,0.00136,-0.0438},
      {-0.0388,-0.0527,-0.0282,-0.0699,-0.0207,-0.042 ,-0.0167,-0.0476,-0.0059,-0.0481},
      {-0.049 ,-0.0297,-0.0193,-0.0609,-0.0472,-0.0248,-0.0419,-0.0355,-0.0397,-0.0221},
      {-0.0091,-0.0044,-0.034 ,-0.0258,-0.0321,-0.0161,0.00438,-0.0562,-0.0427,-0.0432},
      {-0.0462,-0.0138,-0.0317,-0.0384,-0.0141,-0.0432,-0.0197,-0.0414,-0.0129,-0.012 },
      {-0.0156,-0.0094,-0.0353,-0.0469,-0.0269,-0.0534,-0.0332,-0.0242,-0.0322,-0.0261},
      {-0.0358,-0.039 ,-0.037 ,-0.0089,-0.0289,-0.0321,-0.0067,-0.0209,0.0371 ,-0.0328},
      {-0.0065,-0.0313,-0.0152,-0.0532,-0.0333,-0.046 ,-0.0293,-0.0379,-0.0204,-0.0206},
      {-0.0393,-0.024 ,-0.0228,-0.0366,-0.0199,-0.0277,-0.0414,-0.0233,-0.028 ,-0.0425},
      {-0.0355,-0.0402,-0.0389,-0.0201,-0.0212,-0.0141,-0.0547,-0.0561,-0.0592,-0.0274},
      {-0.0499,-0.0364,-0.0316,-0.0114,-0.0331,0.0209 ,-0.0184,-0.0127,-0.0050,-0.0481},
      {0.0119 ,-0.0206,-0.026 ,-0.0369,-0.0367,-0.0433,-0.0171,-0.0185,-0.0426,-0.038 },
      {-0.03  ,-0.0402,-0.0273,-0.0113,-0.0253,-0.0477,-0.0347,-0.0482,0.00285,-0.0209},
      {-0.0188,-0.0107,-0.0302,-0.012 ,-0.0225,-0.0432,-0.0345,-0.0604,-0.0084,-0.0474},
      {-0.0194,-0.0165,-0.0256,-0.0329,-0.0219,-0.031 ,-0.027 ,0.027  ,0      ,-0.0475},
      {-0.0149,-0.0184,-0.0521,-0.0307,-0.0199,-0.0219,0.00491,-0.0225,-0.0546,-0.0098},
      {-0.0041,-0.0057,-0.0126,-0.0272,-0.021 ,0.011  ,-0.0196,-0.0049,0.0183 ,-0.0328},
      {-0.0202,-0.0382,-0.0246,0.0023 ,-0.0258,-0.0441,-0.0079,0.0218 ,-0.0176,-0.0586},
      {-0.0205,-0.0294,-0.0403,-0.033 ,-0.0345,-0.0323,-0.0245,-0.0098,0.0761 ,-0.0153}
    };
  static REAL TFpar1[20][10]= //c1
    { {0.04   ,0.04   ,0.04   ,0.04   ,0.04   ,0.04   ,0.04   ,0.04   ,0.04   ,0.04  },
      {0.0863 ,0.0518 ,0.0828 ,0.0757 ,0.117  ,0.0436 ,0.0762 ,0.0689 ,0.0195 ,0.0987 },
      {0.0884 ,0.124  ,0.0713 ,0.157  ,0.0634 ,0.0952 ,0.0441 ,0.116  ,0.0325 ,0.116  },
      {0.118  ,0.0746 ,0.0484 ,0.142  ,0.113  ,0.0573 ,0.104  ,0.0921 ,0.0988 ,0.0595 },
      {0.0353 ,0.0308 ,0.0873 ,0.065  ,0.0839 ,0.0429 ,0.0123 ,0.133  ,0.0969 ,0.0969 },
      {0.113  ,0.0455 ,0.0974 ,0.0976 ,0.0508 ,0.0974 ,0.0551 ,0.094  ,0.037  ,0.0433 },
      {0.048  ,0.04   ,0.0836 ,0.113  ,0.0749 ,0.126  ,0.0875 ,0.0625 ,0.0789 ,0.0683 },
      {0.085  ,0.0895 ,0.0865 ,0.035  ,0.0713 ,0.0908 ,0.0352 ,0.0641 ,-0.0496,0.0802 },
      {0.0316 ,0.0751 ,0.047  ,0.126  ,0.0804 ,0.103  ,0.0727 ,0.0955 ,0.0551 ,0.0562 },
      {0.0899 ,0.0561 ,0.0536 ,0.102  ,0.0553 ,0.0702 ,0.0934 ,0.0614 ,0.0713 ,0.0942 },
      {0.0859 ,0.092  ,0.0904 ,0.0618 ,0.0583 ,0.045  ,0.13   ,0.132  ,0.137  ,0.0684 },
      {0.122  ,0.0927 ,0.0768 ,0.0413 ,0.0801 ,-0.0177,0.0599 ,0.0476 ,0.0321 ,0.116  },
      {0.00058,0.0568 ,0.0672 ,0.0872 ,0.0859 ,0.109  ,0.0561 ,0.0527 ,0.106  ,0.0961 },
      {0.0753 ,0.0921 ,0.0698 ,0.0389 ,0.0652 ,0.105  ,0.0811 ,0.115  ,0.0151 ,0.061  },
      {0.0521 ,0.0388 ,0.0826 ,0.0421 ,0.0611 ,0.0961 ,0.0917 ,0.14   ,0.0392 ,0.122  },
      {0.0555 ,0.0497 ,0.0725 ,0.0784 ,0.0603 ,0.0855 ,0.0708 ,-0.0311,0      ,0.128  },
      {0.0472 ,0.0642 ,0.124  ,0.0761 ,0.0595 ,0.0584 ,0.0134 ,0.0601 ,0.138  ,0.0413 },
      {0.0294 ,0.0296 ,0.0419 ,0.0694 ,0.0623 ,0.00554,0.0555 ,0.0313 ,-0.0149,0.079  },
      {0.056  ,0.0884 ,0.0575 ,0.0174 ,0.0662 ,0.113  ,0.0352 ,-0.0182,0.0507 ,0.136  },
      {0.0568 ,0.0731 ,0.0919 ,0.0856 ,0.0809 ,0.0786 ,0.0679 ,0.0404 ,-0.136 ,0.0533 }
    };
  static REAL TFpar2[20][10]= //c2
    { {0.62   ,0.62   ,0.62   ,0.62   ,0.62   ,0.62   ,0.62   ,0.62   ,0.62   ,0.62  },
      {0.759  ,0.685  ,0.851  ,0.867  ,0.64   ,0.624  ,0.779  ,1.08   ,0.585  ,0.65  },
      {0.732  ,0.702  ,0.566  ,0.646  ,0.671  ,0.791  ,0.95   ,0.714  ,0.636  ,0.724 },
      {0.824  ,0.662  ,0.807  ,0.891  ,0.769  ,0.793  ,0.799  ,0.69   ,0.787  ,0.683 },
      {0.728  ,0.64   ,0.636  ,0.638  ,0.716  ,0.863  ,0.609  ,0.795  ,0.874  ,0.788 },
      {0.799  ,0.798  ,0.521  ,0.81   ,0.649  ,0.783  ,0.579  ,0.787  ,0.869  ,0.861 },
      {0.913  ,0.618  ,0.873  ,0.745  ,0.718  ,0.724  ,0.651  ,0.698  ,0.955  ,0.808 },
      {0.93   ,0.949  ,0.787  ,0.745  ,0.813  ,0.643  ,0.759  ,0.815  ,0.608  ,1.13  },
      {0.58   ,0.998  ,0.717  ,0.798  ,1.02   ,0.919  ,0.771  ,0.69   ,0.848  ,0.81  },
      {1.07   ,0.809  ,0.837  ,0.698  ,0.857  ,0.702  ,0.788  ,0.953  ,0.747  ,0.806 },
      {1.05   ,1.02   ,0.733  ,0.714  ,0.783  ,0.721  ,0.927  ,0.894  ,0.811  ,0.656 },
      {0.855  ,0.696  ,1.01   ,0.826  ,0.923  ,0.669  ,0.832  ,0.771  ,0.635  ,0.711 },
      {0.572  ,0.766  ,0.941  ,1.05   ,0.914  ,0.986  ,0.693  ,0.702  ,0.796  ,0.66  },
      {0.84   ,0.725  ,0.923  ,0.735  ,0.784  ,0.976  ,0.794  ,0.599  ,0.515  ,0.528 },
      {0.667  ,0.915  ,0.758  ,0.766  ,0.853  ,0.956  ,0.664  ,0.753  ,0.681  ,0.598 },
      {0.869  ,0.654  ,0.741  ,0.78   ,0.737  ,0.96   ,0.896  ,0.583  ,0      ,0.492 },
      {0.678  ,0.618  ,0.937  ,0.911  ,0.634  ,0.673  ,0.621  ,0.686  ,0.621  ,0.656 },
      {0.73   ,0.812  ,0.66   ,0.846  ,0.647  ,0.523  ,0.693  ,0.656  ,0.564  ,0.86  },
      {0.666  ,0.79   ,1.     ,0.661  ,0.771  ,0.581  ,0.735  ,0.794  ,0.591  ,0.713 },
      {0.886  ,0.778  ,0.798  ,0.727  ,0.925  ,0.884  ,0.645  ,0.688  ,0.488  ,0.725 }
    };
  static REAL TFpar3[20][10]= //c3
    { {-0.45 ,-0.45 ,-0.45 ,-0.45 ,-0.45 ,-0.45 ,-0.45 ,-0.45 ,-0.45 ,-0.45 },
      {-0.892,-0.816,-1.04 ,-1.07 ,-0.807,-0.648,-0.928,-1.37 ,-0.739,-0.727},
      {-0.834,-0.893,-0.623,-0.803,-0.853,-0.946,-1.1  ,-0.91 ,-0.818,-0.919},
      {-1.08 ,-0.75 ,-0.898,-1.16 ,-0.999,-0.86 ,-1.04 ,-0.882,-1.03 ,-0.799},
      {-0.906,-0.82 ,-0.808,-0.745,-0.929,-0.982,-0.789,-1.03 ,-1.06 ,-0.944},
      {-1.04 ,-0.981,-0.74 ,-1.06 ,-0.82 ,-0.935,-0.658,-0.947,-1.01 ,-1.07 },
      {-1.16 ,-0.79 ,-1.08 ,-0.968,-0.929,-0.929,-0.829,-0.828,-1.18 ,-0.966},
      {-1.12 ,-1.16 ,-0.954,-0.916,-1.   ,-0.877,-1.   ,-1.07 ,-0.831,-1.42 },
      {-0.664,-1.27 ,-0.859,-1.03 ,-1.29 ,-1.1  ,-0.928,-0.888,-1.05 ,-1.   },
      {-1.35 ,-0.887,-0.932,-0.956,-1.07 ,-0.811,-0.915,-1.2  ,-0.895,-0.946},
      {-1.32 ,-1.25 ,-0.856,-0.924,-0.954,-0.867,-1.21 ,-1.16 ,-1.05 ,-0.777},
      {-1.11 ,-0.888,-1.27 ,-1.01 ,-1.13 ,-0.868,-1.1  ,-1.01 ,-0.808,-0.902},
      {-0.773,-0.925,-1.16 ,-1.31 ,-1.11 ,-1.31 ,-0.894,-0.846,-1.03 ,-0.841},
      {-1.   ,-0.849,-1.13 ,-0.913,-0.955,-1.18 ,-0.961,-0.746,-0.658,-0.658},
      {-0.805,-1.16 ,-0.98 ,-0.934,-1.05 ,-1.19 ,-0.833,-0.962,-0.878,-0.806},
      {-1.08 ,-0.775,-0.962,-0.928,-0.875,-1.28 ,-1.09 ,-0.794,0     ,-0.699},
      {-0.8  ,-0.842,-1.24 ,-1.12 ,-0.811,-0.812,-0.769,-0.811,-0.842,-0.837},
      {-0.957,-1.02 ,-0.779,-1.02 ,-0.828,-0.742,-0.818,-0.847,-0.768,-1.05 },
      {-0.778,-0.934,-1.17 ,-0.864,-0.925,-0.79 ,-0.888,-1.08 ,-0.683,-0.902},
      {-1.1  ,-0.926,-0.96 ,-0.942,-1.16 ,-1.08 ,-0.832,-0.902,-0.691,-0.939}
    };

  REAL dA=0.001,dg=0.1,A_0=0.,g_0=1.0;
  if(w_amp > 0.01) dA=0.01;
  int i,j;

  if(Poisson) i=0,j=0;
  else{
    i=(int) ((w_amp-A_0)/dA);
    if(i<0) i=0; 
    if(w_amp > 0.01) i+=10;
    if(i>19) i=19;
    j=(int) ((w_pow-g_0)/dg)+1;
    if(j<0) j=0;
    if(j>9) j=9;
  }
  alpha=TFpara[i][j];
  beta =TFparb[i][j];
  c0=TFpar0[i][j];
  c1=TFpar1[i][j];
  c2=TFpar2[i][j];
  c3=TFpar3[i][j];
  if(Verbose){
    printf("Background distribution parameters:\n");
    printf(" (alpha,beta): ( %f , %f )\n",alpha,beta);
    printf(" (c0,c1,c2,c3): ( %f , %f , %f , %f )\n",c0,c1,c2,c3);
  }
}

void read_data(void) {
  FILE *fp;
  register int i;
  REAL a, b, c;
  fitsfile *fptr;
  int status=0,hdunum=2,hdutype,frow=1,felem=1,colnum,anynull;
  REAL dnull=0., *tmp1, *tmp2;
  long nrows;
  long long longnull=0;
  if(Fits){ /* reading from fits file */
    fits_open_file(&fptr,file_name,READONLY,&status), cfitsio_error(status);
    fits_movabs_hdu(fptr,hdunum,&hdutype,&status);
    if(hdutype != BINARY_TBL)
      fprintf(stderr,"Error: Binary table expected in HDU #2.\n"),exit(1);
    fits_get_num_rows(fptr, &nrows, &status);
    ndata = (int) nrows;
    if(ndata == 0) fprintf(stderr,"Error: No data points in file.\n"),exit(1);
    in.numberofpoints = ndata;
    in.pointlist = (REAL *) malloc(ndata * 2 * sizeof(REAL));
    tmp1 = (REAL *) malloc(ndata*sizeof(REAL));
    tmp2 = (REAL *) malloc(ndata*sizeof(REAL));
    fits_get_colnum(fptr,CASEINSEN,ttype[1],&colnum,&status); 
    fits_read_col(fptr,TDOUBLE,colnum,frow,felem,ndata,&dnull,tmp1,
		  &anynull,&status), cfitsio_error(status);
    fits_get_colnum(fptr,CASEINSEN,ttype[2],&colnum,&status);
    fits_read_col(fptr,TDOUBLE,colnum,frow,felem,ndata,&dnull,tmp2,
		  &anynull,&status), cfitsio_error(status);
    for(i=0;i<ndata;i++) in.pointlist[2*i]=tmp1[i],in.pointlist[2*i+1]=tmp2[i];
    free(tmp1); free(tmp2);
    number = (long long *) malloc(ndata*sizeof(long long));
    redsh = (REAL *) malloc(ndata*sizeof(REAL));
    redsh1 = (REAL *) malloc(ndata*sizeof(REAL));
    fits_get_colnum(fptr,CASEINSEN,ttype[0],&colnum,&status);
    fits_read_col(fptr,TLONGLONG,colnum,frow,felem,ndata,&longnull,number,
		  &anynull,&status), cfitsio_error(status);
    fits_get_colnum(fptr,CASEINSEN,ttype[3],&colnum,&status);
    fits_read_col(fptr,TDOUBLE,colnum,frow,felem,ndata,&dnull,redsh,
		  &anynull,&status), cfitsio_error(status);
    fits_get_colnum(fptr,CASEINSEN,ttype[3],&colnum,&status);
    fits_read_col(fptr,TDOUBLE,colnum,frow,felem,ndata,&dnull,redsh1,
		  &anynull,&status), cfitsio_error(status);
    fits_close_file(fptr,&status);
  } else { /* reading from ascii file */
    if((fp=fopen(file_name, "r")) == NULL){
      fprintf(stderr, "Error: Failed to open file %s\n", file_name);
      exit(1);
    } else {
      ndata=0;
      skipcomments(fp);
      while((i=getc(fp)) != EOF) if(i == '\n') ndata++;
      rewind(fp);
      in.numberofpoints = ndata;
      in.pointlist = (REAL *) malloc(ndata * 2 * sizeof(REAL));
      redsh = (REAL *) malloc(ndata*sizeof(REAL));
      redsh1 = (REAL *) malloc(ndata*sizeof(REAL));
      number = (long long *) malloc(ndata*sizeof(long long));
      for(i=0;i<ndata;i++) {
	skipcomments(fp);
	fscanf(fp,"%lld %lf %lf %lf", &number[i], &a, &b, &c);
	in.pointlist[2*i]=a, in.pointlist[2*i+1]=b,redsh[i]=c, redsh1[i]=c;
      }
    }
    fclose(fp);
  } 
  /* find the center of the box */
  redshmin=xmin=ymin=1.e10; 
  redshmax=xmax=ymax=-1.e10;
  for(i=0;i<ndata;i++) {
    a=in.pointlist[2*i];
    b=in.pointlist[2*i+1];
    c=redsh[i];
    MAXMIN(a,xmax,xmin); 
    MAXMIN(b,ymax,ymin); 
    MAXMIN(c,redshmax,redshmin);
  }
  xmean=(xmax+xmin)*0.5;
  ymean=(ymax+ymin)*0.5;
  zmean=(redshmax+redshmin)*0.5;
  medredsh=median(redsh1,ndata);
  /* bring the box center to (0,0,zmean) */
  redshmin=xmin=ymin=1.e10; 
  redshmax=xmax=ymax=-1.e10;
  for(i=0;i<ndata;i++) {
    a=in.pointlist[2*i] = (in.pointlist[2*i] -xmean);
    if (! XYcoord) 
      a=in.pointlist[2*i]=in.pointlist[2*i]*cos(in.pointlist[2*i+1]*M_PI/180.);
    b=in.pointlist[2*i+1] -= ymean;
    c=redsh[i];
    MAXMIN(a,xmax,xmin); 
    MAXMIN(b,ymax,ymin); 
    MAXMIN(c,redshmax,redshmin);
  }
  in.numberofpointattributes = in.numberofsegments = 0;
  in.numberofholes = in.numberofregions = 0;
  in.pointmarkerlist = (int *) NULL;
  /* output summary */
  if(Verbose) {
    printf("Input galaxy catalog: %s\n",file_name);
    printf("Number of data points: %d \n",ndata);
    printf("Median of redshifts: %f \n",medredsh);
    printf("Box center: ( %14.10f , %14.10f , %14.10f )\n",xmean,ymean,zmean);
    printf("Box boundaries after shifting the center to (0,0,%f):\n");
    printf("%f < x < %f \n", xmin, xmax);
    printf("%f < y < %f \n", ymin, ymax);
    printf(" %f < z < %f \n",redshmin,redshmax);
  }
} /* read_data() */

REAL median(REAL *mag1, int ndata){
  REAL magmed;
  int n2;
  qsort((REAL *) mag1, ndata, sizeof(REAL), real_compare);
  magmed= mag1[0];
  n2=ndata/2;
  if (2*n2 == ndata) magmed=0.5*(mag1[n2]+mag1[n2+1]);
  else magmed=mag1[n2+1];
  return magmed;
} 

int real_compare(const void *i, const void *j) {
  if ( *(REAL *)i > *(REAL *)j) return 1;
  if ( *(REAL *)i < *(REAL *)j) return -1;
  return 0;
} 

void cfitsio_error(int status){
  if(status) fits_report_error(stderr,status), exit(status);    
  return;
} 
  
void read_holes(){  
  FILE *fp;  
  register int i,j;  
  int i1,i2,i3,i4;
  REAL hfraction=0,harea=0,c;
  fitsfile *fptr;
  int status=0,hdunum=2,hdutype,frow=1,felem=1,colnum,anynull;
  REAL dnull=0., *tmp1;
  long nrows;
  if(Fits){ /* reading from fits file */
    fits_open_file(&fptr,file_holes,READONLY,&status), cfitsio_error(status);
    fits_movabs_hdu(fptr,hdunum,&hdutype,&status);
    if(hdutype != BINARY_TBL)
      fprintf(stderr,"Error: Binary table expected in HDU #2.\n"),exit(1);
    fits_get_num_rows(fptr, &nrows, &status);
    nholes = (int) nrows;
    HolesCoord = (REAL *) malloc(nholes*4*sizeof(REAL)); 
    tmp1 = (REAL *) malloc(nholes*sizeof(REAL));
    for(i=0;i<4;i++){
      fits_get_colnum(fptr,CASEINSEN,ttypeH[i],&colnum,&status); 
      fits_read_col(fptr,TDOUBLE,colnum,frow,felem,nholes,&dnull,tmp1,
		    &anynull,&status), cfitsio_error(status);
      for(j=0;j<nholes;j++)
	if(i<2) HolesCoord[4*j+i] = tmp1[j] - xmean;       
	else    HolesCoord[4*j+i] = tmp1[j] - ymean;       
    }
    free(tmp1);
    fits_close_file(fptr,&status);
  } else { /* reading from ascii file */
    if((fp=fopen(file_holes, "r")) == NULL){  
      fprintf(stderr, "Error: Failed to open file %s\n", file_holes);  
      exit(1);  
    } else {  
      nholes=0; 
      skipcomments(fp);
      while((i=getc(fp)) != EOF) if(i == '\n') nholes++; 
      rewind(fp); 
      HolesCoord = (REAL *) malloc(nholes*4*sizeof(REAL)); 
      for(i=0;i<nholes;i++) { 
	i1=4*i; i2=i1+1; i3=i2+1; i4=i3+1; 
	fscanf(fp,"%lf %lf %lf %lf", &HolesCoord[i1], &HolesCoord[i2],  
	       &HolesCoord[i3], &HolesCoord[i4]); 
	HolesCoord[i1] -= xmean; HolesCoord[i2] -= xmean; 
	HolesCoord[i3] -= ymean; HolesCoord[i4] -= ymean; 
      } 
    } 
    fclose(fp); 
  }
  /* adjust if(!XYcoord) and compute fraction of area in holes */
  for(i=0;i<nholes;i++) {
    if (!XYcoord) {
      HolesCoord[4*i] = HolesCoord[4*i]*
	cos(0.5*(HolesCoord[4*i+2]+HolesCoord[4*i+3])*M_PI/180.); 
      HolesCoord[4*i+1] = HolesCoord[4*i+1]*
	cos(0.5*(HolesCoord[4*i+2]+HolesCoord[4*i+3])*M_PI/180.);
    } 
    harea+=
      (HolesCoord[4*i+1]-HolesCoord[4*i])*
      (HolesCoord[4*i+3]-HolesCoord[4*i+2]);
  }
  hfraction=harea/((xmax-xmin)*(ymax-ymin));
  /* summary */
  if(Verbose) {
    printf("Mask file: %s\n",file_holes);
    printf("Number of holes: %d \n",nholes);
    printf("Masked/Total area fraction: %f \n",hfraction);
  }
} /* read_holes() */   

void adjust_frame(void){
  REAL tmp=1.;
  frame[0]-=xmean;
  frame[1]-=xmean;
  if(! XYcoord){
    tmp=cos(0.5*(frame[2]+frame[3])*M_PI/180);
    frame[0]*=tmp;
    frame[1]*=tmp;
  }
  frame[2]-=ymean;
  frame[3]-=ymean;
  //  printf("%f, %f, %f ,%f\n",frame[0],frame[1],frame[2],frame[3]);
}


REAL poly_area(int n) {
  int ia=edge[n], ib=edge[n+1], /* first and last index of verts */
    nv,                         /* number of verts */
    *idx,                       /* array of ordered verts */
    i_edge, i_edge2,            /* temporary variables */
    i1, i2, i3, i4, i_next;     /* temporary indexes */
  register int i,j,k;
  REAL x1, y1, x2, y2,          /* verts coordinates */
    area=0.;                    /* polygon area */

  if(ia == ib) return 0.;   /* Bug correction ? */
  /* Check if the verts are connected */
  for(i=ia; i<ib; i++) {
    i_edge=2*edge[i];
    if(vor.edgelist[i_edge] == -1 || vor.edgelist[i_edge+1] == -1)
      return 0.;
  }
  /* Write index of verts in (counter)clockwise order */
  nv=ib-ia;
  idx=(int *) malloc(sizeof(int)*nv);
  i_edge=2*edge[ia];
  idx[0]=vor.edgelist[i_edge], i_next=vor.edgelist[i_edge+1];
  for(j=1;j<nv;j++) for(i=ia;i<ib;i++) 
    if((i_edge2=2*edge[i]) != i_edge) {
      i1=vor.edgelist[i_edge2], i2=vor.edgelist[i_edge2+1];
      if (i1 == i_next) {idx[j]=i1, i_edge=i_edge2, i_next=i2; break;}
      else if (i2 == i_next) {idx[j]=i2, i_edge=i_edge2, i_next=i1; break;}
    }
  /* Compute the area of the polygon if it is inside the limits */
  j=2*idx[nv-1];
  x1=vor.pointlist[j]; y1=vor.pointlist[j+1];
  if(x1 > xmax || x1 < xmin || y1 > ymax || y1 < ymin) {
    free(idx); return 0.;
  }
  for(i=0; i<nv; i++) {
    j=2*idx[i];
    x2=vor.pointlist[j]; y2=vor.pointlist[j+1];
    if(x2 > xmax || x2 < xmin || y2 > ymax || y2 < ymin) {
      free(idx); return 0.;
    }
    /* check if (x1,y1)(x2,y2) segment intersects a rectangle */
    if(Holes==1)
      for(k=0;k<nholes;k++) {
        i1=4*k; i2=i1+1; i3=i2+1; i4=i3+1;
        if(intersect(x1,y1,x2,y2,HolesCoord[i1],HolesCoord[i2],
                     HolesCoord[i3],HolesCoord[i4])) {
          free(idx); return 0.;}
      }
    area += x2*y1-y2*x1;
    x1=x2, y1=y2;
  }
  if(area<0.) area *=-1.;
  for(i=ia; i<ib; i++) voronoi_inside[edge[i]]=1;
  free(idx);
  return area*0.5;
} /* poly_area() */

int intersect (REAL x1, REAL y1, REAL x2, REAL y2,
               REAL xrmin, REAL xrmax, REAL yrmin, REAL yrmax) {
  REAL xsmin=x1, xsmax=x2,
    ysmin=y1, ysmax=y2,               /* Limits of the segments */
    a, m;                             /* Intersection, ang. coeff. */

  if (xsmax < xsmin) {xsmin = x2; xsmax = x1;}
  if (ysmax < ysmin) {ysmin = y2; ysmax = y1;}
  if (y1 != y2) {
    m=(x2-x1)/(y2-y1);
    if (ysmin <= yrmin && yrmin <= ysmax) {
      a = x1+(yrmin-y1)*m;
      if (xrmin <= a && a <= xrmax) return 1;
    }
    if (ysmin <= yrmax && yrmax <= ysmax) {
      a = x1+(yrmax-y1)*m;
      if (xrmin <= a && a <= xrmax) return 1;
    }
  }
  if (x1 != x2) {
    m=(y2-y1)/(x2-x1);
    if (xsmin <= xrmin && xrmin <= xsmax) {
      a = y1+(xrmin-x1)*m;
      if (yrmin <= a && a <= yrmax) return 1;
    }
    if (xsmin <= xrmax && xrmax <= xsmax) {
      a = y1+(xrmax-x1)*m;
      if (yrmin <= a && a <= yrmax) return 1;
    }
  }
  return 0;
} /* intersect */

void WriteAreasFile() {
  register int i;
  FILE *fp;
  char * pch;
  char outname[1000]="";
  pch=strrchr(file_name,'.');
  strncat (outname, file_name, pch-file_name);
  strcat(outname,".areas");
  fp=fopen(outname, "w");
  fprintf(fp,"# fiat 1.0\n");
  fprintf(fp,"# Cell areas \n");
  fprintf(fp,"# ttype1=area\n");
  fprintf(fp,"# Original file: %s\n",file_name);
  fprintf(fp,"# Redshift interval: ( %5.2f , %5.2f )\n", redshmin, redshmax);
  fprintf(fp,"# Center of frame: ( %12.8f , %12.8f, %12.8f )\n", 
	  xmean, ymean,zmean);
  for(i=0; i<ndata_ok; i++) fprintf(fp,"%12.8f\n",area[i]);
  fclose(fp);
  if(Verbose) printf("Areas file saved.\n");
} 

int Background(REAL *f, REAL area_tot, int ndata_ok) {
  REAL sum=0., area_mean=0., Nbackground=0.;
  REAL f1, f2, df, sum1, y0, yn, sumold;
  register int i, j, k;

  if(!Poisson) return ndata_ok;

  /* Ramella code */
  sumold=0.;
  for(i=ndata_ok/2;i<=ndata_ok;i++) {
    area_mean = area_tot/(REAL) i;
    sum=0.; j=1;
    f1=f[0]*area_mean;
    y0=TF(f1)*(REAL)i/(REAL) ndata_ok;
    while((f2=f[j]*area_mean) <= 0.8) {
      sum1=y0/2.;
      df=(f2-f1)/6.;
      for(k=0;k<6;k++) {
        yn=fabs(TF(f1)*(REAL)i-(REAL)(j-1))/(REAL) ndata_ok;
        sum1 += yn;
        f1+=df;
      }
      y0=fabs(TF(f1)*(REAL)i-(REAL)(j-1))/(REAL) ndata_ok;
      sum1 += y0*.5;
      sum += df*sum1;
      f1=f2;
      j++;
    }
    /* part revised by WB on 19/07/2000 */
    if (i == ndata_ok/2) sumold=sum; 
    if (i > ndata_ok/2){
      if (sumold > sum) Nbackground=(REAL) i;
    }
    sumold=sum;
  }
  /* Ebeling correction */
  Nbackground *= 1.025;
  if(Nbackground > ndata_ok) return ndata_ok;
  return (int)Nbackground;
}

REAL ComputeThreshold(REAL *f){
  int register i;
  int nresiduals=0;
  REAL *xresiduals, *yresiduals;
  REAL a1=0,a2=0,a3=0,tmp,tmpfmin=0,tmpfmax=3;

  if(Poisson) { // from Ramella's code
    tmpfmin=0.8;
    tmpfmax=2.6;
    }

  if(Ebeling) {
    ARRAY1(REAL, xresiduals, ndata_ok);
    ARRAY1(REAL, yresiduals, ndata_ok);
    for(i=0;i<ndata_ok;i++)
      if(f[i]>tmpfmin && f[i]<tmpfmax) {
	xresiduals[nresiduals]=log10(f[i]);
	yresiduals[nresiduals]=((REAL)(i+1)-TF(f[i])*Nbg)/(REAL) ndata_ok;
	nresiduals++;
      }
    threshold=FitParabol(nresiduals, xresiduals, yresiduals, &a1, &a2, &a3);
    free(xresiduals); free(yresiduals);
    if(Verbose){
      printf("Detection threshold set via Ebeling's method.\n");
      printf("Coefficients of the fit: ( %f , %f , %f) \n",a1,a2,a3);
    }
  } 
  else {
    if(Verbose) printf("Fixed detection threshold for scl = %f\n",conf_lev); 
    for(i=ndata_ok-1,tmp=1.; tmp>=conf_lev;i--) {
      tmp=TF(f[i]);
      threshold=f[i];
    }
  }
  if(Verbose) printf("Density threshold for detection: %f\n", threshold);
  return threshold;
}
 
REAL FitParabol(int ndata, REAL *x, REAL *y, REAL *a, REAL *b, REAL *c) {
  register int i, N=0;
  REAL xx, det,                 /* x coordinate */
    XXY,XY,Y,X4,X3,X2,X;    /* Coeff. of the l.m.s. fit */
/* Make a parabolic fit to the residual between empirical cumulative
   distribution function and the theoretical background function 
   and returns the x-coordinate of the minimum */
  XXY=XY=Y=X4=X3=X2=X=0.;
  for(i=0; i<ndata; i++) {
    xx=x[i];
    N++;
    X+=xx;
    XY+=xx*y[i];
    Y+=y[i];
    xx*=xx;
    X2+=xx; XXY+=xx*y[i]; X3+=xx*x[i];
    X4+=xx*xx;
  }
  det=X4*(X2*N-X*X)+X3*(X2*X-X3*N)+X2*(X3*X-X2*X2);
  *a=XXY*(X2*N-X*X)+X3*(X*Y-XY*N)+X2*(XY*X-X2*Y);
  *b=X4*(XY*N-X*Y)+XXY*(X*X2-X3*N)+X2*(X3*Y-XY*X2);
  *c=X4*(X2*Y-XY*X)+X3*(XY*X2-X3*Y)+XXY*(X3*X-X2*X2);
  *a /= det; *b /= det; *c /=det;
  return pow(10,-0.5*(X4*(XY*N-X*Y)+XXY*(X*X2-X3*N)+X2*(X3*Y-XY*X2))/
	     (XXY*(X2*N-X*X)+X3*(X*Y-XY*N)+X2*(XY*X-X2*Y)));
} /* FitParabol() */

void WriteDistributionsFile(REAL *f) {
  register int i;
  REAL ff;
  FILE *fp;
  char * pch;
  char outname[1000]="";
  pch=strrchr(file_name,'.');
  strncat (outname, file_name, pch-file_name);
  strcat(outname,".distributions");
  fp=fopen(outname, "w");
  fprintf(fp,"# fiat 1.0\n");
  fprintf(fp,"# Distribution of normalized cell densities as well as the \n");
  fprintf(fp,"# background distribution for comparison.\n");  
  fprintf(fp,"# ttype1=delta\n# ttype2=pdf\n# ttype3=pdf_bck\n");
  fprintf(fp,"# Original file: %s\n",file_name);
  fprintf(fp,"# Redshift interval: ( %5.2f , %5.2f )\n", redshmin, redshmax);
  fprintf(fp,"# Center of frame: ( %12.8f , %12.8f , %12.8f )\n", 
	  xmean, ymean,zmean);
  for(i=0; i<ndata_ok; i++) 
    fprintf(fp,"%12.8f %12.8f %12.8f\n",f[i],(REAL) (i+1)/(REAL) ndata_ok,
	    TF(f[i])*(REAL) Nbg/(REAL) ndata_ok);
  fclose(fp);
  if(Verbose) printf("Distributions file saved.\n");
} /* WriteDistributionsFile() */

void WriteCandidatesFile(int *prox, int *mask) {
  register int i,j,k;
  int *list=(int *)NULL, *tmp_list=(int *)NULL, nprox, nl=0,nl0=0, jj,kk,cc=0;
  REAL *sort_area=(REAL *) NULL, area_cluster=0., nback=0, flux_min=0.;
  FILE *fp;
  char * pch;
  char outname[1000]="";
  pch=strrchr(file_name,'.');
  strncat (outname, file_name, pch-file_name);
  strcat(outname,".raw");
  fp=fopen(outname, "w");
  fprintf(fp,"# fiat 1.0\n");
  fprintf(fp,"# Dumped list of candidates. \n");
  fprintf(fp,"# ttype1=n\n# ttype2=flux_min\n# ttype3=Ncorr\n");
  fprintf(fp,"# ttype4=Nlinked\n# ttype5=Nbck\n");
  fprintf(fp,"# Original file: %s\n",file_name);
  fprintf(fp,"# Redshift interval: ( %5.2f , %5.2f )\n", redshmin, redshmax);
  fprintf(fp,"# Center of frame: ( %12.8f , %12.8f  , %12.8f )\n", xmean, 
	  ymean,zmean);
  for(i=0; i<ndata; i++) {
    if(mask[i] == 1) {
      nl=1; ADDLIST(nl,i);
      do { /* while(nprox>0) */
        for(j=0;j<nl;j++) {
          jj=list[j];
          if (mask[jj] != -1) {
            for(k=prox[jj];k<prox[jj+1];k++) {
              kk=prox[k];
              if(mask[kk] == 1) { nl++; ADDLIST(nl,kk); mask[kk]=2; }
            }
            mask[jj]=-1;
          }
        }
        nprox=0;
        for(j=0;j<nl;j++) if(mask[list[j]]==2) nprox++;
      } while (nprox > 0);
      area_cluster=0;
      for(j=0;j<nl;j++) area_cluster+=area[list[j]]; 
      sort_area=(REAL *) realloc(sort_area,nl*sizeof(REAL));
      for(j=0;j<nl;j++) sort_area[j]=area[list[j]];
      qsort((REAL *) sort_area, nl, sizeof(REAL), real_compare);
      nl0=nl;
      nback=Nbg*area_cluster/tot_area_covered;
      flux_min=area_mean/sort_area[nl-1];
      fprintf(fp,"%d %f %d %d %d\n",
	      ++cc,flux_min,nl0-(int)nback,nl0,(int)nback);
    } /* if(mask[i]==1) */
  } /* for(i=0;i<ndata;i++) */
  free(list);
  free(tmp_list);
  free(sort_area);
  fclose(fp);
  if (Verbose) printf("Percolation done.\nNumber of candidates: %d\n",cc);
  if(Verbose) printf("Preliminary list of candidates file saved.\n");
} /* WriteCandidatesFile() */

double erf(double x){
return    1.1283791670955126*x - 
  0.3761263890318375*pow(x,3) + 
  0.11283791670955126*pow(x,5) - 
  0.02686617064513125*pow(x,7) + 
  0.005223977625442187*pow(x,9) - 
  0.0008548327023450852*pow(x,11) + 
  0.00012055332981789664*pow(x,13) - 
  0.000014925650358406252*pow(x,15) + 
  1.6462114365889246e-6*pow(x,17) - 
  1.6365844691234924e-7*pow(x,19) + 
  1.4807192815879218e-8*pow(x,21) - 
  1.2290555301717926e-9*pow(x,23) + 
  9.422759064650411e-11*pow(x,25) - 
  6.7113668551641105e-12*pow(x,27) + 
  4.463224263286477e-13*pow(x,29) - 
  2.7835162072109215e-14*pow(x,31) + 
  1.6342614095367152e-15*pow(x,33) - 
  9.063970842808672e-17*pow(x,35) + 
  4.763348040515068e-18*pow(x,37) - 
  2.3784598852774293e-19*pow(x,39) + 
  1.131218725924631e-20*pow(x,41) - 
  5.136209054585811e-22*pow(x,43) + 
  2.230878680274645e-23*pow(x,45) - 
  9.28672901131906e-25*pow(x,47) + 
  3.71153285316323e-26*pow(x,49) - 
  1.4263930180784176e-27*pow(x,51) + 
  5.279103332510835e-29*pow(x,53) - 
  1.8841244217042036e-30*pow(x,55) + 
  6.492909974544561e-32*pow(x,57) - 
  2.1630383901171243e-33*pow(x,59) + 
  6.973730328792914e-35*pow(x,61) - 
  2.1781748594796097e-36*pow(x,63) + 
  6.597356545539202e-38*pow(x,65) - 
  1.9395213725013484e-39*pow(x,67) + 
  5.539127534424142e-41*pow(x,69) - 
  1.538027363683162e-42*pow(x,71) + 
  4.1552489658106737e-44*pow(x,73) - 
  1.0930925207357809e-45*pow(x,75) + 
  2.8018434400267797e-47*pow(x,77) - 
  7.002335114640116e-49*pow(x,79) + 
  1.7073594878289173e-50*pow(x,81) - 
  4.0639470618319806e-52*pow(x,83) + 
  9.448392328628974e-54*pow(x,85) - 
  2.1467878854142287e-55*pow(x,87) + 
  4.769421502324767e-57*pow(x,89) ;
}

Connect *percola(int *prox, int *mask) {
  register int i,j,k;
  Connect *connection=NULL;
  int *list=(int *)NULL, *tmp_list=(int *)NULL, nprox, nl,nl0, jj, kk, nlist;
  REAL *sort_area=(REAL *) NULL, area_cluster=0., nback=0, flux_min=0.,
    b,Nfluct,Nfluct0,tmp=TF(threshold);
  double erf(double x);

  for(i=0; i<ndata; i++) {
    if(mask[i] == 1) {
      nl=1; ADDLIST(nl,i);
      do { /* while(nprox>0) */
        for(j=0;j<nl;j++) {
          jj=list[j];
          if (mask[jj] != -1) {
            for(k=prox[jj];k<prox[jj+1];k++) {
              kk=prox[k];
              if(mask[kk] == 1) { nl++; ADDLIST(nl,kk); mask[kk]=2; }
            }
            mask[jj]=-1;
          }
        }
        nprox=0;
        for(j=0;j<nl;j++) if(mask[list[j]]==2) nprox++;
      } while (nprox > 0);
      /* Check if fluctuation is significative */
      area_cluster=0;
      for(j=0;j<nl;j++) area_cluster+=area[list[j]]; 
      sort_area=(REAL *) realloc(sort_area,nl*sizeof(REAL));
      for(j=0;j<nl;j++) sort_area[j]=area[list[j]];
      qsort((REAL *) sort_area, nl, sizeof(REAL), real_compare);
      nl0=nl;
      nback=Nbg*area_cluster/tot_area_covered;
      while((REAL)nl0>0.){
  	flux_min=area_mean/sort_area[nl0-1];
	//if (flux_min > 2.2) flux_min=2.2; 
	//Nfluct0=exp(c0*flux_min+c1);
	//b=c2*flux_min+c3;
	//if (b<0) b=0.;
	//if (Nfluct0<0.) Nfluct0=0.;
	//Nfluct=exp(-b*nl0)*Nfluct0;
	Nfluct0=((flux_min/threshold)-1.)*nl0; // signal/noise  
	if (Nfluct0 > 7.) Nfluct0 = 7.; 
	Nfluct=erf(Nfluct0/sqrt(2.0));
	if (Nfluct < 0.) Nfluct = 1.;
	if (Nfluct > 1.) Nfluct = 1.;
	if ((Nfluct >= reject_lev)){  
	  if(nl != nl0) {
            tmp_list=(int *) realloc(tmp_list,nl*sizeof(int));
            for(k=0;k<nl;k++) tmp_list[k]=list[k];
            list=(int *) realloc(list,nl0*sizeof(int));
            nlist=0;
	    for(k=0;k<nl;k++)
	      if(area[tmp_list[k]] < sort_area[nl0])
		list[nlist++]=tmp_list[k];
          }
          connection=AddItem(connection,nl0,Nfluct,list);
          break;
	}
        nl0=nl0-1;
        area_cluster-=sort_area[nl0];
        nback=Nbg*area_cluster/tot_area_covered;
      } /*while()*/
    } /* if(mask[i]==1) */
  } /* for(i=0;i<ndata;i++) */
  free(list);
  free(tmp_list);
  free(sort_area);
  if(Verbose) printf("Percolation done with rejection. cl: %lf\n", reject_lev);
  return connection;
} /* percola() */

Connect *AddItem (Connect *listpointer, int n, REAL sig, int *list) {
  Connect *lp=listpointer;
  register int i;
  if (listpointer != NULL) {
    while (listpointer->link != NULL)
      listpointer=listpointer->link;
    listpointer->link=(Connect *) malloc(sizeof(Connect));
    listpointer=listpointer->link;
    listpointer->link=NULL;
    listpointer->n=n;
    listpointer->sig=sig;
    listpointer->list=(int *) malloc(sizeof(int)*n);
    for(i=0;i<n;i++)
      listpointer->list[i]=list[i];
    return lp;
  } else {
    listpointer=(Connect *) malloc(sizeof(Connect));
    listpointer->link=NULL;
    listpointer->n=n;
    listpointer->sig=sig;
    listpointer->list=(int *) malloc(sizeof(int)*n);
    for(i=0;i<n;i++)
      listpointer->list[i]=list[i];
    return listpointer;
  }
} /* AddItem()*/



int CountItem (Connect *listpointer) {
  register int i=0;
  while (listpointer != NULL) i++, listpointer=listpointer->link;
  if(Verbose) printf("Number of clusters detected: %d\n",i);
  if(i == 0) printf("Warning: No clusters were found.\n");
  return i;
} /* CountItem() */

void WriteVoronoiCatalog() {
  register int i;
  REAL x1,x2,y1,y2;
  int i1,i2;
  FILE *fp;
  char * pch;
  char outname[1000]="";
  pch=strrchr(file_name,'.');
  strncat (outname, file_name, pch-file_name);
  strcat(outname,".voronoi");
  fp=fopen(outname, "w");
  fprintf(fp,"# fiat 1.0\n");
  fprintf(fp,"# List of Voronoi edges. These are line segments linking \n");
  fprintf(fp,"# the points (x1,y1) to (x2,y2).\n");  
  fprintf(fp,"# ttype1=x1\n# ttype2=y1\n# ttype3=x2\n# ttype3=y2\n");
  fprintf(fp,"# Original file: %s\n",file_name);
  fprintf(fp,"# Redshift interval: ( %5.2f , %5.2f )\n", redshmin, redshmax);
  fprintf(fp,"# Center of frame: ( %12.8f , %12.8f  , %12.8f )\n", xmean, 
	  ymean,zmean);
  fprintf(fp,"#    x1            y1          x2            y2\n");
  for(i=0;i<vor.numberofedges;i++) {
    if(voronoi_inside[i]==1) {
      i1=vor.edgelist[2*i]; i2=vor.edgelist[2*i+1];
      x1=vor.pointlist[2*i1]; y1=vor.pointlist[2*i1+1];
      x2=vor.pointlist[2*i2]; y2=vor.pointlist[2*i2+1];
      if (!XYcoord)
        fprintf(fp,"%12.8f %12.8f %12.8f %12.8f\n",
                xmean+x1/cos((y1+ymean)*M_PI/180.), ymean+y1,
                xmean+x2/cos((y2+ymean)*M_PI/180.), ymean+y2);
      else
        fprintf(fp,"%12.8f %12.8f %12.8f %12.8f\n",
                xmean+x1, ymean+y1, xmean+x2, ymean+y2);
    }
  }
  fclose(fp);
  if(Verbose) printf("Voronoi diagram saved.\n");
} /* WriteVoronoiCatalog() */

void WriteDelaunayCatalog() {
  register int i;
  int i1,i2;
  FILE *fp;
  char * pch;
  char outname[1000]="";
  pch=strrchr(file_name,'.');
  strncat (outname, file_name, pch-file_name);
  strcat(outname,".delaunay");
  fp=fopen(outname, "w");
  fprintf(fp,"# fiat 1.0\n");
  fprintf(fp,"# List of Delaunay triangle edges. These are line segments\n");  
  fprintf(fp,"# linking each point (x=i1,y=i2) to the next.\n"); 
  fprintf(fp,"# ttype1=x1\n# ttype2=y1\n# ttype3=x2\n# ttype3=y2\n");
  fprintf(fp,"# Original file: %s\n",file_name);
  fprintf(fp,"# Redshift interval: ( %5.2f , %5.2f )\n", redshmin, redshmax);
  fprintf(fp,"# Center of frame: ( %12.8f , %12.8f  , %12.8f )\n", xmean, 
	  ymean,zmean);
  fprintf(fp,"#    i1            i2\n");
  for(i=0;i<mid.numberofedges;i++) {
    i1=mid.edgelist[2*i]; i2=mid.edgelist[2*i+1];
    fprintf(fp,"%d %d\n", i1, i2);
  }
  fclose(fp);
  if(Verbose) printf("Delaunay triangles saved.\n");
} /* WriteDelaunayCatalog() */

void ComputeRedshifts(Connect *conn) {
  register int i, ii;
  REAL z_mean,*zmemb;
  Connect *p_conn=conn;
  Footprint *p_cir=footprint;
  zmemb=(REAL *) NULL;
  while(p_conn != NULL) {
    zmemb = (REAL *) realloc(zmemb,p_conn->n*sizeof(REAL));
    z_mean=0.;
    for(i=0; i<p_conn->n;i++) {
      ii=p_conn->list[i];
      zmemb[i]=redsh[ii];
      z_mean+=redsh[ii];      
    }
    z_mean/=p_conn->n;
    p_cir->z_mean=z_mean;
    if (p_conn->n<3){
      p_cir->z_median=z_mean;
      p_cir->z_rms=(redshmax-redshmin)*0.5;
    } else {
      p_cir->z_median=median(zmemb,p_conn->n);
      p_cir->z_rms=rms(zmemb,p_conn->n,z_mean);
    }
    /* go to next cluter */
    p_conn=p_conn->link;
    p_cir++;
  } 
  free(zmemb);
  if(Verbose) printf("Cluster redshifts computed.\n");
} /* ComputeRedshifts() */

REAL rms(REAL *data,int n,REAL ave){
  int j;
  REAL s=0.,var=0.,ep=0.,adev=0.;
  if(n<2) return 0;
  for(j=0;j<n;j++){
    adev+=abs(s=data[j]-ave);
    ep+=s;
    var+=s*s;
  }
  adev/=n;
  var=(var-ep*ep/n)/(n-1);
  return sqrt(var);
}

void WriteSimpleMembersList() {
  register int i;
  long long hid;
  Footprint *p_cir=footprint;
  FILE *fp;
  char * pch;
  char outname[1000]="";
  pch=strrchr(file_name,'.');
  strncat (outname, file_name, pch-file_name);
  strcat(outname,".VTmembers.list");
  fp=fopen(outname, "w");
  for(i=0;i<ndata;i++)
    if(edge1[i]==0) fprintf(fp, "%lld %lld \n", number[i], hid=mask1[i]);
  fclose(fp);
  if(Verbose) printf("Members list saved.\n");
} /* WriteSimpleMembersList() */

void WriteGalaxyCatalogAscii() {
  register int i;
  Footprint *p_cir=footprint;
  FILE *fp;
  char * pch;
  char outname[1000]="";
  pch=strrchr(file_name,'.');
  strncat (outname, file_name, pch-file_name);
  strcat(outname,".VTgalaxies.cat");
  fp=fopen(outname, "w");
  fprintf(fp,"# fiat 1.0\n");
  fprintf(fp,"# VT cluster members list. \n");
  fprintf(fp,"# ttype1=id\n# ttype2=ra\n# ttype3=dec\n# ttype4=z\n");
  fprintf(fp,"# ttype5=host_id\n# ttype6=local_dens\n# ttype7=central\n");
  fprintf(fp,"# ttype8=x_c\n# ttype9=y_c\n");
  fprintf(fp,"# Original file: %s\n",file_name);
  fprintf(fp,"# Redshift interval: ( %5.2f , %5.2f )\n", redshmin, redshmax);
  fprintf(fp,"# Center of frame: ( %12.8f , %12.8f  , %12.8f )\n", xmean, 
	  ymean,zmean);
  fprintf(fp,"# Summary of parameters and flags:\n");
  fprintf(fp,"# name       value\n");
  fprintf(fp,"# A          %f\n",w_amp);
  fprintf(fp,"# g          %f\n",w_pow);
  fprintf(fp,"# r          %f\n",reject_lev);
  if(!Ebeling) fprintf(fp,"# s          %f\n",conf_lev);
  fprintf(fp,"# C          %s %s %s %s\n",ttype[0],ttype[1],ttype[2],ttype[3]);
  fprintf(fp,"# B          %s %s %s %s\n",
	  ttypeH[0],ttypeH[1],ttypeH[2],ttypeH[3]);
  fprintf(fp,"# F          %f %f %f %f\n",frame[0],frame[1],frame[2],frame[3]);
  if(Holes) fprintf(fp,"# holesfile  %s\n",file_holes);
  fprintf(fp,"# datafile   %s\n",file_name);
  fprintf(fp,"# Verbose    %d\n",Verbose);
  fprintf(fp,"# Holes      %d\n",Holes);
  fprintf(fp,"# Ebeling    %d\n",Ebeling);
  fprintf(fp,"# Poisson    %d\n",Poisson);
  fprintf(fp,"# Fits       %d\n",Fits);
  fprintf(fp,"# XYcoord    %d\n",XYcoord);
  fprintf(fp,"# NotFind    %d\n",NotFind);
  fprintf(fp,"# RawFind    %d\n",RawFind);

  if (!XYcoord) 
    for(i=0;i<ndata;i++) if(edge1[i]==0)
      fprintf(fp, "%lld %lf %lf %lf %lld %lf %d %lf %lf\n",
              number[i], in.pointlist[2*i], in.pointlist[2*i+1], redsh[i],
	      p_cir->id, 1/area[i], central[i], xmean+
	      in.pointlist[2*i]/cos((in.pointlist[2*i+1]+ymean)*M_PI/180.),
              in.pointlist[2*i+1]+ymean);
  else 
    for(i=0;i<ndata;i++) if(edge1[i]==0)
      fprintf(fp, "%lld %lf %lf %lf %lld %lf %d %lf %lf\n",
              number[i], in.pointlist[2*i], in.pointlist[2*i+1], redsh[i],
	      p_cir->id, 1/area[i], central[i], xmean+in.pointlist[2*i],
	      in.pointlist[2*i+1]+ymean);
  
  fclose(fp);
  if(Verbose) printf("Galaxy catalog saved. Ascii only.\n");
} /* WriteGalaxyCatalogAscii() */

void WriteClusterCatalogAscii() {  
  register int i;  
  FILE *fp;  
  Footprint *p_cir=footprint;  
  char * pch;  
  char outname[1000]="";  
  pch=strrchr(file_name,'.');  
  strncat (outname, file_name, pch-file_name);  
  strcat(outname,".VTclusters.cat");  
  fp=fopen(outname, "w");  
  fprintf(fp,"# fiat 1.0\n");  
  fprintf(fp,"# VT cluster catalog. \n");  
  fprintf(fp,"# ttype1=id\n# ttype2=ra\n# ttype3=dec\n# ttype4=z_median\n"); 
  fprintf(fp,"# ttype5=z_mean\n# ttype6=z_rms\n# ttype7=contrast\n");  
  fprintf(fp,"# ttype8=sig\n# ttype9=nvt\n# ttype10=nbg\n");  
  fprintf(fp,"# ttype11=ncorr\n# ttype12=size\n# ttype13=a\n# ttype14=b\n"); 
  fprintf(fp,"# ttype15=pa\n# ttype16=area\n");  
  fprintf(fp,"# ttype17=min_density\n# ttype18=ra_c\n# ttype19=dec_c\n");  
  fprintf(fp,"# ttype20=z_c\n# ttype21=x_0\n# ttype22=y_0\n");  
  fprintf(fp,"# Original file: %s\n",file_name);  
  fprintf(fp,"# Redshift interval: ( %5.2f , %5.2f )\n", redshmin, redshmax);
  fprintf(fp,"# Center of frame: ( %12.8f , %12.8f , %12.8f )\n",   
  	  xmean, ymean,zmean);  
  fprintf(fp,"# Summary of parameters and flags:\n");  
  fprintf(fp,"# name       value\n");  
  fprintf(fp,"# A          %f\n",w_amp);  
  fprintf(fp,"# g          %f\n",w_pow);  
  fprintf(fp,"# r          %f\n",reject_lev);  
  if(!Ebeling) fprintf(fp,"# s          %f\n",conf_lev);  
  fprintf(fp,"# C          %s %s %s %s\n",
	  ttype[0],ttype[1],ttype[2],ttype[3]);  
  fprintf(fp,"# B          %s %s %s %s\n",  
	  ttypeH[0],ttypeH[1],ttypeH[2],ttypeH[3]);  
  fprintf(fp,"# F          %f %f %f %f\n",
	  frame[0],frame[1],frame[2],frame[3]);  
  if(Holes) fprintf(fp,"# holesfile  %s\n",file_holes);  
  fprintf(fp,"# datafile   %s\n",file_name);  
  fprintf(fp,"# Verbose    %d\n",Verbose);  
  fprintf(fp,"# Holes      %d\n",Holes);  
  fprintf(fp,"# Ebeling    %d\n",Ebeling);  
  fprintf(fp,"# Poisson    %d\n",Poisson);  
  fprintf(fp,"# Fits       %d\n",Fits);  
  fprintf(fp,"# XYcoord    %d\n",XYcoord);  
  fprintf(fp,"# NotFind    %d\n",NotFind);  
  fprintf(fp,"# RawFind    %d\n",RawFind);  
  
  if (!XYcoord)   
    for(i=0; i<nclusters; i++, p_cir++) if(p_cir->edge==0){  
      fprintf(fp, "%lld %lf %lf %lf %lf %lf %lf %lf %d %d %d %lf %lf %lf "  
  	      "%lf %lf %lf %lf %lf %lf %lf %lf",  
  	      p_cir->id,       
  	      xmean+p_cir->xc/cos((p_cir->yc+ymean)*M_PI/180.),  
  	      p_cir->yc+ymean,                 
  	      p_cir->z_median, p_cir->z_mean,    p_cir->z_rms,      
  	      p_cir->contrast, p_cir->sig,  
	      p_cir->nvor,     p_cir->nbg,       p_cir->n,    p_cir->radius, 
  	      p_cir->semi_maj, p_cir->semi_min,  p_cir->pa,     
	      p_cir->area,     p_cir->minflux,       
  	      xmean+p_cir->xgc/cos((p_cir->yc+ymean)*M_PI/180.),  
  	      ymean+p_cir->ygc,  
  	      p_cir->zgc,      p_cir->xc,     p_cir->yc);  
      }  
  else   
    for(i=0; i<nclusters; i++, p_cir++) if(p_cir->edge==0){  
      fprintf(fp, "%lld %lf %lf %lf %lf %lf %lf %lf %d %d %d %lf %lf %lf "  
  	      "%lf %lf %lf %lf %lf %lf %lf %lf",  
  	      p_cir->id,       
  	      xmean+p_cir->xc,  
  	      p_cir->yc+ymean,                 
  	      p_cir->z_median, p_cir->z_mean,    p_cir->z_rms,      
  	      p_cir->contrast, p_cir->sig,  
	      p_cir->nvor,     p_cir->nbg,       p_cir->n,    p_cir->radius, 
  	      p_cir->semi_maj, p_cir->semi_min,  p_cir->pa, 
	      p_cir->area,     p_cir->minflux,       
  	      xmean+p_cir->xgc,  
  	      ymean+p_cir->ygc,  
  	      p_cir->zgc,      p_cir->xc,     p_cir->yc);  
    }  
  fclose(fp);  
  if(Verbose) printf("Cluster catalog saved. Ascii only.\n");
} /* WriteClusterCatalogAscii() */ 

void WriteGalaxyCatalog() {
  register int i;
  Footprint *p_cir=footprint;
  FILE *fp;
  fitsfile *fptr;
  int status=0, hdutype;
  long firstrow, firstelem;
  int tfields   = 9;       
  char extname[] = "VTgalaxies";
  /* col names */
  char *ttype[] = {"id",          "host_id",  "ra",   "dec",      "z",     
		   "local_dens",  "x_0",      "y_0",  "central",             };
  /* data types */
  char *tform[] = { "1K",         "1K",       "1D",   "1D",      "1D",
		    "1D",         "1D",       "1D",   "1I"                   };
  /* units */
  char *tunit[] = { "\0",         "\0",       "deg",  "deg",      "\0",    
		    "\0",         "\0",       "\0",   "\0",       "\0"       };
  if(XYcoord) for(i-0;i<tfields;i++) tunit[i] = "\0";
  /* define some of the col entries  and row vars */
  long long host_id,frow=1,felem=1,nelem=1, colnum=1, **llrow;
  REAL tmp=1.,**Rrow;
  int **irow,j;
  ARRAY2(REAL,Rrow,6,1);
  ARRAY2(long long,llrow,2,1);
  ARRAY2(int,irow,1,1);
  /* open files */
  char * pch;
  char outnameA[1000]="";
  pch=strrchr(file_name,'.');
  strncat (outnameA, file_name, pch-file_name);
  strcat(outnameA,"_VTgalaxies.cat");
  fp=fopen(outnameA, "w");
  char outnameF[1000]="";
  strncat (outnameF, file_name, pch-file_name);
  strcat(outnameF,"_VTgalaxies.fit");
  fits_create_file(&fptr,outnameF,&status);
  fits_create_tbl(fptr,BINARY_TBL,0,tfields,ttype,tform,tunit,extname,&status);
  /* write */
  fprintf(fp,"# fiat 1.0\n");
  fprintf(fp,"# VT cluster members list. \n");
  fprintf(fp,"# ttype1=id\n# ttype2=host_id\n# ttype3=ra\n# ttype4=dec\n");
  fprintf(fp,"# ttype5=z\n# ttype6=local_dens\n# ttype7=x_0\n");
  fprintf(fp,"# ttype8=y_0\n# ttype9=central\n");
  fprintf(fp,"# Original file: %s\n",file_name);
  fprintf(fp,"# Redshift interval: ( %5.2f , %5.2f )\n", redshmin, redshmax);
  fprintf(fp,"# Center of frame: ( %12.8f , %12.8f  , %12.8f )\n", xmean, 
	  ymean,zmean);
  fprintf(fp,"# Summary of parameters and flags:\n");
  fprintf(fp,"# name       value\n");
  fprintf(fp,"# A          %f\n",w_amp);
  fprintf(fp,"# g          %f\n",w_pow);
  fprintf(fp,"# r          %f\n",reject_lev);
  if(!Ebeling) fprintf(fp,"# s          %f\n",conf_lev);
  fprintf(fp,"# C          %s %s %s %s\n",ttype[0],ttype[1],ttype[2],ttype[3]);
  fprintf(fp,"# B          %s %s %s %s\n",
	  ttypeH[0],ttypeH[1],ttypeH[2],ttypeH[3]);
  fprintf(fp,"# F          %f %f %f %f\n",frame[0],frame[1],frame[2],frame[3]);
  if(Holes) fprintf(fp,"# holesfile  %s\n",file_holes);
  fprintf(fp,"# datafile   %s\n",file_name);
  fprintf(fp,"# Verbose    %d\n",Verbose);
  fprintf(fp,"# Holes      %d\n",Holes);
  fprintf(fp,"# Ebeling    %d\n",Ebeling);
  fprintf(fp,"# Poisson    %d\n",Poisson);
  fprintf(fp,"# Fits       %d\n",Fits);
  fprintf(fp,"# XYcoord    %d\n",XYcoord);
  fprintf(fp,"# NotFind    %d\n",NotFind);
  fprintf(fp,"# RawFind    %d\n",RawFind);
  for(i=0;i<ndata;i++) if(edge1[i]==0) {
    tmp=1.;
    if (!XYcoord) tmp=cos((in.pointlist[2*i+1]+ymean)*M_PI/180.);
    fprintf(fp, "%lld %lld %lf %lf %lf %lf %lf %lf %d\n",
	    llrow[0][0]=number[i],
	    llrow[1][0]=mask1[i], 
	    Rrow[0][0]=xmean+in.pointlist[2*i]/tmp,
	    Rrow[1][0]=in.pointlist[2*i+1]+ymean,
	    Rrow[2][0]=redsh[i],
	    Rrow[3][0]=1/area[i],  
	    Rrow[4][0]=in.pointlist[2*i], 
	    Rrow[5][0]=in.pointlist[2*i+1], 
	    irow[0][0]=central[i]);
    fits_write_col_lnglng(fptr, colnum++,frow,felem,nelem,llrow[0],&status);
    fits_write_col_lnglng(fptr, colnum++,frow,felem,nelem,llrow[1],&status);
    for(j=0;j<6;j++)
      fits_write_col_dbl(fptr,colnum++,frow,felem,nelem,Rrow[j],&status);
    fits_write_col_int(fptr,colnum++,frow,felem,nelem,irow[0],&status);
    cfitsio_error(status);
    colnum=1;
    frow++;
  }
  /*
  char comment[70]="";
  sprintf(comment,"VT cluster catalog.");
  fits_write_comment(fptr,comment,&status);
  sprintf(comment,"Original file: %s",file_name);
  fits_write_comment(fptr,comment,&status);
  sprintf(comment,"Redshift interval: ( %5.2f , %5.2f )", 
	  redshmin, redshmax);
  fits_write_comment(fptr,comment,&status);
  sprintf(comment,"Center of frame: ( %12.8f , %12.8f , %12.8f )", 
	  xmean, ymean,zmean);
  fits_write_comment(fptr,comment,&status);
  sprintf(comment,"Summary of parameters and flags:");
  fits_write_comment(fptr,comment,&status);
  sprintf(comment,"name       value");
  fits_write_comment(fptr,comment,&status);
  sprintf(comment,"A          %f",w_amp);
  sprintf(comment,"g          %f\n",w_pow);
  fits_write_comment(fptr,comment,&status);
  sprintf(comment,"r          %f\n",reject_lev);
  fits_write_comment(fptr,comment,&status);
  if(!Ebeling) {
    sprintf(comment,"s          %f\n",conf_lev);
    fits_write_comment(fptr,comment,&status);
  }
  sprintf(comment,"C          %s %s %s %s\n",
	  ttype[0],ttype[1],ttype[2],ttype[3]);
  fits_write_comment(fptr,comment,&status);
  sprintf(comment,"B          %s %s %s %s\n",
	  ttypeH[0],ttypeH[1],ttypeH[2],ttypeH[3]);
  fits_write_comment(fptr,comment,&status);
  sprintf(comment,"F          %f %f %f %f\n",
	  frame[0],frame[1],frame[2],frame[3]);
  if(Holes) {
    sprintf(comment,"holesfile  %s\n",file_holes);
    fits_write_comment(fptr,comment,&status);
  }  
  sprintf(comment,"datafile   %s\n",file_name);
  fits_write_comment(fptr,comment,&status);
  sprintf(comment,"Verbose    %d\n",Verbose);
  fits_write_comment(fptr,comment,&status);
  sprintf(comment,"Holes      %d\n",Holes);
  fits_write_comment(fptr,comment,&status);
  sprintf(comment,"Ebeling    %d\n",Ebeling);
  fits_write_comment(fptr,comment,&status);
  sprintf(comment,"Poisson    %d\n",Poisson);
  fits_write_comment(fptr,comment,&status);
  sprintf(comment,"Fits       %d\n",Fits);
  fits_write_comment(fptr,comment,&status);
  sprintf(comment,"XYcoord    %d\n",XYcoord);
  sprintf(comment,"NotFind    %d\n",NotFind);
  fits_write_comment(fptr,comment,&status);
  sprintf(comment,"RawFind    %d\n",RawFind);
  fits_write_comment(fptr,comment,&status);
  */
  free2(Rrow);
  free2(llrow);
  free2(irow);
  fclose(fp);
  fits_close_file(fptr, &status);
  if(Verbose) printf("Galaxy catalog saved.\n");
} /* WriteGalaxyCatalog() */

void WriteClusterCatalog(){
  register int i;
  Footprint *p_cir=footprint;
  FILE *fp;
  fitsfile *fptr;
  int status=0, hdutype;
  long firstrow, firstelem;
  int tfields   = 24;       
  char extname[] = "VTclusters";
  char *ttype[] = {"id"      ,       
		   "ra"      , "dec" ,  "z"       , "z_mean", "z_rms", 
		   "contrast", "sig" ,  "size"    , "semi_a", "semi_b"    ,       
		   "PA"      , "area",  "min_dens", "ra_c"  , "dec_c",
		   "z_c"     , "x_0" ,  "y_0"     , "A"     , "g"    ,     
		   "nvt"     , "nbg" ,  "ncorr"   };
  char *tform[] = { "1K",        
		    "1D",      "1D",    "1D",        "1D",      "1D",
		    "1D",      "1D",    "1D",        "1D",      "1D",     
		    "1D",      "1D",    "1D",        "1D",      "1D",   
		    "1D",      "1D",    "1D",        "1D",      "1D", 
		    "1I",      "1I",    "1I"};        
  char *tunit[] = { "\0",        
		    "deg",     "deg",   "\0",        "\0",       "\0",
		    "\0",      "\0",    "deg",       "deg",      "deg",    
		    "deg",     "deg2",  "1/deg",     "deg",      "deg",    
		    "\0",      "\0",    "\0",        "\0",      "\0" , 
		    "\0",      "\0",    "\0"};
  if(XYcoord) for(i-0;i<tfields;i++) tunit[i] = "\0";
  /* define the col entries */
  long long id;
  int nvt, nbg, ncorr; 
  REAL ra, dec, z, z_mean, z_median, z_rms, contrast, sig, size, 
    a, b, pa, areac, min_dens, ra_c, dec_c, z_c, x_0, y_0;
  /* and row vars */
  REAL tmp=1.,tmp1=1.,tmp2=1.,**Rrow;
  long long naxis2=0,frow=1,felem=1,nelem=1, colnum=1, **llrow;
  int **irow,j;
  ARRAY2(REAL,Rrow,20,1);
  ARRAY2(long long,llrow,1,1);
  ARRAY2(int,irow,3,1);
  /* open files */
  char * pch;
  char outnameA[1000]="";
  pch=strrchr(file_name,'.');
  strncat (outnameA, file_name, pch-file_name);
  strcat(outnameA,"_VTclusters.cat");
  fp=fopen(outnameA, "w");
  char outnameF[1000]="";
  strncat (outnameF, file_name, pch-file_name);
  strcat(outnameF,"_VTclusters.fit");
  fits_create_file(&fptr,outnameF,&status);
  fits_create_tbl(fptr,BINARY_TBL,naxis2,tfields,ttype,tform,tunit,extname,&status);
  /* write the output files */
  fprintf(fp,"# fiat 1.0\n");
  fprintf(fp,"# VT cluster catalog. \n");
  fprintf(fp,"# ttype1=id\n# ttype2=ra\n# ttype3=dec\n# ttype4=z\n");
  fprintf(fp,"# ttype5=z_mean\n# ttype6=z_rms\n# ttype7=contrast\n");
  fprintf(fp,"# ttype8=sig\n# ttype9=size\n# ttype10=semi_a\n");
  fprintf(fp,"# ttype11=semi_b\n# ttype12=PA\n# ttype13=area\n");
  fprintf(fp,"# ttype14=min_dens\n");
  fprintf(fp,"# ttype15=ra_c\n# ttype16=dec_c\n# ttype17=z_c\n");
  fprintf(fp,"# ttype18=x_0\n# ttype19=y_0\n# ttype20=A\n");
  fprintf(fp,"# ttype21=g\n# ttype22=nvt\n");
  fprintf(fp,"# ttype23=nbg\n# ttype24=ncorr\n");
  fprintf(fp,"# Original file: %s\n",file_name);
  fprintf(fp,"# Redshift interval: ( %5.2f , %5.2f )\n", redshmin, redshmax);
  fprintf(fp,"# Center of frame: ( %12.8f , %12.8f , %12.8f )\n", 
	  xmean, ymean,zmean);
  fprintf(fp,"# Summary of parameters and flags:\n");
  fprintf(fp,"# name       value\n");
  fprintf(fp,"# A          %f\n",w_amp);
  fprintf(fp,"# g          %f\n",w_pow);
  fprintf(fp,"# r          %f\n",reject_lev);
  if(!Ebeling) fprintf(fp,"# s          %f\n",conf_lev);
  fprintf(fp,"# C          %s %s %s %s\n",ttype[0],ttype[1],ttype[2],ttype[3]);
  fprintf(fp,"# B          %s %s %s %s\n",
	  ttypeH[0],ttypeH[1],ttypeH[2],ttypeH[3]);
  fprintf(fp,"# F          %f %f %f %f\n",frame[0],frame[1],frame[2],frame[3]);
  if(Holes) fprintf(fp,"# holesfile  %s\n",file_holes);
  fprintf(fp,"# datafile   %s\n",file_name);
  fprintf(fp,"# Verbose    %d\n",Verbose);
  fprintf(fp,"# Holes      %d\n",Holes);
  fprintf(fp,"# Ebeling    %d\n",Ebeling);
  fprintf(fp,"# Poisson    %d\n",Poisson);
  fprintf(fp,"# Fits       %d\n",Fits);
  fprintf(fp,"# XYcoord    %d\n",XYcoord);
  fprintf(fp,"# NotFind    %d\n",NotFind);
  fprintf(fp,"# RawFind    %d\n",RawFind);
  for(i=0;i<nclusters;i++,p_cir++) if(p_cir->edge==0){
    tmp=1.; tmp1=1.; tmp2=1.;
    if(! XYcoord){ 
      tmp=cos((p_cir->yc+ymean)*M_PI/180.);
      tmp1=cos(p_cir->pa)*cos(p_cir->pa)/tmp+sin(p_cir->pa)*sin(p_cir->pa);
      tmp2=sin(p_cir->pa)*sin(p_cir->pa)/tmp+cos(p_cir->pa)*cos(p_cir->pa);
      tmp1=sqrt(tmp1);
      tmp2=sqrt(tmp2);
    }
    llrow[0][0]=id=p_cir->id;     
    Rrow[0][0]=ra=xmean+p_cir->xc/tmp;
    Rrow[1][0]=dec=p_cir->yc+ymean;
    Rrow[2][0]=z=p_cir->z_median;
    Rrow[3][0]=z_mean=p_cir->z_mean;
    Rrow[4][0]=z_rms=p_cir->z_rms;
    Rrow[5][0]=contrast=p_cir->contrast;
    Rrow[6][0]=sig=p_cir->sig;
    Rrow[7][0]=size=p_cir->radius/tmp;
    Rrow[8][0]=a=p_cir->semi_maj*tmp1;
    Rrow[9][0]=b=p_cir->semi_min*tmp2;
    Rrow[10][0]=pa=p_cir->pa*(180./M_PI)/tmp1;
    Rrow[11][0]=areac=p_cir->area/tmp;
    Rrow[12][0]=min_dens=p_cir->minflux;     
    Rrow[13][0]=ra_c=xmean+p_cir->xgc/tmp;
    Rrow[14][0]=dec_c=ymean+p_cir->ygc;
    Rrow[15][0]=z_c=p_cir->zgc;
    Rrow[16][0]=x_0=p_cir->xc;
    Rrow[17][0]=y_0=p_cir->yc;
    Rrow[18][0]=w_amp;
    Rrow[19][0]=w_pow;
    irow[0][0]=nvt=p_cir->nvor;
    irow[1][0]=nbg=p_cir->nbg;
    irow[2][0]=ncorr=p_cir->n;
    fprintf(fp, "%lld %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf "
	    "%lf %lf %lf %lf %lf %lf %lf %d %d %d\n",
	    id, ra, dec, z, z_mean, z_rms, contrast, sig, size, a, b, pa,
	    areac, min_dens, ra_c, dec_c, z_c, x_0, y_0, w_amp,w_pow,
	    nvt, nbg, ncorr);
    fits_write_col_lnglng(fptr, colnum++,frow,felem,nelem,llrow[0],&status);
    for(j=0;j<20;j++)
      fits_write_col_dbl(fptr,colnum++,frow,felem,nelem,Rrow[j],&status);
    for(j=0;j<3;j++)
      fits_write_col_int(fptr,colnum++,frow,felem,nelem,irow[j],&status);
    cfitsio_error(status);
    colnum=1;
    frow++;
  }
  /*  
  char comment[70]="";
  sprintf(comment,"VT cluster catalog.");
  fits_write_comment(fptr,comment,&status);
  sprintf(comment,"Original file: %s",file_name);
  fits_write_comment(fptr,comment,&status);
  sprintf(comment,"Redshift interval: ( %5.2f , %5.2f )", 
	  redshmin, redshmax);
  fits_write_comment(fptr,comment,&status);
  sprintf(comment,"Center of frame: ( %12.8f , %12.8f , %12.8f )", 
	  xmean, ymean,zmean);
  fits_write_comment(fptr,comment,&status);
  sprintf(comment,"Summary of parameters and flags:");
  fits_write_comment(fptr,comment,&status);
  sprintf(comment,"name       value");
  fits_write_comment(fptr,comment,&status);
  sprintf(comment,"A          %f",w_amp);
  sprintf(comment,"g          %f\n",w_pow);
  fits_write_comment(fptr,comment,&status);
  sprintf(comment,"r          %f\n",reject_lev);
  fits_write_comment(fptr,comment,&status);
  if(!Ebeling) {
    sprintf(comment,"s          %f\n",conf_lev);
    fits_write_comment(fptr,comment,&status);
  }
  sprintf(comment,"C          %s %s %s %s\n",
	  ttype[0],ttype[1],ttype[2],ttype[3]);
  fits_write_comment(fptr,comment,&status);
  sprintf(comment,"B          %s %s %s %s\n",
	  ttypeH[0],ttypeH[1],ttypeH[2],ttypeH[3]);
  fits_write_comment(fptr,comment,&status);
  sprintf(comment,"F          %f %f %f %f\n",
	  frame[0],frame[1],frame[2],frame[3]);
  if(Holes) {
    sprintf(comment,"holesfile  %s\n",file_holes);
    fits_write_comment(fptr,comment,&status);
  }  
  sprintf(comment,"datafile   %s\n",file_name);
  fits_write_comment(fptr,comment,&status);
  sprintf(comment,"Verbose    %d\n",Verbose);
  fits_write_comment(fptr,comment,&status);
  sprintf(comment,"Holes      %d\n",Holes);
  fits_write_comment(fptr,comment,&status);
  sprintf(comment,"Ebeling    %d\n",Ebeling);
  fits_write_comment(fptr,comment,&status);
  sprintf(comment,"Poisson    %d\n",Poisson);
  fits_write_comment(fptr,comment,&status);
  sprintf(comment,"Fits       %d\n",Fits);
  fits_write_comment(fptr,comment,&status);
  sprintf(comment,"XYcoord    %d\n",XYcoord);
  sprintf(comment,"NotFind    %d\n",NotFind);
  fits_write_comment(fptr,comment,&status);
  sprintf(comment,"RawFind    %d\n",RawFind);
  fits_write_comment(fptr,comment,&status);
  */
  free2(Rrow);
  free2(llrow);
  free2(irow);
  fclose(fp);
  fits_close_file(fptr, &status);
  if(Verbose) printf("Cluster catalog saved.\n");
} /* WriteClusterCatalogFits() */

void AssignMembers(Connect *conn){
  register int i;
  Footprint *p_cir=footprint;
  Connect *p_conn=conn;
  for(i=0;i<ndata;i++) mask1[i]=0;  
  while (p_conn != NULL) { 
    for(i=0;i<p_conn->n;i++) mask1[p_conn->list[i]]=p_cir->id; 
    p_conn=p_conn->link; 
    p_cir++;
  } 
  if(Verbose) printf("Cluster members assigned.\n"); 
}

void SetFootprints(Connect *conn) {
  register int i, ii;
  int *boundary, *tmp_list, n_boundary, nlist;
  REAL *poly, *distances, *sort_dist, dist_max=0., tmp=0.,
    total_flux, xc, yc, xp=0., yp=0.,area_max, area_min,
    dxdx, dydy, dxdy, mxx, myy, mxy,t1,semi_maj,semi_min;
  Footprint *p_cir=footprint;
  Connect *p_conn=conn;

  tmp_list=boundary=(int *) NULL;
  poly=distances=sort_dist=(REAL *) NULL;
  ch_points=(coord *) NULL; ch_P=(coord **) NULL;

  while(p_conn != NULL) {
    /* set id */
    p_cir->id=++FieldID;
    /* set significance of detection */
    p_cir->sig=p_conn->sig;
    /* set min cell density */
    area_max=0.; 
    for(i=0; i<p_conn->n;i++)
      if((tmp=area[p_conn->list[i]]) > area_max) area_max=tmp;
    p_cir->minflux=1./area_max;
    /* set coords of central galaxy */
    area_min=area_max+1.; 
    for(i=0; i<p_conn->n;i++)
      if((tmp=area[p_conn->list[i]]) < area_min) {
	area_min=tmp; ii=p_conn->list[i];
      }
    p_cir->xgc=in.pointlist[2*ii]; 
    p_cir->ygc=in.pointlist[2*ii+1]; 
    p_cir->zgc=redsh[ii];
    central[ii]=1;
    /* Find the convex hull (ch2d code) */
    ch_points=(coord *) realloc(ch_points, p_conn->n*sizeof(coord));
    ch_P=(coord **) realloc(ch_P,(p_conn->n+1)*sizeof(coord*));
    for(i=0; i<p_conn->n;i++) {
      ii=2*p_conn->list[i];
      ch_points[i].x=in.pointlist[ii];
      ch_points[i].y=in.pointlist[ii+1];
      ch_P[i]=&ch_points[i];
    }
    n_boundary=ch2d(ch_P, p_conn->n);
    boundary=(int *) NULL; poly=(REAL *) NULL;
    boundary=(int *) realloc(boundary, n_boundary*sizeof(int));
    poly = (REAL *) realloc(poly, 2*n_boundary*sizeof(REAL));
    for(i=0; i<n_boundary; i++) {
      boundary[i]=p_conn->list[ch_P[i]-ch_points];
      poly[2*i]=ch_P[i]->x; poly[2*i+1]=ch_P[i]->y;
    }
    for(i=0;i<ndata;i++)
      if(mask1[i]==0 && 
	 inpoly(poly, n_boundary,in.pointlist[2*i], in.pointlist[2*i+1])){
        (p_conn->n)++;
        p_conn->list=(int *) realloc(p_conn->list, p_conn->n*sizeof(int));
        p_conn->list[p_conn->n-1]=i;
        mask1[i]=1;
      }
    /* set the baricenter, density contrast and total area of the footprint */
    p_cir->area = p_cir->xc=p_cir->yc=0.;
    total_flux=0.;
    for(i=0;i<p_conn->n;i++) {
      ii=p_conn->list[i]; tmp=area[ii];
      if(tmp > 0.) {
        p_cir->area += tmp;
        total_flux += 1./tmp;
        p_cir->xc += in.pointlist[2*ii]/tmp;
        p_cir->yc += in.pointlist[2*ii+1]/tmp;
      }
    }
    p_cir->xc /= total_flux;
    p_cir->yc /= total_flux;
    p_cir->contrast=total_flux/background_density;
    /* set numbers of galaxies: n bbg and nvor */
    p_cir->nvor = p_conn->n;
    p_cir->nbg=rint(p_cir->area*background_density);
    p_cir->n=rint(p_conn->n-p_cir->nbg);
    /* set edge flag to 0 (it will be 1 if any member is on buffer region) */
    p_cir->edge=0;
    /* Compute distances from and weigthed moments about the baricenter */
    xc=p_cir->xc; yc=p_cir->yc;
    distances=sort_dist=(REAL *) NULL;
    distances=(REAL *) realloc(distances,p_conn->n*sizeof(REAL));
    sort_dist=(REAL *) realloc(sort_dist,p_conn->n*sizeof(REAL));
    mxx=myy=mxy=0;
    for(i=0; i<p_conn->n; i++) {
      ii=p_conn->list[i];
      xp=in.pointlist[2*ii]; yp=in.pointlist[2*ii+1];
      dxdx=(xp-xc)*(xp-xc); dydy=(yp-yc)*(yp-yc); dxdy=(xp-xc)*(yp-yc);
      sort_dist[i]=distances[i]=dxdx+dydy;
      mxx+=dxdx/area[ii]; myy+=dydy/area[ii]; mxy+=dxdy/area[ii];
      if(edge1[ii]>0) p_cir->edge=1; /* adjust edge flags */
    }
    mxx/=total_flux; myy/=total_flux; mxy/=total_flux;
    /* set PA, semi-major, semi-minor axis from the moments */
    t1=sqrt((mxx-myy)*(mxx-myy)+4*mxy*mxy);
    semi_maj=sqrt(0.5*(mxx+myy+t1));
    semi_min=sqrt(0.5*abs((mxx+myy-t1)));
    t1=(semi_maj*semi_maj-myy)/mxy;
    p_cir->pa=atan(t1);
    p_cir->semi_maj=semi_maj;
    p_cir->semi_min=semi_min;
    /* set dist_max from center */
    qsort((REAL *) sort_dist, p_conn->n, sizeof(REAL), real_compare);
    nlist=(int) floor(.7*p_conn->n);
    dist_max=sort_dist[nlist-1];
    p_cir->radius = sqrt(dist_max);

    /* go to next cluster */
    p_conn=p_conn->link;
    p_cir++;

  } /* while p_conn != NULL */

  free(ch_points);
  free(ch_P);
  free(boundary);
  free(poly);
  free(distances);
  free(tmp_list);

  if(Verbose) printf("Cluster footprint properties computed.\n");

} /* SetFootprints() */

int ch2d(coord **ch_P, int n) {
  int u = make_chain(ch_P, n, cmpl);            /* make lower hull */
  if (!n) return 0;
  ch_P[n] = ch_P[0];
  return u+make_chain(ch_P+u, n-u+1, cmph);     /* make upper hull */
}

int inpoly(                       /*   1=inside, 0=outside                */
           REAL *poly,            /*   polygon points, [2i]=x, [2i+1]=y   */
           int npoints,           /*   number of points in polygon        */
           REAL xt, REAL yt) {    /*   coordinates of target point        */

/***************************************************************************
 *                                                                         *
 *   INPOLY.C                                                              *
 *                                                                         *
 *   Copyright (c) 1995-1996 Galacticomm, Inc.  Freeware source code.      *
 *                                                                         *
 ***************************************************************************/

/* is target point inside a 2D polygon? */

  REAL xnew,ynew,xold,yold,x1,y1,x2,y2;
  int i, inside=0;

  if (npoints < 3) return 0;
  xold=poly[2*(npoints-1)], yold=poly[2*(npoints-1)+1];
  for (i=0 ; i < npoints ; i++) {
    xnew=poly[2*i], ynew=poly[2*i+1];
    if (xnew > xold) {
      x1=xold; y1=yold;
      x2=xnew; y2=ynew;
    } else {
      x1=xnew; y1=ynew;
      x2=xold; y2=yold;
    }
    if ((xnew < xt) == (xt <= xold)         /* edge "open" at left end */
        && (yt-y1)*(x2-x1)< (y2-y1)*(xt-x1)) inside=!inside;
    xold=xnew, yold=ynew;
  }
  return inside;
} /* inpoly() */

int cmpl(const void *a, const void *b) {
  REAL v;
  CMPX(a,b);
  CMPY(b,a);
  return 0;
} 

int cmph(const void *a, const void *b) {
  return cmpl(b,a);
}


int make_chain(coord** V, int n, int (*cmp)(const void*, const void*)) {
  int i, j, s = 1;
  coord *t;

  qsort(V, n, sizeof(coord*), cmp);
  for (i=2; i<n; i++) {
    for (j=s; j>=1 && ccw(V, i, j, j-1); j--){}
    s = j+1;
    t = V[s]; V[s] = V[i]; V[i] = t;
  }
  return s;
}

int ccw(coord **ch_P, int i, int j, int k) {
  REAL a = ch_P[i]->x - ch_P[j]->x,
    b = ch_P[i]->y - ch_P[j]->y,
    c = ch_P[k]->x - ch_P[j]->x,
    d = ch_P[k]->y - ch_P[j]->y;
  return a*d - b*c <= 0;         /* true if points i, j, k counterclockwise */
}


