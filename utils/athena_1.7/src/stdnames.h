/* ============================================================ *
 * stdnames.h                                                   *
 * Martin Kilbinger 2009                                        *
 * ============================================================ */

#ifndef __STDNAMES_H
#define __STDNAMES_H

/* Default file names */
#define CONFIG_DEF_NAME                        "config_tree"
#define SHEAR_CORR_NAME                        "xi"
#define SHEAR_CORR_REF_NAME                    ""
#define CORR_NAME                              "w"
#define GAL_SHEAR_XCORR_NAME                   "wgl"
#define WEIGHTED_ANGULAR_CORR_NAME             "ww"
#define LOG_NAME                               "log"
#define SHEAR_CORR_NAME2                       ""
#define SHEAR_RESAMPLE_CORR_NAME               SHEAR_CORR_NAME".resample"
#define GAL_SHEAR_RESAMPLE_XCORR_NAME          GAL_SHEAR_XCORR_NAME".resample"
#define SHEAR_ALL_RESAMPLE_CORR_NAME_XIP       ""
#define SHEAR_ALL_RESAMPLE_CORR_NAME_XIM       ""
#define SHEAR_RESAMPLE_COV_NAME                SHEAR_RESAMPLE_CORR_NAME".cov"
#define GAL_SHEAR_COV_XCORR_NAME               GAL_SHEAR_RESAMPLE_XCORR_NAME".cov"
#define ID_GAL_SHEAR_ALL_RESAMPLE_CORR_NAME_GT ""
#define SUFFIX_ASCII                           ""
#define SUFFIX_FITS                            ".fits"

/* Default input column names */
#define COL_NAME_DEFAULT_X                        "ra"
#define COL_NAME_DEFAULT_Y                        "dec"
#define COL_NAME_DEFAULT_Z                        "z"
#define COL_NAME_DEFAULT_E1                       "e1"
#define COL_NAME_DEFAULT_E2                       "e2"
#define COL_NAME_DEFAULT_W                        "w"
#define COL_NAME_DEFAULT_NJK                      "njk"

/* Separator char for pair (column type, column name) */
#define COL_NAME_SEP                              ":"

/* Number of output columns */
#define NCOL_XI                                 8
#define NCOL_XI_RES                             5
#define NCOL_XI_COV_RES                         5
#define NCOL_WGL                                7
#define NCOL_WGL_RES                            3
#define NCOL_COV_COL                            3
#define NCOL_W					2

/* Default output column names */

/* xi */
#define COL_NAME_DISTANCE_2D		"theta"
#define COL_NAME_XI_P			"xi_p"
#define COL_NAME_XI_M			"xi_m"
#define COL_NAME_XI_X			"xi_x"
#define COL_NAME_W_TOT			COL_NAME_DEFAULT_W
#define COL_NAME_SQRT_D			"sqrt_D"		
#define COL_NAME_SQRT_D_COR		"sqrt_Dcor"
#define COL_NAME_N_PAIR			"n_pair"
#define COL_NAME_N_PAIR_RESAMPLE        COL_NAME_N_PAIR"_resample"

/* xi resampled mean and rms */
#define SUFFIX_RESAMPLE			"_resample"
#define COL_NAME_XI_P_RESAMPLE		COL_NAME_XI_P""SUFFIX_RESAMPLE
#define COL_NAME_XI_M_RESAMPLE		COL_NAME_XI_M""SUFFIX_RESAMPLE
#define COL_NAME_RMS_P_RESAMPLE		"rms_p"SUFFIX_RESAMPLE
#define COL_NAME_RMS_M_RESAMPLE		"rms_m"SUFFIX_RESAMPLE

/* xi resampled covariance */
#define COL_NAME_DISTANCE_2D_1          COL_NAME_DISTANCE_2D"_1"
#define COL_NAME_DISTANCE_2D_2          COL_NAME_DISTANCE_2D"_2"
#define COL_NAME_COV_PP                 "Cov_pp"
#define COL_NAME_COV_MM                 "Cov_mm"
#define COL_NAME_COV_PM                 "Cov_pm"

/* wgl */
#define COL_NAME_GT                     "g_t"
#define COL_NAME_GX                     "g_x"

/* wgl resampled mean and rms */
#define COL_NAME_GT_RESAMPLE           COL_NAME_GT""SUFFIX_RESAMPLE
#define COL_NAME_GT_RMS_RESAMPLE       COL_NAME_GT"_rms"SUFFIX_RESAMPLE

/* wgl resampled covariance */
#define COL_NAME_GT_COV                "Cov_"COL_NAME_GT


/* Output file indices */
#define NNAMES 13
#define ID_SHEAR_CORR            	  0
#define ID_SHEAR_CORR_REF        	  1
#define ID_CORR                  	  2
#define ID_GAL_SHEAR_XCORR       	  3
#define ID_WEIGHTED_ANGULAR_CORR 	  4
#define ID_LOG                   	  5
#define ID_SHEAR_CORR2           	  6
#define ID_SHEAR_RESAMPLE_CORR   	  7
#define ID_SHEAR_ALL_RESAMPLE_CORR_XIP    8
#define ID_SHEAR_ALL_RESAMPLE_CORR_XIM    9
#define ID_GAL_SHEAR_RESAMPLE_XCORR       10
#define ID_GAL_SHEAR_ALL_RESAMPLE_CORR_GT 11
#define ID_SUFF                           12

#endif
