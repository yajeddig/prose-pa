/*-------------------------------------------------------------------------------
* 
* SOFTWARE NAME: ProSe-PA
* FILE NAME: matrix_operation.c
* BRANCH NAME: main
* 
* CONTRIBUTORS: Shuaitao WANG, Lauriane VILMIN, Aurélien BORDET, Masihullah HASANYAR, Thomas ROMARY, Nicolas FLIPO
*
* PROJECT MANAGER: Nicolas FLIPO
* 
* SOFTWARE BRIEF DESCRIPTION: The ProSe-PA is a software for simulating the hydro-biogeochemical functioning of rivers, particularly heavily urbanised rivers, and streams. The sofware can
* operate in two modes: direct calculation or data assimilation.
*
* In direct calculation mode, based on a semi-implicit Eulerian numerical scheme, the software simulates the functioning of the water column in contact with a benthic compartment made up of unconsolidated
* sediments and periphyton (librive library). It can be used to simulate the anthropisation of environments, through the explicit representation of developments such as navigation dams, sluice
* gates and river navigation, as well as discharges into the environment, such as those from wastewater treatment plants or combined sewer overflows.
* The software explicitly simulates the growth of micro-organisms in the water column and in the benthic compartment, enabling the carbon, oxygen and nutrient (nitrogen, phosphrus, silica) cycles
* associated with these biological processes to be quantified (librive library). Water temperature is also simulated by the software (libseb library), as are particulate and dissolved exchanges
* between the water column and the benthic compartment. The software can simulate 1D, pseudo-2D hydraulics of river and streams (discharge, water height) using the libhyd library. The advection-dispersion 
* process is simulated using libttc library.
* 
* In data assimilation mode, ProSe-PA includes two filters for assimilating high frequency dissolved oxygen data. These two filters are a particle filter and the ensemble Kalman filter.
*
* ANSI C software developed at the Geosciences and geoengineering Department, joint research center of Mines Paris-PSL and ARMINES, Fontainebleau, France. The code is based on the coupling 
* of 12 libraries developed also at the Geosciences and geoengineering Department, mostly in ANSI C: libprint, libts, libpc, libchronos, libio, libhyd, libtube, libttc, librive, libseb, libmb, scripts.
*
* CITATION: 
* Wang, S., Flipo, N., Romary, T.. (2019). Oxygen data assimilation for estimating micro-organism communities' parameters in river systems. Water Research, 165, 115021. doi:10.1016/j.watres.2019.115021
* Flipo, N., Even, S., Poulin, M., Tusseau-Vuillemin, M-H., Ameziane, T., Dauta, A. (2004). Biogeochemical modelling at the river scale: plankton and periphyton dynamics (Grand Morin case study, France).
*    Ecol. Model. 176(3-4), 333-347. doi:10.1016/j.ecolmodel.2004.01.012. 
* Even S, Poulin M, Garnier J, Billen G, Servais P, Chesterikoff A (1998). River ecosystem modelling: application of the PROSE model to the Seine River (France). Hydrobiologia, 373, pp. 27-45.
*    doi: 10.1023/A:1017045522336
* Vilmin, L, Aissa, N., Garnier, J., Billen, G., Mouchel, J-M., Poulin, M., Flipo, N. (2015). Impact of hydro-sedimentary processes on the dynamics of soluble reactive phosphorus in the Seine River.
*    Biogeochemistry, 122, 229-251. doi:10.1007/s10533-014-0038-3
* Wang, S., Flipo, N., Romary, T. (2023). Which filter for data assimilation in water quality models? Focus on oxygen reaeration and heterotrophic bacteria activity. Journal of Hydrology, 620, 129423. 
*    doi:10.1016/j.jhydrol.2023.129423
* 
* COPYRIGHT: (c) 2023 Contributors to the ProSe-PA software. 
* CONTACT: Nicolas FLIPO <nicolas.flipo@minesparis.psl.eu>
*          
* 
* All rights reserved. This software and the accompanying materials
* are made available under the terms of the Eclipse Public License v2.0
* which accompanies this distribution, and is available at
* http://www.eclipse.org/legal/epl-v20.html
* 
*------------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <malloc.h>
#include <libprint.h>
#include "LA.h"

/*#define MATRIX_TYPE double
enum matrix_types {MATRIX_A, MATRIX_L, MATRIX_U, MATRIX_A_INV, MATRIX_NTYPE};

typedef struct LA_Matrix s_matrix_la;
#define new_matrix_la() ((s_matrix_la *) malloc(sizeof(s_matrix_la)))

MATRIX_TYPE *LA_calc_mean_by_row_matrix(s_matrix_la *, FILE *);
s_matrix_la *LA_subtract_mean_by_row_matrix(s_matrix_la *, FILE *);
s_matrix_la *LA_transpose_matrix(s_matrix_la *, FILE *);
s_matrix_la *LA_calc_mutiply_matrix(s_matrix_la *, s_matrix_la *, FILE *);
s_matrix_la *LA_addition_matrix(s_matrix_la *, s_matrix_la *, FILE *);
void LA_LU_decomposition_L(s_matrix_la *, FILE *);
void LA_print_matrix_type(s_matrix_la *, int , FILE *);
void LA_inverse_matrix(s_matrix_la *, FILE *);
void LA_calc_mutiply_constant(s_matrix_la *, double, FILE *);
void LA_free_matrix(s_matrix_la *, FILE *);
void LA_free_matrix_all(s_matrix_la *, FILE *);

struct LA_Matrix {

    int row_size;
    int col_size;
    
    // matrix data
    MATRIX_TYPE **data;
    
    //inverse matrix data
    //MATRIX_TYPE **data_inv;
    s_matrix_la *mat_inv;

    //LU decomposition lower triangular matrix L
    //MATRIX_TYPE **dataL;
    s_matrix_la *mat_L;

    //LU decomposition upper triangular matrix U, with U[i][i] = 1
    //MATRIX_TYPE **dataU;
    s_matrix_la *mat_U;
}; */


//function to calculate mean values of a matrix by row E(X)
MATRIX_TYPE *LA_calc_mean_by_row_matrix(s_matrix_la *matrix, FILE *fp)
{
    int nr, nc, row_size, col_size;
    MATRIX_TYPE *matrix_mean_by_row;
    MATRIX_TYPE mean_row;
    
    row_size = matrix->row_size;
    col_size = matrix->col_size;
    matrix_mean_by_row = (MATRIX_TYPE *) malloc(row_size * sizeof(MATRIX_TYPE));
    
    for(nr = 0; nr < row_size; nr++)
    {
        mean_row = 0.;
        for(nc = 0; nc < col_size; nc++)
            mean_row += matrix->data[nr][nc]/col_size;

        matrix_mean_by_row[nr] = mean_row;
    }
    return matrix_mean_by_row;
}

// X - E(X)
s_matrix_la *LA_subtract_mean_by_row_matrix(s_matrix_la *matrix, FILE *fp)
{
    int nr, nc;
    MATRIX_TYPE *matrix_mean_by_row;
    s_matrix_la *result;
    
    result = new_matrix_la();
    result->row_size = matrix->row_size;
    result->col_size = matrix->col_size;

    result->data = (MATRIX_TYPE **) malloc(result->row_size * sizeof(MATRIX_TYPE *));
    
    matrix_mean_by_row = LA_calc_mean_by_row_matrix(matrix, fp);

    for(nr = 0; nr < matrix->row_size; nr++)
    {
        result->data[nr] = (MATRIX_TYPE *) malloc(result->col_size * sizeof(MATRIX_TYPE));
        for(nc = 0; nc < matrix->col_size; nc++)
            result->data[nr][nc] = matrix->data[nr][nc] - matrix_mean_by_row[nr];

    }

    return result;
        
}

// X + Y
s_matrix_la *LA_addition_matrix(s_matrix_la *matrixL, s_matrix_la *matrixR, FILE *fp)
{
    int nr, nc;
    s_matrix_la *result;

    result = new_matrix_la();

    if((matrixL->row_size == matrixR->row_size) && (matrixL->col_size == matrixL->col_size))
    {
        result->row_size = matrixL->row_size;
        result->col_size = matrixL->col_size;

        result->data = (MATRIX_TYPE **) malloc(result->row_size * sizeof(MATRIX_TYPE *));
        for(nr = 0; nr < result->row_size; nr++)
        {
            result->data[nr] = (MATRIX_TYPE *) malloc(result->col_size * sizeof(MATRIX_TYPE));
            for(nc = 0; nc < result->col_size; nc++)
                result->data[nr][nc] = matrixL->data[nr][nc] + matrixR->data[nr][nc];

        }
    }
    else
        LP_error(fp, "Matrix dimensions are not coherent, row_sizeL = %d col_sizeL = %d row_siezR = %d col_sizeR = %d in file %s function %s line %d\n", matrixL->row_size, matrixL->col_size, matrixR->row_size, matrixR->col_size, __FILE__, __FUNCTION__, __LINE__);

    return result;
        
}

s_matrix_la *LA_addition_matrix_replace(s_matrix_la *matrixL, s_matrix_la *matrixR, FILE *fp)
{
    int nr, nc;
    s_matrix_la *result;

    if((matrixL->row_size == matrixR->row_size) && (matrixL->col_size == matrixL->col_size))
    {
        
        for(nr = 0; nr < matrixL->row_size; nr++)
        {
            
            for(nc = 0; nc < matrixL->col_size; nc++)
            {
                matrixL->data[nr][nc] = matrixL->data[nr][nc] + matrixR->data[nr][nc];
                //LP_printf(fp,"nr = %d nc = %d  matrixL->data[nr][nc] = %f\n",nr,nc,matrixL->data[nr][nc]);
            }

        }
    }
    else
        LP_error(fp, "Matrix dimensions are not coherent, row_sizeL = %d col_sizeL = %d row_siezR = %d col_sizeR = %d in file %s function %s line %d\n", matrixL->row_size, matrixL->col_size, matrixR->row_size, matrixR->col_size, __FILE__, __FUNCTION__, __LINE__);

    return matrixL;
        
}

// X^T
s_matrix_la *LA_transpose_matrix(s_matrix_la *matrix, FILE *fp)
{
    int nr, nc;
    s_matrix_la *result;
    
    result = new_matrix_la();
    result->row_size = matrix->col_size;
    result->col_size = matrix->row_size;

    result->data = (MATRIX_TYPE **) malloc(result->row_size * sizeof(MATRIX_TYPE *));
    for(nr = 0; nr < result->row_size; nr++)
    {
        result->data[nr] = (MATRIX_TYPE *) malloc(result->col_size * sizeof(MATRIX_TYPE));
        for(nc = 0; nc < result->col_size; nc++)
            result->data[nr][nc] = matrix->data[nc][nr];
    }

    return result;
}

// (X-E(X))(Y - E(Y))^T
s_matrix_la *LA_calc_mutiply_matrix(s_matrix_la *matrixL, s_matrix_la *matrixR, FILE *fp)
{
   int nrl, ncl, ncr;
   s_matrix_la *result;
   
   result = new_matrix_la();
   
   if(matrixL->col_size != matrixR->row_size)
       LP_error(fp, "Matrix dimensions are not coherent, col_sizeL = %d row_siezR = %d in file %s function %s line %d\n", matrixL->col_size, matrixR->col_size, __FILE__, __FUNCTION__, __LINE__);
   else
   {
       result->row_size = matrixL->row_size;
       result->col_size = matrixR->col_size;

       result->data = (MATRIX_TYPE **) malloc(result->row_size * sizeof(MATRIX_TYPE *));

       for(nrl = 0; nrl < matrixL->row_size; nrl++)
       {
           result->data[nrl] = (MATRIX_TYPE *) malloc(result->col_size * sizeof(MATRIX_TYPE));

           for(ncr = 0; ncr < matrixR->col_size; ncr++)
           {
               result->data[nrl][ncr] = 0.;
               for(ncl = 0; ncl < matrixL->col_size; ncl++)
                   result->data[nrl][ncr] += matrixL->data[nrl][ncl] * matrixR->data[ncl][ncr];
           }
       }
    return result;
    }
}

void LA_calc_mutiply_constant(s_matrix_la *matrix, double factor, FILE *fp)
{
    int nr, nc;

    for(nr = 0; nr < matrix->row_size; nr++)
    {
        for(nc = 0; nc < matrix->col_size; nc++)
            matrix->data[nr][nc] *= factor;
    }

}

s_matrix_la *LA_calc_mutiply_constant_noreplace(s_matrix_la *matrix, double factor, FILE *fp)
{
    int nr, nc;
    s_matrix_la *result;

    result = LA_matrix_calloc(matrix->row_size,matrix->col_size);

    for(nr = 0; nr < matrix->row_size; nr++)
    {
        for(nc = 0; nc < matrix->col_size; nc++)
            result->data[nr][nc] = matrix->data[nr][nc] * factor;
    }

    return result;
}

// X = LU, computation of L matrix
void LA_LU_decomposition_L(s_matrix_la *matrix, FILE *fp)
{
    s_matrix_la *L;
    double detL = 1.;
    int nr, nc, k;

    matrix->mat_L = new_matrix_la();
    matrix->mat_U = new_matrix_la();

    matrix->mat_L->row_size = matrix->row_size;
    matrix->mat_U->row_size = matrix->row_size;
    matrix->mat_L->col_size = matrix->col_size;
    matrix->mat_U->col_size = matrix->col_size;

    matrix->mat_L->data = (MATRIX_TYPE **) malloc(matrix->mat_L->row_size * sizeof(MATRIX_TYPE *));
    matrix->mat_U->data = (MATRIX_TYPE **) malloc(matrix->mat_U->row_size * sizeof(MATRIX_TYPE *));
    for(nr = 0; nr < matrix->row_size; nr++)
    {
        matrix->mat_L->data[nr] = (MATRIX_TYPE *) malloc(matrix->mat_L->col_size * sizeof(MATRIX_TYPE));
        memset(matrix->mat_L->data[nr], 0, matrix->mat_L->col_size * sizeof(MATRIX_TYPE));

        matrix->mat_L->data[nr][0] = matrix->data[nr][0];

        matrix->mat_U->data[nr] = (MATRIX_TYPE *) malloc(matrix->mat_U->col_size * sizeof(MATRIX_TYPE));
        
        memset(matrix->mat_U->data[nr], 0, matrix->mat_U->col_size * sizeof(MATRIX_TYPE));
        matrix->mat_U->data[nr][nr] = 1;
    }

    for(nc = 1; nc < matrix->col_size; nc++)
        matrix->mat_U->data[0][nc] = matrix->data[0][nc] / matrix->mat_L->data[0][0];
    for(nr = 1; nr < matrix->row_size; nr++)
    {
        for(nc = 1; nc < matrix->col_size; nc++)
        {
            if(nr >= nc)
            {
                matrix->mat_L->data[nr][nc] = matrix->data[nr][nc];
                for(k = 0; k <= nc - 1; k++)
                    matrix->mat_L->data[nr][nc] -= matrix->mat_L->data[nr][k] * matrix->mat_U->data[k][nc];
            }
            else
            {
                matrix->mat_U->data[nr][nc] = matrix->data[nr][nc] / matrix->mat_L->data[nr][nr];
                for(k=0; k <= nr - 1; k++)
                    matrix->mat_U->data[nr][nc] -= (matrix->mat_L->data[nr][k] * matrix->mat_U->data[k][nc]) / matrix->mat_L->data[nr][nr];
                
            }
            
        }
    }
    
    for(nr = 0; nr < matrix->mat_L->row_size; nr++)
        detL *= matrix->mat_L->data[nr][nr];
    if(detL == 0)
        LP_error(fp, "In FILE %s FUNCTION %s, LINE %d LU decomposition, the determiant of L = %f \n", __FILE__, __FUNCTION__, __LINE__,detL);

}

// inverse A : A*inv(A) = I with A = LU LU*inv(A) = I
// for each column, LUx = b with x is the column of inv(A) and b is the corresponding column I
// so we note Ux = z, Lz=b; we can solve firstly Lz = b; then Ux = z

void LA_inverse_matrix(s_matrix_la *matrix, FILE *fp)
{
   int nr, nc, k, j,i;
   MATRIX_TYPE *x, *b, *z;
  
   if((matrix->mat_L == NULL) || ((matrix->mat_U == NULL)))
       LP_error(fp, "LU decomposition needed, before solving inverse of matrix\n");
   
   matrix->mat_inv = new_matrix_la();

   matrix->mat_inv->row_size = matrix->row_size;
   matrix->mat_inv->col_size = matrix->col_size;

   matrix->mat_inv->data = (MATRIX_TYPE **) malloc(matrix->mat_inv->row_size * sizeof(MATRIX_TYPE *));
   for(nr = 0; nr < matrix->mat_inv->row_size; nr++)
   {
       matrix->mat_inv->data[nr] = (MATRIX_TYPE *) malloc(matrix->mat_inv->col_size * sizeof(MATRIX_TYPE));
       memset(matrix->mat_inv->data[nr], 0, matrix->mat_inv->col_size * sizeof(MATRIX_TYPE));
   }

   // solving Lz = b and Ux = z
   b = (MATRIX_TYPE *) malloc(matrix->mat_inv->row_size * sizeof(MATRIX_TYPE));
   x = (MATRIX_TYPE *) malloc(matrix->mat_inv->row_size * sizeof(MATRIX_TYPE));
   z = (MATRIX_TYPE *) malloc(matrix->mat_inv->row_size * sizeof(MATRIX_TYPE));
   //memset(b, 0, matrix->mat_inv->col_size * sizeof(MATRIX_TYPE));

   for(nc = 0; nc < matrix->mat_inv->col_size; nc++)
   {
       memset(b, 0, matrix->mat_inv->col_size * sizeof(MATRIX_TYPE));
       b[nc] = 1.;

       /*printf("nc = %d:\n",nc);
       printf("b:\n");
       for(i = 0; i < matrix->mat_inv->col_size; i++)
          printf("%f ", b[i]);
       printf("\n");*/

       // solve Lz = b by foward substitution
       z[0] = b[0] / matrix->mat_L->data[0][0];
       for(nr = 1; nr < matrix->mat_inv->row_size; nr++)
       {
          z[nr] = b[nr];
          for(k = 0; k < nr; k++)
             z[nr] -= matrix->mat_L->data[nr][k] * z[k];
          z[nr] /= matrix->mat_L->data[nr][nr];
       }

       /*printf("z:");
       for(i = 0; i < matrix->mat_inv->col_size; i++)
          printf("%f ", z[i]);
       printf("\n");*/


       // solve Ux = z by backward substitution
       x[matrix->mat_inv->row_size -1] = z[matrix->mat_inv->row_size -1];
       for(nr = matrix->mat_inv->row_size - 2; nr >=0; nr--)
       {
           x[nr] = z[nr];
           for(k = nr + 1; k < matrix->mat_inv->col_size; k++)
               x[nr] -= matrix->mat_U->data[nr][k] * x[k];
       }
       
       // x is the nc column of inverse matrix
       for(j = 0; j < matrix->mat_inv->row_size; j++)
           matrix->mat_inv->data[j][nc] = x[j];
   }
   
   free(x);
   free(z);
   free(b);
}

void LA_print_matrix_type(s_matrix_la *matrix, int type, FILE *fp)
{
    int nr, nc; 
    MATRIX_TYPE **mat=NULL;
    switch(type)
    {
        case MATRIX_A :
            mat = matrix->data;
            break;
        case MATRIX_L :
            mat = matrix->mat_L->data;
            break;
        case MATRIX_U :
            mat = matrix->mat_U->data;
            break;
        case MATRIX_A_INV :
            mat = matrix->mat_inv->data;
            break;
        default :
            printf("No matrix type defined\n");

    }
  
    if(mat != NULL)
    {
        for(nr = 0; nr < matrix->row_size; nr++)
        {
            for(nc = 0; nc < matrix->col_size; nc++)
                LP_printf(fp,"%.17f ", mat[nr][nc]);
            LP_printf(fp,"\n");
        }
    }
}

void LA_free_matrix(s_matrix_la *matrix, FILE *fp)
{
    int nr, nc;
   
    for(nr = 0; nr < matrix->row_size; nr++)
    {
        free(matrix->data[nr]);
    }
    free(matrix->data);
    free(matrix);
    matrix = NULL;
}

void LA_free_matrix_all(s_matrix_la *matrix, FILE *fp)
{
    
    if(matrix != NULL)
    {
        if(matrix->mat_L != NULL)
            LA_free_matrix(matrix->mat_L, fp);
        if(matrix->mat_U != NULL)
            LA_free_matrix(matrix->mat_U, fp);
        if(matrix->mat_inv != NULL)
            LA_free_matrix(matrix->mat_inv, fp);
       
        LA_free_matrix(matrix, fp);
    }
    
}


/* void main()
{
   s_matrix_la *mat_A;
   s_matrix_la *LU, *AINVA;
   int nr, nc;
   FILE *fp;
   mat_A = new_matrix_la();
   mat_A->row_size = 4;
   mat_A->col_size = 4;
  
   mat_A->data = (MATRIX_TYPE **) malloc(mat_A->row_size * sizeof(MATRIX_TYPE *));

   for(nr = 0; nr < mat_A->row_size; nr++)
       mat_A->data[nr] = (MATRIX_TYPE *) malloc(mat_A->col_size * sizeof(MATRIX_TYPE));

  if ((fp = fopen("result_mat.out","w")) == 0) 
    LP_error(stderr,"Impossible to open the debug file result_mat.out\n");

   mat_A->data[0][0] = -4;
   mat_A->data[0][1] = 3;
   mat_A->data[0][2] = 4;
   mat_A->data[0][3] = -2;
   mat_A->data[1][0] = 8;
   mat_A->data[1][1] = -9;
   mat_A->data[1][2] = -12;
   mat_A->data[1][3] = 6;
   mat_A->data[2][0] = -20;
   mat_A->data[2][1] = 12;
   mat_A->data[2][2] = 20;
   mat_A->data[2][3] = -12;
   mat_A->data[3][0] = -12;
   mat_A->data[3][1] = -3;
   mat_A->data[3][2] = -4;
   mat_A->data[3][3] = 0;
 
  
  LA_LU_decomposition_L(mat_A, fp);

  LP_printf(fp,"matrix A:\n");
  LA_print_matrix_type(mat_A,  MATRIX_A, fp);
  LP_printf(fp,"matrix L:\n");
  LA_print_matrix_type(mat_A,  MATRIX_L, fp);
  LP_printf(fp,"matrix U:\n");
  LA_print_matrix_type(mat_A,  MATRIX_U, fp);
  
  LU = LA_calc_mutiply_matrix(mat_A->mat_L, mat_A->mat_U, fp);

  LP_printf(fp,"matrix LU:\n");
  LA_print_matrix_type(LU,  MATRIX_A, fp);
  
  LA_inverse_matrix(mat_A, fp);

  LP_printf(fp,"Inverse of matrix A:\n");
  LA_print_matrix_type(mat_A,  MATRIX_A_INV, fp);
 
  AINVA = LA_calc_mutiply_matrix(mat_A, mat_A->mat_inv, fp);
  LP_printf(fp,"A * inv(A):\n");
  LA_print_matrix_type(AINVA,  MATRIX_A, fp);

  return;
  
}*/
