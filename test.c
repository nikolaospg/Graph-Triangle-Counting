//Source file containing functions that help us take an mtx file and take the associated CSC format of the array

#include "mmio.c"
#include <string.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>

/** 
 * Type definition of a struct, which resembles the CSC data structure 
 * (no values vector needed because all values are equal to 1)
 **/

typedef struct{            
    uint32_t* rowVector;     //array containing the row indeices of each nonzero
    uint32_t* colVector;     //array containing the index of the elements which start a column of the sparse matrix
    uint32_t nz;             //number of nonzero elements
    uint32_t M;              //number of columns (=number of rows because the matrix is square)
} CSCArray;

/**
 * Function that frees the memory allocated for the rowVector and colVector of a specific CSCArray structure.
 * Input:
 *      CSCArray* arg: pointer to the CSCArray structure we want to examine
 * Output:
 *      None
 **/
void CSCArrayfree(CSCArray* arg){
    free(arg->colVector);
    free(arg->rowVector);
}

/**
 * This function takes a FILE* associated with an mtx file (with the characteristics that are wanted in this project) and 
 * returns a CSC structure. The .mtx file contains only the coordinates of the nonzero elements of the lower triangular matrix
 * (because the matrix is symmetric). This function takes advantage of the symmetry of the matrix to create a CSC structure containing
 * all the info needed for all nonzero elements of the matrix.
 * Input:
 *      FILE* stream: pointer to the .mtx file cointaining the info of the sparse matrix
 * Output:
 *      CSCArray* retVal: pointer to a CSCArray structure representing the sparse matrix in CSC format
 **/

CSCArray* COOtoCSC(FILE* stream){        

    printf("Started converting mtx file to CSC\n");
    
    //Aquiring data about the sizes
    int32_t M,N,nz;
    
    mm_read_mtx_crd_size(stream, &M, &N, &nz);

    //Finding out how many digits the M number is comprised of (used to create the buffer)
    uint32_t Mdigits=0;
    uint32_t digitCount=M;
    while (digitCount !=0 ) {
        digitCount /= 10;     // n = n/10
        ++Mdigits;
    }

    char buffer[Mdigits+1];                              //buffer used to process the data from each line
    uint32_t* colVector=malloc((M+1)*sizeof(uint32_t));  //index of the elements which start a column of A for the lower triangular part of the matrix
    uint32_t* rowVector=malloc(nz*sizeof(uint32_t));     //row indices of each non zero element for the lower triangular part of the matrix
    colVector[0]=0; 

    if(colVector==NULL){
        printf("Error in COOtoCSC: Couldn't allocate memory for colVector");
        exit(-1);
    }

    if(rowVector==NULL){
        printf("Error in COOtoCSC: Couldn't allocate memory for rowVector");
        exit(-1);
    }        

    uint32_t rowIndiceRead;      //the row indice we read at each line of the .mtx file (first number on the line)
    uint32_t colIndiceRead;      //the column indice we read at each line of the .mtx file (second number on the line)
    uint32_t elemsUntilZeroCols; //temporarily store the number of elements up until the consecutive all-zero columns
    uint32_t colCheck=M+1;       //the integer used to understand whether we have finished converting a column. Initialized with M+1 so that it doesn't go inside the if loop the first time
    uint32_t lowerColElements=0; //how many nonzero elements we have in a specific column of the lower triangular matrix
    uint32_t colIndex=1;         //integer used as the index for the filling of the colVector 
    uint32_t upperColElements=0;  //how many nonzero elements we have in a specific column of the upper triangular matrix

    uint32_t **upperVectors = malloc(M*sizeof(uint32_t*));       //array containing the row indices of each nonzero element in each column for the upper triangular matrix
    
    if(upperVectors==NULL){
        printf("Error in COOtoCSC: Couldn't allocate memory for upperVectors");
        exit(-1);
    }

    //Creating the vectors for the upper triangular part of the matrix. The first element of each vector tells us the number of elements of the vector
    for(uint32_t i=0; i<M; i++){     
        upperVectors[i]=malloc(sizeof(uint32_t));
        if(upperVectors[i]==NULL){
            printf("Error in COOtoCSC: Couldn't allocate memory for upperVectors[%d]", i);
            exit(-1);
        }
        upperVectors[i][0]=1;  
    }
   
    //The loop that will fill out rowVector, colVector and upperVectors
    for(uint32_t i=0; i<nz; i++){ 

        fscanf(stream,"%s",buffer);
        rowIndiceRead=atoi(buffer)-1;
        fscanf(stream,"%s",buffer);
        colIndiceRead=atoi(buffer)-1;

        //Check for nonzero elements in the main diagonal
        if(rowIndiceRead == colIndiceRead){
            printf("There are elements in the main diagonal. Please give me an mtx without elements in the diagonal.\n");
            printf("The row of the first element in the diagonal is the %d\n",rowIndiceRead);
            exit(-1);
        }

        rowVector[i]=rowIndiceRead;
        
        if(colCheck<colIndiceRead){             //Check if the column indice we got is bigger than the previous column indice we examined
            if(colIndiceRead-colCheck>1){       //Checking if one or more consecutive columns have only elements equal to zero
                elemsUntilZeroCols = colVector[colIndex-1] + lowerColElements;
                for (uint32_t k=0; k<(colIndiceRead-colCheck); k++){
                    colVector[colIndex]=elemsUntilZeroCols;   //For all these all-zero columns put in the respective column array the value of the total elements up until that point
                    colIndex++;
                }
                lowerColElements=0;
            }
            //If the current element is in a different column than the previous one
            else{ 
                colVector[colIndex]=colVector[colIndex-1] + lowerColElements;
                colIndex++;
                lowerColElements=0;
            }
        }

        lowerColElements++;
        colCheck=colIndiceRead;

        //Note: the equivalent of a csc down triangular matrix is a crs upper triangular matrix. This is why use the row indices here as column indices and vice versa.
        upperVectors[rowIndiceRead][0]++;        //Increase the element counter of the vector of the specific column
        upperVectors[rowIndiceRead]=realloc(upperVectors[rowIndiceRead], (upperVectors[rowIndiceRead][0])*sizeof(uint32_t));
        if(upperVectors[rowIndiceRead]==NULL){
            printf("Error in COOtoCSC: Couldn't reallocate memory for upperVectors[%d]", rowIndiceRead);
            exit(-1);
        }
        upperColElements=upperVectors[rowIndiceRead][0];
        upperVectors[rowIndiceRead][upperColElements-1]=colIndiceRead; //Add in upperVectors the symmetric element of the one we just read from the file stream
    }

    //Last element of the colVector containing the number of elements of the down triangular matrix
    colVector[colIndex]=colVector[colIndex-1]+lowerColElements;  


    //Filling the last values, in case the last columns are all-zero columns 
    while(colIndex<M){           
        colIndex++;
        colVector[colIndex]=colVector[colIndex-1];
    }  

    uint32_t* finalRowVector = calloc(2*nz,sizeof(uint32_t));   //row indices of each non zero element for the whole matrix
    uint32_t rowVectorCount=0;                                  //shows how many row indices we have added in the finalRowVector
    uint32_t* finalColVector = calloc((M+1),sizeof(uint32_t));  //index of the elements which start a column of the whole matrix
    uint32_t colVectorCount=0;                                  //number of nonzero elements in a particular column for the whole sparse matrix
    uint32_t cscInitialColElems;                                //number of nonzero elements in the lower triangular part of the matrix

    if(finalRowVector==NULL){
        printf("Error in COOtoCSC: Couldn't allocate memory for finalRowVector");
        exit(-1);
    }

    if(finalColVector==NULL){
        printf("Error in COOtoCSC: Couldn't allocate memory for finalColVector");
        exit(-1);
    }

    //The loop that will fill out finalRowVector and finalColVector
    for(uint32_t i=0; i<M; i++){
        
        //Getting the values of the upper triangular half of the matrix
        //Check for nonzero elements in the column
        if(upperVectors[i][0]>1){
            //Add these elements on the final row vector
            for(uint32_t j=0; j<upperVectors[i][0]-1; j++){
                finalRowVector[rowVectorCount] = upperVectors[i][j+1];
                rowVectorCount ++;
                colVectorCount ++;
            }     
            free(upperVectors[i]);
        }
        else{
            free(upperVectors[i]);
        }

        //Getting the values of the lower triangular half of the matrix
        //Check for nonzero elements in the column
        if(i<M-1){
            cscInitialColElems=colVector[i+1] - colVector[i];
            //Add these elements on the final row vector
            for (uint32_t j=0; j<cscInitialColElems; j++){             
                finalRowVector[rowVectorCount] = rowVector[colVector[i] +j];
                rowVectorCount ++;
                colVectorCount ++;
            }
        }
        finalColVector[i+1] = finalColVector[i] + colVectorCount;
        colVectorCount=0;
    }  

    free(upperVectors);
    free(colVector);
    free(rowVector);
    
    //Creating the CSCArray to be returned
    CSCArray* retVal=malloc(sizeof(CSCArray));

    if(retVal==NULL){
        printf("Error in COOtoCSC: Couldn't allocate memory for retVal");
        exit(-1);
    }

    retVal->colVector=finalColVector;
    retVal->rowVector=finalRowVector;
    retVal->M=M;
    retVal->nz=rowVectorCount;    
    
    printf("The conversion is over\n");
    return retVal;  
}