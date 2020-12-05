/**
 * This file includes the implementation of a tester to check that each of the rest of the implementations calculates
 * the correct number of triangles each node is adjacent to.
**/

#include <stdio.h>
#include "test.c"
#include <string.h>
#include <math.h>
#include <stdint.h>

//Structure containing the matrix we are examining in coo format
typedef struct{            
    uint32_t* coo_row;     //array containing the row indices of each nonzero element
    uint32_t* coo_col;     //array containing the column indices of each nonzero element
    uint32_t nz;           //number of nonzero elements
    uint32_t n;            //number of columns (=number of rows because the matrix is square)
} COOArray;


/**
 * Function that checks whether a specific column of the matrix (colNum) has an element in a specific row (wantedRow)
 * This is useful in deciding whether A[wantedRow][colNum] != 0 and calculating the number of triangles adjacen to each node
 * The implementation is based on binary search. Using binary search on rowVector we try to find the element with value wantetRow
 * If we find it we return its position on rowVector, otherwise we return -1.
 * Inputs:
 *      uint32_t* rowVector: the row indices array of the csc format
 *      uint32_t* colVector: the column changes array of the csc format
 *      uint32_t colNum: the index of the column we want to examine
 *      uint32_t wantedRow: the index of the row that we want our element to belong to
 * Outputs:
 *      int32_t result: the position on rowVector if there is an element in (wanted row, colNum), -1 otherwise
 **/

int32_t elementInColumnCheck2(uint32_t* rowVector,uint32_t* colVector, uint32_t colNum, uint32_t wantedRow){
    int32_t result=-1;

    int32_t left = colVector[colNum];       //first element of the sub-array in which we do our binary search
    int32_t right = colVector[colNum+1]-1;  //last element of the sub-array in which we do our binary search
    int32_t middle = (left+right)/2;

    //binary search
    while(left<=right){
        if(rowVector[middle]<wantedRow){
            left = middle+1;
        }
        else if(rowVector[middle]==wantedRow){
            result = middle;
            break;
        }
        else{
            right = middle-1;
        }
        middle = (left+right)/2;    
    }
    return result;
}

/**
 * Function that reads a .mtx file and extracts the nonzero elements written in it in coo format. Note that because 
 * we deal with symmetric matrices the .mtx file (and as a result the coo format of the matrix) contains only
 * the nonzero elements of the lower triangular part of the matrix. We create the whole matrix out of these elements in another function.
 * Input:
 *      FILe* stream: pointer to the .mtx file cointaining the info of the sparse matrix
 * Output:
 *      COOArray* cooArray: pointer to a COOArray structure representing the lower triangular matrix in coo format
**/

COOArray* createCOO(FILE* stream){
    printf("started converting mtx file to COO\n");

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

    char buffer[Mdigits+1];                             //buffer used to process the data from each line
    uint32_t* colVector=malloc(nz*sizeof(uint32_t));    //column indices of each non zero element for the lower triangular part of the matrix
    uint32_t* rowVector=malloc(nz*sizeof(uint32_t));    //row indices of each non zero element for the lower triangular part of the matrix

    if(colVector==NULL){
        printf("Error in COOtoCSC: Couldn't allocate memory for colVector");
        exit(-1);
    }

    if(rowVector==NULL){
        printf("Error in COOtoCSC: Couldn't allocate memory for rowVector");
        exit(-1);
    }

    //Read the row and column indices of each nonzero and store them in rowVector and colVector respectively
    for(int32_t i=0; i<nz; i++){ 
        fscanf(stream,"%s",buffer);
        rowVector[i]=atoi(buffer);
        fscanf(stream,"%s",buffer);
        colVector[i]=atoi(buffer);
    }

    //Creating the COOArray to be returned
    COOArray *cooArray = malloc(sizeof(COOArray));

    if(cooArray==NULL){
        printf("Error in COOtoCSC: Couldn't allocate memory for cooArray");
        exit(-1);
    }

    cooArray->coo_col = colVector;
    cooArray->coo_row = rowVector;
    cooArray->n = M;
    cooArray->nz = nz;

    printf("finished converting mtx file to COO\n");

    return cooArray;
}

/**
 * Function that gets the coo format of the lower triangular part of a sparse matrix A and transforms it to CSC format
 * for the whole matrix. In order to achive that the function takes advantage of the symmetry of the matrix.
 * Because of this symmetry the crs format of the upper triangular part of the matrix is the same as the csc format of 
 * the lower triangular part of the matrix. Knowing that, we can create the colVector and rowVector for the whole matrix.
 * Input:
 *      uint32_t *coo_row: array containing the row indices of each nonzero element for the lower triangular part of the matrix
 *      uint32_t *coo_col: array containing the column indices of each nonzero element for the lower triangular part of the matrix
 *      uint32_t n: number of columns(rows) of the sparse matrix
 *      uint32_t nz: number of nonzero elements of the lower triangular part of the matrix
 * Output:
 *      uint32_t* vector: array containing the number of triangles adjacent to each node
**/

uint32_t* vertexWiseTriangleCounts(uint32_t *coo_row, uint32_t *coo_col, uint32_t n, uint32_t nz){

    uint32_t* vector;               //array containing the number of triangles adjacent to each node

    uint32_t rowIndiceRead;         //the row indice we read at each line of the .mtx file (first number on the line)
    uint32_t colIndiceRead;         //the column indice we read at each line of the .mtx file (second number on the line)
    uint32_t elemsUntilZeroCols;    //temporarily store the number of elements up until the consecutive all-zero columns
    uint32_t colCheck=n+1;          //the integer used to understand whether we have finished converting a column. Initialized with M+1 so that it doesn't go inside the if loop the first time
    uint32_t lowerColElements=0;    //how many nonzero elements we have in a specific column of the lower triangular matrix
    uint32_t colIndex=1;            //integer used as the index for the filling of the colVector 
    uint32_t upperColElements=0;    //how many nonzero elements we have in a specific column of the upper triangular matrix

    uint32_t **upperVectors = malloc(n*sizeof(uint32_t*));       //array containing the row indices of each nonzero element in each column for the upper triangular matrix
    
    if(upperVectors==NULL){
        printf("Error in COOtoCSC: Couldn't allocate memory for upperVectors");
        exit(-1);
    }

    //Creating the vectors for the upper triangular part of the matrix. The first element of each vector tells us the number of elements of the vector
    for(int32_t i=0; i<n; i++){     
        upperVectors[i]=malloc(sizeof(uint32_t));
        if(upperVectors[i]==NULL){
            printf("Error in COOtoCSC: Couldn't allocate memory for upperVectors[%d]", i);
            exit(-1);
        }
        upperVectors[i][0]=1;  
    }

    uint32_t* colVector2=malloc((n+1)*sizeof(uint32_t));   //index of the elements which start a column of A for the lower triangular part of the matrix
    uint32_t* rowVector2=malloc(nz*sizeof(uint32_t));      //row indices of each non zero element for the lower triangular part of the matrix
    colVector2[0]=0; 

    if(colVector2==NULL){
        printf("Error in COOtoCSC: Couldn't allocate memory for colVector");
        exit(-1);
    }

    if(rowVector2==NULL){
        printf("Error in COOtoCSC: Couldn't allocate memory for rowVector");
        exit(-1);
    }

    //The loop that will fill out rowVector2, colVector2 and upperVectors
    for(int32_t i=0; i<nz; i++){ 

        rowIndiceRead=coo_row[i]-1;
        colIndiceRead=coo_col[i]-1;

        //Check for nonzero elements in the main diagonal
        if(rowIndiceRead == colIndiceRead){
            printf("There are elements in the main diagonal. Please give me an mtx without elements in the diagonal.\n");
            printf("The row of the first element in the diagonal is the %d\n",rowIndiceRead);
            exit(-1);
        }

        rowVector2[i]=rowIndiceRead;
        
        if(colCheck<colIndiceRead){             //Check if the column indice we got is bigger than the previous column indice we examined
            if(colIndiceRead-colCheck>1){       //Checking if one or more consecutive columns have only elements equal to zero
                elemsUntilZeroCols = colVector2[colIndex-1] + lowerColElements;
                for (int32_t k=0; k<(colIndiceRead-colCheck); k++){
                    colVector2[colIndex]=elemsUntilZeroCols;   //For all these all-zero columns put in the respective column array the value of the total elements up until that point
                    colIndex++;
                }
                lowerColElements=0;
            }
            //If the current element is in a different column than the previous one
            else{ 
                colVector2[colIndex]=colVector2[colIndex-1] + lowerColElements;
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
    colVector2[colIndex]=colVector2[colIndex-1]+lowerColElements;  


    //Filling the last values, in case the last columns are all-zero columns 
    while(colIndex<n){           
        colIndex++;
        colVector2[colIndex]=colVector2[colIndex-1];
    }  

    uint32_t* finalRowVector = calloc(2*nz,sizeof(uint32_t));   //row indices of each non zero element for the whole matrix
    uint32_t rowVectorCount=0;                                  //shows how many row indices we have added in the finalRowVector
    uint32_t* finalColVector = calloc((n+1),sizeof(uint32_t));  //index of the elements which start a column of the whole matrix
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
    for(int32_t i=0; i<n; i++){
        
        //Getting the values of the upper triangular half of the matrix
        //Check for nonzero elements in the column
        if(upperVectors[i][0]>1){
            //Add these elements on the final row vector
            for(int32_t j=0; j<upperVectors[i][0]-1; j++){
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
        if(i<n-1){
            cscInitialColElems=colVector2[i+1] - colVector2[i];
            //Add these elements on the final row vector
            for (int32_t j=0; j<cscInitialColElems; j++){             
                finalRowVector[rowVectorCount] = rowVector2[colVector2[i] +j];
                rowVectorCount ++;
                colVectorCount ++;
            }
        }
        finalColVector[i+1] = finalColVector[i] + colVectorCount;
        colVectorCount=0;
    } 

    int32_t elemsInCol;                    //Number on nonzero elements in a particular column
    vector = calloc(n, sizeof(uint32_t));   //Each entry contains the number of triangles in which at least an element of this column belongs to
    if(vector==NULL){
        printf("Error in main: Couldn't allocate memory for vector");
        exit(-1);
    }

    uint32_t element1;   //First common idice we investigate
    uint32_t element2;   //Second common indice we investigate

    //We implement the same algorithm used in V3serial
    for(int32_t i=0; i<n; i++){
        elemsInCol = finalColVector[i+1]- finalColVector[i];
        //Check for every pair of the column with this double for loop
        for(int32_t j=0; j<elemsInCol-1; j++){
            element1 = finalRowVector[finalColVector[i]+j];
            for (int32_t k=j+1; k<elemsInCol; k++ ){
                element2 = finalRowVector[finalColVector[i]+k];          
                //Check if the third common indice exists
                if (elementInColumnCheck2(finalRowVector, finalColVector, element1, element2)>=0){
                    vector[i]++;
                }      
            }
        }  
    }

    free(upperVectors);
    free(colVector2);
    free(rowVector2);
    free(finalColVector);
    free(finalRowVector);

    return vector;
}

/**
 * This function is used as a test to see whether the calculation of the vector containing the number of triangles each
 * node is adjacent to is correct. The function calculates this vector with a method that is verifiably correct. Then,
 * each element of the given vector is compared with the correspondent element of the created vector. If all correspondent
 * elements are equal, then the given vector is correct and the calculation of triangles in our serial or parallel
 * method is correct.
 * Input:
 *      uint32_t* vector1: the vector that is created by the serial or parallel method and we want to see if it is correct
 *      char* filename: string containing the name of the .mtx in which the sparse matrix is contained in matrix market format
 * Output:
 *      uint32_t result: 1 if vector1 is correct, 0 otherwise
**/

uint32_t checkCorrectness(uint32_t* vector1, char* filename){
    uint32_t result = 1;    //returned value, 1 if vector1 is correct, 0 otherwise

    FILE *stream;       //file pointer to read the given file
    MM_typecode t;      //the typecode struct  

    //Opening The file as shown in the command line
    stream=fopen(filename, "r");        
    if(stream==NULL){
        printf("Could not open file, pass another file\n");
        exit(-1);
    }

    mm_read_banner(stream,&t);

    COOArray* cooArray = createCOO(stream); //the lower triangular part of the sparse matrix in coo format

    fclose(stream);

    //Create the verifiably correct vector
    uint32_t* vector2 = vertexWiseTriangleCounts(cooArray->coo_row, cooArray->coo_col, cooArray->n, cooArray->nz);
    
    //Check the two vectors
    for(uint32_t i=0;i<cooArray->n;++i){
        if(vector1[i]!=vector2[i]){
            result=0;
            break;
        }
    }

    free(cooArray->coo_row);
    free(cooArray->coo_col);
    free(cooArray);
    free(vector2);

    return result;
}