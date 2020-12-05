/**
 * This program takes a CSC structure and calculates the amount of triangles for each node 
 * implementing the serial version of V3 algorithm.
**/

#include <stdio.h>
#include <time.h>
#include <stdint.h>
#include "tester.c"


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

int32_t elementInColumnCheck(uint32_t* rowVector,uint32_t* colVector, uint32_t colNum, uint32_t wantedRow){
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

int main(int argc, char* argv[]){
    printf("\nStarted V3serial\n");
    FILE *stream;       //file pointer to read the given file
    MM_typecode t;      //the typecode struct
    
    if(argc<2){
        printf("Please pass as argument the .mtx file\n");
        exit(-1);
    }  

    char* s=argv[1];

    //Checking if the argument is .mtx file
    uint32_t nameLength = strlen(s);    //length of the name of the file    
    if(!((s[nameLength-1]=='x') && (s[nameLength-2]=='t') && (s[nameLength-3]=='m') && (s[nameLength-4]=='.'))){
        printf("Your argument is not an .mtx file\n");
        exit(-1);
    }

    //Opening The file as shown in the command line
    stream=fopen(s, "r");        
    if(stream==NULL){
        printf("Could not open file, pass another file\n");
        exit(-1);
    }

    mm_read_banner(stream,&t);

    //Checking if the matrix type is ok
    if (mm_is_sparse(t)==0){
        printf("The array is not sparce. Please give me another matrix market file\n");
        exit(-1);
    }
    if (mm_is_coordinate(t)==0){
        printf("The array is not in coordinate format. Please give me another matrix market file\n");
        exit(-1);
    }
    if (mm_is_symmetric(t)==0){
        printf("The array is not symmetric. Please give me another matrix market file\n");
        exit(-1);
    }

    CSCArray* cscArray = COOtoCSC(stream);  //The sparse array in csc format
    
    fclose(stream);

    uint32_t* rowVector = cscArray->rowVector;
    uint32_t* colVector = cscArray->colVector;
    uint32_t M = cscArray->M;

    int32_t elemsInCol;                                    //Number on nonzero elements in a particular column
    uint32_t* triangleCount = calloc(M, sizeof(uint32_t));  //Each entry contains the number of triangles in which at least an element of this column belongs to
    if(triangleCount==NULL){
        printf("Error in main: Couldn't allocate memory for triangleCount");
        exit(-1);
    }     

    uint32_t element1;   //First common idice we investigate
    uint32_t element2;   //Second common indice we investigate

    //Start timer
    struct timespec init;
    clock_gettime(CLOCK_MONOTONIC, &init); //CLOCK_MONOTONIC might appear underlined as undefined but it's not an actual problem

    /**
     * The following part of the code calculates the number of triangles adjacent to each node and ultimately the total number of triangles of the sparse matrix
     * Algorithm works as follows:
     * First of we take each pair of nonzero elements that belong in the same column i, let's say elements (j,i), (k,i). There is no need to check 
     * the pair (k,i), (j,i) since we would calculate the same triangle twice this way. So we examine only the cases where k>j.
     * If the element (j,k) is a non zero element then we have a triangle. We only increase the value of triangleCount[i] because if we increased
     * triangleCount[element1] and triangleCount[element2] we would add the same triangle three times instead of one.
    **/

    for(int32_t i=0; i<M; i++){
        elemsInCol = colVector[i+1]- colVector[i];
        //Check for every pair of the column with this double for loop
        for(int32_t j=0; j<elemsInCol-1; j++){
            element1 = rowVector[colVector[i]+j];
            for (int32_t k=j+1; k<elemsInCol; k++ ){
                element2 = rowVector[colVector[i]+k];          
                //Check if the third common indice exists
                if (elementInColumnCheck(rowVector, colVector, element1, element2)>=0){
                    triangleCount[i]++;
                }      
            }
        }  
    }

    //End timer
    struct timespec last;   
    clock_gettime(CLOCK_MONOTONIC, &last);

    long ns;
    uint32_t seconds;
    if(last.tv_nsec <init.tv_nsec){
        ns=init.tv_nsec - last.tv_nsec;
        seconds= last.tv_sec - init.tv_sec -1;
    }

    if(last.tv_nsec >init.tv_nsec){
        ns= last.tv_nsec -init.tv_nsec ;
        seconds= last.tv_sec - init.tv_sec ;
    }
    printf("For V3serial the seconds elapsed are %u and the nanoseconds are %ld\n",seconds, ns); 

    CSCArrayfree(cscArray);
    free(cscArray);

    if(checkCorrectness(triangleCount, s)==0){
        printf("Incorrect calculation of triangles\n");
        exit(-1);
    }
    else{
        printf("Correct calculation of triangles\n");
    }
    

    uint32_t totalTriangles=0; //Total number of triangles

    //Compute the total number of triangles
    for (uint32_t i=0; i<M; i++){
        totalTriangles += triangleCount[i];
    }

    //We divide the total number by 3 because now each triangle is added 3 times, 1 for each node in which it is adjacent
    printf("Total triangles = %u\n", totalTriangles/3);

    free(triangleCount);
    
    return 0;
}