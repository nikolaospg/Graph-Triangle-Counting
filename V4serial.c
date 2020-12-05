/**
 * This program takes a CSC structure and calculates the amount of triangles for each node 
 * implementing the serial version of V4 algorithm.
**/

#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include "tester.c"
#include <stdint.h>


/**  
 *  Function that calculates the dot product between two columns.
 *  Every column corresponds to a specific part of the row_vector. By comparing the two parts, we count how many pairs of equal elements we have
 *  in the two parts, and that is the dot product of the two columns.
 *  Useful because the matrices that we deal with are symmetric, so the product between a row1 and col2 is equivalent to the dot product between
 *  col1 and col2.
 *  Calculating the product between a rowi and colj is useful because that is the A[i][j] element of the product of two matrices.
 *  The algorithm is based on binary search and works as follows:
 *  First of all we check which one of the columns has fewer nonzero elements. Then we iterate through the non zero elements of this column
 *  and for each one of them:
 *  We search on the other column whether there is a nonzero element that has the same row indice as this one. If so, we increase the product by 1.
 *  The search on the column is done using binary search to reduce complexity to O(lgn). Moreover, we iterate through the elements of the column
 *  that has fewer nonzero elements than the other in order to decrease the times we loop. In overall, the algorithm runs in O(nlgn) asymptotically
 *  where n is the length of the column of the matrix.
 *  Input:
 *      uint32_t* rowVector: the row indices array of the csc format
 *      uint32_t* colVector: the column changes array of the csc format
 *      uint32_t colNum1: the index of the first column
 *      uint32_t colNum2: the index of the second column
 *  Output:
 *      uint32_t result: the result of the dot product between the column with index colNum1 and the column with index colNum2
 **/

uint32_t product(uint32_t* rowVector, uint32_t* colVector, uint32_t colNum1, uint32_t colNum2){
    uint32_t result = 0;     //the dot product of the two columns

    uint32_t smallLength;    //the number of nonzero elements in the column that has fewer nonzero elements than the other
    uint32_t bigLength;      //the number of nonzero elements in the column that has more nonzero elements than the other

    uint32_t rowIndice1;     //the column index of the element examined. We take advantage here of the symmetry of the matrix, so the csc and crs formats are equivalent
    
    uint32_t smallColIndex;  //the column index of the column that has fewer nonzero elements than the other
    uint32_t bigColIndex;    //the column index of the column that has more nonzero elements than the other
    uint32_t initialLeft;   //initial left, the index in rowVector of the first nonzero element belonging in the "big" column

    int32_t left;       //first element of the sub-array in which we do our binary search
    int32_t right;      //last element of the sub-array in which we do our binary search
    int32_t middle;     //middle element in the binary search
    uint32_t flag;

    uint32_t colLength1 = colVector[colNum1+1] - colVector[colNum1]; //number of nonzero elements in column with index colNum1
    uint32_t colLength2=colVector[colNum2+1] - colVector[colNum2];   //number of nonzero elements in column with index colNum2

    //Check which column has fewer nonzero elements 
    if(colLength1<=colLength2){
        smallLength = colLength1;
        bigLength = colLength2;
        smallColIndex = colNum1;
        bigColIndex = colNum2;
    }
    else{
        smallLength = colLength2;
        bigLength = colLength1;
        smallColIndex = colNum2;
        bigColIndex = colNum1;
    }

    initialLeft = colVector[bigColIndex];
    right = colVector[bigColIndex+1] - 1;

    /**
     *  Loop in which we claculate the dot product. For each nonzero element of the "small" column we execute a binary search
     *  on the "big" column to find out whether there is an element with the same row indice. If so we increase the result by 1.
     *  In order to decrease the number of calculations, we take advantage that both of the elements in these two subarrays are row major
     *  sorted:
     *  1)If we find elememnts in both columns with the same row indice, then, in the next iteration, for the binary search we execute there,
     *    we search on the subarray thst starts right after that element. That is because all the row indices of the following elements are going
     *    to be bigger than those in the previous iteration so there is no reason to search through all nonzero elements of the column
     *  2)If on the sub-array in which we execute the binary search we find an element with smaller row indice than the wanted we immediately
     *    set the sub-array for the next binary search to start at least after that element. That is because, if this element is smaller than the one from
     *    the "small" column, then it will definitely be smaller than the next elements of the "small" columnn (because the lements on the "small" column
     *    are row major sorted).
    **/

   //Loop for each nonzero element on the "small" column
    for(uint32_t i=0; i<smallLength; ++i){
        flag = 0;
        rowIndice1 = rowVector[colVector[smallColIndex]+i]; //wanted row
        left = initialLeft;                     //first element of the subarray in which we execute the bianry search
        right = colVector[bigColIndex+1] - 1;   //last element of the subarray in which we execute the bianry search
        middle = (left+right)/2;
        //binary search
        while(left<=right){
            if(rowVector[middle]==rowIndice1){
                flag = 1;   //element found
                initialLeft = middle+1; //sub-array for the binary search for the next element of the "small" column will begin after this element
                break;
            }
            if(rowVector[middle] < rowIndice1){
                left = middle+1;
                initialLeft = left; //sub-array for the binary search for the next element of the "small" column will definitely begin after this element
            }
            else{
                right = middle-1;
            }
            middle = (left+right)/2;
        }
        //If element is found, increase the result by 1 (because we have a 1x1)
        if(flag == 1){
            result++;
        }
    }

    return result;
}



int main(int argc, char* argv[]){
    printf("\nStarted V4serial\n");
    FILE *stream;       //file pointer to read the given file
    MM_typecode t;      //the typecode struct
    
    if(argc<2){
        printf("Please pass as argument the .mtx file\n");
        exit(-1);
    }  
    char* s=argv[1];

    //Checking if the argument is .mtx file
    uint32_t nameLength = strlen(s);        //length of the name of the file
    if(!((s[nameLength-1]=='x') && (s[nameLength-2]=='t') && (s[nameLength-3]=='m') && (s[nameLength-4]=='.'))){
        printf("Your argument is not an .mtx file\n");
        exit(-1);
    }

    //Opening The file as shown in the command line
    stream=fopen(s, "r");        
    if(stream==NULL){
        printf("could not open file, pass another one\n");
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

    CSCArray* cscArray = COOtoCSC(stream); //The sparse array in csc format

    fclose(stream);

    uint32_t* rowVector = cscArray->rowVector;   //the row vector of the sparse matrix in the csc format
    uint32_t* colVector = cscArray->colVector;   //the column vector of the sparse matrix in the csc format
    uint32_t M = cscArray->M;                    //number of columns/rows of the sparse matrix

    /**
     * The following part of the code calculates the number of triangles adjacent to each node and ultimately the total number of triangles of the sparse matrix
     * Algorithm works as it follows:
     * First of all, we check where we have nonzero elements in A. Let's denote such an element with A[i][j]. Then we calculate the dot product of row with index i
     * with the column of index j. That is because, in the rest of the positions of the matrix there is no need to calculate the dot product due to the fact that we
     * will then multiply it with 0 (since we take the Hadamard product of A with AxA).
     * Then, multiplication of c with vector e means that we end up with a Mx1 matrix where each row is the sum of the elements of the corresponding column.
     * For this reason, in our algorithm, after we have computed the dot product of row j with column i we add this to trianglesArray[i].
     * To compute the total number of triangles we have to add all triangles adjacent to each node i and then divide it by 3 because each triangle is calculated 3 times.
    **/

    uint32_t colLength;      //number of nonzero elements of a column
    uint32_t rowNum;         //the row indice
    uint32_t colNum;         //the column indice
    uint32_t productNum;     //dot product of row with index rowNum with column with index colNum
    
    uint32_t* trianglesArray=calloc(M, sizeof(uint32_t)); //array containing the number of triangles adjacent to each node i (M nodes in total)
    if(trianglesArray==NULL){
        printf("Error in main: Couldn't allocate memory for trianglesArray");
        exit(-1);
    }

    //Start timer
    struct timespec init;
    clock_gettime(CLOCK_MONOTONIC, &init);

    //Loop that calculates the number of triangles adjacent to each node
    //Check for each column
    for(uint32_t i=0; i<M; i++){ 
        colLength = colVector[i+1] - colVector[i];
        colNum=i;
        //Checking for each nonzero element of the column
        for(uint32_t j=0; j<colLength; j++){
            rowNum = rowVector[colVector[i]+j];
            //If this row contains only zeros, skip it. We take advantage of the fact that the row with index rowNum is the same with the column with index rowNum
            if(colVector[rowNum+1]-colVector[rowNum] == 0){
                continue;
            }
            productNum=product(rowVector, colVector, colNum, rowNum);
            trianglesArray[i] += productNum;
        }
        trianglesArray[i] /= 2;
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
    printf("For V4serial the seconds elapsed are %u and the nanoseconds are %ld\n",seconds, ns);

    if(checkCorrectness(trianglesArray, s)==0){
        printf("Incorrect calculation of triangles\n");
        exit(-1);
    }
    else{
        printf("Correct calculation of triangles\n");
    }

    uint32_t totalTriangles=0; //total number of triangles

    //Calculate the total number of triangles
    for(uint32_t i=0;i<M;++i){
        totalTriangles += trianglesArray[i];
    }
    totalTriangles /= 3;

    free(trianglesArray);

    free(colVector);
    free(rowVector);
    free(cscArray);

    printf("Total triangles = %d\n", totalTriangles);

    return 0;
}