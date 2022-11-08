#include <inttypes.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#define MAX 10000

/**
 * Read matrix from file and find its row and col then pass them to
 * some variable.
 *
 * Use fgets and sscanf to scan each line of the file, and count rows
 * and elements in each row (column).
 * Use pointer to pass the row and col that we find from this function
 * for future use.
 *
 * @param  filename is a char array which is the name of file we want to open
 * @param  row is a fixed-size pointer integer, use for storing the row
 * @param  col is a fixed-size pointer integer, use for storing the column
 */
void countRowCol(const char *filename, uint64_t *row, uint64_t *col)
{
    // Read file and store the matrix in double array
    FILE *file_ptr;
    file_ptr = fopen(filename, "r");
    if (file_ptr == NULL)
    {
        perror("Cannot open file");
    }

    // The max line we could read
    char line[MAX];
    double colNum;

    // Record how many space sscanf used
    int buffer;
    while (fgets(line, sizeof(line), file_ptr) != NULL)
    {

        /*
        Scan variables by using %n in sscanf, helps to parse all inputs in the line.
        Modified from:
        https://stackoverflow.com/questions/3975236/how-to-use-sscanf-in-loops
        */
        /*Start*/
        char *tracker = line;
        while (*row == 0 && sscanf(tracker, " %lf%n", &colNum, &buffer) == 1)
        {
            // Count the column from the first row.
            (*col) += 1;
            tracker += buffer;
            // printf("read: %5lf; count = %5" PRId64 "; buffer = %5d\n", colNum, *col, buffer);
        }
        /*End*/
        // If the line is not blank, count the rows
        if (strcmp(line, "\n") != 0)
            (*row)++;
    }

    if (ferror(file_ptr))
        printf("\nEncountered an error while reading the file %s.", filename);

    fclose(file_ptr);
}

/**
 * Read matrix from file and verify if the matrix in the file is valid
 *
 * Verify whether the matrix has same amount of numbers in each row, and
 * whether the matrix contain some non-numerical items.
 *
 * @param  filename is a char array which is the name of file we want to open
 * @return an integer, if it's 1, then the matrix is invalid and will report the
 *         error message. If it's 0, then this matrix is valid.
 */
int verifyMatrix(const char *filename)
{
    uint64_t col = 0;
    uint64_t row = 0;
    countRowCol(filename, &row, &col);

    // "count" array for counting number of double numbers in each line
    uint64_t *count = calloc(row, sizeof(*count));
    // "count2" array for counting number of items in each line, including non-numerical
    uint64_t *count2 = calloc(row, sizeof(*count2));

    // Read file and store the matrix in double array
    FILE *file_ptr;
    file_ptr = fopen(filename, "r");
    if (file_ptr == NULL)
    {
        perror("Cannot open file");
    }

    // The max line we could read
    char line[MAX];
    double colNum;
    // Store the string fragment by strtok
    char *token;

    // Record how many space sscanf used
    int buffer;
    uint64_t i = 0;
    // Check if every row has same elements
    while (fgets(line, sizeof(line), file_ptr) != NULL)
    {
        /*
        Scan variables by using %n in sscanf, helps to parse
        all inputs in the line.
        Modified from:
        https://stackoverflow.com/questions/3975236/how-to-use-sscanf-in-loops
        */
        /*Start*/
        char *tracker = line;
        while (sscanf(tracker, "%lf%n", &colNum, &buffer) == 1)
        {
            // Count the amount of numerical items in row i.
            count[i] += 1;
            tracker += buffer;
            // printf("read: %5lf; count = %5" PRId64 "; buffer = %5d\n", colNum, count[i], buffer);
        }
        /*End*/

        /*
        Using strtok to find out all items separate by whitespace, exclude newline
        Modified from:
        https://www.tutorialspoint.com/c_standard_library/c_function_strtok.htm
        */
        /*Start*/
        token = strtok(line, " \n");
        while (token != NULL)
        {
            // Count the amount of any type items in row i.
            count2[i] += 1;
            token = strtok(NULL, " \n");
        }
        /*End*/

        // If the line is not blank, count the rows
        if (strcmp(line, "\n") != 0)
            i++;
    }

    if (ferror(file_ptr))
        printf("\nEncountered an error while reading the file %s.", filename);

    fclose(file_ptr);

    // Finding out if every row contains same amount of items.
    // Compare each two element in "count",
    // if one of them not equal to each other, return 1.
    for (uint64_t j = 0; j < row - 1; j++)
    {
        if (count[j] != count[j + 1])
        {
            printf("Matrix in %s doesn't have same amount of numbers between row %" PRId64 " and row %" PRId64, filename, j + 1, j + 2);
            free(count);
            free(count2);
            return 1;
        }
    }

    // Finding out whether all items are numerical.
    // Compare total item's quantity with total numerical item's quantity,
    // if mismatch, return 1.
    for (uint64_t j = 0; j < row; j++)
    {
        // printf("Count 1[%" PRId64 "]: % " PRId64 " ", j, count[j]);
        // printf("Count 2[%" PRId64 "]: % " PRId64 "\n", j, count2[j]);
        if (count[j] != count2[j])
        {
            printf("Matrix in %s has invalid input in row %" PRId64, filename, j + 1);
            free(count);
            free(count2);
            return 1;
        }
    }

    free(count);
    free(count2);

    // If the matrix is valid, return 0.
    return 0;
}

/**
 * Read matrices from file 1 and file 2 and add them up.
 *
 * Store matrices in dynamic arrays. If matrices are in the same dimension,
 * add them up and store the result in a new dynamic array.
 *
 * @param filename1 is a char array which is the name of the first matrix file
 * @param filename2 is a char array which is the name of the second matrix file
 * @param storeFile is a char array which is the name of the file we store the result
 */
void add(const char *filename1, const char *filename2, const char *storeFile)
{
    // Check if all the matrices from the files are valid
    if (verifyMatrix(filename1) == 1 || verifyMatrix(filename2) == 1)
    {
        return;
    }

    uint64_t col1 = 0;
    uint64_t row1 = 0;
    uint64_t col2 = 0;
    uint64_t row2 = 0;
    countRowCol(filename1, &row1, &col1);
    countRowCol(filename2, &row2, &col2);
    if (row1 != row2 || col1 != col2)
    {
        printf("Two matrices don't have the same dimension, cannot be added together.\n");
        return;
    }

    // Dynamically create 2d arrays for store matrices
    double *m1 = calloc(row1 * col1, sizeof(*m1));
    if (m1 == NULL)
    {
        printf("Failed to allocate memory for matrix 1!\n");
        return;
    }

    double *m2 = calloc(row2 * col2, sizeof(*m2));
    if (m2 == NULL)
    {
        printf("Failed to allocate memory for matrix 2!\n");
        return;
    }

    double *sum = calloc(row1 * col1, sizeof(*sum));
    if (sum == NULL)
    {
        printf("Failed to allocate memory for sum result!\n");
        return;
    }

    // Read the first file and store the matrix in double array
    FILE *file_ptr;
    file_ptr = fopen(filename1, "r");
    if (file_ptr == NULL)
    {
        perror("Cannot open file");
    }
    for (uint64_t i = 0; i < row1; i++)
    {
        for (uint64_t j = 0; j < col1; j++)
        {
            fscanf(file_ptr, "%lf", &m1[col1 * i + j]);
        }
    }

    if (ferror(file_ptr))
        printf("\nEncountered an error while reading the file %s.", filename1);

    fclose(file_ptr);

    // Read the second file and store the matrix in double array
    file_ptr = fopen(filename2, "r");
    if (file_ptr == NULL)
    {
        perror("Cannot open file");
    }
    for (uint64_t i = 0; i < row1; i++)
    {
        for (uint64_t j = 0; j < col1; j++)
        {
            fscanf(file_ptr, "%lf", &m2[col2 * i + j]);
        }
    }

    if (ferror(file_ptr))
        printf("\nEncountered an error while reading the file %s.", filename2);
    fclose(file_ptr);

    for (uint64_t i = 0; i < row1; i++)
    {
        for (uint64_t j = 0; j < col1; j++)
        {
            sum[col1 * i + j] = m1[col1 * i + j] + m2[col2 * i + j];
            // printf("%f ", sum[i][j]);
        }
    }

    // Write the result in a new file
    FILE *file_ptr_w;
    file_ptr_w = fopen(storeFile, "w");
    if (file_ptr_w == NULL)
    {
        perror("Cannot open file");
    }
    for (uint64_t i = 0; i < row1; i++)
    {
        for (uint64_t j = 0; j < col1; j++)
        {
            fprintf(file_ptr_w, "%lf ", sum[col1 * i + j]);
        }
        fprintf(file_ptr_w, "\n");
    }

    if (ferror(file_ptr_w))
        printf("\nEncountered an error while writing the file %s.", storeFile);
    fclose(file_ptr_w);

    // De-allocate the memory

    free(m1);
    free(m2);
    free(sum);
}

/**
 * Read matrices from file 1 and file 2 and use first matrix subtract the second one.
 *
 * Store matrices in dynamic arrays. If matrices are in the same dimension,
 * subtract each other and store the result in a new dynamic array.
 *
 * @param filename1 is a char array which is the name of the first matrix file
 * @param filename2 is a char array which is the name of the second matrix file
 * @param storeFile is a char array which is the name of the file we store the result
 */
void subtract(const char *filename1, const char *filename2, const char *storeFile)
{
    // Check if all the matrices from the files are valid
    if (verifyMatrix(filename1) == 1 || verifyMatrix(filename2) == 1)
    {
        return;
    }
    uint64_t col1 = 0;
    uint64_t row1 = 0;
    uint64_t col2 = 0;
    uint64_t row2 = 0;
    countRowCol(filename1, &row1, &col1);
    countRowCol(filename2, &row2, &col2);
    if (row1 != row2 || col1 != col2)
    {
        printf("Two matrices don't have the same dimension, cannot be subtracted together.\n");
        return;
    }

    // Dynamically create 2d arrays for store matrices
    double *m1 = calloc(row1 * col1, sizeof(*m1));
    if (m1 == NULL)
    {
        printf("Failed to allocate memory for matrix 1!\n");
        return;
    }

    double *m2 = calloc(row2 * col2, sizeof(*m2));

    if (m2 == NULL)
    {
        printf("Failed to allocate memory for matrix 2!\n");
        return;
    }

    double *result = calloc(row1 * col1, sizeof(*result));
    ;
    if (result == NULL)
    {
        printf("Failed to allocate memory for subtraction result!\n");
        return;
    }

    // Read the first file and store the matrix in double array
    FILE *file_ptr;
    file_ptr = fopen(filename1, "r");
    if (file_ptr == NULL)
    {
        perror("Cannot open file");
    }
    for (uint64_t i = 0; i < row1; i++)
    {
        for (uint64_t j = 0; j < col1; j++)
        {
            fscanf(file_ptr, "%lf", &m1[col1 * i + j]);
        }
    }

    if (ferror(file_ptr))
        printf("\nEncountered an error while reading the file %s.", filename1);
    fclose(file_ptr);

    // Read the second file and store the matrix in double array
    file_ptr = fopen(filename2, "r");
    if (file_ptr == NULL)
    {
        perror("Cannot open file");
    }
    for (uint64_t i = 0; i < row1; i++)
    {
        for (uint64_t j = 0; j < col1; j++)
        {
            fscanf(file_ptr, "%lf", &m2[col2 * i + j]);
        }
    }

    if (ferror(file_ptr))
        printf("\nEncountered an error while reading the file %s.", filename2);
    fclose(file_ptr);

    for (uint64_t i = 0; i < row1; i++)
    {
        for (uint64_t j = 0; j < col1; j++)
        {
            result[col1 * i + j] = m1[col1 * i + j] - m2[col2 * i + j];
            // printf("%f ", result[col1*i+j]);
        }
    }

    // Write the result in a new file
    FILE *file_ptr_w;
    file_ptr_w = fopen(storeFile, "w");
    if (file_ptr_w == NULL)
    {
        perror("Cannot open file");
    }
    for (uint64_t i = 0; i < row1; i++)
    {
        for (uint64_t j = 0; j < col1; j++)
        {
            fprintf(file_ptr_w, "%lf ", result[col1 * i + j]);
        }
        fprintf(file_ptr_w, "\n");
    }

    if (ferror(file_ptr_w))
        printf("\nEncountered an error while writing the file %s.", storeFile);
    fclose(file_ptr_w);

    // De-allocate the memory
    free(m1);
    free(m2);
    free(result);
}

/**
 * Read matrices from file 1 and file 2 and multiply them up.
 *
 * Store matrices in dynamic arrays. If matrix one's column number is same as
 * matrix two's row number, multiply them up and store the result in a new
 * dynamic array. Then write the array in a new file.
 *
 * @param filename1 is a char array which is the name of the first matrix file
 * @param filename2 is a char array which is the name of the second matrix file
 * @param storeFile is a char array which is the name of the file we store the result
 */
void multiply(const char *filename1, const char *filename2, const char *storeFile)
{
    // Check if all the matrices from the files are valid
    if (verifyMatrix(filename1) == 1 || verifyMatrix(filename2) == 1)
    {
        return;
    }
    uint64_t col1 = 0;
    uint64_t row1 = 0;
    uint64_t col2 = 0;
    uint64_t row2 = 0;
    countRowCol(filename1, &row1, &col1);
    countRowCol(filename2, &row2, &col2);
    if (col1 != row2)
    {
        printf("Two matrices don't have the right dimension for doing multiplication.\n");
        return;
    }

    // Dynamically create 2d arrays for store matrices
    double *m1 = calloc(row1 * col1, sizeof(*m1));
    if (m1 == NULL)
    {
        printf("Failed to allocate memory for matrix 1!\n");
        return;
    }

    double *m2 = calloc(row2 * col2, sizeof(*m1));
    if (m2 == NULL)
    {
        printf("Failed to allocate memory for matrix 2!\n");
        return;
    }

    double *result = calloc(row1 * col2, sizeof(*result));
    if (result == NULL)
    {
        printf("Failed to allocate memory for result!\n");
        return;
    }

    // Read the first file and store the matrix in double array
    FILE *file_ptr;
    file_ptr = fopen(filename1, "r");
    if (file_ptr == NULL)
    {
        perror("Cannot open file");
    }
    for (uint64_t i = 0; i < row1; i++)
    {
        for (uint64_t j = 0; j < col1; j++)
        {
            fscanf(file_ptr, "%lf", &m1[col1 * i + j]);
        }
    }

    if (ferror(file_ptr))
        printf("\nEncountered an error while reading the file %s.", filename1);
    fclose(file_ptr);

    // Read the second file and store the matrix in double array
    file_ptr = fopen(filename2, "r");
    if (file_ptr == NULL)
    {
        perror("Cannot open file");
    }
    for (uint64_t i = 0; i < row2; i++)
    {
        for (uint64_t j = 0; j < col2; j++)
        {
            fscanf(file_ptr, "%lf", &m2[col2 * i + j]);
        }
    }

    if (ferror(file_ptr))
        printf("\nEncountered an error while reading the file %s.", filename2);
    fclose(file_ptr);

    for (uint64_t i = 0; i < row1; i++)
    {
        for (uint64_t j = 0; j < col2; j++)
        {
            result[col2 * i + j] = 0;
            for (uint64_t n = 0; n < row2; n++)
            {
                // Dot product rows of first matrix with cols of second matrix and add them up
                result[col2 * i + j] += m1[col1 * i + n] * m2[col2 * n + j];
                // printf("%lf ", result[i][j]);
            }
        }
    }

    // Write the result in a new file
    FILE *file_ptr_w;
    file_ptr_w = fopen(storeFile, "w");
    if (file_ptr_w == NULL)
    {
        perror("Cannot open file");
    }
    for (uint64_t i = 0; i < row1; i++)
    {
        for (uint64_t j = 0; j < col2; j++)
        {
            fprintf(file_ptr_w, "%lf ", result[col2 * i + j]);
        }
        fprintf(file_ptr_w, "\n");
    }

    if (ferror(file_ptr_w))
        printf("\nEncountered an error while writing the file %s.", storeFile);
    fclose(file_ptr_w);

    // De-allocate the memory
    free(m1);
    free(m2);
    free(result);
}

/**
 * Read matrix from file and multiply it by scalar.
 *
 * Store matrix in dynamic arrays. Multiply each element of matrix
 * by the scalar number and store the result in a new dynamic array.
 * Then store the array in a new file.
 *
 * @param filename is a char array which is the name of the first matrix file
 * @param scalar is a char array which is a double number in string form.
 * @param storeFile is a char array which is the name of the file we store the result
 */
void multiply_scalar(const char *filename, char *scalar, const char *storeFile)
{
    // Check if all the matrices from the files are valid
    if (verifyMatrix(filename) == 1)
    {
        return;
    }
    uint64_t col = 0;
    uint64_t row = 0;
    countRowCol(filename, &row, &col);

    // Dynamically create 2d arrays for store matrices
    double *m = calloc(row * col, sizeof(*m));
    if (m == NULL)
    {
        printf("Failed to allocate memory for matrix!\n");
        return;
    }

    double *result = calloc(row * col, sizeof(*result));

    if (result == NULL)
    {
        printf("Failed to allocate memory for result!\n");
        return;
    }

    // Read file and store the matrix in double array
    FILE *file_ptr;
    file_ptr = fopen(filename, "r");
    if (file_ptr == NULL)
    {
        perror("Cannot open file");
    }
    for (uint64_t i = 0; i < row; i++)
    {
        for (uint64_t j = 0; j < col; j++)
        {
            fscanf(file_ptr, "%lf", &m[col * i + j]);
            // printf("i is %" PRId64 " result: %lf ", i, m[col*i+j]);
        }
    }

    if (ferror(file_ptr))
        printf("\nEncountered an error while reading the file %s.", filename);
    fclose(file_ptr);

    // Convert the string argument to double for calculation
    // If the input is not a number, print out the error message and return
    double s;
    if (sscanf(scalar, "%lf", &s) != 1)
    {
        printf("Please use a valid number for the scalar.");
        return;
    }

    for (uint64_t i = 0; i < row; i++)
    {
        for (uint64_t j = 0; j < col; j++)
        {
            result[col * i + j] = m[col * i + j] * s;
        }
    }

    // Write the result in a new file
    FILE *file_ptr_w;
    file_ptr_w = fopen(storeFile, "w");
    if (file_ptr_w == NULL)
    {
        perror("Cannot open file");
    }
    for (uint64_t i = 0; i < row; i++)
    {
        for (uint64_t j = 0; j < col; j++)
        {
            fprintf(file_ptr_w, "%lf ", result[col * i + j]);
        }
        fprintf(file_ptr_w, "\n");
    }

    if (ferror(file_ptr_w))
        printf("\nEncountered an error while writing the file %s.", storeFile);
    fclose(file_ptr_w);

    // De-allocate the memory

    free(m);
    free(result);
}

/**
 * Read matrices from file 1 and file 2 and test whether they are equal to each other.
 *
 * Store matrices in dynamic arrays. If matrices are in the same dimension,
 * compare their elements one by one. Print out the message if or not these
 * two matrices are equal.
 *
 * @param filename1 is a char array which is the name of the first matrix file
 * @param filename2 is a char array which is the name of the second matrix file
 */
void is_equal(const char *filename1, const char *filename2)
{
    // Check if all the matrices from the files are valid
    if (verifyMatrix(filename1) == 1 || verifyMatrix(filename2) == 1)
    {
        return;
    }
    uint64_t col1 = 0;
    uint64_t row1 = 0;
    uint64_t col2 = 0;
    uint64_t row2 = 0;
    countRowCol(filename1, &row1, &col1);
    countRowCol(filename2, &row2, &col2);

    if (row1 != row2 || col1 != col2)
    {
        printf("Two matrices don't have the same dimension, cannot compare equality together.\n");
        return;
    }

    // Dynamically create 2d arrays for store matrices
    double *m1 = calloc(row1 * col1, sizeof(*m1));
    if (m1 == NULL)
    {
        printf("Failed to allocate memory for matrix 1!\n");
        return;
    }
    double *m2 = calloc(row2 * col2, sizeof(*m2));
    if (m2 == NULL)
    {
        printf("Failed to allocate memory for matrix 2!\n");
        return;
    }

    // Read the first file and store the matrix in double array
    FILE *file_ptr;
    file_ptr = fopen(filename1, "r");
    if (file_ptr == NULL)
    {
        perror("Cannot open file");
    }
    for (uint64_t i = 0; i < row1; i++)
    {
        for (uint64_t j = 0; j < col1; j++)
        {
            fscanf(file_ptr, "%lf", &m1[col1 * i + j]);
        }
    }

    if (ferror(file_ptr))
        printf("\nEncountered an error while reading the file %s.", filename1);
    fclose(file_ptr);

    // Read the second file and store the matrix in double array
    file_ptr = fopen(filename2, "r");
    if (file_ptr == NULL)
    {
        perror("Cannot open file");
    }
    for (uint64_t i = 0; i < row1; i++)
    {
        for (uint64_t j = 0; j < col1; j++)
        {
            fscanf(file_ptr, "%lf", &m2[col2 * i + j]);
        }
    }

    if (ferror(file_ptr))
        printf("\nEncountered an error while reading the file %s.", filename2);
    fclose(file_ptr);

    for (uint64_t i = 0; i < row1; i++)
    {
        for (uint64_t j = 0; j < col1; j++)
        {
            if (m1[col1 * i + j] != m2[col2 * i + j])
            {
                printf("False, two matrices are not equal to each other.\n");
                // De-allocate the memory
                free(m1);
                free(m2);

                return;
            }
            else
            {
                continue;
            }
        }
    }
    printf("True, two matrices are equal to each other.\n");

    // De-allocate the memory
    free(m1);
    free(m2);
}

/**
 * Read matrices from file and compute its trace.
 *
 * Store matrix in dynamic arrays. If the matrix is a square matrix,
 * add all the diagonal values. Print the result in the terminal.
 *
 * @param filename is a char array which is the name of the matrix file.
 */
void trace(const char *filename)
{
    // Check if all the matrices from the files are valid
    if (verifyMatrix(filename) == 1)
    {
        return;
    }
    uint64_t col = 0;
    uint64_t row = 0;
    countRowCol(filename, &row, &col);

    // Check if the matrix is a square matrix.
    if (row != col)
    {
        printf("The matrix is not a square matrix, cannot compute its trace.");
        return;
    }

    // The size of the array is the dimension of the first matrix
    // uint64_t size = row * col;

    // Dynamically create 2d arrays for store matrices
    double *m = calloc(row * col, sizeof(*m));

    if (m == NULL)
    {
        printf("Failed to allocate memory for matrix!\n");
        return;
    }

    // Read file and store the matrix in double array
    FILE *file_ptr;
    file_ptr = fopen(filename, "r");
    if (file_ptr == NULL)
    {
        perror("Cannot open file");
    }

    for (uint64_t i = 0; i < row; i++)
    {
        for (uint64_t j = 0; j < col; j++)
        {
            fscanf(file_ptr, "%lf", &m[col * i + j]);
            // printf("i is %" PRId64 " result: %lf ", i, m[col * i + j]);
        }
    }
    if (ferror(file_ptr))
        printf("\nEncountered an error while reading the file %s.", filename);
    else if (feof(file_ptr))
        printf("\nSuccessfully read file %s.", filename);
    fclose(file_ptr);

    double trace = 0.0;
    for (uint64_t i = 0; i < row; i++)
    {
        for (uint64_t j = 0; j < col; j++)
        {
            if (i == j)
            {
                trace += m[col * i + j];
            }
            continue;
        }
    }
    printf("The trace of this matrix is %lf.", trace);

    // De-allocate memory
    free(m);
}

/**
 * Given a nxn matrix and its size n, compute the determinant
 *
 * Use recursion for nxn matrix with n>2. The base cases are when n=1 and n=2.
 *
 * @param m is a double array contains the nxn matrix
 * @param n is a fixed-size integer, it is the size of nxn matrix m.
 * @return A double value, that is the determinant of parameter matrix m.
 */
double determinant(double *m, uint64_t n)
{
    double result = 0;
    if (n == 1)
    {
        return m[0];
    }

    if (n == 2)
    {
        result = m[0] * m[3] - m[1] * m[2];
        return result;
    }

    // Choose the first row as the factor

    else
    {
        // Create dynamic array for minor matrix
        double *minor = calloc((n - 1) * (n - 1), sizeof(*minor));
        if (minor == NULL)
        {
            printf("Failed to allocate memory for matrix!\n");
            return 1;
        }
        // Choose the first row as the factors
        uint64_t row = 0;
        // Rows and cols for minor matrix.
        uint64_t x = 0;
        uint64_t y = 0;

        /*
        Construct minor matrix and do the computation for nxn matrix determinant where n>2
        Modified from:
        https://stackoverflow.com/questions/41384020/c-program-to-calculate-the-determinant-of-a-nxn-matrix
        */
        /*Start*/
        for (uint64_t col = 0; col < n; col++)
        {
            for (uint64_t i = 0; i < n; i++)
            {
                for (uint64_t j = 0; j < n; j++)
                {
                    if (i != row && j != col)
                    {

                        minor[(n - 1) * x + y] = m[n * i + j];
                        // printf("%lf", minor[(n - 1) * x + y]);
                        y++;
                        if (y > n - 2)
                        {
                            x++;
                            y = 0;
                        }
                    }
                }
            }
            // Reset minor matrix row counter
            x = 0;
            // Calculate determinant by recursion.
            result = result + pow(-1, (int)col) * (m[n * row + col] * determinant(minor, n - 1));
        } /*End*/
        free(minor);
        return result;
    }
}

/**
 * Read matrices from file and compute its determinant by function determinant
 *
 * Store matrix in dynamic array. If matrix is a square matrix,
 * run function "determinant" to compute determinant. Print the result in terminal.
 *
 * @param filename is a char array which is the name of the matrix file.
 */
void det(const char *filename)
{
    // Check if all the matrices from the files are valid
    if (verifyMatrix(filename) == 1)
    {
        return;
    }
    uint64_t col = 0;
    uint64_t row = 0;
    countRowCol(filename, &row, &col);

    // Check if the matrix is a square matrix.
    if (row != col)
    {
        printf("The matrix is not a square matrix, cannot compute its determinant.");
        return;
    }

    // Dynamically create 2d arrays for store matrices
    double *m = calloc(row * col, sizeof(*m));
    if (m == NULL)
    {
        printf("Failed to allocate memory for matrix!\n");
        return;
    }

    // Read file and store the matrix in double array
    FILE *file_ptr;
    file_ptr = fopen(filename, "r");
    if (file_ptr == NULL)
    {
        perror("Cannot open file");
    }
    for (uint64_t i = 0; i < row; i++)
    {
        for (uint64_t j = 0; j < col; j++)
        {
            fscanf(file_ptr, "%lf", &m[col * i + j]);
        }
    }

    if (ferror(file_ptr))
        printf("\nEncountered an error while reading the file %s.", filename);
    fclose(file_ptr);

    double result = determinant(m, row);
    free(m);
    printf("The determinant is %lf", result);
}

/**
 * Read matrices from file and compute its power by a given exponent
 *
 * Store matrices in dynamic arrays. If matrix is a square matrix,
 * multiply itself multiple times. Store the result in a file.
 *
 * @param filename is a char array which is the name of the matrix file.
 * @param factor is a char array. It's the exponent, an integer in string form.
 * @param storeFile is a char array which is the name of the file we store the result.
 */
void power(const char *filename, char *factor, const char *storeFile)
{
    // Check if all the matrices from the files are valid
    if (verifyMatrix(filename) == 1)
    {
        return;
    }
    uint64_t col = 0;
    uint64_t row = 0;
    countRowCol(filename, &row, &col);

    // Check if the matrix is a square matrix.
    if (row != col)
    {
        printf("The matrix is not a square matrix, cannot compute its power.");
        return;
    }

    // Dynamically create 2d arrays for store matrices
    double *m = calloc(row * col, sizeof(*m));
    if (m == NULL)
    {
        printf("Failed to allocate memory for matrix!\n");
        return;
    }

    // Save the intermediate result in a matrix.
    double *save = calloc(row * col, sizeof(*save));
    if (save == NULL)
    {
        printf("Failed to allocate memory for matrix!\n");
        return;
    }

    // Double variable for holding the intermediate value.
    double holder;

    double *result = calloc(row * col, sizeof(*result));
    if (result == NULL)
    {
        printf("Failed to allocate memory for result!\n");
        return;
    }

    // Read file and store the matrix in double array
    FILE *file_ptr;
    file_ptr = fopen(filename, "r");
    if (file_ptr == NULL)
    {
        perror("Cannot open file");
    }
    for (uint64_t i = 0; i < row; i++)
    {
        for (uint64_t j = 0; j < col; j++)
        {
            fscanf(file_ptr, "%lf", &m[col * i + j]);
            // printf("i is %" PRId64 " result: %lf \n", i, m[i][j]);
        }
    }

    if (ferror(file_ptr))
        printf("\nEncountered an error while reading the file %s.", filename);
    fclose(file_ptr);

    // Convert the input string to integer.
    // If the input is not an integer, print out error message and return.
    int n;
    if (sscanf(factor, "%d", &n) != 1)
    {
        printf("Please use a valid integer for the exponent.");
        return;
    }

    // If the exponent is 0, return an identity matrix.
    if (n == 0)
    {
        for (uint64_t i = 0; i < row; i++)
        {
            for (uint64_t j = 0; j < col; j++)
            {
                if (i == j)
                {
                    result[col * i + j] = 1;
                }
                else
                {
                    result[col * i + j] = 0;
                }
            }
        }
    }
    // If the exponent is 0, return the same matrix.
    else if (n == 1)
    {
        for (uint64_t i = 0; i < row; i++)
        {
            for (uint64_t j = 0; j < col; j++)
            {
                result[col * i + j] = m[col * i + j];
            }
        }
    }

    // If the exponent is greater than 1, do the power calculation.
    else if (n >= 2)
    {
        for (uint64_t i = 0; i < row; i++)
        {
            for (uint64_t j = 0; j < col; j++)
            {
                result[col * i + j] = m[col * i + j];
                // printf("i is %" PRId64 " result: %lf \n", i, result[col*i+j]);
            }
        }
        while (n >= 2)
        {
            for (uint64_t i = 0; i < row; i++)
            {
                for (uint64_t j = 0; j < col; j++)
                {
                    for (uint64_t k = 0; k < row; k++)
                    {
                        /*
                        Dot product rows of first matrix with cols of
                        second matrix and add them up.
                        Use holder to add up all the values.
                        */
                        holder += result[col * i + k] * m[col * k + j];
                    }
                    // Then store them in the "save" matrix.
                    save[col * i + j] = holder;
                    // Clear holder for the next round of calculation.
                    holder = 0.0;
                }
            }
            for (uint64_t i = 0; i < row; i++)
            {
                for (uint64_t j = 0; j < col; j++)
                {
                    // After all the multiplication, update values in result.
                    result[col * i + j] = save[col * i + j];
                }
            }

            n--;
        }
    }
    else
    {
        printf("The exponent need to be a non-negative integer.");
        return;
    }

    // Write the result in a new file
    FILE *file_ptr_w;
    file_ptr_w = fopen(storeFile, "w");
    if (file_ptr_w == NULL)
    {
        perror("Cannot open file");
    }
    for (uint64_t i = 0; i < row; i++)
    {
        for (uint64_t j = 0; j < col; j++)
        {
            fprintf(file_ptr_w, "%lf ", result[col * i + j]);
        }
        fprintf(file_ptr_w, "\n");
    }

    if (ferror(file_ptr_w))
        printf("\nEncountered an error while writing the file %s.", storeFile);
    fclose(file_ptr_w);

    // De-allocate the memory
    free(m);
    free(result);
}

/**
 * Main function, use commend line argument run every functions from above
 *
 * @param argc, an integer which counts the number of arguments passed to the program
 * @param argv[], an array of argc pointers to strings, contain the actual arguments passed.
 * @return 0 by default, return 1 if some exception happens.
 */
int main(int argc, char *argv[])
{
    // If running the program without putting any commend line argument,
    // print error message.
    if (argc == 1)
    {
        printf("You need to put some arguments in command line.");
        return 1;
    }

    // If the first command line argument is "add",
    // run add function and store the result in a new file.
    if (strcmp(argv[1], "add") == 0)
    {
        // Check if the number of command line argument is correct.
        if (argc != 5)
        {
            printf("The number of arguments is incorrect.");
            return 1;
        }
        add(argv[2], argv[3], argv[4]);
        return 0;
    }
    // If the first command line argument is "subtract",
    // run subtract function and store the result in a new file.
    else if (strcmp(argv[1], "subtract") == 0)
    {
        // Check if the number of command line argument is correct.
        if (argc != 5)
        {
            printf("The number of arguments is incorrect.");
            return 1;
        }
        subtract(argv[2], argv[3], argv[4]);
        return 0;
    }

    // If the first command line argument is "multiply",
    // run multiply function and store the result in a new file.
    else if (strcmp(argv[1], "multiply") == 0)
    {
        // Check if the number of command line argument is correct.
        if (argc != 5)
        {
            printf("The number of arguments is incorrect.");
            return 1;
        }
        multiply(argv[2], argv[3], argv[4]);
        return 0;
    }

    // If the first command line argument is "multiply_scalar",
    // run multiply by scalar function and store the result in a new file.
    else if (strcmp(argv[1], "multiply_scalar") == 0)
    {
        // Check if the number of command line argument is correct.
        if (argc != 5)
        {
            printf("The number of arguments is incorrect.");
            return 1;
        }
        multiply_scalar(argv[2], argv[3], argv[4]);
        return 0;
    }

    // If the first command line argument is "is_equal",
    // compare two matrices from corresponding files is equal or not.
    else if (strcmp(argv[1], "is_equal") == 0)
    {
        // Check if the number of command line argument is correct.
        if (argc != 4)
        {
            printf("The number of arguments is incorrect.");
            return 1;
        }
        is_equal(argv[2], argv[3]);
        return 0;
    }

    // If the first command line argument is "trace",
    // print out the trace of matrix.
    else if (strcmp(argv[1], "trace") == 0)
    {
        // Check if the number of command line argument is correct.
        if (argc != 3)
        {
            printf("The number of arguments is incorrect.");
            return 1;
        }
        trace(argv[2]);
        return 0;
    }

    // If the first command line argument is "det",
    // print out the determinant of matrix.
    else if (strcmp(argv[1], "det") == 0)
    {
        // Check if the number of command line argument is correct.
        if (argc != 3)
        {
            printf("The number of arguments is incorrect.");
            return 1;
        }
        det(argv[2]);
        return 0;
    }

    // If the first command line argument is "power",
    // run power function and store the result in a new file.
    else if (strcmp(argv[1], "power") == 0)
    {
        // Check if the number of command line argument is correct.
        if (argc != 5)
        {
            printf("The number of arguments is incorrect.");
            return 1;
        }
        power(argv[2], argv[3], argv[4]);
        return 0;
    }

    else
    {
        printf("Please use one of the available functions.");
        return 1;
    }
}