#include <stdio.h>
#include <string.h>
#include "parse_csv.h"
#include "return_codes.h"
#include "test_harness.h"


#define NUM_TESTS 7
#define NUM_RECORDS 5
#define NUM_COLUMNS 4


/*Test data that will be written to file and tested against.*/
static char * path = "test-csv";
static char long_token[1031];
static char * tokens[NUM_RECORDS*NUM_COLUMNS] = {
    "columna", "columnb", "columnc", "columnd",
    "23", "5", "0", "1",
    "a", "t", "2", "c",
    "j", "-", "a", ":",
    "alsdjf", "bealdkj", "eokej", "aldkfjd"
};
static char * headers_only_tokens[NUM_RECORDS*NUM_COLUMNS] = {
    "columna", "columnb", "columnc", "columnd",
    NULL, NULL, NULL, NULL,
    NULL, NULL, NULL, NULL,
    NULL, NULL, NULL, NULL,
    NULL, NULL, NULL, NULL
};
static char * blank_file_tokens[NUM_RECORDS*NUM_COLUMNS] = {
    NULL, NULL, NULL, NULL,
    NULL, NULL, NULL, NULL,
    NULL, NULL, NULL, NULL,
    NULL, NULL, NULL, NULL,
    NULL, NULL, NULL, NULL
};
static char * blank_line_tokens[NUM_RECORDS*NUM_COLUMNS] = {
    "columna", "columnb", "columnc", "columnd",
    "23", "5", "0", "1",
    "\n", "t", "2", NULL,
    NULL, NULL, NULL, NULL,
    NULL, NULL, NULL, NULL,
};
static char * too_long_line_tokens[NUM_RECORDS*NUM_COLUMNS] = {
    "columna", "columnb", "columnc", "columnd",
    "23", "5", "0", "1",
    (char *)(&long_token), "t", "2", "c",
    "j", "-", "a", ":",
    "alsdjf", "bealdkj", "eokej", "aldkfjd"
};
static char * wrong_columns_tokens[NUM_RECORDS*NUM_COLUMNS] = {
    "columna", "columnb", "columnc", "columnd",
    "23", "5", "\n", "1",
    NULL, NULL, NULL, NULL,
    NULL, NULL, NULL, NULL,
    NULL, NULL, NULL, NULL,
};
static char * token_too_long_tokens[NUM_RECORDS*NUM_COLUMNS] = {
    "columna", "columnb", "columnc", "columnd",
    "23", "5", "000000000000000000000000000000000000000000000000001", "1",
    "a", "t", "2", "c",
    "j", "-", "a", ":",
    "alsdjf", "bealdkj", "eokej", "aldkfjd"
};


/*Writes a test csv file.*/
static void write_test_file(char ** tokens)
{
    FILE * test_file;
    test_file = fopen(path, "w");
    int i;
    for (i=0; i<NUM_RECORDS; ++i)
    {
        int j;
        for (j=0; j<NUM_COLUMNS; ++j)
        {
            if (tokens[i*NUM_COLUMNS + j] == NULL)
            {
                fclose(test_file);
                return;
            }
            if ((void *)(tokens[i*NUM_COLUMNS + j]) == (void *)(&long_token))
            {
                memset(long_token, 'a', sizeof(long_token));
                fprintf(test_file, "%s", long_token);
            }
            else
            {
                fprintf(test_file, "%s", tokens[i*NUM_COLUMNS + j]);
            }
            if (j == NUM_COLUMNS - 1)
            {
                fprintf(test_file, "\n");
            }
            else
            {
                fprintf(test_file, ",");
            }
        }
    }
    fclose(test_file);
    return;
}


/*Parse a csv file, assuming that the first line in the file contains
  headers for each of the columns.*/
int test_parse_csv_simple()
{
    write_test_file(tokens);
    int num_lines;
    int num_columns;
    int ignore_headers = 0;
    char ** data;
    rc_check(parse_csv(path, &num_lines, &num_columns, ignore_headers, &data));
    if (num_lines != NUM_RECORDS)
    {
        fprintf(stderr, "The number lines %d does not match %d.\n",
                num_lines, NUM_RECORDS);
        return GRTCODE_VALUE_ERR;
    }
    if (num_columns != NUM_COLUMNS)
    {
        fprintf(stderr, "The number of columns %d does not match %d.\n",
                num_columns, NUM_COLUMNS);
        return GRTCODE_VALUE_ERR;
    }
    int i;
    for (i=0; i<NUM_RECORDS; ++i)
    {
        int j;
        for (j=0; j<NUM_COLUMNS; ++j)
        {
            if (strcmp(data[j*NUM_RECORDS + i], tokens[i*NUM_COLUMNS + j]) != 0)
            {
                fprintf(stderr, "Tokens %s and %s do not match.\n",
                        data[j*NUM_RECORDS + i], tokens[i*NUM_COLUMNS + j]);
                return GRTCODE_VALUE_ERR;
            }
        }
    }
    return GRTCODE_SUCCESS;
}


/*Parse a csv file that has a blank line.*/
int test_parse_csv_blank_line()
{
    write_test_file(blank_line_tokens);
    int num_lines;
    int num_columns;
    int ignore_headers = 0;
    char ** data;
    return parse_csv(path, &num_lines, &num_columns, ignore_headers, &data);
}


/*Parse a csv file that has a line that is too long.*/
int test_parse_csv_too_long_line()
{
    write_test_file(too_long_line_tokens);
    int num_lines;
    int num_columns;
    int ignore_headers = 0;
    char ** data;
    return parse_csv(path, &num_lines, &num_columns, ignore_headers, &data);
}


/*Parse a csv file that has a line that has the wrong number of columns.*/
int test_parse_csv_wrong_number_of_columns()
{
    write_test_file(wrong_columns_tokens);
    int num_lines;
    int num_columns;
    int ignore_headers = 0;
    char ** data;
    return parse_csv(path, &num_lines, &num_columns, ignore_headers, &data);
}


/*Parse a csv file that is blank..*/
int test_parse_csv_blank_file()
{
    write_test_file(blank_file_tokens);
    int num_lines;
    int num_columns;
    int ignore_headers = 0;
    char ** data;
    return parse_csv(path, &num_lines, &num_columns, ignore_headers, &data);
}


/*Parse a csv file that only contains headers, while ignoring headers.*/
int test_parse_csv_headers_only_while_being_ignored()
{
    write_test_file(headers_only_tokens);
    int num_lines;
    int num_columns;
    int ignore_headers = 1;
    char ** data;
    return parse_csv(path, &num_lines, &num_columns, ignore_headers, &data);
}


/*Parse a csv file that contains a token that is too long.*/
int test_parse_csv_token_too_long()
{
    write_test_file(token_too_long_tokens);
    int num_lines;
    int num_columns;
    int ignore_headers = 1;
    char ** data;
    return parse_csv(path, &num_lines, &num_columns, ignore_headers, &data);
}


/*Test main function.*/
int main(void)
{
    Test_t tests[NUM_TESTS] = {
        {
            test_parse_csv_simple,
            "test_parse_csv_simple",
            GRTCODE_SUCCESS
        },
        {
            test_parse_csv_blank_line,
            "test_parse_csv_blank_line",
            GRTCODE_VALUE_ERR
        },
        {
            test_parse_csv_too_long_line,
            "test_parse_csv_too_long_line",
            GRTCODE_VALUE_ERR
        },
        {
            test_parse_csv_wrong_number_of_columns,
            "test_parse_csv_wrong_number_of_columns",
            GRTCODE_VALUE_ERR
        },
        {
            test_parse_csv_blank_file,
            "test_parse_csv_blank_file",
            GRTCODE_VALUE_ERR
        },
        {
            test_parse_csv_headers_only_while_being_ignored,
            "test_parse_csv_headers_only_while_being_ignored",
            GRTCODE_VALUE_ERR
        },
        {
            test_parse_csv_token_too_long,
            "test_parse_csv_token_too_long",
            GRTCODE_RANGE_ERR
        }
    };
    return test_harness("test_parse_csv", NUM_TESTS, tests);
}
