#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "argparse.h"
#include "netcdf.h"


#define nc_catch(return_code) { \
    int return_code_ = return_code; \
    if (return_code_ != NC_NOERR) {\
        fprintf(stderr, "[%s: %d] %s\n", __FILE__, __LINE__, nc_strerror(return_code_)); \
        exit(EXIT_FAILURE); \
    } \
}


struct Dimension {
    char axis;
    char name[128];
    size_t length;
    size_t total_length;
    int dimid;
};


struct Variable {
    char axis[128];
    struct Dimension * dimensions[32];
    int dimids[32];
    char name[128];
    int num_attributes;
    int num_dimensions;
    char standard_name[256];
    nc_type type;
    char units[128];
    int varid;
};


struct Dimension * get_dimension(char const * name, struct Dimension * dimensions,
                                 int num_dimensions)
{
    int i;
    for (i=0; i<num_dimensions; ++i)
    {
        if (strcmp(name, dimensions[i].name) == 0)
        {
            return &(dimensions[i]);
        }
    }
    return NULL;
}


struct Variable * get_variable(char const * name, struct Variable * variables,
                               int num_variables)
{
    int i;
    for (i=0; i<num_variables; ++i)
    {
        if (strcmp(name, variables[i].name) == 0)
        {
            return &(variables[i]);
        }
    }
    return NULL;
}


int main(int argc, char **argv)
{
    char *description = "Combines output files from GRTCODE.";
    Parser_t parser = create_parser(argc, argv, description);
    add_argument(&parser, "output_file", NULL, "Output file path.", NULL);
    add_argument(&parser, "input_file", NULL, "Index 0 file path.", NULL);
    add_argument(&parser, "num_files", NULL, "Number of file segments.", NULL);
    parse_args(parser);

    /*Get the number of files.*/
    char buffer[valuelen];
    get_argument(parser, "num_files", buffer);
    int num_files = atoi(buffer);

    /*Get the directory and basename of the input file path.*/
    get_argument(parser, "input_file", buffer);
    int i = valuelen - 1;
    while (buffer[i] != '/')
    {
        i--;
    }
    buffer[i] = '\0';
    char directory[512];
    snprintf(directory, 512, "%s", buffer);

    int starting_point = i + 1;
    i = valuelen - 1;
    while (buffer[i] != '.')
    {
        i--;
    }
    buffer[i] = '\0';
    char basename[512];
    snprintf(basename, 512, "%s", &(buffer[starting_point]));
    fprintf(stderr, "directory: %s\n", directory);
    fprintf(stderr, "input file basename: %s\n", basename);

    /*Gather the necessary metadata about all of the variables.*/
    int num_dimensions = 0;
    struct Dimension dimensions[32];
    int num_variables = 0;
    struct Variable variables[128];

    for (i=1; i<=num_files; ++i)
    {
        /*Construct the path name for and open the specific file segment.*/
        char path[512];
        snprintf(path, 512, "%s/%s.%d", directory, basename, i);
        fprintf(stderr, "Opening dataset: %s\n", path);
        int ncid_segment;
        nc_catch(nc_open(path, NC_NOWRITE, &ncid_segment));

        /*Get the number of dimensions and store the dimension metadata.*/
        int num_dimensions_in_input;
        nc_catch(nc_inq_ndims(ncid_segment, &num_dimensions_in_input));
        int j;
        for (j=0; j<num_dimensions_in_input; ++j)
        {
            char name[128];
            memset(name, 0, 128*sizeof(char));
            size_t length;
            nc_catch(nc_inq_dim(ncid_segment, j, name, &length));

            struct Dimension * dimension = get_dimension(name, dimensions, num_dimensions);
            if (dimension == NULL)
            {
                /*Copy the newly found dimension.*/
                dimensions[num_dimensions].axis = '-';
                dimensions[num_dimensions].length = length;
                memset(dimensions[num_dimensions].name, 0, 128*sizeof(char));
                snprintf(dimensions[num_dimensions].name, 128, "%s", name);
                dimensions[num_dimensions].total_length = 0;
                fprintf(stderr, "Found dimension %d: %s %p.\n",
                        num_dimensions, name, &(dimensions[num_dimensions]));
                num_dimensions++;
                if (num_dimensions > 32)
                {
                    fprintf(stderr, "Error: Too many dimensions.");
                    return EXIT_FAILURE;
                }
            }
        }

        /*Get the number of variables and store the variable metadata.*/
        int num_variables_in_input;
        nc_catch(nc_inq_nvars(ncid_segment, &num_variables_in_input));
        for (j=0; j<num_variables_in_input; ++j)
        {
            /*Get the variables metadata.*/
            char name[128];
            memset(name, 0, 128*sizeof(char));
            nc_type type;
            int num_variable_dimensions;
            int dimids[32];
            int num_attributes;
            nc_catch(nc_inq_var(ncid_segment, j, name, &type, &num_variable_dimensions,
                                dimids, &num_attributes));
            char axis[128];
            memset(axis, 0, 128*sizeof(char));
            nc_get_att_text(ncid_segment, j, "axis", axis);
            char units[128];
            memset(units, 0, 128*sizeof(char));
            nc_get_att_text(ncid_segment, j, "units", units);
            char standard_name[256];
            memset(standard_name, 0, 256*sizeof(char));
            nc_get_att_text(ncid_segment, j, "standard_name", standard_name);

            /*Update the dimensions metadata if this variable is a dimension.*/
            struct Dimension * dimension = get_dimension(name, dimensions, num_dimensions);
            if (dimension != NULL)
            {
                fprintf(stderr, "Found variable for dimension %d: %s %p.\n",
                        num_dimensions, name, dimension);
                size_t length;
                char buffer[128];
                nc_catch(nc_inq_dim(ncid_segment, dimids[0], buffer, &length));
                if (strlen(axis) > 0)
                {
                    dimension->axis = axis[0];
                    if (dimension->axis == 'x' || dimension->axis == 'X' ||
                        dimension->axis == 'y' || dimension->axis == 'Y')
                    {
                        dimension->total_length += length;
                    }
                    else if (dimension->length != length)
                    {
                        fprintf(stderr, "Error: incompatiable dimension lengths for %s.\n",
                                name);
                        return EXIT_FAILURE;
                    }
                }
                if (dimension->total_length == 0)
                {
                    dimension->total_length = length;
                }
            }

            struct Variable * variable = get_variable(name, variables, num_variables);
            if (variable != NULL)
            {
                fprintf(stderr, "Found existing variable %d: %s %p.\n",
                        num_variables, name, variable);
                /*Found existing variable - check that the variable in this file is
                  consistent.*/
                if (variable->type != type)
                {
                    fprintf(stderr, "Error: Incompatible variable types for %s.\n", name);
                    return EXIT_FAILURE;
                }
                if (variable->num_dimensions != num_variable_dimensions)
                {
                    fprintf(stderr, "Error: Incompatible number of dimensions for %s.\n", name);
                    return EXIT_FAILURE;
                }
                int m;
                for (m=0; m<num_variable_dimensions; ++m)
                {
                    char dimension_name[128];
                    memset(dimension_name, 0, 128*sizeof(char));
                    nc_catch(nc_inq_dimname(ncid_segment, dimids[m], dimension_name));
                    struct Dimension * dimension = get_dimension(dimension_name, dimensions,
                                                                 num_dimensions);
                    fprintf(stderr, "Variables: %s  dimid: %d  dimension: %p\n",
                            name, variable->dimensions[m], dimension);
                    if (variable->dimensions[m] != dimension)
                    {
                        fprintf(stderr, "Error: Incompatible dimension order for %s.\n", name);
                        return EXIT_FAILURE;
                    }
                }
            }
            else
            {
                /*Found new variable.  Store its metadata.*/
                fprintf(stderr, "Found new variable %d: %s %p.\n",
                        num_variables, name, variables[num_variables]);
                memset(variables[num_variables].name, 0, 128*sizeof(char));
                snprintf(variables[num_variables].name, 128, "%s", name);
                variables[num_variables].type = type;
                variables[num_variables].num_dimensions = num_variable_dimensions;
                int m;
                for (m=0; m<num_variable_dimensions; ++m)
                {
                    char dimension_name[128];
                    memset(dimension_name, 0, 128*sizeof(char));
                    nc_catch(nc_inq_dimname(ncid_segment, dimids[m], dimension_name));
                    struct Dimension * dimension = get_dimension(dimension_name, dimensions,
                                                                 num_dimensions);
                    variables[num_variables].dimensions[m] = dimension;
                    fprintf(stderr, "%d: %p  == %p?\n", m,
                            variables[num_variables].dimensions[m], dimension);
                }
                variables[num_variables].num_attributes = num_attributes;
                if (num_variables > 128)
                {
                    fprintf(stderr, "Error: Too many variables.");
                    return EXIT_FAILURE;
                }
                memset(variables[num_variables].units, 0, 128*sizeof(char));
                snprintf(variables[num_variables].units, 128, "%s", units);
                memset(variables[num_variables].axis, 0, 128*sizeof(char));
                snprintf(variables[num_variables].axis, 128, "%s", axis);
                memset(variables[num_variables].standard_name, 0, 256*sizeof(char));
                snprintf(variables[num_variables].standard_name, 256, "%s", standard_name);
                num_variables++;
            }
        }
        /*Close the segment file.*/
        nc_catch(nc_close(ncid_segment));
    }

    /*Create the output file.*/
    get_argument(parser, "output_file", buffer);
    int ncid;
    fprintf(stderr, "Creating dataset: %s\n", buffer);
    nc_catch(nc_create(buffer, NC_CLOBBER, &ncid));
    fprintf(stderr, "ncid: %d\n", ncid);

    /*Create the dimensions in the output file.*/
    for (i=0; i<num_dimensions; ++i)
    {
        nc_catch(nc_def_dim(ncid, dimensions[i].name, dimensions[i].total_length,
                 &(dimensions[i].dimid)));
    }

    /*Create the variables in the output file.*/
    for (i=0; i<num_variables; ++i)
    {
        int dimids[32];
        memset(dimids, -1, 32*sizeof(int));
        int j;
        for (j=0; j<variables[i].num_dimensions; ++j)
        {
            dimids[j] = variables[i].dimensions[j]->dimid;
        }
        nc_catch(nc_def_var(ncid, variables[i].name, variables[i].type,
                            variables[i].num_dimensions, dimids,
                            &(variables[i].varid)));
        nc_catch(nc_put_att_text(ncid, variables[i].varid, "units",
                                 strlen(variables[i].units) + 1, variables[i].units));
        nc_catch(nc_put_att_text(ncid, variables[i].varid, "standard_name",
                                 strlen(variables[i].standard_name) + 1,
                                 variables[i].standard_name));
    }

    /*Exit define mode.*/
    nc_catch(nc_enddef(ncid));

    /*Loop through the input files and write their chunks of data.*/
    for (i=1; i<=num_files; ++i)
    {
        /*Construct the path name for and open the specific file segment.*/
        char path[512];
        memset(path, 0, 512*sizeof(char));
        snprintf(path, 512, "%s/%s.%d", directory, basename, i);
        fprintf(stderr, "Opening dataset: %s\n", path);
        int ncid_segment;
        nc_catch(nc_open(path, NC_NOWRITE, &ncid_segment));

        /*Get global attributes.*/
        int x_start = -1;
        nc_catch(nc_get_att_int(ncid_segment, NC_GLOBAL, "x_start", &x_start));
        int x_stop = -1;
        nc_catch(nc_get_att_int(ncid_segment, NC_GLOBAL, "x_stop", &x_stop));

        int j;
        for (j=0; j<num_variables; ++j)
        {
            int varid_segment;
            nc_catch(nc_inq_varid(ncid_segment, variables[j].name, &varid_segment));

            /* Get the variables hypercube parameters.*/
            int num_variable_dimensions;
            int dimids[32];
            nc_type type;
            nc_catch(nc_inq_var(ncid_segment, varid_segment, NULL, &type,
                                &num_variable_dimensions, dimids, NULL));
            size_t start[32];
            memset(start, 0, sizeof(size_t)*32);
            size_t count[32];
            memset(count, 0, sizeof(size_t)*32);
            int num_elements = 1;
            int k;
            for (k=0; k<num_variable_dimensions; ++k)
            {
                char name[128];
                memset(name, 0, 128*sizeof(char));
                size_t length;
                nc_catch(nc_inq_dim(ncid_segment, dimids[k], name, &length));
                struct Dimension * dimension = get_dimension(name, dimensions, num_dimensions);
                count[k] = length;
                num_elements *= length;
                start[k] = 0;
                if (dimension != NULL)
                {
                    if (dimension->axis == 'x' && x_start >= 0 && x_stop >= 0)
                    {
                        start[k] = x_start;
                    }
                }
            }
            if (type == NC_DOUBLE)
            {
                double * buffer = malloc(num_elements*sizeof(buffer));
                nc_catch(nc_get_var_double(ncid_segment, variables[j].varid, buffer));
                nc_catch(nc_put_vara_double(ncid, variables[j].varid, start, count, buffer));
                free(buffer);
            }
        }
        /*Close the segment file.*/
        nc_catch(nc_close(ncid_segment));
    }
    /*Close the output file.*/
    nc_catch(nc_close(ncid));
    return EXIT_SUCCESS;
}
