/* GRTCODE is a GPU-able Radiative Transfer Code
 * Copyright (C) 2016  Garrett Wright
 * Modified in 2019 by Raymond Menzel
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation; version 2.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "debug.h"
#include "parse_csv.h"
#include "return_codes.h"
#include "utilities.h"

#define MAXCHARSPERLINE 1024
#define MAXCHARSPERTOKEN 32
#define offset(col, line, nlines) (col*nlines + line)


/*Copy an input token into an input array.*/
static int copy_token(char * const token, /*Token.*/
                      int const col_index, /*Column index.*/
                      int const line_index, /*Line index.*/
                      int const num_lines, /*Total number of lines.*/
                      char **vals /* Array of tokens.*/
                     )
{
    int s = strlen(token);
    in_range(s, 1, MAXCHARSPERTOKEN-1);
    if (token[s-1] == '\n')
    {
        /*Strip off the trailing new line character.*/
        token[s-1] = '\0';
    }
    strncpy(vals[offset(col_index, line_index, num_lines)], token,
            (size_t)MAXCHARSPERTOKEN);
    return GRTCODE_SUCCESS;
}


/*Parse a csv file, assuming that the first line in the file contains
  headers for each of the columns.*/
int parse_csv(char const * const filepath, int * const num_lines, int * const num_cols,
              int const ignore_headers, char *** out)
{
    not_null(filepath);
    not_null(num_lines);
    not_null(num_cols);
    not_null(out);
    FILE *f = NULL;
    catch(open_file(&f, filepath, "r"));
    char const *mesg = "Reading csv file %s.";
    log_info(mesg, filepath);

    /*Count the number of lines/columns on a line.*/
    char line[MAXCHARSPERLINE];
    *num_lines = 0;
    *num_cols = -1;
    while (fgets(line, MAXCHARSPERLINE, f) != NULL)
    {
        (*num_lines)++;
        char *c = line;
        if (*c == '\n')
        {
            mesg = "line %d in file %s is blank.";
            raise(GRTCODE_VALUE_ERR, mesg, *num_lines, filepath);
        }
        int n = 1;
        int num_chars = 0;
        while (*c != '\n')
        {
            if (*c == ',')
            {
                n++;
            }
            c++;
            num_chars++;
            if (num_chars > MAXCHARSPERLINE)
            {
                mesg = "the number of characters (>=%d) on line %d of file %s"
                       " exceeds the maximum allowed (%d).";
                raise(GRTCODE_VALUE_ERR, mesg, num_chars, *num_lines, filepath,
                      MAXCHARSPERLINE);
            }
        }
        if (*num_cols < 0)
        {
            *num_cols = n;
        }
        else if (n != *num_cols)
        {
            mesg = "the number of columns (%d) on line %d of file"
                   " %s differs from the number of columns (%d) on"
                   " the other lines.";
            raise(GRTCODE_VALUE_ERR, mesg, n, *num_lines, filepath, *num_cols);
        }
    }
    if (*num_lines == 0)
    {
        mesg = "the file %s is empty.";
        raise(GRTCODE_VALUE_ERR, mesg, filepath);
    }

    /*Allocate arrays to hold the values that will be read in.*/
    rewind(f);
    char **vals;
    if (ignore_headers)
    {
        (*num_lines)--;
        if (*num_lines == 0)
        {
            mesg = "the file %s only contains headers, no data.";
            raise(GRTCODE_VALUE_ERR, mesg, filepath);
        }
        fgets(line, MAXCHARSPERLINE, f);
    }
    int num_vals = (*num_cols)*(*num_lines);
    gmalloc(vals, num_vals, HOST_ONLY);
    int i;
    for (i=0; i<num_vals; ++i)
    {
        gmalloc(vals[i], MAXCHARSPERTOKEN, HOST_ONLY);
        snprintf(vals[i], MAXCHARSPERTOKEN, "%c", '\0');
    }

    /*Read in the data.*/
    int line_index = 0;
    while (fgets(line, MAXCHARSPERLINE, f) != NULL)
    {
        int col_index = 0;
        char *token = strtok(line, ",");
        catch(copy_token(token, col_index, line_index, *num_lines, vals));
        while (1)
        {
            token = strtok(NULL, ",");
            if (token == NULL)
            {
                break;
            }
            col_index++;
            catch(copy_token(token, col_index, line_index, *num_lines, vals));
        }
        line_index++;
    }
    *out = vals;

    /*Close the file.*/
    if (0 != fclose(f))
    {
        mesg = "failed to close csv file %s.";
        raise(GRTCODE_VALUE_ERR, mesg, filepath);
    }
    return GRTCODE_SUCCESS;
}
