#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "argparse.h"


#ifdef __cplusplus
#define EXTERNC extern "C"
#else
#define EXTERNC
#endif


/** @brief Check for whitespace and starting dashes in an argument name.*/
static void check_name(char const * const name, /**< Argument name.*/
                       int * const one_dash, /**< Flag telling if name starts with -.*/
                       int * const two_dash, /**< Flag telling if name starts with --.*/
                       int const max_len /**< Maximum allowed name length.*/
                      )
{
    *one_dash = 0;
    *two_dash = 0;
    size_t i;
    if (strlen(name) > (size_t)(max_len-1))
    {
        fprintf(stderr, "[Error]: argument %s must be less than %d characters.\n",
                name, max_len-1);
        exit(EXIT_FAILURE);
    }
    for (i=0; i<strlen(name); ++i)
    {
        if (i == 0 && name[i] == '-')
        {
            *one_dash = 1;
        }
        else if (i == 1 && name[i] == '-' && *one_dash)
        {
            *two_dash = 1;
        }
        else if (isspace((int)(name[i])))
        {
            fprintf(stderr, "[Error]: Argument name cannot have whitespace.\n");
            exit(EXIT_FAILURE);
        }
    }
    return;
}


/** @brief Get the basename of a path.*/
static char * path_basename(char *path, /**< File path.*/
                            char seperator /**< Path seperator.*/
                           )
{
    size_t index = 0;
    size_t i;
    for (i=0; i<strlen(path); ++i)
    {
        if (path[i] == seperator)
        {
            index = i;
        }
    }
    return path + index + 1;
}


/** @brief Print usage message.*/
static void print_usage(Parser_t const p /**< Parser.*/
                       )
{
    printf("\033[3mUsage: %s ", path_basename(p.argv[0], '/'));
    Argument_t *a = p.args;
    while (a != NULL)
    {
        if (a->positional)
        {
            printf("%s ", a->name);
        }
        else if (a->requires_value)
        {
            printf("[%s <value>] ", a->name);
        }
        else
        {
            printf("[%s] ", a->name);
        }
        a = a->head;
    }
    printf("\n\033[0m");
    return;
}


/** @brief Print help message.*/
static void print_help(Parser_t const p /**< Parser.*/
                      )
{
    print_usage(p);
    printf("\033[3m\n%s - %s\n", path_basename(p.argv[0], '/'), p.description);
    Argument_t *a = p.args;
    if (p.num_pos_args > 0)
    {
        printf("\nPositional arguments:\n");
        while (a != NULL)
        {
            if (a->positional)
            {
                char pad[namelen+1];
                size_t i;
                for (i=0; i<(namelen+longnamelen-strlen(a->name)); ++i)
                {
                    pad[i] = ' ';
                }
                pad[i] = '\0';
                printf("%s%s    %s\n", a->name, pad, a->description);
            }
            a = a->head;
        }
    }
    printf("\nOptional arguments:\n");
    a = p.args;
    while (a != NULL)
    {
        if (!(a->positional))
        {
            char pad[namelen+longnamelen+1];
            printf("%s", a->name);
            if (strlen(a->longname) > 0)
            {
                printf(", %s", a->longname);
                size_t i;
                for (i=0;i<(namelen+longnamelen-strlen(a->name)-strlen(a->longname)); ++i)
                {
                    pad[i] = ' ';
                }
                pad[i] = '\0';
            }
            else
            {
                printf("  ");
                size_t i;
                for (i=0;i<(namelen+longnamelen-strlen(a->name)); ++i)
                {
                    pad[i] = ' ';
                }
                pad[i] = '\0';
            }
            printf("%s  %s\n", pad, a->description);
        }
        a = a->head;
    }
    printf("\033[0m");
    return;
}


/*Create a parser and add help option.*/
EXTERNC Parser_t create_parser(int const argc, char **argv, char const * const description)
{
    Parser_t p;
    p.argc = argc;
    p.argv = argv;
    if (description != NULL)
    {
        snprintf(p.description,desclen, "%s", description);
    }
    else
    {
        snprintf(p.description, desclen, "%s", "");
    }
    p.args = NULL;
    p.num_args = 0;
    p.num_pos_args = 0;
    add_argument(&p, "-h", "--help", "Print this help message", 0);
    return p;
}


/*Add argument to parser.*/
EXTERNC void add_argument(Parser_t * const parser, char const * const name,
                          char const * const longname, char const * const description,
                          int const * const requires_value)
{
    Argument_t *arg = (Argument_t *)malloc(sizeof(*arg));
    arg->head = NULL;
    arg->found = 0;
    int one_dash;
    int two_dash;
    check_name(name, &one_dash, &two_dash, namelen);
    snprintf(arg->name, namelen, "%s", name);
    if (two_dash)
    {
        fprintf(stderr, "[Error]: Specify longname using second argument.\n");
        exit(EXIT_FAILURE);
    }
    if (one_dash)
    {
        arg->positional = 0;
        if (requires_value != NULL)
        {
            arg->requires_value = *requires_value;
        }
        else
        {
            arg->requires_value = 0;
        }
        if (longname != NULL)
        {
            check_name(longname, &one_dash, &two_dash, longnamelen);
            if (!two_dash)
            {
                fprintf(stderr, "[Error]: Longname must start with --.\n");
                exit(EXIT_FAILURE);
            }
            snprintf(arg->longname, longnamelen, "%s", longname);
        }
        else
        {
            snprintf(arg->longname, longnamelen, "%s", "");
        }
    }
    else
    {
        arg->requires_value = 1;
        arg->positional = 1;
        (parser->num_pos_args)++;
    }
    if (description != NULL)
    {
        snprintf(arg->description, desclen, "%s", description);
    }
    else
    {
        snprintf(arg->longname, desclen, "%s", "");
    }
    if (parser->args == NULL)
    {
        parser->args = arg;
    }
    else
    {
        Argument_t *p = parser->args;
        while (p->head != NULL)
        {
            p = p->head;
        }
        p->head = arg;
    }
    (parser->num_args)++;
    return;
}


/*Parse arguments.*/
EXTERNC void parse_args(Parser_t const p)
{
    int num_pos_args = 0;
    int i;
    for (i=1; i<p.argc; ++i)
    {
        int one_dash = 0;
        int two_dash = 0;
        check_name(p.argv[i], &one_dash, &two_dash, valuelen);
        if (two_dash || one_dash)
        {
            if (strcmp(p.argv[i], "-h") == 0 || strcmp(p.argv[i], "--help") == 0)
            {
                print_help(p);
                exit(EXIT_SUCCESS);
            }
            Argument_t *a = p.args;
            while (a != NULL)
            {
                char *n;
                if (two_dash)
                {
                    n = a->longname;
                }
                else
                {
                    n = a->name;
                }
                if (strcmp(p.argv[i], n) == 0)
                {
                    if (a->requires_value)
                    {
                        if (i+1 < p.argc)
                        {
                            i++;
                            check_name(p.argv[i], &one_dash, &two_dash, valuelen);
                            if (one_dash || two_dash)
                            {
                                print_usage(p);
                                exit(EXIT_FAILURE);
                            }
                            snprintf(a->value, valuelen, "%s", p.argv[i]);
                        }
                        else
                        {
                            print_usage(p);
                            exit(EXIT_FAILURE);
                        }
                    }
                    a->found = 1;
                    break;
                }
                a = a->head;
            }
        }
        else
        {
            num_pos_args++;
            Argument_t *a = p.args;
            while (a != NULL)
            {
                if (a->positional && !a->found)
                {
                    a->found = 1;
                    snprintf(a->value, valuelen, "%s", p.argv[i]);
                    break;
                }
                a = a->head;
            }
        }
    }
    if (num_pos_args != p.num_pos_args)
    {
        print_usage(p);
        exit(EXIT_FAILURE);
    }
    return;
}


/*Get the value of an argument.*/
EXTERNC int get_argument(Parser_t const p, char const * const name, char buffer[valuelen])
{
    int result = 0;
    Argument_t *a = p.args;
    while (a != NULL)
    {
        if (strcmp(name, a->name) == 0)
        {
            result = a->found;
            if (a->found && a->requires_value)
            {
                snprintf(buffer, valuelen, "%s", "");
                snprintf(buffer, valuelen, "%s", a->value);
            }
            break;
        }
        a = a->head;
    }
    return result;
}


/*Free memory reserved by parser.*/
EXTERNC void destroy_parser(Parser_t * const p)
{
    Argument_t *a = p->args;
    Argument_t *b ;
    while (a != NULL)
    {
        b = a;
        a = a->head;
        free(b);
    }
    return;
}
