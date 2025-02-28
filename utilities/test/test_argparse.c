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


/*Create a parser and add help option.*/
int test_create_parser()
{
    int const argc = 4;
    char * argv[5] = {"test_argparse", "-c", "foo", "-d", "bar"};
    char * description = "Here is a description.";

    Parser_t p = create_parser(argc, argv, description);
    if (p.argc != argc)
    {
        return 
    }
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
