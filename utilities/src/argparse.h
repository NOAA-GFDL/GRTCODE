/** @file*/

#ifndef ARGPARSE_H_
#define ARGPARSE_H_


#ifdef __cplusplus
#define EXTERNC extern "C"
#else
#define EXTERNC
#endif


/** @brief Maximum allowed description length.*/
#define desclen 512


/** @brief Maximum allowed argument "long name" length.*/
#define longnamelen 24


/** @brief Maximum allowed argument name length.*/
#define namelen 18


/** @brief Maximum allowed argument value length.*/
#define valuelen 512


/** @brief Argument.*/
typedef struct Argument
{
    char name[namelen]; /**< Name (must start with a - if optional.*/
    char longname[longnamelen]; /**< Long name (must start with --).*/
    char description[desclen]; /**< Description that will appear in help message.*/
    int found; /**< Flag telling if the arg was passed in.*/
    char value[valuelen]; /**< Buffer to hold input.*/
    int positional; /**< Flag telling if a positional arg.*/
    int requires_value; /**< Flag telling if the argument requires a value.*/
    struct Argument *head; /**< Head of linked list.*/
} Argument_t;


/** @brief Parser.*/
typedef struct Parser
{
    Argument_t *args; /**< Linked list of arguments.*/
    int argc; /**< Number of arguments given to program.*/
    char **argv; /**< Array of args.*/
    char description[desclen]; /**< Description.*/
    int num_args; /**< Number of arguments added..*/
    int num_pos_args; /**< Number of positional arguments.*/
} Parser_t;


/** @brief Create a parser and add help option.
    @return Parser object.*/
EXTERNC Parser_t create_parser(int const argc, /**< Number of arguments.*/
                               char **argv, /**< Array of arguments.*/
                               char const * const description /**< Description.*/
                              );


/** @brief Add argument to parser.*/
EXTERNC void add_argument(Parser_t * const parser, /**< Parser.*/
                          char const * const name, /**< Short name.*/
                          char const * const longname, /**< Long name.*/
                          char const * const description, /**< Description.*/
                          int const * const requires_value /**< Flag telling to look for value.*/
                         );


/** @brief Parse arguments.*/
EXTERNC void parse_args(Parser_t const p /**< Parser.*/
                       );


/**  @brief Get the value of an argument.
     @return 1 if argument was found, else 0.*/
EXTERNC int get_argument(Parser_t const p, /**< Parser.*/
                         char const * const name, /**< Argument name.*/
                         char buffer[valuelen] /**< Buffer to return argument value in.*/
                        );


/** @brief Free memory reserved by parser.*/
EXTERNC void destroy_parser(Parser_t * const p /**< Parser.*/
                           );


#undef EXTERNC
#endif
