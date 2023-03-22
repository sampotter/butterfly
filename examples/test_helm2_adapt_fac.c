#include <bf/fac_streamer.h>

#include "argtable3.h"

int const MAX_NUM_ARG_ERRORS = 20;

typedef struct {
  BfReal k;
} Opts;

static int parseArgs(int argc, char *argv[], Opts *opts) {
  int code = 0, nerrors = 0;
  char const *progname = argv[0];

  struct arg_lit *help;
  struct arg_lit *verbose;
  struct arg_dbl *wavenumber;
  struct arg_end *end;

  void *argtable[] = {
    help = arg_litn(NULL, "help", 0, 1, "display this help and exit"),
    verbose = arg_litn("v", "verbose", 0, 1, "verbose output"),
    wavenumber = arg_dbln("k", "wavenumber", "<k>", 0, 1, "the wavenumber for the problem"),
    end = arg_end(MAX_NUM_ARG_ERRORS),
  };
  nerrors = arg_parse(argc, argv, argtable);

  /* special case: '--help' takes precedence over error reporting */
  if (help->count > 0) {
    printf("Usage: %s", progname);
    arg_print_syntax(stdout, argtable, "\n");
    printf("Demonstrate command-line parsing in argtable3.\n\n");
    arg_print_glossary(stdout, argtable, "  %-25s %s\n");
    code = 0;
    goto exit;
  }

  /* If the parser returned any errors then display them and exit */
  if (nerrors > 0) {
    /* Display the error details contained in the arg_end struct.*/
    arg_print_errors(stdout, end, progname);
    printf("Try '%s --help' for more information.\n", progname);
    code = 1;
    goto exit;
  }

  opts->k = *wavenumber->dval;

exit:
  arg_freetable(argtable, sizeof(argtable) / sizeof(argtable[0]));

  return code;
}

int main(int argc, char *argv[]) {
  Opts opts;
  int code = parseArgs(argc, argv, &opts);
  if (!code) return code;

  printf("k = %g\n", opts.k);
}
