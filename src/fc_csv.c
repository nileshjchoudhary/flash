#include "fc.h"

CSV_LINE * get_CSV_line(char *line, char *delimiter)
{
    CSV_LINE *csv_line;
    char *tok;
    char **value;
    char delimiter_end[10];

    if (delimiter == NULL) {
        delimiter = " ";
    }
    sprintf(delimiter_end, "%s\n", delimiter);

    csv_line = malloc(sizeof(*csv_line));
    csv_line->n = 0;
    csv_line->value = malloc(10240 * sizeof(*value));

    for (tok = strtok(line, delimiter); tok && *tok; tok = strtok(NULL, delimiter_end))
    {
        int nvalue;

        nvalue = csv_line->n;

        csv_line->value[nvalue] = malloc(10240 * sizeof(char));

        strcpy(csv_line->value[nvalue], tok);
        csv_line->n ++;              
    }

    return csv_line;
}

void free_CSV_line(CSV_LINE **csv_line)
{
    int i;
    CSV_LINE *csv_line0 = *csv_line;

    for (i = 0; i < csv_line0->n; i++) {
        free(csv_line0->value[i]);
    }

    free(csv_line0->value);
    free(*csv_line);
}

