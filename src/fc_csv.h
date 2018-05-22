#ifndef FC_CSV
typedef struct CSV_LINE_ {
    int n;
    char **value;
} CSV_LINE;

CSV_LINE * get_CSV_line(char *line, char *delimiter);
void free_CSV_line(CSV_LINE **csv_line);
#define FC_CSV
#endif
