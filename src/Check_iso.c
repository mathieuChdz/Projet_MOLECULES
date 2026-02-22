# include "nauty.h"
# include <stdio.h>

void print_signature(char * path) {
    FILE * f = fopen(path, "r");
    int n, e;
    if (fscanf(f, "%d %d", & n, & e) != 2) return;

    DYNALLSTAT(graph, g, g_sz);
    DYNALLSTAT(int, lab, lab_sz);
    DYNALLSTAT(int, ptn, ptn_sz);
    DYNALLSTAT(int, orbits, orbits_sz);
    static DEFAULTOPTIONS_GRAPH(options);
    statsblk stats;

    DYNALLOC1(int, lab, lab_sz, n, "malloc");
    DYNALLOC1(int, ptn, ptn_sz, n, "malloc");
    DYNALLOC1(int, orbits, orbits_sz, n, "malloc");
    DYNALLOC2(graph, g, g_sz, n, SETWORDSNEEDED(n), "malloc");

    int * colors = malloc(n * sizeof(int));
    for (int i=0; i < n; i++) {
        fscanf(f, "%d", & colors[i]);
    }

    int current_lab = 0;
    int color_list[200];
    int color_count = 0;

    for (int c=0; c < 255; c++) {
        int found = 0;
        for (int i=0; i < n; i++) {
            if (colors[i] == c) {
                lab[current_lab] = i;
                ptn[current_lab] = 0;
                current_lab++;
                found = 1;
            }
        }
        if (found) ptn[current_lab - 1] = 1; 
    }


    options.getcanon = 1;
    options.defaultptn = FALSE;
    EMPTYGRAPH(g, SETWORDSNEEDED(n), n);

    int u, v;
    for (int i = 0; i < e; i++) {
        if (fscanf(f, "%d %d", &u, &v) == 2) {
            ADDONEEDGE(g, u, v, SETWORDSNEEDED(n));
            ADDONEEDGE(g, v, u, SETWORDSNEEDED(n));
        }
    }
    fclose(f);

    DYNALLSTAT(graph, canong, canong_sz);
    DYNALLOC2(graph, canong, canong_sz, n, SETWORDSNEEDED(n), "malloc");
    densenauty(g, lab, ptn, orbits, &options, &stats, SETWORDSNEEDED(n), n, canong);


    for (int i = 0; i < n; i++) {
        printf("%d_", colors[lab[i]]);
    }
    printf("|");

    for (int i = 0; i < SETWORDSNEEDED(n) * n; i++) {
        printf("%lx", canong[i]);
    }
    printf("\n");

    free(colors); 
    DYNFREE(g, g_sz); DYNFREE(canong, canong_sz);
    DYNFREE(lab, lab_sz); DYNFREE(ptn, ptn_sz); DYNFREE(orbits, orbits_sz);
}

int main(int argc, char *argv[]) {
    if (argc < 2) return 1;
    print_signature(argv[1]);
    return 0;
}
