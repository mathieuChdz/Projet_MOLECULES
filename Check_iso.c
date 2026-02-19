#include "nauty.h"
#include <stdio.h>
#include <stdlib.h>

void print_signature(char* path) {
    FILE *f = fopen(path, "r");
    if (!f) return;

    int n;
    if (fscanf(f, "%d", &n) != 1) {
        fclose(f);
        return;
    }

    int m = SETWORDSNEEDED(n);

    DYNALLSTAT(graph, g, g_sz);
    DYNALLSTAT(int, lab, lab_sz);
    DYNALLSTAT(int, ptn, ptn_sz);
    DYNALLSTAT(int, orbits, orbits_sz);
    static DEFAULTOPTIONS_GRAPH(options);
    statsblk stats;

    DYNALLOC1(int, lab, lab_sz, n, "malloc");
    DYNALLOC1(int, ptn, ptn_sz, n, "malloc");
    DYNALLOC1(int, orbits, orbits_sz, n, "malloc");
    DYNALLOC2(graph, g, g_sz, n, m, "malloc");

    int *colors = malloc(n * sizeof(int));
    for (int i = 0; i < n; i++) {
        fscanf(f, "%d", &colors[i]);
    }

    int current_lab = 0;
    for (int c = 0; c < 255; c++) {
        int count_in_color = 0;
        
        for (int i = 0; i < n; i++) {
            if (colors[i] == c) {
                lab[current_lab] = i;
                ptn[current_lab] = 0;
                current_lab++;
                count_in_color++;
            }
        }
        

        if (count_in_color > 0) {
            ptn[current_lab - 1] = 1; 
        }
    }

    // --- 3. Construction du Graphe ---
    EMPTYGRAPH(g, m, n);
    int u, v;
    while (fscanf(f, "%d %d", &u, &v) == 2) {
        ADDONEEDGE(g, u, v, m);
        ADDONEEDGE(g, v, u, m);
    }
    fclose(f);


    options.getcanon = 1;
    options.defaultptn = FALSE;

    DYNALLSTAT(graph, canong, canong_sz);
    DYNALLOC2(graph, canong, canong_sz, n, m, "malloc");


    densenauty(g, lab, ptn, orbits, &options, &stats, m, n, canong);


    for (int i = 0; i < n; i++) {
        printf("%d_", colors[lab[i]]);
    }
    printf("|");

    for (int i = 0; i < m * n; i++) {
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
