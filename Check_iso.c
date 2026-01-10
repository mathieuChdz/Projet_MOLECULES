#include "nauty.h"
#include <stdio.h>

void print_signature(char* path) {
    FILE *f = fopen(path, "r");
    int n, e;
    if (fscanf(f, "%d %d", &n, &e) != 2) return;

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

    options.getcanon = 1;
    EMPTYGRAPH(g, SETWORDSNEEDED(n), n);

    int idx, type;
    for (int i = 0; i < n; i++) {
        fscanf(f, "%d %d", &idx, &type);
        lab[i] = i;
        ptn[i] = 1; 
    }
    ptn[n-1] = 0;

    int u, v;
    for (int i = 0; i < e; i++) {
        fscanf(f, "%d %d", &u, &v);
        ADDONEEDGE(g, u, v, SETWORDSNEEDED(n));
    }
    fclose(f);

    DYNALLSTAT(graph, canong, canong_sz);
    DYNALLOC2(graph, canong, canong_sz, n, SETWORDSNEEDED(n), "malloc");
    densenauty(g, lab, ptn, orbits, &options, &stats, SETWORDSNEEDED(n), n, canong);

    for (int i = 0; i < SETWORDSNEEDED(n) * n; i++) {
        printf("%lx", canong[i]);
    }
    printf("\n");

    DYNFREE(g, g_sz); DYNFREE(canong, canong_sz);
    DYNFREE(lab, lab_sz); DYNFREE(ptn, ptn_sz); DYNFREE(orbits, orbits_sz);
}

int main(int argc, char *argv[]) {
    if (argc < 2) return 1;
    print_signature(argv[1]);
    return 0;
}