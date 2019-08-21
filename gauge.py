import sys
from dwutil.adj import read_adjacency, write_adjacency
from dwutil.graph import Graph
from typing import List


def gauge_transform(graph: Graph, gauge: List) -> Graph:
    if not len(gauge) == len(graph):
        raise ValueError("Gauge spec not the same length as Graph spec")
    gt = []
    for idx, (n, g) in enumerate(zip(graph, gauge)):
        d = {}
        for m, j in n.items():
            if idx != m:
                d[m] = j * gauge[idx] * gauge[m]
            else:
                d[m] = j * gauge[idx]
        gt.append(d)
    return gt


def bin_list(n: int, bitlen: int) -> List:
    li = []
    for _ in range(bitlen):
        li.append(n & 1)
        n = n >> 1
    return li

def bin_to_sz(l : List[int]) -> List[int]:
    li2 = [1 -2 * b for b in l]
    return li2

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Requires 3 arguments")
        exit(0)

    ad_li = read_adjacency(sys.argv[1])
    bit_len = len(ad_li)

    gauge_str = sys.argv[2]
    outfile = sys.argv[3]

    gauge_int = int(gauge_str, 16)
    gauge_li = bin_to_sz(bin_list(gauge_int, bit_len))
    gt = gauge_transform(ad_li, gauge_li)

    write_adjacency(gt, outfile)




