# test reading 2-pt function data and store them in delimited files
using Revise
using DataConverter
using DelimitedFiles

path_dir = "/Users/alessandroconigli/Lattice/data/hq-ssf/3pt/N300"

path_data = filter(x->occursin(".dat", x), readdir(path_dir, join=true))

cdata = read_mesons(path_data, "G5", "G5")

cdata[1]

cdata[1].header

# 1 replicum
write_2pt_data(cdata[4][1], @__DIR__, "test_3pt", "4" )

# multiple replicas
write_2pt_data(cdata[1], @__DIR__, "test_multi", ["4", "5"] )