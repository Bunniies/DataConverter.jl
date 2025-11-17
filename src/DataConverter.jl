module DataConverter

using ADerrors, DelimitedFiles
import Statistics: mean
import ADerrors: err

include("DataConverter_types.jl")
export CData

include("DataConverter_reader.jl")
export read_mesons, read_mesons_multiple_files

include("DataConverter_write.jl")
export write_2pt_data

end # module DataConverter
