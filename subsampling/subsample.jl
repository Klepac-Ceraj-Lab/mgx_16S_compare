# activate and instantiate project
using Pkg
Pkg.activate(@__DIR__)
Pkg.instantiate()

using BioSequences
using FASTX
using GZip
using CodecZlib

# path to raw samples on `ada`
## can set environmental variable to change
if haskey(ENV, "FASTQPATH") 
    rawfastq_path = ENV["FASTQPATH"]    
else
    # default value
    rawfastq_path = "/lovelace/echo/sequencing/mgx/rawfastq/"
end

outpath = joinpath(@__DIR__, "output/")
# create directory if it doesn't exist
isdir(outpath) || mkdir(outpath)
# file provided by Danielle
samples = readlines(joinpath(@__DIR__, "samples.txt"))
fastqs = readdir(rawfastq_path)

"""
    subsample(records, size)

Get subsample of `BioSequence` records with `size` reads.
Returns a view into the original vector,
or if `size` is greater than the total number of reads,
returns the original `records` vector.
"""
function subsample(records, size)
    l = length(records)
    if size < l
        idx = rand(1:l, size)
        return view(records, idx)
    else
        return records
    end
end

"""
    write_subsample(path, ss)

Writes a GZip-compressed vector `ss` of `BioSequences` to `path`
"""
function write_subsample(path, ss)
    io = FASTQ.Writer(GzipCompressorStream(open(path, "w")))
    for s in ss
        write(io, s)
    end
    close(io)
end

function main()
    # Loop through samples, find all fastqs and read in their sequence records
    for sample in samples
        records = FASTQ.Record[]
        for fastq in filter(f-> occursin(replace(sample, "_"=>"-"), f) || occursin(sample, f), fastqs)
            file = joinpath(rawfastq_path, fastq)
            append!(records, FASTQ.Reader(GzipDecompressorStream(open(file))))
        end
        
        # Make 3 replicates each of 10k, 100k, and 1M random subsamples
        # TODO: could make these parameters programmable
        for i in 1:4
            write_subsample(joinpath(outpath, "$(sample)_250k_$i.fastq.gz"), subsample(records, 250_000))
            write_subsample(joinpath(outpath, "$(sample)_500k_$i.fastq.gz"), subsample(records, 500_000))
            write_subsample(joinpath(outpath, "$(sample)_750k_$i.fastq.gz"), subsample(records, 750_000))
        end
    end
end

main()