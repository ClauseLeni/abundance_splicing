using ArgParse
using BioSequences
using XAM
using FASTX

include("countreads.jl")

function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table! s begin
        "--mapQ"
            help = "mapping quality threshold"
            arg_type = Int
            default = 0
        "--single-end"
            help = "flag indicating BAM file contains unpaired reads"
            action = :store_true
        "--paired-end"
            help = "flag indicating BAM file contains paired reads"
            action = :store_true
        "--forwardorientation", "-f"
            help = "flag indicating reads (or first read in a pair) are the forward RNA strand, unlike the usual Illumina libraries"
            action = :store_true
        "--outfile","-o"
            help = "path to results file"
            arg_type = String
        "bam"
            help = "bam file sorted by name (e.g using samtools sort -n) containing reads mapped to spliced and unspliced references"
            required = true
    end

    return parse_args(s)
end

function main()
    parsed_args = parse_commandline()
    reader = open(BAM.Reader, parsed_args["bam"])
    mapQ_threshold = parsed_args["mapQ"]
    #defaults to paired_end
    if parsed_args["single-end"]
        fwd_read_counts, rev_read_counts = unpairedcountreads(f_intervals,r_intervals,reader,mapQ_threshold)
    else parsed_args["paired-end"]
        fwd_read_counts, rev_read_counts = pairedcountreads(reader,mapQ_threshold)
    end
    close(reader)

    out = parsed_args["outfile"]
    if isnothing(out)
        io = Base.stdout
    else
        io = open(out, "w")
    end

    for (ref, (reff, refr)) in enumerate(zip(fwd_read_counts, rev_read_counts))
        for (pos, (countf, countr)) in enumerate(zip(reff, refr))
            write(io, reader.refseqnames[ref], "\t", string(pos), "\t", string(countf), "\t", string(countr), "\n")
        end
    end
    close(io)
end

main()
