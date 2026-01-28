struct Feature
    path::String
    start::Int32
    length::Int32
end

# entire set of Features for one strand of one genome
struct FeatureArray
    genome::String
    genome_length::Int32
    strand::Char
    features::Array{Feature,1}
end

function readFeatures(file::String)
    # let f_strand_features,r_strand_features
    open(file) do f
        header = split(readline(f), '\t')
        seqID = header[1]
        genome_length = parse(Int32, header[2])
        f_strand_features = FeatureArray(seqID, genome_length, '+', Array{Feature,1}(undef, 0))
        r_strand_features = FeatureArray(seqID, genome_length, '-', Array{Feature,1}(undef, 0))
        while !eof(f)
            fields = split(readline(f), '\t')
            feature = Feature(fields[1], parse(Int, fields[3]), parse(Int, fields[4]))
            # println(feature)
            if fields[2][1] == '+'
                push!(f_strand_features.features, feature)
            else
                push!(r_strand_features.features, feature)
            end
        end
        return f_strand_features, r_strand_features
    end
end

function isforwardstrand(record::BAM.Record)
    if BAM.flag(record) & SAM.FLAG_READ1 ≠ 0               # if it's read1
        if BAM.flag(record) & SAM.FLAG_REVERSE ≠ 0         # and if it's reverse complement
            return true
        end
    elseif BAM.flag(record) & SAM.FLAG_READ2 ≠ 0           # if it's read2
        if BAM.flag(record) & SAM.FLAG_REVERSE == 0        # and it's forward
            return true
        end
    elseif BAM.flag(record) & SAM.FLAG_REVERSE ≠ 0        # if it's single-end & reverse complement
            return true
    end
    return false
end

function pairpassesmuster(flags1,flags2)
    flags1 & SAM.FLAG_PROPER_PAIR == 0 && return false
    flags1 & SAM.FLAG_SECONDARY ≠ 0 && return false
    flags1 & SAM.FLAG_QCFAIL ≠ 0 && return false
    flags1 & SAM.FLAG_DUP ≠ 0 && return false
    flags1 & SAM.FLAG_SUPPLEMENTARY ≠ 0 && return false

    flags2 & SAM.FLAG_PROPER_PAIR == 0 && return false
    flags2 & SAM.FLAG_SECONDARY ≠ 0 && return false
    flags2 & SAM.FLAG_QCFAIL ≠ 0 && return false
    flags2 & SAM.FLAG_DUP ≠ 0 && return false
    flags2 & SAM.FLAG_SUPPLEMENTARY ≠ 0 && return false
    return true
end

const ONE = Int32(1)

function pairedcountreads(reader::BAM.Reader,mapQ_threshold::Int)

    fwd_read_counts = Vector{Vector{Int32}}(undef,length(reader.refseqlens))
    rev_read_counts = Vector{Vector{Int32}}(undef,length(reader.refseqlens))

    for (i,refseqlen) in enumerate(reader.refseqlens)
        fwd_read_counts[i] = zeros(Int32,refseqlen)
        rev_read_counts[i] = zeros(Int32,refseqlen)
    end
    
    queue = Vector{BAM.Record}(undef,0)
    num_pairs = 0

    for next_record in reader
        !BAM.isfilled(next_record) && continue
        #println("reading ",BAM.tempname(next_record)," refID: ",BAM.flag(next_record))
        push!(queue,next_record)
        length(queue) != 2 && continue
        if BAM.tempname(queue[1]) != BAM.tempname(queue[2])
            popfirst!(queue)
            continue
        end

        record1 = queue[1]
        record2 = queue[2]
        #println("read1 ",BAM.tempname(record1)," ",BAM.flag(record1) & SAM.FLAG_READ1," ",BAM.flag(record1) & SAM.FLAG_READ2)
        #println("read2 ",BAM.tempname(record2)," ",BAM.flag(record2) & SAM.FLAG_READ1," ",BAM.flag(record2) & SAM.FLAG_READ2)

        # ignore bad read pairs (not properly mapped, duplicates etc)
        if !pairpassesmuster(BAM.flag(record1),BAM.flag(record2)) || BAM.mappingquality(record1) < mapQ_threshold || BAM.mappingquality(record2) < mapQ_threshold
            empty!(queue)
            continue
        end

        #println("counting pair ",BAM.tempname(record1)," ",BAM.tempname(record2))
        refindex = BAM.refid(record1) #both reads should have same reference as they are a proper pair

        # progress tracking
        #num_pairs += 1
        #if num_pairs%1000 == 0; print(num_pairs,"\r");end

        read_pair_interval = min(BAM.position(record1),BAM.position(record2)):max(BAM.rightposition(record1),BAM.rightposition(record2))
        if length(read_pair_interval) > 500
            empty!(queue)
            continue #don't count implausibly long read pairs
        end

        @views if isforwardstrand(record2)
            fwd_read_counts[refindex][read_pair_interval] .+= ONE
        else
            rev_read_counts[refindex][read_pair_interval] .+= ONE
        end
        empty!(queue)
    end
    return fwd_read_counts,rev_read_counts
end

function unpairpassesmuster(flags)
    flags & SAM.FLAG_SECONDARY ≠ 0 && return false
    flags & SAM.FLAG_QCFAIL ≠ 0 && return false
    flags & SAM.FLAG_DUP ≠ 0 && return false
    flags & SAM.FLAG_SUPPLEMENTARY ≠ 0 && return false
    return true
end
