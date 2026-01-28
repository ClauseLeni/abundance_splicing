for i in {1..7}

do 
julia spliced_coverage.jl -o At${i}_splicing.counts At${i}nsorted_chl_cat.bam
done 

