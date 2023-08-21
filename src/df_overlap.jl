import GenomicFeatures;
using GenomicFeatures: Interval, IntervalCollection, Strand
import DataFrames
using DataFrames: DataFrame, DataFrameRow, index


function count_overlaps(row::DataFrameRow, intervals::IntervalCollection; chrom_col::Union{Symbol,String} = :chrom, start_col::Union{Symbol,String} = :start, end_col::Union{Symbol,String} = :end, strand_col::Union{Symbol,String} = :strand)::Integer
    if strand_col in names(row)
        s = Strand(row[strand_col])
    else
        s = '.'
    end
    query = Interval(row[chrom_col],row[start_col], row[end_col], s)
    count = 0 
    for o in GenomicFeatures.eachoverlap(query, intervals)
        count += 1
    end
    return count
end


function count_overlaps(df::DataFrame, intervals::IntervalCollection; chrom_col::Union{Symbol,String} = :chrom, start_col::Union{Symbol,String} = :start, end_col::Union{Symbol,String} = :end, strand_col::Union{Symbol,String} = :strand)::Vector{Integer}
    @assert String(chrom_col) in names(df)
    @assert String(start_col) in names(df)
    @assert String(end_col) in names(df)
    # @assert String(strand_col) in names(df)
    df = deepcopy(df)
    counts = Vector{Int64}(undef, size(df, 1))
    for i in eachindex(counts)
        row = df[i,:]
        counts[i] = count_overlaps(row, intervals; chrom_col = chrom_col, start_col = start_col, end_col = end_col, strand_col = strand_col)
    end
    return counts
end


function count_overlaps(df::DataFrame, df_intervals::DataFrame; chrom_col::Union{Symbol,String} = :chrom, start_col::Union{Symbol,String} = :start, end_col::Union{Symbol,String} = :end, strand_col::Union{Symbol,String} = :strand)::Vector{Integer}
    intervals = to_interval(df_intervals; chrom_col = chrom_col, start_col = start_col, end_col = end_col, strand_col = strand_col)
    return count_overlaps(df, intervals; chrom_col = chrom_col, start_col = start_col, end_col = end_col, strand_col = strand_col)
end