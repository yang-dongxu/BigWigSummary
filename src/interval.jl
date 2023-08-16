import GenomicFeatures;
using GenomicFeatures: Interval, IntervalCollection, Strand
import DataFrames
using DataFrames: DataFrame, DataFrameRow, index

function GenomicFeatures.Interval(interval::Interval{T})::Interval{T} where T<:Any
    return Interval(interval.seqname, interval.first, interval.last, interval.strand, interval.metadata)
end

function shift(interval::Interval, shift::Integer; bystrand::Bool = true)::Interval
    if bystrand
        shift = shift * (isneg(interval) ? -1 : 1)
    end
    # interval.first += shift
    # interval.last += shift
    return Interval(interval.seqname, interval.first + shift, interval.last + shift, interval.strand, interval.metadata)
end



function sort_boundary(interval::Interval)::Interval
    if interval.first > interval.last
        # interval.first, interval.last = interval.last, interval.first
        interval = Interval(interval.seqname, interval.last, interval.first, interval.strand, interval.metadata)
    end
    return interval
end


function resize(interval::Interval, size::Integer; fix::Union{Symbol,Char}=:center, bystrand::Bool = true)::Interval
    # if bystrand
    #     size = size * (isneg(interval) ? 1 : -1)
    # end
    @assert fix in [:center, :left, :right,:c,:l,:r,'c','l','r', :_5, :_3,'5','3']
    if bystrand == false
        if fix == :_5 || fix == '5'
            fix = :left
        elseif fix == :_3 || fix == '3'
            fix = :right
        end
    end
    f, l = interval.first, interval.last

    if fix == :center || fix == :c || fix == 'c'
        center = (interval.first + interval.last) ÷ 2
        f = center - size ÷ 2
        l = center + size ÷ 2
    elseif fix == :left || fix == :l || fix == 'l'
        f = interval.first
        l = interval.first + size 

    elseif fix == :right || fix == :r || fix == 'r'
        f = interval.last - size 
        l = interval.last
    elseif fix == :_5 || fix == '5'
        if isneg(interval)
            f = interval.last - size 
        else
            l = interval.first + size 
        end
    elseif fix == :_3 || fix == '3'
        if isneg(interval)
            l = interval.first + size 
        else
            f = interval.last - size 
        end
    else
        error("fix must be :first, :last or :center")
    end
    f, l = min(f, l), max(f, l)
    interval = Interval(interval.seqname, f, l, interval.strand, interval.metadata)

    return interval
end



function pad(interval::Interval, pad::Integer; direction::Union{Symbol, Char}=:both, bystrand::Bool = true)::Interval
    # if bystrand
    #     pad = pad * (interval.strand == '+' ? 1 : -1)
    # end
    p = [:both, :left, :right,:b,:l,:r,'b','l','r', :_5, :_3,'5','3']
    @assert direction in p
    if bystrand == false
        if direction == :_5 || direction == '5'
            direction = :left
        elseif direction == :_3 || direction == '3'
            direction = :right
        end
    end
    f, l = interval.first, interval.last
    if direction == :both || direction == :b || direction == 'b'
        f -= pad
        l += pad
    elseif direction == :left || direction == :l || direction == 'l'
        f -= pad
    elseif direction == :right || direction == :r || direction == 'r'
        l += pad
    elseif direction == :_5 || direction == '5'
        if ispos(interval)
            f -= pad
        else
            l += pad
        end
    elseif direction == :_3 || direction == '3'
        if ispos(interval)
            l += pad
        else
            f -= pad
        end 
    else
        error("direction must be $p")
    end

    return Interval(interval.seqname, f, l, interval.strand, interval.metadata)
end

function pad(interval::Interval, pad::AbstractFloat; direction::Symbol=:both, bystrand::Bool = true)::Interval
    if pad < -1
        error("pad must be >= -1")
    end
    size = (interval.last - interval.first + 1) * pad
    return pad!(interval, size, direction = direction, bystrand = bystrand)
end


function shift(intervals::IntervalCollection, args...;kwargs...)
    return IntervalCollection([shift(i, args...; kwargs...) for i in intervals])
end

function resize(intervals::IntervalCollection, args...;kwargs...)
    return IntervalCollection([resize(i, args...; kwargs...) for i in intervals])
end

function pad(intervals::IntervalCollection, args...;kwargs...)
    return IntervalCollection([pad(i, args...; kwargs...) for i in intervals])
end


function to_interval(row::DataFrameRow; chrom_col::Union{Symbol,String} = :chrom, start_col::Union{Symbol,String} = :start, end_col::Union{Symbol,String} = :end, strand_col::Union{Symbol,String} = :strand):: Interval{DataFrameRow{DataFrame, DataFrames.Index}}
    if !(String(strand_col) in names(row))
        s = '.'
    else
        s = Strand(row[strand_col])
    end
    return Interval(row[chrom_col],row[start_col], row[end_col], s, row)
end


function to_interval(df::DataFrame; chrom_col::Union{Symbol,String} = :chrom, start_col::Union{Symbol,String} = :start, end_col::Union{Symbol,String} = :end, strand_col::Union{Symbol,String} = :strand)::IntervalCollection
    @assert String(chrom_col) in names(df)
    @assert String(start_col) in names(df)
    @assert String(end_col) in names(df)
    # @assert String(strand_col) in names(df)

    intervals = Vector{Interval{DataFrameRow{DataFrame, DataFrames.Index}}}(undef, size(df, 1))
    for i in 1:size(df, 1)
        intervals[i] = to_interval(df[i, :], chrom_col = chrom_col, start_col = start_col, end_col = end_col, strand_col = strand_col)
    end
    return IntervalCollection(intervals, true) # sorted

end




function to_df(intervals::IntervalCollection; sync::Bool = true,  chrom_col::Union{Symbol,String} = :chrom, start_col::Union{Symbol,String} = :start, end_col::Union{Symbol,String} = :end, strand_col::Union{Symbol,String} = :strand)::DataFrame
    s = Array{DataFrameRow{DataFrame, DataFrames.Index}}(undef, length(intervals))
    for (i,v) in enumerate(intervals)
        s[i] = v.metadata
        if sync
            s[i][chrom_col] = v.seqname
            s[i][start_col] = v.first
            s[i][end_col] = v.last
            s[i][strand_col] = v.strand
        end
    end
    return DataFrame(s)
end

function to_df(
    intervals::Vector{GenomicFeatures.Interval{DataFrameRow{DataFrame, DataFrames.Index}}};kwargs...
    )
    return to_df(IntervalCollection(intervals), kwargs...)
end