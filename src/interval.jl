import GenomicFeatures;
using GenomicFeatures: Interval, IntervalCollection, Strand
import DataFrames
using DataFrames: DataFrame, DataFrameRow, index

function GenomicFeatures.Interval(interval::Interval{T})::Interval{T} where T<:Any
    return Interval(interval.seqname, interval.first, interval.last, interval.strand, interval.metadata)
end

function get_boundary(interval::Interval; side::Union{Symbol,Char}=:center, bystrand::Bool = true)::Tuple{Integer, Integer}
    p = [:center, :left, :right, :both,:c,:l,:r,:b,'c','l','r','b', :_5, :_3,'5','3']
    @assert side in p
    if bystrand == false
        if side == :_5 || side == '5'
            side = :left
        elseif side == :_3 || side == '3'
            side = :right
        end
    end

    if side == :center || side == :c || side == 'c'
        o1 = (interval.first + interval.last) ÷ 2
    elseif side == :left || side == :l || side == 'l'
        o1 = interval.first
    elseif side == :right || side == :r || side == 'r'
        o1 = interval.last
    elseif side == :_5 || side == '5'
        if isneg(interval)
            o1 = interval.last 
        else
            o1 = interval.first 
        end
    elseif side == :_3 || side == '3'
        if isneg(interval)
            o1 = interval.first 
        else
            o1 = interval.last 
        end
    else
        o1 = interval.first
    end
    o2 = (interval.first == o1) ? interval.last : interval.first
    return o1, o2
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


function resize(interval::Interval, size::Integer; fix::Union{Symbol,Char}=:center, bystrand::Bool = true, bydirection::Bool = true)::Interval
    """
    interval: input intervals
    size: final region size, 
    fix: which boundary should be kept. 
        :left, :right will ignore strand of the interval
        :_5, :3 will consider strand of the interval if bystrand is true
    bydirection (works only when fix is '5' or '3'):
        extend direction shoud consider interval or not?
        if false, size >0 means right and <0 means left 
        if true, size >0 means same direction and vice versa.
    """
    if bydirection && fix in [:_5,:_3, '5','3']
        sign = (isneg(interval) ? -1 : 1) * ((fix in [:_5,'5']) ? -1 : 1)
        size = size * sign
    end
    f,l = get_boundary(interval, side = fix, bystrand = bystrand)
    l = max(0, f + size)
    f, l = min(f, l), max(f, l)
    interval = Interval(interval.seqname, f, l, interval.strand, interval.metadata)

    return interval
end

function pad(interval::Interval, size::Integer; side::Union{Symbol, Char}=:both, bystrand::Bool = true, bydirection::Bool = false)::Interval

    f,l = get_boundary(interval, side = side, bystrand = bystrand)
    size = (f < l) ? size : -size

    if side == :both || side == 'b'
        @assert l -f > 2*abs(size) "$(interval.seqname):$f-$l is smaller than 2*pad = $(2*abs(size))"
        l = max(0, l + size)
    end
    f = max(0,f - size)

    f, l = min(f, l), max(f, l)
    interval = Interval(interval.seqname, f, l, interval.strand, interval.metadata)

    return interval
end

function pad(interval::Interval, frac::AbstractFloat; side::Symbol=:both, bystrand::Bool = true)::Interval
    if pad < -1
        error("pad must be >= -1")
    end
    size = (interval.last - interval.first ) * frac
    return pad!(interval, size, side = side, bystrand = bystrand)
end


function shift(intervals::IntervalCollection, args...;kwargs...)
    return IntervalCollection([shift(i, args...; kwargs...) for i in intervals], true)
end

function resize(intervals::IntervalCollection, args...;kwargs...)
    return IntervalCollection([resize(i, args...; kwargs...) for i in intervals], true)
end

function pad(intervals::IntervalCollection, args...;kwargs...)
    return IntervalCollection([pad(i, args...; kwargs...) for i in intervals], true)
end


function to_interval(row::DataFrameRow; chrom_col::Union{Symbol,String} = :chrom, start_col::Union{Symbol,String} = :start, end_col::Union{Symbol,String} = :end, strand_col::Union{Symbol,String} = :strand):: Interval{DataFrameRow{DataFrame, DataFrames.Index}}
    if !(String(strand_col) in names(row))
        s = '.'
    else
        s = Strand(row[strand_col])
    end

    @assert row[start_col] < row[end_col] row
    return Interval(row[chrom_col],row[start_col], row[end_col], s, row)
end


function to_interval(df::DataFrame; chrom_col::Union{Symbol,String} = :chrom, start_col::Union{Symbol,String} = :start, end_col::Union{Symbol,String} = :end, strand_col::Union{Symbol,String} = :strand)::IntervalCollection
    @assert String(chrom_col) in names(df)
    @assert String(start_col) in names(df)
    @assert String(end_col) in names(df)
    # @assert String(strand_col) in names(df)
    df = deepcopy(df)
    intervals = Vector{Interval{DataFrameRow{DataFrame, DataFrames.Index}}}(undef, size(df, 1))
    for i in eachindex(intervals)
        row = df[i,:]
        if !(String(strand_col) in names(row))
            s = '.'
        else
            s = Strand(row[strand_col])
        end
        intervals[i] = Interval(row[chrom_col],row[start_col], row[end_col], s, row)
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

    df = DataFrame()
    for c in names(s[1])
        df[!,Symbol(c)] = map(x -> x[c], s)
    end
    return df
end

function to_df(
    intervals::Vector{GenomicFeatures.Interval{DataFrameRow{DataFrame, DataFrames.Index}}};kwargs...
    )
    return to_df(IntervalCollection(intervals, true), kwargs...)
end