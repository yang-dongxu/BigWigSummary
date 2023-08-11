import BBI
using BigWig 
using LoopVectorization

function read(bwfile::AbstractString)::BigWig.Reader
    return BigWig.Reader(BigWig.open(bwfile))
end


function read_intervals(
    reader::BigWig.Reader,
    chrom::AbstractString, chromstart::T, chromend::T
    )::BigWig.OverlapIterator where T<:Integer
    chromid = reader.chroms[chrom][1]
    return BigWig.OverlapIterator(reader, chromid, chromstart, chromend)
end

function split_n_parts(a::T,b::T,n::Integer)::Array{T,1} where T<:Integer
    @assert n >= 1
    @assert a <= b
    l =  LinRange(a,b,n) 
    return round.(T,l)
end

function coverage2_1base(x::NTuple{2,T1}, y::NTuple{2,T2}) where T1<:Integer where T2<:Integer
    x1, x2 = x
    y1, y2 = y
    return max(0, min(x2, y2) - max(x1, y1)) + T(1) |> T1
end

function coverage2_0base(x::NTuple{2,T1}, y::NTuple{2,T2}) where T1<:Integer where T2<:Integer
    x1, x2 = x
    y1, y2 = y
    return max(0, min(x2, y2) - max(x1, y1)) |> T1
end


function mean_n(
    reader::BigWig.Reader, 
    chrom::AbstractString, chromstart::Integer, chromend::Integer, n::Int = 1;
    strand::Union{Char,Symbol} = '+', zoom::Bool = false
    )::Array{Union{Float32,Missing}}
    @assert n >= 1
    @assert chromstart <= chromend
    @assert n <= chromend - chromstart
    chromid = reader.chroms[chrom][1]

    if n == 1
        v = [BigWig.mean(reader, chrom, chromstart + 1, chromend,usezoom=zoom)]
    else
        if zoom
            zooms = BBI.find_best_zoom(reader.zooms, Int64(chromend - chromstart))
            if zooms !== nothing
                v = mean_n_2(zooms, Int64(chromid), Int64(chromstart), Int64(chromend), n; strand=strand)
                # println("zoom")
            else
                v = mean_n_2(reader, Int64(chromid), Int64(chromstart), Int64(chromend), n; strand=strand)
            end
        else
            v = mean_n_2(reader, Int64(chromid), Int64(chromstart), Int64(chromend), n; strand=strand)
        end
    end
    return v
end

function mean_n_2(
    handler::Union{BigWig.Reader,BBI.Zoom},
    chromid::Int64, chromstart::Int64, chromend::Int64,
    n::Int = 1; strand::Union{Char,Symbol} = '+'
)::Array{Union{Float32,Missing}}
    if typeof(handler) <: BigWig.Reader
        method = BigWig.OverlapIterator
    elseif typeof(handler) <: BBI.Zoom
        method = BBI.find_overlapping_zoomdata
    else
        error("handler must be BigWig.Reader or BBI.ZoomData")
    end

    @assert n >  1
    @assert chromstart <= chromend
    @assert n <= chromend - chromstart
    @assert strand in ['+',:pos,:forward,'-',:neg,:reverse]

    # println("chromstart: $chromstart, chromend: $chromend, n: $n")
    sums = zeros(Union{Float32,Missing},n)
    sizes = zeros(Int,n)

    p = split_n_parts(chromstart,chromend,n+1)
    starts = p[1:end-1] 
    ends = p[2:end]
    # println("starts: $(Int.(starts)), ends: $(Int.(ends))")

    o = method(handler, chromid, chromstart, chromend) |> collect
    sort!(o, by=x->x.chromstart)
    @inbounds for i in 1:n
        for record in o
            cov = coverage2_0base((record.chromstart, record.chromend), (starts[i],ends[i]))
            if cov > 0
                # println("i: $i, cov: $cov")
                sums[i] = sums[i] + record.value * cov
                sizes[i] = sizes[i] + cov
            elseif record.chromstart > ends[i] # no further overlap
                break
            end
        end
        # println("i: $i, sums: $(sums[i]), sizes: $(sizes[i])")

    end

    @inbounds for i in 1:n
        if sizes[i] > 0
            sums[i] = sums[i]/sizes[i]
        else
            sums[i] = missing
        end
    end
    if strand in ['-',:neg,:reverse]
        reverse!(sums)
    end
    return sums    
end