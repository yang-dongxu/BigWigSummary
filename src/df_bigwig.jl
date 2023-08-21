
using DataFrames
summary = Base.summary

function process_scale(reader::BigWig.Reader, chroms::Array{T}, scale::Union{Bool, T, Symbol}, values ::Array{Union{Float32, Missing},2}) where T<:AbstractString
    """
    scale: Union{Bool, String, Symbol}
    true -> "global"
    false -> "none"
    global: scale by global mean
    local: scale by chromosome-wise mean
    none: no scale
    """
    if scale == true
        scale = "global"
    elseif scale == false
        scale = "none"
    end
    scale = Symbol(scale)
    @assert scale in [:none, :global, :local]
    @assert size(values, 1) == length(chroms)
    if scale == :none
        return values
    elseif scale == :global
        global_mean = summary(reader)["all"]
        values ./= global_mean
    elseif scale == :local
        s_dict = summary(reader)
        s = Array{Float32}(undef, length(chroms))
        @inbounds for i in 1:length(chroms)
            s[i] = s_dict[chroms[i]]
        end
        values ./= s

    end
    return values
    
end

function process_inherit(df_region::DataFrame, cols::Union{Symbol, Vector{Symbol}})::Vector{String}

    if cols == :all
        inherit_cols = names(df_region)
    elseif cols == :none
        inherit_cols = String[]
    else
        inherit_cols = [String(i) for i in cols]
    end

    cols_inherit = intersect(names(df_region), inherit_cols)
    # logging
    cols_not_inherit = setdiff(names(df_region), inherit_cols)
    if !isempty(cols_not_inherit)
        println("cols not inherit: ", cols_not_inherit)
    end
    return cols_inherit
end

function get_signal(
    reader::BigWig.Reader,
    df_region::DataFrame, 
    n::Int = 1,
    args...;
    width::Int = -1,
    usezoom::Bool = false,
    bystrand::Bool = true,
    chrom_col::Union{Symbol,T} = :chrom,
    start_col::Union{Symbol,T} = :start,
    end_col::Union{Symbol,T} = :end,
    strand_col::Union{Symbol,T} = :strand,
    inherit_cols::Union{Symbol,Vector{Symbol}} = Array{Symbol}([:chrom, :start, :end, :strand,:name,:id]),
    prefix::String = "bin_",
    missing_as::Union{Float32, Missing} = missing,
    scale::Union{Bool, T, Symbol} = false,
    threads::Bool = true,
)::DataFrame where T <: AbstractString
    @assert length(args) == 0 "you should not pass positional arguments other than reader, df_region and bins, use keyword arguments instead (separated by ; after positional arguments). You passed in: $args"
    @assert n >= 1
    @assert String(chrom_col) in names(df_region)
    @assert String(start_col) in names(df_region)
    @assert String(end_col) in names(df_region)

    chroms = df_region[:, chrom_col]
    starts = df_region[:, start_col]
    ends = df_region[:, end_col]

    if bystrand
        @assert String(strand_col) in names(df_region)
        strands = df_region[:, strand_col]
    else
        strands = fill('+', length(chroms))
    end

    if width > 0
        w = width ÷ 2
        center = (starts .+ ends) .÷ 2
        starts = center .- w
        starts = max.(starts, 0)
        ends = center .+ w
    end

    values = Array{Union{Float32, Missing},2}(missing, (n,length(chroms), ))

    if threads
        Threads.@threads for i in eachindex(chroms)
            values[:, i] = mean_n(
                reader,
                chroms[i],
                starts[i],
                ends[i],            
                strands[i];
                n = n,
                zoom = usezoom,

            )
        end
    else
        for i in eachindex(chroms)
            values[:, i] = mean_n(
                reader,
                chroms[i],
                starts[i],
                ends[i],            
                strands[i];
                n = n,
                zoom = usezoom,

            )
        end
    end
    values = permutedims(values, [2, 1])
    if !(ismissing(missing_as))
        values[ismissing.(values)] .= missing_as
    end

    values = process_scale(reader, chroms, scale, values)
    
    cols_signal = [Symbol("$prefix$i") for i in 1:n]
    df_signal = DataFrame(values, cols_signal)

    
    # filter out cols not exists in inherit_cols
    cols_inherit = process_inherit(df_region, inherit_cols)
    # merge
    df_signal = hcat(df_region[:, cols_inherit], df_signal)


    return df_signal
end

# function get_signal(
#     reader::BigWig.Reader, df_region::DataFrame,
#     args...;
#     inherit_cols::Symbol = :all,
#     kwargs...
    
# )
#     @assert inherit_cols in [:all, :none]
#     if inherit_cols == :all
#         inherit_cols = names(df_region)
#     else
#         inherit_cols = Symbol[]
#     end
#     return get_signal(
#         reader::BigWig.Reader,df_region::DataFrame, args...;
#         inherit_cols = inherit_cols,
#         kwargs...
#     )
# end


function get_signal_flank(
    reader::BigWig.Reader, df_region::DataFrame,
    n5::Int = 1, n::Int=1, n3::Int = 1; 
    left::Int = 2000, right::Int = 2000, bystrand = true, 
    chrom_col::Union{Symbol,T} = :chrom,
    start_col::Union{Symbol,T} = :start,
    end_col::Union{Symbol,T} = :end,
    strand_col::Union{Symbol,T} = :strand,
    prefix::String = "bin_",
    kwargs...
    ) where T <: AbstractString
    @assert n5 >= 0
    @assert n3 >= 0
    @assert n >= 1
    


    kwargs = Dict(kwargs)
    kwargs[:width] = -1

    df_c = get_signal(
        reader, df_region, n; 
        bystrand = bystrand, 
        chrom_col = chrom_col,
        start_col = start_col,
        end_col = end_col,
        strand_col = strand_col,
        prefix = "$(prefix)body_",
        kwargs...
        )

    if n5 == 0 && n3 == 0
        return df_c
    end

    if n5 != 0 
        df_l = begin 
            to_interval(df_region; chrom_col = chrom_col, start_col = start_col, end_col = end_col, strand_col = strand_col) |> 
            x -> GenomicFeatures.IntervalCollection(x, true) |> 
            x-> resize(x, left; fix = :_5, bystrand = bystrand, bydirection = true) |> 
            to_df
        end

        df_5 = get_signal(
            reader, df_l, n5; 
            bystrand = bystrand, 
            chrom_col = chrom_col,
            start_col = start_col,
            end_col = end_col,
            strand_col = strand_col,
            prefix = "$(prefix)5w$(left)_",
            kwargs...
            )

        df_c_n = hcat(df_5, df_c[!,Regex("$(prefix)body_.*")])
        df_c_n[:,chrom_col] = df_c[:,chrom_col]
        df_c_n[:,start_col] = df_c[:,start_col]
        df_c_n[:,end_col] = df_c[:,end_col]
        df_c_n[:,strand_col] = df_c[:,strand_col]
        df_c = df_c_n

    end
    if n3 != 0
        df_r = begin 
            to_interval(df_region; chrom_col = chrom_col, start_col = end_col, end_col = end_col, strand_col = strand_col) |> 
            x-> resize(x, right; fix = :_3, bystrand = bystrand, bydirection = true) |> 
            to_df
        end
        df_3 = get_signal(
            reader, df_r, n3; 
            bystrand = bystrand, 
            chrom_col = chrom_col,
            start_col = start_col,
            end_col = end_col,
            strand_col = strand_col,
            prefix = "$(prefix)3w$(right)_",
            kwargs...
            )
        df_c = hcat(df_c,df_3[!,Regex("$(prefix)3w.*")])
    end
    return df_c

end

