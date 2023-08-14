
using DataFrames
summary = Base.summary

function process_scale(reader::BigWig.Reader, chroms::Array{String}, scale::Union{Bool, String, Symbol}, values ::Array{Union{Float32, Missing},2})
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

function get_signal(
    reader::BigWig.Reader,
    df_region::DataFrame, 
    n::Int = 1,
    args...;
    usezoom::Bool = false,
    bystrand::Bool = true,
    chrom_col::Union{Symbol,String} = :chrom,
    start_col::Union{Symbol,String} = :start,
    end_col::Union{Symbol,String} = :end,
    strand_col::Union{Symbol,String} = :strand,
    inherit_cols::Vector{Symbol} = Array{Symbol}([:chrom, :start, :end, :strand,:name,:id]),
    prefix::String = "bin_",
    missing_as::Union{Float32, Missing} = missing,
    scale::Union{Bool, String, Symbol} = false,
    threads::Bool = true,
)::DataFrame
    @assert length(args) == 0 "you should not pass positional arguments other than reader, df_region and bins, use keyword arguments instead (separated by ; after positional arguments). You passed in: $args"
    @assert n >= 1
    @assert String(chrom_col) in names(df_region)
    @assert String(start_col) in names(df_region)
    @assert String(end_col) in names(df_region)

    chroms = df_region[!, chrom_col]
    starts = df_region[!, start_col]
    ends = df_region[!, end_col]

    if bystrand
        @assert String(strand_col) in names(df_region)
        strands = df_region[!, strand_col]
    else
        strands = fill('+', length(chroms))
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
    cols_inherit = intersect(names(df_region), inherit_cols)
    # logging
    cols_not_inherit = setdiff(names(df_region), inherit_cols)
    if !isempty(cols_not_inherit)
        println("cols not inherit: ", cols_not_inherit)
    end
    # merge
    df_signal = hcat(df_region[:, cols_inherit], df_signal)


    return df_signal
end

function get_signal(
    reader::BigWig.Reader, df_region::DataFrame,
    args...;
    inherit_cols::Union{Symbol,Vector{Symbol}} = :all,
    kwargs...
    
)
    @assert inherit_cols in [:all, :none]
    if inherit_cols == :all
        inherit_cols = names(df_region)
    else
        inherit_cols = Symbol[]
    end
    return get_signal(
        reader::BigWig.Reader,df_region::DataFrame, args...;
        inherit_cols = inherit_cols,
        kwargs...
    )
end