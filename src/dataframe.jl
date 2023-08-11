
using DataFrames

function get_signal(
    reader::BigWig.Reader,
    df_region::DataFrame,
    n::Int = 1;
    usezoom::Bool = false,
    bystrand::Bool = true,
    chrom_col::Symbol = :chrom,
    start_col::Symbol = :start,
    end_col::Symbol = :end,
    strand_col::Symbol = :strand,
    inherit_cols::Vector{Symbol} = Array{Symbol}([:chrom, :start, :end, :strand,:name,:id]),
    prefix::String = "bin_",
    missing_as::Union{Float32, Missing} = missing,
    
)
    @assert n >= 1
    @assert chrom_col in names(df_region)
    @assert start_col in names(df_region)
    @assert end_col in names(df_region)
    if bystrand
        @assert strand_col in names(df_region)
    end

    chroms = df_region[!, chrom_col]
    starts = df_region[!, start_col]
    ends = df_region[!, end_col]
    if bystrand
        strands = df_region[!, strand_col]
    else
        strands = fill('+', length(chroms))
    end

    values = Array{Union{Float32, Missing},2}(missing, (n,length(chroms), ))

    Threads.@threads for i in eachindex(chroms)
        values[:, i] = mean_n(
            reader,
            chroms[i],
            starts[i],
            ends[i],
            n;
            usezoom = usezoom,
            strand = strands[i],
        )
    end
    values = permutedims(values, [2, 1])
    if !(ismissing(missing_as))
        values[ismissing.(values)] .= missing_as
    end
    
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
    reader::BigWig.Reader,
    df_region::DataFrame,
    n::Int = 1;
    usezoom::Bool = false,
    bystrand::Bool = true,
    chrom_col::Symbol = :chrom,
    start_col::Symbol = :start,
    end_col::Symbol = :end,
    strand_col::Symbol = :strand,
    inherit_cols::Symbol = :all,
    prefix::String = "bin_",
    missing_as::Union{Float32, Missing} = missing,
    
)
    @assert inherit_cols in [:all, :none]
    if inherit_cols == :all
        inherit_cols = names(df_region)
    else
        inherit_cols = Symbol[]
    end
    return get_signal(
        reader,
        df_region,
        n;
        usezoom = usezoom,
        bystrand = bystrand,
        chrom_col = chrom_col,
        start_col = start_col,
        end_col = end_col,
        strand_col = strand_col,
        inherit_cols = inherit_cols,
        prefix = prefix,
        missing_as = missing_as,
    )
end