
import GenomicFeatures
using GenomicFeatures: Strand, STRAND_POS, STRAND_NEG, STRAND_NA ,STRAND_BOTH, Interval

function Base.convert(::Type{Strand}, s::String)
    @assert s in ["+","-", "pos", "neg", "forward", "reverse", ".","unk","unknown"]
    if s in ["+","pos", "forward"]
        return Strand('+')
    elseif s in ["-","neg", "reverse"]
        return Strand('-')
    else
        return Strand('.')
    end
end


function Base.convert(::Type{Strand}, s::Symbol)
    @assert s in [:+,:-,:pos,:neg,:forward,:reverse,:.,:unk,:unknown]
    s = String(s)
    return Base.convert(Strand, s::String) 
end

function Base.convert(::Type{Symbol}, s::Strand)
    @assert s in [STRAND_POS, STRAND_NEG, STRAND_NA, STRAND_BOTH]
    if s == STRAND_POS
        return :+
    elseif s == STRAND_NEG
        return :-
    elseif s == STRAND_NA
        return :.
    else
        return :.
    end
end



function ispos(strand::Strand)
    return strand == STRAND_POS
end


function ispos(interval::Interval)
    return ispos(interval.strand)
end

function isneg(strand::Strand)
    return strand == STRAND_NEG
end

function isneg(interval::Interval)
    return isneg(interval.strand)
end

function isstrand(strand::Strand)
    return strand != STRAND_NA || strand != STRAND_BOTH
end

function isstrand(interval::Interval)
    return isstrand(interval.strand)
end

function GenomicFeatures.Strand(strand::String)
    return Base.convert(Strand, strand)
end

function GenomicFeatures.Strand(strand::Symbol)
    return Base.convert(Strand, strand)
end


function  GenomicFeatures.Strand(strand::Strand)
    return strand
end

