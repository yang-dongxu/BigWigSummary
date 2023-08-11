using GenomicFeatures;


function shift!(intervel::Interval, shift::Int; bystrand::Bool = true)::Interval
    if bystrand
        shift = shift * (intervel.strand == '+' ? 1 : -1)
    end
    intervel.first += shift
    intervel.last += shift
    return intervel
end

function shift(intervel::Interval, shift::Int, bystrand::Bool = true)::Interval
    intervel = copy(intervel)
    return shift!(intervel, shift, bystrand = bystrand)
end

function sort_boundary!(interval::Interval)::Interval
    if interval.first > interval.last
        interval.first, interval.last = interval.last, interval.first
    end
    return interval
end

function resize!(intervel::Interval, size::Int; fix::Symbol=:center, bystrand::Bool = true)::Interval
    if bystrand
        size = size * (intervel.strand == '+' ? 1 : -1)
    end
    if fix == :first
        intervel.last = intervel.first + size
    elseif fix == :last
        intervel.first = intervel.last - size
    elseif fix == :center 
        center = (intervel.first + intervel.last) ÷ 2
        intervel.first = center - size ÷ 2
        intervel.last = center + size ÷ 2
    else
        error("fix must be :first, :last or :center")
    end
    if bystrand
        sort_boundary!(intervel)
    end
    return intervel
end


function resize(intervel::Interval, size::Int, fix::Symbol=:center, bystrand::Bool = true)::Interval
    intervel = copy(intervel)
    return resize!(intervel, size, fix = fix, bystrand = bystrand)
end

